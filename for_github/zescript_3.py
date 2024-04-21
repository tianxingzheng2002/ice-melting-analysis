"""
coding: utf-8
Beckermann 1999 model

This notebook will do coupled fluid-solid melting with a phase field method
mixed with volume penalisation
"""
import seawater as sw
from dedalus.core.operators import GeneralFunction
# from dedalus.tools import post
from os.path import join
import logging
import numpy as np
import dedalus.public as de
# from dedalus.extras import flow_tools
import time
import sys

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank, size = comm.rank, comm.size
logger = logging.getLogger(__name__)


class ParityFunction(GeneralFunction):
    """ Seawater equation of state Buoyancy function
    Define GeneralFunction subclass to handle parities """

    def __init__(self, domain, layout, func, args=[], kw={}, out=None,
                 parity={},):
        super().__init__(domain, layout, func, args=[], kw={}, out=None,)
        self._parities = parity

    def meta_parity(self, axis):
        return self._parities.get(axis, 1)  # by default, even parity


def sigmoid(x, a=1):
    return 0.5*(np.tanh(x/a)+1)


def well(x, y, z):
    return 1- np.where(x > 0, 1, 0) * np.where(y > 0, 1, 0) * np.where(z > 0, 1, 0)


d = de.operators.differentiate

# Initial Condition Constants
u_0 = 0
w_0 = 0
p_0 = 0
T_0 = 0
C_0 = 0
ft_0 = 0

# Dimensional parameters

UU = 0
T_B = 20  # C
c_p = 4.2  # J/g*C
L_T = 3.34e2  # J/g
C_B = 30  # g/kg
nu = 1.3e-1  # cm**2/s; physical viscosity is 10 times smaller.
kappa = 1.3e-2  # cm**2/s; NEVER USED (physical value 10 times smaller)
mu = 1.3e-3  # cm**2/s; NEVER USED (physical value 10 times smaller)
em = 0.056  # C/(g/kg)
LL, HH = 10, 5  # cm
ll, hh = 2, 2  # cm
epsilon = 4e-2  # cm

# Non-dimensional parameters

Pr = 7
Sc = 50/4
delta = 1e-4
beta = 4/2.648228  # WHAT'S THIS?

# Save parameters
Nx, Nz = 256, 128
dt = 1e-3
sim_name = "test2023000000"
restart = 0
steps = 25000
save_freq = 500
save_max = 20
print_freq = 200
wall_time = 23*60*60
save_dir = '.'


def main():
    # Define non-dimensional parameters
    Re = 1 / nu
    SS = L_T / (c_p * T_B)
    MM = (em * C_B * C_0) / T_B
    AA = epsilon * (SS * Re * Pr) * (5 / 6)
    GG = epsilon
    eta = 1e-1 * Re * (beta * epsilon) ** 2  # not "optimal"

    rho0 = sw.dens0(s=C_B * C_0, t=T_B)
    # rho0 = sw.dens0(s=5, t=20)

    def buoyancy_func(T, C):
        return -9.8 * 100 * (sw.dens0(C_B * C['g'], T_B * T['g']) - rho0) / rho0

    # Domain
    xbasis = de.Fourier('x', Nx, interval=(0, LL), dealias=3/2)
    zbasis = de.SinCos('z', Nz, interval=(0, HH), dealias=3/2)
    domain = de.Domain([xbasis, zbasis], grid_dtype=np.float64)
    x, z = domain.grids(domain.dealias)
    xx, zz = x + 0*z, 0*x + z
    kx, kz = domain.elements(0), domain.elements(1)
    # Wall penalty boundary
    wall = domain.new_field()
    wall.set_scales(domain.dealias)
    wall.meta['z']['parity'] = 1
    # wall['g'] = sigmoid(-(x-0.02*L), a=2*epsilon) + sigmoid(x-.98*L, a=2*epsilon)
    wall['g'] = 0  # no wall
    wall['c'] *= np.exp(-kx**2/5e6)  # spectral smoothing

    buoyancy = ParityFunction(domain, layout='g', func=buoyancy_func,)

    # Buoyancy multiplier for parity constraints
    par = domain.new_field()
    par.set_scales(domain.dealias)
    par.meta['z']['parity'] = -1
    par['g'] = np.tanh(-(z-HH)/.05)*np.tanh(z/.05)
    par['c'] *= np.exp(-kx**2/5e6)  # spectral smoothing

    # Mathematical problem
    melting = de.IVP(domain, variables=['u', 'w', 'p', 'T', 'C', 'f', 'ft'])
    melting.meta['u', 'p', 'T', 'C', 'f', 'ft']['z']['parity'] = 1
    melting.meta['w']['z']['parity'] = -1

    params = [Nx, Nz, SS, AA, GG, MM, delta, epsilon, Re, Pr,
              Sc, eta, hh, wall, par, UU, LL, HH, buoyancy, rho0]
    param_names = ['Nx', 'Nz', 'S', 'A', 'G', 'M', 'delta', 'epsilon', 'Re', 'Pr',
                   'Sc', 'eta', 'h', 'wall', 'par', 'U', 'L', 'H', 'buoyancy', 'rho0']
    for param, name in zip(params, param_names):
        melting.parameters[name] = param

    melting.substitutions['q'] = 'dz(u) - dx(w)'
    # melting.substitutions['ft'] = (
    #     '(epsilon/A)*((G/epsilon)*(dx(dx(f)) + dz(dz(f)))' +
    #     '- (1/epsilon**2)*f*(1-f)*((G/epsilon)*(1-2*f) + (T+M*C)))')

    print("C_B=", C_B, "C_0=", C_0, "T_B=", T_B)

    # Equations
    melting.add_equation("dx(u) + dz(w) = 0", condition='(nx != 0) or (nz != 0)')
    melting.add_equation("p = 0", condition='(nx == 0) and (nz == 0)')
    melting.add_equation(
        "dt(u) + dx(p) - (1/Re)*dz(q) = - w*q - (f/eta)*u - (wall/eta)*(u-U)")
    melting.add_equation(
        "dt(w) + dz(p) + (1/Re)*dx(q) = " +
        "u*q - (f/eta)*w - (wall/eta)*w + par*buoyancy")
    melting.add_equation(
        "dt(T) - (1/(Pr*Re))*(dx(dx(T)) + dz(dz(T))) - S*dt(f) = " +
        " - (1-f)*(u*dx(T) + w*dz(T)) + T*(u*dx(f) + w*dz(f)) - (wall/eta)*(T-1)")
    melting.add_equation(
        "dt(C) - (1/(Sc*Re))*(dx(dx(C)) + dz(dz(C)))   = " +
        "- u*dx(C) - w*dz(C) + " +
        " (C*ft - (dx(C)*dx(f) + dz(C)*dz(f))/(Sc*Re))/(1-f+delta) " +
        "- (wall/eta)*(C-1)")
    melting.add_equation(
        "(A/epsilon)*dt(f) - (G/epsilon)*(dx(dx(f)) + dz(dz(f))) = " +
        " - (1/epsilon**2)*f*(1-f)*((G/epsilon)*(1-2*f) + (T+M*C))")
    melting.add_equation("ft - dt(f) = 0")

    # Build timestepper and solver
    # CNAB bad, SBDF good, RK only good for BC on acceleration
    ts = de.timesteppers.RK222
    solver = melting.build_solver(ts)

    # Initial conditions
    u, w, p, T, C, f, ft = variables = [
        solver.state[field] for field in melting.variables]
    for field in variables:
        field.set_scales(domain.dealias)
    buoyancy.original_args = buoyancy.args = [T, C]
    if restart == 0:
        u['g'] = u_0
        w['g'] = w_0
        p['g'] = p_0
        f['g'] = (sigmoid(z-(HH-hh), a=2*epsilon) *
                  sigmoid(x-(LL-ll)/2, a=2*epsilon) *
                  sigmoid(-(x-(LL+ll)/2), a=2*epsilon))
        # print(f['g'])
        T['g'] = 1-f['g']
        C['g'] = C_0
        # C['g'] = well(z-(HH-hh), x-(LL-ll)/2, -(x-(LL+ll)/2))
        # C['g'] = abs(1-f['g']) #TZ 20240122: uncomment this to make salt only exist in liquid
        p['g'] = p_0  #pressure fluctuation

        # # NG 230310: unknown module file_tools; commenting out for now, won't
        # # try to restart
        # import file_tools as flts
        # if rank == 0:
        #     flts.save_domain('domain-{}.h5'.format(sim_name), domain)

    else:
        raise NameError("Restart function doesn't work because of file_tools")
        from glob import glob
        import re
        pattern = 'data-{n}-{r:0>2d}'.format(n=sim_name, r=restart-1)
        save_files = glob('{p}/{p}_s*.h5'.format(p=pattern))
        nums = [int(re.search('.*_s(\d+).h5', f).group(1)) for f in save_files]
        last = np.argmax(nums)
        write, _ = solver.load_state(save_files[last], -1)

    # Save configurations
    solver.stop_iteration = steps
    solver.stop_wall_time = wall_time
    solver.stop_sim_time = np.inf

    # CFL # Breaks quickly
    # CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=2, safety=.5,
    #                      max_change=1.5, min_change=0.5, max_dt=η/2,
    # threshold=0.1)
    # CFL.add_velocities(('u', 'w'))

    # # Flow properties # Doesnt work
    # flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
    # flow.add_property("sqrt(u*u + w*w)", name='Re')

    # Save state variables
    analysis = solver.evaluator.add_file_handler(
        join(save_dir, 'data-{}-{:0>2d}'.format(sim_name, restart)),
        iter=save_freq, max_writes=save_max, mode='overwrite')
    for task in melting.variables:
        analysis.add_task(task)
    analysis.add_task("integ(T - S*f,'x','z')", name='energy')
    analysis.add_task("integ((1-f)*C,'x','z')", name='salt')
    analysis.add_task("integ(f, 'x', 'z')", name='volume')
    analysis.add_task("integ(0.5 * (u**2 + w**2) * rho0, 'x', 'z')", name='kinematic energy')
    analysis.add_task("q")
    analysis.add_task("buoyancy")

    # Save parameters
    parameters = solver.evaluator.add_file_handler(
        join(save_dir, 'parameters-{}-{:0>2d}'.format(sim_name, restart)),
        iter=np.inf, max_writes=100, mode='overwrite')
    for task in melting.variables:
        parameters.add_task(task)
    for name in param_names:
        parameters.add_task(name)
    parameters.add_task("q")

    # Main Loop
    start_time = time.time()
    while solver.ok:
        if solver.iteration % print_freq == 0:
            max_speed = u['g'].max()
            logger.info('{:0>6d}, u max {:f}, dt {:.5f}, time {:.2f}'.format(
                solver.iteration, max_speed, dt, (time.time()-start_time)/60))
            if np.isnan(max_speed):
                sys.exit(1)
            if f.integrate('x', 'z')['g'][0, 0] < 1e-3:
                break

        solver.step(dt)
    solver.step(dt)
