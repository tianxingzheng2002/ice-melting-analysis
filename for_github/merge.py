from dedalus.tools import post
import pathlib

sim_name = "data-test2023111201-00"


def main():
    post.merge_process_files(sim_name, cleanup=True)
    set_paths = list(pathlib.Path(sim_name).glob(sim_name+"_s*.h5")) # "data-test2023111201-00_s*.h5"
    post.merge_sets(sim_name+"/"+sim_name+".h5", set_paths, cleanup=True)
    # "data-test2023111201-00/data-test2023111201-00.h5", set_paths, cleanup=True

