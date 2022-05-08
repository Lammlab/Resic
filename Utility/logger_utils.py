import os
import inspect


def get_Pipes_dir():
    # Nasty Hack
    cur_frame = inspect.getfile(inspect.currentframe())
    abs_path = os.path.abspath(cur_frame)
    for i in range(2):  # level of file in realtion to th emain pipes dir
        abs_path = os.path.dirname(abs_path)
    return abs_path
