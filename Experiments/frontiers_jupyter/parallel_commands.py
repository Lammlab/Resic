import logging
from multiprocessing.pool import ThreadPool
import traceback


def apply(*command):
    f = command[0]
    args = command[1:]
    try:
        logging.info(f" started process running {f}({*args,})")
        res = f(*args)
        logging.info(f" finished process running {f}({*args,})")

    except Exception as e:
        logging.info(f" Caught exception in worker process running {f}({*args,})")

        traceback.print_exc()
        raise e

    return res


# recives a list of Command class items
def parallel_commands(list_commands, parallel_limit=2, disable_parallel=True):
    commands = [tuple(com) for com in list_commands]

    if disable_parallel:
        for command in commands:
            apply(*command)

    else:
        pool = ThreadPool(parallel_limit)
        pool.starmap(apply, commands)
        pool.close()
        pool.join()

    return
