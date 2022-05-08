import tempfile
import os
import socket

socket_debug = False


def new_temp_file(dir='/tmp/', *args, **kwargs):
    # TODO fix to check tmp dir with os
    if (dir == '/tmp/' or dir == '/tmp') and socket.gethostname() == 'akalton':
        dir = '/Bigdata/tmp/'

    if socket_debug:
        with open("/tmp/socket_log.log", 'a') as logfp:
            logfp.write("socket is %s, dir is %s \n" % (socket.gethostname(), dir))

    if not os.path.isdir(dir):
        os.makedirs(dir)
    os_handle, name = tempfile.mkstemp(dir=dir, *args, **kwargs)
    os.close(os_handle)
    return name
