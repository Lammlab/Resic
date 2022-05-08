import os


def strip_last_NL(in_file):
    if os.path.getsize(in_file) <= 2:
        return  # file too small
    with open(in_file, "a+b") as fp:
        # get last 2 bytes
        fp.seek(-2, 2)
        bytes = fp.read(2)
        last_chars = [chr(b) for b in bytes]

        # 2 newline chars at end
        if set(last_chars) == {'\n', '\r'}:
            truncate = 2
        # 1 newline char at end
        elif last_chars[1] == '\n' or last_chars[1] == '\r':
            truncate = 1
        # no newline chars at end
        else:
            truncate = 0
        # note we are at end of file here so tell gives as size of file
        fp.truncate(fp.tell() - truncate)
