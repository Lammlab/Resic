import pandas as pd

def column_delta_df(df:pd.DataFrame)->pd.DataFrame:
    """
    returns a DataFrame where each colum j but the last takes the value of col[j]-col[j-1]
    :param df: a pandas df contatining only numeric values
    :return: a pandas df
    """
    delta=pd.DataFrame.copy(df,deep=True)

    delta.iloc[:,1:]=df.values[:,1:] - df.values[:,:-1]
    return  delta

def get_line_count(filename):
    """
    count number of lines in file.
    taken from
    https://stackoverflow.com/a/27518377
    :param filename: file name
    :return: number of lines in file
    """
    def _make_gen(reader):
        b = reader(1024 * 1024)
        while b:
            yield b
            b = reader(1024 * 1024)

    f = open(filename, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum(buf.count(b'\n') for buf in f_gen)
