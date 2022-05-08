from Utility.temp_files import new_temp_file
from Utility.strip_last_newline import strip_last_NL
import os
import operator
import itertools
import logging


def fileGetChunks(file1, size=10000):
    """
    :param file1: the open file
    :param size: the number of lines we want to read each time
    :return: list of lines in the number requested
    """

    chunk = file1.readlines(size)
    while len(chunk) > 0:
        yield chunk
        chunk = file1.readlines(size)
    return


def getLineFromChunk(file1, size=10000, remove_header=False):
    """
    a generator of list of lines
    :param file1:  an open file
    :param size: the number of lines we want to get each time
    :param remove_header: is there a header line
    :return: list of lines that was requested each time, if header==True, skip it
    """

    if remove_header:
        # ignores first line
        file1.readline()

    gen = fileGetChunks(file1, size)

    for chunk in gen:
        for line in chunk:
            yield line
    return


def getKLinesFromChunk(file, size=40000, k=4, remove_header=False):
    """
    generator of k lines at a time
    :param file: the opened file
    :param k: the value of k
    :return: yields an array of 4 lines from the file at a timep
    """
    gen = getLineFromChunk(file, size=size, remove_header=remove_header)

    k_batch = []
    for line in gen:
        # got k lines, yield them
        # TODO use group_generator
        k_batch.append(line)
        if len(k_batch) == k:
            yield k_batch
            k_batch = []

    return


def key_sorted_gen(key, file=None, gen=None, *args, **kwargs):
    is_gen = True  # incdicate wheather gen param was received
    if not gen:
        is_gen = False
        gen = getLineFromChunk(file, *args, **kwargs)

    sorted_file = new_temp_file()
    with open(sorted_file, "w+") as fp:
        for obj in sorted(gen, key=key):
            if not is_gen:
                fp.write(str(obj))
            else:
                fp.write(str(obj) + '\n')

    strip_last_NL(sorted_file)

    with open(sorted_file, "r") as sorted_fp:
        for line in getLineFromChunk(sorted_fp):
            yield line

    os.remove(sorted_file)


def skip_broken_lines_factory(class_name):
    def skip_broken_lines(line):
        try:
            class_name(line)
            return False
        except:
            return True

    return skip_broken_lines


def class_generator(class_type, skip_condition=lambda x: False, file=None, gen=None, number_of_lines=1, class_kwargs={},
                    **kwargs):
    # file should be an open file object
    if not gen:
        if number_of_lines == 1:
            gen = getLineFromChunk(file, **kwargs)
        else:
            gen = getKLinesFromChunk(file, k=number_of_lines)

    for line_num, line in enumerate(gen):
        try:
            if not skip_condition(line):
                yield class_type(line, **class_kwargs)
        except Exception as e:
            logging.exception(
                f"In file {file}, In line {number_of_lines * line_num}, Exception was raised when loading class {class_type}")
            continue

    return


def generates4Lines(file, size=10000):
    yield from getKLinesFromChunk(file, size, k=4)


def window_generator(gate_gen1, gen2, functor_gen1, functor_gen2, closed_interval_first_side=False,
                     closed_interval_second_side=False):
    """
    @param gate_gen1: a generator of gate items (for example primaries)- we look at the items of that gen as a
    border items that we find gen2 items between them
    @param gen2: a generator of items (for example secondaries) - those will be the items we try to match to be in
    between 2 gate items
    @param functor_gen1: gets a gen1 item and returns comparision key (to compare to gen2 comparision key)
    @param functor_gen2: gets a gen2 item and returns comparision key (to compare to gen1 comparision key)
    @param closed_interval_first_side: True if we allow the between gen2 items an equality to the first gate side
    - by default it is not allowed
    @param closed_interval_second_side: True if we allow the between gen2 items an equality to the second gate side
    - by default it is not allowed
    @return: (gate_item1, middle_list, gate_item2) - where the middle list are the items from gen2 there are located
    between the two borders according to the booleans equality allowed by sides
    """

    gate_gen1 = gate1_generator(gate_gen1, functor_gen1)  # run over the received gen with a gen that returns a list
    gate_item1 = list()
    try:
        gate_item2 = next(gate_gen1)
    except StopIteration:
        gate_item2 = None
    try:
        middle_item = next(gen2)
    except StopIteration:
        middle_item = None

    middle_list = list()
    not_finished = True

    gt = operator.ge if closed_interval_first_side else operator.gt  # ">=" if closed_interval_first_side else ">"
    lt = operator.le if closed_interval_second_side else operator.lt  # "<=" if closed_interval_second_side else "<"

    while not_finished:
        if middle_item is None:
            not_finished = False  # exit the loop
            # yield the list we have this far - in between the two curr gate items
            gate_item1, gate_item2, middle_list = yield from yield_window_and_promote(gate_gen1, gate_item1, gate_item2,
                                                                                      middle_list)
            continue

        if gate_item1 is not None and len(gate_item1) != 0:
            if gt(functor_gen2(middle_item), functor_gen1(gate_item1[0])):
                if (len(gate_item2) != 0 and lt(functor_gen2(middle_item), functor_gen1(gate_item2[0]))) \
                        or (len(gate_item2) == 0):
                    # it is in between the two gates or no second border
                    middle_item = add_to_window_and_promote_middle(gen2, middle_item, middle_list)

                elif (len(gate_item2) != 0) and not (lt(functor_gen2(middle_item), functor_gen1(gate_item2[0]))):
                    # when the middle is bigger then the second gate
                    gate_item1, gate_item2, middle_list = yield from yield_window_and_promote(gate_gen1, gate_item1,
                                                                                              gate_item2, middle_list)

            else:  # when the middle is smaller then the first gate
                gate_item1, gate_item2, middle_list = yield from yield_window_and_promote(gate_gen1, gate_item1,
                                                                                          gate_item2, middle_list)

        else:  # when we are at the beginning and have no first boarder
            if gate_item2 is not None:
                if len(gate_item2) != 0 and lt(functor_gen2(middle_item), functor_gen1(gate_item2[0])):
                    middle_item = add_to_window_and_promote_middle(gen2, middle_item, middle_list)
                elif len(gate_item2) != 0:
                    gate_item1, gate_item2, middle_list = yield from yield_window_and_promote(gate_gen1, gate_item1,
                                                                                              gate_item2, middle_list)
            else:  # both gates are None
                return

    # handle the case where there are no more middle, but there are gates

    while gate_item1 is not None and len(gate_item1) != 0:
        yield (gate_item1, middle_list, gate_item2)
        if len(gate_item2) != 0:
            gate_item1, gate_item2, middle_list = promote_window(gate_gen1, gate_item2)
        else:
            return


def yield_window_and_promote(gate_gen1, gate_item1, gate_item2, middle_list):
    """
    yields for the calling func , checks validity of window promoting and promotes the window
    """
    yield (gate_item1, middle_list, gate_item2)
    if len(gate_item2) != 0:
        gate_item1, gate_item2, middle_list = promote_window(gate_gen1, gate_item2)
    return gate_item1, gate_item2, middle_list


def gate1_generator(generator, functor_gen1):
    """
    :param generator: a generator of the gate1 objects
    :param functor_gen1: functor for comparisson of the generator returned items
    :return: yield a list of lines - if there were dup lines - all dup will appear in the list, if no dup -
    one line in the yielded list
    """
    new_item = next(generator)
    r_list = list()
    while new_item:
        r_list.append(new_item)
        try:
            next_item = next(generator)
        except StopIteration:
            yield r_list
            raise
        while operator.eq(functor_gen1(next_item), functor_gen1(new_item)):
            r_list.append(next_item)
            temp = next_item
            try:
                next_item = next(generator)
            except StopIteration:
                yield r_list
                raise
            new_item = temp

        yield r_list
        r_list = list()
        new_item = next_item


def promote_window(gate_gen1, gate_item2):
    """

    :param gate_gen1: the generator of the gates objects
    :param gate_item2: the current gate 2 object
    :return: promoted gate 1 , new middle list , promoted gate 2
    """
    try:
        temp_gate = next(gate_gen1)
    except RuntimeError:
        return gate_item2, [], []

    gate_item1 = gate_item2  # runs over the list
    gate_item2 = temp_gate

    middle_list = list()
    return gate_item1, gate_item2, middle_list


def add_to_window_and_promote_middle(gen2, middle_item, middle_list):
    """

    :param gen2: generator for middle objects
    :param middle_item: the middle item we want to add to the curr window
    :param middle_list: the list we want to add the middle item to
    :return: the updated middle list
    """
    middle_list.append(middle_item)
    try:
        middle_item = next(gen2)
    except StopIteration:
        return None
    return middle_item


def split_on_condition(seq, condition):
    """
    splits gen to two gens based on whether the condition evaluated the items to true or false
    :param seq: iterable
    :parma condition: boolean function
    "returns l1,l2: where l1 is the generator of true element accoridng to condition and 
    l2 is the false elements gen
    """
    l1, l2 = itertools.tee((condition(item), item) for item in seq)
    return (i for p, i in l1 if p), (i for p, i in l2 if not p)


def zip_equal(*iterables):
    """
    zip iterables but raise exception if they are not same length
    :param iterables:
    :return:
    """
    sentinel = object()
    for combo in itertools.zip_longest(*iterables, fillvalue=sentinel):
        if sentinel in combo:
            raise ValueError('Iterables have different lengths')
        yield combo


def group_generator(iterable, n):
    "group_generator(3, 'ABCDEFGHI') --> ABC DEF GHI"
    args = [iter(iterable)] * n
    return zip_equal(*args)
