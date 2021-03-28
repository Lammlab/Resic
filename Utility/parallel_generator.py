
def parallel_generator(generators, functors):
    """
    :param generators: A list of k generators (initialized) repersenting files, 
                                       file are sorted with respect to the functors.
    :param functors: A list of k functors (must be the same size as the generators),
                                       that return a comparable object (comaprable with < ).
    :output: Each yield returns the next equivalence class according to the sorting
                     in a format of a k-list with None if the generator didn't have an
                     object of that equivalence class. Each entry in the list is also
                     a list of the objects in that equivalence class that the 
                     corresponding generator produced.
                     For example output may look like this (where k = 3): [None,[object1], [object2, object3]]
    """

    # Work buffer, lines that aren't in output yet reside here
    #current_lines = [generator.__next__() if generator is not None else None for generator in generators]
    current_lines = []
    for generator in generators:
        try:
            line = next(generator)
        except StopIteration:
            line = None
        current_lines.append(line)

    if len(current_lines) == current_lines.count(None):
        return

    finished_files = 0
    number_of_files = len(generators)
    # Boolean array for which generator finished
    listOfFinished = [0 for index in range(len(generators))]


    # Actual logic
    while (finished_files < number_of_files) :  # As long as not all files ended
        current_intervals = [func(line) if line is not None else None for line, func in zip(
            current_lines, functors)]  # Extracting (str,int)

        # Sort the intervals according to the (default) sorting of the given
        # object
        current_intervals_sorted = sorted(
            [interval for interval in current_intervals if interval is not None])

        min_interval = current_intervals_sorted[0]  # Minimum interval

        # Get all the indexes (in the original list) of the lines with the same
        # minimal equivalence class
        equivalence_class_indexes = [index for index in range(
            len(current_intervals)) if (current_intervals[index] == min_interval)]

        # Number of generators that generated in current eq. class
        num_gen_curr_eq_class = len(equivalence_class_indexes)

        # Init output list with Nones
        next_equivalence_class = [None] * len(generators)
        # Init lists for generators that outputed
        for i in equivalence_class_indexes:
            next_equivalence_class[i] = []

        # For each generator check if he generated an object(line) in the equivalence class,
        # If it did add the object to the output list and try to add next line
        # until he stops generating from this eq. class (else it stays None)
        while(num_gen_curr_eq_class > 0):
            for i in equivalence_class_indexes:
                next_equivalence_class[i].append(current_lines[i])
                try:
                    current_lines[i] = next(generators[i])
                    new_interval = functors[i](current_lines[i])
                except StopIteration:
                    current_lines[i] = None
                    equivalence_class_indexes.remove(i)
                    num_gen_curr_eq_class -= 1
                    break

                if(new_interval > min_interval):
                    equivalence_class_indexes.remove(i)
                    num_gen_curr_eq_class -= 1
                    break

        yield next_equivalence_class  # Yield the equivlanece class

        # Setting up lists for next iteration
        # For each generator that generated an outputted line, get his next
        # line
        for i in equivalence_class_indexes:
            try:
                current_lines[i] = next(generators[i])
            except StopIteration:
                current_lines[i] = None

        # Check if any generator finished yielding eveyrthing
        for index in range(len(current_lines)):
            if current_lines[index] == None:  # Meaning no more output
                if listOfFinished[index] == 0:  # If not already finished
                    finished_files += 1
                    listOfFinished[index] = 1

    return  # When all generators finished
