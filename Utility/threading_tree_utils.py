from contextlib import contextmanager
import multiprocessing
from typing import Dict
from typing import List
from collections import namedtuple
import copy

Fathers_list = List[str]
Specifications = namedtuple('Specifications', ['id', 'flags', 'pre_function', 'post_function'])
Dict_specifications = Dict[Specifications, Fathers_list]

""" this class handles the waiting of the child vertexes to thier fathers, 
handles the dict of filenames that contains the "output" file of each step for it's children ,
and in general handles the tree structure by the given Dict_specifications """


class ThreadingTree:

    def __init__(self, dict_spec: Dict_specifications):
        self.dict_tree = copy.deepcopy(dict_spec)

        if not self.is_not_circular_dependency(dict_spec):
            raise ValueError("The dictionary of specification must have a step without prior dependencies")

        # filling dict of events for each vertex of "tree"
        self.events_by_id_dict = {spec.id: multiprocessing.Event() for spec in dict_spec.keys()}
        self.filenames_dict = multiprocessing.Manager().dict({spec.id: None for spec in dict_spec.keys()})
        self.filenames_dict_lock = multiprocessing.Lock()
        self.fasta_indexing_lock = multiprocessing.Lock()
        self.secondary_fasta_indexing_lock = multiprocessing.Lock()

    def is_not_circular_dependency(self, dict_spec: Dict_specifications):
        """

        :param dict_spec:
        :return: True if there is a head to start the threading from
        """
        hist_dict = dict()
        for tuple in dict_spec.keys():

            id = tuple[0]
            if id not in hist_dict:
                hist_dict[id] = 1
            else:
                hist_dict[id] = hist_dict[id] + 1
        for id, number in hist_dict.items():
            if number == 1:
                return True

        return False

    @contextmanager
    def threading_context(self, id_step):
        """
        this function is used with a with fraze : "with ThreadingTree.threading_context(id): ... "
        when the with begin the part before the yield is executed, then the with scope content is being preformed and
        when the with scope is done the part at the finally is being executed.

        this function handles the waiting of the child vertex to it's fathers.
        :param id_step: the id of the vertex in the tree
        :return:
        """
        # here need to come the enter part - wait for father on the event and so on
        for spec, fathers in self.dict_tree.items():
            if spec.id == id_step:
                for father in fathers:
                    self.events_by_id_dict[father].wait()  # waiting on all fathers

        try:
            yield
        finally:  # here is exit
            self.events_by_id_dict[id_step].set()  # setting the curr event - because the action is done
            pass

    @contextmanager
    def locking_context(self):
        """
        locking for the sake that no two threads will try to do bowtie build at the same time and collapse
        """
        self.fasta_indexing_lock.acquire(True)
        try:
            yield
        finally:  # here is exit
            self.fasta_indexing_lock.release()
            pass

    @contextmanager
    def secondary_locking_context(self):
        """
        locking for the sake that no two threads will try to do bowtie build at the same time and collapse
        """
        self.secondary_fasta_indexing_lock.acquire(True)
        try:
            yield
        finally:  # here is exit
            self.secondary_fasta_indexing_lock.release()
            pass

    def set_filename_of_finished_alignment(self, id_step, name_of_file):
        """
        this sets the filename for the curr step , for the use of its child vertexes
        :param name_of_file: the name of the output file from the id_step vertex
        :return:
        """

        self.filenames_dict_lock.acquire(True)

        self.filenames_dict[id_step] = name_of_file

        self.filenames_dict_lock.release()

    def get_fathers_and_spec(self, step_id):
        """
        :param step_id:
        :return: fathers list and the spec (Specifications) of the curr step_id , like it appeared in the dict_spec
        that was originally given to the class
        """
        for spec, fathers in self.dict_tree.items():
            if spec.id == step_id:
                return copy.deepcopy(fathers), copy.deepcopy(spec)
        return [], None

    def get_previous_steps_file_names(self, id_step):

        filenames = []
        fathers, spec = self.get_fathers_and_spec(id_step)
        for father in fathers:
            filenames.append(copy.deepcopy(self.filenames_dict[father]))
        return filenames

    def get_step_filename(self, id_step):
        return copy.deepcopy(self.filenames_dict[id_step])

    def get_id_of_leafs(self):
        """
        :return: a list of the id's of the leafs children of the graph
        """
        hist_dict = dict()
        for tuple1, fathers in self.dict_tree.items():
            for father in fathers:
                if father not in hist_dict:
                    hist_dict[father] = 1
                else:
                    hist_dict[father] = hist_dict[father] + 1
        leafs_list = []
        for tuple1, fathers in self.dict_tree.items():
            if tuple1[0] not in hist_dict.keys():
                leafs_list.append(copy.deepcopy(tuple1[0]))

        return leafs_list

    def get_filenames_of_leafs(self):
        """
        :return: the filenames the leaf vertex of the tree set as their filenames results
        """
        res = []
        for id1 in self.get_id_of_leafs():
            res.append(copy.deepcopy(self.filenames_dict[id1]))
        return res
