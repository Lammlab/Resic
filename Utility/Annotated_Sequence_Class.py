

from abc import ABC, abstractmethod
from collections.abc import MutableMapping

class SamLineTags(MutableMapping):
    def __init__(self, *strings):
        tags = [tag.split(':') for tag in strings[0]]
        for tag, string in zip(tags, strings[0]):
            if len(tag) is not 3:
                raise ValueError("%s is not a legal sam line tag" % string)
        self.dict = {tag_key: (tag_type, tag_val) for (tag_key, tag_type, tag_val) in tags}

    def __getitem__(self, item):
        return self.dict[item][1]

    def __setitem__(self, key, value):
        tag_type = 'i' if type(value) is int else 'Z'
        self.dict[key] = (tag_type, value)

    def __delitem__(self, key):
        del self.dict[key]

    def __iter__(self):
        for key in self.dict.keys():
            yield key

    def __len__(self):
        return len(self.dict)

    def __repr__(self):
        str_tags = [':'.join([key, type, val]) for key, (type, val) in self.dict.items()]
        return '\t'.join(str_tags)

class Annotated_Sequence(ABC):

    def __init__(self, tags):
        self.tags = SamLineTags(tags)

    @property
    @abstractmethod
    # (reference string,start,end) where start and end are relative to the 1st position in reference
    def position(self):...

    @position.setter
    @abstractmethod
    def position(self):...



    def annotate(self,string):
        # if we havent set an annotation field before
        if not ('annotate_count' in self.__dict__):
            self.annotate_count = sum('annotate' in tag for tag in self.tags)

        new_count = self.annotate_count + 1
        self.tags['annotate_%s' % str(new_count)] = string

        self.annotate_count = new_count
        return

    @abstractmethod
    def flip_strand(self):...
