import re
import tempfile
from collections.abc import MutableMapping
from Utility.generators_utilities import getLineFromChunk
from Utility.temp_files import *


class SamLineTags(MutableMapping):
    def __init__(self, *strings):
        tags = [tag.split(':') for tag in strings]
        for tag, string in zip(tags, strings):
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


class SamLine:
    # TODO - document this class - need to write how we can get the class fields
    fieldNames = {
        "seq_id": 0,
        "flags": 1,
        "reference_id": 2,
        "_gene_pos": 3,
        "mapq": 4,
        "cigar": 5,
        "rnext": 6,
        "pnext": 7,
        "tlen": 8,
        "sequence": 9,
        "quality": 10,
    }

    def __init__(self, line):
        line = line.strip('\n')
        fields = re.split('\t|\n', line)
        # dictionary = {key: fields[index] for (key, index) in self.fieldNames.items()}

        # self.__dict__ = dictionary

        self.seq_id = fields[0]
        self._flags = fields[1]  # int
        self.reference_id = fields[2]
        self._gene_pos = fields[3]  # int
        self._mapq = fields[4]  # int
        self.cigar = fields[5]
        self.rnext = fields[6]
        self._pnext = fields[7]  # int
        self._tlen = fields[8]  # int
        self.sequence = fields[9]
        self.quality = fields[10]

        # create tags that act like dict
        raw_tags = fields[11:]
        self.tags = SamLineTags(*raw_tags)

    #

    def __repr__(self):
        # fields = [self.__dict__[i] for i in self.fieldNames.keys()]
        fields = [self.seq_id, self._flags, self.reference_id, self._gene_pos, self._mapq, self.cigar, self.rnext,
                  self._pnext, self._tlen, self.sequence, self.quality]
        fields.append(str(self.tags))
        return '\t'.join(fields)

    def is_alligned(self):
        # true if the reference is not *
        return self.reference_id is not '*'

    def strand(self):
        """ returns the sense of the sequence """
        # logic: if bit 5 (16) is on, this will become -1
        #        else it will be one
        if int(self.flags) & 16 == 16:
            return -1
        else:
            return 1

    def annotate(self, string):
        # if we havnt set an annotation field before
        if not ('annotate_count' in self.__dict__):
            self.annotate_count = sum('annotate' in tag for tag in self.tags)

        new_count = self.annotate_count + 1
        self.tags['annotate_%s' % str(new_count)] = string

        self.annotate_count = new_count
        return

    # use the position method to get reference_id, sequence start position, and sequence end position
    def position(self):
        return self.reference_id, int(self._gene_pos), int(self._gene_pos) + len(self.sequence)

    def sequence_id(self):
        return self.seq_id

    @classmethod
    def is_header(cls, line):
        return line[0] is '@'

    # TODO move this to Annotated_sequence
    @property
    def gene_pos(self):
        return int(self._gene_pos)

    @gene_pos.setter
    def gene_pos(self, pos):
        self._gene_pos = str(pos)

    @property
    def pnext(self):
        return int(self._pnext)

    @pnext.setter
    def pnext(self, next_pos):
        self._pnext = str(next_pos)

    @property
    def flags(self):
        return int(self._flags)

    @flags.setter
    def flags(self, flag):
        self._flags = str(flag)

    @property
    def mapq(self):
        return int(self._mapq)

    @mapq.setter
    def mapq(self, quality):
        self._mapq = str(quality)

    @property
    def tlen(self):
        return int(self._tlen)

    @tlen.setter
    def tlen(self, length):
        self._tlen = str(length)


"""
:param sam_filename: filename of sam file
:return: 3 output files: 1. header lines file called ${header_filename}
						 2. aligned lines file called ${aligned_filename}
						 3. misaligned lines file called ${misaligned_filename}
"""


def split_sam(sam_filename):
    header_file = open(new_temp_file(),
                       "w+")  # tempfile.NamedTemporaryFile(prefix="sam_splitting", dir="/tmp", mode='w+', delete=False)
    aligned_file = open(new_temp_file(),
                        "w+")  # tempfile.NamedTemporaryFile(prefix="sam_splitting", dir="/tmp", mode='w+', delete=False)
    misaligned_file = open(new_temp_file(),
                           "w+")  # tempfile.NamedTemporaryFile(prefix="sam_splitting", dir="/tmp", mode='w+', delete=False)

    with open(sam_filename, "r") as file:
        gen = getLineFromChunk(file)
        for line in gen:
            if line[0] == "@":  # header sam lines
                header_file.write(line)
            elif re.split('\t|\n', line)[2] == "*":  # If the chr_name field is * than the line is misaligned
                misaligned_file.write(line)
            else:
                aligned_file.write(line)

    header_file.close()
    aligned_file.close()
    misaligned_file.close()

    return (header_file.name, aligned_file.name, misaligned_file.name)
