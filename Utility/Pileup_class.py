from Utility.generators_utilities import getLineFromChunk
import re
from Utility.Annotated_Sequence_Class import Annotated_Sequence
from copy import deepcopy


class Pileup_line(Annotated_Sequence):
    fieldNames = {
        "reference_id": 0,
        "_gene_pos": 1,
        "reference": 2,
        "_base_count": 3,
        "reads_string": 4,
        "quality_string": 5,
    }

    def __init__(self, line , only_tags = False):
        if not only_tags:
            # parsing string
            line = line.strip('\n')
            fields = re.split('[\t\n]', line)
            dictionary = {key: fields[index] for (key, index) in self.fieldNames.items()}

            self.__dict__ = dictionary

            # create tags that act like dict
            raw_tags = fields[6:]
            if raw_tags == [''] and len(raw_tags) == 1:
                raw_tags = []
        else :
            line = line.strip('\n')
            raw_tags = re.split('[\t\n]', line)
            if raw_tags == ['']:
                raw_tags = []
        # We use a cheat here, for further metadata we use the sam format's annotation
        super().__init__(raw_tags)

    @property
    def gene_pos(self):
        return int(self._gene_pos)

    @gene_pos.setter
    def gene_pos(self, pos):
        self._gene_pos = str(pos)

    @property
    def base_count(self):
        return int(self._base_count)

    @base_count.setter
    def base_count(self, count):
        self._base_count = str(count)

    @property
    def reference_id(self):
        return self.__dict__["reference_id"]

    def __repr__(self):
        # get all fields in a tab seperated strings
        fields = [self.__dict__[i] for i in self.fieldNames.keys()]
        if str(self.tags) != "":
            fields.append(str(self.tags))
        return '\t'.join(fields)

    def line_to_csv_with_short_tags(self, list_tags):
        fields = [self.__dict__[i] for i in self.fieldNames.keys()]
        if str(self.tags) != "":
            tags_dict = dict()
            # preparing a tags dict
            for key,(type,value) in self.tags.dict.items():
                tags_dict[key] = value
            # adding the tags requested by the list
            for item in list_tags:
                try:
                    fields.append(tags_dict[item])
                except:
                    fields.append("")
                tags_dict.pop(item, None)
            # adding the misc tags
            if len(tags_dict.values()) == 0: # if we dont have misc tags , we add empty at the end
                fields.append("")
            else :
                for item in tags_dict.values():
                    fields.append(item)

        return '\t'.join(fields)

    def annotate(self, string):
        # We do the same thing for both sam and pileup now
        super().annotate(string)

    def ref(self):
        return self.reference_id

    def pos(self):
        return int(self.gene_pos), int(self.gene_pos)

    def position(self):
        return self.reference_id, int(self.gene_pos), int(self.gene_pos)

    def reads(self):
        return self.reads_string

    def flip_strand(self):
        curr_string = self.fields["reads_string"]
        new_string = ""
        for char in curr_string:
            if (char == '.'):
                new_string += ','
            elif (char == ','):
                new_string += '.'
            elif (char >= 'a' and char <= 'z'):
                new_string += char.upper()
            elif (char >= 'A' and char <= 'Z'):
                new_string += char.lower()
            else:
                pass
        self.fields["reads_string"] = new_string

    def is_with_any_change(self): 

        #return : true if in the reads string there is one or more of a c g t A C G T , else false
        #all the possible IUPAC nucleotide code
        nuc_changes = ["a","c","g","t","r","y","s","w","k","m","b","d","h","v","n"]
        all_nuc_changes = [letter.upper() for letter in nuc_changes] + nuc_changes

        if any(change in self.reads() for change in all_nuc_changes):
            return True
        else : 
            return False


    @property
    def clean_base_string(self):
        """

        :return: the base_string of the pileup lines without indels and ^/$ marks
        """

        # TODO currently assumes no idels. fix this.

        raw_string_iter = self.reads_string.__iter__()
        clean_str = [char for char in raw_string_iter if not char in ('^', '$')]
        clean_str = ''.join(clean_str)
        return clean_str

    @classmethod
    # merge pileup lines
    def merge_pileup_lines(cls, lines):
        # checks that all reference bases are the same
        ids = {line.reference_id for line in lines}
        if len(ids) is not 1:
            raise ValueError("Lines have different references strings\n" + str(lines) + "\n")
        positions = {line.gene_pos for line in lines}
        if len(positions) is not 1:
            raise ValueError("Lines have different positions\n" + str(lines) + "\n")
        ref_bases = {line.reference for line in lines}
        if len(ref_bases) is not 1:
            raise ValueError("Lines in same positions have different reference nucleotide\n" + str(lines) + "\n")

        # merges reads and quality
        merged_reads = ''.join([line.reads_string for line in lines])
        merged_qual = ''.join([line.quality_string for line in lines])
        total_count = sum([line.base_count for line in lines])

        # creates merge pileup
        merged_pileup = deepcopy(lines[0])
        merged_pileup.reads_string = merged_reads
        merged_pileup.quality_string = merged_qual
        merged_pileup.base_count = total_count
        return merged_pileup
