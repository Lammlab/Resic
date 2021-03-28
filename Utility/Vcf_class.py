from Utility.generators_utilities import getLineFromChunk
import re
from Utility.Annotated_Sequence_Class import Annotated_Sequence
from copy import deepcopy


class VcfClass(Annotated_Sequence):
    fieldNames = {
        "reference_id": 0,
        "_gene_pos": 1,
        "ID": 2,
        "REF": 3,
        "ALT": 4,
        "quality": 5,
        "FILTER": 6,
        "INFO": 7
    }


    def __init__(self, line , only_tags = False):
        if not only_tags:
            # parsing string
            line = line.strip('\n')
            fields = re.split('[\t\n]', line)
            dictionary = {key: fields[index] for (key, index) in self.fieldNames.items()}

            self.__dict__ = dictionary

            # create tags that act like dict
            raw_tags = fields[8:]
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

    def flip_strand(self):
        curr_string = self.fields["REF"]
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
        self.fields["REF"] = new_string

    def is_with_any_change(self): 
        pass


    @property
    def clean_base_string(self):
        pass

    @classmethod
    # merge pileup lines
    def merge_pileup_lines(cls, lines):
        pass