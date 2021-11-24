"""This file will generate an appropriate hitpick plan for sgRNA synthesis for the desired guides.
At this time all reactions are done separately.  Pooling for generation of digest reaction conditions happens post-synthesis
to enable easier troubleshooting of the guides."""

class SynRNAPlanner:

    def __init__(self):
        self.template_plate_map = dict()  # dictionary, component name as key, well id as value


    def match_ids_to_wells(self):
        for template in self.template_plate_map:

