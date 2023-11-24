
from aiidalab_qe.common.panel import ResultPanel

from .widget import BandDosPlotsWidget

class Result(ResultPanel):

    title = "SOC"
    workchain_labels = ["soc"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, identifier="soc", **kwargs)

    def _update_view(self):

        #Check if the workchain has the outputs

        bands_dos_widget = BandDosPlotsWidget(node=self.node.outputs.soc)

        self.children = [bands_dos_widget]

        