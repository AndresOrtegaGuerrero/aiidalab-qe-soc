from aiidalab_qe.common.panel import ResultPanel

from aiidalab_qe.common.bandpdoswidget import BandPdosWidget


class Result(ResultPanel):
    title = "SOC"
    workchain_labels = ["soc"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, identifier="soc", **kwargs)

    def _update_view(self):
        # Check if the workchain has the outputs
        try:
            pdos_node = self.node.outputs.soc.pdos
        except AttributeError:
            pdos_node = None

        try:
            bands_node = self.node.outputs.soc.bands
        except AttributeError:
            bands_node = None

        bands_dos_widget = BandPdosWidget(bands=bands_node, pdos=pdos_node)

        self.children = [bands_dos_widget]
