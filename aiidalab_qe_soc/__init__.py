from .setting import Setting
from .workchain import workchain_and_builder
from .result import Result
from aiidalab_qe.common.panel import OutlinePanel
from aiidalab_widgets_base import ComputationalResourcesWidget
from aiidalab_qe.plugins.pdos import dos_code, projwfc_code
class SOCOutline(OutlinePanel):
    title = "Bands/Pdos w SOC"


soc ={
"outline": SOCOutline,
"code": {"dos": dos_code, "projwfc": projwfc_code},
"setting": Setting,
"result": Result,
"workchain": workchain_and_builder,
}