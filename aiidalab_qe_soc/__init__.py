from .setting import Setting
from .workchain import workchain_and_builder
from .result import Result
from aiidalab_qe.common.panel import OutlinePanel
from aiidalab_widgets_base import ComputationalResourcesWidget

class SOCOutline(OutlinePanel):
    title = "Bands/Pdos w SOC"


dos_code = ComputationalResourcesWidget(
    description="dos.x soc",
    default_calc_job_plugin="quantumespresso.dos",
)

projwfc_code = ComputationalResourcesWidget(
    description="projwfc.x soc",
    default_calc_job_plugin="quantumespresso.projwfc",
)

soc ={
"outline": SOCOutline,
"code": {"soc_dos": dos_code, "soc_projwfc": projwfc_code},
"setting": Setting,
"result": Result,
"workchain": workchain_and_builder,
}