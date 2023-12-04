"""Panel for SOC plugin.
"""
import ipywidgets as ipw
import traitlets as tl
from aiida.orm import StructureData
from aiidalab_qe.common.panel import Panel
from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from IPython.display import clear_output, display

class Setting(Panel):
    
    title = "SOC Settings"
    identifier = "soc"

    protocol = tl.Unicode(allow_none=True)
    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Settings</h4></div>"""
        )
        self.calc_options = ipw.Dropdown(
            description="Calculation:",
            options=[
                ("Bands", "bands"),
                ("Pdos", "pdos"),
                ("Bands/Pdos", "bands_pdos"),
            ],
            value="bands",
        )
        self.kpath_2d = ipw.Dropdown(
            description="Lattice:",
            options=[
                ("Hexagonal", "hexagonal"),
                ("Square", "square"),
                ("Rectangular", "rectangular"),
                ("Centered Rectangular", "centered_rectangular"),
                ("Oblique", "oblique"),
            ],
            value="hexagonal",
        )

        self.lattice_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            If your system has periodicity xy. Please select one of the five 2D Bravais lattices corresponding to your system
            </div>"""
        )
        self.pseudos_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Spin-Orbit Coupling calculations utilize Full Relativistic Pseudopotentials from the PseudoDojo library.
            </div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 150%; padding: 0px 10px 5px;">
            Choose the computation type:
            <ul style="margin-top: 5px; margin-bottom: 5px; padding-left: 0; display: flex; flex-wrap: wrap;">
                <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; Band Structure.</li>
                <li style="margin-right: 10px; list-style-type: none; display: inline-block;">&#8226; Projected Density of States.</li>
                <li style="list-style-type: none; display: inline-block;">&#8226; Band Structure and Projected Density of States.</li>
            </ul>
            </div>"""
        )

        self.cutoffs_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Please set the cutoffs for the calculation. The default values are recommended for most systems.
            </div>"""
        )

        self.mesh_grid = ipw.HTML()

        self.lattice_2d_out = ipw.Output()
        self.lattice_2d_box = ipw.VBox([self.lattice_help, self.kpath_2d])

        self.soc_ecutwfc = ipw.BoundedFloatText(
            value=90.0,
            min=1,
            max=2000,
            step=1,
            description="Energy ecutwfc [Ry]:",
            style={"description_width": "initial"},
        )
        self.soc_ecutrho = ipw.BoundedFloatText(
            value=360.0,
            min=1,
            max=4000,
            step=1,
            description="Energy ecutrho [Ry]:",
            style={"description_width": "initial"},
        )

        # nscf kpoints setting widget
        self.nscf_kpoints_distance = ipw.BoundedFloatText(
            min=0.001,
            step=0.01,
            value=0.1,
            description="NSCF K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )

        self.nscf_kpoints_distance.observe(self._display_mesh, "value")

        self.children=[
                self.settings_title,
                self.pseudos_help,
                self.settings_help,
                self.calc_options,
                self.lattice_2d_out,
                self.cutoffs_help,
                ipw.HBox([self.soc_ecutwfc, self.soc_ecutrho,]),
                ipw.HBox([self.nscf_kpoints_distance, self.mesh_grid,]),
            ]
        super().__init__(**kwargs)

    @tl.observe("protocol")
    def _protocol_changed(self, change):
        """Input protocol changed, update the widget values."""
        self.nscf_kpoints_distance.value = PdosWorkChain.get_protocol_inputs(self.protocol)["nscf"]["kpoints_distance"]
        self._display_mesh()

    @tl.observe("input_structure")
    def _update_structure(self, _=None):
        self._display_mesh()

    @tl.observe("input_structure")
    def _update_settings(self, _=None):
        self._display_lattice_box()

    def _display_lattice_box(self):
        with self.lattice_2d_out:
            if self.input_structure:
                if self.input_structure.pbc == (True,True,False):
                    display(self.lattice_2d_box)
                else:
                    clear_output(wait=True)

    def _display_mesh(self, _=None):
        if self.input_structure is None:
            return
        mesh = create_kpoints_from_distance.process_class._func(
            self.input_structure,
            orm.Float(self.nscf_kpoints_distance.value),
            orm.Bool(True),
        )
        self.mesh_grid.value = "Mesh " + str(mesh.get_kpoints_mesh()[0])

    
    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "calc_options": self.calc_options.value,
            "soc_ecutwfc": self.soc_ecutwfc.value,
            "soc_ecutrho": self.soc_ecutrho.value,
            "kpath_2d": self.kpath_2d.value,
            "nscf_kpoints_distance": self.nscf_kpoints_distance.value,
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.kpath_2d.value = input_dict.get("kpath_2d", "hexagonal")
        self.calc_options.value = input_dict.get("calc_options", "bands")
        self.soc_ecutwfc.value = input_dict.get("soc_ecutwfc", 90.0)
        self.soc_ecutrho.value = input_dict.get("soc_ecutrho", 360.0)
        self.nscf_kpoints_distance.value = input_dict.get("nscf_kpoints_distance", 0.1)

    def reset(self):
        """Reset the panel to its default values."""
        self.kpath_2d.value = "hexagonal"
        self.calc_options.value = "bands"
        self.soc_ecutwfc.value = 90.0
        self.soc_ecutrho.value = 360.0
        self.nscf_kpoints_distance.value = 0.1


    