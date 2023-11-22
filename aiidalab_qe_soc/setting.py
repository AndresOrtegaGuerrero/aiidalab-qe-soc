"""Panel for SOC plugin.
"""
import ipywidgets as ipw
from aiidalab_qe.common.panel import Panel

class Setting(Panel):
    

    title = "SOC Settings"
    identifier = "soc"
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

        self.peudos_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Spin Orbit Coupling calculations are conducted using the PseudoDojo library.
            </div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Select if you want a Band structure , Projected density of states or both to be computed.
            </div>"""
        )
        # the widget for DFT functional selection
        self.dft_functional = ipw.Dropdown(
            description="Functional:",
            options=["PBE", "PBEsol"],
            style={"description_width": "initial"},
        )

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

        self.children=[
                self.settings_title,
                self.settings_help,
                self.calc_options,
                self.dft_functional,
                self.kpath_2d,
                self.peudos_help,
                self.soc_ecutwfc,
                self.soc_ecutrho,
                self.nscf_kpoints_distance,
            ]
        super().__init__(**kwargs)

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "calc_options": self.calc_options.value,
            "dft_functional": self.dft_functional.value,
            "soc_ecutwfc": self.soc_ecutwfc.value,
            "soc_ecutrho": self.soc_ecutrho.value,
            "kpath_2d": self.kpath_2d.value,
            "nscf_kpoints_distance": self.nscf_kpoints_distance.value,
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.kpath_2d.value = input_dict.get("kpath_2d", "hexagonal")
        self.calc_options.value = input_dict.get("calc_options", "bands")
        self.dft_functional.value = input_dict.get("dft_functional", "PBE")
        self.soc_ecutwfc.value = input_dict.get("soc_ecutwfc", 90.0)
        self.soc_ecutrho.value = input_dict.get("soc_ecutrho", 360.0)
        self.nscf_kpoints_distance.value = input_dict.get("nscf_kpoints_distance", 0.1)

    def reset(self):
        """Reset the panel to its default values."""
        self.kpath_2d.value = "hexagonal"
        self.calc_options.value = "bands"
        self.dft_functional.value = "PBE"
        self.soc_ecutwfc.value = 90.0
        self.soc_ecutrho.value = 360.0
        self.nscf_kpoints_distance.value = 0.1


    