import ipywidgets as ipw
from pathlib import Path
import json
import base64
import numpy as np
from monty.json import jsanitize
from aiida.orm import ProjectionData
from IPython.display import clear_output, display
from aiidalab_widgets_base.utils import string_range_to_list
import random

class BandDosPlotsWidget(ipw.VBox):
    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Select the style of plotting the projected density of states.
        </div>"""
    )

    def __init__(self, node, **kwargs):
        self.node = node

        self.dos_atoms_group = ipw.Dropdown(
            options=[("Kinds", "kinds"), ("Atoms", "atoms")],
            value="kinds",
        )
        self.dos_plot_group = ipw.Dropdown(
            options=[("Total", "total"), ("Orbital", "orbital")],
            value="total",
        )
        self.selected_atoms = ipw.Text(
            description="Select atoms:",
            value="",
            style={"description_width": "initial"},
        )
        self.wrong_syntax = ipw.HTML(
            value="""<i class="fa fa-times" style="color:red;font-size:2em;" ></i> wrong syntax""",
            layout={"visibility": "hidden"},
        )
        self.update_plot_button = ipw.Button(
            description="Update Plot",
            icon="pencil",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )
        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto", visibility="hidden"),
        )
        self.export_dir = Path.cwd().joinpath("exports")
        self.dos_data = self._get_dos_data()
        self.fermi_energy = self._get_fermi_energy()
        self.bands_data = self._get_bands_data()
        self.band_labels = self.bands_labeling(self.bands_data)
        self.bands_xaxis = self._band_xaxis()
        self.bands_yaxis = self._band_yaxis()
        self.dos_xaxis = self._dos_xaxis()
        self.dos_yaxis = self._dos_yaxis()
        self.bandsplot_widget = self._bandsplot_widget()

        if self.bands_data and not self.dos_data:
            self.update_plot_button.disabled = True
        self.bands_widget = ipw.Output()

        def download_data(_=None):
            file_name_bands = "bands_data.json"
            file_name_dos = "dos_data.json"
            if self.bands_data:
                json_str = json.dumps(self.bands_data)
                b64_str = base64.b64encode(json_str.encode()).decode()
                self._download(payload=b64_str, filename=file_name_bands)
            if self.dos_data:
                json_str = json.dumps(self.dos_data)
                b64_str = base64.b64encode(json_str.encode()).decode()
                self._download(payload=b64_str, filename=file_name_dos)

        self.download_button.on_click(download_data)
        self._initial_view()
        self.update_plot_button.on_click(self._update_plot)

        super().__init__(
            children=[
                self.description,
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Group :",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="80px"
                            ),
                        ),
                        self.dos_atoms_group,
                    ]
                ),
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Plot :",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="80px"
                            ),
                        ),
                        self.dos_plot_group,
                    ]
                ),
                ipw.HBox([self.selected_atoms, self.wrong_syntax]),
                ipw.HBox(
                    children=[
                        self.update_plot_button,
                        self.download_button,
                    ]
                ),
                self.bands_widget,
            ]
        )

    @staticmethod
    def _download(payload, filename):
        """Download payload as a file named as filename."""
        from IPython.display import Javascript

        javas = Javascript(
            """
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """.format(
                payload=payload, filename=filename
            )
        )
        display(javas)

    def _get_dos_data(self):
        if "pdos" in self.node:
            expanded_selection, syntax_ok = string_range_to_list(
                self.selected_atoms.value, shift=-1
            )
            dos = get_pdos_data(
                self.node.pdos,
                group_tag=self.dos_atoms_group.value,
                plot_tag=self.dos_plot_group.value,
                selected_atoms=expanded_selection,
            )
            return dos
        else:
            return None

    def _get_fermi_energy(self):
        fermi_energy = self.dos_data["fermi_energy"] if self.dos_data else 0.0
        return fermi_energy

    def _get_bands_data(self):
        if "bands" in self.node:
            bands = export_bands_data(self.node, self.fermi_energy)
            return bands
        else:
            return None

    def _initial_view(self):
        with self.bands_widget:
            clear_output(wait=True)
            # self.bandsplot_widget.show() #Fix plotly not showing
            display(self.bandsplot_widget)
            self.download_button.layout.visibility = "visible"

    def _update_plot(self, _=None):
        with self.bands_widget:
            expanded_selection, syntax_ok = string_range_to_list(
                self.selected_atoms.value, shift=-1
            )
            if not syntax_ok:
                self.wrong_syntax.layout.visibility = "visible"
                clear_output(wait=True)
            else:
                self.dos_data = self._get_dos_data()
                self.bandsplot_widget = self._bandsplot_widget()
                clear_output(wait=True)
                # self.bandsplot_widget.show() #Fix plotly not showing
                display(self.bandsplot_widget)

    def _bandsplot_widget(self):
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        if self.bands_data and not self.dos_data:
            fig = go.Figure()
            paths = self.bands_data[0].get("paths")

            for band in paths:
                for bands in band["values"]:
                    bands_np = np.array(bands)
                    fig.add_trace(
                        go.Scatter(
                            x=band["x"],
                            y=bands_np - self.fermi_energy,  # substract Fermi Energy
                            mode="lines",
                            line=dict(color="#111111", shape="spline", smoothing=1.3),
                            showlegend=False,
                        )
                    )
            for i in self.band_labels[1]:
                fig.add_vline(x=i, line=dict(color="#111111", width=1))
            if self.fermi_energy != 0:
                fig.add_hline(y=0, line=dict(color="#111111", width=1, dash="dot"))
            fig.update_layout(
                height=600,
                width=950,
                plot_bgcolor="white",
                xaxis=self.bands_xaxis,
                yaxis=self.bands_yaxis,
            )

        elif self.dos_data and not self.bands_data:
            fig = go.Figure()
            for trace in self.dos_data["dos"]:
                if trace["label"] == "Total DOS":
                    my_fill = "tozeroy"
                else:
                    my_fill = None
                dos_np = np.array(trace["x"])
                fig.add_trace(
                    go.Scatter(
                        x=dos_np - self.fermi_energy,
                        y=trace["y"],
                        fill=my_fill,
                        name=trace["label"],
                        line=dict(
                            color=trace["borderColor"], shape="spline", smoothing=1.0
                        ),
                    )
                )
            fig.add_vline(x=0, line=dict(color="#111111", width=1, dash="dot"))
            fig.update_layout(
                height=600,
                width=850,
                plot_bgcolor="white",
                xaxis=self.dos_xaxis,
                yaxis=self.dos_yaxis,
            )

        elif self.dos_data and self.bands_data:
            fig = make_subplots(
                rows=1,
                cols=2,
                shared_yaxes=True,
                column_widths=[0.7, 0.3],
                horizontal_spacing=0.02,
            )
            paths = self.bands_data[0].get("paths")
            for band in paths:
                for bands in band["values"]:
                    bands_np = np.array(bands)
                    fig.add_trace(
                        go.Scatter(
                            x=band["x"],
                            y=bands_np - self.fermi_energy,  # substract Fermi Energy
                            mode="lines",
                            line=dict(color="#111111", shape="spline", smoothing=1.3),
                            showlegend=False,
                        ),
                        row=1,
                        col=1,
                    )

            for trace in self.dos_data["dos"]:
                if trace["label"] == "Total DOS":
                    my_fill = "tozerox"
                else:
                    my_fill = None

                dos_np = np.array(trace["x"])
                fig.add_trace(
                    go.Scatter(
                        x=trace["y"],
                        y=dos_np - self.fermi_energy,
                        fill=my_fill,
                        name=trace["label"],
                        line=dict(
                            color=trace["borderColor"], shape="spline", smoothing=1.3
                        ),
                    ),
                    row=1,
                    col=2,
                )
            for i in self.band_labels[1]:
                fig.add_vline(
                    x=i,
                    line=dict(color="#111111", width=1),
                    row=1,
                    col=1,
                )
            fig.update_xaxes(
                patch=self.bands_xaxis,
                row=1,
                col=1,
            )
            fig.update_yaxes(
                patch=self.bands_yaxis,
                row=1,
                col=1,
            )
            fig.update_xaxes(
                patch=self.dos_xaxis,
                row=1,
                col=2,
            )
            fig.update_yaxes(patch=self.dos_yaxis, row=1, col=2, showticklabels=False)

            fig.add_hline(
                y=0, line=dict(color="#111111", width=1, dash="dot"), row=1, col=1
            )
            fig.add_hline(
                y=0, line=dict(color="#111111", width=1, dash="dot"), row=1, col=2
            )
            fig.update_layout(
                height=600,
                width=850,
                plot_bgcolor="white",
                legend=dict(
                    # yanchor="top",
                    # y=0.99,
                    xanchor="left",
                    x=1.04,
                ),
            )

        else:
            fig = None

        return go.FigureWidget(fig)

    def _band_xaxis(self):
        import plotly.graph_objects as go

        if self.bands_data:
            paths = self.bands_data[0].get("paths")
            # labels = bands_labeling(self.bands_data)
            # labels, labels_values = bands_labeling(self.bands_data)
            slider_bands = go.layout.xaxis.Rangeslider(
                thickness=0.08,
                range=[0, paths[-1]["x"][-1]],
            )

            bandxaxis = go.layout.XAxis(
                title="k-points",
                range=[0, paths[-1]["x"][-1]],
                showgrid=True,
                showline=True,
                tickmode="array",
                rangeslider=slider_bands,
                fixedrange=False,
                tickvals=self.band_labels[1],
                ticktext=self.band_labels[0],
                showticklabels=True,
                linecolor="#111111",
                mirror=True,
                linewidth=2,
                type="linear",
            )

        else:
            bandxaxis = None
        return bandxaxis

    def _band_yaxis(self):
        import plotly.graph_objects as go

        if self.bands_data:
            bandyaxis = go.layout.YAxis(
                # title="Electronic Bands(eV)",
                title=dict(text="Electronic Bands (eV)", standoff=1),
                side="left",
                showgrid=True,
                showline=True,
                zeroline=True,
                fixedrange=False,
                automargin=True,
                mirror="ticks",
                ticks="inside",
                linewidth=2,
                linecolor="#111111",
                tickwidth=2,
                zerolinewidth=2,
            )

        else:
            bandyaxis = None
        return bandyaxis

    def _dos_xaxis(self):
        import plotly.graph_objects as go

        if self.dos_data:
            if self.bands_data:
                dosxaxis = go.layout.XAxis(
                    title="Density of states",
                    side="bottom",
                    showgrid=True,
                    showline=True,
                    linecolor="#111111",
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                )

            else:
                dosxaxis = go.layout.XAxis(
                    title="Density of states (eV)",
                    showgrid=True,
                    showline=True,
                    linecolor="#111111",
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                )

        else:
            dosxaxis = None
        return dosxaxis

    def _dos_yaxis(self):
        import plotly.graph_objects as go

        if self.dos_data:
            if self.bands_data:
                dosyaxis = go.layout.YAxis(
                    # title="Density of states (eV)", #FigureWidget modifies the title position in the meantime lets put it
                    showgrid=True,
                    showline=True,
                    side="right",
                    # position=0.0,
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                    linecolor="#111111",
                    zerolinewidth=2,
                )

            else:
                dosyaxis = go.layout.YAxis(
                    # title="Density of states (eV)",
                    showgrid=True,
                    showline=True,
                    side="left",
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                    linecolor="#111111",
                    zerolinewidth=2,
                )

        else:
            dosyaxis = None
        return dosyaxis

    def bands_labeling(self, bands):
        """Function to return two list containing the labels and values (kpoint) for plotting
        params: bands: dictionary from export_bands_data function
        output: label (list of str), label_values (list of float)
        """
        if bands:
            my_paths = bands[0].get("paths")
            my_labels = []
            for path in my_paths:  # Remove duplicates
                label_a = [path["from"], path["x"][0]]
                label_b = [path["to"], path["x"][-1]]
                if label_a not in my_labels:
                    my_labels.append(label_a)
                if label_b not in my_labels:
                    my_labels.append(label_b)

            my_clean_labels = []  # Format
            for i in my_labels:
                if my_clean_labels:
                    if i not in my_clean_labels:
                        if my_clean_labels[-1][-1] == i[1]:
                            my_clean_labels[-1][0] = my_clean_labels[-1][0] + "|" + i[0]
                        else:
                            my_clean_labels.append(i)
                else:
                    my_clean_labels.append(i)

            labels = [label[0] for label in my_clean_labels]
            labels_values = [label[1] for label in my_clean_labels]
            return [labels, labels_values]
        else:
            return None


def get_pdos_data(work_chain_node, group_tag, plot_tag, selected_atoms):
    if "dos" in work_chain_node:
        _, energy_dos, _ = work_chain_node.dos.output_dos.get_x()
        tdos_values = {f"{n}": v for n, v, _ in work_chain_node.dos.output_dos.get_y()}

        dos = []

        if "projections" in work_chain_node.projwfc:
            # The total dos parsed
            tdos = {
                "label": "Total DOS",
                "x": energy_dos.tolist(),
                "y": tdos_values.get("dos").tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            dos.append(tdos)

            dos += _projections_curated_options(
                work_chain_node.projwfc.projections,
                spin_type="none",
                group_tag=group_tag,
                plot_tag=plot_tag,
                selected_atoms=selected_atoms,
            )

        else:
            # The total dos parsed
            tdos_up = {
                "label": "Total DOS (↑)",
                "x": energy_dos.tolist(),
                "y": tdos_values.get("dos_spin_up").tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            tdos_down = {
                "label": "Total DOS (↓)",
                "x": energy_dos.tolist(),
                "y": (-tdos_values.get("dos_spin_down")).tolist(),  # minus
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "dash",
            }
            dos += [tdos_up, tdos_down]

            # spin-up (↑)
            dos += _projections_curated_options(
                work_chain_node.projwfc.projections_up,
                spin_type="up",
                group_tag=group_tag,
                plot_tag=plot_tag,
                selected_atoms=selected_atoms,
            )

            # spin-dn (↓)
            dos += _projections_curated_options(
                work_chain_node.projwfc.projections_down,
                spin_type="down",
                line_style="dash",
                group_tag=group_tag,
                plot_tag=plot_tag,
                selected_atoms=selected_atoms,
            )

        data_dict = {
            "fermi_energy": work_chain_node.nscf.output_parameters["fermi_energy"],
            "dos": dos,
        }

        return json.loads(json.dumps(data_dict))

    else:
        return None


def _projections_curated_options(
    projections: ProjectionData,
    group_tag,
    plot_tag,
    selected_atoms,
    spin_type="none",
    line_style="solid",
):
    """Collect the data from ProjectionData and parse it as dos list which can be
    understand by bandsplot widget. `curated_tag` is for which tag to be grouped, by atom or by orbital name.
    The spin_type is used to invert all the y values of pdos to be shown as spin down pdos and to set label.
    """
    _pdos = {}
    list_positions = []
    # Setting
    dict_html = {
        "s": "s",
        "pz": "p<sub>z</sub>",
        "px": "p<sub>x</sub>",
        "py": "p<sub>y</sub>",
        "dz2": "d<sub>z<sup>2</sup></sub>",
        "dxy": "d<sub>xy</sub>",
        "dxz": "d<sub>xz</sub>",
        "dyz": "d<sub>yz</sub>",
        "dx2-y2": "d<sub>x<sup>2</sup>-y<sup>2</sup></sub>",
        "fz3": "f<sub>z<sup>3</sup></sub>",
        "fxz2": "f<sub>xz<sup>2</sup></sub>",
        "fyz2": "f<sub>yz<sup>2</sup></sub>",
        "fxyz": "f<sub>xzy</sub>",
        "fx(x2-3y2)": "f<sub>x(x<sup>2</sup>-3y<sup>2</sup>)</sub>",
        "fy(3x2-y2)": "f<sub>y(3x<sup>2</sup>-y<sup>2</sup>)</sub>",
        "fy(x2-z2)": "f<sub>y(x<sup>2</sup>-z<sup>2</sup>)</sub>",
        0.5: "<sup>+1</sup>/<sub>2</sub>",
        -0.5: "<sup>-1</sup>/<sub>2</sub>",
        1.5: "<sup>+3</sup>/<sub>2</sub>",
        -1.5: "<sup>-3</sup>/<sub>2</sub>",
        2.5: "<sup>+5</sup>/<sub>2</sub>",
        -2.5: "<sup>-5</sup>/<sub>2</sub>",
    }
    for orbital, pdos, energy in projections.get_pdos():
        orbital_data = orbital.get_orbital_dict()
        kind_name = orbital_data["kind_name"]
        atom_position = [round(i, 2) for i in orbital_data["position"]]
        if atom_position not in list_positions:
            list_positions.append(atom_position)
        try:
            orbital_name = orbital.get_name_from_quantum_numbers(
                orbital_data["angular_momentum"], orbital_data["magnetic_number"]
            ).lower()
            if orbital_name in dict_html:
                orbital_name_plotly = dict_html[orbital_name]
            else:
                orbital_name_plotly = orbital_name

        except AttributeError:
            orbital_name = "j {j} l {l} m_j{m_j}".format(
                j=orbital_data["total_angular_momentum"],
                l=orbital_data["angular_momentum"],
                m_j=orbital_data["magnetic_number"],
            )
            orbital_name_plotly = "j={j} <i>l</i>={l} m<sub>j</sub>={m_j}".format(
                j=dict_html[orbital_data["total_angular_momentum"]],
                l=orbital_data["angular_momentum"],
                m_j=dict_html[orbital_data["magnetic_number"]],
            )
        if not selected_atoms:
            if group_tag == "atoms" and plot_tag == "total":
                key = r"{var}".format(var=atom_position)
            elif group_tag == "kinds" and plot_tag == "total":
                key = r"{var1}".format(var1=kind_name)
            elif group_tag == "atoms" and plot_tag == "orbital":
                key = r"{var1}<br>{var2}-{var3}".format(
                    var1=atom_position, var2=kind_name, var3=orbital_name_plotly
                )
            elif group_tag == "kinds" and plot_tag == "orbital":
                key = r"{var1}-{var2}".format(var1=kind_name, var2=orbital_name_plotly)
            else:
                key = None

            if key:
                if key in _pdos:
                    _pdos[key][1] += pdos
                else:
                    _pdos[key] = [energy, pdos]

        else:
            try:
                index = list_positions.index(atom_position)
                if index in selected_atoms:
                    if group_tag == "atoms" and plot_tag == "total":
                        key = r"{var}".format(var=atom_position)
                    elif group_tag == "kinds" and plot_tag == "total":
                        key = r"{var1}".format(var1=kind_name)
                    elif group_tag == "atoms" and plot_tag == "orbital":
                        key = r"{var1}<br>{var2}-{var3}".format(
                            var1=atom_position, var2=kind_name, var3=orbital_name_plotly
                        )
                    elif group_tag == "kinds" and plot_tag == "orbital":
                        key = r"{var1}-{var2}".format(
                            var1=kind_name, var2=orbital_name_plotly
                        )
                    else:
                        key = None

                    if key:
                        if key in _pdos:
                            _pdos[key][1] += pdos
                        else:
                            _pdos[key] = [energy, pdos]

            except ValueError:
                pass

    dos = []
    for label, (energy, pdos) in _pdos.items():
        if spin_type == "down":
            # invert y-axis
            pdos = -pdos
            label = f"{label} (↓)"

        if spin_type == "up":
            label = f"{label} (↑)"

        orbital_pdos = {
            "label": label,
            "x": energy.tolist(),
            "y": pdos.tolist(),
            "borderColor": cmap(label),
            "lineStyle": line_style,
        }
        dos.append(orbital_pdos)

    return dos


def export_bands_data(work_chain_node, fermi_energy=None):
    if "band_structure" in work_chain_node.bands:
        data = json.loads(
            work_chain_node.bands.band_structure._exportcontent(
                "json", comments=False
            )[0]
        )
        # The fermi energy from band calculation is not robust.
        data["fermi_level"] = (
            fermi_energy or work_chain_node.bands.band_parameters["fermi_energy"]
        )
        return [
            jsanitize(data),
        ]
    else:
        return None

def cmap(label: str) -> str:
    """Return RGB string of color for given pseudo info
    Hardcoded at the momment.
    """
    # if a unknow type generate random color based on ascii sum
    ascn = sum([ord(c) for c in label])
    random.seed(ascn)

    return "#%06x" % random.randint(0, 0xFFFFFF)