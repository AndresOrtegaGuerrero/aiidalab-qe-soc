import ipywidgets as ipw
from pathlib import Path
import json
import base64
import numpy as np
from monty.json import jsanitize
from aiida.orm import ProjectionData
from IPython.display import clear_output, display
from aiidalab_widgets_base.utils import string_range_to_list



#BandsPdosWidget with Plotly
class BandDosPlotsWidget(ipw.VBox):
    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Select the style of plotting the projected density of states.
        </div>"""
    )

    def __init__(self, bands=None, pdos=None, **kwargs):
        if bands is None and pdos is None:
            raise ValueError("Either bands or pdos must be provided")

        self.bands = bands  # bands node
        self.pdos = pdos  # pdos node

        self.dos_atoms_group = ipw.Dropdown(
            description="Group by:",
            options=[("Kinds", "kinds"), ("Atoms", "atoms")],
            value="kinds",
        )
        self.dos_plot_group = ipw.Dropdown(
            description="Plot contributions:",
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
        )
        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(visibility="hidden"),
        )
        
        self.export_dir = Path.cwd().joinpath("exports")

        #Information for the plot
        self.dos_data = self._get_dos_data()
        self.fermi_energy = self._get_fermi_energy()
        self.bands_data = self._get_bands_data()
        self.band_labels = self._bands_labeling(self.bands_data)

        #Plotly settings
        self.bands_xaxis = self._band_xaxis()
        self.bands_yaxis = self._band_yaxis()
        self.dos_xaxis = self._dos_xaxis()
        self.dos_yaxis = self._dos_yaxis()

        #Plotly widget
        self.bandsplot_widget = self._bandsplot_widget()

        # Output widget to display the bandsplot widget
        self.bands_widget = ipw.Output()

        # Output widget to clear the specific widgets
        self.pdos_options_out = ipw.Output()


        def download_data(_=None):
            """Function to download the data."""
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

        self.pdos_options = ipw.VBox([self.description,
                self.dos_atoms_group,
                self.dos_plot_group,
                ipw.HBox([self.selected_atoms, self.wrong_syntax]),
                self.update_plot_button,])

        self._initial_view()

        # Set the event handlers
        self.download_button.on_click(download_data)
        self.update_plot_button.on_click(self._update_plot)

        super().__init__(
            children=[
                self.pdos_options_out,
                self.download_button,
                self.bands_widget,  # Add the output widget to the VBox
            ],
            **kwargs
        )
        if self.pdos:
            with self.pdos_options_out:
                display(self.pdos_options)

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
        if not self.pdos:
            return None
        expanded_selection, syntax_ok = string_range_to_list(
            self.selected_atoms.value, shift=-1
        )
        dos = get_pdos_data(
            self.pdos,
            group_tag=self.dos_atoms_group.value,
            plot_tag=self.dos_plot_group.value,
            selected_atoms=expanded_selection,
        )
        return dos

    def _get_fermi_energy(self):
        fermi_energy = self.dos_data["fermi_energy"] if self.dos_data else self.bands.band_parameters["fermi_energy"]
        return fermi_energy

    def _get_bands_data(self):
        if not self.bands:
            return None
        
        bands = export_bands_data(self.bands, self.fermi_energy)
        return bands
    
    def _bands_labeling(self, bands):
        """Function to return two lists containing the labels and values (kpoint) for plotting.
        params:
        - bands: dictionary from export_bands_data function
        output: label (list of str), label_values (list of float)
        """
        if not bands:
            return None

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


    def _band_xaxis(self):
        """Function to return the xaxis for the bands plot."""
        import plotly.graph_objects as go
        if not self.bands_data:
            return None
        paths = self.bands_data[0].get("paths")
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
    
        return bandxaxis
        
    def _band_yaxis(self):
        """Function to return the yaxis for the bands plot."""
        import plotly.graph_objects as go
        if not self.bands_data:
            return None

        bandyaxis = go.layout.YAxis(
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

        return bandyaxis

    def _dos_xaxis(self):
        """Function to return the xaxis for the dos plot."""
        import plotly.graph_objects as go

        if not self.dos_data:
            return None

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

        return dosxaxis

    def _dos_yaxis(self):
        """Function to return the yaxis for the dos plot."""
        import plotly.graph_objects as go

        if not self.dos_data:
            return None

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

        return dosyaxis
    
    def _bandsplot_widget(self):
        """Function to return the bands plot widget."""
        conditions = {
            (True, False): self._create_bands_only_plot,
            (False, True): self._create_dos_only_plot,
            (True, True): self._create_combined_plot,
        }

        return conditions.get((bool(self.bands_data), bool(self.dos_data)), None)()

    def _create_bands_only_plot(self):
        """Function to return the bands plot widget."""
        import plotly.graph_objects as go
        fig = go.Figure()
        paths = self.bands_data[0].get("paths")

        self._add_band_traces(fig, paths, "bands_only")

        for i in self.band_labels[1]:
            fig.add_vline(x=i, line=dict(color="#111111", width=1))
        fig.update_layout(xaxis=self.bands_xaxis, yaxis=self.bands_yaxis, plot_bgcolor="white", height=600, width=850,)
        return go.FigureWidget(fig)

    def _create_dos_only_plot(self):
        """Function to return the pdos plot widget."""
        
        import plotly.graph_objects as go
        fig = go.Figure()
        # Extract DOS data
        self._add_dos_traces(fig, plot_type="dos_only")
        # Add a vertical line at zero energy
        fig.add_vline(x=0, line=dict(color="#111111", width=1, dash="dot"))

        # Update the layout of the Figure
        fig.update_layout(
            xaxis=self.dos_xaxis,
            yaxis=self.dos_yaxis,
            plot_bgcolor="white",
            height=600,
            width=850,
        )

        return go.FigureWidget(fig)

    def _create_combined_plot(self):

        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        fig = make_subplots(
            rows=1, cols=2, shared_yaxes=True, column_widths=[0.7, 0.3], horizontal_spacing=0.02
        )
        paths = self.bands_data[0].get("paths")
        self._add_band_traces(fig, paths, plot_type="combined")
        self._add_dos_traces(fig, plot_type="combined")
        for i in self.band_labels[1]:
                fig.add_vline(
                    x=i,
                    line=dict(color="#111111", width=1),
                    row=1,
                    col=1,
                )
        self._customize_combined_layout(fig)
        return go.FigureWidget(fig)

    def _add_band_traces(self, fig, paths, plot_type):
        import plotly.graph_objects as go
        paths = self.bands_data[0].get("paths")
        # Convert paths to a list of Scatter objects
        scatter_objects = []
        for band in paths:
            for bands in band["values"]:
                bands_np = np.array(bands)
                scatter_objects.append(
                    go.Scatter(
                        x=band["x"],
                        y=bands_np - self.fermi_energy,
                        mode="lines",
                        line=dict(color="#111111", shape="spline", smoothing=1.3),
                        showlegend=False,
                    )
                )
        if plot_type == "bands_only":
            fig.add_traces(scatter_objects)
        else:
            rows = [1] * len(scatter_objects)
            cols = [1] * len(scatter_objects)
            fig.add_traces(scatter_objects, rows=rows,cols=cols)

    def _add_dos_traces(self, fig, plot_type):

        import plotly.graph_objects as go
        # Extract DOS data
        dos_data = self.dos_data["dos"]

        # Pre-allocate memory for Scatter objects
        num_traces = len(dos_data)
        scatter_objects = [None] * num_traces

        # Vectorize Scatter object creation
        for i, trace in enumerate(dos_data):
            dos_np = np.array(trace["x"])
            fill_prop = "tozerox" if plot_type == "combined" else "tozeroy"
            my_fill = fill_prop if trace["label"] == "Total DOS" else None
            x_data = trace["y"] if plot_type == "combined" else dos_np - self.fermi_energy
            y_data = dos_np - self.fermi_energy if plot_type == "combined" else trace["y"]
            scatter_objects[i] = go.Scatter(
                x=x_data,
                y=y_data,
                fill=my_fill,
                name=trace["label"],
                line=dict(color=trace["borderColor"], shape="spline", smoothing=1.0),
            )
        if plot_type == "dos_only":
            fig.add_traces(scatter_objects)
        else:
            rows = [1] * len(scatter_objects)
            cols = [2] * len(scatter_objects)
            fig.add_traces(scatter_objects, rows=rows,cols=cols)

    def _customize_combined_layout(self, fig):
        self._customize_layout(fig, self.bands_xaxis, self.bands_yaxis)
        self._customize_layout(fig, self.dos_xaxis, self.dos_yaxis, col=2)
        fig.update_layout(
            legend=dict(xanchor="left", x=1.04),
            height=600,
            width=850,
            plot_bgcolor="white",
        )

    def _customize_layout(self, fig, xaxis, yaxis, row=1, col=1):
        fig.update_xaxes(patch=xaxis, row=row, col=col)
        fig.update_yaxes(patch=yaxis, row=row, col=col, showticklabels=True)
        fig.add_hline(y=0, line=dict(color="#111111", width=1, dash="dot"), row=row, col=col)

    def _initial_view(self):
        with self.bands_widget:
            self._clear_output_and_display(self.bandsplot_widget)
            self.download_button.layout.visibility = "visible"

    def _update_plot(self, _=None):
        with self.bands_widget:
            expanded_selection, syntax_ok = string_range_to_list(
                self.selected_atoms.value, shift=-1
            )
            self._handle_syntax_errors(syntax_ok)
            self.dos_data = self._get_dos_data()
            self.bandsplot_widget = self._bandsplot_widget()
            self._clear_output_and_display(self.bandsplot_widget)

    def _clear_output_and_display(self, widget=None):
        clear_output(wait=True)
        if widget:
            display(widget)

    def _handle_syntax_errors(self, syntax_ok):
        if not syntax_ok:
            self.wrong_syntax.layout.visibility = "visible"
            clear_output(wait=True)

        


def get_pdos_data(pdos, group_tag, plot_tag, selected_atoms):
    dos = []

    if "output_dos" not in pdos.dos:
        return None

    _, energy_dos, _ = pdos.dos.output_dos.get_x()
    tdos_values = {f"{n}": v for n, v, _ in pdos.dos.output_dos.get_y()}

    if "projections" in pdos.projwfc:
        # Total DOS
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
            pdos.projwfc.projections,
            spin_type="none",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )
    else:
        # Total DOS (↑) and Total DOS (↓)
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
            "y": (-tdos_values.get("dos_spin_down")).tolist(),
            "borderColor": "#8A8A8A",  # dark gray
            "backgroundColor": "#999999",  # light gray
            "backgroundAlpha": "40%",
            "lineStyle": "dash",
        }
        dos += [tdos_up, tdos_down]

        # Spin-up (↑) and Spin-down (↓)
        dos += _projections_curated_options(
            pdos.projwfc.projections_up,
            spin_type="up",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )
        dos += _projections_curated_options(
            pdos.projwfc.projections_down,
            spin_type="down",
            line_style="dash",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )

    data_dict = {
        "fermi_energy": pdos.nscf.output_parameters["fermi_energy"],
        "dos": dos,
    }

    return json.loads(json.dumps(data_dict))


def _projections_curated_options(
    projections: ProjectionData,
    group_tag,
    plot_tag,
    selected_atoms,
    spin_type="none",
    line_style="solid",
):
    """Collects data from ProjectionData and parses it as a dos list for the bandsplot widget."""
    
    _pdos = {}
    list_positions = []

    # Setting
    dict_html = {
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
            orbital_name_plotly = dict_html.get(orbital_name, orbital_name)
        except AttributeError:
            # Handling the exception
            orbital_name = "j {j} l {l} m_j{m_j}".format(
                j=orbital_data["total_angular_momentum"],
                l=orbital_data["angular_momentum"],
                m_j=orbital_data["magnetic_number"],
            )
            orbital_name_plotly = "j={j} <i>l</i>={l} m<sub>j</sub>={m_j}".format(
                j=dict_html.get(orbital_data["total_angular_momentum"], ''),
                l=orbital_data["angular_momentum"],
                m_j=dict_html.get(orbital_data["magnetic_number"], ''),
            )

        if not selected_atoms:
            key = _get_key(group_tag, plot_tag, atom_position, kind_name, orbital_name_plotly)

            if key:
                _update_pdos_dict(_pdos, key, pdos, energy)
        else:
            try:
                index = list_positions.index(atom_position)

                if index in selected_atoms:
                    key = _get_key(group_tag, plot_tag, atom_position, kind_name, orbital_name_plotly)

                    if key:
                        _update_pdos_dict(_pdos, key, pdos, energy)
            except ValueError:
                pass

    dos = _format_dos(_pdos, spin_type, cmap, line_style)

    return dos


def _get_key(group_tag, plot_tag, atom_position, kind_name, orbital_name_plotly):
    """Generates the key based on group_tag and plot_tag."""
    
    key_formats = {
        ("atoms", "total"): r"{var}",
        ("kinds", "total"): r"{var1}",
        ("atoms", "orbital"): r"{var1}-{var}<br>{var3}",
        ("kinds", "orbital"): r"{var1}<br>{var3}",
    }

    key = key_formats.get((group_tag, plot_tag))
    if key is not None:
        return key.format(var=atom_position, var1=kind_name, var3=orbital_name_plotly)
    else:
        return None

def _update_pdos_dict(_pdos, key, pdos, energy):
    """Updates the _pdos dictionary with the given key, pdos, and energy."""
    if key in _pdos:
        _pdos[key][1] += pdos
    else:
        _pdos[key] = [energy, pdos]


def _format_dos(_pdos, spin_type, cmap, line_style):
    """Formats the dos list based on _pdos, spin_type, cmap, and line_style."""
    dos = []

    for label, (energy, pdos) in _pdos.items():
        if spin_type == "down":
            # Invert y-axis
            pdos = -pdos
            label = f"{label} (↓)"
        elif spin_type == "up":
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


def export_bands_data(outputs, fermi_energy=None):
    if not "band_structure" in outputs:
        return None
    
    data = json.loads(
        outputs.band_structure._exportcontent(
            "json", comments=False
        )[0]
    )
    # The fermi energy from band calculation is not robust.
    data["fermi_level"] = (
        outputs.band_parameters["fermi_energy"] or fermi_energy 
    )
    return [
        jsanitize(data),
    ]

    
def cmap(label: str) -> str:
    """Return RGB string of color for given pseudo info
    Hardcoded at the momment.
    """
    import random

    # if a unknow type generate random color based on ascii sum
    ascn = sum([ord(c) for c in label])
    random.seed(ascn)

    return "#%06x" % random.randint(0, 0xFFFFFF)