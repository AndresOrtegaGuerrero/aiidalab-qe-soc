from aiida.plugins import WorkflowFactory
from aiida.engine import ToContext, WorkChain, calcfunction, if_
from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida import orm
from aiida.common import AttributeDict

PdosWorkChain = WorkflowFactory("quantumespresso.pdos")
PwBandsWorkChain = WorkflowFactory("quantumespresso.pw.bands")

class SOCWorkChain(WorkChain):
    "WorkChain to compute vibrational property of a crystal."
    label = "soc"

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('structure', valid_type=orm.StructureData,
                   help='The inputs structure.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False),
                   help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.input('properties', valid_type=orm.List, default=lambda: orm.List(),
                   help='The properties to calculate, used to control the logic of SOCWorkChain.')
        spec.expose_inputs(PwBandsWorkChain, namespace='bands',
                           exclude=('clean_workdir', 'structure', 'relax'),
                           namespace_options={'required': False, 'populate_defaults': False,
                                              'help': 'Inputs for the `PwBandsWorkChain`.'})
        spec.expose_inputs(PdosWorkChain, namespace='pdos',
                           exclude=('clean_workdir', 'structure'),
                           namespace_options={'required': False, 'populate_defaults': False,
                                              'help': 'Inputs for the `PdosWorkChain`.'})
        spec.outline(
            cls.setup,
            if_(cls.should_run_bands)(
                cls.run_bands,
                cls.inspect_bands,
            ),
            if_(cls.should_run_pdos)(
                cls.run_pdos,
                cls.inspect_pdos,
            ),
            cls.results,
        )
        spec.expose_outputs(
            PwBandsWorkChain, namespace='bands',
            namespace_options={'required': False, 'help':'Outputs of the `PwBandsWorkChain`.'},
        )
        spec.expose_outputs(
            PdosWorkChain, namespace='pdos',
            namespace_options={'required': False, 'help':'Outputs of the `PdosWorkChain`.'},
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_BANDS',
                       message='the PwBandsWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_PDOS',
                          message='the PdosWorkChain sub process failed')


    @classmethod
    def get_builder_from_protocol(
        cls,
        structure,
        pw_code,
        dos_code,
        projwfc_code,
        protocol,
        properties,
        clean_workdir,
        functional="PBE",
        overrides=None,
        **kwargs):
        """Return a builder prepopulated with inputs selected according to the protocol."""

        overrides = overrides or {}
        builder = cls.get_builder()

        #Use only Fully Relativistic pseudos
        if functional == "PBE":
            family_fr_pseudo = orm.load_group("PseudoDojo/0.4/PBE/FR/stringent/upf")
        else:
            family_fr_pseudo = orm.load_group("PseudoDojo/0.4/PBEsol/FR/stringent/upf")

        #Set spin_orbit coupling
        for key in ["bands", "pdos"]:
            for calc in ["scf", "nscf", "bands"]:
                if calc in overrides[key]:
                    overrides[key][calc]["pw"]["pseudos"] = family_fr_pseudo.get_pseudos(structure=structure)
                    overrides[key][calc]["pw"]["parameters"]["SYSTEM"]["lspinorb"] = True
                    overrides[key][calc]["pw"]["parameters"]["SYSTEM"]["noncolin"] = True
                    overrides[key][calc]["pw"]["metadata"] =  {
                            "options": {"max_wallclock_seconds": 82800}
                        }

        # Set the structure
        builder.structure = structure

        bands_overrides = overrides.pop('bands', {})

        # Bands workchain settings
        soc_bands = PwBandsWorkChain.get_builder_from_protocol(
            structure=structure,
            code=pw_code,
            protocol=protocol,
            overrides=bands_overrides,
            **kwargs
        )

        # pop the inputs that are excluded from the exposed inputs of the bands workchain
        soc_bands.pop('clean_workdir', None)
        soc_bands.pop('structure', None)
        soc_bands.pop('relax', None)
        soc_bands.scf["pw"]["parameters"]["SYSTEM"].pop("nspin", None)
        soc_bands.bands["pw"]["parameters"]["SYSTEM"].pop("nspin", None)

        if structure.pbc != (True, True, True):
            soc_bands.pop("bands_kpoints_distance")
            bands.update({"bands_kpoints": bands_overrides["bands"]["kpoints"]})

        builder.bands = soc_bands

        # Pdos workchain settings

        if dos_code is not None and projwfc_code is not None:
            pdos_overrides = overrides.pop('pdos', {})
            soc_pdos = PdosWorkChain.get_builder_from_protocol(
                structure=structure,
                pw_code = pw_code,
                dos_code=dos_code,
                projwfc_code=projwfc_code,
                protocol=protocol,
                overrides=pdos_overrides,
                **kwargs
            )
            soc_pdos.pop('clean_workdir', None)
            soc_pdos.pop('structure', None)
            soc_pdos.scf["pw"]["parameters"]["SYSTEM"].pop("nspin", None)
            soc_pdos.nscf["pw"]["parameters"]["SYSTEM"].pop("nspin", None)
            builder.pdos = soc_pdos
        
        # Set the properties
        builder.properties = orm.List(list=properties)

        # Set the clean_workdir
        builder.clean_workdir = orm.Bool(clean_workdir)

        return builder

    def setup(self):
        """Define the current structure and the properties to calculate."""
        self.ctx.current_structure = self.inputs.structure
        self.ctx.properties = self.inputs.properties

        #logic to decide if bands should be run
        self.ctx.run_bands = "bands" in self.ctx.properties
        self.ctx.run_pdos = "pdos" in self.ctx.properties
    
    def should_run_bands(self):
        """Return whether a bands calculation should be run."""
        return self.ctx.run_bands
    
    def run_bands(self):
        """Run the bands calculation."""
        inputs = AttributeDict(self.exposed_inputs(PwBandsWorkChain, namespace="bands")) 
        inputs.metadata.call_link_label = "bands"
        inputs.structure = self.ctx.current_structure
        running = self.submit(PwBandsWorkChain, **inputs)
        self.report(f'launching PwBandsWorkChain<{running.pk}>')

        return ToContext(workchain_bands=running)

    def inspect_bands(self):
        """Verify that the bands calculation finished successfully."""
        workchain = self.ctx.workchain_bands

        if not workchain.is_finished_ok:
            self.report(f'bands workchain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_BANDS

        scf = (
            workchain.get_outgoing(orm.WorkChainNode, link_label_filter="scf")
            .one()
            .node
        )
        try:
            self.ctx.current_structure = workchain.outputs.primitive_structure
        except AttributeError:
            self.ctx.current_structure = workchain.inputs.structure
        self.ctx.scf_parent_folder = scf.outputs.remote_folder

    def should_run_pdos(self):
        """Return whether a pdos calculation should be run."""
        return self.ctx.run_pdos
    
    def run_pdos(self):
        """Run the pdos calculation."""
        inputs = AttributeDict(self.exposed_inputs(PdosWorkChain, namespace="pdos")) 
        inputs.metadata.call_link_label = "pdos"
        inputs.structure = self.ctx.current_structure
        inputs.nscf.pw.parameters = inputs.nscf.pw.parameters.get_dict()

        if self.ctx.scf_parent_folder:
            inputs.pop("scf")
            inputs.nscf.pw.parent_folder = self.ctx.scf_parent_folder
        
        inputs = prepare_process_inputs(PdosWorkChain, inputs)
        running = self.submit(PdosWorkChain, **inputs)
        self.report(f'launching PdosWorkChain<{running.pk}>')

        return ToContext(workchain_pdos=running)

    def inspect_pdos(self):
        """Verify that the pdos calculation finished successfully."""
        workchain = self.ctx.workchain_pdos

        if not workchain.is_finished_ok:
            self.report(f'pdos workchain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PDOS

    def results(self):

        if self.ctx.run_bands:
            self.out_many(self.exposed_outputs(self.ctx.workchain_bands, PwBandsWorkChain, namespace="bands"))
        if self.ctx.run_pdos:
            self.out_many(self.exposed_outputs(self.ctx.workchain_pdos, PdosWorkChain, namespace="pdos"))

