from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiidalab_qe.plugins.bands.workchain import generate_kpath_1d, generate_kpath_2d
SOCWorkChain = WorkflowFactory("soc_app.soc")

def get_builder(codes,structure, parameters):

    SOC_PROPERTIES = {
        "bands": ["bands"],
        "pdos": ["pdos"],
        "bands_pdos": ["bands", "pdos"],
    }

    from copy import deepcopy
    protocol = parameters["workchain"].pop("protocol", "moderate")
    pw_code = codes.pop("pw")
    dos_code = codes.pop("dos")
    projwfc_code = codes.pop("projwfc")

    #scf overrides
    scf_overrides = deepcopy(parameters["advanced"])

    #bands overrides
    bands_overrides = deepcopy(parameters["advanced"])
    bands_overrides.pop("kpoints_distance", None)
    bands_overrides["pw"]["parameters"]["SYSTEM"].pop("smearing", None)
    bands_overrides["pw"]["parameters"]["SYSTEM"].pop("degauss", None)

    if structure.pbc != (True, True, True):
        kpoints_distance = parameters["advanced"]["kpoints_distance"]
        if structure.pbc == (True, False, False):
            kpoints = generate_kpath_1d(structure, kpoints_distance)
        elif structure.pbc == (True, True, False):
            kpoints = generate_kpath_2d(
                structure, kpoints_distance, parameters["soc"]["kpath_2d"]
            )
        bands_overrides.pop("bands_kpoints_distance", None)
        bands_overrides.update({"bands_kpoints": kpoints})

    nscf_overrides = deepcopy(parameters["advanced"])
    nscf_overrides["kpoints_distance"] = parameters["soc"]["nscf_kpoints_distance"]

    overrides = {
        "bands": {
            "scf": scf_overrides,
            "bands": bands_overrides,
        },
        "pdos": {
            "scf": scf_overrides,
            "nscf": nscf_overrides,
        },
    }

    builder = SOCWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        dos_code=dos_code,
        projwfc_code=projwfc_code,
        structure=structure,
        protocol=protocol,
        overrides=overrides,
        properties=SOC_PROPERTIES[parameters["soc"]["calc_options"]],
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
        clean_workdir=False,
        )



    return builder


workchain_and_builder = {
    "workchain": SOCWorkChain,
    "get_builder": get_builder,
}
