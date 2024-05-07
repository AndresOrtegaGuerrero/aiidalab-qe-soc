from aiida.plugins import WorkflowFactory
from aiida import orm
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiidalab_qe.plugins.bands.workchain import generate_kpath_1d, generate_kpath_2d

SOCWorkChain = WorkflowFactory("soc_app.soc")


def set_resources(builder_attribute, code_details):
    """
    Set the resources and parallelization for a given attribute of the builder.

    Parameters:
        builder_attribute: The attribute of the builder to update (e.g., builder.scf.pw).
        code_details: A dictionary containing the nodes, ntasks_per_node, and cpus_per_task.
    """
    builder_attribute.metadata.options.resources = {
        "num_machines": code_details["nodes"],
        "num_mpiprocs_per_machine": code_details["ntasks_per_node"],
        "num_cores_per_mpiproc": code_details["cpus_per_task"],
    }
    if "parallelization" in code_details:
        builder_attribute.parallelization = orm.Dict(
            dict=code_details["parallelization"]
        )


def update_resources(builder, codes):
    """
    Update resources and parallelization settings for various components of the builder.

    Parameters:
        builder: The main builder object to be updated.
        codes: A dictionary containing the configuration codes for different components.
    """
    # Update resources for 'scf' and 'nscf' stages using the 'pw' code details
    if "pw" in codes:
        set_resources(builder.pdos.scf.pw, codes["pw"])
        set_resources(builder.pdos.nscf.pw, codes["pw"])
        set_resources(builder.bands.scf.pw, codes["pw"])
        set_resources(builder.bands.bands.pw, codes["pw"])

    # Update resources for 'dos' stage using the 'dos' code details
    if "soc_dos" in codes:
        set_resources(builder.pdos.dos, codes["dos"])

    # Update resources for 'projwfc' stage using the 'projwfc' code details
    if "soc_projwfc" in codes:
        set_resources(builder.pdos.projwfc, codes["projwfc"])


def get_builder(codes, structure, parameters):
    SOC_PROPERTIES = {
        "bands": ["bands"],
        "pdos": ["pdos"],
        "bands_pdos": ["bands", "pdos"],
    }

    from copy import deepcopy

    protocol = parameters["workchain"].pop("protocol", "moderate")
    pseudos = parameters["advanced"].get("pseudo_family")
    pseudo_info = pseudos.split("/")
    functional = pseudo_info[2]
    pw_code = codes.get("pw")["code"]
    dos_code = codes.get("soc_dos")["code"]
    projwfc_code = codes.get("soc_projwfc")["code"]

    # scf overrides
    scf_overrides = deepcopy(parameters["advanced"])
    scf_overrides["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = parameters["soc"][
        "soc_ecutwfc"
    ]
    scf_overrides["pw"]["parameters"]["SYSTEM"]["ecutrho"] = parameters["soc"][
        "soc_ecutrho"
    ]

    # bands overrides
    bands_overrides = deepcopy(parameters["advanced"])
    bands_overrides.pop("kpoints_distance", None)
    bands_overrides["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = parameters["soc"][
        "soc_ecutwfc"
    ]
    bands_overrides["pw"]["parameters"]["SYSTEM"]["ecutrho"] = parameters["soc"][
        "soc_ecutrho"
    ]

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
    nscf_overrides["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = parameters["soc"][
        "soc_ecutwfc"
    ]
    nscf_overrides["pw"]["parameters"]["SYSTEM"]["ecutrho"] = parameters["soc"][
        "soc_ecutrho"
    ]

    dos_overrides = {"parameters": {"DOS": {}}}
    projwfc_overrides = {"parameters": {"PROJWFC": {}}}

    dos_overrides["parameters"]["DOS"] = {"degauss": parameters["soc"]["pdos_smearing"]}
    projwfc_overrides["parameters"]["PROJWFC"] = {
        "degauss": parameters["soc"]["pdos_smearing"]
    }

    overrides = {
        "bands": {
            "scf": scf_overrides,
            "bands": bands_overrides,
        },
        "pdos": {
            "scf": scf_overrides,
            "nscf": nscf_overrides,
            "dos": dos_overrides,
            "projwfc": projwfc_overrides,
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
        functional=functional,
        clean_workdir=False,
    )
    update_resources(builder, codes)
    return builder


workchain_and_builder = {
    "workchain": SOCWorkChain,
    "get_builder": get_builder,
}
