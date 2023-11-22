from subprocess import run

functionals = ['PBE', 'PBEsol']
protocols = ['standard', 'stringent']

for functional in functionals:
    for protocol in protocols:
        run(
            [
                "aiida-pseudo",
                "install",
                "pseudo-dojo",
                "--relativistic",
                "FR",
                "--pseudo-format",
                "upf",
                "--version",
                "0.4",
                "--protocol",
                protocol,
                "--functional",
                functional,
            ],
            check=True,
            capture_output=True,
        )