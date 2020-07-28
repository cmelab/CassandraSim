#!/usr/bin/env python
"""Initialize the project's data space.
Iterates over all defined state points and initializes
the associated job workspace directories."""

import signac
def main():
    project = signac.init_project("MosCas")
    for T in [66]:
        for chem_pot in [299]:
#        for chem_pot in [0,33,66,99,133,166,199,233,266,299]:
                statepoint = dict(
                    chem_pot=chem_pot, T=T,
                )
                job = project.open_job(statepoint)
                job.doc.setdefault("steps_run", int(10000)) #3600549 10000
                job.doc.setdefault("steps_restart", int(50000)) #24097256 50000
                job.init()

if __name__ == "__main__":
    main()