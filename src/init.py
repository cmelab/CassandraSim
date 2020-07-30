#!/usr/bin/env python
"""Initialize the project's data space.
Iterates over all defined state points and initializes
the associated job workspace directories."""

import signac
def main():
    project = signac.init_project("MosCas")
    for T in [ 66, 151, 236, 321, 406, 492, 577, 662, 747, 833]:
#        for chem_pot in [299]:
        for chem_pot in [0]:
                statepoint = dict(
                    chem_pot=chem_pot, T=T,
                )
                job = project.open_job(statepoint)
                job.doc.setdefault("steps_run", int(22000000)) #21603294 3600549 10000
                job.doc.setdefault("steps_restart", int(29000000)) # 28916707 24097256 50000
                job.init()

if __name__ == "__main__":
    main()