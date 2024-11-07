The folder [simulation_scripts](simulation_scripts) contains the scripts to do the simulations (generate_simu_DATASET.R), as well as the script to generate a reference with an additional artificial cell type.
The simulations scripts will load a reference from a folder "references/" and create simulations with the date of the day in a folder called "simulations/"

The reference matrices used for the *in silico* data can be downloaded from Zenodo (DOI 10.5281/zenodo.14024479).

The table below recapitulates where to download *in vitro* and *in vivo* data.

| **Dataset name** | **Source** | **Omic** | **Expression matrix**                                      | **Proportion matrix**                                           | **Reference profiles matrix**                              |
|------------------|------------|----------|------------------------------------------------------------|-----------------------------------------------------------------|------------------------------------------------------------|
| **BlMIX**        | In vitro   | DNAm     | GSE77797                                                   | Figure 4a (DOI: 10.1186/s12859-016-0943-7)                      | GSE35069                                                   |
| **BrMIX**        | In vitro   | RNA      | GSE220605                                                  | "Cell count" column, table S1 (DOI: 10.1186/s13059-023-03016-6) | GSE220605                                                  |
| **PaMIX**        | In vitro   | DNAm     | GSE281305                                         | TODO                                                            | GSE281305                                                                      |
| **PaMIX**        | In vitro   | RNA      | GSE281204                                                  | TODO                                                            | GSE281204                                                  |
| **BlREAL1**      | In vivo    | RNA      | https://github.com/Honchkrow/Deconer_dataset (Linsley.rds) | https://github.com/Honchkrow/Deconer_dataset (Linsley.rds)      | https://github.com/Honchkrow/Deconer_dataset (Linsley.rds) |
| **BlREAL2**      | In vivo    | DNAm     | GSE35069                                                   | DOI: 10.1371/journal.pone.0041361.s004                          | GSE35069                                                   |
| **BlREAL3**      | In vivo    | DNAm     | GSE42861                                                   | Table S2 (DOI: 10.1038/nbt.2487)                                | GSE35069                                                   |
| **BlREAL3**      | In vivo    | RNA      | GSE93722                                                   | Supplementary file 3A (DOI: 10.7554/eLife.26476 )               | EPIC::TRef                                                 |


