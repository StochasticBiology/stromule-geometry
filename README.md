# stromule-geometry
Simulation and statistics for the advantages of different model and observed stromule structures

`model-many.c` is C simulation code for investigating "interaction region" and "plastid access" for different plastid separations and stromule geometries. It produces several output files, which are plotted by `model-plots.R`.

`Data/` contains raw Excel spreadsheets from biological measurements. These are read and processed by `analyse-stats.R` to produce summary statistics of plastid, cellular, and stromule behaviours.
