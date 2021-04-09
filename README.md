# Template_Based_Docking_Project

Tools used to generate docking data via template-based docking with FlexX and HyDE.

1) All pdb files are processed to remove compounds that are not biologically relevant;
2) All compounds are optimized with HyDE, optimized compounds that deviate more than 1 A rmsd are discarded;
3) Pairs of compounds are only used as template:compound-to-dock pairs is they are located in the same pocket;
4) Conformers of compounds-to-dock are generated so that the conformation of the MCS between template and compound-to-dock are as similar as possible;
