from DimensionReduction import HilbertCurve
from ProteinDataBank import PeptideCSV

pdb_big = '/Users/Justin/Desktop/GitHub/AI/Artificial-Intelligence/GEMSEC/PDB_Files/3GFT.pdb'
pdb_small = '/Users/Justin/Desktop/GitHub/AI/Artificial-Intelligence/GEMSEC/PDB_Files/1A3I.pdb'

csv_file = PeptideCSV(pdb_big)
csv_file.pdb_to_hilbert()