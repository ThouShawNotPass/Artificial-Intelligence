from DimensionReduction import HilbertCurve
from ProteinDataBank import PeptideCSV

pdb_big = '/Users/Justin/Desktop/GitHub/AI/Artificial-Intelligence/GEMSEC/PDB_Files/WT_295K_200ns_50ps_0_run.pdb'
pdb_small = '/Users/Justin/Desktop/GitHub/AI/Artificial-Intelligence/GEMSEC/PDB_Files/3GFT.pdb'
pdb_tiny = '/Users/Justin/Desktop/GitHub/AI/Artificial-Intelligence/GEMSEC/PDB_Files/1A3I.pdb'
test = '/Users/Justin/Desktop/GitHub/AI/Artificial-Intelligence/GEMSEC/PDB_Files/365_frame_test.pdb'

# csv_big = 'WT_295K_200ns_50ps_0_run.csv'
# csv_small = '3GFT.csv'
# csv_tiny = '1A3I.csv'

file = pdb_big

def pdb_to_csv():
     csv_file = PeptideCSV(file)
     csv_file.pdb_to_hilbert()

def run_dmd():
     csv_file = PeptideCSV(file)
     csv_file.run_dmd()

def run_mrdmd():
     csv_file = PeptideCSV(file)
     csv_file.run_mrdmd()

def run_test():
     csv_file = PeptideCSV(file)
     csv_file.run_test()

# run_test()
run_mrdmd()
# run_dmd()