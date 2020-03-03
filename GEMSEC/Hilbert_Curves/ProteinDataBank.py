from DimensionReduction import HilbertCurve
from pydmd import DMD
from pydmd import MrDMD
import matplotlib.pyplot as plt
import numpy as np
import math, csv

# This class lets you convert PDB files to CSV format
class PeptideCSV:

     def __init__(self, pdb_file):
          
          # Check that file is a pdb_file
          path = pdb_file.split('/')
          extension = path[len(path) - 1].split('.')
          if len(extension) < 2 or extension[1] != "pdb":
               raise Exception('Filename should end with ".pdb"')
          if len(extension) > 2:
               raise Exception('Invalid filename')

          # Check that file contains proper PDB format
          with open(pdb_file, 'r') as pdb:
               firstLine = pdb.readline().split()
               length = len(firstLine)
               if length < 1:
                    raise Exception('First line of file should not be empty')
               if firstLine[0] == 'HEADER':
                    protein_name = firstLine[length - 1] # name of peptide
               else:
                    protein_name = extension[0] # name of file

          self.data = []
          self.max_value = -1
          self.input = pdb_file
          self.output = protein_name
          self.hilbertCurve = HilbertCurve(1, 1)

          
     # # # # # # # # # # #
     # Multi- PDB Frame  #
     # # # # # # # # # # #


     # Reads in Hilbert Distances from PDB file and writes output to a csv file
     def pdb_to_hilbert(self):
          self._get_hilbert() # updates self.data[]
          with open(self.output + ".csv", "w") as csv_file:
               writer = csv.writer(csv_file)
               writer.writerows(self.data)

     def _get_hilbert(self):

          min_x = 999999999 # arbitrarily
          min_y = 999999999 # large
          min_z = 999999999 # numbers

          max_x = -1 # arbitrarily
          max_y = -1 # small
          max_z = -1 # numbers

          with open(self.input, 'r') as pdb:
               frame = 0
               atoms = []
               for line in pdb:
                    list = line.split() # split line into list of words by spaces
                    id = list[0] # look at the first word in each line
                    if id == "MODEL":
                         frame += 1                         
                    elif id == 'ATOM':
                         if frame == 1:
                              atom = int(list[1])
                              atoms.append(atom) # add atom number to the list of atoms
                         # elif list[2] == N:


                         # calculate XYZ values
                         x = int(float(list[6]) * 100)
                         y = int(float(list[7]) * 100)
                         z = int(float(list[8]) * 100)

                         # calculate minimum values
                         min_x = min(x, min_x)
                         min_y = min(y, min_y)
                         min_z = min(z, min_z)

                         # calculate maximum values
                         max_x = max(x, max_x)
                         max_y = max(y, max_y)
                         max_z = max(z, max_z)
                    # elif id == 'ENDMDL':
                         # process last model

               self.data.append(atoms) # add each atom number as the first row
               self.max_value = max(max_x - min_x, max_y - min_y, max_z - min_z)
               num_dimensions = 3
               num_iterations = math.ceil(math.log(self.max_value, 2))
               self.hilbert_curve = HilbertCurve(num_iterations, num_dimensions)

               pdb.seek(0) # move pointer back to start
               row = 0
               for line in pdb:
                    list = line.split() # split line into words by spaces
                    id = list[0] # look at the first word in each line
                    if id == 'ATOM':
                         # Convert string -> float -> int -> positive int
                         x = int(float(list[6]) * 100) - min_x
                         y = int(float(list[7]) * 100) - min_y
                         z = int(float(list[8]) * 100) - min_z

                         coords = [x, y, z]
                         dist = self.hilbert_curve.distance_from_coordinates(coords)
                         if (len(self.data) == row):
                              empty = []
                              self.data.append(empty)
                         self.data[row].append(dist)

                    elif id == 'MODEL':
                         row += 1
          pdb.close()
          

     # # # # # # # # # # #
     # Single PDB Frame  #
     # # # # # # # # # # #


     # Reads in Hilbert Distances from PDB file and writes output to a csv file
     def pdb_to_hilbert_single_frame(self):
          self._get_hilbert_single_frame() # updates self.data[]
          with open(self.output + ".csv", "w") as csv_file:
               writer = csv.writer(csv_file)
               writer.writerows(self.data)

     # Update self.data[] to include PDB information including xyz coordinates
     def _get_xyz_single_frame(self):

          self.data = [["FRAME", "ATOM", "RESIDUE", "CHAIN", "SEQUENCE", "X", "Y", "Z"]]

          min_x = 999999999 # arbitrarily
          min_y = 999999999 # large
          min_z = 999999999 # numbers

          max_x = -1 # arbitrarily
          max_y = -1 # small
          max_z = -1 # numbers

          with open(self.input, 'r') as pdb:
               for line in pdb:
                    list = line.split() # split line into list of words by spaces
                    id = list[0] # look at the first word in each line
                    if id == 'ATOM':
                         # calculate XYZ values
                         x = int(float(list[6]) * 100)
                         y = int(float(list[7]) * 100)
                         z = int(float(list[8]) * 100)

                         # calculate minimum values
                         min_x = min(x, min_x)
                         min_y = min(y, min_y)
                         min_z = min(z, min_z)

                         # calculate maximum values
                         max_x = max(x, max_x)
                         max_y = max(y, max_y)
                         max_z = max(z, max_z)

               self.max_value = max(max_x - min_x, max_y - min_y, max_z - min_z)
               

               pdb.seek(0) # move pointer back to start
               frame = 0
               for line in pdb:
                    list = line.split() # split line into words by spaces
                    id = list[0] # look at the first word in each line
                    if id == 'ATOM':
                         # Convert string -> float -> int -> positive int
                         x = int(float(list[6]) * 100) - min_x
                         y = int(float(list[7]) * 100) - min_y
                         z = int(float(list[8]) * 100) - min_z

                         # collect data for current row
                         atom      = list[2]
                         residue   = list[3]
                         chain     = list[4]
                         sequence  = list[5]

                         # write a new value to the data array
                         row = [frame, atom, residue, chain, sequence, x, y, z]
                         self.data.append(row)
                    elif id == 'MODEL':
                         frame += 1
          pdb.close()

     # Updates self.data[] to use hilbert distances instead of xyz coordinates
     def _get_hilbert_single_frame(self):

          self._get_xyz_single_frame()

          # Update the column headers
          self.data[0].remove("X")
          self.data[0].remove("Y")
          self.data[0].remove("Z")
          self.data[0].append("HILBERT")

          # pull data needed for hilbert object
          num_iterations = math.ceil(math.log(self.max_value, 2))
          num_dimensions = 3

          # construct a new Hilbert Curve object
          self.hilbert_curve = HilbertCurve(num_iterations, num_dimensions)
          
          # change all the [x, y, z] coordinates to hilbert distances
          for i in range(1, len(self.data)):
               # Read in each coordinate from self.data[]
               x = self.data[i][5]
               y = self.data[i][6]
               z = self.data[i][7]

               # Delete the values from self.data[]
               self.data[i].remove(x)
               self.data[i].remove(y)
               self.data[i].remove(z)

               # Convert coordinates into hilbert distances
               coords = [x, y, z]
               dist = self.hilbert_curve.distance_from_coordinates(coords)
               
               # Update the values in self.data[]
               self.data[i].append(dist)

     # Dynamic Mode Decomposition
     def run_dmd(self):

          def _plot_data():
               atom_axis, time_axis = np.meshgrid(time, atoms)

               plt.figure(figsize=(7, 8))
               plt.subplot(2, 1, 1)
               plt.title("Original PDB data")
               plt.pcolormesh(time_axis, atom_axis, snapshot_matrix)
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.subplot(2, 1, 2)
               plt.title("Reconstructed with DMD")
               plt.pcolormesh(time_axis, atom_axis, dmd.reconstructed_data.real)
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.show()

          def _plot_modes():
               plt.figure(figsize=(8, 8))
               for mode in dmd.modes.T:
                    plt.plot(atoms, mode)
                    plt.title('Modes')
               plt.show()

          def _plot_dynamics():
               plt.figure(figsize=(8, 8))
               for dynamic in dmd.dynamics:
                    plt.plot(time, dynamic)
                    plt.title('Dynamics')
               plt.show()
          
          def _print_eigs():
               for eig in dmd.eigs:
                    dist = np.abs(eig.imag**2 + eig.real**2 - 1)
                    print("Eigenvalue:", eig, " Distance from unit circle:", dist)
          
          def _plot_eigs():
               dmd.plot_eigs(show_axes=True, show_unit_circle=True)

          def _print_error():
               error = np.linalg.norm(snapshot_matrix - dmd.reconstructed_data)
               print("DMD error:", error)

          def _plot_error():
               plt.pcolormesh(time, atoms, (snapshot_matrix - dmd.reconstructed_data).real)
               plt.colorbar()
               plt.show()

          self._get_hilbert() # updates self.data[]

          snapshot_matrix = np.array(self.data[1:365]).transpose()

          dmd = DMD(svd_rank = 2, tlsq_rank = 2, exact = True, opt = True) # create instance of DMD object
          dmd.fit(snapshot_matrix) # populate the matrix with data
          
          num_atoms = len(self.data[0])
          num_frames = len(snapshot_matrix[0])

          atoms = np.linspace(1, num_atoms, num_atoms)
          time = np.linspace(1, num_frames, num_frames)

          _plot_data()
          # _plot_modes()
          # _plot_dynamics()
          # _print_eigs()
          # _plot_eigs()
          # _print_error()
          # _plot_error()

     def run_mrdmd(self):

          '''  MrDMD builds a tree-like structure that is max_levels deep:
                    if current_level < 2**(self.max_level - 1):
                         current_level += 1
               /opt/anaconda3/lib/python3.7/site-packages/pydmd/mrdmd.py
          '''

          def _plot_data(title, time, atoms, data):
               atom_axis, time_axis = np.meshgrid(atoms, time)

               plt.figure(figsize=(8, 8))
               plt.title(title)
               plt.pcolormesh(atom_axis, time_axis, data)
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.show()

          def _print_eigs(dmd):
               print('There are', dmd.eigs.shape[0], 'eigenvalues.')

          def _plot_eigs(dmd):
               dmd.plot_eigs(show_axes=True, show_unit_circle=True, figsize=(8, 8))

          def _plot_partial_modes(dmd, atoms, level):
               partial_modes = dmd.partial_modes(level)
               plt.plot(atoms, partial_modes.real)
               plt.show()

          def _plot_partial_dynamics(dmd, time, level):
               partial_dynamics = dmd.partial_dynamics(level)
               plt.plot(time, partial_dynamics.real.T)
               plt.show()

          def _plot_all_levels(dmd, levels, atoms, time):
               for i in range(1, levels + 1):
                    partial_data = dmd.partial_reconstructed_data(level=i)
                    for j in range(i):
                         partial_data += dmd.partial_reconstructed_data(level=j)
                    _plot_data("DMD Levels 0-" + str(i), time, atoms, partial_data.real.T)

          def _plot_side_by_side(dmd, time, atoms, snapshots):
               plt.figure(figsize=(8, 7))
               plt.subplot(1, 2, 1)
               plt.title("Original PDB data")
               plt.pcolormesh(atoms, time, snapshots.T, cmap='viridis')
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.subplot(1, 2, 2)
               plt.title("Reconstructed with MrDMD")
               plt.pcolormesh(atoms, time, dmd.reconstructed_data.real.T, cmap='viridis')
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.show()

          def _plot_side_by_side_new(dmd, time, atoms, snapshots):
               fig, axs = plt.subplots(1, 2)
               plots = [snapshots.T, dmd.reconstructed_data.real.T]
               for col in range(2):
                    ax = axs[col]
                    pcm = ax.pcolormesh(plots[col], cmap='viridis')
               fig.colorbar(pcm, ax=axs[col])
               plt.show()

          def _run_main():
               self._get_hilbert() # updates self.data[]

               snapshots = np.array(self.data[1:]).transpose() # first 18 nanoseconds
               num_levels = int(np.floor(np.log2(snapshots.shape[1]/8))) + 1 # calc from mrdmd.py
               num_atoms = len(self.data[0])
               num_frames = len(snapshots[0])
               atoms = np.linspace(1, num_atoms, num_atoms)
               time = np.linspace(1, num_frames, num_frames)

               # first_dmd = DMD(svd_rank=-1) # Prints the original data with 
               # first_dmd.fit(snapshots) # a basic DMD algorithm

               dmd = MrDMD(svd_rank=-1, max_level=num_levels, max_cycles=1)
               dmd.fit(snapshots.astype('float64'))

               # _plot_side_by_side(dmd, time, atoms, snapshots)
               # _plot_side_by_side_new(dmd, time, atoms, snapshots)
               # _plot_data("Reconstructed MrDMD", time, atoms, dmd.reconstructed_data.real)
               # _print_eigs(dmd)
               _plot_eigs(dmd)
               # _plot_partial_modes(dmd, atoms, 0)
               # _plot_partial_dynamics(dmd, time, 2)
               # _plot_all_levels(dmd, num_levels, atoms, time)
          
          _run_main()