from DimensionReduction import HilbertCurve
from pydmd import DMD
from pydmd import MrDMD
import matplotlib.pyplot as plt
import numpy as np
import math, csv

class PDB():
     MAX_VALUE    =  999999999
     MIN_VALUE    = -999999999
     ANGSTROM_PRECISION  = 624
     NUM_DIMENSIONS      = 3
     NUM_ATOMS           = 182

# This class lets you convert PDB files to CSV format
class PeptideCSV():

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
          self.max_value = PDB.MIN_VALUE
          self.min_value = PDB.MAX_VALUE
          self.input = pdb_file
          self.output = protein_name
          self.num_iterations = 1
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

          min_x = PDB.MAX_VALUE # arbitrarily
          min_y = PDB.MAX_VALUE # large
          min_z = PDB.MAX_VALUE # numbers

          max_x = PDB.MIN_VALUE # arbitrarily
          max_y = PDB.MIN_VALUE # small
          max_z = PDB.MIN_VALUE # numbers

          with open(self.input, 'r') as pdb:
               # Step 1: go through each frame, make the n-terminus 0 and shift all values
               # Step 2: loop through the whole file and find the gobal minimum (X, Y, Z)
               # Step 3: adjust every (X, Y, Z) to the global minimum and convert to hilbert
               
               frame     = 0
               frame_x   = 0
               frame_y   = 0
               frame_z   = 0
               atoms     = []
               for line in pdb:
                    words = line.split() # split line into list of words by spaces
                    id = words[0] # look at the first word in each line
                    if id == "MODEL":
                         frame += 1                         
                    elif id == 'ATOM':
                         # calculate XYZ values
                         x = int(float(words[6]) * PDB.ANGSTROM_PRECISION)
                         y = int(float(words[7]) * PDB.ANGSTROM_PRECISION)
                         z = int(float(words[8]) * PDB.ANGSTROM_PRECISION)

                         # record position of the n-terminus for each frame
                         if words[1] == '1':
                              self.data.append([]) # add a new frame entry
                              frame_x = x
                              frame_y = y
                              frame_z = z

                         # shift values in frame
                         x -= frame_x
                         y -= frame_y
                         z -= frame_z

                         # store (x, y, z) position in data structure
                         self.data[frame - 1].append([x, y, z])

                         # calculate minimum values
                         min_x = min(x, min_x)
                         min_y = min(y, min_y)
                         min_z = min(z, min_z)

                         # calculate maximum values
                         max_x = max(x, max_x)
                         max_y = max(y, max_y)
                         max_z = max(z, max_z)

               self.min_value = min(min_x, min_y, min_z)
               self.max_value = max(max_x - self.min_value, max_y - self.min_value, max_z - self.min_value)
               self.num_iterations = math.ceil(math.log(self.max_value, 2))
               self.hilbert_curve = HilbertCurve(self.num_iterations, PDB.NUM_DIMENSIONS)

               pdb.seek(0) # move pointer back to start
               for snapshot in self.data:
                    for i in range(PDB.NUM_ATOMS):
                         for j in range(3):
                              snapshot[i][j] -= self.min_value
                         snapshot[i] = self.hilbert_curve.distance_from_coordinates(snapshot[i])
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

          min_x = PDB.MAX_VALUE # arbitrarily
          min_y = PDB.MAX_VALUE # large
          min_z = PDB.MAX_VALUE # numbers

          max_x = PDB.MIN_VALUE # arbitrarily
          max_y = PDB.MIN_VALUE # small
          max_z = PDB.MIN_VALUE # numbers

          with open(self.input, 'r') as pdb:
               for line in pdb:
                    list = line.split() # split line into list of words by spaces
                    id = list[0] # look at the first word in each line
                    if id == 'ATOM':
                         # calculate XYZ values
                         x = int(float(list[6]) * PDB.ANGSTROM_PRECISION)
                         y = int(float(list[7]) * PDB.ANGSTROM_PRECISION)
                         z = int(float(list[8]) * PDB.ANGSTROM_PRECISION)

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
                         x = int(float(list[6]) * PDB.ANGSTROM_PRECISION) - min_x
                         y = int(float(list[7]) * PDB.ANGSTROM_PRECISION) - min_y
                         z = int(float(list[8]) * PDB.ANGSTROM_PRECISION) - min_z

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

          def _plot_future_state():
               print("Shape before manipulation: {}".format(dmd.reconstructed_data.shape))
               dmd.dmd_time['tend'] *= 40
               print("Shape after manipulation: {}".format(dmd.reconstructed_data.shape))
               new_num_frames = dmd.reconstructed_data.shape[1]
               new_time = np.linspace(1, new_num_frames, new_num_frames)

               atom_axis, time_axis = np.meshgrid(new_time, atoms)
               plt.figure(figsize=(7, 8))
               plt.title("Projection with DMD")
               plt.pcolormesh(time_axis, atom_axis, dmd.reconstructed_data.real)
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.show()

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
               # error = np.linalg.norm((snapshot_matrix - dmd.reconstructed_data))
               error = (np.square((snapshot_matrix - dmd.reconstructed_data).real)).mean(axis=None)
               print("DMD error:", error)

          def _plot_error():
               plt.pcolormesh(time, atoms, np.divide((snapshot_matrix - dmd.reconstructed_data).real, snapshot_matrix))
               plt.colorbar()
               plt.show()

          self.data = [] # TODO: DANGEROUS PLEASE REMOVE (TESTING ONLY)
          self._get_hilbert() # updates self.data[]

          snapshot_matrix = np.array(self.data).transpose()

          dmd = DMD(svd_rank = .97, tlsq_rank = 2, exact = True, opt = True) # create instance of DMD object
          dmd.fit(snapshot_matrix) # populate the matrix with data
          
          num_atoms = len(self.data[0])
          num_frames = len(snapshot_matrix[0])

          atoms = np.linspace(1, num_atoms, num_atoms)
          time = np.linspace(1, num_frames, num_frames)

          # _plot_future_state()
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

          def _print_error(dmd, snapshots):
               # error = np.linalg.norm(snapshots - dmd.reconstructed_data)
               # error = (np.square((snapshots - dmd.reconstructed_data).real)).mean(axis=None)
               error = (np.square((snapshots[0] - dmd.reconstructed_data[0]).real)).mean(axis=None)
               print("MrDMD error:", error)

          def _plot_error(snapshots, dmd, time, atoms):
               error = np.divide(abs(snapshots - dmd.reconstructed_data.real), snapshots)
               plt.pcolormesh(time, atoms, error.T)
               plt.colorbar()
               plt.show()

          def _plot_partial_modes(dmd, atoms, level):
               partial_modes = dmd.partial_modes(level)
               plt.plot(atoms, partial_modes.real)
               plt.show()

          def _plot_partial_dynamics(dmd, time, level):
               partial_dynamics = dmd.partial_dynamics(level)
               plt.plot(time, partial_dynamics.real)
               plt.show()

          def _plot_all_levels(dmd, levels, atoms, time):
               for i in range(1, levels + 1):
                    partial_data = dmd.partial_reconstructed_data(level=i)
                    for j in range(i):
                         partial_data += dmd.partial_reconstructed_data(level=j)
                    _plot_data("DMD Levels 0-" + str(i), time, atoms, partial_data.real.T)

          def _plot_side_by_side(reconstructed_data, time, atoms, snapshots):
               plt.figure(figsize=(8, 7))
               plt.subplot(1, 2, 1)
               plt.title("Original PDB data")
               plt.pcolormesh(atoms, time, snapshots, cmap='viridis')
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.subplot(1, 2, 2)
               plt.title("Reconstructed with MrDMD")
               plt.pcolormesh(atoms, time, reconstructed_data, cmap='viridis')
               plt.xlabel("Atom Index")
               plt.ylabel("Frame")
               plt.colorbar()
               plt.show()

          def _plot_side_by_side_new(reconstructed_data, time, atoms, snapshots):
               fig, axs = plt.subplots(1, 2)
               plots = [snapshots, reconstructed_data]
               for col in range(2):
                    ax = axs[col]
                    pcm = ax.pcolormesh(plots[col], cmap='viridis')
               fig.colorbar(pcm, ax=axs[col])
               plt.show()

          def _find_xyz_dist(reconstructed_data, snapshots):
               distances = []
               frames = len(snapshots)
               atoms = len(snapshots[0])
               max_hilbert = 2**(self.num_iterations * PDB.NUM_DIMENSIONS) - 1
               for i in range(frames): # each frame
                    distances.append([])
                    for j in range(atoms): #each atom
                         actual = self.hilbert_curve.coordinates_from_distance(int(snapshots[i][j]))
                         predicted_hilbert = min(max_hilbert, max(0, int(np.rint(reconstructed_data[i][j]))))
                         predicted = self.hilbert_curve.coordinates_from_distance(predicted_hilbert)
                         sum_of_squares = (((actual[0]-predicted[0])**2)+((actual[1]-predicted[1])**2)+((actual[2]-predicted[2])**2))
                         distances[i].append(sum_of_squares**(1/2) / PDB.ANGSTROM_PRECISION)
                    frame_mean = np.mean(distances[i])
                    # print('Average Distance for Frame', i + 1, 'was', frame_mean)
                    distances[i] = frame_mean
               # total_mean = np.mean(distances)
               # print('Level', PDB.ANGSTROM_PRECISION, 'had an average distance of', total_mean)
               return np.mean(distances)

          def _truncate_dmd(dmd):
               result = dmd.reconstructed_data.real
               frames = len(result)
               atoms = len(result[0])
               max_hilbert = 2**(self.num_iterations * PDB.NUM_DIMENSIONS) - 1
               for i in range(frames): # each frame
                    for j in range(atoms): #each atom
                         result[i][j] = min(max_hilbert, max(0, int(np.rint(result[i][j]))))
               return result
          
          def _run_main():
               self.data = [] # TODO: DANGEROUS PLEASE REMOVE (TESTING ONLY)
               self._get_hilbert() # updates self.data[]

               snapshots = np.array(self.data)
               num_levels = int(np.floor(np.log2(snapshots.shape[1]/8))) + 1 # calc from mrdmd.py
               num_atoms = len(snapshots[0])
               num_frames = len(snapshots)
               atoms = np.linspace(1, num_atoms, num_atoms)
               time = np.linspace(1, num_frames, num_frames)

               # first_dmd = DMD(svd_rank=-1) # Prints the original data with 
               # first_dmd.fit(snapshots) # a basic DMD algorithm

               dmd = MrDMD(svd_rank=-1, max_level=num_levels, max_cycles=1)
               dmd.fit(snapshots.astype('float64'))
               reconstructed_data = _truncate_dmd(dmd)

               # return _find_xyz_dist(reconstructed_data, snapshots)
               # _plot_side_by_side(reconstructed_data, time, atoms, snapshots)
               _plot_side_by_side_new(reconstructed_data, time, atoms, snapshots)
               # _plot_data("Reconstructed MrDMD", time, atoms, reconstructed_data)
               # _print_eigs(dmd)
               # _plot_eigs(dmd)
               # print('Mean Error:', _find_xyz_dist(reconstructed_data, snapshots), 'Angstroms')
               # _plot_error(snapshots, dmd, time, atoms)
               # _plot_partial_modes(dmd, atoms, 0)
               # _plot_partial_dynamics(dmd, time, 2)
               # _plot_all_levels(dmd, num_levels, atoms, time)
          # return _run_main()
          _run_main()

     def run_test(self):
          min_dist = PDB.MAX_VALUE
          best_resolution = 1
          for i in range(1000):
               dist = self.run_mrdmd()
               if (dist < min_dist):
                    min_dist = dist
                    best_resolution = i
               PDB.ANGSTROM_PRECISION += 1
               print('Testing resolution', i)
               print('   d =', dist)
          print('Minimum distance:', min_dist)
          print('Found at resolution:', best_resolution)