from DimensionReduction import HilbertCurve
import math
import csv

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
                    protein_name = firstLine[length - 1]
               else:
                    protein_name = extension[0]

          self.data = []
          self.max_value = -1
          self.input = pdb_file
          self.output = protein_name

          
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
                              atoms.append(atom)
                         
                         # calculate XYZ values
                         x = int(float(list[6]) * 1000)
                         y = int(float(list[7]) * 1000)
                         z = int(float(list[8]) * 1000)

                         # calculate minimum values
                         min_x = min(x, min_x)
                         min_y = min(y, min_y)
                         min_z = min(z, min_z)

                         # calculate maximum values
                         max_x = max(x, max_x)
                         max_y = max(y, max_y)
                         max_z = max(z, max_z)

               self.data.append(atoms)
               self.max_value = max(max_x - min_x, max_y - min_y, max_z - min_z)
               num_dimensions = 3
               num_iterations = math.ceil(math.log(self.max_value, 2))
               hilbert_curve = HilbertCurve(num_iterations, num_dimensions)

               pdb.seek(0) # move pointer back to start
               row = 1
               for line in pdb:
                    list = line.split() # split line into words by spaces
                    id = list[0] # look at the first word in each line
                    if id == 'ATOM':
                         # Convert string -> float -> int -> positive int
                         x = int(float(list[6]) * 1000) - min_x
                         y = int(float(list[7]) * 1000) - min_y
                         z = int(float(list[8]) * 1000) - min_z

                         coords = [x, y, z]
                         dist = hilbert_curve.distance_from_coordinates(coords)
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
                         x = int(float(list[6]) * 1000)
                         y = int(float(list[7]) * 1000)
                         z = int(float(list[8]) * 1000)

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
                         x = int(float(list[6]) * 1000) - min_x
                         y = int(float(list[7]) * 1000) - min_y
                         z = int(float(list[8]) * 1000) - min_z

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
          hilbert_curve = HilbertCurve(num_iterations, num_dimensions)
          
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
               dist = hilbert_curve.distance_from_coordinates(coords)
               
               # Update the values in self.data[]
               self.data[i].append(dist)