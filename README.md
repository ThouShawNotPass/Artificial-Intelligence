<p align="center">
  <a href="http://www.gemsec.washington.edu/" target="_blank" >
    <img alt="GEMSEC" src="GEMSEC/img/gemsec_logo.png" width="400" />
  </a>
</p>
<p align="center">
     <a href="https://github.com/ThouShawNotPass/Artificial-Intelligence/blob/master/LICENSE" target="_blank">
        <img alt="Software License" src="https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square">
    </a>
</p>


# Hilbert Curve - Dimension Reduction
A python library for converting PDB files into CSV files where XYZ coordinates can be transformed into Hilbert Distances.

# Get Started

To start using this package, you need to import the HilbertCurve class from the file DimensionReduction.py (be sure they are in the same directory). To construct a HilbertCurve object, you will need to pass two integers to the constructor; the number of iterations, and the number of dimensions.

With the package installed, you can:
1. Convert XYZ coordinates to hilbert distances

          from DimensionReduction import HilbertCurve

          num_iterations = 3
          num_dimensions = 2
          
          hilbert_curve = HilbertCurve(num_iterations, num_dimensions)

          for coords in [[0, 0], [0, 1], [1, 1], [1, 0]]:
               dist = hilbert_curve.distance_from_coordinates(coords)
               print(coords, "-->", dist)

     Output:

          [0, 0] --> 0
          [0, 1] --> 1
          [1, 1] --> 2
          [1, 0] --> 3

2. Calculate coordinates given hilbert distances:

          from DimensionReduction import HilbertCurve

          num_iterations = 3
          num_dimensions = 2

          hilbert_curve = HilbertCurve(num_iterations, num_dimensions)

          for dist in range(4):
               coords = hilbert_curve.coordinates_from_distance(dist)
               print(dist, "-->", coords)

     Output:

          0 --> [0, 0]
          1 --> [0, 1]
          2 --> [1, 1]
          3 --> [1, 0]

# How it works
A typical PDB file describing a protein consists of hundreds to thousands of lines like the [following](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)):

     HEADER    EXTRACELLULAR MATRIX                    22-JAN-98   1A3I
     ...
     ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
     ATOM      2  CA  PRO A   1       7.608  20.729  20.336  1.00 17.44           C
     ATOM      3  C   PRO A   1       8.487  20.707  19.092  1.00 17.44           C
     ATOM      4  O   PRO A   1       9.466  21.457  19.005  1.00 17.44           O
     ATOM      5  CB  PRO A   1       6.460  21.723  20.211  1.00 22.26           C
     ...

The program will create a new CSV file based on the contents of the original PDB file. For example, a CSV file called ```1A3I.CSV``` would be created containing the following rows and columns, based on the official [WWPDB Documentation](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM):

     FRAME    ATOM    RESIDUE    CHAIN    SEQUENCE    HILBERT
     1        N       PRO        A        1           12545533786
     1        CA      PRO        A        1           4421809923
     1        C       PRO        A        1           3987905106
     1        O       PRO        A        1           65327816114
     1        CB      PRO        A        1           8015651838

# TODO

- None.

# Library

- [Introduction to Statistical Learning](http://faculty.marshall.usc.edu/gareth-james/ISL/ISLR%20Seventh%20Printing.pdf)

- [Elements of Statistical Learning](https://web.stanford.edu/~hastie/Papers/ESLII.pdf)

- Data Driven Science and Engineering (Brunton and Kutz)

# References

This module is based on the C code provided in the 2004 article "Programming the Hilbert Curve" by John Skilling found here: http://adsabs.harvard.edu/abs/2004AIPC..707..381S

That code was translated from C into python by github user galtay. Their repository can be found here: https://github.com/galtay/hilbertcurve/blob/develop/hilbertcurve/hilbertcurve.py

There is also an interesting discussion on Stack Overflow about dimension reduction with Hilbert Curves
here: http://stackoverflow.com/questions/499166/mapping-n-dimensional-value-to-a-point-on-hilbert-curve
