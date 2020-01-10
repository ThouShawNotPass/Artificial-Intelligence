# Artificial Intelligence
A repository to organize and track my neural network algorithms.

# TODO
- Convert dataset units from floats to integers (10x, 100x, 1000x)
- Instead of outputting a number between 0 and n^2 - 1, we could constrict the range of values of d to be between 0 and 1 and express the length as a percentage. Then, as you scale up the resolution of the hilbert space, you would approach the same percentage. (UPDATE: current implementation scales resolution by expanding the size of the data-space, rather than increasing the resolution of the current size)


# References

This module is based on the C code provided in the 2004 article "Programming the Hilbert Curve" by John Skilling found here: http://adsabs.harvard.edu/abs/2004AIPC..707..381S

That code was translated from C into python by github user galtay. Their repository can be found here: https://github.com/galtay/hilbertcurve/blob/develop/hilbertcurve/hilbertcurve.py

There is also an interesting discussion on Stack Overflow about dimension reduction with Hilbert Curves
here: http://stackoverflow.com/questions/499166/mapping-n-dimensional-value-to-a-point-on-hilbert-curve
