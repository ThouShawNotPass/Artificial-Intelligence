from DimensionReduction import HilbertCurve

# Notes: Only works with integers and scales up/right rather than increasing detail.
#    The smallest "step" will always be a the integer 1. Need to scale up data 100x.

# Setup:
iterations = 2 # side length is 2**iterations 
num_dimensions = 2

hilbert_curve = HilbertCurve(iterations, num_dimensions)

for coords in [[0, 0], [0, 1], [1, 1], [1, 0]]:
     dist = hilbert_curve.distance_from_coordinates(coords)
     print(coords, "-->", dist)