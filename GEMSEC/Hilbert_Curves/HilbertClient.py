from DimensionReduction import HilbertCurve

iterations = 3 # side length is 2**iterations
num_dimensions = 2 # max value is 2**

hilbert_curve = HilbertCurve(iterations, num_dimensions)

for coords in [[0, 0], [0, 1], [1, 1], [1, 0]]:
     dist = hilbert_curve.distance_from_coordinates(coords)
     print(coords, "-->", dist)

# for dist in range(4):
#      coords = hilbert_curve.coordinates_from_distance(dist)
#      print(dist, "-->", coords)