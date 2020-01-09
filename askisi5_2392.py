import numpy
import matplotlib.pyplot as plt

# Τα σημεία που επέλεξα
x_sin_points = [0.380799109526036, 0.761598219052071, 1.39626340159546, 2.66559376668225,
                2.91945983969961, 3.2367924309713, 3.80799109526036,
                5.20425449685582, 5.83891967939921, 6.09278575241657]

y_sin_points = [0.371662455660328, 0.690079011482112, 0.984807753012208, 0.458226521727411,
                0.220310532786541, -0.0950560433041826, -0.618158986220605,
                -0.881453363447582, -0.429794912089172, -0.189251244360411]


def polyonimiki_prosegisi(x_points, y_points):
    # πλήθος των σημείων
    assert (len(x_points) == len(y_points)), "Το πλήθος των σημέιων πρέπει να είναι ίδιο."
    k = len(x_points)
    x = numpy.poly1d([1, 0])
    print(x)
    N = 0
    for j in range(0, k):
        # print("j", j)
        nj = 1
        for i in range(0, j - 1):
            # print("i", i)
            nj *= numpy.polyadd(x, -x_points[i])


polyonimiki_prosegisi(x_sin_points, y_sin_points)

# p1 = numpy.poly1d([3, 2, 3])
# p2 = numpy.poly1d([2, 4, 19])
# print(numpy.polyadd(p1, p2))

# plt.plot(p1)
# plt.show()
