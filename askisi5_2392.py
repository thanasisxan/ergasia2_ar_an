import numpy
import matplotlib.pyplot as plt

sin = numpy.sin

# Τα σημεία που επέλεξα από την συναρτηση του ημιτόνου
x_sin_points = [0.380799109526036, 0.761598219052071, 1.39626340159546, 2.66559376668225,
                2.91945983969961, 3.2367924309713, 3.80799109526036,
                5.20425449685582, 5.83891967939921, 6.09278575241657]

y_sin_points = [0.371662455660328, 0.690079011482112, 0.984807753012208, 0.458226521727411,
                0.220310532786541, -0.0950560433041826, -0.618158986220605,
                -0.881453363447582, -0.429794912089172, -0.189251244360411]


# Επιστρέφει το πολυώνυμο(προσσέγγιση της συναρτησης) που περνάει απο τα σημέια που περνει σαν όρισμα
# με την μέθοδο του Newton
# Η μορφη της επιστρεφόμενης μεταβλητής είναι: numpy.poly1d
# όπου αποτελεί δομή δεδομένων που προσομοιάζει τα πολυώνυμα
def polyonimiki_prosegisi_newton(x_points, y_points):
    assert (len(x_points) == len(y_points)), "Το πλήθος των σημείων πρέπει να είναι ίδιο."

    # πλήθος των σημείων k
    n = len(x_points)

    # Το πολυώνυμο 1x+0 (μεταβλητή x σε συμβατή μορφη)
    x = numpy.poly1d([1, 0])

    # Διαιρεμένες διαφορές υπολογισμένες στον πίνακα Dij
    Dij = coef(x_points, y_points)

    N = 0
    for j in range(0, n):
        nj = 1
        for i in range(0, j):
            nj = numpy.polymul(numpy.polyadd(x, -x_points[i]), nj)
        N = numpy.polyadd(numpy.polymul(Dij[j], nj), N)
    return N


def coef(x, y):
    n = len(x)
    a = []
    for i in range(n):
        a.append(y[i])

    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            a[i] = (a[i] - a[i - 1]) / (x[i] - x[i - j])

    return numpy.array(a)


point = 2
interp_result = polyonimiki_prosegisi_newton(x_sin_points, y_sin_points)(point)
print(interp_result)
result = sin(point)
print(result)
sfalma = abs(interp_result - result)
print("Σφάλμα προσέγγισης:", sfalma)
