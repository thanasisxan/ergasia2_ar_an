import numpy
import matplotlib.pyplot as plt
import pandas as pd

sin = numpy.sin
pi = numpy.pi

# Τα σημεία που επέλεξα από την συναρτηση του ημιτόνου
x_sin_points = [0.380799109526036, 0.761598219052071, 1.39626340159546, 2.66559376668225,
                2.91945983969961, 3.2367924309713, 3.80799109526036,
                5.20425449685582, 5.83891967939921, 6.09278575241657]

y_sin_points = [sin(i) for i in x_sin_points]


# Επιστρέφει το πολυώνυμο(προσσέγγιση της συναρτησης) που περνάει απο τα σημέια που περνει σαν όρισμα
# με την μέθοδο του Newton
# Η μορφη της επιστρεφόμενης μεταβλητής είναι: numpy.poly1d
# όπου αποτελεί δομή δεδομένων που προσομοιάζει τα πολυώνυμα
def polyonimiki_prosegisi_newton(x_points, y_points):
    assert (len(x_points) == len(y_points)
            ), "Το πλήθος των σημείων πρέπει να είναι ίδιο."

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


# συνάρτηση για τον υπολογισμό των συντελεστών α
def coef(x, y):
    n = len(x)
    a = []
    for i in range(n):
        a.append(y[i])

    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            a[i] = (a[i] - a[i - 1]) / (x[i] - x[i - j])

    return numpy.array(a)


polyonimo_newton = polyonimiki_prosegisi_newton(x_sin_points, y_sin_points)
# print(polyonimo_newton)
# print(polyonimo_newton(pi/2))
# print(sin(pi/2))
points = numpy.linspace(-pi, pi, 200)
sum_sfalma = 0
for point in points:
    interp_result = polyonimo_newton(point)
    # print(interp_result)
    result = sin(point)
    # print(result)
    sfalma = abs(interp_result - result)
    # print("Σφάλμα προσέγγισης:", sfalma)
    sum_sfalma = sum_sfalma + sfalma

avg_sfalma = sum_sfalma / len(points)


print("Μέσο σφάλμα πολυωνυμικής προσσέγισης Newton:", avg_sfalma)


plt.rcParams['figure.figsize'] = (12.0, 9.0)

# Preprocessing Input data
# data = pd.read_csv('data.csv')
# data = [[i, sin(i)] for i in x_sin_points]

X = x_sin_points
Y = y_sin_points

# Building the model
X_mean = numpy.mean(X)
Y_mean = numpy.mean(Y)

num = 0
den = 0
for i in range(len(X)):
    num += (X[i] - X_mean)*(Y[i] - Y_mean)
    den += (X[i] - X_mean)**2
m = num / den
c = Y_mean - m*X_mean

print(m, c)

Y_pred = []
# Making predictions
for i in range(len(X)):
    Y_pred.append(m*X[i]+c)
# Y_pred = m*X + c

plt.plot(X, Y)  # actual
plt.plot(X, Y_pred, color='red')
# plt.plot([min(X), max(X)], [min(Y_pred), max(Y_pred)],
        #  color='red')  # predicted
plt.show()
