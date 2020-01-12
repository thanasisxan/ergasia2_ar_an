import numpy
import matplotlib.pyplot as plt

sin = numpy.sin
pi = numpy.pi

# Τα σημεία που επέλεξα από την συναρτηση του ημιτόνου
x_my_sin_points = [0.380799109526036, 0.761598219052071, 1.39626340159546, 2.66559376668225,
                   2.91945983969961, 3.2367924309713, 3.80799109526036,
                   5.20425449685582, 5.83891967939921, 6.09278575241657]

y_my_sin_points = [sin(i) for i in x_my_sin_points]

# Τα 200 δοκιμαστικά σημεία
x_test_points = numpy.linspace(-pi, pi, 200)
y_test_points = [sin(i) for i in x_test_points]


def polyonimiki_prosegisi_newton(x_points, y_points):
    # Επιστρέφει το πολυώνυμο(προσσέγγιση της συναρτησης) που περνάει απο τα σημέια που περνει σαν όρισμα
    # με την μέθοδο του Newton
    # Η μορφη της επιστρεφόμενης μεταβλητής είναι: numpy.poly1d
    # όπου αποτελεί δομή δεδομένων που προσομοιάζει τα πολυώνυμα
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


def coef(x, y):
    # συνάρτηση για τον υπολογισμό των συντελεστών α για την πολυωνυμική προσσέγιση Newton
    n = len(x)
    a = []
    for i in range(n):
        a.append(y[i])

    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            a[i] = (a[i] - a[i - 1]) / (x[i] - x[i - j])

    return numpy.array(a)


def methodos_elaxistwn_tetragwnwn(x_points, y_points):
    # προσέγγιση των σημειών με πολυώνυμο 1ου βαθμού
    # με την μέθοδο των ελαχίστων τετραγώνων
    x_meso = numpy.mean(x_points)
    y_meso = numpy.mean(y_points)

    sum_arith = 0
    sum_paron = 0
    for i in range(len(x_points)):
        # Δημιουργούμε το άθροισμα του αριθμητή
        sum_arith += (x_points[i] - x_meso)*(y_points[i] - y_meso)
        # Δημιουργούμε το άθροισμα του παρονομαστή
        sum_paron += (x_points[i] - x_meso)**2
    a = sum_arith / sum_paron
    b = y_meso - a*x_meso
    eutheia_prοseggisis = numpy.poly1d([a, b])
    return eutheia_prοseggisis


polyonimo_newton = polyonimiki_prosegisi_newton(
    x_my_sin_points, y_my_sin_points)
# print(polyonimo_newton)
# print(polyonimo_newton(pi/2))
# print(sin(pi/2))


elax_tetr_eutheia = methodos_elaxistwn_tetragwnwn(
    x_my_sin_points, y_my_sin_points)

# z = numpy.polyfit(x_my_sin_points, y_my_sin_points, 1)
# print(a, b)
# print(elax_tetr_eutheia)
# print(z)


# Υπολογισμός του σφάλματος για 200 σημεια μεταξύ -π και π
sum_sfalma_newton = 0
sum_sfalma_leastsq = 0
for point in x_test_points:
    interp_newton = polyonimo_newton(point)
    interp_leastsq = elax_tetr_eutheia(point)
    # print(interp_newton)
    # print(interp_leastsq)
    result = sin(point)
    # print(result)
    sfalma_newton = abs(interp_newton - result)
    sfalma_leastsq = (interp_leastsq - result)**2

    sum_sfalma_newton = sum_sfalma_newton + sfalma_newton
    sum_sfalma_leastsq = sum_sfalma_leastsq + sfalma_leastsq

avg_sfalma_newton = sum_sfalma_newton / len(x_test_points)
final_sfalma_leastsq = numpy.sqrt(sum_sfalma_leastsq)

print("Μέσο σφάλμα πολυωνυμικής προσσέγισης Newton:", avg_sfalma_newton)
print("Σφάλμα προσσέγισης με την μέθοδο των ελάχιστων τετραγώνων:",
      final_sfalma_leastsq)


# print(Y_pred)

# plt.scatter(x_test_points, y_test_points)  # actual
# plt.scatter(x_test_points, Y_pred, color='red')
# # plt.plot([min(X), max(X)], [min(Y_pred), max(Y_pred)],
# #  color='red')  # predicted
# plt.show()
