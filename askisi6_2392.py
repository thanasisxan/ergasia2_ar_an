import numpy

sin = numpy.sin
pi = numpy.pi


def olokliroma_simpson(func, a, b, N):
    sN2_1 = 0  # το άθροισμα μεχρι Ν/2-1
    sN2 = 0  # το άθροισμα μέχρι Ν/2

    def x_i(k):  # υπολογισμός του όρου Xi=X0+κ*(b-a)/N όπου X0=a και κ=0,...,Ν
        return a + k * (b - a) / N

    for i in range(1, N // 2):
        # υπολογισμών των αθροισμάτων
        sN2_1 += func(x_i(2 * i))
        sN2 += func(x_i(2 * i - 1))
    # το τελευταίο στοιχείο που πρέπει να προστεθεί στο άθροισμα που παει μέχρι Ν/2
    sN2 += func(x_i(N - 1))

    E = ((b - a) / (3 * N)) * (func(x_i(0)) + func(x_i(N)) + 2 * sN2_1 + 4 * sN2)
    return E


def custom_olokliroma_simpson(func, a, b, points):
    sN2_1 = 0  # το άθροισμα μεχρι Ν/2-1
    sN2 = 0  # το άθροισμα μέχρι Ν/2
    N = len(points) - 1

    def x_i(k):  # επιστρέφει τα σημεία που εμείς επιλέξαμε για να γίνει η μέθοδος simpson
        return points[k]

    for i in range(1, N // 2):
        # υπολογισμών των αθροισμάτων
        sN2_1 += func(x_i(2 * i))
        sN2 += func(x_i(2 * i - 1))
    # το τελευταίο στοιχείο που πρέπει να προστεθεί στο άθροισμα που παει μέχρι Ν/2
    sN2 += func(x_i(N - 1))

    E = ((b - a) / (3 * N)) * (func(x_i(0)) + func(x_i(N)) + 2 * sN2_1 + 4 * sN2)
    return E


def olokliroma_trapezio(func, a, b, N):
    sumT = 0  # το άθροισμα μέχρι Ν-1

    def x_i(k):  # υπολογισμός του όρου Xi=X0+κ*(b-a)/N όπου X0=a και κ=0,...,Ν
        return a + k * (b - a) / N

    for i in range(1, N):
        # υπολογισμών του αθροίσματος
        sumT += func(x_i(i))
    # το τελευταίο στοιχείο που πρέπει να προστεθεί στο άθροισμα που παει μέχρι Ν/2

    E = ((b - a) / (2 * N)) * (func(x_i(0)) + func(x_i(N)) + 2 * sumT)
    return E


def custom_olokliroma_trapezio(func, a, b, points):
    sumT = 0  # το άθροισμα μέχρι Ν-1
    N = len(points) - 1

    def x_i(k):  # επιστρέφει τα σημεία που εμείς επιλέξαμε για να γίνει η μέθοδος του τραπεζίου
        return points[k]

    for i in range(1, N):
        # υπολογισμών του αθροίσματος
        sumT += func(x_i(i))
    # το τελευταίο στοιχείο που πρέπει να προστεθεί στο άθροισμα που παει μέχρι Ν/2

    E = ((b - a) / (2 * N)) * (func(x_i(0)) + func(x_i(N)) + 2 * sumT)
    return E


print("Υπολογισμός με την μέθοδο Simpson (χωρίς δικά μας σημεία)", olokliroma_simpson(sin, 0, pi / 2, 10))
print("Υπολογισμός με την μέθοδο του τραπεζίου (χωρίς δικά μας σημεία)", olokliroma_trapezio(sin, 0, pi / 2, 10))
mypoints = [0.0158666295635848, 0.237999443453772, 0.333199220835281, 0.571198664289053, 0.650531812106977,
            0.809198107742826, 0.951997773815089, 1.12653069901452, 1.33279688334112, 1.38039677203188, 1.4597299198498]
print("Υπολογισμός με την μέθοδο Simpson (με δικά μας σημεία)", custom_olokliroma_simpson(sin, 0, pi / 2, mypoints))
print("Υπολογισμός με την μέθοδο του τραπεζίου (με δικά μας σημεία)",
      custom_olokliroma_trapezio(sin, 0, pi / 2, mypoints))

# Υπολογισμός του σφάλματος
a1 = 0
b1 = pi / 2
N1 = 11
# Η τέταρτη παράγωγος του ημιτόνου ειναι η συνάρτηση του ημιτόνου
# όπου στο πεδίο [0,π/2] εχει μεγιστη τιμη το 1 για π/2
sfalma_trapezio = ((b - a) ** 3) / (12 * N1 ** 2)
sfalma_simpson = ((b - a) ** 5) / (180 * N1 ** 4)

print()
print("Θεωρητικό σφάλμα της μεθόδου τραπεζίου:", sfalma_trapezio)
print("Θεωρητικό σφάλμα της μεθόδου Simpson:", sfalma_simpson)
