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


print(olokliroma_simpson(sin, 0, pi / 2, 12))
mypoints = [0.0158666295635848, 0.237999443453772, 0.333199220835281, 0.571198664289053, 0.650531812106977,
            0.809198107742826, 0.951997773815089, 1.12653069901452, 1.33279688334112, 1.38039677203188, 1.4597299198498]
print(custom_olokliroma_simpson(sin, 0, pi / 2, mypoints))
