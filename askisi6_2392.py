def olok_simpson(func, a, b, N):
    sN2_1 = 0  # το άθροισμα μεχρι Ν/2-1
    sN2 = 0  # το άθροισμα μέχρι Ν/2

    def x_i(k):  # υπολογισμός του όρου Xi=Xo+κ*(b-a)/N όπου Xo=a και κ=0,...,Ν
        return a + k * (b - a) / N

    for i in range(1, N // 2):
        sN2_1 += func(x_i(2 * i))
        sN2 += func(x_i(2 * i - 1))

    # το τελευταίο στοιχείο που πρέπει να προστεθεί στο άθροισμα που παει μέχρι Ν/2
    sN2 += func(x_i(N - 1))

    E = ((b - a) / (3 * N)) * (func(x_i(0)) + func(x_i(N)) + 2 * sN2_1 + 4 * sN2_1)
    return E


def function(x): return x


print(olok_simpson(function, 0.0, 1.0, 100))
