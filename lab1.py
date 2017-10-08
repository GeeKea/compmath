gamma_0 = 5./3
ro_0 = 0.00271
P_0 = 500000.
U_0 = 0.
gamma_3 = 7./5
ro_3 = 0.0002219
U_3 = 158700.
P_3 = 1000000.
C_0 = (gamma_0 * P_0 / ro_0) ** 0.5
C_3 = (gamma_3 * P_3 / ro_3) ** 0.5
X = P_3 / P_0
alpha_0 = (gamma_0 + 1) / (gamma_0 - 1)
alpha_3 = (gamma_3 + 1) / (gamma_3 - 1)
e_3 = 2 * C_3**2 / (gamma_3 * (gamma_3 - 1) * (U_3 - U_0)**2)
e_0 = 2 * C_0**2 / (gamma_0 * (gamma_0 - 1) * (U_3 - U_0)**2)
a_0 = (alpha_0 * e_3 - alpha_3 * X * e_0)**2
a_1 = 2 * ((alpha_0 * e_3 - alpha_3 * X * e_0) * (e_3 * (1 - 2 * alpha_0 * X) - e_0 * X * (X - 2 * alpha_3)) - alpha_0 * alpha_3 * X * (alpha_0 * e_3 + alpha_3 * X * e_0))
a_2 = e_3 ** 2 * (6 * alpha_0 ** 2 * X ** 2 - 8 * alpha_0 * X + 1) - 2 * e_3 * e_0 * X * (alpha_0 * alpha_3 * (X ** 2 + 4 * X + 1) - 2 * (X + 1) * (alpha_3 + alpha_0 * X) + X) + e_0 ** 2 * X ** 2 * (6 * alpha_3 ** 2 - 8 * alpha_3 * X + X ** 2) + alpha_0 ** 2 * alpha_3 ** 2 * X ** 2 - 2 * alpha_0 * X * e_3 * (alpha_0 * X - 2 * alpha_0 * alpha_3 * X + 2 * alpha_3) - 2 * alpha_3 * X ** 2 * e_0 * (alpha_3 + 2 * alpha_0 * X - 2 * alpha_0 * alpha_3)
a_3 = -2 * X * (2 * e_3**2 * (alpha_0**2 * X**2 - 3 * alpha_0 * X + 1) + e_0 * e_3 * ((alpha_3 + alpha_0 * X) * (X**2 + 4 * X + 1) - 2 * alpha_0 * alpha_3 * X * (X + 1) - 2 * X * (X + 1)) + 2 * e_0 ** 2 * X * (X ** 2 - 3 * alpha_3 * X + alpha_3**2) - alpha_0 * alpha_3 * X * (alpha_0 * X + alpha_3) + e_3 * (alpha_0**2 * alpha_3 * X**2 - 2 * X * (2 * alpha_0 * alpha_3 + alpha_0**2 * X) + (2 * alpha_0 * X + alpha_3)) + e_0 * X * (alpha_0 * alpha_3**2 - 2 * alpha_3 * (alpha_3 + 2 * alpha_0 * X) + 2 * alpha_3 * X + alpha_0 * X**2))
a_4 = X**2 * (e_3**2 * (alpha_0**2 * X**2 - 8 * alpha_0 * X + 6) - 2 * e_0 * e_3 * (alpha_0 * alpha_3 * X - 2 * (X + 1) * (alpha_3 + alpha_0 * X) + X**2 + 4 * X + 1) + e_0**2 * (alpha_3**2 - 8 * alpha_3 * X + 6 * X**2) + (alpha_3**2 + 4*alpha_0*alpha_3*X + alpha_0**2*X**2) - 2 * e_3 * ((alpha_0**2 * X + 2*alpha_0*alpha_3) * X - 2 * (2 * alpha_0 * X + alpha_3) + 1) - 2 * e_0 * (alpha_3 * (2 * alpha_0 * X + alpha_3) - 2 * X * (2 * alpha_3 + alpha_0 * X) + X**2))
a_5 = 2 * X**3 * (e_3**2 * (alpha_0 * X - 2) - e_0 * e_3 * (alpha_0 * X - 2 + alpha_3 - 2 * X) + e_0**2 * (alpha_3 - 2 * X) + (alpha_3 + alpha_0 * X) - e_3 * (2 * alpha_0 * X + alpha_3 - 2) - e_0 * (2 * alpha_3 + alpha_0 * X - 2 * X))
a_6 = X**4 * ((e_3 - e_0)**2 + 1 - 2 * (e_3 + e_0))
coefs = [a_0, a_1, a_2, a_3, a_4, a_5, a_6]
A = max(coefs[1:], key=abs)
B = max(coefs[:-1], key=abs)
R_1 = abs(a_6) / (abs(a_6) + B)
R_2 = 1 + A / abs(a_0)


def f(x, i):
    if i == 0:
        return a_0*x**6 + a_1*x**5 + a_2*x**4 + a_3*x**3 + a_4*x**2 + a_5*x + a_6
    if i == 1:
        return 6*a_0*x**5 + 5*a_1*x**4 + 4*a_2*x**3 + 3*a_3*x**2 + 2*a_4*x + a_5
    if i == 2:
        return 30*a_0*x**4 + 20*a_1*x**3 + 12*a_2*x**2 + 6*a_3*x + 2*a_4
    if i == 3:
        return 120*a_0*x**3 + 60*a_1*x**2 + 24*a_2*x + 6*a_3
    if i == 4:
        return 360*a_0*x**2 + 120*a_1*x + 24*a_2
    if i == 5:
        return 720*a_0*x + 120*a_1
    if i == 6:
        return 720*a_0


def budan_fourier(a, b):
    s_a = 0
    s_b = 0
    for i in range(0, 6):
        if f(a, i) != 0:
            if f(a, i) * f(a, i+1) < 0:
                s_a = s_a + 1
              #  print(a, i, f(a, i), f(a, i + 1), s_a)
        else:
            return -1
        if f(b, i) != 0:
            if f(b, i) * f(b, i+1) < 0:
                s_b = s_b + 1
              #  print(b, i, f(b, i), f(b, i + 1), s_b)
        else:
            return -2
    return s_a - s_b


def newton(a, b):
    x_0 = (a + b)/2
    if f(x_0, 0)*f(x_0, 1) > 0:
        while True:
            x_n = x_0 - f(x_0, 0)/f(x_0, 1)
            if abs(x_n - x_0) < 0.00001:
                return x_n
            x_0 = x_n
    else:
        print("Не выполнены условия сходимости для отрезка", "[", a, ";", b, "]")
        return -1

y = []
r1 = R_1
r2 = R_2
while budan_fourier(r1, r2) > 1:
    r2 = (r1 + r2)/2
r3 = r2
y.append(newton(r1, r2))
r1 = r2
r2 = R_2
while budan_fourier(r1, r2) > 1:
    r1 = (r1 + r2)/2
r4 = r1
y.append(newton(r1, r2))
while budan_fourier(r3, r4) > 1:
   r4 = (r3 + r4)/2
y.append(newton(r3, r4))
for i in range(0, len(y)):
    print("Корень", i+1, ":", y[i])
