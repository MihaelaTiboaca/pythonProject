import math as m
import numpy as np


class Point:
    def __init__(self, id, x, y, X, Y, Z):
        self.x = x
        self.y = y
        self.X = X
        self.Y = Y
        self.Z = Z
        self.id = id

    def __str__(self):
        return "id:{0}, x:{1}, y:{2}, X:{3}, Y:{4}, Z:{5} ".format(self.id, self.x, self.y, self.X, self.Y, self.Z)

    def distance(self, p, ground: bool):
        if ground == True:
            dx = p.X - self.X
            dy = p.Y - self.Y
            return m.sqrt(pow(dx, 2) + pow(dy, 2))
        return m.sqrt(pow(p.x - self.x, 2) + pow(p.y - self.y, 2))


def aprox_height(p1, p2, f):
    AB = p1.distance(p2, True)
    deltax = p2.x - p1.x
    deltay = p2.y - p1.y
    H = (AB * f) / m.sqrt(deltax ** 2 + deltay ** 2)
    return H


def Aprox_Values(points):
    A = []
    T = []
    # product= []
    for i in range(0, len(points)):
        p = points[i]
        A.append([p.x, -p.y, 1, 0])
        A.append([p.x, p.y, 0, 1])
        T.append(p.X)
        T.append(p.Y)
    AT = np.transpose(A)
    N = np.matmul(AT, A)
    invN = np.linalg.inv(N)
    B = np.matmul(AT, T)
    X = np.matmul(invN, B)
    a = X[0]
    b = X[1]
    X0 = X[2]
    Y0 = X[3]
    k = m.atan(a / b)
    return (X0, Y0, k)


def Aprox_Values2(points):
    A = []
    T = []
    # product= []
    for i in range(0, len(points)):
        p = points[i]
        A.append([p.x, -p.y, 1, 0])
        A.append([p.x, p.y, 0, 1])
        T.append(p.X)
        T.append(p.Y)
    X = np.linalg.lstsq(A, T)
    a = X[0][0]
    b = X[0][1]
    X0 = X[0][2]
    Y0 = X[0][3]
    k = m.atan2(a, b)
    return (X0, Y0, k)


def Read_Points(fn):
    # with statement-  de cautat
    file = open(fn)
    lines = file.readlines()
    count = 0
    f = 0
    points = []
    for line in lines:
        info = line.split("|")
        if count == 0:
            f = float(info[0])
        else:
            id = info[0]
            x = float(info[1])
            y = float(info[2])
            X = float(info[3])
            Y = float(info[4])
            Z = float(info[5])
            p = Point(id, x, y, X, Y, Z)
            points.append(p)
        count += 1
    file.close()
    return (f, points)


def f_m11(phi, kappa, omega):
    return m.cos(phi) * m.cos(kappa)


def f_m12(phi, kappa, omega):
    return m.sin(omega) * m.sin(phi) * m.cos(kappa) + m.cos(omega) * m.sin(kappa)


def f_m13(phi, kappa, omega):
    return -m.cos(omega) * m.sin(phi) * m.cos(kappa) + m.sin(omega) * m.sin(kappa)


def f_m21(phi, kappa, omega):
    return -m.cos(phi) * m.sin(kappa)


def f_m22(phi, kappa, omega):
    return -m.sin(omega) * m.sin(phi) * m.sin(kappa) + m.cos(omega) * m.cos(kappa)


def f_m23(phi, kappa, omega):
    return m.cos(omega) * m.sin(phi) * m.sin(kappa) + m.sin(omega) * m.cos(kappa)


def f_m31(phi, kappa, omega):
    return m.sin(phi)


def f_m32(phi, kappa, omega):
    return m.sin(omega) * m.cos(phi)


def f_m33(phi, kappa, omega):
    return m.cos(omega) * m.cos(phi)


def f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL):
    m11 = f_m11(phi, kappa, omega)
    m13 = f_m13(phi, kappa, omega)
    m12 = f_m12(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return m11 * deltaX + m12 * deltaY + m13 * deltaZ


def f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL):
    m31 = f_m31(phi, kappa, omega)
    m32 = f_m32(phi, kappa, omega)
    m33 = f_m33(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return m31 * deltaX + m32 * deltaY + m33 * deltaZ


def f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL):
    m21 = f_m21(phi, kappa, omega)
    m22 = f_m22(phi, kappa, omega)
    m23 = f_m23(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return m21 * deltaX + m22 * deltaY + m23 * deltaZ


def f_b11(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m33 = f_m33(phi, kappa, omega)
    m32 = f_m32(phi, kappa, omega)
    m13 = f_m13(phi, kappa, omega)
    m12 = f_m12(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (r * (-m33 * deltaY + m32 * deltaZ) - q * (-m13 * deltaY + m12 * deltaZ))


def f_b12(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (r * (
                m.cos(phi) * deltaX + m.sin(omega) * m.sin(phi) * deltaY - m.cos(omega) * m.sin(phi) * deltaZ) - q * (
                                         -m.sin(phi) * m.cos(kappa) * deltaX + m.sin(omega) * m.cos(phi) * m.cos(
                                     kappa) * deltaY - m.cos(omega) * m.cos(phi) * m.cos(kappa) * deltaZ))


def f_b13(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m21 = f_m21(phi, kappa, omega)
    m22 = f_m22(phi, kappa, omega)
    m23 = f_m23(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (-f / q) * (m21 * deltaX + m22 * deltaY + m21 * deltaZ)


def f_b14(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m31 = f_m31(phi, kappa, omega)
    m11 = f_m11(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (r * m31 - q * m11)


def f_b15(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m32 = f_m32(phi, kappa, omega)
    m12 = f_m12(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (r * m32 - q * m12)


def f_b16(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m33 = f_m33(phi, kappa, omega)
    m13 = f_m13(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (r * m33 - q * m13)


def f_b21(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    s = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m33 = f_m33(phi, kappa, omega)
    m32 = f_m32(phi, kappa, omega)
    m23 = f_m23(phi, kappa, omega)
    m22 = f_m22(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (s * (-m33 * deltaY + m32 * deltaZ) - q * (-m23 * deltaY + m22 * deltaZ))


def f_b22(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    s = f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (s * (
                m.cos(phi) * deltaX + m.sin(omega) * m.sin(phi) * deltaY - m.cos(omega) * m.sin(phi) * deltaZ) - q * (
                                         m.sin(phi) * m.cos(kappa) * deltaX - m.sin(omega) * m.cos(phi) * m.sin(
                                     kappa) * deltaY + m.cos(omega) * m.cos(phi) * m.sin(kappa) * deltaZ))


def f_b23(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    s = f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m21 = f_m21(phi, kappa, omega)
    m12 = f_m12(phi, kappa, omega)
    m13 = f_m13(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (-f / q) * (m21 * deltaX + m12 * deltaY + m13 * deltaZ)


def f_b24(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    s = f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m31 = f_m31(phi, kappa, omega)
    m21 = f_m21(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (s * m31 - q * m21)


def f_b25(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    s = f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m32 = f_m32(phi, kappa, omega)
    m22 = f_m22(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (s * m32 - q * m22)


def f_b26(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f):
    s = f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    q = f_q(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
    m33 = f_m33(phi, kappa, omega)
    m23 = f_m23(phi, kappa, omega)
    deltaX = XA - XL
    deltaY = YA - YL
    deltaZ = ZA - ZL
    return (f / (q ** 2)) * (s * m33 - q * m23)


def is_small(a, b):
    dif = a - b
    if dif < 0:
        return True
    return False


def Least_Squared(f, points):
    Z0 = aprox_height(points[0], points[1], f)
    (X0, Y0, k) = Aprox_Values(points)
    XL = X0
    YL = Y0
    ZL = Z0
    omega = 0
    phi = 0
    kappa = k
    D = np.zeros((6, 1))
    nr_points = len(points)
    mm = 2 * len(points)
    n = 6
    B = np.zeros((mm, n))
    E = np.zeros((mm, 1))
    N = np.zeros((n, n))
    Q = []
    V = []
    tol = 0.00001
    cond = True
    count = 0
    while cond:
        # calculam coeficientii matricii B
        for i in range(0, len(points)):
            p = points[i]
            XA = p.X
            YA = p.Y
            ZA = p.Z
            XL = XL + D[3]
            YL = YL + D[4]
            ZL = ZL + D[5]
            omega = omega + D[0]
            phi = phi + D[1]
            kappa = kappa + D[2]
            xa = p.x
            ya = p.y
            r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
            q = f_q(phi, kappa, XA, YA, ZA, XL, YL, ZL)
            s = f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
            # coef x
            b11 = f_b11(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b12 = f_b12(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b13 = f_b13(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b14 = f_b14(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b15 = f_b15(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b16 = f_b16(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            # coef y
            b21 = f_b21(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b22 = f_b22(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b23 = f_b23(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b24 = f_b24(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b25 = f_b25(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            b26 = f_b26(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
            # B.append([b11, b12, b13, b14, b15, b16])
            # B.append([b21, b22, b23, b24, b25, b26])
            B[2 * i] = [b11, b12, b13, b14, b15, b16]
            B[2 * i + 1] = [b21, b22, b23, b24, b25, b26]
            J = xa + f * r / q
            K = ya + f * s / q
            E[2 * i] = J
            E[2 * i + 1] = K
        # cele mai mici patrate
        if count > 0:
            V = np.matmul(B, D) - E
        BT = np.transpose(B)
        N = np.matmul(BT, B)
        Q = np.linalg.inv(N)
        L = np.matmul(BT, E)
        D = np.matmul(Q, L)
        cond = (not (is_small(D[0], tol) and is_small(D[1], tol) and is_small(D[2], tol)))
        count += 1
    print('omega:{0},phi:{1}, kappa:{2}, XL:{3},YL:{4},ZL:{5}'.format((omega * 180 / m.pi), (phi * 180 / m.pi),
                                                                      (kappa * 180 / m.pi), XL, YL, ZL))
    VT = np.transpose(V)
    P = np.matmul(VT, V)
    vv = P
    r = float(mm - n)
    sigma = m.sqrt(vv / r)
    Sigmaii = np.zeros((n, 1))
    for i in range(0, n):
        Sigmaii[i] = sigma * Q[i][i]
    print(sigma)
    print(Sigmaii)


def Least_Squared2(f, points):
    Z0 = aprox_height_avg(points, f)
    (X0, Y0, k) = Aprox_Values(points)
    XL = X0
    YL = Y0
    ZL = Z0
    omega = 0
    phi = 0
    kappa = k
    mm = 2 * len(points)
    n = 6
    B = np.zeros((mm, n))
    E = np.zeros((mm, 1))
    # calculam coeficientii matricii B
    for i in range(0, len(points)):
        p = points[i]
        XA = p.X
        YA = p.Y
        ZA = p.Z
        xa = p.x
        ya = p.y
        r = f_r(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
        q = f_q(phi, kappa, XA, YA, ZA, XL, YL, ZL)
        s = f_s(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL)
        # coef x
        b11 = f_b11(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b12 = f_b12(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b13 = f_b13(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b14 = f_b14(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b15 = f_b15(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b16 = f_b16(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        # coef y
        b21 = f_b21(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b22 = f_b22(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b23 = f_b23(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b24 = f_b24(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b25 = f_b25(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        b26 = f_b26(phi, kappa, omega, XA, YA, ZA, XL, YL, ZL, f)
        # B.append([b11, b12, b13, b14, b15, b16])
        # B.append([b21, b22, b23, b24, b25, b26])
        B[2 * i] = [b11, b12, b13, b14, b15, b16]
        B[2 * i + 1] = [b21, b22, b23, b24, b25, b26]
        J = xa + f * r / q
        K = ya + f * s / q
        E[2 * i] = J
        E[2 * i + 1] = K
    Result = np.linalg.lstsq(B, E)
    print(Result[0])


def aprox_height_avg(points, f):
    Havg = 0
    count = 0
    for i in range(0, len(points) - 1):
        p1 = points[i]
        p2 = points[i + 1]
        Havg = Havg + aprox_height(p1, p2, f)
        count += 1
    Havg = Havg / count
    return Havg


(f, points) = Read_Points(r"C:\Users\Admin\Desktop\Facultate\ANUL 2\sem 2\Fotogrammetrie\Date.txt")
Least_Squared(f, points)
# Z0 = aprox_height_avg(points, f)
# (X0, Y0, k) = Aprox_Values2(points)
# print(X0, Y0, Z0, k)
# Least_Squared2(f,points)
