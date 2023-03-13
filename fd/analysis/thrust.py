from math import exp, sqrt, pow


ittmax = 730
Ne = 2
NL = 104
tto = 0
pto = 0
po = 0
t_o = 0
T_isa = 0
rho = 0
nu = 0
mu = 0
a_s = 0
elpcc = 0
r = 0
tel = 0
mh = 0
mc = 0
wh = 0
wc = 0
tt5 = 0
mf = 0
tet = 0
tt1 = 0
tt2 = 0
tt3 = 0
tt4 = 0
tt6 = 0
tt7 = 0
ttb = 0
t8 = 0
t9 = 0
pa = 0
pt1 = 0
pt2 = 0
itt = 0
s = 0
pt3 = 0
pt4 = 0
pt5 = 0
pt6 = 0
pt7 = 0
bpr = 0
dnct = 0
nr = 0
nf = 0
nc = 0
nb = 0
dpt = 0
nt = 0
nm = 0
nnc = 0
nnh = 0
mht3p3 = 0
ehpc = 0
mhd = 0
tt3d = 0
b = 0
thetat = 0
c = 0
d = 0
Ah = 0
Ac = 0
pce = 0
phe = 0
p3krit = 0
p7krit = 0
dncn = 0
dhnr = 0
dtnr = 0
dhmnr = 0
vo = 0
tt4max = 0
itrel = 0
dmhnr = 0
theta = 0
NLcor = 0
NLcort1 = 0
NLcort = 0
dncn1 = 0
elpc = 0
elpc1 = 0
nf1 = 0
fi = 0
deltem = 0
delmh = 0
Tn = 0


def atmos(h, M, T_static):
    global T_isa, po, t_o, tto, pto, a_s
    T_isa = 288.15 - 0.0065 * h
    po = 101325 * pow(T_isa / 288.15, 5.256)
    if h >= 11000:
        T_isa = 216.65
        po = 22631.23 * exp(-9.80665 / 216.65 / 287.05 * (h - 11000))
    t_o = T_static
    # rho=po/287.05/temp
    tto = t_o * (1 + 0.2 * M**2)
    pto = po * (1 + 0.2 * M**2) ** 3.5
    a_s = sqrt(1.4 * 287.05 * t_o)
    # mu=0.000017894*pow((t_o/288.15),1.50)*(288.15+110.4)/(to+110.4)
    # nu=mu/rho
    return T_isa


def stuw(h, M, dtemp, mfi):
    global ittmax, Ne, NL, tto, pto, po, t_o, T_isa, rho, nu, mu, a_s, elpcc, r, tel, mh, mc, wh, wc, tt5, mf, tet, tt1, tt2, tt3, tt4, tt6, tt7, ttb, t8, t9, pa, pt1, pt2, itt, s, pt3, pt4, pt5, pt6, pt7, bpr, dnct, nr, nf, nc, nb, dpt, nt, nm, nnc, nnh, mht3p3, ehpc, mhd, tt3d, b, thetat, c, d, Ah, Ac, pce, phe, p3krit, p7krit, dncn, dhnr, dtnr, dhmnr, vo, tt4max, itrel, dmhnr, theta, NLcor, NLcort1, NLcort, dncn1, elpc, elpc1, nf1, fi, deltem, delmh, Tn

    elpcc = 0
    theta = t_o / 288.15
    thetat = tto / 288.15
    NLcor = NL / sqrt(theta)
    NLcort = NL / sqrt(thetat)
    #    Rendement van de verschillende componenten
    #    van de motor
    #   dnc = 0#
    dncn = -0.1209875 * (NLcort - 104) + 0.0005625 * (pow(NLcort, 2) - pow(104, 2))
    dnct = (4.5 * pow(10, -4) * dtemp + 1.5 * pow(10, -5) * pow(dtemp, 2)) * (
        1 - 0.69263 * h / 2438.4
    )
    dtnr = 0.000032 * pow(dtemp, 2) + 0.000233 * dtemp
    dhnr = -0.02 * (h / 2438.4) * (dtemp / 30)
    dhmnr = (
        -0.001 * pow((h / 3048), 4)
        + 0.0081667 * pow((h / 3048), 3)
        - 0.0085 * pow((h / 3048), 2)
        + 0.0023333 * (h / 3048)
        - 0.001
    )
    dmhnr = 0
    if M >= 0.2:
        dmhnr = (-0.12 * pow(M, 2) + 0.024 * M) * (
            -0.0625 * pow((h / 3048), 2) + 0.5 * h / 3048
        )
    nr = 0.989
    nf = 0.7
    nc = 0.73
    nb = 0.972
    dpt = 0.95
    nt = 0.86
    nm = 0.985
    nnc = 0.925
    nnh = 0.96
    Ac = 0.0779
    Ah = 0.05244
    tt4max = 635
    vo = M * a_s
    pt2 = nr * pto
    tt2 = tto
    tel = 0
    while not mfi == 0:
        #        start iteratieloop voor motorgegevens
        #        bij gegeven brandstofstroom mfi (Tt5=0)
        r = 0
        s = 0
        b = 0
        mhd = 10
        tt3d = 337.466
        c = 0
        while b < 50:
            c = c + 1
            d = 0
            while d < 30:
                d = d + 1
                ttb = nb * 41.865 * pow(10, 6) * mfi / 1147 / mhd
                fi = 0.3452334 * (1 + ttb / tt3d)
                ehpc = pow(
                    1 + fi * (pow(6.5625, (0.4 / 1.4 / nc)) - 1), (nc * 1.4 / 0.4)
                )
                mht3p3 = 0.0011217 * ehpc / sqrt(fi) / 6.5625
                pt3 = mhd * sqrt(tt3d) / mht3p3
                elpc = pt3 / pt2
                tt3 = tt2 * pow(elpc, (0.4 / 1.4 / nf))
                deltem = tt3d / tt3
                if deltem > 1.0005 or deltem < 0.9995:
                    tt3d = tt3d + (tt3 - tt3d) * 1.3
                    if tt3d < 100:
                        break  # inner loop

            if W_end0():
                break

        if W_end1():
            break

    while itt >= ittmax:
        mfi = 1450 * 0.4536 / 3600
        #       start iteratieloop voor motorgegevens bij overschrijding
        #       van ITTmax
        r = 0
        s = 0
        while r < 50:
            b = 0
            mhd = 10
            tt3d = 337.466
            c = 0
            while b < 50:
                c = c + 1
                d = 0
                while d < 30:
                    d = d + 1
                    ttb = nb * 41.865 * pow(10, 6) * mfi / 1147 / mhd
                    fi = 0.3452334 * (1 + ttb / tt3d)
                    ehpc = pow(
                        1 + fi * (pow(6.5625, (0.4 / 1.4 / nc)) - 1), (nc * 1.4 / 0.4)
                    )
                    mht3p3 = 0.0011217 * ehpc / sqrt(fi) / 6.5625
                    pt3 = mhd * sqrt(tt3d) / mht3p3
                    elpc = pt3 / pt2
                    tt3 = tt2 * pow(elpc, (0.4 / 1.4 / nf))
                    deltem = tt3d / tt3
                    if deltem > 1.0005 and deltem < 0.9995:
                        tt3d = tt3d + (tt3 - tt3d) * 1.3
                        if tt3d < 100:
                            break

                if W_end4():
                    break

        if tel < 20:
            if 12231.4 * elpc - 13041.14 < 0:
                NLcort = 30
            else:
                NLcort = 23.198 + sqrt(12231.4 * elpc - 13041.14)
            dncn = (
                46.449
                - 2.681348 * NLcort
                + 0.067482533 * pow(NLcort, 2)
                - 0.00072636 * pow(NLcort, 3)
                + 0.0000027318 * pow(NLcort, 4)
            ) / 100
            dncn = -0.1209875 * (NLcort - 104) + 0.0005625 * (
                pow(NLcort, 2) - pow(104, 2)
            )
            nf1 = (
                -14803
                - 0.2407085
                + 96754 * elpc
                + 0.939903 * elpc
                - 279446 * pow(elpc, 2)
                - 0.26169 * pow(elpc, 2)
                + 468125 * pow(elpc, 3)
                + 0.06834 * pow(elpc, 3)
                - 501269 * pow(elpc, 4)
                - 0.62651 * pow(elpc, 4)
                + 355822 * pow(elpc, 5)
                + 0.2019 * pow(elpc, 5)
                - 167440 * pow(elpc, 6)
                - 0.79086 * pow(elpc, 6)
                + 50370 * pow(elpc, 7)
                + 0.802172 * pow(elpc, 7)
                - 8790 * pow(elpc, 8)
                - 0.376607 * pow(elpc, 8)
                + 678 * pow(elpc, 9)
                + 0.0601726 * pow(elpc, 9)
            )
            nc = nf1 - 0.015
            if abs(elpc1 / elpc - 1) < 0.001:
                break
            nf = nf1
            elpc1 = elpc
            tel = tel + 1
        else:
            break

    return Tn


def W_end0():
    global ittmax, Ne, NL, tto, pto, po, t_o, T_isa, rho, nu, mu, a_s, elpcc, r, tel, mh, mc, wh, wc, tt5, mf, tet, tt1, tt2, tt3, tt4, tt6, tt7, ttb, t8, t9, pa, pt1, pt2, itt, s, pt3, pt4, pt5, pt6, pt7, bpr, dnct, nr, nf, nc, nb, dpt, nt, nm, nnc, nnh, mht3p3, ehpc, mhd, tt3d, b, thetat, c, d, Ah, Ac, pce, phe, p3krit, p7krit, dncn, dhnr, dtnr, dhmnr, vo, tt4max, itrel, dmhnr, theta, NLcor, NLcort1, NLcort, dncn1, elpc, elpc1, nf1, fi, deltem, delmh, Tn

    if pt3 / po < 1:
        if b > 15:
            return True
        b = b + 1
        mhd = mhd - 0.5
        c = 0
    else:
        p3krit = pow((1 - 0.4 / 2.4 / nnc), (-1.4 / 0.4))
        if pt3 / po < p3krit:
            wc = sqrt(2 * nnc * 1005 * tt3 * (1 - pow((po / pt3), (0.4 / 1.4))))
            pce = po
        else:
            pce = pt3 / p3krit
            t9 = tt3 * (1 - nnc * (1 - pow((pce / pt3), (0.4 / 1.4))))
            wc = sqrt(1.4 * 287.05 * t9)
        mc = Ac * pce * wc / 287.05 / t9
        tt4 = tt3 * pow(ehpc, (0.4 / 1.4 / nc))
        tt5 = tt4 + ttb
        tt6 = tt5 - 1005 / 1147 / nm * (tt4 - tt3)
        bpr = mc / mhd
        tt7 = tt6 - 1005 * (1 + bpr) * (tt3 - tt2) / 1147 / nm
        if tt7 < 10:
            if b > 20:
                return True
            mhd = mhd - 0.5
            b = b + 1
            c = 0
        else:
            pt5 = pto * nr * elpc * ehpc * dpt
            pt6 = pt5 * pow((tt6 / tt5), (1.33 / 0.33 / nt))
            pt7 = pt6 * pow((tt7 / tt6), (1.33 / 0.33 / nt))
            if pt7 / po < 1:
                if b > 35:
                    return True
                mhd = mhd - 0.25
                b = b + 1
                c = 0
            else:
                p7krit = pow((1 - 0.33 / 2.33 / nnh), (-1.33 / 0.33))
                if pt7 / po < p7krit:
                    wh = sqrt(
                        2 * nnh * 1147 * tt7 * (1 - pow((po / pt7), (0.33 / 1.33)))
                    )
                    t8 = tt7 * (1 - nnh * (1 - pow((po / pt7), (0.33 / 1.33))))
                    phe = po
                else:
                    phe = pt7 / p7krit
                    t8 = tt7 * (1 - nnh * (1 - pow((phe / pt7), (0.33 / 1.33))))
                    if t8 < 20:
                        return True
                    wh = sqrt(1.33 * 287.05 * t8)
                mh = Ah * wh * phe / 287.05 / t8
                delmh = mh / mhd
                if c >= 40:
                    return
                if delmh < 0.999:
                    mhd = mhd + 0.1 * (mh - mhd)
                elif delmh > 1.001:
                    mhd = mhd + 0.05 * (mh - mhd)
                else:
                    Tn = (
                        mc * (wc - vo)
                        + mh * (wh - vo)
                        + Ac * (pce - po)
                        + Ah * (phe - po)
                    )
                    # echo "Tn (1): ", Tn, "<br>"
                    mf = ttb / nb / 41.875 * pow(10, 6) * 1147 * mh
                    itt = tt7 + 3 * (tt3 - tt2) - 273.15
                    return True
    return False


def W_end1():
    global ittmax, Ne, NL, tto, pto, po, t_o, T_isa, rho, nu, mu, a_s, elpcc, r, tel, mh, mc, wh, wc, tt5, mf, tet, tt1, tt2, tt3, tt4, tt6, tt7, ttb, t8, t9, pa, pt1, pt2, itt, s, pt3, pt4, pt5, pt6, pt7, bpr, dnct, nr, nf, nc, nb, dpt, nt, nm, nnc, nnh, mht3p3, ehpc, mhd, tt3d, b, thetat, c, d, Ah, Ac, pce, phe, p3krit, p7krit, dncn, dhnr, dtnr, dhmnr, vo, tt4max, itrel, dmhnr, theta, NLcor, NLcort1, NLcort, dncn1, elpc, elpc1, nf1, fi, deltem, delmh, Tn

    if c >= 40:
        return True
    if tel < 20:
        if 12231.4 * elpc - 13041.14 < 0:
            NLcort = 30
        else:
            NLcort = 23.198 + sqrt(12231.4 * elpc - 13041.14)
        dncn = (
            46.449
            - 2.681348 * NLcort
            + 0.067482533 * pow(NLcort, 2)
            - 0.00072636 * pow(NLcort, 3)
            + 0.0000027318 * pow(NLcort, 4)
        ) / 100
        dncn = -0.1209875 * (NLcort - 104) + 0.0005625 * (pow(NLcort, 2) - pow(104, 2))
        nf1 = (
            -14803
            - 0.2407085
            + 96754 * elpc
            + 0.939903 * elpc
            - 279446 * pow(elpc, 2)
            - 0.26169 * pow(elpc, 2)
            + 468125 * pow(elpc, 3)
            + 0.06834 * pow(elpc, 3)
            - 501269 * pow(elpc, 4)
            - 0.62651 * pow(elpc, 4)
            + 355822 * pow(elpc, 5)
            + 0.2019 * pow(elpc, 5)
            - 167440 * pow(elpc, 6)
            - 0.79086 * pow(elpc, 6)
            + 50370 * pow(elpc, 7)
            + 0.802172 * pow(elpc, 7)
            - 8790 * pow(elpc, 8)
            - 0.376607 * pow(elpc, 8)
            + 678 * pow(elpc, 9)
            + 0.0601726 * pow(elpc, 9)
        )
        nc = nf1 - 0.015
        elpcc = elpc1 / elpc - 1
        if elpcc < 0:
            elpcc = -elpcc
        if elpcc < 0.001:
            return True
        nf = nf1
        elpc1 = elpc
        tel = tel + 1
    else:
        return True
    return False


def W_end4():
    global ittmax, Ne, NL, tto, pto, po, t_o, T_isa, rho, nu, mu, a_s, elpcc, r, tel, mh, mc, wh, wc, tt5, mf, tet, tt1, tt2, tt3, tt4, tt6, tt7, ttb, t8, t9, pa, pt1, pt2, itt, s, pt3, pt4, pt5, pt6, pt7, bpr, dnct, nr, nf, nc, nb, dpt, nt, nm, nnc, nnh, mht3p3, ehpc, mhd, tt3d, b, thetat, c, d, Ah, Ac, pce, phe, p3krit, p7krit, dncn, dhnr, dtnr, dhmnr, vo, tt4max, itrel, dmhnr, theta, NLcor, NLcort1, NLcort, dncn1, elpc, elpc1, nf1, fi, deltem, delmh, Tn

    if pt3 / po < 1:
        if b > 15:
            return True
        b = b + 1
        mhd = mhd - 0.5
        c = 0
    else:
        p3krit = pow((1 - 0.4 / 2.4 / nnc), (-1.4 / 0.4))
        if pt3 / po < p3krit:
            wc = sqrt(2 * nnc * 1005 * tt3 * (1 - pow((po / pt3), (0.4 / 1.4))))
            t9 = tt3 * (1 - nnc * (1 - pow((po / pt3), (0.4 / 1.4))))
            pce = po
        else:
            pce = pt3 / p3krit
            t9 = tt3 * (1 - nnc * (1 - pow((pce / pt3), (0.4 / 1.4))))
            wc = sqrt(1.4 * 287.05 * t9)
        mc = Ac * pce * wc / 287.05 / t9
        tt4 = tt3 * pow(ehpc, (0.4 / 1.4 / nc))
        tt5 = tt4 + ttb
        tt6 = tt5 - 1005 / 1147 / nm * (tt4 - tt3)
        bpr = mc / mhd
        tt7 = tt6 - 1005 * (1 + bpr) * (tt3 - tt2) / 1147 / nm
        if tt7 < 10:
            if b > 20:
                return True
            mhd = mhd - 0.5
            b = b + 1
            c = 0
        else:
            pt5 = pto * nr * elpc * ehpc * dpt
            pt6 = pt5 * pow((tt6 / tt5), (1.33 / 0.33 / nt))
            pt7 = pt6 * pow((tt7 / tt6), (1.33 / 0.33 / nt))
            if pt7 / po < 1:
                if b > 35:
                    return True
                mhd = mhd - 0.25
                b = b + 1
                c = 0
            else:
                p7krit = pow((1 - 0.33 / 2.33 / nnh), (-1.33 / 0.33))
                if pt7 / po < p7krit:
                    wh = sqrt(
                        2 * nnh * 1147 * tt7 * (1 - pow((po / pt7), (0.33 / 1.33)))
                    )
                    t8 = tt7 * (1 - nnh * (1 - pow((po / pt7), (0.33 / 1.33))))
                    phe = po
                else:
                    phe = pt7 / p7krit
                    t8 = tt7 * (1 - nnh * (1 - pow((phe / pt7), (0.33 / 1.33))))
                    if t8 < 20:
                        return True
                    wh = sqrt(1.33 * 287.05 * t8)
                mh = Ah * wh * phe / 287.05 / t8
                delmh = mh / mhd
                if c < 100:
                    if delmh < 0.999:
                        mhd = mhd + 0.1 * (mh - mhd)
                    elif delmh > 1.001:
                        mhd = mhd + 0.05 * (mh - mhd)
                else:
                    Tn = (
                        mc * (wc - vo)
                        + mh * (wh - vo)
                        + Ac * (pce - po)
                        + Ah * (phe - po)
                    )
                    # echo "Tn (2): ", Tn, "<br>"
                    mf = ttb / nb / 41.875 * pow(10, 6) * 1147 * mh
                    itt = tt7 + 3 * (tt3 - tt2) - 273.15
                    print(itt, ittmax)
                    itrel = itt / ittmax
                    if itrel - 1 < 0:
                        itrel = -(itrel - 1)
                    if itrel > 0.0001:
                        if itrel > 1.5:
                            itrel = 1.5
                        if itrel < 0.5:
                            itrel = 0.5
                        mfi = mfi - (itrel - 1) * mfi
                        r = r + 1
                        if r > 50:
                            return True
                        else:
                            b = 0
                            mhd = 10
                            tt3d = 337.466
                            c = 0
                    else:
                        return True
    return False


def calculate_thrust(h, M, T, fuelflow):
    """
    Calculates thrust based on operating conditions.

    Args:
        h: Altitude [m]
        M: Mach number [-]
        T: Static temperature [K]
        fuelflow: Fuel mass flow [kg/s]

    Returns:
        Thrust [N]
    """
    T_isa = atmos(h, M, T)
    delta_T = T - T_isa
    return stuw(h, M, delta_T, fuelflow)
