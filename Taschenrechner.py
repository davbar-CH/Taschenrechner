import inspect

import numpy as np
from tkinter import *
import tkinter.ttk as ttk
import tkinter as tk
from itertools import permutations

import scipy
import sympy
from sympy import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.patches import Polygon
import re


def zweimalzwei(a1, b1, a2, b2, c1, c2, anzeige):
    anzeige.delete(1.0, END)
    y = ((a2 * c1) - (a1 * c2)) / ((a2 * b1) - (a1 * b2))
    x = (c1 - b1 * y) / a1
    anzeige.insert(END, f"a ist {x} und b ist {y}")


def dreimaldrei(a1, b1, c1, a2, b2, c2, a3, b3, c3, d1, d2, d3, anzeige):
    anzeige.delete(1.0, END)
    rest = (a1 * (b2 * c3 - b3 * c2) - b1 * (a2 * c3 - a3 * c2) + c1 * (a2 * b3 - a3 * b2))
    z = (a1 * (b2 * d3 - b3 * d2) - b1 * (a2 * d3 - a3 * d2) + d1 * (a2 * b3 - a3 * b2)) / rest
    y = (a1 * (d2 * c3 - d3 * c2) - d1 * (a2 * c3 - a3 * c2) + c1 * (a2 * d3 - a3 * d2)) / rest
    x = (d1 * (b2 * c3 - b3 * c2) - b1 * (d2 * c3 - d3 * c2) + c1 * (d2 * b3 - d3 * b2)) / rest
    anzeige.insert(END, f"a ist {x} ,b ist {y} und c ist {z}")


def viermalvier(a1, b1, c1, d1, r1, a2, b2, c2, d2, r2, a3, b3, c3, d3, r3, a4, b4, c4, d4, r4, anzeige):
    anzeige.delete(1.0, END)
    A = np.array([
        [a1, b1, c1, d1],
        [a2, b2, c2, d2],
        [a3, b3, c3, d3],
        [a4, b4, c4, d4]
    ], dtype=float)

    f = np.array([r1, r2, r3, r4], dtype=float)
    x = np.linalg.solve(A, f)

    anzeige.insert(END, f"a ist {x[0]} ,b ist {x[1]}, c ist {x[2]} und d ist {x[3]}")


def kubisch(a, b, c, e):
    p = ((3 * a * c) - (b ** 2)) / (3 * (a ** 2))
    q = ((2 * (b ** 3)) - (9 * a * b * c) + (27 * (a ** 2) * e)) / (27 * (a ** 3))

    diskriminante = ((q / 2) ** 2) + ((p / 3) ** 3)
    if diskriminante > 1 * (10 ** -9):

        u = ((-(q / 2)) + np.sqrt(((q / 2) ** 2) + ((p / 3) ** 3))) ** (1 / 3)
        v = ((-(q / 2)) - np.sqrt(((q / 2) ** 2) + ((p / 3) ** 3))) ** (1 / 3)

        # Die drei Lösungen
        y1 = ((-(q / 2) + (((q / 2) ** 2) + ((p / 3) ** 3)) ** (1 / 2)) ** (1 / 3)) + (
                (-(q / 2) - (((q / 2) ** 2) + ((p / 3) ** 3)) ** (1 / 2)) ** (1 / 3))
        y2 = ((u + v) / -2) + (1j * (np.sqrt(3) / 2) * (u - v))
        y3 = ((u + v) / -2) - (1j * (np.sqrt(3) / 2) * (u - v))

        x1 = y1 + (b / (3 * a))
        x2 = y2 + (b / (3 * a))
        x3 = y3 + (b / (3 * a))
        print(f"Resultat1:{x1}, Resultat2:{x2},Resultat3:{x3}")
    elif abs(diskriminante) < 1 * (10 ** -9):
        if (1 * (10 ** -10) < p < 1 * (10 ** -9)) and (1 * (10 ** -10) < q < 1 * (10 ** -9)):
            print("Resultat ist 0")
        else:
            y1 = 2 * ((-(q / 2)) ** (1 / 3))
            y2 = -(-(q / 2)) ** (1 / 3)
            x1 = y1 + (b / (3 * a))
            x2 = y2 + (b / (3 * a))
            print("Resultat 1 ist", x1, "und Resultat 2 ist", x2)
    else:
        for i in range(0, 3):
            y = 2 * np.sqrt(-(p / 3)) * np.cos(
                (1 / 3) * np.arccos(((3 * q) / (2 * p)) * np.sqrt(-(3 / p))) - ((2 * np.pi * i) / 3))
            x1 = y + (b / (3 * a))
            print(f"Resultat", i + 1, "ist", x1)
    print(p)
    print(q)
    print(diskriminante)


def quadratisch(b, c, e):
    x1 = (-c + np.sqrt((c ** 2) - (4 * b * e))) / (2 * b)
    x2 = (-c - np.sqrt((c ** 2) - (4 * b * e))) / (2 * b)
    print(x1, x2)


def linear(k, c, d):
    if d > 1 * (10 ** -9):
        x = (k - d) / c
        print(x)
    elif d < 1 * (10 ** -9):
        x = (k + d) / c
        print(x)


def ableitung(funktion_expr, n_ableitung_int, anzeige):
    x = Symbol('x')
    abg_funktion = funktion_expr
    for i in range(n_ableitung_int):
        abg_funktion = abg_funktion.diff(x)

    if anzeige == "a":
        pass
    else:
        anzeige.insert(END, f"{n_ableitung_int}.Abgeleitete Funktion:{abg_funktion}" + "\n")

    return abg_funktion


def kurvendiskussion(funktion, anzeige):
    def extremstelle(funktion):
        # extremstelle
        x = Symbol('x')

        abg_funktion_1 = ableitung(funktion, 1, "a")
        x_res_liste = solve(abg_funktion_1, x)

        return x_res_liste

    def hoch_tief(funktion, anzeige):
        x = Symbol("x")
        abg_funktion_2 = ableitung(funktion, 2, "a")
        x_res_liste = extremstelle(funktion)

        for anzahl, x_res in enumerate(x_res_liste):
            y_res = funktion.subs(x, x_res)
            y_res_2 = abg_funktion_2.subs(x, x_res)

            if y_res_2 < 1 * (10 ** -9):
                anzeige.insert(END, f"{anzahl + 1}.Hochpunkt bei ({x_res}und{y_res})" + "\n")
            elif y_res_2 > 1 * (10 ** -9):
                anzeige.insert(END, f"{anzahl + 1}.Tiefpunkt bei ({x_res}und{y_res})" + "\n")

    def wendepunkt(funktion, anzeige):
        x = Symbol("x")
        abg_funktion_2 = ableitung(funktion, 2, "a")

        x_res_liste = solve(abg_funktion_2, x)
        for anzahl, x_res in enumerate(x_res_liste):
            y_res = funktion.subs(x, x_res)
            anzeige.insert(END, f"{anzahl + 1}. Wendepunkt bei ({x_res}und{y_res})" + "\n")

    def symmetrien(funktion, anzeige):
        x = Symbol("x")

        if funktion == funktion.subs(x, -x):
            anzeige.insert(END, "Symmetrisch zur y-Achse" + "\n")
        elif funktion == -1 * (funktion.subs(x, -x)):
            anzeige.insert(END, "Symmetrisch zum Ursprung" + "\n")

    hoch_tief(funktion, anzeige)
    wendepunkt(funktion, anzeige)
    symmetrien(funktion, anzeige)


def aufleitung(funktion_expr, n_aufleitung_int, anzeige):
    anzeige.delete(1.0, END)
    x = Symbol("x")
    aufg_funktion = funktion_expr
    for i in range(n_aufleitung_int):
        aufg_funktion = aufg_funktion.integrate(x)

    if anzeige == "a":
        pass
    else:
        anzeige.insert(END, f"{n_aufleitung_int}.Aufgeleitete Funktion:{aufg_funktion}" + "\n")
    return aufg_funktion


def integrale(funktion, anzeige, a, b):
    anzeige.delete(1.0, END)
    x = Symbol("x")
    area = sympy.integrate(funktion, (x, a, b))
    anzeige.insert(END, f"Die Fläche unter der Funktion{funktion} ist:{area}" + "\n")


# speicher, falls man einen Vektor behalten möchte
speicher_dic = {}

def vektorgeometrie(values, float_values, anzeige, befehl):
    anzeige.delete(1.0, END)

    def vektoren_einlesen(values, float_values):
        a_vektor = b_vektor = c_vektor = d_vektor = None

        if all(v[2] == '' for v in values[:4]):  # nachfolgende Befehle können natürlich vereinfacht werden
            if len(float_values[0]) > 0:  # aber so ist es klar und Einzeiler sind nichts für mich
                a_vektor = np.array(float_values[0])
            if len(float_values[1]) > 0:
                b_vektor = np.array(float_values[1])
            if len(float_values[2]) > 0:
                c_vektor = np.array(float_values[2])
            if len(float_values[3]) > 0:
                d_vektor = np.array(float_values[3])
            return a_vektor, b_vektor, c_vektor, d_vektor

        else:
            if len(float_values[0]) > 0:
                a_vektor = np.array(float_values[0])
            if len(float_values[1]) > 0:
                b_vektor = np.array(float_values[1])
            if len(float_values[2]) > 0:
                c_vektor = np.array(float_values[2])
            if len(float_values[3]) > 0:
                d_vektor = np.array(float_values[3])
            return a_vektor, b_vektor, c_vektor, d_vektor

    a_vektor, b_vektor, c_vektor, d_vektor = vektoren_einlesen(values, float_values)

    orts_vektoren_alle = [
        a_vektor,
        b_vektor,
        c_vektor,
        d_vektor
    ]

    orts_vektoren_richtig = [0 if vektor is None else vektor for vektor in orts_vektoren_alle]

    vektor_namen = [
        "A",
        "B",
        "C",
        "D"
    ]

    differenzen = [(name1, name2, v2 - v1)
                   for (name1, v1), (name2, v2) in permutations(zip(vektor_namen, orts_vektoren_richtig), 2)]

    vektoren_kombi = {f"{name1}-{name2}": diff for name1, name2, diff in differenzen}

    def vektoren_auslesen():
        vektoren_werte = []
        vektoren_namen = []
        resultat = re.search(r"(\b[A-Z]+\b-\b[A-Z]+\b)(.*)(\b[A-Z]+\b-\b[A-Z]+\b)", string=befehl)
        if resultat:
            vektoren = resultat.groups()
            print(vektoren)
            for vektor in vektoren:
                for key in vektoren_kombi:
                    if key == vektor:
                        vektoren_werte.append(vektoren_kombi[key])
                        vektoren_namen.append(key)
            return vektoren_werte, vektoren_namen, "mehrfach"

        einzelner_vektor = re.search(r"(\b[A-D]+\b-\b[A-D]+\b)", string=befehl)
        if einzelner_vektor:
            vektoren = einzelner_vektor.groups()
            print(vektoren)
            for vektor in vektoren:
                for key in vektoren_kombi:
                    if key == vektor:
                        vektoren_werte.append(vektoren_kombi[key])
                        vektoren_namen.append(key)
            return vektoren_werte, vektoren_namen, "einzel"

        anzeige.insert(END, f"Es gibt keinen solchen Vektor: {befehl}\n")
        return None

    def operation_auslesen():
        resultat = re.search(r"(\b[a-z]+\b)", string=befehl)
        if resultat:
            operation = resultat.groups(1)
            operation = operation[0]
            return operation
        else:
            return None

    operation = operation_auslesen()

    def skalarprodukt(vektor1, vektor2):
        skalarprodukt = 0
        for i in range(len(vektor1)):
            skalarprodukt = skalarprodukt + vektor1[i] * vektor2[i]

        return skalarprodukt

    def vektor_produkt(vektor1, vektor2):
        kreuzprodukt = np.cross(vektor1, vektor2)
        return kreuzprodukt

    def winkel(vektor1, vektor2):
        betrag_vektor1 = np.linalg.norm(vektor1)
        betrag_vektor2 = np.linalg.norm(vektor2)

        skalarprodukt = 0
        for i in range(len(vektor1)):
            skalarprodukt = skalarprodukt + vektor1[i] * vektor2[i]

        cos_winkel = skalarprodukt / (betrag_vektor1 * betrag_vektor2)
        cos_winkel = np.clip(cos_winkel, -1.0, 1.0)
        winkel_rad = np.arccos(cos_winkel)
        winkel_deg = np.rad2deg(winkel_rad)

        return winkel_deg

    def ebene(vektor1, vektor2, punkt_wert):
        kreuzprodukt = np.cross(vektor1, vektor2)
        print(type(kreuzprodukt))
        if type(kreuzprodukt) is np.ndarray:
            d = -kreuzprodukt[0] * punkt_wert[0] + kreuzprodukt[1] * punkt_wert[1] + kreuzprodukt[2] * punkt_wert[2]
            ebene_zwei = f"x:{kreuzprodukt[0]} y:{kreuzprodukt[1]} z:{kreuzprodukt[2]} + {d}"
            return ebene_zwei
        else:
            return None

    werte, namen, art = vektoren_auslesen()

    if art == "einzel":
        if operation == "speicher":
            speicher = re.search(r"(\b[S-Z]+\b)", string=befehl)
            speicher = speicher.group()
            speicher_dic.update({speicher: (werte, namen)})
            anzeige.insert(END, f"Vektor {namen} ist jetzt gespeichert als: {speicher}\n")

        if operation is None:
            if not bool(speicher_dic):
                anzeige.insert(END, f"{namen[0]} ist {werte[0]}\n")
            else:
                speicher = re.search(r"(\b[S-Z]+\b)", string=befehl)
                (werte, namen) = speicher_dic[speicher]
                anzeige.insert(END, f"{namen[0]} ist {werte[0]}\n")

    if art == "mehrfach":
        if operation == "sp":
            skalarp = skalarprodukt(werte[0], werte[1])
            anzeige.insert(END, f"Skalarprodukt von {namen[0]} und {namen[1]}: {skalarp}\n")

        if operation == "kp":
            kreuzp = vektor_produkt(werte[0], werte[1])
            anzeige.insert(END, f"Kreuzprodukt von {namen[0]} und {namen[1]}: {kreuzp}\n")

        if operation == "winkel":
            winkel = winkel(werte[0], werte[1])
            anzeige.insert(END, f"Winkel zwischen {namen[0]} und {namen[1]}: {winkel}\n")

        if operation == "ebene":
            punkt = re.search(r"(\b[A-Z]+\b)", string=befehl)
            punkt = punkt.group()

            punkt_wert = None
            if punkt == "A":
                punkt_wert = orts_vektoren_alle[0]
            if punkt == "B":
                punkt_wert = orts_vektoren_alle[1]
            if punkt == "C":
                punkt_wert = orts_vektoren_alle[2]
            if punkt == "D":
                punkt_wert = orts_vektoren_alle[3]

            ebene = ebene(werte[0], werte[1], punkt_wert)
            if ebene is not None:
                anzeige.insert(END, f"Ebene aus {namen[0]} und {namen[1]}: {ebene}\n")
            else:
                anzeige.insert(END, "Entweder ein 2D-Vektor (dasselbe wie Skalarprodukt) oder ein sonstiger Fehler\n")

        if operation == "und":
            for wert, name in zip(werte, namen):
                anzeige.insert(END, f"{name} ist {wert}\n")

def spulen_beschleunigung(werte, anzeige_fehler, felder_dic):
    r_kugel = (werte["r_kugel"], "r_kugel")
    u_permea_vakuum = 1.256637061 * (10**-6)
    u_relativ = (werte["u_relativ"], "u_relativ")

    if u_relativ is None or u_relativ == 0:
        u_relativ = (5150,"u_relativ")

    N_windungen = (werte["N_windungen"], "N_windungen")
    I = (werte["I"], "I")
    r_spule = (werte["r_spule"], "r_spule")
    l = (werte["Länge"], "Länge")

    x = symbols("x")

    f_l = ((-1) * np.pi * (r_kugel[0]**3) * u_permea_vakuum * (u_relativ[0] - 1) *
           (N_windungen[0]**2) * (I[0]**2) * (r_spule[0] ** 4) * (x / ((x**2) + (r_spule[0]**2))**4))

    f_l_abl = ableitung(f_l, 2, "a")
    schritt = 0.1
    beschleunigung_liste = []

    for x_werte in np.arange((-l[0] / 2), 0, schritt):
        beschleunigung = f_l_abl.subs(x, x_werte)
        beschleunigung_liste.append(beschleunigung)
    print(beschleunigung_liste)
    return f_l

def flugbahn_ohne_widerstand(werte, anzeige_fehler, felder_dic):

    # labels = ["Anfangshöhe", "Anfangsgeschwindigkeit", "Abwurfwinkel", "Flugzeit", "Distanz", "maximale Höhe"]
    g = 9.80665

    h = (werte["h"], "h")
    v0 = (werte["v0"], "v0")
    a = (werte["a"], "a")
    t = (werte["t"], "t")
    R = (werte["R"], "R")
    h_max = (werte["h_max"], "h_max")

    pflicht_param_list = [h[0], v0[0], a[0]]
    optional_param_list = [t[0], R[0], h_max[0]]
    alle_param_list = [h, v0, a, t, R, h_max]

    # flugzeit od. reichweite aus dem Stand ist einfach die lange Form mit h = 0
    def flugzeit(v0, a, g, h, t):
        a = np.deg2rad(a)
        x = v0 * np.sin(a)
        return ((x + np.sqrt(x ** 2 + (2 * g * h))) / g) - t

    def reichweite(v0, a, g, h, R):
        t = flugzeit(v0, a, g, h, 0)
        a = np.deg2rad(a)
        return (v0 * np.cos(a) * t) - R

    def maximale_hoehe(v0, a, g, h, h_max):
        a = np.deg2rad(a)
        return (h + (v0 ** 2 * np.sin(a) ** 2) / (2 * g)) - h_max

    def solve(f, args):
        i = args.index(None)
        sig = inspect.signature(f)
        param_name = list(sig.parameters.keys())[i]
        return scipy.optimize.fsolve(lambda x: f(*args[:i], x, *args[i + 1:]), 1), param_name

    optionale_werte = {
        "t": t,
        "R": R,
        "h_max": h_max,
    }

    optionale_formeln = {
        "t": flugzeit,
        "R": reichweite,
        "h_max": maximale_hoehe,
    }

    fehlende_pflicht = pflicht_param_list.count(None)
    fehlende_optional = optional_param_list.count(None)
    try: # was, wenn ich keine parameter habe?
        if (fehlende_pflicht >= 1 and fehlende_optional == 3) or (fehlende_pflicht == 2):
            anzeige_fehler.insert(END, "Es fehlt ein Parameter: v0, a oder h")

        if fehlende_pflicht == 0:
            ergebnisse = [
                solve(flugzeit, (v0[0], a[0], g, h[0], t[0])),
                solve(reichweite, (v0[0], a[0], g, h[0], R[0])),
                solve(maximale_hoehe, (v0[0], a[0], g, h[0], h_max[0]))
            ]

            for res, name in ergebnisse:
                if name in werte:
                    werte[name] = res
                    print(f"{name} berechnet: {res[0]}")
                    felder_dic[name].insert(END, f"{res[0]}")

        if fehlende_pflicht == 1 and fehlende_optional < 3:
            try:
                for name, wert in optionale_werte.items():
                    if wert[0] is not None:
                        res, name = solve(optionale_formeln[name], (v0[0], a[0], g, h[0], wert[0]))

                        for param_name in alle_param_list:
                            if param_name[1] == name:
                                werte[name] = res

                        print(f"{name} berechnet: {res}")
                        felder_dic[name].insert(END, f"{res[0]}")
                        if RuntimeWarning:
                            anzeige_fehler.insert(END, "Physikalisch nicht möglich, heieiei")
            except RuntimeWarning:
                anzeige_fehler.insert(END, "Physikalisch nicht möglich, heieiei")

    except Exception as e:
        print(f"Du hast Mist gebaut, David:{e}")

def flugbahn_mit_widerstand(werte):
    g = 9.80665

    luftdichte = (werte["rho"], "rho")
    cw_wert =  (werte["cw"], "cw")
    querschnitt = (werte["A"], "A")
    masse = (werte["m"], "m")

    h = (werte["h"], "h")
    v0 = (werte["v0"], "v0")
    a = (werte["a"], "a")

    k = (cw_wert[0] * querschnitt[0]  * luftdichte[0]) / (2 * masse[0])
    v_grenz = np.sqrt(g / k)

    h_max = ((v_grenz**2) / (2*g)) * np.log(1+((v0[0]**2) / (v_grenz**2)))

    steigzeit = (v_grenz/g) * np.arctan(v0[0]/v_grenz)
    fallzeit = (v_grenz/g) * np.arccosh(np.exp((h_max * g) / (v_grenz**2)))
# to do vektor save
class tkinterApp(tk.Tk):

    # _init_ function for class tkinterApp
    def __init__(self, *args, **kwargs):
        # _init_ function for class Tk
        tk.Tk.__init__(self, *args, **kwargs)

        # creating a container
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # initializing frames to an empty array
        self.frames = {}

        # iterating through a tuple consisting
        # of the different page layouts
        for F in (StartPage, Page1, Page2, Page3, Page4, Page5, Page6, Page7, Page8, Page9, Page10, Page11, Page12, Page13):
            frame = F(container, self)

            # initializing frame of that object from
            # startpage, page1, page2 respectively with
            # for loop
            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

            self.title("Meth")
            self.geometry("800x400")

        self.show_frame(StartPage)

    # to display the current frame passed as
    # parameter
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


# Startseite
class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        # label of frame Layout 2
        label = ttk.Label(self, text="Übersicht")
        label.grid(row=0, column=4, padx=10, pady=10)

        gleich_sys_button = Button(self, height=2, width=20, text="Gleichungssysteme",
                                   command=lambda: controller.show_frame(Page1))
        gleich_sys_button.grid(row=1, column=1, padx=10, pady=10)

        funktionen_button = Button(self, height=2, width=20, text="Funktionen",
                                   command=lambda: controller.show_frame(Page2))
        funktionen_button.grid(row=2, column=1, padx=10, pady=10)

        vektor_button = Button(self, height=2, width=20, text="Vektoren", command=lambda: controller.show_frame(Page4))
        vektor_button.grid(row=1, column=2, padx=10, pady=10)

        flugbahn_button = Button(self, height=2, width=20, text="Flugbahn",
                                 command=lambda: controller.show_frame(Page11))
        flugbahn_button.grid(row=2, column=2, padx=10, pady=10)


# Gleichungssysteme Übersicht
class Page1(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Gleichungssysteme")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="StartPage", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=1, column=1, padx=10, pady=10)

        zweimalzwei_button = Button(self, height=2, width=20, text="2x2",
                                    command=lambda: controller.show_frame(Page5))
        zweimalzwei_button.grid(row=2, column=1, padx=10, pady=10)

        dreimaldrei_button = Button(self, height=2, width=20, text="3x3",
                                    command=lambda: controller.show_frame(Page6))
        dreimaldrei_button.grid(row=1, column=2, padx=10, pady=10)

        viermalvier_button = Button(self, height=2, width=20, text="4x4",
                                    command=lambda: controller.show_frame(Page7))
        viermalvier_button.grid(row=2, column=2, padx=10, pady=10)


# Funktionen Übersicht
class Page2(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Funktionen")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))

        button1.grid(row=1, column=1, padx=10, pady=10)

        button2 = Button(self, height=2, width=20, text="Kürvendiskussion",
                         command=lambda: controller.show_frame(Page8))

        button2.grid(row=2, column=1, padx=10, pady=10)

        button3 = Button(self, height=2, width=20, text="Integrale", command=lambda: controller.show_frame(Page10))
        button3.grid(row=3, column=1, padx=10, pady=10)


# cheat sheet
class Page3(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Cheat Sheet", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page4))
        button1.grid(row=1, column=1, padx=10, pady=10)

        listbox = Listbox(self, height=7,
                          width=40,
                          bg="light cyan",
                          activestyle='dotbox',
                          font="Arial",
                          fg="black")

        # Define a label for the list.
        label1 = Label(self, text="Übersicht der verschiedenen Operationen")
        label1.grid(row=1, column=4, padx=10, pady=10)
        label2 = Label(self, text="Schreibweise: (grossgeschrieben) von DA-DORT, Operation wird kleingeschrieben")
        label2.grid(row=2, column=4, padx=10, pady=10)

        # insert elements by their
        # index and names.
        listbox.insert(1, "Einzelner Vektor: A-B")
        listbox.insert(2, "Mehrere Vektoren: A-B und B-C")
        listbox.insert(3, "Kreuzprodukt: A-B kp B-C")
        listbox.insert(4, "Skalarprodukt: A-B sp B-C")
        listbox.insert(5, "Ebene: A-B ebene B-C")
        listbox.insert(5, "Winkel: A-B winkel B-C")
        listbox.insert(6, "Vektor speichern: Grossbuchstabe von S bis Z")

        listbox.grid(row=5, column=4, padx=10, pady=10)


# Vektoren
class Page4(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        label = ttk.Label(self, text="Vektoren", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=1, column=1, padx=10, pady=10)

        self.alle_cords = ([
            [],
            [],
            [],
            []
        ])

        # Koeffizienten-Labels für die Vektoren
        labels = ["x", "y", "z", None]

        # vier Punkte Überschrift
        for col in range(4):
            ttk.Label(self,
                      text=f"{["erster Punkt, A", "zweiter Punkt, B", "dritter Punkt, C", "vierter Punkt, D"][col]}").grid(
                row=2, column=2 + col, padx=10, pady=10)
            for row in range(3):
                if labels[col]:  # Falls Text vorhanden, dann labels setzen
                    ttk.Label(self, text=labels[row]).grid(row=row + 3, column=1, padx=10, pady=10, sticky="e")

        # Textfelder für A
        for col in range(4):
            for row in range(3):
                koordinaten = tk.Text(self, height=1, width=5)
                koordinaten.grid(row=row + 3, column=2 + col, padx=10, pady=10, sticky="w")
                einzelner_vektor = self.alle_cords[col]
                einzelner_vektor.append(koordinaten)

        operation = tk.Text(self, height=1, width=15)
        operation.grid(row=3, column=8, padx=10, pady=10, sticky="e")

        ttk.Label(self, text="Operation").grid(row=2, column=8, padx=10, pady=10)

        def solve_and_show():
            # Hilfsfunktion: liest Textfelder -> Strings -> filtert leere -> wandelt in float
            def parse_coords(entries):
                cords = [e.get("1.0", "end-1c").strip() for e in entries]
                nicht_leere = [v for v in cords if v]
                float_werte = [float(v) for v in nicht_leere]
                return cords, float_werte

            # Für alle vier Koordinatengruppen
            results = [parse_coords(group) for group in self.alle_cords]

            # cords = Strings, floats = Floats (jeweils Liste von 4 Listen)
            cords, floats = zip(*results)

            befehl = operation.get("1.0", "end-1c")
            vektorgeometrie(cords, floats, anzeige, befehl)

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=3, column=8, padx=10, pady=10, sticky="w")

        format_button = ttk.Button(self, text="Cheat sheet", command=lambda: controller.show_frame(Page3))
        format_button.grid(row=6, column=1, padx=10, pady=10, sticky="w")

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=4, column=8, padx=10, pady=10)


# 2x2
class Page5(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="2x2", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page1))
        button1.grid(row=1, column=1, padx=10, pady=10)

        # Speichere alle Eingabefelder zur späteren Verwendung (z.B. Auswertung)
        self.entries_2x2 = []

        # Koeffizienten-Labels in der Gleichung
        labels = ["a +", "b =", ""]

        # Zeilen: erste und zweite Gleichung
        for row in range(2):
            ttk.Label(self, text=f"{['Erste', 'Zweite'][row]} Gleichung").grid(row=row + 2, column=1, padx=10, pady=10)

            for col in range(3):
                # Textfeld für Koeffizienten a, b und Ergebnis
                entry = tk.Text(self, height=1, width=5)
                entry.grid(row=row + 2, column=2 + col, padx=10, pady=10, sticky="w")
                self.entries_2x2.append(entry)

                if labels[col]:  # Falls Text vorhanden ist ("" für Ergebnis-Feld)
                    ttk.Label(self, text=labels[col]).grid(row=row + 2, column=2 + col, padx=10, pady=10, sticky="e")

        def solve_and_show():
            values = [e.get("1.0", "end-1c").strip() for e in self.entries_2x2]
            try:
                if len(values) == 6:
                    int_values = [int(v) for v in values]
                    a1, b1, r1, a2, b2, r2 = int_values
                    zweimalzwei(a1, b1, a2, b2, r1, r2, anzeige)
            except ValueError:
                anzeige.insert(END, f"Ein Feld ist noch leer")

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)


# 3x3
class Page6(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="3x3", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page1))
        button1.grid(row=1, column=1, padx=10, pady=10)

        # Speichere alle Eingabefelder zur späteren Verwendung (z.B. Auswertung)
        self.entries_3x3 = []

        # Koeffizienten-Labels in der Gleichung
        labels = ["a +", "b +", "c =", ""]

        # Zeilen: erste und zweite Gleichung
        for row in range(3):
            ttk.Label(self, text=f"{['Erste', 'Zweite', 'Dritte'][row]} Gleichung").grid(row=row + 2, column=1, padx=10,
                                                                                         pady=10)

            for col in range(4):
                # Textfeld für Koeffizienten a, b und Ergebnis
                entry = tk.Text(self, height=1, width=5)
                entry.grid(row=row + 2, column=2 + col, padx=10, pady=10, sticky="w")
                self.entries_3x3.append(entry)

                if labels[col]:  # Falls Text vorhanden ist ("" für Ergebnis-Feld)
                    ttk.Label(self, text=labels[col]).grid(row=row + 2, column=2 + col, padx=10, pady=10, sticky="e")

        def solve_and_show():
            values = [e.get("1.0", "end-1c").strip() for e in self.entries_3x3]
            try:
                if len(values) == 12:
                    int_values = [int(v) for v in values]
                    a1, b1, c1, r1, a2, b2, c2, r2, a3, b3, c3, r3 = int_values
                    dreimaldrei(a1, b1, c1, a2, b2, c2, a3, b3, c3, r1, r2, r3, anzeige)
            except ValueError:
                anzeige.insert(END, f"Ein Feld ist noch leer")

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)


# 4x4
class Page7(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="4x4", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page1))
        button1.grid(row=1, column=1, padx=10, pady=10)

        # Speichere alle Eingabefelder zur späteren Verwendung (z.B. Auswertung)
        self.entries_4x4 = []

        # Koeffizienten-Labels in der Gleichung
        labels = ["a +", "b +", "c +", "d =", ""]

        # Zeilen: erste, zweite usw. Gleichung
        for row in range(4):
            ttk.Label(self, text=f"{['Erste', 'Zweite', 'Dritte', 'Vierte'][row]} Gleichung").grid(row=row + 2,
                                                                                                   column=1, padx=10,
                                                                                                   pady=10)

            for col in range(5):
                # Textfeld für Koeffizienten a, b und Ergebnis
                entry = tk.Text(self, height=1, width=10)
                entry.grid(row=row + 2, column=2 + col, padx=10, pady=10, sticky="w")
                self.entries_4x4.append(entry)

                if labels[col]:  # Falls Text vorhanden ist ("" für Ergebnis-Feld)
                    ttk.Label(self, text=labels[col]).grid(row=row + 2, column=2 + col, padx=10, pady=10, sticky="e")

        def solve_and_show():
            values = [e.get("1.0", "end-1c").strip() for e in self.entries_4x4]
            try:
                if len(values) == 20:
                    int_values = [int(v) for v in values]
                    a1, b1, c1, d1, r1, a2, b2, c2, d2, r2, a3, b3, c3, d3, r3, a4, b4, c4, d4, r4 = int_values
                    viermalvier(a1, b1, c1, d1, r1, a2, b2, c2, d2, r2, a3, b3, c3, d3, r3, a4, b4, c4, d4, r4, anzeige)
            except ValueError:
                anzeige.insert(END, f"Ein Feld ist noch leer")

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)


# kurvendiskussion
class Page8(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Kürven diskutieren", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page2))
        button1.grid(row=1, column=1, padx=10, pady=10)

        funktion_label = ttk.Label(self, text="Funktion")
        funktion_label.grid(row=2, column=1, padx=10, pady=10)
        self.funktion = tk.Text(self, height=3, width=20)
        self.funktion.grid(row=3, column=1, padx=10, pady=10)

        n_ableitung_label = ttk.Label(self, text="Anzahl Ableitungen")
        n_ableitung_label.grid(row=2, column=4, padx=10, pady=10)
        self.n_ableitung = tk.Entry(self)
        self.n_ableitung.grid(row=3, column=4, padx=10, pady=10)

        def solve_and_show():
            funktion = self.funktion.get("1.0", "end-1c")
            n_ableitung_int = int(self.n_ableitung.get())

            funktion_expr = sympify(funktion)
            anzeige.delete(1.0, END)
            ableitung(funktion_expr, n_ableitung_int, anzeige)
            kurvendiskussion(funktion_expr, anzeige)
            return funktion_expr

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=5, padx=10, pady=10)

        anzeige = Text(self, height=10, width=40, bg="light cyan")
        anzeige.grid(row=3, column=5, padx=10, pady=10)

        def anplot():
            try:
                funktion_expr = solve_and_show()
                # Zugriff auf Page8 über den Controller
                page9 = self.controller.frames[Page9]
                self.controller.show_frame(Page9)
                d_oder_i = "diff"
                page9.plot(funktion_expr, d_oder_i)
            except Exception as e:
                print(f"Fehler beim Plotten: {e}")

        plot_button = ttk.Button(self, text="Plot", command=anplot)
        plot_button.grid(row=4, column=1, padx=10, pady=10)


# funktionsgraph
class Page9(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        label = ttk.Label(self, text="Graph", font=("Arial",15))
        label.pack()

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))
        button1.pack()

        # Frame für den Plot
        self.plot_frame = tk.Frame(self)
        self.plot_frame.pack(fill="both", expand=True)

    def plot(self, funktion_expr, d_oder_i, *args, **kwargs):
        # Clear previous plot if it exists
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Create a new figure
        fig = Figure(figsize=(5, 2), dpi=100)
        plot1 = fig.add_subplot(111)

        x_sym = Symbol('x')
        f = lambdify(x_sym, funktion_expr, 'numpy')

        x_vals = np.linspace(-30, 30, 400)
        try:
            y_vals = f(x_vals)
        except:
            # Falls die Funktion nicht überall definiert ist
            y_vals = np.vectorize(lambda x: f(x) if x != 0 else np.nan)(x_vals)

        if d_oder_i == "int":
            ix = np.linspace(args, kwargs)
            iy = f(ix)
            verts = [(args, 0), *zip(ix, iy), (kwargs, 0)]
            poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
            plot1.add_patch(poly)

            plot1.text(0.5 * (args + kwargs), 30, r"$\int_a^b f(x)\mathrm{d}x$",
                       horizontalalignment='center', fontsize=20)
        # Plot the function
        plot1.plot(x_vals, y_vals)
        plot1.grid(True)
        plot1.set_title(f"Plot von {funktion_expr}")

        # Embed the plot in Tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.plot_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


# Integrale
class Page10(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Integrale", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page2))
        button1.grid(row=1, column=1, padx=10, pady=10)

        funktion_label = ttk.Label(self, text="Funktion")
        funktion_label.grid(row=2, column=1, padx=10, pady=10)
        self.funktion = tk.Text(self, height=3, width=20)
        self.funktion.grid(row=3, column=1, padx=10, pady=10)

        n_aufleitung_label = ttk.Label(self, text="Anzahl Integration")
        n_aufleitung_label.grid(row=2, column=3, padx=10, pady=10)
        self.n_aufleitung = tk.Entry(self)
        self.n_aufleitung.grid(row=3, column=3, padx=10, pady=10)

        range_label = ttk.Label(self, text="Von wo bis wo, Format: (a,b)")
        range_label.grid(row=2, column=4, padx=10, pady=10)
        self.range = tk.Entry(self)
        self.range.grid(row=3, column=4, padx=10, pady=10)

        def solve_and_show():
            funktion = self.funktion.get("1.0", "end-1c")
            n_aufleitung_int = int(self.n_aufleitung.get())
            range_str = self.range.get().strip("()")
            a_str, b_str = range_str.split(',')
            a = int(a_str.strip())
            b = int(b_str.strip())

            funktion_expr = sympify(funktion)

            aufleitung(funktion_expr, n_aufleitung_int, anzeige)
            integrale(funktion_expr, anzeige, a, b)
            return funktion_expr, a, b

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=5, padx=10, pady=10)

        anzeige = Text(self, height=10, width=40, bg="light cyan")
        anzeige.grid(row=3, column=5, padx=10, pady=10)

        def anplot():
            try:
                funktion_expr, a, b = solve_and_show()
                # Zugriff auf Page8 über den Controller
                page9 = self.controller.frames[Page9]
                self.controller.show_frame(Page9)
                d_oder_i = "int"
                page9.plot(funktion_expr, d_oder_i, a, b)
            except Exception as e:
                print(f"Fehler beim Plotten: {e}")

        plot_button = ttk.Button(self, text="Plot", command=anplot)
        plot_button.grid(row=4, column=1, padx=10, pady=10)


# Flugbahn
class Page11(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Flugbahn", font=("Arial",15))
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=1, column=1, padx=10, pady=10)

        labels = ["h", "v0", "a", "t", "R", "h_max"]
        felder_dic = {}

        for col in range(1, 7):
            bezeichnung = tk.Label(self, text=labels[col - 1])
            bezeichnung.grid(row=3, column=1 + col, padx=10, pady=10, sticky="w")
            anzeige_res = tk.Text(self, height=1, width=10, bg="light cyan")
            anzeige_res.grid(row=4, column=1 + col, padx=10, sticky="w")
            felder_dic.update({labels[col - 1]: anzeige_res})

        anzeige_fehler = Text(self, height=2, width=20, bg="light cyan")
        anzeige_fehler.grid(row=1, column=4)

        label = ttk.Label(self, text="Fehlermeldung-Box")
        label.grid(row=1, column=3, padx=10, pady=10)

        def solve_and_show():
            werte = {}
            for key, widget in felder_dic.items():
                eintrag_feld = widget.get("1.0", "end-1c")
                if eintrag_feld != '':
                    eintrag_feld = float(eintrag_feld)
                if eintrag_feld == '':
                    eintrag_feld = None
                werte[key] = eintrag_feld
            flugbahn_ohne_widerstand(werte, anzeige_fehler, felder_dic)

        preset_button = ttk.Button(self, text="presets", command=lambda: controller.show_frame(Page12))
        preset_button.grid(row=5, column=1, padx=10, pady=10)

        spulen_flugbahn = Button(self, height=2, width=20, text="Flugbahn aus Spulen", command=lambda: controller.show_frame(Page13))
        spulen_flugbahn.grid(row=6, column=1, padx=10, pady=10)

        solv_button = ttk.Button(self, text="solv", command=solve_and_show)
        solv_button.grid(row=1, column=5, padx=10, pady=10)

        def clear():
            for name in labels:
                felder_dic[name].delete("1.0", END)
            anzeige_fehler.delete("1.0", END)

        clear_button = ttk.Button(self, text="clear", command=clear)
        clear_button.grid(row=1, column=6, padx=10, pady=10)


class Page12(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Presets")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page11))
        button1.grid(row=1, column=1, padx=10, pady=10)

        listbox = Listbox(self, height=9,
                          width=40,
                          bg="light cyan",
                          activestyle='dotbox',
                          font="Arial",
                          fg="black")

        # Define a label for the list.
        label1 = Label(self, text="Übersicht der verschiedenen Einstellungen")
        label1.grid(row=1, column=4, padx=10, pady=10)

        # insert elements by their
        # index and names.
        listbox.insert(1, "Speer") # 0.8kg, Durchmesser 0.03m, cw = 0.5
        listbox.insert(2, "Football") # 0.4252kg, Durchmesser 0.54m, cw = 0.25
        listbox.insert(3, "Die Leiden des jungen Werthers (offen)") #0.204kg, 0.248m (offen), cw = 1.17
        listbox.insert(4, "T-34/85 Schuss")# 9.21 kg, Durchmesser 0.85 m, v0 = 792 m/s, cw = 0.3
        listbox.insert(6, "Medizinball")# 5kg, Durchmesser 0.3m, cw = 0.47
        listbox.insert(8, "Stuhl, Mensa")# 4.9 kg, 0.210375m^2, cw = 1.5
        listbox.insert(9, "Lehrer") # 80kg, 0.95m^2, cw = 1.15

        listbox.grid(row=5, column=4, padx=10, pady=10)

        def selected_item():

            for i in listbox.curselection():
                print(listbox.get(i))

        select_knopf = ttk.Button(self, text="auswählen", command=selected_item)
        select_knopf.grid(row=11, column=1, padx=10, pady=10)
        # to do, dass presets funktionieren

class Page13(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        label = ttk.Label(self, text="Beschleunigte Kugeln aus Spulen und deren Flugbahn", font=("Arial",15))
        label.grid(row=0, column=9)

        label_fehler = ttk.Label(self, text="Fehlermeldung-Box")
        label_fehler.grid(row=2, column=2, padx=10, pady=10, sticky="w")
        self.anzeige_fehler = Text(self, height=2, width=20, bg="light cyan")
        self.anzeige_fehler.grid(row=2, column=3)

        button1 = ttk.Button(self, text="zurück", command=lambda: controller.show_frame(Page11))
        button1.grid(row=2, column=1, padx=10, pady=10)

        self.felder_dic = {}

        # entry-felder
        labels = ["r_kugel", "u_relativ", "N_windungen", "I", "r_spule", "Länge"]

        for col in range(1, 7):
            bezeichnung = tk.Label(self, text=labels[col - 1])
            bezeichnung.grid(row=4, column=col, padx=5, pady=5, sticky="w")
            anzeige_res = tk.Text(self, height=1, width=10, bg="light cyan")
            anzeige_res.grid(row=5, column=col, padx=5, sticky="w")

            anzeige_res.bind("<Return>", self.berechne_und_plotte)
            anzeige_res.bind("<FocusOut>", self.berechne_und_plotte)

            self.felder_dic.update({labels[col - 1]: anzeige_res})
        print(self.felder_dic)

        # ---------- Plot-Bereich ----------

        # ---------- Kraft auf die Kugel Diagramm --------
        label_plot1 = ttk.Label(self, text="Kraft auf die Kugel")
        label_plot1.grid(row=3, column=9)

        self.plot_frame_1 = tk.Frame(self)
        self.plot_frame_1.grid(row=4, column=9, rowspan=5,columnspan=5)

        self.fig_1 = Figure(figsize=(5, 2), dpi=100)
        self.ax_1 = self.fig_1.add_subplot(111)

        self.canvas_1 = FigureCanvasTkAgg(self.fig_1, master=self.plot_frame_1)
        self.canvas_1.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar_1 = NavigationToolbar2Tk(self.canvas_1, self.plot_frame_1)
        self.toolbar_1.update()

        # --------- Flugbahn Diagramm -----------
        label_plot2 = ttk.Label(self, text="Flugbahn der Kugel")
        label_plot2.grid(row=10, column=9)

        self.plot_frame_2 = tk.Frame(self)
        self.plot_frame_2.grid(row=11, column=9, rowspan=5, columnspan=5)

        self.fig_2 = Figure(figsize=(5, 2), dpi=100)
        self.ax_2 = self.fig_2.add_subplot(111)

        self.canvas_2 = FigureCanvasTkAgg(self.fig_2, master=self.plot_frame_2)
        self.canvas_2.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar_2 = NavigationToolbar2Tk(self.canvas_2, self.plot_frame_2)
        self.toolbar_2.update()

    def berechne_und_plotte(self, event=None):
        werte = {}

        for key, widget in self.felder_dic.items():
            text = widget.get("1.0", "end-1c").strip()
            werte[key] = float(text) if text else None

        self.funktion = spulen_beschleunigung(werte, self.anzeige_fehler, self.felder_dic)

        if self.funktion is not None:
            self.update_plot()

        # ---------- Automatisches Update ----------

    def update_plot(self):
        if self.funktion is None:
            return

        self.plot(self.ax_1, self.canvas_1, self.funktion)
        """self.plot(self.ax_2, self.canvas_2, self.funktion)"""

        # ---------- Plot-Funktion ----------

    def plot(self, ax, canvas, funktion_expr):
        ax.clear()
        funktion_expr_sympy = sympy.sympify(funktion_expr)
        print(funktion_expr_sympy)
        x = Symbol("x")
        f = lambdify(x, funktion_expr_sympy, 'numpy')

        x_vals = np.linspace(-30, 30, 400)

        try:
            y_vals = f(x_vals)
        except:
            y_vals = np.full_like(x_vals, np.nan)

        ax.plot(x_vals, y_vals)
        ax.grid(True)
        ax.set_title(f"Plot von {funktion_expr_sympy}")

        canvas.draw()


# Driver Code
app = tkinterApp()
app.mainloop()
