import numpy as np
from tkinter import *
import tkinter.ttk as ttk
import tkinter as tk
from itertools import permutations
import sympy
from sympy import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.patches import Polygon
import re


def zweimalzwei(a1, b1, a2, b2, c1, c2, anzeige):
    y = ((a2 * c1) - (a1 * c2)) / ((a2 * b1) - (a1 * b2))
    x = (c1 - b1 * y) / a1
    anzeige.insert(END, f"a ist {x} und b ist {y}")


def dreimaldrei(a1, b1, c1, a2, b2, c2, a3, b3, c3, d1, d2, d3, anzeige):
    rest = (a1 * (b2 * c3 - b3 * c2) - b1 * (a2 * c3 - a3 * c2) + c1 * (a2 * b3 - a3 * b2))
    z = (a1 * (b2 * d3 - b3 * d2) - b1 * (a2 * d3 - a3 * d2) + d1 * (a2 * b3 - a3 * b2)) / rest
    y = (a1 * (d2 * c3 - d3 * c2) - d1 * (a2 * c3 - a3 * c2) + c1 * (a2 * d3 - a3 * d2)) / rest
    x = (d1 * (b2 * c3 - b3 * c2) - b1 * (d2 * c3 - d3 * c2) + c1 * (d2 * b3 - d3 * b2)) / rest
    anzeige.insert(END, f"a ist {x} ,b ist {y} und c ist {z}")


def viermalvier(a1, b1, c1, d1, r1, a2, b2, c2, d2, r2, a3, b3, c3, d3, r3, a4, b4, c4, d4, r4, anzeige):
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
        print(-1 * funktion)
        if funktion == funktion.subs(x, -x):
            anzeige.insert(END, "Symmetrisch zur y-Achse" + "\n")
        elif funktion == -1 * (funktion.subs(x, -x)):
            anzeige.insert(END, "Symmetrisch zum Ursprung" + "\n")

    hoch_tief(funktion, anzeige)
    wendepunkt(funktion, anzeige)
    symmetrien(funktion, anzeige)


def aufleitung(funktion_expr, n_aufleitung_int, anzeige):
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
    x = Symbol("x")
    area = sympy.integrate(funktion, (x, a, b))
    anzeige.insert(END, f"Die Fläche unter der Funktion{funktion} ist:{area}" + "\n")


def vektorgeometrie(values, float_values, anzeige, befehl):
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

    for key in vektoren_kombi:
        if key == befehl:
            ein_vektor = key
            print(vektoren_kombi[key])
            anzeige_key = "".join(key.split("-"))
            anzeige.insert(END, f"{anzeige_key} ist {vektoren_kombi[key]}\n")

    def operation_auslesen():
        resultat = re.search(r"(\b[a-z]+\b)", string=befehl)
        if resultat:
            operation = resultat.groups(1)
            anzeige.insert(END, f"{operation}\n")

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

    def ebene(vektor1, vektor2):
        kreuzprodukt = np.cross(vektor1, vektor2)
        d = -np.dot(kreuzprodukt, 0)
        try:
            ebene_zwei = [f"x:{kreuzprodukt[0]} y:{kreuzprodukt[1]} z:{kreuzprodukt[2]} + {d}"]
            return ebene_zwei
        except IndexError:
            print("Keine z-Koordinate")



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
        for F in (StartPage, Page1, Page2, Page3, Page4, Page5, Page6, Page7, Page8, Page9):
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

        vektor_button = Button(self, height=2, width=20, text="Vektoren", command=lambda: controller.show_frame(Page3))
        vektor_button.grid(row=1, column=2, padx=10, pady=10)


# Gleichungssysteme Übersicht
class Page1(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Gleichungssysteme")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="StartPage", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=1, column=1, padx=10, pady=10)

        zweimalzwei_button = Button(self, height=2, width=20, text="2x2",
                                    command=lambda: controller.show_frame(Page4))
        zweimalzwei_button.grid(row=2, column=1, padx=10, pady=10)

        dreimaldrei_button = Button(self, height=2, width=20, text="3x3",
                                    command=lambda: controller.show_frame(Page5))
        dreimaldrei_button.grid(row=1, column=2, padx=10, pady=10)

        viermalvier_button = Button(self, height=2, width=20, text="4x4",
                                    command=lambda: controller.show_frame(Page6))
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
                         command=lambda: controller.show_frame(Page7))

        button2.grid(row=2, column=1, padx=10, pady=10)

        button3 = Button(self, height=2, width=20, text="Integrale", command=lambda: controller.show_frame(Page9))
        button3.grid(row=3, column=1, padx=10, pady=10)


# Vektoren
class Page3(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Vektoren")
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

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=4, column=8, padx=10, pady=10)


# 2x2
class Page4(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="2x2")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startseite", command=lambda: controller.show_frame(StartPage))
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
            if len(values) == 6:
                int_values = [int(v) for v in values]
                a1, b1, r1, a2, b2, r2 = int_values
                zweimalzwei(a1, b1, a2, b2, r1, r2, anzeige)

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)


# 3x3
class Page5(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="3x3")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startseite", command=lambda: controller.show_frame(StartPage))
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
            if len(values) == 12:
                int_values = [int(v) for v in values]
                a1, b1, c1, r1, a2, b2, c2, r2, a3, b3, c3, r3 = int_values
                dreimaldrei(a1, b1, c1, a2, b2, c2, a3, b3, c3, r1, r2, r3, anzeige)

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)


# 4x4
class Page6(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="4x4")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))
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
            if len(values) == 20:
                int_values = [int(v) for v in values]
                a1, b1, c1, d1, r1, a2, b2, c2, d2, r2, a3, b3, c3, d3, r3, a4, b4, c4, d4, r4 = int_values
                viermalvier(a1, b1, c1, d1, r1, a2, b2, c2, d2, r2, a3, b3, c3, d3, r3, a4, b4, c4, d4, r4, anzeige)

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=5, width=25, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)


# kurvendiskussion
class Page7(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Kürven diskutieren")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=1, column=1, padx=10, pady=10)

        funktion_label = ttk.Label(self, text="Funktion")
        funktion_label.grid(row=2, column=1, padx=10, pady=10)
        self.funktion = tk.Text(self, height=3, width=20)
        self.funktion.grid(row=3, column=1, padx=10, pady=10)

        n_ableitung_label = ttk.Label(self, text="Anzahl Ableitungen")
        n_ableitung_label.grid(row=2, column=5, padx=10, pady=10)
        self.n_ableitung = tk.Entry(self)
        self.n_ableitung.grid(row=3, column=5, padx=10, pady=10)

        def solve_and_show():
            funktion = self.funktion.get("1.0", "end-1c")
            n_ableitung_int = int(self.n_ableitung.get())

            funktion_expr = sympify(funktion)

            ableitung(funktion_expr, n_ableitung_int, anzeige)
            kurvendiskussion(funktion_expr, anzeige)
            return funktion_expr

        values_button = ttk.Button(self, text="solv", command=solve_and_show)
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=10, width=40, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)

        def anplot():
            try:
                funktion_expr = solve_and_show()
                # Zugriff auf Page8 über den Controller
                page8 = self.controller.frames[Page8]
                self.controller.show_frame(Page8)
                d_oder_i = "diff"
                page8.plot(funktion_expr, d_oder_i)
            except Exception as e:
                print(f"Fehler beim Plotten: {e}")

        plot_button = ttk.Button(self, text="Plot", command=anplot)
        plot_button.grid(row=4, column=2, padx=10, pady=10)


# funktionsgraph
class Page8(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        label = ttk.Label(self, text="Graph")
        label.pack()

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))
        button1.pack()

        # Frame für den Plot
        self.plot_frame = tk.Frame(self)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

    def plot(self, funktion_expr, d_oder_i, *args, **kwargs):
        # Clear previous plot if it exists
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Create a new figure
        fig = Figure(figsize=(5, 2), dpi=100)
        plot1 = fig.add_subplot(111)

        # Convert SymPy function to a NumPy-compatible function
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


class Page9(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Integrale")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=1, column=1, padx=10, pady=10)

        funktion_label = ttk.Label(self, text="Funktion")
        funktion_label.grid(row=2, column=1, padx=10, pady=10)
        self.funktion = tk.Text(self, height=3, width=20)
        self.funktion.grid(row=3, column=1, padx=10, pady=10)

        n_aufleitung_label = ttk.Label(self, text="Anzahl Aufleitungen")
        n_aufleitung_label.grid(row=2, column=5, padx=10, pady=10)
        self.n_aufleitung = tk.Entry(self)
        self.n_aufleitung.grid(row=3, column=5, padx=10, pady=10)

        range_label = ttk.Label(self, text="Von wo bis wo, Format: (a,b)")
        range_label.grid(row=2, column=6, padx=10, pady=10)
        self.range = tk.Entry(self)
        self.range.grid(row=3, column=6, padx=10, pady=10)

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
        values_button.grid(row=2, column=8, padx=10, pady=10)

        anzeige = Text(self, height=10, width=40, bg="light cyan")
        anzeige.grid(row=3, column=8, padx=10, pady=10)

        def anplot():
            try:
                funktion_expr, a, b = solve_and_show()

                # Zugriff auf Page8 über den Controller
                page8 = self.controller.frames[Page8]
                self.controller.show_frame(Page8)
                d_oder_i = "int"
                page8.plot(funktion_expr, d_oder_i, a, b)
            except Exception as e:
                print(f"Fehler beim Plotten: {e}")

        plot_button = ttk.Button(self, text="Plot", command=anplot)
        plot_button.grid(row=4, column=2, padx=10, pady=10)


# Driver Code
app = tkinterApp()
app.mainloop()
