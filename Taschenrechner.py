import numpy as np
from tkinter import *
import tkinter.ttk as ttk
import tkinter as tk
from sympy import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


def zweimalzwei(a1, b1, a2, b2, c1, c2, anzeige):
    y = ((a2 * c1) - (a1 * c2)) / ((a2 * b1) - (a1 * b2))
    x = (c1 - b1 * y) / a1
    anzeige.insert(END, f"a ist {x} und b ist {y}")


def dreimaldrei(a1, b1, c1, a2, b2, c2, a3, b3, c3, d1, d2, d3, anzeige):
    z = (a1 * (b2 * d3 - b3 * d2) - b1 * (a2 * d3 - a3 * d2) + d1 * (a2 * b3 - a3 * b2)) / (
            a1 * (b2 * c3 - b3 * c2) - b1 * (a2 * c3 - a3 * c2) + c1 * (a2 * b3 - a3 * b2))
    y = (a1 * (d2 * c3 - d3 * c2) - d1 * (a2 * c3 - a3 * c2) + c1 * (a2 * d3 - a3 * d2)) / (
            a1 * (b2 * c3 - b3 * c2) - b1 * (a2 * c3 - a3 * c2) + c1 * (a2 * b3 - a3 * b2))
    x = (d1 * (b2 * c3 - b3 * c2) - b1 * (d2 * c3 - d3 * c2) + c1 * (d2 * b3 - d3 * b2)) / (
            a1 * (b2 * c3 - b3 * c2) - b1 * (a2 * c3 - a3 * c2) + c1 * (a2 * b3 - a3 * b2))
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
        abg_funktion_2 = ableitung(funktion,2,"a")
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
        abg_funktion_2 = ableitung(funktion,2,"a")

        x_res_liste = solve(abg_funktion_2, x)
        for anzahl, x_res in enumerate(x_res_liste):
            y_res = funktion.subs(x, x_res)
            anzeige.insert(END, f"{anzahl+1}. Wendepunkt bei ({x_res}und{y_res})" + "\n")

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
        for F in (StartPage, Page1, Page2, Page3, Page4, Page5, Page6, Page7, Page8):
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

        button2 = Button(self, height=2, width=20, text="Kurvendiskussion",
                             command=lambda: controller.show_frame(Page7))

        button2.grid(row=2, column=1, padx=10, pady=10)


# Vektoren Übersicht
class Page3(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Vektoren")
        label.grid(row=0, column=4, padx=10, pady=10)

        button1 = ttk.Button(self, text="Startpage", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=1, column=1, padx=10, pady=10)


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


class Page7(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller  # Speichern des Controllers als Instanzvariable
        label = ttk.Label(self, text="Kurven diskutieren")
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
                page8.plot(funktion_expr)
            except Exception as e:
                print(f"Fehler beim Plotten: {e}")

        plot_button = ttk.Button(self, text="Plot", command=anplot)
        plot_button.grid(row=4, column=2, padx=10, pady=10)


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

    def plot(self, funktion_expr):
        # Clear previous plot if it exists
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Create a new figure
        fig = Figure(figsize=(5, 2), dpi=100)
        plot1 = fig.add_subplot(111)

        # Convert SymPy function to a NumPy-compatible function
        x_sym = Symbol('x')
        f = lambdify(x_sym, funktion_expr, 'numpy')

        # Generate x values
        x_vals = np.linspace(-30, 30, 400)
        try:
            y_vals = f(x_vals)
        except:
            # Falls die Funktion nicht überall definiert ist
            y_vals = np.vectorize(lambda x: f(x) if x != 0 else np.nan)(x_vals)

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


# Driver Code
app = tkinterApp()
app.mainloop()

"""gleichungssystemswitch = input("Ein Gleichungssystem lösen? (ja oder nein)")
if gleichungssystemswitch.lower() == "ja":
    zweimalzweiswitch = input("Ein 2x2 Gleichungssystem lösen? (ja oder nein)")
    if zweimalzweiswitch.lower() == "ja":
        zweimalzwei()
    dreimaldreiswitch = input("Ein 3x3 Gleichungssystem lösen? (ja oder nein)")
    if dreimaldreiswitch.lower() == "ja":
        dreimaldrei()
    viermalvierswitch = input("Ein 4x4 Gleichungssystem lösen? (ja oder nein)")
    if viermalvierswitch.lower() == "ja":
        viermalvier()

nullstellenswitch = input("Die Nullstellen einer Funktion berechnen? (ja oder nein)")
if nullstellenswitch.lower() == "ja":       #bei dem gibt es irgendwo noch Fehler
    a = int(input("Was ist a"))
    b = int(input("Was ist b"))
    c = int(input("Was ist c"))
    d = int(input("Was ist d"))
    k = int(input("Was ergibt diese Funktion"))

    e = d-k


    fig, ax = plt.subplots(figsize=(15, 6))
    ax.set(xlim=(-50, 50), ylim=(-50, 50), aspect='equal')


    # Add coordinate axes
    ax.axhline(0, color='black', linewidth=0.5)  # x-axis
    ax.axvline(0, color='black', linewidth=0.5)  # y-axis



    if a > 0 or a < 0:
        kubisch(a,b,c,e)
        x = np.linspace(-10, 10, 400)
        y = (a * (x ** 3)) + (b * (x ** 2)) + (c * x) + d
    elif 1*(10*-10) < a < 1(10**-9):
        quadratisch(b,c,d)
        x = np.linspace(-10, 10, 400)
        y = (b * (x ** 2)) + (c * x) + e
    elif 1*(10*-10) < a < 1(10*-9) and 1(10*-10) < b < 1(10**-9):
        linear(k, c, d)
        x = np.linspace(-10, 10, 400)
        y = (c * x) + d

    plt.plot(x, y)
    plt.tight_layout()
    plt.show()

ebenen_und_kreuz_switch = input("Willst du die EBENE und/oder das KREUZPRODUKT wissen? (AB und AC Vektoren ebenfalls enthalten), ja oder nein")
if ebenen_und_kreuz_switch.lower() == "ja":
    Ax = int(input("x von a"))
    Ay = int(input("y von a"))
    Az = int(input("z von a"))

    Bx = int(input("x von b"))
    By = int(input("y von b"))
    Bz = int(input("z von b"))

    Cx = int(input("x von c"))
    Cy = int(input("y von c"))
    Cz = int(input("z von c"))

    AB_switch = input("willst du AB")
    AC_switch = input("willst du AC")

    KREUZ_switch = input("willst du das Kreuzprodukt")
    A = (Ax, Ay, By)
    B = (Bx, By, Bz)
    C = (Cx, Cy, Cz)

    ABx = Bx - Ax
    ABy = By - Ay
    ABz = Bz - Az

    ACx = Cx - Ax
    ACy = Cy - Ay
    ACz = Cz - Az

    AB = (ABx, ABy, ABz)
    AC = (ACx, ACy, ACz)

    if AB_switch.lower() == "ja":
        print(AB, "das ist der AB Vektor")
    if AC_switch.lower() == "ja":
        print(AC, "das ist der AC Vektor")


    KREUZ = ((ABy * ACz) - (ABz * ACy), (ABz * ACx) - (ABx * ACz), (ABx * ACy) - (ABy * ACx))

    if KREUZ_switch.lower() == "ja":
        print(KREUZ, "das ist das KREUZPRODUKT")

    KREUZx = (ABy * ACz) - (ABz * ACy)
    KREUZy = (ABz * ACx) - (ABx * ACz)
    KREUZz = (ABx * ACy) - (ABy * ACx)

    d = -((KREUZx*Ax)+(KREUZy*Ay)+(KREUZz*Az))

    print("Ebene ist", KREUZx, "x+", KREUZy, "y+", KREUZz, "z+", d)"""