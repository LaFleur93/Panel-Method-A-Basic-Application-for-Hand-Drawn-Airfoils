import matplotlib.pyplot as plt
from matplotlib import path
import scipy.interpolate as interpolate
from numpy import *

class LineBuilder:
    def __init__(self, line,ax):
        self.line = line
        self.ax = ax
        self.xs = []
        self.ys = []
        self.nodes = []
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        self.counter = 0
        self.precision = 0.02

    def __call__(self, event):
        if event.inaxes != self.line.axes: return
        if self.counter == 0:
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            self.nodes.append((event.xdata, event.ydata))
        if abs(event.xdata-self.xs[0]) <= self.precision and abs(event.ydata-self.ys[0]) <= self.precision and self.counter != 0:
            self.xs.append(self.xs[0])
            self.ys.append(self.ys[0])
            self.nodes.append((self.xs[0], self.ys[0]))
            self.ax.scatter(self.xs,self.ys, s=30, edgecolor='black', linewidth=0, facecolor='black')
            self.ax.scatter(self.xs[0], self.ys[0], s=30, edgecolor='black', linewidth=1, facecolor='red')
            self.ax.plot(self.xs, self.ys, color= 'black')
            self.line.figure.canvas.draw()
            self.counter = 0
        else:
            if self.counter != 0:
                self.xs.append(event.xdata)
                self.ys.append(event.ydata)
                self.nodes.append((event.xdata, event.ydata))
            self.ax.scatter(self.xs, self.ys, s=30, edgecolor='black', linewidth=1, facecolor='gray')
            self.ax.plot(self.xs, self.ys, color='black')
            self.line.figure.canvas.draw()
            self.counter = self.counter + 1