import numpy as np

class pressure:

    def __init__(self, xSteps, dx, dt, c):
        self.Pressure = np.zeros(xSteps)    # P now
        self.PressureOld = np.zeros(xSteps) # P one previous timestep
        self.PressureNew = np.zeros(xSteps) # P in new timestep
        self.PressureSecDer = np.zeros(xSteps) # P second derivative
        self.xSteps = xSteps
        self.dx = dx
        self.dt = dt
        self.c = c

    def partialderivative(self): # 3 point finite difference scheme
        for i in range(1, self.xSteps - 1):
            self.PressureSecDer[i] = (self.Pressure[i + 1] - 2 * self.Pressure[i] + self.Pressure[i - 1]) / self.dx ** 2

    def timeextrapolation(self):
        self.PressureNew = 2* self.Pressure - self.PressureOld + self.c ** 2 * self.dt ** 2 * self.PressureSecDer

