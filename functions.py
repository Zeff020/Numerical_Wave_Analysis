import numpy as np

# Parent class with default constructor
class function:

    def __init__(self, SourceFreq, TimeSteps, dt):

        self.SourceFreq = SourceFreq # Frequency of the source
        self.Shift = 4./SourceFreq # Time at source initialization (shift)
        self.TimeSteps = TimeSteps # Number of timesteps
        self.dt = dt # Timestep size
        self.Source  = np.zeros(TimeSteps + 1)
        self.Time = np.linspace(0 * dt, TimeSteps * dt, TimeSteps)

# First derivative of Gaussian function
class gaussianderivative(function):

    def __init__(self, SourceFreq, TimeSteps, dt):
        super().__init__(SourceFreq, TimeSteps, dt) # Use constructor from parent class
        self.Source = -2. * (self.Time - self.Shift) * (self.SourceFreq ** 2) * (np.exp(-1.0 * (self.SourceFreq ** 2) * (self.Time - self.Shift) ** 2))

    def returnfunction(self):
        return self.Source, self.Time, self.Shift

# Gaussian function
class gaussian(function):

    def __init__(self, SourceFreq, TimeSteps, dt, SourceFunc,ir,isrc,WaveSpeed, dx, xpos):
        super().__init__(SourceFreq, TimeSteps, dt)
        self.WaveSpeed = WaveSpeed
        self.ir = ir
        self.isrc = isrc
        self.SourceFunc =SourceFunc
        self.dx = dx
        self.xpos = xpos
        self.G    = self.Time * 0.

    def returnfunction(self): # Method that calculates and returns the gaussian as an integral of the derivative of the Gaussian using a Green's function
        for it in range(self.TimeSteps):
            if (self.Time[it] - np.abs(self.xpos[self.ir] - self.xpos[self.isrc]) / self.WaveSpeed) >= 0:
                self.G[it] = 1. / (2 * self.WaveSpeed)
        Gc   = np.convolve(self.G, self.SourceFunc * self.dt)
        Gc   = Gc[0:self.TimeSteps]
        lim  = Gc.max() # get limit from maximum amplitude of function for plotting later on

        return Gc, lim
