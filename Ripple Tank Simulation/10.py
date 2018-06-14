import sys
import numpy as np
from matplotlib.cm import jet
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
if sys.version_info[0] >= 3:
    from functools import reduce
from pandas import Series
import random
class Source():
    """
    Sources create perturbations on the ripple tank.
    """
    def __init__(self, rippletank, function, xcorners = (-0.5, 0.5), ycorners=(-0.5, 0.5), freq = 1, phase = 0, amplitude = 0.005):

        self.rippletank = rippletank #: parent tank
        self.X_grid = rippletank.X #: x grid coordinates
        self.Y_grid = rippletank.Y #: y grid coordinates
        self.xcorners = xcorners #: xcorners of the source
        self.ycorners = ycorners #: ycorners of the source

        self.freq = freq #: frequency of the source
        self.period = 1.0/freq #: period of the source
        self.phase = phase #: phase of the source

        self.function = function #: function that describes source behavior
        self.positions = getPositions(self.X_grid, self.Y_grid, self.xcorners, self.ycorners) #: position of the source

        if amplitude < 0 or amplitude > 1:
            raise(Exception("Amplitude is not valid"))
        self.amplitude = amplitude #: relative amplitude of the source

        self.rippletank.addSource(self)

    def evaluate(self, i):
        """
        Receives an int number related with an iterator, evaluates `function` using that number.

        Returns:
            np.ndarray: source values.
        """
        answer = self.function(self, i)
        if type(answer) == type(None):
            return np.zeros_like(self.X_grid)
        return answer*self.amplitude*self.rippletank.deep

class Mask():
    """
    Masks represent obstacules on the ripple tank. Depending on the `rel_deep`
    they can be seen as barriers (`rel_deep=0`) or simply objects that vary the deep of the tank.
    """
    def __init__(self, rippletank, rel_deep = None):
        self.rippletank = rippletank #: parent Tank object
        self.rippletank.addMask(self)

        self.rel_deep = rel_deep #: relative deep of the mask
        self.mask = np.ones_like(self.rippletank.X) #: mask array

    def fromFunc(self, func, args = (), kwargs = {}):
        """
        Mask array can be made by calling a function `func`, the first argument of the
        function needs to be the mask object, other arguments are send to the function
        with the `args` parameter, and the keyword arguments with `kwargs`.
        """
        self.mask = func(self, *args, **kwargs)
        if type(self.rel_deep) != type(None):
            self.mask[self.mask == 0] = self.rel_deep

        self.applyMask()

    def fromArray(self, array, rel_deep = None):
        """
        Mask array can assigned with and incoming 2d array.
        The array must have the same dimensions of the ripple tank, contain ones
        where no effect is wanted and zeros where the deep wants to be modified.
        If the mask object has a `rel_deep` or a value is set by parameter,
        zeros will be changed to `rel_deep`.

        Raises:
            Exception: "Input array does not have the ripple tank dimensions."
            Exception: "Deep must be between 0 and 1."
        """
        if array.shape != self.rippletank.X.shape:
            raise(Exception("Input array does not have the ripple tank dimensions."))
        self.mask = array
        if type(rel_deep) != type(None):
            self.rel_deep = rel_deep
        if type(self.rel_deep) != type(None):
            self.mask[self.mask == 0] = self.rel_deep
        if self.mask.max() > 1 or self.mask.min() < 0:
            raise(Exception("Deep must be between 0 and 1."))
        self.applyMask()

    def applyMask(self):
        """
        Applies the mask object to the `rippletank`.
        """
        self.rippletank.applyMultipleMasks()

    def __sum__(self, other):
        return self.mask + other

    def __sub__(self, other):
        return self.mask - other

    def __mul__(self, other):
        return self.mask * other

    def __rsum__(self, other):
        return self.mask + other

    def __rsub__(self, other):
        return self.mask - other

    def __rmul__(self, other):
        return self.mask * other
    
def getPositions(X_grid, Y_grid, xcorners, ycorners):
    """
    Using the coordinates in `X_grid` and `Y_grid` returns a 2d boolean array
    of a rectangle with `xcorners` and `ycorners`.

    Returns:
        np.ndarray: 2d boolean array.
    """
    xconditions = (X_grid >= min(xcorners)) & (X_grid <= max(xcorners))
    yconditions = (Y_grid >= min(ycorners)) & (Y_grid <= max(ycorners))
    return xconditions*yconditions

def sineSource(source, i):
    """
    Sine function.

    Returns:
        np.ndarray: array with sine values on source positions.
    """
    t = source.rippletank.dt*i
    answer = np.zeros_like(source.X_grid)
    value = np.sin(2*np.pi*source.freq*t + source.phase)
    answer[source.positions] = value
    return answer

def sineSource2(source, i):
    """
    Sine function Plus White Noise.

    Returns:
        np.ndarray: array with sine values on source positions.
    """
    t = source.rippletank.dt*i
    answer = np.zeros_like(source.X_grid)
    value = np.sin(2*np.pi*source.freq*t + source.phase)
    answer[source.positions] = value + [10*random.gauss(0.0, np.sin(i)) for i in range(100)]
    return answer

class RippleTank():
    """
    RippleTank objects are the core of the simulation. They contain both sources and masks.
    The class describes the space in which the waves will move and the force acting on them.
    """
    def __init__(self, xdim = (-15, 15), ydim = (-15, 15), deep = 1.0,
                n_cells_x = 100, n_cells_y = 100, mask = 1.0,
                bc = 'open', alpha = 0.45, units = 'cm'):
        posible_bcs = 'open', 'close'
        if not bc in posible_bcs:
            raise(Exception("'%s' is not a valid boundary condition."%bc))

        posible_units = 'cm', 'm'
        if not units in posible_units:
            raise(Exception("'%s' are not a valid units."%units))

        x = np.linspace(xdim[0], xdim[1], n_cells_x)
        y = np.linspace(ydim[0], ydim[1], n_cells_y)
        self.dx = x[1] - x[0]
        self.dy = y[1] - y[0]

        self.X, self.Y = np.meshgrid(x, y) #: two 2d arrays describing the coordinates of the tank

        self.xdim = xdim #: x dimensions of the grid
        self.ydim = ydim #: y dimensions of the grid

        self.n_cells_x = n_cells_x #: number of cells on x
        self.n_cells_y = n_cells_y #: number of cells on y
        self.units = units #: units used
        self.bc = bc #: boundary conditions, 'open' or 'close'

        self.mask = mask #: mask appplied to the tank
        if self.mask == 1:
            self.mask = np.ones_like(self.X)

        self.deep = deep #: deep of the tank, must be positive
        if self.deep < 0:
            raise(Exception('deep value must be positive.'))
        self.masked_deep = deep*self.mask #: deep on every point

        self.g = 9.8 #: gravity value
        if self.units == 'cm':
            self.g = 980
        self.speed = np.sqrt(self.g*deep) #: speed of propagation on each point

        self.dt = alpha*min(self.dx, self.dy)/self.speed #: dt value

        self.ratiox = (self.speed*self.dt/self.dx)**2 #: finite differences quotient on x
        self.ratioy = (self.speed*self.dt/self.dy)**2 #: finite differences quotient on y

        self.fig = None #: matplotlib figure
        self.ax = None #: matplotlib axes
        self.sources = [] #: stores sources
        self.masks = [] #: stores masks
        self.amplitude = None #: wave amplitude values
        self.complete_values = None #: wave amplitude + deep
        self.forbidden_pos = None #: positions where sources stand

        self.sim_duration = None #: time to simulate
        self.animation_speed = 1.0 #: relative reproduction speed
        self.fps = 24.0 #: frames per second value

        self.extent = [self.xdim[0], self.xdim[1], self.ydim[0], self.ydim[1]] #: matplotlib extent parameter

    def calcSpeed(self, values):
        """
        Calculates the propagation speed depending on the actual height of the wave.
        """
        deep = values + self.masked_deep
        speed = np.sqrt(self.g*deep)
        self.speed = speed

    def solveBorders(self, i):
        """
        Determins the state i+1 of the boundaries using the state i.
        """
        ratiox = self.speed*self.dt/self.dx
        ratioy = self.speed*self.dt/self.dy
        self.amplitude[i+1, 0] = ratioy[0]*(self.amplitude[i, 1] - self.amplitude[i, 0]) + self.amplitude[i, 0]
        self.amplitude[i+1, -1] = -ratioy[-1]*(self.amplitude[i, -1] - self.amplitude[i, -2]) + self.amplitude[i, -1]

        self.amplitude[i+1, :, 0] = ratiox[:, 0]*(self.amplitude[i, :, 1] - self.amplitude[i, :, 0]) + self.amplitude[i, :, 0]
        self.amplitude[i+1, :, -1] = -ratioy[:, -1]*(self.amplitude[i, :, -1] - self.amplitude[i, :, -2]) + self.amplitude[i, :, -1]

    def solveInstant(self, i):
        """
        Solve the differential equation for a single instant of time, from i to i+1.
        """
        self.calcSpeed(self.amplitude[i+1])
        self.amplitude[i+1, 1:-1, 1:-1] = 2*self.amplitude[i, 1:-1, 1:-1] - self.amplitude[i-1, 1:-1, 1:-1]\
                                    + self.getSecondPartEquation(i)[1:-1, 1:-1]
        if self.bc == 'open':
            self.solveBorders(i)

        pos = self.speed == 0
        self.amplitude[i+1, pos] = 0

    def addSource(self, source):
        """
        Includes a source to the ripple tank.
        """
        self.sources += [source]
        if source.period < self.dt:
            self.setdt(0.1*source.period)

    def evaluateSources(self, i):
        """
        Evaluates all sources in the tank.

        Returns:
            np.ndarray: 2d array with the values of the source at the i instant.
        """
        initial = np.zeros_like(self.X)
        if len(self.sources) == 0:
            return initial
        for source in self.sources:
            initial = initial + source.evaluate(i)
        return initial

    def applySources(self, i):
        """
        Sets the sources values in the amplitude array of the waves.
        """
        values = self.evaluateSources(i)
        positions = self.forbidden_pos
        self.amplitude[i, positions] = values[positions]

    def getSourcesPositions(self):
        """
        Gets all the positions of the sources as a boolean array.

        Returns:
            np.ndarray: 2d boolean array.
        """
        positions = np.zeros_like(self.X, dtype=bool)
        if len(self.sources) == 0:
            return positions
        for source in self.sources:
            positions = positions + source.positions
        return positions

    def setdt(self, dt):
        """
        Sets the delta t value.
        """
        self.dt = dt
        self.ratiox = (self.speed*self.dt/self.dx)**2
        self.ratioy = (self.speed*self.dt/self.dy)**2

    def getSecondPartEquation(self, i):
        """
        Evaluates the central differences on x and y for the instant i+1.

        Returns:
            np.ndarray: 2d amplitude values at the instant i+1.
        """
        temp = np.zeros_like(self.X)
        ratiox = (self.speed*self.dt/self.dx)**2
        ratioy = (self.speed*self.dt/self.dy)**2

        if isinstance(self.speed, np.ndarray):
            ratiox = ratiox[1:-1, 1:-1]
            ratioy = ratioy[1:-1, 1:-1]
        temp[1:-1, 1:-1] = ratiox*(self.amplitude[i, 1:-1, :-2] - 2*self.amplitude[i, 1:-1, 1:-1] + self.amplitude[i, 1:-1, 2:])\
                    + ratioy*(self.amplitude[i, :-2, 1:-1] - 2*self.amplitude[i, 1:-1, 1:-1] + self.amplitude[i, 2:, 1:-1])

        return temp

    def simulateTime(self, sim_duration, animation_speed=1.0, fps=24.0):
        """
        Simulates an interval of time, if the animation_speed with the current fps value
        does not match the sim_duration, modifies the `dt` value.

        Returns:
            np.ndarray: 3d array, extra dimension represents time.
        """
        self.sim_duration = sim_duration
        self.animation_speed = animation_speed
        self.fps = fps

        frames = round(fps*sim_duration/animation_speed)
        required_dt = sim_duration/frames
        if required_dt < self.dt:
            self.setdt(required_dt)

        points = round(self.sim_duration/self.dt)
        return self.solvePoints(int(points))

    def solvePoints(self, n_instants):
        """
        Simulates `n_instants` of time.

        Returns:
            np.ndarray: 3d array, extra dimension represents time.
        """
        self.amplitude = np.zeros((n_instants, self.n_cells_y, self.n_cells_x))
        self.amplitude[0] = self.evaluateSources(0)# + self.masked_deep
        self.amplitude[1] = self.amplitude[0] + self.getSecondPartEquation(1)# + self.masked_deep
        self.forbidden_pos = self.getSourcesPositions()

        for i in range(1, n_instants-1):
            self.solveInstant(i)
            self.applySources(i)
            self.applySources(i+1)

        self.complete_values = self.amplitude + self.masked_deep 
        
        return self.complete_values

    def applyMask(self, frame):
        """
        Applies a numpy mask.

        Returns:
            np.ma: masked array.
        """
        if isinstance(self.mask, np.ndarray):
            return np.ma.masked_where(self.mask == 0, frame)
        return frame

    def addMask(self, mask):
        """
        Adds a mask to the ripple tank.
        """
        if isinstance(mask, Mask):
            self.masks += [mask]
        else:
            raise(Exception('Mask type is not valid.'))

    def animate(self, i, values, skip):
        """
        Function used by matplotlib's FuncAnimation.

        Returns:
            matplotlib object: imshow.
            matplotlib object: text.
        """
        i = i*skip
        temp = self.applyMask(values[i])
        self.wave_show.set_array(temp)

        t = i*self.dt
        self.time_label.set_text("%.3f s"%t)
        return self.wave_show, self.time_label,

    def configPlot(self, figsize=(12, 9), xlabel = None, ylabel = None, cmap = jet,
                    vmin = None, vmax = None, cbar_label = None, origin='lower'):
        """
        Configures the plot.

        Returns:
            matplotlib.figure: figure containing the main plot.
            matplotlib.axes: axes containing the imshow.
        """
        self.fig, self.ax = plt.subplots(figsize = figsize)

        if type(cmap) is str:
            exec("from matplotlib.cm import %s as cmap"%cmap, globals())
            exec("cmap.set_bad('black', 1.0)", globals())
        else:
            cmap.set_bad('black', 1.0)

        if xlabel == None:
            xlabel = "$x$ (%s)"%self.units

        if ylabel == None:
            ylabel = "$y$ (%s)"%self.units

        if cbar_label == None:
            cbar_label = "Deep (%s)"%self.units

        if type(self.amplitude) != type(None):
            binary_mask = not ((self.masked_deep > 0) & (self.masked_deep < 1)).any()
            if vmin == None:
                if binary_mask:
                    vmin = self.amplitude.min() + self.deep
                else:
                    vmin = self.complete_values.min()
            if vmax == None:
                if binary_mask:
                    vmax = self.amplitude.max() + self.deep
                else:
                    vmax = self.complete_values.max()


        self.time_label = self.ax.text(self.xdim[0] + 0.1*(self.xdim[1] - self.xdim[0]),
                    self.ydim[1] - 0.1*(self.ydim[1] - self.ydim[0]), "")
        self.wave_show = self.ax.imshow(np.zeros_like(self.X), cmap = 'gray',
                        vmin = 0.95, vmax = 1.05, origin = origin, animated = True,
                        extent = self.extent)

        cbar = self.fig.colorbar(self.wave_show)
        cbar.set_label(cbar_label)

        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)

        return self.fig, self.ax

    def verifyData(self, data):
        """
        Verifies if parameter data is different from None. If no simulation has between
        run, simulates 100 instants of time.

        Returns:
            np.ndarray: 3d array, extra dimension represent time.
        """
        if type(data) == type(None):
            if type(self.complete_values) == type(None):
                return self.solvePoints(100)
            else:
                return self.complete_values
        return data

    def captureFrame(self, data=None, fig = None, frame=-1):
        """
        Uses the rippletank's to plot a single instant of time of `data`.

        Returns:
            matplotlib.figure: figure containing the main plot.
            matplotlib.axes: axes containing the imshow.
        """
        data = self.verifyData(data)[frame]
        if fig == None and self.fig == None:
            self.configPlot()

        self.wave_show.set_array(self.applyMask(data))
        return self.fig, self.ax

    def makeAnimation(self, data=None, fig = None, fps = None, duration = None):
        """
        Makes an animation of `data`, it only uses the required frames depending
        on the duration and fps value.

        Returns:
            matplotlib.animation.FuncAnimation: animation of the data.
        """
        data = self.verifyData(data)

        if fps == None:
            fps = self.fps
        else:
            self.fps = fps

        if duration == None:
            if self.sim_duration != None:
                duration = self.sim_duration
            else:
                duration = 10
        else:
            self.sim_duration = duration

        skip = self.animation_speed*data.shape[0]/(fps*duration)
        skip = round(skip)
        if skip == 0:
            skip = 1

        if fig == None and self.fig == None:
            self.configPlot()

        ani = FuncAnimation(self.fig, self.animate, frames = data.shape[0]//skip,
                interval=50, fargs=(data, skip), blit=True)

        return ani

tank = RippleTank()
Source(tank, sineSource2, xcorners = (-15, -15), ycorners = (-15, 15), freq = 5.0)
tank.simulateTime(10.0, animation_speed=1, fps=60)
ani = tank.makeAnimation()
plt.show()
