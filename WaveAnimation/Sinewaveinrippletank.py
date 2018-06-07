import rippleTank as rt
import matplotlib.pyplot as plt

# creates a ripple tank
tank = rt.RippleTank()

# creates a source on the ripple tank
rt.Source(tank, rt.sineSource, xcorners = (-15, -14.5), ycorners = (-14, 14), freq = 5.0)


tank.simulateTime(3.0, animation_speed=0.5, fps=60)

ani = tank.makeAnimation()

plt.show()
