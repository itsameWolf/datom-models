import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.animation as animation
from MonadModels import Monad, MonadOnMonad
from math import pi

anim_fig = plt.figure()
anim_plot = anim_fig.add_subplot()
steps = 400

monad1 = Monad()
monad2 = Monad()

monad2.Translate([-3*monad2.face_length,0])
monad2.Rotate(-pi/2, monad1.latch1.point1)

monad1_outline = ptc.Polygon(monad1.Getjoints(),color='b')
monad2_outline = ptc.Polygon(monad2.Getjoints(),color='r')

anim_plot.axis('equal')
anim_plot.set_xlim(-0.25, 0.15)
anim_plot.set_ylim(-0.1, 0.25)

anim_plot.add_patch(monad1_outline)
anim_plot.add_patch(monad2_outline)


def init():
    return monad1_outline, monad2_outline

def monadAnimation(i):

    val = abs(200-i)/200
    MonadOnMonad(monad1,monad2,val)

    monad1_outline.set_xy(monad1.Getjoints())
    monad2_outline.set_xy(monad2.Getjoints())

    return monad1_outline, monad2_outline


anim = animation.FuncAnimation(anim_fig, monadAnimation, range(steps), init_func=init)
writergif = animation.PillowWriter(fps=48)
anim.save("monad_on_monad.gif", writer=writergif)