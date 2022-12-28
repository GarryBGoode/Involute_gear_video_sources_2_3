from manim import *
class upd(Scene):
    def construct(self):
        self.totaltime = 0

        def SceneUpdater(dt):
            self.totaltime += dt

        self.add_updater(SceneUpdater)

        mobj = DecimalNumber(0)

        def mobjUpdater(mobj):
            mobj.increment_value(0.01)

        mobj.add_updater(mobjUpdater)

        self.add(mobj)

        self.wait(2)