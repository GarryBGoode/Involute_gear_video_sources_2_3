from manim import *
from manim_gearbox import *
from scipy.optimize import fsolve
from manim_cad_drawing_utils import *


class SlideRoll(MovingCameraScene):
    def construct(self):
        gear1 = Gear(12,module=2,profile_shift=0.6, fill_color=BLUE_A,fill_opacity=1,stroke_opacity=0)
        gear2 = gear1.copy()

        gear1.shift(UP*(gear1.rp+gear1.X*gear1.m/2))
        gear2.shift(DOWN*gear1.rp)
        gear2.mesh_to(gear1)

        # gear_small_1=Gear(12,module=1.5,profile_shift=0.6,
        #                   cutout_teeth_num=8,
        #                   fill_color=BLUE_A,fill_opacity=1,stroke_opacity=0)

        def sin_func(s, a=0.1, f=672):
            return a*(np.cos(s*f*TAU)-1)

        gear_wave1 = Path_Offset_Mobject(gear1,sin_func,num_of_samples=5000, color=RED_E,stroke_width=3)
        gear_wave2 = Path_Offset_Mobject(gear2,sin_func,num_of_samples=5000, color=DARK_BLUE,stroke_width=3)

        geargrp1 = VGroup(gear1,gear_wave1)
        geargrp2 = VGroup(gear2, gear_wave2)

        self.add(geargrp1, geargrp2)
        self.play(Rotate(geargrp1,gear1.pitch_angle*3),
                  Rotate(geargrp2,-gear1.pitch_angle*3),
                  run_time=30,
                  rate_func=linear)
