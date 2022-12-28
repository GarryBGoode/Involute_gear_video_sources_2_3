from manim import *
from manim_gearbox import *
from scipy.optimize import fsolve
from manim_cad_drawing_utils import *


m = 0.5
z = 32
strw=2

class Intro(MovingCameraScene):
    def construct(self):
        rt = ValueTracker(0)
        gear1 = Gear(z, module=m, h_f=1.25, stroke_width=strw)
        gear2 = Gear(8,module=m, h_f=1.25, profile_shift=0.3)
        gear1.reverse_direction()
        # gear1.rotate(PI/2)
        gear2.set_stroke(opacity=0, family=False)
        gear2.set_fill(opacity=1,color=BLUE_B)

        gear3=gear2.copy()
        gear2.shift(RIGHT*(gear1.rp+gear2.rp))
        gear3.shift(LEFT * (gear1.rp + gear2.rp))
        gear2.mesh_to(gear1)
        def gear2_updater(mob:Gear):
            center = rotate_vector(RIGHT*(gear1.rp+gear2.rp),rt.get_value())
            mob.move_to(center)
            mob.mesh_to(gear1)

        def gear3_updater(mob: Gear):
            center = rotate_vector(LEFT * (gear1.rp + gear2.rp), rt.get_value())
            mob.move_to(center)
            mob.mesh_to(gear1)
        gear2.add_updater(gear2_updater)
        gear3.add_updater(gear3_updater)
        gpm = Path_mapper(gear1)

        gear_a = VMobject()
        def gear_anim(mob: Gear):
            mob.pointwise_become_partial(gear1,0,gpm.equalize_alpha(rt.get_value()/TAU))
        gear_a.add_updater(gear_anim)
        gear_b = VMobject()
        gear_b.add_updater(lambda mob: mob.pointwise_become_partial(gear1,0.5,0.5+gpm.equalize_alpha(rt.get_value()/TAU)))


        self.camera.frame.scale((gear2.ra*5 + gear1.ra*2) / 9)
        self.play(FadeIn(gear2),FadeIn(gear3))
        self.add(gear_a,gear_b)
        # point = rotate_vector(UP * (gear1.rp + gear2.rp * 0.5), rt.get_value())
        # self.play(self.camera.frame.animate.move_to(point))
        # self.play(self.camera.frame.animate.scale(gear2.ra * 2.5 / 9))


        self.wait()


        def gear_tracker(mob: Mobject):
            point = rotate_vector(UP * (gear1.rp + gear2.rp*0.5), rt.get_value())
            mob.move_to(point)

        # self.camera.frame.add_updater(gear_tracker)

        self.play(rt.animate.set_value(TAU/2),run_time=10)
        self.wait(2)
        self.play(self.camera.frame.animate.scale(m * 15 / self.camera.frame.width).move_to(gear1.rp * UP),
                  gear_a.animate.set_stroke(width=strw),
                  gear_b.animate.set_stroke(width=strw),
                  run_time=2)
        self.wait()



class Add_height(MovingCameraScene):
    def construct(self):
        # m = 0.5
        # z = 32
        at = ValueTracker(20)
        xt = ValueTracker(0)
        adden = ValueTracker(1)
        deden = ValueTracker(1.25)
        gear1 = Gear(z,module=m, alpha=at.get_value(), h_a=adden.get_value(), h_f=deden.get_value(), stroke_width=strw)
        gear1.add_updater(lambda mob: mob.match_points(Gear(z,
                                                            module=m,
                                                            alpha=at.get_value(),
                                                            h_a=adden.get_value(),
                                                            h_f=deden.get_value(),
                                                            profile_shift=xt.get_value())))

        # the rack has to have reverse params
        rack1 = Rack(6,
                     module=m,
                     alpha=at.get_value(),
                     h_f=adden.get_value(),
                     h_a=deden.get_value(),
                     stroke_width=1,
                     color=RED)
        def rack_updater(mob: Rack):
            mob.match_points(Rack(6,
                                  module=m,
                                  alpha=at.get_value(),
                                  h_f=adden.get_value(),
                                  h_a=deden.get_value(),
                                  stroke_width=1,
                                  make_closed=False))
            mob.rotate(-PI/2)
            mob.shift(UP*gear1.rp)
            mob.shift(UP*xt.get_value()*m)
        rack1.add_updater(rack_updater)

        def get_max_height(alpha):
            return np.pi*0.5/np.tan(alpha*np.pi/180)/2

        pitch_circ_base = Circle(radius=gear1.rp,stroke_width=1)
        pitch_circ = DashDot_mobject(pitch_circ_base,num_dashes=int(z*2), dashed_ratio=0.65, dash_offset=0.35/2)

        self.camera.frame.save_state()
        # self.camera.frame.scale(m * 15 / 16)
        # self.camera.frame.move_to(gear1.rp*UP)
        self.camera.frame.scale(m * 15 / self.camera.frame.width).move_to(gear1.rp * UP)
        self.add(gear1)
        self.play(Create(rack1))
        self.play(Create(pitch_circ))


        dima = Linear_Dimension((gear1.ra)*UP,gear1.rp*UP,
                                text="1m",
                                outside_arrow=True,
                                offset=m*PI/2 * 1,
                                stroke_width=1,
                                tip_len=0.1,
                                color=RED)
        dimb = Linear_Dimension((gear1.rf) * UP, gear1.rp * UP,
                                text="1.25m",
                                outside_arrow=True,
                                offset=m * PI / 2 * 0.8,
                                stroke_width=1,
                                tip_len=0.1,
                                color=RED)

        self.play(Create(dima))
        self.wait()
        self.play(Create(dimb))
        self.wait()
        self.play(Uncreate(dima), Uncreate(dimb))
        self.wait()
        self.play(adden.animate.set_value(get_max_height(at.get_value())),
                  deden.animate.set_value(get_max_height(at.get_value())), run_time=2)
        self.wait()
        hmax = gear1.hmax
        dim = Linear_Dimension((gear1.rb+hmax)*UP,gear1.rp*UP,
                               text=f"{(hmax+gear1.rb-gear1.rp)/m:.3}m",
                               outside_arrow=True,
                               offset=m*PI/2,
                               stroke_width=1,
                               tip_len=0.1,
                               color=RED)

        self.play(Create(dim))
        self.wait()
        self.play(Uncreate(dim))

        # dim_a = Angle_Dimension_3point(start=rack1.get_center()+rack1.h_f*UP+UP*0.1,
        #                                end=rack1.get_center()+rack1.h_f*UP+rotate_vector(UP,at.get_value()*DEGREES)*0.1,
        #                                arc_center=rack1.get_center()+rack1.h_f*UP,
        #                                offset=0.3,
        #                                color=RED,
        #                                outside_arrow=True,
        #                                stroke_width=2)

        # dim_a = MathTex(r'\alpha=',f'{at.get_value():.0f}').move_to(rack1.get_center()+rack1.h_f*UP+UP*0.5)
        dim_a = Text(f'α={at.get_value():.0f}°',color=RED).move_to(rack1.get_center() + rack1.h_f * UP + UP * 0.5).scale(0.5)
        dim_a.add_updater( lambda mob: mob.become(Text(f'α={at.get_value():.1f}°',color=RED).move_to(rack1.get_center() + rack1.h_f * UP + UP * 0.5).scale(0.5)))

        self.play(Create(dim_a))
        self.play(at.animate.set_value(14.5),run_time=2)
        self.wait()

        self.play(adden.animate.set_value(get_max_height(14.5)),
                  deden.animate.set_value(get_max_height(14.5)),
                  run_time=2)
        self.wait()
        self.play(FadeOut(dim_a))

        self.play(at.animate.set_value(20),
                  adden.animate.set_value(1),
                  deden.animate.set_value(1.25),
                  run_time=2)

        self.wait()
        self.play(xt.animate.set_value(0.5))
        self.wait()
        self.play(xt.animate.set_value(-0.5))
        self.wait()
        self.play(xt.animate.set_value(0))
        self.wait()

        self.play(Uncreate(rack1))
        self.play(self.camera.frame.animate.move_to(ORIGIN).scale(gear1.ra*2.5/self.camera.frame.height),
                  gear1.animate.set_stroke(width=4,family=True),
                  pitch_circ.animate.set_stroke(width=4,family=True)
                  )

        dimdp = Linear_Dimension(gear1.rp*UP,
                                 gear1.rp*DOWN,
                                 text="zm",
                                 color=RED,
                                 tip_len=0.5,
                                 offset=gear1.rp+2)

        dimda = Linear_Dimension(gear1.ra * UP,
                                 gear1.ra * DOWN,
                                 text="(z+2)m",
                                 color=RED,
                                 tip_len=0.5,
                                 offset=gear1.rp+4)

        self.play(Create(dimdp))
        self.play(Create(dimda))
        self.wait(2)

# with tempconfig({"quality": "medium_quality", "disable_caching": True, "save_last_frame":True}):
#     scene = Add_height()
#     scene.render()