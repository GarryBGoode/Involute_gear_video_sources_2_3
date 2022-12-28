from manim import *
from manim_cad_drawing_utils import *
from manim.mobject.geometry.tips import *
from manim_gearbox import *
from scipy.optimize import fsolve

class Profile_shift_1(MovingCameraScene):
    def construct(self):
        vt = ValueTracker(0)
        rt = ValueTracker(0)
        m=0.75
        z=8
        P_s = 0.7
        def Gear_constructor(params):
            gear0 = Gear(z,
                         module=m,
                         h_a=1, h_f=1.2,
                         profile_shift=vt.get_value(),
                         fill_color=WHITE,
                         fill_opacity=1,
                         stroke_opacity=0,
                         gen_params=params,
                         make_smooth=False)
            gear0.rotate((1/2+rt.get_value())* gear0.pitch_angle).shift(UP*gear0.rp)
            return gear0
        gear1 = Gear_constructor({})

        gear1.add_updater(lambda mob: mob.match_points(Gear_constructor(gear1.gen_params)))
        rack1 = Rack(5,m,h_a=1.2,h_f=1,fill_color=BLUE_A,fill_opacity=1,stroke_opacity=0).rotate(PI/2)

        Circle_Rb = Circle(radius=gear1.rb,
                           arc_center=gear1.get_center(),
                           color=TEAL,
                           fill_opacity=0)

        # self.add(gear1,rack1,Circle_Rb)
        self.play(FadeIn(gear1),FadeIn(rack1))
        self.play(vt.animate.set_value(P_s),rack1.animate.shift(DOWN*P_s*m))
        self.wait()
        self.play(rt.animate.set_value(0.5), rack1.animate.shift(RIGHT * 0.5 * rack1.pitch))
        self.wait()
        self.play(vt.animate.set_value(0),rack1.animate.shift(-P_s*DOWN*m))
        self.wait()
        # self.play(rt.animate.set_value(0), rack1.animate.shift(RIGHT * (-0.5) * rack1.pitch))

        tmax_0 = gear1.gen_params['tmax_invo']

        def Invo_ref_construc(p, gear_in=None):
            Invo_curve_ref = VMobject()
            if gear_in is None:
                gear0=Gear_constructor(gear1.gen_params)
            else:
                gear0=gear_in
            Invo_curve_ref.points=involute_point_gen(np.linspace(gear0.gen_params['tmin_invo']*0,
                                                                 gear0.gen_params['tmax_invo'],
                                                                 # tmax_0,
                                                                 10),
                                                     gear0.rb)
            Invo_curve_ref.rotate(-gear0.angle_ofs-gear0.pitch_angle/4,about_point=ORIGIN)
            Invo_curve_ref.shift(gear0.get_center())

            Invo_curve_ref.set_color(ORANGE)
            Invo_curve_ref.set_stroke(opacity=0)
            if p:
                Invo_curve_ref.flip(RIGHT,about_point=gear0.get_center())
            Invo_curve_ref.rotate(gear0.get_angle(),about_point=gear0.get_center())
            return Invo_curve_ref

        def Invo_grp_contstruct_right():
            Invo_curve_ref=Invo_ref_construc(0)
            Invo_grp = VGroup()
            Line_grp = VGroup()
            Invo_collect = VDict()
            Invo_collect['invo_curves'] = Invo_grp
            Invo_collect['lines'] = Line_grp
            for k in range(z):
                Invo_grp.add(Invo_curve_ref.copy())
                Line_grp.add(Line(start=gear1.get_center(),
                                  end=Invo_curve_ref.points[0, :]).match_style(Invo_curve_ref))
                Invo_collect.rotate(gear1.pitch_angle,about_point=gear1.get_center())
                Invo_collect.set_color(PURPLE_B)
                Invo_collect.set_stroke(opacity=1,width=5)
            return Invo_collect

        def Invo_grp_contstruct_left():
            Invo_curve_ref = Invo_ref_construc(1)
            Invo_grp = VGroup()
            Line_grp = VGroup()
            Invo_collect = VDict()
            Invo_collect['invo_curves']=Invo_grp
            Invo_collect['lines'] = Line_grp
            for k in range(z):
                Invo_grp.add(Invo_curve_ref.copy())
                Line_grp.add(Line(start=gear1.get_center(),end=Invo_curve_ref.points[0,:]).match_style(Invo_curve_ref))
                Invo_collect.rotate(gear1.pitch_angle, about_point=gear1.get_center())
                Invo_collect.set_stroke(opacity=1,width=5)
            return Invo_collect

        Invo_grp_right = VDict()
        Invo_grp_left = VDict()
        Invo_grp_right.set_color(PURPLE)
        Invo_grp_left.add_updater(lambda mob: mob.become(Invo_grp_contstruct_left()))
        Invo_grp_right.add_updater(lambda mob: mob.become(Invo_grp_contstruct_right()))
        Invo_grp_left.update()
        Invo_grp_right.update()
        # self.add(Invo_curve_ref)
        self.play(Create(Circle_Rb))

        self.play(vt.animate.set_value(P_s), rack1.animate.shift(DOWN * P_s * m))
        self.wait()
        self.play(vt.animate.set_value(0), rack1.animate.shift(-P_s * DOWN * m))
        self.wait()
        Invo_grp_right_2 = Invo_grp_right.copy()
        Invo_grp_right_2.clear_updaters()

        Invo_grp_left_2 = Invo_grp_left.copy()
        Invo_grp_left_2.clear_updaters()

        loc_gear_1 = Gear(z,
                          module=m,
                          h_a=1, h_f=1.2,
                          profile_shift=P_s,
                          make_smooth=False)
        loc_gear_2 = Gear(z,
                          module=m,
                          h_a=1, h_f=1.2,
                          profile_shift=0,
                          make_smooth=False)
        diff_angle = loc_gear_2.angle_ofs-loc_gear_1.angle_ofs

        self.play(Create(Invo_grp_left), Create(Invo_grp_right))
        self.add(Invo_grp_right_2)

        self.play(vt.animate.set_value(P_s),rt.animate.set_value(0.5+diff_angle/gear1.pitch_angle),
                  rack1.animate.shift((DOWN * P_s * m)+ LEFT*(P_s * m * np.tan(gear1.alpha*DEGREES))),run_time=3)

        self.wait()
        self.play(Rotate(Invo_grp_right_2, diff_angle*2,about_point=gear1.get_center()))
        self.wait()
        self.play(Rotate(Invo_grp_right_2, -diff_angle*2, about_point=gear1.get_center()))
        self.wait()





        self.play(FadeOut(rack1),FadeOut(Circle_Rb), FadeOut(Invo_grp_right_2))

        gear2=Gear(30,
                   module=m,
                   h_a=1, h_f=1.2,
                   profile_shift=0,
                   fill_color=WHITE,
                   fill_opacity=1,
                   stroke_opacity=0,
                   make_smooth=False).shift(gear1.rp*UP)
        gear0 = Gear_constructor({})
        gear2.shift(gear2.rp*DOWN)
        gear2.mesh_to(gear0,offset=0.1)

        Invo_ref_3 = Invo_ref_construc(1,gear2)
        Invo_ref_3.set_stroke(opacity=1)
        Invo_ref_3.set_color(ORANGE).rotate(gear2.pitch_angle*8,about_point=gear2.get_center())

        Invo_grp_3 = VGroup()
        for k in range(30):
            Invo_grp_3.add(Invo_ref_3.copy())
            Invo_grp_3.rotate(gear2.pitch_angle,about_point=gear2.get_center())


        # self.camera.frame.move_to(gear2.get_center()).scale(5)

        self.play(FadeIn(gear2))
        self.play(Create(Invo_grp_3,lag_ratio=0))
        gear2_grp = VGroup(gear2,Invo_grp_3)

        self.play(rt.animate.set_value(0.5+diff_angle/gear1.pitch_angle+1),Rotate(gear2_grp,-gear2.pitch_angle,about_point=gear2.get_center()),run_time=2)
        self.wait()
        self.play(rt.animate.set_value(0.5+diff_angle/gear1.pitch_angle),Rotate(gear2_grp,gear2.pitch_angle,about_point=gear2.get_center()),run_time=2)
        self.wait()

        self.play(FadeOut(gear2_grp))
        self.play(FadeIn(rack1),
                  FadeIn(Circle_Rb))

        self.play(vt.animate.set_value(0), rt.animate.set_value(0.5 + 0 * diff_angle / gear1.pitch_angle),
                  rack1.animate.shift(-(DOWN * P_s * m) + LEFT * (-P_s * m * np.tan(gear1.alpha * DEGREES))))
        Pitch_line = DashDot_mobject(Line(start=gear1.get_center()+gear1.rp*DOWN+LEFT*10,
                                          end=gear1.get_center()+gear1.rp*DOWN+RIGHT*10),
                                     num_dashes=30,
                                     dashed_ratio=0.7,
                                     dot_scale=0.5)

        Pitch_line.set_color(RED)
        Pitch_line_shift = Pitch_line.copy().shift(DOWN*P_s*m)
        Pitch_line_shift.add_updater(lambda mob: mob.move_to(Pitch_line).shift(DOWN*vt.get_value()*m))





        Constr_triangle = VMobject()
        Constr_triangle.set_points_as_corners([gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT,
                                               gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT+DOWN*P_s*m+LEFT*P_s*m*np.tan(20*DEGREES),
                                               gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT+LEFT*P_s*m*np.tan(20*DEGREES)*2,
                                               gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT])
        Constr_triangle.set_stroke(color=RED)

        Constr_line_1 = Line(
                             start=gear1.get_center() + gear1.rp * DOWN + gear1.pitch / 4 * LEFT + DOWN * P_s * m + LEFT * P_s * m * np.tan(
                                 20 * DEGREES),
                             end=gear1.get_center() + gear1.rp * DOWN + gear1.pitch / 4 * LEFT + LEFT * P_s * m * np.tan(20 * DEGREES),
                             color=RED,
                             stroke_width=2)
        Invo_grp_right_2.suspend_updating()
        Invo_grp_right.suspend_updating()
        Invo_grp_left.suspend_updating()
        self.play(FadeIn(Pitch_line),
                  FadeIn(Pitch_line_shift),
                  FadeOut(Invo_grp_right_2),
                  FadeOut(Invo_grp_right),
                  FadeOut(Invo_grp_left))
        # self.add(Constr_triangle,Constr_line_1)



        Constr_triangle.set_stroke(width=2)

        # alpha_txt = MathTex(r'\alpha').scale(1/9)
        # alpha_txt.move_to(gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT+DOWN*(P_s*m-0.2)+LEFT*(P_s*m*np.tan(20*DEGREES)-0.04))
        # self.add(alpha_txt)

        h_dim0 = Linear_Dimension(start=Constr_line_1.start, end=Constr_line_1.end, tip_len=0.1,
                                 offset=0.3,
                                 text=MathTex(r'X', color=RED),
                                 stroke_width=2,
                                 color=RED)

        h_dim = Linear_Dimension(start=Constr_line_1.start,end=Constr_line_1.end,tip_len=0.1,
                                 offset=0.3,
                                 text=MathTex(r'xm',color=RED),
                                 stroke_width=2,
                                 color=RED)

        self.play(vt.animate.set_value(P_s), rt.animate.set_value(0.5 + diff_angle / gear1.pitch_angle),
                  rack1.animate.shift((DOWN * P_s * m) + LEFT * (P_s * m * np.tan(gear1.alpha * DEGREES))))

        self.play(Create(h_dim0))
        self.wait()

        self.camera.frame.save_state()
        self.play(self.camera.frame.animate.move_to(Constr_triangle).scale(4 * P_s * m / 9))
        self.wait()
        self.play(ReplacementTransform(h_dim0,h_dim))
        self.play(Create(Constr_triangle))

        ds_dim = Linear_Dimension(start=gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT,
                                 end=gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT+LEFT*P_s*m*np.tan(20*DEGREES)*2,
                                 tip_len=0.1,
                                 offset=-0.1,
                                 text=MathTex(r'\Delta s',color=RED),
                                 stroke_width=2,
                                 color=RED)

        alpha_text = MathTex(r'2\alpha',color=RED).move_to(ds_dim['text'])
        alpha_text.scale(ds_dim['text'].height/alpha_text.height)
        alpha_text.shift(DOWN*P_s*0.6)
        alpha_arc=ArcBetweenPoints(gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT,
                                   gear1.get_center() + gear1.rp * DOWN + gear1.pitch / 4 * LEFT + LEFT * P_s * m * np.tan(
                                       20 * DEGREES) * 2,
                                   arc_center=gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT+DOWN*P_s*m+LEFT*P_s*m*np.tan(20*DEGREES),
                                   color=RED)
        alpha_arc.scale(0.35,about_point=gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT+DOWN*P_s*m+LEFT*P_s*m*np.tan(20*DEGREES))


        self.play(Create(ds_dim),Create(alpha_text),Create(alpha_arc))


        self.wait()

        dS_eq = MathTex(r'\Delta s=2xm \tan(\alpha)',color=RED)

        dS_eq.move_to(ds_dim['text'])
        dS_eq.shift((DOWN*3+RIGHT))
        # phi_eq = MathTex(r'\Delta\phi=\frac{2xm \tan(\alpha)}{r_p}',color=RED)
        # phi_eq.move_to(dS_eq.get_critical_point(DOWN),aligned_edge=UP)
        self.play(self.camera.frame.animate.restore())
        self.play(Create(dS_eq))

        self.play(FadeOut(h_dim),
                  FadeOut(ds_dim),
                  FadeOut(alpha_text),
                  FadeOut(alpha_arc))
        self.wait()

        s_dim = Linear_Dimension(start=gear1.get_center() + gear1.rp * DOWN - gear1.pitch / 4 * LEFT,
                                  end=gear1.get_center() + gear1.rp * DOWN + gear1.pitch / 4 * LEFT + LEFT * P_s * m * np.tan(
                                      20 * DEGREES) * 2,
                                  tip_len=0.1,
                                  offset=2,
                                  text=MathTex(r's', color=RED),
                                  stroke_width=2,
                                  color=RED)
        s_line = Line(start=gear1.get_center() + gear1.rp * DOWN - gear1.pitch / 4 * LEFT,
                      end=gear1.get_center() + gear1.rp * DOWN + gear1.pitch / 4 * LEFT + LEFT * P_s * m * np.tan(
                                      20 * DEGREES) * 2,
                      color=RED)
        S_eq = MathTex(r' s=m \left( \frac{\pi}{2} + 2x\tan(\alpha)\right)', color=RED)
        S_eq.move_to(dS_eq)

        self.play(Create(s_dim),Create(s_line),Transform(dS_eq,S_eq))

        self.wait()
        self.play(FadeOut(Constr_triangle))

        s_arc = Arc(radius=gear1.rp,
                    arc_center=gear1.get_center(),
                    start_angle=-PI/2+gear1.pitch_angle/4,
                    angle=-(gear1.pitch_angle/2+P_s*m*np.tan(20*DEGREES)*2/gear1.rp),
                    color=RED)

        self.wait()
        self.play(self.camera.frame.animate.move_to(s_line).scale(0.3))
        self.play(Transform(s_line,s_arc))
        self.wait()

        # Line_angle_1=Line(start=gear1.get_center(),
        #                   end = gear1.get_center()+gear1.rp*DOWN+gear1.pitch/4*LEFT,
        #                   color=RED)
        # Line_angle_2 = Line_angle_1.copy()
        # self.play(FadeOut(dS_eq))
        # self.play(Create(VGroup(Line_angle_1, Line_angle_2)))
        # self.play(Rotate(Line_angle_2,(-P_s*m*np.tan(20*DEGREES)*2/gear1.rp),about_point=gear1.get_center()))
        # dphi = MathTex(r'\Delta\phi',color=RED).move_to(phi_eq)
        # self.play(Create(dphi))
        # self.wait()

class Profile_shift_calc(MovingCameraScene):
    def construct(self):
        eq_sin2a_1 = MathTex(r"h=r_p \sin^2(\alpha)=r_p-r_d = 1.25m-xm")
        eq_sin2a_2 = MathTex(r"\left(\frac{z}{2} m\right) \sin^2(\alpha)=1.25m-xm")
        eq_sin2a_3 = MathTex(r"x=1.25-\frac{z}{2} \sin^2(\alpha)")
        eq_sin2a_4 = MathTex(r"z_{min}=\frac{2.5}{\sin^2(\alpha)}")
        eq_sin2a_5 = MathTex(r"x=1.25 \left(1-\frac{z}{z_{min}}\right)")

        eq_grp_2 = VGroup(eq_sin2a_1, eq_sin2a_2, eq_sin2a_3,eq_sin2a_4,eq_sin2a_5)
        eq_grp_2.arrange(DOWN)
        for k in range(5):
            self.play(Create(eq_grp_2[k]))
            self.wait()

class Ucut_Pshift_compare(MovingCameraScene):
    def construct(self):
        m0 = 0.08
        teeth_range = range(6,15)
        Undercut_grp = VGroup()
        Pshift_grp = VGroup()
        all_grp = VGroup()
        for k in teeth_range:
            gear1 = Gear(k,
                         module=m0*15/k,
                         h_f=1.25,
                         stroke_opacity=0,
                         fill_opacity=1,
                         fill_color=BLUE_A)

            gear2=Gear(k,
                       module=m0 * 15 / k,
                       stroke_opacity=0,
                       profile_shift=1.25*(1-k/20.5),
                       h_f=1.25,
                       fill_opacity=1,
                       fill_color=BLUE_C)
            num1 = Text(str(k),color=BLACK).scale(0.7).move_to(gear1.get_center())
            num2 = num1.copy().move_to(gear2.get_center())
            # ngear1=VGroup(gear1,num1)
            # ngear2 = VGroup(gear2, num2)
            gear1.add(num1)
            gear2.add(num2)
            loc_grp = VGroup(gear1,gear2).arrange()
            all_grp.add(loc_grp)

        # Undercut_grp.arrange_in_grid(cols=3,rows=3)
        # Pshift_grp.arrange_in_grid(cols=3,rows=3)
        # all_grp = VGroup(Undercut_grp,Pshift_grp)
        all_grp.arrange_in_grid(cols=3,rows=3)
        animations = [FadeIn(mob,lag_ratio=0) for grp in all_grp for mob in grp]
        self.play(AnimationGroup(*animations,lag_ratio=1),run_time=9)
        self.wait()


class Calc_tolerance(MovingCameraScene):
    def construct(self):
        eq_12 = MathTex(r"h_d = 1.25m")
        eq_sin2a_12 = MathTex(r"z_{min}=\frac{2.5}{\sin^2(\alpha)} = 21.3715")
        eq_sin2a_12_2 = MathTex(r" z_{min} \approx 22")
        eq_sin2a_10 = MathTex(r"z_{min}=\frac{2}{\sin^2(\alpha)} = 17.0972")
        eq_sin2a_10_2 = MathTex(r" z_{min} \approx 17")

        eq_10 = MathTex(r"h_d = 1m")
        eq_x_12 = MathTex(r"x=1.25 (1-\frac{z}{z_{min}})")
        eq_x_12_2 = MathTex(r"x \approx 1.25\left(1-\frac{z}{22}\right)")
        eq_x_10 = MathTex(r"x=1-\frac{z}{z_{min}}")
        eq_x_10_2 = MathTex(r"x\approx 1-\frac{z}{17}")

        eq_grp_left = VGroup(Text('with clearance'),
                             eq_12,
                             eq_sin2a_12,
                             eq_sin2a_12_2,
                             eq_x_12,
                             eq_x_12_2).arrange(DOWN)
        eq_grp_right = VGroup(Text('without clearance'),
                              eq_10,
                              eq_sin2a_10,
                              eq_sin2a_10_2,
                              eq_x_10,
                              eq_x_10_2).arrange(DOWN)

        grp_all = VGroup(eq_grp_left,eq_grp_right).arrange(RIGHT,buff=2)

        self.camera.frame.scale(grp_all.width/self.camera.frame.width*1.2)

        vertline = DashedVMobject(Line(start=-10*UP,end=10*UP,color=RED))

        # self.add(grp_all)
        self.play(Create(vertline))
        self.wait()
        self.play(Create(eq_grp_left,lag_ratio=2),run_time=4)
        self.wait()
        self.play(Create(eq_grp_right,lag_ratio=2),run_time=4)
        self.wait()

        self.play(FadeOut(eq_grp_right),FadeOut(eq_grp_left))

        rack1 = Rack(2,module=1,h_a=1.25,h_f=1).rotate(PI/2)
        rack2 = Round_Corners(Rack(2,module=1,h_a=1.25,h_f=1).rotate(PI/2),0.3)

        rack1.shift(LEFT*4)
        rack2.shift(RIGHT*4)
        self.play(Create(rack1),Create(rack2))
        self.wait()
        rack_line_1 = Line(start=rack1.get_center()+3*LEFT,
                           end=rack1.get_center()+3*RIGHT,
                           color=RED)

        rack_line_2 = Line(start=rack2.get_center() + 3 * LEFT,
                           end=rack2.get_center() + 3 * RIGHT,
                           color=RED)

        aux_line_1 = Line(start=rack1.get_center() + rack1.m * UP * 1.249 + LEFT * 12,
                          end=rack1.get_center() + rack1.m * UP * 1.249 + RIGHT * 12)
        aux_line_2 = Line(start=rack1.get_center() + rack1.m * UP + LEFT * 12,
                          end=rack1.get_center() + rack1.m * UP + RIGHT * 12)
        itx_points_0 = curve_intersection(rack_line_1, rack1)
        itx_points_02 = curve_intersection(rack_line_2, rack2)
        itx_points_1 = curve_intersection(aux_line_1,rack1)
        itx_points_2 = curve_intersection(aux_line_1, rack2)
        itx_points_3 = curve_intersection(aux_line_2, rack2)

        dim1 = Linear_Dimension(start=rack_line_1.get_nth_curve_function(0)(itx_points_0[0][0]),
                                end=aux_line_1.get_nth_curve_function(0)(itx_points_1[0][0]),
                                direction=LEFT,
                                offset=0.8,
                                outside_arrow=True,
                                text='1.25 m',
                                color=RED)

        dim2 = Linear_Dimension(start=rack_line_2.get_nth_curve_function(0)(itx_points_0[0][-1]),
                                end=aux_line_2.get_nth_curve_function(0)(itx_points_3[0][-1]),
                                direction=LEFT,
                                offset=-0.8,
                                outside_arrow=True,
                                text=' 1 m ',
                                color=RED)

        self.play(Create(DashedVMobject(rack_line_1)))
        self.play(Create(DashedVMobject(rack_line_2)))
        self.play(Create(dim1))
        self.play(Create(dim2))
        self.wait()
#
# with tempconfig({"quality": "medium_quality", "disable_caching": True, "save_last_frame":True}):
#     scene = Calc_tolerance()
#     scene.render()