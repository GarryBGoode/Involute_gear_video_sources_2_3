import numpy as np
from manim import *
from manim_gearbox import *
from scipy.optimize import fsolve
from manim_cad_drawing_utils import *


class Undercut(MovingCameraScene):
    def construct(self):
        m=1.5
        gear1 = Gear(8,m,h_f=1.25,stroke_opacity=0,fill_opacity=1,fill_color=WHITE,nppc=7,make_smooth=False)
        rz=9


        Invo_ref_pt = involute_point_gen(np.linspace(0,gear1.gen_params['tmax_invo'],20),gear1.rb)
        Invo_ref_mob = VMobject()
        Invo_ref_mob.set_points(Invo_ref_pt)
        Invo_ref_mob.rotate(-gear1.angle_ofs-gear1.pitch_angle/4,about_point=ORIGIN)
        Invo_ref_mob2 = Invo_ref_mob.copy().flip(RIGHT,about_point=ORIGIN)
        Invo_refs = VGroup(Invo_ref_mob,Invo_ref_mob2)

        gear1.shift(LEFT*gear1.rp)
        Invo_refs.shift(LEFT*gear1.rp)
        Invo_refs.set_stroke(color=RED)


        Mask_mob = Invo_ref_mob.copy()
        Mask_mob.reverse_direction()
        # Mask_mob.add_line_to(gear1.get_center())
        Mask_mob.add_line_to(Mask_mob.points[-1,:]+LEFT*gear1.rp)
        Mask_mob.add_line_to(Invo_ref_mob2.points[0, :]+LEFT*gear1.rp)
        Mask_mob.add_line_to(Invo_ref_mob2.points[0, :])
        Mask_mob.append_points(Invo_ref_mob2.points)
        Mask_mob.add_line_to(Invo_ref_mob.points[-1,:])
        Mask_mob.set_stroke(opacity=1,width=4)
        # Ucut_zone = Difference(Mask_mob,gear1)
        # Ucut_zone.set_fill(color=RED,opacity=1)
        # Ucut_zone.set_stroke(opacity=0)
        Mask_mob.set_fill(color=RED,opacity=1)

        ucut_points = involute_point_gen(np.linspace(gear1.gen_params['tmin_ucut'],gear1.gen_params['tmax_ucut'],30),
                                         gear1.rp,
                                         rad_offs=gear1.gen_params['rad_offs'],
                                         tan_offs=gear1.gen_params['tan_offs'])
        Undercut_line1=VMobject(stroke_color=RED)
        Undercut_line1.points=ucut_points
        Undercut_line1.rotate(gear1.alpha*DEGREES+gear1.pitch_angle/4+gear1.angle_ofs,
                              about_point=ORIGIN).shift(gear1.rp*LEFT)
        Undercut_line2 = Undercut_line1.copy().flip(RIGHT,about_point=ORIGIN)

        Invo_refs.add(Undercut_line1,Undercut_line2)
        Invo_refs.set_stroke(width=4)


        rack1 = Rack(rz,m,
                     h_a=1.25,
                     h_f=1,
                     stroke_opacity=0,fill_opacity=1,fill_color=BLUE_A).rotate(PI)
        rack1.shift(rack1.pitch/2*DOWN)



        # Invo_refs.add()

        self.add(self.camera.frame)
        # self.play(self.camera.frame.animate.move_to(rack1.get_center()))
        # self.camera.frame.add_updater(lambda mob: mob.move_to(rack1.get_center()))
        # self.camera.rotate(PI/4)

        Invo_refs.set_stroke(opacity=0)


        all_grp = VGroup(gear1, rack1, Invo_refs)
        all_grp.rotate(-PI / 2, about_point=ORIGIN)

        Ucut_grp = VGroup()
        for k in range(int(gear1.z)):
            Ucut_grp.add(Undercut_line1.copy())
            Ucut_grp.add(Undercut_line2.copy())
            Ucut_grp.rotate(gear1.pitch_angle,about_point=gear1.get_center())
        Ucut_grp.set_stroke(opacity=1,width=6)

        # self.add(gear1)
        self.play(FadeIn(gear1))
        self.play(FadeIn(rack1))
        self.play(FadeIn(Ucut_grp))
        self.add(Invo_refs)
        
        rot_grp = VGroup(Invo_refs,Ucut_grp)

        Ucut_pointer = Pointer_To_Mob(Ucut_grp[0], 0.2, 'Undercut', offset=(UP * 1 + LEFT * 0.3), color=RED)
        # Ucut_pointer.add_updater(
        #     lambda mob: mob.become(Pointer_To_Mob(Ucut_grp[1], 0, 'Undercut', offset=(UP * 0.5 + RIGHT * 0.2), color=RED)))
        self.play(Create(Ucut_pointer))
        self.wait(2)

if 0:

        self.play(Uncreate(Ucut_pointer))
        # rotangle = gear1.pitch_angle / 4 + gear1.angle_ofs + gear1.alpha*DEGREES + gear1.gen_params['tmin_ucut']
        # self.play(Rotate(gear1,rotangle),
        #           Rotate(rot_grp, rotangle,about_point=gear1.get_center()),
        #           rack1.animate.shift(RIGHT*rotangle*gear1.rp),run_time=2)

        # self.wait()
        # self.play(Undercut_line2.animate.set_stroke(opacity=0.5))
        rotangle = gear1.gen_params['tmax_ucut']-gear1.gen_params['tmin_ucut'] + gear1.pitch_angle / 4 + gear1.angle_ofs + gear1.alpha*DEGREES + gear1.gen_params['tmin_ucut']
        self.play(Rotate(gear1,rotangle),
                  Rotate(rot_grp, rotangle,about_point=gear1.get_center()),
                  rack1.animate.shift(RIGHT*rotangle*gear1.rp),run_time=12)



        self.wait()

        self.camera.frame.save_state()
        self.play(self.camera.frame.animate.scale(0.3).move_to(rack1.get_center()+rack1.m*UP),run_time=3)
        self.camera.frame.add_updater(lambda mob:
                                      mob.move_to(rack1.get_center()+rack1.m*UP))

        self.play(Ucut_grp.animate.set_stroke(opacity=0))
        rotangle = -(gear1.gen_params['tmax_ucut'] - gear1.gen_params['tmin_ucut'])
        self.play(Rotate(gear1, rotangle),
                  Rotate(rot_grp, rotangle, about_point=gear1.get_center()),
                  rack1.animate.shift(RIGHT * rotangle * gear1.rp), run_time=5)

        rack2 = Rack(rz,m,
                     h_a=1,
                     h_f=1.25,
                     stroke_opacity=0,fill_opacity=1,fill_color=BLUE_A)

        rack2.shift(rack1.get_center() - rack2.get_center())
        rack2.rotate(PI/2)

        rack1_save = rack1.copy()

        self.play(Transform(rack1,rack2))
        self.wait()

        rotangle = (gear1.gen_params['tmax_ucut'] - gear1.gen_params['tmin_ucut'])
        self.play(Rotate(gear1, rotangle),
                  Rotate(rot_grp, rotangle, about_point=gear1.get_center()),
                  rack1.animate.shift(RIGHT * rotangle * gear1.rp), run_time=5)
        rack1_save.shift(rack1.get_center()-rack1_save.get_center())
        self.play(Transform(rack1,rack1_save))
        self.camera.frame.clear_updaters()
        self.play(self.camera.frame.animate.restore())
        
        GEOM_COLOR=RED



        rotangle = -gear1.gen_params['tmax_ucut']
        self.play(Rotate(gear1, rotangle),
                  Rotate(rot_grp, rotangle, about_point=gear1.get_center()),
                  rack1.animate.shift(RIGHT * rotangle * gear1.rp), run_time=5)
        self.wait()

        Rbase = Circle(arc_center=gear1.get_center(),
                       radius=gear1.rb,
                       fill_opacity=0,
                       stroke_opacity=1,
                       stroke_color=TEAL)
        Rpitch = Circle(arc_center=gear1.get_center(),
                        radius=gear1.rp,
                        fill_opacity=0,
                        stroke_opacity=1,
                        stroke_color=TEAL)

        Rackline_1 = Line(start=gear1.get_center()+(gear1.rf)*DOWN+20*LEFT,
                          end=gear1.get_center()+(gear1.rf)*DOWN+20*RIGHT,
                          color=TEAL)

        RadLine = Line(start=gear1.get_center(),
                       end=gear1.get_center()+gear1.rb*np.array([np.sin(gear1.alpha*DEGREES),
                                                                 -np.cos(gear1.alpha*DEGREES),
                                                                 0]),
                       color=GEOM_COLOR)

        PitchLine = Line(start=gear1.get_center()+(gear1.rp)*DOWN+20*LEFT,
                         end=gear1.get_center()+(gear1.rp)*DOWN+20*RIGHT,
                         color=TEAL)
        VertLine = Line(start=gear1.get_center(),
                         end=gear1.get_center()+(gear1.rp)*DOWN,
                         color=GEOM_COLOR)

        Pitch_point = gear1.get_center()+gear1.rp*DOWN
        Right_angle_line = Line(start=RadLine.get_end(),
                                end=Pitch_point).match_style(RadLine)



        # self.add(Rbase,Rpitch,Rackline_1,RadLine,PitchLine,Right_angle_line,VertLine)
        self.play(FadeIn(Rbase))




        self.play(Invo_ref_mob.animate.set_stroke(opacity=1))
        self.wait()

        self.camera.frame.save_state()
        self.play(self.camera.frame.animate.move_to(RadLine.get_end()).scale(0.3))
        self.wait()
        self.play(self.camera.frame.animate.restore())



        self.play(self.camera.frame.animate.scale(1.5).shift(UP * 2))


        self.play(FadeIn(PitchLine), FadeIn(Rpitch))
        self.play(FadeIn(RadLine), FadeIn(VertLine))
        self.play(FadeIn(Right_angle_line))

        Dim_color=RED

        Alpha_dim_1 = Angle_Dimension_3point(arc_center=gear1.get_center(),
                                             start=gear1.get_center()+gear1.rp/4*DOWN,
                                             end=gear1.get_center()+gear1.rp/4*RadLine.get_unit_vector(),
                                             color=GEOM_COLOR,
                                             text=MathTex(r'\alpha',color=Dim_color))
        Rb_dim = Linear_Dimension(start=RadLine.get_start(),
                                  end=RadLine.get_end(),
                                  color=Dim_color,
                                  text=MathTex(r'r_b=r_p \cos(\alpha)',color=Dim_color),
                                  offset=1)

        Rp_dim = Linear_Dimension(start=VertLine.get_start(),
                                  end=VertLine.get_end(),
                                  color=Dim_color,
                                  text=MathTex(r'r_p', color=Dim_color),
                                  offset=-1.5)

        R_angle_dim = Linear_Dimension(start=Right_angle_line.get_start(),
                                       end=Right_angle_line.get_end(),
                                       color=Dim_color,
                                       text=MathTex(r'r_p \sin(\alpha)', color=Dim_color),
                                       offset=1.5)

        # self.play(Create(VGroup(Alpha_dim_1,Rb_dim,Rp_dim,R_angle_dim)))
        self.play(Create(Alpha_dim_1))
        self.wait()
        self.play(Create(Rp_dim))
        self.wait()
        self.play(Create(Rb_dim))
        self.wait()
        self.play(Create(R_angle_dim))
        self.wait()

        Sin2_Line = Line(start=RadLine.get_end(),
                         end=RadLine.get_end()+ DOWN * gear1.rp*np.sin(gear1.alpha*DEGREES)**2,
                         color=Dim_color
                         )
        Sin2_dim = Linear_Dimension(start=Sin2_Line.get_start(),
                                    end=Sin2_Line.get_end(),
                                    outside_arrow=True,
                                    text=MathTex(r'r_p \sin^2(\alpha)',color=Dim_color),
                                    color=Dim_color,
                                    offset=1)
        Sin2_dim['text'].rotate(-PI/2)
        Sin2_dim['text'].scale(R_angle_dim['text'].height/Sin2_dim['text'].height*0.93)
        Sin2_dim['text'].move_to(Sin2_dim['main_line'], aligned_edge=LEFT).shift(0.1 * RIGHT)

        self.wait()
        self.play(FadeOut(Rp_dim),FadeOut(Rb_dim),FadeOut(Alpha_dim_1))
        self.play(self.camera.frame.animate.scale(0.3).move_to(Sin2_dim['text']).shift(LEFT*2))
        self.play(Create(VGroup(Sin2_Line, Sin2_dim)))
        self.wait()

        self.play(FadeOut(Sin2_dim),
                  FadeOut(Sin2_Line),
                  FadeOut(R_angle_dim),

                  FadeOut(Right_angle_line),
                  FadeOut(R_angle_dim),

                  )

        self.play(self.camera.frame.animate.restore())
        self.play(self.camera.frame.animate.scale(1.5).shift(UP * 2))

        # R_ucut * cos(a) = rp-1.2m
        R_ucut = gear1.rf/np.cos(gear1.alpha*DEGREES)
        VR_ucut = (R_ucut*rotate_vector(DOWN,gear1.alpha*DEGREES))+gear1.get_center()

        # asd = Pointer_Label_Free(VR_ucut,'point',offset_vector=RIGHT+DOWN,color=GEOM_COLOR)

        Invo_anim_lim = 0.7
        Invo_path = VMobject()
        Invo_path.points=involute_point_gen(np.linspace(-Invo_anim_lim,0,10),gear1.rp)
        Invo_path.reverse_direction()
        Invo_path.rotate(-PI/2,about_point=ORIGIN).shift(gear1.get_center()).set_stroke(color=GEOM_COLOR,opacity=1)

        V_ofs_ucut = VR_ucut-(gear1.get_center()+DOWN*gear1.rp)

        Invo_ucut = VMobject()
        Invo_ucut.points = involute_point_gen(np.linspace(-Invo_anim_lim, 0, 10),
                                              gear1.rp,
                                              rad_offs=-V_ofs_ucut[1],
                                              tan_offs=V_ofs_ucut[0])
        Invo_ucut.reverse_direction()
        Invo_ucut.rotate(-PI / 2, about_point=ORIGIN).shift(gear1.get_center()).set_stroke(color=GEOM_COLOR, opacity=1)

        Invo_fi = ValueTracker(0)

        Ofs_v_line = Arrow(start=Invo_path.points[0,:],
                           end=Invo_ucut.points[0,:],
                           buff=0,
                           color=Dim_color)


        def Invo_func_ref(t):
            v =involute_func(t,gear1.rp)
            v = rotate_vector(v,-PI/2)
            v = v+gear1.get_center()
            return v

        def Invo_func_ucut_ref(t):
            v =involute_func(t,gear1.rp,
                             rad_offs=-V_ofs_ucut[1],
                             tan_offs=V_ofs_ucut[0])
            v = rotate_vector(v,-PI/2)
            v = v+gear1.get_center()
            return v

        Ofs_v_line.add_updater(lambda mob: Ofs_v_line.put_start_and_end_on(start=Invo_func_ref(Invo_fi.get_value()),
                                                                           end=Invo_func_ucut_ref(Invo_fi.get_value())))

        rack_bkp = rack1.copy()
        def Invo_update_poz(mob):
            mob.match_points(rack_bkp)
            mob.shift(Invo_func_ref(Invo_fi.get_value())-Invo_func_ref(0))
            mob.rotate(Invo_fi.get_value(),about_point=Invo_func_ref(Invo_fi.get_value()))

        rack1.add_updater(Invo_update_poz)
        PitchLine_bkp = PitchLine.copy()
        Rackline_1_bkp = Rackline_1.copy()
        PitchLine.add_updater(lambda mob: mob.match_points(PitchLine_bkp).rotate(Invo_fi.get_value(),
                                                                                 about_point=gear1.get_center()))
        Rackline_1.add_updater(lambda mob: mob.match_points(Rackline_1_bkp).rotate(Invo_fi.get_value(),
                                                                                   about_point=gear1.get_center()))

        self.wait()
        # self.add(Invo_path,Invo_ucut,Ofs_v_line)
        self.play(Create(Invo_path))
        self.wait()

        self.play(Invo_fi.animate.set_value(-Invo_anim_lim),run_time=2)
        self.wait()
        self.play(Invo_fi.animate.set_value(0), run_time=2)
        self.wait()
        self.play(Invo_fi.animate.set_value(-Invo_anim_lim), run_time=2)
        self.wait()
        self.play(Invo_fi.animate.set_value(0), run_time=2)
        self.wait()

        self.play(Create(Invo_ucut),
                  Create(Ofs_v_line))

        self.play(Invo_fi.animate.set_value(-Invo_anim_lim), run_time=4)
        self.wait()
        self.play(Invo_fi.animate.set_value(0), run_time=4)
        self.wait()


        ru_dim = Linear_Dimension(start=gear1.get_center()+gear1.rf*DOWN,
                                  end=VR_ucut,
                                  text=MathTex(r"r_d \tan(\alpha)",color=Dim_color),
                                  color=Dim_color,
                                  outside_arrow = True,
                                  offset=5)

        Rd_dim = Linear_Dimension(start=gear1.get_center(),
                                  end=gear1.get_center()+gear1.rf*DOWN,
                                  color=Dim_color,
                                  text=MathTex(r'r_d', color=Dim_color),
                                  offset=-0.5)

        # self.add(ru_dim)
        # self.play(Create(Rackline_1),
        #           Create(Rp_dim),
        #           Create(Alpha_dim_1),
        #           Create(Rd_dim),
        #           Create(ru_dim))

        self.play(Create(Rackline_1))
        self.wait()
        self.play(Create(Alpha_dim_1))
        self.wait()
        self.play(Create(Rp_dim))
        self.wait()
        self.play(Create(Rd_dim))
        self.wait()
        self.play(Create(ru_dim))
        self.wait()



        v_ofs_eq = MathTex(r'\begin{bmatrix}v_x\\v_y\end{bmatrix}','=',
                           r'\begin{bmatrix} r_d \tan(\alpha) \\ r_p-r_d\end{bmatrix}',
                           color=RED)
        v_ofs_eq.shift(LEFT*8,UP*5)

        self.play(Create(v_ofs_eq))
        self.wait()

        self.play(FadeOut(Rackline_1),
                  FadeOut(Rp_dim),
                  FadeOut(Alpha_dim_1),
                  FadeOut(Rd_dim),
                  FadeOut(ru_dim),
                  FadeOut(v_ofs_eq))


        HLDot1=Circle(radius=0.05,
                      fill_color=RED,
                      fill_opacity=1,
                      stroke_opacity=1,
                      stroke_color=TEAL)
        HLDot1.move_to(Undercut_line2.points[-1,:])
        HLDot2 = HLDot1.copy()
        HLDot2.move_to(Undercut_line2.points[0,:])

        self.play(Invo_fi.animate.set_value(-Invo_anim_lim), run_time=4)
        self.wait()
        self.play(Invo_fi.animate.set_value(0), run_time=4)
        self.wait()

        # self.play(FadeIn(HLDot1),
        #           FadeIn(HLDot2))
        self.wait()




class Undercut_calc(MovingCameraScene):
    def construct(self):
        # eq_cosa_1a = MathTex(r"r_b = r_p cos(\alpha)<r_d")
        eq_cosa_1 = MathTex(r"r_b = r_p \cos(\alpha) \le r_d")
        eq_cosa_2 = MathTex(r"\left( \frac{z}{2} m \right) \cos(\alpha)=\left( \frac{z}{2} m-1.25m \right)")
        eq_cosa_3 = MathTex(r"\frac{z}{2} m (\cos(\alpha)-1)=-1.25m")
        eq_cosa_4 = MathTex(r"z = \frac{2*1.25m}{m(1-\cos(\alpha))}")
        eq_cosa_5 = MathTex(r"z = 41.4543 \approx 42")
        # eq_cosa_2.move_to(eq_cosa_1,aligned_edge=UP)
        eq_grp = VGroup(eq_cosa_1,eq_cosa_2,eq_cosa_3,eq_cosa_4,eq_cosa_5)
        eq_grp.arrange(DOWN)

        eq_sin2a_1 = MathTex(r"h=r_p \sin^2(\alpha) \ge r_p-r_d")
        eq_sin2a_2 = MathTex(r"\left(\frac{z}{2} m\right) \sin^2(\alpha)=1.25m")
        eq_sin2a_3 = MathTex(r"z=\frac{2.5}{\sin^2(\alpha)} = 21.3715 \approx 22")

        eq_grp_2 = VGroup(eq_sin2a_1,eq_sin2a_2,eq_sin2a_3)
        eq_grp_2.arrange(DOWN)

        for mob in eq_grp:
            self.play(FadeIn(mob))
            self.wait(1)

        self.wait()
        self.play(FadeOut(eq_grp))
        self.wait()

        for mob in eq_grp_2:
            self.play(FadeIn(mob))
            self.wait(1)
        self.wait()
        self.play(FadeOut(eq_grp_2))
        m=1
        gear1 = Gear(50,m,cutout_teeth_num=50-6,make_smooth=False,stroke_width=2).rotate(-PI/2)
        gear1.shift(UP*gear1.rp)
        r_b = Circle(radius=gear1.rb,arc_center=gear1.get_center(),stroke_width=2,color=RED)
        # r_b.add_updater(lambda mob: mob.match_points(Circle(radius=gear1.rb,arc_center=gear1.get_center())))
        r_d = Circle(radius=gear1.rf,arc_center=gear1.get_center(),stroke_width=2,color=TEAL)
        # r_d.add_updater(lambda mob: mob.match_points(Circle(radius=gear1.rb, arc_center=gear1.get_center())))

        self.play(Create(gear1))
        self.add(r_b,r_d)
        self.wait()
        ztext = Text(f'z={gear1.z}').shift(UP*3)
        self.play(Create(ztext))
        gear2 = Gear(30,m,cutout_teeth_num=30-6,make_smooth=False,stroke_width=2).rotate(-PI/2)
        gear2.shift(UP * gear2.rp)
        r_b_2 = Circle(radius=gear2.rb, arc_center=gear2.get_center(), stroke_width=2, color=RED)
        r_d_2 = Circle(radius=gear2.rf, arc_center=gear2.get_center(), stroke_width=2, color=TEAL)

        self.play(Transform(gear1,gear2,make_smooth=False),
                  Transform(r_b,r_b_2),
                  Transform(r_d,r_d_2),
                  Transform(ztext,Text(f'z={gear2.z}').shift(UP*3)))
        self.wait()

        gear3 = Gear(10,m,cutout_teeth_num=10-6,make_smooth=False,stroke_width=2).rotate(-PI/2)
        gear3.shift(UP * gear3.rp)
        r_b_3 = Circle(radius=gear3.rb, arc_center=gear3.get_center(), stroke_width=2, color=RED)
        r_d_3 = Circle(radius=gear3.rf, arc_center=gear3.get_center(), stroke_width=2, color=TEAL)

        self.play(Transform(gear1, gear3),
                  Transform(r_b, r_b_3),
                  Transform(r_d, r_d_3),
                  Transform(ztext,Text(f'z={gear3.z}').shift(UP*3)))
        self.wait()

class Undercut_LOA(MovingCameraScene):
    def construct(self):
        m = 1.5
        a = 20*DEGREES
        v_a = rotate_vector(RIGHT,a)
        rt = ValueTracker(0)

        gear1 = Gear(22, m, h_f=1.25, stroke_opacity=0, fill_opacity=1, fill_color=WHITE, nppc=7, make_smooth=False)
        Invo_ref_pt = involute_point_gen(np.linspace(0, gear1.gen_params['tmax_invo'], 20), gear1.rb)
        Invo_ref_mob = VMobject()
        Invo_ref_mob.set_points(Invo_ref_pt)
        Invo_ref_mob.rotate(-gear1.angle_ofs - gear1.pitch_angle * 3 / 4 - gear1.pitch_angle * 22 / 4,
                            about_point=ORIGIN)
        Invo_ref_mob.set_color(RED)
        rz = 9
        gear2 = Gear(8, m, h_f=1.25, stroke_opacity=0, fill_opacity=1, fill_color=WHITE, nppc=7, make_smooth=False)

        def gear_updater(mob:Gear):
            mob.rotate(-mob.get_angle())
            mob.rotate(rt.get_value()*mob.pitch_angle-PI/2-mob.pitch_angle/2)

        gear1.add_updater(gear_updater)
        gear2.add_updater(gear_updater)

        Invo_ref_pt = involute_point_gen(np.linspace(0, gear1.gen_params['tmax_invo'], 20), gear1.rb)
        Invo_ref_mob = VMobject()
        Invo_ref_mob.set_points(Invo_ref_pt)
        Invo_ref_mob.rotate(-gear1.angle_ofs - gear1.pitch_angle *3 / 4 - gear1.pitch_angle*22/4, about_point=ORIGIN)
        Invo_ref_mob.set_color(RED)

        a0_invo = np.arctan2(Invo_ref_mob.points[0, 1]-gear1.get_center()[1], Invo_ref_mob.points[0, 0]-gear1.get_center()[0])
        def invo_updater(mob):
            angle_invo = np.arctan2(Invo_ref_mob.points[0, 1]-gear1.get_center()[1], Invo_ref_mob.points[0, 0]-gear1.get_center()[0])-a0_invo
            mob.rotate(-angle_invo,about_point=gear1.get_center())
            mob.rotate(rt.get_value()*gear1.pitch_angle,about_point=gear1.get_center())

        Invo_ref_mob.add_updater(invo_updater)


        rz = 9

        rack1 = Rack(rz, m,
                     h_a=1.25,
                     h_f=1,
                     stroke_opacity=0, fill_opacity=1, fill_color=BLUE_A).rotate(PI/2)
        # rack1.shift(rack1.pitch / 2 * DOWN)

        def rack_updater(mob):
            mob.shift(-mob.get_center())
            mob.shift(rt.get_value()*mob.pitch*RIGHT)

        rack1.add_updater(rack_updater)

        gear1.shift(UP * gear1.rp)
        Invo_ref_mob.shift(UP * gear1.rp)
        gear2.shift(UP * gear2.rp)

        LOA_line = Line(start=v_a*1.25*m/np.sin(a),
                        end=-v_a*1*m/np.sin(a),
                        color=RED)

        D_line = Line(start=UP*m*1.25+LEFT*20,
                      end=UP*m*1.25+RIGHT*20,
                      color=TEAL)

        P_line = Line(start=LEFT*20,
                      end=RIGHT*20,
                      color=TEAL)

        Rb_1 = Circle(radius=gear1.rb,
                      arc_center=gear1.get_center(),
                      color=TEAL)

        Rb_2 = Circle(radius=gear2.rb,
                      arc_center=gear2.get_center(),
                      color=TEAL)

        LOA_Dot = Circle(radius=0.1,fill_color=YELLOW,fill_opacity=0.5,stroke_color=RED)

        def LOA_Dot_updater(mob):
            lin_pitch = gear1.pitch*np.cos(a)
            v = v_a*(-0.75+rt.get_value())*lin_pitch
            mob.move_to(v)

        LOA_Dot.add_updater(LOA_Dot_updater)
        LOA_Dot.update()


        self.play(FadeIn(gear1),
                  FadeIn(rack1),
                  FadeIn(D_line),
                  FadeIn(P_line),
                  FadeIn(LOA_line),
                  FadeIn(Rb_1))
        # self.camera.frame.scale(3)
        self.play(Create(Invo_ref_mob))
        self.play(Create(LOA_Dot))
        self.wait()

        h_d_dim = Linear_Dimension(start=UP*1.25*m+LEFT*3,
                                   end=UP*0*m+LEFT*3,
                                   text='1.25m',
                                   offset=0,
                                   color=RED)

        h_dim_1 = Linear_Dimension(start=gear1.get_center()+gear1.rb*rotate_vector(DOWN,a),
                                   end=(gear1.get_center()+gear1.rb*rotate_vector(DOWN,a))[0]*RIGHT,
                                   text='h',
                                   offset=0,
                                   color=RED)

        end_len = 1.25*m/np.sin(a)/np.cos(a)/rack1.pitch
        self.play(rt.animate.set_value(end_len+0.75),run_time=10)
        self.wait()
        self.play(Create(h_dim_1), Create(h_d_dim))
        self.wait()
        self.play(FadeOut(h_dim_1))
        self.play(rt.animate.set_value(0.1), run_time=5)
        self.wait()


        self.play(FadeOut(gear1),FadeOut(Rb_1),FadeOut(Invo_ref_mob))
        gear2.z_index=-1
        Invo_ref_pt = involute_point_gen(np.linspace(0, gear2.gen_params['tmax_invo'], 20), gear2.rb)
        Invo_ref_mob.set_points(Invo_ref_pt)
        Invo_ref_mob.rotate(-gear2.angle_ofs - gear2.pitch_angle * 3 / 4 - gear2.pitch_angle *8/ 4,
                            about_point=ORIGIN)
        Invo_ref_mob.shift(UP*gear2.rp)
        a0_invo = np.arctan2(Invo_ref_mob.points[0, 1] - gear2.get_center()[1],
                             Invo_ref_mob.points[0, 0] - gear2.get_center()[0])
        def invo_updater2(mob):
            angle_invo = np.arctan2(Invo_ref_mob.points[0, 1]-gear2.get_center()[1], Invo_ref_mob.points[0, 0]-gear2.get_center()[0])-a0_invo
            mob.rotate(-angle_invo,about_point=gear2.get_center())
            mob.rotate(rt.get_value()*gear2.pitch_angle,about_point=gear2.get_center())
        Invo_ref_mob.clear_updaters()
        Invo_ref_mob.add_updater(invo_updater2)
        Invo_ref_mob.update()

        # Invo_ref_mob.set_color(RED)

        # self.camera.frame.scale(3)
        self.play(FadeIn(gear2),FadeIn(Rb_2),FadeIn(Invo_ref_mob))

        self.play(rt.animate.set_value(0.75+a/gear2.pitch_angle), run_time=6)
        self.wait()

        h_dim_2 = Linear_Dimension(start=gear2.get_center() + gear2.rb * rotate_vector(DOWN, a),
                                   end=(gear2.get_center() + gear2.rb * rotate_vector(DOWN, a))[0] * RIGHT,
                                   text='h',
                                   offset=1.5,
                                   outside_arrow=True,
                                   color=RED)

        self.play(Create(h_dim_2))
        self.wait()
        self.play(Uncreate(h_dim_2))
        self.play(rt.animate.set_value(end_len + 0.75), run_time=3)
        self.wait()



# with tempconfig({"quality": "medium_quality", "disable_caching": True, "save_last_frame":True}):
#     scene = Undercut_LOA()
#     scene.render()