import numpy as np
from manim import *
from manim_cad_drawing_utils import *
from manim.mobject.geometry.tips import *
from manim_gearbox import *
from scipy.optimize import fsolve
from scipy.optimize import fmin


def Contact_Ratio_Numeric(z1, z2, da_ofs=0, x1=0, x2=0):
    m = 1
    gear1 = Gear(z1, m, profile_shift=x1)
    gear2 = Gear(z2, m, profile_shift=x2)
    gear1.shift(gear1.rp * LEFT)
    gear2.shift((gear2.rp + da_ofs) * RIGHT)
    gear2.mesh_to(gear1,offset=da_ofs)
    da = np.linalg.norm(gear2.get_center()-gear1.get_center())
    r_ucut_1 = np.linalg.norm(involute_func(gear1.gen_params['tmin_invo'], gear1.rb))
    r_ucut_2 = np.linalg.norm(involute_func(gear2.gen_params['tmin_invo'], gear2.rb))
    alpha = np.arccos((gear1.rb + gear2.rb) / da)
    Loa_line = Line(start=gear1.get_center() + rotate_vector(gear1.rb * RIGHT, alpha),
                    end=gear2.get_center() + rotate_vector(gear2.rb * LEFT, alpha))
    ra1 = Circle(gear1.ra, arc_center=gear1.get_center())
    ra2 = Circle(gear2.ra, arc_center=gear2.get_center())
    rb1 = Circle(r_ucut_1, arc_center=gear1.get_center())
    rb2 = Circle(r_ucut_2, arc_center=gear2.get_center())
    ints1 = curve_intersection(Loa_line, ra1)
    ints2 = curve_intersection(Loa_line, ra2)
    intersections = []
    # if it's not empty
    # if ints1[0]:
    #     intersections.append(ints1[0])
    # # if it's not empty
    # if ints2[0]:
    #     intersections.append(ints2[0])
    if r_ucut_1 > gear1.rb * (1 + 1e-5):
        ints_u1 = curve_intersection(Loa_line, rb1)
        # if it's not empty
        if len(ints_u1[0]) > 0:
            if ints2[0]:
                if ints_u1[0][0] > ints2[0][0]:
                    intersections.append(ints_u1[0][0])
                else:
                    intersections.append(ints2[0][0])
            else:
                intersections.append(ints_u1[0][0])
        else:
            if len(ints2[0]) > 0:
                intersections.append(ints2[0][0])
            else:
                intersections.append(0)
    else:
        if len(ints2[0]) > 0:
            intersections.append(ints2[0][0])
        else:
            intersections.append(0)

    if r_ucut_2 > gear2.rb * (1 + 1e-5):
        ints_u2 = curve_intersection(Loa_line, rb2)
        # if it's not empty
        if len(ints_u2[0]) > 0:
            if len(ints1[0]) > 0:
                if ints_u2[0][0] < ints1[0][0]:
                    intersections.append(ints_u2[0][0])
                else:
                    intersections.append(ints1[0][0])
            else:
                intersections.append(ints_u2[0][0])
        else:
            if ints1[0]:
                intersections.append(ints1[0][0])
            else:
                intersections.append(1)
    else:
        if ints1[0]:
            intersections.append(ints1[0])
        else:
            intersections.append([1])
    # if r_ucut_2 > gear2.rb * (1 + 1e-5):
    #     ints_u2 = curve_intersection(Loa_line, rb2)
    #     if ints_u2[0]:
    #         intersections.append(ints_u2[0])

    intersections.sort()

    p1 = Loa_line.get_nth_curve_function(0)(intersections[0])
    p2 = Loa_line.get_nth_curve_function(0)(intersections[1])
    l = np.linalg.norm(p1 - p2)
    pb = gear1.pitch * np.cos(gear1.alpha * DEGREES)
    return l / pb


def Contact_Ratio_Equation(z1, z2, da_ofs=0, m=1):
    rp1 = z1 * m / 2
    ra1 = rp1 + m
    rb1 = rp1 * np.cos(20 * DEGREES)

    rp2 = z2 * m / 2
    ra2 = rp2 + m
    rb2 = rp2 * np.cos(20 * DEGREES)
    da = rp1 + rp2 + da_ofs

    l = np.sqrt(ra2 ** 2 - rb2 ** 2) + np.sqrt(ra1 ** 2 - rb1 ** 2) - np.sqrt(da ** 2 - (rb1 + rb2) ** 2)
    pb = m * PI * np.cos(20 * DEGREES)
    return l / pb


def calc_loa_lengths(gear1:Gear,gear2:Gear):
    da = np.linalg.norm(gear1.get_center()-gear2.get_center())
    alpha_act = np.arccos((gear1.rb + gear2.rb) / da)
    gamma = PI/2-alpha_act
    cosgamma = np.cos(gamma)
    rb1 = gear1.rb
    rb2 = gear2.rb
    ra1 = gear1.ra
    ra2 = gear2.ra
    rp1 = rb1 / (rb1 + rb2) * da
    rp2 = rb2 / (rb1 + rb2) * da
    b1 = 1 / 2 * (-2 * rp1 * cosgamma + np.sqrt((2 * rp1 * cosgamma) ** 2 - 4 * (rp1 ** 2 - ra1 ** 2)))
    b2 = 1 / 2 * (-2 * rp2 * cosgamma + np.sqrt((2 * rp2 * cosgamma) ** 2 - 4 * (rp2 ** 2 - ra2 ** 2)))
    return b1,b2

class LOA_Dots(MovingCameraScene):
    def construct(self):
        m = 1
        d_Ofs = ValueTracker(0.0)
        bt = ValueTracker(1)

        zt = ValueTracker(20)
        pst = ValueTracker(0)
        self.add(zt,pst)
        def gear1_gen():
            shift = pst.get_value()
            gear = Gear( (zt.get_value()),
                         module=m,
                         profile_shift=shift,
                         h_a=1.0,
                         h_f=1.2,
                         nppc=8,
                         make_smooth=False,
                         fill_opacity=1,
                         fill_color=BLUE_A,
                         stroke_opacity=0)
            gear.shift(UP * (gear.rp + (d_Ofs.get_value()+pst.get_value()) * m))
            gear.rotate((zt.get_value()%2)/2*gear.pitch_angle+PI/2)
            return gear
        gear1 = gear1_gen()
        # gear1.add_updater(lambda mob: mob.become(gear1_gen()))
        zt.set_value(20)

        def gear1_update(mob:Gear):
            gear = gear1_gen()
            upds=mob.updaters
            mob.__dict__ = gear.__dict__
            mob.updaters=upds
            return mob



        gear1.add_updater(lambda mob: mob.move_to(UP * (mob.rp + d_Ofs.get_value() * m)))
        gear1.update()
        gear2 = Gear(40, module=m, h_a=1.0, h_f=1.2, nppc=8,make_smooth=False,fill_opacity=1,fill_color=BLUE_B,stroke_opacity=0)

        gear2.shift(gear2.rp * DOWN)
        gear2.mesh_to(gear1)
        dot_op = ValueTracker(0.5)

        LOA = Line(start=LEFT,end=RIGHT,color=RED)
        def get_loa_param():
            bias = 1
            mesh_param = gear2.mesh_param_calc(gear1, offset=d_Ofs.get_value())
            alpha_act = np.arccos((gear1.rb + gear2.rb) / mesh_param['pitch_dist'])
            v12 = normalize( gear2.get_center()-gear1.get_center())
            p1 = gear1.get_center() + rotate_vector(v12, -alpha_act) * gear1.rb
            p2 = gear2.get_center() + rotate_vector(-v12, -alpha_act) * gear2.rb
            return p1, p2

        def get_dot_param():
            bias = 1
            mesh_param = gear1.mesh_param_calc(gear2, offset=d_Ofs.get_value())
            alpha_act = np.arccos((gear1.rb + gear2.rb) / mesh_param['pitch_dist'])
            v12 = normalize(gear2.get_center() - gear1.get_center())
            p1 = gear1.get_center() + rotate_vector(v12, -alpha_act) * gear1.rb
            p2 = gear2.get_center() + rotate_vector(-v12, -alpha_act) * gear2.rb
            loa_len = np.linalg.norm(p2-p1)

            v_alpha = rotate_vector(v12, -alpha_act+PI/2)
            loa_pitch=gear1.pitch*np.cos(gear1.alpha*DEGREES)
            dot0 = gear1.get_center() + mesh_param['rp1'] * v12

            loa_len_1, loa_len_2 = calc_loa_lengths(gear1, gear2)
            mod_ofs = 0.25+mesh_param['mod1']
            ofs_list = np.arange(-10+mod_ofs,11+mod_ofs)*loa_pitch
            dots = []
            for val in ofs_list:
                if val>-loa_len_2 and val<loa_len_1 and val<loa_len:
                    dots.append(dot0+v_alpha*val)

            dots2 = []
            r_ucut_1 = np.linalg.norm(involute_func(gear1.gen_params['tmin_invo'], gear1.rb))
            r_ucut_2 = np.linalg.norm(involute_func(gear2.gen_params['tmin_invo'], gear2.rb))
            for k in range(len(dots)):
                dot=dots[k]
                if not (np.linalg.norm(dot-p1)>loa_len or np.linalg.norm(dot-p2)>loa_len) \
                        and np.linalg.norm(dot-gear1.get_center())>r_ucut_1 \
                        and np.linalg.norm(dot-gear2.get_center())>r_ucut_2:
                    dots2.append(dots[k])



            return dots2

        def gen_dot_list(mob):
            dot_list = [Circle(arc_center=p,
                               radius=0.1,
                               color=RED_E,
                               stroke_opacity=1,
                               stroke_color=RED,
                               fill_opacity=dot_op.get_value()) for p in get_dot_param()]
            mob.become(VGroup(*dot_list))
        dot_grp=VGroup(*[Circle(arc_center=p,
                                radius=0.1,
                                color=RED_E,
                                stroke_opacity=1,
                                stroke_color=RED,
                                fill_opacity=dot_op.get_value()) for p in get_dot_param()])
        dot_grp.add_updater(gen_dot_list)




        LOA.add_updater(lambda mob: LOA.put_start_and_end_on(*get_loa_param()))
        LOA.update()

        # base1 = Arc(start_angle=PI+PI/4,angle=TAU/4, radius=gear1.rb,arc_center=gear1.get_center(),stroke_width=20,num_components=20,color=RED)
        base1 = Circle(color=RED,stroke_width=4)
        base1.add_updater(lambda mob: mob.match_points(Circle(radius=gear1.rb, arc_center=gear1.get_center(),num_components=60)))
        adden1 = Circle(color=BLUE,stroke_width=4)
        adden1.add_updater(lambda mob: mob.match_points(Circle(radius=gear1.ra, arc_center=gear1.get_center())))

        base2 = Circle(radius=gear2.rb, arc_center=gear2.get_center(), stroke_width=4, num_components=60, color=RED)
        base2.add_updater(lambda mob: mob.move_to(gear2.get_center()))
        adden2 = Circle(radius=gear2.ra, arc_center=gear2.get_center(), stroke_width=4,color=BLUE)
        adden2.add_updater(lambda mob: mob.move_to(gear2.get_center()))


        gear2.add_updater(lambda mob: mob.move_to(DOWN))
        gear2.add_updater(lambda mob: mob.mesh_to(gear1,offset=d_Ofs.get_value(),bias=bt.get_value()))
        gear1.update()
        gear2.update()

        d_Ofs.set_value(0.0)
        self.camera.frame.scale(1)
        # self.play(FadeIn(gear1,gear2))
        gear1.suspend_updating()
        gear2.suspend_updating()
        gear1.shift(UP*5)
        gear2.shift(DOWN*5)
        self.add(gear1,gear2)
        self.play(gear1.animate.shift(DOWN*5),gear2.animate.shift(UP*5),rate_func=rate_functions.ease_out_bounce,run_time=3)
        self.wait()

        gear1.resume_updating()
        gear2.resume_updating()

        self.play(Rotate(gear1, angle=gear1.pitch_angle*4), run_time=8)
        self.wait(1)
        self.play(d_Ofs.animate.set_value(0.25))
        self.play(Rotate(gear1, angle=gear1.pitch_angle*1), run_time=2)
        self.wait()
        self.play(bt.animate.set_value(-1))
        self.play(bt.animate.set_value(1))
        self.wait(1)
        self.play(d_Ofs.animate.set_value(0))
        self.wait(1)
        self.play(d_Ofs.animate.set_value(-0.15),bt.animate.set_value(0))
        self.wait(1)

        overlap = Intersection(gear1,gear2)
        overlap.set_fill(color=RED,opacity=1).set_stroke(width=1,color=RED)
        self.play(FadeIn(overlap),run_time=0.25)
        self.wait(0.5)
        self.play(FadeOut(overlap),run_time=0.25)



        self.play(d_Ofs.animate.set_value(0.25), bt.animate.set_value(1))
        self.wait()
        self.play(d_Ofs.animate.set_value(1.5))
        self.wait()
        self.play(d_Ofs.animate.set_value(0))
        #
        self.play(Create(VGroup(base1,base2)))
        LOA.update()
        self.play(Create(LOA))
        self.wait()
        self.play(d_Ofs.animate.set_value(0.0))
        self.wait()

        #

        def dim_gen():
            mesh_data = gear2.mesh_param_calc(gear1,offset=d_Ofs.get_value(),bias=bt.get_value())
            mid_point = gear1.get_center() + mesh_data['rp2']*DOWN
            mid_line = Dashed_line_mobject(Line(start=mid_point+LEFT*10,
                                                end=mid_point+RIGHT*10,
                                                color=RED),
                                           num_dashes=20,
                                           dash_offset=-0.25)
            vert_line = Dashed_line_mobject(Line(start=mid_point+DOWN*10,
                                                end=mid_point+UP*10,
                                                color=RED),
                                            dash_offset=-0.25,
                                            num_dashes=20)
            alpha_dim = Angle_Dimension_3point(start=mid_point-LOA.get_unit_vector()*2.5,
                                               end=mid_point+LEFT*2.5,
                                               arc_center=mid_point,
                                               offset=0.5,
                                               color=RED,
                                               text=MathTex(r'\alpha',color=RED))
            return VGroup(mid_line,vert_line,alpha_dim)

        alpha_dim_grp = dim_gen()
        alpha_dim_grp.add_updater(lambda mob: mob.become(dim_gen()))
        alpha_dim_grp.update()

        self.play(Create(alpha_dim_grp))
        self.wait()

        self.play(d_Ofs.animate.set_value(0))

        self.wait(1)
        self.play(d_Ofs.animate.set_value(1))
        self.wait()
        alpha_dim_grp.suspend_updating()
        self.play(FadeOut(alpha_dim_grp),run_time=2)

        self.wait()
        self.play(d_Ofs.animate.set_value(0))
        self.wait()
        #
        self.play(Create(dot_grp))
        self.wait(1)
        self.play(Rotate(gear1, angle=gear1.pitch_angle), run_time=8)


        self.play(d_Ofs.animate.set_value(0.5))
        self.wait(1)
        self.play(Rotate(gear1, angle=gear1.pitch_angle), run_time=8)
        self.wait(1)

        self.play(d_Ofs.animate.set_value(1.5))
        self.wait(1)
        self.play(Rotate(gear1, angle=gear1.pitch_angle), run_time=8)
        self.wait(1)



        #
        # #
        self.play(Create(VGroup(adden1,adden2)))
        self.wait()

        self.play(d_Ofs.animate.set_value(0.5))
        self.wait(1)
        #
        self.play(Rotate(gear1, angle=gear1.pitch_angle*1.3), run_time=8)
        self.wait(1)

        def loa_act_gen():
            loa_endpoint_1 = LOA.get_nth_curve_function(0)(curve_intersection(LOA,adden1)[0])
            loa_endpoint_2 = LOA.get_nth_curve_function(0)(curve_intersection(LOA,adden2)[0])
            return Linear_Dimension(start=loa_endpoint_1,
                                    end=loa_endpoint_2,
                                    color=RED,
                                    offset=2.5,
                                    text='Line of Contact')

        dim_loa_act = loa_act_gen()
        dim_loa_act.add_updater(lambda mob: mob.become(loa_act_gen()))
        dim_pitch = Linear_Dimension(start=dot_grp[0].get_center(),
                                     end=dot_grp[1].get_center(),
                                     color=RED,
                                     offset=1.25,
                                     text='Base Pitch'
                                     )
        dim_pitch.add_updater(lambda mob: mob.become(Linear_Dimension(start=dot_grp[0].get_center(),
                                                                      end=dot_grp[1].get_center(),
                                                                      offset=1.25,
                                                                      color=RED,
                                                                      text='Base Pitch'
                                                                      )))

        # eq_c_rat = MathTex(r"\epsilon = \frac{Line \  of \  Contact}{Base \  Pitch} > 1", color=RED).shift(UP*3)
        eq_c_rat = MathTex(r"\mbox{Contact Ratio} = \frac{ \mbox{Line  of  Contact}}{\mbox{Base Pitch}} > 1", color=RED).shift(UP * 3)
        self.play(Create(dim_loa_act))
        self.wait()
        self.play(Create(dim_pitch))
        self.wait()

        self.play(Create(eq_c_rat))
        self.wait()

        self.play(d_Ofs.animate.set_value(0),Rotate(gear1,angle=gear1.pitch_angle*(-0.3)))
        self.wait()

        self.play(Uncreate(eq_c_rat))
        self.wait()


        # UNDERCUT




        dim_pitch.suspend_updating()
        dim_loa_act.suspend_updating()
        self.play(FadeOut(dim_loa_act),
                  FadeOut(dim_pitch))

        Ucut_text = Text('!Undercut!',color=RED).shift(UP*2.8).scale(0.9)

        self.play(Create(Ucut_text))
        self.wait()

        gear1.clear_updaters()
        gear1.add_updater(gear1_update)
        gear2.suspend_updating()
        self.play(zt.animate.set_value(6),run_time=8)
        gear2.resume_updating()
        self.wait()

        self.play(Uncreate(Ucut_text))
        self.wait()

        r_ucut = np.linalg.norm(involute_func(gear1.gen_params['tmin_invo'], gear1.rb))
        base1_ucut = Circle(r_ucut,arc_center=gear1.get_center(),color=RED,stroke_width=2)

        self.play(self.camera.frame.animate.scale(0.5))
        self.play(Create(base1_ucut))
        self.wait()

        gear1.remove_updater(gear1_update)
        self.play(Rotate(gear1, angle=gear1.pitch_angle*0.5), run_time=8)
        self.wait()
        self.play(Rotate(gear1, angle=gear1.pitch_angle * 0.5), run_time=8)
        self.wait()
        self.play(Uncreate(base1_ucut))
        # self.play(self.camera.frame.animate.scale(2))
        self.wait()

        gear1.add_updater(gear1_update)
        self.play(pst.animate.set_value(0.65),run_time=4)
        self.wait()
        gear1.remove_updater(gear1_update)
        self.play(Rotate(gear1, angle=gear1.pitch_angle*1.41), run_time=8)
        self.wait(1)
        # self.remove(gear1)
        # self.add(gear3)

class LOA_Calc(MovingCameraScene):
    def construct(self):
        m=0.25
        ofs=0.1
        gear1 = Gear(20,m,fill_opacity=0.3,stroke_opacity=0,nppc=10)
        gear2 = Gear(40,m,fill_opacity=0.3,stroke_opacity=0,nppc=10)
        gear1.shift(RIGHT*(gear1.rp+ofs/2))
        gear2.shift(LEFT*gear2.rp)
        gear2.mesh_to(gear1,offset=+ofs/m,bias=-1)

        # self.add(NumberPlane(background_line_style={
        #                                             "stroke_color": TEAL,
        #                                             "stroke_width": 2,
        #                                             "stroke_opacity": 0
        #                                            },
        #                      faded_line_style={
        #                                        "stroke_color": TEAL,
        #                                        "stroke_width": 2,
        #                                        "stroke_opacity": 0
        #                                       },
        #                      faded_line_ratio=5))
        self.play(FadeIn(VGroup(gear1,gear2)))

        # a=20*DEGREES
        mesh_param=gear2.mesh_param_calc(gear1,offset=ofs/m)
        Rb = gear2.rb
        rb = gear1.rb
        a = np.arccos((Rb+rb)/(np.linalg.norm(gear1.get_center()-gear2.get_center())))
        a2 = PI/2-a

        rb_circ = Circle(arc_center=gear1.get_center(),
                         radius=gear1.rb,
                         num_components=100)
        Rb_circ = Circle(arc_center=gear2.get_center(),
                         radius=gear2.rb,
                         num_components=100)
        ra_circ = Circle(arc_center=gear1.get_center(),
                         radius=gear1.ra,
                         num_components=100)
        Ra_circ = Circle(arc_center=gear2.get_center(),
                         radius=gear2.ra,
                         num_components=100)

        circ_grp = VGroup(Ra_circ,ra_circ,rb_circ,Rb_circ)
        circ_grp.set_stroke(color=PURPLE,opacity=0.3)
        rb_line = Line(start=ORIGIN,
                       end=gear1.rb*rotate_vector(np.array([np.cos(a),np.sin(a),0]),PI)).shift(gear1.get_center())

        Rb_line = Line(start=ORIGIN,
                       end=gear2.rb * rotate_vector(np.array([np.cos(a), np.sin(a), -0]), -PI/2*0)).shift(
            gear2.get_center())

        LOA_line = Line(start=Rb_line.get_end(),
                        end=rb_line.get_end())

        gamma1 = np.arccos(gear1.rb/gear1.ra)
        gamma2 = np.arccos(gear2.rb/gear2.ra)

        ra_line = rb_line.copy().rotate(-gamma1,about_point=gear1.get_center())\
            .scale(gear1.ra/gear1.rb,about_point=gear1.get_center())
        Ra_line = Rb_line.copy().rotate(-gamma2, about_point=gear2.get_center())\
            .scale(gear2.ra / gear2.rb,about_point=gear2.get_center())

        da_line = Line(start=gear1.get_center(),end=gear2.get_center())

        self.play(Create(circ_grp))
        self.play(Create(VGroup(rb_line,Rb_line)))
        self.play(Create(VGroup(ra_line,Ra_line)))
        self.play(Create(LOA_line))
        self.wait()
        shift_v = -(rb_line.get_end()-rb_line.get_start())

        LOA_line_2=LOA_line.copy()
        self.add(LOA_line_2)
        self.play(LOA_line.animate.shift(shift_v),
                  Rb_line.animate.put_start_and_end_on(start=gear2.get_center(), end=Rb_line.get_end() + shift_v),)
                  # rb_line.animate.put_start_and_end_on(start=gear1.get_center(), end=rb_line.get_end() + shift_v))

        dim_rbrb = Linear_Dimension(start=Rb_line.get_start(),
                                    end=Rb_line.get_end(),
                                    text=MathTex(r'r_{b1} + r_{b2}',color=RED),
                                    offset=1,
                                    color=RED)
        dim_da = Linear_Dimension(start=gear1.get_center(),
                                  end=gear2.get_center(),
                                  text=MathTex(r'd_a',color=RED),
                                  color=RED)

        dim_l0 = Linear_Dimension(start=LOA_line.get_start(),
                                  end=LOA_line.get_end(),
                                  text=MathTex(r'l_0',color=RED),
                                  color=RED)
        l0_eq = MathTex(r'l_0=\sqrt{d_a^2-(r_{b1}+r_{b2})^2}')
        l0_eq.shift(RIGHT*4+DOWN*3)


        self.play(Create(dim_rbrb))
        self.play(Create(dim_da))
        self.play(Create(dim_l0))
        self.play(Create(l0_eq))
        self.wait()

        self.play(FadeOut(dim_rbrb),
                  FadeOut(dim_da),
                  FadeOut(dim_l0),
                  FadeOut(l0_eq))
        self.wait()

        self.play(LOA_line.animate.shift(-shift_v),
                  Rb_line.animate.put_start_and_end_on(start=gear2.get_center(), end=Rb_line.get_end() - shift_v) )
        self.remove(LOA_line_2)
        self.wait()

        dim_Rb = Linear_Dimension(start=Rb_line.get_start(),
                                  end=Rb_line.get_end(),
                                  text=MathTex(r'r_{b1}',color=RED),
                                  color=RED)
        dim_Ra = Linear_Dimension(start=Ra_line.get_start(),
                                  end=Ra_line.get_end(),
                                  text=MathTex(r'r_{a1}',color=RED),
                                  color=RED,
                                  offset=-1)

        dim_l1 = Linear_Dimension(start=Rb_line.get_end(),
                                  end=Ra_line.get_end(),
                                  text=MathTex(r'l_1',color=RED),
                                  color=RED,
                                  offset=1)

        l1_eq = MathTex(r'l_1=\sqrt{r_{a1}^2-r_{b1}^2}').shift(RIGHT*4+DOWN*3)

        self.play(Create(dim_Ra))
        self.play(Create(dim_Rb))
        self.play(Create(dim_l1))
        self.play(Create(l1_eq))
        self.wait()

        self.play(FadeOut(dim_Ra),
                  FadeOut(dim_Rb),
                  FadeOut(dim_l1),
                  FadeOut(l1_eq))


        dim_l2 = Linear_Dimension(start=rb_line.get_end(),
                                  end=ra_line.get_end(),
                                  text=MathTex(r'l_2',color=RED),
                                  offset=-3,
                                  color=RED,
                                  outside_arrow=True)

        l2_eq = MathTex(r'l_2=\sqrt{r_{a2}^2-r_{b2}^2}').shift(RIGHT * 4 + DOWN * 3)

        self.play(Create(dim_l2),
                  Create(l2_eq))
        self.wait()

        self.play(FadeOut(dim_l2),
                  FadeOut(l2_eq))

        L1_line = Line(start=rb_line.get_end(),
                       end=ra_line.get_end(),
                       color=RED)
        L2_line = Line(start=Rb_line.get_end(),
                       end=Ra_line.get_end(),
                       color=RED)

        Contact_line = Line(start=Ra_line.get_end(),
                            end=ra_line.get_end(),
                            color=RED)



        lc_eq = MathTex('l=','\sqrt{r_{a2}^2-r_{b2}^2}', '+', '\sqrt{r_{a1}^2-r_{b1}^2}', '-', r' \sqrt{d_a^2 - \left(r_{b1}+r_{b2}\right)^2}')
        # lc_eq = MathTex(r'l=\sqrt{r_a^2-r_b^2} + \sqrt{R_a^2-R_b^2} - \sqrt{d_a^2 - \left(r_b+R_b\right)^2}')
        # l =\sqrt{r_a ^ 2 - r_b ^ 2} + \sqrt{R_a ^ 2 - R_b ^ 2} - \sqrt{d_a ^ 2 - R_b ^ 2 + r_b ^ 2}
        lc_eq.shift(UP*3)

        lc_blackbox = SurroundingRectangle(lc_eq,buff=0.3,fill_color=BLACK,fill_opacity=0.75,stroke_opacity=0)
        lc_blackbox.stretch(3,0)

        self.play(FadeIn(lc_blackbox))
        L1_line.set_stroke(width=8)
        L2_line.set_stroke(width=8)
        self.play( FadeIn(L1_line,lc_eq[1]))
        self.wait()
        self.play(FadeOut(L1_line))
        self.wait()
        self.play(FadeIn(L2_line,lc_eq[2:4]))
        self.wait()
        self.play(FadeOut(L2_line))
        self.wait()
        self.play(LOA_line.animate.set_stroke(color=RED,width=8),FadeIn(lc_eq[4:]))
        self.wait()
        Contact_line.set_stroke(width=8)
        self.play(LOA_line.animate.set_stroke(color=WHITE,width=DEFAULT_STROKE_WIDTH),FadeIn(Contact_line,lc_eq[0]))
        self.wait()
        # self.play(FadeIn(Contact_line,lc_eq[0]))
        self.wait()
        self.play(FadeOut(Contact_line))

class Pitch_loa_calc(MovingCameraScene):
    def construct(self):
        m=0.5
        gear1 = Gear(20,m,make_smooth= False)
        gear1.shift(LEFT*gear1.rp)

        Invo_1 = VMobject()
        Invo_1.points = involute_point_gen(np.linspace(0,1,15),r=gear1.rb)

        Invo_1.rotate(-gear1.angle_ofs-gear1.pitch_angle/4,about_point=ORIGIN)
        Invo_1.shift(gear1.get_center())
        Invo_1.set_stroke(color=RED)

        Invo_2 = Invo_1.copy().rotate(gear1.pitch_angle,about_point=gear1.get_center())
        Circ_rp = Circle(radius=gear1.rp,
                         arc_center=gear1.get_center(),
                         color=TEAL)
        Circ_rb = Circle(radius=gear1.rb,
                         arc_center=gear1.get_center(),
                         color=TEAL)

        Arc_rp = Arc(radius=gear1.rp,
                      arc_center=gear1.get_center(),
                      num_components=40,
                      color=YELLOW,
                      start_angle=-gear1.pitch_angle/4,
                      angle=gear1.pitch_angle)
        Arc_rb = Arc(radius=gear1.rb,
                      arc_center=gear1.get_center(),
                      num_components=40,
                      color=YELLOW,
                      start_angle=-gear1.pitch_angle/4-gear1.angle_ofs,
                      angle=gear1.pitch_angle)


        def gear_invo_point(t,angle):
            p = rotate_vector(involute_func(t,gear1.rb),angle) + gear1.get_center()
            return p

        loa_pitch_line = Line(start=gear_invo_point(20*DEGREES,-gear1.pitch_angle/4-gear1.angle_ofs),
                              end=gear_invo_point(0,gear1.pitch_angle-gear1.pitch_angle/4-gear1.angle_ofs),
                              color=YELLOW)

        pt = ValueTracker(0)
        def loa_update(mob: Line):
            # p1 = gear_invo_point(20*DEGREES+pt.get_value(),-gear1.pitch_angle/4-gear1.angle_ofs)
            # p2 = gear_invo_point(pt.get_value()+gear1.pitch_angle/4+gear1.angle_ofs,gear1.pitch_angle-gear1.pitch_angle/4-gear1.angle_ofs)
            p1 = gear_invo_point(pt.get_value()+20*DEGREES,0)
            p2 = gear_invo_point(pt.get_value(), gear1.pitch_angle)
            mob.put_start_and_end_on(p1,p2)
            mob.rotate(-gear1.pitch_angle/4-gear1.angle_ofs,about_point=gear1.get_center())
        loa_pitch_line.add_updater(loa_update)
        loa_pitch_line.update()

        eq_1 = MathTex(r'p=m\pi')
        eq_2 = MathTex(r'p_b=m\pi \frac{r_b}{r_p}')
        eq_3 = MathTex(r'p_b=m\pi \cos(\alpha_0)')
        eq_grp = VGroup(eq_1, eq_2, eq_3)
        eq_grp.arrange(DOWN).shift(RIGHT * 4 + UP * 1)

        self.play(Create(VGroup(gear1,Circ_rp,Circ_rb,Invo_1,Invo_2)))
        self.wait()
        self.play(Create(Arc_rp),Create(eq_1))
        self.wait()
        self.play(Transform(Arc_rp,Arc_rb))
        self.play(Create(eq_2),Create(eq_3))
        self.wait()
        self.play(ReplacementTransform(Arc_rp, loa_pitch_line))
        self.wait()
        self.play(pt.animate.set_value(1-20*DEGREES))
        self.wait()

        # self.play(Create(eq_grp))
        # self.wait()

class Contact_Ratio(MovingCameraScene):
    def construct(self):
        # eq_contact = MathTex(r'\epsilon = \frac{l}{p_b}=\frac{\sqrt{r_{a2}^2-r_{b2}^2} + \sqrt{r_{a1}^2-r_{b1}^2} - \sqrt{d_a^2 - \left(r_{b1}+r_{b2}\right)^2}}{\frac{m\pi }{\cos(\alpha_0)}}')
        eq_contact = MathTex(
            r'\epsilon = \frac{l}{p_b}=\frac{\sqrt{r_{a2}^2-r_{b2}^2} + \sqrt{r_{a1}^2-r_{b1}^2} - \sqrt{d_a^2 - (r_{b1}+r_{b2})^2}}{m\pi \cos(\alpha_0)}')
        eq_ra = MathTex(r'r_a = m \frac{z+2}{2} + xm')
        eq_rb = MathTex(r'r_b = m \frac{z}{2} \cos(\alpha_0)')
        eq_da = MathTex(r'd_a = r_{p1}+r_{p2}+\delta m = m \frac{z_1+z_2}{2} +\delta m')
        eq_grp = VGroup(eq_ra,eq_rb,eq_da,eq_contact).arrange(DOWN)
        self.play(Create(eq_grp),run_time=4)
        self.wait()


class Contact_table_z(MovingCameraScene):
    def construct(self):
        z1 = 41
        z2 = np.arange(6,24)


        contacts = [Contact_Ratio_Numeric(z1,z2val) for z2val in z2]

        contacts_optim = []
        for zval,c in zip(z2,contacts):
            if zval<25:
                x = fmin(lambda x2: -Contact_Ratio_Numeric(z1,zval,x2=x2),0)
                x = x if x>0 else [0.0]
                contact = Contact_Ratio_Numeric(z1,zval,x2=x[0])
                contacts_optim.append(contact)
            else:
                contacts_optim.append(c)

        # contacts_optim=contacts.copy()

        self.wait()

        chart1 = BarChart(contacts,
                          y_range=[0,2,0.2],
                          bar_colors=[BLUE]*len(contacts_optim))
        chart2 = BarChart(contacts_optim,
                          y_range=[0,2,0.2],
                          bar_colors=[RED]*len(contacts_optim),
                          bar_names=list(map(str,z2)))
        chart2.add(chart1.bars)
        numlines = NumberPlane(x_range=[6,24],
                               x_axis_config=chart2.x_axis_config,
                               x_length=chart2.x_length,
                               y_range=[0,2,0.2],
                               y_length=chart2.y_length,
                               y_axis_config=chart2.y_axis_config)
        chart2.add_to_back(numlines.background_lines)
        chart2.shift(RIGHT*0.5)



        graphtitle = Text('Contact ratio vs tooth number')


        subtext = Text('*Mating with a 41-tooth gear',font_size=34)
        subtext.next_to(chart2,DOWN,buff=0.3)
        subtext2 = VGroup(Text('Profile shifted',color=RED,font_size=34),
                          Text('Standard',color=BLUE,font_size=34)
                          ).arrange(DOWN,buff=0.1)
        subtext2.next_to(chart2,UP,buff=0.15)
        graphtitle.next_to(subtext2, UP, buff=0.1)

        self.camera.frame.scale(1.15)
        self.play(Create(chart2))
        self.play(Create(graphtitle))
        self.play(Create(subtext))
        self.play(Create(subtext2))
        # self.play(Create(chart1))
        self.wait(1)

class Contact_plot_d(MovingCameraScene):
    def construct(self):
        z1 = 41
        z2 = 20
        ax = Axes(x_range=[0, 2, 0.1],
                  y_range=[0, 2, 0.2],
                  tips=False,
                  axis_config={"include_numbers": True},
                  # y_axis_config={"label_constructor": MathTex(r'\epsilon')},
                  )
        numlines = NumberPlane(x_range=[0, 2, 0.1],
                               y_range=[0, 2, 0.2],
                               x_axis_config={"include_numbers": False},
                               x_length=ax.x_length,
                               y_length=ax.y_length,
                               y_axis_config={"include_numbers": False})
        contact_func =  ax.plot(lambda t: Contact_Ratio_Equation(z1, z2, t, m=1),[0,2])
        ref_one_func = ax.plot(lambda t: 1,[0,2],color=RED)
        xlabel = MathTex('\delta').scale(2).next_to(ax,DOWN,buff=0.3)
        ylabel = MathTex('\epsilon').scale(2).next_to(ax, LEFT, buff=0.3)
        # contact_func_n = ax.plot(lambda t: Contact_Ratio_Numeric(z1, z2, t), [0, 2])

        title_tx = Tex(r'Contact ratio of gear pair 20:41 vs axial distance')
        title_tx.scale(1.5).next_to(ax,UP,buff=0.3)
        self.camera.frame.scale(1.3)
        # self.add(numlines,ax,contact_func,ref_one_func,ylabel,xlabel,title_tx)
        self.play(Create(ax),Create(numlines))
        self.play(Create(ylabel),
                  Create(xlabel),
                  Create(title_tx))
        self.play(Create(contact_func))
        self.play(Create(ref_one_func))
        self.wait()
        self.play(FadeOut(ax),
                  FadeOut(numlines),
                  FadeOut(ylabel),
                  FadeOut(xlabel),
                  FadeOut(title_tx),
                  FadeOut(ref_one_func),
                  FadeOut(contact_func))


class Gear_size_comp_2(MovingCameraScene):
    def construct(self):
        k=3
        m=0.2/k
        gear1 = Gear(7*k,m,profile_shift=0.7,fill_color=TEAL,fill_opacity=1,stroke_opacity=0,make_smooth=False)
        gear2 = Gear(70*k,m,fill_color=WHITE,fill_opacity=1,stroke_opacity=0,make_smooth=False)
        gear2.shift(LEFT*gear2.rp*0+RIGHT*self.camera.frame.width/2).rotate(gear2.pitch_angle/2)
        gear1.shift(LEFT*10)
        gear1.mesh_to(gear2)
        self.add(gear1,gear2)


# with tempconfig({"quality": "medium_quality", "disable_caching": True, "save_last_frame":True}):
#     scene = Contact_Ratio()
#     scene.render()