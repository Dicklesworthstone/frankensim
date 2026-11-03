//! Custom section planes in the real native wake owners and mechanical stages.
use fs_airscrew::dynamic_inflow::first_harmonic::HarmonicInduction;
use fs_flyer::aircraft::{OpenLoopDesign, wright_openloop_v1};
use fs_flyer::propcoupling::{CANDIDATE_B, PropCouplingSolverSpec, PropDisk};
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor, GearedState, GearedStep};
use fs_flyer::spine::rotors::propulsion::dynamic::first_harmonic::{
    HarmonicPropulsion, HarmonicState,
};
use fs_flyer::spine::rotors::propulsion::dynamic::{DynamicPropulsion, DynamicState};
use fs_flyer::spine::rotors::propulsion::surfaces::{
    InducedDragGroup, LiftingSurfaces, SectionProperties,
};
use fs_flyer::spine::rotors::propulsion::swirl::SwirlPropulsion;
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};
use fs_wing::nonlinear::continuation::ContinuationControls;
use fs_wing::nonlinear::{
    InfluenceOperator, SectionFrame, StripRegime, StripSpec, solve_nonlinear_with_section_frames,
};
use fs_wing::{Panel, SurfaceId, flat_surface};

const ZERO: [f64; 3] = [0.0; 3];
const ORIGIN: [f64; 3] = [0.4, -0.8, 0.3];
fn vec(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn orientation(state: &GearedState) -> So3 {
    let q = state.rigid.quat;
    So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap()
}
fn inputs() -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.025,
        motor_nm: 80.0,
        gravity_world_mps2: ZERO,
    }
}
fn domain() -> ContinuationControls {
    ContinuationControls {
        maximum_incidence_rad: 0.6,
        maximum_incidence_change_rad: 0.08,
        maximum_row_condition: 16.0,
        minimum_trailing_clearance_chords: 0.02,
        maximum_newton_iterations: 32,
        maximum_residual_evaluations: 600,
        maximum_jacobian_products: 4096,
        section_difference_rad: 1e-6,
    }
}

struct Fixture {
    design: OpenLoopDesign,
    assembly: GearedAssembly<2>,
    panels: Vec<Panel>,
    strips: Vec<StripSpec>,
    properties: Vec<SectionProperties>,
    groups: [InducedDragGroup; 1],
}
impl Fixture {
    fn new() -> Self {
        let mut design = wright_openloop_v1();
        design.cg_m = ZERO;
        design.drag_z_m = 0.0;
        design.coupling = PropCouplingSolverSpec {
            tol: 1e-11,
            cap: 64,
            ..CANDIDATE_B
        };
        design.disks = [
            PropDisk {
                center_m: [-2.0, -1.4, 0.0],
                omega_rad_s: 35.0,
            },
            PropDisk {
                center_m: [-2.0, 1.4, 0.0],
                omega_rad_s: 35.0,
            },
        ];
        let assembly = GearedAssembly::new(
            design.gross_mass_kg,
            Mat3 {
                m: [1000.0, 10.0, 0.0, 10.0, 400.0, 5.0, 0.0, 5.0, 1200.0],
            },
            core::array::from_fn(|k| GearedRotor {
                axis_body: [1.0, 0.0, 0.0],
                axial_inertia_kgm2: [1.5, 2.5][k],
                ratio: [1.0, -1.0][k],
            }),
        )
        .unwrap();
        let rotation = So3::exp(So3Tangent::new(Vec3::new(
            core::f64::consts::FRAC_PI_2,
            0.0,
            0.0,
        )))
        .unwrap();
        let rotate = |value| array(rotation.rotate(vec(value)).unwrap());
        let panels = flat_surface(SurfaceId::Vertical, 4.0, 1.0, -1.0, 0.0, 4, 2)
            .unwrap()
            .into_iter()
            .map(|panel| Panel {
                a: rotate(panel.a),
                b: rotate(panel.b),
                ctrl: rotate(panel.ctrl),
                normal: rotate(panel.normal),
                ..panel
            })
            .collect();
        let strips = (0..4)
            .map(|s| StripSpec {
                panel_indices: vec![s, s + 4],
                chord_m: 1.0,
                twist_rad: 0.01,
            })
            .collect();
        let properties = vec![
            SectionProperties {
                camber_ratio: 0.02,
                attached_cd: 0.015,
                canard_gain: 1.0,
                induced_drag_group: 0
            };
            4
        ];
        Self {
            design,
            assembly,
            panels,
            strips,
            properties,
            groups: [InducedDragGroup {
                span_m: 4.0,
                oswald_e: 0.8,
            }],
        }
    }

    fn model<'a>(&'a self, surfaces: &'a LiftingSurfaces<'a>) -> GearedPropulsion<'a, 2> {
        GearedPropulsion::new_spatial(
            &self.design,
            &self.assembly,
            core::array::from_fn(|k| PropellerMount {
                rotor_index: k,
                thrust_direction_body: [1.0, 0.0, 0.0],
                handedness: if k == 0 {
                    PropellerHandedness::Positive
                } else {
                    PropellerHandedness::Negative
                },
            }),
        )
        .unwrap()
        .with_lifting_surfaces(surfaces)
        .unwrap()
    }

    fn state(&self) -> GearedState {
        let rotation = So3::exp(So3Tangent::new(Vec3::new(0.1, -0.05, 0.2))).unwrap();
        let q = rotation.as_quat();
        self.assembly
            .state_from_speed(
                SixDofState {
                    pos_m: [2.0, -3.0, -10.0],
                    vel_mps: array(rotation.rotate(vec([12.0, -0.55, 0.2])).unwrap()),
                    quat: [q.w, q.x, q.y, q.z],
                    omega_body: [0.03, -0.02, 0.04],
                },
                0.0,
                35.0,
            )
            .unwrap()
    }

    // Independently solve the native section equations at the accepted physical
    // field. This detects a body-xz fallback in any of the actual wake consumers.
    fn assert_wing(&self, state: &GearedState, gamma: &[f64], alpha: &[f64], local: &[[f64; 3]]) {
        assert_eq!(gamma.len(), self.panels.len());
        assert!(alpha.iter().all(|a| a.abs() < 0.30));
        let freestream = array(
            orientation(state)
                .inverse()
                .rotate(vec(state.rigid.vel_mps))
                .unwrap(),
        );
        let op = InfluenceOperator::build(&self.panels, freestream, 1.225).unwrap();
        let mut strips = self.strips.clone();
        for (strip, section) in strips.iter_mut().zip(&self.properties) {
            strip.twist_rad += section.canard_gain * inputs().canard_rad;
        }
        let frames: Vec<_> = strips
            .iter()
            .map(|strip| SectionFrame::from_panel(&self.panels[strip.panel_indices[0]]).unwrap())
            .collect();
        let closure = |s: usize, incidence: f64| {
            (
                core::f64::consts::TAU * (incidence + 2.0 * self.properties[s].camber_ratio),
                StripRegime::Attached,
            )
        };
        let expected = solve_nonlinear_with_section_frames(
            &op,
            &self.panels,
            &strips,
            freestream,
            1.225,
            &closure,
            &frames,
            local,
            None,
        )
        .unwrap();
        for (x, y) in gamma.iter().zip(expected.gamma) {
            assert!((x - y).abs() < 2e-8, "gamma {x:e} != {y:e}");
        }
        for (x, y) in alpha.iter().zip(expected.alpha_eff_rad) {
            assert!((x - y).abs() < 2e-9);
        }
    }

    fn momentum(&self, state: &GearedState) -> Vec3 {
        let mut body = self
            .assembly
            .locked_inertia_kgm2()
            .apply(vec(state.rigid.omega_body));
        for rotor in self.assembly.rotors() {
            body = body
                + vec(rotor.axis_body).scale(
                    rotor.axial_inertia_kgm2 * rotor.ratio * self.assembly.speed(state).unwrap(),
                );
        }
        orientation(state).rotate(body).unwrap()
            + (vec(state.rigid.pos_m) - vec(ORIGIN))
                .cross(vec(state.rigid.vel_mps).scale(self.assembly.mass_kg()))
    }

    fn assert_motion(&self, before: &GearedState, step: &GearedStep<2>) {
        assert!(
            (self.momentum(&step.state)
                - self.momentum(before)
                - vec(step.external_angular_impulse_world))
            .norm()
                < 2e-9
        );
        assert!(vec(step.linear_momentum_residual).norm() < 2e-10);
    }
}

#[test]
fn steady_swirl_dynamic_harmonic_and_continuation_share_the_accepted_section_planes() {
    let fixture = Fixture::new();
    let surfaces = LiftingSurfaces::new(
        &fixture.panels,
        &fixture.strips,
        &fixture.properties,
        &fixture.groups,
    )
    .unwrap();
    let model = fixture.model(&surfaces);
    let state = fixture.state();
    let wind = |_| Ok(ZERO);
    let steady = model.evaluate(&state, inputs(), &wind, None).unwrap();
    let fields = &steady.coupled;
    fixture.assert_wing(
        &state,
        &fields.gamma,
        &fields.strip_alpha_eff_rad,
        fields.strip_velocity_mps.as_ref().unwrap(),
    );
    assert!(steady.aerodynamic_force_body_n[1] > 50.0);
    let swirl = SwirlPropulsion::new(&model, fs_airscrew::swirl::Controls::default())
        .evaluate(&state, inputs(), &wind, None)
        .unwrap();
    let fields = &swirl.coupled.coupled;
    fixture.assert_wing(
        &state,
        &fields.gamma,
        &fields.strip_alpha_eff_rad,
        fields.strip_velocity_mps.as_ref().unwrap(),
    );
    let dynamic = DynamicPropulsion::new(&model, 1.225).unwrap();
    let wake = DynamicState {
        motion: state,
        induced_mps: steady.coupled.w_slip,
    };
    let instantaneous = dynamic.evaluate(&wake, inputs(), &wind).unwrap();
    let fields = &instantaneous.aerodynamics;
    fixture.assert_wing(
        &state,
        &fields.gamma,
        &fields.strip_alpha_eff_rad,
        &fields.strip_velocity_mps,
    );
    let harmonic = HarmonicPropulsion::new(&model, 1.225, 16).unwrap();
    let wake = HarmonicState {
        motion: state,
        induction: core::array::from_fn(|k| HarmonicInduction {
            mean_mps: steady.coupled.w_slip[k],
            first_mps: [0.04, -0.03],
        }),
    };
    let instantaneous = harmonic.evaluate(&wake, inputs(), &wind).unwrap();
    let fields = &instantaneous.aerodynamics;
    fixture.assert_wing(
        &state,
        &fields.gamma,
        &fields.strip_alpha_eff_rad,
        &fields.strip_velocity_mps,
    );
    let continued = harmonic
        .initialize_continuation(wake, inputs(), &wind, domain())
        .unwrap();
    let instantaneous = harmonic
        .evaluate_continued(&continued, inputs(), &wind, domain())
        .unwrap();
    let fields = &instantaneous.aerodynamics;
    fixture.assert_wing(
        &state,
        &fields.gamma,
        &fields.strip_alpha_eff_rad,
        &fields.strip_velocity_mps,
    );
    assert!(fields.wing_continuation.is_some());
}

#[test]
fn actual_vertical_surface_midpoint_loads_close_rigid_and_rotor_momentum_in_every_wake_lane() {
    let fixture = Fixture::new();
    let surfaces = LiftingSurfaces::new(
        &fixture.panels,
        &fixture.strips,
        &fixture.properties,
        &fixture.groups,
    )
    .unwrap();
    let model = fixture.model(&surfaces);
    let state = fixture.state();
    let dt = 0.002;
    let steady = model
        .step(
            &state,
            0.0,
            dt,
            ORIGIN,
            None,
            |_, _| Ok(inputs()),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    fixture.assert_motion(&state, &steady.motion);
    assert!(steady.midpoint.aerodynamic_force_body_n[1] > 50.0);
    let swirl = SwirlPropulsion::new(&model, fs_airscrew::swirl::Controls::default())
        .step(
            &state,
            0.0,
            dt,
            ORIGIN,
            None,
            |_, _| Ok(inputs()),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    fixture.assert_motion(&state, &swirl.motion);
    let dynamic = DynamicPropulsion::new(&model, 1.225).unwrap();
    let wake = DynamicState {
        motion: state,
        induced_mps: steady.initial.coupled.w_slip,
    };
    let advanced = dynamic
        .step(&wake, 0.0, dt, ORIGIN, |_, _| Ok(inputs()), |_, _| Ok(ZERO))
        .unwrap();
    fixture.assert_motion(&state, &advanced.motion);
    assert!(advanced.midpoint.aerodynamic_force_body_n[1] > 50.0);
    let harmonic = HarmonicPropulsion::new(&model, 1.225, 16).unwrap();
    let wake = HarmonicState {
        motion: state,
        induction: core::array::from_fn(|k| HarmonicInduction {
            mean_mps: steady.initial.coupled.w_slip[k],
            first_mps: [0.04, -0.03],
        }),
    };
    let continued = harmonic
        .initialize_continuation(wake, inputs(), &|_| Ok(ZERO), domain())
        .unwrap();
    let advanced = harmonic
        .step_continued(
            &continued,
            domain(),
            0.0,
            dt,
            ORIGIN,
            |_, _| Ok(inputs()),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    fixture.assert_motion(&state, &advanced.step.motion);
    assert!(advanced.step.midpoint.aerodynamic_force_body_n[1] > 50.0);
    assert_ne!(
        advanced.step.initial.aerodynamics.gamma,
        advanced.step.midpoint.aerodynamics.gamma
    );
}
