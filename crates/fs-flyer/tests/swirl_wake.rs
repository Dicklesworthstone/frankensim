//! Real wing/swirl loads, signed shaft ports and complete body/rotor momentum.
use fs_airscrew::swirl::Controls;
use fs_airscrew::{BladeStation, Rotor};
use fs_flyer::aircraft::{OpenLoopDesign, wright_openloop_v1};
use fs_flyer::propcoupling::{CANDIDATE_B, PropCouplingSolverSpec};
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor, GearedState};
use fs_flyer::spine::rotors::propulsion::swirl::SwirlPropulsion;
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};

const ZERO: [f64; 3] = [0.0; 3];
const AXES: [[f64; 3]; 2] = [[1.0, 0.0, 0.0]; 2];
fn vec(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn original_design() -> OpenLoopDesign {
    let mut design = wright_openloop_v1();
    design.coupling = PropCouplingSolverSpec {
        tol: 1e-10,
        cap: 64,
        ..CANDIDATE_B
    };
    design
}
fn assembly(mass: f64, axes: [[f64; 3]; 2], ratios: [f64; 2]) -> GearedAssembly<2> {
    GearedAssembly::new(
        mass,
        Mat3 {
            m: [1000.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0, 0.0, 1200.0],
        },
        core::array::from_fn(|k| GearedRotor {
            axis_body: axes[k],
            ratio: ratios[k],
            axial_inertia_kgm2: [1.5, 2.5][k],
        }),
    )
    .unwrap()
}
fn mounts(axes: [[f64; 3]; 2], ratios: [f64; 2]) -> [PropellerMount; 2] {
    core::array::from_fn(|k| PropellerMount {
        rotor_index: k,
        thrust_direction_body: axes[k],
        handedness: if ratios[k] > 0.0 {
            PropellerHandedness::Positive
        } else {
            PropellerHandedness::Negative
        },
    })
}
fn state(assembly: &GearedAssembly<2>, speed: f64, alpha: f64, omega: f64) -> GearedState {
    assembly
        .state_from_speed(
            SixDofState {
                pos_m: [2.0, -3.0, -10.0],
                vel_mps: [speed * alpha.cos(), 0.0, speed * alpha.sin()],
                quat: [1.0, 0.0, 0.0, 0.0],
                omega_body: ZERO,
            },
            0.0,
            omega,
        )
        .unwrap()
}
fn inputs() -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.03,
        motor_nm: 80.0,
        gravity_world_mps2: ZERO,
    }
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
fn angular_momentum(a: &GearedAssembly<2>, state: &GearedState, origin: [f64; 3]) -> Vec3 {
    let mut spin = a.locked_inertia_kgm2().apply(vec(state.rigid.omega_body));
    let relative = a.speed(state).unwrap();
    for r in a.rotors() {
        spin = spin + vec(r.axis_body).scale(r.axial_inertia_kgm2 * r.ratio * relative);
    }
    orientation(state).rotate(spin).unwrap()
        + (vec(state.rigid.pos_m) - vec(origin)).cross(vec(state.rigid.vel_mps).scale(a.mass_kg()))
}

use fs_flyer::propcoupling::swirl::wake::WakeControls;
use fs_flyer::spine::rotors::propulsion::surfaces::{
    InducedDragGroup, LiftingSurfaces, SectionProperties,
};
use fs_flyer::spine::rotors::propulsion::swirl::wake::WakePropulsion;
use fs_wing::nonlinear::StripSpec;
use fs_wing::{Panel, SurfaceId};

fn fixture() -> (
    OpenLoopDesign,
    Vec<Panel>,
    Vec<StripSpec>,
    Vec<SectionProperties>,
    Vec<InducedDragGroup>,
) {
    let mut d = original_design();
    d.camber_ratio = 0.0;
    d.disks[0].center_m = [0.0, -2.5, 0.0];
    d.disks[1].center_m = [0.0, 2.5, 0.0];
    let panels: Vec<_> = d
        .disks
        .iter()
        .map(|disk| {
            let y = disk.center_m[1] + 0.7;
            Panel {
                surface: SurfaceId::WingLower,
                a: [-0.8, y - 0.15, 0.0],
                b: [-0.8, y + 0.15, 0.0],
                ctrl: [-0.9, y, 0.0],
                normal: [0.0, 0.0, -1.0],
                width_m: 0.3,
            }
        })
        .collect();
    let strips = (0..2)
        .map(|k| StripSpec {
            panel_indices: vec![k],
            chord_m: 0.2,
            twist_rad: 0.0,
        })
        .collect();
    let properties = vec![
        SectionProperties {
            camber_ratio: 0.0,
            attached_cd: 0.01,
            canard_gain: 0.0,
            induced_drag_group: 0
        };
        2
    ];
    let groups = vec![InducedDragGroup {
        span_m: 5.3,
        oswald_e: 0.85,
    }];
    (d, panels, strips, properties, groups)
}

#[test]
fn downstream_handedness_changes_actual_wing_force_and_retains_used_field() {
    let (d, panels, strips, properties, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let mut forces = Vec::new();
    for sign in [1.0, -1.0] {
        let a = assembly(d.gross_mass_kg, AXES, [sign; 2]);
        let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [sign; 2]))
            .unwrap()
            .with_lifting_surfaces(&surfaces)
            .unwrap();
        let st = state(&a, 12.0, 0.0, 35.0);
        let native = SwirlPropulsion::new(&base, Controls::default())
            .evaluate(&st, inputs(), &|_| Ok(ZERO), None)
            .unwrap();
        let model = WakePropulsion::new(&base, Controls::default(), WakeControls::default());
        let result = model.evaluate(&st, inputs(), &|_| Ok(ZERO), None).unwrap();
        let flow = &result.coupled;
        assert!(flow.field_evaluations > 1);
        assert!(flow.field_residual_mps <= WakeControls::default().velocity_tolerance_mps);
        assert!(
            (result.aerodynamic_force_body_n[2] - native.aerodynamic_force_body_n[2]).abs() > 0.1
        );
        assert!(result.aerodynamic_force_body_n[2] * sign > 0.1);
        for k in 0..2 {
            let point = array((vec(panels[k].a) + vec(panels[k].b)).scale(0.5));
            let predicted = vec(flow.wakes[0].relative_velocity_at(point).unwrap())
                + vec(flow.wakes[1].relative_velocity_at(point).unwrap());
            assert!(
                (predicted - vec(flow.used_strip_swirl_velocity_mps[k])).norm()
                    <= flow.field_residual_mps + 1e-14
            );
            let actual = flow.solution.coupled.strip_velocity_mps.as_ref().unwrap()[k];
            near(actual[2], flow.used_strip_swirl_velocity_mps[k][2], 1e-14);
            assert_eq!(
                result.loads.aerodynamic_nm[k],
                -sign * flow.solution.rotors[k].torque_nm
            );
        }
        let continued = model
            .evaluate(&st, inputs(), &|_| Ok(ZERO), Some(&flow.warm_start()))
            .unwrap();
        assert!(continued.coupled.field_evaluations <= 2);
        near(
            continued.aerodynamic_force_body_n[2],
            result.aerodynamic_force_body_n[2],
            2e-6,
        );
        forces.push(result.aerodynamic_force_body_n);
    }
    near(forces[0][2], -forces[1][2], 2e-8);
}

#[test]
fn each_native_annulus_transports_actual_torque_and_resolved_energy_at_every_section() {
    let (d, panels, strips, properties, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXES, [1.0, -1.0]);
    let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0, -1.0]))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let result = WakePropulsion::new(&base, Controls::default(), WakeControls::default())
        .evaluate(&state(&a, 12.0, 0.0, 35.0), inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    for k in 0..2 {
        let wake = &result.coupled.wakes[k];
        let source = &result.coupled.solution.rotors[k];
        let sign = [1.0, -1.0][k];
        for distance in [
            0.05 * d.rotor.radius_m,
            0.8 * d.rotor.radius_m,
            1.99 * d.rotor.radius_m,
        ] {
            let summary = wake.cross_section(distance).unwrap();
            let induction = result.coupled.solution.coupled.w_slip[k];
            let axial = 12.0
                + (1.0
                    + distance
                        / (distance * distance + d.rotor.radius_m * d.rotor.radius_m).sqrt())
                    * induction;
            near(summary.axial_speed_mps, axial, 1e-12);
            let mut total_energy = 0.0;
            for (cell, station) in wake.annuli().iter().zip(&source.stations) {
                let (lo, hi) = (cell.inner_radius_m, cell.outer_radius_m);
                let mut momentum = 0.0;
                let mut energy = 0.0;
                // Independent polar area quadrature of the PHYSICAL velocity,
                // not the reported flux or a re-normalized momentum residual.
                let nr = 256;
                let nt = 16;
                let area = std::f64::consts::PI * (hi * hi - lo * lo) / (nr * nt) as f64;
                for i in 0..nr {
                    let radius =
                        (lo * lo + (hi * hi - lo * lo) * (i as f64 + 0.5) / nr as f64).sqrt();
                    for j in 0..nt {
                        let angle = std::f64::consts::TAU * (j as f64 + 0.5) / nt as f64;
                        let radial = Vec3::new(0.0, radius * angle.cos(), radius * angle.sin());
                        let point =
                            vec(d.disks[k].center_m) + Vec3::new(-distance, 0.0, 0.0) + radial;
                        let fluid =
                            vec(wake.relative_velocity_at(array(point)).unwrap()).scale(-1.0);
                        momentum += 1.225 * axial * radial.cross(fluid).x * area;
                        energy += 0.5 * 1.225 * axial * fluid.dot(fluid) * area;
                    }
                }
                near(momentum, sign * station.torque_nm, 2e-11);
                total_energy += energy;
            }
            near(total_energy, summary.resolved_swirl_power_w, 2e-6);
            near(
                summary.angular_momentum_flux_nm,
                sign * source.torque_nm,
                1e-13,
            );
            near(summary.native_swirl_power_w, source.power.swirl_w, 1e-13);
            near(
                summary.power_defect_w,
                total_energy - source.power.swirl_w,
                2e-6,
            );
            assert!(summary.absolute_cell_power_defect_w <= wake.maximum_power_defect_w() + 1e-12);
        }
        let near_flux = wake.cross_section(0.0).unwrap();
        let far_flux = wake.cross_section(2.0 * d.rotor.radius_m).unwrap();
        near(
            far_flux.axial_speed_mps,
            12.0 + (1.0 + 2.0 / 5.0_f64.sqrt()) * result.coupled.solution.coupled.w_slip[k],
            1e-12,
        );
        assert!(far_flux.mass_flux_kg_s > near_flux.mass_flux_kg_s);
        assert!(wake.maximum_power_defect_w() <= 0.5 * source.power.swirl_w);
        let radius = 0.7;
        assert_eq!(
            wake.relative_velocity_at([0.1, d.disks[k].center_m[1] + radius, 0.0])
                .unwrap(),
            ZERO
        );
        assert_eq!(
            wake.relative_velocity_at([-0.8, d.disks[k].center_m[1] + 2.0 * d.rotor.radius_m, 0.0])
                .unwrap(),
            ZERO
        );
        assert_eq!(
            wake.relative_velocity_at([
                -3.0 * d.rotor.radius_m,
                d.disks[k].center_m[1] + radius,
                0.0
            ])
            .unwrap_err()
            .code,
            "swirl-wake-outside-domain"
        );
    }
}

#[test]
fn actual_wake_midpoint_loads_preserve_global_momentum_and_world_covariance() {
    let (d, panels, strips, properties, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXES, [1.0, 1.0]);
    let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0, 1.0]))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let model = WakePropulsion::new(&base, Controls::default(), WakeControls::default());
    let start = state(&a, 12.0, 0.0, 35.0);
    let origin = [-1.0, 2.0, 3.0];
    let result = model
        .step(
            &start,
            0.0,
            0.001,
            origin,
            None,
            |_, _| {
                Ok(PropulsionInputs {
                    motor_nm: 0.0,
                    ..inputs()
                })
            },
            |_, _| Ok(ZERO),
        )
        .unwrap();
    let dp = (vec(result.motion.state.rigid.vel_mps) - vec(start.rigid.vel_mps)).scale(a.mass_kg());
    let dl =
        angular_momentum(&a, &result.motion.state, origin) - angular_momentum(&a, &start, origin);
    assert!((dp - vec(result.motion.external_linear_impulse_world)).norm() < 2e-10);
    assert!((dl - vec(result.motion.external_angular_impulse_world)).norm() < 2e-9);
    assert!(result.midpoint.loads.carrier.force_n[2] > 0.1);
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.1))).unwrap();
    let mut rotated = start;
    rotated.rigid.pos_m = array(rotation.rotate(vec(start.rigid.pos_m)).unwrap());
    rotated.rigid.vel_mps = array(rotation.rotate(vec(start.rigid.vel_mps)).unwrap());
    let q = rotation.as_quat();
    rotated.rigid.quat = [q.w, q.x, q.y, q.z];
    let evaluation = model
        .evaluate(&rotated, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    assert!(
        (vec(evaluation.aerodynamic_force_body_n) - vec(result.initial.aerodynamic_force_body_n))
            .norm()
            < 1e-8
    );
    assert!(
        (vec(evaluation.loads.carrier.force_n)
            - rotation
                .rotate(vec(result.initial.loads.carrier.force_n))
                .unwrap())
        .norm()
            < 1e-8
    );
}

#[test]
fn physical_mount_axis_and_handedness_determine_swirl_instead_of_gear_sign() {
    let (d, panels, strips, properties, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let negative_axes = [[-1.0, 0.0, 0.0]; 2];
    let a = assembly(d.gross_mass_kg, negative_axes, [-1.0; 2]);
    let mounts = core::array::from_fn(|k| PropellerMount {
        rotor_index: k,
        handedness: PropellerHandedness::Negative,
        thrust_direction_body: AXES[k],
    });
    let base = GearedPropulsion::new_spatial(&d, &a, mounts)
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let result = WakePropulsion::new(&base, Controls::default(), WakeControls::default())
        .evaluate(&state(&a, 12.0, 0.0, 35.0), inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    assert!(result.aerodynamic_force_body_n[2] > 0.1);
    assert!(
        result
            .propeller_absolute_speed_rad_s
            .iter()
            .all(|w| *w < 0.0)
    );
    for k in 0..2 {
        assert!(
            result.coupled.wakes[k]
                .cross_section(0.8)
                .unwrap()
                .angular_momentum_flux_nm
                > 0.0
        );
        assert!(result.loads.aerodynamic_nm[k] > 0.0);
    }
}

#[test]
fn unsupported_disk_swirl_and_unbounded_energy_or_iteration_refuse() {
    let (mut d, panels, strips, properties, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXES, [1.0; 2]);
    let start = state(&a, 12.0, 0.0, 35.0);
    {
        let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0; 2]))
            .unwrap()
            .with_lifting_surfaces(&surfaces)
            .unwrap();
        for (control, code) in [
            (
                WakeControls {
                    relative_power_tolerance: 0.01,
                    ..WakeControls::default()
                },
                "swirl-wake-power-budget",
            ),
            (
                WakeControls {
                    maximum_field_evaluations: 1,
                    ..WakeControls::default()
                },
                "swirl-wake-nonconvergence",
            ),
            (
                WakeControls {
                    relaxation: f64::NAN,
                    ..WakeControls::default()
                },
                "swirl-wake-invalid",
            ),
        ] {
            assert_eq!(
                WakePropulsion::new(&base, Controls::default(), control)
                    .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
                    .unwrap_err()
                    .code,
                code
            );
        }
        let original = start;
        let model = WakePropulsion::new(&base, Controls::default(), WakeControls::default());
        let mut calls = 0;
        let result = model.step(
            &start,
            0.0,
            0.001,
            ZERO,
            None,
            |_, _| {
                calls += 1;
                Ok(PropulsionInputs {
                    motor_nm: if calls == 2 { f64::NAN } else { 0.0 },
                    ..inputs()
                })
            },
            |_, _| Ok(ZERO),
        );
        assert_eq!(result.unwrap_err().code, "propulsion-state-invalid");
        assert_eq!(calls, 2);
        assert_eq!(start, original);
    }
    // Neither a downstream nor a coincident coplanar receiver can consume
    // incident angular momentum through its zero centerline swirl value.
    for center in [[-0.8, -2.5, 0.0], [0.0, -2.5, 0.0]] {
        d.disks[1].center_m = center;
        let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0; 2]))
            .unwrap()
            .with_lifting_surfaces(&surfaces)
            .unwrap();
        assert_eq!(
            WakePropulsion::new(&base, Controls::default(), WakeControls::default())
                .evaluate(&start, inputs(), &|_| Ok(ZERO), None)
                .unwrap_err()
                .code,
            "swirl-wake-disk-interaction"
        );
    }
    // The disks themselves are separated, but their canted cylinders meet at
    // a downstream wing. Separate torque/energy budgets cannot admit mixing.
    let radius = d.rotor.radius_m;
    d.disks[0].center_m = [0.0, -1.5 * radius, 0.0];
    d.disks[1].center_m = [0.0, 1.5 * radius, 0.0];
    let c = 0.5_f64.sqrt();
    let axes = [[c, -c, 0.0], [c, c, 0.0]];
    let crossed_a = assembly(d.gross_mass_kg, axes, [1.0; 2]);
    let crossed_panels = vec![Panel {
        surface: SurfaceId::WingLower,
        a: [-radius, -0.05, 0.0],
        b: [-radius, 0.05, 0.0],
        ctrl: [-radius - 0.1, 0.0, 0.0],
        normal: [0.0, 0.0, -1.0],
        width_m: 0.1,
    }];
    let crossed_strips = vec![StripSpec {
        panel_indices: vec![0],
        chord_m: 0.2,
        twist_rad: 0.0,
    }];
    let crossed_properties = vec![properties[0]];
    let crossed_surfaces = LiftingSurfaces::new(
        &crossed_panels,
        &crossed_strips,
        &crossed_properties,
        &groups,
    )
    .unwrap();
    let crossed = GearedPropulsion::new_spatial(&d, &crossed_a, mounts(axes, [1.0; 2]))
        .unwrap()
        .with_lifting_surfaces(&crossed_surfaces)
        .unwrap();
    let control = WakeControls {
        relative_power_tolerance: 1.0,
        ..WakeControls::default()
    };
    assert_eq!(
        WakePropulsion::new(&crossed, Controls::default(), control)
            .evaluate(
                &state(&crossed_a, 12.0, 0.0, 35.0),
                inputs(),
                &|_| Ok(ZERO),
                None
            )
            .unwrap_err()
            .code,
        "swirl-wake-interaction"
    );
}

/// An attached-flow windmilling rotor with prescribed w=-2 m/s, u=-0.15 m/s
/// at V=10 m/s and Omega=60 rad/s, constructed from both momentum laws.
fn windmill_rotor() -> Rotor {
    let stations = [0.3, 0.5, 0.7, 0.9]
        .into_iter()
        .map(|r| {
            let va = 8.0_f64;
            let vt = 60.0 * r + 0.15;
            let phi = va.atan2(vt);
            let ratio = (phi.cos() + 0.075 * phi.sin()) / (0.075 * phi.cos() - phi.sin());
            let b = std::f64::consts::TAU;
            let alpha = 0.024 * ratio / (b + (b * b - 4.0 * 0.8 * 0.012 * ratio * ratio).sqrt());
            let loss = |distance: f64| {
                2.0 / std::f64::consts::PI * (-1.5 * distance / (r * phi.sin())).exp().acos()
            };
            let f = (loss(1.0 - r) * loss(r - 0.24)).clamp(1e-3, 1.0);
            let target = 4.0 * std::f64::consts::PI * 1.225 * r * f * va * (-2.0);
            let unit_force = 1.5
                * 1.225
                * (va * va + vt * vt)
                * (b * alpha * phi.cos() - (0.012 + 0.8 * alpha * alpha) * phi.sin());
            BladeStation {
                r_over_r: r,
                chord_m: target / unit_force,
                beta_rad: phi + alpha,
            }
        })
        .collect();
    Rotor {
        radius_m: 1.0,
        n_blades: 3,
        camber_ratio: 0.0,
        stations,
    }
}

#[test]
fn windmilling_reverses_angular_transport_and_accelerates_the_real_shaft() {
    let (mut d, panels, strips, properties, groups) = fixture();
    d.rotor = windmill_rotor();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &properties, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXES, [1.0, -1.0]);
    let base = GearedPropulsion::new_spatial(&d, &a, mounts(AXES, [1.0, -1.0]))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let model = WakePropulsion::new(&base, Controls::default(), WakeControls::default());
    let start = state(&a, 10.0, 0.0, 60.0);
    let result = model
        .step(
            &start,
            0.0,
            0.001,
            ZERO,
            None,
            |_, _| {
                Ok(PropulsionInputs {
                    motor_nm: 0.0,
                    canard_rad: 0.0,
                    ..inputs()
                })
            },
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert!(result.motion.speed_rad_s > a.speed(&start).unwrap());
    for stage in [&result.initial, &result.midpoint] {
        for k in 0..2 {
            let source = &stage.coupled.solution.rotors[k];
            let flux = stage.coupled.wakes[k].cross_section(0.8).unwrap();
            assert!(source.torque_nm < 0.0 && source.power.shaft_w < 0.0);
            near(
                flux.angular_momentum_flux_nm,
                [1.0, -1.0][k] * source.torque_nm,
                1e-13,
            );
            assert!(flux.resolved_swirl_power_w > 0.0 && flux.native_swirl_power_w > 0.0);
            let point = [-0.8, d.disks[k].center_m[1] + 0.7, 0.0];
            let relative = stage.coupled.wakes[k].relative_velocity_at(point).unwrap();
            assert!(relative[2] * [1.0, -1.0][k] > 0.0);
        }
    }
}

#[path = "swirl_wake/streamtubes.rs"]
mod streamtubes;
