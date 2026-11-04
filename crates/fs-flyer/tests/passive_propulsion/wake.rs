//! Physical annular wake transport coupled to actual independent shaft states.
use super::*;
use fs_flyer::aircraft::OpenLoopDesign;
use fs_flyer::propcoupling::swirl::wake::WakeTransport;
use fs_flyer::propcoupling::{PropCouplingSolverSpec, swirl::wake::WakeControls};
use fs_flyer::spine::rotors::propulsion::passive::wake::WakeShaftPropulsion;
use fs_flyer::spine::rotors::propulsion::surfaces::{
    InducedDragGroup, LiftingSurfaces, SectionProperties,
};
use fs_flyer::spine::rotors::propulsion::swirl::wake::WakePropulsion;
use fs_wing::nonlinear::StripSpec;
use fs_wing::{Panel, SurfaceId};

const TRANSPORTS: [WakeTransport; 2] = [
    WakeTransport::FixedCylinder,
    WakeTransport::MaterialStreamtubes,
];

fn fixture() -> (
    OpenLoopDesign,
    Vec<Panel>,
    Vec<StripSpec>,
    Vec<SectionProperties>,
    Vec<InducedDragGroup>,
) {
    let mut d = wright_openloop_v1();
    d.coupling = PropCouplingSolverSpec {
        tol: 1e-10,
        cap: 64,
        ..CANDIDATE_B
    };
    d.camber_ratio = 0.0;
    d.disks[0].center_m = [0.0, -2.5, 0.0];
    d.disks[1].center_m = [0.0, 2.5, 0.0];
    let panels = d
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
    let sections = vec![
        SectionProperties {
            camber_ratio: 0.0,
            attached_cd: 0.01,
            canard_gain: 0.0,
            induced_drag_group: 0,
        };
        2
    ];
    let groups = vec![InducedDragGroup {
        span_m: 5.3,
        oswald_e: 0.85,
    }];
    (d, panels, strips, sections, groups)
}
fn start(a: &RotorAssembly<2>, speed: f64, rates: [f64; 2]) -> RotorState<2> {
    a.state_from_relative_spin(
        SixDofState {
            pos_m: [1.0, -2.0, -10.0],
            vel_mps: [speed, 0.0, 0.0],
            quat: [1.0, 0.0, 0.0, 0.0],
            omega_body: ZERO,
        },
        [0.2, -0.1],
        rates,
    )
    .unwrap()
}
fn control() -> ShaftInputs<2> {
    ShaftInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.0,
        motor_nm: [40.0, -55.0],
        gravity_world_mps2: ZERO,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
}
fn kinetic(a: &RotorAssembly<2>, s: &RotorState<2>) -> f64 {
    let omega = v(s.rigid.omega_body);
    0.5 * a.reduced_body.mass_kg * v(s.rigid.vel_mps).dot(v(s.rigid.vel_mps))
        + 0.5 * omega.dot(a.reduced_body.inertia_tensor().matrix().apply(omega))
        + (0..2)
            .map(|i| s.axial_momentum[i].powi(2) / (2.0 * a.rotors[i].axial_inertia_kgm2))
            .sum::<f64>()
}

#[test]
fn g3_annular_fields_match_geared_snapshots_with_physical_axis_signs() {
    let (d, panels, strips, sections, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    for sign in [1.0, -1.0] {
        let axes = [[sign, 0.0, 0.0]; 2];
        let a = assembly(d.gross_mass_kg, axes);
        let mut m = mounts(AXIAL);
        if sign < 0.0 {
            m[0].handedness = PropellerHandedness::Negative;
            m[1].handedness = PropellerHandedness::Positive;
        }
        let mut locked = a.reduced_body.inertia_tensor().matrix();
        for rotor in a.rotors {
            for i in 0..3 {
                for j in 0..3 {
                    locked.m[3 * i + j] +=
                        rotor.axial_inertia_kgm2 * rotor.axis_body[i] * rotor.axis_body[j];
                }
            }
        }
        let gears = GearedAssembly::new(
            d.gross_mass_kg,
            locked,
            core::array::from_fn(|i| GearedRotor {
                axis_body: axes[i],
                axial_inertia_kgm2: a.rotors[i].axial_inertia_kgm2,
                ratio: [sign, -sign][i],
            }),
        )
        .unwrap();
        let passive = ShaftPropulsion::new_spatial(&d, &a, m)
            .unwrap()
            .with_lifting_surfaces(&surfaces)
            .unwrap();
        let geared = GearedPropulsion::new_spatial(&d, &gears, m)
            .unwrap()
            .with_lifting_surfaces(&surfaces)
            .unwrap();
        let state = start(&a, 12.0, [35.0 * sign, -35.0 * sign]);
        let geared_state = gears.state_from_speed(state.rigid, 0.2, 35.0).unwrap();
        let axial = SwirlShaftPropulsion::new(&passive, Controls::default())
            .evaluate(&state, control(), &|_| Ok(ZERO), None)
            .unwrap();
        for transport in TRANSPORTS {
            let model =
                WakeShaftPropulsion::new(&passive, Controls::default(), WakeControls::default())
                    .with_transport(transport);
            let actual = model
                .evaluate(&state, control(), &|_| Ok(ZERO), None)
                .unwrap();
            let expected =
                WakePropulsion::new(&geared, Controls::default(), WakeControls::default())
                    .with_transport(transport)
                    .evaluate(
                        &geared_state,
                        PropulsionInputs {
                            density_kg_m3: control().density_kg_m3,
                            canard_rad: 0.0,
                            motor_nm: 0.0,
                            gravity_world_mps2: ZERO,
                        },
                        &|_| Ok(ZERO),
                        None,
                    )
                    .unwrap();
            assert_eq!(actual.coupled, expected.coupled);
            assert_eq!(actual.carrier, expected.loads.carrier);
            assert_eq!(actual.torques.aerodynamic_nm, expected.loads.aerodynamic_nm);
            assert_eq!(
                actual.propeller_absolute_speed_rad_s,
                [35.0 * sign, -35.0 * sign]
            );
            assert!((v(actual.carrier.moment_nm) - v(axial.carrier.moment_nm)).norm() > 0.1);
            for k in 0..2 {
                let flux = actual.coupled.wakes[k].cross_section(0.8).unwrap();
                near(
                    flux.angular_momentum_flux_nm,
                    -sign * actual.torques.aerodynamic_nm[k],
                    2e-13,
                );
                assert!(actual.coupled.strip_velocity_correction_mps()[k][2].abs() > 0.01);
            }
        }
    }
}

#[test]
fn g1_unequal_shafts_transport_actual_torque_and_preserve_mechanical_work_and_momentum() {
    let (d, panels, strips, sections, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new_spatial(&d, &a, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    // Sustained free rotation with bearing drag; no spring to a fixed body phase.
    let laws = ShaftProperties {
        damping_nms_per_rad: [0.12, 0.18],
        ..ShaftProperties::default()
    };
    let dt = 0.004;
    let origin = [0.3, -0.4, 0.5];
    for transport in TRANSPORTS {
        let model = WakeShaftPropulsion::new(&base, Controls::default(), WakeControls::default())
            .with_transport(transport);
        let mut motion = a
            .shaft_motion(&start(&a, 12.0, [32.0, -39.0]), laws, dt)
            .unwrap();
        motion
            .apply_impulse(
                RotorImpulses {
                    carrier: CarrierImpulse::default(),
                    motor_nms: [0.25, -0.125],
                    aerodynamic_nms: [-0.0625, 0.03125],
                },
                origin,
            )
            .unwrap();
        let before = *motion.state();
        let mut stages = Vec::new();
        let result = model
            .step(
                &mut motion,
                1.0,
                origin,
                None,
                |time, state| {
                    stages.push((time, *state));
                    Ok(control())
                },
                |_, _| Ok(ZERO),
            )
            .unwrap();
        assert_eq!(stages.len(), 2);
        assert_eq!(stages[0], (1.0, before));
        assert_eq!(stages[1].0, 1.0 + 0.5 * dt);
        assert_ne!(stages[1].1.axial_momentum, before.axial_momentum);
        assert_eq!(
            result.midpoint,
            model
                .evaluate(
                    &stages[1].1,
                    control(),
                    &|_| Ok(ZERO),
                    Some(&result.initial.coupled.warm_start())
                )
                .unwrap()
        );
        assert_ne!(result.initial.coupled.wakes, result.midpoint.coupled.wakes);
        assert!(
            (result.midpoint.propeller_absolute_speed_rad_s[0]
                + result.midpoint.propeller_absolute_speed_rad_s[1])
                .abs()
                > 1.0
        );
        for stage in [&result.initial, &result.midpoint] {
            for k in 0..2 {
                let source = &stage.coupled.solution.rotors[k];
                let wake = &stage.coupled.wakes[k];
                near(
                    source.power.shaft_w,
                    -stage.torques.aerodynamic_nm[k] * stage.propeller_absolute_speed_rad_s[k],
                    1e-11,
                );
                let reference = wake.cross_section(0.1 * d.rotor.radius_m).unwrap();
                for distance in [0.8, 1.8].map(|s| s * d.rotor.radius_m) {
                    let section = wake.cross_section(distance).unwrap();
                    near(
                        section.angular_momentum_flux_nm,
                        -stage.torques.aerodynamic_nm[k],
                        2e-13,
                    );
                    if transport == WakeTransport::MaterialStreamtubes {
                        near(section.mass_flux_kg_s, reference.mass_flux_kg_s, 1e-12);
                    }
                }
            }
        }
        let next = result.motion.state;
        let (p0, l0) = momentum(&a, &before, origin);
        let (p1, l1) = momentum(&a, &next, origin);
        let impulse = v(result.midpoint.carrier.force_n).scale(dt);
        let mut torque = v(result.midpoint.carrier.moment_nm);
        for i in 0..2 {
            torque =
                torque + v(a.rotors[i].axis_body).scale(result.midpoint.torques.aerodynamic_nm[i]);
        }
        let angular = (rotation(&before).rotate(torque).unwrap()
            + rotation(&next).rotate(torque).unwrap())
        .scale(0.5 * dt)
            + (v(before.rigid.pos_m) + v(before.rigid.vel_mps).scale(0.5 * dt) - v(origin))
                .cross(impulse);
        assert!((p1 - p0 - impulse).norm() < 1e-10);
        assert!((l1 - l0 - angular).norm() < 1e-9);
        let motor_work: f64 = (0..2)
            .map(|i| {
                control().motor_nm[i] * (next.relative_phase_rad[i] - before.relative_phase_rad[i])
            })
            .sum();
        near(result.motion.motor_work_j, motor_work, 1e-11);
        let heat: f64 = result.motion.native.dissipated_energy_j.into_iter().sum();
        assert!(heat > 0.0);
        let balance = kinetic(&a, &next) - kinetic(&a, &before) + heat
            - motor_work
            - result.motion.external_work_j
            - result.motion.native.motion.drift_energy_defect_j;
        assert!(balance.abs() < 1e-8, "physical work balance {balance}");
        assert_eq!(motion.state(), &next);
        assert_eq!(
            motion.native_state(),
            &result.motion.native.motion.state_after
        );
        eprintln!(
            "{transport:?}: independent speeds {:?}, momentum residual {:.3e}, work residual {balance:.3e}",
            result.midpoint.propeller_absolute_speed_rad_s,
            (l1 - l0 - angular).norm()
        );
    }
}

#[test]
fn g1_windmilling_reverses_source_transport_and_accelerates_each_free_shaft() {
    let (mut d, panels, strips, sections, groups) = fixture();
    d.rotor = windmilling_rotor();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new_spatial(&d, &a, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    for transport in TRANSPORTS {
        let model = WakeShaftPropulsion::new(&base, Controls::default(), WakeControls::default())
            .with_transport(transport);
        let mut motion = a
            .shaft_motion(
                &start(&a, 10.0, [60.0, -60.0]),
                ShaftProperties::default(),
                0.001,
            )
            .unwrap();
        let before = *motion.state();
        let result = model
            .step(
                &mut motion,
                0.0,
                before.rigid.pos_m,
                None,
                |_, _| {
                    Ok(ShaftInputs {
                        motor_nm: [0.0; 2],
                        ..control()
                    })
                },
                |_, _| Ok(ZERO),
            )
            .unwrap();
        for stage in [&result.initial, &result.midpoint] {
            for k in 0..2 {
                let source = &stage.coupled.solution.rotors[k];
                assert!(source.torque_nm < 0.0 && source.power.shaft_w < 0.0);
                let section = stage.coupled.wakes[k].cross_section(0.8).unwrap();
                near(
                    section.angular_momentum_flux_nm,
                    -stage.torques.aerodynamic_nm[k],
                    2e-13,
                );
                assert!(section.angular_momentum_flux_nm * [1.0, -1.0][k] < 0.0);
            }
        }
        for k in 0..2 {
            let delta = result.motion.state.axial_momentum[k] - before.axial_momentum[k];
            near(
                delta,
                0.001 * result.midpoint.torques.aerodynamic_nm[k],
                2e-13,
            );
            assert!(result.motion.state.axial_momentum[k].abs() > before.axial_momentum[k].abs());
        }
        assert_eq!(result.motion.motor_work_j, 0.0);
        assert_eq!(result.motion.native.dissipated_energy_j, [0.0; 2]);
        assert!(motion.state().rigid.vel_mps[0] < before.rigid.vel_mps[0]);
    }
}

#[test]
fn g0_field_budgets_midpoint_and_late_native_refusals_preserve_state_and_retry() {
    let (d, panels, strips, sections, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new_spatial(&d, &a, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let start = start(&a, 12.0, [32.0, -39.0]);
    for transport in TRANSPORTS {
        let model = WakeShaftPropulsion::new(&base, Controls::default(), WakeControls::default())
            .with_transport(transport);
        let mut motion = a
            .shaft_motion(&start, ShaftProperties::default(), 0.001)
            .unwrap();
        let before = *motion.state();
        let canonical = *motion.native_state();
        for (controls, code) in [
            (
                WakeControls {
                    relative_power_tolerance: 0.0,
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
        ] {
            let bounded = WakeShaftPropulsion::new(&base, Controls::default(), controls)
                .with_transport(transport);
            let error = bounded
                .step(
                    &mut motion,
                    0.0,
                    ZERO,
                    None,
                    |_, _| Ok(control()),
                    |_, _| Ok(ZERO),
                )
                .unwrap_err();
            assert_eq!(error.code, code);
            assert_eq!(motion.state(), &before);
            assert_eq!(motion.native_state(), &canonical);
        }
        let seed = model
            .evaluate(&before, control(), &|_| Ok(ZERO), None)
            .unwrap()
            .coupled
            .warm_start();
        let saved_seed = seed.clone();
        for late_native in [false, true] {
            let mut calls = 0;
            let error = model
                .step(
                    &mut motion,
                    0.0,
                    ZERO,
                    Some(&seed),
                    |_, _| {
                        calls += 1;
                        Ok(ShaftInputs {
                            motor_nm: if late_native && calls == 2 {
                                [f64::MAX / 32.0; 2]
                            } else {
                                control().motor_nm
                            },
                            ..control()
                        })
                    },
                    |time, _| {
                        if !late_native && time > 0.0 {
                            Err(fs_flyer::Refusal {
                                code: "wind-domain",
                                message: "midpoint outside material wind domain".into(),
                                ranked_repairs: vec![],
                            })
                        } else {
                            Ok(ZERO)
                        }
                    },
                )
                .unwrap_err();
            assert_eq!(calls, 2);
            assert_eq!(
                error.code,
                if late_native {
                    "rotor-shaft-refused"
                } else {
                    "wind-domain"
                }
            );
            assert_eq!(motion.state(), &before);
            assert_eq!(motion.native_state(), &canonical);
            assert_eq!(seed, saved_seed);
        }
        let mut fresh = a
            .shaft_motion(&start, ShaftProperties::default(), 0.001)
            .unwrap();
        let retry = model
            .step(
                &mut motion,
                0.0,
                ZERO,
                Some(&seed),
                |_, _| Ok(control()),
                |_, _| Ok(ZERO),
            )
            .unwrap();
        let expected = model
            .step(
                &mut fresh,
                0.0,
                ZERO,
                Some(&seed),
                |_, _| Ok(control()),
                |_, _| Ok(ZERO),
            )
            .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(motion.native_state(), fresh.native_state());
    }
}
