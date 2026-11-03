//! Real world wind/wing variation at the native blade points and wrench ports.
use super::*;
use fs_flyer::Refusal;
use fs_flyer::spine::rotors::propulsion::dynamic::BladeIncident;
use std::cell::Cell;

const SAMPLES: usize = 32;

fn posed(assembly: &GearedAssembly<2>, rotated: bool) -> (HarmonicState, So3) {
    let angle: f64 = if rotated { 0.34 } else { 0.0 };
    let axis = unit([0.3, -0.4, 0.5]);
    let q = [
        (0.5 * angle).cos(),
        axis[0] * (0.5 * angle).sin(),
        axis[1] * (0.5 * angle).sin(),
        axis[2] * (0.5 * angle).sin(),
    ];
    let rotation = So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap();
    let mut state = initial(assembly);
    state.motion = assembly
        .state_from_speed(
            SixDofState {
                quat: q,
                vel_mps: a(rotation.rotate(v([12.0, 2.0, 0.8])).unwrap()),
                ..state.motion.rigid
            },
            state.motion.phase_rad,
            assembly.speed(&state.motion).unwrap(),
        )
        .unwrap();
    (state, rotation)
}

fn shear(point_world: [f64; 3]) -> [f64; 3] {
    [
        0.2 + 0.14 * (point_world[1] + 2.0) - 0.06 * (point_world[2] + 10.0),
        -0.1 + 0.025 * (point_world[2] + 10.0),
        0.04 * (point_world[0] - 1.0),
    ]
}

#[test]
fn uniform_world_wind_preserves_exact_hub_loads_and_zero_mode_parity() {
    let mut design = wright_openloop_v1();
    let radius = design.rotor.radius_m;
    // Exercise zero wind increments alongside actual pointwise mutual wake.
    design.disks[0].center_m = [-2.5, 0.0, -0.9];
    design.disks[1].center_m = [-2.5 - radius, 1.5 * radius, -0.9];
    for reverse in [false, true] {
        let assembly = assembly(design.gross_mass_kg, reverse);
        let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(reverse)).unwrap();
        assert!(
            DynamicPropulsion::new(&base, RHO)
                .unwrap()
                .with_blade_incident(BladeIncident::WindField)
                .is_err()
        );
        let hub = DynamicPropulsion::new_oblique(&base, RHO, SAMPLES).unwrap();
        assert_eq!(hub.blade_incident(), BladeIncident::Hub);
        let field = DynamicPropulsion::new_oblique(&base, RHO, SAMPLES)
            .unwrap()
            .with_blade_incident(BladeIncident::WindField)
            .unwrap();
        let harmonic = HarmonicPropulsion::new(&base, RHO, SAMPLES)
            .unwrap()
            .with_blade_incident(BladeIncident::WindField)
            .unwrap();
        let harmonic_hub = HarmonicPropulsion::new(&base, RHO, SAMPLES).unwrap();
        assert_eq!(harmonic.blade_incident(), BladeIncident::WindField);
        for rotated in [false, true] {
            let (mut state, _) = posed(&assembly, rotated);
            let wind = |_: [f64; 3]| Ok([0.35, -0.12, 0.06]);
            assert_eq!(
                harmonic.evaluate(&state, inputs(0.0), &wind).unwrap(),
                harmonic_hub.evaluate(&state, inputs(0.0), &wind).unwrap()
            );
            for mode in &mut state.induction {
                mode.first_mps = [0.0; 2];
            }
            let uniform = DynamicState {
                motion: state.motion,
                induced_mps: state.induction.map(|mode| mode.mean_mps),
            };
            let old = hub.evaluate(&uniform, inputs(0.0), &wind).unwrap();
            let new = field.evaluate(&uniform, inputs(0.0), &wind).unwrap();
            assert_eq!(old, new);
            let modal = harmonic.evaluate(&state, inputs(0.0), &wind).unwrap();
            assert_eq!(new.loads, modal.loads);
            assert_eq!(
                new.aerodynamics.oblique.unwrap(),
                modal.aerodynamics.propellers
            );
        }
    }
}

#[test]
fn world_shear_and_native_wing_field_reach_blades_without_duplicate_motion() {
    let mut design = wright_openloop_v1();
    // Disjoint sources isolate ambient/wing variation from the independently
    // covered mutual-wake field. Both canted, counterrotating shafts remain.
    let radius = design.rotor.radius_m;
    design.disks[0].center_m = [-2.5, -3.0 * radius, -0.9];
    design.disks[1].center_m = [-2.5, 3.0 * radius, -0.9];
    let (panels, _) = design.layout(inputs(0.0).canard_rad).unwrap();
    let rule = fs_airscrew::oblique::AzimuthQuadrature::new(SAMPLES).unwrap();
    for reverse in [false, true] {
        let assembly = assembly(design.gross_mass_kg, reverse);
        let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(reverse)).unwrap();
        for rotated in [false, true] {
            let (state, rotation) = posed(&assembly, rotated);
            let inverse = rotation.inverse();
            let hub = HarmonicPropulsion::new(&base, RHO, SAMPLES)
                .unwrap()
                .evaluate(&state, inputs(0.0), &|p| Ok(shear(p)))
                .unwrap();
            let world_point =
                |point: [f64; 3]| {
                    a(v(state.motion.rigid.pos_m)
                        + rotation.rotate(v(point) - v(design.cg_m)).unwrap())
                };
            let freestream = a(inverse
                .rotate(v(state.motion.rigid.vel_mps) - v(shear(state.motion.rigid.pos_m)))
                .unwrap());
            let mut wind_moments = [ZERO; 2];
            for policy in [BladeIncident::WindField, BladeIncident::WindAndWingField] {
                let model = HarmonicPropulsion::new(&base, RHO, SAMPLES)
                    .unwrap()
                    .with_blade_incident(policy)
                    .unwrap();
                let frames = model.disk_frames();
                let result = model
                    .evaluate(&state, inputs(0.0), &|p| Ok(shear(p)))
                    .unwrap();
                // The field consumes the same accepted wing solution and keeps
                // the same hub incident velocity in wake transport.
                assert_eq!(result.aerodynamics.gamma, hub.aerodynamics.gamma);
                assert_eq!(
                    result.aerodynamics.disk_incident_mps,
                    hub.aerodynamics.disk_incident_mps
                );
                let mut delta_force = Vec3::default();
                let mut delta_moment = Vec3::default();
                let mut shear_couple = 0.0;
                let mut wing_couple = 0.0;
                let mut field_work = 0.0;
                for k in 0..2 {
                    let center = design.disks[k].center_m;
                    let hub_wind = v(shear(world_point(center)));
                    let wing_hub = v(fs_wing::induced_velocity_free(
                        center,
                        &panels,
                        &hub.aerodynamics.gamma,
                        freestream,
                    ));
                    let hub_velocity = inverse
                        .rotate(v(state.motion.rigid.vel_mps) - hub_wind)
                        .unwrap()
                        + v(state.motion.rigid.omega_body).cross(v(center) - v(design.cg_m))
                        + wing_hub;
                    let native = rule
                        .evaluate_with_incident_variation(
                            &design.rotor,
                            RHO,
                            frames[k],
                            a(hub_velocity),
                            state.motion.rigid.omega_body,
                            state.induction[k],
                            result.propeller_absolute_speed_rad_s[k].abs(),
                            |_, arm| {
                                let point = a(v(center) + v(arm));
                                let mut difference = inverse
                                    .rotate(hub_wind - v(shear(world_point(point))))
                                    .unwrap();
                                if policy == BladeIncident::WindAndWingField {
                                    difference = difference
                                        + (v(fs_wing::induced_velocity_free(
                                            point,
                                            &panels,
                                            &hub.aerodynamics.gamma,
                                            freestream,
                                        )) - wing_hub);
                                }
                                Ok(a(difference))
                            },
                        )
                        .unwrap();
                    let actual = result.aerodynamics.propellers[k];
                    assert!((v(actual.force_n) - v(native.force_n)).norm() < 1e-9);
                    assert!((v(actual.hub_moment_nm) - v(native.hub_moment_nm)).norm() < 1e-9);
                    assert!((actual.torque_nm - native.torque_nm).abs() < 1e-10);
                    assert!(
                        (actual.incident_variation_work_w - native.incident_variation_work_w).abs()
                            < 1e-9
                    );
                    assert!(actual.power_residual_w.abs() < 2e-9);
                    field_work += actual.incident_variation_work_w.abs();
                    let before = hub.aerodynamics.propellers[k];
                    shear_couple += (v(actual.hub_moment_nm) - v(before.hub_moment_nm)).norm();
                    if policy == BladeIncident::WindField {
                        wind_moments[k] = actual.hub_moment_nm;
                    } else {
                        wing_couple += (v(actual.hub_moment_nm) - v(wind_moments[k])).norm();
                    }
                    let change = v(actual.force_n) - v(before.force_n);
                    delta_force = delta_force + change;
                    delta_moment = delta_moment
                        + (v(center) - v(design.cg_m)).cross(change)
                        + (v(actual.hub_moment_nm) - v(before.hub_moment_nm));
                    let sign = if reverse { -1.0 } else { 1.0 };
                    assert_eq!(
                        result.loads.aerodynamic_nm[k],
                        sign * [-1.0, 1.0][k] * actual.torque_nm
                    );
                }
                assert!(
                    shear_couple > 0.01,
                    "spatial air must produce a real hub couple"
                );
                assert!(field_work > 1e-4);
                if policy == BladeIncident::WindAndWingField {
                    assert!(
                        wing_couple > 0.01,
                        "native wing variation must reach the blades"
                    );
                }
                assert!(
                    (v(result.aerodynamic_force_body_n)
                        - v(hub.aerodynamic_force_body_n)
                        - delta_force)
                        .norm()
                        < 2e-10
                );
                assert!(
                    (v(result.loads.carrier.moment_nm)
                        - v(hub.loads.carrier.moment_nm)
                        - delta_moment)
                        .norm()
                        < 2e-10
                );
            }
        }
    }
}

#[test]
fn late_blade_wind_refusal_preserves_state_and_exact_retry() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let hub = HarmonicPropulsion::new(&base, RHO, SAMPLES).unwrap();
    let field = HarmonicPropulsion::new(&base, RHO, SAMPLES)
        .unwrap()
        .with_blade_incident(BladeIncident::WindAndWingField)
        .unwrap();
    let (state, _) = posed(&assembly, true);
    let before = state;
    let calls = Cell::new(0);
    hub.evaluate(&state, inputs(0.0), &|p| {
        calls.set(calls.get() + 1);
        Ok(shear(p))
    })
    .unwrap();
    let fail_at = calls.get() + 10;
    let refusal = Refusal {
        code: "blade-material-air-unavailable",
        message: "the tenth native blade sample is unavailable".into(),
        ranked_repairs: vec!["retry with the complete material-air field".into()],
    };
    calls.set(0);
    let error = field
        .evaluate(&state, inputs(0.0), &|p| {
            calls.set(calls.get() + 1);
            if calls.get() == fail_at {
                Err(refusal.clone())
            } else {
                Ok(shear(p))
            }
        })
        .unwrap_err();
    assert_eq!(error, refusal);
    assert_eq!(calls.get(), fail_at);
    assert_eq!(state, before);
    let clean = field
        .step(
            &state,
            0.0,
            0.004,
            ZERO,
            |t, _| Ok(inputs(t)),
            |_, p| Ok(shear(p)),
        )
        .unwrap();
    calls.set(0);
    let error = field
        .step(
            &state,
            0.0,
            0.004,
            ZERO,
            |t, _| Ok(inputs(t)),
            |t, p| {
                if t > 0.0 {
                    calls.set(calls.get() + 1);
                    if calls.get() == fail_at {
                        return Err(refusal.clone());
                    }
                }
                Ok(shear(p))
            },
        )
        .unwrap_err();
    assert_eq!(error, refusal);
    assert_eq!(calls.get(), fail_at);
    assert_eq!(state, before);
    let retry = field
        .step(
            &state,
            0.0,
            0.004,
            ZERO,
            |t, _| Ok(inputs(t)),
            |_, p| Ok(shear(p)),
        )
        .unwrap();
    assert_eq!(clean, retry);
}

#[test]
fn resolved_air_field_trajectory_closes_actual_world_momentum() {
    let design = wright_openloop_v1();
    let assembly = assembly(design.gross_mass_kg, false);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts(false)).unwrap();
    let model = HarmonicPropulsion::new(&base, RHO, SAMPLES)
        .unwrap()
        .with_blade_incident(BladeIncident::WindAndWingField)
        .unwrap();
    let (mut state, _) = posed(&assembly, true);
    let origin = [-4.0, 2.0, 1.0];
    let p0 = v(state.motion.rigid.vel_mps).scale(assembly.mass_kg());
    let l0 = angular_momentum(&assembly, &state.motion, origin);
    let mut dp = Vec3::default();
    let mut dl = Vec3::default();
    let mut air_work = 0.0;
    for tick in 0..12 {
        let step = model
            .step(
                &state,
                tick as f64 * 0.004,
                0.004,
                origin,
                |t, _| Ok(inputs(t)),
                |t, p| {
                    let mut wind = shear(p);
                    wind[0] += 0.1 * t;
                    Ok(wind)
                },
            )
            .unwrap();
        dp = dp + v(step.motion.external_linear_impulse_world);
        dl = dl + v(step.motion.external_angular_impulse_world);
        for blade in step.midpoint.aerodynamics.propellers {
            assert!(blade.power_residual_w.abs() < 2e-9);
            air_work += 0.004 * blade.incident_variation_work_w.abs();
        }
        for wake in step.inflow {
            assert!(v(wake.momentum_residual_ns).norm() < 2e-12);
        }
        state = step.state;
    }
    let p_defect = (v(state.motion.rigid.vel_mps).scale(assembly.mass_kg()) - p0 - dp).norm();
    let l_defect = (angular_momentum(&assembly, &state.motion, origin) - l0 - dl).norm();
    assert!(p_defect < 2e-9);
    assert!(l_defect < 2e-8);
    assert!(air_work > 1e-4);
}
