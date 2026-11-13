//! Original inertia, finite thermal feedback and source-pinned moving beams.

use super::*;
use fs_solid::frame2d::{Rayleigh, laminate_thermal::dynamics::FrameLaminateNewmarkControls};

fn motion_controls() -> FrameLaminateNewmarkControls {
    FrameLaminateNewmarkControls {
        maximum_dofs: 6,
        maximum_members: 1,
        maximum_force_residual_n: 1e-7,
        maximum_moment_residual_n_m: 1e-9,
        maximum_displacement_defect_m: 1e-12,
        maximum_rotation_defect_rad: 1e-11,
        maximum_energy_defect_j: 5e-10,
    }
}

#[test]
fn cold_to_hot_native_motion_matches_the_analytical_inertial_rod() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (mass_form, mass_fraction) in [
            (FrameMassFormulation::Lumped, 0.5),
            (FrameMassFormulation::ConsistentReference, 1.0 / 3.0),
        ] {
            for alpha in [1e-3, -1e-3] {
                for force in [0.0, 0.5] {
                    let p = ply(alpha);
                    let (length, dt, heat) = (0.2, 2e-4, 0.2);
                    let area = WIDTH * p.h;
                    let volume = area * length;
                    let mass = mass_fraction * RHO * volume;
                    fixture_with_mass(
                        &[p],
                        length,
                        &[REFERENCE],
                        &[heat / dt],
                        false,
                        mass_form,
                        cx,
                        |mut frame, thermal, source| {
                            frame.fix(1, [false, true, true]);
                            let load = [0.0, 0.0, 0.0, force, 0.0, 0.0];
                            let model = ReciprocalLaminateFrame::try_new_dynamic(
                                frame,
                                &thermal,
                                source,
                                &mapping(1),
                                &load,
                                LaminateCaloricReference::StressFreeSolid,
                                3,
                                motion_controls(),
                                cx,
                            )
                            .unwrap();
                            let before_mass = model
                                .frame()
                                .laminate_thermal_section(0, 1)
                                .unwrap()
                                .prepared()
                                .reference_mass_moments();
                            let result = model
                                .step_candidate(&thermal, dt, controls(), &[heat], cx)
                                .unwrap();

                            // Solve the ORIGINAL one-free-DOF equations. Cold
                            // a0=force/m; using the hot initial force would add
                            // the thermal forcing twice and fail this oracle.
                            let k = p.e * area / length;
                            let a = dt * dt * k / (4.0 * mass);
                            let r = a / (1.0 + a);
                            let imposed_strain = 2.0 * force * r / (p.e * area);
                            let cp = RHO * p.cp * volume;
                            let b = p.e * alpha * alpha * volume * (1.0 - r);
                            let d = p.e * alpha * volume * imposed_strain;
                            let linear = cp - REFERENCE * b + 0.5 * d;
                            let rhs = heat - REFERENCE * d;
                            let theta =
                                2.0 * rhs / (linear + (linear * linear - 2.0 * b * rhs).sqrt());
                            let strain = imposed_strain + alpha * r * theta;
                            let u = length * strain;
                            let v = 2.0 * u / dt;
                            let kinetic = 0.5 * mass * v * v;
                            let elastic = volume
                                * p.e
                                * (0.5 * strain * strain + REFERENCE * alpha * strain
                                    - alpha * alpha * REFERENCE * theta
                                    - 0.5 * alpha * alpha * theta * theta);
                            let report = result.report();
                            let dynamic = report.mechanical.dynamics.as_ref().unwrap();
                            near(
                                result.thermal().temperature_k(CellId(0)).unwrap(),
                                REFERENCE + theta,
                                2e-8,
                            );
                            near(result.frame().frame().displacements()[3], u, 2e-11);
                            near(result.frame().frame().velocities()[3], v, 2e-7);
                            near(dynamic.kinetic_energy_change_j, kinetic, 2e-10);
                            near(report.mechanical.internal_energy_change_j, elastic, 2e-9);
                            near(report.mechanical.external_work_j, force * u, 2e-10);
                            near(report.mechanical.reversible_heat_j, cp * theta - heat, 2e-9);
                            near(report.external_thermal_heat_j, heat, 2e-11);
                            near(report.mechanical.work_residual_j, 0.0, 5e-10);
                            near(dynamic.native_energy_residual_j, 0.0, 5e-10);
                            near(report.total_energy_residual_j, 0.0, 5e-10);
                            assert_eq!(dynamic.numerical_dissipation_j, 0.0);
                            assert_eq!(dynamic.history.damping_energy, [0.0; 2]);
                            assert_eq!(
                                before_mass,
                                result
                                    .frame()
                                    .frame()
                                    .laminate_thermal_section(0, 1)
                                    .unwrap()
                                    .prepared()
                                    .reference_mass_moments()
                            );
                            assert_eq!(model.frame().displacements(), &[0.0; 6]);
                            assert_eq!(model.frame().velocities(), &[0.0; 6]);
                            assert_eq!(thermal.time_s(), 0.0);
                            assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), REFERENCE);
                        },
                    );
                }
            }
        }
    });
}

// Independent continuous first-law ODE for the homogeneous free-sided rod.
// This scalar reference is test-only; production uses native Newmark and the
// existing implicit caloric owner, not this integrator.
fn rod_reference(p: Ply, length: f64, mass: f64, speed: f64, duration: f64) -> [f64; 3] {
    let area = WIDTH * p.h;
    let volume = area * length;
    let cp = RHO * p.cp * volume;
    let b = p.e * p.alpha * p.alpha * volume;
    let derivative = |[u, v, t]: [f64; 3]| {
        [
            v,
            (-p.e * area * u / length + p.e * area * p.alpha * (t - REFERENCE)) / mass,
            -t * p.e * p.alpha * area * v / (cp - t * b),
        ]
    };
    let mut state = [0.0, speed, REFERENCE];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = derivative(state);
        let b = derivative(std::array::from_fn(|i| state[i] + 0.5 * dt * a[i]));
        let c = derivative(std::array::from_fn(|i| state[i] + 0.5 * dt * b[i]));
        let d = derivative(std::array::from_fn(|i| state[i] + dt * c[i]));
        state =
            std::array::from_fn(|i| state[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn closed_impulsed_motion_keeps_kinetic_storage_and_refines_to_thermoelastic_oscillation() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let p = ply(1e-3);
        let (length, speed) = (0.2, 0.1);
        let area = WIDTH * p.h;
        let volume = area * length;
        let mass = RHO * volume / 3.0;
        let cp = RHO * p.cp * volume;
        let b = p.e * p.alpha * p.alpha * volume;
        let omega = (p.e * area / length / mass * cp / (cp - REFERENCE * b)).sqrt();
        let duration = std::f64::consts::TAU / omega;
        let exact = rod_reference(p, length, mass, speed, duration);
        let displacement_scale = speed / omega;
        let temperature_scale =
            REFERENCE * p.e * p.alpha * area * displacement_scale / (cp - REFERENCE * b);
        let mut errors = Vec::new();
        for count in [12, 24, 48] {
            fixture_with_mass(
                &[p],
                length,
                &[REFERENCE],
                &[0.0],
                false,
                FrameMassFormulation::ConsistentReference,
                cx,
                |mut frame, mut thermal, source| {
                    frame.fix(1, [false, true, true]);
                    let dt = duration / count as f64;
                    frame
                        .apply_impulse(
                            &[0.0, 0.0, 0.0, mass * speed, 0.0, 0.0],
                            &[0.0; 6],
                            Rayleigh { a0: 0.0, a1: 0.0 },
                            dt,
                        )
                        .unwrap();
                    assert!(
                        ReciprocalLaminateFrame::try_new(
                            frame.clone(),
                            &thermal,
                            source,
                            &mapping(1),
                            &[0.0; 6],
                            LaminateCaloricReference::StressFreeSolid,
                            3,
                            cx,
                        )
                        .is_err()
                    );
                    let original_velocity = frame.velocities().to_vec();
                    let mut model = ReciprocalLaminateFrame::try_new_dynamic(
                        frame,
                        &thermal,
                        source,
                        &mapping(1),
                        &[0.0; 6],
                        LaminateCaloricReference::StressFreeSolid,
                        3,
                        motion_controls(),
                        cx,
                    )
                    .unwrap();
                    assert_eq!(model.frame().velocities(), original_velocity);
                    let initial_kinetic = 0.5 * mass * speed * speed;
                    let mut kinetic_change = 0.0;
                    for step in 1..=count {
                        let result = model
                            .step_candidate(&thermal, step as f64 * dt, controls(), &[0.0], cx)
                            .unwrap();
                        let report = result.report();
                        let dynamic = report.mechanical.dynamics.as_ref().unwrap();
                        kinetic_change += dynamic.kinetic_energy_change_j;
                        near(report.external_thermal_heat_j, 0.0, 2e-11);
                        near(report.mechanical.external_work_j, 0.0, 2e-11);
                        near(report.total_energy_residual_j, 0.0, 5e-10);
                        assert_eq!(dynamic.numerical_dissipation_j, 0.0);
                        if step == 1 {
                            assert!(result.thermal().temperature_k(CellId(0)).unwrap() < REFERENCE);
                        }
                        (thermal, model, _) = result.into_parts();
                    }
                    let state = [
                        model.frame().displacements()[3],
                        model.frame().velocities()[3],
                        thermal.temperature_k(CellId(0)).unwrap(),
                    ];
                    let kinetic = 0.5 * mass * state[1] * state[1];
                    near(kinetic - initial_kinetic, kinetic_change, 1e-11);
                    let stored = model
                        .frame()
                        .laminate_thermoelastic_energy(0, cx)
                        .unwrap()
                        .internal_energy_j;
                    near(
                        cp * (state[2] - REFERENCE) + stored + kinetic - initial_kinetic,
                        0.0,
                        2e-8,
                    );
                    errors.push(
                        ((state[0] - exact[0]) / displacement_scale).abs()
                            + ((state[1] - exact[1]) / speed).abs()
                            + ((state[2] - exact[2]) / temperature_scale).abs(),
                    );
                },
            );
        }
        assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.4 * errors[1], "{errors:?}");
        assert!(errors[2] < 0.03, "{errors:?}");
    });
}

#[test]
fn heated_bimetal_moves_and_refreshes_sources_with_atomic_failed_trials() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut bottom = ply(5e-4);
        bottom.h = 0.001;
        bottom.nu = 0.14;
        bottom.cp = 500.0;
        let mut top = bottom;
        top.alpha = 1e-3;
        top.e = 8e7;
        top.slope = -5e5;
        top.nu = 0.34;
        let (dt, heat) = (0.001, 0.05);
        fixture_with_mass(
            &[bottom, top],
            0.05,
            &[REFERENCE; 2],
            &[0.0, heat / dt],
            false,
            FrameMassFormulation::ConsistentReference,
            cx,
            |frame, thermal, source| {
                let model = ReciprocalLaminateFrame::try_new_dynamic(
                    frame,
                    &thermal,
                    source,
                    &mapping(2),
                    &[0.0; 6],
                    LaminateCaloricReference::StressFreeSolid,
                    6,
                    motion_controls(),
                    cx,
                )
                .unwrap();
                let before_mass = model
                    .frame()
                    .laminate_thermal_section(0, 1)
                    .unwrap()
                    .prepared()
                    .reference_mass_moments();
                assert!(
                    model
                        .step_candidate(&thermal, 100.0, controls(), &[0.0, 5000.0], cx)
                        .is_err()
                );
                let gate = CancelGate::new_clock_free();
                gate.request();
                with_cx(&gate, |cancelled| {
                    assert!(matches!(
                        model.step_candidate(&thermal, dt, controls(), &[0.0, heat], cancelled),
                        Err(CircuitEnthalpyError::Cancelled)
                    ));
                });
                let result = model
                    .step_candidate(&thermal, dt, controls(), &[0.0, heat], cx)
                    .unwrap();
                let retry = model
                    .step_candidate(&thermal, dt, controls(), &[0.0, heat], cx)
                    .unwrap();
                assert_eq!(result.report().mechanical, retry.report().mechanical);
                assert_eq!(result.thermal().cells(), retry.thermal().cells());
                let moving = result.frame().frame();
                assert!(moving.displacements()[4] < -1e-8);
                assert!(moving.velocities()[4] < -1e-5);
                assert!(
                    moving
                        .assess_laminate_frame_strength(cx)
                        .unwrap()
                        .utilization()
                        > 0.0
                );
                let section = moving.laminate_thermal_section(0, 1).unwrap();
                assert_eq!(before_mass, section.prepared().reference_mass_moments());
                for layer in 0..2 {
                    near(
                        section
                            .thermoelastic()
                            .unwrap()
                            .layer_temperature_k(layer)
                            .unwrap(),
                        result
                            .thermal()
                            .temperature_k(CellId(layer as u64))
                            .unwrap(),
                        0.0,
                    );
                }
                let report = result.report();
                let dynamic = report.mechanical.dynamics.as_ref().unwrap();
                assert!(dynamic.kinetic_energy_change_j > 0.0);
                assert!(
                    report
                        .mechanical
                        .layers
                        .iter()
                        .any(|layer| layer.source.reversible_heat_j_m != 0.0)
                );
                near(report.external_thermal_heat_j, heat, 2e-11);
                near(report.mechanical.work_residual_j, 0.0, 5e-10);
                near(report.total_energy_residual_j, 0.0, 5e-10);
                assert_eq!(dynamic.numerical_dissipation_j, 0.0);
                assert_eq!(model.frame().displacements(), &[0.0; 6]);
                assert_eq!(model.frame().velocities(), &[0.0; 6]);
                assert_eq!(model.frame().accelerations(), &[0.0; 6]);
                assert_eq!(thermal.time_s(), 0.0);
            },
        );
    });
}
