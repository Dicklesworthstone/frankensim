//! Native unequal-shell inventories, independent diffusion modes and dry fronts.

use super::*;
use fs_couple::moisture::adaptive::*;

#[path = "adaptive/finite.rs"]
mod finite;
#[path = "adaptive/solution.rs"]
mod solution;

fn controls(step_s: f64, tolerance: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: step_s,
        minimum_step_s: 1e-7,
        maximum_step_s: step_s,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: tolerance,
        activity_absolute_tolerance: tolerance,
        reservoir_water_absolute_tolerance_kg: 1.0,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

fn two_shell_parameters(symmetry: RadialMoistureSymmetry) -> ([f64; 2], f64) {
    let radii: [f64; 3] = [0.02, 0.035, 0.1];
    let centers = [0.5 * (radii[0] + radii[1]), 0.5 * (radii[1] + radii[2])];
    let (dimension, factor, resistance) = match symmetry {
        RadialMoistureSymmetry::Cylinder { length_m } => {
            let factor = 2.0 * core::f64::consts::PI * length_m;
            (2, factor, (centers[1] / centers[0]).ln() / factor)
        }
        RadialMoistureSymmetry::Sphere => {
            let factor = 4.0 * core::f64::consts::PI;
            (3, factor, (1.0 / centers[0] - 1.0 / centers[1]) / factor)
        }
    };
    // Independent analytic shell volumes and curved diffusion resistance.
    let mass = core::array::from_fn(|i| {
        500.0 * factor * (radii[i + 1].powi(dimension) - radii[i].powi(dimension))
            / f64::from(dimension)
    });
    let rate = 500.0 * 0.001 / resistance * (1.0 / mass[0] + 1.0 / mass[1]);
    (mass, rate)
}

fn two_shell_mode(mass: [f64; 2], decay: f64) -> [f64; 2] {
    let total = mass[0] + mass[1];
    let mean = (mass[0] * 0.22 + mass[1] * 0.08) / total;
    let difference = (0.22 - 0.08) * decay;
    [
        mean + mass[1] / total * difference,
        mean - mass[0] / total * difference,
    ]
}

fn check_closed<Step>(
    original: &MoistureRadial,
    native: (MoistureRadial, MoistureAdvanceReport<Step, 2>),
    physical_initial_water: f64,
) -> (MoistureRadial, usize) {
    let (next, report) = native;
    assert_eq!(report.interval_s, [0.0, 0.25]);
    assert_eq!(next.time_s(), 0.25);
    assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
    assert_eq!(report.boundary_water_kg, [0.0; 2]);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    assert!(report.rejected_intervals > 0);
    near(report.water_change_kg, 0.0, 1e-11);
    near(report.stored_water_kg[0], physical_initial_water, 1e-12);
    near(report.stored_water_kg[1], physical_initial_water, 1e-11);
    near(next.water_mass_kg(), original.water_mass_kg(), 1e-11);
    near(report.transport_balance_residual_kg, 0.0, 1e-11);
    (next, report.evaluations)
}

#[test]
fn radial_adaptive_refines_both_methods_on_unequal_cylindrical_and_spherical_shells() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
            RadialMoistureSymmetry::Sphere,
        ] {
            let original = specimen(
                symmetry,
                vec![0.02, 0.035, 0.1],
                vec![0.22, 0.08],
                [0.001; 2],
                cx,
            );
            let (mass, rate) = two_shell_parameters(symmetry);
            let physical_water = mass[0] * 0.22 + mass[1] * 0.08;
            let exact = two_shell_mode(mass, (-rate * 0.25).exp());
            let sealed = [MoistureFace::Sealed, MoistureFace::Sealed];
            for second_order in [false, true] {
                let mut errors = Vec::new();
                let mut evaluations = Vec::new();
                for tolerance in [1e-4, 1e-6] {
                    let time = controls(0.25, tolerance);
                    let (next, work) = if second_order {
                        check_closed(
                            &original,
                            original
                                .advance_sdirk2_adaptive_to(
                                    0.25,
                                    &sealed,
                                    config(),
                                    nonlinear(),
                                    time,
                                    cx,
                                )
                                .unwrap(),
                            physical_water,
                        )
                    } else {
                        check_closed(
                            &original,
                            original
                                .advance_adaptive_to(0.25, &sealed, config(), nonlinear(), time, cx)
                                .unwrap(),
                            physical_water,
                        )
                    };
                    evaluations.push(work);
                    errors.push(
                        next.moisture()
                            .iter()
                            .zip(exact)
                            .map(|(&actual, expected)| (actual - expected).abs())
                            .fold(0.0_f64, f64::max),
                    );
                }
                assert!(errors[1] < 0.2 * errors[0], "{symmetry:?}: {errors:?}");
                assert!(errors[1] < if second_order { 1e-5 } else { 1e-4 });
                assert!(evaluations[1] > evaluations[0]);
            }
            assert_eq!(original.moisture(), &[0.22, 0.08]);
            assert_eq!(original.time_s(), 0.0);
        }
    });
}

fn check_water_comparison<Step>(report: &MoistureAdvanceReport<Step, 2>, work: usize) {
    assert_eq!(report.accepted_intervals, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(report.evaluations, work);
    near(report.maximum_accepted_error_ratio, 0.5, 2e-8);
    near(report.water_change_kg, 0.0, 1e-12);
}

#[test]
fn radial_temporal_water_error_weights_each_shell_and_cannot_cancel_interior_exchange() {
    with_cx(|cx| {
        let symmetry = RadialMoistureSymmetry::Sphere;
        let original = specimen(
            symmetry,
            vec![0.02, 0.035, 0.1],
            vec![0.22, 0.08],
            [0.001; 2],
            cx,
        );
        let (mass, rate) = two_shell_parameters(symmetry);
        assert!(mass[1] > 25.0 * mass[0]);
        let duration = 0.05;
        let sealed = [MoistureFace::Sealed, MoistureFace::Sealed];
        for second_order in [false, true] {
            let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
            let stability = |z: f64| {
                if second_order {
                    (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
                } else {
                    1.0 / (1.0 + z)
                }
            };
            let coarse = two_shell_mode(mass, stability(rate * duration));
            let fine = two_shell_mode(mass, stability(0.5 * rate * duration).powi(2));
            let fine_error_factor = if second_order { 1.0 / 3.0 } else { 1.0 };
            // The closed specimen has zero net water error but two nonzero,
            // opposing local errors. The accepted comparison is exactly 1/2
            // only with both physical dry masses and an absolute (L1) sum.
            let physical_difference =
                mass[0] * (coarse[0] - fine[0]).abs() + mass[1] * (coarse[1] - fine[1]).abs();
            let mut time = controls(duration, 1.0);
            time.reservoir_water_absolute_tolerance_kg =
                2.0 * fine_error_factor * physical_difference;
            if second_order {
                let (_, report) = original
                    .advance_sdirk2_adaptive_to(duration, &sealed, config(), nonlinear(), time, cx)
                    .unwrap();
                check_water_comparison(&report, 6);
            } else {
                let (_, report) = original
                    .advance_adaptive_to(duration, &sealed, config(), nonlinear(), time, cx)
                    .unwrap();
                check_water_comparison(&report, 3);
            }
        }
    });
}

#[test]
fn radial_adaptive_nonlinear_wetting_refines_the_dry_front_and_accounts_boundary_water() {
    with_cx(|cx| {
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.0],
            [0.0, 0.03],
            cx,
        );
        let faces = [MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)];
        // Independent nonlinear continuum-in-time one-shell law:
        // M'=15*(0.2^2-M^2), M(0)=0, so M(t)=0.2*tanh(3t).
        let exact = 0.2 * 0.75_f64.tanh();
        let mut errors = Vec::new();
        for tolerance in [1e-4, 1e-7] {
            let (next, report) = original
                .advance_sdirk2_adaptive_to(
                    0.25,
                    &faces,
                    config(),
                    nonlinear(),
                    controls(0.25, tolerance),
                    cx,
                )
                .unwrap();
            errors.push((next.moisture()[0] - exact).abs());
            assert!(next.moisture()[0] > 0.0 && next.moisture()[0] < 0.2);
            assert_eq!(report.boundary_water_kg[0], 0.0);
            assert_eq!(report.stored_water_kg[0], 0.0);
            near(report.stored_water_kg[1], next.water_mass_kg(), 1e-13);
            near(report.boundary_water_kg[1], next.water_mass_kg(), 1e-10);
            near(report.water_change_kg, next.water_mass_kg(), 1e-12);
            assert!(report.maximum_accepted_error_ratio <= 1.0);
            assert!(report.last_substep.endpoint.kirchhoff_source.is_some());
        }
        assert!(errors[1] < 0.05 * errors[0], "{errors:?}");
        assert!(errors[1] < 2e-7, "{errors:?}");
        assert_eq!(original.moisture(), &[0.0]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn radial_adaptive_sdirk_retries_physical_range_failures_without_clipping() {
    with_cx(|cx| {
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.2],
            [0.01; 2],
            cx,
        );
        let faces = [MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.0)];
        assert!(matches!(
            original.advance_sdirk2(1.0, &faces, config(), nonlinear(), cx),
            Err(MoistureError::StageRange { .. })
        ));
        let (next, report) = original
            .advance_sdirk2_adaptive_to(1.0, &faces, config(), nonlinear(), controls(1.0, 1e-6), cx)
            .unwrap();
        assert!(report.rejected_intervals > 0);
        assert!(report.evaluations >= 6 * report.accepted_intervals + 2);
        near(next.moisture()[0], 0.2 * (-3.0_f64).exp(), 1e-5);
        assert!(next.moisture()[0] > 0.0);
        assert!(report.boundary_water_kg[1] < 0.0);
        near(
            report.boundary_water_kg[1],
            next.water_mass_kg() - original.water_mass_kg(),
            1e-10,
        );
        assert_eq!(original.moisture(), &[0.2]);
        assert_eq!(original.time_s(), 0.0);
    });
}

#[test]
fn radial_adaptive_work_failure_discards_an_accepted_prefix() {
    with_cx(|cx| {
        let original = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.1],
            [0.01; 2],
            cx,
        );
        let sealed = [MoistureFace::Sealed, MoistureFace::Sealed];
        let mut time = controls(0.125, 1e-5);
        time.maximum_evaluations = 6;
        assert!(matches!(
            original.advance_sdirk2_adaptive_to(0.25, &sealed, config(), nonlinear(), time, cx,),
            Err(MoistureAdvanceError::WorkLimit {
                time_s: 0.125,
                accepted_intervals: 1,
                evaluations: 6,
            })
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.moisture(), &[0.1]);
    });
}

#[test]
fn radial_adaptive_cancellation_preserves_the_original_checkpoint() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let original = specimen(
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
            vec![0.0, 0.03, 0.1],
            vec![0.2, 0.1],
            [0.01; 2],
            &cx,
        );
        let sealed = [MoistureFace::Sealed, MoistureFace::Sealed];
        gate.request();
        assert!(matches!(
            original.advance_sdirk2_adaptive_to(
                1.0,
                &sealed,
                config(),
                nonlinear(),
                controls(1.0, 1e-5),
                &cx,
            ),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert!(matches!(
            original.advance_adaptive_to(
                1.0,
                &sealed,
                config(),
                nonlinear(),
                controls(1.0, 1e-5),
                &cx,
            ),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert_eq!(original.moisture(), &[0.2, 0.1]);
        assert_eq!(original.time_s(), 0.0);
    });
}
