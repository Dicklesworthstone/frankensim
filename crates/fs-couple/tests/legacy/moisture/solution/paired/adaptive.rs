//! G1/G3/G4 controlled physical durations through the original paired owners.

use super::*;
use fs_couple::moisture::adaptive::{
    KineticPhaseAccuracy, MoistureAdvanceError, MoistureTimeControls,
    PairedSolutionMoistureTimeControls,
};

fn time(span: f64, tolerance: f64) -> PairedSolutionMoistureTimeControls {
    PairedSolutionMoistureTimeControls {
        transport: MoistureTimeControls {
            initial_step_s: span,
            minimum_step_s: 1e-6,
            maximum_step_s: span,
            maximum_evaluations: 100_000,
            moisture_absolute_tolerance: tolerance,
            activity_absolute_tolerance: 1.0,
            reservoir_water_absolute_tolerance_kg: 1.0,
            relative_tolerance: 0.0,
            water_balance_tolerance_kg: 5e-11,
        },
        kinetic_phases: [None, None],
        maximum_field_evaluations: 100_000,
    }
}

fn physical(volumes: [f64; 2]) -> PairedSolutionMoistureSolve {
    let mut controls = pair_controls(volumes);
    for bath in &mut controls.solutions {
        bath.mass_tolerance_kg = 1e-15;
    }
    controls
}

#[test]
fn g1_pure_bath_time_refinement_approaches_the_native_one_cell_exponential() {
    with_cx(|cx| {
        let volume = 1e-3;
        let baths = [reservoir(volume, vec![0.0]), reservoir(volume, vec![0.0])];
        let field = slab(card([0.01; 2], 0.3), 1.0, vec![0.05], cx);
        let exact = 0.2 - 0.15 * (-4.0_f64 * 0.01 * 4.0).exp();
        let mut errors = Vec::new();
        let mut counts = Vec::new();
        for tolerance in [1e-4, 2.5e-5, 6.25e-6] {
            let (next, solid, report) = baths[0]
                .advance_paired_adaptive_to(
                    &baths[1],
                    &field,
                    4.0,
                    &pair_faces(0.2),
                    config(),
                    physical([volume; 2]),
                    time(4.0, tolerance),
                    cx,
                )
                .unwrap();
            errors.push((solid.moisture()[0] - exact).abs());
            counts.push(report.transport.evaluations);
            assert!(report.transport.maximum_accepted_error_ratio <= 1.0);
            assert_eq!(report.field_evaluations, report.transport.evaluations);
            assert_eq!(solid.time_s(), 4.0);
            let gained = field.dry_mass_kg() * (solid.moisture()[0] - 0.05);
            for port in 0..2 {
                near(
                    report.solutions[port].water_to_solid_kg,
                    gained * 0.5,
                    5e-12,
                );
                near(
                    next[port].solvent_mass_kg(),
                    WATER_DENSITY * volume - gained * 0.5,
                    5e-12,
                );
                near(report.solutions[port].water_balance_residual_kg, 0.0, 5e-11);
            }
            near(report.direct_water_balance_residual_kg, 0.0, 5e-11);
        }
        assert!(
            errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
            "{errors:?}"
        );
        assert!(counts[0] < counts[1] && counts[1] < counts[2], "{counts:?}");
        assert_eq!(field.moisture(), &[0.05]);
        assert_eq!(field.time_s(), 0.0);
    });
}

// Independent continuous finite-volume balance for one cell, two ideal baths
// and optional fixed-surface dissolution. No production endpoint/root is called.
// y = [dry-basis moisture, bath water kg, bath water kg, dissolved mol, dissolved mol].
fn reference(
    initial: [f64; 5],
    duration: f64,
    intervals: usize,
    dry_mass: f64,
    face_rate: f64,
    rates: [f64; 2],
    saturation: [f64; 2],
) -> [f64; 5] {
    let derivative = |y: [f64; 5]| {
        let volume = [y[1] / WATER_DENSITY, y[2] / WATER_DENSITY];
        let flux: [f64; 2] = std::array::from_fn(|i| {
            let activity = (-WATER_VOLUME * y[3 + i] / volume[i]).exp();
            dry_mass * face_rate * (0.2 * activity - y[0])
        });
        [
            (flux[0] + flux[1]) / dry_mass,
            -flux[0],
            -flux[1],
            rates[0] * (saturation[0] - y[3] / volume[0]),
            rates[1] * (saturation[1] - y[4] / volume[1]),
        ]
    };
    let dt = duration / intervals as f64;
    let mut y = initial;
    for _ in 0..intervals {
        let k1 = derivative(y);
        let k2 = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]));
        let k3 = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]));
        let k4 = derivative(std::array::from_fn(|i| y[i] + dt * k3[i]));
        y = std::array::from_fn(|i| y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}

#[test]
fn g1_independent_bath_water_accuracy_resolves_nearly_stationary_through_flow() {
    with_cx(|cx| {
        let volume = 2e-6;
        let activity = [0.997_f64, 0.992];
        let amounts = activity.map(|a| -a.ln() * volume / WATER_VOLUME);
        let initial_m = 0.1 * (activity[0] + activity[1]);
        let baths = amounts.map(|n| reservoir(volume, vec![n]));
        let field = slab(card([0.01; 2], 0.3), 1.0, vec![initial_m], cx);
        let initial = [
            initial_m,
            volume * WATER_DENSITY,
            volume * WATER_DENSITY,
            amounts[0],
            amounts[1],
        ];
        let oracle = reference(
            initial,
            4.0,
            4096,
            field.dry_mass_kg(),
            0.02,
            [0.0; 2],
            [0.0; 2],
        );
        let check = reference(
            initial,
            4.0,
            2048,
            field.dry_mass_kg(),
            0.02,
            [0.0; 2],
            [0.0; 2],
        );
        for i in 0..5 {
            near(oracle[i], check[i], 2e-13);
        }
        let loose = time(4.0, 1.0);
        let mut tight = loose;
        tight.transport.reservoir_water_absolute_tolerance_kg = 1e-9;
        let run = |controls| {
            baths[0]
                .advance_paired_adaptive_to(
                    &baths[1],
                    &field,
                    4.0,
                    &pair_faces(0.2),
                    config(),
                    physical([volume; 2]),
                    controls,
                    cx,
                )
                .unwrap()
        };
        let (coarse_baths, _, coarse) = run(loose);
        let (fine_baths, fine_field, fine) = run(tight);
        let error = |b: &[FiniteMoistureSolution; 2]| {
            (0..2)
                .map(|i| (b[i].solvent_mass_kg() - oracle[i + 1]).abs())
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine_baths) < 0.5 * error(&coarse_baths));
        assert!(fine.transport.rejected_intervals > 0);
        assert!(fine.field_evaluations > coarse.field_evaluations);
        assert!(fine.solutions[0].water_to_solid_kg > 0.0);
        assert!(fine.solutions[1].water_to_solid_kg < 0.0);
        assert!(
            fine.transport.water_change_kg.abs() < 0.1 * fine.solutions[0].water_to_solid_kg.abs()
        );
        near(fine_field.moisture()[0], oracle[0], 2e-7);
        for i in 0..2 {
            near(
                fine.solutions[i].water_to_solid_kg,
                fine.transport.boundary_water_kg[i],
                1e-15,
            );
            near(fine.solutions[i].water_balance_residual_kg, 0.0, 5e-11);
            assert_eq!(
                fine_baths[i].solution().solute_moles(),
                baths[i].solution().solute_moles()
            );
        }
        near(fine.direct_water_balance_residual_kg, 0.0, 5e-11);
    });
}

fn kinetic_bath(volume: f64, initial: f64, saturation: f64, rate: f64) -> FiniteMoistureSolution {
    let law = IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("authored-neutral-solute").unwrap(),
        Temperature::new(293.15),
        AmountConcentration::new(saturation),
        0.0,
        [Temperature::new(293.15), Temperature::new(293.15)],
    )
    .unwrap();
    let kinetics = KineticSolubilitySolution::new(
        Temperature::new(293.15),
        WATER_VOLUME,
        Volume::new(volume),
        vec![(
            SolubilityKineticsV1::new(
                law,
                VolumetricFlowRate::new(rate),
                PrecipitationPolicyV1::ExistingSolidOnly,
            )
            .unwrap(),
            SolubilityKineticStateV1::from_phases(Amount::new(initial), Amount::new(0.2 - initial))
                .unwrap(),
        )],
    )
    .unwrap();
    FiniteMoistureSolution::with_kinetic_solubility(kinetics, WATER_MOLAR_MASS).unwrap()
}

#[test]
fn g1_each_kinetic_bath_resolves_actual_phases_against_independent_continuous_balances() {
    with_cx(|cx| {
        let volume = 1e-4;
        let rates = [2.5e-5, 1e-5];
        let saturation = [1000.0, 800.0];
        let initial_moles = [0.02, 0.04];
        let baths: [FiniteMoistureSolution; 2] = std::array::from_fn(|i| {
            kinetic_bath(volume, initial_moles[i], saturation[i], rates[i])
        });
        let field = slab(card([0.001; 2], 0.3), 0.1, vec![0.1], cx);
        let y0 = [0.1, 0.1, 0.1, initial_moles[0], initial_moles[1]];
        let oracle = reference(y0, 2.0, 4096, field.dry_mass_kg(), 0.2, rates, saturation);
        let check = reference(y0, 2.0, 2048, field.dry_mass_kg(), 0.2, rates, saturation);
        for i in 0..5 {
            near(oracle[i], check[i], 2e-13);
        }
        let loose = time(2.0, 1.0);
        let mut tight = loose;
        tight.kinetic_phases = [1e-7, 3e-7].map(|amount| {
            Some(KineticPhaseAccuracy {
                absolute_amount_tolerance_mol: amount,
                relative_tolerance: 0.0,
            })
        });
        let run = |controls| {
            baths[0]
                .advance_paired_adaptive_to(
                    &baths[1],
                    &field,
                    2.0,
                    &pair_faces(0.2),
                    config(),
                    physical([volume; 2]),
                    controls,
                    cx,
                )
                .unwrap()
        };
        let (coarse_baths, _, coarse) = run(loose);
        let (fine_baths, _, fine) = run(tight);
        let error = |b: &[FiniteMoistureSolution; 2]| {
            (0..2)
                .map(|i| {
                    (b[i].kinetic_solubility().unwrap().last_steps()[0]
                        .state()
                        .dissolved()
                        .value()
                        - oracle[3 + i])
                        .abs()
                })
                .fold(0.0_f64, f64::max)
        };
        assert!(error(&fine_baths) < 0.5 * error(&coarse_baths));
        assert!(fine.transport.rejected_intervals > 0);
        assert!(fine.field_evaluations > coarse.field_evaluations);
        for i in 0..2 {
            let before = baths[i].kinetic_solubility().unwrap();
            let after = fine_baths[i].kinetic_solubility().unwrap();
            assert_eq!(before.models(), after.models());
            let phase = after.last_steps()[0].state();
            assert_eq!(
                phase.total_solute(),
                before.last_steps()[0].state().total_solute()
            );
            near(
                phase.dissolved().value() + phase.solid().value(),
                0.2,
                1e-15,
            );
            assert!(phase.solid().value() > 0.0);
        }
        let mut invalid = tight;
        invalid.kinetic_phases[0]
            .as_mut()
            .unwrap()
            .absolute_amount_tolerance_mol = 0.0;
        assert!(matches!(
            baths[0].advance_paired_adaptive_to(
                &baths[1],
                &field,
                2.0,
                &pair_faces(0.2),
                config(),
                physical([volume; 2]),
                invalid,
                cx,
            ),
            Err(MoistureAdvanceError::Invalid(_))
        ));
        let dissolved = reservoir(volume, vec![0.02]);
        assert!(matches!(
            dissolved.advance_paired_adaptive_to(
                &baths[1],
                &field,
                2.0,
                &pair_faces(0.2),
                config(),
                physical([volume; 2]),
                tight,
                cx,
            ),
            Err(MoistureAdvanceError::Invalid(_))
        ));
    });
}

#[test]
fn g1_cartesian_and_curved_shell_durations_keep_their_actual_diffusion_geometry() {
    with_cx(|cx| {
        let volume = 0.01;
        let baths = [reservoir(volume, vec![0.0]), reservoir(volume, vec![0.0])];
        let d = 1e-4;
        let length = 0.1_f64;
        let m0 = 0.1;
        let grid = MoistureCartesian::new(
            card([d; 2], 0.3),
            point(0.1),
            ["moisture_diffusivity"; 3],
            MaterialPropertySelection::SingleClaimOnly,
            MoistureCartesianGeometry {
                lengths_m: [length; 3],
                cells: [1; 3],
                frame: ContentHash([4; 32]),
                maximum_moisture: 0.3,
                maximum_cells: 1,
            },
            vec![m0],
            cx,
        )
        .unwrap();
        let base_faces = pair_faces(0.2);
        let faces = std::array::from_fn(|i| base_faces[i % 2].clone());
        let (_, cart, report) = baths[0]
            .advance_paired_cartesian_adaptive_to(
                &baths[1],
                &grid,
                1.0,
                &faces,
                config(),
                nonlinear(),
                physical([volume; 2]),
                time(1.0, 1e-7),
                cx,
            )
            .unwrap();
        let exact = 0.2 - 0.1 * (-12.0 * d / length.powi(2)).exp();
        near(cart.moisture()[0], exact, 3e-5);
        near(
            report.solutions[0].water_to_solid_kg,
            report.solutions[1].water_to_solid_kg,
            2e-12,
        );
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.08 },
        ] {
            let a = 0.05_f64;
            let b = 0.1_f64;
            let center = 0.5 * (a + b);
            let field = MoistureRadial::new(
                card([d; 2], 0.3),
                point(0.1),
                "moisture_diffusivity",
                MaterialPropertySelection::SingleClaimOnly,
                MoistureRadialGeometry {
                    symmetry,
                    faces_m: vec![a, b],
                    maximum_moisture: 0.3,
                    maximum_cells: 1,
                },
                vec![m0],
                cx,
            )
            .unwrap();
            let conductance = match symmetry {
                RadialMoistureSymmetry::Sphere => {
                    500.0
                        * d
                        * 4.0
                        * core::f64::consts::PI
                        * (1.0 / (1.0 / a - 1.0 / center) + 1.0 / (1.0 / center - 1.0 / b))
                }
                RadialMoistureSymmetry::Cylinder { length_m } => {
                    500.0
                        * d
                        * 2.0
                        * core::f64::consts::PI
                        * length_m
                        * (1.0 / (center / a).ln() + 1.0 / (b / center).ln())
                }
            };
            let exact = 0.2 - 0.1 * (-conductance / field.dry_mass_kg()).exp();
            let (_, final_field, report) = baths[0]
                .advance_paired_radial_adaptive_to(
                    &baths[1],
                    &field,
                    1.0,
                    &base_faces,
                    config(),
                    nonlinear(),
                    physical([volume; 2]),
                    time(1.0, 1e-7),
                    cx,
                )
                .unwrap();
            near(final_field.moisture()[0], exact, 3e-5);
            near(report.direct_water_balance_residual_kg, 0.0, 5e-11);
            assert_eq!(final_field.time_s(), 1.0);
        }
    });
}

#[test]
fn g4_nested_work_includes_rejected_trials_and_late_failure_returns_no_prefix() {
    with_cx(|cx| {
        let volume = 1e-3;
        let baths = [reservoir(volume, vec![0.0]), reservoir(volume, vec![0.0])];
        let field = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let mut controls = time(1.0, 1.0);
        controls.maximum_field_evaluations = 5;
        assert!(matches!(
            baths[0].advance_paired_adaptive_to(
                &baths[1],
                &field,
                2.0,
                &[None, None],
                config(),
                physical([volume; 2]),
                controls,
                cx,
            ),
            Err(MoistureAdvanceError::Step {
                error: SolutionMoistureError::EvaluationBudget {
                    evaluations: 5,
                    maximum: 5
                },
                ..
            })
        ));
        controls = time(2.0, 1e-20);
        controls.maximum_field_evaluations = 3;
        assert!(matches!(
            baths[0].advance_paired_adaptive_to(
                &baths[1],
                &field,
                2.0,
                &pair_faces(0.2),
                config(),
                physical([volume; 2]),
                controls,
                cx,
            ),
            Err(MoistureAdvanceError::Step {
                error: SolutionMoistureError::EvaluationBudget {
                    evaluations: 3,
                    maximum: 3
                },
                ..
            })
        ));
        assert_eq!(field.moisture(), &[0.1]);
        assert_eq!(field.time_s(), 0.0);
        for bath in &baths {
            assert_eq!(bath.solution().solvent_volume_m3(), volume);
        }
    });
}

#[test]
fn g4_cancellation_preserves_both_original_solutions_and_field() {
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
        let volume = 1e-4;
        let baths = [reservoir(volume, vec![0.02]), reservoir(volume, vec![0.04])];
        let field = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], &cx);
        gate.request();
        assert!(matches!(
            baths[0].advance_paired_adaptive_to(
                &baths[1],
                &field,
                1.0,
                &pair_faces(0.2),
                config(),
                physical([volume; 2]),
                time(1.0, 1e-5),
                &cx,
            ),
            Err(MoistureAdvanceError::Cancelled)
        ));
        assert_eq!(field.time_s(), 0.0);
        assert_eq!(field.moisture(), &[0.1]);
        for bath in &baths {
            assert_eq!(bath.solution().solvent_volume_m3(), volume);
        }
    });
}
