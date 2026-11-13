//! G1 temporal convergence of actual geometry; G3/G4 retained accounts and work.
//! Authored numerical sources, not a material calibration or global error proof.

use super::*;
use fs_couple::moisture::adaptive::{
    KineticPhaseAccuracy, MoistureAdvanceError, MoistureSpatialPairedSolutionTimeControls,
    MoistureTimeControls, PairedSolutionMoistureTimeControls,
};

fn time(span: f64, displacement: f64) -> MoistureSpatialPairedSolutionTimeControls {
    MoistureSpatialPairedSolutionTimeControls {
        solutions: PairedSolutionMoistureTimeControls {
            transport: MoistureTimeControls {
                initial_step_s: span,
                minimum_step_s: 1e-6,
                maximum_step_s: span,
                maximum_evaluations: 100_000,
                moisture_absolute_tolerance: 1.0,
                activity_absolute_tolerance: 1.0,
                reservoir_water_absolute_tolerance_kg: 1.0,
                relative_tolerance: 0.0,
                water_balance_tolerance_kg: 5e-11,
            },
            kinetic_phases: [None; 2],
            maximum_field_evaluations: 100_000,
        },
        displacement_absolute_tolerance_m: displacement,
    }
}

fn physical(volumes: [f64; 2]) -> MoistureSpatialPairedSolutionControls {
    let mut controls = coupling(volumes);
    controls.maximum_native_solid_solves = 100_000;
    for bath in &mut controls.solutions.solutions {
        bath.mass_tolerance_kg = 1e-15;
    }
    controls
}

#[test]
fn g1_displacement_alone_controls_swelling_accuracy_and_all_discarded_work_is_charged() {
    with_cx(|cx| {
        let (beta, duration, volume) = (0.1, 0.8, 1e-3);
        let stack = one_cell([beta, 0.0, 0.0], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let baths = [
            aqueous::reservoir(volume, vec![0.0]),
            aqueous::reservoir(volume, vec![0.0]),
        ];
        // Independently integrate dM/dt = 4(.3-M)/(1+beta*(M-.1))².
        // Separation gives this monotone logarithmic primitive. It uses no
        // native diffusion step, solid solve or bath-root implementation.
        let a = 1.0 + 0.2 * beta;
        let exact = root(0.1, 0.3 - 1e-12, |m| {
            let y = 0.3 - m;
            a * a * (0.2 / y).ln() - 2.0 * a * beta * (0.2 - y) + 0.5 * beta * beta * (0.04 - y * y)
                - 4.0 * duration
        });
        let run = |accuracy, limits| {
            solid.advance_spatial_diffusion_paired_solutions_adaptive_to(
                [&baths[0], &baths[1]],
                &stack,
                duration,
                pair_faces(),
                config(),
                nonlinear(),
                limits,
                accuracy,
                cx,
            )
        };
        let mut errors = Vec::new();
        let mut counts = Vec::new();
        for displacement in [1e-5, 2.5e-6, 6.25e-7] {
            let accuracy = time(duration, displacement);
            let (next, result) = run(accuracy, physical([volume; 2])).unwrap();
            let m = result.moisture.layers()[0].slab().moisture()[0];
            let report = &result.transport;
            let summary = &report.transport.summary;
            errors.push((m - exact).abs());
            counts.push(summary.evaluations);
            assert!(summary.rejected_intervals > 0);
            assert!(summary.maximum_accepted_error_ratio <= 1.0);
            assert_eq!(report.field_evaluations, summary.evaluations);
            assert_eq!(result.moisture.time_s(), duration);
            assert_eq!(result.solid.time_s(), duration);
            assert_eq!(result.moisture.dry_mass_kg(), stack.dry_mass_kg());
            assert_eq!(
                result.solid.geometry().identity(),
                summary.last_substep.transport.spatial.geometry_identity
            );
            assert!(result.solid.tangent().is_some());
            near(result.solid.wet_mass_kg(), 0.05 * (1.0 + m), 2e-12);
            for (x, u) in solid
                .reference_nodes_m()
                .iter()
                .zip(result.solid.solution().displacement_m())
            {
                near(u[0], beta * (m - 0.1) * x[0], 2e-10);
                near(u[1], 0.0, 2e-10);
                near(u[2], 0.0, 2e-10);
            }
            for port in 0..2 {
                let received = 0.025 * (m - 0.1);
                near(report.solutions[port].water_to_solid_kg, received, 5e-12);
                near(
                    next[port].solvent_mass_kg(),
                    WATER_DENSITY * volume - received,
                    5e-12,
                );
                near(report.solutions[port].water_balance_residual_kg, 0.0, 5e-11);
            }
            near(report.direct_water_balance_residual_kg, 0.0, 5e-11);
            if displacement == 6.25e-7 {
                // These caps exceed any single endpoint's work but fall one
                // below the complete run, which includes rejected comparisons.
                // A reset per trial or per accepted interval would wrongly pass.
                let mut short = accuracy;
                short.solutions.maximum_field_evaluations = report.field_evaluations - 1;
                assert!(matches!(
                    run(short, physical([volume;2])),
                    Err(MoistureAdvanceError::Step {
                        error: SolutionMoistureError::EvaluationBudget { evaluations, maximum }, ..
                    }) if evaluations == maximum && maximum == report.field_evaluations-1
                ));
                let mut limited = physical([volume; 2]);
                limited.maximum_native_solid_solves = report.native_solid_solves - 1;
                assert!(matches!(
                    run(accuracy, limited),
                    Err(MoistureAdvanceError::Step {
                        error: SolutionMoistureError::Moisture(MoistureError::Invalid(message)), ..
                    }) if message.contains("solid-solve budget")
                ));
            }
        }
        assert!(
            errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
            "{errors:?}"
        );
        assert!(counts[0] < counts[1] && counts[1] < counts[2], "{counts:?}");
        assert_eq!(stack.time_s(), 0.0);
        assert_eq!(stack.layers()[0].slab().moisture(), &[0.1]);
        for bath in &baths {
            assert_eq!(bath.solution().solvent_volume_m3(), volume);
        }
    });
}

#[test]
fn g3_only_original_fine_endpoints_supply_shape_and_opposing_bath_accounts() {
    with_cx(|cx| {
        let (volume, duration) = (2e-5, 0.25);
        let activity = [0.997, 0.99];
        let stack = MoistureLayeredSlab::new(
            vec![layer(
                material(500.0, [0.08, 0.0, 0.0], 0.3, 0.3),
                0.1,
                vec![0.15 * (activity[0] + activity[1])],
                cx,
            )],
            16,
            cx,
        )
        .unwrap();
        let solid = bind_stack(&stack, free_supports(), cx);
        let baths = activity.map(|a| dissolved_at(volume, a));
        let limits = physical([volume; 2]);
        let step = |b: [&FiniteMoistureSolution; 2], field: &MoistureLayeredSlab, dt| {
            solid
                .advance_spatial_diffusion_paired_solutions(
                    b,
                    field,
                    dt,
                    pair_faces(),
                    config(),
                    nonlinear(),
                    limits,
                    cx,
                )
                .unwrap()
        };
        let (coarse, _) = step([&baths[0], &baths[1]], &stack, duration);
        let (first, _) = step([&baths[0], &baths[1]], &stack, duration / 2.0);
        let (fine, native) = step(
            [&first.solutions()[0], &first.solutions()[1]],
            first.moisture(),
            duration / 2.0,
        );
        let (next, result) = solid
            .advance_spatial_diffusion_paired_solutions_adaptive_to(
                [&baths[0], &baths[1]],
                &stack,
                duration,
                pair_faces(),
                config(),
                nonlinear(),
                limits,
                time(duration, 1.0),
                cx,
            )
            .unwrap();
        let report = &result.transport;
        let summary = &report.transport.summary;
        assert_eq!(summary.evaluations, 3);
        assert_eq!(summary.accepted_substeps, 2);
        assert_eq!(summary.rejected_intervals, 0);
        assert_eq!(
            report.field_evaluations,
            coarse.step().field_evaluations
                + first.step().field_evaluations
                + fine.step().field_evaluations
        );
        assert_eq!(
            report.native_solid_solves,
            1 + coarse.step().transport.native_solid_solves
                + first.step().transport.native_solid_solves
                + fine.step().transport.native_solid_solves
        );
        assert_eq!(
            result.moisture.layers()[0].slab().moisture(),
            fine.moisture().layers()[0].slab().moisture()
        );
        assert_eq!(
            result.solid.solution().displacement_m(),
            native.solution().displacement_m()
        );
        assert_eq!(
            result.solid.geometry().identity(),
            native.geometry().identity()
        );
        assert!(report.solutions[0].water_to_solid_kg > 0.0);
        assert!(report.solutions[1].water_to_solid_kg < 0.0);
        assert!(summary.water_change_kg.abs() < 0.1 * report.solutions[0].water_to_solid_kg.abs());
        for port in 0..2 {
            assert_eq!(
                next[port].solution().solvent_volume_m3(),
                fine.solutions()[port].solution().solvent_volume_m3()
            );
            assert_eq!(
                next[port].solution().solute_moles(),
                baths[port].solution().solute_moles()
            );
            near(
                report.solutions[port].water_to_solid_kg,
                first.step().solutions[port].water_to_solid_kg
                    + fine.step().solutions[port].water_to_solid_kg,
                1e-15,
            );
            near(report.solutions[port].water_balance_residual_kg, 0.0, 5e-11);
        }
        near(report.direct_water_balance_residual_kg, 0.0, 5e-11);
    });
}

#[test]
fn g1_sealed_layer_redistribution_and_two_original_kinetic_phase_histories_survive_adaptation() {
    with_cx(|cx| {
        let volumes = [1e-4, 2e-4];
        let rates = [2e-5, 5e-5];
        let stocks = [0.2, 0.3];
        let duration = 0.5;
        let baths: [FiniteMoistureSolution; 2] = std::array::from_fn(|i| {
            let kinetic = KineticSolubilitySolution::new(
                Temperature::new(293.15),
                aqueous::WATER_VOLUME,
                Volume::new(volumes[i]),
                vec![(
                    SolubilityKineticsV1::new(
                        solubility(1000.0),
                        VolumetricFlowRate::new(rates[i]),
                        PrecipitationPolicyV1::ExistingSolidOnly,
                    )
                    .unwrap(),
                    SolubilityKineticStateV1::from_phases(Amount::new(0.0), Amount::new(stocks[i]))
                        .unwrap(),
                )],
            )
            .unwrap();
            FiniteMoistureSolution::with_kinetic_solubility(kinetic, 0.018).unwrap()
        });
        let stack = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(500.0, [0.08, 0.0, 0.0], 0.3, 0.3),
                    0.1,
                    vec![0.2],
                    cx,
                ),
                layer(
                    material(800.0, [0.1, 0.0, 0.0], 0.3, 0.2),
                    0.2,
                    vec![0.1],
                    cx,
                ),
            ],
            16,
            cx,
        )
        .unwrap()
        .with_transfer(
            [MoistureActivityTransfer::Sealed; 2],
            vec![MoistureActivityTransfer::IdealContact],
            cx,
        )
        .unwrap();
        let solid = bind_stack(&stack, free_supports(), cx);
        let mut accuracy = time(duration, 1.0);
        // Local water differences cannot cancel between the unequal layers.
        // Everything else is loose; both sealed phase evolutions are exact
        // fixed-volume exponentials and retain independent original stocks.
        accuracy
            .solutions
            .transport
            .reservoir_water_absolute_tolerance_kg = 1e-7;
        accuracy.solutions.kinetic_phases = [Some(KineticPhaseAccuracy {
            absolute_amount_tolerance_mol: 1e-9,
            relative_tolerance: 0.0,
        }); 2];
        let (next, result) = solid
            .advance_spatial_diffusion_paired_solutions_adaptive_to(
                [&baths[0], &baths[1]],
                &stack,
                duration,
                pair_faces(),
                config(),
                nonlinear(),
                physical(volumes),
                accuracy,
                cx,
            )
            .unwrap();
        let report = &result.transport;
        assert!(report.transport.summary.rejected_intervals > 0);
        assert_eq!(report.transport.summary.boundary_water_kg, [0.0; 2]);
        assert_eq!(
            report.transport.summary.last_substep.face_solutions,
            [None; 2]
        );
        assert!(report.transport.layers[0].water_change_kg < 0.0);
        assert!(report.transport.layers[1].water_change_kg > 0.0);
        near(
            report.transport.layers[0].boundary_water_kg[1],
            -report.transport.layers[1].boundary_water_kg[0],
            1e-15,
        );
        near(
            result.moisture.water_mass_kg(),
            stack.water_mass_kg(),
            5e-11,
        );
        near(result.solid.wet_mass_kg(), stack.wet_mass_kg(), 5e-11);
        for i in 0..2 {
            assert_eq!(next[i].solvent_mass_kg(), baths[i].solvent_mass_kg());
            assert_eq!(report.solutions[i].water_to_solid_kg, 0.0);
            let phase = next[i].kinetic_solubility().unwrap().last_steps()[0].state();
            let dissolved = -1000.0 * volumes[i] * (-rates[i] * duration / volumes[i]).exp_m1();
            near(phase.dissolved().value(), dissolved, 2e-13);
            near(phase.solid().value(), stocks[i] - dissolved, 2e-13);
            assert!(phase.solid().value() > 0.0);
            near(report.solutions[i].water_balance_residual_kg, 0.0, 5e-11);
        }
        assert_eq!(result.solid.time_s(), duration);
        assert!(result.solid.tangent().is_some());
        near(report.direct_water_balance_residual_kg, 0.0, 5e-11);
    });
}
