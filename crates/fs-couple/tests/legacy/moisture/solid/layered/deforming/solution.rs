//! G1 native swelling/bath balances and G3/G4 history and atomic refusal.
//! All coefficients are authored numerical fixtures, not material calibration.

use super::*;
#[path = "solution/adaptive.rs"]
mod adaptive;
use crate::solution as aqueous;
use fs_couple::moisture::solution::{
    FiniteMoistureSolution, MoistureSolutionPort, PairedSolutionMoistureSolve,
    SolutionMoistureError,
};
use fs_material::osmosis::solubility::{
    FiniteSolubilitySolution, kinetics::KineticSolubilitySolution,
};
use fs_qty::{Amount, AmountConcentration, Temperature, Time, Volume, VolumetricFlowRate};
use fs_thermochem::{
    SpeciesId,
    solution::solubility::{
        IdealSolubilityV1,
        kinetics::{PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1},
    },
};

const WATER_DENSITY: f64 = 0.018 / aqueous::WATER_VOLUME;

fn coupling(volumes: [f64; 2]) -> MoistureSpatialPairedSolutionControls {
    MoistureSpatialPairedSolutionControls {
        solutions: PairedSolutionMoistureSolve {
            solutions: volumes.map(aqueous::controls),
            maximum_evaluations: 1000,
            mass_tolerance_kg: 2e-11,
        },
        maximum_geometry_residual_m: 1e-9,
        maximum_native_solid_solves: 20_000,
    }
}

fn pair_faces() -> [Option<MoistureSolutionPort>; 2] {
    [
        Some(MoistureSolutionPort::First),
        Some(MoistureSolutionPort::Second),
    ]
}

fn dissolved_at(volume: f64, activity: f64) -> FiniteMoistureSolution {
    aqueous::reservoir(
        volume,
        vec![-activity.ln() * volume / aqueous::WATER_VOLUME],
    )
}

fn solubility(saturation: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("authored-neutral-solute").unwrap(),
        Temperature::new(293.15),
        AmountConcentration::new(saturation),
        0.0,
        [Temperature::new(293.15), Temperature::new(293.15)],
    )
    .unwrap()
}

#[test]
fn g1_actual_bath_composition_area_and_length_close_one_native_solid_endpoint() {
    with_cx(|cx| {
        let beta = [0.1, 0.08, 0.04];
        let initial_m = 0.295;
        let target = 0.2976;
        let activity = [0.997_f64, 0.99];
        let final_volume = [1e-4, 2e-4];
        let h = 0.4;
        let stack = MoistureLayeredSlab::new(
            vec![layer(
                material(500.0, beta, 0.3, 0.3),
                0.1,
                vec![initial_m],
                cx,
            )],
            16,
            cx,
        )
        .unwrap()
        .with_transfer(
            [MoistureActivityTransfer::ConductanceKgM2S(h); 2],
            vec![],
            cx,
        )
        .unwrap();
        let solid = bind_stack(&stack, free_supports(), cx);
        let before = solid.resolve(&stack, cx).unwrap();
        // Manufacture the endpoint from independent current-length and current-
        // area resistance equations, then infer each ORIGINAL finite inventory.
        let stretch = beta.map(|b| 1.0 + b * (target - 0.1));
        let area = 0.001 * stretch[1] * stretch[2];
        let half = 2.0 * 0.05 * 0.01 * 0.3 / (0.1_f64 * stretch[0]).powi(2);
        let g = 1.0 / (1.0 / half + 1.0 / (h * area));
        let flux = activity.map(|a| g * (a - target / 0.3));
        let dt = 0.05 * (target - initial_m) / (flux[0] + flux[1]);
        let transfer = flux.map(|q| q * dt);
        assert!(transfer[0] > 0.0 && transfer[1] < 0.0);
        let volumes: [f64; 2] =
            std::array::from_fn(|i| final_volume[i] + transfer[i] / WATER_DENSITY);
        let baths: [FiniteMoistureSolution; 2] = std::array::from_fn(|i| {
            aqueous::reservoir(
                volumes[i],
                vec![-activity[i].ln() * final_volume[i] / aqueous::WATER_VOLUME],
            )
        });
        let (candidate, state) = solid
            .advance_spatial_diffusion_paired_solutions(
                [&baths[0], &baths[1]],
                &stack,
                dt,
                pair_faces(),
                config(),
                nonlinear(),
                coupling(volumes),
                cx,
            )
            .unwrap();
        let actual = candidate.moisture().layers()[0].slab().moisture()[0];
        near(actual, target, 5e-10);
        let step = candidate.step();
        let spatial = &step.transport.spatial;
        near(spatial.axial_stretch[0], stretch[0], 5e-10);
        near(spatial.current_cross_section_m2, area, 2e-12);
        for i in 0..2 {
            near(
                candidate.solutions()[i].solution().solvent_volume_m3(),
                final_volume[i],
                5e-15,
            );
            near(
                candidate.solutions()[i]
                    .solution()
                    .water_activity()
                    .unwrap(),
                activity[i],
                5e-11,
            );
            near(step.solutions[i].water_to_solid_kg, transfer[i], 5e-12);
            near(step.solutions[i].water_balance_residual_kg, 0.0, 2e-11);
            assert_eq!(
                candidate.solutions()[i].solution().solute_moles(),
                baths[i].solution().solute_moles()
            );
            let surface = spatial.transport.surface_activities[i].unwrap();
            near(
                spatial.transport.boundary_water_kg[i],
                dt * h * area * (activity[i] - surface),
                5e-12,
            );
        }
        near(step.water_change_kg, 0.05 * (target - initial_m), 5e-12);
        near(step.direct_water_balance_residual_kg, 0.0, 2e-11);
        assert_eq!(candidate.moisture().dry_mass_kg(), stack.dry_mass_kg());
        near(state.wet_mass_kg(), 0.05 * (1.0 + target), 5e-12);
        assert_eq!(state.time_s(), candidate.moisture().time_s());
        assert_eq!(state.geometry().identity(), spatial.geometry_identity);
        for (x, u) in solid
            .reference_nodes_m()
            .iter()
            .zip(state.solution().displacement_m())
        {
            for axis in 0..3 {
                near(u[axis], beta[axis] * (target - 0.1) * x[axis], 2e-10);
            }
        }
        // Uniform free swelling has zero prestress. The actual source E(M) and
        // wet reference mass therefore scale every generalized modal eigenvalue
        // by E_after/E_before * (1+M_before)/(1+M_after). Inspect a native pencil
        // Rayleigh quotient, not a separately assembled oscillator.
        let quotient = |state: &MoistureSolidState| {
            let pencil = state.tangent().unwrap().pencil();
            pencil.stiffness.get(0, 0) / pencil.mass.get(0, 0)
        };
        let ratio = (1.0 - target) / (1.0 - initial_m) * (1.0 + initial_m) / (1.0 + target);
        near(quotient(&state) / quotient(&before), ratio, 1e-7);
        assert!(quotient(&state) < quotient(&before));
        // Hold the same accepted bath activities to isolate geometry feedback:
        // the ordinary reference-metric transport reaches a different endpoint.
        let fixed = stack
            .advance(
                dt,
                &activity.map(MoistureLayeredFace::PrescribedRelativeHumidity),
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        assert!((fixed.0.layers()[0].slab().moisture()[0] - actual).abs() > 1e-6);
        assert!(step.field_evaluations > 1);
        assert!(step.transport.native_solid_solves > spatial.native_solid_solves);
        assert_eq!(stack.layers()[0].slab().moisture(), &[initial_m]);
    });
}

#[test]
fn g1_unlike_stretched_layers_pass_water_between_independent_finite_crystal_baths() {
    with_cx(|cx| {
        let source_activity = [0.997, 0.99];
        let slopes = [0.3, 0.2];
        let m = [
            slopes[0] * source_activity[0],
            slopes[1] * source_activity[1],
        ];
        let beta = [0.1, 0.2];
        let lengths = [0.1_f64, 0.2];
        let dry = [0.05, 0.16];
        let stretch: [f64; 2] = std::array::from_fn(|i| 1.0 + beta[i] * (m[i] - 0.1));
        let halves: [f64; 2] = std::array::from_fn(|i| {
            2.0 * dry[i] * 0.01 * slopes[i] / (lengths[i] * stretch[i]).powi(2)
        });
        let contact_h = 0.5;
        let flux = (source_activity[0] - source_activity[1])
            / (1.0 / halves[0] + 1.0 / (contact_h * 0.001) + 1.0 / halves[1]);
        let bath_activity = [
            source_activity[0] + flux / halves[0],
            source_activity[1] - flux / halves[1],
        ];
        let dt = 1.0;
        let final_volume = [1e-4; 2];
        let transfer = [dt * flux, -dt * flux];
        let volumes: [f64; 2] =
            std::array::from_fn(|i| final_volume[i] + transfer[i] / WATER_DENSITY);
        let saturation = bath_activity.map(|a| -a.ln() / aqueous::WATER_VOLUME);
        let totals: [f64; 2] =
            std::array::from_fn(|i| 2.0 * saturation[i] * volumes[i].max(final_volume[i]));
        let baths: [FiniteMoistureSolution; 2] = std::array::from_fn(|i| {
            FiniteMoistureSolution::with_solubility(
                FiniteSolubilitySolution::new(
                    Temperature::new(293.15),
                    aqueous::WATER_VOLUME,
                    Volume::new(volumes[i]),
                    vec![(solubility(saturation[i]), Amount::new(totals[i]))],
                )
                .unwrap(),
                0.018,
            )
            .unwrap()
        });
        let stack = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(500.0, [beta[0], 0.0, 0.0], 0.3, slopes[0]),
                    lengths[0],
                    vec![m[0]],
                    cx,
                ),
                layer(
                    material(800.0, [beta[1], 0.0, 0.0], 0.3, slopes[1]),
                    lengths[1],
                    vec![m[1]],
                    cx,
                ),
            ],
            16,
            cx,
        )
        .unwrap()
        .with_transfer(
            [MoistureActivityTransfer::IdealContact; 2],
            vec![MoistureActivityTransfer::ConductanceKgM2S(contact_h)],
            cx,
        )
        .unwrap();
        let solid = bind_stack(&stack, free_supports(), cx);
        let (candidate, state) = solid
            .advance_spatial_diffusion_paired_solutions(
                [&baths[0], &baths[1]],
                &stack,
                dt,
                pair_faces(),
                config(),
                nonlinear(),
                coupling(volumes),
                cx,
            )
            .unwrap();
        let step = candidate.step();
        near(step.water_change_kg, 0.0, 5e-12);
        near(
            step.transport.spatial.transport.interfaces[0].water_left_to_right_kg,
            dt * flux,
            5e-12,
        );
        for i in 0..2 {
            near(
                candidate.moisture().layers()[i].slab().moisture()[0],
                m[i],
                5e-10,
            );
            near(step.transport.spatial.axial_stretch[i], stretch[i], 5e-10);
            near(step.solutions[i].water_to_solid_kg, transfer[i], 5e-12);
            near(
                candidate.solutions()[i].solution().solvent_volume_m3(),
                final_volume[i],
                5e-15,
            );
            let original = baths[i].solubility().unwrap();
            let final_bath = candidate.solutions()[i].solubility().unwrap();
            assert_eq!(original.components(), final_bath.components());
            let phase = final_bath.phase_partitions()[0];
            near(
                phase.dissolved.value(),
                saturation[i] * final_volume[i],
                5e-12,
            );
            near(
                phase.dissolved.value() + phase.solid.value(),
                totals[i],
                1e-14,
            );
            assert!(phase.solid.value() > 0.0);
            assert_eq!(
                state.elastic_states()[i].card_identity(),
                stack.layers()[i].slab().card().content_hash()
            );
        }
        near(state.wet_mass_kg(), stack.wet_mass_kg(), 5e-12);
        near(step.direct_water_balance_residual_kg, 0.0, 2e-11);
    });
}

#[test]
fn g3_zero_swelling_reproduces_original_finite_solution_transport() {
    with_cx(|cx| {
        let volume = 1e-4;
        let baths = [dissolved_at(volume, 0.99), dissolved_at(volume, 0.98)];
        let stack = one_cell([0.0; 3], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let original = baths[0]
            .advance_layered(
                &stack,
                0.25,
                [true; 2],
                config(),
                nonlinear(),
                aqueous::controls(volume),
                cx,
            )
            .unwrap();
        let (candidate, state) = solid
            .advance_spatial_diffusion_paired_solutions(
                [&baths[0], &baths[1]],
                &stack,
                0.25,
                [Some(MoistureSolutionPort::First); 2],
                config(),
                nonlinear(),
                coupling([volume; 2]),
                cx,
            )
            .unwrap();
        near(
            candidate.moisture().activity()[0],
            original.moisture().activity()[0],
            2e-12,
        );
        for i in 0..2 {
            near(
                candidate
                    .step()
                    .transport
                    .spatial
                    .transport
                    .boundary_water_kg[i],
                original.step().transport.boundary_water_kg[i],
                2e-12,
            );
        }
        near(
            candidate.solutions()[0].solution().solvent_volume_m3(),
            original.reservoir().solution().solvent_volume_m3(),
            2e-15,
        );
        // At the identical accepted activity, the zero-swelling metric is
        // exactly one and the original native transport arithmetic is retained.
        let activity = candidate.solutions()[0]
            .solution()
            .water_activity()
            .unwrap();
        let fixed = stack
            .advance(
                0.25,
                &[MoistureLayeredFace::PrescribedRelativeHumidity(activity); 2],
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.moisture().activity(), fixed.0.activity());
        assert_eq!(
            candidate
                .step()
                .transport
                .spatial
                .transport
                .boundary_water_kg,
            fixed.1.boundary_water_kg
        );
        assert_eq!(
            candidate.solutions()[1].solution().solvent_volume_m3(),
            volume
        );
        assert_eq!(candidate.step().solutions[1].water_to_solid_kg, 0.0);
        assert_eq!(candidate.step().transport.spatial.axial_stretch, vec![1.0]);
        near(state.stress().elastic_energy_j(), 0.0, 1e-14);
        near(
            state.wet_mass_kg(),
            candidate.moisture().wet_mass_kg(),
            1e-13,
        );
    });
}

#[test]
fn g3_actual_wet_dry_bath_history_changes_shape_without_resetting_either_inventory() {
    with_cx(|cx| {
        let volume = 5e-4;
        let beta = 0.05;
        let mut baths = [dissolved_at(volume, 0.9999), dissolved_at(volume, 0.98)];
        let mut stack = MoistureLayeredSlab::new(
            vec![layer(
                material(500.0, [beta, 0.0, 0.0], 0.3, 0.3),
                0.1,
                vec![0.297],
                cx,
            )],
            16,
            cx,
        )
        .unwrap();
        let solid = bind_stack(&stack, free_supports(), cx);
        let initial_water = baths
            .iter()
            .map(FiniteMoistureSolution::solvent_mass_kg)
            .sum::<f64>()
            + stack.water_mass_kg();
        let stocks = baths
            .each_ref()
            .map(|b| b.solution().solute_moles().to_vec());
        for cycle in 0..4 {
            let side = cycle % 2;
            let old_m = stack.layers()[0].slab().moisture()[0];
            let ports = [Some(if side == 0 {
                MoistureSolutionPort::First
            } else {
                MoistureSolutionPort::Second
            }); 2];
            let (candidate, state) = solid
                .advance_spatial_diffusion_paired_solutions(
                    [&baths[0], &baths[1]],
                    &stack,
                    0.25,
                    ports,
                    config(),
                    nonlinear(),
                    coupling([volume; 2]),
                    cx,
                )
                .unwrap();
            let (next_baths, next, step) = candidate.into_parts();
            let m = next.layers()[0].slab().moisture()[0];
            assert_eq!(m > old_m, side == 0);
            near(
                step.transport.spatial.axial_stretch[0],
                1.0 + beta * (m - 0.1),
                2e-10,
            );
            assert_eq!(step.solutions[1 - side].water_to_solid_kg, 0.0);
            near(state.wet_mass_kg(), next.wet_mass_kg(), 2e-12);
            assert!(state.tangent().is_some());
            baths = next_baths;
            stack = next;
            near(
                baths
                    .iter()
                    .map(FiniteMoistureSolution::solvent_mass_kg)
                    .sum::<f64>()
                    + stack.water_mass_kg(),
                initial_water,
                5e-11,
            );
            for i in 0..2 {
                assert_eq!(baths[i].solution().solute_moles(), stocks[i]);
            }
            assert_eq!(stack.time_s(), 0.25 * (cycle + 1) as f64);
        }
    });
}

#[test]
fn g1_sealed_baths_retain_water_while_phase_kinetics_and_spatial_redistribution_continue() {
    with_cx(|cx| {
        let volume = 1e-4;
        let kinetic = KineticSolubilitySolution::new(
            Temperature::new(293.15),
            aqueous::WATER_VOLUME,
            Volume::new(volume),
            vec![(
                SolubilityKineticsV1::new(
                    solubility(1000.0),
                    VolumetricFlowRate::new(2e-5),
                    PrecipitationPolicyV1::ExistingSolidOnly,
                )
                .unwrap(),
                SolubilityKineticStateV1::from_phases(Amount::new(0.0), Amount::new(0.2)).unwrap(),
            )],
        )
        .unwrap();
        let mut expected = kinetic.clone();
        expected.advance(Time::new(0.2), || false).unwrap();
        let baths = [
            dissolved_at(volume, 0.99),
            FiniteMoistureSolution::with_kinetic_solubility(kinetic, 0.018).unwrap(),
        ];
        let stack = MoistureLayeredSlab::new(
            vec![
                layer(
                    material(500.0, [0.1, 0.0, 0.0], 0.3, 0.3),
                    0.1,
                    vec![0.2],
                    cx,
                ),
                layer(
                    material(800.0, [0.2, 0.0, 0.0], 0.3, 0.2),
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
        let (candidate, state) = solid
            .advance_spatial_diffusion_paired_solutions(
                [&baths[0], &baths[1]],
                &stack,
                0.2,
                pair_faces(),
                config(),
                nonlinear(),
                coupling([volume; 2]),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.step().face_solutions, [None; 2]);
        assert_eq!(candidate.step().field_evaluations, 1);
        assert_eq!(
            candidate
                .step()
                .transport
                .spatial
                .transport
                .boundary_water_kg,
            [0.0; 2]
        );
        assert!(candidate.moisture().layers()[0].slab().moisture()[0] < 0.2);
        assert!(candidate.moisture().layers()[1].slab().moisture()[0] > 0.1);
        near(state.wet_mass_kg(), stack.wet_mass_kg(), 2e-12);
        let actual = candidate.solutions()[1]
            .kinetic_solubility()
            .unwrap()
            .last_steps()[0]
            .state();
        let expected = expected.last_steps()[0].state();
        assert_eq!(actual.dissolved(), expected.dissolved());
        assert_eq!(actual.solid(), expected.solid());
        for bath in candidate.solutions() {
            assert_eq!(bath.solution().solvent_volume_m3(), volume);
        }
    });
}

#[test]
fn g4_cumulative_work_late_source_and_geometry_refusals_preserve_every_original_owner() {
    with_cx(|cx| {
        let volume = 1e-4;
        let baths = [dissolved_at(volume, 0.99), dissolved_at(volume, 0.98)];
        let stack = one_cell([0.1, 0.0, 0.0], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let run = |controls| {
            solid.advance_spatial_diffusion_paired_solutions(
                [&baths[0], &baths[1]],
                &stack,
                0.25,
                pair_faces(),
                config(),
                nonlinear(),
                controls,
                cx,
            )
        };
        let (accepted, _) = run(coupling([volume; 2])).unwrap();
        let mut bounded = coupling([volume; 2]);
        bounded.solutions.maximum_evaluations = 1;
        assert!(matches!(
            run(bounded),
            Err(SolutionMoistureError::EvaluationBudget {
                evaluations: 1,
                maximum: 1
            })
        ));
        bounded = coupling([volume; 2]);
        bounded.maximum_native_solid_solves = accepted.step().transport.native_solid_solves - 1;
        assert!(
            matches!(run(bounded), Err(SolutionMoistureError::Moisture(MoistureError::Invalid(message))) if message.contains("solid-solve budget"))
        );
        let restricted = MoistureLayeredSlab::new(
            vec![layer(
                material(500.0, [0.1, 0.0, 0.0], 0.1001, 0.3),
                0.1,
                vec![0.1],
                cx,
            )],
            16,
            cx,
        )
        .unwrap();
        assert!(
            bind_stack(&restricted, free_supports(), cx)
                .advance_spatial_diffusion_paired_solutions(
                    [&baths[0], &baths[1]],
                    &restricted,
                    0.25,
                    pair_faces(),
                    config(),
                    nonlinear(),
                    coupling([volume; 2]),
                    cx,
                )
                .is_err()
        );
        let warped = one_cell([0.01, 0.08, 0.04], cx);
        assert!(
            matches!(bind_stack(&warped, (0..12).collect(), cx).advance_spatial_diffusion_paired_solutions(
            [&baths[0], &baths[1]], &warped, 0.25, pair_faces(), config(), nonlinear(), coupling([volume; 2]), cx,
        ), Err(SolutionMoistureError::Moisture(MoistureError::Invalid(message))) if message.contains("axis-aligned"))
        );
        let (retry, _) = run(coupling([volume; 2])).unwrap();
        assert_eq!(retry.moisture().activity(), accepted.moisture().activity());
        assert_eq!(stack.time_s(), 0.0);
        assert_eq!(stack.layers()[0].slab().moisture(), &[0.1]);
        for bath in &baths {
            assert_eq!(bath.solution().solvent_volume_m3(), volume);
        }
    });
}

#[test]
fn g4_cancelled_spatial_pair_keeps_baths_water_and_solid_at_their_start() {
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
        let baths = [dissolved_at(volume, 0.99), dissolved_at(volume, 0.98)];
        let stack = one_cell([0.1, 0.0, 0.0], &cx);
        let solid = bind_stack(&stack, free_supports(), &cx);
        gate.request();
        assert!(matches!(
            solid.advance_spatial_diffusion_paired_solutions(
                [&baths[0], &baths[1]],
                &stack,
                0.25,
                pair_faces(),
                config(),
                nonlinear(),
                coupling([volume; 2]),
                &cx,
            ),
            Err(SolutionMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(stack.time_s(), 0.0);
        assert_eq!(stack.layers()[0].slab().moisture(), &[0.1]);
        for bath in &baths {
            assert_eq!(bath.solution().solvent_volume_m3(), volume);
        }
    });
}
