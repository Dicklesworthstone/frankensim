//! Original kinetic phases on the native slab/Cartesian/shared water root.

use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::cartesian::{MoistureCartesian, MoistureCartesianGeometry};
use fs_material::osmosis::solubility::kinetics::KineticSolubilitySolution;
use fs_qty::VolumetricFlowRate;
use fs_thermochem::solution::solubility::kinetics::{
    PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1,
};

#[path = "kinetics/adaptive.rs"]
mod adaptive;

fn kinetic_bath(
    volume: f64,
    dissolved: f64,
    solid: f64,
    conductance: f64,
    policy: PrecipitationPolicyV1,
) -> FiniteMoistureSolution {
    let law = IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("authored-neutral-solute").unwrap(),
        Temperature::new(293.15),
        AmountConcentration::new(1000.0),
        0.0,
        [Temperature::new(293.15), Temperature::new(293.15)],
    )
    .unwrap();
    FiniteMoistureSolution::with_kinetic_solubility(
        KineticSolubilitySolution::new(
            Temperature::new(293.15),
            WATER_VOLUME,
            Volume::new(volume),
            vec![(
                SolubilityKineticsV1::new(law, VolumetricFlowRate::new(conductance), policy)
                    .unwrap(),
                SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))
                    .unwrap(),
            )],
        )
        .unwrap(),
        WATER_MOLAR_MASS,
    )
    .unwrap()
}

fn phase(bath: &FiniteMoistureSolution) -> SolubilityKineticStateV1 {
    *bath.kinetic_solubility().unwrap().last_steps()[0].state()
}

fn grid(length: f64, diffusion: f64, moisture: f64, cx: &Cx<'_>) -> MoistureCartesian {
    MoistureCartesian::new(
        card([diffusion; 2], 0.3),
        point(0.1),
        ["moisture_diffusivity"; 3],
        MaterialPropertySelection::SingleClaimOnly,
        MoistureCartesianGeometry {
            lengths_m: [length; 3],
            cells: [1; 3],
            frame: ContentHash([4; 32]),
            maximum_moisture: 0.3,
            maximum_cells: 256,
        },
        vec![moisture],
        cx,
    )
    .unwrap()
}

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

#[test]
fn initially_undissolved_crystals_change_native_uptake_and_do_not_take_the_pure_water_shortcut() {
    with_cx(|cx| {
        let final_volume = 5e-5;
        let conductance = 2.5e-5;
        let dt = 2.0;
        // At the manufactured endpoint volume G*dt/V=1. This independent
        // kinetic trajectory starts from zero dissolved solute, not saturation.
        let concentration = 1000.0 * (1.0 - (-1.0_f64).exp());
        let activity = (-WATER_VOLUME * concentration).exp();
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        let alpha = 4.0 * dt * 0.01;
        let expected_moisture = (0.1 + alpha * 0.2 * activity) / (1.0 + alpha);
        let water = initial.dry_mass_kg() * (expected_moisture - 0.1);
        let initial_volume = final_volume + water / WATER_DENSITY;
        let bath = kinetic_bath(
            initial_volume,
            0.0,
            0.2,
            conductance,
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let end = bath
            .advance(&initial, dt, &faces, config(), controls(initial_volume), cx)
            .unwrap();
        assert!(end.step().iterations > 0);
        assert_eq!(end.step().water_activity[0], 1.0);
        near(end.step().water_activity[1], activity, 2e-13);
        near(end.moisture().moisture()[0], expected_moisture, 2e-13);
        near(
            end.reservoir().solution().solvent_volume_m3(),
            final_volume,
            2e-16,
        );
        let actual = phase(end.reservoir());
        near(
            actual.dissolved().value(),
            concentration * final_volume,
            2e-13,
        );
        near(
            actual.solid().value(),
            0.2 - concentration * final_volume,
            2e-13,
        );
        assert_eq!(actual.total_solute().value(), 0.2);
        near(end.step().water_balance_residual_kg, 0.0, 1e-13);
        for source in end.step().transport.surface_sorption.iter().flatten() {
            near(source.query().axes()["relative_humidity"], activity, 2e-13);
        }

        let closed = kinetic_bath(
            initial_volume,
            0.0,
            0.2,
            0.0,
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let no_dissolution = closed
            .advance(&initial, dt, &faces, config(), controls(initial_volume), cx)
            .unwrap();
        assert_eq!(no_dissolution.step().water_activity, [1.0; 2]);
        assert_eq!(phase(no_dissolution.reservoir()).solid().value(), 0.2);
        assert!(end.step().water_to_solid_kg < no_dissolution.step().water_to_solid_kg);
        assert_eq!(phase(&bath).dissolved().value(), 0.0);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn sealed_reservoir_matches_fixed_volume_kinetics_and_continues_actual_phases() {
    with_cx(|cx| {
        let volume = 5e-5;
        let bath = kinetic_bath(
            volume,
            0.0,
            0.2,
            2.5e-5,
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.25, 0.1, 0.05], cx);
        let sealed = [None, None];
        let end = bath
            .advance_shared(
                &[MoistureParticipant {
                    id: 7,
                    slab: &initial,
                    faces: &sealed,
                    transport: config(),
                }],
                2.0,
                3,
                controls(volume),
                cx,
            )
            .unwrap();
        near(
            phase(end.reservoir()).dissolved().value(),
            0.05 * (1.0 - (-1.0_f64).exp()),
            2e-15,
        );
        assert_eq!(end.reservoir().solution().solvent_volume_m3(), volume);
        assert_eq!(end.step().iterations, 0);
        near(end.step().water_to_solids_kg, 0.0, 1e-13);
        assert!(end.solids()[0].moisture.moisture()[0] < 0.25);
        let next = end
            .reservoir()
            .advance(
                &end.solids()[0].moisture,
                2.0,
                &sealed,
                config(),
                controls(volume),
                cx,
            )
            .unwrap();
        let state = phase(next.reservoir());
        near(
            state.dissolved().value(),
            0.05 * (1.0 - (-2.0_f64).exp()),
            2e-15,
        );
        near(
            state.dissolved().value() + state.solid().value(),
            0.2,
            2e-16,
        );
        assert_eq!(next.moisture().time_s(), 4.0);
        near(next.step().water_balance_residual_kg, 0.0, 1e-13);
        assert_eq!(phase(&bath).solid().value(), 0.2);
    });
}

#[test]
fn closed_moisture_field_preserves_seed_policy_and_finite_precipitation_rate() {
    with_cx(|cx| {
        let volume = 5e-5;
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.1], cx);
        for policy in [
            PrecipitationPolicyV1::ExistingSolidOnly,
            PrecipitationPolicyV1::SeededSurface,
        ] {
            let bath = kinetic_bath(volume, 0.1, 0.0, 2.5e-5, policy);
            let end = bath
                .advance(&initial, 2.0, &[None, None], config(), controls(volume), cx)
                .unwrap();
            let expected = if policy == PrecipitationPolicyV1::SeededSurface {
                0.05 + 0.05 * (-1.0_f64).exp()
            } else {
                0.1
            };
            let state = phase(end.reservoir());
            near(state.dissolved().value(), expected, 2e-15);
            near(state.solid().value(), 0.1 - expected, 2e-15);
            assert_eq!(state.total_solute().value(), 0.1);
            assert_eq!(end.reservoir().solution().solvent_volume_m3(), volume);
            near(end.step().water_to_solid_kg, 0.0, 1e-14);
            near(
                end.step().water_activity[1],
                (-WATER_VOLUME * expected / volume).exp(),
                2e-14,
            );
        }
    });
}

#[test]
fn native_cartesian_and_shared_specimens_use_one_kinetic_endpoint_activity() {
    with_cx(|cx| {
        let volume = 5e-5;
        let dt = 2.0;
        let concentration = 1000.0 * (1.0 - (-1.0_f64).exp());
        let activity = (-WATER_VOLUME * concentration).exp();
        let length = [0.04_f64, 0.03];
        let diffusion = [1e-5, 2e-5];
        let old = [0.19, 0.21];
        let expected: [f64; 2] = std::array::from_fn(|i| {
            let alpha = 12.0 * dt * diffusion[i] / length[i].powi(2);
            (old[i] + alpha * 0.2 * activity) / (1.0 + alpha)
        });
        let water: [f64; 2] =
            std::array::from_fn(|i| 500.0 * length[i].powi(3) * (expected[i] - old[i]));
        assert!(water[0] > 0.0 && water[1] < 0.0);
        let initial_volume = volume + (water[0] + water[1]) / WATER_DENSITY;
        let bath = kinetic_bath(
            initial_volume,
            0.0,
            0.2,
            2.5e-5,
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let a = grid(length[0], diffusion[0], old[0], cx);
        let b = grid(length[1], diffusion[1], old[1], cx);
        let faces = std::array::from_fn(|_| Some(atmosphere::sorption(0.2)));
        let participants = [
            CartesianMoistureParticipant {
                id: 20,
                moisture: &b,
                faces: &faces,
                transport: config(),
                nonlinear: nonlinear(),
            },
            CartesianMoistureParticipant {
                id: 10,
                moisture: &a,
                faces: &faces,
                transport: config(),
                nonlinear: nonlinear(),
            },
        ];
        let end = bath
            .advance_shared_cartesian(&participants, dt, 2, controls(initial_volume), cx)
            .unwrap();
        near(
            end.reservoir().solution().solvent_volume_m3(),
            volume,
            2e-16,
        );
        near(
            phase(end.reservoir()).dissolved().value(),
            concentration * volume,
            2e-13,
        );
        for (i, body) in end.solids().iter().enumerate() {
            assert_eq!(body.id, 10 * (i as u64 + 1));
            near(body.moisture.moisture()[0], expected[i], 2e-13);
            near(
                body.transport.boundary_water_kg.iter().sum(),
                water[i],
                2e-13,
            );
            for source in body.transport.surface_sorption.iter().flatten() {
                near(source.query().axes()["relative_humidity"], activity, 2e-13);
            }
        }
        near(end.step().water_balance_residual_kg, 0.0, 1e-13);
        near(end.step().direct_water_balance_residual_kg, 0.0, 1e-13);
        let [second, first] = participants;
        let reordered = [first, second];
        let replay = bath
            .advance_shared_cartesian(&reordered, dt, 2, controls(initial_volume), cx)
            .unwrap();
        assert_eq!(
            end.reservoir().kinetic_solubility(),
            replay.reservoir().kinetic_solubility()
        );

        let solo_volume = volume + water[0] / WATER_DENSITY;
        let solo = kinetic_bath(
            solo_volume,
            0.0,
            0.2,
            2.5e-5,
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let single = solo
            .advance_cartesian(
                &a,
                dt,
                &faces,
                config(),
                nonlinear(),
                controls(solo_volume),
                cx,
            )
            .unwrap();
        near(
            single.reservoir().solution().solvent_volume_m3(),
            volume,
            2e-16,
        );
        near(single.moisture().moisture()[0], expected[0], 2e-13);
        assert_eq!(phase(&bath).dissolved().value(), 0.0);
    });
}

#[test]
fn rejected_roots_and_cancellation_leave_original_phase_and_moisture_stocks_untouched() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = kinetic_bath(
            volume,
            0.0,
            0.2,
            5e-5,
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let saved = bath.kinetic_solubility().unwrap().clone();
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let mut narrow = controls(volume);
        narrow.solvent_volume_bracket_m3 = [0.9999 * volume, 1.0001 * volume];
        assert!(matches!(
            bath.advance(&initial, 0.25, &faces, config(), narrow, cx),
            Err(SolutionMoistureError::Unbracketed { .. })
        ));
        let mut limited = controls(volume);
        limited.maximum_iterations = 1;
        assert!(matches!(
            bath.advance(&initial, 0.25, &faces, config(), limited, cx),
            Err(SolutionMoistureError::NotConverged { .. })
        ));
        assert_eq!(bath.kinetic_solubility(), Some(&saved));
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
        let accepted = bath
            .advance(&initial, 0.25, &faces, config(), controls(volume), cx)
            .unwrap();
        assert!(phase(accepted.reservoir()).dissolved().value() > 0.0);
        assert_eq!(phase(accepted.reservoir()).total_solute().value(), 0.2);
        near(accepted.step().water_balance_residual_kg, 0.0, 1e-13);
        assert_eq!(bath.kinetic_solubility(), Some(&saved));
    });

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
        let bath = kinetic_bath(
            1e-4,
            0.0,
            0.2,
            5e-5,
            PrecipitationPolicyV1::ExistingSolidOnly,
        );
        let saved = bath.kinetic_solubility().unwrap().clone();
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], &cx);
        gate.request();
        assert!(matches!(
            bath.advance(&initial, 0.25, &[None, None], config(), controls(1e-4), &cx),
            Err(SolutionMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(bath.kinetic_solubility(), Some(&saved));
        assert_eq!(initial.time_s(), 0.0);
    });
}
