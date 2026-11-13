//! G1 finite solution balances and G3/G4 native material coupling/refusals.
//! Solution and specimen coefficients are authored references, not calibration.
use super::*;
use fs_couple::moisture::solution::*;
use fs_material::osmosis::{
    AqueousSolution, IdealDiluteSolution, OsmosisError,
    nonideal::{NonidealSolutionError, VirialSolutionModel},
    solubility::FiniteSolubilitySolution,
};
use fs_qty::{Amount, AmountConcentration, Temperature, Volume};
use fs_thermochem::{SpeciesId, solution::solubility::IdealSolubilityV1};
use std::sync::Arc;

#[path = "solution/electrolyte.rs"]
mod electrolyte;
#[path = "solution/kinetics.rs"]
mod kinetics;
#[path = "solution/paired.rs"]
mod paired;
#[path = "solution/shared.rs"]
mod shared;
#[path = "solution/slab_sdirk.rs"]
mod slab_sdirk;

pub(super) const WATER_VOLUME: f64 = 18e-6;
const WATER_MOLAR_MASS: f64 = 0.018;
const WATER_DENSITY: f64 = WATER_MOLAR_MASS / WATER_VOLUME;

pub(super) fn reservoir(volume: f64, moles: Vec<f64>) -> FiniteMoistureSolution {
    FiniteMoistureSolution::new(
        IdealDiluteSolution::new(293.15, WATER_VOLUME, volume, moles).unwrap(),
        WATER_MOLAR_MASS,
    )
    .unwrap()
}
pub(super) fn controls(volume: f64) -> SolutionMoistureSolve {
    SolutionMoistureSolve {
        solvent_volume_bracket_m3: [0.5 * volume, 1.5 * volume],
        maximum_iterations: 100,
        mass_tolerance_kg: 1e-13,
    }
}

fn soluble_reservoir(volume: f64, total_moles: f64, saturation: f64) -> FiniteMoistureSolution {
    let law = IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new("authored-neutral-solute").unwrap(),
        Temperature::new(293.15),
        AmountConcentration::new(saturation),
        0.0,
        [Temperature::new(293.15), Temperature::new(293.15)],
    )
    .unwrap();
    FiniteMoistureSolution::with_solubility(
        FiniteSolubilitySolution::new(
            Temperature::new(293.15),
            WATER_VOLUME,
            Volume::new(volume),
            vec![(law, Amount::new(total_moles))],
        )
        .unwrap(),
        WATER_MOLAR_MASS,
    )
    .unwrap()
}

#[test]
fn solubility_buffers_actual_moisture_activity_and_preserves_finite_crystals() {
    with_cx(|cx| {
        let final_volume = 5e-5;
        let saturation = 100.0;
        let activity = (-WATER_VOLUME * saturation).exp();
        let slope = 0.2;
        let diffusion = 0.01;
        let dt = 2.0;
        let alpha = 4.0 * dt * diffusion;
        for initial_moisture in [0.198, 0.201] {
            let initial = slab(card([diffusion; 2], 0.3), 1.0, vec![initial_moisture], cx);
            let final_moisture = (initial_moisture + alpha * slope * activity) / (1.0 + alpha);
            let water = initial.dry_mass_kg() * (final_moisture - initial_moisture);
            let initial_volume = final_volume + water / WATER_DENSITY;
            let total = 2.0 * saturation * initial_volume.max(final_volume);
            let bath = soluble_reservoir(initial_volume, total, saturation);
            let before = bath.solubility().unwrap().clone();
            let faces = [
                Some(atmosphere::sorption(slope)),
                Some(atmosphere::sorption(slope)),
            ];
            let endpoint = bath
                .advance(&initial, dt, &faces, config(), controls(initial_volume), cx)
                .unwrap();
            near(endpoint.moisture().moisture()[0], final_moisture, 2e-13);
            near(
                endpoint.reservoir().solution().solvent_volume_m3(),
                final_volume,
                2e-16,
            );
            for actual in endpoint.step().water_activity {
                near(actual, activity, 2e-14);
            }
            let phase = endpoint
                .reservoir()
                .solubility()
                .unwrap()
                .phase_partitions()[0];
            near(phase.dissolved.value(), saturation * final_volume, 2e-13);
            near(
                phase.solid.value(),
                total - saturation * final_volume,
                2e-13,
            );
            near(phase.dissolved.value() + phase.solid.value(), total, 1e-16);
            near(
                phase.solid.value() - before.phase_partitions()[0].solid.value(),
                saturation * water / WATER_DENSITY,
                2e-13,
            );
            near(endpoint.step().water_balance_residual_kg, 0.0, 1e-13);
            for source in endpoint.step().transport.surface_sorption.iter().flatten() {
                near(source.query().axes()["relative_humidity"], activity, 2e-14);
            }
            let continued = endpoint
                .reservoir()
                .advance(
                    endpoint.moisture(),
                    dt,
                    &faces,
                    config(),
                    controls(final_volume),
                    cx,
                )
                .unwrap();
            let continued_phase = continued
                .reservoir()
                .solubility()
                .unwrap()
                .phase_partitions()[0];
            near(
                continued_phase.dissolved.value() + continued_phase.solid.value(),
                total,
                1e-16,
            );
            assert_eq!(bath.solubility(), Some(&before));
            assert_eq!(initial.moisture(), &[initial_moisture]);
        }
    });
}

#[test]
fn native_water_exchange_crosses_saturation_and_exhausts_or_creates_solid() {
    with_cx(|cx| {
        let final_volume = 5e-5;
        let saturation = 1000.0;
        let slope = 0.2;
        let diffusion = 0.01;
        let dt = 2.0;
        let alpha = 4.0 * dt * diffusion;
        for (initial_moisture, final_concentration) in [(0.1, saturation), (0.26, 980.0)] {
            let activity = (-WATER_VOLUME * final_concentration).exp();
            let initial = slab(card([diffusion; 2], 0.3), 1.0, vec![initial_moisture], cx);
            let final_moisture = (initial_moisture + alpha * slope * activity) / (1.0 + alpha);
            let water = initial.dry_mass_kg() * (final_moisture - initial_moisture);
            let initial_volume = final_volume + water / WATER_DENSITY;
            let total = if water > 0.0 {
                saturation * 0.5 * (initial_volume + final_volume)
            } else {
                final_concentration * final_volume
            };
            let bath = soluble_reservoir(initial_volume, total, saturation);
            let faces = [
                Some(atmosphere::sorption(slope)),
                Some(atmosphere::sorption(slope)),
            ];
            let endpoint = bath
                .advance(&initial, dt, &faces, config(), controls(initial_volume), cx)
                .unwrap();
            near(
                endpoint.reservoir().solution().solvent_volume_m3(),
                final_volume,
                2e-16,
            );
            near(endpoint.moisture().moisture()[0], final_moisture, 2e-13);
            near(endpoint.step().water_activity[1], activity, 2e-13);
            let start_solid = bath.solubility().unwrap().phase_partitions()[0]
                .solid
                .value();
            let end = endpoint
                .reservoir()
                .solubility()
                .unwrap()
                .phase_partitions()[0];
            if water > 0.0 {
                assert_eq!(start_solid, 0.0);
                assert!(end.solid.value() > 0.0);
            } else {
                assert!(start_solid > 0.0);
                assert_eq!(end.solid.value(), 0.0);
                assert_eq!(end.dissolved.value(), total);
            }
            near(end.dissolved.value() + end.solid.value(), total, 1e-15);
            near(endpoint.step().water_balance_residual_kg, 0.0, 1e-13);
        }
    });
}

#[test]
fn precipitated_solute_changes_native_wet_plate_and_failure_retains_both_stocks() {
    with_cx(|cx| {
        let volume = 1e-4;
        let total = 0.2;
        let bath = soluble_reservoir(volume, total, 1000.0);
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1; 4], cx);
        let plate = plate(&initial, PlateMaterialModel::Isotropic, cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let solve = controls(volume);
        let endpoint = bath
            .advance(&initial, 0.25, &faces, config(), solve, cx)
            .unwrap();
        let all_dissolved = reservoir(volume, vec![total])
            .advance(&initial, 0.25, &faces, config(), solve, cx)
            .unwrap();
        assert!(endpoint.step().water_to_solid_kg > all_dissolved.step().water_to_solid_kg);
        let wet = plate.resolve(endpoint.moisture(), cx).unwrap();
        let overdissolved = plate.resolve(all_dissolved.moisture(), cx).unwrap();
        assert!(wet.wet_mass_kg() > overdissolved.wet_mass_kg());
        near(wet.wet_mass_kg(), endpoint.moisture().wet_mass_kg(), 1e-13);
        let before = bath.solubility().unwrap().clone();
        let too_small = SolutionMoistureSolve {
            solvent_volume_bracket_m3: [0.999999999 * volume, 1.000000001 * volume],
            ..solve
        };
        assert!(matches!(
            bath.advance(&initial, 0.25, &faces, config(), too_small, cx),
            Err(SolutionMoistureError::Unbracketed { .. }),
        ));
        assert_eq!(bath.solubility(), Some(&before));
        assert_eq!(initial.moisture(), &[0.1; 4]);
    });
}

#[test]
fn finite_solution_endpoint_matches_manufactured_uptake_and_desorption() {
    with_cx(|cx| {
        let final_volume = 5e-5;
        let slope = 0.2;
        let diffusion = 0.01;
        let dt = 2.0;
        let alpha = 4.0 * dt * diffusion; // Two half-cell faces; unit slab length.
        for (final_activity, initial_moisture) in [(0.998_f64, 0.198), (0.995_f64, 0.201)] {
            let initial = slab(card([diffusion; 2], 0.3), 1.0, vec![initial_moisture], cx);
            let final_moisture =
                (initial_moisture + alpha * slope * final_activity) / (1.0 + alpha);
            let water = initial.dry_mass_kg() * (final_moisture - initial_moisture);
            let initial_volume = final_volume + water / WATER_DENSITY;
            let solute = -final_activity.ln() * final_volume / WATER_VOLUME;
            let bath = reservoir(initial_volume, vec![0.6 * solute, 0.4 * solute]);
            let faces = [
                Some(atmosphere::sorption(slope)),
                Some(atmosphere::sorption(slope)),
            ];
            let solve = controls(initial_volume);
            let endpoint = bath
                .advance(&initial, dt, &faces, config(), solve, cx)
                .unwrap();
            near(
                endpoint.reservoir().solution().solvent_volume_m3(),
                final_volume,
                2e-16,
            );
            near(
                endpoint.reservoir().solution().water_activity().unwrap(),
                final_activity,
                1e-12,
            );
            near(endpoint.moisture().moisture()[0], final_moisture, 2e-13);
            near(endpoint.step().water_to_solid_kg, water, 2e-13);
            near(
                endpoint.step().water_balance_residual_kg,
                0.0,
                solve.mass_tolerance_kg,
            );
            near(
                endpoint.reservoir().solvent_mass_kg() + endpoint.moisture().water_mass_kg(),
                bath.solvent_mass_kg() + initial.water_mass_kg(),
                2e-13,
            );
            assert_eq!(
                endpoint.reservoir().solution().solute_moles(),
                bath.solution().solute_moles()
            );
            assert_eq!(
                endpoint.step().water_to_solid_kg.is_sign_positive(),
                water.is_sign_positive()
            );
            for source in endpoint.step().transport.surface_sorption.iter().flatten() {
                near(
                    source.query().axes()["relative_humidity"],
                    final_activity,
                    1e-12,
                );
            }
            let repeat = bath
                .advance(&initial, dt, &faces, config(), solve, cx)
                .unwrap();
            assert_eq!(repeat.moisture().moisture(), endpoint.moisture().moisture());
            assert_eq!(repeat.step().water_activity, endpoint.step().water_activity);
            assert_eq!(initial.moisture(), &[initial_moisture]);
            assert_eq!(initial.time_s(), 0.0);
            assert_eq!(bath.solution().solvent_volume_m3(), initial_volume);
        }
    });
}

#[test]
fn nonideal_activity_reaches_manufactured_native_moisture_endpoint_and_source_queries() {
    with_cx(|cx| {
        let final_activity = 0.95_f64;
        let final_volume = 5e-5;
        let initial_moisture = 0.18;
        let slope = 0.2;
        let diffusion = 0.01;
        let dt = 2.0;
        let interaction = 0.002;
        let initial = slab(card([diffusion; 2], 0.3), 1.0, vec![initial_moisture], cx);
        let alpha = 4.0 * dt * diffusion;
        let final_moisture = (initial_moisture + alpha * slope * final_activity) / (1.0 + alpha);
        let water = initial.dry_mass_kg() * (final_moisture - initial_moisture);
        let initial_volume = final_volume + water / WATER_DENSITY;
        // Independent inversion of -ln(a_w)/v_water = c + B*c²/2.
        let osmotic_concentration = -final_activity.ln() / WATER_VOLUME;
        let concentration = 2.0 * osmotic_concentration
            / (1.0 + (1.0 + 2.0 * interaction * osmotic_concentration).sqrt());
        let amount = concentration * final_volume;
        let model = Arc::new(
            VirialSolutionModel::new(293.15, WATER_VOLUME, 1.0, vec![3000.0], vec![interaction])
                .unwrap(),
        );
        let bath = FiniteMoistureSolution::new(
            AqueousSolution::nonideal(model.clone(), initial_volume, vec![amount]).unwrap(),
            WATER_MOLAR_MASS,
        )
        .unwrap();
        let original = bath.solution().clone();
        let faces = [
            Some(atmosphere::sorption(slope)),
            Some(atmosphere::sorption(slope)),
        ];
        let solve = controls(initial_volume);
        let endpoint = bath
            .advance(&initial, dt, &faces, config(), solve, cx)
            .unwrap();
        near(
            endpoint.reservoir().solution().solvent_volume_m3(),
            final_volume,
            2e-16,
        );
        near(
            endpoint.reservoir().solution().water_activity().unwrap(),
            final_activity,
            1e-12,
        );
        near(endpoint.moisture().moisture()[0], final_moisture, 2e-13);
        near(endpoint.step().water_to_solid_kg, water, 2e-13);
        near(
            endpoint.step().water_balance_residual_kg,
            0.0,
            solve.mass_tolerance_kg,
        );
        near(
            endpoint.reservoir().solvent_mass_kg() + endpoint.moisture().water_mass_kg(),
            bath.solvent_mass_kg() + initial.water_mass_kg(),
            2e-13,
        );
        assert_eq!(endpoint.reservoir().solution().solute_moles(), &[amount]);
        assert!(std::ptr::eq(
            endpoint.reservoir().solution().nonideal_model().unwrap(),
            model.as_ref(),
        ));
        assert_eq!(
            endpoint
                .step()
                .transport
                .surface_sorption
                .iter()
                .flatten()
                .count(),
            2
        );
        for source in endpoint.step().transport.surface_sorption.iter().flatten() {
            near(
                source.query().axes()["relative_humidity"],
                final_activity,
                1e-12,
            );
        }
        let ideal_activity =
            IdealDiluteSolution::new(293.15, WATER_VOLUME, final_volume, vec![amount])
                .unwrap()
                .water_activity()
                .unwrap();
        assert!(ideal_activity - final_activity > 0.02);
        assert_eq!(bath.solution(), &original);
        assert_eq!(initial.moisture(), &[initial_moisture]);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn nonideal_concentration_domain_refusal_preserves_bath_and_native_slab() {
    with_cx(|cx| {
        let volume = 1e-4;
        let model = Arc::new(
            VirialSolutionModel::new(293.15, WATER_VOLUME, 1.0, vec![1100.0], vec![0.002]).unwrap(),
        );
        let bath = FiniteMoistureSolution::new(
            AqueousSolution::nonideal(model, volume, vec![0.1]).unwrap(),
            WATER_MOLAR_MASS,
        )
        .unwrap();
        let original = bath.solution().clone();
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        // The broad numerical bracket proposes 2000 mol/m³ at its lower end;
        // the original virial source admits at most 1100. No extrapolation.
        assert!(matches!(
            bath.advance(&initial, 0.25, &faces, config(), controls(volume), cx),
            Err(SolutionMoistureError::Solution(OsmosisError::Nonideal(
                NonidealSolutionError::OutsideDomain { .. }
            )))
        ));
        assert_eq!(bath.solution(), &original);
        assert_eq!(initial.moisture(), &[0.1]);
        assert_eq!(initial.time_s(), 0.0);
        let mut supported = controls(volume);
        supported.solvent_volume_bracket_m3 = [0.95 * volume, 1.05 * volume];
        let endpoint = bath
            .advance(&initial, 0.25, &faces, config(), supported, cx)
            .unwrap();
        assert!(endpoint.step().water_to_solid_kg > 0.0);
        assert_eq!(endpoint.reservoir().solution().solute_moles(), &[0.1]);
        assert_eq!(bath.solution(), &original);
        assert_eq!(initial.moisture(), &[0.1]);
    });
}

#[test]
fn solute_reduces_uptake_and_changes_native_wet_plate_without_losing_solute() {
    with_cx(|cx| {
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1; 4], cx);
        let plate = plate(&initial, PlateMaterialModel::Isotropic, cx);
        let volume = 1e-4;
        let pure = reservoir(volume, vec![0.0]);
        let saline = reservoir(volume, vec![300.0 * volume]);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let pure_end = pure
            .advance(&initial, 0.25, &faces, config(), controls(volume), cx)
            .unwrap();
        assert_eq!(pure_end.step().iterations, 0);
        assert_eq!(pure_end.step().water_activity, [1.0, 1.0]);
        let salt_end = saline
            .advance(&initial, 0.25, &faces, config(), controls(volume), cx)
            .unwrap();
        assert!(salt_end.step().water_to_solid_kg < pure_end.step().water_to_solid_kg);
        assert!(salt_end.step().water_to_solid_kg > 0.0);
        assert!(salt_end.step().osmotic_pressure_pa[1] > salt_end.step().osmotic_pressure_pa[0]);
        assert!(salt_end.step().water_activity[1] < salt_end.step().water_activity[0]);
        assert_eq!(
            salt_end.reservoir().solution().solute_moles(),
            saline.solution().solute_moles()
        );
        let pure_plate = plate.resolve(pure_end.moisture(), cx).unwrap();
        let salt_plate = plate.resolve(salt_end.moisture(), cx).unwrap();
        assert!(salt_plate.wet_mass_kg() < pure_plate.wet_mass_kg());
        near(
            salt_plate.wet_mass_kg(),
            salt_end.moisture().wet_mass_kg(),
            1e-13,
        );
        let first_e = |state: &MoisturePlateState| {
            state.elastic_states()[0]
                .property("young_modulus")
                .unwrap()
                .value_si()
        };
        assert!(first_e(&salt_plate) > first_e(&pure_plate));
        let (next_bath, next_slab, _) = salt_end.into_parts();
        let next = next_bath
            .advance(
                &next_slab,
                0.25,
                &faces,
                config(),
                controls(next_bath.solution().solvent_volume_m3()),
                cx,
            )
            .unwrap();
        assert_eq!(
            next.reservoir().solution().solute_moles(),
            saline.solution().solute_moles()
        );
        assert_eq!(next.moisture().time_s(), 0.5);
        near(
            next.reservoir().solvent_mass_kg() + next.moisture().water_mass_kg(),
            saline.solvent_mass_kg() + initial.water_mass_kg(),
            3e-13,
        );
    });
}

#[test]
fn equilibrium_and_sealed_faces_do_not_spend_solution_inventory() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = reservoir(volume, vec![0.01]);
        let activity = bath.solution().water_activity().unwrap();
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let uniform = slab(card([0.01; 2], 0.3), 1.0, vec![0.2 * activity], cx);
        let endpoint = bath
            .advance(&uniform, 1.0, &faces, config(), controls(volume), cx)
            .unwrap();
        assert_eq!(endpoint.step().iterations, 0);
        assert_eq!(endpoint.reservoir().solution().solvent_volume_m3(), volume);
        near(endpoint.step().water_to_solid_kg, 0.0, 1e-14);

        let nonuniform = slab(card([0.01; 2], 0.3), 1.0, vec![0.25, 0.1, 0.05], cx);
        let sealed = bath
            .advance(
                &nonuniform,
                1.0,
                &[None, None],
                config(),
                controls(volume),
                cx,
            )
            .unwrap();
        assert_eq!(sealed.reservoir().solution().solvent_volume_m3(), volume);
        assert_eq!(
            sealed.reservoir().solution().solute_moles(),
            bath.solution().solute_moles()
        );
        assert!(sealed.moisture().moisture()[0] < nonuniform.moisture()[0]);
        assert!(sealed.moisture().moisture()[2] > nonuniform.moisture()[2]);
        near(sealed.step().water_to_solid_kg, 0.0, 1e-14);
    });
}

#[test]
fn unsupported_volume_temperature_and_iteration_budget_preserve_originals() {
    with_cx(|cx| {
        let volume = 1e-4;
        let bath = reservoir(volume, vec![0.02]);
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
        let mut too_few = controls(volume);
        too_few.maximum_iterations = 1;
        assert!(matches!(
            bath.advance(&initial, 0.25, &faces, config(), too_few, cx),
            Err(SolutionMoistureError::NotConverged { .. })
        ));
        let hot = FiniteMoistureSolution::new(
            IdealDiluteSolution::new(303.15, WATER_VOLUME, volume, vec![0.02]).unwrap(),
            WATER_MOLAR_MASS,
        )
        .unwrap();
        assert!(matches!(
            hot.advance(&initial, 0.25, &faces, config(), controls(volume), cx),
            Err(SolutionMoistureError::Invalid(_))
        ));
        let good = bath
            .advance(&initial, 0.25, &faces, config(), controls(volume), cx)
            .unwrap();
        assert!(good.step().water_to_solid_kg > 0.0);
        assert_eq!(bath.solution().solvent_volume_m3(), volume);
        assert_eq!(bath.solution().solute_moles(), &[0.02]);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[0.1]);
    });
}

#[test]
fn pure_solvent_refuses_exhaustion_and_accepts_a_smaller_physical_interval() {
    with_cx(|cx| {
        let volume = 1e-7;
        let bath = reservoir(volume, Vec::new());
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        // At dt=.25, native one-cell uptake is .0025 kg but the complete
        // solvent inventory is only .0001 kg. No negative volume or clipping.
        assert!(matches!(
            bath.advance(&initial, 0.25, &faces, config(), controls(volume), cx),
            Err(SolutionMoistureError::Unbracketed { .. })
        ));
        let dt = 0.001;
        let endpoint = bath
            .advance(&initial, dt, &faces, config(), controls(volume), cx)
            .unwrap();
        let alpha = 4.0 * dt * 0.01 / 0.1_f64.powi(2);
        let expected_m = (0.1 + alpha * 0.2) / (1.0 + alpha);
        let expected_water = initial.dry_mass_kg() * (expected_m - 0.1);
        near(endpoint.moisture().moisture()[0], expected_m, 1e-13);
        near(
            endpoint.reservoir().solution().solvent_volume_m3(),
            volume - expected_water / WATER_DENSITY,
            1e-18,
        );
        assert_eq!(endpoint.step().iterations, 0);
        assert_eq!(endpoint.step().water_activity, [1.0, 1.0]);
        assert!(endpoint.reservoir().solution().solute_moles().is_empty());
        assert_eq!(bath.solution().solvent_volume_m3(), volume);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[0.1]);
    });
}

#[test]
fn cancellation_cannot_publish_either_water_inventory() {
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
        let bath = reservoir(volume, vec![0.02]);
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], &cx);
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        gate.request();
        assert!(matches!(
            bath.advance(&initial, 0.25, &faces, config(), controls(volume), &cx),
            Err(SolutionMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(bath.solution().solvent_volume_m3(), volume);
        assert_eq!(bath.solution().solute_moles(), &[0.02]);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[0.1]);
    });
}
