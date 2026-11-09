//! G1/G3 kinetic balances through original caloric owners, including returned gas.

use super::super::tests::{interface, near, pure};
use super::*;

const AREA: f64 = 1e-10;
const DIAMETER: f64 = 1e-9;

fn channel(gas: &ClosedGasParcel, length: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        interface(gas.model().clone(), 1.0),
        CylindricalMolecularGeometry::try_new(
            DIAMETER,
            length,
            1e-7,
            "synthetic collisionless diffuse bore qualification",
        )
        .unwrap(),
    )
    .unwrap()
}

fn incident_moles(pressure_pa: f64, temperature_k: f64, molar_mass: f64) -> f64 {
    AREA * pressure_pa / (2.0 * core::f64::consts::PI * molar_mass * R * temperature_k).sqrt()
}

#[test]
fn equal_gases_and_wall_are_in_equilibrium_despite_gross_returned_molecules() {
    let gas = pure(0.028, 3.5, 0.001, 0.0, 600.0, 100_000.0);
    let law = channel(&gas, DIAMETER);
    let flux = law.sample(&gas, &gas, AREA, 600.0, || false).unwrap();
    assert_eq!(flux.transmitted.net_molar_rate_left_to_right_mol_s, 0.0);
    assert_eq!(flux.transmitted.net_mass_rate_left_to_right_kg_s, 0.0);
    assert_eq!(flux.transmitted.energy_rate_left_to_right_w, 0.0);
    assert_eq!(flux.wall.heat_into_gas_w, [0.0; 2]);
    assert_eq!(flux.wall.heat_into_wall_w, 0.0);
    assert_eq!(flux.wall.reflected_species.len(), 1);
    let reflected = flux.wall.reflected_species[0].mass_rate_kg_s;
    assert!(reflected[0] > 0.0);
    assert_eq!(reflected[0], reflected[1]);
    near(
        flux.transmitted.species[0].left_to_right.molar_rate_mol_s,
        incident_moles(100_000.0, 600.0, 0.028) * law.geometry().transmission_probability(),
    );
}

#[test]
fn wall_changes_both_returned_and_transmitted_energy_with_a_complete_balance() {
    let mass = 0.004;
    let temperature = 400.0;
    let pressures = [100_000.0, 70_000.0];
    let left = pure(mass, 2.5, 0.0, 0.0, temperature, pressures[0]);
    let right = pure(mass, 2.5, 0.0, 0.0, temperature, pressures[1]);
    let law = channel(&left, DIAMETER);
    // Independent one-diameter geometry and monatomic molecular energies.
    let tau = 18.0 / 35.0;
    let direct = 0.171_572_875_253_809_9;
    let q = pressures.map(|p| incident_moles(p, temperature, mass));
    let source_energy = 2.0 * R * temperature;
    for wall_temperature in [300.0, 900.0] {
        let flux = law
            .sample(&left, &right, AREA, wall_temperature, || false)
            .unwrap();
        let wall_energy = 2.0 * R * wall_temperature;
        let change = wall_energy - source_energy;
        for side in 0..2 {
            near(
                flux.wall.reflected_species[0].mass_rate_kg_s[side],
                q[side] * mass * (1.0 - tau),
            );
            near(
                flux.wall.heat_into_gas_w[side],
                (q[side] * (1.0 - tau) + q[1 - side] * (tau - direct)) * change,
            );
            // Directly account for every entering, returning and arriving
            // molecule rather than reusing the sampler's heat decomposition.
            let expected_gas_change = -q[side] * source_energy
                + q[side] * (1.0 - tau) * wall_energy
                + q[1 - side] * direct * source_energy
                + q[1 - side] * (tau - direct) * wall_energy;
            let transported_energy = if side == 0 {
                -flux.transmitted.energy_rate_left_to_right_w
            } else {
                flux.transmitted.energy_rate_left_to_right_w
            };
            near(
                transported_energy + flux.wall.heat_into_gas_w[side],
                expected_gas_change,
            );
        }
        near(
            flux.wall.heat_into_wall_w,
            -(q[0] + q[1]) * (1.0 - direct) * change,
        );
        near(
            flux.wall.heat_into_wall_w,
            -flux.wall.heat_into_gas_w.iter().sum::<f64>(),
        );
    }
}

#[test]
fn relative_wall_heat_survives_large_formation_references_and_nonlinear_calorics() {
    let mass = 0.028;
    let cp_r = 3.5;
    let slope = 0.001;
    let temperatures = [400.0, 1600.0];
    let pressures = [100_000.0, 200_000.0];
    let wall_temperature = 900.0;
    let mut baseline: Option<(f64, [f64; 2], f64)> = None;
    for reference in [0.0, 1e24, -1e24] {
        let left = pure(mass, cp_r, slope, reference, temperatures[0], pressures[0]);
        let right = pure(mass, cp_r, slope, reference, temperatures[1], pressures[1]);
        let law = channel(&left, DIAMETER);
        let flux = law
            .sample(&left, &right, AREA, wall_temperature, || false)
            .unwrap();
        assert_eq!(flux.transmitted.net_mass_rate_left_to_right_kg_s, 0.0);
        let tau = 18.0 / 35.0;
        let direct = 0.171_572_875_253_809_9;
        let source_heat = [0, 1].map(|side| {
            let t = temperatures[side];
            let change =
                R * (wall_temperature - t) * (cp_r - 0.5 + 0.5 * slope * (wall_temperature + t));
            incident_moles(pressures[side], t, mass) * change
        });
        for side in 0..2 {
            near(
                flux.wall.heat_into_gas_w[side],
                (1.0 - tau) * source_heat[side] + (tau - direct) * source_heat[1 - side],
            );
        }
        near(
            flux.wall.heat_into_wall_w,
            -(1.0 - direct) * (source_heat[0] + source_heat[1]),
        );
        if let Some((energy, heat, wall)) = baseline {
            near(flux.transmitted.energy_rate_left_to_right_w, energy);
            for (actual, expected) in flux.wall.heat_into_gas_w.into_iter().zip(heat) {
                near(actual, expected);
            }
            near(flux.wall.heat_into_wall_w, wall);
        } else {
            baseline = Some((
                flux.transmitted.energy_rate_left_to_right_w,
                flux.wall.heat_into_gas_w,
                flux.wall.heat_into_wall_w,
            ));
        }
    }
}

#[test]
fn zero_length_recovers_the_original_aperture_without_querying_a_wall() {
    let left = pure(0.028, 3.5, 0.001, 0.0, 400.0, 110_000.0);
    let right = pure(0.028, 3.5, 0.001, 0.0, 900.0, 80_000.0);
    let law = channel(&left, 0.0);
    let original = left
        .effusive_flux_to(&right, law.mouth(), AREA, || false)
        .unwrap();
    let sampled = law.sample(&left, &right, AREA, f64::NAN, || false).unwrap();
    assert_eq!(sampled.transmitted, original);
    assert!(sampled.wall.reflected_species.is_empty());
    assert_eq!(sampled.wall.heat_into_gas_w, [0.0; 2]);
    assert_eq!(sampled.wall.heat_into_wall_w, 0.0);
}

#[test]
fn zero_area_leaves_inactive_wall_and_foreign_gas_unqueried_but_observes_cancellation() {
    let original = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let foreign = pure(0.004, 2.5, 0.0, 0.0, 1000.0, 100_000.0);
    let law = channel(&original, DIAMETER);
    let dormant = law
        .sample(&foreign, &foreign, 0.0, f64::NAN, || false)
        .unwrap();
    assert!(dormant.transmitted.species.is_empty());
    assert!(dormant.wall.reflected_species.is_empty());
    assert_eq!(dormant.wall.heat_into_gas_w, [0.0; 2]);
    assert_eq!(dormant.wall.heat_into_wall_w, 0.0);
    assert!(
        law.sample(&foreign, &foreign, AREA, 600.0, || false)
            .is_err()
    );
    for invalid_temperature in [f64::NAN, 199.0, 2001.0] {
        assert!(
            law.sample(&original, &original, AREA, invalid_temperature, || false)
                .is_err()
        );
    }
    assert!(matches!(
        law.sample(&original, &original, 0.0, f64::NAN, || true),
        Err(GasEnergyError::Cancelled)
    ));
}

#[test]
fn reversing_reservoirs_swaps_gas_heat_and_reflections_and_keeps_wall_heat() {
    let left = pure(0.028, 3.5, 0.001, 0.0, 400.0, 110_000.0);
    let right = pure(0.028, 3.5, 0.001, 0.0, 900.0, 80_000.0);
    let law = channel(&left, 2.0 * DIAMETER);
    let forward = law.sample(&left, &right, AREA, 650.0, || false).unwrap();
    let reverse = law.sample(&right, &left, AREA, 650.0, || false).unwrap();
    near(
        forward.transmitted.energy_rate_left_to_right_w,
        -reverse.transmitted.energy_rate_left_to_right_w,
    );
    near(
        forward.transmitted.net_mass_rate_left_to_right_kg_s,
        -reverse.transmitted.net_mass_rate_left_to_right_kg_s,
    );
    for side in 0..2 {
        near(
            forward.wall.heat_into_gas_w[side],
            reverse.wall.heat_into_gas_w[1 - side],
        );
        near(
            forward.wall.reflected_species[0].mass_rate_kg_s[side],
            reverse.wall.reflected_species[0].mass_rate_kg_s[1 - side],
        );
    }
    near(forward.wall.heat_into_wall_w, reverse.wall.heat_into_wall_w);
}

#[test]
fn channel_requires_matching_unit_transmission_mouths() {
    let gas = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let geometry = channel(&gas, DIAMETER).geometry().clone();
    assert!(
        DiffuseMolecularChannel::try_new(interface(gas.model().clone(), 0.5), geometry).is_err()
    );
    for (diameter, mean_free_path) in [(2.0 * DIAMETER, 1e-7), (DIAMETER, 1e-6)] {
        let geometry = CylindricalMolecularGeometry::try_new(
            diameter,
            DIAMETER,
            mean_free_path,
            "synthetic mismatched mouth",
        )
        .unwrap();
        assert!(
            DiffuseMolecularChannel::try_new(interface(gas.model().clone(), 1.0), geometry)
                .is_err()
        );
    }
}

#[test]
fn resolved_geometry_drives_original_species_and_accommodation_energy() {
    let mass = 0.004;
    let temperature = 400.0;
    let pressures = [100_000.0, 70_000.0];
    let left = pure(mass, 2.5, 0.0, 0.0, temperature, pressures[0]);
    let right = pure(mass, 2.5, 0.0, 0.0, temperature, pressures[1]);
    let geometry = channel(&left, DIAMETER)
        .geometry()
        .clone()
        .resolve_transmission(
            ClausingControls {
                initial_rings: 32,
                maximum_rings: 1024,
                maximum_iterations: 256,
                absolute_tolerance: 3e-5,
                relative_tolerance: 0.0,
            },
            || false,
        )
        .unwrap();
    assert!((geometry.transmission_probability() - 18.0 / 35.0).abs() > 4e-5);
    let law =
        DiffuseMolecularChannel::try_new(interface(left.model().clone(), 1.0), geometry.clone())
            .unwrap();
    assert_eq!(law.geometry().resolution(), geometry.resolution());
    let wall_temperature = 900.0;
    let flux = law
        .sample(&left, &right, AREA, wall_temperature, || false)
        .unwrap();
    let q = pressures.map(|p| incident_moles(p, temperature, mass));
    let change = 2.0 * R * (wall_temperature - temperature);
    near(
        flux.transmitted.species[0].left_to_right.molar_rate_mol_s,
        q[0] * geometry.transmission_probability(),
    );
    near(
        flux.transmitted.species[0].right_to_left.molar_rate_mol_s,
        q[1] * geometry.transmission_probability(),
    );
    for side in 0..2 {
        near(
            flux.wall.reflected_species[0].mass_rate_kg_s[side],
            q[side] * mass * geometry.reflected_probability(),
        );
        near(
            flux.wall.heat_into_gas_w[side],
            (q[side] * geometry.reflected_probability()
                + q[1 - side] * geometry.wall_transmitted_probability())
                * change,
        );
    }
    near(
        flux.wall.heat_into_wall_w,
        -(q[0] + q[1]) * geometry.wall_hit_probability() * change,
    );
    near(
        flux.wall.heat_into_wall_w,
        -flux.wall.heat_into_gas_w.iter().sum::<f64>(),
    );
}

fn accommodation_controls() -> ClausingControls {
    ClausingControls {
        initial_rings: 8,
        maximum_rings: 512,
        maximum_iterations: 256,
        absolute_tolerance: 1e-6,
        relative_tolerance: 1e-4,
    }
}

#[test]
fn elastic_diffuse_walls_keep_all_gross_paths_but_query_no_wall_calorics() {
    let left = pure(0.028, 3.5, 0.001, 0.0, 400.0, 110_000.0);
    let right = pure(0.028, 3.5, 0.001, 0.0, 900.0, 80_000.0);
    let original = channel(&left, 2.0 * DIAMETER);
    let elastic = original
        .clone()
        .with_thermal_accommodation(
            0.0,
            accommodation_controls(),
            "synthetic diffuse elastic wall",
            || false,
        )
        .unwrap();
    let full = original
        .sample(&left, &right, AREA, 650.0, || false)
        .unwrap();
    let flux = elastic
        .sample(&left, &right, AREA, f64::NAN, || false)
        .unwrap();
    assert_eq!(flux.transmitted, full.transmitted);
    assert_eq!(flux.wall.reflected_species, full.wall.reflected_species);
    assert!(flux.wall.reflected_species[0].mass_rate_kg_s[0] > 0.0);
    assert_eq!(flux.wall.heat_into_gas_w, [0.0; 2]);
    assert_eq!(flux.wall.heat_into_wall_w, 0.0);
    elastic
        .admit_wall_for_parcel(&left, f64::NAN, || false)
        .unwrap();
    let vacuum = elastic
        .sample_to_vacuum(&left, AREA, f64::NAN, || false)
        .unwrap();
    assert!(vacuum.wall.reflected_species[0].mass_rate_kg_s > 0.0);
    assert_eq!(vacuum.wall.heat_into_gas_w, 0.0);
    assert_eq!(vacuum.wall.heat_into_escaped_w, 0.0);
    assert_eq!(vacuum.wall.heat_into_wall_w, 0.0);
    assert_eq!(
        vacuum.escaped_energy_rate_w,
        vacuum.transmitted.energy_rate_w
    );
}

#[test]
fn unit_collision_accommodation_preserves_the_original_full_wall_results() {
    let left = pure(0.028, 3.5, 0.001, 0.0, 400.0, 110_000.0);
    let right = pure(0.028, 3.5, 0.001, 0.0, 900.0, 80_000.0);
    let original = channel(&left, 2.0 * DIAMETER);
    let full = original
        .clone()
        .with_thermal_accommodation(
            1.0,
            accommodation_controls(),
            "synthetic fully accommodating wall",
            || false,
        )
        .unwrap();
    assert_eq!(original.geometry(), full.geometry());
    assert_eq!(full.thermal_accommodation().unwrap().coefficient(), 1.0);
    assert_eq!(
        full.thermal_accommodation_provenance(),
        Some("synthetic fully accommodating wall")
    );
    assert_eq!(
        original
            .sample(&left, &right, AREA, 650.0, || false)
            .unwrap(),
        full.sample(&left, &right, AREA, 650.0, || false).unwrap(),
    );
    assert_eq!(
        original
            .sample_to_vacuum(&left, AREA, 650.0, || false)
            .unwrap(),
        full.sample_to_vacuum(&left, AREA, 650.0, || false).unwrap(),
    );
}

#[test]
fn repeated_collision_heat_consumes_original_nonlinear_calorics_without_changing_mass() {
    let temperatures = [400.0, 600.0];
    let pressures = [110_000.0, 80_000.0];
    let tw = 900.0;
    let alpha = 0.25;
    let mut baseline = None;
    for datum in [0.0, 1e24] {
        let left = pure(0.028, 3.5, 0.001, datum, temperatures[0], pressures[0]);
        let right = pure(0.028, 3.5, 0.001, datum, temperatures[1], pressures[1]);
        let partial = channel(&left, 3.0 * DIAMETER)
            .with_thermal_accommodation(
                alpha,
                accommodation_controls(),
                "synthetic common per-hit thermalization",
                || false,
            )
            .unwrap();
        let full =
            DiffuseMolecularChannel::try_new(partial.mouth().clone(), partial.geometry().clone())
                .unwrap();
        let flux = partial.sample(&left, &right, AREA, tw, || false).unwrap();
        let whole = full.sample(&left, &right, AREA, tw, || false).unwrap();
        assert_eq!(flux.transmitted, whole.transmitted);
        assert_eq!(flux.wall.reflected_species, whole.wall.reflected_species);
        let model = partial.thermal_accommodation().unwrap();
        let hr = model.reflected_heat_probability();
        let ht = model.transmitted_heat_probability();
        assert!(hr > alpha * partial.geometry().reflected_probability());
        assert!(ht > alpha * partial.geometry().wall_transmitted_probability());
        // Independent NASA fixture: e_cross/R = 3T + 0.0005 T² + datum.
        let source_heat = [0, 1].map(|side| {
            let t = temperatures[side];
            incident_moles(pressures[side], t, 0.028) * R * (tw - t) * (3.0 + 0.0005 * (tw + t))
        });
        for side in 0..2 {
            near(
                flux.wall.heat_into_gas_w[side],
                hr * source_heat[side] + ht * source_heat[1 - side],
            );
            assert!(flux.wall.heat_into_gas_w[side] > alpha * whole.wall.heat_into_gas_w[side]);
            assert!(flux.wall.heat_into_gas_w[side] < whole.wall.heat_into_gas_w[side]);
        }
        near(
            flux.wall.heat_into_wall_w,
            -(hr + ht) * source_heat.iter().sum::<f64>(),
        );
        if let Some(old) = baseline {
            assert_eq!(flux.wall, old);
        }
        baseline = Some(flux.wall);
    }
}

#[test]
fn partial_accommodation_requires_explicit_physical_input_and_preserves_refusals() {
    let gas = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let original = channel(&gas, DIAMETER);
    for alpha in [-0.1, 1.1, f64::NAN] {
        assert!(
            original
                .clone()
                .with_thermal_accommodation(
                    alpha,
                    accommodation_controls(),
                    "synthetic surface",
                    || false,
                )
                .is_err()
        );
    }
    assert!(
        original
            .clone()
            .with_thermal_accommodation(0.25, accommodation_controls(), " ", || false,)
            .is_err()
    );
    assert!(matches!(
        original.clone().with_thermal_accommodation(
            0.25,
            accommodation_controls(),
            "synthetic surface",
            || true,
        ),
        Err(GasEnergyError::Cancelled)
    ));
    let partial = original
        .with_thermal_accommodation(0.25, accommodation_controls(), "synthetic surface", || {
            false
        })
        .unwrap();
    assert!(
        partial
            .sample(&gas, &gas, AREA, f64::NAN, || false)
            .is_err()
    );
    assert!(
        partial
            .admit_wall_for_parcel(&gas, 199.0, || false)
            .is_err()
    );
}
