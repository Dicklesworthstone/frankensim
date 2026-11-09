//! Mixture-specific energy exchange through the original caloric/source owners.
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    tests::{point, species_in_domain, transport},
};
use fs_qty::Amount;
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

fn mixture(fractions: [f64; 2], state: [f64; 2], reference: f64, a_upper: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_in_domain("A", 0.028, 3.5, 0.001, reference, [200.0, a_upper]),
                species_in_domain("B", 0.004, 2.5, 0.0, -reference, [200.0, 2000.0]),
            ],
            Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, point(state[0], state[1]), 0.001).unwrap()
}

fn qualification(name: &str, alpha: f64) -> SpeciesThermalAccommodation {
    SpeciesThermalAccommodation::try_new(
        SpeciesId::new(name).unwrap(),
        alpha,
        format!("synthetic {name}/surface per-collision energy qualification"),
    )
    .unwrap()
}

#[test]
fn distinct_species_weights_preserve_gross_transport_and_close_mixture_heat() {
    let temperatures = [400.0, 600.0];
    let pressures = [110_000.0, 80_000.0];
    let fractions = [[0.75, 0.25], [0.25, 0.75]];
    let tw = 900.0;
    let mut baseline = None;
    for reference in [0.0, 1e24] {
        let left = mixture(
            fractions[0],
            [temperatures[0], pressures[0]],
            reference,
            2000.0,
        );
        let right = mixture(
            fractions[1],
            [temperatures[1], pressures[1]],
            reference,
            2000.0,
        );
        let path = channel(&left, 3.0 * DIAMETER)
            .with_species_thermal_accommodation(
                vec![qualification("B", 0.75), qualification("A", 0.15)],
                accommodation_controls(),
                || false,
            )
            .unwrap();
        let full = DiffuseMolecularChannel::try_new(path.mouth().clone(), path.geometry().clone())
            .unwrap();
        let sampled = path.sample(&left, &right, AREA, tw, || false).unwrap();
        let full_sampled = full.sample(&left, &right, AREA, tw, || false).unwrap();
        assert_eq!(sampled.transmitted, full_sampled.transmitted);
        assert_eq!(
            sampled.wall.reflected_species,
            full_sampled.wall.reflected_species
        );
        assert!(path.thermal_accommodation().is_none());
        assert!(path.thermal_accommodation_provenance().is_none());
        let rows = path.species_thermal_accommodation().unwrap();
        assert_eq!(rows[0].qualification().species().as_str(), "A");
        assert_eq!(rows[1].qualification().species().as_str(), "B");
        assert_eq!(rows[0].qualification().coefficient(), 0.15);
        assert!(rows[0].qualification().provenance().contains("A/surface"));
        assert!(
            rows[0].accommodation().reflected_heat_probability()
                < rows[1].accommodation().reflected_heat_probability()
        );

        let mut gas_heat = [0.0; 2];
        let (mut vacuum_return, mut vacuum_escape) = (0.0, 0.0);
        for (index, mass, cp_r, slope) in [(0, 0.028, 3.5, 0.001), (1, 0.004, 2.5, 0.0)] {
            let response = rows[index].accommodation();
            let hr = response.reflected_heat_probability();
            let ht = response.transmitted_heat_probability();
            // Independent NASA polynomial difference: each incident mole's
            // crossing energy changes by R*dT*(cp/R - 1/2 + slope*(T+Tw)/2).
            let heat = [0, 1].map(|side| {
                let t = temperatures[side];
                incident_moles(pressures[side] * fractions[side][index], t, mass)
                    * R
                    * (tw - t)
                    * (cp_r - 0.5 + 0.5 * slope * (tw + t))
            });
            for side in 0..2 {
                gas_heat[side] += hr * heat[side] + ht * heat[1 - side];
            }
            vacuum_return += hr * heat[0];
            vacuum_escape += ht * heat[0];
        }
        for side in 0..2 {
            near(sampled.wall.heat_into_gas_w[side], gas_heat[side]);
        }
        near(sampled.wall.heat_into_wall_w, -gas_heat.iter().sum::<f64>());
        let vacuum = path.sample_to_vacuum(&left, AREA, tw, || false).unwrap();
        let full_vacuum = full.sample_to_vacuum(&left, AREA, tw, || false).unwrap();
        assert_eq!(vacuum.transmitted, full_vacuum.transmitted);
        assert_eq!(
            vacuum.wall.reflected_species,
            full_vacuum.wall.reflected_species
        );
        near(vacuum.wall.heat_into_gas_w, vacuum_return);
        near(vacuum.wall.heat_into_escaped_w, vacuum_escape);
        near(vacuum.wall.heat_into_wall_w, -vacuum_return - vacuum_escape);
        if let Some(old) = baseline {
            assert_eq!(sampled.wall, old);
        }
        baseline = Some(sampled.wall);
    }
}

#[test]
fn equal_species_coefficients_recover_common_law_and_builders_replace_modes() {
    let gas = mixture([0.4, 0.6], [400.0, 100_000.0], 0.0, 2000.0);
    let original = channel(&gas, 2.0 * DIAMETER);
    for alpha in [0.0, 0.25, 1.0] {
        let common = original
            .clone()
            .with_thermal_accommodation(
                alpha,
                accommodation_controls(),
                "synthetic common qualification",
                || false,
            )
            .unwrap();
        let by_species = original
            .clone()
            .with_species_thermal_accommodation(
                vec![qualification("B", alpha), qualification("A", alpha)],
                accommodation_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(common.geometry(), by_species.geometry());
        assert_eq!(
            common.sample(&gas, &gas, AREA, 900.0, || false).unwrap(),
            by_species
                .sample(&gas, &gas, AREA, 900.0, || false)
                .unwrap(),
        );
        let restored = by_species
            .with_thermal_accommodation(
                alpha,
                accommodation_controls(),
                "synthetic common qualification",
                || false,
            )
            .unwrap();
        assert!(restored.species_thermal_accommodation().is_none());
        assert_eq!(
            restored.thermal_accommodation().unwrap().coefficient(),
            alpha
        );
        let replaced = common
            .with_species_thermal_accommodation(
                vec![qualification("A", 0.0), qualification("B", 1.0)],
                accommodation_controls(),
                || false,
            )
            .unwrap();
        assert!(replaced.thermal_accommodation().is_none());
        assert_eq!(replaced.species_thermal_accommodation().unwrap().len(), 2);
    }
}

#[test]
fn species_wall_calorics_are_queried_only_for_present_nonelastic_species() {
    let gas = mixture([0.5, 0.5], [400.0, 100_000.0], 0.0, 600.0);
    let original = channel(&gas, DIAMETER);
    let elastic_a = original
        .clone()
        .with_species_thermal_accommodation(
            vec![qualification("A", 0.0), qualification("B", 1.0)],
            accommodation_controls(),
            || false,
        )
        .unwrap();
    // A is actually present, but elastic: its 600 K caloric limit must not
    // prevent B exchanging energy with a 900 K wall.
    elastic_a
        .admit_wall_for_parcel(&gas, 900.0, || false)
        .unwrap();
    assert!(
        elastic_a
            .sample(&gas, &gas, AREA, 900.0, || false)
            .unwrap()
            .wall
            .heat_into_gas_w[0]
            > 0.0
    );
    assert!(
        elastic_a
            .sample_to_vacuum(&gas, AREA, 900.0, || false)
            .unwrap()
            .wall
            .heat_into_gas_w
            > 0.0
    );
    assert!(
        elastic_a
            .admit_wall_for_parcel(&gas, f64::NAN, || false)
            .is_err()
    );

    let active_a = original
        .with_species_thermal_accommodation(
            vec![qualification("A", 1.0), qualification("B", 0.0)],
            accommodation_controls(),
            || false,
        )
        .unwrap();
    assert!(
        active_a
            .admit_wall_for_parcel(&gas, 900.0, || false)
            .is_err()
    );
    assert!(active_a.sample(&gas, &gas, AREA, 900.0, || false).is_err());
    assert!(
        active_a
            .sample_to_vacuum(&gas, AREA, 900.0, || false)
            .is_err()
    );

    let only_a = ClosedGasParcel::new(
        gas.model()
            .with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
            .unwrap(),
        point(400.0, 100_000.0),
        0.001,
    )
    .unwrap();
    // Absent B must not activate a wall temperature/source query either.
    elastic_a
        .admit_wall_for_parcel(&only_a, f64::NAN, || false)
        .unwrap();
    let closed = elastic_a
        .sample(&only_a, &only_a, AREA, f64::NAN, || false)
        .unwrap();
    let vacuum = elastic_a
        .sample_to_vacuum(&only_a, AREA, f64::NAN, || false)
        .unwrap();
    assert_eq!(closed.wall.heat_into_gas_w, [0.0; 2]);
    assert_eq!(vacuum.wall.heat_into_wall_w, 0.0);
    assert_eq!(closed.wall.reflected_species.len(), 1);
    assert!(closed.wall.reflected_species[0].mass_rate_kg_s[0] > 0.0);
}

#[test]
fn species_table_refuses_incomplete_foreign_or_substituted_original_sources() {
    let gas = mixture([0.5, 0.5], [400.0, 100_000.0], 0.0, 2000.0);
    let original = channel(&gas, DIAMETER);
    for rows in [
        vec![],
        vec![qualification("A", 0.2)],
        vec![qualification("A", 0.2), qualification("A", 0.8)],
        vec![qualification("A", 0.2), qualification("C", 0.8)],
        vec![
            qualification("A", 0.2),
            qualification("B", 0.8),
            qualification("C", 0.4),
        ],
    ] {
        assert!(
            original
                .clone()
                .with_species_thermal_accommodation(rows, accommodation_controls(), || false,)
                .is_err()
        );
    }
    for alpha in [-0.1, 1.1, f64::INFINITY, f64::NAN] {
        assert!(
            SpeciesThermalAccommodation::try_new(
                SpeciesId::new("A").unwrap(),
                alpha,
                "synthetic qualification",
            )
            .is_err()
        );
    }
    assert!(SpeciesThermalAccommodation::try_new(SpeciesId::new("A").unwrap(), 0.2, " ",).is_err());
    assert!(matches!(
        original.clone().with_species_thermal_accommodation(
            vec![qualification("A", 0.2), qualification("B", 0.8)],
            accommodation_controls(),
            || true,
        ),
        Err(GasEnergyError::Cancelled)
    ));
    let path = original
        .with_species_thermal_accommodation(
            vec![qualification("A", 0.0), qualification("B", 0.0)],
            accommodation_controls(),
            || false,
        )
        .unwrap();
    let substituted = mixture([0.5, 0.5], [400.0, 100_000.0], 1.0, 2000.0);
    // Elastic energy exchange does not waive immutable source identity.
    assert!(
        path.admit_wall_for_parcel(&substituted, f64::NAN, || false)
            .is_err()
    );
    assert!(
        path.sample(&gas, &substituted, AREA, f64::NAN, || false)
            .is_err()
    );
    assert!(
        path.sample_to_vacuum(&substituted, AREA, f64::NAN, || false)
            .is_err()
    );
}
