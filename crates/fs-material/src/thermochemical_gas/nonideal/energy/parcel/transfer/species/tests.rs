//! Actual upstream partial enthalpy, finite inventories and original PR storage.
use super::super::RealGasMassTransferSolve;
use super::*;
use crate::fluid_transport::FluidTransport;
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::nonideal::PengRobinsonSpeciesBinding;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_qty::{Amount, Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

fn source(x_a: f64, references: [f64; 2]) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.001, references[0]),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, references[1]),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![x_a, 1.0 - x_a]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        None,
    )
    .unwrap()
}

fn real_model(source: ThermochemicalGasModel, interaction: f64) -> ThermochemicalPengRobinsonModel {
    let bindings: Vec<_> = source
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            let (tc, pc, omega) = if species.species() == &SpeciesId::new("A").unwrap() {
                (190.564, 4.5992e6, 0.01142)
            } else {
                (150.687, 4.863e6, -0.00219)
            };
            PengRobinsonSpeciesBinding {
                species: species.species().clone(),
                critical_temperature: Temperature::new(tc),
                critical_pressure: Pressure::new(pc),
                acentric_factor: omega,
            }
        })
        .collect();
    let interactions = if bindings.len() == 1 {
        vec![]
    } else {
        vec![0.0, interaction, interaction, 0.0]
    };
    ThermochemicalPengRobinsonModel::new(source, bindings, interactions).unwrap()
}
fn at(
    model: ThermochemicalPengRobinsonModel,
    temperature: f64,
    pressure: f64,
    owner: f64,
) -> IsochoricRealGasParcel {
    IsochoricRealGasParcel::new(
        model,
        point(temperature, pressure).with("owner", owner).unwrap(),
        0.001,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}
fn pair(references: [f64; 2]) -> [IsochoricRealGasParcel; 2] {
    [
        at(real_model(source(0.75, references), 0.03), 300.0, 5e6, 11.0),
        at(
            real_model(source(0.35, references), 0.03),
            350.0,
            3.8e6,
            22.0,
        ),
    ]
}
fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn controls() -> GasSpeciesExchangeSolve {
    GasSpeciesExchangeSolve {
        gas: [energy(); 2],
        maximum_species_fraction: 0.25,
        mass_tolerance_kg: 1e-14,
        energy_tolerance_j: 1e-7,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.15e} != {b:.15e}, tolerance {tolerance:.3e}"
    );
}
fn row(name: &str, mass: f64) -> SpeciesMassTransfer {
    SpeciesMassTransfer {
        species: SpeciesId::new(name).unwrap(),
        mass_left_to_right_kg: mass,
    }
}
fn counterflow() -> [SpeciesMassTransfer; 2] {
    [row("A", 0.0002), row("B", -0.0002)]
}

// Independent conversion of the published mole inventory into extensive mass.
fn species_mass(gas: &IsochoricRealGasParcel, id: &SpeciesId) -> f64 {
    gas.model()
        .source()
        .thermodynamics()
        .components()
        .iter()
        .zip(gas.model().mechanics().mole_fractions())
        .find(|(species, _)| species.species() == id)
        .map_or(0.0, |(species, x)| {
            gas.mass_kg() * x * species.molar_mass().value()
                / gas.model().mechanics().molar_mass().value()
        })
}
fn partial_enthalpy(gas: &IsochoricRealGasParcel, id: &SpeciesId) -> f64 {
    let index = gas
        .model()
        .source()
        .thermodynamics()
        .components()
        .iter()
        .position(|species| species.species() == id)
        .unwrap();
    gas.model()
        .partial_calorics_at_density(gas.density_kg_m3(), gas.state().temperature_k())
        .unwrap()
        .partial_specific_enthalpies_j_kg()[index]
}

#[test]
fn dense_counterflow_transports_actual_partial_enthalpies_with_zero_net_mass() {
    let [left, right] = pair([0.0; 2]);
    let rows = counterflow();
    // Each complete mixture recovers its original bulk h as sum(Y_i hbar_i).
    // Nonideal departures must be material in this dense fixture, not silently
    // replaced by the pure NASA species enthalpy or by the mixture mean.
    for gas in [&left, &right] {
        let partial = gas
            .model()
            .partial_calorics_at_density(gas.density_kg_m3(), gas.state().temperature_k())
            .unwrap();
        assert!(std::ptr::eq(partial.model(), gas.model()));
        let mut weighted_h = 0.0;
        for (i, species) in gas
            .model()
            .source()
            .thermodynamics()
            .components()
            .iter()
            .enumerate()
        {
            let mass_fraction = species_mass(gas, species.species()) / gas.mass_kg();
            let h = partial.partial_specific_enthalpies_j_kg()[i];
            weighted_h += mass_fraction * h;
            let ideal = species
                .evaluate_properties(Temperature::new(gas.state().temperature_k()))
                .unwrap()
                .h()
                .value()
                / species.molar_mass().value();
            let departure = partial.departure().components()[i].residual_enthalpy_j_mol()
                / species.molar_mass().value();
            close(h - ideal, departure, 1e-7);
            assert!(departure.abs() > 1.0);
        }
        close(weighted_h, gas.state().specific_enthalpy_j_kg(), 2e-7);
    }
    let expected_heat: f64 = rows
        .iter()
        .enumerate()
        .map(|(i, stream)| {
            stream.mass_left_to_right_kg * partial_enthalpy([&left, &right][i], &stream.species)
        })
        .sum();
    let candidate = left
        .exchange_species_with(&right, &rows, controls(), || false)
        .unwrap();
    let step = candidate.step();
    let endpoints = [candidate.left(), candidate.right()];
    close(step.enthalpy_left_to_right_j, expected_heat, 1e-9);
    assert_eq!(step.net_mass_left_to_right_kg, 0.0);
    assert_eq!(step.species.len(), 2);
    for (i, account) in step.species.iter().enumerate() {
        assert_eq!(account.transfer, rows[i]);
        close(
            account.source_enthalpy_j_kg,
            partial_enthalpy([&left, &right][i], &rows[i].species),
            1e-8,
        );
        close(
            account.enthalpy_left_to_right_j,
            rows[i].mass_left_to_right_kg * account.source_enthalpy_j_kg,
            1e-9,
        );
    }
    for (side, (after, before)) in endpoints.into_iter().zip([&left, &right]).enumerate() {
        assert_eq!(after.mass_kg(), before.mass_kg());
        assert_eq!(after.volume_m3(), before.volume_m3());
        assert_eq!(after.density_kg_m3(), before.density_kg_m3());
        let sign = if side == 0 { -1.0 } else { 1.0 };
        let actual_u = after.mass_kg() * after.state().specific_internal_energy_j_kg()
            - before.mass_kg() * before.state().specific_internal_energy_j_kg();
        close(
            actual_u,
            sign * expected_heat,
            controls().energy_tolerance_j,
        );
        close(
            step.energy_change_j[side],
            actual_u,
            controls().energy_tolerance_j,
        );
        for stream in &rows {
            close(
                species_mass(after, &stream.species) - species_mass(before, &stream.species),
                sign * stream.mass_left_to_right_kg,
                controls().mass_tolerance_kg,
            );
        }
        assert!((after.state().temperature_k() - before.state().temperature_k()).abs() > 0.1);
        assert!((after.state().pressure_pa() - before.state().pressure_pa()).abs() > 100.0);
        let state = after.state();
        assert_eq!(state.density_kg_m3(), after.density_kg_m3());
        assert_eq!(
            state.transport_state().unwrap().density_kg_m3(),
            after.density_kg_m3()
        );
        assert_eq!(
            state.transport_state().unwrap().specific_heat_j_kg_k(),
            state.specific_heat_cp_j_kg_k()
        );
        assert!(state.specific_heat_cv_j_kg_k() > 0.0 && state.sound_speed_m_s() > 0.0);
        for component in state.source().components() {
            assert_eq!(
                component.query().axes().get("owner"),
                Some(&[11.0, 22.0][side])
            );
        }
    }
    assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(step.mass_residual_kg.abs() <= controls().mass_tolerance_kg);
    assert!(step.maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
}

#[test]
fn proportional_partial_species_streams_recover_the_original_bulk_transfer() {
    let [left, right] = pair([0.0; 2]);
    let dm = 0.001;
    let rows: Vec<_> = left
        .model()
        .source()
        .thermodynamics()
        .components()
        .iter()
        .map(|species| SpeciesMassTransfer {
            species: species.species().clone(),
            mass_left_to_right_kg: dm * species_mass(&left, species.species()) / left.mass_kg(),
        })
        .collect();
    let selected = left
        .exchange_species_with(&right, &rows, controls(), || false)
        .unwrap();
    let bulk = left
        .transfer_mass_to(
            &right,
            dm,
            RealGasMassTransferSolve {
                donor: energy(),
                receiver: energy(),
                mass_tolerance_kg: controls().mass_tolerance_kg,
                energy_tolerance_j: controls().energy_tolerance_j,
            },
            || false,
        )
        .unwrap();
    close(selected.step().net_mass_left_to_right_kg, dm, 1e-17);
    close(
        selected.step().enthalpy_left_to_right_j,
        bulk.step().transferred_enthalpy_j,
        1e-8,
    );
    for (a, b) in [selected.left(), selected.right()]
        .into_iter()
        .zip([bulk.donor(), bulk.receiver()])
    {
        close(a.mass_kg(), b.mass_kg(), 1e-14);
        close(a.state().temperature_k(), b.state().temperature_k(), 2e-8);
        close(a.state().pressure_pa(), b.state().pressure_pa(), 5e-4);
        for id in [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()] {
            close(species_mass(a, &id), species_mass(b, &id), 1e-14);
        }
    }
}

#[test]
fn balanced_counterflow_retains_accepted_density_across_an_inexact_mass_volume_roundtrip() {
    let [original_left, right] = pair([0.0; 2]);
    let density = original_left.density_kg_m3();
    // Choose a deterministic nearby volume that makes rho*V/V differ from rho.
    // Searching only scalar arithmetic avoids depending on the last EOS bit or
    // constructing a synthetic parcel with inconsistent private fields.
    let volume = (1..=128)
        .map(|i| 0.001 * (1.0 + f64::from(i) / 127.0))
        .find(|&volume| ((density * volume) / volume).to_bits() != density.to_bits())
        .expect("dense fixture must expose a non-exact mass/volume roundtrip");
    let left = IsochoricRealGasParcel::new(
        original_left.model().clone(),
        original_left.point().clone(),
        volume,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap();
    assert_eq!(left.density_kg_m3().to_bits(), density.to_bits());
    assert_ne!(
        (left.mass_kg() / left.volume_m3()).to_bits(),
        left.density_kg_m3().to_bits()
    );
    let rows = counterflow();
    let expected_heat = rows[0].mass_left_to_right_kg * partial_enthalpy(&left, &rows[0].species)
        + rows[1].mass_left_to_right_kg * partial_enthalpy(&right, &rows[1].species);
    let candidate = left
        .exchange_species_with(&right, &rows, controls(), || false)
        .unwrap();
    assert_eq!(candidate.step().net_mass_left_to_right_kg, 0.0);
    for (side, (after, before)) in [candidate.left(), candidate.right()]
        .into_iter()
        .zip([&left, &right])
        .enumerate()
    {
        let retained = before.density_kg_m3().to_bits();
        assert_eq!(after.mass_kg().to_bits(), before.mass_kg().to_bits());
        assert_eq!(after.volume_m3().to_bits(), before.volume_m3().to_bits());
        assert_eq!(after.density_kg_m3().to_bits(), retained);
        assert_eq!(after.state().density_kg_m3().to_bits(), retained);
        assert_eq!(
            after
                .state()
                .transport_state()
                .unwrap()
                .density_kg_m3()
                .to_bits(),
            retained
        );
        let actual_u = after.mass_kg() * after.state().specific_internal_energy_j_kg()
            - before.mass_kg() * before.state().specific_internal_energy_j_kg();
        let sign = if side == 0 { -1.0 } else { 1.0 };
        close(
            actual_u,
            sign * expected_heat,
            controls().energy_tolerance_j,
        );
        close(
            candidate.step().energy_change_j[side],
            actual_u,
            controls().energy_tolerance_j,
        );
    }
    assert!(candidate.step().energy_residual_j.abs() <= controls().energy_tolerance_j);
}

#[test]
fn independent_huge_species_datums_preserve_selective_exchange_physics() {
    let [left, right] = pair([0.0; 2]);
    let expected = left
        .exchange_species_with(&right, &counterflow(), controls(), || false)
        .unwrap();
    for references in [[1e24, -2e24], [-3e24, 4e24]] {
        let [left, right] = pair(references);
        let actual = left
            .exchange_species_with(&right, &counterflow(), controls(), || false)
            .unwrap();
        for (a, b) in [actual.left(), actual.right()]
            .into_iter()
            .zip([expected.left(), expected.right()])
        {
            assert_eq!(a.mass_kg(), b.mass_kg());
            assert_eq!(a.density_kg_m3(), b.density_kg_m3());
            assert_eq!(
                a.model().mechanics().mole_fractions(),
                b.model().mechanics().mole_fractions()
            );
            close(a.state().temperature_k(), b.state().temperature_k(), 2e-8);
            close(a.state().pressure_pa(), b.state().pressure_pa(), 5e-4);
        }
        // Absolute reported energies carry the chosen references; only the
        // datum-cancelled original residual decides local/global acceptance.
        assert!(actual.step().enthalpy_left_to_right_j.abs() > 1e19);
        assert!(actual.step().energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(actual.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    }
}

#[test]
fn stream_permutation_reversal_and_zero_exchange_preserve_original_owners() {
    let [left, right] = pair([0.0; 2]);
    let rows = counterflow();
    let expected = left
        .exchange_species_with(&right, &rows, controls(), || false)
        .unwrap();
    let shuffled = left
        .exchange_species_with(
            &right,
            &[rows[1].clone(), rows[0].clone()],
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(shuffled, expected);
    let reversed_rows: Vec<_> = rows
        .iter()
        .map(|stream| SpeciesMassTransfer {
            species: stream.species.clone(),
            mass_left_to_right_kg: -stream.mass_left_to_right_kg,
        })
        .collect();
    let reversed = right
        .exchange_species_with(&left, &reversed_rows, controls(), || false)
        .unwrap();
    for (a, b) in [reversed.left(), reversed.right()]
        .into_iter()
        .zip([expected.right(), expected.left()])
    {
        close(a.mass_kg(), b.mass_kg(), 1e-14);
        close(a.state().temperature_k(), b.state().temperature_k(), 2e-8);
        close(a.state().pressure_pa(), b.state().pressure_pa(), 5e-4);
        assert_eq!(a.model(), b.model());
    }
    close(
        reversed.step().enthalpy_left_to_right_j,
        -expected.step().enthalpy_left_to_right_j,
        1e-9,
    );
    for rows in [vec![], vec![row("A", 0.0), row("B", -0.0)]] {
        let zero = left
            .exchange_species_with(&right, &rows, controls(), || false)
            .unwrap();
        assert_eq!(zero.left(), &left);
        assert_eq!(zero.right(), &right);
        assert_eq!(zero.step().iterations, 0);
        assert!(zero.step().species.is_empty());
        assert_eq!(zero.step().energy_change_j, [0.0; 2]);
        let (a, b, step) = zero.into_parts();
        assert_eq!([a, b], [left.clone(), right.clone()]);
        assert_eq!(step.energy_residual_j, 0.0);
    }
}

#[test]
fn selective_composition_requires_original_complete_consistent_pr_sources() {
    let [left, right] = pair([0.0; 2]);
    for source in [
        real_model(right.model().source().clone(), 0.04),
        real_model(source(0.35, [10.0, 0.0]), 0.03),
    ] {
        let incompatible = at(source, 350.0, 3.8e6, 22.0);
        let saved = [left.clone(), incompatible.clone()];
        assert!(
            left.exchange_species_with(&incompatible, &counterflow(), controls(), || false)
                .is_err()
        );
        assert_eq!([&left, &incompatible], [&saved[0], &saved[1]]);
    }
    let base = source(0.75, [0.0; 2]);
    let pure_a = base
        .with_mole_amounts(&[Amount::new(1.0), Amount::new(0.0)])
        .unwrap();
    let pure_b = base
        .with_mole_amounts(&[Amount::new(0.0), Amount::new(1.0)])
        .unwrap();
    let right = at(real_model(pure_b, 0.03), 350.0, 3.8e6, 22.0);
    // The mixed left source knows AB and can introduce A into initially pure B.
    let known = left
        .exchange_species_with(&right, &counterflow(), controls(), || false)
        .unwrap();
    assert_eq!(known.right().model().mechanics().components().len(), 2);
    assert_eq!(
        known.right().model().mechanics().binary_interactions(),
        &[0.0, 0.03, 0.03, 0.0]
    );
    close(
        species_mass(known.right(), &SpeciesId::new("A").unwrap()),
        0.0002,
        1e-14,
    );
    // Two pure catalogs do not supply a cross interaction, even when each has
    // an empty matrix and its own diagonal defaults to zero.
    let left = at(real_model(pure_a, 0.03), 300.0, 5e6, 11.0);
    let saved = [left.clone(), right.clone()];
    assert!(matches!(
        left.exchange_species_with(&right, &counterflow(), controls(), || false),
        Err(GasEnergyError::Invalid(
            "mixed gas needs explicit Peng-Robinson interaction data for every new species pair"
        ))
    ));
    assert_eq!([left, right], saved);
}

#[test]
fn malformed_withdrawals_and_late_cancellation_publish_neither_parcel() {
    let [left, right] = pair([0.0; 2]);
    let saved = [left.clone(), right.clone()];
    for rows in [
        vec![row("A", species_mass(&left, &SpeciesId::new("A").unwrap()))],
        vec![row("unknown", 0.0001)],
        vec![row("A", 0.0001), row("A", 0.0002)],
        vec![row("A", f64::NAN)],
    ] {
        assert!(
            left.exchange_species_with(&right, &rows, controls(), || false)
                .is_err()
        );
        assert_eq!([&left, &right], [&saved[0], &saved[1]]);
    }
    let mut polls = 0;
    let expected = left
        .exchange_species_with(&right, &counterflow(), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        left.exchange_species_with(&right, &counterflow(), controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(seen, polls);
    let mut exhausted = controls();
    exhausted.gas[1].maximum_iterations = 1;
    assert!(matches!(
        left.exchange_species_with(&right, &counterflow(), exhausted, || false),
        Err(GasEnergyError::NotConverged { .. })
    ));
    assert_eq!([&left, &right], [&saved[0], &saved[1]]);
    assert_eq!(
        left.exchange_species_with(&right, &counterflow(), controls(), || false)
            .unwrap(),
        expected
    );
}
