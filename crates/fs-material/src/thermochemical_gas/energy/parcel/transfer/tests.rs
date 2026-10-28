//! Synthetic source fixtures exercise actual material/energy consumers.
use super::*;
use crate::fluid_transport::FluidTransport;
use crate::thermochemical_gas::tests::{model, point, species, species_with_energy_reference, transport};

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn controls() -> GasMassTransferSolve {
    GasMassTransferSolve { donor: energy(), receiver: energy(),
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }
}
fn parcel(fraction: f64, slope: f64, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![fraction, 1.0-fraction], slope),
        point(t, p).with("local-lot", p).unwrap(), 0.001).unwrap()
}
fn pure(name: &str, molar_mass: f64, cp_r: f64, slope: f64, t: f64, p: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species(name, molar_mass, cp_r, slope)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport(name, 1.7e-5, 0.025)], None).unwrap();
    ClosedGasParcel::new(model, point(t, p), 0.001).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol*b.abs().max(1.0), "{a} != {b}");
}

pub(super) fn reference_parcel(x: f64, slope: f64, t: f64, p: f64, reference: [f64; 2]) -> ClosedGasParcel {
    let source = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A", 0.028, 3.5, slope, reference[0]),
            species_with_energy_reference("B", 0.004, 2.5, 0.0, reference[1])],
        Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0-x]).unwrap()).unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap();
    ClosedGasParcel::new(source, point(t, p), 0.001).unwrap()
}

#[test]
fn g1_bulk_filling_and_blowdown_retain_flow_work_under_arbitrary_formation_datums() {
    let ordinary = [reference_parcel(0.75, 0.0, 500.0, 200000.0, [0.0; 2]),
        reference_parcel(0.75, 0.0, 300.0, 100000.0, [0.0; 2])];
    let shifted = [reference_parcel(0.75, 0.0, 500.0, 200000.0, [1e24, -2e24]),
        reference_parcel(0.75, 0.0, 300.0, 100000.0, [1e24, -2e24])];
    let dm = ordinary[0].mass_kg() * 0.1;
    let expected = ordinary[0].transfer_mass_to(&ordinary[1], dm, controls(), || false).unwrap();
    let actual = shifted[0].transfer_mass_to(&shifted[1], dm, controls(), || false).unwrap();
    let gas = ordinary[0].state().state();
    let cv = gas.specific_heat_cp - gas.specific_gas_constant;
    let donor_t = 500.0 * (1.0 - dm / (ordinary[0].mass_kg() - dm) * gas.specific_gas_constant / cv);
    let receiver_t = (ordinary[1].mass_kg() * cv * 300.0 + dm * gas.specific_heat_cp * 500.0)
        / ((ordinary[1].mass_kg() + dm) * cv);
    near(actual.donor().state().state().temperature, donor_t, 1e-11);
    near(actual.receiver().state().state().temperature, receiver_t, 1e-11);
    for (a, b) in [(actual.donor(), expected.donor()), (actual.receiver(), expected.receiver())] {
        assert_eq!(a.state().state(), b.state().state());
        assert_eq!(a.mass_kg(), b.mass_kg());
    }
    assert_eq!(actual.step().energy_residual_j, expected.step().energy_residual_j);
    assert!(actual.step().energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert_eq!(shifted[0].state().specific_enthalpy_j_kg()
        - shifted[0].state().specific_internal_energy_j_kg(), 0.0,
        "the old absolute h-u path cannot resolve this physical flow work");
}

#[test]
fn withdrawal_uses_enthalpy_not_internal_energy_and_closes_both_volumes() {
    let donor = parcel(0.75, 0.0, 500.0, 200000.0);
    let receiver = parcel(0.75, 0.0, 300.0, 100000.0);
    let dm = 0.1*donor.mass_kg();
    let proposal = donor.transfer_mass_to(&receiver, dm, controls(), || false).unwrap();
    let gas = donor.state().state();
    let cv = gas.specific_heat_cp-gas.specific_gas_constant;
    let td = 500.0*(1.0-gas.specific_gas_constant/cv*dm/(donor.mass_kg()-dm));
    let tr = (receiver.mass_kg()*cv*300.0+dm*gas.specific_heat_cp*500.0)
        / ((receiver.mass_kg()+dm)*cv);
    near(proposal.donor().state().state().temperature, td, 1e-10);
    near(proposal.receiver().state().state().temperature, tr, 1e-10);
    assert!(proposal.donor().state().state().temperature < 500.0);
    near(proposal.step().transferred_enthalpy_j, dm*donor.state().specific_enthalpy_j_kg(), 1e-12);
    assert!((proposal.step().transferred_enthalpy_j-dm*donor.state().specific_internal_energy_j_kg()).abs() > 1.0);
    assert!(proposal.step().energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(proposal.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    assert_eq!(proposal.donor().volume_m3(), donor.volume_m3());
    assert_eq!(proposal.receiver().volume_m3(), receiver.volume_m3());
    assert_eq!(proposal.donor().point().axes().get("local-lot"), Some(&200000.0));
    assert_eq!(proposal.receiver().point().axes().get("local-lot"), Some(&100000.0));
}

#[test]
fn new_species_rebuilds_shared_thermal_transport_and_acoustic_medium_from_inventory() {
    let donor = pure("B", 0.004, 2.5, 0.0, 500.0, 250000.0);
    let receiver = pure("A", 0.028, 3.5, 0.002, 400.0, 100000.0);
    let dm = 0.1*donor.mass_kg();
    let trial = donor.transfer_mass_to(&receiver, dm, controls(), || false).unwrap();
    let next = trial.receiver();
    let n_a = receiver.mass_kg()/0.028;
    let n_b = dm/0.004;
    let thermo = next.parcel.model.thermodynamics();
    near(thermo.mole_composition().fractions()[0], n_a/(n_a+n_b), 1e-12);
    near(thermo.mole_composition().fractions()[1], n_b/(n_a+n_b), 1e-12);
    assert_eq!(&thermo.components()[0], &receiver.parcel.model.thermodynamics().components()[0]);
    assert_eq!(&thermo.components()[1], &donor.parcel.model.thermodynamics().components()[0]);
    let state = next.state().state();
    near(state.density, next.mass_kg()/next.volume_m3(), 1e-12);
    near(state.sound_speed.powi(2), state.gamma*state.pressure/state.density, 1e-12);
    near(next.state().transport_state().unwrap().density_kg_m3(), state.density, 1e-12);
    assert_ne!(state.sound_speed, receiver.state().state().sound_speed);
    assert_ne!(state.specific_heat_cp, receiver.state().state().specific_heat_cp);
    assert_ne!(next.state().components(), receiver.state().components());
    assert_eq!(trial.step().species.len(), 1);
    assert_eq!(trial.step().species[0].species, SpeciesId::new("B").unwrap());
    near(trial.step().species[0].mass_kg, dm, 1e-14);
}

#[test]
fn shared_species_sources_cannot_be_selected_by_name_or_silently_overwritten() {
    let donor = parcel(0.75, 0.001, 500.0, 200000.0);
    let conflict = parcel(0.25, 0.002, 400.0, 100000.0);
    let dm = donor.mass_kg()*0.01;
    assert!(donor.transfer_mass_to(&conflict, dm, controls(), || false).is_err());
    let mut recipient = parcel(0.25, 0.001, 400.0, 100000.0);
    let model = &recipient.parcel.model;
    let replaced = ThermochemicalGasModel::new(model.thermodynamics().clone(),
        vec![transport("A", 1.8e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap();
    recipient = ClosedGasParcel::new(replaced, point(400.0, 100000.0), 0.001).unwrap();
    assert!(donor.transfer_mass_to(&recipient, dm, controls(), || false).is_err());
    // No transfer does not claim an active mixed state or require compatibility.
    let noop = donor.transfer_mass_to(&conflict, 0.0, controls(), || false).unwrap();
    assert_eq!(noop.donor(), &donor);
    assert_eq!(noop.receiver(), &conflict);
    assert_eq!(noop.step().iterations, 0);
}

#[test]
fn donor_blowdown_refines_to_the_adiabatic_law_without_faking_exact_large_steps() {
    let initial = parcel(0.75, 0.0, 500.0, 300000.0);
    let reservoir = parcel(0.25, 0.0, 350.0, 100000.0);
    let exact = 500.0*0.8_f64.powf(initial.state().state().gamma-1.0);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let (mut a, mut b) = (initial.clone(), reservoir.clone());
        let dm = 0.2*initial.mass_kg()/f64::from(n);
        for _ in 0..n {
            let (next_a, next_b, report) = a.transfer_mass_to(&b, dm, controls(), || false).unwrap().into_parts();
            assert!(report.energy_residual_j.abs() <= controls().energy_tolerance_j);
            a = next_a; b = next_b;
        }
        errors.push((a.state().state().temperature-exact).abs());
        near(a.mass_kg(), 0.8*initial.mass_kg(), 1e-14);
        near(stored_change(&initial, &a)+stored_change(&reservoir, &b), 0.0, 1e-7);
    }
    assert!(errors[0]/errors[1] > 1.9);
    assert!(errors[1]/errors[2] > 1.9);
}

#[test]
fn transferred_endpoints_continue_through_existing_heat_work_and_force_owners() {
    let mut a = parcel(0.75, 0.001, 500.0, 200000.0);
    let b = parcel(0.25, 0.001, 400.0, 100000.0);
    let initial = a.clone();
    a.add_heat(5.0, energy(), || false).unwrap();
    let (donor, receiver, _) = a.transfer_mass_to(&b, 0.02*a.mass_kg(), controls(), || false).unwrap().into_parts();
    let motion = receiver.midpoint_heat_volume(0.0011, 3.0, energy(), || false).unwrap();
    let work = motion.step().work_by_gas_j;
    let load = motion.discrete_mechanical_load(0.0011, &[0.01], 100000.0).unwrap();
    near(load.generalized_forces().unwrap()[0], 0.01*(motion.midpoint().state().pressure-100000.0), 1e-12);
    let (next, _) = motion.into_parts();
    near(stored_change(&initial, &donor)+stored_change(&b, &next)+work, 8.0, 1e-8);
}

#[test]
fn late_source_budget_and_cancellation_failures_never_publish_half_a_transfer() {
    let donor = parcel(0.75, 0.001, 500.0, 200000.0);
    let receiver = parcel(0.25, 0.001, 400.0, 100000.0);
    let saved = (donor.clone(), receiver.clone());
    for mass in [-1.0, f64::NAN, f64::INFINITY, donor.mass_kg(), controls().mass_tolerance_kg*0.5] {
        assert!(donor.transfer_mass_to(&receiver, mass, controls(), || false).is_err());
    }
    let dm = 0.01*donor.mass_kg();
    let narrow = GasMassTransferSolve { receiver: GasEnergySolve { maximum_temperature_k: 250.0, ..energy() }, ..controls() };
    assert!(donor.transfer_mass_to(&receiver, dm, narrow, || false).is_err());
    let limited = GasMassTransferSolve { donor: GasEnergySolve { maximum_iterations: 1, ..energy() }, ..controls() };
    assert!(donor.transfer_mass_to(&receiver, dm, limited, || false).is_err());
    let mut total = 0;
    donor.transfer_mass_to(&receiver, dm, controls(), || { total += 1; false }).unwrap();
    for stop in 1..=total {
        let mut calls = 0;
        assert!(donor.transfer_mass_to(&receiver, dm, controls(), || { calls += 1; calls == stop }).is_err());
        assert_eq!((&donor, &receiver), (&saved.0, &saved.1));
    }
}
