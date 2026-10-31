use super::*;
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use crate::thermochemical_gas::tests::{model, point};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

fn bound() -> ThermochemicalPengRobinsonModel {
    ThermochemicalPengRobinsonModel::new(
        model(vec![0.75, 0.25], 0.001),
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    ).unwrap()
}
fn parcel(model: ThermochemicalPengRobinsonModel, pressure: f64, volume: f64) -> IsochoricRealGasParcel {
    IsochoricRealGasParcel::new(
        model, point(300.0, pressure), volume, PengRobinsonPhaseV1::Vapor,
    ).unwrap()
}
fn energy() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn controls() -> RealGasMassTransferSolve {
    RealGasMassTransferSolve {
        donor: energy(),
        receiver: energy(),
        mass_tolerance_kg: 1e-14,
        energy_tolerance_j: 1e-7,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol, "{a:.15e} != {b:.15e}, tol={tol:.3e}");
}

#[test]
fn dense_fixed_volumes_transport_upstream_enthalpy_once_and_refresh_pressures() {
    let m = bound();
    let donor = parcel(m.clone(), 5.0e6, 0.01);
    let receiver = parcel(m, 3.0e6, 0.02);
    let dm = 0.04 * donor.mass_kg();
    let candidate = donor.transfer_mass_to(&receiver, dm, controls(), || false).unwrap();
    let step = candidate.step();
    close(candidate.donor().mass_kg(), donor.mass_kg()-dm, 1e-14);
    close(candidate.receiver().mass_kg(), receiver.mass_kg()+dm, 1e-14);
    assert_eq!(candidate.donor().volume_m3(), donor.volume_m3());
    assert_eq!(candidate.receiver().volume_m3(), receiver.volume_m3());
    close(step.transferred_enthalpy_j, dm * donor.state().specific_enthalpy_j_kg(), 1e-9);
    close(step.donor_energy_change_j + step.receiver_energy_change_j, step.energy_residual_j, 1e-9);
    assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(step.mass_residual_kg.abs() <= controls().mass_tolerance_kg);
    assert!(candidate.donor().state().pressure_pa() < donor.state().pressure_pa());
    assert!(candidate.receiver().state().pressure_pa() > receiver.state().pressure_pa());
}

#[test]
fn zero_transfer_is_identity_and_late_cancellation_is_atomic() {
    let m = bound();
    let donor = parcel(m.clone(), 5.0e6, 0.01);
    let receiver = parcel(m, 3.0e6, 0.02);
    let zero = donor.transfer_mass_to(&receiver, 0.0, controls(), || false).unwrap();
    assert_eq!(zero.donor(), &donor);
    assert_eq!(zero.receiver(), &receiver);
    assert_eq!(zero.step().iterations, 0);

    let dm = 0.02 * donor.mass_kg();
    let mut polls = 0;
    donor.transfer_mass_to(&receiver, dm, controls(), || { polls += 1; false }).unwrap();
    let mut seen = 0;
    assert!(matches!(
        donor.transfer_mass_to(&receiver, dm, controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(seen, polls);
}

#[test]
fn incompatible_pr_source_or_unrepresentable_transfer_refuses() {
    let m = bound();
    let donor = parcel(m.clone(), 5.0e6, 0.01);
    let receiver = parcel(m, 3.0e6, 0.02);
    let other_source: ThermochemicalGasModel = model(vec![0.75, 0.25], 0.002);
    let other = ThermochemicalPengRobinsonModel::new(
        other_source,
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    ).unwrap();
    let incompatible = parcel(other, 3.0e6, 0.02);
    assert!(donor.transfer_mass_to(&incompatible, 0.01*donor.mass_kg(), controls(), || false).is_err());
    assert!(donor.transfer_mass_to(&receiver, donor.mass_kg(), controls(), || false).is_err());
    assert!(donor.transfer_mass_to(&receiver, f64::NAN, controls(), || false).is_err());
}
