use super::*;
use crate::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use crate::thermochemical_gas::tests::{model, point};
use crate::thermochemical_gas::energy::{GasEnergySolve, entropy::GasEntropySolve};
use crate::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
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
        minimum_temperature_k: 220.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn controls() -> RealGasNozzleStepSolve {
    RealGasNozzleStepSolve {
        entropy: GasEntropySolve {
            minimum_temperature_k: 220.0,
            maximum_temperature_k: 450.0,
            entropy_tolerance_j_kg_k: 1e-9,
            maximum_iterations: 128,
        },
        transfer: RealGasMassTransferSolve {
            donor: energy(),
            receiver: energy(),
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-7,
        },
        minimum_density_fraction: 0.2,
        pressure_tolerance_pa: 1e-5,
        sonic_tolerance_m2_s2: 1e-5,
        maximum_iterations: 128,
        maximum_donor_fraction: 0.1,
    }
}
fn nozzle() -> GasNozzle {
    GasNozzle { area_m2: 1e-6, discharge_coefficient: 0.8 }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol, "{a:.15e} != {b:.15e}, tol={tol:.3e}");
}

#[test]
fn modest_pressure_drop_resolves_an_unchoked_pr_throat() {
    let gas = parcel(bound(), 5.0e6, 0.01);
    let flow = gas.nozzle_flow_to(4.8e6, nozzle(), controls(), || false).unwrap();
    assert!(!flow.choked);
    assert!(flow.mass_flow_kg_s > 0.0);
    assert!(flow.velocity_m_s > 0.0);
    assert!(flow.mach_number > 0.0 && flow.mach_number < 1.0);
    close(flow.throat().pressure_pa(), 4.8e6, controls().pressure_tolerance_pa * 2.0);
    assert!(flow.entropy_residual_j_kg_k.abs()
        <= controls().entropy.entropy_tolerance_j_kg_k);
    assert_eq!(
        flow.upstream_enthalpy_j_kg,
        gas.state().specific_enthalpy_j_kg()
    );
}

#[test]
fn deep_backpressure_selects_the_source_backed_sonic_state() {
    let gas = parcel(bound(), 5.0e6, 0.01);
    let flow = gas.nozzle_flow_to(1.0e6, nozzle(), controls(), || false).unwrap();
    assert!(flow.choked);
    close(flow.mach_number, 1.0, 1e-5);
    assert!(flow.throat().pressure_pa() > 1.0e6);
    assert!(flow.throat().temperature_k() < gas.state().temperature_k());
    assert!(flow.throat().density_kg_m3() < gas.density_kg_m3());
}

#[test]
fn finite_reservoir_step_transfers_real_enthalpy_once_and_preserves_pair_energy() {
    let m = bound();
    let donor = parcel(m.clone(), 5.0e6, 0.01);
    let receiver = parcel(m, 4.6e6, 0.02);
    let (candidate, flow) = donor
        .through_nozzle_to(&receiver, nozzle(), 2e-5, controls(), || false)
        .unwrap();
    let step = candidate.step();
    close(
        step.transferred_mass_kg,
        flow.mass_flow_kg_s * 2e-5,
        controls().transfer.mass_tolerance_kg,
    );
    close(
        step.transferred_enthalpy_j,
        step.transferred_mass_kg * donor.state().specific_enthalpy_j_kg(),
        1e-9,
    );
    close(
        step.donor_energy_change_j + step.receiver_energy_change_j,
        step.energy_residual_j,
        1e-9,
    );
    assert!(step.energy_residual_j.abs() <= controls().transfer.energy_tolerance_j);
    assert!(candidate.donor().state().pressure_pa()
        >= candidate.receiver().state().pressure_pa());
}

#[test]
fn direction_crossing_work_limits_and_late_cancellation_refuse_atomically() {
    let m = bound();
    let donor = parcel(m.clone(), 5.0e6, 0.01);
    let receiver = parcel(m, 4.9e6, 0.02);
    assert!(donor.nozzle_flow_to(5.1e6, nozzle(), controls(), || false).is_err());
    assert!(donor
        .through_nozzle_to(&receiver, nozzle(), 1.0, controls(), || false)
        .is_err());

    let mut polls = 0;
    donor.nozzle_flow_to(4.8e6, nozzle(), controls(), || {
        polls += 1;
        false
    }).unwrap();
    assert!(polls > 2);
    let mut seen = 0;
    assert!(matches!(
        donor.nozzle_flow_to(4.8e6, nozzle(), controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(seen, polls);
}
