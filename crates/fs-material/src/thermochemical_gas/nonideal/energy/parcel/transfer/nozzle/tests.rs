use super::*;
use crate::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use crate::thermochemical_gas::nonideal::energy::tests::bound_with_reference;
use crate::thermochemical_gas::tests::point;
use crate::thermochemical_gas::energy::{GasEnergySolve, entropy::GasEntropySolve};
use crate::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

fn bound() -> ThermochemicalPengRobinsonModel {
    bound_with_reference(0.0)
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
        // The original source only admits T >= 200 K. At rho/rho0=.5,
        // this isentrope reaches 220.6928 K and brackets its sonic state.
        minimum_density_fraction: 0.5,
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
    // Independent 65-digit PR/NASA isentrope. The manufactured mixture has
    // x=[.75,.25], M=.022, Cp/R=3.25+.00075*T (not constant Cp).
    close(flow.throat().temperature_k(), 296.327_413_254_839_35, 2e-8);
    close(flow.throat().density_kg_m3(), 46.586_933_930_034_394, 2e-8);
    close(flow.velocity_m_s, 92.021_757_539_019_58, 2e-7);
    close(flow.mass_flow_kg_s, 0.003_429_609_230_876_759_5, 2e-11);
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
    // Same independent variable-Cp PR/NASA oracle, solving c^2=-2*Delta h
    // separately from pressure matching before selecting the sonic throat.
    close(flow.throat().temperature_k(), 246.138_609_574_328_88, 2e-8);
    close(flow.throat().pressure_pa(), 2_606_775.478_614_349_7, 5e-4);
    close(flow.throat().density_kg_m3(), 30.672_132_687_941_348, 2e-8);
    close(flow.velocity_m_s, 350.440_664_370_255_6, 2e-7);
    close(flow.mass_flow_kg_s, 0.008_599_010_045_451_84, 2e-11);
}

#[test]
fn formation_references_cannot_erase_choked_unchoked_or_small_pressure_drives() {
    let plain = parcel(bound(), 5e6, 0.01);
    for back_pressure in [4.8e6, 1e6, 5e6 - 1.0] {
        let expected = plain.nozzle_flow_to(back_pressure, nozzle(), controls(), || false).unwrap();
        assert!(expected.mass_flow_kg_s > 0.0);
        for reference in [-1e24, 1e24] {
            let shifted = parcel(bound_with_reference(reference), 5e6, 0.01);
            let actual = shifted.nozzle_flow_to(back_pressure, nozzle(), controls(), || false).unwrap();
            assert_eq!(actual.choked, expected.choked);
            assert_eq!(actual.velocity_m_s, expected.velocity_m_s);
            assert_eq!(actual.mass_flow_kg_s, expected.mass_flow_kg_s);
            assert_eq!(actual.mach_number, expected.mach_number);
            assert_eq!(actual.throat().density_kg_m3(), expected.throat().density_kg_m3());
            assert_eq!(actual.throat().temperature_k(), expected.throat().temperature_k());
            assert_eq!(actual.throat().pressure_pa(), expected.throat().pressure_pa());
            assert_eq!(actual.upstream_enthalpy_j_kg, shifted.state().specific_enthalpy_j_kg());
            // Absolute subtraction really is unusable at this reference,
            // including the represented 1-Pa drive next to the shut boundary.
            let naive_drop = actual.upstream_enthalpy_j_kg - actual.throat().specific_enthalpy_j_kg();
            let physical_drop = 0.5 * actual.velocity_m_s.powi(2);
            assert!((naive_drop - physical_drop).abs() > 0.5 * physical_drop);
        }
    }
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
