//! Independent dense-fluid endpoints, work and transactional volume changes.
//! The binary fixture is manufactured, not an empirical gas calibration.

use super::*;
use crate::fluid_acoustics::FluidAcoustics;
use crate::fluid_transport::ThermalFluidTransport;
use crate::thermochemical_gas::nonideal::PengRobinsonSpeciesBinding;
use crate::thermochemical_gas::tests::{model, point};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::SpeciesId;

fn parcel_at(pressure_pa: f64) -> ClosedRealGasParcel {
    let model = ThermochemicalPengRobinsonModel::new(
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
    )
    .unwrap();
    ClosedRealGasParcel::new(
        model,
        point(300.0, pressure_pa)
            .with("material-lot", 12.0)
            .unwrap(),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

fn energy_controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1.0e-7,
        maximum_iterations: 128,
    }
}

fn entropy_controls() -> GasEntropySolve {
    GasEntropySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 450.0,
        entropy_tolerance_j_kg_k: 1.0e-9,
        maximum_iterations: 128,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_compression_matches_independent_pr_work_and_shared_material_views() {
    let mut gas = parcel_at(5.0e6);
    let initial = gas.clone();
    let step = gas
        .set_volume_isentropic(0.008, entropy_controls(), || false)
        .unwrap();
    // Independent 75-digit PR Helmholtz/NASA integration and entropy root.
    // x=[.75,.25], M=.022, cp/R=3.25+.00075*T, k_AB=.03.
    near(gas.mass_kg(), 0.478_890_031_676_333_95, 5.0e-13);
    near(gas.state().temperature_k(), 331.606_449_235_535_6, 2.0e-9);
    near(gas.state().pressure_pa(), 6_981_183.016_655_789, 2.0e-4);
    near(step.work_by_gas_j, -11_791.630_070_073_763, 2.0e-6);
    assert_eq!(step.stored_energy_change_j, -step.work_by_gas_j);
    assert!(
        step.entropy_change_j_k.abs()
            <= gas.mass_kg() * entropy_controls().entropy_tolerance_j_kg_k
    );
    assert_eq!(gas.mass_kg(), initial.mass_kg());
    assert_eq!(gas.density_kg_m3(), gas.mass_kg() / gas.volume_m3());
    assert_eq!(gas.model(), initial.model());
    assert_eq!(gas.point().axes().get("material-lot"), Some(&12.0));
    assert_eq!(
        gas.point().axes().get("pressure"),
        Some(&gas.state().pressure_pa())
    );
    let thermal = gas.state().thermal_transport_state().unwrap();
    let acoustic = gas.state().acoustic_state().unwrap();
    assert_eq!(thermal.transport().density_kg_m3(), gas.density_kg_m3());
    assert_eq!(acoustic.density_kg_m3(), gas.density_kg_m3());
    assert_eq!(thermal.temperature_k(), gas.state().temperature_k());
    near(
        gas.model()
            .specific_entropy_change_between_densities(
                initial.density_kg_m3(),
                initial.state().temperature_k(),
                gas.density_kg_m3(),
                gas.state().temperature_k(),
            )
            .unwrap(),
        0.0,
        entropy_controls().entropy_tolerance_j_kg_k,
    );
}

#[test]
fn g1_joule_expansion_cools_without_inventing_heat_or_boundary_work() {
    let initial = parcel_at(5.0e6);
    let (expanded, energy) = initial
        .at_volume_internal_energy_increment(0.0125, 0.0, energy_controls(), || false)
        .unwrap();
    // The endpoint closure holds total internal energy; the consumer owns
    // the free-expansion protocol. An ideal caloric shortcut would leave T=300.
    near(
        expanded.state().temperature_k(),
        294.872_058_857_231_76,
        2.0e-9,
    );
    near(expanded.state().pressure_pa(), 3_975_006.042_948_88, 2.0e-4);
    near(
        energy,
        0.0,
        initial.mass_kg() * energy_controls().energy_tolerance_j_kg,
    );
    assert_eq!(expanded.mass_kg(), initial.mass_kg());
    assert_eq!(expanded.density_kg_m3(), initial.mass_kg() / 0.0125);
    let entropy = initial
        .model()
        .specific_entropy_change_between_densities(
            initial.density_kg_m3(),
            300.0,
            expanded.density_kg_m3(),
            expanded.state().temperature_k(),
        )
        .unwrap();
    near(entropy, 78.087_722_432_685_2, 2.0e-8);
    let (restored, residual) = expanded
        .at_volume_internal_energy_increment(0.01, 0.0, energy_controls(), || false)
        .unwrap();
    near(restored.state().temperature_k(), 300.0, 2.0e-9);
    near(restored.state().pressure_pa(), 5.0e6, 2.0e-4);
    near(
        residual,
        0.0,
        initial.mass_kg() * energy_controls().energy_tolerance_j_kg,
    );
}

#[test]
fn g3_partitioned_work_roundtrip_and_post_compression_heat_preserve_inventory() {
    let initial = parcel_at(5.0e6);
    let mut whole = initial.clone();
    let mut partitioned = initial.clone();
    let report = whole
        .set_volume_isentropic(0.008, entropy_controls(), || false)
        .unwrap();
    let mut work = 0.0;
    for volume in [0.0095, 0.009, 0.0085, 0.008] {
        work += partitioned
            .set_volume_isentropic(volume, entropy_controls(), || false)
            .unwrap()
            .work_by_gas_j;
    }
    near(
        partitioned.state().temperature_k(),
        whole.state().temperature_k(),
        1.0e-8,
    );
    near(work, report.work_by_gas_j, 5.0e-6);
    let returned = partitioned
        .set_volume_isentropic(0.01, entropy_controls(), || false)
        .unwrap();
    near(returned.work_by_gas_j + work, 0.0, 5.0e-6);
    near(partitioned.state().temperature_k(), 300.0, 1.0e-8);
    near(partitioned.state().pressure_pa(), 5.0e6, 5.0e-4);
    let density = whole.density_kg_m3();
    let volume = whole.volume_m3();
    let heat = whole.add_heat(100.0, energy_controls(), || false).unwrap();
    assert_eq!(whole.mass_kg(), initial.mass_kg());
    assert_eq!(whole.volume_m3(), volume);
    assert_eq!(whole.density_kg_m3(), density);
    assert!(whole.state().temperature_k() > 331.606);
    near(
        heat.stored_energy_change_j,
        100.0,
        whole.mass_kg() * energy_controls().energy_tolerance_j_kg,
    );
    let held = whole.clone().into_isochoric();
    assert_eq!(ClosedRealGasParcel::from_isochoric(held), whole);
}

#[test]
fn g4_volume_refusal_source_limits_and_late_cancellation_are_atomic() {
    let mut gas = parcel_at(5.0e6);
    let saved = gas.clone();
    for volume in [0.0, -0.01, f64::NAN, f64::INFINITY, 1.0e10] {
        assert!(
            gas.set_volume_isentropic(volume, entropy_controls(), || false)
                .is_err()
        );
        assert_eq!(gas, saved);
    }
    let limited = GasEntropySolve {
        maximum_iterations: 1,
        ..entropy_controls()
    };
    assert!(matches!(
        gas.set_volume_isentropic(0.008, limited, || false),
        Err(GasEntropyError::NotConverged { .. })
    ));
    assert_eq!(gas, saved);
    let mut probe = gas.clone();
    let mut total = 0;
    probe
        .set_volume_isentropic(0.008, entropy_controls(), || {
            total += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        gas.set_volume_isentropic(0.008, entropy_controls(), || {
            polls += 1;
            polls == total
        }),
        Err(GasEntropyError::Cancelled)
    ));
    assert_eq!(gas, saved);
    // The final state would exceed the original transport card's 10 MPa domain.
    let mut dense = parcel_at(9.0e6);
    let dense_saved = dense.clone();
    assert!(
        dense
            .set_volume_isentropic(0.008, entropy_controls(), || false)
            .is_err()
    );
    assert_eq!(dense, dense_saved);
    assert!(
        gas.at_volume_internal_energy_increment(0.008, 1.0e20, energy_controls(), || false)
            .is_err()
    );
    assert!(matches!(
        gas.at_volume_internal_energy_increment(0.0125, 0.0, energy_controls(), || true),
        Err(GasEnergyError::Cancelled)
    ));
}

#[test]
fn g0_unchanged_volume_is_exact_but_still_checks_controls_and_cancellation() {
    let mut gas = parcel_at(5.0e6);
    let saved = gas.clone();
    let step = gas
        .set_volume_isentropic(gas.volume_m3(), entropy_controls(), || false)
        .unwrap();
    assert_eq!(gas, saved);
    assert_eq!(step.work_by_gas_j, 0.0);
    assert_eq!(step.iterations, 0);
    let invalid = GasEntropySolve {
        maximum_iterations: 0,
        ..entropy_controls()
    };
    assert!(
        gas.set_volume_isentropic(gas.volume_m3(), invalid, || false)
            .is_err()
    );
    assert!(matches!(
        gas.set_volume_isentropic(gas.volume_m3(), entropy_controls(), || true),
        Err(GasEntropyError::Cancelled)
    ));
    assert_eq!(gas, saved);
    let (same, energy) = gas
        .at_volume_internal_energy_increment(gas.volume_m3(), 0.0, energy_controls(), || false)
        .unwrap();
    assert_eq!(same.density_kg_m3(), gas.density_kg_m3());
    assert_eq!(same.state().temperature_k(), gas.state().temperature_k());
    assert_eq!(energy, 0.0);
}
