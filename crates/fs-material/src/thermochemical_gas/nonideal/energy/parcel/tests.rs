//! G1/G3/G4 vessel energy, material-view and atomic-update regressions.
//! The mixture is manufactured; these are not empirical gas-validation claims.

use super::*;
use crate::fluid_acoustics::FluidAcoustics;
use crate::fluid_transport::{FluidTransport, ThermalFluidTransport};
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::nonideal::PengRobinsonSpeciesBinding;
use crate::thermochemical_gas::tests::{model, point, species_with_energy_reference, transport};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

mod volume;

fn bound(source: ThermochemicalGasModel) -> ThermochemicalPengRobinsonModel {
    ThermochemicalPengRobinsonModel::new(
        source,
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
    .unwrap()
}

fn parcel() -> IsochoricRealGasParcel {
    IsochoricRealGasParcel::new(
        bound(model(vec![0.75, 0.25], 0.001)),
        point(300.0, 5.0e6).with("material-lot", 12.0).unwrap(),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 280.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1.0e-7,
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
fn g1_dense_vessel_heat_matches_independent_pr_energy_and_refreshes_all_views() {
    let mut vessel = parcel();
    let before = vessel.state().clone();
    let mass = vessel.mass_kg();
    // Independent 75-digit Decimal PR cubic/departure calculation, using
    // NASA cp/R = .75*(3.5+.001*T)+.25*2.5 and M=.022 kg/mol.
    // Delta u = 56632.82833319377 ideal + 1361.9215822505413 residual.
    near(vessel.density_kg_m3(), 47.889_003_167_633_395, 5.0e-11);
    near(
        before.specific_heat_cv_j_kg_k(),
        959.158_977_263_523_2,
        2.0e-9,
    );
    assert_ne!(mass, before.source().state().density * vessel.volume_m3());
    let heat_j = mass * 57_994.749_915_444_315;
    let report = vessel.add_heat(heat_j, controls(), || false).unwrap();
    let after = vessel.state();
    near(after.temperature_k(), 360.0, 2.0e-9);
    near(after.pressure_pa(), 6_226_134.941_447_334, 1.0e-4);
    near(after.density_kg_m3(), before.density_kg_m3(), 2.0e-12);
    assert_eq!(vessel.mass_kg(), mass);
    assert_eq!(vessel.volume_m3(), 0.01);
    near(
        report.stored_energy_change_j,
        heat_j,
        mass * controls().energy_tolerance_j_kg,
    );
    assert!(report.residual_j.abs() <= mass * controls().energy_tolerance_j_kg);
    // This vessel does not obey the ideal fixed-volume p/T=constant shortcut.
    assert!((after.pressure_pa() - 5.0e6 * 360.0 / 300.0).abs() > 200_000.0);
    let fluid = after.transport_state().unwrap();
    let thermal = after.thermal_transport_state().unwrap();
    let acoustic = after.acoustic_state().unwrap();
    assert_eq!(fluid.density_kg_m3(), after.density_kg_m3());
    assert_eq!(
        fluid.specific_heat_j_kg_k(),
        after.specific_heat_cp_j_kg_k()
    );
    assert_eq!(thermal.transport(), fluid);
    assert_eq!(thermal.temperature_k(), after.temperature_k());
    assert_eq!(acoustic.density_kg_m3(), after.density_kg_m3());
    assert_eq!(acoustic.sound_speed_m_s(), after.sound_speed_m_s());
    assert_eq!(
        acoustic.characteristic_impedance_pa_s_m(),
        after.density_kg_m3() * after.sound_speed_m_s()
    );
    assert!(
        fluid.dynamic_viscosity_pa_s() > before.transport_state().unwrap().dynamic_viscosity_pa_s()
    );
    assert_ne!(after.sound_speed_m_s(), before.sound_speed_m_s());
    assert_eq!(
        vessel.point().axes().get("pressure"),
        Some(&after.pressure_pa())
    );
    assert_eq!(vessel.point().axes().get("material-lot"), Some(&12.0));
}

#[test]
fn g3_heat_partition_and_cooling_preserve_mass_and_recover_the_source_state() {
    let mut whole = parcel();
    let initial = whole.clone();
    let mut partitioned = initial.clone();
    let heat_j = whole.mass_kg() * 40_000.0;
    let report = whole.add_heat(heat_j, controls(), || false).unwrap();
    let mut accumulated = 0.0;
    for _ in 0..8 {
        accumulated += partitioned
            .add_heat(heat_j / 8.0, controls(), || false)
            .unwrap()
            .stored_energy_change_j;
    }
    near(
        partitioned.state().temperature_k(),
        whole.state().temperature_k(),
        2.0e-8,
    );
    near(
        partitioned.state().pressure_pa(),
        whole.state().pressure_pa(),
        1.0e-3,
    );
    near(
        accumulated,
        report.stored_energy_change_j,
        9.0 * whole.mass_kg() * controls().energy_tolerance_j_kg,
    );
    assert_eq!(partitioned.mass_kg(), initial.mass_kg());
    assert_eq!(partitioned.volume_m3(), initial.volume_m3());
    assert_eq!(partitioned.density_kg_m3(), initial.density_kg_m3());
    let inverse = whole.add_heat(-heat_j, controls(), || false).unwrap();
    near(
        whole.state().temperature_k(),
        initial.state().temperature_k(),
        2.0e-8,
    );
    near(
        whole.state().pressure_pa(),
        initial.state().pressure_pa(),
        1.0e-3,
    );
    assert!(inverse.residual_j.abs() <= whole.mass_kg() * controls().energy_tolerance_j_kg);
    assert_eq!(whole.point().axes().get("material-lot"), Some(&12.0));
    assert_eq!(whole.model(), initial.model());
}

#[test]
fn g4_refusal_and_late_cancellation_leave_the_complete_parcel_unchanged() {
    let mut vessel = parcel();
    let saved = vessel.clone();
    let heat_j = vessel.mass_kg() * 40_000.0;
    for heat in [f64::NAN, f64::INFINITY, 1.0e12] {
        assert!(vessel.add_heat(heat, controls(), || false).is_err());
        assert_eq!(vessel, saved);
    }
    let limited = GasEnergySolve {
        maximum_iterations: 1,
        ..controls()
    };
    assert!(matches!(
        vessel.add_heat(heat_j, limited, || false),
        Err(GasEnergyError::NotConverged { .. })
    ));
    assert_eq!(vessel, saved);
    // Discover the final precommit cancellation point on the same deterministic
    // operation, then cancel there after all constitutive evaluation completed.
    let mut probe = vessel.clone();
    let mut total_polls = 0;
    probe
        .add_heat(heat_j, controls(), || {
            total_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        vessel.add_heat(heat_j, controls(), || {
            polls += 1;
            polls == total_polls
        }),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(vessel, saved);

    // The final EOS pressure exceeds the original transport card's 10 MPa
    // support. Thermodynamic convergence cannot authorize that source query.
    let mut dense = IsochoricRealGasParcel::new(
        bound(model(vec![0.75, 0.25], 0.001)),
        point(300.0, 9.0e6),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap();
    let saved_dense = dense.clone();
    let heat = dense.mass_kg()
        * dense
            .model()
            .specific_internal_energy_change_at_density(dense.density_kg_m3(), 300.0, 440.0)
            .unwrap();
    assert!(dense.add_heat(heat, controls(), || false).is_err());
    assert_eq!(dense, saved_dense);
}

#[test]
fn g3_small_heat_is_independent_of_large_formation_energy_references() {
    let source = |reference| {
        ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![
                    species_with_energy_reference("A", 0.028, 3.5, 0.001, reference),
                    species_with_energy_reference("B", 0.004, 2.5, 0.0, reference),
                ],
                Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
            )
            .unwrap(),
            vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
            None,
        )
        .unwrap()
    };
    let initialize = |reference| {
        IsochoricRealGasParcel::new(
            bound(source(reference)),
            point(300.0, 5.0e6),
            0.01,
            PengRobinsonPhaseV1::Vapor,
        )
        .unwrap()
    };
    let mut reference = initialize(0.0);
    let mut shifted = initialize(1.0e20);
    let heat_j = reference.mass_kg() * 1.0e-3;
    let first = reference.add_heat(heat_j, controls(), || false).unwrap();
    let second = shifted.add_heat(heat_j, controls(), || false).unwrap();
    assert!(reference.state().temperature_k() > 300.0);
    assert_eq!(
        reference.state().temperature_k(),
        shifted.state().temperature_k()
    );
    assert_eq!(
        reference.state().pressure_pa(),
        shifted.state().pressure_pa()
    );
    assert_eq!(first, second);
    assert!(second.residual_j.abs() <= shifted.mass_kg() * controls().energy_tolerance_j_kg);
}

#[test]
fn g0_zero_heat_preserves_state_and_does_not_bypass_invalid_controls() {
    let mut vessel = parcel();
    let saved = vessel.clone();
    let report = vessel.add_heat(0.0, controls(), || false).unwrap();
    assert_eq!(vessel, saved);
    assert_eq!(report.stored_energy_change_j, 0.0);
    assert_eq!(report.iterations, 0);
    assert!(
        vessel
            .add_heat(
                0.0,
                GasEnergySolve {
                    maximum_iterations: 0,
                    ..controls()
                },
                || false
            )
            .is_err()
    );
    assert_eq!(vessel, saved);
    assert!(matches!(
        vessel.add_heat(0.0, controls(), || true),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(vessel, saved);
}
