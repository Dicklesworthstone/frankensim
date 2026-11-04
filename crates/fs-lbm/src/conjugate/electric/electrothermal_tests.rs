//! Independent lumped-resistance references for coupled contact heating.

use super::*;
use crate::conjugate::{
    EnergyConfig, FlowField, FluidProperties, SolidMaterial, ThermalFace, solve_energy,
};

fn fixture(alpha: f64, drive: ElectrodeDrive) -> (VoxelDomain, ElectricSetup, ThermalSetup) {
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let electric = ElectricSetup {
        conductors: vec![Some(
            Conductor::new(2.0).with_temperature_coefficient(alpha, 300.0),
        )],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 0],
                hi: [1, 1, 1],
                drive,
            },
            Electrode {
                lo: [2, 0, 0],
                hi: [3, 1, 1],
                drive: ElectrodeDrive::Voltage(0.0),
            },
        ],
    };
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    faces[1] = ThermalFace::Temperature(300.0);
    (domain, electric, ThermalSetup::new(faces))
}

fn heat(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    k: f64,
    gate: &CancelGate,
) -> Result<EnergySolution, ChtError> {
    solve_energy(
        domain,
        &FluidProperties {
            density_kg_m3: 1.0,
            specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: 1.0,
            kinematic_viscosity_m2_s: 1.0,
        },
        &[SolidMaterial::new("analytic resistor", k)],
        &FlowField::quiescent(domain),
        setup,
        &EnergyConfig {
            tolerance: 1e-13,
            ..EnergyConfig::default()
        },
        gate,
    )
}

fn contact() -> ElectricalContact {
    ElectricalContact {
        cells: [0, 1],
        resistance_ohm_m2: 4.0,
        heat_fraction_to_first: 0.25,
    }
}

#[test]
fn electrothermal_contact_heat_enters_electrodes_and_bulk_once() {
    let gate = CancelGate::new();
    let (domain, electric, thermal) = fixture(0.0, ElectrodeDrive::Current(1.0));
    let (temperature, field, report) = solve_electrothermal_with_contacts(
        &domain,
        &electric,
        &thermal,
        &[contact()],
        |source| heat(&domain, source, 1.0, &gate),
        &gate,
    )
    .unwrap();
    // One unit-area, unit-length bulk resistor (2 ohm), plus a 4 ohm
    // contact. At 1 A the electrode gets 1 W, and the bulk gets 2+3 W.
    assert_eq!(report.passes, 1);
    assert_eq!(report.joule_source_change_w, 0.0);
    for (actual, expected) in field.joule_w.iter().zip([1.0, 5.0, 0.0]) {
        assert!((actual - expected).abs() < 1e-12);
    }
    assert!((field.electrodes[0].voltage_v - 6.0).abs() < 1e-12);
    assert!((field.total_joule_w - 6.0).abs() < 1e-12);
    assert!((field.delivered_w - 6.0).abs() < 1e-12);
    assert!((temperature.report.balance.source_w - 6.0).abs() < 1e-12);
    // Solve [3,-1,0; -1,2,-1; 0,-1,3] theta = [1,5,0].
    for (actual, rise) in temperature
        .temperature
        .iter()
        .zip([5.0 / 3.0, 4.0, 4.0 / 3.0])
    {
        assert!((actual - (300.0 + rise)).abs() < 1e-10);
    }
}

#[test]
fn electrothermal_contact_feedback_matches_nonlinear_series_resistance() {
    let gate = CancelGate::new();
    let (domain, electric, thermal) = fixture(0.05, ElectrodeDrive::Voltage(3.0));
    let (temperature, field, report) = solve_electrothermal_with_contacts(
        &domain,
        &electric,
        &thermal,
        &[contact()],
        |source| heat(&domain, source, 1.0, &gate),
        &gate,
    )
    .unwrap();
    // R_bulk = 2 + 0.1 theta; I = 3 / (6 + 0.1 theta).
    // The central thermal inverse row is [1/4,3/4,1/4], so
    // theta = I^2 (4 + 0.075 theta). Bracket the independent cubic.
    let (mut lo, mut hi): (f64, f64) = (0.0, 4.0);
    for _ in 0..80 {
        let theta = 0.5 * (lo + hi);
        if theta * (6.0 + 0.1 * theta).powi(2) > 9.0 * (4.0 + 0.075 * theta) {
            hi = theta;
        } else {
            lo = theta;
        }
    }
    let rise = 0.5 * (lo + hi);
    let current = 3.0 / (6.0 + 0.1 * rise);
    assert!(report.passes > 2);
    assert!((temperature.temperature[1] - 300.0 - rise).abs() < 1e-9);
    assert!((field.electrodes[0].current_a - current).abs() < 1e-10);
    assert!((temperature.report.balance.source_w - field.total_joule_w).abs() < 1e-9);
    // Re-solving the electrical problem at the actual published temperature
    // must reproduce the published current and local heat, not a lagged state.
    let endpoint = electric
        .solve_with_contacts(
            &domain,
            &temperature.temperature,
            &[contact()],
            |_, _, t| Ok(electric.conductors[0].unwrap().resistivity(t)),
            &gate,
        )
        .unwrap();
    assert_eq!(field, endpoint);
}

#[test]
fn electrothermal_small_temperature_updates_cannot_hide_joule_mismatch() {
    let gate = CancelGate::new();
    let (domain, electric, thermal) = fixture(1e6, ElectrodeDrive::Current(0.1));
    let (temperature, field, report) = solve_electrothermal(
        &domain,
        &electric,
        &thermal,
        |source| heat(&domain, source, 1e6, &gate),
        &gate,
    )
    .unwrap();
    // A highly cooled, sensitive resistor: the second temperature update is
    // already below 1e-9 K while its Joule-source mismatch is still resolved.
    // theta = (3/4) I^2 rho0 (1 + alpha theta) / k.
    let source_rise = 0.75 * 0.01 * 2.0 / 1e6;
    let exact = 300.0 + source_rise / (1.0 - 1e6 * source_rise);
    assert!(report.passes > 2, "{report:?}");
    assert!((temperature.temperature[1] - exact).abs() < 1e-12);
    let endpoint =
        solve_electric(&domain, &electric, Some(&temperature.temperature), &gate).unwrap();
    assert_eq!(field, endpoint);
    let mismatch = (temperature.report.balance.source_w - field.total_joule_w).abs();
    assert!(mismatch <= 1e-12 + 1e-9 * field.total_joule_w, "{mismatch}");
    assert!((report.joule_source_change_w - mismatch).abs() < 1e-15);
}
