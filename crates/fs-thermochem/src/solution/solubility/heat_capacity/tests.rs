//! G0/G1/G3: independent Kirchhoff/van't Hoff oracles and source-domain checks.

use super::super::{IdealSolubilityV1, R};
use super::*;
use crate::MolarThermalQuantityV1;
use fs_qty::AmountConcentration;

fn law(enthalpy: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("explicit-solvent").unwrap(),
        SpeciesId::new("neutral-solute").unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(2.0),
        enthalpy,
        [Temperature::new(250.0), Temperature::new(400.0)],
    )
    .unwrap()
}

fn curve(datum: f64, knots: &[(f64, f64)]) -> MolarHeatCapacityCurveV1 {
    MolarHeatCapacityCurveV1::new(
        SpeciesId::new("neutral-solute").unwrap(),
        Temperature::new(datum),
        knots
            .iter()
            .map(|&(t, cp)| (Temperature::new(t), MolarThermalQuantityV1::new(cp)))
            .collect(),
    )
    .unwrap()
}

fn constant(cp: f64) -> MolarHeatCapacityCurveV1 {
    curve(300.0, &[(250.0, cp), (400.0, cp)])
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-100),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn constant_phase_cp_difference_matches_integrated_kirchhoff_law() {
    let model = law(12_000.0)
        .with_phase_heat_capacities(constant(25.0), constant(85.0))
        .unwrap();
    for t in [250.0, 275.0, 300.0, 340.0, 400.0] {
        let delta = t - 300.0;
        relative(
            model.dissolution_enthalpy_at(Temperature::new(t)).unwrap(),
            12_000.0 + 60.0 * delta,
            3e-15,
        );
        let expected_log =
            12_000.0 / R * delta / 300.0 / t + 60.0 / R * ((t / 300.0).ln() - delta / t);
        relative(
            model
                .saturation_concentration(Temperature::new(t))
                .unwrap()
                .value(),
            2.0 * expected_log.exp(),
            2e-13,
        );
    }
    assert_eq!(model.dissolution_enthalpy_j_mol(), 12_000.0);
    assert_eq!(model.phase_heat_capacities().unwrap()[0], constant(25.0));
    assert_eq!(model.phase_heat_capacities().unwrap()[1], constant(85.0));
}

#[test]
fn direct_thermal_increments_keep_neighboring_changes_and_original_source_domain() {
    let model = law(1e16)
        .with_phase_heat_capacities(constant(25.0), constant(85.0))
        .unwrap();
    let from = Temperature::new(300.0);
    let to = Temperature::new(300.0_f64.next_up());
    assert_eq!(
        model.dissolution_enthalpy_at(from).unwrap(),
        model.dissolution_enthalpy_at(to).unwrap()
    );
    let expected = 60.0 * (to.value() - from.value());
    relative(
        model.dissolution_enthalpy_change_j_mol(from, to).unwrap(),
        expected,
        2e-15,
    );
    relative(
        model.dissolution_enthalpy_change_j_mol(to, from).unwrap(),
        -expected,
        2e-15,
    );
    let log_change = (1e16 / R) * ((to.value() - from.value()) / from.value() / to.value());
    relative(
        model.log_saturation_ratio(from, to).unwrap(),
        log_change,
        2e-15,
    );
    assert!(
        model
            .log_saturation_ratio(Temperature::new(249.0), to)
            .is_err()
    );
    assert!(
        model
            .dissolution_enthalpy_change_j_mol(from, Temperature::new(401.0))
            .is_err()
    );

    let ordinary = law(12_000.0)
        .with_phase_heat_capacities(constant(25.0), constant(85.0))
        .unwrap();
    for (a, b) in [(270.0, 350.0), (335.0, 290.0), (300.0, 300.0)] {
        let (a, b) = (Temperature::new(a), Temperature::new(b));
        let increment = ordinary.log_saturation_ratio(a, b).unwrap();
        let reverse = ordinary.log_saturation_ratio(b, a).unwrap();
        assert!((increment + reverse).abs() < 5e-14);
        relative(
            ordinary.saturation_concentration(a).unwrap().value() * increment.exp(),
            ordinary.saturation_concentration(b).unwrap().value(),
            2e-13,
        );
    }
}

#[test]
fn unequal_knots_and_independent_datums_preserve_linear_cp_analytic_solution() {
    // Cp_s=20+0.1T, Cp_d=-5+0.3T: ΔCp=-25+0.2T. Different knot locations
    // and sensible datums must neither re-fit the law nor shift ΔH_ref.
    let solid = curve(275.0, &[(250.0, 45.0), (285.0, 48.5), (400.0, 60.0)]);
    let dissolved = curve(
        375.0,
        &[(250.0, 70.0), (310.0, 88.0), (340.0, 97.0), (400.0, 115.0)],
    );
    let model = law(-250.0)
        .with_phase_heat_capacities(solid, dissolved)
        .unwrap();
    for t in [
        250.0, 267.0, 285.0, 299.0, 300.0, 310.0, 340.0, 371.0, 400.0,
    ] {
        let delta = t - 300.0;
        let expected_enthalpy = -250.0 + (-25.0 + 0.1 * (t + 300.0)) * delta;
        relative(
            model.dissolution_enthalpy_at(Temperature::new(t)).unwrap(),
            expected_enthalpy,
            2e-14,
        );
        let expected_log = -250.0 / R * delta / 300.0 / t
            - 25.0 / R * ((t / 300.0).ln() - delta / t)
            + 0.1 / R * delta * delta / t;
        relative(
            model
                .saturation_concentration(Temperature::new(t))
                .unwrap()
                .value(),
            2.0 * expected_log.exp(),
            2e-13,
        );
    }
}

#[test]
fn piecewise_phase_sources_match_independent_enthalpy_quadrature_and_temperature_derivative() {
    let solid = curve(
        300.0,
        &[(250.0, 25.0), (290.0, 40.0), (350.0, 20.0), (400.0, 60.0)],
    );
    let dissolved = curve(
        315.0,
        &[(250.0, 60.0), (270.0, 20.0), (315.0, 70.0), (400.0, 45.0)],
    );
    let model = law(150.0)
        .with_phase_heat_capacities(solid.clone(), dissolved.clone())
        .unwrap();
    let enthalpy = |t: f64| {
        150.0
            + dissolved
                .enthalpy_change_j_per_mol(Temperature::new(300.0), Temperature::new(t))
                .unwrap()
            - solid
                .enthalpy_change_j_per_mol(Temperature::new(300.0), Temperature::new(t))
                .unwrap()
    };
    let breaks: [f64; 7] = [250.0, 270.0, 290.0, 300.0, 315.0, 350.0, 400.0];
    for t in [250.0_f64, 280.0, 301.0, 330.0, 400.0] {
        // Independently integrate ΔH(T)/(R*T²) with Simpson quadrature, using
        // the existing ORIGINAL phase enthalpy integrals rather than the
        // prepared ΔCp curve or its analytic log kernels.
        let mut integral = 0.0;
        for limits in breaks.windows(2) {
            let a = limits[0].max(t.min(300.0));
            let b = limits[1].min(t.max(300.0));
            if b <= a {
                continue;
            }
            let h = (b - a) / 256.0;
            let f = |u| enthalpy(u) / R / u / u;
            let mut partial = f(a) + f(b);
            for i in 1..256 {
                partial += if i % 2 == 0 { 2.0 } else { 4.0 } * f(a + f64::from(i) * h);
            }
            integral += partial * h / 3.0;
        }
        if t < 300.0 {
            integral = -integral;
        }
        relative(
            model.dissolution_enthalpy_at(Temperature::new(t)).unwrap(),
            enthalpy(t),
            4e-14,
        );
        relative(
            model
                .saturation_concentration(Temperature::new(t))
                .unwrap()
                .value(),
            2.0 * integral.exp(),
            3e-12,
        );
    }
    for t in [265.0, 285.0, 305.0, 325.0, 365.0] {
        let h = 0.001;
        let log_c = |u| {
            model
                .saturation_concentration(Temperature::new(u))
                .unwrap()
                .value()
                .ln()
        };
        let derivative = (log_c(t + h) - log_c(t - h)) / (2.0 * h);
        relative(derivative, enthalpy(t) / R / t / t, 3e-8);
    }
}

#[test]
fn zero_reference_enthalpy_retains_adjacent_temperature_curvature() {
    // Large but finite declared Cp exposes cancellation that would be hidden
    // in ordinary concentration roundoff. Adjacent temperatures have a real
    // O((T-Tref)^2) contribution even though ΔH_ref is exactly zero.
    let model = law(0.0)
        .with_phase_heat_capacities(constant(1.0), constant(1e32))
        .unwrap();
    for t in [
        f64::from_bits(300.0_f64.to_bits() - 1),
        f64::from_bits(300.0_f64.to_bits() + 1),
    ] {
        let x = (t - 300.0) / 300.0;
        let curvature = x * x * (0.5 + x * (-2.0 / 3.0 + x * (0.75 - 0.8 * x)));
        let expected = 2.0 * (1e32 / R * curvature).exp();
        let actual = model
            .saturation_concentration(Temperature::new(t))
            .unwrap()
            .value();
        assert!(actual > 2.0);
        relative(actual, expected, 5e-15);
        relative(
            model.dissolution_enthalpy_at(Temperature::new(t)).unwrap(),
            1e32 * (t - 300.0),
            3e-15,
        );
    }
}

#[test]
fn zero_heat_capacity_difference_preserves_constant_law_exactly_and_shares_sources() {
    for reference_enthalpy in [0.0, 12_000.0, -12_000.0] {
        let original = law(reference_enthalpy);
        let solid = curve(275.0, &[(240.0, 35.0), (400.0, 35.0)]);
        let dissolved = curve(350.0, &[(250.0, 35.0), (290.0, 35.0), (425.0, 35.0)]);
        let model = original
            .clone()
            .with_phase_heat_capacities(solid, dissolved)
            .unwrap();
        let cloned = model.clone();
        assert!(core::ptr::eq(
            model.phase_heat_capacities().unwrap(),
            cloned.phase_heat_capacities().unwrap()
        ));
        assert!(model.heap_bytes() > original.heap_bytes());
        for t in [250.0, 280.0, 300.0, 310.0, 400.0] {
            assert_eq!(
                model.saturation_concentration(Temperature::new(t)).unwrap(),
                original
                    .saturation_concentration(Temperature::new(t))
                    .unwrap()
            );
            assert_eq!(
                model.dissolution_enthalpy_at(Temperature::new(t)).unwrap(),
                reference_enthalpy
            );
        }
    }
}

#[test]
fn phase_sources_refuse_mismatched_identity_truncated_domains_and_unrepresentable_heat() {
    let wrong_species = MolarHeatCapacityCurveV1::new(
        SpeciesId::new("different-solute").unwrap(),
        Temperature::new(300.0),
        vec![
            (Temperature::new(250.0), MolarThermalQuantityV1::new(30.0)),
            (Temperature::new(400.0), MolarThermalQuantityV1::new(30.0)),
        ],
    )
    .unwrap();
    assert!(
        law(100.0)
            .with_phase_heat_capacities(wrong_species, constant(40.0))
            .is_err()
    );
    for limited in [
        curve(300.0, &[(251.0, 30.0), (400.0, 30.0)]),
        curve(300.0, &[(250.0, 30.0), (399.0, 30.0)]),
    ] {
        assert!(
            law(100.0)
                .with_phase_heat_capacities(constant(20.0), limited)
                .is_err()
        );
    }
    let model = law(100.0)
        .with_phase_heat_capacities(constant(20.0), constant(40.0))
        .unwrap();
    for temperature in [249.0, 401.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            model.dissolution_enthalpy_at(Temperature::new(temperature)),
            Err(SolubilityErrorV1::OutsideTemperatureDomain)
        );
        assert_eq!(
            model.saturation_concentration(Temperature::new(temperature)),
            Err(SolubilityErrorV1::OutsideTemperatureDomain)
        );
    }
    let extreme = law(0.0)
        .with_phase_heat_capacities(constant(1.0), constant(1e307))
        .unwrap();
    assert!(matches!(
        extreme.dissolution_enthalpy_at(Temperature::new(400.0)),
        Err(SolubilityErrorV1::Unrepresentable(_))
    ));
}
