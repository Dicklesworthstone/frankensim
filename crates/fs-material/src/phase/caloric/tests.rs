use super::*;
use crate::phase::{HeatCapacityKnot, PhaseIntervalDensityLaw, SolidLiquidPhase};
use fs_blake3::ContentHash;

fn source(cp: [f64; 2]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([0x96; 32]),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: cp[0],
                bulk_density_kg_m3: 2000.0,
            },
            HeatCapacityKnot {
                temperature_k: 500.0,
                specific_heat_capacity_j_kg_k: cp[1],
                bulk_density_kg_m3: 1000.0,
            },
        ],
        1.0e6,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn coarse_source_boundaries_retain_exact_capacity_inverse_tangent_and_density() {
    for cp in [[1000.0, 3000.0], [3000.0, 1000.0], [2000.0, 2000.0]] {
        let curve = source(cp);
        assert_eq!(curve.knots().len(), 2);
        assert_eq!(curve.sensible_heat_capacities_on_interval(0), Some(cp));
        assert_eq!(
            curve.sensible_heat_capacity_bounds_j_kg_k(0).unwrap(),
            [cp[0].min(cp[1]), cp[0].max(cp[1])]
        );
        assert_eq!(curve.sensible_heat_capacities_on_interval(1), None);
        for t in [300.0, 300.125, 329.125, 400.0, 474.75, 499.875, 500.0] {
            let dt = t - 300.0;
            let capacity = cp[0] + (cp[1] - cp[0]) * dt / 200.0;
            let h = (0.5 * cp[0] + 0.5 * capacity) * dt;
            let state = curve.values_at_specific_enthalpy(h).unwrap();
            close(state.temperature_k, t, 8.0 * f64::EPSILON * t);
            close(state.bulk_density_kg_m3, 2000.0 - 5.0 * dt, 2e-12);
            close(
                curve
                    .temperature_derivative_at_specific_enthalpy(h)
                    .unwrap(),
                capacity.recip(),
                8.0 * f64::EPSILON / capacity,
            );
            close(
                curve.specific_enthalpy_on_sensible_interval(0, t).unwrap(),
                h,
                8.0 * f64::EPSILON * h.abs().max(1.0),
            );
        }
        assert!(
            curve
                .specific_enthalpy_on_sensible_interval(0, 299.0)
                .is_err()
        );
        assert!(
            curve
                .specific_enthalpy_on_sensible_interval(1, 400.0)
                .is_err()
        );
        assert!(
            curve
                .specific_enthalpy_on_sensible_interval(usize::MAX, 400.0)
                .is_err()
        );
        assert!(
            curve
                .specific_enthalpy_on_sensible_interval(0, f64::NAN)
                .is_err()
        );
        assert!(curve.sensible_heat_capacity_bounds_j_kg_k(1).is_err());
        assert!(
            curve
                .sensible_heat_capacity_bounds_j_kg_k(usize::MAX)
                .is_err()
        );
    }
}

#[test]
fn caloric_identity_distinguishes_laws_with_identical_boundary_chords() {
    let varying = source([1000.0, 3000.0]);
    let constant = source([2000.0, 2000.0]);
    assert_eq!(varying.knots(), constant.knots());
    assert_ne!(varying.identity(), constant.identity());
    let authored = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        varying.material_card_identity(),
        SolidLiquidPhase::Solid,
        varying.knots().to_vec(),
    )
    .unwrap();
    assert_eq!(authored.sensible_heat_capacities_on_interval(0), None);
    assert_eq!(
        authored.sensible_heat_capacity_bounds_j_kg_k(0).unwrap(),
        [2000.0; 2]
    );
    assert_eq!(
        authored
            .temperature_at_specific_enthalpy(200_000.0)
            .unwrap(),
        400.0
    );
    assert_eq!(
        authored
            .specific_enthalpy_on_sensible_interval(0, 400.0)
            .unwrap(),
        200_000.0
    );
    assert!(varying.temperature_at_specific_enthalpy(200_000.0).unwrap() > 420.0);
    let mut unresolved = authored.knots().to_vec();
    unresolved[1].temperature_k = 300.0_f64.next_up();
    unresolved[1].specific_enthalpy_j_kg = 1e300;
    let unresolved = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        authored.material_card_identity(),
        SolidLiquidPhase::Solid,
        unresolved,
    )
    .unwrap();
    assert!(unresolved.sensible_heat_capacity_bounds_j_kg_k(0).is_err());
}

#[test]
fn density_mixing_remainder_includes_the_original_capacity_curvature() {
    let curve = source([1000.0, 3000.0]);
    assert_eq!(
        curve.density_law_on_interval(0).unwrap(),
        PhaseIntervalDensityLaw::AffineTemperatureDensity
    );
    for (receiver_t, donor_t) in [(335.0, 470.0), (470.0, 335.0)] {
        let hr = curve
            .specific_enthalpy_on_sensible_interval(0, receiver_t)
            .unwrap();
        let hd = curve
            .specific_enthalpy_on_sensible_interval(0, donor_t)
            .unwrap();
        let rho_r = 2000.0 - 5.0 * (receiver_t - 300.0);
        let rho_d = 2000.0 - 5.0 * (donor_t - 300.0);
        let cp_r = 1000.0 + 10.0 * (receiver_t - 300.0);
        let slope = 5.0 / (rho_r * rho_r * cp_r);
        let expected = rho_r.recip() + slope * (hd - hr) - rho_d.recip();
        let (actual, intervals) = curve
            .mixing_specific_volume_m3_kg(hr, hd, 1, || false)
            .unwrap();
        close(actual, expected, 4e-19);
        assert_eq!(intervals, 1);
        assert!(
            curve
                .mixing_specific_volume_m3_kg(hr, hd, 0, || false)
                .is_err()
        );
    }
    assert_eq!(
        curve
            .mixing_specific_volume_m3_kg(0.0, 0.0, 0, || false)
            .unwrap(),
        (0.0, 0)
    );
    assert!(
        curve
            .mixing_specific_volume_m3_kg(0.0, 10.0, 1, || true)
            .is_err()
    );
}

#[test]
fn sensible_capacities_do_not_smooth_or_double_count_a_latent_plateau() {
    let knot = |h, t, f, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: rho,
    };
    let mut curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x97; 32]),
        vec![
            knot(0.0, 300.0, 0.0, 1000.0),
            knot(100_000.0, 400.0, 0.0, 900.0),
            knot(300_000.0, 400.0, 1.0, 800.0),
            knot(600_000.0, 500.0, 1.0, 700.0),
        ],
    )
    .unwrap();
    curve.retain_sensible_density_law();
    curve
        .retain_sensible_heat_capacities(vec![Some([800.0, 1200.0]), None, Some([2000.0, 4000.0])])
        .unwrap();
    assert_eq!(curve.sensible_heat_capacities_on_interval(1), None);
    assert!(curve.sensible_heat_capacity_bounds_j_kg_k(1).is_err());
    let solid = curve.values_at_specific_enthalpy(45_000.0).unwrap();
    assert_eq!(solid.temperature_k, 350.0);
    assert_eq!(solid.bulk_density_kg_m3, 950.0);
    let mixed = curve.values_at_specific_enthalpy(200_000.0).unwrap();
    assert_eq!(mixed.temperature_k, 400.0);
    assert_eq!(mixed.liquid_mass_fraction, 0.5);
    close(
        mixed.bulk_density_kg_m3,
        (0.5 / 900.0 + 0.5 / 800.0_f64).recip(),
        2e-13,
    );
    assert_eq!(
        curve
            .temperature_derivative_at_specific_enthalpy(100_000.0)
            .unwrap(),
        0.0
    );
    assert_eq!(
        curve
            .temperature_derivative_at_specific_enthalpy(300_000.0)
            .unwrap(),
        1.0 / 2000.0
    );
    assert_eq!(
        curve.temperature_at_specific_enthalpy(425_000.0).unwrap(),
        450.0
    );
    assert!(
        curve
            .specific_enthalpy_on_sensible_interval(1, 400.0)
            .is_err()
    );
}
