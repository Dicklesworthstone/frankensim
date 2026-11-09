//! G0/G1: chemical work, equilibrium, finite kinetics and thermal consistency.

use super::*;
use crate::SpeciesId;
use crate::solution::solubility::kinetics::{
    PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticsV1,
};
use fs_qty::{Time, VolumetricFlowRate};

fn law(enthalpy: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("solvent").unwrap(),
        SpeciesId::new("neutral-solute").unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(2.0),
        enthalpy,
        [Temperature::new(250.0), Temperature::new(400.0)],
    )
    .unwrap()
}

fn change(law: &IdealSolubilityV1, a: f64, b: f64) -> f64 {
    law.solid_relative_free_energy_change_j(
        Temperature::new(300.0),
        Volume::new(1.0),
        Amount::new(a),
        Amount::new(b),
    )
    .unwrap()
}

#[test]
fn chemical_potential_follows_activity_and_retains_absent_limit() {
    let model = law(0.0);
    let chemical = |c| {
        model
            .solid_relative_chemical_potential_j_mol(
                Temperature::new(300.0),
                AmountConcentration::new(c),
            )
            .unwrap()
    };
    assert_eq!(chemical(0.0), None);
    assert_eq!(chemical(2.0), Some(0.0));
    for c in [1e-300_f64, 0.5, 1.0, 4.0, 1e300] {
        let expected = (R * 300.0) * (c.ln() - 2.0_f64.ln());
        assert!((chemical(c).unwrap() - expected).abs() < 3e-10);
    }
    let adjacent = 2.0_f64.next_up();
    let expected = (R * 300.0) * ((adjacent - 2.0) / 2.0).ln_1p();
    assert!((chemical(adjacent).unwrap() / expected - 1.0).abs() < 3e-15);
}

#[test]
fn potential_derivative_is_the_same_solid_relative_chemical_potential() {
    let model = law(17_000.0);
    let t = Temperature::new(330.0);
    let v = Volume::new(0.7);
    for n in [0.2, 1.5, 6.0] {
        let h = 1e-5 * n;
        let derivative = model
            .solid_relative_free_energy_change_j(t, v, Amount::new(n - h), Amount::new(n + h))
            .unwrap()
            / (2.0 * h);
        let chemical = model
            .solid_relative_chemical_potential_j_mol(t, AmountConcentration::new(n / v.value()))
            .unwrap()
            .unwrap();
        assert!((derivative - chemical).abs() < 2e-6);
    }
}

#[test]
fn phase_equilibrium_minimizes_the_correct_total_inventory_potential() {
    let model = law(0.0);
    for total in [0.5, 2.0, 5.0] {
        let equilibrium = model
            .equilibrate(
                Temperature::new(300.0),
                Volume::new(1.0),
                Amount::new(total),
            )
            .unwrap();
        for j in 0..=100 {
            let aqueous = total * f64::from(j) / 100.0;
            assert!(change(&model, aqueous, equilibrium.dissolved.value()) <= 1e-12);
        }
    }
    // The saturation reference cannot be cancelled as if aqueous moles were
    // conserved: dissolution can increase dissolved-only ideal mixing energy.
    let high_saturation = IdealSolubilityV1::new(
        model.solvent().clone(),
        model.solute().clone(),
        Temperature::new(300.0),
        AmountConcentration::new(100.0),
        0.0,
        model.temperature_domain(),
    )
    .unwrap();
    assert!(10.0 * (10.0_f64.ln() - 1.0) - (1.0_f64.ln() - 1.0) > 0.0);
    assert!(change(&high_saturation, 1.0, 10.0) < 0.0);
}

#[test]
fn tiny_saturation_departures_preserve_work_and_reversal() {
    let model = law(0.0);
    for endpoint in [2.0_f64.next_down(), 2.0_f64.next_up(), 2.0 + 1e-8] {
        let delta = endpoint - 2.0;
        let expected_leading = (R * 300.0) * (delta * delta / 4.0);
        let work = change(&model, 2.0, endpoint);
        assert!(work > 0.0);
        assert!((work / expected_leading - 1.0).abs() < 2e-8);
        assert_eq!(change(&model, endpoint, 2.0), -work);
    }
    let delta = 2.0_f64.powi(-18);
    // Independent series of the symmetric antiderivative difference:
    // F(2+d)-F(2-d) = -RT*d^3/12 + O(d^5).
    let expected = -(R * 300.0) * delta.powi(3) / 12.0;
    assert!((change(&model, 2.0 - delta, 2.0 + delta) / expected - 1.0).abs() < 2e-12);

    // The reduced molar work is below f64 range; the final joules are not.
    let capacity = 1e-293_f64;
    let next = capacity.next_up();
    let delta = next - capacity;
    let leading = ((R * 300.0) * delta) * (0.5 * (delta / capacity));
    let work = model
        .solid_relative_free_energy_change_j(
            Temperature::new(300.0),
            Volume::new(5e-294),
            Amount::new(capacity),
            Amount::new(next),
        )
        .unwrap();
    assert!(work > 0.0);
    assert_eq!(work, leading);
}

#[test]
fn fixed_area_and_particle_phase_updates_dissipate_the_same_potential() {
    let model = law(0.0);
    let fixed = SolubilityKineticsV1::new(
        model.clone(),
        VolumetricFlowRate::new(0.4),
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
    .unwrap();
    let particles = fixed
        .clone()
        .with_particle_surface(Amount::new(1.0))
        .unwrap();
    for kinetics in [&fixed, &particles] {
        for (aqueous, solid) in [(0.0, 5.0), (0.5, 0.1), (1.0, 4.0), (4.0, 1.0), (4.0, 0.0)] {
            let initial =
                SolubilityKineticStateV1::from_phases(Amount::new(aqueous), Amount::new(solid))
                    .unwrap();
            for duration in [0.001, 1.0, 100.0] {
                let final_state = kinetics
                    .advance(
                        Temperature::new(300.0),
                        Volume::new(1.0),
                        &initial,
                        Time::new(duration),
                    )
                    .unwrap();
                assert!(change(&model, aqueous, final_state.state().dissolved().value()) <= 1e-12);
            }
        }
    }
}

#[test]
fn chemical_work_obeys_the_original_gibbs_helmholtz_enthalpy() {
    let model = law(21_000.0);
    let reduced_mu = |t| {
        model
            .solid_relative_chemical_potential_j_mol(
                Temperature::new(t),
                AmountConcentration::new(0.3),
            )
            .unwrap()
            .unwrap()
            / t
    };
    let (t, h) = (320.0, 0.001);
    let inferred_enthalpy = -t * t * (reduced_mu(t + h) - reduced_mu(t - h)) / (2.0 * h);
    assert!((inferred_enthalpy / 21_000.0 - 1.0).abs() < 1e-9);
}

#[test]
fn zero_inventory_limits_scaling_and_domains_are_explicit() {
    let model = law(0.0);
    assert_eq!(change(&model, 0.0, 0.0), 0.0);
    assert_eq!(change(&model, 0.0, 2.0), -(R * 300.0) * 2.0);
    for scale in [1e-280, 1.0, 1e280] {
        let work = model
            .solid_relative_free_energy_change_j(
                Temperature::new(300.0),
                Volume::new(scale),
                Amount::new(0.3 * scale),
                Amount::new(1.1 * scale),
            )
            .unwrap();
        assert!((work / scale / change(&model, 0.3, 1.1) - 1.0).abs() < 3e-13);
    }
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .solid_relative_free_energy_j(
                    Temperature::new(300.0),
                    Volume::new(1.0),
                    Amount::new(invalid),
                )
                .is_err()
        );
    }
    assert!(
        model
            .solid_relative_free_energy_j(
                Temperature::new(300.0),
                Volume::new(0.0),
                Amount::new(0.0),
            )
            .is_err()
    );
    assert_eq!(
        model.solid_relative_free_energy_j(
            Temperature::new(401.0),
            Volume::new(1.0),
            Amount::new(1.0),
        ),
        Err(SolubilityErrorV1::OutsideTemperatureDomain)
    );
}
