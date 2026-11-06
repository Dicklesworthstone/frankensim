//! G0/G1/G3/G4: reacting work, geometry tangents, reference invariance and retry.
use super::*;
use crate::thermochemical_gas::energy::equilibrium::tests::{fixture, isomers};
use crate::thermochemical_gas::tests::{
    point, species_with_energy_reference, transport, transport_in_domain,
};
use fs_thermochem::equilibrium::composition::EquilibriumSolveV1;

fn controls() -> EquilibriumGasParcelSolve {
    let mut equilibrium = EquilibriumSolveV1::default();
    equilibrium.relative_balance_tolerance = 1e-13;
    equilibrium.dimensionless_potential_tolerance = 1e-12;
    EquilibriumGasParcelSolve {
        thermochemistry: EquilibriumEntropySolveV1 {
            equilibrium,
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            entropy_tolerance_j_k: 1e-8,
            maximum_temperature_iterations: 128,
        },
        density: EquilibriumDensitySolveV1 {
            minimum_pressure_pa: 1e4,
            maximum_pressure_pa: 8e6,
            relative_density_tolerance: 1e-11,
            maximum_pressure_iterations: 96,
        },
        mass_tolerance_kg: 1e-12,
        energy_tolerance_j: 1e-5,
    }
}
fn feed(amount: f64) -> Vec<(SpeciesId, Amount)> {
    vec![(SpeciesId::new("A").unwrap(), Amount::new(amount))]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn fraction(t: f64) -> f64 {
    1.0 / (1.0 + (-24000.0 / (R * t)).exp())
}
fn mix(t: f64) -> f64 {
    let y = fraction(t);
    -y * y.ln() - (1.0 - y) * (1.0 - y).ln()
}
fn reference_temperature(ratio: f64) -> f64 {
    let (mut lo, mut hi) = (200.0_f64, 1500.0_f64);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        let residual = 3.0 * (t / 900.0).ln() + ratio.ln() + mix(t) - mix(900.0);
        if residual < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    0.5 * (lo + hi)
}
fn parcel<'a>(
    bank: &'a ThermochemicalGasModel,
    chemistry: &'a IdealGasEquilibriumV1,
) -> EquilibriumGasParcel<'a> {
    EquilibriumGasParcel::new(
        bank,
        chemistry,
        &point(300.0, 1e5).with("sample_id", 7.0).unwrap(),
        &feed(2.0),
        Temperature::new(900.0),
        Pressure::new(2.0 * R * 900.0 / 0.02),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_reacting_moving_parcel_matches_independent_isomer_work_and_round_trip() {
    let (bank, chemistry) = isomers();
    let mut gas = parcel(&bank, &chemistry);
    let v0 = gas.volume_m3();
    let mass = gas.mass_kg();
    let initial_fraction = gas.chamber().equilibrium().equilibrium().species()[1].mole_fraction;
    let mut prior_work = 0.0;
    for ratio in [0.8, 1.3, 0.5, 1.6, 1.0, 0.7, 1.0] {
        let step = gas
            .set_volume_isentropic(Volume::new(v0 * ratio), || false)
            .unwrap();
        let expected = reference_temperature(ratio);
        let expected_work = -2.0
            * (3.0 * R * (expected - 900.0) - 24000.0 * (fraction(expected) - fraction(900.0)));
        close(
            gas.chamber().parcel().state().state().temperature,
            expected,
            3e-6,
        );
        close(step.work_by_gas_since_reference_j, expected_work, 1e-4);
        close(step.work_by_gas_j, expected_work - prior_work, 2e-4);
        assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(step.reference_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(step.entropy_change_since_reference_j_k.abs() <= 1e-8);
        assert_eq!(gas.mass_kg(), mass);
        assert_eq!(gas.volume_m3(), v0 * ratio);
        assert_eq!(gas.feed(), feed(2.0));
        assert_eq!(gas.chamber().parcel().point().axes()["sample_id"], 7.0);
        assert_eq!(
            gas.chamber().parcel().model().transport_materials(),
            bank.transport_materials()
        );
        close(
            gas.chamber().parcel().state().state().pressure * gas.volume_m3(),
            2.0 * R * expected,
            1e-4,
        );
        if ratio == 0.5 {
            assert!(
                gas.chamber().equilibrium().equilibrium().species()[1].mole_fraction
                    < initial_fraction
            );
        }
        prior_work = expected_work;
    }
    close(gas.work_by_gas_since_reference_j(), 0.0, 1e-4);
}

#[test]
fn g1_reacting_geometric_tangent_matches_a_re_equilibrated_curved_volume() {
    let (bank, chemistry) = isomers();
    let gas = parcel(&bank, &chemistry);
    let v = gas.volume_m3();
    let p = gas.chamber().parcel().state().state().pressure;
    let outside = 2e5;
    let gradient = [v]; // V(q)=V0*exp(q), evaluated at q=0.
    let load = gas.mechanical_load(v, &gradient, outside).unwrap();
    let action = load.tangent_action(&[1.0], &[v]).unwrap()[0];
    let response = gas.pressure_response().unwrap();
    let y = fraction(900.0);
    let cv = 3.0 * R + 24000.0_f64.powi(2) * y * (1.0 - y) / (R * 900.0_f64.powi(2));
    close(response.bulk_modulus_pa, p * (1.0 + R / cv), 2e-5);
    assert!(
        response.bulk_modulus_pa
            < gas
                .chamber()
                .parcel()
                .pressure_response()
                .unwrap()
                .bulk_modulus_pa
    );
    close(
        load.generalized_forces().unwrap()[0],
        (p - outside) * v,
        1e-10,
    );
    close(
        load.mechanical_power_w(&[0.03]).unwrap(),
        (p - outside) * v * 0.03,
        1e-10,
    );
    let h = 2e-4_f64;
    let mut forces = [0.0; 2];
    let mut pressures = [0.0; 2];
    for (i, q) in [-h, h].into_iter().enumerate() {
        let mut perturbed = parcel(&bank, &chemistry);
        let volume = v * q.exp();
        perturbed
            .set_volume_isentropic(Volume::new(volume), || false)
            .unwrap();
        pressures[i] = perturbed.pressure_response().unwrap().pressure_pa;
        forces[i] = (pressures[i] - outside) * volume;
    }
    close((forces[1] - forces[0]) / (2.0 * h), action, 1e-3);
    let dp = (pressures[1] - pressures[0]) / (v * (h.exp() - (-h).exp()));
    close(dp / response.pressure_volume_derivative_pa_m3, 1.0, 2e-6);
    assert!(gas.mechanical_load(v * 1.01, &gradient, outside).is_err());
    assert!(
        gas.mechanical_load_with_limit(v, &gradient, outside, 0, || false)
            .is_err()
    );
}

#[test]
fn g0_g4_noop_cancellation_and_pressure_refusal_preserve_accepted_state() {
    let (bank, chemistry) = isomers();
    let mut gas = parcel(&bank, &chemistry);
    let initial = gas.chamber().clone();
    let sample = gas.sample();
    let v = gas.volume_m3();
    let noop = gas.set_volume_isentropic(Volume::new(v), || false).unwrap();
    assert_eq!(noop.work_by_gas_j, 0.0);
    assert_eq!(noop.pressure_queries, 0);
    assert_eq!(gas.chamber(), &initial);
    assert_eq!(gas.sample(), sample);
    let mut polls = 0;
    gas.set_volume_isentropic(Volume::new(v * 0.7), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut trial = parcel(&bank, &chemistry);
        let old = trial.chamber().clone();
        let mut count = 0;
        assert!(matches!(
            trial.set_volume_isentropic(Volume::new(v * 0.7), || {
                count += 1;
                count >= stop
            }),
            Err(Error::Cancelled)
        ));
        assert_eq!(trial.chamber(), &old);
        assert_eq!(trial.work_by_gas_since_reference_j(), 0.0);
        trial
            .set_volume_isentropic(Volume::new(v * 0.7), || false)
            .unwrap();
        close(
            trial.chamber().parcel().state().state().temperature,
            reference_temperature(0.7),
            3e-6,
        );
    }
    let accepted = gas.chamber().clone();
    let work = gas.work_by_gas_since_reference_j();
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e12] {
        assert!(
            gas.set_volume_isentropic(Volume::new(volume), || false)
                .is_err()
        );
        assert_eq!(gas.chamber(), &accepted);
        assert_eq!(gas.work_by_gas_since_reference_j(), work);
    }
    gas.set_volume_isentropic(Volume::new(v), || false).unwrap();
    close(gas.work_by_gas_since_reference_j(), 0.0, 1e-4);
}

#[test]
fn g0_transport_domain_escape_cannot_publish_a_partially_reacted_parcel() {
    let (original, chemistry) = isomers();
    let limited = ThermochemicalGasModel::new(
        original.thermodynamics().clone(),
        vec![
            transport_in_domain("A", 1.7e-5, 0.025, [850.0, 950.0], [1000.0, 1e7]),
            transport_in_domain("B", 2.0e-5, 0.15, [850.0, 950.0], [1000.0, 1e7]),
        ],
        None,
    )
    .unwrap();
    let mut gas = parcel(&limited, &chemistry);
    let before = gas.chamber().clone();
    let sample = gas.sample();
    assert!(matches!(
        gas.set_volume_isentropic(Volume::new(gas.volume_m3() * 0.5), || false),
        Err(Error::Material(_))
    ));
    assert_eq!(gas.chamber(), &before);
    assert_eq!(gas.sample(), sample);
    assert_eq!(gas.work_by_gas_since_reference_j(), 0.0);
    gas.set_volume_isentropic(Volume::new(gas.volume_m3() * 0.99), || false)
        .unwrap();
}

#[test]
fn g3_extensive_scaling_retains_pressure_temperature_and_work_per_mass() {
    let (bank, chemistry) = isomers();
    let mut baseline: Option<(f64, f64, f64)> = None;
    for scale in [0.25, 1.0, 4.0] {
        let mut c = controls();
        c.mass_tolerance_kg *= scale;
        c.energy_tolerance_j *= scale;
        c.thermochemistry.entropy_tolerance_j_k *= scale;
        let mut gas = EquilibriumGasParcel::new(
            &bank,
            &chemistry,
            &point(300.0, 1e5),
            &feed(2.0 * scale),
            Temperature::new(900.0),
            Pressure::new(2.0 * R * 900.0 / 0.02),
            c,
            || false,
        )
        .unwrap();
        let step = gas
            .set_volume_isentropic(Volume::new(gas.volume_m3() * 0.7), || false)
            .unwrap();
        let actual = (
            gas.sample().temperature.value(),
            gas.pressure_response().unwrap().pressure_pa,
            step.work_by_gas_j / gas.mass_kg(),
        );
        if let Some(reference) = baseline {
            close(actual.0, reference.0, 3e-6);
            close(actual.1 / reference.1, 1.0, 1e-8);
            close(actual.2 / reference.2, 1.0, 1e-7);
        } else {
            baseline = Some(actual);
        }
    }
}

#[test]
fn g3_small_inert_work_retains_large_common_formation_reference() {
    let ratio = 1.0 + 1e-5_f64;
    let expected = -2.0 * 3.0 * R * 900.0 * ((-ratio.ln() / 3.0).exp_m1());
    for offset in [0.0, 1e11] {
        let (bank, chemistry) = fixture(
            vec![species_with_energy_reference("A", 0.01, 4.0, 0.0, offset)],
            vec![transport("A", 1.7e-5, 0.025)],
            &["H"],
            vec![vec![1]],
        );
        let mut gas = parcel(&bank, &chemistry);
        let result = gas
            .set_volume_isentropic(Volume::new(gas.volume_m3() * ratio), || false)
            .unwrap();
        close(result.work_by_gas_j, expected, 2e-7);
        assert!(result.reference_energy_residual_j.abs() <= controls().energy_tolerance_j);
        gas.set_volume_isentropic(Volume::new(gas.volume_m3() / ratio), || false)
            .unwrap();
        close(gas.work_by_gas_since_reference_j(), 0.0, 2e-7);
    }
}
