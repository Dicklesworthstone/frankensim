//! Electron-impact ionization and three-body recombination on the existing
//! finite-rate reactor. The forward rate is explicitly synthetic: this checks
//! numerical kinetics, conservation and detailed balance, not a measured
//! hydrogen ionization time or a non-Maxwellian discharge model.
//!
//! Ideal monatomic source entropy follows Sackur--Tetrode; the resulting
//! equilibrium is independently given by Saha. Ionization energy is NIST H I,
//! 13.598434599702 eV, <https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=H+I>.
//! NASA RP-1311 section 3.7 describes the ideal charge-balanced equilibrium:
//! <https://ntrs.nasa.gov/citations/19950013764>.

use super::*;

const K_B: f64 = 1.380_649e-23;
const PLANCK: f64 = 6.626_070_15e-34;
const N_A: f64 = 6.022_140_76e23;
const M_E: f64 = 9.109_383_713_9e-31;
const M_P: f64 = 1.672_621_925_95e-27;
const M_H: f64 = M_P + M_E;
const CHI: f64 = 13.598_434_599_702 * 1.602_176_634e-19;
const CHI_MOL: f64 = CHI * N_A;
const FORWARD_RATE: f64 = 20.0;

fn atom(name: &str, mass: f64, degeneracy: f64, formation: f64) -> Nasa9StandardStateModelV1 {
    let entropy = 2.5
        + (degeneracy * K_B / 100000.0
            * (2.0 * core::f64::consts::PI * mass * K_B / PLANCK.powi(2)).powf(1.5))
        .ln();
    let model = source(name, mass * N_A, 2.5, formation * N_A, entropy);
    let mut card = model.regions()[0].card().clone();
    card.validity = ValidityDomain::unconstrained().with("T", 1000.0, 30000.0);
    card.provenance.source =
        "ideal monatomic statistical mechanics; NIST H ionization energy".into();
    Nasa9StandardStateModelV1::new(
        model.species().clone(),
        model.molar_mass(),
        model.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn hydrogen_ionization() -> IdealGasKineticsV1 {
    let ids = vec![id("H"), id("H+"), id("e-")];
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        vec![
            atom("H", M_H, 2.0, 0.0),
            atom("H+", M_P, 1.0, CHI),
            atom("e-", M_E, 2.0, 0.0),
        ],
    )
    .unwrap();
    IdealGasKineticsV1::new(
        source,
        "synthetic-electron-impact-rate",
        1,
        vec![reaction(
            "H-ionization",
            &[("H", 1), ("e-", 1)],
            &[("H+", 1), ("e-", 2)],
            FORWARD_RATE,
        )],
        1e-12,
    )
    .unwrap()
}

// Saha in mol/m3, including the actual translational-mass correction. One
// mole of nuclei occupies one cubic metre throughout these tests.
fn saha_molar(t: f64) -> f64 {
    (2.0 * core::f64::consts::PI * (M_E * M_P / M_H) * K_B * t / PLANCK.powi(2)).powf(1.5)
        * (-CHI / (K_B * t)).exp()
        / N_A
}

fn equilibrium_fraction(t: f64) -> f64 {
    let r = saha_molar(t);
    2.0 * r.sqrt() / ((r + 4.0).sqrt() + r.sqrt())
}

fn energy(t: f64, x: f64) -> f64 {
    1.5 * R * t * (1.0 + x) + x * CHI_MOL
}

fn temperature(u: f64, x: f64) -> f64 {
    (u - x * CHI_MOL) / (1.5 * R * (1.0 + x))
}

fn feed(x: f64) -> Vec<(SpeciesId, Amount)> {
    vec![
        (id("H"), Amount::new(1.0 - x)),
        (id("H+"), Amount::new(x)),
        (id("e-"), Amount::new(x)),
    ]
}

fn controls() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        minimum_temperature_k: 1000.0,
        maximum_temperature_k: 30000.0,
        energy_tolerance_j: 1e-6,
        relative_tolerance: 1e-5,
        absolute_amount_tolerance_mol: 1e-12,
        entropy_tolerance_j_k: 1e-8,
        initial_step_s: 1e-4,
        maximum_step_s: 0.02,
        maximum_steps: 10000,
        ..IsochoricKineticsSolve::default()
    }
}

fn amount(result: &IsochoricKineticsResult, name: &str) -> f64 {
    let index = result
        .mechanism()
        .thermodynamics()
        .elemental_matrix()
        .species()
        .binary_search(&id(name))
        .unwrap();
    result.amounts()[index].value()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_ionization_forward_reverse_rates_balance_at_independent_saha_state() {
    let mechanism = hydrogen_ionization();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    for t in [7000.0, 12000.0, 20000.0] {
        let x = equilibrium_fraction(t);
        let state = mechanism
            .evaluate(Temperature::new(t), &[1.0 - x, x, x], &mut work, || false)
            .unwrap();
        relative(
            work.log_equilibrium_constants()[0].exp(),
            saha_molar(t),
            2e-12,
        );
        let forward = FORWARD_RATE * (1.0 - x) * x;
        assert!(work.reaction_rates()[0].abs() <= 2e-12 * forward);
        assert!(state.adiabatic_temperature_rate_k_s.abs() <= 2e-12 * forward * CHI_MOL / R);
    }
}

#[test]
fn g0_ionization_collisions_do_not_invent_an_electron_seed() {
    let mechanism = hydrogen_ionization();
    let result = advance_adiabatic(
        &mechanism,
        &[(id("H"), Amount::new(1.0))],
        Volume::new(1.0),
        Temperature::new(16000.0),
        Time::new(2.0),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(amount(&result, "H"), 1.0);
    assert_eq!(amount(&result, "H+"), 0.0);
    assert_eq!(amount(&result, "e-"), 0.0);
    relative(result.temperature().value(), 16000.0, 1e-12);
    assert_eq!(result.balance().entropy_change_j_k, 0.0);
}

#[test]
fn g1_ionization_reactor_freezes_actual_endpoint_without_equilibrium_projection() {
    let mechanism = hydrogen_ionization();
    let initial = [Amount::new(0.99), Amount::new(0.01), Amount::new(0.01)];
    let mut reactor = IsochoricReactor::new(mechanism, controls(), || false).unwrap();
    let endpoint = reactor
        .advance_adiabatic(
            &initial,
            Volume::new(1.0),
            Temperature::new(12000.0),
            Time::new(1e-4),
            || false,
        )
        .unwrap();
    let frozen = endpoint.frozen_state(1e-10, || false).unwrap();
    let owned = endpoint.to_owned(|| false).unwrap();
    let owned_frozen = owned.frozen_state(1e-10, || false).unwrap();
    assert_eq!(frozen, owned_frozen);
    assert_eq!(frozen.amounts(), owned.amounts());
    assert_eq!(frozen.heavy_temperature(), owned.temperature());
    assert_eq!(frozen.electron_temperature(), owned.temperature());
    relative(frozen.pressure().value(), owned.pressure().value(), 3e-15);
    relative(
        frozen.internal_energy().value(),
        owned.internal_energy().value(),
        3e-15,
    );
    relative(frozen.mass().value(), owned.mass().value(), 3e-15);
    // A deliberately short finite-rate interval stays far from Saha. Freezing
    // must retain its small electron density rather than replacing it by LTE.
    assert!(frozen.amounts()[2].value() < 0.02);
    assert!(equilibrium_fraction(owned.temperature().value()) > 5.0 * frozen.amounts()[2].value());
}

#[test]
fn g1_ionization_finite_rate_reactor_conserves_energy_mass_and_charge() {
    let mechanism = hydrogen_ionization();
    let initial_x = 0.01;
    let initial_t = 16000.0;
    let u = energy(initial_t, initial_x);
    let elapsed = 0.01;
    // Independent scalar reduction: charge and nuclei conservation leave x as
    // the only degree of freedom; the first law gives T(x) explicitly. Fixed
    // small-step RK4 is independent of the production implicit extent solver.
    let derivative =
        |x: f64| FORWARD_RATE * ((1.0 - x) * x - x.powi(3) / saha_molar(temperature(u, x)));
    let mut reference = initial_x;
    let dt = elapsed / 1000.0;
    for _ in 0..1000 {
        let a = derivative(reference);
        let b = derivative(reference + 0.5 * dt * a);
        let c = derivative(reference + 0.5 * dt * b);
        let d = derivative(reference + dt * c);
        reference += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let result = advance_adiabatic(
        &mechanism,
        &feed(initial_x),
        Volume::new(1.0),
        Temperature::new(initial_t),
        Time::new(elapsed),
        controls(),
        || false,
    )
    .unwrap();
    let x = amount(&result, "H+");
    relative(x, reference, 1e-3); // Backward Euler with bounded local error.
    relative(amount(&result, "e-"), x, 1e-12);
    relative(amount(&result, "H") + x, 1.0, 1e-12);
    relative(result.mass().value(), M_H * N_A, 1e-12);
    relative(result.temperature().value(), temperature(u, x), 1e-11);
    relative(
        result.pressure().value(),
        (1.0 + x) * R * result.temperature().value(),
        1e-12,
    );
    assert!(x > initial_x && result.temperature().value() < initial_t);
    assert!(result.balance().energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(
        result.balance().maximum_relative_charge_residual
            <= controls().relative_conservation_tolerance
    );
    assert!(result.balance().entropy_change_j_k > 0.0);
    assert_eq!(result.elapsed_time(), Time::new(elapsed));
    assert!(
        result.balance().accepted_steps + result.balance().rejected_steps
            <= controls().maximum_steps
    );
}

#[test]
fn g1_ionization_recombination_relaxes_to_adiabatic_saha_equilibrium() {
    let mechanism = hydrogen_ionization();
    let initial_x = 0.01;
    let initial_t = 16000.0;
    let u = energy(initial_t, initial_x);
    let (mut lo, mut hi) = (1000.0, initial_t);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        if energy(t, equilibrium_fraction(t)) < u {
            lo = t;
        } else {
            hi = t;
        }
    }
    let expected_t = 0.5 * (lo + hi);
    let expected_x = equilibrium_fraction(expected_t);
    let result = advance_adiabatic(
        &mechanism,
        &feed(initial_x),
        Volume::new(1.0),
        Temperature::new(initial_t),
        Time::new(2.0),
        controls(),
        || false,
    )
    .unwrap();
    relative(amount(&result, "H+"), expected_x, 2e-5);
    relative(result.temperature().value(), expected_t, 2e-5);
    relative(amount(&result, "e-"), amount(&result, "H+"), 1e-12);
    assert!(result.balance().energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(
        result.balance().maximum_relative_charge_residual
            <= controls().relative_conservation_tolerance
    );
    assert!(result.balance().entropy_change_j_k > 0.0);
}

#[test]
fn g0_ionization_charge_audit_rejects_trace_carrier_drift() {
    let mechanism = hydrogen_ionization();
    for scale in [1e-250, 1.0, 1e250] {
        let before = [scale, 1e-18 * scale, 1e-18 * scale];
        let after = [scale, 1e-18 * scale, 1.1e-18 * scale];
        let mut balance = IsochoricKineticsBalance::default();
        assert!(matches!(
            audit_charge(&mechanism, |i| before[i], |i| after[i], 1e-10, &mut balance),
            Err(IsochoricKineticsError::Balance(_))
        ));
        // A 1e-19 change relative to neutral gas is a 1/21 charge error.
        relative(balance.maximum_relative_charge_residual, 1.0 / 21.0, 1e-12);
    }
    let mut balance = IsochoricKineticsBalance::default();
    audit_charge(
        &mechanism,
        |i| [1.0, 0.0, 0.0][i],
        |i| [0.99, 0.01, 0.01][i],
        1e-10,
        &mut balance,
    )
    .unwrap();
    assert_eq!(balance.maximum_relative_charge_residual, 0.0);
}

#[test]
fn g1_ionization_kinetics_preserves_authored_net_charge_without_neutralizing() {
    let mechanism = hydrogen_ionization();
    let result = advance_adiabatic(
        &mechanism,
        &[
            (id("H"), Amount::new(0.98)),
            (id("H+"), Amount::new(0.02)),
            (id("e-"), Amount::new(0.01)),
        ],
        Volume::new(1.0),
        Temperature::new(16000.0),
        Time::new(0.01),
        controls(),
        || false,
    )
    .unwrap();
    relative(amount(&result, "H+") - amount(&result, "e-"), 0.01, 1e-11);
    assert!(amount(&result, "e-") > 0.01);
    assert!(
        result.balance().maximum_relative_charge_residual
            <= controls().relative_conservation_tolerance
    );
}
