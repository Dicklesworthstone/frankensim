use super::*;
use crate::transport::maxwell_stefan::{
    MaxwellStefanMolarSolveV1, MaxwellStefanSolveV1, MaxwellStefanV1,
};
use fs_qty::{AmountConcentration, MolarMass};

fn model(a: f64, b: f64) -> BinaryRegularSolutionV1 {
    BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        a,
        b,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}
fn eq(law: &BinaryRegularSolutionV1, t: f64, x: f64) -> BinaryPhaseEquilibriumV1 {
    law.equilibrate(
        Temperature::new(t),
        x,
        BinaryCoexistenceOptionsV1::default(),
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g0_ideal_pure_and_trace_limits_have_no_fictitious_solute() {
    let law = model(0.0, 0.0);
    for x in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let state = law.evaluate(Temperature::new(300.0), x).unwrap();
        let independent = if x == 0.0 || x == 1.0 {
            0.0
        } else {
            R * 300.0 * (x * x.ln() + (1.0 - x) * (-x).ln_1p())
        };
        near(state.mixing_gibbs_j_per_mol(), independent, 5e-15, 0.0);
        assert_eq!(state.thermodynamic_factor(), 1.0);
        assert_eq!(state.log_activity_coefficients(), [0.0, 0.0]);
        assert_eq!(state.driving_forces(0.2).unwrap(), [0.2, -0.2]);
        assert_eq!(eq(&law, 300.0, x).phases().len(), 1);
    }
    let pure = law.evaluate(Temperature::new(300.0), 0.0).unwrap();
    assert_eq!(pure.reduced_chemical_potentials(), [None, Some(0.0)]);
    let x = 1e-100;
    let trace = law.evaluate(Temperature::new(300.0), x).unwrap();
    near(
        trace.mixing_gibbs_j_per_mol(),
        R * 300.0 * x * (x.ln() - 1.0),
        5e-15,
        0.0,
    );
    near(
        trace.reduced_chemical_potentials()[1].unwrap(),
        -x,
        1e-15,
        0.0,
    );
    assert!(trace.reduced_chemical_potentials()[0].unwrap().is_finite());
}

#[test]
fn g1_free_energy_derivatives_and_calorics_share_one_temperature_law() {
    let law = model(-0.4, 900.0);
    let (t, x, dx, dt) = (400.0, 0.3, 1e-6, 1e-3);
    let state = law.evaluate(Temperature::new(t), x).unwrap();
    let g = |temperature, composition| {
        law.evaluate(Temperature::new(temperature), composition)
            .unwrap()
            .mixing_gibbs_j_per_mol()
    };
    let mu = state.reduced_chemical_potentials().map(Option::unwrap);
    near(
        (g(t, x + dx) - g(t, x - dx)) / (2.0 * dx),
        R * t * (mu[0] - mu[1]),
        2e-9,
        1e-7,
    );
    near(
        (g(t + dt, x) - g(t - dt, x)) / (2.0 * dt),
        -state.mixing_entropy_j_per_mol_k(),
        3e-10,
        1e-9,
    );
    near(
        state.mixing_gibbs_j_per_mol(),
        state.mixing_enthalpy_j_per_mol() - t * state.mixing_entropy_j_per_mol_k(),
        2e-15,
        1e-12,
    );
    near(
        x * mu[0] + (1.0 - x) * mu[1],
        state.mixing_gibbs_j_per_mol() / (R * t),
        2e-15,
        1e-15,
    );
    let chemical = |composition| {
        law.evaluate(Temperature::new(t), composition)
            .unwrap()
            .reduced_chemical_potentials()
            .map(Option::unwrap)
    };
    let a = chemical(x - dx);
    let b = chemical(x + dx);
    near(
        x * (b[0] - a[0]) / (2.0 * dx) + (1.0 - x) * (b[1] - a[1]) / (2.0 * dx),
        0.0,
        0.0,
        1e-10,
    );
}

#[test]
fn g1_common_tangent_binodal_spinodal_and_molar_lever_rule() {
    let law = model(3.0, 0.0);
    let result = eq(&law, 300.0, 0.3);
    assert_eq!(result.phases().len(), 2);
    let coexistence = result.coexistence().unwrap();
    let lo = coexistence.lower_a().mole_fraction_a();
    near(lo, 0.0707201816799448, 3e-12, 1e-15);
    near(
        coexistence.spinodal_mole_fractions()[0],
        0.21132486540518713,
        2e-15,
        0.0,
    );
    near(
        (lo / (1.0 - lo)).ln() + 3.0 * (1.0 - 2.0 * lo),
        0.0,
        0.0,
        1e-12,
    );
    let a = coexistence
        .lower_a()
        .reduced_chemical_potentials()
        .map(Option::unwrap);
    let b = coexistence
        .upper_a()
        .reduced_chemical_potentials()
        .map(Option::unwrap);
    for i in 0..2 {
        near(a[i], b[i], 0.0, 1e-12);
    }
    let phases = result.phases();
    near(
        phases.iter().map(|p| p.fraction_of_total_moles()).sum(),
        1.0,
        0.0,
        3e-16,
    );
    for (i, expected) in [0.3, 0.7].into_iter().enumerate() {
        near(
            phases
                .iter()
                .map(|p| p.fraction_of_total_moles() * p.mole_fractions()[i])
                .sum(),
            expected,
            2e-15,
            0.0,
        );
    }
    assert!(result.mixing_gibbs_j_per_mol() < result.feed().mixing_gibbs_j_per_mol());
    assert_eq!(
        eq(&law, 300.0, 0.01).feed_stability(),
        BinaryMixingStabilityV1::Stable
    );
    assert_eq!(
        eq(&law, 300.0, lo).feed_stability(),
        BinaryMixingStabilityV1::Binodal
    );
    assert_eq!(
        eq(&law, 300.0, 0.15).feed_stability(),
        BinaryMixingStabilityV1::Metastable
    );
    assert_eq!(
        eq(&law, 300.0, 0.5).feed_stability(),
        BinaryMixingStabilityV1::Spinodal
    );
}

#[test]
fn g3_temperature_crosses_upper_and_lower_critical_solution_regimes() {
    let upper = model(0.0, 900.0);
    assert_eq!(eq(&upper, 300.0, 0.5).phases().len(), 2);
    assert_eq!(
        eq(&upper, 450.0, 0.5).feed_stability(),
        BinaryMixingStabilityV1::Critical
    );
    assert_eq!(eq(&upper, 600.0, 0.5).phases().len(), 1);
    let lower = model(3.0, -300.0);
    assert_eq!(eq(&lower, 200.0, 0.5).phases().len(), 1);
    assert_eq!(eq(&lower, 600.0, 0.5).phases().len(), 2);
}

#[test]
fn g1_near_critical_root_does_not_collapse_to_central_stationarity() {
    for chi in [
        f64::from_bits(2.0_f64.to_bits() + 1),
        2.0 + 1e-10,
        2.0 + 1e-6,
    ] {
        let coexistence = model(chi, 0.0)
            .coexistence(
                Temperature::new(300.0),
                BinaryCoexistenceOptionsV1::default(),
                || false,
            )
            .unwrap()
            .unwrap();
        let m = 1.0 - 2.0 * coexistence.lower_a().mole_fraction_a();
        assert!(m > 0.0);
        // Independent Landau leading term: m² ~ 3(chi-2)/2.
        near(m * m, 1.5 * (chi - 2.0), 2e-6, 0.0);
        let bracket = coexistence.log_odds_bracket();
        assert!(bracket[0] > 0.0 && bracket[1] > bracket[0]);
    }
}

#[test]
fn g0_dilute_binodal_preserves_both_components_and_refuses_unrepresentable_gap() {
    let law = model(30.0, 0.0);
    let coexistence = law
        .coexistence(
            Temperature::new(300.0),
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap()
        .unwrap();
    let lo = coexistence.lower_a().mole_fractions();
    let hi = coexistence.upper_a().mole_fractions();
    assert_eq!(lo[0], hi[1]);
    assert_eq!(lo[1], hi[0]);
    near(lo[0], (-30.0_f64).exp(), 1e-10, 0.0);
    assert!(lo[0] > 0.0 && lo[0] < 1e-12);
    assert!(lo[1] < 1.0);
    assert!(matches!(
        model(100.0, 0.0).coexistence(
            Temperature::new(300.0),
            BinaryCoexistenceOptionsV1::default(),
            || false
        ),
        Err(BinarySolutionErrorV1::Unrepresentable(_))
    ));
    assert_eq!(eq(&model(100.0, 0.0), 300.0, 0.0).phases().len(), 1);
}

#[test]
fn g3_component_exchange_reflects_phases_and_preserves_energy() {
    let law = model(3.0, 0.0);
    let a = eq(&law, 300.0, 0.2);
    let b = eq(&law, 300.0, 0.8);
    near(
        a.mixing_gibbs_j_per_mol(),
        b.mixing_gibbs_j_per_mol(),
        2e-15,
        0.0,
    );
    near(
        a.phases()[0].fraction_of_total_moles(),
        b.phases()[1].fraction_of_total_moles(),
        2e-15,
        0.0,
    );
    assert_eq!(
        a.phases()[0].mole_fractions()[0],
        b.phases()[1].mole_fractions()[1]
    );
}

#[test]
fn g1_shared_law_drives_both_existing_maxwell_stefan_reference_frames() {
    let masses = [0.018, 0.046];
    let owner = MaxwellStefanV1::new(vec![
        (SpeciesId::new("A").unwrap(), MolarMass::new(masses[0])),
        (SpeciesId::new("B").unwrap(), MolarMass::new(masses[1])),
    ])
    .unwrap();
    let mut workspace = owner.workspace().unwrap();
    let (t, c, d, gradient) = (350.0, 18000.0, 8e-10, 0.2);
    for (chi, x) in [(0.0, 0.25), (1.5, 0.25), (3.0, 0.5), (1.5, 0.0)] {
        let state = model(chi, 0.0).evaluate(Temperature::new(t), x).unwrap();
        let fractions = state.mole_fractions();
        let forces = state.driving_forces(gradient).unwrap();
        assert_eq!(forces[0] + forces[1], 0.0);
        let diffusivity = [0.0, d, d, 0.0];
        let expected_molar = -c * d * state.thermodynamic_factor() * gradient;
        let molar = owner
            .evaluate_molar_driving_forces(
                &fractions,
                &forces,
                AmountConcentration::new(c),
                &diffusivity,
                MaxwellStefanMolarSolveV1::default(),
                &mut workspace,
                || false,
            )
            .unwrap();
        near(molar.molar_fluxes()[0], expected_molar, 3e-14, 0.0);
        near(molar.molar_fluxes()[1], -expected_molar, 3e-14, 0.0);
        if x > 0.0 {
            let mu_gradient = [
                R * t * (1.0 / x - 2.0 * chi * (1.0 - x)) * gradient,
                -R * t * (1.0 / (1.0 - x) - 2.0 * chi * x) * gradient,
            ];
            let entropy = -(molar.molar_fluxes()[0] * mu_gradient[0]
                + molar.molar_fluxes()[1] * mu_gradient[1])
                / t;
            assert!(entropy >= 0.0);
            near(
                molar.drag_entropy_production_w_m3_k().unwrap(),
                entropy,
                5e-14,
                0.0,
            );
        }
        let mass = owner
            .evaluate_driving_forces(
                &fractions,
                &forces,
                AmountConcentration::new(c),
                &diffusivity,
                MaxwellStefanSolveV1::default(),
                &mut workspace,
                || false,
            )
            .unwrap();
        let expected_mass =
            expected_molar * masses[0] * masses[1] / (x * masses[0] + (1.0 - x) * masses[1]);
        near(mass.mass_fluxes()[0], expected_mass, 4e-14, 0.0);
        near(mass.mass_fluxes()[1], -expected_mass, 4e-14, 0.0);
    }
}

#[test]
fn g4_invalid_input_work_bounds_and_cancellation_do_not_publish_partial_states() {
    let law = model(3.0, 0.0);
    for x in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
        assert!(law.evaluate(Temperature::new(300.0), x).is_err());
    }
    for t in [0.0, 99.0, 1001.0, f64::NAN] {
        assert_eq!(
            law.evaluate(Temperature::new(t), 0.5),
            Err(BinarySolutionErrorV1::OutsideTemperatureDomain)
        );
    }
    assert_eq!(
        law.equilibrate(
            Temperature::new(300.0),
            0.5,
            BinaryCoexistenceOptionsV1::default(),
            || true
        ),
        Err(BinarySolutionErrorV1::Cancelled)
    );
    let mut calls = 0;
    assert_eq!(
        law.equilibrate(
            Temperature::new(300.0),
            0.5,
            BinaryCoexistenceOptionsV1::default(),
            || {
                calls += 1;
                calls == 8
            }
        ),
        Err(BinarySolutionErrorV1::Cancelled)
    );
    assert_eq!(
        law.coexistence(
            Temperature::new(300.0),
            BinaryCoexistenceOptionsV1 {
                maximum_iterations: 1,
                ..Default::default()
            },
            || false
        ),
        Err(BinarySolutionErrorV1::IterationLimit)
    );
    assert!(
        law.equilibrate(
            Temperature::new(300.0),
            0.0,
            BinaryCoexistenceOptionsV1 {
                relative_tolerance: f64::NAN,
                ..Default::default()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(eq(&law, 300.0, 0.5), eq(&law, 300.0, 0.5));
}
