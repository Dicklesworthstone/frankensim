use super::*;
use crate::chemistry::isobaric::{self as constant_pressure, IsobaricReactor};

fn pressure_controls(method: IsochoricKineticsMethod) -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        method,
        initial_step_s: 0.002,
        maximum_step_s: 0.01,
        ..controls(2e-7)
    }
}

fn isothermal_dimer(recombine: bool, rate: f64) -> IdealGasKineticsV1 {
    mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, if recombine { 0.0 } else { 10.0 }),
            source("A2", 0.02, 8.0, 0.0, if recombine { 20.0 } else { 0.0 }),
        ],
        vec![1, 2],
        vec![if recombine {
            reaction("recombine", &[("A", 2)], &[("A2", 1)], rate).irreversible()
        } else {
            reaction("dissociate", &[("A2", 1)], &[("A", 2)], rate).irreversible()
        }],
    )
}

fn check_first_law(result: &constant_pressure::IsobaricKineticsResult, initial_u: f64) {
    let work =
        result.pressure().value() * (result.volume().value() - result.initial_volume().value());
    assert!((result.boundary_work().value() - work).abs() < 1e-12);
    assert!((result.internal_energy().value() - initial_u + work).abs() < 2e-8);
    assert!((result.enthalpy().value() - result.target_enthalpy().value()).abs() < 2e-8);
    assert!(result.balance().energy_residual_j.abs() <= 1e-8);
    assert!(result.balance().maximum_relative_element_residual < 1e-13);
    assert!(result.balance().entropy_change_j_k > 0.0);
    assert!(result.pressure_residual_pa().abs() < 1e-9);
    let moles: f64 = result.amounts().iter().map(|n| n.value()).sum();
    assert!(
        (moles / result.volume().value()
            - result.pressure().value() / (R * result.temperature().value()))
        .abs()
            < 1e-12
    );
}

#[test]
fn isothermal_dissociation_changes_volume_and_internal_energy_at_fixed_enthalpy() {
    let mechanism = isothermal_dimer(false, 2.0);
    for method in [
        IsochoricKineticsMethod::BackwardEuler,
        IsochoricKineticsMethod::Sdirk2,
    ] {
        let result = constant_pressure::advance_adiabatic(
            &mechanism,
            &[(id("A2"), Amount::new(1.0))],
            Pressure::new(R * 500.0),
            Temperature::new(500.0),
            Time::new(0.2),
            pressure_controls(method),
            || false,
        )
        .unwrap();
        let remaining = (-0.4_f64).exp();
        assert!((result.amounts()[1].value() - remaining).abs() < 3e-4);
        assert!((result.amounts()[0].value() - 2.0 * (1.0 - remaining)).abs() < 6e-4);
        assert!((result.temperature().value() - 500.0).abs() < 2e-9);
        assert!((result.volume().value() - (2.0 - result.amounts()[1].value())).abs() < 2e-13);
        assert!(result.volume().value() > result.initial_volume().value());
        check_first_law(&result, 7.0 * R * 500.0);
    }
}

#[test]
fn bimolecular_recombination_uses_actual_dilution_in_the_shared_stiff_controller() {
    let mechanism = isothermal_dimer(true, 2.0);
    let duration = 0.15;
    // y=n_A/n_A0: exact fixed-p law 1/y-1-ln(y)=4*k*x_A0*t.
    // Cp(A2)=2Cp(A) gives constant T, V/V0=(1+y)/2 and x_A0=1.
    let (mut lower, mut upper) = (0.01_f64, 1.0_f64);
    for _ in 0..80 {
        let middle = 0.5 * (lower + upper);
        if 1.0 / middle - 1.0 - middle.ln() > 8.0 * duration {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let expected = 0.5 * (lower + upper);
    for method in [
        IsochoricKineticsMethod::BackwardEuler,
        IsochoricKineticsMethod::Sdirk2,
    ] {
        let result = constant_pressure::advance_adiabatic(
            &mechanism,
            &[(id("A"), Amount::new(1.0))],
            Pressure::new(R * 500.0),
            Temperature::new(500.0),
            Time::new(duration),
            pressure_controls(method),
            || false,
        )
        .unwrap();
        let y = result.amounts()[0].value();
        assert!((y - expected).abs() < 3e-4, "{y} != {expected}");
        assert!(
            (y - 1.0 / (1.0 + 4.0 * duration)).abs() > 0.01,
            "constant-volume concentration law was substituted"
        );
        assert!((result.volume().value() - 0.5 * (1.0 + y)).abs() < 2e-13);
        assert!((result.temperature().value() - 500.0).abs() < 2e-9);
        assert!(result.boundary_work().value() < 0.0);
        check_first_law(&result, 3.0 * R * 500.0);
    }
}

fn exothermic_recombination(reference: f64, binding: f64, reactions: usize) -> IdealGasKineticsV1 {
    // Apply the exact elemental shift to original NASA a7 coefficients. A
    // multiply/divide by R would alter the represented 1024 K binding at
    // +2^40, turning a datum metamorphic into a physically different source.
    let coefficient_source = |name, mass, cp, formation_r, entropy| {
        let original = source(name, mass, cp, 0.0, entropy);
        let mut card = original.regions()[0].card().clone();
        card.parameters.get_mut("a7").unwrap().value = formation_r;
        Nasa9StandardStateModelV1::new(
            original.species().clone(),
            original.molar_mass(),
            original.convention().clone(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        )
        .unwrap()
    };
    mechanism(
        vec![
            coefficient_source("A", 0.01, 4.0, reference, 0.0),
            coefficient_source("A2", 0.02, 6.0, 2.0 * reference - binding, 20.0),
        ],
        vec![1, 2],
        (0..reactions)
            .map(|i| {
                reaction(
                    &format!("recombine-{i}"),
                    &[("A", 2)],
                    &[("A2", 1)],
                    20.0 / reactions as f64,
                )
                .irreversible()
            })
            .collect(),
    )
}

#[test]
fn reaction_enthalpy_preserves_elemental_references_and_real_binding_energy() {
    for binding in [512.0, 1024.0] {
        let mut baseline: Option<(f64, f64)> = None;
        for reference in [0.0, 2.0_f64.powi(40), -2.0_f64.powi(40)] {
            // Five reaction columns also exercise the original smaller ns+1
            // factorization; they share one physical stoichiometric direction.
            let mechanism = exothermic_recombination(reference, binding, 5);
            let coefficients: Vec<_> = mechanism
                .thermodynamics()
                .models()
                .iter()
                .map(|model| model.regions()[0].coefficients()[7])
                .collect();
            let retained_binding = 2.0 * coefficients[0] - coefficients[1];
            assert_eq!(retained_binding, binding);
            let result = constant_pressure::advance_adiabatic(
                &mechanism,
                &[(id("A"), Amount::new(1.0))],
                Pressure::new(R * 300.0),
                Temperature::new(300.0),
                Time::new(0.03),
                pressure_controls(IsochoricKineticsMethod::Sdirk2),
                || false,
            )
            .unwrap();
            let (a, a2) = (result.amounts()[0].value(), result.amounts()[1].value());
            let expected = (4.0 * 300.0 + retained_binding * a2) / (4.0 * a + 6.0 * a2);
            assert!((result.temperature().value() - expected).abs() < 3e-9);
            assert!((a + 2.0 * a2 - 1.0).abs() < 2e-14);
            assert!(result.temperature().value() > 300.0);
            assert!(result.balance().energy_residual_j.abs() <= 1e-8);
            if let Some((temperature, amount)) = baseline {
                assert!((result.temperature().value() - temperature).abs() < 1e-6);
                assert!((a2 - amount).abs() < 1e-9);
            } else {
                baseline = Some((result.temperature().value(), a2));
            }
        }
    }
}

#[test]
fn isobaric_constrained_rate_jacobian_matches_independent_finite_differences() {
    for mechanism in [
        isothermal_dimer(true, 2.0),
        exothermic_recombination(0.0, 512.0, 5),
    ] {
        let nr = mechanism.reactions().len();
        let mut reactor = solve::Reactor::new(
            mechanism,
            pressure_controls(IsochoricKineticsMethod::Sdirk2),
        )
        .unwrap();
        reactor.prepare_pressure().unwrap();
        reactor
            .initialize(
                &[Amount::new(1.0), Amount::new(0.0)],
                Volume::new(1.0),
                Some(Pressure::new(R * 500.0)),
                Temperature::new(500.0),
                None,
                &mut || false,
            )
            .unwrap();
        let extent = vec![0.08 / nr as f64; nr];
        let (_, analytic) = reactor
            .constrained_rate_linearization(&extent, 500.0)
            .unwrap();
        for column in 0..nr {
            let epsilon = 1e-6;
            let mut plus = extent.clone();
            plus[column] += epsilon;
            let mut minus = extent.clone();
            minus[column] -= epsilon;
            let (plus, _) = reactor
                .constrained_rate_linearization(&plus, 500.0)
                .unwrap();
            let (minus, _) = reactor
                .constrained_rate_linearization(&minus, 500.0)
                .unwrap();
            for row in 0..nr {
                let numeric = (plus[row] - minus[row]) / (2.0 * epsilon);
                let expected = analytic[row * nr + column];
                assert!(
                    (numeric - expected).abs() < 2e-7 * expected.abs().max(1.0),
                    "column {column}, row {row}: {numeric} != {expected}"
                );
            }
        }
    }
}

#[test]
fn isobaric_explicit_enthalpy_and_cancelled_prepared_calls_retain_atomic_history() {
    let mechanism = isothermal_dimer(false, 2.0);
    let controls = pressure_controls(IsochoricKineticsMethod::Sdirk2);
    let mut reactor = IsobaricReactor::new(mechanism, controls, || false).unwrap();
    let feed = [Amount::new(0.0), Amount::new(1.0)];
    let storage = reactor.workspace_bytes();
    let target = Energy::new(8.0 * R * 600.0);
    let zero = reactor
        .advance_adiabatic_enthalpy(
            &feed,
            target,
            Pressure::new(R * 500.0),
            Temperature::new(500.0),
            Time::new(0.0),
            || false,
        )
        .unwrap();
    assert!((zero.temperature().value() - 600.0).abs() < 2e-10);
    assert!((zero.initial_volume().value() - 1.2).abs() < 1e-13);
    assert_eq!(zero.volume(), zero.initial_volume());
    assert_eq!(zero.boundary_work().value(), 0.0);
    assert_eq!(zero.amounts(), &feed);
    let mut polls = 0;
    let refused = reactor.advance_adiabatic(
        &feed,
        Pressure::new(R * 500.0),
        Temperature::new(500.0),
        Time::new(0.2),
        || {
            polls += 1;
            polls > 100
        },
    );
    assert!(matches!(refused, Err(IsochoricKineticsError::Cancelled)));
    let replay = reactor
        .advance_adiabatic(
            &feed,
            Pressure::new(R * 500.0),
            Temperature::new(500.0),
            Time::new(0.2),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    let again = reactor
        .advance_adiabatic(
            &feed,
            Pressure::new(R * 500.0),
            Temperature::new(500.0),
            Time::new(0.2),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(replay, again);
    assert_eq!(reactor.workspace_bytes(), storage);
    assert_eq!(feed, [Amount::new(0.0), Amount::new(1.0)]);
    assert!(matches!(
        reactor.advance_adiabatic_enthalpy(
            &feed,
            Energy::new(1e20),
            Pressure::new(R * 500.0),
            Temperature::new(500.0),
            Time::new(0.0),
            || false
        ),
        Err(IsochoricKineticsError::EnergyNotBracketed { .. })
    ));
}
