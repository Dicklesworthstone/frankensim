//! Independent Troe equations on the existing material, kinetic and gas owners.
use super::*;

const TIMES: [f64; 3] = [200.0, 2000.0, 2500.0];
const ALPHA_SLOPE: f64 = 0.0001;

fn troe_card(
    alpha: f64,
    times: [f64; 3],
    include_t2: bool,
    low: f64,
    efficiencies: [f64; 3],
    upper: f64,
) -> MaterialCard {
    let original = collider_card(200.0, low, efficiencies);
    let mut claims = original.claims().clone();
    let mut values = vec![
        (TROE_ALPHA_PROPERTY, Dims([0; 6]), alpha),
        (TROE_T3_PROPERTY, Temperature::DIMS, times[0]),
        (TROE_T1_PROPERTY, Temperature::DIMS, times[1]),
    ];
    if include_t2 {
        values.push((TROE_T2_PROPERTY, Temperature::DIMS, times[2]));
    }
    for (name, dims, value) in values {
        let (value, interpolation) = if name == TROE_ALPHA_PROPERTY {
            (
                PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (200.0, value - 400.0 * ALPHA_SLOPE),
                        (upper, value + (upper - 600.0) * ALPHA_SLOPE),
                    ],
                    dims,
                },
                InterpolationPolicy::LinearInside,
            )
        } else {
            (
                PropertyValue::Scalar { value, dims },
                InterpolationPolicy::ConstantWithinValidity,
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        upper,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn bind(
    card: MaterialCard,
    regime: ColliderRegime,
    activated: bool,
) -> Result<MaterialRateScale, GasReactionError> {
    let base = if activated {
        MaterialRateScale::modified_arrhenius(
            card,
            MaterialPropertySelection::SingleClaimOnly,
            Pressure::new(100000.0),
            arrhenius::ArrheniusRateBasis::PressureActivity,
        )
    } else {
        MaterialRateScale::new(
            card,
            MaterialPropertySelection::SingleClaimOnly,
            Pressure::new(100000.0),
        )
    }?;
    base.with_colliders(bindings(regime), MaterialPropertySelection::SingleClaimOnly)
}

fn sample(
    law: &GasReaction,
    gas: &ClosedGasParcel,
) -> DeclaredStoichiometricMassActionEvaluationV1 {
    let amounts = law.inventory(gas).unwrap();
    law.evaluate_channel(
        &law.kinetics,
        law.rate_material.as_ref(),
        gas,
        &amounts,
        &mut || false,
    )
    .unwrap()
    .0
}

fn center(temperature: f64, alpha: f64, four: bool) -> f64 {
    (1.0 - alpha) * (-temperature / TIMES[0]).exp()
        + alpha * (-temperature / TIMES[1]).exp()
        + if four {
            (-TIMES[2] / temperature).exp()
        } else {
            0.0
        }
}

// Literal independent equations from the Cantera reference, using direct Pr
// and base-10 logarithms. Production instead uses natural logs and scaled
// squared coordinates, including where this direct center would underflow.
fn hand_factor(reduced: f64, log_center10: f64) -> f64 {
    let c = -0.4 - 0.67 * log_center10;
    let n = 0.75 - 1.27 * log_center10;
    let x = reduced.log10() + c;
    let y = n - 0.14 * x;
    let broadening = if y == 0.0 {
        1.0
    } else {
        10.0_f64.powf(log_center10 / (1.0 + (x / y).powi(2)))
    };
    reduced / (1.0 + reduced) * broadening
}

#[test]
fn troe_forms_match_independent_pressure_equations_and_continuous_singular_limit() {
    for (regime, four) in [
        (ColliderRegime::Troe3, false),
        (ColliderRegime::Troe4, true),
    ] {
        for low in [2e-100, 20.0, 2e100] {
            let law = reaction(4000.0, 3.5, 1.0)
                .with_material_rate(
                    bind(
                        troe_card(0.35, TIMES, four, low, [1.0, 2.0, 4.0], 1800.0),
                        regime,
                        false,
                    )
                    .unwrap(),
                )
                .unwrap();
            for pressure in [50000.0, 100000.0, 200000.0] {
                let gas = parcel(&law, pressure, [0.25, 0.25, 0.5]);
                let rate = sample(&law, &gas);
                let q = 2.75 * pressure / (R * 600.0 * C0);
                let reduced = low * q / 200.0;
                let expected = hand_factor(reduced, center(600.0, 0.35, four).log10());
                near(rate.kinetic_multiplier().value() / expected, 1.0, 2e-12);
                // Common broadening preserves the original equilibrium ratio.
                let equilibrium = rate
                    .reverse_rate()
                    .equilibrium()
                    .equilibrium_constant()
                    .unwrap()
                    .value();
                near(
                    rate.forward_progress_rate().unwrap().value()
                        / rate.reverse_progress_rate().unwrap().value(),
                    equilibrium,
                    1e-10,
                );
                assert!(rate.kinetic_multiplier().value() > 0.0);
                assert!(rate.kinetic_multiplier().value() <= 1.0);
            }
        }
        // Choose Pr from the independent y=n-.14*(log10Pr+c)=0 equation.
        // The continuous broadening value is exactly one, not a refusal.
        let lc = center(600.0, 0.35, four).log10();
        let c = -0.4 - 0.67 * lc;
        let n = 0.75 - 1.27 * lc;
        let reduced = 10.0_f64.powf(n / 0.14 - c);
        let q = 1.75 * 100000.0 / (R * 600.0 * C0);
        let law = reaction(0.0, 3.5, 1.0)
            .with_material_rate(
                bind(
                    troe_card(
                        0.35,
                        TIMES,
                        four,
                        200.0 * reduced / q,
                        [1.0, 2.0, 4.0],
                        1800.0,
                    ),
                    regime,
                    false,
                )
                .unwrap(),
            )
            .unwrap();
        let gas = parcel(&law, 100000.0, [0.75, 0.0, 0.25]);
        near(
            sample(&law, &gas).kinetic_multiplier().value(),
            reduced / (1.0 + reduced),
            1e-12,
        );
    }
    // exp(-1000) is not a represented center, but its logarithm and the
    // resulting high-pressure broadened rate remain positive and finite.
    let law = reaction(0.0, 3.5, 1.0)
        .with_material_rate(
            bind(
                troe_card(
                    0.35,
                    [0.6, 0.6, TIMES[2]],
                    false,
                    1e300,
                    [1.0, 2.0, 4.0],
                    1800.0,
                ),
                ColliderRegime::Troe3,
                false,
            )
            .unwrap(),
        )
        .unwrap();
    let gas = parcel(&law, 100000.0, [0.75, 0.0, 0.25]);
    let reduced = 1e300 * (1.75 * 100000.0 / (R * 600.0 * C0)) / 200.0;
    let expected = hand_factor(reduced, -1000.0 / core::f64::consts::LN_10);
    let rate = sample(&law, &gas);
    assert!(rate.forward_progress_rate().unwrap().value() > 0.0);
    near(rate.kinetic_multiplier().value() / expected, 1.0, 2e-10);
}

#[test]
fn source_only_troe_broadening_reaches_implicit_piston_motion_and_conserved_heat() {
    let run = |alpha| {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(
                bind(
                    troe_card(alpha, TIMES, true, 2000.0, [1.0, 2.0, 4.0], 1800.0),
                    ColliderRegime::Troe4,
                    true,
                )
                .unwrap(),
            )
            .unwrap()
            .with_implicit(ImplicitReactionSolve {
                extent_bracket_mol: [-1.0, 1.0],
                extent_tolerance_mol: 1e-12,
                maximum_iterations: 128,
            })
            .unwrap();
        let mut piston = MaterialGasPiston::new(
            parcel(&law, 100000.0, [0.75, 0.0, 0.25]),
            body(8000.0),
            0.01,
            0.0,
        )
        .unwrap();
        let initial = piston.clone();
        let energy0 =
            initial.gas().mass_kg() * initial.gas().state().specific_internal_energy_j_kg();
        let mut pins = None;
        for _ in 0..12 {
            let volume = piston.gas().volume_m3();
            let candidate = piston
                .reactive_candidate(
                    &law,
                    drive(),
                    chemical_controls(),
                    mechanical_controls(),
                    || false,
                )
                .unwrap();
            let report = candidate.step();
            let reaction = &report.reaction;
            let resolved = reaction.rate_material.as_ref().unwrap();
            let temperature = resolved.query().axes()["temperature"];
            let a = alpha + ALPHA_SLOPE * (temperature - 600.0);
            near(
                resolved.property(TROE_ALPHA_PROPERTY).unwrap().value_si(),
                a,
                1e-12,
            );
            let n: Vec<_> = reaction.species.iter().map(|s| s.final_mol).collect();
            let q = (n[0] + 2.0 * temperature / 600.0 * n[1] + 4.0 * n[2]) / (volume * C0);
            let high = 2000.0 * (-5000.0 / (R * temperature)).exp();
            near(
                reaction.rate.kinetic_multiplier().value(),
                hand_factor(2000.0 * q / high, center(temperature, a, true).log10()),
                1e-10,
            );
            assert!(temperature > 600.0);
            assert!(reaction.implicit.as_ref().unwrap().residual_mol.abs() <= 1e-12);
            assert!(
                report.total_energy_residual_j.abs() <= mechanical_controls().energy_tolerance_j
            );
            let selected: Vec<_> = resolved
                .properties()
                .iter()
                .map(|p| p.answer().receipt.selected.clone())
                .collect();
            if let Some(previous) = &pins {
                assert_eq!(&selected, previous);
            } else {
                pins = Some(selected);
            }
            piston = candidate.into_parts().0;
        }
        let final_energy = piston.gas().mass_kg()
            * piston.gas().state().specific_internal_energy_j_kg()
            + 0.5 * piston.body().inertia().mass() * piston.piston().velocity_m_s().powi(2)
            + drive().exterior_pressure_pa * (piston.gas().volume_m3() - initial.gas().volume_m3());
        near(
            final_energy,
            energy0,
            12.0 * mechanical_controls().energy_tolerance_j,
        );
        assert_eq!(piston.gas().mass_kg(), initial.gas().mass_kg());
        assert_eq!(piston.body(), initial.body());
        piston
    };
    let narrow = run(0.7);
    let broad = run(0.2);
    assert!(narrow.piston().velocity_m_s() > broad.piston().velocity_m_s());
    assert!(broad.piston().velocity_m_s() > 0.0);
    assert!(narrow.gas().state().state().temperature > broad.gas().state().state().temperature);
}

#[test]
fn troe_source_domains_exact_absence_and_numerical_refusals_preserve_accepted_state() {
    assert!(
        bind(
            troe_card(0.35, TIMES, false, 20.0, [1.0; 3], 1800.0),
            ColliderRegime::Troe4,
            false
        )
        .is_err()
    );
    for (regime, four) in [
        (ColliderRegime::Troe3, false),
        (ColliderRegime::Troe4, true),
    ] {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(
                bind(
                    troe_card(0.35, TIMES, four, 20.0, [0.0; 3], 1800.0),
                    regime,
                    false,
                )
                .unwrap(),
            )
            .unwrap();
        let before = parcel(&law, 100000.0, [0.25, 0.25, 0.5]);
        let stopped = law
            .parcel_candidate(&before, drive().duration_s, chemical_controls(), || false)
            .unwrap();
        assert_eq!(stopped.gas(), &before);
        assert_eq!(stopped.step().stored_energy_change_j, 0.0);
        assert_eq!(stopped.step().rate.kinetic_multiplier().value(), 0.0);
    }
    for (alpha, times, regime, efficiencies) in [
        (1.1, TIMES, ColliderRegime::Troe3, [0.0; 3]),
        (0.35, [0.0, 2000.0, 2500.0], ColliderRegime::Troe3, [1.0; 3]),
        (
            0.35,
            [1e-320, 2000.0, 2500.0],
            ColliderRegime::Troe3,
            [1.0; 3],
        ),
        (0.35, [0.6, 0.6, 2500.0], ColliderRegime::Troe3, [1.0; 3]),
        (
            0.35,
            [10000.0, 10000.0, 1.0],
            ColliderRegime::Troe4,
            [1.0; 3],
        ),
    ] {
        let law = reaction(4000.0, 3.5, 1.0)
            .with_material_rate(
                bind(
                    troe_card(
                        alpha,
                        times,
                        regime == ColliderRegime::Troe4,
                        20.0,
                        efficiencies,
                        1800.0,
                    ),
                    regime,
                    false,
                )
                .unwrap(),
            )
            .unwrap();
        let before = parcel(&law, 100000.0, [0.75, 0.0, 0.25]);
        let original = before.clone();
        assert!(
            law.parcel_candidate(&before, drive().duration_s, chemical_controls(), || false)
                .is_err()
        );
        assert_eq!(before, original);
    }
    // A pinned Troe source must not freeze when reaction heating leaves its
    // support, even though every existing Arrhenius/collider claim still fits.
    let law = reaction(4000.0, 3.5, 1.0)
        .with_material_rate(
            bind(
                troe_card(0.35, TIMES, true, 2000.0, [1.0, 2.0, 4.0], 600.000001),
                ColliderRegime::Troe4,
                true,
            )
            .unwrap(),
        )
        .unwrap();
    let before = parcel(&law, 100000.0, [0.75, 0.0, 0.25]);
    let original = before.clone();
    let mut polls = 0;
    let heated = law
        .parcel_candidate(&before, drive().duration_s, chemical_controls(), || {
            polls += 1;
            false
        })
        .unwrap()
        .into_parts()
        .0;
    assert!(heated.state().state().temperature > 600.000001);
    let accepted = heated.clone();
    assert!(matches!(
        law.parcel_candidate(&heated, drive().duration_s, chemical_controls(), || false),
        Err(GasReactionError::RateMaterial(_))
    ));
    assert_eq!(heated, accepted);
    let mut seen = 0;
    assert!(matches!(
        law.parcel_candidate(&before, drive().duration_s, chemical_controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(GasReactionError::Cancelled)
    ));
    assert_eq!(before, original);
}
