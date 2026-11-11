//! G3: the derivative follows the admitted piecewise enthalpy chart exactly.
use fs_blake3::ContentHash;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, PhaseStateError, SolidLiquidPhase,
};

fn knot(h: f64, temperature: f64, liquid: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: temperature,
        liquid_mass_fraction: liquid,
        bulk_density_kg_m3: 1000.0,
    }
}

fn melting_curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x93; 32]),
        vec![
            knot(-100.0, 200.0, 0.0),
            knot(100.0, 300.0, 0.0),
            knot(500.0, 300.0, 1.0),
            knot(1300.0, 500.0, 1.0),
        ],
    )
    .unwrap()
}

fn slope(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> f64 {
    curve
        .temperature_derivative_at_specific_enthalpy(h)
        .unwrap()
}

#[test]
fn solid_latent_and_liquid_slopes_match_both_directional_differences() {
    let curve = melting_curve();
    for (enthalpy, expected) in [(0.0, 0.5), (300.0, 0.0), (900.0, 0.25)] {
        let derivative = slope(&curve, enthalpy);
        assert_eq!(derivative, expected);
        let temperature = curve
            .state_at_specific_enthalpy(enthalpy)
            .unwrap()
            .temperature_k();
        for increment in [-8.0, 8.0] {
            let shifted = curve
                .state_at_specific_enthalpy(enthalpy + increment)
                .unwrap()
                .temperature_k();
            assert_eq!(derivative, (shifted - temperature) / increment);
        }
    }
    for enthalpy in [101.0, 299.0, 499.0] {
        assert_eq!(slope(&curve, enthalpy).to_bits(), 0.0f64.to_bits());
        assert_eq!(
            curve
                .state_at_specific_enthalpy(enthalpy)
                .unwrap()
                .temperature_k(),
            300.0
        );
    }
}

#[test]
fn internal_knots_use_the_right_slope_and_upper_endpoint_uses_the_left_slope() {
    let curve = melting_curve();
    for (enthalpy, expected) in [
        (-100.0, 0.5),
        (100.0f64.next_down(), 0.5),
        (100.0, 0.0),
        (500.0f64.next_down(), 0.0),
        (500.0, 0.25),
        (1300.0, 0.25),
    ] {
        assert_eq!(slope(&curve, enthalpy), expected);
    }
    for (phase, liquid) in [
        (SolidLiquidPhase::Solid, 0.0),
        (SolidLiquidPhase::Liquid, 1.0),
    ] {
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x94; 32]),
            phase,
            vec![
                knot(-20.0, 200.0, liquid),
                knot(40.0, 230.0, liquid),
                knot(100.0, 320.0, liquid),
            ],
        )
        .unwrap();
        for (enthalpy, expected) in [
            (-20.0, 0.5),
            (10.0, 0.5),
            (40.0, 1.5),
            (70.0, 1.5),
            (100.0, 1.5),
        ] {
            assert_eq!(slope(&curve, enthalpy), expected);
        }
    }
}

#[test]
fn nonfinite_and_outside_queries_return_the_existing_domain_refusals() {
    let curve = melting_curve();
    for enthalpy in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            curve.temperature_derivative_at_specific_enthalpy(enthalpy),
            Err(PhaseStateError::NonFiniteSpecificEnthalpy)
        );
    }
    for enthalpy in [-101.0, 1301.0] {
        assert_eq!(
            curve.temperature_derivative_at_specific_enthalpy(enthalpy),
            Err(PhaseStateError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg: enthalpy,
                lower_j_kg: -100.0,
                upper_j_kg: 1300.0,
            })
        );
    }
}

#[test]
fn unrepresentable_positive_slopes_refuse_without_erasing_plateaus_or_subnormals() {
    let smallest = f64::from_bits(1);
    let single = |h1, t0, t1| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x95; 32]),
            SolidLiquidPhase::Solid,
            vec![knot(0.0, t0, 0.0), knot(h1, t1, 0.0)],
        )
        .unwrap()
    };
    let overflow = single(smallest, 1.0, 2.0);
    let underflow = single(f64::MAX, smallest, f64::from_bits(2));
    for (curve, upper) in [(&overflow, smallest), (&underflow, f64::MAX)] {
        for enthalpy in [0.0, upper] {
            assert_eq!(
                curve.temperature_derivative_at_specific_enthalpy(enthalpy),
                Err(PhaseStateError::UnrepresentableTemperatureDerivative)
            );
        }
    }
    let subnormal = single(1.0, smallest, f64::from_bits(2));
    assert_eq!(slope(&subnormal, 0.5).to_bits(), smallest.to_bits());
    let plateau = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x96; 32]),
        vec![knot(0.0, 1.0, 0.0), knot(smallest, 1.0, 1.0)],
    )
    .unwrap();
    assert_eq!(slope(&plateau, 0.0).to_bits(), 0.0f64.to_bits());
}

#[test]
fn signed_zero_queries_resolve_the_same_boundary_state_and_right_slope() {
    for lower in [-0.0, 0.0] {
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x97; 32]),
            vec![knot(lower, 300.0, 0.0), knot(100.0, 350.0, 1.0)],
        )
        .unwrap();
        let negative = curve.state_at_specific_enthalpy(-0.0).unwrap();
        let positive = curve.state_at_specific_enthalpy(0.0).unwrap();
        assert_eq!(negative.temperature_k(), positive.temperature_k());
        assert_eq!(negative.temperature_k(), 300.0);
        assert_eq!(
            negative.solid_mass_fraction(),
            positive.solid_mass_fraction()
        );
        assert_eq!(
            negative.liquid_mass_fraction(),
            positive.liquid_mass_fraction()
        );
        assert_eq!(negative.bulk_density_kg_m3(), positive.bulk_density_kg_m3());
        assert_eq!(negative.phase(), positive.phase());
        assert_eq!(slope(&curve, -0.0), 0.5);
        assert_eq!(slope(&curve, 0.0), 0.5);
    }
}

/// The same chart as a Component of the composition contract: a body of
/// fixed mass whose heat port returns the chart's entropy-conjugate mean.
#[cfg(feature = "component")]
mod body {
    use super::{ContentHash, EquilibriumEnthalpyPhaseCurve, PhaseStateError, knot, melting_curve};
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_material::phase::EquilibriumPhaseState;
    use fs_material::phase::component::UniformEnthalpyBody;
    use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
    use fs_ports::component::{Class, Component, PortTangent, Proposal};

    /// An eighth of a kilogram over half a second: a watt moves the specific
    /// enthalpy by 4 J/kg, exactly.
    const MASS_KG: f64 = 0.125;
    const DT_S: f64 = 0.5;

    fn in_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let key = StreamKey { seed: 31, kernel_id: 3, tile: 0, iteration: 0 };
            f(&Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic))
        })
    }

    /// Equal to the bit.
    #[track_caller]
    fn same(actual: f64, expected: f64) {
        assert_eq!(actual.to_bits(), expected.to_bits(), "{actual} vs {expected}");
    }

    /// One step from enthalpy `from` under `heat_w`: the state reached and
    /// the temperature returned.
    fn stepped(
        body: &UniformEnthalpyBody,
        from: f64,
        heat_w: f64,
        cx: &Cx<'_>,
    ) -> (EquilibriumPhaseState, f64) {
        let old = body.curve().state_at_specific_enthalpy(from).unwrap();
        let mut next = old;
        let mut returned = [f64::NAN];
        let proposal = body.propose(&old, DT_S, &[heat_w], &mut next, &mut returned, cx).unwrap();
        assert!(matches!(proposal, Proposal::Step(_)));
        (next, returned[0])
    }

    fn slope_k_w(body: &UniformEnthalpyBody, from: f64, heat_w: f64, cx: &Cx<'_>) -> f64 {
        let old = body.curve().state_at_specific_enthalpy(from).unwrap();
        let (mut forward, mut back) = ([f64::NAN], [f64::NAN]);
        body.jvp(&old, DT_S, &[heat_w], &[1.0], &mut forward, cx).unwrap();
        body.vjp(&old, DT_S, &[heat_w], &[1.0], &mut back, cx).unwrap();
        assert_eq!(forward[0].to_bits(), back[0].to_bits());
        forward[0]
    }

    #[test]
    fn a_uniform_body_conforms_and_returns_the_charts_mean_on_and_across_its_corners() {
        // Solid from 200 K to 300 K over 200 J/kg, a plateau at 300 K from
        // 100 to 500 J/kg, liquid to 500 K at 1300 J/kg.
        let curve = melting_curve();
        let body = UniformEnthalpyBody::new(curve.clone(), MASS_KG).unwrap();
        let ports: Vec<&str> = body.ports().iter().map(|port| port.id().as_str()).collect();
        assert_eq!(ports, ["heat"]);
        assert_eq!(body.class(), Class::Monotone1d);
        let state = |h: f64| curve.state_at_specific_enthalpy(h).unwrap();
        let stale = state(900.0);
        in_cx(&CancelGate::new_clock_free(), |cx| {
            for (from, heat_w, what) in [
                (0.0, 12.5, "within the solid"),
                (0.0, 25.0, "ending on the lower corner of the plateau"),
                (100.0, 25.0, "starting on that corner"),
                (0.0, 50.0, "from the solid onto the plateau"),
                (200.0, 50.0, "along the plateau"),
                (300.0, 50.0, "ending on its upper corner"),
                (500.0, 50.0, "starting on that corner"),
                (50.0, 200.0, "across the whole plateau into the liquid"),
                (900.0, -225.0, "back down across both corners"),
                (900.0, -100.0, "from the liquid onto the upper corner"),
            ] {
                let old = state(from);
                let reached = from + DT_S * heat_w / MASS_KG;
                same(reached, from + 4.0 * heat_w);
                let case = Case { state: &old, stale: &stale, dt_s: DT_S, inputs: &[heat_w] };
                let Ok(Observed::Stepped { receipt, outputs, .. }) = check_component(&body, &case, cx)
                else {
                    panic!("no conforming step {what}");
                };
                // The absolute allowance is for the central difference, not
                // for the slope: where a step ends on a corner the slope is
                // continuous but its own derivative jumps, and a difference
                // across the corner is off by a quarter of that jump times
                // its step, about 2e-6 K/W here.
                check_port_tangent(&body, &case, 2.0e-5, 1.0e-5, cx)
                    .unwrap_or_else(|failure| panic!("{what}: {failure}"));
                // The state and the temperature are the chart's own answers
                // to the same two enthalpies, to the bit.
                let mean = curve.thermodynamic_mean_temperature(from, reached, 3, || false).unwrap();
                let (next, returned) = stepped(&body, from, heat_w, cx);
                assert_eq!(next, state(reached), "{what}");
                assert_eq!(returned.to_bits(), mean.temperature_k.to_bits(), "{what}");
                assert_eq!(outputs[0].to_bits(), returned.to_bits(), "{what}");
                assert!(old.temperature_k().min(next.temperature_k()) <= returned, "{what}");
                assert!(returned <= old.temperature_k().max(next.temperature_k()), "{what}");
                // The heat over the returned temperature is the chart's
                // entropy change: the body generates none of its own.
                let entropy = MASS_KG * curve.specific_entropy_change_j_kg_k(from, reached).unwrap();
                let carried = DT_S * heat_w / returned;
                assert!((carried - entropy).abs() <= 1.0e-13 * entropy.abs(), "{what}: {carried} vs {entropy}");
                assert_eq!(receipt.entropy_generated_j_per_k.to_bits(), 0.0_f64.to_bits());
                // All of the heat is stored; nothing rounds at these values.
                same(receipt.stored_change_j, DT_S * heat_w);
                same(receipt.numerical_loss_j, 0.0);
                same(body.stored_energy_j(&next).unwrap(), MASS_KG * reached);
                // The slope is the chart's derivative of its mean in the
                // final enthalpy, per watt.
                let slope = slope_k_w(&body, from, heat_w, cx);
                let expected = mean.d_temperature_d_final_enthalpy_kg_k_j * (DT_S / MASS_KG);
                assert_eq!(slope.to_bits(), expected.to_bits(), "{what}");
            }
            // Along the plateau the mean is the transition temperature and
            // does not move with the heat.
            assert_eq!(stepped(&body, 200.0, 50.0, cx).1.to_bits(), 300.0_f64.to_bits());
            assert_eq!(slope_k_w(&body, 200.0, 50.0, cx).to_bits(), 0.0_f64.to_bits());

            // It rises with the heat: from one state in the solid, heats
            // that reach from near the cold end of the solid to near the
            // hot end of the liquid.
            let heats = [-20.0, -5.0, 0.0, 5.0, 25.0, 50.0, 125.0, 200.0, 320.0];
            let returned: Vec<f64> = heats.iter().map(|&heat_w| stepped(&body, 0.0, heat_w, cx).1).collect();
            assert!(returned.windows(2).all(|pair| pair[0] < pair[1]), "{returned:?}");

            // No heat at all returns the temperature of the state. On a
            // corner of the chart the slope there is the chart's one-sided
            // choice: half its right derivative. At the upper corner of the
            // plateau that is the liquid's, 0.5 * 0.25 K kg/J * 4 J/kg/W,
            // and a little heat confirms it from the right; from the left,
            // on the plateau, the temperature does not move at all.
            for (at, kelvin, slope) in [(0.0, 250.0, 1.0), (300.0, 300.0, 0.0), (100.0, 300.0, 0.0), (500.0, 300.0, 0.5)] {
                assert_eq!(stepped(&body, at, 0.0, cx), (state(at), kelvin));
                same(slope_k_w(&body, at, 0.0, cx), slope);
            }
            let little = 1.0e-3;
            let right = (stepped(&body, 500.0, little, cx).1 - 300.0) / little;
            let left = (300.0 - stepped(&body, 500.0, -little, cx).1) / little;
            assert!((right - 0.5).abs() <= 1.0e-6, "{right}");
            same(left, 0.0);

            // Heat too small to move the enthalpy is not refused: the state
            // stays, and the receipt books that heat as lost to rounding.
            let old = state(900.0);
            let faint = 1.0e-30;
            let case = Case { state: &old, stale: &stale, dt_s: DT_S, inputs: &[faint] };
            let Ok(Observed::Stepped { receipt, outputs, .. }) = check_component(&body, &case, cx)
            else {
                panic!("no conforming step under a faint heat");
            };
            same(outputs[0], 400.0);
            same(receipt.stored_change_j, 0.0);
            same(receipt.numerical_loss_j, DT_S * faint);
            assert_eq!(stepped(&body, 900.0, faint, cx).0, old);
            // A state that does not move is kept as it is, sign of zero
            // and identity included, as the chart's own advance keeps it.
            let (kept, _) = stepped(&body, -0.0, 0.0, cx);
            assert_eq!(kept, state(-0.0));
            same(kept.specific_enthalpy_j_kg(), -0.0);
        });
    }

    #[test]
    fn a_uniform_body_refuses_what_its_chart_refuses_and_leaves_the_candidate_alone() {
        let curve = melting_curve();
        let refused = |result: Result<UniformEnthalpyBody, PhaseStateError>| {
            matches!(result, Err(PhaseStateError::InvalidIsobaricHeatStep { .. }))
        };
        for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(refused(UniformEnthalpyBody::new(curve.clone(), mass)), "{mass}");
        }
        let body = UniformEnthalpyBody::new(curve.clone(), MASS_KG).unwrap();
        let old = curve.state_at_specific_enthalpy(0.0).unwrap();
        let stale = curve.state_at_specific_enthalpy(900.0).unwrap();
        // The same knots under another card are another chart.
        let other = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x94; 32]),
            vec![
                knot(-100.0, 200.0, 0.0),
                knot(100.0, 300.0, 0.0),
                knot(500.0, 300.0, 1.0),
                knot(1300.0, 500.0, 1.0),
            ],
        )
        .unwrap();
        let foreign = other.state_at_specific_enthalpy(0.0).unwrap();
        let gate = CancelGate::new_clock_free();
        in_cx(&gate, |cx| {
            let mut next = stale;
            let (mut one, mut two) = ([7.0], [7.0, 8.0]);
            // Two outputs, two seeds, two slots.
            assert!(body.propose(&old, DT_S, &[1.0], &mut next, &mut two, cx).is_err());
            assert!(body.jvp(&old, DT_S, &[1.0], &[1.0, 0.0], &mut one, cx).is_err());
            assert!(body.vjp(&old, DT_S, &[1.0], &[1.0], &mut two, cx).is_err());
            assert_eq!(two.map(f64::to_bits), [7.0_f64.to_bits(), 8.0_f64.to_bits()]);
            let mut attempt = |state: &EquilibriumPhaseState, dt_s: f64, inputs: &[f64]| {
                let error = body.propose(state, dt_s, inputs, &mut next, &mut one, cx).unwrap_err();
                assert!(body.jvp(state, dt_s, inputs, &[1.0], &mut one, cx).is_err());
                assert_eq!(next, stale);
                assert_eq!(one[0].to_bits(), 7.0_f64.to_bits());
                error
            };
            let step = |what| PhaseStateError::InvalidIsobaricHeatStep { what };
            assert_eq!(attempt(&old, DT_S, &[]), step("a uniform enthalpy body requires one heat rate"));
            assert_eq!(attempt(&old, DT_S, &[1.0, 2.0]), step("a uniform enthalpy body requires one heat rate"));
            assert_eq!(attempt(&old, DT_S, &[f64::NAN]), PhaseStateError::NonFiniteSpecificEnergyIncrement);
            for dt_s in [0.0, -1.0, f64::NAN] {
                assert_eq!(
                    attempt(&old, dt_s, &[1.0]),
                    step("a uniform enthalpy body requires a positive finite timestep")
                );
            }
            assert_eq!(attempt(&foreign, DT_S, &[1.0]), PhaseStateError::StateCurveMismatch);
            assert_eq!(body.stored_energy_j(&foreign), Err(PhaseStateError::StateCurveMismatch));
            // A mass so large that its enthalpy is not a number of joules:
            // the stored energy and the step both say so.
            let vast = UniformEnthalpyBody::new(curve.clone(), 1.0e306).unwrap();
            assert_eq!(
                vast.stored_energy_j(&stale),
                Err(step("a uniform enthalpy body's stored energy is unrepresentable"))
            );
            let (mut kept, mut slot) = (old, [7.0]);
            assert_eq!(
                vast.propose(&stale, DT_S, &[1.0], &mut kept, &mut slot, cx),
                Err(step("a uniform enthalpy body's energy account is unrepresentable"))
            );
            assert_eq!((kept, slot[0].to_bits()), (old, 7.0_f64.to_bits()));
            // Heat that would take the body off either end of its chart is
            // the chart's own refusal.
            for heat_w in [330.0, -30.0] {
                assert!(
                    matches!(attempt(&old, DT_S, &[heat_w]), PhaseStateError::OutsideEnthalpyDomain { .. }),
                    "{heat_w}"
                );
            }
            // Cancellation, last: the gate stays requested.
            cx.cancel_gate().request();
            assert_eq!(attempt(&old, DT_S, &[1.0]), PhaseStateError::EntropyEvaluationCancelled);
        });
    }
}
