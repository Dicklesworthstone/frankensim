//! A card's Ohmic conductor behind resistor ports, its resistance re-read at
//! the heat port's temperature. A manufactured resistivity table with three
//! knots, whose resistance and slope the tests write out by hand; no measured
//! material is claimed.
//!
//! Claims, written before the first run:
//! - Both terminal forms answer the card's `rho(T) L / A` at the heat port's
//!   temperature: outputs within 1e-12 of the hand values, the current form's
//!   equal to the node's own graph evaluation to the bit, receipts of the
//!   Joule heat, and the conformance kit and its tangent check (relative
//!   2e-5) at points inside the segments. At the interior knot the
//!   temperature column of the tangent is the right segment's slope, by hand.
//! - A heater on a fixed voltage through the kernel, its heat into a 1 J/K
//!   body, follows a hand-written backward step that bisects the body's end
//!   temperature with the resistance at the logarithmic mean temperature the
//!   body's port publishes: within 2e-9 K at each of 60 steps, across the
//!   knot. The current falls as the resistance rises. The boundary work, the
//!   Joule heat and the stored energy agree within 1e-9 J.
//! - In a current loop with an inductance, the inductance's lost energy is
//!   the body's gain, within 1e-9 J, and the current decays.
//! - Refusals publish nothing and are the node's `ConductorError`: a
//!   temperature off the card's interval (the node's own refusal), a value
//!   that is not finite, a timestep that is not positive, a wrong slot count,
//!   and cancellation, each named. Through the kernel a
//!   heater that would leave the interval is refused and the body keeps its
//!   state to the bit.
//! - The identity binds the card, the geometry and the terminal form.
//!
//! Added after the first scene run found every seed card with a resistivity
//! refused, its provenance flags (`engineering_reference_model = 1` and the
//! like) counted as unowned coordinates:
//! - A validity flag declared dimensionless and bounded to one value is held
//!   at its query value: the node answers the unflagged card's values to the
//!   bit. An untyped flag, and a dimensionless range, still refuse.

#![cfg(feature = "component")]
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::component::{OhmicResistor, Terminals};
use fs_material::conductor::graph::MaterialOhmicNode;
use fs_material::conductor::{ConductorError, ElectricalPropertyBasis, ResolvedConductorState};
use fs_material::graph::LawNode;
use fs_material::phase::component::UniformEnthalpyBody;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortTangent};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_system::elements::{Mass, Motion};
use fs_system::{SystemBuilder, SystemError};

/// Resistivity [ohm m] at 300, 350 and 400 K.
const KNOTS: [(f64, f64); 3] = [(300.0, 2.0e-7), (350.0, 3.0e-7), (400.0, 3.5e-7)];
const LENGTH_M: f64 = 2.0;
const AREA_M2: f64 = 1.0e-6;

fn in_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 61,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + absolute,
        "actual {actual:e}, expected {expected:e}, relative {relative:e}, absolute {absolute:e}"
    );
}

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

/// The manufactured card. With `flag`, the resistivity's validity also bounds
/// `engineering_reference_model` to the given interval, which the query
/// meets at 1, as the seed cards' provenance flags do: declared with the
/// given quantity, or untyped for `None`.
fn card(
    knots: &[(f64, f64)],
    flag: Option<(Option<QuantitySpec>, [f64; 2])>,
) -> (MaterialCard, ResolvedConductorState) {
    const FLAG: &str = "engineering_reference_model";
    let basis = ElectricalPropertyBasis::Resistivity;
    let key = PropertyKey::new("electrical-resistivity", basis.dimensions());
    let (low, high) = (knots[0].0, knots[knots.len() - 1].0);
    let mut validity = ValidityDomain::unconstrained().with_quantity("T", absolute(), low, high);
    let mut query = QueryPoint::new()
        .with_quantity("T", absolute(), low)
        .unwrap();
    match flag {
        Some((Some(quantity), [from, to])) => {
            validity = validity.with_quantity(FLAG, quantity, from, to);
            query = query.with_quantity(FLAG, quantity, 1.0).unwrap();
        }
        Some((None, [from, to])) => {
            validity = validity.with(FLAG, from, to);
            query = query.with(FLAG, 1.0).unwrap();
        }
        None => {}
    }
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: knots.to_vec(),
                dims: basis.dimensions(),
            },
            interpolation: InterpolationPolicy::LinearInside,
            validity,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: "manufactured resistivity table".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductor".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &query,
        &[basis.requirement_for(&key).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original =
        ResolvedConductorState::try_from_material_with_property(material, basis, &key).unwrap();
    (card, original)
}

fn node_of(knots: &[(f64, f64)], length_m: f64, area_m2: f64) -> MaterialOhmicNode {
    let (card, original) = card(knots, None);
    let range = [knots[0].0, knots[knots.len() - 1].0];
    MaterialOhmicNode::try_new(&card, original, "T", range, length_m, area_m2).unwrap()
}

fn resistor(terminals: Terminals) -> OhmicResistor {
    OhmicResistor::new(node_of(&KNOTS, LENGTH_M, AREA_M2), terminals).unwrap()
}

/// The segment a temperature reads: the right one at an interior knot.
fn segment(kelvin: f64) -> usize {
    KNOTS
        .iter()
        .rposition(|&(knot, _)| knot <= kelvin)
        .unwrap()
        .min(KNOTS.len() - 2)
}

/// `rho(T) L / A` by hand.
fn resistance(kelvin: f64) -> f64 {
    let i = segment(kelvin);
    let ((t0, r0), (t1, r1)) = (KNOTS[i], KNOTS[i + 1]);
    (r0 + (r1 - r0) * (kelvin - t0) / (t1 - t0)) * LENGTH_M / AREA_M2
}

/// `dR/dT` by hand [ohm/K].
fn resistance_slope(kelvin: f64) -> f64 {
    let i = segment(kelvin);
    let ((t0, r0), (t1, r1)) = (KNOTS[i], KNOTS[i + 1]);
    (r1 - r0) / (t1 - t0) * LENGTH_M / AREA_M2
}

fn observed(
    component: &OhmicResistor,
    inputs: &[f64],
    cx: &Cx<'_>,
) -> (Vec<f64>, fs_ports::component::Receipt) {
    let case = Case {
        state: &(),
        stale: &(),
        dt_s: 0.25,
        inputs,
    };
    let Observed::Stepped {
        outputs, receipt, ..
    } = check_component(component, &case, cx).unwrap()
    else {
        panic!("a resistor does not shorten time");
    };
    (outputs, receipt)
}

#[test]
fn both_terminal_forms_answer_the_cards_resistance_at_the_heat_port() {
    let law = node_of(&KNOTS, LENGTH_M, AREA_M2);
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for terminals in [Terminals::Current, Terminals::Voltage] {
            let component = resistor(terminals);
            for (a, b, kelvin) in [
                (1.5, -0.5, 300.0),
                (-2.0, 0.25, 325.0),
                (0.8, 0.8, 330.0),
                (3.0, 1.0, 350.0),
                (0.5, 2.5, 399.0),
            ] {
                let inputs = [a, b, kelvin];
                let (outputs, receipt) = observed(&component, &inputs, cx);
                let r = resistance(kelvin);
                let (output, power) = match terminals {
                    Terminals::Current => ((a - b) * r, (a - b) * (a - b) * r),
                    Terminals::Voltage => ((a - b) / r, (a - b) * (a - b) / r),
                };
                close(outputs[0], output, 1e-12, 1e-300);
                assert_eq!(outputs[1], -outputs[0]);
                close(outputs[2], -power, 1e-12, 1e-300);
                assert_eq!(receipt.dissipated_j, 0.25 * -outputs[2]);
                assert_eq!(
                    receipt.entropy_generated_j_per_k,
                    receipt.dissipated_j / kelvin
                );
                if terminals == Terminals::Current {
                    let own = law.evaluate(&[], &[kelvin - 300.0, a - b]).unwrap().outputs;
                    assert_eq!(outputs[0].to_bits(), own[0].to_bits());
                    assert_eq!((-outputs[2]).to_bits(), own[1].to_bits());
                }
                // The temperature column of the tangent, by hand (the right
                // segment at the knot).
                let mut pushed = [0.0; 3];
                component
                    .jvp(&(), 0.25, &inputs, &[0.0, 0.0, 1.0], &mut pushed, cx)
                    .unwrap();
                let (current, slope) = match terminals {
                    Terminals::Current => (a - b, resistance_slope(kelvin)),
                    Terminals::Voltage => ((a - b) / r, resistance_slope(kelvin)),
                };
                let expected = match terminals {
                    Terminals::Current => [
                        current * slope,
                        -current * slope,
                        -current * current * slope,
                    ],
                    Terminals::Voltage => [
                        -current * slope / r,
                        current * slope / r,
                        current * current * slope,
                    ],
                };
                for (got, want) in pushed.iter().zip(expected) {
                    close(*got, want, 1e-12, 1e-300);
                }
                // Inside a segment the kit's difference quotients agree.
                if kelvin != 350.0 && kelvin != 300.0 {
                    let case = Case {
                        state: &(),
                        stale: &(),
                        dt_s: 0.25,
                        inputs: &inputs,
                    };
                    check_port_tangent(&component, &case, 2e-5, 1e-10, cx).unwrap();
                }
            }
        }
    });
}

/// One backward step of a 1 J/K body heated by `V^2 / R` at the logarithmic
/// mean of its two end temperatures, which its heat port publishes: the end
/// temperature by bisection, and the Joule heat [J].
fn heater_oracle(kelvin: f64, voltage: f64, dt: f64) -> (f64, f64) {
    let power = |next: f64| {
        let increment = next - kelvin;
        let mean = if increment == 0.0 {
            kelvin
        } else {
            increment / (increment / kelvin).ln_1p()
        };
        voltage * voltage / resistance(mean)
    };
    // The power falls as the body warms, so the step at the start's power
    // bounds the end temperature from above.
    let (mut lower, mut upper) = (kelvin, kelvin + dt * power(kelvin));
    for _ in 0..80 {
        let trial = 0.5 * lower + 0.5 * upper;
        if trial - kelvin < dt * power(trial) {
            lower = trial;
        } else {
            upper = trial;
        }
    }
    let next = 0.5 * lower + 0.5 * upper;
    (next, dt * power(next))
}

/// A solid chart of 10 J/(kg K) from 280 K to 4000 K, `h = 0` at 300 K: wide
/// enough that the body never refuses before the wire does.
fn chart() -> EquilibriumEnthalpyPhaseCurve {
    let knot = |kelvin: f64| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: 10.0 * (kelvin - 300.0),
        temperature_k: kelvin,
        liquid_mass_fraction: 0.0,
        bulk_density_kg_m3: 8000.0,
    };
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x5a; 32]),
        SolidLiquidPhase::Solid,
        vec![knot(280.0), knot(4000.0)],
    )
    .unwrap()
}

#[test]
fn a_heater_on_a_fixed_voltage_follows_its_rising_resistance_as_the_hand_step_does() {
    let chart = chart();
    let voltage = 1.0;
    let mut builder = SystemBuilder::new();
    let wire = builder
        .add_with_tangent("wire", resistor(Terminals::Voltage), ())
        .unwrap();
    // 0.1 kg at 10 J/(kg K): 1 J/K.
    let body = builder
        .add_with_tangent(
            "body",
            UniformEnthalpyBody::new(chart.clone(), 0.1).unwrap(),
            chart.state_at_specific_enthalpy(0.0).unwrap(),
        )
        .unwrap();
    let (a, b, heat) = (
        builder.port(wire, "a").unwrap(),
        builder.port(wire, "b").unwrap(),
        builder.port(wire, "heat").unwrap(),
    );
    builder
        .boundary(a, Box::new(move |_, _, values| values[0] = voltage))
        .unwrap();
    builder
        .boundary(b, Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .common_effort(builder.port(body, "heat").unwrap(), &[heat])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));
    let (mut kelvin, mut joule) = (300.0, 0.0);
    let mut currents = Vec::new();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..60 {
            let (next, step_heat) = heater_oracle(kelvin, voltage, 0.5);
            kelvin = next;
            joule += step_heat;
            system.step(0.5, cx).unwrap();
            let state = system.state::<EquilibriumPhaseState>(body).unwrap();
            close(state.temperature_k(), kelvin, 0.0, 2e-9);
            currents.push(system.last_output(a)[0]);
        }
    });
    assert!(
        kelvin > 355.0,
        "the run must cross the 350 K knot: {kelvin}"
    );
    assert!(
        currents.windows(2).all(|pair| pair[1] < pair[0]),
        "the current falls as the wire warms"
    );
    assert!(currents[59] < 0.75 * currents[0], "{currents:?}");
    let ledger = system.ledger();
    let stored = ledger.nodes[body.index()].stored_change_j.value();
    let dissipated = ledger.nodes[wire.index()].dissipated_j.value();
    let work: f64 = ledger
        .boundaries
        .iter()
        .map(|boundary| boundary.energy_in_j.value())
        .sum();
    close(stored, joule, 0.0, 1e-9);
    close(dissipated, stored, 0.0, 1e-9);
    close(work, stored, 0.0, 1e-9);
    close(ledger.energy_residual_j(), 0.0, 0.0, 1e-9);
    eprintln!(
        "heater: 60 steps of 0.5 s at {voltage} V, body {kelvin} K, current {} A -> {} A, stored {stored} J, boundary work {work} J",
        currents[0], currents[59]
    );
}

#[test]
fn in_a_current_loop_the_inductances_lost_energy_is_the_bodys_gain() {
    let chart = chart();
    let inductance = 0.5;
    let mut builder = SystemBuilder::new();
    let coil = builder
        .add_with_tangent(
            "coil",
            Mass::in_domain(PortKind::ElectricalVoltageCurrent, inductance).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 4.0,
            },
        )
        .unwrap();
    let wire = builder
        .add_with_tangent("wire", resistor(Terminals::Current), ())
        .unwrap();
    let body = builder
        .add_with_tangent(
            "body",
            UniformEnthalpyBody::new(chart.clone(), 0.1).unwrap(),
            chart.state_at_specific_enthalpy(0.0).unwrap(),
        )
        .unwrap();
    builder
        .bond(
            builder.port(coil, "p").unwrap(),
            builder.port(wire, "a").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(wire, "b").unwrap(),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    let heat = builder.port(wire, "heat").unwrap();
    builder
        .common_effort(builder.port(body, "heat").unwrap(), &[heat])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..40 {
            system.step(0.05, cx).unwrap();
        }
    });
    let current = system.state::<Motion>(coil).unwrap().velocity_m_s;
    let lost = 0.5 * inductance * (4.0 * 4.0 - current * current);
    let ledger = system.ledger();
    let gained = ledger.nodes[body.index()].stored_change_j.value();
    assert!(
        current > 0.0 && current < 1.0,
        "the current decays: {current}"
    );
    close(gained, lost, 0.0, 1e-9);
    close(ledger.energy_residual_j(), 0.0, 0.0, 1e-9);
    let kelvin = system
        .state::<EquilibriumPhaseState>(body)
        .unwrap()
        .temperature_k();
    eprintln!(
        "loop: 40 steps of 0.05 s, current 4 A -> {current} A, body {kelvin} K, lost {lost} J, gained {gained} J"
    );
}

#[test]
fn refusals_publish_nothing_and_the_kernel_keeps_the_body() {
    let component = resistor(Terminals::Voltage);
    in_cx(&CancelGate::new_clock_free(), |cx| {
        let attempt = |dt: f64, inputs: &[f64], slots: usize| {
            let mut outputs = vec![f64::NAN; slots];
            let result = component.propose(&(), dt, inputs, &mut (), &mut outputs, cx);
            assert!(outputs.iter().all(|value| value.is_nan()), "{outputs:?}");
            result
        };
        let refused = |quantity: &'static str| Err(ConductorError::InvalidInput { quantity });
        for kelvin in [299.0, 400.5] {
            assert_eq!(
                attempt(0.25, &[1.0, 0.0, kelvin], 3),
                refused("Ohmic graph temperature outside its admitted source interval")
            );
        }
        for inputs in [[f64::NAN, 0.0, 320.0], [1.0, f64::INFINITY, 320.0]] {
            assert_eq!(
                attempt(0.25, &inputs, 3),
                refused("Ohmic resistor needs three finite port values")
            );
        }
        for dt in [0.0, -1.0, f64::NAN] {
            assert_eq!(
                attempt(dt, &[1.0, 0.0, 320.0], 3),
                refused("Ohmic resistor timestep must be positive and finite")
            );
        }
        assert_eq!(
            attempt(0.25, &[1.0, 0.0, 320.0], 2),
            refused("Ohmic resistor needs three output slots")
        );
        cx.cancel_gate().request();
        assert_eq!(
            attempt(0.25, &[1.0, 0.0, 320.0], 3),
            refused("Ohmic resistor cancelled")
        );
    });
    // 10 V across at most the card's 0.7 ohm is at least 143 W into 1 J/K:
    // one 5 s step ends above 1000 K, where the logarithmic mean the body's
    // port publishes is past 586 K. No resistance in the card's interval
    // closes that step, so the kernel refuses it and names the wire. (A first
    // draft drove 3 V; the rising resistance closes that step near 370 K, a
    // mean of about 333 K, and the kernel took it.)
    let chart = chart();
    let mut builder = SystemBuilder::new();
    let wire = builder
        .add_with_tangent("wire", resistor(Terminals::Voltage), ())
        .unwrap();
    let body = builder
        .add_with_tangent(
            "body",
            UniformEnthalpyBody::new(chart.clone(), 0.1).unwrap(),
            chart.state_at_specific_enthalpy(0.0).unwrap(),
        )
        .unwrap();
    let heat = builder.port(wire, "heat").unwrap();
    builder
        .boundary(
            builder.port(wire, "a").unwrap(),
            Box::new(|_, _, values| values[0] = 10.0),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(wire, "b").unwrap(),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    builder
        .common_effort(builder.port(body, "heat").unwrap(), &[heat])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));
    let before = *system.state::<EquilibriumPhaseState>(body).unwrap();
    let error = in_cx(&CancelGate::new_clock_free(), |cx| system.step(5.0, cx)).unwrap_err();
    eprintln!("refused: {error}");
    assert!(
        matches!(&error, SystemError::Component { node, .. } if node == "wire"),
        "{error:?}"
    );
    let after = *system.state::<EquilibriumPhaseState>(body).unwrap();
    assert_eq!(
        after.specific_enthalpy_j_kg().to_bits(),
        before.specific_enthalpy_j_kg().to_bits()
    );
    assert_eq!(system.steps(), 0);
}

#[test]
fn the_identity_binds_the_card_the_geometry_and_the_terminal_form() {
    let words = |component: &OhmicResistor| {
        let mut words = Vec::new();
        component.identity_words(&mut words);
        words
    };
    let base = words(&resistor(Terminals::Current));
    assert_eq!(base, words(&resistor(Terminals::Current)));
    for other in [
        resistor(Terminals::Voltage),
        OhmicResistor::new(node_of(&KNOTS, 2.5, AREA_M2), Terminals::Current).unwrap(),
        OhmicResistor::new(node_of(&KNOTS, LENGTH_M, 2.0e-6), Terminals::Current).unwrap(),
        OhmicResistor::new(
            node_of(
                &[(300.0, 2.0e-7), (350.0, 3.1e-7), (400.0, 3.5e-7)],
                LENGTH_M,
                AREA_M2,
            ),
            Terminals::Current,
        )
        .unwrap(),
    ] {
        assert_ne!(base, words(&other));
    }
    assert!(matches!(
        OhmicResistor::STATE_CHARTER,
        Some(charter) if charter.owner == "fs-material::conductor::component::OhmicResistor"
    ));
    let component = resistor(Terminals::Current);
    let mut encoded = Vec::new();
    assert!(component.encode_state(&(), &mut encoded));
    assert!(encoded.is_empty());
    assert_eq!(component.decode_state(&[]), Some(()));
    assert_eq!(component.decode_state(&[0]), None);
}

#[test]
fn a_dimensionless_validity_flag_is_held_and_an_untyped_or_ranged_one_still_refuses() {
    let range = [KNOTS[0].0, KNOTS[2].0];
    let dimensionless = Some(QuantitySpec::dimensional(Dims([0; 6])));
    let flagged = |quantity, bound| {
        let (card, original) = card(&KNOTS, Some((quantity, bound)));
        MaterialOhmicNode::try_new(&card, original, "T", range, LENGTH_M, AREA_M2)
    };
    let plain = node_of(&KNOTS, LENGTH_M, AREA_M2);
    let held = flagged(dimensionless, [1.0, 1.0]).unwrap();
    for kelvin in [300.0, 325.0, 350.0, 399.0] {
        let inputs = [kelvin - 300.0, 1.5];
        let want = plain.evaluate(&[], &inputs).unwrap().outputs;
        let got = held.evaluate(&[], &inputs).unwrap().outputs;
        assert_eq!(
            got.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            want.iter().map(|value| value.to_bits()).collect::<Vec<_>>()
        );
    }
    // An untyped bound has no known dimension: it could be a frequency or a
    // field, which the conductor must not freeze. A dimensionless range, such
    // as a phase fraction over [0, 1], is a coordinate, not a flag.
    let error = flagged(None, [1.0, 1.0]).unwrap_err();
    eprintln!("untyped flag refused: {error}");
    let error = flagged(dimensionless, [0.0, 1.0]).unwrap_err();
    eprintln!("ranged flag refused: {error}");
}
