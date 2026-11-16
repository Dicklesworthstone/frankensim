//! A card's Newtonian liquid in a laminar circular pipe behind two hydraulic
//! ports and one heat port, its viscosity re-read at the heat port's
//! temperature. Manufactured viscosity and density tables, whose conductance
//! and slopes the tests write out by hand; no measured liquid is claimed.
//!
//! Claims, written before the first run:
//! - With a dynamic-viscosity card the outputs are `Q = pi r^4 dp / (8 mu L)`,
//!   `-Q` and `-Q dp` at the heat port's temperature within 1e-12 of the hand
//!   values, receipts of the viscous heat, and the conformance kit and its
//!   tangent check (relative 2e-5) at points inside the segments. The
//!   temperature column of the tangent is `-Q mu'/mu`, `Q mu'/mu`,
//!   `Q dp mu'/mu` by hand, with the right segment's slope at the interior
//!   knot and the left one's at the last knot.
//! - With a kinematic-viscosity card and a density on other knots, the
//!   viscosity is `rho nu` and its slope `rho' nu + rho nu'`, each on its own
//!   line's segment, within 1e-12 by hand, at a knot of either line too.
//! - A pipe on a fixed pressure drop through the kernel, its heat into a
//!   1 J/K body, follows a hand-written backward step that bisects the body's
//!   end temperature with the viscosity at the logarithmic mean temperature
//!   the body's port publishes: within 2e-9 K at each of 60 steps, across the
//!   knot. The flow rises as the liquid thins. The boundary work, the viscous
//!   heat and the stored energy agree within 1e-9 J.
//! - Preparation refuses: a dynamic viscosity with a density, a kinematic one
//!   without, a property that is neither, a range that does not hold the
//!   original temperature, a geometry that is not positive and finite, and a
//!   viscosity that is not positive at a knot inside the range (the card's
//!   own admissibility, named).
//! - Refusals publish nothing and are the pipe's `MaterialStatePointError`:
//!   a temperature off the range, a value that is not finite, a timestep that
//!   is not positive, a wrong slot count and cancellation, each named.
//!   Through the kernel a pipe that would heat its body past the range is
//!   refused and the body keeps its state to the bit.
//! - The identity binds the card, the geometry and the viscosity basis.

#![cfg(feature = "component")]
#![allow(clippy::float_cmp)]

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::laminar_pipe::{
    LaminarPipeComponent, MaterialLaminarPipe, ViscosityBasis,
};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY,
};
use fs_material::phase::component::UniformEnthalpyBody;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortTangent};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_system::{NodeId, PortRef, System, SystemBuilder, SystemError};

/// Dynamic viscosity [Pa s] at 300, 350 and 400 K: a liquid thinning as it
/// warms.
const MU: [(f64, f64); 3] = [(300.0, 1.0e-3), (350.0, 5.0e-4), (400.0, 4.0e-4)];
/// Kinematic viscosity [m²/s] and density [kg/m³], on different knots.
const NU: [(f64, f64); 3] = [(300.0, 1.2e-6), (350.0, 6.0e-7), (400.0, 4.5e-7)];
const RHO: [(f64, f64); 3] = [(300.0, 1000.0), (325.0, 990.0), (400.0, 960.0)];
const RADIUS_M: f64 = 1.0e-3;
const LENGTH_M: f64 = 1.0;
const DT: f64 = 0.25;

/// A named curve of the stated dimensions over `T`: `(T [K], value)` knots.
type Curve<'a> = (&'a str, Dims, &'a [(f64, f64)]);

fn in_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 67,
                kernel_id: 29,
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

/// A manufactured liquid card of the named curves over `T` and the state
/// point of all of them at `reference_k`, each strictly positive.
fn card(
    curves: &[Curve<'_>],
    reference_k: f64,
) -> (
    Arc<MaterialCard>,
    fs_material::state_point::ResolvedMaterialStatePoint,
) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, knots) in curves {
        let (low, high) = (knots[0].0, knots[knots.len() - 1].0);
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: knots.to_vec(),
                    dims,
                },
                interpolation: InterpolationPolicy::LinearInside,
                validity: ValidityDomain::unconstrained().with_quantity("T", absolute(), low, high),
                uncertainty: UncertaintyModel::Unstated,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured liquid table".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured liquid".into(),
            phase: "liquid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("T", absolute(), reference_k)
        .unwrap();
    let original = resolve_material_state_point(
        &card,
        &query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (Arc::new(card), original)
}

fn dynamic_curves(mu: &[(f64, f64)]) -> Vec<Curve<'_>> {
    vec![(DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu)]
}

fn kinematic_curves() -> Vec<Curve<'static>> {
    vec![
        (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS, &NU),
        ("density", Density::DIMS, &RHO),
    ]
}

fn pipe_of(
    curves: &[Curve<'_>],
    radius_m: f64,
    length_m: f64,
) -> Result<MaterialLaminarPipe, MaterialStatePointError> {
    let (card, original) = card(curves, 300.0);
    MaterialLaminarPipe::new(card, original, "T", [300.0, 400.0], radius_m, length_m)
}

fn dynamic() -> LaminarPipeComponent {
    LaminarPipeComponent::new(pipe_of(&dynamic_curves(&MU), RADIUS_M, LENGTH_M).unwrap())
}

fn kinematic() -> LaminarPipeComponent {
    LaminarPipeComponent::new(pipe_of(&kinematic_curves(), RADIUS_M, LENGTH_M).unwrap())
}

/// A line's value and slope at `kelvin`, by hand: the right segment at an
/// interior knot, the left one at the last knot.
fn on_line(knots: &[(f64, f64)], kelvin: f64) -> (f64, f64) {
    let i = knots
        .iter()
        .rposition(|&(knot, _)| knot <= kelvin)
        .unwrap()
        .min(knots.len() - 2);
    let ((t0, v0), (t1, v1)) = (knots[i], knots[i + 1]);
    let slope = (v1 - v0) / (t1 - t0);
    (v0 + slope * (kelvin - t0), slope)
}

/// `pi r^4 / (8 L)` by hand [m³].
fn geometry() -> f64 {
    core::f64::consts::PI * RADIUS_M.powi(4) / (8.0 * LENGTH_M)
}

fn observed(
    component: &LaminarPipeComponent,
    inputs: &[f64],
    cx: &Cx<'_>,
) -> (Vec<f64>, fs_ports::component::Receipt) {
    let case = Case {
        state: &(),
        stale: &(),
        dt_s: DT,
        inputs,
    };
    let Observed::Stepped {
        outputs, receipt, ..
    } = check_component(component, &case, cx).unwrap()
    else {
        panic!("a pipe does not shorten time");
    };
    (outputs, receipt)
}

/// Both bases against `viscosity(T) -> (mu, mu')` by hand, at the knots and
/// inside the segments.
fn answers_by_hand(component: &LaminarPipeComponent, viscosity: impl Fn(f64) -> (f64, f64)) {
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for (a, b, kelvin) in [
            (2.0e4, 0.0, 300.0),
            (-1.5e4, 5.0e3, 320.0),
            (8.0e3, 8.0e3, 330.0),
            (3.0e4, 1.0e4, 325.0),
            (3.0e4, -1.0e4, 350.0),
            (1.0e3, 2.5e4, 380.0),
            (5.0e4, 1.0e4, 400.0),
        ] {
            let inputs = [a, b, kelvin];
            let (outputs, receipt) = observed(component, &inputs, cx);
            let (mu, mu_slope) = viscosity(kelvin);
            let (conductance, conductance_slope) =
                component.pipe().conductance_and_slope(kelvin).unwrap();
            close(conductance, geometry() / mu, 1e-12, 0.0);
            close(
                conductance_slope,
                -geometry() * mu_slope / (mu * mu),
                1e-12,
                0.0,
            );
            let drop = a - b;
            let flow = geometry() / mu * drop;
            close(outputs[0], flow, 1e-12, 1e-300);
            assert_eq!(outputs[1], -outputs[0]);
            close(outputs[2], -flow * drop, 1e-12, 1e-300);
            assert_eq!(receipt.dissipated_j, DT * -outputs[2]);
            assert_eq!(receipt.stored_change_j, 0.0);
            assert_eq!(
                receipt.entropy_generated_j_per_k,
                receipt.dissipated_j / kelvin
            );
            // The temperature column of the tangent, by hand.
            let mut pushed = [0.0; 3];
            component
                .jvp(&(), DT, &inputs, &[0.0, 0.0, 1.0], &mut pushed, cx)
                .unwrap();
            let relative = mu_slope / mu;
            for (got, want) in
                pushed
                    .iter()
                    .zip([-flow * relative, flow * relative, flow * drop * relative])
            {
                close(*got, want, 1e-12, 1e-300);
            }
            // The pressure columns: G, -G and -2 G dp.
            component
                .jvp(&(), DT, &inputs, &[1.0, 0.0, 0.0], &mut pushed, cx)
                .unwrap();
            let g = geometry() / mu;
            for (got, want) in pushed.iter().zip([g, -g, -2.0 * g * drop]) {
                close(*got, want, 1e-12, 1e-300);
            }
            // Inside the segments of every line the kit's difference
            // quotients agree.
            if ![300.0, 325.0, 350.0, 400.0].contains(&kelvin) {
                let case = Case {
                    state: &(),
                    stale: &(),
                    dt_s: DT,
                    inputs: &inputs,
                };
                check_port_tangent(component, &case, 2e-5, 1e-10, cx).unwrap();
            }
        }
    });
}

#[test]
fn a_dynamic_viscosity_card_answers_hagen_poiseuille_at_the_heat_port() {
    let component = dynamic();
    assert_eq!(component.pipe().basis(), ViscosityBasis::Dynamic);
    answers_by_hand(&component, |kelvin| on_line(&MU, kelvin));
    let (mu, slope) = component.pipe().viscosity_and_slope(350.0).unwrap();
    assert_eq!(mu, 5.0e-4);
    // The right segment at the interior knot, the left one at the last.
    close(slope, (4.0e-4 - 5.0e-4) / 50.0, 1e-15, 0.0);
    let (_, last) = component.pipe().viscosity_and_slope(400.0).unwrap();
    close(last, (4.0e-4 - 5.0e-4) / 50.0, 1e-15, 0.0);
    let (_, first) = component.pipe().viscosity_and_slope(300.0).unwrap();
    close(first, (5.0e-4 - 1.0e-3) / 50.0, 1e-15, 0.0);
}

#[test]
fn a_kinematic_viscosity_card_answers_rho_nu_with_the_product_rules_slope() {
    let component = kinematic();
    assert_eq!(component.pipe().basis(), ViscosityBasis::Kinematic);
    let viscosity = |kelvin: f64| {
        let (nu, nu_slope) = on_line(&NU, kelvin);
        let (rho, rho_slope) = on_line(&RHO, kelvin);
        (rho * nu, rho_slope * nu + rho * nu_slope)
    };
    for kelvin in [300.0, 310.0, 325.0, 340.0, 350.0, 375.0, 400.0] {
        let (mu, slope) = component.pipe().viscosity_and_slope(kelvin).unwrap();
        let (want, want_slope) = viscosity(kelvin);
        close(mu, want, 1e-12, 0.0);
        close(slope, want_slope, 1e-12, 0.0);
    }
    answers_by_hand(&component, viscosity);
}

/// One backward step of a 1 J/K body heated by `G dp^2` at the logarithmic
/// mean of its two end temperatures, which its heat port publishes: the end
/// temperature by bisection, and the viscous heat [J].
fn heater_oracle(kelvin: f64, drop: f64, dt: f64) -> (f64, f64) {
    let power = |next: f64| {
        let increment = next - kelvin;
        let mean = if increment == 0.0 {
            kelvin
        } else {
            increment / (increment / kelvin).ln_1p()
        };
        geometry() / on_line(&MU, mean).0 * drop * drop
    };
    // The power rises as the liquid thins, so the step at the start's power
    // bounds the end temperature from below, and twice it from above.
    let (mut lower, mut upper) = (
        kelvin + dt * power(kelvin),
        kelvin + 2.0 * dt * power(kelvin),
    );
    assert!(
        upper - kelvin > dt * power(upper),
        "the bracket holds the step"
    );
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
/// enough that the body never refuses before the pipe does.
fn chart() -> EquilibriumEnthalpyPhaseCurve {
    let knot = |kelvin: f64| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: 10.0 * (kelvin - 300.0),
        temperature_k: kelvin,
        liquid_mass_fraction: 0.0,
        bulk_density_kg_m3: 8000.0,
    };
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x5b; 32]),
        SolidLiquidPhase::Solid,
        vec![knot(280.0), knot(4000.0)],
    )
    .unwrap()
}

/// A pipe on a fixed pressure drop, its heat into a 1 J/K body.
fn heated(drop: f64) -> (System, NodeId, NodeId, PortRef) {
    let chart = chart();
    let mut builder = SystemBuilder::new();
    let pipe = builder.add_with_tangent("pipe", dynamic(), ()).unwrap();
    // 0.1 kg at 10 J/(kg K): 1 J/K.
    let body = builder
        .add_with_tangent(
            "body",
            UniformEnthalpyBody::new(chart.clone(), 0.1).unwrap(),
            chart.state_at_specific_enthalpy(0.0).unwrap(),
        )
        .unwrap();
    let (a, b, heat) = (
        builder.port(pipe, "a").unwrap(),
        builder.port(pipe, "b").unwrap(),
        builder.port(pipe, "heat").unwrap(),
    );
    builder
        .boundary(a, Box::new(move |_, _, values| values[0] = drop))
        .unwrap();
    builder
        .boundary(b, Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .common_effort(builder.port(body, "heat").unwrap(), &[heat])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));
    (system, pipe, body, a)
}

#[test]
fn a_pipe_on_a_fixed_drop_follows_its_thinning_liquid_as_the_hand_step_does() {
    // 6e4 Pa across: 1.41 W at 300 K, 2.83 W at 350 K, into 1 J/K.
    let drop = 6.0e4;
    let (mut system, pipe, body, a) = heated(drop);
    let (mut kelvin, mut heat) = (300.0, 0.0);
    let mut flows = Vec::new();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..60 {
            let (next, step_heat) = heater_oracle(kelvin, drop, 0.5);
            kelvin = next;
            heat += step_heat;
            system.step(0.5, cx).unwrap();
            let state = system.state::<EquilibriumPhaseState>(body).unwrap();
            close(state.temperature_k(), kelvin, 0.0, 2e-9);
            flows.push(system.last_output(a)[0]);
        }
    });
    assert!(
        (352.0..400.0).contains(&kelvin),
        "the run must cross the 350 K knot and stay in the card: {kelvin}"
    );
    assert!(
        flows.windows(2).all(|pair| pair[1] > pair[0]),
        "the flow rises as the liquid thins"
    );
    assert!(flows[59] > 1.5 * flows[0], "{flows:?}");
    let ledger = system.ledger();
    let stored = ledger.nodes[body.index()].stored_change_j.value();
    let dissipated = ledger.nodes[pipe.index()].dissipated_j.value();
    let work: f64 = ledger
        .boundaries
        .iter()
        .map(|boundary| boundary.energy_in_j.value())
        .sum();
    close(stored, heat, 0.0, 1e-9);
    close(dissipated, stored, 0.0, 1e-9);
    close(work, stored, 0.0, 1e-9);
    close(ledger.energy_residual_j(), 0.0, 0.0, 1e-9);
    eprintln!(
        "pipe: 60 steps of 0.5 s at {drop} Pa, body {kelvin} K, flow {} -> {} m3/s, stored {stored} J, boundary work {work} J",
        flows[0], flows[59]
    );
}

#[test]
fn preparation_refuses_what_the_pipe_cannot_answer() {
    let refused =
        |quantity: &'static str| Err(MaterialStatePointError::InvalidDerived { quantity });
    let mut with_density = dynamic_curves(&MU);
    with_density.push(("density", Density::DIMS, &RHO));
    assert_eq!(
        pipe_of(&with_density, RADIUS_M, LENGTH_M).map(|_| ()),
        refused("laminar pipe with a dynamic viscosity takes no density")
    );
    assert_eq!(
        pipe_of(
            &[(KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS, &NU)],
            RADIUS_M,
            LENGTH_M
        )
        .map(|_| ()),
        refused("laminar pipe with a kinematic viscosity needs the card's density")
    );
    let mut with_other = dynamic_curves(&MU);
    with_other.push((
        "thermal_conductivity",
        Dims([1, 1, -3, -1, 0, 0]),
        &[(300.0, 0.6), (400.0, 0.7)],
    ));
    assert_eq!(
        pipe_of(&with_other, RADIUS_M, LENGTH_M).map(|_| ()),
        refused("laminar pipe source holds a property that is neither a viscosity nor a density")
    );
    assert_eq!(
        pipe_of(&[("density", Density::DIMS, &RHO)], RADIUS_M, LENGTH_M).map(|_| ()),
        refused("laminar pipe source holds no viscosity")
    );
    for (radius, length) in [
        (0.0, LENGTH_M),
        (-1.0e-3, LENGTH_M),
        (RADIUS_M, 0.0),
        (RADIUS_M, f64::INFINITY),
        (f64::NAN, LENGTH_M),
    ] {
        assert_eq!(
            pipe_of(&dynamic_curves(&MU), radius, length).map(|_| ()),
            refused("laminar pipe radius and length must be positive and finite")
        );
    }
    let (card, original) = card(&dynamic_curves(&MU), 300.0);
    for range in [[310.0, 400.0], [300.0, 300.0], [0.0, 400.0]] {
        assert_eq!(
            MaterialLaminarPipe::new(
                card.clone(),
                original.clone(),
                "T",
                range,
                RADIUS_M,
                LENGTH_M
            )
            .map(|_| ()),
            refused("laminar pipe range must contain the original positive temperature")
        );
    }
    assert_eq!(
        MaterialLaminarPipe::new(
            card.clone(),
            original.clone(),
            "p",
            [300.0, 400.0],
            RADIUS_M,
            LENGTH_M
        )
        .map(|_| ()),
        refused("laminar pipe source temperature is absent")
    );
    // A range past the card's own support is the card's refusal.
    assert!(
        MaterialLaminarPipe::new(card, original, "T", [300.0, 450.0], RADIUS_M, LENGTH_M).is_err()
    );
    // A viscosity that reaches zero inside the range: the original point is
    // admitted, the knot is not.
    let thinning = [(300.0, 1.0e-3), (350.0, 5.0e-4), (400.0, 0.0)];
    assert_eq!(
        pipe_of(&dynamic_curves(&thinning), RADIUS_M, LENGTH_M).map(|_| ()),
        Err(MaterialStatePointError::OutsideConsumerDomain {
            property: DYNAMIC_VISCOSITY_PROPERTY.to_owned(),
            value: 0.0,
            admissibility: ScalarAdmissibility::StrictlyPositive,
        })
    );
}

#[test]
fn refusals_publish_nothing_and_the_kernel_keeps_the_body() {
    let component = dynamic();
    in_cx(&CancelGate::new_clock_free(), |cx| {
        let attempt = |dt: f64, inputs: &[f64], slots: usize| {
            let mut outputs = vec![f64::NAN; slots];
            let result = component.propose(&(), dt, inputs, &mut (), &mut outputs, cx);
            assert!(outputs.iter().all(|value| value.is_nan()), "{outputs:?}");
            result.map(|_| ())
        };
        let refused =
            |quantity: &'static str| Err(MaterialStatePointError::InvalidDerived { quantity });
        for kelvin in [299.0, 400.5] {
            assert_eq!(
                attempt(DT, &[1.0e4, 0.0, kelvin], 3),
                refused("laminar pipe temperature outside its declared range")
            );
        }
        for inputs in [[f64::NAN, 0.0, 320.0], [1.0e4, f64::INFINITY, 320.0]] {
            assert_eq!(
                attempt(DT, &inputs, 3),
                refused("laminar pipe needs three finite port values")
            );
        }
        for dt in [0.0, -1.0, f64::NAN] {
            assert_eq!(
                attempt(dt, &[1.0e4, 0.0, 320.0], 3),
                refused("laminar pipe timestep must be positive and finite")
            );
        }
        assert_eq!(
            attempt(DT, &[1.0e4, 0.0, 320.0], 2),
            refused("laminar pipe needs three output slots")
        );
        let mut slots = [f64::NAN; 2];
        assert_eq!(
            component.jvp(
                &(),
                DT,
                &[1.0e4, 0.0, 320.0],
                &[1.0, 0.0, 0.0],
                &mut slots,
                cx
            ),
            refused("laminar pipe needs three tangent output slots")
        );
        assert_eq!(
            component.vjp(
                &(),
                DT,
                &[1.0e4, 0.0, 320.0],
                &[1.0, 0.0, 0.0],
                &mut slots,
                cx
            ),
            refused("laminar pipe needs three tangent input slots")
        );
        assert!(slots.iter().all(|value| value.is_nan()));
        cx.cancel_gate().request();
        assert_eq!(
            attempt(DT, &[1.0e4, 0.0, 320.0], 3),
            refused("laminar pipe cancelled")
        );
    });
    // 3e5 Pa across is at least 35 W into 1 J/K: one 5 s step ends far past
    // 400 K, where the logarithmic mean the body's port publishes leaves the
    // card. No viscosity in the range closes that step, so the kernel
    // refuses it and names the pipe.
    let (mut system, _, body, _) = heated(3.0e5);
    let before = *system.state::<EquilibriumPhaseState>(body).unwrap();
    let error = in_cx(&CancelGate::new_clock_free(), |cx| system.step(5.0, cx)).unwrap_err();
    eprintln!("refused: {error}");
    assert!(
        matches!(&error, SystemError::Component { node, .. } if node == "pipe"),
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
fn the_identity_binds_the_card_the_geometry_and_the_basis() {
    let words = |component: &LaminarPipeComponent| {
        let mut words = Vec::new();
        component.identity_words(&mut words);
        words
    };
    let base = words(&dynamic());
    assert_eq!(base, words(&dynamic()));
    let other_card = [(300.0, 1.0e-3), (350.0, 5.1e-4), (400.0, 4.0e-4)];
    for other in [
        kinematic(),
        LaminarPipeComponent::new(pipe_of(&dynamic_curves(&MU), 1.1e-3, LENGTH_M).unwrap()),
        LaminarPipeComponent::new(pipe_of(&dynamic_curves(&MU), RADIUS_M, 2.0).unwrap()),
        LaminarPipeComponent::new(
            pipe_of(&dynamic_curves(&other_card), RADIUS_M, LENGTH_M).unwrap(),
        ),
    ] {
        assert_ne!(base, words(&other));
    }
    assert!(matches!(
        LaminarPipeComponent::STATE_CHARTER,
        Some(charter)
            if charter.owner == "fs-material::fluid_transport::laminar_pipe::LaminarPipeComponent"
    ));
    let component = dynamic();
    let mut encoded = Vec::new();
    assert!(component.encode_state(&(), &mut encoded));
    assert!(encoded.is_empty());
    assert_eq!(component.decode_state(&[]), Some(()));
    assert_eq!(component.decode_state(&[0]), None);
}
