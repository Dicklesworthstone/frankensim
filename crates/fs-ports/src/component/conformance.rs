//! Conformance checks every [`Component`] implementation runs.
//!
//! One call per component replaces the self-consistency tests each pairing
//! used to write for itself. A failure names the law that was broken.
//!
//! Covered here: port layout, step-limit bounds, purity of `propose`
//! (including independence from the stale candidate slot), that `propose`
//! refuses a step asked under a context cancelled before the call (any
//! error: whether it is the owner's cancellation error is not checked, nor
//! that `cx` is polled again later), the step audit (first law, endpoint
//! energy, non-negative dissipation and entropy production, residual),
//! available mass and entropy inventories, `admit`, and for
//! [`PortTangent`] the JVP against
//! central differences and the VJP against the JVP by the dot-product test,
//! and for a component that declares [`Class::ExactLinear`]
//! ([`check_exact_linear`]) that its slopes are the same bits at another
//! state and at other inputs.
//! For [`LoadLine`]: the line's ports, purity, the audit at the flows on the
//! line, and that a zero compliance gives `propose`. For a component a test
//! claims reversible, [`check_time_reversal`]: a step there and back, each
//! audited and dissipating nothing, whose outputs flip with their parity and
//! whose end returns to the start within bounds the test derives. For a
//! component a test builds in any units, [`check_unit_rescaling`]: the same
//! step in units changed by powers of two, which must give the same numbers
//! in the new units, to the bit.
//!
//! For [`StageResidual`]: [`check_stage_residual`] checks unknown-count
//! consistency, finite initial guess, cancellation polling of `residual`,
//! root satisfaction at the proposed solution, agreement of `accept` with
//! `propose` and `audit_step`, `jvp_z` against central finite differences of
//! both residual and outputs, and `vjp_z` against `jvp_z` by the dot-product
//! test.
//!
//! Temporal convergence, seeded shrinking, and zero-allocation gates:
//! [`check_temporal_order`] runs a refinement ladder through `fs-mms`
//! ([`RefinementLadder`], [`OrderGate`]); [`check_component_sweep`] runs a
//! deterministic `fs-propcheck` seeded sweep and shrinks any counterexample
//! with [`fs_propcheck::minimize_with_budget`]; [`check_zero_alloc_propose`]
//! verifies zero heap allocations after warm-up via an injected counting
//! allocator hook.

use core::fmt;

use fs_exec::{CancelGate, Cx};
pub use fs_mms::{LadderSide, OrderFit, OrderGate, OrderVerdict, RefinementLadder};
pub use fs_propcheck::{MinimizeBudget, Shrink, Stream};
use fs_qty::{DIMENSION_COUNT, Dims, Energy, Mass, Power, Temperature, Time};

use super::{
    Carrier, Class, Component, ContractViolation, LoadLine, PortLayout, PortRole, PortSpec,
    PortTangent, Proposal, Receipt, StageResidual, StepAudit, audit_carried_step,
    audit_entropy_step, audit_mass_step, audit_step,
};
use crate::Entropy;

/// One operating point at which a component is checked.
#[derive(Debug, Clone, Copy)]
pub struct Case<'a, S> {
    /// Accepted state the step starts from.
    pub state: &'a S,
    /// Any other valid state, used to fill the candidate slot with stale
    /// content; `propose` must give the same result whatever the slot holds.
    pub stale: &'a S,
    /// Step [s].
    pub dt_s: f64,
    /// Held inputs, in port-layout order.
    pub inputs: &'a [f64],
}

/// Why a component failed conformance.
#[derive(Debug)]
pub enum Failure<E> {
    /// The component refused the call with its own error.
    Owner(E),
    /// The step broke the contract's audit.
    Contract(ContractViolation),
    /// A law of the contract that the audit does not cover.
    Law {
        /// The law, stated as what must hold.
        law: &'static str,
        /// What was observed instead.
        detail: String,
    },
}

impl<E: fmt::Debug> fmt::Display for Failure<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Owner(error) => write!(f, "the component refused: {error:?}"),
            Self::Contract(violation) => write!(f, "contract violation: {violation}"),
            Self::Law { law, detail } => write!(f, "law broken: {law} ({detail})"),
        }
    }
}

/// What a passing check observed.
#[derive(Debug, Clone, PartialEq)]
pub enum Observed {
    /// The component took the step; the audit is attached.
    Stepped {
        /// The component's receipt.
        receipt: Receipt,
        /// Port energies and defects.
        audit: StepAudit,
        /// Outputs, in port-layout order.
        outputs: Vec<f64>,
    },
    /// The component asked for a shorter step, which the contract allows.
    Shortened {
        /// The step it offered [s].
        dt_s: f64,
    },
}

fn law<E>(law: &'static str, detail: String) -> Failure<E> {
    Failure::Law { law, detail }
}

/// Check the mandatory tier at one operating point.
///
/// # Errors
/// The first [`Failure`] found.
#[allow(clippy::too_many_lines)]
pub fn check_component<C>(
    component: &C,
    case: &Case<'_, C::State>,
    cx: &Cx<'_>,
) -> Result<Observed, Failure<C::Error>>
where
    C: Component,
    C::State: PartialEq + fmt::Debug,
{
    let ports = component.ports();
    let layout = PortLayout::new(ports).map_err(Failure::Contract)?;
    if case.inputs.len() != layout.input_len() {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "inputs",
            expected: layout.input_len(),
            actual: case.inputs.len(),
        }));
    }

    let limited = component
        .limit_dt_s(case.state, case.dt_s)
        .map_err(Failure::Owner)?;
    if !(limited > 0.0 && limited <= case.dt_s) {
        return Err(law(
            "limit_dt_s returns a positive step not above the proposal",
            format!("proposed {} s, got {limited} s", case.dt_s),
        ));
    }
    if let Some(event) = component.next_event_s(case.state)
        && !(event.is_finite() && event > 0.0)
    {
        return Err(law(
            "next_event_s is positive and finite when present",
            format!("got {event} s"),
        ));
    }

    let before = case.state.clone();
    // The contract has `propose` poll `cx` and return the owner's
    // cancellation error. Asked under a context built from `cx` on a gate
    // whose cancellation is requested before the call, it must refuse. Any
    // error passes: each owner's cancellation error is its own.
    let gate = CancelGate::new_clock_free();
    gate.request();
    let cancelled = Cx::new(&gate, cx.arena(), cx.stream_key(), cx.budget(), cx.mode());
    let mut cancelled_state = case.stale.clone();
    let mut cancelled_out = vec![f64::NAN; layout.output_len()];
    if let Ok(proposal) = component.propose(
        case.state,
        case.dt_s,
        case.inputs,
        &mut cancelled_state,
        &mut cancelled_out,
        &cancelled,
    ) {
        return Err(law(
            "propose polls cx: a step asked under a cancelled context is refused",
            format!("a context cancelled before the call gave {proposal:?}"),
        ));
    }
    let mut first_state = case.state.clone();
    let mut second_state = case.stale.clone();
    let mut first_out = vec![0.0; layout.output_len()];
    let mut second_out = vec![f64::NAN; layout.output_len()];
    let first = component
        .propose(
            case.state,
            case.dt_s,
            case.inputs,
            &mut first_state,
            &mut first_out,
            cx,
        )
        .map_err(Failure::Owner)?;
    let second = component
        .propose(
            case.state,
            case.dt_s,
            case.inputs,
            &mut second_state,
            &mut second_out,
            cx,
        )
        .map_err(Failure::Owner)?;
    if *case.state != before {
        return Err(law(
            "propose does not change the accepted state",
            format!("before {before:?}, after {:?}", case.state),
        ));
    }
    // The two calls differ only in what the candidate slot held, so a
    // difference in the candidates is reported as that before anything else.
    let receipt = match (first, second) {
        (Proposal::Step(receipt), Proposal::Step(_)) => receipt,
        (Proposal::Shorten { dt_s }, Proposal::Shorten { dt_s: again })
            if dt_s.to_bits() == again.to_bits() =>
        {
            if !(dt_s > 0.0 && dt_s < case.dt_s) {
                return Err(law(
                    "a shortened step is positive and below the proposal",
                    format!("proposed {} s, offered {dt_s} s", case.dt_s),
                ));
            }
            return Ok(Observed::Shortened { dt_s });
        }
        _ => {
            return Err(law(
                "propose is pure: equal arguments give equal proposals",
                format!("first {first:?}, second {second:?}"),
            ));
        }
    };
    if first_state != second_state {
        return Err(law(
            "propose reads nothing from the candidate slot",
            format!("clean slot gave {first_state:?}, stale slot gave {second_state:?}"),
        ));
    }
    let same_bits = first_out
        .iter()
        .zip(&second_out)
        .all(|(a, b)| a.to_bits() == b.to_bits());
    if !same_bits {
        return Err(law(
            "propose is pure: outputs are bitwise equal on repetition",
            format!("first {first_out:?}, second {second_out:?}"),
        ));
    }
    if first != second {
        return Err(law(
            "propose is pure: equal arguments give equal proposals",
            format!("first {first:?}, second {second:?}"),
        ));
    }

    let stored_old = component
        .stored_energy_j(case.state)
        .map_err(Failure::Owner)?;
    let stored_new = component
        .stored_energy_j(&first_state)
        .map_err(Failure::Owner)?;
    let audit = audit_step(
        ports,
        &layout,
        case.dt_s,
        case.inputs,
        &first_out,
        &receipt,
        stored_old,
        stored_new,
    )
    .map_err(Failure::Contract)?;
    let mass_old = component
        .stored_mass_kg(case.state)
        .map_err(Failure::Owner)?;
    let mass_new = component
        .stored_mass_kg(&first_state)
        .map_err(Failure::Owner)?;
    audit_mass_step(
        ports,
        &layout,
        case.dt_s,
        case.inputs,
        &first_out,
        mass_old,
        mass_new,
    )
    .map_err(Failure::Contract)?;
    let entropy_old = component
        .stored_entropy_j_per_k(case.state)
        .map_err(Failure::Owner)?;
    let entropy_new = component
        .stored_entropy_j_per_k(&first_state)
        .map_err(Failure::Owner)?;
    audit_entropy_step(
        ports,
        &layout,
        case.dt_s,
        case.inputs,
        &first_out,
        receipt.entropy_generated_j_per_k,
        entropy_old,
        entropy_new,
    )
    .map_err(Failure::Contract)?;
    for carried in 0..ports
        .iter()
        .map(PortSpec::carried_quantities)
        .max()
        .unwrap_or(0)
    {
        let old = component
            .stored_carried_quantity(case.state, carried)
            .map_err(Failure::Owner)?;
        let new = component
            .stored_carried_quantity(&first_state, carried)
            .map_err(Failure::Owner)?;
        audit_carried_step(
            ports,
            &layout,
            case.dt_s,
            case.inputs,
            &first_out,
            carried,
            old,
            new,
        )
        .map_err(Failure::Contract)?;
    }
    component
        .admit(case.state, &first_state, &receipt)
        .map_err(Failure::Owner)?;
    Ok(Observed::Stepped {
        receipt,
        audit,
        outputs: first_out,
    })
}

/// What a passing [`check_load_line`] observed.
#[derive(Debug, Clone, PartialEq)]
pub struct OnLine {
    /// The step, audited at the flows on the line.
    pub observed: Observed,
    /// Those flows, `released - compliance * effort`, in the order of the
    /// line's coefficients. Empty if the component asked to shorten.
    pub flows: Vec<f64>,
}

/// One call of a step function into a fresh copy of `slot`.
fn stepped<S, E, F>(
    slot: &S,
    len: usize,
    fill: f64,
    step: &F,
) -> Result<(Proposal, S, Vec<f64>), Failure<E>>
where
    S: Clone,
    F: Fn(&mut S, &mut [f64]) -> Result<Proposal, E>,
{
    let mut state = slot.clone();
    let mut outputs = vec![fill; len];
    let proposal = step(&mut state, &mut outputs).map_err(Failure::Owner)?;
    Ok((proposal, state, outputs))
}

/// Check the load-line tier at one operating point.
///
/// `case.inputs` holds the released flows at the line's ports, and
/// `compliance` is the line's matrix, as [`LoadLine::propose_on_line`] takes
/// them.
///
/// Covered: the line's ports are power ports that receive a flow, listed in
/// increasing order; the size of `compliance`; purity of `propose_on_line`
/// (including independence from the stale candidate slot); the step audit
/// at the flows on the line; `admit`; and that with a zero compliance the
/// call is `propose`, to the bit. Not covered: that the efforts satisfy the
/// component's own relation, which only the component can state.
///
/// # Errors
/// The first [`Failure`] found.
#[allow(clippy::too_many_lines)]
pub fn check_load_line<C>(
    component: &C,
    case: &Case<'_, C::State>,
    compliance: &[f64],
    cx: &Cx<'_>,
) -> Result<OnLine, Failure<C::Error>>
where
    C: LoadLine,
    C::State: PartialEq + fmt::Debug,
{
    let ports = component.ports();
    let layout = PortLayout::new(ports).map_err(Failure::Contract)?;
    if case.inputs.len() != layout.input_len() {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "inputs",
            expected: layout.input_len(),
            actual: case.inputs.len(),
        }));
    }

    // Per coefficient of the line: where its flow is received and where its
    // effort is returned.
    let line = component.line_ports();
    if line.is_empty() {
        return Err(law(
            "the line spans at least one port",
            "it lists none".to_owned(),
        ));
    }
    let mut places: Vec<(usize, usize)> = Vec::new();
    let mut last = None;
    for &port in line {
        if port >= ports.len() || last.is_some_and(|before| port <= before) {
            return Err(law(
                "the line's ports are indices of the component's ports, in increasing order",
                format!("{line:?} with {} ports", ports.len()),
            ));
        }
        let spec = &ports[port];
        if !matches!(spec.carrier(), Carrier::Power(_)) || spec.role() != PortRole::FlowIn {
            return Err(law(
                "the line's ports are power ports that receive a flow",
                format!(
                    "port {port} is {:?} with role {:?}",
                    spec.carrier(),
                    spec.role()
                ),
            ));
        }
        places.extend(layout.input_range(port).zip(layout.output_range(port)));
        last = Some(port);
    }
    let width = places.len();
    // A component that names end inputs is handed one more row per
    // coefficient of them.
    let end_rows: usize = component
        .line_ends()
        .iter()
        .flatten()
        .map(|&port| layout.input_range(port).len())
        .sum();
    if compliance.len() != (width + end_rows) * width {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "compliance",
            expected: (width + end_rows) * width,
            actual: compliance.len(),
        }));
    }

    let before = case.state.clone();
    let given = |state: &mut C::State, outputs: &mut [f64]| {
        component.propose_on_line(
            case.state,
            case.dt_s,
            case.inputs,
            compliance,
            state,
            outputs,
            cx,
        )
    };
    let (first, first_state, first_out) = stepped(case.state, layout.output_len(), 0.0, &given)?;
    let (second, second_state, second_out) =
        stepped(case.stale, layout.output_len(), f64::NAN, &given)?;
    if *case.state != before {
        return Err(law(
            "propose_on_line does not change the accepted state",
            format!("before {before:?}, after {:?}", case.state),
        ));
    }
    let receipt = match (first, second) {
        (Proposal::Step(receipt), Proposal::Step(_)) => receipt,
        (Proposal::Shorten { dt_s }, Proposal::Shorten { dt_s: again })
            if dt_s.to_bits() == again.to_bits() =>
        {
            if !(dt_s > 0.0 && dt_s < case.dt_s) {
                return Err(law(
                    "a shortened step is positive and below the proposal",
                    format!("proposed {} s, offered {dt_s} s", case.dt_s),
                ));
            }
            return Ok(OnLine {
                observed: Observed::Shortened { dt_s },
                flows: Vec::new(),
            });
        }
        _ => {
            return Err(law(
                "propose_on_line is pure: equal arguments give equal proposals",
                format!("first {first:?}, second {second:?}"),
            ));
        }
    };
    if first_state != second_state {
        return Err(law(
            "propose_on_line reads nothing from the candidate slot",
            format!("clean slot gave {first_state:?}, stale slot gave {second_state:?}"),
        ));
    }
    let same_bits = |a: &[f64], b: &[f64]| a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits());
    if !same_bits(&first_out, &second_out) || first != second {
        return Err(law(
            "propose_on_line is pure: outputs and proposals are equal on repetition",
            format!("first {first:?} {first_out:?}, second {second:?} {second_out:?}"),
        ));
    }

    // The step is audited where it was taken: at the flows on the line.
    let mut at_line = case.inputs.to_vec();
    let mut flows = Vec::with_capacity(width);
    for (row, &(input, _)) in places.iter().enumerate() {
        let mut flow = case.inputs[input];
        for (column, &(_, output)) in places.iter().enumerate() {
            flow -= compliance[row * width + column] * first_out[output];
        }
        at_line[input] = flow;
        flows.push(flow);
    }
    let stored_old = component
        .stored_energy_j(case.state)
        .map_err(Failure::Owner)?;
    let stored_new = component
        .stored_energy_j(&first_state)
        .map_err(Failure::Owner)?;
    let audit = audit_step(
        ports, &layout, case.dt_s, &at_line, &first_out, &receipt, stored_old, stored_new,
    )
    .map_err(Failure::Contract)?;
    component
        .admit(case.state, &first_state, &receipt)
        .map_err(Failure::Owner)?;

    // With nothing taken away by an effort the line is the released flows
    // themselves, and the call is the mandatory tier's.
    let zero = vec![0.0; (width + end_rows) * width];
    let unloaded = |state: &mut C::State, outputs: &mut [f64]| {
        component.propose_on_line(
            case.state,
            case.dt_s,
            case.inputs,
            &zero,
            state,
            outputs,
            cx,
        )
    };
    let (flat, flat_state, flat_out) = stepped(case.state, layout.output_len(), 0.0, &unloaded)?;
    let plain = |state: &mut C::State, outputs: &mut [f64]| {
        component.propose(case.state, case.dt_s, case.inputs, state, outputs, cx)
    };
    let (held, held_state, held_out) = stepped(case.state, layout.output_len(), 0.0, &plain)?;
    if flat != held || flat_state != held_state || !same_bits(&flat_out, &held_out) {
        return Err(law(
            "with a zero compliance propose_on_line is propose, to the bit",
            format!(
                "on the line {flat:?} {flat_state:?} {flat_out:?}, held {held:?} {held_state:?} \
                 {held_out:?}"
            ),
        ));
    }

    Ok(OnLine {
        observed: Observed::Stepped {
            receipt,
            audit,
            outputs: first_out,
        },
        flows,
    })
}

/// Deterministic direction with entries in (-1, 1), different for each seed.
fn direction(len: usize, seed: u64) -> Vec<f64> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1_u64 << 53) as f64 * 2.0 - 1.0
        })
        .collect()
}

fn outputs_at<C>(
    component: &C,
    case: &Case<'_, C::State>,
    inputs: &[f64],
    len: usize,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, Failure<C::Error>>
where
    C: Component,
{
    let mut candidate = case.state.clone();
    let mut outputs = vec![0.0; len];
    match component
        .propose(
            case.state,
            case.dt_s,
            inputs,
            &mut candidate,
            &mut outputs,
            cx,
        )
        .map_err(Failure::Owner)?
    {
        Proposal::Step(_) => Ok(outputs),
        Proposal::Shorten { dt_s } => Err(law(
            "the tangent is checked where the component takes the step",
            format!("the component asked to shorten to {dt_s} s"),
        )),
    }
}

/// Check the tangent tier at one operating point.
///
/// The JVP is compared with a central difference of `propose` along three
/// fixed directions; `relative_tolerance` applies to the larger of the two
/// norms plus `absolute_tolerance`. The VJP is compared with the JVP through
/// `<bar, J d> = <J^T bar, d>`.
///
/// # Errors
/// The first [`Failure`] found.
pub fn check_port_tangent<C>(
    component: &C,
    case: &Case<'_, C::State>,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<(), Failure<C::Error>>
where
    C: PortTangent,
{
    let layout = PortLayout::new(component.ports()).map_err(Failure::Contract)?;
    let (n_in, n_out) = (layout.input_len(), layout.output_len());
    if case.inputs.len() != n_in {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "inputs",
            expected: n_in,
            actual: case.inputs.len(),
        }));
    }
    let scale = case.inputs.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
    let step = scale * f64::EPSILON.cbrt();
    for seed in 1..=3_u64 {
        let d = direction(n_in, seed);
        let plus: Vec<f64> = case
            .inputs
            .iter()
            .zip(&d)
            .map(|(x, d)| x + step * d)
            .collect();
        let minus: Vec<f64> = case
            .inputs
            .iter()
            .zip(&d)
            .map(|(x, d)| x - step * d)
            .collect();
        let f_plus = outputs_at(component, case, &plus, n_out, cx)?;
        let f_minus = outputs_at(component, case, &minus, n_out, cx)?;
        let mut jvp = vec![0.0; n_out];
        component
            .jvp(case.state, case.dt_s, case.inputs, &d, &mut jvp, cx)
            .map_err(Failure::Owner)?;
        let (mut error, mut size) = (0.0_f64, 0.0_f64);
        for ((p, m), j) in f_plus.iter().zip(&f_minus).zip(&jvp) {
            let fd = (p - m) / (2.0 * step);
            error = error.max((fd - j).abs());
            size = size.max(fd.abs()).max(j.abs());
        }
        if !(error <= relative_tolerance * size + absolute_tolerance) {
            return Err(law(
                "jvp equals the derivative of the port map",
                format!("direction {seed}: max difference {error}, derivative size {size}"),
            ));
        }
        let bar = direction(n_out, seed + 100);
        let mut vjp = vec![0.0; n_in];
        component
            .vjp(case.state, case.dt_s, case.inputs, &bar, &mut vjp, cx)
            .map_err(Failure::Owner)?;
        let left: f64 = bar.iter().zip(&jvp).map(|(a, b)| a * b).sum();
        let right: f64 = vjp.iter().zip(&d).map(|(a, b)| a * b).sum();
        let magnitude: f64 = bar
            .iter()
            .zip(&jvp)
            .map(|(a, b)| (a * b).abs())
            .sum::<f64>()
            + vjp.iter().zip(&d).map(|(a, b)| (a * b).abs()).sum::<f64>();
        if !((left - right).abs() <= 64.0 * f64::EPSILON * magnitude + f64::MIN_POSITIVE) {
            return Err(law(
                "vjp is the transpose of jvp",
                format!("<bar, J d> = {left}, <J^T bar, d> = {right}"),
            ));
        }
    }
    if component.class() == Class::ExactLinear {
        // The class is a claim about every state and every input. It is
        // held here to the two this case can name without the test's help:
        // its own point, and its other valid state at inputs moved off it.
        let moved: Vec<f64> = case.inputs.iter().map(|x| 1.5 * x + 0.25).collect();
        let other = Case {
            state: case.stale,
            stale: case.state,
            dt_s: case.dt_s,
            inputs: &moved,
        };
        check_exact_linear(component, case, &other, cx)?;
    }
    Ok(())
}

/// Check the claim of [`Class::ExactLinear`] between two operating points.
///
/// A component of that class says that the slopes of its port map depend on
/// the step length alone. A solver may then keep what it derived from them
/// at one step and use it at a later step of the same length, and a run
/// that kept it must take the steps of a run that derived it again. So the
/// slopes must be the same numbers at any two states and any two inputs, to
/// the bit, a zero's sign included. The JVP along every unit input
/// direction is taken at `case` and at `other`, which share `dt_s`, and
/// compared bit for bit. The caller chooses two points that differ;
/// [`check_port_tangent`] asks this of every component of the class, at its
/// case and at the case's stale state with other inputs.
///
/// # Errors
/// The first [`Failure`] found: a component of another class, two steps, an
/// input vector of the wrong length, or a slope that differs, named by its
/// output and its input with both values and their bits.
pub fn check_exact_linear<C>(
    component: &C,
    case: &Case<'_, C::State>,
    other: &Case<'_, C::State>,
    cx: &Cx<'_>,
) -> Result<(), Failure<C::Error>>
where
    C: PortTangent,
{
    let class = component.class();
    if class != Class::ExactLinear {
        return Err(law(
            "the exact-linear law is asked of a component of that class",
            format!("the component declares {class:?}"),
        ));
    }
    if case.dt_s.to_bits() != other.dt_s.to_bits() {
        return Err(law(
            "the two points of the exact-linear law share a step",
            format!("{} s and {} s", case.dt_s, other.dt_s),
        ));
    }
    let layout = PortLayout::new(component.ports()).map_err(Failure::Contract)?;
    let (n_in, n_out) = (layout.input_len(), layout.output_len());
    for point in [case, other] {
        if point.inputs.len() != n_in {
            return Err(Failure::Contract(ContractViolation::VectorLength {
                vector: "inputs",
                expected: n_in,
                actual: point.inputs.len(),
            }));
        }
    }
    let mut unit = vec![0.0; n_in];
    let (mut here, mut there) = (vec![0.0; n_out], vec![0.0; n_out]);
    for input in 0..n_in {
        unit[input] = 1.0;
        component
            .jvp(case.state, case.dt_s, case.inputs, &unit, &mut here, cx)
            .map_err(Failure::Owner)?;
        component
            .jvp(other.state, other.dt_s, other.inputs, &unit, &mut there, cx)
            .map_err(Failure::Owner)?;
        unit[input] = 0.0;
        for (output, (a, b)) in here.iter().zip(&there).enumerate() {
            if a.to_bits() != b.to_bits() {
                return Err(law(
                    "the slopes of an exact-linear component depend on the step alone",
                    format!(
                        "d output {output} / d input {input} is {a:e} ({:#018x}) at the first \
                         point and {b:e} ({:#018x}) at the second",
                        a.to_bits(),
                        b.to_bits()
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// Absolute bounds for [`check_time_reversal`], one per compared number.
///
/// The test derives them from the rounding of the component's own
/// arithmetic, evaluated at the states of the round trip; the kit holds no
/// tolerance of its own.
#[derive(Debug, Clone, PartialEq)]
pub struct ReversalAllowance {
    /// One bound per output, in port-layout order.
    pub outputs: Vec<f64>,
    /// One bound per number the test's `observe` returns.
    pub observed: Vec<f64>,
    /// Bound on the stored energy's change over the round trip [J].
    pub energy_j: f64,
}

/// The test's bounds for [`check_time_reversal`]: from the start `x0`, the
/// end `x1` of the forward step and the end `x2` of the backward step, in
/// that order, the [`ReversalAllowance`] of that round trip. The lifetime
/// lets the closure borrow, as a `&dyn Fn` written in place would.
pub type ReversalBounds<'a, S> = dyn Fn(&S, &S, &S) -> ReversalAllowance + 'a;

/// What a passing [`check_time_reversal`] measured, beside what it allowed.
#[derive(Debug, Clone, PartialEq)]
pub struct Reversal {
    /// Per output: the backward output less the forward one times the
    /// output's parity, in magnitude.
    pub output_gaps: Vec<f64>,
    /// Per observed number: `R` of the backward end against the start.
    pub observed_gaps: Vec<f64>,
    /// Stored energy of `R` of the backward end against the start's [J].
    pub energy_gap_j: f64,
    /// The bounds the gaps were held to.
    pub allowance: ReversalAllowance,
}

/// Parity under time reversal of every input and output coefficient: `1.0`
/// for an even one, `-1.0` for an odd one. Effort and temperature are even,
/// flow and heat rate odd; a port that is neither power nor heat, or that
/// pairs a field through a metric, has no parity here.
fn parities<E>(
    ports: &[PortSpec],
    layout: &PortLayout,
) -> Result<(Vec<f64>, Vec<f64>), Failure<E>> {
    let mut inputs = vec![1.0; layout.input_len()];
    let mut outputs = vec![1.0; layout.output_len()];
    for (index, port) in ports.iter().enumerate() {
        let carried = matches!(port.carrier(), Carrier::Power(_) | Carrier::Heat)
            && port.field_metric().is_none();
        let (input, output) = match port.role() {
            PortRole::EffortIn if carried => (1.0, -1.0),
            PortRole::FlowIn if carried => (-1.0, 1.0),
            role => {
                let metric = port
                    .field_metric()
                    .map_or("", |_| " through a field metric");
                return Err(law(
                    "time reversal is checked on power and heat ports only",
                    format!(
                        "port `{}` carries {:?} as {role:?}{metric}",
                        port.id().as_str(),
                        port.carrier(),
                    ),
                ));
            }
        };
        inputs[layout.input_range(index)].fill(input);
        outputs[layout.output_range(index)].fill(output);
    }
    Ok((inputs, outputs))
}

/// Where a leg of the round trip ends, and the outputs of its step.
type LegEnd<S> = (S, Vec<f64>);

/// One leg of the round trip: the mandatory tier's check at `case`, a step
/// taken that dissipates nothing, and the candidate and outputs of that step.
fn reversible_leg<C>(
    component: &C,
    case: &Case<'_, C::State>,
    leg: &str,
    cx: &Cx<'_>,
) -> Result<LegEnd<C::State>, Failure<C::Error>>
where
    C: Component,
    C::State: PartialEq + fmt::Debug,
{
    let receipt = match check_component(component, case, cx)? {
        Observed::Stepped { receipt, .. } => receipt,
        Observed::Shortened { dt_s } => {
            return Err(law(
                "time reversal is checked where the component takes the step",
                format!("the {leg} step asked to shorten to {dt_s} s"),
            ));
        }
    };
    if receipt.dissipated_j != 0.0 || receipt.entropy_generated_j_per_k != 0.0 {
        return Err(law(
            "a reversible step dissipates nothing and generates no entropy",
            format!(
                "the {leg} step dissipated {} J and generated {} J/K",
                receipt.dissipated_j, receipt.entropy_generated_j_per_k
            ),
        ));
    }
    // check_component has shown that propose is pure, so this third call
    // returns what it audited.
    let layout = PortLayout::new(component.ports()).map_err(Failure::Contract)?;
    let plain = |state: &mut C::State, outputs: &mut [f64]| {
        component.propose(case.state, case.dt_s, case.inputs, state, outputs, cx)
    };
    let (_, candidate, outputs) = stepped(case.state, layout.output_len(), 0.0, &plain)?;
    Ok((candidate, outputs))
}

/// Each gap against its bound, the counts first.
fn within<E>(
    law_text: &'static str,
    what: &str,
    gaps: &[f64],
    bounds: &[f64],
) -> Result<(), Failure<E>> {
    if gaps.len() != bounds.len() {
        return Err(law(
            "the allowance gives one bound per compared number",
            format!("{} {what} numbers, {} bounds", gaps.len(), bounds.len()),
        ));
    }
    for (index, (gap, bound)) in gaps.iter().zip(bounds).enumerate() {
        if !(gap <= bound) {
            return Err(law(
                law_text,
                format!("{what} {index}: gap {gap:e}, allowed {bound:e}"),
            ));
        }
    }
    Ok(())
}

/// Check time reversal at one operating point, for a component the test
/// claims reversible.
///
/// `reverse` is the test's `R`: the state with its odd parts (velocities,
/// momenta, fluxes) negated and its even parts kept. From the start `x0` the
/// component steps to `x1` at the case's inputs `u`; from `R(x1)` it steps to
/// `x2` at `P u`, the inputs with their odd coefficients negated. The parity
/// `P` comes from the port roles, not from the test: effort and temperature
/// are even, flow and heat rate odd, and a component with another kind of
/// port is refused. The laws, in order:
/// 1. `R` is an involution: `R(R(x0)) == x0`.
/// 2. Every port is a power or a heat port without a field metric.
/// 3. The forward step passes [`check_component`], is taken, and its receipt
///    dissipates nothing and generates no entropy, exactly.
/// 4. So does the backward step.
/// 5. Each backward output is the forward one times its parity, within its
///    bound.
/// 6. `R(x2)` returns to `x0` as `observe` sees a state, within a bound per
///    number, and holds `x0`'s stored energy, within `energy_j`.
///
/// `allowance` receives `x0`, `x1` and `x2` and returns absolute bounds the
/// test derives from the rounding of the component's arithmetic.
///
/// # Errors
/// The first [`Failure`] found.
pub fn check_time_reversal<C>(
    component: &C,
    case: &Case<'_, C::State>,
    reverse: &dyn Fn(&C::State) -> C::State,
    observe: &dyn Fn(&C::State) -> Vec<f64>,
    allowance: &ReversalBounds<'_, C::State>,
    cx: &Cx<'_>,
) -> Result<Reversal, Failure<C::Error>>
where
    C: Component,
    C::State: PartialEq + fmt::Debug,
{
    let start = case.state;
    let twice = reverse(&reverse(start));
    if twice != *start {
        return Err(law(
            "time reversal is an involution: R(R(x)) == x",
            format!("x {start:?}, R(R(x)) {twice:?}"),
        ));
    }
    let ports = component.ports();
    let layout = PortLayout::new(ports).map_err(Failure::Contract)?;
    let (input_parity, output_parity) = parities(ports, &layout)?;

    let (forward_end, forward_outputs) = reversible_leg(component, case, "forward", cx)?;
    let back_start = reverse(&forward_end);
    let back_inputs: Vec<f64> = case
        .inputs
        .iter()
        .zip(&input_parity)
        .map(|(input, parity)| parity * input)
        .collect();
    let back_case = Case {
        state: &back_start,
        stale: case.stale,
        dt_s: case.dt_s,
        inputs: &back_inputs,
    };
    let (back_end, back_outputs) = reversible_leg(component, &back_case, "backward", cx)?;
    let bounds = allowance(start, &forward_end, &back_end);

    let output_gaps: Vec<f64> = back_outputs
        .iter()
        .zip(&forward_outputs)
        .zip(&output_parity)
        .map(|((back, forward), parity)| (back - parity * forward).abs())
        .collect();
    within(
        "the backward outputs are the forward outputs times their parity",
        "output",
        &output_gaps,
        &bounds.outputs,
    )?;

    let returned = reverse(&back_end);
    let (seen, expected) = (observe(&returned), observe(start));
    if seen.len() != expected.len() {
        return Err(law(
            "observe gives as many numbers for every state",
            format!("{} for R(x2), {} for x0", seen.len(), expected.len()),
        ));
    }
    let observed_gaps: Vec<f64> = seen
        .iter()
        .zip(&expected)
        .map(|(seen, expected)| (seen - expected).abs())
        .collect();
    within(
        "R of the backward end returns to the start",
        "observed",
        &observed_gaps,
        &bounds.observed,
    )?;
    let returned_j = component
        .stored_energy_j(&returned)
        .map_err(Failure::Owner)?;
    let start_j = component.stored_energy_j(start).map_err(Failure::Owner)?;
    let energy_gap_j = (returned_j - start_j).abs();
    if !(energy_gap_j <= bounds.energy_j) {
        return Err(law(
            "R of the backward end holds the start's stored energy",
            format!("gap {energy_gap_j:e} J, allowed {:e} J", bounds.energy_j),
        ));
    }
    Ok(Reversal {
        output_gaps,
        observed_gaps,
        energy_gap_j,
        allowance: bounds,
    })
}

/// Base-2 exponents of a change of units, one per SI base dimension in the
/// order of [`Dims`], `[m, kg, s, K, A, mol]`. In the new units a number of
/// dimensions `d` is `2^(d . exponents)` times the old one: exponent `3` on
/// metres measures a length in eighths of a metre.
pub type UnitExponents = [i32; DIMENSION_COUNT];

/// The largest exponent [`check_unit_rescaling`] takes, in size.
pub const UNIT_EXPONENT_LIMIT: i32 = 8;

/// `2^(d . exponents)`: the factor by which a number of dimensions `d`
/// changes in the units `exponents` name. Exact while it is a normal `f64`;
/// infinity or zero beyond that.
#[must_use]
pub fn unit_factor(dimensions: Dims, exponents: &UnitExponents) -> f64 {
    let power: i64 = dimensions
        .0
        .iter()
        .zip(exponents)
        .map(|(dimension, exponent)| i64::from(*dimension) * i64::from(*exponent))
        .sum();
    match u64::try_from(power + 1023) {
        // The biased exponent of a normal f64, over a zero fraction.
        Ok(biased @ 1..=2046) => f64::from_bits(biased << 52),
        _ if power > 0 => f64::INFINITY,
        _ => 0.0,
    }
}

const RANGE_LAW: &str = "the rescaled numbers stay in the normal range of f64";
const STEP_LAW: &str = "the rescaled step is taken or shortened as the base step is";
const START_LAW: &str =
    "the rescaled start's energy, mass, entropy and step limits are the base ones in the new units";

/// Zero or a normal `f64`: a number a power of two scales exactly.
fn scalable(value: f64) -> bool {
    value == 0.0 || value.is_normal()
}

/// Dimensions of every input and output coefficient, from the ports alone:
/// a power port's from its kind, a heat port's kelvin and watts.
fn coefficient_dimensions<E>(
    ports: &[PortSpec],
    layout: &PortLayout,
) -> Result<(Vec<Dims>, Vec<Dims>), Failure<E>> {
    let mut inputs = vec![Dims::NONE; layout.input_len()];
    let mut outputs = vec![Dims::NONE; layout.output_len()];
    for (index, port) in ports.iter().enumerate() {
        let pair = match port.carrier() {
            Carrier::Power(kind) if port.field_metric().is_none() => Some((
                kind.canonical_effort_dimensions(),
                kind.canonical_flow_dimensions(),
            )),
            Carrier::Heat => Some((Temperature::DIMS, Power::DIMS)),
            _ => None,
        };
        let (input, output) = match (pair, port.role()) {
            (Some((effort, flow)), PortRole::EffortIn) => (effort, flow),
            (Some((effort, flow)), PortRole::FlowIn) => (flow, effort),
            (_, role) => {
                let metric = port
                    .field_metric()
                    .map_or("", |_| " through a field metric");
                return Err(law(
                    "unit rescaling is checked on power and heat ports only",
                    format!(
                        "port `{}` carries {:?} as {role:?}{metric}",
                        port.id().as_str(),
                        port.carrier(),
                    ),
                ));
            }
        };
        inputs[layout.input_range(index)].fill(input);
        outputs[layout.output_range(index)].fill(output);
    }
    Ok((inputs, outputs))
}

/// Each rescaled number against the base number times its factor, to the
/// bit: `(name, rescaled, base, factor)`.
fn in_new_units<E>(
    law_text: &'static str,
    numbers: impl IntoIterator<Item = (String, f64, f64, f64)>,
) -> Result<(), Failure<E>> {
    for (name, rescaled, base, factor) in numbers {
        let expected = base * factor;
        if !scalable(expected) {
            return Err(law(
                RANGE_LAW,
                format!("{name}: {base:e} x {factor:e} = {expected:e}"),
            ));
        }
        if rescaled.to_bits() != expected.to_bits() {
            return Err(law(
                law_text,
                format!("{name}: {rescaled:e}, not {base:e} x {factor:e} = {expected:e}"),
            ));
        }
    }
    Ok(())
}

/// The candidate and outputs of one more call: `check_component` has shown
/// that `propose` is pure, so this is the step it checked.
fn end_of_step<C>(
    component: &C,
    case: &Case<'_, C::State>,
    outputs: usize,
    cx: &Cx<'_>,
) -> Result<LegEnd<C::State>, Failure<C::Error>>
where
    C: Component,
{
    let plain = |state: &mut C::State, out: &mut [f64]| {
        component.propose(case.state, case.dt_s, case.inputs, state, out, cx)
    };
    let (_, candidate, outputs) = stepped(case.state, outputs, 0.0, &plain)?;
    Ok((candidate, outputs))
}

/// Check, at one operating point, that a component keeps its law in other
/// units.
///
/// `build` makes the component in the units `exponents` name, each parameter
/// times the factor of its own dimensions ([`unit_factor`]); the base
/// component is `build` of all-zero exponents. `scale_state` gives a state in
/// those units. The kit scales the step by the factor of seconds and each
/// input and output by the factor of its port's dimensions: a power port's
/// from its kind, a heat port's kelvin and watts. Every factor is a power of
/// two, so the change of units is exact in binary, and the laws hold to the
/// bit, with no tolerance. In order:
/// 1. Every exponent is at most [`UNIT_EXPONENT_LIMIT`] in size.
/// 2. Every port is a power port without a field metric, or a heat port.
/// 3. The rescaled component has the base component's ports.
/// 4. The base case passes [`check_component`]; so does the rescaled case:
///    the state and the stale slot through `scale_state`, the step and the
///    inputs times their factors.
/// 5. The rescaled step is taken where the base step is taken, and shortened
///    where it is shortened, to the base step times the factor of seconds.
/// 6. The rescaled outputs are the base outputs times their factors.
/// 7. The rescaled step ends at `scale_state` of the base end.
/// 8. The rescaled receipt's stored change, dissipation, numerical loss and
///    energy tolerance are the base ones times the factor of joules, and its
///    entropy generated the base one times the factor of J/K.
/// 9. At the start, the stored energy and mass, `limit_dt_s` and
///    `next_event_s` are the base ones times the factors of joules,
///    kilograms and seconds.
///
/// Every number scaled or compared stays zero or normal, or the case is
/// refused as out of range, before any law about the component.
///
/// A failure is a number in the component's arithmetic that does not carry
/// units: a constant with dimensions (an absolute tolerance, a floor in
/// metres or kelvin, a default time), a parameter `build` left in the old
/// units, or a function other than a power applied to a dimensional value
/// (`exp(T / 300)`, `ln(p)`). The receipt's `residual` and `tolerance` are
/// in the owner's units and are not compared.
///
/// # Errors
/// The first [`Failure`] found.
#[allow(clippy::too_many_lines)]
pub fn check_unit_rescaling<C>(
    build: &dyn Fn(&UnitExponents) -> C,
    case: &Case<'_, C::State>,
    scale_state: &dyn Fn(&C::State, &UnitExponents) -> C::State,
    exponents: &UnitExponents,
    cx: &Cx<'_>,
) -> Result<Observed, Failure<C::Error>>
where
    C: Component,
    C::State: PartialEq + fmt::Debug,
{
    if exponents
        .iter()
        .any(|exponent| exponent.abs() > UNIT_EXPONENT_LIMIT)
    {
        return Err(law(
            "unit exponents are at most UNIT_EXPONENT_LIMIT in size",
            format!("{exponents:?}"),
        ));
    }
    let base = build(&[0; DIMENSION_COUNT]);
    let rescaled = build(exponents);
    let ports = base.ports();
    let layout = PortLayout::new(ports).map_err(Failure::Contract)?;
    let (input_dimensions, output_dimensions) = coefficient_dimensions(ports, &layout)?;
    if rescaled.ports() != ports {
        return Err(law(
            "the rescaled component has the base component's ports",
            format!("base {ports:?}, rescaled {:?}", rescaled.ports()),
        ));
    }
    let base_seen = check_component(&base, case, cx)?;

    let factor = |dimensions: Dims| unit_factor(dimensions, exponents);
    let (seconds, joules) = (factor(Time::DIMS), factor(Energy::DIMS));
    let dt_s = case.dt_s * seconds;
    let inputs: Vec<f64> = case
        .inputs
        .iter()
        .zip(&input_dimensions)
        .map(|(input, dimensions)| input * factor(*dimensions))
        .collect();
    if let Some(value) = inputs
        .iter()
        .chain([&dt_s])
        .find(|value| !scalable(**value))
    {
        return Err(law(
            RANGE_LAW,
            format!("step {dt_s:e} s, inputs {inputs:?}: {value:e}"),
        ));
    }
    let start = scale_state(case.state, exponents);
    let stale = scale_state(case.stale, exponents);
    let rescaled_case = Case {
        state: &start,
        stale: &stale,
        dt_s,
        inputs: &inputs,
    };
    let seen = check_component(&rescaled, &rescaled_case, cx)?;

    match (&base_seen, &seen) {
        (
            Observed::Stepped {
                receipt: before, ..
            },
            Observed::Stepped { receipt: after, .. },
        ) => {
            let (base_end, base_outputs) = end_of_step(&base, case, layout.output_len(), cx)?;
            let (end, outputs) = end_of_step(&rescaled, &rescaled_case, layout.output_len(), cx)?;
            in_new_units(
                "the rescaled outputs are the base outputs in the new units",
                (0..outputs.len()).map(|index| {
                    let scale = factor(output_dimensions[index]);
                    let name = format!("output {index}");
                    (name, outputs[index], base_outputs[index], scale)
                }),
            )?;
            let expected = scale_state(&base_end, exponents);
            if end != expected {
                return Err(law(
                    "the rescaled step ends at the base end in the new units",
                    format!("{end:?}, not {expected:?}"),
                ));
            }
            let entropy = factor(Entropy::DIMS);
            let receipt = [
                (
                    "stored_change_j",
                    after.stored_change_j,
                    before.stored_change_j,
                    joules,
                ),
                (
                    "dissipated_j",
                    after.dissipated_j,
                    before.dissipated_j,
                    joules,
                ),
                (
                    "numerical_loss_j",
                    after.numerical_loss_j,
                    before.numerical_loss_j,
                    joules,
                ),
                (
                    "energy_tolerance_j",
                    after.energy_tolerance_j,
                    before.energy_tolerance_j,
                    joules,
                ),
                (
                    "entropy_generated_j_per_k",
                    after.entropy_generated_j_per_k,
                    before.entropy_generated_j_per_k,
                    entropy,
                ),
            ];
            in_new_units(
                "the rescaled receipt is the base receipt in the new units",
                receipt.map(|(name, after, before, scale)| (name.to_owned(), after, before, scale)),
            )?;
        }
        (Observed::Shortened { dt_s: before }, Observed::Shortened { dt_s: after }) => {
            in_new_units(
                STEP_LAW,
                [("shortened step".to_owned(), *after, *before, seconds)],
            )?;
        }
        _ => {
            return Err(law(
                STEP_LAW,
                format!("base {base_seen:?}, rescaled {seen:?}"),
            ));
        }
    }

    let owner = Failure::Owner;
    let mut numbers = vec![
        (
            "stored_energy_j".to_owned(),
            rescaled.stored_energy_j(&start).map_err(owner)?,
            base.stored_energy_j(case.state).map_err(owner)?,
            joules,
        ),
        (
            "limit_dt_s".to_owned(),
            rescaled.limit_dt_s(&start, dt_s).map_err(owner)?,
            base.limit_dt_s(case.state, case.dt_s).map_err(owner)?,
            seconds,
        ),
    ];
    let optional = [
        (
            "stored_mass_kg",
            rescaled.stored_mass_kg(&start).map_err(owner)?,
            base.stored_mass_kg(case.state).map_err(owner)?,
            factor(Mass::DIMS),
        ),
        (
            "stored_entropy_j_per_k",
            rescaled.stored_entropy_j_per_k(&start).map_err(owner)?,
            base.stored_entropy_j_per_k(case.state).map_err(owner)?,
            factor(Entropy::DIMS),
        ),
        (
            "next_event_s",
            rescaled.next_event_s(&start),
            base.next_event_s(case.state),
            seconds,
        ),
    ];
    for (name, after, before, scale) in optional {
        match (after, before) {
            (Some(after), Some(before)) => numbers.push((name.to_owned(), after, before, scale)),
            (None, None) => {}
            _ => {
                return Err(law(
                    START_LAW,
                    format!("{name}: {after:?}, base {before:?}"),
                ));
            }
        }
    }
    in_new_units(START_LAW, numbers)?;
    Ok(seen)
}

/// Check the [`StageResidual`] tier at one operating point and its converged
/// stage unknowns `z_solution`.
///
/// Covered:
/// 1. Layout and vector lengths (`inputs`, `z_solution` of length `unknowns()`,
///    including zero for algebraic owners).
/// 2. `initial_guess` fills a finite starting vector.
/// 3. `residual` polls `cx` and refuses under a cancelled context.
/// 4. At `z_solution`, `residual` vanishes within `absolute_tolerance` and its
///    outputs agree with `propose` within `relative_tolerance * scale + absolute_tolerance`.
/// 5. `accept` builds a candidate state and receipt passing `audit_step` and
///    `admit`, and matching `propose`.
/// 6. `jvp_z` matches central finite differences of `residual` (both `dr` and
///    `d_outputs`) along three deterministic directions at a perturbed point
///    `z_probe` (so nonzero residuals and curvature are exercised).
/// 7. `vjp_z` is the exact transpose of `jvp_z`:
///    `<bar_r, dr> + <bar_outputs, d_outputs> == <bar_z, dz>`.
///
/// # Errors
/// The first [`Failure`] found.
#[allow(clippy::too_many_lines)]
pub fn check_stage_residual<C>(
    component: &C,
    case: &Case<'_, C::State>,
    z_solution: &[f64],
    relative_tolerance: f64,
    absolute_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<(), Failure<C::Error>>
where
    C: StageResidual,
    C::State: PartialEq + fmt::Debug,
{
    let ports = component.ports();
    let layout = PortLayout::new(ports).map_err(Failure::Contract)?;
    let (n_in, n_out) = (layout.input_len(), layout.output_len());
    if case.inputs.len() != n_in {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "inputs",
            expected: n_in,
            actual: case.inputs.len(),
        }));
    }
    let n_z = component.unknowns();
    if z_solution.len() != n_z {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "z",
            expected: n_z,
            actual: z_solution.len(),
        }));
    }

    // 1. Finite initial guess.
    let mut guess = vec![f64::NAN; n_z];
    component
        .initial_guess(case.state, case.dt_s, case.inputs, &mut guess)
        .map_err(Failure::Owner)?;
    if !guess.iter().all(|v| v.is_finite()) {
        return Err(law(
            "initial_guess fills every unknown with a finite value",
            format!("got {guess:?}"),
        ));
    }

    // 2. Cancellation polling on residual.
    let gate = CancelGate::new_clock_free();
    gate.request();
    let cancelled = Cx::new(&gate, cx.arena(), cx.stream_key(), cx.budget(), cx.mode());
    let mut cancelled_r = vec![f64::NAN; n_z];
    let mut cancelled_out = vec![f64::NAN; n_out];
    if component
        .residual(
            case.state,
            case.dt_s,
            case.inputs,
            z_solution,
            &mut cancelled_r,
            &mut cancelled_out,
            &cancelled,
        )
        .is_ok()
    {
        return Err(law(
            "StageResidual::residual polls cx: a call under a cancelled context is refused",
            "residual returned Ok under a pre-cancelled context".to_owned(),
        ));
    }

    // 3. Residual vanishes at z_solution and outputs match propose.
    let mut r_sol = vec![f64::NAN; n_z];
    let mut out_sol = vec![f64::NAN; n_out];
    component
        .residual(
            case.state,
            case.dt_s,
            case.inputs,
            z_solution,
            &mut r_sol,
            &mut out_sol,
            cx,
        )
        .map_err(Failure::Owner)?;
    for (idx, &val) in r_sol.iter().enumerate() {
        if !(val.is_finite() && val.abs() <= absolute_tolerance) {
            return Err(law(
                "residual vanishes at the converged stage unknowns",
                format!("r[{idx}] = {val:e}, allowed {absolute_tolerance:e}"),
            ));
        }
    }

    let mut prop_state = case.stale.clone();
    let mut prop_out = vec![f64::NAN; n_out];
    let Proposal::Step(_) = component
        .propose(
            case.state,
            case.dt_s,
            case.inputs,
            &mut prop_state,
            &mut prop_out,
            cx,
        )
        .map_err(Failure::Owner)?
    else {
        return Err(law(
            "StageResidual is checked where propose takes the step",
            "propose asked to shorten".to_owned(),
        ));
    };
    for (idx, (&a, &b)) in out_sol.iter().zip(&prop_out).enumerate() {
        let scale = a.abs().max(b.abs());
        if !((a - b).abs() <= relative_tolerance * scale + absolute_tolerance) {
            return Err(law(
                "residual outputs at z_solution match propose",
                format!("output {idx}: residual gave {a:e}, propose gave {b:e}"),
            ));
        }
    }

    // 4. accept matches audit_step, admit, and propose candidate state.
    let mut accepted_state = case.stale.clone();
    let receipt = component
        .accept(
            case.state,
            case.dt_s,
            case.inputs,
            z_solution,
            &mut accepted_state,
            cx,
        )
        .map_err(Failure::Owner)?;
    if accepted_state != prop_state {
        return Err(law(
            "accept at z_solution produces the same candidate state as propose",
            format!("accept gave {accepted_state:?}, propose gave {prop_state:?}"),
        ));
    }
    let stored_old = component
        .stored_energy_j(case.state)
        .map_err(Failure::Owner)?;
    let stored_new = component
        .stored_energy_j(&accepted_state)
        .map_err(Failure::Owner)?;
    audit_step(
        ports,
        &layout,
        case.dt_s,
        case.inputs,
        &out_sol,
        &receipt,
        stored_old,
        stored_new,
    )
    .map_err(Failure::Contract)?;
    component
        .admit(case.state, &accepted_state, &receipt)
        .map_err(Failure::Owner)?;

    // 5. jvp_z against central differences and vjp_z transpose dot-product test.
    let z_scale = z_solution.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
    let probe_shift = direction(n_z, 999);
    let z_probe: Vec<f64> = z_solution
        .iter()
        .zip(&probe_shift)
        .map(|(z, s)| z + 1.0e-3 * z_scale * s)
        .collect();
    let step = z_scale * f64::EPSILON.cbrt();
    let mut r_plus = vec![0.0; n_z];
    let mut r_minus = vec![0.0; n_z];
    let mut out_plus = vec![0.0; n_out];
    let mut out_minus = vec![0.0; n_out];
    let mut dr = vec![0.0; n_z];
    let mut d_out = vec![0.0; n_out];
    let mut bar_z = vec![0.0; n_z];

    for seed in 1..=3_u64 {
        let dz = direction(n_z, seed);
        let z_plus: Vec<f64> = z_probe.iter().zip(&dz).map(|(z, d)| z + step * d).collect();
        let z_minus: Vec<f64> = z_probe.iter().zip(&dz).map(|(z, d)| z - step * d).collect();
        component
            .residual(
                case.state,
                case.dt_s,
                case.inputs,
                &z_plus,
                &mut r_plus,
                &mut out_plus,
                cx,
            )
            .map_err(Failure::Owner)?;
        component
            .residual(
                case.state,
                case.dt_s,
                case.inputs,
                &z_minus,
                &mut r_minus,
                &mut out_minus,
                cx,
            )
            .map_err(Failure::Owner)?;
        component
            .jvp_z(
                case.state,
                case.dt_s,
                case.inputs,
                &z_probe,
                &dz,
                &mut dr,
                &mut d_out,
                cx,
            )
            .map_err(Failure::Owner)?;

        let (mut r_err, mut r_size) = (0.0_f64, 0.0_f64);
        for ((p, m), j) in r_plus.iter().zip(&r_minus).zip(&dr) {
            let fd = (p - m) / (2.0 * step);
            r_err = r_err.max((fd - j).abs());
            r_size = r_size.max(fd.abs()).max(j.abs());
        }
        if !(r_err <= relative_tolerance * r_size + absolute_tolerance) {
            return Err(law(
                "jvp_z residual derivative matches finite differences",
                format!("direction {seed}: max diff {r_err:e}, derivative size {r_size:e}"),
            ));
        }

        let (mut o_err, mut o_size) = (0.0_f64, 0.0_f64);
        for ((p, m), j) in out_plus.iter().zip(&out_minus).zip(&d_out) {
            let fd = (p - m) / (2.0 * step);
            o_err = o_err.max((fd - j).abs());
            o_size = o_size.max(fd.abs()).max(j.abs());
        }
        if !(o_err <= relative_tolerance * o_size + absolute_tolerance) {
            return Err(law(
                "jvp_z output derivative matches finite differences",
                format!("direction {seed}: max diff {o_err:e}, derivative size {o_size:e}"),
            ));
        }

        let bar_r = direction(n_z, seed + 200);
        let bar_out = direction(n_out, seed + 300);
        component
            .vjp_z(
                case.state,
                case.dt_s,
                case.inputs,
                &z_probe,
                &bar_r,
                &bar_out,
                &mut bar_z,
                cx,
            )
            .map_err(Failure::Owner)?;
        let left: f64 = bar_r.iter().zip(&dr).map(|(a, b)| a * b).sum::<f64>()
            + bar_out.iter().zip(&d_out).map(|(a, b)| a * b).sum::<f64>();
        let right: f64 = bar_z.iter().zip(&dz).map(|(a, b)| a * b).sum();
        let magnitude: f64 = bar_r
            .iter()
            .zip(&dr)
            .map(|(a, b)| (a * b).abs())
            .sum::<f64>()
            + bar_out
                .iter()
                .zip(&d_out)
                .map(|(a, b)| (a * b).abs())
                .sum::<f64>()
            + bar_z
                .iter()
                .zip(&dz)
                .map(|(a, b)| (a * b).abs())
                .sum::<f64>();
        if !((left - right).abs() <= 64.0 * f64::EPSILON * magnitude + f64::MIN_POSITIVE) {
            return Err(law(
                "vjp_z is the transpose of jvp_z",
                format!("<bar_r, dr> + <bar_out, d_out> = {left}, <bar_z, dz> = {right}"),
            ));
        }
    }
    Ok(())
}

/// Output of a passing [`check_temporal_order`] convergence ladder.
#[derive(Debug, Clone, PartialEq)]
pub struct TemporalOrderReport {
    /// Passing `fs-mms` convergence verdict (`OrderGate`).
    pub verdict: OrderVerdict,
    /// Step sizes `dt = horizon_s / n` for each rung [s], strictly decreasing.
    pub dts_s: Vec<f64>,
    /// Final-time error norms corresponding to each rung.
    pub errors: Vec<f64>,
}

/// Check a component's temporal convergence order over a fixed `horizon_s` using
/// `fs-mms` ([`RefinementLadder`] and [`OrderGate`]).
///
/// For each step count `n` in `step_counts` (at least three strictly increasing
/// positive integers), the component is marched `n` steps of length
/// `dt = horizon_s / n` from `initial_state` under `inputs_at(t, dt, inputs)`.
/// Every step is audited through [`audit_step`] and [`Component::admit`]; the
/// final state is scored by `error_of(&final_state, horizon_s)`. The resulting
/// `(dt, error)` ladder is gated against `expected_order` within
/// [`fs_mms::ORDER_GATE_TOLERANCE`] (0.2).
///
/// # Errors
/// Refuses invalid horizon/step counts, shortened steps, audit/owner failures,
/// non-positive errors, or a fitted order outside the `fs-mms` gate.
#[allow(clippy::too_many_arguments)]
pub fn check_temporal_order<C, F, G>(
    component: &C,
    case_name: &str,
    initial_state: &C::State,
    horizon_s: f64,
    step_counts: &[usize],
    expected_order: f64,
    mut inputs_at: F,
    mut error_of: G,
    cx: &Cx<'_>,
) -> Result<TemporalOrderReport, Failure<C::Error>>
where
    C: Component,
    F: FnMut(f64, f64, &mut [f64]),
    G: FnMut(&C::State, f64) -> f64,
{
    let ports = component.ports();
    let layout = PortLayout::new(ports).map_err(Failure::Contract)?;
    if !(horizon_s.is_finite() && horizon_s > 0.0) {
        return Err(law(
            "temporal order requires a positive finite horizon",
            format!("got horizon_s = {horizon_s:e}"),
        ));
    }
    if step_counts.len() < 3
        || step_counts.contains(&0)
        || !step_counts.windows(2).all(|w| w[1] > w[0])
    {
        return Err(law(
            "temporal order requires at least three strictly increasing positive step counts",
            format!("got step_counts = {step_counts:?}"),
        ));
    }

    let mut dts_s = Vec::with_capacity(step_counts.len());
    let mut errors = Vec::with_capacity(step_counts.len());
    let mut inputs = vec![0.0; layout.input_len()];
    let mut outputs = vec![0.0; layout.output_len()];

    for &steps in step_counts {
        let dt_s = horizon_s / (steps as f64);
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(law(
                "temporal order requires representable positive step sizes",
                format!("horizon_s={horizon_s:e} / steps={steps} produced dt_s={dt_s:e}"),
            ));
        }
        let mut state = initial_state.clone();
        let mut candidate = initial_state.clone();
        for step in 0..steps {
            let t_s = (step as f64) * dt_s;
            inputs.fill(0.0);
            inputs_at(t_s, dt_s, &mut inputs);
            let proposal = component
                .propose(&state, dt_s, &inputs, &mut candidate, &mut outputs, cx)
                .map_err(Failure::Owner)?;
            let receipt = match proposal {
                Proposal::Step(receipt) => receipt,
                Proposal::Shorten { dt_s: shorter } => {
                    return Err(law(
                        "temporal order is checked where the component takes every step",
                        format!("rung {steps} step {step} asked to shorten to {shorter:e} s"),
                    ));
                }
            };
            let stored_old = component.stored_energy_j(&state).map_err(Failure::Owner)?;
            let stored_new = component
                .stored_energy_j(&candidate)
                .map_err(Failure::Owner)?;
            audit_step(
                ports, &layout, dt_s, &inputs, &outputs, &receipt, stored_old, stored_new,
            )
            .map_err(Failure::Contract)?;
            let mass_old = component.stored_mass_kg(&state).map_err(Failure::Owner)?;
            let mass_new = component
                .stored_mass_kg(&candidate)
                .map_err(Failure::Owner)?;
            audit_mass_step(ports, &layout, dt_s, &inputs, &outputs, mass_old, mass_new)
                .map_err(Failure::Contract)?;
            component
                .admit(&state, &candidate, &receipt)
                .map_err(Failure::Owner)?;
            core::mem::swap(&mut state, &mut candidate);
        }
        let err = error_of(&state, horizon_s);
        dts_s.push(dt_s);
        errors.push(err);
    }

    let ladder = RefinementLadder::new(dts_s.clone(), errors.clone()).map_err(|e| {
        law(
            "temporal order matches theoretical convergence rate",
            e.to_string(),
        )
    })?;
    let verdict = OrderGate {
        theoretical: expected_order,
    }
    .check(case_name, LadderSide::Primal, &ladder)
    .map_err(|e| {
        law(
            "temporal order matches theoretical convergence rate",
            e.to_string(),
        )
    })?;

    Ok(TemporalOrderReport {
        verdict,
        dts_s,
        errors,
    })
}

/// Owned operating point built from a generated/shrunk property-test input.
#[derive(Debug, Clone)]
pub struct OwnedCase<S> {
    /// Accepted state the step starts from.
    pub state: S,
    /// Stale candidate state used for candidate-slot purity verification.
    pub stale: S,
    /// Step length [s].
    pub dt_s: f64,
    /// Held inputs in port-layout order.
    pub inputs: Vec<f64>,
}

impl<S> OwnedCase<S> {
    /// Borrow this owned operating point as a [`Case`].
    #[must_use]
    pub fn as_case(&self) -> Case<'_, S> {
        Case {
            state: &self.state,
            stale: &self.stale,
            dt_s: self.dt_s,
            inputs: &self.inputs,
        }
    }
}

/// Configuration for [`check_component_sweep`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SweepConfig {
    /// Deterministic `fs-propcheck` suite seed.
    pub suite_seed: u64,
    /// Number of generated cases to check.
    pub cases: u64,
    /// Maximum accepted shrink descents on a failing input.
    pub max_shrink_steps: usize,
    /// Maximum property evaluations during shrinking.
    pub max_shrink_evaluations: usize,
}

impl Default for SweepConfig {
    fn default() -> Self {
        Self {
            suite_seed: 0x05EE_D511,
            cases: 64,
            max_shrink_steps: 64,
            max_shrink_evaluations: 512,
        }
    }
}

/// Minimized counterexample returned when [`check_component_sweep`] catches a
/// conformance violation.
#[derive(Debug)]
pub struct ShrunkFailure<T, E> {
    /// Suite seed used to run the sweep.
    pub suite_seed: u64,
    /// Index of the first failing generated case.
    pub case_index: u64,
    /// Original unshrunk failing input.
    pub initial_input: T,
    /// Minimized failing input produced by `fs-propcheck` shrinking.
    pub minimized: T,
    /// Accepted shrink steps taken.
    pub shrink_steps: usize,
    /// Total evaluations spent during shrinking.
    pub shrink_evaluations: usize,
    /// Whether shrinking reached a fixpoint within budget.
    pub converged: bool,
    /// Exact [`Failure`] observed at `minimized`.
    pub failure: Failure<E>,
}

impl<T: fmt::Debug, E: fmt::Debug> fmt::Display for ShrunkFailure<T, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "conformance sweep failed on case {} (suite_seed={:#x}, shrunk in {} steps / {} evals, converged={}): \
             minimized={:?}; {}",
            self.case_index,
            self.suite_seed,
            self.shrink_steps,
            self.shrink_evaluations,
            self.converged,
            self.minimized,
            self.failure
        )
    }
}

/// Discriminant tag of a [`Failure`] used to keep shrinking locked onto the
/// same violated law or contract class (preventing shrinks toward zero from
/// escaping into unrelated out-of-domain owner refusals).
fn same_failure_kind<E>(a: &Failure<E>, b: &Failure<E>) -> bool {
    match (a, b) {
        (Failure::Owner(_), Failure::Owner(_)) => true,
        (Failure::Law { law: la, .. }, Failure::Law { law: lb, .. }) => la == lb,
        (Failure::Contract(ca), Failure::Contract(cb)) => {
            core::mem::discriminant(ca) == core::mem::discriminant(cb)
        }
        _ => false,
    }
}

/// Run [`check_component`] across `config.cases` deterministic `fs-propcheck`
/// seeded operating points, shrinking the first counterexample via
/// [`fs_propcheck::minimize_with_budget`] while staying locked onto the same
/// failure class.
///
/// Returns the number of passing cases on success, or [`ShrunkFailure`] naming
/// the replay key, minimized input, and violated law.
///
/// # Errors
/// Returns [`ShrunkFailure`] on the first conformance violation.
pub fn check_component_sweep<C, T, G, B>(
    component: &C,
    config: SweepConfig,
    mut generate: G,
    build_case: B,
    cx: &Cx<'_>,
) -> Result<u64, Box<ShrunkFailure<T, C::Error>>>
where
    C: Component,
    C::State: PartialEq + fmt::Debug,
    T: Shrink + fmt::Debug,
    G: FnMut(&mut Stream) -> T,
    B: Fn(&T) -> OwnedCase<C::State>,
{
    for case_index in 0..config.cases {
        let mut stream = Stream::for_case(config.suite_seed, case_index);
        let initial = generate(&mut stream);
        let owned = build_case(&initial);
        if let Err(initial_err) = check_component(component, &owned.as_case(), cx) {
            let budget = MinimizeBudget {
                max_steps: config.max_shrink_steps,
                max_evaluations: config.max_shrink_evaluations.max(1),
                max_candidates_per_step: fs_propcheck::DEFAULT_MAX_SHRINK_CANDIDATES_PER_STEP,
            };
            let report = fs_propcheck::minimize_with_budget(
                initial.clone(),
                |cand| {
                    let built = build_case(cand);
                    match check_component(component, &built.as_case(), cx) {
                        Ok(_) => true,
                        Err(err) => !same_failure_kind(&err, &initial_err),
                    }
                },
                budget,
            );
            let (minimized, shrink_steps, shrink_evaluations, converged) = match report {
                Ok(r) => (r.minimized, r.steps, r.tried, r.converged),
                Err(_) => (initial.clone(), 0, 1, false),
            };
            let min_owned = build_case(&minimized);
            let failure = match check_component(component, &min_owned.as_case(), cx) {
                Err(err) => err,
                Ok(_) => law(
                    "propose is deterministic across shrink replay",
                    format!("minimized candidate {minimized:?} passed on replay"),
                ),
            };
            return Err(Box::new(ShrunkFailure {
                suite_seed: config.suite_seed,
                case_index,
                initial_input: initial,
                minimized,
                shrink_steps,
                shrink_evaluations,
                converged,
                failure,
            }));
        }
    }
    Ok(config.cases)
}

/// Conformance-kit zero-allocation gate for [`Component::propose`] (bead `.15.1`):
/// runs `warmup_steps` calls and then asserts zero heap allocations across
/// `measured_steps` repeated calls, using a caller-supplied counting-allocator
/// measurement closure `measure_allocations`.
///
/// Keeping allocation counting injected avoids any `unsafe` global allocator
/// inside `fs-ports` while letting every test binary enforce the exact same law.
///
/// # Errors
/// Fails if layout/input sizes mismatch, if `measured_steps == 0`, if `propose`
/// returns an error, or if `measure_allocations` reports any nonzero allocation.
pub fn check_zero_alloc_propose<C, M>(
    component: &C,
    case: &Case<'_, C::State>,
    warmup_steps: usize,
    measured_steps: usize,
    mut measure_allocations: M,
    cx: &Cx<'_>,
) -> Result<(), Failure<C::Error>>
where
    C: Component,
    M: FnMut(&mut dyn FnMut()) -> usize,
{
    let layout = PortLayout::new(component.ports()).map_err(Failure::Contract)?;
    if case.inputs.len() != layout.input_len() {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "inputs",
            expected: layout.input_len(),
            actual: case.inputs.len(),
        }));
    }
    if measured_steps == 0 {
        return Err(law(
            "zero-allocation check requires at least one measured call",
            "measured_steps was 0".to_owned(),
        ));
    }
    let mut candidate = case.stale.clone();
    let mut outputs = vec![0.0; layout.output_len()];
    for _ in 0..warmup_steps {
        component
            .propose(
                case.state,
                case.dt_s,
                case.inputs,
                &mut candidate,
                &mut outputs,
                cx,
            )
            .map_err(Failure::Owner)?;
    }
    let mut owner_err: Option<C::Error> = None;
    let allocs = measure_allocations(&mut || {
        for _ in 0..measured_steps {
            if let Err(err) = component.propose(
                case.state,
                case.dt_s,
                case.inputs,
                &mut candidate,
                &mut outputs,
                cx,
            ) {
                owner_err = Some(err);
                break;
            }
        }
    });
    if let Some(err) = owner_err {
        return Err(Failure::Owner(err));
    }
    if allocs > 0 {
        return Err(law(
            "propose allocates no heap memory after warm-up",
            format!("{allocs} heap allocations across {measured_steps} warm calls"),
        ));
    }
    Ok(())
}
