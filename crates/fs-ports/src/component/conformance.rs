//! Conformance checks every [`Component`] implementation runs.
//!
//! One call per component replaces the self-consistency tests each pairing
//! used to write for itself. A failure names the law that was broken.
//!
//! Covered here: port layout, step-limit bounds, purity of `propose`
//! (including independence from the stale candidate slot), the step audit
//! (first law, endpoint energy, non-negative dissipation and entropy
//! production, residual), available mass inventories, `admit`, and for
//! [`PortTangent`] the JVP against
//! central differences and the VJP against the JVP by the dot-product test.
//! For [`LoadLine`]: the line's ports, purity, the audit at the flows on the
//! line, and that a zero compliance gives `propose`.
//!
//! Not covered yet (bead `frankensim-rc-1109-root-3gc12.5.11`):
//! [`StageResidual`](super::StageResidual) checks, unit rescaling, time
//! reversal of lossless components, and the temporal order hook.

use core::fmt;

use fs_exec::Cx;

use super::{
    Carrier, Component, ContractViolation, LoadLine, PortLayout, PortRole, PortTangent, Proposal,
    Receipt, StepAudit, audit_mass_step, audit_step,
};

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
    if compliance.len() != width * width {
        return Err(Failure::Contract(ContractViolation::VectorLength {
            vector: "compliance",
            expected: width * width,
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
    let zero = vec![0.0; width * width];
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
    Ok(())
}
