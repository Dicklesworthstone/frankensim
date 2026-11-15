//! Exact held-input port-Hamiltonian owners with complete port derivatives.
//!
//! The original exact quadratic owner prepares every requested duration and
//! supplies the state, mean ports and independently integrated physical heat.
//! A bounded cache reuses two durations, such as adaptive coarse/half trials.
//! The complete map is smooth, including its quadratic heat outlet.
//!
//! The source types this owner is built from are re-exported, so a caller
//! can assemble one without depending on `fs-phs` itself.
//!
//! The end of a step. The step is exact for a held input, not the midpoint
//! rule, so a port's output does not end its step at twice its mean less its
//! start. A port answers [`EndOfStep::Untold`], and a system does not join it
//! to a load line, unless the owner publishes its ends
//! ([`ExactPortHamiltonianComponent::publishing_ends`]): one signal output per
//! power port, `<id>.end`, returning `y1 = G^T grad H(x1)` at the step's end
//! state, which the port then names ([`EndOfStep::Told`]).

use core::cell::RefCell;
use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_phs::linear_exact::{ExactZohWorkspace, PreparedExactZoh};
use fs_phs::{PreparedStepError, arithmetic_rounding_floor};
use fs_ports::StableId;
use fs_ports::component::{
    Carrier, Class, Component, EndOfStep, PortLayout, PortRole, PortSpec, PortTangent, Proposal,
    Receipt,
};

pub use fs_phs::linear_exact::LinearSplitBudget;
pub use fs_phs::{PortHamiltonian, QuadraticStorage};

use crate::phs::PhsComponentError;

const CACHED_DURATIONS: usize = 2;

fn invalid(what: &'static str) -> PhsComponentError {
    PhsComponentError::Invalid { what }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), PhsComponentError> {
    cx.checkpoint()
        .map_err(|_| PhsComponentError::Owner(PreparedStepError::Cancelled))
}

struct Scratch {
    // Only immutable constitutive callbacks are shared. No accepted state or
    // clock lives in this prototype, a cache entry, or the native workspace.
    source: PortHamiltonian,
    cached: Vec<PreparedExactZoh<'static>>,
    workspace: ExactZohWorkspace,
    next: Vec<f64>,
    output: Vec<f64>,
    heat_gradient: Vec<f64>,
    tangent: Vec<f64>,
    // The published ends: the end state's effort, the end outputs, and the
    // row-major end rows of the tangent with the rest state, unit input and
    // step buffers they are found with.
    effort: Vec<f64>,
    end: Vec<f64>,
    end_rows: Vec<f64>,
    rest: Vec<f64>,
    unit: Vec<f64>,
    unit_next: Vec<f64>,
    unit_output: Vec<f64>,
}

impl Scratch {
    fn prepare(
        &mut self,
        dt_s: f64,
        budget: LinearSplitBudget,
        cx: &Cx<'_>,
    ) -> Result<usize, PhsComponentError> {
        checkpoint(cx)?;
        if let Some(index) = self
            .cached
            .iter()
            .position(|prepared| prepared.duration().to_bits() == dt_s.to_bits())
        {
            let last = self.cached.len() - 1;
            self.cached.swap(index, last);
            return Ok(last);
        }
        let prepared = PreparedExactZoh::new_owned_controlled(
            self.source.share_storage(),
            dt_s,
            budget,
            || cx.is_cancel_requested(),
        )
        .map_err(PhsComponentError::Owner)?;
        // Cold preparation must succeed before an old reusable entry is lost.
        if self.cached.len() == CACHED_DURATIONS {
            self.cached.remove(0);
        }
        self.cached.push(prepared);
        Ok(self.cached.len() - 1)
    }
}

/// A complete exact-quadratic owner, its physical heat outlet and port tangent.
///
/// Each original input has one caller-declared ordinary scalar power port in
/// the original order. The last port, `heat`, receives positive temperature
/// and emits `-D/dt`, where D is the original integrated resistance work.
/// Temperature affects entropy only, so its column in the port derivative is
/// exactly zero. The heat row includes the full analytic input derivative.
///
/// Cold preparation retains the source storage through `share_storage` and
/// the original dense capacity limits. Two exact-duration caches are retained
/// by least-recent use, keyed by the duration's IEEE bits. Warm proposals and
/// tangents reuse native scratch; accepted state always belongs to the caller.
/// Cancellation or refusal leaves caller outputs unchanged. This adapter
/// neither infers quadratic storage nor substitutes a discrete-gradient step.
///
/// Its ports answer [`EndOfStep::Untold`] unless it was built with
/// [`Self::publishing_ends`], which adds one signal output per power port
/// after `heat` carrying the port's output at the end of the step.
pub struct ExactPortHamiltonianComponent {
    system: PortHamiltonian,
    ports: Vec<PortSpec>,
    budget: LinearSplitBudget,
    tells_ends: bool,
    scratch: RefCell<Scratch>,
}

impl ExactPortHamiltonianComponent {
    /// Admit the original complete quadratic model and its power-port metadata.
    ///
    /// # Errors
    /// Refuses nonquadratic storage, native structure/capacity refusals, wrong
    /// port count, non-scalar or non-power ports, field metrics and repeated
    /// identifiers. `heat` is reserved for the mandatory thermal outlet.
    /// The constructor performs cold zero-duration admission; positive trial
    /// durations are prepared under the proposal's cancellation capability.
    pub fn new(
        mut system: PortHamiltonian,
        mut power_ports: Vec<PortSpec>,
        budget: LinearSplitBudget,
    ) -> Result<Self, PhsComponentError> {
        if power_ports.len() != system.port_dim()
            || power_ports.iter().any(|port| {
                !matches!(port.carrier(), Carrier::Power(_))
                    || !matches!(port.role(), PortRole::EffortIn | PortRole::FlowIn)
                    || port.width() != 1
                    || port.field_metric().is_some()
                    || port.id().as_str() == "heat"
            })
        {
            return Err(invalid(
                "one ordinary scalar power port per native input; heat is reserved",
            ));
        }
        power_ports.push(
            PortSpec::heat(
                StableId::new("heat").map_err(PhsComponentError::Identifier)?,
                PortRole::EffortIn,
            )
            .map_err(PhsComponentError::Port)?,
        );
        PortLayout::new(&power_ports).map_err(PhsComponentError::Port)?;
        let workspace = PreparedExactZoh::new(&system, 0.0, budget)?.workspace();
        let (n, m) = (system.state_dim(), system.port_dim());
        let scratch = Scratch {
            source: system.share_storage(),
            cached: Vec::with_capacity(CACHED_DURATIONS),
            workspace,
            next: vec![0.0; n],
            output: vec![0.0; m],
            heat_gradient: vec![0.0; m],
            tangent: vec![0.0; m + 1],
            effort: vec![0.0; n],
            end: vec![0.0; m],
            end_rows: Vec::new(),
            rest: vec![0.0; n],
            unit: vec![0.0; m],
            unit_next: vec![0.0; n],
            unit_output: vec![0.0; m],
        };
        Ok(Self {
            system,
            ports: power_ports,
            budget,
            tells_ends: false,
            scratch: RefCell::new(scratch),
        })
    }

    /// The same owner telling the end of each power port's output: one signal
    /// output per power port after `heat`, `<id>.end`, of width one, returning
    /// the port's output `y1 = G^T grad H(x1)` at the step's end state. A port
    /// that returns a flow then answers [`EndOfStep::Told`] with that output, so
    /// a system may join it to a load line. Its tangent rows are exact: the end
    /// output is linear in the held inputs. Asking twice adds them once.
    ///
    /// # Errors
    /// [`PhsComponentError::Port`] when an end's name repeats a port (a port
    /// named `a.end` beside one named `a`).
    pub fn publishing_ends(mut self) -> Result<Self, PhsComponentError> {
        if self.tells_ends {
            return Ok(self);
        }
        let m = self.system.port_dim();
        for port in 0..m {
            let id = StableId::new(format!("{}.end", self.ports[port].id().as_str()))
                .map_err(PhsComponentError::Identifier)?;
            self.ports.push(PortSpec::signal_out(id, NonZeroUsize::MIN));
        }
        PortLayout::new(&self.ports).map_err(PhsComponentError::Port)?;
        let scratch = self.scratch.get_mut();
        scratch.tangent = vec![0.0; 2 * m + 1];
        scratch.end_rows = vec![0.0; m * m];
        self.tells_ends = true;
        Ok(self)
    }

    /// The unchanged physical operators and original shared storage law.
    #[must_use]
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    fn energy(&self, state: &[f64]) -> Result<f64, PhsComponentError> {
        if state.len() != self.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite original state required"));
        }
        let energy = self.system.hamiltonian(state);
        if !energy.is_finite() {
            return Err(invalid("original Hamiltonian is not finite"));
        }
        Ok(energy)
    }

    fn inputs(&self, dt_s: f64, inputs: &[f64]) -> Result<(), PhsComponentError> {
        let m = self.system.port_dim();
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != m + 1
            || inputs.iter().any(|v| !v.is_finite())
            || inputs[m] <= 0.0
        {
            return Err(invalid(
                "positive duration, complete finite ports and positive temperature required",
            ));
        }
        Ok(())
    }

    fn energy_error_bound(&self, state: &[f64]) -> Result<f64, PhsComponentError> {
        let bound = self
            .system
            .energy_evaluation_error_bound_j(state)
            .unwrap_or(0.0);
        if !(bound.is_finite() && bound >= 0.0) {
            return Err(invalid(
                "source energy evaluation error bound must be finite and nonnegative",
            ));
        }
        Ok(bound)
    }

    /// Proposals and tangents share the complete native physical admission.
    /// Scratch may change; no caller state, derivative or port is published here.
    fn trial(
        &self,
        old: &[f64],
        dt_s: f64,
        inputs: &[f64],
        scratch: &mut Scratch,
        cx: &Cx<'_>,
    ) -> Result<(usize, Receipt), PhsComponentError> {
        self.inputs(dt_s, inputs)?;
        let energy_old = self.energy(old)?;
        let index = scratch.prepare(dt_s, self.budget, cx)?;
        let m = self.system.port_dim();
        let prepared = &scratch.cached[index];
        let record = prepared
            .step_with_workspace_into_controlled(
                &mut scratch.workspace,
                old,
                &inputs[..m],
                &mut scratch.next,
                &mut scratch.output,
                || cx.is_cancel_requested(),
            )
            .map_err(PhsComponentError::Owner)?;
        let energy_new = self.energy(&scratch.next)?;
        // Native work contracts both words of each port integral before the
        // wire rounds its mean. Retain their absolute magnitude even when
        // large opposing port works leave a tiny signed net supply.
        let magnitude = energy_old.abs()
            + energy_new.abs()
            + record.delta_h.abs()
            + record.supplied.abs()
            + record.dissipated.abs()
            + prepared.port_work_magnitude(&scratch.workspace, &inputs[..m])?;
        let rounding = arithmetic_rounding_floor(
            32_usize
                .saturating_add(self.system.state_dim().saturating_mul(8))
                .saturating_add(m.saturating_mul(8)),
            magnitude,
        );
        // A source may evaluate H as large cancelling terms. Its declared
        // endpoint evaluation bounds concern the Component remeasurement,
        // never the native independent work/heat/energy balance admission.
        let energy_tolerance =
            rounding + self.energy_error_bound(old)? + self.energy_error_bound(&scratch.next)?;
        let numerical_loss = record.supplied - record.dissipated - record.delta_h;
        let heat_rate = -record.dissipated / dt_s;
        let entropy = record.dissipated / inputs[m];
        if ![energy_tolerance, numerical_loss, heat_rate, entropy]
            .iter()
            .all(|value| value.is_finite())
            || record.dissipated < 0.0
        {
            return Err(invalid("unrepresentable physical heat or energy account"));
        }
        checkpoint(cx)?;
        Ok((
            index,
            Receipt {
                stored_change_j: record.delta_h,
                dissipated_j: record.dissipated,
                numerical_loss_j: numerical_loss,
                entropy_generated_j_per_k: entropy,
                energy_tolerance_j: energy_tolerance,
                residual: record.solver_residual,
                tolerance: 0.0,
            },
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn derivative(
        &self,
        old: &[f64],
        dt_s: f64,
        inputs: &[f64],
        direction: &[f64],
        result: &mut [f64],
        transpose: bool,
        cx: &Cx<'_>,
    ) -> Result<(), PhsComponentError> {
        checkpoint(cx)?;
        let m = self.system.port_dim();
        let ends = if self.tells_ends { m } else { 0 };
        // Inputs are the ports' and the temperature; outputs add the heat and
        // the published ends. A JVP maps inputs to outputs, a VJP back.
        let (inputs_len, outputs_len) = (m + 1, m + 1 + ends);
        let (direction_len, result_len) = if transpose {
            (outputs_len, inputs_len)
        } else {
            (inputs_len, outputs_len)
        };
        if direction.len() != direction_len
            || result.len() != result_len
            || direction.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite port derivative direction required",
            ));
        }
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| invalid("reentrant use of numerical workspace"))?;
        let (index, _) = self.trial(old, dt_s, inputs, &mut scratch, cx)?;
        if ends > 0 {
            self.end_rows(index, &mut scratch, cx)?;
        }
        let Scratch {
            cached,
            workspace,
            heat_gradient,
            tangent,
            end_rows,
            ..
        } = &mut *scratch;
        let prepared = &cached[index];
        prepared
            .heat_input_gradient_into_controlled(
                workspace,
                old,
                &inputs[..m],
                heat_gradient,
                || cx.is_cancel_requested(),
            )
            .map_err(PhsComponentError::Owner)?;
        tangent.fill(0.0);
        for row in 0..m {
            checkpoint(cx)?;
            let (_, input_row) = prepared.average_output_map(row)?;
            if transpose {
                for (result, coefficient) in tangent[..m].iter_mut().zip(input_row) {
                    *result += coefficient * direction[row];
                }
            } else {
                tangent[row] = input_row
                    .iter()
                    .zip(&direction[..m])
                    .map(|(coefficient, delta)| coefficient * delta)
                    .sum();
            }
        }
        if transpose {
            for (result, gradient) in tangent[..m].iter_mut().zip(heat_gradient) {
                *result -= (*gradient / dt_s) * direction[m];
            }
            // Every output is independent of the thermal input. In particular
            // the entropy account is not itself an output of the port map.
            tangent[m] = 0.0;
        } else {
            tangent[m] = 0.0
                - heat_gradient
                    .iter()
                    .zip(&direction[..m])
                    .map(|(gradient, delta)| (gradient / dt_s) * delta)
                    .sum::<f64>();
        }
        for row in 0..ends {
            checkpoint(cx)?;
            let coefficients = &end_rows[row * m..(row + 1) * m];
            if transpose {
                for (result, coefficient) in tangent[..m].iter_mut().zip(coefficients) {
                    *result += coefficient * direction[m + 1 + row];
                }
            } else {
                tangent[m + 1 + row] = coefficients
                    .iter()
                    .zip(&direction[..m])
                    .map(|(coefficient, delta)| coefficient * delta)
                    .sum();
            }
        }
        let tangent = &tangent[..result_len];
        if tangent.iter().any(|value| !value.is_finite()) {
            return Err(invalid("port derivative is not finite"));
        }
        checkpoint(cx)?;
        result.copy_from_slice(tangent);
        Ok(())
    }

    /// The end output `y1 = G^T grad H(x1)` of the state the last trial ended
    /// at, into `scratch.end`.
    fn end_output(&self, scratch: &mut Scratch) -> Result<(), PhsComponentError> {
        let Scratch {
            next, effort, end, ..
        } = scratch;
        self.system.storage().gradient(next, effort);
        self.system
            .output_from_effort_into(effort, end)
            .map_err(|_| invalid("unrepresentable end output"))?;
        if end.iter().any(|value| !value.is_finite()) {
            return Err(invalid("end output is not finite"));
        }
        Ok(())
    }

    /// The end rows of the tangent, `d y1 / d u`, row-major into
    /// `scratch.end_rows`. The step is affine in its held inputs and the
    /// output linear in the state, so column `j` is the end output of a step
    /// of the same duration from rest under a unit input on port `j`; the
    /// temperature does not enter.
    fn end_rows(
        &self,
        index: usize,
        scratch: &mut Scratch,
        cx: &Cx<'_>,
    ) -> Result<(), PhsComponentError> {
        let m = self.system.port_dim();
        let Scratch {
            cached,
            workspace,
            effort,
            end_rows,
            rest,
            unit,
            unit_next,
            unit_output,
            ..
        } = scratch;
        let prepared = &cached[index];
        for column in 0..m {
            checkpoint(cx)?;
            unit.fill(0.0);
            unit[column] = 1.0;
            prepared
                .step_with_workspace_into_controlled(
                    workspace,
                    rest,
                    unit,
                    unit_next,
                    unit_output,
                    || cx.is_cancel_requested(),
                )
                .map_err(PhsComponentError::Owner)?;
            self.system.storage().gradient(unit_next, effort);
            self.system
                .output_from_effort_into(effort, unit_output)
                .map_err(|_| invalid("unrepresentable end-output derivative"))?;
            for (row, value) in unit_output.iter().enumerate() {
                end_rows[row * m + column] = *value;
            }
        }
        Ok(())
    }
}

impl Component for ExactPortHamiltonianComponent {
    type State = Vec<f64>;
    type Error = PhsComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, PhsComponentError> {
        self.energy(state)
    }

    /// Exact for a held input, so never the midpoint relation: the end is
    /// told when the owner publishes its ends and untold otherwise.
    fn end_of_step(&self, port: usize) -> EndOfStep {
        let m = self.system.port_dim();
        if self.tells_ends && port < m {
            EndOfStep::Told(m + 1 + port)
        } else {
            EndOfStep::Untold
        }
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, PhsComponentError> {
        checkpoint(cx)?;
        let m = self.system.port_dim();
        let ends = if self.tells_ends { m } else { 0 };
        if outputs.len() != m + 1 + ends {
            return Err(invalid("complete port output buffer required"));
        }
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| invalid("reentrant use of numerical workspace"))?;
        let (_, receipt) = self.trial(old, dt_s, inputs, &mut scratch, cx)?;
        if ends > 0 {
            self.end_output(&mut scratch)?;
        }
        checkpoint(cx)?;
        new.clone_from(&scratch.next);
        outputs[..m].copy_from_slice(&scratch.output);
        outputs[m] = -receipt.dissipated_j / dt_s;
        if ends > 0 {
            outputs[m + 1..].copy_from_slice(&scratch.end);
        }
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for ExactPortHamiltonianComponent {
    fn jvp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), PhsComponentError> {
        self.derivative(old, dt_s, inputs, d_inputs, d_outputs, false, cx)
    }

    fn vjp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), PhsComponentError> {
        self.derivative(old, dt_s, inputs, bar_outputs, bar_inputs, true, cx)
    }
}
