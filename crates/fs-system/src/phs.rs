//! Original port-Hamiltonian owners on the composition kernel.
//!
//! [`PortHamiltonianComponent`] delegates every duration to the owner's
//! variable-step discrete-gradient solver. It retains all coordinates, source
//! storage laws and power ports, and sends the actual resistance heat through
//! a mandatory thermal port. No integration equation is defined here.
//!
//! This is not the exact-exponential tier: nonlinear heat makes the complete
//! port map [`Class::Smooth`], and no port tangent is claimed. Source storage
//! must be differentiable on its admitted support. Switching and contact laws
//! need their original event treatment and are outside this smooth adapter.
//! The native Newton budget and source-domain refusals remain authoritative.

use core::cell::RefCell;
use core::fmt;

use fs_exec::Cx;
use fs_phs::{
    PhsError, PortHamiltonian, PreparedStepError, SolveBudget, StepWorkspace,
    arithmetic_rounding_floor,
};
use fs_ports::component::{
    Carrier, Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{CoupleError, StableId};

/// Refusal from the original owner or its Component boundary.
#[derive(Debug)]
pub enum PhsComponentError {
    /// The unchanged numerical-owner refusal, including cooperative cancellation.
    Owner(PreparedStepError),
    /// Invalid or repeated port metadata.
    Port(ContractViolation),
    /// An invalid port identifier.
    Identifier(CoupleError),
    /// A shape, value or reentrant workspace request was refused.
    Invalid {
        /// The rejected quantity or operation.
        what: &'static str,
    },
    /// Native work, heat and storage do not close within equation and arithmetic error.
    Balance {
        /// Supplied work minus physical heat minus stored-energy change [J].
        defect_j: f64,
        /// Allowed equation-work and arithmetic error [J].
        tolerance_j: f64,
    },
}

impl fmt::Display for PhsComponentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Owner(error) => write!(f, "port-Hamiltonian owner refused: {error}"),
            Self::Port(error) => write!(f, "port-Hamiltonian port refused: {error:?}"),
            Self::Identifier(error) => write!(f, "port identifier refused: {error:?}"),
            Self::Invalid { what } => write!(f, "port-Hamiltonian component refused: {what}"),
            Self::Balance {
                defect_j,
                tolerance_j,
            } => {
                write!(
                    f,
                    "port-Hamiltonian balance defect {defect_j} J exceeds {tolerance_j} J"
                )
            }
        }
    }
}

impl std::error::Error for PhsComponentError {}

impl From<PhsError> for PhsComponentError {
    fn from(error: PhsError) -> Self {
        Self::Owner(PreparedStepError::Solver(error))
    }
}

fn invalid(what: &'static str) -> PhsComponentError {
    PhsComponentError::Invalid { what }
}

struct Scratch {
    solver: StepWorkspace,
    next: Vec<f64>,
    output: Vec<f64>,
}

/// A complete original port-Hamiltonian model and its mandatory heat outlet.
///
/// State is the original coordinate vector. The caller supplies one scalar
/// power port per original input, in the same order: `EffortIn` interprets
/// the native input as effort and output as flow; `FlowIn` does the reverse.
/// No sign, coordinate or unit conversion is inferred. The final `heat` port
/// receives temperature [K] and returns `-D/dt` [W], where D is the native
/// resistance work, never an energy remainder.
///
/// Numerical scratch is private and contains no accepted state or clock.
/// The owner resets its iteration and re-reads the supplied state, inputs,
/// duration and operators on every call. Refused trials therefore require no
/// cache rollback. Reentrant use refuses instead of panicking. Warm proposals
/// reuse buffers when the candidate has the original shape; allocation inside
/// user-provided storage callbacks remains the owner's responsibility.
pub struct PortHamiltonianComponent {
    system: PortHamiltonian,
    ports: Vec<PortSpec>,
    scratch: RefCell<Scratch>,
}

impl PortHamiltonianComponent {
    /// Wrap the original model without changing its operators or storage.
    ///
    /// # Errors
    /// Requires exactly one scalar, ordinary power port per native input,
    /// unique identifiers other than `heat`, and admitted native workspace
    /// dimensions. Field, stream, heat and signal ports are not native power
    /// inputs. Native structural and storage validity remain owner obligations.
    /// The source must define a smooth port map on its admitted support;
    /// switching and contact laws requiring event treatment are unsupported.
    pub fn new(
        system: PortHamiltonian,
        mut power_ports: Vec<PortSpec>,
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
        let heat = PortSpec::heat(
            StableId::new("heat").map_err(PhsComponentError::Identifier)?,
            PortRole::EffortIn,
        )
        .map_err(PhsComponentError::Port)?;
        power_ports.push(heat);
        PortLayout::new(&power_ports).map_err(PhsComponentError::Port)?;
        let scratch = Scratch {
            solver: StepWorkspace::new(&system)?,
            next: vec![0.0; system.state_dim()],
            output: vec![0.0; system.port_dim()],
        };
        Ok(Self {
            system,
            ports: power_ports,
            scratch: RefCell::new(scratch),
        })
    }

    /// Original immutable numerical owner, including its complete storage law.
    #[must_use]
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Set the native residual and Newton-update limits before composition.
    /// Explicit budgets retain the owner's below-rounding-floor refusal.
    pub fn set_solve_budget(&mut self, budget: SolveBudget) -> Result<(), PhsComponentError> {
        self.scratch.get_mut().solver.set_solve_budget(budget)?;
        Ok(())
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
}

impl Component for PortHamiltonianComponent {
    type State = Vec<f64>;
    type Error = PhsComponentError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        // Even with quadratic storage, the heat outlet is quadratic in input.
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, PhsComponentError> {
        self.energy(state)
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
        let cancelled = || PhsComponentError::Owner(PreparedStepError::Cancelled);
        cx.checkpoint().map_err(|_| cancelled())?;
        let m = self.system.port_dim();
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != m + 1
            || outputs.len() != m + 1
            || inputs.iter().any(|v| !v.is_finite())
            || inputs[m] <= 0.0
        {
            return Err(invalid(
                "positive duration, complete finite ports and positive temperature required",
            ));
        }
        let energy_old = self.energy(old)?;
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| invalid("reentrant use of numerical workspace"))?;
        let Scratch {
            solver,
            next,
            output,
        } = &mut *scratch;
        let record = solver
            .step_into_controlled(&self.system, old, &inputs[..m], dt_s, next, output, || {
                cx.is_cancel_requested()
            })
            .map_err(PhsComponentError::Owner)?;
        let energy_new = self.energy(next)?;
        let effort = solver
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("native accepted discrete effort is unavailable"))?;
        let certificate = record
            .solver_certificate
            .ok_or_else(|| invalid("native equation certificate is unavailable"))?;

        // Preserve weak work under cancellation in the original G map. The
        // wire uses rounded outputs, while native work retains low words.
        let (_, _, map) = self.system.structure();
        let mut magnitude = energy_old.abs()
            + energy_new.abs()
            + record.delta_h.abs()
            + record.supplied.abs()
            + record.dissipated.abs();
        for (row, &e) in effort.iter().enumerate() {
            cx.checkpoint().map_err(|_| cancelled())?;
            for (column, &u) in inputs[..m].iter().enumerate() {
                magnitude += ((map[row * m + column] * e) * u * dt_s).abs();
            }
        }
        let rounding = arithmetic_rounding_floor(
            32_usize
                .saturating_add(effort.len().saturating_mul(8))
                .saturating_add(m.saturating_mul(8)),
            magnitude,
        );
        // A residual in coordinate i contributes effort_i * residual_i
        // to the first law. Nothing beyond this equation-work allowance is
        // silently relabelled as numerical loss.
        let residual_work = effort
            .iter()
            .map(|e| e.abs() * record.solver_residual)
            .sum::<f64>();
        let tolerance = rounding + residual_work;
        let numerical_loss = record.supplied - record.dissipated - record.delta_h;
        let heat_rate = -record.dissipated / dt_s;
        let entropy = record.dissipated / inputs[m];
        if ![rounding, tolerance, numerical_loss, heat_rate, entropy]
            .iter()
            .all(|v| v.is_finite())
            || record.dissipated < 0.0
        {
            return Err(invalid(
                "physical heat and energy account must be finite and nonnegative",
            ));
        }
        if numerical_loss.abs() > tolerance {
            return Err(PhsComponentError::Balance {
                defect_j: numerical_loss,
                tolerance_j: tolerance,
            });
        }
        cx.checkpoint().map_err(|_| cancelled())?;
        new.clone_from(next);
        outputs[..m].copy_from_slice(output);
        outputs[m] = heat_rate;
        Ok(Proposal::Step(Receipt {
            stored_change_j: record.delta_h,
            dissipated_j: record.dissipated,
            numerical_loss_j: numerical_loss,
            entropy_generated_j_per_k: entropy,
            energy_tolerance_j: rounding,
            residual: record.solver_residual,
            tolerance: certificate.effective_budget,
        }))
    }
}
