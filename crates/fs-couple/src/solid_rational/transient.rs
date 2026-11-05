//! Full spatial solid motion with the original rational acoustic load history.
//!
//! The native solid advances midpoint inertia and its exact material strain
//! ramp. Its actual port displacement increment divided by the unchanged clock
//! is held at the original MIMO source for that interval. The source's exact
//! mean effort enters the SAME mechanical solve through its complete mutual
//! mean-direct matrix. There is no explicit delayed reaction, force fit, new
//! rational integrator, structural mode truncation or direct-term inverse.
//!
//! This is a work-consistent time discretization, not the continuous joint
//! exponential used by the smaller assembly owner. Time refinement remains
//! necessary for loaded resonance and decay accuracy. Source exchange is signed
//! boundary work; arbitrary rational coordinates have no inferred heat/storage.

use fs_exec::Cx;
use fs_solid::linear3::visco::{
    ViscoDynamicPortBudget, ViscoDynamicPortOperator, ViscoDynamicPortStep, ViscoDynamicState,
    ViscoElementError,
};
use fs_vfit::mimo::{
    MimoPassivityReport,
    load::{MimoLoadError, PreparedMimoLoad},
};

/// Separate limits for the original spatial connection and original load.
#[derive(Clone, Copy, Debug)]
pub struct SolidMimoTransientBudget {
    /// Full native port map, cached sparse responses and small Schur factor.
    pub ports: ViscoDynamicPortBudget,
    /// Original admitted load allowance plus 2048 bytes of fixed mean-effort /
    /// history scratch and 32 bytes per port. Preview reuses the source's own
    /// candidate state; its model and reports are never cloned per interval.
    /// The body, constitutive history and their native candidate allowance are
    /// separate.
    pub maximum_load_payload_bytes: usize,
    /// Relative effort and signed exchange-work gate, in (0, 1e-8]. Scales use
    /// the uncancelled source terms; no absolute force/work floor is inserted.
    pub maximum_relative_exchange_error: f64,
}

/// A refused interval publishes neither solid nor original source history.
#[derive(Debug)]
pub enum SolidMimoTransientError {
    Invalid(&'static str),
    Structure(ViscoElementError),
    Source(MimoLoadError),
    Exchange {
        quantity: &'static str,
        residual: f64,
        allowance: f64,
    },
    Cancelled,
}

impl core::fmt::Display for SolidMimoTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SOLID-MIMO-TRANSIENT: {self:?}")
    }
}
impl std::error::Error for SolidMimoTransientError {}
impl From<ViscoElementError> for SolidMimoTransientError {
    fn from(error: ViscoElementError) -> Self {
        Self::Structure(error)
    }
}
impl From<MimoLoadError> for SolidMimoTransientError {
    fn from(error: MimoLoadError) -> Self {
        Self::Source(error)
    }
}
type Result<T> = core::result::Result<T, SolidMimoTransientError>;

fn checkpoint(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| SolidMimoTransientError::Cancelled)
}

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolidMimoTransientError::Invalid(
            "unrepresentable original source exchange",
        ))
    }
}

fn check_exchange(quantity: &'static str, residual: f64, scale: f64, tolerance: f64) -> Result<()> {
    let allowance = finite(tolerance * finite(scale)?)?;
    let residual = finite(residual)?;
    if residual.abs() > allowance {
        Err(SolidMimoTransientError::Exchange {
            quantity,
            residual,
            allowance,
        })
    } else {
        Ok(())
    }
}

/// Original native step and independently checked rational boundary exchange.
pub struct SolidMimoTransientStep {
    /// Complete original constitutive forces, stored energy, physical heat and
    /// native port observations. `structure.body.external_work_j` is NET nodal
    /// work: applied work minus the resisting boundary's work.
    pub structure: ViscoDynamicPortStep,
    /// Exact original MIMO mean effort dotted with its held physical flow,
    /// times the actual source duration [J]. Existing history can release work.
    pub source_work_j: f64,
    /// Source work minus the independently computed native boundary work [J].
    pub source_work_residual_j: f64,
    /// Applied work plus the native support/eigenstrain work, minus the native
    /// kinetic/stored-energy changes, material heat and original source work [J].
    pub balance_residual_j: f64,
}

/// One immutable spatial connection, accepted solid state and original load.
///
/// Attachment consumes existing states without resetting material memory,
/// displacement, velocity, elapsed time, accepted step count or MIMO history.
/// The source's physical clock and continuous/sampled reports are retained.
/// Stationary supports and the original material-source domain remain those of
/// the native step. Source-band validity for an authored transient and physical
/// radiation qualification remain separate from this discrete exchange check.
pub struct PreparedSolidMimoTransient<'p, 'a> {
    trajectory: ViscoDynamicState<'p, 'a>,
    source: PreparedMimoLoad,
    ports: ViscoDynamicPortOperator<'p, 'a>,
    exchange_tolerance: f64,
}

impl<'p, 'a> PreparedSolidMimoTransient<'p, 'a> {
    /// Attach every original source port through a real row-major
    /// `free_dofs * source.ports()` force map in `equilibrium().free_dofs`
    /// order. Flow is its transpose times the actual midpoint nodal velocity;
    /// the structure receives the opposing mapped mean source effort.
    ///
    /// Original sparse inverse actions and the small port factor are prepared
    /// once. Each interval needs one ordinary sparse structural solve and a
    /// small Schur solve, not a dense spatial matrix or one solve per port.
    /// Full mutual, zero-direct, singular-direct and signed port maps survive.
    pub fn new(
        trajectory: ViscoDynamicState<'p, 'a>,
        source: PreparedMimoLoad,
        force_map: &[f64],
        budget: SolidMimoTransientBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        checkpoint(cx)?;
        let operator = trajectory.operator();
        if operator.duration_s().to_bits() != source.time_step_s().to_bits() {
            return Err(SolidMimoTransientError::Invalid(
                "solid and original MIMO load require exactly the same physical clock",
            ));
        }
        let bytes = source.limits().max_memory_bytes as u128 + 2048 + 32 * source.ports() as u128;
        if bytes > budget.maximum_load_payload_bytes as u128
            || !budget.maximum_relative_exchange_error.is_finite()
            || budget.maximum_relative_exchange_error <= 0.0
            || budget.maximum_relative_exchange_error > 1e-8
        {
            return Err(SolidMimoTransientError::Invalid(
                "complete original-load candidate allowance and finite exchange tolerance required",
            ));
        }
        let ports = operator.prepare_port_impedance(
            source.ports(),
            force_map,
            source.mean_direct_term(),
            budget.ports,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(Self {
            trajectory,
            source,
            ports,
            exchange_tolerance: budget.maximum_relative_exchange_error,
        })
    }

    pub fn duration_s(&self) -> f64 {
        self.source.time_step_s()
    }

    /// Original accepted full displacement/velocity and material quadrature histories.
    pub fn trajectory(&self) -> &ViscoDynamicState<'p, 'a> {
        &self.trajectory
    }

    /// Original admitted source, including its unchanged band and accepted pole states.
    pub fn source(&self) -> &PreparedMimoLoad {
        &self.source
    }

    /// Original continuous and sampled classifications; GridOnly remains GridOnly.
    pub fn source_reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        self.source.reports()
    }

    pub fn connection(&self) -> &ViscoDynamicPortOperator<'p, 'a> {
        &self.ports
    }

    /// Recover both accepted original owners, e.g. for a native impulse or
    /// explicit material continuation before preparing its new connection.
    pub fn into_parts(self) -> (ViscoDynamicState<'p, 'a>, PreparedMimoLoad) {
        (self.trajectory, self.source)
    }

    /// Advance one full node-major mean applied force. All physical states and
    /// native work/heat accounts commit together only after the source returns
    /// the same mean effort used in the mechanics and both exchange checks pass.
    pub fn advance(&mut self, mean_force_n: &[f64], cx: &Cx<'_>) -> Result<SolidMimoTransientStep> {
        checkpoint(cx)?;
        let p = self.source.ports();
        // The source and native port owners both admit at most 128 ports.
        let mut history = [0.0; 128];
        self.source.mean_history_output_into(&mut history[..p])?;
        let mut body_candidate = self.trajectory.clone();
        let structure = body_candidate.advance_with_port_impedance(
            &self.ports,
            mean_force_n,
            &history[..p],
            cx,
        )?;
        checkpoint(cx)?;
        let duration = self.duration_s();
        let source_candidate = self.source.prepare_mean_step(&structure.port.mean_flows)?;
        let effort = source_candidate.mean_effort();
        let mut source_work_j = 0.0;
        let mut work_correction = 0.0;
        let mut absolute_work_j = 0.0;
        for a in 0..p {
            checkpoint(cx)?;
            let direct_scale = self.ports.mean_impedance()[a * p..(a + 1) * p]
                .iter()
                .zip(&structure.port.mean_flows)
                .map(|(d, u)| (d * u).abs())
                .sum::<f64>();
            check_exchange(
                "original source mean effort",
                effort[a] - structure.port.mean_efforts[a],
                history[a].abs() + direct_scale + structure.port.mean_efforts[a].abs(),
                self.exchange_tolerance,
            )?;
            let work = finite(duration * structure.port.mean_flows[a] * effort[a])?;
            // Coherent ports may exchange large opposite works. Preserve the
            // remaining small transfer instead of erasing it by cancellation.
            let (sum, error) = fs_math::eft::two_sum(source_work_j, work);
            source_work_j = sum;
            work_correction += error;
            absolute_work_j += work.abs();
        }
        let source_work_j = finite(source_work_j + work_correction)?;
        let source_work_residual_j = finite(source_work_j - structure.port.boundary_work_j)?;
        check_exchange(
            "original source boundary work",
            source_work_residual_j,
            absolute_work_j + structure.port.boundary_work_j.abs(),
            self.exchange_tolerance,
        )?;
        let exchange =
            finite(structure.port.applied_work_j - structure.body.external_work_j - source_work_j)?;
        check_exchange(
            "applied nodal versus net solid and source work",
            exchange,
            structure.port.applied_work_j.abs()
                + structure.body.external_work_j.abs()
                + absolute_work_j,
            self.exchange_tolerance,
        )?;
        let balance_residual_j = finite(structure.body.balance_residual_j + exchange)?;
        checkpoint(cx)?;
        // Every fallible operation is complete. The exclusive source preview
        // now publishes its already checked pole state without a payload copy.
        source_candidate.commit();
        self.trajectory = body_candidate;
        Ok(SolidMimoTransientStep {
            structure,
            source_work_j,
            source_work_residual_j,
            balance_residual_j,
        })
    }
}

#[cfg(test)]
mod tests;
