//! Continuous physical forcing through the existing complete-coordinate owner.

use super::{
    ExactLinearSolid, ExactLinearSolidControls, ExactLinearSolidError, ExactLinearSolidStep,
    invalid, poll,
};
use crate::linear3::TetElasticAssembly;
use fs_exec::Cx;
use fs_phs::PortHamiltonian;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, PreparedExactTones,
};
use fs_sparse::Csr;

pub(super) enum SolidFlow {
    Polynomial(PreparedQuadraticFlow),
    Harmonic(HarmonicFlow),
}

pub(super) struct HarmonicFlow {
    plan: Box<PreparedExactTones<'static>>,
    average_output: Vec<f64>,
}

impl SolidFlow {
    pub(super) fn prepare_harmonic(
        system: PortHamiltonian,
        frequencies: &[f64],
        dt_s: f64,
        split_budget: LinearSplitBudget,
        flow_budget: DenseExponentialBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ExactLinearSolidError> {
        let ports = system.m;
        let plan = PreparedExactTones::new_owned_with_resistance_channels_controlled(
            system,
            frequencies,
            dt_s,
            split_budget,
            flow_budget,
            &[],
            || cx.checkpoint().is_err(),
        )?;
        Ok(Self::Harmonic(HarmonicFlow {
            plan: Box::new(plan),
            average_output: vec![0.0; ports],
        }))
    }

    pub(super) fn into_polynomial(self) -> Result<PreparedQuadraticFlow, ExactLinearSolidError> {
        match self {
            Self::Polynomial(flow) => Ok(flow),
            Self::Harmonic(_) => Err(invalid("held/affine-force preparation required")),
        }
    }

    fn into_harmonic(self) -> Result<HarmonicFlow, ExactLinearSolidError> {
        match self {
            Self::Harmonic(flow) => Ok(flow),
            Self::Polynomial(_) => Err(invalid("harmonic-force preparation required")),
        }
    }
}

// Match the shared tone owner's source, flow and independent power dimensions
// before constructing any physical factor. Reserve its numerical-payload bound
// twice, including the original M/K/C, inverse, pHS matrices and adapter state.
// Repeated frequencies use this conservative caller-count admission; the native
// owner combines them coherently before constructing its actual source map.
pub(super) fn preparation_shape(
    coordinates: usize,
    tones: usize,
) -> Result<(usize, usize), ExactLinearSolidError> {
    let state = coordinates
        .checked_mul(2)
        .ok_or_else(|| invalid("harmonic state dimension overflow"))?;
    let source = state
        .checked_mul(tones)
        .ok_or_else(|| invalid("harmonic source dimension overflow"))?;
    let flow = state
        .checked_add(source)
        .and_then(|v| v.checked_add(1))
        .ok_or_else(|| invalid("harmonic flow dimension overflow"))?;
    let power = state
        .checked_mul(2)
        .and_then(|v| v.checked_add(source))
        .ok_or_else(|| invalid("harmonic power dimension overflow"))?;
    let dimension = flow.max(power);
    if dimension > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
        return Err(ExactLinearSolidError::Budget {
            required: dimension,
            maximum: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        });
    }
    let native_bytes = 8 * (flow * flow * (2 * coordinates + 25) + power * power * 80)
        + 128 * (state + coordinates + source + tones);
    let bytes = 2 * native_bytes + if coordinates > 32 { 2 * 1024 * 1024 } else { 0 };
    Ok((dimension, bytes))
}

/// Complete structural motion under caller-phased continuous tones.
///
/// Physical source admission, original matrix energy and transactional state
/// acceptance are shared with `ExactLinearSolid`. The distinct owner retains
/// the generic tone service's original storage without changing the existing
/// held/affine owner's `Send + Sync` numerical representation.
pub struct ExactHarmonicSolid {
    inner: ExactLinearSolid<HarmonicFlow>,
}

impl ExactHarmonicSolid {
    /// Prepare continuous harmonic loads on every original free coordinate.
    ///
    /// Frequencies are finite nonnegative angular frequencies [rad/s]; zero is
    /// an explicit DC source. They need not be integer multiples or distinct.
    /// All original mass/stiffness entries and initial transients remain. Exact
    /// resonance retains its growing transient; no steady-state solution,
    /// artificial damping, modal truncation or Nyquist retuning is substituted.
    ///
    /// The shared complete exponential/Gramian owner prepares the motion and
    /// independent signed work once. Its flow/power dimension and numerical
    /// payload are admitted before physical factorization under `controls`.
    /// This remains a bounded dense owner, not a sparse-scale harmonic solver.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        assembly: &TetElasticAssembly,
        angular_frequencies_rad_s: &[f64],
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        controls: ExactLinearSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, ExactLinearSolidError> {
        let source = ExactLinearSolid::prepare_source(
            assembly,
            dt_s,
            initial_displacement_m,
            initial_velocity_mps,
            None,
            None,
            Some(angular_frequencies_rad_s),
            controls,
            cx,
        )?;
        Ok(Self {
            inner: source.change_flow(SolidFlow::into_harmonic)?,
        })
    }

    /// Continuous tones with the same full physical symmetric PSD viscosity C
    /// [N s/m] as `ExactLinearSolid::prepare_viscous`. Nonproportional damping and consistent mass
    /// are not diagonalized away. Native original resistance-square integrals
    /// account for physical heat independently of force work and stored energy.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_viscous(
        assembly: &TetElasticAssembly,
        damping_ns_m: &Csr,
        angular_frequencies_rad_s: &[f64],
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        controls: ExactLinearSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self, ExactLinearSolidError> {
        let source = ExactLinearSolid::prepare_source(
            assembly,
            dt_s,
            initial_displacement_m,
            initial_velocity_mps,
            None,
            Some(damping_ns_m),
            Some(angular_frequencies_rad_s),
            controls,
            cx,
        )?;
        Ok(Self {
            inner: source.change_flow(SolidFlow::into_harmonic)?,
        })
    }

    /// Caller-ordered angular frequencies [rad/s].
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.inner.flow.plan.angular_frequencies_rad_s()
    }

    pub fn free_dofs(&self) -> &[usize] {
        self.inner.free_dofs()
    }
    pub fn displacement_m(&self) -> &[f64] {
        self.inner.displacement_m()
    }
    pub fn velocity_mps(&self) -> &[f64] {
        self.inner.velocity_mps()
    }
    pub fn time_s(&self) -> f64 {
        self.inner.time_s()
    }
    pub fn duration_s(&self) -> f64 {
        self.inner.duration_s()
    }
    pub fn accepted_steps(&self) -> u64 {
        self.inner.accepted_steps()
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.inner.stored_energy_j()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.inner.kinetic_energy_j()
    }
    pub fn cumulative_work_j(&self) -> f64 {
        self.inner.cumulative_work_j()
    }
    pub fn cumulative_dissipated_energy_j(&self) -> f64 {
        self.inner.cumulative_dissipated_energy_j()
    }

    /// Rotate one source coefficient pair to the next interval's local origin.
    ///
    /// The convention is `Fc*cos(omega*tau) + Fs*sin(omega*tau)` (equivalently
    /// peak phasor `Fc+i*Fs` under `exp(-i*omega*tau)`). Use the returned pair
    /// only after an accepted step. This uses cached sine/cosine coefficients;
    /// it neither advances the accepted clock nor mutates caller source phase.
    /// A checkpoint must retain source coefficients together with physical state.
    pub fn advance_force_phase(
        &self,
        tone: usize,
        cosine_n: f64,
        sine_n: f64,
    ) -> Result<(f64, f64), ExactLinearSolidError> {
        Ok(self
            .inner
            .flow
            .plan
            .advance_force_phase(tone, cosine_n, sine_n)?)
    }

    /// Advance the complete coupled trajectory under simultaneous continuous tones.
    ///
    /// Arrays are free-coordinate-major: `i*K+j` is tone j on free DOF i, in N.
    /// Their cosine/sine phases refer to THIS interval's start (`tau=0`). A
    /// continuing source uses `advance_force_phase` after success; an
    /// explicitly changed next coefficient pair represents a finite force jump,
    /// not an impulse or a phase inferred from the global clock.
    ///
    /// Signed work and physical viscosity include every cross-frequency term
    /// and the initial transient. They are not sums of isolated cycle averages.
    /// Original M/K energy and per-step/cumulative balance gate acceptance.
    /// Warm calls reuse all exponentials and work/heat factors; only O(N) result
    /// vectors are allocated. Cancellation is polled before and after the bounded
    /// native contraction. Refusal preserves accepted state, clock and accounts.
    pub fn advance(
        &mut self,
        cosine_force_n: &[f64],
        sine_force_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep, ExactLinearSolidError> {
        poll(cx)?;
        let (step, time_s) = self.inner.next_interval()?;
        let mut candidate = vec![0.0; self.inner.state.len()];
        let record = self.inner.flow.plan.step_into(
            &self.inner.state,
            cosine_force_n,
            sine_force_n,
            &mut candidate,
            &mut self.inner.flow.average_output,
        )?;
        self.inner.accept(
            candidate,
            record.supplied,
            record.dissipated,
            step,
            time_s,
            cx,
        )
    }
}
