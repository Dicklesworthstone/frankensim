//! Reusable physical source, periodic tangent and carrier factor ownership.

use super::super::{
    forced_problem,
    tangent::StorageTangent,
    variational_feedback::{FeedbackVariationalProblem, feedback_tangent_problem},
};
use super::{
    AutonomousPortFeedback, C64, HbOrbit, PeriodicPortInput, PhsFeedbackStateJacobian,
    PhsOrbitProblem, PhsPortSidebandControls, PhsPortSidebandError, PhsPortSidebandResponse,
    PhsStorageHessianAction, PortHamiltonian, Result, SidebandResponse, SidebandScales, admit,
    admit_carrier, admit_forcing, admit_source, admit_storage, input_forcing, observe, poll,
};
use fs_orbit::{
    IslandStateJacobian, OrbitProblem,
    sideband::{
        PreparedPeriodicLinearization, PreparedSidebandResponse,
        refinement::{
            SidebandRefinementCheck, SidebandRefinementControls, solve_sideband_adaptive,
        },
    },
};

mod adaptive;
mod frequency;
mod parameter;
mod resonance;
pub(super) mod sensitivity;
mod stability;
pub use adaptive::{
    AdaptivePhsPortResonanceResponse, AdaptivePhsPortSidebandCarrierDerivativeResponse,
    PhsPortAdaptiveResonanceRequest, PhsPortResonanceRefinementCheck,
};
pub use frequency::PhsForcedFrequencyTangent;
pub use parameter::{
    PhsForcedParameterDirection, PhsForcedParameterTangent, PhsForcedStorageParameterTangent,
    PhsParameterControls, PhsStorageParameterControls,
};
pub use resonance::{PhsPortResonanceRequest, PhsPortResonanceResponse};
pub use sensitivity::PhsPortSidebandCarrierDerivativeResponse;
pub use stability::{
    PhsBaseStabilityControls, PhsBaseStabilityError, PhsBaseStabilityLevel, PhsBaseStabilityReport,
    PhsBaseTransientControls, PhsBaseTransientError, PhsBaseTransientLevel, PhsBaseTransientReport,
};

/// A complete physical response with independently measured sideband convergence.
/// The state order is `response.state_harmonics`; the output order remains the
/// separately authored `response.output_harmonics`. Neither finite-grid trace
/// is an infinite-tail, forward-error or Floquet-stability certificate.
#[derive(Clone, Debug)]
pub struct AdaptivePhsPortSidebandResponse {
    pub response: PhsPortSidebandResponse,
    /// Every actually solved state order and its original-equation assessment.
    pub refinement_trace: Vec<SidebandRefinementCheck>,
}

// This dispatch only selects existing ORIGINAL source owners. In particular,
// feedback uses its total derivative and never the storage-only Jacobian.
enum SourceProblem<'a> {
    Forced(PhsOrbitProblem<'a>),
    Feedback(FeedbackVariationalProblem<'a>),
}
impl SourceProblem<'_> {
    fn original(&self) -> &PhsOrbitProblem<'_> {
        match self {
            Self::Forced(problem) => problem,
            Self::Feedback(problem) => problem.original_problem(),
        }
    }

    fn state_jacobian(&self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        match self {
            Self::Forced(problem) => problem
                .island_state_jacobian()
                .is_some_and(|action| action(problem, time, state, out)),
            Self::Feedback(problem) => problem
                .island_state_jacobian()
                .is_some_and(|action| action(problem, time, state, out)),
        }
    }
}
impl OrbitProblem for SourceProblem<'_> {
    fn dim(&self) -> usize {
        self.original().dim()
    }
    fn autonomous(&self) -> bool {
        match self {
            Self::Forced(problem) => problem.autonomous(),
            Self::Feedback(problem) => problem.autonomous(),
        }
    }
    fn port(&self, frequency: C64) -> Vec<C64> {
        match self {
            Self::Forced(problem) => problem.port(frequency),
            Self::Feedback(problem) => problem.port(frequency),
        }
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        match self {
            Self::Forced(problem) => problem.island(time, state, out),
            Self::Feedback(problem) => problem.island(time, state, out),
        }
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.original()
            .tangent
            .available()
            .then_some(Self::state_jacobian)
    }
}

/// Original physical source and base waveform for repeated resonance probes.
///
/// Constructors bind the immutable original system, input/feedback, complete
/// storage tangent and base orbit. All material and controller memory stays in
/// its original state coordinates. Sources and authored callbacks must remain
/// deterministic and unchanged for this owner's lifetime, including through
/// interior mutability. Changing the model or base requires a new source.
///
/// Call [`Self::prepare`] once, then factor each signed carrier once and solve
/// as many physical excitations as needed. The owners borrow each other; there
/// is no self-reference, copied physical model or inferred output law. Mutable
/// preparation/factor borrows allow one live cached tangent and one live factor
/// per source. Every solve returns a complete independently observed response.
///
/// Controls are copied and admitted before adapter scratch is allocated. The
/// generic payload cap and adapter payload cap are independent: their SUM
/// bounds live library numerical storage, including observation while the
/// tangent and factor remain cached. Work/action limits apply per operation,
/// not to an unbounded sequence of calls. Previously returned responses are
/// caller storage and are not charged to later operations.
pub struct PhsPortSidebandSource<'model> {
    problem: SourceProblem<'model>,
    orbit: &'model HbOrbit,
    controls: PhsPortSidebandControls,
}

impl<'model> PhsPortSidebandSource<'model> {
    /// Bind a periodically driven source using its complete native storage tangent.
    /// # Errors
    /// Invalid original input/orbit, unresolved forcing, or adapter resource caps.
    pub fn forced(
        system: &'model PortHamiltonian,
        input: &'model PeriodicPortInput,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
    ) -> Result<Self> {
        Self::from_forced(system, input, orbit, controls, None)
    }

    /// Bind the same driven source with an authored ORIGINAL complete Hessian.
    /// The action supplies both vector-field tangents and physical observation;
    /// refusal never substitutes a native or finite-difference action.
    /// # Errors
    /// The native constructor's admission failures; action failures occur at use.
    pub fn forced_analytic(
        system: &'model PortHamiltonian,
        input: &'model PeriodicPortInput,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
        hessian: &'model PhsStorageHessianAction<'model>,
    ) -> Result<Self> {
        Self::from_forced(system, input, orbit, controls, Some(hessian))
    }

    /// Bind a self-excited source with its TOTAL original feedback derivative.
    /// The derivative includes the physical flow/Hessian chain exactly once.
    /// # Errors
    /// Invalid complete ports/base, or adapter/tangent-scratch resource caps.
    pub fn autonomous_feedback(
        system: &'model PortHamiltonian,
        feedback: &'model dyn AutonomousPortFeedback,
        feedback_jacobian: &'model PhsFeedbackStateJacobian<'model>,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
    ) -> Result<Self> {
        Self::from_feedback(
            system,
            None,
            feedback,
            feedback_jacobian,
            orbit,
            controls,
            None,
        )
    }

    /// Bind self-excited feedback and an authored complete original Hessian.
    /// # Errors
    /// The native feedback constructor's admission failures; no action fallback.
    pub fn autonomous_feedback_analytic(
        system: &'model PortHamiltonian,
        feedback: &'model dyn AutonomousPortFeedback,
        feedback_jacobian: &'model PhsFeedbackStateJacobian<'model>,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
        hessian: &'model PhsStorageHessianAction<'model>,
    ) -> Result<Self> {
        Self::from_feedback(
            system,
            None,
            feedback,
            feedback_jacobian,
            orbit,
            controls,
            Some(hessian),
        )
    }

    /// Bind periodic physical driving AND instantaneous nonlinear port feedback.
    /// The sustaining drive fixes the base clock; feedback supplies its TOTAL
    /// state derivative. Probe inputs are additional original-port efforts.
    /// # Errors
    /// Invalid complete ports/base, mismatched drive clock/bandwidth, or budgets.
    pub fn forced_feedback(
        system: &'model PortHamiltonian,
        input: &'model PeriodicPortInput,
        feedback: &'model dyn AutonomousPortFeedback,
        feedback_jacobian: &'model PhsFeedbackStateJacobian<'model>,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
    ) -> Result<Self> {
        Self::from_feedback(
            system,
            Some(input),
            feedback,
            feedback_jacobian,
            orbit,
            controls,
            None,
        )
    }

    /// Bind the driven feedback source with an authored complete storage Hessian.
    /// The feedback derivative must use this same original effort law.
    /// # Errors
    /// The native driven-feedback constructor's admission failures; no fallback.
    pub fn forced_feedback_analytic(
        system: &'model PortHamiltonian,
        input: &'model PeriodicPortInput,
        feedback: &'model dyn AutonomousPortFeedback,
        feedback_jacobian: &'model PhsFeedbackStateJacobian<'model>,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
        hessian: &'model PhsStorageHessianAction<'model>,
    ) -> Result<Self> {
        Self::from_feedback(
            system,
            Some(input),
            feedback,
            feedback_jacobian,
            orbit,
            controls,
            Some(hessian),
        )
    }

    pub(super) fn from_forced(
        system: &'model PortHamiltonian,
        input: &'model PeriodicPortInput,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
        hessian: Option<&'model PhsStorageHessianAction<'model>>,
    ) -> Result<Self> {
        admit_source(system, orbit, controls)?;
        admit_input(system, input, orbit, controls)?;
        let controls = source_controls(system, controls, false)?;
        let mut problem = forced_problem(system, input).map_err(|_| {
            PhsPortSidebandError::Invalid("complete original physical port input required")
        })?;
        if let Some(action) = hessian {
            problem.tangent = StorageTangent::hessian(system.state_dim(), action);
        }
        Ok(Self {
            problem: SourceProblem::Forced(problem),
            orbit,
            controls,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn from_feedback(
        system: &'model PortHamiltonian,
        input: Option<&'model PeriodicPortInput>,
        feedback: &'model dyn AutonomousPortFeedback,
        feedback_jacobian: &'model PhsFeedbackStateJacobian<'model>,
        orbit: &'model HbOrbit,
        controls: &PhsPortSidebandControls,
        hessian: Option<&'model PhsStorageHessianAction<'model>>,
    ) -> Result<Self> {
        admit_source(system, orbit, controls)?;
        if feedback.port_count() != system.port_dim() {
            return Err(PhsPortSidebandError::Invalid(
                "feedback must cover every original physical port",
            ));
        }
        if let Some(input) = input {
            admit_input(system, input, orbit, controls)?;
        }
        let controls = source_controls(system, controls, true)?;
        let mut problem = feedback_tangent_problem(system, feedback, feedback_jacobian, hessian)
            .map_err(|_| {
                PhsPortSidebandError::Invalid("complete original feedback source required")
            })?;
        problem.external = input;
        Ok(Self {
            problem: SourceProblem::Feedback(problem),
            orbit,
            controls,
        })
    }

    /// Verify the complete base balance and cache its full periodic tangent once.
    /// The same native generic preparation owns all sampling and scaling checks.
    /// An exclusive borrow prevents another live preparation for this source.
    /// # Errors
    /// Generic admission/budget, absent/refused tangents, unresolved source
    /// balance or tangent sampling, and cancellation. No partial cache escapes.
    ///
    /// A second preparation requires dropping the first, even if its last
    /// calculation has finished; its numerical allocations are still live.
    /// ```compile_fail
    /// use fs_couple::phs_orbit::sideband::PhsPortSidebandSource;
    /// use fs_orbit::sideband::SidebandScales;
    /// fn overlapping(source: &mut PhsPortSidebandSource<'_>, scales: &SidebandScales) {
    ///     let _first = source.prepare(scales, || false).unwrap();
    ///     let _second = source.prepare(scales, || false).unwrap();
    /// }
    /// ```
    pub fn prepare(
        &mut self,
        scales: &SidebandScales,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsPortSidebands<'_, 'model>> {
        let linearization = PreparedPeriodicLinearization::new(
            &self.problem,
            self.orbit,
            self.controls.sideband,
            scales,
            cancelled,
        )?;
        Ok(PreparedPhsPortSidebands {
            source: self,
            linearization,
        })
    }

    /// Resolve the complete physical response by growing the retained state order.
    ///
    /// `forcing` contains all original-port efforts at the SOURCE'S initial
    /// sideband order. Map G once, then let the shared generic owner preserve
    /// that exact support while solving and independently assessing larger full
    /// state systems. A small change between two truncated solves is insufficient
    /// without the original-equation and omitted-band checks.
    ///
    /// After state acceptance, observe the ENTIRE accepted state response through
    /// the original `G^T Hxx` action. The output harmonic order stays independently
    /// authored; it is never silently raised or replaced by the state order.
    /// Observation grids grow enough to resolve the accepted state, within the
    /// unchanged maximum grid, Hessian-action, work and payload allowances.
    ///
    /// The generic `maximum_work` is cumulative across every attempted order in
    /// this call. Adapter work covers one physical input map and the final output
    /// assessment. Cloned tolerances, mapped input and retained refinement trace
    /// are admitted explicitly; the generic and adapter live payload caps add.
    /// This method recomputes the adaptive solve; use [`Self::prepare`] to reuse
    /// one chosen fixed-order tangent/factor across subsequent drives.
    ///
    /// # Errors
    /// Invalid original forcing/scales/refinement, generic resolution or numerical
    /// refusal, insufficient final observation resources, unresolved physical
    /// output, source refusal or cancellation. A generic state solution alone
    /// never escapes as an accepted physical response. Source controls/base are
    /// unchanged after success or refusal, so the same source can be reused.
    pub fn solve_adaptive(
        &mut self,
        carrier_angular_frequency: f64,
        forcing: &[C64],
        scales: &SidebandScales,
        refinement: &SidebandRefinementControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptivePhsPortSidebandResponse> {
        poll(&mut cancelled)?;
        let system = self.problem.original().system;
        admit_carrier(self.orbit, carrier_angular_frequency, &self.controls)?;
        admit_forcing(system, forcing, &self.controls)?;
        // The generic 0..=64 order ladder has at most eight checks, including
        // both endpoints. Keep that returned history alive during observation,
        // and charge this second owned tolerance vector before either allocation.
        let extra = 8 * system.port_dim() as u128
            + 8 * core::mem::size_of::<SidebandRefinementCheck>() as u128;
        if extra > self.controls.max_workspace_bytes as u128 {
            return Err(PhsPortSidebandError::Budget(
                "adaptive physical controls and trace exceed the adapter payload cap",
            ));
        }
        let mut observer = self.controls.clone();
        observer.max_workspace_bytes -= extra as usize;
        // The existing 32*n*(2*N+1) term includes both mapped input and response;
        // the mapped input is dropped before any output-grid allocation.
        admit_storage(
            system.state_dim(),
            system.port_dim(),
            &observer,
            observer.initial_samples,
        )?;
        let resolved = {
            let rhs = input_forcing(system, forcing, observer.sideband.harmonics, &mut cancelled)?;
            solve_sideband_adaptive(
                &self.problem,
                self.orbit,
                carrier_angular_frequency,
                &rhs,
                observer.sideband,
                scales,
                refinement,
                &mut cancelled,
            )?
        };
        observer.sideband.harmonics = resolved.harmonics;
        let largest = (self.orbit.coeffs[0].len() - 1)
            .max(resolved.harmonics)
            .max(observer.output_harmonics);
        // These extents were bounded by source and generic admission. The FFT
        // grid must be strictly larger than twice every synthesized order.
        observer.initial_samples = observer
            .initial_samples
            .max((2 * largest + 1).next_power_of_two());
        if observer.initial_samples > observer.max_samples / 4 {
            return Err(PhsPortSidebandError::Budget(
                "refined physical response needs three resolved observation grids within the maximum",
            ));
        }
        admit_source(system, self.orbit, &observer)?;
        let response = self.observe_response(
            carrier_angular_frequency,
            resolved.response,
            &observer,
            resolved.base_residual,
            resolved.tangent_sampling_defect,
            &mut cancelled,
        )?;
        Ok(AdaptivePhsPortSidebandResponse {
            response,
            refinement_trace: resolved.refinement_trace,
        })
    }

    fn observe_response(
        &self,
        carrier: f64,
        state_response: SidebandResponse,
        controls: &PhsPortSidebandControls,
        base_residual: f64,
        tangent_sampling_defect: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhsPortSidebandResponse> {
        let original = self.problem.original();
        let (port_coefficients, observation_trace, hessian_actions) = observe(
            original,
            self.orbit,
            &state_response.coefficients,
            controls,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(PhsPortSidebandResponse {
            omega_rad_s: self.orbit.omega,
            carrier_angular_frequency: carrier,
            state_harmonics: controls.sideband.harmonics,
            output_harmonics: controls.output_harmonics,
            state_dimension: original.system.state_dim(),
            port_count: original.system.port_dim(),
            state_response,
            port_coefficients,
            base_residual,
            tangent_sampling_defect,
            observation_trace,
            hessian_actions,
        })
    }
}

fn admit_input(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    orbit: &HbOrbit,
    controls: &PhsPortSidebandControls,
) -> Result<()> {
    if system.port_dim() != input.port_count()
        || orbit.omega.to_bits() != input.omega_rad_s().to_bits()
        || controls.sideband.samples_per_period as u128
            <= 2 * input.highest_nonzero_harmonic() as u128
    {
        return Err(PhsPortSidebandError::Invalid(
            "matching complete resolved periodic input required",
        ));
    }
    Ok(())
}

fn source_controls(
    system: &PortHamiltonian,
    controls: &PhsPortSidebandControls,
    feedback: bool,
) -> Result<PhsPortSidebandControls> {
    // Charge the owned tolerances, plus complete feedback derivative/flow
    // scratch when present, before cloning controls or constructing the source.
    let columns = if feedback {
        system.state_dim() as u128 + 2
    } else {
        1
    };
    let extra = 8 * system.port_dim() as u128 * columns;
    if extra > controls.max_workspace_bytes as u128 {
        return Err(PhsPortSidebandError::Budget(
            "physical source scratch exceeds the adapter payload cap",
        ));
    }
    let mut remaining = controls.clone();
    remaining.max_workspace_bytes -= extra as usize;
    admit_storage(
        system.state_dim(),
        system.port_dim(),
        &remaining,
        remaining.initial_samples,
    )?;
    Ok(remaining)
}

/// One independently checked complete periodic tangent, reusable across carriers.
/// This borrows a separate physical source; it cannot outlive or replace its base.
pub struct PreparedPhsPortSidebands<'source, 'model> {
    source: &'source PhsPortSidebandSource<'model>,
    linearization: PreparedPeriodicLinearization<'source, SourceProblem<'model>>,
}

// Explicit Drop keeps the source borrow live until its cached allocations are
// actually destroyed. Without it, NLL can end reference borrows at last use
// while sibling Vec allocations still exist, admitting another live cache.
impl Drop for PreparedPhsPortSidebands<'_, '_> {
    fn drop(&mut self) {}
}

impl<'model> PreparedPhsPortSidebands<'_, 'model> {
    /// Complete original base-equation residual measured on both preparation grids.
    pub fn base_residual(&self) -> f64 {
        self.linearization.base_residual
    }
    /// Independently sampled difference of the complete normalized tangent arrays.
    pub fn tangent_sampling_defect(&self) -> f64 {
        self.linearization.tangent_sampling_defect
    }

    /// Factor one signed carrier using the cached full periodic tangent.
    /// The exclusive borrow permits only one live dense factor per preparation.
    /// No base or tangent callbacks are repeated. Signed shifted frequencies,
    /// singularities and residual semantics remain those of fs-orbit.
    /// # Errors
    /// Invalid physical frequencies, singular original response, source-port
    /// refusal, generic numerical failure or cancellation. No damping floor.
    ///
    /// Drop one factor before asking for another; merely finishing its last
    /// solve does not release the retained operator and LU payload.
    /// ```compile_fail
    /// use fs_couple::phs_orbit::sideband::PreparedPhsPortSidebands;
    /// fn overlapping(prepared: &mut PreparedPhsPortSidebands<'_, '_>) {
    ///     let _first = prepared.factor(0.3, || false).unwrap();
    ///     let _second = prepared.factor(0.7, || false).unwrap();
    /// }
    /// ```
    pub fn factor(
        &mut self,
        carrier_angular_frequency: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FactoredPhsPortSidebands<'_, 'model>> {
        poll(&mut cancelled)?;
        admit_carrier(
            self.source.orbit,
            carrier_angular_frequency,
            &self.source.controls,
        )?;
        let base_residual = self.base_residual();
        let tangent_sampling_defect = self.tangent_sampling_defect();
        let factor = self
            .linearization
            .factor(carrier_angular_frequency, cancelled)?;
        Ok(FactoredPhsPortSidebands {
            source: self.source,
            factor,
            base_residual,
            tangent_sampling_defect,
        })
    }
}

/// One carrier's original coupled operator/LU, reusable across physical drives.
/// Each solve independently observes `G^T Hxx(x_base) dx` with the same source.
pub struct FactoredPhsPortSidebands<'prepared, 'model> {
    source: &'prepared PhsPortSidebandSource<'model>,
    factor: PreparedSidebandResponse<'prepared>,
    base_residual: f64,
    tangent_sampling_defect: f64,
}

// Retain the preparation borrow until the dense operator and factor are
// destroyed, including when the caller no longer reads this factor's fields.
impl Drop for FactoredPhsPortSidebands<'_, '_> {
    fn drop(&mut self) {}
}

impl FactoredPhsPortSidebands<'_, '_> {
    /// Signed probe carrier in radians per unit time, on the original base clock.
    pub fn carrier_angular_frequency(&self) -> f64 {
        self.factor.carrier_angular_frequency()
    }

    /// Solve a new physical-port excitation without resampling A(t) or refactoring.
    /// `forcing[(k+N)*ports+p]` is the peak complex effort at the original port.
    /// Observation still uses the complete original Hessian; it is not inferred
    /// from the cached vector-field tangent or copied from state coefficients.
    /// Per-response work/actions and independent output refinement start afresh.
    /// # Errors
    /// Malformed forcing, original input/output refusal, unresolved response or
    /// output, exhausted per-operation budgets, or cancellation. No partial
    /// physical result escapes; failure does not consume or change the factor.
    pub fn solve(
        &mut self,
        forcing: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsPortSidebandResponse> {
        poll(&mut cancelled)?;
        let original = self.source.problem.original();
        let system = original.system;
        let controls = &self.source.controls;
        admit_forcing(system, forcing, controls)?;
        let state_response = {
            let rhs = input_forcing(system, forcing, controls.sideband.harmonics, &mut cancelled)?;
            self.factor.solve(&rhs, &mut cancelled)?
        };
        self.source.observe_response(
            self.carrier_angular_frequency(),
            state_response,
            controls,
            self.base_residual,
            self.tangent_sampling_defect,
            &mut cancelled,
        )
    }
}

/// External physical probes about a periodically driven nonlinear feedback loop.
///
/// Base balance uses `(J-R) grad H + G (u_periodic + u_feedback)`; its exact
/// tangent uses `(J-R) Hxx + G du_total/dx`. The authored feedback derivative
/// includes the original flow/Hessian chain exactly once. The sustaining drive
/// and all original memory coordinates remain present; the additional probe
/// enters the original G ports. Outputs use the same original storage Hessian.
/// The external drive fixes the frequency/phase, with no autonomous phase pin.
///
/// # Errors
/// Invalid complete source/clock, unavailable or refusing exact derivatives,
/// unresolved base/response/output, budgets, singularity, or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_feedback_port_sideband_response(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    orbit: &HbOrbit,
    carrier_angular_frequency: f64,
    forcing: &[C64],
    controls: &PhsPortSidebandControls,
    scales: &SidebandScales,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    solve_driven_feedback(
        system,
        input,
        feedback,
        feedback_jacobian,
        orbit,
        carrier_angular_frequency,
        forcing,
        controls,
        scales,
        None,
        cancelled,
    )
}

/// The driven-feedback probe with an authored complete ORIGINAL storage Hessian.
/// Both vector-field and physical-output tangents use this action, while the
/// TOTAL feedback derivative includes the chain through the same effort law.
/// # Errors
/// The native driven-feedback wrapper's failures; declared actions never fall back.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_feedback_port_sideband_response_analytic(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    orbit: &HbOrbit,
    carrier_angular_frequency: f64,
    forcing: &[C64],
    controls: &PhsPortSidebandControls,
    scales: &SidebandScales,
    hessian: &PhsStorageHessianAction<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    solve_driven_feedback(
        system,
        input,
        feedback,
        feedback_jacobian,
        orbit,
        carrier_angular_frequency,
        forcing,
        controls,
        scales,
        Some(hessian),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn solve_driven_feedback(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    orbit: &HbOrbit,
    carrier: f64,
    forcing: &[C64],
    controls: &PhsPortSidebandControls,
    scales: &SidebandScales,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    poll(&mut cancelled)?;
    admit(system, orbit, carrier, forcing, controls)?;
    let mut source = PhsPortSidebandSource::from_feedback(
        system,
        Some(input),
        feedback,
        feedback_jacobian,
        orbit,
        controls,
        hessian,
    )?;
    let mut prepared = source.prepare(scales, &mut cancelled)?;
    let mut factor = prepared.factor(carrier, &mut cancelled)?;
    factor.solve(forcing, cancelled)
}
