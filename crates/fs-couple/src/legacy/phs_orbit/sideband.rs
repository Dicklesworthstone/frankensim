//! Physical small-probe response about a driven or self-excited nonlinear pHS orbit.
//!
//! fs-orbit solves the complete periodic tangent. This adapter supplies the
//! ORIGINAL port forcing `G du` and observes `dy = G^T Hxx(x_base) dx`.
//! Nonlinear stiffness, reciprocal port signs and every explicit loss-memory
//! coordinate remain in their original units and order. State coefficients
//! alone are not physical output coefficients when the storage is nonlinear.
//!
//! Input, state and output sidebands all use peak complex amplitudes in
//! `Re(exp(i*nu*t) sum_k Z_k exp(i*k*omega*t))`. No conjugation, doubling or
//! positive-frequency folding is performed. At coincident signed physical
//! frequencies contributions must be combined coherently by the consumer.
//! This is an infinitesimal-probe response, not a new finite-amplitude orbit
//! or a claim of Floquet stability. Its finite-grid checks are not tail bounds.

use super::{
    AutonomousPortFeedback, PeriodicPortInput, PhsFeedbackStateJacobian, PhsOrbitProblem,
    PhsStorageHessianAction,
};
use fs_fft::{C64 as FftComplex, Fft};
use fs_math::{c64::C64, det};
use fs_orbit::{
    HbOrbit,
    sideband::{SidebandControls, SidebandError, SidebandResponse, SidebandScales},
};
use fs_phs::PortHamiltonian;

mod prepared;
pub(super) use prepared::sensitivity::{admit_port_tolerances, observe_original_pair};
pub use prepared::{
    AdaptivePhsPortResonanceResponse, AdaptivePhsPortSidebandCarrierDerivativeResponse,
    AdaptivePhsPortSidebandResponse, FactoredPhsPortSidebands, PhsBaseStabilityControls,
    PhsBaseStabilityError, PhsBaseStabilityLevel, PhsBaseStabilityReport, PhsBaseTransientControls,
    PhsBaseTransientError, PhsBaseTransientLevel, PhsBaseTransientReport,
    PhsForcedFrequencyTangent, PhsForcedParameterDirection, PhsForcedParameterTangent,
    PhsForcedStorageParameterTangent, PhsParameterControls, PhsPortAdaptiveResonanceRequest,
    PhsPortResonanceRefinementCheck, PhsPortResonanceRequest, PhsPortResonanceResponse,
    PhsPortSidebandCarrierDerivativeResponse, PhsPortSidebandSource, PhsStorageParameterControls,
    PreparedPhsPortSidebands, solve_forced_feedback_port_sideband_response,
    solve_forced_feedback_port_sideband_response_analytic,
};

/// Complete derivative of the ORIGINAL storage effort with respect to one
/// caller-declared real parameter, at fixed physical state. Write every entry
/// of `partial_p grad H(x,p)` in original effort units per parameter unit.
/// Include all nonlinear and explicit material-memory terms. This is distinct
/// from the Hessian action, which differentiates the state. The callback must
/// describe the same immutable source, be deterministic/read-only, and refuse
/// unsupported states or unrepresentable arithmetic. No differencing fallback
/// is supplied for refusal, incomplete writes or nonfinite outputs.
pub type PhsStorageParameterAction<'a> = dyn Fn(&[f64], &mut [f64]) -> bool + 'a;

/// Independent state-response and physical-output accuracy/resource limits.
/// Work and action caps apply to each operation, not cumulatively across a
/// reusable source's lifetime. Returned responses become caller-owned storage.
#[derive(Clone, Debug)]
pub struct PhsPortSidebandControls {
    pub sideband: SidebandControls,
    /// Retain output indices -output_harmonics through +output_harmonics.
    /// This may exceed the state order: Hxx(t) can generate more output bands.
    pub output_harmonics: usize,
    /// First power-of-two observation grid, at least 32 and strictly above
    /// twice the base-state, response-state and output orders.
    pub initial_samples: usize,
    /// Largest power-of-two grid, at most 65536 and at least four times the
    /// first. Two successive grid comparisons must pass independently.
    pub max_samples: usize,
    /// Maximum observation Hessian actions per operation; two per sample, for
    /// the real and imaginary envelopes. A response-with-carrier-derivative
    /// operation shares this cap across BOTH observations. Generic tangent
    /// preparation has its own budget.
    pub max_hessian_actions: usize,
    /// Extra adapter numerical payload, including mapped input, retained state
    /// response, tangent scratch, FFT plans/samples and both output spectra.
    /// The generic solver's separate payload cap still applies. A prepared
    /// response keeps its tangent and one factor alive during observation:
    /// their combined live bound is the SUM of this cap and the generic cap.
    /// Original model/caller storage, callback-internal allocation and allocator
    /// metadata are excluded. Adapter source/tolerance scratch is included.
    pub max_workspace_bytes: usize,
    /// Conservative observer/port-map arithmetic cap per operation, including
    /// dense quadratic Hessian actions. A response/derivative pair shares this
    /// cap across both observations and its exact physical port derivative.
    /// Other source callbacks' internal work is source-defined.
    pub max_observation_work: u64,
    /// Absolute complex-envelope L2 tolerance per original output port, in its
    /// physical units. Retained sampling change and observed omitted-band norm
    /// SHARE this allowance. All signed bands count once, without sqrt(2).
    pub envelope_tolerances: Vec<f64>,
}

/// Independent finite-grid output evidence in the observed envelope's units:
/// original physical port units, or those units per angular frequency for a
/// carrier-derivative observation, or original equation units per parameter
/// for a storage-parameter forcing projection.
#[derive(Clone, Debug, PartialEq)]
pub struct PhsPortSidebandCheck {
    pub samples: usize,
    pub omitted_envelope_norm: Vec<f64>,
    /// None on the initial grid; later grids use the same original phase clock.
    pub sampling_change_norm: Option<Vec<f64>>,
}

/// Complete state response and its physical work-conjugate port observation.
#[derive(Clone, Debug)]
pub struct PhsPortSidebandResponse {
    pub omega_rad_s: f64,
    pub carrier_angular_frequency: f64,
    pub state_harmonics: usize,
    pub output_harmonics: usize,
    pub state_dimension: usize,
    pub port_count: usize,
    pub state_response: SidebandResponse,
    /// Sideband-major, then original physical port: `(k+output_harmonics)*m+p`.
    pub port_coefficients: Vec<C64>,
    pub base_residual: f64,
    pub tangent_sampling_defect: f64,
    pub observation_trace: Vec<PhsPortSidebandCheck>,
    pub hessian_actions: usize,
}

/// A refusal publishes no partial physical response and never changes the source.
#[derive(Clone, Debug, PartialEq)]
pub enum PhsPortSidebandError {
    Invalid(&'static str),
    Budget(&'static str),
    SourceRefused,
    Cancelled,
    Sideband(SidebandError),
    ObservationUnresolved {
        last: PhsPortSidebandCheck,
    },
    /// The primal output resolved, but its carrier derivative did not. The
    /// last norms are in physical output units per angular frequency.
    DerivativeObservationUnresolved {
        last: PhsPortSidebandCheck,
    },
    /// Maximum common state order was reached before all physical pair and
    /// independently measured peak/linewidth changes resolved.
    ResonanceUnresolved {
        last: PhsPortResonanceRefinementCheck,
    },
}
impl core::fmt::Display for PhsPortSidebandError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical pHS sideband response: {self:?}")
    }
}
impl std::error::Error for PhsPortSidebandError {}
impl From<SidebandError> for PhsPortSidebandError {
    fn from(value: SidebandError) -> Self {
        if value == SidebandError::Cancelled {
            Self::Cancelled
        } else {
            Self::Sideband(value)
        }
    }
}
type Result<T> = core::result::Result<T, PhsPortSidebandError>;
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PhsPortSidebandError::Cancelled)
    } else {
        Ok(())
    }
}

/// Probe a driven finite-amplitude orbit through all original physical ports.
///
/// `forcing[(k+N)*m+p]` is the PEAK complex perturbation at port p and signed
/// frequency `carrier_angular_frequency + k*orbit.omega`. Missing outside bands
/// are exactly zero by the input definition. The original forcing that sustains
/// the base orbit remains in `input`; it is not replaced by the probe.
///
/// Exact quadratic storage and native complete continuous tangents are discovered
/// automatically. The full source balance and tangent are independently checked
/// by fs-orbit. Physical outputs have a separate shifted-grid refinement because
/// Hxx components hidden by a singular J-R can still be visible at the ports.
/// FFT synthesis/projection costs O((states+ports) S log S); the storage supplies
/// two Hessian actions per sample. No finite-difference tangent is substituted.
///
/// # Errors
/// Invalid inputs or budgets, unavailable/refused source tangents, unresolved
/// base/tangent/response/output bands, singular response, or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_port_sideband_response(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    orbit: &HbOrbit,
    carrier_angular_frequency: f64,
    forcing: &[C64],
    controls: &PhsPortSidebandControls,
    scales: &SidebandScales,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    solve(
        system,
        input,
        orbit,
        carrier_angular_frequency,
        forcing,
        controls,
        scales,
        None,
        cancelled,
    )
}

/// The same complete physical response with an explicitly authored continuous
/// Hessian action of the ORIGINAL storage, including every memory coordinate.
/// The action supplies both the vector-field tangent and physical observation.
/// A declared action that refuses or leaves output unwritten cannot fall back.
///
/// # Errors
/// The native wrapper's refusals, including incomplete/nonfinite Hessian output.
#[allow(clippy::too_many_arguments)]
pub fn solve_forced_port_sideband_response_analytic(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    orbit: &HbOrbit,
    carrier_angular_frequency: f64,
    forcing: &[C64],
    controls: &PhsPortSidebandControls,
    scales: &SidebandScales,
    hessian: &PhsStorageHessianAction<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    solve(
        system,
        input,
        orbit,
        carrier_angular_frequency,
        forcing,
        controls,
        scales,
        Some(hessian),
        cancelled,
    )
}

/// External physical probes about an original self-excited closed-loop orbit.
///
/// The unchanged closed-loop vector field supplies the base balance. Its exact
/// tangent includes `(J-R) Hxx + G du_total/dx`; the feedback Jacobian MUST
/// include the original flow/Hessian chain exactly once. The additional probe
/// enters the same physical G ports without replacing the sustaining source.
/// Every original storage, controller and material-memory coordinate remains.
///
/// The supplied orbit fixes the phase clock for all complex amplitudes. No
/// phase pin, extra damping or frequency detuning regularizes neutral modes;
/// singular shifted carriers remain generic solver refusals. A resolved
/// particular response does not establish transverse Floquet stability.
///
/// # Errors
/// The forced wrapper's numerical/refinement refusals, mismatched feedback
/// ports, an unavailable storage tangent or a refused TOTAL feedback derivative.
#[allow(clippy::too_many_arguments)]
pub fn solve_autonomous_feedback_port_sideband_response(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    orbit: &HbOrbit,
    carrier_angular_frequency: f64,
    forcing: &[C64],
    controls: &PhsPortSidebandControls,
    scales: &SidebandScales,
    cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    solve_autonomous(
        system,
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

/// The same external probe with an authored complete original storage Hessian.
/// It supplies both the storage-flow tangent and physical port observation;
/// the TOTAL feedback derivative must use the same original effort law.
/// # Errors
/// All native autonomous wrapper refusals; declared actions never fall back.
#[allow(clippy::too_many_arguments)]
pub fn solve_autonomous_feedback_port_sideband_response_analytic(
    system: &PortHamiltonian,
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
    solve_autonomous(
        system,
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
fn solve_autonomous(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
    orbit: &HbOrbit,
    carrier: f64,
    forcing: &[C64],
    c: &PhsPortSidebandControls,
    scales: &SidebandScales,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    poll(&mut cancelled)?;
    admit(system, orbit, carrier, forcing, c)?;
    let mut source = PhsPortSidebandSource::from_feedback(
        system,
        None,
        feedback,
        feedback_jacobian,
        orbit,
        c,
        hessian,
    )?;
    let mut prepared = source.prepare(scales, &mut cancelled)?;
    let mut factor = prepared.factor(carrier, &mut cancelled)?;
    factor.solve(forcing, cancelled)
}

#[allow(clippy::too_many_arguments)]
fn solve(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    orbit: &HbOrbit,
    carrier: f64,
    forcing: &[C64],
    c: &PhsPortSidebandControls,
    scales: &SidebandScales,
    hessian: Option<&PhsStorageHessianAction<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhsPortSidebandResponse> {
    poll(&mut cancelled)?;
    admit(system, orbit, carrier, forcing, c)?;
    let mut source = PhsPortSidebandSource::from_forced(system, input, orbit, c, hessian)?;
    let mut prepared = source.prepare(scales, &mut cancelled)?;
    let mut factor = prepared.factor(carrier, &mut cancelled)?;
    factor.solve(forcing, cancelled)
}

fn admit(
    system: &PortHamiltonian,
    orbit: &HbOrbit,
    carrier: f64,
    forcing: &[C64],
    c: &PhsPortSidebandControls,
) -> Result<()> {
    admit_source(system, orbit, c)?;
    admit_carrier(orbit, carrier, c)?;
    admit_forcing(system, forcing, c)
}

pub(super) fn admit_carrier(
    orbit: &HbOrbit,
    carrier: f64,
    c: &PhsPortSidebandControls,
) -> Result<()> {
    if !carrier.is_finite()
        || !((c.output_harmonics as f64).mul_add(orbit.omega, carrier)).is_finite()
        || !(-(c.output_harmonics as f64))
            .mul_add(orbit.omega, carrier)
            .is_finite()
    {
        return Err(PhsPortSidebandError::Invalid(
            "finite signed physical output frequencies required",
        ));
    }
    Ok(())
}

pub(super) fn admit_forcing(
    system: &PortHamiltonian,
    forcing: &[C64],
    c: &PhsPortSidebandControls,
) -> Result<()> {
    let bands = 2 * c.sideband.harmonics + 1;
    if system.port_dim().checked_mul(bands) != Some(forcing.len())
        || forcing.iter().any(|&z| !finite(z))
    {
        return Err(PhsPortSidebandError::Invalid(
            "finite sideband-major original port forcing required",
        ));
    }
    Ok(())
}

pub(super) fn admit_source(
    system: &PortHamiltonian,
    orbit: &HbOrbit,
    c: &PhsPortSidebandControls,
) -> Result<()> {
    let n = system.state_dim();
    let m = system.port_dim();
    let width = orbit.coeffs.first().map_or(0, Vec::len);
    if n == 0
        || m == 0
        || width == 0
        || width > 513
        || c.sideband.harmonics > 64
        || orbit.coeffs.len() != n
        || !orbit.omega.is_finite()
        || orbit.omega <= 0.0
        || orbit
            .coeffs
            .iter()
            .any(|r| r.len() != width || r[0].im != 0.0 || r.iter().any(|&z| !finite(z)))
        || c.initial_samples < 32
        || !c.initial_samples.is_power_of_two()
        || !c.max_samples.is_power_of_two()
        || c.max_samples > 65_536
        || c.initial_samples > c.max_samples / 4
        || c.output_harmonics >= c.initial_samples / 2
        || (width - 1).max(c.sideband.harmonics) >= c.initial_samples / 2
        || c.envelope_tolerances.len() != m
        || c.envelope_tolerances
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(PhsPortSidebandError::Invalid(
            "complete original orbit, resolved grids and positive physical tolerances required",
        ));
    }
    let bands = 2 * c.sideband.harmonics + 1;
    if n as u128 * bands as u128 > c.sideband.maximum_unknowns.min(1024) as u128 {
        return Err(PhsPortSidebandError::Budget(
            "full original state exceeds the sideband dimension cap",
        ));
    }
    if 14 * c.initial_samples > c.max_hessian_actions {
        return Err(PhsPortSidebandError::Budget(
            "two independent output comparisons exceed the Hessian-action cap",
        ));
    }
    let minimum_work =
        input_work(n, m, c.sideband.harmonics) + minimum_observation_work(n, m, c.initial_samples);
    if minimum_work > c.max_observation_work as u128 {
        return Err(PhsPortSidebandError::Budget(
            "two independent output comparisons exceed the arithmetic cap",
        ));
    }
    admit_storage(n, m, c, c.initial_samples)?;
    Ok(())
}

fn admit_storage(n: usize, m: usize, c: &PhsPortSidebandControls, samples: usize) -> Result<()> {
    admit_storage_with_retained(n, m, c, samples, 0)
}

pub(super) fn admit_storage_with_retained(
    n: usize,
    m: usize,
    c: &PhsPortSidebandControls,
    samples: usize,
    retained_bytes: u128,
) -> Result<()> {
    let (n, m, s) = (n as u128, m as u128, samples as u128);
    let bytes = (24 * n + 16 * m + 128) * s
        + 32 * m * (2 * c.output_harmonics + 1) as u128
        + 32 * n * (2 * c.sideband.harmonics + 1) as u128
        + 24 * n * n
        + 64 * n
        + 296 * m
        + retained_bytes;
    if bytes > c.max_workspace_bytes as u128 {
        Err(PhsPortSidebandError::Budget(
            "physical sideband adapter exceeds its extra numeric payload cap",
        ))
    } else {
        Ok(())
    }
}

pub(super) fn input_forcing(
    system: &PortHamiltonian,
    forcing: &[C64],
    harmonics: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>> {
    let (n, m) = (system.state_dim(), system.port_dim());
    let mut rhs = vec![C64::ZERO; n * (2 * harmonics + 1)];
    let mut ports = vec![0.0; m];
    let mut flow = vec![0.0; n];
    for (band, row) in forcing.chunks_exact(m).enumerate() {
        for imaginary in [false, true] {
            poll(cancelled)?;
            for (p, z) in ports.iter_mut().zip(row) {
                *p = if imaginary { z.im } else { z.re };
            }
            system
                .input_flow_into(&ports, &mut flow)
                .map_err(|_| PhsPortSidebandError::SourceRefused)?;
            for (z, &f) in rhs[band * n..(band + 1) * n].iter_mut().zip(&flow) {
                if imaginary {
                    z.im = f;
                } else {
                    z.re = f;
                }
            }
        }
    }
    Ok(rhs)
}

fn observation_work(n: usize, m: usize, samples: usize) -> u128 {
    let (n, m, s) = (n as u128, m as u128, samples as u128);
    s * (64 * (2 * n + m) * samples.ilog2() as u128 + 64 * n * (2 * n + m) + 64 * (n + m))
}

pub(super) fn input_work(n: usize, m: usize, harmonics: usize) -> u128 {
    128 * n as u128 * m as u128 * (2 * harmonics + 1) as u128
}

pub(super) fn minimum_observation_work(n: usize, m: usize, initial_samples: usize) -> u128 {
    [1, 2, 4]
        .iter()
        .map(|factor| observation_work(n, m, initial_samples * factor))
        .sum()
}

// One operation can observe several related envelopes. Retained results are
// charged while the next grid is live; work and actions never restart between
// the primal and its derivative. Single-response calls retain their old costs.
pub(super) struct ObservationBudget {
    pub(super) actions: usize,
    pub(super) work: u128,
    pub(super) retained_bytes: u128,
}

// The same synthesis/projection/refinement owner observes a constitutive
// parameter derivative and the full effort response. No second FFT observer
// or finite-grid error policy is introduced for parameter sensitivities.
#[derive(Clone, Copy)]
pub(super) enum EffortObservation<'a> {
    State,
    Parameter(&'a PhsStorageParameterAction<'a>),
    StateAndParameter(&'a PhsStorageParameterAction<'a>),
}
impl<'a> EffortObservation<'a> {
    fn state(self) -> bool {
        !matches!(self, Self::Parameter(_))
    }
    fn parameter(self) -> Option<&'a PhsStorageParameterAction<'a>> {
        match self {
            Self::State => None,
            Self::Parameter(action) | Self::StateAndParameter(action) => Some(action),
        }
    }
}

fn observe(
    problem: &PhsOrbitProblem<'_>,
    orbit: &HbOrbit,
    response: &[C64],
    c: &PhsPortSidebandControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<PhsPortSidebandCheck>, usize)> {
    let (n, m) = (problem.system.state_dim(), problem.system.port_dim());
    let mut budget = ObservationBudget {
        actions: 0,
        work: input_work(n, m, c.sideband.harmonics),
        retained_bytes: 0,
    };
    observe_with_budget(
        problem,
        orbit,
        response,
        c,
        &c.envelope_tolerances,
        &mut budget,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn observe_with_budget(
    problem: &PhsOrbitProblem<'_>,
    orbit: &HbOrbit,
    response: &[C64],
    c: &PhsPortSidebandControls,
    envelope_tolerances: &[f64],
    budget: &mut ObservationBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<PhsPortSidebandCheck>, usize)> {
    observe_mapped_with_budget(
        problem,
        orbit,
        response,
        c,
        envelope_tolerances,
        problem.system.port_dim(),
        0,
        &|effort, output| {
            problem
                .system
                .output_from_effort_into(effort, output)
                .is_ok()
        },
        budget,
        cancelled,
    )
}

// One original Hessian evaluation can observe several work-conjugate maps.
// Internal callers must admit the full output count, retained payload and
// callback arithmetic before entering this shared refinement owner.
#[allow(clippy::too_many_arguments)]
pub(super) fn observe_mapped_with_budget(
    problem: &PhsOrbitProblem<'_>,
    orbit: &HbOrbit,
    response: &[C64],
    c: &PhsPortSidebandControls,
    envelope_tolerances: &[f64],
    output_count: usize,
    extra_work_per_sample: u128,
    map_effort: &impl Fn(&[f64], &mut [f64]) -> bool,
    budget: &mut ObservationBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<PhsPortSidebandCheck>, usize)> {
    observe_effort_with_budget(
        problem,
        orbit,
        response,
        c,
        envelope_tolerances,
        output_count,
        extra_work_per_sample,
        map_effort,
        EffortObservation::State,
        budget,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn observe_effort_with_budget(
    problem: &PhsOrbitProblem<'_>,
    orbit: &HbOrbit,
    response: &[C64],
    c: &PhsPortSidebandControls,
    envelope_tolerances: &[f64],
    output_count: usize,
    extra_work_per_sample: u128,
    map_effort: &impl Fn(&[f64], &mut [f64]) -> bool,
    effort_observation: EffortObservation<'_>,
    budget: &mut ObservationBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<PhsPortSidebandCheck>, usize)> {
    observe_sampled_effort_with_budget(
        problem,
        orbit,
        response,
        c,
        envelope_tolerances,
        output_count,
        extra_work_per_sample,
        &|_, _, effort, _, out| map_effort(effort, out),
        effort_observation,
        budget,
        cancelled,
    )
}

// The original phase and state are available to a changing physical map.
// A real parameter's direct terms occur only in the real envelope component.
// Legacy maps ignore these arguments and retain their original observer path.
#[allow(clippy::too_many_arguments)]
pub(super) fn observe_sampled_effort_with_budget(
    problem: &PhsOrbitProblem<'_>,
    orbit: &HbOrbit,
    response: &[C64],
    c: &PhsPortSidebandControls,
    envelope_tolerances: &[f64],
    output_count: usize,
    extra_work_per_sample: u128,
    map_effort: &impl Fn(f64, &[f64], &[f64], bool, &mut [f64]) -> bool,
    effort_observation: EffortObservation<'_>,
    budget: &mut ObservationBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<PhsPortSidebandCheck>, usize)> {
    if output_count == 0
        || envelope_tolerances.len() != output_count
        || envelope_tolerances
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(PhsPortSidebandError::Invalid(
            "complete positive mapped-output tolerances required",
        ));
    }
    let (n, m) = (problem.system.state_dim(), output_count);
    let starting_actions = budget.actions;
    let mut previous: Option<Vec<C64>> = None;
    let mut trace = Vec::new();
    let mut samples = c.initial_samples;
    let mut consecutive = 0;
    loop {
        poll(cancelled)?;
        let parameter_scratch = if effort_observation.parameter().is_some() {
            8 * n as u128
        } else {
            0
        };
        admit_storage_with_retained(n, m, c, samples, budget.retained_bytes + parameter_scratch)?;
        budget.actions = budget
            .actions
            .checked_add(if effort_observation.state() {
                2 * samples
            } else {
                0
            })
            .ok_or(PhsPortSidebandError::Budget(
                "physical output Hessian-action count overflow",
            ))?;
        if budget.actions > c.max_hessian_actions {
            return Err(PhsPortSidebandError::Budget(
                "physical output Hessian-action cap exhausted",
            ));
        }
        budget.work = extra_work_per_sample
            .checked_mul(samples as u128)
            .and_then(|extra| observation_work(n, m, samples).checked_add(extra))
            .and_then(|work| budget.work.checked_add(work))
            .ok_or(PhsPortSidebandError::Budget(
                "mapped physical output arithmetic count overflow",
            ))?;
        if budget.work > c.max_observation_work as u128 {
            return Err(PhsPortSidebandError::Budget(
                "physical output arithmetic cap exhausted",
            ));
        }
        let shift = core::f64::consts::TAU * (trace.len() % 3) as f64 / (3.0 * samples as f64);
        let (coefficients, omitted) = project(
            problem,
            orbit,
            response,
            c,
            m,
            map_effort,
            effort_observation,
            samples,
            shift,
            cancelled,
        )?;
        let change = previous.as_ref().map(|old| {
            let mut norms = vec![0.0; m];
            for (i, (&a, &b)) in coefficients.iter().zip(old).enumerate() {
                norms[i % m] = det::hypot(norms[i % m], (a - b).abs());
            }
            norms
        });
        let resolved = change.as_ref().is_some_and(|delta| {
            delta
                .iter()
                .zip(&omitted)
                .zip(envelope_tolerances)
                .all(|((&d, &tail), &limit)| d + tail <= limit)
        });
        consecutive = if resolved { consecutive + 1 } else { 0 };
        let check = PhsPortSidebandCheck {
            samples,
            omitted_envelope_norm: omitted,
            sampling_change_norm: change,
        };
        trace.push(check.clone());
        if consecutive >= 2 {
            return Ok((coefficients, trace, budget.actions - starting_actions));
        }
        if samples == c.max_samples {
            return Err(PhsPortSidebandError::ObservationUnresolved { last: check });
        }
        previous = Some(coefficients);
        samples *= 2;
    }
}

#[allow(clippy::too_many_arguments)]
fn project(
    problem: &PhsOrbitProblem<'_>,
    orbit: &HbOrbit,
    response: &[C64],
    c: &PhsPortSidebandControls,
    output_count: usize,
    map_effort: &impl Fn(f64, &[f64], &[f64], bool, &mut [f64]) -> bool,
    effort_observation: EffortObservation<'_>,
    samples: usize,
    shift: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<f64>)> {
    let (n, m, h) = (
        problem.system.state_dim(),
        output_count,
        c.sideband.harmonics,
    );
    let fft = Fft::new(samples);
    let mut data = vec![FftComplex::default(); samples];
    let mut scratch = data.clone();
    let mut states = vec![0.0; n * samples];
    let mut envelopes = if effort_observation.state() {
        vec![C64::ZERO; n * samples]
    } else {
        Vec::new()
    };
    let mut outputs = vec![C64::ZERO; m * samples];
    for (coordinate, row) in orbit.coeffs.iter().enumerate() {
        poll(cancelled)?;
        data.fill(FftComplex::default());
        data[0].re = row[0].re;
        for (k, &coefficient) in row.iter().enumerate().skip(1) {
            let angle = k as f64 * shift;
            let z = coefficient * C64::new(det::cos(angle), det::sin(angle));
            data[k] = FftComplex::new(z.re, -z.im);
            data[samples - k] = FftComplex::new(z.re, z.im);
        }
        // Reversed Fourier bins synthesize with forward FFT, avoiding an
        // overflow-prone multiply-by-S/inverse-divide-by-S round trip.
        fft.forward(&mut data, &mut scratch);
        for (j, z) in data.iter().enumerate() {
            if !z.re.is_finite() || !z.im.is_finite() {
                return Err(PhsPortSidebandError::SourceRefused);
            }
            states[j * n + coordinate] = z.re;
        }
        if !effort_observation.state() {
            continue;
        }
        poll(cancelled)?;
        data.fill(FftComplex::default());
        for band in 0..2 * h + 1 {
            let k = band as isize - h as isize;
            let angle = k as f64 * shift;
            let z = response[band * n + coordinate] * C64::new(det::cos(angle), det::sin(angle));
            let index = if k > 0 {
                samples - k as usize
            } else {
                k.unsigned_abs()
            };
            data[index] = FftComplex::new(z.re, z.im);
        }
        fft.forward(&mut data, &mut scratch);
        for (j, z) in data.iter().enumerate() {
            if !z.re.is_finite() || !z.im.is_finite() {
                return Err(PhsPortSidebandError::SourceRefused);
            }
            envelopes[j * n + coordinate] = C64::new(z.re, z.im);
        }
    }
    let mut direction = vec![0.0; n];
    let mut effort = vec![0.0; n];
    let mut parameter_effort = if effort_observation.parameter().is_some() {
        vec![0.0; n]
    } else {
        Vec::new()
    };
    let mut flows = vec![0.0; m];
    for j in 0..samples {
        let phase = core::f64::consts::TAU * j as f64 / samples as f64 + shift;
        for imaginary in [false, true] {
            if imaginary && !effort_observation.state() {
                continue;
            }
            poll(cancelled)?;
            if effort_observation.state() {
                for (v, z) in direction.iter_mut().zip(&envelopes[j * n..(j + 1) * n]) {
                    *v = if imaginary { z.im } else { z.re };
                }
                if !problem.tangent.apply_hessian(
                    &states[j * n..(j + 1) * n],
                    &direction,
                    &mut effort,
                ) {
                    return Err(PhsPortSidebandError::SourceRefused);
                }
            } else {
                effort.fill(0.0);
            }
            if !imaginary {
                if let Some(action) = effort_observation.parameter() {
                    parameter_effort.fill(f64::NAN);
                    if !action(&states[j * n..(j + 1) * n], &mut parameter_effort)
                        || parameter_effort.iter().any(|value| !value.is_finite())
                    {
                        return Err(PhsPortSidebandError::SourceRefused);
                    }
                    for (value, direct) in effort.iter_mut().zip(&parameter_effort) {
                        *value += direct;
                    }
                }
            }
            flows.fill(f64::NAN);
            if effort.iter().any(|value| !value.is_finite())
                || !map_effort(
                    phase,
                    &states[j * n..(j + 1) * n],
                    &effort,
                    imaginary,
                    &mut flows,
                )
                || flows.iter().any(|x| !x.is_finite())
            {
                return Err(PhsPortSidebandError::SourceRefused);
            }
            for (p, &flow) in flows.iter().enumerate() {
                if imaginary {
                    outputs[p * samples + j].im = flow;
                } else {
                    outputs[p * samples + j].re = flow;
                }
            }
        }
    }
    let mut coefficients = vec![C64::ZERO; m * (2 * c.output_harmonics + 1)];
    let mut omitted = vec![0.0; m];
    for (p, values) in outputs.chunks_exact(samples).enumerate() {
        poll(cancelled)?;
        let scale = values
            .iter()
            .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
        let scale = if scale == 0.0 { 1.0 } else { scale };
        for (z, &value) in data.iter_mut().zip(values) {
            *z = FftComplex::new(
                (value.re / scale) / samples as f64,
                (value.im / scale) / samples as f64,
            );
            if (value.re != 0.0 && z.re == 0.0) || (value.im != 0.0 && z.im == 0.0) {
                return Err(PhsPortSidebandError::Invalid(
                    "physical output Fourier normalization is not representable",
                ));
            }
        }
        fft.forward(&mut data, &mut scratch);
        for (bin, &value) in data.iter().enumerate() {
            let k = if bin <= samples / 2 {
                bin as isize
            } else {
                bin as isize - samples as isize
            };
            let angle = -(k as f64) * shift;
            let z = C64::new(value.re, value.im).scale(scale)
                * C64::new(det::cos(angle), det::sin(angle));
            if !finite(z) {
                return Err(PhsPortSidebandError::SourceRefused);
            }
            if k.unsigned_abs() <= c.output_harmonics {
                coefficients[(k + c.output_harmonics as isize) as usize * m + p] = z;
            } else {
                // The envelope is complex. Each measured signed bin, including
                // Nyquist, counts once; there is no inferred conjugate partner.
                omitted[p] = det::hypot(omitted[p], z.abs());
            }
        }
        if !omitted[p].is_finite() {
            return Err(PhsPortSidebandError::SourceRefused);
        }
    }
    Ok((coefficients, omitted))
}
