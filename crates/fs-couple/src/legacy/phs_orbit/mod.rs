//! Finite-amplitude periodic regimes of complete nonlinear pHS ODEs.
//!
//! This is an adapter, not a second orbit solver. fs-orbit owns harmonic balance,
//! shooting and continuation; fs-phs owns energy, interconnection, dissipation
//! and material-memory coordinates. The adapter evaluates the ORIGINAL continuous
//! vector field at every AFT sample. The basic ODE adapter exposes the default sI
//! frequency port; [`impedance`] adds physical frequency-dependent back-reaction.
//! Descriptor systems are intentionally excluded because their algebraic rows need
//! an explicit descriptor OrbitProblem rather than an inverted mass matrix.
//!
//! The path is useful after small-signal pole analysis: linear poles seed the
//! first periodic branch, while finite amplitude keeps nonlinear storage and every
//! explicit Maxwell coordinate. Harmonic truncation remains disclosed structure;
//! time-domain confirmation is still required for attacks and aperiodic regimes.

use core::cell::{Cell, RefCell};
use fs_math::{c64::C64, det};
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbScales, OrbitError,
    OrbitProblem, ScaledHbOrbit, solve_hb, solve_hb_adaptive, solve_hb_adaptive_scaled_seeded,
    solve_hb_scaled_seeded,
};
use fs_phs::PortHamiltonian;

mod tangent;
pub use tangent::{
    PhsStorageHessianAction, solve_forced_harmonic_balance_analytic_adaptive_scaled_seeded,
    solve_forced_harmonic_balance_analytic_scaled_seeded,
};

mod sensitivity;
pub use sensitivity::{
    PhsForcedParameterDirection, PhsForcedSensitivityError,
    PreparedPhsForcedAdjoint, PreparedPhsForcedSensitivity,
    PhsPowerObjective, PreparedPhsForcedPowerAdjoint,
};

mod backbone_sensitivity;
pub use backbone_sensitivity::{
    PhsBackboneParameterDirection, PhsBackboneSensitivityError, PreparedPhsBackboneSensitivity,
};

/// Periodically modulated original Hamiltonians and exact physical pump work.
pub mod parametric;

/// Finite-amplitude internal nonlinear dissipation in the original Hamiltonian.
pub mod dissipative;

/// Physical port harmonics evaluated through the original nonlinear storage.
pub mod spectrum;

/// Local finite-amplitude resonance surveys of the original forced response.
pub mod forced_resonance;

/// Original harmonic and sideband outputs through physical BEM receivers.
pub mod acoustic;

/// Full physical impedance back-reaction on forced nonlinear periodic regimes.
pub mod impedance;

/// Original causal rational boundary states in nonlinear HB and variational shooting.
pub mod rational;

/// Physical probe response and frequency conversion about driven or self-excited periodic orbits.
pub mod sideband;

mod continuation;
pub use continuation::{
    PhsForcedBranchStabilityControls, PhsForcedPortBranchControls, PhsForcedPortBranchError,
    PhsForcedPortBranchPoint, continue_forced_port_response_analytic_scaled_seeded,
    continue_forced_port_response_scaled_seeded,
    PhsFrequencyResponsePoint, PhsScaledFrequencyResponsePoint,
    continue_forced_frequency_response, continue_forced_frequency_response_scaled_seeded,
};
mod stability;
pub use stability::{
    PhsFloquetControls, PhsFloquetLevel, PhsForcedStability, PhsNumericalStability,
    PhsVariationalFloquetLevel, PhsVariationalForcedStability,
    analyze_forced_orbit_stability, analyze_forced_orbit_stability_variational,
    analyze_forced_orbit_stability_variational_analytic,
};
mod shooting;
pub use shooting::{
    solve_forced_shooting_variational, solve_forced_shooting_variational_analytic,
};
mod variational_feedback;
pub use variational_feedback::{
    PhsFeedbackParameterAction, PhsFeedbackStateJacobian,
    solve_autonomous_feedback_harmonic_balance_adaptive_scaled_seeded,
    solve_autonomous_feedback_harmonic_balance_analytic_adaptive_scaled_seeded,
    solve_autonomous_feedback_harmonic_balance_analytic_scaled_seeded,
    solve_autonomous_feedback_harmonic_balance_scaled_seeded,
    solve_autonomous_feedback_shooting_variational,
    solve_autonomous_feedback_shooting_variational_analytic,
};
mod forced_feedback;
pub use forced_feedback::{
    PhsForcedFeedbackAudit, audit_forced_feedback_orbit,
    solve_forced_feedback_harmonic_balance_adaptive_scaled_seeded,
    solve_forced_feedback_harmonic_balance_analytic_adaptive_scaled_seeded,
    solve_forced_feedback_harmonic_balance_analytic_scaled_seeded,
    solve_forced_feedback_harmonic_balance_scaled_seeded,
    solve_forced_feedback_shooting_variational,
    solve_forced_feedback_shooting_variational_analytic,
};
mod autonomous;
/// Full-state parameter branches of the original self-excited pHS.
pub use autonomous::continuation as autonomous_continuation;
pub use autonomous::{
    AutonomousPortFeedback, DissipativeAutonomousHamiltonian, DissipativeAutonomousParameterDirection,
    DissipativeAutonomousPowerObjective, DissipativeAutonomousSensitivityError,
    PreparedDissipativeAutonomousAdjoint, PreparedDissipativeAutonomousPowerAdjoint,
    PreparedDissipativeAutonomousSensitivity, PhsAutonomousControls, PhsAutonomousLevel,
    PhsAutonomousParameterDirection, PhsAutonomousSensitivityError,
    PhsAutonomousStability, PreparedPhsAutonomousAdjoint, PreparedPhsAutonomousSensitivity,
    PreparedPhsAutonomousPowerAdjoint,
    analyze_autonomous_feedback_stability,
    audit_autonomous_feedback_orbit, solve_autonomous_feedback_harmonic_balance,
    solve_autonomous_feedback_harmonic_balance_adaptive,
};

const MAX_INPUT_HARMONICS: usize = 256;

/// Why a periodic regime was not produced.
#[derive(Debug, Clone, PartialEq)]
pub enum PhsOrbitError {
    /// An argument or input row was malformed (non-finite, non-positive or out of range).
    Invalid {
        /// The violated requirement.
        what: &'static str,
    },
    /// The original storage/input evaluation refused at an AFT or time-domain sample.
    SourceRefused,
    /// Shooting moved farther from the supplied phase-zero state than permitted.
    SeedDisagreement {
        /// Original pHS state coordinate.
        coordinate: usize,
        /// Absolute difference, in this coordinate's physical units.
        absolute_change: f64,
        /// Caller-authored limit in the same units.
        tolerance: f64,
    },
    /// fs-orbit's underlying periodic solve refused.
    Orbit(OrbitError),
    /// Exact discrete variational shooting, including cancellation and budgets.
    Variational(fs_orbit::VariationalShootError),
}
impl core::fmt::Display for PhsOrbitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pHS periodic regime refused: {self:?}")
    }
}
impl std::error::Error for PhsOrbitError {}
impl From<OrbitError> for PhsOrbitError {
    fn from(value: OrbitError) -> Self { Self::Orbit(value) }
}
fn bad(what: &'static str) -> PhsOrbitError { PhsOrbitError::Invalid { what } }

/// Peak-periodic physical port input in fs-orbit's exp(+i*n*omega*t) chart.
/// coefficients[port][0] is the real DC value; n>=1 stores the positive-frequency
/// complex Fourier coefficient, so a real cosine of peak amplitude A uses A/2.
#[derive(Debug, Clone)]
pub struct PeriodicPortInput {
    omega_rad_s: f64,
    coefficients: Vec<Vec<C64>>,
}
impl PeriodicPortInput {
    /// Validate a port input: a positive finite fundamental, at most 256 ports, and
    /// equal-width finite Fourier rows (at most 256 harmonics) with real DC terms.
    ///
    /// # Errors
    /// [`PhsOrbitError::Invalid`] naming the violated requirement.
    pub fn new(omega_rad_s: f64, coefficients: Vec<Vec<C64>>) -> Result<Self, PhsOrbitError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 || coefficients.len() > 256 {
            return Err(bad("positive finite frequency and at most 256 ports required"));
        }
        let width = coefficients.first().map_or(1, Vec::len);
        if width == 0 || width > MAX_INPUT_HARMONICS + 1
            || coefficients.iter().any(|row| row.len() != width
                || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
                || row[0].im != 0.0)
        {
            return Err(bad("complete finite Fourier rows with real DC terms required"));
        }
        Ok(Self { omega_rad_s, coefficients })
    }
    /// Fundamental angular frequency [rad/s].
    pub fn omega_rad_s(&self) -> f64 { self.omega_rad_s }
    /// Reuse the same DC values, complex harmonic amplitudes and phase origin
    /// at another positive finite fundamental, for example a continuation point.
    pub fn at_frequency(&self, omega_rad_s: f64) -> Result<Self, PhsOrbitError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(bad("positive finite forcing frequency required"));
        }
        Ok(Self { omega_rad_s, coefficients: self.coefficients.clone() })
    }
    /// Number of driven ports (Fourier rows).
    pub fn port_count(&self) -> usize { self.coefficients.len() }
    /// Declared harmonics per port, excluding the DC term, including zero padding.
    pub fn harmonic_count(&self) -> usize {
        self.coefficients.first().map_or(0, |row| row.len().saturating_sub(1))
    }
    /// Highest harmonic with a nonzero coefficient on any physical port.
    /// Exact-zero padding does not increase this bandwidth; DC-only and zero
    /// inputs both return zero. No amplitude threshold is applied.
    pub fn highest_nonzero_harmonic(&self) -> usize {
        self.coefficients.iter().map(|row|
            row.iter().rposition(|c| c.re != 0.0 || c.im != 0.0).unwrap_or(0)
        ).max().unwrap_or(0)
    }
    fn sample_into(&self, time_s: f64, out: &mut [f64]) -> bool {
        self.sample_at_frequency_into(time_s, self.omega_rad_s, out)
    }
    fn sample_at_frequency_into(&self, time_s: f64, omega_rad_s: f64, out: &mut [f64]) -> bool {
        if out.len() != self.coefficients.len() || !time_s.is_finite()
            || !omega_rad_s.is_finite() || omega_rad_s <= 0.0
        { return false; }
        for (port, row) in self.coefficients.iter().enumerate() {
            let mut value = row[0].re;
            let mut correction = 0.0;
            for (n, &coefficient) in row.iter().enumerate().skip(1) {
                if coefficient.re == 0.0 && coefficient.im == 0.0 { continue; }
                let phase = n as f64 * omega_rad_s * time_s;
                let term = 2.0 * (coefficient.re * det::cos(phase)
                    - coefficient.im * det::sin(phase));
                let adjusted = term - correction;
                let next = value + adjusted;
                correction = (next - value) - adjusted;
                value = next;
            }
            if !value.is_finite() { return false; }
            out[port] = value;
        }
        true
    }
}

/// Private adapter so source-domain refusals can be promoted by the public solve
/// wrappers instead of escaping as an undifferentiated Newton stall.
struct PhsOrbitProblem<'a> {
    system: &'a PortHamiltonian,
    input: InputKind<'a>,
    tangent: tangent::StorageTangent<'a>,
    effort: RefCell<Vec<f64>>,
    ports: RefCell<Vec<f64>>,
    source_refused: Cell<bool>,
}
enum InputKind<'a> {
    Forced { periodic: &'a PeriodicPortInput, omega_rad_s: f64 },
    Constant(&'a [f64]),
    Feedback { law: &'a dyn AutonomousPortFeedback, flows: RefCell<Vec<f64>> },
}
impl PhsOrbitProblem<'_> {
    fn source_refused(&self) -> bool { self.source_refused.get() }
    fn fail(&self, out: &mut [f64]) {
        self.source_refused.set(true);
        out.fill(f64::NAN);
    }
}
impl OrbitProblem for PhsOrbitProblem<'_> {
    fn dim(&self) -> usize { self.system.state_dim() }
    fn autonomous(&self) -> bool { !matches!(self.input, InputKind::Forced { .. }) }
    fn island_state_jacobian(&self) -> Option<fs_orbit::IslandStateJacobian<Self>> {
        self.tangent.available().then_some(Self::state_jacobian)
    }
    fn time_domain_rhs(&self) -> Option<fs_orbit::TimeDomainRhs<Self>> {
        // This adapter's frequency port is exactly sI: island already evaluates
        // the complete original ODE, including all explicit memory coordinates.
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<fs_orbit::TimeDomainStateJacobian<Self>> {
        // Only this complete ODE adapter has the exact sI port. Its declared
        // storage tangent differentiates the whole flow for state-independent
        // forcing/bias. An arbitrary feedback law also needs its own derivative.
        (!matches!(self.input, InputKind::Feedback { .. }) && self.tangent.available())
            .then_some(Self::state_jacobian)
    }

    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        self.island_with_external_input(t, x, out, None);
    }
}

impl PhsOrbitProblem<'_> {
    fn island_with_external_input(
        &self, t: f64, x: &[f64], out: &mut [f64], external: Option<&PeriodicPortInput>,
    ) {
        let n = self.system.state_dim();
        let m = self.system.port_dim();
        if x.len() != n || out.len() != n || x.iter().any(|v| !v.is_finite()) {
            self.fail(out); return;
        }
        let mut effort = self.effort.borrow_mut();
        effort.fill(f64::NAN);
        self.system.storage().gradient(x, &mut effort);
        if effort.iter().any(|v| !v.is_finite()) { self.fail(out); return; }

        let mut input = self.ports.borrow_mut();
        let ok = match &self.input {
            InputKind::Forced { periodic, omega_rad_s } =>
                periodic.sample_at_frequency_into(t, *omega_rad_s, &mut input),
            InputKind::Constant(bias) => {
                if bias.len() != m || bias.iter().any(|v| !v.is_finite()) { false }
                else { input.copy_from_slice(bias); true }
            }
            InputKind::Feedback { law, flows } => {
                let mut flows = flows.borrow_mut();
                input.fill(f64::NAN);
                port_flows_into(self.system, &effort, &mut flows)
                    && law.input(x, &flows, &mut input)
                    && input.iter().all(|value| value.is_finite())
            }
        };
        if !ok { self.fail(out); return; }

        if let Some(periodic) = external {
            // Feedback has consumed the conjugate flows. Reuse that scratch
            // for the state-independent drive, then apply the complete force
            // through the original compensated J/R/G owner in one contraction.
            let InputKind::Feedback { flows, .. } = &self.input else {
                self.fail(out); return;
            };
            let mut drive = flows.borrow_mut();
            if !periodic.sample_into(t, &mut drive) { self.fail(out); return; }
            for (value, imposed) in input.iter_mut().zip(drive.iter()) {
                *value += imposed;
            }
            if input.iter().any(|v| !v.is_finite()) { self.fail(out); return; }
        }

        if self.system.vector_field_from_effort_into(&effort, &input, out).is_err() {
            self.fail(out);
        }
    }
}

fn forced_problem<'a>(
    system: &'a PortHamiltonian, input: &'a PeriodicPortInput,
) -> Result<PhsOrbitProblem<'a>, PhsOrbitError> {
    if system.state_dim() == 0 || system.port_dim() != input.port_count() {
        return Err(bad("periodic input must cover every original physical port"));
    }
    Ok(PhsOrbitProblem {
        system, input: InputKind::Forced { periodic: input, omega_rad_s: input.omega_rad_s() },
        tangent: tangent::StorageTangent::from_storage(system),
        effort: RefCell::new(vec![0.0; system.state_dim()]),
        ports: RefCell::new(vec![0.0; system.port_dim()]),
        source_refused: Cell::new(false),
    })
}
fn autonomous_problem<'a>(
    system: &'a PortHamiltonian, bias: &'a [f64],
) -> Result<PhsOrbitProblem<'a>, PhsOrbitError> {
    if system.state_dim() == 0 || bias.len() != system.port_dim()
        || bias.iter().any(|v| !v.is_finite())
    {
        return Err(bad("constant bias must cover every original physical port"));
    }
    Ok(PhsOrbitProblem {
        system, input: InputKind::Constant(bias),
        tangent: tangent::StorageTangent::from_storage(system),
        effort: RefCell::new(vec![0.0; system.state_dim()]),
        ports: RefCell::new(bias.to_vec()),
        source_refused: Cell::new(false),
    })
}

/// Forced finite-amplitude harmonic balance using the original nonlinear pHS.
/// The forcing frequency is taken from input, preventing an inconsistent anchor.
/// Every nonzero input harmonic must fit the retained HB truncation. Admission
/// precedes storage/source evaluation, so unresolved forcing cannot alias into a
/// lower harmonic. Trailing zero coefficients do not consume this bandwidth.
/// Nonlinear response harmonics still require truncation refinement; use
/// [`solve_forced_harmonic_balance_adaptive`] to request automatic refinement.
pub fn solve_forced_harmonic_balance(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    guess_amplitude: f64,
    budget: &HbBudget,
) -> Result<HbOrbit, PhsOrbitError> {
    if !guess_amplitude.is_finite() || guess_amplitude <= 0.0 {
        return Err(bad("positive finite first-state seed amplitude required"));
    }
    if input.highest_nonzero_harmonic() > budget.harmonics {
        return Err(bad("nonzero forcing harmonic exceeds the retained HB truncation"));
    }
    let problem = forced_problem(system, input)?;
    let result = solve_hb(&problem, HbAnchor::Forced { omega: input.omega_rad_s() },
        guess_amplitude, budget);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

/// Forced physical pHS response with explicit mixed-coordinate Newton scales.
///
/// The packed seed uses every original state, including material memory, in
/// fs-orbit's canonical physical Fourier ordering. State and equation scales
/// respectively carry the original coordinate and time-derivative units. They
/// change numerical conditioning and accuracy admission, never the Hamiltonian,
/// physical port input or dissipation. The initial order must retain every
/// nonzero authored force harmonic; unlike an amplitude-only seed it is not
/// silently resized. See [`solve_hb_scaled_seeded`] for convergence semantics.
pub fn solve_forced_harmonic_balance_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    if input.highest_nonzero_harmonic() > budget.harmonics {
        return Err(bad(
            "nonzero forcing harmonic exceeds the retained HB truncation",
        ));
    }
    let problem = forced_problem(system, input)?;
    let result = solve_hb_scaled_seeded(
        &problem,
        HbAnchor::Forced {
            omega: input.omega_rad_s(),
        },
        seed,
        budget,
        scales,
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// The same physically scaled pHS response with automatic harmonic refinement.
///
/// Retains the original storage, complete state, signed multiport input and
/// source-refusal behavior. The seed's initial order must include the input
/// bandwidth; per-state and per-equation refinement tolerances stay in physical
/// units at every larger solve. Finite-band agreement is not an infinite-tail
/// certificate or a time-domain stability claim.
pub fn solve_forced_harmonic_balance_adaptive_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    if input.highest_nonzero_harmonic() > budget.harmonics {
        return Err(bad(
            "nonzero forcing harmonic exceeds the retained HB truncation",
        ));
    }
    let problem = forced_problem(system, input)?;
    let result = solve_hb_adaptive_scaled_seeded(
        &problem,
        HbAnchor::Forced {
            omega: input.omega_rad_s(),
        },
        seed,
        budget,
        refinement,
        scales,
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// Autonomous finite-amplitude orbit about a constant physical port bias.
pub fn solve_autonomous_harmonic_balance(
    system: &PortHamiltonian,
    bias: &[f64],
    omega_guess: f64,
    guess_amplitude: f64,
    budget: &HbBudget,
) -> Result<HbOrbit, PhsOrbitError> {
    if !omega_guess.is_finite() || omega_guess <= 0.0
        || !guess_amplitude.is_finite() || guess_amplitude <= 0.0
    {
        return Err(bad("positive finite frequency and seed amplitude required"));
    }
    let problem = autonomous_problem(system, bias)?;
    let result = solve_hb(&problem, HbAnchor::Autonomous { omega_guess },
        guess_amplitude, budget);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

/// Conservative/backbone point with first-state amplitude pinned by fs-orbit.
/// The source system is not required to be lossless; if dissipation prevents a
/// periodic backbone the original balance refuses rather than deleting the loss.
pub fn solve_backbone_harmonic_balance(
    system: &PortHamiltonian, bias: &[f64], amplitude: f64, omega_guess: f64,
    budget: &HbBudget,
) -> Result<HbOrbit, PhsOrbitError> {
    if !amplitude.is_finite() || amplitude <= 0.0
        || !omega_guess.is_finite() || omega_guess <= 0.0
    {
        return Err(bad("positive finite backbone amplitude and frequency required"));
    }
    let problem = autonomous_problem(system, bias)?;
    let result = solve_hb(&problem, HbAnchor::Backbone { amplitude, omega_guess },
        amplitude, budget);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}


/// Forced pHS response with automatic retained-harmonic refinement.
///
/// The first solve retains at least the highest nonzero authored input harmonic.
/// The maximum order must leave room for a second, larger solve; unresolved
/// forcing is refused before storage evaluation. Zero padding does not change
/// either the starting order or the result. State tolerances use the original
/// pHS coordinates (for example displacement and momentum, including memory);
/// residual tolerances use their time-derivative units for this ODE adapter.
///
/// The returned evidence checks finite spectral bands and successive solutions,
/// not an infinite-spectrum bound or a proof of periodic-orbit existence.
pub fn solve_forced_harmonic_balance_adaptive(
    system: &PortHamiltonian, input: &PeriodicPortInput,
    guess_amplitude: f64, budget: &HbBudget, refinement: &HbRefinementBudget,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    if !guess_amplitude.is_finite() || guess_amplitude <= 0.0 {
        return Err(bad("positive finite first-state seed amplitude required"));
    }
    if budget.harmonics == 0 {
        return Err(bad("initial harmonic count must be positive"));
    }
    let mut initial = *budget;
    initial.harmonics = initial.harmonics.max(input.highest_nonzero_harmonic());
    if initial.harmonics >= refinement.max_harmonics {
        return Err(bad("adaptive maximum must exceed the starting and forcing harmonic orders"));
    }
    let problem = forced_problem(system, input)?;
    let result = solve_hb_adaptive(&problem,
        HbAnchor::Forced { omega: input.omega_rad_s() },
        guess_amplitude, &initial, refinement);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

/// Autonomous pHS orbit with retained-harmonic, waveform, and frequency checks.
///
/// Uses the original storage, dissipation, physical bias, and memory coordinates.
/// Per-coordinate tolerances have the units described by
/// [`solve_forced_harmonic_balance_adaptive`]. A non-isolated conservative family
/// should instead be given an explicit amplitude through the backbone entry.
pub fn solve_autonomous_harmonic_balance_adaptive(
    system: &PortHamiltonian, bias: &[f64], omega_guess: f64,
    guess_amplitude: f64, budget: &HbBudget, refinement: &HbRefinementBudget,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    if !omega_guess.is_finite() || omega_guess <= 0.0
        || !guess_amplitude.is_finite() || guess_amplitude <= 0.0
    {
        return Err(bad("positive finite frequency and seed amplitude required"));
    }
    let problem = autonomous_problem(system, bias)?;
    let result = solve_hb_adaptive(&problem,
        HbAnchor::Autonomous { omega_guess }, guess_amplitude, budget, refinement);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

/// Amplitude-constrained pHS backbone with automatic harmonic refinement.
///
/// The pinned quantity remains the first-state fundamental amplitude. Nonlinear
/// storage can generate higher harmonics and change the physical turning point.
/// Every original balance row remains enforced; no loss or memory law is removed.
pub fn solve_backbone_harmonic_balance_adaptive(
    system: &PortHamiltonian, bias: &[f64], amplitude: f64, omega_guess: f64,
    budget: &HbBudget, refinement: &HbRefinementBudget,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    if !amplitude.is_finite() || amplitude <= 0.0
        || !omega_guess.is_finite() || omega_guess <= 0.0
    {
        return Err(bad("positive finite backbone amplitude and frequency required"));
    }
    let problem = autonomous_problem(system, bias)?;
    let result = solve_hb_adaptive(&problem,
        HbAnchor::Backbone { amplitude, omega_guess }, amplitude, budget, refinement);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

/// Physical whole-cycle observation of an HB orbit against the ORIGINAL pHS.
/// This is an independent quadrature audit, not part of fs-orbit's Newton
/// convergence test. Nonlinear storage may generate harmonics above the retained
/// HB truncation, so callers should repeat with a sample ladder before treating
/// a small residual as converged evidence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhsOrbitAudit {
    /// Quadrature samples over one period.
    pub samples: usize,
    /// Cycle mean of the original Hamiltonian [J].
    pub mean_stored_energy_j: f64,
    /// Smallest stored energy at the quadrature samples [J].
    pub minimum_stored_energy_j: f64,
    /// Largest stored energy at the quadrature samples [J].
    pub maximum_stored_energy_j: f64,
    /// Cycle mean of the power the port input supplies [W].
    pub mean_supplied_power_w: f64,
    /// Cycle mean of the native resistive dissipation [W].
    pub mean_resistance_loss_w: f64,
    /// supply - native R loss; a periodic exact solution should approach zero.
    pub cycle_power_residual_w: f64,
    /// |residual| / mean(|instantaneous supply| + native loss), or zero at rest.
    pub relative_power_residual: f64,
}

/// Audit one FORCED orbit using its exact authored Fourier input and original
/// pHS energy/resistance. No damping or energy correction is applied. The grid
/// must resolve every authored drive and retained state harmonic; nonlinear
/// work can need finer sampling and still requires successive-grid comparison.
pub fn audit_forced_orbit(
    system:&PortHamiltonian, input:&PeriodicPortInput, orbit:&HbOrbit, samples:usize,
) -> Result<PhsOrbitAudit,PhsOrbitError> {
    if system.port_dim()!=input.port_count() || orbit.coeffs.len()!=system.state_dim()
        || orbit.omega.to_bits()!=input.omega_rad_s().to_bits()
        || samples<32 || samples>65_536
        || input.highest_nonzero_harmonic() > samples.saturating_sub(1) / 2
        || orbit.coeffs.iter().any(|row|row.is_empty())
    { return Err(bad("cycle audit needs matching forced orbit and 32..=65536 samples")); }
    audit_orbit_with_input(system, orbit, samples,
        |time, _, _, out| input.sample_into(time, out))
}

fn port_flows_into(system: &PortHamiltonian, effort: &[f64], out: &mut [f64]) -> bool {
    system.output_from_effort_into(effort, out).is_ok()
}

fn audit_orbit_with_input(
    system: &PortHamiltonian, orbit: &HbOrbit, samples: usize,
    mut input: impl FnMut(f64, &[f64], &[f64], &mut [f64]) -> bool,
) -> Result<PhsOrbitAudit, PhsOrbitError> {
    audit_orbit_with_input_effort(system, orbit, samples,
        |time, state, _, flows, out| input(time, state, flows, out))
}

fn audit_orbit_with_input_effort(
    system: &PortHamiltonian, orbit: &HbOrbit, samples: usize,
    mut input: impl FnMut(f64, &[f64], &[f64], &[f64], &mut [f64]) -> bool,
) -> Result<PhsOrbitAudit, PhsOrbitError> {
    let width=orbit.coeffs.first().map_or(0,Vec::len);
    if system.state_dim()==0 || orbit.coeffs.len()!=system.state_dim() || width==0
        || !orbit.omega.is_finite() || orbit.omega<=0.0
        || !(core::f64::consts::TAU/orbit.omega).is_finite()
        || !(32..=65_536).contains(&samples)
        || width.saturating_sub(1) > samples.saturating_sub(1) / 2
        || orbit.coeffs.iter().any(|row|row.len()!=width || row[0].im!=0.0
            || row.iter().any(|z|!z.re.is_finite() || !z.im.is_finite()))
    { return Err(bad("cycle audit needs finite original-state Fourier rows and 32..=65536 samples")); }
    let n=system.state_dim(); let m=system.port_dim();
    let mut x=vec![0.0;n]; let mut e=vec![0.0;n]; let mut u=vec![0.0;m];
    let mut flows=vec![0.0;m];
    let (mut energy_sum,mut supply_sum,mut loss_sum,mut scale_sum)=(0.0,0.0,0.0,0.0);
    let (mut ce,mut cs,mut cl,mut cx)=(0.0,0.0,0.0,0.0);
    let (mut emin,mut emax)=(f64::INFINITY,f64::NEG_INFINITY);
    let add=|sum:&mut f64,c:&mut f64,value:f64| {let y=value-*c;let next=*sum+y;*c=(next-*sum)-y;*sum=next;};
    for j in 0..samples {
        let theta=core::f64::consts::TAU*j as f64/samples as f64;
        for (k,value) in x.iter_mut().enumerate() {*value=orbit.sample(k,theta);}
        e.fill(f64::NAN); system.storage().gradient(&x,&mut e);
        if e.iter().any(|v|!v.is_finite()) {return Err(PhsOrbitError::SourceRefused);}
        let time=theta/orbit.omega;
        u.fill(f64::NAN);
        if !port_flows_into(system,&e,&mut flows) || !input(time,&x,&e,&flows,&mut u)
            || u.iter().any(|value|!value.is_finite())
        {return Err(PhsOrbitError::SourceRefused);}
        let h=system.hamiltonian(&x);
        // Use the original effort and map: individually rounded output ports
        // can hide real work when large opposite source powers cancel.
        let supply=system.supplied_power_from_effort(&e,&u)
            .map_err(|_|PhsOrbitError::SourceRefused)?;
        let loss=system.resistance_power_from_effort(&e)
            .map_err(|_|PhsOrbitError::SourceRefused)?;
        if ![h,supply,loss].iter().all(|v|v.is_finite()) {
            return Err(bad("cycle audit produced nonfinite physical work or storage"));
        }
        add(&mut energy_sum,&mut ce,h); add(&mut supply_sum,&mut cs,supply);
        add(&mut loss_sum,&mut cl,loss); add(&mut scale_sum,&mut cx,supply.abs()+loss.abs());
        emin=emin.min(h); emax=emax.max(h);
    }
    let inv=1.0/samples as f64;
    let mean_h=energy_sum*inv; let mean_supply=supply_sum*inv; let mean_loss=loss_sum*inv;
    let residual=mean_supply-mean_loss; let scale=scale_sum*inv;
    let relative=if scale==0.0 {0.0}else{residual.abs()/scale};
    if ![mean_h,mean_supply,mean_loss,residual,relative,emin,emax].iter().all(|v|v.is_finite()) {
        return Err(bad("cycle audit reduction overflow"));
    }
    Ok(PhsOrbitAudit {samples,mean_stored_energy_j:mean_h,minimum_stored_energy_j:emin,
        maximum_stored_energy_j:emax,mean_supplied_power_w:mean_supply,
        mean_resistance_loss_w:mean_loss,cycle_power_residual_w:residual,
        relative_power_residual:relative})
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage};

    fn hb(n: usize) -> HbBudget {
        HbBudget { harmonics:n, max_newton:60, tolerance:1e-11, deflate_equilibrium:false }
    }

    fn refinement(dimension: usize, max_harmonics: usize) -> HbRefinementBudget {
        HbRefinementBudget {
            max_harmonics,
            residual_tolerances: vec![1e-8; dimension],
            state_tolerances: vec![1e-7; dimension],
            frequency_tolerance: 1e-9,
        }
    }

    #[test]
    fn g1_coherent_port_forcing_retains_the_weak_original_drive() {
        let system=PortHamiltonian::new(1,3,vec![0.0],vec![0.0],vec![1.0;3],
            Box::new(QuadraticStorage::new(vec![1.0],1).unwrap())).unwrap();
        let large=2.0_f64.powi(54);
        let input=[large,1.0,-large];
        let problem=autonomous_problem(&system,&input).unwrap();
        let mut derivative=[f64::NAN];
        problem.island(0.0,&[0.0],&mut derivative);
        assert!(!problem.source_refused());
        assert_eq!(derivative,[1.0]);
    }

    #[test]
    fn g1_forcing_and_internal_dynamics_cancel_before_final_rounding() {
        let system=PortHamiltonian::new(1,2,vec![0.0],vec![1.0],vec![1.0;2],
            Box::new(QuadraticStorage::new(vec![1.0],1).unwrap())).unwrap();
        let large=2.0_f64.powi(54);
        let input=[large,1.0];
        let problem=autonomous_problem(&system,&input).unwrap();
        let mut derivative=[f64::NAN];
        // -R*effort + G*u = -L + L + 1. Rounding Gu first would lose 1.
        problem.island(0.0,&[large],&mut derivative);
        assert!(!problem.source_refused());
        assert_eq!(derivative,[1.0]);
    }

    #[test]
    fn g1_cycle_work_uses_original_map_even_when_observed_ports_round_equal() {
        let system=PortHamiltonian::new(2,2,vec![0.0;4],
            vec![0.0,0.0,0.0,1.0],vec![1.0,1.0,1.0,0.0],
            Box::new(QuadraticStorage::new(vec![1.0,0.0,0.0,1.0],2).unwrap())).unwrap();
        let large=2.0_f64.powi(54);
        let state=[large,1.0];
        let input=PeriodicPortInput::new(1.0,
            vec![vec![C64::ONE],vec![C64::from_re(-1.0)]]).unwrap();
        let problem=forced_problem(&system,&input).unwrap();
        let mut derivative=[f64::NAN;2];
        problem.island(0.0,&state,&mut derivative);
        assert_eq!(derivative,[0.0,0.0]);
        // This is an exact constant equilibrium with Gu=[0,1], supply=loss=1.
        // Its observed outputs [L+1,L] both round to L; even a compensated
        // input dot those rounded outputs would report zero supplied power.
        assert_eq!(system.output(&state),[large,large]);
        let orbit=HbOrbit {
            omega:input.omega_rad_s(),
            coeffs:state.into_iter().map(|x|vec![C64::from_re(x)]).collect(),
            residual_trace:Vec::new(),residual:f64::NAN,
            aft_samples:0,sampling_defect:f64::NAN,
        };
        for samples in [32,64] {
            let audit=audit_forced_orbit(&system,&input,&orbit,samples).unwrap();
            assert_eq!(audit.mean_supplied_power_w,1.0);
            assert_eq!(audit.mean_resistance_loss_w,1.0);
            assert_eq!(audit.cycle_power_residual_w,0.0);
            assert_eq!(audit.relative_power_residual,0.0);
        }
    }

    #[test]
    fn g1_cycle_loss_retains_a_weak_collective_resistance_direction() {
        let identity=vec![1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0];
        let system=PortHamiltonian::new(3,3,vec![0.0;9],vec![1.0;9],identity.clone(),
            Box::new(QuadraticStorage::new(identity,3).unwrap())).unwrap();
        let large=2.0_f64.powi(54);
        let input=PeriodicPortInput::new(1.0,vec![vec![C64::ONE];3]).unwrap();
        for state in [[large,1.0,-large],[1.0,large,-large],[-large,1.0,large]] {
            // R=b*b^T with b=[1,1,1]. Independently, (b dot state)^2=1,
            // while the authored input balances each resistance flow exactly.
            let orbit=HbOrbit {
                omega:1.0,
                coeffs:state.into_iter().map(|x|vec![C64::from_re(x)]).collect(),
                residual_trace:Vec::new(),residual:0.0,aft_samples:0,sampling_defect:0.0,
            };
            for samples in [32,64,128] {
                let audit=audit_forced_orbit(&system,&input,&orbit,samples).unwrap();
                assert_eq!(audit.mean_supplied_power_w,1.0);
                assert_eq!(audit.mean_resistance_loss_w,1.0);
                assert_eq!(audit.cycle_power_residual_w,0.0);
                assert_eq!(audit.relative_power_residual,0.0);
            }
        }
    }


    #[test]
    fn forced_cycle_audit_closes_original_source_and_resistance_work() {
        let (m,k,d,w,f)=(2.0_f64,1000.0_f64,4.0_f64,15.0_f64,3.0_f64);
        let system=fs_phs::mass_spring_damper(m,k,d).unwrap();
        let input=PeriodicPortInput::new(w,vec![vec![C64::ZERO,C64::from_re(0.5*f)]]).unwrap();
        let expected=f/det::sqrt((k-m*w*w).powi(2)+(d*w).powi(2));
        let orbit=solve_forced_harmonic_balance(&system,&input,expected,&hb(3)).unwrap();
        let a=audit_forced_orbit(&system,&input,&orbit,256).unwrap();
        assert!(a.mean_supplied_power_w>0.0 && a.mean_resistance_loss_w>0.0);
        assert!(a.relative_power_residual<2e-11,"{a:?}");
        let b=audit_forced_orbit(&system,&input,&orbit,512).unwrap();
        assert!((a.mean_supplied_power_w-b.mean_supplied_power_w).abs()<1e-12);
        assert!(a.maximum_stored_energy_j>a.minimum_stored_energy_j);
    }

    #[test]
    fn maxwell_cycle_audit_counts_internal_relaxation_loss_without_extra_law() {
        let (k,kb,tau,w,f)=(60_000.0_f64,50_000.0_f64,1.0/300.0,240.0_f64,2.0_f64);
        let system=fs_phs::mass_spring_damper(1.0,k,0.0).unwrap()
            .with_relaxation_branches(vec![RelaxationBranch {
                projection:vec![1.0,0.0],stiffness:kb,relaxation_time_s:tau,
            }]).unwrap();
        let input=PeriodicPortInput::new(w,vec![vec![C64::ZERO,C64::from_re(0.5*f)]]).unwrap();
        let orbit=solve_forced_harmonic_balance(&system,&input,1e-4,&hb(3)).unwrap();
        let audit=audit_forced_orbit(&system,&input,&orbit,512).unwrap();
        assert!(audit.mean_resistance_loss_w>0.0);
        assert!(audit.relative_power_residual<5e-10,"{audit:?}");
    }

    #[test]
    fn forced_linear_oscillator_matches_closed_form_amplitude() {
        let (m,k,d,w,f) = (2.0_f64, 1000.0_f64, 4.0_f64, 15.0_f64, 3.0_f64);
        let system = fs_phs::mass_spring_damper(m,k,d).unwrap();
        let input = PeriodicPortInput::new(w, vec![vec![C64::ZERO,C64::from_re(0.5*f)]]).unwrap();
        let expected = f / det::sqrt((k-m*w*w).powi(2)+(d*w).powi(2));
        let orbit = solve_forced_harmonic_balance(&system,&input,expected,&hb(3)).unwrap();
        assert!((orbit.first_harmonic_amplitude(0)-expected).abs() < 2e-10*expected);
        assert!((orbit.first_harmonic_amplitude(1)-w*expected*m).abs() < 2e-10*w*expected*m);
        assert!(orbit.coeffs.iter().all(|row| row.iter().skip(2).all(|z|z.abs()<1e-10)));
    }

    #[test]
    fn g0_unresolved_forcing_refuses_before_storage_evaluation() {
        use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};

        struct CountedStorage(Arc<AtomicUsize>);
        impl Storage for CountedStorage {
            fn hamiltonian(&self, x: &[f64]) -> f64 {
                self.0.fetch_add(1, Ordering::Relaxed);
                0.5 * x[0] * x[0]
            }
            fn gradient(&self, x: &[f64], out: &mut [f64]) {
                self.0.fetch_add(1, Ordering::Relaxed);
                out[0] = x[0];
            }
        }
        let calls = Arc::new(AtomicUsize::new(0));
        let system = PortHamiltonian::new(1, 2, vec![0.0], vec![1.0], vec![1.0, 1.0],
            Box::new(CountedStorage(calls.clone()))).unwrap();
        calls.store(0, Ordering::Relaxed);
        // At N=1 the former eight-point AFT grid folded harmonics 7 and 9
        // into the fundamental. Check every port and both complex components;
        // even a subnormal authored coefficient must not be silently removed.
        for (harmonic, coefficient) in [
            (7, C64::from_re(0.5)),
            (9, C64::new(0.0, -0.5)),
            (MAX_INPUT_HARMONICS, C64::from_re(f64::from_bits(1))),
        ] {
            let mut second = vec![C64::ZERO; harmonic + 1];
            second[harmonic] = coefficient;
            let input = PeriodicPortInput::new(1.0,
                vec![vec![C64::ZERO; harmonic + 1], second]).unwrap();
            assert_eq!(
                solve_forced_harmonic_balance(&system, &input, 1.0, &hb(1)).unwrap_err(),
                PhsOrbitError::Invalid {
                    what: "nonzero forcing harmonic exceeds the retained HB truncation",
                },
            );
            assert_eq!(calls.load(Ordering::Relaxed), 0);
            assert_eq!(
                solve_forced_harmonic_balance_adaptive(
                    &system, &input, 1.0, &hb(1), &refinement(1, 7),
                ).unwrap_err(),
                PhsOrbitError::Invalid {
                    what: "adaptive maximum must exceed the starting and forcing harmonic orders",
                },
            );
            assert_eq!(calls.load(Ordering::Relaxed), 0);
        }
    }

    #[test]
    fn g3_zero_padding_preserves_the_same_forced_solution() {
        let system = fs_phs::mass_spring_damper(2.0, 1000.0, 4.0).unwrap();
        let compact = vec![C64::from_re(2.0), C64::new(1.5, -0.5)];
        let mut padded = compact.clone();
        padded.resize(MAX_INPUT_HARMONICS + 1, C64::new(-0.0, 0.0));
        let a = PeriodicPortInput::new(15.0, vec![compact]).unwrap();
        let b = PeriodicPortInput::new(15.0, vec![padded]).unwrap();
        assert_eq!(b.harmonic_count(), MAX_INPUT_HARMONICS);
        let a = solve_forced_harmonic_balance(&system, &a, 0.01, &hb(1)).unwrap();
        let b = solve_forced_harmonic_balance(&system, &b, 0.01, &hb(1)).unwrap();
        assert_eq!(a.coeffs, b.coeffs);
        assert_eq!(a.residual_trace, b.residual_trace);
        assert_eq!(a.omega.to_bits(), b.omega.to_bits());
    }

    #[test]
    fn g1_dc_and_higher_harmonics_match_each_analytic_phasor_and_power() {
        let (m, k, d, w) = (2.0, 1000.0, 4.0, 15.0);
        let system = fs_phs::mass_spring_damper(m, k, d).unwrap();
        // No fundamental: the real response contains a static displacement,
        // second and fourth harmonics with independently authored phases.
        let coefficients = vec![C64::from_re(2.0), C64::ZERO,
            C64::new(1.5, -0.75), C64::ZERO, C64::new(-0.4, 0.3)];
        let input = PeriodicPortInput::new(w, vec![coefficients.clone()]).unwrap();
        let orbit = solve_forced_harmonic_balance(&system, &input, 0.01, &hb(4)).unwrap();
        let mut expected_loss = 0.0;
        for (n, force) in coefficients.iter().enumerate() {
            let wn = n as f64 * w;
            let displacement = *force / C64::new(k - m * wn * wn, d * wn);
            let momentum = displacement * C64::new(0.0, m * wn);
            assert!((orbit.coeffs[0][n] - displacement).abs()
                < 1e-10 * displacement.abs().max(1e-4));
            assert!((orbit.coeffs[1][n] - momentum).abs()
                < 1e-10 * momentum.abs().max(1e-4));
            if n > 0 { expected_loss += 2.0 * d * wn * wn * displacement.abs().powi(2); }
        }
        let audit = audit_forced_orbit(&system, &input, &orbit, 256).unwrap();
        assert!((audit.mean_resistance_loss_w - expected_loss).abs() < 1e-9 * expected_loss);
        assert!(audit.relative_power_residual < 1e-9);
    }

    #[test]
    fn adaptive_forcing_resolves_authored_bandwidth_and_ignores_zero_padding() {
        let (m, k, d, w) = (2.0, 1000.0, 4.0, 15.0);
        let system = fs_phs::mass_spring_damper(m, k, d).unwrap();
        let coefficients = vec![C64::from_re(2.0), C64::ZERO,
            C64::new(1.5, -0.75), C64::ZERO, C64::new(-0.4, 0.3)];
        let input = PeriodicPortInput::new(w, vec![coefficients.clone()]).unwrap();
        let mut padded = coefficients.clone();
        padded.resize(MAX_INPUT_HARMONICS + 1, C64::new(-0.0, 0.0));
        let padded = PeriodicPortInput::new(w, vec![padded]).unwrap();
        assert_eq!(input.highest_nonzero_harmonic(), 4);
        assert_eq!(padded.highest_nonzero_harmonic(), 4);
        assert_eq!(padded.harmonic_count(), MAX_INPUT_HARMONICS);
        // N=1 is deliberately below the authored fourth harmonic. The adaptive
        // wrapper admits its bandwidth, then independently solves a larger order.
        let result = solve_forced_harmonic_balance_adaptive(
            &system, &input, 0.01, &hb(1), &refinement(2, 17),
        ).unwrap();
        let repeated = solve_forced_harmonic_balance_adaptive(
            &system, &padded, 0.01, &hb(1), &refinement(2, 17),
        ).unwrap();
        let orbit = &result.orbit;
        assert!(orbit.coeffs[0].len() > coefficients.len());
        assert_eq!(orbit.coeffs, repeated.orbit.coeffs);
        assert_eq!(orbit.residual_trace, repeated.orbit.residual_trace);
        assert_eq!(orbit.omega.to_bits(), repeated.orbit.omega.to_bits());
        for (n, force) in coefficients.iter().enumerate() {
            let wn = n as f64 * w;
            let displacement = *force / C64::new(k - m * wn * wn, d * wn);
            let momentum = displacement * C64::new(0.0, m * wn);
            assert!((orbit.coeffs[0][n] - displacement).abs() < 1e-12);
            assert!((orbit.coeffs[1][n] - momentum).abs() < 1e-11);
        }
        let audit = audit_forced_orbit(&system, &input, orbit, 512).unwrap();
        assert!(audit.relative_power_residual < 1e-9);
    }

    #[test]
    fn scaled_physical_forcing_preserves_tiny_multiharmonic_phasors_and_power() {
        let (mass, stiffness, damping, omega) = (2.0, 1000.0, 4.0, 15.0);
        let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
        let force = vec![
            C64::ZERO,
            C64::ZERO,
            C64::new(1.5e-12, -0.75e-12),
            C64::new(-0.4e-12, 0.3e-12),
        ];
        let input = PeriodicPortInput::new(omega, vec![force.clone()]).unwrap();
        let scales = HbScales {
            states: vec![1e-14, 1e-12],
            equations: vec![1e-12, 1e-10],
            angular_frequency: omega,
        };
        let budget = HbBudget {
            harmonics: 3,
            tolerance: 1e-11,
            ..HbBudget::default()
        };
        let fixed = solve_forced_harmonic_balance_scaled_seeded(
            &system,
            &input,
            vec![0.0; 14],
            &budget,
            &scales,
        )
        .unwrap();
        let adaptive = solve_forced_harmonic_balance_adaptive_scaled_seeded(
            &system,
            &input,
            vec![0.0; 14],
            &budget,
            &HbRefinementBudget {
                max_harmonics: 7,
                residual_tolerances: vec![1e-20, 1e-18],
                state_tolerances: vec![1e-22, 1e-20],
                frequency_tolerance: 1e-10,
            },
            &scales,
        )
        .unwrap();
        assert!(fixed.scaled_residual < budget.tolerance);
        assert_eq!(
            adaptive
                .refinement_trace
                .iter()
                .map(|r| r.harmonics)
                .collect::<Vec<_>>(),
            vec![3, 7]
        );
        for orbit in [&fixed.orbit, &adaptive.orbit] {
            let mut expected_loss = 0.0;
            for (n, force) in force.iter().enumerate() {
                let w = n as f64 * omega;
                let displacement = *force / C64::new(stiffness - mass * w * w, damping * w);
                let momentum = displacement * C64::new(0.0, mass * w);
                assert!((orbit.coeffs[0][n] - displacement).abs() < 1e-22);
                assert!((orbit.coeffs[1][n] - momentum).abs() < 1e-20);
                expected_loss += 2.0 * damping * w * w * displacement.abs().powi(2);
            }
            let audit = audit_forced_orbit(&system, &input, orbit, 256).unwrap();
            assert!((audit.mean_resistance_loss_w / expected_loss - 1.0).abs() < 1e-7);
            assert!(audit.relative_power_residual < 1e-7);
        }
    }

    #[test]
    fn adaptive_backbone_matches_independent_energy_period_integral() {
        let (k, alpha, mass, amplitude) = (4.0, 10.0, 2.0, 0.6);
        let system = PortHamiltonian::new(
            2,
            0,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0; 4],
            vec![],
            Box::new(Duffing {
                k,
                alpha,
                inv_mass: 1.0 / mass,
            }),
        )
        .unwrap();
        let omega_guess = det::sqrt((k + 0.75 * alpha * amplitude * amplitude) / mass);
        let result = solve_backbone_harmonic_balance_adaptive(
            &system,
            &[],
            amplitude,
            omega_guess,
            &hb(1),
            &refinement(2, 31),
        )
        .unwrap();
        let orbit = &result.orbit;
        assert!(orbit.coeffs[0].len() > 4);
        assert!(orbit.coeffs[0][3].abs() > 1e-4);
        assert!((orbit.first_harmonic_amplitude(0) - amplitude).abs() < 1e-9);
        let initial = [orbit.sample(0, 0.0), orbit.sample(1, 0.0)];
        let energy = system.hamiltonian(&initial);
        let mut energy_defect: f64 = 0.0;
        for j in 0..256 {
            let phase = core::f64::consts::TAU * j as f64 / 256.0;
            let state = [orbit.sample(0, phase), orbit.sample(1, phase)];
            energy_defect = energy_defect.max((system.hamiltonian(&state) - energy).abs());
        }
        assert!(energy_defect < 1e-7 * energy, "energy defect {energy_defect}");
        // Independent quadrature of dt = sqrt(m) dphi /
        // sqrt(k + alpha*A^2*(1 + sin(phi)^2)/2), with A inferred from
        // the conserved energy rather than equated to the fundamental amplitude.
        let turning_squared = 4.0 * energy
            / (k + det::sqrt(k * k + 4.0 * alpha * energy));
        let panels = 4096;
        let step = core::f64::consts::FRAC_PI_2 / panels as f64;
        let mut integral = 0.0;
        for j in 0..=panels {
            let sine = det::sin(j as f64 * step);
            let value = 1.0 / det::sqrt(k + 0.5 * alpha * turning_squared * (1.0 + sine * sine));
            let weight = if j == 0 || j == panels { 1.0 }
                else if j % 2 == 0 { 2.0 } else { 4.0 };
            integral += weight * value;
        }
        let period = 4.0 * det::sqrt(mass) * step * integral / 3.0;
        let reference_omega = core::f64::consts::TAU / period;
        assert!((orbit.omega - reference_omega).abs() / reference_omega < 1e-8);
    }

    #[test]
    fn maxwell_memory_is_retained_in_the_periodic_balance() {
        let (k,kb,tau,w,f) = (60_000.0_f64,50_000.0_f64,1.0/300.0,240.0_f64,2.0_f64);
        let system = fs_phs::mass_spring_damper(1.0,k,0.0).unwrap()
            .with_relaxation_branches(vec![RelaxationBranch {
                projection:vec![1.0,0.0], stiffness:kb, relaxation_time_s:tau,
            }]).unwrap();
        let input = PeriodicPortInput::new(w,vec![vec![C64::ZERO,C64::from_re(0.5*f)]]).unwrap();
        let iw = C64::new(0.0,w);
        let modulus = (iw.scale(tau)/(C64::ONE+iw.scale(tau))).scale(kb);
        let denominator = C64::from_re(k-w*w)+modulus;
        let expected = f/denominator.abs();
        let orbit = solve_forced_harmonic_balance(&system,&input,expected,&hb(3)).unwrap();
        assert_eq!(orbit.coeffs.len(),3);
        assert!((orbit.first_harmonic_amplitude(0)-expected).abs()<5e-9*expected);
        assert!(orbit.first_harmonic_amplitude(2)>0.0);
    }

    struct Duffing { k:f64, alpha:f64, inv_mass:f64 }
    impl Storage for Duffing {
        fn hamiltonian(&self,x:&[f64])->f64 {
            0.5*self.k*x[0]*x[0]+0.25*self.alpha*x[0].powi(4)+0.5*self.inv_mass*x[1]*x[1]
        }
        fn gradient(&self,x:&[f64],out:&mut[f64]) {
            out[0]=self.k*x[0]+self.alpha*x[0].powi(3);
            out[1]=self.inv_mass*x[1];
        }
    }
    #[test]
    fn nonlinear_backbone_keeps_the_original_cubic_storage() {
        let (k,alpha,m,a)=(400.0_f64,100.0_f64,1.0_f64,0.1_f64);
        let system=PortHamiltonian::new(2,0,vec![0.,1.,-1.,0.],vec![0.;4],vec![],
            Box::new(Duffing{k,alpha,inv_mass:1.0/m})).unwrap();
        let first_order=det::sqrt((k+0.75*alpha*a*a)/m);
        let orbit=solve_backbone_harmonic_balance(&system,&[],a,first_order,&hb(7)).unwrap();
        assert!((orbit.first_harmonic_amplitude(0)-a).abs()<1e-10);
        assert!((orbit.omega-first_order).abs()/first_order<2e-4);
        assert!(orbit.coeffs[0][3].abs()>0.0, "cubic storage must generate a third harmonic");
    }

    #[test]
    fn input_shape_and_source_failure_are_refusals_not_partial_orbits() {
        assert!(PeriodicPortInput::new(10.0,vec![vec![C64::new(0.0,1.0)]]).is_err());
        let system=fs_phs::mass_spring_damper(1.0,100.0,0.2).unwrap();
        let wrong=PeriodicPortInput::new(10.0,vec![]).unwrap();
        assert!(solve_forced_harmonic_balance(&system,&wrong,0.1,&hb(3)).is_err());
        let ok=PeriodicPortInput::new(10.0,vec![vec![C64::ZERO,C64::from_re(0.5)]]).unwrap();
        let a=solve_forced_harmonic_balance(&system,&ok,0.1,&hb(3)).unwrap();
        let b=solve_forced_harmonic_balance(&system,&ok,0.1,&hb(3)).unwrap();
        assert_eq!(a.omega.to_bits(),b.omega.to_bits());
        assert_eq!(a.coeffs,b.coeffs);
    }
}
