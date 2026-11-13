//! Full-state Floquet evidence for the SAME original periodic response base.

use super::{
    C64, HbOrbit, PhsPortSidebandError, PhsPortSidebandSource, PreparedPhsPortSidebands,
    SidebandScales, SourceProblem,
};
use core::cell::RefCell;
use fs_la::gemm::gemm_scalar_checked;
use fs_math::{dd::Dd, det};
use fs_orbit::{
    OrbitProblem, ShootBudget, TimeDomainRhs, TimeDomainStateJacobian, VariationalShootControls,
    VariationalShootError, solve_shooting_variational,
};
use fs_spectral::floquet::{
    FloquetClassification, FloquetControls, FloquetError, FloquetSpectrum,
    assess_monodromy_with_multipliers,
};

mod transient;
pub use transient::{
    PhsBaseTransientControls, PhsBaseTransientError, PhsBaseTransientLevel, PhsBaseTransientReport,
};

/// Independent temporal, spectral and TOTAL stability-operation limits.
#[derive(Clone, Debug)]
pub struct PhsBaseStabilityControls {
    /// Three grids use N, 2N and 4N steps over the ORIGINAL HB base period.
    /// The first must resolve the complete retained base waveform at Nyquist.
    pub initial_steps_per_period: usize,
    /// Frobenius norm allowance for adjacent complete S^-1*M*S maps. This is a
    /// numerical comparison in the disclosed dimensionless state chart.
    pub monodromy_change_tolerance: f64,
    /// Absolute change allowance for both full and effective spectral radii.
    pub spectral_radius_change_tolerance: f64,
    /// Original state scales and per-map native limits. Each internal linear
    /// query uses exactly one map, a zero seed, and a Newton cap of one. The
    /// autonomous excursion floor applies to the ORIGINAL HB waveform.
    pub variational: VariationalShootControls,
    /// Per-spectrum complete basis, residual and autonomous phase controls.
    /// Multipliers come from the same map's existing QR; no second QR runs.
    pub floquet: FloquetControls,
    /// TOTAL new live numerical storage for this call: stage-time cache, native
    /// workspaces, all three retained maps/spectra, controls and phase scratch.
    /// Existing prepared source/tangent storage and caller arrays are excluded;
    /// their already admitted caps ADD to this bound. Allocator metadata and
    /// arbitrary source-callback internal allocations remain excluded.
    pub maximum_workspace_bytes: usize,
    /// TOTAL conservative arithmetic admission for the three one-map solves,
    /// original-base/tangent sampling and three spectral allowances. The latter
    /// reserves 3*floquet.spectrum.maximum_work; this is not measured FLOPs.
    /// Native callback internals remain source-defined. No cap resets per grid.
    pub maximum_work: u64,
}

/// One temporal resolution of the original fixed-base variational equations.
#[derive(Clone, Debug)]
pub struct PhsBaseStabilityLevel {
    pub steps_per_period: usize,
    /// Complete row-major S^-1*M*S in ORIGINAL state order, including memory.
    /// This is the RK4 map of y'=A_original(t,x_HB(t))*y, not a corrected orbit.
    pub monodromy: Vec<f64>,
    /// Every multiplier/right vector and original-map residual; vectors use
    /// the same disclosed scaled state chart. An identified phase is retained.
    pub spectrum: FloquetSpectrum,
    /// Full Frobenius change divided by its separate declared allowance.
    pub monodromy_change_ratio: Option<f64>,
    pub spectral_radius_change_ratio: Option<f64>,
    /// Transverse for a resolved autonomous phase, otherwise full for forcing.
    /// None when this is the first grid or either phase assessment is unresolved.
    pub effective_radius_change_ratio: Option<f64>,
    /// Maximum original HB state excursion from phase zero, divided coordinate
    /// by coordinate by state_scales, measured at this map's actual RK4 stages.
    pub maximum_scaled_base_excursion: f64,
    /// The shared linear-map owner requests four RHS and Jacobian stages per
    /// step. A same-time cache avoids repeated ORIGINAL Hessian/feedback work.
    pub linear_rhs_evaluations: usize,
    pub linear_jacobian_evaluations: usize,
    pub original_jacobian_evaluations: usize,
}

/// Local numerical stability evidence bound to the prepared physical base.
#[derive(Clone, Debug)]
pub struct PhsBaseStabilityReport {
    pub omega_rad_s: f64,
    pub period_s: f64,
    pub autonomous: bool,
    pub state_scales: Vec<f64>,
    /// S^-1*f_original(x_HB(0)), complete and in original state ordering.
    /// None for forcing, even at zero source amplitude or with a unit multiplier.
    pub scaled_phase_flow: Option<Vec<C64>>,
    /// Original sampled balance/tangent checks from the existing preparation.
    pub base_residual: f64,
    pub tangent_sampling_defect: f64,
    /// Requires all three spectral classifications and BOTH complete-map/radius
    /// comparisons to agree. Failed numerical resolution remains Unresolved.
    pub classification: FloquetClassification,
    pub levels: Vec<PhsBaseStabilityLevel>,
    /// Finest measured full/transverse log-radius per second, when defined.
    /// Its availability does not override an Unresolved classification.
    pub dominant_growth_rate_per_s: Option<f64>,
    pub original_jacobian_evaluations: usize,
    /// Conservative admitted operation bounds, useful for repeatable budgets.
    pub admitted_workspace_bytes: usize,
    pub admitted_work: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PhsBaseStabilityError {
    Invalid(&'static str),
    Budget(&'static str),
    Preparation(PhsPortSidebandError),
    /// The original physical flow or declared complete tangent refused.
    SourceRefused {
        time: f64,
    },
    Variational(VariationalShootError),
    Spectrum(FloquetError),
    Cancelled,
}
impl core::fmt::Display for PhsBaseStabilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical periodic-base stability refused: {self:?}")
    }
}
impl std::error::Error for PhsBaseStabilityError {}
impl From<PhsPortSidebandError> for PhsBaseStabilityError {
    fn from(error: PhsPortSidebandError) -> Self {
        match error {
            PhsPortSidebandError::Cancelled => Self::Cancelled,
            other => Self::Preparation(other),
        }
    }
}
impl From<VariationalShootError> for PhsBaseStabilityError {
    fn from(error: VariationalShootError) -> Self {
        match error {
            VariationalShootError::Cancelled => Self::Cancelled,
            other => Self::Variational(other),
        }
    }
}
impl From<FloquetError> for PhsBaseStabilityError {
    fn from(error: FloquetError) -> Self {
        match error {
            FloquetError::Cancelled => Self::Cancelled,
            other => Self::Spectrum(other),
        }
    }
}
type Result<T> = core::result::Result<T, PhsBaseStabilityError>;

impl PhsPortSidebandSource<'_> {
    /// Prepare this original response base once, then assess its fixed-base map.
    /// The preparation's existing core/adapter allowances are ADDITIONAL to the
    /// stability operation's total limits. This does not perform a carrier solve.
    /// # Errors
    /// Original preparation refusals and [`PreparedPhsPortSidebands::analyze_base_stability`]
    /// errors. The source and its original base remain unchanged and reusable.
    pub fn analyze_base_stability(
        &mut self,
        scales: &SidebandScales,
        controls: &PhsBaseStabilityControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhsBaseStabilityReport> {
        let mut cancelled = cancelled;
        let mut prepared = self.prepare(scales, &mut cancelled)?;
        prepared.analyze_base_stability(controls, cancelled)
    }
}

impl PreparedPhsPortSidebands<'_, '_> {
    /// Assess the stability of the SAME fixed original HB base used by probes.
    ///
    /// This complete pHS source has an explicitly realized first-order ODE,
    /// including every constitutive/controller-memory coordinate. Its exact
    /// original tangent supplies the linear periodic problem
    /// y'=A_original(t,x_HB(t))*y. The existing variational RK4 owner integrates
    /// its zero primal and consistent tangent together. Its exact zero closure
    /// needs one period map; no physical waveform, period or phase is corrected.
    /// A same-time cache shares original A between the RHS/Jacobian and repeated
    /// RK4 stages. No cached Fourier truncation replaces the original tangent.
    ///
    /// Every actual map retains all coordinates. Existing QR multipliers go to
    /// the shared spectral owner for full original-map eigenpair residuals and
    /// complete-basis conditioning, without repeating QR. For autonomous cycles,
    /// the ORIGINAL full scaled phase flow must satisfy M*g=g and identify one
    /// isolated near-+1 multiplier before a transverse radius is available.
    /// The phase mode stays in the spectrum. Multiple neutral modes remain
    /// unresolved; no coordinate deletion, amplitude constraint or damping floor
    /// is introduced. Forced bases never remove a phase, including zero input.
    ///
    /// Maps at N, 2N and 4N share the same physical period, base and state chart.
    /// Both adjacent full-map and full/effective-radius comparisons must pass,
    /// and all three guarded spectral classifications must agree. A large or
    /// even identical finite forced-response magnitude alone says nothing about
    /// this classification. Analyze once per base and reuse its carrier factors.
    ///
    /// The prepared base balance and this measured temporal ladder are numerical
    /// evidence, not a continuous-orbit existence proof, infinite-tail enclosure,
    /// forward spectral-error bound, global attraction or hybrid-event analysis.
    /// A response about an unstable base remains a mathematical particular
    /// response; this method exposes its local growth instead of relabeling it.
    ///
    /// # Errors
    /// Invalid complete controls/scales, exhausted TOTAL budgets, stationary or
    /// unresolved autonomous excursion, original source/tangent refusal, native
    /// map/spectral refusal, or cancellation. No partial report escapes. The
    /// prepared source is unchanged and remains reusable after any refusal.
    pub fn analyze_base_stability(
        &mut self,
        controls: &PhsBaseStabilityControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsBaseStabilityReport> {
        poll(&mut cancelled)?;
        let source = self.source;
        let n = source.problem.dim();
        let (grids, period, admitted_bytes, admitted_work) = admit(source, controls)?;
        let scales = &controls.variational.state_scales;
        let initial: Vec<_> = (0..n).map(|i| source.orbit.sample(i, 0.0)).collect();
        if initial.iter().any(|x| !x.is_finite()) {
            return Err(PhsBaseStabilityError::Invalid(
                "original phase-zero state is not representable",
            ));
        }
        let autonomous = source.problem.autonomous();
        let phase = if autonomous {
            let mut flow = vec![f64::NAN; n];
            source.problem.island(0.0, &initial, &mut flow);
            poll(&mut cancelled)?;
            if flow.iter().any(|x| !x.is_finite()) {
                return Err(PhsBaseStabilityError::SourceRefused { time: 0.0 });
            }
            let mut phase = Vec::with_capacity(n);
            for (&value, &scale) in flow.iter().zip(scales) {
                let normalized = value / scale;
                if !normalized.is_finite() || (value != 0.0 && normalized == 0.0) {
                    return Err(PhsBaseStabilityError::Invalid(
                        "complete scaled physical phase flow is not representable",
                    ));
                }
                phase.push(C64::from_re(normalized));
            }
            if phase.iter().all(|z| *z == C64::ZERO) {
                return Err(PhsBaseStabilityError::Invalid(
                    "autonomous base must have a nonzero original phase flow",
                ));
            }
            Some(phase)
        } else {
            None
        };
        let problem = LinearizedBase {
            source: &source.problem,
            base: source.orbit,
            scales,
            cache: RefCell::new(TangentCache {
                time: None,
                state: vec![0.0; n],
                initial,
                tangent: vec![0.0; n * n],
                evaluations: 0,
                excursion: 0.0,
                refused: None,
            }),
        };
        let zero = vec![0.0; n];
        let mut variational = controls.variational.clone();
        variational.maximum_period_maps = 1;
        let mut levels: Vec<PhsBaseStabilityLevel> = Vec::with_capacity(3);
        let mut total_jacobians = 0;
        for steps in grids {
            poll(&mut cancelled)?;
            problem.reset();
            let budget = ShootBudget {
                steps_per_period: steps,
                max_newton: 1,
                tolerance: f64::EPSILON,
            };
            let solved = solve_shooting_variational(
                &problem,
                &zero,
                period,
                &budget,
                &variational,
                &mut cancelled,
            );
            poll(&mut cancelled)?;
            let (original_jacobians, excursion, refused) = {
                let cache = problem.cache.borrow();
                (cache.evaluations, cache.excursion, cache.refused)
            };
            if let Some(time) = refused {
                return Err(PhsBaseStabilityError::SourceRefused { time });
            }
            let solved = solved?;
            if solved.period_maps != 1
                || solved.normalized_residual != 0.0
                || solved.orbit.x0.iter().any(|&x| x != 0.0)
                || solved.orbit.period.to_bits() != period.to_bits()
            {
                return Err(PhsBaseStabilityError::Invalid(
                    "linearized zero solution must retain exactly one original-period map",
                ));
            }
            if autonomous && excursion < variational.minimum_autonomous_excursion {
                return Err(PhsBaseStabilityError::Invalid(
                    "original autonomous base excursion is unresolved in physical state scales",
                ));
            }
            let matrix: Vec<_> = solved.monodromy.iter().copied().map(C64::from_re).collect();
            let multipliers: Vec<_> = solved
                .orbit
                .multipliers
                .iter()
                .map(|&(re, im)| C64::new(re, im))
                .collect();
            let spectrum = assess_monodromy_with_multipliers(
                &matrix,
                n,
                period,
                &multipliers,
                phase.as_deref(),
                controls.floquet,
                &mut cancelled,
            )?;
            let previous = levels.last();
            let map_change = previous.map(|old| {
                matrix_change(
                    &old.monodromy,
                    &solved.monodromy,
                    controls.monodromy_change_tolerance,
                )
            });
            let radius_change = previous.map(|old| {
                scalar_change(
                    old.spectrum.spectral_radius,
                    spectrum.spectral_radius,
                    controls.spectral_radius_change_tolerance,
                )
            });
            let effective_change = previous.and_then(|old| {
                old.spectrum
                    .effective_spectral_radius
                    .zip(spectrum.effective_spectral_radius)
                    .map(|(a, b)| scalar_change(a, b, controls.spectral_radius_change_tolerance))
            });
            total_jacobians += original_jacobians;
            levels.push(PhsBaseStabilityLevel {
                steps_per_period: steps,
                monodromy: solved.monodromy,
                spectrum,
                monodromy_change_ratio: map_change,
                spectral_radius_change_ratio: radius_change,
                effective_radius_change_ratio: effective_change,
                maximum_scaled_base_excursion: excursion,
                linear_rhs_evaluations: solved.rhs_evaluations,
                linear_jacobian_evaluations: solved.jacobian_evaluations,
                original_jacobian_evaluations: original_jacobians,
            });
        }
        let temporal = levels.iter().skip(1).all(|level| {
            level.monodromy_change_ratio.is_some_and(|r| r <= 1.0)
                && level.spectral_radius_change_ratio.is_some_and(|r| r <= 1.0)
                && level
                    .effective_radius_change_ratio
                    .is_some_and(|r| r <= 1.0)
        });
        let classification = if temporal
            && levels
                .iter()
                .all(|level| level.spectrum.classification == FloquetClassification::Stable)
        {
            FloquetClassification::Stable
        } else if temporal
            && levels
                .iter()
                .all(|level| level.spectrum.classification == FloquetClassification::Unstable)
        {
            FloquetClassification::Unstable
        } else {
            FloquetClassification::Unresolved
        };
        let growth = levels[2].spectrum.dominant_growth_rate_per_s;
        poll(&mut cancelled)?;
        Ok(PhsBaseStabilityReport {
            omega_rad_s: source.orbit.omega,
            period_s: period,
            autonomous,
            state_scales: scales.clone(),
            scaled_phase_flow: phase,
            base_residual: self.base_residual(),
            tangent_sampling_defect: self.tangent_sampling_defect(),
            classification,
            levels,
            dominant_growth_rate_per_s: growth,
            original_jacobian_evaluations: total_jacobians,
            admitted_workspace_bytes: admitted_bytes,
            admitted_work,
        })
    }
}

struct TangentCache {
    time: Option<u64>,
    state: Vec<f64>,
    initial: Vec<f64>,
    tangent: Vec<f64>,
    evaluations: usize,
    excursion: f64,
    refused: Option<f64>,
}
struct LinearizedBase<'a, 'model> {
    source: &'a SourceProblem<'model>,
    base: &'a HbOrbit,
    scales: &'a [f64],
    cache: RefCell<TangentCache>,
}
impl LinearizedBase<'_, '_> {
    fn reset(&self) {
        let mut cache = self.cache.borrow_mut();
        cache.time = None;
        cache.evaluations = 0;
        cache.excursion = 0.0;
        cache.refused = None;
    }
    fn tangent(&self, time: f64) -> bool {
        let mut cache = self.cache.borrow_mut();
        if cache.time == Some(time.to_bits()) {
            return true;
        }
        let theta = self.base.omega * time;
        if !time.is_finite() || !theta.is_finite() {
            return false;
        }
        for i in 0..cache.state.len() {
            let value = self.base.sample(i, theta);
            let change =
                (Dd::from_f64(value) - Dd::from_f64(cache.initial[i])).to_f64() / self.scales[i];
            if !value.is_finite() || !change.is_finite() {
                return false;
            }
            cache.state[i] = value;
            cache.excursion = cache.excursion.max(change.abs());
        }
        cache.evaluations += 1;
        let TangentCache { state, tangent, .. } = &mut *cache;
        tangent.fill(f64::NAN);
        if !self.source.state_jacobian(time, state, tangent)
            || tangent.iter().any(|x| !x.is_finite())
        {
            cache.refused = Some(time);
            return false;
        }
        cache.time = Some(time.to_bits());
        true
    }
    fn rhs(&self, time: f64, state: &[f64], out: &mut [f64]) {
        if !self.tangent(time) {
            out.fill(f64::NAN);
            return;
        }
        let cache = self.cache.borrow();
        out.fill(0.0);
        if gemm_scalar_checked(
            self.dim(),
            1,
            self.dim(),
            1.0,
            &cache.tangent,
            state,
            0.0,
            out,
        )
        .is_err()
        {
            out.fill(f64::NAN);
        }
    }
    fn jacobian(&self, time: f64, _: &[f64], out: &mut [f64]) -> bool {
        if out.len() != self.dim() * self.dim() || !self.tangent(time) {
            return false;
        }
        out.copy_from_slice(&self.cache.borrow().tangent);
        true
    }
}
impl OrbitProblem for LinearizedBase<'_, '_> {
    fn dim(&self) -> usize {
        self.source.dim()
    }
    // The linear coefficient clock is prescribed, even when its physical base
    // is autonomous. Phase neutrality is checked on the COMPLETE resulting map.
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        self.rhs(time, state, out);
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::rhs)
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        Some(Self::jacobian)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PhsBaseStabilityError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn admit(
    source: &PhsPortSidebandSource<'_>,
    controls: &PhsBaseStabilityControls,
) -> Result<([usize; 3], f64, usize, u64)> {
    let n = source.problem.dim();
    let v = &controls.variational;
    let finest =
        controls
            .initial_steps_per_period
            .checked_mul(4)
            .ok_or(PhsBaseStabilityError::Budget(
                "stability time-grid extent overflows",
            ))?;
    let period = core::f64::consts::TAU / source.orbit.omega;
    let width = source.orbit.coeffs[0].len();
    if n == 0
        || n > 256
        || n > v.maximum_dimension
        || v.state_scales.len() != n
        || v.state_scales.iter().any(|&x| !positive(x))
        || !positive(v.period_scale)
        || !positive(v.minimum_autonomous_excursion)
        || v.maximum_period_maps == 0
        || controls.initial_steps_per_period <= 2 * (width - 1)
        || !positive(period)
        || !positive(controls.monodromy_change_tolerance)
        || !positive(controls.spectral_radius_change_tolerance)
    {
        return Err(PhsBaseStabilityError::Invalid(
            "complete finite scales, resolved base grids and positive stability controls required",
        ));
    }
    if source.problem.autonomous()
        && source
            .orbit
            .coeffs
            .iter()
            .all(|row| row.iter().skip(1).all(|&z| z == C64::ZERO))
    {
        return Err(PhsBaseStabilityError::Invalid(
            "a constant autonomous state is not a nonstationary periodic base",
        ));
    }
    let f = controls.floquet;
    let s = f.spectrum;
    if n > s.maximum_dimension
        || s.svd_max_sweeps == 0
        || !positive(s.svd_tolerance)
        || s.svd_tolerance <= f64::EPSILON
        || s.svd_tolerance >= 1.0
        || !positive(s.maximum_relative_residual)
        || s.maximum_relative_residual >= 1.0
        || !positive(s.minimum_eigenbasis_reciprocal_condition)
        || s.minimum_eigenbasis_reciprocal_condition >= 1.0
        || !positive(f.unit_circle_guard)
        || f.unit_circle_guard >= 1.0
        || !positive(f.phase_multiplier_tolerance)
        || f.phase_multiplier_tolerance >= f.unit_circle_guard
        || !positive(f.maximum_phase_residual)
        || f.maximum_phase_residual >= 1.0
    {
        return Err(PhsBaseStabilityError::Invalid(
            "complete guarded Floquet and eigenbasis controls required",
        ));
    }
    let grids = [controls.initial_steps_per_period, finest / 2, finest];
    let (d, m) = (
        n as u128,
        source.problem.original().system.port_dim() as u128,
    );
    // Existing variational owner's exact ONE-map, ONE-Newton admission formula.
    let map_bytes = 384 * d * d + 1024 * d + 128;
    let mut work = 0_u128;
    for steps in grids {
        let requested = steps as u128 * (16 * d * d * d + 128 * d * d + 128 * d) + 8320 * d * d * d;
        if map_bytes > v.maximum_payload_bytes as u128 || requested > v.maximum_work as u128 {
            return Err(PhsBaseStabilityError::Budget(
                "one complete stability map exceeds its native per-map allowance",
            ));
        }
        work += requested;
    }
    let stages = 28 * controls.initial_steps_per_period as u128;
    if stages > usize::MAX as u128 {
        return Err(PhsBaseStabilityError::Budget(
            "stability evaluation counts are not representable",
        ));
    }
    // Count every stage conservatively although the same-time cache skips many
    // original callbacks. Include physical base synthesis, full Hessian/J-R/G
    // contractions, total feedback derivative, and normalized-map comparisons.
    let original_stage = 4 * d * d * d
        + 64 * d * d
        + 16 * d * d * m
        + 64 * d * m
        + 128 * d * width as u128
        + 128 * d;
    work = (stages + 1)
        .checked_mul(original_stage)
        .and_then(|extra| work.checked_add(extra))
        .and_then(|total| total.checked_add(3 * controls.floquet.spectrum.maximum_work as u128))
        .and_then(|total| total.checked_add(1024 * d * d + 4096 * d))
        .ok_or(PhsBaseStabilityError::Budget(
            "total stability arithmetic envelope overflows",
        ))?;
    // Three complete maps/right bases, one original-time cache, conversions,
    // zero seed, phase flow, cloned scales, bounded traces and level headers.
    let wrapper_bytes = 128 * d * d + 512 * d + 4096;
    let bytes =
        wrapper_bytes + map_bytes.max(controls.floquet.spectrum.maximum_dense_bytes as u128);
    if work > controls.maximum_work as u128 || bytes > controls.maximum_workspace_bytes as u128 {
        return Err(PhsBaseStabilityError::Budget(
            "three complete stability levels exceed the TOTAL operation allowance",
        ));
    }
    Ok((grids, period, bytes as usize, work as u64))
}

fn matrix_change(old: &[f64], new: &[f64], tolerance: f64) -> f64 {
    let mut norm = 0.0;
    for (&a, &b) in old.iter().zip(new) {
        let difference = (Dd::from_f64(a) - Dd::from_f64(b)).to_f64() / tolerance;
        norm = det::hypot(norm, difference);
    }
    if norm.is_finite() {
        norm
    } else {
        f64::INFINITY
    }
}
fn scalar_change(a: f64, b: f64, tolerance: f64) -> f64 {
    let change = ((Dd::from_f64(a) - Dd::from_f64(b)).to_f64() / tolerance).abs();
    if change.is_finite() {
        change
    } else {
        f64::INFINITY
    }
}
