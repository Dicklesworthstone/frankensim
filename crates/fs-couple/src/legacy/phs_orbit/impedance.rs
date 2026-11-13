//! Frequency-dependent physical loads on the complete nonlinear pHS orbit.
//!
//! With original port flow y=G^T grad H, solve x'=(J-R)grad H+G(f-Z*y).
//! Augment the HB unknown by y: P(s)=[[sI,G Z(s)],[0,I]] and the nonlinear
//! island is [(J-R)grad H+Gf, G^T grad H]. No storage or memory coordinate is
//! condensed. The load's exp(-i omega t) convention is conjugated exactly once
//! for fs-orbit's positive-frequency Fourier coefficients. A table supplies
//! only the declared harmonics; it is neither a causal fit nor a time ODE.

use super::{
    PeriodicPortInput, PhsOrbitError, PhsOrbitProblem, PhsStorageHessianAction,
    audit_orbit_with_input_effort, forced_problem, tangent::StorageTangent,
};
use core::cell::{Cell, RefCell};
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd, det};
use fs_modal::harmonic::{HarmonicResponseError, admit_passive_port_impedance};
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbRefinementCheck, HbScales, OrbitProblem,
    solve_hb_adaptive_scaled_seeded, solve_hb_scaled_seeded,
};
use fs_phs::PortHamiltonian;

pub mod sideband;

/// Bounds for the table, prepared physical maps and independent work audit.
#[derive(Clone, Copy, Debug)]
pub struct PortImpedanceLimits {
    /// Original pHS states, including every material-memory coordinate.
    pub maximum_states: usize,
    /// Complete original physical ports; no subset is silently discarded.
    pub maximum_ports: usize,
    /// Largest supplied harmonic. Additional hard caps are 256 states, ports
    /// and positive harmonics; the real DC matrix is supplied separately.
    pub maximum_harmonics: usize,
    /// Largest augmented real Fourier system, capped additionally at 4096.
    /// HB's Newton, AFT and refinement budgets remain independently applicable.
    pub maximum_hb_unknowns: usize,
    /// Preparation score L*(8*p^3+n*p^2+(n+p)^2), or audit sample products.
    /// The shared fs-la PSD owner additionally caps each dense admission.
    pub maximum_products: u64,
    /// Additional numeric table/map/tangent/audit payload, excluding HB solver
    /// workspace, source storage, allocator metadata and the caller's seed.
    /// A 256 MiB hard cap also applies.
    pub maximum_workspace_bytes: usize,
}

/// A refusal publishes neither a partial orbit nor an altered physical source.
#[derive(Clone, Debug, PartialEq)]
pub enum PortImpedanceError {
    /// Inconsistent clocks, dimensions, coverage or nonfinite physical data.
    Invalid(&'static str),
    /// A complete operation exceeds an explicit resource cap.
    Budget(&'static str),
    /// Cancellation at a bounded source/map/audit boundary.
    Cancelled,
    /// The shared physical impedance admission refused this harmonic (zero is DC).
    Admission {
        /// Supplied integer harmonic; zero names the real DC matrix.
        harmonic: usize,
        /// Unaltered refusal from the shared physical impedance owner.
        source: HarmonicResponseError,
    },
    /// Pointwise physical admission refused an explicitly supplied probe frequency.
    FrequencyAdmission {
        angular_frequency: f64,
        source: HarmonicResponseError,
    },
    /// Original generic sideband or independently refined physical observation refusal.
    Sideband(super::sideband::PhsPortSidebandError),
    /// Retained auxiliary flow disagrees with the original physical observation.
    FlowMismatch { derivative: bool, norms: Vec<f64> },
    /// Original source, tangent or shared HB/AFT/refinement refusal.
    Regime(PhsOrbitError),
}
impl core::fmt::Display for PortImpedanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear physical impedance regime: {self:?}")
    }
}
impl std::error::Error for PortImpedanceError {}
impl From<PhsOrbitError> for PortImpedanceError {
    fn from(value: PhsOrbitError) -> Self {
        Self::Regime(value)
    }
}
impl From<super::sideband::PhsPortSidebandError> for PortImpedanceError {
    fn from(value: super::sideband::PhsPortSidebandError) -> Self {
        if value == super::sideband::PhsPortSidebandError::Cancelled {
            Self::Cancelled
        } else {
            Self::Sideband(value)
        }
    }
}
impl From<fs_orbit::sideband::SidebandError> for PortImpedanceError {
    fn from(value: fs_orbit::sideband::SidebandError) -> Self {
        super::sideband::PhsPortSidebandError::from(value).into()
    }
}
type Result<T> = core::result::Result<T, PortImpedanceError>;
fn invalid(message: &'static str) -> PortImpedanceError {
    PortImpedanceError::Invalid(message)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| PortImpedanceError::Cancelled)
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn resources(limits: PortImpedanceLimits, work: u128, bytes: u128) -> Result<()> {
    if work > u128::from(limits.maximum_products) {
        return Err(PortImpedanceError::Budget(
            "physical impedance work exceeds its limit",
        ));
    }
    if bytes > limits.maximum_workspace_bytes as u128 || bytes > 256 * 1024 * 1024 {
        return Err(PortImpedanceError::Budget(
            "physical impedance workspace exceeds its limit",
        ));
    }
    Ok(())
}

/// Complete physical impedance at DC and every integer multiple of one clock.
///
/// Matrices are row-major, effort=Z*flow, in the original power-conjugate port
/// units. Positive frequencies use exp(-i omega t) peak phasors; impedance is
/// unchanged by the peak/Fourier factor of two. DC must be finite and real.
/// Full mutual and nonreciprocal passive matrices are retained. Admission is
/// pointwise only: the caller owns the physical source, validity band, causality
/// and any distinct heat/radiation interpretation of its absorbed work.
#[derive(Clone, Debug)]
pub struct PeriodicPortImpedance {
    omega: f64,
    ports: usize,
    matrices: Vec<Vec<C64>>,
    limits: PortImpedanceLimits,
}
impl PeriodicPortImpedance {
    /// Admit an explicitly evaluated DC..N table before any storage evaluation.
    /// There is no interpolation, missing-harmonic default, DC limit inference,
    /// reciprocity assumption, pole fit, or artificial resistance.
    ///
    /// # Errors
    /// Invalid/complex DC, incomplete/nonfinite rows, nonpassive directions,
    /// count/work/payload limits or cancellation. Singular DC loads must use an
    /// explicit physical realization instead of an invented finite matrix.
    pub fn new(
        omega_rad_s: f64,
        ports: usize,
        matrices: Vec<Vec<C64>>,
        limits: PortImpedanceLimits,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        if ports == 0
            || ports > 256
            || ports > limits.maximum_ports
            || matrices.is_empty()
            || matrices.len() > 257
            || matrices.len() - 1 > limits.maximum_harmonics
        {
            return Err(invalid(
                "complete bounded physical ports and DC..N table required",
            ));
        }
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || !(core::f64::consts::TAU / omega_rad_s).is_finite()
            || !(omega_rad_s * (matrices.len() - 1) as f64).is_finite()
            || matrices
                .iter()
                .any(|row| row.len() != ports * ports || row.iter().any(|z| !finite(*z)))
            || matrices[0].iter().any(|z| z.im != 0.0)
        {
            return Err(invalid(
                "finite physical harmonic matrices, real DC and representable clock required",
            ));
        }
        let (p, l) = (ports as u128, matrices.len() as u128);
        resources(
            limits,
            8 * l * p * p * p,
            16 * l * p * p + 256 * p * p + 128 * p,
        )?;
        for (harmonic, matrix) in matrices.iter().enumerate() {
            poll(cx)?;
            admit_passive_port_impedance(ports, matrix, cx).map_err(|source| match source {
                HarmonicResponseError::Cancelled => PortImpedanceError::Cancelled,
                source => PortImpedanceError::Admission { harmonic, source },
            })?;
        }
        Ok(Self {
            omega: omega_rad_s,
            ports,
            matrices,
            limits,
        })
    }
    /// Largest explicitly supplied harmonic, excluding DC.
    #[must_use]
    pub fn harmonic_count(&self) -> usize {
        self.matrices.len() - 1
    }
    /// Fundamental physical angular frequency [rad/s].
    #[must_use]
    pub fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    /// Borrow an admitted original-convention matrix; absent bands stay absent.
    #[must_use]
    pub fn matrix(&self, harmonic: usize) -> Option<&[C64]> {
        self.matrices.get(harmonic).map(Vec::as_slice)
    }

    /// Prepare full original-state back-reaction and cache each harmonic port.
    /// The drive, load and original physical port ordering must agree exactly.
    /// Preparation evaluates neither the source storage nor any load callback.
    pub fn prepare<'a, 's>(
        &'a self,
        system: &'a PortHamiltonian,
        input: &'a PeriodicPortInput,
        cx: &'a Cx<'s>,
    ) -> Result<PreparedPhsImpedance<'a, 's>> {
        poll(cx)?;
        let (n, p) = (system.state_dim(), self.ports);
        if n == 0
            || n > 256
            || n > self.limits.maximum_states
            || system.port_dim() != p
            || input.port_count() != p
            || input.omega_rad_s().to_bits() != self.omega.to_bits()
        {
            return Err(invalid(
                "load and drive must match every original state/port and clock",
            ));
        }
        let d = n + p;
        let (nu, pu, du, l) = (n as u128, p as u128, d as u128, self.matrices.len() as u128);
        resources(
            self.limits,
            l * (8 * pu * pu * pu + nu * pu * pu + du * du),
            16 * l * (du * du + pu * pu) + 80 * du * du + 128 * du,
        )?;
        let (_, _, g) = system.structure();
        if g.len() != n * p || g.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "finite complete original physical port map required",
            ));
        }
        let mut operators = Vec::with_capacity(self.matrices.len());
        for (h, impedance) in self.matrices.iter().enumerate() {
            poll(cx)?;
            operators.push(physical_operator(
                system,
                h as f64 * self.omega,
                impedance,
                false,
                cx,
            )?);
        }
        Ok(PreparedPhsImpedance {
            original: forced_problem(system, input)?,
            load: self,
            input,
            cx,
            operators,
            cancelled: Cell::new(false),
            port_refused: Cell::new(false),
            tangent_cache: RefCell::new(None),
            tangent_scratch: RefCell::new((vec![0.0; n], vec![0.0; n], vec![0.0; d])),
            zero_input: vec![0.0; p],
        })
    }
}

// One complete original physical map serves base harmonics, signed probe
// sidebands and their declared frequency derivatives. Inputs use exp(-iwt).
fn physical_operator(
    system: &PortHamiltonian,
    frequency: f64,
    impedance: &[C64],
    derivative: bool,
    cx: &Cx<'_>,
) -> Result<Vec<C64>> {
    let (n, p) = (system.state_dim(), system.port_dim());
    let d = n + p;
    let (_, _, g) = system.structure();
    let mut operator = vec![C64::ZERO; d * d];
    for i in 0..n {
        operator[i * d + i] = C64::new(0.0, if derivative { 1.0 } else { frequency });
    }
    if !derivative {
        for a in 0..p {
            operator[(n + a) * d + n + a] = C64::ONE;
        }
    }
    for i in 0..n {
        poll(cx)?;
        for b in 0..p {
            let mut sum = ComplexSum::default();
            for a in 0..p {
                sum.product(C64::from_re(g[i * p + a]), impedance[a * p + b].conj());
            }
            operator[i * d + n + b] = sum.value()?;
        }
    }
    Ok(operator)
}

/// Numerical evidence remains that of the original shared HB implementation.
#[derive(Clone, Debug)]
pub enum PhsImpedanceEvidence {
    /// Fixed-order physical scaled residual and independent quadrature change.
    Fixed {
        /// Full augmented equation norm in the authored physical row scales.
        scaled_residual: f64,
        /// Independent AFT grid change in the same normalized chart.
        scaled_sampling_defect: f64,
    },
    /// Successive solved orders, omitted-equation checks and waveform changes.
    Adaptive(Vec<HbRefinementCheck>),
}
/// Loaded periodic response with the original pHS state order intact.
#[derive(Clone, Debug)]
pub struct PhsImpedanceOrbit {
    /// Only original states, including memory. Residual diagnostics still measure
    /// the COMPLETE augmented balance, including the physical flow constraints.
    pub state: HbOrbit,
    /// Retained physical flow Fourier coefficients, exp(+i n omega t).
    /// DC is real; positive-frequency peak exp(-i omega t) phasors are 2*conj(Y).
    pub flow_coefficients: Vec<Vec<C64>>,
    /// Opposing boundary effort coefficients in the same Fourier chart.
    /// Applied net input is authored drive minus this effort, not plus it.
    pub effort_coefficients: Vec<Vec<C64>>,
    /// Fixed-order or adaptive evidence for the full augmented problem.
    pub evidence: PhsImpedanceEvidence,
}

/// Prepared forced physical impedance problem; no time-domain RHS is supplied.
///
/// Packed seeds/scales have original x coordinates followed by physical y.
/// State-equation scales have x/time units, whereas added equations have flow
/// units. Use HbOrbit::packed_seed on that augmented coordinate ordering.
pub struct PreparedPhsImpedance<'a, 's> {
    original: PhsOrbitProblem<'a>,
    load: &'a PeriodicPortImpedance,
    input: &'a PeriodicPortInput,
    cx: &'a Cx<'s>,
    operators: Vec<Vec<C64>>,
    cancelled: Cell<bool>,
    port_refused: Cell<bool>,
    tangent_cache: RefCell<Option<Vec<f64>>>,
    tangent_scratch: RefCell<(Vec<f64>, Vec<f64>, Vec<f64>)>,
    zero_input: Vec<f64>,
}
impl<'a, 's> PreparedPhsImpedance<'a, 's> {
    /// Override with the COMPLETE original continuous storage Hessian action.
    /// Native exact quadratic/continuous tangents are discovered automatically.
    /// A declared action's refusal never falls back to finite differences.
    #[must_use]
    pub fn with_storage_hessian(mut self, hessian: &'a PhsStorageHessianAction<'a>) -> Self {
        self.original.tangent = StorageTangent::hessian(self.original.system.state_dim(), hessian);
        self.tangent_cache = RefCell::new(None);
        self
    }
    fn admit_solve(&self, budget: &HbBudget, maximum: usize) -> Result<()> {
        poll(self.cx)?;
        if self.input.highest_nonzero_harmonic() > budget.harmonics
            || maximum > self.load.harmonic_count()
            || maximum < budget.harmonics
        {
            return Err(invalid(
                "drive must fit initial order and load must cover every retained/refinement harmonic",
            ));
        }
        let unknowns = (self.dim() as u128) * (2 * maximum as u128 + 1);
        if unknowns > self.load.limits.maximum_hb_unknowns as u128 || unknowns > 4096 {
            return Err(PortImpedanceError::Budget(
                "augmented HB Fourier system exceeds its limit",
            ));
        }
        self.cancelled.set(false);
        self.port_refused.set(false);
        self.original.source_refused.set(false);
        Ok(())
    }
    fn after_solve(&self) -> Result<()> {
        if self.cancelled.get() {
            return Err(PortImpedanceError::Cancelled);
        }
        poll(self.cx)?;
        if self.port_refused.get() {
            return Err(invalid("HB queried an unprovided physical load frequency"));
        }
        if self.original.source_refused() {
            return Err(PhsOrbitError::SourceRefused.into());
        }
        Ok(())
    }
    /// Fixed-order loaded HB using the existing scaled Newton/AFT solver.
    /// Nonlinear output harmonics require separate refinement; no infinite tail
    /// or stability claim follows from convergence at one retained order.
    pub fn solve_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        scales: &HbScales,
    ) -> Result<PhsImpedanceOrbit> {
        self.admit_solve(budget, budget.harmonics)?;
        let result = solve_hb_scaled_seeded(
            self,
            HbAnchor::Forced {
                omega: self.load.omega,
            },
            seed,
            budget,
            scales,
        );
        self.after_solve()?;
        let result = result.map_err(PhsOrbitError::Orbit)?;
        self.finish(
            result.orbit,
            PhsImpedanceEvidence::Fixed {
                scaled_residual: result.scaled_residual,
                scaled_sampling_defect: result.scaled_sampling_defect,
            },
        )
    }
    /// Refine the same full loaded problem with original fs-orbit diagnostics.
    /// The table must cover DC..max_harmonics before any source evaluation.
    /// The 2N+1 omitted-equation assessment needs no unprovided load values:
    /// omitted state and flow coefficients are exactly zero, as in fs-orbit.
    pub fn solve_adaptive_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        refinement: &HbRefinementBudget,
        scales: &HbScales,
    ) -> Result<PhsImpedanceOrbit> {
        self.admit_solve(budget, refinement.max_harmonics)?;
        let result = solve_hb_adaptive_scaled_seeded(
            self,
            HbAnchor::Forced {
                omega: self.load.omega,
            },
            seed,
            budget,
            refinement,
            scales,
        );
        self.after_solve()?;
        let result = result.map_err(PhsOrbitError::Orbit)?;
        self.finish(
            result.orbit,
            PhsImpedanceEvidence::Adaptive(result.refinement_trace),
        )
    }
    fn finish(
        &self,
        mut state: HbOrbit,
        evidence: PhsImpedanceEvidence,
    ) -> Result<PhsImpedanceOrbit> {
        let flow_coefficients = state.coeffs.split_off(self.original.system.state_dim());
        let p = self.load.ports;
        let width = flow_coefficients[0].len();
        let mut effort_coefficients = vec![vec![C64::ZERO; width]; p];
        for h in 0..width {
            poll(self.cx)?;
            for a in 0..p {
                let mut sum = ComplexSum::default();
                for b in 0..p {
                    sum.product(
                        self.load.matrices[h][a * p + b].conj(),
                        flow_coefficients[b][h],
                    );
                }
                effort_coefficients[a][h] = sum.value()?;
            }
        }
        Ok(PhsImpedanceOrbit {
            state,
            flow_coefficients,
            effort_coefficients,
            evidence,
        })
    }
    fn source_poll(&self) -> bool {
        if self.cx.checkpoint().is_err() {
            self.cancelled.set(true);
            false
        } else {
            true
        }
    }
    fn state_jacobian(&self, _time: f64, state: &[f64], out: &mut [f64]) -> bool {
        let n = self.original.system.state_dim();
        let d = self.dim();
        if !self.source_poll()
            || state.len() != d
            || out.len() != d * d
            || state.iter().any(|v| !v.is_finite())
        {
            return false;
        }
        if let Some(cached) = self.tangent_cache.borrow().as_ref() {
            out.copy_from_slice(cached);
            return true;
        }
        out.fill(0.0);
        let mut scratch = self.tangent_scratch.borrow_mut();
        let (direction, effort, flow) = &mut *scratch;
        for column in 0..n {
            if !self.source_poll() {
                return false;
            }
            direction.fill(0.0);
            direction[column] = 1.0;
            if !self
                .original
                .tangent
                .apply_hessian(&state[..n], direction, effort)
                || self
                    .original
                    .system
                    .vector_field_from_effort_into(effort, &self.zero_input, &mut flow[..n])
                    .is_err()
                || self
                    .original
                    .system
                    .output_from_effort_into(effort, &mut flow[n..])
                    .is_err()
            {
                return false;
            }
            for row in 0..d {
                out[row * d + column] = flow[row];
            }
        }
        if matches!(self.original.tangent, StorageTangent::Quadratic { .. }) {
            *self.tangent_cache.borrow_mut() = Some(out.to_vec());
        }
        true
    }
}
impl OrbitProblem for PreparedPhsImpedance<'_, '_> {
    fn dim(&self) -> usize {
        self.original.system.state_dim() + self.load.ports
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island_state_jacobian(&self) -> Option<fs_orbit::IslandStateJacobian<Self>> {
        self.original
            .tangent
            .available()
            .then_some(Self::state_jacobian)
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        let n = self.original.system.state_dim();
        if !self.source_poll()
            || state.len() != self.dim()
            || out.len() != self.dim()
            || state.iter().any(|v| !v.is_finite())
        {
            self.original.fail(out);
            return;
        }
        self.original.island(time, &state[..n], &mut out[..n]);
        if self.original.source_refused()
            || self
                .original
                .system
                .output_from_effort_into(&self.original.effort.borrow(), &mut out[n..])
                .is_err()
        {
            self.original.fail(out);
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let harmonic = (s.im / self.load.omega).round();
        if !self.source_poll()
            || !finite(s)
            || s.re != 0.0
            || harmonic < 0.0
            || harmonic > self.load.harmonic_count() as f64
            || (harmonic * self.load.omega).to_bits() != s.im.to_bits()
        {
            self.port_refused.set(true);
            return vec![C64::new(f64::NAN, f64::NAN); self.dim() * self.dim()];
        }
        self.operators[harmonic as usize].clone()
    }
    // Inherit time_domain_rhs=None: a frequency table is not a state realization.
}

#[derive(Default)]
struct ComplexSum {
    re: Dd,
    im: Dd,
}
impl ComplexSum {
    fn product(&mut self, a: C64, b: C64) {
        let (ar, ai, br, bi) = (
            Dd::from_f64(a.re),
            Dd::from_f64(a.im),
            Dd::from_f64(b.re),
            Dd::from_f64(b.im),
        );
        self.re = self.re + ar * br - ai * bi;
        self.im = self.im + ar * bi + ai * br;
    }
    fn value(self) -> Result<C64> {
        let value = C64::new(self.re.to_f64(), self.im.to_f64());
        if finite(value) {
            Ok(value)
        } else {
            Err(invalid("unrepresentable physical impedance contraction"))
        }
    }
}

/// Independent original-storage work and physical-output checks over one cycle.
/// A single sample grid is measured evidence; nonlinear power/flow needs a
/// successive-grid comparison, even after harmonic-order refinement.
#[derive(Clone, Debug)]
pub struct PhsImpedanceAudit {
    /// Number of uniformly spaced physical cycle samples.
    pub samples: usize,
    /// Original cycle-mean stored energy, including all material memory [J].
    pub mean_stored_energy_j: f64,
    /// Authored external source work rate through original e^T G f [W].
    pub mean_source_power_w: f64,
    /// Original e^T R e [W], kept distinct from absorbed boundary work.
    /// Its physical owner determines which part represents structural heat.
    pub mean_native_resistance_loss_w: f64,
    /// Original e^T G*(Z*y) [W], sampled using the original nonlinear effort.
    /// No split into radiation, boundary heat or stored load energy is inferred.
    pub mean_boundary_absorbed_power_w: f64,
    /// DC plus Parseval harmonic work using the retained load-flow unknowns [W].
    /// Its difference from the preceding field exposes physical-flow mismatch.
    pub retained_boundary_absorbed_power_w: f64,
    /// External source minus native resistance minus original boundary work [W].
    pub cycle_power_residual_w: f64,
    /// RMS of original G^T grad H(x(t)) minus retained flow, per physical port.
    /// Includes nonlinear output above the retained harmonic order.
    pub physical_flow_rms_defect: Vec<f64>,
}
impl PreparedPhsImpedance<'_, '_> {
    /// Audit the full physical source without treating the augmented flow as the
    /// original nonlinear output. The grid must resolve retained state/drive
    /// support; nonlinear work still requires successive-grid convergence.
    #[allow(clippy::too_many_lines)] // one transactional original-source cycle audit
    pub fn audit(&self, result: &PhsImpedanceOrbit, samples: usize) -> Result<PhsImpedanceAudit> {
        poll(self.cx)?;
        let p = self.load.ports;
        let width = result.state.coeffs.first().map_or(0, Vec::len);
        if result.state.omega.to_bits() != self.load.omega.to_bits()
            || result.flow_coefficients.len() != p
            || result.effort_coefficients.len() != p
            || width == 0
            || width > self.load.matrices.len()
            || result
                .flow_coefficients
                .iter()
                .chain(&result.effort_coefficients)
                .any(|row| row.len() != width)
            || !(32..=65_536).contains(&samples)
            || self.input.highest_nonzero_harmonic() > (samples - 1) / 2
        {
            return Err(invalid(
                "audit needs matching complete physical spectra and a resolving bounded grid",
            ));
        }
        let d = self.dim() as u128;
        resources(
            self.load.limits,
            samples as u128 * (d * d + 4 * p as u128 * width as u128)
                + 4 * width as u128 * p as u128 * p as u128,
            64 * p as u128 * width as u128 + 128 * d,
        )?;
        let flows = PeriodicPortInput::new(self.load.omega, result.flow_coefficients.clone())?;
        let efforts = PeriodicPortInput::new(self.load.omega, result.effort_coefficients.clone())?;
        // Public results may have been edited or come from another load. Do not
        // audit their effort as this source's boundary unless Z*y still agrees.
        for h in 0..width {
            poll(self.cx)?;
            for a in 0..p {
                let mut sum = ComplexSum::default();
                for b in 0..p {
                    sum.product(
                        self.load.matrices[h][a * p + b].conj(),
                        result.flow_coefficients[b][h],
                    );
                }
                if sum.value()? != result.effort_coefficients[a][h] {
                    return Err(invalid(
                        "audited boundary effort differs from the prepared physical load",
                    ));
                }
            }
        }
        let mut flow = vec![0.0; p];
        let mut effort = vec![0.0; p];
        let mut defects = vec![0.0; p];
        let mut boundary = Dd::ZERO;
        let mut cancelled = false;
        let base = audit_orbit_with_input_effort(
            self.original.system,
            &result.state,
            samples,
            |time, _, original_effort, original_flow, out| {
                if self.cx.checkpoint().is_err() {
                    cancelled = true;
                    return false;
                }
                if !flows.sample_into(time, &mut flow)
                    || !efforts.sample_into(time, &mut effort)
                    || !self.input.sample_into(time, out)
                {
                    return false;
                }
                let Ok(work) = self
                    .original
                    .system
                    .supplied_power_from_effort(original_effort, &effort)
                else {
                    return false;
                };
                boundary = boundary + Dd::from_f64(work);
                for a in 0..p {
                    defects[a] = det::hypot(defects[a], original_flow[a] - flow[a]);
                }
                true
            },
        );
        if cancelled {
            return Err(PortImpedanceError::Cancelled);
        }
        let base = base?;
        let boundary = (boundary / Dd::from_f64(samples as f64)).to_f64();
        let mut retained = Dd::ZERO;
        for a in 0..p {
            poll(self.cx)?;
            for h in 0..width {
                let y = result.flow_coefficients[a][h];
                let e = result.effort_coefficients[a][h];
                let work = Dd::from_f64(y.re) * Dd::from_f64(e.re)
                    + Dd::from_f64(y.im) * Dd::from_f64(e.im);
                retained = retained + work * Dd::from_f64(if h == 0 { 1.0 } else { 2.0 });
            }
            defects[a] /= (samples as f64).sqrt();
        }
        let residual = (Dd::from_f64(base.mean_supplied_power_w)
            - Dd::from_f64(base.mean_resistance_loss_w)
            - Dd::from_f64(boundary))
        .to_f64();
        if !boundary.is_finite()
            || !retained.to_f64().is_finite()
            || !residual.is_finite()
            || defects.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable original physical work audit"));
        }
        Ok(PhsImpedanceAudit {
            samples,
            mean_stored_energy_j: base.mean_stored_energy_j,
            mean_source_power_w: base.mean_supplied_power_w,
            mean_native_resistance_loss_w: base.mean_resistance_loss_w,
            mean_boundary_absorbed_power_w: boundary,
            retained_boundary_absorbed_power_w: retained.to_f64(),
            cycle_power_residual_w: residual,
            physical_flow_rms_defect: defects,
        })
    }
}
