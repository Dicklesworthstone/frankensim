//! Periodically modulated quadratic energy on an original nonlinear pHS.
//!
//! H(t,x) = H0(x) + x^T Q(t) x/2, with unchanged original J, R and G.
//! Q(t) is an authored physical energy addition, not a fitted modal frequency.
//! Its explicit pump power is x^T Qdot(t) x/2. Symmetry does not establish
//! positivity, passivity, bounded motion or Floquet stability.

use super::{PeriodicPortInput, PhsStorageHessianAction, tangent::StorageTangent};
use core::cell::{Cell, RefCell};
use fs_math::{c64::C64, dd::Dd, det};
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbScales, OrbitError,
    OrbitProblem, ScaledHbOrbit, ShootBudget, VariationalShootControls, VariationalShootError,
    VariationalShootOrbit,
    sideband::{
        PreparedPeriodicLinearization, SidebandControls, SidebandError, SidebandResponse,
        SidebandScales,
    },
    solve_hb_adaptive_scaled_seeded, solve_hb_scaled_seeded, solve_shooting_variational,
};
use fs_phs::PortHamiltonian;
pub use fs_phs::harmonic::equilibrium::NonlinearDissipation;

mod dissipation;
mod periods;

mod sensitivity;
pub use sensitivity::{
    ParametricParameterDirection, ParametricPowerObjective, PreparedParametricPowerAdjoint,
    PreparedParametricSensitivity,
};

/// Hard source limits, independent of every solver's own numerical budgets.
pub const MAX_PARAMETRIC_STATES: usize = 256;
/// Maximum positive harmonic index of the authored pump.
pub const MAX_PARAMETRIC_HARMONICS: usize = 256;
/// Maximum supplied complex entries, including all full symmetric matrices.
pub const MAX_PARAMETRIC_COEFFICIENTS: usize = 1_048_576;

/// Named refusals; no partial numerical response is published.
#[derive(Clone, Debug, PartialEq)]
pub enum ParametricError {
    Invalid(&'static str),
    SourceRefused,
    Orbit(OrbitError),
    Shooting(VariationalShootError),
    Sideband(SidebandError),
    Physical(fs_phs::PhsError),
    Dissipation(fs_phs::harmonic::equilibrium::LinearizationError),
    Sensitivity(fs_orbit::sensitivity::ForcedSensitivityError),
}
impl core::fmt::Display for ParametricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pHS parametric excitation refused: {self:?}")
    }
}
impl std::error::Error for ParametricError {}
type Result<T> = core::result::Result<T, ParametricError>;

struct Entry {
    row: usize,
    column: usize,
    coefficients: Vec<C64>,
}

/// Real symmetric Fourier matrix in the ORIGINAL complete state coordinates.
///
/// `matrices[n][row*d+column]` is the positive-frequency Fourier coefficient
/// Q_n in Q(t)=Q_0+2 Re sum Q_n exp(+i*n*omega*t). Q_0 is real. Each complex
/// matrix is TRANSPOSE-symmetric, not Hermitian: both real and imaginary parts
/// are symmetric. A cosine pump of peak matrix amplitude A therefore uses A/2.
/// Coefficients have energy/(state_row*state_column) units; state coordinates
/// can carry different units. No mass normalization or parameter law is inferred.
///
/// Exactly absent entries/zero padding are skipped; every nonzero coefficient
/// and its phase is retained. Signed energy additions are permitted, including
/// stiffness modulation through zero. The author owns physical applicability
/// and stability/energy support; this source is not a passivity certificate.
pub struct PeriodicQuadraticModulation {
    dimension: usize,
    omega: f64,
    harmonics: usize,
    entries: Vec<Entry>,
}
impl PeriodicQuadraticModulation {
    pub fn new(omega: f64, dimension: usize, matrices: Vec<Vec<C64>>) -> Result<Self> {
        let square = dimension
            .checked_mul(dimension)
            .ok_or(ParametricError::Invalid(
                "modulation matrix extent overflows",
            ))?;
        let count = square
            .checked_mul(matrices.len())
            .ok_or(ParametricError::Invalid(
                "modulation coefficient extent overflows",
            ))?;
        if dimension == 0
            || dimension > MAX_PARAMETRIC_STATES
            || matrices.is_empty()
            || matrices.len() > MAX_PARAMETRIC_HARMONICS + 1
            || count > MAX_PARAMETRIC_COEFFICIENTS
            || !positive(omega)
            || !positive(core::f64::consts::TAU / omega)
            || matrices
                .iter()
                .any(|matrix| matrix.len() != square || matrix.iter().any(|z| !finite(*z)))
            || matrices[0].iter().any(|z| z.im != 0.0)
        {
            return Err(ParametricError::Invalid(
                "finite complete bounded Fourier matrices and real DC required",
            ));
        }
        for matrix in &matrices {
            for row in 0..dimension {
                for column in 0..row {
                    if matrix[row * dimension + column] != matrix[column * dimension + row] {
                        return Err(ParametricError::Invalid(
                            "every Fourier energy matrix must be exactly transpose-symmetric",
                        ));
                    }
                }
            }
        }
        let mut harmonics = 0;
        let mut entries = Vec::new();
        for row in 0..dimension {
            for column in row..dimension {
                let index = row * dimension + column;
                let Some(last) = matrices
                    .iter()
                    .rposition(|matrix| matrix[index] != C64::ZERO)
                else {
                    continue;
                };
                harmonics = harmonics.max(last);
                entries.push(Entry {
                    row,
                    column,
                    coefficients: matrices[..=last]
                        .iter()
                        .map(|matrix| matrix[index])
                        .collect(),
                });
            }
        }
        if harmonics > 0 && !positive(harmonics as f64 * omega) {
            return Err(ParametricError::Invalid(
                "pump harmonic frequencies must be representable",
            ));
        }
        Ok(Self {
            dimension,
            omega,
            harmonics,
            entries,
        })
    }
    pub fn state_dimension(&self) -> usize {
        self.dimension
    }
    pub fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn highest_nonzero_harmonic(&self) -> usize {
        self.harmonics
    }

    fn sample(&self, time: f64, scratch: &mut Scratch, rate: bool) -> Result<()> {
        if !time.is_finite() {
            return Err(ParametricError::Invalid("finite physical time required"));
        }
        if scratch.time == Some(time.to_bits()) && (!rate || scratch.has_rate) {
            return Ok(());
        }
        scratch.time = None;
        scratch.q.fill(0.0);
        if rate {
            scratch.q_rate.fill(0.0);
        }
        for n in 1..=self.harmonics {
            let phase = n as f64 * self.omega * time;
            if !phase.is_finite() {
                return Err(ParametricError::SourceRefused);
            }
            let (sin, cos) = det::sin_cos(phase);
            scratch.phases[n] = C64::new(cos, sin);
        }
        for entry in &self.entries {
            let mut value = Dd::from_f64(entry.coefficients[0].re);
            let mut derivative = Dd::ZERO;
            for (n, coefficient) in entry.coefficients.iter().enumerate().skip(1) {
                let phase = scratch.phases[n];
                let re = Dd::from_f64(coefficient.re);
                let im = Dd::from_f64(coefficient.im);
                value = value
                    + (re * Dd::from_f64(phase.re) - im * Dd::from_f64(phase.im))
                        * Dd::from_f64(2.0);
                if rate {
                    derivative = derivative
                        - (re * Dd::from_f64(phase.im) + im * Dd::from_f64(phase.re))
                            * Dd::from_f64(n as f64 * self.omega)
                            * Dd::from_f64(2.0);
                }
            }
            let (value, derivative) = (value.to_f64(), derivative.to_f64());
            if !value.is_finite() || (rate && !derivative.is_finite()) {
                return Err(ParametricError::SourceRefused);
            }
            for index in [
                entry.row * self.dimension + entry.column,
                entry.column * self.dimension + entry.row,
            ] {
                scratch.q[index] = value;
                if rate {
                    scratch.q_rate[index] = derivative;
                }
            }
        }
        scratch.time = Some(time.to_bits());
        scratch.has_rate = rate;
        Ok(())
    }
}

struct Scratch {
    q: Vec<f64>,
    q_rate: Vec<f64>,
    phases: Vec<C64>,
    time: Option<u64>,
    has_rate: bool,
    effort: Vec<f64>,
    ports: Vec<f64>,
    direction: Vec<f64>,
    flow: Vec<f64>,
    dissipation: Option<dissipation::Scratch>,
}

impl Scratch {
    fn new(n: usize, m: usize, harmonics: usize) -> Self {
        Self {
            q: vec![0.0; n * n],
            q_rate: vec![0.0; n * n],
            phases: vec![C64::ZERO; harmonics + 1],
            time: None,
            has_rate: false,
            effort: vec![0.0; n],
            ports: vec![0.0; m],
            direction: vec![0.0; n],
            flow: vec![0.0; n],
            dissipation: None,
        }
    }
}

/// Instantaneous power identity of the original explicitly driven Hamiltonian.
/// Positive pump/port power enters the system; resistance power leaves it.
#[derive(Clone, Copy, Debug)]
pub struct ParametricPower {
    pub energy: f64,
    /// Explicit partial derivative of H at fixed state: x^T Qdot x/2.
    pub pump_power: f64,
    pub port_power: f64,
    pub resistance_power: f64,
    /// Original nonlinear internal loss e^T D(x,e), with the COMPLETE pumped
    /// effort e=H0_x+Q(t)x. Zero when no nonlinear resistance is attached.
    /// Roundoff-sized signed values admitted by the original law are retained.
    pub nonlinear_power: f64,
    /// Total dH/dt = gradient H dot xdot + partial H/partial t.
    pub energy_rate: f64,
}

/// Prepared original source, periodic pump and drive, with reusable scratch.
///
/// The source is smooth in state wherever H0 declares that capability. Every
/// original coordinate, including material-memory states, remains present. The
/// private OrbitProblem adapter prevents unadmitted pump/drive bandwidth from
/// reaching the shared HB/sideband owners through these entry points. No new
/// Newton, Fourier, time integrator, linear solver or eigensolver is supplied.
pub struct ParametricHamiltonian<'a> {
    system: &'a PortHamiltonian,
    input: &'a PeriodicPortInput,
    modulation: &'a PeriodicQuadraticModulation,
    periods: usize,
    tangent: StorageTangent<'a>,
    dissipation: Option<NonlinearDissipation<'a>>,
    scratch: RefCell<Scratch>,
    source_refused: Cell<bool>,
    failure: RefCell<Option<ParametricError>>,
}
impl<'a> ParametricHamiltonian<'a> {
    pub fn new(
        system: &'a PortHamiltonian,
        input: &'a PeriodicPortInput,
        modulation: &'a PeriodicQuadraticModulation,
    ) -> Result<Self> {
        Self::build(system, input, modulation, None)
    }

    /// Supply the complete original H0 Hessian for a storage without a native
    /// action, or an explicitly authored override. Q(t) is added exactly once.
    pub fn new_analytic(
        system: &'a PortHamiltonian,
        input: &'a PeriodicPortInput,
        modulation: &'a PeriodicQuadraticModulation,
        hessian: &'a PhsStorageHessianAction<'a>,
    ) -> Result<Self> {
        Self::build(system, input, modulation, Some(hessian))
    }

    fn build(
        system: &'a PortHamiltonian,
        input: &'a PeriodicPortInput,
        modulation: &'a PeriodicQuadraticModulation,
        hessian: Option<&'a PhsStorageHessianAction<'a>>,
    ) -> Result<Self> {
        let (n, m) = (system.state_dim(), system.port_dim());
        let (j, r, g) = system.structure();
        if n != modulation.dimension
            || m != input.port_count()
            || input.omega_rad_s() != modulation.omega
            || j.len() != n * n
            || r.len() != n * n
            || g.len() != n * m
            || j.iter().chain(r).chain(g).any(|v| !v.is_finite())
        {
            return Err(ParametricError::Invalid(
                "original source, pump and drive must share dimensions and clock",
            ));
        }
        let tangent = hessian.map_or_else(
            || StorageTangent::from_storage(system),
            |action| StorageTangent::hessian(n, action),
        );
        Ok(Self {
            system,
            input,
            modulation,
            periods: 1,
            tangent,
            dissipation: None,
            scratch: RefCell::new(Scratch::new(n, m, modulation.harmonics)),
            source_refused: Cell::new(false),
            failure: RefCell::new(None),
        })
    }

    fn bandwidth(&self) -> usize {
        self.modulation
            .harmonics
            .max(self.input.highest_nonzero_harmonic())
            * self.periods
    }
    fn admit_order(&self, harmonics: usize) -> Result<()> {
        if harmonics < self.bandwidth() {
            Err(ParametricError::Invalid(
                "retained order must resolve every nonzero drive and pump harmonic",
            ))
        } else {
            Ok(())
        }
    }
    fn effort(&self, time: f64, state: &[f64], scratch: &mut Scratch, rate: bool) -> Result<()> {
        if state.len() != self.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(ParametricError::Invalid(
                "complete finite original state required",
            ));
        }
        self.modulation.sample(time, scratch, rate)?;
        scratch.effort.fill(f64::NAN);
        self.system.storage().gradient(state, &mut scratch.effort);
        add_matrix_action(&scratch.q, state, &mut scratch.effort)?;
        if !self.input.sample_into(time, &mut scratch.ports) {
            return Err(ParametricError::SourceRefused);
        }
        Ok(())
    }

    /// Full physical xdot=(J-R)e-D(x,e)+G u(t), e=gradient H0+Q(t)x.
    /// D is zero unless attached through `with_dissipation`.
    /// The original compensated field owner sees the combined effort and drive
    /// in ONE contraction; small motion is not lost between opposing fields.
    pub fn vector_field_into(&self, time: f64, state: &[f64], out: &mut [f64]) -> Result<()> {
        let mut scratch = self.scratch.borrow_mut();
        self.effort(time, state, &mut scratch, false)?;
        self.evaluate_dissipation(state, &mut scratch)?;
        self.field_from_effort(&scratch, out)
    }

    /// Original state-dependent storage, explicit pump work and physical port/
    /// resistance powers at one state. This is an instantaneous identity, not
    /// an integrated work certificate or a passivity/stability determination.
    pub fn power(&self, time: f64, state: &[f64]) -> Result<ParametricPower> {
        let mut scratch = self.scratch.borrow_mut();
        self.effort(time, state, &mut scratch, true)?;
        let nonlinear_power = self.evaluate_dissipation(state, &mut scratch)?;
        let Scratch {
            effort,
            ports,
            flow,
            q,
            q_rate,
            dissipation,
            ..
        } = &mut *scratch;
        match dissipation {
            Some(loss) => self
                .system
                .vector_field_with_dissipation_from_effort_into(effort, ports, &loss.flow, flow),
            None => self
                .system
                .vector_field_from_effort_into(effort, ports, flow),
        }
        .map_err(|_| ParametricError::SourceRefused)?;
        let energy = (Dd::from_f64(self.system.hamiltonian(state))
            + quadratic(q, state) * Dd::from_f64(0.5))
        .to_f64();
        let pump_power = (quadratic(q_rate, state) * Dd::from_f64(0.5)).to_f64();
        let port_power = self
            .system
            .supplied_power_from_effort(effort, ports)
            .map_err(|_| ParametricError::SourceRefused)?;
        let resistance_power = quadratic(self.system.structure().1, effort).to_f64();
        let energy_rate = (dot(effort, flow) + Dd::from_f64(pump_power)).to_f64();
        if [
            energy,
            pump_power,
            port_power,
            resistance_power,
            nonlinear_power,
            energy_rate,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(ParametricError::SourceRefused);
        }
        Ok(ParametricPower {
            energy,
            pump_power,
            port_power,
            resistance_power,
            nonlinear_power,
            energy_rate,
        })
    }

    /// Complete exact continuous tangent; used by the SAME HB, variational
    /// shooting and periodic-sideband owners. No finite-difference replacement
    /// is used when H0 declares a tangent which subsequently refuses.
    pub fn state_jacobian_into(&self, time: f64, state: &[f64], out: &mut [f64]) -> Result<()> {
        let n = self.system.state_dim();
        if state.len() != n || out.len() != n * n || state.iter().any(|v| !v.is_finite()) {
            return Err(ParametricError::Invalid(
                "complete original state and square Jacobian required",
            ));
        }
        if self.dissipation.is_some() {
            return self.dissipative_state_jacobian_into(time, state, out);
        }
        let mut scratch = self.scratch.borrow_mut();
        self.modulation.sample(time, &mut scratch, false)?;
        let Scratch {
            q,
            effort,
            ports,
            direction,
            flow,
            ..
        } = &mut *scratch;
        ports.fill(0.0);
        direction.fill(0.0);
        for column in 0..n {
            direction[column] = 1.0;
            if !self.tangent.apply_hessian(state, direction, effort) {
                return Err(ParametricError::SourceRefused);
            }
            add_matrix_action(q, direction, effort)?;
            self.system
                .vector_field_from_effort_into(effort, ports, flow)
                .map_err(|_| ParametricError::SourceRefused)?;
            for row in 0..n {
                out[row * n + column] = flow[row];
            }
            direction[column] = 0.0;
        }
        Ok(())
    }

    pub fn solve_harmonic_balance_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        scales: &HbScales,
    ) -> Result<ScaledHbOrbit> {
        self.admit_order(budget.harmonics)?;
        self.reset_failure();
        let result = solve_hb_scaled_seeded(
            &Adapter(self),
            HbAnchor::Forced {
                omega: self.angular_frequency(),
            },
            seed,
            budget,
            scales,
        );
        self.finish(result.map_err(ParametricError::Orbit))
    }

    /// The shared adaptive owner increases actual solved orders, checks retained
    /// and omitted original equations and compares successive full waveforms.
    /// The starting seed/order must already represent the entire authored pump
    /// and drive bandwidth; it is never silently truncated or resampled.
    pub fn solve_harmonic_balance_adaptive_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        refinement: &HbRefinementBudget,
        scales: &HbScales,
    ) -> Result<AdaptiveHbOrbit> {
        self.admit_order(budget.harmonics)?;
        self.reset_failure();
        let result = solve_hb_adaptive_scaled_seeded(
            &Adapter(self),
            HbAnchor::Forced {
                omega: self.angular_frequency(),
            },
            seed,
            budget,
            refinement,
            scales,
        );
        self.finish(result.map_err(ParametricError::Orbit))
    }

    /// Independent original-ODE period solve and exact discrete RK4 monodromy.
    /// Resolve pump and drive on the declared time grid, then refine that grid
    /// separately. The returned finite-grid multipliers imply no stability claim.
    /// Shared controls bound the solver; this already prepared source's bounded
    /// coefficient/storage work and retained scratch are excluded from those caps.
    pub fn solve_shooting_variational(
        &self,
        seed: &[f64],
        budget: &ShootBudget,
        controls: &VariationalShootControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VariationalShootOrbit> {
        if cancelled() {
            return Err(ParametricError::Shooting(VariationalShootError::Cancelled));
        }
        if budget.steps_per_period == 0 || self.bandwidth() > (budget.steps_per_period - 1) / 2 {
            return Err(ParametricError::Invalid(
                "shooting grid must resolve every pump and drive harmonic",
            ));
        }
        self.reset_failure();
        let result = solve_shooting_variational(
            &Adapter(self),
            seed,
            core::f64::consts::TAU / self.angular_frequency(),
            budget,
            controls,
            &mut cancelled,
        );
        if matches!(result, Err(VariationalShootError::Cancelled)) {
            return Err(ParametricError::Shooting(VariationalShootError::Cancelled));
        }
        self.finish(result.map_err(ParametricError::Shooting))
    }

    /// Physical-port probe around the complete periodically pumped base orbit.
    /// Index (k+N)*ports+p is a PEAK complex port coefficient at carrier+k*omega.
    /// Signed converted sidebands are retained. The original input map forms
    /// forcing and the shared periodic tangent/factor owns frequency conversion.
    /// Adapter forcing scratch and its explicit mapping score are charged to
    /// the supplied caps before the shared owner runs. Already prepared source
    /// storage and source callbacks' internal work remain excluded, as in that
    /// owner's contract. The mapping score is 64*width*d*max(ports,1).
    pub fn solve_port_sideband_response(
        &self,
        orbit: &HbOrbit,
        carrier: f64,
        forcing: &[C64],
        controls: SidebandControls,
        scales: &SidebandScales,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SidebandResponse> {
        if cancelled() {
            return Err(ParametricError::Sideband(SidebandError::Cancelled));
        }
        self.admit_order(
            orbit
                .coeffs
                .first()
                .map_or(0, |row| row.len().saturating_sub(1)),
        )?;
        self.admit_order(controls.harmonics)?;
        let (d, m) = (self.system.state_dim(), self.system.port_dim());
        let width = controls
            .harmonics
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .ok_or(ParametricError::Invalid("sideband extent overflows"))?;
        if orbit.omega != self.angular_frequency()
            || controls.harmonics > 64
            || width.checked_mul(m) != Some(forcing.len())
            || forcing.iter().any(|z| !finite(*z))
        {
            return Err(ParametricError::Invalid(
                "complete finite probe and matching base clock required",
            ));
        }
        let overhead = 16 * d * width + 8 * (m + d);
        let mapping_work = 64 * width as u64 * d as u64 * (m as u64).max(1);
        if overhead > controls.maximum_payload_bytes {
            return Err(ParametricError::Sideband(SidebandError::Budget {
                quantity: "parametric port forcing payload",
                requested: overhead as u128,
                cap: controls.maximum_payload_bytes as u128,
            }));
        }
        if mapping_work > controls.maximum_work {
            return Err(ParametricError::Sideband(SidebandError::Budget {
                quantity: "parametric port forcing work",
                requested: mapping_work as u128,
                cap: controls.maximum_work as u128,
            }));
        }
        let shared_controls = SidebandControls {
            maximum_payload_bytes: controls.maximum_payload_bytes - overhead,
            maximum_work: controls.maximum_work - mapping_work,
            ..controls
        };
        self.reset_failure();
        let adapter = Adapter(self);
        let prepared = PreparedPeriodicLinearization::new(
            &adapter,
            orbit,
            shared_controls,
            scales,
            &mut cancelled,
        );
        if matches!(prepared, Err(SidebandError::Cancelled)) {
            return Err(ParametricError::Sideband(SidebandError::Cancelled));
        }
        let prepared = self.finish(prepared.map_err(ParametricError::Sideband))?;
        let factor = prepared
            .factor(carrier, &mut cancelled)
            .map_err(ParametricError::Sideband)?;
        let mut physical = vec![C64::ZERO; d * width];
        let mut input = vec![0.0; m];
        let mut field = vec![0.0; d];
        for band in 0..width {
            if cancelled() {
                return Err(ParametricError::Sideband(SidebandError::Cancelled));
            }
            for imaginary in [false, true] {
                for port in 0..m {
                    input[port] = if imaginary {
                        forcing[band * m + port].im
                    } else {
                        forcing[band * m + port].re
                    };
                }
                self.system
                    .input_flow_into(&input, &mut field)
                    .map_err(|_| ParametricError::SourceRefused)?;
                for state in 0..d {
                    if imaginary {
                        physical[band * d + state].im = field[state];
                    } else {
                        physical[band * d + state].re = field[state];
                    }
                }
            }
        }
        factor
            .solve(&physical, cancelled)
            .map_err(ParametricError::Sideband)
    }
}

struct Adapter<'s, 'a>(&'s ParametricHamiltonian<'a>);
impl Adapter<'_, '_> {
    fn record(&self, result: Result<()>) -> bool {
        if let Err(error) = result {
            self.0.source_refused.set(true);
            if self.0.dissipation.is_some() {
                let mut failure = self.0.failure.borrow_mut();
                if failure.is_none() {
                    *failure = Some(error);
                }
            }
            false
        } else {
            true
        }
    }
}
impl OrbitProblem for Adapter<'_, '_> {
    fn dim(&self) -> usize {
        self.0.system.state_dim()
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        if !self.record(self.0.vector_field_into(time, state, out)) {
            out.fill(f64::NAN);
        }
    }
    fn island_state_jacobian(&self) -> Option<fs_orbit::IslandStateJacobian<Self>> {
        self.0
            .tangent
            .available()
            .then_some(|problem, time, state, out| {
                let result = problem.0.state_jacobian_into(time, state, out);
                if problem.0.dissipation.is_some() {
                    problem.record(result)
                } else {
                    // Keep the original no-D native tangent-refusal semantics.
                    result.is_ok()
                }
            })
    }
    fn time_domain_rhs(&self) -> Option<fs_orbit::TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<fs_orbit::TimeDomainStateJacobian<Self>> {
        self.island_state_jacobian()
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn add_matrix_action(matrix: &[f64], state: &[f64], out: &mut [f64]) -> Result<()> {
    let n = state.len();
    for row in 0..n {
        let mut value = Dd::from_f64(out[row]);
        for column in 0..n {
            if matrix[row * n + column] != 0.0 {
                value =
                    value + Dd::from_f64(matrix[row * n + column]) * Dd::from_f64(state[column]);
            }
        }
        out[row] = value.to_f64();
        if !out[row].is_finite() {
            return Err(ParametricError::SourceRefused);
        }
    }
    Ok(())
}
fn dot(a: &[f64], b: &[f64]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (a, b)| {
        sum + Dd::from_f64(*a) * Dd::from_f64(*b)
    })
}
fn quadratic(matrix: &[f64], x: &[f64]) -> Dd {
    let n = x.len();
    let mut result = Dd::ZERO;
    for row in 0..n {
        let mut value = Dd::ZERO;
        for column in 0..n {
            if matrix[row * n + column] != 0.0 {
                value = value + Dd::from_f64(matrix[row * n + column]) * Dd::from_f64(x[column]);
            }
        }
        result = result + Dd::from_f64(x[row]) * value;
    }
    result
}

#[cfg(test)]
mod tests;
