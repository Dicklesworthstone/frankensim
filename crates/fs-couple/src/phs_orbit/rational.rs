//! Causal original-state rational loading of nonlinear pHS periodic responses.
//!
//! For the supplied proper source z'=Az+By, e=Cz+Dy, retain all [x,z] and use
//! y=S^T grad H(x), x'=(J-R)grad H+Gf-S e. The SAME physical S carries flow and
//! reaction; external drive ports G remain independent. fs-orbit owns HB, RK4,
//! its exact discrete tangent and Floquet spectrum. No realization, damping,
//! state reduction, direct-term inverse or additional integrator is inferred.
//!
//! This is exact continuous evolution of the SUPPLIED rational model, including
//! its original DC response. Positive retained harmonics must lie in its source
//! band. That band and the original continuous/sampled passivity reports remain
//! unchanged: physical validity of DC, transient content and unresolved higher
//! harmonics needs separate source evidence. Arbitrary rational states have no
//! inferred stored energy or heat. Exchanged boundary work remains signed.

use super::{
    PeriodicPortInput, PhsOrbitError, PhsOrbitProblem, PhsStorageHessianAction, forced_problem,
    tangent::StorageTangent,
};
use core::cell::{Cell, RefCell};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbScales, OrbitProblem,
    ScaledHbOrbit, ShootBudget, VariationalShootControls, VariationalShootOrbit,
    solve_hb_adaptive_scaled_seeded, solve_hb_scaled_seeded, solve_shooting_variational,
};
use fs_phs::PortHamiltonian;
use fs_vfit::mimo::{
    MimoPassivityReport,
    load::{MimoLoadError, PreparedMimoLoad},
};

mod operator;
pub mod sideband;
use operator::{Boundary, dot};

/// Separate adapter limits. Shared HB/shooting workspace and work caps also apply.
#[derive(Clone, Copy, Debug)]
pub struct PhsMimoLimits {
    /// Complete original host plus source states, additionally capped at 256.
    pub maximum_states: usize,
    /// Original boundary ports, additionally capped at 128.
    pub maximum_ports: usize,
    /// Authored and retained positive harmonics, additionally capped at 256.
    pub maximum_harmonics: usize,
    /// Complete real Fourier unknown count, additionally capped at 4096.
    pub maximum_hb_unknowns: usize,
    /// Adapter contraction score consumed by one solve/audit, not solver FLOPs.
    /// Preparation separately requires host_states*source_states*ports <= cap.
    /// Original constitutive/Hessian callback internals remain caller-owned.
    pub maximum_products: u64,
    /// Matrix, tangent and scratch payload, excluding original owners, caller
    /// seeds and shared numerical solvers. A 256 MiB hard limit also applies.
    pub maximum_workspace_bytes: usize,
}

/// A refused query leaves both original source owners and their history untouched.
#[derive(Debug)]
pub enum PhsMimoError {
    /// Incomplete/nonfinite original data, unsupported derivative or grid.
    Invalid(&'static str),
    /// An explicitly bounded matrix, Fourier system or source contraction.
    Budget(&'static str),
    /// Execution scope cancelled before publishing the candidate.
    Cancelled,
    /// Unaltered original realization, history or source-band refusal.
    Source(MimoLoadError),
    /// Original storage/tangent or shared numerical solver refusal.
    Regime(PhsOrbitError),
}
impl core::fmt::Display for PhsMimoError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear original rational boundary: {self:?}")
    }
}
impl std::error::Error for PhsMimoError {}
impl From<PhsOrbitError> for PhsMimoError {
    fn from(error: PhsOrbitError) -> Self {
        Self::Regime(error)
    }
}
impl From<MimoLoadError> for PhsMimoError {
    fn from(error: MimoLoadError) -> Self {
        Self::Source(error)
    }
}
type Result<T> = core::result::Result<T, PhsMimoError>;

struct Scratch {
    host_effort: Vec<f64>,
    direction: Vec<f64>,
    host_rate: Vec<f64>,
    source_rate: Vec<f64>,
    flow: Vec<f64>,
    boundary_effort: Vec<f64>,
}

/// Prepared complete [original host, original real source] ODE.
///
/// Source pair coordinates use the original MimoModel chart: all real parts,
/// then all MINUS imaginary parts for each conjugate pole pair. No accepted
/// history is silently inserted into a steady periodic seed. Use `initial_state`
/// explicitly to carry existing source memory into a time-domain seed.
pub struct PreparedPhsMimo<'a, 's> {
    original: PhsOrbitProblem<'a>,
    source: &'a PreparedMimoLoad,
    input: &'a PeriodicPortInput,
    boundary: Boundary,
    source_reaction: Vec<f64>,
    limits: PhsMimoLimits,
    cx: &'a Cx<'s>,
    scratch: RefCell<Scratch>,
    zero_input: Vec<f64>,
    tangent_cache: RefCell<Option<Vec<f64>>>,
    products: Cell<u64>,
    budget_refused: Cell<bool>,
    cancelled: Cell<bool>,
}

impl<'a, 's> PreparedPhsMimo<'a, 's> {
    /// Attach every original rational boundary port using finite row-major
    /// host-state-by-source-port columns. These need not be host drive columns.
    /// Complete source A/B/C/D/E comes from its existing continuous owner.
    /// E MUST be exactly zero: a differentiating load needs an explicit
    /// descriptor/mass owner and cannot enter this bounded explicit ODE path.
    pub fn new(
        system: &'a PortHamiltonian,
        source: &'a PreparedMimoLoad,
        input: &'a PeriodicPortInput,
        columns: &[f64],
        limits: PhsMimoLimits,
        cx: &'a Cx<'s>,
    ) -> Result<Self> {
        cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
        let view = source.continuous_operator();
        let (n, q, p, m) = (
            system.state_dim(),
            view.dimension(),
            view.ports(),
            system.port_dim(),
        );
        let d = n
            .checked_add(q)
            .ok_or(PhsMimoError::Invalid("joint state extent overflow"))?;
        if n == 0
            || d > 256
            || d > limits.maximum_states
            || p == 0
            || p > 128
            || p > limits.maximum_ports
            || m > 256
            || input.port_count() != m
            || n.checked_mul(p) != Some(columns.len())
            || columns.iter().any(|v| !v.is_finite())
            || input.harmonic_count() > limits.maximum_harmonics
            || input.harmonic_count() > 256
            || !(core::f64::consts::TAU / input.omega_rad_s()).is_finite()
        {
            return Err(PhsMimoError::Invalid(
                "bounded complete original states, physical ports and forcing required",
            ));
        }
        // Includes source matrices, original tangent metadata, a complete cached
        // joint Jacobian, all reusable scratch and an independent audit workspace.
        let bytes =
            8 * (q * q + 2 * q * p + 2 * p * p + n * p) + 32 * d * d + 128 * d + 64 * p + 32 * m;
        if bytes > limits.maximum_workspace_bytes || bytes > 256 * 1024 * 1024 {
            return Err(PhsMimoError::Budget(
                "original-state boundary payload exceeds its cap",
            ));
        }
        let state = view.state_space(limits.maximum_workspace_bytes, &mut || {
            cx.checkpoint().is_err()
        });
        cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
        let state = state?;
        if state.n != q
            || state.ports != p
            || state.a.len() != q * q
            || state.b.len() != q * p
            || state.c.len() != p * q
            || state.d.len() != p * p
            || state.e.len() != p * p
            || state.e.iter().any(|&e| e != 0.0)
        {
            return Err(PhsMimoError::Invalid(
                "complete proper original rational realization required; nonzero E needs a descriptor owner",
            ));
        }
        if n as u64 * q as u64 * p as u64 > limits.maximum_products {
            return Err(PhsMimoError::Budget(
                "constant original source tangent exceeds its preparation cap",
            ));
        }
        let boundary = Boundary {
            host: n,
            states: q,
            ports: p,
            columns: columns.to_vec(),
            a: state.a,
            b: state.b,
            c: state.c,
            d: state.d,
        };
        let mut source_reaction = vec![0.0; n * q];
        for row in 0..n {
            cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
            if !boundary.state_reaction_row(row, &mut source_reaction[row * q..(row + 1) * q]) {
                return Err(PhsMimoError::Invalid(
                    "unrepresentable original source reaction tangent",
                ));
            }
        }
        Ok(Self {
            original: forced_problem(system, input)?,
            source,
            input,
            boundary,
            source_reaction,
            limits,
            cx,
            scratch: RefCell::new(Scratch {
                host_effort: vec![0.0; n],
                direction: vec![0.0; n],
                host_rate: vec![0.0; n],
                source_rate: vec![0.0; q],
                flow: vec![0.0; p],
                boundary_effort: vec![0.0; p],
            }),
            zero_input: vec![0.0; m],
            tangent_cache: RefCell::new(None),
            products: Cell::new(0),
            budget_refused: Cell::new(false),
            cancelled: Cell::new(false),
        })
    }

    /// Original admitted source, including unchanged band, clock and history.
    pub fn source(&self) -> &'a PreparedMimoLoad {
        self.source
    }
    /// Original continuous and sampled reports; neither becomes a storage witness.
    pub fn reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        self.source.reports()
    }
    /// Complete original host state count, including material memory.
    pub fn host_dimension(&self) -> usize {
        self.boundary.host
    }
    /// Complete original real rational state count, including unobserved modes.
    pub fn source_dimension(&self) -> usize {
        self.boundary.states
    }
    /// State/row extent for every shared numerical query, in [host, source] order.
    pub fn dimension(&self) -> usize {
        self.boundary.host + self.boundary.states
    }

    /// Supply a complete original-storage Hessian action when the source does
    /// not already expose one. Refusal never silently switches to differences.
    #[must_use]
    pub fn with_storage_hessian(mut self, hessian: &'a PhsStorageHessianAction<'a>) -> Self {
        self.original.tangent = StorageTangent::hessian(self.boundary.host, hessian);
        self.tangent_cache = RefCell::new(None);
        self
    }

    /// Explicitly append the SAME source's accepted real pole coordinates.
    pub fn initial_state(&self, host_state: &[f64]) -> Result<Vec<f64>> {
        if host_state.len() != self.boundary.host || host_state.iter().any(|v| !v.is_finite()) {
            return Err(PhsMimoError::Invalid(
                "complete finite original host state required",
            ));
        }
        self.cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
        let mut state = Vec::with_capacity(self.dimension());
        state.extend_from_slice(host_state);
        state.extend(self.source.continuous_operator().real_coordinates());
        self.cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
        Ok(state)
    }

    /// Candidate of the original source at a caller-accepted joint endpoint.
    /// This uses the existing source history-conversion owner. No source is
    /// mutated; the enclosing trajectory owns acceptance and publication.
    pub fn source_candidate(&self, joint_state: &[f64]) -> Result<PreparedMimoLoad> {
        self.reset()?;
        self.observe(joint_state)?;
        let scratch = self.scratch.borrow();
        let result = self
            .source
            .continuous_operator()
            .candidate_from_coordinates(&joint_state[self.boundary.host..], &scratch.flow)
            .map_err(PhsMimoError::Source);
        self.after()?;
        result
    }

    fn reset(&self) -> Result<()> {
        self.cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
        self.products.set(0);
        self.budget_refused.set(false);
        self.cancelled.set(false);
        self.original.source_refused.set(false);
        Ok(())
    }
    fn poll(&self, products: u64) -> bool {
        if self.cx.checkpoint().is_err() {
            self.cancelled.set(true);
            return false;
        }
        let Some(total) = self.products.get().checked_add(products) else {
            self.budget_refused.set(true);
            return false;
        };
        if total > self.limits.maximum_products {
            self.budget_refused.set(true);
            return false;
        }
        self.products.set(total);
        true
    }
    fn after(&self) -> Result<()> {
        if self.cancelled.get() || self.cx.checkpoint().is_err() {
            return Err(PhsMimoError::Cancelled);
        }
        if self.budget_refused.get() {
            return Err(PhsMimoError::Budget(
                "original-state boundary contraction cap exhausted",
            ));
        }
        if self.original.source_refused() {
            return Err(PhsOrbitError::SourceRefused.into());
        }
        Ok(())
    }
    fn source_work(&self) -> u64 {
        let (n, q, p) = (
            self.boundary.host as u64,
            self.boundary.states as u64,
            self.boundary.ports as u64,
        );
        3 * n * p + q * q + 2 * q * p + p * p + 4 * (n + q + p)
    }
    fn rhs_work(&self) -> u64 {
        let (n, m) = (
            self.boundary.host as u64,
            self.original.system.port_dim() as u64,
        );
        self.source_work()
            + 4 * n * n
            + 2 * n * m
            + 32 * m * (self.input.harmonic_count() as u64 + 1)
    }
    fn tangent_work(&self) -> u64 {
        let (n, q, p, m) = (
            self.boundary.host as u64,
            self.boundary.states as u64,
            self.boundary.ports as u64,
            self.original.system.port_dim() as u64,
        );
        4 * n * n + 2 * n * m + 3 * n * p + q * p + p * p + 8 * (n + q + p)
    }
    fn admit_harmonics(&self, maximum: usize) -> Result<()> {
        if maximum == 0 || maximum > 256 || maximum > self.limits.maximum_harmonics {
            return Err(PhsMimoError::Invalid(
                "positive bounded retained harmonic order required",
            ));
        }
        for h in 1..=maximum {
            self.cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
            self.source
                .continuous_operator()
                .admit_frequency(self.input.omega_rad_s() * h as f64)?;
        }
        Ok(())
    }
    fn admit_hb(&self, budget: &HbBudget, maximum: usize) -> Result<()> {
        self.reset()?;
        if budget.harmonics > maximum || self.input.highest_nonzero_harmonic() > budget.harmonics {
            return Err(PhsMimoError::Invalid(
                "drive must fit the initial retained harmonic order",
            ));
        }
        self.admit_harmonics(maximum)?;
        let unknowns = self.dimension() * (2 * maximum + 1);
        if unknowns > 4096 || unknowns > self.limits.maximum_hb_unknowns {
            return Err(PhsMimoError::Budget(
                "complete original-state Fourier system exceeds its cap",
            ));
        }
        Ok(())
    }

    /// Fixed-order HB of the full original [host, source] ODE. Fourier seeds
    /// and scales retain all source states; no accepted history is added as an
    /// extra periodic forcing. Original fs-orbit owns Newton/AFT convergence.
    pub fn solve_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        scales: &HbScales,
    ) -> Result<ScaledHbOrbit> {
        self.admit_hb(budget, budget.harmonics)?;
        let result = solve_hb_scaled_seeded(
            self,
            HbAnchor::Forced {
                omega: self.input.omega_rad_s(),
            },
            seed,
            budget,
            scales,
        );
        self.after()?;
        result.map_err(|e| PhsOrbitError::Orbit(e).into())
    }

    /// Shared adaptive harmonic refinement with unchanged source law and every
    /// original state. All retained orders must fit the source's positive band;
    /// sampled omitted nonlinear content remains a numerical refinement measure.
    pub fn solve_adaptive_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        refinement: &HbRefinementBudget,
        scales: &HbScales,
    ) -> Result<AdaptiveHbOrbit> {
        self.admit_hb(budget, refinement.max_harmonics)?;
        let result = solve_hb_adaptive_scaled_seeded(
            self,
            HbAnchor::Forced {
                omega: self.input.omega_rad_s(),
            },
            seed,
            budget,
            refinement,
            scales,
        );
        self.after()?;
        result.map_err(|e| PhsOrbitError::Orbit(e).into())
    }

    /// Shoot the actual full causal model and retain every original Floquet
    /// multiplier, including source modes invisible at the selected boundary.
    /// Requires the complete native/caller storage Hessian. Shared controls
    /// bound the numerical map; PhsMimoLimits separately bounds adapter work.
    /// One RK4 grid is numerical evidence, not a temporal convergence certificate.
    pub fn solve_shooting_variational(
        &self,
        seed: &[f64],
        budget: &ShootBudget,
        controls: &VariationalShootControls,
    ) -> Result<VariationalShootOrbit> {
        self.reset()?;
        if budget.steps_per_period == 0
            || self.input.highest_nonzero_harmonic() > (budget.steps_per_period - 1) / 2
        {
            return Err(PhsMimoError::Invalid(
                "shooting grid must resolve every authored drive harmonic",
            ));
        }
        let highest = self.input.highest_nonzero_harmonic();
        if highest > 0 {
            self.admit_harmonics(highest)?;
        }
        let result = solve_shooting_variational(
            self,
            seed,
            core::f64::consts::TAU / self.input.omega_rad_s(),
            budget,
            controls,
            || self.cx.checkpoint().is_err(),
        );
        self.after()?;
        result.map_err(|e| PhsOrbitError::Variational(e).into())
    }

    fn observe(&self, state: &[f64]) -> Result<()> {
        if !self.poll(self.source_work()) {
            self.after()?;
        }
        if state.len() != self.dimension() || state.iter().any(|v| !v.is_finite()) {
            return Err(PhsMimoError::Invalid(
                "complete finite joint state required",
            ));
        }
        let n = self.boundary.host;
        let mut scratch = self.scratch.borrow_mut();
        let Scratch {
            host_effort,
            flow,
            boundary_effort,
            ..
        } = &mut *scratch;
        host_effort.fill(f64::NAN);
        self.original
            .system
            .storage()
            .gradient(&state[..n], host_effort);
        if host_effort.iter().any(|v| !v.is_finite())
            || !self.boundary.flow_into(host_effort, flow)
            || !self
                .boundary
                .effort_into(&state[n..], flow, boundary_effort)
        {
            self.original.source_refused.set(true);
            return Err(PhsOrbitError::SourceRefused.into());
        }
        Ok(())
    }

    fn state_jacobian(&self, _time: f64, state: &[f64], out: &mut [f64]) -> bool {
        let (n, q, d) = (self.boundary.host, self.boundary.states, self.dimension());
        if !self.poll((d * d) as u64)
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
        let mut scratch = self.scratch.borrow_mut();
        let Scratch {
            host_effort,
            direction,
            host_rate,
            source_rate,
            flow,
            boundary_effort,
        } = &mut *scratch;
        out.fill(0.0);
        for column in 0..n {
            if !self.poll(self.tangent_work()) {
                return false;
            }
            direction.fill(0.0);
            direction[column] = 1.0;
            if !self
                .original
                .tangent
                .apply_hessian(&state[..n], direction, host_effort)
                || !self.boundary.flow_into(host_effort, flow)
                || !self.boundary.direct_effort_into(flow, boundary_effort)
                || self
                    .original
                    .system
                    .vector_field_with_reaction_from_effort_into(
                        host_effort,
                        &self.zero_input,
                        &self.boundary.columns,
                        boundary_effort,
                        host_rate,
                    )
                    .is_err()
                || !self.boundary.driven_rate_into(flow, source_rate)
            {
                self.original.source_refused.set(true);
                return false;
            }
            for row in 0..n {
                out[row * d + column] = host_rate[row];
            }
            for row in 0..q {
                out[(n + row) * d + column] = source_rate[row];
            }
        }
        // Copy constant -S C and A: O(d*q), retaining all dark source modes.
        // Neither block needs q dense source-rate sweeps at every RK4 stage.
        for row in 0..n {
            out[row * d + n..(row + 1) * d]
                .copy_from_slice(&self.source_reaction[row * q..(row + 1) * q]);
        }
        for row in 0..q {
            out[(n + row) * d + n..(n + row + 1) * d]
                .copy_from_slice(&self.boundary.a[row * q..(row + 1) * q]);
        }
        if matches!(self.original.tangent, StorageTangent::Quadratic { .. }) {
            *self.tangent_cache.borrow_mut() = Some(out.to_vec());
        }
        true
    }
}

impl OrbitProblem for PreparedPhsMimo<'_, '_> {
    fn dim(&self) -> usize {
        self.dimension()
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
    fn time_domain_rhs(&self) -> Option<fs_orbit::TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<fs_orbit::TimeDomainStateJacobian<Self>> {
        self.original
            .tangent
            .available()
            .then_some(Self::state_jacobian)
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        let n = self.boundary.host;
        if !time.is_finite()
            || !self.poll(self.rhs_work())
            || state.len() != self.dimension()
            || out.len() != self.dimension()
            || state.iter().any(|v| !v.is_finite())
        {
            self.original.fail(out);
            return;
        }
        let mut effort = self.original.effort.borrow_mut();
        effort.fill(f64::NAN);
        self.original
            .system
            .storage()
            .gradient(&state[..n], &mut effort);
        let mut input = self.original.ports.borrow_mut();
        if effort.iter().any(|v| !v.is_finite()) || !self.input.sample_into(time, &mut input) {
            self.original.fail(out);
            return;
        }
        let mut scratch = self.scratch.borrow_mut();
        let Scratch {
            flow,
            boundary_effort,
            ..
        } = &mut *scratch;
        if !self.boundary.flow_into(&effort, flow)
            || !self
                .boundary
                .effort_into(&state[n..], flow, boundary_effort)
            || self
                .original
                .system
                .vector_field_with_reaction_from_effort_into(
                    &effort,
                    &input,
                    &self.boundary.columns,
                    boundary_effort,
                    &mut out[..n],
                )
                .is_err()
            || !self.boundary.rate_into(&state[n..], flow, &mut out[n..])
        {
            self.original.fail(out);
        }
    }
}

/// Independent sampled original-host energy and signed rational boundary work.
#[derive(Clone, Copy, Debug)]
pub struct PhsMimoAudit {
    /// Independent uniform quadrature points over the physical forcing period.
    pub samples: usize,
    /// Original host Hamiltonian only; no rational storage is inferred.
    pub mean_host_stored_energy_j: f64,
    /// Minimum original host Hamiltonian on this quadrature grid.
    pub minimum_host_stored_energy_j: f64,
    /// Maximum original host Hamiltonian on this quadrature grid.
    pub maximum_host_stored_energy_j: f64,
    /// Authored external effort supplied through the original host G.
    pub mean_supplied_power_w: f64,
    /// Original host gradient dotted with its original R gradient.
    pub mean_host_resistance_loss_w: f64,
    /// Original source effort dotted with its actual flow; may be negative.
    pub mean_source_exchange_power_w: f64,
    /// Independent host-side physical reaction work.
    pub mean_host_boundary_power_w: f64,
    /// Source-side exchange minus independently contracted host boundary work.
    pub mean_exchange_residual_w: f64,
    /// Applied supply minus original host loss minus signed source work.
    pub cycle_power_residual_w: f64,
    /// Absolute cycle residual divided by the mean sum of absolute supply,
    /// host loss and both independently evaluated boundary works; zero at rest.
    pub relative_power_residual: f64,
}

impl PreparedPhsMimo<'_, '_> {
    /// Sample the SAME original H/R/G and source C/D over a full joint orbit.
    /// The rational states carry no invented heat/storage. This finite quadrature
    /// is independent of Newton's residual and needs a sample ladder for nonlinear
    /// storage. It does not certify the source outside its original physical band.
    pub fn audit(&self, orbit: &HbOrbit, samples: usize) -> Result<PhsMimoAudit> {
        self.reset()?;
        let width = orbit.coeffs.first().map_or(0, Vec::len);
        if orbit.omega.to_bits() != self.input.omega_rad_s().to_bits()
            || orbit.coeffs.len() != self.dimension()
            || width == 0
            || !(32..=65_536).contains(&samples)
            || width - 1 > (samples - 1) / 2
            || self.input.highest_nonzero_harmonic() > (samples - 1) / 2
            || orbit.coeffs.iter().any(|row| {
                row.len() != width
                    || row[0].im != 0.0
                    || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
        {
            return Err(PhsMimoError::Invalid(
                "audit needs complete finite joint Fourier rows and a resolved 32..=65536 grid",
            ));
        }
        if width > 1 {
            self.admit_harmonics(width - 1)?;
        }
        let n = self.boundary.host;
        let mut state = vec![0.0; self.dimension()];
        let mut input = vec![0.0; self.input.port_count()];
        let (_, r, _) = self.original.system.structure();
        let mut sums = [Dd::ZERO; 6];
        let (mut minimum, mut maximum) = (f64::INFINITY, f64::NEG_INFINITY);
        for sample in 0..samples {
            let audit_work = 4 * n * n
                + 2 * n * self.input.port_count()
                + 32 * self.dimension() * width
                + 32 * self.input.port_count() * (self.input.harmonic_count() + 1);
            if !self.poll(audit_work as u64) {
                self.after()?;
            }
            let phase = core::f64::consts::TAU * sample as f64 / samples as f64;
            for (coordinate, x) in state.iter_mut().enumerate() {
                *x = orbit.sample(coordinate, phase);
            }
            self.observe(&state)?;
            let scratch = self.scratch.borrow();
            let effort = &scratch.host_effort;
            if !self.input.sample_into(phase / orbit.omega, &mut input) {
                return Err(PhsOrbitError::SourceRefused.into());
            }
            let stored = self.original.system.hamiltonian(&state[..n]);
            let supplied = self
                .original
                .system
                .supplied_power_from_effort(effort, &input)
                .map_err(|_| PhsOrbitError::SourceRefused)?;
            let mut loss = Dd::ZERO;
            for i in 0..n {
                for j in 0..n {
                    loss = loss
                        + Dd::from_f64(effort[i])
                            * Dd::from_f64(r[i * n + j])
                            * Dd::from_f64(effort[j]);
                }
            }
            let loss = loss.to_f64();
            let source = dot(&scratch.flow, &scratch.boundary_effort);
            let host = self.boundary.host_power(effort, &scratch.boundary_effort);
            let scale = supplied.abs() + loss.abs() + source.abs() + host.abs();
            let values = [stored, supplied, loss, source, host, scale];
            if values.iter().any(|v| !v.is_finite()) {
                return Err(PhsMimoError::Invalid(
                    "unrepresentable original storage or boundary work",
                ));
            }
            for (sum, value) in sums.iter_mut().zip(values) {
                *sum = *sum + Dd::from_f64(value);
            }
            minimum = minimum.min(stored);
            maximum = maximum.max(stored);
        }
        let mean = sums.map(|sum| sum.to_f64() / samples as f64);
        let residual =
            (Dd::from_f64(mean[1]) - Dd::from_f64(mean[2]) - Dd::from_f64(mean[3])).to_f64();
        let relative = if mean[5] == 0.0 {
            0.0
        } else {
            residual.abs() / mean[5]
        };
        let exchange = (Dd::from_f64(mean[3]) - Dd::from_f64(mean[4])).to_f64();
        if mean
            .iter()
            .chain([residual, relative, exchange].iter())
            .any(|v| !v.is_finite())
        {
            return Err(PhsMimoError::Invalid(
                "unrepresentable cycle work reduction",
            ));
        }
        self.after()?;
        Ok(PhsMimoAudit {
            samples,
            mean_host_stored_energy_j: mean[0],
            minimum_host_stored_energy_j: minimum,
            maximum_host_stored_energy_j: maximum,
            mean_supplied_power_w: mean[1],
            mean_host_resistance_loss_w: mean[2],
            mean_source_exchange_power_w: mean[3],
            mean_host_boundary_power_w: mean[4],
            mean_exchange_residual_w: exchange,
            cycle_power_residual_w: residual,
            relative_power_residual: relative,
        })
    }
}
