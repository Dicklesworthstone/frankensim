//! Complete native frame and attached-memory histories on the existing exact flow.
//!
//! Only an ordinary differential descriptor is admitted here. Its original
//! J/R/G and storage are retained; no impedance, modal damping or new generator
//! is inferred. The common quadratic flow integrates every heat and effort
//! observation on the same held-force trajectory.

use super::{FrameCoupledError, PreparedFrameResonance, limit, poll};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, QuadraticPart, split_quadratic,
};
use fs_phs::{DescriptorPortHamiltonian, PortHamiltonian, PreparedStepError, Storage};
use std::ops::Range;
use std::sync::Arc;

type Result<T> = core::result::Result<T, FrameCoupledError>;

/// Complete fixed-clock preparation bounds; no coordinate is discarded to fit.
#[derive(Clone, Copy, Debug)]
pub struct FrameTransientControls {
    /// Includes the adapter's numerical payload and preparation work before
    /// passing the remaining allowance to the original exponential/Gramian owner.
    /// The already prepared source, caller inputs and allocator overhead are excluded.
    pub exact: DenseExponentialBudget,
    /// Relative energy and original-equation gate, finite in (0, 1e-8].
    pub maximum_relative_balance_error: f64,
}

/// One accepted full native endpoint and this interval's independent accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameTransientReport {
    pub time_s: f64,
    /// Original node-major [ux, uy, theta]; stationary supported entries are zero.
    pub displacement: Vec<f64>,
    pub velocity: Vec<f64>,
    pub kinetic_energy_j: f64,
    pub equilibrium_energy_j: f64,
    /// All original attached states, including constitutive memory when present.
    pub attached_energy_j: f64,
    pub total_energy_j: f64,
    /// Polarized change of the original complete storage, independent of heat/work.
    pub energy_change_j: f64,
    /// Held original nodal forces/moments times actual displacement increments.
    pub applied_work_j: f64,
    /// Supports in this source are stationary and do exactly zero work.
    /// Dynamic support reactions are not reconstructed from the reduced mass.
    pub stationary_support_work_j: f64,
    /// Integral of the original full v^T C v, including explicitly authored
    /// Rayleigh loss and unselected members' material viscosity exactly once.
    pub frame_dissipation_j: f64,
    /// Integral of the original attached effort^T R effort. This includes any
    /// material memory; mechanism-specific reports are subaccounts of this heat.
    pub attached_dissipation_j: f64,
    pub energy_balance_residual_j: f64,
    /// Original native integrated q/p equation residual, scaled componentwise.
    pub maximum_relative_equation_residual: f64,
}

impl FrameTransientReport {
    fn zero(ndof: usize) -> Self {
        Self {
            time_s: 0.0,
            displacement: vec![0.0; ndof],
            velocity: vec![0.0; ndof],
            kinetic_energy_j: 0.0,
            equilibrium_energy_j: 0.0,
            attached_energy_j: 0.0,
            total_energy_j: 0.0,
            energy_change_j: 0.0,
            applied_work_j: 0.0,
            stationary_support_work_j: 0.0,
            frame_dissipation_j: 0.0,
            attached_dissipation_j: 0.0,
            energy_balance_residual_j: 0.0,
            maximum_relative_equation_residual: 0.0,
        }
    }

    fn clear_interval(&mut self) {
        self.energy_change_j = 0.0;
        self.applied_work_j = 0.0;
        self.stationary_support_work_j = 0.0;
        self.frame_dissipation_j = 0.0;
        self.attached_dissipation_j = 0.0;
        self.energy_balance_residual_j = 0.0;
        self.maximum_relative_equation_residual = 0.0;
    }
}

// Forward the SAME immutable storage callbacks, including its exact declaration
// and stable energy remainder. This adapter does not rebuild a material law.
struct OriginalStorage(Arc<DescriptorPortHamiltonian>);
impl Storage for OriginalStorage {
    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(Self(Arc::clone(&self.0))))
    }
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.0.storage().hamiltonian(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.0.storage().gradient(x, out);
    }
    fn continuous_tangent(&self) -> Option<&dyn fs_phs::ContinuousStorageTangent> {
        self.0.storage().continuous_tangent()
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn fs_phs::ConstitutiveDiscreteGradient> {
        self.0.storage().discrete_gradient_rule()
    }
    fn coordinate_scale(&self, coordinate: usize) -> Option<f64> {
        self.0.storage().coordinate_scale(coordinate)
    }
    fn quadratic_part(&self) -> Option<QuadraticPart> {
        self.0.storage().quadratic_part()
    }
    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.0.storage().exact_quadratic_matrix()
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        self.0.storage().midpoint_energy_remainder(a, b, mid)
    }
}

/// Fixed-clock complete native history. No successful warm step allocates or
/// prepares an exponential; accepted state and observations publish together.
/// Original reference geometry/inertia, small strain and stationary supports
/// apply. No thermal deposition, evolving material or moving support is inferred.
pub struct PreparedFrameTransient<'a> {
    source: &'a PreparedFrameResonance,
    system: PortHamiltonian,
    flow: PreparedQuadraticFlow,
    duration_s: f64,
    tolerance: f64,
    state: Vec<f64>,
    candidate: Vec<f64>,
    input: Vec<f64>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    integrals: Vec<f64>,
    effort: Vec<f64>,
    old_effort: Vec<f64>,
    heat: Vec<f64>,
    staged_heat: Vec<f64>,
    report: FrameTransientReport,
    staged: FrameTransientReport,
}

impl PreparedFrameResonance {
    /// Retain every native q/p and original attached state in one exact
    /// held-force evolution. Material sources must explicitly support the
    /// broadband continuous law [0, f64::MAX] rad/s. Boundary applicability
    /// remains with its original owner; exact propagation supplies no new data.
    /// Algebraic boundary constraints refuse rather than acquiring fake inertia.
    pub fn prepare_transient(
        &self,
        duration_s: f64,
        controls: FrameTransientControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedFrameTransient<'_>> {
        let nf = self.free_dofs.len();
        PreparedFrameTransient::prepare(
            self,
            duration_s,
            controls,
            vec![nf..2 * nf, 2 * nf..self.dimension()],
            cx,
        )
    }
}

impl<'a> PreparedFrameTransient<'a> {
    pub(in crate::frame2d) fn prepare(
        source: &'a PreparedFrameResonance,
        duration_s: f64,
        controls: FrameTransientControls,
        groups: Vec<Range<usize>>,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let n = source.dimension();
        let nf = source.free_dofs.len();
        if source.system.differential_dim() != n {
            return Err(invalid(
                "exact native frame time evolution requires differential boundary states",
            ));
        }
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || !controls.maximum_relative_balance_error.is_finite()
            || controls.maximum_relative_balance_error <= 0.0
            || controls.maximum_relative_balance_error > 1e-8
            || source
                .material_band
                .is_some_and(|band| band != (0.0, f64::MAX))
            || groups.len() < 2
            || groups.len() > n + 1
            || groups[0] != (nf..2 * nf)
        {
            return Err(invalid(
                "positive clock, explicit broadband material support and bounded original heat groups required",
            ));
        }
        let d = n as u128 + nf as u128 + 1;
        limit(
            "exact frame augmented states",
            2 * d,
            controls.exact.max_dimension.min(256) as u128,
        )?;
        // Original operators/storage remain borrowed. Charge copied structure,
        // full source heat forms, split scratch, physical endpoints and staging.
        let bytes = 8
            * ((groups.len() as u128 + 24) * (n as u128).pow(2)
                + 4 * n as u128 * nf as u128
                + 24 * n as u128
                + 4 * source.ndof as u128);
        let work = (4096 + 4 * groups.len() as u128) * (n as u128).pow(3);
        limit(
            "exact frame adapter bytes",
            bytes,
            controls.exact.max_dense_bytes as u128,
        )?;
        limit(
            "exact frame adapter work",
            work,
            controls.exact.max_multiply_terms as u128,
        )?;
        let cache_budget = DenseExponentialBudget {
            max_dimension: controls.exact.max_dimension,
            max_dense_bytes: controls.exact.max_dense_bytes - bytes as usize,
            max_multiply_terms: controls.exact.max_multiply_terms - work as u64,
        };
        let (j, r, g) = source.system.structure();
        let mut membership = vec![usize::MAX; n];
        for (group, range) in groups.iter().enumerate() {
            if range.start > range.end || range.end > n {
                return Err(invalid("heat group is outside the original frame state"));
            }
            for row in range.clone() {
                if membership[row] != usize::MAX {
                    return Err(invalid(
                        "original resistance cannot be counted by two heat groups",
                    ));
                }
                membership[row] = group;
            }
        }
        for i in 0..n {
            poll(cx)?;
            for k in 0..n {
                if r[i * n + k] != 0.0
                    && (membership[i] == usize::MAX || membership[i] != membership[k])
                {
                    return Err(invalid(
                        "heat partition cannot discard an original resistance cross term",
                    ));
                }
            }
        }
        let system = PortHamiltonian::new(
            n,
            nf,
            j.to_vec(),
            r.to_vec(),
            g.to_vec(),
            Box::new(OriginalStorage(Arc::clone(&source.system))),
        )
        .map_err(FrameCoupledError::Physical)?;
        poll(cx)?;
        let split = split_quadratic(
            &system,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: core::f64::consts::PI,
            },
        )
        .map_err(FrameCoupledError::Physical)?;
        poll(cx)?;
        let q = split.quadratic_matrix();
        let forms = groups
            .iter()
            .map(|range| resistance_form(q, r, n, range.clone(), cx))
            .collect::<Result<Vec<_>>>()?;
        // Integrate all original efforts, retaining material/boundary reaction
        // terms in the independently checked native momentum equation.
        let observations: Vec<_> = q.chunks_exact(n).map(<[f64]>::to_vec).collect();
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            duration_s,
            &forms,
            &observations,
            cache_budget,
            || cx.checkpoint().is_err(),
        )
        .map_err(step_error)?;
        let report = FrameTransientReport::zero(source.ndof);
        let staged = report.clone();
        let mut result = Self {
            source,
            system,
            flow,
            duration_s,
            tolerance: controls.maximum_relative_balance_error,
            state: vec![0.0; n],
            candidate: vec![0.0; n],
            input: vec![0.0; nf],
            augmented: vec![0.0; n + nf + 1],
            balanced: vec![0.0; n],
            integrals: vec![0.0; groups.len() + n],
            effort: vec![0.0; n],
            old_effort: vec![0.0; n],
            heat: vec![0.0; groups.len()],
            staged_heat: vec![0.0; groups.len()],
            report,
            staged,
        };
        result.stage_endpoint(0.0, cx)?;
        poll(cx)?;
        result.publish_staged();
        Ok(result)
    }

    pub fn source(&self) -> &PreparedFrameResonance {
        self.source
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Original [q, p=M*v, attached state] including every material history.
    pub fn canonical_state(&self) -> &[f64] {
        &self.state
    }
    pub fn report(&self) -> &FrameTransientReport {
        &self.report
    }

    /// Restore the complete state and physical clock on this unchanged source.
    /// Neither history nor accepted observations change if admission fails.
    pub fn restore_state(&mut self, time_s: f64, state: &[f64], cx: &Cx<'_>) -> Result<()> {
        self.stage_restore(time_s, state, cx)?;
        poll(cx)?;
        self.publish_staged();
        Ok(())
    }

    /// Original full nodal displacements/velocities and complete attached state.
    /// Fixed entries must be zero; no unsupported prescribed support motion enters.
    pub fn initialize(
        &mut self,
        time_s: f64,
        displacement: &[f64],
        velocity: &[f64],
        attached_state: &[f64],
        cx: &Cx<'_>,
    ) -> Result<()> {
        self.stage_initialize(time_s, displacement, velocity, attached_state, cx)?;
        poll(cx)?;
        self.publish_staged();
        Ok(())
    }

    /// Advance one held full nodal force/moment. Fixed loads do no work. Every
    /// off-diagonal inertial, elastic, viscous and attached-memory term survives.
    pub fn step(&mut self, nodal_force: &[f64], cx: &Cx<'_>) -> Result<&FrameTransientReport> {
        self.stage_step(nodal_force, cx)?;
        poll(cx)?;
        self.publish_staged();
        Ok(&self.report)
    }

    pub(in crate::frame2d) fn stage_restore(
        &mut self,
        time_s: f64,
        state: &[f64],
        cx: &Cx<'_>,
    ) -> Result<()> {
        poll(cx)?;
        if state.len() != self.state.len() || state.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "frame restart requires every finite original state",
            ));
        }
        self.candidate.copy_from_slice(state);
        self.stage_endpoint(time_s, cx)
    }

    pub(in crate::frame2d) fn stage_initialize(
        &mut self,
        time_s: f64,
        displacement: &[f64],
        velocity: &[f64],
        attached_state: &[f64],
        cx: &Cx<'_>,
    ) -> Result<()> {
        poll(cx)?;
        let nf = self.source.free_dofs.len();
        if displacement.len() != self.source.ndof
            || velocity.len() != self.source.ndof
            || attached_state.len() != self.source.load_dimension
            || displacement
                .iter()
                .chain(velocity)
                .chain(attached_state)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "initial frame motion and attached history must match original coordinates",
            ));
        }
        let mut free = 0;
        for dof in 0..self.source.ndof {
            if self.source.free_dofs.get(free) == Some(&dof) {
                free += 1;
            } else if displacement[dof] != 0.0 || velocity[dof] != 0.0 {
                return Err(invalid(
                    "stationary frame supports cannot carry initial motion",
                ));
            }
        }
        for (i, &dof) in self.source.free_dofs.iter().enumerate() {
            poll(cx)?;
            self.candidate[i] = displacement[dof];
            self.candidate[nf + i] = sum(self
                .source
                .free_dofs
                .iter()
                .enumerate()
                .map(|(k, &dof)| self.source.mass[i * nf + k] * velocity[dof]))?;
        }
        self.candidate[2 * nf..].copy_from_slice(attached_state);
        self.stage_endpoint(time_s, cx)
    }

    fn stage_endpoint(&mut self, time_s: f64, cx: &Cx<'_>) -> Result<()> {
        poll(cx)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(invalid("frame history requires a nonnegative finite time"));
        }
        self.system
            .storage()
            .gradient(&self.candidate, &mut self.effort);
        if self.effort.iter().any(|x| !x.is_finite()) {
            return Err(invalid("original frame effort is unrepresentable"));
        }
        let nf = self.source.free_dofs.len();
        for (i, &dof) in self.source.free_dofs.iter().enumerate() {
            poll(cx)?;
            self.staged.displacement[dof] = self.candidate[i];
            self.staged.velocity[dof] = self.effort[nf + i];
        }
        self.staged.time_s = time_s;
        self.staged.equilibrium_energy_j = 0.5 * dot(&self.candidate[..nf], &self.effort[..nf])?;
        self.staged.kinetic_energy_j =
            0.5 * dot(&self.candidate[nf..2 * nf], &self.effort[nf..2 * nf])?;
        self.staged.attached_energy_j =
            0.5 * dot(&self.candidate[2 * nf..], &self.effort[2 * nf..])?;
        self.staged.total_energy_j = finite(self.system.hamiltonian(&self.candidate))?;
        self.staged.clear_interval();
        self.staged_heat.fill(0.0);
        poll(cx)
    }

    pub(in crate::frame2d) fn stage_step(
        &mut self,
        nodal_force: &[f64],
        cx: &Cx<'_>,
    ) -> Result<()> {
        poll(cx)?;
        if nodal_force.len() != self.source.ndof || nodal_force.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "one finite held force/moment per original frame coordinate required",
            ));
        }
        let time_s = finite(self.report.time_s + self.duration_s)?;
        if time_s <= self.report.time_s {
            return Err(invalid("frame clock cannot represent this interval"));
        }
        for (i, &dof) in self.source.free_dofs.iter().enumerate() {
            self.input[i] = nodal_force[dof];
        }
        self.flow
            .evaluate_into(
                &self.state,
                &self.input,
                &mut self.augmented,
                &mut self.balanced,
                &mut self.candidate,
                &mut self.integrals,
            )
            .map_err(FrameCoupledError::Physical)?;
        self.stage_endpoint(time_s, cx)?;
        let groups = self.staged_heat.len();
        for (heat, &integral) in self.staged_heat.iter_mut().zip(&self.integrals[..groups]) {
            *heat = finite(integral)?;
            if *heat < 0.0 {
                return Err(invalid("original frame heat integral lost nonnegativity"));
            }
        }
        self.system
            .storage()
            .gradient(&self.state, &mut self.old_effort);
        let mut change = Dd::ZERO;
        for i in 0..self.state.len() {
            change = change
                + (Dd::from_f64(self.candidate[i]) - Dd::from_f64(self.state[i]))
                    * (Dd::from_f64(self.effort[i]) + Dd::from_f64(self.old_effort[i]))
                    * Dd::from_f64(0.5);
        }
        self.staged.energy_change_j = finite(change.to_f64())?;
        self.staged.applied_work_j = sum(self
            .source
            .free_dofs
            .iter()
            .enumerate()
            .map(|(i, &dof)| nodal_force[dof] * (self.candidate[i] - self.state[i])))?;
        self.staged.frame_dissipation_j = self.staged_heat[0];
        self.staged.attached_dissipation_j = sum(self.staged_heat[1..].iter().copied())?;
        self.staged.energy_balance_residual_j = sum([
            self.staged.energy_change_j,
            self.staged.frame_dissipation_j,
            self.staged.attached_dissipation_j,
            -self.staged.applied_work_j,
        ])?;
        let scale = self
            .report
            .total_energy_j
            .abs()
            .max(self.staged.total_energy_j.abs())
            .max(self.staged.energy_change_j.abs())
            .max(self.staged.applied_work_j.abs())
            .max(self.staged.frame_dissipation_j)
            .max(self.staged.attached_dissipation_j);
        if self.staged.energy_balance_residual_j.abs() > self.tolerance * scale {
            return Err(invalid(
                "complete native frame storage, physical heat and applied work do not close",
            ));
        }
        self.staged.maximum_relative_equation_residual = self.equation_residual(cx)?;
        if self.staged.maximum_relative_equation_residual > self.tolerance {
            return Err(FrameCoupledError::Residual {
                relative: self.staged.maximum_relative_equation_residual,
                tolerance: self.tolerance,
            });
        }
        poll(cx)
    }

    fn equation_residual(&self, cx: &Cx<'_>) -> Result<f64> {
        let nf = self.source.free_dofs.len();
        let effort = &self.integrals[self.heat.len()..];
        let mut maximum = 0.0_f64;
        for i in 0..nf {
            poll(cx)?;
            let displacement = self.candidate[i] - self.state[i];
            maximum = maximum.max(relative(
                displacement - effort[nf + i],
                self.candidate[i].abs() + self.state[i].abs() + effort[nf + i].abs(),
            )?);
            let delta_p = self.candidate[nf + i] - self.state[nf + i];
            let applied = self.duration_s * self.input[i];
            let mut residual =
                Dd::from_f64(delta_p) - Dd::from_f64(applied) + Dd::from_f64(effort[i]);
            let mut scale = self.candidate[nf + i].abs()
                + self.state[nf + i].abs()
                + applied.abs()
                + effort[i].abs();
            for k in 0..nf {
                let term = self.source.damping[i * nf + k] * effort[nf + k];
                residual = residual + Dd::from_f64(term);
                scale += term.abs();
            }
            for p in 0..self.source.ports {
                let load = sum((0..self.source.load_dimension)
                    .map(|k| self.source.load_g[k * self.source.ports + p] * effort[2 * nf + k]))?;
                let term = self.source.force_map[i * self.source.ports + p] * load;
                residual = residual + Dd::from_f64(term);
                scale += term.abs();
            }
            maximum = maximum.max(relative(residual.to_f64(), scale)?);
        }
        Ok(maximum)
    }

    pub(in crate::frame2d) fn staged_state(&self) -> &[f64] {
        &self.candidate
    }
    pub(in crate::frame2d) fn staged_report(&self) -> &FrameTransientReport {
        &self.staged
    }
    pub(in crate::frame2d) fn staged_heat(&self) -> &[f64] {
        &self.staged_heat
    }
    pub(in crate::frame2d) fn publish_staged(&mut self) {
        core::mem::swap(&mut self.state, &mut self.candidate);
        core::mem::swap(&mut self.report, &mut self.staged);
        core::mem::swap(&mut self.heat, &mut self.staged_heat);
    }
}

fn resistance_form(
    q: &[f64],
    r: &[f64],
    n: usize,
    range: Range<usize>,
    cx: &Cx<'_>,
) -> Result<Vec<f64>> {
    let mut rq = vec![0.0; n * n];
    for i in range.clone() {
        poll(cx)?;
        for k in 0..n {
            rq[i * n + k] = sum(range.clone().map(|j| r[i * n + j] * q[j * n + k]))?;
        }
    }
    let mut form = vec![0.0; n * n];
    for i in 0..n {
        poll(cx)?;
        for k in 0..n {
            form[i * n + k] = sum(range.clone().map(|j| q[j * n + i] * rq[j * n + k]))?;
        }
    }
    Ok(form)
}

pub(in crate::frame2d) fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "native frame transient arithmetic is unrepresentable",
        ))
    }
}
pub(in crate::frame2d) fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let value = values
        .into_iter()
        .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value));
    finite(value.to_f64())
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    let value = a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    });
    finite(value.to_f64())
}
fn relative(residual: f64, scale: f64) -> Result<f64> {
    finite(residual)?;
    finite(scale)?;
    Ok(if scale == 0.0 {
        if residual == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        residual.abs() / scale
    })
}
fn invalid(message: &'static str) -> FrameCoupledError {
    FrameCoupledError::Invalid(message)
}
fn step_error(error: PreparedStepError) -> FrameCoupledError {
    match error {
        PreparedStepError::Cancelled => FrameCoupledError::Cancelled,
        PreparedStepError::Solver(error) => FrameCoupledError::Physical(error),
    }
}
