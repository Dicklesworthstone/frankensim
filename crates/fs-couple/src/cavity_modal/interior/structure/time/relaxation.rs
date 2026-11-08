//! Exact held-input evolution of the full coupled cavity with Maxwell histories.
//!
//! The original coordinates are q=[u,psi], q_dot=[v,p]. The native mechanical
//! owner retains M, C, RELAXED K and every supplied strain map; its canonical
//! memory state represents sigma_dot=k*b^T*q_dot-sigma/tau. No fitted damping,
//! mass lumping, modal truncation or extra differential-equation solver is used.
//! Pressure is always psi_dot. Original interface work and material work are
//! integrated on the same trajectory, separately from positive relaxation heat.

use super::{StructuralInterior, StructuralInteriorTimeEnergy};
use fs_exec::Cx;
pub use fs_modal::relaxation_modes::MaxwellDampingBranch;
use fs_phs::PhsError;
use fs_phs::linear_exact::DenseExponentialBudget;

/// Bounds on original coordinates, histories and the complete exact cache.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationTimeLimits {
    pub maximum_coordinates: usize,
    pub maximum_branches: usize,
    /// Conservative adapter payload, including construction and one workspace.
    /// Native exponential caches have the separate `preparation` budget. Borrowed
    /// source arrays, other accepted states and allocator overhead are excluded.
    pub maximum_adapter_bytes: u128,
    /// Conservative multiply/add score for adapter maps and physical accounts.
    /// Native mechanical admission has its own state cap; native exponential
    /// work is covered by `preparation`.
    pub maximum_adapter_work: u128,
    pub preparation: DenseExponentialBudget,
    pub energy_balance_tolerance: f64,
}
impl Default for StructuralInteriorRelaxationTimeLimits {
    fn default() -> Self {
        let ordinary = super::StructuralInteriorTimeLimits::default();
        Self {
            maximum_coordinates: ordinary.maximum_coordinates,
            maximum_branches: 32,
            maximum_adapter_bytes: 64 * 1024 * 1024,
            maximum_adapter_work: 1_000_000_000,
            preparation: ordinary.preparation,
            energy_balance_tolerance: ordinary.energy_balance_tolerance,
        }
    }
}

/// Source, native dynamics, work-admission or cancellation refusal.
#[derive(Debug)]
pub enum StructuralInteriorRelaxationTimeError {
    Invalid(&'static str),
    Limit {
        quantity: &'static str,
        requested: u128,
        maximum: u128,
    },
    Dynamics(PhsError),
    Cancelled,
    EnergyBalance {
        account: &'static str,
        residual_j: f64,
        limit_j: f64,
    },
}
impl core::fmt::Display for StructuralInteriorRelaxationTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "structural interior relaxation time: {what}"),
            Self::Limit {
                quantity,
                requested,
                maximum,
            } => write!(
                f,
                "structural interior relaxation time {quantity}: {requested} exceeds {maximum}"
            ),
            Self::Dynamics(error) => error.fmt(f),
            Self::Cancelled => f.write_str("structural interior relaxation time cancelled"),
            Self::EnergyBalance {
                account,
                residual_j,
                limit_j,
            } => write!(
                f,
                "structural interior {account} work defect {residual_j} J exceeds {limit_j} J"
            ),
        }
    }
}
impl std::error::Error for StructuralInteriorRelaxationTimeError {}
impl From<PhsError> for StructuralInteriorRelaxationTimeError {
    fn from(error: PhsError) -> Self {
        Self::Dynamics(error)
    }
}

/// Spatial storage uses the original relaxed stiffness. Memory is additional.
#[derive(Clone, Copy, Debug, Default)]
pub struct StructuralInteriorRelaxationTimeEnergy {
    pub spatial: StructuralInteriorTimeEnergy,
    pub relaxation_j: f64,
}
impl StructuralInteriorRelaxationTimeEnergy {
    pub fn total_j(self) -> f64 {
        self.spatial.total_j() + self.relaxation_j
    }
}

/// Independently integrated work of one original branch, in caller source order.
#[derive(Clone, Copy, Debug, Default)]
pub struct StructuralInteriorRelaxationTimeMemoryWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    /// Bilinear change, avoiding subtraction of two large stored energies.
    pub energy_change_j: f64,
    pub structural_work_j: f64,
    pub acoustic_work_j: f64,
    /// Integral sigma^2/(k*tau), never inferred from an energy deficit.
    pub dissipation_j: f64,
}
impl StructuralInteriorRelaxationTimeMemoryWork {
    pub fn balance_residual_j(self) -> f64 {
        self.energy_change_j + self.dissipation_j - self.structural_work_j - self.acoustic_work_j
    }
}

/// Separate original source, reciprocal interface, storage and heat accounts.
#[derive(Clone, Copy, Debug, Default)]
pub struct StructuralInteriorRelaxationTimeWork {
    pub initial: StructuralInteriorRelaxationTimeEnergy,
    pub final_energy: StructuralInteriorRelaxationTimeEnergy,
    pub structural_energy_change_j: f64,
    pub acoustic_energy_change_j: f64,
    pub relaxation_energy_change_j: f64,
    pub structural_input_work_j: f64,
    pub acoustic_input_work_j: f64,
    pub structural_dissipation_j: f64,
    pub wall_dissipation_j: f64,
    pub relaxation_dissipation_j: f64,
    /// Positive from fluid into structure; explicit material work is separate.
    pub fluid_to_structure_work_j: f64,
    pub structural_to_relaxation_work_j: f64,
    pub acoustic_to_relaxation_work_j: f64,
    pub balance_tolerance_j: f64,
}
impl StructuralInteriorRelaxationTimeWork {
    pub fn structural_balance_residual_j(self) -> f64 {
        self.structural_energy_change_j
            + self.structural_dissipation_j
            + self.structural_to_relaxation_work_j
            - self.structural_input_work_j
            - self.fluid_to_structure_work_j
    }
    pub fn acoustic_balance_residual_j(self) -> f64 {
        self.acoustic_energy_change_j + self.wall_dissipation_j + self.acoustic_to_relaxation_work_j
            - self.acoustic_input_work_j
            + self.fluid_to_structure_work_j
    }
    pub fn relaxation_balance_residual_j(self) -> f64 {
        self.relaxation_energy_change_j + self.relaxation_dissipation_j
            - self.structural_to_relaxation_work_j
            - self.acoustic_to_relaxation_work_j
    }
    pub fn balance_residual_j(self) -> f64 {
        self.structural_energy_change_j
            + self.acoustic_energy_change_j
            + self.relaxation_energy_change_j
            + self.structural_dissipation_j
            + self.wall_dissipation_j
            + self.relaxation_dissipation_j
            - self.structural_input_work_j
            - self.acoustic_input_work_j
    }
}

#[derive(Clone, Debug)]
pub struct StructuralInteriorRelaxationTimeFrame {
    pub time_s: f64,
    pub duration_s: f64,
    pub work: StructuralInteriorRelaxationTimeWork,
    pub memory: Vec<StructuralInteriorRelaxationTimeMemoryWork>,
    /// Original B^T*p, excluding explicit branch reaction.
    pub fluid_force: Vec<f64>,
    /// Original B*v, excluding explicit branch contribution.
    pub displaced_flow_m3_s: Vec<f64>,
    pub integrated_structural_velocity: Vec<f64>,
    pub integrated_pressure_pa_s: Vec<f64>,
}

/// Original fields and material histories, bound to one physical source.
#[derive(Clone)]
pub struct StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    branches: &'m [MaxwellDampingBranch],
    inner: kernel::State,
}
impl StructuralInteriorRelaxationTimeState<'_, '_, '_, '_> {
    pub fn time_s(&self) -> f64 {
        self.inner.time_s
    }
    pub fn structural_displacement(&self) -> &[f64] {
        &self.inner.canonical[..self.model.structural_dofs]
    }
    pub fn structural_velocity(&self) -> &[f64] {
        &self.inner.velocity[..self.model.structural_dofs]
    }
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        &self.inner.canonical[self.model.structural_dofs..self.inner.velocity.len()]
    }
    pub fn pressure_pa(&self) -> &[f64] {
        &self.inner.velocity[self.model.structural_dofs..]
    }
    pub fn branch_stresses(&self) -> &[f64] {
        &self.inner.stress
    }
    pub fn branch_viscous_strains(&self) -> &[f64] {
        &self.inner.viscous_strain
    }
}

/// Fixed scratch for allocation-free candidate evaluation and physical checks.
pub struct StructuralInteriorRelaxationTimeWorkspace<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    branches: &'m [MaxwellDampingBranch],
    inner: kernel::Workspace,
}
pub struct StructuralInteriorRelaxationTimeStep<'m, 'a, 'f, 'c> {
    pub state: StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c>,
    pub frame: StructuralInteriorRelaxationTimeFrame,
}

/// One native whole-system propagation cache and original physical source.
pub struct PreparedStructuralInteriorRelaxationTime<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    branches: &'m [MaxwellDampingBranch],
    inner: kernel::Prepared<'m>,
}

impl<'a, 'f, 'c> StructuralInterior<'a, 'f, 'c> {
    /// Exact held generalized forces and projected inward-flow loads, with
    /// explicit Maxwell channels in original [u,psi] order. K must be RELAXED.
    /// Structural, potential-only and caller-declared mixed strain maps retain
    /// their units and source order. A zero strain map retains pure relaxation.
    ///
    /// Preparation preserves the complete 2*n+b state, including acoustic gauge
    /// coordinates and dark material modes. The native doubled augmented bound
    /// `2*(3*n+b+1)<=256` applies without truncation. Cancellation is polled
    /// around bounded native mechanical admission and throughout trajectory
    /// preparation. Native SPD/PSD/passivity admission remains authoritative.
    pub fn prepare_relaxation_exact_time<'m>(
        &'m self,
        branches: &'m [MaxwellDampingBranch],
        duration_s: f64,
        limits: StructuralInteriorRelaxationTimeLimits,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedStructuralInteriorRelaxationTime<'m, 'a, 'f, 'c>,
        StructuralInteriorRelaxationTimeError,
    > {
        let mut cancelled = || cx.checkpoint().is_err();
        let source = kernel::Source {
            structural: self.structural_dofs,
            coordinates: self.structural_dofs + self.pressure_dofs(),
            mass: &self.mass,
            damping: &self.damping,
            stiffness: &self.stiffness,
            coupling: &self.coupling,
        };
        kernel::admit(
            source.coordinates,
            branches.len(),
            duration_s,
            limits,
            &mut cancelled,
        )?;
        let mut native = Vec::with_capacity(branches.len());
        for branch in branches {
            kernel::poll(&mut cancelled)?;
            if branch.strain.len() != source.coordinates {
                return Err(kernel::invalid(
                    "every branch needs the complete original strain map",
                ));
            }
            native.push(fs_phs::RelaxationBranch {
                projection: branch.strain.clone(),
                stiffness: branch.stiffness,
                relaxation_time_s: branch.relaxation_time_s,
            });
        }
        let inner = kernel::Prepared::new(source, native, duration_s, limits, &mut cancelled)?;
        Ok(PreparedStructuralInteriorRelaxationTime {
            model: self,
            branches,
            inner,
        })
    }
}

impl<'m, 'a, 'f, 'c> PreparedStructuralInteriorRelaxationTime<'m, 'a, 'f, 'c> {
    pub fn duration_s(&self) -> f64 {
        self.inner.duration_s
    }

    /// Initialize the ORIGINAL stresses, not an assumed relaxed material.
    /// The conversion to native viscous memory is z=sqrt(k)*b^T*q-sigma/sqrt(k).
    #[allow(clippy::too_many_arguments)]
    pub fn state(
        &self,
        time_s: f64,
        displacement: &[f64],
        velocity: &[f64],
        potential_pa_s: &[f64],
        pressure_pa: &[f64],
        stresses: &[f64],
        cx: &Cx<'_>,
    ) -> Result<
        StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c>,
        StructuralInteriorRelaxationTimeError,
    > {
        if displacement.len() != self.model.structural_dofs
            || velocity.len() != self.model.structural_dofs
            || potential_pa_s.len() != self.model.pressure_dofs()
            || pressure_pa.len() != self.model.pressure_dofs()
        {
            return Err(kernel::invalid(
                "complete original structural and acoustic state required",
            ));
        }
        let q: Vec<_> = displacement.iter().chain(potential_pa_s).copied().collect();
        let v: Vec<_> = velocity.iter().chain(pressure_pa).copied().collect();
        let inner = self
            .inner
            .state(time_s, &q, &v, stresses, &mut || cx.checkpoint().is_err())?;
        Ok(StructuralInteriorRelaxationTimeState {
            model: self.model,
            branches: self.branches,
            inner,
        })
    }

    pub fn workspace(&self) -> StructuralInteriorRelaxationTimeWorkspace<'m, 'a, 'f, 'c> {
        StructuralInteriorRelaxationTimeWorkspace {
            model: self.model,
            branches: self.branches,
            inner: self.inner.workspace(),
        }
    }
    /// Explicitly initialize every branch at zero stress on the supplied fields.
    pub fn relaxed_state(
        &self,
        time_s: f64,
        displacement: &[f64],
        velocity: &[f64],
        potential_pa_s: &[f64],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<
        StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c>,
        StructuralInteriorRelaxationTimeError,
    > {
        self.state(
            time_s,
            displacement,
            velocity,
            potential_pa_s,
            pressure_pa,
            &vec![0.0; self.branches.len()],
            cx,
        )
    }
    pub fn energy(
        &self,
        state: &StructuralInteriorRelaxationTimeState<'_, '_, '_, '_>,
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorRelaxationTimeEnergy, StructuralInteriorRelaxationTimeError> {
        self.validate(state.model, state.branches)?;
        self.inner
            .energy(&state.inner, &mut || cx.checkpoint().is_err())
            .map(|v| v.0)
    }
    pub fn step(
        &self,
        initial: &StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        inward_flow_load: &[f64],
        cx: &Cx<'_>,
    ) -> Result<
        StructuralInteriorRelaxationTimeStep<'m, 'a, 'f, 'c>,
        StructuralInteriorRelaxationTimeError,
    > {
        self.validate(initial.model, initial.branches)?;
        let mut state = initial.clone();
        let mut workspace = self.workspace();
        self.step_into(
            &mut state,
            structural_force,
            inward_flow_load,
            &mut workspace,
            cx,
        )?;
        Ok(StructuralInteriorRelaxationTimeStep {
            state,
            frame: workspace.inner.frame,
        })
    }
    /// Cached contractions only: no allocation, mass solve or exponential.
    /// State/history/clock publish together after every independent balance and
    /// the final cancellation poll. Refusal may alter private scratch only.
    pub fn step_into<'w>(
        &self,
        state: &mut StructuralInteriorRelaxationTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        inward_flow_load: &[f64],
        workspace: &'w mut StructuralInteriorRelaxationTimeWorkspace<'m, 'a, 'f, 'c>,
        cx: &Cx<'_>,
    ) -> Result<&'w StructuralInteriorRelaxationTimeFrame, StructuralInteriorRelaxationTimeError>
    {
        self.validate(state.model, state.branches)?;
        self.validate(workspace.model, workspace.branches)?;
        self.inner.step_into(
            &mut state.inner,
            structural_force,
            inward_flow_load,
            &mut workspace.inner,
            &mut || cx.checkpoint().is_err(),
        )
    }
    fn validate(
        &self,
        model: &StructuralInterior<'_, '_, '_>,
        branches: &[MaxwellDampingBranch],
    ) -> Result<(), StructuralInteriorRelaxationTimeError> {
        if !core::ptr::eq(self.model, model) || !core::ptr::eq(self.branches, branches) {
            return Err(kernel::invalid(
                "state/workspace belongs to another original source or branch array",
            ));
        }
        Ok(())
    }
}

// Pure source-to-native mapping and physical accounting. Kept independent of
// FE assembly and executor state so the actual runtime can be checked directly
// against original-coordinate analytic trajectories.
mod kernel {
    use super::{
        StructuralInteriorRelaxationTimeEnergy as Energy,
        StructuralInteriorRelaxationTimeError as Error,
        StructuralInteriorRelaxationTimeFrame as Frame,
        StructuralInteriorRelaxationTimeLimits as Limits,
        StructuralInteriorRelaxationTimeMemoryWork as MemoryWork,
        StructuralInteriorRelaxationTimeWork as Work, StructuralInteriorTimeEnergy,
    };
    use fs_math::{dd::Dd, det};
    use fs_phs::linear_exact::{
        LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, prepared::PreparedQuadraticFlow,
        split_quadratic,
    };
    use fs_phs::{PreparedStepError, RelaxationBranch, mechanical::LinearMechanicalSystem};

    #[derive(Clone, Copy)]
    pub(super) struct Source<'s> {
        pub structural: usize,
        pub coordinates: usize,
        pub mass: &'s [f64],
        pub damping: &'s [f64],
        pub stiffness: &'s [f64],
        pub coupling: &'s [f64],
    }
    #[derive(Clone)]
    pub(super) struct State {
        pub time_s: f64,
        pub canonical: Vec<f64>,
        pub velocity: Vec<f64>,
        pub stress: Vec<f64>,
        pub viscous_strain: Vec<f64>,
        elastic: Vec<f64>,
    }
    pub(super) struct Workspace {
        candidate: State,
        input: Vec<f64>,
        augmented: Vec<f64>,
        balanced: Vec<f64>,
        integrals: Vec<f64>,
        pub frame: Frame,
    }
    pub(super) struct Prepared<'s> {
        source: Source<'s>,
        mechanical: LinearMechanicalSystem,
        flow: PreparedQuadraticFlow,
        velocity_map: Vec<f64>,
        elastic_map: Vec<f64>,
        roots: Vec<f64>,
        pub duration_s: f64,
        limits: Limits,
    }

    pub(super) fn invalid(what: &'static str) -> Error {
        Error::Invalid(what)
    }
    pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
        if cancelled() {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
    pub(super) fn admit(
        n: usize,
        b: usize,
        dt: f64,
        limits: Limits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if n == 0
            || !dt.is_finite()
            || dt <= 0.0
            || !limits.energy_balance_tolerance.is_finite()
            || limits.energy_balance_tolerance <= 0.0
            || limits.energy_balance_tolerance >= 1.0
        {
            return Err(invalid(
                "positive finite interval, nonempty coordinates and work tolerance in (0,1) required",
            ));
        }
        let (n, b) = (n as u128, b as u128);
        for (quantity, requested, maximum) in [
            ("coordinates", n, limits.maximum_coordinates as u128),
            ("branches", b, limits.maximum_branches as u128),
            (
                "doubled augmented dimension",
                2 * (3 * n + b + 1),
                MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128,
            ),
        ] {
            if requested > maximum {
                return Err(Error::Limit {
                    quantity,
                    requested,
                    maximum,
                });
            }
        }
        // Dimension admission precedes cubic scores, so even hostile extents
        // cannot overflow a budget calculation before they are refused.
        let d = 2 * n + b;
        let forms = 3 + 3 * b;
        let bytes =
            8 * ((forms + 24) * d * d + (8 * n + 4 * b + 32) * d + 16 * n * n + 32 * b + 32 * n);
        let work = 64 * d * d * d + (8 * forms + 32) * d * d;
        for (quantity, requested, maximum) in [
            ("adapter bytes", bytes, limits.maximum_adapter_bytes),
            ("adapter work", work, limits.maximum_adapter_work),
        ] {
            if requested > maximum {
                return Err(Error::Limit {
                    quantity,
                    requested,
                    maximum,
                });
            }
        }
        Ok(())
    }

    impl<'s> Prepared<'s> {
        pub(super) fn new(
            source: Source<'s>,
            branches: Vec<RelaxationBranch>,
            duration_s: f64,
            limits: Limits,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<Self, Error> {
            let (n, ns, b) = (source.coordinates, source.structural, branches.len());
            admit(n, b, duration_s, limits, cancelled)?;
            if ns == 0
                || ns >= n
                || source.coupling.len() != ns * (n - ns)
                || source.coupling.iter().any(|v| !v.is_finite())
            {
                return Err(invalid("complete original reciprocal interface required"));
            }
            let mut identity = vec![0.0; n * n];
            for i in 0..n {
                identity[i * n + i] = 1.0;
            }
            let mechanical = LinearMechanicalSystem::new(
                source.mass,
                source.damping,
                source.stiffness,
                &identity,
                n,
                n,
                branches,
            )?;
            poll(cancelled)?;
            let d = 2 * n + b;
            let split = split_quadratic(
                mechanical.system(),
                LinearSplitBudget {
                    max_modes: d,
                    max_omega_dt: 1.0,
                },
            )?;
            let velocity_map = split.quadratic_matrix()[n * d..2 * n * d].to_vec();
            let mut elastic_map = vec![0.0; b * d];
            let mut roots = Vec::with_capacity(b);
            for (j, branch) in mechanical.branches().iter().enumerate() {
                poll(cancelled)?;
                let root = det::sqrt(branch.stiffness);
                roots.push(root);
                for i in 0..n {
                    elastic_map[j * d + i] = root * branch.projection[i];
                }
                elastic_map[j * d + 2 * n + j] = -1.0;
            }
            let mut forms = Vec::with_capacity(3 + 3 * b);
            let mut physical = vec![0.0; n * n];
            for (start, end) in [(0, ns), (ns, n)] {
                physical.fill(0.0);
                for i in start..end {
                    for j in start..end {
                        physical[i * n + j] =
                            source.damping[i * n + j].midpoint(source.damping[j * n + i]);
                    }
                }
                forms.push(pullback(&physical, &velocity_map, n, d, cancelled)?);
            }
            physical.fill(0.0);
            for a in 0..n - ns {
                for s in 0..ns {
                    let value = 0.5 * source.coupling[a * ns + s];
                    physical[s * n + ns + a] = value;
                    physical[(ns + a) * n + s] = value;
                }
            }
            forms.push(pullback(&physical, &velocity_map, n, d, cancelled)?);
            for (j, branch) in mechanical.branches().iter().enumerate() {
                poll(cancelled)?;
                let elastic = &elastic_map[j * d..(j + 1) * d];
                let loss: Vec<_> = elastic
                    .iter()
                    .map(|v| v / det::sqrt(branch.relaxation_time_s))
                    .collect();
                forms.push(product_form(&loss, &loss)?);
                let stress: Vec<_> = elastic.iter().map(|v| v * roots[j]).collect();
                for (start, end) in [(0, ns), (ns, n)] {
                    let mut rate = vec![0.0; d];
                    for k in 0..d {
                        let mut sum = Dd::from_f64(0.0);
                        for i in start..end {
                            sum = sum
                                + Dd::from_f64(branch.projection[i])
                                    * Dd::from_f64(velocity_map[i * d + k]);
                        }
                        rate[k] = sum.to_f64();
                    }
                    forms.push(product_form(&stress, &rate)?);
                }
            }
            let linear: Vec<_> = velocity_map.chunks_exact(d).map(|v| v.to_vec()).collect();
            let flow = PreparedQuadraticFlow::new_controlled(
                mechanical.system(),
                &split,
                duration_s,
                &forms,
                &linear,
                limits.preparation,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                PreparedStepError::Cancelled => Error::Cancelled,
                PreparedStepError::Solver(error) => Error::Dynamics(error),
            })?;
            poll(cancelled)?;
            Ok(Self {
                source,
                mechanical,
                flow,
                velocity_map,
                elastic_map,
                roots,
                duration_s,
                limits,
            })
        }

        fn blank(&self) -> State {
            let (n, b) = (self.source.coordinates, self.roots.len());
            State {
                time_s: 0.0,
                canonical: vec![0.0; 2 * n + b],
                velocity: vec![0.0; n],
                stress: vec![0.0; b],
                viscous_strain: vec![0.0; b],
                elastic: vec![0.0; b],
            }
        }
        pub(super) fn state(
            &self,
            time_s: f64,
            q: &[f64],
            velocity: &[f64],
            stresses: &[f64],
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<State, Error> {
            poll(cancelled)?;
            let (n, b) = (self.source.coordinates, self.roots.len());
            if !time_s.is_finite() || stresses.len() != b || stresses.iter().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "finite physical time and every original branch stress required",
                ));
            }
            let mut state = self.blank();
            state.canonical = self.mechanical.state(q, velocity, &vec![0.0; b])?;
            state.time_s = time_s;
            let d = state.canonical.len();
            for j in 0..b {
                poll(cancelled)?;
                let mut strain = Dd::from_f64(0.0);
                for (coefficient, coordinate) in self.elastic_map[j * d..j * d + n].iter().zip(q) {
                    strain = strain + Dd::from_f64(*coefficient) * Dd::from_f64(*coordinate);
                }
                state.canonical[2 * n + j] =
                    (strain - Dd::from_f64(stresses[j] / self.roots[j])).to_f64();
            }
            self.recover(&mut state, cancelled)?;
            for (j, (&expected, &actual)) in stresses.iter().zip(&state.stress).enumerate() {
                // A requested nonzero history must survive the native chart.
                // An explicitly relaxed arm can retain only roundoff from
                // representing its projected equilibrium strain in one word.
                let scale = if expected == 0.0 {
                    let (_, absolute) = work_dot(&self.elastic_map[j * d..j * d + n], q)?;
                    finite(self.roots[j] * (absolute + state.canonical[2 * n + j].abs()))?
                } else {
                    actual.abs().max(expected.abs())
                };
                if (actual - expected).abs() > 512.0 * f64::EPSILON * d as f64 * scale {
                    return Err(invalid(
                        "supplied initial stress cannot be retained in the native memory coordinates",
                    ));
                }
            }
            self.energy(&state, cancelled)?;
            Ok(state)
        }
        fn recover(
            &self,
            state: &mut State,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<(), Error> {
            let n = self.source.coordinates;
            let d = state.canonical.len();
            if state.canonical.iter().any(|v| !v.is_finite()) {
                return Err(invalid("finite complete native state required"));
            }
            for i in 0..n {
                poll(cancelled)?;
                state.velocity[i] = dot(&self.velocity_map[i * d..(i + 1) * d], &state.canonical)?;
            }
            for j in 0..self.roots.len() {
                poll(cancelled)?;
                state.elastic[j] = dot(&self.elastic_map[j * d..(j + 1) * d], &state.canonical)?;
                state.stress[j] = self.roots[j] * state.elastic[j];
                state.viscous_strain[j] = state.canonical[2 * n + j] / self.roots[j];
                finite(state.stress[j])?;
                finite(state.viscous_strain[j])?;
            }
            Ok(())
        }
        pub(super) fn workspace(&self) -> Workspace {
            let (n, ns, b) = (
                self.source.coordinates,
                self.source.structural,
                self.roots.len(),
            );
            Workspace {
                candidate: self.blank(),
                input: vec![0.0; n],
                augmented: vec![0.0; 3 * n + b + 1],
                balanced: vec![0.0; 2 * n + b],
                integrals: vec![0.0; 3 + 3 * b + n],
                frame: Frame {
                    time_s: 0.0,
                    duration_s: self.duration_s,
                    work: Work::default(),
                    memory: vec![MemoryWork::default(); b],
                    fluid_force: vec![0.0; ns],
                    displaced_flow_m3_s: vec![0.0; n - ns],
                    integrated_structural_velocity: vec![0.0; ns],
                    integrated_pressure_pa_s: vec![0.0; n - ns],
                },
            }
        }
        pub(super) fn energy(
            &self,
            state: &State,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<(Energy, [f64; 3]), Error> {
            let (n, ns) = (self.source.coordinates, self.source.structural);
            let (sk, rsk) = quadratic(self.source.mass, n, 0, ns, &state.velocity, cancelled)?;
            let (ss, rss) = quadratic(
                self.source.stiffness,
                n,
                0,
                ns,
                &state.canonical[..n],
                cancelled,
            )?;
            let (fk, rfk) = quadratic(
                self.source.stiffness,
                n,
                ns,
                n,
                &state.canonical[..n],
                cancelled,
            )?;
            let (fc, rfc) = quadratic(self.source.mass, n, ns, n, &state.velocity, cancelled)?;
            let mut memory = Dd::from_f64(0.0);
            for &e in &state.elastic {
                memory = memory + Dd::from_f64(0.5 * e) * Dd::from_f64(e);
            }
            let memory = finite(memory.to_f64())?;
            Ok((
                Energy {
                    spatial: StructuralInteriorTimeEnergy {
                        structural_kinetic_j: sk,
                        structural_strain_j: ss,
                        fluid_kinetic_j: fk,
                        fluid_compression_j: fc,
                    },
                    relaxation_j: memory,
                },
                [
                    rsk + rss,
                    rfk + rfc,
                    128.0 * f64::EPSILON * self.roots.len() as f64 * memory.abs(),
                ],
            ))
        }

        pub(super) fn step_into<'w>(
            &self,
            state: &mut State,
            force: &[f64],
            inward: &[f64],
            workspace: &'w mut Workspace,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<&'w Frame, Error> {
            poll(cancelled)?;
            let (n, ns, b) = (
                self.source.coordinates,
                self.source.structural,
                self.roots.len(),
            );
            if force.len() != ns
                || inward.len() != n - ns
                || force.iter().chain(inward).any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "complete finite original held forces and inward-flow loads required",
                ));
            }
            let next_time = state.time_s + self.duration_s;
            if !next_time.is_finite() || next_time <= state.time_s {
                return Err(invalid("interval must advance the finite physical clock"));
            }
            workspace.input[..ns].copy_from_slice(force);
            workspace.input[ns..].copy_from_slice(inward);
            let (initial, r0) = self.energy(state, cancelled)?;
            self.flow.evaluate_into(
                &state.canonical,
                &workspace.input,
                &mut workspace.augmented,
                &mut workspace.balanced,
                &mut workspace.candidate.canonical,
                &mut workspace.integrals,
            )?;
            self.recover(&mut workspace.candidate, cancelled)?;
            workspace.candidate.time_s = next_time;
            let (final_energy, r1) = self.energy(&workspace.candidate, cancelled)?;
            let (ds, rs) = self.spatial_change(state, &workspace.candidate, 0, ns, cancelled)?;
            let (da, ra) = self.spatial_change(state, &workspace.candidate, ns, n, cancelled)?;
            let offset = 3 + 3 * b;
            let (ws, aws) = work_dot(force, &workspace.integrals[offset..offset + ns])?;
            let (wa, awa) = work_dot(inward, &workspace.integrals[offset + ns..])?;
            let mut work = Work {
                initial,
                final_energy,
                structural_energy_change_j: ds,
                acoustic_energy_change_j: da,
                structural_input_work_j: ws,
                acoustic_input_work_j: wa,
                structural_dissipation_j: workspace.integrals[0],
                wall_dissipation_j: workspace.integrals[1],
                fluid_to_structure_work_j: workspace.integrals[2],
                ..Work::default()
            };
            let mut memory_change = Dd::from_f64(0.0);
            let mut memory_heat = Dd::from_f64(0.0);
            let mut branch_structural = Dd::from_f64(0.0);
            let mut branch_acoustic = Dd::from_f64(0.0);
            let mut absolute_structural = 0.0;
            let mut absolute_acoustic = 0.0;
            let mut memory_roundoff = 0.0;
            for j in 0..b {
                poll(cancelled)?;
                let (a, z) = (state.elastic[j], workspace.candidate.elastic[j]);
                let change = (Dd::from_f64(0.5)
                    * (Dd::from_f64(z) - Dd::from_f64(a))
                    * (Dd::from_f64(z) + Dd::from_f64(a)))
                .to_f64();
                let item = MemoryWork {
                    initial_energy_j: 0.5 * a * a,
                    final_energy_j: 0.5 * z * z,
                    energy_change_j: change,
                    dissipation_j: workspace.integrals[3 + 3 * j],
                    structural_work_j: workspace.integrals[4 + 3 * j],
                    acoustic_work_j: workspace.integrals[5 + 3 * j],
                };
                let roundoff = 128.0
                    * f64::EPSILON
                    * (item.initial_energy_j.abs() + item.final_energy_j.abs() + change.abs());
                let limit = self.limit(
                    &[
                        item.initial_energy_j,
                        item.final_energy_j,
                        change,
                        item.dissipation_j,
                        item.structural_work_j,
                        item.acoustic_work_j,
                    ],
                    roundoff,
                );
                balance("material memory", item.balance_residual_j(), limit)?;
                nonnegative("relaxation heat", item.dissipation_j, limit)?;
                workspace.frame.memory[j] = item;
                work.balance_tolerance_j = work.balance_tolerance_j.max(limit);
                memory_change = memory_change + Dd::from_f64(change);
                memory_heat = memory_heat + Dd::from_f64(item.dissipation_j);
                branch_structural = branch_structural + Dd::from_f64(item.structural_work_j);
                branch_acoustic = branch_acoustic + Dd::from_f64(item.acoustic_work_j);
                absolute_structural += item.structural_work_j.abs();
                absolute_acoustic += item.acoustic_work_j.abs();
                memory_roundoff += roundoff;
            }
            work.relaxation_energy_change_j = finite(memory_change.to_f64())?;
            work.relaxation_dissipation_j = finite(memory_heat.to_f64())?;
            work.structural_to_relaxation_work_j = finite(branch_structural.to_f64())?;
            work.acoustic_to_relaxation_work_j = finite(branch_acoustic.to_f64())?;
            let ls = self.limit(
                &[
                    initial.spatial.structural_j(),
                    final_energy.spatial.structural_j(),
                    ds,
                    aws,
                    work.structural_dissipation_j,
                    work.fluid_to_structure_work_j,
                    absolute_structural,
                ],
                r0[0] + r1[0] + rs,
            );
            let la = self.limit(
                &[
                    initial.spatial.acoustic_j(),
                    final_energy.spatial.acoustic_j(),
                    da,
                    awa,
                    work.wall_dissipation_j,
                    work.fluid_to_structure_work_j,
                    absolute_acoustic,
                ],
                r0[1] + r1[1] + ra,
            );
            let lm = self.limit(
                &[
                    initial.relaxation_j,
                    final_energy.relaxation_j,
                    work.relaxation_energy_change_j,
                    work.relaxation_dissipation_j,
                    absolute_structural,
                    absolute_acoustic,
                ],
                r0[2] + r1[2] + memory_roundoff,
            );
            let lt = self.limit(
                &[
                    initial.total_j(),
                    final_energy.total_j(),
                    ds,
                    da,
                    work.relaxation_energy_change_j,
                    aws + awa,
                    work.structural_dissipation_j,
                    work.wall_dissipation_j,
                    work.relaxation_dissipation_j,
                ],
                r0.iter().sum::<f64>() + r1.iter().sum::<f64>() + rs + ra + memory_roundoff,
            );
            for (name, residual, limit) in [
                ("structure", work.structural_balance_residual_j(), ls),
                ("acoustic field", work.acoustic_balance_residual_j(), la),
                (
                    "all material memory",
                    work.relaxation_balance_residual_j(),
                    lm,
                ),
                ("complete system", work.balance_residual_j(), lt),
            ] {
                balance(name, residual, limit)?;
            }
            nonnegative("structural heat", work.structural_dissipation_j, ls)?;
            nonnegative("wall heat", work.wall_dissipation_j, la)?;
            work.balance_tolerance_j = work.balance_tolerance_j.max(ls).max(la).max(lm).max(lt);
            workspace.frame.fluid_force.fill(0.0);
            workspace.frame.displaced_flow_m3_s.fill(0.0);
            for s in 0..ns {
                poll(cancelled)?;
                let mut value = Dd::from_f64(0.0);
                for a in 0..n - ns {
                    value = value
                        + Dd::from_f64(self.source.coupling[a * ns + s])
                            * Dd::from_f64(workspace.candidate.velocity[ns + a]);
                }
                workspace.frame.fluid_force[s] = finite(value.to_f64())?;
            }
            for a in 0..n - ns {
                workspace.frame.displaced_flow_m3_s[a] = dot(
                    &self.source.coupling[a * ns..(a + 1) * ns],
                    &workspace.candidate.velocity[..ns],
                )?;
            }
            workspace
                .frame
                .integrated_structural_velocity
                .copy_from_slice(&workspace.integrals[offset..offset + ns]);
            workspace
                .frame
                .integrated_pressure_pa_s
                .copy_from_slice(&workspace.integrals[offset + ns..]);
            workspace.frame.time_s = next_time;
            workspace.frame.duration_s = self.duration_s;
            workspace.frame.work = work;
            poll(cancelled)?;
            core::mem::swap(state, &mut workspace.candidate);
            Ok(&workspace.frame)
        }

        fn spatial_change(
            &self,
            a: &State,
            z: &State,
            start: usize,
            end: usize,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<(f64, f64), Error> {
            let n = self.source.coordinates;
            let mut sum = Dd::from_f64(0.0);
            let mut absolute = 0.0;
            for (matrix, before, after) in [
                (
                    self.source.mass,
                    a.velocity.as_slice(),
                    z.velocity.as_slice(),
                ),
                (self.source.stiffness, &a.canonical[..n], &z.canonical[..n]),
            ] {
                for i in start..end {
                    poll(cancelled)?;
                    for j in start..end {
                        let term = Dd::from_f64(0.5 * matrix[i * n + j])
                            * (Dd::from_f64(after[i]) - Dd::from_f64(before[i]))
                            * (Dd::from_f64(after[j]) + Dd::from_f64(before[j]));
                        sum = sum + term;
                        absolute += term.to_f64().abs();
                    }
                }
            }
            Ok((
                finite(sum.to_f64())?,
                finite(128.0 * f64::EPSILON * n as f64 * absolute)?,
            ))
        }
        fn limit(&self, values: &[f64], roundoff: f64) -> f64 {
            let scale = values.iter().map(|v| v.abs()).fold(0.0, f64::max);
            self.limits.energy_balance_tolerance * scale
                + roundoff
                + 128.0
                    * f64::EPSILON
                    * (2 * self.source.coordinates + self.roots.len()) as f64
                    * values.iter().map(|v| v.abs()).sum::<f64>()
        }
    }

    fn dot(a: &[f64], b: &[f64]) -> Result<f64, Error> {
        work_dot(a, b).map(|v| v.0)
    }
    fn work_dot(a: &[f64], b: &[f64]) -> Result<(f64, f64), Error> {
        let mut sum = Dd::from_f64(0.0);
        let mut absolute = 0.0;
        for (&a, &b) in a.iter().zip(b) {
            let term = Dd::from_f64(a) * Dd::from_f64(b);
            sum = sum + term;
            absolute += term.to_f64().abs();
        }
        Ok((finite(sum.to_f64())?, finite(absolute)?))
    }
    fn finite(value: f64) -> Result<f64, Error> {
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid(
                "original physical state or account is unrepresentable",
            ))
        }
    }
    fn balance(account: &'static str, residual: f64, limit: f64) -> Result<(), Error> {
        finite(residual)?;
        finite(limit)?;
        if residual.abs() > limit {
            Err(Error::EnergyBalance {
                account,
                residual_j: residual,
                limit_j: limit,
            })
        } else {
            Ok(())
        }
    }
    fn nonnegative(account: &'static str, value: f64, limit: f64) -> Result<(), Error> {
        finite(value)?;
        if value < -limit {
            Err(Error::EnergyBalance {
                account,
                residual_j: value,
                limit_j: limit,
            })
        } else {
            Ok(())
        }
    }
    fn quadratic(
        matrix: &[f64],
        n: usize,
        start: usize,
        end: usize,
        x: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), Error> {
        let mut sum = Dd::from_f64(0.0);
        let mut absolute = 0.0;
        for i in start..end {
            poll(cancelled)?;
            for j in start..end {
                let term =
                    Dd::from_f64(0.5 * matrix[i * n + j]) * Dd::from_f64(x[i]) * Dd::from_f64(x[j]);
                sum = sum + term;
                absolute += term.to_f64().abs();
            }
        }
        let (value, roundoff) = (
            finite(sum.to_f64())?,
            finite(128.0 * f64::EPSILON * n as f64 * absolute)?,
        );
        nonnegative("spatial storage", value, roundoff)?;
        Ok((value, roundoff))
    }
    fn product_form(a: &[f64], b: &[f64]) -> Result<Vec<f64>, Error> {
        let d = a.len();
        let mut form = vec![0.0; d * d];
        for i in 0..d {
            for j in i..d {
                let value = (Dd::from_f64(0.5 * a[i]) * Dd::from_f64(b[j])
                    + Dd::from_f64(0.5 * b[i]) * Dd::from_f64(a[j]))
                .to_f64();
                form[i * d + j] = finite(value)?;
                form[j * d + i] = value;
            }
        }
        Ok(form)
    }
    fn pullback(
        physical: &[f64],
        velocity: &[f64],
        n: usize,
        d: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, Error> {
        let mut right = vec![0.0; n * d];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..d {
                let mut sum = Dd::from_f64(0.0);
                for k in 0..n {
                    sum =
                        sum + Dd::from_f64(physical[i * n + k]) * Dd::from_f64(velocity[k * d + j]);
                }
                right[i * d + j] = finite(sum.to_f64())?;
            }
        }
        let mut form = vec![0.0; d * d];
        for i in 0..d {
            poll(cancelled)?;
            for j in i..d {
                let mut sum = Dd::from_f64(0.0);
                for k in 0..n {
                    sum = sum + Dd::from_f64(velocity[k * d + i]) * Dd::from_f64(right[k * d + j]);
                }
                let value = finite(sum.to_f64())?;
                form[i * d + j] = value;
                form[j * d + i] = value;
            }
        }
        Ok(form)
    }
}

#[cfg(test)]
#[path = "relaxation_tests.rs"]
mod tests;
