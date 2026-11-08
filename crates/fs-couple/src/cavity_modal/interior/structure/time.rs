//! Exact held-input evolution of the original structural/interior pencil.
//!
//! The state contains every original structural displacement/velocity and every
//! pressure-potential/pressure coefficient. A numerical diagonal coordinate
//! change and the original complete mass inverse give canonical momenta; the
//! physical matrices are unchanged. Zero-stiffness structural coordinates and
//! the constant fluid potential remain in the state without a small frequency.
//!
//! The existing pHS owner advances the full system and independently integrates
//! structural loss, wall loss and SIGNED fluid-to-structure work on that same
//! trajectory. Held-force work uses integrated structural velocity; held inward
//! flow work uses integrated pressure. Interface work is a quadratic trajectory
//! integral, not a product of averaged port values or an energy deficit.
//!
//! Exactness concerns this frozen semidiscrete linear model under held inputs,
//! up to numerical exponential error. Spatial convergence, changing materials or
//! geometry, source bandwidth and output antialiasing remain caller obligations.
//! Dense preparation is bounded; oversized original spaces refuse without
//! truncation. Warm `step_into` uses fixed scratch and cached contractions, with
//! no factorization, exponential, material query or allocation.

use super::StructuralInterior;
use fs_exec::Cx;
use fs_la::factor::{FactorError, FactorRunError, cholesky_cancellable};
use fs_math::{dd::Dd, det};
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    prepared::PreparedQuadraticFlow, split_quadratic,
};
use fs_phs::{PhsError, PortHamiltonian, PreparedStepError, QuadraticStorage};

/// Complete original coordinates with explicit causal Maxwell histories.
pub mod relaxation;

/// Admission of all original coordinates and the native exact trajectory cache.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorTimeLimits {
    /// Structural plus pressure coordinates, each retaining two state variables.
    /// The native doubled augmented-state cap also applies: with one held input
    /// per coordinate, `2 * (3*n + 1) <= 256`. No coordinate is discarded.
    pub maximum_coordinates: usize,
    /// Native exponential/Gramian cache, memory and multiply-work budget.
    /// Original model and O(n^2) coordinate/observation maps are additional.
    pub preparation: DenseExponentialBudget,
    /// Relative physical work defect, plus evaluated contraction roundoff.
    pub energy_balance_tolerance: f64,
}

impl Default for StructuralInteriorTimeLimits {
    fn default() -> Self {
        Self {
            maximum_coordinates: 42,
            preparation: DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: 64 * 1024 * 1024,
                max_multiply_terms: 32_000_000_000,
            },
            energy_balance_tolerance: 5e-10,
        }
    }
}

/// Original factorization, propagation and physical-work refusals.
#[derive(Debug)]
pub enum StructuralInteriorTimeError {
    Invalid(&'static str),
    Limit { requested: usize, maximum: usize },
    Mass(FactorError),
    Dynamics(PhsError),
    Cancelled,
    EnergyBalance {
        account: &'static str,
        residual_j: f64,
        limit_j: f64,
    },
}

impl core::fmt::Display for StructuralInteriorTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "structural interior time: {what}"),
            Self::Limit { requested, maximum } => write!(f,
                "structural interior time requires {requested} coordinates, maximum {maximum}"),
            Self::Mass(error) => error.fmt(f),
            Self::Dynamics(error) => error.fmt(f),
            Self::Cancelled => f.write_str("structural interior time cancelled"),
            Self::EnergyBalance { account, residual_j, limit_j } => write!(f,
                "structural interior {account} work defect {residual_j} J exceeds {limit_j} J"),
        }
    }
}
impl std::error::Error for StructuralInteriorTimeError {}
impl From<PhsError> for StructuralInteriorTimeError {
    fn from(error: PhsError) -> Self { Self::Dynamics(error) }
}

/// Complete accepted physical coordinates and one accepted clock.
#[derive(Clone)]
pub struct StructuralInteriorTimeState<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    time_s: f64,
    displacement: Vec<f64>,
    velocity: Vec<f64>,
}

impl StructuralInteriorTimeState<'_, '_, '_, '_> {
    #[must_use]
    pub const fn time_s(&self) -> f64 { self.time_s }
    #[must_use]
    pub fn structural_displacement(&self) -> &[f64] {
        &self.displacement[..self.model.structural_dofs]
    }
    #[must_use]
    pub fn structural_velocity(&self) -> &[f64] {
        &self.velocity[..self.model.structural_dofs]
    }
    #[must_use]
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        &self.displacement[self.model.structural_dofs..]
    }
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.velocity[self.model.structural_dofs..]
    }
}

/// Original physical storage; small signed summation roundoff is retained.
#[derive(Clone, Copy, Debug, Default)]
pub struct StructuralInteriorTimeEnergy {
    pub structural_kinetic_j: f64,
    pub structural_strain_j: f64,
    pub fluid_kinetic_j: f64,
    pub fluid_compression_j: f64,
}

impl StructuralInteriorTimeEnergy {
    #[must_use]
    pub fn structural_j(self) -> f64 { self.structural_kinetic_j + self.structural_strain_j }
    #[must_use]
    pub fn acoustic_j(self) -> f64 { self.fluid_kinetic_j + self.fluid_compression_j }
    #[must_use]
    pub fn total_j(self) -> f64 { self.structural_j() + self.acoustic_j() }
}

/// Independent source, storage, loss and reciprocal interface accounts [J].
#[derive(Clone, Copy, Debug, Default)]
pub struct StructuralInteriorTimeWork {
    pub initial: StructuralInteriorTimeEnergy,
    pub final_energy: StructuralInteriorTimeEnergy,
    /// Evaluated bilinearly from original matrices, without subtracting large energies.
    pub structural_energy_change_j: f64,
    pub acoustic_energy_change_j: f64,
    pub structural_input_work_j: f64,
    pub acoustic_input_work_j: f64,
    pub structural_dissipation_j: f64,
    pub wall_dissipation_j: f64,
    /// Signed work, positive from the fluid into the structure.
    pub fluid_to_structure_work_j: f64,
    /// Largest allowance used by the three separate physical balance checks.
    pub balance_tolerance_j: f64,
}

impl StructuralInteriorTimeWork {
    #[must_use]
    pub fn structural_balance_residual_j(self) -> f64 {
        self.structural_energy_change_j + self.structural_dissipation_j
            - self.structural_input_work_j - self.fluid_to_structure_work_j
    }
    #[must_use]
    pub fn acoustic_balance_residual_j(self) -> f64 {
        self.acoustic_energy_change_j + self.wall_dissipation_j
            - self.acoustic_input_work_j + self.fluid_to_structure_work_j
    }
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        self.structural_energy_change_j + self.acoustic_energy_change_j
            + self.structural_dissipation_j + self.wall_dissipation_j
            - self.structural_input_work_j - self.acoustic_input_work_j
    }
}

/// Observations from one accepted common trajectory.
#[derive(Clone, Debug)]
pub struct StructuralInteriorTimeFrame {
    pub time_s: f64,
    pub duration_s: f64,
    pub work: StructuralInteriorTimeWork,
    /// Endpoint generalized force on the structure, `B^T p`.
    pub fluid_force: Vec<f64>,
    /// Endpoint original pressure-basis integrals of outward structural velocity.
    pub displaced_flow_m3_s: Vec<f64>,
    /// Exact trajectory integrals, in original structural coordinate units.
    pub integrated_structural_velocity: Vec<f64>,
    /// Exact pressure integrals, in original pressure-basis order [Pa s].
    pub integrated_pressure_pa_s: Vec<f64>,
}

/// Candidate-returning convenience API; caller adopts the complete state.
pub struct StructuralInteriorTimeStep<'m, 'a, 'f, 'c> {
    pub state: StructuralInteriorTimeState<'m, 'a, 'f, 'c>,
    pub frame: StructuralInteriorTimeFrame,
}

/// Reusable candidate and native contraction buffers. Accepted state stays
/// untouched on refusal; scratch contents are not exposed as accepted output.
pub struct StructuralInteriorTimeWorkspace<'m, 'a, 'f, 'c> {
    candidate: StructuralInteriorTimeState<'m, 'a, 'f, 'c>,
    original: Vec<f64>,
    input: Vec<f64>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    integrals: Vec<f64>,
    frame: StructuralInteriorTimeFrame,
}

/// Complete original operator and one fixed held-input interval.
pub struct PreparedStructuralInteriorTime<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    duration_s: f64,
    displacement_scale: Vec<f64>,
    mass_scale: Vec<f64>,
    normalized_mass: Vec<f64>,
    velocity_map: Vec<f64>,
    flow: PreparedQuadraticFlow,
    limits: StructuralInteriorTimeLimits,
}

impl<'a, 'f, 'c> StructuralInterior<'a, 'f, 'c> {
    /// Prepare one whole-system exact transition using the existing pHS owner.
    /// All original forces and projected inward flows are independent held inputs.
    pub fn prepare_exact_time(
        &self,
        duration_s: f64,
        limits: StructuralInteriorTimeLimits,
        cx: &Cx<'_>,
    ) -> Result<PreparedStructuralInteriorTime<'_, 'a, 'f, 'c>, StructuralInteriorTimeError> {
        poll(cx)?;
        let n = self.structural_dofs + self.pressure_dofs();
        let maximum = limits.maximum_coordinates.min((MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2 - 1) / 3);
        if n > maximum {
            return Err(StructuralInteriorTimeError::Limit { requested: n, maximum });
        }
        if !duration_s.is_finite() || duration_s <= 0.0
            || !limits.energy_balance_tolerance.is_finite()
            || limits.energy_balance_tolerance <= 0.0
            || limits.energy_balance_tolerance >= 1.0
        {
            return Err(invalid("positive finite duration and work tolerance in (0,1) required"));
        }
        let states = 2 * n;
        let mut displacement_scale = Vec::with_capacity(n);
        let mut mass_scale = Vec::with_capacity(n);
        for i in 0..n {
            poll(cx)?;
            let m = self.mass[i * n + i];
            let k = self.stiffness[i * n + i];
            if !m.is_finite() || m <= 0.0 || !k.is_finite() || k < 0.0 {
                return Err(invalid("positive mass and nonnegative stiffness diagonals required"));
            }
            mass_scale.push(det::sqrt(m));
            // A free displacement keeps its original unit coordinate. This is
            // a change of numerical scale, never an added physical stiffness.
            displacement_scale.push(if k > 0.0 { det::sqrt(k) } else { 1.0 });
            for j in 0..n {
                if !self.mass[i * n + j].is_finite()
                    || !self.stiffness[i * n + j].is_finite()
                    || !self.damping[i * n + j].is_finite()
                    || self.mass[i * n + j] != self.mass[j * n + i]
                    || self.stiffness[i * n + j] != self.stiffness[j * n + i]
                {
                    return Err(invalid("finite original pencil with exactly symmetric mass/stiffness required"));
                }
            }
        }
        let mut normalized_mass = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let value = self.mass[i * n + j] / mass_scale[i] / mass_scale[j];
                normalized_mass[i * n + j] = value;
                normalized_mass[j * n + i] = value;
            }
        }
        let factor = cholesky_cancellable(&normalized_mass, n, || cx.checkpoint().is_err())
            .map_err(|error| match error {
                FactorRunError::Factor(error) => StructuralInteriorTimeError::Mass(error),
                FactorRunError::Cancelled => StructuralInteriorTimeError::Cancelled,
            })?;
        let mut inverse = vec![0.0; n * n];
        let mut rhs = vec![0.0; n];
        for column in 0..n {
            poll(cx)?;
            rhs.fill(0.0);
            rhs[column] = 1.0;
            factor.solve(&mut rhs);
            for row in 0..n { inverse[row * n + column] = rhs[row]; }
        }
        verify_inverse(&normalized_mass, &inverse, n, cx)?;
        let mut velocity_map = inverse.clone();
        for i in 0..n {
            for j in 0..n { velocity_map[i * n + j] /= mass_scale[i]; }
        }
        let mut q = vec![0.0; states * states];
        let mut j = vec![0.0; states * states];
        let mut r = vec![0.0; states * states];
        let mut g = vec![0.0; states * n];
        // Physical velocity forms: structural loss, relative-wall loss,
        // and signed fluid-to-structure work. Only these observations are split;
        // the generator retains the complete original velocity operator.
        let mut forms = vec![vec![0.0; n * n]; 3];
        for a in 0..n {
            poll(cx)?;
            let canonical = displacement_scale[a] / mass_scale[a];
            j[a * states + n + a] = canonical;
            j[(n + a) * states + a] = -canonical;
            g[(n + a) * n + a] = 1.0 / mass_scale[a];
            for b in a..n {
                let k = self.stiffness[a * n + b] / displacement_scale[a] / displacement_scale[b];
                q[a * states + b] = k;
                q[b * states + a] = k;
                let symmetric = self.damping[a * n + b].midpoint(self.damping[b * n + a]);
                let skew = self.damping[a * n + b].midpoint(-self.damping[b * n + a]);
                let resistance = symmetric / mass_scale[a] / mass_scale[b];
                let gyroscopic = -skew / mass_scale[a] / mass_scale[b];
                r[(n + a) * states + n + b] = resistance;
                r[(n + b) * states + n + a] = resistance;
                j[(n + a) * states + n + b] = gyroscopic;
                j[(n + b) * states + n + a] = -gyroscopic;
                if b < self.structural_dofs {
                    forms[0][a * n + b] = symmetric;
                    forms[0][b * n + a] = symmetric;
                } else if a >= self.structural_dofs {
                    forms[1][a * n + b] = symmetric;
                    forms[1][b * n + a] = symmetric;
                }
            }
            for b in 0..n {
                q[(n + a) * states + n + b] = inverse[a * n + b];
            }
        }
        for a in 0..self.pressure_dofs() {
            for s in 0..self.structural_dofs {
                let b = 0.5 * self.coupling[a * self.structural_dofs + s];
                forms[2][s * n + self.structural_dofs + a] = b;
                forms[2][(self.structural_dofs + a) * n + s] = b;
            }
        }
        let quadratic: Vec<_> = forms.iter().map(|form| {
            velocity_form(form, &velocity_map, n, cx)
        }).collect::<Result<_, _>>()?;
        let linear: Vec<_> = (0..n).map(|a| {
            let mut observation = vec![0.0; states];
            observation[n..].copy_from_slice(&velocity_map[a * n..(a + 1) * n]);
            observation
        }).collect();
        if q.iter().chain(&j).chain(&r).chain(&g).chain(&velocity_map).any(|v| !v.is_finite()) {
            return Err(invalid("canonical coordinate transformation overflow"));
        }
        let storage = QuadraticStorage::new(q, states)?;
        let system = PortHamiltonian::new(states, n, j, r, g, Box::new(storage))?;
        poll(cx)?;
        let split = split_quadratic(&system, LinearSplitBudget {
            max_modes: states,
            max_omega_dt: 1.0,
        })?;
        let flow = PreparedQuadraticFlow::new_controlled(
            &system, &split, duration_s, &quadratic, &linear, limits.preparation,
            || cx.checkpoint().is_err(),
        ).map_err(propagation)?;
        poll(cx)?;
        Ok(PreparedStructuralInteriorTime {
            model: self, duration_s, displacement_scale, mass_scale,
            normalized_mass, velocity_map, flow, limits,
        })
    }
}

impl<'m, 'a, 'f, 'c> PreparedStructuralInteriorTime<'m, 'a, 'f, 'c> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 { self.duration_s }

    /// Initialize all original coordinates at one finite physical time.
    #[allow(clippy::too_many_arguments)] // original structural and pressure state, clock and cancellation
    pub fn state(
        &self,
        time_s: f64,
        structural_displacement: &[f64],
        structural_velocity: &[f64],
        pressure_potential_pa_s: &[f64],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorTimeState<'m, 'a, 'f, 'c>, StructuralInteriorTimeError> {
        poll(cx)?;
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        if !time_s.is_finite()
            || structural_displacement.len() != ns || structural_velocity.len() != ns
            || pressure_potential_pa_s.len() != na || pressure_pa.len() != na
            || structural_displacement.iter().chain(structural_velocity)
                .chain(pressure_potential_pa_s).chain(pressure_pa).any(|v| !v.is_finite())
        {
            return Err(invalid("complete finite original structural/pressure state and clock required"));
        }
        let state = StructuralInteriorTimeState {
            model: self.model, time_s,
            displacement: structural_displacement.iter().chain(pressure_potential_pa_s).copied().collect(),
            velocity: structural_velocity.iter().chain(pressure_pa).copied().collect(),
        };
        self.energy(&state, cx)?;
        Ok(state)
    }

    /// Allocate all scratch once, outside repeated physical stepping.
    #[must_use]
    pub fn workspace(&self) -> StructuralInteriorTimeWorkspace<'m, 'a, 'f, 'c> {
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        let n = ns + na;
        StructuralInteriorTimeWorkspace {
            candidate: StructuralInteriorTimeState {
                model: self.model, time_s: 0.0,
                displacement: vec![0.0; n], velocity: vec![0.0; n],
            },
            original: vec![0.0; 2 * n],
            input: vec![0.0; n],
            augmented: vec![0.0; 3 * n + 1],
            balanced: vec![0.0; 2 * n],
            endpoint: vec![0.0; 2 * n],
            integrals: vec![0.0; n + 3],
            frame: StructuralInteriorTimeFrame {
                time_s: 0.0, duration_s: self.duration_s,
                work: StructuralInteriorTimeWork::default(),
                fluid_force: vec![0.0; ns],
                displaced_flow_m3_s: vec![0.0; na],
                integrated_structural_velocity: vec![0.0; ns],
                integrated_pressure_pa_s: vec![0.0; na],
            },
        }
    }

    /// Original matrix energy, with no coordinate normalization in the result.
    pub fn energy(
        &self,
        state: &StructuralInteriorTimeState<'_, '_, '_, '_>,
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorTimeEnergy, StructuralInteriorTimeError> {
        self.validate(state)?;
        self.energy_with_roundoff(state, cx).map(|value| value.0)
    }

    /// Return a complete candidate without modifying the accepted input.
    pub fn step(
        &self,
        initial: &StructuralInteriorTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        inward_flow_load: &[f64],
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorTimeStep<'m, 'a, 'f, 'c>, StructuralInteriorTimeError> {
        self.validate(initial)?;
        let mut next = initial.clone();
        let mut workspace = self.workspace();
        self.step_into(&mut next, structural_force, inward_flow_load, &mut workspace, cx)?;
        Ok(StructuralInteriorTimeStep { state: next, frame: workspace.frame })
    }

    /// Advance every physical state and the clock in one transaction.
    ///
    /// Inputs are held original generalized forces and projected inward flow
    /// `s_a = integral(phi_a v_in dA)` [m^3/s], in the same order as the harmonic
    /// owner. No area, density or frequency normalization is introduced.
    /// Success publishes the endpoint only after source, structural, acoustic
    /// and total work checks and a final cancellation poll. A failure can alter
    /// private scratch, but leaves the accepted physical state and clock intact.
    pub fn step_into<'w>(
        &self,
        state: &mut StructuralInteriorTimeState<'m, 'a, 'f, 'c>,
        structural_force: &[f64],
        inward_flow_load: &[f64],
        workspace: &'w mut StructuralInteriorTimeWorkspace<'m, 'a, 'f, 'c>,
        cx: &Cx<'_>,
    ) -> Result<&'w StructuralInteriorTimeFrame, StructuralInteriorTimeError> {
        poll(cx)?;
        self.validate(state)?;
        self.validate(&workspace.candidate)?;
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        let n = ns + na;
        if structural_force.len() != ns || inward_flow_load.len() != na
            || structural_force.iter().chain(inward_flow_load).any(|v| !v.is_finite())
        {
            return Err(invalid("complete finite held structural forces and inward flow loads required"));
        }
        let next_time = state.time_s + self.duration_s;
        if !next_time.is_finite() || next_time <= state.time_s {
            return Err(invalid("held interval must advance the finite physical clock"));
        }
        workspace.input[..ns].copy_from_slice(structural_force);
        workspace.input[ns..].copy_from_slice(inward_flow_load);
        let (initial_energy, initial_roundoff) = self.energy_with_roundoff(state, cx)?;
        for i in 0..n {
            poll(cx)?;
            workspace.original[i] = self.displacement_scale[i] * state.displacement[i];
            let mut momentum = Dd::from_f64(0.0);
            for j in 0..n {
                momentum = momentum + Dd::from_f64(self.normalized_mass[i * n + j])
                    * Dd::from_f64(self.mass_scale[j]) * Dd::from_f64(state.velocity[j]);
            }
            workspace.original[n + i] = momentum.to_f64();
        }
        self.flow.evaluate_into(
            &workspace.original, &workspace.input, &mut workspace.augmented,
            &mut workspace.balanced, &mut workspace.endpoint, &mut workspace.integrals,
        )?;
        poll(cx)?;
        for i in 0..n {
            workspace.candidate.displacement[i] = workspace.endpoint[i] / self.displacement_scale[i];
            let mut velocity = Dd::from_f64(0.0);
            for j in 0..n {
                velocity = velocity + Dd::from_f64(self.velocity_map[i * n + j])
                    * Dd::from_f64(workspace.endpoint[n + j]);
            }
            workspace.candidate.velocity[i] = velocity.to_f64();
        }
        workspace.candidate.time_s = next_time;
        if workspace.candidate.displacement.iter().chain(&workspace.candidate.velocity)
            .any(|v| !v.is_finite())
        {
            return Err(invalid("physical exact endpoint is unrepresentable"));
        }
        let (final_energy, final_roundoff) = self.energy_with_roundoff(&workspace.candidate, cx)?;
        let (structural_change, structural_roundoff) = self.energy_change(state, &workspace.candidate, 0, ns, cx)?;
        let (acoustic_change, acoustic_roundoff) = self.energy_change(state, &workspace.candidate, ns, n, cx)?;
        let (structural_input, structural_absolute) = input_work(structural_force, &workspace.integrals[3..3 + ns])?;
        let (acoustic_input, acoustic_absolute) = input_work(inward_flow_load, &workspace.integrals[3 + ns..])?;
        let mut work = StructuralInteriorTimeWork {
            initial: initial_energy, final_energy,
            structural_energy_change_j: structural_change,
            acoustic_energy_change_j: acoustic_change,
            structural_input_work_j: structural_input,
            acoustic_input_work_j: acoustic_input,
            structural_dissipation_j: workspace.integrals[0],
            wall_dissipation_j: workspace.integrals[1],
            fluid_to_structure_work_j: workspace.integrals[2],
            balance_tolerance_j: 0.0,
        };
        let roundoff = [initial_roundoff[0] + final_roundoff[0], initial_roundoff[1] + final_roundoff[1]];
        let structural_limit = self.work_limit(
            work.initial.structural_j(), work.final_energy.structural_j(), structural_change,
            structural_absolute, work.structural_dissipation_j, work.fluid_to_structure_work_j,
            roundoff[0] + structural_roundoff,
        );
        let acoustic_limit = self.work_limit(
            work.initial.acoustic_j(), work.final_energy.acoustic_j(), acoustic_change,
            acoustic_absolute, work.wall_dissipation_j, work.fluid_to_structure_work_j,
            roundoff[1] + acoustic_roundoff,
        );
        let total_limit = self.work_limit(
            work.initial.total_j(), work.final_energy.total_j(), structural_change + acoustic_change,
            structural_absolute + acoustic_absolute,
            work.structural_dissipation_j + work.wall_dissipation_j, 0.0,
            roundoff[0] + roundoff[1] + structural_roundoff + acoustic_roundoff,
        );
        check_balance("structural", work.structural_balance_residual_j(), structural_limit)?;
        check_balance("acoustic", work.acoustic_balance_residual_j(), acoustic_limit)?;
        check_balance("combined", work.balance_residual_j(), total_limit)?;
        for (account, loss, limit) in [
            ("structural dissipation", work.structural_dissipation_j, structural_limit),
            ("wall dissipation", work.wall_dissipation_j, acoustic_limit),
        ] {
            if loss < -limit {
                return Err(StructuralInteriorTimeError::EnergyBalance { account, residual_j: loss, limit_j: limit });
            }
        }
        work.balance_tolerance_j = structural_limit.max(acoustic_limit).max(total_limit);
        workspace.frame.fluid_force.fill(0.0);
        workspace.frame.displaced_flow_m3_s.fill(0.0);
        for a in 0..na {
            poll(cx)?;
            for s in 0..ns {
                let b = self.model.coupling[a * ns + s];
                workspace.frame.fluid_force[s] += b * workspace.candidate.velocity[ns + a];
                workspace.frame.displaced_flow_m3_s[a] += b * workspace.candidate.velocity[s];
            }
        }
        if workspace.frame.fluid_force.iter().chain(&workspace.frame.displaced_flow_m3_s).any(|v| !v.is_finite()) {
            return Err(invalid("endpoint interface observation overflow"));
        }
        workspace.frame.integrated_structural_velocity.copy_from_slice(&workspace.integrals[3..3 + ns]);
        workspace.frame.integrated_pressure_pa_s.copy_from_slice(&workspace.integrals[3 + ns..]);
        workspace.frame.time_s = next_time;
        workspace.frame.duration_s = self.duration_s;
        workspace.frame.work = work;
        poll(cx)?;
        core::mem::swap(state, &mut workspace.candidate);
        Ok(&workspace.frame)
    }

    fn validate(&self, state: &StructuralInteriorTimeState<'_, '_, '_, '_>)
        -> Result<(), StructuralInteriorTimeError>
    {
        if !core::ptr::eq(self.model, state.model) {
            return Err(invalid("time state/workspace belongs to another original structural-interior operator"));
        }
        Ok(())
    }

    fn energy_with_roundoff(
        &self, state: &StructuralInteriorTimeState<'_, '_, '_, '_>, cx: &Cx<'_>,
    ) -> Result<(StructuralInteriorTimeEnergy, [f64; 2]), StructuralInteriorTimeError> {
        let ns = self.model.structural_dofs;
        let n = ns + self.model.pressure_dofs();
        let (sk, rsk) = energy_form(&self.model.mass, n, 0, ns, &state.velocity, cx)?;
        let (ss, rss) = energy_form(&self.model.stiffness, n, 0, ns, &state.displacement, cx)?;
        let (fk, rfk) = energy_form(&self.model.stiffness, n, ns, n, &state.displacement, cx)?;
        let (fc, rfc) = energy_form(&self.model.mass, n, ns, n, &state.velocity, cx)?;
        Ok((StructuralInteriorTimeEnergy {
            structural_kinetic_j: sk, structural_strain_j: ss,
            fluid_kinetic_j: fk, fluid_compression_j: fc,
        }, [rsk + rss, rfk + rfc]))
    }

    fn energy_change(
        &self, initial: &StructuralInteriorTimeState<'_, '_, '_, '_>,
        next: &StructuralInteriorTimeState<'_, '_, '_, '_>,
        start: usize, end: usize, cx: &Cx<'_>,
    ) -> Result<(f64, f64), StructuralInteriorTimeError> {
        let n = self.model.structural_dofs + self.model.pressure_dofs();
        let mut change = Dd::from_f64(0.0);
        let mut absolute = 0.0;
        for (matrix, before, after) in [
            (&self.model.mass, &initial.velocity, &next.velocity),
            (&self.model.stiffness, &initial.displacement, &next.displacement),
        ] {
            for i in start..end {
                poll(cx)?;
                for j in start..end {
                    let term = Dd::from_f64(0.5 * matrix[i * n + j])
                        * (Dd::from_f64(after[i]) - Dd::from_f64(before[i]))
                        * (Dd::from_f64(after[j]) + Dd::from_f64(before[j]));
                    change = change + term;
                    absolute += term.to_f64().abs();
                }
            }
        }
        finite_account(change.to_f64(), absolute)?;
        Ok((change.to_f64(), 128.0 * f64::EPSILON * n as f64 * absolute))
    }

    #[allow(clippy::too_many_arguments)] // independent original physical accounts and evaluated roundoff
    fn work_limit(&self, initial: f64, final_energy: f64, change: f64,
        absolute_input: f64, loss: f64, exchange: f64, roundoff: f64) -> f64
    {
        let physical_scale = initial.abs().max(final_energy.abs()).max(change.abs())
            .max(absolute_input).max(loss.abs()).max(exchange.abs());
        self.limits.energy_balance_tolerance * physical_scale + roundoff
            + 128.0 * f64::EPSILON * (self.model.structural_dofs + self.model.pressure_dofs()) as f64
                * (absolute_input + loss.abs() + exchange.abs())
    }
}

fn invalid(what: &'static str) -> StructuralInteriorTimeError {
    StructuralInteriorTimeError::Invalid(what)
}
fn poll(cx: &Cx<'_>) -> Result<(), StructuralInteriorTimeError> {
    cx.checkpoint().map_err(|_| StructuralInteriorTimeError::Cancelled)
}
fn propagation(error: PreparedStepError) -> StructuralInteriorTimeError {
    match error {
        PreparedStepError::Cancelled => StructuralInteriorTimeError::Cancelled,
        PreparedStepError::Solver(error) => StructuralInteriorTimeError::Dynamics(error),
    }
}

fn verify_inverse(matrix: &[f64], inverse: &[f64], n: usize, cx: &Cx<'_>)
    -> Result<(), StructuralInteriorTimeError>
{
    for i in 0..n {
        poll(cx)?;
        for j in 0..n {
            let mut value = Dd::from_f64(0.0);
            let mut bound = 0.0;
            for k in 0..n {
                let term = Dd::from_f64(matrix[i * n + k]) * Dd::from_f64(inverse[k * n + j]);
                value = value + term;
                bound += term.to_f64().abs();
            }
            let expected: f64 = if i == j { 1.0 } else { 0.0 };
            if !value.to_f64().is_finite() || !bound.is_finite()
                || (value.to_f64() - expected).abs() > 512.0 * f64::EPSILON * n as f64 * bound.max(expected)
            {
                return Err(invalid("original complete mass inverse failed its residual"));
            }
        }
    }
    Ok(())
}

// Pull back the physical v^T L v observation through v = velocity_map * momentum.
// One product per symmetric pair retains exact observation symmetry.
fn velocity_form(form: &[f64], velocity: &[f64], n: usize, cx: &Cx<'_>)
    -> Result<Vec<f64>, StructuralInteriorTimeError>
{
    let mut right = vec![0.0; n * n];
    for i in 0..n {
        poll(cx)?;
        for j in 0..n {
            let mut value = Dd::from_f64(0.0);
            for k in 0..n {
                value = value + Dd::from_f64(form[i * n + k]) * Dd::from_f64(velocity[k * n + j]);
            }
            right[i * n + j] = value.to_f64();
        }
    }
    let states = 2 * n;
    let mut result = vec![0.0; states * states];
    for i in 0..n {
        poll(cx)?;
        for j in i..n {
            let mut value = Dd::from_f64(0.0);
            for k in 0..n {
                value = value + Dd::from_f64(velocity[k * n + i]) * Dd::from_f64(right[k * n + j]);
            }
            result[(n + i) * states + n + j] = value.to_f64();
            result[(n + j) * states + n + i] = value.to_f64();
        }
    }
    if result.iter().any(|v| !v.is_finite()) { return Err(invalid("physical work observation overflow")); }
    Ok(result)
}

fn energy_form(matrix: &[f64], n: usize, start: usize, end: usize,
    state: &[f64], cx: &Cx<'_>) -> Result<(f64, f64), StructuralInteriorTimeError>
{
    let mut energy = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for i in start..end {
        poll(cx)?;
        for j in start..end {
            let term = Dd::from_f64(0.5 * matrix[i * n + j])
                * Dd::from_f64(state[i]) * Dd::from_f64(state[j]);
            energy = energy + term;
            absolute += term.to_f64().abs();
        }
    }
    let value = energy.to_f64();
    finite_account(value, absolute)?;
    let roundoff = 128.0 * f64::EPSILON * n as f64 * absolute;
    if value < -roundoff { return Err(invalid("original physical storage is materially negative")); }
    Ok((value, roundoff))
}
fn input_work(input: &[f64], integral: &[f64]) -> Result<(f64, f64), StructuralInteriorTimeError> {
    let mut value = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for (&input, &integral) in input.iter().zip(integral) {
        let term = Dd::from_f64(input) * Dd::from_f64(integral);
        value = value + term;
        absolute += term.to_f64().abs();
    }
    finite_account(value.to_f64(), absolute)?;
    Ok((value.to_f64(), absolute))
}
fn finite_account(value: f64, absolute: f64) -> Result<(), StructuralInteriorTimeError> {
    if !value.is_finite() || !absolute.is_finite() { return Err(invalid("physical energy/work is unrepresentable")); }
    Ok(())
}
fn check_balance(account: &'static str, residual: f64, limit: f64)
    -> Result<(), StructuralInteriorTimeError>
{
    if !residual.is_finite() || !limit.is_finite() {
        return Err(invalid("physical work balance is unrepresentable"));
    }
    if residual.abs() > limit {
        return Err(StructuralInteriorTimeError::EnergyBalance { account, residual_j: residual, limit_j: limit });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
