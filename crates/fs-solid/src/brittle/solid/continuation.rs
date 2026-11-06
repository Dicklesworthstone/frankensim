//! Proportional-force equilibrium paths through native tensile softening.
//!
//! The existing Keller bordered solver advances dimensionless free coordinates
//! `q = u/L` and load factor `lambda`; its equation is `(f_int-lambda*f_ref)/F`.
//! Every force and consistent tangent evaluates the ORIGINAL crack-band trial
//! from the same last ACCEPTED histories. Trial damage is allowed to evolve;
//! no Newton iterate, failed corrector or numerical predictor commits damage.
//! Native P1 reference geometry is prepared once per endpoint and reused by
//! every residual/Jacobian assembly. Fixed supports are eliminated exactly.
//!
//! At an accepted loading envelope the material deliberately selects its
//! unloading tangent. Reusing that tangent as the next physical predictor can
//! turn a descending fracture branch into unloading. Subsequent steps therefore
//! retain the actual accepted secant direction; the corrector still uses the
//! native current trial tangent. An explicit direction can start a damaged
//! branch or request a reversal. Numerical arc length is not physical time.
//!
//! Publication shares the original solid's work, irreversibility and released-
//! face gates, plus independent free-force, damage-increment and work budgets.
//! A sharp constitutive corner may require a smaller step or explicit direction.
//! This is a path-following equilibrium model, not a stable dynamic trajectory,
//! unique branch selection, crack-path certificate, moving-support method or
//! automatic tensile-event localization. Source conditions and topology stay
//! fixed. Endpoint work error remains distinct from physical fracture loss.

use std::cell::{Cell, RefCell};
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::{BrittleEnergy, BrittleSmallStrain};
use fs_sparse::Csr;

use super::{BrittleTetControls, BrittleTetEnergy, BrittleTetError, BrittleTetSolid, checkpoint};
use crate::SolidError;
use crate::continuation::{self, ArcSettings, PathResidual, PathState};
use crate::plastic3::{
    IncrementRecord, IncrementSettings, PreparedSmallStrainTet, SmallStrainLawField,
    SmallStrainTetError, SmallStrainTetProblem, SmallStrainTetSolution,
};

/// Explicit SI scales for the dimensionless force equation and arc metric.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetPathScales {
    pub displacement_m: f64,
    pub force_n: f64,
}

/// Retained trace bound and native assembly allowance per proposed endpoint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BrittleTetArcBudget {
    pub maximum_steps: usize,
    /// Counts complete native force/tangent assemblies, including independent
    /// initial and final force checks. Local energy recovery is separately
    /// bounded by the original element count and runs only after convergence.
    pub maximum_operator_evaluations: usize,
}

/// Physical gates applied in addition to the dimensionless solver tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetArcControls {
    pub work: BrittleTetControls,
    pub maximum_damage_increment: f64,
}

/// Optional physical direction for a first damaged step or deliberate reversal.
/// Supports must have exactly zero direction. Magnitude is immaterial after
/// normalizing by the retained displacement scale and the full arc norm.
#[derive(Clone, Copy, Debug)]
pub struct BrittleTetArcDirection<'a> {
    pub displacement_m: &'a [[f64; 3]],
    pub load_factor: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetContinuationError {
    Invalid(&'static str),
    Solid(BrittleTetError),
    Path(SolidError),
    InitialEquilibrium {
        scaled_residual_norm: f64,
        tolerance: f64,
    },
    ForceResidual {
        scaled_residual_norm: f64,
        tolerance: f64,
    },
    DamageIncrement {
        observed: f64,
        maximum: f64,
    },
    StepBudget,
    OperatorBudget {
        evaluations: usize,
    },
}

impl core::fmt::Display for BrittleTetContinuationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "brittle arclength endpoint refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetContinuationError {}
impl From<BrittleTetError> for BrittleTetContinuationError {
    fn from(error: BrittleTetError) -> Self {
        Self::Solid(error)
    }
}
impl From<SmallStrainTetError> for BrittleTetContinuationError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error.into())
    }
}

#[derive(Debug)]
struct Binding {
    prescribed_m: Vec<(usize, f64)>,
    free_dofs: Vec<usize>,
    full_to_free: Vec<usize>,
    reference_force_n: Vec<f64>,
    reference_force_scaled: Vec<f64>,
    scales: BrittleTetPathScales,
    tangent_scale: f64,
    budget: BrittleTetArcBudget,
}

/// One immutable accepted mechanical state and its original path definition.
/// Source owners, physical history and continuation coordinates cannot diverge
/// through public mutation. Clone/retry is performed by each returned step.
#[derive(Debug)]
pub struct BrittleTetContinuation<'a> {
    solid: BrittleTetSolid<'a>,
    binding: Arc<Binding>,
    state: PathState,
}

/// A complete accepted equilibrium; no partially advanced solid is exposed.
#[derive(Debug)]
pub struct BrittleTetArcStep<'a> {
    pub continuation: BrittleTetContinuation<'a>,
    pub energy: BrittleTetEnergy,
    pub element_energy: Vec<BrittleEnergy>,
    pub maximum_damage_increment: f64,
    pub operator_evaluations: usize,
    /// Independent physical Euclidean norm over FREE force residuals, N.
    pub force_residual_norm_n: f64,
    /// Sign change of consecutive accepted load increments; a bracketed
    /// observation, not an exact extremum or stability certificate.
    pub load_direction_changed: bool,
}

impl<'a> BrittleTetContinuation<'a> {
    /// Begin at an existing equilibrium under `lambda * reference_force_n`.
    ///
    /// Fixed support coordinates must equal the original accepted displacement
    /// exactly, including a nonzero retained preload. Repeated force entries
    /// accumulate in input order; support forces remain in reactions. The first
    /// automatic predictor follows increasing load. Starting continued damage
    /// or unloading from a pre-damaged state may need an explicit direction.
    #[allow(clippy::too_many_arguments)]
    pub fn try_from_solid(
        solid: &BrittleTetSolid<'a>,
        prescribed_m: &[(usize, f64)],
        reference_force_n: &[(usize, f64)],
        load_factor: f64,
        scales: BrittleTetPathScales,
        budget: BrittleTetArcBudget,
        settings: &ArcSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetContinuationError> {
        checkpoint(cx)?;
        let old = solid.solution.displacement_m.as_flattened();
        let dofs = old.len();
        let tangent_scale = scales.displacement_m / scales.force_n;
        if ![scales.displacement_m, scales.force_n, tangent_scale]
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
            || !load_factor.is_finite()
            || prescribed_m.len() > dofs
            || reference_force_n.len() > dofs
            || budget.maximum_steps == 0
            || budget.maximum_steps.checked_add(1).is_none()
            || budget.maximum_operator_evaluations == 0
        {
            return Err(BrittleTetContinuationError::Invalid(
                "positive finite scales, bounded original loads and positive work limits required",
            ));
        }
        let mut full_to_free = vec![0; dofs];
        for &(dof, value) in prescribed_m {
            checkpoint(cx)?;
            if dof >= dofs
                || !value.is_finite()
                || value != old[dof]
                || full_to_free[dof] == usize::MAX
            {
                return Err(BrittleTetContinuationError::Invalid(
                    "distinct fixed supports must equal retained physical coordinates",
                ));
            }
            full_to_free[dof] = usize::MAX;
        }
        let mut free_dofs = Vec::new();
        for (dof, free) in full_to_free.iter_mut().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if *free != usize::MAX {
                *free = free_dofs.len();
                free_dofs.push(dof);
            }
        }
        if free_dofs.is_empty() || free_dofs.len() > solid.budget.maximum_free_dofs {
            return Err(BrittleTetContinuationError::Invalid(
                "bounded nonempty free coordinates required",
            ));
        }
        let mut forces = vec![0.0; dofs];
        for &(dof, force) in reference_force_n {
            checkpoint(cx)?;
            if dof >= dofs || !force.is_finite() {
                return Err(BrittleTetContinuationError::Invalid(
                    "invalid reference force",
                ));
            }
            forces[dof] += force;
            if !forces[dof].is_finite() || !(load_factor * forces[dof]).is_finite() {
                return Err(BrittleTetContinuationError::Invalid(
                    "unrepresentable accumulated force",
                ));
            }
        }
        let reference_force_scaled: Vec<_> = free_dofs
            .iter()
            .map(|&dof| forces[dof] / scales.force_n)
            .collect();
        if reference_force_scaled
            .iter()
            .any(|value| !value.is_finite())
            || !reference_force_scaled.iter().any(|&value| value != 0.0)
        {
            return Err(BrittleTetContinuationError::Invalid(
                "finite nonzero scaled free force required",
            ));
        }
        let mut state = PathState::start(free_dofs.len(), settings);
        state.u = free_dofs
            .iter()
            .map(|&dof| old[dof] / scales.displacement_m)
            .collect();
        state.lambda = load_factor;
        state.trace = vec![(load_factor, max_abs(&state.u))];
        let owner = Self {
            solid: retain(solid),
            binding: Arc::new(Binding {
                prescribed_m: prescribed_m.to_vec(),
                free_dofs,
                full_to_free,
                reference_force_n: forces,
                reference_force_scaled,
                scales,
                tangent_scale,
                budget,
            }),
            state,
        };
        let laws: Vec<_> = owner
            .solid
            .materials
            .iter()
            .map(|m| m.small_strain_law())
            .collect();
        let native = owner.problem(&laws);
        let prepared = native.prepare_evaluator(cx)?;
        let adapter = Adapter::new(&owner, &prepared, cx);
        let mut check = owner.state.clone();
        continuation::advance(&adapter, &mut check, settings, 0)
            .map_err(|error| adapter.recover(error))?;
        let residual = adapter.residual_at(&owner.state.u, owner.state.lambda)?;
        let residual_norm = norm(&residual);
        if residual_norm > settings.tol {
            return Err(BrittleTetContinuationError::InitialEquilibrium {
                scaled_residual_norm: residual_norm,
                tolerance: settings.tol,
            });
        }
        owner
            .solid
            .require_released_faces(&owner.solid.solution.displacement_m, cx)?;
        checkpoint(cx)?;
        Ok(owner)
    }

    #[must_use]
    pub const fn solid(&self) -> &BrittleTetSolid<'a> {
        &self.solid
    }
    #[must_use]
    pub const fn load_factor(&self) -> f64 {
        self.state.lambda
    }
    #[must_use]
    pub const fn accepted_steps(&self) -> usize {
        self.state.step
    }
    #[must_use]
    pub const fn next_step_size(&self) -> f64 {
        self.state.ds
    }
    #[must_use]
    pub fn free_dofs(&self) -> &[usize] {
        &self.binding.free_dofs
    }
    #[must_use]
    pub fn reference_force_n(&self) -> &[f64] {
        &self.binding.reference_force_n
    }
    /// `(load factor, maximum |u/L|)` observations, including the initial state.
    #[must_use]
    pub fn trace(&self) -> &[(f64, f64)] {
        &self.state.trace
    }

    fn problem<'p>(
        &'p self,
        laws: &'p [BrittleSmallStrain<'a>],
    ) -> SmallStrainTetProblem<'p, BrittleSmallStrain<'a>> {
        SmallStrainTetProblem {
            nodes_m: self.solid.geometry.nodes_m(),
            tetrahedra: self.solid.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(laws),
            prescribed_m: &self.binding.prescribed_m,
            nodal_forces_n: &[],
            load_path: &[1.0],
            budget: self.solid.budget,
            // Evaluation uses no Newton solve or numerical constitutive fallback.
            settings: IncrementSettings::default(),
        }
    }

    /// Propose one complete arc endpoint. `None` preserves the previous physical
    /// secant (or increasing load on the first step). A work/damage refusal
    /// leaves this owner unchanged; retry with smaller `settings.ds_max`.
    pub fn advance(
        &self,
        settings: &ArcSettings,
        controls: BrittleTetArcControls,
        direction: Option<BrittleTetArcDirection<'_>>,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetArcStep<'a>, BrittleTetContinuationError> {
        checkpoint(cx)?;
        if self.state.step >= self.binding.budget.maximum_steps {
            return Err(BrittleTetContinuationError::StepBudget);
        }
        if !controls.maximum_damage_increment.is_finite()
            || !(0.0..=1.0).contains(&controls.maximum_damage_increment)
            || [
                controls.work.material_energy_tolerance_j,
                controls.work.equilibrium_work_tolerance_j,
                controls.work.maximum_endpoint_work_defect_j,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(BrittleTetContinuationError::Invalid(
                "finite nonnegative work and damage gates required",
            ));
        }
        let laws: Vec<_> = self
            .solid
            .materials
            .iter()
            .map(|m| m.small_strain_law())
            .collect();
        let native = self.problem(&laws);
        let prepared = native.prepare_evaluator(cx)?;
        let adapter = Adapter::new(self, &prepared, cx);
        let mut next = self.state.clone();
        // Lowering the explicit maximum permits a retry after a physical gate.
        next.ds = next.ds.min(settings.ds_max);
        if let Some(direction) = direction {
            self.set_direction(&mut next, direction)?;
        }
        continuation::advance(&adapter, &mut next, settings, 0)
            .map_err(|error| adapter.recover(error))?;
        let initial = norm(&adapter.residual_at(&self.state.u, self.state.lambda)?);
        if initial > settings.tol {
            return Err(BrittleTetContinuationError::InitialEquilibrium {
                scaled_residual_norm: initial,
                tolerance: settings.tol,
            });
        }
        continuation::advance(&adapter, &mut next, settings, 1)
            .map_err(|error| adapter.recover(error))?;
        let displacement = self.displacement(&next.u, next.lambda)?;
        let internal = adapter.internal(&next.u, next.lambda)?;
        let mut residual = Vec::with_capacity(internal.len());
        for (dof, (&force, &reference)) in internal
            .iter()
            .zip(&self.binding.reference_force_n)
            .enumerate()
        {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            let value = force - next.lambda * reference;
            if !value.is_finite() {
                return Err(BrittleTetContinuationError::Invalid(
                    "nonfinite final physical residual",
                ));
            }
            residual.push(value);
        }
        let physical_norm = self
            .binding
            .free_dofs
            .iter()
            .fold(0.0_f64, |n, &dof| n.hypot(residual[dof]));
        let scaled_norm = physical_norm / self.binding.scales.force_n;
        if !scaled_norm.is_finite() || scaled_norm > settings.tol {
            return Err(BrittleTetContinuationError::ForceResidual {
                scaled_residual_norm: scaled_norm,
                tolerance: settings.tol,
            });
        }
        let strains = prepared.strain_field(&displacement, cx)?;
        let mut states = Vec::with_capacity(strains.len());
        let mut stresses = Vec::with_capacity(strains.len());
        let mut maximum_damage_increment = 0.0_f64;
        let mut evolving_elements = 0;
        for (element, (source, strain)) in self.solid.materials.iter().zip(&strains).enumerate() {
            checkpoint(cx)?;
            let old = &self.solid.solution.element_states[element];
            let trial = source
                .law()
                .trial(strain, old)
                .map_err(|source| BrittleTetError::Material { element, source })?;
            maximum_damage_increment =
                maximum_damage_increment.max(trial.state.damage() - old.damage());
            evolving_elements += usize::from(trial.state != *old);
            states.push(trial.state);
            stresses.push(trial.stress_pa);
        }
        if maximum_damage_increment > controls.maximum_damage_increment {
            return Err(BrittleTetContinuationError::DamageIncrement {
                observed: maximum_damage_increment,
                maximum: controls.maximum_damage_increment,
            });
        }
        let solution = SmallStrainTetSolution {
            displacement_m: displacement.as_chunks::<3>().0.to_vec(),
            element_strain: strains,
            element_stress: stresses,
            element_states: states,
            increments: vec![IncrementRecord {
                load_factor: next.lambda,
                residual_history: vec![initial, scaled_norm],
                reactions_n: self
                    .binding
                    .prescribed_m
                    .iter()
                    .map(|&(dof, _)| residual[dof])
                    .collect(),
                evolving_elements,
                indefinite_fallbacks: 0,
            }],
        };
        let forces: Vec<_> = self
            .binding
            .reference_force_n
            .iter()
            .enumerate()
            .filter(|(_, force)| **force != 0.0)
            .map(|(dof, &force)| (dof, next.lambda * force))
            .collect();
        let accepted = self.solid.accept_equilibrium(
            solution,
            &self.binding.prescribed_m,
            &forces,
            controls.work,
            cx,
        )?;
        let mut secant: Vec<_> = next
            .u
            .iter()
            .zip(&self.state.u)
            .map(|(b, a)| b - a)
            .collect();
        let load_increment = next.lambda - self.state.lambda;
        let length = norm(&secant).hypot(load_increment);
        if !length.is_finite() || length == 0.0 {
            return Err(BrittleTetContinuationError::Invalid(
                "arclength endpoint made no representable progress",
            ));
        }
        let load_direction_changed =
            self.state.step > 0 && load_increment * self.state.tangent.1 < 0.0;
        for value in &mut secant {
            *value /= length;
        }
        next.tangent = (secant, load_increment / length);
        next.pending_switch = true;
        checkpoint(cx)?;
        Ok(BrittleTetArcStep {
            continuation: Self {
                solid: accepted.solid,
                binding: Arc::clone(&self.binding),
                state: next,
            },
            energy: accepted.energy,
            element_energy: accepted.element_energy,
            maximum_damage_increment,
            operator_evaluations: adapter.evaluations.get(),
            force_residual_norm_n: physical_norm,
            load_direction_changed,
        })
    }

    fn set_direction(
        &self,
        state: &mut PathState,
        direction: BrittleTetArcDirection<'_>,
    ) -> Result<(), BrittleTetContinuationError> {
        let full = direction.displacement_m.as_flattened();
        if full.len() != self.binding.full_to_free.len()
            || full.iter().any(|value| !value.is_finite())
            || !direction.load_factor.is_finite()
            || self
                .binding
                .prescribed_m
                .iter()
                .any(|&(dof, _)| full[dof] != 0.0)
        {
            return Err(BrittleTetContinuationError::Invalid(
                "finite complete direction with immobile supports required",
            ));
        }
        let mut q: Vec<_> = self
            .binding
            .free_dofs
            .iter()
            .map(|&dof| full[dof] / self.binding.scales.displacement_m)
            .collect();
        let length = norm(&q).hypot(direction.load_factor);
        if !length.is_finite() || length == 0.0 {
            return Err(BrittleTetContinuationError::Invalid(
                "nonzero representable arc direction required",
            ));
        }
        for value in &mut q {
            *value /= length;
        }
        state.tangent = (q, direction.load_factor / length);
        state.pending_switch = true;
        Ok(())
    }

    fn displacement(
        &self,
        q: &[f64],
        lambda: f64,
    ) -> Result<Vec<f64>, BrittleTetContinuationError> {
        if q.len() != self.binding.free_dofs.len()
            || !lambda.is_finite()
            || q.iter().any(|value| !value.is_finite())
            || self
                .binding
                .reference_force_n
                .iter()
                .any(|force| !(lambda * force).is_finite())
        {
            return Err(BrittleTetContinuationError::Invalid(
                "finite dimension-matched arc coordinates and force required",
            ));
        }
        let mut full = vec![0.0; self.binding.full_to_free.len()];
        for &(dof, value) in &self.binding.prescribed_m {
            full[dof] = value;
        }
        for (&dof, &value) in self.binding.free_dofs.iter().zip(q) {
            full[dof] = value * self.binding.scales.displacement_m;
            if !full[dof].is_finite() {
                return Err(BrittleTetContinuationError::Invalid(
                    "unrepresentable physical displacement",
                ));
            }
        }
        Ok(full)
    }
}

struct Adapter<'s, 'p, 'a, 'c, 'cx> {
    owner: &'s BrittleTetContinuation<'a>,
    prepared: &'s PreparedSmallStrainTet<'p, 's, BrittleSmallStrain<'a>>,
    cx: &'c Cx<'cx>,
    evaluations: Cell<usize>,
    native_error: RefCell<Option<BrittleTetContinuationError>>,
}

impl<'s, 'p, 'a, 'c, 'cx> Adapter<'s, 'p, 'a, 'c, 'cx> {
    fn new(
        owner: &'s BrittleTetContinuation<'a>,
        prepared: &'s PreparedSmallStrainTet<'p, 's, BrittleSmallStrain<'a>>,
        cx: &'c Cx<'cx>,
    ) -> Self {
        Self {
            owner,
            prepared,
            cx,
            evaluations: Cell::new(0),
            native_error: RefCell::new(None),
        }
    }
    fn charge(&self) -> Result<(), BrittleTetContinuationError> {
        checkpoint(self.cx)?;
        if self.evaluations.get() >= self.owner.binding.budget.maximum_operator_evaluations {
            return Err(BrittleTetContinuationError::OperatorBudget {
                evaluations: self.evaluations.get(),
            });
        }
        self.evaluations.set(self.evaluations.get() + 1);
        Ok(())
    }
    fn internal(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, BrittleTetContinuationError> {
        self.charge()?;
        let full = self.owner.displacement(q, lambda)?;
        Ok(self.prepared.internal_force(
            &full,
            &self.owner.solid.solution.element_states,
            self.cx,
        )?)
    }
    fn residual_at(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, BrittleTetContinuationError> {
        let internal = self.internal(q, lambda)?;
        let result: Vec<_> = self
            .owner
            .binding
            .free_dofs
            .iter()
            .zip(&self.owner.binding.reference_force_scaled)
            .map(|(&dof, &force)| {
                internal[dof] / self.owner.binding.scales.force_n - lambda * force
            })
            .collect();
        if result.iter().any(|value| !value.is_finite()) {
            return Err(BrittleTetContinuationError::Invalid(
                "nonfinite normalized free-force residual",
            ));
        }
        Ok(result)
    }
    fn tangent_at(&self, q: &[f64], lambda: f64) -> Result<Csr, BrittleTetContinuationError> {
        self.charge()?;
        let full = self.owner.displacement(q, lambda)?;
        let matrix = self.prepared.tangent_matrix(
            &full,
            &self.owner.solid.solution.element_states,
            self.cx,
        )?;
        let binding = &self.owner.binding;
        let n = binding.free_dofs.len();
        let mut rows = Vec::with_capacity(n + 1);
        let mut columns = Vec::new();
        let mut values = Vec::new();
        rows.push(0);
        for &dof in &binding.free_dofs {
            checkpoint(self.cx)?;
            let (source_columns, source_values) = matrix.row(dof);
            for (&full_column, &value) in source_columns.iter().zip(source_values) {
                let column = binding.full_to_free[full_column];
                if column != usize::MAX {
                    let value = value * binding.tangent_scale;
                    if !value.is_finite() {
                        return Err(BrittleTetContinuationError::Invalid(
                            "nonfinite normalized tangent",
                        ));
                    }
                    columns.push(column);
                    values.push(value);
                }
            }
            rows.push(columns.len());
        }
        // Native sorted columns and ordered free DOFs preserve canonical CSR.
        Ok(Csr::from_parts(n, n, rows, columns, values))
    }
    fn bridge<T>(&self, result: Result<T, BrittleTetContinuationError>) -> Result<T, SolidError> {
        *self.native_error.borrow_mut() = None;
        result.map_err(|error| {
            *self.native_error.borrow_mut() = Some(error.clone());
            match error {
                BrittleTetContinuationError::Solid(BrittleTetError::Cancelled)
                | BrittleTetContinuationError::OperatorBudget { .. } => SolidError::Cancelled,
                other => SolidError::InvalidInput {
                    what: other.to_string(),
                },
            }
        })
    }
    fn recover(&self, error: SolidError) -> BrittleTetContinuationError {
        self.native_error
            .borrow_mut()
            .take()
            .unwrap_or(BrittleTetContinuationError::Path(error))
    }
}

impl PathResidual for Adapter<'_, '_, '_, '_, '_> {
    fn checkpoint(&self) -> Result<(), SolidError> {
        self.bridge(checkpoint(self.cx).map_err(Into::into))
    }
    fn ndof(&self) -> usize {
        self.owner.binding.free_dofs.len()
    }
    fn residual(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.bridge(self.residual_at(q, lambda))
    }
    fn tangent(&self, q: &[f64], lambda: f64) -> Result<Csr, SolidError> {
        self.bridge(self.tangent_at(q, lambda))
    }
    fn load_vector(&self) -> Vec<f64> {
        self.owner.binding.reference_force_scaled.clone()
    }
}

fn norm(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |norm, &value| norm.hypot(value))
}
fn max_abs(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |max, &value| max.max(value.abs()))
}
fn retain<'a>(solid: &BrittleTetSolid<'a>) -> BrittleTetSolid<'a> {
    BrittleTetSolid {
        geometry: solid.geometry,
        materials: solid.materials,
        bands: Arc::clone(&solid.bands),
        released_faces: Arc::clone(&solid.released_faces),
        budget: solid.budget,
        solution: solid.solution.clone(),
        stored_energy_j: solid.stored_energy_j,
        fracture_dissipation_j: solid.fracture_dissipation_j,
    }
}
