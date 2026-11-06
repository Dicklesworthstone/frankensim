//! Source-bound inertial evolution of the original retained brittle solid.
//!
//! The native P1 backward-Euler owner advances displacement, velocity and the
//! original fracture histories together. A consistent reference mass is
//! assembled once from the same source bundles and original geometry, then
//! retained across every step. Tensile failure never deletes mass or resets
//! damage. A supported static limit is not imposed on an inertial trajectory.
//!
//! External and support work are balanced against elastic storage, physical
//! fracture loss, kinetic change, constitutive endpoint quadrature error and
//! backward-Euler numerical kinetic loss as SEPARATE accounts. Numerical
//! damping is not fracture dissipation or heat. First-order time and fixed-plane
//! small-strain limitations remain. The [`time`] path refines native timesteps
//! against local state comparisons, damage jumps and numerical-loss budgets.
//! The [`fracture`] path preserves this same motion through explicitly completed
//! traction-free cuts. Closing released faces requires a contact response;
//! arbitrary fragment contact and finite-strain mass mapping are absent.

/// Analytic material path work with centered inertia and no numerical kinetic loss.
pub mod midpoint;

/// Temperature-dependent strength with unchanged reference mass and retained motion.
pub mod temperature;

/// Native step doubling with bounded damage, state error and numerical losses.
pub mod time;

/// Conservative continuation after complete traction-free fracture cuts.
pub mod fracture;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::BrittleEnergy;
use fs_material::brittle::source::BrittleSourceError;

use super::{BrittleTetControls, BrittleTetError, BrittleTetSolid, checkpoint, energy_gate};
use crate::inertia::P1ReferenceMass;
use crate::plastic3::dynamics::InertialEnergy;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
};

/// Independent physical and numerical work gates, all in joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetDynamicControls {
    /// Material/cumulative-fracture and force-equilibrium gates. The endpoint
    /// allowance bounds constitutive quadrature error, excluding kinetic loss.
    /// The material tolerance also gates the combined full energy balance.
    pub work: BrittleTetControls,
    /// Native kinetic-work identity and retained kinetic-storage continuity.
    pub kinetic_work_tolerance_j: f64,
    /// Maximum backward-Euler numerical kinetic loss for this physical step.
    /// Reducing dt is a remedy; this is not a physical damping coefficient.
    pub maximum_numerical_kinetic_loss_j: f64,
}

/// Whole-solid accounts for one accepted physical-time increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetDynamicEnergy {
    /// Original consistent-mass kinetic/work receipt.
    pub inertia: InertialEnergy,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub material_work_j: f64,
    pub fracture_dissipation_j: f64,
    pub material_balance_error_j: f64,
    pub internal_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    /// Full native reactions, including consistent-mass support coupling.
    pub support_endpoint_work_j: f64,
    /// External work minus material internal endpoint work and inertial work.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus exact native constitutive path work.
    pub constitutive_endpoint_defect_j: f64,
    /// External work minus elastic change, fracture, kinetic change and the
    /// two separately reported numerical accounts.
    pub total_balance_error_j: f64,
}

impl BrittleTetDynamicEnergy {
    #[must_use]
    pub fn external_endpoint_work_j(&self) -> f64 {
        self.dead_load_endpoint_work_j + self.support_endpoint_work_j
    }
}

/// Immutable accepted motion; retrying a rejected step preserves this owner.
#[derive(Debug)]
pub struct BrittleTetMotion<'a> {
    solid: BrittleTetSolid<'a>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<P1ReferenceMass>,
    kinetic_energy_j: f64,
    time_s: f64,
    maximum_steps: usize,
    accepted_steps: usize,
}

/// Complete accepted candidate and original per-element constitutive work.
#[derive(Debug)]
pub struct BrittleTetMotionStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    pub element_energy: Vec<BrittleEnergy>,
    pub energy: BrittleTetDynamicEnergy,
}

/// Source, native solid/inertia, time budget or coupled-work refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetDynamicError {
    Invalid(&'static str),
    Source {
        element: usize,
        source: BrittleSourceError,
    },
    Solid(BrittleTetError),
    StepBudget {
        accepted_steps: usize,
        maximum_steps: usize,
    },
}
impl core::fmt::Display for BrittleTetDynamicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "brittle physical-time endpoint refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetDynamicError {}
impl From<BrittleTetError> for BrittleTetDynamicError {
    fn from(error: BrittleTetError) -> Self {
        Self::Solid(error)
    }
}
impl From<SmallStrainTetError> for BrittleTetDynamicError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Start motion from an existing accepted solid, retaining all fracture
    /// history and storage. The caller supplies actual initial nodal velocity.
    ///
    /// Each material must include positive source-resolved reference density.
    /// This explicitly binds that density snapshot to the supplied reference
    /// geometry: no conversion from a current/hot volume is inferred. Its mass
    /// inventory is frozen for the entire motion. Starting motion does not
    /// perform an impulse, claim a heat transfer or erase the supplied kinetic
    /// energy. The original static owner remains usable on every refusal.
    pub fn try_from_solid(
        solid: &BrittleTetSolid<'a>,
        initial_velocity_mps: &[[f64; 3]],
        initial_time_s: f64,
        maximum_steps: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetDynamicError> {
        checkpoint(cx)?;
        if !initial_time_s.is_finite()
            || initial_time_s < 0.0
            || maximum_steps == 0
            || initial_velocity_mps.len() != solid.geometry.nodes_m().len()
            || initial_velocity_mps
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(BrittleTetDynamicError::Invalid(
                "finite full initial velocity, nonnegative time and positive step budget required",
            ));
        }
        let mut densities = Vec::with_capacity(solid.materials.len());
        for (element, material) in solid.materials.iter().enumerate() {
            checkpoint(cx)?;
            densities.push(
                material
                    .reference_density_kg_m3()
                    .map_err(|source| BrittleTetDynamicError::Source { element, source })?,
            );
        }
        let mass = P1ReferenceMass::new(
            solid
                .geometry
                .tetrahedra()
                .iter()
                .zip(solid.geometry.element_volumes_m3())
                .map(|(&nodes, &volume)| (nodes, volume)),
            &densities,
            initial_velocity_mps.len() * 3,
            cx,
        )
        .map_err(SmallStrainTetError::from)?;
        let momentum = mass
            .action(initial_velocity_mps.as_flattened(), cx)
            .map_err(SmallStrainTetError::from)?;
        let mut kinetic = 0.0;
        for (dof, (&v, &p)) in initial_velocity_mps
            .as_flattened()
            .iter()
            .zip(&momentum)
            .enumerate()
        {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            kinetic = v.mul_add(p, kinetic);
        }
        kinetic *= 0.5;
        if !kinetic.is_finite() || kinetic < 0.0 {
            return Err(BrittleTetDynamicError::Invalid(
                "initial consistent kinetic storage must be finite and nonnegative",
            ));
        }
        let retained = BrittleTetSolid {
            geometry: solid.geometry,
            materials: solid.materials,
            bands: Arc::clone(&solid.bands),
            released_faces: Arc::clone(&solid.released_faces),
            budget: solid.budget,
            solution: solid.solution.clone(),
            stored_energy_j: solid.stored_energy_j,
            fracture_dissipation_j: solid.fracture_dissipation_j,
        };
        checkpoint(cx)?;
        Ok(Self {
            solid: retained,
            velocity_mps: initial_velocity_mps.to_vec(),
            reference_density_kg_m3: densities.into(),
            mass: Arc::new(mass),
            kinetic_energy_j: kinetic,
            time_s: initial_time_s,
            maximum_steps,
            accepted_steps: 0,
        })
    }

    #[must_use]
    pub const fn solid(&self) -> &BrittleTetSolid<'a> {
        &self.solid
    }
    #[must_use]
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    #[must_use]
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    #[must_use]
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }
    #[must_use]
    pub const fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_energy_j
    }
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    #[must_use]
    pub const fn accepted_steps(&self) -> usize {
        self.accepted_steps
    }
    #[must_use]
    pub const fn remaining_steps(&self) -> usize {
        self.maximum_steps - self.accepted_steps
    }

    /// Advance one actual physical interval with native backward-Euler inertia.
    ///
    /// Loads and prescribed displacements are absolute endpoint values. Free
    /// bodies need no artificial supports. Reactions include acceleration of
    /// prescribed coordinates and their consistent-mass coupling to free DOFs.
    /// All nonlinear trials borrow the same old fracture history. Any solve,
    /// domain, cancellation or late work refusal leaves the original motion
    /// unchanged; only a complete candidate advances the physical clock.
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetDynamicControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMotionStep<'a>, BrittleTetDynamicError> {
        checkpoint(cx)?;
        if self.accepted_steps >= self.maximum_steps {
            return Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: self.accepted_steps,
                maximum_steps: self.maximum_steps,
            });
        }
        let next_time = self.time_s + dt_s;
        let tolerances = [
            controls.work.material_energy_tolerance_j,
            controls.work.equilibrium_work_tolerance_j,
            controls.work.maximum_endpoint_work_defect_j,
            controls.kinetic_work_tolerance_j,
            controls.maximum_numerical_kinetic_loss_j,
        ];
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !next_time.is_finite()
            || next_time <= self.time_s
            || tolerances.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "representable positive physical time and finite nonnegative work gates required",
            ));
        }
        let old_u = self.solid.solution.displacement_m.as_flattened();
        if prescribed_m.len() > old_u.len() || nodal_forces_n.len() > old_u.len() {
            return Err(BrittleTetDynamicError::Invalid(
                "endpoint load entries exceed the original DOF count",
            ));
        }
        let mut laws = Vec::with_capacity(self.solid.materials.len());
        for material in self.solid.materials {
            checkpoint(cx)?;
            laws.push(material.small_strain_law());
        }
        let problem = SmallStrainTetProblem {
            nodes_m: self.solid.geometry.nodes_m(),
            tetrahedra: self.solid.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.solid.budget,
            settings,
        };
        let step = problem.advance_backward_euler_with_mass(
            old_u,
            &self.velocity_mps,
            &self.solid.solution.element_states,
            &self.mass,
            dt_s,
            cx,
        )?;
        let (material, element_energy, cumulative_loss) =
            self.solid.material_energy(&step.solution, cx)?;
        let new_u = step.solution.displacement_m.as_flattened();
        let mut dead_work = 0.0;
        for &(dof, force) in nodal_forces_n {
            checkpoint(cx)?;
            dead_work = force.mul_add(new_u[dof] - old_u[dof], dead_work);
        }
        let mut support_work = 0.0;
        for (&(dof, _), &reaction) in prescribed_m
            .iter()
            .zip(&step.solution.increments[0].reactions_n)
        {
            checkpoint(cx)?;
            support_work = reaction.mul_add(new_u[dof] - old_u[dof], support_work);
        }
        let external = dead_work + support_work;
        let storage_change = material.stored_after_j - material.stored_before_j;
        let constitutive_defect = material.internal_endpoint_work_j - material.material_work_j;
        let equilibrium_defect =
            (external - material.internal_endpoint_work_j) - step.energy.inertial_work_j;
        let total_defect = ((((external - storage_change) - material.fracture_dissipation_j)
            - step.energy.kinetic_change_j)
            - step.energy.numerical_loss_j)
            - constitutive_defect;
        let energy = BrittleTetDynamicEnergy {
            inertia: step.energy,
            stored_before_j: material.stored_before_j,
            stored_after_j: material.stored_after_j,
            material_work_j: material.material_work_j,
            fracture_dissipation_j: material.fracture_dissipation_j,
            material_balance_error_j: (material.material_work_j - storage_change)
                - material.fracture_dissipation_j,
            internal_endpoint_work_j: material.internal_endpoint_work_j,
            dead_load_endpoint_work_j: dead_work,
            support_endpoint_work_j: support_work,
            equilibrium_work_residual_j: equilibrium_defect,
            constitutive_endpoint_defect_j: constitutive_defect,
            total_balance_error_j: total_defect,
        };
        for (account, defect, tolerance) in [
            (
                "dynamic material work/storage/fracture",
                energy.material_balance_error_j,
                controls.work.material_energy_tolerance_j,
            ),
            (
                "dynamic cumulative fracture loss",
                (cumulative_loss - self.solid.fracture_dissipation_j)
                    - energy.fracture_dissipation_j,
                controls.work.material_energy_tolerance_j,
            ),
            (
                "dynamic support/dead-load equilibrium work",
                equilibrium_defect,
                controls.work.equilibrium_work_tolerance_j,
            ),
            (
                "dynamic constitutive endpoint integration defect",
                constitutive_defect,
                controls.work.maximum_endpoint_work_defect_j,
            ),
            (
                "native kinetic work identity",
                step.energy.balance_error_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "retained kinetic storage",
                step.energy.kinetic_before_j - self.kinetic_energy_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "backward-Euler numerical kinetic loss",
                step.energy.numerical_loss_j,
                controls.maximum_numerical_kinetic_loss_j,
            ),
            (
                "complete brittle kinetic work balance",
                total_defect,
                controls.work.material_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        if [
            material.stored_after_j,
            material.fracture_dissipation_j,
            cumulative_loss,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "finite nonnegative physical storage and fracture loss required",
            ));
        }
        self.solid
            .require_released_faces(&step.solution.displacement_m, cx)?;
        let solid = BrittleTetSolid {
            geometry: self.solid.geometry,
            materials: self.solid.materials,
            bands: Arc::clone(&self.solid.bands),
            released_faces: Arc::clone(&self.solid.released_faces),
            budget: self.solid.budget,
            solution: step.solution,
            stored_energy_j: material.stored_after_j,
            fracture_dissipation_j: cumulative_loss,
        };
        checkpoint(cx)?;
        Ok(BrittleTetMotionStep {
            motion: Self {
                solid,
                velocity_mps: step.velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
                kinetic_energy_j: step.energy.kinetic_after_j,
                time_s: next_time,
                maximum_steps: self.maximum_steps,
                accepted_steps: self.accepted_steps + 1,
            },
            element_energy,
            energy,
        })
    }
}
