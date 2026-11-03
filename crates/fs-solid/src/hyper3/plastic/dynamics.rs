//! Finite plastic motion with consistent reference mass and backward Euler.
//!
//! The endpoint momentum equation and its tangent enter the original hyper3
//! Newton solve. Each trial uses the same accepted plastic history; motion and
//! full constitutive history are returned together only after every gate passes.
//! Free bodies need no artificial restraints. Support reactions include the
//! consistent-mass coupling to prescribed displacement components.
//!
//! The original material, geometry and shared inertia owners account for
//! `external work = stored change + plastic dissipation + signed material
//! endpoint defect + kinetic change + backward-Euler kinetic loss + residual`.
//! The nonnegative kinetic loss is numerical, never material heat. Its explicit
//! budget and the separate signed material-defect budget require physical
//! timestep refinement; a work identity alone does not establish time accuracy.
//!
//! Reference density is an explicit original mass inventory, admitted once and
//! retained as the body deforms. Supplying it does not create a material-source
//! or density-calibration claim. This path holds dead nodal/body loads over one
//! interval. Thermal evolution, follower loads, contact and mass transfer are
//! outside this path; the parent module's plastic-locking limitations remain.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Material};

use super::super::{HyperTetSettings, HyperTetSolution, inertia_error};
use super::{
    FinitePlasticStepControls, FinitePlasticTetEnergy, FinitePlasticTetError,
    FinitePlasticTetSolid, FinitePlasticTetState, checkpoint, energy_gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, P1ReferenceMass};

/// Shared physical kinetic-work account; numerical loss is never plastic heat.
pub use crate::inertia::InertialEnergy;
pub use crate::inertia::{SolidImpulseControls, SolidImpulseReport, SolidMomentum};
pub use crate::inertia::{SolidPointKinematics, SolidPointMobility, TetImpulsePoint};

/// Separate material-work and inertial time-discretization tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePlasticMotionControls {
    /// Original material storage, dissipation, yield and signed-defect gates.
    pub material: FinitePlasticStepControls,
    /// Maximum admitted nonnegative backward-Euler kinetic loss for this step.
    pub maximum_numerical_kinetic_loss_j: f64,
}

/// Complete accepted displacement, velocity, material history and original mass.
#[derive(Debug)]
pub struct FinitePlasticTetMotion<'a, L: FiniteJ2Material = FiniteJ2> {
    solid: FinitePlasticTetSolid<'a, L>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<P1ReferenceMass>,
}

/// A complete proposed interval; the borrowed accepted motion is unchanged.
#[derive(Debug)]
pub struct FinitePlasticMotionStep<'a, L: FiniteJ2Material = FiniteJ2> {
    /// Accepted full next displacement, velocity, source binding and history.
    pub motion: FinitePlasticTetMotion<'a, L>,
    /// Native equilibrium evidence, including inertial support reactions.
    /// `strain_energy_j` remains the material's reduced incremental potential.
    pub equilibrium: HyperTetSolution,
    /// Physical material storage/dissipation and signed endpoint-work defect.
    pub material_energy: FinitePlasticTetEnergy,
    /// Independent kinetic change, inertial work and numerical kinetic loss.
    pub inertia: InertialEnergy,
    /// Sum of the material/momentum work residual and kinetic identity residual.
    pub balance_error_j: f64,
}

/// Instantaneous motion candidate retaining complete plastic constitutive history.
#[derive(Debug)]
pub struct FinitePlasticMotionImpulse<'a, L: FiniteJ2Material = FiniteJ2> {
    pub motion: FinitePlasticTetMotion<'a, L>,
    /// Applied and support impulses, actual momentum change and kinetic work.
    pub impulse: SolidImpulseReport,
}

impl<'a, L: FiniteJ2Material> FinitePlasticTetMotion<'a, L> {
    /// Bind an admitted static configuration to initial velocity and fixed mass.
    ///
    /// The supplied solid's complete history and original source owners are
    /// retained. Density has one value per original tetrahedron and never uses
    /// current volume; velocity has one finite vector per original node.
    ///
    /// # Errors
    /// Invalid density/velocity, unrepresentable reference mass, or cancellation.
    pub fn try_new(
        solid: FinitePlasticTetSolid<'a, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePlasticTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.state.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity vector per original node required",
            ));
        }
        let mass = P1ReferenceMass::new(
            solid
                .state
                .elements
                .iter()
                .map(|element| (element.nodes, element.volume)),
            reference_density_kg_m3,
            solid.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(reference_density_kg_m3),
            reference_nodal_density_kg_m3: None,
            mass: Arc::new(mass),
        })
    }

    /// Bind positive nodal reference density to the SAME native momentum owner.
    ///
    /// Density is interpolated with the original P1 material coordinates.
    /// The consistent matrix integrates rho*N_a*N_b exactly, rather than
    /// replacing the field by its cell mean. The immutable nodal values and
    /// mass operator survive every accepted deformation and plastic step.
    /// This is an authored reference inventory, not a source/calibration claim.
    ///
    /// # Errors
    /// Invalid nodal density/velocity, unrepresentable mass, or cancellation.
    pub fn try_new_with_nodal_density(
        solid: FinitePlasticTetSolid<'a, L>,
        reference_nodal_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePlasticTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.state.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity vector per original node required",
            ));
        }
        let (mass, means) = P1ReferenceMass::from_nodal_density(
            solid
                .state
                .elements
                .iter()
                .map(|element| (element.nodes, element.volume)),
            reference_nodal_density_kg_m3,
            solid.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(means),
            reference_nodal_density_kg_m3: Some(Arc::from(reference_nodal_density_kg_m3)),
            mass: Arc::new(mass),
        })
    }

    /// Accepted finite deformation, source binding and full plastic history.
    pub fn solid(&self) -> &FinitePlasticTetSolid<'a, L> {
        &self.solid
    }

    /// Accepted velocity in original node order, m/s.
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }

    /// Reference cell means, in original tetrahedron order. The cellwise
    /// constructor retains its exact supplied values. A nodal-density body
    /// also retains its complete field through `reference_nodal_density_kg_m3`;
    /// these means alone do not reproduce its mass operator or mass distribution.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Complete original nodal field, or None for the cellwise-density model.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_nodal_density_kg_m3.as_deref()
    }

    /// Total original mass, unchanged by current deformation or plastic flow.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Current physical material-point position and work-conjugate velocity.
    pub fn point_kinematics(
        &self,
        point: TetImpulsePoint,
        cx: &Cx<'_>,
    ) -> Result<SolidPointKinematics, FinitePlasticTetError> {
        let state = &self.solid.state;
        crate::inertia::point::kinematics(
            point,
            state.tetrahedra,
            state.nodes_m,
            &state.displacement_m,
            &self.velocity_mps,
            cx,
        )
        .map_err(|error| FinitePlasticTetError::from(inertia_error(error)))
    }

    /// Supported directional inverse mass at a material point [1/kg].
    /// Uses the original consistent mass and a unit world direction. No material
    /// evaluation, history clone, dense inverse or lumped mass is needed.
    pub fn point_directional_inverse_mass_per_kg(
        &self,
        held: &[usize],
        point: TetImpulsePoint,
        direction: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobility, FinitePlasticTetError> {
        crate::inertia::point::directional_mobility(
            &self.mass,
            self.solid.state.tetrahedra,
            held,
            point,
            direction,
            controls,
            cx,
        )
        .map_err(|error| FinitePlasticTetError::from(inertia_error(error)))
    }

    /// Apply a world impulse at the current material point through B^T.
    /// This is the exact transpose of `point_kinematics`' velocity trace, so
    /// point work and torque use the same position and supported mass response.
    pub fn point_impulse_candidate(
        &self,
        held: &[usize],
        point: TetImpulsePoint,
        impulse_world: [f64; 3],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePlasticMotionImpulse<'a, L>, FinitePlasticTetError> {
        let entries = crate::inertia::point::distribution(
            point,
            self.solid.state.tetrahedra,
            self.velocity_mps.len(),
            impulse_world,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Apply an instantaneous nodal impulse through the original consistent mass.
    ///
    /// Impulses are full nodal DOFs and world components in kg m/s; repeated
    /// entries accumulate. Held DOFs retain their existing velocity, including
    /// nonzero moving supports. Support impulses and their work are reported.
    /// Current deformed positions define angular momentum about `origin_m`.
    ///
    /// The event has zero duration and changes only velocity. Deformation,
    /// plastic metrics, accumulated flow, first Piola stress, stored energy,
    /// source laws and original density are retained exactly. Plastic flow
    /// resumes through the native `advance` owner from the changed velocity;
    /// an impulse itself is not plastic dissipation. The caller supplies the
    /// contact/restitution law and any equal-and-opposite partner impulse.
    ///
    /// # Errors
    /// Invalid input or controls, mass-solve/work/momentum refusal, or cancellation.
    /// No part of the accepted motion changes on failure.
    pub fn impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePlasticMotionImpulse<'a, L>, FinitePlasticTetError> {
        checkpoint(cx)?;
        let state = &self.solid.state;
        let mut positions_m = Vec::with_capacity(self.velocity_mps.len());
        for (node, (reference, displacement)) in
            state.nodes_m.iter().zip(&state.displacement_m).enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            positions_m.push(core::array::from_fn(|axis| {
                reference[axis] + displacement[axis]
            }));
        }
        let (velocity_mps, impulse) = self
            .mass
            .impulse_candidate(
                &positions_m,
                &self.velocity_mps,
                held_velocity_dofs,
                nodal_impulses_kg_m_s,
                origin_m,
                controls,
                cx,
            )
            .map_err(inertia_error)?;
        let solid = FinitePlasticTetSolid {
            material: self.solid.material,
            state: FinitePlasticTetState {
                nodes_m: state.nodes_m,
                tetrahedra: state.tetrahedra,
                elements: Arc::clone(&state.elements),
                budget: state.budget,
                displacement_m: state.displacement_m.clone(),
                element_deformation: state.element_deformation.clone(),
                element_states: state.element_states.clone(),
                element_first_piola_pa: state.element_first_piola_pa.clone(),
                stored_energy_j: state.stored_energy_j,
            },
        };
        checkpoint(cx)?;
        Ok(FinitePlasticMotionImpulse {
            motion: Self {
                solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                reference_nodal_density_kg_m3: self
                    .reference_nodal_density_kg_m3
                    .as_ref()
                    .map(Arc::clone),
                mass: Arc::clone(&self.mass),
            },
            impulse,
        })
    }

    /// Project an acceleration field with this body's actual consistent mass.
    ///
    /// Supply one world acceleration vector [m/s²] per original node. P1
    /// interpolation and the retained reference density give the exact nodal
    /// force projection M*a [N], including uniform gravity on a graded body.
    /// Returned nonzero DOF/force pairs feed `advance` as held nodal loads.
    /// Do not also add the same acceleration as a body-force density. This
    /// assembles a load; position-dependent field evolution and its potential
    /// or tangent remain with the external field owner.
    ///
    /// # Errors
    /// Wrong shape, nonfinite acceleration/force, or cancellation.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FinitePlasticTetError> {
        checkpoint(cx)?;
        if acceleration_m_s2.len() != self.velocity_mps.len()
            || acceleration_m_s2
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite acceleration vector per original node required",
            ));
        }
        let forces = self
            .mass
            .action(acceleration_m_s2.as_flattened(), cx)
            .map_err(inertia_error)?;
        let mut loads = Vec::with_capacity(forces.len());
        for (dof, force) in forces.into_iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if force != 0.0 {
                loads.push((dof, force));
            }
        }
        checkpoint(cx)?;
        Ok(loads)
    }

    /// Advance one physical interval with held dead loads and absolute supports.
    ///
    /// Material history stays frozen throughout the one native endpoint solve.
    /// Force convergence uses physical free-row forces, including inertia;
    /// stiffness times rigid travel cannot hide missing momentum. The accepted
    /// motion is borrowed and no part changes on any error or budget rejection.
    ///
    /// # Errors
    /// Invalid timestep/controls, native solve/material refusal, numerical-loss
    /// budget exceedance, independent work-balance failure, or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        settings: HyperTetSettings,
        controls: FinitePlasticMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePlasticMotionStep<'a, L>, FinitePlasticTetError> {
        checkpoint(cx)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(invalid(
                "finite nonnegative numerical kinetic-loss budget required",
            ));
        }
        let inertia = BackwardEulerInertia::new(
            &self.mass,
            self.solid.state.displacement_m.as_flattened(),
            self.velocity_mps.as_flattened(),
            dt_s,
            cx,
        )
        .map_err(inertia_error)?;
        let (step, motion) = self.solid.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let (velocity_mps, inertia) = motion.expect("inertial endpoint returns native motion");
        energy_gate(
            "backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        let balance_error_j = step.energy.equilibrium_work_residual_j + inertia.balance_error_j;
        energy_gate(
            "total finite inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        let motion = Self {
            solid: step.solid,
            velocity_mps,
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            reference_nodal_density_kg_m3: self
                .reference_nodal_density_kg_m3
                .as_ref()
                .map(Arc::clone),
            mass: Arc::clone(&self.mass),
        };
        Ok(FinitePlasticMotionStep {
            motion,
            equilibrium: step.equilibrium,
            material_energy: step.energy,
            inertia,
            balance_error_j,
        })
    }
}
