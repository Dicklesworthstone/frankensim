//! Native momentum evolution with complete finite-strain viscoelastic memory.
//!
//! The original consistent reference mass and backward-Euler inertia enter the
//! SAME hyper3 Newton solve as the frozen-history material potential. Physical
//! relaxation heat, signed material endpoint error and backward-Euler numerical
//! kinetic loss remain separate. Free bodies require no artificial supports.
//! Density is explicit original mass per reference volume and never follows
//! current deformation. This is an isothermal, closed-mass, held-dead-load path;
//! contact, follower loads and thermal evolution retain their separate owners.

use std::sync::Arc;

use fs_exec::Cx;

use super::super::{HyperTetSettings, HyperTetSolution, inertia_error};
use super::{
    FiniteViscoStepControls, FiniteViscoTetEnergy, FiniteViscoTetError, FiniteViscoTetSolid,
    checkpoint, energy_gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, P1ReferenceMass};

pub use crate::inertia::InertialEnergy;
pub use crate::inertia::{SolidImpulseControls, SolidImpulseReport, SolidMomentum};
pub use crate::inertia::{
    SolidPointKinematics, SolidPointMobility, SolidPointMobilityMatrix, TetImpulsePoint,
};

/// Independent budgets for material work and numerical kinetic damping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteViscoMotionControls {
    pub material: FiniteViscoStepControls,
    /// Nonnegative BE kinetic loss is a discretization error, never heat [J].
    pub maximum_numerical_kinetic_loss_j: f64,
}

/// Accepted full material history, motion and original reference inventory.
#[derive(Debug)]
pub struct FiniteViscoTetMotion<'a> {
    solid: FiniteViscoTetSolid<'a>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<P1ReferenceMass>,
}

/// One complete candidate; previous displacement, velocity and memory are intact.
#[derive(Debug)]
pub struct FiniteViscoMotionStep<'a> {
    pub motion: FiniteViscoTetMotion<'a>,
    /// Includes native inertial support reactions; the potential is material only.
    pub equilibrium: HyperTetSolution,
    pub material_energy: FiniteViscoTetEnergy,
    pub inertia: InertialEnergy,
    /// Momentum work residual plus the independent kinetic identity residual.
    pub balance_error_j: f64,
}

/// Instantaneous velocity change with the complete original material memory.
#[derive(Debug)]
pub struct FiniteViscoMotionImpulse<'a> {
    pub motion: FiniteViscoTetMotion<'a>,
    /// Applied and support impulses, actual momentum change and kinetic work.
    pub impulse: SolidImpulseReport,
}

impl<'a> FiniteViscoTetMotion<'a> {
    /// Bind an accepted solid to original density and initial velocity. Complete
    /// relaxed/deformed histories and elapsed material time are retained.
    pub fn try_new(
        solid: FiniteViscoTetSolid<'a>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteViscoTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity per original node required",
            ));
        }
        let mass = P1ReferenceMass::new(
            solid
                .elements
                .iter()
                .map(|element| (element.nodes, element.volume)),
            reference_density_kg_m3,
            solid.displacement_m.as_flattened().len(),
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

    /// Bind positive P1 nodal reference density to the original momentum owner.
    ///
    /// The shared consistent mass integrates `rho*N_a*N_b` exactly. Cell means
    /// remain inventory diagnostics; they do not replace the graded mass matrix.
    /// The accepted solid's deformation, elapsed time and every Maxwell branch
    /// history are retained. Density is an authored reference inventory, not a
    /// new material-source claim, and stays fixed through subsequent motion.
    ///
    /// # Errors
    /// Invalid density/velocity, unrepresentable reference mass, or cancellation.
    pub fn try_new_with_nodal_density(
        solid: FiniteViscoTetSolid<'a>,
        reference_nodal_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteViscoTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity per original node required",
            ));
        }
        let (mass, means) = P1ReferenceMass::from_nodal_density(
            solid
                .elements
                .iter()
                .map(|element| (element.nodes, element.volume)),
            reference_nodal_density_kg_m3,
            solid.displacement_m.as_flattened().len(),
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

    pub fn solid(&self) -> &FiniteViscoTetSolid<'a> {
        &self.solid
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    /// Reference cell means in original tetrahedron order. Cellwise bindings
    /// retain their supplied values; graded bindings retain their complete
    /// field separately because means alone do not reproduce their inertia.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    /// Complete original nodal field, or None for the cellwise-density model.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_nodal_density_kg_m3.as_deref()
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Current physical material-point position and work-conjugate velocity.
    pub fn point_kinematics(
        &self,
        point: TetImpulsePoint,
        cx: &Cx<'_>,
    ) -> Result<SolidPointKinematics, FiniteViscoTetError> {
        crate::inertia::point::kinematics(
            point,
            self.solid.tetrahedra,
            self.solid.nodes_m,
            &self.solid.displacement_m,
            &self.velocity_mps,
            cx,
        )
        .map_err(|error| FiniteViscoTetError::from(inertia_error(error)))
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
    ) -> Result<SolidPointMobility, FiniteViscoTetError> {
        crate::inertia::point::directional_mobility(
            &self.mass,
            self.solid.tetrahedra,
            held,
            point,
            direction,
            controls,
            cx,
        )
        .map_err(|error| FiniteViscoTetError::from(inertia_error(error)))
    }

    /// Complete simultaneous cross-point/direction mobility [1/kg].
    /// Each row uses a unit world direction and the original supported mass.
    /// Only sparse point rows, the bounded result and one full response at a
    /// time are retained; no constitutive history is copied by this query.
    pub fn point_mobility_matrix_per_kg(
        &self,
        held: &[usize],
        rows: &[(TetImpulsePoint, [f64; 3])],
        maximum_rows: usize,
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<SolidPointMobilityMatrix, FiniteViscoTetError> {
        crate::inertia::point::mobility_matrix(
            &self.mass,
            self.solid.tetrahedra,
            held,
            rows,
            maximum_rows,
            controls,
            cx,
        )
        .map_err(|error| FiniteViscoTetError::from(inertia_error(error)))
    }

    /// Apply all material-point impulses in one native velocity update.
    /// The original mass accounts for every cross term, support reaction and
    /// signed work; compensated accumulation retains small net nodal impulses.
    pub fn point_impulses_candidate(
        &self,
        held: &[usize],
        impulses: &[(TetImpulsePoint, [f64; 3])],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoMotionImpulse<'a>, FiniteViscoTetError> {
        let entries = crate::inertia::point::distributions(
            self.solid.tetrahedra,
            self.velocity_mps.len(),
            impulses,
            controls,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Apply a world impulse through the transpose of the material-point trace.
    /// Current point work and torque, original mass and all relaxation memory
    /// are retained by the same native nodal impulse account.
    pub fn point_impulse_candidate(
        &self,
        held: &[usize],
        point: TetImpulsePoint,
        impulse_world: [f64; 3],
        origin: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoMotionImpulse<'a>, FiniteViscoTetError> {
        let entries = crate::inertia::point::distribution(
            point,
            self.solid.tetrahedra,
            self.velocity_mps.len(),
            impulse_world,
            cx,
        )
        .map_err(inertia_error)?;
        self.impulse_candidate(held, &entries, origin, controls, cx)
    }

    /// Apply an instantaneous nodal impulse through the original consistent mass.
    ///
    /// Each impulse is a full nodal DOF and a world component in kg m/s;
    /// repeated impulse entries accumulate. Held DOFs keep their current
    /// velocity, including nonzero moving-support velocity. Their reaction
    /// impulses and work remain in the returned physical balance.
    ///
    /// The current deformed positions define angular momentum about `origin_m`.
    /// The event has zero duration: displacement, every Maxwell branch, stress,
    /// stored energy, elapsed material time and original density are unchanged.
    /// Subsequent `advance` calls evolve this retained memory from the new
    /// velocity. The caller owns the contact/restitution law and any opposite
    /// body's impulse; this operation adds no relaxation heat.
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
    ) -> Result<FiniteViscoMotionImpulse<'a>, FiniteViscoTetError> {
        checkpoint(cx)?;
        let mut positions_m = Vec::with_capacity(self.velocity_mps.len());
        for (node, (reference, displacement)) in self
            .solid
            .nodes_m
            .iter()
            .zip(&self.solid.displacement_m)
            .enumerate()
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
        let solid = FiniteViscoTetSolid {
            nodes_m: self.solid.nodes_m,
            tetrahedra: self.solid.tetrahedra,
            material: self.solid.material,
            elements: Arc::clone(&self.solid.elements),
            budget: self.solid.budget,
            displacement_m: self.solid.displacement_m.clone(),
            states: self.solid.states.clone(),
            first_piola_pa: self.solid.first_piola_pa.clone(),
            stored_energy_j: self.solid.stored_energy_j,
            time_s: self.solid.time_s,
        };
        checkpoint(cx)?;
        Ok(FiniteViscoMotionImpulse {
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

    /// Project a P1 acceleration field through the actual reference mass, M*a.
    ///
    /// One world acceleration [m/s²] per original node gives the exact nodal
    /// load [N] for the retained density, including uniform gravity on a graded
    /// body. Pass these held loads to `advance` without adding the same force
    /// again as a body-force density. Evolving field values, potential and
    /// tangent remain the external field owner's responsibility.
    ///
    /// # Errors
    /// Wrong shape, nonfinite acceleration/force, or cancellation.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteViscoTetError> {
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

    /// Advance one positive physical interval with absolute supports and held
    /// dead loads. Native inertia supplies the free-row momentum stopping gate,
    /// consistent mass tangent and support coupling. Every material trial uses
    /// the same accepted history and the SAME interval as the momentum equation.
    ///
    /// All candidates are transactional, including post-solve kinetic budgets.
    /// Reduce the physical timestep when numerical-loss or work gates refuse;
    /// a closed work identity alone does not establish time accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        controls: FiniteViscoMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoMotionStep<'a>, FiniteViscoTetError> {
        checkpoint(cx)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(invalid(
                "finite nonnegative numerical kinetic-loss allowance required",
            ));
        }
        let inertia = BackwardEulerInertia::new(
            &self.mass,
            self.solid.displacement_m.as_flattened(),
            self.velocity_mps.as_flattened(),
            dt_s,
            cx,
        )
        .map_err(inertia_error)?;
        let (step, motion) = self.solid.advance_loaded(
            dt_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let (velocity_mps, inertia) = motion.expect("inertial endpoint supplies native motion");
        energy_gate(
            "backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        let balance_error_j = step.energy.equilibrium_work_residual_j + inertia.balance_error_j;
        energy_gate(
            "total viscoelastic inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteViscoMotionStep {
            motion: Self {
                solid: step.solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                reference_nodal_density_kg_m3: self
                    .reference_nodal_density_kg_m3
                    .as_ref()
                    .map(Arc::clone),
                mass: Arc::clone(&self.mass),
            },
            equilibrium: step.equilibrium,
            material_energy: step.energy,
            inertia,
            balance_error_j,
        })
    }
}
