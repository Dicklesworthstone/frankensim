//! Finite porous motion through the original momentum and compaction owners.
//!
//! The consistent P1 reference mass and backward-Euler tangent enter the same
//! native finite-element Newton solve as the porous incremental potential.
//! Every trial starts from the accepted non-isochoric plastic metric. A complete
//! endpoint retains displacement, velocity, plastic volume, hardening and mass
//! only after all material, momentum and work gates pass.
//!
//! Original mass is independent of plastic volume: current density is rho0/J,
//! never rho0/Jp. Physical plastic dissipation, signed material endpoint-work
//! defect and backward-Euler numerical kinetic loss have separate accounts and
//! budgets. Numerical kinetic loss is not heat; timestep refinement is still
//! necessary for time accuracy even when the work identity closes.
//!
//! Density is an explicit original mass inventory. `FinitePorousMaterial`
//! currently exposes no sourced density; borrowing a source-resolved law keeps
//! its card and context but makes no density-calibration claim. Conditions stay
//! fixed: this path adds no evolving porosity, densification, thermal coupling,
//! contact, follower loads, mass transfer or pressure stabilization.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::plastic::finite::porous::{FinitePorous, FinitePorousMaterial};

use super::{FinitePorousTetError, FinitePorousTetSolid};
use crate::hyper3::plastic::{FinitePlasticTetEnergy, checkpoint, energy_gate, invalid};
use crate::hyper3::{HyperTetSettings, HyperTetSolution, inertia_error};
use crate::inertia::{BackwardEulerInertia, P1ReferenceMass};

/// Original separate material and numerical-inertia admission controls.
pub use crate::hyper3::plastic::dynamics::FinitePlasticMotionControls;
/// Original physical kinetic-work account, separate from material heat.
pub use crate::inertia::InertialEnergy;

/// Accepted porous configuration, complete history, velocity and original mass.
#[derive(Debug)]
pub struct FinitePorousTetMotion<'a, L: FinitePorousMaterial = FinitePorous> {
    solid: FinitePorousTetSolid<'a, L>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<P1ReferenceMass>,
}

/// Complete admitted physical interval; the previously accepted body is intact.
#[derive(Debug)]
pub struct FinitePorousMotionStep<'a, L: FinitePorousMaterial = FinitePorous> {
    /// Complete next state, including non-isochoric plastic metric and mass.
    pub motion: FinitePorousTetMotion<'a, L>,
    /// Native solve evidence and support reactions including consistent inertia.
    /// `strain_energy_j` remains the material's reduced incremental potential.
    pub equilibrium: HyperTetSolution,
    /// Physical storage/yield loss and signed material integration defect.
    pub material_energy: FinitePlasticTetEnergy,
    /// Kinetic change, momentum work and separate numerical kinetic loss.
    pub inertia: InertialEnergy,
    /// Combined material/momentum and independent kinetic-identity residual.
    pub balance_error_j: f64,
}

fn validate_nodal_vectors(
    node_count: usize,
    values: &[[f64; 3]],
    invalid_input: &'static str,
    cx: &Cx<'_>,
) -> Result<(), FinitePorousTetError> {
    checkpoint(cx)?;
    if values.len() != node_count {
        return Err(invalid(invalid_input).into());
    }
    for (node, vector) in values.iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        if vector.iter().any(|value| !value.is_finite()) {
            return Err(invalid(invalid_input).into());
        }
    }
    Ok(())
}

impl<'a, L: FinitePorousMaterial> FinitePorousTetMotion<'a, L> {
    /// Bind a virgin or previously compacted solid to fixed mass and velocity.
    ///
    /// One finite positive density per original tetrahedron defines reference
    /// mass, regardless of accepted total/plastic volume. The solid retains its
    /// original material owners and all accepted histories. The original
    /// consistent mass matrix is prepared once and shared by later endpoints.
    ///
    /// # Errors
    /// Invalid density/velocity, unrepresentable reference mass or cancellation.
    pub fn try_new(
        solid: FinitePorousTetSolid<'a, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePorousTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.state.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(
                invalid("one finite initial velocity vector per original node required").into(),
            );
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

    /// Bind positive nodal density through exact P1 reference-mass integration.
    ///
    /// The native mass owner integrates `rho_ref * N_a * N_b`, retaining the
    /// original field and its full matrix across every compaction step. Cell
    /// means are available for volume/density queries but do not replace this
    /// mass distribution or its coupling to moving supports.
    ///
    /// # Errors
    /// Invalid nodal density/velocity, unrepresentable mass, or cancellation.
    pub fn try_new_with_nodal_density(
        solid: FinitePorousTetSolid<'a, L>,
        reference_nodal_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePorousTetError> {
        validate_nodal_vectors(
            solid.state.nodes_m.len(),
            velocity_mps,
            "one finite initial velocity vector per original node required",
            cx,
        )?;
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

    /// Accepted finite configuration, source binding and full compaction history.
    pub fn solid(&self) -> &FinitePorousTetSolid<'a, L> {
        &self.solid
    }

    /// Accepted velocity in original node order, m/s.
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }

    /// Original cell density, or exact cell means of a retained P1 nodal field.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Complete original nodal density field, if that constructor was selected.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_nodal_density_kg_m3.as_deref()
    }

    /// Original mass, unchanged by deformation or permanent plastic compaction.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Density over each current cell: retained reference density divided by J.
    /// The plastic Jacobian never rescales the mass inventory.
    ///
    /// # Errors
    /// Unrepresentable current density or cancellation.
    pub fn current_density_kg_m3(&self, cx: &Cx<'_>) -> Result<Vec<f64>, FinitePorousTetError> {
        self.solid
            .current_density_kg_m3(&self.reference_density_kg_m3, cx)
    }

    /// Project a nodal acceleration field through this body's consistent mass.
    ///
    /// Returned original DOF/force pairs are `M*a` in newtons, including uniform
    /// gravity on graded material. Feed them to `advance` as held nodal loads;
    /// do not also apply the same acceleration as a body-force density. Evolving
    /// external-field potentials and tangents remain with their owning model.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FinitePorousTetError> {
        validate_nodal_vectors(
            self.velocity_mps.len(),
            acceleration_m_s2,
            "one finite acceleration vector per original node required",
            cx,
        )?;
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

    /// Advance one physical interval under held dead loads and absolute supports.
    ///
    /// Body force is force per original reference volume. Free bodies require
    /// no artificial restraints. Prescribed motion participates in the full
    /// consistent mass action and support reactions. Every trial borrows the
    /// same accepted history; any failed solve or gate leaves `self` unchanged.
    ///
    /// # Errors
    /// Invalid timestep/controls, native solve/material refusal, work or kinetic
    /// loss budget exceedance, unrepresentable motion or cancellation.
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
    ) -> Result<FinitePorousMotionStep<'a, L>, FinitePorousTetError> {
        checkpoint(cx)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(
                invalid("finite nonnegative numerical kinetic-loss budget required").into(),
            );
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
            "porous backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        let balance_error_j = step.energy.equilibrium_work_residual_j + inertia.balance_error_j;
        energy_gate(
            "total finite porous inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FinitePorousMotionStep {
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
