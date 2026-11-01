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
    FinitePlasticTetSolid, checkpoint, energy_gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, P1ReferenceMass};

/// Shared physical kinetic-work account; numerical loss is never plastic heat.
pub use crate::inertia::InertialEnergy;

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

    /// Exact originally admitted reference density, in tetrahedron order.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Total original mass, unchanged by current deformation or plastic flow.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
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
