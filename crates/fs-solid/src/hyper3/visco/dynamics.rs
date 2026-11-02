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
            mass: Arc::new(mass),
        })
    }

    pub fn solid(&self) -> &FiniteViscoTetSolid<'a> {
        &self.solid
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
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
                mass: Arc::clone(&self.mass),
            },
            equilibrium: step.equilibrium,
            material_energy: step.energy,
            inertia,
            balance_error_j,
        })
    }
}
