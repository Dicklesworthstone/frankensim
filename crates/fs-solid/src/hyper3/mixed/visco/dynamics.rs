//! Consistent quadratic inertia and retained mixed viscoelastic motion.
//!
//! The native backward-Euler momentum equation enters the SAME P2/P1 solve as
//! the frozen-history material. Pressure has no mass. Original reference
//! density and the signed consistent P2 mass are immutable; no row-sum lumping
//! or deformation-dependent mass update occurs. Material loss, signed material
//! endpoint-work defect and numerical kinetic loss remain separate.
//!
//! Exact pressure constrains endpoint volume weakly; velocity is the native
//! backward-Euler displacement chord, not a separate tangent projection. A
//! closed work identity does not establish time accuracy. Numerical-loss and
//! work allowances can reject the complete candidate without changing history.
//! This first-order path permits budgeted numerical damping; it is not the
//! instrument-render doctrine's zero-algorithmic-damping/discrete-gradient integrator.

use std::sync::Arc;

use fs_exec::Cx;

use super::super::{MixedUpSettings, MixedUpSolution, inertia_error};
use super::{
    MixedViscoTetEnergy, MixedViscoTetError, MixedViscoTetSolid, checkpoint, gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, ReferenceMass};

pub use crate::hyper3::visco::dynamics::FiniteViscoMotionControls;
pub use crate::inertia::InertialEnergy;

/// Immutable accepted displacement, continuous pressure, memory and velocity.
#[derive(Debug)]
pub struct MixedViscoTetMotion<'a> {
    solid: MixedViscoTetSolid<'a>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<ReferenceMass>,
}

/// Complete proposed motion; the previous endpoint remains unchanged.
#[derive(Debug)]
pub struct MixedViscoMotionStep<'a> {
    pub motion: MixedViscoTetMotion<'a>,
    /// Support reactions include native inertia; potential_j is material only.
    pub equilibrium: MixedUpSolution,
    pub material_energy: MixedViscoTetEnergy,
    pub inertia: InertialEnergy,
    /// Momentum work residual plus the independent kinetic identity residual.
    pub balance_error_j: f64,
}

impl<'a> MixedViscoTetMotion<'a> {
    /// Bind an accepted solid to original reference densities and initial
    /// velocities. Deformation, pressure, relaxation history and time persist.
    /// Mass uses the admitted cached quadrature: exact for affine P2 tetrahedra;
    /// curved geometry retains the caller's explicit quadrature approximation.
    pub fn try_new(
        solid: MixedViscoTetSolid<'a>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        maximum_mass_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedViscoTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.displacement_m.len()
            || velocity_mps.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity per original quadratic node required",
            )
            .into());
        }
        let mass = ReferenceMass::from_tet10_quadrature(
            solid.mesh.elements.iter().map(|element| {
                (
                    element.connectivity,
                    element.points.iter().map(|q| (q.shapes, q.weight_m3)),
                )
            }),
            reference_density_kg_m3,
            solid.displacement_m.as_flattened().len(),
            maximum_mass_assembly_bytes,
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

    pub fn solid(&self) -> &MixedViscoTetSolid<'a> {
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

    /// Advance one positive interval using absolute support targets and held
    /// dead loads. Native inertia supplies M/dt² only in the displacement block;
    /// the original pressure equations and frozen material histories are kept.
    /// The material interval equals the momentum interval for every trial.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        settings: MixedUpSettings,
        controls: FiniteViscoMotionControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoMotionStep<'a>, MixedViscoTetError> {
        checkpoint(cx)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(
                invalid("finite nonnegative numerical kinetic-loss allowance required").into(),
            );
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
            pressure_gauge,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let (velocity_mps, inertia) = motion.expect("inertial endpoint supplies native motion");
        gate(
            "backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        let balance_error_j = step.energy.equilibrium_work_residual_j + inertia.balance_error_j;
        gate(
            "total mixed viscoelastic inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(MixedViscoMotionStep {
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
