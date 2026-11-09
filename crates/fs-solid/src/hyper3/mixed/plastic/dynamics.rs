//! Consistent quadratic inertia with complete retained mixed inelastic history.
//!
//! Backward Euler adds the original reference mass only to the displacement
//! block of the same native P2/P1 solve. Continuous pressure has no inertia.
//! Every Newton trial borrows the accepted plastic, creep and recovery state;
//! the complete candidate publishes motion and history together after the
//! original material, momentum and work gates. A kinetic material's frozen
//! duration must equal the physical momentum interval.
//!
//! Physical inelastic heat, signed material endpoint-work error and numerical
//! kinetic loss remain distinct. Backward Euler permits numerical damping
//! under its explicit allowance; the work identity is not a time-accuracy
//! claim, and physical timestep refinement remains necessary. Exact
//! incompressibility enforces accepted endpoint volume weakly, while
//! velocity is the displacement chord, without a separate tangent projection.
//!
//! The original reference density and signed consistent P2 mass are retained
//! across deformation. No mass lumping, current-volume update, thermal law,
//! contact response or follower load is introduced by this motion wrapper.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::plastic::finite::FiniteJ2;
use fs_material::plastic::finite::mixed::MixedFiniteInelasticMaterial;

use super::super::{MixedUpSettings, MixedUpSolution, inertia_error};
use super::{
    MixedPlasticTetEnergy, MixedPlasticTetError, MixedPlasticTetSolid, checkpoint, gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, ReferenceMass};

/// Independent material-work and numerical kinetic-loss allowances.
pub use crate::hyper3::plastic::dynamics::FinitePlasticMotionControls;
/// Original-mass kinetic and momentum accounts; neither advances material state.
pub use crate::inertia::{InertialEnergy, SolidMomentum};

/// Immutable accepted displacement, pressure, complete material memory and velocity.
#[derive(Debug)]
pub struct MixedPlasticTetMotion<'a, L: MixedFiniteInelasticMaterial = FiniteJ2> {
    solid: MixedPlasticTetSolid<'a, L>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<ReferenceMass>,
}

/// A complete proposed interval; any refusal leaves the borrowed motion unchanged.
#[derive(Debug)]
pub struct MixedPlasticMotionStep<'a, L: MixedFiniteInelasticMaterial = FiniteJ2> {
    /// Complete next displacement, pressure, velocity and original material state.
    pub motion: MixedPlasticTetMotion<'a, L>,
    /// Native support reactions include inertia; potential_j is material only.
    pub equilibrium: MixedUpSolution,
    /// Physical storage and inelastic heat, with separate signed endpoint error.
    pub material_energy: MixedPlasticTetEnergy,
    /// Kinetic change, inertial work and backward-Euler numerical kinetic loss.
    pub inertia: InertialEnergy,
    /// Momentum work residual plus the independent kinetic identity residual.
    pub balance_error_j: f64,
}

impl<'a, L: MixedFiniteInelasticMaterial> MixedPlasticTetMotion<'a, L> {
    /// Bind an accepted mixed configuration to its fixed original mass inventory.
    ///
    /// One positive finite density [kg/m³] belongs to each original tetrahedron;
    /// one finite velocity belongs to each quadratic node. The existing cached
    /// quadrature integrates the signed P2 shape products without lumping. It
    /// is exact for affine geometry and retains the declared quadrature
    /// approximation on curved elements. All source bindings and history persist.
    ///
    /// # Errors
    /// Invalid density/velocity, mass assembly budget, unrepresentable mass,
    /// or cancellation refuses before returning any motion state.
    pub fn try_new(
        solid: MixedPlasticTetSolid<'a, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        maximum_mass_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedPlasticTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity vector per original quadratic node required",
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

    /// Accepted mixed configuration and complete original constitutive history.
    pub fn solid(&self) -> &MixedPlasticTetSolid<'a, L> {
        &self.solid
    }

    /// Accepted world velocities, in original quadratic-node order [m/s].
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }

    /// Immutable original reference densities, in tetrahedron order [kg/m³].
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Original inventory [kg], independent of accepted deformation and pressure.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Original finite compliance or explicitly selected exact volume constraint.
    pub fn volumetric_treatment(&self) -> VolumetricTreatment {
        self.solid.volumetric
    }

    /// Physical momentum at the accepted deformed geometry about a finite origin.
    /// Reuses the original mass and compensated inventory without evolving motion.
    pub fn momentum(
        &self,
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<SolidMomentum, MixedPlasticTetError> {
        let account = self
            .mass
            .constant_velocity_drift(
                self.solid.nodes_m,
                &self.solid.displacement_m,
                &self.solid.displacement_m,
                &self.velocity_mps,
                origin_m,
                cx,
            )
            .map_err(inertia_error)?;
        Ok(account.before)
    }

    /// Physical kinetic energy [J] through the same original consistent mass.
    pub fn kinetic_energy_j(&self, cx: &Cx<'_>) -> Result<f64, MixedPlasticTetError> {
        let account = self
            .mass
            .constant_velocity_drift(
                self.solid.nodes_m,
                &self.solid.displacement_m,
                &self.solid.displacement_m,
                &self.velocity_mps,
                [0.0; 3],
                cx,
            )
            .map_err(inertia_error)?;
        Ok(account.kinetic_j)
    }

    /// Project a P2 reference acceleration [m/s²] through the retained mass.
    ///
    /// Return original nodal dead loads [N], including signed corner/edge
    /// weights and supported components. Feed them to `advance` once; adding
    /// the same acceleration again as body-force density would double the load.
    /// This does not introduce a position-dependent gravity potential.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, MixedPlasticTetError> {
        checkpoint(cx)?;
        if acceleration_m_s2.len() != self.solid.nodes_m.len()
            || acceleration_m_s2
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(
                invalid("one finite acceleration vector per original node required").into(),
            );
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
    /// Original inertia enters only the displacement block of the same mixed
    /// solve. Pressure, full constitutive history and velocity publish together.
    /// Every selected kinetic point must use this exact duration; a numerical
    /// Newton iteration cannot consume another material interval.
    ///
    /// # Errors
    /// Invalid or mismatched duration, invalid controls, native solve/material
    /// refusal, numerical kinetic-loss or work budget excess, or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        settings: MixedUpSettings,
        controls: FinitePlasticMotionControls,
        cx: &Cx<'_>,
    ) -> Result<MixedPlasticMotionStep<'a, L>, MixedPlasticTetError> {
        checkpoint(cx)?;
        if self
            .solid
            .duration_s()
            .is_some_and(|duration| duration != dt_s)
        {
            return Err(invalid(
                "mixed material kinetics and inertia must use the same physical duration",
            )
            .into());
        }
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
            "total mixed inelastic inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(MixedPlasticMotionStep {
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
