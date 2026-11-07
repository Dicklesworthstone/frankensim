//! In-plane momentum evolution with native finite-strain viscoelastic memory.
//!
//! The same physical interval drives the native material update and consistent
//! backward-Euler inertia in the original plane-stress Newton solve. Original
//! density times original area and thickness defines a closed mass inventory.
//! The transverse column remains locally quasistatic; no transverse kinetic
//! energy, shell bending or rotational inertia is inferred.
//!
//! Physical branch dissipation, signed material endpoint error and nonnegative
//! backward-Euler kinetic loss are separate accounts with explicit gates. This
//! is a first-order motion discretization; refining time controls its numerical
//! loss and phase error. It does not supply a zero-algorithmic-damping or exact
//! modal-frequency integrator.

use std::sync::Arc;

use crate::inertia::{BackwardEulerInertia, ReferenceMass};
use crate::{NewtonReport, NewtonSettings, Patch};
use fs_exec::Cx;

use super::{
    FiniteViscoPlaneEnergy, FiniteViscoPlaneError, FiniteViscoPlaneReferenceEnergy,
    FiniteViscoPlaneSolid, FiniteViscoStepControls, ViscoelasticLawField, checkpoint, gate,
    inertia_error, invalid,
};

pub use crate::hyper3::visco::dynamics::FiniteViscoMotionControls;
pub use crate::inertia::InertialEnergy;

/// Accepted material state, velocity and immutable original consistent mass.
#[derive(Debug)]
pub struct FiniteViscoPlaneMotion<'a> {
    solid: FiniteViscoPlaneSolid<'a>,
    velocity_mps: Arc<[[f64; 2]]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<ReferenceMass>,
}

/// Complete proposed interval. Dropping it preserves the previous motion.
#[derive(Debug)]
pub struct FiniteViscoPlaneMotionStep<'a> {
    pub motion: FiniteViscoPlaneMotion<'a>,
    pub equilibrium: NewtonReport,
    pub material_energy: FiniteViscoPlaneEnergy,
    pub inertia: InertialEnergy,
    /// Support/dead-load equilibrium work plus the native kinetic identity
    /// residual. Material endpoint error and numerical kinetic loss have their
    /// own explicit accounts and are not hidden in this balance.
    pub balance_error_j: f64,
    pub reactions_n: Vec<(usize, f64)>,
}

/// Reference change at fixed in-plane displacement, velocity and physical time.
/// Reference mass is shared unchanged; current thickness/temperature never
/// creates or removes mass.
#[derive(Debug)]
pub struct FiniteViscoPlaneMotionReferenceChange<'a> {
    pub motion: FiniteViscoPlaneMotion<'a>,
    pub energy: FiniteViscoPlaneReferenceEnergy,
}

impl<'a> FiniteViscoPlaneMotion<'a> {
    /// Bind an already accepted solid to explicit original per-element density
    /// and one velocity vector per original node. Existing material memory,
    /// displacement and age survive. The mass assembly has an explicit byte cap.
    pub fn try_new(
        solid: FiniteViscoPlaneSolid<'a>,
        reference_density_kg_m3: &[f64],
        initial_velocity_mps: &[[f64; 2]],
        maximum_mass_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteViscoPlaneError> {
        checkpoint(cx)?;
        if initial_velocity_mps.len() != solid.reference_mesh().nodes.len()
            || initial_velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite velocity vector per original plane node required",
            ));
        }
        let mass = ReferenceMass::from_plane(
            solid.reference_mesh(),
            solid.reference_thickness_m(),
            reference_density_kg_m3,
            maximum_mass_assembly_bytes,
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: initial_velocity_mps.into(),
            reference_density_kg_m3: reference_density_kg_m3.into(),
            reference_nodal_density_kg_m3: None,
            mass: Arc::new(mass),
        })
    }

    /// Bind a positive original P1/Q1 nodal density [kg/m³] to native momentum.
    /// The consistent mass integrates `rho_0*N_a*N_b` exactly, including a
    /// bilinearly mapped Q1 element. Cell means are diagnostics, not replacement
    /// densities. The complete field and prepared mass stay fixed through all
    /// subsequent deformation, relaxation and material-reference changes.
    /// This authored reference inventory does not imply a material calibration.
    ///
    /// # Errors
    /// Invalid density, velocity, original geometry, mass-assembly byte cap,
    /// unrepresentable mass, or cancellation.
    pub fn try_new_with_nodal_density(
        solid: FiniteViscoPlaneSolid<'a>,
        reference_nodal_density_kg_m3: &[f64],
        initial_velocity_mps: &[[f64; 2]],
        maximum_mass_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteViscoPlaneError> {
        checkpoint(cx)?;
        if initial_velocity_mps.len() != solid.reference_mesh().nodes.len()
            || initial_velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite velocity vector per original plane node required",
            ));
        }
        let (mass, means) = ReferenceMass::from_plane_nodal_density(
            solid.reference_mesh(),
            solid.reference_thickness_m(),
            reference_nodal_density_kg_m3,
            maximum_mass_assembly_bytes,
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: initial_velocity_mps.into(),
            reference_density_kg_m3: means.into(),
            reference_nodal_density_kg_m3: Some(reference_nodal_density_kg_m3.into()),
            mass: Arc::new(mass),
        })
    }

    pub fn solid(&self) -> &FiniteViscoPlaneSolid<'a> {
        &self.solid
    }
    pub fn velocity_mps(&self) -> &[[f64; 2]] {
        &self.velocity_mps
    }
    /// Reference cell means in original element order. Graded bindings use
    /// original-volume weighting and retain their full nodal field separately.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    /// Complete original nodal field, or None for a cellwise-density binding.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_nodal_density_kg_m3.as_deref()
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Preserve closed inventory and actual velocity through a prescribed
    /// material-reference event. Transverse strain is re-equilibrated by the
    /// original solid owner, including exact pressure-volume reference input.
    pub fn change_reference<'b>(
        &self,
        target_material: ViscoelasticLawField<'b>,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoPlaneMotionReferenceChange<'b>, FiniteViscoPlaneError>
    where
        'a: 'b,
    {
        let changed = self.solid.change_reference(target_material, controls, cx)?;
        checkpoint(cx)?;
        Ok(FiniteViscoPlaneMotionReferenceChange {
            motion: FiniteViscoPlaneMotion {
                solid: changed.solid,
                velocity_mps: Arc::clone(&self.velocity_mps),
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                reference_nodal_density_kg_m3: self.reference_nodal_density_kg_m3.clone(),
                mass: Arc::clone(&self.mass),
            },
            energy: changed.energy,
        })
    }

    /// Project an explicitly supplied reference nodal acceleration through the
    /// SAME consistent mass used by momentum. Returns physical nodal forces,
    /// including supported rows; no current-density resampling or lumping occurs.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 2]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteViscoPlaneError> {
        checkpoint(cx)?;
        if acceleration_m_s2.len() != self.velocity_mps.len()
            || acceleration_m_s2
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite acceleration vector per original plane node required",
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

    /// Advance one positive physical interval with absolute endpoint supports
    /// [m], held nominal boundary tractions [Pa], and explicit nodal forces [N].
    /// Duplicate nodal entries accumulate in original order. Consistent mass
    /// enters force, tangent, potential, normalization and support reactions.
    ///
    /// Every material trial uses the same accepted memory and physical interval.
    /// All material, equilibrium and kinetic gates pass before returning the
    /// successor. An unsupported free body can translate without artificial pins.
    /// Reduce the physical interval when numerical-loss or work gates refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        dirichlet: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        traction: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        nodal_forces_n: &[(usize, f64)],
        settings: NewtonSettings,
        controls: FiniteViscoMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteViscoPlaneMotionStep<'a>, FiniteViscoPlaneError> {
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
            self.solid.displacement_m().as_flattened(),
            self.velocity_mps.as_flattened(),
            dt_s,
            cx,
        )
        .map_err(inertia_error)?;
        let (step, motion) = self.solid.advance_loaded(
            dt_s,
            dirichlet,
            traction,
            nodal_forces_n,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let (velocity, inertia) =
            motion.ok_or_else(|| invalid("missing admitted inertial endpoint"))?;
        if velocity.len() != self.velocity_mps.len() * 2 {
            return Err(invalid(
                "inertial velocity does not cover original plane DOFs",
            ));
        }
        gate(
            "backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        let balance_error_j = step.energy.equilibrium_work_residual_j + inertia.balance_error_j;
        gate(
            "total thin viscoelastic inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteViscoPlaneMotionStep {
            motion: Self {
                solid: step.solid,
                velocity_mps: velocity.as_chunks::<2>().0.to_vec().into(),
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                reference_nodal_density_kg_m3: self.reference_nodal_density_kg_m3.clone(),
                mass: Arc::clone(&self.mass),
            },
            equilibrium: step.equilibrium,
            material_energy: step.energy,
            inertia,
            balance_error_j,
            reactions_n: step.reactions_n,
        })
    }
}
