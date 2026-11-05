//! Consistent quadratic inertia and retained mixed viscoelastic motion.
//!
//! The native backward-Euler momentum equation enters the SAME P2/P1 solve as
//! the frozen-history material. Pressure has no mass. Original reference
//! density and the signed consistent P2 mass are immutable; no row-sum lumping
//! or deformation-dependent mass update occurs. Material loss, signed material
//! endpoint-work defect and numerical kinetic loss remain separate.
//! Finite-compliance motion also admits instantaneous nodal impulses through
//! this same consistent mass, retaining pressure and every material history.
//! Exact incompressibility has a separate pressure-constrained impulse entry
//! using the original weak-volume derivative at the accepted geometry.
//!
//! Exact pressure constrains endpoint volume weakly; velocity is the native
//! backward-Euler displacement chord, not a separate tangent projection. A
//! closed work identity does not establish time accuracy. Numerical-loss and
//! work allowances can reject the complete candidate without changing history.
//! This first-order path permits budgeted numerical damping; it is not the
//! instrument-render doctrine's zero-algorithmic-damping/discrete-gradient integrator.
//!
//! The separate finite-compliance force/contact/drift entry uses accepted Piola
//! forces, the same consistent mass and a pressure-only native material endpoint.
//! It retains post-contact velocity during drift and reports numerical splitting
//! work independently of relaxation. Complete path admission belongs to the caller.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::VolumetricTreatment;

use super::super::{MixedUpSettings, MixedUpSolution, inertia_error};
use super::{
    MixedViscoLawField, MixedViscoReferenceChangeControls, MixedViscoTetEnergy, MixedViscoTetError,
    MixedViscoTetSolid, checkpoint, gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, ReferenceMass};
use crate::linear3::tet10::TET10_EDGES;

pub use crate::hyper3::mixed::impulse::{MixedImpulseControls, MixedImpulseReport};
pub use crate::hyper3::visco::dynamics::FiniteViscoMotionControls;
pub use crate::hyper3::visco::dynamics::{
    FiniteViscoContactDrift, FiniteViscoContactDriftControls,
};
pub use crate::inertia::{
    InertialEnergy, SolidImpulseControls, SolidImpulseReport, SolidMomentum, SolidMotionMomentum,
    SolidPointKinematics, SolidPointMobility, SolidPointMobilityMatrix, Tet10ImpulsePoint,
};

mod impulse;
pub use impulse::{
    MixedViscoConstrainedImpulse, MixedViscoConstrainedResponse, MixedViscoMotionImpulse,
};

mod split;
pub use split::MixedViscoForceKick;

/// Immutable accepted displacement, continuous pressure, memory and velocity.
#[derive(Debug)]
pub struct MixedViscoTetMotion<'a> {
    solid: MixedViscoTetSolid<'a>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_corner_density_kg_m3: Option<Arc<[[f64; 4]]>>,
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

/// Admitted reference/pressure change with unchanged original-mass motion.
#[derive(Debug)]
pub struct MixedViscoMotionReferenceChange<'a> {
    pub motion: MixedViscoTetMotion<'a>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub remainder_reference_work_j: f64,
    pub bulk_reference_work_j: f64,
    /// Signed material-reference exchange; pressure-constraint work is reported
    /// separately below. Velocity, kinetic energy and viscous loss are unchanged.
    pub reference_work_j: f64,
    /// Signed weak-constraint contribution retained by the native body account:
    /// `stored_after - stored_before = reference_work - constraint_defect`.
    pub pressure_constraint_defect_j: f64,
    pub relative_volume_residual: f64,
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
            reference_corner_density_kg_m3: None,
            mass: Arc::new(mass),
        })
    }

    /// Bind a positive P1 reference-density field to the original P2 inertia.
    ///
    /// Four densities [kg/m³] per element follow its corner connectivity. P1
    /// interpolation remains positive between corners; unconstrained quadratic
    /// density interpolation would not. Material, pressure and viscous memory
    /// are retained, and density is never re-queried on the deformed geometry.
    /// Four Duffy points per direction integrate affine graded P2 mass exactly;
    /// a varying profile on a curved map requires at least six. The cached
    /// reference quadrature supplies both inertia and acceleration loads.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite or incomplete density, insufficient curved
    /// quadrature, unrepresentable mass, assembly-budget excess or cancellation.
    pub fn try_new_with_reference_density(
        solid: MixedViscoTetSolid<'a>,
        corner_density_kg_m3: &[[f64; 4]],
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
        if corner_density_kg_m3.len() != solid.mesh.elements.len()
            || corner_density_kg_m3
                .iter()
                .flatten()
                .any(|rho| !rho.is_finite() || *rho <= 0.0)
        {
            return Err(
                invalid("four positive reference densities per original Tet10 required").into(),
            );
        }
        for (element, density) in solid.mesh.elements.iter().zip(corner_density_kg_m3) {
            checkpoint(cx)?;
            if density.iter().all(|rho| *rho == density[0]) {
                continue;
            }
            let nodes = &element.connectivity;
            let affine = TET10_EDGES.iter().enumerate().all(|(edge, &[a, b])| {
                (0..3).all(|c| {
                    solid.nodes_m[nodes[edge + 4]][c]
                        == f64::midpoint(solid.nodes_m[nodes[a]][c], solid.nodes_m[nodes[b]][c])
                })
            });
            if !affine && solid.quadrature_points_per_dir < 6 {
                return Err(invalid(
                    "curved Tet10 graded reference density requires at least six Duffy points",
                )
                .into());
            }
        }
        let mass = ReferenceMass::from_tet10_density_quadrature(
            solid
                .mesh
                .elements
                .iter()
                .zip(corner_density_kg_m3)
                .map(|(element, &density)| {
                    (
                        element.connectivity,
                        element.points.iter().map(move |q| {
                            (
                                q.shapes,
                                q.weight_m3,
                                corner_density_at(density, q.pressure_shapes),
                            )
                        }),
                    )
                }),
            solid.displacement_m.as_flattened().len(),
            maximum_mass_assembly_bytes,
            cx,
        )
        .map_err(inertia_error)?;
        let mut means = Vec::with_capacity(corner_density_kg_m3.len());
        for (element, &density) in solid.mesh.elements.iter().zip(corner_density_kg_m3) {
            checkpoint(cx)?;
            let mean = if density.iter().all(|rho| *rho == density[0]) {
                density[0]
            } else {
                element
                    .points
                    .iter()
                    .map(|q| {
                        corner_density_at(density, q.pressure_shapes)
                            * (q.weight_m3 / element.reference_volume_m3)
                    })
                    .sum::<f64>()
            };
            if !mean.is_finite() || mean <= 0.0 {
                return Err(invalid("reference-volume mean density is unrepresentable").into());
            }
            means.push(mean);
        }
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: means.into(),
            reference_corner_density_kg_m3: Some(Arc::from(corner_density_kg_m3)),
            mass: Arc::new(mass),
        })
    }

    pub fn solid(&self) -> &MixedViscoTetSolid<'a> {
        &self.solid
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    /// Original reference-volume mean densities. Reporting only: graded
    /// inertia and acceleration loads never substitute these for the field.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    /// Original positive corner profiles, retained by all motion transitions.
    pub fn reference_corner_density_kg_m3(&self) -> Option<&[[f64; 4]]> {
        self.reference_corner_density_kg_m3.as_deref()
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Admit a prescribed-reference change at the held total deformation and
    /// caller-supplied pressure, without changing physical motion or inventory.
    ///
    /// The native solid re-admits its weak volume constraint and accounts for
    /// remainder and reconstructed bulk energy separately. In particular, a
    /// pressure choice cannot repair an incompatible exact-volume change at
    /// fixed displacement. The original corner density and signed consistent
    /// P2 mass are shared directly; no mass assembly or averaging occurs.
    ///
    /// Velocity and time are retained, and the new accepted body identity
    /// invalidates old prepared contact/force/drift responses. This successor
    /// does not borrow the predecessor wrapper. Any refusal or cancellation
    /// leaves the previous motion, pressure and histories unchanged.
    pub fn change_reference<'b>(
        &self,
        target_material: MixedViscoLawField<'b>,
        pressure_pa: &[f64],
        controls: MixedViscoReferenceChangeControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoMotionReferenceChange<'b>, MixedViscoTetError>
    where
        'a: 'b,
    {
        let change = self
            .solid
            .change_reference(target_material, pressure_pa, controls, cx)?;
        let mut velocity_mps = Vec::with_capacity(self.velocity_mps.len());
        for chunk in self.velocity_mps.chunks(256) {
            checkpoint(cx)?;
            velocity_mps.extend_from_slice(chunk);
        }
        checkpoint(cx)?;
        Ok(MixedViscoMotionReferenceChange {
            motion: MixedViscoTetMotion {
                solid: change.solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                reference_corner_density_kg_m3: self
                    .reference_corner_density_kg_m3
                    .as_ref()
                    .map(Arc::clone),
                mass: Arc::clone(&self.mass),
            },
            stored_before_j: change.stored_before_j,
            stored_after_j: change.stored_after_j,
            remainder_reference_work_j: change.remainder_reference_work_j,
            bulk_reference_work_j: change.bulk_reference_work_j,
            reference_work_j: change.reference_work_j,
            pressure_constraint_defect_j: change.pressure_constraint_defect_j,
            relative_volume_residual: change.relative_volume_residual,
        })
    }
    /// Project a P2 reference acceleration [m/s²] through the retained mass.
    ///
    /// Supply the returned original nodal forces [N] as dead loads exactly
    /// once. Signed corner and edge contributions, supported coordinates and
    /// density-weighted moments are preserved. A constant acceleration gives
    /// constant free-body acceleration for graded and curved elements alike.
    /// This does not introduce a current-position-dependent gravity potential.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, MixedViscoTetError> {
        checkpoint(cx)?;
        if acceleration_m_s2.len() != self.solid.nodes_m.len()
            || acceleration_m_s2.iter().flatten().any(|a| !a.is_finite())
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
    /// Original pressure treatment retained by every accepted native interval.
    pub fn volumetric_treatment(&self) -> VolumetricTreatment {
        self.solid.volumetric
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
        self.advance_accounted(
            dt_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            pressure_gauge,
            settings,
            controls,
            cx,
        )
        .map(|(step, _)| step)
    }

    /// Advance the same native interval and report actual original-mass momentum.
    /// All signed quadratic load weights and original support reactions enter
    /// the receipt at the accepted deformed geometry. Endpoint applied/support
    /// torque and backward-Euler angular drift remain separate; consumers must
    /// bound physical residuals AND numerical drift. The original load assembly
    /// is reused without another mass assembly or material update.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_momentum(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        settings: MixedUpSettings,
        controls: FiniteViscoMotionControls,
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<(MixedViscoMotionStep<'a>, SolidMotionMomentum), MixedViscoTetError> {
        if origin_m.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite momentum origin required").into());
        }
        let (step, force) = self.advance_accounted(
            dt_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            pressure_gauge,
            settings,
            controls,
            cx,
        )?;
        let momentum = self
            .mass
            .backward_euler_momentum(
                self.solid.nodes_m,
                &self.solid.displacement_m,
                &step.motion.solid.displacement_m,
                &self.velocity_mps,
                &step.motion.velocity_mps,
                &force,
                &step.equilibrium.reactions_n,
                dt_s,
                origin_m,
                cx,
            )
            .map_err(inertia_error)?;
        checkpoint(cx)?;
        Ok((step, momentum))
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_accounted(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        settings: MixedUpSettings,
        controls: FiniteViscoMotionControls,
        cx: &Cx<'_>,
    ) -> Result<(MixedViscoMotionStep<'a>, Vec<f64>), MixedViscoTetError> {
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
        let (step, motion, force) = self.solid.advance_loaded(
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
        Ok((
            MixedViscoMotionStep {
                motion: Self {
                    solid: step.solid,
                    velocity_mps,
                    reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                    reference_corner_density_kg_m3: self.reference_corner_density_kg_m3.clone(),
                    mass: Arc::clone(&self.mass),
                },
                equilibrium: step.equilibrium,
                material_energy: step.energy,
                inertia,
                balance_error_j,
            },
            force,
        ))
    }
}

fn corner_density_at(density: [f64; 4], lambda: [f64; 4]) -> f64 {
    if density.iter().all(|rho| *rho == density[0]) {
        density[0]
    } else {
        (0..4).map(|a| lambda[a] * density[a]).sum()
    }
}
