//! Momentum and progressive finite ductility on the original reference mass.
//!
//! Backward Euler supplies the original consistent momentum force, derivative,
//! endpoint velocity and kinetic-work identity. The parent retains the SAME
//! nonpotential Newton solve, complete finite plastic/damage histories, fixed
//! reference bands and separate physical plastic/damage dissipation accounts.
//! One declared physical interval drives both inertia and ductile accumulation.
//!
//! Material endpoint-work defects and backward-Euler kinetic loss have separate
//! absolute budgets. Numerical loss never becomes material heat or damage
//! release. Every input is borrowed; no candidate displacement, velocity or
//! constitutive history escapes a failed material, solve, work or cancellation
//! gate. This path adds no contact, impulse, mass transfer or adaptive stepping.

use super::*;
use crate::inertia::P1ReferenceMass;

/// Original kinetic change, endpoint momentum work and numerical kinetic loss.
pub use crate::inertia::InertialEnergy;

/// Independent material and inertial accuracy budgets for one physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctileMotionControls {
    /// Parent material/yield/work gates, including summed absolute cell defects.
    pub material: FinitePlasticStepControls,
    /// Maximum nonnegative backward-Euler kinetic loss in this interval [J].
    pub maximum_numerical_kinetic_loss_j: f64,
}

/// Complete accepted finite solid, physical velocity and original mass inventory.
#[derive(Debug)]
pub struct FiniteDuctileTetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    solid: FiniteDuctileTetSolid<'g, 'm, L>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<P1ReferenceMass>,
}

/// One complete accepted interval with unchanged material and kinetic accounts.
#[derive(Debug)]
pub struct FiniteDuctileMotionStep<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub motion: FiniteDuctileTetMotion<'g, 'm, L>,
    /// Native force evidence, including consistent inertia on moving supports.
    pub equilibrium: IncrementRecord,
    pub material_energy: FiniteDuctileTetEnergy,
    /// Unchanged original per-element reference-volume material work receipts.
    pub element_energy: Vec<FiniteDuctileEnergy>,
    pub inertia: InertialEnergy,
    /// External/support work minus material, kinetic and numerical work terms.
    pub balance_error_j: f64,
}

impl<'g, 'm, L> FiniteDuctileTetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Bind an admitted finite specimen to one physical velocity per node and
    /// positive reference density per tetrahedron. The reference mass is
    /// assembled once and retained unchanged as the body deforms. Supplying
    /// density declares original inventory; it does not infer source calibration.
    pub fn try_new(
        solid: FiniteDuctileTetSolid<'g, 'm, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteDuctileTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.state.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity vector per original finite ductile node required",
            )
            .into());
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
        .map_err(SmallStrainTetError::from)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(reference_density_kg_m3),
            reference_nodal_density_kg_m3: None,
            mass: Arc::new(mass),
        })
    }

    /// Bind a positive P1 nodal density on the ORIGINAL material tetrahedra.
    ///
    /// The shared reference-mass owner integrates rho_0*N_a*N_b exactly.
    /// Neither cell means nor damage-scaled stiffness replace that mass.
    /// Its immutable nodal inventory and operator survive every accepted
    /// deformation and damage step. With no matter transfer, current density
    /// is rho_0/J; it is not another multiplier on the reference mass.
    /// Density is an authored inventory, not a material-calibration claim.
    ///
    /// # Errors
    /// Invalid nodal density/velocity, unrepresentable mass, or cancellation.
    pub fn try_new_with_nodal_density(
        solid: FiniteDuctileTetSolid<'g, 'm, L>,
        reference_nodal_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteDuctileTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.state.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite initial velocity vector per original finite ductile node required",
            )
            .into());
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
        .map_err(SmallStrainTetError::from)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(means),
            reference_nodal_density_kg_m3: Some(Arc::from(reference_nodal_density_kg_m3)),
            mass: Arc::new(mass),
        })
    }

    pub fn solid(&self) -> &FiniteDuctileTetSolid<'g, 'm, L> {
        &self.solid
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    /// Original cell-average density, in reference tetrahedron order. Scalar
    /// inputs retain their exact values. Nodal inputs also retain their complete
    /// field through [`Self::reference_nodal_density_kg_m3`]; these reporting means
    /// never replace the density field in momentum or acceleration loads.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    /// Complete original nodal density, or `None` for cellwise density inputs.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_nodal_density_kg_m3.as_deref()
    }
    /// Original total mass, independent of current deformation and damage.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Project a prescribed P1 world acceleration [m/s²] through the SAME
    /// original consistent mass used by the momentum solve.
    ///
    /// The returned nonzero DOF/force pairs are held nodal loads for `advance`.
    /// Uniform gravity retains the graded mass distribution and its moment,
    /// including supported nodes. Plastic flow, damage and current volume do
    /// not rescale this closed material inventory. Do not also supply the same
    /// acceleration as `body_force_n_m3`. Position-dependent field evolution,
    /// its potential and its tangent remain with the external field owner.
    ///
    /// # Errors
    /// Wrong shape, nonfinite acceleration/force, or cancellation.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteDuctileTetError> {
        checkpoint(cx)?;
        if acceleration_m_s2.len() != self.velocity_mps.len()
            || acceleration_m_s2
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "one finite acceleration vector per original finite ductile node required",
            )
            .into());
        }
        let forces = self
            .mass
            .action(acceleration_m_s2.as_flattened(), cx)
            .map_err(SmallStrainTetError::from)?;
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

    /// Advance one original physical interval with absolute displacement targets
    /// and held nodal/body loads. Body force is per ORIGINAL reference volume;
    /// use [`Self::reference_acceleration_load_n`] for gravity on graded inventory.
    /// Free bodies require no artificial supports. The same
    /// `conditions.duration_s` sets momentum and initiation rate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FiniteDuctileMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileMotionStep<'g, 'm, L>, FiniteDuctileTetError> {
        checkpoint(cx)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(
                invalid("finite nonnegative kinetic numerical-loss budget required").into(),
            );
        }
        let inertia = BackwardEulerInertia::new(
            &self.mass,
            self.solid.state.displacement_m.as_flattened(),
            self.velocity_mps.as_flattened(),
            conditions.duration_s,
            cx,
        )
        .map_err(SmallStrainTetError::from)?;
        let (step, motion) = self.solid.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            conditions,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let (velocity_mps, inertia) =
            motion.expect("native inertial solve returns complete motion");
        gate(
            "backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        gate(
            "native kinetic work identity",
            inertia.balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        let balance_error_j = step.energy.balance_error_j + inertia.balance_error_j;
        gate(
            "complete finite ductile inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteDuctileMotionStep {
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
            element_energy: step.element_energy,
            inertia,
            balance_error_j,
        })
    }
}
