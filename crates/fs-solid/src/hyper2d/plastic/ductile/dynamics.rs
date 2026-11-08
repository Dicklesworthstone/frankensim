//! Backward-Euler motion of progressive finite ductile P1/Q1 sheets.
//!
//! The original consistent reference mass retains each element's physical
//! thickness and volumetric density. Neither damage, transverse contraction
//! nor current volume changes this closed material inventory. One physical
//! interval drives the material return, initiation and momentum equation.
//! Material plastic/damage loss, signed material endpoint-work defect and
//! numerical kinetic loss stay separate; none is silently converted to heat.

use super::*;
use crate::inertia::{BackwardEulerInertia, ReferenceMass};

/// Original kinetic change, endpoint momentum work and numerical kinetic loss.
pub use crate::inertia::InertialEnergy;

/// Independent material and inertial accuracy bounds for one physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctileSheetMotionControls {
    pub material: FinitePlasticStepControls,
    /// Maximum backward-Euler numerical kinetic loss [J], never material heat.
    pub maximum_numerical_kinetic_loss_j: f64,
}

/// Complete accepted sheet, physical velocity and immutable reference mass.
#[derive(Clone, Debug)]
pub struct FiniteDuctileSheetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    sheet: FiniteDuctileSheet<'g, 'm, L>,
    velocity_mps: Vec<[f64; 2]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<ReferenceMass>,
}

/// Complete candidate, published only after equilibrium and every work gate.
#[derive(Debug)]
pub struct FiniteDuctileSheetMotionStep<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub motion: FiniteDuctileSheetMotion<'g, 'm, L>,
    /// Original force evidence, including inertial reactions on moving supports.
    pub equilibrium: IncrementRecord,
    pub material_energy: FiniteDuctileSheetEnergy,
    /// Unchanged material receipts in original element/quadrature order.
    pub point_energy: Vec<FiniteDuctileEnergy>,
    pub inertia: InertialEnergy,
    pub balance_error_j: f64,
}

impl<'g, 'm, L> FiniteDuctileSheetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Bind an admitted sheet to one finite velocity per node and one positive
    /// original volumetric density per element. Mass is assembled once with
    /// the original physical thicknesses and shared native mass quadrature.
    /// The byte cap bounds consistent-mass assembly before allocation.
    pub fn try_new(
        sheet: FiniteDuctileSheet<'g, 'm, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 2]],
        maximum_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteDuctileSheetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != sheet.mesh.nodes.len()
            || velocity_mps.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite velocity pair per original ductile sheet node required",
            )
            .into());
        }
        let mass = ReferenceMass::from_plane_layers(
            sheet.mesh,
            sheet.thickness_m,
            reference_density_kg_m3,
            maximum_assembly_bytes,
            cx,
        )
        .map_err(SmallStrainTetError::from)?;
        checkpoint(cx)?;
        Ok(Self {
            sheet,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(reference_density_kg_m3),
            mass: Arc::new(mass),
        })
    }

    pub fn sheet(&self) -> &FiniteDuctileSheet<'g, 'm, L> {
        &self.sheet
    }
    pub fn velocity_mps(&self) -> &[[f64; 2]] {
        &self.velocity_mps
    }
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Project an original nodal acceleration [m/s²] through the SAME retained
    /// consistent mass. This preserves the distribution of a heterogeneous
    /// inventory, including supported nodes. Returned pairs are held nodal
    /// forces for `advance`; do not also apply that gravity as a body force.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 2]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteDuctileSheetError> {
        checkpoint(cx)?;
        if acceleration_m_s2.len() != self.velocity_mps.len()
            || acceleration_m_s2.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite acceleration pair per original ductile sheet node required",
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

    fn with_body_force(
        &self,
        nodal: &[(usize, f64)],
        body: [f64; 2],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteDuctileSheetError> {
        let mut forces = vec![0.0; self.velocity_mps.len() * 2];
        for (index, &(dof, value)) in nodal.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if dof >= forces.len() || !value.is_finite() {
                return Err(invalid("invalid ductile sheet dead nodal force").into());
            }
            forces[dof] += value;
        }
        for element in self.sheet.elements.iter() {
            checkpoint(cx)?;
            for point in &element.points {
                for a in 0..element.arity {
                    let volume = point.shape[a] * point.volume_m3;
                    for i in 0..2 {
                        let dof = 2 * element.nodes[a] + i;
                        forces[dof] = body[i].mul_add(volume, forces[dof]);
                    }
                }
            }
        }
        let mut loads = Vec::with_capacity(forces.len());
        for (dof, force) in forces.into_iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if !force.is_finite() {
                return Err(invalid("ductile sheet body/nodal load is not representable").into());
            }
            if force != 0.0 {
                loads.push((dof, force));
            }
        }
        Ok(loads)
    }

    /// Advance one physical interval with absolute support targets and held
    /// forces. `body_force_n_m3` is per ORIGINAL reference volume; gravity on
    /// varying density should use `reference_acceleration_load_n` instead.
    ///
    /// The same `conditions.duration_s` drives momentum, plastic kinetics and
    /// initiation. Free translation requires no artificial supports. The old
    /// sheet, velocity and reference inventory remain intact on any refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 2],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FiniteDuctileSheetMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileSheetMotionStep<'g, 'm, L>, FiniteDuctileSheetError> {
        checkpoint(cx)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
            || body_force_n_m3.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite body force and nonnegative kinetic numerical-loss budget required",
            )
            .into());
        }
        let combined;
        let loads = if body_force_n_m3 == [0.0; 2] {
            nodal_forces_n
        } else {
            combined = self.with_body_force(nodal_forces_n, body_force_n_m3, cx)?;
            &combined
        };
        let inertia = BackwardEulerInertia::new(
            &self.mass,
            self.sheet.displacement_m.as_flattened(),
            self.velocity_mps.as_flattened(),
            conditions.duration_s,
            cx,
        )
        .map_err(SmallStrainTetError::from)?;
        let (step, motion) = self.sheet.advance_loaded(
            prescribed_m,
            loads,
            conditions,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let (velocity_mps, inertia) = motion.expect("native inertial sheet solve returns motion");
        energy_gate(
            "backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        energy_gate(
            "native kinetic work identity",
            inertia.balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        let balance_error_j = step.energy.balance_error_j + inertia.balance_error_j;
        energy_gate(
            "complete finite ductile sheet inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteDuctileSheetMotionStep {
            motion: Self {
                sheet: step.sheet,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
            },
            equilibrium: step.equilibrium,
            material_energy: step.energy,
            point_energy: step.point_energy,
            inertia,
            balance_error_j,
        })
    }
}
