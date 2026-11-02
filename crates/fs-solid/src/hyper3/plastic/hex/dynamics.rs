//! Retained finite Hex8 motion with original consistent mass and backward Euler.
//!
//! The parent solid supplies the complete variational finite-element force,
//! tangent and eight independent constitutive histories per cell. The shared
//! inertia owner supplies the momentum residual and exact discrete kinetic-work
//! identity. Every Newton trial uses the same accepted history and velocity;
//! all proposed material, thermal and motion state is returned atomically.
//!
//! Mass integrates the original trilinear shapes over the original quadrature
//! volumes. Neither current volume nor the assumed mean dilatation changes that
//! inventory. Numerical kinetic loss remains separate from plastic/creep heat
//! and from the signed material/projection endpoint integration defects.
//!
//! Density is explicitly supplied reference data; this binding does not create
//! a calibration or material-card claim. Temperature is prescribed per cell.
//! Conductive heat transport, follower loads, contact and mass transfer are
//! outside this motion owner; the parent formulation's stability limits remain.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::plastic::finite::thermal::{AffineThermalFiniteJ2, FiniteThermalJ2Material};

use super::{
    FiniteThermoplasticHexEnergy, FiniteThermoplasticHexError, FiniteThermoplasticHexSolid,
    checkpoint, energy_gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, InertiaError, ReferenceMass};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

/// The same independent material and kinetic-loss gates as the native tetrahedra.
pub use super::super::thermal::FiniteThermoplasticControls;
use super::super::thermal::FiniteThermoplasticElementEnergy;
/// Shared physical kinetic-work account, including separately reported numerical loss.
pub use crate::inertia::InertialEnergy;

/// Complete accepted finite solid, velocity and immutable original mass.
#[derive(Debug)]
pub struct FiniteThermoplasticHexMotion<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    solid: FiniteThermoplasticHexSolid<'a, L>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<ReferenceMass>,
}

/// Fully admitted candidate; the borrowed previous motion remains unchanged.
#[derive(Debug)]
pub struct FiniteThermoplasticHexMotionStep<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2>
{
    pub motion: FiniteThermoplasticHexMotion<'a, L>,
    /// Native physical equilibrium, including inertial support reactions.
    pub equilibrium: IncrementRecord,
    /// Original-volume material accounts in element-major, eight-point order.
    pub quadrature_energy: Vec<FiniteThermoplasticElementEnergy>,
    /// Native storage, physical heat and separate material/projection defects.
    pub material_energy: FiniteThermoplasticHexEnergy,
    pub inertia: InertialEnergy,
    /// Combined external/material/thermal/kinetic endpoint balance residual.
    pub balance_error_j: f64,
}

fn inertia_error(error: InertiaError) -> FiniteThermoplasticHexError {
    match error {
        InertiaError::InvalidInput { what } => invalid(what).into(),
        InertiaError::Cancelled => SmallStrainTetError::Cancelled.into(),
    }
}

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticHexMotion<'a, L> {
    /// Bind an already admitted complete solid to fixed reference mass and
    /// initial velocity. No constitutive history or source identity is rebuilt.
    ///
    /// Density has one positive finite value per original cell and velocity
    /// one finite vector per original node. The original eight-point geometry
    /// supplies the consistent trilinear mass, exact on affine reference bricks.
    /// The retained tangent-entry budget includes both material and mass entries
    /// before sparse duplicate assembly: at most 768 entries per original cell.
    pub fn try_new(
        solid: FiniteThermoplasticHexSolid<'a, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticHexError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "finite hex motion needs one finite initial velocity per original node",
            )
            .into());
        }
        // The material tangent contributes 24² entries per cell. The consistent
        // mass adds at most three 8² blocks, even after shared-node assembly.
        // Admit their combined bound before allocating the additional mass.
        let entries = solid
            .elements
            .len()
            .checked_mul(24 * 24 + 3 * 8 * 8)
            .ok_or_else(|| invalid("finite hex inertial tangent count overflow"))?;
        if entries > solid.budget.maximum_tangent_entries {
            return Err(invalid(format!(
                "finite hex inertial tangent entries {entries} exceed budget {}",
                solid.budget.maximum_tangent_entries
            ))
            .into());
        }
        let mass = ReferenceMass::from_hex8_quadrature(
            solid.elements.iter().map(|element| {
                (
                    element.nodes,
                    element.points.map(|point| (point.shapes, point.volume_m3)),
                )
            }),
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

    /// Accepted displacement, temperature, original families and all histories.
    pub fn solid(&self) -> &FiniteThermoplasticHexSolid<'a, L> {
        &self.solid
    }

    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }

    /// Original reference density in cell order, retained through deformation.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Original integrated mass; thermal/plastic volume changes cannot alter it.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Advance one actual physical interval with absolute supports and held
    /// dead nodal loads. Material kinetics and momentum use this same duration.
    /// Free bodies require no artificial constraints. Inertia couples the
    /// prescribed coordinates through the original consistent mass matrix.
    ///
    /// Every refusal, including a numerical kinetic-loss budget rejection after
    /// equilibrium, leaves all accepted motion and constitutive state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        dt_s: f64,
        settings: IncrementSettings,
        controls: FiniteThermoplasticControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticHexMotionStep<'a, L>, FiniteThermoplasticHexError> {
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
            self.solid.displacement_m.as_flattened(),
            self.velocity_mps.as_flattened(),
            dt_s,
            cx,
        )
        .map_err(inertia_error)?;
        let (step, receipt) = self.solid.advance_with_inertia(
            new_temperatures_k,
            prescribed_m,
            nodal_forces_n,
            dt_s,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let (velocity_mps, inertia) = receipt.expect("inertial hex endpoint retains its motion");
        energy_gate(
            "backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        // The parent removes inertial endpoint work from its material balance.
        // Restore the independent kinetic identity to obtain the full physical
        // balance, without adding its numerical loss to caloric heat.
        let balance_error_j = step.energy.balance_error_j + inertia.balance_error_j;
        energy_gate(
            "total finite hex inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FiniteThermoplasticHexMotionStep {
            motion: Self {
                solid: step.solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
            },
            equilibrium: step.equilibrium,
            quadrature_energy: step.quadrature_energy,
            material_energy: step.energy,
            inertia,
            balance_error_j,
        })
    }
}
