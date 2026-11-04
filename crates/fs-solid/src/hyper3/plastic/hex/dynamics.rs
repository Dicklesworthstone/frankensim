//! Retained finite Hex8 motion with original consistent mass and backward Euler.
//!
//! The parent solid supplies the complete variational finite-element force,
//! tangent and eight independent constitutive histories per cell. The shared
//! inertia owner supplies the momentum residual and exact discrete kinetic-work
//! identity. Every Newton trial uses the same accepted history and velocity;
//! all proposed material, thermal and motion state is returned atomically.
//!
//! Mass integrates the original trilinear shapes with a separate 3×3×3 rule,
//! exact on warped reference bricks. Neither current volume nor assumed mean
//! dilatation changes that inventory. Numerical kinetic loss remains separate
//! from plastic/creep heat and the signed material/projection endpoint defects.
//!
//! Density is explicitly supplied reference data; this binding does not create
//! a calibration or material-card claim. Temperature is prescribed per cell.
//! Conductive heat transport, follower loads, contact and mass transfer are
//! outside this motion owner; the parent formulation's stability limits remain.

use std::sync::Arc;

use fs_exec::Cx;
use fs_feec::highorder::quad1d::gauss_legendre;
use fs_material::plastic::finite::thermal::{AffineThermalFiniteJ2, FiniteThermalJ2Material};

use super::{
    FiniteThermoplasticHexEnergy, FiniteThermoplasticHexError, FiniteThermoplasticHexSolid,
    checkpoint, energy_gate, invalid,
};
use crate::inertia::{BackwardEulerInertia, InertiaError, ReferenceMass};
use crate::linear3::hex::hex8_point;
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

/// The same independent material and kinetic-loss gates as the native tetrahedra.
pub use super::super::thermal::FiniteThermoplasticControls;
use super::super::thermal::FiniteThermoplasticElementEnergy;
/// Shared physical kinetic-work account, including separately reported numerical loss.
pub use crate::inertia::InertialEnergy;

/// Impulses on the same full/mean-dilatation thermal-plastic solid and mass.
pub mod impulse;
pub use crate::inertia::{
    HexImpulsePoint, SolidImpulseControls, SolidImpulseReport, SolidMomentum, SolidPointKinematics,
    SolidPointMobility, SolidPointMobilityMatrix,
};
pub use impulse::FiniteThermoplasticHexMotionImpulse;

/// Complete accepted finite solid, velocity and immutable original mass.
#[derive(Debug)]
pub struct FiniteThermoplasticHexMotion<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    solid: FiniteThermoplasticHexSolid<'a, L>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_corner_density_kg_m3: Option<Arc<[[f64; 8]]>>,
    mass: Arc<ReferenceMass>,
}

enum ReferenceHexDensity<'a> {
    Cell(&'a [f64]),
    Corners(&'a [[f64; 8]]),
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
    /// one finite vector per original node. A separate 27-point reference rule
    /// supplies consistent mass, exact for warped trilinear bricks as well as
    /// affine ones. The eight constitutive histories remain unchanged.
    /// The retained tangent-entry budget includes both material and mass entries
    /// before sparse duplicate assembly: at most 768 entries per original cell.
    pub fn try_new(
        solid: FiniteThermoplasticHexSolid<'a, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticHexError> {
        Self::with_density(
            solid,
            ReferenceHexDensity::Cell(reference_density_kg_m3),
            velocity_mps,
            cx,
        )
    }

    /// Bind a positive trilinear reference-density field, with eight corner
    /// values [kg/m³] per cell in original connectivity order. Cell profiles
    /// may be discontinuous. Three Gauss points per direction integrate
    /// `rho_0 N_a N_b det(J_0)` exactly, up to rounding, even on warped Hex8s.
    /// The eight original constitutive points and their histories are unchanged.
    ///
    /// This profile shares the exact mass rule with the scalar constructor.
    /// An elementwise constant profile preserves the scalar mass bits on both
    /// affine and warped bricks.
    /// Density, reference mass and its moments persist through thermal/plastic
    /// deformation and impulses. No current-temperature density is resampled.
    ///
    /// # Errors
    /// Invalid density, original geometry at a mass point, motion or assembly
    /// budget, unrepresentable mass, or cancellation before publication.
    pub fn try_new_with_reference_density(
        solid: FiniteThermoplasticHexSolid<'a, L>,
        corner_density_kg_m3: &[[f64; 8]],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticHexError> {
        Self::with_density(
            solid,
            ReferenceHexDensity::Corners(corner_density_kg_m3),
            velocity_mps,
            cx,
        )
    }

    fn with_density(
        solid: FiniteThermoplasticHexSolid<'a, L>,
        density: ReferenceHexDensity<'_>,
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
        let ndofs = solid.displacement_m.as_flattened().len();
        // det(J) is quadratic in each coordinate on a warped Q1 map. Constant
        // density needs degree-four mass integration; a Q1 density needs five.
        // One shared 27-point rule supplies both without additional histories.
        let (coordinates, weights) = gauss_legendre(3);
        let reference_points = |index, nodes: [usize; 8]| {
            let positions = nodes.map(|node| solid.nodes_m[node]);
            let mut points = [([0.0; 8], 0.0); 27];
            for (q, entry) in points.iter_mut().enumerate() {
                cx.checkpoint().map_err(|_| InertiaError::Cancelled)?;
                let [i, j, k] = [q / 9, (q / 3) % 3, q % 3];
                let point = hex8_point(
                    positions,
                    [coordinates[i], coordinates[j], coordinates[k]],
                    solid.budget.minimum_scaled_jacobian,
                    index,
                    q,
                )
                .map_err(|error| InertiaError::InvalidInput {
                    what: error.to_string(),
                })?;
                *entry = (
                    point.shapes,
                    point.determinant_m3 * weights[i] * weights[j] * weights[k],
                );
            }
            Ok(points)
        };
        let (mass, reference_density_kg_m3, reference_corner_density_kg_m3) = match density {
            ReferenceHexDensity::Cell(densities) => (
                ReferenceMass::from_hex8_quadrature(
                    solid.elements.iter().enumerate().map(|(index, element)| {
                        Ok((element.nodes, reference_points(index, element.nodes)?))
                    }),
                    densities,
                    ndofs,
                    cx,
                )
                .map_err(inertia_error)?,
                Arc::from(densities),
                None,
            ),
            ReferenceHexDensity::Corners(densities) => {
                if densities.len() != solid.elements.len() {
                    return Err(invalid(
                        "one eight-corner reference density profile per hex required",
                    )
                    .into());
                }
                for rho in densities {
                    checkpoint(cx)?;
                    if rho.iter().any(|value| !value.is_finite() || *value <= 0.0) {
                        return Err(invalid(
                            "positive finite Hex8 corner reference densities required",
                        )
                        .into());
                    }
                }
                let mut means = Vec::with_capacity(densities.len());
                let geometry = solid.elements.iter().zip(densities).enumerate().map(
                    |(index, (element, rho))| {
                        let uniform = rho.iter().all(|value| *value == rho[0]);
                        let mut points = [([0.0; 8], 0.0, 0.0); 27];
                        let (mut volume, mut inventory) = (0.0, 0.0);
                        for (entry, (shapes, weight)) in points
                            .iter_mut()
                            .zip(reference_points(index, element.nodes)?)
                        {
                            let density = if uniform {
                                rho[0]
                            } else {
                                (0..8).map(|a| shapes[a] * rho[a]).sum()
                            };
                            *entry = (shapes, weight, density);
                            volume += weight;
                            inventory += weight * density;
                        }
                        let mean = if uniform { rho[0] } else { inventory / volume };
                        if !mean.is_finite() || mean <= 0.0 {
                            return Err(InertiaError::InvalidInput {
                                what: "unrepresentable Hex8 mean reference density".into(),
                            });
                        }
                        means.push(mean);
                        Ok((element.nodes, points))
                    },
                );
                let mass = ReferenceMass::from_hex8_density_quadrature(geometry, ndofs, cx)
                    .map_err(inertia_error)?;
                (mass, Arc::from(means), Some(Arc::from(densities)))
            }
        };
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3,
            reference_corner_density_kg_m3,
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

    /// Original volume-weighted mean reference density per cell. These means
    /// report inventory; they never replace the graded consistent mass.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Original corner profiles, if explicitly supplied, retained unchanged.
    pub fn reference_corner_density_kg_m3(&self) -> Option<&[[f64; 8]]> {
        self.reference_corner_density_kg_m3.as_deref()
    }

    /// Original integrated mass; thermal/plastic volume changes cannot alter it.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Project a trilinear reference acceleration through the retained mass, M*a.
    ///
    /// One acceleration vector [m/s²] is supplied per original node. The full
    /// nodal forces [N] retain supported coordinates and the original density's
    /// force and moment; a uniform field therefore gives uniform free-body
    /// acceleration even on a graded or warped brick. Cell-mean density and
    /// nodal mass lumping cannot replace this consistent projection.
    ///
    /// Supply the result as held dead loads to [`Self::advance`] exactly once.
    /// Density is not resampled at the current temperature or geometry. A
    /// changing acceleration field, its potential and its tangent remain the
    /// caller's responsibility.
    ///
    /// # Errors
    /// Refuses a wrong node count, nonfinite acceleration or force, and
    /// cancellation before returning the complete load vector.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteThermoplasticHexError> {
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
                reference_corner_density_kg_m3: self.reference_corner_density_kg_m3.clone(),
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
