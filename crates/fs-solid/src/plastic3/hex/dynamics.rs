//! Backward-Euler motion of the native Standard or B-bar inelastic Hex8 solid.
//!
//! The original eight constitutive histories per cell remain frozen throughout
//! Newton and backtracking. The same consistent tangent, full support reactions
//! and converged endpoint stress are used before history is committed once.
//!
//! Consistent mass uses the original trilinear shapes and a separate 3x3x3
//! reference rule, exact for warped bricks with constant or trilinear density.
//! B-bar projection changes neither reference inventory nor inertial coupling.
//! Endpoint inertial work equals kinetic change plus the separately reported
//! backward-Euler kinetic loss; that numerical loss is never material heat.
//!
//! This is a first-order motion adapter, not a caloric or constitutive energy
//! owner. A timed SmallStrainLaw must be configured by its caller for this same
//! physical interval: the generic trait does not expose a duration to check.
//! Density is immutable reference inventory, not a current-temperature update.
//! Finite deformation, follower loads, contact and mass transfer are not added.

use super::{Element, SmallStrainHexProblem, SmallStrainHexSolution, checkpoint, invalid};
use crate::inertia::{BackwardEulerInertia, InertiaError, ReferenceMass};
use crate::linear3::hex::hex8_point;
use crate::plastic3::SmallStrainTetError;
use fs_exec::Cx;
use fs_feec::highorder::quad1d::gauss_legendre;
use fs_material::SmallStrainLaw;

/// Original kinetic-work receipt, with numerical loss separate from heat.
pub use crate::inertia::InertialEnergy;

/// Complete motion proposal; all borrowed old displacement/history is retained.
#[derive(Clone, Debug)]
pub struct BackwardEulerHexStep<S> {
    /// Original constitutive endpoint and inertial support reactions.
    pub solution: SmallStrainHexSolution<S>,
    /// Velocity from the actual accepted displacement increment divided by dt.
    pub velocity_mps: Vec<[f64; 3]>,
    /// Consistent-mass kinetic storage, work and numerical kinetic loss.
    pub energy: InertialEnergy,
}

#[derive(Clone, Copy)]
pub(super) enum ReferenceDensity<'a> {
    Cell(&'a [f64]),
    Corners(&'a [[f64; 8]]),
}

impl<L: SmallStrainLaw> SmallStrainHexProblem<'_, L> {
    /// Advance one physical interval from the accepted complete nodal motion
    /// and eight histories per element. Exactly one load factor, 1.0, is
    /// required: forces and prescribed displacements are absolute endpoints.
    ///
    /// Each positive density belongs to an original cell. The caller retains
    /// that reference inventory across motion and thermal intervals; changing
    /// it requires a separate matter/momentum-transfer law. Free bodies need
    /// no artificial restraints. All refusals leave the inputs unchanged.
    pub fn advance_backward_euler(
        &self,
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        quadrature_states: &[L::State],
        reference_density_kg_m3: &[f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<BackwardEulerHexStep<L::State>, SmallStrainTetError> {
        self.advance_with_density(
            displacement_m,
            velocity_mps,
            quadrature_states,
            ReferenceDensity::Cell(reference_density_kg_m3),
            dt_s,
            cx,
        )
    }

    /// The same motion with a positive original density at each cell corner,
    /// in the element's connectivity order. Trilinear interpolation occurs
    /// only inside that original cell; distinct neighboring material profiles
    /// are not averaged. The separate mass rule integrates their full spatial
    /// distribution, including translation/deformation and support coupling.
    pub fn advance_backward_euler_with_corner_density(
        &self,
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        quadrature_states: &[L::State],
        reference_corner_density_kg_m3: &[[f64; 8]],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<BackwardEulerHexStep<L::State>, SmallStrainTetError> {
        self.advance_with_density(
            displacement_m,
            velocity_mps,
            quadrature_states,
            ReferenceDensity::Corners(reference_corner_density_kg_m3),
            dt_s,
            cx,
        )
    }

    fn advance_with_density(
        &self,
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        quadrature_states: &[L::State],
        density: ReferenceDensity<'_>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<BackwardEulerHexStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        if self.load_path != [1.0] {
            return Err(invalid(
                "one endpoint load at factor one is required per hex inertial interval",
            ));
        }
        let elements = self.prepare(cx)?;
        self.advance_prepared(
            &elements,
            displacement_m,
            velocity_mps,
            quadrature_states,
            density,
            dt_s,
            cx,
        )
    }

    /// Reuse geometry already admitted for this exact one-endpoint problem.
    /// Thermal history admission and the inertial solve share that preparation.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_prepared(
        &self,
        elements: &[Element],
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        quadrature_states: &[L::State],
        density: ReferenceDensity<'_>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<BackwardEulerHexStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        self.check(displacement_m, quadrature_states)?;
        let (count, values) = match density {
            ReferenceDensity::Cell(values) => (values.len(), values),
            ReferenceDensity::Corners(values) => (values.len(), values.as_flattened()),
        };
        if count != elements.len()
            || values
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid(
                "one positive finite reference density or corner profile per hex required",
            ));
        }
        // Material contributes at most 24^2 entries; consistent inertia adds
        // three 8^2 blocks per original cell before shared-node assembly.
        let entries = elements
            .len()
            .checked_mul(24 * 24 + 3 * 8 * 8)
            .ok_or_else(|| invalid("hex inertial tangent count overflow"))?;
        if entries > self.budget.maximum_tangent_entries {
            return Err(invalid(format!(
                "hex inertial tangent entries {entries} exceed budget {}",
                self.budget.maximum_tangent_entries
            )));
        }
        let (coordinate, weight) = gauss_legendre(3);
        let geometry = elements.iter().enumerate().map(|(index, element)| {
            let positions = element.nodes.map(|node| self.nodes_m[node]);
            let (uniform, profile) = match density {
                ReferenceDensity::Cell(values) => (Some(values[index]), [values[index]; 8]),
                ReferenceDensity::Corners(values) => {
                    let profile = values[index];
                    (
                        profile
                            .iter()
                            .all(|value| *value == profile[0])
                            .then_some(profile[0]),
                        profile,
                    )
                }
            };
            let mut points = [([0.0; 8], 0.0, 0.0); 27];
            for (q, entry) in points.iter_mut().enumerate() {
                cx.checkpoint().map_err(|_| InertiaError::Cancelled)?;
                let [i, j, k] = [q / 9, (q / 3) % 3, q % 3];
                let point = hex8_point(
                    positions,
                    [coordinate[i], coordinate[j], coordinate[k]],
                    self.budget.minimum_scaled_jacobian,
                    index,
                    q,
                )
                .map_err(|error| InertiaError::InvalidInput {
                    what: error.to_string(),
                })?;
                let rho =
                    uniform.unwrap_or_else(|| (0..8).map(|a| point.shapes[a] * profile[a]).sum());
                *entry = (
                    point.shapes,
                    point.determinant_m3 * weight[i] * weight[j] * weight[k],
                    rho,
                );
            }
            Ok((element.nodes, points))
        });
        let mass = ReferenceMass::from_hex8_density_quadrature(geometry, displacement_m.len(), cx)?;
        let inertia = BackwardEulerInertia::new(
            &mass,
            displacement_m,
            velocity_mps.as_flattened(),
            dt_s,
            cx,
        )?;
        let solution = self.solve_prepared(
            elements,
            Some((displacement_m, quadrature_states)),
            Some(&inertia),
            cx,
        )?;
        let (velocity_mps, energy) = inertia.motion(solution.displacement_m.as_flattened(), cx)?;
        checkpoint(cx)?;
        Ok(BackwardEulerHexStep {
            solution,
            velocity_mps,
            energy,
        })
    }
}
