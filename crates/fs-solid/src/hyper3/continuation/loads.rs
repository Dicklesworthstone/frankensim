//! Original configuration-dependent loading for the native continuation owner.

use fs_exec::Cx;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_sparse::Csr;

use super::{HyperTetContinuationError, invalid};
use crate::SolidError;
use crate::hyper3::gravity::{GravityField3, ReferenceGravityDensity, ReferenceGravityWeights};
use crate::hyper3::prepared::PreparedHyperTetProblem;
use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use crate::hyper3::{HyperTetError, HyperTetProblem, checkpoint};

/// Optional native loads in addition to the problem's original dead forces.
///
/// The path parameter multiplies their forces at the CURRENT configuration.
/// Pressure directions, areas, spatial origins and free-surface planes remain
/// native; a negative parameter reverses the held force, not its spatial datum.
/// Gravity retains fixed reference mass and the original stationary field.
#[derive(Clone, Copy, Default)]
pub struct HyperTetPathLoads<'a> {
    /// Exact reference surface and held current-area pressure field.
    /// Nonuniform face patches retain their nonsymmetric follower tangent.
    pub pressure: Option<(&'a ClosedTetPressureSurface, TetPressure<'a>)>,
    /// Fixed reference density basis and stationary external gravity field.
    pub gravity: Option<(ReferenceGravityDensity<'a>, &'a GravityField3<'a>)>,
}

impl core::fmt::Debug for HyperTetPathLoads<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("HyperTetPathLoads")
            .field("pressure", &self.pressure)
            .field("gravity_density", &self.gravity.map(|(density, _)| density))
            .finish_non_exhaustive()
    }
}

pub(super) struct PreparedPathLoads<'a> {
    pub(super) input: HyperTetPathLoads<'a>,
    unloaded_pressure: bool,
    gravity_weights: Option<ReferenceGravityWeights>,
}

impl core::fmt::Debug for PreparedPathLoads<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedPathLoads")
            .field("input", &self.input)
            .field("unloaded_pressure", &self.unloaded_pressure)
            .finish_non_exhaustive()
    }
}

impl<'a> PreparedPathLoads<'a> {
    pub(super) fn prepare<L: HyperelasticMaterial>(
        problem: &HyperTetProblem<'a, L>,
        prepared: &PreparedHyperTetProblem<'a, L>,
        input: HyperTetPathLoads<'a>,
        cx: &Cx<'_>,
    ) -> Result<Option<Self>, HyperTetError> {
        checkpoint(cx)?;
        if input.pressure.is_none() && input.gravity.is_none() {
            return Ok(None);
        }
        let unloaded_pressure = if let Some((surface, pressure)) = input.pressure {
            if !surface.matches_reference(problem.nodes_m, problem.tetrahedra) {
                return Err(HyperTetError::InvalidInput {
                    what: "continuation pressure requires the exact reference surface".into(),
                });
            }
            surface.validate_pressure(pressure, cx)?
        } else {
            true
        };
        let gravity_weights = input
            .gravity
            .map(|(density, _)| prepared.bind_gravity_weights(density, cx))
            .transpose()?;
        checkpoint(cx)?;
        Ok(Some(Self {
            input,
            unloaded_pressure,
            gravity_weights,
        }))
    }

    pub(super) fn add_force<L: HyperelasticMaterial>(
        &self,
        prepared: &PreparedHyperTetProblem<'_, L>,
        displacement: &[f64],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        if let Some((surface, pressure)) = self.input.pressure {
            if !self.unloaded_pressure {
                add_force(force, &surface.force_n(displacement, pressure, cx)?, cx)?;
            }
        }
        if let Some((_, field)) = self.input.gravity {
            let gravity = prepared.gravity_evaluate_with_weights(
                self.gravity_weights
                    .as_ref()
                    .expect("bound gravity inventory"),
                displacement,
                field,
                false,
                cx,
            )?;
            add_force(force, &gravity.force_n, cx)?;
        }
        checkpoint(cx)
    }

    pub(super) fn tangents<L: HyperelasticMaterial>(
        &self,
        prepared: &PreparedHyperTetProblem<'_, L>,
        displacement: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Option<Csr>, Option<Csr>), HyperTetError> {
        let pressure = self
            .input
            .pressure
            .filter(|_| !self.unloaded_pressure)
            .map(|(surface, pressure)| surface.force_tangent(displacement, pressure, cx))
            .transpose()?;
        let gravity = if let Some((_, field)) = self.input.gravity {
            prepared
                .gravity_evaluate_with_weights(
                    self.gravity_weights
                        .as_ref()
                        .expect("bound gravity inventory"),
                    displacement,
                    field,
                    true,
                    cx,
                )?
                .tangent
        } else {
            None
        };
        checkpoint(cx)?;
        Ok((pressure, gravity))
    }
}

fn add_force(target: &mut [f64], load: &[f64], cx: &Cx<'_>) -> Result<(), HyperTetError> {
    for (dof, (force, addition)) in target.iter_mut().zip(load).enumerate() {
        if dof % 1024 == 0 {
            checkpoint(cx)?;
        }
        *force += addition;
        if !force.is_finite() {
            return Err(HyperTetError::InvalidInput {
                what: "combined continuation load is not representable".into(),
            });
        }
    }
    Ok(())
}

/// Merge canonical native rows directly; preserve the nonsymmetric Jacobian.
/// No global COO buffer, sort or dense intermediate is needed for load stiffness.
#[allow(clippy::too_many_arguments)]
pub(super) fn restrict_tangent(
    material: &Csr,
    pressure: Option<&Csr>,
    gravity: Option<&Csr>,
    lambda: f64,
    scale: f64,
    free: &[usize],
    reduced: &[usize],
    cx: &Cx<'_>,
) -> Result<Csr, HyperTetContinuationError> {
    let n = free.len();
    let mut rows = Vec::with_capacity(n + 1);
    let mut columns = Vec::new();
    let mut values = Vec::new();
    rows.push(0);
    for (row, &dof) in free.iter().enumerate() {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        let (ac, av) = material.row(dof);
        let (pc, pv) = pressure.map_or((&[][..], &[][..]), |matrix| matrix.row(dof));
        let (gc, gv) = gravity.map_or((&[][..], &[][..]), |matrix| matrix.row(dof));
        let (mut a, mut p, mut g) = (0, 0, 0);
        while a < ac.len() || p < pc.len() || g < gc.len() {
            let column = ac
                .get(a)
                .copied()
                .unwrap_or(usize::MAX)
                .min(pc.get(p).copied().unwrap_or(usize::MAX))
                .min(gc.get(g).copied().unwrap_or(usize::MAX));
            let mut value = if ac.get(a) == Some(&column) {
                let value = av[a];
                a += 1;
                value
            } else {
                0.0
            };
            if pc.get(p) == Some(&column) {
                value -= lambda * pv[p];
                p += 1;
            }
            if gc.get(g) == Some(&column) {
                value -= lambda * gv[g];
                g += 1;
            }
            let free_column = reduced[column];
            if free_column != usize::MAX {
                let scaled = value * scale;
                if !scaled.is_finite() {
                    return Err(invalid("normalized loaded tangent is not finite"));
                }
                columns.push(free_column);
                values.push(scaled);
            }
        }
        rows.push(values.len());
    }
    let tangent =
        Csr::try_from_parts_with_checkpoint(n, n, rows, columns, values, || checkpoint(cx))?
            .ok_or_else(|| SolidError::InternalInvariant {
                what: "canonical loaded tangent lost ordering during free-DOF restriction".into(),
            })?;
    checkpoint(cx)?;
    Ok(tangent)
}
