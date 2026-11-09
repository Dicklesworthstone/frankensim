//! Full consistent reference-mass action for physical acceleration loads.
//!
//! Reuse the element owner's admitted geometry, positive density and mass rule.
//! Interpolate acceleration before testing against each native shape: this is
//! the same integral as M*a without assembling stiffness or a global matrix.

use super::Point;
use crate::linear3::{TetElasticError, TetMaterialField};
use fs_exec::Cx;

/// Prescribed reference body acceleration [m/s²] in the mesh's Cartesian axes.
///
/// This is a held dead-load field. It does not infer a gravity, centrifugal or
/// follower-force law or resample acceleration after deformation. The original
/// reference density supplies inertia independently of the acceleration field.
#[derive(Debug, Clone, Copy)]
pub enum ReferenceAccelerationField<'a> {
    /// Spatially uniform acceleration, including ordinary uniform gravity.
    Uniform([f64; 3]),
    /// One vector per original geometry node, interpolated by the native
    /// displacement shapes. This includes every supported node; a reduced
    /// free-DOF field would lose the full consistent mass's supported columns.
    PerNode(&'a [[f64; 3]]),
}

impl ReferenceAccelerationField<'_> {
    fn validate(self, node_count: usize, cx: &Cx<'_>) -> Result<(), TetElasticError> {
        checkpoint(cx)?;
        match self {
            Self::Uniform(value) => {
                if value.iter().any(|value| !value.is_finite()) {
                    return Err(invalid("uniform acceleration must be finite"));
                }
            }
            Self::PerNode(values) => {
                if values.len() != node_count {
                    return Err(invalid(
                        "one acceleration vector per original geometry node required",
                    ));
                }
                for (node, value) in values.iter().enumerate() {
                    if node % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    if value.iter().any(|value| !value.is_finite()) {
                        return Err(invalid("nodal acceleration must be finite"));
                    }
                }
            }
        }
        Ok(())
    }

    fn at<const N: usize>(self, nodes: &[usize; N], shapes: &[f64; N]) -> [f64; 3] {
        match self {
            Self::Uniform(value) => value,
            Self::PerNode(values) => {
                // Partition of unity expressed relative to one local sample
                // retains an exactly uniform acceleration without summing
                // signed quadratic shape values against its common datum.
                let anchor = values[nodes[0]];
                let mut value = anchor;
                for a in 1..N {
                    for axis in 0..3 {
                        value[axis] += shapes[a] * (values[nodes[a]][axis] - anchor[axis]);
                    }
                }
                value
            }
        }
    }
}

/// Private mass action on already admitted native connectivity and density.
/// The caller's original node/element/free-DOF budgets remain the work envelope.
/// Full output rows include supports and are sorted for the existing dead-load
/// equilibrium/reaction owners; no lumping or constraint projection occurs.
#[allow(clippy::too_many_arguments)]
pub(in crate::linear3) fn assemble<const N: usize, R: Copy>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    acceleration: ReferenceAccelerationField<'_>,
    mass_rule: &[(R, f64)],
    evaluate: impl Fn([[f64; 3]; N], R, usize, usize) -> Result<Point<N>, TetElasticError>,
    reference_density: impl Fn(usize, R) -> Option<f64>,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, f64)>, TetElasticError> {
    acceleration.validate(nodes_m.len(), cx)?;
    let full_dofs = nodes_m
        .len()
        .checked_mul(3)
        .ok_or(TetElasticError::SizeOverflow)?;
    let mut force = Vec::with_capacity(full_dofs);
    for node in 0..nodes_m.len() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for axis in 0..3 {
            force.push((3 * node + axis, 0.0));
        }
    }
    let mut total_mass = 0.0;
    for (element, connectivity) in elements.iter().enumerate() {
        checkpoint(cx)?;
        let nodes = connectivity.map(|node| nodes_m[node]);
        let material = materials.at(element);
        for (point, &(reference, weight)) in mass_rule.iter().enumerate() {
            checkpoint(cx)?;
            let p = evaluate(nodes, reference, element, point)?;
            let density = reference_density(element, reference).unwrap_or(material.density_kg_m3);
            let dv = weight * p.determinant;
            let dm = density * dv;
            if !density.is_finite() || density <= 0.0 || !dm.is_finite() || dm <= 0.0 {
                return Err(TetElasticError::InvalidMaterial {
                    what: "reference acceleration requires representable positive quadrature mass",
                });
            }
            total_mass += dm;
            if !total_mass.is_finite() {
                return Err(TetElasticError::InvalidMaterial {
                    what: "reference acceleration total mass is not representable",
                });
            }
            let value = acceleration.at(connectivity, &p.shapes);
            if value.iter().any(|value| !value.is_finite()) {
                return Err(invalid("interpolated reference acceleration is not finite"));
            }
            for a in 0..N {
                let nodal_mass = dm * p.shapes[a];
                for axis in 0..3 {
                    let dof = 3 * connectivity[a] + axis;
                    force[dof].1 += nodal_mass * value[axis];
                    if !force[dof].1.is_finite() {
                        return Err(invalid(
                            "consistent reference-acceleration force is not finite",
                        ));
                    }
                }
            }
        }
    }
    checkpoint(cx)?;
    Ok(force)
}

fn invalid(what: &'static str) -> TetElasticError {
    TetElasticError::InvalidBodyAcceleration { what }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TetElasticError> {
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)
}
