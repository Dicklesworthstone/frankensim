//! Moving transverse forces, bending moments, and point probes on a flat plate.
//!
//! A located stencil uses the containing triangle's P1 barycentric weights for
//! all three nodal fields `(w, wx, wy)`. This is an explicit interface
//! interpolation, not the DKT element's quadratic interior rotation field.
//! Applying the transpose of that same interpolation preserves virtual work
//! and instantaneous power. Physical rotations obey `theta_x = wy` and
//! `theta_y = -wx`, so a Cartesian load `[Fz, Mx, My]` contributes generalized
//! nodal loads `[Fz, -My, Mx]` times each node's weight.
//!
//! Locate at each current path position, then reuse the stencil for load,
//! probe, and load-parameter pullback. Mesh coordinates and support indexing
//! are snapshots: changing either requires a new stencil. Positions are in
//! the undeformed xy plane, with zero prescribed motion at eliminated DOFs.
//! This module does not differentiate the moving path or model contact.

use crate::{PlateMesh, PlateModel};

/// Admission bounds, checked before scanning mesh geometry or the DOF map.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlateLoadBudget {
    /// Maximum number of mesh nodes to inspect.
    pub max_nodes: usize,
    /// Maximum number of triangles to inspect.
    pub max_triangles: usize,
}

/// A bounded point-transfer refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PlateLoadError {
    /// The caller's bound cannot cover admission and location work.
    Budget {
        /// Bounded item: nodes or triangles.
        what: &'static str,
        /// Required count.
        required: usize,
        /// Caller-supplied maximum.
        limit: usize,
    },
    /// Empty mesh, invalid connectivity, coordinates, or support mapping.
    InvalidMesh(&'static str),
    /// Degenerate, inverted, or unrepresentably ill-conditioned triangle.
    InvalidTriangle {
        /// Index in `PlateMesh::tris`.
        triangle: usize,
    },
    /// No triangle contains the point, including the mesh boundary.
    Outside {
        /// Requested undeformed point in metres.
        point: [f64; 2],
    },
    /// A reduced state, seed, or load buffer has the wrong length.
    Dimension {
        /// Reduced model dimension.
        expected: usize,
        /// Supplied length.
        actual: usize,
    },
    /// Nonfinite input or an unrepresentable accumulation.
    NonFinite(&'static str),
}

impl core::fmt::Display for PlateLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-PLATE-LOAD: {self:?}")
    }
}

impl std::error::Error for PlateLoadError {}

/// The applied resultant and the portion on eliminated DOFs.
#[derive(Debug, Clone, PartialEq)]
pub struct PlateAppliedLoad {
    /// `(full_global_dof, generalized_load)` pairs on fixed DOFs only.
    /// At most nine entries; these are applied loads, not complete support
    /// reactions (which also require inertia, damping, and internal forces).
    pub fixed_dof_loads: Vec<(usize, f64)>,
    /// Total Cartesian force `[0, 0, Fz]` in newtons, including fixed DOFs.
    pub force: [f64; 3],
    /// Total Cartesian moment about the undeformed global origin, in N m:
    /// `[Mx + y*Fz, My - x*Fz, 0]`, including fixed DOFs.
    pub moment: [f64; 3],
}

/// Reconstructed point motion, using the same weights as load transfer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlatePointSample {
    /// Transverse displacement in metres.
    pub displacement: f64,
    /// The interpolated nodal slope fields `[wx, wy]` in radians.
    /// They need not equal the gradient of the interpolated P1 displacement.
    pub slopes: [f64; 2],
    /// Transverse velocity in m/s.
    pub velocity: f64,
    /// Physical Cartesian angular velocity `[v_wy, -v_wx, 0]` in rad/s.
    pub angular_velocity: [f64; 3],
}

/// One admitted triangle location and its immutable reduced-index snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct PlatePointStencil {
    triangle: usize,
    nodes: [usize; 3],
    weights: [f64; 3],
    point: [f64; 2],
    reduced: [[Option<usize>; 3]; 3],
    free: usize,
}

impl PlatePointStencil {
    /// Locate a point without nearest-node snapping.
    ///
    /// All geometry and support-map entries are checked within the explicit
    /// bounds. CCW triangles must have normalized doubled area greater than
    /// `64*EPSILON`. Boundary barycentric coordinates admit roundoff within
    /// `64*EPSILON`, clamp that roundoff, and renormalize. On an edge or shared
    /// vertex, the lowest containing triangle index wins deterministically.
    /// Mesh overlap is not inferred or repaired.
    ///
    /// # Errors
    /// Refuses exhausted bounds, invalid mesh/index maps, nonfinite positions,
    /// degenerate/inverted triangles, or a point outside the mesh.
    pub fn locate(
        mesh: &PlateMesh,
        model: &PlateModel,
        point: [f64; 2],
        budget: PlateLoadBudget,
    ) -> Result<Self, PlateLoadError> {
        admit(mesh, model, point, budget)?;
        let mut located = None;
        for (triangle, &nodes) in mesh.tris.iter().enumerate() {
            let geometry = TriangleGeometry::new(mesh, nodes, triangle)?;
            if located.is_none()
                && let Some(weights) = geometry.weights(point)
            {
                located = Some(Self {
                    triangle,
                    nodes,
                    weights,
                    point,
                    reduced: nodes.map(|node| {
                        std::array::from_fn(|component| model.dof_map[3 * node + component])
                    }),
                    free: model.free,
                });
            }
        }
        located.ok_or(PlateLoadError::Outside { point })
    }

    /// Index of the containing triangle.
    #[must_use]
    pub fn triangle(&self) -> usize {
        self.triangle
    }

    /// P1 weights in the containing triangle's node order.
    #[must_use]
    pub fn weights(&self) -> [f64; 3] {
        self.weights
    }

    /// Requested point in the undeformed global xy plane, in metres.
    #[must_use]
    pub fn point(&self) -> [f64; 2] {
        self.point
    }

    /// Add `[Fz, Mx, My]` to an existing reduced force vector.
    ///
    /// Load units are N, N m, N m. The reduced vector is unchanged on any
    /// refusal. Repeated calls accumulate loads, so a caller can represent a
    /// distributed patch with its own explicitly bounded quadrature stencil.
    /// Fixed-DOF contributions are returned separately and include zeros.
    ///
    /// # Errors
    /// Refuses a dimension mismatch, nonfinite load, or overflow in the
    /// resultant or any touched reduced-vector entry.
    pub fn add_load(
        &self,
        load: [f64; 3],
        reduced: &mut [f64],
    ) -> Result<PlateAppliedLoad, PlateLoadError> {
        self.check_dimension(reduced)?;
        if load.iter().any(|value| !value.is_finite()) {
            return Err(PlateLoadError::NonFinite("point load"));
        }
        let [fz, mx, my] = load;
        let moment = [mx + self.point[1] * fz, my - self.point[0] * fz, 0.0];
        if moment.iter().any(|value| !value.is_finite()) {
            return Err(PlateLoadError::NonFinite("global moment resultant"));
        }
        let generalized = [fz, -my, mx];
        let mut applied = PlateAppliedLoad {
            fixed_dof_loads: Vec::new(),
            force: [0.0, 0.0, fz],
            moment,
        };
        let mut updates = [[0.0; 3]; 3];
        for (local, (&weight, mapping)) in self.weights.iter().zip(&self.reduced).enumerate() {
            for (component, (&effort, &index)) in generalized.iter().zip(mapping).enumerate() {
                let contribution = weight * effort;
                if let Some(index) = index {
                    updates[local][component] = reduced[index] + contribution;
                    if !updates[local][component].is_finite() {
                        return Err(PlateLoadError::NonFinite("reduced load accumulation"));
                    }
                } else {
                    applied
                        .fixed_dof_loads
                        .push((3 * self.nodes[local] + component, contribution));
                }
            }
        }
        for (mapping, updates) in self.reduced.iter().zip(updates) {
            for (&index, value) in mapping.iter().zip(updates) {
                if let Some(index) = index {
                    reduced[index] = value;
                }
            }
        }
        Ok(applied)
    }

    /// Interpolate displacement, slopes, transverse and angular velocity.
    ///
    /// # Errors
    /// Refuses dimension mismatch or nonfinite contributing/reconstructed
    /// values. Eliminated DOFs have zero prescribed displacement and velocity.
    pub fn sample(&self, q: &[f64], v: &[f64]) -> Result<PlatePointSample, PlateLoadError> {
        let q = self.interpolate(q)?;
        let v = self.interpolate(v)?;
        Ok(PlatePointSample {
            displacement: q[0],
            slopes: [q[1], q[2]],
            velocity: v[0],
            angular_velocity: [v[2], -v[1], 0.0],
        })
    }

    /// Pull a reduced-force seed back to physical `[Fz, Mx, My]`.
    ///
    /// This is the exact load derivative at this fixed point and support
    /// mapping. For load-amplitude parameters, dot it with the physical-load
    /// derivative; it does not include derivatives of the path coordinates.
    ///
    /// # Errors
    /// Refuses a dimension mismatch or nonfinite contributing/returned values.
    pub fn load_vjp(&self, seed: &[f64]) -> Result<[f64; 3], PlateLoadError> {
        let field = self.interpolate(seed)?;
        Ok([field[0], field[2], -field[1]])
    }

    fn check_dimension(&self, field: &[f64]) -> Result<(), PlateLoadError> {
        if field.len() != self.free {
            return Err(PlateLoadError::Dimension {
                expected: self.free,
                actual: field.len(),
            });
        }
        Ok(())
    }

    fn interpolate(&self, field: &[f64]) -> Result<[f64; 3], PlateLoadError> {
        self.check_dimension(field)?;
        let mut values = [0.0; 3];
        for (&weight, mapping) in self.weights.iter().zip(&self.reduced) {
            for (value, &index) in values.iter_mut().zip(mapping) {
                if let Some(index) = index {
                    if !field[index].is_finite() {
                        return Err(PlateLoadError::NonFinite("point field"));
                    }
                    *value += weight * field[index];
                }
            }
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(PlateLoadError::NonFinite("interpolated point field"));
        }
        Ok(values)
    }
}

fn admit(
    mesh: &PlateMesh,
    model: &PlateModel,
    point: [f64; 2],
    budget: PlateLoadBudget,
) -> Result<(), PlateLoadError> {
    for (what, required, limit) in [
        ("nodes", mesh.nodes.len(), budget.max_nodes),
        ("triangles", mesh.tris.len(), budget.max_triangles),
    ] {
        if required > limit {
            return Err(PlateLoadError::Budget {
                what,
                required,
                limit,
            });
        }
    }
    if point.iter().any(|value| !value.is_finite()) {
        return Err(PlateLoadError::NonFinite("point coordinates"));
    }
    if mesh.nodes.len() < 3 || mesh.tris.is_empty() {
        return Err(PlateLoadError::InvalidMesh("empty or undersized mesh"));
    }
    if mesh
        .nodes
        .iter()
        .any(|&(x, y)| !x.is_finite() || !y.is_finite())
    {
        return Err(PlateLoadError::InvalidMesh("nonfinite node coordinates"));
    }
    let full = mesh
        .nodes
        .len()
        .checked_mul(3)
        .ok_or(PlateLoadError::InvalidMesh("full DOF count overflows"))?;
    if model.dof_map.len() != full || model.free > full {
        return Err(PlateLoadError::InvalidMesh("support map dimension"));
    }
    let mut seen = vec![false; model.free];
    for &index in model.dof_map.iter().flatten() {
        if index >= model.free || seen[index] {
            return Err(PlateLoadError::InvalidMesh("support map is not one-to-one"));
        }
        seen[index] = true;
    }
    if seen.iter().any(|&present| !present) {
        return Err(PlateLoadError::InvalidMesh(
            "support map has missing free DOFs",
        ));
    }
    Ok(())
}

const GEOMETRY_TOLERANCE: f64 = 64.0 * f64::EPSILON;

struct TriangleGeometry {
    origin: [f64; 2],
    edge1: [f64; 2],
    edge2: [f64; 2],
    scale: f64,
    determinant: f64,
}

impl TriangleGeometry {
    fn new(mesh: &PlateMesh, nodes: [usize; 3], triangle: usize) -> Result<Self, PlateLoadError> {
        if nodes.iter().any(|&node| node >= mesh.nodes.len()) {
            return Err(PlateLoadError::InvalidTriangle { triangle });
        }
        let (x, y) = mesh.nodes[nodes[0]];
        let (x1, y1) = mesh.nodes[nodes[1]];
        let (x2, y2) = mesh.nodes[nodes[2]];
        let scale = (x1 - x)
            .abs()
            .max((y1 - y).abs())
            .max((x2 - x).abs())
            .max((y2 - y).abs());
        if !(scale > 0.0 && scale.is_finite()) {
            return Err(PlateLoadError::InvalidTriangle { triangle });
        }
        let edge1 = [(x1 - x) / scale, (y1 - y) / scale];
        let edge2 = [(x2 - x) / scale, (y2 - y) / scale];
        let determinant = cross(edge1, edge2);
        if determinant <= GEOMETRY_TOLERANCE {
            return Err(PlateLoadError::InvalidTriangle { triangle });
        }
        Ok(Self {
            origin: [x, y],
            edge1,
            edge2,
            scale,
            determinant,
        })
    }

    fn weights(&self, point: [f64; 2]) -> Option<[f64; 3]> {
        let relative = [
            (point[0] - self.origin[0]) / self.scale,
            (point[1] - self.origin[1]) / self.scale,
        ];
        let weight1 = cross(relative, self.edge2) / self.determinant;
        let weight2 = cross(self.edge1, relative) / self.determinant;
        let mut weights = [1.0 - weight1 - weight2, weight1, weight2];
        if weights
            .iter()
            .any(|&value| !value.is_finite() || value < -GEOMETRY_TOLERANCE)
        {
            return None;
        }
        for weight in &mut weights {
            *weight = weight.max(0.0);
        }
        let sum: f64 = weights.iter().sum();
        for weight in &mut weights {
            *weight /= sum;
        }
        Some(weights)
    }
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
