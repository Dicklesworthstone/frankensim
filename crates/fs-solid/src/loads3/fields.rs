//! Consistent loads from piecewise-linear fields on the retained P1 mesh.
//!
//! These are field coefficients, not point-load forces: pressure/traction
//! samples are in Pa and body-force-density samples in N/m^3. Each selected
//! triangle retains its own three samples, allowing discontinuous patch data.
//! Body data may be continuous nodal values or four independent values per tet.
//! Integration uses the exact affine-simplex mass matrix, including off-diagonal
//! terms. Replacing it with equal nodal shares preserves force but generally
//! destroys first moments and work against an affine displacement.

use super::*;

/// P1 surface fields in the canonical node order of each selected `LoadFace3`.
/// The outer value array follows `faces`; shared vertices need not share values.
#[derive(Debug, Clone, Copy)]
pub enum SurfaceField3<'a> {
    /// Three Cartesian traction vectors per selected face [Pa].
    Traction { faces: &'a [usize], nodal_pascals: &'a [[[f64; 3]; 3]] },
    /// Three inward-positive scalar pressures per selected face [Pa].
    Pressure { faces: &'a [usize], nodal_pascals: &'a [[f64; 3]] },
}

/// P1 body-force density [N/m^3], never a mass density or acceleration.
#[derive(Debug, Clone, Copy)]
pub enum VolumeField3<'a> {
    None,
    /// One vector per geometry node; the interpolant is continuous across faces.
    PerNode(&'a [[f64; 3]]),
    /// Four vectors per tetrahedron in its original connectivity order.
    /// Separate element samples permit material-interface jumps without averaging.
    PerElement(&'a [[[f64; 3]; 4]]),
}

/// Integral N_i N_j = measure * (1 + delta_ij) / (N * (N + 1)), where
/// N is the number of vertices (three on triangles, four on tetrahedra).
fn scatter<const N: usize>(
    forces: &mut [f64], nodes: [usize; N], values: [[f64; 3]; N], measure: f64,
) -> Result<(), LoadError3> {
    for value in values { vector(value)?; }
    let unit = finite(measure / (N * (N + 1)) as f64)?;
    for i in 0..N {
        for j in 0..N {
            let weight = if i == j { 2.0 * unit } else { unit };
            for c in 0..3 {
                let contribution = finite(weight * values[j][c])?;
                let dof = 3 * nodes[i] + c;
                forces[dof] = finite(forces[dof] + contribution)?;
            }
        }
    }
    Ok(())
}

impl<'m> TetLoadGeometry3<'m> {
    /// Integrate the declared P1 pressure, traction and body fields into the
    /// same `TetLoadCase3` consumed by native equilibrium/stress/reaction recovery.
    /// The interpolation is exactly the supplied affine-element field; no
    /// reconstruction of an unknown source function or quadrature-error bound
    /// is implied. It is not a follower force on deformed geometry.
    ///
    /// Surface selection must be nonempty, exterior and assigned at most once.
    /// Count mismatch, nonfinite coefficients/arithmetic or cancellation returns
    /// no partial case. Work is bounded by the admitted geometry and checked
    /// per face/element; each 3x3 or 4x4 field-mass contraction is indivisible.
    /// Static nodal forces are superposed in their usual sorted, unique order.
    /// For a time-dependent combination, pass this case's `nodal_forces_n()`
    /// as the static nodal component of `history::TetLoadHistory3`.
    pub fn assemble_p1_fields<'g>(
        &'g self, surfaces: &[SurfaceField3<'_>], body: VolumeField3<'_>,
        nodal: &[(usize, f64)], cx: &Cx<'_>,
    ) -> Result<TetLoadCase3<'g, 'm>, LoadError3> {
        poll(cx)?;
        let full_dofs = self.nodes.len() * 3;
        if surfaces.len() > self.faces.len() || nodal.len() > full_dofs {
            return Err(LoadError3::Budget);
        }
        match body {
            VolumeField3::PerNode(values) if values.len() != self.nodes.len() => {
                return Err(LoadError3::Invalid("one P1 body-force-density vector per mesh node required"));
            }
            VolumeField3::PerElement(values) if values.len() != self.tets.len() => {
                return Err(LoadError3::Invalid("one four-vector P1 body field per tetrahedron required"));
            }
            _ => {}
        }
        let mut used = vec![false; self.faces.len()];
        // Validate the entire surface partition before any force accumulation.
        for row in surfaces {
            poll(cx)?;
            let (faces, count) = match *row {
                SurfaceField3::Traction { faces, nodal_pascals } => (faces, nodal_pascals.len()),
                SurfaceField3::Pressure { faces, nodal_pascals } => (faces, nodal_pascals.len()),
            };
            if faces.is_empty() || faces.len() > self.faces.len() || faces.len() != count {
                return Err(LoadError3::Invalid("one P1 field per selected exterior face required"));
            }
            for &face in faces {
                poll(cx)?;
                let taken = used.get_mut(face).ok_or(LoadError3::Invalid("unknown exterior face index"))?;
                if std::mem::replace(taken, true) {
                    return Err(LoadError3::Invalid("duplicate surface load assignment"));
                }
            }
        }
        let mut forces = vec![0.0; full_dofs];
        for row in surfaces {
            let faces = match *row {
                SurfaceField3::Traction { faces, .. } | SurfaceField3::Pressure { faces, .. } => faces,
            };
            for (index, &id) in faces.iter().enumerate() {
                poll(cx)?;
                let face = self.faces[id];
                let values = match *row {
                    SurfaceField3::Traction { nodal_pascals, .. } => nodal_pascals[index],
                    SurfaceField3::Pressure { nodal_pascals, .. } => {
                        nodal_pascals[index].map(|p| face.outward_normal.map(|n| -p * n))
                    }
                };
                scatter(&mut forces, face.nodes, values, face.area_m2)?;
            }
        }
        if !matches!(body, VolumeField3::None) {
            for (index, &tet) in self.tets.iter().enumerate() {
                poll(cx)?;
                let values = match body {
                    VolumeField3::None => unreachable!("no body field is excluded above"),
                    VolumeField3::PerNode(values) => tet.map(|node| values[node]),
                    VolumeField3::PerElement(values) => values[index],
                };
                scatter(&mut forces, tet, values, self.volumes[index])?;
            }
        }
        self.finish_forces(forces, nodal, cx)
    }
}
