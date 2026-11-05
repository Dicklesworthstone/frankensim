//! Tangential Fourier transport in thermally thin, perfectly bonded coatings.
//!
//! The film shares the original solid trace temperature. Its declared sheet
//! conductance is k*t [W/K]; its heat flux is per edge length [W/m]. It adds
//! neither storage nor through-thickness resistance. Keep any film storage in
//! the original surface-capacity owner and omit a separately meshed layer.

use std::collections::BTreeSet;

use fs_exec::Cx;

use super::ThermalBoundary;
use crate::{ConductionError, ConductionMesh};

/// Temperature-independent sheet conductance on one original face.
/// Coefficients follow that face's vertex order and are interpolated with P1
/// basis functions. Tensor laws use an explicitly supplied orthonormal tangent
/// frame; they are intrinsic sheet laws, not unqualified bulk tensors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceHeatConduction {
    boundary_face: usize,
    law: SheetLaw,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum SheetLaw {
    Isotropic([f64; 3]),
    Tensor {
        axes: [[f64; 3]; 2],
        nodal: [[[f64; 2]; 2]; 3],
    },
}

impl SurfaceHeatConduction {
    /// Declare uniform k*t [W/K], finite and nonnegative. Zero adds no film.
    pub fn uniform(
        boundary_face: usize,
        conductance_w_per_k: f64,
    ) -> Result<Self, ConductionError> {
        Self::nodal(boundary_face, [conductance_w_per_k; 3])
    }

    /// Declare a P1 k*t field in the original face's three-vertex order.
    /// These are sheet conductances, not separately interpolated k and t.
    pub fn nodal(boundary_face: usize, nodal_w_per_k: [f64; 3]) -> Result<Self, ConductionError> {
        for value in nodal_w_per_k {
            crate::require_finite("sheet conductance", value)?;
            if value < 0.0 {
                return Err(invalid("sheet conductance must be nonnegative"));
            }
        }
        Ok(Self {
            boundary_face,
            law: SheetLaw::Isotropic(nodal_w_per_k),
        })
    }

    /// Declare a uniform symmetric positive-semidefinite 2-D sheet tensor
    /// [W/K] in two orthonormal world-space tangent axes. Zero and rank-one
    /// sheets are admitted. The axes must lie in the original face plane;
    /// binding to a boundary or evaluating flux checks that geometric fact.
    ///
    /// This is an already reduced in-plane law. It does not discard the normal
    /// coupling of a bulk anisotropic tensor or infer its thin-layer reduction.
    pub fn uniform_tensor(
        boundary_face: usize,
        tangent_axes: [[f64; 3]; 2],
        conductance_w_per_k: [[f64; 2]; 2],
    ) -> Result<Self, ConductionError> {
        Self::nodal_tensor(boundary_face, tangent_axes, [conductance_w_per_k; 3])
    }

    /// Reduce a uniform world-space bulk tensor [W/(m K)] under an explicitly
    /// declared zero-normal-heat-flux constraint, with positive thickness [m].
    ///
    /// In the orthonormal tangent/normal frame, `K = [[A,b],[b^T,c]]` gives
    /// `grad_n(T) = -b^T grad_s(T)/c` and `K_s = thickness*(A-b*b^T/c)`.
    /// The bulk tensor must be finite, symmetric and positive semidefinite,
    /// with strictly positive normal conductivity. Its Schur complement is
    /// checked without replacing negative eigenvalues by zero. Zero/rank-one
    /// reduced sheets remain admissible when the reduction is representable.
    ///
    /// Zero normal flux can require a NONZERO normal temperature gradient.
    /// The caller must justify neglecting the resulting temperature drop
    /// across this thin, bonded film. Neither perfect bonding nor thickness
    /// alone establishes that regime. A prescribed nonzero normal flux would
    /// add `thickness*b*q_n/c` to the tangential flux; this law does not include
    /// that coupling or solve a through-thickness field. It adds no storage,
    /// normal resistance or material-source authority.
    pub fn uniform_from_bulk_zero_normal_flux(
        boundary_face: usize,
        tangent_axes: [[f64; 3]; 2],
        bulk_w_per_m_k: [[f64; 3]; 3],
        thickness_m: f64,
    ) -> Result<Self, ConductionError> {
        check_axes(tangent_axes)?;
        crate::require_finite("sheet thickness", thickness_m)?;
        if thickness_m <= 0.0 {
            return Err(invalid("bulk sheet reduction needs positive thickness"));
        }
        let mut scale = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                let value = bulk_w_per_m_k[i][j];
                crate::require_finite("bulk sheet conductivity", value)?;
                if value != bulk_w_per_m_k[j][i] {
                    return Err(invalid("bulk sheet conductivity must be symmetric"));
                }
                scale = scale.max(value.abs());
            }
        }
        if scale == 0.0 {
            return Err(invalid(
                "bulk sheet reduction needs positive normal conductivity",
            ));
        }
        let scaled = bulk_w_per_m_k.map(|row| row.map(|v| v / scale));
        for (original, normalized) in bulk_w_per_m_k
            .into_iter()
            .flatten()
            .zip(scaled.into_iter().flatten())
        {
            if original != 0.0 && normalized == 0.0 {
                return Err(invalid("bulk sheet normalization underflowed"));
            }
        }
        let normal = cross(tangent_axes[0], tangent_axes[1]);
        let action = |v| scaled.map(|row| dot(row, v));
        let normal_action = action(normal);
        let normal_k = dot(normal, normal_action);
        if !normal_k.is_finite() || normal_k <= 0.0 {
            return Err(invalid(
                "bulk sheet reduction needs positive normal conductivity",
            ));
        }
        // Dividing by sqrt(c) BEFORE multiplying avoids an underflowed b_i*b_j
        // numerator when its quotient by c is still physically significant.
        let root_normal_k = fs_math::det::sqrt(normal_k);
        let coupling = tangent_axes.map(|axis| dot(axis, normal_action) / root_normal_k);
        for value in coupling {
            crate::require_finite("relaxed bulk sheet coupling", value)?;
        }
        let a = dot(tangent_axes[0], action(tangent_axes[0]));
        let b = dot(tangent_axes[0], action(tangent_axes[1]));
        let c = dot(tangent_axes[1], action(tangent_axes[1]));
        let off_diagonal = (-coupling[0]).mul_add(coupling[1], b);
        let reduced = [
            [(-coupling[0]).mul_add(coupling[0], a), off_diagonal],
            [off_diagonal, (-coupling[1]).mul_add(coupling[1], c)],
        ];
        check_tensor(reduced)?;
        let mut sheet = [[0.0; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                sheet[i][j] = sheet_product(reduced[i][j], scale, thickness_m)?;
            }
        }
        Self::uniform_tensor(boundary_face, tangent_axes, sheet)
    }

    /// Declare three P1 sheet tensors [W/K] in one fixed tangent frame, in the
    /// original face's vertex order. Interpolating PSD tensors preserves local
    /// nonnegative dissipation. The coefficient, rather than its eigenvalues
    /// or separately varying thickness and conductivity, is interpolated.
    pub fn nodal_tensor(
        boundary_face: usize,
        tangent_axes: [[f64; 3]; 2],
        nodal_w_per_k: [[[f64; 2]; 2]; 3],
    ) -> Result<Self, ConductionError> {
        check_axes(tangent_axes)?;
        for tensor in nodal_w_per_k {
            check_tensor(tensor)?;
        }
        Ok(Self {
            boundary_face,
            law: SheetLaw::Tensor {
                axes: tangent_axes,
                nodal: nodal_w_per_k,
            },
        })
    }

    /// Slot in the original mesh boundary table.
    #[must_use]
    pub const fn boundary_face(self) -> usize {
        self.boundary_face
    }

    /// Original isotropic P1 sheet conductances [W/K], or `None` for a tensor
    /// law. A directional tensor is never silently replaced by a scalar mean.
    #[must_use]
    pub const fn nodal_conductance_w_per_k(self) -> Option<[f64; 3]> {
        match self.law {
            SheetLaw::Isotropic(values) => Some(values),
            SheetLaw::Tensor { .. } => None,
        }
    }

    /// Declared world-space tangent axes, or `None` for an isotropic law.
    #[must_use]
    pub const fn tangent_axes(self) -> Option<[[f64; 3]; 2]> {
        match self.law {
            SheetLaw::Tensor { axes, .. } => Some(axes),
            SheetLaw::Isotropic(_) => None,
        }
    }

    /// Original P1 sheet tensors [W/K] in the declared tangent frame.
    #[must_use]
    pub const fn nodal_tensor_conductance_w_per_k(self) -> Option<[[[f64; 2]; 2]; 3]> {
        match self.law {
            SheetLaw::Tensor { nodal, .. } => Some(nodal),
            SheetLaw::Isotropic(_) => None,
        }
    }

    fn active(self) -> bool {
        match self.law {
            SheetLaw::Isotropic(values) => values.iter().any(|&value| value > 0.0),
            SheetLaw::Tensor { nodal, .. } => nodal
                .iter()
                .any(|tensor| tensor[0][0] > 0.0 || tensor[1][1] > 0.0),
        }
    }

    fn mean_action(self, gradient: [f64; 3]) -> [f64; 3] {
        match self.law {
            SheetLaw::Isotropic(values) => gradient.map(|g| mean(values) * g),
            SheetLaw::Tensor { axes, nodal } => {
                let a = mean(nodal.map(|k| k[0][0]));
                let b = mean(nodal.map(|k| k[0][1]));
                let c = mean(nodal.map(|k| k[1][1]));
                let g = axes.map(|axis| dot(axis, gradient));
                let q = [a.mul_add(g[0], b * g[1]), b.mul_add(g[0], c * g[1])];
                std::array::from_fn(|i| axes[0][i].mul_add(q[0], axes[1][i] * q[1]))
            }
        }
    }

    /// Exact face-average tangential flux -<K_s> grad_s(T_h), in W/m.
    /// A constant temperature has exactly zero gradient, independently of its
    /// Kelvin baseline. This is not the bulk flux density in W/m².
    pub fn mean_heat_flux_w_per_m(
        self,
        mesh: &ConductionMesh,
        temperature: &[f64],
    ) -> Result<[f64; 3], ConductionError> {
        if temperature.len() != mesh.vertex_count() {
            return Err(ConductionError::FieldLength {
                field: "surface temperature",
                expected: mesh.vertex_count(),
                found: temperature.len(),
            });
        }
        let (vertices, gradients, _) = self.geometry(mesh)?;
        let values = vertices.map(|v| temperature[v]);
        for value in values {
            crate::require_finite("surface temperature", value)?;
        }
        let gradient = std::array::from_fn(|axis| {
            (values[1] - values[0]).mul_add(
                gradients[1][axis],
                (values[2] - values[0]) * gradients[2][axis],
            )
        });
        let flux = self.mean_action(gradient).map(|q| -q);
        for value in flux {
            crate::require_finite("surface heat flux", value)?;
        }
        Ok(flux)
    }

    /// Exact P1 integral A*grad_s(N_a).<K_s> grad_s(N_b). Off-diagonal
    /// pairs determine the diagonal, retaining the constant-temperature null
    /// mode and symmetric internal transfer without adding a reservoir load.
    pub(crate) fn matrix(
        self,
        mesh: &ConductionMesh,
    ) -> Result<([usize; 3], [[f64; 3]; 3]), ConductionError> {
        let (vertices, gradients, area) = self.geometry(mesh)?;
        // Preserve the scalar arithmetic path; tensor faces add only fixed-size
        // stack work and the same nine COO entries/unchanged sparse graph.
        let isotropic_scale = match self.law {
            SheetLaw::Isotropic(values) => {
                let scale = area * mean(values);
                crate::require_finite("surface conductivity moment", scale)?;
                if self.active() && scale == 0.0 {
                    return Err(invalid(
                        "positive sheet conductance has an unrepresentable area moment",
                    ));
                }
                Some(scale)
            }
            SheetLaw::Tensor { .. } => None,
        };
        let actions = if isotropic_scale.is_none() {
            gradients.map(|gradient| self.mean_action(gradient))
        } else {
            [[0.0; 3]; 3]
        };
        let mut block = [[0.0; 3]; 3];
        for a in 0..3 {
            for b in (a + 1)..3 {
                let value = isotropic_scale.map_or_else(
                    || area * dot(gradients[a], actions[b]),
                    |scale| scale * dot(gradients[a], gradients[b]),
                );
                block[a][b] = value;
                block[b][a] = value;
                block[a][a] -= value;
                block[b][b] -= value;
            }
        }
        for value in block.into_iter().flatten() {
            crate::require_finite("surface conduction stiffness", value)?;
        }
        Ok((vertices, block))
    }

    fn geometry(
        self,
        mesh: &ConductionMesh,
    ) -> Result<([usize; 3], [[f64; 3]; 3], f64), ConductionError> {
        let face = mesh.boundary().get(self.boundary_face).ok_or_else(|| {
            invalid("surface conduction face is outside the original boundary table")
        })?;
        let vertices = face.vertices.map(|v| v as usize);
        let points = vertices.map(|v| mesh.positions()[v]);
        let first = std::array::from_fn(|i| points[1][i] - points[0][i]);
        let second = std::array::from_fn(|i| points[2][i] - points[0][i]);
        let inverse_twice_area = 0.5 / face.area;
        // The canonical vertex order, not the outward-normal sign, determines
        // this orientation. Gradients therefore remain correct on either side.
        let normal = cross(first, second).map(|v| v * inverse_twice_area);
        if let SheetLaw::Tensor { axes, .. } = self.law {
            if axes
                .iter()
                .any(|&axis| dot(axis, normal).abs() > FRAME_TOLERANCE)
            {
                return Err(invalid(
                    "sheet tensor axes must lie in the original face plane",
                ));
            }
        }
        let g1 = cross(second, normal).map(|v| v * inverse_twice_area);
        let g2 = cross(normal, first).map(|v| v * inverse_twice_area);
        let gradients = [std::array::from_fn(|i| -g1[i] - g2[i]), g1, g2];
        for value in gradients.into_iter().flatten() {
            crate::require_finite("surface basis gradient", value)?;
        }
        Ok((vertices, gradients, face.area))
    }
}

impl ThermalBoundary {
    /// Add perfectly bonded film conduction to this original face partition.
    /// It coexists with convection, prescribed flux/temperature, radiation and
    /// thermal contact, and remains in every shared spatial residual/tangent.
    ///
    /// Active faces must be unique across all additions; sum coincident layers'
    /// sheet conductances explicitly. Empty/zero additions preserve the original
    /// operator. Keep this boundary with the same mesh; no remapping is inferred.
    pub fn with_surface_heat_conduction(
        mut self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        surfaces: &[SurfaceHeatConduction],
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if self.face_count() != mesh.boundary().len() || surfaces.len() > mesh.boundary().len() {
            return Err(invalid(
                "surface conduction must use bounded original boundary-face slots",
            ));
        }
        let mut seen: BTreeSet<_> = self
            .surface_conduction
            .iter()
            .map(|s| s.boundary_face)
            .collect();
        let mut active = 0usize;
        for (index, &surface) in surfaces.iter().enumerate() {
            if index % crate::assemble::ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            if surface.boundary_face >= mesh.boundary().len() {
                return Err(invalid(
                    "surface conduction face is outside the original boundary table",
                ));
            }
            if !surface.active() {
                continue;
            }
            if !seen.insert(surface.boundary_face) {
                return Err(invalid(
                    "surface conduction face is declared more than once",
                ));
            }
            surface.matrix(mesh)?;
            active += 1;
        }
        self.surface_conduction
            .try_reserve(active)
            .map_err(|_| invalid("surface conduction allocation refused"))?;
        self.surface_conduction
            .extend(surfaces.iter().copied().filter(|s| s.active()));
        poll(cx, surfaces.len())?;
        Ok(self)
    }

    /// Active original face laws, with zero coatings omitted.
    #[must_use]
    pub fn surface_heat_conduction(&self) -> &[SurfaceHeatConduction] {
        &self.surface_conduction
    }
}

const FRAME_TOLERANCE: f64 = 1e-12;

fn check_axes(axes: [[f64; 3]; 2]) -> Result<(), ConductionError> {
    for axis in axes {
        for value in axis {
            crate::require_finite("sheet tangent axis", value)?;
        }
        if (dot(axis, axis) - 1.0).abs() > FRAME_TOLERANCE {
            return Err(invalid("sheet tangent axes must have unit length"));
        }
    }
    if dot(axes[0], axes[1]).abs() > FRAME_TOLERANCE {
        return Err(invalid("sheet tangent axes must be orthogonal"));
    }
    Ok(())
}

fn check_tensor(tensor: [[f64; 2]; 2]) -> Result<(), ConductionError> {
    for value in tensor.into_iter().flatten() {
        crate::require_finite("sheet conductance tensor", value)?;
    }
    if tensor[0][1] != tensor[1][0] {
        return Err(invalid("sheet conductance tensor must be symmetric"));
    }
    // The square-root product avoids overflow/underflow in a*c-b*b.
    if tensor[0][0] < 0.0
        || tensor[1][1] < 0.0
        || tensor[0][1].abs() > fs_math::det::sqrt(tensor[0][0]) * fs_math::det::sqrt(tensor[1][1])
    {
        return Err(invalid(
            "sheet conductance tensor must be positive semidefinite",
        ));
    }
    Ok(())
}

fn sheet_product(value: f64, bulk_scale: f64, thickness: f64) -> Result<f64, ConductionError> {
    if value == 0.0 {
        return Ok(value);
    }
    // Multiply the smallest and largest magnitudes first. A common k*t scale
    // can overflow even though each final reduced coefficient is representable.
    let mut factors = [value.abs(), bulk_scale, thickness];
    factors.sort_by(f64::total_cmp);
    let magnitude = (factors[0] * factors[2]) * factors[1];
    crate::require_finite("reduced sheet conductance", magnitude)?;
    if magnitude == 0.0 {
        return Err(invalid("nonzero reduced sheet conductance underflowed"));
    }
    Ok(magnitude.copysign(value))
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn mean(values: [f64; 3]) -> f64 {
    let scale = values.into_iter().map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 {
        0.0
    } else {
        scale * (values.into_iter().map(|v| v / scale).sum::<f64>() / 3.0)
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}
fn invalid(what: &str) -> ConductionError {
    ConductionError::Config {
        parameter: "surface heat conduction",
        what: what.into(),
    }
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "surface-conduction",
        at,
    })
}
