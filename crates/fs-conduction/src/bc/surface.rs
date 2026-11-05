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

/// Isotropic, temperature-independent sheet conductance on one original face.
/// Coefficients follow that face's vertex order and are interpolated with P1
/// basis functions. No material provenance or tensor isotropy is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceHeatConduction {
    boundary_face: usize,
    nodal_w_per_k: [f64; 3],
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
            nodal_w_per_k,
        })
    }

    /// Slot in the original mesh boundary table.
    #[must_use]
    pub const fn boundary_face(self) -> usize {
        self.boundary_face
    }

    /// Original P1 sheet conductances [W/K].
    #[must_use]
    pub const fn nodal_conductance_w_per_k(self) -> [f64; 3] {
        self.nodal_w_per_k
    }

    fn active(self) -> bool {
        self.nodal_w_per_k.iter().any(|&value| value > 0.0)
    }

    fn mean(self) -> f64 {
        let scale = self.nodal_w_per_k.into_iter().fold(0.0, f64::max);
        if scale == 0.0 {
            0.0
        } else {
            scale
                * (self
                    .nodal_w_per_k
                    .into_iter()
                    .map(|v| v / scale)
                    .sum::<f64>()
                    / 3.0)
        }
    }

    /// Exact face-average tangential flux -<k*t> grad_s(T_h), in W/m.
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
        let mean = self.mean();
        let flux = std::array::from_fn(|axis| {
            -mean
                * (values[1] - values[0]).mul_add(
                    gradients[1][axis],
                    (values[2] - values[0]) * gradients[2][axis],
                )
        });
        for value in flux {
            crate::require_finite("surface heat flux", value)?;
        }
        Ok(flux)
    }

    /// Exact P1 integral A*<k*t>*grad_s(N_a).grad_s(N_b). Off-diagonal
    /// pairs determine the diagonal, retaining the constant-temperature null
    /// mode and symmetric internal transfer without adding a reservoir load.
    pub(crate) fn matrix(
        self,
        mesh: &ConductionMesh,
    ) -> Result<([usize; 3], [[f64; 3]; 3]), ConductionError> {
        let (vertices, gradients, area) = self.geometry(mesh)?;
        let scale = area * self.mean();
        crate::require_finite("surface conductivity moment", scale)?;
        if self.active() && scale == 0.0 {
            return Err(invalid(
                "positive sheet conductance has an unrepresentable area moment",
            ));
        }
        let mut block = [[0.0; 3]; 3];
        for a in 0..3 {
            for b in (a + 1)..3 {
                let value = scale
                    * gradients[a][0].mul_add(
                        gradients[b][0],
                        gradients[a][1].mul_add(gradients[b][1], gradients[a][2] * gradients[b][2]),
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
