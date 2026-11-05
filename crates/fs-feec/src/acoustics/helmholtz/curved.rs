//! Original curved geometry in the existing complete acoustic pressure space.

use super::geometry::{GeometryFacePoint, GeometryPoint, QuadraticGeometry};
use super::{AcousticFeError, HelmholtzFe, Medium, QuadraticGeometryControls};
use crate::highorder::simplex::duffy_quadrature;
use fs_exec::Cx;
use fs_rep_mesh::TetComplex;
use fs_sparse::Coo;

/// Largest number of Gauss points on each curved-geometry Duffy axis.
pub const MAX_GEOMETRY_QUADRATURE_ORDER: usize = 16;

impl<'c> HelmholtzFe<'c> {
    /// Assemble the original H1 pressure space on a conforming quadratic map.
    ///
    /// `edge_midpoints_m` follows `complex.edges` exactly. Each supplied point
    /// is the physical position of the reference edge midpoint, not a Bezier
    /// control point. Shared edges therefore share geometry on every incident
    /// face. All ten geometric nodes are retained independently of pressure
    /// order. The existing interval/Bernstein query must prove a nonsingular
    /// orientation over each complete closed element before any form is built.
    ///
    /// `quadrature_order` is an independent convergence axis in
    /// `order+2..=MAX_GEOMETRY_QUADRATURE_ORDER`. It is also the minimum used by
    /// every default curved face form. Curved stiffness and surface measures
    /// are generally nonpolynomial; this count is not an accuracy certificate.
    /// Local positive orientation does not prove global mesh injectivity.
    /// Cancellation is checked during proof and every volume quadrature point.
    #[allow(clippy::too_many_arguments)]
    pub fn new_quadratic(
        complex: &'c TetComplex,
        positions: &[[f64; 3]],
        edge_midpoints_m: &[[f64; 3]],
        order: usize,
        quadrature_order: usize,
        medium: Medium,
        controls: QuadraticGeometryControls,
        cx: &Cx<'_>,
    ) -> Result<Self, AcousticFeError> {
        cx.checkpoint().map_err(|_| AcousticFeError::Cancelled)?;
        if quadrature_order < order.saturating_add(2)
            || quadrature_order > MAX_GEOMETRY_QUADRATURE_ORDER
        {
            return Err(AcousticFeError::InvalidInput {
                what: "curved acoustic quadrature order outside admitted range",
            });
        }
        let mut model = Self::initialize(complex, positions, order, medium)?;
        model.geometry = Some(QuadraticGeometry::new(
            complex,
            positions,
            edge_midpoints_m,
            controls,
            cx,
        )?);
        model.geometry_quadrature_order = quadrature_order;
        model.assemble_curved(cx)?;
        Ok(model)
    }

    /// Whether the owner retains a quadratic geometric map.
    #[must_use]
    pub const fn is_curved(&self) -> bool {
        self.geometry.is_some()
    }

    /// Declared volume quadrature count; a minimum for curved surface forms.
    #[must_use]
    pub const fn geometry_quadrature_order(&self) -> usize {
        self.geometry_quadrature_order
    }

    /// Physical map, volume measure and exact mapped barycentric gradients.
    /// `barycentric` follows the original tetrahedron's stored corner order.
    /// Finite affine coordinates summing to one may be sampled outside the
    /// simplex, but whole-element orientation evidence applies only inside it.
    /// Physical probes impose their own closed-simplex location admission.
    pub fn geometry_point(
        &self,
        element: usize,
        barycentric: [f64; 4],
    ) -> Result<GeometryPoint, AcousticFeError> {
        match &self.geometry {
            Some(geometry) => {
                geometry.point(self.space.complex, &self.positions, element, barycentric)
            }
            None => super::geometry::affine_point(
                self.space.complex,
                &self.positions,
                element,
                barycentric,
                &self.affine_geometry,
            ),
        }
    }

    /// Pointwise surface measure and normal outward from one original cell.
    ///
    /// `omitted` is the local corner opposite the face and its barycentric
    /// coordinate must be zero. The returned surface Jacobian multiplies the
    /// reference triangle quadrature (whose weights sum to one half).
    pub fn geometry_face_point(
        &self,
        element: usize,
        omitted: usize,
        barycentric: [f64; 4],
    ) -> Result<GeometryFacePoint, AcousticFeError> {
        if omitted >= 4 || barycentric[omitted] != 0.0 {
            return Err(AcousticFeError::InvalidInput {
                what: "surface sample must lie on the selected original face",
            });
        }
        self.geometry_point(element, barycentric)?.face(omitted)
    }

    /// The affine constructor retains its original quadrature arithmetic.
    pub(super) fn face_quadrature_order(&self, minimum: usize) -> usize {
        if self.is_curved() {
            minimum.max(self.geometry_quadrature_order)
        } else {
            minimum
        }
    }

    fn assemble_curved(&mut self, cx: &Cx<'_>) -> Result<(), AcousticFeError> {
        let quad = duffy_quadrature(self.geometry_quadrature_order);
        let mut stiffness = Coo::new(self.ndof(), self.ndof());
        let mut mass = Coo::new(self.ndof(), self.ndof());
        for element in 0..self.tetrahedra().len() {
            cx.checkpoint().map_err(|_| AcousticFeError::Cancelled)?;
            let dofs = self.space.element_dofs(element);
            let n = dofs.len();
            let mut k = vec![0.0; n * n];
            let mut m = vec![0.0; n * n];
            let mut values = vec![0.0; n];
            let mut gradients = vec![[0.0; 3]; n];
            for &(barycentric, weight) in &quad {
                cx.checkpoint().map_err(|_| AcousticFeError::Cancelled)?;
                let point = self.geometry_point(element, barycentric)?;
                let dv = weight * point.determinant_m3;
                if !dv.is_finite() || dv <= 0.0 {
                    return Err(AcousticFeError::DegenerateTet { tet: element });
                }
                for (i, (_, shape)) in dofs.iter().enumerate() {
                    values[i] = shape.eval(barycentric, self.order);
                    let derivatives = shape.d_lambda(barycentric, self.order);
                    gradients[i] = std::array::from_fn(|axis| {
                        (0..4)
                            .map(|a| derivatives[a] * point.barycentric_gradients[a][axis])
                            .sum()
                    });
                }
                for i in 0..n {
                    for j in i..n {
                        let product: f64 = (0..3)
                            .map(|axis| gradients[i][axis] * gradients[j][axis])
                            .sum();
                        k[i * n + j] += dv * product;
                        m[i * n + j] += dv * values[i] * values[j];
                    }
                }
            }
            for (i, &(row, _)) in dofs.iter().enumerate() {
                for (j, &(col, _)) in dofs.iter().enumerate().skip(i) {
                    let (ki, mi) = (k[i * n + j], m[i * n + j]);
                    if !ki.is_finite() || !mi.is_finite() {
                        return Err(AcousticFeError::InvalidInput {
                            what: "curved acoustic assembly overflow",
                        });
                    }
                    stiffness.push(row, col, ki);
                    mass.push(row, col, mi);
                    if i != j {
                        stiffness.push(col, row, ki);
                        mass.push(col, row, mi);
                    }
                }
            }
        }
        cx.checkpoint().map_err(|_| AcousticFeError::Cancelled)?;
        let stiffness = stiffness.assemble();
        let mass = mass.assemble();
        for row in 0..self.ndof() {
            cx.checkpoint().map_err(|_| AcousticFeError::Cancelled)?;
            if stiffness
                .row(row)
                .1
                .iter()
                .chain(mass.row(row).1)
                .any(|v| !v.is_finite())
            {
                return Err(AcousticFeError::InvalidInput {
                    what: "curved acoustic global assembly overflow",
                });
            }
        }
        self.stiffness = stiffness;
        self.mass = mass;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
