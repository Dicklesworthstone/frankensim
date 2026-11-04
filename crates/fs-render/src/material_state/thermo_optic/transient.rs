//! Native P1 thermal snapshots on the original conduction tetrahedra.
//!
//! The temperature-based theta and backward-Euler solvers publish one nodal
//! temperature per original mesh vertex. Their spatial reconstruction is P1:
//! T = sum_a lambda_a T_a, grad(T) = sum_a T_a grad(lambda_a), and the thermal
//! Hessian is zero inside each element. This adapter borrows that reconstruction
//! without smoothing nodal values, averaging elements, or solving another heat
//! problem. A nonlinear thermo-optic law still supplies the nonzero optical
//! Hessian n_TT grad(T) grad(T)^T through the existing material owner.
//!
//! An element is a bounded, one-sided field. Its native four-face oracle stops
//! curved propagation before a gradient jump or a thermal-contact discontinuity.
//! The optional [`mesh`] owner continues primal rays across original shared
//! vertex-ID faces in one continuous temperature/material field. Distinct
//! contact nodes still stop propagation. No global mesh locator, geometry
//! deformation, or interpolation between thermal times is inferred. The caller
//! selects one accepted immutable thermal state and its original mesh, retaining
//! its snapshot identity when
//! constructing `MaterialThermoOpticField`.
//!
//! A pointwise enthalpy reconstruction T(h(x)) is generally not P1 T. Its nodal
//! temperatures alone do not authorize this adapter, nor do cell averages.

use super::{
    TemperatureGradientField, TemperatureGradientSample, TemperatureJet,
    TemperatureSecondOrderField,
};
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::thermal::TemperatureField;
use fs_conduction::ConductionMesh;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

mod boundary;
pub use boundary::ConductionElementBoundary;
/// Bounded native primal-ray traversal across conforming temperature elements.
pub mod mesh;

/// Immutable original mesh and complete P1 nodal temperature field. Preparation
/// validates once in O(vertices), without allocating or copying the state.
#[derive(Clone, Copy)]
pub struct P1ConductionTemperatureSnapshot<'a> {
    mesh: &'a ConductionMesh,
    temperature: &'a [f64],
}

impl<'a> P1ConductionTemperatureSnapshot<'a> {
    /// Borrow the original temperature-based spatial reconstruction, normally
    /// `&accepted.temperature` from the native transient solver. The original
    /// mesh and state pairing is explicit; equal array lengths alone establish
    /// no physical correspondence. Every node must be finite and above 0 K.
    pub fn try_new(
        cx: &Cx<'_>,
        mesh: &'a ConductionMesh,
        temperature: &'a [f64],
    ) -> Result<Self, TransientTemperatureFieldError> {
        if temperature.len() != mesh.vertex_count() {
            return Err(TransientTemperatureFieldError::FieldLength {
                expected: mesh.vertex_count(),
                found: temperature.len(),
            });
        }
        for (vertex, &value) in temperature.iter().enumerate() {
            if vertex % 512 == 0 {
                cx.checkpoint()
                    .map_err(|_| TransientTemperatureFieldError::Cancelled)?;
            }
            if !value.is_finite() || value <= 0.0 {
                return Err(TransientTemperatureFieldError::InvalidTemperature { vertex });
            }
        }
        Ok(Self { mesh, temperature })
    }

    pub const fn mesh(&self) -> &'a ConductionMesh {
        self.mesh
    }

    pub const fn nodal_temperature_k(&self) -> &'a [f64] {
        self.temperature
    }

    /// Prepare one original tetrahedron in constant work and storage. The
    /// original fs-feec barycentric gradients are reused; no matrix is inverted
    /// here. Subsequent temperature, derivative and geometry queries allocate
    /// nothing. Element ordering is the conduction mesh's original ordering.
    pub fn element(
        &self,
        element: usize,
    ) -> Result<P1ConductionElement<'a>, TransientTemperatureFieldError> {
        let vertices = self
            .mesh
            .complex()
            .tets
            .get(element)
            .ok_or(TransientTemperatureFieldError::InvalidElement { element })?;
        let positions = vertices.map(|v| {
            let p = self.mesh.positions()[v as usize];
            Point3::new(p[0], p[1], p[2])
        });
        let values = vertices.map(|v| self.temperature[v as usize]);
        let gradients = self.mesh.geometry().grads[element].map(|g| Vec3::new(g[0], g[1], g[2]));
        // Temperature differences preserve an exactly uniform field, avoiding
        // cancellation between large absolute-temperature gradient terms.
        let differences = [
            values[1] - values[0],
            values[2] - values[0],
            values[3] - values[0],
        ];
        let component = |g: [f64; 3]| {
            differences[0].mul_add(g[0], differences[1].mul_add(g[1], differences[2] * g[2]))
        };
        let gradient = Vec3::new(
            component([gradients[1].x, gradients[2].x, gradients[3].x]),
            component([gradients[1].y, gradients[2].y, gradients[3].y]),
            component([gradients[1].z, gradients[2].z, gradients[3].z]),
        );
        if ![gradient.x, gradient.y, gradient.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(TransientTemperatureFieldError::NumericalRange);
        }
        Ok(P1ConductionElement {
            snapshot: *self,
            element,
            positions,
            values,
            gradients,
            gradient,
        })
    }
}

/// Exact original P1 thermal reconstruction on one closed tetrahedron. Its
/// derivatives are one-sided on a face; neighboring gradients are not averaged.
#[derive(Clone, Copy)]
pub struct P1ConductionElement<'a> {
    snapshot: P1ConductionTemperatureSnapshot<'a>,
    element: usize,
    positions: [Point3; 4],
    values: [f64; 4],
    gradients: [Vec3; 4],
    gradient: Vec3,
}

impl<'a> P1ConductionElement<'a> {
    pub const fn snapshot(&self) -> P1ConductionTemperatureSnapshot<'a> {
        self.snapshot
    }

    pub const fn element_index(&self) -> usize {
        self.element
    }

    pub const fn vertex_positions_m(&self) -> &[Point3; 4] {
        &self.positions
    }

    pub const fn vertex_temperature_k(&self) -> &[f64; 4] {
        &self.values
    }

    /// Actual original face coordinates, outward ordered. Local face i is
    /// opposite local vertex i, independently of the tet's signed orientation.
    pub fn faces_m(&self) -> [[[f64; 3]; 3]; 4] {
        let mut indices = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];
        if self.snapshot.mesh.geometry().vol_signed[self.element] < 0.0 {
            for face in &mut indices {
                face.swap(1, 2);
            }
        }
        indices.map(|face| {
            face.map(|v| {
                let p = self.positions[v];
                [p.x, p.y, p.z]
            })
        })
    }

    /// The source's barycentric coordinates, with no domain extension or
    /// clamping. Stored vertices retain their exact nodal values. Floating-point
    /// uncertainty at an arbitrary boundary point can refuse; native event
    /// propagation approaches from the element interior to its declared metric
    /// tolerance instead of teleporting onto a face.
    pub fn barycentric(
        &self,
        position_m: Point3,
    ) -> Result<[f64; 4], TransientTemperatureFieldError> {
        if ![position_m.x, position_m.y, position_m.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(TransientTemperatureFieldError::OutsideElement);
        }
        if let Some(vertex) = self.positions.iter().position(|p| *p == position_m) {
            let mut value = [0.0; 4];
            value[vertex] = 1.0;
            return Ok(value);
        }
        let origin = self.positions[0];
        let displacement = Vec3::new(
            position_m.x - origin.x,
            position_m.y - origin.y,
            position_m.z - origin.z,
        );
        let mut lambda = [0.0; 4];
        for (entry, gradient) in lambda.iter_mut().zip(self.gradients).skip(1) {
            *entry = gradient.dot(displacement);
        }
        lambda[0] = ((1.0 - lambda[1]) - lambda[2]) - lambda[3];
        if lambda
            .iter()
            .any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
        {
            return Err(TransientTemperatureFieldError::OutsideElement);
        }
        Ok(lambda)
    }

    pub fn sample_temperature_jet(
        &self,
        position_m: Point3,
    ) -> Result<TemperatureJet, TransientTemperatureFieldError> {
        let lambda = self.barycentric(position_m)?;
        let vertex = lambda.iter().enumerate().find_map(|(index, &value)| {
            (value == 1.0
                && lambda
                    .iter()
                    .enumerate()
                    .all(|(i, &v)| i == index || v == 0.0))
            .then_some(index)
        });
        let temperature = if let Some(vertex) = vertex {
            self.values[vertex]
        } else {
            lambda[1].mul_add(
                self.values[1] - self.values[0],
                lambda[2].mul_add(
                    self.values[2] - self.values[0],
                    lambda[3].mul_add(self.values[3] - self.values[0], self.values[0]),
                ),
            )
        };
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(TransientTemperatureFieldError::NumericalRange);
        }
        Ok(TemperatureJet {
            sample: TemperatureGradientSample {
                temperature_k: temperature,
                gradient_k_per_m: self.gradient,
            },
            diagonal_k_per_m2: [0.0; 3],
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}

impl TemperatureField for P1ConductionElement<'_> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_temperature_jet(p)
            .map(|v| v.sample.temperature_k)
            .map_err(TransientTemperatureFieldError::into_gradient_error)
    }
}
impl TemperatureGradientField for P1ConductionElement<'_> {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        self.sample_temperature_jet(p)
            .map(|v| v.sample)
            .map_err(TransientTemperatureFieldError::into_gradient_error)
    }
}
impl TemperatureSecondOrderField for P1ConductionElement<'_> {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        self.sample_temperature_jet(p)
            .map_err(TransientTemperatureFieldError::into_gradient_error)
    }
}

/// Snapshot admission and native P1 evaluation failures.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransientTemperatureFieldError {
    FieldLength { expected: usize, found: usize },
    InvalidTemperature { vertex: usize },
    InvalidElement { element: usize },
    OutsideElement,
    NumericalRange,
    Cancelled,
}
impl TransientTemperatureFieldError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::OutsideElement => GradientIndexError::OutsideFieldDomain,
            Self::NumericalRange => GradientIndexError::NumericalRange,
            Self::Cancelled => GradientIndexError::Cancelled,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl core::fmt::Display for TransientTemperatureFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient P1 temperature field: {self:?}")
    }
}
impl core::error::Error for TransientTemperatureFieldError {}

#[cfg(test)]
mod tests;
