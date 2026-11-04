//! Original P1 pressure and temperature on one shared native geometry.
//!
//! Pressure is a separately supplied, accepted nodal field in absolute Pa.
//! The conduction mesh is borrowed only for its original topology, coordinates
//! and P1 gradients; no pressure solve or temperature interpretation is inferred.
//! The temperature snapshot remains the original thermal owner's reconstruction.
//! Cell averages, DG traces and acoustic pressure perturbations do not implicitly
//! become absolute nodal pressure through this adapter.
//!
//! Each element uses one native barycentric admission for both fields. Their
//! in-element Hessians are zero; the original IF97 density jet and IAPWS optical
//! chain rule still supply the complete nonlinear optical Hessian, including
//! mixed pressure/temperature terms. Source preparation and sampling allocate
//! nothing, and the original borrowed buffers and accepted identities remain
//! available to the caller.

use crate::gradient_index::GradientIndexError;
use crate::material_state::thermo_optic::transient::{
    P1ConductionElement, P1ConductionTemperatureSnapshot, TransientTemperatureFieldError,
};
use crate::material_state::water_refraction::pressure_temperature::{
    PressureTemperatureGradientField, PressureTemperatureGradientSample, PressureTemperatureJet,
    PressureTemperatureSecondOrderField, PressureTemperatureWaterError, PressureTemperatureWaterField,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_conduction::ConductionMesh;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_qty::Pressure;
use fs_thermochem::water::if97::If97Region;

/// A borrowed SI pressure view; neither representation copies or converts a buffer.
/// Both variants mean absolute thermodynamic pressure, never gauge pressure.
#[derive(Clone, Copy)]
pub enum NodalPressureValues<'a> {
    Quantities(&'a [Pressure]),
    Pascals(&'a [f64]),
}
impl NodalPressureValues<'_> {
    pub const fn len(self) -> usize {
        match self {
            Self::Quantities(values) => values.len(),
            Self::Pascals(values) => values.len(),
        }
    }
    pub const fn is_empty(self) -> bool {
        self.len() == 0
    }
    pub fn get(self, vertex: usize) -> Option<Pressure> {
        match self {
            Self::Quantities(values) => values.get(vertex).copied(),
            Self::Pascals(values) => values.get(vertex).map(|&value| Pressure::new(value)),
        }
    }
    fn value_at(self, vertex: usize) -> f64 {
        match self {
            Self::Quantities(values) => values[vertex].value(),
            Self::Pascals(values) => values[vertex],
        }
    }
}

/// One immutable original nodal pressure field with explicit native geometry.
/// The accepted identity must name that field and its original P1 reconstruction.
#[derive(Clone, Copy)]
pub struct P1PressureSnapshot<'a> {
    mesh: &'a ConductionMesh,
    pressure: NodalPressureValues<'a>,
    snapshot_identity: ContentHash,
}
impl<'a> P1PressureSnapshot<'a> {
    /// Bind an accepted pressure field in the original mesh's vertex ordering.
    /// Equal lengths alone do not establish physical correspondence. This is
    /// a caller-declared original nodal reconstruction, not a resampling step.
    pub fn try_new(
        cx: &Cx<'_>,
        mesh: &'a ConductionMesh,
        pressure: NodalPressureValues<'a>,
        snapshot_identity: ContentHash,
    ) -> Result<Self, PressureTemperatureSnapshotError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(PressureTemperatureSnapshotError::MissingPressureIdentity);
        }
        if pressure.len() != mesh.vertex_count() {
            return Err(PressureTemperatureSnapshotError::PressureLength {
                expected: mesh.vertex_count(), found: pressure.len(),
            });
        }
        for vertex in 0..pressure.len() {
            if vertex % 512 == 0 {
                cx.checkpoint().map_err(|_| PressureTemperatureSnapshotError::Cancelled)?;
            }
            let value = pressure.value_at(vertex);
            if !value.is_finite() || value <= 0.0 {
                return Err(PressureTemperatureSnapshotError::InvalidPressure { vertex });
            }
        }
        Ok(Self { mesh, pressure, snapshot_identity })
    }
    pub const fn mesh(self) -> &'a ConductionMesh {
        self.mesh
    }
    pub const fn nodal_pressure_values(self) -> NodalPressureValues<'a> {
        self.pressure
    }
    pub const fn snapshot_identity(self) -> ContentHash {
        self.snapshot_identity
    }

    fn polynomial(self, element: usize) -> Result<PressurePolynomial, PressureTemperatureSnapshotError> {
        let vertices = self.mesh.complex().tets.get(element)
            .ok_or(PressureTemperatureSnapshotError::InvalidElement { element })?;
        let values_pa = vertices.map(|vertex| self.pressure.value_at(vertex as usize));
        let gradients = self.mesh.geometry().grads[element];
        let differences = [
            values_pa[1] - values_pa[0],
            values_pa[2] - values_pa[0],
            values_pa[3] - values_pa[0],
        ];
        // Differences retain exactly uniform pressure without cancellation
        // between large absolute-pressure terms.
        let component = |axis: usize| {
            differences[0].mul_add(gradients[1][axis],
                differences[1].mul_add(gradients[2][axis], differences[2] * gradients[3][axis]))
        };
        let gradient_pa_per_m = Vec3::new(component(0), component(1), component(2));
        if ![gradient_pa_per_m.x, gradient_pa_per_m.y, gradient_pa_per_m.z]
            .into_iter().all(f64::is_finite)
        {
            return Err(PressureTemperatureSnapshotError::NumericalRange);
        }
        Ok(PressurePolynomial { values_pa, gradient_pa_per_m })
    }
}

/// The original temperature and pressure fields on exactly the same mesh.
/// The combined identity includes BOTH accepted sources and the P1 reconstruction.
#[derive(Clone, Copy)]
pub struct P1PressureTemperatureSnapshot<'a> {
    temperature: P1ConductionTemperatureSnapshot<'a>,
    temperature_identity: ContentHash,
    pressure: P1PressureSnapshot<'a>,
    identity: ContentHash,
}
impl<'a> P1PressureTemperatureSnapshot<'a> {
    pub fn try_new(
        temperature: P1ConductionTemperatureSnapshot<'a>,
        temperature_identity: ContentHash,
        pressure: P1PressureSnapshot<'a>,
    ) -> Result<Self, PressureTemperatureSnapshotError> {
        if temperature_identity == ContentHash([0; 32]) {
            return Err(PressureTemperatureSnapshotError::MissingTemperatureIdentity);
        }
        if !core::ptr::eq(temperature.mesh(), pressure.mesh()) {
            return Err(PressureTemperatureSnapshotError::ForeignGeometry);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.p1-pressure-temperature.v1");
        hash.update(temperature_identity.as_bytes());
        hash.update(pressure.snapshot_identity().as_bytes());
        Ok(Self { temperature, temperature_identity, pressure, identity: hash.finalize() })
    }
    pub const fn temperature_snapshot(self) -> P1ConductionTemperatureSnapshot<'a> {
        self.temperature
    }
    pub const fn temperature_snapshot_identity(self) -> ContentHash {
        self.temperature_identity
    }
    pub const fn pressure_snapshot(self) -> P1PressureSnapshot<'a> {
        self.pressure
    }
    pub const fn content_identity(self) -> ContentHash {
        self.identity
    }
    pub fn element(self, element: usize) -> Result<P1PressureTemperatureElement<'a>, PressureTemperatureSnapshotError> {
        Ok(P1PressureTemperatureElement {
            snapshot: self,
            temperature: self.temperature.element(element)?,
            pressure: self.pressure.polynomial(element)?,
        })
    }
}

/// One original joint P1 element with native geometric and source admission.
#[derive(Clone, Copy)]
pub struct P1PressureTemperatureElement<'a> {
    snapshot: P1PressureTemperatureSnapshot<'a>,
    temperature: P1ConductionElement<'a>,
    pressure: PressurePolynomial,
}
impl<'a> P1PressureTemperatureElement<'a> {
    pub const fn snapshot(self) -> P1PressureTemperatureSnapshot<'a> {
        self.snapshot
    }
    /// Original geometry/boundary owner; it remains a temperature field, not a
    /// pressure surrogate. Pressure is retained separately by this joint element.
    pub const fn temperature_element(self) -> P1ConductionElement<'a> {
        self.temperature
    }
    pub fn sample_jet(&self, position_m: Point3) -> Result<PressureTemperatureJet, PressureTemperatureSnapshotError> {
        let lambda = self.temperature.barycentric(position_m)?;
        joint_jet(&self.temperature, self.pressure, lambda)
    }
    /// Bind the existing IF97 owner using both original source identities.
    pub fn bind_water(
        &self,
        region: If97Region,
    ) -> Result<PressureTemperatureWaterField<'_, Self>, PressureTemperatureWaterError> {
        PressureTemperatureWaterField::try_new(self, region, self.snapshot.content_identity())
    }
}
impl PressureTemperatureGradientField for P1PressureTemperatureElement<'_> {
    fn pressure_temperature_gradient(
        &self,
        position_m: Point3,
    ) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        self.sample_jet(position_m).map(|jet| jet.sample())
            .map_err(PressureTemperatureSnapshotError::into_gradient_error)
    }
}
impl PressureTemperatureSecondOrderField for P1PressureTemperatureElement<'_> {
    fn pressure_temperature_jet(&self, position_m: Point3) -> Result<PressureTemperatureJet, GradientIndexError> {
        self.sample_jet(position_m).map_err(PressureTemperatureSnapshotError::into_gradient_error)
    }
}

#[derive(Clone, Copy)]
struct PressurePolynomial {
    values_pa: [f64; 4],
    gradient_pa_per_m: Vec3,
}
impl PressurePolynomial {
    fn sample(self, lambda: [f64; 4]) -> Result<f64, PressureTemperatureSnapshotError> {
        let vertex = lambda.iter().enumerate().find_map(|(index, &value)| {
            (value == 1.0 && lambda.iter().enumerate().all(|(i, &v)| i == index || v == 0.0))
                .then_some(index)
        });
        let pressure_pa = if let Some(vertex) = vertex {
            self.values_pa[vertex]
        } else {
            lambda[1].mul_add(
                self.values_pa[1] - self.values_pa[0],
                lambda[2].mul_add(
                    self.values_pa[2] - self.values_pa[0],
                    lambda[3].mul_add(self.values_pa[3] - self.values_pa[0], self.values_pa[0]),
                ),
            )
        };
        if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
            return Err(PressureTemperatureSnapshotError::NumericalRange);
        }
        Ok(pressure_pa)
    }
}
fn joint_jet(
    temperature: &P1ConductionElement<'_>,
    pressure: PressurePolynomial,
    lambda: [f64; 4],
) -> Result<PressureTemperatureJet, PressureTemperatureSnapshotError> {
    Ok(PressureTemperatureJet {
        temperature: temperature.temperature_jet_at_barycentric(lambda)?,
        pressure_pa: pressure.sample(lambda)?,
        pressure_gradient_pa_per_m: pressure.gradient_pa_per_m,
        pressure_diagonal_pa_per_m2: [0.0; 3],
        pressure_off_diagonal_pa_per_m2: [0.0; 3],
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PressureTemperatureSnapshotError {
    PressureLength { expected: usize, found: usize },
    InvalidPressure { vertex: usize },
    InvalidElement { element: usize },
    MissingPressureIdentity,
    MissingTemperatureIdentity,
    ForeignGeometry,
    Temperature(TransientTemperatureFieldError),
    NumericalRange,
    Cancelled,
}
impl PressureTemperatureSnapshotError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Temperature(error) => error.into_gradient_error(),
            Self::NumericalRange => GradientIndexError::NumericalRange,
            Self::Cancelled => GradientIndexError::Cancelled,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl From<TransientTemperatureFieldError> for PressureTemperatureSnapshotError {
    fn from(error: TransientTemperatureFieldError) -> Self {
        Self::Temperature(error)
    }
}
impl core::fmt::Display for PressureTemperatureSnapshotError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original P1 pressure/temperature field: {self:?}")
    }
}
impl core::error::Error for PressureTemperatureSnapshotError {}

#[cfg(test)]
mod tests;
