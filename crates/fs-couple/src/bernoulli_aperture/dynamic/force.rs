//! Physical mechanical loads on the original plate, projected by virtual work.
//! Positive force acts in positive plate displacement. A patch force is TOTAL N
//! spread uniformly over its selected triangles, not a pressure or fitted gain.
use super::{AcousticRealizeError, invalid};
use crate::bernoulli_aperture::plate::PlateApertureReduction;

/// Explicit application site on the retained mesh; no nearest-node inference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PlateForceFootprint {
    /// Concentrated transverse force at an original node (a point-load idealization).
    Node(usize),
    /// Uniform transverse traction over unique original triangles; force is total N.
    Patch(Vec<usize>),
}

/// The SAME signed coefficient maps force to the opening coordinate and maps
/// opening velocity to physical actuator velocity: Q=B F, v_port=B v_opening.
/// This is a fixed linear projection, not an additional dynamic/contact body.
#[derive(Clone, Debug)]
pub struct PlateForcePort {
    footprint: PlateForceFootprint,
    coefficient: f64,
    area_m2: Option<f64>,
}
impl PlateForcePort {
    /// Original node or triangle addresses, in declared order.
    #[must_use]
    pub const fn footprint(&self) -> &PlateForceFootprint { &self.footprint }
    /// Dimensionless signed force/velocity projection. Zero is a stationary site
    /// in this retained basis, not a missing-data error.
    #[must_use]
    pub const fn coefficient(&self) -> f64 { self.coefficient }
    /// Actual patch area [m²], absent for an ideal point force.
    #[must_use]
    pub const fn area_m2(&self) -> Option<f64> { self.area_m2 }
    /// Opening-coordinate force [N]. No pressure area or PCM scaling enters it.
    ///
    /// # Errors
    /// Nonfinite input or overflow in the projected force.
    pub fn generalized_force_n(&self, force_n: f64) -> Result<f64, AcousticRealizeError> {
        let value = self.coefficient * force_n;
        if !force_n.is_finite() || !value.is_finite() {
            return Err(invalid("plate actuator force must project to finite newtons"));
        }
        Ok(value)
    }
    /// Physical actuator velocity [m/s]; patch ports return the area mean.
    ///
    /// # Errors
    /// Nonfinite input or overflow in the projected velocity.
    pub fn velocity_m_s(&self, opening_velocity_m_s: f64) -> Result<f64, AcousticRealizeError> {
        let value = self.coefficient * opening_velocity_m_s;
        if !opening_velocity_m_s.is_finite() || !value.is_finite() {
            return Err(invalid("plate actuator velocity must project to finite metres per second"));
        }
        Ok(value)
    }
}
impl PlateApertureReduction {
    /// Project a physical load on THIS retained specimen. Patches use the same
    /// linear displacement trace as uniform-pressure work on the source plate.
    /// Supports, thickness/material changes and the actual mode shape therefore
    /// affect the actuator coupling; no effective lever arm is prescribed.
    ///
    /// # Errors
    /// Missing node, empty/duplicate/absent triangles, or nonfinite projection.
    pub fn force_port(&self, footprint: PlateForceFootprint) -> Result<PlateForcePort, AcousticRealizeError> {
        let (coefficient, area_m2) = match &footprint {
            PlateForceFootprint::Node(node) => {
                let shape = self.shape_per_opening().get(*node)
                    .ok_or_else(|| invalid("plate actuator node is outside the original mesh"))?;
                (shape[0], None)
            }
            PlateForceFootprint::Patch(triangles) => {
                let mesh = &self.chart().mesh;
                if triangles.is_empty() || triangles.len() > mesh.tris.len() {
                    return Err(invalid("plate actuator patch needs a nonempty bounded triangle set"));
                }
                let mut seen = std::collections::BTreeSet::new();
                let (mut area, mut integral) = (0.0, 0.0);
                for &index in triangles {
                    let nodes = mesh.tris.get(index)
                        .ok_or_else(|| invalid("plate actuator triangle is outside the original mesh"))?;
                    if !seen.insert(index) { return Err(invalid("plate actuator patch repeats a triangle")); }
                    let [a, b, c] = nodes.map(|i| mesh.nodes[i]);
                    let da = 0.5 * ((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0));
                    let mean = nodes.iter().map(|&i| self.shape_per_opening()[i][0]).sum::<f64>() / 3.0;
                    area += da;
                    integral += da * mean;
                }
                if !area.is_finite() || area <= 0.0 { return Err(invalid("plate actuator patch has no finite positive area")); }
                (integral / area, Some(area))
            }
        };
        if !coefficient.is_finite() { return Err(invalid("plate actuator projection overflowed")); }
        Ok(PlateForcePort { footprint, coefficient, area_m2 })
    }
}
