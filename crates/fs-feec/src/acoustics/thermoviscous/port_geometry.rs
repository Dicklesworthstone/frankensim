//! Geometry-only access to the original work-conjugate pressure ports.

use super::ThermoviscousFe;

/// One original planar boundary triangle, without a new fluid discretization.
/// Vertex ordering is canonical and need not encode its outward orientation;
/// the retained normal carries that orientation explicitly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousPortTriangle {
    pub vertices_m: [[f64; 3]; 3],
    pub outward_normal: [f64; 3],
    pub area_m2: f64,
}

impl ThermoviscousFe<'_> {
    /// Original triangles of a unit, planar, exclusively owned pressure port.
    /// Admission is exactly `plane_wave_port_area_m2`: scaled, signed, curved
    /// or overlapping traction modes return None. No area normalization or
    /// replacement boundary/volume operator is constructed.
    ///
    /// This establishes scalar traction/net-flow coordinates. It does not
    /// assert that the resolved opening velocity is spatially uniform; an
    /// exterior constant-velocity lift remains a declared port approximation.
    pub fn plane_wave_port_triangles(
        &self,
        port: usize,
    ) -> Option<impl ExactSizeIterator<Item = ThermoviscousPortTriangle> + '_> {
        self.plane_wave_port_area_m2(port)?;
        Some(self.physical_modes[port].iter().map(|trace| {
            ThermoviscousPortTriangle {
                vertices_m: trace.vertices.map(|vertex| vertex.map(f64::from_bits)),
                outward_normal: trace.outward_normal,
                area_m2: trace.area_m2,
            }
        }))
    }
}
