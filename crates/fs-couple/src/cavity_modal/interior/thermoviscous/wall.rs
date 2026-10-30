//! Source-evaluated finite solids coupled to resolved gas temperature.
//!
//! Every patch keeps its original material sample, geometry, depth and backing.
//! The original mixed FE owns the gas heat equation and the normal solid heat
//! admittance. No gas boundary-layer loss is added by this adapter.

use super::{ThermoviscousInterior, ThermoviscousSourceError, invalid};
use crate::cavity_modal::interior::conductive_wall_transport;
use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalWallDepth,
};
use fs_feec::acoustics::thermoviscous::{ThermoviscousError, wall::ConductiveWallRegion};
use fs_material::thermal::transport::ThermalTransportSample;

/// A disjoint patch of original no-slip boundary faces and its evaluated solid.
/// All dimensions are SI. The sample's temperature must equal the source gas
/// bitwise. The original solid owner admits its local normal-conduction model
/// at each requested frequency; gas thermal layers remain spatially resolved.
pub struct MaterialWallRegion<'a> {
    pub faces: &'a [[u32; 3]],
    pub sample: &'a ThermalTransportSample,
    pub depth: ThermalWallDepth,
    pub minimum_tangential_scale_m: f64,
    pub geometry: BoundaryLayerGeometry,
}

impl<'c> ThermoviscousInterior<'c> {
    /// Attach all finite wall patches once. Faces cannot overlap another
    /// patch, a pressure port or an exact isothermal wall. Unassigned no-slip
    /// walls remain insulated. Source samples are retained in patch order.
    ///
    /// This consumes the unprepared model; no borrowed frequency factorization
    /// can outlive a material change. Solid transport/profile admission is
    /// frequency-local in the original owner, and its errors remain typed.
    pub fn with_conductive_walls(
        mut self,
        regions: &[MaterialWallRegion<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousSourceError> {
        if cancelled() {
            return Err(ThermoviscousError::Cancelled.into());
        }
        if regions.is_empty() || !self.wall_sources.is_empty() {
            return Err(invalid(
                "source conductive wall attachment requires a nonempty first patch set",
            ));
        }
        let mut numerical = Vec::with_capacity(regions.len());
        let mut sources = Vec::with_capacity(regions.len());
        for region in regions {
            if cancelled() {
                return Err(ThermoviscousError::Cancelled.into());
            }
            let transport = conductive_wall_transport(&self.gas, region.sample)
                .map_err(ThermoviscousError::from)?;
            numerical.push(ConductiveWallRegion {
                faces: region.faces.to_vec(),
                wall: ConductiveThermalWall {
                    transport,
                    depth: region.depth,
                    minimum_tangential_scale_m: region.minimum_tangential_scale_m,
                },
                geometry: region.geometry,
                temperature_k: region.sample.temperature_k(),
            });
            sources.push(*region.sample);
        }
        self.fe = self.fe.with_conductive_walls(&numerical, cancelled)?;
        self.wall_sources = sources;
        Ok(self)
    }
}
