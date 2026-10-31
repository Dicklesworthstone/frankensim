//! Source-evaluated finite solids coupled to resolved fluid temperature.
//!
//! Every patch keeps its original material sample, geometry, depth and backing.
//! The original mixed FE owns the gas heat equation and the normal solid heat
//! admittance. No gas boundary-layer loss is added by this adapter.

use super::{ThermoviscousInterior, ThermoviscousSourceError, invalid};
use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalWallDepth, WallThermalTransport,
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

impl<'c, S> ThermoviscousInterior<'c, S> {
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
            // The complete fluid source was admitted at construction. Keep
            // the SAME physical temperature for the original wall transport;
            // do not require liquid/real-fluid coefficients to impersonate gas.
            if self.fe.medium().temperature_k.to_bits() != region.sample.temperature_k().to_bits() {
                return Err(invalid(
                    "conductive acoustic wall and fluid require one equilibrium temperature",
                ));
            }
            let transport = WallThermalTransport {
                conductivity_w_m_k: region.sample.conductivity_w_m_k(),
                volumetric_heat_capacity_j_m3_k: region.sample.volumetric_heat_capacity_j_m3_k(),
            };
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


/// Source paths for one explicitly transversely isotropic planar solid.
/// `normal` must be the exact prepared path that produced the attached normal
/// wall sample. `tangential` declares isotropy in the plane and zero coupling to
/// the normal, using its OWN conductivity source. Both paths must select the
/// same material card, state path, density and specific-heat law/constraint.
pub struct MaterialTangentialWallRegion<'a> {
    pub normal: &'a fs_material::thermal::transport::PreparedThermalTransport,
    pub tangential: &'a fs_material::thermal::transport::PreparedThermalTransport,
}

impl<'c, S> ThermoviscousInterior<'c, S> {
    /// Attach in-plane conductivity in complete original wall-region order.
    /// This is an explicit material-symmetry declaration, not an inference from
    /// a scalar normal measurement. Source identities are retained beside the
    /// original normal samples. Failure leaves no partially attached owner.
    ///
    /// `prepare_time_with_wall_diffusion` now includes spatial solid conduction.
    /// Harmonic normal-only response methods refuse this material attachment.
    /// Actual coplanarity, outward normals and finite thickness are admitted by
    /// the original causal field preparation; no curved extrusion is invented.
    pub fn with_tangentially_isotropic_wall_transport(
        mut self,
        regions: &[MaterialTangentialWallRegion<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousSourceError> {
        if cancelled() {
            return Err(ThermoviscousError::Cancelled.into());
        }
        if regions.is_empty() || regions.len() != self.wall_sources.len()
            || !self.wall_tangential_sources.is_empty()
        {
            return Err(invalid("complete first tangential wall source set required"));
        }
        let mut sources = Vec::with_capacity(regions.len());
        for (binding, normal_sample) in regions.iter().zip(&self.wall_sources) {
            if cancelled() {
                return Err(ThermoviscousError::Cancelled.into());
            }
            if binding.normal.identity() != normal_sample.path_identity()
                || binding.normal.material_card_identity() != binding.tangential.material_card_identity()
                || binding.normal.query_endpoints() != binding.tangential.query_endpoints()
                || binding.normal.source_lines()[..2] != binding.tangential.source_lines()[..2]
            {
                return Err(invalid(
                    "normal and tangential transport require the same original card, state and caloric sources",
                ));
            }
            let temperature = normal_sample.temperature_k();
            let expected = binding.normal.at_temperature(temperature)
                .map_err(|_| invalid("normal wall temperature outside original material path"))?;
            let tangent = binding.tangential.at_temperature(temperature)
                .map_err(|_| invalid("tangential wall temperature outside original material path"))?;
            if expected != *normal_sample
                || tangent.density_kg_m3().to_bits() != normal_sample.density_kg_m3().to_bits()
                || tangent.specific_heat_j_kg_k().to_bits() != normal_sample.specific_heat_j_kg_k().to_bits()
            {
                return Err(invalid("spatial wall transport differs from its original thermal state"));
            }
            sources.push(tangent);
        }
        self.wall_tangential_sources = sources;
        Ok(self)
    }

    /// Original in-plane source samples, empty on the original normal-only path.
    pub fn wall_tangential_sources(&self) -> &[ThermalTransportSample] {
        &self.wall_tangential_sources
    }
}
