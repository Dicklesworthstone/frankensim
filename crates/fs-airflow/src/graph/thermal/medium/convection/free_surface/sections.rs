//! One physical leading edge across the original liquid surface's axial strips.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceGasFace;

/// One original surface patch occupying an authored streamwise plate strip.
///
/// The patch's actual area, distances, liquid/gas identities and outward normal
/// must match every trial. `start` and `end` are distances from the SAME physical
/// leading edge, measured along `downstream_direction`, not distances from this
/// patch's own upstream side. They describe its streamwise extent; its original
/// area includes the transverse width. Several patches at the same axial extent
/// can share these distances. The caller owns this geometry correspondence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarFreeSurfaceSection {
    pub face: FreeSurfaceGasFace,
    pub outward_normal: [f64; 3],
    pub start: Length,
    pub end: Length,
    /// Nonzero downstream direction tangent to the actual surface. It need not
    /// be normalized. Relative flow must follow this direction, including sign.
    pub downstream_direction: [f64; 3],
    /// Actual gas velocity RELATIVE to the liquid [m/s], held during the split.
    pub relative_velocity_m_s: [f64; 3],
    pub regime: PlanarForcedConvectionRegime,
}

/// Section-dependent convection in original call-local surface-contact order.
///
/// Lookup is O(1) and allocates no geometry or replacement heat state. Each
/// trial uses the original gas at its actual surface/bulk mean temperature.
/// The lower owner still solves liquid half-cell resistance, radiation and the
/// reciprocal finite-gas energy exchange, with the same timestep and rollback.
///
/// Only a stationary set of authored patches is admitted. Newly exposed PLIC
/// patches need an explicit physical section binding; index zero cannot stand
/// in for all streamwise positions. Likewise, flow reversal needs a rebind to
/// the opposite leading edge. This adapter does not resolve a non-isothermal
/// gas boundary layer, curvature, gas momentum or transition between cards.
#[derive(Clone, Copy, Debug)]
pub struct PlanarFreeSurfaceSections<'a> {
    pub sections: &'a [PlanarFreeSurfaceSection],
}

impl FreeSurfaceGasConvection for PlanarFreeSurfaceSections<'_> {
    fn coefficient(
        &self,
        query: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        let section = self.sections.get(query.patch_index).ok_or_else(|| {
            invalid("surface convection has no authored streamwise section for this patch")
        })?;
        if section.face != query.surface.face {
            return Err(invalid(
                "streamwise convection must retain its original patch geometry and owner",
            ));
        }
        let normal = admit_plate_query(&query, section.outward_normal, section.end)?;
        let (downstream, _) = direction(section.downstream_direction)?;
        let (flow, speed) = physical_vector(section.relative_velocity_m_s)?;
        if !tangent(normal, downstream) || !tangent(normal, flow) || !aligned(flow, downstream) {
            return Err(invalid(
                "plate section distances require tangent relative flow in their declared downstream direction",
            ));
        }
        let correlation = match section.regime {
            PlanarForcedConvectionRegime::Laminar => CorrelationId::FlatPlateLaminarAverage,
            PlanarForcedConvectionRegime::Turbulent => CorrelationId::FlatPlateTurbulentAverage,
        };
        let result = evaluate_flat_plate_section_from_material(
            correlation,
            query.film,
            speed,
            section.start,
            section.end,
        )
        .map_err(original)?;
        Ok(Some(result.into_coefficient()))
    }
}
