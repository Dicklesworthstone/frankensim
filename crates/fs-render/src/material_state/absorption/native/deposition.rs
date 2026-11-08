//! Accepted material absorption drives the original conservative beam owner.

use super::{MaterialAbsorptionField, checkpoint};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::deposition::{BeamPowerDeposition, deposit_prepared};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

impl<'g, 'field, 'material, const N: usize, Model>
    MaterialAbsorptionField<'g, 'field, 'material, N, Model>
{
    /// Deposit a collimated beam through the exact prepared material spectra.
    ///
    /// Each incident lane is watts already integrated over the beam's area and
    /// spectral quadrature or line/bin width. Wavelengths are this field's fixed
    /// nodes; no coefficient interpolation, material evaluation or RGB-to-power
    /// conversion occurs during traversal. The narrow beam must see the same
    /// piecewise-constant material across its cross-section.
    ///
    /// The original native DDA and stable Beer power partition assign absorbed
    /// watts to the original density cells. Incident = transmitted + absorbed
    /// to floating-point roundoff; scattering is exactly zero in this model.
    /// Finite endpoints exclude downstream material. The geometry's traversal
    /// ceiling and the independent output-cell limit remain binding. An
    /// admitted vacuum traverses no cells, while invalid rays and powers still
    /// refuse. Cancellation or any refusal publishes no partial report.
    ///
    /// These are absorbed INCIDENT watts, even in accepted LTE emission mode;
    /// thermal radiation/cooling is a separate source. The report retains this
    /// field's identity, original source-band intersection and exact nodes. It
    /// can be passed directly to the existing optical-heating observer with an
    /// explicit optical-to-thermal cell map, within that owner's admission.
    /// Reprepare the material field and deposition after a source/state change;
    /// no optical feedback or thermal material update is inferred here.
    pub fn deposit_beam_power(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
        incident_power_w: &[f64; N],
        max_deposition_cells: usize,
    ) -> Result<BeamPowerDeposition<N>, VolumeTransportError> {
        checkpoint(cx)?;
        // The native owner validates origin and direction below. Validate the
        // original distance before shortening a proven vacuum to an empty ray.
        if max_distance_m.is_nan() || max_distance_m < 0.0 {
            return Err(VolumeTransportError::Ray);
        }
        let distance = if self.vacuum { 0.0 } else { max_distance_m };
        deposit_prepared(
            &mut || checkpoint(cx),
            self.geometry,
            (origin, unit_direction, distance),
            &self.wavelengths_nm,
            incident_power_w,
            max_deposition_cells,
            (
                self.identity,
                (self.source_bounds_nm[0], self.source_bounds_nm[1]),
            ),
            |cell| {
                let absorption =
                    self.models[self.materials.cell_indices()[cell] as usize].absorption;
                Ok((absorption, [0.0; N], absorption))
            },
        )
    }
}

#[cfg(test)]
mod tests;
