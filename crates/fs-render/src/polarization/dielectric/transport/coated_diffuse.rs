//! Polarized smooth coating reflection over an ideal depolarizing substrate.

use super::{BeamError, BeamSample, PolarizationFrame, PolarizedBeam, cross, unit};
use super::connection::{SurfaceConnectionResponse, checked_product, normalized_camera, reciprocal_contraction};
use crate::coated_diffuse::SmoothDiffuseCoating;
use crate::dielectric::DielectricEvent;
use crate::polarization::{Mueller, Stokes};
use crate::polarization::dielectric::fresnel_dielectric_polarized;
use fs_geom::Vec3;

impl PolarizedBeam {
    /// Sample the mirror delta or diffuse return from a smooth dielectric coat.
    ///
    /// The first two uniform coordinates drive the original cosine-hemisphere
    /// proposal; the third selects the Fresnel branch. All lie in `[0,1)`.
    /// The diffuse density includes its branch mass. Delta `pdf` is a discrete
    /// mass, and neither branch changes the surrounding medium or ray origin.
    ///
    /// Entry retains its full s/p response, the substrate depolarizes at each
    /// encounter, all internal returns are summed, and escape can repolarize
    /// the light. Optional gray Beer absorption retains both internal path
    /// lengths and every return. The geometrically collapsed coat adds no
    /// displacement, flight time, coherent-film phase or inferred index.
    pub fn sample_coated_diffuse(
        self,
        normal: Vec3,
        coating: SmoothDiffuseCoating,
        reflectance: f64,
        samples: [f64; 3],
    ) -> Result<BeamSample, BeamError> {
        if !(0.0..=1.0).contains(&reflectance)
            || samples.iter().any(|value| !(0.0..1.0).contains(value)) {
            return Err(BeamError::InvalidDiffuseSample);
        }
        let normal = self.coated_normal(normal, coating)?;
        if coating.is_matched() && coating.normal_optical_depth() == 0.0 {
            let sample = self.sample_depolarizing_diffuse(normal, reflectance, [samples[0], samples[1]])?;
            return Ok(BeamSample { beam: sample.beam, pdf: sample.pdf_solid_angle,
                event: DielectricEvent::Reflection, delta: false });
        }
        let cosine = (-normal.dot(self.frame.direction)).clamp(0.0, 1.0);
        let fresnel = fresnel_dielectric_polarized(cosine, coating.exterior_index(), coating.coating_index())?;
        let probability = fresnel.scalar.event_probabilities();
        // Sample the smaller branch next to zero. This retains weak diffuse
        // support even when the mirror probability rounds to one.
        let mirror = if probability[0] <= probability[1] {
            samples[2] < probability[0]
        } else { samples[2] >= probability[1] };
        if mirror {
            let direction = Vec3::new(self.frame.direction.x + 2.0*cosine*normal.x,
                self.frame.direction.y + 2.0*cosine*normal.y,
                self.frame.direction.z + 2.0*cosine*normal.z);
            let matrix = Mueller(fresnel.reflection.0.map(|row| row.map(|value| value / probability[0])));
            let beam = self.scatter(normal, direction, matrix, DielectricEvent::Reflection, self.medium)?;
            return Ok(BeamSample { beam, event: DielectricEvent::Reflection,
                pdf: probability[0], delta: true });
        }
        let proposal = self.sample_depolarizing_diffuse(normal, 1.0,
            [samples[0] * coating.diffuse_cone_fraction(), samples[1]])?;
        let beam = self.coated_diffuse_response(normal, proposal.beam.frame.direction,
            coating, reflectance, probability[1], true)?;
        let pdf = checked_product(probability[1], coating.diffuse_direction_pdf(
            normal.dot(proposal.beam.frame.direction).clamp(0.0, 1.0))?)?;
        Ok(BeamSample { beam, event: DielectricEvent::Reflection, pdf, delta: false })
    }

    /// Evaluate the continuous coated-substrate radiance response and its
    /// original mixed proposal at an independently supplied source direction.
    /// The mirror delta is absent from this solid-angle response.
    pub fn evaluate_coated_diffuse_connection(
        self,
        normal: Vec3,
        coating: SmoothDiffuseCoating,
        reflectance: f64,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
    ) -> Result<SurfaceConnectionResponse, BeamError> {
        if !(0.0..=1.0).contains(&reflectance) { return Err(BeamError::InvalidDiffuseSample); }
        let camera = normalized_camera(self, normalized_source)?;
        let normal = camera.coated_normal(normal, coating)?;
        let direction = source_frame.direction().scale(-1.0);
        let outgoing_cosine = normal.dot(direction);
        if outgoing_cosine <= 0.0 {
            return Ok(SurfaceConnectionResponse::zero(DielectricEvent::Reflection));
        }
        let response = camera.coated_diffuse_response(normal, direction, coating, reflectance,
            core::f64::consts::PI, false)?;
        let factor = reciprocal_contraction(response, source_frame, normalized_source)?;
        let probability = coating.branch_probabilities((-normal.dot(camera.frame.direction)).clamp(0.0, 1.0))?[1];
        let pdf = checked_product(probability,
            coating.diffuse_direction_pdf(outgoing_cosine.clamp(0.0, 1.0))?)?;
        SurfaceConnectionResponse::checked(factor, pdf, DielectricEvent::Reflection)
    }

    fn coated_normal(self, normal: Vec3, coating: SmoothDiffuseCoating) -> Result<Vec3, BeamError> {
        if self.medium.phase_index_at(self.origin, self.wavelength_nm)? != coating.exterior_index() {
            return Err(BeamError::MediumMismatch);
        }
        let normal = self.incident_normal(normal)?;
        if normal.dot(self.frame.direction) >= 0.0 { return Err(BeamError::InvalidFrame); }
        Ok(normal)
    }

    fn coated_diffuse_response(
        self, normal: Vec3, direction: Vec3, coating: SmoothDiffuseCoating,
        reflectance: f64, divisor: f64, sampled: bool,
    ) -> Result<Self, BeamError> {
        let incident_cosine = (-normal.dot(self.frame.direction)).clamp(0.0, 1.0);
        let outgoing_cosine = normal.dot(direction).clamp(0.0, 1.0);
        let incident_axis = interface_axis(self.frame.direction, normal)?;
        let incident = self.with_reference_axis(incident_axis)?;
        let entry = fresnel_dielectric_polarized(incident_cosine, coating.exterior_index(), coating.coating_index())?;
        let exit = fresnel_dielectric_polarized(outgoing_cosine, coating.exterior_index(), coating.coating_index())?;
        // An external TIR direction has exactly zero diffuse response. This
        // is a physical support boundary, not an underflow of positive power.
        if entry.scalar.total_internal_reflection || exit.scalar.total_internal_reflection {
            let frame = PolarizationFrame::try_new(direction, interface_axis(direction, normal)?)?;
            return Self::try_new(self.origin, frame, self.wavelength_nm, self.medium, Stokes::unpolarized(0.0));
        }
        let weight = if sampled { coating.sampled_diffuse_weight(reflectance)? }
            else { coating.diffuse_weight(reflectance)? };
        let weight = checked_product(weight, coating.path_pair_transmittance(
            entry.scalar.transmitted_cosine, exit.scalar.transmitted_cosine))?;
        // A sampled transmission is conditional on its branch. Normalize the
        // native coefficients before applying them: weight/probability can
        // overflow at grazing, and T*input can underflow before /probability
        // restores a finite conditional power. The connection's pi divisor
        // has no small sampling probability and stays in the bounded weight.
        let (transmission, weight) = if divisor < 1.0 {
            (Mueller(entry.transmission.0.map(|row| row.map(|value| value / divisor))), weight)
        } else { (entry.transmission, weight / divisor) };
        let entered = transmission.apply_diagonal_jones(incident.power)?.i;
        let base = checked_product(entered, weight)?;
        // Power transmission is reciprocal at paired Snell angles. Evaluating
        // from the external direction retains the exact same native s/p law
        // without subtracting a critical-angle cosine inside the escape cone.
        let power = exit
            .transmission
            .apply_diagonal_jones(Stokes::unpolarized(base))?;
        if base > 0.0 && outgoing_cosine > 0.0 && power.i == 0.0 { return Err(BeamError::NumericalRange); }
        let axis = interface_axis(direction, normal)?;
        let frame = PolarizationFrame::try_new(direction, axis)?;
        Self::try_new(self.origin, frame, self.wavelength_nm, self.medium, power)
    }
}

fn interface_axis(direction: Vec3, normal: Vec3) -> Result<Vec3, BeamError> {
    let axis = cross(direction, normal);
    if axis.x == 0.0 && axis.y == 0.0 && axis.z == 0.0 {
        // At normal incidence transmission is isotropic. The least-aligned
        // Cartesian reference has a nondegenerate transverse projection.
        let reference = if direction.x.abs() <= direction.y.abs() && direction.x.abs() <= direction.z.abs() {
            Vec3::new(1.0, 0.0, 0.0)
        } else if direction.y.abs() <= direction.z.abs() { Vec3::new(0.0, 1.0, 0.0) }
        else { Vec3::new(0.0, 0.0, 1.0) };
        return Ok(reference);
    }
    unit(axis)
}

#[cfg(test)]
mod tests;
