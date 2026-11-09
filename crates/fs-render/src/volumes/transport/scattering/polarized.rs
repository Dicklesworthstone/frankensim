//! Molecular and spherical-particle Stokes transport on the original flights.
//!
//! The Hansen--Travis Rayleigh matrix (doi:10.1007/BF00168069) uses the same
//! unpolarized 90-degree depolarization ratio and normalization as `PhaseFunction`.
//! Here +Q is perpendicular to the scattering plane; (s,p,travel) is right-handed,
//! U=2 Re(Es conj(Ep)), V=2 Im(Es conj(Ep)), and both directions are forward travel.
//! No scalar HG, isotropic or intensity-only particle table determines a Mueller
//! matrix: those laws refuse instead of silently erasing or preserving Q/U/V.

use super::{PhaseFunction, ScatteringEvent, ScatteringGrid, checkpoint, unit_variate, validate_ray};
use crate::polarization::{Mueller, Stokes};
use crate::polarization::dielectric::transport::{
    BeamError, BeamMedium, PolarizationFrame, PolarizedBeam,
};
use crate::volumes::transport::VolumeTransportError;
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

// Reciprocal scalar camera/source connections retain the full Stokes state.
mod connection;

/// A refusal is not absorption, a rejected sample, or a zero detector reading.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PolarizedScatteringError {
    /// The scalar angular law does not retain its polarized response.
    UnsupportedPhase,
    /// Native grids are index-one. Their carrier must be homogeneous, index-one
    /// and clear at the transported wavelength; the grid owns ALL extinction.
    CarrierMedium,
    Beam(BeamError),
    Volume(VolumeTransportError),
}
impl From<BeamError> for PolarizedScatteringError {
    fn from(value: BeamError) -> Self { Self::Beam(value) }
}
impl From<VolumeTransportError> for PolarizedScatteringError {
    fn from(value: VolumeTransportError) -> Self { Self::Volume(value) }
}
impl core::fmt::Display for PolarizedScatteringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized scattering: {self:?}")
    }
}
impl core::error::Error for PolarizedScatteringError {}

impl PhaseFunction {
    /// Differential Stokes scattering matrix per steradian in the common s/p
    /// frames. M00 equals the existing scalar phase PDF; its sphere integral
    /// is one. The intensity integral is one for EVERY incident polarization.
    /// At 90 degrees unpolarized output obeys I_parallel/I_perpendicular=delta.
    /// Molecular anisotropy also changes U/V; it is not an isotropic depolarizer.
    pub fn rayleigh_mueller(&self, cosine: f64) -> Result<Mueller, PolarizedScatteringError> {
        let delta = depolarization(self)?;
        let scalar = self.eval_cos(cosine)?;
        let c = cosine.clamp(-1.0, 1.0);
        let scale = 3.0 / (8.0 * core::f64::consts::PI * (2.0 + delta));
        let linear = 1.0 - delta;
        let off_diagonal = scale * linear * ((1.0 - c) * (1.0 + c));
        Ok(Mueller([
            [scalar, off_diagonal, 0.0, 0.0],
            [off_diagonal, scale * linear * (1.0 + c*c), 0.0, 0.0],
            [0.0, 0.0, scale * 2.0 * linear * c, 0.0],
            [0.0, 0.0, 0.0, scale * 2.0 * (1.0 - 2.0*delta) * c],
        ]))
    }

    /// Reuse the exact native angular CDF and solid-angle PDF. The returned
    /// Stokes power is M/PDF times the incoming power, NOT M00 times it again.
    /// A polarized sample can weigh up to twice its input; conservation holds
    /// in expectation, so clamping individual weights to one would bias paths.
    pub fn sample_polarized(
        &self,
        beam: PolarizedBeam,
        samples: [f64; 2],
    ) -> Result<PolarizedPhaseSample, PolarizedScatteringError> {
        self.admit_polarized()?;
        let incoming = beam.frame().direction();
        let sample = self.sample_packet(&[beam.wavelength_nm()], incoming, samples[0], samples[1])?;
        if !sample.pdf.is_finite() || sample.pdf <= 0.0 {
            return Err(VolumeTransportError::Unrepresentable.into());
        }
        let transverse = cross(incoming, sample.direction);
        let axis = if transverse.x == 0.0 && transverse.y == 0.0 && transverse.z == 0.0 {
            // Forward/backward scattering has no unique plane. Retain the
            // incoming horizontal axis; the signed U/V entries handle reversal.
            beam.frame().horizontal()
        } else { unit(transverse)? };
        let incident = beam.with_reference_axis(axis)?;
        let power = if matches!(self, Self::Rayleigh | Self::RayleighDepolarized { .. }) {
            // Retain the admitted molecular sampler's original arithmetic.
            let mut matrix=self.rayleigh_mueller(incoming.dot(sample.direction))?;
            for row in &mut matrix.0 { for entry in row { *entry/=sample.pdf; } }
            let out=transform(matrix,incident.power())?;
            if out.i==0.0 && incident.power().i>0.0 {
                // A rounded M00-M01 can erase a representable near-dipole-null
                // channel. Use the factored physical response only in that
                // case; ordinary molecular samples retain their original bits.
                sampled_response(self,beam.wavelength_nm(),incoming.dot(sample.direction),incident.power(),sample.pdf)?
            } else { out }
        } else if incident.power().i==0.0 {
            Stokes::unpolarized(0.0)
        } else {
            sampled_response(self,beam.wavelength_nm(),incoming.dot(sample.direction),incident.power(),sample.pdf)?
        };
        let frame = PolarizationFrame::try_new(sample.direction, axis)?;
        let beam = PolarizedBeam::try_new(
            beam.ray().origin, frame, beam.wavelength_nm(), beam.medium(), power,
        )?;
        Ok(PolarizedPhaseSample { beam, pdf: sample.pdf })
    }

    /// Sample the physical polarized angular density without changing packet
    /// power. This is the matched continuation needed by conservative heating
    /// and momentum histories; `sample_polarized` instead returns importance
    /// weights that can exceed the incoming power.
    ///
    /// Rejection reuses the original scalar CDF and full Stokes response.
    /// Its envelope is `1 + degree_of_linear_polarization`: the differential
    /// intensity is M00 I + M01 Q in the scattering frame, and |M01| <= M00.
    /// Thus at most two proposals are needed on average. Circular/unpolarized
    /// inputs need one proposal and no acceptance draw. Every trial polls the
    /// caller's checkpoint; exhausting the explicit cap refuses the operation.
    /// The returned PDF is the actual polarized density, not the scalar proposal.
    pub fn sample_polarized_power_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        beam: PolarizedBeam,
        max_proposals: u32,
        random: &mut impl FnMut() -> f64,
    ) -> Result<PolarizedPhaseSample, PolarizedScatteringError> {
        check()?;
        self.admit_polarized()?;
        if max_proposals == 0 {
            return Err(VolumeTransportError::Budget.into());
        }
        let input = beam.power();
        if input.i <= 0.0 {
            return Err(BeamError::InvalidStokes.into());
        }
        // Unit intensity avoids overflowing a large physical packet while
        // evaluating the at-most-two proposal weight.
        let normalized = Stokes::try_new(
            1.0, input.q / input.i, input.u / input.i, input.v / input.i,
        ).ok_or(BeamError::InvalidStokes)?;
        let unit_beam = PolarizedBeam::try_new(
            beam.ray().origin, beam.frame(), beam.wavelength_nm(), beam.medium(), normalized,
        )?;
        let linear = det::hypot(normalized.q, normalized.u).min(1.0);
        let envelope = 1.0 + linear;
        for _ in 0..max_proposals {
            check()?;
            let proposal = self.sample_polarized(unit_beam, [random(), random()])?;
            let response = proposal.beam.power();
            let weight = response.i;
            if !weight.is_finite() || weight > envelope * (1.0 + 128.0*f64::EPSILON) {
                return Err(BeamError::NumericalRange.into());
            }
            if linear != 0.0 {
                let accept = random();
                if !unit_variate(accept) { return Err(VolumeTransportError::Sample.into()); }
                if accept * envelope >= weight { continue; }
            }
            if weight <= 0.0 {
                return Err(BeamError::NumericalRange.into());
            }
            let pdf = proposal.pdf * weight;
            if !pdf.is_finite() || pdf <= 0.0 {
                return Err(VolumeTransportError::Unrepresentable.into());
            }
            let power = Stokes::try_new(
                input.i,
                (response.q / weight) * input.i,
                (response.u / weight) * input.i,
                (response.v / weight) * input.i,
            ).ok_or(BeamError::NumericalRange)?;
            let beam = PolarizedBeam::try_new(
                beam.ray().origin, proposal.beam.frame(), beam.wavelength_nm(), beam.medium(), power,
            )?;
            check()?;
            return Ok(PolarizedPhaseSample { beam, pdf });
        }
        Err(VolumeTransportError::Budget.into())
    }
}

/// Sampled phase continuation; PDF is a solid-angle density, never a delta.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedPhaseSample {
    pub beam: PolarizedBeam,
    pub pdf: f64,
}

/// Borrowed, pre-admitted polarized transport. Keeps original cell identities,
/// spectra, density, source bounds, and native traversal/work limits. Preparation
/// scans the optical palette once; path samples allocate and clone nothing.
#[derive(Clone, Copy)]
pub struct RayleighGrid<'grid, 'data> {
    grid: &'grid ScatteringGrid<'data>,
}

impl<'grid, 'data> RayleighGrid<'grid, 'data> {
    /// Admit all active optical models before a path starts. Pure absorption
    /// requires no polarization law. An indexed cell retains its own molecular
    /// depolarization or spherical response, including the original spectral
    /// scattering weights of independent material mixtures.
    pub fn try_new(
        cx: &Cx<'_>, grid: &'grid ScatteringGrid<'data>,
    ) -> Result<Self, PolarizedScatteringError> {
        checkpoint(cx)?;
        if grid.has_scattering() {
            if let Some(palette) = grid.scattering_palette() {
                for optical in palette.values() {
                    checkpoint(cx)?;
                    if !optical.scattering.is_zero() { optical.phase.admit_polarized()?; }
                }
            } else {
                grid.phase_ref().ok_or(PolarizedScatteringError::UnsupportedPhase)?.admit_polarized()?;
            }
        }
        checkpoint(cx)?;
        Ok(Self { grid })
    }

    /// Original immutable physical snapshot, not a retabulated optical medium.
    pub const fn grid(&self) -> &'grid ScatteringGrid<'data> { self.grid }

    /// Transport the incident source's power, using one native monochromatic
    /// free flight followed by a full Stokes collision when present. Inputs are
    /// [distance survival variate, phase cosine variate, phase azimuth variate].
    /// All three must be in [0,1), even when the selected segment escapes.
    ///
    /// Grid extinction replaces carrier absorption, never adds to it. Do not
    /// call `beam.advance` or apply transmittance again to the returned beam.
    /// Native rebased event positions retain cell-face ownership. Escape ends
    /// at the requested physical endpoint. Source-only forward power contains
    /// no LTE radiance term; thermal emission requires separately launched rays.
    pub fn sample_beam_segment(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        max_distance_m: f64,
        samples: [f64; 3],
    ) -> Result<PolarizedFreeFlight, PolarizedScatteringError> {
        checkpoint(cx)?;
        validate_carrier(beam)?;
        if !max_distance_m.is_finite() { return Err(VolumeTransportError::Ray.into()); }
        validate_ray(beam.ray().origin, beam.ray().dir, max_distance_m)?;
        if !samples.into_iter().all(unit_variate) {
            return Err(VolumeTransportError::Sample.into());
        }
        let wavelengths = [beam.wavelength_nm()];
        let (a, s, t) = self.grid.coefficients(&wavelengths)?;
        let flight = self.grid.sample_segment_prepared(
            cx, beam.ray().origin, beam.ray().dir, max_distance_m, &wavelengths,
            0.5, samples[0], (&a, &s, &t), false,
        )?;
        let (position, distance_m) = flight.event.map_or(
            (beam.ray().at(max_distance_m), max_distance_m),
            |event| (event.position, event.distance_m),
        );
        let power = scale_power(beam.power(), flight.weight[0])?;
        let mut continuation = PolarizedBeam::try_new(
            position, beam.frame(), beam.wavelength_nm(), beam.medium(), power,
        )?;
        let mut phase_pdf = None;
        if let Some(event) = flight.event
            && continuation.power().i > 0.0
        {
            let phase = self.grid.phase_at_cell(event.cell_index)?;
            let scattered = phase.sample_polarized(continuation, [samples[1], samples[2]])?;
            continuation = scattered.beam;
            phase_pdf = Some(scattered.pdf);
        }
        checkpoint(cx)?;
        Ok(PolarizedFreeFlight {
            beam: continuation, event: flight.event, distance_m, phase_pdf,
            escape_probability: flight.escape_probability,
        })
    }
}

/// One physical segment. Event=None reaches the requested endpoint. Extinction
/// and importance weights are ALREADY in beam.power(), including absorption-only
/// segments and vacuum clipping. A collision with zero surviving power requires
/// no angular response; phase_pdf is present only when a phase was sampled.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedFreeFlight {
    pub beam: PolarizedBeam,
    pub event: Option<ScatteringEvent>,
    pub distance_m: f64,
    pub phase_pdf: Option<f64>,
    pub escape_probability: f64,
}

pub(crate) fn validate_carrier(beam: PolarizedBeam) -> Result<(), PolarizedScatteringError> {
    let BeamMedium::Homogeneous(material) = beam.medium() else {
        return Err(PolarizedScatteringError::CarrierMedium);
    };
    let wavelength = beam.wavelength_nm();
    if material.ior().eval(wavelength).map_err(BeamError::from)? != 1.0
        || material.absorption_ref().extinction_per_m(wavelength).map_err(BeamError::from)? != 0.0 {
        return Err(PolarizedScatteringError::CarrierMedium);
    }
    Ok(())
}

fn depolarization(phase: &PhaseFunction) -> Result<f64, PolarizedScatteringError> {
    match phase {
        PhaseFunction::Rayleigh => Ok(0.0),
        PhaseFunction::RayleighDepolarized { depolarization: d }
            if d.is_finite() && (0.0..6.0/7.0).contains(d) => Ok(*d),
        PhaseFunction::RayleighDepolarized { .. } => Err(VolumeTransportError::Phase.into()),
        _ => Err(PolarizedScatteringError::UnsupportedPhase),
    }
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
}
fn unit(v: Vec3) -> Result<Vec3, BeamError> {
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if !scale.is_finite() || scale == 0.0 { return Err(BeamError::InvalidFrame); }
    let v = Vec3::new(v.x/scale, v.y/scale, v.z/scale);
    Ok(v.scale(1.0/det::sqrt(v.dot(v))))
}
fn scale_power(input: Stokes, weight: f64) -> Result<Stokes, BeamError> {
    if !weight.is_finite() || weight < 0.0 { return Err(BeamError::NumericalRange); }
    let out = Stokes { i:input.i*weight, q:input.q*weight, u:input.u*weight, v:input.v*weight };
    if ![out.i,out.q,out.u,out.v].into_iter().all(f64::is_finite)
        || (input.i>0.0 && weight>0.0 && out.i==0.0) {
        return Err(BeamError::NumericalRange);
    }
    Ok(out)
}

fn sampled_response(phase:&PhaseFunction,wavelength_nm:f64,cosine:f64,input:Stokes,pdf:f64) -> Result<Stokes,PolarizedScatteringError> {
    let response=phase.apply_polarized_at(wavelength_nm,cosine,
        Stokes { i:1.0,q:input.q/input.i,u:input.u/input.i,v:input.v/input.i })?;
    Ok(scale_power(Stokes { i:response.i/pdf,q:response.q/pdf,u:response.u/pdf,v:response.v/pdf },input.i)?)
}
fn transform(matrix: Mueller, input: Stokes) -> Result<Stokes, BeamError> {
    if input.i == 0.0 { return Ok(Stokes::unpolarized(0.0)); }
    let mut out = matrix.apply(Stokes {
        i:1.0, q:input.q/input.i, u:input.u/input.i, v:input.v/input.i,
    });
    let norm = det::hypot(det::hypot(out.q,out.u),out.v);
    let tolerance = 128.0*f64::EPSILON*out.i.abs().max(1.0);
    if ![out.i,out.q,out.u,out.v,norm].into_iter().all(f64::is_finite)
        || out.i < -tolerance || norm > out.i.max(0.0)+tolerance {
        return Err(BeamError::NumericalRange);
    }
    out.i = out.i.max(0.0);
    if norm > out.i {
        let factor = out.i/norm;
        out.q *= factor; out.u *= factor; out.v *= factor;
    }
    scale_power(out,input.i)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod particles;
