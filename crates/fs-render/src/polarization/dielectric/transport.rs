//! Forward polarized ray-power transport through native dielectric materials.
//!
//! A ray carries its transverse frame, vacuum wavelength, Stokes power and
//! current glass together. Intersections still belong to the existing chart/
//! scene backends: intersect `beam.ray()`, advance by the returned metric
//! distance, then sample the interface with the medium on the other side.
//! Reflection keeps the incident medium; transmission replaces it. The caller
//! owns nested-medium topology and geometric self-intersection offsets.
//!
//! These are emitter-to-scene POWER samples, not camera-path radiance. The
//! native sampler's radiance eta-squared factor is removed. Per-sample power
//! can exceed input power because it is divided by a proposal probability;
//! conservation applies to the probability-weighted expectation. This opt-in
//! path does not alter the legacy scalar scene integrators.

use fs_geom::{Point3, Vec3};
use fs_math::det;

use crate::charts::Ray;
use crate::dielectric::{DielectricError, DielectricEvent, DielectricGlass};
use crate::polarization::{Mueller, Stokes};

use super::{attenuate_stokes, sample_smooth_dielectric_polarized};

/// Admission or representability refusal. Failed operations leave the beam unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BeamError {
    /// A direction/reference axis is nonfinite, zero or parallel to its partner.
    InvalidFrame,
    /// Negative intensity, nonfinite components or a nonphysical Stokes state.
    InvalidStokes,
    /// A position, distance or derived ray point is invalid/unrepresentable.
    InvalidPosition,
    /// A native material, wavelength or interface refused the request.
    Dielectric(DielectricError),
    /// A finite physical result cannot be represented in binary64.
    NumericalRange,
}

impl From<DielectricError> for BeamError {
    fn from(error: DielectricError) -> Self {
        Self::Dielectric(error)
    }
}

impl core::fmt::Display for BeamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized beam: {self:?}")
    }
}

impl core::error::Error for BeamError {}

/// Right-handed `(horizontal, vertical, propagation)` Stokes reference frame.
/// Fields are private so direction and polarization axes cannot drift apart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizationFrame {
    direction: Vec3,
    horizontal: Vec3,
    vertical: Vec3,
}

impl PolarizationFrame {
    /// Normalize a finite nonzero direction and project a reference axis onto
    /// its transverse plane. A parallel reference axis has no orientation and
    /// is refused. Scaling either supplied vector does not change the frame.
    pub fn try_new(direction: Vec3, reference_axis: Vec3) -> Result<Self, BeamError> {
        let direction = unit(direction)?;
        let reference = unit(reference_axis)?;
        let vertical = unit(cross(direction, reference))?;
        let horizontal = unit(cross(vertical, direction))?;
        Ok(Self { direction, horizontal, vertical })
    }

    /// Unit physical direction of light propagation, not an away-pointing `wo`.
    #[must_use]
    pub const fn direction(self) -> Vec3 { self.direction }

    /// Unit axis for positive Q.
    #[must_use]
    pub const fn horizontal(self) -> Vec3 { self.horizontal }

    /// Unit axis for negative Q; `horizontal x vertical = direction`.
    #[must_use]
    pub const fn vertical(self) -> Vec3 { self.vertical }

    // Re-express a Stokes state for the SAME propagation direction. Dot
    // products give cos/sin directly: no atan or platform trigonometry.
    fn express_in(self, state: Stokes, target: Self) -> Stokes {
        let c = self.horizontal.dot(target.horizontal);
        let s = self.vertical.dot(target.horizontal);
        let length2 = c * c + s * s;
        let c2 = (c * c - s * s) / length2;
        let s2 = 2.0 * c * s / length2;
        Stokes {
            i: state.i,
            q: c2 * state.q + s2 * state.u,
            u: -s2 * state.q + c2 * state.u,
            v: state.v,
        }
    }
}

/// A monochromatic ray and its transported Stokes POWER weight.
///
/// Wavelengths are vacuum nanometres, coordinates and traveled distances are
/// metres. Power may be W, W/nm, or a dimensionless ray weight, but callers must
/// use the same convention for a whole path. Different rays add incoherently;
/// this representation does not model coherent path interference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedBeam {
    origin: Point3,
    frame: PolarizationFrame,
    wavelength_nm: f64,
    medium: DielectricGlass,
    power: Stokes,
}

/// A sampled interface continuation with its original native proposal measure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BeamSample {
    /// Propagating continuation, including its new frame, power and medium.
    pub beam: PolarizedBeam,
    /// Reflection or transmission, independent of the sign of the input normal.
    pub event: DielectricEvent,
    /// Discrete event probability if `delta`, otherwise solid-angle density.
    pub pdf: f64,
    /// Delta probabilities must not be used as solid-angle MIS densities.
    pub delta: bool,
}

impl PolarizedBeam {
    /// Admit a beam without inferring its surrounding medium or polarization.
    pub fn try_new(
        origin: Point3,
        frame: PolarizationFrame,
        wavelength_nm: f64,
        medium: DielectricGlass,
        power: Stokes,
    ) -> Result<Self, BeamError> {
        if !finite_point(origin) { return Err(BeamError::InvalidPosition); }
        medium.ior().eval(wavelength_nm)?;
        if Stokes::try_new(power.i, power.q, power.u, power.v).is_none() {
            return Err(BeamError::InvalidStokes);
        }
        Ok(Self { origin, frame, wavelength_nm, medium,
            power: physical_roundoff(power, power.i)? })
    }

    /// Native chart ray with unit direction: its intersection parameter is metres.
    #[must_use]
    pub const fn ray(self) -> Ray {
        Ray { origin: self.origin, dir: self.frame.direction }
    }

    /// The current physical polarization reference frame.
    #[must_use]
    pub const fn frame(self) -> PolarizationFrame { self.frame }

    /// Stokes power expressed in `frame()`.
    #[must_use]
    pub const fn power(self) -> Stokes { self.power }

    /// Current material, retained on reflection and replaced on transmission.
    #[must_use]
    pub const fn medium(self) -> DielectricGlass { self.medium }

    /// Vacuum wavelength in nanometres.
    #[must_use]
    pub const fn wavelength_nm(self) -> f64 { self.wavelength_nm }

    /// Change coordinates, not the physical beam. Circular polarization and
    /// intensity are invariant; Q/U undergo the required double-angle rotation.
    pub fn with_reference_axis(self, reference_axis: Vec3) -> Result<Self, BeamError> {
        let frame = PolarizationFrame::try_new(self.frame.direction, reference_axis)?;
        let power = physical_roundoff(self.frame.express_in(self.power, frame), self.power.i)?;
        Ok(Self { frame, power, ..self })
    }

    /// Travel within the CURRENT medium using its original spectral absorption
    /// law. This changes both position and all Stokes components, not polarization
    /// degree. No geometric ray-origin epsilon or untraveled distance is inferred.
    pub fn advance(self, distance_m: f64) -> Result<Self, BeamError> {
        if !distance_m.is_finite() || distance_m < 0.0 {
            return Err(BeamError::InvalidPosition);
        }
        let origin = self.ray().at(distance_m);
        if !finite_point(origin) { return Err(BeamError::InvalidPosition); }
        let power = attenuate_stokes(
            self.power, self.medium.absorption(), self.wavelength_nm, distance_m,
        )?;
        Ok(Self { origin, power, ..self })
    }

    /// Sample a smooth interface at the current ray origin. The geometric
    /// normal is normalized and face-forwarded; `target` is explicitly the
    /// medium across this boundary. Uses the native directions/probabilities.
    ///
    /// Input/output s/p frames are constructed automatically. Exactly normal
    /// incidence uses the beam's existing horizontal axis, not an arbitrary
    /// world-axis discontinuity. The returned beam is ready for the next hit.
    pub fn sample_smooth_interface(
        self,
        normal: Vec3,
        target: DielectricGlass,
        event_sample: f64,
    ) -> Result<BeamSample, BeamError> {
        let normal = self.incident_normal(normal)?;
        let eta_i = self.medium.ior().eval(self.wavelength_nm)?;
        let eta_t = target.ior().eval(self.wavelength_nm)?;
        let sample = sample_smooth_dielectric_polarized(
            normal, self.frame.direction.scale(-1.0), eta_i, eta_t, event_sample,
        )?;
        // This is a forward ray-power path. The native camera/radiance
        // proposal is useful, but its eta Jacobian is not a power loss.
        let factor = if sample.event == DielectricEvent::Transmission {
            let ratio = eta_t / eta_i;
            ratio * ratio
        } else { 1.0 };
        let matrix = Mueller(sample.radiance_weight.0.map(|row| row.map(|v| v * factor)));
        let beam = self.scatter(normal, sample.direction, matrix, sample.event, target)?;
        Ok(BeamSample { beam, event: sample.event, pdf: sample.probability, delta: true })
    }

    fn incident_normal(self, normal: Vec3) -> Result<Vec3, BeamError> {
        let normal = unit(normal)?;
        Ok(if normal.dot(self.frame.direction) > 0.0 { normal.scale(-1.0) } else { normal })
    }

    fn scatter(
        self, normal: Vec3, direction: Vec3, matrix: Mueller,
        event: DielectricEvent, target: DielectricGlass,
    ) -> Result<Self, BeamError> {
        let transverse = cross(self.frame.direction, normal);
        let s = if transverse.x == 0.0 && transverse.y == 0.0 && transverse.z == 0.0 {
            self.frame.horizontal
        } else { unit(transverse)? };
        let incident_frame = PolarizationFrame::try_new(self.frame.direction, s)?;
        let input = self.frame.express_in(self.power, incident_frame);
        let output = matrix.apply(input);
        // Bounds the four-term matrix/vector roundoff, including extinction
        // where the correct result is zero. Refuse meaningful violations.
        let magnitude = matrix.0.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
        let reference = self.power.i * magnitude;
        if !reference.is_finite() { return Err(BeamError::NumericalRange); }
        let power = physical_roundoff(output, reference)?;
        let frame = PolarizationFrame::try_new(direction, s)?;
        let medium = if event == DielectricEvent::Transmission { target } else { self.medium };
        Ok(Self { frame, power, medium, ..self })
    }
}

fn unit(v: Vec3) -> Result<Vec3, BeamError> {
    if ![v.x, v.y, v.z].iter().all(|x| x.is_finite()) { return Err(BeamError::InvalidFrame); }
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if scale == 0.0 { return Err(BeamError::InvalidFrame); }
    // Divide componentwise: the reciprocal of a subnormal scale can overflow.
    let v = Vec3::new(v.x / scale, v.y / scale, v.z / scale);
    Ok(v.scale(1.0 / det::sqrt(v.dot(v))))
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
}

fn finite_point(p: Point3) -> bool { [p.x, p.y, p.z].iter().all(|v| v.is_finite()) }

fn physical_roundoff(mut s: Stokes, reference: f64) -> Result<Stokes, BeamError> {
    if ![s.i, s.q, s.u, s.v].iter().all(|x| x.is_finite()) {
        return Err(BeamError::NumericalRange);
    }
    let tolerance = 128.0 * f64::EPSILON * reference;
    let scale = s.q.abs().max(s.u.abs()).max(s.v.abs());
    let norm = if scale == 0.0 { 0.0 } else {
        let q = s.q / scale;
        let u = s.u / scale;
        let v = s.v / scale;
        scale * det::sqrt(q * q + u * u + v * v)
    };
    if !norm.is_finite() || s.i < -tolerance || norm - s.i > tolerance {
        return Err(BeamError::NumericalRange);
    }
    s.i = s.i.max(0.0);
    if s.i == 0.0 { return Ok(Stokes::unpolarized(0.0)); }
    if norm > s.i {
        let ratio = s.i / norm;
        s.q *= ratio; s.u *= ratio; s.v *= ratio;
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};

    fn glass(eta: f64, extinction: f64) -> DielectricGlass {
        DielectricGlass::new(CauchyIor::try_constant(eta).unwrap(),
            BeerLambertAbsorption::try_constant(extinction).unwrap(), GlassProvenance::Custom)
    }
    fn beam(k: Vec3, x: Vec3, medium: DielectricGlass, power: Stokes) -> PolarizedBeam {
        PolarizedBeam::try_new(Point3::new(0.0, 0.0, 0.0),
            PolarizationFrame::try_new(k, x).unwrap(), 550.0, medium, power).unwrap()
    }
    fn near(a: f64, b: f64) { assert!((a - b).abs() <= 3e-12, "{a} != {b}"); }
    fn same(a: Stokes, b: Stokes) {
        near(a.i, b.i); near(a.q, b.q); near(a.u, b.u); near(a.v, b.v);
    }

    #[test]
    fn g0_frames_admit_large_and_subnormal_axes_without_trig() {
        for scale in [f64::MAX, 1.0, f64::from_bits(1)] {
            let f = PolarizationFrame::try_new(Vec3::new(0.0, 0.0, scale),
                Vec3::new(scale, 0.0, 0.0)).unwrap();
            near(f.direction().norm(), 1.0); near(f.horizontal().norm(), 1.0);
            near(f.vertical().dot(f.direction()), 0.0);
            assert_eq!(cross(f.horizontal(), f.vertical()), f.direction());
        }
        assert!(PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 2.0)).is_err());
        assert!(PolarizationFrame::try_new(Vec3::new(f64::NAN, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).is_err());
    }

    #[test]
    fn g3_reframing_changes_coordinates_not_physical_scattering() {
        let b = beam(Vec3::new(0.6, 0.0, -0.8), Vec3::new(0.0, 1.0, 0.0),
            glass(1.0, 0.0), Stokes::try_new(1.0, 0.4, 0.5, -0.2).unwrap());
        let c = b.with_reference_axis(Vec3::new(1.0, 2.0, 3.0)).unwrap();
        for event in [0.0, 0.99] {
            let a = b.sample_smooth_interface(Vec3::new(0.0, 0.0, 1.0), glass(1.5, 0.0), event).unwrap();
            let d = c.sample_smooth_interface(Vec3::new(0.0, 0.0, -1.0), glass(1.5, 0.0), event).unwrap();
            same(a.beam.power(), d.beam.power());
            assert_eq!(a.beam.ray().dir, d.beam.ray().dir);
            assert_eq!(a.pdf, d.pdf);
        }
    }

    #[test]
    fn g1_two_crossed_brewster_reflections_extinguish_the_ray() {
        let ci = 1.0 / 1.5_f64.hypot(1.0);
        let si = 1.5 * ci;
        let b = beam(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0),
            glass(1.0, 0.0), Stokes::unpolarized(1.0));
        let first = b.sample_smooth_interface(Vec3::new(si, 0.0, ci), glass(1.5, 0.0), 0.0).unwrap();
        near(first.beam.power().q, first.beam.power().i);
        let k = first.beam.frame().direction();
        let x = first.beam.frame().horizontal();
        let n = Vec3::new(-ci*k.x + si*x.x, -ci*k.y + si*x.y, -ci*k.z + si*x.z);
        let second = first.beam.sample_smooth_interface(n, glass(1.5, 0.0), 0.0).unwrap();
        assert_eq!(second.event, DielectricEvent::Reflection);
        near(second.beam.power().i, 0.0);
    }

    #[test]
    fn g1_tir_retains_phase_and_normal_reflection_flips_handedness() {
        let t = 0.5_f64.sqrt();
        let b = beam(Vec3::new(t, 0.0, -t), Vec3::new(0.0, -1.0, 0.0),
            glass(1.5, 0.0), Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap());
        let out = b.sample_smooth_interface(Vec3::new(0.0, 0.0, 1.0), glass(1.0, 0.0), 0.5).unwrap();
        same(out.beam.power(), Stokes { i: 1.0, q: 0.0, u: 0.8, v: 0.6 });
        assert_eq!(out.beam.medium(), b.medium());
        let b = beam(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0),
            glass(1.0, 0.0), Stokes::try_new(2.0, 0.0, 1.0, 1.0).unwrap());
        let out = b.sample_smooth_interface(Vec3::new(0.0, 0.0, 1.0), glass(1.5, 0.0), 0.0).unwrap();
        same(out.beam.power(), Stokes { i: 2.0, q: 0.0, u: -1.0, v: -1.0 });
        assert_eq!(out.beam.frame().horizontal(), b.frame().horizontal());
    }

    #[test]
    fn g1_transmission_switches_material_and_uses_power_not_radiance() {
        let b = beam(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0),
            glass(1.0, 0.0), Stokes::try_new(1.0, 0.3, 0.2, 0.1).unwrap());
        let target = glass(1.5, 2.0);
        let t = b.sample_smooth_interface(Vec3::new(0.0, 0.0, 1.0), target, 0.99).unwrap();
        same(t.beam.power(), b.power());
        assert_eq!(t.beam.medium(), target);
        let moved = t.beam.advance(0.25).unwrap();
        near(moved.power().i, det::exp(-0.5));
        near(moved.power().degree(), b.power().degree());
        near(moved.ray().origin.z, -0.25);
        let r = b.sample_smooth_interface(Vec3::new(0.0, 0.0, 1.0), target, 0.0).unwrap();
        near(r.beam.advance(0.25).unwrap().power().i, 1.0);
        near(r.pdf * r.beam.power().i + t.pdf * t.beam.power().i, b.power().i);
    }

    #[test]
    fn g0_invalid_operations_do_not_mutate_the_original_beam() {
        let b = beam(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0),
            glass(1.0, 0.0), Stokes::unpolarized(1.0));
        let original = b;
        assert!(b.advance(-1.0).is_err());
        assert!(b.advance(f64::INFINITY).is_err());
        assert!(b.sample_smooth_interface(Vec3::new(1.0, 0.0, 0.0), glass(1.5, 0.0), 0.5).is_err());
        assert!(b.sample_smooth_interface(Vec3::new(0.0, 0.0, 1.0), glass(1.5, 0.0), 1.0).is_err());
        assert!(PolarizedBeam::try_new(b.ray().origin, b.frame(), 10.0, b.medium(), b.power()).is_err());
        assert_eq!(b, original);
    }
}
