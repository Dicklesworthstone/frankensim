//! Forward polarized ray-power transport through native dielectric materials.
//!
//! A ray carries its transverse frame, vacuum wavelength, Stokes power and
//! current medium together. Intersections still belong to the existing chart/
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
//!
//! Spatially varying media use the native adaptive eikonal propagator, including
//! parallel transport of the transverse frame. The field is supplied explicitly
//! at each advance; a sampled medium cannot silently take a straight homogeneous
//! step. Boundary scattering uses its local phase index, never a nominal glass
//! value or the group index. Curved-ray intersection/event location remains the
//! caller's responsibility: a propagation interval must stay within one smooth
//! field domain.

use fs_geom::{Point3, Vec3};
use fs_math::det;

use crate::charts::Ray;
use crate::dielectric::{
    DielectricError, DielectricEvent, DielectricGlass, DielectricSurface, sample_rough_dielectric,
};
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexPropagation,
    GradientIndexSample, propagate_with_control,
};
use crate::polarization::{Mueller, Stokes};

use super::{
    attenuate_stokes, fresnel_dielectric_polarized, sample_smooth_dielectric_polarized,
};

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
    /// A spatially varying field or its bounded integration refused the request.
    GradientIndex(GradientIndexError),
    /// A gradient medium requires curved propagation with its original field.
    GradientFieldRequired,
    /// The supplied field, position or wavelength differs from the retained local state.
    MediumMismatch,
    /// A finite physical result cannot be represented in binary64.
    NumericalRange,
}

impl From<DielectricError> for BeamError {
    fn from(error: DielectricError) -> Self {
        Self::Dielectric(error)
    }
}

impl From<GradientIndexError> for BeamError {
    fn from(error: GradientIndexError) -> Self {
        Self::GradientIndex(error)
    }
}

impl core::fmt::Display for BeamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized beam: {self:?}")
    }
}

impl core::error::Error for BeamError {}

/// A field sample tied to the exact point and vacuum wavelength where it was
/// evaluated. Construct through [`BeamMedium::try_gradient_at`]; it is not a
/// homogeneous material and cannot be used to advance a straight ray.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalGradientMedium {
    position: Point3,
    wavelength_nm: f64,
    sample: GradientIndexSample,
}

impl LocalGradientMedium {
    /// Local optical properties, including the phase gradient and true absorption.
    #[must_use]
    pub const fn sample(self) -> GradientIndexSample {
        self.sample
    }

    /// The point in metres at which this state applies.
    #[must_use]
    pub const fn position(self) -> Point3 {
        self.position
    }

    /// Vacuum wavelength in nanometres at which this state applies.
    #[must_use]
    pub const fn wavelength_nm(self) -> f64 {
        self.wavelength_nm
    }
}

/// The physical medium at a polarized beam's current position.
///
/// A homogeneous law applies throughout a segment. A gradient field sample is
/// local: subsequent propagation requires the same immutable physical field.
/// The owner of the accepted simulation state retains the field and its source
/// identity. Equality at one point does not prove that two whole fields agree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BeamMedium {
    /// Native spectral dispersion and absorption in a homogeneous dielectric.
    Homogeneous(DielectricGlass),
    /// Validated optical data at exactly one position and wavelength.
    GradientIndex(LocalGradientMedium),
}

impl From<DielectricGlass> for BeamMedium {
    fn from(material: DielectricGlass) -> Self {
        Self::Homogeneous(material)
    }
}

impl BeamMedium {
    /// Sample the incident or target medium for a beam launch or interface.
    /// This does not perform a physical boundary crossing: use the existing
    /// smooth/rough sampling operation to apply Snell/Fresnel at an interface.
    pub fn try_gradient_at<F: GradientIndexField + ?Sized>(
        field: &F,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<Self, BeamError> {
        if !finite_point(position) {
            return Err(BeamError::InvalidPosition);
        }
        if !(crate::spectral::LAMBDA_MIN..=crate::spectral::LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(DielectricError::InvalidWavelength.into());
        }
        let sample = field.sample(position, wavelength_nm)?;
        Ok(Self::GradientIndex(LocalGradientMedium {
            position,
            wavelength_nm,
            sample,
        }))
    }

    fn phase_index_at(self, position: Point3, wavelength_nm: f64) -> Result<f64, BeamError> {
        match self {
            Self::Homogeneous(material) => Ok(material.ior().eval(wavelength_nm)?),
            Self::GradientIndex(local) => {
                if position != local.position || wavelength_nm != local.wavelength_nm {
                    return Err(BeamError::MediumMismatch);
                }
                Ok(local.sample.phase_index())
            }
        }
    }
}

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
    medium: BeamMedium,
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

/// Curved polarized continuation and the actual integrated optical quantities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BeamPropagation {
    /// Endpoint position, transported frame, attenuated power and local medium.
    pub beam: PolarizedBeam,
    /// Native error-controlled propagation result and work counters.
    pub propagation: GradientIndexPropagation,
}

impl PolarizedBeam {
    /// Admit a beam without inferring its surrounding medium or polarization.
    pub fn try_new(
        origin: Point3,
        frame: PolarizationFrame,
        wavelength_nm: f64,
        medium: impl Into<BeamMedium>,
        power: Stokes,
    ) -> Result<Self, BeamError> {
        if !finite_point(origin) { return Err(BeamError::InvalidPosition); }
        let medium = medium.into();
        medium.phase_index_at(origin, wavelength_nm)?;
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
    pub const fn medium(self) -> BeamMedium {
        self.medium
    }

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
        let BeamMedium::Homogeneous(material) = self.medium else {
            return Err(BeamError::GradientFieldRequired);
        };
        let origin = self.ray().at(distance_m);
        if !finite_point(origin) { return Err(BeamError::InvalidPosition); }
        let power = attenuate_stokes(
            self.power, material.absorption(), self.wavelength_nm, distance_m,
        )?;
        Ok(Self { origin, power, ..self })
    }

    /// Advance through the same smooth, isotropic field used to admit this
    /// beam. The Stokes state is expressed in the frame parallel-transported
    /// along the entire curved path; endpoint projection is not a substitute.
    /// All four power components receive the integrated Beer--Lambert factor.
    ///
    /// The local index changes geometrical focusing, not the power of an
    /// individual forward ray. No radiance eta-squared factor is applied here.
    /// This does not integrate scattering or locate a material discontinuity.
    pub fn advance_gradient<F: GradientIndexField + ?Sized>(
        self,
        field: &F,
        distance_m: f64,
        options: GradientIndexOptions,
    ) -> Result<BeamPropagation, BeamError> {
        self.advance_gradient_with_control(field, distance_m, options, || true)
    }

    /// Curved advance with a bounded-work cancellation check before every
    /// attempted adaptive step, including rejected steps. Returning false
    /// refuses the whole operation and leaves the original beam available.
    pub fn advance_gradient_with_control<F: GradientIndexField + ?Sized>(
        self,
        field: &F,
        distance_m: f64,
        options: GradientIndexOptions,
        control: impl FnMut() -> bool,
    ) -> Result<BeamPropagation, BeamError> {
        let BeamMedium::GradientIndex(local) = self.medium else {
            return Err(BeamError::MediumMismatch);
        };
        if field.sample(self.origin, self.wavelength_nm)? != local.sample {
            return Err(BeamError::MediumMismatch);
        }
        let propagation = propagate_with_control(
            field,
            self.ray(),
            self.frame.horizontal,
            self.wavelength_nm,
            distance_m,
            options,
            control,
        )?;
        self.consume_propagation(self.ray(), self.wavelength_nm, propagation)
    }

    /// Internal continuation for the event owner, which already solved this
    /// segment. Keep source coordinates/wavelength and local medium checks at
    /// this boundary; callers cannot publish arbitrary propagation records.
    pub(super) fn consume_propagation(
        self,
        initial_ray: Ray,
        wavelength_nm: f64,
        propagation: GradientIndexPropagation,
    ) -> Result<BeamPropagation, BeamError> {
        if initial_ray != self.ray() || wavelength_nm != self.wavelength_nm {
            return Err(BeamError::MediumMismatch);
        }
        let medium = match self.medium {
            BeamMedium::Homogeneous(material) => {
                let field = crate::gradient_index::HomogeneousIndexField::new(&material);
                if field.sample(self.origin, wavelength_nm)? != propagation.initial_sample
                    || propagation.final_sample != propagation.initial_sample
                {
                    return Err(BeamError::MediumMismatch);
                }
                self.medium
            }
            BeamMedium::GradientIndex(local) => {
                if local.position != self.origin
                    || local.wavelength_nm != wavelength_nm
                    || local.sample != propagation.initial_sample
                {
                    return Err(BeamError::MediumMismatch);
                }
                BeamMedium::GradientIndex(LocalGradientMedium {
                    position: propagation.ray.origin,
                    wavelength_nm,
                    sample: propagation.final_sample,
                })
            }
        };
        let frame = PolarizationFrame::try_new(propagation.ray.dir, propagation.transverse_axis)?;
        let transmission = propagation.power_transmittance();
        let power = physical_roundoff(
            Stokes {
                i: self.power.i * transmission,
                q: self.power.q * transmission,
                u: self.power.u * transmission,
                v: self.power.v * transmission,
            },
            self.power.i,
        )?;
        let origin = propagation.ray.origin;
        Ok(BeamPropagation {
            beam: Self {
                origin,
                frame,
                power,
                medium,
                ..self
            },
            propagation,
        })
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
        target: impl Into<BeamMedium>,
        event_sample: f64,
    ) -> Result<BeamSample, BeamError> {
        let normal = self.incident_normal(normal)?;
        let target = target.into();
        let eta_i = self
            .medium
            .phase_index_at(self.origin, self.wavelength_nm)?;
        let eta_t = target.phase_index_at(self.origin, self.wavelength_nm)?;
        let sample = sample_smooth_dielectric_polarized(
            normal, self.frame.direction.scale(-1.0), eta_i, eta_t, event_sample,
        )?;
        self.scatter_smooth_sample(normal, target, eta_i, eta_t, sample)
    }

    /// Sample one original transparent coating on a smooth dielectric boundary.
    /// Reflection/transmission use this beam's actual spectral/local index and
    /// the supplied target medium, with coherent s/p phase and normal power.
    /// Only transmission changes the medium. Film thickness is optical; both
    /// local film faces are collapsed onto the current geometric hit.
    pub fn sample_coated_interface(
        self,
        normal: Vec3,
        target: impl Into<BeamMedium>,
        coating: &crate::conductor::CoherentDielectricCoating,
        event_sample: f64,
    ) -> Result<BeamSample, BeamError> {
        let normal = self.incident_normal(normal)?;
        let target = target.into();
        let eta_i = self.medium.phase_index_at(self.origin, self.wavelength_nm)?;
        let eta_t = target.phase_index_at(self.origin, self.wavelength_nm)?;
        let sample = super::sample_smooth_coated_dielectric_polarized(
            normal, self.frame.direction.scale(-1.0), eta_i, eta_t,
            self.wavelength_nm, coating, event_sample,
        )?;
        self.scatter_smooth_sample(normal, target, eta_i, eta_t, sample)
    }

    fn scatter_smooth_sample(
        self,
        normal: Vec3,
        target: BeamMedium,
        eta_i: f64,
        eta_t: f64,
        sample: super::PolarizedDielectricSample,
    ) -> Result<BeamSample, BeamError> {
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

    /// Sample the native smooth or isotropic-GGX surface without losing its
    /// polarization frame. `samples` is `[microfacet_u, microfacet_v, event]`,
    /// each in `[0,1)`. Smooth surfaces ignore the first two admitted samples.
    ///
    /// The original sampler owns directions, visibility/masking, branch
    /// selection and density. Its scalar weight multiplies the normalized
    /// Fresnel Mueller response of the reconstructed microfacet. Thus an
    /// unpolarized input reproduces the native forward-power sample, while
    /// subsequent noncoplanar interactions retain linear/circular polarization.
    ///
    /// `None` is the native zero-contribution microfacet rejection, not an
    /// error. This keeps the native single-scattering microfacet approximation;
    /// unresolved inter-facet multiple scattering is not synthesized here.
    pub fn sample_interface(
        self,
        normal: Vec3,
        target: impl Into<BeamMedium>,
        surface: DielectricSurface,
        samples: [f64; 3],
    ) -> Result<Option<BeamSample>, BeamError> {
        if samples.iter().any(|u| !u.is_finite() || !(0.0..1.0).contains(u)) {
            return Err(DielectricError::InvalidInterface.into());
        }
        let target = target.into();
        let Some(alpha) = surface.roughness_alpha() else {
            return self.sample_smooth_interface(normal, target, samples[2]).map(Some);
        };
        let normal = self.incident_normal(normal)?;
        let wo = self.frame.direction.scale(-1.0);
        let eta_i = self
            .medium
            .phase_index_at(self.origin, self.wavelength_nm)?;
        let eta_t = target.phase_index_at(self.origin, self.wavelength_nm)?;
        let Some(sample) = sample_rough_dielectric(
            normal, wo, eta_i, eta_t, alpha, samples[0], samples[1], samples[2],
        )? else { return Ok(None); };
        if sample.delta {
            // Equal phase index is a geometrically null interface, but not
            // necessarily the same absorption law. Preserve the exact frame
            // and power while installing the material on the other side.
            return Ok(Some(BeamSample {
                beam: Self { medium: target, ..self },
                event: sample.event, pdf: sample.pdf, delta: true,
            }));
        }
        // Match the native evaluator's normalized half-vector convention.
        let wo = wo.scale(1.0 / wo.norm());
        let wi = sample.direction.scale(1.0 / sample.direction.norm());
        let ratio = if sample.event == DielectricEvent::Transmission { eta_t / eta_i } else { 1.0 };
        let half = Vec3::new(wo.x + ratio * wi.x, wo.y + ratio * wi.y, wo.z + ratio * wi.z);
        let mut micro_normal = unit(half)?;
        if normal.dot(micro_normal) < 0.0 { micro_normal = micro_normal.scale(-1.0); }
        let fresnel = fresnel_dielectric_polarized(wo.dot(micro_normal).clamp(0.0, 1.0), eta_i, eta_t)?;
        let matrix = if sample.event == DielectricEvent::Reflection {
            fresnel.reflection
        } else { fresnel.transmission };
        let mean = matrix.0[0][0];
        if !mean.is_finite() || mean <= 0.0 { return Err(BeamError::NumericalRange); }
        let weight = sample.radiance_weight * ratio * ratio;
        if !weight.is_finite() { return Err(BeamError::NumericalRange); }
        // Normalize before multiplying, avoiding weight/mean overflow for
        // extremely weak reflectors. Each normalized entry is bounded by two.
        let matrix = Mueller(matrix.0.map(|row| row.map(|v| (v / mean) * weight)));
        let beam = self.scatter(micro_normal, sample.direction, matrix, sample.event, target)?;
        Ok(Some(BeamSample { beam, event: sample.event, pdf: sample.pdf, delta: false }))
    }

    fn incident_normal(self, normal: Vec3) -> Result<Vec3, BeamError> {
        let normal = unit(normal)?;
        Ok(if normal.dot(self.frame.direction) > 0.0 { normal.scale(-1.0) } else { normal })
    }

    fn scatter(
        self,
        normal: Vec3,
        direction: Vec3,
        matrix: Mueller,
        event: DielectricEvent,
        target: BeamMedium,
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

#[cfg(feature = "photometry")]
impl crate::polarization::photometry::PolarizedSource {
    /// Emit a monochromatic polarized ray-power weight from the calibrated
    /// fixture's native position, orientation and angular intensity profile.
    ///
    /// `bin` selects the existing midpoint wavelength. The result is W/nm:
    /// surviving spectral intensity [W/(sr nm)] times the explicit angular
    /// quadrature weight [sr], or `1/pdf` for a Monte Carlo direction sample.
    /// Average Monte Carlo samples and integrate wavelength bins separately.
    /// No implicit PDF, wavelength bandwidth or replacement flux is inferred.
    ///
    /// The local +X fixture axis projected transverse to the outgoing ray
    /// defines positive Q. This is an explicit off-axis extension of the source
    /// optics' paraxial reference; directions parallel to +X/-X are refused.
    /// The caller supplies the actual surrounding medium and owns intersections
    /// and nested-medium topology after emission. Scalar camera paths remain
    /// unchanged; this method feeds the forward polarized beam path.
    pub fn emit_beam(
        &self,
        world_direction: Vec3,
        bin: usize,
        solid_angle_weight_sr: f64,
        medium: impl Into<BeamMedium>,
    ) -> Result<PolarizedBeam, crate::photometry::source::SourceError> {
        use crate::photometry::source::{SPECTRAL_BINS, SourceError, wavelength_nm};
        if bin >= SPECTRAL_BINS { return Err(SourceError::InvalidSpectrum); }
        if !solid_angle_weight_sr.is_finite() || solid_angle_weight_sr < 0.0 {
            return Err(SourceError::InvalidFlux);
        }
        let source = self.intensity_source();
        let axis = source.transform().transform_vector(Vec3::new(1.0, 0.0, 0.0));
        let frame = PolarizationFrame::try_new(world_direction, axis)
            .map_err(|_| SourceError::InvalidDirection)?;
        let wavelength = wavelength_nm(bin);
        let intensity = source.spectral_intensity_w_sr_nm(frame.direction(), wavelength)?
            * solid_angle_weight_sr;
        if !intensity.is_finite() { return Err(SourceError::NumericalRange); }
        let state = self.relative_stokes()[bin];
        // The native intensity source ALREADY includes all optical-chain loss.
        // Restore only normalized Q/U/V, never multiply by state.i again.
        // A dark bin needs neither division nor an invented polarization axis.
        let power = if state.i == 0.0 || intensity == 0.0 { Stokes::unpolarized(0.0) } else {
            Stokes { i: intensity, q: intensity * (state.q / state.i),
                u: intensity * (state.u / state.i), v: intensity * (state.v / state.i) }
        };
        PolarizedBeam::try_new(source.position_m(), frame, wavelength, medium, power)
            .map_err(|_| SourceError::NumericalRange)
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
        assert_eq!(t.beam.medium(), BeamMedium::Homogeneous(target));
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

    #[test]
    fn g1_rough_glass_keeps_native_proposals_directions_and_unpolarized_power() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let mut contributing = [0; 2];
        let mut rejected = 0;
        for (eta_i, eta_t) in [(1.0, 1.5), (1.5, 1.0)] {
            for cosine in [0.4_f64, 0.8, 1.0] {
                let b = beam(Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, -cosine),
                    Vec3::new(0.0, 1.0, 0.0), glass(eta_i, 0.0), Stokes::unpolarized(1.0));
                for alpha in [0.05, 0.4, 1.0] {
                    let surface = DielectricSurface::try_rough(alpha).unwrap();
                    for samples in [[0.0, 0.0, 0.0], [0.17, 0.63, 0.999], [0.9, 0.2, 0.0]] {
                        let native = sample_rough_dielectric(normal, b.frame().direction().scale(-1.0),
                            eta_i, eta_t, alpha, samples[0], samples[1], samples[2]).unwrap();
                        let actual = b.sample_interface(normal, glass(eta_t, 0.0), surface, samples).unwrap();
                        match (native, actual) {
                            (None, None) => rejected += 1,
                            (Some(native), Some(actual)) => {
                                assert_eq!(actual.event, native.event);
                                assert_eq!(actual.pdf, native.pdf);
                                assert_eq!(actual.delta, native.delta);
                                near(actual.beam.ray().dir.x, native.direction.x);
                                near(actual.beam.ray().dir.y, native.direction.y);
                                near(actual.beam.ray().dir.z, native.direction.z);
                                let transmitted = native.event == DielectricEvent::Transmission;
                                contributing[usize::from(transmitted)] += 1;
                                let ratio = if transmitted { eta_t / eta_i } else { 1.0 };
                                near(actual.beam.power().i, native.radiance_weight * ratio * ratio);
                            }
                            other => panic!("native and polarized support differ: {other:?}"),
                        }
                    }
                }
            }
        }
        assert!(contributing.iter().all(|n| *n > 0));
        assert!(rejected > 0);
    }

    #[test]
    fn g3_rough_scattering_is_invariant_under_incoming_reference_changes() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let b = beam(Vec3::new(0.6, 0.0, -0.8), Vec3::new(0.0, 1.0, 0.0),
            glass(1.0, 0.0), Stokes::try_new(1.0, 0.3, -0.4, 0.2).unwrap());
        let reframed = b.with_reference_axis(Vec3::new(1.0, 2.0, 3.0)).unwrap();
        for event in [0.0, 0.999] {
            let samples = [0.17, 0.63, event];
            let a = b.sample_interface(normal, glass(1.5, 0.0), DielectricSurface::POLISHED, samples)
                .unwrap().unwrap();
            let c = reframed.sample_interface(normal, glass(1.5, 0.0), DielectricSurface::POLISHED, samples)
                .unwrap().unwrap();
            same(a.beam.power(), c.beam.power());
            assert_eq!(a.pdf, c.pdf);
            let s = a.beam.power();
            assert!(Stokes::try_new(s.i, s.q, s.u, s.v).is_some());
        }
    }

    #[test]
    fn g0_null_rough_boundary_changes_absorption_without_changing_frame_or_power() {
        let b = beam(Vec3::new(0.6, 0.0, -0.8), Vec3::new(0.0, 1.0, 0.0),
            glass(1.5, 0.0), Stokes::try_new(2.0, 0.4, 0.5, -0.7).unwrap());
        let target = glass(1.5, 3.0);
        let out = b.sample_interface(Vec3::new(0.0, 0.0, 1.0), target,
            DielectricSurface::POLISHED, [0.3, 0.7, 0.5]).unwrap().unwrap();
        assert!(out.delta);
        assert_eq!(out.pdf, 1.0);
        assert_eq!(out.beam.frame(), b.frame());
        assert_eq!(out.beam.power(), b.power());
        assert_eq!(out.beam.medium(), BeamMedium::Homogeneous(target));
        near(out.beam.advance(0.2).unwrap().power().i, 2.0 * det::exp(-0.6));
        assert!(b.sample_interface(Vec3::new(0.0, 0.0, 1.0), target,
            DielectricSurface::SMOOTH, [f64::NAN, 0.0, 0.5]).is_err());
    }

    struct BendingField {
        absorption_per_m: f64,
    }

    impl GradientIndexField for BendingField {
        fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            let n = 1.4 + 0.3 * p.y;
            GradientIndexSample::try_new(
                n,
                Vec3::new(0.0, 0.3, 0.0),
                self.absorption_per_m,
                Some(n),
            )
        }
    }

    fn gradient_beam(field: &impl GradientIndexField) -> PolarizedBeam {
        let origin = Point3::new(0.0, 0.0, 0.0);
        PolarizedBeam::try_new(
            origin,
            PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 1.0)).unwrap(),
            550.0,
            BeamMedium::try_gradient_at(field, origin, 550.0).unwrap(),
            Stokes::try_new(2.0, 0.6, -0.8, 0.4).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn g1_curved_beam_transports_polarization_and_uses_local_snell_index() {
        let field = BendingField {
            absorption_per_m: 0.7,
        };
        let initial = gradient_beam(&field);
        let result = initial
            .advance_gradient(&field, 1.0, GradientIndexOptions::default())
            .unwrap();
        let moved = result.beam;
        let k = moved.frame().direction();
        // In an affine y-gradient, n t_x is invariant and n t_y = g s.
        let n = (1.4_f64 * 1.4 + 0.3 * 0.3).sqrt();
        assert!((k.x - 1.4 / n).abs() < 2e-10);
        assert!((k.y - 0.3 / n).abs() < 2e-10);
        let h = moved.frame().horizontal();
        let inv_sqrt_two = core::f64::consts::FRAC_1_SQRT_2;
        // The initial 45-degree axis combines a rotating in-plane vector and
        // a fixed out-of-plane vector. Projecting only at the endpoint gives
        // a different (incorrect) relative angle between these components.
        assert!((h.x + k.y * inv_sqrt_two).abs() < 2e-10);
        assert!((h.y - k.x * inv_sqrt_two).abs() < 2e-10);
        assert!((h.z - inv_sqrt_two).abs() < 2e-10);
        let projection = PolarizationFrame::try_new(k, initial.frame().horizontal()).unwrap();
        assert!((projection.horizontal().z - h.z).abs() > 0.001);
        let factor = det::exp(-0.7);
        same(
            moved.power(),
            Stokes {
                i: 2.0 * factor,
                q: 0.6 * factor,
                u: -0.8 * factor,
                v: 0.4 * factor,
            },
        );
        near(moved.power().degree(), initial.power().degree());
        assert!(result.propagation.optical_path_length_m > 1.4);

        let normal = Vec3::new(-1.0, 0.0, 0.0);
        let exit = moved
            .sample_smooth_interface(normal, glass(1.0, 0.0), 0.999)
            .unwrap();
        assert_eq!(exit.event, DielectricEvent::Transmission);
        assert!((exit.beam.ray().dir.y - 0.3).abs() < 3e-10);
        assert!((exit.beam.ray().dir.x - 0.91_f64.sqrt()).abs() < 3e-10);
        let after_air = exit.beam.advance(0.5).unwrap();
        same(after_air.power(), exit.beam.power());
        let reflection = moved
            .sample_smooth_interface(normal, glass(1.0, 0.0), 0.0)
            .unwrap();
        assert_eq!(reflection.beam.medium(), moved.medium());
        assert_eq!(
            reflection.beam.advance(0.1),
            Err(BeamError::GradientFieldRequired)
        );
    }

    #[test]
    fn g3_curved_beam_segment_composition_preserves_the_physical_frame() {
        let field = BendingField {
            absorption_per_m: 0.2,
        };
        let initial = gradient_beam(&field);
        let options = GradientIndexOptions::default();
        let whole = initial.advance_gradient(&field, 1.0, options).unwrap();
        let first = initial.advance_gradient(&field, 0.4, options).unwrap();
        let split = first.beam.advance_gradient(&field, 0.6, options).unwrap();
        assert!(
            whole
                .beam
                .ray()
                .origin
                .delta_from(split.beam.ray().origin)
                .norm()
                < 2e-10
        );
        for (a, b) in [
            (
                whole.beam.frame().horizontal(),
                split.beam.frame().horizontal(),
            ),
            (
                whole.beam.frame().direction(),
                split.beam.frame().direction(),
            ),
        ] {
            assert!((a.x - b.x).abs() < 2e-10);
            assert!((a.y - b.y).abs() < 2e-10);
            assert!((a.z - b.z).abs() < 2e-10);
        }
        same(whole.beam.power(), split.beam.power());
        assert!(
            (whole.propagation.optical_path_length_m
                - first.propagation.optical_path_length_m
                - split.propagation.optical_path_length_m)
                .abs()
                < 2e-10
        );
    }

    #[test]
    fn g0_gradient_beam_requires_matching_local_state_and_transactional_control() {
        let field = BendingField {
            absorption_per_m: 0.0,
        };
        let initial = gradient_beam(&field);
        let original = initial;
        let options = GradientIndexOptions::default();
        assert_eq!(initial.advance(0.1), Err(BeamError::GradientFieldRequired));
        let changed = BendingField {
            absorption_per_m: 0.5,
        };
        assert_eq!(
            initial.advance_gradient(&changed, 0.5, options),
            Err(BeamError::MediumMismatch)
        );
        assert_eq!(
            initial.advance_gradient_with_control(&field, 1.0, options, || false),
            Err(BeamError::GradientIndex(GradientIndexError::Cancelled))
        );
        let displaced =
            BeamMedium::try_gradient_at(&field, Point3::new(0.0, 0.1, 0.0), 550.0).unwrap();
        assert_eq!(
            initial.sample_smooth_interface(Vec3::new(-1.0, 0.0, 0.0), displaced, 0.9),
            Err(BeamError::MediumMismatch)
        );
        assert_eq!(initial, original);
    }

    #[test]
    fn g1_native_surface_sampling_enters_a_gradient_medium_without_power_rescaling() {
        let field = BendingField {
            absorption_per_m: 0.0,
        };
        let origin = Point3::new(0.0, 0.0, 0.0);
        let initial = beam(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            glass(1.0, 0.0),
            Stokes::unpolarized(1.0),
        );
        let target = BeamMedium::try_gradient_at(&field, origin, 550.0).unwrap();
        let normal = Vec3::new(-1.0, 0.0, 0.0);
        for surface in [DielectricSurface::SMOOTH, DielectricSurface::POLISHED] {
            let samples = [0.0, 0.0, 0.999];
            let actual = initial
                .sample_interface(normal, target, surface, samples)
                .unwrap()
                .unwrap();
            let homogeneous = initial
                .sample_interface(normal, glass(1.4, 0.0), surface, samples)
                .unwrap()
                .unwrap();
            assert_eq!(actual.event, DielectricEvent::Transmission);
            assert_eq!(actual.pdf, homogeneous.pdf);
            assert_eq!(actual.beam.ray(), homogeneous.beam.ray());
            same(actual.beam.power(), homogeneous.beam.power());
            assert_eq!(actual.beam.medium(), target);
            let moved = actual
                .beam
                .advance_gradient(&field, 0.5, GradientIndexOptions::default())
                .unwrap();
            assert!(moved.beam.ray().dir.y > 0.0);
            same(moved.beam.power(), actual.beam.power());
        }
    }

    #[cfg(feature = "photometry")]
    fn source() -> crate::photometry::source::PhotometricSource {
        use crate::instances::RigidTransform;
        use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum, SPECTRAL_BINS};
        PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0; 3]).unwrap(),
            1000.0, 1.0, AngularDistribution::isotropic(),
            SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
        ).unwrap()
    }

    #[test]
    #[cfg(feature = "photometry")]
    fn g1_source_emission_keeps_circular_power_and_filters_exactly_once() {
        use crate::photometry::source::{SPECTRAL_BINS, wavelength_nm};
        use crate::polarization::photometry::SpectralOptics;
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let original = source();
        let circular = original.clone().with_unpolarized_optics()
            .through(&SpectralOptics::linear_polarizer(core::f64::consts::FRAC_PI_4).unwrap()).unwrap()
            .through(&SpectralOptics::linear_retarder(0.0,
                [core::f64::consts::FRAC_PI_2; SPECTRAL_BINS]).unwrap()).unwrap();
        let expected = original.spectral_intensity_w_sr_nm(direction, wavelength_nm(30)).unwrap()
            * 0.125 * 0.5;
        let emitted = circular.emit_beam(direction, 30, 0.125, glass(1.5, 2.0)).unwrap();
        near(emitted.power().i, expected);
        near(emitted.power().v, -expected);
        near(emitted.power().degree(), 1.0);
        assert_eq!(emitted.ray().origin, original.position_m());
        near(emitted.advance(0.3).unwrap().power().i, expected * det::exp(-0.6));
        let dark = circular.through(&SpectralOptics::neutral([0.0; SPECTRAL_BINS]).unwrap()).unwrap();
        assert_eq!(dark.emit_beam(direction, 30, 0.125, glass(1.0, 0.0)).unwrap().power(),
            Stokes::unpolarized(0.0));
    }

    #[test]
    #[cfg(feature = "photometry")]
    fn g3_emission_follows_fixture_rotation_and_refuses_undefined_axes() {
        use crate::instances::RigidTransform;
        use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum, SPECTRAL_BINS};
        use crate::polarization::photometry::SpectralOptics;
        let h = 0.5_f64.sqrt();
        let polarized = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, h, h], [1.0, 2.0, 3.0]).unwrap(),
            1000.0, 1.0, AngularDistribution::isotropic(),
            SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
        ).unwrap().with_unpolarized_optics()
            .through(&SpectralOptics::linear_polarizer(0.0).unwrap()).unwrap();
        let out = polarized.emit_beam(Vec3::new(0.0, 0.0, 1.0), 30, 0.2, glass(1.0, 0.0)).unwrap();
        near(out.frame().horizontal().x, 0.0); near(out.frame().horizontal().y, 1.0);
        near(out.power().q, out.power().i);
        assert_eq!(out.ray().origin, Point3::new(1.0, 2.0, 3.0));
        let original = source().with_unpolarized_optics();
        assert!(original.emit_beam(Vec3::new(1.0, 0.0, 0.0), 30, 0.2, glass(1.0, 0.0)).is_err());
        assert!(original.emit_beam(Vec3::new(0.0, 0.0, 1.0), SPECTRAL_BINS, 0.2, glass(1.0, 0.0)).is_err());
        assert!(original.emit_beam(Vec3::new(0.0, 0.0, 1.0), 30, -0.2, glass(1.0, 0.0)).is_err());
    }

    #[test]
    #[cfg(feature = "photometry")]
    fn g1_calibrated_source_through_two_spatial_brewster_reflections() {
        let b = source().with_unpolarized_optics()
            .emit_beam(Vec3::new(0.0, 0.0, -1.0), 30, 0.2, glass(1.0, 0.0)).unwrap();
        let ci = 1.0 / 1.5_f64.hypot(1.0);
        let si = 1.5 * ci;
        let first = b.advance(0.4).unwrap().sample_interface(Vec3::new(si, 0.0, ci),
            glass(1.5, 0.0), DielectricSurface::SMOOTH, [0.0; 3]).unwrap().unwrap();
        assert!(first.beam.power().i > 0.0);
        let k = first.beam.frame().direction();
        let x = first.beam.frame().horizontal();
        let n = Vec3::new(-ci*k.x + si*x.x, -ci*k.y + si*x.y, -ci*k.z + si*x.z);
        let second = first.beam.advance(0.3).unwrap().sample_interface(n,
            glass(1.5, 0.0), DielectricSurface::SMOOTH, [0.0; 3]).unwrap().unwrap();
        near(second.beam.power().i, 0.0);
    }
}
