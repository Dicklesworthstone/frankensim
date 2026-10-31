//! Scalar incoherent imaging through uniformly illuminated camera pupils.
//!
//! In the focused image plane, `alpha = pi*D/(lambda*q)` and the unit-integral
//! PSF is `alpha^2/(4*pi) * [2*J1(alpha*r)/(alpha*r)]^2`. Here `q` is image
//! distance, not focal length: at finite conjugates `1/f = 1/s + 1/q`.
//! In this camera's angular raster coordinates `r = q*theta`, so `q` cancels
//! and the scale is `lambda/(pi*D)`. No focal-length approximation is needed.
//! The normalization follows from `integral_0^inf J1(x)^2/x dx = 1/2`.
//!
//! References: TU Delft, *Scalar Diffraction Optics*, “Imaging by a lens”
//! (Airy pattern and finite image coordinates), and *Advanced Optical Imaging*,
//! “Incoherent imaging” (intensity convolution).
//! <https://interactivetextbooks.tudelft.nl/interactive-optics/content/Chap6_Diffraction/DiffractiveOptics_2022_01Clean.html>
//! <https://qiweb.tudelft.nl/aoi/incoherentimaging/incoherentimaging/>
//!
//! This is not a defocused-pupil model. A declared upper bound on the
//! quadratic pupil phase, `pi*a^2/lambda * |1/z - 1/s|`, gates visible depths.
//! Geometric pupil blur is not added to the Airy pattern. Coherence between
//! scene points, polarization, aberrations and pupil-dependent occlusion are
//! outside this scalar, paraxial, isoplanatic model.

use super::{ApertureKind, CameraError, LensSample, PhysicalCamera, direct_unit};
use crate::charts::Ray;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::{det, special::bessel_j1};

mod polygon;
use polygon::PolygonPupil;

const PI: f64 = core::f64::consts::PI;
const PROPOSAL_RADIUS: f64 = 2.0;

/// Explicit in-focus scalar diffraction for circular and regular bladed pupils.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FocusedDiffraction {
    maximum_defocus_phase_rad: f64,
}

impl FocusedDiffraction {
    /// Admit a maximum pupil-rim defocus phase in `[0, 0.1]` radians.
    /// Zero requires exact focus. This is a phase budget, not a certified image
    /// error: omitted phase can also interact with occlusion or aberrations.
    pub fn try_new(maximum_defocus_phase_rad: f64) -> Result<Self, CameraError> {
        if !maximum_defocus_phase_rad.is_finite()
            || !(0.0..=0.1).contains(&maximum_defocus_phase_rad)
        {
            return Err(CameraError::InvalidDiffraction);
        }
        Ok(Self {
            maximum_defocus_phase_rad,
        })
    }

    /// Admitted pupil-rim phase budget in radians.
    #[must_use]
    pub const fn maximum_defocus_phase_rad(self) -> f64 {
        self.maximum_defocus_phase_rad
    }
}

/// One shared geometric ray with wavelength-specific pupil/proposal weights.
/// These weights multiply incoming spectral radiance before color integration.
#[derive(Clone, Copy, Debug)]
pub struct SpectralCameraRay<const N: usize> {
    /// Incoming ray through the optical centre for focused diffraction, or
    /// through the sampled geometric pupil when diffraction is not selected.
    pub ray: Ray,
    /// Importance weights, in the same order as the supplied wavelengths.
    pub weights: [f64; N],
}

impl PhysicalCamera {
    /// Select in-focus diffraction for a nonzero circular or bladed pupil. The rendered
    /// raster retains the camera's declared angular field of view. The pupil
    /// sets optical resolution; it is not also sampled for geometric defocus.
    /// A declared physical focal length requires focus beyond that length,
    /// so the conjugate image plane lies at a positive distance.
    /// Existing geometric-camera behavior is unchanged until this is called.
    pub fn with_in_focus_diffraction(
        mut self,
        policy: FocusedDiffraction,
    ) -> Result<Self, CameraError> {
        if self.aperture.is_pinhole() {
            return Err(CameraError::InvalidDiffraction);
        }
        self.admit_diffraction_projection()?;
        self.diffraction = Some(policy);
        Ok(self)
    }

    /// Selected focused diffraction policy, if any.
    #[must_use]
    pub const fn in_focus_diffraction(&self) -> Option<FocusedDiffraction> {
        self.diffraction
    }

    /// Trace one sample of the wavelength-dependent, normalized pupil PSF.
    ///
    /// A shared two-dimensional Cauchy proposal has the Airy tail's `r^-3`
    /// decay. It requires two independent unit variates, has no rejection loop
    /// or cutoff, and retains the rings at all wavelengths through separate
    /// importance weights. Bladed pupils use an analytic polygon Fourier
    /// transform and an edge-aligned product-Cauchy mixture, retaining their
    /// directional spikes without a pupil raster or an image convolution.
    /// The lens and diffraction draws must be independent
    /// of pixel, wavelength, BSDF and shutter draws. Geometric cameras ignore
    /// diffraction samples and wavelengths and return unit weights.
    ///
    /// Wavelengths are in nanometres. Scalar admission requires `lambda/D <=
    /// 0.01` for every wavelength, using the inscribed pupil diameter for a
    /// regular polygon. Call [`Self::admit_diffraction_depth`] on the
    /// primary visible surface (or infinite environment) before accepting the
    /// path. Primary surfaces must be opaque diffuse: a mirror/refractor's
    /// virtual-image depth is not its first-hit depth. Volumetric primary
    /// emission/scattering is not covered.
    #[allow(clippy::too_many_arguments)] // independent optical/sample coordinates
    pub fn generate_spectral_ray_from_tangent_offsets<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        diffraction_sample: LensSample,
        wavelengths_nm: &[f64; N],
    ) -> Result<SpectralCameraRay<N>, CameraError> {
        if self.diffraction.is_none() {
            return Ok(SpectralCameraRay {
                ray: self.generate_ray_from_tangent_offsets(cx, x_tan, y_tan, lens_sample)?,
                weights: [1.0; N],
            });
        }
        cx.checkpoint()?;
        if !x_tan.is_finite() || !y_tan.is_finite() {
            return Err(CameraError::NonFiniteInput);
        }
        let (minimum, maximum) = self.diffraction_wavelength_range(wavelengths_nm)?;
        let reference_nm = 0.5 * minimum + 0.5 * maximum;
        let scale = (reference_nm * 1.0e-9 / self.aperture.radius_m()) / (2.0 * PI);
        let (x, y, weights) = if let ApertureKind::RegularPolygon { vertices, .. } =
            &self.aperture.kind
        {
            let pupil = PolygonPupil::new(self.aperture.radius_m(), vertices);
            let (q, proposal) = pupil.sample([diffraction_sample.u(), diffraction_sample.v()]);
            let mut weights = [0.0; N];
            for (weight, &wavelength_nm) in weights.iter_mut().zip(wavelengths_nm) {
                cx.checkpoint()?;
                let ratio = reference_nm / wavelength_nm;
                *weight = ratio * ratio * pupil.density(q.map(|value| ratio * value)) / proposal;
                if !weight.is_finite() || *weight < 0.0 {
                    return Err(CameraError::InvalidDiffraction);
                }
            }
            (x_tan + scale * q[0], y_tan + scale * q[1], weights)
        } else {
            // Preserve the original circular sampling arithmetic and streams.
            let u = diffraction_sample.u();
            let remainder = 1.0 - u;
            // Inverse radial CDF: F(r) = 1 - 1/sqrt(1 + (r/s)^2).
            let radius = PROPOSAL_RADIUS * det::sqrt(u * (2.0 - u)) / remainder;
            let (sin_phi, cos_phi) = det::sin_cos(2.0 * PI * diffraction_sample.v());
            let x = x_tan + scale * radius * cos_phi;
            let y = y_tan + scale * radius * sin_phi;
            let mut weights = [0.0; N];
            for (weight, &wavelength_nm) in weights.iter_mut().zip(wavelengths_nm) {
                let ratio = reference_nm / wavelength_nm;
                let amplitude = jinc(ratio * radius);
                // Divide the normalized PSF by q(r) = (1-u)^3/(2*pi*s^2).
                *weight =
                    0.5 * PROPOSAL_RADIUS * PROPOSAL_RADIUS * ratio * ratio * amplitude * amplitude
                        / remainder
                        / remainder
                        / remainder;
                if !weight.is_finite() || *weight < 0.0 {
                    return Err(CameraError::InvalidDiffraction);
                }
            }
            (x, y, weights)
        };
        if !scale.is_finite() || scale <= 0.0 || !x.is_finite() || !y.is_finite() {
            return Err(CameraError::InvalidDiffraction);
        }
        let ray = Ray {
            origin: self.eye,
            dir: direct_unit(Vec3::new(
                self.forward.x + x * self.right.x + y * self.up.x,
                self.forward.y + x * self.right.y + y * self.up.y,
                self.forward.z + x * self.right.z + y * self.up.z,
            ))?,
        };
        cx.checkpoint()?;
        Ok(SpectralCameraRay { ray, weights })
    }

    /// Evaluate the selected focused PSF per unit angular-raster area.
    ///
    /// `tangent_offset` is the displacement in the camera's two tangent
    /// coordinates, and wavelength is in nanometres. The integral over those
    /// two coordinates is one. The on-axis density is physical pupil area
    /// divided by wavelength squared; it is not normalized to unit peak.
    /// This bounded analytic query allocates nothing. The same focus and
    /// primary-surface conditions as spectral ray generation still apply.
    pub fn in_focus_diffraction_density(
        &self,
        tangent_offset: [f64; 2],
        wavelength_nm: f64,
    ) -> Result<f64, CameraError> {
        if self.diffraction.is_none() || tangent_offset.iter().any(|v| !v.is_finite()) {
            return Err(CameraError::InvalidDiffraction);
        }
        self.diffraction_wavelength_range(&[wavelength_nm])?;
        let scale = (wavelength_nm * 1.0e-9 / self.aperture.radius_m()) / (2.0 * PI);
        let q = tangent_offset.map(|value| value / scale);
        let radius_squared = q[0] * q[0] + q[1] * q[1];
        if !scale.is_finite() || scale <= 0.0 || !radius_squared.is_finite() {
            return Err(CameraError::InvalidDiffraction);
        }
        let density = match &self.aperture.kind {
            ApertureKind::Circular { .. } => {
                let amplitude = jinc(det::sqrt(radius_squared));
                amplitude * amplitude / (4.0 * PI)
            }
            ApertureKind::RegularPolygon { vertices, .. } => {
                PolygonPupil::new(self.aperture.radius_m(), vertices).density(q)
            }
            ApertureKind::Pinhole => return Err(CameraError::InvalidDiffraction),
        } / scale
            / scale;
        if !density.is_finite() || density < 0.0 {
            return Err(CameraError::InvalidDiffraction);
        }
        Ok(density)
    }

    /// Admit the actual first visible depth under the declared phase budget.
    /// `None` denotes an infinitely distant environment. The smallest packet
    /// wavelength supplies the strongest test. The distance is axial, matching
    /// the camera's focus definition. Geometric cameras accept immediately.
    pub fn admit_diffraction_depth(
        &self,
        point: Option<Point3>,
        wavelengths_nm: &[f64],
    ) -> Result<(), CameraError> {
        let Some(policy) = self.diffraction else {
            return Ok(());
        };
        let (minimum, _) = self.diffraction_wavelength_range(wavelengths_nm)?;
        let inverse_depth = match point {
            Some(point) => {
                let depth = point.delta_from(self.eye).dot(self.forward);
                if !depth.is_finite() || depth <= 0.0 {
                    return Err(CameraError::DefocusOutsideDiffractionModel);
                }
                1.0 / depth
            }
            None => 0.0,
        };
        let phase = PI
            * (self.aperture.radius_m() / (minimum * 1.0e-9))
            * self.aperture.radius_m()
            * (inverse_depth - 1.0 / self.focus_distance_m).abs();
        if !phase.is_finite() || phase > policy.maximum_defocus_phase_rad {
            return Err(CameraError::DefocusOutsideDiffractionModel);
        }
        Ok(())
    }

    fn diffraction_wavelength_range(
        &self,
        wavelengths_nm: &[f64],
    ) -> Result<(f64, f64), CameraError> {
        self.admit_diffraction_projection()?;
        let radius_m = match &self.aperture.kind {
            ApertureKind::RegularPolygon {
                radius_m, blades, ..
            } => radius_m * det::cos(PI / f64::from(*blades)),
            _ => self.aperture.radius_m(),
        };
        let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
        for &wavelength in wavelengths_nm {
            let lambda_m = wavelength * 1.0e-9;
            let relative = 0.5 * (lambda_m / radius_m);
            if !wavelength.is_finite()
                || wavelength <= 0.0
                || lambda_m <= 0.0
                || !relative.is_finite()
                || relative <= 0.0
                || relative > 0.01
            {
                return Err(CameraError::InvalidDiffraction);
            }
            minimum = minimum.min(wavelength);
            maximum = maximum.max(wavelength);
        }
        if !minimum.is_finite() {
            return Err(CameraError::InvalidDiffraction);
        }
        Ok((minimum, maximum))
    }

    fn admit_diffraction_projection(&self) -> Result<(), CameraError> {
        if self
            .projection
            .focal_length_m()
            .is_some_and(|focal_length| self.focus_distance_m <= focal_length)
        {
            // A finite positive image plane exists only for an object beyond
            // the front focal plane. FOV-only cameras declare angular imaging.
            return Err(CameraError::InvalidDiffraction);
        }
        Ok(())
    }
}

fn jinc(x: f64) -> f64 {
    if x.abs() < 1.0e-4 {
        1.0 - x * x / 8.0
    } else {
        2.0 * bessel_j1(x) / x
    }
}
