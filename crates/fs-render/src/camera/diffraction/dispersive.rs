//! Longitudinal dispersion of coincident, paraxial thin elements.
//!
//! In index one, each element contributes `P(lambda) = (n(lambda)-1)*(Cf-Cb)` in inverse
//! metres and phase `-pi*P(lambda)*r^2/lambda`. Powers add because their
//! transmission functions multiply at the same pupil plane. The sensor is
//! at `1/t = P(reference)-1/focus`, so propagation from source depth z leaves
//! `pi*a^2/lambda * (1/z-1/focus+P(reference)-P(lambda))` at the pupil rim.
//! Compute the chromatic difference element by element, before summation.
//!
//! Reference: Kirchauer, TU Wien, section 4.1.2, equations (4.15)–(4.18):
//! <https://www.iue.tuwien.ac.at/phd/kirchauer/node51.html>.
//! These are ideal quadratic phase elements, not finite-thickness spherical
//! lenses. Surface aberrations, coatings, absorption and polarization are
//! deliberately not inferred from the supplied material phase-index law.
//! Geometric cameras may additionally retain an explicit surrounding material:
//! their inverse focal length is `(n_lens/n_exterior-1)*(Cf-Cb)`. The same
//! medium occupies both sides, so the usual object/sensor conjugacy survives.
//! See OpenStax, *University Physics III*, section 2.4, lens maker's equation:
//! <https://openstax.org/books/university-physics-volume-3/pages/2-4-thin-lenses>.

use std::sync::Arc;

use super::super::{CameraError, PhysicalCamera};
use super::{LensSample, PolygonPupil, SpectralCameraRay, annular};
use crate::dielectric::DielectricIor;
use fs_exec::Cx;
use fs_math::{dd::Dd, det};

/// One ideal thin optical element with a native material phase-index law.
///
/// Signed curvatures are reciprocals of surface radii, in inverse metres:
/// a converging biconvex element has positive front and negative back
/// curvature when light travels from the object toward the sensor. Zero is
/// a plane surface. The prescription separately selects the surrounding medium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinLensElement {
    ior: DielectricIor,
    curvatures_per_m: [f64; 2],
}

impl ThinLensElement {
    /// Admit finite signed surface curvatures and a validated native IOR law.
    /// Pupil-dependent paraxial admission happens when attaching the lens.
    pub fn try_new(
        ior: DielectricIor,
        front_curvature_per_m: f64,
        back_curvature_per_m: f64,
    ) -> Result<Self, CameraError> {
        if !front_curvature_per_m.is_finite()
            || !back_curvature_per_m.is_finite()
            || !(front_curvature_per_m - back_curvature_per_m).is_finite()
        {
            return Err(CameraError::InvalidDispersiveLens);
        }
        Ok(Self {
            ior,
            curvatures_per_m: [front_curvature_per_m, back_curvature_per_m],
        })
    }

    /// Native material phase index; group index never determines lens power.
    #[must_use]
    pub const fn ior(&self) -> &DielectricIor {
        &self.ior
    }

    /// Signed front/back surface curvatures, in inverse metres.
    #[must_use]
    pub const fn curvatures_per_m(self) -> [f64; 2] {
        self.curvatures_per_m
    }

    fn curvature_difference(self) -> f64 {
        self.curvatures_per_m[0] - self.curvatures_per_m[1]
    }
}

#[derive(Debug, PartialEq)]
struct Prescription {
    reference_wavelength_nm: f64,
    reference_power_per_m: f64,
    elements: Vec<ThinLensElement>,
    reference_indices: Vec<f64>,
    exterior_ior: Option<DielectricIor>,
    reference_exterior_index: f64,
}

/// Material-driven longitudinal chromatic aberration of a thin lens stack.
///
/// One to eight coincident ideal elements may combine positive and negative
/// powers, including achromatic doublets. Every wavelength uses the exact
/// native Cauchy, Sellmeier or tabulated phase index within their common source
/// support; neither a fitted Abbe number nor group index is used. Source-bound
/// infrared laws retain their own wavelength limits and knot-gap policies.
/// The fixed sensor is conjugate to the declared focus at the reference
/// wavelength. Its raster magnification is common to every spectral lane.
///
/// This unit-transmission phase model preserves the camera's normalized
/// radiance convention. It does not infer reflection, absorption, finite
/// thickness, element separation, lateral color or spherical aberration.
#[derive(Clone, Debug, PartialEq)]
pub struct DispersiveThinLens {
    prescription: Arc<Prescription>,
}

impl DispersiveThinLens {
    /// Hard work bound on coincident optical elements.
    pub const MAX_ELEMENTS: usize = 8;
    /// Maximum pupil-radius/surface-radius ratio in this paraxial model.
    pub const MAXIMUM_SURFACE_SLOPE: f64 = 0.05;
    /// Maximum nominal chief-ray or object/image pupil tangent slope.
    /// This does not truncate the mathematical scalar PSF's infinite tails.
    pub const MAXIMUM_PARAXIAL_SLOPE: f64 = 0.1;

    /// Construct a positive-power reference lens in index one from native materials.
    /// Only construction allocates; camera interpolation shares the immutable
    /// bounded prescription, and spectral evaluation allocates nothing.
    pub fn try_new(
        reference_wavelength_nm: f64,
        elements: &[ThinLensElement],
    ) -> Result<Self, CameraError> {
        Self::prepare(reference_wavelength_nm, elements, None)
    }

    /// Retain an explicit homogeneous surrounding phase-index law on both
    /// sides of a geometric thin-lens prescription. Native absolute indices
    /// determine the relative contrast at each original vacuum wavelength;
    /// neither an air-relative fit nor a group index is substituted.
    ///
    /// The camera's sampled exterior index must equal this law at every lane.
    /// Source bands and table gap limits remain intact. This is a transparent,
    /// unit-transmission ideal lens; exterior absorption/scattering along scene
    /// paths remains owned by the original transport. Wave diffraction refuses
    /// this explicit-medium prescription until its pupil propagator supports it.
    pub fn try_new_in_medium(
        reference_wavelength_nm: f64,
        elements: &[ThinLensElement],
        exterior_ior: DielectricIor,
    ) -> Result<Self, CameraError> {
        Self::prepare(reference_wavelength_nm, elements, Some(exterior_ior))
    }

    fn prepare(
        reference_wavelength_nm: f64,
        elements: &[ThinLensElement],
        exterior_ior: Option<DielectricIor>,
    ) -> Result<Self, CameraError> {
        if elements.is_empty() || elements.len() > Self::MAX_ELEMENTS {
            return Err(CameraError::InvalidDispersiveLens);
        }
        let reference_exterior_index = exterior_ior
            .map_or(Ok(1.0), |ior| ior.eval(reference_wavelength_nm))
            .map_err(|_| CameraError::InvalidDispersiveLens)?;
        let reference_indices = elements
            .iter()
            .map(|element| {
                element
                    .ior
                    .eval(reference_wavelength_nm)
                    .map_err(|_| CameraError::InvalidDispersiveLens)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let reference_power_per_m = compensated_power(
            elements
                .iter()
                .zip(&reference_indices)
                .map(|(element, &index)| {
                    Ok(
                        (index - reference_exterior_index) / reference_exterior_index
                            * element.curvature_difference(),
                    )
                }),
        )?;
        if reference_power_per_m <= 0.0 || !(1.0 / reference_power_per_m).is_finite() {
            return Err(CameraError::InvalidDispersiveLens);
        }
        Ok(Self {
            prescription: Arc::new(Prescription {
                reference_wavelength_nm,
                reference_power_per_m,
                elements: elements.to_vec(),
                reference_indices,
                exterior_ior,
                reference_exterior_index,
            }),
        })
    }

    /// Wavelength [nm] at which the camera's focus distance is sharp.
    #[must_use]
    pub fn reference_wavelength_nm(&self) -> f64 {
        self.prescription.reference_wavelength_nm
    }

    /// Coincident elements, retaining the exact material laws and curvatures.
    #[must_use]
    pub fn elements(&self) -> &[ThinLensElement] {
        &self.prescription.elements
    }

    /// Original surrounding material law; `None` is the index-one constructor.
    #[must_use]
    pub fn exterior_ior(&self) -> Option<&DielectricIor> {
        self.prescription.exterior_ior.as_ref()
    }

    /// Exterior phase index at an original vacuum wavelength [nm].
    pub fn exterior_phase_index(&self, wavelength_nm: f64) -> Result<f64, CameraError> {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        self.exterior_ior().map_or(Ok(1.0), |ior| {
            ior.eval(wavelength_nm)
                .map_err(|_| CameraError::DispersiveLensOutsideModel)
        })
    }

    /// Bind camera geometry to the actual local transport index. Equality is
    /// exact: the same native material evaluation can be shared by both owners.
    /// A different medium is not silently substituted into the lens maker law.
    pub fn admit_exterior_phase_index(
        &self,
        wavelength_nm: f64,
        actual_phase_index: f64,
    ) -> Result<(), CameraError> {
        if actual_phase_index != self.exterior_phase_index(wavelength_nm)? {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        Ok(())
    }

    /// Inclusive common support [vacuum nm] of the original material laws.
    /// Construction requires the reference wavelength to belong to every law.
    #[must_use]
    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        let exterior = self
            .exterior_ior()
            .map_or((0.0, f64::INFINITY), DielectricIor::wavelength_bounds_nm);
        self.elements()
            .iter()
            .fold(exterior, |(minimum, maximum), element| {
                let (low, high) = element.ior.wavelength_bounds_nm();
                (minimum.max(low), maximum.min(high))
            })
    }

    /// Positive reference optical power in inverse metres (dioptres).
    #[must_use]
    pub fn reference_power_per_m(&self) -> f64 {
        self.prescription.reference_power_per_m
    }

    /// Positive reference focal length in metres.
    #[must_use]
    pub fn reference_focal_length_m(&self) -> f64 {
        1.0 / self.reference_power_per_m()
    }

    /// Optical power [1/m] at a supported vacuum wavelength [nm].
    /// Signed powers at other wavelengths remain meaningful to the phase
    /// model; the camera separately admits its total pupil phase.
    pub fn power_per_m(&self, wavelength_nm: f64) -> Result<f64, CameraError> {
        let exterior = self.exterior_phase_index(wavelength_nm)?;
        compensated_power(self.elements().iter().map(|element| {
            let index = element
                .ior
                .eval(wavelength_nm)
                .map_err(|_| CameraError::DispersiveLensOutsideModel)?;
            Ok((index - exterior) / exterior * element.curvature_difference())
        }))
    }

    /// `P(reference)-P(wavelength)` [1/m], subtracting each material's
    /// phase indices before multiplying its curvature. This retains small
    /// chromatic differences in nearly cancelling doublet prescriptions.
    pub fn power_difference_per_m(&self, wavelength_nm: f64) -> Result<f64, CameraError> {
        if wavelength_nm == self.reference_wavelength_nm() {
            return Ok(0.0);
        }
        let exterior = self.exterior_phase_index(wavelength_nm)?;
        let reference_exterior = self.prescription.reference_exterior_index;
        compensated_power(
            self.elements()
                .iter()
                .zip(&self.prescription.reference_indices)
                .map(|(element, &reference_index)| {
                    let index = element
                        .ior
                        .eval(wavelength_nm)
                        .map_err(|_| CameraError::DispersiveLensOutsideModel)?;
                    let difference = if self.exterior_ior().is_none() {
                        // Preserve the original index-one arithmetic exactly.
                        reference_index - index
                    } else {
                        // Subtract cross products before division. Double-double
                        // products retain weak chromatic contrast even when
                        // both ordinary index ratios round to the same value.
                        let n0 = Dd::from_f64(reference_index);
                        let n = Dd::from_f64(index);
                        let e0 = Dd::from_f64(reference_exterior);
                        let e = Dd::from_f64(exterior);
                        ((n0 * e - n * e0) / (e0 * e)).to_f64()
                    };
                    Ok(difference * element.curvature_difference())
                }),
        )
    }

    /// Positive sensor conjugate [m] for a finite reference focus [m].
    /// The sensor plane remains common to all wavelengths.
    pub fn sensor_distance_m(&self, focus_distance_m: f64) -> Result<f64, CameraError> {
        let inverse_distance = self.reference_power_per_m() - 1.0 / focus_distance_m;
        let distance = 1.0 / inverse_distance;
        if !focus_distance_m.is_finite()
            || focus_distance_m <= 0.0
            || !inverse_distance.is_finite()
            || inverse_distance <= 0.0
            || !distance.is_finite()
            || distance <= 0.0
        {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        Ok(distance)
    }

    fn admit_camera(&self, camera: &PhysicalCamera) -> Result<(), CameraError> {
        let focal = camera
            .projection
            .focal_length_m()
            .ok_or(CameraError::DispersiveLensOutsideModel)?;
        let sensor = camera
            .projection
            .sensor_height_m()
            .ok_or(CameraError::DispersiveLensOutsideModel)?;
        let expected_focal = self.reference_focal_length_m();
        let tolerance = 64.0 * f64::EPSILON * focal.max(expected_focal);
        let radius = camera.aperture.radius_m();
        let distance = self.sensor_distance_m(camera.focus_distance_m)?;
        let half_tan = super::super::half_ratio(sensor, distance);
        if camera.aperture.is_pinhole()
            || camera.diffraction.is_some()
            || (self.exterior_ior().is_some() && camera.has_diffraction())
            || (focal - expected_focal).abs() > tolerance
            || radius / camera.focus_distance_m > Self::MAXIMUM_PARAXIAL_SLOPE
            || radius / distance > Self::MAXIMUM_PARAXIAL_SLOPE
            || !half_tan.is_finite()
            || half_tan <= 0.0
            || half_tan > Self::MAXIMUM_PARAXIAL_SLOPE
            || self.elements().iter().any(|element| {
                element
                    .curvatures_per_m
                    .iter()
                    .any(|curvature| radius * curvature.abs() > Self::MAXIMUM_SURFACE_SLOPE)
            })
        {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        Ok(())
    }
}

fn compensated_power(
    terms: impl Iterator<Item = Result<f64, CameraError>>,
) -> Result<f64, CameraError> {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for term in terms {
        let term = term?;
        let next = sum + term;
        if !term.is_finite() || !next.is_finite() {
            return Err(CameraError::InvalidDispersiveLens);
        }
        correction += if sum.abs() >= term.abs() {
            (sum - next) + term
        } else {
            (term - next) + sum
        };
        sum = next;
    }
    let result = sum + correction;
    if !result.is_finite() {
        return Err(CameraError::InvalidDispersiveLens);
    }
    Ok(result)
}

#[cfg(test)]
mod immersion_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dielectric::{CauchyIor, DielectricIorSample, SellmeierIor, TabulatedIor};

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn native_lensmaker_power_and_reference_sensor_conjugate() {
        let lens = DispersiveThinLens::try_new(
            550.0,
            &[
                ThinLensElement::try_new(CauchyIor::try_constant(1.5).unwrap().into(), 20.0, -20.0)
                    .unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(lens.reference_power_per_m(), 20.0);
        assert_eq!(lens.reference_focal_length_m(), 0.05);
        for wavelength in [380.0, 550.0, 780.0] {
            assert_eq!(lens.power_per_m(wavelength).unwrap(), 20.0);
            assert_eq!(lens.power_difference_per_m(wavelength).unwrap(), 0.0);
        }
        close(lens.sensor_distance_m(1.0).unwrap(), 1.0 / 19.0, 1.0e-17);
        assert!(lens.sensor_distance_m(0.05).is_err());
        assert!(lens.sensor_distance_m(0.04).is_err());
    }

    #[test]
    fn native_cauchy_sellmeier_and_table_keep_phase_index_and_units() {
        let cauchy = CauchyIor::try_new(1.5, 0.004, 0.0002).unwrap();
        let strengths = [1.03961212, 0.231792344, 1.01046945];
        let resonances = [0.00600069867, 0.0200179144, 103.560653];
        let sellmeier = SellmeierIor::try_new(strengths, resonances).unwrap();
        let indices = [1.53, 1.52, 1.51, 1.55, 1.50, 1.49, 1.48, 1.47, 1.46];
        let table = TabulatedIor::try_new(core::array::from_fn(|i| {
            DielectricIorSample::try_new(380.0 + 50.0 * i as f64, indices[i]).unwrap()
        }))
        .unwrap();
        for law in [cauchy.into(), sellmeier.into(), table.into()] {
            let lens = DispersiveThinLens::try_new(
                550.0,
                &[ThinLensElement::try_new(law, 30.0, -10.0).unwrap()],
            )
            .unwrap();
            for wavelength in [380.0, 405.0, 530.0, 550.0, 780.0] {
                let micrometres = wavelength * 0.001;
                let squared = micrometres * micrometres;
                let expected_index = match law {
                    DielectricIor::Cauchy(_) => 1.5 + 0.004 / squared + 0.0002 / squared / squared,
                    DielectricIor::Sellmeier(_) => {
                        let n2 = 1.0
                            + strengths
                                .into_iter()
                                .zip(resonances)
                                .map(|(b, c)| b * squared / (squared - c))
                                .sum::<f64>();
                        det::sqrt(n2)
                    }
                    DielectricIor::Tabulated(_) => {
                        let interval = (((wavelength - 380.0) / 50.0) as usize).min(7);
                        let t = (wavelength - (380.0 + 50.0 * interval as f64)) / 50.0;
                        (1.0 - t) * indices[interval] + t * indices[interval + 1]
                    }
                    #[cfg(feature = "material-state")]
                    DielectricIor::Water(_) => unreachable!("fixture contains no water law"),
                };
                close(
                    lens.power_per_m(wavelength).unwrap(),
                    40.0 * (expected_index - 1.0),
                    3.0e-13,
                );
                close(
                    lens.power_difference_per_m(wavelength).unwrap(),
                    40.0 * (law.eval(550.0).unwrap() - expected_index),
                    3.0e-13,
                );
            }
            assert!(lens.power_per_m(379.0).is_err());
            assert!(lens.power_difference_per_m(781.0).is_err());
        }
        assert_ne!(
            cauchy.eval(550.0).unwrap(),
            cauchy.group_index(550.0).unwrap()
        );
    }

    #[test]
    fn signed_native_doublet_cancels_two_colors_and_retains_secondary_spectrum() {
        let crown = CauchyIor::try_new(1.50, 0.004, 0.0001).unwrap();
        let flint = CauchyIor::try_new(1.62, 0.012, 0.0005).unwrap();
        let (blue, red, reference) = (486.1, 656.3, 550.0);
        let ratio = -(crown.eval(blue).unwrap() - crown.eval(red).unwrap())
            / (flint.eval(blue).unwrap() - flint.eval(red).unwrap());
        let crown_curvature = 20.0
            / ((crown.eval(reference).unwrap() - 1.0)
                + ratio * (flint.eval(reference).unwrap() - 1.0));
        let flint_curvature = ratio * crown_curvature;
        assert!(crown_curvature > 0.0 && flint_curvature < 0.0);
        let lens = DispersiveThinLens::try_new(
            reference,
            &[
                ThinLensElement::try_new(
                    crown.into(),
                    0.5 * crown_curvature,
                    -0.5 * crown_curvature,
                )
                .unwrap(),
                ThinLensElement::try_new(
                    flint.into(),
                    0.5 * flint_curvature,
                    -0.5 * flint_curvature,
                )
                .unwrap(),
            ],
        )
        .unwrap();
        close(lens.reference_power_per_m(), 20.0, 2.0e-14);
        close(
            lens.power_per_m(blue).unwrap(),
            lens.power_per_m(red).unwrap(),
            2.0e-13,
        );
        close(
            lens.power_difference_per_m(blue).unwrap(),
            lens.power_difference_per_m(red).unwrap(),
            2.0e-13,
        );
        assert!(lens.power_difference_per_m(blue).unwrap().abs() > 1.0e-5);
        assert_eq!(
            lens.power_difference_per_m(reference).unwrap().to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            compensated_power([Ok(1.0e16), Ok(1.0), Ok(-1.0e16)].into_iter()).unwrap(),
            1.0
        );
    }

    #[test]
    fn source_band_intersection_retains_infrared_phase_laws_and_knot_gaps() {
        let wide_samples = core::array::from_fn(|i| {
            let wavelength = 8000.0 + 750.0 * i as f64;
            DielectricIorSample::try_new_in_band(wavelength, 1.5 + (11000.0 - wavelength) * 2e-5)
                .unwrap()
        });
        assert!(TabulatedIor::try_new_in_band(wide_samples, 749.0).is_err());
        let wide = TabulatedIor::try_new_in_band(wide_samples, 750.0).unwrap();
        let narrow = TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(9000.0 + 500.0 * i as f64, 1.6).unwrap()
            }),
            500.0,
        )
        .unwrap();
        let lens = DispersiveThinLens::try_new(
            11000.0,
            &[
                ThinLensElement::try_new(wide.into(), 20.0, -20.0).unwrap(),
                ThinLensElement::try_new(narrow.into(), 10.0, -10.0).unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(lens.wavelength_bounds_nm(), (9000.0, 13000.0));
        assert_eq!(lens.elements()[0].ior(), &DielectricIor::Tabulated(wide));
        for wavelength in [9000.0, 9750.0, 11000.0, 12250.0, 13000.0] {
            close(
                lens.power_difference_per_m(wavelength).unwrap(),
                (wavelength - 11000.0) * 8e-4,
                2e-14,
            );
        }
        for wavelength in [8999.0, 13001.0] {
            assert!(lens.power_per_m(wavelength).is_err());
            assert!(lens.power_difference_per_m(wavelength).is_err());
        }
    }

    #[test]
    fn thin_element_prescription_refuses_unsupported_inputs() {
        let law = CauchyIor::try_constant(1.5).unwrap().into();
        let converging = ThinLensElement::try_new(law, 20.0, -20.0).unwrap();
        assert!(ThinLensElement::try_new(law, f64::NAN, 0.0).is_err());
        assert!(ThinLensElement::try_new(law, f64::MAX, -f64::MAX).is_err());
        assert!(DispersiveThinLens::try_new(550.0, &[]).is_err());
        assert!(DispersiveThinLens::try_new(550.0, &[converging; 9]).is_err());
        assert!(DispersiveThinLens::try_new(550.0, &[converging; 8]).is_ok());
        assert!(DispersiveThinLens::try_new(379.0, &[converging]).is_err());
        assert!(
            DispersiveThinLens::try_new(
                550.0,
                &[ThinLensElement::try_new(law, 0.0, 0.0).unwrap(),]
            )
            .is_err()
        );
        assert!(
            DispersiveThinLens::try_new(
                550.0,
                &[ThinLensElement::try_new(law, -20.0, 20.0).unwrap(),]
            )
            .is_err()
        );
    }
}

impl PhysicalCamera {
    /// Enable material-driven longitudinal dispersion for an explicit
    /// finite-defocus or layered diffraction camera. The focal/sensor
    /// projection must name this lens's reference focal length (within 64
    /// relative machine epsilons). Its sensor height is fixed; focus changes
    /// move the sensor conjugate and therefore change the active raster scale.
    ///
    /// Geometric and focused-Airy policies refuse this attachment because they
    /// cannot propagate the wavelength-dependent pupil phase. For a source
    /// sharp at the reference wavelength, declare its plane with
    /// [`super::DefocusedDiffraction`] at the camera's focus distance.
    pub fn with_dispersive_lens(mut self, lens: DispersiveThinLens) -> Result<Self, CameraError> {
        if self.defocused_diffraction.is_none() && self.layered_diffraction.is_none() {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        lens.admit_camera(&self)?;
        self.dispersive_lens = Some(lens);
        Ok(self)
    }

    /// Explicitly select geometric longitudinal chromatic aberration using
    /// the original thin-element prescription and aperture sampler.
    ///
    /// The fixed sensor is conjugate to the reference focus. Each wavelength
    /// therefore needs its own outgoing ray through the sampled pupil. This
    /// geometric-optics approximation omits diffraction; it never superposes
    /// geometric blur on the wave PSF. The original diffraction attachment
    /// remains separate. The scene's phase index at each sampled pupil point
    /// must match the prescription's original exterior law at that wavelength.
    pub fn with_geometric_dispersive_lens(
        mut self,
        lens: DispersiveThinLens,
    ) -> Result<Self, CameraError> {
        if self.has_diffraction() {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        lens.admit_camera(&self)?;
        self.dispersive_lens = Some(lens);
        Ok(self)
    }

    /// Generate separate geometric rays at supplied vacuum wavelengths [nm].
    /// A shared original pupil draw couples the lanes without changing their
    /// aperture-average radiance measure; every lane has unit camera weight.
    /// Achromatic cameras reuse their original ray bit for bit.
    ///
    /// For pupil coordinate `p`, raster tangent `u`, sensor distance `t` and
    /// native phase power `P(lambda)`, the backward ray has tangent
    /// `u + p*(1/t-P(lambda))`. Since `1/t=P(reference)-1/focus`, the stable
    /// material index difference supplies the chromatic change. Group index
    /// does not determine this direction. All wavelengths retain the original
    /// source bands/gaps and optical-opening admission
    /// `lambda_vac/(n_exterior*width) <= 0.01`.
    /// The entire pupil's outgoing paraxial slope must remain at most 0.1.
    /// No wavelength extrapolation, diffraction, or lateral color is inferred.
    pub fn generate_geometric_spectral_rays_from_tangent_offsets<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
    ) -> Result<[crate::charts::Ray; N], CameraError> {
        cx.checkpoint()?;
        if self.has_diffraction() {
            return Err(CameraError::DiffractionRequiresSpectralPath);
        }
        let Some(lens) = &self.dispersive_lens else {
            return Ok([self.generate_ray_from_tangent_offsets(cx, x_tan, y_tan, lens_sample)?; N]);
        };
        self.admit_dispersive_raster(x_tan, y_tan)?;
        self.diffraction_wavelength_range(wavelengths_nm)?;
        let reference_ray =
            self.generate_geometric_ray_from_tangent_offsets(cx, x_tan, y_tan, lens_sample)?;
        let [pupil_x, pupil_y] = self.aperture.sample(lens_sample);
        let raster_slope = det::hypot(x_tan, y_tan);
        let mut rays = [reference_ray; N];
        for (ray, &wavelength_nm) in rays.iter_mut().zip(wavelengths_nm) {
            cx.checkpoint()?;
            let difference = lens.power_difference_per_m(wavelength_nm)?;
            let inverse_conjugate = difference - 1.0 / self.focus_distance_m;
            let maximum_slope = raster_slope + self.aperture.radius_m() * inverse_conjugate.abs();
            if !maximum_slope.is_finite()
                || maximum_slope > DispersiveThinLens::MAXIMUM_PARAXIAL_SLOPE
            {
                return Err(CameraError::DispersiveLensOutsideModel);
            }
            if difference != 0.0 {
                let x = x_tan + pupil_x * inverse_conjugate;
                let y = y_tan + pupil_y * inverse_conjugate;
                ray.dir = super::super::direct_unit(fs_geom::Vec3::new(
                    self.forward.x + x * self.right.x + y * self.up.x,
                    self.forward.y + x * self.right.y + y * self.up.y,
                    self.forward.z + x * self.right.z + y * self.up.z,
                ))?;
            }
        }
        cx.checkpoint()?;
        Ok(rays)
    }

    /// Selected material-based optical prescription, if explicitly enabled.
    #[must_use]
    pub const fn dispersive_lens(&self) -> Option<&DispersiveThinLens> {
        self.dispersive_lens.as_ref()
    }

    /// Actual sensor distance [m] for an enabled material-based thin lens.
    /// Other cameras retain their declared angular projection and return None.
    #[must_use]
    pub fn dispersive_sensor_distance_m(&self) -> Option<f64> {
        self.dispersive_lens
            .as_ref()
            .map(|lens| 1.0 / (lens.reference_power_per_m() - 1.0 / self.focus_distance_m))
    }

    pub(in crate::camera) fn admit_dispersive_lens(&self) -> Result<(), CameraError> {
        if let Some(lens) = &self.dispersive_lens {
            lens.admit_camera(self)?;
        }
        Ok(())
    }

    pub(crate) fn admit_dispersive_raster(&self, x: f64, y: f64) -> Result<(), CameraError> {
        if self.dispersive_lens.is_some()
            && fs_math::det::hypot(x, y) > DispersiveThinLens::MAXIMUM_PARAXIAL_SLOPE
        {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        Ok(())
    }

    pub(super) fn dispersive_proposal_phase(
        &self,
        wavelengths_nm: &[f64],
        reference_nm: f64,
    ) -> Result<f64, CameraError> {
        let mut maximum = 0.0_f64;
        for &wavelength in wavelengths_nm {
            let phase = if self.layered_diffraction.is_some() {
                self.layered_maximum_phase(wavelength)?
            } else {
                self.defocus_phase(
                    self.defocused_diffraction
                        .ok_or(CameraError::InvalidDiffraction)?,
                    wavelength,
                )?
                .abs()
            };
            // Geometric blur radius is 2|phase| in this wavelength's q.
            // Convert to the proposal's reference-wavelength coordinates.
            maximum = maximum.max(phase * (wavelength / reference_nm));
        }
        Ok(maximum)
    }

    pub(super) fn generate_dispersive_spectral_ray<const N: usize>(
        &self,
        cx: &Cx<'_>,
        raster_tangents: [f64; 2],
        sample: LensSample,
        wavelengths_nm: &[f64; N],
        reference_nm: f64,
        scale: f64,
    ) -> Result<SpectralCameraRay<N>, CameraError> {
        self.admit_dispersive_raster(raster_tangents[0], raster_tangents[1])?;
        let maximum_phase = self.dispersive_proposal_phase(wavelengths_nm, reference_nm)?;
        let (q, proposal) = match &self.aperture.kind {
            super::super::ApertureKind::RegularPolygon { vertices, .. } => {
                PolygonPupil::new(self.aperture.radius_m(), vertices)
                    .sample_phase_range([sample.u(), sample.v()], maximum_phase)
            }
            _ => {
                let (radius, proposal) = annular::sample_radius(
                    self.aperture.obscuration_ratio(),
                    maximum_phase,
                    sample.u(),
                );
                let (sin, cos) = det::sin_cos(2.0 * super::PI * sample.v());
                ([radius * cos, radius * sin], proposal)
            }
        };
        let policy = self
            .defocused_diffraction
            .ok_or(CameraError::InvalidDiffraction)?;
        let mut weights = [0.0; N];
        for (weight, &wavelength) in weights.iter_mut().zip(wavelengths_nm) {
            cx.checkpoint()?;
            let ratio = reference_nm / wavelength;
            let density = self.layered_pupil_density(
                q.map(|value| value * ratio),
                self.defocus_phase(policy, wavelength)?,
            )?;
            *weight = ratio * ratio * density / proposal;
            if !weight.is_finite() || *weight < 0.0 {
                return Err(CameraError::InvalidDiffraction);
            }
        }
        let x = raster_tangents[0] + scale * q[0];
        let y = raster_tangents[1] + scale * q[1];
        if !scale.is_finite() || scale <= 0.0 || !x.is_finite() || !y.is_finite() {
            return Err(CameraError::InvalidDiffraction);
        }
        let ray = crate::charts::Ray {
            origin: self.eye,
            dir: super::direct_unit(fs_geom::Vec3::new(
                self.forward.x + x * self.right.x + y * self.up.x,
                self.forward.y + x * self.right.y + y * self.up.y,
                self.forward.z + x * self.right.z + y * self.up.z,
            ))?,
        };
        cx.checkpoint()?;
        Ok(SpectralCameraRay { ray, weights })
    }
}
