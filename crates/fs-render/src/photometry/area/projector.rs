//! An ideal paraxial projector expressed as an exit-pupil light field.
//!
//! The supplied clear flux is uniform in pupil area and in coordinates of the
//! rectangular gobo plane. This is an explicit source model, not a recovered
//! lamp/condenser distribution or a solved lens train. For gobo distance `s`,
//! focal length `f`, pupil coordinate `a` and outgoing slope `u = w.xy/w.z`,
//! the thin-lens map is `g = a*(1-s/f) - s*u`.
//!
//! With gobo area `G`, pupil area `A` and passive transmission `M`, radiance is
//! `L = clear_flux/A * s²/G * M(g)/cos(theta)^4`. The factor follows from
//! `dOmega = cos(theta)^3 du.x du.y`, `dg.x dg.y = s² du.x du.y` and the extra
//! projected-area cosine in radiometric flux. Consequently
//! `L*cos(theta)*dA*dOmega = clear_flux/(A*G)*M(g)*dA*dg`, and total transmitted
//! flux is exactly clear flux times mean gobo transmission in the represented
//! piecewise-constant model, apart from floating-point arithmetic.
//!
//! An optional ideal beam splitter adds outgoing slope offsets `d_k` with
//! absolute power fractions `t_k`. The mask becomes
//! `sum_k t_k * M(a*(1-s/f) - s*(u-d_k))`; the outgoing `cos(theta)^4`
//! denominator remains shared. Slope translation has unit Jacobian, so total
//! flux is clear flux times mean gobo transmission times `sum_k t_k`.
//! Overlapping images add; at focus each copy shifts by `q*d_k` metres.
//!
//! At `1/q = 1/f - 1/s`, screen coordinates obey `g = -s*p/q` independently
//! of pupil position. At another screen distance `z`, the pupil contributes
//! the true geometric defocus displacement `(1-z/q)*a`. Diffraction,
//! aberrations and nonuniform condenser illumination are not modeled. Lens
//! transmission losses must already be included in the supplied clear-output
//! lumens; that calibration does not recover omitted spatial/angular physics.
//!
//! References: <https://openstax.org/books/university-physics-volume-3/pages/2-4-thin-lenses>
//! and <https://www.pbr-book.org/4ed/Light_Sources/Point_Lights>.

use std::sync::{Arc, OnceLock};

use fs_geom::{Point3, Vec3};

use crate::instances::RigidTransform;
use crate::photometry::source::SourceError;

mod beam_splitter;
pub use beam_splitter::{MAX_PROJECTOR_BEAMS, ParaxialBeamSplitter, ProjectorBeam};
mod diffuser;
pub use diffuser::ParaxialBoxDiffuser;
use diffuser::PreparedDiffuser;

/// Maximum admitted samples in one passive gobo transmission map.
pub const MAX_GOBO_TEXELS: usize = 1_048_576;

/// A passive, achromatic transmission map on a centered rectangular gobo.
///
/// Values are piecewise constant, row-major with X changing fastest. The
/// support is `[-width/2,width/2) × [-height/2,height/2)` in metres; everything
/// outside is opaque. X points right and Y points up in the optical frame.
/// Clones share the immutable sample storage and any prepared frost integral,
/// retained until the final clone of that data is dropped.
/// Equality compares the represented map, independently of preparation state.
#[derive(Debug, Clone)]
pub struct GoboTransmission {
    data: Arc<GoboData>,
}

#[derive(Debug)]
struct GoboData {
    width_m: f64,
    height_m: f64,
    area_m2: f64,
    columns: usize,
    rows: usize,
    values: Vec<f64>,
    mean: f64,
    maximum: f64,
    diffusion_integral: OnceLock<Arc<diffuser::PositiveIntegral>>,
}

impl PartialEq for GoboTransmission {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
            || (self.data.width_m == other.data.width_m
                && self.data.height_m == other.data.height_m
                && self.data.columns == other.data.columns
                && self.data.rows == other.data.rows
                && self.data.values == other.data.values)
    }
}

impl GoboTransmission {
    /// Admit positive finite dimensions and one finite `[0,1]` transmission
    /// per texel. A wholly opaque map is valid and emits zero flux.
    pub fn try_new(
        width_m: f64,
        height_m: f64,
        columns: usize,
        rows: usize,
        transmission: Vec<f64>,
    ) -> Result<Self, SourceError> {
        let count = columns
            .checked_mul(rows)
            .ok_or(SourceError::DistributionTooLarge)?;
        if count > MAX_GOBO_TEXELS {
            return Err(SourceError::DistributionTooLarge);
        }
        if columns == 0
            || rows == 0
            || transmission.len() != count
            || !width_m.is_finite()
            || !height_m.is_finite()
            || width_m <= 0.0
            || height_m <= 0.0
        {
            return Err(SourceError::InvalidDistribution);
        }
        let area_m2 = width_m * height_m;
        if !area_m2.is_finite() || area_m2 <= 0.0 || width_m * 0.5 == 0.0 || height_m * 0.5 == 0.0 {
            return Err(SourceError::NumericalRange);
        }
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        let mut maximum = 0.0_f64;
        for &value in &transmission {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(SourceError::InvalidFilter);
            }
            maximum = maximum.max(value);
            let next = sum + value;
            correction += if sum >= value {
                (sum - next) + value
            } else {
                (value - next) + sum
            };
            sum = next;
        }
        let mean = (sum + correction) / count as f64;
        if !mean.is_finite() || !(0.0..=1.0).contains(&mean) || (maximum > 0.0 && mean == 0.0) {
            return Err(SourceError::NumericalRange);
        }
        Ok(Self {
            data: Arc::new(GoboData {
                width_m,
                height_m,
                area_m2,
                columns,
                rows,
                values: transmission,
                mean,
                maximum,
                diffusion_integral: OnceLock::new(),
            }),
        })
    }

    /// Full width of the represented gobo plane in metres.
    #[must_use]
    pub fn width_m(&self) -> f64 {
        self.data.width_m
    }

    /// Full height of the represented gobo plane in metres.
    #[must_use]
    pub fn height_m(&self) -> f64 {
        self.data.height_m
    }

    /// Physical area of the rectangular gobo support in square metres.
    #[must_use]
    pub fn area_m2(&self) -> f64 {
        self.data.area_m2
    }

    /// Samples along the optical X axis.
    #[must_use]
    pub fn columns(&self) -> usize {
        self.data.columns
    }

    /// Samples along the optical Y axis.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.data.rows
    }

    /// Area-mean passive transmission, without replacement of blocked flux.
    #[must_use]
    pub fn mean_transmission(&self) -> f64 {
        self.data.mean
    }

    /// Evaluate the represented piecewise-constant transmission in metres.
    /// Finite positions outside the rectangle return zero.
    pub fn transmission_at(&self, x_m: f64, y_m: f64) -> Result<f64, SourceError> {
        if !x_m.is_finite() || !y_m.is_finite() {
            return Err(SourceError::NumericalRange);
        }
        let half_width = 0.5 * self.data.width_m;
        let half_height = 0.5 * self.data.height_m;
        if x_m < -half_width || x_m >= half_width || y_m < -half_height || y_m >= half_height {
            return Ok(0.0);
        }
        let x = ((x_m / self.data.width_m + 0.5) * self.data.columns as f64) as usize;
        let y = ((y_m / self.data.height_m + 0.5) * self.data.rows as f64) as usize;
        // A point immediately below the excluded upper edge can round onto
        // the final integer boundary. It remains in the final admitted texel.
        Ok(self.data.values
            [y.min(self.data.rows - 1) * self.data.columns + x.min(self.data.columns - 1)])
    }
}

/// A calibrated ideal thin-lens light field on a native planar mesh aperture.
///
/// `frame` maps optical coordinates into the emitter primitive's local frame;
/// native rigid-instance placement is composed during scene admission. The
/// aperture lies in optical Z=0 and emits along +Z. Lens and gobo distances are
/// positive metres. A gobo beyond the focal plane has a real forward image;
/// a gobo at the focal plane is collimated; a nearer gobo has a virtual image.
#[derive(Debug, Clone)]
pub struct ThinLensProjector {
    frame: RigidTransform,
    gobo_distance_m: f64,
    focal_length_m: f64,
    pupil_factor: f64,
    focus_distance_m: Option<f64>,
    radiance_factor: f64,
    gobo: GoboTransmission,
    beam_splitter: Option<ParaxialBeamSplitter>,
    diffuser: Option<PreparedDiffuser>,
    transmitted_fraction: f64,
}

impl ThinLensProjector {
    /// Admit an ideal geometric projector with the supplied passive gobo.
    pub fn try_new(
        frame: RigidTransform,
        gobo_distance_m: f64,
        focal_length_m: f64,
        gobo: GoboTransmission,
    ) -> Result<Self, SourceError> {
        if !gobo_distance_m.is_finite()
            || !focal_length_m.is_finite()
            || gobo_distance_m <= 0.0
            || focal_length_m <= 0.0
        {
            return Err(SourceError::InvalidDistribution);
        }
        let ratio = gobo_distance_m / focal_length_m;
        let pupil_factor = 1.0 - ratio;
        let radiance_factor =
            (gobo_distance_m / gobo.width_m()) * (gobo_distance_m / gobo.height_m());
        if !ratio.is_finite()
            || ratio <= 0.0
            || !pupil_factor.is_finite()
            || !radiance_factor.is_finite()
            || radiance_factor <= 0.0
            || (pupil_factor == 0.0 && gobo_distance_m != focal_length_m)
        {
            return Err(SourceError::NumericalRange);
        }
        let focus_distance_m = if pupil_factor == 0.0 {
            None
        } else {
            let distance = -gobo_distance_m / pupil_factor;
            if !distance.is_finite() || distance == 0.0 {
                return Err(SourceError::NumericalRange);
            }
            Some(distance)
        };
        Ok(Self {
            frame,
            gobo_distance_m,
            focal_length_m,
            pupil_factor,
            focus_distance_m,
            radiance_factor,
            transmitted_fraction: gobo.mean_transmission(),
            gobo,
            beam_splitter: None,
            diffuser: None,
        })
    }

    /// Change the positive gobo separation and focal length in metres.
    ///
    /// The original optical frame, gobo, splitter, diffuser law and passive
    /// transmission are retained. Native lens admission recomputes focus and
    /// radiance; changing the gobo distance refreshes only the diffuser's
    /// metric widths and shares its prepared integral. Focal-length-only
    /// changes retain the complete prepared diffuser. Equal distances return
    /// an unchanged clone, and refusal leaves this snapshot intact.
    pub fn with_lens_distances(
        &self,
        gobo_distance_m: f64,
        focal_length_m: f64,
    ) -> Result<Self, SourceError> {
        self.with_optical_state(
            gobo_distance_m,
            focal_length_m,
            self.gobo.clone(),
            self.beam_splitter.clone(),
            self.diffuser().copied(),
        )
    }

    /// Admit a complete optical state while retaining this projector's frame.
    ///
    /// `None` removes the corresponding splitter or diffuser. The final gobo,
    /// distances and passive components are admitted together: a removed frost
    /// filter is never prepared at an obsolete gobo size or lens distance.
    /// Unchanged components retain their transmitted-fraction bits, including
    /// a lens-only update. A changed gobo recomputes the native radiance factor;
    /// its opaque exterior and physical dimensions remain part of the source.
    ///
    /// Frost integrals belong to shared immutable gobo data, independently of
    /// diffuser angles, fractions and lens distances. Preparation is fallible
    /// and introduces no per-ray cache access. This method changes the optical
    /// transmission, not a source's clear-output photometric calibration.
    pub fn with_optical_state(
        &self,
        gobo_distance_m: f64,
        focal_length_m: f64,
        gobo: GoboTransmission,
        beam_splitter: Option<ParaxialBeamSplitter>,
        diffuser: Option<ParaxialBoxDiffuser>,
    ) -> Result<Self, SourceError> {
        let same_gobo = Arc::ptr_eq(&gobo.data, &self.gobo.data);
        let same_splitter = beam_splitter == self.beam_splitter;
        let same_diffuser = diffuser.as_ref() == self.diffuser();
        let same_components = same_gobo && same_splitter && same_diffuser;
        let same_distances =
            gobo_distance_m == self.gobo_distance_m && focal_length_m == self.focal_length_m;
        if same_components && same_distances {
            return Ok(self.clone());
        }
        // A distinct supplied map is adopted even if its values equal the
        // original. Shared fixture controls can thereby consolidate separately
        // allocated originals without scanning their texels per emitter.
        let mut result = if same_gobo && same_distances {
            self.clone()
        } else {
            Self::try_new(self.frame, gobo_distance_m, focal_length_m, gobo.clone())?
        };
        let transmitted_fraction = if same_components {
            self.transmitted_fraction
        } else {
            let splitter = beam_splitter
                .as_ref()
                .map_or(1.0, |value| value.transmitted_fraction());
            let frost = diffuser.map_or(1.0, |law| law.transmission());
            if gobo.mean_transmission() == 0.0 || splitter == 0.0 || frost == 0.0 {
                0.0
            } else {
                diffuser::product(
                    diffuser::product(gobo.mean_transmission(), splitter)?,
                    frost,
                )?
            }
        };
        result.diffuser = if same_gobo && same_diffuser {
            if gobo_distance_m == self.gobo_distance_m {
                self.diffuser.clone()
            } else {
                self.diffuser
                    .as_ref()
                    .map(|prepared| prepared.with_distance(&gobo, gobo_distance_m))
                    .transpose()?
            }
        } else {
            diffuser
                .map(|law| PreparedDiffuser::try_new(law, &gobo, gobo_distance_m))
                .transpose()?
        };
        result.beam_splitter = if same_splitter {
            self.beam_splitter.clone()
        } else {
            beam_splitter
        };
        result.transmitted_fraction = transmitted_fraction;
        Ok(result)
    }

    /// Replace outgoing beam replication with this ideal passive splitter.
    ///
    /// Gobo geometry, focus and the clear-output calibration are retained.
    /// Replacing a splitter does not compound the previous splitter's losses.
    /// The full aperture's radiance bound is checked during scene admission.
    pub fn with_beam_splitter(self, splitter: ParaxialBeamSplitter) -> Result<Self, SourceError> {
        self.with_optical_state(
            self.gobo_distance_m,
            self.focal_length_m,
            self.gobo.clone(),
            Some(splitter),
            self.diffuser().copied(),
        )
    }

    /// Optional outgoing beam replication, or the original unsplit light field.
    #[must_use]
    pub const fn beam_splitter(&self) -> Option<&ParaxialBeamSplitter> {
        self.beam_splitter.as_ref()
    }

    /// Replace outgoing slope redistribution with an ideal passive diffuser.
    ///
    /// Gobo, focus, splitter and clear-output calibration are retained. The
    /// scattered fraction changes the angular distribution; only the absolute
    /// transmission changes total flux. A shared positive-sum range tree is
    /// prepared once per shared gobo for nonzero scattering, with at most
    /// `4*MAX_GOBO_TEXELS` f64 storage slots including power-of-two padding
    /// (32 MiB at the current cap). Concurrent first preparations may build
    /// redundant temporary trees; completed projectors share the retained tree.
    /// Too much padding or unavailable storage refuses. Source evaluation has
    /// no allocation or additional random draws. Replacing the diffuser does
    /// not compound its previous loss.
    pub fn with_diffuser(self, diffuser: ParaxialBoxDiffuser) -> Result<Self, SourceError> {
        self.with_optical_state(
            self.gobo_distance_m,
            self.focal_length_m,
            self.gobo.clone(),
            self.beam_splitter.clone(),
            Some(diffuser),
        )
    }

    /// Declared diffuser state, or the original exit-pupil light field.
    #[must_use]
    pub fn diffuser(&self) -> Option<&ParaxialBoxDiffuser> {
        self.diffuser.as_ref().map(|prepared| &prepared.law)
    }

    /// Surviving fraction of clear-output flux after gobo, splitter and diffuser.
    #[must_use]
    pub const fn transmitted_fraction(&self) -> f64 {
        self.transmitted_fraction
    }

    /// Optical-to-primitive rigid frame; optical +Z is the emitting direction.
    #[must_use]
    pub const fn frame(&self) -> RigidTransform {
        self.frame
    }

    /// Positive axial separation between the gobo and lens planes in metres.
    #[must_use]
    pub const fn gobo_distance_m(&self) -> f64 {
        self.gobo_distance_m
    }

    /// Positive thin-lens focal length in metres.
    #[must_use]
    pub const fn focal_length_m(&self) -> f64 {
        self.focal_length_m
    }

    /// Signed image distance in metres, or `None` for exact collimation.
    /// Positive values are real forward images; negative values are virtual.
    #[must_use]
    pub const fn focus_distance_m(&self) -> Option<f64> {
        self.focus_distance_m
    }

    /// Original immutable gobo transmission data and physical dimensions.
    #[must_use]
    pub const fn gobo(&self) -> &GoboTransmission {
        &self.gobo
    }

    /// Luminance divided by clear flux per pupil area, in inverse steradians.
    ///
    /// Inputs are in the optical frame. The point supplies aperture X/Y;
    /// scene admission establishes its plane. A finite nonzero outgoing
    /// direction is normalized. Backward rays and rays outside the bounded
    /// gobo return zero before the grazing-angle Jacobian is evaluated.
    pub fn luminance_per_clear_flux_per_area(
        &self,
        optical_point: Point3,
        optical_dir: Vec3,
    ) -> Result<f64, SourceError> {
        if [optical_point.x, optical_point.y, optical_point.z]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(SourceError::NumericalRange);
        }
        let direction = normalized_direction(optical_dir)?;
        if direction.z <= 0.0 || self.transmitted_fraction == 0.0 {
            return Ok(0.0);
        }
        let pupil_x = optical_point.x * self.pupil_factor;
        let pupil_y = optical_point.y * self.pupil_factor;
        if !pupil_x.is_finite() || !pupil_y.is_finite() {
            return Err(SourceError::NumericalRange);
        }
        if let Some(diffuser) = &self.diffuser
            && !diffuser.is_identity()
        {
            return self.diffused_luminance(diffuser, [pupil_x, pupil_y], direction);
        }
        if let Some(splitter) = &self.beam_splitter {
            return self.split_luminance(splitter, [pupil_x, pupil_y], direction);
        }
        // Compare homogeneous coordinates first. Forming w.xy/w.z before
        // rejecting an unilluminated grazing ray could overflow needlessly.
        let x = pupil_x * direction.z - self.gobo_distance_m * direction.x;
        let y = pupil_y * direction.z - self.gobo_distance_m * direction.y;
        let half_width = 0.5 * self.gobo.width_m() * direction.z;
        let half_height = 0.5 * self.gobo.height_m() * direction.z;
        if x < -half_width || x >= half_width || y < -half_height || y >= half_height {
            return Ok(0.0);
        }
        let transmission = self
            .gobo
            .transmission_at(x / direction.z, y / direction.z)?;
        if transmission == 0.0 {
            return Ok(0.0);
        }
        let cosine_squared = direction.z * direction.z;
        let value = self.radiance_factor * transmission / cosine_squared / cosine_squared;
        if !value.is_finite() || value <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        Ok(value)
    }

    fn diffused_luminance(
        &self,
        diffuser: &PreparedDiffuser,
        pupil: [f64; 2],
        direction: Vec3,
    ) -> Result<f64, SourceError> {
        let mask_at = |offset: [f64; 2]| -> Result<f64, SourceError> {
            let x = pupil[0] * direction.z
                - self.gobo_distance_m * (direction.x - offset[0] * direction.z);
            let y = pupil[1] * direction.z
                - self.gobo_distance_m * (direction.y - offset[1] * direction.z);
            if !x.is_finite() || !y.is_finite() {
                return Err(SourceError::NumericalRange);
            }
            // Retain possible boundary mass even when a very narrow kernel's
            // expanded support rounds back onto the original gobo boundary.
            let half = [self.gobo.width_m(), self.gobo.height_m()].map(|size| 0.5 * size);
            let bound_x = ((half[0] + diffuser.half_width_m[0]).next_up() * direction.z).next_up();
            let bound_y = ((half[1] + diffuser.half_width_m[1]).next_up() * direction.z).next_up();
            if x < -bound_x || x > bound_x || y < -bound_y || y > bound_y {
                return Ok(0.0);
            }
            diffuser.mask(&self.gobo, x / direction.z, y / direction.z)
        };
        let transmission = if let Some(splitter) = &self.beam_splitter {
            let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
            for beam in splitter.beams() {
                if beam.throughput() == 0.0 {
                    continue;
                }
                let value = diffuser::product(beam.throughput(), mask_at(beam.slope_offset())?)?;
                let next = sum + value;
                correction += if sum >= value {
                    (sum - next) + value
                } else {
                    (value - next) + sum
                };
                sum = next;
            }
            sum + correction
        } else {
            mask_at([0.0; 2])?
        };
        if transmission == 0.0 {
            return Ok(0.0);
        }
        let cosine_squared = direction.z * direction.z;
        let value = self.radiance_factor * transmission / cosine_squared / cosine_squared;
        if !value.is_finite() || value <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        Ok(value)
    }

    fn split_luminance(
        &self,
        splitter: &ParaxialBeamSplitter,
        pupil: [f64; 2],
        direction: Vec3,
    ) -> Result<f64, SourceError> {
        let half_width = 0.5 * self.gobo.width_m() * direction.z;
        let half_height = 0.5 * self.gobo.height_m() * direction.z;
        let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
        for beam in splitter.beams() {
            if beam.throughput() == 0.0 {
                continue;
            }
            let offset = beam.slope_offset();
            let x = pupil[0] * direction.z
                - self.gobo_distance_m * (direction.x - offset[0] * direction.z);
            let y = pupil[1] * direction.z
                - self.gobo_distance_m * (direction.y - offset[1] * direction.z);
            if !x.is_finite() || !y.is_finite() {
                return Err(SourceError::NumericalRange);
            }
            // Test bounded homogeneous coordinates before any grazing divide.
            if x < -half_width || x >= half_width || y < -half_height || y >= half_height {
                continue;
            }
            let mask = self
                .gobo
                .transmission_at(x / direction.z, y / direction.z)?;
            let value = beam.throughput() * mask;
            if mask > 0.0 && value == 0.0 {
                return Err(SourceError::NumericalRange);
            }
            let next = sum + value;
            correction += if sum >= value {
                (sum - next) + value
            } else {
                (value - next) + sum
            };
            sum = next;
        }
        let transmission = sum + correction;
        if transmission == 0.0 {
            return Ok(0.0);
        }
        let cosine_squared = direction.z * direction.z;
        let value = self.radiance_factor * transmission / cosine_squared / cosine_squared;
        if !value.is_finite() || value <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        Ok(value)
    }

    /// Finite upper bound on the represented radiance factor over an aperture
    /// with the supplied maximum absolute optical X/Y coordinates in metres.
    /// This bounds the source law, not floating-point or physical model error.
    pub fn maximum_luminance_factor(&self, pupil_abs_xy_m: [f64; 2]) -> Result<f64, SourceError> {
        if pupil_abs_xy_m
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(SourceError::InvalidDistribution);
        }
        if self.transmitted_fraction == 0.0 {
            return Ok(0.0);
        }
        let mut maximum_slope = [
            (pupil_abs_xy_m[0] * self.pupil_factor.abs() + 0.5 * self.gobo.width_m())
                / self.gobo_distance_m,
            (pupil_abs_xy_m[1] * self.pupil_factor.abs() + 0.5 * self.gobo.height_m())
                / self.gobo_distance_m,
        ];
        let mut maximum_transmission = self.gobo.data.maximum;
        if let Some(splitter) = &self.beam_splitter {
            for (maximum, offset) in maximum_slope.iter_mut().zip(splitter.maximum_abs_offset()) {
                *maximum = (*maximum + offset).next_up();
            }
            maximum_transmission *= splitter.transmitted_fraction();
        }
        if let Some(diffuser) = &self.diffuser {
            if diffuser.law.scattered_fraction() > 0.0 {
                for (maximum, spread) in maximum_slope
                    .iter_mut()
                    .zip(diffuser.law.half_slope_widths())
                {
                    *maximum = (*maximum + spread).next_up();
                }
            }
            maximum_transmission =
                diffuser::product(maximum_transmission, diffuser.law.transmission())?;
        }
        let inverse_cosine_squared =
            1.0 + maximum_slope[0] * maximum_slope[0] + maximum_slope[1] * maximum_slope[1];
        let value = self.radiance_factor
            * maximum_transmission
            * inverse_cosine_squared
            * inverse_cosine_squared;
        if !value.is_finite() || value <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        Ok(value)
    }

    /// Check one optical-frame aperture vertex against Z=0. The caller supplies
    /// the scale of original coordinates/translation used in its rigid transform
    /// so the rounding tolerance does not grow with unrelated scene geometry.
    pub(crate) fn validate_aperture_point(
        &self,
        point: Point3,
        roundoff_scale_m: f64,
    ) -> Result<(), SourceError> {
        if [point.x, point.y, point.z, roundoff_scale_m]
            .iter()
            .any(|value| !value.is_finite())
            || roundoff_scale_m < 0.0
        {
            return Err(SourceError::InvalidDistribution);
        }
        let scale = roundoff_scale_m
            .max(point.x.abs())
            .max(point.y.abs())
            .max(point.z.abs());
        if point.z.abs() > 128.0 * f64::EPSILON * scale {
            return Err(SourceError::InvalidDistribution);
        }
        Ok(())
    }

    /// Admit the front winding of an optical-frame planar aperture face.
    pub(crate) fn validate_aperture_normal(&self, normal: Vec3) -> Result<(), SourceError> {
        let normal = normalized_direction(normal)?;
        if normal.z <= 0.0
            || normal.x.abs() > 128.0 * f64::EPSILON
            || normal.y.abs() > 128.0 * f64::EPSILON
        {
            return Err(SourceError::InvalidDistribution);
        }
        Ok(())
    }
}

fn normalized_direction(direction: Vec3) -> Result<Vec3, SourceError> {
    if [direction.x, direction.y, direction.z]
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(SourceError::InvalidDirection);
    }
    let scale = direction
        .x
        .abs()
        .max(direction.y.abs())
        .max(direction.z.abs());
    if scale == 0.0 {
        return Err(SourceError::InvalidDirection);
    }
    let scaled = Vec3::new(
        direction.x / scale,
        direction.y / scale,
        direction.z / scale,
    );
    Ok(scaled.scale(1.0 / scaled.norm()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_gobo_coordinates_passivity_and_admission_are_explicit() {
        let gobo = GoboTransmission::try_new(2.0, 4.0, 2, 2, vec![0.0, 0.25, 0.5, 1.0]).unwrap();
        assert_eq!(gobo.mean_transmission(), 0.4375);
        assert_eq!(gobo.transmission_at(-0.5, -1.0), Ok(0.0));
        assert_eq!(gobo.transmission_at(0.5, -1.0), Ok(0.25));
        assert_eq!(gobo.transmission_at(-0.5, 1.0), Ok(0.5));
        assert_eq!(gobo.transmission_at(0.5, 1.0), Ok(1.0));
        assert_eq!(gobo.transmission_at(-1.0, 1.0), Ok(0.5));
        assert_eq!(gobo.transmission_at(1.0, 1.0), Ok(0.0));
        assert_eq!(gobo.transmission_at(0.5, 2.0), Ok(0.0));
        assert_eq!(gobo.transmission_at(1.0e300, 1.0), Ok(0.0));
        assert!(gobo.transmission_at(f64::NAN, 0.0).is_err());
        assert!(GoboTransmission::try_new(0.0, 1.0, 1, 1, vec![1.0]).is_err());
        assert!(GoboTransmission::try_new(1.0, 1.0, 1, 1, vec![1.01]).is_err());
        assert!(GoboTransmission::try_new(1.0, 1.0, 1, 1, vec![f64::NAN]).is_err());
        assert!(GoboTransmission::try_new(1.0, 1.0, 2, 1, vec![1.0]).is_err());
        assert!(matches!(
            GoboTransmission::try_new(1.0, 1.0, MAX_GOBO_TEXELS + 1, 1, Vec::new()),
            Err(SourceError::DistributionTooLarge)
        ));
        assert_eq!(
            GoboTransmission::try_new(1.0, 1.0, 1, 1, vec![0.0])
                .unwrap()
                .mean_transmission(),
            0.0
        );
    }

    #[test]
    fn g0_thin_lens_backprojection_normalizes_directions_and_rejects_dark_grazing_rays() {
        let gobo = GoboTransmission::try_new(0.2, 0.2, 2, 1, vec![0.0, 1.0]).unwrap();
        let projector =
            ThinLensProjector::try_new(RigidTransform::identity(), 1.0, 0.5, gobo.clone()).unwrap();
        assert_eq!(projector.focus_distance_m(), Some(1.0));
        let pupil = Point3::new(-0.1, 0.07, 0.0);
        let direction = Vec3::new(0.07, -0.06, 1.0);
        // This ray reaches screen (-0.03, 0.01) at q=1, hence g=(0.03,-0.01)
        // in the transmitting half. The exact law is s²/G * (1+u²+v²)².
        let expected = 25.0 * (1.0_f64 + 0.07 * 0.07 + 0.06 * 0.06).powi(2);
        let actual = projector
            .luminance_per_clear_flux_per_area(pupil, direction)
            .unwrap();
        assert!((actual - expected).abs() < 2.0e-14 * expected);
        let scaled = projector
            .luminance_per_clear_flux_per_area(pupil, direction.scale(1.0e100))
            .unwrap();
        assert!((scaled - actual).abs() < 2.0e-14 * actual);
        assert!(actual <= projector.maximum_luminance_factor([0.1, 0.07]).unwrap());
        for direction in [
            Vec3::new(0.13, -0.06, 1.0),
            Vec3::new(1.0, 0.0, f64::MIN_POSITIVE),
            Vec3::new(0.0, 0.0, -1.0),
        ] {
            assert_eq!(
                projector.luminance_per_clear_flux_per_area(pupil, direction),
                Ok(0.0)
            );
        }
        assert!(
            projector
                .luminance_per_clear_flux_per_area(pupil, Vec3::new(0.0, 0.0, 0.0))
                .is_err()
        );
        assert_eq!(
            ThinLensProjector::try_new(RigidTransform::identity(), 1.0, 1.0, gobo.clone())
                .unwrap()
                .focus_distance_m(),
            None
        );
        assert_eq!(
            ThinLensProjector::try_new(RigidTransform::identity(), 1.0, 2.0, gobo)
                .unwrap()
                .focus_distance_m(),
            Some(-2.0)
        );
    }

    #[test]
    fn g1_lens_refresh_preserves_assets_and_updates_native_focus_and_radiance() {
        let frame =
            RigidTransform::try_new([0.0, 0.0, 0.2_f64.sin(), 0.2_f64.cos()], [0.3, -0.2, 0.1])
                .unwrap();
        let gobo = GoboTransmission::try_new(0.2, 0.4, 2, 1, vec![0.25, 1.0]).unwrap();
        let splitter =
            ParaxialBeamSplitter::try_new(vec![ProjectorBeam::try_new([0.0; 2], 0.7).unwrap()])
                .unwrap();
        let base = ThinLensProjector::try_new(frame, 1.0, 0.5, gobo)
            .unwrap()
            .with_beam_splitter(splitter)
            .unwrap();
        let point = Point3::new(0.0, 0.0, 0.0);
        let direction = Vec3::new(-0.02, 0.01, 1.0);
        let original = base
            .luminance_per_clear_flux_per_area(point, direction)
            .unwrap();
        let identity = base.with_lens_distances(1.0, 0.5).unwrap();
        assert_eq!(
            identity
                .luminance_per_clear_flux_per_area(point, direction)
                .unwrap()
                .to_bits(),
            original.to_bits()
        );

        let updated = base.with_lens_distances(2.0, 1.0).unwrap();
        assert_eq!(updated.frame(), frame);
        assert_eq!(updated.focus_distance_m(), Some(2.0));
        assert!(Arc::ptr_eq(&base.gobo.data, &updated.gobo.data));
        assert!(std::ptr::eq(
            base.beam_splitter().unwrap().beams(),
            updated.beam_splitter().unwrap().beams(),
        ));
        assert_eq!(
            updated.transmitted_fraction().to_bits(),
            base.transmitted_fraction().to_bits()
        );
        // The ray remains in the fully transmitting cell. Its independent law
        // is s² / gobo_area * splitter_power / cos(theta)^4.
        let expected = 4.0 / (0.2 * 0.4) * 0.7 * (1.0_f64 + 0.02 * 0.02 + 0.01 * 0.01).powi(2);
        let actual = updated
            .luminance_per_clear_flux_per_area(point, direction)
            .unwrap();
        assert!((actual - expected).abs() <= 2.0e-14 * expected);
        assert!((actual - 4.0 * original).abs() <= 2.0e-14 * actual);
        let focused = base.with_lens_distances(1.0, 0.75).unwrap();
        assert!((focused.focus_distance_m().unwrap() - 3.0).abs() < 2.0e-14);
        assert_eq!(
            base.with_lens_distances(1.0, 1.0)
                .unwrap()
                .focus_distance_m(),
            None
        );
        assert_eq!(
            base.with_lens_distances(1.0, 2.0)
                .unwrap()
                .focus_distance_m(),
            Some(-2.0)
        );
        assert_eq!(base.gobo_distance_m(), 1.0);
        assert_eq!(base.focal_length_m(), 0.5);
        assert_eq!(base.focus_distance_m(), Some(1.0));
        assert_eq!(
            base.luminance_per_clear_flux_per_area(point, direction)
                .unwrap()
                .to_bits(),
            original.to_bits()
        );
    }

    #[test]
    fn g0_lens_refresh_refuses_invalid_and_unrepresentable_distances() {
        let base = ThinLensProjector::try_new(
            RigidTransform::identity(),
            1.0,
            0.5,
            GoboTransmission::try_new(0.2, 0.4, 1, 1, vec![1.0]).unwrap(),
        )
        .unwrap();
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                base.with_lens_distances(bad, 0.5),
                Err(SourceError::InvalidDistribution)
            ));
            assert!(matches!(
                base.with_lens_distances(1.0, bad),
                Err(SourceError::InvalidDistribution)
            ));
        }
        for (gobo_distance, focal_length) in [
            (f64::MAX, f64::MAX),
            (f64::from_bits(1), f64::from_bits(1)),
            (1.0, f64::from_bits(1)),
        ] {
            assert!(matches!(
                base.with_lens_distances(gobo_distance, focal_length),
                Err(SourceError::NumericalRange)
            ));
        }
        assert_eq!(base.gobo_distance_m(), 1.0);
        assert_eq!(base.focal_length_m(), 0.5);
    }

    #[test]
    fn g0_atomic_optical_state_replaces_and_removes_without_compounding_loss() {
        use crate::photometry::area::{EmissionSidedness, PhotometricAreaEmitter};
        use crate::photometry::source::{SPECTRAL_BINS, SourceSpectrum};

        let splitter = |power| {
            ParaxialBeamSplitter::try_new(vec![
                ProjectorBeam::try_new([0.1, -0.05], power).unwrap(),
            ])
            .unwrap()
        };
        let base = ThinLensProjector::try_new(
            RigidTransform::identity(),
            1.0,
            0.5,
            GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![0.5]).unwrap(),
        )
        .unwrap()
        .with_beam_splitter(splitter(0.5))
        .unwrap()
        .with_diffuser(ParaxialBoxDiffuser::try_new([0.1; 2], 0.75, 0.5).unwrap())
        .unwrap();
        let source = PhotometricAreaEmitter::try_new(
            7,
            320.0,
            0.5,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()
        .filtered([0.5; SPECTRAL_BINS])
        .unwrap()
        .with_projector(base.clone())
        .unwrap();
        let clear = source.clear_flux_lm();
        assert_eq!(base.transmitted_fraction(), 0.125);
        assert_eq!(source.flux_lm(), clear * 0.125);

        let updated = base
            .with_optical_state(
                1.0,
                0.5,
                GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![0.25]).unwrap(),
                Some(splitter(0.75)),
                Some(ParaxialBoxDiffuser::try_new([0.2; 2], 0.25, 0.8).unwrap()),
            )
            .unwrap();
        let updated_source = source.clone().with_projector(updated.clone()).unwrap();
        assert_eq!(updated.transmitted_fraction(), 0.25 * 0.75 * 0.8);
        assert_eq!(
            updated_source.flux_lm(),
            clear * updated.transmitted_fraction()
        );
        assert_eq!(updated_source.clear_flux_lm().to_bits(), clear.to_bits());
        assert_eq!(
            updated_source.spectrum().bins_w_per_lumen_nm(),
            source.spectrum().bins_w_per_lumen_nm()
        );
        let removed = updated
            .with_optical_state(1.0, 0.5, updated.gobo().clone(), None, None)
            .unwrap();
        assert!(removed.beam_splitter().is_none());
        assert!(removed.diffuser().is_none());
        assert_eq!(removed.transmitted_fraction(), 0.25);
        assert_eq!(
            source
                .clone()
                .with_projector(removed.clone())
                .unwrap()
                .flux_lm(),
            clear * 0.25
        );
        let opened = removed
            .with_optical_state(
                1.0,
                0.5,
                GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![1.0]).unwrap(),
                None,
                None,
            )
            .unwrap();
        assert_eq!(
            source.clone().with_projector(opened).unwrap().flux_lm(),
            clear
        );
        assert_eq!(source.clear_flux_lm(), clear);
        assert_eq!(source.flux_lm(), clear * 0.125);
        assert_eq!(base.transmitted_fraction(), 0.125);
    }

    #[test]
    fn g1_gobo_replacement_updates_native_support_and_radiance() {
        let frame =
            RigidTransform::try_new([0.0, 0.0, 0.2_f64.sin(), 0.2_f64.cos()], [0.3, -0.2, 0.1])
                .unwrap();
        let base = ThinLensProjector::try_new(
            frame,
            1.0,
            0.5,
            GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![1.0]).unwrap(),
        )
        .unwrap();
        let updated = base
            .with_optical_state(
                1.0,
                0.5,
                GoboTransmission::try_new(4.0, 2.0, 1, 1, vec![0.5]).unwrap(),
                None,
                None,
            )
            .unwrap();
        assert_eq!(updated.frame(), frame);
        assert_eq!(updated.gobo_distance_m(), base.gobo_distance_m());
        assert_eq!(updated.focal_length_m(), base.focal_length_m());
        assert_eq!(updated.focus_distance_m(), Some(1.0));
        assert_eq!(updated.transmitted_fraction(), 0.5);
        let pupil = Point3::new(0.0, 0.0, 0.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        assert_eq!(
            base.luminance_per_clear_flux_per_area(pupil, axis),
            Ok(0.25)
        );
        assert_eq!(
            updated.luminance_per_clear_flux_per_area(pupil, axis),
            Ok(0.5 / 8.0)
        );
        let direction = Vec3::new(1.5, 0.0, 1.0);
        assert_eq!(
            base.luminance_per_clear_flux_per_area(pupil, direction),
            Ok(0.0)
        );
        // g=(-1.5,0) is inside only the wider replacement. Its independent
        // law is s²/G * transmission * (1 + slope_x² + slope_y²)².
        let expected = 0.5 / 8.0 * (1.0_f64 + 1.5 * 1.5).powi(2);
        let actual = updated
            .luminance_per_clear_flux_per_area(pupil, direction)
            .unwrap();
        assert!((actual - expected).abs() <= 2.0e-14 * expected);
        assert!(actual <= updated.maximum_luminance_factor([0.0; 2]).unwrap());
        assert_eq!(
            updated.luminance_per_clear_flux_per_area(pupil, Vec3::new(2.1, 0.0, 1.0)),
            Ok(0.0)
        );
    }

    #[test]
    fn g0_atomic_optics_admit_only_final_preparation_and_representable_power() {
        let gobo = GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![1.0]).unwrap();
        let base = ThinLensProjector::try_new(RigidTransform::identity(), 1.0, 0.5, gobo.clone())
            .unwrap()
            .with_diffuser(ParaxialBoxDiffuser::try_new([1.0e-200; 2], 1.0, 1.0).unwrap())
            .unwrap();
        // Retaining the old frost would underflow its metric width. The final
        // state removes it, so no obsolete intermediate preparation is needed.
        assert!(matches!(
            base.with_lens_distances(1.0e-150, 0.5e-150),
            Err(SourceError::NumericalRange)
        ));
        let removed = base
            .with_optical_state(1.0e-150, 0.5e-150, gobo.clone(), None, None)
            .unwrap();
        assert!(removed.diffuser().is_none());
        assert_eq!(removed.transmitted_fraction(), 1.0);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                base.with_optical_state(bad, 0.5, gobo.clone(), None, None),
                Err(SourceError::InvalidDistribution)
            ));
        }
        let faint_gobo =
            GoboTransmission::try_new(2.0, 2.0, 1, 1, vec![f64::MIN_POSITIVE]).unwrap();
        let faint_splitter = ParaxialBeamSplitter::try_new(vec![
            ProjectorBeam::try_new([0.0; 2], f64::MIN_POSITIVE).unwrap(),
        ])
        .unwrap();
        assert!(matches!(
            base.with_optical_state(
                1.0,
                0.5,
                faint_gobo.clone(),
                Some(faint_splitter.clone()),
                None
            ),
            Err(SourceError::NumericalRange)
        ));
        let opaque = base
            .with_optical_state(
                1.0,
                0.5,
                faint_gobo,
                Some(faint_splitter),
                Some(ParaxialBoxDiffuser::try_new([0.1; 2], 1.0, 0.0).unwrap()),
            )
            .unwrap();
        assert_eq!(opaque.transmitted_fraction(), 0.0);
        assert_eq!(
            opaque.luminance_per_clear_flux_per_area(
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0)
            ),
            Ok(0.0)
        );
        assert_eq!(base.gobo_distance_m(), 1.0);
        assert!(base.diffuser().is_some());
        assert_eq!(base.transmitted_fraction(), 1.0);
    }
}
