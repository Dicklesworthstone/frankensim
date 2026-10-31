//! Calibrated point emitters in SI metres, radians, lumens and candela.
//!
//! Angular distributions integrate to one over solid angle. Tabulated data is
//! interpolated in `cos(theta)` and azimuth, whose integration measure is
//! `d(cos(theta)) d(phi)`; normalization therefore integrates the represented
//! distribution exactly, apart from floating-point rounding. Tables are plain
//! Rust data, so file-format conversion remains the caller's responsibility.
//!
//! Spectra use 80-bin, 380--780 nm quadrature and the published CIE 1931
//! 2-degree observer, linearly interpolated at bin midpoints (see [`super::observer`]).
//! The factor 683 has units lm/W. Normalization is with respect to that declared
//! quadrature, not certified CIE metrology. Narrower spectral features,
//! omitted wavelengths and measured-source uncertainty are not bounded here.
//! Angular interpolation is a declared
//! source model, not a bound on an unknown manufacturer's intensity function.
//! Photometric unit reference: <https://www.bipm.org/en/si-base-units/candela>.

use core::f64::consts::{FRAC_PI_2, PI, TAU};
use core::fmt;
use std::sync::Arc;

use fs_geom::{Point3, Vec3};
use fs_math::det;

use crate::instances::RigidTransform;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, QUAD_BINS};

mod ray_partition;
pub mod validation;

/// Spectral bins shared with the existing spectral renderer.
pub const SPECTRAL_BINS: usize = QUAD_BINS;
/// Width of each spectral bin, in nanometres.
pub const BIN_WIDTH_NM: f64 = (LAMBDA_MAX - LAMBDA_MIN) / SPECTRAL_BINS as f64;
/// Photometric conversion constant, in lumens per watt.
pub const LUMINOUS_EFFICACY_LM_W: f64 = 683.0;
/// Maximum admitted nodal entries of a tabulated angular distribution.
pub const MAX_ANGULAR_ENTRIES: usize = 1_048_576;

/// Midpoint wavelength of spectral bin `bin`, in nanometres.
///
/// The caller supplies `bin < SPECTRAL_BINS`.
#[must_use]
pub fn wavelength_nm(bin: usize) -> f64 {
    LAMBDA_MIN + (bin as f64 + 0.5) * BIN_WIDTH_NM
}

/// Photopic quadrature weight `y_bar(lambda) * delta_lambda`, in nanometres.
/// Multiply by 683 and spectral radiance/irradiance to obtain cd/m²/lux.
#[must_use]
pub fn photopic_weight_nm(bin: usize) -> f64 {
    super::observer::xyz_at_bin(bin)[1] * BIN_WIDTH_NM
}

/// A finite-input or representability refusal for a calibrated source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceError {
    /// A direction is zero or contains a non-finite component.
    InvalidDirection,
    /// An angular law, grid, seam, pole, or sample is invalid.
    InvalidDistribution,
    /// A table exceeds the bounded admission size.
    DistributionTooLarge,
    /// A spectrum is non-finite, negative, or has zero photopic response.
    InvalidSpectrum,
    /// Filter transmittance must be finite and in `[0,1]` in every bin.
    InvalidFilter,
    /// Flux must be finite and nonnegative; dimmer must be in `[0,1]`.
    InvalidFlux,
    /// Evaluation or normalization would produce a non-finite value.
    NumericalRange,
}

impl fmt::Display for SourceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "photometric source: {self:?}")
    }
}

impl std::error::Error for SourceError {}

/// A normalized intensity distribution about the local **+Z** optical axis.
/// Azimuth increases from +X toward +Y in a right-handed fixture frame.
#[derive(Debug, Clone)]
pub struct AngularDistribution {
    law: AngularLaw,
    peak_per_sr: f64,
}

#[derive(Debug, Clone)]
enum AngularLaw {
    Isotropic,
    Cone {
        cosine_cutoff: f64,
    },
    CosinePower {
        exponent: f64,
        cosine_cutoff: f64,
        e_fold_cosine: f64,
    },
    Table(Arc<AngularTable>),
}

#[derive(Debug)]
struct AngularTable {
    // Descending cos(theta), including +1 and -1.
    cosine: Vec<f64>,
    phi: Vec<f64>,
    // Polar-major rows, azimuth-minor columns, already in sr^-1.
    values: Vec<f64>,
}

impl AngularDistribution {
    /// Uniform intensity over the full sphere: `1/(4*pi)` per steradian.
    #[must_use]
    pub fn isotropic() -> Self {
        Self {
            law: AngularLaw::Isotropic,
            peak_per_sr: 1.0 / (4.0 * PI),
        }
    }

    /// Uniform intensity inside a circular cone of half-angle `(0, pi]`.
    /// Intensity is zero outside the cone, with the boundary included.
    pub fn uniform_cone(half_angle_rad: f64) -> Result<Self, SourceError> {
        if !half_angle_rad.is_finite() || !(0.0 < half_angle_rad && half_angle_rad <= PI) {
            return Err(SourceError::InvalidDistribution);
        }
        let cosine_cutoff = if half_angle_rad == PI {
            -1.0
        } else {
            det::cos(half_angle_rad)
        };
        if cosine_cutoff >= 1.0 {
            return Err(SourceError::NumericalRange);
        }
        // Integrate the admitted cosine boundary itself. A different
        // trigonometric identity can disagree after rounding for a narrow cone.
        let peak_per_sr = 1.0 / (TAU * (1.0 - cosine_cutoff));
        if !peak_per_sr.is_finite() {
            return Err(SourceError::NumericalRange);
        }
        Ok(Self {
            law: AngularLaw::Cone { cosine_cutoff },
            peak_per_sr,
        })
    }

    /// `cos(theta)^exponent` within a forward cone, normalized to unit flux.
    /// The exponent is nonnegative and half-angle is in `(0, pi/2]`.
    /// Exponent one over a hemisphere is a Lambertian angular source.
    /// A positive exponent whose natural e-fold cosine rounds to one returns
    /// `NumericalRange`: its angular width is not representable for quadrature.
    pub fn cosine_power(exponent: f64, half_angle_rad: f64) -> Result<Self, SourceError> {
        if !exponent.is_finite()
            || exponent < 0.0
            || !half_angle_rad.is_finite()
            || !(0.0 < half_angle_rad && half_angle_rad <= FRAC_PI_2)
        {
            return Err(SourceError::InvalidDistribution);
        }
        let cosine_cutoff = if half_angle_rad == FRAC_PI_2 {
            0.0
        } else {
            det::cos(half_angle_rad)
        };
        if cosine_cutoff >= 1.0 {
            return Err(SourceError::NumericalRange);
        }
        let e_fold_cosine = if exponent == 0.0 {
            0.0
        } else {
            det::exp(-1.0 / exponent)
        };
        if !e_fold_cosine.is_finite() || e_fold_cosine >= 1.0 {
            return Err(SourceError::NumericalRange);
        }
        let exponent_plus_one = exponent + 1.0;
        let numerator = if cosine_cutoff == 0.0 {
            1.0
        } else {
            -det::expm1(exponent_plus_one * det::ln(cosine_cutoff))
        };
        let peak_per_sr = exponent_plus_one / (TAU * numerator);
        if !peak_per_sr.is_finite() || peak_per_sr <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        Ok(Self {
            law: AngularLaw::CosinePower {
                exponent,
                cosine_cutoff,
                e_fold_cosine,
            },
            peak_per_sr,
        })
    }

    /// An axisymmetric beam/field profile with explicit tail cutoff.
    ///
    /// All three arguments are **full** opening angles in radians. Relative
    /// intensity is 1 on axis, 0.5 at half the beam angle, 0.1 at half the field
    /// angle, and 0 at half the cutoff angle. Between knots it is linear in
    /// cos(theta), with exactly integrated flux normalization. Angles satisfy
    /// `0 < beam < field < cutoff <= 2*pi`.
    ///
    /// The cutoff is deliberately supplied: beam and field angles alone do
    /// not determine the unmeasured tails or total flux. This interpolation
    /// model is a declared approximation of the fixture, not recovered optics.
    pub fn beam_field(
        beam_full_angle_rad: f64,
        field_full_angle_rad: f64,
        cutoff_full_angle_rad: f64,
    ) -> Result<Self, SourceError> {
        if !beam_full_angle_rad.is_finite()
            || !field_full_angle_rad.is_finite()
            || !cutoff_full_angle_rad.is_finite()
            || !(0.0 < beam_full_angle_rad
                && beam_full_angle_rad < field_full_angle_rad
                && field_full_angle_rad < cutoff_full_angle_rad
                && cutoff_full_angle_rad <= TAU)
        {
            return Err(SourceError::InvalidDistribution);
        }
        let mut theta = vec![
            0.0,
            0.5 * beam_full_angle_rad,
            0.5 * field_full_angle_rad,
            0.5 * cutoff_full_angle_rad,
        ];
        let mut values = vec![1.0, 1.0, 0.5, 0.5, 0.1, 0.1, 0.0, 0.0];
        if cutoff_full_angle_rad < TAU {
            theta.push(PI);
            values.extend_from_slice(&[0.0, 0.0]);
        }
        Self::tabulated(theta, vec![0.0, TAU], values)
    }

    /// A full-sphere intensity table, normalized by its represented flux.
    ///
    /// `theta_rad` strictly increases from 0 to pi; `phi_rad` from 0 to 2*pi.
    /// `values[theta_index * phi_rad.len() + phi_index]` is a nonnegative
    /// relative intensity. The repeated 0/2*pi seam must agree, and each pole
    /// must have one value independent of azimuth. Interpolation is bilinear in
    /// cos(theta) and phi, rather than in theta; this is the explicit model.
    /// Each cell's integrated flux is its four-corner mean times solid angle.
    /// The arrays are shared by clones, so a rig can reuse one fixture profile.
    pub fn tabulated(
        mut theta_rad: Vec<f64>,
        phi_rad: Vec<f64>,
        mut values: Vec<f64>,
    ) -> Result<Self, SourceError> {
        let count = theta_rad
            .len()
            .checked_mul(phi_rad.len())
            .ok_or(SourceError::DistributionTooLarge)?;
        if count > MAX_ANGULAR_ENTRIES {
            return Err(SourceError::DistributionTooLarge);
        }
        if theta_rad.len() < 2
            || phi_rad.len() < 2
            || values.len() != count
            || theta_rad[0] != 0.0
            || theta_rad[theta_rad.len() - 1] != PI
            || phi_rad[0] != 0.0
            || phi_rad[phi_rad.len() - 1] != TAU
            || theta_rad.iter().chain(&phi_rad).any(|v| !v.is_finite())
            || theta_rad.windows(2).any(|w| w[0] >= w[1])
            || phi_rad.windows(2).any(|w| w[0] >= w[1])
            || values.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(SourceError::InvalidDistribution);
        }
        let columns = phi_rad.len();
        for row in values.chunks_exact(columns) {
            if row[0] != row[columns - 1] {
                return Err(SourceError::InvalidDistribution);
            }
        }
        for pole in [0, theta_rad.len() - 1] {
            let row = &values[pole * columns..(pole + 1) * columns];
            if row.iter().any(|v| *v != row[0]) {
                return Err(SourceError::InvalidDistribution);
            }
        }
        // Scale first so normalization handles very large or tiny tabulated
        // candela values without overflowing a sum or discarding their ratios.
        let peak = values.iter().copied().fold(0.0_f64, f64::max);
        if peak == 0.0 {
            return Err(SourceError::InvalidDistribution);
        }
        for value in &mut values {
            *value /= peak;
        }
        let last = theta_rad.len() - 1;
        for (index, theta) in theta_rad.iter_mut().enumerate() {
            *theta = if index == 0 {
                1.0
            } else if index == last {
                -1.0
            } else {
                det::cos(*theta)
            };
        }
        if theta_rad.windows(2).any(|w| w[0] <= w[1]) {
            return Err(SourceError::NumericalRange);
        }
        let mut integral = 0.0;
        for (row, mu) in theta_rad.windows(2).enumerate() {
            for (column, phi) in phi_rad.windows(2).enumerate() {
                let i = row * columns + column;
                let mean = 0.25
                    * (values[i] + values[i + 1] + values[i + columns] + values[i + columns + 1]);
                integral += mean * (mu[0] - mu[1]) * (phi[1] - phi[0]);
            }
        }
        let peak_per_sr = 1.0 / integral;
        if !peak_per_sr.is_finite() || peak_per_sr <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        for value in &mut values {
            *value *= peak_per_sr;
        }
        Ok(Self {
            law: AngularLaw::Table(Arc::new(AngularTable {
                cosine: theta_rad,
                phi: phi_rad,
                values,
            })),
            peak_per_sr,
        })
    }

    /// Largest intensity per unit flux, in inverse steradians.
    #[must_use]
    pub const fn peak_per_steradian(&self) -> f64 {
        self.peak_per_sr
    }

    /// Intensity per unit luminous flux, in inverse steradians.
    /// Finite nonzero directions are normalized before evaluation.
    pub fn intensity_per_steradian(&self, local_direction: Vec3) -> Result<f64, SourceError> {
        let direction = unit(local_direction)?;
        let mu = direction.z.clamp(-1.0, 1.0);
        Ok(match &self.law {
            AngularLaw::Isotropic => self.peak_per_sr,
            AngularLaw::Cone { cosine_cutoff } => {
                if mu >= *cosine_cutoff {
                    self.peak_per_sr
                } else {
                    0.0
                }
            }
            AngularLaw::CosinePower {
                exponent,
                cosine_cutoff,
                ..
            } => {
                if mu < *cosine_cutoff {
                    0.0
                } else if *exponent == 0.0 {
                    self.peak_per_sr
                } else {
                    self.peak_per_sr * det::pow(mu.max(0.0), *exponent)
                }
            }
            AngularLaw::Table(table) => {
                let raw_phi = det::atan2(direction.y, direction.x);
                let phi = if raw_phi < 0.0 {
                    raw_phi + TAU
                } else {
                    raw_phi
                };
                let row = table
                    .cosine
                    .partition_point(|value| *value > mu)
                    .saturating_sub(1)
                    .min(table.cosine.len() - 2);
                let column = table
                    .phi
                    .partition_point(|value| *value <= phi)
                    .saturating_sub(1)
                    .min(table.phi.len() - 2);
                let t = ((table.cosine[row] - mu) / (table.cosine[row] - table.cosine[row + 1]))
                    .clamp(0.0, 1.0);
                let u = ((phi - table.phi[column]) / (table.phi[column + 1] - table.phi[column]))
                    .clamp(0.0, 1.0);
                let i = row * table.phi.len() + column;
                let upper = (1.0 - u) * table.values[i] + u * table.values[i + 1];
                let lower = (1.0 - u) * table.values[i + table.phi.len()]
                    + u * table.values[i + table.phi.len() + 1];
                (1.0 - t) * upper + t * lower
            }
        })
    }
}

/// A nonnegative visible spectrum normalized to one photopic lumen.
/// Values are W/(lm nm), constant inside each of the 80 spectral bins.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceSpectrum {
    bins: [f64; SPECTRAL_BINS],
}

impl SourceSpectrum {
    /// Equal radiant power at every visible wavelength, per emitted lumen.
    #[must_use]
    pub fn equal_energy() -> Self {
        Self::try_from_bins([1.0; SPECTRAL_BINS]).expect("finite positive observer integral")
    }

    /// Normalize arbitrary nonnegative relative spectral samples to one lumen
    /// using the calibrated transport's explicit CIE photopic quadrature.
    pub fn try_from_bins(mut bins: [f64; SPECTRAL_BINS]) -> Result<Self, SourceError> {
        if bins.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(SourceError::InvalidSpectrum);
        }
        let peak = bins.iter().copied().fold(0.0_f64, f64::max);
        if peak == 0.0 {
            return Err(SourceError::InvalidSpectrum);
        }
        let mut integral = 0.0;
        for (index, value) in bins.iter_mut().enumerate() {
            *value /= peak;
            integral += *value * photopic_weight_nm(index);
        }
        let normalization = 1.0 / (LUMINOUS_EFFICACY_LM_W * integral);
        if !normalization.is_finite() || normalization <= 0.0 {
            return Err(SourceError::NumericalRange);
        }
        for value in &mut bins {
            *value *= normalization;
        }
        Ok(Self { bins })
    }

    /// Normalized 80-bin spectrum in W/(lm nm), without copying its buffer.
    #[must_use]
    pub const fn bins_w_per_lumen_nm(&self) -> &[f64; SPECTRAL_BINS] {
        &self.bins
    }

    /// Spectrum value in W/(lm nm); zero outside the admitted visible band.
    /// The upper endpoint belongs to the last bin. A non-finite wavelength
    /// returns zero, since it does not identify any admitted spectral bin.
    #[must_use]
    pub fn value_w_per_lumen_nm(&self, lambda_nm: f64) -> f64 {
        if !lambda_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&lambda_nm) {
            return 0.0;
        }
        let bin = ((lambda_nm - LAMBDA_MIN) / BIN_WIDTH_NM) as usize;
        self.bins[bin.min(SPECTRAL_BINS - 1)]
    }

    /// Shared passive filtering for calibrated point and finite sources.
    /// Return the surviving normalized spectrum and its photopic transmission;
    /// callers must apply that fraction to their luminous-flux calibration.
    /// An opaque result retains a valid spectrum but has zero transmission.
    pub(crate) fn filtered(
        &self,
        transmittance: [f64; SPECTRAL_BINS],
    ) -> Result<(Self, f64), SourceError> {
        if transmittance
            .iter()
            .any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
        {
            return Err(SourceError::InvalidFilter);
        }
        if transmittance.iter().all(|value| *value == transmittance[0]) {
            return Ok((self.clone(), transmittance[0]));
        }
        let mut surviving = [0.0; SPECTRAL_BINS];
        let mut transmission = 0.0;
        for (index, value) in surviving.iter_mut().enumerate() {
            *value = self.bins[index] * transmittance[index];
            transmission += *value * LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(index);
        }
        let spectrum = if transmission > 0.0 {
            Self::try_from_bins(surviving)?
        } else {
            self.clone()
        };
        Ok((spectrum, transmission.min(1.0)))
    }
}

/// A far-field point fixture or one emitter of a multi-emitter fixture.
/// The rigid transform places the optical origin in metres and rotates +Z
/// toward the beam direction. No source size or near-field optics is inferred.
#[derive(Debug, Clone)]
pub struct PhotometricSource {
    transform: RigidTransform,
    world_to_local_rotation: RigidTransform,
    emitted_flux_lm: f64,
    distribution: AngularDistribution,
    spectrum: SourceSpectrum,
}

impl PhotometricSource {
    /// Admit a fixture at its supplied transform and emitted luminous flux.
    /// `lumens` is full-output flux, `dimmer` is a linear factor in `[0,1]`.
    pub fn try_new(
        transform: RigidTransform,
        lumens: f64,
        dimmer: f64,
        distribution: AngularDistribution,
        spectrum: SourceSpectrum,
    ) -> Result<Self, SourceError> {
        if !lumens.is_finite()
            || lumens < 0.0
            || !dimmer.is_finite()
            || !(0.0..=1.0).contains(&dimmer)
        {
            return Err(SourceError::InvalidFlux);
        }
        let emitted_flux_lm = lumens * dimmer;
        let peak_cd = emitted_flux_lm * distribution.peak_per_steradian();
        if !peak_cd.is_finite() || spectrum.bins.iter().any(|v| !(peak_cd * *v).is_finite()) {
            return Err(SourceError::NumericalRange);
        }
        let [x, y, z, w] = transform.rotation_xyzw();
        let world_to_local_rotation = RigidTransform::try_new([-x, -y, -z, w], [0.0; 3])
            .map_err(|_| SourceError::NumericalRange)?;
        Ok(Self {
            transform,
            world_to_local_rotation,
            emitted_flux_lm,
            distribution,
            spectrum,
        })
    }

    /// Source position in world metres.
    #[must_use]
    pub fn position_m(&self) -> Point3 {
        let [x, y, z] = self.transform.translation_m();
        Point3::new(x, y, z)
    }

    /// Original proper-rigid local-to-world fixture placement.
    #[must_use]
    pub const fn transform(&self) -> RigidTransform {
        self.transform
    }

    /// Emitted flux in lumens after dimming and any applied filters.
    #[must_use]
    pub const fn flux_lm(&self) -> f64 {
        self.emitted_flux_lm
    }

    /// The normalized emitted spectrum, in W/(lm nm).
    #[must_use]
    pub const fn spectrum(&self) -> &SourceSpectrum {
        &self.spectrum
    }

    /// Luminous intensity in candela along a world-space outgoing direction.
    pub fn intensity_cd(&self, world_direction: Vec3) -> Result<f64, SourceError> {
        // Normalize before rotating to avoid overflow for finite large inputs.
        let local = self
            .world_to_local_rotation
            .transform_vector(unit(world_direction)?);
        Ok(self.emitted_flux_lm * self.distribution.intensity_per_steradian(local)?)
    }

    /// Split a camera segment at source-distance scales and represented
    /// source-profile/phase boundaries and smooth power-profile scales.
    /// Knots are sorted, exactly deduplicated fractions in `(0,1)`, shared
    /// across fixtures. They guide quadrature but do not certify geometry,
    /// resolve shadow boundaries, or bound phase/smooth-profile variation.
    pub(crate) fn append_scattering_breakpoints(
        &self,
        ray: &crate::charts::Ray,
        distance_m: f64,
        phase: &super::medium::PhaseFunction,
        max_panels: u32,
        knots: &mut Vec<f64>,
        cx: &fs_exec::Cx<'_>,
    ) -> Result<(), super::PhotometryError> {
        use super::{PhotometryError, medium::PhaseFunction};
        use ray_partition::{
            KnotError, azimuth_crossing, cone_crossings, cosine_power_levels, distance_breakpoints,
            insert_knot, phase_crossing,
        };

        cx.checkpoint()?;
        if self.emitted_flux_lm == 0.0 || distance_m == 0.0 {
            return Ok(());
        }
        let limit = max_panels
            .checked_sub(1)
            .ok_or(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels })?
            as usize;
        if knots.len() > limit {
            return Err(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels });
        }
        let offset = ray.origin.delta_from(self.position_m());
        let scale = offset
            .x
            .abs()
            .max(offset.y.abs())
            .max(offset.z.abs())
            .max(distance_m);
        if !scale.is_finite() || scale <= 0.0 || !distance_m.is_finite() || distance_m < 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        // Rotate after scaling: both the displacement and offset remain
        // representable even for very large metric input coordinates.
        let origin = self.world_to_local_rotation.transform_vector(Vec3::new(
            offset.x / scale,
            offset.y / scale,
            offset.z / scale,
        ));
        let delta = self
            .world_to_local_rotation
            .transform_vector(ray.dir.scale(distance_m / scale));
        let origin = [origin.x, origin.y, origin.z];
        let delta = [delta.x, delta.y, delta.z];
        if origin.iter().chain(&delta).any(|value| !value.is_finite()) {
            return Err(PhotometryError::NonFiniteResult);
        }
        let mut insert = |value: f64| {
            insert_knot(knots, value, limit).map_err(|error| match error {
                KnotError::Budget => {
                    PhotometryError::ScatteringPartitionBudgetExceeded { max_panels }
                }
                KnotError::Allocation => PhotometryError::AllocationRefused,
            })
        };
        // Even an isotropic source has a potentially narrow inverse-square
        // peak. Resolve its metric scales independently of fixture rotation.
        // Off-segment collinear sources use their endpoint distance, not a
        // zero line impact parameter. Scene admission already excludes every
        // active source inside the ray's singularity radius.
        let camera_to_source = [-offset.x / scale, -offset.y / scale, -offset.z / scale];
        let direction = [ray.dir.x, ray.dir.y, ray.dir.z];
        let distance_knots = distance_breakpoints(camera_to_source, direction, distance_m / scale)
            .ok_or(PhotometryError::NonFiniteResult)?;
        for knot in distance_knots {
            cx.checkpoint()?;
            insert(knot)?;
        }
        match &self.distribution.law {
            AngularLaw::Isotropic => {}
            AngularLaw::Cone { cosine_cutoff } if *cosine_cutoff == -1.0 => {}
            AngularLaw::Cone { cosine_cutoff } => {
                for root in cone_crossings(origin, delta, *cosine_cutoff)
                    .into_iter()
                    .flatten()
                {
                    insert(root)?;
                }
            }
            AngularLaw::CosinePower {
                cosine_cutoff,
                e_fold_cosine,
                ..
            } => {
                let levels = cosine_power_levels(*e_fold_cosine, *cosine_cutoff)
                    .ok_or(PhotometryError::NonFiniteResult)?;
                for cosine in std::iter::once(*cosine_cutoff).chain(levels) {
                    cx.checkpoint()?;
                    for root in cone_crossings(origin, delta, cosine).into_iter().flatten() {
                        insert(root)?;
                    }
                }
            }
            AngularLaw::Table(table) => {
                for cosine in &table.cosine[1..table.cosine.len() - 1] {
                    cx.checkpoint()?;
                    for root in cone_crossings(origin, delta, *cosine).into_iter().flatten() {
                        insert(root)?;
                    }
                }
                // A two-column seam table is axisymmetric. Otherwise include
                // the 0/2pi seam once and retain only each actual azimuth ray.
                if table.phi.len() > 2 {
                    for phi in &table.phi[..table.phi.len() - 1] {
                        cx.checkpoint()?;
                        let (sine, cosine) = det::sin_cos(*phi);
                        if let Some(root) = azimuth_crossing(origin, delta, sine, cosine) {
                            insert(root)?;
                        }
                    }
                }
            }
        }
        if let PhaseFunction::Tabulated(table) = phase {
            // The phase angle belongs to the source-camera geometry, not the
            // fixture's optical frame. Scale before forming dot/cross products.
            let phase_knots = table.knots_per_sr();
            for &(cosine, _) in &phase_knots[1..phase_knots.len() - 1] {
                cx.checkpoint()?;
                if let Some(root) =
                    phase_crossing(camera_to_source, direction, distance_m / scale, cosine)
                {
                    insert(root)?;
                }
            }
        }
        cx.checkpoint()?;
        Ok(())
    }

    /// Spectral radiant intensity in W/(sr nm). Use `I*cos(theta)/r²` for
    /// receiver spectral irradiance, with distances in metres.
    pub fn spectral_intensity_w_sr_nm(
        &self,
        world_direction: Vec3,
        lambda_nm: f64,
    ) -> Result<f64, SourceError> {
        Ok(self.intensity_cd(world_direction)? * self.spectrum.value_w_per_lumen_nm(lambda_nm))
    }

    /// Apply a passive spectral filter without replacing lost flux.
    ///
    /// Surviving spectral shape is normalized, and emitted lumens decrease by
    /// the photopic transmission of the original shape. A fully opaque filter
    /// returns a zero-flux source. Repeated filters multiply transmission.
    pub fn filtered(&self, transmittance: [f64; SPECTRAL_BINS]) -> Result<Self, SourceError> {
        let (spectrum, transmission) = self.spectrum.filtered(transmittance)?;
        let mut result = self.clone();
        result.emitted_flux_lm *= transmission;
        result.spectrum = spectrum;
        Ok(result)
    }
}

fn unit(direction: Vec3) -> Result<Vec3, SourceError> {
    if !direction.x.is_finite() || !direction.y.is_finite() || !direction.z.is_finite() {
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
    Ok(scaled.scale(1.0 / det::sqrt(scaled.dot(scaled))))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual}, expected={expected}"
        );
    }

    #[test]
    fn analytic_distributions_have_known_candela_and_flux() {
        let isotropic = AngularDistribution::isotropic();
        close(
            isotropic
                .intensity_per_steradian(Vec3::new(1.0, 2.0, 3.0))
                .unwrap()
                * 4.0
                * PI,
            1.0,
            1e-15,
        );
        let lambertian = AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap();
        close(
            lambertian
                .intensity_per_steradian(Vec3::new(0.0, 0.0, 1.0))
                .unwrap(),
            1.0 / PI,
            1e-15,
        );
        close(
            lambertian
                .intensity_per_steradian(Vec3::new(3.0_f64.sqrt(), 0.0, 1.0))
                .unwrap(),
            0.5 / PI,
            1e-15,
        );
        assert_eq!(
            lambertian
                .intensity_per_steradian(Vec3::new(0.0, 0.0, -1.0))
                .unwrap(),
            0.0
        );
        let cone = AngularDistribution::uniform_cone(PI / 3.0).unwrap();
        close(cone.peak_per_steradian(), 1.0 / PI, 1e-15);
        assert_eq!(
            cone.intensity_per_steradian(Vec3::new(1.0, 0.0, 0.0))
                .unwrap(),
            0.0
        );
    }

    #[test]
    fn cosine_power_refuses_an_unrepresentable_natural_width() {
        assert!(matches!(
            AngularDistribution::cosine_power(1e20, FRAC_PI_2),
            Err(SourceError::NumericalRange)
        ));
        for exponent in [0.0, f64::MIN_POSITIVE, 1.0, 20_000.0] {
            assert!(AngularDistribution::cosine_power(exponent, FRAC_PI_2).is_ok());
        }
    }

    #[test]
    fn tabulated_linear_cosine_law_integrates_to_unit_flux() {
        let table = AngularDistribution::tabulated(
            vec![0.0, FRAC_PI_2, PI],
            vec![0.0, TAU],
            vec![0.0, 0.0, 1.0, 1.0, 2.0, 2.0],
        )
        .unwrap();
        let mut integral = 0.0;
        for i in 0..1024 {
            let mu = -1.0 + 2.0 * (f64::from(i) + 0.5) / 1024.0;
            let value = table
                .intensity_per_steradian(Vec3::new((1.0 - mu * mu).sqrt(), 0.0, mu))
                .unwrap();
            close(value, (1.0 - mu) / (4.0 * PI), 2e-15);
            integral += value * (4.0 * PI / 1024.0);
        }
        close(integral, 1.0, 2e-14);
    }

    #[test]
    fn beam_and_field_thresholds_survive_flux_normalization() {
        let profile = AngularDistribution::beam_field(0.2, 0.4, 0.6).unwrap();
        let peak = profile.peak_per_steradian();
        for (theta, fraction) in [(0.0, 1.0), (0.1, 0.5), (0.2, 0.1), (0.4, 0.0)] {
            let (sine, cosine) = det::sin_cos(theta);
            let intensity = profile
                .intensity_per_steradian(Vec3::new(sine, 0.0, cosine))
                .unwrap();
            close(intensity / peak, fraction, 3e-14);
        }
        let expected_integral = TAU
            * (0.75 * (1.0 - det::cos(0.1))
                + 0.3 * (det::cos(0.1) - det::cos(0.2))
                + 0.05 * (det::cos(0.2) - det::cos(0.3)));
        close(peak * expected_integral, 1.0, 2e-15);
        assert!(AngularDistribution::beam_field(0.4, 0.2, 0.6).is_err());
    }

    #[test]
    fn azimuth_table_wraps_and_normalizes_in_solid_angle() {
        let table = AngularDistribution::tabulated(
            vec![0.0, FRAC_PI_2, PI],
            vec![0.0, PI, TAU],
            vec![1.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0],
        )
        .unwrap();
        close(
            table
                .intensity_per_steradian(Vec3::new(-1.0, 0.0, 0.0))
                .unwrap(),
            3.0 / (6.0 * PI),
            1e-15,
        );
        close(
            table
                .intensity_per_steradian(Vec3::new(1.0, 0.0, 0.0))
                .unwrap(),
            1.0 / (6.0 * PI),
            1e-15,
        );
        let (mut integral, mut compensation) = (0.0, 0.0);
        for i in 0..64 {
            let mu = -1.0 + 2.0 * (f64::from(i) + 0.5) / 64.0;
            let radial = (1.0 - mu * mu).sqrt();
            for j in 0..128 {
                let phi = TAU * (f64::from(j) + 0.5) / 128.0;
                let (sine, cosine) = det::sin_cos(phi);
                let value = table
                    .intensity_per_steradian(Vec3::new(radial * cosine, radial * sine, mu))
                    .unwrap()
                    * 4.0
                    * PI
                    / (64.0 * 128.0);
                // The independent 8192-cell quadrature should measure the
                // angular law, not lose its tolerance to naive summation.
                let adjusted = value - compensation;
                let next = integral + adjusted;
                compensation = (next - integral) - adjusted;
                integral = next;
            }
        }
        close(integral, 1.0, 3e-14);
    }

    #[test]
    fn transform_dimmer_and_spectrum_preserve_declared_flux() {
        let q = 0.5_f64.sqrt();
        let transform = RigidTransform::try_new([0.0, q, 0.0, q], [2.0, 3.0, 4.0]).unwrap();
        let source = PhotometricSource::try_new(
            transform,
            2000.0,
            0.5,
            AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        assert_eq!(source.position_m(), Point3::new(2.0, 3.0, 4.0));
        assert_eq!(source.flux_lm(), 1000.0);
        close(
            source.intensity_cd(Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            1000.0 / PI,
            1e-12,
        );
        let mut integrated_cd = 0.0;
        for i in 0..SPECTRAL_BINS {
            integrated_cd += source
                .spectral_intensity_w_sr_nm(Vec3::new(1.0, 0.0, 0.0), wavelength_nm(i))
                .unwrap()
                * LUMINOUS_EFFICACY_LM_W
                * photopic_weight_nm(i);
        }
        close(integrated_cd, 1000.0 / PI, 1e-12);
        let scale = source.intensity_cd(Vec3::new(f64::MAX, 0.0, 0.0)).unwrap();
        close(scale, 1000.0 / PI, 1e-12);
    }

    #[test]
    fn passive_filter_reduces_flux_and_preserves_spectral_transmission() {
        let source = PhotometricSource::try_new(
            RigidTransform::identity(),
            1000.0,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let filter = std::array::from_fn(|i| if i < 40 { 0.1 } else { 0.8 });
        let filtered = source.filtered(filter).unwrap();
        assert!(filtered.flux_lm() < source.flux_lm());
        let direction = Vec3::new(0.0, 0.0, 1.0);
        for (i, transmission) in filter.iter().enumerate() {
            let reference = source
                .spectral_intensity_w_sr_nm(direction, wavelength_nm(i))
                .unwrap()
                * transmission;
            close(
                filtered
                    .spectral_intensity_w_sr_nm(direction, wavelength_nm(i))
                    .unwrap(),
                reference,
                1e-15,
            );
        }
        close(
            source.filtered([0.25; SPECTRAL_BINS]).unwrap().flux_lm(),
            250.0,
            1e-12,
        );
        assert_eq!(
            source.filtered([0.0; SPECTRAL_BINS]).unwrap().flux_lm(),
            0.0
        );
    }

    #[test]
    fn invalid_and_singular_inputs_refuse() {
        assert!(AngularDistribution::uniform_cone(0.0).is_err());
        assert!(AngularDistribution::cosine_power(-1.0, 0.3).is_err());
        assert!(
            AngularDistribution::tabulated(vec![0.0, PI], vec![0.0, TAU], vec![1.0, 2.0, 1.0, 1.0])
                .is_err()
        );
        assert!(
            AngularDistribution::tabulated(vec![0.0, PI], vec![0.0, TAU], vec![0.0; 4]).is_err()
        );
        assert!(SourceSpectrum::try_from_bins([0.0; SPECTRAL_BINS]).is_err());
        assert!(SourceSpectrum::try_from_bins([f64::NAN; SPECTRAL_BINS]).is_err());
        assert!(
            AngularDistribution::isotropic()
                .intensity_per_steradian(Vec3::new(0.0, 0.0, 0.0))
                .is_err()
        );
        assert!(
            PhotometricSource::try_new(
                RigidTransform::identity(),
                1.0,
                2.0,
                AngularDistribution::isotropic(),
                SourceSpectrum::equal_energy()
            )
            .is_err()
        );
    }
}
