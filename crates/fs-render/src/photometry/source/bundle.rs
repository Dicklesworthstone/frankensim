//! Point-source irradiance through a smooth, continuously refracting medium.
//!
//! The native fixture supplies spectral radiant intensity per launch solid
//! angle. The differentiated ray supplies receiving area per that same solid
//! angle. Their ratio, with the central ray's integrated absorption, gives
//! spectral irradiance. The projected area already includes receiving-plane
//! obliquity; no additional cosine or refractive-index factor is applied.
//!
//! Each supplied bundle must actually end at the declared receiving point,
//! within an explicit numerical position tolerance, and face the receiver.
//! This owner does not search connecting rays, test visibility, differentiate
//! discrete interfaces, or enumerate multiple optical paths. An 80-bin reading
//! describes one supplied connecting branch per wavelength, not a completeness
//! certificate for all paths from an extended or multiply imaged source.

use super::{
    BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, PhotometricSource, SPECTRAL_BINS, SourceError,
    photopic_weight_nm, unit, wavelength_nm,
};
use crate::gradient_index::bundle::{BundleError, ProjectionOptions, RayBundlePropagation};
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Admission failures leave an irradiance accumulator unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BundlePhotometryError {
    InvalidReceiver,
    SourceOriginMismatch,
    ReceiverPositionMismatch,
    BackFacingReceiver,
    InvalidBin,
    WavelengthMismatch,
    DuplicateBin,
    IncompleteSpectrum,
    NumericalRange,
    Source(SourceError),
    Bundle(BundleError),
}

impl From<SourceError> for BundlePhotometryError {
    fn from(error: SourceError) -> Self {
        Self::Source(error)
    }
}
impl From<BundleError> for BundlePhotometryError {
    fn from(error: BundleError) -> Self {
        Self::Bundle(error)
    }
}
impl core::fmt::Display for BundlePhotometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved bundle photometry: {self:?}")
    }
}
impl core::error::Error for BundlePhotometryError {}

/// A one-sided local irradiance receiver, with its normal toward arriving light.
/// Its tolerance admits numerical endpoint localization, not a finite-area
/// aperture or a bound on optical-model or irradiance error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarIrradianceReceiver {
    position_m: Point3,
    normal: Vec3,
    position_tolerance_m: f64,
}

impl PlanarIrradianceReceiver {
    pub fn try_new(
        position_m: Point3,
        normal: Vec3,
        position_tolerance_m: f64,
    ) -> Result<Self, BundlePhotometryError> {
        if ![position_m.x, position_m.y, position_m.z]
            .iter()
            .all(|v| v.is_finite())
            || !position_tolerance_m.is_finite()
            || position_tolerance_m < 0.0
        {
            return Err(BundlePhotometryError::InvalidReceiver);
        }
        let normal = unit(normal).map_err(|_| BundlePhotometryError::InvalidReceiver)?;
        Ok(Self {
            position_m,
            normal,
            position_tolerance_m,
        })
    }

    pub const fn position_m(self) -> Point3 {
        self.position_m
    }
    pub const fn normal(self) -> Vec3 {
        self.normal
    }
    pub const fn position_tolerance_m(self) -> f64 {
        self.position_tolerance_m
    }
}

impl PhotometricSource {
    /// Spectral irradiance in W/(m² nm) at an explicit receiver. The launch
    /// origin must match this source exactly. Native angular calibration is
    /// evaluated along the retained launch direction and at the traced vacuum
    /// wavelength; absorption is applied once by the bundle projection.
    pub fn spectral_irradiance_from_bundle(
        &self,
        bundle: &RayBundlePropagation,
        receiver: PlanarIrradianceReceiver,
        projection: ProjectionOptions,
    ) -> Result<f64, BundlePhotometryError> {
        if bundle.initial_ray().origin != self.position_m() {
            return Err(BundlePhotometryError::SourceOriginMismatch);
        }
        let endpoint = bundle.propagation().ray;
        let delta = endpoint.origin.delta_from(receiver.position_m);
        let offset = det::hypot(det::hypot(delta.x, delta.y), delta.z);
        if !offset.is_finite() || offset > receiver.position_tolerance_m {
            return Err(BundlePhotometryError::ReceiverPositionMismatch);
        }
        if receiver.normal.dot(endpoint.dir) >= 0.0 {
            return Err(BundlePhotometryError::BackFacingReceiver);
        }
        let intensity =
            self.spectral_intensity_w_sr_nm(bundle.initial_ray().dir, bundle.wavelength_nm())?;
        Ok(bundle
            .project_at_endpoint(receiver.normal, projection)?
            .point_source_irradiance(intensity)?)
    }
}

/// Fixed-size spectral quadrature over one supplied connecting bundle in each
/// native wavelength bin. The borrowed fixture and receiver remain unchanged.
/// Repeated bins refuse instead of silently duplicating one optical path.
pub struct BundleIrradianceMeter<'a> {
    source: &'a PhotometricSource,
    receiver: PlanarIrradianceReceiver,
    projection: ProjectionOptions,
    bins_w_m2_nm: [f64; SPECTRAL_BINS],
    present: u128,
}

impl<'a> BundleIrradianceMeter<'a> {
    pub fn new(
        source: &'a PhotometricSource,
        receiver: PlanarIrradianceReceiver,
        projection: ProjectionOptions,
    ) -> Self {
        Self {
            source,
            receiver,
            projection,
            bins_w_m2_nm: [0.0; SPECTRAL_BINS],
            present: 0,
        }
    }

    /// Record one completed central ray and its differentiated angular bundle.
    /// Every accepted sample uses the native midpoint wavelength exactly.
    pub fn record(
        &mut self,
        bin: usize,
        bundle: &RayBundlePropagation,
    ) -> Result<(), BundlePhotometryError> {
        if bin >= SPECTRAL_BINS {
            return Err(BundlePhotometryError::InvalidBin);
        }
        if bundle.wavelength_nm() != wavelength_nm(bin) {
            return Err(BundlePhotometryError::WavelengthMismatch);
        }
        let bit = 1_u128 << bin;
        if self.present & bit != 0 {
            return Err(BundlePhotometryError::DuplicateBin);
        }
        let value =
            self.source
                .spectral_irradiance_from_bundle(bundle, self.receiver, self.projection)?;
        if !value.is_finite() || value < 0.0 {
            return Err(BundlePhotometryError::NumericalRange);
        }
        self.bins_w_m2_nm[bin] = value;
        self.present |= bit;
        Ok(())
    }

    pub fn spectral_irradiance_w_m2_nm(&self, bin: usize) -> Option<f64> {
        (bin < SPECTRAL_BINS && self.present & (1_u128 << bin) != 0).then(|| self.bins_w_m2_nm[bin])
    }

    /// Integrate all 80 admitted bins with the original CIE observer. Missing
    /// wavelengths refuse; successful quadrature adds no visibility or global
    /// path-completeness claim to the supplied smooth-medium bundles.
    pub fn reading(&self) -> Result<BundleIrradianceReading, BundlePhotometryError> {
        if self.present != (1_u128 << SPECTRAL_BINS) - 1 {
            return Err(BundlePhotometryError::IncompleteSpectrum);
        }
        let mut radiant = Sum::default();
        let mut luminous = Sum::default();
        for (bin, value) in self.bins_w_m2_nm.iter().enumerate() {
            radiant.add(value * BIN_WIDTH_NM);
            luminous.add(value * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W);
        }
        let radiant_irradiance_w_m2 = radiant.value();
        let illuminance_lux = luminous.value();
        if !radiant_irradiance_w_m2.is_finite() || !illuminance_lux.is_finite() {
            return Err(BundlePhotometryError::NumericalRange);
        }
        Ok(BundleIrradianceReading {
            radiant_irradiance_w_m2,
            illuminance_lux,
        })
    }
}

/// Visible-band irradiance and photopic illuminance for the supplied branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleIrradianceReading {
    pub radiant_irradiance_w_m2: f64,
    pub illuminance_lux: f64,
}

#[derive(Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.total + value;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
    }
    fn value(self) -> f64 {
        self.total + self.correction
    }
}

#[cfg(test)]
mod tests;
