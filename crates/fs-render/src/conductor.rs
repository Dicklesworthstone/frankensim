//! Validated spectral conductor optics for the frontier path tracer.
//!
//! Wavelengths are vacuum nanometres.  Each table stores the absolute real
//! refractive index `eta` and extinction coefficient `k`; Fresnel evaluation
//! divides both by the incident dielectric index before solving the exact
//! complex interface equations.  Built-in metal tables are representative
//! look-development assets, not measurements of a specimen or finish.

use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;

use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};

/// Bit semantics of table validation, interpolation, and Fresnel evaluation.
pub const CONDUCTOR_OPTICS_SEMANTICS_VERSION: u32 = 1;
/// Bit semantics of the single-scattering isotropic-GGX conductor BSDF.
pub const CONDUCTOR_BSDF_SEMANTICS_VERSION: u32 = 1;
/// Fixed number of knots in a conductor optical table.
pub const CONDUCTOR_IOR_SAMPLE_COUNT: usize = 9;
/// Largest accepted interval between adjacent optical knots.
pub const MAX_CONDUCTOR_KNOT_SPAN_NM: f64 = 100.0;

const MAX_IOR_COMPONENT: f64 = 100.0;
const MIN_ROUGHNESS_ALPHA: f64 = 1.0e-4;
const MAX_ROUGHNESS_ALPHA: f64 = 1.0;
const OPTICS_IDENTITY_DOMAIN: &str = "org.frankensim.fs-render.conductor-optics.v1";
const SOURCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-render.conductor-source.v1";

/// Refusal from conductor table construction or analytic evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConductorError {
    /// One knot contained a non-finite or out-of-domain value.
    InvalidSample,
    /// Knots were not canonical, did not cover the visible interval, or left
    /// a gap larger than the declared interpolation envelope.
    InvalidTable,
    /// A wavelength was outside the table's inclusive support.
    InvalidWavelength,
    /// The incident dielectric index was non-finite or non-positive.
    InvalidIncidentMedium,
    /// An incidence cosine was non-finite or outside `[0, 1]`.
    InvalidIncidence,
    /// Isotropic GGX roughness was outside the admitted interval.
    InvalidRoughness,
    /// A source identity was zero.
    InvalidSource,
    /// Valid inputs nevertheless produced an invalid complex intermediate.
    NumericalFailure,
}

impl core::fmt::Display for ConductorError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str(match self {
            Self::InvalidSample => "invalid conductor optical sample",
            Self::InvalidTable => "invalid conductor optical table",
            Self::InvalidWavelength => "wavelength is outside the conductor table",
            Self::InvalidIncidentMedium => "invalid incident-medium refractive index",
            Self::InvalidIncidence => "invalid conductor incidence cosine",
            Self::InvalidRoughness => "invalid conductor GGX roughness",
            Self::InvalidSource => "invalid conductor source identity",
            Self::NumericalFailure => "conductor Fresnel evaluation was non-finite",
        })
    }
}

impl core::error::Error for ConductorError {}

/// Whether a table is an artistic representative or a caller-retained claim
/// of measurement.  The latter is metadata only; this module does not verify
/// calibration, licensing, specimen identity, or measurement uncertainty.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ConductorDataStatus {
    /// Representative visual data with no specimen claim.
    Representative,
    /// Caller assertion that the source contains measured optical constants.
    CallerAssertedMeasured,
}

/// Caller-supplied source identifier and declared status of admitted optical
/// data. Built-in identifiers bind their complete canonical tables; custom
/// identifiers remain caller assertions. In all cases,
/// [`ConductorOptics::content_identity`] binds the actual admitted knots.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ConductorSource {
    identity: ContentHash,
    status: ConductorDataStatus,
}

impl ConductorSource {
    /// Attach a nonzero caller-supplied source identifier and honest status.
    pub fn try_new(
        identity: ContentHash,
        status: ConductorDataStatus,
    ) -> Result<Self, ConductorError> {
        if identity.as_bytes().iter().all(|byte| *byte == 0) {
            return Err(ConductorError::InvalidSource);
        }
        Ok(Self { identity, status })
    }

    /// Source identifier. Built-in values bind their canonical table bytes;
    /// custom values are retained caller assertions and confer no authority.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }

    /// Declared representative-versus-measured status.
    #[must_use]
    pub const fn status(self) -> ConductorDataStatus {
        self.status
    }
}

/// Provenance class attached to a validated conductor table.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ConductorProvenance {
    /// Caller-supplied knots and source binding.
    Custom,
    /// Built-in representative tungsten-like visual table v1.
    RepresentativeTungstenV1,
    /// Built-in representative stainless-steel-like visual table v1.
    RepresentativeStainlessSteelV1,
}

/// One absolute complex-index knot at a vacuum wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductorIorSample {
    wavelength_nm: f64,
    eta: f64,
    k: f64,
}

impl ConductorIorSample {
    /// Construct one visible-band optical knot.
    pub fn try_new(wavelength_nm: f64, eta: f64, k: f64) -> Result<Self, ConductorError> {
        let sample = Self {
            wavelength_nm: canonical_zero(wavelength_nm),
            eta: canonical_zero(eta),
            k: canonical_zero(k),
        };
        if !sample_is_valid(sample) {
            return Err(ConductorError::InvalidSample);
        }
        Ok(sample)
    }

    /// Vacuum wavelength in nanometres.
    #[must_use]
    pub const fn wavelength_nm(self) -> f64 {
        self.wavelength_nm
    }

    /// Absolute real refractive index.
    #[must_use]
    pub const fn eta(self) -> f64 {
        self.eta
    }

    /// Absolute extinction coefficient.
    #[must_use]
    pub const fn k(self) -> f64 {
        self.k
    }
}

/// Interpolated absolute complex refractive index.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexIor {
    eta: f64,
    k: f64,
}

impl ComplexIor {
    /// Absolute real refractive index.
    #[must_use]
    pub const fn eta(self) -> f64 {
        self.eta
    }

    /// Absolute extinction coefficient.
    #[must_use]
    pub const fn k(self) -> f64 {
        self.k
    }
}

/// Validated fixed-band complex-IOR table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductorOptics {
    samples: [ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT],
    source: ConductorSource,
    provenance: ConductorProvenance,
}

impl ConductorOptics {
    /// Admit a caller-supplied table.  Knots must be strictly increasing,
    /// cover exactly 380--780 nm, and leave no interval above the declared
    /// interpolation ceiling.
    pub fn try_new(
        samples: [ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT],
        source: ConductorSource,
    ) -> Result<Self, ConductorError> {
        Self::from_samples(samples, source, ConductorProvenance::Custom)
    }

    /// Representative tungsten-like optical table.  These values are an
    /// uncalibrated visual starting point, not measured stock.
    #[must_use]
    pub fn representative_tungsten() -> Self {
        Self::from_samples(
            TUNGSTEN_SAMPLES,
            representative_source(b"representative-tungsten-v1", &TUNGSTEN_SAMPLES),
            ConductorProvenance::RepresentativeTungstenV1,
        )
        .expect("built-in tungsten table is valid")
    }

    /// Representative stainless-steel-like optical table.  Alloy, passive
    /// film, machining, and finish are deliberately not inferred.
    #[must_use]
    pub fn representative_stainless_steel() -> Self {
        Self::from_samples(
            STAINLESS_STEEL_SAMPLES,
            representative_source(
                b"representative-stainless-steel-v1",
                &STAINLESS_STEEL_SAMPLES,
            ),
            ConductorProvenance::RepresentativeStainlessSteelV1,
        )
        .expect("built-in stainless-steel table is valid")
    }

    fn from_samples(
        samples: [ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT],
        source: ConductorSource,
        provenance: ConductorProvenance,
    ) -> Result<Self, ConductorError> {
        validate_table(&samples)?;
        Ok(Self {
            samples,
            source,
            provenance,
        })
    }

    /// Evaluate linearly interpolated absolute `eta + i k`.
    pub fn complex_ior(self, wavelength_nm: f64) -> Result<ComplexIor, ConductorError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(ConductorError::InvalidWavelength);
        }
        Ok(self.complex_ior_unchecked(wavelength_nm))
    }

    pub(crate) fn complex_ior_unchecked(self, wavelength_nm: f64) -> ComplexIor {
        if wavelength_nm.to_bits() == LAMBDA_MAX.to_bits() {
            let sample = self.samples[CONDUCTOR_IOR_SAMPLE_COUNT - 1];
            return ComplexIor {
                eta: sample.eta,
                k: sample.k,
            };
        }
        let upper = self
            .samples
            .partition_point(|sample| sample.wavelength_nm <= wavelength_nm)
            .clamp(1, CONDUCTOR_IOR_SAMPLE_COUNT - 1);
        let lower = self.samples[upper - 1];
        let upper = self.samples[upper];
        if wavelength_nm.to_bits() == lower.wavelength_nm.to_bits() {
            return ComplexIor {
                eta: lower.eta,
                k: lower.k,
            };
        }
        let fraction =
            (wavelength_nm - lower.wavelength_nm) / (upper.wavelength_nm - lower.wavelength_nm);
        ComplexIor {
            eta: lower.eta + fraction * (upper.eta - lower.eta),
            k: lower.k + fraction * (upper.k - lower.k),
        }
    }

    /// Exact unpolarized Fresnel reflectance at this table's wavelength.
    pub fn fresnel(
        self,
        wavelength_nm: f64,
        incident_eta: f64,
        incidence_cosine: f64,
    ) -> Result<f64, ConductorError> {
        let ior = self.complex_ior(wavelength_nm)?;
        fresnel_conductor(incident_eta, ior.eta, ior.k, incidence_cosine)
    }

    /// Canonical knots in ascending wavelength order.
    #[must_use]
    pub const fn samples(self) -> [ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT] {
        self.samples
    }

    /// Bound source identity and declared status.
    #[must_use]
    pub const fn source(self) -> ConductorSource {
        self.source
    }

    /// Built-in or caller-supplied provenance class.
    #[must_use]
    pub const fn provenance(self) -> ConductorProvenance {
        self.provenance
    }

    /// Identity of the exact knots, conventions, interpolation policy, and
    /// source binding.  Surface roughness is intentionally separate.
    #[must_use]
    pub fn content_identity(self) -> ContentHash {
        let mut hasher = DomainHasher::new(OPTICS_IDENTITY_DOMAIN);
        hasher.update(&CONDUCTOR_OPTICS_SEMANTICS_VERSION.to_le_bytes());
        hasher.update(&[match self.provenance {
            ConductorProvenance::Custom => 0,
            ConductorProvenance::RepresentativeTungstenV1 => 1,
            ConductorProvenance::RepresentativeStainlessSteelV1 => 2,
        }]);
        hasher.update(&[match self.source.status {
            ConductorDataStatus::Representative => 0,
            ConductorDataStatus::CallerAssertedMeasured => 1,
        }]);
        hasher.update(self.source.identity.as_bytes());
        hasher.update(&LAMBDA_MIN.to_bits().to_le_bytes());
        hasher.update(&LAMBDA_MAX.to_bits().to_le_bytes());
        hasher.update(&MAX_CONDUCTOR_KNOT_SPAN_NM.to_bits().to_le_bytes());
        for sample in self.samples {
            hasher.update(&sample.wavelength_nm.to_bits().to_le_bytes());
            hasher.update(&sample.eta.to_bits().to_le_bytes());
            hasher.update(&sample.k.to_bits().to_le_bytes());
        }
        hasher.finalize()
    }
}

/// Validated isotropic single-scattering GGX surface for conductor optics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductorSurface {
    roughness_alpha: f64,
}

impl ConductorSurface {
    /// Construct an isotropic GGX surface.  Roughness changes the microfacet
    /// distribution but never rewrites optical constants.
    pub fn try_rough(roughness_alpha: f64) -> Result<Self, ConductorError> {
        if !roughness_alpha.is_finite()
            || !(MIN_ROUGHNESS_ALPHA..=MAX_ROUGHNESS_ALPHA).contains(&roughness_alpha)
        {
            return Err(ConductorError::InvalidRoughness);
        }
        Ok(Self {
            roughness_alpha: canonical_zero(roughness_alpha),
        })
    }

    /// Isotropic GGX alpha.
    #[must_use]
    pub const fn roughness_alpha(self) -> f64 {
        self.roughness_alpha
    }
}

/// Exact unpolarized conductor Fresnel for absolute transmitted `eta`, `k`
/// and a real incident dielectric index.
pub fn fresnel_conductor(
    incident_eta: f64,
    transmitted_eta: f64,
    transmitted_k: f64,
    incidence_cosine: f64,
) -> Result<f64, ConductorError> {
    if !incident_eta.is_finite() || incident_eta <= 0.0 || incident_eta > MAX_IOR_COMPONENT {
        return Err(ConductorError::InvalidIncidentMedium);
    }
    if !transmitted_eta.is_finite()
        || !transmitted_k.is_finite()
        || transmitted_eta <= 0.0
        || transmitted_k <= 0.0
        || transmitted_eta > MAX_IOR_COMPONENT
        || transmitted_k > MAX_IOR_COMPONENT
    {
        return Err(ConductorError::InvalidSample);
    }
    if !incidence_cosine.is_finite() || !(0.0..=1.0).contains(&incidence_cosine) {
        return Err(ConductorError::InvalidIncidence);
    }
    if incidence_cosine == 0.0 {
        return Ok(1.0);
    }

    let cosine = incidence_cosine;
    let relative = Complex::new(
        transmitted_eta / incident_eta,
        -transmitted_k / incident_eta,
    );
    let one = Complex::new(1.0, 0.0);
    let sine_squared = 1.0 - cosine * cosine;
    let transmitted_cosine =
        complex_sqrt(one - Complex::new(sine_squared, 0.0) / relative.square());
    let parallel =
        (relative * cosine - transmitted_cosine) / (relative * cosine + transmitted_cosine);
    let perpendicular = (Complex::new(cosine, 0.0) - relative * transmitted_cosine)
        / (Complex::new(cosine, 0.0) + relative * transmitted_cosine);
    let reflectance = 0.5 * (parallel.norm_squared() + perpendicular.norm_squared());
    if !reflectance.is_finite() || !(-1.0e-12..=1.0 + 1.0e-12).contains(&reflectance) {
        return Err(ConductorError::NumericalFailure);
    }
    Ok(reflectance.clamp(0.0, 1.0))
}

#[derive(Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    fn square(self) -> Self {
        self * self
    }

    fn norm_squared(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

impl core::ops::Add for Complex {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl core::ops::Sub for Complex {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl core::ops::Mul for Complex {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl core::ops::Mul<f64> for Complex {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl core::ops::Div for Complex {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let denominator = rhs.norm_squared();
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denominator,
            (self.im * rhs.re - self.re * rhs.im) / denominator,
        )
    }
}

fn complex_sqrt(value: Complex) -> Complex {
    let magnitude = det::sqrt(value.norm_squared());
    let real = det::sqrt(((magnitude + value.re) * 0.5).max(0.0));
    let imaginary_magnitude = det::sqrt(((magnitude - value.re) * 0.5).max(0.0));
    let imaginary = if value.im < 0.0 {
        -imaginary_magnitude
    } else {
        imaginary_magnitude
    };
    Complex::new(real, imaginary)
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn sample_is_valid(sample: ConductorIorSample) -> bool {
    sample.wavelength_nm.is_finite()
        && (LAMBDA_MIN..=LAMBDA_MAX).contains(&sample.wavelength_nm)
        && sample.eta.is_finite()
        && sample.eta > 0.0
        && sample.eta <= MAX_IOR_COMPONENT
        && sample.k.is_finite()
        && sample.k > 0.0
        && sample.k <= MAX_IOR_COMPONENT
}

fn validate_table(
    samples: &[ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT],
) -> Result<(), ConductorError> {
    if samples
        .iter()
        .copied()
        .any(|sample| !sample_is_valid(sample))
    {
        return Err(ConductorError::InvalidSample);
    }
    if samples[0].wavelength_nm.to_bits() != LAMBDA_MIN.to_bits()
        || samples[CONDUCTOR_IOR_SAMPLE_COUNT - 1]
            .wavelength_nm
            .to_bits()
            != LAMBDA_MAX.to_bits()
    {
        return Err(ConductorError::InvalidTable);
    }
    for pair in samples.windows(2) {
        let span = pair[1].wavelength_nm - pair[0].wavelength_nm;
        if !span.is_finite() || span <= 0.0 || span > MAX_CONDUCTOR_KNOT_SPAN_NM {
            return Err(ConductorError::InvalidTable);
        }
    }
    Ok(())
}

fn representative_source(
    label: &[u8],
    samples: &[ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT],
) -> ConductorSource {
    let mut hasher = DomainHasher::new(SOURCE_IDENTITY_DOMAIN);
    hasher.update(label);
    for sample in samples {
        hasher.update(&sample.wavelength_nm.to_bits().to_le_bytes());
        hasher.update(&sample.eta.to_bits().to_le_bytes());
        hasher.update(&sample.k.to_bits().to_le_bytes());
    }
    ConductorSource {
        identity: hasher.finalize(),
        status: ConductorDataStatus::Representative,
    }
}

const fn sample(wavelength_nm: f64, eta: f64, k: f64) -> ConductorIorSample {
    ConductorIorSample {
        wavelength_nm,
        eta,
        k,
    }
}

// Representative visual curves.  They intentionally encode no grade, oxide,
// roughness, temperature, machining, or specimen claim.
const TUNGSTEN_SAMPLES: [ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT] = [
    sample(380.0, 3.45, 2.62),
    sample(430.0, 3.58, 2.74),
    sample(480.0, 3.72, 2.88),
    sample(530.0, 3.90, 3.04),
    sample(580.0, 4.12, 3.23),
    sample(630.0, 4.36, 3.45),
    sample(680.0, 4.58, 3.70),
    sample(730.0, 4.78, 3.96),
    sample(780.0, 4.96, 4.22),
];

const STAINLESS_STEEL_SAMPLES: [ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT] = [
    sample(380.0, 2.48, 3.54),
    sample(430.0, 2.55, 3.66),
    sample(480.0, 2.62, 3.78),
    sample(530.0, 2.69, 3.90),
    sample(580.0, 2.76, 4.02),
    sample(630.0, 2.83, 4.14),
    sample(680.0, 2.90, 4.26),
    sample(730.0, 2.97, 4.38),
    sample(780.0, 3.04, 4.50),
];

#[cfg(test)]
mod tests {
    use super::*;

    fn source(label: &[u8], status: ConductorDataStatus) -> ConductorSource {
        ConductorSource::try_new(
            fs_blake3::hash_domain("org.frankensim.test.conductor-source.v1", label),
            status,
        )
        .expect("test source")
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64, context: &str) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{context}: actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    fn real_form_oracle(incident_eta: f64, eta: f64, k: f64, cosine: f64) -> f64 {
        let eta = eta / incident_eta;
        let k = k / incident_eta;
        let cosine_squared = cosine * cosine;
        let sine_squared = 1.0 - cosine_squared;
        let eta_squared = eta * eta;
        let k_squared = k * k;
        let t0 = eta_squared - k_squared - sine_squared;
        let a2_plus_b2 = (t0 * t0 + 4.0 * eta_squared * k_squared).sqrt();
        let a = (0.5 * (a2_plus_b2 + t0)).sqrt();
        let t1 = a2_plus_b2 + cosine_squared;
        let t2 = 2.0 * cosine * a;
        let perpendicular = (t1 - t2) / (t1 + t2);
        let t3 = cosine_squared * a2_plus_b2 + sine_squared * sine_squared;
        let t4 = t2 * sine_squared;
        let parallel = perpendicular * (t3 - t4) / (t3 + t4);
        0.5 * (parallel + perpendicular)
    }

    #[test]
    fn exact_fresnel_matches_known_complex_ior_values() {
        for (incident, eta, k, cosine, expected) in [
            (1.0, 0.2, 3.0, 1.0, 0.923_371_647_509_578_6),
            (1.0, 0.2, 3.0, 0.5, 0.918_411_084_659_368_9),
            (
                1.0,
                0.2,
                3.0,
                0.173_648_177_666_930_41,
                0.938_558_369_711_149_9,
            ),
            (
                1.33,
                0.2,
                3.0,
                core::f64::consts::FRAC_1_SQRT_2,
                0.905_565_341_833_169,
            ),
            (
                1.0,
                3.0,
                4.0,
                core::f64::consts::FRAC_1_SQRT_2,
                0.617_383_883_096_731_8,
            ),
        ] {
            assert_close(
                fresnel_conductor(incident, eta, k, cosine).expect("valid Fresnel"),
                expected,
                2.0e-14,
                "complex Fresnel known answer",
            );
        }
    }

    #[test]
    fn complex_fresnel_matches_independent_real_form_grid() {
        for incident_eta in [1.0, 1.33, 1.52] {
            for (eta, k) in [(0.2, 3.0), (2.7, 4.1), (4.8, 2.6)] {
                for cosine_step in 0..=32 {
                    let cosine = f64::from(cosine_step) / 32.0;
                    let expected = real_form_oracle(incident_eta, eta, k, cosine);
                    let actual = fresnel_conductor(incident_eta, eta, k, cosine).unwrap();
                    assert_close(actual, expected, 5.0e-14, "independent real-form Fresnel");
                }
            }
        }
    }

    #[test]
    fn normal_grazing_and_common_scale_relations_hold() {
        for (incident, eta, k) in [(1.0, 0.2, 3.0), (1.33, 2.7, 4.1), (1.0, 4.0, 2.5)] {
            let expected = ((eta - incident) * (eta - incident) + k * k)
                / ((eta + incident) * (eta + incident) + k * k);
            assert_close(
                fresnel_conductor(incident, eta, k, 1.0).unwrap(),
                expected,
                2.0e-14,
                "normal incidence",
            );
            assert_eq!(
                fresnel_conductor(incident, eta, k, 0.0).unwrap().to_bits(),
                1.0_f64.to_bits()
            );
            let base = fresnel_conductor(incident, eta, k, 0.37).unwrap();
            let scaled = fresnel_conductor(incident * 2.0, eta * 2.0, k * 2.0, 0.37).unwrap();
            assert_close(scaled, base, 2.0e-14, "common IOR scale");
        }
    }

    #[test]
    fn tables_interpolate_exact_knots_and_refuse_outside_support() {
        let optics = ConductorOptics::representative_tungsten();
        let knots = optics.samples();
        for knot in knots {
            let evaluated = optics.complex_ior(knot.wavelength_nm()).unwrap();
            assert_eq!(evaluated.eta().to_bits(), knot.eta().to_bits());
            assert_eq!(evaluated.k().to_bits(), knot.k().to_bits());
        }
        let midpoint = optics.complex_ior(405.0).unwrap();
        assert_close(
            midpoint.eta(),
            (knots[0].eta() + knots[1].eta()) * 0.5,
            0.0,
            "eta midpoint",
        );
        assert_close(
            midpoint.k(),
            (knots[0].k() + knots[1].k()) * 0.5,
            0.0,
            "k midpoint",
        );
        assert_eq!(
            optics.complex_ior(379.999),
            Err(ConductorError::InvalidWavelength)
        );
        assert_eq!(
            optics.complex_ior(780.001),
            Err(ConductorError::InvalidWavelength)
        );
    }

    #[test]
    fn malformed_samples_tables_sources_and_surfaces_refuse() {
        assert_eq!(
            ConductorIorSample::try_new(500.0, f64::NAN, 3.0),
            Err(ConductorError::InvalidSample)
        );
        assert_eq!(
            ConductorIorSample::try_new(500.0, 2.0, 0.0),
            Err(ConductorError::InvalidSample)
        );
        assert_eq!(
            ConductorSource::try_new(
                ContentHash([0; 32]),
                ConductorDataStatus::CallerAssertedMeasured,
            ),
            Err(ConductorError::InvalidSource)
        );
        let mut duplicate = TUNGSTEN_SAMPLES;
        duplicate[2].wavelength_nm = duplicate[1].wavelength_nm;
        assert_eq!(
            ConductorOptics::try_new(
                duplicate,
                source(b"duplicate", ConductorDataStatus::Representative),
            ),
            Err(ConductorError::InvalidTable)
        );
        let mut incomplete = TUNGSTEN_SAMPLES;
        incomplete[0].wavelength_nm = 381.0;
        assert_eq!(
            ConductorOptics::try_new(
                incomplete,
                source(b"incomplete", ConductorDataStatus::Representative),
            ),
            Err(ConductorError::InvalidTable)
        );
        for alpha in [0.0, 1.0e-5, 1.000_001, f64::INFINITY] {
            assert_eq!(
                ConductorSurface::try_rough(alpha),
                Err(ConductorError::InvalidRoughness)
            );
        }
    }

    #[test]
    fn presets_and_custom_tables_have_sensitive_replay_identities() {
        let tungsten = ConductorOptics::representative_tungsten();
        let replay = ConductorOptics::representative_tungsten();
        let stainless = ConductorOptics::representative_stainless_steel();
        assert_eq!(tungsten.content_identity(), replay.content_identity());
        assert_ne!(tungsten.content_identity(), stainless.content_identity());
        assert_eq!(
            tungsten.source().status(),
            ConductorDataStatus::Representative
        );
        assert_eq!(
            stainless.source().status(),
            ConductorDataStatus::Representative
        );

        let custom_source = source(b"custom", ConductorDataStatus::CallerAssertedMeasured);
        let custom = ConductorOptics::try_new(TUNGSTEN_SAMPLES, custom_source).unwrap();
        assert_eq!(custom.provenance(), ConductorProvenance::Custom);
        assert_ne!(custom.content_identity(), tungsten.content_identity());
        let mut perturbed = TUNGSTEN_SAMPLES;
        perturbed[4].eta = f64::from_bits(perturbed[4].eta.to_bits() + 1);
        let changed = ConductorOptics::try_new(perturbed, custom_source).unwrap();
        assert_ne!(changed.content_identity(), custom.content_identity());
        assert_ne!(
            representative_source(b"source-byte-binding", &TUNGSTEN_SAMPLES).identity(),
            representative_source(b"source-byte-binding", &perturbed).identity(),
            "built-in source identity must bind the canonical optical table"
        );
    }

    #[test]
    fn spectral_packet_permutation_does_not_change_lane_values() {
        let optics = ConductorOptics::representative_stainless_steel();
        let wavelengths = [410.0, 510.0, 610.0, 710.0];
        let original = wavelengths.map(|wavelength| optics.fresnel(wavelength, 1.0, 0.63).unwrap());
        let permuted_wavelengths = [
            wavelengths[2],
            wavelengths[0],
            wavelengths[3],
            wavelengths[1],
        ];
        let permuted =
            permuted_wavelengths.map(|wavelength| optics.fresnel(wavelength, 1.0, 0.63).unwrap());
        assert_eq!(original[0].to_bits(), permuted[1].to_bits());
        assert_eq!(original[1].to_bits(), permuted[3].to_bits());
        assert_eq!(original[2].to_bits(), permuted[0].to_bits());
        assert_eq!(original[3].to_bits(), permuted[2].to_bits());
    }

    #[test]
    fn admitted_grid_is_finite_and_energy_bounded() {
        for optics in [
            ConductorOptics::representative_tungsten(),
            ConductorOptics::representative_stainless_steel(),
        ] {
            for wavelength_step in 0..=40 {
                let wavelength =
                    LAMBDA_MIN + (LAMBDA_MAX - LAMBDA_MIN) * f64::from(wavelength_step) / 40.0;
                for cosine_step in 0..=64 {
                    let value = optics
                        .fresnel(wavelength, 1.0, f64::from(cosine_step) / 64.0)
                        .unwrap();
                    assert!(value.is_finite() && (0.0..=1.0).contains(&value));
                }
            }
        }
    }
}
