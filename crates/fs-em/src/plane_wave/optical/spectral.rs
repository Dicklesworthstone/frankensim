//! Native finite-band material optics in the original Maxwell response owner.

use super::{ReflectionError, VACUUM_SPEED_M_S, dispersion, positive};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::{
    field_storage::{FieldStorageError, linear_magnetic_permeability},
    state_point::{
        ResolvedMaterialStatePoint, SpectralComplexIndexStatePoint, SpectralConductorStatePoint,
        SpectralOpticalStatePoint, SpectralTabulatedDielectricStatePoint,
    },
};
use fs_math::{c64::C64, det};

const FREQUENCY_WAVELENGTH_M: f64 = core::f64::consts::TAU * VACUUM_SPEED_M_S;

/// Borrowed original spectral source, including its wavelength interpretation.
/// A source-selected SI conductor grid refreshes its wavelength properties;
/// the generic complex-index projection retains its explicitly fixed bindings.
/// Both use the SAME canonical complex-index interpolation owner.
#[derive(Clone, Copy, Debug)]
pub enum SpectralOpticalSource<'a> {
    ComplexIndex(&'a SpectralComplexIndexStatePoint),
    Conductor(&'a SpectralConductorStatePoint),
    TabulatedDielectric(&'a SpectralTabulatedDielectricStatePoint),
}

impl<'a> From<&'a SpectralOpticalStatePoint> for SpectralOpticalSource<'a> {
    fn from(source: &'a SpectralOpticalStatePoint) -> Self {
        match source {
            SpectralOpticalStatePoint::Conductor(value) => Self::Conductor(value),
            SpectralOpticalStatePoint::TabulatedDielectric(value) => {
                Self::TabulatedDielectric(value)
            }
        }
    }
}

impl<'a> SpectralOpticalSource<'a> {
    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        match self {
            Self::ComplexIndex(value) => value.resolved(),
            Self::Conductor(value) => value.resolved(),
            Self::TabulatedDielectric(value) => value.resolved(),
        }
    }

    fn complex_index(self) -> Option<&'a SpectralComplexIndexStatePoint> {
        match self {
            Self::ComplexIndex(value) => Some(value),
            Self::Conductor(value) => Some(value.complex_index()),
            Self::TabulatedDielectric(_) => None,
        }
    }
}

/// Piecewise source dispersion and absorption, on its original finite band.
///
/// Scalar absolute permeability MUST belong to the same complete snapshot.
/// The positive maximum knot span is authored in vacuum metres; it admits
/// interpolation, not unresolved resonances or an error certificate. There is
/// no nonmagnetic default, visible clamp, extrapolation or additional Ohmic loss.
/// These monochromatic data do not establish a causal time-domain reservoir
/// or a two-sided group derivative at a piecewise-linear source knot.
#[derive(Clone, Copy, Debug)]
pub struct SpectralOpticalResponse<'a> {
    source: SpectralOpticalSource<'a>,
    maximum_knot_span_m: f64,
    permeability_h_m: f64,
    index_to_sqrt_permittivity: f64,
}

impl<'a> SpectralOpticalResponse<'a> {
    pub fn try_new(
        source: SpectralOpticalSource<'a>,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        if !maximum_knot_span_m.is_finite() || maximum_knot_span_m <= 0.0 {
            return Err(ReflectionError::Invalid(
                "spectral optical interpolation needs a positive finite maximum knot span in metres",
            ));
        }
        let admit_gap = |gap: f64| {
            if gap.is_finite() && gap > 0.0 && gap <= maximum_knot_span_m {
                Ok(())
            } else {
                Err(ReflectionError::Invalid(
                    "spectral optical source exceeds its interpolation-gap or SI representability limit",
                ))
            }
        };
        if let Some(index) = source.complex_index() {
            for pair in index.samples().windows(2) {
                admit_gap((pair[1].wavelength_nm - pair[0].wavelength_nm) * 1e-9)?;
            }
        } else if let SpectralOpticalSource::TabulatedDielectric(index) = source {
            for pair in index.samples().windows(2) {
                admit_gap(pair[1].wavelength_m - pair[0].wavelength_m)?;
            }
        }
        let permeability_h_m = linear_magnetic_permeability(source.material())?;
        let index_to_sqrt_permittivity =
            positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        Ok(Self {
            source,
            maximum_knot_span_m,
            permeability_h_m,
            index_to_sqrt_permittivity,
        })
    }

    pub const fn source(self) -> SpectralOpticalSource<'a> {
        self.source
    }
    pub const fn material(self) -> &'a ResolvedMaterialStatePoint {
        self.source.material()
    }
    pub const fn maximum_knot_span_m(self) -> f64 {
        self.maximum_knot_span_m
    }
    pub const fn permeability_h_m(self) -> f64 {
        self.permeability_h_m
    }

    /// Original n-i*k at an admitted angular frequency [rad/s]. Evaluation
    /// allocates nothing. Dielectric n and intensity alpha are interpolated
    /// before k=alpha*lambda_vac/(4*pi), retaining the original loss spectrum.
    pub fn complex_index(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        if let Some(index) = self.source.complex_index() {
            let wavelength = dispersion::wavelength_from_frequency(
                omega_rad_s, super::FREQUENCY_WAVELENGTH_NM, index.wavelength_bounds_nm(),
            ).ok_or(ReflectionError::Invalid("wave frequency is outside the original complex-index source band or is unrepresentable"))?;
            let value = index
                .at_wavelength(wavelength)
                .map_err(FieldStorageError::from)?;
            return Ok(C64::new(value.eta, -value.k));
        }
        let SpectralOpticalSource::TabulatedDielectric(source) = self.source else {
            unreachable!("complex-index sources returned above")
        };
        let (low, high) = source.wavelength_bounds_m();
        let wavelength = dispersion::wavelength_from_frequency(
            omega_rad_s,
            FREQUENCY_WAVELENGTH_M,
            [low, high],
        )
        .ok_or(ReflectionError::Invalid(
            "wave frequency is outside the original dielectric source band or is unrepresentable",
        ))?;
        let samples = source.samples();
        let upper = samples.partition_point(|sample| sample.wavelength_m < wavelength);
        let (n, alpha) = if samples[upper].wavelength_m == wavelength {
            (samples[upper].phase_index, samples[upper].absorption_per_m)
        } else {
            let (left, right) = (samples[upper - 1], samples[upper]);
            let t = (wavelength - left.wavelength_m) / (right.wavelength_m - left.wavelength_m);
            let mix = |a: f64, b: f64| {
                if t <= 0.5 {
                    a + t * (b - a)
                } else {
                    b + (1.0 - t) * (a - b)
                }
            };
            (
                mix(left.phase_index, right.phase_index),
                mix(left.absorption_per_m, right.absorption_per_m),
            )
        };
        dispersion::attenuating_index_m(n, alpha, wavelength).ok_or(ReflectionError::Invalid(
            "spectral dielectric phase index or original intensity loss is unrepresentable",
        ))
    }

    pub fn permittivity_f_m(self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        dispersion::permittivity(
            self.complex_index(omega_rad_s)?,
            self.index_to_sqrt_permittivity,
        )
        .ok_or(ReflectionError::Invalid(
            "spectral optical permittivity or original loss is unrepresentable",
        ))
    }

    /// Refresh the exact original source claims, including same-state mu and
    /// source-selected SI wavelength coordinates. Re-admit the original gap
    /// policy before publishing an owned result; the input remains unchanged.
    /// The returned owner exposes a borrowed response without self-reference.
    pub fn at_query_point(
        self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<RefreshedSpectralOpticalSource, ReflectionError> {
        let source = match self.source {
            SpectralOpticalSource::ComplexIndex(value) => OwnedSource::ComplexIndex(
                value
                    .at_query_point(card, point)
                    .map_err(FieldStorageError::from)?,
            ),
            SpectralOpticalSource::Conductor(value) => OwnedSource::Conductor(
                value
                    .at_query_point(card, point)
                    .map_err(FieldStorageError::from)?,
            ),
            SpectralOpticalSource::TabulatedDielectric(value) => OwnedSource::TabulatedDielectric(
                value
                    .at_query_point(card, point)
                    .map_err(FieldStorageError::from)?,
            ),
        };
        let refreshed = RefreshedSpectralOpticalSource {
            source,
            maximum_knot_span_m: self.maximum_knot_span_m,
        };
        refreshed.response()?;
        Ok(refreshed)
    }

    pub(crate) fn same_source(self, other: Self) -> bool {
        self.maximum_knot_span_m == other.maximum_knot_span_m
            && match (self.source, other.source) {
                (
                    SpectralOpticalSource::ComplexIndex(a),
                    SpectralOpticalSource::ComplexIndex(b),
                ) => core::ptr::eq(a, b),
                (SpectralOpticalSource::Conductor(a), SpectralOpticalSource::Conductor(b)) => {
                    core::ptr::eq(a, b)
                }
                (
                    SpectralOpticalSource::TabulatedDielectric(a),
                    SpectralOpticalSource::TabulatedDielectric(b),
                ) => core::ptr::eq(a, b),
                _ => false,
            }
    }
}

#[derive(Clone, Debug)]
enum OwnedSource {
    ComplexIndex(SpectralComplexIndexStatePoint),
    Conductor(SpectralConductorStatePoint),
    TabulatedDielectric(SpectralTabulatedDielectricStatePoint),
}

/// A fully admitted source refresh with its original interpolation policy.
/// Source ownership is separate from any response borrowing it. This is an
/// endpoint constitutive observation, not a thermal evolution or path proof.
#[derive(Clone, Debug)]
pub struct RefreshedSpectralOpticalSource {
    source: OwnedSource,
    maximum_knot_span_m: f64,
}

impl RefreshedSpectralOpticalSource {
    pub fn source(&self) -> SpectralOpticalSource<'_> {
        match &self.source {
            OwnedSource::ComplexIndex(value) => SpectralOpticalSource::ComplexIndex(value),
            OwnedSource::Conductor(value) => SpectralOpticalSource::Conductor(value),
            OwnedSource::TabulatedDielectric(value) => {
                SpectralOpticalSource::TabulatedDielectric(value)
            }
        }
    }

    pub fn response(&self) -> Result<SpectralOpticalResponse<'_>, ReflectionError> {
        SpectralOpticalResponse::try_new(self.source(), self.maximum_knot_span_m)
    }
}
