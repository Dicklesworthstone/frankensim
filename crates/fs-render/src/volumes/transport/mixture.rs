//! Exact finite sums of the retained constituent optical models.

use super::{MAX_KNOTS, MassAbsorptionSpectrum, VolumeTransportError, checkpoint};
use fs_blake3::DomainHasher;
use fs_exec::Cx;

/// Hard bound on retained constituents and work at a mixture scattering vertex.
pub const MAX_OPTICAL_MIXTURE_COMPONENTS: usize = 64;

/// Caller-owned preparation and per-event constituent-work limit.
#[derive(Clone, Copy, Debug)]
pub struct OpticalMixtureBudget {
    /// Maximum supplied components, including declared zero-weight entries.
    /// The native hard bound is [`MAX_OPTICAL_MIXTURE_COMPONENTS`].
    pub max_components: usize,
}

impl MassAbsorptionSpectrum {
    /// Sum original mass coefficients with explicit nonnegative mass weights.
    ///
    /// A weight is kg of constituent per kg of the density field's declared
    /// reference mass. Fractions for a TOTAL mixture density sum to one; carrier
    /// loadings need not. This numerical operation never normalizes them or
    /// decides which density basis the caller has supplied.
    ///
    /// Every source must have identical wavelength endpoints. The exact union
    /// of active linear knots must fit the native 16-knot representation;
    /// otherwise admission refuses rather than removing spectral structure.
    /// No per-ray component sum or allocation is needed after preparation.
    pub fn weighted_sum(
        cx: &Cx<'_>,
        components: &[(f64, Self)],
        budget: OpticalMixtureBudget,
    ) -> Result<Self, VolumeTransportError> {
        validate_count(cx, components.len(), budget)?;
        let bounds = components[0].1.wavelength_bounds_nm();
        let mut wavelengths = [0.0; MAX_KNOTS];
        wavelengths[0] = bounds.0;
        wavelengths[1] = bounds.1;
        let mut count = 2;
        let mut hash = DomainHasher::new("org.frankensim.render.absorption-mixture.v1");
        hash.update(&(components.len() as u64).to_le_bytes());
        for &(weight, spectrum) in components {
            checkpoint(cx)?;
            validate_weight(weight)?;
            if spectrum.wavelength_bounds_nm() != bounds {
                return Err(VolumeTransportError::Wavelength);
            }
            hash.update(&weight.to_bits().to_le_bytes());
            hash_absorption(&spectrum, &mut hash);
            if weight > 0.0 {
                for &(wavelength, _) in &spectrum.knots[..spectrum.count] {
                    insert_knot(&mut wavelengths, &mut count, wavelength)?;
                }
            }
        }
        if components.len() == 1 && components[0].0 == 1.0 {
            return Ok(components[0].1);
        }
        let mut knots = [(0.0, 0.0); MAX_KNOTS];
        for (index, &wavelength) in wavelengths[..count].iter().enumerate() {
            checkpoint(cx)?;
            let mut total = PositiveSum::default();
            for &(weight, spectrum) in components {
                if weight > 0.0 {
                    total.add(weighted_product(weight, spectrum.at(wavelength)?)?)?;
                }
            }
            knots[index] = (wavelength, total.finish()?);
        }
        Self::source_table(&knots[..count], hash.finalize())
    }
}

pub(super) fn validate_count(
    cx: &Cx<'_>,
    count: usize,
    budget: OpticalMixtureBudget,
) -> Result<(), VolumeTransportError> {
    checkpoint(cx)?;
    if count == 0 {
        return Err(VolumeTransportError::Mixture);
    }
    if count > budget.max_components || count > MAX_OPTICAL_MIXTURE_COMPONENTS {
        return Err(VolumeTransportError::Budget);
    }
    Ok(())
}

pub(super) fn validate_weight(weight: f64) -> Result<(), VolumeTransportError> {
    if !weight.is_finite() || weight < 0.0 {
        return Err(VolumeTransportError::Mixture);
    }
    Ok(())
}

pub(super) fn weighted_product(weight: f64, coefficient: f64) -> Result<f64, VolumeTransportError> {
    let result = weight * coefficient;
    if !result.is_finite() || (weight > 0.0 && coefficient > 0.0 && result == 0.0) {
        return Err(VolumeTransportError::Unrepresentable);
    }
    Ok(result)
}

pub(super) fn hash_absorption(spectrum: &MassAbsorptionSpectrum, hash: &mut DomainHasher) {
    hash.update(&(spectrum.count as u64).to_le_bytes());
    for &(wavelength, coefficient) in &spectrum.knots[..spectrum.count] {
        hash.update(&wavelength.to_bits().to_le_bytes());
        hash.update(&coefficient.to_bits().to_le_bytes());
    }
    if let Some(source) = spectrum.source_identity {
        hash.update(b"source-table");
        hash.update(source.as_bytes());
    }
}

fn insert_knot(
    wavelengths: &mut [f64; MAX_KNOTS],
    count: &mut usize,
    wavelength: f64,
) -> Result<(), VolumeTransportError> {
    let index = wavelengths[..*count].partition_point(|&value| value < wavelength);
    if index < *count && wavelengths[index] == wavelength {
        return Ok(());
    }
    if *count == MAX_KNOTS {
        return Err(VolumeTransportError::Budget);
    }
    wavelengths.copy_within(index..*count, index + 1);
    wavelengths[index] = wavelength;
    *count += 1;
    Ok(())
}

/// Fixed-order Neumaier sum for nonnegative physical contributions.
#[derive(Default)]
pub(super) struct PositiveSum {
    sum: f64,
    correction: f64,
}

impl PositiveSum {
    pub(super) fn add(&mut self, value: f64) -> Result<(), VolumeTransportError> {
        let next = self.sum + value;
        if !value.is_finite() || value < 0.0 || !next.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        self.correction += if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        Ok(())
    }

    pub(super) fn finish(self) -> Result<f64, VolumeTransportError> {
        let result = self.sum + self.correction;
        if result.is_finite() && result >= 0.0 {
            Ok(result)
        } else {
            Err(VolumeTransportError::Unrepresentable)
        }
    }
}
