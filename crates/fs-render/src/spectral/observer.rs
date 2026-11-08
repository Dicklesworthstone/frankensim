//! Allocation-free visible projection shared by spectral volume observers.
//!
//! Transport owners supply their own units and quadrature. Keeping emission
//! and transmitted illumination separate avoids overflowing a spectral sum
//! whose weighted XYZ contribution is still representable.

use super::{cie_x, cie_y, cie_z};
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XyzObserverError {
    Wavelength,
    Weight,
    Radiance,
    Unrepresentable,
}

/// One prepared observer node. The original factors survive range loss in
/// their cached product, so even very small bin widths remain usable.
#[derive(Clone, Copy)]
pub(crate) struct XyzObserverLane {
    cmf: [f64; 3],
    weight: f64,
    weighted: [f64; 3],
}

impl XyzObserverLane {
    pub(crate) const ZERO: Self = Self {
        cmf: [0.0; 3],
        weight: 0.0,
        weighted: [0.0; 3],
    };

    pub(crate) fn new(
        wavelength_nm: f64,
        spectral_weight: f64,
        normalization: f64,
    ) -> Result<Self, XyzObserverError> {
        if !wavelength_nm.is_finite() || !(380.0..=780.0).contains(&wavelength_nm) {
            return Err(XyzObserverError::Wavelength);
        }
        if !spectral_weight.is_finite()
            || spectral_weight < 0.0
            || !normalization.is_finite()
            || normalization <= 0.0
        {
            return Err(XyzObserverError::Weight);
        }
        let cmf = [
            cie_x(wavelength_nm),
            cie_y(wavelength_nm),
            cie_z(wavelength_nm),
        ]
        .map(|value| value * normalization);
        Ok(Self {
            cmf,
            weight: spectral_weight,
            weighted: cmf.map(|value| value * spectral_weight),
        })
    }

    fn project(self, value: f64, channel: usize, spatial_weight: f64) -> f64 {
        if value == 0.0 || self.weight == 0.0 || self.cmf[channel] == 0.0 || spatial_weight == 0.0 {
            return 0.0;
        }
        let weighted = self.weighted[channel];
        let product = value * weighted * spatial_weight;
        if weighted.is_normal() && product.is_normal() {
            product
        } else {
            det::exp(
                det::ln(value)
                    + det::ln(self.cmf[channel])
                    + det::ln(self.weight)
                    + det::ln(spatial_weight),
            )
        }
    }

    fn transmit(
        self,
        background: f64,
        transmission: f64,
        depth: f64,
        channel: usize,
        spatial_weight: f64,
    ) -> f64 {
        if background == 0.0
            || self.weight == 0.0
            || self.cmf[channel] == 0.0
            || spatial_weight == 0.0
        {
            return 0.0;
        }
        let weighted = self.weighted[channel];
        let attenuated = background * transmission;
        let product = attenuated * weighted * spatial_weight;
        if transmission.is_normal()
            && attenuated.is_normal()
            && weighted.is_normal()
            && product.is_normal()
        {
            product
        } else {
            // A stored zero/subnormal transmission can have finite log survival.
            // Use it before applying bright illumination or a large bin weight.
            let log_transmission = if transmission < f64::MIN_POSITIVE {
                -depth
            } else {
                det::ln(transmission)
            };
            det::exp(
                det::ln(background)
                    + log_transmission
                    + det::ln(self.cmf[channel])
                    + det::ln(self.weight)
                    + det::ln(spatial_weight),
            )
        }
    }
}

/// Nonnegative compensated XYZ accumulation. An overflowing channel cannot be
/// cancelled by a later physical contribution and is refused at its source.
#[derive(Default)]
pub(crate) struct XyzAccumulator {
    value: [f64; 3],
    correction: [f64; 3],
}

impl XyzAccumulator {
    pub(crate) fn add(
        &mut self,
        observer: &XyzObserverLane,
        emission: f64,
        transmission: f64,
        depth: f64,
        background: f64,
    ) -> Result<(), XyzObserverError> {
        self.add_scaled(observer, emission, transmission, depth, background, 1.0)
    }

    /// Include a nonnegative spatial quadrature weight before rounding any
    /// spectral contribution. This preserves finite pixel integrals when an
    /// unweighted column overflows or an intermediate radiance underflows.
    pub(crate) fn add_scaled(
        &mut self,
        observer: &XyzObserverLane,
        emission: f64,
        transmission: f64,
        depth: f64,
        background: f64,
        spatial_weight: f64,
    ) -> Result<(), XyzObserverError> {
        if !spatial_weight.is_finite() || spatial_weight < 0.0 {
            return Err(XyzObserverError::Weight);
        }
        if !emission.is_finite()
            || emission < 0.0
            || !background.is_finite()
            || background < 0.0
            || !transmission.is_finite()
            || !(0.0..=1.0).contains(&transmission)
            || depth.is_nan()
            || depth < 0.0
        {
            return Err(XyzObserverError::Radiance);
        }
        for channel in 0..3 {
            for term in [
                observer.project(emission, channel, spatial_weight),
                observer.transmit(background, transmission, depth, channel, spatial_weight),
            ] {
                if term == 0.0 {
                    continue;
                }
                let corrected = term - self.correction[channel];
                let next = self.value[channel] + corrected;
                if !term.is_finite() || !next.is_finite() {
                    return Err(XyzObserverError::Unrepresentable);
                }
                self.correction[channel] = (next - self.value[channel]) - corrected;
                self.value[channel] = next;
            }
        }
        Ok(())
    }

    pub(crate) const fn finish(self) -> Result<[f64; 3], XyzObserverError> {
        Ok(self.value)
    }
}

#[cfg(test)]
mod tests;
