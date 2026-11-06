use super::MAX_ALPHA_BINS;
use crate::fusion::{DtFusionError, poll, positive};
use fs_math::det::{exp, ln, sqrt};
use fs_qty::{Amount, Energy};
use std::sync::Arc;

/// Original positive kinetic-energy cell centers and their integration widths.
///
/// There are 2–256 strictly increasing centers. Interior faces bisect adjacent
/// centers; the lower face is zero and the upper face lies half a final center
/// spacing above the last center. A numerical collision solve imposes zero
/// flux at both faces. Energies are joules **per alpha**, amounts are moles.
#[derive(Clone, Debug, PartialEq)]
pub struct AlphaEnergyGrid {
    energies: Arc<[Energy]>,
    widths_j: Arc<[f64]>,
}

impl AlphaEnergyGrid {
    /// Retain an explicitly authored energy grid without sorting or repair.
    pub fn new(
        centers: Vec<Energy>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        poll(&mut cancelled)?;
        if !(2..=MAX_ALPHA_BINS).contains(&centers.len()) {
            return Err(DtFusionError::Invalid(
                "alpha energy grid requires 2..=256 centers",
            ));
        }
        for (i, center) in centers.iter().enumerate() {
            poll(&mut cancelled)?;
            if !positive(center.value()) || (i > 0 && *center <= centers[i - 1]) {
                return Err(DtFusionError::Invalid(
                    "alpha energy centers must be positive, finite and strictly increasing",
                ));
            }
        }
        let mut widths_j = Vec::with_capacity(centers.len());
        let mut left_face = 0.0;
        for i in 0..centers.len() {
            poll(&mut cancelled)?;
            let e = centers[i].value();
            let right_face = if i + 1 < centers.len() {
                e + 0.5 * (centers[i + 1].value() - e)
            } else {
                e + 0.5 * (e - centers[i - 1].value())
            };
            let width = right_face - left_face;
            if !positive(right_face) || !positive(width) {
                return Err(DtFusionError::Arithmetic(
                    "finite positive alpha energy cell widths",
                ));
            }
            widths_j.push(width);
            left_face = right_face;
        }
        poll(&mut cancelled)?;
        Ok(Self {
            energies: centers.into(),
            widths_j: widths_j.into(),
        })
    }

    /// Author logarithmically spaced positive centers, retaining both endpoints.
    /// The finite lower center is not a thermalization threshold; the lower
    /// cell still extends to zero energy with a reflecting face.
    pub fn logarithmic(
        minimum: Energy,
        maximum: Energy,
        count: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        poll(&mut cancelled)?;
        if !(2..=MAX_ALPHA_BINS).contains(&count)
            || !positive(minimum.value())
            || !positive(maximum.value())
            || minimum >= maximum
        {
            return Err(DtFusionError::Invalid(
                "positive ordered alpha energies and 2..=256 logarithmic centers required",
            ));
        }
        let logarithm = log_ratio(maximum.value(), minimum.value());
        let mut centers = Vec::with_capacity(count);
        for i in 0..count {
            poll(&mut cancelled)?;
            let energy = if i == 0 {
                minimum
            } else if i + 1 == count {
                maximum
            } else {
                let (anchor, exponent) = if i <= (count - 1) / 2 {
                    (minimum.value(), (i as f64 / (count - 1) as f64) * logarithm)
                } else {
                    (
                        maximum.value(),
                        -((count - 1 - i) as f64 / (count - 1) as f64) * logarithm,
                    )
                };
                let factor = exp(exponent);
                let product = anchor * factor;
                Energy::new(if factor.is_normal() && positive(product) {
                    product
                } else {
                    exp(ln(anchor) + exponent)
                })
            };
            centers.push(energy);
        }
        Self::new(centers, cancelled)
    }

    /// Author centers equally spaced in particle speed, equivalently sqrt(E).
    /// This usually resolves a thermal core and a fast birth with substantially
    /// less collision stiffness than extremely narrow low-energy log cells.
    /// Both energy endpoints are retained. For continuum refinement, also move
    /// the first center toward zero; a fixed positive first center retains a
    /// finite first-cell quadrature approximation even as more cells are added.
    pub fn uniform_speed(
        minimum: Energy,
        maximum: Energy,
        count: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        poll(&mut cancelled)?;
        if !(2..=MAX_ALPHA_BINS).contains(&count)
            || !positive(minimum.value())
            || !positive(maximum.value())
            || minimum >= maximum
        {
            return Err(DtFusionError::Invalid(
                "positive ordered alpha energies and 2..=256 speed-spaced centers required",
            ));
        }
        let low = sqrt(minimum.value());
        let high = sqrt(maximum.value());
        let gap = high - low;
        let mut centers = Vec::with_capacity(count);
        for i in 0..count {
            poll(&mut cancelled)?;
            let energy = if i == 0 {
                minimum
            } else if i + 1 == count {
                maximum
            } else {
                let root = if i <= (count - 1) / 2 {
                    low + gap * (i as f64 / (count - 1) as f64)
                } else {
                    high - gap * ((count - 1 - i) as f64 / (count - 1) as f64)
                };
                Energy::new(root * root)
            };
            centers.push(energy);
        }
        Self::new(centers, cancelled)
    }

    /// Retained cell centers, in joules per alpha.
    pub fn energies(&self) -> &[Energy] {
        &self.energies
    }

    /// Midpoint-face integration widths, in joules per alpha.
    pub fn widths_j(&self) -> &[f64] {
        &self.widths_j
    }

    /// Number of authored cells, at least two.
    pub fn len(&self) -> usize {
        self.energies.len()
    }

    /// Admitted grids are never empty.
    pub fn is_empty(&self) -> bool {
        self.energies.is_empty()
    }

    /// Place a monoenergetic birth between adjacent centers, preserving its
    /// number and first energy moment to floating-point precision. Refusal
    /// leaves the input population array unchanged; no edge clipping occurs.
    pub(crate) fn inject(
        &self,
        amounts: &mut [f64; MAX_ALPHA_BINS],
        amount: Amount,
        particle_energy: Energy,
    ) -> Result<(), DtFusionError> {
        self.admit_amounts(amounts)?;
        let (n, e) = (amount.value(), particle_energy.value());
        if !n.is_finite() || n < 0.0 || !positive(e) {
            return Err(DtFusionError::Invalid(
                "nonnegative finite alpha birth amount and positive finite particle energy",
            ));
        }
        if particle_energy < self.energies[0] || particle_energy > self.energies[self.len() - 1] {
            return Err(DtFusionError::OutsideDomain);
        }
        if n == 0.0 {
            return Ok(());
        }
        let high = self
            .energies
            .partition_point(|value| *value < particle_energy);
        let mut next = *amounts;
        if self.energies[high] == particle_energy {
            next[high] += n;
        } else {
            let low = high - 1;
            let from_low = e - self.energies[low].value();
            let from_high = self.energies[high].value() - e;
            let gap = self.energies[high].value() - self.energies[low].value();
            // Compute the smaller population directly; its complement retains
            // the full birth amount without subtracting nearby energy ratios.
            let (lower, upper) = if from_low <= from_high {
                let upper = n * (from_low / gap);
                (n - upper, upper)
            } else {
                let lower = n * (from_high / gap);
                (lower, n - lower)
            };
            next[low] += lower;
            next[high] += upper;
        }
        if next[..self.len()].iter().any(|n| !n.is_finite()) {
            return Err(DtFusionError::Arithmetic("finite alpha birth populations"));
        }
        *amounts = next;
        Ok(())
    }

    pub(super) fn admit_amounts(
        &self,
        amounts: &[f64; MAX_ALPHA_BINS],
    ) -> Result<(), DtFusionError> {
        if amounts[..self.len()]
            .iter()
            .any(|&n| !n.is_finite() || n < 0.0)
            || amounts[self.len()..].iter().any(|&n| n != 0.0)
        {
            return Err(DtFusionError::Invalid(
                "finite nonnegative alpha populations and zero unused cells required",
            ));
        }
        Ok(())
    }
}

pub(super) fn log_ratio(high: f64, low: f64) -> f64 {
    let ratio = high / low;
    if ratio.is_normal() {
        ln(ratio)
    } else {
        ln(high) - ln(low)
    }
}

#[cfg(test)]
mod tests;
