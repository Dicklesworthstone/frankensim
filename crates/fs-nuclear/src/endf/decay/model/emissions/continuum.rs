//! Closed positive moments of ENDF-102 (2023), §0.5.2.1, INT=3/4/5.
//!
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>
//!
//! INT=3 interpolates density in log energy, INT=4 log density in energy,
//! and INT=5 log density in log energy. Coefficients are formed once per source
//! segment, then reused across its group intersections. The formulas integrate
//! the supplied law without adaptive quadrature or a linearized replacement.
//! Positive endpoint bases avoid subtracting large count/energy primitives.

use fs_math::{
    decay::{attenuate, decayed_fraction},
    det,
};

use super::{DecaySpectrumPoint, EndfError, FissionError, PositiveIntegral, convex, invalid};

pub(super) struct LogSegment {
    a: DecaySpectrumPoint,
    b: DecaySpectrumPoint,
    law: u8,
    log_energy_span: f64,
    log_density_ratio: f64,
}

impl LogSegment {
    pub(super) fn new(
        a: DecaySpectrumPoint,
        b: DecaySpectrumPoint,
        law: u8,
    ) -> Result<Self, EndfError> {
        if matches!(law, 3 | 5) && a.energy_j <= 0.0 {
            return Err(invalid(
                "log-energy decay interpolation requires positive energies",
            ));
        }
        if matches!(law, 4 | 5) && (a.density_per_j <= 0.0 || b.density_per_j <= 0.0) {
            return Err(invalid(
                "log-density decay interpolation requires positive endpoint densities",
            ));
        }
        Ok(Self {
            a,
            b,
            law,
            log_energy_span: if matches!(law, 3 | 5) {
                log_ratio(b.energy_j, a.energy_j)
            } else {
                0.0
            },
            log_density_ratio: if matches!(law, 4 | 5) {
                if b.density_per_j >= a.density_per_j {
                    log_ratio(b.density_per_j, a.density_per_j)
                } else {
                    -log_ratio(a.density_per_j, b.density_per_j)
                }
            } else {
                0.0
            },
        })
    }

    pub(super) fn add_moments(
        &self,
        left: f64,
        right: f64,
        normalization: f64,
        count: &mut PositiveIntegral,
        energy: &mut PositiveIntegral,
    ) -> Result<(), EndfError> {
        let p = self.density_at(left)?;
        let q = self.density_at(right)?;
        let width = right - left;
        match self.law {
            3 => {
                let span = log_ratio(right, left);
                // For the measure E^n dE, the lower log-linear endpoint has
                // normalized weight 1/z - 1/expm1(z), z=(n+1) ln(right/left).
                // ∫E dE = width*(left+right)/2, retaining both energy factors
                // separately so their sum cannot overflow.
                let lower = lower_weight(span);
                count.add_triple_product(width, normalization, p, lower, 1.0)?;
                count.add_triple_product(width, normalization, q, 1.0 - lower, 1.0)?;
                let lower = lower_weight(2.0 * span);
                for (density, weight) in [(p, lower), (q, 1.0 - lower)] {
                    energy.add_triple_product(width, normalization, density, left, 2.0 / weight)?;
                    energy.add_triple_product(
                        width,
                        normalization,
                        density,
                        right,
                        2.0 / weight,
                    )?;
                }
            }
            4 => {
                let contrast =
                    self.log_density_ratio.abs() * (width / (self.b.energy_j - self.a.energy_j));
                let area = exponential_average(contrast);
                let toward_low_density = lower_weight(contrast);
                let (density, high_energy, low_energy) = if self.log_density_ratio >= 0.0 {
                    (q, right, left)
                } else {
                    (p, left, right)
                };
                count.add_product(width, normalization, density, 1.0 / area)?;
                energy.add_triple_product(
                    width,
                    normalization,
                    density,
                    high_energy,
                    1.0 / (area * (1.0 - toward_low_density)),
                )?;
                energy.add_triple_product(
                    width,
                    normalization,
                    density,
                    low_energy,
                    1.0 / (area * toward_low_density),
                )?;
            }
            5 => {
                let span = log_ratio(right, left);
                let fraction = span / self.log_energy_span;
                // In log energy, ∫E^n y(E)dE is an exponential integral.
                // Anchor to the larger E^(n+1)*y endpoint, identified from its
                // log contrast without ever forming those endpoint products.
                // This also handles the logarithmic limits y∝1/E and y∝1/E².
                let count_contrast = (self.log_density_ratio + self.log_energy_span) * fraction;
                let (density, endpoint) = if count_contrast >= 0.0 {
                    (q, right)
                } else {
                    (p, left)
                };
                count.add_triple_product(
                    span,
                    normalization,
                    density,
                    endpoint,
                    1.0 / exponential_average(count_contrast.abs()),
                )?;
                let energy_contrast =
                    (self.log_density_ratio + 2.0 * self.log_energy_span) * fraction;
                let (density, endpoint) = if energy_contrast >= 0.0 {
                    (q, right)
                } else {
                    (p, left)
                };
                energy.add_four_product(
                    span,
                    normalization,
                    density,
                    endpoint,
                    endpoint,
                    1.0 / exponential_average(energy_contrast.abs()),
                )?;
            }
            _ => unreachable!("only logarithmic interpolation enters LogSegment"),
        }
        Ok(())
    }

    fn density_at(&self, energy: f64) -> Result<f64, EndfError> {
        if energy == self.a.energy_j {
            return Ok(self.a.density_per_j);
        }
        if energy == self.b.energy_j {
            return Ok(self.b.density_per_j);
        }
        if self.a.density_per_j == self.b.density_per_j {
            return Ok(self.a.density_per_j);
        }
        let (lower, upper) = if self.law == 4 {
            let span = self.b.energy_j - self.a.energy_j;
            (
                (self.b.energy_j - energy) / span,
                (energy - self.a.energy_j) / span,
            )
        } else {
            let below = log_ratio(energy, self.a.energy_j);
            let above = log_ratio(self.b.energy_j, energy);
            // ln(E/a)+ln(b/E)=ln(b/a). Forming the two positive pieces keeps
            // each interpolation coordinate in [0,1] under rounding.
            (above / (below + above), below / (below + above))
        };
        if self.law == 3 {
            return convex(self.a.density_per_j, self.b.density_per_j, lower, upper);
        }
        let value = if self.log_density_ratio >= 0.0 {
            attenuate(self.b.density_per_j, self.log_density_ratio * lower)
        } else {
            attenuate(self.a.density_per_j, -self.log_density_ratio * upper)
        };
        // An interpolated subnormal would lose relative precision before the
        // normalization/energy factors can rescale it. Original endpoint data
        // remain admitted; a rounded intermediate density refuses explicitly.
        if value.is_normal() {
            Ok(value)
        } else {
            Err(FissionError::Arithmetic.into())
        }
    }
}

// ln(high/low), high>=low>0. A close-ratio log1p correction preserves the
// difference of neighboring values; distant ratios never overflow or underflow.
fn log_ratio(high: f64, low: f64) -> f64 {
    if high * 0.5 <= low {
        let delta = (high - low) / low;
        let shifted = 1.0 + delta;
        if shifted == 1.0 {
            delta
        } else {
            det::ln(shifted) * (delta / (shifted - 1.0))
        }
    } else {
        let ratio = high / low;
        if ratio.is_finite() {
            det::ln(ratio)
        } else {
            det::ln(high) - det::ln(low)
        }
    }
}

fn exponential_average(contrast: f64) -> f64 {
    if contrast == 0.0 {
        1.0
    } else {
        decayed_fraction(contrast) / contrast
    }
}

// ∫₀¹ t exp(-z t) dt / ∫₀¹ exp(-z t) dt. The Bernoulli expansion removes
// cancellation at z≈0; exp(-z) avoids overflow for large positive contrasts.
fn lower_weight(contrast: f64) -> f64 {
    if contrast <= 0.5 {
        let square = contrast * contrast;
        0.5 - contrast
            * (1.0 / 12.0
                - square
                    * (1.0 / 720.0
                        - square
                            * (1.0 / 30_240.0
                                - square
                                    * (1.0 / 1_209_600.0
                                        - square
                                            * (1.0 / 47_900_160.0
                                                - square
                                                    * (691.0 / 1_307_674_368_000.0
                                                        - square / 74_724_249_600.0))))))
    } else {
        1.0 / contrast - det::exp(-contrast) / decayed_fraction(contrast)
    }
}
