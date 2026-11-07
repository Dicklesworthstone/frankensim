//! Positive coherency interpolation for randomly located spherical scatterers.
//!
//! The two diagonal Jones powers and their complex cross product are retained
//! before conversion to I/Q. Convex mixtures of these positive semidefinite
//! 2x2 matrices describe independent spheres and angular/spectral interpolation;
//! averaging complex amplitudes would introduce nonexistent coherent particles.

use super::{DomainHasher, TabulatedPhase, VolumeTransportError, det};
use crate::polarization::{Mueller, Stokes};
use fs_math::c64::C64;
use std::sync::Arc;

/// Differential spherical-particle response in common right-handed s/p/travel
/// frames. +Q is perpendicular to the scattering plane, V=2 Im(Es conj(Ep)).
/// Powers are per steradian; the cross term is Es conj(Ep) in the same units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphericalPhaseValue {
    powers: [f64; 2],
    cross: C64,
}

impl SphericalPhaseValue {
    pub const ZERO: Self = Self {
        powers: [0.0; 2],
        cross: C64::ZERO,
    };

    /// Admit a positive semidefinite coherency response. Only roundoff-sized
    /// Cauchy excess is projected; a materially unphysical matrix refuses.
    pub fn try_new(powers: [f64; 2], mut cross: C64) -> Result<Self, VolumeTransportError> {
        if powers.iter().any(|p| !p.is_finite() || *p < 0.0)
            || !cross.re.is_finite()
            || !cross.im.is_finite()
        {
            return Err(VolumeTransportError::Phase);
        }
        let limit = det::sqrt(powers[0]) * det::sqrt(powers[1]);
        let magnitude = det::hypot(cross.re, cross.im);
        if limit == 0.0 {
            if magnitude != 0.0 {
                return Err(VolumeTransportError::Phase);
            }
        } else if magnitude > limit {
            if magnitude / limit > 1.0 + 128.0 * f64::EPSILON {
                return Err(VolumeTransportError::Phase);
            }
            cross = cross.scale(limit / magnitude);
        }
        Ok(Self { powers, cross })
    }

    pub const fn powers_per_sr(self) -> [f64; 2] {
        self.powers
    }
    pub const fn cross_per_sr(self) -> C64 {
        self.cross
    }
    pub fn unpolarized_per_sr(self) -> f64 {
        half_sum(self.powers[0], self.powers[1])
    }

    pub fn mueller(self) -> Mueller {
        let mean = self.unpolarized_per_sr();
        let difference = 0.5 * self.powers[0] - 0.5 * self.powers[1];
        Mueller([
            [mean, difference, 0.0, 0.0],
            [difference, mean, 0.0, 0.0],
            [0.0, 0.0, self.cross.re, -self.cross.im],
            [0.0, 0.0, self.cross.im, self.cross.re],
        ])
    }

    /// Apply retained positive channels, including exact pure-p or pure-s
    /// inputs whose weak output would disappear in a rounded M00-M01.
    pub fn apply(self, input: Stokes) -> Result<Stokes, VolumeTransportError> {
        if Stokes::try_new(input.i, input.q, input.u, input.v).is_none() {
            return Err(VolumeTransportError::Phase);
        }
        let s = channel_power(self.powers[0], input.i, input.q)?;
        let p = channel_power(self.powers[1], input.i, -input.q)?;
        let out = Stokes {
            i: s + p,
            q: s - p,
            u: self.cross.re * input.u - self.cross.im * input.v,
            v: self.cross.im * input.u + self.cross.re * input.v,
        };
        if Stokes::try_new(out.i, out.q, out.u, out.v).is_none() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(out)
    }

    pub(in crate::volumes::transport::scattering) fn scaled(
        self,
        scale: f64,
    ) -> Result<Self, VolumeTransportError> {
        if !scale.is_finite() || scale < 0.0 {
            return Err(VolumeTransportError::Phase);
        }
        Self::try_new(
            [
                product(self.powers[0], scale)?,
                product(self.powers[1], scale)?,
            ],
            C64::new(
                product(self.cross.re, scale)?,
                product(self.cross.im, scale)?,
            ),
        )
    }

    pub(in crate::volumes::transport::scattering) fn plus(
        self,
        rhs: Self,
    ) -> Result<Self, VolumeTransportError> {
        Self::try_new(
            [
                self.powers[0] + rhs.powers[0],
                self.powers[1] + rhs.powers[1],
            ],
            self.cross + rhs.cross,
        )
    }

    pub(in crate::volumes::transport::scattering) fn interpolate(
        self,
        rhs: Self,
        t: f64,
    ) -> Result<Self, VolumeTransportError> {
        if t == 0.0 {
            return Ok(self);
        }
        if t == 1.0 {
            return Ok(rhs);
        }
        self.scaled(1.0 - t)?.plus(rhs.scaled(t)?)
    }
}

/// Piecewise-linear positive coherency with the original scalar inverse CDF.
/// The same cosine knots describe all four retained real components. Convex
/// interpolation preserves realizability and includes incoherent depolarization.
#[derive(Clone, Debug, PartialEq)]
pub struct SphericalPhaseTable {
    scalar: TabulatedPhase,
    values: Arc<[SphericalPhaseValue]>,
}

impl SphericalPhaseTable {
    /// Normalize the complete response with its unpolarized solid-angle
    /// integral. The scalar owner enforces cosine support, knot/work bounds,
    /// positive sampling support and representability. This is an explicitly
    /// discretized angular law, not an inference from an intensity-only table.
    pub fn try_new(knots: &[(f64, SphericalPhaseValue)]) -> Result<Self, VolumeTransportError> {
        if !(2..=super::MAX_PHASE_KNOTS).contains(&knots.len()) {
            return Err(VolumeTransportError::Phase);
        }
        let scalar_knots: Vec<_> = knots
            .iter()
            .map(|(mu, v)| (*mu, v.unpolarized_per_sr()))
            .collect();
        let scalar = TabulatedPhase::try_new(&scalar_knots)?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(knots.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        for ((_, value), &(_, normalized)) in knots.iter().zip(scalar.knots_per_sr()) {
            let mean = value.unpolarized_per_sr();
            if mean == 0.0 {
                values.push(SphericalPhaseValue::ZERO);
                continue;
            }
            // Divide by the local mean before multiplying the normalized PDF:
            // channel/mean <= 2, even for tiny raw cross sections.
            let rescale = |x: f64| product_ratio(x, normalized, mean);
            values.push(SphericalPhaseValue::try_new(
                [rescale(value.powers[0])?, rescale(value.powers[1])?],
                C64::new(rescale(value.cross.re)?, rescale(value.cross.im)?),
            )?);
        }
        Ok(Self {
            scalar,
            values: values.into(),
        })
    }

    pub fn scalar_phase(&self) -> &TabulatedPhase {
        &self.scalar
    }
    pub fn values_per_sr(&self) -> &[SphericalPhaseValue] {
        &self.values
    }

    pub fn evaluate(&self, cosine: f64) -> Result<SphericalPhaseValue, VolumeTransportError> {
        self.scalar.eval_cos(cosine)?;
        let cosine = cosine.clamp(-1.0, 1.0);
        let knots = self.scalar.knots_per_sr();
        let upper = knots.partition_point(|&(mu, _)| mu <= cosine);
        let index = upper.saturating_sub(1).min(knots.len() - 2);
        let fraction =
            ((cosine - knots[index].0) / (knots[index + 1].0 - knots[index].0)).clamp(0.0, 1.0);
        self.values[index].interpolate(self.values[index + 1], fraction)
    }

    pub(super) fn hash_into(&self, hash: &mut DomainHasher) {
        self.scalar.hash_into(hash);
        for value in self.values.iter() {
            for x in [
                value.powers[0],
                value.powers[1],
                value.cross.re,
                value.cross.im,
            ] {
                hash.update(&x.to_bits().to_le_bytes());
            }
        }
    }
}

fn half_sum(a: f64, b: f64) -> f64 {
    let sum = a + b;
    if sum.is_finite() {
        0.5 * sum
    } else {
        0.5 * a + 0.5 * b
    }
}

fn channel_power(power: f64, intensity: f64, signed_q: f64) -> Result<f64, VolumeTransportError> {
    let half = half_sum(intensity, signed_q).max(0.0);
    if half >= f64::MIN_POSITIVE || intensity == -signed_q {
        return product(power, half);
    }
    // Halving a subnormal channel first can erase or substantially round a
    // representable result. Keep its exponent until after the phase power acts.
    product_ratio(power, (intensity + signed_q).max(0.0), 2.0)
}

// a*b/c with ordinary division-first arithmetic in the normal range. The rare
// exponent path avoids losing a weak channel before a large normalized angular
// density rescales it. It uses exact powers of two, without logarithmic error.
fn product_ratio(a: f64, b: f64, c: f64) -> Result<f64, VolumeTransportError> {
    if a == 0.0 || b == 0.0 {
        return Ok(a * b);
    }
    let quotient = a / c;
    let direct = quotient * b;
    if quotient.is_normal() && direct.is_finite() && direct != 0.0 {
        return Ok(direct);
    }
    let (a, ea) = binary_parts(a);
    let (b, eb) = binary_parts(b);
    let (c, ec) = binary_parts(c);
    let (mantissa, em) = binary_parts(a * b / c);
    let exponent = ea + eb - ec + em;
    let fraction_sign = mantissa.to_bits() & 0x800f_ffff_ffff_ffff;
    let value = if exponent > 1023 {
        return Err(VolumeTransportError::Unrepresentable);
    } else if exponent >= -1022 {
        f64::from_bits(fraction_sign | (((exponent + 1023) as u64) << 52))
    } else if exponent >= -1075 {
        // Round into the subnormal range once, including the half-ULP boundary.
        let normal = f64::from_bits(fraction_sign | (1_u64 << 52));
        normal * f64::from_bits(((exponent + 2045) as u64) << 52)
    } else {
        return Err(VolumeTransportError::Unrepresentable);
    };
    if !value.is_finite() || value == 0.0 {
        Err(VolumeTransportError::Unrepresentable)
    } else {
        Ok(value)
    }
}

// Finite nonzero x = mantissa * 2^exponent, 1 <= |mantissa| < 2.
fn binary_parts(x: f64) -> (f64, i32) {
    let bits = x.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent == 0 {
        let (mantissa, exponent) = binary_parts(x * 18_014_398_509_481_984.0);
        (mantissa, exponent - 54)
    } else {
        (
            f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52)),
            exponent - 1023,
        )
    }
}

impl super::PhaseFunction {
    /// Retained spherical coherency at one vacuum wavelength. Spectral tables
    /// interpolate the positive matrices; non-spherical laws refuse this API.
    pub fn spherical_value_at(
        &self,
        wavelength_nm: f64,
        cosine: f64,
    ) -> Result<SphericalPhaseValue, VolumeTransportError> {
        self.validate_wavelength(wavelength_nm)?;
        match self {
            Self::SphericalTable(table) => table.evaluate(cosine),
            Self::SpectralTable(table) => {
                let (a, b, t) = table.interval(wavelength_nm);
                if t == 0.0 {
                    return a.spherical_value_at(wavelength_nm, cosine);
                }
                if t == 1.0 {
                    return b.spherical_value_at(wavelength_nm, cosine);
                }
                a.spherical_value_at(wavelength_nm, cosine)?
                    .interpolate(b.spherical_value_at(wavelength_nm, cosine)?, t)
            }
            _ => Err(VolumeTransportError::Phase),
        }
    }
}

fn product(a: f64, b: f64) -> Result<f64, VolumeTransportError> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        Err(VolumeTransportError::Unrepresentable)
    } else {
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
