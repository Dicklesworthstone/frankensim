//! Positive pressure covariance for stationary two-sided porous excitation.
//! C_ij = E[u_i conj(u_j)] uses RMS complex amplitudes: a deterministic peak
//! phasor p has C = p p*/2. C_LR is the LEFT amplitude times conjugate RIGHT.
//! Components represent one frequency, never extra independent frequencies.

use super::{DuctError, PorousLayerStack, PorousStackExcitation, bad, finite};
use fs_math::{c64::C64, dd::Dd, det};

/// Admitted Hermitian positive-semidefinite 2×2 RMS pressure covariance [Pa²].
/// A deterministic positive factor supplies at most two original peak phasor
/// pairs. Their native cycle-mean fields and losses add without a subtraction
/// of large boundary powers or an eigenvalue clamp.
#[derive(Clone, Copy, Debug)]
pub struct PressureCovariance {
    matrix_pa2: [C64; 4],
    components: [[C64; 2]; 2],
    count: usize,
}

// An exact binary scale near sqrt(x), including subnormal x. Dividing x by
// this scale twice normalizes its diagonal to [1,4) without under/overflow.
fn root_scale(x: f64) -> f64 {
    let bits = x.to_bits();
    let encoded = ((bits >> 52) & 0x7ff) as i32;
    let exponent = if encoded == 0 {
        -1074 + (63 - bits.leading_zeros()) as i32
    } else {
        encoded - 1023
    };
    f64::from_bits(((exponent.div_euclid(2) + 1023) as u64) << 52)
}
fn real(value: Dd) -> Result<f64, DuctError> {
    let value = value.to_f64();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("pressure covariance or factor is unrepresentable"))
    }
}

// Exact sign of a*b-re²-im² for the ORIGINAL encoded f64 values. A product
// has at most 106 significand bits and exponent -2148..1942; 66 fixed limbs
// cover the entire range including the carry from two squared components.
// This matters at rank one: even an imaginary part whose square underflows
// f64 makes [1,1+i*epsilon;1-i*epsilon,1] strictly indefinite.
fn determinant(diagonal: [f64; 2], cross: C64, scales: [f64; 2]) -> Result<Dd, DuctError> {
    fn parts(value: f64) -> (u64, i32) {
        let bits = value.to_bits() & 0x7fff_ffff_ffff_ffff;
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        let mantissa = bits & 0x000f_ffff_ffff_ffff;
        if exponent == 0 {
            (mantissa, -1074)
        } else {
            (mantissa | (1_u64 << 52), exponent - 1075)
        }
    }
    fn add_product(words: &mut [u64; 66], a: f64, b: f64) {
        let (a, ea) = parts(a);
        let (b, eb) = parts(b);
        let product = u128::from(a) * u128::from(b);
        let offset = (ea + eb + 2148) as usize;
        let shift = offset % 64;
        let low = product as u64;
        let high = (product >> 64) as u64;
        let chunks = if shift == 0 {
            [low, high, 0]
        } else {
            [
                low << shift,
                (high << shift) | (low >> (64 - shift)),
                high >> (64 - shift),
            ]
        };
        let mut carry = 0_u128;
        for (i, word) in words.iter_mut().enumerate().skip(offset / 64) {
            let chunk = chunks.get(i - offset / 64).copied().unwrap_or(0);
            let sum = u128::from(*word) + u128::from(chunk) + carry;
            *word = sum as u64;
            carry = sum >> 64;
            if i >= offset / 64 + 2 && carry == 0 {
                break;
            }
        }
    }
    let mut lhs = [0_u64; 66];
    let mut rhs = [0_u64; 66];
    add_product(&mut lhs, diagonal[0], diagonal[1]);
    add_product(&mut rhs, cross.re, cross.re);
    add_product(&mut rhs, cross.im, cross.im);
    match lhs.iter().rev().cmp(rhs.iter().rev()) {
        core::cmp::Ordering::Less => {
            return Err(bad("pressure covariance is not positive semidefinite"));
        }
        core::cmp::Ordering::Equal => return Ok(Dd::ZERO),
        core::cmp::Ordering::Greater => {}
    }
    let mut borrow = false;
    for (a, b) in lhs.iter_mut().zip(rhs) {
        let (difference, first) = a.overflowing_sub(b);
        let (difference, second) = difference.overflowing_sub(u64::from(borrow));
        *a = difference;
        borrow = first || second;
    }
    let top = lhs
        .iter()
        .rposition(|&v| v != 0)
        .expect("positive exact determinant");
    let scale_exponent = scales
        .map(|s| ((s.to_bits() >> 52) as i32) - 1023)
        .into_iter()
        .sum::<i32>()
        * 2;
    let mut value = Dd::ZERO;
    // Three leading limbs retain at least 129 bits, beyond DD precision. All
    // terms are positive, so conversion cannot erase a small Schur remainder.
    for index in (top.saturating_sub(2)..=top).rev() {
        let exponent = -2148 + 64 * index as i32 - scale_exponent;
        let high = (lhs[index] >> 32) as f64;
        let low = (lhs[index] & 0xffff_ffff) as f64;
        value = value
            + Dd::from_f64(high * 2.0_f64.powi(exponent + 32))
            + Dd::from_f64(low * 2.0_f64.powi(exponent));
    }
    Ok(value)
}

impl PressureCovariance {
    /// Admit a row-major matrix [C_LL,C_LR,C_RL,C_RR]. Diagonals must be finite
    /// real nonnegative values; off-diagonals must be exact conjugates.
    /// An indefinite encoded matrix refuses, including a negative near-zero
    /// determinant. No correlation is projected or negative eigenvalue clipped.
    /// Use `from_peak_pressure` when the known source is exactly coherent;
    /// its positive witness avoids separately rounded outer-product entries.
    pub fn try_new(matrix_pa2: [C64; 4]) -> Result<Self, DuctError> {
        if !matrix_pa2.into_iter().all(finite)
            || matrix_pa2[0].im != 0.0
            || matrix_pa2[3].im != 0.0
            || matrix_pa2[0].re < 0.0
            || matrix_pa2[3].re < 0.0
            || matrix_pa2[1] != matrix_pa2[2].conj()
        {
            return Err(bad(
                "pressure covariance must be finite and Hermitian with nonnegative diagonal",
            ));
        }
        let diagonal = [matrix_pa2[0].re, matrix_pa2[3].re];
        let cross = matrix_pa2[1];
        let mut components = [[C64::ZERO; 2]; 2];
        if diagonal.contains(&0.0) {
            if cross != C64::ZERO {
                return Err(bad("zero pressure variance requires zero cross covariance"));
            }
            for i in 0..2 {
                components[0][i] = C64::from_re(det::sqrt(diagonal[i]) * core::f64::consts::SQRT_2);
            }
            return Ok(Self {
                matrix_pa2,
                components,
                count: 1,
            });
        }
        let scales = diagonal.map(root_scale);
        let normalized = std::array::from_fn::<_, 2, _>(|i| diagonal[i] / scales[i] / scales[i]);
        let d = normalized.map(Dd::from_f64);
        let determinant = determinant(diagonal, cross, scales)?;
        // Pivot on the larger PHYSICAL variance. The scaled DD determinant
        // retains a small positive second component beside a strong coherent
        // component; direct b-|c|²/a subtraction would lose that information.
        let pivot = usize::from(diagonal[1] > diagonal[0]);
        let other = 1 - pivot;
        let root = d[pivot].sqrt() * Dd::from_f64(scales[pivot]);
        let sqrt_two = Dd::from_f64(2.0).sqrt();
        components[0][pivot] = C64::from_re(real(root * sqrt_two)?);
        let correlated = if pivot == 0 { cross.conj() } else { cross };
        components[0][other] = C64::new(
            real((Dd::from_f64(correlated.re) / root) * sqrt_two)?,
            real((Dd::from_f64(correlated.im) / root) * sqrt_two)?,
        );
        let remainder =
            (determinant / d[pivot] * Dd::from_f64(2.0)).sqrt() * Dd::from_f64(scales[other]);
        components[1][other] = C64::from_re(real(remainder)?);
        let count = if determinant == Dd::ZERO { 1 } else { 2 };
        Ok(Self {
            matrix_pa2,
            components,
            count,
        })
    }

    /// Construct a known coherent source without testing the determinant of
    /// separately rounded p p*/2 entries. The ORIGINAL peak phasors remain the
    /// exact positive witness, preserving their phase and native arithmetic.
    pub fn from_peak_pressure(pressure_pa: [C64; 2]) -> Result<Self, DuctError> {
        if !pressure_pa.into_iter().all(finite) {
            return Err(bad("coherent covariance pressure must be finite"));
        }
        let p = pressure_pa;
        let h = p.map(|value| value.scale(0.5));
        let mul = |a, b| Dd::from_f64(a) * Dd::from_f64(b);
        let a = real(mul(h[0].re, p[0].re) + mul(h[0].im, p[0].im))?;
        let b = real(mul(h[1].re, p[1].re) + mul(h[1].im, p[1].im))?;
        let c = C64::new(
            real(mul(h[0].re, p[1].re) + mul(h[0].im, p[1].im))?,
            real(mul(h[0].im, p[1].re) - mul(h[0].re, p[1].im))?,
        );
        Ok(Self {
            matrix_pa2: [C64::from_re(a), c, c.conj(), C64::from_re(b)],
            components: [pressure_pa, [C64::ZERO; 2]],
            count: 1,
        })
    }

    /// One-sided cross-spectral density per rad/s, with C_bin = S(omega) * w.
    /// `quadrature_weight_rad_s` includes the positive authored integration
    /// weight. Mean-square pressure is integral_0^infinity S(omega) d omega;
    /// there is no extra 1/2 or 2 pi. Per-Hz matrices require S_omega=S_Hz/(2 pi)
    /// and w_omega=2 pi*w_Hz. PSD density is admitted BEFORE positive scaling,
    /// then its positive factor is retained through ordinary rounding.
    pub fn from_angular_density(
        density_pa2_per_rad_s: [C64; 4],
        quadrature_weight_rad_s: f64,
    ) -> Result<Self, DuctError> {
        if !quadrature_weight_rad_s.is_finite() || quadrature_weight_rad_s <= 0.0 {
            return Err(bad(
                "pressure cross-spectrum needs a positive finite angular weight",
            ));
        }
        let mut source = Self::try_new(density_pa2_per_rad_s)?;
        source.matrix_pa2 = source
            .matrix_pa2
            .map(|value| value.scale(quadrature_weight_rad_s));
        let scale = det::sqrt(quadrature_weight_rad_s);
        source.components = source
            .components
            .map(|pair| pair.map(|value| value.scale(scale)));
        if !source.matrix_pa2.into_iter().all(finite)
            || !source.components.into_iter().flatten().all(finite)
        {
            return Err(bad("weighted pressure cross-spectrum is unrepresentable"));
        }
        Ok(source)
    }

    /// RMS covariance in row-major order. Entries are rounded; constructors
    /// with a known coherent/density witness retain that positive witness.
    pub const fn matrix_pa2(&self) -> [C64; 4] {
        self.matrix_pa2
    }
    /// At most two mutually uncorrelated PEAK phasor pairs at THIS frequency.
    /// Each pair's two ports remain coherent. Sum native powers/mean squares,
    /// never the complex pressures. Zero source uses one all-zero pair.
    pub fn peak_components(&self) -> &[[C64; 2]] {
        &self.components[..self.count]
    }
}

impl<'l, 'm> PorousLayerStack<'l, 'm> {
    /// Prepare one covariance on the already-composed original stack.
    /// At most two original coherent excitations reuse its same incident maps;
    /// no second composition or material evaluation is performed. The returned
    /// component fields/integrals add statistically at this one frequency.
    pub fn prepare_covariance(
        &self,
        covariance: PressureCovariance,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<PorousStackExcitation<'_, 'l, 'm>>, DuctError> {
        let mut components = Vec::with_capacity(covariance.peak_components().len());
        for &pressure in covariance.peak_components() {
            components.push(self.prepare_excitation(pressure, &mut cancelled)?);
        }
        Ok(components)
    }
}

#[cfg(test)]
mod tests;
