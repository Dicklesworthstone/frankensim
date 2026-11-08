//! Spatial corrections composed with the existing real harmonic-product kernel.
//!
//! For pressure polynomials C_a and C_b, H = <C_a conj(C_b)> and
//! S = <C_a C_b>. The four real spatial products are
//! (Re H + Re S, Im S - Im H, Im S + Im H, Re H - Re S)/2.
//! Only separation corrections are cached: their colocated limits are applied
//! directly through the temporal Gram kernel, preserving small sine weights.

use super::{
    AeroacError, Atom, ComplexDd, Dd, RadiatingMultipole, finite_dd, nonfinite, poll,
    window_error,
};
use super::super::{PairMoment, TensorPairMoment};
use super::super::momentum::{FifthMoment, correction};
use fs_fft::harmonic_window::HarmonicProduct;
use fs_math::special::spherical_bessel_moments;

#[derive(Clone, Copy, Default)]
pub(super) struct SpatialCorrection {
    conjugate: [ComplexDd; 4],
    unconjugated: [ComplexDd; 4],
}

impl SpatialCorrection {
    pub(super) fn apply(self, temporal: HarmonicProduct) -> Result<[Dd; 4], AeroacError> {
        let real = [Dd::ONE, Dd::ZERO];
        let imag = [Dd::ZERO, Dd::ONE];
        let rr = temporal.mean_real_product_parts_dd(real, real).map_err(window_error)?;
        let ri = temporal.mean_real_product_parts_dd(real, imag).map_err(window_error)?;
        let ir = temporal.mean_real_product_parts_dd(imag, real).map_err(window_error)?;
        let ii = temporal.mean_real_product_parts_dd(imag, imag).map_err(window_error)?;
        let half = Dd::from_f64(0.5);
        let weights = [(rr + ii) * half, (ir - ri) * half, (rr - ii) * half, (ir + ri) * half];
        Ok(core::array::from_fn(|i| {
            self.conjugate[i].re * weights[0] + self.conjugate[i].im * weights[1]
                + self.unconjugated[i].re * weights[2] + self.unconjugated[i].im * weights[3]
        }))
    }
}

pub(super) fn prepare_correction(
    a: &[Atom],
    b: &[Atom],
    same_frequency: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialCorrection, AeroacError> {
    let mut result = SpatialCorrection::default();
    for (i, left) in a.iter().enumerate() {
        let count = if same_frequency { i + 1 } else { b.len() };
        for (j, right) in b[..count].iter().enumerate() {
            poll(cancelled)?;
            let weight = Dd::from_f64(if same_frequency && i != j { 2.0 } else { 1.0 });
            let conjugate = if same_frequency && i == j {
                [ComplexDd::default(); 4]
            } else {
                angular_correction(left.source, right.source,
                    core::array::from_fn(|axis| left.wavevector[axis] - right.wavevector[axis]))?
            };
            let unconjugated = angular_correction(
                left.source,
                transform(right.source, |v| ComplexDd { re: v.re, im: -v.im }),
                core::array::from_fn(|axis| left.wavevector[axis] + right.wavevector[axis]),
            )?;
            for axis in 0..4 {
                result.conjugate[axis].re =
                    result.conjugate[axis].re + weight * conjugate[axis].re;
                // Within one frequency, reverse-order Hermitian terms are
                // exact conjugates. Their imaginary parts cancel analytically.
                if !same_frequency {
                    result.conjugate[axis].im =
                        result.conjugate[axis].im + conjugate[axis].im;
                }
                result.unconjugated[axis].re =
                    result.unconjugated[axis].re + weight * unconjugated[axis].re;
                result.unconjugated[axis].im =
                    result.unconjugated[axis].im + weight * unconjugated[axis].im;
            }
        }
    }
    if result.conjugate.iter().chain(&result.unconjugated).any(|v| !v.finite()) {
        return Err(nonfinite("radiated-energy angular interference"));
    }
    Ok(result)
}

pub(super) struct AngularKernel {
    pair: PairMoment,
    tensor: TensorPairMoment,
    fifth: FifthMoment,
}

impl AngularKernel {
    fn apply(&self, a: RadiatingMultipole, b: RadiatingMultipole) -> [Dd; 4] {
        let power = self.tensor.correction(self.pair, a, b);
        let moment = correction(self.pair, self.tensor, self.fifth, a, b);
        [power, moment[0], moment[1], moment[2]]
    }
}

fn angular_correction(
    a: RadiatingMultipole,
    b: RadiatingMultipole,
    phase: [Dd; 3],
) -> Result<[ComplexDd; 4], AeroacError> {
    if phase.iter().any(|v| !finite_dd(*v)) {
        return Err(nonfinite("radiated-energy pair phase"));
    }
    let scale = phase.iter().fold(0.0_f64, |s, v| s.max(v.hi.abs()).max(v.lo.abs()));
    if scale == 0.0 {
        return Ok([ComplexDd::default(); 4]);
    }
    let scaled = phase.map(|v| v / Dd::from_f64(scale));
    let norm = Dd::from_f64(scale)
        * scaled.into_iter().fold(Dd::ZERO, |sum, v| sum + v * v).sqrt();
    if !finite_dd(norm) || norm.to_f64() <= 0.0 || !norm.to_f64().is_finite() {
        return Err(nonfinite("radiated-energy phase separation"));
    }
    let x = norm.to_f64();
    let pair = PairMoment {
        direction: phase.map(|v| (v / norm).to_f64()),
        moments: spherical_bessel_moments(x)
            .ok_or_else(|| nonfinite("radiated-energy spherical moments"))?,
    };
    let tensor = TensorPairMoment::new(x, pair.moments)?;
    let kernel = AngularKernel { pair, tensor, fifth: FifthMoment::new(x, tensor)? };
    let real = kernel.apply(a, b);
    // Re[C_a conj(i C_b) exp(-i n.v)] = Im[H]. The identical
    // operation also recovers the unconjugated product after conjugating b.
    let imaginary = kernel.apply(a, transform(b, |v| ComplexDd { re: -v.im, im: v.re }));
    Ok(core::array::from_fn(|i| ComplexDd { re: real[i], im: imaginary[i] }))
}

pub(super) fn transform(
    mut source: RadiatingMultipole,
    map: impl Fn(ComplexDd) -> ComplexDd + Copy,
) -> RadiatingMultipole {
    source.compact.volume = map(source.compact.volume);
    source.compact.force = source.compact.force.map(map);
    source.tensor = source.tensor.map(map);
    source
}

pub(super) fn colocated_product(
    kernel: HarmonicProduct,
    a: RadiatingMultipole,
    b: RadiatingMultipole,
) -> Result<[Dd; 4], AeroacError> {
    let product = |a: ComplexDd, b: ComplexDd| {
        kernel.mean_real_product_parts_dd([a.re, a.im], [b.re, b.im]).map_err(window_error)
    };
    let third = Dd::ONE / Dd::from_f64(3.0);
    let two_fifteenths = Dd::from_f64(2.0) / Dd::from_f64(15.0);
    let mut result = [Dd::ZERO; 4];
    result[0] = product(a.compact.volume, b.compact.volume)?;
    for axis in 0..3 {
        result[0] = result[0] + third * product(a.compact.force[axis], b.compact.force[axis])?;
    }
    for axis in 0..6 {
        let value = product(a.tensor[axis], b.tensor[axis])?;
        result[0] = result[0] + two_fifteenths * if axis < 3 { value } else { value + value };
    }
    for (axis, row) in [[0, 3, 4], [3, 1, 5], [4, 5, 2]].into_iter().enumerate() {
        result[axis + 1] = third * (product(a.compact.volume, b.compact.force[axis])?
            + product(a.compact.force[axis], b.compact.volume)?);
        for (j, tensor) in row.into_iter().enumerate() {
            result[axis + 1] = result[axis + 1] + two_fifteenths
                * (product(a.tensor[tensor], b.compact.force[j])?
                    + product(a.compact.force[j], b.tensor[tensor])?);
        }
    }
    Ok(result)
}
