//! Original binary-factored Gram channels over a physical harmonic interval.

use super::{Sum, binary_scale, split};
use fs_fft::harmonic_window::{HarmonicWindow, WindowError};
use fs_math::{c64::C64, dd::Dd, det};

/// Shared analytic time weights, reused for all material channels and points.
/// Fields retain the original source clock and the separately authored duration.
pub(crate) struct HarmonicEnergyWindow {
    duration: f64,
    omega: f64,
    gram: [Dd; 2],
    start: [Dd; 2],
    center: [Dd; 2],
    end: [Dd; 2],
}

fn finite(value: f64) -> Result<f64, WindowError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WindowError::Unrepresentable)
    }
}

fn positive(value: f64) -> Result<f64, WindowError> {
    if value >= 0.0 {
        finite(value)
    } else {
        Err(WindowError::Unrepresentable)
    }
}

fn phase(angle: Dd) -> Result<[Dd; 2], WindowError> {
    finite(angle.hi)?;
    finite(angle.lo)?;
    let (s, c) = det::sin_cos_dd(angle);
    for part in [s.hi, s.lo, c.hi, c.lo] {
        finite(part)?;
    }
    Ok([c, -s])
}

fn phase_product(a: [Dd; 2], b: [Dd; 2]) -> [Dd; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

// Keep each channel's binary exponent outside products, including the
// temporal weights and original tau. A rounded spring strain, squared strain,
// strain rate, or reciprocal tau is never an intermediate physical output.
fn normalized(value: Sum) -> Sum {
    if value.fraction == Dd::ZERO {
        return Sum::default();
    }
    let (_, shift) = split(if value.fraction.hi != 0.0 {
        value.fraction.hi
    } else {
        value.fraction.lo
    });
    Sum {
        fraction: Dd::from_pair(
            binary_scale(value.fraction.hi, -shift),
            binary_scale(value.fraction.lo, -shift),
        ),
        exponent: value.exponent + shift,
    }
}

fn product(a: Sum, b: Sum) -> Sum {
    if a.fraction == Dd::ZERO || b.fraction == Dd::ZERO {
        return Sum::default();
    }
    let (a, b) = (normalized(a), normalized(b));
    Sum {
        fraction: a.fraction * b.fraction,
        exponent: a.exponent + b.exponent,
    }
}

fn dd(value: Dd) -> Sum {
    normalized(Sum {
        fraction: value,
        exponent: 0,
    })
}

fn scaled(value: Sum, factor: Dd) -> Sum {
    product(value, dd(factor))
}

fn add(target: &mut Sum, value: Sum) {
    if value.fraction != Dd::ZERO {
        target.accumulate(value.fraction, value.exponent);
    }
}

fn rotate(z: [Sum; 2], phase: [Dd; 2]) -> [Sum; 2] {
    let (mut re, mut im) = (Sum::default(), Sum::default());
    add(&mut re, scaled(z[0], phase[0]));
    add(&mut re, scaled(z[1], -phase[1]));
    add(&mut im, scaled(z[0], phase[1]));
    add(&mut im, scaled(z[1], phase[0]));
    [re, im]
}

fn weighted(target: &mut Sum, a: Sum, b: Sum, weight: Dd, numerator: &[f64], denominator: &[f64]) {
    let mut factor = Sum::default();
    factor.add_ratio(numerator, denominator);
    add(target, product(product(a, b), product(dd(weight), factor)));
}

fn project(row: &[f64; 6], strain: &[C64; 6], spectrum: [(f64, i32); 2]) -> [Sum; 2] {
    let [a, b] = spectrum;
    let (mut re, mut im) = (Sum::default(), Sum::default());
    for (&map, value) in row.iter().zip(strain) {
        re.add_binary_ratio(&[map, value.re, a.0], &[], a.1);
        re.add_binary_ratio(&[map, value.im, b.0], &[], b.1);
        im.add_binary_ratio(&[map, value.im, a.0], &[], a.1);
        im.add_binary_ratio(&[-map, value.re, b.0], &[], b.1);
    }
    [re, im]
}

impl HarmonicEnergyWindow {
    pub(crate) fn new(omega: f64, window: HarmonicWindow) -> Result<Self, WindowError> {
        let gram = window.product_kernel(omega, omega)?.centered_gram_dd();
        let w = Dd::from_f64(omega);
        let start = phase(w * Dd::from_f64(window.start_s()))?;
        // Apply half before multiplying the larger factor, as the shared
        // temporal owner does. Keep start/width rotations separate even when
        // adding the duration to a large timestamp would round it away.
        let half_width = if omega.abs() >= window.duration_s() {
            (w * Dd::from_f64(0.5)) * Dd::from_f64(window.duration_s())
        } else {
            w * (Dd::from_f64(window.duration_s()) * Dd::from_f64(0.5))
        };
        let half = phase(half_width)?;
        Ok(Self {
            duration: window.duration_s(),
            omega,
            gram,
            start,
            center: phase_product(start, half),
            end: phase_product(start, phase_product(half, half)),
        })
    }

    /// Independent original spring storage at both endpoints, integrated
    /// storage rate, and (when tau is supplied) physical Maxwell heat.
    pub(crate) fn spring(
        &self,
        gram: &[[f64; 6]; 6],
        strain: &[C64; 6],
        spectrum: [(f64, i32); 2],
        tau: Option<f64>,
    ) -> Result<[f64; 4], WindowError> {
        let mut totals = [Sum::default(); 4];
        for row in gram {
            let z = project(row, strain, spectrum);
            let start = rotate(z, self.start)[0];
            let end = rotate(z, self.end)[0];
            weighted(&mut totals[0], start, start, Dd::ONE, &[0.5], &[]);
            weighted(&mut totals[1], end, end, Dd::ONE, &[0.5], &[]);
            let [re, im] = rotate(z, self.center);
            // Integral q*(dq/dt) is the spring energy change. Evaluating
            // the original storage derivative avoids subtracting two rounded
            // endpoint energies on an extremely short interval.
            weighted(
                &mut totals[2],
                re,
                im,
                self.gram[0] - self.gram[1],
                &[self.omega, self.duration],
                &[],
            );
            if let Some(tau) = tau {
                weighted(
                    &mut totals[3],
                    re,
                    re,
                    self.gram[0],
                    &[self.duration],
                    &[tau],
                );
                weighted(
                    &mut totals[3],
                    im,
                    im,
                    self.gram[1],
                    &[self.duration],
                    &[tau],
                );
            }
        }
        Ok([
            positive(totals[0].value())?,
            positive(totals[1].value())?,
            finite(totals[2].value())?,
            positive(totals[3].value())?,
        ])
    }

    pub(crate) fn kelvin_heat(
        &self,
        gram: &[[f64; 6]; 6],
        strain: &[C64; 6],
    ) -> Result<f64, WindowError> {
        let mut total = Sum::default();
        for row in gram {
            let rate = project(row, strain, [(0.0, 0), (self.omega, 0)]);
            let [re, im] = rotate(rate, self.center);
            weighted(&mut total, re, re, self.gram[0], &[self.duration], &[]);
            weighted(&mut total, im, im, self.gram[1], &[self.duration], &[]);
        }
        positive(total.value())
    }

    /// Accumulate one original stress channel times strain rate. Original
    /// coefficients, spectrum and rate stay factored until all channels have
    /// supplied their work; even an unrepresentable complex-modulus entry can
    /// yield representable interval work. No storage or heat enforces closure.
    pub(crate) fn add_stress_work(
        &self,
        matrix: &[[f64; 6]; 6],
        strain: &[C64; 6],
        spectrum: [(f64, i32); 2],
        total: &mut Sum,
    ) {
        for (i, row) in matrix.iter().enumerate() {
            let stress = rotate(project(row, strain, spectrum), self.center);
            let mut rate = [Sum::default(); 2];
            rate[0].add([self.omega, strain[i].im, 1.0, 1.0]);
            rate[1].add([-self.omega, strain[i].re, 1.0, 1.0]);
            let rate = rotate(rate, self.center);
            for j in 0..2 {
                weighted(
                    total,
                    stress[j],
                    rate[j],
                    self.gram[j],
                    &[self.duration],
                    &[],
                );
            }
        }
    }
}
