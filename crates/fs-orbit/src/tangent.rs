//! Exact coefficient derivatives of the original sampled HB balance.
//!
//! If A_l is the Fourier coefficient of the time-local island Jacobian,
//! a real X_k column contributes A_(n-k) + A_(n+k), and an imaginary
//! column contributes i(A_(n-k) - A_(n+k)). A negative index conjugates.
//! DC has only one real coefficient. These are derivatives of the discrete
//! AFT sum, not an assumption that a nonlinear Jacobian is band-limited.

use super::{AftWorkspace, IslandStateJacobian, OrbitError, OrbitProblem, Packing, TAU};
use fs_math::{c64::C64, det};

fn invalid(what: &'static str) -> OrbitError {
    OrbitError::BadParameter { what }
}

fn zeros<T: Clone>(count: usize, zero: T) -> Result<Vec<T>, OrbitError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("analytic AFT workspace is unrepresentable or unavailable"))?;
    values.resize(count, zero);
    Ok(values)
}

fn product(a: usize, b: usize) -> Result<usize, OrbitError> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("analytic AFT workspace extent overflows"))
}

/// Reused at every Newton iterate; phase tables change only on AFT refinement.
pub(super) struct Workspace {
    square: usize,
    orders: usize,
    samples: usize,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    local: Vec<f64>,
    coefficients: Vec<C64>,
    correction: Vec<C64>,
    ports: Vec<C64>,
}

impl Workspace {
    pub(super) fn new(pack: &Packing) -> Result<Self, OrbitError> {
        let square = product(pack.d, pack.d)?;
        let orders = pack
            .n_harm
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| invalid("analytic AFT order overflows"))?;
        let entries = product(square, orders)?;
        let port_entries = product(square, pack.n_harm + 1)?;
        // Check the combined numerical payload before allocating any component.
        let complex = entries
            .checked_mul(2)
            .and_then(|n| n.checked_add(port_entries))
            .and_then(|n| n.checked_mul(core::mem::size_of::<C64>()));
        let real = product(square, core::mem::size_of::<f64>())?;
        if complex
            .and_then(|n| n.checked_add(real))
            .is_none_or(|n| n > isize::MAX as usize)
        {
            return Err(invalid("analytic AFT numerical payload overflows"));
        }
        Ok(Self {
            square,
            orders,
            samples: 0,
            cosine: Vec::new(),
            sine: Vec::new(),
            local: zeros(square, 0.0)?,
            coefficients: zeros(entries, C64::ZERO)?,
            correction: zeros(entries, C64::ZERO)?,
            ports: zeros(port_entries, C64::ZERO)?,
        })
    }

    fn prepare_basis(&mut self, samples: usize) -> Result<(), OrbitError> {
        if self.samples == samples {
            return Ok(());
        }
        let count = product(samples, self.orders)?;
        let bytes = product(count, 2 * core::mem::size_of::<f64>())?;
        if bytes > isize::MAX as usize {
            return Err(invalid("analytic AFT phase-table payload overflows"));
        }
        let mut cosine = zeros(count, 0.0)?;
        let mut sine = zeros(count, 0.0)?;
        for sample in 0..samples {
            let theta = TAU * sample as f64 / samples as f64;
            cosine[sample * self.orders] = 1.0;
            for order in 1..self.orders {
                let (s, c) = det::sin_cos(order as f64 * theta);
                sine[sample * self.orders + order] = s;
                cosine[sample * self.orders + order] = c;
            }
        }
        self.cosine = cosine;
        self.sine = sine;
        self.samples = samples;
        Ok(())
    }

    /// `aft.xs` must still contain this iterate's original physical waveform.
    pub(super) fn prepare<P: OrbitProblem>(
        &mut self,
        problem: &P,
        callback: IslandStateJacobian<P>,
        pack: &Packing,
        omega: f64,
        aft: &AftWorkspace,
    ) -> Result<(), OrbitError> {
        self.prepare_basis(aft.samples)?;
        self.coefficients.fill(C64::ZERO);
        self.correction.fill(C64::ZERO);
        let period = TAU / omega;
        for sample in 0..aft.samples {
            self.local.fill(f64::NAN);
            let time = period * sample as f64 / aft.samples as f64;
            if !callback(
                problem,
                time,
                &aft.xs[sample * pack.d..(sample + 1) * pack.d],
                &mut self.local,
            ) || self.local.iter().any(|value| !value.is_finite())
            {
                return Err(OrbitError::TangentRefused { sample });
            }
            // Normalize before accumulation to avoid a gratuitous sample-count
            // overflow. Fixed sample order and compensated sums preserve weak
            // harmonic entries beside a much larger constant tangent.
            for order in 0..self.orders {
                let index = sample * self.orders + order;
                let c = self.cosine[index] / aft.samples as f64;
                let s = -self.sine[index] / aft.samples as f64;
                for entry in 0..self.square {
                    let index = order * self.square + entry;
                    let term = C64::new(c * self.local[entry], s * self.local[entry]);
                    let adjusted = term - self.correction[index];
                    let next = self.coefficients[index] + adjusted;
                    self.correction[index] = (next - self.coefficients[index]) - adjusted;
                    self.coefficients[index] = next;
                }
            }
        }
        if self
            .coefficients
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid("analytic AFT Fourier Jacobian must be finite"));
        }
        for order in 0..=pack.n_harm {
            let port = problem.port(C64::new(0.0, order as f64 * omega));
            if port.len() != self.square
                || port.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            {
                return Err(invalid(
                    "analytic AFT port must be a complete finite matrix",
                ));
            }
            self.ports[order * self.square..(order + 1) * self.square].copy_from_slice(&port);
        }
        Ok(())
    }

    fn fourier(&self, order: usize, entry: usize) -> C64 {
        self.coefficients[order * self.square + entry]
    }

    /// One original physical derivative, before row/column conditioning.
    pub(super) fn coefficient_entry(&self, row: usize, col: usize, pack: &Packing) -> f64 {
        let (n, output, imaginary_row) = coordinate(row, pack.d);
        let (k, input, imaginary_col) = coordinate(col, pack.d);
        let entry = output * pack.d + input;
        let nonlinear = if k == 0 {
            self.fourier(n, entry)
        } else {
            let difference = if n >= k {
                self.fourier(n - k, entry)
            } else {
                self.fourier(k - n, entry).conj()
            };
            let sum = self.fourier(n + k, entry);
            if imaginary_col {
                let value = difference - sum;
                C64::new(-value.im, value.re)
            } else {
                difference + sum
            }
        };
        let mut linear = if n == k {
            self.ports[n * self.square + entry]
        } else {
            C64::ZERO
        };
        if imaginary_col {
            linear = C64::new(-linear.im, linear.re);
        }
        let derivative = linear - nonlinear;
        if imaginary_row {
            derivative.im
        } else {
            derivative.re
        }
    }
}

fn coordinate(index: usize, d: usize) -> (usize, usize, bool) {
    if index < d {
        (0, index, false)
    } else {
        let block = (index - d) / d;
        (block / 2 + 1, index % d, block % 2 == 1)
    }
}

#[cfg(test)]
mod tests;
