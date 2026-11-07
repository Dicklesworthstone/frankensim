//! Original sparse pressure equations and their frozen-coefficient tangent.
//!
//! For `A = K - omega^2 M - omega B + i omega G` and `A p = i omega b`,
//! `A p' = i b + i omega b' + 2 omega M p + B p - i G p`.
//! Derivatives are per rad/s. Material, geometry and the complete complex
//! admittance `G + i B` are held fixed. Both equations use the retained LU.

use fs_la::eigen_complex::{LuComplex, LuComplexError};
use fs_math::{c64::C64, dd::Dd};
use fs_sparse::Csr;

#[derive(Debug)]
pub(super) enum PressureError {
    Invalid(&'static str),
    Factor(LuComplexError),
    Residual { observed: f64, limit: f64 },
    Cancelled,
}

impl From<LuComplexError> for PressureError {
    fn from(error: LuComplexError) -> Self {
        match error {
            LuComplexError::Cancelled => Self::Cancelled,
            _ => Self::Factor(error),
        }
    }
}

#[derive(Debug)]
pub(super) struct PressureSolution {
    pub(super) pressure: Vec<C64>,
    pub(super) backward_error: f64,
}

/// Borrowed from one prepared field. Its factor is for `S A_ff S`, with
/// `scales = diag(S)` in the same ascending `free_dofs` order. The original
/// preparation owner establishes this binding; it is not a public raw API.
pub(super) struct PressurePencil<'a> {
    pub(super) omega: f64,
    pub(super) stiffness: &'a Csr,
    pub(super) mass: &'a Csr,
    pub(super) conductance: &'a Csr,
    pub(super) susceptance: &'a Csr,
    pub(super) factor: &'a LuComplex,
    pub(super) scales: &'a [f64],
    pub(super) free_dofs: &'a [usize],
    pub(super) residual_tolerance: f64,
}

// Retain unrounded tangent forcing for original-equation residual correction.
// The magnitude sums its separate physical terms, including coherent sources.
struct Forcing {
    values: Vec<[Dd; 2]>,
    magnitude: Vec<f64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PressureError> {
    if cancelled() {
        Err(PressureError::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

impl PressurePencil<'_> {
    fn admit(
        &self,
        load: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PressureError> {
        poll(cancelled)?;
        let n = self.stiffness.nrows();
        if n == 0
            || load.len() != n
            || [
                self.stiffness,
                self.mass,
                self.conductance,
                self.susceptance,
            ]
            .iter()
            .any(|matrix| matrix.nrows() != n || matrix.ncols() != n)
            || self.free_dofs.is_empty()
            || self.free_dofs.len() != self.scales.len()
            || self.free_dofs.iter().any(|&row| row >= n)
            || self.free_dofs.windows(2).any(|pair| pair[0] >= pair[1])
            || self.scales.iter().any(|s| !s.is_finite() || *s <= 0.0)
            || !self.omega.is_finite()
            || self.omega <= 0.0
            || !(self.omega * self.omega).is_finite()
            || self.omega * self.omega == 0.0
            || !self.residual_tolerance.is_finite()
            || self.residual_tolerance <= 0.0
            || load.iter().any(|&value| !finite(value))
        {
            return Err(PressureError::Invalid(
                "invalid prepared pressure pencil or source",
            ));
        }
        Ok(())
    }

    pub(super) fn solve(
        &self,
        load: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        self.admit(load, cancelled)?;
        let mut forcing = Forcing {
            values: Vec::with_capacity(load.len()),
            magnitude: Vec::with_capacity(load.len()),
        };
        for &value in load {
            poll(cancelled)?;
            // Preserve the primary source's original f64 arithmetic.
            let rhs = C64::new(-self.omega * value.im, self.omega * value.re);
            if !finite(rhs) {
                return Err(PressureError::Invalid("harmonic velocity source overflow"));
            }
            forcing
                .values
                .push([Dd::from_f64(rhs.re), Dd::from_f64(rhs.im)]);
            forcing.magnitude.push(rhs.abs());
        }
        self.solve_forcing(&forcing, cancelled)
    }

    /// `pressure` is the admitted primary solution of this same pencil.
    pub(super) fn frequency_derivative(
        &self,
        pressure: &[C64],
        load: &[C64],
        load_derivative: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        self.admit(load, cancelled)?;
        if pressure.len() != load.len()
            || load_derivative.len() != load.len()
            || pressure
                .iter()
                .chain(load_derivative)
                .any(|&value| !finite(value))
        {
            return Err(PressureError::Invalid(
                "invalid harmonic pressure or source tangent",
            ));
        }
        let omega = Dd::from_f64(self.omega);
        let mut forcing = Forcing {
            values: vec![[Dd::from_f64(0.0); 2]; load.len()],
            magnitude: vec![0.0; load.len()],
        };
        // Eliminated rows are reactions, not equations; exact zero pressure
        // and zero tangent are lifted into those coordinates after solving.
        for &row in self.free_dofs {
            poll(cancelled)?;
            let mut real =
                -Dd::from_f64(load[row].im) - omega * Dd::from_f64(load_derivative[row].im);
            let mut imaginary =
                Dd::from_f64(load[row].re) + omega * Dd::from_f64(load_derivative[row].re);
            let mut magnitude = load[row].abs() + self.omega * load_derivative[row].abs();
            for (matrix, multiplier, rotate) in [
                (self.mass, Dd::from_f64(2.0) * omega, false),
                (self.susceptance, Dd::from_f64(1.0), false),
                (self.conductance, Dd::from_f64(-1.0), true),
            ] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    let coefficient = multiplier * Dd::from_f64(value);
                    let p = pressure[column];
                    if rotate {
                        real = real - coefficient * Dd::from_f64(p.im);
                        imaginary = imaginary + coefficient * Dd::from_f64(p.re);
                    } else {
                        real = real + coefficient * Dd::from_f64(p.re);
                        imaginary = imaginary + coefficient * Dd::from_f64(p.im);
                    }
                    magnitude += coefficient.to_f64().abs() * p.abs();
                }
            }
            if !finite(C64::new(real.to_f64(), imaginary.to_f64())) || !magnitude.is_finite() {
                return Err(PressureError::Invalid(
                    "harmonic frequency tangent forcing overflow",
                ));
            }
            forcing.values[row] = [real, imaginary];
            forcing.magnitude[row] = magnitude;
        }
        self.solve_forcing(&forcing, cancelled)
    }

    fn solve_forcing(
        &self,
        forcing: &Forcing,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        let mut free_pressure = self
            .free_dofs
            .iter()
            .zip(self.scales)
            .map(|(&row, &scale)| {
                C64::new(
                    forcing.values[row][0].to_f64(),
                    forcing.values[row][1].to_f64(),
                )
                .scale(scale)
            })
            .collect::<Vec<_>>();
        if free_pressure.iter().any(|&value| !finite(value)) {
            return Err(PressureError::Invalid(
                "unrepresentable equilibrated pressure source",
            ));
        }
        self.factor
            .solve_with_cancel(&mut free_pressure, &mut *cancelled)?;
        let mut pressure = vec![C64::ZERO; self.stiffness.nrows()];
        for ((&row, &p), &scale) in self.free_dofs.iter().zip(&free_pressure).zip(self.scales) {
            pressure[row] = p.scale(scale);
        }
        for correction in 0..=2 {
            poll(cancelled)?;
            if pressure.iter().any(|&z| !finite(z)) {
                return Err(PressureError::Invalid(
                    "nonfinite harmonic pressure or tangent",
                ));
            }
            let (mut residual, backward_error) = self.residual(&pressure, forcing, cancelled)?;
            if correction == 2 || backward_error < 1e-18 {
                if backward_error > self.residual_tolerance {
                    return Err(PressureError::Residual {
                        observed: backward_error,
                        limit: self.residual_tolerance,
                    });
                }
                poll(cancelled)?;
                return Ok(PressureSolution {
                    pressure,
                    backward_error,
                });
            }
            for (r, &scale) in residual.iter_mut().zip(self.scales) {
                *r = r.scale(scale);
            }
            self.factor
                .solve_with_cancel(&mut residual, &mut *cancelled)?;
            for (i, &row) in self.free_dofs.iter().enumerate() {
                pressure[row] = pressure[row] + residual[i].scale(self.scales[i]);
            }
        }
        unreachable!("the bounded correction loop returns on its last residual")
    }

    fn residual(
        &self,
        pressure: &[C64],
        forcing: &Forcing,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, f64), PressureError> {
        let omega = Dd::from_f64(self.omega);
        let omega2 = omega * omega;
        let mut residual = vec![C64::ZERO; self.free_dofs.len()];
        let mut worst: f64 = 0.0;
        for (free_row, &row) in self.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let [mut real, mut imaginary] = forcing.values[row];
            let mut scale = forcing.magnitude[row];
            for (matrix, multiplier, rotate) in [
                (self.stiffness, Dd::from_f64(1.0), false),
                (self.mass, -omega2, false),
                (self.susceptance, -omega, false),
                (self.conductance, omega, true),
            ] {
                let (cols, vals) = matrix.row(row);
                for (&col, &value) in cols.iter().zip(vals) {
                    let coefficient = multiplier * Dd::from_f64(value);
                    let p = pressure[col];
                    if rotate {
                        real = real + coefficient * Dd::from_f64(p.im);
                        imaginary = imaginary - coefficient * Dd::from_f64(p.re);
                    } else {
                        real = real - coefficient * Dd::from_f64(p.re);
                        imaginary = imaginary - coefficient * Dd::from_f64(p.im);
                    }
                    scale += coefficient.to_f64().abs() * p.abs();
                }
            }
            let value = C64::new(real.to_f64(), imaginary.to_f64());
            if !finite(value) || !scale.is_finite() {
                return Err(PressureError::Invalid(
                    "unrepresentable physical pressure residual",
                ));
            }
            worst = worst.max(value.abs() / scale.max(f64::MIN_POSITIVE));
            residual[free_row] = value;
        }
        Ok((residual, worst))
    }
}

#[cfg(test)]
mod tests;
