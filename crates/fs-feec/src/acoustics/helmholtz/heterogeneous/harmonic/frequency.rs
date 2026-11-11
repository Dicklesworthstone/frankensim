//! Original sparse pressure equations and their analytic frequency tangent.
//!
//! For `A = K - omega^2 M - omega B + i omega G` and `A p = i omega b`,
//! `A p' = i b + i omega b' + 2 omega M p + (B + omega B') p
//!          - i (G + omega G') p`.
//! Derivatives are per rad/s. Bulk material and geometry remain fixed;
//! boundary derivatives are explicit inputs. Both equations use the retained LU.

use crate::acoustics::backward_error::row_scale;
use fs_la::eigen_complex::{LuComplex, LuComplexError};
use fs_math::{c64::C64, dd::Dd};
use fs_sparse::Csr;

pub(crate) mod coupled;

#[derive(Debug)]
pub(super) enum PressureError {
    Invalid(&'static str),
    Factor(LuComplexError),
    Residual {
        observed: f64,
        limit: f64,
    },
    Limit {
        what: &'static str,
        requested: usize,
        maximum: usize,
    },
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

/// Borrowed from one original physical pencil. The coupled wall owner also
/// uses these real polynomial blocks, without assuming symmetric or positive
/// blocks outside the fluid equations. Its factor is for `S A_ff S`, with
/// `scales = diag(S)` in the same ascending `free_dofs` order. The original
/// preparation owner establishes this binding; it is not a public raw API.
pub(super) struct PressurePencil<'a> {
    pub(super) omega: f64,
    pub(super) stiffness: &'a Csr,
    pub(super) mass: &'a Csr,
    pub(super) conductance: &'a Csr,
    pub(super) susceptance: &'a Csr,
    pub(super) boundary_frequency_derivative: Option<(&'a Csr, &'a Csr)>,
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

/// Real part of a complete original complex bilinear form and the magnitude
/// of its separate summands. Peak-phasor work/energy factors belong to callers.
fn real_form(
    matrix: &Csr,
    left: &[C64],
    right: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), PressureError> {
    if left.len() != matrix.nrows() || right.len() != matrix.ncols() {
        return Err(PressureError::Invalid("physical form dimensions disagree"));
    }
    let mut sum = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for (row, p) in left.iter().enumerate() {
        poll(cancelled)?;
        let (cols, values) = matrix.row(row);
        for (&col, &coefficient) in cols.iter().zip(values) {
            let q = right[col];
            sum = sum
                + Dd::from_f64(coefficient)
                    * (Dd::from_f64(p.re) * Dd::from_f64(q.re)
                        + Dd::from_f64(p.im) * Dd::from_f64(q.im));
            absolute += coefficient.abs() * p.abs() * q.abs();
        }
    }
    let value = sum.to_f64();
    if !value.is_finite() || !absolute.is_finite() {
        return Err(PressureError::Invalid(
            "unrepresentable physical energy or work form",
        ));
    }
    Ok((value, absolute))
}

pub(super) fn quadratic_nonnegative(
    matrix: &Csr,
    field: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, PressureError> {
    let (value, absolute) = real_form(matrix, field, field, cancelled)?;
    if value < -64.0 * f64::EPSILON * absolute {
        return Err(PressureError::Invalid("nonpositive physical energy form"));
    }
    Ok(value.max(0.0))
}

/// Harmonic image of the original surface mass, equilibrium spring, Kelvin
/// resistance and positive Maxwell arms. The adapter supplies the coefficients
/// of its existing wall owner, without fitting or changing their units.
/// `Z = r + i(omega sigma - K/omega) + sum k/(1/tau + i omega)`.
/// The reciprocal is admitted only for finite, nonzero impedance; a lossless
/// wall pole requires the uncondensed coupled-state formulation instead.
pub(super) fn mechanical_admittance(
    omega: f64,
    surface_density: f64,
    stiffness_per_area: f64,
    resistance: f64,
    relaxation: impl IntoIterator<Item = (f64, f64)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(C64, C64), PressureError> {
    poll(cancelled)?;
    if ![omega, surface_density, stiffness_per_area]
        .iter()
        .all(|v| v.is_finite() && *v > 0.0)
        || !resistance.is_finite()
        || resistance < 0.0
    {
        return Err(PressureError::Invalid(
            "physical mechanical wall coefficients required",
        ));
    }
    let w = Dd::from_f64(omega);
    let mass = Dd::from_f64(surface_density);
    let spring = Dd::from_f64(stiffness_per_area);
    // Retain spring/mass cancellation at a high-Q wall resonance.
    let mut real = Dd::from_f64(resistance);
    let mut imaginary = w * mass - spring / w;
    let mut real_derivative = Dd::from_f64(0.0);
    let mut imaginary_derivative = mass + spring / w / w;
    for (stiffness, tau) in relaxation {
        poll(cancelled)?;
        if !stiffness.is_finite()
            || stiffness <= 0.0
            || !tau.is_finite()
            || tau <= 0.0
            || !(1.0 / tau).is_finite()
        {
            return Err(PressureError::Invalid(
                "positive representable Maxwell wall arm required",
            ));
        }
        // The rate form avoids forming omega*tau or tau^2. Complex division
        // also avoids squaring a very large/small relaxation denominator.
        let denominator = C64::new(1.0 / tau, omega);
        let arm = C64::from_re(stiffness) / denominator;
        let derivative = C64::new(arm.im, -arm.re) / denominator;
        if !finite(arm) || !finite(derivative) {
            return Err(PressureError::Invalid(
                "mechanical wall relaxation overflow",
            ));
        }
        real = real + Dd::from_f64(arm.re);
        imaginary = imaginary + Dd::from_f64(arm.im);
        real_derivative = real_derivative + Dd::from_f64(derivative.re);
        imaginary_derivative = imaginary_derivative + Dd::from_f64(derivative.im);
    }
    let impedance = C64::new(real.to_f64(), imaginary.to_f64());
    let impedance_derivative = C64::new(real_derivative.to_f64(), imaginary_derivative.to_f64());
    if !finite(impedance) || !finite(impedance_derivative) || impedance == C64::ZERO {
        return Err(PressureError::Invalid(
            "finite nonzero mechanical wall impedance required",
        ));
    }
    let admittance = impedance.recip();
    let derivative = -(impedance_derivative / impedance) / impedance;
    if !finite(admittance) || !finite(derivative) {
        return Err(PressureError::Invalid(
            "unrepresentable mechanical wall admittance or slope",
        ));
    }
    poll(cancelled)?;
    Ok((admittance, derivative))
}

impl PressurePencil<'_> {
    fn admit(
        &self,
        load_len: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PressureError> {
        poll(cancelled)?;
        let n = self.stiffness.nrows();
        if n == 0
            || load_len != n
            || [
                self.stiffness,
                self.mass,
                self.conductance,
                self.susceptance,
            ]
            .iter()
            .any(|matrix| matrix.nrows() != n || matrix.ncols() != n)
            || self.boundary_frequency_derivative.is_some_and(|(g, b)| {
                [g, b]
                    .iter()
                    .any(|matrix| matrix.nrows() != n || matrix.ncols() != n)
            })
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
        self.admit(load.len(), cancelled)?;
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
        self.admit(load.len(), cancelled)?;
        if pressure.len() != load.len()
            || load_derivative.len() != load.len()
            || pressure
                .iter()
                .chain(load)
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
        for &row in self.free_dofs {
            poll(cancelled)?;
            forcing.values[row] = [
                -Dd::from_f64(load[row].im) - omega * Dd::from_f64(load_derivative[row].im),
                Dd::from_f64(load[row].re) + omega * Dd::from_f64(load_derivative[row].re),
            ];
            forcing.magnitude[row] = load[row].abs() + self.omega * load_derivative[row].abs();
        }
        self.frequency_derivative_rhs(pressure, forcing, cancelled)
    }

    /// Actual physical RHS, already including any i*omega factor belonging to
    /// fluid velocity sources. Mechanical momentum forces have no such factor.
    fn solve_rhs(
        &self,
        forcing: &Forcing,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        self.admit_rhs(forcing, cancelled)?;
        self.solve_forcing(forcing, cancelled)
    }

    fn admit_rhs(
        &self,
        forcing: &Forcing,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PressureError> {
        self.admit(forcing.values.len(), cancelled)?;
        if forcing.magnitude.len() != forcing.values.len()
            || forcing
                .values
                .iter()
                .flatten()
                .any(|v| !v.hi.is_finite() || !v.lo.is_finite())
            || forcing.magnitude.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(PressureError::Invalid("invalid complete harmonic RHS"));
        }
        Ok(())
    }

    /// Solve A*x' = rhs' - A'*x, retaining DD physical forcing through the
    /// original-equation residual correction and the same prepared LU.
    fn frequency_derivative_rhs(
        &self,
        pressure: &[C64],
        mut forcing: Forcing,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSolution, PressureError> {
        self.admit_rhs(&forcing, cancelled)?;
        if pressure.len() != forcing.values.len() || pressure.iter().any(|&v| !finite(v)) {
            return Err(PressureError::Invalid(
                "invalid complete harmonic tangent state",
            ));
        }
        let omega = Dd::from_f64(self.omega);
        // Eliminated rows are reactions, not equations; exact zero pressure
        // and zero tangent are lifted into those coordinates after solving.
        for &row in self.free_dofs {
            poll(cancelled)?;
            let [mut real, mut imaginary] = forcing.values[row];
            let mut magnitude = forcing.magnitude[row];
            for (matrix, multiplier, rotate) in [
                (self.mass, Dd::from_f64(2.0) * omega, false),
                (self.susceptance, Dd::from_f64(1.0), false),
                (self.conductance, Dd::from_f64(-1.0), true),
            ]
            .into_iter()
            .chain(
                self.boundary_frequency_derivative
                    .into_iter()
                    .flat_map(|(g, b)| [(b, omega, false), (g, -omega, true)]),
            ) {
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
            let (mut residual, backward_error, componentwise) =
                self.residual(&pressure, forcing, cancelled)?;
            // The refinement stops early only when the strict componentwise error is exact to double-double, as before; the
            // mixed error is smaller and would stop it at the first iterate. The mixed one decides acceptance.
            if correction == 2 || componentwise < 1e-18 {
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
    ) -> Result<(Vec<C64>, f64, f64), PressureError> {
        let omega = Dd::from_f64(self.omega);
        let omega2 = omega * omega;
        let mut residual = vec![C64::ZERO; self.free_dofs.len()];
        let mut worst: f64 = 0.0;
        let mut componentwise: f64 = 0.0;
        let pressure_max = pressure.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        for (free_row, &row) in self.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let [mut real, mut imaginary] = forcing.values[row];
            let mut scale = forcing.magnitude[row];
            let mut absolute_row = 0.0;
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
                    absolute_row += coefficient.to_f64().abs();
                }
            }
            let value = C64::new(real.to_f64(), imaginary.to_f64());
            if !finite(value) || !scale.is_finite() {
                return Err(PressureError::Invalid(
                    "unrepresentable physical pressure residual",
                ));
            }
            let norm = absolute_row * pressure_max + forcing.magnitude[row];
            let denominator = row_scale(scale, norm, self.residual_tolerance);
            worst = worst.max(value.abs() / denominator.max(f64::MIN_POSITIVE));
            componentwise = componentwise.max(value.abs() / scale.max(f64::MIN_POSITIVE));
            residual[free_row] = value;
        }
        Ok((residual, worst, componentwise))
    }
}

#[cfg(test)]
mod tests;
