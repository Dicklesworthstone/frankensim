//! Analytic original-source impedance jets for coupled resonance calculations.
//!
//! Shared poles, original port order, both conjugate members and mutual terms
//! remain in the source chart. Directions are authored physical coefficient
//! derivatives, not derivatives of a fit or a pole-matching procedure.

use super::MimoModel;
use crate::model::PoleTerm;
use fs_math::{c64::C64, dd::Dd};

/// Invalid source/direction data, a true pole, or an unrepresentable jet.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MimoDerivativeError(
    /// Static source, direction or arithmetic refusal reason.
    pub &'static str,
);

impl core::fmt::Display for MimoDerivativeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for MimoDerivativeError {}

/// Local derivatives per unit of one real parameter in the original chart.
///
/// None means exactly zero. Pole slots follow the stored groups; a real pole
/// has a real derivative, while a pair supplies only its stored nonreal
/// member and implies the conjugate derivative. Either pair orientation is
/// accepted and retained. Residue matrices follow the
/// same groups and row-major port order as the original source. Real-pole
/// residues and their derivatives must be real. Signed tangents need not be
/// passive, stable, or positive. No finite parameter interval is admitted.
#[derive(Clone, Copy, Debug, Default)]
pub struct MimoParameterDirection<'a> {
    /// One stored-pole derivative per original group.
    pub poles: Option<&'a [C64]>,
    /// One row-major residue derivative matrix per original group.
    pub residues: Option<&'a [Vec<C64>]>,
    /// Row-major derivative of the direct matrix D.
    pub direct: Option<&'a [f64]>,
    /// Row-major derivative of the improper matrix E.
    pub improper: Option<&'a [f64]>,
}

/// Original impedance and its first two analytic Laplace derivatives.
///
/// Entries are row-major. On the structural/acoustic exp(-i omega t) axis,
/// evaluate at s=-i*omega, then use Z_omega=-i*Z_s and Z_omegaomega=-Z_ss.
/// These derivatives retain dispersion; they are not a damping approximation.
#[derive(Clone, Debug, PartialEq)]
pub struct MimoLaplaceJet {
    /// Original Z(s), in row-major port order.
    pub impedance: Vec<C64>,
    /// First partial derivative dZ/ds at fixed coefficients.
    pub laplace_derivative: Vec<C64>,
    /// Second partial derivative d²Z/ds² at fixed coefficients.
    pub laplace_second_derivative: Vec<C64>,
}

/// Complete local Laplace/parameter jet of the same original impedance.
///
/// The frequency fields hold the parameter fixed. The parameter fields hold
/// s fixed. On s=-i*omega, Z_omega_theta=-i*Z_s_theta.
#[derive(Clone, Debug, PartialEq)]
pub struct MimoLaplaceParameterJet {
    /// The original value and first two Laplace derivatives.
    pub frequency: MimoLaplaceJet,
    /// Partial derivative dZ/dtheta at fixed s.
    pub parameter_derivative: Vec<C64>,
    /// Mixed partial derivative d²Z/(ds dtheta).
    pub laplace_parameter_derivative: Vec<C64>,
}

impl MimoModel {
    /// Evaluate Z, dZ/ds and d²Z/ds² without fitting or forming another model.
    ///
    /// Binary-scaled complex quotient chains use the shared fs-math double-double
    /// arithmetic without explicit squared/cubed pole distances or a rounded
    /// reciprocal. Weak components survive until the final output rounding. All
    /// conjugate terms,
    /// mutual coefficients and the improper E term remain. The deterministic
    /// operation costs O(ports² * stored_groups) with three output matrices and
    /// constant scratch. Admission precedes allocation; the caller owns its
    /// physical frequency band and any enclosing work/memory/cancellation limits.
    ///
    /// This pointwise analytic query admits neither passivity nor causality.
    /// It does not alter runtime history or provide a between-frequency bound.
    pub fn eval_jet(&self, s: C64) -> Result<MimoLaplaceJet, MimoDerivativeError> {
        self.validate_jet(s, None)?;
        let n = self.ports * self.ports;
        let mut result = MimoLaplaceJet {
            impedance: Vec::with_capacity(n),
            laplace_derivative: Vec::with_capacity(n),
            laplace_second_derivative: Vec::with_capacity(n),
        };
        for entry in 0..n {
            let mut sums = [Sum::default(); 3];
            sums[0].add(C64::from_re(self.d[entry]))?;
            sums[0].add(Sum::from(s).multiply(C64::from_re(self.e[entry])))?;
            sums[1].add(C64::from_re(self.e[entry]))?;
            for (group, residues) in self.poles.iter().zip(&self.residues) {
                let (pole, paired) = pole(*group);
                let residue = residues[entry];
                add_frequency(&mut sums, residue, distance(s, pole))?;
                if paired {
                    add_frequency(&mut sums, residue.conj(), distance(s, pole.conj()))?;
                }
            }
            result.impedance.push(sums[0].value()?);
            result.laplace_derivative.push(sums[1].value()?);
            result.laplace_second_derivative.push(sums[2].value()?);
        }
        Ok(result)
    }

    /// Evaluate the full analytic jet, including authored coefficient motion.
    ///
    /// For each member R/(s-p), the parameter derivative is
    /// R_theta/(s-p)+R*p_theta/(s-p)² and its Laplace derivative is
    /// -R_theta/(s-p)²-2*R*p_theta/(s-p)³. Direct and improper tangents are
    /// D_theta+s*E_theta and E_theta. Both members of every pair participate;
    /// no frequency-dependent impedance is frozen while differentiating.
    ///
    /// The fixed grouping, port chart and parameter meaning belong to the
    /// caller. Every constitutive or source-identification chain rule must be
    /// supplied; a refit, mode crossing or topology change is not inferred.
    /// Five output matrices and constant scratch are used, with the same
    /// asymptotic work as eval_jet. Refusal publishes no partial jet.
    pub fn eval_parameter_jet(
        &self,
        s: C64,
        direction: MimoParameterDirection<'_>,
    ) -> Result<MimoLaplaceParameterJet, MimoDerivativeError> {
        self.validate_jet(s, Some(direction))?;
        let n = self.ports * self.ports;
        let mut result = MimoLaplaceParameterJet {
            frequency: MimoLaplaceJet {
                impedance: Vec::with_capacity(n),
                laplace_derivative: Vec::with_capacity(n),
                laplace_second_derivative: Vec::with_capacity(n),
            },
            parameter_derivative: Vec::with_capacity(n),
            laplace_parameter_derivative: Vec::with_capacity(n),
        };
        for entry in 0..n {
            let mut sums = [Sum::default(); 5];
            sums[0].add(C64::from_re(self.d[entry]))?;
            sums[0].add(Sum::from(s).multiply(C64::from_re(self.e[entry])))?;
            sums[1].add(C64::from_re(self.e[entry]))?;
            sums[3].add(C64::from_re(direction.direct.map_or(0.0, |v| v[entry])))?;
            let de = direction.improper.map_or(0.0, |v| v[entry]);
            sums[3].add(Sum::from(s).multiply(C64::from_re(de)))?;
            sums[4].add(C64::from_re(de))?;
            for (group, (term, residues)) in self.poles.iter().zip(&self.residues).enumerate() {
                let (pole, paired) = pole(*term);
                let residue = residues[entry];
                let dp = direction.poles.map_or(C64::ZERO, |v| v[group]);
                let dr = direction.residues.map_or(C64::ZERO, |v| v[group][entry]);
                add_parameter(&mut sums, residue, distance(s, pole), dp, dr)?;
                if paired {
                    add_parameter(
                        &mut sums,
                        residue.conj(),
                        distance(s, pole.conj()),
                        dp.conj(),
                        dr.conj(),
                    )?;
                }
            }
            result.frequency.impedance.push(sums[0].value()?);
            result.frequency.laplace_derivative.push(sums[1].value()?);
            result
                .frequency
                .laplace_second_derivative
                .push(sums[2].value()?);
            result.parameter_derivative.push(sums[3].value()?);
            result.laplace_parameter_derivative.push(sums[4].value()?);
        }
        Ok(result)
    }

    fn validate_jet(
        &self,
        s: C64,
        direction: Option<MimoParameterDirection<'_>>,
    ) -> Result<(), MimoDerivativeError> {
        let entries = self
            .ports
            .checked_mul(self.ports)
            .ok_or(MimoDerivativeError("MIMO jet matrix size overflow"))?;
        let groups = self.poles.len();
        if self.ports == 0
            || !finite(s)
            || self.d.len() != entries
            || self.e.len() != entries
            || self.residues.len() != groups
            || self.d.iter().chain(&self.e).any(|v| !v.is_finite())
            || self
                .residues
                .iter()
                .any(|r| r.len() != entries || r.iter().any(|&z| !finite(z)))
        {
            return Err(MimoDerivativeError(
                "MIMO jet requires finite source data with exact shapes",
            ));
        }
        if let Some(d) = direction {
            if d.poles
                .is_some_and(|v| v.len() != groups || v.iter().any(|&z| !finite(z)))
                || d.residues.is_some_and(|v| {
                    v.len() != groups
                        || v.iter()
                            .any(|r| r.len() != entries || r.iter().any(|&z| !finite(z)))
                })
                || d.direct
                    .is_some_and(|v| v.len() != entries || v.iter().any(|z| !z.is_finite()))
                || d.improper
                    .is_some_and(|v| v.len() != entries || v.iter().any(|z| !z.is_finite()))
            {
                return Err(MimoDerivativeError(
                    "MIMO jet requires finite direction data with exact shapes",
                ));
            }
        }
        for (i, term) in self.poles.iter().enumerate() {
            let (p, paired) = pole(*term);
            if !finite(p)
                || (paired && p.im == 0.0)
                || (!paired && self.residues[i].iter().any(|z| z.im != 0.0))
                || direction.is_some_and(|d| {
                    !paired
                        && (d.poles.is_some_and(|v| v[i].im != 0.0)
                            || d.residues.is_some_and(|v| v[i].iter().any(|z| z.im != 0.0)))
                })
            {
                return Err(MimoDerivativeError(
                    "MIMO jet must preserve real poles and conjugate-pair grouping",
                ));
            }
            if s == p || (paired && s == p.conj()) {
                return Err(MimoDerivativeError(
                    "MIMO jet is undefined at a stored pole",
                ));
            }
        }
        Ok(())
    }
}

fn pole(term: PoleTerm) -> (C64, bool) {
    match term {
        PoleTerm::Real { pole, .. } => (C64::from_re(pole), false),
        PoleTerm::Pair { pole, .. } => (pole, true),
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

// Preserve the authored source coordinates before subtraction rounds. This
// matters in R_theta/d + R*p_theta/d² when the physical direction nearly
// cancels, and when a finite s-p lies outside the standalone f64 range.
fn distance(s: C64, p: C64) -> Sum {
    Sum {
        re: Number::from(s.re).add(Number::from(p.re).negative()),
        im: Number::from(s.im).add(Number::from(p.im).negative()),
    }
}

fn frequency(residue: C64, distance: Sum) -> [Sum; 3] {
    let value = Sum::from(residue).divide(distance);
    let first = value.divide(distance).negative();
    let second = first.divide(distance).negative().twice();
    [value, first, second]
}

fn add_frequency(
    sums: &mut [Sum; 3],
    residue: C64,
    distance: Sum,
) -> Result<(), MimoDerivativeError> {
    for (sum, value) in sums.iter_mut().zip(frequency(residue, distance)) {
        sum.add(value)?;
    }
    Ok(())
}

fn add_parameter(
    sums: &mut [Sum; 5],
    residue: C64,
    distance: Sum,
    pole_derivative: C64,
    residue_derivative: C64,
) -> Result<(), MimoDerivativeError> {
    let values = frequency(residue, distance);
    for (sum, &value) in sums[..3].iter_mut().zip(&values) {
        sum.add(value)?;
    }
    let residue_term = Sum::from(residue_derivative).divide(distance);
    sums[3].add(residue_term)?;
    sums[3].add(values[1].multiply(pole_derivative).negative())?;
    sums[4].add(residue_term.divide(distance).negative())?;
    sums[4].add(values[2].multiply(pole_derivative).negative())?;
    Ok(())
}

/// Each component retains its own exponent. Scaling an entire complex number
/// together could erase a weak quadrature before a later division or parameter
/// coefficient restores it. All significand arithmetic belongs to fs-math::Dd.
#[derive(Clone, Copy, Default)]
struct Number {
    fraction: Dd,
    exponent: i32,
}

impl Number {
    fn normalized(fraction: Dd, exponent: i32) -> Self {
        if fraction.hi == 0.0 {
            return Self::default();
        }
        let (hi, shift) = binary_parts(fraction.hi);
        Self {
            fraction: Dd::from_pair(hi, scale_binary(fraction.lo, -shift)),
            exponent: exponent + shift,
        }
    }

    fn add(self, other: Self) -> Self {
        if self.fraction.hi == 0.0 {
            return other;
        }
        if other.fraction.hi == 0.0 {
            return self;
        }
        let exponent = self.exponent.max(other.exponent);
        let rescale = |value: Self| {
            Dd::from_pair(
                scale_binary(value.fraction.hi, value.exponent - exponent),
                scale_binary(value.fraction.lo, value.exponent - exponent),
            )
        };
        Self::normalized(rescale(self) + rescale(other), exponent)
    }

    fn negative(mut self) -> Self {
        self.fraction = -self.fraction;
        self
    }

    fn multiply(self, other: Self) -> Self {
        Self::normalized(
            self.fraction * other.fraction,
            self.exponent + other.exponent,
        )
    }

    fn divide(self, other: Self) -> Self {
        Self::normalized(
            self.fraction / other.fraction,
            self.exponent - other.exponent,
        )
    }

    fn value(self) -> f64 {
        scale_binary(self.fraction.to_f64(), self.exponent)
    }
}

impl From<f64> for Number {
    fn from(value: f64) -> Self {
        let (fraction, exponent) = binary_parts(value);
        Self {
            fraction: Dd::from_f64(fraction),
            exponent,
        }
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    re: Number,
    im: Number,
}

impl From<C64> for Sum {
    fn from(z: C64) -> Self {
        Self {
            re: z.re.into(),
            im: z.im.into(),
        }
    }
}

impl Sum {
    fn add(&mut self, z: impl Into<Self>) -> Result<(), MimoDerivativeError> {
        let z = z.into();
        self.re = self.re.add(z.re);
        self.im = self.im.add(z.im);
        Ok(())
    }

    fn negative(self) -> Self {
        Self {
            re: self.re.negative(),
            im: self.im.negative(),
        }
    }

    fn twice(mut self) -> Self {
        self.re.exponent += 1;
        self.im.exponent += 1;
        self
    }

    fn multiply(self, z: C64) -> Self {
        let z = Self::from(z);
        Self {
            re: self
                .re
                .multiply(z.re)
                .add(self.im.multiply(z.im).negative()),
            im: self.re.multiply(z.im).add(self.im.multiply(z.re)),
        }
    }

    fn divide(self, z: impl Into<Self>) -> Self {
        let z = z.into();
        let denominator = z.re.multiply(z.re).add(z.im.multiply(z.im));
        Self {
            re: self
                .re
                .multiply(z.re)
                .add(self.im.multiply(z.im))
                .divide(denominator),
            im: self
                .im
                .multiply(z.re)
                .add(self.re.multiply(z.im).negative())
                .divide(denominator),
        }
    }

    fn value(self) -> Result<C64, MimoDerivativeError> {
        let result = C64::new(self.re.value(), self.im.value());
        if finite(result) {
            Ok(result)
        } else {
            Err(MimoDerivativeError("MIMO jet is not representable"))
        }
    }
}

// Exact decomposition, including subnormals. Bounded derivative chains use at
// most a handful of products, so their separated exponents cannot overflow i32.
fn binary_parts(value: f64) -> (f64, i32) {
    if value == 0.0 {
        return (value, 0);
    }
    let mut bits = value.to_bits();
    let mut exponent = ((bits >> 52) & 0x7ff) as i32;
    let correction = if exponent == 0 {
        bits = (value * f64::from_bits(1077_u64 << 52)).to_bits();
        exponent = ((bits >> 52) & 0x7ff) as i32;
        -54
    } else {
        0
    };
    (
        f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52)),
        exponent - 1023 + correction,
    )
}

fn scale_binary(value: f64, exponent: i32) -> f64 {
    if value == 0.0 {
        return value;
    }
    let (fraction, own) = binary_parts(value);
    let exponent = exponent + own;
    if exponent > 1023 {
        f64::INFINITY.copysign(value)
    } else if exponent < -1075 {
        0.0_f64.copysign(value)
    } else if exponent >= -1022 {
        fraction * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent == -1075 {
        (fraction * 0.5) * f64::from_bits(1)
    } else {
        fraction * f64::from_bits(1_u64 << (exponent + 1074))
    }
}
