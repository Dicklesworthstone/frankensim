//! Authored physical mechanical directions in the existing harmonic state chart.
//!
//! The complete state is `[q,p=M*v,z=sqrt(k)*viscous_strain]`. Differentiate
//! the original coefficient laws, not a mass factorization, modal fit or solve:
//! `Q_pp,theta=-M^-1 M_theta M^-1`, `J_pp,theta=-skew(C_theta)`,
//! `R_pp,theta=sym(C_theta)` and `G_p,theta=P_theta`. A Maxwell coupling
//! `c=sqrt(k)*b` contributes `c_theta*c^T+c*c_theta^T` to Q_qq,
//! `-c_theta` to Q_qz and `-(tau_theta/tau)/tau` to R_zz.

use super::{LinearMechanicalSystem, MAX_MECHANICAL_STATES};
use crate::harmonic::{HarmonicError, HarmonicParameterDirection};
use crate::interconnection::accumulation::Sum;
use fs_math::det;

/// One signed local direction of the original dimensional source matrices.
///
/// None holds that entire law fixed. M and relaxed K directions are exactly
/// symmetric n-by-n arrays; C is a full n-by-n array, including gyroscopic
/// terms, and P is n-by-ports. Included branch arrays have one entry per
/// ORIGINAL arm, including equal-time and zero-projection arms. Derivatives
/// need not be positive: the base source is already admitted, and no finite
/// parameter step or material-validity interval is inferred.
#[derive(Clone, Copy, Debug, Default)]
pub struct MechanicalHarmonicDirection<'a> {
    pub mass: Option<&'a [f64]>,
    pub damping: Option<&'a [f64]>,
    pub relaxed_stiffness: Option<&'a [f64]>,
    pub force_map: Option<&'a [f64]>,
    pub branches: Option<&'a [MechanicalRelaxationDirection<'a>]>,
}

/// Derivatives per unit of the SAME real parameter as the matrix direction.
#[derive(Clone, Copy, Debug, Default)]
pub struct MechanicalRelaxationDirection<'a> {
    /// n original displacement entries of db/dtheta; None holds b fixed.
    pub projection: Option<&'a [f64]>,
    /// dk/dtheta of the excess (not relaxed) stiffness of this arm.
    pub stiffness: f64,
    /// d tau/dtheta in seconds per parameter unit, not a decay-rate slope.
    pub relaxation_time_s: f64,
}

/// Explicit cold-preparation limits. No direction or arm is dropped to fit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MechanicalHarmonicDirectionControls {
    /// Complete `[q,p,z]` dimension, including every material memory.
    pub max_state_dimension: usize,
    /// Conservative peak internal numerical payload, excluding the source.
    pub max_dense_bytes: usize,
    /// Conservative normalized scalar product/contraction work score.
    pub max_multiply_terms: u64,
}

/// Immutable coefficient direction, prepared once and reused across frequencies.
///
/// `as_harmonic_direction` borrows the owned complete J/R/Q/G arrays. The
/// original state/port order and number of arms remain fixed; do not apply it
/// to another mechanical source. No inverse, square root or matrix product is
/// repeated by that accessor. Physical stress slopes can use the stored
/// `d sqrt(k)/dtheta` alongside native memory-effort slopes.
pub struct PreparedMechanicalHarmonicDirection {
    coordinates: usize,
    n: usize,
    m: usize,
    interconnection: Vec<f64>,
    dissipation: Vec<f64>,
    storage: Vec<f64>,
    ports: Vec<f64>,
    branch_roots: Vec<f64>,
    branch_rates: Vec<f64>,
}

impl PreparedMechanicalHarmonicDirection {
    pub const fn coordinates(&self) -> usize {
        self.coordinates
    }
    pub const fn state_dim(&self) -> usize {
        self.n
    }
    pub const fn port_dim(&self) -> usize {
        self.m
    }

    pub fn as_harmonic_direction(&self) -> HarmonicParameterDirection<'_> {
        HarmonicParameterDirection {
            interconnection: Some(&self.interconnection),
            dissipation: Some(&self.dissipation),
            storage: Some(&self.storage),
            ports: Some(&self.ports),
            resisting_flow: None,
        }
    }

    /// d sqrt(k_j)/dtheta in original arm order; zero for fixed stiffness.
    pub fn branch_root_derivatives(&self) -> &[f64] {
        &self.branch_roots
    }
    /// d(1/tau_j)/dtheta in original arm order; zero for fixed time constants.
    pub fn branch_rate_derivatives(&self) -> &[f64] {
        &self.branch_rates
    }
}

impl LinearMechanicalSystem {
    /// Prepare all physical coefficient chain rules through the cached mass inverse.
    ///
    /// Two compensated dense products form `-M^-1 M_theta M^-1`; both words
    /// of the first product survive the second. Binary numerical scales are
    /// retained until a complete coefficient is rounded, avoiding overflowing
    /// or vanishing intermediate products when coordinate units change. The
    /// original complete physical arrays are never normalized or replaced.
    ///
    /// With N complete states, n displacements, m ports and b arms, the payload
    /// ceiling is `8*(3*N*N+N*m+2*b)+24*(n*n+2*n)` bytes, and the work score is
    /// `16*n^3+32*b*n^2+16*N^2+8*N*m+32*b`. Limits and source layouts are
    /// checked before allocation. Cancellation is polled per dense row and arm
    /// and before publication. No partial result or physical state is published.
    pub fn prepare_harmonic_parameter_direction(
        &self,
        direction: MechanicalHarmonicDirection<'_>,
        controls: MechanicalHarmonicDirectionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalHarmonicDirection, HarmonicError> {
        poll(&mut cancelled)?;
        let n = self.coordinates;
        let total = self.system.state_dim();
        let ports = self.system.port_dim();
        let arms = self.branches.len();
        admit_capacity(n, total, ports, arms, controls)?;
        validate(direction, n, ports, arms, &mut cancelled)?;
        poll(&mut cancelled)?;
        let mut result = PreparedMechanicalHarmonicDirection {
            coordinates: n,
            n: total,
            m: ports,
            interconnection: vec![0.0; total * total],
            dissipation: vec![0.0; total * total],
            storage: vec![0.0; total * total],
            ports: vec![0.0; total * ports],
            branch_roots: vec![0.0; arms],
            branch_rates: vec![0.0; arms],
        };
        if let Some(dm) = direction.mass {
            let mut first = vec![ScaledSum::default(); n * n];
            for i in 0..n {
                poll(&mut cancelled)?;
                for j in 0..n {
                    for k in 0..n {
                        first[i * n + j].product(self.inverse_mass[i * n + k], dm[k * n + j])?;
                    }
                }
            }
            for i in 0..n {
                poll(&mut cancelled)?;
                for j in 0..=i {
                    let mut entry = ScaledSum::default();
                    for k in 0..n {
                        entry.add_scaled(first[i * n + k], -self.inverse_mass[k * n + j])?;
                    }
                    let value = entry.value()?;
                    // Symmetry is an algebraic identity of this admitted dM;
                    // install both entries from the SAME completed contraction.
                    result.storage[(n + i) * total + n + j] = value;
                    result.storage[(n + j) * total + n + i] = value;
                }
            }
        }
        for i in 0..n {
            poll(&mut cancelled)?;
            if let Some(dc) = direction.damping {
                result.dissipation[(n + i) * total + n + i] = dc[i * n + i];
                for j in 0..i {
                    let a = dc[i * n + j];
                    let b = dc[j * n + i];
                    let symmetric = a.midpoint(b);
                    let skew = b.midpoint(-a);
                    if (symmetric == 0.0 && a != -b) || (skew == 0.0 && a != b) {
                        return Err(arithmetic());
                    }
                    result.dissipation[(n + i) * total + n + j] = symmetric;
                    result.dissipation[(n + j) * total + n + i] = symmetric;
                    result.interconnection[(n + i) * total + n + j] = skew;
                    result.interconnection[(n + j) * total + n + i] = -skew;
                }
            }
            if let Some(dp) = direction.force_map {
                result.ports[(n + i) * ports..(n + i + 1) * ports]
                    .copy_from_slice(&dp[i * ports..(i + 1) * ports]);
            }
        }
        // Keep every arm contribution unrounded until it has been combined
        // with the relaxed stiffness. Opposing material directions can cancel.
        let mut elastic = vec![ScaledSum::default(); n * n];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..=i {
                elastic[i * n + j].product(at(direction.relaxed_stiffness, i * n + j), 1.0)?;
            }
        }
        let mut coupling = vec![0.0; n];
        let mut coupling_derivative = vec![ScaledSum::default(); n];
        for (arm, branch) in self.branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let d = direction
                .branches
                .map_or(MechanicalRelaxationDirection::default(), |a| a[arm]);
            let root = det::sqrt(branch.stiffness);
            let root_derivative = quotient(d.stiffness, root, 2.0)?;
            let rate_derivative = quotient(
                -d.relaxation_time_s,
                branch.relaxation_time_s,
                branch.relaxation_time_s,
            )?;
            result.branch_roots[arm] = root_derivative;
            result.branch_rates[arm] = rate_derivative;
            let z = 2 * n + arm;
            result.dissipation[z * total + z] = rate_derivative;
            for i in 0..n {
                // Same represented coupling as the admitted source storage.
                coupling[i] = root * branch.projection[i];
                let mut change = ScaledSum::default();
                change.product(root_derivative, branch.projection[i])?;
                change.product(root, at(d.projection, i))?;
                let value = -change.value()?;
                result.storage[i * total + z] = value;
                result.storage[z * total + i] = value;
                coupling_derivative[i] = change;
            }
            for i in 0..n {
                poll(&mut cancelled)?;
                for j in 0..=i {
                    elastic[i * n + j].add_scaled(coupling_derivative[i], coupling[j])?;
                    elastic[i * n + j].add_scaled(coupling_derivative[j], coupling[i])?;
                }
            }
        }
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..=i {
                let value = elastic[i * n + j].value()?;
                result.storage[i * total + j] = value;
                result.storage[j * total + i] = value;
            }
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

fn validate(
    d: MechanicalHarmonicDirection<'_>,
    n: usize,
    ports: usize,
    arms: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), HarmonicError> {
    for matrix in [d.mass, d.damping, d.relaxed_stiffness]
        .into_iter()
        .flatten()
    {
        poll(cancelled)?;
        if matrix.len() != n * n || matrix.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete finite original n*n mechanical direction required",
            ));
        }
    }
    for matrix in [d.mass, d.relaxed_stiffness].into_iter().flatten() {
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..i {
                if matrix[i * n + j] != matrix[j * n + i] {
                    return Err(invalid(
                        "original mass and relaxed stiffness directions must be symmetric",
                    ));
                }
            }
        }
    }
    if let Some(p) = d.force_map {
        if p.len() != n * ports || p.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete finite n*ports mechanical force-map direction required",
            ));
        }
    }
    if let Some(branches) = d.branches {
        if branches.len() != arms {
            return Err(invalid(
                "one physical direction per original Maxwell arm required",
            ));
        }
        for branch in branches {
            poll(cancelled)?;
            if !branch.stiffness.is_finite()
                || !branch.relaxation_time_s.is_finite()
                || branch
                    .projection
                    .is_some_and(|p| p.len() != n || p.iter().any(|v| !v.is_finite()))
            {
                return Err(invalid(
                    "finite original Maxwell stiffness/time/projection directions required",
                ));
            }
        }
    }
    Ok(())
}

fn admit_capacity(
    n: usize,
    total: usize,
    ports: usize,
    arms: usize,
    controls: MechanicalHarmonicDirectionControls,
) -> Result<(), HarmonicError> {
    // Source constructor already bounds each dimension; checked arithmetic
    // still makes admission independent of target pointer width.
    let square = total.checked_mul(total).ok_or_else(capacity)?;
    let port_extent = total.checked_mul(ports).ok_or_else(capacity)?;
    let scratch = n
        .checked_mul(n)
        .and_then(|v| v.checked_add(2 * n))
        .and_then(|v| v.checked_mul(core::mem::size_of::<ScaledSum>()))
        .ok_or_else(capacity)?;
    let bytes = square
        .checked_mul(3)
        .and_then(|v| v.checked_add(port_extent))
        .and_then(|v| v.checked_add(2 * arms))
        .and_then(|v| v.checked_mul(8))
        .and_then(|v| v.checked_add(scratch))
        .ok_or_else(capacity)?;
    let n64 = n as u64;
    let work = n64
        .checked_pow(3)
        .and_then(|v| v.checked_mul(16))
        .and_then(|v| {
            (arms as u64)
                .checked_mul(32 * n64 * n64)
                .and_then(|b| v.checked_add(b))
        })
        .and_then(|v| v.checked_add(16 * square as u64))
        .and_then(|v| v.checked_add(8 * port_extent as u64))
        .and_then(|v| v.checked_add(32 * arms as u64))
        .ok_or_else(capacity)?;
    if total > MAX_MECHANICAL_STATES
        || total > controls.max_state_dimension
        || bytes > controls.max_dense_bytes
        || work > controls.max_multiply_terms
    {
        return Err(capacity());
    }
    Ok(())
}

fn at(values: Option<&[f64]>, index: usize) -> f64 {
    values.map_or(0.0, |a| a[index])
}
fn invalid(what: &'static str) -> HarmonicError {
    HarmonicError::Invalid { what }
}
fn arithmetic() -> HarmonicError {
    invalid("representable complete physical mechanical coefficient direction")
}
fn capacity() -> HarmonicError {
    invalid("complete mechanical direction dimension/payload/work capacity")
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HarmonicError> {
    if cancelled() {
        Err(HarmonicError::Cancelled)
    } else {
        Ok(())
    }
}

// Local binary scaling around the native two-word Sum. This is preparation
// arithmetic, not another solver or precision owner. No reciprocal of a mass
// or time constant is squared before its authored numerator can restore range.
#[derive(Clone, Copy, Default)]
struct ScaledSum {
    sum: Sum,
    exponent: i32,
}

impl ScaledSum {
    fn product(&mut self, a: f64, b: f64) -> Result<(), HarmonicError> {
        self.product_at(a, b, 0)
    }
    fn product_at(&mut self, a: f64, b: f64, exponent: i32) -> Result<(), HarmonicError> {
        if a == 0.0 || b == 0.0 {
            return Ok(());
        }
        let (a, ae) = binary_parts(a);
        let (b, be) = binary_parts(b);
        let mut term = Sum::default();
        term.product(a, b).ok_or_else(arithmetic)?;
        self.add(term, exponent + ae + be)
    }
    fn add_scaled(&mut self, source: Self, coefficient: f64) -> Result<(), HarmonicError> {
        let (high, low) = source.sum.components().ok_or_else(arithmetic)?;
        self.product_at(high, coefficient, source.exponent)?;
        self.product_at(low, coefficient, source.exponent)
    }
    fn add(&mut self, term: Sum, exponent: i32) -> Result<(), HarmonicError> {
        let (a, b) = self.sum.components().ok_or_else(arithmetic)?;
        let (c, d) = term.components().ok_or_else(arithmetic)?;
        if a == 0.0 && b == 0.0 {
            self.sum = term;
            self.exponent = exponent;
            return Ok(());
        }
        let common = self.exponent.max(exponent);
        let mut sum = Sum::default();
        for (value, shift) in [
            (a, self.exponent - common),
            (b, self.exponent - common),
            (c, exponent - common),
            (d, exponent - common),
        ] {
            let scaled = scale_binary(value, shift);
            if value != 0.0 && scaled == 0.0 {
                return Err(arithmetic());
            }
            sum.product(scaled, 1.0).ok_or_else(arithmetic)?;
        }
        self.sum = sum;
        self.exponent = common;
        Ok(())
    }
    fn value(self) -> Result<f64, HarmonicError> {
        let sum = self.sum.value().ok_or_else(arithmetic)?;
        let value = scale_binary(sum, self.exponent);
        if !value.is_finite() || (value == 0.0 && sum != 0.0) {
            Err(arithmetic())
        } else {
            Ok(value)
        }
    }
}

fn quotient(numerator: f64, first: f64, second: f64) -> Result<f64, HarmonicError> {
    if numerator == 0.0 {
        return Ok(numerator);
    }
    let (a, ae) = binary_parts(numerator);
    let (b, be) = binary_parts(first);
    let (c, ce) = binary_parts(second);
    let result = scale_binary((a / b) / c, ae - be - ce);
    if !result.is_finite() || result == 0.0 {
        Err(arithmetic())
    } else {
        Ok(result)
    }
}

// Finite nonzero x = mantissa*2^exponent, 1 <= |mantissa| < 2.
fn binary_parts(mut value: f64) -> (f64, i32) {
    let mut bits = value.to_bits();
    let mut exponent = ((bits >> 52) & 0x7ff) as i32;
    let correction = if exponent == 0 {
        value *= f64::from_bits((1023_u64 + 54) << 52);
        bits = value.to_bits();
        exponent = ((bits >> 52) & 0x7ff) as i32;
        -54
    } else {
        0
    };
    (
        f64::from_bits((bits & (1_u64 << 63 | ((1_u64 << 52) - 1))) | (1023_u64 << 52)),
        exponent - 1023 + correction,
    )
}

fn scale_binary(value: f64, exponent: i32) -> f64 {
    if value == 0.0 {
        return value;
    }
    let (mantissa, own_exponent) = binary_parts(value);
    let exponent = exponent + own_exponent;
    if exponent > 1023 {
        f64::INFINITY.copysign(value)
    } else if exponent < -1075 {
        0.0_f64.copysign(value)
    } else if exponent >= -1022 {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent == -1075 {
        (0.5 * mantissa) * f64::from_bits(1)
    } else {
        mantissa * f64::from_bits(1_u64 << (exponent + 1074))
    }
}
