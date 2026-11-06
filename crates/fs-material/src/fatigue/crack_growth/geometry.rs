//! Declared crack geometry with analytical growth through log-log Y(a) pieces.
//!
//! The original table is retained without sorting, discarding or fitting knots.
//! Inside one interval Y=Y0*(a/a0)^s, so the Paris rate is proportional to
//! a^[m*(s+1/2)]. The same analytical power-law integral handles increasing,
//! constant and decreasing intensity. Each interval's extrema are its endpoints;
//! checking every crossed knot also catches interior peaks of the whole table.
//! Interpolation/model error is not a statistical confidence bound. The caller
//! must supply a resolved geometry solution with the declared stress/length
//! conventions; a table alone does not establish small-scale yielding or change
//! the uncracked mechanical stiffness.

use std::sync::Arc;

use fs_math::det;

use super::{
    ConstantGeometryCrack, CrackGrowthError, ParisCrackGrowth, StressCycle, positive,
    positive_log_ratio, power_cycles, power_cycles_to_boundary, power_length,
};

mod centre_crack;
pub use centre_crack::FiniteWidthCentreCrack;

/// One original geometry-factor sample. Length and Y must be strictly positive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrackGeometryPoint {
    pub crack_length_m: f64,
    pub geometry_factor: f64,
}

/// Applicability of an original Y(a) solution, separate from the material fit.
/// These are caller declarations, not independently validated geometry evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct CrackGeometryTableContext {
    pub source_reference: String,
    /// For example: half-length of a central through crack, in metres.
    pub crack_length_convention: String,
    /// For example: remote gross-section tensile stress in an uncracked plate.
    pub nominal_stress_basis: String,
}

/// Immutable bounded log-log interpolation of the original Y(a) points.
/// No extrapolation or monotonicity repair is performed. The retained point
/// count bounds per-block work; no interval is traversed more than once.
#[derive(Debug, Clone, PartialEq)]
pub struct TabulatedCrackGeometry {
    context: CrackGeometryTableContext,
    points: Vec<CrackGeometryPoint>,
    log_slopes: Vec<f64>,
    centre_crack: Option<FiniteWidthCentreCrack>,
}

impl TabulatedCrackGeometry {
    pub fn new(
        context: CrackGeometryTableContext,
        points: Vec<CrackGeometryPoint>,
        maximum_points: usize,
    ) -> Result<Self, CrackGrowthError> {
        Self::new_with_cancel(context, points, maximum_points, &mut || false)
    }

    /// Retain the table only after every point and interval has been admitted.
    pub fn new_with_cancel(
        context: CrackGeometryTableContext,
        points: Vec<CrackGeometryPoint>,
        maximum_points: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CrackGrowthError> {
        poll(cancelled)?;
        if points.len() > maximum_points {
            return Err(CrackGrowthError::Budget);
        }
        if points.len() < 2
            || context.source_reference.trim().is_empty()
            || context.crack_length_convention.trim().is_empty()
            || context.nominal_stress_basis.trim().is_empty()
        {
            return Err(CrackGrowthError::Invalid(
                "declare geometry source, crack-length and nominal-stress conventions with at least two points",
            ));
        }
        let mut log_slopes = Vec::with_capacity(points.len() - 1);
        for (index, point) in points.iter().enumerate() {
            poll(cancelled)?;
            if !positive(point.crack_length_m) || !positive(point.geometry_factor) {
                return Err(CrackGrowthError::Invalid(
                    "geometry points need finite positive crack length and factor",
                ));
            }
            if index > 0 {
                let before = points[index - 1];
                if point.crack_length_m <= before.crack_length_m {
                    return Err(CrackGrowthError::Invalid(
                        "original geometry lengths must be strictly increasing",
                    ));
                }
                let slope = positive_log_ratio(before.geometry_factor, point.geometry_factor)
                    / positive_log_ratio(before.crack_length_m, point.crack_length_m);
                if !slope.is_finite() {
                    return Err(CrackGrowthError::Arithmetic);
                }
                log_slopes.push(slope);
            }
        }
        poll(cancelled)?;
        Ok(Self {
            context,
            points,
            log_slopes,
            centre_crack: None,
        })
    }

    pub fn context(&self) -> &CrackGeometryTableContext {
        &self.context
    }
    pub fn points(&self) -> &[CrackGeometryPoint] {
        &self.points
    }
    /// Original physical geometry, retained only by the typed centre-crack
    /// constructor. A caller's table description cannot establish plate width.
    pub fn finite_width_centre_crack(&self) -> Option<FiniteWidthCentreCrack> {
        self.centre_crack
    }
    pub fn length_bounds(&self) -> (f64, f64) {
        (
            self.points[0].crack_length_m,
            self.points[self.points.len() - 1].crack_length_m,
        )
    }

    /// Original values at knots; log-log interpolation strictly inside support.
    pub fn factor_at(&self, length_m: f64) -> Result<f64, CrackGrowthError> {
        check_length(length_m, self.length_bounds())?;
        let index = self.segment_index(length_m);
        if index == self.points.len() - 1 {
            return Ok(self.points[index].geometry_factor);
        }
        self.factor_in_segment(index, length_m)
    }

    fn segment_index(&self, length_m: f64) -> usize {
        self.points
            .partition_point(|point| point.crack_length_m <= length_m)
            .saturating_sub(1)
    }

    fn factor_in_segment(&self, index: usize, length_m: f64) -> Result<f64, CrackGrowthError> {
        let left = self.points[index];
        let right = self.points[index + 1];
        if length_m == left.crack_length_m || self.log_slopes[index] == 0.0 {
            return Ok(left.geometry_factor);
        }
        if length_m == right.crack_length_m {
            return Ok(right.geometry_factor);
        }
        let log_change = self.log_slopes[index] * positive_log_ratio(left.crack_length_m, length_m);
        let factor = if log_change.abs() < 0.5 {
            left.geometry_factor + left.geometry_factor * det::expm1(log_change)
        } else {
            det::exp(det::ln(left.geometry_factor) + log_change)
        };
        if !positive(factor) {
            return Err(CrackGrowthError::Arithmetic);
        }
        Ok(factor)
    }
}

/// One shared geometry owner for material, quasistatic and inertial consumers.
/// Tables are immutable and shared by checkpoint clones.
#[derive(Debug, Clone, PartialEq)]
pub enum CrackGeometry {
    Constant(ConstantGeometryCrack),
    Tabulated(Arc<TabulatedCrackGeometry>),
}

impl From<ConstantGeometryCrack> for CrackGeometry {
    fn from(value: ConstantGeometryCrack) -> Self {
        Self::Constant(value)
    }
}
impl From<TabulatedCrackGeometry> for CrackGeometry {
    fn from(value: TabulatedCrackGeometry) -> Self {
        Self::Tabulated(Arc::new(value))
    }
}
impl From<Arc<TabulatedCrackGeometry>> for CrackGeometry {
    fn from(value: Arc<TabulatedCrackGeometry>) -> Self {
        Self::Tabulated(value)
    }
}

impl CrackGeometry {
    pub fn length_bounds(&self) -> (f64, f64) {
        match self {
            Self::Constant(value) => (value.minimum_length_m, value.maximum_length_m),
            Self::Tabulated(table) => table.length_bounds(),
        }
    }

    pub fn factor_at(&self, length_m: f64) -> Result<f64, CrackGrowthError> {
        self.validate()?;
        self.check_length(length_m)?;
        match self {
            Self::Constant(value) => Ok(value.geometry_factor),
            Self::Tabulated(table) => table.factor_at(length_m),
        }
    }

    pub(super) fn validate(&self) -> Result<(), CrackGrowthError> {
        match self {
            Self::Constant(value) => value.validate(),
            Self::Tabulated(_) => Ok(()), // Construction already checked every immutable point.
        }
    }
    pub(super) fn check_length(&self, length_m: f64) -> Result<(), CrackGrowthError> {
        match self {
            Self::Constant(value) => value.check_length(length_m),
            Self::Tabulated(table) => check_length(length_m, table.length_bounds()),
        }
    }
    pub(super) fn intensity(&self, length_m: f64, stress_pa: f64) -> Result<f64, CrackGrowthError> {
        match self {
            Self::Constant(value) => value.intensity(length_m, stress_pa),
            Self::Tabulated(table) => intensity(length_m, stress_pa, table.factor_at(length_m)?),
        }
    }
}

fn check_length(length_m: f64, (minimum, maximum): (f64, f64)) -> Result<(), CrackGrowthError> {
    if !length_m.is_finite() || length_m < minimum || length_m > maximum {
        return Err(CrackGrowthError::CrackLengthOutsideDomain {
            length_m,
            minimum,
            maximum,
        });
    }
    Ok(())
}

fn intensity(length_m: f64, stress_pa: f64, factor: f64) -> Result<f64, CrackGrowthError> {
    let value = (stress_pa * det::sqrt(length_m)) * (factor * det::sqrt(std::f64::consts::PI));
    if !value.is_finite() || (stress_pa > 0.0 && value <= 0.0) {
        return Err(CrackGrowthError::Arithmetic);
    }
    Ok(value)
}

impl ParisCrackGrowth {
    /// Evaluate one interval endpoint, retaining the source's exact knot value.
    /// Intensity is monotone within a log-log geometry interval, even when it
    /// decreases, so admission at both ends covers that entire interval.
    fn table_endpoint(
        &self,
        table: &TabulatedCrackGeometry,
        index: usize,
        length_m: f64,
        range_pa: f64,
        maximum_pa: f64,
    ) -> Result<f64, CrackGrowthError> {
        let factor = table.factor_in_segment(index, length_m)?;
        let delta_k = intensity(length_m, range_pa, factor)?;
        self.law.check_delta_k(delta_k)?;
        self.law
            .check_maximum_k(intensity(length_m, maximum_pa, factor)?)?;
        Ok(delta_k)
    }

    pub(super) fn tabulated_length(
        &self,
        table: &TabulatedCrackGeometry,
        cycle: StressCycle,
        mut rate: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, CrackGrowthError> {
        let mut a = self.state.crack_length_m;
        let mut spent = 0.0;
        let mut compensation = 0.0;
        let mut index = table.segment_index(a);
        let range = cycle.maximum_pa - cycle.minimum_pa;
        let ratio = cycle.minimum_pa / cycle.maximum_pa;
        while index < table.log_slopes.len() {
            poll(cancelled)?;
            let end = table.points[index + 1].crack_length_m;
            let power = self.law.parameters.exponent * (0.5 + table.log_slopes[index]);
            let to_boundary = power_cycles_to_boundary(a, end, rate, power)?;
            // Use the same compensated prefix as the inverse. Repeatedly
            // subtracting interval counts can falsely leave exposure beyond a
            // final knot even when the input is that knot's exact forecast.
            let (boundary_total, boundary_compensation) =
                exposure_prefix(spent, compensation, to_boundary)?;
            if cycle.cycles < boundary_total {
                let remaining = cycle.cycles - spent;
                let next = power_length(a, rate, power, remaining)?;
                if !positive(next) || next <= a || next > end {
                    return Err(CrackGrowthError::Arithmetic);
                }
                self.table_endpoint(table, index, next, range, cycle.maximum_pa)?;
                poll(cancelled)?;
                return Ok(next);
            }
            // A boundary reached within this block must itself be representable;
            // never round an unresolved finite-cycle singularity into a knot.
            power_cycles(a, end, rate, power)?;
            let delta_k = self.table_endpoint(table, index, end, range, cycle.maximum_pa)?;
            poll(cancelled)?;
            if cycle.cycles == boundary_total {
                return Ok(end);
            }
            a = end;
            spent = boundary_total;
            compensation = boundary_compensation;
            index += 1;
            if index < table.log_slopes.len() {
                rate = self.law.growth_rate(delta_k, ratio)?;
            }
        }
        Err(CrackGrowthError::GeometrySupportExhausted {
            length_m: a,
            remaining_cycles: cycle.cycles - spent,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn tabulated_cycles(
        &self,
        table: &TabulatedCrackGeometry,
        minimum_pa: f64,
        maximum_pa: f64,
        target: f64,
        mut rate: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, CrackGrowthError> {
        let mut a = self.state.crack_length_m;
        let mut index = table.segment_index(a);
        let range = maximum_pa - minimum_pa;
        let ratio = minimum_pa / maximum_pa;
        let mut total = 0.0;
        let mut compensation = 0.0;
        while index < table.log_slopes.len() {
            poll(cancelled)?;
            let end = target.min(table.points[index + 1].crack_length_m);
            let power = self.law.parameters.exponent * (0.5 + table.log_slopes[index]);
            let cycles = power_cycles(a, end, rate, power)?;
            let delta_k = self.table_endpoint(table, index, end, range, maximum_pa)?;
            // Compensated summation keeps small interval exposures when one
            // geometry table spans substantially different local growth rates.
            (total, compensation) = exposure_prefix(total, compensation, cycles)?;
            if !positive(total) {
                return Err(CrackGrowthError::Arithmetic);
            }
            poll(cancelled)?;
            if end == target {
                return Ok(total);
            }
            a = end;
            index += 1;
            rate = self.law.growth_rate(delta_k, ratio)?;
        }
        Err(CrackGrowthError::Arithmetic) // Caller admitted target inside table support.
    }
}

/// One deterministic compensated prefix shared by both directions. A distinct
/// geometry endpoint needs a distinct representable cycle coordinate. Infinity
/// is useful only to advancement: a finite requested block cannot reach it.
fn exposure_prefix(
    total: f64,
    compensation: f64,
    interval: f64,
) -> Result<(f64, f64), CrackGrowthError> {
    let adjusted = interval - compensation;
    let next = total + adjusted;
    if next.is_nan() || next <= total {
        return Err(CrackGrowthError::Arithmetic);
    }
    if next.is_infinite() {
        return Ok((next, 0.0));
    }
    Ok((next, (next - total) - adjusted))
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CrackGrowthError> {
    if cancelled() {
        Err(CrackGrowthError::Cancelled)
    } else {
        Ok(())
    }
}
