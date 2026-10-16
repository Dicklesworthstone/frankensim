//! Measured geometry authoring for the extensible assembly-v1 Parametric table.
//! The versioned vocabulary below retains values, absolute coordinate bounds,
//! row identities and methods in ordinary dimensioned parameters. It changes
//! neither the assembly wire layout nor scenario IR. Rows are citations, not
//! resolved authority. No meshing, material defaults or specimen inference.
use super::{DimensionBound, GeometrySource, ParameterKey, ReferenceState};
use crate::entity::ToleranceSource;
use fs_ivl::interval::Interval;
use fs_qty::{Dims, QtyAny};

mod arching;
mod bore;
mod outline;
pub use arching::ArchingSurface;
pub use bore::{BoreInterpolant, BoreProfile, BoreSampling, GeneratedSegment};
pub use outline::OutlineWithHoles;

pub const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
const VOLUME: Dims = Dims([3, 0, 0, 0, 0, 0]);
/// Per-dimension source row in the M8-W1 ledger; resolved above L3.
pub type LedgerRowId = String;

#[derive(Debug, Clone, PartialEq)]
pub enum GeneratorError {
    MissingSource { dimension: String },
    InvalidInput(String),
    FixtureRequired,
    NonMonotoneStations,
    SelfIntersection { loop_index: usize },
    HoleOutsideOutline { hole: usize },
    OverlappingHoles,
    BudgetExceeded { resource: &'static str, cap: usize },
    Cancelled,
}
impl std::fmt::Display for GeneratorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "geometry generator: {self:?}")
    }
}
impl std::error::Error for GeneratorError {}
fn bad(s: impl Into<String>) -> GeneratorError {
    GeneratorError::InvalidInput(s.into())
}
fn cap(n: usize, limit: usize, resource: &'static str) -> Result<(), GeneratorError> {
    if n > limit {
        Err(GeneratorError::BudgetExceeded {
            resource,
            cap: limit,
        })
    } else {
        Ok(())
    }
}
fn fixture(label: &str) -> bool {
    label
        .strip_prefix("test-fixture:")
        .is_some_and(|s| !s.trim().is_empty())
}
fn label_valid(label: &str) -> Result<(), GeneratorError> {
    if label.trim().is_empty() || label.len() > 256 || label.chars().any(char::is_control) {
        return Err(bad("nonempty bounded geometry label required"));
    }
    if label.starts_with("test-fixture:") && !fixture(label) {
        return Err(bad("fixture label needs a nonempty suffix"));
    }
    Ok(())
}

/// A nominal SI observation and its *absolute* lower/upper measurement bounds.
/// A missing ledger row is admissible only for explicitly labelled fixtures.
#[derive(Debug, Clone, PartialEq)]
pub struct MeasuredDimension {
    pub value: QtyAny,
    pub bound: DimensionBound,
    pub ledger_row: Option<LedgerRowId>,
}
impl MeasuredDimension {
    fn validate(&self, label: &str, dims: Dims) -> Result<(), GeneratorError> {
        let b = &self.bound;
        if self.value.dims != dims
            || b.property.dims != dims
            || b.lower.dims != dims
            || b.upper.dims != dims
            || !self.value.value.is_finite()
            || !b.lower.value.is_finite()
            || !b.upper.value.is_finite()
            || b.lower.value > self.value.value
            || self.value.value > b.upper.value
            || b.property.name.is_empty()
            || b.property.name.len() > 256
        {
            return Err(bad(format!("invalid dimension/bound {}", b.property.name)));
        }
        match (&self.ledger_row, &b.source) {
            (Some(row), ToleranceSource::Measured { ledger_row, method })
            | (
                Some(row),
                ToleranceSource::Reconstructed {
                    ledger_row,
                    route: method,
                },
            ) if !row.trim().is_empty()
                && row == ledger_row
                && !method.trim().is_empty()
                && row.len() <= 256
                && method.len() <= 256
                && !row.chars().any(char::is_control)
                && !method.chars().any(char::is_control) =>
            {
                if row.starts_with("test-fixture:") && !fixture(label) {
                    return Err(GeneratorError::FixtureRequired);
                }
            }
            (None, ToleranceSource::Assumed { rationale })
                if fixture(label) && !rationale.trim().is_empty() && rationale.len() <= 256 => {}
            _ => {
                return Err(GeneratorError::MissingSource {
                    dimension: b.property.name.clone(),
                });
            }
        }
        Ok(())
    }
    fn interval(&self) -> Interval {
        Interval::new(self.bound.lower.value, self.bound.upper.value)
    }
}

/// A derived dimension carries all its original observations and bounds.
/// Intervals enclose the declared geometry, not unmeasured model discrepancy.
#[derive(Debug, Clone, PartialEq)]
pub struct DerivedDimension {
    pub value: QtyAny,
    pub lower: QtyAny,
    pub upper: QtyAny,
    pub inputs: Vec<MeasuredDimension>,
}
impl DerivedDimension {
    fn new(
        value: f64,
        enclosure: Interval,
        dims: Dims,
        inputs: Vec<MeasuredDimension>,
    ) -> Result<Self, GeneratorError> {
        if !value.is_finite()
            || !enclosure.lo().is_finite()
            || !enclosure.hi().is_finite()
            || !enclosure.contains(value)
        {
            return Err(bad(
                "derived quantity is unrepresentable or outside its bound",
            ));
        }
        Ok(Self {
            value: QtyAny::new(value, dims),
            lower: QtyAny::new(enclosure.lo(), dims),
            upper: QtyAny::new(enclosure.hi(), dims),
            inputs,
        })
    }
}

/// Explicit caps apply before source traversal and output allocation.
#[derive(Debug, Clone, Copy)]
pub struct GeneratorBudget {
    pub max_dimensions: usize,
    pub max_samples: usize,
    pub max_work: usize,
}
/// Refinement inserts samples; it never removes a measured station.
#[derive(Debug, Clone, Copy)]
pub struct SamplingRung {
    pub subdivisions_per_span: usize,
}

#[derive(Debug)]
pub enum Generator {
    BoreProfile(BoreProfile),
    OutlineWithHoles(OutlineWithHoles),
    ArchingSurface(ArchingSurface),
}
impl Generator {
    /// Reconstruct and re-admit a versioned table, retaining every source bit.
    pub fn from_source(
        source: &GeometrySource,
        budget: GeneratorBudget,
    ) -> Result<Self, GeneratorError> {
        let GeometrySource::Parametric { generator, params } = source else {
            return Err(bad("expected Parametric geometry"));
        };
        cap(
            params.len(),
            budget.max_dimensions.saturating_mul(5).saturating_add(16),
            "parameter records",
        )?;
        let mut r = TableReader::new(params);
        let label = r.label()?;
        let result = match generator.as_str() {
            "BoreProfile.v1" => Self::BoreProfile(BoreProfile::read(label, &mut r, budget)?),
            "OutlineWithHoles.v1" => {
                Self::OutlineWithHoles(OutlineWithHoles::read(label, &mut r, budget)?)
            }
            "ArchingSurface.v1" => {
                Self::ArchingSurface(ArchingSurface::read(label, &mut r, budget)?)
            }
            _ => return Err(bad("unknown measured generator version")),
        };
        r.finish()?;
        Ok(result)
    }
    /// Fixture labels survive serialization and must refuse specimen claims.
    pub fn require_specimen(&self) -> Result<(), GeneratorError> {
        let label = match self {
            Self::BoreProfile(v) => &v.label,
            Self::OutlineWithHoles(v) => &v.label,
            Self::ArchingSurface(v) => &v.label,
        };
        if fixture(label) {
            Err(GeneratorError::FixtureRequired)
        } else {
            Ok(())
        }
    }
}
pub(super) fn validate_source(source: &GeometrySource, limit: usize) -> Result<(), GeneratorError> {
    if matches!(source,GeometrySource::Parametric{generator,..} if matches!(generator.as_str(),"BoreProfile.v1"|"OutlineWithHoles.v1"|"ArchingSurface.v1"))
    {
        Generator::from_source(
            source,
            GeneratorBudget {
                max_dimensions: limit,
                max_samples: limit,
                max_work: limit.saturating_mul(limit),
            },
        )?;
    }
    Ok(())
}

// Metadata strings are losslessly hex-encoded in legal ParameterKey tokens;
// dimensioned numeric values remain ordinary QtyAny records in assembly v1.
fn hex(s: &str) -> String {
    s.bytes().map(|b| format!("{b:02x}")).collect()
}
fn unhex(s: &str) -> Result<String, GeneratorError> {
    if s.len() % 2 != 0 || s.len() > 512 || !s.is_ascii() {
        return Err(bad("invalid metadata token"));
    }
    let bytes = (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).map_err(|_| bad("invalid metadata hex")))
        .collect::<Result<Vec<_>, _>>()?;
    String::from_utf8(bytes).map_err(|_| bad("invalid metadata UTF8"))
}
struct TableWriter {
    params: Vec<(ParameterKey, QtyAny)>,
    next: usize,
}
impl TableWriter {
    fn new(label: &str) -> Self {
        Self {
            params: vec![(
                ParameterKey::new(format!("label.{}", hex(label)), Dims::NONE),
                QtyAny::dimensionless(0.0),
            )],
            next: 0,
        }
    }
    fn count(&mut self, name: &str, n: usize) {
        self.params.push((
            ParameterKey::new(name, Dims::NONE),
            QtyAny::dimensionless(n as f64),
        ));
    }
    fn dimension(&mut self, d: &MeasuredDimension) {
        let p = format!("d{}", self.next);
        self.next += 1;
        for (suffix, q) in [
            ("value", d.value),
            ("lower", d.bound.lower),
            ("upper", d.bound.upper),
        ] {
            self.params
                .push((ParameterKey::new(format!("{p}.{suffix}"), q.dims), q));
        }
        let (kind, method) = match &d.bound.source {
            ToleranceSource::Measured { method, .. } => ("m", method),
            ToleranceSource::Reconstructed { route, .. } => ("r", route),
            ToleranceSource::Assumed { rationale } => ("a", rationale),
            _ => unreachable!("admitted dimension"),
        };
        // Assembly's quantity codec canonicalizes signed zero. Preserve those
        // three source sign bits in metadata without changing its frozen wire.
        let zeros = [d.value, d.bound.lower, d.bound.upper]
            .iter()
            .enumerate()
            .fold(0u8, |bits, (i, q)| {
                bits | (u8::from(q.value.to_bits() == (-0.0f64).to_bits()) << i)
            });
        self.params.push((
            ParameterKey::new(
                format!(
                    "{p}.meta.{}.{}.{}.{}.{zeros}",
                    hex(&d.bound.property.name),
                    kind,
                    hex(d.ledger_row.as_deref().unwrap_or("")),
                    hex(method)
                ),
                Dims::NONE,
            ),
            QtyAny::dimensionless(0.0),
        ));
    }
    fn source(self, id: &str) -> GeometrySource {
        GeometrySource::Parametric {
            generator: id.into(),
            params: self.params,
        }
    }
}
struct TableReader<'a> {
    params: &'a [(ParameterKey, QtyAny)],
    at: usize,
    next: usize,
}
impl<'a> TableReader<'a> {
    fn new(params: &'a [(ParameterKey, QtyAny)]) -> Self {
        Self {
            params,
            at: 0,
            next: 0,
        }
    }
    fn take(&mut self) -> Result<(&'a ParameterKey, QtyAny), GeneratorError> {
        let (k, v) = self
            .params
            .get(self.at)
            .ok_or_else(|| bad("incomplete generator table"))?;
        self.at += 1;
        if k.dims != v.dims || !v.value.is_finite() {
            return Err(bad("invalid parameter dimensions/value"));
        }
        Ok((k, *v))
    }
    fn label(&mut self) -> Result<String, GeneratorError> {
        let (k, v) = self.take()?;
        if v != QtyAny::dimensionless(0.0) {
            return Err(bad("invalid label marker"));
        }
        let label = unhex(
            k.name
                .strip_prefix("label.")
                .ok_or_else(|| bad("missing label"))?,
        )?;
        label_valid(&label)?;
        Ok(label)
    }
    fn count(&mut self, name: &str) -> Result<usize, GeneratorError> {
        let (k, v) = self.take()?;
        if k.name != name
            || v.dims != Dims::NONE
            || v.value < 0.0
            || v.value.fract() != 0.0
            || v.value > u32::MAX as f64
        {
            return Err(bad("invalid generator count"));
        }
        Ok(v.value as usize)
    }
    fn dimension(&mut self) -> Result<MeasuredDimension, GeneratorError> {
        let p = format!("d{}", self.next);
        self.next += 1;
        let mut q = Vec::new();
        for suffix in ["value", "lower", "upper"] {
            let (k, v) = self.take()?;
            if k.name != format!("{p}.{suffix}") {
                return Err(bad("unordered dimension table"));
            }
            q.push(v);
        }
        let (k, v) = self.take()?;
        if v != QtyAny::dimensionless(0.0) {
            return Err(bad("invalid metadata marker"));
        }
        let m = k
            .name
            .strip_prefix(&format!("{p}.meta."))
            .ok_or_else(|| bad("missing dimension metadata"))?
            .split('.')
            .collect::<Vec<_>>();
        if m.len() != 5 {
            return Err(bad("invalid dimension metadata"));
        }
        let property = unhex(m[0])?;
        let row = unhex(m[2])?;
        let method = unhex(m[3])?;
        let zeros = m[4]
            .parse::<u8>()
            .map_err(|_| bad("invalid signed-zero metadata"))?;
        if zeros > 7 {
            return Err(bad("invalid signed-zero mask"));
        }
        for (i, value) in q.iter_mut().enumerate() {
            if zeros & (1 << i) != 0 {
                if value.value != 0.0 {
                    return Err(bad("signed-zero mask on nonzero quantity"));
                }
                value.value = -0.0;
            }
        }
        let source = match m[1] {
            "m" => ToleranceSource::Measured {
                ledger_row: row.clone(),
                method,
            },
            "r" => ToleranceSource::Reconstructed {
                ledger_row: row.clone(),
                route: method,
            },
            "a" => ToleranceSource::Assumed { rationale: method },
            _ => return Err(bad("invalid source kind")),
        };
        Ok(MeasuredDimension {
            value: q[0],
            bound: DimensionBound {
                property: ParameterKey::new(property, q[0].dims),
                lower: q[1],
                upper: q[2],
                source,
            },
            ledger_row: if row.is_empty() { None } else { Some(row) },
        })
    }
    fn finish(&self) -> Result<(), GeneratorError> {
        if self.at == self.params.len() {
            Ok(())
        } else {
            Err(bad("trailing generator parameters"))
        }
    }
}

fn linear_knots(axis: &[f64]) -> Result<fs_rep_nurbs::KnotVector<f64>, GeneratorError> {
    let mut knots = vec![axis[0]];
    knots.extend_from_slice(axis);
    knots.push(axis[axis.len() - 1]);
    fs_rep_nurbs::KnotVector::new(knots, 1).map_err(|e| bad(e.to_string()))
}
fn checked_axis(axis: &[MeasuredDimension], label: &str) -> Result<(), GeneratorError> {
    if axis.len() < 2 {
        return Err(bad("axis needs at least two stations"));
    }
    for d in axis {
        d.validate(label, LENGTH)?;
    }
    if axis
        .windows(2)
        .any(|w| w[0].bound.upper.value >= w[1].bound.lower.value)
    {
        return Err(GeneratorError::NonMonotoneStations);
    }
    Ok(())
}
