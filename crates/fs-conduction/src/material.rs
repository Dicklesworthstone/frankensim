//! Conductivity models, and the receipt discipline that makes property
//! provenance travel with the solve.
//!
//! # Why a SAMPLED table and not a live query per element
//!
//! `fs-matdb` answers one query at one point and hands back a
//! [`PropertyUsageReceipt`] describing exactly which claims were
//! considered, which was selected under which policy, and how it was
//! evaluated. Querying per element per Newton iteration would produce a
//! receipt volume nobody can audit and would make the receipt set a
//! function of the iteration path.
//!
//! So a [`ConductivityTable`] is built ONCE from a DECLARED temperature
//! grid: one matdb query per grid point, every receipt retained. The
//! solve then reads the table. That splits the claim cleanly:
//!
//! - the KNOT VALUES are matdb's claim, and each carries its receipt;
//! - the values BETWEEN knots are THIS crate's claim, and the
//!   interpolation is declared: piecewise linear in T.
//!
//! Outside the sampled span the table REFUSES
//! ([`crate::ConductionError::OutsideTemperatureSpan`]). Extrapolating
//! material data is exactly the move that turns a solve into a confident
//! wrong answer, so it is not available.

use fs_matdb::{ClaimSet, PropertyUsageReceipt, QueryPoint, SelectionPolicy};
use fs_qty::Dims;

use crate::ConductionError;

/// SI exponents (m, kg, s, K, A, mol) of a thermal conductivity,
/// W/(m·K) = kg·m·s⁻³·K⁻¹.
pub const CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);

/// The temperature axis name used when querying `fs-matdb`. It matches
/// the `fs_evidence::ValidityDomain` axis convention used by material
/// cards.
pub const TEMPERATURE_AXIS: &str = "T";

/// Where a conductivity number came from. A model is never silently
/// provenance-free: a declared constant SAYS it is declared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProvenanceClass {
    /// Every knot value came from an `fs-matdb` query and its receipt is
    /// retained.
    MatdbReceipts,
    /// The value was supplied inline by the caller. There is no material
    /// provenance and this crate does not invent any.
    Declared,
}

impl ProvenanceClass {
    /// A stable tag for receipts and structured logs.
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self {
            ProvenanceClass::MatdbReceipts => "matdb-receipts",
            ProvenanceClass::Declared => "declared",
        }
    }
}

/// The temperature interval a model is usable over.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TemperatureSpan {
    /// The model has no temperature dependence, so every temperature is
    /// admissible.
    Unbounded,
    /// Sampled over `[low, high]` K; outside it the model refuses.
    Sampled {
        /// Lowest sampled temperature, K.
        low: f64,
        /// Highest sampled temperature, K.
        high: f64,
    },
}

impl TemperatureSpan {
    /// Refuse a temperature outside the span.
    ///
    /// # Errors
    /// [`ConductionError::OutsideTemperatureSpan`].
    pub fn check(self, temperature: f64) -> Result<(), ConductionError> {
        match self {
            TemperatureSpan::Unbounded => Ok(()),
            TemperatureSpan::Sampled { low, high } => {
                if temperature >= low && temperature <= high {
                    Ok(())
                } else {
                    Err(ConductionError::OutsideTemperatureSpan {
                        temperature,
                        low,
                        high,
                    })
                }
            }
        }
    }

    /// Intersection of two spans (the composition law: never wider).
    #[must_use]
    pub fn intersect(self, other: TemperatureSpan) -> TemperatureSpan {
        match (self, other) {
            (TemperatureSpan::Unbounded, s) | (s, TemperatureSpan::Unbounded) => s,
            (
                TemperatureSpan::Sampled { low: a, high: b },
                TemperatureSpan::Sampled { low: c, high: d },
            ) => TemperatureSpan::Sampled {
                low: a.max(c),
                high: b.min(d),
            },
        }
    }
}

/// One scalar conductivity component as a function of temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductivityTable {
    property: String,
    knots: Vec<(f64, f64)>,
    receipts: Vec<PropertyUsageReceipt>,
    provenance: ProvenanceClass,
}

impl ConductivityTable {
    /// A caller-declared constant. Carries NO material provenance and
    /// says so ([`ProvenanceClass::Declared`]).
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] for a non-finite or
    /// non-positive value.
    pub fn declared(value: f64) -> Result<ConductivityTable, ConductionError> {
        if !(value.is_finite() && value > 0.0) {
            return Err(ConductionError::Conductivity {
                what: format!("declared conductivity {value} must be finite and positive"),
            });
        }
        Ok(ConductivityTable {
            property: "declared".to_string(),
            knots: vec![(0.0, value)],
            receipts: Vec::new(),
            provenance: ProvenanceClass::Declared,
        })
    }

    /// A caller-declared `k(T)` curve. Carries NO material provenance
    /// ([`ProvenanceClass::Declared`]) and is usable only inside
    /// `[knots.first().0, knots.last().0]`.
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] for fewer than two knots,
    /// a non-increasing abscissa, or a non-positive conductivity.
    pub fn declared_curve(knots: Vec<(f64, f64)>) -> Result<ConductivityTable, ConductionError> {
        if knots.len() < 2 {
            return Err(ConductionError::Conductivity {
                what: format!(
                    "a declared curve needs at least two knots, got {}",
                    knots.len()
                ),
            });
        }
        for w in knots.windows(2) {
            if !(w[0].0.is_finite() && w[1].0.is_finite() && w[0].0 < w[1].0) {
                return Err(ConductionError::Conductivity {
                    what: "curve abscissae must be finite and strictly increasing".to_string(),
                });
            }
        }
        for &(t, k) in &knots {
            if !(k.is_finite() && k > 0.0) {
                return Err(ConductionError::Conductivity {
                    what: format!("declared conductivity {k} at T = {t} K must be positive"),
                });
            }
        }
        Ok(ConductivityTable {
            property: "declared-curve".to_string(),
            knots,
            receipts: Vec::new(),
            provenance: ProvenanceClass::Declared,
        })
    }

    /// Sample an `fs-matdb` property over a declared temperature grid,
    /// retaining one receipt per grid point.
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] for a malformed grid;
    /// [`ConductionError::MaterialQuery`] wrapping any upstream matdb
    /// refusal (unknown property, out of validity, ambiguous selection);
    /// [`ConductionError::Dimensions`] when the answered sample does not
    /// carry [`CONDUCTIVITY_DIMS`].
    pub fn from_claims(
        claims: &ClaimSet,
        property: &str,
        grid: &[f64],
        policy: SelectionPolicy,
    ) -> Result<ConductivityTable, ConductionError> {
        if grid.len() < 2 {
            return Err(ConductionError::Conductivity {
                what: format!(
                    "a sampled conductivity table needs at least two grid points, got {}",
                    grid.len()
                ),
            });
        }
        for w in grid.windows(2) {
            if !(w[0].is_finite() && w[1].is_finite() && w[0] < w[1]) {
                return Err(ConductionError::Conductivity {
                    what: "the temperature grid must be finite and strictly increasing".to_string(),
                });
            }
        }
        let mut knots = Vec::with_capacity(grid.len());
        let mut receipts = Vec::with_capacity(grid.len());
        for &t in grid {
            let point = QueryPoint::new().with(TEMPERATURE_AXIS, t).map_err(|e| {
                ConductionError::MaterialQuery {
                    property: property.to_string(),
                    temperature: t,
                    upstream: e.to_string(),
                }
            })?;
            let answer = claims.query(property, &point, policy).map_err(|e| {
                ConductionError::MaterialQuery {
                    property: property.to_string(),
                    temperature: t,
                    upstream: e.to_string(),
                }
            })?;
            let sample = &answer.evidence.value;
            if sample.dims != CONDUCTIVITY_DIMS {
                return Err(ConductionError::Dimensions {
                    context: format!("fs-matdb property {property:?} at T = {t} K"),
                    expected: CONDUCTIVITY_DIMS.0,
                    found: sample.dims.0,
                });
            }
            if !(sample.value.is_finite() && sample.value > 0.0) {
                return Err(ConductionError::Conductivity {
                    what: format!(
                        "fs-matdb answered {} W/(m K) for {property:?} at T = {t} K; \
                         conduction requires a finite positive conductivity",
                        sample.value
                    ),
                });
            }
            knots.push((t, sample.value));
            receipts.push(answer.receipt);
        }
        Ok(ConductivityTable {
            property: property.to_string(),
            knots,
            receipts,
            provenance: ProvenanceClass::MatdbReceipts,
        })
    }

    /// The property name this table was sampled from (`"declared"` for a
    /// caller-supplied constant).
    #[must_use]
    pub fn property(&self) -> &str {
        &self.property
    }

    /// The retained `fs-matdb` receipts, one per sampled grid point.
    #[must_use]
    pub fn receipts(&self) -> &[PropertyUsageReceipt] {
        &self.receipts
    }

    /// Where the numbers came from.
    #[must_use]
    pub const fn provenance(&self) -> ProvenanceClass {
        self.provenance
    }

    /// The sampled `(T, k)` knots.
    #[must_use]
    pub fn knots(&self) -> &[(f64, f64)] {
        &self.knots
    }

    /// True when the table actually varies with temperature.
    #[must_use]
    pub fn is_temperature_dependent(&self) -> bool {
        self.knots.len() > 1
    }

    /// The temperature interval this table is usable over.
    #[must_use]
    pub fn span(&self) -> TemperatureSpan {
        if self.knots.len() > 1 {
            TemperatureSpan::Sampled {
                low: self.knots[0].0,
                high: self.knots[self.knots.len() - 1].0,
            }
        } else {
            TemperatureSpan::Unbounded
        }
    }

    fn segment(&self, temperature: f64) -> Result<usize, ConductionError> {
        self.span().check(temperature)?;
        // Deterministic tie-break: a temperature landing exactly on an
        // interior knot uses the segment to its RIGHT, so `eval` and
        // `derivative` agree on which linear piece is in force.
        let last = self.knots.len() - 2;
        let mut lo = 0usize;
        let mut hi = last;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if temperature < self.knots[mid + 1].0 {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        Ok(lo)
    }

    /// `k(T)` by piecewise-linear interpolation between sampled knots.
    ///
    /// # Errors
    /// [`ConductionError::OutsideTemperatureSpan`] outside the sampled
    /// span — this crate never extrapolates material data.
    pub fn eval(&self, temperature: f64) -> Result<f64, ConductionError> {
        if self.knots.len() == 1 {
            return Ok(self.knots[0].1);
        }
        let s = self.segment(temperature)?;
        let (x0, y0) = self.knots[s];
        let (x1, y1) = self.knots[s + 1];
        let w = (temperature - x0) / (x1 - x0);
        Ok(w.mul_add(y1 - y0, y0))
    }

    /// `dk/dT` — the slope of the linear piece in force at `T`.
    ///
    /// # Errors
    /// [`ConductionError::OutsideTemperatureSpan`] outside the span.
    pub fn derivative(&self, temperature: f64) -> Result<f64, ConductionError> {
        if self.knots.len() == 1 {
            return Ok(0.0);
        }
        let s = self.segment(temperature)?;
        let (x0, y0) = self.knots[s];
        let (x1, y1) = self.knots[s + 1];
        Ok((y1 - y0) / (x1 - x0))
    }
}

#[derive(Debug, Clone, PartialEq)]
enum ConductivityKind {
    /// A general constant SPD tensor, supplied inline.
    ConstantTensor(Box<[[f64; 3]; 3]>),
    /// `k(T)·I`.
    Isotropic(Box<ConductivityTable>),
    /// `Σ_i k_i(T)·e_i e_iᵀ` over an orthonormal principal frame.
    Orthotropic {
        axes: Box<[[f64; 3]; 3]>,
        tables: Box<[ConductivityTable; 3]>,
    },
}

/// The conductivity model a solve reads: a symmetric positive-definite
/// 3×3 tensor field that may depend on temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductivityModel {
    kind: ConductivityKind,
}

const ORTHONORMALITY_TOL: f64 = 1e-10;
const SYMMETRY_TOL: f64 = 1e-12;

fn check_spd(k: &[[f64; 3]; 3]) -> Result<(), ConductionError> {
    let mut scale = 0.0f64;
    for row in k {
        for &v in row {
            if !v.is_finite() {
                return Err(ConductionError::Conductivity {
                    what: format!("tensor entry {v} is not finite"),
                });
            }
            scale = scale.max(v.abs());
        }
    }
    if scale == 0.0 {
        return Err(ConductionError::Conductivity {
            what: "the conductivity tensor is identically zero".to_string(),
        });
    }
    for (i, row) in k.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if (value - k[j][i]).abs() > SYMMETRY_TOL * scale {
                return Err(ConductionError::Conductivity {
                    what: format!(
                        "the conductivity tensor is not symmetric: k[{i}][{j}] = {value} vs \
                         k[{j}][{i}] = {}",
                        k[j][i]
                    ),
                });
            }
        }
    }
    // Sylvester's criterion on the leading principal minors.
    let m1 = k[0][0];
    let m2 = k[0][0].mul_add(k[1][1], -(k[0][1] * k[1][0]));
    let m3 = k[0][0].mul_add(
        k[1][1].mul_add(k[2][2], -(k[1][2] * k[2][1])),
        -k[0][1].mul_add(
            k[1][0].mul_add(k[2][2], -(k[1][2] * k[2][0])),
            -k[0][2] * k[1][0].mul_add(k[2][1], -(k[1][1] * k[2][0])),
        ),
    );
    if !(m1 > 0.0 && m2 > 0.0 && m3 > 0.0) {
        return Err(ConductionError::Conductivity {
            what: format!(
                "the conductivity tensor is not positive definite (leading principal \
                 minors {m1:e}, {m2:e}, {m3:e} must all be positive)"
            ),
        });
    }
    Ok(())
}

fn check_orthonormal(axes: &[[f64; 3]; 3]) -> Result<(), ConductionError> {
    for i in 0..3 {
        for j in 0..3 {
            let d = axes[i][0].mul_add(
                axes[j][0],
                axes[i][1].mul_add(axes[j][1], axes[i][2] * axes[j][2]),
            );
            let want = if i == j { 1.0 } else { 0.0 };
            if (d - want).abs() > ORTHONORMALITY_TOL {
                return Err(ConductionError::Conductivity {
                    what: format!(
                        "principal axes are not orthonormal: e{i}·e{j} = {d}, expected {want}"
                    ),
                });
            }
        }
    }
    Ok(())
}

impl ConductivityModel {
    /// An isotropic model from one scalar table (declared or sampled).
    #[must_use]
    pub fn isotropic(table: ConductivityTable) -> ConductivityModel {
        ConductivityModel {
            kind: ConductivityKind::Isotropic(Box::new(table)),
        }
    }

    /// A caller-declared constant isotropic conductivity, W/(m·K).
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] for a non-positive value.
    pub fn isotropic_declared(k: f64) -> Result<ConductivityModel, ConductionError> {
        Ok(ConductivityModel::isotropic(ConductivityTable::declared(
            k,
        )?))
    }

    /// A general constant anisotropic tensor, W/(m·K). Checked symmetric
    /// and positive definite; carries no material provenance.
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] when the tensor is non-finite,
    /// asymmetric, or not positive definite.
    pub fn constant_tensor(k: [[f64; 3]; 3]) -> Result<ConductivityModel, ConductionError> {
        check_spd(&k)?;
        Ok(ConductivityModel {
            kind: ConductivityKind::ConstantTensor(Box::new(k)),
        })
    }

    /// An orthotropic model: one table per principal direction, over an
    /// orthonormal frame whose rows are the principal axes.
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] when the frame is not
    /// orthonormal.
    pub fn orthotropic(
        axes: [[f64; 3]; 3],
        tables: [ConductivityTable; 3],
    ) -> Result<ConductivityModel, ConductionError> {
        check_orthonormal(&axes)?;
        Ok(ConductivityModel {
            kind: ConductivityKind::Orthotropic {
                axes: Box::new(axes),
                tables: Box::new(tables),
            },
        })
    }

    /// The tensor at a temperature.
    ///
    /// # Errors
    /// [`ConductionError::OutsideTemperatureSpan`] when a table refuses.
    pub fn tensor_at(&self, temperature: f64) -> Result<[[f64; 3]; 3], ConductionError> {
        match &self.kind {
            ConductivityKind::ConstantTensor(k) => Ok(**k),
            ConductivityKind::Isotropic(table) => {
                let k = table.eval(temperature)?;
                Ok([[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]])
            }
            ConductivityKind::Orthotropic { axes, tables } => {
                let mut out = [[0.0f64; 3]; 3];
                for (a, table) in axes.iter().zip(tables.iter()) {
                    let k = table.eval(temperature)?;
                    for i in 0..3 {
                        for j in 0..3 {
                            out[i][j] = (k * a[i]).mul_add(a[j], out[i][j]);
                        }
                    }
                }
                Ok(out)
            }
        }
    }

    /// `dk/dT` as a tensor — the Newton Jacobian's material term.
    ///
    /// # Errors
    /// [`ConductionError::OutsideTemperatureSpan`] when a table refuses.
    pub fn tensor_derivative_at(&self, temperature: f64) -> Result<[[f64; 3]; 3], ConductionError> {
        match &self.kind {
            ConductivityKind::ConstantTensor(_) => Ok([[0.0f64; 3]; 3]),
            ConductivityKind::Isotropic(table) => {
                let d = table.derivative(temperature)?;
                Ok([[d, 0.0, 0.0], [0.0, d, 0.0], [0.0, 0.0, d]])
            }
            ConductivityKind::Orthotropic { axes, tables } => {
                let mut out = [[0.0f64; 3]; 3];
                for (a, table) in axes.iter().zip(tables.iter()) {
                    let d = table.derivative(temperature)?;
                    for i in 0..3 {
                        for j in 0..3 {
                            out[i][j] = (d * a[i]).mul_add(a[j], out[i][j]);
                        }
                    }
                }
                Ok(out)
            }
        }
    }

    /// True when any component varies with temperature (i.e. the solve
    /// is nonlinear).
    #[must_use]
    pub fn is_temperature_dependent(&self) -> bool {
        match &self.kind {
            ConductivityKind::ConstantTensor(_) => false,
            ConductivityKind::Isotropic(table) => table.is_temperature_dependent(),
            ConductivityKind::Orthotropic { tables, .. } => tables
                .iter()
                .any(ConductivityTable::is_temperature_dependent),
        }
    }

    /// The temperature interval every component is usable over.
    #[must_use]
    pub fn temperature_span(&self) -> TemperatureSpan {
        match &self.kind {
            ConductivityKind::ConstantTensor(_) => TemperatureSpan::Unbounded,
            ConductivityKind::Isotropic(table) => table.span(),
            ConductivityKind::Orthotropic { tables, .. } => tables
                .iter()
                .fold(TemperatureSpan::Unbounded, |acc, t| acc.intersect(t.span())),
        }
    }

    /// Every retained `fs-matdb` receipt behind this model, in declared
    /// component order then grid order.
    #[must_use]
    pub fn receipts(&self) -> Vec<&PropertyUsageReceipt> {
        match &self.kind {
            ConductivityKind::ConstantTensor(_) => Vec::new(),
            ConductivityKind::Isotropic(table) => table.receipts().iter().collect(),
            ConductivityKind::Orthotropic { tables, .. } => {
                tables.iter().flat_map(|t| t.receipts().iter()).collect()
            }
        }
    }

    /// The model's provenance class: [`ProvenanceClass::MatdbReceipts`]
    /// only when EVERY component carries receipts.
    #[must_use]
    pub fn provenance(&self) -> ProvenanceClass {
        let all_backed = match &self.kind {
            ConductivityKind::ConstantTensor(_) => false,
            ConductivityKind::Isotropic(table) => {
                table.provenance() == ProvenanceClass::MatdbReceipts
            }
            ConductivityKind::Orthotropic { tables, .. } => tables
                .iter()
                .all(|t| t.provenance() == ProvenanceClass::MatdbReceipts),
        };
        if all_backed {
            ProvenanceClass::MatdbReceipts
        } else {
            ProvenanceClass::Declared
        }
    }
}
