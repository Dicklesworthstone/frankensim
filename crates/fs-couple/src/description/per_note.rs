//! Bounded model-only threshold search (G1).
//!
//! A model supplies a predeclared, monotone onset predicate and its identity,
//! including horizon/observable/tolerance. No recording or spectral reference
//! is accepted. The receipt certifies the evaluated bracket, not global
//! monotonicity or an infinite-time bifurcation. Those remain model obligations.
use fs_blake3::{ContentHash, hash_bytes};
use fs_qty::{Dims, QtyAny};
use fs_scenario::performance::{DeclaredObjective, PerNoteControl, Performance, PerformanceError};

const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);

/// Deliberately has no audio/reference argument. Model identity must bind all
/// initial conditions and the predeclared physical onset observable.
pub trait ThresholdModel {
    fn identity(&self) -> ContentHash;
    fn oscillates(&mut self, pressure_pa: f64) -> Result<bool, String>;
}
#[derive(Clone, Copy, Debug)]
pub struct ThresholdBudget {
    pub pressure_tolerance_pa: f64,
    /// Includes the two endpoint evaluations.
    pub max_evaluations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub enum ResolutionError {
    Invalid(&'static str),
    UnsupportedObjective,
    Unbracketed,
    BudgetExceeded,
    Model(String),
    Performance(PerformanceError),
}
impl std::fmt::Display for ResolutionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "per-note resolution: {self:?}")
    }
}
impl std::error::Error for ResolutionError {}
#[derive(Clone, Debug, PartialEq)]
pub struct ThresholdReceipt {
    pub model: ContentHash,
    pub target: fs_scenario::assembly::ParameterPath,
    pub note_index: usize,
    pub domain_pa: (f64, f64),
    pub bracket_pa: (f64, f64),
    pub tolerance_pa: f64,
    pub max_evaluations: usize,
    pub probes: Vec<(f64, bool)>,
}
impl ThresholdReceipt {
    /// Domain-separated deterministic evaluation preimage, including failures
    /// to oscillate, the objective, search budget and exact returned bracket.
    pub fn content_hash(&self) -> ContentHash {
        let mut b = b"fs-minimum-threshold-pressure-v1".to_vec();
        b.extend_from_slice(self.model.as_bytes());
        let target = format!(
            "{}:{}",
            self.target.entity.token(),
            self.target.parameter.name
        );
        b.extend_from_slice(&(target.len() as u64).to_be_bytes());
        b.extend_from_slice(target.as_bytes());
        b.extend(self.target.parameter.dims.0.iter().map(|v| *v as u8));
        b.extend_from_slice(&(self.note_index as u64).to_be_bytes());
        b.extend_from_slice(&(self.max_evaluations as u64).to_be_bytes());
        for v in [
            self.domain_pa.0,
            self.domain_pa.1,
            self.bracket_pa.0,
            self.bracket_pa.1,
            self.tolerance_pa,
        ] {
            b.extend_from_slice(&v.to_bits().to_be_bytes());
        }
        for (v, on) in &self.probes {
            b.extend_from_slice(&v.to_bits().to_be_bytes());
            b.push(u8::from(*on));
        }
        hash_bytes(&b)
    }
    pub fn resolved(&self) -> QtyAny {
        QtyAny::new(self.bracket_pa.1, PRESSURE)
    }
}

/// Deterministic bisection, retaining a non-oscillating lower and oscillating
/// upper bound. Returns the upper bound within the declared absolute tolerance.
pub fn resolve_threshold<M: ThresholdModel>(
    model: &mut M,
    control: &PerNoteControl,
    budget: ThresholdBudget,
) -> Result<ThresholdReceipt, ResolutionError> {
    if !matches!(
        control.objective,
        DeclaredObjective::MinimumThresholdPressure
    ) {
        return Err(ResolutionError::UnsupportedObjective);
    }
    let (a, b) = control.domain;
    if a.dims != PRESSURE
        || b.dims != PRESSURE
        || control.target.parameter.dims != PRESSURE
        || !a.value.is_finite()
        || !b.value.is_finite()
        || a.value < 0.0
        || a.value >= b.value
        || !budget.pressure_tolerance_pa.is_finite()
        || budget.pressure_tolerance_pa <= 0.0
        || budget.max_evaluations < 2
    {
        return Err(ResolutionError::Invalid("pressure domain or budget"));
    }
    let mut receipt = ThresholdReceipt {
        model: model.identity(),
        target: control.target.clone(),
        note_index: control.note_index,
        domain_pa: (a.value, b.value),
        bracket_pa: (a.value, b.value),
        tolerance_pa: budget.pressure_tolerance_pa,
        max_evaluations: budget.max_evaluations,
        probes: Vec::new(),
    };
    let mut probe = |p| -> Result<bool, ResolutionError> {
        if receipt.probes.len() >= budget.max_evaluations {
            return Err(ResolutionError::BudgetExceeded);
        }
        let on = model.oscillates(p).map_err(ResolutionError::Model)?;
        receipt.probes.push((p, on));
        Ok(on)
    };
    if probe(a.value)? || !probe(b.value)? {
        return Err(ResolutionError::Unbracketed);
    }
    let (mut lo, mut hi) = (a.value, b.value);
    while hi - lo > budget.pressure_tolerance_pa {
        let mid = lo + (hi - lo) * 0.5;
        if mid <= lo || mid >= hi {
            return Err(ResolutionError::Invalid(
                "pressure tolerance below float resolution",
            ));
        }
        if probe(mid)? {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    receipt.bracket_pa = (lo, hi);
    Ok(receipt)
}

/// Attach the actual evaluation receipt to the performance, so lowering copies
/// it into render provenance. A receipt for another target/note/domain refuses.
pub fn record_resolution(
    performance: &mut Performance,
    index: usize,
    receipt: &ThresholdReceipt,
) -> Result<(), ResolutionError> {
    let c = performance
        .per_note()
        .get(index)
        .ok_or(ResolutionError::Invalid("per-note index"))?;
    if c.target != receipt.target
        || c.note_index != receipt.note_index
        || !matches!(c.objective, DeclaredObjective::MinimumThresholdPressure)
        || c.domain.0.value.to_bits() != receipt.domain_pa.0.to_bits()
        || c.domain.1.value.to_bits() != receipt.domain_pa.1.to_bits()
        || c.domain.0.dims != PRESSURE
        || c.domain.1.dims != PRESSURE
    {
        return Err(ResolutionError::Invalid(
            "receipt does not belong to control",
        ));
    }
    performance
        .resolve(index, receipt.resolved(), receipt.content_hash())
        .map_err(ResolutionError::Performance)
}
