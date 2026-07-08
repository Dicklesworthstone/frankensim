//! fs-report — automatic lab notebooks + semantic design diffs. Layer: L6.
//!
//! Reproducibility should be a SIDE EFFECT of running a study, not a virtue you
//! remember to practice. A [`LabNotebook`] is the automatic lab notebook: every
//! study emits a deterministic, human-readable report — provenance, prose,
//! Qty-labelled metrics (units on every value, P10), AND THE EXACT IR TO
//! REPRODUCE IT ([`LabNotebook::repro_ir`]). Because the render is deterministic
//! it is CONTENT-ADDRESSED ([`LabNotebook::content_hash`]), so replaying the IR
//! and re-rendering yields the same hash — the reproducibility loop closes by
//! construction.
//!
//! [`semantic_diff`] is the other half: a diff between two designs that is a
//! GEOMETRIC attribution ("lip curvature −18%, wall thinned 0.4 mm"), ranked by
//! significance — not a file diff. Deterministic; no dependencies.

use core::fmt::Write as _;
use std::collections::BTreeMap;

/// A dimensioned quantity — a value with its unit (units on every value).
#[derive(Debug, Clone, PartialEq)]
pub struct Quantity {
    /// The numeric value.
    pub value: f64,
    /// The unit label (e.g. `"mm"`, `"kg"`, `"1/mm"`).
    pub unit: String,
}

impl Quantity {
    /// A quantity.
    #[must_use]
    pub fn new(value: f64, unit: impl Into<String>) -> Quantity {
        Quantity {
            value,
            unit: unit.into(),
        }
    }
}

/// One replayable operation of the reproducibility IR.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReproStep {
    /// The operation name.
    pub op: String,
    /// Its serialized arguments.
    pub args: Vec<String>,
}

/// A notebook block.
#[derive(Debug, Clone, PartialEq)]
pub enum Block {
    /// Free prose.
    Prose(String),
    /// A named, dimensioned metric.
    Metric {
        /// The metric name.
        name: String,
        /// The value + unit.
        quantity: Quantity,
    },
    /// A reproducibility step (part of the replay IR).
    Step(ReproStep),
}

/// An automatic lab notebook for a study.
#[derive(Debug, Clone, PartialEq)]
pub struct LabNotebook {
    /// The study title.
    pub title: String,
    /// The RNG seed (provenance).
    pub seed: u64,
    /// The toolchain / crate version (provenance).
    pub version: String,
    /// The report body.
    pub blocks: Vec<Block>,
}

impl LabNotebook {
    /// A new notebook with provenance.
    #[must_use]
    pub fn new(title: impl Into<String>, seed: u64, version: impl Into<String>) -> LabNotebook {
        LabNotebook {
            title: title.into(),
            seed,
            version: version.into(),
            blocks: Vec::new(),
        }
    }

    /// Append prose.
    pub fn prose(&mut self, text: impl Into<String>) -> &mut LabNotebook {
        self.blocks.push(Block::Prose(text.into()));
        self
    }

    /// Append a dimensioned metric.
    pub fn metric(
        &mut self,
        name: impl Into<String>,
        value: f64,
        unit: impl Into<String>,
    ) -> &mut LabNotebook {
        self.blocks.push(Block::Metric {
            name: name.into(),
            quantity: Quantity::new(value, unit),
        });
        self
    }

    /// Append a reproducibility step.
    pub fn step(&mut self, op: impl Into<String>, args: Vec<String>) -> &mut LabNotebook {
        self.blocks.push(Block::Step(ReproStep {
            op: op.into(),
            args,
        }));
        self
    }

    /// The metrics recorded (name + quantity).
    #[must_use]
    pub fn metrics(&self) -> Vec<(&str, &Quantity)> {
        self.blocks
            .iter()
            .filter_map(|b| match b {
                Block::Metric { name, quantity } => Some((name.as_str(), quantity)),
                _ => None,
            })
            .collect()
    }

    /// THE EXACT IR TO REPRODUCE the study — the ordered replay steps.
    #[must_use]
    pub fn repro_ir(&self) -> Vec<ReproStep> {
        self.blocks
            .iter()
            .filter_map(|b| match b {
                Block::Step(s) => Some(s.clone()),
                _ => None,
            })
            .collect()
    }

    /// The report rendered to Markdown (deterministic).
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "# {}", self.title);
        let _ = writeln!(s);
        let _ = writeln!(s, "_seed: {} · version: {}_", self.seed, self.version);
        let _ = writeln!(s);
        for block in &self.blocks {
            match block {
                Block::Prose(t) => {
                    let _ = writeln!(s, "{t}");
                    let _ = writeln!(s);
                }
                Block::Metric { name, quantity } => {
                    let _ = writeln!(s, "- **{}**: {} {}", name, quantity.value, quantity.unit);
                }
                Block::Step(step) => {
                    let _ = writeln!(s, "- repro: `{}({})`", step.op, step.args.join(", "));
                }
            }
        }
        s
    }

    /// A content hash of the rendered report — a report is as content-addressed
    /// as any other ledger artifact (FNV-1a of the Markdown).
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        fnv1a(self.render_markdown().as_bytes())
    }
}

/// A per-feature semantic difference between two designs.
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureDelta {
    /// The feature name.
    pub name: String,
    /// The value before.
    pub before: f64,
    /// The value after.
    pub after: f64,
    /// The absolute change (`after − before`).
    pub abs_change: f64,
    /// The relative change (`abs_change / before`; `0` if `before == 0`).
    pub rel_change: f64,
    /// The unit.
    pub unit: String,
}

impl FeatureDelta {
    /// A human attribution string, e.g. `"wall_thickness: 2 mm → 1.6 mm (−20.0%)"`.
    #[must_use]
    pub fn describe(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{}: {} {} → {} {} ({:+.1}%)",
            self.name,
            self.before,
            self.unit,
            self.after,
            self.unit,
            self.rel_change * 100.0
        );
        s
    }
}

/// A SEMANTIC (per-feature) diff between two designs described as
/// `feature → Quantity` maps: the changed features with absolute + relative
/// deltas, ranked by significance (largest relative change first). Not a file
/// diff — a geometric attribution.
#[must_use]
pub fn semantic_diff(
    before: &BTreeMap<String, Quantity>,
    after: &BTreeMap<String, Quantity>,
) -> Vec<FeatureDelta> {
    let mut deltas: Vec<FeatureDelta> = before
        .iter()
        .filter_map(|(name, b)| {
            after.get(name).map(|a| {
                let abs_change = a.value - b.value;
                let rel_change = if b.value == 0.0 {
                    0.0
                } else {
                    abs_change / b.value
                };
                FeatureDelta {
                    name: name.clone(),
                    before: b.value,
                    after: a.value,
                    abs_change,
                    rel_change,
                    unit: b.unit.clone(),
                }
            })
        })
        .collect();
    // rank by significance (largest |relative change| first); name as tiebreak.
    deltas.sort_by(|x, y| {
        y.rel_change
            .abs()
            .partial_cmp(&x.rel_change.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| x.name.cmp(&y.name))
    });
    deltas
}

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}
