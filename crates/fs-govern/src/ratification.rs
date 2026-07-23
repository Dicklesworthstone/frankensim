//! The vertical ratification decision record (bead f85xj.1.4): the program's
//! ONE ratified vertical, recorded the way this project records claims — with
//! inputs, alternatives, a kill criterion bound to the MEASURED cycle-time
//! baseline, named mechanically-evaluable falsifiers, and a review date — so
//! a future re-evaluation is a data update, not an archaeology project.
//!
//! The record is fail-closed: [`ratified_vertical`] re-validates everything
//! it cites before handing the record out, including recomputing the measured
//! comparison in `fs-wedge` and checking the baseline's provenance class. A
//! ratification whose scoring table drifted, whose baseline is a placeholder,
//! or whose falsifiers are incomplete refuses with a typed error instead of
//! standing on stale authority (P8, Governance Rule 2).

use fs_wedge::{
    BaselineProvenance, CHT_BASELINE, CycleTimeBaseline, ScoringError, comparison_candidates,
    default_recommendation,
};

use crate::json_escape;
use core::fmt::Write as _;

/// One named falsifier: an observed fact that would force re-selection,
/// stated so a quarterly review can evaluate it mechanically.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Falsifier {
    /// Stable falsifier id.
    pub id: &'static str,
    /// The observed fact that forces re-selection.
    pub statement: &'static str,
    /// The mechanical evaluation a quarterly review performs.
    pub measurement: &'static str,
    /// The numeric or categorical trigger.
    pub threshold: &'static str,
}

impl Falsifier {
    /// Is every field populated?
    #[must_use]
    pub fn is_complete(self) -> bool {
        !self.id.trim().is_empty()
            && !self.statement.trim().is_empty()
            && !self.measurement.trim().is_empty()
            && !self.threshold.trim().is_empty()
    }
}

/// One recorded candidate total from the ratifying comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordedTotal {
    /// Candidate slug.
    pub candidate: &'static str,
    /// Weighted total recorded at ratification time.
    pub weighted_total: u16,
}

/// The program-level vertical ratification decision record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VerticalRatification {
    /// Stable record id, cited by downstream beads.
    pub id: &'static str,
    /// Decision date.
    pub decided_on: &'static str,
    /// The ratified vertical (candidate slug from the measured comparison).
    pub chosen_vertical: &'static str,
    /// Rank-2 candidate whose case is retained.
    pub runner_up: &'static str,
    /// The runner-up's retained strongest case, verbatim from the comparison.
    pub minority_report: &'static str,
    /// Git revision of the code inventory the scoring table rests on.
    pub scoring_inventory_revision: &'static str,
    /// Candidate totals recorded at ratification, in rank order.
    pub recorded_totals: &'static [RecordedTotal],
    /// Required cycle-time reduction factor, bound to the measured baseline.
    pub kill_target_reduction: f64,
    /// Quarters after GA to meet it or re-select, bound to the baseline.
    pub kill_within_quarters: u8,
    /// The named falsifiers.
    pub falsifiers: &'static [Falsifier],
    /// Next scheduled review.
    pub review_due: &'static str,
    /// Design principles this record operationalizes.
    pub principles: &'static [&'static str],
    /// Downstream beads gated on this record.
    pub downstream_gates: &'static [&'static str],
}

/// Why a ratification record refuses to stand.
#[derive(Debug, Clone, PartialEq)]
pub enum RatificationError {
    /// A required text field is empty.
    EmptyField {
        /// Which field.
        field: &'static str,
    },
    /// The record names no falsifier at all.
    MissingFalsifiers,
    /// A named falsifier is structurally incomplete.
    IncompleteFalsifier {
        /// The offending falsifier id (or its index when the id is empty).
        id: &'static str,
    },
    /// The kill criterion does not match the measured baseline record.
    KillCriterionUnbound {
        /// Which parameter disagrees.
        field: &'static str,
    },
    /// The cycle-time baseline in the decision path is not a measured record.
    BaselineNotMeasured {
        /// The refused provenance class label.
        provenance: &'static str,
    },
    /// The measured comparison cannot be recomputed.
    ScoringUnavailable {
        /// The underlying comparison refusal.
        source: ScoringError,
    },
    /// The recomputed comparison disagrees with the recorded decision.
    ScoringDrift {
        /// Which recorded fact drifted.
        field: &'static str,
        /// The recorded value.
        recorded: String,
        /// The recomputed value.
        recomputed: String,
    },
}

/// The ratification record for the 0.1 product vertical.
///
/// The chosen vertical is the measured comparison's rank-1 candidate:
/// thermal design assurance — the electronics-cooling/thermal family limited
/// to conduction, interfaces, radiation, and fan/correlation rungs, deferring
/// only the RANS rung. This is the "thermal design assurance" compromise the
/// program's default expectation named as an acceptable ratified outcome, and
/// it is selected BY the data, not despite it. The kill-criterion denominator
/// is the measured incumbent electronics-cooling workflow envelope in
/// `fs-wedge` (`CHT_BASELINE`), which spans the whole thermal family.
pub const VERTICAL_RATIFICATION_V1: VerticalRatification = VerticalRatification {
    id: "frankensim-vertical-ratification-v1",
    decided_on: "2026-07-22",
    chosen_vertical: "thermal-design-assurance",
    runner_up: "sdf-structural-topology-assurance",
    minority_report: "Structural/topology assurance already has the deepest verified kernel \
        stack in the workspace (fs-solid, fs-topopt, fs-truss-e2e, fs-frame) and public \
        validation data with retained raw records; if thermal data acquisition stalls, it is \
        the strongest fallback.",
    scoring_inventory_revision: "b3b5f2c1c809eec06cde1e40cbc916d6995469b5",
    recorded_totals: &[
        RecordedTotal {
            candidate: "thermal-design-assurance",
            weighted_total: 638,
        },
        RecordedTotal {
            candidate: "sdf-structural-topology-assurance",
            weighted_total: 623,
        },
        RecordedTotal {
            candidate: "full-electronics-cooling-cht",
            weighted_total: 502,
        },
    ],
    kill_target_reduction: 3.0,
    kill_within_quarters: 2,
    falsifiers: &[
        Falsifier {
            id: "level-c-thermal-data-unobtainable",
            statement: "Level-C electronics-thermal experimental datasets with retained raw \
                records and stated uncertainties prove unobtainable, so the vertical cannot \
                ever earn an L4 validation claim.",
            measurement: "Count fs-vvreg corpus rows at Level C for the thermal family whose \
                provenance includes original raw data and uncertainty statements (the current \
                Martin-Moyce curve is derived-only and does not count).",
            threshold: "Zero qualifying rows by 2027-01-22 (two quarterly reviews) forces \
                re-selection or a recorded scope change.",
        },
        Falsifier {
            id: "kill-criterion-not-met",
            statement: "The measured FrankenSim iteration cycle time fails the required \
                reduction against the measured incumbent envelope.",
            measurement: "fs_wedge::CHT_BASELINE.evaluate_kill_criterion(measured_days) on the \
                instrumented per-iteration cycle time of the reference cooling project.",
            threshold: "Verdict NotMet at GA plus kill_within_quarters ends the wedge; verdict \
                Indeterminate at that point without an ExecutedRun baseline upgrade in flight \
                escalates to a re-selection review. The reduction claim is never marketed \
                until the verdict is Met.",
        },
        Falsifier {
            id: "scoring-reversal-on-reinventory",
            statement: "A quarterly re-inventory of the measured comparison inputs reverses \
                the ranking that ratified this vertical.",
            measurement: "Refresh the fs-wedge measured inputs, then recompute \
                fs_wedge::default_recommendation() under the default weights.",
            threshold: "recommended != thermal-design-assurance on two consecutive quarterly \
                reviews forces re-ratification with a new record id.",
        },
    ],
    review_due: "2026-10-22",
    principles: &["P7", "P8"],
    downstream_gates: &[
        "frankensim-extreal-program-f85xj.6.1",
        "frankensim-extreal-program-f85xj.4.4",
        "frankensim-extreal-program-f85xj.10.4",
        "frankensim-extreal-program-f85xj.5.8",
    ],
};

impl VerticalRatification {
    /// Validate this record against the measured world it cites, using the
    /// supplied cycle-time baseline as the kill-criterion denominator.
    ///
    /// Field checks run first, then falsifiers, then the baseline binding,
    /// then the scoring recomputation, so the first refusal names the
    /// earliest broken layer.
    pub fn validate_against(&self, baseline: &CycleTimeBaseline) -> Result<(), RatificationError> {
        self.validate_fields()?;
        self.validate_falsifiers()?;
        self.validate_kill_binding(baseline)?;
        self.validate_scoring()
    }

    /// Validate against the workspace's actual decision-path baseline.
    pub fn validate(&self) -> Result<(), RatificationError> {
        self.validate_against(&CHT_BASELINE)
    }

    fn validate_fields(&self) -> Result<(), RatificationError> {
        let fields = [
            ("id", self.id),
            ("decided_on", self.decided_on),
            ("chosen_vertical", self.chosen_vertical),
            ("runner_up", self.runner_up),
            ("minority_report", self.minority_report),
            (
                "scoring_inventory_revision",
                self.scoring_inventory_revision,
            ),
            ("review_due", self.review_due),
        ];
        for (field, value) in fields {
            if value.trim().is_empty() {
                return Err(RatificationError::EmptyField { field });
            }
        }
        if self.recorded_totals.is_empty() {
            return Err(RatificationError::EmptyField {
                field: "recorded_totals",
            });
        }
        if self.principles.is_empty() {
            return Err(RatificationError::EmptyField {
                field: "principles",
            });
        }
        if self.downstream_gates.is_empty() {
            return Err(RatificationError::EmptyField {
                field: "downstream_gates",
            });
        }
        Ok(())
    }

    fn validate_falsifiers(&self) -> Result<(), RatificationError> {
        if self.falsifiers.is_empty() {
            return Err(RatificationError::MissingFalsifiers);
        }
        for falsifier in self.falsifiers {
            if !falsifier.is_complete() {
                return Err(RatificationError::IncompleteFalsifier { id: falsifier.id });
            }
        }
        Ok(())
    }

    fn validate_kill_binding(&self, baseline: &CycleTimeBaseline) -> Result<(), RatificationError> {
        if baseline.provenance == BaselineProvenance::Placeholder || !baseline.is_complete() {
            return Err(RatificationError::BaselineNotMeasured {
                provenance: baseline.provenance.label(),
            });
        }
        if (self.kill_target_reduction - baseline.target_reduction).abs() >= f64::EPSILON {
            return Err(RatificationError::KillCriterionUnbound {
                field: "kill_target_reduction",
            });
        }
        if self.kill_within_quarters != baseline.kill_within_quarters {
            return Err(RatificationError::KillCriterionUnbound {
                field: "kill_within_quarters",
            });
        }
        Ok(())
    }

    fn validate_scoring(&self) -> Result<(), RatificationError> {
        let recommendation = default_recommendation()
            .map_err(|source| RatificationError::ScoringUnavailable { source })?;
        if recommendation.recommended != self.chosen_vertical {
            return Err(RatificationError::ScoringDrift {
                field: "chosen_vertical",
                recorded: self.chosen_vertical.to_string(),
                recomputed: recommendation.recommended.to_string(),
            });
        }
        if recommendation.runner_up != self.runner_up {
            return Err(RatificationError::ScoringDrift {
                field: "runner_up",
                recorded: self.runner_up.to_string(),
                recomputed: recommendation.runner_up.to_string(),
            });
        }
        if recommendation.ranked.len() != self.recorded_totals.len() {
            return Err(RatificationError::ScoringDrift {
                field: "recorded_totals.len",
                recorded: self.recorded_totals.len().to_string(),
                recomputed: recommendation.ranked.len().to_string(),
            });
        }
        for (recorded, recomputed) in self.recorded_totals.iter().zip(&recommendation.ranked) {
            if recorded.candidate != recomputed.candidate
                || recorded.weighted_total != recomputed.weighted_total
            {
                return Err(RatificationError::ScoringDrift {
                    field: "recorded_totals",
                    recorded: format!("{}={}", recorded.candidate, recorded.weighted_total),
                    recomputed: format!("{}={}", recomputed.candidate, recomputed.weighted_total),
                });
            }
        }
        for candidate in comparison_candidates() {
            if candidate.inventory_revision != self.scoring_inventory_revision {
                return Err(RatificationError::ScoringDrift {
                    field: "scoring_inventory_revision",
                    recorded: self.scoring_inventory_revision.to_string(),
                    recomputed: candidate.inventory_revision.to_string(),
                });
            }
        }
        Ok(())
    }

    /// Render the record as one deterministic JSON object. Callers wanting a
    /// fail-closed render should use [`ratification_json`].
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::new();
        write!(
            out,
            "{{\"id\":\"{}\",\"decided_on\":\"{}\",\"chosen_vertical\":\"{}\",\"runner_up\":\"{}\",\"minority_report\":\"{}\",\"scoring_inventory_revision\":\"{}\",\"recorded_totals\":[",
            json_escape(self.id),
            json_escape(self.decided_on),
            json_escape(self.chosen_vertical),
            json_escape(self.runner_up),
            json_escape(self.minority_report),
            json_escape(self.scoring_inventory_revision),
        )
        .expect("write to String");
        for (index, total) in self.recorded_totals.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            write!(
                out,
                "{{\"candidate\":\"{}\",\"weighted_total\":{}}}",
                json_escape(total.candidate),
                total.weighted_total
            )
            .expect("write to String");
        }
        write!(
            out,
            "],\"kill_target_reduction\":{},\"kill_within_quarters\":{},\"falsifiers\":[",
            self.kill_target_reduction, self.kill_within_quarters
        )
        .expect("write to String");
        for (index, falsifier) in self.falsifiers.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            write!(
                out,
                "{{\"id\":\"{}\",\"statement\":\"{}\",\"measurement\":\"{}\",\"threshold\":\"{}\"}}",
                json_escape(falsifier.id),
                json_escape(falsifier.statement),
                json_escape(falsifier.measurement),
                json_escape(falsifier.threshold)
            )
            .expect("write to String");
        }
        write!(
            out,
            "],\"review_due\":\"{}\",\"principles\":[",
            json_escape(self.review_due)
        )
        .expect("write to String");
        for (index, principle) in self.principles.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            write!(out, "\"{}\"", json_escape(principle)).expect("write to String");
        }
        out.push_str("],\"downstream_gates\":[");
        for (index, gate) in self.downstream_gates.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            write!(out, "\"{}\"", json_escape(gate)).expect("write to String");
        }
        out.push_str("]}");
        out
    }
}

/// Every program-level decision record, in decision order.
#[must_use]
pub fn decision_records() -> &'static [VerticalRatification] {
    &[VERTICAL_RATIFICATION_V1]
}

/// Fail-closed accessor: the ratified vertical, only if the record still
/// validates against the measured comparison and the measured baseline.
pub fn ratified_vertical() -> Result<&'static VerticalRatification, RatificationError> {
    VERTICAL_RATIFICATION_V1.validate()?;
    Ok(&VERTICAL_RATIFICATION_V1)
}

/// Fail-closed JSON render of every decision record.
pub fn ratification_json() -> Result<String, RatificationError> {
    let mut out = String::from("{\"decision_records\":[");
    for (index, record) in decision_records().iter().enumerate() {
        record.validate()?;
        if index > 0 {
            out.push(',');
        }
        out.push_str(&record.to_json());
    }
    out.push_str("]}");
    Ok(out)
}
