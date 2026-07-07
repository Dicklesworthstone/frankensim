//! The `RegimeReport`: groups + dominant balance + model verdicts +
//! recommended scaling + similarity match, returned as
//! `Evidence<RegimeReport>` (exact numerics — the report is arithmetic
//! over declared inputs — with content-addressed provenance).

use crate::RegimeError;
use crate::cards::{distance_to_validity, flux_model_cards};
use crate::groups::{RoleInput, standard_groups};
use crate::pi::{Input, pi_groups};
use crate::scaling::ScalingMap;
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate,
};
use fs_math::det;
use std::collections::BTreeMap;
use std::fmt::Write as _;

fn log10(x: f64) -> f64 {
    det::ln(x) / core::f64::consts::LN_10
}

/// A canonical benchmark neighbor ("your Re=98 cylinder is near the
/// Re=100 benchmark; expected Cd ≈ 1.33").
#[derive(Debug, Clone, PartialEq)]
pub struct BenchmarkMatch {
    /// Benchmark name.
    pub name: &'static str,
    /// Log-space distance over shared groups.
    pub distance: f64,
    /// What to expect there.
    pub expectation: &'static str,
    /// "info" (close), "warning" (nearby, extrapolating), "far".
    pub grade: &'static str,
}

struct Benchmark {
    name: &'static str,
    groups: &'static [(&'static str, f64)],
    expectation: &'static str,
}

const BENCHMARKS: &[Benchmark] = &[
    Benchmark {
        name: "cylinder-crossflow-Re100",
        groups: &[("Re", 100.0)],
        expectation: "expected Cd ~ 1.33, Karman shedding at St ~ 0.16",
    },
    Benchmark {
        name: "stokes-sphere-Re0.1",
        groups: &[("Re", 0.1)],
        expectation: "Stokes drag Cd = 24/Re within ~2%",
    },
    Benchmark {
        name: "lid-driven-cavity-Re1000",
        groups: &[("Re", 1000.0)],
        expectation: "primary vortex center at (0.531, 0.565) (Ghia et al.)",
    },
    Benchmark {
        name: "dam-break-Fr1",
        groups: &[("Fr", 1.0)],
        expectation: "surge front celerity ~ 2*sqrt(g*h0) (Ritter)",
    },
];

/// The regime assessment.
#[derive(Debug, Clone, PartialEq)]
pub struct RegimeReport {
    /// Named dimensionless groups.
    pub groups: BTreeMap<String, f64>,
    /// Rank of the dimension matrix over the inputs.
    pub pi_rank: usize,
    /// Number of independent Pi groups (Buckingham count).
    pub pi_count: usize,
    /// One-line dominant-balance reading of the group values.
    pub dominant_balance: String,
    /// Models whose validity boxes contain the group point.
    pub valid_models: Vec<String>,
    /// Refused models with the violated-bound reasons.
    pub invalid_models: Vec<(String, String)>,
    /// Recommended nondimensionalization scales.
    pub recommended_scaling: ScalingMap,
    /// log10 spread of the inputs' scale factors — a conditioning
    /// red-flag when large.
    pub conditioning_risk: f64,
    /// Nearest canonical benchmark, when one is close enough to teach.
    pub nearest_benchmark: Option<BenchmarkMatch>,
}

fn dominant_balance(groups: &BTreeMap<String, f64>) -> String {
    if let Some(&re) = groups.get("Re") {
        if re < 1.0 {
            return format!("viscous-dominated (creeping, Re = {re:.3e})");
        }
        if re < 2300.0 {
            return format!("laminar inertial-viscous balance (Re = {re:.3e})");
        }
        return format!("turbulent inertia-dominated (Re = {re:.3e})");
    }
    if let Some(&sl) = groups.get("slenderness") {
        if sl >= 20.0 {
            return format!("flexure-dominated slender member (L/r = {sl:.1})");
        }
        return format!("shear-deformation-significant member (L/r = {sl:.1})");
    }
    "no dominant-balance rule matched the available groups".to_string()
}

fn nearest_benchmark(groups: &BTreeMap<String, f64>) -> Option<BenchmarkMatch> {
    let mut best: Option<BenchmarkMatch> = None;
    for b in BENCHMARKS {
        let mut d = 0.0f64;
        let mut shared = 0usize;
        for &(name, target) in b.groups {
            if let Some(&v) = groups.get(name)
                && v > 0.0
                && target > 0.0
            {
                d += log10(v / target).abs();
                shared += 1;
            }
        }
        if shared == 0 {
            continue;
        }
        let grade = if d < 0.05 {
            "info"
        } else if d < 0.7 {
            "warning"
        } else {
            "far"
        };
        let m = BenchmarkMatch {
            name: b.name,
            distance: d,
            expectation: b.expectation,
            grade,
        };
        if best.as_ref().is_none_or(|cur| m.distance < cur.distance) {
            best = Some(m);
        }
    }
    best
}

/// Assess the regime of a role-tagged problem: groups, Pi count,
/// model verdicts, scaling, similarity — as `Evidence<RegimeReport>`.
///
/// # Errors
/// [`RegimeError`] on dimensional inconsistencies, degenerate inputs, or
/// a missing Length role (scaling needs one).
pub fn assess(inputs: &[RoleInput]) -> Result<Evidence<RegimeReport>, RegimeError> {
    let named = standard_groups(inputs)?;
    let mut groups = BTreeMap::new();
    for g in &named {
        groups.insert(g.name.to_string(), g.value);
    }
    let pi_inputs: Vec<Input> = inputs
        .iter()
        .map(|i| Input {
            name: i.role.tag().to_string(),
            qty: i.qty,
        })
        .collect();
    let basis = pi_groups(&pi_inputs)?;
    let scaling = ScalingMap::recommend(inputs)?;
    // Conditioning risk: spread of the inputs' nondimensionalization
    // factors in decades. Large spread = raw-SI assembly mixes scales.
    let mut min_f = f64::INFINITY;
    let mut max_f = 0.0f64;
    for i in inputs {
        let f = scaling.factor(i.qty.dims.0).abs().max(f64::MIN_POSITIVE);
        min_f = min_f.min(f);
        max_f = max_f.max(f);
    }
    let conditioning_risk = log10(max_f / min_f);
    let registry = flux_model_cards();
    let mut valid = Vec::new();
    let mut invalid = Vec::new();
    for card in &registry {
        if card.validity.contains(&groups) {
            valid.push(card.name.clone());
        } else {
            let d = distance_to_validity(card, &groups);
            invalid.push((
                card.name.clone(),
                format!("distance {d:.3} decades to validity box"),
            ));
        }
    }
    let report = RegimeReport {
        dominant_balance: dominant_balance(&groups),
        nearest_benchmark: nearest_benchmark(&groups),
        pi_rank: basis.rank,
        pi_count: basis.groups.len(),
        valid_models: valid,
        invalid_models: invalid,
        recommended_scaling: scaling,
        conditioning_risk,
        groups,
    };
    let mut canon = String::from("regime-report");
    for (k, v) in &report.groups {
        let _ = write!(canon, ";{k}={v}");
    }
    let _ = write!(canon, ";rank={};count={}", report.pi_rank, report.pi_count);
    let qoi = report.groups.get("Re").copied().unwrap_or(f64::from(
        u32::try_from(report.pi_count).unwrap_or(u32::MAX),
    ));
    Ok(Evidence {
        value: report,
        qoi,
        numerical: NumericalCertificate::exact(qoi),
        statistical: StatisticalCertificate::None,
        model: ModelEvidence::none(),
        sensitivity: SensitivitySummary::default(),
        provenance: ProvenanceHash::of_bytes(canon.as_bytes()),
        adjoint_ref: None,
    })
}
