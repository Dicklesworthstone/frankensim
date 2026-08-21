//! ABComparisonReceiptV1 (bead wf-root-guzez.7.1.3, E6.1-iii):
//! SameInputTrace + HumanRefly ghost semantics on the engine.
//!
//! SameInputTrace: run B replays run A's APPLIED input trace under a
//! modified axis (scenario scalar, assist, model mode). HumanRefly:
//! two different traces (the live re-fly vs the ghost). Either way the
//! receipt records the COMMON-PREFIX tick — the first tick whose
//! frozen 12-float snapshot payloads differ BITWISE — plus both
//! terminals (terminal-divergence cases included) and both final
//! chained digests. A comparison that cannot fail is not a
//! comparison: the battery drives identical twins (full-length
//! prefix), a perturbed-axis pair (receipted divergence tick), and a
//! terminal-divergence pair (different TerminalEvent kinds).

use crate::Refusal;
use crate::simloop::{ControlInput, Phase, PilotMode, SNAPSHOT_LEN, ScenarioSpec, SimLoop};
use fs_blake3::hash_domain;

/// Receipt schema id.
pub const AB_SCHEMA: &str = "org.frankensim.wf.ab-comparison-receipt.v1";

/// Trace cap (absurd-input guard; matches MAX_TICKS).
pub const MAX_TRACE_LEN: usize = 72_000;

/// One run's summary inside the receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct RunSummary {
    /// RunIntentId.
    pub run_intent_id: String,
    /// Terminal tick.
    pub terminal_tick: u64,
    /// Terminal phase code (snapshot payload convention, 2..=5).
    pub terminal_code: u8,
    /// Final chained digest.
    pub final_digest: String,
}

/// The A/B receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct AbComparisonReceipt {
    /// Schema id.
    pub schema: &'static str,
    /// Comparison mode word.
    pub mode: &'static str,
    /// Run A summary.
    pub a: RunSummary,
    /// Run B summary.
    pub b: RunSummary,
    /// Ticks whose payloads matched BITWISE from the start.
    pub common_prefix_ticks: u64,
    /// First divergent tick (None = identical through the shorter run
    /// AND equal terminals).
    pub first_divergence_tick: Option<u64>,
    /// The two terminal kinds differ (terminal-divergence case).
    pub terminal_divergence: bool,
    /// Receipt digest.
    pub receipt_digest: String,
}

/// March one run, feeding the trace when the mode is Human; collect
/// per-tick payloads. Deterministic modes take `trace = None`.
fn run_collect(
    spec: &ScenarioSpec,
    trace: Option<&[ControlInput]>,
) -> Result<(RunSummary, Vec<[f64; SNAPSHOT_LEN]>), Refusal> {
    if matches!(spec.pilot_mode, PilotMode::Human) != trace.is_some() {
        return Err(Refusal {
            code: "ab-trace-mode-mismatch",
            message: "Human runs need a trace; deterministic runs must not carry one".into(),
            ranked_repairs: vec!["bind the applied input trace to the Human run".into()],
        });
    }
    if let Some(t) = trace {
        if t.is_empty() || t.len() > MAX_TRACE_LEN {
            return Err(Refusal {
                code: "ab-trace-invalid",
                message: format!("trace length {} outside [1, {MAX_TRACE_LEN}]", t.len()),
                ranked_repairs: vec!["record the trace from a real session".into()],
            });
        }
    }
    let mut sim = SimLoop::init(spec.clone())?;
    let run_intent_id = sim.run_intent_id.clone();
    let mut payloads = Vec::new();
    let (terminal_tick, terminal_code) = loop {
        let input = trace.map(|t| {
            let k = payloads.len().min(t.len() - 1);
            t[k]
        });
        let out = sim.step(input)?;
        payloads.push(sim.snapshot_payload(&out));
        if let Phase::Ended(_) = out.phase {
            let code = payloads.last().expect("pushed")[SNAPSHOT_LEN - 1] as u8;
            break (out.tick, code);
        }
    };
    Ok((
        RunSummary {
            run_intent_id,
            terminal_tick,
            terminal_code,
            final_digest: sim.digest_hex(),
        },
        payloads,
    ))
}

/// Run the A/B comparison and emit the receipt.
///
/// # Errors
/// `ab-trace-mode-mismatch`, `ab-trace-invalid`; engine refusals pass
/// through.
pub fn ab_compare(
    mode: &'static str,
    spec_a: &ScenarioSpec,
    trace_a: Option<&[ControlInput]>,
    spec_b: &ScenarioSpec,
    trace_b: Option<&[ControlInput]>,
) -> Result<AbComparisonReceipt, Refusal> {
    let (a, pa) = run_collect(spec_a, trace_a)?;
    let (b, pb) = run_collect(spec_b, trace_b)?;
    let shorter = pa.len().min(pb.len());
    let mut prefix = shorter as u64;
    let mut first_div = None;
    'outer: for k in 0..shorter {
        for slot in 0..SNAPSHOT_LEN {
            if pa[k][slot].to_bits() != pb[k][slot].to_bits() {
                prefix = k as u64;
                first_div = Some((k + 1) as u64); // ticks are 1-based
                break 'outer;
            }
        }
    }
    if first_div.is_none() && pa.len() != pb.len() {
        // Bitwise-equal through the shorter run but different lengths:
        // the divergence is the first tick past the shorter terminal.
        first_div = Some((shorter + 1) as u64);
    }
    let terminal_divergence = a.terminal_code != b.terminal_code;
    let mut bytes = Vec::new();
    bytes.extend_from_slice(mode.as_bytes());
    for r in [&a, &b] {
        bytes.extend_from_slice(r.run_intent_id.as_bytes());
        bytes.extend_from_slice(&r.terminal_tick.to_le_bytes());
        bytes.push(r.terminal_code);
        bytes.extend_from_slice(r.final_digest.as_bytes());
    }
    bytes.extend_from_slice(&prefix.to_le_bytes());
    bytes.extend_from_slice(&first_div.unwrap_or(0).to_le_bytes());
    bytes.push(u8::from(terminal_divergence));
    let receipt_digest = hash_domain(AB_SCHEMA, &bytes).to_hex();
    Ok(AbComparisonReceipt {
        schema: AB_SCHEMA,
        mode,
        a,
        b,
        common_prefix_ticks: prefix,
        first_divergence_tick: first_div,
        terminal_divergence,
        receipt_digest,
    })
}
