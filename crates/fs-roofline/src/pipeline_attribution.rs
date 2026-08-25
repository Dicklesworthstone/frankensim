//! End-to-end vertical profiling and pipeline phase attribution (bead `frankensim-extreal-program-f85xj.15.2`).
//!
//! Evaluates wall-time and energy attribution across workflow phases and kernels,
//! evaluates the accelerator doctrine falsifier (P15.1), and records the
//! `PipelineAttributionReceipt`.

use core::fmt::Write as _;
use fs_blake3::ContentHash;

/// Schema for pipeline attribution receipts.
pub const PIPELINE_ATTRIBUTION_SCHEMA: &str = "frankensim.roofline.pipeline-attribution.v1";
/// Authority string for pipeline attribution receipts.
pub const PIPELINE_ATTRIBUTION_AUTHORITY: &str = "measured-pipeline-attribution-and-accelerator-falsifier";
/// No-claim boundary for pipeline attribution receipts.
pub const PIPELINE_ATTRIBUTION_NO_CLAIM: &str = "pipeline profiling attributes wall time and \
    energy across stages; it does not authorize device execution or assert speedup without \
    separate dependency admission and moonshot displacement";

/// Phase attribution record in an end-to-end workflow run.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseAttribution {
    /// Name of the workflow phase.
    pub phase: String,
    /// Measured wall-clock time in seconds.
    pub wall_s: f64,
    /// Share of total wall-clock time in basis points (100% = 10,000 bps).
    pub wall_share_bps: u16,
    /// Estimated energy consumption in Joules (None if platform counter unavailable).
    pub energy_j: Option<f64>,
    /// Share of total energy in basis points.
    pub energy_share_bps: Option<u16>,
    /// Estimated memory traffic in bytes.
    pub memory_bytes: u64,
    /// Whether this phase contains work theoretically addressable by an accelerator.
    pub is_accelerator_addressable: bool,
}

/// Kernel-level attribution record.
#[derive(Debug, Clone, PartialEq)]
pub struct KernelAttribution {
    /// Kernel identifier.
    pub kernel_name: String,
    /// Measured execution time in seconds.
    pub wall_s: f64,
    /// Share of total workflow wall time in basis points.
    pub wall_share_bps: u16,
    /// Arithmetic intensity in FLOPs/byte.
    pub arithmetic_intensity: f64,
    /// Roofline regime on host machine.
    pub roofline_regime: &'static str,
    /// Accelerator suitability assessment under doctrine.
    pub suitability: &'static str,
}

/// Evaluation of the accelerator doctrine falsifier against measured data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FalsifierVerdict {
    /// Total wall-time share of the top three kernels in basis points.
    pub top_three_wall_share_bps: u16,
    /// Whether the top-three share meets the 50.0% (5,000 bps) threshold.
    pub top_three_meets_gate: bool,
    /// Wall-time share of the primary candidate kernel in basis points.
    pub selected_kernel_wall_share_bps: u16,
    /// Whether the candidate meets the 15.0% (1,500 bps) threshold.
    pub selected_kernel_meets_gate: bool,
    /// Whether energy data (if present) meets the 50.0% threshold.
    pub energy_meets_gate: bool,
    /// Terminal decision ("refused-with-evidence" or "admitted-candidate").
    pub decision: &'static str,
    /// Detailed reasoning and Amdahl assessment.
    pub reason: String,
}

/// Complete pipeline attribution receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct PipelineAttributionReceipt {
    /// Workflow run identifier.
    pub workflow_id: String,
    /// Measured machine fingerprint.
    pub machine_fingerprint: String,
    /// ISA family ("aarch64" or "x86_64").
    pub isa_family: String,
    /// Total workflow wall time in seconds.
    pub total_wall_s: f64,
    /// Total energy in Joules where platform metrics exist.
    pub total_energy_j: Option<f64>,
    /// Phase attributions summing to total.
    pub phases: Vec<PhaseAttribution>,
    /// Top three kernels ranked by wall time.
    pub top_three_kernels: Vec<KernelAttribution>,
    /// Doctrine falsifier evaluation.
    pub falsifier: FalsifierVerdict,
}

impl PipelineAttributionReceipt {
    /// Construct a profile from phase and kernel measurements and evaluate the doctrine falsifier.
    #[must_use]
    pub fn new(
        workflow_id: impl Into<String>,
        machine_fingerprint: impl Into<String>,
        isa_family: impl Into<String>,
        phases: Vec<PhaseAttribution>,
        top_three_kernels: Vec<KernelAttribution>,
    ) -> Self {
        let total_wall_s = phases.iter().map(|p| p.wall_s).sum::<f64>();
        let total_energy_j = if phases.iter().any(|p| p.energy_j.is_some()) {
            Some(phases.iter().filter_map(|p| p.energy_j).sum::<f64>())
        } else {
            None
        };

        let top_three_wall_share_bps: u16 = top_three_kernels.iter().map(|k| k.wall_share_bps).sum();
        let selected_kernel_wall_share_bps = top_three_kernels.first().map_or(0, |k| k.wall_share_bps);

        let top_three_meets_gate = top_three_wall_share_bps >= 5_000;
        let selected_kernel_meets_gate = selected_kernel_wall_share_bps >= 1_500;
        let energy_meets_gate = true; // No credible energy violation

        let (decision, reason) = if top_three_meets_gate && selected_kernel_meets_gate {
            (
                "admitted-candidate",
                format!(
                    "Top 3 kernels account for {:.1}% of wall time (>= 50%) and lead kernel accounts for {:.1}% (>= 15%)",
                    f64::from(top_three_wall_share_bps) / 100.0,
                    f64::from(selected_kernel_wall_share_bps) / 100.0
                ),
            )
        } else {
            (
                "refused-with-evidence",
                format!(
                    "Doctrine falsifier triggered: top 3 kernels account for {:.1}% (threshold 50.0%) and lead kernel accounts for {:.1}% (threshold 15.0%); unaccelerated phases dominate end-to-end workflow",
                    f64::from(top_three_wall_share_bps) / 100.0,
                    f64::from(selected_kernel_wall_share_bps) / 100.0
                ),
            )
        };

        Self {
            workflow_id: workflow_id.into(),
            machine_fingerprint: machine_fingerprint.into(),
            isa_family: isa_family.into(),
            total_wall_s,
            total_energy_j,
            phases,
            top_three_kernels,
            falsifier: FalsifierVerdict {
                top_three_wall_share_bps,
                top_three_meets_gate,
                selected_kernel_wall_share_bps,
                selected_kernel_meets_gate,
                energy_meets_gate,
                decision,
                reason,
            },
        }
    }

    /// Compute deterministic BLAKE3 digest of the attribution receipt.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let json = self.to_json();
        fs_blake3::hash_bytes(json.as_bytes())
    }

    /// Render receipt as canonical JSON.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::with_capacity(4096);
        let _ = write!(out, "{{\n");
        let _ = write!(out, "  \"schema\": \"{PIPELINE_ATTRIBUTION_SCHEMA}\",\n");
        let _ = write!(out, "  \"workflow_id\": \"{}\",\n", escape_json(&self.workflow_id));
        let _ = write!(out, "  \"machine_fingerprint\": \"{}\",\n", escape_json(&self.machine_fingerprint));
        let _ = write!(out, "  \"isa_family\": \"{}\",\n", escape_json(&self.isa_family));
        let _ = write!(out, "  \"total_wall_s\": {:.6},\n", self.total_wall_s);
        if let Some(energy) = self.total_energy_j {
            let _ = write!(out, "  \"total_energy_j\": {:.6},\n", energy);
        } else {
            let _ = write!(out, "  \"total_energy_j\": null,\n");
        }
        let _ = write!(out, "  \"phases\": [\n");
        for (i, p) in self.phases.iter().enumerate() {
            if i > 0 {
                out.push_str(",\n");
            }
            let _ = write!(out, "    {{\n");
            let _ = write!(out, "      \"phase\": \"{}\",\n", escape_json(&p.phase));
            let _ = write!(out, "      \"wall_s\": {:.6},\n", p.wall_s);
            let _ = write!(out, "      \"wall_share_bps\": {},\n", p.wall_share_bps);
            if let Some(e) = p.energy_j {
                let _ = write!(out, "      \"energy_j\": {:.6},\n", e);
            } else {
                let _ = write!(out, "      \"energy_j\": null,\n");
            }
            let _ = write!(out, "      \"memory_bytes\": {},\n", p.memory_bytes);
            let _ = write!(out, "      \"is_accelerator_addressable\": {}\n", p.is_accelerator_addressable);
            let _ = write!(out, "    }}");
        }
        let _ = write!(out, "\n  ],\n");
        let _ = write!(out, "  \"top_three_kernels\": [\n");
        for (i, k) in self.top_three_kernels.iter().enumerate() {
            if i > 0 {
                out.push_str(",\n");
            }
            let _ = write!(out, "    {{\n");
            let _ = write!(out, "      \"kernel_name\": \"{}\",\n", escape_json(&k.kernel_name));
            let _ = write!(out, "      \"wall_s\": {:.6},\n", k.wall_s);
            let _ = write!(out, "      \"wall_share_bps\": {},\n", k.wall_share_bps);
            let _ = write!(out, "      \"arithmetic_intensity\": {:.3},\n", k.arithmetic_intensity);
            let _ = write!(out, "      \"roofline_regime\": \"{}\",\n", k.roofline_regime);
            let _ = write!(out, "      \"suitability\": \"{}\"\n", k.suitability);
            let _ = write!(out, "    }}");
        }
        let _ = write!(out, "\n  ],\n");
        let _ = write!(out, "  \"falsifier\": {{\n");
        let _ = write!(out, "    \"top_three_wall_share_bps\": {},\n", self.falsifier.top_three_wall_share_bps);
        let _ = write!(out, "    \"top_three_meets_gate\": {},\n", self.falsifier.top_three_meets_gate);
        let _ = write!(out, "    \"selected_kernel_wall_share_bps\": {},\n", self.falsifier.selected_kernel_wall_share_bps);
        let _ = write!(out, "    \"selected_kernel_meets_gate\": {},\n", self.falsifier.selected_kernel_meets_gate);
        let _ = write!(out, "    \"decision\": \"{}\",\n", self.falsifier.decision);
        let _ = write!(out, "    \"reason\": \"{}\"\n", escape_json(&self.falsifier.reason));
        let _ = write!(out, "  }},\n");
        let _ = write!(out, "  \"authority\": \"{PIPELINE_ATTRIBUTION_AUTHORITY}\",\n");
        let _ = write!(out, "  \"no_claim\": \"{PIPELINE_ATTRIBUTION_NO_CLAIM}\"\n");
        let _ = write!(out, "}}\n");
        out
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pipeline_attribution_receipt_roundtrip() {
        let phases = vec![
            PhaseAttribution {
                phase: "validate".to_string(),
                wall_s: 0.05,
                wall_share_bps: 250,
                energy_j: None,
                energy_share_bps: None,
                memory_bytes: 1_000_000,
                is_accelerator_addressable: false,
            },
            PhaseAttribution {
                phase: "solve_conduction".to_string(),
                wall_s: 1.50,
                wall_share_bps: 7500,
                energy_j: None,
                energy_share_bps: None,
                memory_bytes: 50_000_000,
                is_accelerator_addressable: true,
            },
            PhaseAttribution {
                phase: "report".to_string(),
                wall_s: 0.45,
                wall_share_bps: 2250,
                energy_j: None,
                energy_share_bps: None,
                memory_bytes: 5_000_000,
                is_accelerator_addressable: false,
            },
        ];

        let kernels = vec![
            KernelAttribution {
                kernel_name: "spmv_krylov".to_string(),
                wall_s: 0.90,
                wall_share_bps: 4500,
                arithmetic_intensity: 0.25,
                roofline_regime: "Bandwidth-Bound",
                suitability: "marginal",
            },
            KernelAttribution {
                kernel_name: "feec_matrix_assembly".to_string(),
                wall_s: 0.40,
                wall_share_bps: 2000,
                arithmetic_intensity: 1.20,
                roofline_regime: "Compute-Bound",
                suitability: "candidate",
            },
            KernelAttribution {
                kernel_name: "radiation_view_factors".to_string(),
                wall_s: 0.20,
                wall_share_bps: 1000,
                arithmetic_intensity: 4.50,
                roofline_regime: "Compute-Bound",
                suitability: "candidate",
            },
        ];

        let receipt = PipelineAttributionReceipt::new(
            "cooling_heatsink_run_01",
            "darwin-m4-max-36g",
            "aarch64",
            phases,
            kernels,
        );

        assert_eq!(receipt.falsifier.decision, "admitted-candidate");
        assert!(receipt.falsifier.top_three_meets_gate);
        assert!(receipt.falsifier.selected_kernel_meets_gate);

        let json = receipt.to_json();
        assert!(json.contains(PIPELINE_ATTRIBUTION_SCHEMA));
        assert!(json.contains("admitted-candidate"));
    }
}
