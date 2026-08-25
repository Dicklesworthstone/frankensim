//! Deterministic simulation/software reports and honest comparison placeholders
//! (bead `frankensim-euler-disc-emergent-flagship-t6314.8.7`).

use core::fmt;

/// Schema identifier for the Euler disc simulation report.
pub const EULER_SIMULATION_REPORT_SCHEMA_V1: &str =
    "org.frankensim.euler-disc.simulation-report.v1";

/// Physical comparison section authority state.
#[derive(Debug, Clone, PartialEq)]
pub enum PhysicalComparisonSection {
    /// No physical data is bound; typed placeholder with explicit next evidence requirements.
    NoData {
        /// Explanation of what physical evidence is required.
        required_next_evidence: &'static str,
    },
    /// Separately admitted physical comparison with content-bound scorecard and disposition.
    AdmittedPhysicalScorecard {
        /// Retained scorecard digest.
        scorecard_digest: String,
        /// Admitted disposition authority identifier.
        disposition_authority: String,
        /// Observed versus predicted spin duration difference.
        spin_time_relative_error: f64,
    },
}

impl PhysicalComparisonSection {
    /// Default honest placeholder for simulation-only reports.
    #[must_use]
    pub const fn no_data() -> Self {
        Self::NoData {
            required_next_evidence: "Certified high-speed optical metrology and calibrated physical telemetry",
        }
    }
}

/// A structured, deterministic simulation report for an Euler disc campaign.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerSimulationReport {
    /// Report schema version.
    pub schema_version: &'static str,
    /// Campaign identifier.
    pub campaign_id: String,
    /// Specimen identifier.
    pub specimen_id: String,
    /// Numerical model fidelity description.
    pub model_ladder: &'static str,
    /// Duration of the simulated run in seconds.
    pub duration_s: f64,
    /// Initial inclination angle in degrees.
    pub initial_inclination_deg: f64,
    /// Final terminal inclination angle in degrees.
    pub terminal_inclination_deg: f64,
    /// Initial mechanical energy in Joules.
    pub initial_energy_j: f64,
    /// Final mechanical energy in Joules.
    pub final_energy_j: f64,
    /// Energy conservation / dissipation defect in Joules.
    pub energy_defect_j: f64,
    /// Physical comparison section state.
    pub physical_section: PhysicalComparisonSection,
    /// Explicit no-claim disclosure.
    pub no_claim_disclosure: &'static str,
}

impl EulerSimulationReport {
    /// Create a new simulation report with honest NO-DATA physical placeholders.
    #[must_use]
    pub fn new_simulation_only(
        campaign_id: impl Into<String>,
        specimen_id: impl Into<String>,
        duration_s: f64,
        initial_inclination_deg: f64,
        terminal_inclination_deg: f64,
        initial_energy_j: f64,
        final_energy_j: f64,
    ) -> Self {
        let energy_defect = (initial_energy_j - final_energy_j).abs();
        Self {
            schema_version: EULER_SIMULATION_REPORT_SCHEMA_V1,
            campaign_id: campaign_id.into(),
            specimen_id: specimen_id.into(),
            model_ladder: "Coupled reduced rigid-body + unilateral contact + 1-mode base + exterior air",
            duration_s,
            initial_inclination_deg,
            terminal_inclination_deg,
            initial_energy_j,
            final_energy_j,
            energy_defect_j: energy_defect,
            physical_section: PhysicalComparisonSection::no_data(),
            no_claim_disclosure: "Simulation report only; does not establish physical validation or experimental ground truth without independent metrology.",
        }
    }

    /// Render human-readable Markdown summary.
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut out = String::new();
        out.push_str("# Euler Disc Simulation Report\n\n");
        out.push_str(&format!("- **Campaign ID**: `{}`\n", self.campaign_id));
        out.push_str(&format!("- **Specimen ID**: `{}`\n", self.specimen_id));
        out.push_str(&format!("- **Model Ladder**: {}\n", self.model_ladder));
        out.push_str(&format!("- **Simulated Duration**: {:.3} s\n", self.duration_s));
        out.push_str(&format!(
            "- **Inclination Range**: {:.2}° -> {:.2}°\n",
            self.initial_inclination_deg, self.terminal_inclination_deg
        ));
        out.push_str(&format!(
            "- **Energy (Initial / Final / Defect)**: {:.4} J / {:.4} J / {:.6} J\n\n",
            self.initial_energy_j, self.final_energy_j, self.energy_defect_j
        ));

        out.push_str("## Physical Validation Status\n\n");
        match &self.physical_section {
            PhysicalComparisonSection::NoData {
                required_next_evidence,
            } => {
                out.push_str("> [!NOTE]\n");
                out.push_str("> **NO DATA**: No physical comparison dataset bound.\n");
                out.push_str(&format!("> Required next evidence: {}\n\n", required_next_evidence));
            }
            PhysicalComparisonSection::AdmittedPhysicalScorecard {
                scorecard_digest,
                disposition_authority,
                spin_time_relative_error,
            } => {
                out.push_str(&format!("- **Scorecard Digest**: `{}`\n", scorecard_digest));
                out.push_str(&format!("- **Authority**: `{}`\n", disposition_authority));
                out.push_str(&format!(
                    "- **Spin Time Relative Error**: {:.2}%\n\n",
                    spin_time_relative_error * 100.0
                ));
            }
        }

        out.push_str("## No-Claim Disclosure\n\n");
        out.push_str(self.no_claim_disclosure);
        out.push('\n');
        out
    }
}

impl fmt::Display for EulerSimulationReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render_markdown())
    }
}

impl fmt::Display for PhysicalComparisonSection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoData { required_next_evidence } => {
                write!(f, "NO DATA (Required: {})", required_next_evidence)
            }
            Self::AdmittedPhysicalScorecard {
                scorecard_digest,
                disposition_authority,
                spin_time_relative_error,
            } => {
                write!(
                    f,
                    "Admitted [digest: {}, authority: {}, error: {:.2}%]",
                    scorecard_digest,
                    disposition_authority,
                    spin_time_relative_error * 100.0
                )
            }
        }
    }
}
