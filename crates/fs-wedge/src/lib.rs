//! fs-wedge — go-to-market wedge selection as data (plan addendum,
//! Proposal 7). Layer: UTIL (pure data + audit; no dependencies).
//!
//! The wedge is the beachhead. The load-bearing DOCTRINE is a NEGATIVE one:
//!
//! > DO NOT SELL AGAINST PEAK SINGLE-PHYSICS FIDELITY ANYWHERE.
//!
//! Unification at the solver level loses to specialized codes on every
//! individual physics; nobody buys a beautifully glued assembly of second-rate
//! solvers, and nobody needs FrankenSim where one mature code owns the whole
//! problem. The wedge must be a WORKFLOW that is today three tools, lossy
//! handoffs, and week-long iteration — where certified seams (the sheaf),
//! incremental re-solve of variants (Proposal 2), and autonomous gradient
//! exploration (Proposal 1) dominate EVEN WITH merely-decent kernels.
//!
//! This crate preserves the plan's original ranking for replay, but those
//! judgment-only scores are [`ScoreUse::SupersededForDecisionUse`]. Current
//! decisions use [`MeasuredWedgeInputs`]: source-inventory readiness,
//! validation-data access, CAD burden, and static compute envelopes, each with
//! a method and evidence pointer. An explicit nine-factor comparison then ranks
//! full CHT, SDF structural/topology assurance, and a narrower thermal-design-
//! assurance candidate. Ratings remain separate from evidence authority, and
//! both rating and weight flip sensitivities are recorded. The
//! [`CycleTimeBaseline`] is a per-step sourced envelope assembled under the
//! dossier protocol in `data/cycle-time-baseline-dossier.md`; the retired flat
//! placeholder survives only as [`RETIRED_PLACEHOLDER_BASELINE`], which the
//! kill-criterion evaluation refuses with a typed error.

/// The load-bearing negative doctrine of wedge selection.
pub const WEDGE_DOCTRINE: &str = "Do not sell against peak single-physics fidelity anywhere; the wedge is a \
     multi-tool workflow with lossy handoffs where certified seams + incremental \
     re-solve + autonomous gradients win even with merely-decent kernels.";

/// The four criteria a wedge vertical is scored on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WedgeCriterion {
    /// Kernels are individually MATURE AND MODEST (no peak-fidelity arms race);
    /// correlation-based bottom rungs make the fidelity ladder immediately real.
    KernelMaturity,
    /// The cross-team iteration loop is the ACKNOWLEDGED, quantified pain today.
    IterationPain,
    /// ROI is QUANTIFIABLE per design cycle.
    QuantifiableRoi,
    /// Regulatory friction is LOW (the evidence-package story matures on
    /// friendly ground before facing the FAA).
    LowRegulatoryFriction,
}

impl WedgeCriterion {
    /// All four criteria, in order.
    pub const ALL: [WedgeCriterion; 4] = [
        WedgeCriterion::KernelMaturity,
        WedgeCriterion::IterationPain,
        WedgeCriterion::QuantifiableRoi,
        WedgeCriterion::LowRegulatoryFriction,
    ];

    /// A stable slug.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            WedgeCriterion::KernelMaturity => "kernel-maturity",
            WedgeCriterion::IterationPain => "iteration-pain",
            WedgeCriterion::QuantifiableRoi => "quantifiable-roi",
            WedgeCriterion::LowRegulatoryFriction => "low-regulatory-friction",
        }
    }
}

/// A criterion score (`0..=10`) with its rationale.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriterionScore {
    /// Which criterion.
    pub criterion: WedgeCriterion,
    /// The score, `0..=10`.
    pub score: u8,
    /// Why.
    pub rationale: &'static str,
}

/// Whether a historical criterion score may drive a current wedge decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreUse {
    /// The plan score is retained for replay, but measured inputs supersede it.
    SupersededForDecisionUse,
}

impl ScoreUse {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            ScoreUse::SupersededForDecisionUse => "superseded-for-decision-use",
        }
    }

    /// Historical scores never authorize a current selection.
    #[must_use]
    pub const fn permits_decision(self) -> bool {
        false
    }
}

const fn s(criterion: WedgeCriterion, score: u8, rationale: &'static str) -> CriterionScore {
    CriterionScore {
        criterion,
        score,
        rationale,
    }
}

/// A candidate vertical with its rank, four-criteria scores, the proposals it
/// exercises, and a rationale.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Vertical {
    /// A stable slug.
    pub name: &'static str,
    /// A human name.
    pub display: &'static str,
    /// Historical plan rank (1 = proposed beachhead, then 2, 3).
    pub rank: u8,
    /// The historical four-criteria scores (in `WedgeCriterion::ALL` order).
    pub scores: [CriterionScore; 4],
    /// Whether those historical scores can drive a current decision.
    pub score_use: ScoreUse,
    /// The proposals this vertical progressively exercises.
    pub exercises: &'static [&'static str],
    /// Why this vertical, at this rank.
    pub rationale: &'static str,
}

impl Vertical {
    /// This vertical's historical score for a criterion.
    #[must_use]
    pub fn score(&self, criterion: WedgeCriterion) -> u8 {
        self.scores
            .iter()
            .find(|s| s.criterion == criterion)
            .map_or(0, |s| s.score)
    }

    /// The minimum score across all four criteria (a wedge is only as good as
    /// its weakest criterion).
    #[must_use]
    pub fn weakest_criterion_score(&self) -> u8 {
        self.scores.iter().map(|s| s.score).min().unwrap_or(0)
    }

    /// A current decision score, if this record has measured authority.
    ///
    /// The retained plan scores are deliberately never promoted by this API.
    #[must_use]
    pub const fn decision_score(&self, _criterion: WedgeCriterion) -> Option<u8> {
        match self.score_use {
            ScoreUse::SupersededForDecisionUse => None,
        }
    }
}

use WedgeCriterion::{IterationPain, KernelMaturity, LowRegulatoryFriction, QuantifiableRoi};

/// The ranked verticals: V1 conjugate heat transfer, then aeroelastic
/// screening, then additive-manufacturing distortion.
const VERTICALS: [Vertical; 3] = [
    Vertical {
        name: "conjugate-heat-transfer",
        display: "Conjugate heat transfer for electronics cooling",
        rank: 1,
        scores: [
            s(
                KernelMaturity,
                8,
                "conduction FEM + forced-convection CFD with correlation-based Nusselt rungs (the fs-ladder cht() bottom rung — makes Proposal 3 real)",
            ),
            s(
                IterationPain,
                9,
                "the thermal<->mechanical/layout iteration loop is the acknowledged pain: today 3 tools, lossy handoffs, week-long cycles",
            ),
            s(
                QuantifiableRoi,
                9,
                "ROI is quantifiable per design cycle (cycle-time reduction directly measurable)",
            ),
            s(
                LowRegulatoryFriction,
                9,
                "low regulatory friction — the evidence-package story matures on friendly ground before the FAA",
            ),
        ],
        score_use: ScoreUse::SupersededForDecisionUse,
        exercises: &["2", "1", "3", "12"],
        rationale: "the beachhead: modest mature kernels, acknowledged cross-team pain, quantifiable ROI, friendly regulatory ground",
    },
    Vertical {
        name: "aeroelastic-screening",
        display: "Aeroelastic screening",
        rank: 2,
        scores: [
            s(
                KernelMaturity,
                6,
                "structural + aerodynamic kernels are mature but the coupling is where handoffs hurt",
            ),
            s(
                IterationPain,
                8,
                "flutter/divergence screening iterates across structures and aero teams",
            ),
            s(
                QuantifiableRoi,
                7,
                "ROI via faster screening of the design envelope",
            ),
            s(
                LowRegulatoryFriction,
                5,
                "moderate friction — closer to certification-sensitive aerospace",
            ),
        ],
        score_use: ScoreUse::SupersededForDecisionUse,
        exercises: &["1"],
        rationale: "second vertical: progressively exercises Proposal 1 (autonomous gradient exploration across the coupled loop)",
    },
    Vertical {
        name: "additive-manufacturing-distortion",
        display: "Additive-manufacturing distortion",
        rank: 3,
        scores: [
            s(
                KernelMaturity,
                6,
                "thermo-mechanical distortion kernels exist but validation against builds is the pain",
            ),
            s(
                IterationPain,
                8,
                "print-measure-recompensate loops are slow and physical",
            ),
            s(
                QuantifiableRoi,
                7,
                "ROI via fewer scrapped builds / compensation iterations",
            ),
            s(
                LowRegulatoryFriction,
                6,
                "moderate friction depending on the end-use part",
            ),
        ],
        score_use: ScoreUse::SupersededForDecisionUse,
        exercises: &["11", "4"],
        rationale: "third vertical: exercises Proposal 11 (reality as another chart — registration against scans) and Proposal 4 (extend the complex into time)",
    },
];

/// One measured input axis for wedge selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputAxis {
    /// Executable kernel inventory, including explicitly missing seams.
    KernelReadiness,
    /// Public validation data, raw-data access, and reuse terms.
    ValidationDataAccess,
    /// Required geometry semantics compared with native `fs-io` admission.
    CadBurden,
    /// Static work envelope for one fidelity rung.
    ComputeCost,
}

impl InputAxis {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            InputAxis::KernelReadiness => "kernel-readiness",
            InputAxis::ValidationDataAccess => "validation-data-access",
            InputAxis::CadBurden => "cad-burden",
            InputAxis::ComputeCost => "compute-cost",
        }
    }
}

/// Inventory status of one measured input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Readiness {
    /// An executable or directly obtainable input exists for the stated scope.
    Present,
    /// A narrower component exists, but a required seam or authority is absent.
    Partial,
    /// No decision-usable implementation or data package was found.
    Absent,
}

impl Readiness {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Readiness::Present => "present",
            Readiness::Partial => "partial",
            Readiness::Absent => "absent",
        }
    }

    /// Highest readiness score this status is allowed to carry.
    #[must_use]
    pub const fn score_ceiling(self) -> u8 {
        match self {
            Readiness::Present => 10,
            Readiness::Partial => 7,
            Readiness::Absent => 2,
        }
    }
}

/// How a decision input was measured.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasurementMethod {
    /// Direct symbol/module inventory in the tracked Rust workspace.
    WorkspaceInventory,
    /// A crate contract's explicit no-claim boundary.
    ContractBoundaryReview,
    /// Review of an official publisher's data-access record.
    OfficialDatasetReview,
    /// Static operation-count or algorithmic-complexity analysis.
    StaticComplexityAnalysis,
    /// Review of a declared planning assumption whose empirical measurement is
    /// still pending. This method can record uncertainty; it cannot erase it.
    DecisionAssumptionReview,
}

impl MeasurementMethod {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            MeasurementMethod::WorkspaceInventory => "workspace-inventory",
            MeasurementMethod::ContractBoundaryReview => "contract-boundary-review",
            MeasurementMethod::OfficialDatasetReview => "official-dataset-review",
            MeasurementMethod::StaticComplexityAnalysis => "static-complexity-analysis",
            MeasurementMethod::DecisionAssumptionReview => "decision-assumption-review",
        }
    }
}

/// Kind of durable evidence pointer attached to a measured input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceKind {
    /// A tracked workspace-relative file plus a required text marker.
    WorkspacePath,
    /// A Beads issue identifier.
    Bead,
    /// A primary-source URL published by the data owner.
    OfficialSource,
}

impl EvidenceKind {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            EvidenceKind::WorkspacePath => "workspace-path",
            EvidenceKind::Bead => "bead",
            EvidenceKind::OfficialSource => "official-source",
        }
    }
}

/// Durable evidence for one measurement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidencePointer {
    /// Pointer kind.
    pub kind: EvidenceKind,
    /// Workspace-relative path, Bead ID, or official URL.
    pub reference: &'static str,
    /// Required source marker, Bead scope, or dataset locator.
    pub locator: &'static str,
}

impl EvidencePointer {
    /// Is the pointer structurally complete?
    #[must_use]
    pub fn is_complete(self) -> bool {
        !self.reference.trim().is_empty() && !self.locator.trim().is_empty()
    }
}

/// Common measured fields carried by every wedge input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Measurement {
    /// Measured inventory/access status.
    pub readiness: Readiness,
    /// Decision-readiness score, constrained by [`Readiness::score_ceiling`].
    pub score: u8,
    /// How the status was established.
    pub method: MeasurementMethod,
    /// One or more replay pointers.
    pub evidence: &'static [EvidencePointer],
    /// Concise measured result and its no-claim boundary.
    pub finding: &'static str,
}

impl Measurement {
    /// Does this measurement satisfy the structural evidence contract?
    #[must_use]
    pub fn is_complete(self) -> bool {
        self.score <= self.readiness.score_ceiling()
            && !self.finding.trim().is_empty()
            && !self.evidence.is_empty()
            && self.evidence.iter().all(|pointer| pointer.is_complete())
    }
}

/// One capability required by a candidate vertical.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KernelReadinessEntry {
    /// Stable capability label.
    pub capability: &'static str,
    /// Inventory result.
    pub measurement: Measurement,
}

/// Access assessment for one published validation source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValidationDataEntry {
    /// Dataset or benchmark name.
    pub dataset: &'static str,
    /// What raw data are directly obtainable.
    pub raw_data: &'static str,
    /// Explicit reuse/license terms found, or an explicit missing-terms note.
    pub license_terms: &'static str,
    /// Access assessment.
    pub measurement: Measurement,
}

/// Geometry burden compared with native `fs-io` admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CadBurdenEntry {
    /// Geometry semantics the candidate requires.
    pub required_geometry: &'static str,
    /// Geometry semantics admitted in the current workspace.
    pub admitted_geometry: &'static str,
    /// Decision-relevant missing semantics.
    pub missing_semantics: &'static str,
    /// Burden assessment.
    pub measurement: Measurement,
}

/// Static work envelope for one fidelity rung.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComputeCostEntry {
    /// Stable rung label.
    pub rung: &'static str,
    /// Variables in the work model.
    pub variables: &'static str,
    /// Static operation-count or complexity envelope; never a wall-time claim.
    pub work_envelope: &'static str,
    /// Availability and evidence for this envelope.
    pub measurement: Measurement,
}

/// Replayable, measured inputs for one wedge candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MeasuredWedgeInputs {
    /// Candidate vertical slug.
    pub vertical: &'static str,
    /// Inventory date; URLs and source paths are the replay handles.
    pub measured_on: &'static str,
    /// Required kernel capabilities.
    pub kernels: &'static [KernelReadinessEntry],
    /// Published validation-data access.
    pub validation_data: &'static [ValidationDataEntry],
    /// Required geometry compared with current native admission.
    pub cad_burden: &'static [CadBurdenEntry],
    /// Static compute envelopes by fidelity rung.
    pub compute_cost: &'static [ComputeCostEntry],
}

impl MeasuredWedgeInputs {
    /// Iterate over every common measurement in stable axis order.
    pub fn measurements(&self) -> impl Iterator<Item = &Measurement> {
        self.kernels
            .iter()
            .map(|entry| &entry.measurement)
            .chain(self.validation_data.iter().map(|entry| &entry.measurement))
            .chain(self.cad_burden.iter().map(|entry| &entry.measurement))
            .chain(self.compute_cost.iter().map(|entry| &entry.measurement))
    }

    /// Are all four axes populated and structurally evidence-complete?
    #[must_use]
    pub fn is_complete(&self) -> bool {
        !self.vertical.is_empty()
            && !self.measured_on.is_empty()
            && !self.kernels.is_empty()
            && !self.validation_data.is_empty()
            && !self.cad_burden.is_empty()
            && !self.compute_cost.is_empty()
            && self
                .measurements()
                .all(|measurement| measurement.is_complete())
            && self
                .kernels
                .iter()
                .all(|entry| !entry.capability.trim().is_empty())
            && self.validation_data.iter().all(|entry| {
                !entry.dataset.trim().is_empty()
                    && !entry.raw_data.trim().is_empty()
                    && !entry.license_terms.trim().is_empty()
            })
            && self.cad_burden.iter().all(|entry| {
                !entry.required_geometry.trim().is_empty()
                    && !entry.admitted_geometry.trim().is_empty()
                    && !entry.missing_semantics.trim().is_empty()
            })
            && self.compute_cost.iter().all(|entry| {
                !entry.rung.trim().is_empty()
                    && !entry.variables.trim().is_empty()
                    && !entry.work_envelope.trim().is_empty()
            })
    }
}

const fn evidence(
    kind: EvidenceKind,
    reference: &'static str,
    locator: &'static str,
) -> EvidencePointer {
    EvidencePointer {
        kind,
        reference,
        locator,
    }
}

const fn measured(
    readiness: Readiness,
    score: u8,
    method: MeasurementMethod,
    evidence: &'static [EvidencePointer],
    finding: &'static str,
) -> Measurement {
    Measurement {
        readiness,
        score,
        method,
        evidence,
        finding,
    }
}

const CHT_KERNELS: [KernelReadinessEntry; 6] = [
    KernelReadinessEntry {
        capability: "steady-conduction-fem",
        measurement: measured(
            Readiness::Present,
            8,
            MeasurementMethod::WorkspaceInventory,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-conduction/src/solve.rs",
                "pub struct ConductionSolver",
            )],
            "A tracked tetrahedral P1 steady-conduction solver now exists with anisotropic and temperature-dependent conductivity, typed boundary conditions, residual remeasurement, and energy-balance evidence; it is not a CHT coupling or flow solver.",
        ),
    },
    KernelReadinessEntry {
        capability: "thermal-natural-convection-lbm",
        measurement: measured(
            Readiness::Present,
            8,
            MeasurementMethod::WorkspaceInventory,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-lbm/src/thermal.rs",
                "pub struct ThermalLbm",
            )],
            "D2Q9 flow plus D2Q5 temperature with Boussinesq forcing and Nusselt reporting is executable, but it is a natural-convection slab rather than electronics forced flow.",
        ),
    },
    KernelReadinessEntry {
        capability: "forced-convection-correlations-and-fan-curve",
        measurement: measured(
            Readiness::Absent,
            1,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-ladder/src/lib.rs",
                "cheap bottom rung: forced-convection Nusselt correlation",
            )],
            "The CHT ladder names a correlation rung but executes only the generic Refine1d transfer; no Nusselt correlation catalog, pressure-drop curve, or fan operating-point solve is present.",
        ),
    },
    KernelReadinessEntry {
        capability: "time-dependent-heat-adjoint",
        measurement: measured(
            Readiness::Partial,
            6,
            MeasurementMethod::WorkspaceInventory,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-adjoint/src/timedep.rs",
                "pub struct HeatAdjoint",
            )],
            "A backward-Euler adjoint over caller-assembled mass and stiffness matrices exists; it does not assemble or differentiate a CHT model.",
        ),
    },
    KernelReadinessEntry {
        capability: "temperature-dependent-material-properties",
        measurement: measured(
            Readiness::Partial,
            6,
            MeasurementMethod::WorkspaceInventory,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-matdb/src/lib.rs",
                "conductivity(T)",
            )],
            "The material database can represent temperature-indexed conductivity claims with provenance; application-specific electronics materials and coolant coverage remain dataset dependent.",
        ),
    },
    KernelReadinessEntry {
        capability: "solid-fluid-thermal-coupling-and-contact-resistance",
        measurement: measured(
            Readiness::Absent,
            1,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-couple/CONTRACT.md",
                "The scalar evaluator does not execute vector/tensor/field",
            )],
            "Typed coupling metadata exists, but no field transfer closes solid conduction, fluid temperature, interface heat flux, or contact resistance.",
        ),
    },
];

const CHT_VALIDATION: [ValidationDataEntry; 1] = [ValidationDataEntry {
    dataset: "Sandia transient forced-to-natural convection vertical-plate benchmark",
    raw_data: "The publisher states that boundary-condition and system-response data are downloadable, including velocity profiles, wall heat flux, and wall shear stress.",
    license_terms: "The official record exposes the paper/data download but does not state an explicit dataset reuse license; electronics-package applicability is not established.",
    measurement: measured(
        Readiness::Partial,
        6,
        MeasurementMethod::OfficialDatasetReview,
        &[evidence(
            EvidenceKind::OfficialSource,
            "https://www.sandia.gov/research/publications/details/experimental-validation-benchmark-data-for-cfd-of-transient-convection-from-2016-06-23/",
            "SAND2016-4201J publisher record and data-download statement",
        )],
        "Raw benchmark quantities and uncertainty-qualified boundary conditions are identified, but reuse terms and direct electronics-CHT coverage are not pinned.",
    ),
}];

const CHT_CAD: [CadBurdenEntry; 1] = [CadBurdenEntry {
    required_geometry: "electronics assemblies, material regions, thin interfaces, internal flow passages, and declared units",
    admitted_geometry: "bounded strict triangular faceted STEP resource closure and estimated tessellation-to-SDF handoff",
    missing_semantics: "assemblies, product/material linkage, units/context, NURBS and general B-rep topology, and interface/thickness identity",
    measurement: measured(
        Readiness::Partial,
        3,
        MeasurementMethod::ContractBoundaryReview,
        &[evidence(
            EvidenceKind::WorkspacePath,
            "crates/fs-io/CONTRACT.md",
            "Full native STEP CAD semantics remain STAGED",
        )],
        "A faceted handoff can carry a prepared mesh, but native CAD admission cannot reconstruct the multi-material assembly semantics the vertical needs.",
    ),
}];

const CHT_COMPUTE: [ComputeCostEntry; 3] = [
    ComputeCostEntry {
        rung: "correlation-Nu",
        variables: "C = number of components or thermal interfaces",
        work_envelope: "Target envelope O(C), but no executable correlation/fan kernel exists, so neither operation count nor wall time is measured.",
        measurement: measured(
            Readiness::Absent,
            1,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-ladder/src/lib.rs",
                "correlation-Nu",
            )],
            "Only a rung label and advisory relative cost exist.",
        ),
    },
    ComputeCostEntry {
        rung: "thermal-lbm-slab",
        variables: "N = lattice cells; S = time steps; Qf = 9 flow populations; Qt = 5 temperature populations",
        work_envelope: "O(S*N*(Qf+Qt)) work and O(N*(Qf+Qt)) state for the implemented two-dimensional slab; no coupled solid mesh is included.",
        measurement: measured(
            Readiness::Present,
            7,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-lbm/src/thermal.rs",
                "One coupled step",
            )],
            "Loop structure gives a static linear-in-cells-per-step envelope, not a wall-time or electronics accuracy claim.",
        ),
    },
    ComputeCostEntry {
        rung: "coupled-RANS-or-LES",
        variables: "Nf = fluid degrees of freedom; Ns = solid degrees of freedom; I = nonlinear/coupling iterations; S = time steps",
        work_envelope: "Unmeasured: the declared RANS and LES rungs have no executable solver or transfer, so no defensible bound beyond symbolic variables is available.",
        measurement: measured(
            Readiness::Absent,
            0,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-ladder/CONTRACT.md",
                "The ladder does not run solves",
            )],
            "A relative-cost hint is not compute evidence.",
        ),
    },
];

const AERO_KERNELS: [KernelReadinessEntry; 4] = [
    KernelReadinessEntry {
        capability: "wing-shell-structure-and-modes",
        measurement: measured(
            Readiness::Partial,
            5,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-solid/CONTRACT.md",
                "3D and shells; higher-order families",
            )],
            "Two-dimensional solid, rod, and stability kernels exist, while three-dimensional shell elements needed for a wing model are explicitly staged.",
        ),
    },
    KernelReadinessEntry {
        capability: "unsteady-aerodynamic-loads",
        measurement: measured(
            Readiness::Partial,
            5,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-vpm/CONTRACT.md",
                "INVISCID 2-D core by DIRECT",
            )],
            "A deterministic two-dimensional inviscid vortex-particle core exists; three-dimensional filaments and the BEM+VPM airfoil credential are staged.",
        ),
    },
    KernelReadinessEntry {
        capability: "nonlinear-field-fsi",
        measurement: measured(
            Readiness::Absent,
            2,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-couple/CONTRACT.md",
                "The FSI fixture is the classic LINEARIZED",
            )],
            "The coupling crate demonstrates a scalar added-mass map and Aitken relaxation, not an interface transfer between aerodynamic and structural fields.",
        ),
    },
    KernelReadinessEntry {
        capability: "coupled-flutter-gradient",
        measurement: measured(
            Readiness::Absent,
            0,
            MeasurementMethod::WorkspaceInventory,
            &[evidence(
                EvidenceKind::Bead,
                "frankensim-extreal-program-f85xj.1.1",
                "inventory found no executable coupled flutter objective or adjoint chain",
            )],
            "No admitted objective connects aerodynamic load, structural dynamics, flutter detection, and a verified coupled gradient.",
        ),
    },
];

const AERO_VALIDATION: [ValidationDataEntry; 1] = [ValidationDataEntry {
    dataset: "NASA/AGARD Wing 445.6 flutter benchmark",
    raw_data: "A public NASA report contains geometry descriptions and flutter plots/points; a later NASA assessment says many sets lack sufficient geometric and modal information.",
    license_terms: "NASA marks the report public; no separate raw numeric package and explicit dataset license are pinned by this record.",
    measurement: measured(
        Readiness::Partial,
        5,
        MeasurementMethod::OfficialDatasetReview,
        &[
            evidence(
                EvidenceKind::OfficialSource,
                "https://ntrs.nasa.gov/api/citations/19890009875/downloads/19890009875.pdf",
                "AGARD standard configuration report, Wing 445.6",
            ),
            evidence(
                EvidenceKind::OfficialSource,
                "https://c3.ndc.nasa.gov/dashlink/static/media/other/Aeroelasticity_Benchmark_Assessment_InterimReport.pdf",
                "NASA assessment documents limited geometry and modal information",
            ),
        ],
        "The benchmark is publicly inspectable but not a pinned, machine-readable, license-explicit raw corpus sufficient for end-to-end validation.",
    ),
}];

const AERO_CAD: [CadBurdenEntry; 1] = [CadBurdenEntry {
    required_geometry: "three-dimensional wing surfaces, shell midsurfaces/thickness, material axes, control surfaces, and modal correspondence",
    admitted_geometry: "strict triangular faceted STEP subset without product/context semantics",
    missing_semantics: "NURBS and general B-rep, shell thickness/material axes, assemblies/control-surface joints, and modal mesh correspondence",
    measurement: measured(
        Readiness::Absent,
        2,
        MeasurementMethod::ContractBoundaryReview,
        &[evidence(
            EvidenceKind::WorkspacePath,
            "crates/fs-io/CONTRACT.md",
            "does not fit NURBS",
        )],
        "The native import subset does not preserve the structural/aerodynamic geometry semantics required to construct a flutter model.",
    ),
}];

const AERO_COMPUTE: [ComputeCostEntry; 3] = [
    ComputeCostEntry {
        rung: "structural-mode-screen",
        variables: "Nd = free structural degrees of freedom; K = requested modes; I = eigensolver iterations",
        work_envelope: "Problem-dependent sparse operator work approximately O(I*K*operator(Nd)); the current dense reduction is fixture-gated at Nd <= 4096.",
        measurement: measured(
            Readiness::Partial,
            5,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-solid/CONTRACT.md",
                "dense reduction fixture-gated at n ≤ 4096",
            )],
            "A structural stability envelope exists for current 2-D fixtures, not production wing shells.",
        ),
    },
    ComputeCostEntry {
        rung: "direct-vpm-aerodynamic-screen",
        variables: "N = vortex particles; S = RK4 steps",
        work_envelope: "Exactly 4*S*N^2 attempted source-target contributions on the checked direct kernel, plus O(S*N) stage work.",
        measurement: measured(
            Readiness::Present,
            7,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-vpm/CONTRACT.md",
                "exactly `4 S N²`",
            )],
            "The exact work receipt applies to the two-dimensional inviscid core only.",
        ),
    },
    ComputeCostEntry {
        rung: "coupled-flutter-boundary",
        variables: "Na = aerodynamic state; Ns = structural state; I = coupling iterations; F = frequency or flight-condition samples",
        work_envelope: "Unmeasured because no executable aero-structure transfer, coupled residual, or flutter-boundary driver exists.",
        measurement: measured(
            Readiness::Absent,
            0,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-couple/CONTRACT.md",
                "nonlinear FSI solve over real fluid/structure subsystems",
            )],
            "Component complexity cannot be promoted into a coupled cost claim.",
        ),
    },
];

const AM_KERNELS: [KernelReadinessEntry; 5] = [
    KernelReadinessEntry {
        capability: "moving-heat-source-and-phase-change",
        measurement: measured(
            Readiness::Absent,
            0,
            MeasurementMethod::WorkspaceInventory,
            &[evidence(
                EvidenceKind::Bead,
                "frankensim-extreal-program-f85xj.1.1",
                "inventory found no laser path, melt-pool, phase-change, or powder-bed kernel",
            )],
            "No production kernel represents the AM process heat source, melt pool, or phase evolution.",
        ),
    },
    KernelReadinessEntry {
        capability: "three-dimensional-inelastic-distortion",
        measurement: measured(
            Readiness::Absent,
            1,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-solid/CONTRACT.md",
                "J2 continuum element wiring",
            )],
            "The solid crate explicitly stages three-dimensional elements and J2 continuum wiring required for residual-stress distortion.",
        ),
    },
    KernelReadinessEntry {
        capability: "layer-activation-time-sequencing",
        measurement: measured(
            Readiness::Partial,
            3,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-time/CONTRACT.md",
                "`slabs` module",
            )],
            "Generic feature-gated time slabs and activation reporting exist, but no AM layer birth/death or process-state adapter is implemented.",
        ),
    },
    KernelReadinessEntry {
        capability: "manufacturing-constraint-screen",
        measurement: measured(
            Readiness::Partial,
            3,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-fab/CONTRACT.md",
                "evaluators (overhang via surface-normal fields",
            )],
            "Scalar overhang and minimum-feature constraints exist; geometry-derived additive checks are explicitly staged.",
        ),
    },
    KernelReadinessEntry {
        capability: "as-built-registration",
        measurement: measured(
            Readiness::Partial,
            4,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-asbuilt/CONTRACT.md",
                "v1 is 2-D rigid registration",
            )],
            "Two-dimensional known-correspondence registration exists; three-dimensional Kabsch/ICP and full CT or point-cloud admission are staged.",
        ),
    },
];

const AM_VALIDATION: [ValidationDataEntry; 1] = [ValidationDataEntry {
    dataset: "NIST Additive Manufacturing Benchmark Test Series (AM Bench)",
    raw_data: "Public measurement data and metadata are stored in the NIST Public Data Repository; datasets include thermography, residual strain/stress, and part deflection, with some datasets larger than 1 TB and mirrored to SciServer.",
    license_terms: "NIST directs users to the dataset DOI and Fair Use citation guidance; this inventory does not promote that guidance into a software-style license or cover unpublished challenge keys.",
    measurement: measured(
        Readiness::Partial,
        7,
        MeasurementMethod::OfficialDatasetReview,
        &[
            evidence(
                EvidenceKind::OfficialSource,
                "https://www.nist.gov/ambench/am-bench-data-and-challenge-problems-0",
                "direct links to AM Bench measurement data",
            ),
            evidence(
                EvidenceKind::OfficialSource,
                "https://www.nist.gov/ambench/am-bench-data-management-systems",
                "PDR raw-data access, DOI citation guidance, and SciServer size/access notes",
            ),
        ],
        "This is the strongest data-access candidate, but a specific distortion case, version, files, checksum, and dataset-specific reuse terms still must be pinned before validation execution.",
    ),
}];

const AM_CAD: [CadBurdenEntry; 1] = [CadBurdenEntry {
    required_geometry: "build orientation, supports, scan regions, powder layers, material/process zones, and pre/post-build correspondence",
    admitted_geometry: "strict faceted STEP import and write-only minimal 3MF",
    missing_semantics: "3MF import, build/support/process metadata, assembly/material linkage, and three-dimensional scan correspondence",
    measurement: measured(
        Readiness::Partial,
        3,
        MeasurementMethod::ContractBoundaryReview,
        &[evidence(
            EvidenceKind::WorkspacePath,
            "crates/fs-io/CONTRACT.md",
            "3MF/GLB are WRITE-ONLY",
        )],
        "Prepared faceted meshes can enter, but native ingestion does not preserve the process and as-built semantics needed by a distortion workflow.",
    ),
}];

const AM_COMPUTE: [ComputeCostEntry; 3] = [
    ComputeCostEntry {
        rung: "process-thermal",
        variables: "Ne = active thermal elements; L = layers; St = thermal steps per layer; In = nonlinear iterations",
        work_envelope: "Unmeasured because the moving heat source, phase change, activation, and thermal process solver are absent.",
        measurement: measured(
            Readiness::Absent,
            0,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::Bead,
                "frankensim-extreal-program-f85xj.1.1",
                "no executable AM process rung found in tracked inventory",
            )],
            "Symbolic variables are retained without fabricating an operation or wall-time estimate.",
        ),
    },
    ComputeCostEntry {
        rung: "thermo-mechanical-distortion",
        variables: "Nt = thermal degrees of freedom; Nu = displacement degrees of freedom; L = layers; I = nonlinear/coupling iterations",
        work_envelope: "Unmeasured because no three-dimensional inelastic element, activation adapter, or coupled distortion driver exists.",
        measurement: measured(
            Readiness::Absent,
            0,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-solid/CONTRACT.md",
                "contact-law coefficient wiring. Plasticity flow remains successor",
            )],
            "No single-physics structural solve is used as a proxy for the missing coupled rung.",
        ),
    },
    ComputeCostEntry {
        rung: "as-built-registration-screen",
        variables: "N = known fiducial correspondences",
        work_envelope: "Exactly 6*N point visits for the current two-dimensional registration preflight/solve path; full 3-D scan registration is outside the envelope.",
        measurement: measured(
            Readiness::Present,
            7,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-asbuilt/CONTRACT.md",
                "exactly `6n` point visits",
            )],
            "The exact count covers a narrow inspection rung, not process simulation or three-dimensional compensation.",
        ),
    },
];

const MEASURED_INPUTS: [MeasuredWedgeInputs; 3] = [
    MeasuredWedgeInputs {
        vertical: "conjugate-heat-transfer",
        measured_on: "2026-07-22",
        kernels: &CHT_KERNELS,
        validation_data: &CHT_VALIDATION,
        cad_burden: &CHT_CAD,
        compute_cost: &CHT_COMPUTE,
    },
    MeasuredWedgeInputs {
        vertical: "aeroelastic-screening",
        measured_on: "2026-07-22",
        kernels: &AERO_KERNELS,
        validation_data: &AERO_VALIDATION,
        cad_burden: &AERO_CAD,
        compute_cost: &AERO_COMPUTE,
    },
    MeasuredWedgeInputs {
        vertical: "additive-manufacturing-distortion",
        measured_on: "2026-07-22",
        kernels: &AM_KERNELS,
        validation_data: &AM_VALIDATION,
        cad_burden: &AM_CAD,
        compute_cost: &AM_COMPUTE,
    },
];

/// Measured decision inputs for all candidate verticals.
#[must_use]
pub fn measured_wedge_inputs() -> &'static [MeasuredWedgeInputs] {
    &MEASURED_INPUTS
}

/// Measured decision inputs for one candidate vertical.
#[must_use]
pub fn measured_inputs_for(vertical: &str) -> Option<&'static MeasuredWedgeInputs> {
    MEASURED_INPUTS
        .iter()
        .find(|inputs| inputs.vertical == vertical)
}

/// One factor in the explicit wedge-comparison function.
///
/// Every factor is oriented so that a larger rating is better. In particular,
/// `CadBurden`, `ComputeCost`, and `RegulatoryRisk` rate LOW burden/cost/risk.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ScoringFactor {
    /// Evidence that the workflow is costly or slow for a real customer.
    CustomerPain,
    /// Readiness of the required executable kernels and their seams.
    KernelReadiness,
    /// Practicality of constructing a discriminating validation program.
    ValidationTractability,
    /// Access to raw, reusable, uncertainty-bearing validation data.
    DataAccess,
    /// LOW native-geometry admission and preparation burden.
    CadBurden,
    /// SHORT path to a defensible customer decision.
    TimeToDecision,
    /// Product differentiation from certified composition and iteration.
    Differentiation,
    /// LOW compute burden at decision-useful fidelity.
    ComputeCost,
    /// LOW regulatory and certification friction for the beachhead.
    RegulatoryRisk,
}

impl ScoringFactor {
    /// Canonical factor order used by scoring, reports, and tie-breaking.
    pub const ALL: [ScoringFactor; 9] = [
        ScoringFactor::CustomerPain,
        ScoringFactor::KernelReadiness,
        ScoringFactor::ValidationTractability,
        ScoringFactor::DataAccess,
        ScoringFactor::CadBurden,
        ScoringFactor::TimeToDecision,
        ScoringFactor::Differentiation,
        ScoringFactor::ComputeCost,
        ScoringFactor::RegulatoryRisk,
    ];

    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            ScoringFactor::CustomerPain => "customer-pain",
            ScoringFactor::KernelReadiness => "kernel-readiness",
            ScoringFactor::ValidationTractability => "validation-tractability",
            ScoringFactor::DataAccess => "data-access",
            ScoringFactor::CadBurden => "low-cad-burden",
            ScoringFactor::TimeToDecision => "short-time-to-decision",
            ScoringFactor::Differentiation => "differentiation",
            ScoringFactor::ComputeCost => "low-compute-cost",
            ScoringFactor::RegulatoryRisk => "low-regulatory-risk",
        }
    }

    const fn index(self) -> usize {
        match self {
            ScoringFactor::CustomerPain => 0,
            ScoringFactor::KernelReadiness => 1,
            ScoringFactor::ValidationTractability => 2,
            ScoringFactor::DataAccess => 3,
            ScoringFactor::CadBurden => 4,
            ScoringFactor::TimeToDecision => 5,
            ScoringFactor::Differentiation => 6,
            ScoringFactor::ComputeCost => 7,
            ScoringFactor::RegulatoryRisk => 8,
        }
    }
}

/// One normalized weight in the wedge-comparison function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FactorWeight {
    /// Factor being weighted.
    pub factor: ScoringFactor,
    /// Integer percentage points. All nine weights must sum to 100.
    pub weight: u8,
}

/// Recorded default weights. Their sum is exactly 100.
pub const DEFAULT_FACTOR_WEIGHTS: [FactorWeight; 9] = [
    FactorWeight {
        factor: ScoringFactor::CustomerPain,
        weight: 15,
    },
    FactorWeight {
        factor: ScoringFactor::KernelReadiness,
        weight: 18,
    },
    FactorWeight {
        factor: ScoringFactor::ValidationTractability,
        weight: 12,
    },
    FactorWeight {
        factor: ScoringFactor::DataAccess,
        weight: 10,
    },
    FactorWeight {
        factor: ScoringFactor::CadBurden,
        weight: 10,
    },
    FactorWeight {
        factor: ScoringFactor::TimeToDecision,
        weight: 12,
    },
    FactorWeight {
        factor: ScoringFactor::Differentiation,
        weight: 10,
    },
    FactorWeight {
        factor: ScoringFactor::ComputeCost,
        weight: 8,
    },
    FactorWeight {
        factor: ScoringFactor::RegulatoryRisk,
        weight: 5,
    },
];

/// A 0..=10 factor rating attached to the measurement that justifies it.
///
/// `rating` is the comparative decision input. `measurement.score` separately
/// records the readiness/authority of the evidence, so an attractive but
/// assumption-heavy factor cannot masquerade as a well-validated one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FactorRating {
    /// Which factor this input supplies.
    pub factor: ScoringFactor,
    /// Comparative value, 0 (worst) through 10 (best).
    pub rating: u8,
    /// Measurement method, authority, finding, and replay pointers.
    pub measurement: Measurement,
    /// Why the evidence maps to this comparative rating.
    pub rationale: &'static str,
}

impl FactorRating {
    /// Is this a structurally valid, evidence-bearing decision input?
    #[must_use]
    pub fn is_complete(self) -> bool {
        self.rating <= 10 && self.measurement.is_complete() && !self.rationale.trim().is_empty()
    }
}

/// One candidate in the measured comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComparisonCandidate {
    /// Stable candidate slug.
    pub name: &'static str,
    /// Human-readable candidate name.
    pub display: &'static str,
    /// Date on which the comparison inputs were reviewed.
    pub measured_on: &'static str,
    /// Git revision of the code inventory reviewed before this comparison.
    pub inventory_revision: &'static str,
    /// Exactly one measured input for every scoring factor.
    pub factors: &'static [FactorRating],
    /// Strongest case for this candidate if it finishes second.
    pub minority_case: &'static str,
}

/// Refusals from the pure scoring and sensitivity routines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringError {
    /// The supplied weights do not sum to 100.
    WeightsNotNormalized {
        /// Observed sum of the supplied integer weights.
        sum: u16,
    },
    /// A factor appears more than once in the weight vector.
    DuplicateWeight {
        /// Repeated factor.
        factor: ScoringFactor,
    },
    /// A factor is absent from the weight vector.
    MissingWeight {
        /// Missing factor.
        factor: ScoringFactor,
    },
    /// No candidates were supplied.
    NoCandidates,
    /// Candidate slugs must be unique.
    DuplicateCandidate {
        /// Repeated candidate slug.
        candidate: &'static str,
    },
    /// A candidate does not carry exactly one input for every factor.
    IncompleteCandidate {
        /// Incomplete candidate slug.
        candidate: &'static str,
    },
    /// A candidate repeats one factor.
    DuplicateFactor {
        /// Candidate containing the repeated factor.
        candidate: &'static str,
        /// Repeated factor.
        factor: ScoringFactor,
    },
    /// A factor rating, rationale, or measurement is incomplete.
    InvalidFactorInput {
        /// Candidate containing the invalid input.
        candidate: &'static str,
        /// Invalid factor.
        factor: ScoringFactor,
    },
    /// A requested weight tilt does not increase the named factor.
    InvalidWeightTilt {
        /// Factor whose weight was to increase.
        factor: ScoringFactor,
        /// Requested replacement weight.
        requested: u8,
    },
}

/// One scored row in deterministic rank order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateScore {
    /// Candidate slug.
    pub candidate: &'static str,
    /// Weighted sum in integer points, in the closed range 0..=1000.
    pub weighted_total: u16,
}

/// Minimum improvement to one challenger rating that would win the ranking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RatingFlipSensitivity {
    /// Candidate being tested against the baseline winner.
    pub challenger: &'static str,
    /// Single rating allowed to change.
    pub factor: ScoringFactor,
    /// Smallest winning rating, or `None` if 10 still cannot flip the result.
    pub required_rating: Option<u8>,
    /// Resulting challenger total at the required rating.
    pub resulting_total: Option<u16>,
}

/// Minimum upward tilt of one factor weight that would win the ranking.
///
/// The boosted factor gets the recorded weight; all other factors are scaled
/// proportionally back to a 100-point total with deterministic largest-
/// remainder allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WeightFlipSensitivity {
    /// Candidate being tested against the baseline winner.
    pub challenger: &'static str,
    /// Single weight allowed to increase.
    pub factor: ScoringFactor,
    /// Smallest winning factor weight, or `None` if even 100 cannot flip it.
    pub required_weight: Option<u8>,
    /// Resulting challenger total under the tilted weights.
    pub resulting_total: Option<u16>,
}

/// Complete ranked recommendation and its one-factor sensitivity tables.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RankedRecommendation {
    /// Candidate scores in descending order, then slug order for exact ties.
    pub ranked: Vec<CandidateScore>,
    /// Rank-1 candidate.
    pub recommended: &'static str,
    /// Rank-2 candidate.
    pub runner_up: &'static str,
    /// The runner-up's retained strongest case.
    pub minority_report: &'static str,
    /// Rating-change sensitivity for every challenger/factor pair.
    pub rating_sensitivities: Vec<RatingFlipSensitivity>,
    /// Weight-change sensitivity for every challenger/factor pair.
    pub weight_sensitivities: Vec<WeightFlipSensitivity>,
}

const fn factor_rating(
    factor: ScoringFactor,
    rating: u8,
    measurement: Measurement,
    rationale: &'static str,
) -> FactorRating {
    FactorRating {
        factor,
        rating,
        measurement,
        rationale,
    }
}

const CHT_FULL_FACTORS: [FactorRating; 9] = [
    factor_rating(
        ScoringFactor::CustomerPain,
        6,
        measured(
            Readiness::Absent,
            1,
            MeasurementMethod::DecisionAssumptionReview,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md",
                    "iteration loop between thermal and mechanical/layout teams is the acknowledged pain",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.1.3",
                    "measured cycle-time baseline protocol is still successor work",
                ),
            ],
            "Customer pain is declared in the plan, but no retained interview, workflow trace, or measured cycle-time baseline exists yet.",
        ),
        "The hypothesized multi-tool cooling loop is plausible and specific, but the rating is deliberately middling until the baseline protocol produces observations.",
    ),
    factor_rating(
        ScoringFactor::KernelReadiness,
        5,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::WorkspaceInventory,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-conduction/src/solve.rs",
                    "pub struct ConductionSolver",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-ladder/CONTRACT.md",
                    "The ladder does not run solves",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.8",
                    "low-Re forced-convection RANS remains gated on ratification",
                ),
            ],
            "Steady anisotropic conduction is real, while correlation, fan-network, solid-fluid transfer, and RANS/LES authority remain incomplete.",
        ),
        "One load-bearing thermal kernel now exists, but full electronics CHT still depends on several absent seams and its most expensive flow rung.",
    ),
    factor_rating(
        ScoringFactor::ValidationTractability,
        5,
        measured(
            Readiness::Partial,
            5,
            MeasurementMethod::OfficialDatasetReview,
            &[
                evidence(
                    EvidenceKind::OfficialSource,
                    "https://www.sandia.gov/research/publications/details/experimental-validation-benchmark-data-for-cfd-of-transient-convection-from-2016-06-23/",
                    "SAND2016-4201J boundary-condition and response data",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.4.4",
                    "electronics-thermal Level-C acquisition is not complete",
                ),
            ],
            "A convection benchmark is obtainable, but a complete electronics-package CHT validation chain and retained uncertainties are not pinned.",
        ),
        "The physics admits staged tests, but coupled geometry, material, interface, and airflow effects make the end-to-end validation program substantial.",
    ),
    factor_rating(
        ScoringFactor::DataAccess,
        6,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::OfficialDatasetReview,
            &[evidence(
                EvidenceKind::OfficialSource,
                "https://www.sandia.gov/research/publications/details/experimental-validation-benchmark-data-for-cfd-of-transient-convection-from-2016-06-23/",
                "publisher states boundary-condition and response data are downloadable",
            )],
            "Raw benchmark quantities are identified, but explicit reuse terms and direct electronics-package coverage remain unpinned.",
        ),
        "Public convection data prevent a zero, but the decision still lacks a versioned electronics-cooling corpus with explicit reuse terms.",
    ),
    factor_rating(
        ScoringFactor::CadBurden,
        3,
        measured(
            Readiness::Partial,
            3,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-io/CONTRACT.md",
                "broader CAD/EXPRESS",
            )],
            "Native admission remains a strict faceted STEP subset without assemblies, material linkage, general B-rep/NURBS, or interface identity.",
        ),
        "A full enclosure CHT model is unusually sensitive to assembly, material-region, thin-interface, and passage semantics that native import does not retain.",
    ),
    factor_rating(
        ScoringFactor::TimeToDecision,
        3,
        measured(
            Readiness::Absent,
            2,
            MeasurementMethod::DecisionAssumptionReview,
            &[
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.7",
                    "real solid-conduction to airflow transfer remains successor work",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.8",
                    "validated RANS rung remains successor work",
                ),
            ],
            "No measured schedule exists; two decision-critical coupled-flow steps remain unimplemented.",
        ),
        "The full candidate carries the longest critical path because it requires both the low-fidelity thermal stack and a ratified higher-fidelity flow rung.",
    ),
    factor_rating(
        ScoringFactor::Differentiation,
        8,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::DecisionAssumptionReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md",
                "The differentiation is composition, not physics",
            )],
            "Certified cross-tool composition is a declared product thesis; no customer comparison has measured willingness to pay for it.",
        ),
        "CHT strongly exercises the composition thesis, although the rating remains a product hypothesis rather than market evidence.",
    ),
    factor_rating(
        ScoringFactor::ComputeCost,
        2,
        measured(
            Readiness::Absent,
            0,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-ladder/CONTRACT.md",
                "The ladder does not run solves",
            )],
            "The decisive coupled RANS/LES rung has no executable operator count, memory envelope, or wall-time evidence.",
        ),
        "The full candidate contains the largest unmeasured state and iteration envelope, so low compute cost cannot be credited prospectively.",
    ),
    factor_rating(
        ScoringFactor::RegulatoryRisk,
        8,
        measured(
            Readiness::Partial,
            5,
            MeasurementMethod::DecisionAssumptionReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md",
                "regulatory friction is low",
            )],
            "The plan declares low friction for electronics cooling; no segment-specific legal or certification study is retained.",
        ),
        "The beachhead can target design assurance rather than safety-of-flight authority, but the score is an explicit planning assumption.",
    ),
];

const SDF_STRUCTURAL_FACTORS: [FactorRating; 9] = [
    factor_rating(
        ScoringFactor::CustomerPain,
        4,
        measured(
            Readiness::Absent,
            1,
            MeasurementMethod::DecisionAssumptionReview,
            &[evidence(
                EvidenceKind::Bead,
                "frankensim-extreal-program-f85xj.1.2",
                "candidate is motivated by technical depth, not a retained customer-pain measurement",
            )],
            "No retained interview, workflow trace, or cycle-time baseline establishes demand for SDF structural/topology assurance.",
        ),
        "The candidate has a credible engineering workflow but weaker market evidence than the declared thermal family.",
    ),
    factor_rating(
        ScoringFactor::KernelReadiness,
        8,
        measured(
            Readiness::Partial,
            7,
            MeasurementMethod::WorkspaceInventory,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-rep-sdf/CONTRACT.md",
                    "Signed-distance-field charts",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-cutfem/CONTRACT.md",
                    "vector Q1, plane-strain elasticity frontend",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-topopt/CONTRACT.md",
                    "SIMP with the modern hygiene stack",
                ),
            ],
            "SDF charts, 2-D CutFEM elasticity, density topology optimization, and verified gradients exist, but 3-D octree elasticity and production benchmark envelopes do not.",
        ),
        "This candidate aligns best with current integrated technical depth; the partial authority prevents its 2-D fixtures from being scored as a complete product kernel.",
    ),
    factor_rating(
        ScoringFactor::ValidationTractability,
        6,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::ContractBoundaryReview,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-cutfem/CONTRACT.md",
                    "fixed three-level log-log Q1 convergence",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-topopt/CONTRACT.md",
                    "FULL-CHAIN sensitivity vs FD",
                ),
            ],
            "Manufactured convergence, adjoint checks, deterministic goldens, and optimization laws exist; an independent product-scale structural corpus does not.",
        ),
        "The internal falsification ladder is unusually mature, but external validation would still need a scoped benchmark and acceptance envelope.",
    ),
    factor_rating(
        ScoringFactor::DataAccess,
        4,
        measured(
            Readiness::Absent,
            2,
            MeasurementMethod::DecisionAssumptionReview,
            &[evidence(
                EvidenceKind::Bead,
                "frankensim-extreal-program-f85xj.1.2",
                "no versioned external structural/topology corpus is named by the comparison task",
            )],
            "The current evidence is generated internally; no raw, reusable, uncertainty-bearing customer or public corpus is pinned for this candidate.",
        ),
        "Canonical fixtures are easy to generate, but they do not substitute for accessible external validation data.",
    ),
    factor_rating(
        ScoringFactor::CadBurden,
        6,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::ContractBoundaryReview,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-rep-sdf/CONTRACT.md",
                    "narrow-band level sets",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-io/CONTRACT.md",
                    "strict native triangular faceted-resource",
                ),
            ],
            "Native SDF design paths avoid some remeshing burden, but production CAD import still loses broad B-rep, assembly, and material semantics.",
        ),
        "SDF-native optimization lowers internal geometry burden relative to CHT assemblies, without pretending the external CAD bridge is complete.",
    ),
    factor_rating(
        ScoringFactor::TimeToDecision,
        7,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::WorkspaceInventory,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-topopt/CONTRACT.md",
                    "whole runs replay bitwise",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-cutfem/CONTRACT.md",
                    "3D octree instantiation",
                ),
            ],
            "A deterministic 2-D decision loop exists now; expansion to production 3-D remains outside current authority.",
        ),
        "A narrow 2-D assurance product can be demonstrated quickly, though production relevance depends on resisting a premature 3-D claim.",
    ),
    factor_rating(
        ScoringFactor::Differentiation,
        9,
        measured(
            Readiness::Present,
            8,
            MeasurementMethod::ContractBoundaryReview,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md",
                    "differentiability a routing requirement",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-adjoint/CONTRACT.md",
                    "density_pullback",
                ),
            ],
            "The candidate directly combines SDF routing, CutFEM, topology optimization, and checked adjoint paths that are unusual as one typed workflow.",
        ),
        "This is the strongest architectural fit and the runner-up's central minority argument.",
    ),
    factor_rating(
        ScoringFactor::ComputeCost,
        6,
        measured(
            Readiness::Partial,
            5,
            MeasurementMethod::StaticComplexityAnalysis,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-topopt/CONTRACT.md",
                "Per-iteration wall times are DEBUG-build measurements",
            )],
            "Deterministic fixture-scale loops exist, but production performance and 3-D memory envelopes are not certified.",
        ),
        "Existing sparse/adjoint structure supports a moderate rating, capped because debug fixture timing is not a product cost model.",
    ),
    factor_rating(
        ScoringFactor::RegulatoryRisk,
        5,
        measured(
            Readiness::Partial,
            4,
            MeasurementMethod::DecisionAssumptionReview,
            &[evidence(
                EvidenceKind::Bead,
                "frankensim-extreal-program-f85xj.1.2",
                "structural assurance segment and authority boundary are not yet ratified",
            )],
            "Regulatory exposure depends on the selected structural segment; no segment-specific assessment is retained.",
        ),
        "Structural decisions can approach certification-sensitive territory, so the comparison does not inherit electronics cooling's declared low-friction score.",
    ),
];

const THERMAL_ASSURANCE_FACTORS: [FactorRating; 9] = [
    factor_rating(
        ScoringFactor::CustomerPain,
        6,
        measured(
            Readiness::Absent,
            1,
            MeasurementMethod::DecisionAssumptionReview,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md",
                    "iteration loop between thermal and mechanical/layout teams is the acknowledged pain",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.1.3",
                    "measured baseline is still pending",
                ),
            ],
            "The thermal workflow pain remains a declared hypothesis pending the measured baseline protocol.",
        ),
        "The compromise addresses the same hypothesized customer loop as full CHT, without taking extra credit for unmeasured pain.",
    ),
    factor_rating(
        ScoringFactor::KernelReadiness,
        6,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::WorkspaceInventory,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-conduction/src/solve.rs",
                    "pub struct ConductionSolver",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.3",
                    "contact resistance remains successor work",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.4",
                    "surface radiation remains successor work",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.5",
                    "fan curve and flow network remain successor work",
                ),
            ],
            "Steady conduction with anisotropy and nonlinear material data exists; contact, radiation, fan networks, and calibrated forced-convection correlations remain incomplete.",
        ),
        "The compromise starts from a real core and has a bounded missing-kernel list, while explicitly deferring RANS rather than claiming it.",
    ),
    factor_rating(
        ScoringFactor::ValidationTractability,
        7,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::OfficialDatasetReview,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-conduction/CONTRACT.md",
                    "Every test prints a JSON-lines verdict",
                ),
                evidence(
                    EvidenceKind::OfficialSource,
                    "https://www.sandia.gov/research/publications/details/experimental-validation-benchmark-data-for-cfd-of-transient-convection-from-2016-06-23/",
                    "uncertainty-qualified boundary conditions and response data",
                ),
            ],
            "Analytic/manufactured conduction checks and a public convection benchmark support staged validation, but a retained electronics Level-C corpus is pending.",
        ),
        "Separating conduction, interfaces, radiation, and calibrated correlations makes falsification more tractable than validating a coupled RANS stack at once.",
    ),
    factor_rating(
        ScoringFactor::DataAccess,
        6,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::OfficialDatasetReview,
            &[evidence(
                EvidenceKind::OfficialSource,
                "https://www.sandia.gov/research/publications/details/experimental-validation-benchmark-data-for-cfd-of-transient-convection-from-2016-06-23/",
                "publisher states boundary-condition and response data are downloadable",
            )],
            "Public thermal data are identifiable, but a versioned electronics corpus and explicit dataset reuse terms remain to be pinned.",
        ),
        "The compromise can use the same accessible low-fidelity thermal evidence without waiting for a RANS validation corpus.",
    ),
    factor_rating(
        ScoringFactor::CadBurden,
        4,
        measured(
            Readiness::Partial,
            3,
            MeasurementMethod::ContractBoundaryReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "crates/fs-io/CONTRACT.md",
                "broader CAD/EXPRESS",
            )],
            "Prepared faceted models can enter, but board stackups, assemblies, materials, thin interfaces, and units are not natively reconstructed.",
        ),
        "Deferring resolved airflow reduces passage-detail burden, but thermal regions and interfaces still require semantics beyond native import.",
    ),
    factor_rating(
        ScoringFactor::TimeToDecision,
        6,
        measured(
            Readiness::Partial,
            5,
            MeasurementMethod::DecisionAssumptionReview,
            &[
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.2",
                    "forced-convection correlation rung is in progress",
                ),
                evidence(
                    EvidenceKind::Bead,
                    "frankensim-extreal-program-f85xj.5.8",
                    "RANS is explicitly deferred by this candidate",
                ),
            ],
            "No measured delivery schedule exists, but the candidate removes the gated RANS rung from its critical path.",
        ),
        "A bounded analytic/correlation stack should reach decisions sooner than full CHT, while remaining honest that four supporting kernels are unfinished.",
    ),
    factor_rating(
        ScoringFactor::Differentiation,
        8,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::DecisionAssumptionReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md",
                "The differentiation is composition, not physics",
            )],
            "The product thesis is an evidence-preserving ladder across thermal models; market differentiation has not been independently measured.",
        ),
        "The candidate showcases the architecture's evidence and fidelity-ladder strengths without entering a single-physics RANS fidelity race.",
    ),
    factor_rating(
        ScoringFactor::ComputeCost,
        8,
        measured(
            Readiness::Partial,
            6,
            MeasurementMethod::StaticComplexityAnalysis,
            &[
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-conduction/CONTRACT.md",
                    "continuous P₁ Lagrange",
                ),
                evidence(
                    EvidenceKind::WorkspacePath,
                    "crates/fs-ladder/CONTRACT.md",
                    "The ladder does not run solves",
                ),
            ],
            "The implemented sparse conduction core and planned algebraic/correlation rungs avoid a resolved-flow state, but no wall-time envelope is yet retained.",
        ),
        "Deferring RANS gives this candidate the smallest expected decision-useful state among the thermal options; the rating is not a performance claim.",
    ),
    factor_rating(
        ScoringFactor::RegulatoryRisk,
        8,
        measured(
            Readiness::Partial,
            5,
            MeasurementMethod::DecisionAssumptionReview,
            &[evidence(
                EvidenceKind::WorkspacePath,
                "COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md",
                "regulatory friction is low",
            )],
            "The plan declares low friction for electronics cooling; no segment-specific legal or certification study is retained.",
        ),
        "The assurance framing supports internal design decisions and evidence packaging without initially claiming safety-of-flight authority.",
    ),
];

const COMPARISON_CANDIDATES: [ComparisonCandidate; 3] = [
    ComparisonCandidate {
        name: "full-electronics-cooling-cht",
        display: "Full electronics-cooling CHT",
        measured_on: "2026-07-22",
        inventory_revision: "e5c8061f4faed986b831b8978d0c8d1812e960fb",
        factors: &CHT_FULL_FACTORS,
        minority_case: "Full CHT most completely exercises the original multi-tool coupling thesis and could become the strongest long-run moat. Its present score is held down by missing flow/coupling kernels, unmeasured compute, and a longer validation path—not by evidence that the market need is false.",
    },
    ComparisonCandidate {
        name: "sdf-structural-topology-assurance",
        display: "SDF structural/topology optimization assurance",
        measured_on: "2026-07-22",
        inventory_revision: "e5c8061f4faed986b831b8978d0c8d1812e960fb",
        factors: &SDF_STRUCTURAL_FACTORS,
        minority_case: "SDF structural assurance is the lowest technical-risk route to a differentiated demo: SDF charts, CutFEM elasticity, topology optimization, deterministic replay, and checked adjoints already compose. It should win if customer-pain evidence appears or if kernel readiness/CAD burden receive materially more weight; its current weakness is external market and validation-data evidence, plus honest 2-D scope.",
    },
    ComparisonCandidate {
        name: "thermal-design-assurance",
        display: "Thermal design assurance",
        measured_on: "2026-07-22",
        inventory_revision: "e5c8061f4faed986b831b8978d0c8d1812e960fb",
        factors: &THERMAL_ASSURANCE_FACTORS,
        minority_case: "Thermal design assurance preserves the declared thermal market thesis while limiting the first product to conduction, interfaces, radiation, fan/correlation rungs, and evidence-led decisions. Its remaining risk is that the customer-pain baseline and several supporting kernels are still pending.",
    },
];

/// The three candidates in the explicit measured comparison.
#[must_use]
pub fn comparison_candidates() -> &'static [ComparisonCandidate] {
    &COMPARISON_CANDIDATES
}

fn canonical_weight_values(weights: &[FactorWeight]) -> Result<[u8; 9], ScoringError> {
    let mut values = [0_u8; 9];
    let mut seen = [false; 9];
    let mut sum = 0_u16;
    for weight in weights {
        let index = weight.factor.index();
        if seen[index] {
            return Err(ScoringError::DuplicateWeight {
                factor: weight.factor,
            });
        }
        seen[index] = true;
        values[index] = weight.weight;
        sum += u16::from(weight.weight);
    }
    if sum != 100 {
        return Err(ScoringError::WeightsNotNormalized { sum });
    }
    for factor in ScoringFactor::ALL {
        if !seen[factor.index()] {
            return Err(ScoringError::MissingWeight { factor });
        }
    }
    Ok(values)
}

fn validate_candidate(candidate: &ComparisonCandidate) -> Result<(), ScoringError> {
    if candidate.name.trim().is_empty()
        || candidate.display.trim().is_empty()
        || candidate.measured_on.trim().is_empty()
        || candidate.inventory_revision.trim().is_empty()
        || candidate.minority_case.trim().is_empty()
        || candidate.factors.len() != ScoringFactor::ALL.len()
    {
        return Err(ScoringError::IncompleteCandidate {
            candidate: candidate.name,
        });
    }
    let mut seen = [false; 9];
    for input in candidate.factors {
        let index = input.factor.index();
        if seen[index] {
            return Err(ScoringError::DuplicateFactor {
                candidate: candidate.name,
                factor: input.factor,
            });
        }
        seen[index] = true;
        if !input.is_complete() {
            return Err(ScoringError::InvalidFactorInput {
                candidate: candidate.name,
                factor: input.factor,
            });
        }
    }
    if seen.iter().all(|present| *present) {
        Ok(())
    } else {
        Err(ScoringError::IncompleteCandidate {
            candidate: candidate.name,
        })
    }
}

/// Score and rank candidates with integer arithmetic.
///
/// Higher totals rank first. Exact ties are broken by ascending stable slug,
/// making the result independent of candidate input order.
pub fn score_candidates(
    weights: &[FactorWeight],
    candidates: &[ComparisonCandidate],
) -> Result<Vec<CandidateScore>, ScoringError> {
    let canonical = canonical_weight_values(weights)?;
    if candidates.is_empty() {
        return Err(ScoringError::NoCandidates);
    }
    for (index, candidate) in candidates.iter().enumerate() {
        validate_candidate(candidate)?;
        if candidates[..index]
            .iter()
            .any(|prior| prior.name == candidate.name)
        {
            return Err(ScoringError::DuplicateCandidate {
                candidate: candidate.name,
            });
        }
    }

    let mut scored = Vec::with_capacity(candidates.len());
    for candidate in candidates {
        let weighted_total = candidate.factors.iter().fold(0_u16, |total, input| {
            total + u16::from(input.rating) * u16::from(canonical[input.factor.index()])
        });
        scored.push(CandidateScore {
            candidate: candidate.name,
            weighted_total,
        });
    }
    scored.sort_by(|left, right| {
        right
            .weighted_total
            .cmp(&left.weighted_total)
            .then_with(|| left.candidate.cmp(right.candidate))
    });
    Ok(scored)
}

/// Increase one factor weight and deterministically renormalize all others.
pub fn tilted_weights(
    weights: &[FactorWeight],
    factor: ScoringFactor,
    requested: u8,
) -> Result<[FactorWeight; 9], ScoringError> {
    let original = canonical_weight_values(weights)?;
    let factor_index = factor.index();
    if requested <= original[factor_index] || requested > 100 {
        return Err(ScoringError::InvalidWeightTilt { factor, requested });
    }

    let old_other_total = 100_u16 - u16::from(original[factor_index]);
    let new_other_total = 100_u16 - u16::from(requested);
    let mut values = [0_u8; 9];
    let mut remainders = [0_u16; 9];
    let mut awarded = [false; 9];
    values[factor_index] = requested;
    let mut allocated = 0_u16;
    for other in ScoringFactor::ALL {
        let index = other.index();
        if index == factor_index {
            continue;
        }
        let numerator = u16::from(original[index]) * new_other_total;
        let quotient = numerator / old_other_total;
        values[index] = u8::try_from(quotient).expect("renormalized weight is at most 100");
        remainders[index] = numerator % old_other_total;
        allocated += quotient;
    }
    let mut leftover = new_other_total - allocated;
    while leftover > 0 {
        let mut best_index = None;
        for other in ScoringFactor::ALL {
            let index = other.index();
            if index == factor_index || awarded[index] {
                continue;
            }
            if best_index.is_none_or(|best| remainders[index] > remainders[best]) {
                best_index = Some(index);
            }
        }
        let index = best_index.expect("at most eight largest-remainder awards");
        values[index] += 1;
        awarded[index] = true;
        leftover -= 1;
    }

    Ok(ScoringFactor::ALL.map(|candidate_factor| FactorWeight {
        factor: candidate_factor,
        weight: values[candidate_factor.index()],
    }))
}

fn rating_sensitivities(
    weights: &[FactorWeight],
    candidates: &[ComparisonCandidate],
    baseline: &[CandidateScore],
) -> Result<Vec<RatingFlipSensitivity>, ScoringError> {
    let canonical = canonical_weight_values(weights)?;
    let winner = baseline[0];
    let mut rows = Vec::new();
    for challenger in baseline.iter().skip(1) {
        let candidate = candidates
            .iter()
            .find(|candidate| candidate.name == challenger.candidate)
            .expect("validated ranked candidate has source data");
        for factor in ScoringFactor::ALL {
            let input = candidate
                .factors
                .iter()
                .find(|input| input.factor == factor)
                .expect("validated candidate has every factor");
            let mut required_rating = None;
            let mut resulting_total = None;
            for rating in input.rating.saturating_add(1)..=10 {
                let total = challenger.weighted_total
                    + u16::from(rating - input.rating) * u16::from(canonical[factor.index()]);
                if total > winner.weighted_total
                    || (total == winner.weighted_total && challenger.candidate < winner.candidate)
                {
                    required_rating = Some(rating);
                    resulting_total = Some(total);
                    break;
                }
            }
            rows.push(RatingFlipSensitivity {
                challenger: challenger.candidate,
                factor,
                required_rating,
                resulting_total,
            });
        }
    }
    Ok(rows)
}

fn weight_sensitivities(
    weights: &[FactorWeight],
    candidates: &[ComparisonCandidate],
    baseline: &[CandidateScore],
) -> Result<Vec<WeightFlipSensitivity>, ScoringError> {
    let canonical = canonical_weight_values(weights)?;
    let winner = baseline[0].candidate;
    let mut rows = Vec::new();
    for challenger in baseline.iter().skip(1) {
        for factor in ScoringFactor::ALL {
            let mut required_weight = None;
            let mut resulting_total = None;
            for requested in canonical[factor.index()].saturating_add(1)..=100 {
                let tilted = tilted_weights(weights, factor, requested)?;
                let ranking = score_candidates(&tilted, candidates)?;
                if ranking[0].candidate == challenger.candidate {
                    required_weight = Some(requested);
                    resulting_total = ranking
                        .iter()
                        .find(|row| row.candidate == challenger.candidate)
                        .map(|row| row.weighted_total);
                    break;
                }
            }
            rows.push(WeightFlipSensitivity {
                challenger: challenger.candidate,
                factor,
                required_weight,
                resulting_total,
            });
        }
    }
    debug_assert_eq!(baseline[0].candidate, winner);
    Ok(rows)
}

/// Build the ranked recommendation, minority report, and complete sensitivity
/// tables for supplied normalized weights.
pub fn ranked_recommendation(
    weights: &[FactorWeight],
    candidates: &[ComparisonCandidate],
) -> Result<RankedRecommendation, ScoringError> {
    let ranked = score_candidates(weights, candidates)?;
    if ranked.len() < 2 {
        return Err(ScoringError::IncompleteCandidate {
            candidate: ranked[0].candidate,
        });
    }
    let recommended = ranked[0].candidate;
    let runner_up = ranked[1].candidate;
    let minority_report = candidates
        .iter()
        .find(|candidate| candidate.name == runner_up)
        .expect("validated runner-up has source data")
        .minority_case;
    let rating_sensitivities = rating_sensitivities(weights, candidates, &ranked)?;
    let weight_sensitivities = weight_sensitivities(weights, candidates, &ranked)?;
    Ok(RankedRecommendation {
        ranked,
        recommended,
        runner_up,
        minority_report,
        rating_sensitivities,
        weight_sensitivities,
    })
}

/// The default measured recommendation.
pub fn default_recommendation() -> Result<RankedRecommendation, ScoringError> {
    ranked_recommendation(&DEFAULT_FACTOR_WEIGHTS, comparison_candidates())
}

/// Render the full scoring table and both one-factor sensitivity tables as a
/// deterministic, verbose report artifact.
pub fn render_comparison_report() -> Result<String, ScoringError> {
    use core::fmt::Write as _;
    let recommendation = default_recommendation()?;
    let weights = canonical_weight_values(&DEFAULT_FACTOR_WEIGHTS)?;
    let mut out = String::from("FS-WEDGE-COMPARISON\tv1\n");
    out.push_str("WEIGHTS\tfactor\tweight\n");
    for factor in ScoringFactor::ALL {
        writeln!(
            out,
            "WEIGHT\t{}\t{}",
            factor.label(),
            weights[factor.index()]
        )
        .expect("write to String");
    }
    for candidate in comparison_candidates() {
        writeln!(
            out,
            "CANDIDATE\t{}\t{}\t{}\t{}",
            candidate.name, candidate.display, candidate.measured_on, candidate.inventory_revision
        )
        .expect("write to String");
        for input in candidate.factors {
            let contribution = u16::from(input.rating) * u16::from(weights[input.factor.index()]);
            writeln!(
                out,
                "FACTOR\t{}\t{}\trating={}\tweight={}\tcontribution={}\tauthority={}\tauthority_score={}\tmethod={}\trationale={}\tfinding={}",
                candidate.name,
                input.factor.label(),
                input.rating,
                weights[input.factor.index()],
                contribution,
                input.measurement.readiness.label(),
                input.measurement.score,
                input.measurement.method.label(),
                input.rationale,
                input.measurement.finding,
            )
            .expect("write to String");
            for pointer in input.measurement.evidence {
                writeln!(
                    out,
                    "EVIDENCE\t{}\t{}\t{}\t{}\t{}",
                    candidate.name,
                    input.factor.label(),
                    pointer.kind.label(),
                    pointer.reference,
                    pointer.locator,
                )
                .expect("write to String");
            }
        }
    }
    for (index, row) in recommendation.ranked.iter().enumerate() {
        writeln!(
            out,
            "RANK\t{}\t{}\t{}",
            index + 1,
            row.candidate,
            row.weighted_total
        )
        .expect("write to String");
    }
    writeln!(out, "RECOMMENDED\t{}", recommendation.recommended).expect("write to String");
    writeln!(out, "RUNNER_UP\t{}", recommendation.runner_up).expect("write to String");
    writeln!(out, "MINORITY_REPORT\t{}", recommendation.minority_report).expect("write to String");
    for row in &recommendation.rating_sensitivities {
        writeln!(
            out,
            "RATING_FLIP\t{}\t{}\t{}\t{}",
            row.challenger,
            row.factor.label(),
            row.required_rating
                .map_or("none".to_string(), |value| value.to_string()),
            row.resulting_total
                .map_or("none".to_string(), |value| value.to_string()),
        )
        .expect("write to String");
    }
    for row in &recommendation.weight_sensitivities {
        writeln!(
            out,
            "WEIGHT_FLIP\t{}\t{}\t{}\t{}",
            row.challenger,
            row.factor.label(),
            row.required_weight
                .map_or("none".to_string(), |value| value.to_string()),
            row.resulting_total
                .map_or("none".to_string(), |value| value.to_string()),
        )
        .expect("write to String");
    }
    Ok(out)
}

/// The ranked verticals.
#[must_use]
pub fn verticals() -> &'static [Vertical] {
    &VERTICALS
}

/// The four wedge-selection criteria.
#[must_use]
pub fn four_criteria() -> [WedgeCriterion; 4] {
    WedgeCriterion::ALL
}

/// The plan's historical rank-1 beachhead: conjugate heat transfer.
///
/// This accessor preserves the original proposal for replay; it does not
/// override [`ScoreUse::SupersededForDecisionUse`].
#[must_use]
pub fn chosen_wedge() -> &'static Vertical {
    VERTICALS
        .iter()
        .find(|v| v.rank == 1)
        .expect("a rank-1 wedge")
}

/// The six incumbent workflow steps of the measured cycle-time baseline
/// protocol, in fixed pipeline order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IncumbentStep {
    /// CAD import, cleanup, and defeaturing.
    CadPreparation,
    /// Grid/mesh generation for the thermal model.
    Meshing,
    /// Materials, boundary conditions, and power maps.
    SolverSetup,
    /// The solve itself (wall-clock; a running solve blocks the iteration).
    Solve,
    /// Result review and visualization.
    PostProcessing,
    /// Engineering report assembly.
    ReportAssembly,
}

impl IncumbentStep {
    /// Fixed pipeline order.
    pub const ALL: [IncumbentStep; 6] = [
        IncumbentStep::CadPreparation,
        IncumbentStep::Meshing,
        IncumbentStep::SolverSetup,
        IncumbentStep::Solve,
        IncumbentStep::PostProcessing,
        IncumbentStep::ReportAssembly,
    ];

    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            IncumbentStep::CadPreparation => "cad-preparation",
            IncumbentStep::Meshing => "meshing",
            IncumbentStep::SolverSetup => "solver-setup",
            IncumbentStep::Solve => "solve",
            IncumbentStep::PostProcessing => "post-processing",
            IncumbentStep::ReportAssembly => "report-assembly",
        }
    }
}

/// Authority class of one sourced duration figure, in descending authority.
/// Vendor material never becomes load-bearing on its own.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DurationSourceClass {
    /// An incumbent-workflow run executed and timed by this project.
    ExecutedIncumbentRun,
    /// A government-lab time-and-motion study.
    GovLabTimeStudy,
    /// A peer-reviewed publication reporting its own timings.
    PeerReviewedRun,
    /// A practitioner survey (self-reported; sponsorship flagged in caveat).
    PractitionerSurvey,
    /// A vendor or vendor-channel case study (marketing bias flagged).
    VendorCaseStudy,
}

impl DurationSourceClass {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            DurationSourceClass::ExecutedIncumbentRun => "executed-incumbent-run",
            DurationSourceClass::GovLabTimeStudy => "gov-lab-time-study",
            DurationSourceClass::PeerReviewedRun => "peer-reviewed-run",
            DurationSourceClass::PractitionerSurvey => "practitioner-survey",
            DurationSourceClass::VendorCaseStudy => "vendor-case-study",
        }
    }

    /// May a step bound rest on this class alone?
    #[must_use]
    pub const fn is_load_bearing(self) -> bool {
        !matches!(self, DurationSourceClass::VendorCaseStudy)
    }
}

/// One published figure backing a step bound, with its verbatim quote and
/// bias caveat. The figure is reproduced as published, never re-rounded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SourcedFigure {
    /// The figure exactly as published, with units.
    pub figure: &'static str,
    /// Verbatim sentence containing the figure.
    pub quote: &'static str,
    /// Author/organization, title, venue, year.
    pub citation: &'static str,
    /// URL verified reachable when the dossier was assembled.
    pub url: &'static str,
    /// Authority class.
    pub class: DurationSourceClass,
    /// Bias and applicability caveat.
    pub caveat: &'static str,
}

impl SourcedFigure {
    /// Is the provenance row structurally complete?
    #[must_use]
    pub fn is_complete(self) -> bool {
        !self.figure.trim().is_empty()
            && !self.quote.trim().is_empty()
            && !self.citation.trim().is_empty()
            && !self.url.trim().is_empty()
            && !self.caveat.trim().is_empty()
    }
}

/// The sourced duration envelope for one incumbent workflow step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaselineStepEstimate {
    /// Which step.
    pub step: IncumbentStep,
    /// Lower bound, engineer/wall hours (see the step's basis).
    pub low_hours: f64,
    /// Upper bound, hours.
    pub high_hours: f64,
    /// How the bounds were derived from the sources.
    pub basis: &'static str,
    /// The published figures behind the bounds.
    pub sources: &'static [SourcedFigure],
}

impl BaselineStepEstimate {
    /// Is the estimate structurally complete: finite positive bounds in
    /// order, a stated basis, and at least one complete load-bearing source?
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.low_hours.is_finite()
            && self.high_hours.is_finite()
            && self.low_hours > 0.0
            && self.high_hours >= self.low_hours
            && !self.basis.trim().is_empty()
            && !self.sources.is_empty()
            && self.sources.iter().all(|source| source.is_complete())
            && self
                .sources
                .iter()
                .any(|source| source.class.is_load_bearing())
    }
}

/// Where a baseline's numbers came from. Only non-placeholder records may
/// enter the kill-criterion decision path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaselineProvenance {
    /// A declared scoping figure with no protocol behind it. Refused.
    Placeholder,
    /// Assembled from published sources under the dossier protocol. This is
    /// an ESTIMATED quantity: nobody on this project timed the incumbent run.
    PublishedSourceDerived,
    /// An incumbent-workflow run executed and timed by this project.
    ExecutedRun,
}

impl BaselineProvenance {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            BaselineProvenance::Placeholder => "placeholder",
            BaselineProvenance::PublishedSourceDerived => "published-source-derived",
            BaselineProvenance::ExecutedRun => "executed-run",
        }
    }
}

/// The baseline that makes the cycle-time kill criterion measurable: a
/// per-step sourced envelope, never a point value. Uncertainty is carried in
/// the open, and the decision path refuses placeholders.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CycleTimeBaseline {
    /// Which vertical.
    pub vertical: &'static str,
    /// When the record was assembled (dossier date).
    pub assembled_on: &'static str,
    /// The representative incumbent task the envelope describes.
    pub representative_task: &'static str,
    /// The admissibility protocol the record was assembled under.
    pub protocol: &'static str,
    /// The provenance dossier backing every figure.
    pub dossier: EvidencePointer,
    /// Where the numbers came from.
    pub provenance: BaselineProvenance,
    /// Per-step sourced envelopes.
    pub steps: &'static [BaselineStepEstimate],
    /// Explicit hours-to-working-day conversion.
    pub hours_per_working_day: f64,
    /// The cycle-time reduction factor the kill criterion demands (`3.0`).
    pub target_reduction: f64,
    /// The window (quarters after GA) to hit it or re-select the wedge.
    pub kill_within_quarters: u8,
}

impl CycleTimeBaseline {
    /// Sum of per-step lower bounds, hours.
    #[must_use]
    pub fn low_hours(&self) -> f64 {
        self.steps.iter().map(|step| step.low_hours).sum()
    }

    /// Sum of per-step upper bounds, hours.
    #[must_use]
    pub fn high_hours(&self) -> f64 {
        self.steps.iter().map(|step| step.high_hours).sum()
    }

    /// Envelope lower bound in working days.
    #[must_use]
    pub fn baseline_days_low(&self) -> f64 {
        self.low_hours() / self.hours_per_working_day
    }

    /// Envelope upper bound in working days.
    #[must_use]
    pub fn baseline_days_high(&self) -> f64 {
        self.high_hours() / self.hours_per_working_day
    }

    /// Is the record structurally complete: one estimate per protocol step in
    /// pipeline order, every estimate complete, a dossier pointer, and a
    /// usable day conversion?
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.steps.len() == IncumbentStep::ALL.len()
            && self
                .steps
                .iter()
                .zip(IncumbentStep::ALL)
                .all(|(estimate, step)| estimate.step == step)
            && self.steps.iter().all(BaselineStepEstimate::is_complete)
            && !self.vertical.trim().is_empty()
            && !self.assembled_on.trim().is_empty()
            && !self.representative_task.trim().is_empty()
            && !self.protocol.trim().is_empty()
            && self.dossier.is_complete()
            && self.hours_per_working_day > 0.0
            && self.hours_per_working_day.is_finite()
    }

    /// Render the complete audit derivation: record identity, per-step
    /// breakdown, envelope, conversion, target, measured value, and the
    /// verdict rule.
    fn render_derivation(
        &self,
        measured_days: f64,
        reduction_low: f64,
        reduction_high: f64,
        verdict: KillVerdict,
    ) -> String {
        use core::fmt::Write as _;
        let mut derivation = String::new();
        writeln!(
            derivation,
            "kill-criterion evaluation for `{}` (baseline assembled {}, provenance {}, dossier {} # {})",
            self.vertical,
            self.assembled_on,
            self.provenance.label(),
            self.dossier.reference,
            self.dossier.locator,
        )
        .expect("write to String");
        writeln!(
            derivation,
            "  representative task: {}",
            self.representative_task
        )
        .expect("write to String");
        for step in self.steps {
            writeln!(
                derivation,
                "  step {:<16} {:>6.1}..{:>6.1} h  ({} sources) — {}",
                step.step.label(),
                step.low_hours,
                step.high_hours,
                step.sources.len(),
                step.basis,
            )
            .expect("write to String");
        }
        writeln!(
            derivation,
            "  baseline envelope: {:.1}..{:.1} h = {:.2}..{:.2} working days at {:.0} h/day",
            self.low_hours(),
            self.high_hours(),
            self.baseline_days_low(),
            self.baseline_days_high(),
            self.hours_per_working_day,
        )
        .expect("write to String");
        writeln!(
            derivation,
            "  measured FrankenSim cycle time: {measured_days:.3} days",
        )
        .expect("write to String");
        writeln!(
            derivation,
            "  implied reduction: {:.2}x..{:.2}x against target {:.1}x",
            reduction_low, reduction_high, self.target_reduction,
        )
        .expect("write to String");
        writeln!(
            derivation,
            "  verdict: {} (met only if the LOW bound clears the target; not-met only if the HIGH bound cannot)",
            verdict.label(),
        )
        .expect("write to String");
        derivation
    }

    /// Evaluate the `>= target_reduction x` kill criterion against a measured
    /// FrankenSim cycle time, refusing placeholders and non-measurable
    /// inputs, and returning the full derivation for audit.
    ///
    /// The verdict is conservative across the envelope: `Met` only when even
    /// the LOW baseline bound clears the target; `NotMet` when even the HIGH
    /// bound cannot; `Indeterminate` when the envelope straddles the target.
    pub fn evaluate_kill_criterion(
        &self,
        measured_days: f64,
    ) -> Result<KillCriterionEvaluation, KillCriterionError> {
        if self.provenance == BaselineProvenance::Placeholder {
            return Err(KillCriterionError::PlaceholderBaseline {
                vertical: self.vertical,
            });
        }
        if !self.is_complete() {
            return Err(KillCriterionError::IncompleteBaseline {
                vertical: self.vertical,
            });
        }
        if !measured_days.is_finite() || measured_days <= 0.0 {
            return Err(KillCriterionError::NonMeasurableCycleTime { measured_days });
        }

        let baseline_days_low = self.baseline_days_low();
        let baseline_days_high = self.baseline_days_high();
        let reduction_low = baseline_days_low / measured_days;
        let reduction_high = baseline_days_high / measured_days;
        let verdict = if reduction_low >= self.target_reduction {
            KillVerdict::Met
        } else if reduction_high < self.target_reduction {
            KillVerdict::NotMet
        } else {
            KillVerdict::Indeterminate
        };

        let derivation =
            self.render_derivation(measured_days, reduction_low, reduction_high, verdict);

        Ok(KillCriterionEvaluation {
            vertical: self.vertical,
            assembled_on: self.assembled_on,
            provenance: self.provenance,
            baseline_days_low,
            baseline_days_high,
            target_reduction: self.target_reduction,
            measured_days,
            reduction_low,
            reduction_high,
            verdict,
            derivation,
        })
    }
}

/// Typed refusals from the kill-criterion evaluation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KillCriterionError {
    /// The baseline is a declared placeholder; the decision path refuses it.
    PlaceholderBaseline {
        /// Vertical whose placeholder was offered.
        vertical: &'static str,
    },
    /// The baseline record is structurally incomplete.
    IncompleteBaseline {
        /// Vertical whose record is incomplete.
        vertical: &'static str,
    },
    /// The measured cycle time is not a positive finite number of days.
    NonMeasurableCycleTime {
        /// The refused value.
        measured_days: f64,
    },
}

/// Conservative three-way kill-criterion verdict across the envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KillVerdict {
    /// Even the LOW baseline bound clears the target reduction.
    Met,
    /// Even the HIGH baseline bound cannot clear the target reduction.
    NotMet,
    /// The envelope straddles the target; marketing may not claim "met".
    Indeterminate,
}

impl KillVerdict {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            KillVerdict::Met => "met",
            KillVerdict::NotMet => "not-met",
            KillVerdict::Indeterminate => "indeterminate",
        }
    }
}

/// One audited kill-criterion evaluation, logging which baseline record was
/// used and the full derivation.
#[derive(Debug, Clone, PartialEq)]
pub struct KillCriterionEvaluation {
    /// Vertical evaluated.
    pub vertical: &'static str,
    /// Baseline record date used.
    pub assembled_on: &'static str,
    /// Baseline provenance used (never `Placeholder`).
    pub provenance: BaselineProvenance,
    /// Envelope lower bound, working days.
    pub baseline_days_low: f64,
    /// Envelope upper bound, working days.
    pub baseline_days_high: f64,
    /// Required reduction factor.
    pub target_reduction: f64,
    /// The measured FrankenSim cycle time, days.
    pub measured_days: f64,
    /// `baseline_days_low / measured_days` — the claim's worst case.
    pub reduction_low: f64,
    /// `baseline_days_high / measured_days`.
    pub reduction_high: f64,
    /// Conservative verdict.
    pub verdict: KillVerdict,
    /// Full printed derivation: record identity, per-step breakdown,
    /// envelope, conversion, target, measured value, and verdict rule.
    pub derivation: String,
}

/// The retired scoping placeholder (`baseline_days = 5.0`, declared without a
/// protocol on 2026-07-07). Retained as data so the refusal path is proven
/// against the genuine legacy object; unreachable in the decision path.
pub const RETIRED_PLACEHOLDER_BASELINE: CycleTimeBaseline = CycleTimeBaseline {
    vertical: "conjugate-heat-transfer",
    assembled_on: "2026-07-07",
    representative_task: "an unspecified week-long loop (the retired flat 5.0-working-day figure)",
    protocol: "none: declared as a scoping placeholder without measurement or sources",
    dossier: EvidencePointer {
        kind: EvidenceKind::WorkspacePath,
        reference: "crates/fs-wedge/src/lib.rs",
        locator: "RETIRED_PLACEHOLDER_BASELINE",
    },
    provenance: BaselineProvenance::Placeholder,
    steps: &[],
    hours_per_working_day: 8.0,
    target_reduction: 3.0,
    kill_within_quarters: 2,
};

/// Dossier-backed measured cycle-time baseline for the cooling vertical:
/// one first-pass design iteration of a representative electronics-cooling
/// task (PCB + components + heat sink, forced convection / CHT) in a
/// commercial incumbent tool, dirty CAD to engineering report.
///
/// This record is ESTIMATED (published-source-derived): its authority ends at
/// the cited sources. It is NOT customer-pain evidence and must not upgrade
/// the comparison's customer-pain factor; an executed incumbent run
/// supersedes it under the dossier protocol.
pub const CHT_BASELINE: CycleTimeBaseline = CycleTimeBaseline {
    vertical: "conjugate-heat-transfer",
    assembled_on: "2026-07-22",
    representative_task: "thermal analysis of an electronics assembly (PCB + components + heat \
        sink, forced convection / conjugate heat transfer) taken from dirty CAD to an \
        engineering report in a commercial tool, one first-pass design iteration",
    protocol: "crates/fs-wedge/data/cycle-time-baseline-dossier.md section 1: six-step \
        decomposition; sources admitted in descending authority (executed run, gov-lab or \
        peer-reviewed time study, practitioner survey, vendor case with bias caveat); every \
        figure carries citation, verbatim quote, URL, and bias class; bounds are ranges, \
        never points; quarterly review alongside the ratification record",
    dossier: EvidencePointer {
        kind: EvidenceKind::WorkspacePath,
        reference: "crates/fs-wedge/data/cycle-time-baseline-dossier.md",
        locator: "fs-wedge-cycle-time-baseline-dossier-v1",
    },
    provenance: BaselineProvenance::PublishedSourceDerived,
    steps: &[
        BaselineStepEstimate {
            step: IncumbentStep::CadPreparation,
            low_hours: 1.0,
            high_hours: 16.0,
            basis: "floor: clean-CAD import plus minimal cleanup (survey fast quartile); \
                ceiling: dirty assembly at two working days, consistent with the 30%-over-one-day \
                model-build survey row while staying below defense-lab scale",
            sources: &[
                SourcedFigure {
                    figure: "39% of thermal engineers spend >1 h importing CAD; 24% <10 min",
                    quote: "39 percent spend over an hour importing CAD data, while 24 percent \
                        claimed they can run the same task in under 10 minutes.",
                    citation: "6SigmaET (Future Facilities), State of Thermal survey (n>170), \
                        via engineering.com, 2018",
                    url: "https://www.engineering.com/70-of-engineers-left-out-in-the-cold-by-thermal-simulations/",
                    class: DurationSourceClass::PractitionerSurvey,
                    caveat: "vendor-run survey; self-reported; sample skews to the vendor's \
                        contact list",
                },
                SourcedFigure {
                    figure: "30% spend >1 day building models; 37% <1 h",
                    quote: "only 37 percent of engineers spend less than an hour building their \
                        models, while 30 percent spend more than a day.",
                    citation: "6SigmaET State of Thermal survey via engineering.com, 2018",
                    url: "https://www.engineering.com/70-of-engineers-left-out-in-the-cold-by-thermal-simulations/",
                    class: DurationSourceClass::PractitionerSurvey,
                    caveat: "model build includes geometry idealization; same vendor-survey bias",
                },
                SourcedFigure {
                    figure: "preprocessing = 38% of total simulation time",
                    quote: "Preprocessing is the most time consuming part of the simulation \
                        process, taking up 38% of total simulation time.",
                    citation: "Tech-Clarity, Addressing the Bottlenecks of FEA Simulation \
                        (n>160), 2016",
                    url: "https://tech-clarity.com/simulation_bottlenecks/5467",
                    class: DurationSourceClass::PractitionerSurvey,
                    caveat: "all-FEA population, not thermal-specific; analyst research \
                        typically vendor-sponsored",
                },
                SourcedFigure {
                    figure: "defeaturing 48.8 h + decomposition 36.8 h initial engineer-hours \
                        (heat-transfer problems)",
                    quote: "Analysis Solid Model Creation/Edit 48.8; Geometry Decomposition 36.8 \
                        (Table 5, heat-transfer column, engineering hours)",
                    citation: "Hardwick, Clay, Boggs, Walsh, Larzelere, Altshuler, DART System \
                        Analysis, SAND2005-4647, Sandia National Laboratories, 2005",
                    url: "https://www.osti.gov/servlets/purl/876325",
                    class: DurationSourceClass::GovLabTimeStudy,
                    caveat: "defense-lab problem scale bounds context, not the representative \
                        task; 2005; heat-transfer n=5 with high variance",
                },
            ],
        },
        BaselineStepEstimate {
            step: IncumbentStep::Meshing,
            low_hours: 0.5,
            high_hours: 24.0,
            basis: "floor: automated gridding in a CAD-embedded tool; ceiling: the 10%-over-one-\
                day survey tail for electronics-scale models",
            sources: &[
                SourcedFigure {
                    figure: "41% spend >1 h gridding; 10% >1 day",
                    quote: "41 percent also said that they typically spend over an hour gridding \
                        their designs, while 10 percent spend more than a day.",
                    citation: "6SigmaET State of Thermal survey via engineering.com, 2018",
                    url: "https://www.engineering.com/70-of-engineers-left-out-in-the-cold-by-thermal-simulations/",
                    class: DurationSourceClass::PractitionerSurvey,
                    caveat: "vendor-run survey; meshing automation is the sponsor's key \
                        differentiator",
                },
                SourcedFigure {
                    figure: "1.9M-cell mesh in <3 h wall-clock (detailed package + test board)",
                    quote: "1.9-million-cell mesh in less than three hours of computation on a \
                        12-core Intel Xeon processor.",
                    citation: "Thales Corporate Engineering / Flotherm XT case via ETS Solution \
                        Asia (Siemens reseller)",
                    url: "https://www.etssolution-asia.com/blog/simcenter-flotherm-accurate-simulations-for-electronics-thermal-design",
                    class: DurationSourceClass::VendorCaseStudy,
                    caveat: "vendor-channel case; favorable framing; one detailed package, not a \
                        full assembly",
                },
                SourcedFigure {
                    figure: "meshing 19.6 h + mesh manipulation 11.3 h initial engineer-hours \
                        (heat transfer)",
                    quote: "Meshing 19.62; Mesh Manipulation 11.30 (Table 5, heat-transfer \
                        column, engineering hours)",
                    citation: "Sandia SAND2005-4647 DART System Analysis, 2005",
                    url: "https://www.osti.gov/servlets/purl/876325",
                    class: DurationSourceClass::GovLabTimeStudy,
                    caveat: "defense-lab scale; thermal is DART's own exception where meshing \
                        does NOT dominate",
                },
            ],
        },
        BaselineStepEstimate {
            step: IncumbentStep::SolverSetup,
            low_hours: 0.5,
            high_hours: 16.0,
            basis: "floor: existing material/component libraries; ceiling: board-level power \
                map and interface cards assembled by hand",
            sources: &[
                SourcedFigure {
                    figure: "39% spend >1 h defining material properties",
                    quote: "39 percent spend more than an hour defining properties",
                    citation: "6SigmaET State of Thermal survey via engineering.com, 2018",
                    url: "https://www.engineering.com/70-of-engineers-left-out-in-the-cold-by-thermal-simulations/",
                    class: DurationSourceClass::PractitionerSurvey,
                    caveat: "vendor-run survey; self-reported",
                },
                SourcedFigure {
                    figure: "23% spend >1 h on boundary conditions",
                    quote: "23 percent must spend over an hour setting boundary conditions.",
                    citation: "6SigmaET State of Thermal survey via engineering.com, 2018",
                    url: "https://www.engineering.com/70-of-engineers-left-out-in-the-cold-by-thermal-simulations/",
                    class: DurationSourceClass::PractitionerSurvey,
                    caveat: "vendor-run survey; self-reported",
                },
                SourcedFigure {
                    figure: "parameter assignment 19.4 h; simulation-model assembly 55.4 h \
                        (heat transfer, initial)",
                    quote: "Model Parameter Assignment 19.40; Simulation Model Assembly 55.40 \
                        (Table 5, heat-transfer column, engineering hours)",
                    citation: "Sandia SAND2005-4647 DART System Analysis, 2005",
                    url: "https://www.osti.gov/servlets/purl/876325",
                    class: DurationSourceClass::GovLabTimeStudy,
                    caveat: "defense-lab scale bounds context only",
                },
            ],
        },
        BaselineStepEstimate {
            step: IncumbentStep::Solve,
            low_hours: 0.5,
            high_hours: 48.0,
            basis: "floor: module-scale steady case on current hardware; ceiling: turbulence-\
                model-heavy chassis CHT counted 1:1 into cycle time because a running solve \
                blocks the iteration regardless of engineer attention",
            sources: &[
                SourcedFigure {
                    figure: "66% spend up to a day or more per solve; 14% <30 min",
                    quote: "66 percent of thermal engineers spend up to a day or more solving \
                        their simulations—and yet 14 percent said that they can perform this \
                        task in under 30 minutes.",
                    citation: "6SigmaET State of Thermal survey via engineering.com, 2018",
                    url: "https://www.engineering.com/70-of-engineers-left-out-in-the-cold-by-thermal-simulations/",
                    class: DurationSourceClass::PractitionerSurvey,
                    caveat: "vendor-run survey; 'up to a day or more' is loose phrasing",
                },
                SourcedFigure {
                    figure: "<4.5 h to full convergence, 1.9M cells, 12-core Xeon, 10 GB RAM",
                    quote: "The powerful solver needed less than 4.5 hours to reach full \
                        convergence on the same processor using 10 Gb of memory.",
                    citation: "Thales / Flotherm XT case via ETS Solution Asia",
                    url: "https://www.etssolution-asia.com/blog/simcenter-flotherm-accurate-simulations-for-electronics-thermal-design",
                    class: DurationSourceClass::VendorCaseStudy,
                    caveat: "vendor-channel case; single detailed package plus board",
                },
                SourcedFigure {
                    figure: "Icepak transient runs 730-4800 s (~12-80 min), 0.3-3M cells",
                    quote: "Table 3 reports Ansys Icepak run times of 730-4800 s across \
                        natural/forced-convection transient cases (Intel i7-10850H, 6 parallel \
                        processes).",
                    citation: "Padmanabhan (ZF), A Scalable Approach for Transient Thermal \
                        Modeling of Automotive Power Electronics, SAE 2025-01-5073 (author \
                        manuscript), 2025",
                    url: "https://arxiv.org/pdf/2606.11226",
                    class: DurationSourceClass::PeerReviewedRun,
                    caveat: "module scale on laptop-class hardware",
                },
                SourcedFigure {
                    figure: "data-hall CHT steady solve 0.44-5.95 h (32 to 1 cores), 10M cells",
                    quote: "CFD model solving time: 5.95 h (1 core) ... 0.44 h (32 cores) \
                        (Table 2, 6SigmaDCX, one production data hall).",
                    citation: "Wang, Zhou, Dong, Wen, Tan, Chen, Wang, Zeng, Kalibre, ACM \
                        BuildSys 2020",
                    url: "https://arxiv.org/pdf/2001.10681",
                    class: DurationSourceClass::PeerReviewedRun,
                    caveat: "data-center scale, larger than the representative task",
                },
                SourcedFigure {
                    figure: "RNG k-epsilon chassis CHT: two days of continuous runs (900k \
                        cells, 2008 hardware)",
                    quote: "Using the RNG k-e model, which is a two-equation model, doubles the \
                        solution time. This corresponds to two days of continuous runs.",
                    citation: "Ozturk, Tari, Forced Air Cooling of CPUs with Heat Sinks: a \
                        Numerical Study, IEEE Trans. Components and Packaging Technologies \
                        31(3), 2008",
                    url: "https://users.metu.edu.tr/itari/OzturkTariIEEETCAPT2008.pdf",
                    class: DurationSourceClass::PeerReviewedRun,
                    caveat: "2008 workstation hardware; modern hosts are faster, retained as \
                        the turbulence-model ceiling",
                },
            ],
        },
        BaselineStepEstimate {
            step: IncumbentStep::PostProcessing,
            low_hours: 1.0,
            high_hours: 16.0,
            basis: "floor: standard temperature-map review; ceiling: scaled down from the DART \
                heat-transfer post-processing share (37% of engineer-time, the largest single \
                step) to representative-task size",
            sources: &[
                SourcedFigure {
                    figure: "post-processing 123.4 h initial / 245.3 h expected = 37% of \
                        engineer-time for heat-transfer problems",
                    quote: "Results Post-Processing 123.40 (Table 5); 245.3 h of a 659.2 h \
                        total (Table 7); normalized 0.37 (Table 8) — the largest single \
                        fraction for heat-transfer problems.",
                    citation: "Sandia SAND2005-4647 DART System Analysis, 2005",
                    url: "https://www.osti.gov/servlets/purl/876325",
                    class: DurationSourceClass::GovLabTimeStudy,
                    caveat: "defense-lab scale; strongest available signal that post-processing, \
                        not meshing, dominates thermal analysis",
                },
                SourcedFigure {
                    figure: "thermal-to-structural results handoff: half a day of manual \
                        scripting per transfer",
                    quote: "the manual baseline (MATLAB scripts, multiple engineers, about half \
                        a day per transfer)",
                    citation: "Siemens Mynaric case study",
                    url: "https://resources.sw.siemens.com/en-US/case-study-mynaric/",
                    class: DurationSourceClass::VendorCaseStudy,
                    caveat: "vendor case; space thermal, not forced-convection PCB cooling",
                },
            ],
        },
        BaselineStepEstimate {
            step: IncumbentStep::ReportAssembly,
            low_hours: 1.0,
            high_hours: 8.0,
            basis: "weakest row: no independent electronics-specific figure exists (explicit \
                evidence gap); bounds set from the DART post-processing share net of \
                visualization review, flagged for the quarterly review",
            sources: &[SourcedFigure {
                figure: "report generation folded into DART post-processing step I \
                        (123.4 h initial, heat transfer)",
                quote: "Results Post-Processing ... explicitly includes 'generate reports' \
                        in step I of the DART process model.",
                citation: "Sandia SAND2005-4647 DART System Analysis, 2005",
                url: "https://www.osti.gov/servlets/purl/876325",
                class: DurationSourceClass::GovLabTimeStudy,
                caveat: "report time is not separated from visualization in the primary; \
                        this row is the envelope's weakest and says so",
            }],
        },
    ],
    hours_per_working_day: 8.0,
    target_reduction: 3.0,
    kill_within_quarters: 2,
};

/// One named go-to-market audit check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuditCheck {
    /// The check name.
    pub name: &'static str,
    /// Did it pass?
    pub passed: bool,
}

/// The go-to-market audit result.
#[derive(Debug, Clone, PartialEq)]
pub struct WedgeAudit {
    /// Named checks for supersession, evidence completeness, score/status
    /// consistency, explicit comparison/ranking/sensitivity shape,
    /// rank/proposal shape, and the cycle-time criterion.
    pub checks: Vec<AuditCheck>,
    /// Any gaps (human-readable).
    pub gaps: Vec<String>,
}

impl WedgeAudit {
    /// Is the go-to-market story complete and self-consistent?
    #[must_use]
    pub fn ok(&self) -> bool {
        self.gaps.is_empty()
    }

    /// Did a named check pass?
    #[must_use]
    pub fn passed(&self, name: &str) -> bool {
        self.checks.iter().any(|c| c.name == name && c.passed)
    }
}

/// The historical threshold retained for replay and the forbidden lower bound
/// for an absent capability.
pub const STRONG_THRESHOLD: u8 = 8;

fn audit_comparison(gaps: &mut Vec<String>) -> [AuditCheck; 4] {
    let inputs_complete = COMPARISON_CANDIDATES.len() == 3
        && COMPARISON_CANDIDATES
            .iter()
            .all(|candidate| validate_candidate(candidate).is_ok());
    if !inputs_complete {
        gaps.push("the explicit comparison lacks a complete measured factor table".to_string());
    }

    let weights_normalized = canonical_weight_values(&DEFAULT_FACTOR_WEIGHTS).is_ok();
    if !weights_normalized {
        gaps.push("the default comparison weights do not normalize to 100".to_string());
    }

    let recommendation = default_recommendation();
    let ranking_complete = recommendation.as_ref().is_ok_and(|record| {
        record.ranked.len() == COMPARISON_CANDIDATES.len()
            && record.recommended == "thermal-design-assurance"
            && !record.runner_up.is_empty()
            && !record.minority_report.is_empty()
    });
    if !ranking_complete {
        gaps.push(
            "the default explicit comparison does not produce its recorded ranking".to_string(),
        );
    }

    let sensitivity_complete = recommendation.as_ref().is_ok_and(|record| {
        let expected = (COMPARISON_CANDIDATES.len() - 1) * ScoringFactor::ALL.len();
        record.rating_sensitivities.len() == expected
            && record.weight_sensitivities.len() == expected
            && record
                .rating_sensitivities
                .iter()
                .any(|row| row.challenger == record.runner_up && row.required_rating.is_some())
            && record
                .weight_sensitivities
                .iter()
                .any(|row| row.challenger == record.runner_up && row.required_weight.is_some())
    });
    if !sensitivity_complete {
        gaps.push("the comparison lacks complete one-factor flip sensitivities".to_string());
    }

    [
        AuditCheck {
            name: "comparison-inputs-complete",
            passed: inputs_complete,
        },
        AuditCheck {
            name: "default-weights-normalized",
            passed: weights_normalized,
        },
        AuditCheck {
            name: "comparison-ranking-complete",
            passed: ranking_complete,
        },
        AuditCheck {
            name: "comparison-sensitivity-complete",
            passed: sensitivity_complete,
        },
    ]
}

/// Audit the wedge input ledger.
///
/// Historical scores must be superseded, every candidate needs complete
/// measurements on all four axes, absent inputs may not carry strong authority
/// scores, the comparison must be complete and normalized, ranks/proposal
/// mappings must remain complete, and the kill-criterion shape must remain
/// `>= 3×`.
#[must_use]
pub fn audit() -> WedgeAudit {
    let mut gaps = Vec::new();

    let historic_scores_superseded = VERTICALS
        .iter()
        .all(|vertical| !vertical.score_use.permits_decision());
    if !historic_scores_superseded {
        gaps.push("a historical plan score still permits decision use".to_string());
    }

    let measured_inputs_complete = MEASURED_INPUTS.len() == VERTICALS.len()
        && VERTICALS.iter().all(|vertical| {
            measured_inputs_for(vertical.name).is_some_and(MeasuredWedgeInputs::is_complete)
        });
    if !measured_inputs_complete {
        gaps.push("a candidate lacks a complete four-axis measured-input record".to_string());
    }

    let no_absent_strong_scores = MEASURED_INPUTS.iter().all(|inputs| {
        inputs.measurements().all(|measurement| {
            measurement.readiness != Readiness::Absent || measurement.score < STRONG_THRESHOLD
        })
    }) && COMPARISON_CANDIDATES.iter().all(|candidate| {
        candidate.factors.iter().all(|input| {
            input.measurement.readiness != Readiness::Absent
                || input.measurement.score < STRONG_THRESHOLD
        })
    });
    if !no_absent_strong_scores {
        gaps.push(format!(
            "an absent input carries a score at or above {STRONG_THRESHOLD}"
        ));
    }

    let comparison_checks = audit_comparison(&mut gaps);

    let mut ranks: Vec<u8> = VERTICALS.iter().map(|v| v.rank).collect();
    ranks.sort_unstable();
    let ranks_complete = ranks == vec![1, 2, 3];
    if !ranks_complete {
        gaps.push("verticals are not ranked exactly 1, 2, 3".to_string());
    }

    let all_exercise_proposals = VERTICALS.iter().all(|v| !v.exercises.is_empty());
    if !all_exercise_proposals {
        gaps.push("a vertical names no exercised proposal".to_string());
    }

    let kill_criterion_measurable = (CHT_BASELINE.target_reduction - 3.0).abs() < f64::EPSILON;
    if !kill_criterion_measurable {
        gaps.push("cycle-time kill criterion is not the required 3x".to_string());
    }

    let baseline_measured =
        CHT_BASELINE.provenance != BaselineProvenance::Placeholder && CHT_BASELINE.is_complete();
    if !baseline_measured {
        gaps.push("the cycle-time baseline is not a complete measured record".to_string());
    }

    let placeholder_refused = matches!(
        RETIRED_PLACEHOLDER_BASELINE.evaluate_kill_criterion(1.0),
        Err(KillCriterionError::PlaceholderBaseline { .. })
    );
    if !placeholder_refused {
        gaps.push(
            "the retired placeholder baseline is not refused by the kill criterion".to_string(),
        );
    }

    let mut checks = vec![
        AuditCheck {
            name: "historic-scores-superseded",
            passed: historic_scores_superseded,
        },
        AuditCheck {
            name: "measured-inputs-complete",
            passed: measured_inputs_complete,
        },
        AuditCheck {
            name: "no-absent-strong-scores",
            passed: no_absent_strong_scores,
        },
    ];
    checks.extend(comparison_checks);
    checks.extend([
        AuditCheck {
            name: "ranks-complete",
            passed: ranks_complete,
        },
        AuditCheck {
            name: "all-exercise-proposals",
            passed: all_exercise_proposals,
        },
        AuditCheck {
            name: "kill-criterion-measurable",
            passed: kill_criterion_measurable,
        },
        AuditCheck {
            name: "cycle-time-baseline-measured",
            passed: baseline_measured,
        },
        AuditCheck {
            name: "placeholder-baseline-refused",
            passed: placeholder_refused,
        },
    ]);

    WedgeAudit { checks, gaps }
}

fn push_json_string(out: &mut String, value: &str) {
    use core::fmt::Write as _;
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            character if character.is_control() => {
                write!(out, "\\u{:04x}", u32::from(character)).expect("write to String");
            }
            character => out.push(character),
        }
    }
    out.push('"');
}

fn write_measurement_fields(out: &mut String, axis: InputAxis, measurement: Measurement) {
    use core::fmt::Write as _;
    out.push_str("\"axis\":");
    push_json_string(out, axis.label());
    out.push_str(",\"readiness\":");
    push_json_string(out, measurement.readiness.label());
    write!(out, ",\"score\":{}", measurement.score).expect("write to String");
    out.push_str(",\"method\":");
    push_json_string(out, measurement.method.label());
    out.push_str(",\"finding\":");
    push_json_string(out, measurement.finding);
    out.push_str(",\"evidence\":[");
    for (index, pointer) in measurement.evidence.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"kind\":");
        push_json_string(out, pointer.kind.label());
        out.push_str(",\"reference\":");
        push_json_string(out, pointer.reference);
        out.push_str(",\"locator\":");
        push_json_string(out, pointer.locator);
        out.push('}');
    }
    out.push(']');
}

/// Emit the historical ranking and complete measured-input ledger as
/// deterministic machine-readable JSON.
#[must_use]
pub fn to_json() -> String {
    use core::fmt::Write as _;
    let mut out = String::from("{\"verticals\":[");
    for (index, vertical) in VERTICALS.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"name\":");
        push_json_string(&mut out, vertical.name);
        write!(
            out,
            ",\"rank\":{},\"historic_weakest_score\":{}",
            vertical.rank,
            vertical.weakest_criterion_score()
        )
        .expect("write to String");
        out.push_str(",\"score_use\":");
        push_json_string(&mut out, vertical.score_use.label());
        out.push_str(",\"exercises\":[");
        for (proposal_index, proposal) in vertical.exercises.iter().enumerate() {
            if proposal_index > 0 {
                out.push(',');
            }
            push_json_string(&mut out, proposal);
        }
        out.push_str("]}");
    }

    out.push_str("],\"measured_inputs\":[");
    for (index, inputs) in MEASURED_INPUTS.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str("{\"vertical\":");
        push_json_string(&mut out, inputs.vertical);
        out.push_str(",\"measured_on\":");
        push_json_string(&mut out, inputs.measured_on);

        out.push_str(",\"kernel_readiness\":[");
        for (entry_index, entry) in inputs.kernels.iter().enumerate() {
            if entry_index > 0 {
                out.push(',');
            }
            out.push_str("{\"capability\":");
            push_json_string(&mut out, entry.capability);
            out.push(',');
            write_measurement_fields(&mut out, InputAxis::KernelReadiness, entry.measurement);
            out.push('}');
        }

        out.push_str("],\"validation_data\":[");
        for (entry_index, entry) in inputs.validation_data.iter().enumerate() {
            if entry_index > 0 {
                out.push(',');
            }
            out.push_str("{\"dataset\":");
            push_json_string(&mut out, entry.dataset);
            out.push_str(",\"raw_data\":");
            push_json_string(&mut out, entry.raw_data);
            out.push_str(",\"license_terms\":");
            push_json_string(&mut out, entry.license_terms);
            out.push(',');
            write_measurement_fields(&mut out, InputAxis::ValidationDataAccess, entry.measurement);
            out.push('}');
        }

        out.push_str("],\"cad_burden\":[");
        for (entry_index, entry) in inputs.cad_burden.iter().enumerate() {
            if entry_index > 0 {
                out.push(',');
            }
            out.push_str("{\"required_geometry\":");
            push_json_string(&mut out, entry.required_geometry);
            out.push_str(",\"admitted_geometry\":");
            push_json_string(&mut out, entry.admitted_geometry);
            out.push_str(",\"missing_semantics\":");
            push_json_string(&mut out, entry.missing_semantics);
            out.push(',');
            write_measurement_fields(&mut out, InputAxis::CadBurden, entry.measurement);
            out.push('}');
        }

        out.push_str("],\"compute_cost\":[");
        for (entry_index, entry) in inputs.compute_cost.iter().enumerate() {
            if entry_index > 0 {
                out.push(',');
            }
            out.push_str("{\"rung\":");
            push_json_string(&mut out, entry.rung);
            out.push_str(",\"variables\":");
            push_json_string(&mut out, entry.variables);
            out.push_str(",\"work_envelope\":");
            push_json_string(&mut out, entry.work_envelope);
            out.push(',');
            write_measurement_fields(&mut out, InputAxis::ComputeCost, entry.measurement);
            out.push('}');
        }
        out.push_str("]}");
    }

    out.push_str("],\"cycle_time_baseline\":");
    write_baseline_json(&mut out);
    out.push('}');
    out
}

fn write_baseline_json(out: &mut String) {
    use core::fmt::Write as _;
    out.push_str("{\"vertical\":");
    push_json_string(out, CHT_BASELINE.vertical);
    out.push_str(",\"assembled_on\":");
    push_json_string(out, CHT_BASELINE.assembled_on);
    out.push_str(",\"provenance\":");
    push_json_string(out, CHT_BASELINE.provenance.label());
    out.push_str(",\"dossier\":");
    push_json_string(out, CHT_BASELINE.dossier.reference);
    write!(
        out,
        ",\"baseline_days_low\":{},\"baseline_days_high\":{},\"hours_per_working_day\":{},\"target_reduction\":{},\"kill_within_quarters\":{},\"steps\":[",
        CHT_BASELINE.baseline_days_low(),
        CHT_BASELINE.baseline_days_high(),
        CHT_BASELINE.hours_per_working_day,
        CHT_BASELINE.target_reduction,
        CHT_BASELINE.kill_within_quarters
    )
    .expect("write to String");
    for (step_index, step) in CHT_BASELINE.steps.iter().enumerate() {
        if step_index > 0 {
            out.push(',');
        }
        out.push_str("{\"step\":");
        push_json_string(out, step.step.label());
        write!(
            out,
            ",\"low_hours\":{},\"high_hours\":{},\"sources\":{}}}",
            step.low_hours,
            step.high_hours,
            step.sources.len()
        )
        .expect("write to String");
    }
    out.push_str("]}");
}
