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
//! a method and evidence pointer. The [`CycleTimeBaseline`] remains a
//! separately identified placeholder until a customer measurement replaces it.

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
            Readiness::Absent,
            2,
            MeasurementMethod::WorkspaceInventory,
            &[evidence(
                EvidenceKind::Bead,
                "frankensim-extreal-program-f85xj.5.1",
                "steady heat-conduction implementation is successor work, not shipped authority",
            )],
            "No tracked production crate owned steady conduction at this snapshot; an in-progress Bead is not decision-usable kernel evidence.",
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

/// The baseline that makes the cycle-time kill criterion measurable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CycleTimeBaseline {
    /// Which vertical.
    pub vertical: &'static str,
    /// Today's baseline design-cycle time (days) for the acknowledged loop.
    pub baseline_days: f64,
    /// The cycle-time reduction factor the kill criterion demands (`3.0`).
    pub target_reduction: f64,
    /// The window (quarters after GA) to hit it or re-select the wedge.
    pub kill_within_quarters: u8,
}

impl CycleTimeBaseline {
    /// Does a measured cycle time meet the `>=target_reduction×` kill
    /// criterion? (`baseline / measured >= target_reduction`.)
    #[must_use]
    pub fn meets_kill_criterion(&self, measured_days: f64) -> bool {
        measured_days > 0.0 && self.baseline_days / measured_days >= self.target_reduction
    }
}

/// The conjugate-heat-transfer cycle-time baseline (a week-long loop today).
pub const CHT_BASELINE: CycleTimeBaseline = CycleTimeBaseline {
    vertical: "conjugate-heat-transfer",
    baseline_days: 5.0,
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
    /// consistency, rank/proposal shape, and the cycle-time criterion.
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

/// Audit the wedge input ledger.
///
/// Historical scores must be superseded, every candidate needs complete
/// measurements on all four axes, absent inputs may not carry strong scores,
/// ranks/proposal mappings must remain complete, and the kill-criterion shape
/// must remain `>= 3×`.
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
    });
    if !no_absent_strong_scores {
        gaps.push(format!(
            "an absent input carries a score at or above {STRONG_THRESHOLD}"
        ));
    }

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

    WedgeAudit {
        checks: vec![
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
        ],
        gaps,
    }
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

    write!(
        out,
        "],\"baseline_days\":{},\"target_reduction\":{}}}",
        CHT_BASELINE.baseline_days, CHT_BASELINE.target_reduction
    )
    .expect("write to String");
    out
}
