//! Reference encoding of a REAL in-repo V&V fixture, for the corpus schema
//! (bead `frankensim-extreal-program-f85xj.4.1`).
//!
//! # Status: not wired into `lib.rs`
//!
//! This module is deliberately NOT declared from `lib.rs`, so it compiles
//! only if someone wires it. It exists as a handoff artifact: it records,
//! field by field, what the one genuinely measured fixture already in this
//! repository actually says, so that whoever extends
//! [`crate::corpus`](crate::corpus) can seed it without re-deriving the
//! provenance. It deliberately depends on NO schema type, so it cannot rot
//! against a schema still under construction.
//!
//! # Why this fixture
//!
//! `data/reference/martin-moyce-1952.jsonl` is 946 bytes of digitized
//! surge-front coordinates from Martin & Moyce (1952), and
//! `crates/fs-lbm/tests/d3q19_freesurface3.rs::lbm3_105_martin_moyce_front`
//! already consumes it. It is therefore an EXISTING fixture (the bead's
//! acceptance wording) and a published physical experiment (portfolio level
//! C), as opposed to a fixture authored to demonstrate the schema.
//!
//! # Why it cannot be seeded today
//!
//! Every completeness axis except the payload itself is honestly ABSENT,
//! and a corpus schema can only hold this row if absence is representable:
//!
//! | schema field | honest state for this fixture |
//! | --- | --- |
//! | raw sensor payload | ABSENT — digitized figure coordinates only; no cine frames, no timing records, no length-scale frames were ever retained here |
//! | calibration certificate | ABSENT — none retained, none reproduced in the secondary sources |
//! | instrument identity | ABSENT — unrecorded in the retained artifact |
//! | sensor placement + placement uncertainty | ABSENT — measurement station and camera geometry unrecorded |
//! | as-built geometry | ABSENT — nominal only (see [`GEOMETRY_NOMINAL`]) |
//! | environmental conditions | ABSENT — ambient/fluid temperature and surface state unrecorded |
//! | measurement uncertainty | ABSENT — see [`UNCERTAINTY_NOTE`]; a qualitative phrase is not a bound |
//! | preprocessing lineage | ABSENT — figure scan, digitizer tool and operator all unknown |
//! | acquisition window | ABSENT — 1952 is a publication year, not an acquisition window |
//! | acceptance envelope | UNPINNABLE as a scalar rule — see [`ACCEPTANCE_BASIS`] |
//!
//! A schema that requires each of those to be present forces an author to
//! invent them, which inverts the purpose of the corpus: the only rows it
//! can hold are the ones nobody measured. Representing each absence as an
//! explicit typed variant (as
//! `crate::corpus::MeasurementUncertainty::Unstated` already does for
//! uncertainty) is what lets this row exist AND caps it at `Estimated`.

/// Proposed dataset id.
pub const DATASET_ID: &str = "martin-moyce-1952-square-column";

/// Repo-relative path of the retained payload.
pub const PAYLOAD_PATH: &str = "data/reference/martin-moyce-1952.jsonl";

/// Byte length of the retained payload at the time of writing.
pub const PAYLOAD_BYTES: u64 = 946;

/// Media type of the retained payload.
pub const PAYLOAD_MEDIA_TYPE: &str = "application/x-ndjson";

/// The existing in-repo consumer of this fixture.
pub const EXISTING_CONSUMER: &str =
    "crates/fs-lbm/tests/d3q19_freesurface3.rs::lbm3_105_martin_moyce_front";

/// Portfolio evidence level: C, published experiment.
pub const EVIDENCE_LEVEL: &str = "C";

/// Declared partition role.
pub const PARTITION: &str = "validation";

/// Why the partition role is validation and not calibration.
pub const PARTITION_RATIONALE: &str = "the curve is compared against free-surface solver output only; no FrankenSim model \
     parameter, closure coefficient, or lattice setting is fitted to it";

/// Human citation.
pub const CITATION: &str = "J. C. Martin & W. J. Moyce (1952), 'Part IV. An experimental study of \
                            the collapse of liquid columns on a rigid horizontal plane', Phil. \
                            Trans. R. Soc. Lond. A 244, 312-324";

/// Exact locator.
pub const LOCATOR: &str = "Phil. Trans. R. Soc. Lond. A 244:312-324 (1952)";

/// Who measured it, as far as the retained artifact records.
pub const MEASURED_BY: &str = "J. C. Martin and W. J. Moyce (original experiment); the digitizer \
                               of these coordinates is unrecorded";

/// Raw-retention state. The payload is post-processed only.
pub const RETENTION_NOTE: &str = "only digitized figure coordinates survive in this repository: no cine frames, no raw timing \
     records, and no length-scale calibration frames were ever acquired or retained here";

/// Calibration state.
pub const CALIBRATION_NOTE: &str = "no calibration certificate for the 1952 imaging, timing, or length-scale chain is retained \
     here, and none is reproduced in the secondary sources this curve is taken from";

/// Nominal geometry, as stated in the retained artifact.
pub const GEOMETRY_NOMINAL: &str = "square-based water column collapsing on a rigid horizontal plane; base a = 2.25 in, initial \
     aspect ratio n^2 = 2";

/// Why the geometry is nominal-only.
pub const GEOMETRY_NOTE: &str = "as-built tank dimensions, plate flatness, and gate-release \
                                 geometry are not recorded in the retained artifact";

/// THE crux field. The source states a qualitative phrase, not a bound.
pub const UNCERTAINTY_NOTE: &str = "the retained artifact states only that 'digitization uncertainty is a few percent' — a \
     qualitative phrase with no half-width and no confidence level. No measurement covariance, \
     repeatability record, or original figure-reading uncertainty survives. Converting that \
     phrase into a number would invent a bound the source never stated, so it is Unstated and \
     caps use at Estimated.";

/// Why the preprocessing lineage cannot be recorded as complete.
pub const LINEAGE_NOTE: &str = "the chain from the published figure to these coordinates is not replayable: the source \
     figure scan, the digitizing tool and its version, and the operator are all unrecorded. \
     Exactly one transform is known to have occurred (published figure -> coordinate pairs) and \
     none of its parameters are retained.";

/// The decision this dataset is evidence for.
pub const CONTEXT_DECISION: &str = "coarse-lattice free-surface dam-break front-position comparison for the FrankenSim LBM/SPH \
     free-surface batteries";

/// Context axis name. `T = t*sqrt(2g/a)` in the retained artifact.
pub const CONTEXT_AXIS: &str = "t_star";

/// Inclusive lower end of the span the retained points actually cover.
pub const CONTEXT_LO: f64 = 0.41;

/// Inclusive upper end of the span the retained points actually cover.
pub const CONTEXT_HI: f64 = 2.95;

/// Uses this dataset does NOT support.
pub const CONTEXT_EXCLUSIONS: &[&str] = &[
    "pressure, impact load, or force prediction",
    "splash, air entrainment, or fragmentation metrics",
    "quantitative central-band acceptance at any lattice resolution",
    "any initial geometry other than the square-based n^2 = 2 column",
    "any working fluid other than the water column of the cited experiment",
];

/// The metric the acceptance record would name. `Z = x/a`, dimensionless.
pub const ACCEPTANCE_METRIC: &str = "surge-front-position-z";

/// Inclusive lower end of the regime the in-repo consumer actually gates in.
pub const ACCEPTANCE_REGIME_LO: f64 = 0.5;

/// Inclusive upper end of the regime the in-repo consumer actually gates in.
pub const ACCEPTANCE_REGIME_HI: f64 = 2.0;

/// Why no scalar acceptance envelope is defensible for this metric.
pub const ACCEPTANCE_BASIS: &str = "the in-repo consumer applies a monotone-advance check plus the broad upper envelope \
     z <= 2.2*t_star + 1 for 0.5 < t_star < 2, and compares this curve REPORT-ONLY. That gate is \
     a function of t_star, which a scalar tolerance/interval algebra cannot express, and no \
     quantitative central band is defensible while the digitization uncertainty is only a \
     qualitative phrase. The envelope should therefore be recorded as unpinned-with-a-basis, \
     never invented.";

/// License state of the retained bytes.
pub const LICENSE_TERMS: &str = "the underlying figure was published in Phil. Trans. R. Soc. Lond. A 244 (1952) and is Royal \
     Society copyright; the retained bytes are numeric coordinate pairs only — no figure, plate, \
     table image, or text is reproduced";

/// Redistribution state. Deliberately unresolved rather than asserted.
pub const REDISTRIBUTION_NOTE: &str = "redistribution terms for these coordinates have NOT been established. Numeric coordinate \
     values reproduced across the LBM/SPH validation literature are commonly treated as facts, \
     but this repository has obtained no determination; packaging must resolve terms before any \
     release ships this dataset.";
