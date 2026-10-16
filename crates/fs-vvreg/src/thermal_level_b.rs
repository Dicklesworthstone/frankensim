//! Level-B thermal cross-code frozen references (EXTREAL E04, bead
//! `frankensim-extreal-program-f85xj.4.3`).
//!
//! Each case is a steady-conduction problem solved by an EXTERNAL,
//! independently implemented FEM code (scikit-fem on numpy/scipy with
//! SuperLU solves; pinned by `tools/vvref/uv.lock`) on the exact
//! Kuhn/Freudenthal mesh `fs_conduction::fixtures` builds, with the same
//! declared discretization (P1, element-mean `k(T)`, consistent
//! source/Robin mass rules). The frozen probe temperatures in the retained
//! manifest are therefore SAME-DISCRETIZATION cross-code parity
//! references: they check independent assembly/boundary-condition/solver
//! implementations against each other, not discretization error, and they
//! are NOT physical validation. Two codes agreeing is not truth — the
//! corpus registers every row as an Estimated-colour `CrossCode` dataset.
//!
//! Fail-closed binding runs in three layers:
//! 1. the manifest parser refuses malformed rows, missing mandatory
//!    metadata, and a failed external self-check;
//! 2. [`verify_spec_echo`] compares every load-bearing number the external
//!    solver echoed back against this catalog's typed constants
//!    BIT-EXACTLY, so the frozen values cannot drift from the case they
//!    claim to solve;
//! 3. [`verify_probe_grid`] recomputes every probe's mesh position with
//!    the fixture arithmetic (`index * (extent / count)`) and requires
//!    bit-identity, witnessing the cross-language mesh-parity claim
//!    without running the external stack.
//!
//! Reproducibility contract: `data/vv-corpus/thermal-level-b/README.md`.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::LazyLock;

/// Retained tab-separated manifest backing every Level-B thermal row.
pub(crate) const THERMAL_LEVEL_B_MANIFEST: &[u8] =
    include_bytes!("../../../data/vv-corpus/thermal-level-b/thermal-level-b-references-v1.tsv");

/// Repository-relative locator of the retained manifest.
pub const THERMAL_LEVEL_B_MANIFEST_LOCATOR: &str =
    "data/vv-corpus/thermal-level-b/thermal-level-b-references-v1.tsv";

/// Committed external-deck bytes, keyed by case id, for hash binding.
pub(crate) const THERMAL_LEVEL_B_DECKS: [(&str, &[u8]); 4] = [
    (
        "thermal-b-block-source-robin-v1",
        include_bytes!(
            "../../../data/vv-corpus/thermal-level-b/thermal-b-block-source-robin-v1.case.json"
        ),
    ),
    (
        "thermal-b-orthotropic-rotated-v1",
        include_bytes!(
            "../../../data/vv-corpus/thermal-level-b/thermal-b-orthotropic-rotated-v1.case.json"
        ),
    ),
    (
        "thermal-b-kt-nonlinear-slab-v1",
        include_bytes!(
            "../../../data/vv-corpus/thermal-level-b/thermal-b-kt-nonlinear-slab-v1.case.json"
        ),
    ),
    (
        "thermal-b-fin-film-v1",
        include_bytes!("../../../data/vv-corpus/thermal-level-b/thermal-b-fin-film-v1.case.json"),
    ),
];

/// The committed external-deck bytes for a case, if the case exists.
#[must_use]
pub fn thermal_level_b_deck_bytes(case_id: &str) -> Option<&'static [u8]> {
    THERMAL_LEVEL_B_DECKS
        .iter()
        .find(|(id, _)| *id == case_id)
        .map(|(_, bytes)| *bytes)
}

/// Conductivity model of a Level-B case, matching the committed deck.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalLevelBMaterial {
    /// Constant isotropic conductivity, W/(m K).
    Isotropic {
        /// Conductivity value.
        k: f64,
    },
    /// Constant symmetric positive-definite conductivity tensor, W/(m K).
    Tensor {
        /// Row-major 3x3 tensor.
        k: [[f64; 3]; 3],
    },
    /// Piecewise-linear `k(T)` evaluated once per element at the
    /// element-mean temperature (the declared discretization).
    LinearKt {
        /// `(T, k)` knots with strictly increasing temperatures.
        knots: &'static [(f64, f64)],
    },
}

/// Volumetric source of a Level-B case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalLevelBSource {
    /// No volumetric source.
    None,
    /// `q(x, y) = q0 * (x / ex) * (1.0 - y / ey)` W/m^3, evaluated at
    /// mesh nodes in exactly this IEEE-754 operation order by both codes.
    PolyXy {
        /// Source amplitude, W/m^3.
        q0: f64,
    },
}

/// Boundary condition attached to one axis-aligned box face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalLevelBBcKind {
    /// Prescribed temperature, K.
    Dirichlet {
        /// Face temperature.
        t_k: f64,
    },
    /// Convective (film) boundary.
    Robin {
        /// Heat-transfer coefficient, W/(m^2 K).
        h: f64,
        /// Ambient temperature, K.
        t_inf_k: f64,
    },
}

/// One named face boundary condition. Unlisted faces are adiabatic
/// (natural) in both codes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalLevelBBc {
    /// Stable name matching the committed deck.
    pub name: &'static str,
    /// Face axis: 0 = x, 1 = y, 2 = z.
    pub axis: usize,
    /// Whether the face sits at `extent[axis]` rather than `0.0`.
    pub at_max: bool,
    /// The condition on that face.
    pub kind: ThermalLevelBBcKind,
}

/// Picard controls for the nonlinear `k(T)` case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalLevelBPicard {
    /// Relative step tolerance.
    pub tol: f64,
    /// Iteration budget; exhausting it is a refusal in both codes.
    pub max_iterations: u64,
}

/// One Level-B cross-code case definition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalLevelBCase {
    /// Stable corpus dataset id.
    pub id: &'static str,
    /// Human-readable title.
    pub title: &'static str,
    /// Kuhn grid cell counts per axis.
    pub mesh_counts: [usize; 3],
    /// Box extent per axis, m.
    pub mesh_extent: [f64; 3],
    /// Conductivity model.
    pub material: ThermalLevelBMaterial,
    /// Volumetric source.
    pub source: ThermalLevelBSource,
    /// Named face boundary conditions, in deck order.
    pub bcs: &'static [ThermalLevelBBc],
    /// Picard controls; present exactly for `LinearKt` materials.
    pub picard: Option<ThermalLevelBPicard>,
    /// Probe locations as vertex grid indices `(i, j, k)`.
    pub probes: &'static [[usize; 3]],
    /// Declared absolute agreement envelope on probe temperatures, K.
    ///
    /// A violation opens an investigation bead; it is never silently
    /// widened (golden-bump discipline).
    pub acceptance_atol_k: f64,
    /// Explicit boundary on what agreement means.
    pub no_claim_reason: &'static str,
}

const CROSS_CODE_ONLY: &str = "same-discretization cross-code parity reference; agreement checks \
     independent implementations, not discretization error, and is not physical validation";

static THERMAL_LEVEL_B_CASES: [ThermalLevelBCase; 4] = [
    ThermalLevelBCase {
        id: "thermal-b-block-source-robin-v1",
        title: "Convectively cooled block with a bilinear volumetric source",
        mesh_counts: [8, 5, 3],
        mesh_extent: [0.08, 0.05, 0.03],
        material: ThermalLevelBMaterial::Isotropic { k: 15.0 },
        source: ThermalLevelBSource::PolyXy { q0: 200_000.0 },
        bcs: &[
            ThermalLevelBBc {
                name: "x-lo",
                axis: 0,
                at_max: false,
                kind: ThermalLevelBBcKind::Robin {
                    h: 25.0,
                    t_inf_k: 300.0,
                },
            },
            ThermalLevelBBc {
                name: "x-hi",
                axis: 0,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 25.0,
                    t_inf_k: 300.0,
                },
            },
            ThermalLevelBBc {
                name: "y-lo",
                axis: 1,
                at_max: false,
                kind: ThermalLevelBBcKind::Robin {
                    h: 25.0,
                    t_inf_k: 300.0,
                },
            },
            ThermalLevelBBc {
                name: "y-hi",
                axis: 1,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 25.0,
                    t_inf_k: 300.0,
                },
            },
            ThermalLevelBBc {
                name: "z-lo",
                axis: 2,
                at_max: false,
                kind: ThermalLevelBBcKind::Robin {
                    h: 25.0,
                    t_inf_k: 300.0,
                },
            },
            ThermalLevelBBc {
                name: "z-hi",
                axis: 2,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 25.0,
                    t_inf_k: 300.0,
                },
            },
        ],
        picard: None,
        probes: &[[0, 0, 0], [4, 2, 2], [8, 5, 3], [2, 4, 1], [6, 1, 3]],
        acceptance_atol_k: 1.0e-6,
        no_claim_reason: CROSS_CODE_ONLY,
    },
    ThermalLevelBCase {
        id: "thermal-b-orthotropic-rotated-v1",
        title: "Anisotropic block with a declared rotated-frame conductivity tensor",
        mesh_counts: [6, 4, 4],
        mesh_extent: [0.06, 0.04, 0.04],
        material: ThermalLevelBMaterial::Tensor {
            k: [
                [23.75, 10.825_317_547_305_483, 0.0],
                [10.825_317_547_305_483, 11.25, 0.0],
                [0.0, 0.0, 1.0],
            ],
        },
        source: ThermalLevelBSource::None,
        bcs: &[
            ThermalLevelBBc {
                name: "hot",
                axis: 0,
                at_max: false,
                kind: ThermalLevelBBcKind::Dirichlet { t_k: 320.0 },
            },
            ThermalLevelBBc {
                name: "cold",
                axis: 0,
                at_max: true,
                kind: ThermalLevelBBcKind::Dirichlet { t_k: 300.0 },
            },
            ThermalLevelBBc {
                name: "cooled-top",
                axis: 1,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 15.0,
                    t_inf_k: 290.0,
                },
            },
        ],
        picard: None,
        probes: &[[3, 2, 2], [1, 1, 3], [5, 3, 1], [3, 0, 0], [2, 4, 2]],
        acceptance_atol_k: 1.0e-6,
        no_claim_reason: CROSS_CODE_ONLY,
    },
    ThermalLevelBCase {
        id: "thermal-b-kt-nonlinear-slab-v1",
        title: "Temperature-dependent k(T) slab with Dirichlet base and Robin tip",
        mesh_counts: [10, 4, 4],
        mesh_extent: [0.05, 0.02, 0.02],
        material: ThermalLevelBMaterial::LinearKt {
            knots: &[(250.0, 8.0), (500.0, 18.0)],
        },
        source: ThermalLevelBSource::None,
        bcs: &[
            ThermalLevelBBc {
                name: "hot",
                axis: 0,
                at_max: false,
                kind: ThermalLevelBBcKind::Dirichlet { t_k: 450.0 },
            },
            ThermalLevelBBc {
                name: "cooled",
                axis: 0,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 40.0,
                    t_inf_k: 300.0,
                },
            },
        ],
        picard: Some(ThermalLevelBPicard {
            tol: 1.0e-13,
            max_iterations: 60,
        }),
        probes: &[[5, 2, 2], [10, 2, 2], [2, 1, 1], [8, 3, 0], [10, 0, 0]],
        acceptance_atol_k: 1.0e-5,
        no_claim_reason: CROSS_CODE_ONLY,
    },
    ThermalLevelBCase {
        id: "thermal-b-fin-film-v1",
        title: "Three-dimensional rectangular fin with film on every exposed face",
        mesh_counts: [16, 2, 4],
        mesh_extent: [0.04, 0.004, 0.01],
        material: ThermalLevelBMaterial::Isotropic { k: 170.0 },
        source: ThermalLevelBSource::None,
        bcs: &[
            ThermalLevelBBc {
                name: "base",
                axis: 0,
                at_max: false,
                kind: ThermalLevelBBcKind::Dirichlet { t_k: 350.0 },
            },
            ThermalLevelBBc {
                name: "tip",
                axis: 0,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 60.0,
                    t_inf_k: 295.0,
                },
            },
            ThermalLevelBBc {
                name: "y-lo",
                axis: 1,
                at_max: false,
                kind: ThermalLevelBBcKind::Robin {
                    h: 60.0,
                    t_inf_k: 295.0,
                },
            },
            ThermalLevelBBc {
                name: "y-hi",
                axis: 1,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 60.0,
                    t_inf_k: 295.0,
                },
            },
            ThermalLevelBBc {
                name: "z-lo",
                axis: 2,
                at_max: false,
                kind: ThermalLevelBBcKind::Robin {
                    h: 60.0,
                    t_inf_k: 295.0,
                },
            },
            ThermalLevelBBc {
                name: "z-hi",
                axis: 2,
                at_max: true,
                kind: ThermalLevelBBcKind::Robin {
                    h: 60.0,
                    t_inf_k: 295.0,
                },
            },
        ],
        picard: None,
        probes: &[[8, 1, 2], [16, 1, 2], [4, 0, 4], [12, 2, 0], [16, 2, 4]],
        acceptance_atol_k: 1.0e-6,
        no_claim_reason: CROSS_CODE_ONLY,
    },
];

/// Complete, stable Level-B thermal cross-code case catalog.
#[must_use]
pub fn thermal_level_b_cases() -> &'static [ThermalLevelBCase] {
    &THERMAL_LEVEL_B_CASES
}

/// The catalog case with the given id.
#[must_use]
pub fn thermal_level_b_case(case_id: &str) -> Option<&'static ThermalLevelBCase> {
    THERMAL_LEVEL_B_CASES.iter().find(|case| case.id == case_id)
}

/// One frozen probe row from the external reference manifest.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalLevelBProbe {
    /// Zero-based probe index, dense per case.
    pub index: usize,
    /// Vertex grid indices `(i, j, k)`.
    pub grid: [usize; 3],
    /// Node position recorded by the external code, m.
    pub position_m: [f64; 3],
    /// Frozen external probe temperature, K.
    pub temperature_k: f64,
}

/// One parsed per-case block of the external reference manifest.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalLevelBReference {
    /// Case id the block claims to solve.
    pub case_id: String,
    /// BLAKE3 of the committed deck bytes the external run consumed.
    pub deck_blake3: String,
    /// BLAKE3 of the canonical mesh bytes the external run assembled.
    pub mesh_blake3: String,
    /// External code identity string, e.g. `scikit-fem 11.0.0`.
    pub external_code: String,
    /// External linear-solver identity string.
    pub linear_solver: String,
    /// Picard iterations the external run needed (1 for linear cases).
    pub picard_iterations: u64,
    /// Minimum nodal temperature of the external field, K.
    pub t_min_k: f64,
    /// Maximum nodal temperature of the external field, K.
    pub t_max_k: f64,
    /// Every load-bearing spec value the external solver echoed back.
    pub spec_echo: BTreeMap<String, f64>,
    /// Frozen probe rows, dense by index.
    pub probes: Vec<ThermalLevelBProbe>,
}

/// Fail-closed manifest parsing and binding refusals.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermalLevelBError {
    /// A data line does not have the column count its kind requires.
    Columns {
        /// One-based manifest line number.
        line: usize,
        /// The row kind.
        kind: String,
    },
    /// A row kind outside `case_meta`/`spec_echo`/`probe`.
    UnknownKind {
        /// One-based manifest line number.
        line: usize,
        /// The unrecognized kind token.
        kind: String,
    },
    /// A numeric field failed strict parsing.
    Number {
        /// One-based manifest line number.
        line: usize,
        /// The offending field text.
        field: String,
    },
    /// A metadata or spec-echo key appeared twice for one case.
    Duplicate {
        /// Case id.
        case_id: String,
        /// Duplicated key.
        key: String,
    },
    /// A mandatory metadata key is missing for a case.
    MissingMeta {
        /// Case id.
        case_id: String,
        /// The absent key.
        key: String,
    },
    /// The external harness self-check did not record `pass`.
    SelfCheck {
        /// Case id.
        case_id: String,
        /// The recorded value.
        recorded: String,
    },
    /// Probe indices are not dense and ordered from zero.
    ProbeOrder {
        /// Case id.
        case_id: String,
        /// Expected next index.
        expected: usize,
        /// Observed index.
        observed: usize,
    },
    /// The manifest names a case this catalog does not define.
    UnknownCase {
        /// The unmatched case id.
        case_id: String,
    },
    /// A catalog case has no manifest block.
    MissingCase {
        /// The unreferenced case id.
        case_id: String,
    },
    /// The echoed spec-key set differs from the catalog expectation.
    EchoKeys {
        /// Case id.
        case_id: String,
        /// A key present on exactly one side.
        key: String,
    },
    /// An echoed spec value is not bit-identical to the catalog constant.
    EchoValue {
        /// Case id.
        case_id: String,
        /// Spec-echo key.
        key: String,
        /// Catalog value.
        expected: f64,
        /// Echoed value.
        observed: f64,
    },
    /// Probe rows disagree with the catalog probe list.
    ProbeGrid {
        /// Case id.
        case_id: String,
        /// Probe index.
        index: usize,
        /// Human-readable disagreement.
        detail: String,
    },
}

impl fmt::Display for ThermalLevelBError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Columns { line, kind } => {
                write!(
                    f,
                    "manifest line {line}: wrong column count for kind {kind}"
                )
            }
            Self::UnknownKind { line, kind } => {
                write!(f, "manifest line {line}: unknown row kind {kind}")
            }
            Self::Number { line, field } => {
                write!(f, "manifest line {line}: unparsable number {field:?}")
            }
            Self::Duplicate { case_id, key } => {
                write!(f, "case {case_id}: duplicate key {key}")
            }
            Self::MissingMeta { case_id, key } => {
                write!(f, "case {case_id}: missing mandatory metadata {key}")
            }
            Self::SelfCheck { case_id, recorded } => write!(
                f,
                "case {case_id}: external self-check recorded {recorded:?}, not \"pass\""
            ),
            Self::ProbeOrder {
                case_id,
                expected,
                observed,
            } => write!(
                f,
                "case {case_id}: probe index {observed} where {expected} was required"
            ),
            Self::UnknownCase { case_id } => {
                write!(f, "manifest names unknown case {case_id}")
            }
            Self::MissingCase { case_id } => {
                write!(f, "catalog case {case_id} has no manifest block")
            }
            Self::EchoKeys { case_id, key } => {
                write!(f, "case {case_id}: spec-echo key set mismatch at {key}")
            }
            Self::EchoValue {
                case_id,
                key,
                expected,
                observed,
            } => write!(
                f,
                "case {case_id}: spec echo {key} is {observed:?}, catalog holds {expected:?}"
            ),
            Self::ProbeGrid {
                case_id,
                index,
                detail,
            } => write!(f, "case {case_id} probe {index}: {detail}"),
        }
    }
}

const MANDATORY_META: [&str; 8] = [
    "deck_blake3",
    "external_code",
    "linear_solver",
    "mesh_blake3",
    "picard_iterations",
    "self_check",
    "t_max_k",
    "t_min_k",
];

fn parse_f64(line: usize, text: &str) -> Result<f64, ThermalLevelBError> {
    text.parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or(ThermalLevelBError::Number {
            line,
            field: text.to_string(),
        })
}

fn parse_usize(line: usize, text: &str) -> Result<usize, ThermalLevelBError> {
    text.parse::<usize>()
        .map_err(|_| ThermalLevelBError::Number {
            line,
            field: text.to_string(),
        })
}

#[derive(Default)]
struct ReferenceBuilder {
    meta: BTreeMap<String, String>,
    spec_echo: BTreeMap<String, f64>,
    probes: Vec<ThermalLevelBProbe>,
}

impl ReferenceBuilder {
    fn finish(mut self, case_id: String) -> Result<ThermalLevelBReference, ThermalLevelBError> {
        for key in MANDATORY_META {
            if !self.meta.contains_key(key) {
                return Err(ThermalLevelBError::MissingMeta {
                    case_id: case_id.clone(),
                    key: key.to_string(),
                });
            }
        }
        // Every MANDATORY_META key was just checked present.
        let mut take = |key: &str| self.meta.remove(key).expect("mandatory key checked above");
        let self_check = take("self_check");
        if self_check != "pass" {
            return Err(ThermalLevelBError::SelfCheck {
                case_id: case_id.clone(),
                recorded: self_check,
            });
        }
        let finite = |field: String| -> Result<f64, ThermalLevelBError> {
            field
                .parse::<f64>()
                .ok()
                .filter(|value| value.is_finite())
                .ok_or(ThermalLevelBError::Number { line: 0, field })
        };
        let picard_text = take("picard_iterations");
        let picard_iterations =
            picard_text
                .parse::<u64>()
                .map_err(|_| ThermalLevelBError::Number {
                    line: 0,
                    field: picard_text,
                })?;
        let t_min_k = finite(take("t_min_k"))?;
        let t_max_k = finite(take("t_max_k"))?;
        Ok(ThermalLevelBReference {
            deck_blake3: take("deck_blake3"),
            mesh_blake3: take("mesh_blake3"),
            external_code: take("external_code"),
            linear_solver: take("linear_solver"),
            picard_iterations,
            t_min_k,
            t_max_k,
            spec_echo: self.spec_echo,
            probes: self.probes,
            case_id,
        })
    }
}

/// Per-manifest accumulation state shared by the row handlers.
#[derive(Default)]
struct ManifestBuilder {
    order: Vec<String>,
    blocks: BTreeMap<String, ReferenceBuilder>,
}

impl ManifestBuilder {
    fn block(&mut self, case_id: &str) -> &mut ReferenceBuilder {
        if !self.blocks.contains_key(case_id) {
            self.order.push(case_id.to_string());
        }
        self.blocks.entry(case_id.to_string()).or_default()
    }

    fn keyed_row(
        &mut self,
        line: usize,
        columns: &[&str],
        is_echo: bool,
    ) -> Result<(), ThermalLevelBError> {
        if columns.len() != 4 {
            return Err(ThermalLevelBError::Columns {
                line,
                kind: columns[0].to_string(),
            });
        }
        let (case_id, key) = (columns[1], columns[2]);
        let duplicate = if is_echo {
            let value = parse_f64(line, columns[3])?;
            self.block(case_id)
                .spec_echo
                .insert(key.to_string(), value)
                .is_some()
        } else {
            self.block(case_id)
                .meta
                .insert(key.to_string(), columns[3].to_string())
                .is_some()
        };
        if duplicate {
            return Err(ThermalLevelBError::Duplicate {
                case_id: case_id.to_string(),
                key: key.to_string(),
            });
        }
        Ok(())
    }

    fn probe_row(&mut self, line: usize, columns: &[&str]) -> Result<(), ThermalLevelBError> {
        if columns.len() != 10 {
            return Err(ThermalLevelBError::Columns {
                line,
                kind: columns[0].to_string(),
            });
        }
        let case_id = columns[1];
        let probe_index = parse_usize(line, columns[2])?;
        let grid = [
            parse_usize(line, columns[3])?,
            parse_usize(line, columns[4])?,
            parse_usize(line, columns[5])?,
        ];
        let position_m = [
            parse_f64(line, columns[6])?,
            parse_f64(line, columns[7])?,
            parse_f64(line, columns[8])?,
        ];
        let temperature_k = parse_f64(line, columns[9])?;
        let block = self.block(case_id);
        if block.probes.len() != probe_index {
            return Err(ThermalLevelBError::ProbeOrder {
                case_id: case_id.to_string(),
                expected: block.probes.len(),
                observed: probe_index,
            });
        }
        block.probes.push(ThermalLevelBProbe {
            index: probe_index,
            grid,
            position_m,
            temperature_k,
        });
        Ok(())
    }
}

/// Parse manifest bytes into per-case reference blocks, fail-closed.
///
/// # Errors
/// Every malformed line, duplicate or missing key, failed external
/// self-check, and unknown/unreferenced case id is a typed refusal.
pub fn parse_thermal_level_b_manifest(
    bytes: &[u8],
) -> Result<Vec<ThermalLevelBReference>, ThermalLevelBError> {
    let text = String::from_utf8_lossy(bytes);
    let mut builder = ManifestBuilder::default();
    for (index, raw) in text.lines().enumerate() {
        let line = index + 1;
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = raw.split('\t').collect();
        match columns[0] {
            "case_meta" => builder.keyed_row(line, &columns, false)?,
            "spec_echo" => builder.keyed_row(line, &columns, true)?,
            "probe" => builder.probe_row(line, &columns)?,
            other => {
                return Err(ThermalLevelBError::UnknownKind {
                    line,
                    kind: other.to_string(),
                });
            }
        }
    }
    let mut references = Vec::with_capacity(builder.order.len());
    for case_id in builder.order {
        let block = builder
            .blocks
            .remove(&case_id)
            .expect("every ordered id was inserted");
        references.push(block.finish(case_id)?);
    }
    Ok(references)
}

/// Every load-bearing spec value the external run must echo, mirroring
/// `tools/vvref/solve_skfem.py::spec_echo_rows` exactly.
#[must_use]
pub fn expected_spec_echo(case: &ThermalLevelBCase) -> BTreeMap<String, f64> {
    fn as_f64(value: usize) -> f64 {
        value as f64
    }
    let mut rows = BTreeMap::new();
    for (axis, name) in ["x", "y", "z"].iter().enumerate() {
        rows.insert(
            format!("mesh.counts.{name}"),
            as_f64(case.mesh_counts[axis]),
        );
        rows.insert(format!("mesh.extent.{name}"), case.mesh_extent[axis]);
    }
    match case.material {
        ThermalLevelBMaterial::Isotropic { k } => {
            rows.insert("material.isotropic.k".to_string(), k);
        }
        ThermalLevelBMaterial::Tensor { k } => {
            for (i, row) in k.iter().enumerate() {
                for (j, value) in row.iter().enumerate() {
                    rows.insert(format!("material.tensor.k{i}{j}"), *value);
                }
            }
        }
        ThermalLevelBMaterial::LinearKt { knots } => {
            for (n, (t, k)) in knots.iter().enumerate() {
                rows.insert(format!("material.linear_kt.knot{n}.t"), *t);
                rows.insert(format!("material.linear_kt.knot{n}.k"), *k);
            }
            let picard = case
                .picard
                .expect("catalog invariant: LinearKt cases declare picard controls");
            rows.insert("picard.tol".to_string(), picard.tol);
            rows.insert(
                "picard.max_iterations".to_string(),
                picard.max_iterations as f64,
            );
        }
    }
    if let ThermalLevelBSource::PolyXy { q0 } = case.source {
        rows.insert("source.poly_xy.q0".to_string(), q0);
    }
    for bc in case.bcs {
        let prefix = format!("bc.{}", bc.name);
        rows.insert(format!("{prefix}.axis"), as_f64(bc.axis));
        let target = if bc.at_max {
            case.mesh_extent[bc.axis]
        } else {
            0.0
        };
        rows.insert(format!("{prefix}.target"), target);
        match bc.kind {
            ThermalLevelBBcKind::Dirichlet { t_k } => {
                rows.insert(format!("{prefix}.dirichlet.t"), t_k);
            }
            ThermalLevelBBcKind::Robin { h, t_inf_k } => {
                rows.insert(format!("{prefix}.robin.h"), h);
                rows.insert(format!("{prefix}.robin.t_inf"), t_inf_k);
            }
        }
    }
    rows
}

/// Require the echoed spec to match the catalog constants bit-exactly.
///
/// # Errors
/// [`ThermalLevelBError::EchoKeys`] on any key-set difference and
/// [`ThermalLevelBError::EchoValue`] on any non-bit-identical value.
pub fn verify_spec_echo(
    case: &ThermalLevelBCase,
    reference: &ThermalLevelBReference,
) -> Result<(), ThermalLevelBError> {
    let expected = expected_spec_echo(case);
    for key in expected.keys() {
        if !reference.spec_echo.contains_key(key) {
            return Err(ThermalLevelBError::EchoKeys {
                case_id: case.id.to_string(),
                key: key.clone(),
            });
        }
    }
    for key in reference.spec_echo.keys() {
        if !expected.contains_key(key) {
            return Err(ThermalLevelBError::EchoKeys {
                case_id: case.id.to_string(),
                key: key.clone(),
            });
        }
    }
    for (key, expected_value) in &expected {
        let observed = reference.spec_echo[key];
        if expected_value.to_bits() != observed.to_bits() {
            return Err(ThermalLevelBError::EchoValue {
                case_id: case.id.to_string(),
                key: key.clone(),
                expected: *expected_value,
                observed,
            });
        }
    }
    Ok(())
}

/// Require every frozen probe row to sit exactly on the catalog probe
/// list, with positions bit-identical to the fixture arithmetic
/// `index * (extent / count)`.
///
/// # Errors
/// [`ThermalLevelBError::ProbeGrid`] naming the first disagreement.
pub fn verify_probe_grid(
    case: &ThermalLevelBCase,
    reference: &ThermalLevelBReference,
) -> Result<(), ThermalLevelBError> {
    let fail = |index: usize, detail: String| ThermalLevelBError::ProbeGrid {
        case_id: case.id.to_string(),
        index,
        detail,
    };
    if reference.probes.len() != case.probes.len() {
        return Err(fail(
            reference.probes.len(),
            format!(
                "manifest froze {} probes, catalog declares {}",
                reference.probes.len(),
                case.probes.len()
            ),
        ));
    }
    for (probe, declared) in reference.probes.iter().zip(case.probes) {
        if probe.grid != *declared {
            return Err(fail(
                probe.index,
                format!("grid {:?} differs from catalog {declared:?}", probe.grid),
            ));
        }
        for axis in 0..3 {
            if probe.grid[axis] > case.mesh_counts[axis] {
                return Err(fail(
                    probe.index,
                    format!(
                        "grid index {} exceeds vertex range 0..={}",
                        probe.grid[axis], case.mesh_counts[axis]
                    ),
                ));
            }
            let spacing = case.mesh_extent[axis] / case.mesh_counts[axis] as f64;
            let expected = probe.grid[axis] as f64 * spacing;
            if expected.to_bits() != probe.position_m[axis].to_bits() {
                return Err(fail(
                    probe.index,
                    format!(
                        "axis {axis} position {:?} is not bit-identical to fixture arithmetic {expected:?}",
                        probe.position_m[axis]
                    ),
                ));
            }
        }
        if !(reference.t_min_k..=reference.t_max_k).contains(&probe.temperature_k) {
            return Err(fail(
                probe.index,
                format!(
                    "temperature {} outside the recorded field range [{}, {}]",
                    probe.temperature_k, reference.t_min_k, reference.t_max_k
                ),
            ));
        }
    }
    Ok(())
}

/// Fully verified references: parsed fail-closed, every catalog case
/// present, no unknown cases, spec echo and probe grid bound bit-exactly.
///
/// # Errors
/// The first refusal from parsing or binding, by case order.
pub fn verified_thermal_level_b_references()
-> Result<Vec<ThermalLevelBReference>, ThermalLevelBError> {
    let references = parse_thermal_level_b_manifest(THERMAL_LEVEL_B_MANIFEST)?;
    for reference in &references {
        let case = thermal_level_b_case(&reference.case_id).ok_or_else(|| {
            ThermalLevelBError::UnknownCase {
                case_id: reference.case_id.clone(),
            }
        })?;
        verify_spec_echo(case, reference)?;
        verify_probe_grid(case, reference)?;
    }
    for case in thermal_level_b_cases() {
        if !references.iter().any(|r| r.case_id == case.id) {
            return Err(ThermalLevelBError::MissingCase {
                case_id: case.id.to_string(),
            });
        }
    }
    Ok(references)
}

static VERIFIED: LazyLock<Result<Vec<ThermalLevelBReference>, ThermalLevelBError>> =
    LazyLock::new(verified_thermal_level_b_references);

/// Shared verified references.
///
/// # Errors
/// The cached refusal when the committed manifest fails verification.
pub fn thermal_level_b_references()
-> Result<&'static [ThermalLevelBReference], &'static ThermalLevelBError> {
    match &*VERIFIED {
        Ok(references) => Ok(references.as_slice()),
        Err(error) => Err(error),
    }
}

/// The verified reference block for one case id.
///
/// # Errors
/// The cached manifest refusal; a missing id yields `Ok(None)` only for
/// ids outside the catalog (catalog members are guaranteed present by
/// verification).
pub fn thermal_level_b_reference(
    case_id: &str,
) -> Result<Option<&'static ThermalLevelBReference>, &'static ThermalLevelBError> {
    Ok(thermal_level_b_references()?
        .iter()
        .find(|reference| reference.case_id == case_id))
}
