//! Acoustic/music validation-corpus rows (music bead
//! `frankensim-music-v8-root-3ez8g.1.1`, program root
//! `frankensim-music-v8-root-3ez8g`).
//!
//! The music program's claims registry (`instrument-claims.json`) may cite
//! corpora only through rows registered HERE, under the same discipline as
//! the thermal Level-A/B/C catalogs: licensing resolved BEFORE ingestion,
//! two-source verification for transcribed values, provenance named per row,
//! and absence recorded as a named hunt rather than silence. Every
//! value-bearing row stays reference-only until a consuming kernel binds it
//! in a test — registration is not validation, and no corpus query
//! manufactures a physical claim.
//!
//! Licensing decisions recorded here (the bead's fourth DONE-WHEN):
//!
//! - **Ernoult et al. 2021** (Acta Acustica 5:47): the paper is CC-BY-4.0
//!   and its Table 1 geometry is retained verbatim. The measured impedance
//!   CURVES are published through openwind under GPLv3 and are NOT retained
//!   in this MIT+rider repository (see [`acoustic_absences`] row
//!   `acoustic-refused-openwind-curves`). The five scalar first-peak
//!   frequencies are retained as facts with dual citation; the three
//!   measurement sessions agree within ~3 Hz, and the paper states ±2-cent
//!   peak accuracy.
//! - **Olson & Hazell 1977** values are transcribed from TWO independent
//!   secondary sources (Srivastava/Datta/Sheikh 2004 Table 4;
//!   Thinh/Binh/Tu 2013 Table 1) per the two-source transcription law; the
//!   committed `fs-plate` literature test is the retained transcription.
//! - **Carcagno et al. 2018** (JASA 144(6):3533) is CC-BY; the tabulated
//!   BR-rosewood F/Q values were pdftotext-verified in-session and the
//!   committed `fs-modalid` benchmark test is the retained transcription.
//!   Mobility MAGNITUDES are figures-only across the CC-BY literature and
//!   are refused (never digitized from figures).
//! - **Analytic rows** (free-free bar ratios, Leissa's clamped-square
//!   eigenvalue) are formula-documented definitions; the bar ratios are
//!   re-derivable in-fixture from the cosh·cos = 1 characteristic equation
//!   and need no transcription trust at all.

use fs_qty::Dims;

/// Tracked tab-separated manifest backing every acoustic corpus row. The
/// conformance test renders the expected bytes FROM [`acoustic_cases`] and
/// [`acoustic_absences`] and byte-compares, so the tracked file cannot
/// drift from the code catalog.
pub const ACOUSTIC_MANIFEST_LOCATOR: &str = "data/vv-corpus/acoustic/acoustic-v1.tsv";

/// Coherent-SI dimensions of frequency (s^-1).
pub const FREQUENCY_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);
/// Coherent-SI dimensions of length.
pub const LENGTH_DIMS: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// Coverage family for acoustic rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AcousticFamily {
    /// Wind-instrument fingering ladder (measured first impedance peaks).
    FingeringLadder,
    /// Plate/panel modal references (analytic and literature).
    PlateReference,
    /// Coupled guitar-body modal parameters.
    GuitarBody,
    /// Analytic bar/beam definitions.
    BarAnalytic,
}

impl AcousticFamily {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::FingeringLadder => "fingering-ladder",
            Self::PlateReference => "plate-reference",
            Self::GuitarBody => "guitar-body",
            Self::BarAnalytic => "bar-analytic",
        }
    }
}

/// Evidence level of the row's reference value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AcousticLevel {
    /// Analytic/textbook definition (Level-A analog): re-derivable, no
    /// transcription trust required beyond the cited formula.
    AnalyticDefinition,
    /// Published-experiment record (Level-C analog): derived/tabulated
    /// values with named licensing and retention decisions.
    PublishedExperiment,
}

impl AcousticLevel {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::AnalyticDefinition => "analytic-definition",
            Self::PublishedExperiment => "published-experiment",
        }
    }
}

/// Acceptance envelope a consuming binding must satisfy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AcousticAcceptance {
    /// Pitch-class envelope in cents: `outer` is the hard gate, `inner`
    /// the post-refinement band a modern model is expected to hold.
    CentsEnvelope {
        /// Hard gate, cents.
        outer: f64,
        /// Expected refined band, cents.
        inner: f64,
    },
    /// Plain relative tolerance on the metric.
    Relative {
        /// `|model - reference| / |reference|` ceiling.
        rtol: f64,
    },
}

impl AcousticAcceptance {
    /// Stable manifest spelling.
    #[must_use]
    pub fn render(self) -> String {
        match self {
            Self::CentsEnvelope { outer, inner } => {
                format!("cents-outer={outer},cents-inner={inner}")
            }
            Self::Relative { rtol } => format!("rtol={rtol}"),
        }
    }
}

/// One inclusive query-context range (point ranges pin exact geometry).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcousticContext {
    /// Stable axis name.
    pub name: &'static str,
    /// Coherent-SI dimensions.
    pub dims: Dims,
    /// Inclusive lower bound.
    pub lo: f64,
    /// Inclusive upper bound.
    pub hi: f64,
}

/// One value-bearing acoustic corpus row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcousticCase {
    /// Stable corpus row id (`vvreg:` prefix in registry refs).
    pub id: &'static str,
    /// Human-readable title.
    pub title: &'static str,
    /// Coverage family.
    pub family: AcousticFamily,
    /// Evidence level.
    pub level: AcousticLevel,
    /// Stable scalar metric name.
    pub metric: &'static str,
    /// Metric dimensions.
    pub metric_dims: Dims,
    /// Frozen reference value in coherent SI.
    pub reference_value_si: f64,
    /// Exact formula/measurement semantics.
    pub formula: &'static str,
    /// Acceptance envelope for consuming bindings.
    pub acceptance: AcousticAcceptance,
    /// Complete query context (geometry as point ranges where exact).
    pub context: &'static [AcousticContext],
    /// Citation of the source of the value.
    pub source: &'static str,
    /// License of the source and of what is retained.
    pub license: &'static str,
    /// Exactly what bytes/values this repository retains.
    pub retention: &'static str,
    /// Why this row carries no solver evidence by itself.
    pub no_claim_reason: &'static str,
}

/// One recorded absence or refusal — the population signal, named.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcousticAbsence {
    /// Stable row id.
    pub id: &'static str,
    /// Human-readable title.
    pub title: &'static str,
    /// What is missing or refused.
    pub what: &'static str,
    /// Why (license refusal, no compatible source found, not yet hunted).
    pub why: &'static str,
    /// Which work the row unblocks when filled.
    pub unblocks: &'static str,
}

const REFERENCE_ONLY: &str =
    "reference row only; a consuming kernel test must bind it before any gate cites it";

const fn ctx(name: &'static str, dims: Dims, value: f64) -> AcousticContext {
    AcousticContext {
        name,
        dims,
        lo: value,
        hi: value,
    }
}

/// Ernoult 2021 Table 1 caliper-measured four-hole cylinder geometry
/// (CC-BY-4.0). Hole positions are distances from the input plane; the
/// bore continues 47.5 mm past the last hole to L = 287.5 mm.
pub const ERNOULT_GEOMETRY: &[AcousticContext] = &[
    ctx("bore-radius-m", LENGTH_DIMS, 2.0e-3),
    ctx("bore-length-m", LENGTH_DIMS, 0.2875),
    ctx("hole-1-position-m", LENGTH_DIMS, 0.100),
    ctx("hole-1-radius-m", LENGTH_DIMS, 1.5e-3),
    ctx("hole-1-chimney-m", LENGTH_DIMS, 1.7e-3),
    ctx("hole-2-position-m", LENGTH_DIMS, 0.130),
    ctx("hole-2-radius-m", LENGTH_DIMS, 1.75e-3),
    ctx("hole-2-chimney-m", LENGTH_DIMS, 1.3e-3),
    ctx("hole-3-position-m", LENGTH_DIMS, 0.180),
    ctx("hole-3-radius-m", LENGTH_DIMS, 1.75e-3),
    ctx("hole-3-chimney-m", LENGTH_DIMS, 1.5e-3),
    ctx("hole-4-position-m", LENGTH_DIMS, 0.240),
    ctx("hole-4-radius-m", LENGTH_DIMS, 1.25e-3),
    ctx("hole-4-chimney-m", LENGTH_DIMS, 1.4e-3),
];

const ERNOULT_SOURCE: &str = "Ernoult, Chabassier, Rodriguez, Humeau, Acta Acustica 5:47 (2021), \
     Table 1 geometry; measured first-peak values from the paper's two-microphone method \
     (sessions Measure1-3 agree within ~3 Hz; stated peak accuracy +-2 cents)";
const ERNOULT_LICENSE: &str = "paper CC-BY-4.0; measured curves openwind GPLv3 (curves NOT \
     retained - see acoustic-refused-openwind-curves); scalar peak values retained as facts \
     with dual citation";
const ERNOULT_RETENTION: &str =
    "Table 1 geometry + five scalar first-peak frequencies; no curve bytes";
const ERNOULT_ACCEPT: AcousticAcceptance = AcousticAcceptance::CentsEnvelope {
    outer: 30.0,
    inner: 10.0,
};

const fn ernoult(
    id: &'static str,
    title: &'static str,
    measured_hz: f64,
    formula: &'static str,
) -> AcousticCase {
    AcousticCase {
        id,
        title,
        family: AcousticFamily::FingeringLadder,
        level: AcousticLevel::PublishedExperiment,
        metric: "first-impedance-peak-frequency",
        metric_dims: FREQUENCY_DIMS,
        reference_value_si: measured_hz,
        formula,
        acceptance: ERNOULT_ACCEPT,
        context: ERNOULT_GEOMETRY,
        source: ERNOULT_SOURCE,
        license: ERNOULT_LICENSE,
        retention: ERNOULT_RETENTION,
        no_claim_reason: REFERENCE_ONLY,
    }
}

const OH_SOURCE: &str = "Olson & Hazell (1977) clamped stiffened panel; values reproduced in \
     Srivastava, Datta & Sheikh, Shock and Vibration 11 (2004) 9-19, Table 4 AND Thinh, Binh \
     & Tu (2013), Table 1 (two independent secondary sources)";
const OH_LICENSE: &str = "factual frequency values transcribed under the two-source law; the \
     committed fs-plate literature test is the retained transcription";
const OH_RETENTION: &str = "five theory frequencies (Hz); experiment column noted in formulas";

const fn olson_hazell(
    id: &'static str,
    title: &'static str,
    theory_hz: f64,
    formula: &'static str,
) -> AcousticCase {
    AcousticCase {
        id,
        title,
        family: AcousticFamily::PlateReference,
        level: AcousticLevel::PublishedExperiment,
        metric: "panel-mode-frequency",
        metric_dims: FREQUENCY_DIMS,
        reference_value_si: theory_hz,
        formula,
        acceptance: AcousticAcceptance::Relative { rtol: 0.08 },
        context: &[],
        source: OH_SOURCE,
        license: OH_LICENSE,
        retention: OH_RETENTION,
        no_claim_reason: REFERENCE_ONLY,
    }
}

const CARCAGNO_SOURCE: &str = "Carcagno, Bucknall, Woodhouse, Fritz, Plack, JASA 144(6):3533 \
     (2018), tabulated modal parameters of the BR-rosewood guitar";
const CARCAGNO_LICENSE: &str = "CC-BY; F/Q values pdftotext-verified in-session; the committed \
     fs-modalid published-parameter benchmark is the retained transcription; mobility \
     MAGNITUDES are figures-only across the CC-BY literature and are refused";
const CARCAGNO_RETENTION: &str = "three (F, Q) scalar pairs; no figure digitization";

const fn carcagno(
    id: &'static str,
    title: &'static str,
    frequency_hz: f64,
    q_context: &'static [AcousticContext],
    formula: &'static str,
) -> AcousticCase {
    AcousticCase {
        id,
        title,
        family: AcousticFamily::GuitarBody,
        level: AcousticLevel::PublishedExperiment,
        metric: "body-mode-frequency",
        metric_dims: FREQUENCY_DIMS,
        reference_value_si: frequency_hz,
        formula,
        acceptance: AcousticAcceptance::Relative { rtol: 0.10 },
        context: q_context,
        source: CARCAGNO_SOURCE,
        license: CARCAGNO_LICENSE,
        retention: CARCAGNO_RETENTION,
        no_claim_reason: REFERENCE_ONLY,
    }
}

const CARCAGNO_Q1: &[AcousticContext] = &[ctx("quality-factor", Dims::NONE, 34.0)];
const CARCAGNO_Q2: &[AcousticContext] = &[ctx("quality-factor", Dims::NONE, 18.0)];
const CARCAGNO_Q3: &[AcousticContext] = &[ctx("quality-factor", Dims::NONE, 36.0)];

/// Free-free Euler-Bernoulli bar frequency ratios from the cosh·cos = 1
/// characteristic roots (betaL = 4.730040745..., 7.853204624...,
/// 10.995607838...). Ratios are (beta_n/beta_1)^2; consumers re-derive the
/// roots in-fixture by deterministic Newton (self-verified — the analytic
/// row exists so the registry can NAME the definition, not to be trusted
/// blindly).
const BAR_SOURCE: &str = "free-free Euler-Bernoulli characteristic equation cosh(bL)cos(bL)=1; \
     roots re-derivable in-fixture (self-verified pin technique)";
const BAR_LICENSE: &str = "analytic; no transcription trust required";
const BAR_RETENTION: &str = "ratio values documented; consumers recompute roots";

const fn bar_ratio(
    id: &'static str,
    title: &'static str,
    ratio: f64,
    formula: &'static str,
) -> AcousticCase {
    AcousticCase {
        id,
        title,
        family: AcousticFamily::BarAnalytic,
        level: AcousticLevel::AnalyticDefinition,
        metric: "partial-frequency-ratio",
        metric_dims: Dims::NONE,
        reference_value_si: ratio,
        formula,
        acceptance: AcousticAcceptance::Relative { rtol: 1.0e-3 },
        context: &[],
        source: BAR_SOURCE,
        license: BAR_LICENSE,
        retention: BAR_RETENTION,
        no_claim_reason: REFERENCE_ONLY,
    }
}

const CASES: &[AcousticCase] = &[
    ernoult(
        "acoustic-ernoult-2021-xxxx",
        "Ernoult 2021 fingering xxxx (all holes closed)",
        283.0,
        "measured first impedance peak, all four holes closed",
    ),
    ernoult(
        "acoustic-ernoult-2021-xxxo",
        "Ernoult 2021 fingering xxxo (hole 4 open)",
        332.0,
        "measured first impedance peak, hole 4 open",
    ),
    ernoult(
        "acoustic-ernoult-2021-xxox",
        "Ernoult 2021 fingering xxox (hole 3 open)",
        449.0,
        "measured first impedance peak, hole 3 open",
    ),
    ernoult(
        "acoustic-ernoult-2021-xoxx",
        "Ernoult 2021 fingering xoxx (hole 2 open)",
        619.0,
        "measured first impedance peak, hole 2 open",
    ),
    ernoult(
        "acoustic-ernoult-2021-oxxx",
        "Ernoult 2021 fingering oxxx (hole 1 open)",
        770.0,
        "measured first impedance peak, hole 1 open",
    ),
    AcousticCase {
        id: "acoustic-leissa-clamped-square-lambda1",
        title: "Leissa clamped-square fundamental eigenvalue",
        family: AcousticFamily::PlateReference,
        level: AcousticLevel::AnalyticDefinition,
        metric: "clamped-square-lambda1",
        metric_dims: Dims::NONE,
        reference_value_si: 35.992,
        formula: "lambda = omega a^2 sqrt(rho h / D) = 35.992 (Leissa, NASA SP-160); next \
                  clamped mode at 73.41/35.99 ~ 2.04x in omega",
        acceptance: AcousticAcceptance::Relative { rtol: 0.02 },
        context: &[],
        source: "Leissa, Vibration of Plates, NASA SP-160 (1969), clamped-square table",
        license: "US government publication; tabulated eigenvalue is a factual constant",
        retention: "one dimensionless eigenvalue",
        no_claim_reason: REFERENCE_ONLY,
    },
    olson_hazell(
        "acoustic-olson-hazell-1977-mode1",
        "Olson-Hazell stiffened panel mode 1",
        718.1,
        "theory 718.1 Hz; experiment 689 Hz",
    ),
    olson_hazell(
        "acoustic-olson-hazell-1977-mode2",
        "Olson-Hazell stiffened panel mode 2",
        751.4,
        "theory 751.4 Hz; experiment 725 Hz",
    ),
    olson_hazell(
        "acoustic-olson-hazell-1977-mode3",
        "Olson-Hazell stiffened panel mode 3",
        997.4,
        "theory 997.4 Hz; experiment 961 Hz",
    ),
    olson_hazell(
        "acoustic-olson-hazell-1977-mode4",
        "Olson-Hazell stiffened panel mode 4",
        1007.1,
        "theory 1007.1 Hz; experiment 986 Hz",
    ),
    olson_hazell(
        "acoustic-olson-hazell-1977-mode5",
        "Olson-Hazell stiffened panel mode 5",
        1419.8,
        "theory 1419.8 Hz; experiment 1376 Hz",
    ),
    carcagno(
        "acoustic-carcagno-2018-mode1",
        "Carcagno 2018 BR-rosewood guitar body mode 1",
        97.0,
        CARCAGNO_Q1,
        "F1 = 97 Hz, Q1 = 34 (tabulated)",
    ),
    carcagno(
        "acoustic-carcagno-2018-mode2",
        "Carcagno 2018 BR-rosewood guitar body mode 2",
        177.0,
        CARCAGNO_Q2,
        "F2 = 177 Hz, Q2 = 18 (tabulated)",
    ),
    carcagno(
        "acoustic-carcagno-2018-mode3",
        "Carcagno 2018 BR-rosewood guitar body mode 3",
        336.0,
        CARCAGNO_Q3,
        "F3 = 336 Hz, Q3 = 36 (tabulated)",
    ),
    bar_ratio(
        "acoustic-bar-free-free-f2f1",
        "Free-free bar partial ratio f2/f1",
        2.756_538_507_099_962,
        "(7.853204624095838/4.730040744862704)^2",
    ),
    bar_ratio(
        "acoustic-bar-free-free-f3f1",
        "Free-free bar partial ratio f3/f1",
        5.403_917_632_383_322_5,
        "(10.995607838001671/4.730040744862704)^2",
    ),
];

const ABSENCES: &[AcousticAbsence] = &[
    AcousticAbsence {
        id: "acoustic-refused-openwind-curves",
        title: "openwind measured impedance curves",
        what: "full measured input-impedance curves for the Ernoult 2021 fingerings",
        why: "GPLv3; retention refused in this MIT+rider repository - scalar peak facts \
              retained instead (see the acoustic-ernoult-2021-* rows)",
        unblocks: "nothing (deliberate refusal, not a gap); curve-level comparisons would \
                   need an independently licensed measurement",
    },
    AcousticAbsence {
        id: "acoustic-absent-felt-coupon",
        title: "wool-felt force-compression coupon",
        what: "cited felt loading/unloading force-compression data with measurement protocol",
        why: "not yet hunted under the PD/CC-BY two-source recipe",
        unblocks: "frankensim-music-t-piano-felt-87zbd (the Uniaxial fit) and the piano \
                   tilt-vs-velocity gate (3ez8g.5.3)",
    },
    AcousticAbsence {
        id: "acoustic-absent-glottal-lf",
        title: "published glottal-flow waveform parameters",
        what: "LF-model-class glottal flow parameter sets with license",
        why: "not yet hunted; parameters must register as model-card rows, distinct from \
              constitutive tissue data",
        unblocks: "vowel gates (3ez8g.8.3) source-shape QoIs",
    },
    AcousticAbsence {
        id: "acoustic-absent-measured-frf-trace",
        title: "license-compatible measured FRF trace",
        what: "a published, license-compatible instrument FRF TRACE (not just parameters)",
        why: "fs-modalid's recorded no-claim: figures-only across the CC-BY literature; the \
              published-parameter round-trip remains the honest benchmark",
        unblocks: "fs-modalid trace-level benchmarks and simulate-vs-measure board gates \
                   (3ez8g.5.3, 3ez8g.7.4)",
    },
];

/// Complete, stable value-bearing acoustic catalog.
#[must_use]
pub fn acoustic_cases() -> &'static [AcousticCase] {
    CASES
}

/// Recorded absences and refusals (the population signal).
#[must_use]
pub fn acoustic_absences() -> &'static [AcousticAbsence] {
    ABSENCES
}

/// Render the canonical manifest bytes from the catalog. The tracked TSV at
/// [`ACOUSTIC_MANIFEST_LOCATOR`] must byte-equal this rendering (enforced by
/// the conformance test), so the file is a projection that cannot drift.
#[must_use]
pub fn render_acoustic_manifest() -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    out.push_str(
        "schema_version\tcase_id\tfamily\tlevel\tmetric\treference_value_si\tacceptance\t\
         formula\tsource\tlicense\tretention\tstatus\n",
    );
    for case in CASES {
        let _ = writeln!(
            out,
            "1\t{}\t{}\t{}\t{}\t{:e}\t{}\t{}\t{}\t{}\t{}\treference-only",
            case.id,
            case.family.name(),
            case.level.name(),
            case.metric,
            case.reference_value_si,
            case.acceptance.render(),
            case.formula,
            case.source,
            case.license,
            case.retention,
        );
    }
    for absence in ABSENCES {
        let _ = writeln!(
            out,
            "1\t{}\tabsence\tabsence\tnone\t\t\t{}\t{}\t{}\t\t{}",
            absence.id,
            absence.what,
            absence.why,
            absence.unblocks,
            if absence.id.contains("refused") {
                "refused-retention"
            } else {
                "absent-hunt"
            },
        );
    }
    out
}
