//! Valve/crook delta-L extraction (music bead
//! `frankensim-music-v8-root-3ez8g.3.2`): each valve crook's INSERTED
//! LENGTH is minted from its lumen geometry by the bore extractor
//! ([`crate::extract_bore`] — this module applies that machinery
//! per-crook), certified against the CAD's own declared centerline arc
//! length, and emitted as a control-rate chart delta the performance
//! images consume. A crook that measures 10% long is a DATA ERROR, not
//! a detuned trumpet: the band certificate REFUSES by name instead of
//! silently detuning the horn.
//!
//! Junction honesty: the crook's port areas are checked against the
//! main bore's area. A mismatch beyond band is EITHER a real acoustic
//! feature (a deliberate step the CAD intent declares) or an extraction
//! error — the receipt distinguishes by comparing against the declared
//! intent, refusing only the undeclared mismatch.
//!
//! Authority: everything here inherits the bore extractor's Estimate
//! authority; the band certificate is a CONSISTENCY certificate against
//! CAD intent, not a measurement of a physical horn.

use fs_exec::Cx;
use fs_geom::Chart;
use fs_rep_mesh::Soup;

use crate::{BoreConfig, BoreError, BoreExtraction, extract_bore};

/// CAD-declared intent for one crook.
#[derive(Debug, Clone, Copy)]
pub struct CrookCadIntent {
    /// Declared centerline arc length of the inserted tube [m].
    pub centerline_length_m: f64,
    /// Declared port radius where the crook meets the main bore [m].
    pub port_radius_m: f64,
    /// Declared DELIBERATE port-area step ratio (crook area / main bore
    /// area); 1.0 = flush junction intended.
    pub declared_area_step: f64,
}

/// Configuration for one crook extraction.
#[derive(Debug, Clone)]
pub struct CrookConfig {
    /// Bore-extractor knobs for the lumen.
    pub bore: BoreConfig,
    /// Refuse when `|extracted - cad| / cad` exceeds this.
    pub length_band: f64,
    /// Refuse when the port area deviates from the DECLARED step by more
    /// than this relative band.
    pub junction_band: f64,
}

/// One valve's minted chart delta — the control-rate object a
/// performance image consumes (the numbers; the D17 switching lift
/// lives with the runtime, per the bead's own boundary).
#[derive(Debug, Clone)]
pub struct ValveChartDelta {
    /// Caller label (e.g. `valve-1`).
    pub label: String,
    /// Extracted inserted length [m].
    pub delta_l_m: f64,
    /// CAD-declared length [m] (the certificate's reference).
    pub cad_length_m: f64,
    /// `|extracted - cad| / cad` — inside the band by construction.
    pub length_deviation: f64,
    /// Mean extracted port area over the two end stations [m^2].
    pub port_area_m2: f64,
    /// Extracted port-area step ratio vs the main bore.
    pub area_step: f64,
    /// The full bore receipt (provenance chain: boundary digest,
    /// stations, closure).
    pub bore: BoreExtraction,
}

/// A continuous slide's certified range (tuning slide, trombone slide):
/// two certified endpoint extractions; the runtime interpolates between
/// them (the sub-sample delay question belongs to the wind epic's
/// fractional-delay bead — a dependency of claim, not of code).
#[derive(Debug, Clone)]
pub struct SlideRange {
    /// Label (e.g. `tuning-slide`).
    pub label: String,
    /// Shortest inserted length [m].
    pub min_delta_l_m: f64,
    /// Longest inserted length [m].
    pub max_delta_l_m: f64,
}

/// Typed refusals from crook extraction.
#[derive(Debug)]
pub enum CrookError {
    /// The underlying bore extraction refused.
    Bore(BoreError),
    /// A CAD-intent or config parameter is unusable.
    Invalid {
        /// Diagnosis.
        what: &'static str,
    },
    /// The extracted length disagrees with the CAD centerline beyond the
    /// authored band — a data error, refused by name.
    DeltaLengthOutOfBand {
        /// Extracted length [m].
        extracted_m: f64,
        /// CAD-declared length [m].
        cad_m: f64,
        /// Measured relative deviation.
        deviation: f64,
        /// The authored band.
        band: f64,
    },
    /// The port area deviates from the DECLARED junction step beyond
    /// band — an undeclared mismatch (extraction error or undocumented
    /// geometry), refused rather than silently voicing a step.
    JunctionAreaMismatch {
        /// Extracted step ratio (crook area / main bore area).
        extracted_step: f64,
        /// Declared step ratio.
        declared_step: f64,
        /// The authored band.
        band: f64,
    },
}

impl core::fmt::Display for CrookError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CrookError::Bore(e) => write!(f, "FS-QUERY-CROOK bore: {e}"),
            CrookError::Invalid { what } => write!(f, "FS-QUERY-CROOK: {what}"),
            CrookError::DeltaLengthOutOfBand {
                extracted_m,
                cad_m,
                deviation,
                band,
            } => write!(
                f,
                "FS-QUERY-CROOK-LENGTH-BAND: extracted {extracted_m:.6} m vs CAD \
                 {cad_m:.6} m (deviation {deviation:.4} > band {band:.4}) — a crook that \
                 measures wrong is a data error, not a detuned horn"
            ),
            CrookError::JunctionAreaMismatch {
                extracted_step,
                declared_step,
                band,
            } => write!(
                f,
                "FS-QUERY-CROOK-JUNCTION: extracted area step {extracted_step:.4} vs \
                 declared {declared_step:.4} beyond band {band:.4} — an undeclared port \
                 mismatch is refused, never silently voiced"
            ),
        }
    }
}

impl core::error::Error for CrookError {}

impl From<BoreError> for CrookError {
    fn from(e: BoreError) -> Self {
        CrookError::Bore(e)
    }
}

/// Extract one crook's delta-L receipt from its lumen geometry.
///
/// `main_bore_area_m2` is the area the crook must match at its ports
/// (up to the declared step).
///
/// # Errors
/// [`CrookError`] — band and junction refusals by name, plus any
/// underlying bore refusal.
pub fn extract_crook_delta(
    chart: &dyn Chart,
    boundary: &Soup,
    config: &CrookConfig,
    cad: &CrookCadIntent,
    main_bore_area_m2: f64,
    label: &str,
    cx: &Cx<'_>,
) -> Result<ValveChartDelta, CrookError> {
    if !(cad.centerline_length_m > 0.0 && cad.centerline_length_m.is_finite()) {
        return Err(CrookError::Invalid {
            what: "CAD centerline length must be positive finite",
        });
    }
    if !(config.length_band > 0.0 && config.junction_band > 0.0) {
        return Err(CrookError::Invalid {
            what: "bands must be positive",
        });
    }
    if !(main_bore_area_m2 > 0.0 && main_bore_area_m2.is_finite()) {
        return Err(CrookError::Invalid {
            what: "main bore area must be positive finite",
        });
    }
    if !(cad.declared_area_step > 0.0 && cad.declared_area_step.is_finite()) {
        return Err(CrookError::Invalid {
            what: "declared area step must be positive finite",
        });
    }
    let bore = extract_bore(chart, boundary, &config.bore, label, cx)?;
    let extracted_m = bore.total_length_m;
    let deviation = (extracted_m - cad.centerline_length_m).abs() / cad.centerline_length_m;
    if deviation > config.length_band {
        return Err(CrookError::DeltaLengthOutOfBand {
            extracted_m,
            cad_m: cad.centerline_length_m,
            deviation,
            band: config.length_band,
        });
    }
    // Port areas are read at near-end INTERIOR stations: the bore
    // extractor's cut-face end stations carry the disclosed tangent-tilt
    // inflation/clipping (its own CONTRACT row), so the honest port
    // measurement sits a couple of stations inboard.
    let n_st = bore.stations.len();
    if n_st < 7 {
        return Err(CrookError::Invalid {
            what: "too few stations to read interior port areas",
        });
    }
    let port_area_m2 = 0.5 * (bore.stations[2].area_m2 + bore.stations[n_st - 3].area_m2);
    let area_step = port_area_m2 / main_bore_area_m2;
    let step_deviation = (area_step - cad.declared_area_step).abs() / cad.declared_area_step;
    if step_deviation > config.junction_band {
        return Err(CrookError::JunctionAreaMismatch {
            extracted_step: area_step,
            declared_step: cad.declared_area_step,
            band: config.junction_band,
        });
    }
    Ok(ValveChartDelta {
        label: label.to_string(),
        delta_l_m: extracted_m,
        cad_length_m: cad.centerline_length_m,
        length_deviation: deviation,
        port_area_m2,
        area_step,
        bore,
    })
}

impl ValveChartDelta {
    /// JSON log row: extracted vs CAD, bands, junction continuity — the
    /// bead's logging clause.
    #[must_use]
    pub fn debug_line(&self) -> String {
        format!(
            "{{\"suite\":\"fs-query\",\"case\":\"crook-delta\",\"label\":\"{}\",\
             \"delta_l_m\":{:.6e},\"cad_length_m\":{:.6e},\"length_deviation\":{:.4e},\
             \"port_area_m2\":{:.6e},\"area_step\":{:.4},\"boundary_digest\":{},\
             \"stations\":{}}}",
            self.label,
            self.delta_l_m,
            self.cad_length_m,
            self.length_deviation,
            self.port_area_m2,
            self.area_step,
            self.bore.boundary_digest,
            self.bore.stations.len()
        )
    }
}

/// Certify a continuous slide range from its two endpoint extractions.
///
/// # Errors
/// [`CrookError::Invalid`] when the endpoints are inverted or labels
/// disagree with the intent.
pub fn certify_slide_range(
    label: &str,
    shortest: &ValveChartDelta,
    longest: &ValveChartDelta,
) -> Result<SlideRange, CrookError> {
    if longest.delta_l_m <= shortest.delta_l_m {
        return Err(CrookError::Invalid {
            what: "slide endpoints inverted (longest must exceed shortest)",
        });
    }
    Ok(SlideRange {
        label: label.to_string(),
        min_delta_l_m: shortest.delta_l_m,
        max_delta_l_m: longest.delta_l_m,
    })
}
