//! Declared scalar-flux uncertainty on the original direct transport model.
//!
//! For source i, let E_i be its exact retained-model direct lux and let g_i
//! multiply its retained emitted lumens, AFTER dimming and filtering. A supplied
//! [l_i,u_i] asserts l_i <= g_i <= u_i at the declared operating state. Positivity
//! gives sum(l_i E_i) <= E(g) <= sum(u_i E_i). No distribution, independence,
//! confidence level or root-sum-square combination is assumed. Sources sharing
//! an instrument or calibration error may be arbitrarily correlated.
//!
//! The bounds are conditional on the caller's declarations. An evidence hash
//! identifies the calibration/assumption record; it does not authenticate it.
//! Statistical standard deviations or coverage intervals are NOT deterministic
//! bounds. Changes to angular shape, spectrum, placement, medium or material are
//! outside this scalar model. In particular a zero-flux source stays off: stray
//! light or a dimmer's additive black-level offset is NOT bounded here.
//!
//! The nominal numerical certificate is unchanged. Numerical-plus-flux ranges,
//! scalar-flux perturbation bounds and unknown total physical error stay separate.

use super::{
    CalculationGrid, ContentHash, Cx, DomainHasher, EnclosedIlluminance,
    GridEnclosureBudget, IlluminanceProbe, Interval, PhotometricScene,
    PhotometryError, Prepared, Receivers, evaluate_with_contributions, finite,
    hash_float, nonnegative, prepare,
};
use std::io::{self, Write};

/// Infer a common gain from bounded lux measurements using the original transport.
pub mod calibration;

/// Version of the scalar-gain declaration, accumulation and output contract.
pub const FLUX_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

/// Invalid caller-supplied flux-uncertainty declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluxBoundError {
    /// Endpoints must be finite, ordered and nonnegative.
    InvalidRange,
    /// A nonzero identity must name the calibration or assumption record.
    MissingEvidence,
}

impl core::fmt::Display for FluxBoundError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "flux-bound declaration refused: {self:?}")
    }
}
impl std::error::Error for FluxBoundError {}

/// Immutable multiplicative bound on one retained emitted-flux coefficient.
///
/// The interval need not contain one: a documented calibration bias may place
/// the entire admissible flux above or below the original nominal value. It
/// never modifies that nominal source or any accepted renderer sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxScaleBound {
    scale: Interval,
    evidence: ContentHash,
}

impl FluxScaleBound {
    /// Admit an explicitly identified deterministic range, not a standard error.
    /// A synthetic/assumed range must be identified as such by its evidence record.
    pub fn try_new(lower: f64, upper: f64, evidence: ContentHash) -> Result<Self, FluxBoundError> {
        if !lower.is_finite() || !upper.is_finite() || lower < 0.0 || lower > upper {
            return Err(FluxBoundError::InvalidRange);
        }
        if evidence == ContentHash([0; 32]) {
            return Err(FluxBoundError::MissingEvidence);
        }
        Ok(Self { scale: Interval::new(lower.max(0.0), upper.max(0.0)), evidence })
    }

    /// Admitted dimensionless multiplier range on emitted lumens.
    #[must_use]
    pub const fn scale(self) -> Interval { self.scale }

    /// Caller-supplied calibration/assumption record identity, not a signature.
    #[must_use]
    pub const fn evidence_identity(self) -> ContentHash { self.evidence }
}

/// Original numerical result plus a separate conditional scalar-flux account.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxBoundedIlluminance {
    /// Unchanged original nominal value, numerical certificate and total claim.
    pub nominal: EnclosedIlluminance,
    /// Enclosure of direct lux with ALL supplied gains applied to the retained
    /// model. Includes transport arithmetic/visibility uncertainty, not physical
    /// model discrepancy. It may exclude the uncorrected nominal value.
    pub flux_conditioned_lux: Option<Interval>,
    /// Upper bound on |E(g)-E(1)| for exact retained-model source contributions.
    /// This is only the scalar-flux perturbation, not numerical or total error.
    pub flux_perturbation_bound_lux: Option<f64>,
    /// Potentially contributing sources without a supplied gain bound. An
    /// exactly off or provably blocked source is irrelevant to this scalar model.
    /// `None` means the underlying transport enclosure was unavailable.
    pub unbounded_flux_sources: Option<usize>,
    /// Ordered declarations, including missing entries and evidence identities.
    /// This alone does not identify a scene or certify an input measurement.
    pub flux_declaration_identity: ContentHash,
    /// Binds declarations to this probe's original numerical-model identity.
    /// Absent when the numerical transport model could not be enclosed.
    pub combined_identity: Option<ContentHash>,
}

impl FluxBoundedIlluminance {
    /// Bound on |nominal binary64 lux - E(g)|, conditional on the supplied
    /// scalar gains and retained transport model. This is NOT total physical
    /// error. Unavailable/nonrepresentable bounds return `None`, never zero.
    #[must_use]
    pub fn numerical_and_flux_error_bound_lux(self) -> Option<f64> {
        let range = self.flux_conditioned_lux?;
        let value = self.nominal.estimate.lux;
        if range.lo() == value && range.hi() == value { return Some(0.0); }
        let upper = (Interval::point(value) - range).abs().hi();
        upper.is_finite().then_some(upper)
    }

    /// Self-contained conditional-flux CSV; all physical-model totals are unknown.
    pub fn write_csv_header(writer: &mut impl Write) -> io::Result<()> {
        writer.write_all(concat!(
            "index,nominal_lux,numerical_lower_lux,numerical_upper_lux,",
            "numerical_error_bound_lux,flux_conditioned_lower_lux,flux_conditioned_upper_lux,",
            "flux_perturbation_bound_lux,numerical_and_flux_error_bound_lux,",
            "unbounded_flux_sources,flux_status,flux_declaration_identity,combined_identity,",
            "source_construction_error_bound_lux,physical_model_error_bound_lux,",
            "total_error_bound_lux,flux_enclosure_version\n"
        ).as_bytes())
    }

    /// Stream round-trip values without substituting zero for missing evidence.
    /// Writer failures propagate; atomic file publication belongs to the caller.
    pub fn write_csv_row(self, index: usize, writer: &mut impl Write) -> io::Result<()> {
        write!(writer, "{index},{}", self.nominal.estimate.lux)?;
        for value in [
            self.nominal.numerical_error_bound_lux.map(|_| self.nominal.estimate.numerical.lo),
            self.nominal.numerical_error_bound_lux.map(|_| self.nominal.estimate.numerical.hi),
            self.nominal.numerical_error_bound_lux,
            self.flux_conditioned_lux.map(|r| r.lo()),
            self.flux_conditioned_lux.map(|r| r.hi()),
            self.flux_perturbation_bound_lux,
            self.numerical_and_flux_error_bound_lux(),
        ] {
            writer.write_all(b",")?;
            match value { Some(v) => write!(writer, "{v}")?, None => writer.write_all(b"unknown")? }
        }
        writer.write_all(b",")?;
        match self.unbounded_flux_sources {
            Some(count) => write!(writer, "{count}")?,
            None => writer.write_all(b"unknown")?,
        }
        let status = if self.nominal.unavailable.is_some() { "no-numerical-enclosure" }
            else if self.unbounded_flux_sources.is_some_and(|n| n > 0) { "missing-flux-bounds" }
            else if self.flux_conditioned_lux.is_none() { "arithmetic-range" }
            else { "conditional-model-bound" };
        write!(writer, ",{status},{},", self.flux_declaration_identity)?;
        match self.combined_identity {
            Some(identity) => write!(writer, "{identity}")?,
            None => writer.write_all(b"unknown")?,
        }
        writeln!(writer, ",unknown,unknown,unknown,{FLUX_ENCLOSURE_SEMANTICS_VERSION}")
    }
}

impl PhotometricScene<'_> {
    /// Enclose one receiver under caller-declared scalar fixture gains. Missing
    /// entries stay unknown. Use [`FluxEnclosureJob`] to amortize preparation.
    pub fn illuminance_with_flux_bounds(
        &self, cx: &Cx<'_>, probe: &IlluminanceProbe, bounds: &[Option<FluxScaleBound>],
    ) -> Result<FluxBoundedIlluminance, PhotometryError> {
        cx.checkpoint()?;
        if bounds.len() != self.sources.len() { return Err(PhotometryError::InvalidOptions); }
        let identity = declaration_identity(bounds, cx)?;
        evaluate_flux(self, cx, probe, &prepare(self, cx)?, bounds, identity)
    }
}

/// Progressive probes/grids with immutable, source-order flux declarations.
///
/// `bounds[i]` applies to `scene.sources[i]`; `None` means unknown, not exact.
/// The scene, source declarations and probes stay borrowed. Each receiver uses
/// ONE original numerical-enclosure traversal; scalar accumulation adds O(1)
/// work per contributing source, with no source-sized temporary per receiver.
/// The native nominal evaluation remains separate exactly as in the original API.
/// Cancellation/refusal never publishes partially accumulated flux accounts.
pub struct FluxEnclosureJob<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    receivers: Receivers<'a>,
    bounds: &'a [Option<FluxScaleBound>],
    identity: ContentHash,
    prepared: Prepared<'scene>,
    count: usize,
    values: Vec<FluxBoundedIlluminance>,
}

impl<'a, 'scene> FluxEnclosureJob<'a, 'scene> {
    /// Admit original probes and exactly one optional declaration per source.
    pub fn new(
        scene: &'a PhotometricScene<'scene>, probes: &'a [IlluminanceProbe],
        bounds: &'a [Option<FluxScaleBound>], cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Probes(probes), probes.len(), bounds, cx)
    }

    /// Admit a native point-sampled calculation grid without copying probes.
    pub fn for_grid(
        scene: &'a PhotometricScene<'scene>, grid: CalculationGrid,
        bounds: &'a [Option<FluxScaleBound>], cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::begin(scene, Receivers::Grid(grid), grid.len(), bounds, cx)
    }

    fn begin(
        scene: &'a PhotometricScene<'scene>, receivers: Receivers<'a>, count: usize,
        bounds: &'a [Option<FluxScaleBound>], cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if bounds.len() != scene.sources.len() { return Err(PhotometryError::InvalidOptions); }
        if count > scene.options.max_outputs { return Err(PhotometryError::OutputBudgetExceeded); }
        let identity = declaration_identity(bounds, cx)?;
        let prepared = prepare(scene, cx)?;
        let mut values = Vec::new();
        values.try_reserve_exact(count).map_err(|_| PhotometryError::AllocationRefused)?;
        cx.checkpoint()?;
        Ok(Self { scene, receivers, bounds, identity, prepared, count, values })
    }

    /// Change optional medium work only before publication. The native grid
    /// budget remains binding; exhausted bounds never erase successful nominal lux.
    pub fn set_grid_enclosure_budget(&mut self, budget: GridEnclosureBudget) -> Result<(), PhotometryError> {
        if !self.values.is_empty() || budget.max_nodes_per_segment == 0 {
            return Err(PhotometryError::InvalidOptions);
        }
        if let Prepared::Available { medium, .. } = &mut self.prepared { medium.set_grid_budget(budget); }
        Ok(())
    }

    /// Publish at most this many complete receivers, retaining the accepted prefix.
    pub fn advance(&mut self, cx: &Cx<'_>, max_points: usize) -> Result<usize, PhotometryError> {
        cx.checkpoint()?;
        let start = self.values.len();
        let end = start.saturating_add(max_points).min(self.count);
        for index in start..end {
            let probe = match &self.receivers {
                Receivers::Probes(probes) => probes[index],
                Receivers::Grid(grid) => grid.probe(index)?,
            };
            let report = evaluate_flux(
                self.scene, cx, &probe, &self.prepared, self.bounds, self.identity,
            )?;
            cx.checkpoint()?;
            self.values.push(report);
        }
        Ok(self.values.len() - start)
    }

    /// Completed receiver prefix; uncomputed values are absent, not zero.
    #[must_use]
    pub fn results(&self) -> &[FluxBoundedIlluminance] { &self.values }

    /// Every admitted receiver has been published.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.values.len() == self.count }
}

fn evaluate_flux(
    scene: &PhotometricScene<'_>, cx: &Cx<'_>, probe: &IlluminanceProbe,
    prepared: &Prepared<'_>, bounds: &[Option<FluxScaleBound>], identity: ContentHash,
) -> Result<FluxBoundedIlluminance, PhotometryError> {
    let mut account = FluxAccount::new();
    let nominal = evaluate_with_contributions(scene, cx, probe, prepared,
        &mut |source, value| account.add(value, bounds[source]))?;
    let report = account.finish(nominal, identity);
    cx.checkpoint()?;
    Ok(report)
}

fn declaration_identity(bounds: &[Option<FluxScaleBound>], cx: &Cx<'_>) -> Result<ContentHash, PhotometryError> {
    let mut hash = DomainHasher::new("org.frankensim.photometry.flux-declarations.v1");
    hash.update(&FLUX_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
    hash.update(&(bounds.len() as u64).to_le_bytes());
    for bound in bounds {
        cx.checkpoint()?;
        match bound {
            None => { hash.update(&[0]); }
            Some(bound) => {
                hash.update(&[1]);
                hash_float(&mut hash, bound.scale.lo());
                hash_float(&mut hash, bound.scale.hi());
                hash.update(&bound.evidence.0);
            }
        }
    }
    cx.checkpoint()?;
    Ok(hash.finalize())
}

// Private receiver-local accumulator. A failed transport evaluation discards
// this entire account, including callbacks from earlier successful sources.
struct FluxAccount {
    lux: Interval,
    perturbation: f64,
    unknown: usize,
}

impl FluxAccount {
    fn new() -> Self { Self { lux: Interval::point(0.0), perturbation: 0.0, unknown: 0 } }

    fn add(&mut self, value: Interval, bound: Option<FluxScaleBound>) {
        if value.hi() == 0.0 { return; }
        let Some(bound) = bound else { self.unknown += 1; return; };
        let gain = bound.scale;
        let scaled = if gain.hi() == 0.0 { Interval::point(0.0) }
            else if gain == Interval::point(1.0) { value }
            else { nonnegative(value * gain) };
        if scaled.hi() > 0.0 { self.lux = nonnegative(self.lux + scaled); }
        if gain != Interval::point(1.0) {
            let deviation = (gain - Interval::point(1.0)).abs().hi();
            let error = (Interval::point(value.hi()) * Interval::point(deviation)).hi();
            self.perturbation = (Interval::point(self.perturbation) + Interval::point(error)).hi();
        }
    }

    fn finish(self, nominal: EnclosedIlluminance, identity: ContentHash) -> FluxBoundedIlluminance {
        let numerical_available = nominal.unavailable.is_none();
        let supported = numerical_available && self.unknown == 0;
        let combined_identity = nominal.numerical_model_identity.map(|model| {
            let mut hash = DomainHasher::new("org.frankensim.photometry.flux-probe.v1");
            hash.update(&model.0);
            hash.update(&identity.0);
            hash.finalize()
        });
        FluxBoundedIlluminance {
            nominal,
            flux_conditioned_lux: (supported && finite(self.lux)).then_some(self.lux),
            flux_perturbation_bound_lux: (supported && self.perturbation.is_finite()).then_some(self.perturbation),
            unbounded_flux_sources: numerical_available.then_some(self.unknown),
            flux_declaration_identity: identity,
            combined_identity,
        }
    }
}

#[cfg(test)]
mod tests;
