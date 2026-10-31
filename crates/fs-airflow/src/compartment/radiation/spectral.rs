//! Nongray, nonscattering transfer through the original reciprocal enclosure.
//!
//! Molecular absorption is supplied as constant coefficients on explicit,
//! contiguous wavelength bands. Optional soot uses the original Rayleigh
//! `kappa_lambda` at each wavelength; it is not replaced by one Planck mean.
//! Every surface/layer emits its own Planck spectrum. The original chord
//! expansion, diffuse reflection equations and prepared LU own the transfer.
//!
//! Integration uses positive 8/16-point Gauss--Legendre rules in log wavelength,
//! split at every supplied band boundary. Their difference is an **estimate**,
//! not a certified error bound. Initial panels span at most one logarithmic
//! wavelength unit, so a broad window cannot skip the smooth Planck peak.
//! Both the work limit and absolute per-output tolerances are explicit.
//! Constant band optics commute with wavelength integration: their Planck
//! sources are integrated analytically and one prepared solve is performed per
//! band. Adaptive wavelength quadrature is reserved for genuinely varying soot
//! opacity, where the transfer operator itself changes with wavelength.
//!
//! The Planck shape follows NIST SP 250-1, section 2:
//! <https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nbsspecialpublication250-1.pdf>.
//! Vacuum wavelength and the exact SI `h c / k_B` are used. Its normalization
//! is the original conduction owner's Stefan--Boltzmann constant, retaining
//! exactly that owner's gray energy convention. Tail evaluation uses the
//! exponential series of `1/(exp(x)-1)` and, for small x, the Bernoulli series
//! in <https://dlmf.nist.gov/24.2.E1>.
//!
//! This is a sourced band-absorption approximation, not a molecular line database
//! or a RADCAL replacement. The caller resolves absorption at the actual gas
//! composition, pressure and temperature and justifies the band's resolution.
//! Surfaces use either explicitly gray total emissivities or independent,
//! source-bounded spectral-band emissivities. At every wavelength an opaque
//! diffuse surface's absorptivity equals its emissivity (Kirchhoff's law), and
//! its reflectivity is their complement. Total emissivity is never relabeled
//! as a measured spectral property. Optical geometry, omitted wavelengths,
//! band-data error, scattering and particle/gas nonequilibrium remain separate
//! from the reported quadrature and algebraic errors.
//! Booking soot emission into a gas layer's thermal account assumes that its
//! particles and gas share the supplied layer temperature.

use super::soot::RayleighSoot;
use super::{
    LayerAbsorption, MAX_RADIATION_SURFACES, ParticipatingEnclosure, RadiationError,
    RadiationSolve, RadiationWorkspace, checkpoint, finite, nonnegative, positive,
};
use fs_conduction::radiation::{BlackbodyBand, STEFAN_BOLTZMANN_W_M2_K4, SurfaceEmissivity};
use fs_gas::multicomponent::compartment::{Layer, LayerRates, TwoZone};
use fs_qty::{Length, Power, Temperature};
use std::sync::Arc;

/// Fixed bound on the source-supplied wavelength bands.
pub const MAX_SPECTRAL_BANDS: usize = 256;
/// Hard ceiling on wavelength evaluations in one returned result.
pub const MAX_SPECTRAL_EVALUATIONS: usize = 65_536;
const MAX_DEPTH: usize = 24;
const MAX_CHANNELS: usize = 2 * MAX_RADIATION_SURFACES + 3;
/// Explicit molecular/background absorption over a finite wavelength interval.
/// Adjacent bands must touch exactly; gaps are never assumed transparent.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralBand {
    /// Increasing vacuum wavelengths [m], both strictly positive and finite.
    pub wavelengths: [Length; 2],
    /// Frozen `[lower, upper]` absorption [1/m], with each original source.
    /// Zero is an explicitly transparent molecular/background contribution.
    pub absorption: [LayerAbsorption; 2],
}

/// Source-supplied diffuse, opaque spectral emissivity, constant within one
/// enclosing wavelength band and the declared surface-temperature range.
/// This is distinct from a material card's hemispherical TOTAL emissivity.
#[derive(Clone, Debug, PartialEq)]
pub struct BandSurfaceEmissivity {
    value: f64,
    temperature_range: [Temperature; 2],
    source: String,
}
impl BandSurfaceEmissivity {
    /// Bind the band value to its measured/model source and temperature domain.
    /// Zero explicitly denotes a perfect reflector in this band. Each connected
    /// optical component must still reach an emitting surface or absorbing gas
    /// when evaluated, so its reflected radiation has a unique thermal anchor.
    pub fn new(
        value: f64,
        temperature_range: [Temperature; 2],
        source: impl Into<String>,
    ) -> Result<Self, RadiationError> {
        let source = source.into();
        let [lo, hi] = temperature_range.map(Temperature::value);
        if !nonnegative(value)
            || value > 1.0
            || !positive(lo)
            || !positive(hi)
            || lo > hi
            || source.trim().is_empty()
        {
            return Err(RadiationError::Invalid(
                "sourced band emissivity in [0,1] and a positive temperature domain required",
            ));
        }
        Ok(Self {
            value,
            temperature_range,
            source,
        })
    }
    /// Both emissivity and absorptivity within this opaque diffuse band.
    pub const fn value(&self) -> f64 {
        self.value
    }
    /// Closed admitted surface-temperature interval.
    pub const fn temperature_range(&self) -> [Temperature; 2] {
        self.temperature_range
    }
    /// Original spectral source, including the applicability of a constant band.
    pub fn source(&self) -> &str {
        &self.source
    }
}

/// Original Rayleigh optical model and actual particulate loading.
/// The finite wavelength window is declared by the enclosing band's endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SootPopulation<'a> {
    /// Original constant-index optical model, including its source and domain.
    pub model: &'a RayleighSoot,
    /// Condensed soot volume per mixture volume. Never ideal-gas moles.
    pub volume_fraction: f64,
}

/// Immutable sourced band optics, optionally augmented by spectral soot.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralOptics<'a> {
    bands: Vec<SpectralBand>,
    soot: [Option<SootPopulation<'a>>; 2],
    gray_surface_assumption: String,
    // Immutable source bank is shared by reports and nonlinear wall trials;
    // copying a report must not reallocate every band's per-surface strings.
    surface_bands: Option<Arc<[Vec<BandSurfaceEmissivity>]>>,
    varying: bool,
}
impl<'a> SpectralOptics<'a> {
    /// Admit finite contiguous bands and both original soot populations once.
    /// Surface card values remain unchanged; the supplied text records why a
    /// wavelength-independent boundary emissivity is acceptable for this use.
    pub fn new(
        bands: Vec<SpectralBand>,
        soot: [Option<SootPopulation<'a>>; 2],
        gray_surface_assumption: impl Into<String>,
    ) -> Result<Self, RadiationError> {
        let gray_surface_assumption = gray_surface_assumption.into();
        if gray_surface_assumption.trim().is_empty() {
            return Err(RadiationError::Invalid(
                "bounded spectral bands and explicit gray surface assumption required",
            ));
        }
        let varying = Self::admit_bands(&bands, &soot)?;
        Ok(Self {
            bands,
            soot,
            gray_surface_assumption,
            surface_bands: None,
            varying,
        })
    }

    /// Use genuine spectral surface data instead of extending total emissivity
    /// over wavelength. The outer vector matches `bands` exactly; every inner
    /// vector follows the original enclosure surface axis. Gas and surface
    /// discontinuities must both appear in this common band partition.
    /// No gray emissivity cards are needed by `evaluate_spectral` in this mode.
    pub fn new_with_surface_bands(
        bands: Vec<SpectralBand>,
        soot: [Option<SootPopulation<'a>>; 2],
        surface_bands: Vec<Vec<BandSurfaceEmissivity>>,
    ) -> Result<Self, RadiationError> {
        let varying = Self::admit_bands(&bands, &soot)?;
        let n = surface_bands.first().map_or(0, Vec::len);
        if surface_bands.len() != bands.len()
            || !(2..=MAX_RADIATION_SURFACES).contains(&n)
            || surface_bands.iter().any(|row| row.len() != n)
        {
            return Err(RadiationError::Invalid(
                "spectral emissivity bands must match the common band and surface axes",
            ));
        }
        Ok(Self {
            bands,
            soot,
            gray_surface_assumption: String::new(),
            surface_bands: Some(surface_bands.into()),
            varying,
        })
    }

    fn admit_bands(
        bands: &[SpectralBand],
        soot: &[Option<SootPopulation<'_>>; 2],
    ) -> Result<bool, RadiationError> {
        if bands.is_empty() || bands.len() > MAX_SPECTRAL_BANDS {
            return Err(RadiationError::Invalid("bounded spectral bands required"));
        }
        for (i, band) in bands.iter().enumerate() {
            let [a, b] = band.wavelengths.map(Length::value);
            if !positive(a)
                || !positive(b)
                || a >= b
                || a.ln() >= b.ln()
                || (i > 0 && bands[i - 1].wavelengths[1] != band.wavelengths[0])
            {
                return Err(RadiationError::Invalid(
                    "spectral bands must be finite, increasing, contiguous and resolvable",
                ));
            }
        }
        let mut varying = false;
        for population in soot.iter().flatten() {
            let k = population
                .model
                .spectral_coefficient(population.volume_fraction, bands[0].wavelengths[0])?;
            varying |= k.value() > 0.0;
        }
        Ok(varying)
    }

    /// Unmodified original band coefficients and source identities.
    pub fn bands(&self) -> &[SpectralBand] {
        &self.bands
    }
    /// Original particle models and actual loadings.
    pub const fn soot(&self) -> &[Option<SootPopulation<'a>>; 2] {
        &self.soot
    }
    /// Actual finite spectral window. No opacity is assumed outside this range.
    pub fn wavelengths(&self) -> [Length; 2] {
        [
            self.bands[0].wavelengths[0],
            self.bands[self.bands.len() - 1].wavelengths[1],
        ]
    }
    /// Named assumption extending original total emissivities over wavelength.
    /// Empty when original spectral-band surface inputs replace that assumption.
    pub fn gray_surface_assumption(&self) -> &str {
        &self.gray_surface_assumption
    }

    /// Original spectral surface inputs in band-then-surface order, if selected.
    pub fn surface_emissivity_bands(&self) -> Option<&[Vec<BandSurfaceEmissivity>]> {
        self.surface_bands.as_deref()
    }

    fn coefficients(&self, band: usize, wavelength: Length) -> Result<[f64; 2], RadiationError> {
        let mut k = self.bands[band]
            .absorption
            .each_ref()
            .map(|a| a.coefficient().value());
        for (layer, population) in self.soot.iter().enumerate() {
            if let Some(population) = population {
                k[layer] = finite(
                    k[layer]
                        + population
                            .model
                            .spectral_coefficient(population.volume_fraction, wavelength)?
                            .value(),
                )?;
            }
        }
        Ok(k)
    }

    fn participating(&self, layer: usize) -> bool {
        self.bands
            .iter()
            .any(|b| b.absorption[layer].coefficient().value() > 0.0)
            || self.soot[layer]
                .is_some_and(|p| p.volume_fraction > 0.0 && p.model.refractive_index()[1] > 0.0)
    }
}

/// Numerical controls; source, band and optical-geometry errors are separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralSolve {
    /// Original algebraic gates, also applied to the final integrated balance.
    pub algebraic: RadiationSolve,
    /// Absolute estimated quadrature error per signed gas/surface power [W].
    pub heat_tolerance_w: f64,
    /// Absolute estimated quadrature error per integrated radiosity [W/m²].
    pub radiosity_tolerance_w_m2: f64,
    /// Largest allowed omitted blackbody fraction for any participating node.
    /// This gates spectral coverage, not the unknown omitted net heat exchange.
    pub maximum_omitted_blackbody_fraction: f64,
    /// Actual evaluation ceiling, including rejected quadrature panels.
    pub maximum_evaluations: usize,
}
impl Default for SpectralSolve {
    fn default() -> Self {
        Self {
            algebraic: RadiationSolve::default(),
            heat_tolerance_w: 1e-4,
            radiosity_tolerance_w_m2: 1e-4,
            maximum_omitted_blackbody_fraction: 1e-4,
            maximum_evaluations: 8192,
        }
    }
}

/// Blackbody emission outside the declared spectral window at one actual T.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralCoverage {
    /// Actual boundary/layer temperature used in the emission law.
    pub temperature: Temperature,
    /// Whether this node participates in any admitted optical exchange.
    pub participating: bool,
    /// Blackbody fraction at wavelengths shorter than the lower endpoint.
    pub shortwave_fraction: f64,
    /// Blackbody fraction at wavelengths longer than the upper endpoint.
    pub longwave_fraction: f64,
    /// Omitted blackbody emission density [W/m²], not a net-heat error bound.
    pub omitted_blackbody_w_m2: f64,
}

/// One integrated, frozen-state nongray enclosure evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralRadiationReport<'a> {
    /// Watts INTO the original boundary patches, integrated over this window.
    pub surface_heat_power: Vec<Power>,
    /// Watts INTO the original `[lower, upper]` thermal gas/particle accounts.
    pub gas_heat_power: [Power; 2],
    /// Outgoing surface radiosities integrated over the wavelength window [W/m²].
    pub surface_radiosity_w_m2: Vec<f64>,
    /// Estimated absolute quadrature errors for each surface heat power [W].
    pub surface_heat_error_estimate_w: Vec<f64>,
    /// Estimated absolute quadrature errors for lower/upper gas powers [W].
    pub gas_heat_error_estimate_w: [f64; 2],
    /// Estimated absolute quadrature errors for each radiosity [W/m²].
    pub radiosity_error_estimate_w_m2: Vec<f64>,
    /// Applied original gray cards; empty when `optics` retains the actual
    /// band-resolved surface inputs instead. No unused gray source is reported.
    pub surface_emissivities: Vec<SurfaceEmissivity>,
    /// Unmodified optical model and assumptions used in this calculation.
    pub optics: SpectralOptics<'a>,
    /// Coverage diagnostics in original surface order, then lower and upper gas.
    pub coverage: Vec<SpectralCoverage>,
    /// Integrated absolute surface-equation defect estimate [W/m²].
    pub surface_residual_integral_w_m2: f64,
    /// Actual signed sum of all integrated powers [W], never redistributed.
    pub energy_residual_w: f64,
    /// Actual point-wavelength evaluations, including rejected panels. Fixed
    /// band optics integrate Planck sources analytically and contribute zero.
    pub wavelength_evaluations: usize,
    /// Original LU factorizations performed; fixed band optics use one per band.
    pub factorizations: usize,
    /// Accepted adaptive 16-point panels. Fixed bands contribute zero panels.
    pub accepted_panels: usize,
}
impl SpectralRadiationReport<'_> {
    /// Add the identical gas powers atomically to the existing source account.
    /// Apply corresponding surface powers to their wall owners exactly once.
    pub fn add_gas_heat(&self, rates: &mut [LayerRates; 2]) -> Result<(), RadiationError> {
        let mut candidate = *rates;
        for (rate, &heat) in candidate.iter_mut().zip(&self.gas_heat_power) {
            rate.add_heat(heat)?;
        }
        *rates = candidate;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Panel {
    lower: f64,
    upper: f64,
    depth: usize,
}

/// Reusable original radiation buffers and bounded spectral accumulators.
/// No storage is allocated per wavelength or per adaptive refinement.
pub struct SpectralWorkspace {
    radiation: RadiationWorkspace,
    totals: [f64; MAX_CHANNELS],
    corrections: [f64; MAX_CHANNELS],
    errors: [f64; MAX_CHANNELS],
    pending: Vec<Panel>,
}

impl ParticipatingEnclosure {
    /// Prepare storage once; it can be reused for any admitted spectral optics
    /// on this original boundary count. Failed calls only overwrite scratch.
    pub fn spectral_workspace(&self) -> Result<SpectralWorkspace, RadiationError> {
        Ok(SpectralWorkspace {
            radiation: self.workspace()?,
            totals: [0.0; MAX_CHANNELS],
            corrections: [0.0; MAX_CHANNELS],
            errors: [0.0; MAX_CHANNELS],
            pending: Vec::with_capacity(MAX_DEPTH + 1),
        })
    }

    /// Use original room temperatures, refusing paths through absent layers.
    /// Paths must still describe the actual accepted interface geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_room_spectral<'a>(
        &self,
        room: &TwoZone,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        optics: &SpectralOptics<'a>,
        controls: SpectralSolve,
        work: &mut SpectralWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpectralRadiationReport<'a>, RadiationError> {
        for path in self.paths() {
            checkpoint(&mut cancelled)?;
            for segment in &path.segments {
                if segment.length.value() > 0.0 && room.layer(segment.layer).volume().value() == 0.0
                {
                    return Err(RadiationError::Invalid(
                        "optical path crosses an absent gas layer",
                    ));
                }
            }
        }
        self.evaluate_spectral(
            surface_temperatures,
            emissivities,
            [
                room.layer(Layer::Lower).temperature(),
                room.layer(Layer::Upper).temperature(),
            ],
            optics,
            controls,
            work,
            cancelled,
        )
    }

    /// Resolve each node's emission spectrum and wavelength-dependent transfer.
    /// Surface data is either the original gray cards or `optics`' independent
    /// band sources. Pass an empty gray slice with band sources; supplying both
    /// refuses ambiguous boundary ownership. This call publishes no state or
    /// time-order claim.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_spectral<'a>(
        &self,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        gas_temperatures: [Temperature; 2],
        optics: &SpectralOptics<'a>,
        controls: SpectralSolve,
        work: &mut SpectralWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpectralRadiationReport<'a>, RadiationError> {
        let gray = optics.surface_bands.is_none();
        if !gray && !emissivities.is_empty() {
            return Err(RadiationError::Invalid(
                "spectral surface bands and gray emissivity cards cannot both own a boundary",
            ));
        }
        self.admit_state(
            surface_temperatures,
            gray.then_some(emissivities),
            gas_temperatures,
            controls.algebraic,
            &work.radiation,
            &mut cancelled,
        )?;
        if !positive(controls.heat_tolerance_w)
            || !positive(controls.radiosity_tolerance_w_m2)
            || !nonnegative(controls.maximum_omitted_blackbody_fraction)
            || controls.maximum_omitted_blackbody_fraction > 1.0
            || !(24..=MAX_SPECTRAL_EVALUATIONS).contains(&controls.maximum_evaluations)
        {
            return Err(RadiationError::Invalid("invalid spectral solve controls"));
        }
        let n = self.geometry().len();
        let m = n + 2;
        if let Some(bands) = &optics.surface_bands {
            for band in bands.iter() {
                checkpoint(&mut cancelled)?;
                if band.len() != n {
                    return Err(RadiationError::Invalid(
                        "spectral surface axis differs from enclosure",
                    ));
                }
                for (source, temperature) in band.iter().zip(surface_temperatures) {
                    let [lo, hi] = source.temperature_range.map(Temperature::value);
                    if !(lo..=hi).contains(&temperature.value()) {
                        return Err(RadiationError::Invalid(
                            "surface temperature is outside spectral emissivity validity",
                        ));
                    }
                }
            }
        }
        let channels = m + n + 1;
        let wavelengths = optics.wavelengths();
        let log_range = finite(wavelengths[1].value().ln() - wavelengths[0].value().ln())?;
        let mut temperatures = [0.0; MAX_RADIATION_SURFACES + 2];
        let mut coverage = Vec::with_capacity(m);
        let mut reference = f64::INFINITY;
        let mut maximum_omitted = 0.0_f64;
        for (i, t) in surface_temperatures
            .iter()
            .chain(gas_temperatures.iter())
            .enumerate()
        {
            checkpoint(&mut cancelled)?;
            temperatures[i] = t.value();
            // A perfect reflector has no thermal spectrum of its own. Its
            // temperature cannot choose the subtracted background or reject
            // the wavelength window for emission that is identically zero.
            let participating = if i < n {
                optics
                    .surface_bands
                    .as_ref()
                    .is_none_or(|bands| bands.iter().any(|band| band[i].value > 0.0))
            } else {
                optics.participating(i - n)
                    && self.paths().iter().any(|p| {
                        p.segments
                            .iter()
                            .any(|s| s.layer.index() == i - n && s.length.value() > 0.0)
                    })
            };
            let record = coverage_at(wavelengths, *t, participating)?;
            if participating {
                reference = reference.min(t.value());
                maximum_omitted =
                    maximum_omitted.max(record.shortwave_fraction + record.longwave_fraction);
            }
            coverage.push(record);
        }
        if !reference.is_finite() {
            return Err(RadiationError::Invalid(
                "closed lossless enclosure has no thermal radiation anchor",
            ));
        }
        if maximum_omitted > controls.maximum_omitted_blackbody_fraction {
            return Err(RadiationError::SpectralTruncation {
                maximum_fraction: maximum_omitted,
                allowed_fraction: controls.maximum_omitted_blackbody_fraction,
            });
        }
        work.totals.fill(0.0);
        work.corrections.fill(0.0);
        work.errors.fill(0.0);
        work.pending.clear();
        let mut evaluations = 0;
        let mut factorizations = 0;
        let mut accepted_panels = 0;
        for (band, source) in optics.bands.iter().enumerate() {
            checkpoint(&mut cancelled)?;
            let mut values = [0.0; MAX_RADIATION_SURFACES];
            if let Some(bands) = &optics.surface_bands {
                for (value, source) in values.iter_mut().zip(&bands[band]) {
                    *value = source.value;
                }
            } else {
                for (value, source) in values.iter_mut().zip(emissivities) {
                    *value = source.value();
                }
            }
            let begin = source.wavelengths[0].value().ln();
            let end = source.wavelengths[1].value().ln();
            if !optics.varying {
                self.prepare_optics(
                    &values[..n],
                    optics.coefficients(band, source.wavelengths[0])?,
                    &mut work.radiation,
                    &mut cancelled,
                )?;
                factorizations += 1;

                // With wavelength-independent opacity and emissivity, the
                // complete diffuse/reflection operator is linear and constant
                // across this band. Integrate each Planck source first, then
                // solve the identical radiosity system once. This is both more
                // accurate and substantially cheaper than sampling wavelength.
                let blackbody = BlackbodyBand::new(
                    source.wavelengths[0].value(),
                    source.wavelengths[1].value(),
                )
                .map_err(|_| RadiationError::Invalid("shared blackbody spectral window refused"))?;
                let reference_emission = blackbody
                    .evaluate(reference)
                    .map_err(|_| RadiationError::Invalid("shared blackbody band evaluation refused"))?
                    .power_w_m2;
                for (i, (value, &temperature)) in work
                    .radiation
                    .emission
                    .iter_mut()
                    .zip(&temperatures[..m])
                    .enumerate()
                {
                    *value = if i < n && values[i] == 0.0 {
                        0.0
                    } else {
                        blackbody
                            .difference(temperature, reference)
                            .map_err(|_| RadiationError::Invalid("shared blackbody band evaluation refused"))?
                    };
                }
                let (residual, _) = self.solve_emission(
                    &values[..n],
                    controls.algebraic,
                    &mut work.radiation,
                    &mut cancelled,
                )?;
                for i in 0..m {
                    let y = finite(work.radiation.powers[i] - work.corrections[i])?;
                    let next = finite(work.totals[i] + y)?;
                    work.corrections[i] = (next - work.totals[i]) - y;
                    work.totals[i] = next;
                }
                for i in 0..n {
                    let radiosity = finite(reference_emission + work.radiation.radiosity[i])?;
                    if radiosity < 0.0 {
                        return Err(RadiationError::Invalid("negative spectral radiosity"));
                    }
                    let channel = m + i;
                    let y = finite(radiosity - work.corrections[channel])?;
                    let next = finite(work.totals[channel] + y)?;
                    work.corrections[channel] = (next - work.totals[channel]) - y;
                    work.totals[channel] = next;
                }
                work.totals[m + n] = finite(work.totals[m + n] + residual)?;
                continue;
            }
            let mut sample = SpectralSample {
                enclosure: self,
                emissivities: &values[..n],
                optics,
                temperatures: &temperatures[..m],
                reference,
                controls: RadiationSolve {
                    surface_residual_w_m2: controls.algebraic.surface_residual_w_m2 / log_range,
                    energy_residual_w: controls.algebraic.energy_residual_w / log_range,
                },
            };
            let mut left = begin;
            while left < end {
                let right = (left + 1.0).min(end);
                if right <= left {
                    return Err(RadiationError::Invalid("unresolvable spectral panel"));
                }
                work.pending.push(Panel {
                    lower: left,
                    upper: right,
                    depth: 0,
                });
                while let Some(panel) = work.pending.pop() {
                    checkpoint(&mut cancelled)?;
                    if evaluations + 24 > controls.maximum_evaluations {
                        return Err(RadiationError::SpectralBudget {
                            evaluations,
                            maximum_evaluations: controls.maximum_evaluations,
                        });
                    }
                    let mut coarse = [0.0; MAX_CHANNELS];
                    let mut fine = [0.0; MAX_CHANNELS];
                    sample.rule(
                        band,
                        panel,
                        &GL8,
                        &mut coarse[..channels],
                        &mut work.radiation,
                        &mut evaluations,
                        &mut factorizations,
                        &mut cancelled,
                    )?;
                    sample.rule(
                        band,
                        panel,
                        &GL16,
                        &mut fine[..channels],
                        &mut work.radiation,
                        &mut evaluations,
                        &mut factorizations,
                        &mut cancelled,
                    )?;
                    let share = (panel.upper - panel.lower) / log_range;
                    let resolved = (0..m + n).all(|i| {
                        let tolerance = if i < m {
                            controls.heat_tolerance_w
                        } else {
                            controls.radiosity_tolerance_w_m2
                        };
                        (fine[i] - coarse[i]).abs() <= tolerance * share
                    });
                    if !resolved {
                        let middle = panel.lower + 0.5 * (panel.upper - panel.lower);
                        if panel.depth == MAX_DEPTH
                            || middle <= panel.lower
                            || middle >= panel.upper
                        {
                            return Err(RadiationError::Invalid(
                                "spectral refinement reached numerical resolution",
                            ));
                        }
                        work.pending.push(Panel {
                            lower: middle,
                            upper: panel.upper,
                            depth: panel.depth + 1,
                        });
                        work.pending.push(Panel {
                            lower: panel.lower,
                            upper: middle,
                            depth: panel.depth + 1,
                        });
                        continue;
                    }
                    accepted_panels += 1;
                    for i in 0..channels {
                        let y = finite(fine[i] - work.corrections[i])?;
                        let next = finite(work.totals[i] + y)?;
                        work.corrections[i] = (next - work.totals[i]) - y;
                        work.totals[i] = next;
                        work.errors[i] = finite(work.errors[i] + (fine[i] - coarse[i]).abs())?;
                    }
                }
                left = right;
            }
        }
        let balance = finite(work.totals[..m].iter().sum())?;
        let surface_residual = work.totals[m + n];
        if balance.abs() > controls.algebraic.energy_residual_w
            || surface_residual > controls.algebraic.surface_residual_w_m2
        {
            return Err(RadiationError::Residual {
                surface_w_m2: surface_residual,
                energy_w: balance,
            });
        }
        checkpoint(&mut cancelled)?;
        Ok(SpectralRadiationReport {
            surface_heat_power: work.totals[..n].iter().copied().map(Power::new).collect(),
            gas_heat_power: [Power::new(work.totals[n]), Power::new(work.totals[n + 1])],
            surface_radiosity_w_m2: work.totals[m..m + n].to_vec(),
            surface_heat_error_estimate_w: work.errors[..n].to_vec(),
            gas_heat_error_estimate_w: [work.errors[n], work.errors[n + 1]],
            radiosity_error_estimate_w_m2: work.errors[m..m + n].to_vec(),
            surface_emissivities: emissivities.to_vec(),
            optics: optics.clone(),
            coverage,
            surface_residual_integral_w_m2: surface_residual,
            energy_residual_w: balance,
            wavelength_evaluations: evaluations,
            factorizations,
            accepted_panels,
        })
    }
}

struct SpectralSample<'a, 'b> {
    enclosure: &'a ParticipatingEnclosure,
    emissivities: &'a [f64],
    optics: &'a SpectralOptics<'b>,
    temperatures: &'a [f64],
    reference: f64,
    controls: RadiationSolve,
}
impl SpectralSample<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    fn rule(
        &mut self,
        band: usize,
        panel: Panel,
        rule: &[(f64, f64)],
        out: &mut [f64],
        work: &mut RadiationWorkspace,
        evaluations: &mut usize,
        factorizations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RadiationError> {
        let half = 0.5 * (panel.upper - panel.lower);
        let midpoint = panel.lower + half;
        let n = self.emissivities.len();
        let m = n + 2;
        for &(node, weight) in rule {
            for sign in [-1.0, 1.0] {
                checkpoint(cancelled)?;
                let wavelength = finite((midpoint + sign * half * node).exp())?;
                if self.optics.varying {
                    self.enclosure.prepare_optics(
                        self.emissivities,
                        self.optics.coefficients(band, Length::new(wavelength))?,
                        work,
                        cancelled,
                    )?;
                    *factorizations += 1;
                }
                let reference_emission = planck_log_emission(wavelength, self.reference)?;
                for (i, (value, &temperature)) in
                    work.emission.iter_mut().zip(self.temperatures).enumerate()
                {
                    *value = if i < n && self.emissivities[i] == 0.0 {
                        0.0
                    } else {
                        planck_log_difference(wavelength, temperature, self.reference)?
                    };
                }
                let (residual, _) = self.enclosure.solve_emission(
                    self.emissivities,
                    self.controls,
                    work,
                    cancelled,
                )?;
                *evaluations += 1;
                let measure = half * weight;
                for i in 0..m {
                    out[i] = finite(measure.mul_add(work.powers[i], out[i]))?;
                }
                for i in 0..n {
                    let radiosity = finite(reference_emission + work.radiosity[i])?;
                    if radiosity < 0.0 {
                        return Err(RadiationError::Invalid("negative spectral radiosity"));
                    }
                    out[m + i] = finite(measure.mul_add(radiosity, out[m + i]))?;
                }
                out[m + n] = finite(measure.mul_add(residual, out[m + n]))?;
            }
        }
        Ok(())
    }
}

// Positive nodes/weights for unit-interval symmetric Gauss--Legendre rules.
// Full rules use both signs. Their measure is d(log lambda), not d lambda.
const GL8: [(f64, f64); 4] = [
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];
const GL16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_5),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_7),
    (0.755_404_408_355_003, 0.124_628_971_255_533_9),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_09),
];

fn planck_log_emission(wavelength: f64, temperature: f64) -> Result<f64, RadiationError> {
    BlackbodyBand::spectral_power_per_log_wavelength(wavelength, temperature)
        .map_err(|_| RadiationError::Invalid("shared blackbody spectral evaluation refused"))
}

fn planck_log_difference(
    wavelength: f64,
    temperature: f64,
    reference: f64,
) -> Result<f64, RadiationError> {
    BlackbodyBand::spectral_difference_per_log_wavelength(wavelength, temperature, reference)
        .map_err(|_| RadiationError::Invalid("shared blackbody spectral evaluation refused"))
}

fn coverage_at(
    wavelengths: [Length; 2],
    temperature: Temperature,
    participating: bool,
) -> Result<SpectralCoverage, RadiationError> {
    let window = BlackbodyBand::new(wavelengths[0].value(), wavelengths[1].value())
        .map_err(|_| RadiationError::Invalid("shared blackbody spectral window refused"))?;
    let emission = window
        .evaluate(temperature.value())
        .map_err(|_| RadiationError::Invalid("shared blackbody spectral coverage refused"))?;
    let omitted = finite(
        (emission.short_wave_fraction + emission.long_wave_fraction)
            * STEFAN_BOLTZMANN_W_M2_K4
            * temperature.value().powi(4),
    )?;
    Ok(SpectralCoverage {
        temperature,
        participating,
        shortwave_fraction: emission.short_wave_fraction,
        longwave_fraction: emission.long_wave_fraction,
        omitted_blackbody_w_m2: omitted,
    })
}

#[cfg(test)]
mod tests;
