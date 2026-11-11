//! Prompt fluorescence and local heat from the original absorbed beam power.
//!
//! Quantum yield counts emitted photons per absorbed photon, not emitted joules
//! per absorbed joule (IUPAC Gold Book, doi:10.1351/goldbook.Q04991). A photon
//! emitted at lambda_e carries lambda_p/lambda_e of the pump photon's energy.
//! Nonradiative loss and the Stokes energy difference heat the original cell.
//!
//! This restricted consumer requires the SAME non-scattering, index-one medium
//! to have exactly zero absorption at every emitted wavelength, throughout all
//! occupied cells. Its exterior is open vacuum: no surfaces, trapping, bath or
//! incoming fluorescence is inferred. Thus every emitted photon escapes; a
//! self-absorbing source refuses instead of silently losing or heating photons.
//! Emission is prompt on the resolved thermal clock, with no excited-state
//! storage, lifetime, saturation, bleaching, anti-Stokes process or cascade.
//! The uniform authored yield/spectrum and optical snapshot remain fixed over
//! their declared temperature interval. These assumptions need source support;
//! a supplied identity does not confer measurement authority.

use super::{
    CheckedSum, OpticalEnergyBudget, OpticalHeating, OpticalHeatingError, OpticalPowerBudget,
    finite, invalid, nonnegative,
};
use crate::volumes::transport::deposition::BeamPowerDeposition;
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError};
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;

/// A supplied discrete photon spectrum and constant photon quantum yield.
/// Emission weights describe photon counts, NOT watts or spectral density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluorescenceSource<const M: usize> {
    source_identity: ContentHash,
    excitation_bounds_nm: (f64, f64),
    temperature_bounds_k: (f64, f64),
    quantum_yield: f64,
    wavelengths_nm: [f64; M],
    photon_fractions: [f64; M],
}

impl<const M: usize> FluorescenceSource<M> {
    /// Admit 1–16 strictly increasing emission nodes and nonnegative relative
    /// photon weights. The weights are explicitly normalized here. Each node
    /// must be at or beyond the ENTIRE declared excitation band (Stokes only).
    /// The temperature interval binds this fixed source, not a fitted law.
    pub fn try_new(
        source_identity: ContentHash,
        excitation_bounds_nm: (f64, f64),
        temperature_bounds_k: (f64, f64),
        quantum_yield: f64,
        wavelengths_nm: [f64; M],
        photon_weights: [f64; M],
    ) -> Result<Self, OpticalHeatingError> {
        let valid_interval =
            |(lo, hi): (f64, f64)| lo.is_finite() && hi.is_finite() && lo > 0.0 && hi >= lo;
        if !(1..=16).contains(&M)
            || source_identity == ContentHash([0; 32])
            || !valid_interval(excitation_bounds_nm)
            || !valid_interval(temperature_bounds_k)
            || !quantum_yield.is_finite()
            || !(0.0..=1.0).contains(&quantum_yield)
            || wavelengths_nm
                .iter()
                .any(|v| !v.is_finite() || *v < excitation_bounds_nm.1)
            || wavelengths_nm.windows(2).any(|v| v[0] >= v[1])
            || !photon_weights.into_iter().all(nonnegative)
        {
            return Err(invalid(
                "invalid sourced Stokes fluorescence spectrum or domain",
            ));
        }
        let maximum = photon_weights.into_iter().fold(0.0, f64::max);
        if maximum == 0.0 {
            return Err(invalid(
                "fluorescence requires positive total photon weight",
            ));
        }
        let mut total = CheckedSum::default();
        for weight in photon_weights {
            total.add(weight / maximum)?;
        }
        let total = total.value()?;
        let photon_fractions = photon_weights.map(|w| (w / maximum) / total);
        if photon_weights
            .iter()
            .zip(photon_fractions)
            .any(|(w, f)| *w > 0.0 && f == 0.0)
        {
            return Err(invalid(
                "a positive fluorescence photon weight is below representation",
            ));
        }
        Ok(Self {
            source_identity,
            excitation_bounds_nm,
            temperature_bounds_k,
            quantum_yield,
            wavelengths_nm,
            photon_fractions,
        })
    }

    /// Original declared source/condition identity.
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
    /// Emitted photons per absorbed pump photon.
    pub const fn quantum_yield(&self) -> f64 {
        self.quantum_yield
    }
    /// Vacuum wavelength of each discrete emission node, nm.
    pub const fn wavelengths_nm(&self) -> &[f64; M] {
        &self.wavelengths_nm
    }
    /// Normalized photon fractions; these are not energy fractions.
    pub const fn photon_fractions(&self) -> &[f64; M] {
        &self.photon_fractions
    }
    /// Inclusive source-supported pump band, vacuum nm.
    pub const fn excitation_bounds_nm(&self) -> (f64, f64) {
        self.excitation_bounds_nm
    }
    /// Inclusive temperature interval of the fixed fluorescence source, K.
    pub const fn temperature_bounds_k(&self) -> (f64, f64) {
        self.temperature_bounds_k
    }
}

/// Fixed escaping emission; each bin is already integrated power, in watts.
/// `residual_w` is absorbed pump power minus heat and escaping radiation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluorescencePower<const M: usize> {
    /// Escaping radiant watts in the source's emission-node order.
    pub escaped_w: [f64; M],
    /// Pump watts thermalized without photon emission.
    pub nonradiative_heat_w: f64,
    /// Positive absorbed-minus-emitted photon energy, W.
    pub stokes_heat_w: f64,
    /// Remaining absorbed-power partition roundoff, W.
    pub residual_w: f64,
}

/// Escaping radiation over the SAME accepted duration as the native heat audit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluorescenceEnergy<const M: usize> {
    /// Actual accepted native thermal duration, s.
    pub duration_s: f64,
    /// Escaped radiation at the supplied emission nodes, J.
    pub escaped_j: [f64; M],
    /// Nonradiative heat, already included once in the native heat audit, J.
    pub nonradiative_heat_j: f64,
    /// Stokes heat, already included once in the native heat audit, J.
    pub stokes_heat_j: f64,
    /// Integrated absorbed-power partition roundoff, J.
    pub residual_j: f64,
}

/// Prompt fluorescent heating on the original conservative enthalpy timeline.
/// Absorbed pump watts remain in `pump_power`; only the positive Stokes and
/// nonradiative losses enter its `thermal_input_w` and the native cell sources.
/// Candidate cloning rolls back heat and both radiation accounts together.
#[derive(Clone, Debug)]
pub struct FluorescentOpticalHeating<const N: usize, const M: usize> {
    source: FluorescenceSource<M>,
    heating: OpticalHeating<N>,
    power: FluorescencePower<M>,
    energy: FluorescenceEnergy<M>,
    cell_emission_w: Vec<(usize, [f64; M])>,
}

impl<const N: usize, const M: usize> FluorescentOpticalHeating<N, M> {
    /// Bind an original exact `AbsorbingGrid::deposit_beam_power` report to
    /// prompt fluorescence. Every occupied cell is checked for transparent
    /// emitted wavelengths, including cells outside the illuminated ray.
    /// Calling this constructor declares the supplied absorption-only grid is
    /// the COMPLETE optical environment with an open index-one boundary. Do not
    /// pass an absorption view stripped from a scattering field or a scene with
    /// reflecting/refracting surfaces: zero absorption alone cannot prove escape.
    /// `max_grid_cells` bounds that scan and the retained deposition/bindings.
    /// Public deposition reports undergo
    /// the same original conservation/mapping admission as `OpticalHeating`.
    /// The absolute watt tolerance also admits the fluorescence partition.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        medium: &AbsorbingGrid<'_>,
        source: FluorescenceSource<M>,
        deposition: BeamPowerDeposition<N>,
        cell_bindings: &[(usize, CellId)],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        power_tolerance_w: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        validate_fluorescent_medium(
            medium,
            &source,
            &deposition,
            cell_bindings.len(),
            max_grid_cells,
            power_tolerance_w,
            true,
            &mut checkpoint,
        )?;
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        let bindings: BTreeMap<_, _> = cell_bindings.iter().copied().collect();
        let mut escaped = [CheckedSum::default(); M];
        let mut nonradiative = CheckedSum::default();
        let mut stokes = CheckedSum::default();
        let mut total_heat = CheckedSum::default();
        let mut cell_emission_w = Vec::new();
        cell_emission_w
            .try_reserve_exact(deposition.cells.len())
            .map_err(|_| invalid("fluorescence cell-emission allocation refused"))?;
        for cell in &deposition.cells {
            checkpoint()?;
            if cell.cell_index >= medium.cell_count() {
                return Err(invalid(
                    "fluorescent deposition names a cell outside its medium",
                ));
            }
            let mut heat = CheckedSum::default();
            let mut cell_emission = [CheckedSum::default(); M];
            for lane in 0..N {
                let absorbed = cell.absorbed_power_w[lane];
                let wavelength = deposition.wavelengths_nm[lane];
                let dark = finite(absorbed * (1.0 - source.quantum_yield))?;
                nonradiative.add(dark)?;
                heat.add(dark)?;
                for bin in 0..M {
                    let emitted_nm = source.wavelengths_nm[bin];
                    let weight = absorbed * source.quantum_yield * source.photon_fractions[bin];
                    let emitted = ratio_power(weight, wavelength, emitted_nm)?;
                    escaped[bin].add(emitted)?;
                    cell_emission[bin].add(emitted)?;
                    // Positive loss preserves weak Stokes heating even when
                    // the escaped/absorbed power ratio rounds to exactly one.
                    let loss = ratio_power(weight, emitted_nm - wavelength, emitted_nm)?;
                    stokes.add(loss)?;
                    heat.add(loss)?;
                }
            }
            let mut emitted = [0.0; M];
            for (value, sum) in emitted.iter_mut().zip(cell_emission) {
                *value = sum.value()?;
            }
            cell_emission_w.push((cell.cell_index, emitted));
            let heat = heat.value()?;
            total_heat.add(heat)?;
            if let Some(&id) = bindings.get(&cell.cell_index) {
                sources.entry(id).or_default().add(heat)?;
            }
        }
        // Original admission owns report integrity, positive powers, complete
        // mappings and pump conservation; it publishes nothing on refusal.
        let mut heating = OpticalHeating::try_new(
            deposition,
            cell_bindings,
            thermal,
            power_tolerance_w,
            &mut checkpoint,
        )?;
        let mut power = FluorescencePower {
            escaped_w: [0.0; M],
            nonradiative_heat_w: nonradiative.value()?,
            stokes_heat_w: stokes.value()?,
            residual_w: 0.0,
        };
        let mut balance = CheckedSum::default();
        for absorbed in heating.power().absorbed_w {
            balance.add(absorbed)?;
        }
        for (value, sum) in power.escaped_w.iter_mut().zip(escaped) {
            *value = sum.value()?;
            balance.add(-*value)?;
        }
        balance.add(-total_heat.value()?)?;
        power.residual_w = balance.value()?;
        if power.residual_w.abs() > power_tolerance_w {
            return Err(invalid(
                "fluorescent heat and escaped radiation do not conserve absorbed power",
            ));
        }
        heating.sources_w = sources
            .into_iter()
            .map(|(id, sum)| Ok((id, sum.value()?)))
            .collect::<Result<_, OpticalHeatingError>>()?;
        heating.power.thermal_input_w = total_heat.value()?;
        heating.energy = heating.power.integrate(0.0)?;
        let result = Self {
            source,
            heating,
            power,
            energy: power.integrate(0.0)?,
            cell_emission_w,
        };
        result.check_temperatures(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }

    /// Immutable source and its declared photon/domain conventions.
    pub const fn source(&self) -> &FluorescenceSource<M> {
        &self.source
    }
    /// Original pump absorption and the reduced actual thermal-input watts.
    pub fn pump_power(&self) -> &OpticalPowerBudget<N> {
        self.heating.power()
    }
    /// Pump and native heat accounts integrated over accepted substeps only.
    pub fn pump_energy(&self) -> &OpticalEnergyBudget<N> {
        self.heating.energy()
    }
    /// Prepared escaping fluorescence and positive heat partition, W.
    pub const fn power(&self) -> &FluorescencePower<M> {
        &self.power
    }
    /// Fluorescence and heat partition on the accepted native clock, J.
    pub const fn energy(&self) -> &FluorescenceEnergy<M> {
        &self.energy
    }
    /// Original native thermal cells receiving the positive local heat, W.
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        self.heating.sources_w()
    }

    /// Isotropic emitted watts per original optical cell and discrete photon
    /// node. These are integrated powers, not a per-nanometre density. The
    /// spatial observer may distribute each power uniformly over that cell;
    /// this adds an explicit cell-uniform source approximation to the beam.
    pub fn cell_emission_w(&self) -> &[(usize, [f64; M])] {
        &self.cell_emission_w
    }

    fn check_temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalHeatingError> {
        let (lo, hi) = self.source.temperature_bounds_k;
        for &id in self.heating.sources_w().keys() {
            checkpoint()?;
            if !(lo..=hi).contains(&thermal.temperature_k(id)?) {
                return Err(invalid(
                    "native temperature exceeds the fixed fluorescence source domain",
                ));
            }
        }
        Ok(())
    }
}

impl<const N: usize, const M: usize> DepositionObserver for FluorescentOpticalHeating<N, M> {
    type Error = OpticalHeatingError;

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_temperatures(thermal, checkpoint)?;
        self.heating.sources(thermal, checkpoint)
    }
    // visit_sources reaches the heating owner's own, as sources() does. The step path stays
    // all-default while the owner has no step sources; its three hooks are forwarded together.
    fn visit_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.check_temperatures(thermal, checkpoint)?;
        self.heating.visit_sources(thermal, emit, checkpoint)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.check_temperatures(thermal, &mut || Ok(()))?;
        self.heating.event(event, thermal)
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.check_temperatures(thermal, checkpoint)?;
        self.heating
            .step(start_s, duration_s, old_temperatures_k, thermal, checkpoint)?;
        self.energy = self.power.integrate(self.heating.energy().duration_s)?;
        Ok(())
    }
}

impl<const M: usize> FluorescencePower<M> {
    fn integrate(&self, duration_s: f64) -> Result<FluorescenceEnergy<M>, OpticalHeatingError> {
        let mut escaped_j = [0.0; M];
        for (energy, power) in escaped_j.iter_mut().zip(self.escaped_w) {
            *energy = finite(power * duration_s)?;
        }
        Ok(FluorescenceEnergy {
            duration_s,
            escaped_j,
            nonradiative_heat_j: finite(self.nonradiative_heat_w * duration_s)?,
            stokes_heat_j: finite(self.stokes_heat_w * duration_s)?,
            residual_j: finite(self.residual_w * duration_s)?,
        })
    }
}

fn optical_error(_: VolumeTransportError) -> OpticalHeatingError {
    invalid("fluorescence requires the original source-supported absorption law")
}

/// Common pump/source admission. Only an owner with explicit secondary photon
/// transport may set `require_transparent_escape` to false; the prompt owner
/// always requires transparent escape. This does not construct an unphysical
/// prompt heating report for a nontransparent medium.
#[allow(clippy::too_many_arguments)]
pub(super) fn validate_fluorescent_medium<const N: usize, const M: usize>(
    medium: &AbsorbingGrid<'_>,
    source: &FluorescenceSource<M>,
    deposition: &BeamPowerDeposition<N>,
    binding_count: usize,
    max_grid_cells: usize,
    power_tolerance_w: f64,
    require_transparent_escape: bool,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), OpticalHeatingError> {
    checkpoint()?;
    if medium.cell_count() > max_grid_cells
        || deposition.cells.len() > medium.cell_count()
        || binding_count > deposition.cells.len()
        || deposition.medium_identity != medium.content_identity()
        || deposition.wavelength_bounds_nm != medium.wavelength_bounds_nm()
        || !nonnegative(power_tolerance_w)
        || deposition.scattered_power_w.iter().any(|p| *p != 0.0)
        || deposition.wavelengths_nm.iter().any(|l| {
            !l.is_finite()
                || !(source.excitation_bounds_nm.0..=source.excitation_bounds_nm.1).contains(l)
        })
    {
        return Err(invalid(
            "fluorescence requires the bound non-scattering pump and source domain",
        ));
    }
    for cell in &deposition.cells {
        checkpoint()?;
        if cell.cell_index >= medium.cell_count() {
            return Err(invalid(
                "fluorescent deposition names a cell outside its medium",
            ));
        }
    }
    let bounds = medium.wavelength_bounds_nm();
    if source
        .wavelengths_nm
        .iter()
        .any(|l| !(bounds.0..=bounds.1).contains(l))
    {
        return Err(invalid(
            "fluorescence wavelengths exceed the original absorption support",
        ));
    }
    if require_transparent_escape {
        // A uniform law is evaluated once. Indexed fields only repeat spectral
        // work when their occupied optical model changes.
        let mut previous = None;
        for cell in 0..medium.cell_count() {
            checkpoint()?;
            if medium.cell_density_kg_m3(cell) == Some(0.0) {
                continue;
            }
            let law = medium.absorption_at_cell(cell).map_err(optical_error)?;
            if previous.is_some_and(|old| core::ptr::eq(old, law)) {
                continue;
            }
            for &wavelength in &source.wavelengths_nm {
                if law.at(wavelength).map_err(optical_error)? != 0.0 {
                    return Err(invalid(
                        "self-absorbing fluorescence requires secondary transport",
                    ));
                }
            }
            previous = Some(law);
        }
    }
    Ok(())
}

// A tiny wavelength fraction can still carry representable watts. Evaluate
// ordinary ratios directly, preserving the near-zero Stokes difference; use
// the native deterministic log owner only when that ratio is subnormal.
pub(super) fn ratio_power(
    power: f64,
    numerator: f64,
    denominator: f64,
) -> Result<f64, OpticalHeatingError> {
    let ratio = numerator / denominator;
    if power == 0.0 || numerator == 0.0 || ratio >= f64::MIN_POSITIVE {
        finite(power * ratio)
    } else {
        finite(fs_math::det::exp(
            fs_math::det::ln(power) + fs_math::det::ln(numerator) - fs_math::det::ln(denominator),
        ))
    }
}
