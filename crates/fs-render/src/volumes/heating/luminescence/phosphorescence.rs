//! Singlet/triplet crossing on the existing optical/thermal timeline.
//!
//! The pump excites the original singlet source. Intersystem crossing transfers
//! actual populations to a separate triplet store; reverse crossing takes the
//! state-energy gap from the bath. Each state keeps its own photon spectrum.
//! The base crossing mechanism is dilute and first order. Photon recycling and
//! pair annihilation require separate explicit bindings; neither establishes
//! trap-mediated persistence, ground-state depletion or measured spin rates.

use super::{
    LuminescenceSource, LuminescentOpticalHeating, finite, invalid, partition, rate_power,
};
use crate::volumes::heating::{CheckedSum, OpticalHeatingError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_thermal::luminescence::ArrheniusQuenching;
use fs_thermal::luminescence::transfer::{
    SingletTripletDecay, SingletTripletError, SingletTripletPlan,
};
use std::sync::Arc;
#[cfg(feature = "material-state")]
use std::{collections::BTreeMap, sync::Mutex};

mod annihilation;
#[cfg(feature = "material-state")]
mod material;
use annihilation::AnnihilationFeedback;
pub use annihilation::{AnnihilationLimits, PhosphorescenceAnnihilation, TripletAnnihilation};
mod photon_recycling;
pub use photon_recycling::PhosphorescenceReceiver;
use photon_recycling::PhotonRecyclingFeedback;

/// The thermal owner's existing positive-network work and accuracy limits.
pub use fs_thermal::luminescence::recycling::RecyclingLimits as SingletTripletLimits;

const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;

#[derive(Clone, Debug)]
enum RateLaw {
    Intrinsic {
        forward_crossing_s: f64,
        reverse_crossing: ArrheniusQuenching,
    },
    #[cfg(feature = "material-state")]
    Material(Arc<material::MaterialTemperaturePhosphor>),
}

/// Explicit triplet spectrum and competing singlet/triplet transition laws.
/// Authored scalar lifetimes exclude crossing. Original material laws instead
/// supply all six native rates directly, including exactly stable states.
#[derive(Clone, Debug)]
pub struct SingletTripletSource<const M: usize> {
    triplet: LuminescenceSource<M>,
    rates: RateLaw,
    expected_singlet: Option<LuminescenceSource<M>>,
    limits: SingletTripletLimits,
}

impl<const M: usize> SingletTripletSource<M> {
    /// Supply a triplet spectrum/lifetime, a fixed forward rate [1/s], and
    /// an authored Arrhenius reverse rate. A zero barrier supplies a constant
    /// reverse rate; a zero prefactor disables it. No detailed balance or
    /// oxygen concentration is inferred from these inputs.
    pub fn try_new(
        triplet: LuminescenceSource<M>,
        forward_crossing_s: f64,
        reverse_crossing: ArrheniusQuenching,
        limits: SingletTripletLimits,
    ) -> Result<Self, OpticalHeatingError> {
        if !forward_crossing_s.is_finite() || forward_crossing_s < 0.0 || triplet.decay.is_none() {
            return Err(invalid(
                "authored crossing requires a finite intrinsic triplet decay",
            ));
        }
        // Keep numerical admission with the existing evolution owner. Actual
        // source energies and temperature-dependent rates are checked on bind.
        SingletTripletDecay::try_new([2.0, 1.0], [1.0; 2], [0.0; 2], [0.0; 2])
            .and_then(|model| model.prepare(0.0, limits, || false))
            .map_err(|_| invalid("invalid singlet/triplet evolution limits"))?;
        Ok(Self {
            triplet,
            rates: RateLaw::Intrinsic {
                forward_crossing_s,
                reverse_crossing,
            },
            expected_singlet: None,
            limits,
        })
    }

    /// Original triplet spectrum; native material rates need no scalar lifetime.
    pub const fn triplet(&self) -> &LuminescenceSource<M> {
        &self.triplet
    }

    /// Authored fixed singlet-to-triplet rate, inverse seconds. Native
    /// material rates refresh at the actual thermal temperature instead.
    pub fn forward_crossing_s(&self) -> Option<f64> {
        match &self.rates {
            RateLaw::Intrinsic {
                forward_crossing_s, ..
            } => Some(*forward_crossing_s),
            #[cfg(feature = "material-state")]
            RateLaw::Material(_) => None,
        }
    }

    /// Authored reverse-crossing activation law, absent for native material
    /// rates. A sampled material curve is never replaced with an Arrhenius fit.
    pub fn reverse_crossing(&self) -> Option<ArrheniusQuenching> {
        match &self.rates {
            RateLaw::Intrinsic {
                reverse_crossing, ..
            } => Some(*reverse_crossing),
            #[cfg(feature = "material-state")]
            RateLaw::Material(_) => None,
        }
    }

    /// Numerical work and accuracy limits for the native excitation evolution.
    pub const fn evolution_limits(&self) -> SingletTripletLimits {
        self.limits
    }

    /// Identity of the authored crossing law or original material selection.
    /// Spectral descriptors and intrinsic decay remain the caller's separate
    /// source inputs; this identity alone does not confer measurement authority.
    pub fn rate_source_identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.render.phosphorescence-rate-source.v1");
        match &self.rates {
            RateLaw::Intrinsic {
                forward_crossing_s,
                reverse_crossing,
            } => {
                hash.update(&[0]);
                for value in [
                    *forward_crossing_s,
                    reverse_crossing.prefactor_rate_s(),
                    reverse_crossing.activation_energy_j(),
                ] {
                    hash.update(&value.to_bits().to_le_bytes());
                }
            }
            #[cfg(feature = "material-state")]
            RateLaw::Material(source) => {
                hash.update(&[1]);
                hash.update(source.identity.as_bytes());
            }
        }
        hash.finalize()
    }

    /// Resolve the shared constitutive rates at an actual thermal temperature.
    /// Both dilute pumping and finite-ground-population pumping use this same
    /// source, including original material claims when explicitly bound.
    /// This query changes no excitation history; coefficients are fixed over
    /// the caller's accepted interval.
    pub fn model_at_temperature(
        &self,
        singlet: &LuminescenceSource<M>,
        temperature_k: f64,
    ) -> Result<SingletTripletDecay, OpticalHeatingError> {
        self.admit_singlet(singlet)?;
        singlet.check_temperature(temperature_k)?;
        self.triplet.check_temperature(temperature_k)?;
        #[cfg(feature = "material-state")]
        if let RateLaw::Material(source) = &self.rates {
            return source.model_at_temperature(temperature_k);
        }
        let (forward_crossing_s, reverse_crossing) = match &self.rates {
            RateLaw::Intrinsic {
                forward_crossing_s,
                reverse_crossing,
            } => (*forward_crossing_s, *reverse_crossing),
            #[cfg(feature = "material-state")]
            RateLaw::Material(_) => unreachable!("material rate dispatch returned above"),
        };
        let decays = [
            singlet.decay_at_temperature(temperature_k)?,
            self.triplet.decay_at_temperature(temperature_k)?,
        ];
        let mut radiative = [0.0; 2];
        let mut nonradiative = [0.0; 2];
        for (i, decay) in decays.iter().enumerate() {
            radiative[i] = rate_power(1.0, decay.quantum_yield(), decay.lifetime_s())?;
            nonradiative[i] = rate_power(1.0, decay.nonradiative_yield(), decay.lifetime_s())?;
            if (decay.quantum_yield() > 0.0 && radiative[i] == 0.0)
                || (decay.nonradiative_yield() > 0.0 && nonradiative[i] == 0.0)
            {
                return Err(invalid(
                    "positive singlet/triplet removal rate is below representation",
                ));
            }
        }
        let reverse = reverse_crossing
            .rate_at_temperature(temperature_k)
            .map_err(|_| invalid("unrepresentable reverse intersystem crossing rate"))?;
        SingletTripletDecay::try_new(
            [
                HC_J_NM / singlet.excitation_wavelength_nm(),
                HC_J_NM / self.triplet.excitation_wavelength_nm(),
            ],
            radiative,
            nonradiative,
            [forward_crossing_s, reverse],
        )
        .map_err(|_| invalid("invalid singlet/triplet rates or energies"))
    }

    fn admit_singlet(&self, singlet: &LuminescenceSource<M>) -> Result<(), OpticalHeatingError> {
        if let Some(expected) = &self.expected_singlet {
            if singlet != expected {
                return Err(invalid(
                    "native phosphorescence requires its original singlet spectrum and energy, without additional quenching",
                ));
            }
        } else if singlet.decay.is_none() {
            return Err(invalid(
                "authored crossing requires an intrinsic singlet decay",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct Prepared {
    model: SingletTripletDecay,
    duration_s: f64,
    plan: SingletTripletPlan,
}

#[cfg(feature = "material-state")]
#[derive(Debug, Default)]
struct TemperatureCache {
    interval_start_bits: Option<u64>,
    models: BTreeMap<u64, SingletTripletDecay>,
}

#[derive(Debug)]
pub(super) struct SingletTripletFeedback<const M: usize> {
    pub(super) source: SingletTripletSource<M>,
    prepared: Option<Arc<Prepared>>,
    photon_recycling: Option<PhotonRecyclingFeedback>,
    annihilation: Option<AnnihilationFeedback>,
    #[cfg(feature = "material-state")]
    temperatures: Mutex<TemperatureCache>,
    #[cfg(feature = "material-state")]
    maximum_cached_temperatures: usize,
}

impl<const M: usize> Clone for SingletTripletFeedback<M> {
    fn clone(&self) -> Self {
        Self {
            source: self.source.clone(),
            prepared: self.prepared.clone(),
            photon_recycling: self.photon_recycling.clone(),
            annihilation: self.annihilation.clone(),
            // Only immutable source/propagator data are shared. A candidate
            // starts with its own empty cache, independent of accepted stores.
            #[cfg(feature = "material-state")]
            temperatures: Mutex::new(TemperatureCache::default()),
            #[cfg(feature = "material-state")]
            maximum_cached_temperatures: self.maximum_cached_temperatures,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct CoupledStep<const M: usize> {
    pub stored_energy_j: f64,
    pub state_stored_energy_j: [f64; 2],
    pub emitted_energy_j: [f64; M],
    pub nonradiative_heat_j: f64,
    pub relaxation_heat_j: f64,
    pub crossing_heat_j: f64,
    pub reverse_crossing_heat_absorbed_j: f64,
    pub annihilation_heat_j: f64,
}

impl<const N: usize, const M: usize> LuminescentOpticalHeating<N, M> {
    /// Bind a distinct triplet before advancing. The original source remains
    /// the singlet and the external pump populates it. Emission nodes must be
    /// identical so the existing transport operator can carry both spectra;
    /// photon weights may differ, including zero-weight nodes in either state.
    /// State energies must obey Es >= Et. Passive emitted-light absorption is
    /// supported; explicit spin-state photon recycling can be bound afterward.
    pub fn with_singlet_triplet(
        mut self,
        coupling: SingletTripletSource<M>,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        if self.pending.is_some()
            || self.energy.duration_s != 0.0
            || self.singlet_triplet.is_some()
            || self.recycling.is_some()
            || self.source.fluorescence().wavelengths_nm()
                != coupling.triplet.fluorescence().wavelengths_nm()
            || self.source.excitation_wavelength_nm() > coupling.triplet.excitation_wavelength_nm()
        {
            return Err(invalid(
                "singlet/triplet coupling requires an unadvanced supported common-grid source",
            ));
        }
        coupling.admit_singlet(&self.source)?;
        self.singlet_triplet = Some(SingletTripletFeedback {
            source: coupling,
            prepared: None,
            photon_recycling: None,
            annihilation: None,
            #[cfg(feature = "material-state")]
            temperatures: Mutex::new(TemperatureCache::default()),
            #[cfg(feature = "material-state")]
            maximum_cached_temperatures: self.cells.len().saturating_mul(2).max(1),
        });
        self.check_state(thermal, true, &mut checkpoint)?;
        Ok(self)
    }

    /// Bound triplet spectrum and crossing/rate source.
    pub fn singlet_triplet_source(&self) -> Option<&SingletTripletSource<M>> {
        self.singlet_triplet
            .as_ref()
            .map(|feedback| &feedback.source)
    }

    /// Attribution for the complete emitting source, including both spectra,
    /// state energies, intrinsic lifetimes, quenching and crossing laws.
    /// This identifies authored data; it does not confer measurement authority.
    pub fn emission_source_identity(&self) -> ContentHash {
        let Some(coupling) = &self.singlet_triplet else {
            return self.source.fluorescence().source_identity();
        };
        let mut hash = DomainHasher::new("org.frankensim.render.singlet-triplet-source.v2");
        for source in [&self.source, &coupling.source.triplet] {
            hash.update(source.fluorescence().source_identity().as_bytes());
            let (pump_lo, pump_hi) = source.fluorescence().excitation_bounds_nm();
            let (temp_lo, temp_hi) = source.fluorescence().temperature_bounds_k();
            for value in [
                source.excitation_wavelength_nm(),
                pump_lo,
                pump_hi,
                temp_lo,
                temp_hi,
            ] {
                hash.update(&value.to_bits().to_le_bytes());
            }
            match source.decay {
                Some(decay) => {
                    hash.update(&[1]);
                    hash.update(&decay.lifetime_s().to_bits().to_le_bytes());
                    hash.update(&decay.quantum_yield().to_bits().to_le_bytes());
                }
                None => hash.update(&[0]),
            }
            for value in source
                .fluorescence()
                .wavelengths_nm()
                .iter()
                .chain(source.fluorescence().photon_fractions())
            {
                hash.update(&value.to_bits().to_le_bytes());
            }
            match source.quenching {
                Some(law) => {
                    hash.update(&[1]);
                    hash.update(&law.prefactor_rate_s().to_bits().to_le_bytes());
                    hash.update(&law.activation_energy_j().to_bits().to_le_bytes());
                }
                None => hash.update(&[0]),
            }
        }
        match &coupling.source.rates {
            RateLaw::Intrinsic {
                forward_crossing_s,
                reverse_crossing,
            } => {
                hash.update(&[0]);
                hash.update(&forward_crossing_s.to_bits().to_le_bytes());
                hash.update(&reverse_crossing.prefactor_rate_s().to_bits().to_le_bytes());
                hash.update(
                    &reverse_crossing
                        .activation_energy_j()
                        .to_bits()
                        .to_le_bytes(),
                );
            }
            #[cfg(feature = "material-state")]
            RateLaw::Material(source) => {
                hash.update(&[1]);
                hash.update(source.identity.as_bytes());
            }
        }
        let source_identity = hash.finalize();
        let source_identity = match &coupling.photon_recycling {
            Some(recycling) => {
                let mut hash =
                    DomainHasher::new("org.frankensim.render.singlet-triplet-recycling-source.v1");
                hash.update(source_identity.as_bytes());
                hash.update(recycling.identity().as_bytes());
                hash.finalize()
            }
            None => source_identity,
        };
        match &coupling.annihilation {
            Some(annihilation) => {
                let mut hash = DomainHasher::new(
                    "org.frankensim.render.singlet-triplet-annihilation-source.v1",
                );
                hash.update(source_identity.as_bytes());
                hash.update(annihilation.identity().as_bytes());
                hash.finalize()
            }
            None => source_identity,
        }
    }
}

impl<const M: usize> SingletTripletFeedback<M> {
    /// Memoized constitutive data have no physical state and need no rollback.
    /// Retrying at the same accepted start retains already admitted temperatures.
    pub(super) fn begin_interval(&mut self, start_s: f64) {
        #[cfg(feature = "material-state")]
        {
            if self.temperatures.is_poisoned() {
                self.temperatures = Mutex::new(TemperatureCache::default());
            }
            let cache = self
                .temperatures
                .get_mut()
                .expect("unpoisoned source cache");
            if cache.interval_start_bits != Some(start_s.to_bits()) {
                cache.models.clear();
                cache.interval_start_bits = Some(start_s.to_bits());
            }
        }
        #[cfg(not(feature = "material-state"))]
        let _ = start_s;
    }

    pub(super) fn model_at_temperature(
        &self,
        singlet: &LuminescenceSource<M>,
        temperature_k: f64,
    ) -> Result<SingletTripletDecay, OpticalHeatingError> {
        self.source.admit_singlet(singlet)?;
        singlet.check_temperature(temperature_k)?;
        self.source.triplet.check_temperature(temperature_k)?;
        #[cfg(feature = "material-state")]
        if let RateLaw::Material(source) = &self.source.rates {
            let mut cache = self
                .temperatures
                .lock()
                .map_err(|_| invalid("material temperature cache was interrupted"))?;
            let key = temperature_k.to_bits();
            if let Some(model) = cache.models.get(&key) {
                return Ok(*model);
            }
            let model = source.model_at_temperature(temperature_k)?;
            if cache.models.len() == self.maximum_cached_temperatures {
                cache.models.clear();
            }
            cache.models.insert(key, model);
            return Ok(model);
        }
        self.source.model_at_temperature(singlet, temperature_k)
    }

    pub(super) fn advance(
        &mut self,
        singlet: &LuminescenceSource<M>,
        grid_cell: usize,
        stores_j: [f64; 2],
        pump_w: f64,
        duration_s: f64,
        temperature_k: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CoupledStep<M>, OpticalHeatingError> {
        let model = self.model_at_temperature(singlet, temperature_k)?;
        if let Some(annihilation) = &mut self.annihilation
            && let Some(step) = annihilation.advance(
                grid_cell,
                singlet,
                &self.source.triplet,
                model,
                stores_j,
                pump_w,
                duration_s,
                self.source.limits,
                checkpoint,
            )?
        {
            return Ok(step);
        }
        if self
            .prepared
            .as_ref()
            .is_none_or(|cached| cached.model != model || cached.duration_s != duration_s)
        {
            let plan = preserve_checkpoint(checkpoint, |cancelled| {
                model.prepare(duration_s, self.source.limits, cancelled)
            })?;
            self.prepared = Some(Arc::new(Prepared {
                model,
                duration_s,
                plan,
            }));
        }
        let energies = model.excitation_energies_j();
        let populations = [
            finite(stores_j[0] / energies[0])?,
            finite(stores_j[1] / energies[1])?,
        ];
        let pump = [finite(pump_w / energies[0])?, 0.0];
        let prepared = self
            .prepared
            .as_ref()
            .expect("prepared singlet/triplet interval");
        let step = preserve_checkpoint(checkpoint, |cancelled| {
            prepared
                .plan
                .advance_with_cancel(populations, pump, cancelled)
        })?;
        let (singlet_emission, singlet_heat) =
            partition(singlet, step.radiative_decay_energy_j[0])?;
        let (triplet_emission, triplet_heat) =
            partition(&self.source.triplet, step.radiative_decay_energy_j[1])?;
        let mut emitted = [0.0; M];
        for i in 0..M {
            emitted[i] = finite(singlet_emission[i] + triplet_emission[i])?;
        }
        Ok(CoupledStep {
            stored_energy_j: sum(step.stored_energy_j)?,
            state_stored_energy_j: step.stored_energy_j,
            emitted_energy_j: emitted,
            nonradiative_heat_j: sum(step.nonradiative_heat_j)?,
            relaxation_heat_j: finite(singlet_heat + triplet_heat)?,
            crossing_heat_j: step.crossing_heat_j,
            reverse_crossing_heat_absorbed_j: step.reverse_crossing_heat_absorbed_j,
            annihilation_heat_j: 0.0,
        })
    }

    pub(super) fn emission_power_w(
        &self,
        singlet: &LuminescenceSource<M>,
        stores_j: [f64; 2],
        temperature_k: f64,
    ) -> Result<[f64; M], OpticalHeatingError> {
        let rates = self
            .model_at_temperature(singlet, temperature_k)?
            .radiative_rates_s();
        let (spectrum_s, _) = partition(singlet, finite(stores_j[0] * rates[0])?)?;
        let (spectrum_t, _) = partition(&self.source.triplet, finite(stores_j[1] * rates[1])?)?;
        let mut emitted = [0.0; M];
        for i in 0..M {
            emitted[i] = finite(spectrum_s[i] + spectrum_t[i])?;
        }
        Ok(emitted)
    }
}

fn sum(values: [f64; 2]) -> Result<f64, OpticalHeatingError> {
    let mut total = CheckedSum::default();
    for value in values {
        total.add(value)?;
    }
    total.value()
}

pub(super) fn preserve_checkpoint<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, SingletTripletError>,
) -> Result<T, OpticalHeatingError> {
    let mut refused = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            refused = Some(error);
            true
        }
    });
    if let Some(error) = refused {
        return Err(error.into());
    }
    result.map_err(OpticalHeatingError::SingletTriplet)
}
