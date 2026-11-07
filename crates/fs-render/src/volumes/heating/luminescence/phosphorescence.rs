//! Authored singlet/triplet crossing on the existing optical/thermal timeline.
//!
//! The pump excites the original singlet source. Intersystem crossing transfers
//! actual populations to a separate triplet store; reverse crossing takes the
//! state-energy gap from the bath. Each state keeps its own photon spectrum.
//! This dilute first-order mechanism does not imply trap-mediated persistent
//! luminescence, ground-state depletion, annihilation or measured spin rates.

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

/// The thermal owner's existing positive-network work and accuracy limits.
pub use fs_thermal::luminescence::recycling::RecyclingLimits as SingletTripletLimits;

const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;

/// Explicit triplet source and forward/reverse intersystem-crossing laws.
/// Both scalar lifetimes describe intrinsic radiative plus dark removal;
/// crossing is added as a competing channel, never counted inside either one.
#[derive(Clone, Copy, Debug)]
pub struct SingletTripletSource<const M: usize> {
    triplet: LuminescenceSource<M>,
    forward_crossing_s: f64,
    reverse_crossing: ArrheniusQuenching,
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
        if !forward_crossing_s.is_finite() || forward_crossing_s < 0.0 {
            return Err(invalid("invalid singlet-to-triplet crossing rate"));
        }
        // Keep numerical admission with the existing evolution owner. Actual
        // source energies and temperature-dependent rates are checked on bind.
        SingletTripletDecay::try_new([2.0, 1.0], [1.0; 2], [0.0; 2], [0.0; 2])
            .and_then(|model| model.prepare(0.0, limits, || false))
            .map_err(|_| invalid("invalid singlet/triplet evolution limits"))?;
        Ok(Self {
            triplet,
            forward_crossing_s,
            reverse_crossing,
            limits,
        })
    }

    /// Original sourced triplet emission and intrinsic decay law.
    pub const fn triplet(&self) -> &LuminescenceSource<M> {
        &self.triplet
    }

    /// Authored fixed singlet-to-triplet rate, inverse seconds.
    pub const fn forward_crossing_s(&self) -> f64 {
        self.forward_crossing_s
    }

    /// Authored reverse-crossing thermal activation law.
    pub const fn reverse_crossing(&self) -> ArrheniusQuenching {
        self.reverse_crossing
    }
}

#[derive(Clone, Debug)]
struct Prepared {
    model: SingletTripletDecay,
    duration_s: f64,
    plan: SingletTripletPlan,
}

#[derive(Clone, Debug)]
pub(super) struct SingletTripletFeedback<const M: usize> {
    pub(super) source: SingletTripletSource<M>,
    prepared: Option<Arc<Prepared>>,
}

pub(super) struct CoupledStep<const M: usize> {
    pub stored_energy_j: f64,
    pub state_stored_energy_j: [f64; 2],
    pub emitted_energy_j: [f64; M],
    pub nonradiative_heat_j: f64,
    pub relaxation_heat_j: f64,
    pub crossing_heat_j: f64,
    pub reverse_crossing_heat_absorbed_j: f64,
}

impl<const N: usize, const M: usize> LuminescentOpticalHeating<N, M> {
    /// Bind a distinct triplet before advancing. The original source remains
    /// the singlet and the external pump populates it. Emission nodes must be
    /// identical so the existing transport operator can carry both spectra;
    /// photon weights may differ, including zero-weight nodes in either state.
    /// State energies must obey Es >= Et. Passive emitted-light absorption is
    /// supported; resonant recycling would require a different coupled network.
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
        self.singlet_triplet = Some(SingletTripletFeedback {
            source: coupling,
            prepared: None,
        });
        self.check_state(thermal, true, &mut checkpoint)?;
        Ok(self)
    }

    /// Authored triplet/crossing source, when explicitly bound.
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
        let mut hash = DomainHasher::new("org.frankensim.render.singlet-triplet-source.v1");
        for source in [&self.source, &coupling.source.triplet] {
            hash.update(source.fluorescence().source_identity().as_bytes());
            let (pump_lo, pump_hi) = source.fluorescence().excitation_bounds_nm();
            let (temp_lo, temp_hi) = source.fluorescence().temperature_bounds_k();
            for value in [
                source.excitation_wavelength_nm(),
                source.decay.lifetime_s(),
                source.decay.quantum_yield(),
                pump_lo,
                pump_hi,
                temp_lo,
                temp_hi,
            ] {
                hash.update(&value.to_bits().to_le_bytes());
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
        hash.update(&coupling.source.forward_crossing_s.to_bits().to_le_bytes());
        hash.update(
            &coupling
                .source
                .reverse_crossing
                .prefactor_rate_s()
                .to_bits()
                .to_le_bytes(),
        );
        hash.update(
            &coupling
                .source
                .reverse_crossing
                .activation_energy_j()
                .to_bits()
                .to_le_bytes(),
        );
        hash.finalize()
    }
}

impl<const M: usize> SingletTripletFeedback<M> {
    pub(super) fn model_at_temperature(
        &self,
        singlet: &LuminescenceSource<M>,
        temperature_k: f64,
    ) -> Result<SingletTripletDecay, OpticalHeatingError> {
        let decays = [
            singlet.decay_at_temperature(temperature_k)?,
            self.source.triplet.decay_at_temperature(temperature_k)?,
        ];
        let mut radiative = [0.0; 2];
        let mut nonradiative = [0.0; 2];
        for (i, decay) in decays.iter().enumerate() {
            radiative[i] = rate_power(1.0, decay.quantum_yield(), decay.lifetime_s())?;
            nonradiative[i] = rate_power(1.0, 1.0 - decay.quantum_yield(), decay.lifetime_s())?;
            if (decay.quantum_yield() > 0.0 && radiative[i] == 0.0)
                || (decay.quantum_yield() < 1.0 && nonradiative[i] == 0.0)
            {
                return Err(invalid(
                    "positive singlet/triplet removal rate is below representation",
                ));
            }
        }
        let reverse = self
            .source
            .reverse_crossing
            .rate_at_temperature(temperature_k)
            .map_err(|_| invalid("unrepresentable reverse intersystem crossing rate"))?;
        SingletTripletDecay::try_new(
            [
                HC_J_NM / singlet.excitation_wavelength_nm(),
                HC_J_NM / self.source.triplet.excitation_wavelength_nm(),
            ],
            radiative,
            nonradiative,
            [self.source.forward_crossing_s, reverse],
        )
        .map_err(|_| invalid("invalid singlet/triplet rates or energies"))
    }

    pub(super) fn advance(
        &mut self,
        singlet: &LuminescenceSource<M>,
        stores_j: [f64; 2],
        pump_w: f64,
        duration_s: f64,
        temperature_k: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CoupledStep<M>, OpticalHeatingError> {
        let model = self.model_at_temperature(singlet, temperature_k)?;
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

fn preserve_checkpoint<T>(
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
    result.map_err(|_| invalid("coupled singlet/triplet evolution refused"))
}
