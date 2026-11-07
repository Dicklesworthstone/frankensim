//! Explicit triplet-pair annihilation on the existing excitation timeline.
//!
//! Each authored law owns its actual active population volume. The renderer
//! never substitutes a thermal cell's volume or mass. The native composer
//! retains exact isolated pair kinetics and controls nonlinear evolution through
//! the original positive linear/photon owner; the optical owner retains
//! acceptance, cancellation and heat deposition.

use super::super::{Cell, LuminescenceSource, LuminescentOpticalHeating};
#[cfg(feature = "material-state")]
use super::RateLaw;
use super::{
    CoupledStep, HC_J_NM, SingletTripletFeedback, SingletTripletLimits, finite, invalid,
    preserve_checkpoint, sum,
};
use crate::volumes::heating::{CheckedSum, OpticalHeatingError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::DepositionError;
use fs_thermal::luminescence::transfer::SingletTripletDecay;
use fs_thermal::luminescence::transfer::recycling::{EmissionChannel, SingletTripletRecyclingPlan};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

pub use fs_thermal::luminescence::transfer::annihilation::TripletAnnihilation;
pub use fs_thermal::luminescence::transfer::recycling::annihilation::AnnihilationLimits;

/// A pair-event law for one existing optical excitation store.
/// The coefficient counts pair EVENTS: each event consumes two triplets.
/// A coefficient authored for triplet LOSS must be converted by its source;
/// this binding never guesses that convention or an active population volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhosphorescenceAnnihilation {
    pub grid_cell: usize,
    pub law: TripletAnnihilation,
}

#[derive(Debug)]
struct Settings {
    laws: BTreeMap<usize, TripletAnnihilation>,
    limits: AnnihilationLimits,
    has_nonzero_law: bool,
    identity: ContentHash,
}

#[derive(Debug)]
struct Prepared {
    duration_s: f64,
    plan: SingletTripletRecyclingPlan,
    channel_bins: Arc<[usize]>,
}

#[derive(Clone, Debug)]
pub(super) struct AnnihilationFeedback {
    settings: Arc<Settings>,
    prepared: Option<Arc<Prepared>>,
}

pub(super) struct NetworkAnnihilation {
    pub laws: Vec<Option<TripletAnnihilation>>,
    pub limits: AnnihilationLimits,
}

impl<const N: usize, const M: usize> LuminescentOpticalHeating<N, M> {
    /// Bind authored triplet-pair annihilation before the first advance.
    /// Singlet/triplet coupling must already exist. Bind photon recycling first
    /// when it must add unilluminated active receivers; this method never adds
    /// excitation stores and refuses duplicate or inactive grid-cell bindings.
    ///
    /// Active volume, pair-event coefficient and productive singlet probability
    /// come only from each supplied law. A productive pair must be able to pay
    /// the singlet energy event by event, even for a very small singlet yield.
    ///
    /// Limits apply to each independently evolved cell, or to the complete
    /// coupled photon network when recycling is bound. A zero coefficient uses
    /// the original linear fast path. Nonzero laws use the native controlled
    /// composition; exact isolated pair kinetics are retained by that owner.
    pub fn with_triplet_annihilation(
        mut self,
        bindings: &[PhosphorescenceAnnihilation],
        limits: AnnihilationLimits,
    ) -> Result<Self, OpticalHeatingError> {
        if self.pending.is_some()
            || self.energy.duration_s != 0.0
            || self
                .singlet_triplet
                .as_ref()
                .is_none_or(|coupling| coupling.annihilation.is_some())
        {
            return Err(invalid(
                "triplet annihilation requires an unadvanced singlet/triplet source",
            ));
        }
        if bindings.len() > self.cells.len() {
            return Err(invalid(
                "annihilation bindings exceed the existing excitation-store count",
            ));
        }
        let coupling = self
            .singlet_triplet
            .as_ref()
            .expect("admitted singlet/triplet source");
        let energies = source_energies(&self.source, coupling);
        let mut laws = BTreeMap::new();
        for binding in bindings {
            binding
                .law
                .advance(energies, [0.0; 2], 0.0)
                .map_err(|_| invalid("annihilation pair cannot supply its authored singlet"))?;
            if laws.insert(binding.grid_cell, binding.law).is_some() {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell));
            }
        }
        let mut missing: BTreeSet<_> = laws.keys().copied().collect();
        for cell in &self.cells {
            missing.remove(&cell.grid_cell);
        }
        if let Some(&grid_cell) = missing.first() {
            return Err(OpticalHeatingError::UnmappedCell(grid_cell));
        }
        // Admission belongs to the same native nonlinear owner used at runtime.
        // A stable, empty zero-duration model validates controls without making
        // spectra, nonlinear substeps or a finite lifetime necessary for a no-op.
        let admission = SingletTripletDecay::try_new(energies, [0.0; 2], [0.0; 2], [0.0; 2])
            .and_then(|model| {
                SingletTripletRecyclingPlan::prepare(
                    &[model],
                    &[],
                    &[],
                    0.0,
                    0.0,
                    coupling.source.limits,
                    || false,
                )
            })
            .and_then(|plan| {
                plan.advance_with_annihilation(&[[0.0; 2]], &[[0.0; 2]], &[None], limits, || false)
            });
        admission.map_err(|_| invalid("invalid triplet-annihilation solver limits"))?;
        let has_nonzero_law = laws.values().any(|law| law.pair_rate_m3_s() > 0.0);
        let identity = binding_identity(&laws, limits);
        self.singlet_triplet
            .as_mut()
            .expect("admitted singlet/triplet source")
            .annihilation = Some(AnnihilationFeedback {
            settings: Arc::new(Settings {
                laws,
                limits,
                has_nonzero_law,
                identity,
            }),
            prepared: None,
        });
        Ok(self)
    }

    /// Original explicit pair-event law for an active optical grid cell.
    pub fn triplet_annihilation_law(&self, grid_cell: usize) -> Option<TripletAnnihilation> {
        self.singlet_triplet
            .as_ref()?
            .annihilation
            .as_ref()?
            .settings
            .laws
            .get(&grid_cell)
            .copied()
    }

    /// Original nonlinear accuracy and work controls, if explicitly bound.
    pub fn triplet_annihilation_limits(&self) -> Option<AnnihilationLimits> {
        Some(
            self.singlet_triplet
                .as_ref()?
                .annihilation
                .as_ref()?
                .settings
                .limits,
        )
    }
}

impl AnnihilationFeedback {
    pub(super) fn identity(&self) -> ContentHash {
        self.settings.identity
    }

    pub(super) fn network_parameters<const N: usize>(
        &self,
        cells: &[Cell<N>],
    ) -> Option<NetworkAnnihilation> {
        if !self.settings.has_nonzero_law {
            return None;
        }
        Some(NetworkAnnihilation {
            laws: cells
                .iter()
                .map(|cell| self.settings.laws.get(&cell.grid_cell).copied())
                .collect(),
            limits: self.settings.limits,
        })
    }

    /// One independent cell uses two native states regardless of the total
    /// number of optical cells. Only immutable linear plans are cached here.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance<const M: usize>(
        &mut self,
        grid_cell: usize,
        singlet: &LuminescenceSource<M>,
        triplet: &LuminescenceSource<M>,
        model: SingletTripletDecay,
        stores_j: [f64; 2],
        pump_w: f64,
        duration_s: f64,
        linear_limits: SingletTripletLimits,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<CoupledStep<M>>, OpticalHeatingError> {
        let Some(law) = self
            .settings
            .laws
            .get(&grid_cell)
            .copied()
            .filter(|law| law.pair_rate_m3_s() > 0.0)
        else {
            return Ok(None);
        };
        if self.prepared.as_ref().is_none_or(|prepared| {
            prepared.plan.models() != core::slice::from_ref(&model)
                || prepared.duration_s != duration_s
        }) {
            let fraction_roundoff = spectrum_roundoff(singlet, triplet)?;
            let (plan, channel_bins) = match &self.prepared {
                Some(prepared) => {
                    let plan = preserve_checkpoint(checkpoint, |cancelled| {
                        if prepared.plan.models() == core::slice::from_ref(&model) {
                            prepared.plan.with_duration(duration_s, cancelled)
                        } else {
                            SingletTripletRecyclingPlan::prepare(
                                &[model],
                                prepared.plan.channels(),
                                &[],
                                duration_s,
                                fraction_roundoff,
                                linear_limits,
                                cancelled,
                            )
                        }
                    })?;
                    (plan, prepared.channel_bins.clone())
                }
                None => {
                    let mut channels = Vec::new();
                    let mut bins = Vec::new();
                    for (state, source) in [singlet, triplet].into_iter().enumerate() {
                        for (lane, &weight) in
                            source.fluorescence().photon_fractions().iter().enumerate()
                        {
                            checkpoint()?;
                            if weight == 0.0 {
                                continue;
                            }
                            if channels.len() >= linear_limits.maximum_edges {
                                return Err(invalid("annihilation photon-channel budget exceeded"));
                            }
                            channels.push(EmissionChannel {
                                source_cell: 0,
                                source_state: state,
                                photon_energy_j: finite(
                                    HC_J_NM / source.fluorescence().wavelengths_nm()[lane],
                                )?,
                                photon_fraction: weight,
                                escape_fraction: 1.0,
                                unresolved_fraction: 0.0,
                            });
                            bins.push(lane);
                        }
                    }
                    let plan = preserve_checkpoint(checkpoint, |cancelled| {
                        SingletTripletRecyclingPlan::prepare(
                            &[model],
                            &channels,
                            &[],
                            duration_s,
                            fraction_roundoff,
                            linear_limits,
                            cancelled,
                        )
                    })?;
                    (plan, Arc::from(bins))
                }
            };
            self.prepared = Some(Arc::new(Prepared {
                duration_s,
                plan,
                channel_bins,
            }));
        }
        let energies = model.excitation_energies_j();
        let populations = [[
            finite(stores_j[0] / energies[0])?,
            finite(stores_j[1] / energies[1])?,
        ]];
        let pump_rates = [[finite(pump_w / energies[0])?, 0.0]];
        let prepared = self
            .prepared
            .as_ref()
            .expect("prepared independent annihilation interval");
        let step = preserve_checkpoint(checkpoint, |cancelled| {
            prepared.plan.advance_with_annihilation(
                &populations,
                &pump_rates,
                &[Some(law)],
                self.settings.limits,
                cancelled,
            )
        })?;
        let mut emitted_energy_j = [0.0; M];
        for (index, &lane) in prepared.channel_bins.iter().enumerate() {
            checkpoint()?;
            emitted_energy_j[lane] =
                finite(emitted_energy_j[lane] + step.linear.emitted_energy_j[index])?;
        }
        Ok(Some(CoupledStep {
            stored_energy_j: sum(step.linear.stored_energy_j[0])?,
            state_stored_energy_j: step.linear.stored_energy_j[0],
            emitted_energy_j,
            nonradiative_heat_j: sum(step.linear.nonradiative_heat_j[0])?,
            relaxation_heat_j: sum(step.linear.relaxation_heat_j[0])?,
            crossing_heat_j: step.linear.crossing_heat_j[0],
            reverse_crossing_heat_absorbed_j: step.linear.reverse_crossing_heat_absorbed_j[0],
            annihilation_heat_j: step.annihilation_heat_j[0],
        }))
    }
}

fn source_energies<const M: usize>(
    singlet: &LuminescenceSource<M>,
    coupling: &SingletTripletFeedback<M>,
) -> [f64; 2] {
    #[cfg(feature = "material-state")]
    if let RateLaw::Material(source) = &coupling.source.rates {
        return source.original.decay().excitation_energies_j();
    }
    [
        HC_J_NM / singlet.excitation_wavelength_nm(),
        HC_J_NM / coupling.source.triplet.excitation_wavelength_nm(),
    ]
}

/// Only the already normalized source weights need a representation allowance:
/// every synthetic no-return transport row is exactly one. Retain the computed
/// normalization defect plus one unit-rounding allowance, without rescaling any
/// source fraction or reusing a nonlinear integration tolerance as photon data.
fn spectrum_roundoff<const M: usize>(
    singlet: &LuminescenceSource<M>,
    triplet: &LuminescenceSource<M>,
) -> Result<f64, OpticalHeatingError> {
    let mut error: f64 = 0.0;
    for source in [singlet, triplet] {
        let mut weights = CheckedSum::default();
        for &weight in source.fluorescence().photon_fractions() {
            weights.add(weight)?;
        }
        error = error.max((weights.value()? - 1.0).abs());
    }
    finite(error + f64::EPSILON)
}

fn binding_identity(
    laws: &BTreeMap<usize, TripletAnnihilation>,
    limits: AnnihilationLimits,
) -> ContentHash {
    let mut hash = DomainHasher::new("org.frankensim.render.triplet-annihilation.v1");
    hash.update(&(laws.len() as u64).to_le_bytes());
    hash.update(&limits.absolute_energy_tolerance_j.to_bits().to_le_bytes());
    hash.update(&limits.relative_tolerance.to_bits().to_le_bytes());
    hash.update(&(limits.maximum_substeps as u64).to_le_bytes());
    for (&grid_cell, law) in laws {
        hash.update(&(grid_cell as u64).to_le_bytes());
        for value in [
            law.pair_rate_m3_s(),
            law.active_volume_m3(),
            law.singlet_yield(),
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
    }
    hash.finalize()
}
