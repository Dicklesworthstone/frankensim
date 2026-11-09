//! Authored or original-material pair kinetics on the finite-inventory clock.
//!
//! The active volume belongs to the supplied law. Neither thermal mass nor
//! beam area substitutes for it. Only immutable constitutive snapshots are
//! cached; every pair event acts on the same accepted ground/singlet/triplet
//! inventory as pumping, crossing and decay.

use super::{
    DepositionError, DepositionSimulation, DomainHasher, OpticalHeatingError,
    SaturableLuminescentHeating, SaturableSingletTripletPlan, SingletTripletDecay, invalid,
    preserve_checkpoint,
};
use std::sync::Arc;

pub use super::super::phosphorescence::{PhosphorescenceAnnihilation, TripletAnnihilation};
pub use fs_thermal::luminescence::transfer::saturation::annihilation::AnnihilationLimits;

#[derive(Clone, Debug)]
enum PairSource {
    Authored(TripletAnnihilation),
    #[cfg(feature = "material-state")]
    Material {
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::annihilation::TripletAnnihilationStatePoint,
        active_volume_m3: f64,
    },
}

/// One original emitter cell's pair-event coefficient, active volume and yield.
/// Material sources refresh at that cell's actual native temperature while
/// retaining their exact claims and all other external conditions.
#[derive(Clone, Debug)]
pub struct EmitterTripletAnnihilation {
    grid_cell: usize,
    source: PairSource,
}

impl EmitterTripletAnnihilation {
    pub fn authored(grid_cell: usize, law: TripletAnnihilation) -> Self {
        Self {
            grid_cell,
            source: PairSource::Authored(law),
        }
    }

    /// Bind original molecular pair data with an explicitly owned absolute
    /// temperature coordinate. The supplied physical active volume contains
    /// the complete cell's triplet count; it is never inferred from density.
    #[cfg(feature = "material-state")]
    pub fn from_material(
        grid_cell: usize,
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::annihilation::TripletAnnihilationStatePoint,
        active_volume_m3: f64,
    ) -> Result<Self, OpticalHeatingError> {
        if source.temperature_axis().is_none() {
            return Err(invalid(
                "thermal pair source requires an owned absolute-temperature axis",
            ));
        }
        source
            .at_query_point(&card, &source.resolved().query())
            .map_err(|_| invalid("triplet pair data differ from the original material card"))?;
        source
            .law(active_volume_m3)
            .map_err(|_| invalid("invalid triplet-pair active volume or source"))?;
        Ok(Self {
            grid_cell,
            source: PairSource::Material {
                card,
                source,
                active_volume_m3,
            },
        })
    }

    pub const fn grid_cell(&self) -> usize {
        self.grid_cell
    }

    fn evaluate(&self, temperature_k: f64) -> Result<TripletAnnihilation, OpticalHeatingError> {
        match &self.source {
            PairSource::Authored(law) => {
                let _ = temperature_k;
                Ok(*law)
            }
            #[cfg(feature = "material-state")]
            PairSource::Material {
                card,
                source,
                active_volume_m3,
            } => source
                .at_temperature(card, Some(temperature_k))
                .and_then(|point| point.law(*active_volume_m3))
                .map_err(|_| {
                    invalid("native temperature exceeds original triplet-pair source support")
                }),
        }
    }

    fn hash_into(&self, hash: &mut DomainHasher) {
        hash.update(&(self.grid_cell as u64).to_le_bytes());
        match &self.source {
            PairSource::Authored(law) => {
                hash.update(&[0]);
                for value in [
                    law.pair_rate_m3_s(),
                    law.active_volume_m3(),
                    law.singlet_yield(),
                ] {
                    hash.update(&value.to_bits().to_le_bytes());
                }
            }
            #[cfg(feature = "material-state")]
            PairSource::Material {
                source,
                active_volume_m3,
                ..
            } => {
                hash.update(&[1]);
                hash.update(source.identity().as_bytes());
                hash.update(&active_volume_m3.to_bits().to_le_bytes());
            }
        }
    }
}

#[derive(Debug)]
struct Settings {
    /// Same index order as the original owner's unique emitter inventory.
    sources: Vec<Option<EmitterTripletAnnihilation>>,
    limits: AnnihilationLimits,
}

#[derive(Clone, Debug)]
pub(super) struct AnnihilationFeedback {
    settings: Arc<Settings>,
    prepared: Option<Arc<PreparedAnnihilation>>,
}

#[derive(Debug)]
pub(super) struct PreparedAnnihilation {
    temperatures: Vec<Option<u64>>,
    pub laws: Vec<Option<TripletAnnihilation>>,
    pub limits: AnnihilationLimits,
}

impl<const M: usize> SaturableLuminescentHeating<M> {
    /// Add explicit pair-event laws to the SAME finite inventories before the
    /// first advance. Unbound cells retain their existing linear kinetics.
    /// A productive pair must pay its singlet energy event by event. Pump
    /// depletion and pair recovery are solved together, including in the dark.
    ///
    /// The original law's k counts pair events, each consuming two triplets.
    /// Its active volume and productive yield are retained without conversion.
    /// Zero coefficients preserve the existing linear fast path exactly.
    pub fn with_triplet_annihilation(
        self,
        bindings: &[PhosphorescenceAnnihilation],
        limits: AnnihilationLimits,
    ) -> Result<Self, OpticalHeatingError> {
        if bindings.len() > self.emitters.len() {
            return Err(invalid("pair bindings exceed the finite emitter count"));
        }
        let sources: Vec<_> = bindings
            .iter()
            .map(|binding| EmitterTripletAnnihilation::authored(binding.grid_cell, binding.law))
            .collect();
        self.install_annihilation(&sources, limits, &mut || Ok(()))
    }

    /// Bind authored and/or source-bound pair laws, resolving every admitted
    /// source at its original thermal cell's current temperature. Both source
    /// domains and eventwise energy admission are checked again at every
    /// accepted endpoint before any population or energy account is published.
    pub fn with_triplet_annihilation_sources(
        self,
        bindings: &[EmitterTripletAnnihilation],
        thermal: &DepositionSimulation<'_>,
        limits: AnnihilationLimits,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        self.check_state(thermal, true, &mut checkpoint)?;
        let mut owner = self.install_annihilation(bindings, limits, &mut checkpoint)?;
        let prepared = owner
            .annihilation
            .as_ref()
            .expect("bound pair laws")
            .prepare(&owner, thermal, &mut checkpoint)?;
        checkpoint()?;
        owner
            .annihilation
            .as_mut()
            .expect("bound pair laws")
            .accept(prepared);
        Ok(owner)
    }

    /// Last accepted constitutive law for this original cell. Authored laws
    /// are available immediately; material laws are resolved at binding and
    /// replaced only after a successful native endpoint.
    pub fn triplet_annihilation_law(&self, grid_cell: usize) -> Option<TripletAnnihilation> {
        let index = self
            .emitters
            .binary_search_by_key(&grid_cell, |e| e.binding.grid_cell)
            .ok()?;
        let feedback = self.annihilation.as_ref()?;
        if let Some(prepared) = &feedback.prepared {
            return prepared.laws[index];
        }
        match &feedback.settings.sources[index].as_ref()?.source {
            PairSource::Authored(law) => Some(*law),
            #[cfg(feature = "material-state")]
            PairSource::Material { .. } => None,
        }
    }

    pub fn triplet_annihilation_limits(&self) -> Option<AnnihilationLimits> {
        self.annihilation
            .as_ref()
            .map(|feedback| feedback.settings.limits)
    }

    fn install_annihilation(
        mut self,
        bindings: &[EmitterTripletAnnihilation],
        limits: AnnihilationLimits,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        if self.pending.is_some()
            || self.energy.duration_s != 0.0
            || self.annihilation.is_some()
            || bindings.len() > self.emitters.len()
        {
            return Err(invalid(
                "triplet pairs require unadvanced, uniquely bound finite inventories",
            ));
        }
        // Delegate control admission to the original nonlinear owner even for
        // empty/zero-law bindings; no source lifetime or numerical rate is invented.
        preserve_checkpoint(checkpoint, |cancelled| {
            let model = SingletTripletDecay::try_new([1.0; 2], [0.0; 2], [0.0; 2], [0.0; 2])?;
            let plan = SaturableSingletTripletPlan::prepare(
                model,
                [0.0; 2],
                [0.0; 2],
                0.0,
                self.triplet.evolution_limits(),
                &mut *cancelled,
            )?;
            plan.advance_with_annihilation([0.0; 3], None, limits, cancelled)
        })?;
        let mut sources = vec![None; self.emitters.len()];
        for binding in bindings {
            checkpoint()?;
            let index = self
                .emitters
                .binary_search_by_key(&binding.grid_cell, |e| e.binding.grid_cell)
                .map_err(|_| OpticalHeatingError::UnmappedCell(binding.grid_cell))?;
            if sources[index].is_some() {
                return Err(OpticalHeatingError::DuplicateCell(binding.grid_cell));
            }
            match &binding.source {
                PairSource::Authored(law) => {
                    validate_pair(
                        *law,
                        self.emitters[index]
                            .prepared
                            .as_ref()
                            .expect("prepared initial finite inventory")
                            .model()
                            .excitation_energies_j(),
                    )?;
                }
                #[cfg(feature = "material-state")]
                PairSource::Material { .. } => {}
            }
            sources[index] = Some(binding.clone());
        }
        let mut hash = DomainHasher::new("org.frankensim.render.saturable-triplet-pairs.v1");
        hash.update(self.source_identity.as_bytes());
        hash.update(&(bindings.len() as u64).to_le_bytes());
        for source in sources.iter().flatten() {
            source.hash_into(&mut hash);
        }
        for value in [
            limits.absolute_energy_tolerance_j,
            limits.relative_tolerance,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&(limits.maximum_substeps as u64).to_le_bytes());
        checkpoint()?;
        self.source_identity = hash.finalize();
        self.annihilation = Some(AnnihilationFeedback {
            settings: Arc::new(Settings { sources, limits }),
            prepared: None,
        });
        Ok(self)
    }
}

impl AnnihilationFeedback {
    pub(super) fn prepare<const M: usize>(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Arc<PreparedAnnihilation>, OpticalHeatingError> {
        if let Some(prepared) = &self.prepared {
            let mut matches = true;
            for (emitter, expected) in owner.emitters.iter().zip(&prepared.temperatures) {
                checkpoint()?;
                if let Some(expected) = expected
                    && thermal
                        .temperature_k(emitter.binding.thermal_cell)?
                        .to_bits()
                        != *expected
                {
                    matches = false;
                    break;
                }
            }
            if matches {
                return Ok(Arc::clone(prepared));
            }
        }
        let mut temperatures = vec![None; owner.emitters.len()];
        let mut laws = vec![None; owner.emitters.len()];
        let mut cached_model = None;
        for (index, source) in self.settings.sources.iter().enumerate() {
            checkpoint()?;
            let Some(source) = source else {
                continue;
            };
            let temperature = thermal.temperature_k(owner.emitters[index].binding.thermal_cell)?;
            let key = Some(temperature.to_bits());
            if let Some(prepared) = &self.prepared
                && prepared.temperatures[index] == key
            {
                temperatures[index] = key;
                laws[index] = prepared.laws[index];
                continue;
            }
            let law = source.evaluate(temperature)?;
            let model = super::resolve_model(
                &owner.singlet,
                &owner.triplet,
                temperature,
                &mut cached_model,
            )?;
            validate_pair(law, model.excitation_energies_j())?;
            temperatures[index] = key;
            laws[index] = Some(law);
        }
        checkpoint()?;
        Ok(Arc::new(PreparedAnnihilation {
            temperatures,
            laws,
            limits: self.settings.limits,
        }))
    }

    pub(super) fn accept(&mut self, prepared: Arc<PreparedAnnihilation>) {
        self.prepared = Some(prepared);
    }
}

fn validate_pair(law: TripletAnnihilation, energies: [f64; 2]) -> Result<(), OpticalHeatingError> {
    law.advance(energies, [0.0; 2], 0.0)
        .map_err(|_| invalid("triplet pair cannot pay its authored singlet energy"))?;
    Ok(())
}
