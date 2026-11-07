//! Evaluated radioactive data assembled into the original physical model owners.
//!
//! ENDF MF8/MT457 reports average emitted energies per parent transformation,
//! summed over decay modes. These are not branch-conditional spectra. Every
//! branch therefore receives the same explicitly routed parent-average energy:
//! expected parent calorimetry is preserved, while correlations between a
//! selected daughter and emitted energy are not supplied. Q values are never
//! substituted for emitted energy or heat.
//!
//! Daughter identities follow ENDF-102 (2023), §8.4: RTYP describes an initial
//! transformation followed by immediate particle emissions, and RFS selects the
//! final isomer. Only modes with one represented residual nucleus are lowered.
//! Spontaneous fission needs separately authored product yields. Unrepresented
//! daughters are refused; no missing evaluation is assumed stable. Burnup may
//! explicitly retain emitted alpha particles and protons as supplied He4/H1
//! inventories. This does not infer stopping, heat, or a nuclear mass/charge ledger.
//!
//! This is an import-time assembly adapter, not a second inventory or solver.
//! Decay admission/evolution and irradiation retain their existing owners.

use std::{collections::BTreeMap, sync::Arc};

use fs_la::decay::{DecayLimits, MAX_DECAY_EDGES, MAX_DECAY_NODES};

use super::DecayEvaluation;
use crate::burnup::{BurnupData, BurnupNuclide, BurnupRate, BurnupReaction};
use crate::decay::{DecayBranch, DecayError, DecayModel, DecayNuclide};
use crate::irradiation::{NeutronGroups, NuclideId};
use crate::{Accumulator, EnergyChannels, FissionError, nonnegative, positive, sum};

/// Energy-resolved emissions driven by the original integrated decay events.
pub mod emissions;

mod products;
pub use products::DecayProductRetention;
use products::{ChargedMultiplicity, RetentionBinding, resolve_mode};

/// Explicit fractions of one category's emitted energy sent to each destination.
/// Fractions must be finite, nonnegative, and sum to one within rounding; they
/// are never renormalized. The routing is a caller-supplied deposition model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayEnergyFractions {
    pub deposited: f64,
    pub radiation: f64,
    pub escaped: f64,
}

impl DecayEnergyFractions {
    fn values(self) -> [f64; 3] {
        [self.deposited, self.radiation, self.escaped]
    }

    fn validate(self) -> Result<(), DecayAssemblyError> {
        let values = self.values();
        if !values
            .iter()
            .all(|&value| nonnegative(value) && value <= 1.0)
            || (sum(&values) - 1.0).abs() > 32.0 * f64::EPSILON
        {
            return Err(DecayAssemblyError::Invalid(
                "energy destination fractions must be nonnegative and sum to one",
            ));
        }
        Ok(())
    }
}

/// Explicit treatment of the three aggregate ENDF decay-heat energy categories.
/// Photon escape/absorption and particle stopping are not inferred from isotope
/// identity. The heavy-particle category includes neutrons and nuclear recoil;
/// a common category fraction cannot resolve their different transport laws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayEnergyRouting {
    pub light_particles: DecayEnergyFractions,
    pub electromagnetic: DecayEnergyFractions,
    pub heavy_particles: DecayEnergyFractions,
    /// Include the seventeenth ENDF average, when explicitly reported, in the
    /// escaped channel only. Refuses an unknown/absent average. When false, only
    /// the three non-neutrino aggregates are accounted; Q supplies no remainder.
    pub include_reported_neutrinos: bool,
}

/// One evaluated nuclide paired with its caller-authored energy routing.
#[derive(Clone, Debug)]
pub struct DecayModelEntry {
    pub evaluation: Arc<DecayEvaluation>,
    pub routing: DecayEnergyRouting,
}

/// Refusal before a partial physical model is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum DecayAssemblyError {
    Invalid(&'static str),
    Budget(&'static str),
    DuplicateNuclide(NuclideId),
    UnsupportedDecayMode {
        parent: NuclideId,
        rtyp: f64,
    },
    InvalidDaughter {
        parent: NuclideId,
        rtyp: f64,
    },
    MissingDaughter {
        parent: NuclideId,
        daughter: NuclideId,
    },
    /// Retaining a charged coproduct requires its explicitly supplied evaluation.
    MissingRetainedProduct(NuclideId),
    /// Ground-state He4/H1 must have an explicitly stable evaluated law.
    UnstableRetainedProduct(NuclideId),
    /// Component names are ENDF ELP, EEM, EHP and neutrino energy.
    UnknownEnergy {
        nuclide: NuclideId,
        component: &'static str,
    },
    Physics(FissionError),
    Decay(DecayError),
}

impl From<FissionError> for DecayAssemblyError {
    fn from(error: FissionError) -> Self {
        Self::Physics(error)
    }
}

impl From<DecayError> for DecayAssemblyError {
    fn from(error: DecayError) -> Self {
        Self::Decay(error)
    }
}

impl core::fmt::Display for DecayAssemblyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "evaluated decay assembly refused: {self:?}")
    }
}

impl std::error::Error for DecayAssemblyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Physics(error) => Some(error),
            Self::Decay(error) => Some(error),
            _ => None,
        }
    }
}

/// Immutable import result in the caller's original nuclide order.
/// The original model is shared directly with inventories and thermal sources;
/// evaluation uncertainties remain observable metadata and are not propagated
/// into numerical error budgets. Source labels are declarations, not validation.
#[derive(Clone, Debug)]
pub struct EvaluatedDecayAssembly {
    model: Arc<DecayModel>,
    evaluations: Vec<Arc<DecayEvaluation>>,
    routing: Vec<DecayEnergyRouting>,
    nuclides: Vec<NuclideId>,
    indices: BTreeMap<NuclideId, usize>,
    // For each admitted physical branch, retain its original record/mode index.
    mode_sources: Vec<(usize, usize)>,
    charged_multiplicities: Vec<ChargedMultiplicity>,
    product_retention: RetentionBinding,
}

impl EvaluatedDecayAssembly {
    /// Resolve a complete closed set of evaluated decay products. The original
    /// DecayModel admits positive branching, branch sums, topology and numerical
    /// limits. No probability normalization, stable endpoint or isotope is added.
    /// Zero-probability source entries create no physical branch.
    pub fn new(
        data_id: String,
        entries: Vec<DecayModelEntry>,
        limits: DecayLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DecayAssemblyError> {
        checkpoint()?;
        if data_id.trim().is_empty() || entries.is_empty() {
            return Err(DecayAssemblyError::Invalid(
                "an identified nonempty evaluated decay assembly is required",
            ));
        }
        if entries.len() > limits.maximum_nodes || entries.len() > MAX_DECAY_NODES {
            return Err(DecayAssemblyError::Budget("evaluated decay nuclides"));
        }
        let mut mode_count = 0_usize;
        let mut indices = BTreeMap::new();
        for (index, entry) in entries.iter().enumerate() {
            checkpoint()?;
            let id = entry.evaluation.metadata().nuclide;
            if indices.insert(id, index).is_some() {
                return Err(DecayAssemblyError::DuplicateNuclide(id));
            }
            mode_count = mode_count
                .checked_add(entry.evaluation.modes().len())
                .ok_or(DecayAssemblyError::Budget("evaluated decay modes"))?;
            if mode_count > limits.maximum_edges || mode_count > MAX_DECAY_EDGES {
                return Err(DecayAssemblyError::Budget("evaluated decay modes"));
            }
            entry.routing.light_particles.validate()?;
            entry.routing.electromagnetic.validate()?;
            entry.routing.heavy_particles.validate()?;
        }
        let mut nuclides = Vec::with_capacity(entries.len());
        let mut nodes = Vec::with_capacity(entries.len());
        let mut branches = Vec::with_capacity(mode_count);
        let mut mode_sources = Vec::with_capacity(mode_count);
        let mut charged_multiplicities = Vec::with_capacity(mode_count);
        for (parent, entry) in entries.iter().enumerate() {
            checkpoint()?;
            let evaluation = &entry.evaluation;
            let id = evaluation.metadata().nuclide;
            nuclides.push(id);
            nodes.push(DecayNuclide {
                id: format!(
                    "Z{}-A{}-m{}",
                    id.atomic_number(),
                    id.mass_number(),
                    id.isomer()
                ),
                law: evaluation.law(),
            });
            let energy = routed_energy(evaluation, entry.routing)?;
            for (mode_index, mode) in evaluation.modes().iter().enumerate() {
                checkpoint()?;
                if mode.probability == 0.0 {
                    continue;
                }
                let resolved = resolve_mode(id, mode.rtyp, mode.daughter_isomer)?;
                let daughter_id = resolved.residual;
                let &daughter =
                    indices
                        .get(&daughter_id)
                        .ok_or(DecayAssemblyError::MissingDaughter {
                            parent: id,
                            daughter: daughter_id,
                        })?;
                branches.push(DecayBranch {
                    parent,
                    daughter,
                    probability: mode.probability,
                    energy_per_decay: energy,
                });
                mode_sources.push((parent, mode_index));
                charged_multiplicities.push(resolved.charged);
            }
        }
        checkpoint()?;
        let model = Arc::new(DecayModel::new(data_id, nodes, branches, limits)?);
        let (evaluations, routing) = entries
            .into_iter()
            .map(|entry| (entry.evaluation, entry.routing))
            .unzip();
        checkpoint()?;
        Ok(Self {
            model,
            evaluations,
            routing,
            nuclides,
            indices,
            mode_sources,
            charged_multiplicities,
            product_retention: RetentionBinding::default(),
        })
    }

    /// The original admitted physical model, directly reusable by decay owners.
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
    }

    /// Explicit identities in inventory/population order.
    pub fn nuclides(&self) -> &[NuclideId] {
        &self.nuclides
    }

    /// Original immutable evaluations, including their uncertainties and Q values.
    pub fn evaluations(&self) -> &[Arc<DecayEvaluation>] {
        &self.evaluations
    }

    /// Admitted routing assumptions in the same original nuclide order.
    pub fn routing(&self) -> &[DecayEnergyRouting] {
        &self.routing
    }

    /// Resolve one identity without assuming a sorted input order.
    pub fn index(&self, nuclide: NuclideId) -> Option<usize> {
        self.indices.get(&nuclide).copied()
    }

    /// Produce ordinary burnup data with the same order, half-lives and decay
    /// channels and the explicitly selected charged-coproduct retention policy.
    /// Append independently supplied capture/fission reactions using
    /// index(), then admit through BurnupModel::new with its original work limits.
    /// No cross sections, fission yields or material-temperature range are inferred.
    pub fn burnup_data(
        &self,
        groups: Arc<NeutronGroups>,
        temperature_min_k: f64,
        temperature_max_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupData, DecayAssemblyError> {
        checkpoint()?;
        if !positive(temperature_min_k)
            || !positive(temperature_max_k)
            || temperature_max_k < temperature_min_k
        {
            return Err(DecayAssemblyError::Invalid(
                "an explicit finite positive burnup temperature range is required",
            ));
        }
        let mut nuclides = Vec::with_capacity(self.nuclides.len());
        for (&id, node) in self.nuclides.iter().zip(self.model.nuclides()) {
            checkpoint()?;
            nuclides.push(BurnupNuclide { id, law: node.law });
        }
        let mut reactions = Vec::with_capacity(self.mode_sources.len());
        for (index, (branch, &(record, mode))) in self
            .model
            .branches()
            .iter()
            .zip(&self.mode_sources)
            .enumerate()
        {
            checkpoint()?;
            let source = &self.evaluations[record];
            reactions.push(BurnupReaction {
                data_id: format!(
                    "{} | {} | MAT{} MF8/MT457 mode{}",
                    self.model.data_id(),
                    source.data_id(),
                    source.metadata().material_number,
                    mode
                ),
                parent: branch.parent,
                rate: BurnupRate::Decay {
                    branch_fraction: branch.probability,
                },
                products: self.burnup_products(index).as_slice().to_vec(),
                energy_per_event: branch.energy_per_decay,
            });
        }
        checkpoint()?;
        Ok(BurnupData {
            data_id: self.model.data_id().to_owned(),
            groups,
            temperature_min_k,
            temperature_max_k,
            nuclides,
            reactions,
        })
    }
}

fn routed_energy(
    evaluation: &DecayEvaluation,
    routing: DecayEnergyRouting,
) -> Result<EnergyChannels, DecayAssemblyError> {
    if evaluation.law().decay_constant().value() == 0.0 {
        return Ok(EnergyChannels::default());
    }
    let mut totals = [Accumulator::new(0.0); 3];
    for (index, name, fractions) in [
        (0, "ELP", routing.light_particles),
        (1, "EEM", routing.electromagnetic),
        (2, "EHP", routing.heavy_particles),
    ] {
        let energy = average_energy(evaluation, index, name)?;
        for (total, fraction) in totals.iter_mut().zip(fractions.values()) {
            let value = energy * fraction;
            if !nonnegative(value) || (energy > 0.0 && fraction > 0.0 && value == 0.0) {
                return Err(FissionError::Arithmetic.into());
            }
            total.add(value)?;
        }
    }
    if routing.include_reported_neutrinos {
        totals[2].add(average_energy(evaluation, 16, "neutrino energy")?)?;
    }
    let channels = EnergyChannels {
        deposited_j: totals[0].value(),
        radiation_j: totals[1].value(),
        escaped_j: totals[2].value(),
    };
    channels.validate()?;
    Ok(channels)
}

fn average_energy(
    evaluation: &DecayEvaluation,
    index: usize,
    component: &'static str,
) -> Result<f64, DecayAssemblyError> {
    evaluation
        .average_energies()
        .get(index)
        .and_then(|energy| energy.value_j)
        .ok_or(DecayAssemblyError::UnknownEnergy {
            nuclide: evaluation.metadata().nuclide,
            component,
        })
}

#[cfg(test)]
mod tests;
