//! Neutron kinetic and prompt-reaction energy on the multigroup reaction operator.
//!
//! The original diffusion owner advances neutron populations, continuous births,
//! reaction counts and face histories. This companion attaches an explicitly
//! authored representative energy to each group. Downscattering removes the
//! group-energy difference; capture removes the incident kinetic energy.
//! Prompt fission and other neutron production retain their distinct authored
//! multiplicities and event counters. Explicit signed prompt Q and threshold
//! declarations give the energy available
//! to heat/radiation as incident energy + Q - outgoing neutron energy. Nuclear
//! conversion is retained separately, including endothermic consumption; leakage
//! retains its outgoing neutron energy. [`activation`] binds finite target and
//! product inventories to the same accepted history. No Q, threshold, secondary
//! spectrum, stopping law or photon transport is inferred.
//! Fission Q includes only prompt conversion; radioactive product decays retain
//! their own event energies. Explicit precursor storage reserves each original
//! family's future neutron kinetic energy. Its native decay transfers that
//! energy into neutrons; no additional beta/gamma release is inferred.
//! Explicit induced-fission funding reserves separately declared future-neutron
//! nuclear conversion at each family birth. Prompt Q retains its original
//! meaning; delayed funding cannot pay a prompt energy deficit or threshold.
//!
//! A fusion interval retains separate D-T and D-D-neutron counts and energies.
//! `FusionNeutronSource::project_interval_mean` explicitly replaces each channel
//! by its mean energy, projects conservatively to neighboring representative
//! energies, and assumes uniform emission during that interval. This is a
//! numerical source closure, not a resolved thermal/Doppler spectrum or birth
//! history. Split/refine fusion intervals when birth timing affects transport.
//!
//! Diffusion remains an angular approximation, inappropriate for resolved void
//! streaming or strongly anisotropic sources. Its source equations and limits:
//! <https://arfc.github.io/moltres/getting_started/theory.html>. Distinct capture,
//! scattering and outgoing-particle energy follow the physical channels in
//! <https://docs.openmc.org/en/stable/methods/neutron_physics.html>.

use std::sync::Arc;

/// Accepted fusion intervals projected onto one original spatial neutron model.
pub mod history;

use fs_thermochem::fusion::neutrons::FusionNeutronEmissionV1;

use crate::diffusion::{
    DiffusionBudget, DiffusionError, DiffusionModel, DiffusionReactionKind, DiffusionState,
    DiffusionStep,
};
use crate::{Accumulator, EnergyChannels, FissionError, nonnegative, positive, sum};

/// Finite target activation and reaction energy on the same neutron history.
pub mod activation;

mod precursors;

/// Authored prompt nuclear mass-energy conversion for one neutron reaction.
///
/// Positive Q releases rest energy; negative Q consumes incident kinetic energy.
/// Q is the prompt nuclear conversion before subtracting outgoing neutron
/// kinetic energy. Incident kinetic energy and later product decays are separate
/// accounts. The laboratory threshold is an independent physical input:
/// `-Q` alone does not include target recoil or determine a reaction threshold.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NeutronReactionEnergy {
    pub nuclear_q_j: f64,
    pub threshold_j: f64,
}

/// Explicit funding for the future neutron energy of induced-fission families.
///
/// This declaration supplements, and never changes, the original prompt Q.
/// The actual reserve is the original delayed yield times each family's neutron
/// energy. Unused allowance releases no energy. A finite inventory consumer
/// must additionally admit ownership of the unrepresented fragment cohort.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InducedFissionPrecursorFunding {
    /// Separate nuclear conversion stored for later neutron emission [J/fission].
    /// No delayed beta/gamma energy or local heat is inferred from this cap.
    AdditionalDelayedConversion { maximum_energy_j_per_fission: f64 },
}

/// Immutable energy interpretation of an original neutron-diffusion model.
#[derive(Clone, Debug, PartialEq)]
pub struct NeutronKineticModel {
    data_id: String,
    diffusion: Arc<DiffusionModel>,
    energies_j: Vec<f64>,
    precursor_energies_j: Vec<f64>,
    fission_precursor_funding: Vec<Option<InducedFissionPrecursorFunding>>,
    fission_precursor_energy_j: Vec<f64>,
    induced_fission_reactions: Vec<usize>,
    local_heat_fractions: Vec<f64>,
    reaction_energies: Vec<NeutronReactionEnergy>,
    available_energy_j: Vec<f64>,
}

impl NeutronKineticModel {
    /// Bind one positive representative energy per original ascending group,
    /// and one local heat fraction per original reaction. The remaining removed
    /// kinetic energy is outgoing radiation, never silently local heat.
    ///
    /// Group energies must lie inside their original group edges. Source-owned
    /// effective speeds and cross sections remain unchanged; their averaging
    /// and physical validity must be compatible with this energy discretization.
    /// Only capture and nonincreasing-energy scattering are admitted. Explicit
    /// reaction fractions are required even for an exactly elastic self-scatter.
    pub fn new(
        data_id: String,
        diffusion: Arc<DiffusionModel>,
        energies_j: Vec<f64>,
        local_heat_fractions: Vec<f64>,
    ) -> Result<Self, NeutronKineticError> {
        for reaction in diffusion.reactions() {
            match reaction.kind {
                DiffusionReactionKind::Capture | DiffusionReactionKind::Scatter { .. } => {}
                _ => {
                    return Err(NeutronKineticError::Invalid(
                        "neutron production requires explicit reaction energies and thresholds",
                    ));
                }
            }
        }
        let reaction_energies = vec![NeutronReactionEnergy::default(); diffusion.reactions().len()];
        let model = Self::new_with_reaction_energies(
            data_id,
            diffusion,
            energies_j,
            local_heat_fractions,
            reaction_energies,
        )?;
        if model.diffusion.reactions().iter().any(|reaction| {
            matches!(reaction.kind, DiffusionReactionKind::Scatter { outgoing_group }
                if model.energies_j[outgoing_group] > model.energies_j[reaction.incident_group])
        }) {
            return Err(NeutronKineticError::Invalid(
                "legacy kinetic-only scattering must not increase group energy",
            ));
        }
        Ok(model)
    }

    /// Admit explicit prompt fission, other neutron production and nuclear energy.
    ///
    /// Every original reaction supplies signed Q, a nonnegative laboratory
    /// threshold, and an independent local heat fraction. Secondary yields are
    /// the original diffusion declaration and are never normalized. Each active
    /// reaction must lie above threshold and satisfy E_in + Q >= sum(yield E_out).
    /// Zero cross sections may cover below-threshold groups; composition- or
    /// temperature-driven rebinding must still pass these checks before use.
    /// These are representative-energy checks, not within-group kinematics;
    /// source data must be consistently collapsed and thresholds resolved.
    /// Fission Q excludes later radioactive-product release. Only prompt fission
    /// yields are admitted: abstract delayed-precursor energy storage is not
    /// inferred. Explicit radioactive products use the original activation owner.
    pub fn new_with_reaction_energies(
        data_id: String,
        diffusion: Arc<DiffusionModel>,
        energies_j: Vec<f64>,
        local_heat_fractions: Vec<f64>,
        reaction_energies: Vec<NeutronReactionEnergy>,
    ) -> Result<Self, NeutronKineticError> {
        Self::prepare(
            data_id,
            diffusion,
            energies_j,
            local_heat_fractions,
            reaction_energies,
            false,
            None,
        )
    }

    /// Admit original delayed families with explicitly reserved neutron energy.
    ///
    /// Each precursor stores `sum(spectrum[g] * group_energy[g])` joules until
    /// its original decay emits that spectrum. Initial precursor inventories
    /// and external precursor feeds therefore supply this stored energy; bound
    /// radioactive births must debit it from their original event energy.
    /// Decay emits no additional heat or radiation here. Material outflow
    /// exports the reserve separately from neutron leakage and local heat.
    ///
    /// Positive induced-fission delayed yields remain unsupported: prompt Q
    /// does not declare funding for later precursor release. Matching all-zero
    /// delayed yields are admitted without inventing such a contribution.
    pub fn new_with_precursor_storage(
        data_id: String,
        diffusion: Arc<DiffusionModel>,
        energies_j: Vec<f64>,
        local_heat_fractions: Vec<f64>,
        reaction_energies: Vec<NeutronReactionEnergy>,
    ) -> Result<Self, NeutronKineticError> {
        Self::prepare(
            data_id,
            diffusion,
            energies_j,
            local_heat_fractions,
            reaction_energies,
            true,
            None,
        )
    }

    /// Admit explicit additional future-neutron conversion for induced fission.
    ///
    /// Funding follows the exact original reaction order. Positive delayed
    /// yields require a declaration even while their cross section is zero;
    /// nonfission reactions require None. Prompt Q, thresholds and prompt
    /// outgoing-energy admission remain unchanged. The native family reserve
    /// must fit its per-fission cap; unused cap is never generated or deposited.
    /// A finite activation consumer must separately establish exclusive
    /// ownership of each induced-fission fragment cohort.
    pub fn new_with_fission_precursor_funding(
        data_id: String,
        diffusion: Arc<DiffusionModel>,
        energies_j: Vec<f64>,
        local_heat_fractions: Vec<f64>,
        reaction_energies: Vec<NeutronReactionEnergy>,
        funding: Vec<Option<InducedFissionPrecursorFunding>>,
    ) -> Result<Self, NeutronKineticError> {
        Self::prepare(
            data_id,
            diffusion,
            energies_j,
            local_heat_fractions,
            reaction_energies,
            true,
            Some(funding),
        )
    }

    fn prepare(
        data_id: String,
        diffusion: Arc<DiffusionModel>,
        energies_j: Vec<f64>,
        local_heat_fractions: Vec<f64>,
        reaction_energies: Vec<NeutronReactionEnergy>,
        precursor_storage: bool,
        funding: Option<Vec<Option<InducedFissionPrecursorFunding>>>,
    ) -> Result<Self, NeutronKineticError> {
        if data_id.trim().is_empty()
            || energies_j.len() != diffusion.groups().len()
            || local_heat_fractions.len() != diffusion.reactions().len()
            || reaction_energies.len() != diffusion.reactions().len()
            || funding
                .as_ref()
                .is_some_and(|values| values.len() != diffusion.reactions().len())
            || (!precursor_storage && !diffusion.delayed_groups().is_empty())
        {
            return Err(NeutronKineticError::Invalid(
                "identified original groups, matching reaction energies/fractions and explicit precursor storage required",
            ));
        }
        for (group, &energy) in energies_j.iter().enumerate() {
            let edges = diffusion.groups().edges_j();
            if !positive(energy)
                || energy < edges[group]
                || energy > edges[group + 1]
                || (group > 0 && energy <= energies_j[group - 1])
            {
                return Err(NeutronKineticError::Invalid(
                    "ascending positive representative energies inside original groups required",
                ));
            }
        }
        let precursor_energies_j = precursors::energies(&diffusion, &energies_j)?;
        let fission_precursor_funding =
            funding.unwrap_or_else(|| vec![None; diffusion.reactions().len()]);
        let fission_precursor_energy_j = precursors::funding(
            &diffusion,
            &precursor_energies_j,
            &fission_precursor_funding,
        )?;
        // The hot accounting path visits only reactions that can create a
        // funded family; ordinary capture/scatter models pay no repeated scan.
        let induced_fission_reactions = fission_precursor_energy_j
            .iter()
            .enumerate()
            .filter_map(|(index, &energy)| (energy > 0.0).then_some(index))
            .collect();
        let mut available_energy_j = Vec::with_capacity(diffusion.reactions().len());
        for ((reaction, &fraction), energy) in diffusion
            .reactions()
            .iter()
            .zip(&local_heat_fractions)
            .zip(&reaction_energies)
        {
            if !fraction.is_finite()
                || !(0.0..=1.0).contains(&fraction)
                || !energy.nuclear_q_j.is_finite()
                || !nonnegative(energy.threshold_j)
                || energy.threshold_j < (-energy.nuclear_q_j).max(0.0)
            {
                return Err(NeutronKineticError::Invalid(
                    "finite signed Q, threshold at least max(0,-Q), and heat fraction in [0,1] required",
                ));
            }
            let outgoing = match &reaction.kind {
                DiffusionReactionKind::Capture => 0.0,
                DiffusionReactionKind::Scatter { outgoing_group } => energies_j[*outgoing_group],
                DiffusionReactionKind::NeutronProduction { neutron_yields } => {
                    weighted_sum(neutron_yields, &energies_j)?
                }
                DiffusionReactionKind::Fission {
                    prompt_yields,
                    delayed_yields,
                } => {
                    if !precursor_storage && !delayed_yields.is_empty() {
                        return Err(NeutronKineticError::Invalid(
                            "induced-fission delayed production requires an explicit physical energy funding declaration",
                        ));
                    }
                    weighted_sum(prompt_yields, &energies_j)?
                }
            };
            available_energy_j.push(checked(sum(&[
                energies_j[reaction.incident_group],
                energy.nuclear_q_j,
                -outgoing,
            ]))?);
        }
        let model = Self {
            data_id,
            diffusion,
            energies_j,
            precursor_energies_j,
            fission_precursor_funding,
            fission_precursor_energy_j,
            induced_fission_reactions,
            local_heat_fractions,
            reaction_energies,
            available_energy_j,
        };
        model.admit_active_reactions(&model.diffusion)?;
        Ok(model)
    }

    /// Original declaration of group energies and reaction energy partition.
    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    /// The unchanged source-bound neutron operator.
    pub fn diffusion(&self) -> &Arc<DiffusionModel> {
        &self.diffusion
    }
    /// Authored representative neutron kinetic energies [J/particle].
    pub fn group_energies_j(&self) -> &[f64] {
        &self.energies_j
    }
    /// Reserved future neutron energy in original delayed-family order [J].
    pub fn precursor_energies_j(&self) -> &[f64] {
        &self.precursor_energies_j
    }
    /// Explicit future-neutron funding in unchanged reaction order.
    pub fn fission_precursor_funding(&self) -> &[Option<InducedFissionPrecursorFunding>] {
        &self.fission_precursor_funding
    }
    /// Actual reserved future-neutron conversion per induced event [J].
    /// Zero for nonfission reactions and zero-yield families; unused cap is absent.
    pub fn fission_precursor_energy_j(&self) -> &[f64] {
        &self.fission_precursor_energy_j
    }
    /// Local heat fractions in original reaction order.
    pub fn local_heat_fractions(&self) -> &[f64] {
        &self.local_heat_fractions
    }

    /// Signed prompt Q and threshold in the unchanged reaction ordering.
    pub fn reaction_energies(&self) -> &[NeutronReactionEnergy] {
        &self.reaction_energies
    }

    fn admit_active_reactions(&self, current: &DiffusionModel) -> Result<(), NeutronKineticError> {
        if current.reactions().len() != self.reaction_energies.len() {
            return Err(NeutronKineticError::Invalid(
                "changed reaction energy layout",
            ));
        }
        for (index, reaction) in current.reactions().iter().enumerate() {
            if reaction.macroscopic_cross_section_m_inv > 0.0
                && (self.energies_j[reaction.incident_group]
                    < self.reaction_energies[index].threshold_j
                    || !nonnegative(self.available_energy_j[index]))
            {
                return Err(NeutronKineticError::Invalid(
                    "active neutron reaction violates its threshold or outgoing energy budget",
                ));
            }
        }
        Ok(())
    }
}

/// Explicit interval-uniform, channel-mean source projection to one spatial cell.
#[derive(Clone, Debug, PartialEq)]
pub struct FusionNeutronSource {
    model: Arc<NeutronKineticModel>,
    emission: FusionNeutronEmissionV1,
    cell: usize,
    rates_s: Vec<f64>,
    number_projection_residual: f64,
    energy_projection_residual_j: f64,
}

impl FusionNeutronSource {
    /// Preserve total particle count and kinetic energy by linear interpolation
    /// between neighboring representative energies. No extrapolation, nearest
    /// group clipping, extra branch fraction or source renormalization occurs.
    ///
    /// This intentionally named approximation collapses each original channel
    /// to its mean and its interval to a constant rate. The immutable original
    /// emission remains attached. A zero source may have zero duration.
    pub fn project_interval_mean(
        model: Arc<NeutronKineticModel>,
        emission: FusionNeutronEmissionV1,
        cell: usize,
    ) -> Result<Self, NeutronKineticError> {
        if cell >= model.diffusion.cell_count() {
            return Err(NeutronKineticError::Invalid(
                "fusion source cell outside original geometry",
            ));
        }
        let energies = &model.energies_j;
        let groups = energies.len();
        let mut counts = vec![Accumulator::new(0.0); model.diffusion.cell_count() * groups];
        for channel in emission.channels() {
            let Some(mean) = channel.mean_particle_energy() else {
                continue;
            };
            let energy = mean.value();
            let count = channel.particle_count();
            if energy < energies[0] || energy > energies[groups - 1] {
                return Err(NeutronKineticError::Invalid(
                    "fusion channel mean outside representative neutron energies",
                ));
            }
            let upper = energies.partition_point(|&value| value < energy);
            if energies[upper] == energy {
                counts[cell * groups + upper].add(count)?;
            } else {
                let lower = upper - 1;
                let span = energies[upper] - energies[lower];
                let low = checked(count * ((energies[upper] - energy) / span))?;
                let high = checked(count * ((energy - energies[lower]) / span))?;
                // Underflowed positive weights cannot silently discard a source.
                if low <= 0.0 || high <= 0.0 {
                    return Err(NeutronKineticError::Arithmetic(
                        "fusion source projection underflow",
                    ));
                }
                counts[cell * groups + lower].add(low)?;
                counts[cell * groups + upper].add(high)?;
            }
        }
        let duration = emission.duration().value();
        let mut number = Accumulator::new(0.0);
        let mut energy = Accumulator::new(0.0);
        let mut rates_s = Vec::with_capacity(counts.len());
        for (node, count) in counts.iter().enumerate() {
            let count = count.value();
            let rate = if count == 0.0 {
                0.0
            } else {
                checked(count / duration)?
            };
            if count > 0.0 && rate == 0.0 {
                return Err(NeutronKineticError::Arithmetic(
                    "fusion source rate underflow",
                ));
            }
            let represented = checked(rate * duration)?;
            number.add(represented)?;
            energy.add(checked(represented * energies[node % groups])?)?;
            rates_s.push(rate);
        }
        let number_projection_residual = checked(number.value() - emission.total_particle_count())?;
        let energy_projection_residual_j =
            checked(energy.value() - emission.total_kinetic_energy().value())?;
        // A projection is accepted only at floating-point accounting accuracy;
        // the caller separately chooses transport/time/model error budgets.
        if number_projection_residual.abs() > 256.0 * f64::EPSILON * emission.total_particle_count()
            || energy_projection_residual_j.abs()
                > 256.0 * f64::EPSILON * emission.total_kinetic_energy().value()
        {
            return Err(NeutronKineticError::Balance(
                "fusion projection does not preserve count and energy",
            ));
        }
        Ok(Self {
            model,
            emission,
            cell,
            rates_s,
            number_projection_residual,
            energy_projection_residual_j,
        })
    }

    /// Exact retained original D-T and D-D-neutron interval accounts.
    pub const fn emission(&self) -> &FusionNeutronEmissionV1 {
        &self.emission
    }
    /// Immutable transport/energy interpretation used by this projection.
    pub fn model(&self) -> &Arc<NeutronKineticModel> {
        &self.model
    }
    /// Original spatial injection cell.
    pub const fn cell_index(&self) -> usize {
        self.cell
    }
    /// Cell-major constant rates [particles/s]; validity is the authored interval.
    pub fn neutron_rates_s(&self) -> &[f64] {
        &self.rates_s
    }
    /// Projected count minus original emitted count [particles].
    pub const fn number_projection_residual(&self) -> f64 {
        self.number_projection_residual
    }
    /// Projected kinetic energy minus original emitted energy [J].
    pub const fn energy_projection_residual_j(&self) -> f64 {
        self.energy_projection_residual_j
    }
}

/// Original diffusion work limits plus an absolute kinetic-energy account limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeutronKineticBudget {
    /// Original count/propagation/work limits, unmodified by this adapter.
    pub diffusion: DiffusionBudget,
    /// Noncancelling energy-account allowance [J], not physical uncertainty.
    pub energy_tolerance_j: f64,
}

/// Neutron storage and cumulative energy exports on one accepted history.
#[derive(Clone, Debug)]
pub struct NeutronKineticState {
    model: Arc<NeutronKineticModel>,
    neutrons: DiffusionState,
    released: Vec<[Accumulator; 3]>,
    nuclear_reaction_energy: Vec<Accumulator>,
    induced_precursor_energy: Vec<Accumulator>,
    exported_precursor_energy: Vec<Accumulator>,
}

/// The same interval's particle, energy, and per-material-cell heat histories.
#[derive(Clone, Debug, PartialEq)]
pub struct NeutronKineticStep {
    /// Accepted native reaction/face histories; no neutron solver is duplicated.
    pub diffusion: DiffusionStep,
    /// Local heat, outgoing reaction radiation, and neutron leakage by cell [J].
    pub cell_released: Vec<EnergyChannels>,
    /// Aggregated disjoint exports [J]. Only deposited_j is local heat.
    pub released: EnergyChannels,
    /// Signed prompt nuclear conversion by cell [J]. Negative values are energy
    /// consumed into nuclear rest mass, not negative deposited heat. This Q is
    /// already included in `cell_released` after outgoing neutron kinetic energy.
    pub cell_nuclear_reaction_energy_j: Vec<f64>,
    /// Sum of signed prompt nuclear conversion on this neutron history [J].
    pub nuclear_reaction_energy_j: f64,
    /// Additional nuclear conversion reserved by induced-fission family births,
    /// in original spatial-cell order [J]. Neither prompt heat nor injection.
    pub cell_induced_precursor_energy_j: Vec<f64>,
    /// Sum of additional future-neutron conversion reserved at induced births [J].
    pub induced_precursor_energy_j: f64,
    /// Initial stored neutron kinetic energy [J].
    pub initial_kinetic_energy_j: f64,
    /// Final stored neutron kinetic energy [J].
    pub final_kinetic_energy_j: f64,
    /// Actual integrated energy brought by external neutron injection [J].
    pub injected_kinetic_energy_j: f64,
    /// Initial reserved future-neutron energy in original precursor stock [J].
    pub initial_precursor_energy_j: f64,
    /// Final reserved future-neutron energy in original precursor stock [J].
    pub final_precursor_energy_j: f64,
    /// Energy reserved by actual integrated precursor births [J]. Native
    /// internal radioactive births are included; activation reconciles them.
    pub injected_precursor_energy_j: f64,
    /// Internal reserve transferred into delayed neutron kinetic energy [J].
    /// This is neither an external source nor an additional energy release.
    pub precursor_conversion_energy_j: f64,
    /// Reserved energy carried out by original precursor material flows [J].
    /// Separate from `cell_released`; no future exterior decay is credited.
    pub cell_exported_precursor_energy_j: Vec<f64>,
    /// Sum of reserved precursor energy exported to the exterior [J].
    pub exported_precursor_energy_j: f64,
    /// Final - initial neutron and precursor storage + disjoint exports
    /// - injected neutron and precursor energy - prompt nuclear Q
    /// - additional induced-fission precursor conversion [J].
    pub energy_balance_residual_j: f64,
    /// Sum of absolute per-cell energy/account defects [J].
    pub energy_roundoff_j: f64,
    /// Original fusion intervals retained when trial_with_fusion was used.
    pub fusion_sources: Vec<FusionNeutronSource>,
    /// Native injection minus original fusion emission, separately from projection.
    pub fusion_number_residual: f64,
    /// Native injected energy minus original fusion emitted energy [J].
    pub fusion_energy_residual_j: f64,
}

impl NeutronKineticState {
    /// Bind initial extensive neutron counts to the unchanged diffusion model.
    pub fn new(
        model: Arc<NeutronKineticModel>,
        neutron_counts: Vec<f64>,
    ) -> Result<Self, NeutronKineticError> {
        let precursors = vec![0.0; model.diffusion.cell_count() * model.precursor_energies_j.len()];
        Self::new_with_precursors(model, neutron_counts, precursors)
    }

    /// Bind original cell-major neutron and precursor counts. Initial precursor
    /// stock supplies its declared future-neutron energy as initial storage.
    pub fn new_with_precursors(
        model: Arc<NeutronKineticModel>,
        neutron_counts: Vec<f64>,
        precursor_counts: Vec<f64>,
    ) -> Result<Self, NeutronKineticError> {
        let neutrons = DiffusionState::new(
            Arc::clone(&model.diffusion),
            neutron_counts,
            precursor_counts,
        )?;
        let released = vec![[Accumulator::new(0.0); 3]; model.diffusion.cell_count()];
        let nuclear_reaction_energy = vec![Accumulator::new(0.0); model.diffusion.cell_count()];
        let induced_precursor_energy = vec![Accumulator::new(0.0); model.diffusion.cell_count()];
        let exported_precursor_energy = vec![Accumulator::new(0.0); model.diffusion.cell_count()];
        let state = Self {
            model,
            neutrons,
            released,
            nuclear_reaction_energy,
            induced_precursor_energy,
            exported_precursor_energy,
        };
        state.stored_energy_j()?;
        Ok(state)
    }

    /// Retained energy/source interpretation.
    pub fn model(&self) -> &Arc<NeutronKineticModel> {
        &self.model
    }
    /// Actual underlying neutron inventory and original group fluxes.
    pub const fn neutrons(&self) -> &DiffusionState {
        &self.neutrons
    }
    /// Total retained neutron kinetic energy [J].
    pub fn kinetic_energy_j(&self) -> Result<f64, NeutronKineticError> {
        weighted_sum(self.neutrons.neutron_counts(), &self.model.energies_j)
    }
    /// Total reserved future-neutron energy in the retained precursor stock [J].
    pub fn precursor_energy_j(&self) -> Result<f64, NeutronKineticError> {
        weighted_sum(
            self.neutrons.precursor_counts(),
            &self.model.precursor_energies_j,
        )
    }
    /// Neutron kinetic energy plus reserved precursor energy [J].
    pub fn stored_energy_j(&self) -> Result<f64, NeutronKineticError> {
        checked(sum(&[self.kinetic_energy_j()?, self.precursor_energy_j()?]))
    }
    /// Reserved energy exported by precursor outflow from this cell since
    /// binding [J]. It is not part of the cell's released heat or radiation.
    pub fn cell_exported_precursor_energy_j(&self, cell: usize) -> Option<f64> {
        self.exported_precursor_energy
            .get(cell)
            .copied()
            .map(Accumulator::value)
    }
    /// Lifetime reserved energy exported by all precursor outflows [J].
    pub fn cumulative_exported_precursor_energy_j(&self) -> Result<f64, NeutronKineticError> {
        let mut total = Accumulator::new(0.0);
        for value in &self.exported_precursor_energy {
            total.add(value.value())?;
        }
        Ok(total.value())
    }
    /// Accepted lifetime exports for each cell, suitable for the thermal owner.
    pub fn cumulative_cell_released(&self) -> Vec<EnergyChannels> {
        self.released
            .iter()
            .copied()
            .map(crate::channel_values)
            .collect()
    }
    /// Allocation-free lifetime exports for one original spatial cell.
    pub fn cell_released(&self, cell: usize) -> Option<EnergyChannels> {
        self.released.get(cell).copied().map(crate::channel_values)
    }

    /// Signed prompt nuclear conversion in one cell since this owner was bound.
    /// This separate account must not be added to the already combined heat.
    pub fn cell_nuclear_reaction_energy_j(&self, cell: usize) -> Option<f64> {
        self.nuclear_reaction_energy
            .get(cell)
            .copied()
            .map(Accumulator::value)
    }

    /// Additional future-neutron conversion reserved in this cell since binding.
    /// It is separate from prompt Q and from already released material heat.
    pub fn cell_induced_precursor_energy_j(&self, cell: usize) -> Option<f64> {
        self.induced_precursor_energy
            .get(cell)
            .copied()
            .map(Accumulator::value)
    }

    /// Lifetime additional conversion reserved by induced-fission families [J].
    pub fn cumulative_induced_precursor_energy_j(&self) -> Result<f64, NeutronKineticError> {
        let mut total = Accumulator::new(0.0);
        for value in &self.induced_precursor_energy {
            total.add(value.value())?;
        }
        Ok(total.value())
    }

    /// Native continuous source integration with explicit constant group rates.
    /// This entrypoint is also the seam for an independently scheduled source.
    /// Zero rates continue slowing/absorption/leakage after fusion switches off.
    pub fn trial_with_rates(
        &self,
        duration_s: f64,
        neutron_rates_s: &[f64],
        budget: NeutronKineticBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, NeutronKineticStep), NeutronKineticError> {
        let precursor_rates = vec![0.0; self.neutrons.precursor_counts().len()];
        self.trial_with_sources(
            duration_s,
            neutron_rates_s,
            &precursor_rates,
            budget,
            checkpoint,
        )
    }

    /// Native uniform neutron and precursor feeds. Precursor injection supplies
    /// its explicitly declared reserve; conversion follows each family's native
    /// within-interval residence time, transport and original decay constant.
    pub fn trial_with_sources(
        &self,
        duration_s: f64,
        neutron_rates_s: &[f64],
        precursor_rates_s: &[f64],
        budget: NeutronKineticBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, NeutronKineticStep), NeutronKineticError> {
        checkpoint()?;
        if !nonnegative(budget.energy_tolerance_j) {
            return Err(NeutronKineticError::Invalid(
                "finite nonnegative energy allowance required",
            ));
        }
        let (neutrons, diffusion) = self.neutrons.trial_with_sources(
            duration_s,
            neutron_rates_s,
            precursor_rates_s,
            budget.diffusion,
            &mut checkpoint,
        )?;
        self.account_accepted_diffusion(neutrons, diffusion, budget.energy_tolerance_j, checkpoint)
    }

    // Attach energy to an already accepted native history. The activation owner
    // supplies its coupled neutron/burnup step here; no neutron propagator is
    // run again and no fixed macroscopic coefficient replaces its finite stock.
    pub(crate) fn account_accepted_diffusion(
        &self,
        neutrons: DiffusionState,
        diffusion: DiffusionStep,
        energy_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, NeutronKineticStep), NeutronKineticError> {
        self.model.admit_active_reactions(neutrons.model())?;
        let groups = self.model.energies_j.len();
        let cells = self.model.diffusion.cell_count();
        let mut cell_exports = vec![[Accumulator::new(0.0); 3]; cells];
        let mut cell_balance = vec![Accumulator::new(0.0); cells];
        let mut defect = Accumulator::new(0.0);
        let mut injected = Accumulator::new(0.0);
        let mut initial_energy = Accumulator::new(0.0);
        let mut final_energy = Accumulator::new(0.0);
        let mut cell_nuclear = vec![Accumulator::new(0.0); cells];
        let mut nuclear_total = Accumulator::new(0.0);
        for (node, &energy) in self
            .model
            .energies_j
            .iter()
            .cycle()
            .take(self.neutrons.neutron_counts().len())
            .enumerate()
        {
            checkpoint()?;
            let before = checked(self.neutrons.neutron_counts()[node] * energy)?;
            let after = checked(neutrons.neutron_counts()[node] * energy)?;
            let input = checked(diffusion.injected_neutron_counts[node] * energy)?;
            initial_energy.add(before)?;
            final_energy.add(after)?;
            injected.add(input)?;
            for value in [after, -before, -input] {
                defect.add(cell_balance[node / groups].add(value)?.abs())?;
            }
        }
        let precursors = precursors::account(
            &self.model,
            &self.neutrons,
            &neutrons,
            &diffusion,
            &mut cell_balance,
            &mut defect,
            &mut checkpoint,
        )?;
        for (index, reaction) in self.model.diffusion.reactions().iter().enumerate() {
            checkpoint()?;
            let count = diffusion.reaction_counts[index];
            let available = self.model.available_energy_j[index];
            let q = self.model.reaction_energies[index].nuclear_q_j;
            let nuclear = checked(count * q)?;
            let energy = checked(count * available)?;
            let heat = checked(energy * self.model.local_heat_fractions[index])?;
            let radiation = checked(energy - heat)?;
            defect.add(count.mul_add(q, -nuclear).abs())?;
            defect.add(count.mul_add(available, -energy).abs())?;
            defect.add(cell_nuclear[reaction.cell].add(nuclear)?.abs())?;
            defect.add(nuclear_total.add(nuclear)?.abs())?;
            defect.add(cell_exports[reaction.cell][0].add(heat)?.abs())?;
            defect.add(cell_exports[reaction.cell][1].add(radiation)?.abs())?;
            defect.add(cell_balance[reaction.cell].add(energy)?.abs())?;
            defect.add(cell_balance[reaction.cell].add(-nuclear)?.abs())?;
        }
        for (face, boundary) in self.model.diffusion.leakage().iter().enumerate() {
            checkpoint()?;
            for (group, &energy) in self.model.energies_j.iter().enumerate() {
                let lost = checked(diffusion.leakage_counts[face * groups + group] * energy)?;
                defect.add(cell_exports[boundary.cell][2].add(lost)?.abs())?;
                defect.add(cell_balance[boundary.cell].add(lost)?.abs())?;
            }
        }
        // Reciprocal transfers cancel globally but each cell must also close.
        for (face, link) in self.model.diffusion.links().iter().enumerate() {
            checkpoint()?;
            for (group, &energy) in self.model.energies_j.iter().enumerate() {
                let ab = checked(diffusion.link_a_to_b_counts[face * groups + group] * energy)?;
                let ba = checked(diffusion.link_b_to_a_counts[face * groups + group] * energy)?;
                for (cell, values) in [(link.cell_a, [ab, -ba]), (link.cell_b, [ba, -ab])] {
                    for value in values {
                        defect.add(cell_balance[cell].add(value)?.abs())?;
                    }
                }
            }
        }
        let mut released = self.released.clone();
        let mut nuclear_reaction_energy = self.nuclear_reaction_energy.clone();
        let mut induced_precursor_energy = self.induced_precursor_energy.clone();
        let mut exported_precursor_energy = self.exported_precursor_energy.clone();
        let mut totals = [Accumulator::new(0.0); 3];
        let mut cell_released = Vec::with_capacity(cells);
        for cell in 0..cells {
            checkpoint()?;
            defect.add(cell_balance[cell].value().abs())?;
            defect.add(
                induced_precursor_energy[cell]
                    .add(precursors.cell_induced_j[cell])?
                    .abs(),
            )?;
            defect.add(
                exported_precursor_energy[cell]
                    .add(precursors.cell_exported_j[cell])?
                    .abs(),
            )?;
            defect.add(
                nuclear_reaction_energy[cell]
                    .add(cell_nuclear[cell].value())?
                    .abs(),
            )?;
            for channel in 0..3 {
                let energy = cell_exports[cell][channel].value();
                defect.add(released[cell][channel].add(energy)?.abs())?;
                defect.add(totals[channel].add(energy)?.abs())?;
            }
            cell_released.push(crate::channel_values(cell_exports[cell]));
        }
        let exports = crate::channel_values(totals);
        let energy_balance_residual_j = checked(sum(&[
            final_energy.value(),
            precursors.final_j,
            -initial_energy.value(),
            -precursors.initial_j,
            exports.deposited_j,
            exports.radiation_j,
            exports.escaped_j,
            precursors.exported_j,
            -injected.value(),
            -precursors.injected_j,
            -nuclear_total.value(),
            -precursors.induced_j,
        ]))?;
        defect.add(energy_balance_residual_j.abs())?;
        if defect.value() > energy_tolerance_j {
            return Err(NeutronKineticError::Balance(
                "neutron and precursor energy account exceeds allowance",
            ));
        }
        let step = NeutronKineticStep {
            diffusion,
            cell_released,
            released: exports,
            cell_nuclear_reaction_energy_j: cell_nuclear
                .into_iter()
                .map(Accumulator::value)
                .collect(),
            nuclear_reaction_energy_j: nuclear_total.value(),
            cell_induced_precursor_energy_j: precursors.cell_induced_j,
            induced_precursor_energy_j: precursors.induced_j,
            initial_kinetic_energy_j: initial_energy.value(),
            final_kinetic_energy_j: final_energy.value(),
            injected_kinetic_energy_j: injected.value(),
            initial_precursor_energy_j: precursors.initial_j,
            final_precursor_energy_j: precursors.final_j,
            injected_precursor_energy_j: precursors.injected_j,
            precursor_conversion_energy_j: precursors.converted_j,
            cell_exported_precursor_energy_j: precursors.cell_exported_j,
            exported_precursor_energy_j: precursors.exported_j,
            energy_balance_residual_j,
            energy_roundoff_j: defect.value(),
            fusion_sources: vec![],
            fusion_number_residual: 0.0,
            fusion_energy_residual_j: 0.0,
        };
        checkpoint()?;
        Ok((
            Self {
                model: Arc::clone(&self.model),
                neutrons,
                released,
                nuclear_reaction_energy,
                induced_precursor_energy,
                exported_precursor_energy,
            },
            step,
        ))
    }

    /// Evolve one or more original fusion emissions over their shared duration.
    /// Timing is explicitly interval-uniform; native continuous injection means
    /// newborn neutrons only react/transport for their own residence time.
    pub fn trial_with_fusion(
        &self,
        duration_s: f64,
        sources: &[FusionNeutronSource],
        budget: NeutronKineticBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, NeutronKineticStep), NeutronKineticError> {
        checkpoint()?;
        let mut rates = vec![Accumulator::new(0.0); self.neutrons.neutron_counts().len()];
        let mut source_count = Accumulator::new(0.0);
        let mut source_energy = Accumulator::new(0.0);
        for source in sources {
            checkpoint()?;
            if *source.model != *self.model || source.emission.duration().value() != duration_s {
                return Err(NeutronKineticError::Invalid(
                    "fusion source model and complete interval must match transport",
                ));
            }
            source_count.add(source.emission.total_particle_count())?;
            source_energy.add(source.emission.total_kinetic_energy().value())?;
            for (rate, &value) in rates.iter_mut().zip(&source.rates_s) {
                rate.add(value)?;
            }
        }
        let rates: Vec<_> = rates.into_iter().map(Accumulator::value).collect();
        let (next, mut step) =
            self.trial_with_rates(duration_s, &rates, budget, &mut checkpoint)?;
        step.fusion_number_residual =
            checked(step.diffusion.injected_neutrons - source_count.value())?;
        step.fusion_energy_residual_j =
            checked(step.injected_kinetic_energy_j - source_energy.value())?;
        step.energy_roundoff_j =
            checked(step.energy_roundoff_j + step.fusion_energy_residual_j.abs())?;
        if step.fusion_number_residual.abs() + step.diffusion.population_roundoff
            > budget.diffusion.population_tolerance
            || step.energy_roundoff_j > budget.energy_tolerance_j
        {
            return Err(NeutronKineticError::Balance(
                "fusion-to-neutron source account exceeds allowance",
            ));
        }
        step.fusion_sources = sources.to_vec();
        checkpoint()?;
        Ok((next, step))
    }
}

/// Atomic source/domain/work/energy refusal; input histories remain unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum NeutronKineticError {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Balance(&'static str),
    Diffusion(DiffusionError),
    Source(FissionError),
}
impl From<DiffusionError> for NeutronKineticError {
    fn from(error: DiffusionError) -> Self {
        Self::Diffusion(error)
    }
}
impl From<FissionError> for NeutronKineticError {
    fn from(error: FissionError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for NeutronKineticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "neutron kinetic energy: {self:?}")
    }
}
impl std::error::Error for NeutronKineticError {}

fn checked(value: f64) -> Result<f64, NeutronKineticError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NeutronKineticError::Arithmetic(
            "unrepresentable energy or rate",
        ))
    }
}
fn weighted_sum(counts: &[f64], energies: &[f64]) -> Result<f64, NeutronKineticError> {
    let mut total = Accumulator::new(0.0);
    for (&count, &energy) in counts.iter().zip(energies.iter().cycle()) {
        total.add(checked(count * energy)?)?;
    }
    Ok(total.value())
}

#[cfg(test)]
mod tests;

#[cfg(test)]
pub(crate) use tests::material_induced_fission_fixture;
