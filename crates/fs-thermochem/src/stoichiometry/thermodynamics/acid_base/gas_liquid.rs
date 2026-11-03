//! Finite headspace partition coupled to all aqueous acid/base families.
//!
//! A neutral volatile family member obeys the original source-bound relation
//! `a_neutral = K p_partial/p_ref`. At a trial proton activity, let alpha be
//! that member's fraction of the dissolved family. The ideal gas then holds
//! `n_g/n_aqueous = alpha V_g p_ref/(n_ref R T K)`. Normalizing the gas and
//! aqueous weights closes each finite analytical total before evaluating the
//! original electroneutrality residual. Thus ionization changes total aqueous
//! uptake and gas pressure while pH responds to that uptake, in ONE solve.
//!
//! The gas member has zero charge; including it in the family normalization
//! preserves monotonic charge versus pH (the derivative is minus a charge
//! variance). No extra outer iteration or independent gas/pH approximation is
//! needed. Both phase tails are represented without subtracting from a total.
//!
//! This is a fixed-temperature, fixed-volume ideal gas and ideal aqueous model.
//! It adds no nonideal activity, gas EOS, salt effect, finite solvent balance,
//! gas-transfer kinetics, thermal duty, redox or competing mineral chemistry.
//! Fixed-volume and fixed-pressure gas phases are distinct formulations; see
//! USGS PHREEQC, gas-component mass action and fixed-volume inventories:
//! <https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-10.html>.

use super::super::gas_liquid::IdealGasLiquidPartition;
use super::*;
use crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use fs_math::det;
use fs_qty::{Pressure, Volume};

#[derive(Clone, Debug)]
struct Binding {
    family: usize,
    neutral: usize,
    partition: IdealGasLiquidPartition,
}

/// Original aqueous and neutral gas owners, bound by exact source identities.
#[derive(Clone, Debug)]
pub struct VolatileAcidBaseSystem {
    aqueous: IdealAcidBaseSystem,
    bindings: Vec<Binding>,
    network: ReactionNetwork,
}

/// One volatile analytical family across the finite gas and every aqueous form.
#[derive(Clone, Debug)]
pub struct VolatileFamilySpeciation {
    /// Original family identity, in the aqueous model's canonical order.
    pub family: String,
    /// Original neutral gas phase identity.
    pub gas_species: SpeciesId,
    /// Original neutral aqueous member to which Henry partition applies.
    pub neutral_aqueous_species: SpeciesId,
    /// Supplied total of gas plus ALL dissolved protonation forms [mol].
    pub component_total: Amount,
    /// Retained gas inventory [mol], including representable trace tails.
    pub gas_amount: Amount,
    /// Sum of all aqueous family forms [mol], not just the neutral molecule.
    pub aqueous_amount: Amount,
    /// This species' ideal partial pressure in the common headspace [Pa].
    pub partial_pressure: Pressure,
    /// Independently reduced gas plus aqueous minus supplied total [mol].
    pub component_balance_residual: Amount,
    /// Recomputed ln(a_neutral / (K p_partial/p_ref)); absent for an empty total.
    pub log_equilibrium_residual: Option<Dimensionless>,
}

/// One simultaneous finite-gas and aqueous equilibrium at the source temperature.
#[derive(Clone, Debug)]
pub struct VolatileAcidBaseSpeciation {
    model: VolatileAcidBaseSystem,
    /// Original complete aqueous result. Family totals here are DISSOLVED only.
    pub aqueous: AcidBaseSpeciation,
    /// Common fixed ideal headspace volume [m3].
    pub gas_volume: Volume,
    /// Canonically ordered volatile-family inventories and independent balances.
    pub volatile_families: Vec<VolatileFamilySpeciation>,
}
impl VolatileAcidBaseSpeciation {
    /// Original equations, source receipts and gas/aqueous bindings.
    pub fn model(&self) -> &VolatileAcidBaseSystem {
        &self.model
    }
}

impl VolatileAcidBaseSystem {
    /// Bind one or more original gas partitions to their neutral aqueous family
    /// members. The member may be an acid, a base or an intermediate form; its
    /// charge and composition come from the original exact reaction network.
    /// Each family has at most one neutral volatile member and one gas phase.
    /// All shared source answers/receipts and standard references must match.
    pub fn new(
        aqueous: IdealAcidBaseSystem,
        partitions: Vec<IdealGasLiquidPartition>,
    ) -> Result<Self> {
        if partitions.is_empty() || partitions.len() > aqueous.families().len() {
            return Err(AcidBaseError::Invalid(
                "require one neutral gas partition per selected distinct aqueous family",
            ));
        }
        let mut bindings = Vec::with_capacity(partitions.len());
        let mut gases = BTreeSet::new();
        let mut families = BTreeSet::new();
        for partition in partitions {
            let (family, neutral) = aqueous
                .families()
                .iter()
                .enumerate()
                .find_map(|(index, family)| {
                    family
                        .species()
                        .iter()
                        .position(|id| id == partition.aqueous())
                        .map(|neutral| (index, neutral))
                })
                .ok_or(AcidBaseError::Invalid(
                    "gas partition must name an original aqueous family member",
                ))?;
            if aqueous.families()[family].steps[0].acid_charge - neutral as i32 != 0
                || !families.insert(family)
                || !gases.insert(partition.gas().clone())
                || aqueous
                    .reaction_network()
                    .species()
                    .binary_search(partition.gas())
                    .is_ok()
            {
                return Err(AcidBaseError::Invalid(
                    "require distinct gas identities and one neutral member per volatile family",
                ));
            }
            bindings.push(Binding {
                family,
                neutral,
                partition,
            });
        }
        bindings.sort_by_key(|binding| binding.family);
        let standards = core::iter::once(aqueous.water().standard())
            .chain(
                aqueous
                    .families()
                    .iter()
                    .flat_map(|f| f.steps().iter().map(|s| s.standard())),
            )
            .chain(bindings.iter().map(|b| b.partition.standard()));
        consistent_sources(standards)?;
        let reactions = aqueous
            .reaction_network()
            .reactions()
            .iter()
            .cloned()
            .chain(bindings.iter().map(|b| b.partition.reaction().clone()))
            .collect();
        let network = ReactionNetwork::new(reactions)?;
        Ok(Self {
            aqueous,
            bindings,
            network,
        })
    }

    /// Original aqueous-family owner, including all unmodified nonvolatile pools.
    pub fn aqueous_model(&self) -> &IdealAcidBaseSystem {
        &self.aqueous
    }
    /// Original gas transfer owners in canonical aqueous-family order.
    pub fn partitions(&self) -> impl ExactSizeIterator<Item = &IdealGasLiquidPartition> {
        self.bindings.iter().map(|binding| &binding.partition)
    }
    /// Exact conserved equations with shared atom and charge identities checked.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.network
    }

    /// Solve one common proton activity, gas pressure and complete dissolved
    /// speciation. Supply every family total exactly once: a bound family's
    /// total includes its gas plus ALL aqueous forms; other totals are aqueous.
    /// The declared solvent activity basis is common and independent of volume.
    /// Spectator counterions remain explicit; no background ions are invented.
    ///
    /// This uses the original bounded scalar pH solve, with work linear in the
    /// number of family members per evaluation and no per-iteration allocation.
    /// Its charge, source-domain, bracket, numerical-range and cancellation
    /// refusals remain atomic. Residuals are numerical diagnostics, not source
    /// uncertainty propagation or experimental validation of Henry's law.
    #[allow(clippy::too_many_arguments)]
    pub fn speciate(
        &self,
        totals: &[(String, Amount)],
        spectators: &[SpectatorIon],
        amount_at_unit_activity: Amount,
        gas_volume: Volume,
        options: SpeciationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VolatileAcidBaseSpeciation> {
        if cancelled() {
            return Err(AcidBaseError::Cancelled);
        }
        let basis = amount_at_unit_activity.value();
        let volume = gas_volume.value();
        if !basis.is_finite() || basis <= 0.0 || !volume.is_finite() || volume <= 0.0 {
            return Err(AcidBaseError::Invalid(
                "positive aqueous activity basis and finite headspace volume required",
            ));
        }
        if spectators.iter().any(|ion| {
            self.bindings
                .iter()
                .any(|binding| &ion.species == binding.partition.gas())
        }) {
            return Err(AcidBaseError::Invalid(
                "a neutral gas identity cannot also be an aqueous spectator ion",
            ));
        }
        let mut weights = vec![None; self.aqueous.families().len()];
        for binding in &self.bindings {
            let standard = binding.partition.standard();
            let mut log_ratio = Sum::default();
            for term in [
                det::ln(volume),
                det::ln(standard.reference_pressure().value()),
                -det::ln(basis),
                -det::ln(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
                -det::ln(standard.temperature().value()),
                -standard.log_equilibrium_constant().value(),
            ] {
                log_ratio.add(term)?;
            }
            weights[binding.family] = Some(GasWeight {
                neutral: binding.neutral,
                log_gas_to_neutral: log_ratio.finish()?,
            });
        }
        let aqueous = self.aqueous.speciate_partitioned(
            totals,
            spectators,
            amount_at_unit_activity,
            options,
            &weights,
            &mut cancelled,
        )?;
        let mut volatile_families = Vec::with_capacity(self.bindings.len());
        for binding in &self.bindings {
            if cancelled() {
                return Err(AcidBaseError::Cancelled);
            }
            let family = &self.aqueous.families()[binding.family];
            // The original solve has admitted each total and canonical family.
            let total = totals.iter().find(|(id, _)| id == family.id()).unwrap().1;
            let fractions =
                solve::fractions(family, -aqueous.ph.value() * core::f64::consts::LN_10)?;
            let [_, gas] = weights[binding.family]
                .unwrap()
                .partition(total.value(), fractions[binding.neutral])?;
            let dissolved = &aqueous.families[binding.family];
            let pressure = binding
                .partition
                .ideal_partial_pressure(Amount::new(gas), gas_volume)?;
            let mut balance = Sum::default();
            balance.add(gas)?;
            for form in &dissolved.forms {
                balance.add(form.amount.value())?;
            }
            balance.add(-total.value())?;
            let residual = if total.value() == 0.0 {
                None
            } else {
                let mut value = Sum::default();
                for term in [
                    det::ln(dissolved.forms[binding.neutral].activity.value()),
                    -det::ln(pressure.value()),
                    det::ln(binding.partition.standard().reference_pressure().value()),
                    -binding
                        .partition
                        .standard()
                        .log_equilibrium_constant()
                        .value(),
                ] {
                    value.add(term)?;
                }
                Some(Dimensionless::new(value.finish()?))
            };
            volatile_families.push(VolatileFamilySpeciation {
                family: family.id().into(),
                gas_species: binding.partition.gas().clone(),
                neutral_aqueous_species: binding.partition.aqueous().clone(),
                component_total: total,
                gas_amount: Amount::new(gas),
                aqueous_amount: dissolved.total,
                partial_pressure: pressure,
                component_balance_residual: Amount::new(balance.finish()?),
                log_equilibrium_residual: residual,
            });
        }
        if cancelled() {
            return Err(AcidBaseError::Cancelled);
        }
        Ok(VolatileAcidBaseSpeciation {
            model: self.clone(),
            aqueous,
            gas_volume,
            volatile_families,
        })
    }
}

/// Fixed neutral gas weight; the aqueous ionization fraction varies with pH.
#[derive(Clone, Copy, Debug)]
pub(super) struct GasWeight {
    pub neutral: usize,
    log_gas_to_neutral: f64,
}
impl GasWeight {
    /// Return [aqueous, gas] without cancellation of a positive trace stock.
    pub(super) fn partition(self, total: f64, neutral_fraction: f64) -> Result<[f64; 2]> {
        if total == 0.0 {
            return Ok([0.0; 2]);
        }
        let log_ratio = self.log_gas_to_neutral + det::ln(neutral_fraction);
        if !log_ratio.is_finite() {
            return Err(AcidBaseError::Arithmetic(
                "gas/aqueous log partition loses range",
            ));
        }
        let small = det::exp(-log_ratio.abs());
        let major = positive(total / (1.0 + small), "dominant volatile stock loses range")?;
        let mut minor = major * small;
        if minor == 0.0 {
            minor = det::exp(det::ln(total) - log_ratio.abs() - det::ln(1.0 + small));
        }
        positive(minor, "positive volatile trace stock loses range")?;
        Ok(if log_ratio >= 0.0 {
            [minor, major]
        } else {
            [major, minor]
        })
    }
}
