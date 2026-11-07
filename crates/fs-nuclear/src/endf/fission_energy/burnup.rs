//! Evaluated prompt fission products on the original burnup and heat owners.
//!
//! Only fragment kinetic energy, prompt-neutron kinetic energy and prompt-photon
//! energy enter these fission-event coefficients. Delayed neutrons, photons,
//! betas and neutrinos remain separate MT458 observations: represented product
//! decays already release their own energy through the original reaction model.
//!
//! The coefficients below describe OUTGOING PRODUCT ENERGY. That includes the
//! incident neutron's kinetic input; it is not all newly generated nuclear Q.
//! Incident energy and prompt Q are returned separately. No arbitrary debit is
//! taken from a material or radiation channel. A coupled neutron kinetic-energy
//! owner must retire the same incident energy, and must not credit transported
//! outgoing neutrons a second time. A population-only prescribed-flux model has
//! no such kinetic-energy ledger; its external neutron supply carries that input.

use super::FissionEnergyComponent;
use super::model::CollapsedFissionEnergy;
use crate::burnup::{BurnupData, BurnupLimits};
use crate::cross_sections::NeutronReaction;
use crate::endf::fission_yields::model::FissionYieldCollapse;
use crate::endf::{EndfError, budget};
use crate::irradiation::NuclideId;
use crate::{Accumulator, EnergyChannels, FissionError, finite, nonnegative, sum};

/// Explicit destination of an entire prompt energy component. Fractional
/// transport can subsequently use the original `RadiationResponse` owner.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FissionEnergyDestination {
    /// Thermalized at the original event's material cell.
    Deposited,
    /// Leaves the event cell for separately owned particle/radiation transport.
    Radiation,
    /// Leaves the represented system without subsequent modeled recovery.
    Escaped,
}

impl FissionEnergyDestination {
    const fn index(self) -> usize {
        match self {
            Self::Deposited => 0,
            Self::Radiation => 1,
            Self::Escaped => 2,
        }
    }
}

/// The caller's treatment of the three prompt MT458 components. Neither local
/// absorption nor neutron/photon transport is inferred from the evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PromptFissionRouting {
    pub fragments: FissionEnergyDestination,
    pub prompt_neutrons: FissionEnergyDestination,
    pub prompt_photons: FissionEnergyDestination,
}

/// One reacting neutron group's energy per fission, in joules. These are
/// reaction-weighted means for the original within-group spectrum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PromptFissionEnergy {
    /// The outgoing prompt products, each assigned once to its declared route.
    pub channels: EnergyChannels,
    pub outgoing_product_energy_j: f64,
    pub incident_neutron_energy_j: f64,
    /// Outgoing prompt products minus incident kinetic energy. This diagnostic
    /// is signed; delayed release is deliberately excluded.
    pub prompt_q_j: f64,
    /// Noncancelling aggregation and scalar-export defects in this routing.
    /// This does not include evaluated-data uncertainty or upstream quadrature.
    pub roundoff_j: f64,
    /// Sum of the prepared prompt component interpolation-error estimates.
    /// This is a coefficient approximation, separate from runtime energy closure.
    pub interpolation_error_j: f64,
}

impl CollapsedFissionEnergy {
    /// Route prompt products without changing the evaluated source components.
    /// A zero-reaction group stays `None`; no mean energy is invented for it.
    pub fn prompt_products(
        &self,
        routing: PromptFissionRouting,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<Option<PromptFissionEnergy>>, EndfError> {
        checkpoint()?;
        let components = [
            FissionEnergyComponent::Fragments,
            FissionEnergyComponent::PromptNeutrons,
            FissionEnergyComponent::PromptPhotons,
        ];
        let destinations = [
            routing.fragments,
            routing.prompt_neutrons,
            routing.prompt_photons,
        ];
        let weighted = self.reaction_weighted();
        let groups = weighted.flux.groups().len();
        if weighted.reaction != NeutronReaction::Fission
            || weighted.cross_sections_m2.len() != groups
            || self.incident_energy_j().len() != groups
            || components
                .iter()
                .any(|&component| self.component_j(component).len() != groups)
        {
            return Err(invalid(
                "matching fission component and neutron-group shape required",
            ));
        }
        let mut interpolation_error = NonnegativeBound::default();
        for component in components {
            let error = self.interpolation_error_j(component);
            if !nonnegative(error) {
                return Err(invalid(
                    "finite nonnegative component interpolation error required",
                ));
            }
            interpolation_error.add(error)?;
        }
        let interpolation_error_j = interpolation_error.0;
        let mut output = Vec::new();
        output
            .try_reserve_exact(groups)
            .map_err(|_| budget("prompt fission energy groups"))?;
        for (group, &cross_section) in weighted.cross_sections_m2.iter().enumerate() {
            checkpoint()?;
            if !nonnegative(cross_section) {
                return Err(invalid("finite nonnegative fission cross section required"));
            }
            let values = components.map(|component| self.component_j(component)[group]);
            let incident = self.incident_energy_j()[group];
            if cross_section == 0.0 {
                if incident.is_some() || values.iter().any(Option::is_some) {
                    return Err(invalid(
                        "zero-reaction group cannot carry a fission energy mean",
                    ));
                }
                output.push(None);
                continue;
            }
            let [Some(fragment), Some(neutron), Some(photon)] = values else {
                return Err(invalid(
                    "reacting group requires all three prompt energy components",
                ));
            };
            let Some(incident_neutron_energy_j) = incident else {
                return Err(invalid(
                    "reacting group requires its incident neutron energy mean",
                ));
            };
            let prompt = [fragment, neutron, photon];
            if !prompt.into_iter().all(nonnegative) || !nonnegative(incident_neutron_energy_j) {
                return Err(invalid(
                    "finite nonnegative prompt and incident energies required",
                ));
            }
            let mut channels = [Accumulator::new(0.0); 3];
            let mut outgoing = Accumulator::new(0.0);
            let mut q = Accumulator::new(0.0);
            let mut roundoff = NonnegativeBound::default();
            for (energy, destination) in prompt.into_iter().zip(destinations) {
                roundoff.add(channels[destination.index()].add(energy)?.abs())?;
                roundoff.add(outgoing.add(energy)?.abs())?;
                roundoff.add(q.add(energy)?.abs())?;
            }
            roundoff.add(q.add(-incident_neutron_energy_j)?.abs())?;
            for value in channels.into_iter().chain([outgoing, q]) {
                roundoff.add(export_error(value)?)?;
            }
            let channels = EnergyChannels {
                deposited_j: finite(channels[0].value())?,
                radiation_j: finite(channels[1].value())?,
                escaped_j: finite(channels[2].value())?,
            };
            channels.validate()?;
            let outgoing_product_energy_j = finite(outgoing.value())?;
            let prompt_q_j = finite(q.value())?;
            // Charge the independently exported scalar balances as well as the
            // retained low parts. Opposite defects cannot cancel this account.
            roundoff.add(
                finite(sum(&[
                    channels.deposited_j,
                    channels.radiation_j,
                    channels.escaped_j,
                    -outgoing_product_energy_j,
                ]))?
                .abs(),
            )?;
            roundoff.add(
                finite(sum(&[
                    prompt_q_j,
                    incident_neutron_energy_j,
                    -outgoing_product_energy_j,
                ]))?
                .abs(),
            )?;
            output.push(Some(PromptFissionEnergy {
                channels,
                outgoing_product_energy_j,
                incident_neutron_energy_j,
                prompt_q_j,
                roundoff_j: roundoff.0,
                interpolation_error_j,
            }));
        }
        checkpoint()?;
        Ok(output)
    }

    /// Append evaluated OUTGOING prompt product energy to the original disjoint
    /// fission-product reactions. The original input data are never mutated.
    ///
    /// Both condensations must use the same original parent, processed metadata,
    /// group grid, cross-section/spectrum identities, original pointwise knots
    /// and interpolation laws, numerical cross sections and absolute group
    /// fluxes. Matching group means alone cannot establish matching spectral
    /// shapes. These checks do not certify evaluated-data accuracy.
    ///
    /// `maximum_routing_roundoff_j` is an explicit J/fission allowance for the
    /// coefficient aggregation above. It is not a whole-run energy-error bound:
    /// actual fission counts scale coefficient errors, including the separately
    /// reported interpolation error. Runtime accounts retain their own budgets.
    ///
    /// Only EFR, ENP and EGP are inserted. END, EGD, EB and ENU are not released
    /// immediately and their lifetimes are not inferred from MT458. Product
    /// decay and shutdown heat retain the original inventory/thermal owners.
    #[allow(clippy::too_many_arguments)]
    pub fn append_prompt_products_to_burnup(
        &self,
        yields: &FissionYieldCollapse,
        data: &BurnupData,
        routing: PromptFissionRouting,
        exported_products: &[NuclideId],
        limits: BurnupLimits,
        maximum_routing_roundoff_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupData, EndfError> {
        let channels = self.prompt_channels_for_burnup(
            yields,
            routing,
            limits,
            maximum_routing_roundoff_j,
            &mut checkpoint,
        )?;
        let mut assembled =
            yields.append_to_burnup(data, &channels, exported_products, limits, &mut checkpoint)?;
        for reaction in &mut assembled.reactions[data.reactions.len()..] {
            checkpoint()?;
            reaction.data_id = format!(
                "{} | {} | MF1/MT458 outgoing prompt products",
                reaction.data_id,
                self.evaluation().data_id()
            );
        }
        checkpoint()?;
        Ok(assembled)
    }

    /// Shared validation and routing for fixed-state and temperature-series
    /// assembly. Neither caller may weaken original pointwise source matching.
    pub(super) fn prompt_channels_for_burnup(
        &self,
        yields: &FissionYieldCollapse,
        routing: PromptFissionRouting,
        limits: BurnupLimits,
        maximum_routing_roundoff_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<EnergyChannels>, EndfError> {
        checkpoint()?;
        if !nonnegative(maximum_routing_roundoff_j) {
            return Err(invalid(
                "finite nonnegative prompt routing allowance required",
            ));
        }
        let energy = self.reaction_weighted();
        let products = yields.reaction_weighted();
        if energy.flux.groups().len() > limits.max_neutron_groups {
            return Err(budget("prompt fission energy groups"));
        }
        let energy_cross = self.source_cross_section();
        let product_cross = yields.source_cross_section();
        let energy_spectrum = self.source_spectrum();
        let product_spectrum = yields.source_spectrum();
        if yields.parent() != self.evaluation().metadata().nuclide
            || energy.reaction != NeutronReaction::Fission
            || products.reaction != NeutronReaction::Fission
            || energy.metadata != products.metadata
            || energy.cross_section_data_id != products.cross_section_data_id
            || energy.spectrum_data_id != products.spectrum_data_id
            || energy.flux.groups().data_id() != products.flux.groups().data_id()
            || !same_values(
                energy.flux.groups().edges_j(),
                products.flux.groups().edges_j(),
                &mut checkpoint,
            )?
            || !same_values(
                &energy.cross_sections_m2,
                &products.cross_sections_m2,
                &mut checkpoint,
            )?
            || !same_values(
                energy.flux.values_m2_s(),
                products.flux.values_m2_s(),
                &mut checkpoint,
            )?
            || energy_cross.data_id() != product_cross.data_id()
            || energy_cross.metadata() != product_cross.metadata()
            || energy_cross.reaction() != product_cross.reaction()
            || energy_cross.mass_difference_q_j() != product_cross.mass_difference_q_j()
            || energy_cross.reaction_q_j() != product_cross.reaction_q_j()
            || !same_values(
                energy_cross.energies_j(),
                product_cross.energies_j(),
                &mut checkpoint,
            )?
            || !same_values(
                energy_cross.cross_sections_m2(),
                product_cross.cross_sections_m2(),
                &mut checkpoint,
            )?
            || !same_values(
                energy_cross.interpolation_regions(),
                product_cross.interpolation_regions(),
                &mut checkpoint,
            )?
            || energy_spectrum.data_id() != product_spectrum.data_id()
            || !same_values(
                energy_spectrum.energies_j(),
                product_spectrum.energies_j(),
                &mut checkpoint,
            )?
            || !same_values(
                energy_spectrum.flux_density_m2_s_j(),
                product_spectrum.flux_density_m2_s_j(),
                &mut checkpoint,
            )?
            || !same_values(
                energy_spectrum.interpolation_regions(),
                product_spectrum.interpolation_regions(),
                &mut checkpoint,
            )?
        {
            return Err(invalid(
                "fission energy and products require the same parent and processed reaction/spectrum condensation",
            ));
        }
        let prompt = self.prompt_products(routing, &mut checkpoint)?;
        let mut channels = Vec::new();
        channels
            .try_reserve_exact(prompt.len())
            .map_err(|_| budget("prompt fission energy channels"))?;
        for group in prompt {
            checkpoint()?;
            match group {
                Some(group) => {
                    if group.roundoff_j > maximum_routing_roundoff_j {
                        return Err(budget("prompt fission coefficient routing roundoff"));
                    }
                    channels.push(group.channels);
                }
                // The original bridge omits this zero-cross-section reaction.
                None => channels.push(EnergyChannels::default()),
            }
        }
        checkpoint()?;
        Ok(channels)
    }
}

pub(super) fn same_values<T: PartialEq>(
    left: &[T],
    right: &[T],
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    checkpoint()?;
    if left.len() != right.len() {
        return Ok(false);
    }
    // Original source clones share immutable Arc payloads, so the usual path
    // requires no linear scan. Distinct payloads remain cancellable in batches.
    if std::ptr::eq(left, right) {
        return Ok(true);
    }
    for (left, right) in left.chunks(64).zip(right.chunks(64)) {
        checkpoint()?;
        if left != right {
            return Ok(false);
        }
    }
    Ok(true)
}

fn export_error(value: Accumulator) -> Result<f64, FissionError> {
    Ok(finite(sum(&[value.hi, value.lo, -finite(value.value())?]))?.abs())
}

/// A nonnegative error account rounds outwards after every addition. Merely
/// rounding a many-term expansion at its final export could miss earlier
/// positive defects below that expansion's retained precision.
#[derive(Default)]
struct NonnegativeBound(f64);

impl NonnegativeBound {
    fn add(&mut self, value: f64) -> Result<(), FissionError> {
        if !nonnegative(value) {
            return Err(FissionError::Arithmetic);
        }
        // Two input scalars need at most the original two-term expansion.
        let mut pair = Accumulator::new(self.0);
        pair.add(value)?;
        self.0 = nonnegative_upper_bound(pair)?;
        Ok(())
    }
}

fn nonnegative_upper_bound(value: Accumulator) -> Result<f64, FissionError> {
    let rounded = finite(value.value())?;
    if finite(sum(&[value.hi, value.lo, -rounded]))? > 0.0 {
        finite(rounded.next_up())
    } else {
        Ok(rounded)
    }
}

fn invalid(reason: &'static str) -> EndfError {
    FissionError::InvalidInput(reason).into()
}

#[cfg(test)]
#[path = "burnup_tests.rs"]
mod tests;
