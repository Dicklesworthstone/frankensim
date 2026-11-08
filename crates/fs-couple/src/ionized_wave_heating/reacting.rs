//! Held Maxwell work on the original finite-rate atomic energy/population owner.

use super::{IonizedWaveHeatingError, MagnetizedIonizedGasHeating, poll};
use fs_gas::chemistry::collisional_radiative::CollisionalRadiativeThermal;
use fs_gas::chemistry::energy::ElectronHeavyCollisionSources;
use fs_gas::chemistry::ionization_network::{
    IonizationNetworkError, IonizationNetworkResult, IonizationNetworkSolve,
    advance_ionization_network_with_thermal,
};
use fs_matdb::MaterialCard;
use fs_material::ionized_gas::{
    IonizedGasConductivity, IonizedGasMomentumSource, IonizedGasRefreshError, IonizedGasTransportBank,
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::{Power, Time};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;

/// Native reaction controls and explicit limits on one held-field interval.
///
/// Temperature and amount limits compare the accepted endpoints. They are
/// modeling limits, not error estimates or proof that intermediate transport
/// and optical coefficients remain close to the initial ones. A new Maxwell
/// solve is required after every accepted interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizedGasReactingHeatControls {
    /// Original reaction/radiation source, time-error and work allowances.
    pub network: IonizationNetworkSolve,
    /// Original heavy-species caloric support [K], checked by the gas owner.
    pub minimum_heavy_temperature_k: f64,
    /// Upper bound of the same original heavy-species support [K].
    pub maximum_heavy_temperature_k: f64,
    /// Work bound for each source heavy-caloric inversion, in 8..=256.
    pub maximum_heavy_caloric_iterations: usize,
    /// Largest absolute change of either pool temperature [K].
    pub maximum_temperature_increment_k: f64,
    /// Absolute allowance in each species' amount change [mol], nonnegative.
    /// This explicitly controls creation of initially absent species.
    pub maximum_absolute_population_increment_mol: f64,
    /// Additional allowance times that species' ORIGINAL amount, in (0,1].
    /// Trace carrier changes are not normalized by the much larger neutral pool.
    pub maximum_relative_population_increment: f64,
}

impl<'w, 'a> MagnetizedIonizedGasHeating<'w, 'a> {
    /// Deposit the wave's two signed powers once while evolving the original
    /// atomic charge stages, electron energy, and heavy temperature together.
    /// Impact ionization spends source level energy; three-body recombination
    /// returns it; spontaneous recombination radiation escapes. The existing
    /// source-bound network retains its temporal and energy error controls.
    ///
    /// `transport_bank` must reproduce every current wave carrier's original
    /// momentum source, optional Coulomb law and transport limits. Its complete
    /// neutral cross sections also supply thermal exchange for targets the
    /// chemistry may create. Neither cached collision rates nor caller-invented
    /// heat partitions enter the gas evolution. Coulomb drift heat is ALREADY
    /// in the held powers; the gas owner adds only baseline thermal exchange.
    ///
    /// The coherent wave is a collisional work source only. Its photons are
    /// not separately counted as a photoionization beam. Optical recombination
    /// transport, line/free-free radiation and continuous field feedback are
    /// outside this interval model. Bound-charge background and ion momentum
    /// cards must remain applicable to the resulting mixture, as authored;
    /// endpoint T/p refresh alone does not establish that physical validity.
    ///
    /// All source/bound checks, chemistry, and next material evaluations succeed
    /// before a result is published. Refusal leaves the original wave, gas and
    /// complete transport bank unchanged. An exactly neutral endpoint returns
    /// no free-carrier transport; no conductivity floor is introduced.
    pub fn advance_reacting<'h>(
        &'h self,
        duration: Time,
        sources: &[AtomicRadiativeRecombinationV1],
        transport_bank: &'h IonizedGasTransportBank,
        controls: IonizedGasReactingHeatControls,
        original_background_card: &MaterialCard,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MagnetizedIonizedGasReactingHeatStep<'h, 'w, 'a>, IonizedWaveHeatingError> {
        use IonizedWaveHeatingError::Invalid;
        poll(&mut cancelled)?;
        if !duration.value().is_finite()
            || duration.value() <= 0.0
            || !controls.maximum_temperature_increment_k.is_finite()
            || controls.maximum_temperature_increment_k <= 0.0
            || !controls.maximum_absolute_population_increment_mol.is_finite()
            || controls.maximum_absolute_population_increment_mol < 0.0
            || !controls.maximum_relative_population_increment.is_finite()
            || controls.maximum_relative_population_increment <= 0.0
            || controls.maximum_relative_population_increment > 1.0
        {
            return Err(Invalid(
                "positive duration and explicit held-field increment limits required",
            ));
        }
        let source = self.source();
        let initial = source.state();
        let background = self.absorption().wave().interface().background();
        if original_background_card.content_hash() != background.card_identity() {
            return Err(Invalid(
                "bound-charge refresh requires its original immutable material card",
            ));
        }
        let bound = transport_bank
            .evaluate(initial, &mut cancelled)
            .map_err(transport_error)?
            .ok_or(Invalid(
                "transport bank cannot remove the original wave carriers",
            ))?;
        if bound.carriers() != source.carriers()
            || bound.limits() != source.limits()
            || bound.coulomb_model() != source.coulomb_model()
        {
            return Err(Invalid(
                "transport bank must reproduce the original wave carrier laws and limits",
            ));
        }
        let electron_index = initial
            .electron_index()
            .ok_or(Invalid("structural electron source required"))?;
        let electron = initial.thermodynamics().models()[electron_index].species();
        let electron_source = transport_bank
            .sources()
            .find(|(species, _)| *species == electron)
            .map(|(_, source)| source)
            .ok_or(Invalid("complete original electron momentum law required"))?;
        let IonizedGasMomentumSource::ElectronNeutral {
            cross_sections,
            controls: rate_controls,
            maximum_target_thermal_speed_ratio,
        } = electron_source
        else {
            return Err(Invalid(
                "reacting thermal exchange requires the original target-resolved cross sections",
            ));
        };
        let electron_power = self.absorption().electron_power_w();
        let heavy_power = self.absorption().heavy_power_w();
        let absorbed_energy_j = self.absorption().absorbed_power_w() * duration.value();
        for power in [
            electron_power,
            heavy_power,
            self.absorption().absorbed_power_w(),
        ] {
            let energy = power * duration.value();
            if !energy.is_finite() || (power != 0.0 && energy == 0.0) {
                return Err(Invalid("physical wave heating energy is unrepresentable"));
            }
        }
        let thermal = CollisionalRadiativeThermal {
            collisions: Some(ElectronHeavyCollisionSources {
                momentum_sources: cross_sections.clone(),
                coulomb: source.coulomb_model().cloned(),
                maximum_charged_mole_fraction: source.limits().maximum_charged_mole_fraction,
                maximum_target_thermal_speed_ratio: *maximum_target_thermal_speed_ratio,
                rate_controls: *rate_controls,
            }),
            electron_power: Power::new(electron_power),
            heavy_power: Power::new(heavy_power),
            minimum_heavy_temperature_k: controls.minimum_heavy_temperature_k,
            maximum_heavy_temperature_k: controls.maximum_heavy_temperature_k,
            maximum_iterations: controls.maximum_heavy_caloric_iterations,
        };
        let chemistry = advance_ionization_network_with_thermal(
            initial,
            sources,
            &[],
            &thermal,
            duration,
            controls.network,
            &mut cancelled,
        )
        .map_err(|error| match error {
            IonizationNetworkError::Cancelled => IonizedWaveHeatingError::Cancelled,
            error => IonizedWaveHeatingError::Reaction(error),
        })?;
        let next = chemistry.state();
        let heavy_k =
            (next.heavy_temperature().value() - initial.heavy_temperature().value()).abs();
        let electron_k =
            (next.electron_temperature().value() - initial.electron_temperature().value()).abs();
        if heavy_k.max(electron_k) > controls.maximum_temperature_increment_k {
            return Err(IonizedWaveHeatingError::TemperatureIncrement {
                heavy_k,
                electron_k,
                maximum_k: controls.maximum_temperature_increment_k,
            });
        }
        for (i, (new, old)) in next.amounts().iter().zip(initial.amounts()).enumerate() {
            poll(&mut cancelled)?;
            let change_mol = (new.value() - old.value()).abs();
            let maximum_mol = controls.maximum_absolute_population_increment_mol
                + controls.maximum_relative_population_increment * old.value();
            if !maximum_mol.is_finite() {
                return Err(Invalid("population increment allowance is unrepresentable"));
            }
            if change_mol > maximum_mol {
                return Err(IonizedWaveHeatingError::PopulationIncrement {
                    species: initial.thermodynamics().models()[i].species().clone(),
                    change_mol,
                    maximum_mol,
                });
            }
        }
        let refreshed_transport = transport_bank
            .evaluate(next, &mut cancelled)
            .map_err(transport_error)?;
        let mut query = background.query();
        for (axis, _) in background.query_point() {
            let value = match axis.as_str() {
                "T" | "T_h" => next.heavy_temperature().value(),
                "T_e" => next.electron_temperature().value(),
                "p" => next.pressure().value(),
                _ => return Err(Invalid("unmapped original plasma background condition")),
            };
            query = query
                .with(axis.clone(), value)
                .map_err(|_| Invalid("refreshed bound-charge coordinate is inadmissible"))?;
        }
        let refreshed_background = background
            .at_query_point(original_background_card, &query)
            .map_err(IonizedWaveHeatingError::Background)?;
        poll(&mut cancelled)?;
        Ok(MagnetizedIonizedGasReactingHeatStep {
            heating: self,
            transport_bank,
            chemistry,
            refreshed_transport,
            refreshed_background,
            absorbed_energy_j,
            controls,
        })
    }
}

/// One atomic endpoint and the original laws for the next Maxwell preparation.
pub struct MagnetizedIonizedGasReactingHeatStep<'h, 'w, 'a> {
    heating: &'h MagnetizedIonizedGasHeating<'w, 'a>,
    transport_bank: &'h IonizedGasTransportBank,
    chemistry: IonizationNetworkResult,
    refreshed_transport: Option<IonizedGasConductivity>,
    refreshed_background: ResolvedMaterialStatePoint,
    absorbed_energy_j: f64,
    controls: IonizedGasReactingHeatControls,
}
impl<'h, 'w, 'a> MagnetizedIonizedGasReactingHeatStep<'h, 'w, 'a> {
    /// Original retained Maxwell absorption and gas snapshot.
    pub const fn heating(&self) -> &'h MagnetizedIonizedGasHeating<'w, 'a> {
        self.heating
    }
    /// Complete original momentum laws, including currently absent species.
    pub const fn transport_bank(&self) -> &'h IonizedGasTransportBank {
        self.transport_bank
    }
    /// Native population, chemical/sensible energy, radiation and exchange accounts.
    pub const fn chemistry(&self) -> &IonizationNetworkResult {
        &self.chemistry
    }
    /// Actual endpoint carriers; `None` denotes exact free-carrier extinction.
    pub const fn refreshed_transport(&self) -> Option<&IonizedGasConductivity> {
        self.refreshed_transport.as_ref()
    }
    /// Original bound-charge card and claims at the endpoint temperatures/pressure.
    pub const fn refreshed_background(&self) -> &ResolvedMaterialStatePoint {
        &self.refreshed_background
    }
    /// Original Maxwell friction work, distinct from escaped recombination light.
    pub const fn absorbed_energy_j(&self) -> f64 {
        self.absorbed_energy_j
    }
    /// Requested source support, work budgets and held-field limits.
    pub const fn controls(&self) -> IonizedGasReactingHeatControls {
        self.controls
    }
}

fn transport_error(error: IonizedGasRefreshError) -> IonizedWaveHeatingError {
    match error {
        IonizedGasRefreshError::Cancelled => IonizedWaveHeatingError::Cancelled,
        error => IonizedWaveHeatingError::Transport(error),
    }
}
