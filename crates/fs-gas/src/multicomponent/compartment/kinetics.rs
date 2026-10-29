//! Finite-rate chemistry on the original two-zone gas inventories.
//!
//! Each selected, present layer first uses the original prepared
//! [`IsochoricReactor`] at its actual fixed volume and physical internal energy.
//! Both resulting species increments are then applied together through the
//! original common-pressure/interface-work solve. This is a first-order
//! homogeneous chemistry / mechanical-relaxation composition. The reactor's
//! entropy diagnostics concern its fixed-U,V stage, not a fully discrete entropy
//! theorem for the complete moving-interface or vent-coupled update.
//!
//! No formation energy, molecular data, Arrhenius law, reverse rate, nonlinear
//! solver, or heat of combustion is supplied a second time. Chemical rates carry
//! only the numerical-reference correction `-sum(dm_i/dt * e_ref_i)`; physical
//! internal energy has no chemical source. In a network pressure-stage callback,
//! call [`ZoneKinetics::apply_split`] on the transported trial and add each
//! returned chemical rate to the network's additional source account once.
//! Every pressure trial must start from the original transported inventories.
//!
//! Construction allocates the original reactor and one room candidate. Accepted
//! evaluations and updates reuse those buffers without allocating. An absent
//! layer remains exactly absent, and disabled layers receive no chemical rates.

use super::{Layer, LayerRates, MAX_SPECIES, TwoZone, ZoneError, ZoneSolve, ZoneStepReport, sum};
use crate::chemistry::{
    IsochoricKineticsBalance, IsochoricKineticsError, IsochoricKineticsSolve, IsochoricReactor,
};
use crate::multicomponent::MixtureGas;
use fs_qty::{Amount, Energy, Pressure, Temperature, Time, Volume};
use fs_thermochem::kinetics::IdealGasKineticsV1;

/// Refusal leaves the supplied room and its physical clock unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum ZoneKineticsError {
    /// Invalid interval, arithmetic, or original-domain relationship.
    Invalid(&'static str),
    /// Complete species/source bank or original energy references differ.
    SourceMismatch,
    /// Original prepared homogeneous reactor refused a layer endpoint.
    Chemistry {
        layer: Layer,
        source: IsochoricKineticsError,
    },
    /// Original reactor preparation refused the source or controls.
    Preparation(IsochoricKineticsError),
    /// Original pressure/interface solve or candidate identity refused.
    Zone(ZoneError),
    /// Rate arithmetic or accepted physical balance exceeded its numerical gate.
    Balance(&'static str),
    /// A selected layer leaves the admitted kinetic temperature interval.
    TemperatureDomain { layer: Layer, temperature_k: f64 },
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for ZoneKineticsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "two-zone finite-rate chemistry: {self:?}")
    }
}
impl std::error::Error for ZoneKineticsError {}
impl From<ZoneError> for ZoneKineticsError {
    fn from(error: ZoneError) -> Self {
        match error {
            ZoneError::Cancelled => Self::Cancelled,
            other => Self::Zone(other),
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ZoneKineticsError> {
    if cancelled() {
        Err(ZoneKineticsError::Cancelled)
    } else {
        Ok(())
    }
}

/// Original homogeneous-reactor diagnostics for one present, selected layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerKineticsReport {
    /// Fixed volume used by the chemistry substep, before interface relaxation.
    pub volume: Volume,
    /// One original physical internal-energy target, including formation terms.
    pub physical_internal_energy: Energy,
    /// Homogeneous endpoint temperature, before interface work.
    pub isochoric_temperature: Temperature,
    /// Homogeneous endpoint pressure, before mechanical equilibration.
    pub isochoric_pressure: Pressure,
    /// Original fixed-U,V conservation, entropy, local error and work report.
    pub reactor: IsochoricKineticsBalance,
    /// Largest actual rate-reconstructed versus reactor species-mass difference [kg].
    pub maximum_species_roundoff_kg: f64,
    /// Actual rate-integrated net molecular mass [kg], retained without repair.
    pub mass_source_kg: f64,
    /// Actual rate-integrated physical energy [J], nominally zero.
    pub physical_energy_source_j: f64,
}

/// Prepared finite-duration chemistry-only source; no room has been changed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZoneKineticsReport {
    /// Original lower/upper species rates and numerical energy-reference change.
    pub chemical_rates: [LayerRates; 2],
    /// None means absent or explicitly disabled, with exactly zero chemical rates.
    pub layers: [Option<LayerKineticsReport>; 2],
    /// Interval represented by these time-averaged rates.
    pub duration: Time,
}

/// One accepted chemistry/pressure-relaxation stage without another clock tick.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitKineticsReport {
    /// Original homogeneous reactor endpoints and actual chemical source rates.
    pub chemistry: ZoneKineticsReport,
    /// Original simultaneous common-pressure/interface-work solve.
    pub relaxation: ZoneStepReport,
    /// Final minus initial actual whole-room molecular mass [kg].
    pub mass_residual_kg: f64,
    /// Final minus initial actual whole-room physical internal energy [J].
    pub physical_energy_residual_j: f64,
}

/// Reusable chemistry workspace bound to one original gas and room geometry.
/// Supplied physical state is authoritative on every call, including retries.
pub struct ZoneKinetics {
    reactor: IsochoricReactor,
    candidate: TwoZone,
    amounts: [Amount; MAX_SPECIES],
}
impl ZoneKinetics {
    /// Bind a complete original mechanism and reserve the original reactor and
    /// room candidate. Kinetic temperatures must fit the gas's admitted domain.
    pub fn new(
        room: &TwoZone,
        mechanism: IdealGasKineticsV1,
        controls: IsochoricKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ZoneKineticsError> {
        poll(&mut cancelled)?;
        if room.gas().source() != mechanism.thermodynamics() {
            return Err(ZoneKineticsError::SourceMismatch);
        }
        let domain = room.gas().domain();
        if controls.minimum_temperature_k < domain.minimum_temperature
            || controls.maximum_temperature_k > domain.maximum_temperature
        {
            return Err(ZoneKineticsError::Invalid(
                "kinetic interval exceeds gas caloric domain",
            ));
        }
        let reactor =
            IsochoricReactor::new(mechanism, controls, &mut cancelled).map_err(|e| match e {
                IsochoricKineticsError::Cancelled => ZoneKineticsError::Cancelled,
                other => ZoneKineticsError::Preparation(other),
            })?;
        let result = Self {
            reactor,
            candidate: room.clone(),
            amounts: [Amount::new(0.0); MAX_SPECIES],
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Original complete source bank, caloric domain and energy-reference gauge.
    pub fn gas(&self) -> &MixtureGas {
        self.candidate.gas()
    }
    /// Exact original source mechanism, including reversible rate construction.
    pub fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.reactor.mechanism()
    }
    /// Original integration tolerances and bounded adaptive work controls.
    pub fn controls(&self) -> IsochoricKineticsSolve {
        self.reactor.controls()
    }

    /// Evaluate finite-duration chemistry from the actual current inventory.
    /// The layer selection permits chemistry in either/both layers; absent
    /// layers are skipped exactly. The room and accepted clock are not mutated.
    pub fn evaluate(
        &mut self,
        room: &TwoZone,
        duration: Time,
        selected: [bool; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ZoneKineticsReport, ZoneKineticsError> {
        poll(&mut cancelled)?;
        if room.gas() != self.gas() {
            return Err(ZoneKineticsError::SourceMismatch);
        }
        let dt = duration.value();
        if !dt.is_finite() || dt <= 0.0 {
            return Err(ZoneKineticsError::Invalid(
                "positive finite chemical duration required",
            ));
        }
        let controls = self.controls();
        let gas = room.gas();
        let n = gas.species_count();
        let mut result = ZoneKineticsReport {
            chemical_rates: [LayerRates::default(); 2],
            layers: [None; 2],
            duration,
        };
        for layer in [Layer::Lower, Layer::Upper] {
            poll(&mut cancelled)?;
            let k = layer.index();
            let state = room.layer(layer);
            if !selected[k] || state.mass_kg() == 0.0 {
                continue;
            }
            if !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
                .contains(&state.temperature().value())
            {
                return Err(ZoneKineticsError::TemperatureDomain {
                    layer,
                    temperature_k: state.temperature().value(),
                });
            }
            let initial = state.species_masses();
            for (i, amount) in self.amounts[..n].iter_mut().enumerate() {
                let value = initial[i] / gas.source().models()[i].molar_mass().value();
                if !value.is_finite() || (initial[i] > 0.0 && value == 0.0) {
                    return Err(ZoneKineticsError::Invalid(
                        "unrepresentable original species amount",
                    ));
                }
                *amount = Amount::new(value);
            }
            let reference = sum(initial
                .iter()
                .zip(gas.energy_references())
                .map(|(m, e)| m * e));
            let physical_energy = state.sensible_internal_energy().value() + reference;
            if !physical_energy.is_finite() {
                return Err(ZoneKineticsError::Invalid(
                    "unrepresentable physical layer energy",
                ));
            }
            let endpoint = self
                .reactor
                .advance_adiabatic_energy(
                    &self.amounts[..n],
                    Energy::new(physical_energy),
                    state.volume(),
                    state.temperature(),
                    duration,
                    &mut cancelled,
                )
                .map_err(|e| match e {
                    IsochoricKineticsError::Cancelled => ZoneKineticsError::Cancelled,
                    source => ZoneKineticsError::Chemistry { layer, source },
                })?;
            let rates = &mut result.chemical_rates[k];
            let mut maximum_roundoff = 0.0_f64;
            for i in 0..n {
                poll(&mut cancelled)?;
                let mass =
                    endpoint.amounts()[i].value() * gas.source().models()[i].molar_mass().value();
                rates.species_kg_s[i] = (mass - initial[i]) / dt;
                let actual = initial[i] + dt * rates.species_kg_s[i];
                let roundoff = (actual - mass).abs();
                let tolerance = controls.relative_conservation_tolerance * mass.max(initial[i]);
                if !mass.is_finite()
                    || !rates.species_kg_s[i].is_finite()
                    || !actual.is_finite()
                    || actual < 0.0
                    || !roundoff.is_finite()
                    || roundoff > tolerance
                {
                    return Err(ZoneKineticsError::Balance(
                        "species rate cannot represent original kinetic endpoint",
                    ));
                }
                maximum_roundoff = maximum_roundoff.max(roundoff);
            }
            rates.referenced_energy_w =
                -sum((0..n).map(|i| rates.species_kg_s[i] * gas.energy_references()[i]));
            let mass_source = sum(rates.species_kg_s[..n].iter().map(|r| dt * r));
            let energy_source = dt * rates.referenced_energy_w
                + sum((0..n).map(|i| dt * rates.species_kg_s[i] * gas.energy_references()[i]));
            if !rates.referenced_energy_w.is_finite()
                || !mass_source.is_finite()
                || mass_source.abs() > controls.relative_conservation_tolerance * state.mass_kg()
                || !energy_source.is_finite()
                || energy_source.abs() > controls.energy_tolerance_j
            {
                return Err(ZoneKineticsError::Balance(
                    "chemical rate mass or physical energy",
                ));
            }
            result.layers[k] = Some(LayerKineticsReport {
                volume: state.volume(),
                physical_internal_energy: Energy::new(physical_energy),
                isochoric_temperature: endpoint.temperature(),
                isochoric_pressure: endpoint.pressure(),
                reactor: endpoint.balance(),
                maximum_species_roundoff_kg: maximum_roundoff,
                mass_source_kg: mass_source,
                physical_energy_source_j: energy_source,
            });
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Apply both chemical increments and original pressure/interface work on
    /// one unpublished candidate. The physical clock keeps its pre-call bits:
    /// the caller's transport or enclosing time owner already counts duration.
    /// The returned rates belong in its additional source account exactly once.
    pub fn apply_split(
        &mut self,
        room: &mut TwoZone,
        duration: Time,
        selected: [bool; 2],
        zone_controls: ZoneSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SplitKineticsReport, ZoneKineticsError> {
        let chemistry = self.evaluate(room, duration, selected, &mut cancelled)?;
        self.candidate.copy_state_from(room)?;
        let initial_energy = room.physical_internal_energy().value();
        let initial_mass = room.layers[0].mass_kg() + room.layers[1].mass_kg();
        let relaxation = self.candidate.advance(
            duration,
            chemistry.chemical_rates,
            zone_controls,
            &mut cancelled,
        )?;
        self.candidate.time = room.time;
        let controls = self.controls();
        for layer in [Layer::Lower, Layer::Upper] {
            let state = self.candidate.layer(layer);
            if selected[layer.index()]
                && state.mass_kg() > 0.0
                && !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
                    .contains(&state.temperature().value())
            {
                return Err(ZoneKineticsError::TemperatureDomain {
                    layer,
                    temperature_k: state.temperature().value(),
                });
            }
        }
        let mass_residual =
            self.candidate.layers[0].mass_kg() + self.candidate.layers[1].mass_kg() - initial_mass;
        let energy_residual = self.candidate.physical_internal_energy().value() - initial_energy;
        if !mass_residual.is_finite()
            || mass_residual.abs() > controls.relative_conservation_tolerance * initial_mass
            || !energy_residual.is_finite()
            || energy_residual.abs()
                > 2.0 * controls.energy_tolerance_j + zone_controls.absolute_energy_tolerance_j
        {
            return Err(ZoneKineticsError::Balance(
                "relaxed whole-room physical mass or energy",
            ));
        }
        poll(&mut cancelled)?;
        room.copy_state_from(&self.candidate)?;
        Ok(SplitKineticsReport {
            chemistry,
            relaxation,
            mass_residual_kg: mass_residual,
            physical_energy_residual_j: energy_residual,
        })
    }
}

#[cfg(test)]
mod tests;
