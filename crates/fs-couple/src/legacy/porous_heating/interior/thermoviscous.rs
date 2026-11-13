//! Resolved viscous acoustic heating on the original thermal timeline.
//!
//! The full mixed acoustic owner independently integrates shear and bulk
//! viscous heating on its original curved or affine tetrahedra. Those two
//! positive forms supply the native caloric cells. Conductive entropy production
//! and molecular/caloric availability destruction remain separate observations:
//! neither is a closed mean caloric source. No signed port work, pressure-work
//! difference, first-order heat excitation or numerical residual enters heat.
//!
//! This is a held, cycle-mean source, with peak phasors in `exp(+i omega t)`.
//! The caller owns the spatial and caloric correspondence and must justify
//! fixed acoustic coefficients and amplitudes over the thermal advance. No
//! acoustic startup, streaming, thermal expansion or material refresh is inferred.

use super::{
    CellBinding, InteriorHeatCell, InteriorHeatingError, MAX_INTERIOR_HEATING_CELLS, finite,
    same_measure,
};
use fs_feec::acoustics::thermoviscous::{
    PreparedThermoviscous, ThermoviscousElementPower, ThermoviscousError,
};
use fs_math::{c64::C64, dd::Dd};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig,
};
use std::collections::BTreeMap;

/// One original tetrahedron, its native physical powers and caloric receiver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousCellHeating {
    pub cell: CellId,
    pub element: usize,
    /// Already integrated, cycle-mean powers; no second volume or half factor.
    pub power: ThermoviscousElementPower,
}

/// Native acoustic, thermal, or original receiver-binding refusal.
#[derive(Debug)]
pub enum ThermoviscousHeatingError {
    Acoustic(ThermoviscousError),
    Thermal(DepositionError),
    Binding(InteriorHeatingError),
}

impl core::fmt::Display for ThermoviscousHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved viscous heating refused: {self:?}")
    }
}
impl std::error::Error for ThermoviscousHeatingError {}
impl From<ThermoviscousError> for ThermoviscousHeatingError {
    fn from(value: ThermoviscousError) -> Self {
        Self::Acoustic(value)
    }
}
impl From<DepositionError> for ThermoviscousHeatingError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<InteriorHeatingError> for ThermoviscousHeatingError {
    fn from(value: InteriorHeatingError) -> Self {
        Self::Binding(value)
    }
}

fn invalid(message: &'static str) -> ThermoviscousHeatingError {
    InteriorHeatingError::Invalid(message).into()
}

struct BoundThermoviscousCell<'t> {
    heat: ThermoviscousCellHeating,
    binding: CellBinding<'t>,
}

/// Frozen full-field sources retaining the original acoustic preparation.
///
/// Source/material quadrature is evaluated once during construction. Native
/// thermal substeps reuse cached watts; observer clones copy only bounded work
/// accounts and retain this same immutable owner and original field borrows.
pub struct ThermoviscousViscousHeating<'p, 'a, 'c, 't> {
    prepared: &'p PreparedThermoviscous<'a, 'c>,
    velocity_m_s: &'p [[C64; 3]],
    temperature_k: &'p [C64],
    cells: Vec<BoundThermoviscousCell<'t>>,
    sources_w: BTreeMap<CellId, f64>,
}

/// Accumulation over accepted native thermal durations only.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousHeatingEnergy {
    pub shear_heating_j: f64,
    pub bulk_heating_j: f64,
    /// Conductive entropy production [J/K], observed at the frozen source state.
    pub conductive_entropy_j_k: f64,
    /// Acoustic availability destruction [J]; not included in deposited_j.
    pub thermal_availability_loss_j: f64,
    /// Molecular/caloric availability destruction [J]; not deposited as heat.
    pub caloric_availability_loss_j: f64,
    /// Shear plus bulk heat, included once in native external-heat accounting.
    pub deposited_j: f64,
}

/// The accepted native thermal candidate and its separate acoustic account.
pub struct ThermoviscousHeatingStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    /// All original authored thermal exchanges plus resolved viscous heat.
    pub heat_report: EnergyAudit,
    pub acoustic: ThermoviscousHeatingEnergy,
    /// Viscous heat only, in the original tetrahedron/destination order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

impl<'p, 'a, 'c, 't> ThermoviscousViscousHeating<'p, 'a, 'c, 't> {
    /// Bind every original fluid tetrahedron to one distinct native heat cell.
    ///
    /// Full P2 velocity and temperature coefficients use the native hierarchical
    /// basis. They must belong to this exact frequency/material/geometry owner;
    /// as in its loss observation, finite supplied fields are not recertified
    /// solutions of the acoustic equations. Moving-wall lifts must already be
    /// present in the resolved velocity. Structural and solid-wall losses are
    /// outside this fluid-volume source.
    ///
    /// Native volume and mass must match each original acoustic element to
    /// roundoff. The original caloric object, storage basis, material identity
    /// and body/segment bindings are retained. Geometry or mass equality alone
    /// does not establish the caller's spatial/material correspondence.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        prepared: &'p PreparedThermoviscous<'a, 'c>,
        velocity_m_s: &'p [[C64; 3]],
        temperature_k: &'p [C64],
        thermal: &DepositionSimulation<'_>,
        destinations: &[InteriorHeatCell<'t>],
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermoviscousHeatingError> {
        checkpoint()?;
        let count = prepared.element_count();
        if maximum_cells == 0 || maximum_cells > MAX_INTERIOR_HEATING_CELLS || count > maximum_cells
        {
            return Err(DepositionError::Budget("resolved viscous heating cells").into());
        }
        if count == 0 || destinations.len() != count {
            return Err(invalid(
                "resolved viscous heat needs every original tetrahedron exactly once",
            ));
        }
        let mut bindings = Vec::with_capacity(count);
        let mut sources_w = BTreeMap::new();
        for route in destinations {
            checkpoint()?;
            if sources_w.insert(route.cell, 0.0).is_some() {
                return Err(invalid(
                    "resolved viscous heat destinations must be distinct",
                ));
            }
            bindings.push(CellBinding::try_new(thermal, route.cell, route.caloric)?);
        }
        let mut refusal = None;
        let powers =
            prepared.element_dissipation(velocity_m_s, temperature_k, &mut || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    refusal = Some(error);
                    true
                }
            });
        if let Some(error) = refusal {
            return Err(error.into());
        }
        let mut cells = Vec::with_capacity(count);
        for (element, (binding, power)) in bindings.into_iter().zip(powers?).enumerate() {
            checkpoint()?;
            if !same_measure(binding.geometry.volume_m3, power.volume_m3)
                || !same_measure(binding.mass_kg, power.mass_kg)
            {
                return Err(invalid(
                    "resolved viscous heat needs the original fluid volume and mass",
                ));
            }
            for value in [
                power.shear_heating_w,
                power.bulk_heating_w,
                power.conductive_entropy_w_k,
                power.thermal_availability_loss_w,
                power.caloric_availability_loss_w,
            ] {
                if !value.is_finite() || value < 0.0 {
                    return Err(invalid(
                        "resolved physical loss forms must be finite and passive",
                    ));
                }
            }
            let watts = finite(power.viscous_heat_w())?;
            sources_w.insert(binding.cell, watts);
            cells.push(BoundThermoviscousCell {
                heat: ThermoviscousCellHeating {
                    cell: binding.cell,
                    element,
                    power,
                },
                binding,
            });
        }
        checkpoint()?;
        Ok(Self {
            prepared,
            velocity_m_s,
            temperature_k,
            cells,
            sources_w,
        })
    }

    #[must_use]
    pub const fn prepared(&self) -> &'p PreparedThermoviscous<'a, 'c> {
        self.prepared
    }

    #[must_use]
    pub const fn velocity_m_s(&self) -> &'p [[C64; 3]] {
        self.velocity_m_s
    }

    #[must_use]
    pub const fn temperature_k(&self) -> &'p [C64] {
        self.temperature_k
    }

    pub fn cells(&self) -> impl ExactSizeIterator<Item = &ThermoviscousCellHeating> {
        self.cells.iter().map(|cell| &cell.heat)
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermoviscousHeatingError> {
        for cell in &self.cells {
            checkpoint()?;
            cell.binding.validate(thermal)?;
        }
        Ok(())
    }

    /// Start a native composable observer with a zero accepted-work account.
    /// Retain the returned accepted observer to continue an existing account.
    pub fn observer(&self) -> ThermoviscousHeatingObserver<'_, 'p, 'a, 'c, 't> {
        ThermoviscousHeatingObserver {
            owner: self,
            account: ThermoviscousHeatAccount {
                cells: vec![Dd::ZERO; self.cells.len()],
                mechanisms: [Dd::ZERO; 5],
            },
        }
    }

    /// Advance with held native viscous heat. Original links, reservoirs,
    /// authored sources, events, caloric domains and thermal timestep admission
    /// remain in DepositionSimulation. Do not also author these watts as thermal
    /// sources. A refusal publishes neither candidate state nor new heat debit.
    pub fn advance_to<'s>(
        &self,
        initial: &DepositionSimulation<'s>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThermoviscousHeatingStep<'s>, ThermoviscousHeatingError> {
        self.validate(initial, &mut checkpoint)?;
        let (thermal, observer, heat_report) =
            initial.advance_observed(&self.observer(), target_s, config, &mut checkpoint)?;
        self.validate(&thermal, &mut checkpoint)?;
        let acoustic = observer.energy()?;
        let cell_energy_j = observer.cell_energy_j()?;
        checkpoint()?;
        Ok(ThermoviscousHeatingStep {
            thermal,
            heat_report,
            acoustic,
            cell_energy_j,
        })
    }
}

#[derive(Clone)]
struct ThermoviscousHeatAccount {
    cells: Vec<Dd>,
    mechanisms: [Dd; 5],
}

/// Native observer retaining the immutable source and only its accepted work.
#[derive(Clone)]
pub struct ThermoviscousHeatingObserver<'h, 'p, 'a, 'c, 't> {
    owner: &'h ThermoviscousViscousHeating<'p, 'a, 'c, 't>,
    account: ThermoviscousHeatAccount,
}

impl ThermoviscousHeatingObserver<'_, '_, '_, '_, '_> {
    pub fn energy(&self) -> Result<ThermoviscousHeatingEnergy, ThermoviscousHeatingError> {
        let deposited = self
            .account
            .cells
            .iter()
            .copied()
            .fold(Dd::ZERO, |a, b| a + b);
        Ok(ThermoviscousHeatingEnergy {
            shear_heating_j: finite(self.account.mechanisms[0].to_f64())?,
            bulk_heating_j: finite(self.account.mechanisms[1].to_f64())?,
            conductive_entropy_j_k: finite(self.account.mechanisms[2].to_f64())?,
            thermal_availability_loss_j: finite(self.account.mechanisms[3].to_f64())?,
            caloric_availability_loss_j: finite(self.account.mechanisms[4].to_f64())?,
            deposited_j: finite(deposited.to_f64())?,
        })
    }

    pub fn cell_energy_j(&self) -> Result<Vec<(CellId, f64)>, ThermoviscousHeatingError> {
        self.owner
            .cells
            .iter()
            .zip(&self.account.cells)
            .map(|(cell, energy)| Ok((cell.heat.cell, finite(energy.to_f64())?)))
            .collect()
    }
}

impl DepositionObserver for ThermoviscousHeatingObserver<'_, '_, '_, '_, '_> {
    type Error = ThermoviscousHeatingError;

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        Ok(self.owner.sources_w.clone())
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.owner.validate(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        for (cell, energy) in self.owner.cells.iter().zip(&mut self.account.cells) {
            checkpoint()?;
            let power = cell.heat.power;
            *energy = *energy + Dd::from_f64(finite(power.viscous_heat_w() * dt)?);
            for (sum, watts) in self.account.mechanisms.iter_mut().zip([
                power.shear_heating_w,
                power.bulk_heating_w,
                power.conductive_entropy_w_k,
                power.thermal_availability_loss_w,
                power.caloric_availability_loss_w,
            ]) {
                *sum = *sum + Dd::from_f64(finite(watts * dt)?);
            }
        }
        Ok(())
    }
}
