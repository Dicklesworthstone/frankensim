//! Temperature-dependent adhesion on the actual deposited-cell thermal state.
//!
//! The original cohesive owner changes coefficients at fixed jump and retained
//! damage, then relaxes onto the new envelope. Thermal failure updates the SAME
//! contact closure used by subsequent heat transport. Cooling never rebonds a
//! failed patch. Parameter work and damage release remain separate mechanical
//! energy accounts: neither is silently added to cell enthalpy.
//!
//! This is an explicit staggered endpoint update. Complete the preceding
//! mechanical segment before changing its law; otherwise a later timed replay
//! would incorrectly apply the new coefficients to the entire earlier interval.
//! The supplied weighted cell temperature is a declared interface approximation,
//! not a resolved through-thickness temperature or a calibrated cure model.

use super::{ContactReport, Contacts, FractureCooling, HeatSum, Interface, positive};
use crate::transient::deposition::{ContactUpdate, DepositionError};
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_material::state_point::MaterialStatePointError;
use std::collections::BTreeMap;

/// Original source and its explicitly resolved thermal successor for one patch.
/// Use `interface_temperature_k` and the original source's `at_temperature` to
/// prepare the successor. One update is required for every retained interface.
#[derive(Clone, Copy, Debug)]
pub struct TemperatureUpdate<'old, 'new> {
    pub interface: u64,
    pub previous: &'old ResolvedBilinearCohesive,
    pub successor: &'new ResolvedBilinearCohesive,
    /// T_interface = weight_a T_a + (1 - weight_a) T_b, in [0, 1].
    pub temperature_weight_a: f64,
}

/// Independent constitutive and area-integrated energy acceptance gates.
#[derive(Clone, Copy, Debug)]
pub struct TemperatureConfig {
    pub local_energy_tolerance_j_m2: f64,
    pub total_energy_tolerance_j: f64,
}

/// Work and loss at fixed mechanical displacement, integrated on original area.
#[derive(Clone, Copy, Debug)]
pub struct TemperatureReport {
    /// New damage release and any thermal-contact opening/conductance changes.
    /// Viscous dissipation is zero: no displacement or time passes here.
    pub contact: ContactReport,
    /// Signed work of changing coefficients at fixed jump and retained damage.
    pub parameter_work_j: f64,
    pub stored_change_j: f64,
    /// Parameter work minus storage change minus new damage dissipation.
    pub defect_j: f64,
}

#[derive(Debug, Clone)]
pub enum Error {
    Contact(super::Error),
    Source {
        interface: u64,
        error: MaterialStatePointError,
    },
    Invalid(&'static str),
    Temperature {
        interface: u64,
        expected_k: f64,
        supplied_k: Option<f64>,
    },
    Energy {
        defect_j: f64,
        tolerance_j: f64,
    },
}
impl From<super::Error> for Error {
    fn from(error: super::Error) -> Self {
        Self::Contact(error)
    }
}
impl From<DepositionError> for Error {
    fn from(error: DepositionError) -> Self {
        Self::Contact(error.into())
    }
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal bond transfer refused: {self:?}")
    }
}
impl std::error::Error for Error {}

impl<'t> FractureCooling<'t, '_> {
    /// Evaluate the declared interface temperature from the current caloric
    /// cells, including any latent-heat plateau owned by their original law.
    /// The same arithmetic is used to verify a proposed source successor.
    pub fn interface_temperature_k(&self, interface: u64, weight_a: f64) -> Result<f64, Error> {
        if !weight_a.is_finite() || !(0.0..=1.0).contains(&weight_a) {
            return Err(Error::Invalid(
                "interface temperature weight must be in [0,1]",
            ));
        }
        let contact = self
            .contacts
            .0
            .get(&interface)
            .ok_or(Error::Invalid("unknown cohesive interface"))?;
        let a = self.thermal.temperature_k(contact.reference.a)?;
        let b = self.thermal.temperature_k(contact.reference.b)?;
        let temperature = weight_a * a + (1.0 - weight_a) * b;
        if !positive(temperature) {
            return Err(Error::Invalid("unrepresentable interface temperature"));
        }
        Ok(temperature)
    }

    /// Change the original source-defined bond at the CURRENT synchronized
    /// mechanical/thermal endpoint. Geometry, area, jump, damage history, source
    /// claims, nonthermal coordinates and both clocks remain bound.
    ///
    /// The original cohesive temperature law supplies every point and its work;
    /// source temperatures must equal the declared actual cell average. Failure
    /// may open a patch or select its explicitly declared failed-contact thermal
    /// conductance immediately. Existing open/failed patches may not heal.
    ///
    /// Neither parameter work nor damage release becomes heat automatically.
    /// A thermodynamic coupling must supply the former and partition the latter
    /// explicitly. This operation does not infer events inside the last thermal
    /// step, thermal expansion, mechanical equilibrium or a heat-source closure.
    /// Every refusal, including cancellation, leaves the accepted input intact.
    pub fn transfer_temperature<'new>(
        &self,
        time_s: f64,
        updates: &[TemperatureUpdate<'_, 'new>],
        temperature_axis: &str,
        config: TemperatureConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(FractureCooling<'t, 'new>, TemperatureReport), Error> {
        checkpoint()?;
        if time_s != self.thermal.time_s() || time_s != self.last_synchronized_time_s {
            return Err(Error::Invalid(
                "complete the mechanical segment at the current thermal time before changing its law",
            ));
        }
        if updates.len() != self.contacts.0.len()
            || !positive(config.local_energy_tolerance_j_m2)
            || !positive(config.total_energy_tolerance_j)
        {
            return Err(Error::Invalid(
                "complete temperature updates and positive energy gates required",
            ));
        }
        let mut ordered = BTreeMap::new();
        for update in updates {
            checkpoint()?;
            if !self.contacts.0.contains_key(&update.interface)
                || ordered.insert(update.interface, update).is_some()
            {
                return Err(Error::Invalid("unknown or duplicate temperature interface"));
            }
        }
        let mut contacts = BTreeMap::new();
        let mut replacements = Vec::with_capacity(updates.len());
        let mut parameter = HeatSum::default();
        let mut storage = HeatSum::default();
        let mut damage = HeatSum::default();
        let mut contact = ContactReport::default();
        for (id, update) in ordered {
            checkpoint()?;
            let old = self.contacts.0[&id];
            if !core::ptr::eq(old.point.law(), update.previous.law()) {
                return Err(Error::Invalid(
                    "previous temperature source does not own the accepted bond",
                ));
            }
            update
                .previous
                .require_temperature_successor(update.successor, temperature_axis)
                .map_err(|error| Error::Source {
                    interface: id,
                    error,
                })?;
            let expected_k = self.interface_temperature_k(id, update.temperature_weight_a)?;
            let supplied_k = update
                .successor
                .resolved()
                .state_coordinate(temperature_axis);
            if supplied_k != Some(expected_k) {
                return Err(Error::Temperature {
                    interface: id,
                    expected_k,
                    supplied_k,
                });
            }
            let step = old
                .point
                .transfer_temperature(update.successor.law(), config.local_energy_tolerance_j_m2)
                .map_err(|error| super::Error::Material {
                    interface: id,
                    error,
                })?;
            let closed = old.next_closed(step.point)?;
            let next = Interface {
                reference: old.reference,
                point: step.point,
                closed,
                policy: old.policy,
            };
            let expected = old.thermal_link();
            let replacement = next.thermal_link();
            if old.closed && !closed {
                contact.opened += 1;
            }
            if !old.closed && closed {
                contact.reclosed += 1;
            }
            if old.closed && closed && expected != replacement {
                contact.changed_conductances += 1;
            }
            let area = old.reference.area_m2;
            parameter.add(area * step.parameter_work_j_m2)?;
            storage.add(area * step.stored_change_j_m2)?;
            damage.add(area * step.damage_dissipated_j_m2)?;
            replacements.push(ContactUpdate {
                id,
                expected,
                replacement,
            });
            contacts.insert(id, next);
        }
        let parameter_work_j = parameter.value()?;
        let stored_change_j = storage.value()?;
        contact.fracture_dissipation_j = damage.value()?;
        let mut defect = HeatSum::default();
        defect.add(parameter_work_j)?;
        defect.add(-stored_change_j)?;
        defect.add(-contact.fracture_dissipation_j)?;
        let defect_j = defect.value()?;
        if defect_j.abs() > config.total_energy_tolerance_j {
            return Err(Error::Energy {
                defect_j,
                tolerance_j: config.total_energy_tolerance_j,
            });
        }
        let thermal = self
            .thermal
            .replace_contacts(&replacements, &mut checkpoint)?;
        checkpoint()?;
        Ok((
            FractureCooling {
                thermal,
                contacts: Contacts(contacts),
                last_synchronized_time_s: self.last_synchronized_time_s,
            },
            TemperatureReport {
                contact,
                parameter_work_j,
                stored_change_j,
                defect_j,
            },
        ))
    }
}

#[cfg(test)]
mod tests;
