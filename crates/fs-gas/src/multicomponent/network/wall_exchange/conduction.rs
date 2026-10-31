//! Source-qualified solid resistance on the original gas/wall heat transaction.

use super::{NetworkWallOperator, WallContact, WallExchangeError, WallThermal, poll};
use crate::multicomponent::{allocated, positive};
use fs_material::thermal::transport::PreparedThermalTransport;
use fs_qty::{Area, Length};
use std::sync::Arc;

/// Fixed-geometry conduction between the gas-facing surface and the retained
/// wall temperature. The original material supplies directional k(T).
///
/// `length / area` is the planar resistance geometry. For a finite uniform
/// wall whose state represents its centre, supply the distance from surface to
/// centre, not the complete wall thickness. The existing contact conductance
/// represents ONLY the gas-side film. This adds their resistances in series;
/// it never adds a second wall mass, heat supply, or caloric law.
///
/// Conductivity is sampled at the actual wall temperature at the beginning of
/// each source half and frozen for that half's conservative implicit exchange.
/// The final wall temperature must remain in the original source domain. This
/// is a first-order local constitutive approximation, not an exact nonlinear
/// through-wall profile. Refine the interval and wall discretization when its
/// spatial or temporal temperature variation matters.
///
/// This fixed-context source admits temperature-coordinate contacts only.
/// Latent/phase-changing energy contacts require their separate phase-aware
/// conductivity owner and refuse here; a melting plateau must not reuse a
/// solid conductivity merely because its temperature has not changed.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialWallConduction {
    source: Arc<PreparedThermalTransport>,
    length: Length,
    area: Area,
    length_over_area: f64,
}

impl MaterialWallConduction {
    /// Retain the admitted original material, direction and all source receipts.
    /// Density and heat capacity stay available on the source but are not
    /// converted to additional storage by this resistance binding.
    pub fn new(
        source: Arc<PreparedThermalTransport>,
        length: Length,
        area: Area,
    ) -> Result<Self, WallExchangeError> {
        let length_over_area = length.value() / area.value();
        if !positive(length.value()) || !positive(area.value()) || !positive(length_over_area) {
            return Err(WallExchangeError::Invalid("wall conduction geometry"));
        }
        Ok(Self {
            source,
            length,
            area,
            length_over_area,
        })
    }

    /// Original material-card identity, source lines, fixed context and k(T).
    pub fn source(&self) -> &PreparedThermalTransport {
        &self.source
    }
    /// Distance from gas-facing surface to the retained wall-temperature state.
    pub const fn length(&self) -> Length {
        self.length
    }
    /// Actual exposed conduction area, counted once for this contact.
    pub const fn area(&self) -> Area {
        self.area
    }

    fn conductance(&self, temperature_k: f64) -> Result<f64, WallExchangeError> {
        let sample = self
            .source
            .at_temperature(temperature_k)
            .map_err(WallExchangeError::Conduction)?;
        let conductance = sample.conductivity_w_m_k() / self.length_over_area;
        if !positive(conductance) {
            return Err(WallExchangeError::Invalid(
                "wall conduction is unrepresentable",
            ));
        }
        Ok(conductance)
    }
}

#[derive(Clone)]
pub(super) struct BoundWallConduction {
    layer: MaterialWallConduction,
    // Keep the authored film separate from the combined coefficient stored in
    // accepted contacts, so repeated steps cannot accumulate the same resistance.
    fixed_film_w_k: f64,
}

impl NetworkWallOperator<'_> {
    /// Bind solid k(T) and geometry to any branch/plenum thermal contact.
    ///
    /// The original contact's conductance becomes the gas-side film coefficient.
    /// A bound circular gas-transport model still refreshes that film from its
    /// actual gas composition, T and p. Other contacts retain the initial film;
    /// the resulting series conductance is never mistaken for a new film on the
    /// next half. `None` leaves the original contact path unchanged.
    ///
    /// Finite thermal storage remains exclusively in the original contact and
    /// optional material caloric binding. The gas/wall heat debit is performed
    /// exactly once by that same exchange owner. Bind before evolution; source,
    /// geometry and topology refusals publish no gas or wall state.
    pub fn with_material_conduction(
        mut self,
        branches: Vec<Vec<Option<MaterialWallConduction>>>,
        plenums: Vec<Option<MaterialWallConduction>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if self.last_report.is_some()
            || self.trial_active
            || self.conduction.is_some()
            || branches.len() != self.branches.len()
            || plenums.len() != self.plenum_volumes.len()
            || !branches
                .iter()
                .flatten()
                .chain(&plenums)
                .any(Option::is_some)
        {
            return Err(WallExchangeError::Invalid(
                "material wall conduction state or topology",
            ));
        }
        for (input, shape) in branches.iter().zip(&self.branches) {
            poll(&mut cancelled)?;
            if input.len() != shape.mesh.cell_count() {
                return Err(WallExchangeError::Invalid(
                    "material conduction branch cell count",
                ));
            }
        }
        let mut bindings = allocated(self.contacts.len(), None)?;
        for (index, layer) in branches.into_iter().flatten().chain(plenums).enumerate() {
            poll(&mut cancelled)?;
            if let Some(layer) = layer {
                let (temperature, film) = thermal_state(self.contacts[index].thermal)?;
                layer.conductance(temperature)?;
                bindings[index] = Some(BoundWallConduction {
                    layer,
                    fixed_film_w_k: film,
                });
            }
        }
        poll(&mut cancelled)?;
        self.conduction = Some(bindings);
        Ok(self)
    }

    /// Original conduction source/geometry for a particular branch cell.
    pub fn branch_conduction(&self, branch: usize, cell: usize) -> Option<&MaterialWallConduction> {
        let shape = self.branches.get(branch)?;
        if cell >= shape.mesh.cell_count() {
            return None;
        }
        Some(
            &self
                .conduction
                .as_ref()?
                .get(shape.start + cell)?
                .as_ref()?
                .layer,
        )
    }

    /// Original conduction source/geometry for a compact plenum contact.
    pub fn plenum_conduction(&self, plenum: usize) -> Option<&MaterialWallConduction> {
        if plenum >= self.plenum_volumes.len() {
            return None;
        }
        Some(
            &self
                .conduction
                .as_ref()?
                .get(self.plenum_start + plenum)?
                .as_ref()?
                .layer,
        )
    }

    pub(super) fn film_contact(&self, index: usize, mut contact: WallContact) -> WallContact {
        if let Some(binding) = self.conduction.as_ref().and_then(|all| all[index].as_ref()) {
            set_conductance(&mut contact, binding.fixed_film_w_k);
        }
        contact
    }

    pub(super) fn conduction_contact(
        &self,
        index: usize,
        mut contact: WallContact,
    ) -> Result<WallContact, WallExchangeError> {
        if let Some(binding) = self.conduction.as_ref().and_then(|all| all[index].as_ref()) {
            let (temperature, film) = thermal_state(contact.thermal)?;
            let solid = binding.layer.conductance(temperature)?;
            // Stable harmonic combination without overflowing Gfilm * Gsolid
            // or taking the reciprocal of an exactly disabled gas film.
            let smaller = film.min(solid);
            let combined = smaller / (1.0 + smaller / film.max(solid));
            if !combined.is_finite() || combined < 0.0 || film > 0.0 && combined == 0.0 {
                return Err(WallExchangeError::Invalid(
                    "series wall conductance is unrepresentable",
                ));
            }
            set_conductance(&mut contact, combined);
        }
        Ok(contact)
    }

    pub(super) fn admit_conduction_endpoint(
        &self,
        index: usize,
        contact: WallContact,
    ) -> Result<(), WallExchangeError> {
        if let Some(binding) = self.conduction.as_ref().and_then(|all| all[index].as_ref()) {
            let (temperature, _) = thermal_state(contact.thermal)?;
            binding.layer.conductance(temperature)?;
        }
        Ok(())
    }
}

fn thermal_state(thermal: WallThermal) -> Result<(f64, f64), WallExchangeError> {
    match thermal {
        WallThermal::Lumped {
            temperature_k,
            conductance_w_per_k,
            ..
        }
        | WallThermal::Reservoir {
            temperature_k,
            conductance_w_per_k,
        } => {
            if !positive(temperature_k)
                || !conductance_w_per_k.is_finite()
                || conductance_w_per_k < 0.0
            {
                return Err(WallExchangeError::Invalid("wall conduction contact state"));
            }
            Ok((temperature_k, conductance_w_per_k))
        }
        WallThermal::Adiabatic => Err(WallExchangeError::Invalid(
            "adiabatic wall has no conduction temperature",
        )),
        WallThermal::Caloric { .. } => Err(WallExchangeError::Invalid(
            "energy-coordinate wall requires phase-aware conductivity",
        )),
    }
}

fn set_conductance(contact: &mut WallContact, conductance: f64) {
    match &mut contact.thermal {
        WallThermal::Lumped {
            conductance_w_per_k,
            ..
        }
        | WallThermal::Reservoir {
            conductance_w_per_k,
            ..
        } => *conductance_w_per_k = conductance,
        // Only admitted temperature-coordinate contacts have a binding.
        WallThermal::Adiabatic | WallThermal::Caloric { .. } => {}
    }
}
