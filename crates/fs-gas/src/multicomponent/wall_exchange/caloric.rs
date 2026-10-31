//! Original material storage in the existing gas/wall scalar exchange.

use super::*;
#[cfg(feature = "material")]
use fs_material::thermal::transport::storage::PreparedHeatStorage;
use fs_thermal::enthalpy::{DepositedCell, ThermalLaw, ThermalStorageBasis};
use std::sync::Arc;

mod thermal_law;
use thermal_law::EnergyStorage;

/// Required physical interpretation of the source capacity. No Cp-to-Cv
/// conversion or unrepresented pressure/deformation work is inferred.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WallCaloricConstraint {
    /// The supplied capacity is dU/dT under the fixed wall-volume constraint.
    IsochoricInternalEnergy,
    /// Stationary solid approximation: thermal expansion, pressure-volume work
    /// and the Cp-Cv difference are explicitly neglected over this source span.
    StationarySolidNegligibleExpansion,
    /// Specific enthalpy from an original thermal chart. The source owns its
    /// pressure-work convention; this storage is not rigid-vessel U.
    SourceEnthalpy,
}

/// Immutable source and fixed constituent amounts, shared without copying tables.
/// Temperature or specific h/u remains exclusively in the original contact.
/// A bound contact's capacity field is a refreshed tangent cache; its storage
/// change is ALWAYS the source integral, never that cache times delta T.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialWallCaloric<'law> {
    storage: Storage<'law>,
    constraint: WallCaloricConstraint,
    #[cfg(feature = "material")]
    phase_curve: Option<Arc<fs_material::phase::EquilibriumEnthalpyPhaseCurve>>,
}

#[derive(Clone, Debug, PartialEq)]
enum Storage<'law> {
    #[cfg(feature = "material")]
    Capacity(Arc<PreparedHeatStorage>),
    Energy(EnergyStorage<'law>),
}

struct StorageUpdate {
    contact: WallContact,
    temperature_k: f64,
    energy_change_j: f64,
    temperature_slope_k_per_j: Option<f64>,
}

impl<'law> MaterialWallCaloric<'law> {
    /// Bind an admitted source with an explicit compatible work constraint.
    /// This retains fixed amounts and source context; it adds no phase change,
    /// wall deformation, reaction, spatial wall conduction or external work.
    #[cfg(feature = "material")]
    pub fn new(storage: Arc<PreparedHeatStorage>, constraint: WallCaloricConstraint) -> Self {
        Self {
            storage: Storage::Capacity(storage),
            constraint,
            phase_curve: None,
        }
    }
    /// Retain the original immutable T(h/u), including a latent plateau.
    /// The cell's identity, mass and reference geometry are admitted by that
    /// law. Its deposition boundary is not applied: this gas contact owns heat.
    /// The caller declares a continuous supported specific-energy interval.
    pub fn from_thermal_law<L: ThermalLaw + Send + Sync + 'law>(
        law: L,
        initial: DepositedCell,
        energy_range_j_kg: [f64; 2],
    ) -> Result<Self, WallExchangeError> {
        let constraint = match law.storage_basis() {
            ThermalStorageBasis::InternalEnergy => WallCaloricConstraint::IsochoricInternalEnergy,
            ThermalStorageBasis::Enthalpy => WallCaloricConstraint::SourceEnthalpy,
        };
        Ok(Self {
            storage: Storage::Energy(EnergyStorage::new(law, initial, energy_range_j_kg)?),
            constraint,
            #[cfg(feature = "material")]
            phase_curve: None,
        })
    }

    /// Bind the original phase chart simultaneously as the caloric law and as
    /// the material phase identity used by wall conduction. Both consumers
    /// share one immutable Arc; temperature alone never identifies latent h.
    ///
    /// This stationary, fixed-volume reduction requires the chart's constant
    /// density to agree with the reference mass/volume. Expansion, shrinkage,
    /// motion and phase-front geometry require their separate physical owners.
    #[cfg(feature = "material")]
    pub fn from_phase_curve(
        curve: Arc<fs_material::phase::EquilibriumEnthalpyPhaseCurve>,
        initial: DepositedCell,
        energy_range_j_kg: [f64; 2],
    ) -> Result<Self, WallExchangeError> {
        let density = initial.mass_kg / initial.geometry.volume_m3;
        if !positive(density)
            || curve.knots().iter().any(|knot| {
                (knot.bulk_density_kg_m3 - density).abs()
                    > 64.0 * f64::EPSILON * density.max(knot.bulk_density_kg_m3)
            })
        {
            return Err(WallExchangeError::Invalid(
                "phase wall requires the chart's constant reference mass density",
            ));
        }
        Ok(Self {
            storage: Storage::Energy(EnergyStorage::from_shared(
                curve.clone(), initial, energy_range_j_kg,
            )?),
            constraint: WallCaloricConstraint::SourceEnthalpy,
            phase_curve: Some(curve),
        })
    }

    /// The exact chart also held by the native h-storage law, when phase-aware
    /// material wall conduction was explicitly admitted.
    #[cfg(feature = "material")]
    pub fn phase_curve(&self) -> Option<&Arc<fs_material::phase::EquilibriumEnthalpyPhaseCurve>> {
        self.phase_curve.as_ref()
    }
    /// Original fixed-amount capacity source, when that representation is bound.
    #[cfg(feature = "material")]
    pub fn storage(&self) -> Option<&PreparedHeatStorage> {
        match &self.storage {
            Storage::Capacity(source) => Some(source),
            Storage::Energy(_) => None,
        }
    }
    /// Original thermal law for an energy-coordinate binding.
    pub fn thermal_law(&self) -> Option<&dyn ThermalLaw> {
        match &self.storage {
            Storage::Energy(source) => Some(source.law.as_ref()),
            #[cfg(feature = "material")]
            Storage::Capacity(_) => None,
        }
    }
    /// Original source identities, mass, reference geometry and initial energy.
    pub fn initial_cell(&self) -> Option<&DepositedCell> {
        match &self.storage {
            Storage::Energy(source) => Some(&source.initial),
            #[cfg(feature = "material")]
            Storage::Capacity(_) => None,
        }
    }
    /// The selected source's physical storage coordinate.
    pub fn storage_basis(&self) -> ThermalStorageBasis {
        match self.constraint {
            WallCaloricConstraint::IsochoricInternalEnergy => ThermalStorageBasis::InternalEnergy,
            WallCaloricConstraint::StationarySolidNegligibleExpansion
            | WallCaloricConstraint::SourceEnthalpy => ThermalStorageBasis::Enthalpy,
        }
    }
    pub const fn constraint(&self) -> WallCaloricConstraint {
        self.constraint
    }

    /// Construct the original finite-wall state from material data alone.
    pub fn contact(
        &self,
        temperature_k: f64,
        conductance_w_per_k: f64,
        drag: WallDrag,
    ) -> Result<WallContact, WallExchangeError> {
        #[cfg(feature = "material")]
        if let Storage::Capacity(storage) = &self.storage {
            let contact = WallContact {
                drag,
                thermal: WallThermal::Lumped {
                    temperature_k,
                    conductance_w_per_k,
                    heat_capacity_j_per_k: storage.heat_capacity_j_k(temperature_k)?,
                },
            };
            contact.validate()?;
            return Ok(contact);
        }
        let _ = (temperature_k, conductance_w_per_k, drag);
        Err(WallExchangeError::Invalid(
            "thermal-law wall must be initialized by energy",
        ))
    }

    /// Initialize h/u from the original reference cell. A temperature alone
    /// cannot specify the inventory inside a latent plateau.
    pub fn initial_contact(
        &self,
        conductance_w_per_k: f64,
        drag: WallDrag,
    ) -> Result<WallContact, WallExchangeError> {
        match &self.storage {
            Storage::Energy(source) => {
                source.contact(source.initial.specific_h_j_kg, conductance_w_per_k, drag)
            }
            #[cfg(feature = "material")]
            Storage::Capacity(_) => Err(WallExchangeError::Invalid(
                "capacity wall requires its initial temperature",
            )),
        }
    }

    /// Resolve the original source at one accepted contact state.
    pub fn temperature_k(&self, state: WallThermal) -> Result<f64, WallExchangeError> {
        let contact = self.bind_contact(WallContact {
            drag: WallDrag::None,
            thermal: state,
        })?;
        match &self.storage {
            Storage::Energy(source) => {
                let (energy, _) = source.state(contact)?;
                source.at(energy)
            }
            #[cfg(feature = "material")]
            Storage::Capacity(_) => match contact.thermal {
                WallThermal::Lumped { temperature_k, .. } => Ok(temperature_k),
                _ => unreachable!("binding admitted temperature storage"),
            },
        }
    }

    pub(crate) fn bind_contact(
        &self,
        contact: WallContact,
    ) -> Result<WallContact, WallExchangeError> {
        match &self.storage {
            Storage::Energy(source) => {
                source.state(contact)?;
                Ok(contact)
            }
            #[cfg(feature = "material")]
            Storage::Capacity(_) => {
                let WallThermal::Lumped {
                    temperature_k,
                    conductance_w_per_k,
                    ..
                } = contact.thermal
                else {
                    return Err(WallExchangeError::Invalid(
                        "material storage requires a finite wall",
                    ));
                };
                self.contact(temperature_k, conductance_w_per_k, contact.drag)
            }
        }
    }

    /// Signed heat into gas that reaches each end of the wall source domain.
    fn heat_limits(
        &self,
        contact: WallContact,
        gas_temperature: f64,
    ) -> Result<[f64; 2], WallExchangeError> {
        match &self.storage {
            Storage::Energy(source) => {
                let _ = gas_temperature;
                let (energy, _) = source.state(contact)?;
                Ok([
                    source.initial.mass_kg * (energy - source.range_j_kg[1]),
                    source.initial.mass_kg * (energy - source.range_j_kg[0]),
                ])
            }
            #[cfg(feature = "material")]
            Storage::Capacity(source) => {
                let temperature = self.temperature_k(contact.thermal)?;
                let [low, high] = source.temperature_span_k();
                let limit =
                    -source.energy_change_j(temperature, gas_temperature.clamp(low, high))?;
                Ok([limit.min(0.0), limit.max(0.0)])
            }
        }
    }

    fn after(&self, contact: WallContact, heat: f64) -> Result<StorageUpdate, WallExchangeError> {
        match &self.storage {
            Storage::Energy(source) => source.after(contact, heat),
            #[cfg(feature = "material")]
            Storage::Capacity(source) => {
                let WallThermal::Lumped {
                    temperature_k,
                    conductance_w_per_k,
                    ..
                } = contact.thermal
                else {
                    return Err(WallExchangeError::Invalid(
                        "capacity wall requires temperature storage",
                    ));
                };
                let next = source.temperature_after_energy_change(temperature_k, heat)?;
                Ok(StorageUpdate {
                    contact: self.contact(next, conductance_w_per_k, contact.drag)?,
                    temperature_k: next,
                    energy_change_j: source.energy_change_j(temperature_k, next)?,
                    temperature_slope_k_per_j: Some(1.0 / source.heat_capacity_j_k(next)?),
                })
            }
        }
    }
}

/// Source-backed nonlinear wall storage with the SAME gas calorics, radial
/// drag and safeguarded scalar owner as `exchange_cell`.
///
/// Solve Q = dt G (Tw(Q)-Tg(Q)), where each original energy owner supplies its
/// temperature and an exact tangent when available. The bracket intersects
/// both source domains and the passive temperature interval. Actual changes and
/// endpoint heat flux must pass the original budget before anything returns.
/// Heat integration is backward Euler; coefficients and geometry are held for
/// this source interval. The gas retains drag heating and all original species.
/// Original h/u walls retain their chart's pressure-work convention in the
/// storage balance; enthalpy is not relabelled as rigid-vessel internal energy.
#[allow(clippy::too_many_arguments)]
pub fn exchange_cell_with_heat_storage(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    wall: &MaterialWallCaloric<'_>,
    control: WallControl,
    cancelled: impl FnMut() -> bool,
) -> Result<WallUpdate, WallExchangeError> {
    exchange_cell_with_heat_storage_and_conductance(
        gas, state, volume_m3, duration, contact, wall, None, control, cancelled,
    )
}

/// The same finite-storage owner with an optional source-owned endpoint
/// conductance. Its callback reads a proposed wall state, whose conductance
/// field always retains the original film. No candidate is published here.
/// The bracket uses passive storage limits, not a frozen-flux bound that could
/// exclude roots when conductivity increases during heating or melting.
#[allow(clippy::too_many_arguments)]
pub(crate) fn exchange_cell_with_heat_storage_and_conductance(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    wall: &MaterialWallCaloric<'_>,
    mut conductance: Option<&mut dyn FnMut(WallContact) -> Result<f64, WallExchangeError>>,
    control: WallControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WallUpdate, WallExchangeError> {
    poll(&mut cancelled)?;
    let contact = wall.bind_contact(contact)?;
    let initial_wall = wall.temperature_k(contact.thermal)?;
    let authored_conductance = match contact.thermal {
        WallThermal::Lumped {
            conductance_w_per_k,
            ..
        }
        | WallThermal::Caloric {
            conductance_w_per_k,
            ..
        } => conductance_w_per_k,
        _ => unreachable!("binding admitted finite storage"),
    };
    let state_dependent = conductance.is_some();
    let mut conductance_at = |candidate: WallContact| {
        let value = match &mut conductance {
            Some(resolve) => resolve(candidate)?,
            None => authored_conductance,
        };
        if !nonnegative(value) {
            return Err(WallExchangeError::Invalid("material wall endpoint conductance"));
        }
        Ok(value)
    };
    let conductance_w_per_k = conductance_at(contact)?;
    // Reuse the original radial loss and gas sensible-heating owner exactly.
    let mut update = exchange_cell(
        gas,
        state,
        volume_m3,
        duration,
        WallContact {
            drag: contact.drag,
            thermal: WallThermal::Adiabatic,
        },
        control,
        &mut cancelled,
    )?;
    update.contact = with_conductance(contact, conductance_w_per_k);
    let gas_at = |energy, cancelled: &mut dyn FnMut() -> bool| {
        gas.primitive(
            Cell {
                partial_densities: state.partial_densities,
                momentum: update.momentum,
                energy,
            },
            cancelled,
        )
    };
    let initial_gas = gas_at(update.energy, &mut cancelled)?;
    let integrated = conductance_w_per_k * duration.value();
    if !nonnegative(integrated) {
        return Err(WallExchangeError::Invalid(
            "integrated material wall conductance",
        ));
    }
    if integrated == 0.0 || initial_wall == initial_gas.temperature {
        poll(&mut cancelled)?;
        return Ok(update);
    }
    let initial_flux = integrated * (initial_wall - initial_gas.temperature);
    if !initial_flux.is_finite() {
        return Err(WallExchangeError::Invalid("material wall heat increment"));
    }
    if initial_flux.abs() <= control.tolerance(initial_flux.abs()) {
        update.report.thermal_residual_j = -initial_flux;
        poll(&mut cancelled)?;
        return Ok(update);
    }
    let [wall_low, wall_high] = wall.heat_limits(contact, initial_gas.temperature)?;
    let domain = gas.domain();
    let kinetic = 0.5 * update.momentum * (update.momentum / initial_gas.density);
    let target = update.energy - kinetic;
    let gas_change_at = |temperature| -> Result<f64, WallExchangeError> {
        let energy = gas.caloric(state.partial_densities, temperature)?.0;
        let change = (energy - target) * volume_m3;
        if !change.is_finite() {
            return Err(WallExchangeError::Invalid("material wall gas energy range"));
        }
        Ok(change)
    };
    // Source-energy bounds also retain exact latent intervals, whose inventory
    // is not determined by temperature. A held-flux bound is valid only for
    // fixed conductance: an endpoint material law can increase that flux.
    let flux_limit = if state_dependent { f64::INFINITY } else { initial_flux.abs() };
    let (mut low, mut high) = if initial_wall > initial_gas.temperature {
        (
            0.0,
            gas_change_at(initial_wall.min(domain.maximum_temperature))?
                .min(wall_high)
                .min(flux_limit),
        )
    } else {
        (
            gas_change_at(initial_wall.max(domain.minimum_temperature))?
                .max(wall_low)
                .max(-flux_limit),
            0.0,
        )
    };
    // Reconstructed endpoint energy can round past the source's exact bound.
    // Move the trial bound one represented value toward zero; accepted storage
    // is never clipped or extrapolated.
    if low < 0.0 {
        low = low.next_up();
    }
    if high > 0.0 {
        high = high.next_down();
    }
    if !low.is_finite() || !high.is_finite() || low >= high {
        return Err(WallExchangeError::Invalid(
            "material wall has no supported heat interval",
        ));
    }
    let equation = |heat: f64, cancelled: &mut dyn FnMut() -> bool| {
        let gas_state = gas_at(update.energy + heat / volume_m3, cancelled)?;
        let wall_state = wall.after(contact, -heat)?;
        let gas_capacity = gas
            .caloric(state.partial_densities, gas_state.temperature)?
            .1
            * volume_m3;
        let endpoint_integrated = conductance_at(wall_state.contact)? * duration.value();
        let flux = endpoint_integrated * (wall_state.temperature_k - gas_state.temperature);
        let residual = heat - flux;
        // A material source supplies values, not an exact derivative of G.
        // Do not pass the frozen-G tangent off as the nonlinear Jacobian.
        let slope = if state_dependent {
            None
        } else {
            wall_state.temperature_slope_k_per_j
                .map(|wall_slope| 1.0 + integrated * (wall_slope + 1.0 / gas_capacity))
        };
        if !positive(gas_capacity)
            || !nonnegative(endpoint_integrated)
            || slope.is_some_and(|slope| !positive(slope))
            || !residual.is_finite()
            || !flux.is_finite()
        {
            return Err(WallExchangeError::Invalid("material wall heat equation"));
        }
        Ok((residual, slope, heat.abs().max(flux.abs())))
    };
    let (heat, iterations) = solve_scalar(low, high, 0.0, control, equation, &mut cancelled)?;
    let energy = update.energy + heat / volume_m3;
    let actual_heat = (energy - update.energy) * volume_m3;
    if !actual_heat.is_finite() || (heat != 0.0 && actual_heat == 0.0) {
        return Err(WallExchangeError::Invalid(
            "material wall heat is unrepresentable",
        ));
    }
    let actual_gas = gas_at(energy, &mut cancelled)?;
    let actual_wall = wall.after(contact, -actual_heat)?;
    let actual_conductance = conductance_at(actual_wall.contact)?;
    let actual_integrated = actual_conductance * duration.value();
    let wall_change = actual_wall.energy_change_j;
    let energy_residual = actual_heat + wall_change;
    let thermal_residual =
        actual_heat - actual_integrated * (actual_wall.temperature_k - actual_gas.temperature);
    let tolerance = control.tolerance(actual_heat.abs().max(wall_change.abs()));
    if !nonnegative(actual_integrated)
        || !energy_residual.is_finite()
        || !thermal_residual.is_finite()
        || energy_residual.abs() > tolerance
        || thermal_residual.abs() > tolerance
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j: energy_residual,
            thermal_residual_j: thermal_residual,
        });
    }
    update.contact = with_conductance(actual_wall.contact, actual_conductance);
    update.energy = energy;
    update.report.heat_to_gas_j = actual_heat;
    update.report.wall_energy_change_j = wall_change;
    update.report.energy_residual_j = energy_residual;
    update.report.thermal_residual_j = thermal_residual;
    update.report.iterations = iterations;
    poll(&mut cancelled)?;
    Ok(update)
}

fn with_conductance(mut contact: WallContact, value: f64) -> WallContact {
    match &mut contact.thermal {
        WallThermal::Lumped { conductance_w_per_k, .. }
        | WallThermal::Caloric { conductance_w_per_k, .. } => *conductance_w_per_k = value,
        _ => unreachable!("binding admitted finite storage"),
    }
    contact
}
