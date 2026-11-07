//! Original material storage in the existing gas/wall scalar exchange.

use super::*;
#[cfg(feature = "material")]
use fs_material::thermal::transport::storage::PreparedHeatStorage;
use fs_thermal::enthalpy::{DepositedCell, ThermalLaw, ThermalStorageBasis};
use std::sync::Arc;

mod thermal_law;
use thermal_law::EnergyStorage;
#[cfg(feature = "material")]
mod reservoir;
#[cfg(feature = "material")]
pub(crate) use reservoir::{
    exchange_caloric_with_reservoir_conduction, exchange_cell_with_reservoir_conduction,
};

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
    /// Original sensible-capacity integral without an established Cp/Cv or
    /// pressure/deformation-work convention. Heat-only exchange retains this
    /// coordinate; it is neither stream enthalpy nor rigid-body internal energy.
    SourceSensibleStorage,
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

pub(crate) struct StorageUpdate {
    pub(crate) contact: WallContact,
    pub(crate) temperature_k: f64,
    pub(crate) energy_change_j: f64,
    pub(crate) temperature_slope_k_per_j: Option<f64>,
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
            ThermalStorageBasis::DeclaredSensibleStorage => {
                WallCaloricConstraint::SourceSensibleStorage
            }
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
                curve.clone(),
                initial,
                energy_range_j_kg,
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
            WallCaloricConstraint::SourceSensibleStorage => {
                ThermalStorageBasis::DeclaredSensibleStorage
            }
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

    /// Actual finite-store endpoint difference in this original source's h/u
    /// coordinate. A mechanical composing owner must admit its work convention
    /// before using the result as internal energy. Energy gauges cancel here.
    pub(crate) fn energy_change_j(
        &self,
        before: WallContact,
        after: WallContact,
    ) -> Result<f64, WallExchangeError> {
        match &self.storage {
            Storage::Energy(source) => {
                let (initial, _) = source.state(before)?;
                let (final_energy, _) = source.state(after)?;
                Ok(source.initial.mass_kg * (final_energy - initial))
            }
            #[cfg(feature = "material")]
            Storage::Capacity(source) => source
                .energy_change_j(
                    self.temperature_k(before.thermal)?,
                    self.temperature_k(after.thermal)?,
                )
                .map_err(WallExchangeError::from),
        }
    }

    /// Signed heat into gas that reaches each end of the wall source domain.
    pub(crate) fn heat_limits(
        &self,
        contact: WallContact,
        gas_temperature: f64,
    ) -> Result<[f64; 2], WallExchangeError> {
        match &self.storage {
            Storage::Energy(source) => {
                let _ = gas_temperature;
                let (energy, _) = source.state(contact)?;
                let limit = |endpoint: f64| {
                    let mut heat = source.initial.mass_kg * (energy - endpoint);
                    // Form only bracket proposals here. Multiplication by
                    // mass and reconstruction through division/addition can
                    // round a nominal endpoint outside its native h/u domain.
                    // Contract extensive heat toward zero until that SAME
                    // represented state is supported; never clip stored h/u.
                    for _ in 0..8 {
                        let candidate = energy + (-heat) / source.initial.mass_kg;
                        if !heat.is_finite()
                            || (source.range_j_kg[0]..=source.range_j_kg[1]).contains(&candidate)
                        {
                            return Ok(heat);
                        }
                        heat = if heat > 0.0 {
                            heat.next_down()
                        } else {
                            heat.next_up()
                        };
                    }
                    Err(WallExchangeError::Invalid(
                        "unrepresentable caloric wall heat bound",
                    ))
                };
                Ok([limit(source.range_j_kg[1])?, limit(source.range_j_kg[0])?])
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

    pub(crate) fn after(
        &self,
        contact: WallContact,
        heat: f64,
    ) -> Result<StorageUpdate, WallExchangeError> {
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

/// One original endpoint heat law. Derivatives, when supplied, are of POWER
/// with respect to [Twall,Tgas], including every material and film resistance.
pub(crate) struct WallHeatResponse {
    pub(crate) conductance_w_per_k: f64,
    pub(crate) power_to_gas_w: f64,
    pub(crate) temperature_derivative_w_per_k: Option<[f64; 2]>,
    pub(crate) flux_error_w: f64,
}

impl WallHeatResponse {
    fn validate(&self) -> Result<(), WallExchangeError> {
        if !nonnegative(self.conductance_w_per_k)
            || !self.power_to_gas_w.is_finite()
            || !nonnegative(self.flux_error_w)
            || self
                .temperature_derivative_w_per_k
                .is_some_and(|d| !nonnegative(d[0]) || !nonnegative(-d[1]))
        {
            return Err(WallExchangeError::Invalid(
                "material wall endpoint heat law",
            ));
        }
        Ok(())
    }
}

pub(crate) struct WallEndpointConduction<'a> {
    // Some: the heat law depends on both temperatures; refine gas against this
    // global slope bound before querying it. None: the original wall-only k law
    // supplies its coefficient first, retaining its original caloric precision.
    pub(crate) gas_sensitivity_bound_w_per_k: Option<f64>,
    pub(crate) resolve: &'a mut dyn FnMut(
        WallContact,
        [f64; 2],
        &mut dyn FnMut() -> bool,
    ) -> Result<WallHeatResponse, WallExchangeError>,
}

/// The same finite-storage owner with an optional source-owned endpoint heat
/// law. Its callback reads original wall storage and both actual temperatures;
/// the contact's conductance field retains the original film. No candidate is
/// published here. An exact integrated profile can supply its power tangents.
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
    conductance: Option<WallEndpointConduction<'_>>,
    control: WallControl,
    cancelled: impl FnMut() -> bool,
) -> Result<WallUpdate, WallExchangeError> {
    exchange_caloric_with_heat_storage(
        &MixtureCaloric {
            gas,
            partial_densities: state.partial_densities,
        },
        state.momentum,
        state.energy,
        volume_m3,
        duration,
        contact,
        wall,
        conductance,
        control,
        cancelled,
    )
}

/// Borrow the unchanged gas energy owner; wall source and implicit equation
/// are shared with the original species-carrying entry points above.
#[allow(clippy::too_many_arguments)]
pub(crate) fn exchange_caloric_with_heat_storage(
    gas: &impl WallGasCaloric,
    initial_momentum: f64,
    initial_energy: f64,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    wall: &MaterialWallCaloric<'_>,
    mut conductance: Option<WallEndpointConduction<'_>>,
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
    let gas_sensitivity_bound = conductance
        .as_ref()
        .and_then(|law| law.gas_sensitivity_bound_w_per_k);
    let mut response_at = |candidate: WallContact,
                           temperatures: [f64; 2],
                           cancelled: &mut dyn FnMut() -> bool| {
        let response = match &mut conductance {
            Some(law) => (law.resolve)(candidate, temperatures, cancelled)?,
            None => WallHeatResponse {
                conductance_w_per_k: authored_conductance,
                power_to_gas_w: authored_conductance * (temperatures[0] - temperatures[1]),
                temperature_derivative_w_per_k: Some([authored_conductance, -authored_conductance]),
                flux_error_w: 0.0,
            },
        };
        response.validate()?;
        Ok(response)
    };
    // Reuse the original radial loss and gas sensible-heating owner exactly.
    let mut update = exchange_caloric(
        gas,
        initial_momentum,
        initial_energy,
        volume_m3,
        duration,
        WallContact {
            drag: contact.drag,
            thermal: WallThermal::Adiabatic,
        },
        control,
        &mut cancelled,
    )?;
    let gas_at = |energy,
                  wall_temperature,
                  precision_integrated,
                  heat_scale,
                  cancelled: &mut dyn FnMut() -> bool| {
        let mut state = gas.primitive(update.momentum, energy, &mut *cancelled)?;
        let target = energy - 0.5 * update.momentum * state.velocity;
        state.temperature = if gas_sensitivity_bound.is_some() {
            refine_caloric_temperature_with_allowance(
                gas,
                state.temperature,
                target,
                volume_m3,
                precision_integrated,
                control.tolerance(heat_scale),
                control,
                cancelled,
            )?
        } else {
            refine_caloric_temperature(
                gas,
                state.temperature,
                target,
                volume_m3,
                precision_integrated,
                wall_temperature,
                heat_scale,
                control,
                cancelled,
            )?
        };
        Ok::<_, WallExchangeError>(state)
    };
    let mut endpoint_at = |candidate: WallContact,
                           wall_temperature: f64,
                           energy: f64,
                           heat_scale: f64,
                           cancelled: &mut dyn FnMut() -> bool| {
        // A wall-only coefficient does not depend on gas T. Resolve it before
        // the original gas inverse, then form its power at the actual endpoint.
        let sampled = if gas_sensitivity_bound.is_none() {
            Some(response_at(candidate, [wall_temperature; 2], cancelled)?)
        } else {
            None
        };
        let bound = sampled.as_ref().map_or_else(
            || gas_sensitivity_bound.expect("two-temperature law has a slope bound"),
            |response| response.conductance_w_per_k,
        );
        let integrated_bound = bound * duration.value();
        if !nonnegative(integrated_bound) {
            return Err(WallExchangeError::Invalid(
                "integrated material wall conductance",
            ));
        }
        let gas_state = gas_at(
            energy,
            wall_temperature,
            integrated_bound,
            heat_scale,
            cancelled,
        )?;
        let response = if let Some(mut response) = sampled {
            response.power_to_gas_w =
                response.conductance_w_per_k * (wall_temperature - gas_state.temperature);
            response
        } else {
            response_at(
                candidate,
                [wall_temperature, gas_state.temperature],
                cancelled,
            )?
        };
        Ok::<_, WallExchangeError>((gas_state, response))
    };
    let (initial_gas, initial_response) =
        endpoint_at(contact, initial_wall, update.energy, 0.0, &mut cancelled)?;
    update.contact = with_conductance(contact, initial_response.conductance_w_per_k);
    let integrated = initial_response.conductance_w_per_k * duration.value();
    if !nonnegative(integrated) {
        return Err(WallExchangeError::Invalid(
            "integrated material wall conductance",
        ));
    }
    if integrated == 0.0 || initial_wall == initial_gas.temperature {
        poll(&mut cancelled)?;
        return Ok(update);
    }
    let initial_flux = duration.value() * initial_response.power_to_gas_w;
    if !initial_flux.is_finite() {
        return Err(WallExchangeError::Invalid("material wall heat increment"));
    }
    if initial_flux.abs() + duration.value() * initial_response.flux_error_w
        <= control.tolerance(initial_flux.abs())
    {
        update.report.thermal_residual_j = -initial_flux;
        poll(&mut cancelled)?;
        return Ok(update);
    }
    let [wall_low, wall_high] = wall.heat_limits(contact, initial_gas.temperature)?;
    let domain = gas.domain();
    let kinetic = 0.5 * update.momentum * (update.momentum / initial_gas.density);
    let target = update.energy - kinetic;
    let gas_change_at = |temperature| -> Result<f64, WallExchangeError> {
        let energy = gas.caloric(temperature)?.0;
        let change = (energy - target) * volume_m3;
        if !change.is_finite() {
            return Err(WallExchangeError::Invalid("material wall gas energy range"));
        }
        Ok(change)
    };
    // Source-energy bounds also retain exact latent intervals, whose inventory
    // is not determined by temperature. A held-flux bound is valid only for
    // fixed conductance: an endpoint material law can increase that flux.
    let flux_limit = if gas_sensitivity_bound.is_some() {
        // Positive k(T) makes the original integrated material/film power
        // decrease along passive exchange. Include constitutive and gas-inverse
        // error in its initial-flux bound; point-sampled k has no such guarantee.
        initial_flux.abs()
            + duration.value() * initial_response.flux_error_w
            + control.tolerance(initial_flux.abs())
    } else if state_dependent {
        f64::INFINITY
    } else {
        initial_flux.abs()
    };
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
    let mut equation = |heat: f64, cancelled: &mut dyn FnMut() -> bool| {
        let wall_state = wall.after(contact, -heat)?;
        let (gas_state, response) = endpoint_at(
            wall_state.contact,
            wall_state.temperature_k,
            update.energy + heat / volume_m3,
            heat.abs(),
            cancelled,
        )?;
        let gas_capacity = gas.caloric(gas_state.temperature)?.1 * volume_m3;
        let flux = duration.value() * response.power_to_gas_w;
        let residual = heat - flux;
        // Only exact power derivatives from the original law enter Newton.
        // Point-sampled state-dependent k still supplies no invented dG/dh.
        let slope = response
            .temperature_derivative_w_per_k
            .zip(wall_state.temperature_slope_k_per_j)
            .map(|(derivatives, wall_slope)| {
                1.0 + duration.value()
                    * (derivatives[0] * wall_slope - derivatives[1] / gas_capacity)
            });
        if !positive(gas_capacity)
            || slope.is_some_and(|slope| !positive(slope))
            || !residual.is_finite()
            || !flux.is_finite()
        {
            return Err(WallExchangeError::Invalid("material wall heat equation"));
        }
        Ok((residual, slope, heat.abs().max(flux.abs())))
    };
    // Reserve part of the physical budget for a nested constitutive flux solve.
    // The final check below includes that law's reported power-error bound.
    let mut solve_control =
        if initial_response.temperature_derivative_w_per_k.is_some() && state_dependent {
            WallControl {
                absolute_energy_tolerance_j: 0.5 * control.absolute_energy_tolerance_j,
                relative_energy_tolerance: 0.5 * control.relative_energy_tolerance,
                ..control
            }
        } else {
            control
        };
    let bracket_iterations = if gas_sensitivity_bound.is_some() {
        let far = if initial_flux > 0.0 { high } else { low };
        let (supported, iterations) =
            supported_profile_heat_bracket(0.0, far, solve_control, &mut equation, &mut cancelled)?;
        if initial_flux > 0.0 {
            high = supported;
        } else {
            low = supported;
        }
        solve_control.maximum_iterations -= iterations;
        iterations
    } else {
        0
    };
    let (heat, iterations) = solve_scalar(low, high, 0.0, solve_control, equation, &mut cancelled)?;
    let iterations = iterations + bracket_iterations;
    let energy = update.energy + heat / volume_m3;
    let actual_heat = (energy - update.energy) * volume_m3;
    if !actual_heat.is_finite() || (heat != 0.0 && actual_heat == 0.0) {
        return Err(WallExchangeError::Invalid(
            "material wall heat is unrepresentable",
        ));
    }
    let actual_wall = wall.after(contact, -actual_heat)?;
    let (_, actual_response) = endpoint_at(
        actual_wall.contact,
        actual_wall.temperature_k,
        energy,
        actual_heat.abs(),
        &mut cancelled,
    )?;
    let actual_conductance = actual_response.conductance_w_per_k;
    let actual_integrated = actual_conductance * duration.value();
    let wall_change = actual_wall.energy_change_j;
    let energy_residual = update.report.energy_residual_j + actual_heat + wall_change;
    let thermal_residual = actual_heat - duration.value() * actual_response.power_to_gas_w;
    let scale = actual_heat.abs().max(wall_change.abs());
    if !nonnegative(actual_integrated)
        || !energy_residual.is_finite()
        || !thermal_residual.is_finite()
        || energy_residual.abs()
            > control.tolerance(scale.max(update.report.mechanical_work_to_gas_j.abs()))
        || thermal_residual.abs() + duration.value() * actual_response.flux_error_w
            > control.tolerance(scale)
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

/// Find a supported far bracket for a passive integrated profile. Storage
/// bounds need not support the material FACE, even when the physical heat root
/// does. Search between the last supported wrong-sign point and the unsupported
/// endpoint; never clip a material temperature, invent k, or accept a wrong sign.
pub(super) fn supported_profile_heat_bracket(
    anchor: f64,
    far: f64,
    control: WallControl,
    equation: &mut impl FnMut(
        f64,
        &mut dyn FnMut() -> bool,
    ) -> Result<(f64, Option<f64>, f64), WallExchangeError>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<(f64, usize), WallExchangeError> {
    let direction = (far - anchor).signum();
    let mut supported = anchor;
    let mut unsupported = far;
    let mut trial = unsupported;
    let mut domain_error = None;
    for contractions in 0..control.maximum_iterations {
        if cancelled() {
            return Err(WallExchangeError::Cancelled);
        }
        match equation(trial, cancelled) {
            Ok((residual, _, scale)) => {
                if direction * residual >= -control.tolerance(scale) || domain_error.is_none() {
                    return Ok((trial, contractions));
                }
                supported = trial;
            }
            Err(
                error @ WallExchangeError::Thermal(
                    fs_thermal::enthalpy::DepositionError::MaterialDomain,
                ),
            ) => {
                domain_error = Some(error);
                unsupported = trial;
            }
            Err(error) => return Err(error),
        }
        let next = 0.5 * supported + 0.5 * unsupported;
        if next <= supported.min(unsupported) || next >= supported.max(unsupported) {
            return Err(domain_error.unwrap_or(WallExchangeError::TemperatureIterations));
        }
        trial = next;
    }
    Err(WallExchangeError::TemperatureIterations)
}

fn with_conductance(mut contact: WallContact, value: f64) -> WallContact {
    match &mut contact.thermal {
        WallThermal::Lumped {
            conductance_w_per_k,
            ..
        }
        | WallThermal::Caloric {
            conductance_w_per_k,
            ..
        } => *conductance_w_per_k = value,
        _ => unreachable!("binding admitted finite storage"),
    }
    contact
}
