//! Source-bearing molecular escape from a freely translating finite vessel.
//!
//! All separately exposed outlets spend the same original species inventory
//! and share one caloric reconstruction. Their actual accepted species masses
//! carry qualified half-Maxwellian momentum and molecular energy into the
//! native variable-mass mechanics. Exhaust momentum and energy are transformed
//! into the same inertial frame as the vessel, so recoil is funded by the gas.
//!
//! The gas is locally Maxwellian in the vessel frame. Outward normals remain
//! fixed in the declared inertial frame; rotation, wall accommodation, gravity,
//! external forces, nozzle pressure thrust and finite receiver filling are not
//! part of this free-translation model. Fluxes are held for each first-order
//! finite-gas step. The caller refines time as temperature/composition changes.

use super::GasRocketState;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_TRANSFERS,
    species::{
        effusion::momentum::EffusiveMomentumInterface,
        network::transport::{
            GasTransportError, GasTransportPaths, GasTransportSolve,
            effusion::vacuum::{self, VacuumEffusionLink, VacuumEffusionStep},
        },
    },
};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_mbd::{
    Vec3,
    rocket::{
        RocketError,
        exhaust::{ExhaustForcing, ExhaustStep},
    },
};
use std::collections::BTreeMap;

/// One physical aperture group; groups must not count the same exposed area twice.
#[derive(Clone, Copy, Debug)]
pub struct EffusionOutlet<'a> {
    pub id: u64,
    /// Original thermochemical and uncollimated angular qualification.
    pub interface: &'a EffusiveMomentumInterface,
    /// Exposed area [m²]; zero closes this outlet without querying its source.
    pub area_m2: f64,
    /// Finite unit vector pointing OUT of the vessel in its declared frame.
    pub outward_normal_world: Vec3,
}

/// Physical budgets for the shared material and mechanical transaction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionRocketSolve {
    pub transport: GasTransportSolve,
    /// Absolute componentwise momentum residual limit [N s].
    pub momentum_tolerance_n_s: f64,
}

/// A refusal publishes neither the thermodynamic nor mechanical endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum EffusionRocketError {
    Invalid(&'static str),
    Gas(GasTransportError),
    Source {
        outlet: u64,
        error: GasEnergyError,
    },
    Dynamics(RocketError),
    InventoryBalance {
        remaining_residual_kg: f64,
        consumed_residual_kg: f64,
    },
    MomentumBalance {
        residual_n_s: Vec3,
    },
    EnergyBalance {
        residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for EffusionRocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled effusive vacuum recoil: {self:?}")
    }
}
impl std::error::Error for EffusionRocketError {}
impl From<GasTransportError> for EffusionRocketError {
    fn from(error: GasTransportError) -> Self {
        Self::Gas(error)
    }
}
impl From<RocketError> for EffusionRocketError {
    fn from(error: RocketError) -> Self {
        Self::Dynamics(error)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EffusionRocketError> {
    if cancelled() {
        Err(EffusionRocketError::Cancelled)
    } else {
        Ok(())
    }
}

// Preserve small components between opposed outlets in canonical ID/species order.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), EffusionRocketError> {
        let next = self.value + term;
        self.correction += if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else {
            (term - next) + self.value
        };
        self.value = next;
        if self.total().is_finite() {
            Ok(())
        } else {
            Err(EffusionRocketError::Invalid(
                "unrepresentable aggregate exhaust",
            ))
        }
    }
    fn total(&self) -> f64 {
        self.value + self.correction
    }
}

/// Original accepted gas and mechanics plus their actual cross-owner balances.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionRocketStep {
    pub duration_s: f64,
    pub gas: VacuumEffusionStep,
    /// Native mechanics and inertial-frame exhaust from the held physical rates.
    /// Its raw endpoint mass is retained; the joint endpoint uses the actual
    /// gas inventory after the explicit mass-chart correction below is admitted.
    pub exhaust: ExhaustStep,
    pub expelled_mass_kg: f64,
    /// Accepted outward species impulse in the vessel frame [kg m/s].
    pub relative_momentum_kg_m_s: Vec3,
    /// Accepted molecular energy in the original species energy datum [J].
    pub relative_energy_j: f64,
    /// Actual accepted escaped momentum in the declared inertial frame.
    /// Includes source-to-rate rounding; `exhaust` retains the raw native report.
    pub exhaust_momentum_world_kg_m_s: Vec3,
    /// Accepted molecular energy plus its inertial kinetic transport [J].
    /// May be signed when original species use signed formation datums.
    pub exhaust_energy_world_j: f64,
    /// Native remaining propellant minus actual accepted gas inventory [kg].
    pub remaining_mass_residual_kg: f64,
    pub consumed_mass_residual_kg: f64,
    /// Inertial momentum correction from accepted masses/impulses to native rates.
    pub source_momentum_correction_n_s: Vec3,
    /// Corresponding kinetic transport correction, without formation datums.
    pub source_kinetic_correction_j: f64,
    /// Momentum change when the native endpoint uses the actual gas inventory.
    pub inventory_momentum_correction_n_s: Vec3,
    /// Corresponding kinetic-energy change, evaluated without bulk subtraction.
    pub inventory_kinetic_correction_j: f64,
    pub momentum_residual_kg_m_s: Vec3,
    /// Original source-relative gas closure plus native kinetic/exhaust closure
    /// and the admitted inventory correction. No formation energies are subtracted.
    pub energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionRocketCandidate {
    endpoint: GasRocketState,
    step: EffusionRocketStep,
}
impl EffusionRocketCandidate {
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    #[must_use]
    pub const fn step(&self) -> &EffusionRocketStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(self) -> (GasRocketState, EffusionRocketStep) {
        (self.endpoint, self.step)
    }
}

impl GasRocketState {
    /// Tentatively cool, fractionate and recoil through all specified vacuum outlets.
    ///
    /// `inertial_frame` must match this vessel; it declares the normals' frame,
    /// not a coordinate transformation. There is no incident vacuum gas and no
    /// added pressure-thrust term. Opposed outlets can cancel momentum while
    /// still removing species and molecular energy. Closed outlets or an empty
    /// list coast with the original gas intact. Duration must be positive.
    ///
    /// The original material owner's withdrawal, source-domain and final-state
    /// checks apply. Native mechanics consumes the accepted masses, never a
    /// resampled endpoint rate. Gas remains the sole propellant inventory:
    /// remapping floating-point native mass to that inventory is explicitly
    /// included in the momentum and energy budgets before either owner publishes.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion(
        &self,
        inertial_frame: u64,
        outlets: &[EffusionOutlet<'_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        controls: EffusionRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionRocketCandidate, EffusionRocketError> {
        poll(&mut cancelled)?;
        if inertial_frame != self.inertial_frame
            || !controls.momentum_tolerance_n_s.is_finite()
            || controls.momentum_tolerance_n_s < 0.0
            || outlets.len() > MAX_SNAPSHOT_TRANSFERS
        {
            return Err(EffusionRocketError::Invalid(
                "effusive recoil frame, momentum budget or outlet count",
            ));
        }
        let mut ordered = BTreeMap::new();
        for outlet in outlets {
            poll(&mut cancelled)?;
            let normal = outlet.outward_normal_world;
            let squared = normal.norm_squared();
            if !normal.is_finite()
                || !squared.is_finite()
                || (squared - 1.0).abs() > 64.0 * f64::EPSILON
                || ordered.insert(outlet.id, outlet).is_some()
            {
                return Err(EffusionRocketError::Invalid(
                    "effusive outlets need unique IDs and finite unit outward normals",
                ));
            }
        }
        let parcels = BTreeMap::from([(0, self.reservoir.clone())]);
        let endpoints = BTreeMap::from([(0, endpoint)]);
        let links: Vec<_> = ordered
            .values()
            .map(|outlet| VacuumEffusionLink {
                id: outlet.id,
                parcel: 0,
                interface: outlet.interface.interface(),
                area_m2: outlet.area_m2,
            })
            .collect();
        let gas = vacuum::step_candidate(
            &parcels,
            GasTransportPaths::default(),
            &links,
            duration_s,
            &endpoints,
            controls.transport,
            &mut cancelled,
        )?;
        let (mut parcels, gas) = gas.into_parts();
        let mut mass = Sum::default();
        let mut energy = Sum::default();
        let mut impulse: [Sum; 3] = core::array::from_fn(|_| Sum::default());
        let mut boundaries = gas.exchange.boundaries.iter().peekable();
        for outlet in &gas.outlets {
            poll(&mut cancelled)?;
            let authored = ordered[&outlet.id];
            let source = |error| EffusionRocketError::Source {
                outlet: outlet.id,
                error,
            };
            let momentum = authored
                .interface
                .vacuum_flux(&outlet.flux, &mut cancelled)
                .map_err(source)?;
            while boundaries.peek().is_some_and(|row| row.id == outlet.id) {
                poll(&mut cancelled)?;
                let row = boundaries.next().ok_or(EffusionRocketError::Invalid(
                    "accepted boundary stream vanished",
                ))?;
                let index = momentum
                    .species()
                    .binary_search_by(|beam| beam.species.cmp(&row.transfer.species))
                    .map_err(|_| {
                        EffusionRocketError::Invalid(
                            "accepted species lost its original effusive beam",
                        )
                    })?;
                let dm = row.transfer.outgoing_mass_kg;
                let emitted = momentum.species()[index]
                    .outward
                    .impulse_for_mass_kg(dm)
                    .map_err(source)?;
                mass.add(dm)?;
                energy.add(-row.energy_into_gas_j)?;
                for (sum, direction) in impulse.iter_mut().zip([
                    authored.outward_normal_world.x,
                    authored.outward_normal_world.y,
                    authored.outward_normal_world.z,
                ]) {
                    sum.add(direction * emitted)?;
                }
            }
        }
        if boundaries.next().is_some() {
            return Err(EffusionRocketError::Invalid(
                "accepted escape lost its physical outlet",
            ));
        }
        let expelled_mass_kg = mass.total();
        let relative_energy_j = energy.total();
        let relative_momentum_kg_m_s =
            Vec3::new(impulse[0].total(), impulse[1].total(), impulse[2].total());
        let forcing = ExhaustForcing {
            mass_flow_kg_s: expelled_mass_kg / duration_s,
            relative_momentum_rate_n: Vec3::new(
                relative_momentum_kg_m_s.x / duration_s,
                relative_momentum_kg_m_s.y / duration_s,
                relative_momentum_kg_m_s.z / duration_s,
            ),
            relative_energy_rate_w: relative_energy_j / duration_s,
        };
        let exhaust = self
            .vehicle
            .advance_with_exhaust(forcing, duration_s, &mut cancelled)?;
        let reservoir = parcels.remove(&0).ok_or(EffusionRocketError::Invalid(
            "accepted gas inventory vanished",
        ))?;
        let mut vehicle = exhaust.mechanics.state_after;
        let remaining_mass_residual_kg = vehicle.propellant_mass_kg - reservoir.mass_kg();
        let consumed_mass_residual_kg =
            exhaust.mechanics.burned_propellant_mass_kg - expelled_mass_kg;
        let initial_mass_residual = self.vehicle.propellant_mass_kg - self.reservoir.mass_kg();
        let mass_budget = controls.transport.exchange.mass_tolerance_kg;
        if !remaining_mass_residual_kg.is_finite()
            || !consumed_mass_residual_kg.is_finite()
            || remaining_mass_residual_kg.abs() > mass_budget
            || consumed_mass_residual_kg.abs() > mass_budget
            || initial_mass_residual.abs() > mass_budget
            || (expelled_mass_kg > 0.0 && exhaust.mechanics.coast_duration_s != 0.0)
        {
            return Err(EffusionRocketError::InventoryBalance {
                remaining_residual_kg: remaining_mass_residual_kg,
                consumed_residual_kg: consumed_mass_residual_kg,
            });
        }
        let delta_v = vehicle.velocity_m_s.sub(self.vehicle.velocity_m_s);
        let mass_chart_change = initial_mass_residual - remaining_mass_residual_kg;
        let inventory_momentum_correction_n_s = self
            .vehicle
            .velocity_m_s
            .scale(mass_chart_change)
            .sub(delta_v.scale(remaining_mass_residual_kg));
        let inventory_kinetic_correction_j =
            0.5 * mass_chart_change * self.vehicle.velocity_m_s.norm_squared()
                - remaining_mass_residual_kg
                    * (self.vehicle.velocity_m_s.dot(delta_v) + 0.5 * delta_v.norm_squared());
        let mass_rate_correction = expelled_mass_kg - exhaust.mechanics.burned_propellant_mass_kg;
        let impulse_rate_correction = relative_momentum_kg_m_s.sub(exhaust.relative_momentum_n_s);
        let source_momentum_correction_n_s = exhaust
            .mean_velocity_world_m_s
            .scale(mass_rate_correction)
            .add(impulse_rate_correction);
        let source_kinetic_correction_j = impulse_rate_correction
            .dot(exhaust.mean_velocity_world_m_s)
            + 0.5 * mass_rate_correction * exhaust.mean_speed_squared_m2_s2;
        let exhaust_momentum_world_kg_m_s = exhaust
            .exhaust_momentum_world_kg_m_s
            .add(source_momentum_correction_n_s);
        let exhaust_energy_world_j = relative_energy_j
            + exhaust.relative_kinetic_transport_energy_j
            + source_kinetic_correction_j;
        if !exhaust_momentum_world_kg_m_s.is_finite() || !exhaust_energy_world_j.is_finite() {
            return Err(EffusionRocketError::Invalid(
                "unrepresentable accepted inertial exhaust",
            ));
        }
        let momentum_residual_kg_m_s = exhaust
            .momentum_residual_kg_m_s
            .add(inventory_momentum_correction_n_s)
            .add(source_momentum_correction_n_s);
        if !momentum_residual_kg_m_s.is_finite()
            || [
                momentum_residual_kg_m_s.x,
                momentum_residual_kg_m_s.y,
                momentum_residual_kg_m_s.z,
            ]
            .into_iter()
            .any(|value| value.abs() > controls.momentum_tolerance_n_s)
        {
            return Err(EffusionRocketError::MomentumBalance {
                residual_n_s: momentum_residual_kg_m_s,
            });
        }
        let energy_residual_j = gas.exchange.total_energy_residual_j
            + exhaust.kinetic_energy_residual_j
            + inventory_kinetic_correction_j
            + source_kinetic_correction_j;
        if !energy_residual_j.is_finite()
            || energy_residual_j.abs() > controls.transport.exchange.energy_tolerance_j
        {
            return Err(EffusionRocketError::EnergyBalance {
                residual_j: energy_residual_j,
            });
        }
        vehicle.propellant_mass_kg = reservoir.mass_kg();
        poll(&mut cancelled)?;
        let endpoint = Self {
            reservoir,
            vehicle,
            inertial_frame: self.inertial_frame,
        };
        Ok(EffusionRocketCandidate {
            endpoint,
            step: EffusionRocketStep {
                duration_s,
                gas,
                exhaust,
                expelled_mass_kg,
                relative_momentum_kg_m_s,
                relative_energy_j,
                exhaust_momentum_world_kg_m_s,
                exhaust_energy_world_j,
                remaining_mass_residual_kg,
                consumed_mass_residual_kg,
                source_momentum_correction_n_s,
                source_kinetic_correction_j,
                inventory_momentum_correction_n_s,
                inventory_kinetic_correction_j,
                momentum_residual_kg_m_s,
                energy_residual_j,
            },
        })
    }
}

#[cfg(test)]
mod tests;
