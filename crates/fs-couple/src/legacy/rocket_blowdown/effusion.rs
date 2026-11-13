//! Source-bearing molecular escape and explicitly forced vessel translation.
//!
//! All separately exposed outlets spend the same original species inventory
//! and share one caloric reconstruction. Their actual accepted species masses
//! carry qualified half-Maxwellian momentum and molecular energy into the
//! native variable-mass mechanics. Exhaust momentum and energy are transformed
//! into the same inertial frame as the vessel, so recoil is funded by the gas.
//!
//! The gas is locally Maxwellian in the vessel frame. Outward normals remain
//! fixed in the declared inertial frame. Uniform gravity and another held force
//! supply their original impulse and mechanical work. Rotation, wall
//! accommodation, nozzle pressure thrust and finite receiver filling are not
//! modeled. Fluxes and forcing are held for each first-order finite-gas step;
//! changing temperature, composition or loads require caller time refinement.

/// Bounded adaptive trajectories using the same native half-step transactions.
pub mod trajectory;

/// Bidirectional exchange with explicit comoving molecular-source ports.
pub mod exchange;

/// Finite material-wall heat feedback on the same freely recoiling gas inventory.
pub mod thermal;

/// Finite diffuse channels with material-wall feedback and resolved exit momentum.
pub mod channel;

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
use fs_math::det;
use fs_mbd::{
    Vec3,
    rocket::{
        RocketError,
        exhaust::{ExhaustEnvironment, ExhaustForcing, ExhaustStep},
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

/// Loads held in the same inertial frame as the vessel and outlet normals.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionRocketEnvironment {
    /// Frame identity, without an inferred coordinate transformation.
    pub inertial_frame: u64,
    /// Uniform gravitational acceleration [m/s²]; its force follows gas inventory.
    pub gravity_world_m_s2: Vec3,
    /// Other held force [N], excluding gravity and all molecular particle momentum.
    pub external_force_world_n: Vec3,
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
    StreamInventoryBalance {
        remaining_residual_kg: f64,
        incoming_residual_kg: f64,
        outgoing_residual_kg: f64,
    },
    MomentumBalance {
        residual_n_s: Vec3,
    },
    EnergyBalance {
        residual_j: f64,
    },
    /// Remaining physical budgets cannot cover numerical inventory uncertainty.
    RoundingBudget {
        momentum_residual_n_s: Vec3,
        energy_residual_j: f64,
        momentum_bound_n_s: Vec3,
        energy_bound_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for EffusionRocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled effusive vessel recoil: {self:?}")
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

// Gravity is a mass-dependent force. Do not invent when within the interval a
// finite gas reconstruction's inventory rounding occurred. Its impulse is
// bounded by epsilon_m*|g|*dt; Cauchy--Schwarz bounds work using the native
// powered trajectory's mean speed squared. No mass escapes during a dry coast.
fn nonnegative_product3(mut factors: [f64; 3]) -> Result<f64, EffusionRocketError> {
    if factors
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(EffusionRocketError::Invalid(
            "gravity rounding reserve needs finite nonnegative factors",
        ));
    }
    if factors.contains(&0.0) {
        return Ok(0.0);
    }
    factors.sort_by(f64::total_cmp);
    // For positive factors this pairing cannot overflow or vanish unless the
    // final triple product is also outside the positive representable range.
    let product = (factors[0] * factors[2]) * factors[1];
    if !product.is_finite() || product == 0.0 {
        return Err(EffusionRocketError::Invalid(
            "positive gravity rounding reserve vanished or overflowed",
        ));
    }
    Ok(product)
}

fn gravity_rounding_bounds(
    mass_uncertainty_kg: f64,
    gravity: Vec3,
    duration_s: f64,
    mean_speed_squared_m2_s2: f64,
) -> Result<(Vec3, f64), EffusionRocketError> {
    if !mass_uncertainty_kg.is_finite() || mass_uncertainty_kg < 0.0 {
        return Err(EffusionRocketError::Invalid(
            "gravity mass-history uncertainty is unrepresentable",
        ));
    }
    if mass_uncertainty_kg == 0.0 || gravity == Vec3::ZERO {
        return Ok((Vec3::ZERO, 0.0));
    }
    let component =
        |value: f64| nonnegative_product3([mass_uncertainty_kg, duration_s, value.abs()]);
    let impulse = Vec3::new(
        component(gravity.x)?,
        component(gravity.y)?,
        component(gravity.z)?,
    );
    let scale = impulse.x.max(impulse.y).max(impulse.z);
    let work = if scale == 0.0 {
        0.0
    } else {
        let unit = Vec3::new(impulse.x / scale, impulse.y / scale, impulse.z / scale);
        nonnegative_product3([
            scale,
            det::sqrt(unit.norm_squared()),
            det::sqrt(mean_speed_squared_m2_s2),
        ])?
    };
    if !impulse.is_finite() || !work.is_finite() || work < 0.0 {
        return Err(EffusionRocketError::Invalid(
            "gravity mass-history rounding bound is unrepresentable",
        ));
    }
    Ok((impulse, work))
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
    /// Original native propellant minus authoritative starting gas inventory [kg].
    /// A local native starting copy uses that gas inventory; no gas is added.
    pub initial_mass_residual_kg: f64,
    /// Numerical change in the original native representation when canonicalized.
    /// Its magnitude is reserved in the physical budget, not supplied by a force.
    pub initial_inventory_momentum_correction_n_s: Vec3,
    pub initial_inventory_kinetic_correction_j: f64,
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
    /// Absolute inventory uncertainty behind the gravity reserve [kg]; zero for coast.
    pub gravity_mass_history_uncertainty_kg: f64,
    /// Componentwise nonnegative gravity impulse reserve [N s], separate from
    /// the signed residual and the native nominal external impulse.
    pub gravity_impulse_uncertainty_n_s: Vec3,
    /// Nonnegative gravity work reserve [J]; no exact rounded mass history is claimed.
    pub gravity_work_uncertainty_j: f64,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionRocketCandidate, EffusionRocketError> {
        self.advance_effusion_with_environment(
            EffusionRocketEnvironment {
                inertial_frame,
                gravity_world_m_s2: Vec3::ZERO,
                external_force_world_n: Vec3::ZERO,
            },
            outlets,
            duration_s,
            endpoint,
            controls,
            cancelled,
        )
    }

    /// Advance the same source-bearing escape under explicit gravity and force.
    ///
    /// The native trajectory integrates the force against the changing mass,
    /// and its external impulse/work enter the same momentum/energy balances.
    /// Closed outlets retain their gas while the vessel follows loaded coast.
    /// Gravity-induced internal stratification is outside the compact mixed-gas
    /// model. Geometry, loads and molecular fluxes are fixed over this interval.
    ///
    /// Gas is the authoritative initial and final propellant inventory. Initial
    /// native inventory remapping and unresolved gravity mass-history rounding
    /// reserve space in the existing physical budgets; neither is hidden as an
    /// exact force/work correction. All endpoints remain tentative until the
    /// final balance, original-source and cancellation checks succeed.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_with_environment(
        &self,
        environment: EffusionRocketEnvironment,
        outlets: &[EffusionOutlet<'_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        controls: EffusionRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionRocketCandidate, EffusionRocketError> {
        poll(&mut cancelled)?;
        if environment.inertial_frame != self.inertial_frame
            || !environment.gravity_world_m_s2.is_finite()
            || !environment.external_force_world_n.is_finite()
            || !controls.momentum_tolerance_n_s.is_finite()
            || controls.momentum_tolerance_n_s < 0.0
            || outlets.len() > MAX_SNAPSHOT_TRANSFERS
        {
            return Err(EffusionRocketError::Invalid(
                "effusive recoil frame, finite loads, momentum budget or outlet count",
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
        let initial_mass_residual_kg = self.vehicle.propellant_mass_kg - self.reservoir.mass_kg();
        let mass_budget = controls.transport.exchange.mass_tolerance_kg;
        if !initial_mass_residual_kg.is_finite() || initial_mass_residual_kg.abs() > mass_budget {
            return Err(EffusionRocketError::InventoryBalance {
                remaining_residual_kg: initial_mass_residual_kg,
                consumed_residual_kg: 0.0,
            });
        }
        let mut native_initial = self.vehicle;
        native_initial.propellant_mass_kg = self.reservoir.mass_kg();
        let initial_inventory_momentum_correction_n_s =
            self.vehicle.velocity_m_s.scale(-initial_mass_residual_kg);
        let initial_inventory_kinetic_correction_j = if initial_mass_residual_kg == 0.0 {
            0.0
        } else {
            (-0.5 * initial_mass_residual_kg) * self.vehicle.velocity_m_s.norm_squared()
        };
        let exhaust = native_initial.advance_with_exhaust_forced(
            forcing,
            ExhaustEnvironment {
                gravity_m_s2: environment.gravity_world_m_s2,
                external_force_n: environment.external_force_world_n,
            },
            duration_s,
            &mut cancelled,
        )?;
        let reservoir = parcels.remove(&0).ok_or(EffusionRocketError::Invalid(
            "accepted gas inventory vanished",
        ))?;
        let mut vehicle = exhaust.mechanics.state_after;
        let remaining_mass_residual_kg = vehicle.propellant_mass_kg - reservoir.mass_kg();
        let consumed_mass_residual_kg =
            exhaust.mechanics.burned_propellant_mass_kg - expelled_mass_kg;
        if !remaining_mass_residual_kg.is_finite()
            || !consumed_mass_residual_kg.is_finite()
            || remaining_mass_residual_kg.abs() > mass_budget
            || consumed_mass_residual_kg.abs() > mass_budget
            || (expelled_mass_kg > 0.0 && exhaust.mechanics.coast_duration_s != 0.0)
        {
            return Err(EffusionRocketError::InventoryBalance {
                remaining_residual_kg: remaining_mass_residual_kg,
                consumed_residual_kg: consumed_mass_residual_kg,
            });
        }
        let inventory_momentum_correction_n_s =
            vehicle.velocity_m_s.scale(-remaining_mass_residual_kg);
        let inventory_kinetic_correction_j = if remaining_mass_residual_kg == 0.0 {
            0.0
        } else {
            (-0.5 * remaining_mass_residual_kg) * vehicle.velocity_m_s.norm_squared()
        };
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
        let gravity_mass_history_uncertainty_kg =
            if expelled_mass_kg == 0.0 || environment.gravity_world_m_s2 == Vec3::ZERO {
                0.0
            } else {
                consumed_mass_residual_kg.abs()
                    + exhaust.mass_residual_kg.abs()
                    + remaining_mass_residual_kg.abs()
                    + gas.exchange.total_mass_residual_kg.abs()
            };
        let (gravity_impulse_uncertainty_n_s, gravity_work_uncertainty_j) =
            gravity_rounding_bounds(
                gravity_mass_history_uncertainty_kg,
                environment.gravity_world_m_s2,
                exhaust.mechanics.powered_duration_s,
                exhaust.mean_speed_squared_m2_s2,
            )?;
        let momentum_bound = gravity_impulse_uncertainty_n_s.add(Vec3::new(
            initial_inventory_momentum_correction_n_s.x.abs(),
            initial_inventory_momentum_correction_n_s.y.abs(),
            initial_inventory_momentum_correction_n_s.z.abs(),
        ));
        let energy_bound =
            gravity_work_uncertainty_j + initial_inventory_kinetic_correction_j.abs();
        if !momentum_bound.is_finite()
            || !energy_bound.is_finite()
            || [
                (momentum_residual_kg_m_s.x, momentum_bound.x),
                (momentum_residual_kg_m_s.y, momentum_bound.y),
                (momentum_residual_kg_m_s.z, momentum_bound.z),
            ]
            .into_iter()
            .any(|(residual, bound)| bound > controls.momentum_tolerance_n_s - residual.abs())
            || energy_bound
                > controls.transport.exchange.energy_tolerance_j - energy_residual_j.abs()
        {
            return Err(EffusionRocketError::RoundingBudget {
                momentum_residual_n_s: momentum_residual_kg_m_s,
                energy_residual_j,
                momentum_bound_n_s: momentum_bound,
                energy_bound_j: energy_bound,
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
                initial_mass_residual_kg,
                initial_inventory_momentum_correction_n_s,
                initial_inventory_kinetic_correction_j,
                remaining_mass_residual_kg,
                consumed_mass_residual_kg,
                source_momentum_correction_n_s,
                source_kinetic_correction_j,
                inventory_momentum_correction_n_s,
                inventory_kinetic_correction_j,
                gravity_mass_history_uncertainty_kg,
                gravity_impulse_uncertainty_n_s,
                gravity_work_uncertainty_j,
                momentum_residual_kg_m_s,
                energy_residual_j,
            },
        })
    }
}

#[cfg(test)]
mod tests;
