//! Geometry-resolved molecular escape, finite channel walls and vessel recoil.
//!
//! Every outlet spends one original gas inventory. The material channel owner
//! retains gross reflected molecules, pools accommodation heat by physical wall
//! identity and reconstructs gas and finite material storage together. Its
//! accepted escaping masses acquire axial momentum from the SAME qualified
//! Clausing geometry, without a second sample or thin-aperture approximation.
//!
//! The translating vessel includes the gas AND channel walls. Incident and
//! returned impulses are internal exchanges, so only escaping momentum drives
//! its recoil. Local wall traction remains available in each channel receipt;
//! adding it to that recoil would count the same molecules twice. Complete
//! escaped molecular energy includes accommodation before its inertial-frame
//! transform along the native variable-mass trajectory.
//!
//! Distinct wall IDs denote disjoint inventories already included in vehicle
//! dry mass. Channels can share a wall ID; its mass and storage count once.
//! The quasistatic, first-order stage holds geometry, normals, fluxes and loads;
//! repeat with accepted gas/walls to resolve thermal and composition feedback.
//! Prepared material conduction can connect distinct finite wall inventories;
//! its equal and opposite heat joins molecular heat before one wall inversion.
//! Internal conduction changes subsequent escape without adding external thrust.
//! Rotation, transit delay, spatial gas gradients and continuum pressure thrust
//! are outside this point-vessel owner. Native source and gross withdrawal
//! limits remain active, and refusal publishes none of the physical endpoints.

/// Adaptive finite-channel gas, material-wall and vessel evolution.
pub mod trajectory;

use super::{
    EffusionRocketEnvironment, EffusionRocketError, EffusionRocketSolve, GasRocketState, Sum,
    gravity_rounding_bounds, poll,
};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::{
            snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
            species::{
                effusion::channel::momentum::{
                    ChannelMomentumStep, DiffuseChannelMomentumInterface,
                },
                network::transport::{
                    GasTransportPaths,
                    effusion::channel::thermal::{
                        FiniteChannelError, FixedVolumeChannelWall,
                        conduction::ChannelWallConduction,
                        vacuum::{self, FiniteVacuumChannelLink, FiniteVacuumChannelStep},
                    },
                },
            },
        },
    },
};
use fs_math::dd::Dd;
use fs_mbd::{
    Vec3,
    rocket::exhaust::{ExhaustEnvironment, ExhaustForcing, ExhaustStep},
};
use std::collections::BTreeMap;

/// A separately exposed circular-channel array connected to a finite wall.
#[derive(Clone, Copy, Debug)]
pub struct EffusionChannelOutlet<'a> {
    pub id: u64,
    /// Original geometry, species, accommodation and speed-distribution laws.
    pub interface: &'a DiffuseChannelMomentumInterface,
    /// Sum of physical mouth areas [m²]; zero closes this outlet.
    pub area_m2: f64,
    /// Stable identity in the finite wall map, shared by all its channels.
    pub wall: u64,
    /// Finite unit bore axis pointing OUT of the vessel in its inertial frame.
    pub outward_normal_world: Vec3,
}

#[derive(Clone, Debug, PartialEq)]
pub enum EffusionChannelError {
    Invalid(&'static str),
    Channel(FiniteChannelError),
    Recoil(EffusionRocketError),
    /// Signed defects cannot cancel their absolute shares of the fixed budget.
    EnergyBalance {
        residual_j: f64,
        absolute_residual_sum_j: f64,
        uncertainty_bound_j: f64,
    },
}
impl core::fmt::Display for EffusionChannelError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-channel effusive recoil: {self:?}")
    }
}
impl std::error::Error for EffusionChannelError {}
impl From<FiniteChannelError> for EffusionChannelError {
    fn from(error: FiniteChannelError) -> Self {
        Self::Channel(error)
    }
}
impl From<EffusionRocketError> for EffusionChannelError {
    fn from(error: EffusionRocketError) -> Self {
        Self::Recoil(error)
    }
}

/// Actual accepted species impulses; the native receipt uses the stated axis.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionChannelMomentumStep {
    pub outward_normal_world: Vec3,
    /// `exterior_impulses_n_s[1]` escapes. The left exterior and local wall
    /// traction together are its opposite; they are not extra external loads.
    pub channel: ChannelMomentumStep,
}

/// Original variable-mass mechanics and explicit accepted-inventory corrections.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionChannelRecoil {
    pub exhaust: ExhaustStep,
    pub expelled_mass_kg: f64,
    pub relative_momentum_kg_m_s: Vec3,
    /// Full escaped energy, including channel accommodation and source datums.
    pub relative_energy_j: f64,
    pub exhaust_momentum_world_kg_m_s: Vec3,
    pub exhaust_energy_world_j: f64,
    /// Original native propellant minus authoritative initial gas mass [kg].
    pub initial_mass_residual_kg: f64,
    pub initial_inventory_momentum_correction_n_s: Vec3,
    pub initial_inventory_kinetic_correction_j: f64,
    pub remaining_mass_residual_kg: f64,
    pub consumed_mass_residual_kg: f64,
    /// Accepted amounts minus rate*duration, transported on the native path.
    pub source_momentum_correction_n_s: Vec3,
    pub source_kinetic_correction_j: f64,
    /// Native endpoint remapping to the actual material gas inventory.
    pub inventory_momentum_correction_n_s: Vec3,
    pub inventory_kinetic_correction_j: f64,
    pub gravity_mass_history_uncertainty_kg: f64,
    pub gravity_impulse_uncertainty_n_s: Vec3,
    pub gravity_work_uncertainty_j: f64,
    pub momentum_residual_kg_m_s: Vec3,
    /// Actual gas + wall + kinetic + escaped energy minus external work [J],
    /// assembled in the original source-relative chart without datum subtraction.
    pub energy_residual_j: f64,
    pub absolute_energy_residual_sum_j: f64,
    /// Initial mass remapping and unresolved gravity mass-history work [J].
    pub energy_uncertainty_bound_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionChannelStep {
    pub duration_s: f64,
    /// Original finite gas/wall transaction, including actual escaped energy.
    pub thermal: FiniteVacuumChannelStep,
    /// Canonical outlet-ID order, independent of caller ordering.
    pub momentum: Vec<EffusionChannelMomentumStep>,
    pub recoil: EffusionChannelRecoil,
}

/// All gas, translation and finite-wall endpoints are accepted together.
#[derive(Clone, Debug)]
pub struct EffusionChannelCandidate<'l, 'p> {
    endpoint: GasRocketState,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: EffusionChannelStep,
}
impl<'l, 'p> EffusionChannelCandidate<'l, 'p> {
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    pub fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> {
        &self.walls
    }
    pub const fn step(&self) -> &EffusionChannelStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        GasRocketState,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        EffusionChannelStep,
    ) {
        (self.endpoint, self.walls, self.step)
    }
}

impl GasRocketState {
    /// Advance finite diffuse channels, their material walls and vessel recoil
    /// over one physical interval using the original accepted molecular beams.
    ///
    /// Every distinct wall is already part of this vessel's dry mass; their sum
    /// must fit it. Walls sharing a material law may still be distinct physical
    /// inventories. To share one inventory, use the same wall ID on all outlets.
    /// The existing energy budget bounds the sum of absolute material/mechanical
    /// defects and inventory/gravity reserves, not only their signed sum.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_channels<'l, 'p>(
        &self,
        environment: EffusionRocketEnvironment,
        walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        outlets: &[EffusionChannelOutlet<'_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        controls: EffusionRocketSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionChannelCandidate<'l, 'p>, EffusionChannelError> {
        self.advance_effusion_channels_with_conduction(
            environment,
            walls,
            outlets,
            &[],
            duration_s,
            endpoint,
            controls,
            cancelled,
        )
    }

    /// Advance molecular escape and heat spreading through connected walls.
    ///
    /// Each prepared conduction link retains the original wall caloric and
    /// conductivity sources, aligned transfer distances, area and contact
    /// resistance. All conductive and molecular heat uses the same initial
    /// snapshot and one material inversion per affected wall. Wall inventories
    /// without an outlet may conduct heat and still belong to the dry mass.
    /// Conductive heat is internal to this vessel; only escaping molecular
    /// momentum and the explicit environment supply external mechanical loads.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_channels_with_conduction<'l, 'p>(
        &self,
        environment: EffusionRocketEnvironment,
        walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        outlets: &[EffusionChannelOutlet<'_>],
        conduction: &[ChannelWallConduction<'_, '_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        controls: EffusionRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionChannelCandidate<'l, 'p>, EffusionChannelError> {
        poll(&mut cancelled)?;
        if environment.inertial_frame != self.inertial_frame
            || !environment.gravity_world_m_s2.is_finite()
            || !environment.external_force_world_n.is_finite()
            || !controls.momentum_tolerance_n_s.is_finite()
            || controls.momentum_tolerance_n_s < 0.0
            || walls.len() > MAX_SNAPSHOT_PARCELS
            || outlets.len() > MAX_SNAPSHOT_TRANSFERS
        {
            return Err(EffusionChannelError::Invalid(
                "finite-channel frame, loads, momentum budget or work count",
            ));
        }
        let mut wall_mass = Dd::ZERO;
        for wall in walls.values() {
            poll(&mut cancelled)?;
            wall_mass = wall_mass + Dd::from_f64(wall.mass_kg());
        }
        if !wall_mass.to_f64().is_finite()
            || (wall_mass - Dd::from_f64(self.vehicle.dry_mass_kg)).to_f64() > 0.0
        {
            return Err(EffusionChannelError::Invalid(
                "finite channel walls must be included once within existing dry mass",
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
                return Err(EffusionChannelError::Invalid(
                    "finite channel outlets need unique IDs and unit outward axes",
                ));
            }
        }
        let links: Vec<_> = ordered
            .values()
            .map(|outlet| FiniteVacuumChannelLink {
                id: outlet.id,
                parcel: 0,
                channel: outlet.interface.channel(),
                area_m2: outlet.area_m2,
                wall: outlet.wall,
            })
            .collect();
        let candidate = vacuum::step_candidate_with_conduction(
            &BTreeMap::from([(0, self.reservoir.clone())]),
            walls,
            GasTransportPaths::default(),
            &links,
            conduction,
            duration_s,
            &BTreeMap::from([(0, endpoint)]),
            controls.transport,
            &mut cancelled,
        )?;
        let (mut parcels, walls, thermal) = candidate.into_parts();
        let mut momentum = Vec::with_capacity(ordered.len());
        let mut impulse: [Sum; 3] = core::array::from_fn(|_| Sum::default());
        for outlet in &thermal.gas.outlets {
            poll(&mut cancelled)?;
            let authored = ordered[&outlet.id];
            let channel = authored
                .interface
                .vacuum_outlet(outlet, &mut cancelled)
                .map_err(|error| EffusionRocketError::Source {
                    outlet: outlet.id,
                    error,
                })?;
            let axis = authored.outward_normal_world;
            for (sum, direction) in impulse.iter_mut().zip([axis.x, axis.y, axis.z]) {
                sum.add(direction * channel.exterior_impulses_n_s[1])?;
            }
            momentum.push(EffusionChannelMomentumStep {
                outward_normal_world: axis,
                channel,
            });
        }
        let relative_momentum =
            Vec3::new(impulse[0].total(), impulse[1].total(), impulse[2].total());
        let reservoir = parcels.remove(&0).ok_or(EffusionChannelError::Invalid(
            "finite-channel candidate lost its original gas inventory",
        ))?;
        let (endpoint, recoil) = translate(
            self,
            reservoir,
            &thermal,
            relative_momentum,
            environment,
            duration_s,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(EffusionChannelCandidate {
            endpoint,
            walls,
            step: EffusionChannelStep {
                duration_s,
                thermal,
                momentum,
                recoil,
            },
        })
    }
}

// Only accepted escaping moments cross into the original mechanics owner. Its
// explicit inventory/rate rounding follows the existing thin-aperture owner;
// material heat and local channel wall forces are never added as mechanical loads.
#[allow(clippy::too_many_arguments)]
fn translate(
    original: &GasRocketState,
    reservoir: ClosedGasParcel,
    thermal: &FiniteVacuumChannelStep,
    relative_momentum_kg_m_s: Vec3,
    environment: EffusionRocketEnvironment,
    duration_s: f64,
    controls: EffusionRocketSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(GasRocketState, EffusionChannelRecoil), EffusionChannelError> {
    poll(cancelled)?;
    let expelled_mass_kg = thermal.gas.audit.escaped_mass_kg;
    let relative_energy_j = thermal.gas.audit.escaped_energy_j;
    let forcing = ExhaustForcing {
        mass_flow_kg_s: expelled_mass_kg / duration_s,
        relative_momentum_rate_n: Vec3::new(
            relative_momentum_kg_m_s.x / duration_s,
            relative_momentum_kg_m_s.y / duration_s,
            relative_momentum_kg_m_s.z / duration_s,
        ),
        relative_energy_rate_w: relative_energy_j / duration_s,
    };
    let initial_mass_residual_kg =
        original.vehicle.propellant_mass_kg - original.reservoir.mass_kg();
    let mass_budget = controls.transport.exchange.mass_tolerance_kg;
    if !initial_mass_residual_kg.is_finite() || initial_mass_residual_kg.abs() > mass_budget {
        return Err(EffusionRocketError::InventoryBalance {
            remaining_residual_kg: initial_mass_residual_kg,
            consumed_residual_kg: 0.0,
        }
        .into());
    }
    let mut native_initial = original.vehicle;
    native_initial.propellant_mass_kg = original.reservoir.mass_kg();
    let initial_inventory_momentum_correction_n_s = original
        .vehicle
        .velocity_m_s
        .scale(-initial_mass_residual_kg);
    let initial_inventory_kinetic_correction_j = if initial_mass_residual_kg == 0.0 {
        0.0
    } else {
        (-0.5 * initial_mass_residual_kg) * original.vehicle.velocity_m_s.norm_squared()
    };
    let exhaust = native_initial
        .advance_with_exhaust_forced(
            forcing,
            ExhaustEnvironment {
                gravity_m_s2: environment.gravity_world_m_s2,
                external_force_n: environment.external_force_world_n,
            },
            duration_s,
            &mut *cancelled,
        )
        .map_err(EffusionRocketError::from)?;
    let mut vehicle = exhaust.mechanics.state_after;
    let remaining_mass_residual_kg = vehicle.propellant_mass_kg - reservoir.mass_kg();
    let consumed_mass_residual_kg = exhaust.mechanics.burned_propellant_mass_kg - expelled_mass_kg;
    if !remaining_mass_residual_kg.is_finite()
        || !consumed_mass_residual_kg.is_finite()
        || remaining_mass_residual_kg.abs() > mass_budget
        || consumed_mass_residual_kg.abs() > mass_budget
        || (expelled_mass_kg > 0.0 && exhaust.mechanics.coast_duration_s != 0.0)
    {
        return Err(EffusionRocketError::InventoryBalance {
            remaining_residual_kg: remaining_mass_residual_kg,
            consumed_residual_kg: consumed_mass_residual_kg,
        }
        .into());
    }
    let inventory_momentum_correction_n_s = vehicle.velocity_m_s.scale(-remaining_mass_residual_kg);
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
    let source_kinetic_correction_j = impulse_rate_correction.dot(exhaust.mean_velocity_world_m_s)
        + 0.5 * mass_rate_correction * exhaust.mean_speed_squared_m2_s2;
    let exhaust_momentum_world_kg_m_s = exhaust
        .exhaust_momentum_world_kg_m_s
        .add(source_momentum_correction_n_s);
    let exhaust_energy_world_j = (Dd::from_f64(relative_energy_j)
        + Dd::from_f64(exhaust.relative_kinetic_transport_energy_j)
        + Dd::from_f64(source_kinetic_correction_j))
    .to_f64();
    if !exhaust_momentum_world_kg_m_s.is_finite() || !exhaust_energy_world_j.is_finite() {
        return Err(EffusionChannelError::Invalid(
            "unrepresentable accepted inertial channel exhaust",
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
        .any(|v| v.abs() > controls.momentum_tolerance_n_s)
    {
        return Err(EffusionRocketError::MomentumBalance {
            residual_n_s: momentum_residual_kg_m_s,
        }
        .into());
    }
    let energy_terms = [
        thermal.energy_residual_j,
        exhaust.kinetic_energy_residual_j,
        inventory_kinetic_correction_j,
        source_kinetic_correction_j,
    ];
    let energy_residual_j = energy_terms
        .into_iter()
        .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value))
        .to_f64();
    let absolute_energy_residual_sum_j = energy_terms
        .into_iter()
        .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value.abs()))
        .to_f64();
    let gravity_mass_history_uncertainty_kg =
        if expelled_mass_kg == 0.0 || environment.gravity_world_m_s2 == Vec3::ZERO {
            0.0
        } else {
            consumed_mass_residual_kg.abs()
                + exhaust.mass_residual_kg.abs()
                + remaining_mass_residual_kg.abs()
                + thermal.gas.audit.total_mass_residual_kg.abs()
        };
    let (gravity_impulse_uncertainty_n_s, gravity_work_uncertainty_j) = gravity_rounding_bounds(
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
    let energy_uncertainty_bound_j =
        gravity_work_uncertainty_j + initial_inventory_kinetic_correction_j.abs();
    if !energy_residual_j.is_finite()
        || !absolute_energy_residual_sum_j.is_finite()
        || !energy_uncertainty_bound_j.is_finite()
        || energy_uncertainty_bound_j
            > controls.transport.exchange.energy_tolerance_j - absolute_energy_residual_sum_j
    {
        return Err(EffusionChannelError::EnergyBalance {
            residual_j: energy_residual_j,
            absolute_residual_sum_j: absolute_energy_residual_sum_j,
            uncertainty_bound_j: energy_uncertainty_bound_j,
        });
    }
    if !momentum_bound.is_finite()
        || [
            (momentum_residual_kg_m_s.x, momentum_bound.x),
            (momentum_residual_kg_m_s.y, momentum_bound.y),
            (momentum_residual_kg_m_s.z, momentum_bound.z),
        ]
        .into_iter()
        .any(|(r, bound)| bound > controls.momentum_tolerance_n_s - r.abs())
    {
        return Err(EffusionRocketError::RoundingBudget {
            momentum_residual_n_s: momentum_residual_kg_m_s,
            energy_residual_j,
            momentum_bound_n_s: momentum_bound,
            energy_bound_j: energy_uncertainty_bound_j,
        }
        .into());
    }
    vehicle.propellant_mass_kg = reservoir.mass_kg();
    poll(cancelled)?;
    Ok((
        GasRocketState {
            reservoir,
            vehicle,
            inertial_frame: original.inertial_frame,
        },
        EffusionChannelRecoil {
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
            absolute_energy_residual_sum_j,
            energy_uncertainty_bound_j,
        },
    ))
}
