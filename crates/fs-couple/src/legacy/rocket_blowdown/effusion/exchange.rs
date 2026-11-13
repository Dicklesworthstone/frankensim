//! Bidirectional molecular ports on the original finite translating vessel.
//!
//! A source can remain instantaneously COMOVING with the aperture or have a
//! prescribed velocity in the declared inertial frame. Ambient distributions
//! are sampled relative to the actual vessel velocity at each step. Sources
//! supply and receive the reported particles and their inertial momentum/energy;
//! they are not second finite vessels or an exterior wall-pressure law.
//! The original molecular crossing energy enters the gas once; no additional
//! invented thermalization heat or continuum pressure thrust is added.

use super::{
    EffusionOutlet, EffusionRocketEnvironment, EffusionRocketError, EffusionRocketSolve,
    GasRocketState, Sum, gravity_rounding_bounds, poll,
};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::{
            snapshot::MAX_SNAPSHOT_TRANSFERS,
            species::{
                effusion::{
                    drift::{DriftingEffusiveSpeciesFlux, EffusiveDriftInterface},
                    momentum::{EffusiveBeamMomentum, EffusiveMomentumInterface},
                },
                network::transport::{
                    GasTransportPaths,
                    effusion::{
                        reservoir::{
                            self, EffusionDriftingReservoirLink, EffusionReservoirLink,
                            ReservoirEffusionStep,
                        },
                        vacuum::VacuumEffusionLink,
                    },
                },
            },
        },
    },
};
use fs_mbd::{
    Vec3,
    rocket::exhaust::{
        ExhaustEnvironment, ParticleExchangeForcing, ParticleExchangeStep, ParticleStreamIntegral,
        ParticleStreamMoments,
    },
};
use std::collections::BTreeMap;

/// A prescribed thermochemical source on a separately exposed physical aperture.
#[derive(Clone, Copy, Debug)]
pub struct EffusionReservoirPort<'a> {
    pub id: u64,
    /// Thermochemical snapshot held throughout this interval. The external
    /// source follows the vessel velocity; its finite mass is NOT depleted.
    pub source: &'a ClosedGasParcel,
    pub interface: &'a EffusiveMomentumInterface,
    pub area_m2: f64,
    /// Unit vector OUT of the vessel in its declared inertial frame.
    pub outward_normal_world: Vec3,
}

/// A Maxwellian source with fixed thermochemical state and inertial bulk velocity.
/// Its nominal volume does not impose a finite donor inventory. The separate
/// drift qualification admits the complete relative velocity and aperture law.
#[derive(Clone, Copy, Debug)]
pub struct EffusionAmbientPort<'a> {
    pub id: u64,
    pub source: &'a ClosedGasParcel,
    /// Prescribed bulk gas velocity [m/s], fixed during the requested interval.
    pub velocity_world_m_s: Vec3,
    /// Must match the vessel/environment; no coordinate transform is inferred.
    pub inertial_frame: u64,
    pub interface: &'a EffusiveDriftInterface,
    /// Separately exposed area [m²]; zero leaves this source and drift inactive.
    pub area_m2: f64,
    /// Unit vector OUT of the vessel in the same declared inertial frame.
    pub outward_normal_world: Vec3,
}

/// One actual gross particle stream. An incoming relative moment normally
/// points inward; its positive mass and complete energy are never negated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionExchangeStream {
    pub mass_kg: f64,
    pub relative_momentum_kg_m_s: Vec3,
    pub relative_energy_j: f64,
    pub momentum_world_kg_m_s: Vec3,
    pub energy_world_j: f64,
}

/// An external owner can accept the opposite of this source's own receipt.
/// IDs are canonical and shared with the original material boundary account.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionExchangePortStep {
    pub id: u64,
    pub incoming: EffusionExchangeStream,
    pub outgoing: EffusionExchangeStream,
    /// Authoritative signed inventory increment, accumulated before gross sums.
    pub net_mass_into_gas_kg: f64,
    /// Original relative-caloric evaluation of the COMPLETE net molecular energy.
    /// Unlike subtracting gross energies, it retains small thermal counterflow
    /// beside arbitrary large species formation-energy references.
    pub net_relative_energy_into_gas_j: f64,
    pub net_momentum_into_vessel_world_kg_m_s: Vec3,
    pub net_energy_into_vessel_world_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionExchangeStep {
    pub duration_s: f64,
    pub gas: ReservoirEffusionStep,
    /// Raw native held-rate mechanics; the accepted endpoint below uses the
    /// authoritative gas inventory after admitting the explicit correction.
    pub mechanics: ParticleExchangeStep,
    pub incoming: EffusionExchangeStream,
    pub outgoing: EffusionExchangeStream,
    pub ports: Vec<EffusionExchangePortStep>,
    pub net_mass_into_gas_kg: f64,
    pub net_relative_energy_into_gas_j: f64,
    pub net_momentum_into_vessel_world_kg_m_s: Vec3,
    pub net_energy_into_vessel_world_j: f64,
    pub initial_mass_residual_kg: f64,
    pub initial_inventory_momentum_correction_n_s: Vec3,
    pub initial_inventory_kinetic_correction_j: f64,
    pub remaining_mass_residual_kg: f64,
    pub incoming_mass_residual_kg: f64,
    pub outgoing_mass_residual_kg: f64,
    /// Correction to NET OUTWARD native momentum from the original accepted
    /// boundary inventories and moments, including gross-sum rate rounding.
    pub source_momentum_correction_n_s: Vec3,
    pub source_kinetic_correction_j: f64,
    pub inventory_momentum_correction_n_s: Vec3,
    pub inventory_kinetic_correction_j: f64,
    pub gravity_mass_history_uncertainty_kg: f64,
    pub gravity_impulse_uncertainty_n_s: Vec3,
    pub gravity_work_uncertainty_j: f64,
    pub momentum_residual_kg_m_s: Vec3,
    /// Gas source-relative closure plus kinetic/source/inventory closure;
    /// gross formation-energy reports are never subtracted to evaluate it.
    pub energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionExchangeCandidate {
    endpoint: GasRocketState,
    step: EffusionExchangeStep,
}
impl EffusionExchangeCandidate {
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    pub const fn step(&self) -> &EffusionExchangeStep {
        &self.step
    }
    pub fn into_parts(self) -> (GasRocketState, EffusionExchangeStep) {
        (self.endpoint, self.step)
    }
}

#[derive(Default)]
struct StreamSum {
    mass: Sum,
    energy: Sum,
    impulse: [Sum; 3],
}
impl StreamSum {
    fn add(
        &mut self,
        mass: f64,
        beam: &EffusiveBeamMomentum<'_>,
        direction: Vec3,
        id: u64,
    ) -> Result<(), EffusionRocketError> {
        if mass == 0.0 {
            return Ok(());
        }
        let impulse = beam
            .impulse_for_mass_kg(mass)
            .map_err(|error| EffusionRocketError::Source { outlet: id, error })?;
        self.add_moments(
            mass,
            beam.beam().crossing_energy_j_kg,
            [direction.x, direction.y, direction.z].map(|component| impulse * component),
        )
    }
    fn add_drifting(
        &mut self,
        mass: f64,
        beam: &DriftingEffusiveSpeciesFlux,
        id: u64,
    ) -> Result<(), EffusionRocketError> {
        if mass == 0.0 {
            return Ok(());
        }
        let impulse = beam
            .impulse_for_mass_kg(mass)
            .map_err(|error| EffusionRocketError::Source { outlet: id, error })?;
        self.add_moments(mass, beam.incoming.crossing_energy_j_kg, impulse)
    }
    fn add_moments(
        &mut self,
        mass: f64,
        specific: Option<f64>,
        impulse: [f64; 3],
    ) -> Result<(), EffusionRocketError> {
        let specific = specific.ok_or(EffusionRocketError::Invalid(
            "accepted particle stream lost original molecular energy",
        ))?;
        let energy = mass * specific;
        if !energy.is_finite() || (specific != 0.0 && energy == 0.0) {
            return Err(EffusionRocketError::Invalid(
                "accepted gross molecular energy is unrepresentable",
            ));
        }
        self.mass.add(mass)?;
        self.energy.add(energy)?;
        for (sum, component) in self.impulse.iter_mut().zip(impulse) {
            sum.add(component)?;
        }
        Ok(())
    }
    fn add_stream(&mut self, stream: &EffusionExchangeStream) -> Result<(), EffusionRocketError> {
        self.mass.add(stream.mass_kg)?;
        self.energy.add(stream.relative_energy_j)?;
        for (sum, component) in self.impulse.iter_mut().zip([
            stream.relative_momentum_kg_m_s.x,
            stream.relative_momentum_kg_m_s.y,
            stream.relative_momentum_kg_m_s.z,
        ]) {
            sum.add(component)?;
        }
        Ok(())
    }
    fn finish(self) -> EffusionExchangeStream {
        EffusionExchangeStream {
            mass_kg: self.mass.total(),
            relative_momentum_kg_m_s: Vec3::new(
                self.impulse[0].total(),
                self.impulse[1].total(),
                self.impulse[2].total(),
            ),
            relative_energy_j: self.energy.total(),
            momentum_world_kg_m_s: Vec3::ZERO,
            energy_world_j: 0.0,
        }
    }
}

impl EffusionExchangeStream {
    fn rates(&self, duration_s: f64) -> ParticleStreamMoments {
        ParticleStreamMoments {
            mass_flow_kg_s: self.mass_kg / duration_s,
            relative_momentum_rate_n: Vec3::new(
                self.relative_momentum_kg_m_s.x / duration_s,
                self.relative_momentum_kg_m_s.y / duration_s,
                self.relative_momentum_kg_m_s.z / duration_s,
            ),
            relative_energy_rate_w: self.relative_energy_j / duration_s,
        }
    }
    fn frame(&mut self, mean: Vec3, mean_squared: f64) -> Result<(), EffusionRocketError> {
        self.momentum_world_kg_m_s = mean.scale(self.mass_kg).add(self.relative_momentum_kg_m_s);
        self.energy_world_j = self.relative_energy_j
            + self.relative_momentum_kg_m_s.dot(mean)
            + 0.5 * self.mass_kg * mean_squared;
        self.admit_frame()
    }
    fn corrected_frame(
        &mut self,
        native: &ParticleStreamIntegral,
        mean: Vec3,
        mean_squared: f64,
    ) -> Result<(), EffusionRocketError> {
        let dm = self.mass_kg - native.mass_kg;
        let dj = self
            .relative_momentum_kg_m_s
            .sub(native.relative_momentum_n_s);
        self.momentum_world_kg_m_s = native.momentum_world_kg_m_s.add(mean.scale(dm)).add(dj);
        self.energy_world_j = self.relative_energy_j
            + native.kinetic_transport_energy_j
            + dj.dot(mean)
            + 0.5 * dm * mean_squared;
        self.admit_frame()
    }
    fn admit_frame(&self) -> Result<(), EffusionRocketError> {
        if self.momentum_world_kg_m_s.is_finite() && self.energy_world_j.is_finite() {
            Ok(())
        } else {
            Err(EffusionRocketError::Invalid(
                "accepted inertial particle stream is unrepresentable",
            ))
        }
    }
}

// Combine source-coordinate terms before transforming to a possibly boosted
// inertial frame. In particular, gross mass/energy sums never define net heat.
fn net_frame(
    incoming: &EffusionExchangeStream,
    outgoing: &EffusionExchangeStream,
    net_mass: f64,
    net_energy: f64,
    mean: Vec3,
    mean_squared: f64,
) -> Result<(Vec3, f64), EffusionRocketError> {
    let relative = incoming
        .relative_momentum_kg_m_s
        .sub(outgoing.relative_momentum_kg_m_s);
    let momentum = mean.scale(net_mass).add(relative);
    let energy = net_energy + relative.dot(mean) + 0.5 * net_mass * mean_squared;
    if momentum.is_finite() && energy.is_finite() {
        Ok((momentum, energy))
    } else {
        Err(EffusionRocketError::Invalid(
            "net inertial particle exchange is unrepresentable",
        ))
    }
}

impl GasRocketState {
    /// Exchange molecules with prescribed comoving sources and true vacuum
    /// outlets on ONE original species budget and ONE caloric reconstruction.
    ///
    /// Both incoming and outgoing molecules push opposite the outward aperture
    /// normal. Equal opposing masses therefore retain their physical stress.
    /// External wall pressure, supporting apparatus and any source acceleration
    /// remain external physical owners; this method reports their particle
    /// exchange, not a stationary ambient or an implicit pressure-difference force.
    /// All fluxes, normals and loads are held for this first-order finite-gas
    /// interval. Source changes require resampling between accepted calls.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_exchange(
        &self,
        environment: EffusionRocketEnvironment,
        vacuum_outlets: &[EffusionOutlet<'_>],
        reservoir_ports: &[EffusionReservoirPort<'_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        controls: EffusionRocketSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionExchangeCandidate, EffusionRocketError> {
        self.advance_effusion_exchange_with_ambient(
            environment,
            vacuum_outlets,
            reservoir_ports,
            &[],
            duration_s,
            endpoint,
            controls,
            cancelled,
        )
    }

    /// Exchange with vacuum, comoving sources and fixed inertial Maxwellians in
    /// one original gas transaction and one native mechanical trajectory.
    ///
    /// Each ambient's velocity minus the ORIGINAL vessel velocity selects its
    /// incoming half-space mass, full vector momentum and crossing energy. The
    /// native step holds these relative moments for the finite interval; this
    /// is a first-order frozen-moment approximation to the fixed ambient.
    /// Repeated steps or the adaptive owner resample actual relative drift.
    /// Final actual relative velocities must remain inside their original
    /// qualification, without resampling rates or spending inventory twice.
    /// Admission covers sampled and endpoint speeds; it is not an enclosure
    /// of the continuous speed throughout the native held-moment trajectory.
    ///
    /// Tangential momentum and its energy enter through the same particles.
    /// No additional thermalization heat, pressure thrust or exterior wall load
    /// is inferred. All port IDs share one namespace and one withdrawal budget.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_exchange_with_ambient(
        &self,
        environment: EffusionRocketEnvironment,
        vacuum_outlets: &[EffusionOutlet<'_>],
        reservoir_ports: &[EffusionReservoirPort<'_>],
        ambient_ports: &[EffusionAmbientPort<'_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        controls: EffusionRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionExchangeCandidate, EffusionRocketError> {
        poll(&mut cancelled)?;
        if environment.inertial_frame != self.inertial_frame
            || !environment.gravity_world_m_s2.is_finite()
            || !environment.external_force_world_n.is_finite()
            || !controls.momentum_tolerance_n_s.is_finite()
            || controls.momentum_tolerance_n_s < 0.0
            || vacuum_outlets
                .len()
                .checked_add(reservoir_ports.len())
                .and_then(|n| n.checked_add(ambient_ports.len()))
                .is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
        {
            return Err(EffusionRocketError::Invalid(
                "molecular exchange frame, loads, budget or port count",
            ));
        }
        for port in ambient_ports {
            poll(&mut cancelled)?;
            if port.inertial_frame != self.inertial_frame
                || (port.area_m2 != 0.0 && !port.velocity_world_m_s.is_finite())
            {
                return Err(EffusionRocketError::Invalid(
                    "ambient velocity needs the vessel frame and finite active components",
                ));
            }
        }
        let mut ordered = BTreeMap::new();
        for (id, interface, normal) in vacuum_outlets
            .iter()
            .map(|p| (p.id, p.interface, p.outward_normal_world))
            .chain(
                reservoir_ports
                    .iter()
                    .map(|p| (p.id, p.interface, p.outward_normal_world)),
            )
            .chain(ambient_ports.iter().map(|p| {
                (
                    p.id,
                    p.interface.momentum_interface(),
                    p.outward_normal_world,
                )
            }))
        {
            poll(&mut cancelled)?;
            let squared = normal.norm_squared();
            if !normal.is_finite()
                || !squared.is_finite()
                || (squared - 1.0).abs() > 64.0 * f64::EPSILON
                || ordered.insert(id, (interface, normal)).is_some()
            {
                return Err(EffusionRocketError::Invalid(
                    "molecular ports need globally unique IDs and unit outward normals",
                ));
            }
        }
        let parcels = BTreeMap::from([(0, self.reservoir.clone())]);
        let endpoints = BTreeMap::from([(0, endpoint)]);
        let vacuum: Vec<_> = vacuum_outlets
            .iter()
            .map(|p| VacuumEffusionLink {
                id: p.id,
                parcel: 0,
                interface: p.interface.interface(),
                area_m2: p.area_m2,
            })
            .collect();
        let reservoirs: Vec<_> = reservoir_ports
            .iter()
            .map(|p| EffusionReservoirLink {
                id: p.id,
                parcel: 0,
                reservoir: p.source,
                interface: p.interface.interface(),
                area_m2: p.area_m2,
            })
            .collect();
        let drifting: Vec<_> = ambient_ports
            .iter()
            .map(|p| {
                let relative = if p.area_m2 == 0.0 {
                    Vec3::ZERO
                } else {
                    p.velocity_world_m_s.sub(self.vehicle.velocity_m_s)
                };
                let normal = p.outward_normal_world;
                EffusionDriftingReservoirLink {
                    id: p.id,
                    parcel: 0,
                    reservoir: p.source,
                    interface: p.interface,
                    area_m2: p.area_m2,
                    relative_bulk_velocity_m_s: [relative.x, relative.y, relative.z],
                    outward_unit_normal: [normal.x, normal.y, normal.z],
                }
            })
            .collect();
        let candidate = reservoir::step_candidate_with_drift(
            &parcels,
            GasTransportPaths::default(),
            &vacuum,
            &reservoirs,
            &drifting,
            duration_s,
            &endpoints,
            controls.transport,
            &mut cancelled,
        )?;
        let (mut parcels, gas) = candidate.into_parts();
        let mut incoming = StreamSum::default();
        let mut outgoing = StreamSum::default();
        let mut net_mass = Sum::default();
        let mut net_energy = Sum::default();
        let mut ports = Vec::with_capacity(ordered.len());
        let mut boundaries = gas.exchange.boundaries.iter().peekable();
        let mut vacuum = gas.vacuum_outlets.iter().peekable();
        let mut reservoirs = gas.reservoirs.iter().peekable();
        let mut drifting = gas.drifting_reservoirs.iter().peekable();
        for (id, (interface, normal)) in ordered {
            poll(&mut cancelled)?;
            let source = |error| EffusionRocketError::Source { outlet: id, error };
            let vacuum_momentum = if vacuum.peek().is_some_and(|p| p.id == id) {
                Some(
                    interface
                        .vacuum_flux(&vacuum.next().unwrap().flux, &mut cancelled)
                        .map_err(source)?,
                )
            } else {
                None
            };
            let reciprocal_momentum = if reservoirs.peek().is_some_and(|p| p.id == id) {
                Some(
                    interface
                        .reciprocal_flux(&reservoirs.next().unwrap().flux, &mut cancelled)
                        .map_err(source)?,
                )
            } else {
                None
            };
            let drifting_source = if drifting.peek().is_some_and(|p| p.id == id) {
                drifting.next()
            } else {
                None
            };
            let drifting_outgoing = drifting_source
                .map(|row| {
                    interface
                        .vacuum_flux(&row.outgoing, &mut cancelled)
                        .map_err(source)
                })
                .transpose()?;
            if usize::from(vacuum_momentum.is_some())
                + usize::from(reciprocal_momentum.is_some())
                + usize::from(drifting_source.is_some())
                != 1
            {
                return Err(EffusionRocketError::Invalid(
                    "accepted port lost its unique original source",
                ));
            }
            let (mut pin, mut pout) = (StreamSum::default(), StreamSum::default());
            let (mut pmass, mut penergy) = (Sum::default(), Sum::default());
            while boundaries.peek().is_some_and(|row| row.id == id) {
                poll(&mut cancelled)?;
                let row = boundaries.next().unwrap();
                let missing = || {
                    EffusionRocketError::Invalid(
                        "accepted species lost its original molecular beam",
                    )
                };
                if let Some(momentum) = &vacuum_momentum {
                    let index = momentum
                        .species()
                        .binary_search_by(|b| b.species.cmp(&row.transfer.species))
                        .map_err(|_| missing())?;
                    pout.add(
                        row.transfer.outgoing_mass_kg,
                        &momentum.species()[index].outward,
                        normal,
                        id,
                    )?;
                    if row.transfer.incoming_mass_kg != 0.0 {
                        return Err(EffusionRocketError::Invalid(
                            "vacuum supplied incoming molecules",
                        ));
                    }
                } else if let Some(momentum) = &reciprocal_momentum {
                    let index = momentum
                        .species()
                        .binary_search_by(|b| b.species.cmp(&row.transfer.species))
                        .map_err(|_| missing())?;
                    let beam = &momentum.species()[index];
                    pin.add(
                        row.transfer.incoming_mass_kg,
                        &beam.right_to_left,
                        normal.scale(-1.0),
                        id,
                    )?;
                    pout.add(
                        row.transfer.outgoing_mass_kg,
                        &beam.left_to_right,
                        normal,
                        id,
                    )?;
                } else if let (Some(ambient), Some(momentum)) =
                    (drifting_source, &drifting_outgoing)
                {
                    // The source and vessel can have distinct active species.
                    // Each accepted direction uses its own original sample.
                    if row.transfer.incoming_mass_kg > 0.0 {
                        let index = ambient
                            .incoming
                            .species
                            .binary_search_by(|b| b.species.cmp(&row.transfer.species))
                            .map_err(|_| missing())?;
                        pin.add_drifting(
                            row.transfer.incoming_mass_kg,
                            &ambient.incoming.species[index],
                            id,
                        )?;
                    }
                    if row.transfer.outgoing_mass_kg > 0.0 {
                        let index = momentum
                            .species()
                            .binary_search_by(|b| b.species.cmp(&row.transfer.species))
                            .map_err(|_| missing())?;
                        pout.add(
                            row.transfer.outgoing_mass_kg,
                            &momentum.species()[index].outward,
                            normal,
                            id,
                        )?;
                    }
                }
                pmass.add(row.transfer.net_mass_into_gas_kg)?;
                penergy.add(row.energy_into_gas_j)?;
            }
            let (pin, pout) = (pin.finish(), pout.finish());
            incoming.add_stream(&pin)?;
            outgoing.add_stream(&pout)?;
            net_mass.add(pmass.total())?;
            net_energy.add(penergy.total())?;
            ports.push(EffusionExchangePortStep {
                id,
                incoming: pin,
                outgoing: pout,
                net_mass_into_gas_kg: pmass.total(),
                net_relative_energy_into_gas_j: penergy.total(),
                net_momentum_into_vessel_world_kg_m_s: Vec3::ZERO,
                net_energy_into_vessel_world_j: 0.0,
            });
        }
        if boundaries.next().is_some()
            || vacuum.next().is_some()
            || reservoirs.next().is_some()
            || drifting.next().is_some()
        {
            return Err(EffusionRocketError::Invalid(
                "accepted molecular exchange lost its physical port",
            ));
        }
        let (mut incoming, mut outgoing) = (incoming.finish(), outgoing.finish());
        let net_mass_into_gas_kg = net_mass.total();
        let net_relative_energy_into_gas_j = net_energy.total();
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
        let mechanics = native_initial.advance_with_particle_exchange(
            ParticleExchangeForcing {
                incoming: incoming.rates(duration_s),
                outgoing: outgoing.rates(duration_s),
            },
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
        let mut vehicle = mechanics.state_after;
        let remaining_mass_residual_kg = vehicle.propellant_mass_kg - reservoir.mass_kg();
        let incoming_mass_residual_kg = mechanics.incoming.mass_kg - incoming.mass_kg;
        let outgoing_mass_residual_kg = mechanics.outgoing.mass_kg - outgoing.mass_kg;
        if [
            remaining_mass_residual_kg,
            incoming_mass_residual_kg,
            outgoing_mass_residual_kg,
        ]
        .into_iter()
        .any(|r| !r.is_finite() || r.abs() > mass_budget)
        {
            return Err(EffusionRocketError::StreamInventoryBalance {
                remaining_residual_kg: remaining_mass_residual_kg,
                incoming_residual_kg: incoming_mass_residual_kg,
                outgoing_residual_kg: outgoing_mass_residual_kg,
            });
        }
        let mean = mechanics.mean_velocity_world_m_s;
        let mean_squared = mechanics.mean_speed_squared_m2_s2;
        incoming.corrected_frame(&mechanics.incoming, mean, mean_squared)?;
        outgoing.corrected_frame(&mechanics.outgoing, mean, mean_squared)?;
        for port in &mut ports {
            poll(&mut cancelled)?;
            port.incoming.frame(mean, mean_squared)?;
            port.outgoing.frame(mean, mean_squared)?;
            (
                port.net_momentum_into_vessel_world_kg_m_s,
                port.net_energy_into_vessel_world_j,
            ) = net_frame(
                &port.incoming,
                &port.outgoing,
                port.net_mass_into_gas_kg,
                port.net_relative_energy_into_gas_j,
                mean,
                mean_squared,
            )?;
        }
        let (net_momentum_into_vessel_world_kg_m_s, net_energy_into_vessel_world_j) = net_frame(
            &incoming,
            &outgoing,
            net_mass_into_gas_kg,
            net_relative_energy_into_gas_j,
            mean,
            mean_squared,
        )?;
        // Use the original boundary's signed inventory before gross sums and
        // division/multiplication can erase a small counterflow difference.
        let dm = -net_mass_into_gas_kg - (mechanics.outgoing.mass_kg - mechanics.incoming.mass_kg);
        let dj = outgoing
            .relative_momentum_kg_m_s
            .sub(incoming.relative_momentum_kg_m_s)
            .sub(
                mechanics
                    .outgoing
                    .relative_momentum_n_s
                    .sub(mechanics.incoming.relative_momentum_n_s),
            );
        let source_momentum_correction_n_s = mean.scale(dm).add(dj);
        let source_kinetic_correction_j = dj.dot(mean) + 0.5 * dm * mean_squared;
        let inventory_momentum_correction_n_s =
            vehicle.velocity_m_s.scale(-remaining_mass_residual_kg);
        let inventory_kinetic_correction_j = if remaining_mass_residual_kg == 0.0 {
            0.0
        } else {
            (-0.5 * remaining_mass_residual_kg) * vehicle.velocity_m_s.norm_squared()
        };
        let momentum_residual_kg_m_s = mechanics
            .momentum_residual_kg_m_s
            .add(source_momentum_correction_n_s)
            .add(inventory_momentum_correction_n_s);
        if !momentum_residual_kg_m_s.is_finite()
            || [
                momentum_residual_kg_m_s.x,
                momentum_residual_kg_m_s.y,
                momentum_residual_kg_m_s.z,
            ]
            .into_iter()
            .any(|r| r.abs() > controls.momentum_tolerance_n_s)
        {
            return Err(EffusionRocketError::MomentumBalance {
                residual_n_s: momentum_residual_kg_m_s,
            });
        }
        let energy_residual_j = gas.exchange.total_energy_residual_j
            + mechanics.kinetic_energy_residual_j
            + source_kinetic_correction_j
            + inventory_kinetic_correction_j;
        if !energy_residual_j.is_finite()
            || energy_residual_j.abs() > controls.transport.exchange.energy_tolerance_j
        {
            return Err(EffusionRocketError::EnergyBalance {
                residual_j: energy_residual_j,
            });
        }
        let gravity_mass_history_uncertainty_kg = if (incoming.mass_kg == 0.0
            && outgoing.mass_kg == 0.0)
            || environment.gravity_world_m_s2 == Vec3::ZERO
        {
            0.0
        } else {
            incoming_mass_residual_kg.abs()
                + outgoing_mass_residual_kg.abs()
                + mechanics.mass_residual_kg.abs()
                + remaining_mass_residual_kg.abs()
                + gas.exchange.total_mass_residual_kg.abs()
                + dm.abs()
        };
        let (gravity_impulse_uncertainty_n_s, gravity_work_uncertainty_j) =
            gravity_rounding_bounds(
                gravity_mass_history_uncertainty_kg,
                environment.gravity_world_m_s2,
                duration_s,
                mean_squared,
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
            .any(|(r, b)| b > controls.momentum_tolerance_n_s - r.abs())
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
        for port in ambient_ports {
            poll(&mut cancelled)?;
            if port.area_m2 > 0.0 {
                let relative = port.velocity_world_m_s.sub(vehicle.velocity_m_s);
                port.interface
                    .admit_relative_velocity([relative.x, relative.y, relative.z])
                    .map_err(|error| EffusionRocketError::Source {
                        outlet: port.id,
                        error,
                    })?;
            }
        }
        vehicle.propellant_mass_kg = reservoir.mass_kg();
        poll(&mut cancelled)?;
        Ok(EffusionExchangeCandidate {
            endpoint: Self {
                reservoir,
                vehicle,
                inertial_frame: self.inertial_frame,
            },
            step: EffusionExchangeStep {
                duration_s,
                gas,
                mechanics,
                incoming,
                outgoing,
                ports,
                net_mass_into_gas_kg,
                net_relative_energy_into_gas_j,
                net_momentum_into_vessel_world_kg_m_s,
                net_energy_into_vessel_world_j,
                initial_mass_residual_kg,
                initial_inventory_momentum_correction_n_s,
                initial_inventory_kinetic_correction_j,
                remaining_mass_residual_kg,
                incoming_mass_residual_kg,
                outgoing_mass_residual_kg,
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
