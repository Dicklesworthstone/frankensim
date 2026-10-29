//! Conservative, pressure-coupled networks of original-source two-zone rooms.
//!
//! The original stratified opening owner computes every directional slab;
//! fs-gas transports those species and enthalpies and accounts for moving
//! layer-interface work. A bounded Newton solve makes opening pressure implicit
//! so low-Mach compartment compressibility does not impose an acoustic-size
//! explicit step. Temperature, composition and interface geometry in the vent
//! flux are frozen at the beginning of a step. This is a first-order estimated
//! zone-network integrator, not a field fire or general compressible duct solver.

use crate::stratified::{self, TwoLayerColumn, VentError, WallOpening};
use fs_gas::multicomponent::MAX_SPECIES;
use fs_gas::multicomponent::compartment::{Layer, LayerRates, TwoZone, ZoneError, ZoneSolve};
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_qty::{Acceleration, Density, Length, MassFlowRate, Power, Pressure, Temperature, Time};

/// Original enthalpy-kernel, finite-volume wall conduction and gas contacts.
pub mod wall;

/// Original-source plume entrainment and explicitly prescribed heat release.
pub mod fire;

#[cfg(test)]
mod tests;

/// Maximum number of coupled rooms in the bounded dense pressure solve.
pub const MAX_ROOMS: usize = 32;
/// Maximum number of authored wall and mechanical connections in total.
pub const MAX_CONNECTIONS: usize = 512;

/// One room and its floor elevation in the shared vertical coordinate system.
#[derive(Clone, Debug, PartialEq)]
pub struct Room {
    /// Source-owned gas inventories and moving interface.
    pub gas: TwoZone,
    /// Floor height relative to the common vent datum.
    pub floor_elevation: Length,
}

/// Source-matched, fixed ambient reservoir. Fractions are in the rooms' common
/// original species order. No default atmosphere or oxygen content is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct Ambient {
    /// Absolute pressure at the supplied reference height.
    pub pressure: Pressure,
    /// Height at which pressure is declared.
    pub reference_height: Length,
    /// Reservoir temperature.
    pub temperature: Temperature,
    /// Canonical mass fractions.
    pub mass_fractions: Vec<f64>,
}

/// One endpoint of an opening or prescribed mechanical transfer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Endpoint {
    /// An index into the original ordered room list.
    Room(usize),
    /// The declared external reservoir.
    Ambient,
}

/// A bidirectional wall opening between rooms or between one room and ambient.
#[derive(Clone, Debug, PartialEq)]
pub struct WallVent {
    /// Left side of the original signed opening law.
    pub left: Endpoint,
    /// Right side of the original signed opening law.
    pub right: Endpoint,
    /// World-coordinate rectangular opening and source-bearing discharge law.
    pub opening: WallOpening,
}

/// Prescribed mechanical ventilation. The caller explicitly supplies any fan
/// work delivered to the receiving gas; zero means no deposited mechanical heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MechanicalVent {
    /// Donor compartment or ambient.
    pub from: Endpoint,
    /// Donor layer (ignored for a homogeneous ambient).
    pub from_layer: Layer,
    /// Receiver compartment or ambient.
    pub to: Endpoint,
    /// Receiver layer (ignored for ambient).
    pub to_layer: Layer,
    /// Nonnegative directed mass transfer.
    pub mass_rate: MassFlowRate,
    /// External mechanical power delivered to the receiver, nonnegative.
    pub deposited_power: Power,
}

/// Explicit nonlinear, source-domain, and duration budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkControl {
    /// Maximum Newton iterations of the coupled pressure endpoint equation.
    pub maximum_iterations: usize,
    /// Maximum whole-step halvings after a failed nonlinear/source attempt.
    pub maximum_halvings: usize,
    /// Smallest admitted attempted duration [s].
    pub minimum_step_s: f64,
    /// Maximum absolute pressure-equation residual [Pa].
    pub pressure_tolerance_pa: f64,
    /// Maximum nominal forward pressure difference [Pa]. It shrinks with the
    /// pressure residual, with a representable-increment floor at the trial pressure.
    pub pressure_difference_pa: f64,
    /// Original per-room moving-interface solve controls.
    pub zone: ZoneSolve,
}
impl Default for NetworkControl {
    fn default() -> Self {
        Self {
            maximum_iterations: 24,
            maximum_halvings: 20,
            minimum_step_s: 1e-8,
            pressure_tolerance_pa: 1e-5,
            pressure_difference_pa: 0.01,
            zone: ZoneSolve::default(),
        }
    }
}

/// Balance of one fully accepted network step. Every room advances exactly the
/// same duration; all report values exclude discarded pressure trials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkStepReport {
    /// Common accepted physical duration.
    pub elapsed: Time,
    /// Whole-step duration reductions.
    pub halvings: usize,
    /// Accepted pressure iterations.
    pub iterations: usize,
    /// Largest pressure equation residual [Pa] before an optional split source.
    pub maximum_pressure_residual_pa: f64,
    /// Largest pressure change [Pa] caused by an optional split source. This is
    /// not part of the preceding transport pressure-equation residual.
    pub maximum_split_pressure_change_pa: f64,
    /// Integrated net authored source and ambient exchange in canonical
    /// species order [kg], including any declared split-reaction sources.
    pub external_species_change_kg: [f64; MAX_SPECIES],
    /// Integrated net physical energy supplied to gas [J], including ambient
    /// enthalpy, authored heat, wall exchange and deposited mechanical power.
    pub external_physical_energy_j: f64,
    /// Species residuals against authored sources, ambient exchange and fans [kg].
    pub species_residual_kg: [f64; MAX_SPECIES],
    /// Whole-network mass residual [kg].
    pub mass_residual_kg: f64,
    /// Whole-network physical energy residual [J].
    pub physical_energy_residual_j: f64,
}

/// No public room is modified on any failed attempt or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum NetworkError {
    /// Invalid dimensions, sources, room identity, controls, or geometry.
    Invalid(&'static str),
    /// An original two-zone owner refused an endpoint.
    Zone(ZoneError),
    /// The original opening owner refused a profile or flow.
    Vent(VentError),
    /// The shared prepared dense solver refused a pressure Jacobian.
    Linear(LuWorkspaceError),
    /// The original wall enthalpy/conduction adapter refused a candidate.
    Wall(wall::WallError),
    /// The original prescribed-fire/plume source refused its regime.
    Fire(fire::PrescribedFireError),
    /// The original inventory-limited reaction adapter refused a split source.
    Burning(fs_gas::multicomponent::compartment::burning::BurningError),
    /// Pressure closure exhausted its bounded nonlinear budget.
    PressureIterations,
    /// All bounded complete-network attempts failed.
    RetryLimit {
        attempts: usize,
        last: Box<NetworkError>,
    },
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for NetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "compartment network: {self:?}")
    }
}
impl std::error::Error for NetworkError {}
impl From<ZoneError> for NetworkError {
    fn from(e: ZoneError) -> Self {
        match e {
            ZoneError::Cancelled => Self::Cancelled,
            other => Self::Zone(other),
        }
    }
}
impl From<VentError> for NetworkError {
    fn from(e: VentError) -> Self {
        match e {
            VentError::Cancelled => Self::Cancelled,
            other => Self::Vent(other),
        }
    }
}
impl From<LuWorkspaceError> for NetworkError {
    fn from(e: LuWorkspaceError) -> Self {
        Self::Linear(e)
    }
}
impl From<wall::WallError> for NetworkError {
    fn from(e: wall::WallError) -> Self {
        match e {
            wall::WallError::Cancelled => Self::Cancelled,
            other => Self::Wall(other),
        }
    }
}
impl From<fire::PrescribedFireError> for NetworkError {
    fn from(e: fire::PrescribedFireError) -> Self {
        match e {
            fire::PrescribedFireError::Cancelled => Self::Cancelled,
            other => Self::Fire(other),
        }
    }
}
impl From<fs_gas::multicomponent::compartment::burning::BurningError> for NetworkError {
    fn from(e: fs_gas::multicomponent::compartment::burning::BurningError) -> Self {
        use fs_gas::multicomponent::compartment::burning::BurningError;
        match e {
            BurningError::Cancelled => Self::Cancelled,
            other => Self::Burning(other),
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), NetworkError> {
    if cancelled() {
        Err(NetworkError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn max_abs(x: &[f64]) -> f64 {
    x.iter().fold(0.0_f64, |a, b| a.max(b.abs()))
}

/// Reusable multi-room owner, including pressure Jacobian/LU storage and a
/// complete unpublished candidate. Runtime pressure evaluations allocate no
/// successful-path vectors, matrices, gas banks, or flow slabs.
pub struct ZoneNetwork {
    rooms: Vec<Room>,
    trial: Vec<Room>,
    ambient: Ambient,
    ambient_density: f64,
    walls: Vec<WallVent>,
    mechanical: Vec<MechanicalVent>,
    gravity: Acceleration,
    rates: Vec<[LayerRates; 2]>,
    external: LayerRates,
    pressures: Vec<f64>,
    perturbation: Vec<f64>,
    residual: Vec<f64>,
    evaluated: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    delta: Vec<f64>,
    solver: LuWorkspace,
}
impl ZoneNetwork {
    /// Admit a bounded topology over identical original species/caloric banks.
    pub fn new(
        rooms: Vec<Room>,
        ambient: Ambient,
        walls: Vec<WallVent>,
        mechanical: Vec<MechanicalVent>,
        gravity: Acceleration,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, NetworkError> {
        poll(&mut cancelled)?;
        let n = rooms.len();
        if !(1..=MAX_ROOMS).contains(&n)
            || walls.len() + mechanical.len() > MAX_CONNECTIONS
            || !gravity.value().is_finite()
            || gravity.value() < 0.0
            || !positive(ambient.pressure.value())
            || !ambient.reference_height.value().is_finite()
        {
            return Err(NetworkError::Invalid(
                "network size, gravity or ambient pressure",
            ));
        }
        for room in &rooms {
            if !room.floor_elevation.value().is_finite()
                || room.gas.gas() != rooms[0].gas.gas()
                || room.gas.time() != rooms[0].gas.time()
            {
                return Err(NetworkError::Invalid(
                    "room source, clock or elevation mismatch",
                ));
            }
        }
        rooms[0].gas.inflow(
            MassFlowRate::new(1.0),
            &ambient.mass_fractions,
            ambient.temperature,
        )?;
        let inverse_mass: f64 = ambient
            .mass_fractions
            .iter()
            .zip(rooms[0].gas.gas().source().models())
            .map(|(y, m)| y / m.molar_mass().value())
            .sum();
        let ambient_density = ambient.pressure.value()
            / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * ambient.temperature.value()
                * inverse_mass);
        if !positive(ambient_density) {
            return Err(NetworkError::Invalid("ambient density"));
        }
        for vent in &walls {
            admit_wall(&rooms, vent)?;
        }
        for vent in &mechanical {
            admit_endpoints(n, vent.from, vent.to)?;
            if !vent.mass_rate.value().is_finite()
                || vent.mass_rate.value() < 0.0
                || !vent.deposited_power.value().is_finite()
                || vent.deposited_power.value() < 0.0
            {
                return Err(NetworkError::Invalid(
                    "mechanical mass rate or deposited power",
                ));
            }
        }
        let trial = rooms.clone();
        Ok(Self {
            rooms,
            trial,
            ambient,
            ambient_density,
            walls,
            mechanical,
            gravity,
            rates: vec![[LayerRates::default(); 2]; n],
            external: LayerRates::default(),
            pressures: vec![0.0; n],
            perturbation: vec![0.0; n],
            residual: vec![0.0; n],
            evaluated: vec![0.0; n],
            matrix: vec![0.0; n * n],
            rhs: vec![0.0; n],
            delta: vec![0.0; n],
            solver: LuWorkspace::new(n)?,
        })
    }
    /// Accepted compartments in original topology order.
    pub fn rooms(&self) -> &[Room] {
        &self.rooms
    }
    /// Original ambient reservoir declaration.
    pub const fn ambient(&self) -> &Ambient {
        &self.ambient
    }
    /// Change one authored opening between steps, including exact zero-area
    /// closure. Existing state, other connections and original loss laws persist.
    pub fn set_wall_opening(
        &mut self,
        index: usize,
        opening: WallOpening,
    ) -> Result<(), NetworkError> {
        let old = self
            .walls
            .get(index)
            .ok_or(NetworkError::Invalid("wall index"))?;
        let next = WallVent {
            left: old.left,
            right: old.right,
            opening,
        };
        admit_wall(&self.rooms, &next)?;
        self.walls[index] = next;
        Ok(())
    }

    /// Advance with explicit source rates frozen over the accepted duration.
    /// Opening pressure is solved implicitly; all other supplied source/closure
    /// states must be refreshed by the caller after each accepted step.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        sources: &[[LayerRates; 2]],
        controls: NetworkControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, NetworkError> {
        poll(&mut cancelled)?;
        if sources.len() != self.rooms.len() {
            return Err(NetworkError::Invalid("network source axis"));
        }
        self.step_with_sources(
            maximum_duration,
            controls,
            |_, _, output, _| {
                output.copy_from_slice(sources);
                Ok(())
            },
            cancelled,
        )
        .map(|(report, ())| report)
    }

    /// Recompute duration-dependent sources from accepted rooms once per
    /// attempted interval. The provider returns an unpublished payload, such
    /// as wall candidates. It must not modify an external physical owner.
    /// Failed attempts drop their payload; only the accepted payload is returned
    /// alongside the atomically accepted network. The caller may then install
    /// that payload by an infallible move, without another cancellation point.
    pub fn step_with_sources<T>(
        &mut self,
        maximum_duration: Time,
        controls: NetworkControl,
        prepare: impl FnMut(
            &[Room],
            Time,
            &mut [[LayerRates; 2]],
            &mut dyn FnMut() -> bool,
        ) -> Result<T, NetworkError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(NetworkStepReport, T), NetworkError> {
        self.step_with_stages(
            maximum_duration,
            controls,
            prepare,
            |_, _, _, _, _| Ok(()),
            cancelled,
        )
    }

    /// Compose transport and a subsequent source operator on one unpublished
    /// candidate. The finish provider sees transported inventories, so reactions
    /// can consume oxygen remaining after outflow. It must preserve identity,
    /// geometry and clocks and add its actual aggregate species/energy source
    /// rates to additional for the balance report. Use the reaction owner's
    /// apply_split to avoid counting elapsed physical time twice.
    ///
    /// Both providers run again after every whole-step halving. Their payloads
    /// and all rooms remain unpublished on any failure. A source stage changes
    /// pressure after the implicit vent solve; the report distinguishes that
    /// change from the transport pressure residual. This is a first-order
    /// operator composition, without a second-order or entropy certificate.
    pub fn step_with_stages<T>(
        &mut self,
        maximum_duration: Time,
        controls: NetworkControl,
        mut prepare: impl FnMut(
            &[Room],
            Time,
            &mut [[LayerRates; 2]],
            &mut dyn FnMut() -> bool,
        ) -> Result<T, NetworkError>,
        mut finish: impl FnMut(
            &mut [Room],
            Time,
            &mut T,
            &mut LayerRates,
            &mut dyn FnMut() -> bool,
        ) -> Result<(), NetworkError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(NetworkStepReport, T), NetworkError> {
        poll(&mut cancelled)?;
        let n = self.rooms.len();
        let mut dt = maximum_duration.value();
        if !positive(dt)
            || !(1..=128).contains(&controls.maximum_iterations)
            || controls.maximum_halvings > 32
            || !positive(controls.minimum_step_s)
            || !positive(controls.pressure_tolerance_pa)
            || !positive(controls.pressure_difference_pa)
        {
            return Err(NetworkError::Invalid("network duration or controls"));
        }
        let mut last = NetworkError::PressureIterations;
        let mut attempts = 0;
        // Bounded stack storage, cleared before each duration-dependent source.
        let mut source_storage = [[LayerRates::default(); 2]; MAX_ROOMS];
        for halving in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            if dt < controls.minimum_step_s {
                break;
            }
            attempts += 1;
            let sources = &mut source_storage[..n];
            sources.fill([LayerRates::default(); 2]);
            for i in 0..n {
                self.pressures[i] = self.rooms[i].gas.pressure().value();
            }
            let candidate: Result<_, NetworkError> = (|| {
                let duration = Time::new(dt);
                let mut payload = prepare(&self.rooms, duration, sources, &mut cancelled)?;
                let (iterations, residual) = self.attempt(dt, sources, controls, &mut cancelled)?;
                let mut before_split = [0.0; MAX_ROOMS];
                for (i, room) in self.trial.iter().enumerate() {
                    before_split[i] = room.gas.pressure().value();
                }
                let mut additional = LayerRates::default();
                finish(
                    &mut self.trial,
                    duration,
                    &mut payload,
                    &mut additional,
                    &mut cancelled,
                )?;
                let mut split_pressure_change = 0.0_f64;
                for (i, room) in self.trial.iter().enumerate() {
                    let original = &self.rooms[i];
                    if room.floor_elevation != original.floor_elevation
                        || room.gas.gas() != original.gas.gas()
                        || room.gas.floor_area() != original.gas.floor_area()
                        || room.gas.height() != original.gas.height()
                        || room.gas.time().value() != original.gas.time().value() + dt
                    {
                        return Err(NetworkError::Invalid(
                            "split source changed room identity or clock",
                        ));
                    }
                    split_pressure_change = split_pressure_change
                        .max((room.gas.pressure().value() - before_split[i]).abs());
                }
                if additional.species_kg_s[self.rooms[0].gas.gas().species_count()..]
                    .iter()
                    .any(|v| *v != 0.0)
                {
                    return Err(NetworkError::Invalid("split source species axis"));
                }
                self.external.add(&additional)?;
                Ok((iterations, residual, split_pressure_change, payload))
            })();
            match candidate {
                Ok((iterations, residual, split_pressure_change, payload)) => {
                    let ns = self.rooms[0].gas.gas().species_count();
                    let mut species = [0.0; MAX_SPECIES];
                    let mut energy = 0.0;
                    for i in 0..n {
                        energy += self.trial[i].gas.sensible_internal_energy().value()
                            - self.rooms[i].gas.sensible_internal_energy().value();
                        for layer in [Layer::Lower, Layer::Upper] {
                            for s in 0..ns {
                                species[s] += self.trial[i].gas.layer(layer).species_masses()[s]
                                    - self.rooms[i].gas.layer(layer).species_masses()[s];
                            }
                        }
                    }
                    energy -= dt * self.external.referenced_energy_w;
                    for s in 0..ns {
                        species[s] -= dt * self.external.species_kg_s[s];
                        energy += species[s] * self.rooms[0].gas.gas().energy_references()[s];
                    }
                    let mass = species[..ns].iter().sum();
                    let mut external_species = [0.0; MAX_SPECIES];
                    let mut external_energy = self.external.referenced_energy_w;
                    for s in 0..ns {
                        external_species[s] = dt * self.external.species_kg_s[s];
                        external_energy += self.external.species_kg_s[s]
                            * self.rooms[0].gas.gas().energy_references()[s];
                    }
                    external_energy *= dt;
                    if !energy.is_finite()
                        || !external_energy.is_finite()
                        || species[..ns].iter().any(|v| !v.is_finite())
                    {
                        return Err(NetworkError::Invalid("network extensive balance"));
                    }
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.rooms, &mut self.trial);
                    return Ok((
                        NetworkStepReport {
                            elapsed: Time::new(dt),
                            halvings: halving,
                            iterations,
                            maximum_pressure_residual_pa: residual,
                            maximum_split_pressure_change_pa: split_pressure_change,
                            external_species_change_kg: external_species,
                            external_physical_energy_j: external_energy,
                            species_residual_kg: species,
                            mass_residual_kg: mass,
                            physical_energy_residual_j: energy,
                        },
                        payload,
                    ));
                }
                Err(NetworkError::Cancelled) => return Err(NetworkError::Cancelled),
                Err(error) => {
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(NetworkError::RetryLimit {
            attempts,
            last: Box::new(last),
        })
    }

    fn attempt(
        &mut self,
        dt: f64,
        sources: &[[LayerRates; 2]],
        controls: NetworkControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(usize, f64), NetworkError> {
        let n = self.rooms.len();
        for iteration in 0..controls.maximum_iterations {
            poll(cancelled)?;
            self.perturbation.copy_from_slice(&self.pressures);
            self.evaluate(dt, sources, controls.zone, cancelled)?;
            self.residual.copy_from_slice(&self.evaluated);
            let norm = max_abs(&self.residual);
            if norm <= controls.pressure_tolerance_pa {
                return Ok((iteration + 1, norm));
            }
            for j in 0..n {
                poll(cancelled)?;
                self.perturbation.copy_from_slice(&self.pressures);
                // A fixed pressure difference can span a nearby neutral-plane
                // transition and leave a biased Jacobian even at convergence.
                // Tighten the linearization with the actual closure residual,
                // while retaining a representable difference at the pressure scale.
                let difference = controls
                    .pressure_difference_pa
                    .min(0.1 * norm)
                    .max(64.0 * f64::EPSILON * self.pressures[j].abs());
                self.perturbation[j] += difference;
                let actual_difference = self.perturbation[j] - self.pressures[j];
                if !positive(actual_difference) {
                    return Err(NetworkError::Invalid("unrepresentable pressure difference"));
                }
                self.evaluate(dt, sources, controls.zone, cancelled)?;
                for i in 0..n {
                    self.matrix[i * n + j] =
                        (self.evaluated[i] - self.residual[i]) / actual_difference;
                }
            }
            for i in 0..n {
                self.rhs[i] = -self.residual[i];
            }
            poll(cancelled)?;
            self.solver
                .solve_into(&self.matrix, &self.rhs, &mut self.delta)?;
            let mut scale = 1.0;
            let mut accepted = false;
            for _ in 0..24 {
                poll(cancelled)?;
                for i in 0..n {
                    self.perturbation[i] = self.pressures[i] + scale * self.delta[i];
                }
                if self.perturbation.iter().all(|p| positive(*p)) {
                    match self.evaluate(dt, sources, controls.zone, cancelled) {
                        Ok(()) if max_abs(&self.evaluated) < norm => {
                            accepted = true;
                            break;
                        }
                        Err(NetworkError::Cancelled) => return Err(NetworkError::Cancelled),
                        _ => {}
                    }
                }
                scale *= 0.5;
            }
            if !accepted {
                return Err(NetworkError::PressureIterations);
            }
            self.pressures.copy_from_slice(&self.perturbation);
        }
        Err(NetworkError::PressureIterations)
    }

    fn evaluate(
        &mut self,
        dt: f64,
        sources: &[[LayerRates; 2]],
        zone: ZoneSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), NetworkError> {
        self.rates.copy_from_slice(sources);
        self.external = LayerRates::default();
        for source in sources {
            self.external.add(&source[0])?;
            self.external.add(&source[1])?;
        }
        for index in 0..self.walls.len() {
            poll(cancelled)?;
            let vent = &self.walls[index];
            let left = column(
                &self.rooms,
                &self.ambient,
                self.ambient_density,
                &self.perturbation,
                vent.left,
            )?;
            let right = column(
                &self.rooms,
                &self.ambient,
                self.ambient_density,
                &self.perturbation,
                vent.right,
            )?;
            let flow = stratified::wall_vent_flow(
                &vent.opening,
                left,
                right,
                self.gravity,
                &mut *cancelled,
            )?;
            for slab in flow.slabs() {
                let (from, from_layer, to, to_layer, rate) = if slab.mass_rate.value() >= 0.0 {
                    (
                        vent.left,
                        gas_layer(slab.left_layer),
                        vent.right,
                        gas_layer(slab.right_layer),
                        slab.mass_rate,
                    )
                } else {
                    (
                        vent.right,
                        gas_layer(slab.right_layer),
                        vent.left,
                        gas_layer(slab.left_layer),
                        MassFlowRate::new(-slab.mass_rate.value()),
                    )
                };
                let carried = carried(&self.rooms, &self.ambient, from, from_layer, rate)?;
                transfer(
                    &mut self.rates,
                    &mut self.external,
                    from,
                    from_layer,
                    to,
                    to_layer,
                    carried,
                )?;
            }
        }
        for vent in &self.mechanical {
            poll(cancelled)?;
            if vent.mass_rate.value() > 0.0 {
                let carried = carried(
                    &self.rooms,
                    &self.ambient,
                    vent.from,
                    vent.from_layer,
                    vent.mass_rate,
                )?;
                transfer(
                    &mut self.rates,
                    &mut self.external,
                    vent.from,
                    vent.from_layer,
                    vent.to,
                    vent.to_layer,
                    carried,
                )?;
            }
            if let Endpoint::Room(i) = vent.to {
                self.rates[i][vent.to_layer.index()].add_heat(vent.deposited_power)?;
                self.external.add_heat(vent.deposited_power)?;
            }
        }
        for i in 0..self.rooms.len() {
            poll(cancelled)?;
            self.trial[i].gas.copy_state_from(&self.rooms[i].gas)?;
            self.trial[i]
                .gas
                .advance(Time::new(dt), self.rates[i], zone, &mut *cancelled)?;
            self.evaluated[i] = self.trial[i].gas.pressure().value() - self.perturbation[i];
        }
        Ok(())
    }
}

fn gas_layer(layer: stratified::Layer) -> Layer {
    match layer {
        stratified::Layer::Lower => Layer::Lower,
        stratified::Layer::Upper => Layer::Upper,
    }
}
fn admit_endpoints(n: usize, left: Endpoint, right: Endpoint) -> Result<(), NetworkError> {
    if left == right {
        return Err(NetworkError::Invalid(
            "a connection must join distinct endpoints",
        ));
    }
    for endpoint in [left, right] {
        if let Endpoint::Room(i) = endpoint {
            if i >= n {
                return Err(NetworkError::Invalid("room endpoint index"));
            }
        }
    }
    Ok(())
}
fn admit_wall(rooms: &[Room], vent: &WallVent) -> Result<(), NetworkError> {
    admit_endpoints(rooms.len(), vent.left, vent.right)?;
    for endpoint in [vent.left, vent.right] {
        if let Endpoint::Room(i) = endpoint {
            let floor = rooms[i].floor_elevation.value();
            let ceiling = floor + rooms[i].gas.height().value();
            if vent.opening.bottom().value() < floor || vent.opening.top().value() > ceiling {
                return Err(NetworkError::Invalid("opening outside room height"));
            }
        }
    }
    Ok(())
}
fn column(
    rooms: &[Room],
    ambient: &Ambient,
    ambient_density: f64,
    pressures: &[f64],
    endpoint: Endpoint,
) -> Result<TwoLayerColumn, NetworkError> {
    Ok(match endpoint {
        Endpoint::Ambient => TwoLayerColumn::new(
            ambient.reference_height,
            ambient.pressure,
            ambient.reference_height,
            Density::new(ambient_density),
            Density::new(ambient_density),
        )?,
        Endpoint::Room(i) => {
            let room = &rooms[i];
            let lower = room.gas.layer(Layer::Lower).density_kg_m3();
            let upper = room.gas.layer(Layer::Upper).density_kg_m3();
            let lower = lower
                .or(upper)
                .ok_or(NetworkError::Invalid("empty compartment"))?;
            let upper = upper.unwrap_or(lower);
            let ratio = pressures[i] / room.gas.pressure().value();
            TwoLayerColumn::new(
                room.floor_elevation,
                Pressure::new(pressures[i]),
                Length::new(room.floor_elevation.value() + room.gas.interface_height().value()),
                Density::new(lower * ratio),
                Density::new(upper * ratio),
            )?
        }
    })
}
fn carried(
    rooms: &[Room],
    ambient: &Ambient,
    from: Endpoint,
    layer: Layer,
    rate: MassFlowRate,
) -> Result<LayerRates, NetworkError> {
    Ok(match from {
        Endpoint::Room(i) => rooms[i].gas.carried_rate(layer, rate)?,
        Endpoint::Ambient => {
            rooms[0]
                .gas
                .inflow(rate, &ambient.mass_fractions, ambient.temperature)?
        }
    })
}
fn transfer(
    rates: &mut [[LayerRates; 2]],
    external: &mut LayerRates,
    from: Endpoint,
    from_layer: Layer,
    to: Endpoint,
    to_layer: Layer,
    carried: LayerRates,
) -> Result<(), NetworkError> {
    match from {
        Endpoint::Room(i) => rates[i][from_layer.index()].add(&carried.opposite())?,
        Endpoint::Ambient => external.add(&carried)?,
    }
    match to {
        Endpoint::Room(i) => rates[i][to_layer.index()].add(&carried)?,
        Endpoint::Ambient => external.add(&carried.opposite())?,
    }
    Ok(())
}
