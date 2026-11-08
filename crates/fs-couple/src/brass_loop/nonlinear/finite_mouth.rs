//! Finite upstream gas storage on the bore's existing two-stage gas clock.
//!
//! A rigid, acoustically compact Plenum receives an optional source-owned
//! nozzle feed and supplies the lip face. Its mass and internal energy use
//! exactly the same signed face exchanges as the bore. Incoming jet kinetic
//! energy thermalizes in the compact volume; its fixed support does no work.
//!
//! The opt-in noise rows (`jet_noise`, `constriction_noise`) are SERIES pressure sources between the stored gas and the lip
//! face. A trial carries their summed pressure `series_pa`; both stages of the step then offer the lip the stored gas at its
//! own temperature and the pressure raised by `series_pa` (an isothermal modulation: the specific enthalpy of the gas that
//! passes is the stored gas's own, so the storage is debited exactly what it delivers), while the feed nozzle and the stored
//! inventory never see the offset. The mechanical work such a source would do, `series_pa x` the volume flow into the bore, is
//! reported as `series_work_j`; the isothermal modulation rejects it as heat, so it is neither drawn from the storage nor added
//! to the bore. With `series_pa == 0.0` every expression is the one used before the offset existed, bit for bit.

use fs_gas::dynamics::network::Plenum;
use fs_gas::dynamics::{Conserved, Error, FaceFlux, FrozenGas, PortStage, Primitive, StepReport};

/// Finite, rigid, adiabatic mouth storage and its upstream feed opening.
///
/// The blowing-pressure control sets the stationary feed reservoir pressure,
/// relative to the room. Its temperature and composition remain the room's.
/// A zero feed area gives a sealed, initially charged supply. The source-owned
/// converging nozzle has unit discharge coefficient; no neck inertia, fitted
/// resistance, wall heat, moving volume or anatomical lung model is implied.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteMouthSpec {
    /// Acoustically compact storage volume [m^3].
    pub volume_m3: f64,
    /// Effective feed opening area [m^2], including any authored discharge factor.
    pub feed_area_m2: f64,
    /// Initial absolute storage pressure [Pa].
    pub initial_pressure_pa: f64,
    /// Initial storage temperature [K], within the original caloric source.
    pub initial_temperature_k: f64,
}

/// Accepted finite-mouth storage, decoded by the original frozen gas source.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteMouthState {
    /// Fixed storage volume [m^3].
    pub volume_m3: f64,
    /// Remaining mass [kg].
    pub mass_kg: f64,
    /// Internal energy [J] in the bore's original sensible-energy reference.
    /// Physical energy adds the source reference energy times the stored mass.
    pub energy_j: f64,
    /// Absolute static pressure [Pa].
    pub pressure_pa: f64,
    /// Static temperature [K].
    pub temperature_k: f64,
}

/// Accepted gas-step exchange; positive feed enters the mouth, positive bore
/// exchange leaves it. Negative values are physical reverse flow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteMouthReport {
    /// Common accepted step [s].
    pub elapsed_s: f64,
    /// Storage before this step.
    pub before: FiniteMouthState,
    /// Storage after this step.
    pub after: FiniteMouthState,
    /// Net mass delivered by the held feed reservoir [kg].
    pub feed_mass_kg: f64,
    /// Net feed stagnation enthalpy delivered [J], in the common energy reference.
    pub feed_energy_j: f64,
    /// Net mass delivered to the bore, from its accepted native boundary report [kg].
    pub bore_mass_kg: f64,
    /// Net total enthalpy delivered to the bore, from that same report [J].
    pub bore_energy_j: f64,
    /// Storage change minus feed plus bore exchange [kg].
    pub mass_residual_kg: f64,
    /// Storage change minus feed plus bore exchange [J].
    pub energy_residual_j: f64,
    /// Largest compact-volume face-speed CFL over the two accepted stages.
    pub maximum_cfl: f64,
    /// Mechanical work [J] a series noise source of pressure `series_pa` would do on the volume that crossed the lip face
    /// (the trapezoid over the two stages of `series_pa x inlet_area x face mass flux / offered density`); exactly zero without
    /// a series offset. It is rejected as heat by the isothermal modulation and is in no mass or energy account of the storage.
    pub series_work_j: f64,
}

pub(super) struct FiniteMouth {
    plenum: Plenum,
    feed: Primitive,
    feed_area_m2: f64,
    report: Option<FiniteMouthReport>,
}

impl FiniteMouth {
    pub(super) fn new(
        gas: &FrozenGas,
        spec: FiniteMouthSpec,
        feed: Primitive,
    ) -> Result<(Self, Primitive), Error> {
        if !(spec.volume_m3.is_finite() && spec.volume_m3 > 0.0)
            || !(spec.feed_area_m2.is_finite() && spec.feed_area_m2 >= 0.0)
            || !(spec.initial_pressure_pa.is_finite() && spec.initial_pressure_pa > 0.0)
            || !(spec.initial_temperature_k.is_finite() && spec.initial_temperature_k > 0.0)
        {
            return Err(Error::Invalid("finite mouth volume, feed area, pressure or temperature"));
        }
        let density = spec.initial_pressure_pa / (gas.gas_constant() * spec.initial_temperature_k);
        let state = gas.conserved(density, 0.0, spec.initial_pressure_pa)?;
        let primitive = gas.primitive(state, || false)?;
        let plenum = Plenum::new(spec.volume_m3, state)?;
        Ok((Self { plenum, feed, feed_area_m2: spec.feed_area_m2, report: None }, primitive))
    }

    pub(super) fn set_feed(&mut self, feed: Primitive) {
        self.feed = feed;
    }

    pub(super) fn state(&self, primitive: Primitive) -> FiniteMouthState {
        inventory(&self.plenum, primitive)
    }

    pub(super) const fn report(&self) -> Option<&FiniteMouthReport> {
        self.report.as_ref()
    }

    pub(super) fn trial(
        &self,
        primitive: Primitive,
        dt: f64,
        inlet_area_m2: f64,
        series_pa: f64,
    ) -> FiniteMouthTrial {
        FiniteMouthTrial {
            series_pa,
            initial: self.plenum.clone(),
            initial_primitive: primitive,
            feed: self.feed,
            feed_area_m2: self.feed_area_m2,
            inlet_area_m2,
            dt,
            predictor: None,
            first: None,
            candidate: None,
        }
    }

    pub(super) fn commit(
        &mut self,
        trial: FiniteMouthTrial,
        bore: &StepReport,
    ) -> Result<Primitive, Error> {
        let next = trial.candidate.ok_or(Error::Invalid("finite mouth has no accepted corrector"))?;
        let before = inventory(&trial.initial, trial.initial_primitive);
        let after = inventory(&next.plenum, next.primitive);
        let bore_mass_kg = bore.left_boundary_change.mass;
        let bore_energy_j = bore.left_boundary_change.energy;
        let report = FiniteMouthReport {
            elapsed_s: bore.elapsed.value(),
            before,
            after,
            feed_mass_kg: next.feed_mass_kg,
            feed_energy_j: next.feed_energy_j,
            bore_mass_kg,
            bore_energy_j,
            mass_residual_kg: after.mass_kg - before.mass_kg - next.feed_mass_kg + bore_mass_kg,
            energy_residual_j: after.energy_j - before.energy_j - next.feed_energy_j + bore_energy_j,
            maximum_cfl: next.maximum_cfl,
            series_work_j: next.series_work_j,
        };
        self.plenum = next.plenum;
        self.report = Some(report);
        Ok(next.primitive)
    }
}

fn inventory(plenum: &Plenum, primitive: Primitive) -> FiniteMouthState {
    FiniteMouthState {
        volume_m3: plenum.volume_m3(),
        mass_kg: plenum.volume_m3() * plenum.state().density,
        energy_j: plenum.volume_m3() * plenum.state().energy,
        pressure_pa: primitive.pressure,
        temperature_k: primitive.temperature,
    }
}

#[derive(Clone, Copy)]
struct Exchange {
    feed_mass_rate: f64,
    feed_energy_rate: f64,
    net_mass_rate: f64,
    net_energy_rate: f64,
    series_work_rate: f64,
    cfl: f64,
}

struct Candidate {
    plenum: Plenum,
    primitive: Primitive,
    feed_mass_kg: f64,
    feed_energy_j: f64,
    maximum_cfl: f64,
    series_work_j: f64,
}

/// Stack-only candidate storage. Every repeated stage zero starts again from
/// the accepted inventory, including the bore's MUSCL-to-first-order retry.
pub(super) struct FiniteMouthTrial {
    /// The summed series noise pressure of this gas step [Pa]; zero when no noise source speaks.
    series_pa: f64,
    initial: Plenum,
    initial_primitive: Primitive,
    feed: Primitive,
    feed_area_m2: f64,
    inlet_area_m2: f64,
    dt: f64,
    predictor: Option<Primitive>,
    first: Option<Exchange>,
    candidate: Option<Candidate>,
}

impl FiniteMouthTrial {
    pub(super) fn flux(
        &mut self,
        gas: &FrozenGas,
        stage: &PortStage,
        outlet: impl FnOnce(Primitive) -> Result<FaceFlux, Error>,
    ) -> Result<FaceFlux, Error> {
        if stage.dt != self.dt || !(self.dt.is_finite() && self.dt > 0.0) || stage.stage > 1 {
            return Err(Error::Invalid("finite mouth requires the bore's fixed two-stage clock"));
        }
        let mouth = if stage.stage == 0 {
            self.predictor = None;
            self.first = None;
            self.candidate = None;
            self.initial_primitive
        } else {
            self.predictor.ok_or(Error::Invalid("finite mouth corrector needs its predictor"))?
        };
        let (feed_mass, feed_energy, feed_signal) =
            feed_flux(gas, self.feed, mouth, self.feed_area_m2)?;
        // The feed nozzle and the inventory use the stored gas as it is; only the lip face is offered the series offset.
        let offered = if self.series_pa == 0.0 { mouth } else { offer(gas, mouth, self.series_pa)? };
        let face = outlet(offered)?;
        let series_work_rate = if self.series_pa == 0.0 {
            0.0
        } else {
            self.series_pa * self.inlet_area_m2 * face.mass / offered.density
        };
        let cfl = self.dt
            * (self.feed_area_m2 * feed_signal + self.inlet_area_m2 * face.signal_speed)
            / self.initial.volume_m3();
        if !cfl.is_finite() || cfl > 0.45 {
            return Err(Error::CflExceeded { observed: cfl, maximum: 0.45 });
        }
        let exchange = Exchange {
            feed_mass_rate: feed_mass,
            feed_energy_rate: feed_energy,
            net_mass_rate: feed_mass - self.inlet_area_m2 * face.mass,
            net_energy_rate: feed_energy - self.inlet_area_m2 * face.energy,
            series_work_rate,
            cfl,
        };
        let initial = self.initial.state();
        let volume = self.initial.volume_m3();
        let (mass_change, energy_change) = if stage.stage == 0 {
            (self.dt * exchange.net_mass_rate, self.dt * exchange.net_energy_rate)
        } else {
            let first = self.first.ok_or(Error::Invalid("finite mouth lost its first-stage exchange"))?;
            (
                (0.5 * self.dt) * (first.net_mass_rate + exchange.net_mass_rate),
                (0.5 * self.dt) * (first.net_energy_rate + exchange.net_energy_rate),
            )
        };
        let state = Conserved {
            density: initial.density + mass_change / volume,
            momentum: 0.0,
            energy: initial.energy + energy_change / volume,
        };
        let primitive = gas.primitive(state, || false)?;
        let plenum = Plenum::new(volume, state)?;
        if stage.stage == 0 {
            self.predictor = Some(primitive);
            self.first = Some(exchange);
        } else {
            let first = self.first.ok_or(Error::Invalid("finite mouth lost its first-stage exchange"))?;
            self.candidate = Some(Candidate {
                plenum,
                primitive,
                feed_mass_kg: (0.5 * self.dt) * (first.feed_mass_rate + exchange.feed_mass_rate),
                feed_energy_j: (0.5 * self.dt) * (first.feed_energy_rate + exchange.feed_energy_rate),
                maximum_cfl: first.cfl.max(exchange.cfl),
                series_work_j: (0.5 * self.dt) * (first.series_work_rate + exchange.series_work_rate),
            });
        }
        Ok(face)
    }
}

// The stored gas at rest, raised by a series pressure `series_pa` at its own temperature (the isothermal modulation the
// original reservoir applies at the room temperature). Its specific enthalpy is the stored gas's own.
fn offer(gas: &FrozenGas, mouth: Primitive, series_pa: f64) -> Result<Primitive, Error> {
    let pressure = mouth.pressure + series_pa;
    if !(series_pa.is_finite() && pressure.is_finite() && pressure > 0.0) {
        return Err(Error::Invalid("finite mouth series pressure offset leaves no positive pressure"));
    }
    let density = pressure / (gas.gas_constant() * mouth.temperature);
    gas.primitive(gas.conserved(density, 0.0, pressure)?, || false)
}

// Both compact owners are stationary. The donor supplies its own stagnation
// enthalpy; jet kinetic energy is already included and is not added twice.
fn feed_flux(
    gas: &FrozenGas,
    feed: Primitive,
    mouth: Primitive,
    area: f64,
) -> Result<(f64, f64, f64), Error> {
    if area == 0.0 {
        return Ok((0.0, 0.0, 0.0));
    }
    if feed.pressure == mouth.pressure {
        return Ok((0.0, 0.0, feed.sound_speed.max(mouth.sound_speed)));
    }
    let (donor, backpressure, sign) = if feed.pressure > mouth.pressure {
        (feed, mouth.pressure, 1.0)
    } else {
        (mouth, feed.pressure, -1.0)
    };
    let (jet, _) = gas.nozzle_exit(donor, backpressure, || false)?;
    let mass = sign * area * jet.density * jet.velocity;
    let energy = mass * gas.stagnation_enthalpy(donor.temperature)?;
    if !mass.is_finite() || !energy.is_finite() {
        return Err(Error::Invalid("finite mouth feed mass or enthalpy overflow"));
    }
    Ok((mass, energy, jet.velocity + jet.sound_speed))
}

#[cfg(test)]
mod tests;
