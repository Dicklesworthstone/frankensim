//! Actual deposition events and thermal substeps drive a growing clamped road.
//!
//! One thermal cell corresponds to one single-layer mechanical cell. Source
//! centreline coordinates describe the horizontal material REFERENCE domain;
//! the beam exposes the separately solved current geometry. New material is
//! assumed perfectly joined to the current tip, not forced onto a nozzle path.
//! This consumer uses the existing thermal observer and beam owners,
//! not another time integrator, heat-capacity law or constitutive model.

use std::collections::BTreeMap;
use fs_conduction::transient::deposition::{CellId, DepositedCell, DepositionError,
    DepositionObserver, DepositionSimulation, EnergyAudit, Event, MaterialId, StepConfig};
use fs_material::thermovisco::ThermoMaxwell;
use fs_solid::bridge::{self, Bridge};
use fs_solid::deposition::LayerGeometry;
use fs_solid::bridge::growth::CellBirth;

/// Explicit thermal/material correspondence for one future road cell.
#[derive(Debug, Clone, Copy)]
pub struct Binding<'a> {
    pub cell: CellId,
    pub thermal_material: MaterialId,
    pub mechanical_law: &'a ThermoMaxwell,
}

#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub mechanics: bridge::Config,
    pub gravity_m_s2: f64,
    /// Admission tolerances for the initial seed and incoming reference cells.
    pub coordinate_tolerance_m: f64,
    pub volume_tolerance_m3: f64,
    pub mass_tolerance_kg: f64,
    pub temperature_tolerance_k: f64,
    /// Seed plus all declared future bindings, not just currently active cells.
    pub max_bound_cells: usize,
    /// One shared allowance for initialized and evolved Gauss points per call.
    pub max_material_point_updates: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Thermal(DepositionError),
    Mechanical(bridge::Error),
    Invalid(&'static str),
    Budget,
}
impl From<DepositionError> for Error { fn from(e: DepositionError) -> Self { Self::Thermal(e) } }
impl From<bridge::Error> for Error { fn from(e: bridge::Error) -> Self { Self::Mechanical(e) } }
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "growing thermal road refused: {self:?}")
    }
}
impl std::error::Error for Error {}

/// Complete per-call accounting, not a coupled thermodynamic certificate.
#[derive(Debug, Clone)]
pub struct Report {
    pub thermal: EnergyAudit,
    pub born_cells: usize,
    pub material_point_updates: usize,
    /// Last actual equilibrium in this call; absent when no thermal step ran.
    pub last_mechanics: Option<bridge::Report>,
    /// The seed or latest birth has not yet received a positive-time solve here.
    pub birth_equilibrium_pending: bool,
}

#[derive(Clone)]
struct History<'a> {
    beam: Bridge<'a>,
    body: u64,
    layer: u64,
    geometry: LayerGeometry,
    bindings: BTreeMap<CellId, Binding<'a>>,
    config: Config,
    updates: usize,
    births: usize,
    last: Option<bridge::Report>,
    pending: bool,
}

impl History<'_> {
    fn reserve(&mut self, updates: usize) -> Result<(), Error> {
        self.updates = self.updates.checked_add(updates).ok_or(Error::Budget)?;
        if self.updates > self.config.max_material_point_updates { return Err(Error::Budget); }
        Ok(())
    }
}
impl DepositionObserver for History<'_> {
    type Error = Error;
    fn begin(&mut self) -> Result<(), Error> {
        self.updates = 0; self.births = 0; self.last = None;
        Ok(())
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Error> {
        let Event::Birth(cell) = event else { return Ok(()); };
        if cell.body_id != self.body {
            if self.bindings.contains_key(&cell.id) {
                return Err(Error::Invalid("bound birth belongs to a different physical body"));
            }
            return Ok(());
        }
        let binding = *self.bindings.get(&cell.id)
            .ok_or(Error::Invalid("road birth lacks an explicit mechanical binding"))?;
        if binding.thermal_material != cell.material_id || thermal.time_s() != self.beam.time_s() {
            return Err(Error::Invalid("birth material or physical clock mismatch"));
        }
        let length = cell.geometry.end_m[0]-cell.geometry.start_m[0];
        check_geometry(cell, self.beam.length_m(), length, self.geometry, self.config)?;
        self.reserve(2)?;
        let birth = CellBirth { id: cell.id.0, length_m: length, mass_kg: cell.mass_kg,
            temperature_k: thermal.temperature_k(cell.id)?, law: binding.mechanical_law };
        // The thermal owner polls immediately around each event. This bounded
        // birth/clone operation has no callback in the existing observer trait.
        let (beam, _) = self.beam.append_cantilever_controlled(&[birth],
            self.config.mechanics, &mut || false)?;
        self.beam = beam;
        self.bindings.remove(&cell.id);
        self.births += 1;
        self.pending = true;
        self.last = None;
        Ok(())
    }
    fn step(&mut self, start_s: f64, _: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Error> {
        checkpoint()?;
        if self.beam.time_s() != start_s {
            return Err(Error::Invalid("thermal substep and mechanical history are out of phase"));
        }
        self.reserve(self.beam.cells().len().checked_mul(2).ok_or(Error::Budget)?)?;
        let mut temperatures = BTreeMap::new();
        for cell in self.beam.cells() {
            checkpoint()?;
            temperatures.insert((cell.id(), self.layer), thermal.temperature_k(CellId(cell.id()))?);
        }
        let mut fault = None;
        let result = self.beam.advance_controlled(thermal.time_s(), &temperatures,
            self.config.gravity_m_s2, &vec![0.0;self.beam.cells().len()], self.config.mechanics,
            &mut || match checkpoint() {
                Ok(()) => false,
                Err(error) => { fault = Some(error); true }
            });
        if let Some(error) = fault { return Err(error.into()); }
        let (beam, report) = result?;
        self.beam = beam;
        self.last = Some(report);
        self.pending = false;
        Ok(())
    }
}

/// Accepted heat and growing mechanics on one physical timeline. Clone is an
/// in-memory checkpoint of both, not durable serialization or a full history log.
#[derive(Clone)]
pub struct GrowingRoad<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    history: History<'m>,
}
impl<'t, 'm> GrowingRoad<'t, 'm> {
    /// Bind an already activated single-layer clamped seed and future births.
    /// Cell IDs, body, reference centreline/volume, masses, temperatures and
    /// clocks must correspond. Other thermal bodies may exist and exchange heat;
    /// their mechanics are NOT supplied by this owner. Every later birth of the
    /// selected body must have an explicit binding or the whole advance refuses.
    pub fn new(thermal: DepositionSimulation<'t>, beam: Bridge<'m>, body: u64,
        future: Vec<Binding<'m>>, config: Config,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>) -> Result<Self, Error> {
        checkpoint()?;
        for value in [config.coordinate_tolerance_m, config.volume_tolerance_m3,
            config.mass_tolerance_kg, config.temperature_tolerance_k] {
            if !value.is_finite() || value < 0.0 { return Err(Error::Invalid("finite nonnegative binding tolerances required")); }
        }
        if !config.gravity_m_s2.is_finite() || config.gravity_m_s2 < 0.0
            || config.max_material_point_updates == 0 || thermal.time_s() != beam.time_s() {
            return Err(Error::Invalid("invalid gravity, work allowance or physical clock"));
        }
        if beam.cells().len().checked_add(future.len()).ok_or(Error::Budget)? > config.max_bound_cells {
            return Err(Error::Budget);
        }
        // Validate the growth model through the physical owner, without aging.
        beam.append_cantilever_controlled(&[], config.mechanics, &mut || false)?;
        let (&layer, first) = beam.cells()[0].section().layers().iter().next().expect("one born layer");
        let geometry = first.geometry;
        let mut active = std::collections::BTreeSet::new();
        for cell in beam.cells() {
            checkpoint()?;
            let id = CellId(cell.id());
            let source = thermal.cells().get(&id).ok_or(Error::Invalid("mechanical seed cell is not thermally active"))?;
            if source.body_id != body || !near(source.mass_kg, cell.mass_kg(), config.mass_tolerance_kg) {
                return Err(Error::Invalid("seed mass or physical body differs"));
            }
            check_geometry(source, cell.start_m(), cell.length_m(), geometry, config)?;
            let temperature = thermal.temperature_k(id)?;
            for point in cell.section().layers()[&layer].material_points() {
                if !near(point.temperature_k(), temperature, config.temperature_tolerance_k) {
                    return Err(Error::Invalid("seed material and thermal temperatures differ"));
                }
            }
            active.insert(id);
        }
        for (id, cell) in thermal.cells() {
            checkpoint()?;
            if cell.body_id == body && !active.contains(id) {
                return Err(Error::Invalid("thermal seed contains unrepresented road material"));
            }
        }
        let mut bindings = BTreeMap::new();
        for binding in future {
            checkpoint()?;
            if thermal.cells().contains_key(&binding.cell) || bindings.insert(binding.cell, binding).is_some() {
                return Err(Error::Invalid("future binding is active or duplicated"));
            }
        }
        checkpoint()?;
        Ok(Self { thermal, history: History { beam, body, layer, geometry, bindings, config,
            updates: 0, births: 0, last: None, pending: true } })
    }
    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn beam(&self) -> &Bridge<'m> { &self.history.beam }

    /// Every actual CFL-limited thermal step supplies the mechanical endpoint
    /// temperatures; every ordered birth supplies its real incoming mass and
    /// enthalpy-derived temperature. No final-temperature reconstruction or
    /// independent mechanics resampling is used. Births never age in advance.
    ///
    /// Thermal, material, geometry, work-budget and cancellation failures drop
    /// BOTH proposed states. Exact-time births may leave pending equilibrium;
    /// a positive-time step is needed before calling the enlarged road solved.
    /// The thermal ledger and last mechanical ledger remain separate: no
    /// viscous heating, fracture, moving contact or coupled passivity is claimed.
    pub fn advance_to(&self, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>) -> Result<(Self, Report), Error> {
        let (thermal, history, audit) = self.thermal.advance_observed(&self.history, target_s, config, checkpoint)?;
        let report = Report { thermal: audit, born_cells: history.births,
            material_point_updates: history.updates, last_mechanics: history.last.clone(),
            birth_equilibrium_pending: history.pending };
        Ok((Self { thermal, history }, report))
    }
}

fn near(a: f64, b: f64, tolerance: f64) -> bool {
    a.is_finite() && b.is_finite() && (a-b).abs() <= tolerance
}
fn check_geometry(cell: &DepositedCell, start: f64, length: f64,
    section: LayerGeometry, config: Config) -> Result<(), Error> {
    if !length.is_finite() || length <= 0.0 {
        return Err(Error::Invalid("positive finite reference cell length required"));
    }
    let a = [start, 0.0, section.centroid_z_m];
    let b = [start+length, 0.0, section.centroid_z_m];
    if !cell.geometry.start_m.into_iter().zip(a).chain(cell.geometry.end_m.into_iter().zip(b))
        .all(|(x,y)| near(x,y,config.coordinate_tolerance_m))
        || !near(cell.geometry.volume_m3, length*section.area_m2, config.volume_tolerance_m3) {
        return Err(Error::Invalid("thermal reference cell does not continue the mechanical cross section"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
