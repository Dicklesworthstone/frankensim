//! Atomic self-heating intervals and scheduled stress-free layer births.
//!
//! The thermal timeline remains the only event scheduler. This owner applies
//! its real endpoint events AFTER the original coupled solve, on private clones.
//! Births never age in Newton trials. A failure in the last birth discards the
//! whole interval, including already successful thermal and mechanical trials.
//!
//! Callers split intervals at event times and load discontinuities. This is a
//! growing, axially uniform bonded section, not a spatial bridge or impact law.

use core::convert::Infallible;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use fs_solid::deposition::{Control, LayerGeometry, LayeredStrip};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionObserver,
    DepositionSimulation, EnergyAudit, Event, StepConfig};
use super::{LayerThermalBinding, StripEnthalpyConfig, StripEnthalpyError,
    StripEnthalpyReport, ThermalCoefficients};

/// A future perfectly bonded layer, activated only by its thermal cell's birth.
#[derive(Clone, Copy, Debug)]
pub struct LayerBirth<'m> {
    pub binding: LayerThermalBinding<'m>,
    /// Reference section geometry, not a shape inferred from cell proximity.
    pub geometry: LayerGeometry,
}

/// Source/event refusal or the original coupled interval failure.
#[derive(Debug)]
pub enum GrowthError<E> {
    Events(StripEnthalpyError<Infallible>),
    Coupled(StripEnthalpyError<E>),
}
impl<E: core::fmt::Debug> core::fmt::Display for GrowthError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "growing self-heating strip: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for GrowthError<E> {}

/// The interval and births have separate accounts; neither is counted twice.
pub struct GrowthReport<S> {
    /// Original left-limit coupled interval, before endpoint events.
    pub interval: StripEnthalpyReport<S>,
    /// Actual endpoint event mass and enthalpy, including passive-body births.
    pub events: EnergyAudit,
    pub born_layers: usize,
    /// Newborn material has not undergone an additional equilibrium/impact solve.
    pub birth_equilibrium_pending: bool,
    /// Full storage change minus imposed work, other heat and incoming enthalpy,
    /// plus numerical loss, measured across BOTH the interval and its events.
    pub total_energy_residual_j: f64,
}

#[derive(Clone)]
struct LayerEvents<'m> {
    strip: LayeredStrip<'m>,
    bindings: Vec<LayerThermalBinding<'m>>,
    future: BTreeMap<CellId, LayerBirth<'m>>,
    body: u64,
    born: usize,
    pending: bool,
}

impl LayerEvents<'_> {
    fn admit(&self, thermal: &DepositionSimulation<'_>,
        mut check: impl FnMut() -> Result<(), DepositionError>)
        -> Result<(), StripEnthalpyError<Infallible>> {
        super::solve::admit_bindings::<Infallible>(thermal, &self.strip, &self.bindings, &mut check)?;
        let mut represented = BTreeSet::new();
        for binding in &self.bindings {
            check()?;
            if thermal.cells()[&binding.cell].body_id != self.body {
                return Err(StripEnthalpyError::Invalid("bound cell belongs to another body"));
            }
            represented.insert(binding.cell);
        }
        for (id, cell) in thermal.cells() {
            check()?;
            if cell.body_id == self.body && !represented.contains(id) {
                return Err(StripEnthalpyError::Invalid("selected thermal body has unrepresented mechanical material"));
            }
        }
        check()?;
        Ok(())
    }
}

impl DepositionObserver for LayerEvents<'_> {
    type Error = StripEnthalpyError<Infallible>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.born = 0;
        Ok(())
    }

    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        let Event::Birth(cell) = event else { return Ok(()); };
        if cell.body_id != self.body {
            if self.future.contains_key(&cell.id) {
                return Err(StripEnthalpyError::Invalid("registered birth belongs to another body"));
            }
            return Ok(());
        }
        if self.strip.time_s() != thermal.time_s() {
            return Err(StripEnthalpyError::Invalid("birth and mechanical history clocks differ"));
        }
        let birth = self.future.remove(&cell.id)
            .ok_or(StripEnthalpyError::Invalid("selected-body birth has no material binding"))?;
        // Existing deposit initializes both points on the CURRENT strain/curvature.
        // Mass, caloric object, source context and physical limits are checked by
        // the same admission function as the nonlinear driver, after the batch.
        self.strip = self.strip.deposit(birth.binding.layer, birth.geometry,
            birth.binding.material.law(), thermal.temperature_k(cell.id)?)?;
        self.bindings.push(birth.binding);
        self.born += 1;
        self.pending = true;
        Ok(())
    }

    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>, _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error> {
        Err(StripEnthalpyError::Invalid("birth synchronization must not advance thermal time"))
    }
}

/// One accepted thermal/mechanical history and bounded future material bindings.
/// Clone is an in-memory checkpoint; laws are borrowed immutable source objects.
#[derive(Clone)]
pub struct GrowingStrip<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    layers: LayerEvents<'m>,
}

impl<'t, 'm> GrowingStrip<'t, 'm> {
    /// Bind an existing common-time state and process events due exactly now.
    /// Empty mechanical/thermal seeds are allowed when current-time births create
    /// a nonempty selected body. Initial bindings plus future births must fit the
    /// explicit registry limit; the original section/event budgets still apply.
    /// Returns the actual initialization event audit, including incoming enthalpy.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(thermal: DepositionSimulation<'t>, strip: LayeredStrip<'m>,
        body: u64, bindings: Vec<LayerThermalBinding<'m>>, future: Vec<LayerBirth<'m>>,
        maximum_registered_layers: usize, event_energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool)
        -> Result<(Self, EnergyAudit), StripEnthalpyError<Infallible>> {
        let mut check = || if cancelled() { Err(DepositionError::Cancelled) } else { Ok(()) };
        check()?;
        let count = bindings.len().checked_add(future.len())
            .ok_or(StripEnthalpyError::Invalid("layer registry size overflow"))?;
        if maximum_registered_layers == 0 || count > maximum_registered_layers {
            return Err(StripEnthalpyError::Invalid("layer registry exceeds its declared bound"));
        }
        let mut ids = BTreeSet::new();
        for binding in &bindings {
            check()?;
            if !ids.insert(binding.layer) {
                return Err(StripEnthalpyError::Invalid("duplicate active layer binding"));
            }
        }
        let mut registered = BTreeMap::new();
        for birth in future {
            check()?;
            if !ids.insert(birth.binding.layer) || thermal.cells().contains_key(&birth.binding.cell)
                || registered.insert(birth.binding.cell, birth).is_some() {
                return Err(StripEnthalpyError::Invalid("future cell/layer is active or duplicated"));
            }
        }
        let layers = LayerEvents { strip, bindings, future: registered, body, born: 0, pending: true };
        layers.admit(&thermal, &mut check)?;
        let (thermal, layers, events) = synchronize(&thermal, &layers, event_energy_tolerance_j, &mut check)?;
        if layers.strip.layers().is_empty() {
            return Err(StripEnthalpyError::Invalid("initialization requires a born mechanical seed"));
        }
        check()?;
        Ok((Self { thermal, layers }, events))
    }

    pub const fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn strip(&self) -> &LayeredStrip<'m> { &self.layers.strip }
    pub fn bindings(&self) -> &[LayerThermalBinding<'m>] { &self.layers.bindings }
    pub fn remaining_births(&self) -> usize { self.layers.future.len() }
    pub fn birth_equilibrium_pending(&self) -> bool { self.layers.pending }

    /// Solve one interval, then apply all actual target-time events atomically.
    /// Interior events are refused by the original thermal owner: split at their
    /// times rather than skipping them. No automatic timestep/error estimate is
    /// claimed. The returned state includes the right-limit endpoint births.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<S, E>(&self, target_s: f64, control: Control,
        config: StripEnthalpyConfig, initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<ThermalCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool)
        -> Result<(Self, GrowthReport<S>), GrowthError<E>> {
        let cancel = RefCell::new(cancelled);
        let mut check = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
        let candidate = super::step_candidate(&self.thermal, &self.layers.strip,
            &self.layers.bindings, target_s, control, config, initial_energy_increments_j,
            resolve, || (*cancel.borrow_mut())()).map_err(GrowthError::Coupled)?;
        let (thermal, strip, interval) = candidate.into_parts();
        let mut layers = self.layers.clone();
        layers.strip = strip;
        layers.pending = false;
        let (thermal, layers, events) = synchronize(&thermal, &layers,
            config.thermal_energy_tolerance_j, &mut check).map_err(GrowthError::Events)?;
        // Native births retain every existing Gauss-point history and add
        // stress-free points (Xi=0), so their mechanical entropy change is zero.
        // The actual full-interval spring storage is free energy; the interval's
        // original Delta(T*s) completes internal energy through the birth batch.
        let residual = (layers.strip.stored_energy_j() - self.layers.strip.stored_energy_j())
            + interval.mechanical_entropy_energy_change_j
            + (thermal.energy_j() - self.thermal.energy_j())
            + interval.mechanical.energy_j.numerical_dissipation
            - interval.external_mechanical_work_j - interval.non_mechanical_heat_j
            - events.incoming_energy_j - events.external_heat_j;
        if !residual.is_finite() || residual.abs() > config.total_energy_tolerance_j {
            return Err(GrowthError::Events(StripEnthalpyError::Invalid(
                "combined interval and birth energy account does not close")));
        }
        check().map_err(|e| GrowthError::Events(e.into()))?;
        let report = GrowthReport { interval, events, born_layers: layers.born,
            birth_equilibrium_pending: layers.pending, total_energy_residual_j: residual };
        Ok((Self { thermal, layers }, report))
    }
}

fn synchronize<'t, 'm>(thermal: &DepositionSimulation<'t>, layers: &LayerEvents<'m>,
    energy_tolerance_j: f64, mut check: impl FnMut() -> Result<(), DepositionError>)
    -> Result<(DepositionSimulation<'t>, LayerEvents<'m>, EnergyAudit), StripEnthalpyError<Infallible>> {
    // Same time and zero allowed substeps: max_dt is never used to age physics.
    let config = StepConfig { max_dt_s: 1.0, energy_tolerance_j };
    let (thermal, layers, audit) = thermal.advance_observed_limited(layers, thermal.time_s(),
        config, 0, &mut check)?;
    layers.admit(&thermal, &mut check)?;
    check()?;
    Ok((thermal, layers, audit))
}
