//! Retain the original thermal network while observing impact pressure.
//!
//! The existing coupled steppers own mechanical work deposition, thermal
//! subcycling, contact refresh and liquid-film viscosity. This adapter keeps
//! their inputs alive across audio callbacks and invokes that step exactly once
//! per mechanical sample. It does not advance the adiabatic source first or
//! deposit the same work again in an observer.

use super::super::material::MaterialContact;
use super::super::material::heating::liquid_film::LiquidFilmImpactSystem;
use super::super::material::heating::network::{NetworkImpactFrame, ReservoirCell};
use super::super::material::heating::{CaloricReservoir, HeatedImpactSystem};
use super::{ImpactError, ImpactSource};
use fs_contact::normal_patch::source::SourcePointContact;
use fs_exec::CancelGate;
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

/// An existing impact image with an atomic spatial thermal-network step.
///
/// Implementations retain the same mechanical coordinates, clock and budgets
/// exposed by [`ImpactSource`]. A refusal must preserve source, network and
/// observer together. External effects captured by a material resolver remain
/// the caller's responsibility, as in the original coupled step.
pub trait ThermalNetworkImpact: ImpactSource {
    /// Invoke the original network-coupled step, including its material gates.
    fn advance_with_network<'n, O, R>(
        &mut self,
        thermal: &mut DepositionSimulation<'n>,
        observer: &mut O,
        bindings: &[ReservoirCell],
        config: StepConfig,
        forces: &[f64],
        gate: &CancelGate,
        resolve: &mut R,
    ) -> Result<NetworkImpactFrame, ImpactError>
    where
        O: DepositionObserver,
        O::Error: core::fmt::Display,
        R: for<'m> FnMut(
            &MaterialContact,
            &[CaloricReservoir<'m>],
        ) -> Result<Option<SourcePointContact>, ImpactError>;
}

impl ThermalNetworkImpact for HeatedImpactSystem<'_> {
    fn advance_with_network<'n, O, R>(
        &mut self,
        thermal: &mut DepositionSimulation<'n>,
        observer: &mut O,
        bindings: &[ReservoirCell],
        config: StepConfig,
        forces: &[f64],
        gate: &CancelGate,
        resolve: &mut R,
    ) -> Result<NetworkImpactFrame, ImpactError>
    where
        O: DepositionObserver,
        O::Error: core::fmt::Display,
        R: for<'m> FnMut(
            &MaterialContact,
            &[CaloricReservoir<'m>],
        ) -> Result<Option<SourcePointContact>, ImpactError>,
    {
        self.step_with_thermal_network(thermal, observer, bindings, config, forces, gate, resolve)
    }
}

impl ThermalNetworkImpact for LiquidFilmImpactSystem<'_> {
    fn advance_with_network<'n, O, R>(
        &mut self,
        thermal: &mut DepositionSimulation<'n>,
        observer: &mut O,
        bindings: &[ReservoirCell],
        config: StepConfig,
        forces: &[f64],
        gate: &CancelGate,
        resolve: &mut R,
    ) -> Result<NetworkImpactFrame, ImpactError>
    where
        O: DepositionObserver,
        O::Error: core::fmt::Display,
        R: for<'m> FnMut(
            &MaterialContact,
            &[CaloricReservoir<'m>],
        ) -> Result<Option<SourcePointContact>, ImpactError>,
    {
        self.step_with_thermal_network(thermal, observer, bindings, config, forces, gate, resolve)
    }
}

/// One retained mechanical source, spatial heat network and thermal observer.
///
/// This is an [`ImpactSource`] for [`super::ImpactPressureRenderer`], so the
/// existing radiation, propagation, decimation, pressure mixing and WAV paths
/// all observe the same network-coupled trajectory. No physical clock, state,
/// material law, heat capacity or dissipation partition is reconstructed here.
///
/// Each advance delegates complete binding, clock, material-domain, energy and
/// cancellation admission to the original owner. Its accepted frame retains
/// both the network audit and the combined mechanical/thermal account. The
/// owner's thermal splitting and unloaded-only contact refresh limitations
/// still apply; audio observation does not improve their physical accuracy.
pub struct ThermalNetworkSource<'n, M, O, R> {
    source: M,
    thermal: DepositionSimulation<'n>,
    observer: O,
    bindings: Vec<ReservoirCell>,
    config: StepConfig,
    resolve: R,
}

impl<'n, M, O, R> ThermalNetworkSource<'n, M, O, R>
where
    M: ThermalNetworkImpact,
    O: DepositionObserver,
    O::Error: core::fmt::Display,
    R: for<'m> FnMut(
        &MaterialContact,
        &[CaloricReservoir<'m>],
    ) -> Result<Option<SourcePointContact>, ImpactError>,
{
    /// Retain caller-supplied coupled state without advancing either clock.
    ///
    /// The original step checks the complete reservoir/cell binding on every
    /// advance. The observer supplies only its original environmental/material
    /// inputs; mechanical dissipation is already routed by the coupled owner.
    /// Return `Ok(None)` from `resolve` to retain the original contact sources.
    pub fn new(
        source: M,
        thermal: DepositionSimulation<'n>,
        observer: O,
        bindings: Vec<ReservoirCell>,
        config: StepConfig,
        resolve: R,
    ) -> Self {
        Self {
            source,
            thermal,
            observer,
            bindings,
            config,
            resolve,
        }
    }

    /// Accepted original mechanical and finite caloric histories.
    #[must_use]
    pub const fn source(&self) -> &M {
        &self.source
    }

    /// Accepted spatial thermal cells, links and physical clock.
    #[must_use]
    pub const fn thermal(&self) -> &DepositionSimulation<'n> {
        &self.thermal
    }

    /// Observer state published by the last accepted coupled step.
    #[must_use]
    pub const fn observer(&self) -> &O {
        &self.observer
    }

    /// Original one-to-one reservoir/cell identities.
    #[must_use]
    pub fn bindings(&self) -> &[ReservoirCell] {
        &self.bindings
    }
}

impl<M, O, R> ImpactSource for ThermalNetworkSource<'_, M, O, R>
where
    M: ThermalNetworkImpact,
    O: DepositionObserver,
    O::Error: core::fmt::Display,
    R: for<'m> FnMut(
        &MaterialContact,
        &[CaloricReservoir<'m>],
    ) -> Result<Option<SourcePointContact>, ImpactError>,
{
    type Frame = NetworkImpactFrame;

    fn state(&self) -> &[f64] {
        self.source.state()
    }
    fn mode_count(&self) -> usize {
        self.source.mode_count()
    }
    fn samples_rendered(&self) -> u64 {
        self.source.samples_rendered()
    }
    fn sample_period_s(&self) -> f64 {
        self.source.sample_period_s()
    }
    fn remaining_steps(&self) -> u64 {
        self.source.remaining_steps()
    }
    fn maximum_generalized_force(&self) -> f64 {
        self.source.maximum_generalized_force()
    }

    fn advance(&mut self, forces: &[f64], gate: &CancelGate) -> Result<Self::Frame, ImpactError> {
        self.source.advance_with_network(
            &mut self.thermal,
            &mut self.observer,
            &self.bindings,
            self.config,
            forces,
            gate,
            &mut self.resolve,
        )
    }
}
