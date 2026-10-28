//! Native nodal enthalpy drives the existing continuing thermal-valve motion.
//!
//! The spatial heat owner supplies local material states, the common valve
//! reducer supplies equilibrium and tangent modes, and fs-modal transfers the
//! same physical displacement and momentum. Thermal discretization does not
//! select a different mechanics, observation or state-transfer implementation.

use std::convert::Infallible;

use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_modal::transfer::{ModalTransfer, ModalTransferControls};
use fs_solid::reduce::thermal::ThermalPressureLoadedValve;
use fs_solid::reduce::thermal::motion::{
    ThermalValveModalMotion, ThermalValveMotionError, ThermalValveMotionFrame,
};

use super::{
    EnthalpyStepConfig, EnthalpyStepSolution, ScalarField, SpatialSolidEnthalpy,
    SpatialSolidEnthalpyError, SpatialSolidState, ThermalBoundary, ThermalInterfaces,
};
use crate::solid_enthalpy::valve::{HeatedValveError, HeatedValveReduction};

/// A failed proposal changes neither nodal heat nor accepted motion.
#[derive(Debug)]
pub enum SpatialValveMotionError<E = Infallible> {
    /// Original native heat, source or solid-state refusal.
    Spatial(SpatialSolidEnthalpyError),
    /// Original equilibrium, pressure projection or modal reduction refusal.
    Valve(HeatedValveError),
    /// Original modal transition or physical state-transfer refusal.
    Motion(ThermalValveMotionError),
    /// A new source-basis observation could not be formed.
    Observation(E),
    /// Cancellation before joint publication.
    Cancelled,
}

impl<E> From<SpatialSolidEnthalpyError> for SpatialValveMotionError<E> {
    fn from(error: SpatialSolidEnthalpyError) -> Self {
        Self::Spatial(error)
    }
}
impl<E> From<HeatedValveError> for SpatialValveMotionError<E> {
    fn from(error: HeatedValveError) -> Self {
        Self::Valve(error)
    }
}
impl<E> From<ThermalValveMotionError> for SpatialValveMotionError<E> {
    fn from(error: ThermalValveMotionError) -> Self {
        Self::Motion(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for SpatialValveMotionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Spatial(error) => write!(f, "spatial valve motion: {error}"),
            Self::Valve(error) => write!(f, "spatial valve motion: {error}"),
            Self::Motion(error) => write!(f, "spatial valve motion: {error}"),
            Self::Observation(error) => write!(f, "spatial valve observation: {error}"),
            Self::Cancelled => f.write_str("spatial valve motion cancelled"),
        }
    }
}
impl<E: core::fmt::Debug + core::fmt::Display> std::error::Error for SpatialValveMotionError<E> {}

fn poll<E>(cx: &Cx<'_>) -> Result<(), SpatialValveMotionError<E>> {
    cx.checkpoint()
        .map_err(|_| SpatialValveMotionError::Cancelled)
}

impl SpatialSolidState<'_> {
    /// Lower the exact nodal-heat-derived material field through the SAME
    /// equilibrium, initial-stress, aperture and modal reducer as finite cells.
    /// Reference density, local thermal strain and frame remain source-owned.
    pub fn reduce_pressure_valve(
        &self,
        reduction: &HeatedValveReduction<'_>,
        cx: &Cx<'_>,
    ) -> Result<ThermalPressureLoadedValve, HeatedValveError> {
        reduction.reduce_state(
            &self.problem(reduction.fixed_dofs),
            self.thermal_strains(),
            self.mesh_frame,
            self.maximum_absolute_mandel_strain,
            cx,
        )
    }
}

/// One accepted nodal thermal field and its continuing native modal trajectory.
///
/// Reference geometry, material sources, mass, supports and pressure footprint
/// stay fixed. Each interval freezes the accepted mechanical law, advances its
/// exact held-pressure motion, advances native backward-Euler heat, then
/// transfers physical motion to the new equilibrium and tangent basis. The
/// complete heat field is retained; it is never replaced with a mean-temperature
/// thermal store. Only the original solid quadrature uses enthalpy centroids.
///
/// This first-order partition has explicit parameter work and projection loss.
/// Neither is thermalized. Reciprocal thermoelastic heat, radiation reaction,
/// moving geometry and phase changes require their own admitted owners. The
/// observation callback must be a pure proposal, with no external publication.
pub struct SpatialHeatedValveMotion<'s, 'a, 'm, 'c, 'r> {
    source: &'s SpatialSolidEnthalpy<'a, 'm, 'c>,
    specific_enthalpy_j_kg: Vec<f64>,
    reduction: HeatedValveReduction<'r>,
    motion: ThermalValveModalMotion,
}

/// Original thermal and mechanical accounts for one jointly accepted interval.
#[derive(Clone, Debug)]
pub struct SpatialValveMotionStep {
    /// Elapsed time from the supplied initial heat/motion state [s].
    pub time_s: f64,
    /// Held-pressure motion using the accepted start-of-step coefficients.
    pub pressure_step: ThermalValveMotionFrame,
    /// Complete native nodal heat result, including its original energy audit.
    pub thermal: EnthalpyStepSolution,
    /// Existing physical shape-transfer work and discarded-state accounting.
    pub transfer: ModalTransfer,
    /// Diagnostic pressure on the NEW source basis and transferred state [Pa].
    pub endpoint_observer_pressure_pa: f64,
}

impl<'s, 'a, 'm: 'a, 'c: 'a, 'r> SpatialHeatedValveMotion<'s, 'a, 'm, 'c, 'r> {
    /// Bind the initial vibration to the reduction of this exact nodal heat
    /// field. One explicit pressure transfer is required per actual mode; no
    /// loss or radiation model is inferred from a material name.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new<E>(
        source: &'s SpatialSolidEnthalpy<'a, 'm, 'c>,
        specific_enthalpy_j_kg: Vec<f64>,
        reduction: HeatedValveReduction<'r>,
        sample_rate_hz: u32,
        initial: &[ModalAcousticState],
        observe: impl FnOnce(&ThermalPressureLoadedValve) -> Result<Vec<C64>, E>,
        budget: ModalAcousticTimeBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, SpatialValveMotionError<E>> {
        poll(cx)?;
        let state = source.resolve(&specific_enthalpy_j_kg, cx)?;
        let valve = state.reduce_pressure_valve(&reduction, cx)?;
        let transfers = observe(&valve).map_err(SpatialValveMotionError::Observation)?;
        let motion =
            ThermalValveModalMotion::new(valve, sample_rate_hz, initial, &transfers, budget, cx)?;
        poll(cx)?;
        Ok(Self {
            source,
            specific_enthalpy_j_kg,
            reduction,
            motion,
        })
    }

    /// Accepted nodal heat in the original conduction mesh's vertex order.
    pub fn specific_enthalpy_j_kg(&self) -> &[f64] {
        &self.specific_enthalpy_j_kg
    }

    /// Accepted pressure port, thermal tangent, modal state and sample clock.
    pub const fn motion(&self) -> &ThermalValveModalMotion {
        &self.motion
    }

    /// Resolve source receipts from the same accepted heat field as mechanics.
    pub fn state(&self, cx: &Cx<'_>) -> Result<SpatialSolidState<'_>, SpatialSolidEnthalpyError> {
        self.source.resolve(&self.specific_enthalpy_j_kg, cx)
    }

    /// Advance heat, material operators, vibration and observation atomically.
    ///
    /// Boundary, contact and volumetric generation stay with native conduction;
    /// conductivity stays with the original material binding. A source-domain,
    /// phase, equilibrium, observation, transfer, budget or cancellation refusal
    /// publishes no heat, motion, pressure history or sample increment. Refine
    /// time, spatial mesh and retained modes to assess the declared split.
    #[allow(clippy::too_many_arguments)]
    pub fn step<E>(
        &mut self,
        pressure_pa: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        source: &ScalarField,
        thermal_config: EnthalpyStepConfig,
        controls: ModalTransferControls,
        observe: impl FnOnce(&ThermalPressureLoadedValve) -> Result<Vec<C64>, E>,
        cx: &Cx<'_>,
    ) -> Result<SpatialValveMotionStep, SpatialValveMotionError<E>> {
        poll(cx)?;
        let mut motion = self.motion.clone();
        let pressure_step = motion.step_pressure(pressure_pa, cx)?;
        let thermal = self.source.advance_thermal(
            &self.specific_enthalpy_j_kg,
            motion.model().sample_period_s(),
            boundary,
            interfaces,
            source,
            thermal_config,
            cx,
        )?;
        let state = self.source.resolve(&thermal.specific_enthalpy_j_kg, cx)?;
        let valve = state.reduce_pressure_valve(&self.reduction, cx)?;
        let transfers = observe(&valve).map_err(SpatialValveMotionError::Observation)?;
        let (motion, transfer) = motion.rebind_candidate(valve, &transfers, controls, cx)?;
        let endpoint_observer_pressure_pa = motion
            .model()
            .observer_pressure_with_transfers(&transfers)
            .map_err(ThermalValveMotionError::from)?;
        // Allocate the retained snapshot before the final poll. No fallible
        // operation separates publication of the two accepted physical owners.
        let enthalpy = thermal.specific_enthalpy_j_kg.clone();
        poll(cx)?;
        self.specific_enthalpy_j_kg = enthalpy;
        self.motion = motion;
        Ok(SpatialValveMotionStep {
            time_s: pressure_step.time_s,
            pressure_step,
            thermal,
            transfer,
            endpoint_observer_pressure_pa,
        })
    }
}
