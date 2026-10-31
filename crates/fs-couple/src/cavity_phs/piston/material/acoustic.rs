//! Original material inertia retained through reciprocal acoustic motion and heat.
//!
//! These entry points bind the complete native candidates back to their original
//! source-resolved rigid body. No gas, acoustic, thermal or mechanical equation is
//! repeated here, and neither gas snapshots nor port clocks are reconstructed.

use super::{MaterialDoubleActingGasPiston, MaterialDoubleThermalPistonCandidate};
use crate::cavity_phs::piston::double_acting::{
    DoublePistonInput, DoublePistonSolve, DoublePistonStep,
    acoustic::{AcousticMotionError, AcousticMotionFace, AcousticMotionSolve, AcousticMotionStep},
};
use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use crate::cavity_phs::piston::{ClosedGasParcel, PistonError, PistonGas};
use crate::dg_port::DgPortLoad;
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;

/// Source-bound shared mechanism and every original acoustic endpoint together.
#[derive(Clone, Debug)]
pub struct MaterialAcousticPistonCandidate<G: PistonGas = ClosedGasParcel> {
    piston: MaterialDoubleActingGasPiston<G>,
    loads: Vec<DgPortLoad>,
    step: DoublePistonStep<G::State>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}
impl<G: PistonGas> MaterialAcousticPistonCandidate<G> {
    /// Complete mechanism with its original density, geometry and inertia records.
    pub const fn piston(&self) -> &MaterialDoubleActingGasPiston<G> {
        &self.piston
    }
    /// Tentative original DG and terminal-memory owners in authored face order.
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }
    /// Original gas/momentum report retaining both actual material midpoints.
    pub const fn piston_step(&self) -> &DoublePistonStep<G::State> {
        &self.step
    }
    /// Original reciprocal acoustic work, storage, losses and terminal exchange.
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }
    /// Original enclosing gas, wall and acoustic energy defect [J].
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }
    /// Publish the source-bound mechanism and all acoustic histories together.
    pub fn into_parts(
        self,
    ) -> (
        MaterialDoubleActingGasPiston<G>,
        Vec<DgPortLoad>,
        DoublePistonStep<G::State>,
        AcousticMotionStep,
        f64,
    ) {
        (
            self.piston,
            self.loads,
            self.step,
            self.motion,
            self.total_energy_residual_j,
        )
    }
}

/// Original material-bound finite-body candidate and all acoustic endpoints.
#[derive(Debug)]
pub struct MaterialAcousticThermalPistonCandidate<S, G: PistonGas = ClosedGasParcel> {
    native: MaterialDoubleThermalPistonCandidate<S, G>,
    loads: Vec<DgPortLoad>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}
impl<S, G: PistonGas> MaterialAcousticThermalPistonCandidate<S, G> {
    /// Both gases, one source-bound moving body and the finite thermal store.
    pub const fn native(&self) -> &MaterialDoubleThermalPistonCandidate<S, G> {
        &self.native
    }
    /// Original acoustic and terminal-memory candidates in authored face order.
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }
    /// Original complete reciprocal exchange and acoustic energy account.
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }
    /// Original gas, single inertia, thermal body and acoustic defect [J].
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }
    /// Publish the returned mechanism, thermal state and all loads together.
    /// The native body heat is already exchanged and must not be applied twice.
    pub fn into_parts(
        self,
    ) -> (
        MaterialDoubleThermalPistonCandidate<S, G>,
        Vec<DgPortLoad>,
        AcousticMotionStep,
        f64,
    ) {
        (
            self.native,
            self.loads,
            self.motion,
            self.total_energy_residual_j,
        )
    }
}

impl<G: PistonGas> MaterialDoubleActingGasPiston<G> {
    /// Apply original reciprocal DG reactions to the same source-bound actuator.
    /// Its canonical solid inertia and both complete gas/port histories survive
    /// the native solve unchanged. Ideal and real-gas material owners are admitted.
    /// Native geometry, clocks, source and energy requirements remain in force.
    /// Any refusal or final cancellation publishes neither mechanism nor loads.
    pub fn acoustic_candidate(
        &self,
        input: DoublePistonInput,
        controls: DoublePistonSolve,
        faces: &[AcousticMotionFace<'_>],
        motion: AcousticMotionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialAcousticPistonCandidate<G>, AcousticMotionError<PistonError>> {
        let candidate =
            self.piston
                .acoustic_candidate(input, controls, faces, motion, &mut cancelled)?;
        let (piston, loads, step, motion, total_energy_residual_j) = candidate.into_parts();
        let piston = self
            .bind_piston_candidate(piston)
            .map_err(AcousticMotionError::Mechanical)?;
        let candidate = MaterialAcousticPistonCandidate {
            piston,
            loads,
            step,
            motion,
            total_energy_residual_j,
        };
        if cancelled() {
            return Err(AcousticMotionError::Cancelled);
        }
        Ok(candidate)
    }

    /// Retain original solid provenance while finite-body heat drives both gases,
    /// their common motion and additional native acoustic domains together.
    /// Each callback sees both actual gas midpoints and the shared body midpoint.
    /// The original thermal/acoustic owner prepares the reaction once and marches
    /// each load once after convergence; this wrapper only retains the rigid body.
    /// The thermal store may differ from the inertial specimen. Heating it never
    /// silently changes the specimen's density, geometry or canonical inertia.
    #[allow(clippy::too_many_arguments)]
    pub fn acoustic_thermal_candidate<S, E>(
        &self,
        faces: &[AcousticMotionFace<'_>],
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        exchange: DoubleHeatExchangeSolve,
        acoustic: AcousticMotionSolve,
        resolve: impl FnMut(
            &[G::State; 2],
            EquilibriumPhaseState,
        ) -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MaterialAcousticThermalPistonCandidate<S, G>,
        AcousticMotionError<ThermalPistonError<E>>,
    > {
        let candidate = self.piston.acoustic_thermal_candidate(
            faces,
            input,
            body,
            controls,
            exchange,
            acoustic,
            resolve,
            &mut cancelled,
        )?;
        let (native, loads, motion, total_energy_residual_j) = candidate.into_parts();
        let (piston, thermal, step) = native.into_parts();
        let piston = self.bind_piston_candidate(piston).map_err(|error| {
            AcousticMotionError::Mechanical(ThermalPistonError::<E>::from(error))
        })?;
        let native = MaterialDoubleThermalPistonCandidate {
            piston,
            thermal,
            step,
        };
        let candidate = MaterialAcousticThermalPistonCandidate {
            native,
            loads,
            motion,
            total_energy_residual_j,
        };
        if cancelled() {
            return Err(AcousticMotionError::Cancelled);
        }
        Ok(candidate)
    }

    /// Use the existing source-pinned conduction layers on the complete loaded
    /// thermal candidate, retaining the original material body and port history.
    #[allow(clippy::too_many_arguments)]
    pub fn acoustic_thermal_layers_candidate(
        &self,
        faces: &[AcousticMotionFace<'_>],
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        exchange: DoubleHeatExchangeSolve,
        acoustic: AcousticMotionSolve,
        layers: [MaterialConductionLayers<'_>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        MaterialAcousticThermalPistonCandidate<ResolvedMaterialStatePoint, G>,
        AcousticMotionError<ThermalPistonError<LayerConductionError>>,
    >
    where
        G::State: FluidTransport,
    {
        self.acoustic_thermal_candidate(
            faces,
            input,
            body,
            controls,
            exchange,
            acoustic,
            |gases, wall| {
                Ok([
                    layers[0].resolve(&gases[0], wall)?,
                    layers[1].resolve(&gases[1], wall)?,
                ])
            },
            cancelled,
        )
    }
}
