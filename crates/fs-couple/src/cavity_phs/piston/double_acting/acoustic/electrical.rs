//! Reciprocal circuit, gas, mechanical and acoustic work on one moving body.
//!
//! The original electrical velocity residual sees the exact affine acoustic
//! reaction in every native momentum proposal. Acoustic histories advance only
//! once, after current, back emf and loaded motion have converged together.

/// Finite material heat, electrical feedback and acoustic pressure together.
pub mod network;

use super::{
    AcousticMotionError, AcousticMotionFace, AcousticMotionSolve, AcousticMotionStep, finish,
    prepare,
};
use crate::cavity_phs::piston::PistonError;
use crate::cavity_phs::piston::double_acting::electrical::{
    ElectricalMotionError, ElectricalMotionInput, ElectricalMotionSolve, ElectricalPistonCandidate,
};
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
};
use crate::dg_port::DgPortLoad;

/// Original circuit/two-gas candidate and all continuing acoustic endpoints.
/// Resistor heat remains the native electrical output; use the finite-network
/// path when that heat must enter material cells and change their properties.
pub struct AcousticElectricalCandidate {
    native: ElectricalPistonCandidate,
    loads: Vec<DgPortLoad>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}
impl AcousticElectricalCandidate {
    pub const fn native(&self) -> &ElectricalPistonCandidate {
        &self.native
    }
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }
    /// Circuit, gases, one kinetic store and acoustic storage/loss balance [J].
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }
    /// Publish every endpoint together; neither circuit nor duct is stepped again.
    pub fn into_parts(
        self,
    ) -> (
        ElectricalPistonCandidate,
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

impl DoubleActingGasPiston {
    /// Solve the original electrical transducer against gas AND acoustic loads.
    ///
    /// Each original mean-pressure relation is prepared once. Its history force
    /// and exact velocity resistance enter every mechanical trial in the native
    /// electrical solve, so pressure changes back emf, current and motion within
    /// the same interval. No acoustic marching occurs in electrical stencils.
    /// Each complete DG owner advances once at the accepted displacement flow.
    ///
    /// `input.external_force_n` excludes both new reciprocal reactions. The
    /// native circuit source work, resistor heat, gas storage and actual wall
    /// work are retained; only the two internal transfer pairs are cancelled.
    /// Source refusal, exhausted work, acoustic admission failure or cancellation
    /// returns no physical endpoints. Callers publish all returned owners together.
    ///
    /// The original fixed C/L/topology, frozen linear transducer and uniform gas
    /// assumptions remain. Acoustic faces belong to additional fixed-geometry
    /// perturbation domains, not second copies of either thermodynamic chamber.
    /// This does not infer moving-mesh sound, magnetic field dynamics or heat
    /// deposition into an unauthored thermal body.
    #[allow(clippy::too_many_arguments)]
    pub fn electroacoustic_candidate(
        &self,
        faces: &[AcousticMotionFace<'_>],
        electrical: ElectricalMotionInput<'_>,
        input: DoublePistonInput,
        piston_controls: DoublePistonSolve,
        electrical_controls: ElectricalMotionSolve,
        acoustic_controls: AcousticMotionSolve,
        initial_midpoint_velocity_m_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticElectricalCandidate, AcousticMotionError<ElectricalMotionError<PistonError>>>
    {
        let prepared = prepare(
            faces,
            input.duration_s,
            self.ports.next_tick() - 1,
            acoustic_controls,
            &mut cancelled,
        )?;
        let force = if prepared.history_force_n == 0.0 {
            input.external_force_n
        } else {
            input.external_force_n + prepared.history_force_n
        };
        if !force.is_finite() {
            return Err(AcousticMotionError::Invalid(
                "acoustic history and prescribed force overflowed",
            ));
        }
        let native = match self.electrical_candidate_with_midpoint_resistance(
            electrical,
            DoublePistonInput {
                external_force_n: force,
                ..input
            },
            piston_controls,
            electrical_controls,
            initial_midpoint_velocity_m_s,
            prepared.resistance_n_s_m,
            &mut cancelled,
        ) {
            Ok(candidate) => candidate,
            Err(ElectricalMotionError::Cancelled) => return Err(AcousticMotionError::Cancelled),
            Err(error) => return Err(AcousticMotionError::Mechanical(error)),
        };
        // Finish subtracts every non-acoustic force from actual native work.
        // The accepted electrical force belongs here, never in the DG source.
        let non_acoustic_force = input.external_force_n + native.motion().force_n;
        let (loads, motion) = finish(
            prepared,
            native.piston_step(),
            non_acoustic_force,
            &mut cancelled,
        )?;
        let total = (native.total_energy_residual_j() + motion.acoustic_energy_residual_j)
            + motion.transfer_residual_j;
        let tolerance = electrical_controls
            .energy_tolerance_j
            .min(acoustic_controls.energy_tolerance_j);
        if !total.is_finite() || total.abs() > tolerance {
            return Err(AcousticMotionError::Invalid(
                "circuit, both gases, shared motion and acoustic storage do not close external supply",
            ));
        }
        if cancelled() {
            return Err(AcousticMotionError::Cancelled);
        }
        Ok(AcousticElectricalCandidate {
            native,
            loads,
            motion,
            total_energy_residual_j: total,
        })
    }
}

#[cfg(test)]
mod tests;
