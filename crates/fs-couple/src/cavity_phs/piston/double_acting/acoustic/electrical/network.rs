//! Acoustic reaction, reciprocal electrical work and finite material heat on
//! the original shared mechanism. No acoustic owner advances during a thermal
//! or circuit trial; its affine mean load is eliminated before those solves.

use super::super::{
    AcousticMotionError, AcousticMotionFace, AcousticMotionSolve, AcousticMotionStep, check,
    finish, prepare,
};
use crate::cavity_phs::piston::double_acting::electrical::network::{
    DrivenNetworkCandidate, DrivenNetworkError, ElectromagneticNetworkDrive,
    MaterialDrivenNetworkCandidate,
};
use crate::cavity_phs::piston::double_acting::{DoubleActingGasPiston, DoublePistonInput};
use crate::cavity_phs::piston::material::MaterialDoubleActingGasPiston;
use crate::cavity_phs::piston::thermal::boundary::double_acting::network::DoubleNetworkDrive;
use crate::cavity_phs::piston::thermal::boundary::network::NetworkGasSolve;
use crate::circuit_enthalpy::CircuitEnthalpyError;
use crate::circuit_enthalpy::material::MaterialNetwork;
use crate::dg_port::DgPortLoad;
use crate::material_circuit::MaterialCircuit;
use fs_thermal::enthalpy::DepositionSimulation;

/// One immutable proposal containing the original finite circuit/heat/gas
/// network and every additional acoustic domain. Consume all endpoints once.
pub struct AcousticDrivenNetworkCandidate<'m> {
    native: DrivenNetworkCandidate<'m>,
    loads: Vec<DgPortLoad>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}

impl<'m> AcousticDrivenNetworkCandidate<'m> {
    /// Original cell enthalpy, charges, fluxes, source evaluations and motion.
    pub const fn native(&self) -> &DrivenNetworkCandidate<'m> {
        &self.native
    }

    /// Complete acoustic endpoints in the authored face order.
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }

    /// Actual acoustic pressure reaction, native storage, losses and work.
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }

    /// Whole-system balance against only external electrical, mechanical and
    /// thermal supplies. Joule heat is already inside finite cell enthalpy.
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }

    /// Publish the original network and acoustic histories together. Neither
    /// the native transducer work nor resistor heat may be applied again.
    pub fn into_parts(
        self,
    ) -> (
        DrivenNetworkCandidate<'m>,
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
    /// Close acoustic pressure, native back-emf/current and finite material
    /// heat on the same two-gas trajectory. Each thermal trial resolves the
    /// original temperature-dependent resistance; its current supplies BOTH
    /// electromagnetic force and Joule heat. The acoustic load changes that
    /// same trial motion and therefore its back-emf and gas/interface heats.
    ///
    /// Original circuit, gas and thermal work budgets apply to every nested
    /// proposal. Acoustic mean laws are prepared once and the original DG loads
    /// advance once after the complete native network converges. Nothing is
    /// published on a late source, energy, budget or cancellation refusal.
    ///
    /// The finite network owns resistor losses and both gas-contact transfers;
    /// acoustic losses and terminal work remain in their original accounts.
    /// The global gate uses the smaller native/acoustic energy tolerance.
    /// Zero face area retains the original network path and free acoustic
    /// continuation. Fixed acoustic geometry, held interval flow, frozen linear
    /// transducer factor and the native first-order thermal model still apply.
    /// Acoustic domains are additional physical faces, not either gas chamber.
    #[allow(clippy::too_many_arguments)]
    pub fn electroacoustic_network_candidate<'m>(
        &self,
        faces: &[AcousticMotionFace<'_>],
        binding: &MaterialNetwork<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        inputs: &[f64],
        drive: ElectromagneticNetworkDrive,
        controls: NetworkGasSolve,
        acoustic: AcousticMotionSolve,
        guess_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticDrivenNetworkCandidate<'m>, AcousticMotionError<DrivenNetworkError>> {
        let input = drive.gas.input;
        let prepared = prepare(
            faces,
            input.duration_s,
            self.ports.next_tick() - 1,
            acoustic,
            &mut cancelled,
        )?;
        let loaded_force = if prepared.history_force_n == 0.0 {
            input.external_force_n
        } else {
            input.external_force_n + prepared.history_force_n
        };
        if !loaded_force.is_finite() {
            return Err(AcousticMotionError::Invalid(
                "acoustic history and prescribed force overflowed",
            ));
        }
        let loaded_drive = ElectromagneticNetworkDrive {
            gas: DoubleNetworkDrive {
                input: DoublePistonInput {
                    external_force_n: loaded_force,
                    ..input
                },
                ..drive.gas
            },
            ..drive
        };
        let native = self
            .electromagnetic_network_candidate_with_midpoint_resistance(
                binding,
                thermal,
                target_s,
                circuit,
                inputs,
                loaded_drive,
                controls,
                guess_j,
                prepared.resistance_n_s_m,
                &mut cancelled,
            )
            .map_err(|error| match error {
                CircuitEnthalpyError::Cancelled => AcousticMotionError::Cancelled,
                other => AcousticMotionError::Mechanical(other),
            })?;
        let source = &native.network().report().source;
        let (loads, motion) = finish(
            prepared,
            &source.boundary.step().piston,
            input.external_force_n + source.motion.force_n,
            &mut cancelled,
        )?;
        let total = native.coupled().total_energy_residual_j
            + motion.acoustic_energy_residual_j
            + motion.transfer_residual_j;
        if !total.is_finite()
            || total.abs()
                > controls
                    .total_energy_tolerance_j
                    .min(acoustic.energy_tolerance_j)
        {
            return Err(AcousticMotionError::Invalid(
                "electrical, finite thermal, gas, kinetic and acoustic storage do not close external supply",
            ));
        }
        check(&mut cancelled)?;
        Ok(AcousticDrivenNetworkCandidate {
            native,
            loads,
            motion,
            total_energy_residual_j: total,
        })
    }
}

/// The same complete network and acoustic proposal with original specimen
/// geometry, density evidence and canonical inertia attached to live motion.
pub struct MaterialAcousticDrivenNetworkCandidate<'m> {
    native: MaterialDrivenNetworkCandidate<'m>,
    loads: Vec<DgPortLoad>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}

impl<'m> MaterialAcousticDrivenNetworkCandidate<'m> {
    /// Original finite cells, circuit and source-bound shared mechanism.
    pub const fn native(&self) -> &MaterialDrivenNetworkCandidate<'m> {
        &self.native
    }

    /// Complete acoustic endpoints in the authored face order.
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }

    /// Actual original acoustic pressure and reciprocal work account.
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }

    /// Electrical, thermal, gas, kinetic and acoustic enclosing balance [J].
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }

    /// Publish all original state owners together without reconstructing the
    /// piston or resetting its accepted clock and cumulative port history.
    pub fn into_parts(
        self,
    ) -> (
        MaterialDrivenNetworkCandidate<'m>,
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

impl MaterialDoubleActingGasPiston {
    /// Use the original source-resolved rigid specimen in the same coupled
    /// electrical, finite-heat, gas and acoustic solve. The complete proposed
    /// dynamics are retained, including port clocks required by the next
    /// acoustic interval. No numeric replacement mass is accepted.
    ///
    /// The translating specimen remains rigid with frozen source state. The
    /// authored finite thermal cells may represent separate physical bodies;
    /// this adapter does not infer thermal expansion or changing body inertia.
    #[allow(clippy::too_many_arguments)]
    pub fn electroacoustic_network_candidate<'m>(
        &self,
        faces: &[AcousticMotionFace<'_>],
        binding: &MaterialNetwork<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        inputs: &[f64],
        drive: ElectromagneticNetworkDrive,
        controls: NetworkGasSolve,
        acoustic: AcousticMotionSolve,
        guess_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialAcousticDrivenNetworkCandidate<'m>, AcousticMotionError<DrivenNetworkError>>
    {
        let candidate = self.piston().electroacoustic_network_candidate(
            faces,
            binding,
            thermal,
            target_s,
            circuit,
            inputs,
            drive,
            controls,
            acoustic,
            guess_j,
            &mut cancelled,
        )?;
        let (native, loads, motion, total_energy_residual_j) = candidate.into_parts();
        let native = MaterialDrivenNetworkCandidate::from_native(self, native)
            .map_err(AcousticMotionError::Mechanical)?;
        check(&mut cancelled)?;
        Ok(MaterialAcousticDrivenNetworkCandidate {
            native,
            loads,
            motion,
            total_energy_residual_j,
        })
    }
}

#[cfg(test)]
mod tests;
