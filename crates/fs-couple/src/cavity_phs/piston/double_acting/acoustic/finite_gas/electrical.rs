//! Finite regional gas transport on the original reciprocal electrical drive.
//!
//! The existing circuit/momentum/thermal candidate resolves against the initial
//! acoustic medium once. The common finite-gas finisher then pays the actual
//! parameter work of the already-marched loads. No circuit, piston or acoustic
//! history is advanced inside regional gas-work trials.

use super::{FiniteGasAcousticCandidate, FiniteGasAcousticError, FiniteGasAcousticFace};
use super::{finish, prepare_faces};
use crate::cavity_phs::piston::PistonError;
use crate::cavity_phs::piston::double_acting::acoustic::AcousticMotionSolve;
use crate::cavity_phs::piston::double_acting::electrical::network::{
    DrivenNetworkCandidate, DrivenNetworkError, ElectromagneticNetworkDrive,
    MaterialDrivenNetworkCandidate,
};
use crate::cavity_phs::piston::double_acting::electrical::{
    ElectricalMotionError, ElectricalMotionInput, ElectricalMotionSolve, ElectricalPistonCandidate,
};
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
};
use crate::cavity_phs::piston::material::MaterialDoubleActingGasPiston;
use crate::cavity_phs::piston::thermal::boundary::network::NetworkGasSolve;
use crate::circuit_enthalpy::material::MaterialNetwork;
use crate::dg_port::finite_gas::transport::{FiniteGasDgTransportInput, FiniteGasDgTransportSolve};
use crate::material_circuit::MaterialCircuit;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportPaths;
use fs_thermal::enthalpy::DepositionSimulation;
use std::collections::BTreeMap;

impl DoubleActingGasPiston {
    /// Resolve current, back emf, both piston gases and all acoustic reactions,
    /// then transport the separate finite acoustic-region inventories once.
    ///
    /// The current interval uses each original acoustic medium. Its complete
    /// marched load is rebound to the transported/heated finite gas with native
    /// parameter work debited from that gas. The next interval's original
    /// pressure law therefore feeds the changed medium back into current and
    /// motion. This is the native first-order frozen-medium partition.
    ///
    /// Every caller gas input has empty ports: the solved displacement supplies
    /// the only acoustic inlet flow. Native nozzle, diffusion, porous and
    /// membrane transport remain distinct material flows. The regional parcels
    /// are fixed-volume inventories separate from either moving piston chamber;
    /// uncompensated mean-pressure forces still require an authored external
    /// force or balanced reference. Native resistor heat remains an exported
    /// electrical loss here; the finite-network method owns its deposition.
    ///
    /// Source, budget, transport endpoint, work or cancellation refusal returns
    /// no endpoints. Publish the circuit, piston, regional gases and complete
    /// loads together by consuming the returned candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn transported_electroacoustic_candidate<C: FnMut() -> bool>(
        &self,
        faces: &[FiniteGasAcousticFace<'_>],
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        gas_inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
        paths: GasTransportPaths<'_>,
        electrical: ElectricalMotionInput<'_>,
        input: DoublePistonInput,
        piston_controls: DoublePistonSolve,
        electrical_controls: ElectricalMotionSolve,
        acoustic_controls: AcousticMotionSolve,
        transport_controls: FiniteGasDgTransportSolve,
        initial_midpoint_velocity_m_s: f64,
        cancelled: &mut C,
    ) -> Result<
        FiniteGasAcousticCandidate<ElectricalPistonCandidate>,
        FiniteGasAcousticError<ElectricalMotionError<PistonError>>,
    > {
        let acoustic_faces =
            prepare_faces(faces, parcels, gas_inputs, acoustic_controls, cancelled)?;
        let candidate = self
            .electroacoustic_candidate(
                &acoustic_faces,
                electrical,
                input,
                piston_controls,
                electrical_controls,
                acoustic_controls,
                initial_midpoint_velocity_m_s,
                &mut *cancelled,
            )
            .map_err(FiniteGasAcousticError::from)?;
        let (native, marched, motion, original_total) = candidate.into_parts();
        finish(
            faces,
            parcels,
            gas_inputs,
            paths,
            input.duration_s,
            transport_controls,
            AcousticMotionSolve {
                energy_tolerance_j: acoustic_controls
                    .energy_tolerance_j
                    .min(electrical_controls.energy_tolerance_j)
                    .min(piston_controls.piston.energy_tolerance_j),
                ..acoustic_controls
            },
            native,
            marched,
            motion,
            original_total,
            cancelled,
        )
    }

    /// Add the same conservative regional gas evolution to the original finite
    /// electrical/material-heat network. The native thermal solve retains its
    /// temperature-dependent resistance, Joule heat, both gas contacts and one
    /// mechanical trajectory; only its already-marched acoustic endpoints are
    /// refreshed afterward.
    ///
    /// Regional supplied heat is an independent authored supply. The caller
    /// must not count the network's finite cells or either piston gas again in
    /// this parcel graph. Acoustic parameter work and material-transfer enthalpy
    /// are internal accounts, never a second electrical or thermal supply.
    /// All native source histories and the strictest enclosing energy tolerance
    /// are retained. The fixed-volume and first-order partition assumptions of
    /// the ordinary transported electroacoustic candidate apply unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn transported_electroacoustic_network_candidate<'m, C: FnMut() -> bool>(
        &self,
        faces: &[FiniteGasAcousticFace<'_>],
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        gas_inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
        paths: GasTransportPaths<'_>,
        binding: &MaterialNetwork<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        circuit_inputs: &[f64],
        drive: ElectromagneticNetworkDrive,
        controls: NetworkGasSolve,
        acoustic_controls: AcousticMotionSolve,
        transport_controls: FiniteGasDgTransportSolve,
        guess_j: &[f64],
        cancelled: &mut C,
    ) -> Result<
        FiniteGasAcousticCandidate<DrivenNetworkCandidate<'m>>,
        FiniteGasAcousticError<DrivenNetworkError>,
    > {
        let acoustic_faces =
            prepare_faces(faces, parcels, gas_inputs, acoustic_controls, cancelled)?;
        let candidate = self
            .electroacoustic_network_candidate(
                &acoustic_faces,
                binding,
                thermal,
                target_s,
                circuit,
                circuit_inputs,
                drive,
                controls,
                acoustic_controls,
                guess_j,
                &mut *cancelled,
            )
            .map_err(FiniteGasAcousticError::from)?;
        let (native, marched, motion, original_total) = candidate.into_parts();
        finish(
            faces,
            parcels,
            gas_inputs,
            paths,
            drive.gas.input.duration_s,
            transport_controls,
            AcousticMotionSolve {
                energy_tolerance_j: acoustic_controls
                    .energy_tolerance_j
                    .min(controls.total_energy_tolerance_j),
                ..acoustic_controls
            },
            native,
            marched,
            motion,
            original_total,
            cancelled,
        )
    }
}

impl MaterialDoubleActingGasPiston {
    /// Preserve the original source-resolved specimen, canonical inertia and
    /// complete accepted motion clock through finite regional gas transport.
    ///
    /// This uses the same finite circuit/heat network and gas-work finisher as
    /// the numeric-inertia owner. The rigid body's material state remains frozen;
    /// neither regional heat nor network heat implies a second thermal body,
    /// changing geometry or changing inertia.
    #[allow(clippy::too_many_arguments)]
    pub fn transported_electroacoustic_network_candidate<'m, C: FnMut() -> bool>(
        &self,
        faces: &[FiniteGasAcousticFace<'_>],
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        gas_inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
        paths: GasTransportPaths<'_>,
        binding: &MaterialNetwork<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        circuit_inputs: &[f64],
        drive: ElectromagneticNetworkDrive,
        controls: NetworkGasSolve,
        acoustic_controls: AcousticMotionSolve,
        transport_controls: FiniteGasDgTransportSolve,
        guess_j: &[f64],
        cancelled: &mut C,
    ) -> Result<
        FiniteGasAcousticCandidate<MaterialDrivenNetworkCandidate<'m>>,
        FiniteGasAcousticError<DrivenNetworkError>,
    > {
        let acoustic_faces =
            prepare_faces(faces, parcels, gas_inputs, acoustic_controls, cancelled)?;
        let candidate = self
            .electroacoustic_network_candidate(
                &acoustic_faces,
                binding,
                thermal,
                target_s,
                circuit,
                circuit_inputs,
                drive,
                controls,
                acoustic_controls,
                guess_j,
                &mut *cancelled,
            )
            .map_err(FiniteGasAcousticError::from)?;
        let (native, marched, motion, original_total) = candidate.into_parts();
        finish(
            faces,
            parcels,
            gas_inputs,
            paths,
            drive.gas.input.duration_s,
            transport_controls,
            AcousticMotionSolve {
                energy_tolerance_j: acoustic_controls
                    .energy_tolerance_j
                    .min(controls.total_energy_tolerance_j),
                ..acoustic_controls
            },
            native,
            marched,
            motion,
            original_total,
            cancelled,
        )
    }
}

#[cfg(test)]
mod tests;
