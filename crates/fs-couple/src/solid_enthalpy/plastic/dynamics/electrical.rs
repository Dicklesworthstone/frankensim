//! One electrical, thermal, pressure-contact and inertial material trajectory.
//!
//! This adapter selects inertia in the original coupled residual. Circuit
//! charge/flux, endpoint material conductivity, native solid history and local
//! caloric energy keep their original owners and publish as one candidate.
//! Joule, plastic and reversible heat are deposited once; neither constitutive
//! nor inertial numerical loss heats the material.

use super::super::contact::ThermalSupportContact;
use super::super::electrical::{
    ElectricalSolidError, ElectricalSolidReport, ElectricalSolidSolve, ElectricalSupportContact,
};
use super::{ThermoplasticDynamics, ThermoplasticDynamicsReport};
use crate::circuit_enthalpy::contacts::heat::{HeatAssembly, HeatNetworkSample};
use crate::material_circuit::MaterialCircuit;
use fs_exec::Cx;
use fs_solid::plastic3::IncrementSettings;
use fs_thermal::enthalpy::DepositionSimulation;

/// Accepted field storage, finite heat storage and solid motion together.
pub struct DynamicElectricalCandidate<'m, 't, S = HeatNetworkSample> {
    thermal: DepositionSimulation<'t>,
    dynamics: ThermoplasticDynamics<'m>,
    circuit: MaterialCircuit,
    report: ThermoplasticDynamicsReport<ElectricalSolidReport<S>>,
}
impl<'m, 't, S> DynamicElectricalCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn dynamics(&self) -> &ThermoplasticDynamics<'m> {
        &self.dynamics
    }
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    pub const fn report(&self) -> &ThermoplasticDynamicsReport<ElectricalSolidReport<S>> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticDynamics<'m>,
        MaterialCircuit,
        ThermoplasticDynamicsReport<ElectricalSolidReport<S>>,
    ) {
        (self.thermal, self.dynamics, self.circuit, self.report)
    }
}

impl<'m> ThermoplasticDynamics<'m> {
    /// Finite charge/flux and source-resolved Joule heating on the same interval
    /// as actual nodal inertia and thermal/plastic feedback. The original
    /// thermal geometry remains fixed under the small-strain model. Bulk paths
    /// bound by `ThermoplasticSolid::with_deforming_conductors` follow this same
    /// inertial endpoint; this adds no motional EMF, Lorentz force or changing C/L.
    /// Complete heat assemblies retain finite coolant storage, reservoir powers
    /// and their material evaluations in the same residual and accepted report.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_electrical_candidate<'t, A: HeatAssembly + ?Sized>(
        &self,
        thermal: &DepositionSimulation<'t>,
        circuit: &MaterialCircuit,
        assembly: &A,
        ports: &[f64],
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: ElectricalSolidSolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicElectricalCandidate<'m, 't, A::Sample>, ElectricalSolidError<A::Error>> {
        self.advance_electrical_with_pressure_contacts_candidate(
            thermal,
            circuit,
            assembly,
            &[],
            &[],
            ports,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
            guess_j,
            cx,
        )
    }

    /// Add pressure-dependent thermal supports while retaining electrical and
    /// inertial histories. Bulk and interface resistance use the original
    /// series composition and the stress from this same mechanical endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_electrical_with_supports_candidate<'t, A: HeatAssembly + ?Sized>(
        &self,
        thermal: &DepositionSimulation<'t>,
        circuit: &MaterialCircuit,
        assembly: &A,
        supports: &[ThermalSupportContact<'_>],
        ports: &[f64],
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: ElectricalSolidSolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicElectricalCandidate<'m, 't, A::Sample>, ElectricalSolidError<A::Error>> {
        self.advance_electrical_with_pressure_contacts_candidate(
            thermal,
            circuit,
            assembly,
            supports,
            &[],
            ports,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
            guess_j,
            cx,
        )
    }

    /// Jointly solve inertia, local temperature, source resistivity and solved
    /// support-pressure heat/current transfer. Original source pins, finite
    /// reference mass, complete prescribed support faces and heat partitions
    /// remain enforced. Initialize pressure-contact circuit sources through
    /// self.solid().resolve_electrical_supports before the first interval.
    ///
    /// Each nonlinear trial restarts from accepted u/v/J2/charge/flux/enthalpy.
    /// The closed account includes kinetic change and separate numerical loss.
    /// Source exit, cancellation or failed balance returns no partial endpoint.
    /// Resume with all THREE returned owners to preserve their common clock.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_electrical_with_pressure_contacts_candidate<'t, A: HeatAssembly + ?Sized>(
        &self,
        thermal: &DepositionSimulation<'t>,
        circuit: &MaterialCircuit,
        assembly: &A,
        supports: &[ThermalSupportContact<'_>],
        electrical_supports: &[ElectricalSupportContact],
        ports: &[f64],
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: ElectricalSolidSolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicElectricalCandidate<'m, 't, A::Sample>, ElectricalSolidError<A::Error>> {
        let (candidate, motion) = self.solid.advance_electrical_with_motion_candidate(
            thermal,
            circuit,
            assembly,
            supports,
            electrical_supports,
            ports,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
            guess_j,
            Some(&self.velocity_mps),
            cx,
        )?;
        let motion = motion.ok_or(ElectricalSolidError::Invalid(
            "inertial electrical candidate must retain native velocity and kinetic account",
        ))?;
        let (thermal, solid, circuit, material) = candidate.into_parts();
        cx.checkpoint()
            .map_err(|_| ElectricalSolidError::Cancelled)?;
        Ok(DynamicElectricalCandidate {
            thermal,
            circuit,
            dynamics: Self {
                solid,
                velocity_mps: motion.velocity_mps,
            },
            report: ThermoplasticDynamicsReport {
                material,
                inertia: motion.energy,
            },
        })
    }
}
