//! Actual spatial current, reference-mass heat and native inertial solid motion.
//!
//! The original coupled endpoint assembles current on its proposed X+u mesh.
//! Only the existing optional-velocity mechanical selection changes; conductor
//! source refresh, per-tet Joule heat, native plastic history and the complete
//! finite heat assembly still share the same enthalpy residual and candidate.

use super::super::{ThermoplasticDynamics, ThermoplasticDynamicsReport};
use crate::circuit_enthalpy::{
    conduction::spatial::{SpatialConductionDrive, SpatialConductionThermalStep},
    contacts::heat::{HeatAssembly, HeatNetworkSample},
};
use crate::solid_enthalpy::plastic::electrical::{
    ElectricalSolidError, ElectricalSolidReport, spatial::SpatialElectricalSolidSolve,
};
use fs_em::spatial::conduction::SpatialConduction;
use fs_exec::Cx;
use fs_solid::plastic3::IncrementSettings;
use fs_thermal::enthalpy::DepositionSimulation;

/// Accepted heat, full nodal motion and actual current operator on one clock.
pub struct DynamicSpatialElectricalCandidate<'m, 't, S = HeatNetworkSample> {
    thermal: DepositionSimulation<'t>,
    dynamics: ThermoplasticDynamics<'m>,
    problem: SpatialConduction,
    report: ThermoplasticDynamicsReport<ElectricalSolidReport<S, SpatialConductionThermalStep>>,
}
impl<'m, 't, S> DynamicSpatialElectricalCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn dynamics(&self) -> &ThermoplasticDynamics<'m> {
        &self.dynamics
    }
    pub const fn problem(&self) -> &SpatialConduction {
        &self.problem
    }
    pub const fn report(
        &self,
    ) -> &ThermoplasticDynamicsReport<ElectricalSolidReport<S, SpatialConductionThermalStep>> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticDynamics<'m>,
        SpatialConduction,
        ThermoplasticDynamicsReport<ElectricalSolidReport<S, SpatialConductionThermalStep>>,
    ) {
        (self.thermal, self.dynamics, self.problem, self.report)
    }
}

impl<'m> ThermoplasticDynamics<'m> {
    /// Advance actual native current and finite heat with retained u/v/J2 history.
    /// Every nonlinear trial starts from the same accepted mechanical state,
    /// solves inertia at the proposed local temperatures, and rebuilds the
    /// original tetrahedral electrical operator on THAT trial's X+u geometry.
    /// Original E.J heats the reference-mass cells once. Kinetic energy and
    /// inertial backward-Euler loss stay in their independent accounts; neither
    /// becomes fictitious material heating.
    ///
    /// The supplied voltage trace follows its original material vertices;
    /// prescribed current density acts on the actual deformed face geometry.
    /// Shared total-current control instead holds the signed TOTAL current
    /// through moving equipotential terminal faces. Each source/geometry trial
    /// solves their voltage and actual current distribution; supplied work is
    /// I times the solved voltage drop, independently of the Joule integral.
    /// This is quasistatic fully-solid isotropic DC transport, within the native
    /// small-strain/rotation limits. It adds no motional EMF, Lorentz force,
    /// electrical contact motion, field storage or evolving thermal geometry.
    /// Boundary work remains with the native electrical owner. Resolve load
    /// pulses and refine intervals independently of thermal nonlinear tolerance.
    ///
    /// Retain all returned owners for continuation. Source/domain, budget,
    /// energy, or cancellation refusal publishes none of the proposed histories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_spatial_electrical_candidate<'t, 'd, A: HeatAssembly + ?Sized>(
        &self,
        thermal: &DepositionSimulation<'t>,
        problem: &SpatialConduction,
        assembly: &A,
        drive: impl Into<SpatialConductionDrive<'d>>,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: SpatialElectricalSolidSolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicSpatialElectricalCandidate<'m, 't, A::Sample>, ElectricalSolidError<A::Error>>
    {
        let (candidate, motion) = self
            .solid
            .advance_spatial_electrical_with_motion_candidate(
                thermal,
                problem,
                assembly,
                drive,
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
            "inertial spatial-current candidate must retain native velocity and kinetic account",
        ))?;
        let (thermal, solid, problem, material) = candidate.into_parts();
        cx.checkpoint()
            .map_err(|_| ElectricalSolidError::Cancelled)?;
        Ok(DynamicSpatialElectricalCandidate {
            thermal,
            problem,
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
