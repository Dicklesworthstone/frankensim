//! Sustained rigid-tool forming and reciprocal heat on the original owners.
//!
//! Each thermal proposal solves the original normal-contact complementarity
//! problem through the same finite thermoplastic equilibrium and accepted
//! material/momentum histories. Tool motion is affine over the actual physical
//! interval. There is no penalty stiffness or impulse-to-force conversion.
//!
//! Only native material caloric transfer enters the thermal cells. Prescribed
//! tool work, body-contact work and the signed normal endpoint-work defect stay
//! explicit; the latter and numerical kinetic loss are never thermal sources.
//! Native full-plane and affine tetrahedron path admission remains mandatory.
//! Friction, tool rotation, compliant tooling and evolving thermal conductance
//! require their own admitted models and are not inferred by this adapter.

use super::feedback::{self, MechanicalEnergy, MechanicalProposal};
use super::{
    FiniteThermoplasticCandidate, FiniteThermoplasticConfig, FiniteThermoplasticControls,
    FiniteThermoplasticError, FiniteThermoplasticMaterial, FiniteThermoplasticMaterialFamily,
    FiniteThermoplasticSolid, HyperTetSettings, HyperTetSolution, finite,
};
use crate::iqn_ils::driver::CouplingMethod;
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_solid::hyper3::plastic::thermal::contact::{
    PlaneContactControls, PlaneContactReport, RigidPlaneTool, TetPlaneContact,
};
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticTetEnergy, FiniteThermoplasticTetSolid,
};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

/// Original material storage/heat and the separately admitted tool work.
#[derive(Clone, Debug)]
pub struct FiniteThermoplasticContactEnergy {
    /// Native account: its held nodal-load work includes body-contact work.
    /// `contact.dead_load_work_j` removes that contribution exactly once for
    /// the combined external balance.
    pub material: FiniteThermoplasticTetEnergy,
    pub contact: PlaneContactReport,
    /// Material endpoint defect plus the signed tool-minus-body work defect.
    /// Both original absolute gates must pass before this sum is formed.
    pub external_work_defect_j: f64,
}

impl MechanicalEnergy for FiniteThermoplasticContactEnergy {
    fn internal_energy_change_j(&self) -> f64 {
        self.material.internal_energy_change_j
    }
    fn endpoint_work_defect_j(&self) -> f64 {
        self.external_work_defect_j
    }
    fn dead_load_endpoint_work_j(&self) -> f64 {
        self.contact.dead_load_work_j
    }
    fn support_endpoint_work_j(&self) -> f64 {
        self.material.support_endpoint_work_j
    }
    fn prescribed_tool_work_j(&self) -> f64 {
        self.contact.tool_work_j
    }
    fn nonthermal_loss_j(&self) -> f64 {
        0.0
    }
}

/// Atomic solid, thermal and observer endpoint with its original contact solve.
pub type FiniteThermoplasticContactCandidate<'a, 't, O, M = FiniteThermoplasticMaterial<'a>> =
    FiniteThermoplasticCandidate<'a, 't, O, FiniteThermoplasticContactEnergy, M>;

impl<'a, M: FiniteThermoplasticMaterialFamily> FiniteThermoplasticSolid<'a, M> {
    /// Form against explicitly translating rigid planes while solving heat.
    ///
    /// Contact row scales condition the native natural map in N/m; they are
    /// not material stiffnesses. The caller supplies complete contact rows and
    /// numerical force, gap, path and work allowances. Refine physical time to
    /// resolve contact onset/release and material-rate effects independently of
    /// the within-interval thermal iteration.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_plane_contact<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        contact_controls: PlaneContactControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticContactCandidate<'a, 't, O, M>, FiniteThermoplasticError<O::Error>>
    {
        self.advance_with_plane_contact_and_coupling(
            thermal,
            observer,
            target_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            tools,
            contacts,
            settings,
            mechanical_controls,
            contact_controls,
            thermal_config,
            config,
            CouplingMethod::RelaxedPicard,
            cx,
        )
    }

    /// Use the original Picard or IQN-ILS thermal iteration during contact.
    ///
    /// Every force probe and every temperature replay starts at the same
    /// accepted deformation, plastic metric, temperature and optional velocity.
    /// Final mechanics is re-solved at the actual returned caloric temperature;
    /// temperature, per-cell heat and the complete energy balance all gate the
    /// candidate. Failure or cancellation publishes none of its three owners.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_plane_contact_and_coupling<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        contact_controls: PlaneContactControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticContactCandidate<'a, 't, O, M>, FiniteThermoplasticError<O::Error>>
    {
        let next = feedback::advance_with_coupling(
            self.binding(),
            self.time_s,
            thermal,
            observer,
            target_s,
            thermal_config,
            config,
            method,
            cx,
            |temperatures, duration_s| {
                self.propose_plane_contact(
                    temperatures,
                    duration_s,
                    prescribed_m,
                    nodal_forces_n,
                    body_force_n_m3,
                    tools,
                    contacts,
                    settings,
                    mechanical_controls,
                    contact_controls,
                    cx,
                )
            },
        )?;
        Ok(FiniteThermoplasticCandidate {
            thermal: next.thermal,
            solid: Self {
                solid: next.state,
                cells: next.cells,
                time_s: target_s,
            },
            observer: next.observer,
            report: next.report,
        })
    }

    /// Shared proposal for fixed-observer and current-conduction compositions.
    /// Only the original material heat is transferred, whichever transport is selected.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn propose_plane_contact<E>(
        &self,
        temperatures: &[f64],
        duration_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        contact_controls: PlaneContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        MechanicalProposal<
            FiniteThermoplasticTetSolid<'a, M>,
            HyperTetSolution,
            FiniteThermoplasticContactEnergy,
        >,
        FiniteThermoplasticError<E>,
    > {
        let formed = self
            .solid
            .advance_with_plane_contact(
                temperatures,
                prescribed_m,
                nodal_forces_n,
                body_force_n_m3,
                duration_s,
                tools,
                contacts,
                settings,
                mechanical_controls,
                contact_controls,
                cx,
            )
            .map_err(FiniteThermoplasticError::Contact)?;
        let step = formed.step;
        let heat_j = step
            .element_energy
            .iter()
            .enumerate()
            .map(|(i, energy)| {
                finite(
                    self.solid
                        .reference_volume_m3(i)
                        .expect("retained tetrahedron index")
                        * energy.caloric_heat_j_m3(),
                )
            })
            .collect::<Result<_, FiniteThermoplasticError<E>>>()?;
        let external_work_defect_j = finite(
            (Dd::from_f64(step.energy.endpoint_work_defect_j)
                + Dd::from_f64(formed.contact.normal_work_defect_j))
            .to_f64(),
        )?;
        Ok(MechanicalProposal {
            state: step.solid,
            equilibrium: step.equilibrium,
            energy: FiniteThermoplasticContactEnergy {
                material: step.energy,
                contact: formed.contact,
                external_work_defect_j,
            },
            inertia: step.inertia,
            heat_j,
        })
    }
}
