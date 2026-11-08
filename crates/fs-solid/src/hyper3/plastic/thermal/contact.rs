//! Frictionless rigid tooling on the original finite thermal plastic response.
//!
//! Shared native complementarity, original vertex/plane gaps, swept tetrahedron
//! orientation and contact work are owned by the common plane-contact module.
//! Each force probe calls this original thermal/material/momentum endpoint with
//! the same accepted history, temperatures and physical interval.
//!
//! Held contact work on the solid and prescribed tool work remain distinct.
//! Their signed endpoint defect and its noncancelling bound never become heat.
//! This supports translating planes and explicit candidates, not self-contact,
//! generated fracture surfaces, impulses or arbitrary curved tools.

use super::super::plane_contact::{
    PlaneContactError, TetPlaneEndpoint, TetPlaneGeometry, solve_tet_planes,
};
use super::{
    AffineThermalFiniteJ2, FiniteThermalJ2Material, FiniteThermoplasticControls,
    FiniteThermoplasticTetError, FiniteThermoplasticTetSolid, FiniteThermoplasticTetStep,
    HyperTetSettings,
};
use fs_exec::Cx;
use fs_math::dd::Dd;

/// Original tooling primitives, retained at their existing public import path.
pub use super::super::plane_contact::{
    PlaneContactControls, PlaneContactReport, RigidPlaneTool, TetPlaneContact,
};

/// Native thermal/material errors remain intact through contact rejection.
pub type FiniteThermoplasticPlaneContactError = PlaneContactError<FiniteThermoplasticTetError>;

/// Complete native candidate together with its contact and work evidence.
#[derive(Debug)]
pub struct FiniteThermoplasticPlaneContactStep<
    'a,
    L: FiniteThermalJ2Material = AffineThermalFiniteJ2,
> {
    /// Native dead-load work includes applied contact forces. Count tool work
    /// with contact.dead_load_work_j to remove that contribution once.
    pub step: FiniteThermoplasticTetStep<'a, L>,
    pub contact: PlaneContactReport,
}

impl<L: FiniteThermalJ2Material> TetPlaneEndpoint for FiniteThermoplasticTetStep<'_, L> {
    fn displacement_m(&self) -> &[[f64; 3]] {
        self.solid.displacement_m()
    }
    fn dead_load_work_j(&self) -> f64 {
        self.energy.dead_load_endpoint_work_j
    }
    fn support_work_j(&self) -> f64 {
        self.energy.support_endpoint_work_j
    }
    fn balance(&self, external_work: Dd, normal_defect: Dd) -> Dd {
        let mut balance = external_work
            - Dd::from_f64(self.energy.internal_energy_change_j)
            - Dd::from_f64(self.energy.caloric_heat_j)
            - Dd::from_f64(self.energy.endpoint_work_defect_j)
            - normal_defect;
        if let Some(inertia) = &self.inertia {
            balance = balance
                - Dd::from_f64(inertia.kinetic_change_j)
                - Dd::from_f64(inertia.numerical_loss_j);
        }
        balance
    }
}

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticTetSolid<'a, L> {
    /// Propose one physical interval against frictionless translating tools.
    ///
    /// All probes retain the original history, source families, physical time
    /// and reference mass. The final response must match the original contact
    /// receipt before full-mesh gap, swept-orientation and work admission.
    /// Contact arrival/release and temporal accuracy require physical refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_plane_contact(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        settings: HyperTetSettings,
        material_controls: FiniteThermoplasticControls,
        controls: PlaneContactControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticPlaneContactStep<'a, L>, FiniteThermoplasticPlaneContactError>
    {
        let (step, contact) = solve_tet_planes(
            TetPlaneGeometry {
                nodes_m: self.state.nodes_m,
                tetrahedra: self.state.tetrahedra,
                elements: &self.state.elements,
                displacement_m: self.displacement_m(),
            },
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            tools,
            contacts,
            controls,
            |applied| {
                self.advance(
                    new_temperatures_k,
                    prescribed_m,
                    applied,
                    body_force_n_m3,
                    dt_s,
                    settings,
                    material_controls,
                    cx,
                )
            },
            cx,
        )?;
        Ok(FiniteThermoplasticPlaneContactStep { step, contact })
    }
}
