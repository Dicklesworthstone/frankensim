//! Translating frictionless tools on progressive finite ductile tetrahedra.
//!
//! Native normal complementarity probes the original nonpotential equilibrium
//! or momentum owner at a fixed accepted history and physical interval. Plastic
//! flow, damage, physical losses and kinetic losses remain with those owners.
//! Contact arrival/release has its own bounded endpoint-work defect. It is not
//! plastic heat, fracture energy, an impact law or a complete-failure event.

use super::super::plane_contact::{
    PlaneContactError, TetPlaneEndpoint, TetPlaneGeometry, solve_tet_planes,
};
use super::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileMotionStep, FiniteDuctileTetMotion,
};
use super::*;
use fs_math::dd::Dd;

pub use super::super::plane_contact::{
    PlaneContactControls, PlaneContactReport, RigidPlaneTool, TetPlaneContact,
};

/// Original ductile refusal stays typed, including a rejected damage increment.
pub type FiniteDuctilePlaneContactError = PlaneContactError<FiniteDuctileTetError>;

/// Complete native endpoint and original contact reaction/work record.
#[derive(Debug)]
pub struct FiniteDuctilePlaneContactStep<S> {
    /// Native dead-load work includes contact. Use contact.dead_load_work_j
    /// when separately accounting for tool work and normal contact defects.
    pub step: S,
    pub contact: PlaneContactReport,
}

fn material_balance(energy: &FiniteDuctileTetEnergy, external: Dd, normal_defect: Dd) -> Dd {
    external
        - Dd::from_f64(energy.stored_after_j - energy.stored_before_j)
        - Dd::from_f64(energy.plastic_dissipation_j)
        - Dd::from_f64(energy.damage_release_j)
        - Dd::from_f64(energy.endpoint_work_defect_j)
        - normal_defect
}

impl<L> TetPlaneEndpoint for FiniteDuctileTetStep<'_, '_, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    fn displacement_m(&self) -> &[[f64; 3]] {
        self.solid.displacement_m()
    }
    fn dead_load_work_j(&self) -> f64 {
        self.energy.dead_load_endpoint_work_j
    }
    fn support_work_j(&self) -> f64 {
        self.energy.support_endpoint_work_j
    }
    fn balance(&self, external: Dd, normal_defect: Dd) -> Dd {
        material_balance(&self.energy, external, normal_defect)
    }
}
impl<L> TetPlaneEndpoint for FiniteDuctileMotionStep<'_, '_, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    fn displacement_m(&self) -> &[[f64; 3]] {
        self.motion.solid().displacement_m()
    }
    fn dead_load_work_j(&self) -> f64 {
        self.material_energy.dead_load_endpoint_work_j
    }
    fn support_work_j(&self) -> f64 {
        self.material_energy.support_endpoint_work_j
    }
    fn balance(&self, external: Dd, normal_defect: Dd) -> Dd {
        material_balance(&self.material_energy, external, normal_defect)
            - Dd::from_f64(self.inertia.kinetic_change_j)
            - Dd::from_f64(self.inertia.numerical_loss_j)
    }
}

impl<'g, 'm, L> FiniteDuctileTetSolid<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    fn plane_geometry(&self) -> TetPlaneGeometry<'_> {
        TetPlaneGeometry {
            nodes_m: self.state.nodes_m,
            tetrahedra: self.state.tetrahedra,
            elements: &self.state.elements,
            displacement_m: self.displacement_m(),
        }
    }

    /// Resolve fixed-normal rigid tools while this specimen flows and softens.
    ///
    /// Every native force probe starts from this same accepted specimen and
    /// conditions.duration_s. Candidate rows retain original node/plane IDs;
    /// every vertex is checked against every plane, including omitted rows.
    /// Rejection publishes no changed material, contact or geometry endpoint.
    /// A locally unstable quasistatic branch may require displacement control,
    /// inertia or physical interval refinement; no arclength method is implied.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_plane_contact(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        settings: IncrementSettings,
        material_controls: FinitePlasticStepControls,
        controls: PlaneContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteDuctilePlaneContactStep<FiniteDuctileTetStep<'g, 'm, L>>,
        FiniteDuctilePlaneContactError,
    > {
        let (step, contact) = solve_tet_planes(
            self.plane_geometry(),
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            tools,
            contacts,
            controls,
            |applied| {
                self.advance(
                    prescribed_m,
                    applied,
                    body_force_n_m3,
                    conditions,
                    settings,
                    material_controls,
                    cx,
                )
            },
            cx,
        )?;
        Ok(FiniteDuctilePlaneContactStep { step, contact })
    }
}

impl<'g, 'm, L> FiniteDuctileTetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// The same native reference mass and physical duration drive momentum and
    /// ductile evolution at every contact probe. Numerical kinetic loss remains
    /// independent of material plastic/damage loss and contact-work defects.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_plane_contact(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        settings: IncrementSettings,
        material_controls: FiniteDuctileMotionControls,
        controls: PlaneContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        FiniteDuctilePlaneContactStep<FiniteDuctileMotionStep<'g, 'm, L>>,
        FiniteDuctilePlaneContactError,
    > {
        let (step, contact) = solve_tet_planes(
            self.solid().plane_geometry(),
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            tools,
            contacts,
            controls,
            |applied| {
                self.advance(
                    prescribed_m,
                    applied,
                    body_force_n_m3,
                    conditions,
                    settings,
                    material_controls,
                    cx,
                )
            },
            cx,
        )?;
        Ok(FiniteDuctilePlaneContactStep { step, contact })
    }
}
