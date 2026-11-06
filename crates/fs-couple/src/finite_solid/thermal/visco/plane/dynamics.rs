//! Thermal source trajectories with native in-plane momentum and retained memory.
//!
//! The prepared source field remains bound to the accepted native motion.
//! Original density, consistent mass and velocity survive prescribed-reference
//! events. Temperature updates never refresh the original mass inventory.
//! Reference input, material relaxation loss, endpoint error and backward-Euler
//! kinetic loss retain their distinct native reports and acceptance gates.

use super::{
    ExpansionViscoPlaneError, ExpansionViscoPlaneField, ExpansionViscoPlaneSolid, checkpoint,
};
use fs_exec::Cx;
use fs_solid::hyper2d::visco::dynamics::{
    FiniteViscoMotionControls, FiniteViscoPlaneMotion, InertialEnergy,
};
use fs_solid::hyper2d::visco::{
    FiniteViscoPlaneEnergy, FiniteViscoPlaneReferenceEnergy, FiniteViscoStepControls,
};
use fs_solid::{NewtonReport, NewtonSettings, Patch};

/// Native material histories and momentum bound to the original source owners.
#[derive(Debug)]
pub struct ExpansionViscoPlaneMotion<'a> {
    source: &'a ExpansionViscoPlaneField,
    native: FiniteViscoPlaneMotion<'a>,
}

/// Complete physical interval; no accepted source, motion or history is mutated.
#[derive(Debug)]
pub struct ExpansionViscoPlaneMotionStep<'a> {
    pub motion: ExpansionViscoPlaneMotion<'a>,
    pub equilibrium: NewtonReport,
    pub material_energy: FiniteViscoPlaneEnergy,
    pub inertia: InertialEnergy,
    pub balance_error_j: f64,
    pub reactions_n: Vec<(usize, f64)>,
}

/// Zero-duration source event with unchanged in-plane velocity and closed mass.
#[derive(Debug)]
pub struct ExpansionViscoPlaneMotionReferenceChange<'a> {
    pub motion: ExpansionViscoPlaneMotion<'a>,
    pub energy: FiniteViscoPlaneReferenceEnergy,
}

/// Atomic reference event and subsequent physical interval. Reference work is
/// independent of the material/kinetic accounts of the mechanical interval.
#[derive(Debug)]
pub struct ExpansionViscoPlaneThermalMotionStep<'a> {
    pub motion: ExpansionViscoPlaneMotion<'a>,
    pub reference: FiniteViscoPlaneReferenceEnergy,
    pub material_energy: FiniteViscoPlaneEnergy,
    pub inertia: InertialEnergy,
    pub balance_error_j: f64,
    pub equilibrium: NewtonReport,
    pub reactions_n: Vec<(usize, f64)>,
}

impl ExpansionViscoPlaneThermalMotionStep<'_> {
    /// Count fixed-deformation reference input and exact pressure actuation once.
    pub fn reference_input_j(&self) -> f64 {
        self.reference.reference_work_j + self.reference.constraint_reference_work_j
    }
}

impl<'a> ExpansionViscoPlaneMotion<'a> {
    /// Bind the accepted source solid to explicit original density [kg/m³],
    /// nodal velocity [m/s] and a consistent-mass assembly byte cap. Existing
    /// displacement, material age and every quadrature history are retained.
    /// The density is the original inventory; it is not queried at endpoint T.
    pub fn try_new(
        solid: ExpansionViscoPlaneSolid<'a>,
        reference_density_kg_m3: &[f64],
        initial_velocity_mps: &[[f64; 2]],
        maximum_mass_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, ExpansionViscoPlaneError> {
        let native = FiniteViscoPlaneMotion::try_new(
            solid.native,
            reference_density_kg_m3,
            initial_velocity_mps,
            maximum_mass_assembly_bytes,
            cx,
        )?;
        Ok(Self {
            source: solid.source,
            native,
        })
    }

    pub fn native(&self) -> &FiniteViscoPlaneMotion<'a> {
        &self.native
    }

    pub fn source_field(&self) -> &'a ExpansionViscoPlaneField {
        self.source
    }

    /// Project nodal reference acceleration through the retained consistent
    /// mass. The result is a physical nodal force field [N], including supports.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 2]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, ExpansionViscoPlaneError> {
        self.native
            .reference_acceleration_load_n(acceleration_m_s2, cx)
            .map_err(Into::into)
    }

    /// Advance the native physical interval under absolute endpoint supports
    /// [m], nominal reference tractions [Pa] and explicit nodal forces [N].
    /// Physical relaxation and numerical kinetic loss remain separate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        dirichlet: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        traction: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        nodal_forces_n: &[(usize, f64)],
        settings: NewtonSettings,
        controls: FiniteViscoMotionControls,
        cx: &Cx<'_>,
    ) -> Result<ExpansionViscoPlaneMotionStep<'a>, ExpansionViscoPlaneError> {
        let step = self.native.advance(
            dt_s,
            dirichlet,
            traction,
            nodal_forces_n,
            settings,
            controls,
            cx,
        )?;
        Ok(ExpansionViscoPlaneMotionStep {
            motion: Self {
                source: self.source,
                native: step.motion,
            },
            equilibrium: step.equilibrium,
            material_energy: step.material_energy,
            inertia: step.inertia,
            balance_error_j: step.balance_error_j,
            reactions_n: step.reactions_n,
        })
    }

    /// Admit each original element's source/frame/datum/calibration once, then
    /// restore native transverse equilibrium without changing physical time,
    /// in-plane velocity, original density or consistent mass. The native
    /// reference report retains exact pressure-volume actuation separately.
    pub fn change_reference<'b>(
        &self,
        target: &'b ExpansionViscoPlaneField,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<ExpansionViscoPlaneMotionReferenceChange<'b>, ExpansionViscoPlaneError>
    where
        'a: 'b,
    {
        checkpoint(cx)?;
        target.check_transition(
            self.source,
            self.native.solid().reference_mesh().elems.len(),
            cx,
        )?;
        let changed = self
            .native
            .change_reference(target.native_laws(), controls, cx)?;
        Ok(ExpansionViscoPlaneMotionReferenceChange {
            motion: ExpansionViscoPlaneMotion {
                source: target,
                native: changed.motion,
            },
            energy: changed.energy,
        })
    }

    /// Return a successor only after both the source event and actual motion
    /// interval pass every gate, including the late numerical kinetic-loss
    /// allowance. A refusal leaves the original source, mass, velocity and all
    /// material histories intact. The caller refines the prescribed temperature
    /// path and physical interval; no unseen thermal trajectory is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn change_reference_and_advance<'b>(
        &self,
        target: &'b ExpansionViscoPlaneField,
        dt_s: f64,
        dirichlet: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        traction: &[(Patch, &dyn Fn(f64, f64) -> [f64; 2])],
        nodal_forces_n: &[(usize, f64)],
        settings: NewtonSettings,
        reference_controls: FiniteViscoStepControls,
        motion_controls: FiniteViscoMotionControls,
        cx: &Cx<'_>,
    ) -> Result<ExpansionViscoPlaneThermalMotionStep<'b>, ExpansionViscoPlaneError>
    where
        'a: 'b,
    {
        let changed = self.change_reference(target, reference_controls, cx)?;
        let step = changed.motion.advance(
            dt_s,
            dirichlet,
            traction,
            nodal_forces_n,
            settings,
            motion_controls,
            cx,
        )?;
        Ok(ExpansionViscoPlaneThermalMotionStep {
            motion: step.motion,
            reference: changed.energy,
            material_energy: step.material_energy,
            inertia: step.inertia,
            balance_error_j: step.balance_error_j,
            equilibrium: step.equilibrium,
            reactions_n: step.reactions_n,
        })
    }
}
