//! Original scalar constitutive closure in a prescribed moving conductor.
//!
//! This uses the same nonlinear bulk/contact iteration and contact-refinement
//! owner as stationary response. Material laws see E'=E+u cross B; native
//! laboratory work, conductor-frame heat and signed Lorentz output stay in
//! distinct fields. No accepted nonlinear operating point is rescaled.

use super::*;
use crate::spatial::conduction::motion::{
    SpatialConductionMotion, SpatialMovingConductionSnapshot,
};

/// Actual-amplitude held electrical controls for a moving material snapshot.
/// Current controls include zero-current floating equipotential electrodes.
/// A regulated-power law requires a separate signed motional operating policy.
#[derive(Clone, Copy, Debug)]
pub enum SpatialMovingConductionResponseDrive<'a> {
    Boundary(&'a SpatialConductionBoundary),
    Terminals(ConductionTerminalDrive<'a>),
}

impl<'a> From<&'a SpatialConductionBoundary> for SpatialMovingConductionResponseDrive<'a> {
    fn from(value: &'a SpatialConductionBoundary) -> Self {
        Self::Boundary(value)
    }
}

impl<'a> From<ConductionTerminalDrive<'a>> for SpatialMovingConductionResponseDrive<'a> {
    fn from(value: ConductionTerminalDrive<'a>) -> Self {
        Self::Terminals(value)
    }
}

/// Accepted native motion diagnostics and complete original source closure.
/// Only `field.joule_power_w()` is irreversible heat. Signed terminal input
/// equals that heat plus signed Lorentz mechanical output.
pub type SpatialMovingConductionResponseReport =
    SpatialConductionResponseReport<SpatialMovingConductionSnapshot>;

/// Original nonlinear source owner paired with its accepted moving snapshot.
/// The caller continues any thermal/mechanical history atomically with this
/// candidate; prescribed velocities and B remain external physical inputs.
pub type SpatialMovingConductionResponseCandidate =
    SpatialConductionResponseCandidate<SpatialMovingConductionSnapshot>;

impl SpatialMovingConductionResponseCandidate {
    /// Integrate the Lorentz force from this accepted nonlinear current field
    /// on its own original mesh. No new electrical solve or source query occurs.
    /// Shared vertices accumulate once through the native P1 load integrator.
    /// These loads have not yet advanced a mechanical history.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: crate::spatial::conduction::motion::SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<
        crate::spatial::conduction::motion::SpatialLorentzLoads,
        crate::spatial::conduction::motion::SpatialLorentzLoadError,
    > {
        self.problem.native.moving_snapshot_lorentz_loads(
            &self.report.field,
            selected_cells,
            None,
            limits,
            cancelled,
        )
    }

    /// Require the receiving P1 nodal velocity to reproduce the original
    /// conductor-frame solve's prescribed mean velocity in every selected cell.
    /// The local tolerance and virtual-work semantics are the native moving
    /// field's; an unrelated mesh or reconstructed force is never substituted.
    pub fn lorentz_loads_for_velocity(
        &self,
        selected_cells: &[usize],
        nodal_velocity_m_s: &[[f64; 3]],
        velocity_tolerance_m_s: f64,
        limits: crate::spatial::conduction::motion::SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<
        crate::spatial::conduction::motion::SpatialLorentzLoads,
        crate::spatial::conduction::motion::SpatialLorentzLoadError,
    > {
        self.problem.native.moving_snapshot_lorentz_loads(
            &self.report.field,
            selected_cells,
            Some((nodal_velocity_m_s, velocity_tolerance_m_s)),
            limits,
            cancelled,
        )
    }
}

impl SpatialConductionResponse {
    /// Close the original scalar conductivity/resistivity at the actual
    /// conductor-frame field and current under prescribed motion and B.
    ///
    /// Each constitutive trial solves original charge continuity with its
    /// motional weak load. Direct terminal drives cost one actual-amplitude
    /// sparse solve per trial, including floating and mixed voltage/current
    /// electrodes. Native current, heat, Hodge and signed-work gates remain
    /// binding; positive material values alone do not prove convergence.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_motion_candidate<'d>(
        &self,
        drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: SpatialConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionResponseCandidate> {
        self.solve_with_motion_and_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            drive,
            motion,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// The same moving source closure with externally refreshed original
    /// temperature/phase context, such as a finite-enthalpy endpoint trial.
    ///
    /// Source pins, axes, fixed coordinates, geometry, contact order and
    /// relaxation/refinement policy keep their original owners. Every discarded
    /// nonlinear or integration trial consumes the shared work limits and starts
    /// from this same external context. Only the final driver/field proposal is
    /// returned; no numerical residual is deposited as heat.
    ///
    /// Finite contacts are stationary surface laws: both contact-adjacent
    /// original cells must have zero prescribed velocity. Motion in other
    /// cells is admitted. No sliding-interface work, moving-mesh transport,
    /// induction, field storage or evolution of prescribed mechanics is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_motion_and_materials_and_contacts_candidate<'d>(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<SpatialConductionContact>,
        drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: SpatialConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionResponseCandidate> {
        check(cancelled)?;
        let drive = drive.into();
        let (mut candidate, terminals) = self.solve_coupled_candidate(
            seeds,
            contacts,
            1,
            assembly,
            field_controls,
            controls,
            cancelled,
            |native, cancelled| {
                let (field, terminals) = match drive {
                    SpatialMovingConductionResponseDrive::Boundary(boundary) => (
                        native.solve_with_motion(boundary, motion, field_controls, &mut || {
                            cancelled()
                        })?,
                        None,
                    ),
                    SpatialMovingConductionResponseDrive::Terminals(drive) => {
                        let (field, report) = native.solve_terminals_with_motion(
                            drive.terminals,
                            drive.control,
                            drive.limits,
                            motion,
                            field_controls,
                            &mut || cancelled(),
                        )?;
                        (field, Some(report))
                    }
                };
                // Transfer native buffers into the generic source owner; no
                // secondary field reconstruction or motion diagnostic copy.
                Ok::<_, SpatialConductionResponseError>((field.into_snapshot(), terminals))
            },
        )?;
        candidate.report.terminals = terminals;
        Ok(candidate)
    }
}

/// Private observation seam shared by the two native result types. All field
/// creation remains inside the original stationary or moving solver; external
/// callers cannot substitute an arbitrary constitutive-field adapter.
pub(super) trait ResponseField {
    fn electric_field_v_m(&self) -> &[[f64; 3]];
    fn conductor_field_v_m(&self) -> &[[f64; 3]];
    fn current_density_a_m2(&self) -> &[[f64; 3]];
    fn cell_power_w(&self) -> &[f64];
    fn contacts(&self) -> &[SpatialConductionContactResponse];
    fn nodal_boundary_current_a(&self) -> &[Option<f64>];
    fn power_w(&self) -> f64;
    fn supplied_power_w(&self) -> f64;
    fn cell_mechanical_power_w(&self) -> Option<&[f64]> {
        None
    }
    fn mechanical_power_w(&self) -> Option<f64> {
        None
    }
}

impl ResponseField for SpatialConductionField {
    fn electric_field_v_m(&self) -> &[[f64; 3]] {
        &self.electric_field_v_m
    }
    fn conductor_field_v_m(&self) -> &[[f64; 3]] {
        &self.electric_field_v_m
    }
    fn current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.current_density_a_m2
    }
    fn cell_power_w(&self) -> &[f64] {
        &self.cell_power_w
    }
    fn contacts(&self) -> &[SpatialConductionContactResponse] {
        &self.contacts
    }
    fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.nodal_boundary_current_a
    }
    fn power_w(&self) -> f64 {
        self.power_w
    }
    fn supplied_power_w(&self) -> f64 {
        self.supplied_power_w
    }
}

impl ResponseField for SpatialMovingConductionSnapshot {
    fn electric_field_v_m(&self) -> &[[f64; 3]] {
        SpatialMovingConductionSnapshot::electric_field_v_m(self)
    }
    fn conductor_field_v_m(&self) -> &[[f64; 3]] {
        SpatialMovingConductionSnapshot::conductor_field_v_m(self)
    }
    fn current_density_a_m2(&self) -> &[[f64; 3]] {
        SpatialMovingConductionSnapshot::current_density_a_m2(self)
    }
    fn cell_power_w(&self) -> &[f64] {
        self.cell_joule_power_w()
    }
    fn contacts(&self) -> &[SpatialConductionContactResponse] {
        self.contact_response()
    }
    fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        SpatialMovingConductionSnapshot::nodal_boundary_current_a(self)
    }
    fn power_w(&self) -> f64 {
        self.joule_power_w()
    }
    fn supplied_power_w(&self) -> f64 {
        self.terminal_power_w()
    }
    fn cell_mechanical_power_w(&self) -> Option<&[f64]> {
        Some(SpatialMovingConductionSnapshot::cell_mechanical_power_w(
            self,
        ))
    }
    fn mechanical_power_w(&self) -> Option<f64> {
        Some(SpatialMovingConductionSnapshot::mechanical_power_w(self))
    }
}
