//! Spatial caloric state -> compatible thermal equilibrium -> pressure valve.
//!
//! This adapter joins the existing heat, material, initial-stress and valve
//! owners. It does not introduce a thermal or mechanical integrator. The
//! returned valve can enter `material_assembly::valve` without copying or
//! overriding its material-derived mass, opening, stiffness or pressure port.
//!
//! `motion` carries this finite thermal state into the existing continuing
//! modal trajectory. Reversible thermoelastic heat and moving the pressure
//! footprint still require separate admitted models.

/// Finite spatial heating and continuing native modal motion on one clock.
pub mod motion;

use fs_exec::Cx;
use fs_solid::linear3::{TetMaterialField, TetStaticSolveConfig};
use fs_solid::reduce::thermal::{
    ThermalPressureLoadedValve, ThermalValveControls, reduce_pressure_loaded_valve,
};
use fs_solid::reduce::{OrificeSpec, ReduceError, ValveCardRequest, ValveElementOrder};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit, StepConfig};

use super::{SolidEnthalpy, SolidEnthalpyError, SolidEnthalpyState, SolidSourceObserver, poll};

/// Explicit mechanics and reduction choices on the original P1 reference mesh.
#[derive(Clone, Debug)]
pub struct HeatedValveReduction<'a> {
    /// Sorted homogeneous support DOFs, shared by equilibrium and vibration.
    pub fixed_dofs: &'a [usize],
    /// Opening observation on the accepted displaced nodes, with a fixed plane.
    pub orifice: OrificeSpec<'a>,
    /// Sorted reference dead-force DOFs [N], independent of acoustic pressure.
    pub equilibrium_nodal_forces_n: &'a [(usize, f64)],
    /// Fixed nodal force per pascal [m²], work-conjugate to swept volume.
    pub nodal_load_per_pa_m2: &'a [[f64; 3]],
    /// Requested frequency window [Hz], with the original modal admission.
    pub window_hz: (f64, f64),
    /// Required fraction of the exact opening compliance in the retained modes.
    pub retained_compliance_floor: f64,
    /// Explicitly authored structural loss factor, not measured source loss.
    pub loss_factor: f64,
    /// Provenance of the original geometry/support/pressure declaration.
    pub source_id: &'a str,
    /// Original bounded equilibrium solver controls.
    pub static_solve: TetStaticSolveConfig,
    /// Maximum Euclidean nodal displacement [m] before aperture admission.
    pub maximum_nodal_displacement_m: f64,
}

/// Preserve the actual owner of a failed thermal/material or reduction step.
#[derive(Debug)]
pub enum HeatedValveError {
    Solid(SolidEnthalpyError),
    Reduction(ReduceError),
}
impl core::fmt::Display for HeatedValveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => write!(f, "heated valve: {error}"),
            Self::Reduction(error) => write!(f, "heated valve: {error}"),
        }
    }
}
impl std::error::Error for HeatedValveError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Reduction(error) => Some(error),
        }
    }
}
impl From<SolidEnthalpyError> for HeatedValveError {
    fn from(error: SolidEnthalpyError) -> Self {
        Self::Solid(error)
    }
}
impl From<ReduceError> for HeatedValveError {
    fn from(error: ReduceError) -> Self {
        Self::Reduction(error)
    }
}

impl SolidEnthalpyState<'_> {
    /// Reduce this exact accepted material/temperature state into its actual
    /// pressure valve. Reference mass, per-tet eigenstrain and frame come from
    /// the retained specimen, never from independently supplied temperatures.
    ///
    /// The original reducer checks equilibrium, displaced aperture, complete
    /// tangent stability, modal compliance and pressure projection. It solves
    /// equilibrium once. Loss remains explicit research input. P1 cells do not
    /// imply a quadratic thermal interpolation or a finite-strain formulation.
    pub fn reduce_pressure_valve(
        &self,
        input: &HeatedValveReduction<'_>,
        cx: &Cx<'_>,
    ) -> Result<ThermalPressureLoadedValve, HeatedValveError> {
        poll(cx)?;
        let request = ValveCardRequest {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            element_order: ValveElementOrder::Tet4,
            fixed_dofs: input.fixed_dofs,
            material: TetMaterialField::PerElement(&self.materials),
            orifice: input.orifice.clone(),
            window_hz: input.window_hz,
            retained_compliance_floor: input.retained_compliance_floor,
            loss_factor: input.loss_factor,
            source_id: input.source_id,
        };
        let valve = reduce_pressure_loaded_valve(
            &request,
            self.thermal_strains(),
            input.equilibrium_nodal_forces_n,
            input.nodal_load_per_pa_m2,
            ThermalValveControls {
                solve: input.static_solve,
                mesh_frame: self.mesh_frame,
                maximum_absolute_mandel_strain: self.maximum_absolute_mandel_strain,
                maximum_nodal_displacement_m: input.maximum_nodal_displacement_m,
            },
            cx,
        )?;
        poll(cx)?;
        Ok(valve)
    }
}

/// A complete candidate: no accepted heat debit can outlive a refused valve.
pub struct HeatedValveCandidate<'m, 't> {
    pub thermal: DepositionSimulation<'t>,
    /// Caloric heat/storage balance only; elastic energy is reported separately
    /// by the retained valve equilibrium and is not fed back as fictitious heat.
    pub thermal_audit: EnergyAudit,
    pub state: SolidEnthalpyState<'m>,
    pub valve: ThermalPressureLoadedValve,
}

impl<'m> SolidEnthalpy<'m> {
    /// Advance finite local heat storage, refresh pinned material sources and
    /// reduce the compatible prestressed valve before publishing any candidate.
    /// Every accepted thermal substep/event passes the existing source/solid
    /// phase observer, so a hidden intermediate melt or source exit still fails.
    ///
    /// The original reference nodes, mass, support pattern and dead pressure
    /// footprint remain fixed. Only the admitted equilibrium geometry observes
    /// the aperture. Cancellation, source/domain exit, instability, closed gap
    /// and inadequate reduction leave the borrowed accepted thermal state intact.
    pub fn advance_pressure_valve_to<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        thermal_config: StepConfig,
        input: &HeatedValveReduction<'_>,
        cx: &Cx<'_>,
    ) -> Result<HeatedValveCandidate<'m, 't>, HeatedValveError> {
        poll(cx)?;
        self.resolve(thermal, cx)?;
        let observer = SolidSourceObserver { source: self, cx };
        let (thermal, _, thermal_audit) = thermal
            .advance_observed(&observer, target_s, thermal_config, || {
                cx.checkpoint().map_err(|_| DepositionError::Cancelled)
            })
            .map_err(SolidEnthalpyError::from)?;
        let state = self.resolve(&thermal, cx)?;
        let valve = state.reduce_pressure_valve(input, cx)?;
        poll(cx)?;
        Ok(HeatedValveCandidate {
            thermal,
            thermal_audit,
            state,
            valve,
        })
    }
}
