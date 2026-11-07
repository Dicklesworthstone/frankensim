//! Deformation-dependent normal conduction through the original heat timeline.
//!
//! The geometric callback sees the SAME proposed native finite solid that
//! supplies physical plastic/reversible heat. Its complete current areas and
//! normal distances refresh the existing conduction observer before every
//! thermal replay. Reference cells, original source laws and mass never move.
//!
//! This is an explicit, first-order endpoint-geometry composition of finite
//! mechanics with the existing normal two-point resistance model. The caller
//! must derive admissible faces/paths from solved displacement and justify the
//! normal-flow reduction. No arbitrary distorted-tet Fourier discretization,
//! contact opening, remeshing, area reconstruction or thermal anisotropy is
//! inferred by supplying scalar metrics. All other boundary/source channels
//! retain their authored original meanings.

use fs_exec::Cx;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetSolid,
};
use fs_solid::hyper3::{HyperTetSettings, HyperTetSolution};
use fs_thermal::enthalpy::conduction::{ConductionObserver, CurrentConductionGeometry};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, StepConfig};

use crate::iqn_ils::driver::CouplingMethod;

#[cfg(feature = "experimental-solid-contact-heat")]
use super::contact::FiniteThermoplasticContactCandidate;
#[cfg(feature = "experimental-solid-contact-heat")]
use fs_solid::hyper3::plastic::thermal::contact::{
    PlaneContactControls, RigidPlaneTool, TetPlaneContact,
};

use super::{
    FiniteThermoplasticCandidate, FiniteThermoplasticConfig, FiniteThermoplasticError,
    FiniteThermoplasticMaterial, FiniteThermoplasticSolid, feedback, finite,
};

impl<'a> FiniteThermoplasticSolid<'a> {
    /// Solve material history, current normal conductances and caloric feedback.
    ///
    /// `geometry` must be a pure evaluation of the supplied proposed native
    /// solid and immutable geometric bindings. Every call begins from the same
    /// accepted history. Geometry errors, source-domain exits and cancellation
    /// withhold the entire candidate; external callback side effects cannot be
    /// rolled back. Current metrics must cover every original observer channel.
    /// The original reference cell shape and mass stay unchanged and inspectable.
    ///
    /// The native feedback owner re-solves mechanics at the actual returned
    /// caloric temperatures. Besides its existing per-cell heat, temperature
    /// and whole-energy gates, ANY area/distance discrepancy between the used
    /// transport and that actual mechanical endpoint must fit the explicit
    /// relative allowance in `[0,1)`. The returned second value is that largest
    /// discrepancy. It is a coupling residual, not a continuum-error bound.
    ///
    /// The returned observer retains the metrics actually used for heat. Its
    /// conductivity laws still refresh on every native thermal substep; actual
    /// bath events, source watts, caloric clocks and CFL remain original-owned.
    /// Advancing another physical interval invokes `geometry` afresh. Physical
    /// interval, thermal substep and finite-cell resolution need independent
    /// refinement. The plane-contact companion also composes the original
    /// sustained rigid-tool work with these same current thermal metrics.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_deformed_conduction<'t, 'k>(
        &self,
        thermal: &DepositionSimulation<'t>,
        conduction: &ConductionObserver<'k>,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        mechanical_controls: FiniteThermoplasticControls,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        maximum_geometry_relative_defect: f64,
        geometry: impl FnMut(
            &FiniteThermoplasticTetSolid<'a, FiniteThermoplasticMaterial<'a>>,
            &Cx<'_>,
        ) -> Result<CurrentConductionGeometry, DepositionError>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            FiniteThermoplasticCandidate<'a, 't, ConductionObserver<'k>>,
            f64,
        ),
        FiniteThermoplasticError<DepositionError>,
    > {
        self.advance_with_deformed_transport(
            thermal,
            conduction,
            target_s,
            thermal_config,
            config,
            method,
            maximum_geometry_relative_defect,
            |temperatures, duration_s| {
                let step = self
                    .solid
                    .advance(
                        temperatures,
                        prescribed_m,
                        nodal_forces_n,
                        body_force_n_m3,
                        duration_s,
                        settings,
                        mechanical_controls,
                        cx,
                    )
                    .map_err(FiniteThermoplasticError::Solid)?;
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
                    .collect::<Result<_, FiniteThermoplasticError<DepositionError>>>()?;
                Ok(feedback::MechanicalProposal {
                    state: step.solid,
                    equilibrium: step.equilibrium,
                    energy: step.energy,
                    inertia: step.inertia,
                    heat_j,
                })
            },
            geometry,
            cx,
        )
    }

    /// Form against translating rigid planes while current normal heat paths
    /// respond to that SAME solved deformation. Native complementarity, affine
    /// path admission, actual tool work, retained material histories and all
    /// thermal/geometry gates compose through their existing owners. Friction
    /// and changing contact topology are not inferred from conductive metrics.
    #[cfg(feature = "experimental-solid-contact-heat")]
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_plane_contact_and_deformed_conduction<'t, 'k>(
        &self,
        thermal: &DepositionSimulation<'t>,
        conduction: &ConductionObserver<'k>,
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
        maximum_geometry_relative_defect: f64,
        geometry: impl FnMut(
            &FiniteThermoplasticTetSolid<'a, FiniteThermoplasticMaterial<'a>>,
            &Cx<'_>,
        ) -> Result<CurrentConductionGeometry, DepositionError>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            FiniteThermoplasticContactCandidate<'a, 't, ConductionObserver<'k>>,
            f64,
        ),
        FiniteThermoplasticError<DepositionError>,
    > {
        self.advance_with_deformed_transport(
            thermal,
            conduction,
            target_s,
            thermal_config,
            config,
            method,
            maximum_geometry_relative_defect,
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
            geometry,
            cx,
        )
    }

    // Only the original mechanical proposal differs between ordinary forming
    // and contact forming; their geometry preparation and heat replay are one.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn advance_with_deformed_transport<'t, 'k, A: feedback::MechanicalEnergy>(
        &self,
        thermal: &DepositionSimulation<'t>,
        conduction: &ConductionObserver<'k>,
        target_s: f64,
        thermal_config: StepConfig,
        config: FiniteThermoplasticConfig,
        method: CouplingMethod,
        maximum_geometry_relative_defect: f64,
        mechanics: impl FnMut(
            &[f64],
            f64,
        ) -> Result<
            feedback::MechanicalProposal<
                FiniteThermoplasticTetSolid<'a, FiniteThermoplasticMaterial<'a>>,
                HyperTetSolution,
                A,
            >,
            FiniteThermoplasticError<DepositionError>,
        >,
        mut geometry: impl FnMut(
            &FiniteThermoplasticTetSolid<'a, FiniteThermoplasticMaterial<'a>>,
            &Cx<'_>,
        ) -> Result<CurrentConductionGeometry, DepositionError>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            FiniteThermoplasticCandidate<'a, 't, ConductionObserver<'k>, A>,
            f64,
        ),
        FiniteThermoplasticError<DepositionError>,
    > {
        let (next, geometry_relative_defect) = feedback::advance_with_transport(
            self.binding(),
            self.time_s,
            thermal,
            conduction,
            target_s,
            thermal_config,
            config,
            method,
            Some(maximum_geometry_relative_defect),
            cx,
            mechanics,
            |original, solid| {
                let metrics = geometry(solid, cx)?;
                original
                    .with_current_geometry(thermal, metrics, || {
                        cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                    })
                    .map_err(FiniteThermoplasticError::from)
            },
            |used, actual| {
                let used = used
                    .current_geometry()
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "missing current conduction geometry",
                    ))?;
                let actual = actual
                    .current_geometry()
                    .ok_or(FiniteThermoplasticError::Invalid(
                        "missing actual-endpoint conduction geometry",
                    ))?;
                used.relative_difference(actual, || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })
                .map_err(FiniteThermoplasticError::from)
            },
        )?;
        Ok((
            FiniteThermoplasticCandidate {
                thermal: next.thermal,
                solid: Self {
                    solid: next.state,
                    cells: next.cells,
                    time_s: target_s,
                },
                observer: next.observer,
                report: next.report,
            },
            geometry_relative_defect,
        ))
    }
}
