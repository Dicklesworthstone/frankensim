//! Reciprocal spatial heat on the original full tetrahedral trajectory.
//!
//! Native average-acceleration Newmark keeps displacement, velocity, acceleration,
//! reference mass and the accepted step budget. Original nodal backward Euler
//! retains heat transport, source charts, surface storage and boundary reactions.
//! Every material exchange uses the same enthalpy-centroid quadrature and
//! polarized source potential as the quasistatic owner. Kinetic energy is storage;
//! neither a numerical residual nor supplied Rayleigh loss becomes nodal heat.
//!
//! The endpoint coupling is generally first order and needs time/mesh refinement.
//! Centroid material temperatures do not resolve within-element thermoelastic
//! gradients or establish exact nodal entropy production. Supports remain fixed
//! and zero, transport geometry remains the original reference geometry, and
//! every corner must stay in the admitted solid source interval.

use super::*;
use fs_solid::loads3::transient::response::TransientResponse3;

/// Physical destination of the independently supplied native Rayleigh loss.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpatialRayleighDestination {
    /// Energy leaves this modeled solid/thermal system. No nodal deposition is
    /// inferred from a global Rayleigh matrix. The supplied damping must exclude
    /// thermoelastic loss already resolved by this heat/mechanics coupling.
    ExternalSink,
}

/// The accepted nodal heat and the SAME retained native inertial checkpoint.
pub struct DynamicReciprocalSpatialSolid<'s, 'a, 'm, 'c, 'g, 'mesh> {
    source: &'s SpatialSolidEnthalpy<'a, 'm, 'c>,
    state: SpatialSolidState<'s>,
    laws: Arc<[AffineThermoelastic]>,
    specific_h: Vec<f64>,
    fixed_dofs: Vec<usize>,
    minimum_nodal_capacity_j_k: Arc<[f64]>,
    motion: TetTransientCheckpoint3<'g, 'mesh>,
    response: TransientResponse3,
}

/// Native force/mass audit and independent constitutive/caloric accounts.
#[derive(Debug)]
pub struct DynamicSpatialThermoelasticReport {
    pub iterations: usize,
    pub native: ThermalTransientStep3,
    pub elements: Vec<ThermoelasticExchange>,
    pub nodal_reversible_heat_j: Vec<Dd>,
    pub kinetic_energy_change_j: f64,
    pub stress_work_j: f64,
    pub external_work_j: f64,
    /// Supplied physical Rayleigh loss exported to the declared external sink.
    pub damping_export_j: f64,
    /// Delta K + endpoint stress work + exported damping - actual native work.
    pub mechanical_residual_j: f64,
    pub internal_energy_change_j: f64,
    pub reversible_heat_j: f64,
    /// Original volume/boundary/radiative heat and prescribed-T reactions.
    pub external_thermal_heat_j: f64,
    pub heat_defect_j: f64,
    /// Delta caloric storage + Delta U + Delta K + damping - external inputs.
    pub total_energy_residual_j: f64,
}

/// A failed or cancelled proposal returns neither heat nor a motion checkpoint.
/// The caller publishes this complete pair; the supplied accepted owner remains
/// unchanged and can be retried or retained as the accepted prefix.
pub struct DynamicSpatialThermoelasticCandidate<'s, 'a, 'm, 'c, 'g, 'mesh> {
    pub solid: DynamicReciprocalSpatialSolid<'s, 'a, 'm, 'c, 'g, 'mesh>,
    pub thermal: NodalHeatStep,
    pub report: DynamicSpatialThermoelasticReport,
}

impl<'s, 'a: 's, 'm: 'a, 'c: 'a> ReciprocalSpatialSolid<'s, 'a, 'm, 'c> {
    /// Launch real inertia from the already admitted loaded equilibrium and
    /// supplied nodal velocity. Initialization does not simulate a heat step.
    /// The original native P1 density profile supplies consistent mass; neither
    /// current constitutive density nor a scalar average replaces that inventory.
    ///
    /// Loads must be the same assembled nodal forces as this equilibrium.
    /// Subsequent load changes use the native endpoint work convention.
    pub fn begin_transient<'g, 'mesh>(
        self,
        loads: &TetLoadCase3<'g, 'mesh>,
        velocity_m_s: &[[f64; 3]],
        controls: TransientConfig3,
        _damping_destination: SpatialRayleighDestination,
        cx: &Cx<'_>,
    ) -> Result<DynamicReciprocalSpatialSolid<'s, 'a, 'm, 'c, 'g, 'mesh>> {
        check(cx)?;
        if loads.nodal_forces_n() != self.nodal_forces_n {
            return Err(invalid("initial spatial motion load must match the admitted equilibrium"));
        }
        let problem = self.state.problem(&self.fixed_dofs);
        let run = self.state.begin_transient(
            &problem, loads, self.equilibrium.solution.displacement.displacement_m(),
            velocity_m_s, controls, cx,
        )?;
        let response = run.recover_response(
            self.source.mesh_frame, self.source.maximum_absolute_mandel_strain, cx,
        ).map_err(SpatialSolidEnthalpyError::from)?;
        let motion = run.into_checkpoint(cx).map_err(SpatialSolidEnthalpyError::from)?;
        check(cx)?;
        Ok(DynamicReciprocalSpatialSolid {
            source: self.source, state: self.state, laws: self.laws,
            specific_h: self.specific_h, fixed_dofs: self.fixed_dofs,
            minimum_nodal_capacity_j_k: self.minimum_nodal_capacity_j_k,
            motion, response,
        })
    }
}

impl<'s, 'a: 's, 'm: 'a, 'c: 'a, 'g, 'mesh>
    DynamicReciprocalSpatialSolid<'s, 'a, 'm, 'c, 'g, 'mesh>
{
    pub fn state(&self) -> &SpatialSolidState<'s> { &self.state }
    pub fn response(&self) -> &TransientResponse3 { &self.response }
    pub fn specific_enthalpy_j_kg(&self) -> &[f64] { &self.specific_h }
    pub fn time_s(&self) -> f64 { self.motion.state().time_s }
    pub fn step(&self) -> usize { self.motion.state().step }
    pub fn kinetic_energy_j(&self) -> f64 { self.motion.state().kinetic_energy_j }
    /// Original conservative fixed-strain admission, not a dynamic effective Cp.
    pub fn minimum_nodal_heat_capacities_j_k(&self) -> &[f64] {
        &self.minimum_nodal_capacity_j_k
    }

    /// Actual accepted displacement, velocity and acceleration in native order.
    pub fn nodal_motion(&self, cx: &Cx<'_>) -> Result<[Vec<[f64; 3]>; 3]> {
        let problem = self.state.problem(&self.fixed_dofs);
        let run = self.motion.bind(&problem, cx).map_err(SpatialSolidEnthalpyError::from)?;
        Ok(run.nodal_motion(cx).map_err(SpatialSolidEnthalpyError::from)?)
    }

    /// Advance the original nodal caloric and native inertial histories together.
    /// No static solve replaces motion. Every heat trial starts from exactly the
    /// accepted h field and mechanical checkpoint, retaining the native step cap.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: EnthalpyStepConfig,
        config: ThermoelasticSolve,
        cx: &Cx<'_>,
    ) -> Result<DynamicSpatialThermoelasticCandidate<'s, 'a, 'm, 'c, 'g, 'mesh>> {
        self.advance_with_radiation(
            dt_s, boundary, interfaces, volumetric_source, loads,
            thermal_config, config, EnthalpyRadiation::None, cx,
        )
    }

    /// Compose native pointwise radiation with the same reciprocal trajectory.
    /// The final heat residual is freshly evaluated with actual material heat
    /// and the original radiation owner. Original row, global heat, native
    /// mechanics and combined-energy gates all remain independent.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_radiation(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        interfaces: Option<&ThermalInterfaces>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: EnthalpyStepConfig,
        config: ThermoelasticSolve,
        radiation: EnthalpyRadiation<'_>,
        cx: &Cx<'_>,
    ) -> Result<DynamicSpatialThermoelasticCandidate<'s, 'a, 'm, 'c, 'g, 'mesh>> {
        self.advance_heat_boundary(
            dt_s, boundary, HeatBoundary::Fixed(interfaces, thermal_config),
            volumetric_source, loads, config, radiation, cx,
        )
    }

    /// Live original contact resistance, spatial heat and real inertial motion.
    /// Every native thermal trial queries both independent contact traces;
    /// final source receipts and fluxes belong to the actual accepted endpoint.
    /// Original contact geometry and all undeclared nonthermal coordinates stay
    /// fixed: this does not infer pressure, gap, friction or contact-area motion.
    ///
    /// Contact transfer stays internal. Native kinetic storage, externally
    /// exported Rayleigh loss and pointwise radiation keep their separate
    /// accounts; a refused contact query publishes neither heat nor motion.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_temperature_contacts(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        contacts: &TemperatureContacts<'_, '_>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        contact_config: ContactEnthalpyNewtonConfig,
        config: ThermoelasticSolve,
        radiation: EnthalpyRadiation<'_>,
        cx: &Cx<'_>,
    ) -> Result<DynamicSpatialThermoelasticCandidate<'s, 'a, 'm, 'c, 'g, 'mesh>> {
        self.advance_heat_boundary(
            dt_s, boundary, HeatBoundary::Temperature(contacts, contact_config),
            volumetric_source, loads, config, radiation, cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_heat_boundary(
        &self,
        dt_s: f64,
        boundary: &ThermalBoundary,
        heat_boundary: HeatBoundary<'_, '_, '_>,
        volumetric_source: &ScalarField,
        loads: &TetLoadCase3<'_, '_>,
        config: ThermoelasticSolve,
        radiation: EnthalpyRadiation<'_>,
        cx: &Cx<'_>,
    ) -> Result<DynamicSpatialThermoelasticCandidate<'s, 'a, 'm, 'c, 'g, 'mesh>> {
        check(cx)?;
        let thermal_config = heat_boundary.step_config();
        admit_controls(config)?;
        let end_s = self.time_s() + dt_s;
        if !dt_s.is_finite() || dt_s <= 0.0 || !end_s.is_finite() || end_s <= self.time_s()
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(invalid("representable positive physical interval and compatible native heat tolerance required"));
        }
        let source = self.source;
        let accepted_problem = self.state.problem(&self.fixed_dofs);
        let accepted = self.motion.bind(&accepted_problem, cx)
            .map_err(SpatialSolidEnthalpyError::from)?;
        let mut applied = vec![Dd::ZERO; self.specific_h.len()];
        let mut defects = [f64::INFINITY; 3];
        for iteration in 1..=config.maximum_iterations {
            check(cx)?;
            let mut thermal = heat_boundary.advance(
                source, boundary, volumetric_source, &self.specific_h, dt_s, &applied, radiation, cx,
            )?;
            let state = source.resolve_at(&thermal.specific_enthalpy_j_kg, end_s, cx)?;
            let problem = state.problem(&self.fixed_dofs);
            // Preserve the supplied dt exactly: subtracting rounded absolute
            // clocks can otherwise give mechanics a different heat interval.
            let (run, native) = state.advance_transient_interval(
                &self.state, &accepted, &problem, loads, dt_s, cx,
            )?;
            let response = run.recover_response(
                source.mesh_frame, source.maximum_absolute_mandel_strain, cx,
            ).map_err(SpatialSolidEnthalpyError::from)?;
            let ElementExchange { elements, nodal_heat, work, internal, heat } = element_exchange(
                source, &self.laws, &self.state, &self.response.stress,
                &state, &response.stress, config, cx,
            )?;
            let kinetic = native.native.retained_kinetic_energy_change_j();
            let external = native.native.retained_external_work_j();
            let damping = native.native.retained_damping_loss_j();
            let mechanical_residual = kinetic + work + damping - external;
            let mechanical_residual_j = finite(mechanical_residual)?;
            if !within(mechanical_residual, config.energy_tolerance_j)? {
                return Err(invalid("native inertia, actual stress work, exported damping and applied work disagree"));
            }
            let endpoint = heat_boundary.evaluate(
                source, boundary, volumetric_source, &self.specific_h,
                &thermal.specific_enthalpy_j_kg, dt_s, &nodal_heat, radiation, cx,
            )?;
            let external_heat = endpoint.energy.volumetric_source_j
                + endpoint.energy.dirichlet_in_j - endpoint.energy.neumann_out_j
                - endpoint.energy.robin_out_j - endpoint.energy.radiation_out_j;
            let total = endpoint.energy.stored_energy_change_j + internal + kinetic
                + damping - external - external_heat;
            let mut heat_defect = Dd::ZERO;
            for (&actual, &used) in nodal_heat.iter().zip(&applied) {
                heat_defect = heat_defect + (actual - used).abs();
            }
            let heat_defect_j = finite(heat_defect)?;
            let total_energy_residual_j = finite(total)?;
            defects = [heat_defect_j, endpoint.maximum_free_residual_j, total_energy_residual_j];
            if within(heat_defect, config.heat_tolerance_j)?
                && heat_boundary.admits_contact_residual(endpoint.residual_norm_j)
                && endpoint.maximum_free_residual_j <= thermal_config.energy_tolerance_j
                && within(endpoint.energy.residual_j, thermal_config.energy_tolerance_j)?
                && within(total, config.energy_tolerance_j)?
            {
                let motion = run.into_checkpoint(cx).map_err(SpatialSolidEnthalpyError::from)?;
                thermal.endpoint = endpoint;
                let report = DynamicSpatialThermoelasticReport {
                    iterations: iteration, native, elements, nodal_reversible_heat_j: nodal_heat,
                    kinetic_energy_change_j: finite(kinetic)?, stress_work_j: finite(work)?,
                    external_work_j: finite(external)?, damping_export_j: finite(damping)?,
                    mechanical_residual_j, internal_energy_change_j: finite(internal)?,
                    reversible_heat_j: finite(heat)?, external_thermal_heat_j: finite(external_heat)?,
                    heat_defect_j, total_energy_residual_j,
                };
                let solid = Self {
                    source, state, laws: Arc::clone(&self.laws),
                    specific_h: thermal.specific_enthalpy_j_kg.clone(),
                    fixed_dofs: self.fixed_dofs.clone(),
                    minimum_nodal_capacity_j_k: Arc::clone(&self.minimum_nodal_capacity_j_k),
                    motion, response,
                };
                check(cx)?;
                return Ok(DynamicSpatialThermoelasticCandidate { solid, thermal, report });
            }
            for (used, actual) in applied.iter_mut().zip(nodal_heat) {
                *used = *used + dd(config.heat_relaxation) * (actual - *used);
                finite(*used)?;
            }
        }
        Err(SpatialThermoelasticError::NotConverged {
            iterations: config.maximum_iterations,
            heat_defect_j: defects[0], maximum_nodal_residual_j: defects[1],
            combined_residual_j: defects[2],
        })
    }
}
