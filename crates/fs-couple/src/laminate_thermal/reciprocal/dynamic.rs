//! Reciprocal finite heat on the actual enriched shell trajectory.
//!
//! The two internal director coordinates per triangle remain independent
//! mechanical unknowns. In the native basis y = T q + B r, T is only a change
//! of coordinates: r is evolved with its original mass, cross inertia and
//! stiffness. It is never reset to its static thermal offset after a step.
//! This avoids a prescribed internal lift with an unaccounted actuator port.
//!
//! Implicit midpoint dynamics and the original spatial heat solver propose
//! from the same accepted endpoint on every coupling iteration. Each original
//! heat cell receives W_stress - delta(Psi + T S), including stress work on the
//! actual internal motion. Kinetic energy, native corrected shear storage and
//! numerical drilling storage remain separate from thermal energy. Changing
//! elastic sources retain the INITIAL full coordinate basis and physical mass,
//! while the material operator, thermal loading and nonlocal shear entropy
//! follow each actual temperature. The numerical drilling stiffness stays
//! fixed at launch; it is not a material heating or inertia mechanism.
//! The parent's full-temperature-box capacity and strain admission apply.
//! Explicit source-bound in-plane Kelvin viscosity may add physical damping
//! and return its independently evaluated heat to each original thermal cell.
//! It is not finite-rotation dynamics, moving thermal geometry, hereditary
//! material memory or exact discrete isentropy.

mod viscosity;
pub use viscosity::LaminateViscosityControls;

use super::*;
use fs_plate::shell::mitc::thermal::transient::{
    EnrichedMitcShell, EnrichedThermalResponse, MitcTransientReport, MitcTransientState,
    PreparedMitcMidpoint,
};

/// Accepted thermal material, compatible physical strain and inertial motion.
///
/// Observation reads this state directly; it never re-solves equilibrium or
/// substitutes the statically condensed internal-director displacement.
pub struct DynamicLaminateState {
    thermomechanical: ReciprocalLaminateState,
    motion: MitcTransientState,
    bubble_material_force: Vec<[f64; 2]>,
}

impl DynamicLaminateState {
    pub fn time_s(&self) -> f64 {
        self.motion.time_s()
    }
    pub fn temperatures_k(&self) -> &[f64] {
        self.thermomechanical.shell.temperatures_k()
    }
    pub fn material(&self) -> &ResolvedLaminateThermal {
        self.thermomechanical.shell.material()
    }
    /// Actual native material response at the inertial endpoint.
    pub fn response(&self) -> &ThermalShellResponse {
        self.thermomechanical.shell.response()
    }
    /// Native polynomials with the evolved internal director coordinates.
    pub fn strain_fields(&self) -> &[RealElementStrain] {
        self.thermomechanical.shell.strain_fields()
    }
    pub fn cell_potentials(&self) -> &[LaminateCellPotential] {
        self.thermomechanical.cell_potentials()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.motion.kinetic_energy_j()
    }
    pub fn transverse_shear_energy_j(&self) -> f64 {
        self.thermomechanical.transverse_shear_energy_j()
    }
    pub fn transverse_shear_internal_energy_j(&self) -> f64 {
        self.thermomechanical.transverse_shear_internal_energy_j()
    }
    pub fn transverse_shear_entropy_j_k(&self) -> &[f64] {
        self.thermomechanical.transverse_shear_entropy_j_k()
    }
    pub fn drilling_gauge_energy_j(&self) -> f64 {
        self.thermomechanical.drilling_gauge_energy_j()
    }
    pub fn nodal_velocity_6n(&self) -> &[f64] {
        self.motion.nodal_velocity_6n()
    }
    /// Native generalized coordinates and velocities, including every retained
    /// internal mode; reference mass and stationary supports are unchanged.
    pub const fn motion(&self) -> &MitcTransientState {
        &self.motion
    }
    pub fn external_load_6n(&self) -> &[f64] {
        self.thermomechanical.external_load_6n()
    }
    pub fn reduced_elastic_strain_upper_bound(&self) -> f64 {
        self.thermomechanical
            .shell
            .reduced_elastic_strain_upper_bound()
    }
}

/// Combined account over one accepted physical interval.
#[derive(Clone, Debug, PartialEq)]
pub struct DynamicLaminateReport {
    /// Original per-cell heat and source-potential account. Its mechanical and
    /// combined residuals INCLUDE the kinetic-energy change reported below.
    pub thermoelastic: LaminateThermoelasticReport,
    pub kinetic_energy_change_j: f64,
    /// Elastic stress work on independent residual director coordinates.
    pub bubble_material_work_j: f64,
    /// Elastic plus viscous work on nodal and internal coordinates, excluding
    /// numerical drilling. Viscous heat is also reported separately per cell.
    pub generalized_material_work_j: f64,
    /// Native linear mechanical identity; its load includes thermal forcing.
    /// It is an independent check, not an additional heat or energy credit.
    pub native: MitcTransientReport,
}

pub struct DynamicLaminateStep<'m> {
    pub thermal: DepositionSimulation<'m>,
    pub report: DynamicLaminateReport,
}

/// Finite heat coupled to retained native displacements and velocities.
pub struct DynamicLayeredThermalShell<'s, 'b, 'a> {
    source: &'s LayeredThermalShell<'b, 'a>,
    laws: Vec<AffineLaminaThermoelastic>,
    native: EnrichedMitcShell,
    fixed_elasticity: bool,
    midpoint: Option<PreparedMitcMidpoint>,
    viscosity: Option<viscosity::LayerViscosity>,
    state: DynamicLaminateState,
    quadrature: Vec<([f64; 3], f64)>,
}

impl<'s, 'b, 'a> ReciprocalLayeredThermalShell<'s, 'b, 'a> {
    /// Launch inertia from this admitted thermal equilibrium.
    ///
    /// Nodal velocities are physical world velocities. The two velocities per
    /// facet are residual internal director rates in y=Tq+Br; zero means the
    /// initial internal velocity follows T times the nodal velocity. They are
    /// explicit initial conditions, never an inferred temperature rate. The
    /// original consistent mass is required; a lumped model without internal
    /// inertia cannot establish this full dynamic state. Free massless drilling
    /// directions must satisfy the native homogeneous gauge in both initial
    /// displacement and velocity; director-axis gauge efforts are unsupported.
    /// The launch basis, full physical reference mass and numerical drilling
    /// stiffness stay fixed. Later source stiffness changes act on those SAME
    /// coordinates, retaining every internal motion and cross-inertia term.
    pub fn begin_transient(
        self,
        nodal_velocity_6n: &[f64],
        residual_bubble_velocity_rad_s: &[[f64; 2]],
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DynamicLayeredThermalShell<'s, 'b, 'a>, LaminateThermoelasticError> {
        check(&mut cancelled)?;
        if !energy_tolerance_j.is_finite() || energy_tolerance_j <= 0.0 {
            return Err(invalid("positive finite launch energy tolerance required"));
        }
        let fixed_elasticity = self
            .laws
            .iter()
            .all(AffineLaminaThermoelastic::is_temperature_independent_elasticity);
        let native =
            EnrichedMitcShell::try_new(Arc::clone(&self.state.shell.model), &mut cancelled)
                .map_err(LaminateThermoelasticError::Native)?;
        let loading = native
            .model()
            .prepare_laminate_thermal(&self.state.shell.thermal_fields, &mut cancelled)
            .map_err(LaminateThermoelasticError::Native)?;
        let seed = native
            .equilibrium_state(
                self.time_s(),
                &loading,
                &self.state.shell.response.displacement_6n,
                &mut cancelled,
            )
            .map_err(LaminateThermoelasticError::Native)?;
        let motion = native
            .state(
                self.time_s(),
                &self.state.shell.response.displacement_6n,
                nodal_velocity_6n,
                seed.bubble_coordinates(),
                residual_bubble_velocity_rad_s,
                &mut cancelled,
            )
            .map_err(LaminateThermoelasticError::Native)?;
        let recovered = native
            .recover(&loading, &motion, &mut cancelled)
            .map_err(LaminateThermoelasticError::Native)?;
        gate(
            recovered.response.strain_energy_j - self.state.shell.response.strain_energy_j,
            energy_tolerance_j,
            "dynamic launch must retain the original equilibrium material storage",
        )?;
        gate(
            recovered.drilling_gauge_energy_j - self.state.gauge_energy_j,
            energy_tolerance_j,
            "dynamic launch must preserve the original numerical gauge",
        )?;
        let ReciprocalLaminateState { shell, loads, .. } = self.state;
        let prepared = PreparedLaminateSnapshot {
            temperatures_k: shell.temperatures_k,
            material: shell.material,
            model: shell.model,
            fields: shell.thermal_fields,
        };
        let thermal_load = shell.thermal_load_6n;
        let (shell, bubble_material_force, native_gauge) = admit_recovered(
            self.source,
            motion.time_s(),
            prepared,
            recovered,
            thermal_load,
            &mut cancelled,
        )?;
        let proposal = assemble_state(
            self.source,
            &self.laws,
            &self.quadrature,
            shell,
            &loads,
            Some(native_gauge),
            None,
            energy_tolerance_j,
            &mut cancelled,
        )?;
        check(&mut cancelled)?;
        Ok(DynamicLayeredThermalShell {
            source: self.source,
            laws: self.laws,
            native,
            fixed_elasticity,
            midpoint: None,
            viscosity: None,
            state: DynamicLaminateState {
                thermomechanical: proposal.state,
                motion,
                bubble_material_force,
            },
            quadrature: self.quadrature,
        })
    }
}

impl<'s, 'b, 'a> DynamicLayeredThermalShell<'s, 'b, 'a> {
    pub const fn state(&self) -> &DynamicLaminateState {
        &self.state
    }
    pub fn time_s(&self) -> f64 {
        self.state.time_s()
    }

    /// Propose and accept one reciprocal heat/motion interval.
    ///
    /// Constant elastic sources reuse one effective sparse factor for the exact
    /// requested dt. Varying sources prepare the current physical stiffness in
    /// the original reference basis and use both accepted and proposed endpoint
    /// forces. No trial advances the accepted motion
    /// or caloric history. Range exits, source mismatches, nonconvergence,
    /// work/energy refusals and cancellation leave both physical histories
    /// untouched; the caller's thermal workspace may retain scratch values.
    /// Supports are stationary and loads are endpoint world nodal loads.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'m>(
        &mut self,
        initial: &DepositionSimulation<'m>,
        target_s: f64,
        external_load_6n: &[f64],
        conduction: LayerConductionConfig,
        heat_controls: LayeredThermalSolveControls,
        controls: LaminateThermoelasticControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DynamicLaminateStep<'m>, LaminateThermoelasticError> {
        check(&mut cancelled)?;
        admit_controls(controls)?;
        let dt = target_s - self.time_s();
        if !target_s.is_finite() || !dt.is_finite() || dt <= 0.0 {
            return Err(invalid(
                "positive finite inertial thermal interval required",
            ));
        }
        let initial_t = self
            .source
            .thermal
            .cell_temperatures(initial, || thermal_poll(&mut cancelled))
            .map_err(LaminateThermoelasticError::Thermal)?;
        if initial.time_s().to_bits() != self.time_s().to_bits()
            || initial_t.len() != self.state.temperatures_k().len()
            || initial_t
                .iter()
                .zip(self.state.temperatures_k())
                .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err(invalid(
                "accepted inertial shell and original thermal history must remain exactly paired",
            ));
        }
        let old = &self.state.thermomechanical;
        let old_thermal = self
            .native
            .model()
            .prepare_laminate_thermal(&old.shell.thermal_fields, &mut cancelled)
            .map_err(LaminateThermoelasticError::Native)?;
        let old_load = self
            .native
            .generalized_load(&old_thermal, &old.loads, &mut cancelled)
            .map_err(LaminateThermoelasticError::Native)?;
        // Validate nodal load length/finiteness before the expensive factor.
        self.native
            .model()
            .project_load(external_load_6n)
            .map_err(LaminateThermoelasticError::Native)?;
        if self.fixed_elasticity
            && self
                .midpoint
                .as_ref()
                .is_none_or(|factor| factor.dt_s().to_bits() != dt.to_bits())
        {
            self.midpoint = Some(
                self.native
                    .prepare_midpoint(dt, &mut cancelled)
                    .map_err(LaminateThermoelasticError::Native)?,
            );
        }
        let mut applied = vec![0.0; self.laws.len()];
        let mut last = (f64::INFINITY, f64::INFINITY);
        for iteration in 0..controls.maximum_iterations {
            check(&mut cancelled)?;
            let step =
                LayeredThermalStep::new(self.source.thermal, initial, target_s, conduction, || {
                    thermal_poll(&mut cancelled)
                })
                .map_err(LaminateThermoelasticError::Thermal)?
                .with_cell_energy_j(&applied, || thermal_poll(&mut cancelled))
                .map_err(LaminateThermoelasticError::Thermal)?;
            let (thermal, _, heat) = step
                .solve(heat_controls, workspace, || thermal_poll(&mut cancelled))
                .map_err(LaminateThermoelasticError::HeatSolve)?;
            let prepared = self
                .source
                .prepare_snapshot(
                    &thermal,
                    self.fixed_elasticity.then(|| self.native.shared_model()),
                    &mut cancelled,
                )
                .map_err(LaminateThermoelasticError::Shell)?;
            // The basis and mass belong to launch. Only the physical material
            // coefficients change; the two internal directors remain actual
            // unknowns in the same coordinates through EVERY thermal trial.
            let trial_viscosity = self
                .viscosity
                .as_ref()
                .map(|loss| {
                    loss.prepare_midpoint(
                        self.source,
                        self.state.temperatures_k(),
                        &prepared.temperatures_k,
                        &mut cancelled,
                    )
                })
                .transpose()?;
            let same_viscosity = match (&self.viscosity, &trial_viscosity) {
                (Some(loss), Some(prepared)) => prepared.same_operator(&loss.accepted),
                _ => true,
            };
            // One frozen midpoint law applies to BOTH endpoint mechanical
            // owners. Rebinding changes no state, mass or elastic operator.
            let interval_native = match (&self.viscosity, &trial_viscosity) {
                (Some(loss), Some(prepared)) => {
                    loss.interval_native(&self.native, prepared, &mut cancelled)?
                }
                _ => self.native.clone(),
            };
            let trial_native = if self.fixed_elasticity {
                interval_native.clone()
            } else {
                interval_native
                    .with_material_model(Arc::clone(&prepared.model), &mut cancelled)
                    .map_err(LaminateThermoelasticError::Native)?
            };
            let trial_factor = if self.fixed_elasticity && same_viscosity {
                None
            } else {
                Some(
                    trial_native
                        .prepare_midpoint(dt, &mut cancelled)
                        .map_err(LaminateThermoelasticError::Native)?,
                )
            };
            let mechanical = trial_factor
                .as_ref()
                .or(self.midpoint.as_ref())
                .ok_or_else(|| invalid("missing admitted mechanical factor"))?;
            let loading = prepared
                .model
                .prepare_laminate_thermal(&prepared.fields, &mut cancelled)
                .map_err(LaminateThermoelasticError::Native)?;
            let new_load = trial_native
                .generalized_load(&loading, external_load_6n, &mut cancelled)
                .map_err(LaminateThermoelasticError::Native)?;
            let candidate = mechanical
                .advance_from(
                    &interval_native,
                    &self.state.motion,
                    &old_load,
                    &new_load,
                    &mut cancelled,
                )
                .map_err(LaminateThermoelasticError::Native)?;
            if candidate.state.time_s().to_bits() != thermal.time_s().to_bits() {
                return Err(invalid("thermal and mechanical endpoint times differ"));
            }
            gate(
                candidate.report.energy_residual_j,
                controls.energy_tolerance_j,
                "native midpoint kinetic, quadratic and parameter work do not close",
            )?;
            let recovered = trial_native
                .recover(&loading, &candidate.state, &mut cancelled)
                .map_err(LaminateThermoelasticError::Native)?;
            let thermal_load = loading.load_6n().to_vec();
            let (shell, bubble_force, native_gauge) = admit_recovered(
                self.source,
                thermal.time_s(),
                prepared,
                recovered,
                thermal_load,
                &mut cancelled,
            )?;
            let proposal = assemble_state(
                self.source,
                &self.laws,
                &self.quadrature,
                shell,
                external_load_6n,
                Some(native_gauge),
                Some(old),
                controls.energy_tolerance_j,
                &mut cancelled,
            )?;
            let mut bubble_work = 0.0;
            for (element, (a, b)) in self
                .state
                .motion
                .bubble_coordinates()
                .iter()
                .zip(candidate.state.bubble_coordinates())
                .enumerate()
            {
                check(&mut cancelled)?;
                for axis in 0..2 {
                    bubble_work = finite(
                        bubble_work
                            + 0.5
                                * (self.state.bubble_material_force[element][axis]
                                    + bubble_force[element][axis])
                                * (b[axis] - a[axis]),
                    )?;
                }
            }
            // The original midpoint owner evaluates the accepted endpoint
            // increment on the unchanged full reference mass. Re-subtracting
            // rounded totals here would discard small physical stress work
            // beside an unchanged large kinetic inventory.
            let kinetic_change = finite(candidate.report.kinetic_energy_change_j)?;
            let cell_dissipated = if let Some(loss) = &trial_viscosity {
                let rates = trial_native
                    .midpoint_strain_rate_fields(
                        &self.state.motion,
                        &candidate.state,
                        &mut cancelled,
                    )
                    .map_err(LaminateThermoelasticError::Native)?;
                loss.cell_heat(
                    self.source,
                    &rates,
                    &self.quadrature,
                    dt,
                    candidate.report.dissipated_energy_j,
                    &mut cancelled,
                )?
            } else {
                Vec::new()
            };
            let thermoelastic = account(
                old,
                &proposal,
                &applied,
                heat,
                iteration + 1,
                controls.energy_tolerance_j,
                Some(InertialWork {
                    kinetic_energy_change_j: kinetic_change,
                    bubble_material_work_j: bubble_work,
                    cell_dissipated_energy_j: &cell_dissipated,
                }),
                &mut cancelled,
            )?;
            last = (
                thermoelastic.maximum_heat_defect_j,
                thermoelastic.combined_energy_residual_j,
            );
            if thermoelastic.maximum_heat_defect_j <= controls.cell_heat_tolerance_j
                && thermoelastic.total_absolute_heat_defect_j <= controls.energy_tolerance_j
                && thermoelastic.combined_energy_residual_j.abs() <= controls.energy_tolerance_j
            {
                let generalized_material_work_j = finite(
                    proposal.mechanical_work.iter().sum::<f64>()
                        + proposal.shear_work
                        + thermoelastic.dissipated_energy_j,
                )?;
                let report = DynamicLaminateReport {
                    thermoelastic,
                    kinetic_energy_change_j: kinetic_change,
                    bubble_material_work_j: bubble_work,
                    generalized_material_work_j,
                    native: candidate.report,
                };
                check(&mut cancelled)?;
                if self.fixed_elasticity
                    && let Some(factor) = trial_factor
                {
                    self.midpoint = Some(factor);
                }
                if let (Some(loss), Some(prepared)) = (&mut self.viscosity, trial_viscosity) {
                    loss.accepted = prepared;
                }
                self.native = trial_native;
                self.state = DynamicLaminateState {
                    thermomechanical: proposal.state,
                    motion: candidate.state,
                    bubble_material_force: bubble_force,
                };
                return Ok(DynamicLaminateStep { thermal, report });
            }
            for (input, cell) in applied.iter_mut().zip(&thermoelastic.cells) {
                check(&mut cancelled)?;
                *input = finite(*input + controls.heat_relaxation * cell.heat_defect_j)?;
            }
        }
        Err(LaminateThermoelasticError::NotConverged {
            iterations: controls.maximum_iterations,
            maximum_heat_defect_j: last.0,
            combined_residual_j: last.1,
        })
    }
}

fn admit_recovered(
    source: &LayeredThermalShell<'_, '_>,
    time_s: f64,
    prepared: PreparedLaminateSnapshot,
    recovered: EnrichedThermalResponse,
    thermal_load: Vec<f64>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(LaminateShellSnapshot, Vec<[f64; 2]>, NumericalGaugeState), LaminateThermoelasticError>
{
    envelope(
        "complete director motion",
        recovered.maximum_director_motion_rad,
        source.limits.maximum_rotation_rad,
    )
    .map_err(LaminateThermoelasticError::Shell)?;
    let nodal = prepared.model.pencil().free;
    if recovered.numerical_gauge_force.len() != nodal + 2 * recovered.bubble_material_force.len()
        || recovered.numerical_gauge_force[nodal..]
            .iter()
            .any(|force| *force != 0.0)
    {
        return Err(invalid(
            "reference numerical drilling gauge cannot act on physical interior directors",
        ));
    }
    let gauge = NumericalGaugeState {
        energy_j: recovered.drilling_gauge_energy_j,
        force: recovered.numerical_gauge_force[..nodal].to_vec(),
    };
    let shell = source
        .admit_snapshot(
            time_s,
            prepared,
            recovered.response,
            recovered.strain_fields,
            thermal_load,
            cancelled,
        )
        .map_err(LaminateThermoelasticError::Shell)?;
    Ok((shell, recovered.bubble_material_force, gauge))
}
