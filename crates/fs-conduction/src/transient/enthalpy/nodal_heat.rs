//! Signed interval heat on the original spatial enthalpy residual.
//!
//! Each entry is already an assembled nodal heat amount [J], not a volumetric
//! density or a power. It is added once to the original P1 source load. The
//! source is held fixed during one native Newton solve; an endpoint-dependent
//! physical owner must iterate it and re-evaluate the final endpoint.
//!
//! This reuses the original caloric charts, reference inventory, conductivity,
//! contact, boundary assembly, Newton tangent and cancellation. Prescribed
//! temperature nodes retain their supplied heat in the physical reaction.
//! Local and global energy checks use retained arithmetic independently of the
//! Newton stopping criterion. No residual becomes a physical heat source.
//!
//! This entry has the original backward-Euler temporal order and stationary
//! reference geometry. The radiation entry retains the original pointwise
//! ambient/enclosure owners and their material response. Neither entry adds an
//! entropy-flow model or a derivative of the supplied heat with respect to state.

/// Live source-backed contact and signed nodal heat in the same Newton solve.
pub mod contact;
/// Native pointwise radiation and signed nodal heat in the same Newton solve.
pub mod radiation;

use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_solver::{NewtonReport, norm2};

use super::interference::{checked_dd, product};
use super::{
    ASSEMBLY_TILE, DofMap, EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig,
    StepContext, energy_balance, finite, poll,
};
use crate::{
    ConductionProblem, InterfaceFlux, PointwiseEnclosureReport,
    PointwiseRadiationReport, RobinFlux, ThermalInterfaces,
};
use super::pointwise::PointwiseBoundary;

/// Independently integrated physical interval energies; positive input heats
/// the original inventory. Boundary heat is positive outward except Dirichlet.
#[derive(Clone, Debug)]
pub struct NodalHeatEnergy {
    pub nodal_heat_j: Dd,
    pub volumetric_source_j: Dd,
    pub stored_energy_change_j: Dd,
    pub neumann_out_j: Dd,
    pub robin_out_j: Dd,
    /// Original pointwise radiation, positive outward; zero without radiation.
    pub radiation_out_j: Dd,
    pub dirichlet_in_j: Dd,
    pub residual_j: Dd,
}

/// Measurements at a supplied endpoint. Evaluation alone does not admit it.
#[derive(Clone, Debug)]
pub struct NodalHeatEndpoint {
    pub temperature: Vec<f64>,
    /// Full physical joule residual in vertex order, zero on prescribed rows.
    pub residual_j: Vec<f64>,
    pub residual_norm_j: f64,
    pub maximum_free_residual_j: f64,
    pub free_dofs: usize,
    pub energy: NodalHeatEnergy,
    pub robin_fluxes: Vec<RobinFlux>,
    pub contact_fluxes: Vec<InterfaceFlux>,
    /// Actual endpoint source queries for a live temperature/phase contact set.
    /// Geometry and all nonthermal source coordinates remain held by that owner.
    pub temperature_contacts: Option<crate::interface::temperature::TemperatureContactEvaluation>,
    /// Original ambient radiation report; convection remains separate.
    pub pointwise_radiation: Option<PointwiseRadiationReport>,
    /// Original closed/open enclosure report with reflected exchange.
    pub pointwise_enclosure_radiation: Option<PointwiseEnclosureReport>,
}

/// An admitted step from the original Newton solve and original inventory.
#[derive(Clone, Debug)]
pub struct NodalHeatStep {
    pub specific_enthalpy_j_kg: Vec<f64>,
    pub liquid_mass_fraction: Vec<f64>,
    pub endpoint: NodalHeatEndpoint,
    pub newton: NewtonReport,
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Add signed assembled interval joules to the original volumetric input.
    /// Each source entry belongs to exactly one original vertex, including
    /// prescribed vertices. Neither input history nor the supplied heat changes.
    ///
    /// The original Newton tangent is exact for this held source. Both every
    /// free-row residual and the independently integrated global balance must
    /// satisfy the unchanged step energy tolerance before returning a candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_nodal_heat(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[Dd],
    ) -> Result<NodalHeatStep, EnthalpyError> {
        self.advance_nodal_heat(
            cx, problem, interfaces, old_specific_h, dt_s, config, nodal_heat_j, None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_nodal_heat(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[Dd],
        pointwise_radiation: Option<PointwiseBoundary<'_>>,
    ) -> Result<NodalHeatStep, EnthalpyError> {
        self.admit_step(cx, problem, old_specific_h, dt_s, config)?;
        let heat = admit_heat(cx, self.masses.len(), nodal_heat_j)?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_specific_h,
            dt: dt_s,
            pointwise_radiation,
        };
        let (h, newton) = context.solve_newton(config, Some(&heat))?;
        let endpoint = evaluate(&context, &h, &heat)?;
        self.finish_nodal_heat(cx, h, newton, endpoint, config)
    }

    fn finish_nodal_heat(
        &self,
        cx: &Cx<'_>,
        h: Vec<f64>,
        newton: NewtonReport,
        endpoint: NodalHeatEndpoint,
        config: EnthalpyStepConfig,
    ) -> Result<NodalHeatStep, EnthalpyError> {
        let tolerance = Dd::from_f64(config.energy_tolerance_j);
        if endpoint.maximum_free_residual_j > config.energy_tolerance_j {
            return Err(EnthalpyError::EndpointResidual {
                residual_j: endpoint.maximum_free_residual_j,
                tolerance_j: config.energy_tolerance_j,
            });
        }
        if positive(absolute(endpoint.energy.residual_j) - tolerance) {
            return Err(EnthalpyError::EnergyBalance {
                residual_j: endpoint.energy.residual_j.to_f64(),
                tolerance_j: config.energy_tolerance_j,
            });
        }
        let mut liquid_mass_fraction = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            liquid_mass_fraction.push(
                self.curve_for_vertex(vertex)
                    .values_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?
                    .liquid_mass_fraction.clamp(0.0, 1.0),
            );
        }
        poll(cx, h.len())?;
        Ok(NodalHeatStep {
            specific_enthalpy_j_kg: h,
            liquid_mass_fraction,
            endpoint,
            newton,
        })
    }

    /// Measure the original physical endpoint with a freshly supplied heat
    /// vector, without taking a Newton step or changing accepted history.
    /// This permits an outer constitutive coupling to reject equal-and-opposite
    /// local defects even when its total heat happens to balance.
    ///
    /// Prescribed enthalpies, source/domain admission and finite arithmetic
    /// remain mandatory. Residual tolerances belong to the calling owner.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_nodal_heat(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[Dd],
    ) -> Result<NodalHeatEndpoint, EnthalpyError> {
        self.evaluate_nodal_heat(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, nodal_heat_j, None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_nodal_heat(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[Dd],
        pointwise_radiation: Option<PointwiseBoundary<'_>>,
    ) -> Result<NodalHeatEndpoint, EnthalpyError> {
        self.admit_step(cx, problem, old_specific_h, dt_s, config)?;
        let heat = admit_heat(cx, self.masses.len(), nodal_heat_j)?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_specific_h,
            dt: dt_s,
            pointwise_radiation,
        };
        evaluate(&context, candidate_specific_h, &heat)
    }
}

fn admit_heat(cx: &Cx<'_>, n: usize, heat: &[Dd]) -> Result<Vec<Dd>, EnthalpyError> {
    if heat.len() != n {
        return Err(EnthalpyError::InvalidInput(
            "one signed interval heat amount is required per original vertex",
        ));
    }
    let mut admitted = Vec::with_capacity(n);
    for (vertex, value) in heat.iter().enumerate() {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(cx, vertex)?;
        }
        admitted.push(checked_dd(Dd::from_pair(value.hi, value.lo))?);
    }
    poll(cx, n)?;
    Ok(admitted)
}

fn positive(value: Dd) -> bool {
    value.hi > 0.0 || (value.hi == 0.0 && value.lo > 0.0)
}
fn absolute(value: Dd) -> Dd {
    if positive(-value) { -value } else { value }
}
fn represented(value: Dd) -> Result<f64, EnthalpyError> {
    let rounded = finite(value.to_f64())?;
    if rounded == 0.0 && (value.hi != 0.0 || value.lo != 0.0) {
        Err(EnthalpyError::NonFiniteArithmetic)
    } else {
        Ok(rounded)
    }
}

fn evaluate(
    context: &StepContext<'_, '_, '_, '_>,
    h: &[f64],
    heat: &[Dd],
) -> Result<NodalHeatEndpoint, EnthalpyError> {
    if h.len() != context.storage.masses.len() {
        return Err(EnthalpyError::InvalidInput(
            "candidate enthalpy must match the original reference inventory",
        ));
    }
    context.check_prescribed_values(h)?;
    // Endpoint measurement needs no tangent, but both original caloric states
    // must still be valid even for a stationary or completely prescribed mesh.
    context.temperatures(context.old)?;
    let temperature = context.temperatures(h)?;
    let scales = context.phase_scales(h)?;
    let system = context.assemble(&temperature, scales.as_ref())?;
    let radiation = context.pointwise_radiation
        .map(|owner| owner.load_and_report(context.cx, &temperature, h)).transpose()?;
    let dt = Dd::from_f64(context.dt);
    let mut stored = Dd::ZERO;
    let mut supplied = Dd::ZERO;
    let mut dirichlet = Dd::ZERO;
    let mut residual_j = Vec::with_capacity(h.len());
    let mut maximum = Dd::ZERO;
    for (vertex, &value) in h.iter().enumerate() {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(context.cx, vertex)?;
        }
        let change = context.storage.storage_change_dd(
            vertex, value, context.old[vertex], temperature[vertex],
        )?;
        stored = checked_dd(stored + change)?;
        supplied = checked_dd(supplied + heat[vertex])?;
        let (transport, load) = system.heat_row_terms(vertex, &temperature);
        let radiative_load = match &radiation {
            Some((load, _)) => Dd::from_f64(finite(load[vertex])?),
            None => Dd::ZERO,
        };
        let transport = product(
            dt, checked_dd(Dd::from_f64(transport) - Dd::from_f64(load) + radiative_load)?,
        )?;
        let residual = checked_dd(change + transport - heat[vertex])?;
        if context.is_prescribed(vertex) {
            dirichlet = checked_dd(dirichlet + residual)?;
            residual_j.push(0.0);
        } else {
            let magnitude = absolute(residual);
            if positive(magnitude - maximum) {
                maximum = magnitude;
            }
            residual_j.push(represented(residual)?);
        }
    }
    let dofs = DofMap::new_allow_fully_prescribed(context.problem.boundary, h.len())?;
    let (boundary, robin_fluxes) = energy_balance(
        context.storage.mesh,
        context.problem.boundary,
        context.problem.source,
        &system,
        &dofs,
        &temperature,
    );
    let volumetric = product(dt, Dd::from_f64(finite(boundary.source_w)?))?;
    let neumann = product(dt, Dd::from_f64(finite(boundary.neumann_out_w)?))?;
    let robin = product(dt, Dd::from_f64(finite(boundary.robin_out_w)?))?;
    let radiation_out = match &radiation {
        Some((_, report)) => product(dt, Dd::from_f64(finite(report.outward_heat_w())?))?,
        None => Dd::ZERO,
    };
    let residual = checked_dd(
        stored - (volumetric + supplied + dirichlet - neumann - robin - radiation_out),
    )?;
    let contact_fluxes = match context.interfaces {
        Some(interfaces) => interfaces.fluxes(&temperature)?,
        None => Vec::new(),
    };
    let residual_norm_j = finite(norm2(&residual_j))?;
    // Round a positive maximum outward so an accepted local gate cannot hide
    // a retained low component just beyond the caller's f64 tolerance.
    let rounded_maximum = represented(maximum)?;
    let maximum_free_residual_j = finite(
        if positive(maximum - Dd::from_f64(rounded_maximum)) {
            rounded_maximum.next_up()
        } else {
            rounded_maximum
        },
    )?;
    let (pointwise_radiation, pointwise_enclosure_radiation) = radiation
        .map_or((None, None), |(_, report)| report.into_reports());
    poll(context.cx, h.len())?;
    Ok(NodalHeatEndpoint {
        temperature,
        residual_j,
        residual_norm_j,
        maximum_free_residual_j,
        free_dofs: dofs.n(),
        energy: NodalHeatEnergy {
            nodal_heat_j: supplied,
            volumetric_source_j: volumetric,
            stored_energy_change_j: stored,
            neumann_out_j: neumann,
            robin_out_j: robin,
            radiation_out_j: radiation_out,
            dirichlet_in_j: dirichlet,
            residual_j: residual,
        },
        robin_fluxes,
        contact_fluxes,
        temperature_contacts: None,
        pointwise_radiation,
        pointwise_enclosure_radiation,
    })
}
