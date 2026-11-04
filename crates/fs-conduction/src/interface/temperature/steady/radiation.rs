//! Ambient patch radiation around the original steady contact/bulk owners.

use fs_exec::Cx;
use fs_solver::norm2;

use crate::assemble::{assemble_operator_scaled_with_interfaces, full_residual};
use crate::radiation::{AmbientRadiationTrial, solve_ambient_radiation_with};
use crate::solve::energy_balance;
use crate::{
    AmbientRadiationConfig, AmbientRadiationPatch, AmbientRadiationPatchReport,
    AmbientRadiationReport, ConductionError, ConductionProblem, DofMap, EnergyBalance,
    RobinFlux, ScalarField, SolveConfig, ThermalBc, ThermalBoundary,
};

use super::super::{TemperatureContactEvaluation, finite, invalid, poll};
use super::{ContactSteadyNewtonSolution, ContactSteadySolution, SteadyContactConfig,
    SteadyContactNewtonConfig, TemperatureContacts};

/// One unchanged field accepted against both constitutive laws.
#[derive(Debug)]
pub struct ContactSteadyNewtonRadiationSolution {
    /// Physical endpoint rows use freshly queried contacts and radiation at
    /// this same field. The numerical history remains the final Newton solve
    /// with its applied radiative secant, not a monolithic radiation Newton.
    pub contact: ContactSteadyNewtonSolution,
    /// Exact radiative secant boundary used by the final numerical response.
    pub applied_boundary: ThermalBoundary,
    /// Exact radiative boundary reassembled at the unchanged accepted field.
    pub physical_boundary: ThermalBoundary,
    /// Original convection only; radiative heat never becomes heat into air.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Sum of original convective heat, W outward.
    pub convective_out_w: f64,
    /// Original applied/physical patch evidence and all Newton/Krylov work.
    pub radiation: AmbientRadiationReport,
    /// Residual callbacks across every attempted contact/bulk response.
    pub residual_evaluations: usize,
    /// Complete contact/bulk tangents across every attempted response.
    pub jacobian_assemblies: usize,
}

/// Bounded contact Picard responses inside the same ambient-radiation owner.
#[derive(Debug)]
pub struct ContactSteadyRadiationSolution {
    /// Last contact solve with freshly audited physical endpoint rows. Its
    /// inner numerical history still describes its original applied operator.
    pub contact: ContactSteadySolution,
    /// Combined radiative boundary used by the last thermal response.
    pub applied_boundary: ThermalBoundary,
    /// Combined boundary freshly evaluated at the unchanged accepted field.
    pub physical_boundary: ThermalBoundary,
    /// Original convective heat only, available to transported-air consumers.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Sum of original convective heat, W outward.
    pub convective_out_w: f64,
    /// Radiation callbacks and cumulative nonlinear/Krylov work.
    pub radiation: AmbientRadiationReport,
    /// All inner spatial solves, including discarded contact/radiation trials.
    /// This differs from the number of ambient-radiation callbacks.
    pub thermal_responses: usize,
}

struct Trial<T> {
    contact: T,
    physical_boundary: Option<ThermalBoundary>,
}

struct Endpoint {
    physical_contacts: TemperatureContactEvaluation,
    free_residual_w: Vec<(usize, f64)>,
    residual_norm_w: f64,
    energy: EnergyBalance,
    robin_fluxes: Vec<RobinFlux>,
    dirichlet_reactions_w: Vec<(usize, f64)>,
    boundary: ThermalBoundary,
}

impl TemperatureContacts<'_, '_> {
    /// Compose source-backed contact and heterogeneous k(T) with the existing
    /// area-mean gray radiation model, including radiation-only zero-flux rows.
    ///
    /// Each ambient-radiation callback uses the original complete bulk/contact
    /// Newton owner and its analytic tangent. Only its warm start changes;
    /// source selection, pressure, geometry and the physical boundary remain
    /// fixed. The original radiation owner retains patch admission, raw
    /// temperature/watt gates and the separation of convection from radiation.
    ///
    /// Once those gates pass, acceptance independently requeries every contact
    /// and rebuilds the radiation coefficient at the returned field. Both the
    /// complete free residual and global heat closure must satisfy the caller's
    /// original absolute watt budgets. The field is never replaced by a hidden
    /// final solve. No radiation tangent, arbitrary-law convergence, uniqueness,
    /// continuum bound or material authority is supplied by this fixed point.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_steady_newton_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        initial: &[f64],
        solver: SteadyContactNewtonConfig,
        patches: &[AmbientRadiationPatch],
        radiation: AmbientRadiationConfig,
    ) -> Result<ContactSteadyNewtonRadiationSolution, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        radiation.validate()?;
        for allowance in [solver.max_iterations, solver.max_residual_evaluations] {
            allowance
                .checked_mul(radiation.max_iterations)
                .ok_or_else(|| {
                    invalid(
                        "<steady-contact-radiation>",
                        "combined radiation/contact work allowance overflows",
                    )
                })?;
        }
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        let mut next_initial = initial.to_vec();
        let mut residual_evaluations = 0usize;
        let mut jacobian_assemblies = 0usize;
        let accepted = solve_ambient_radiation_with(
            cx,
            problem,
            patches,
            radiation,
            |combined_problem| {
                let contact =
                    self.solve_steady_newton(cx, combined_problem, &next_initial, solver)?;
                residual_evaluations = count(residual_evaluations, contact.residual_evaluations)?;
                jacobian_assemblies = count(jacobian_assemblies, contact.jacobian_assemblies)?;
                let krylov_iterations = contact
                    .newton
                    .history
                    .iter()
                    .try_fold(0usize, |total, row| count(total, row.linear_iterations))?;
                next_initial.clone_from(&contact.temperature);
                Ok::<_, ConductionError>(AmbientRadiationTrial {
                    robin_fluxes: contact.robin_fluxes.clone(),
                    robin_out_w: contact.energy.robin_out_w,
                    solid_iterations: contact.newton.iterations,
                    krylov_iterations,
                    solution: Trial {
                        contact,
                        physical_boundary: None,
                    },
                })
            },
            |trial, rows, convective_out_w| {
                let Some(endpoint) = self.radiation_endpoint(cx, problem, &dofs,
                    &trial.contact.temperature, rows, convective_out_w,
                    solver.residual_tolerance_w, solver.energy_tolerance_w)?
                    else { return Ok(false); };
                trial.contact.physical_contacts = endpoint.physical_contacts;
                trial.contact.free_residual_w = endpoint.free_residual_w;
                trial.contact.residual_norm_w = endpoint.residual_norm_w;
                trial.contact.energy = endpoint.energy;
                trial.contact.robin_fluxes = endpoint.robin_fluxes;
                trial.contact.dirichlet_reactions_w = endpoint.dirichlet_reactions_w;
                trial.physical_boundary = Some(endpoint.boundary);
                Ok(true)
            },
        )?;
        let physical_boundary = accepted.solution.physical_boundary.ok_or_else(|| {
            invalid(
                "<steady-contact-radiation>",
                "accepted radiation lost its physical endpoint audit",
            )
        })?;
        Ok(ContactSteadyNewtonRadiationSolution {
            contact: accepted.solution.contact,
            applied_boundary: accepted.combined_boundary,
            physical_boundary,
            convective_robin_fluxes: accepted.convective_robin_fluxes,
            convective_out_w: accepted.convective_out_w,
            radiation: accepted.radiation,
            residual_evaluations,
            jacobian_assemblies,
        })
    }

    /// Compose the original bounded contact fixed point with ambient radiation.
    /// Source surfaces need no unique gradient for this primal solve. Every
    /// callback retains its original source queries and nonlinear bulk solve;
    /// both this route and Newton share the identical physical endpoint audit.
    /// Geometry/nonthermal context stay fixed. This nested fixed point supplies
    /// no convergence, uniqueness or contact/radiation derivative certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_steady_with_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        solver: SolveConfig, coupling: SteadyContactConfig,
        patches: &[AmbientRadiationPatch], radiation: AmbientRadiationConfig,
    ) -> Result<ContactSteadyRadiationSolution, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        radiation.validate()?;
        let (nonlinear, krylov) = coupling.admit(&solver)?;
        for allowance in [coupling.max_iterations, nonlinear, krylov] {
            allowance.checked_mul(radiation.max_iterations).ok_or_else(||
                invalid("<steady-contact-radiation>", "combined radiation/contact work allowance overflows"))?;
        }
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        let mut next_initial = initial.to_vec();
        let mut thermal_responses = 0usize;
        let accepted = solve_ambient_radiation_with(cx, problem, patches, radiation,
            |combined_problem| {
                let contact = self.solve_steady(cx, combined_problem, &next_initial, solver.clone(), coupling)?;
                thermal_responses = count(thermal_responses, contact.iterations)?;
                next_initial.clone_from(&contact.conduction.temperature);
                Ok::<_, ConductionError>(AmbientRadiationTrial {
                    robin_fluxes: contact.robin_fluxes.clone(), robin_out_w: contact.energy.robin_out_w,
                    solid_iterations: contact.nonlinear_iterations, krylov_iterations: contact.krylov_iterations,
                    solution: Trial { contact, physical_boundary: None },
                })
            },
            |trial, rows, convective_out_w| {
                let Some(endpoint) = self.radiation_endpoint(cx, problem, &dofs,
                    &trial.contact.conduction.temperature, rows, convective_out_w,
                    coupling.residual_tolerance_w, coupling.energy_tolerance_w)?
                    else { return Ok(false); };
                trial.contact.physical_contacts = endpoint.physical_contacts;
                trial.contact.free_residual_w = endpoint.free_residual_w;
                trial.contact.residual_norm_w = endpoint.residual_norm_w;
                trial.contact.energy = endpoint.energy;
                trial.contact.robin_fluxes = endpoint.robin_fluxes;
                trial.contact.dirichlet_reactions_w = endpoint.dirichlet_reactions_w;
                trial.physical_boundary = Some(endpoint.boundary);
                Ok(true)
            },
        )?;
        let physical_boundary = accepted.solution.physical_boundary.ok_or_else(||
            invalid("<steady-contact-radiation>", "accepted radiation lost its physical endpoint audit"))?;
        Ok(ContactSteadyRadiationSolution {
            contact: accepted.solution.contact, applied_boundary: accepted.combined_boundary,
            physical_boundary, convective_robin_fluxes: accepted.convective_robin_fluxes,
            convective_out_w: accepted.convective_out_w, radiation: accepted.radiation, thermal_responses,
        })
    }

    /// Shared endpoint check: requery both physical laws without another solve.
    #[allow(clippy::too_many_arguments)]
    fn radiation_endpoint(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, dofs: &DofMap,
        temperature: &[f64], rows: &[AmbientRadiationPatchReport], convective_out_w: f64,
        residual_tolerance_w: f64, energy_tolerance_w: f64,
    ) -> Result<Option<Endpoint>, ConductionError> {
        let boundary = boundary_at_endpoint(problem, rows)?;
        let physical_contacts = self.evaluate(cx, &boundary, temperature)?;
        let system = assemble_operator_scaled_with_interfaces(cx, problem.mesh, &boundary,
            problem.material, problem.source, temperature, None,
            Some(physical_contacts.frozen_interfaces()), problem.element_materials)?;
        let residual = full_residual(&system, temperature);
        for (index, &value) in residual.iter().enumerate() {
            if index % 256 == 0 { poll(cx, index)?; }
            finite(value)?;
        }
        let free_residual_w: Vec<_> = dofs.free().iter().map(|&vertex| (vertex, residual[vertex])).collect();
        let free: Vec<_> = free_residual_w.iter().map(|&(_, value)| value).collect();
        let residual_norm_w = finite(norm2(&free))?;
        let (energy, robin_fluxes) = energy_balance(problem.mesh, &boundary,
            problem.source, &system, dofs, temperature);
        for value in [energy.source_w, energy.neumann_out_w, energy.robin_out_w,
            energy.dirichlet_in_w, energy.closure_w, energy.scale_w] { finite(value)?; }
        let radiation_out_w = rows.iter().try_fold(0.0, |total, row| finite(total + row.nonlinear_heat_w))?;
        // Independent external reservoirs supplement the full nodal residual;
        // contact transfer remains entirely internal on both numerical routes.
        let closure_w = finite(energy.source_w + energy.dirichlet_in_w
            - energy.neumann_out_w - convective_out_w - radiation_out_w)?;
        let split_w = finite(energy.robin_out_w - convective_out_w - radiation_out_w)?;
        if residual_norm_w > residual_tolerance_w || energy.closure_w.abs() > energy_tolerance_w
            || closure_w.abs() > energy_tolerance_w || split_w.abs() > energy_tolerance_w
        { return Ok(None); }
        let dirichlet_reactions_w = dofs.fixed().iter().map(|&vertex| (vertex, residual[vertex])).collect();
        poll(cx, rows.len())?;
        Ok(Some(Endpoint { physical_contacts, free_residual_w, residual_norm_w,
            energy, robin_fluxes, dirichlet_reactions_w, boundary }))
    }
}

/// The shared radiation owner has already measured the actual P1 trace mean.
/// Reuse those physical patch rows rather than introducing another quadrature.
fn boundary_at_endpoint(
    problem: ConductionProblem<'_>,
    rows: &[AmbientRadiationPatchReport],
) -> Result<ThermalBoundary, ConductionError> {
    let mut replacements = Vec::with_capacity(rows.len());
    for row in rows {
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == row.patch.region())
            .ok_or_else(|| {
                invalid(
                    "<steady-contact-radiation>",
                    "accepted patch lost its original boundary",
                )
            })?;
        let (h, reference) = match &problem.boundary.conditions()[region] {
            ThermalBc::Robin {
                htc: ScalarField::Uniform(h),
                t_ref: ScalarField::Uniform(t),
            } => (*h, *t),
            ThermalBc::Neumann {
                outward_flux: ScalarField::Uniform(flux),
            } if *flux == 0.0 => (0.0, row.patch.ambient_temperature_k()),
            _ => {
                return Err(invalid(
                    "<steady-contact-radiation>",
                    "radiation requires its original uniform Robin or zero-flux trace",
                ));
            }
        };
        let radiative = row
            .patch
            .secant_coefficient_w_m2_k(row.mean_surface_temperature_k)?;
        let combined = finite(h + radiative)?;
        let reference = finite(
            (h / combined) * reference + (radiative / combined) * row.patch.ambient_temperature_k(),
        )?;
        replacements.push((region, combined, reference));
    }
    problem
        .boundary
        .with_radiative_robin_replacements(&replacements)
}

fn count(total: usize, increment: usize) -> Result<usize, ConductionError> {
    total.checked_add(increment).ok_or_else(|| {
        invalid(
            "<steady-contact-radiation>",
            "counted coupled work overflowed",
        )
    })
}
