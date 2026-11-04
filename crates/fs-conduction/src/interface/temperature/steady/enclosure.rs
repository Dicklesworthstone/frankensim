//! Live source contacts inside the original closed or open radiosity iteration.

use fs_exec::Cx;
use crate::{ConductionError, ConductionProblem, RobinFlux, SolveConfig, ThermalBoundary};
use crate::radiation::{CoupledRadiationConfig, CoupledRadiationReport, GrayDiffuseEnclosure,
    OpenDiffuseEnclosure, RadiationSurface, RadiosityReport, enclosure_exchange as exchange, equilibrium,
    require_enclosure_solid_residual, solve_radiating_surfaces_with};
use super::super::{TemperatureContactEvaluation, finite, invalid, poll};
use super::{ContactSteadyNewtonSolution, ContactSteadySolution, SteadyContactConfig,
    SteadyContactNewtonConfig, TemperatureContacts};

/// Actual work across the initial conduction response and every radiative trial.
#[derive(Debug, Clone, Copy, Default)]
pub struct ContactEnclosureWork {
    /// Complete spatial responses, including every inner contact Picard solve.
    pub thermal_responses: usize,
    /// Complete bulk/contact nonlinear updates, including discarded trials.
    pub nonlinear_iterations: usize,
    /// Actual Krylov iterations in those responses.
    pub krylov_iterations: usize,
    /// Complete contact Newton residual calls; zero for the Picard route.
    pub residual_evaluations: usize,
    /// Complete contact Newton tangents; zero for the Picard route.
    pub jacobian_assemblies: usize,
}

impl ContactEnclosureWork {
    fn add(&mut self, responses: usize, nonlinear: usize, krylov: usize,
        residuals: usize, jacobians: usize) -> Result<(), ConductionError>
    {
        for (total, amount) in [(&mut self.thermal_responses, responses),
            (&mut self.nonlinear_iterations, nonlinear), (&mut self.krylov_iterations, krylov),
            (&mut self.residual_evaluations, residuals), (&mut self.jacobian_assemblies, jacobians)]
        { *total = count(*total, amount)?; }
        Ok(())
    }
}

/// An unchanged accepted field, its physical endpoint and original work history.
#[derive(Debug)]
pub struct ContactSteadyEnclosureSolution<T> {
    /// The final contact response with freshly audited physical rows. Numerical
    /// history still describes its applied load, not a monolithic radiation solve.
    pub contact: T,
    /// Exact shifted Robin/Neumann load used by the final numerical response.
    pub applied_boundary: ThermalBoundary,
    /// Original convection alone; internal radiation is never heat into air.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Original convective outward heat at the accepted field, W.
    pub convective_out_w: f64,
    /// Original reflected radiosity iteration and its fresh equilibrium audit.
    pub radiation: CoupledRadiationReport,
    /// All work, including discarded radiation and contact trials.
    pub work: ContactEnclosureWork,
}

impl TemperatureContacts<'_, '_> {
    /// Compose complete source-contact/bulk Newton with the existing closed,
    /// area-mean gray enclosure. The same prepared radiosity operator supplies
    /// each applied and fresh endpoint load. Geometry, source selection and
    /// manufactured context stay fixed; no contact or radiation law is frozen
    /// across responses. The initial response retains the original solid anchor.
    ///
    /// Acceptance requeries contacts at the unchanged returned field and adds
    /// actual radiative nodal loads to the original conductive/convection rows,
    /// including prescribed reactions. Independent free-residual, heat closure,
    /// patch mismatch and radiosity gates must pass the original watt budgets.
    /// This bounded composition supplies no radiation tangent, arbitrary-law
    /// convergence, uniqueness or continuum error certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_steady_newton_with_gray_diffuse_enclosure(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        solver: SteadyContactNewtonConfig, enclosure: &GrayDiffuseEnclosure,
        radiation: CoupledRadiationConfig,
    ) -> Result<ContactSteadyEnclosureSolution<ContactSteadyNewtonSolution>, ConductionError> {
        poll(cx, 0)?;
        let prepared = enclosure.prepare_radiosity(cx)?;
        self.solve_steady_newton_with_enclosure(cx, problem, initial, solver,
            enclosure.surfaces(), radiation, largest_area(enclosure.view_factors().areas_m2()),
            &|t| prepared.solve(cx, t))
    }

    /// Compose the same complete contact Newton response with an open gray
    /// cavity. Its explicit black reservoir has no solid degree of freedom;
    /// only physical-surface radiation enters the solid heat/reaction account.
    /// Original convection alone remains available to transported air.
    /// Source support, original thermal anchoring and all closed-path physical
    /// acceptance gates remain required. This supplies no coupled adjoint.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_steady_newton_with_open_gray_diffuse_enclosure(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        solver: SteadyContactNewtonConfig, enclosure: &OpenDiffuseEnclosure,
        ambient_temperature_k: f64, radiation: CoupledRadiationConfig,
    ) -> Result<ContactSteadyEnclosureSolution<ContactSteadyNewtonSolution>, ConductionError> {
        poll(cx, 0)?;
        require_ambient(ambient_temperature_k)?;
        self.solve_steady_newton_with_enclosure(cx, problem, initial, solver,
            enclosure.surfaces(), radiation, largest_area(enclosure.view_factors().areas_m2()),
            &|t| enclosure.solve(cx, t, ambient_temperature_k).map(|r| r.exchange))
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_steady_newton_with_enclosure<R>(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        solver: SteadyContactNewtonConfig, surfaces: &[RadiationSurface],
        radiation: CoupledRadiationConfig, largest_area: f64, evaluate: &R,
    ) -> Result<ContactSteadyEnclosureSolution<ContactSteadyNewtonSolution>, ConductionError>
    where R: Fn(&[f64]) -> Result<RadiosityReport, ConductionError>,
    {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        let responses = count(radiation.max_iterations, 1)?;
        for allowance in [solver.max_iterations, solver.max_residual_evaluations] {
            bound(allowance, responses)?;
        }
        let krylov_allowance = solver.newton.linear_restart
            .checked_mul(solver.newton.max_linear_cycles)
            .and_then(|n| n.checked_mul(solver.max_iterations))
            .ok_or_else(|| invalid("<steady-contact-enclosure>", "combined Newton/Krylov allowance overflows"))?;
        bound(krylov_allowance, responses)?;
        let patches = exchange::bind_surface_traces(cx, problem, surfaces)?;
        let boundary = |radiosity: &RadiosityReport| exchange::shifted_boundary(problem, &patches, radiosity);
        let mut work = ContactEnclosureWork::default();
        let mut convection = None;
        let (contact, radiation, applied_boundary) = solve_radiating_surfaces_with(cx, problem,
            surfaces, radiation, &boundary, evaluate,
            |combined, warm| {
                let contact = self.solve_steady_newton(cx, combined, warm.unwrap_or(initial), solver)?;
                let krylov = contact.newton.history.iter()
                    .try_fold(0usize, |n, row| count(n, row.linear_iterations))?;
                work.add(1, contact.newton.iterations, krylov,
                    contact.residual_evaluations, contact.jacobian_assemblies)?;
                Ok(contact)
            },
            |contact| &contact.temperature,
            |contact, nonlinear, applied| {
                let (physical, endpoint, accepted) = self.enclosure_endpoint(cx, problem,
                    &contact.temperature, surfaces, nonlinear, applied, &boundary,
                    solver.residual_tolerance_w, solver.energy_tolerance_w, largest_area)?;
                if accepted {
                    convection = Some((endpoint.robin_fluxes.clone(), endpoint.convective_out_w));
                    contact.physical_contacts = physical;
                    contact.free_residual_w = endpoint.free_residual_w;
                    contact.residual_norm_w = endpoint.report.residual_w;
                    contact.energy = endpoint.energy;
                    contact.robin_fluxes = endpoint.robin_fluxes;
                    contact.dirichlet_reactions_w = endpoint.dirichlet_reactions_w;
                }
                Ok((endpoint.report, accepted))
            })?;
        let (convective_robin_fluxes, convective_out_w) = convection.ok_or_else(missing_audit)?;
        Ok(ContactSteadyEnclosureSolution { contact, applied_boundary, convective_robin_fluxes,
            convective_out_w, radiation, work })
    }

    /// Use the original bounded contact Picard responses in the same enclosure
    /// owner. Source surfaces need no unique gradient for this primal route.
    /// Its final physical gates are identical to the complete Newton route;
    /// every inner spatial response remains counted, including rejected trials.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_steady_with_gray_diffuse_enclosure(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        solver: SolveConfig, coupling: SteadyContactConfig,
        enclosure: &GrayDiffuseEnclosure, radiation: CoupledRadiationConfig,
    ) -> Result<ContactSteadyEnclosureSolution<ContactSteadySolution>, ConductionError> {
        poll(cx, 0)?;
        let prepared = enclosure.prepare_radiosity(cx)?;
        self.solve_steady_with_enclosure(cx, problem, initial, solver, coupling,
            enclosure.surfaces(), radiation, largest_area(enclosure.view_factors().areas_m2()),
            &|t| prepared.solve(cx, t))
    }

    /// Retain source-surface Picard contact responses inside the same open
    /// cavity owner. The black reservoir is an external heat channel, never a
    /// fictitious solid trace or a load to air. Every inner response is counted
    /// and the unchanged final field passes fresh contact/radiation watt gates.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_steady_with_open_gray_diffuse_enclosure(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        solver: SolveConfig, coupling: SteadyContactConfig,
        enclosure: &OpenDiffuseEnclosure, ambient_temperature_k: f64,
        radiation: CoupledRadiationConfig,
    ) -> Result<ContactSteadyEnclosureSolution<ContactSteadySolution>, ConductionError> {
        poll(cx, 0)?;
        require_ambient(ambient_temperature_k)?;
        self.solve_steady_with_enclosure(cx, problem, initial, solver, coupling,
            enclosure.surfaces(), radiation, largest_area(enclosure.view_factors().areas_m2()),
            &|t| enclosure.solve(cx, t, ambient_temperature_k).map(|r| r.exchange))
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_steady_with_enclosure<R>(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, initial: &[f64],
        solver: SolveConfig, coupling: SteadyContactConfig,
        surfaces: &[RadiationSurface], radiation: CoupledRadiationConfig,
        largest_area: f64, evaluate: &R,
    ) -> Result<ContactSteadyEnclosureSolution<ContactSteadySolution>, ConductionError>
    where R: Fn(&[f64]) -> Result<RadiosityReport, ConductionError>,
    {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        let responses = count(radiation.max_iterations, 1)?;
        let (nonlinear, krylov) = coupling.admit(&solver)?;
        for allowance in [coupling.max_iterations, nonlinear, krylov] { bound(allowance, responses)?; }
        let patches = exchange::bind_surface_traces(cx, problem, surfaces)?;
        let boundary = |radiosity: &RadiosityReport| exchange::shifted_boundary(problem, &patches, radiosity);
        let mut work = ContactEnclosureWork::default();
        let mut convection = None;
        let (contact, radiation, applied_boundary) = solve_radiating_surfaces_with(cx, problem,
            surfaces, radiation, &boundary, evaluate,
            |combined, warm| {
                let contact = self.solve_steady(cx, combined, warm.unwrap_or(initial), solver.clone(), coupling)?;
                require_enclosure_solid_residual(cx, &contact.conduction)?;
                work.add(contact.iterations, contact.nonlinear_iterations, contact.krylov_iterations, 0, 0)?;
                Ok(contact)
            },
            |contact| &contact.conduction.temperature,
            |contact, nonlinear, applied| {
                let (physical, endpoint, accepted) = self.enclosure_endpoint(cx, problem,
                    &contact.conduction.temperature, surfaces, nonlinear, applied, &boundary,
                    coupling.residual_tolerance_w, coupling.energy_tolerance_w, largest_area)?;
                if accepted {
                    convection = Some((endpoint.robin_fluxes.clone(), endpoint.convective_out_w));
                    contact.physical_contacts = physical;
                    contact.free_residual_w = endpoint.free_residual_w;
                    contact.residual_norm_w = endpoint.report.residual_w;
                    contact.energy = endpoint.energy;
                    contact.robin_fluxes = endpoint.robin_fluxes;
                    contact.dirichlet_reactions_w = endpoint.dirichlet_reactions_w;
                }
                Ok((endpoint.report, accepted))
            })?;
        let (convective_robin_fluxes, convective_out_w) = convection.ok_or_else(missing_audit)?;
        Ok(ContactSteadyEnclosureSolution { contact, applied_boundary, convective_robin_fluxes,
            convective_out_w, radiation, work })
    }

    #[allow(clippy::too_many_arguments)]
    fn enclosure_endpoint<F>(&self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        temperature: &[f64], surfaces: &[RadiationSurface], nonlinear: RadiosityReport,
        applied: &RadiosityReport, boundary: &F, residual_tolerance_w: f64, energy_tolerance_w: f64,
        largest_area: f64,
    ) -> Result<(TemperatureContactEvaluation, equilibrium::Endpoint, bool), ConductionError>
    where F: Fn(&RadiosityReport) -> Result<ThermalBoundary, ConductionError>,
    {
        let physical = self.evaluate(cx, problem.boundary, temperature)?;
        let endpoint = equilibrium::endpoint(cx, problem, Some(physical.frozen_interfaces()),
            temperature, residual_tolerance_w, nonlinear, applied, boundary, Some(surfaces))?;
        let energy = &endpoint.energy;
        for value in [energy.source_w, energy.neumann_out_w, energy.robin_out_w,
            energy.dirichlet_in_w, energy.closure_w, energy.scale_w, endpoint.convective_out_w]
        { finite(value)?; }
        let mut accepted = endpoint.report.residual_w <= residual_tolerance_w
            && energy.closure_w.abs() <= energy_tolerance_w
            && endpoint.report.max_patch_heat_mismatch_w <= energy_tolerance_w;
        for radiosity in [applied, &endpoint.report.nonlinear_radiosity] {
            accepted &= finite(radiosity.enclosure_energy_closure_w)?.abs() <= energy_tolerance_w
                && finite(radiosity.linear_residual_max_w_m2 * largest_area)? <= energy_tolerance_w;
        }
        Ok((physical, endpoint, accepted))
    }
}

fn largest_area(areas: &[f64]) -> f64 {
    areas.iter().copied().fold(0.0, f64::max)
}

fn require_ambient(temperature_k: f64) -> Result<(), ConductionError> {
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(invalid("<steady-contact-enclosure>", "black reservoir temperature must be positive and finite"));
    }
    Ok(())
}

fn count(total: usize, increment: usize) -> Result<usize, ConductionError> {
    total.checked_add(increment).ok_or_else(|| invalid("<steady-contact-enclosure>", "counted work overflowed"))
}

fn bound(allowance: usize, responses: usize) -> Result<(), ConductionError> {
    allowance.checked_mul(responses).map(|_| ()).ok_or_else(||
        invalid("<steady-contact-enclosure>", "combined enclosure/contact work allowance overflows"))
}

fn missing_audit() -> ConductionError {
    invalid("<steady-contact-enclosure>", "accepted enclosure lost its physical endpoint audit")
}
