//! Two-field goal comparisons for the actual reflected enclosure equations.
//!
//! Keep original convection and prescribed fluxes in the solid assembly. The
//! physical residual owner adds radiative loads directly, and the complete
//! thermal transpose adds G (dq/dmean) W on every physical radiation patch.

use super::{
    Enclosure, EnclosureRadiationGradient, GrayDiffuseEnclosure, OpenDiffuseEnclosure,
    RadiosityReport, checked, invalid, poll, require_temperature,
};
use crate::adjoint::robin::admit_linear;
use crate::adjoint::{DiscreteGoalComparison, MeanFluxFeedback, RobinResponse};
use crate::assemble::assemble_operator_scaled_with_interfaces;
use crate::interface::temperature::TemperatureContacts;
use crate::radiation::{enclosure_exchange, equilibrium};
use crate::{
    ConductionError, ConductionProblem, DofMap, LinearConfig, ScalarField, ThermalBc,
    ThermalInterfaces,
};
use fs_exec::Cx;

enum Contacts<'a, 'm, 'c> {
    Fixed(Option<&'a ThermalInterfaces>),
    Live(&'a TemperatureContacts<'m, 'c>, usize),
}

impl EnclosureRadiationGradient {
    /// Compare a linear nodal goal under closed reflected radiation and actual
    /// area-mean-dependent convection, with fixed or absent contact resistance.
    ///
    /// In selected original Robin order, reference rows are `[h, r, dh/dmean,
    /// dr/dmean]`, approximate rows are `[h, r]`, and the row-major reference
    /// map is ADDITIONAL `dr_i/dmean_j`. The law owner supplies the actual values,
    /// domains and derivatives. Unselected convection is fixed; the selected
    /// list may be empty. Every radiating patch retains its exact original
    /// whole Robin or Neumann region, including authored nonzero fluxes.
    ///
    /// Both fields receive fresh physical radiosity, material assembly and
    /// direct radiative nodal loads. The reference independently passes the
    /// supplied absolute free-node residual and global energy targets in W,
    /// even for a zero objective. Each radiosity solve also passes the explicit
    /// heat allowance. A shifted reference `r-q/h` is never used for admission,
    /// residuals or the tangent, so rounding cannot hide physical radiation.
    ///
    /// One complete thermal transpose contains K'(T), consistent convection
    /// mass, local/reference feedback and all reflected cross-patch radiation
    /// derivatives. Only the reference requires interior emissivity support
    /// and a smooth material tangent; the approximation needs physical values.
    /// Geometry, source loads, fixed contact and prescribed temperatures stay
    /// unchanged. Original conduction must retain a thermal anchor.
    ///
    /// The result retains signed physical nodal contributions and the observed
    /// finite remainder. It is a discrete refinement diagnostic, not a continuum
    /// error bound. At most 64 Robin controls and 64 radiation patches are
    /// admitted; `4*n*controls + 2*n*patches` factor entries must fit the given
    /// allowance. `dual_iterations` includes the radiosity transpose work used
    /// to build the full patch Jacobian and the complete thermal goal solve.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_goal_at(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: &[f64],
        enclosure: &GrayDiffuseEnclosure,
        heat_tolerance_w: f64,
        full_nodal_weights: &[f64],
        max_feedback_entries: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        compare_goal_fields(
            cx,
            problem,
            Contacts::Fixed(interfaces),
            linear,
            reference_temperature,
            approximate_temperature,
            residual_tolerance_w,
            energy_tolerance_w,
            regions,
            reference_points,
            approximate_values,
            reference_feedback,
            Enclosure::Closed(enclosure),
            heat_tolerance_w,
            full_nodal_weights,
            max_feedback_entries,
        )
    }

    /// Open-enclosure form of [`Self::compare_goal_at`]. The black bath keeps
    /// its original absolute temperature and reflected paths at both fields.
    /// It is an algebraic radiation reservoir, never a solid degree of freedom
    /// or a second contribution to the solid's external heat account.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_open_goal_at(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: &[f64],
        enclosure: &OpenDiffuseEnclosure,
        ambient_temperature_k: f64,
        heat_tolerance_w: f64,
        full_nodal_weights: &[f64],
        max_feedback_entries: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        compare_goal_fields(
            cx,
            problem,
            Contacts::Fixed(interfaces),
            linear,
            reference_temperature,
            approximate_temperature,
            residual_tolerance_w,
            energy_tolerance_w,
            regions,
            reference_points,
            approximate_values,
            reference_feedback,
            Enclosure::Open(enclosure, ambient_temperature_k),
            heat_tolerance_w,
            full_nodal_weights,
            max_feedback_entries,
        )
    }
}

impl TemperatureContacts<'_, '_> {
    /// Compare fields with the original live contact source laws and closed
    /// reflected radiation. Robin rows, additional reference feedback, physical
    /// watt gates and factor budgets have the same meaning as
    /// [`EnclosureRadiationGradient::compare_goal_at`].
    ///
    /// Both residuals query the original selected R(T_A,T_B) at their own
    /// matching or nonmatching traces. The reference adds the complete contact
    /// Jacobian to the material/convection/radiation transpose, preserving
    /// independent A/B partials. The approximation needs physical source
    /// support only; source corners or domain endpoints there do not demand
    /// derivatives. No frozen-contact substitute or extra primal solve enters.
    /// `max_face_pairs` separately bounds original contact trace work.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_steady_goal_with_gray_diffuse_enclosure(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        linear: LinearConfig,
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: &[f64],
        enclosure: &GrayDiffuseEnclosure,
        heat_tolerance_w: f64,
        full_nodal_weights: &[f64],
        max_feedback_entries: usize,
        max_face_pairs: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        compare_goal_fields(
            cx,
            problem,
            Contacts::Live(self, max_face_pairs),
            linear,
            reference_temperature,
            approximate_temperature,
            residual_tolerance_w,
            energy_tolerance_w,
            regions,
            reference_points,
            approximate_values,
            reference_feedback,
            Enclosure::Closed(enclosure),
            heat_tolerance_w,
            full_nodal_weights,
            max_feedback_entries,
        )
    }

    /// Open-enclosure form of
    /// [`Self::compare_steady_goal_with_gray_diffuse_enclosure`], retaining the
    /// original black-bath temperature and all multiple-reflection paths.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_steady_goal_with_open_gray_diffuse_enclosure(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        linear: LinearConfig,
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: &[f64],
        enclosure: &OpenDiffuseEnclosure,
        ambient_temperature_k: f64,
        heat_tolerance_w: f64,
        full_nodal_weights: &[f64],
        max_feedback_entries: usize,
        max_face_pairs: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        compare_goal_fields(
            cx,
            problem,
            Contacts::Live(self, max_face_pairs),
            linear,
            reference_temperature,
            approximate_temperature,
            residual_tolerance_w,
            energy_tolerance_w,
            regions,
            reference_points,
            approximate_values,
            reference_feedback,
            Enclosure::Open(enclosure, ambient_temperature_k),
            heat_tolerance_w,
            full_nodal_weights,
            max_feedback_entries,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn compare_goal_fields(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    contacts: Contacts<'_, '_, '_>,
    linear: LinearConfig,
    reference_temperature: &[f64],
    approximate_temperature: &[f64],
    residual_tolerance_w: f64,
    energy_tolerance_w: f64,
    regions: &[&str],
    reference_points: &[[f64; 4]],
    approximate_values: &[[f64; 2]],
    reference_feedback: &[f64],
    enclosure: Enclosure<'_>,
    heat_tolerance_w: f64,
    full_nodal_weights: &[f64],
    max_feedback_entries: usize,
) -> Result<DiscreteGoalComparison, ConductionError> {
    poll(cx, 0)?;
    admit_linear(linear)?;
    let n = problem.mesh.vertex_count();
    let m = regions.len();
    let p = enclosure.surfaces().len();
    let robin_entries = n.checked_mul(4).and_then(|count| count.checked_mul(m));
    let entries = n
        .checked_mul(2)
        .and_then(|count| count.checked_mul(p))
        .and_then(|flux| robin_entries.and_then(|robin| robin.checked_add(flux)));
    if linear.restart == 0
        || !(heat_tolerance_w.is_finite() && heat_tolerance_w > 0.0)
        || !(residual_tolerance_w.is_finite() && residual_tolerance_w > 0.0)
        || !(energy_tolerance_w.is_finite() && energy_tolerance_w > 0.0)
        || m > 64
        || p == 0
        || p > 64
        || reference_points.len() != m
        || approximate_values.len() != m
        || m.checked_mul(m) != Some(reference_feedback.len())
        || entries.is_none_or(|count| {
            count > max_feedback_entries
                || count
                    .checked_mul(8)
                    .is_none_or(|bytes| bytes > isize::MAX as usize)
        })
    {
        return Err(invalid(
            "enclosure goal needs positive watt targets/restart, matching Robin rows/map and at most 64 controls/patches within its combined factor allowance",
        ));
    }
    for (field, values, temperature) in [
        (
            "enclosure goal reference temperature",
            reference_temperature,
            true,
        ),
        (
            "enclosure goal approximate temperature",
            approximate_temperature,
            true,
        ),
        ("enclosure goal nodal weights", full_nodal_weights, false),
    ] {
        if values.len() != n {
            return Err(ConductionError::FieldLength {
                field,
                expected: n,
                found: values.len(),
            });
        }
        for (i, &value) in values.iter().enumerate() {
            if i % 256 == 0 {
                poll(cx, i)?;
            }
            checked(value)?;
            if temperature {
                require_temperature("enclosure-goal", field, value)?;
            }
        }
    }
    for &(vertex, prescribed) in problem.boundary.dirichlet() {
        poll(cx, vertex)?;
        if vertex >= n
            || reference_temperature[vertex].to_bits() != prescribed.to_bits()
            || approximate_temperature[vertex].to_bits() != prescribed.to_bits()
        {
            return Err(invalid(
                "both enclosure goal fields must retain exact prescribed temperatures",
            ));
        }
    }
    if let Contacts::Live(owner, max_face_pairs) = &contacts {
        owner.admit_goal_comparison(
            cx,
            problem,
            reference_temperature,
            approximate_temperature,
            linear,
            residual_tolerance_w,
            energy_tolerance_w,
            full_nodal_weights,
            *max_face_pairs,
        )?;
    }
    for (i, &value) in reference_feedback.iter().enumerate() {
        if i % 256 == 0 {
            poll(cx, i)?;
        }
        checked(value)?;
    }
    let mut reference_rows = [(0, 0.0, 0.0); 64];
    let mut approximate_rows = [(0, 0.0, 0.0); 64];
    let mut h_slopes = [0.0; 64];
    let mut local_reference_slopes = [0.0; 64];
    for (i, &name) in regions.iter().enumerate() {
        poll(cx, i)?;
        for &value in reference_points[i].iter().chain(&approximate_values[i]) {
            checked(value)?;
        }
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|value| value == name)
            .ok_or_else(|| invalid("unknown original Robin control region"))?;
        let [h, r, dh, dr] = reference_points[i];
        let [ha, ra] = approximate_values[i];
        if regions[..i].contains(&name)
            || h <= 0.0
            || ha <= 0.0
            || !matches!(
                problem.boundary.conditions().get(region),
                Some(ThermalBc::Robin {
                    htc: ScalarField::Uniform(_),
                    t_ref: ScalarField::Uniform(_),
                })
            )
        {
            return Err(invalid(
                "unique existing uniform Robin controls and positive actual coefficients are required",
            ));
        }
        reference_rows[i] = (region, h, r);
        approximate_rows[i] = (region, ha, ra);
        h_slopes[i] = dh;
        local_reference_slopes[i] = dr;
    }
    let reference_boundary = problem
        .boundary
        .with_uniform_robin_replacements(&reference_rows[..m])?;
    let approximate_boundary = problem
        .boundary
        .with_uniform_robin_replacements(&approximate_rows[..m])?;
    let reference_problem = ConductionProblem {
        boundary: &reference_boundary,
        ..problem
    };
    let approximate_problem = ConductionProblem {
        boundary: &approximate_boundary,
        ..problem
    };
    // Shared primal trace owner: exact mesh, whole named regions and original
    // uniform Robin/Neumann laws. No radiation is encoded into either boundary.
    enclosure_exchange::bind_surface_traces(cx, reference_problem, enclosure.surfaces())?;
    enclosure_exchange::bind_surface_traces(cx, approximate_problem, enclosure.surfaces())?;
    let reference_means = enclosure_exchange::surface_means(
        cx,
        reference_problem,
        enclosure.surfaces(),
        reference_temperature,
    )?;
    let approximate_means = enclosure_exchange::surface_means(
        cx,
        approximate_problem,
        enclosure.surfaces(),
        approximate_temperature,
    )?;
    let reference_radiosity =
        checked_radiosity(cx, &enclosure, &reference_means, heat_tolerance_w)?;
    let approximate_radiosity =
        checked_radiosity(cx, &enclosure, &approximate_means, heat_tolerance_w)?;
    // The approximate field is only a physical query. A source-support end
    // there must not acquire the reference tangent's two-sided requirement.
    for (i, (surface, &mean)) in enclosure
        .surfaces()
        .iter()
        .zip(&reference_means)
        .enumerate()
    {
        poll(cx, i)?;
        if let Some((lower, upper)) = surface.emissivity().temperature_validity_k()
            && (mean <= lower || mean >= upper)
        {
            return Err(invalid(
                "enclosure reference derivative requires interior emissivity temperature validity",
            ));
        }
    }
    let fixed_interfaces = match &contacts {
        Contacts::Fixed(value) => *value,
        Contacts::Live(_, _) => None,
    };
    let reference_contact = match &contacts {
        Contacts::Live(owner, max_face_pairs) => Some(owner.linearize_temperature(
            cx,
            &reference_boundary,
            reference_temperature,
            *max_face_pairs,
        )?),
        Contacts::Fixed(_) => None,
    };
    let approximate_contact = match &contacts {
        Contacts::Live(owner, _) => {
            Some(owner.evaluate(cx, &approximate_boundary, approximate_temperature)?)
        }
        Contacts::Fixed(_) => None,
    };
    let reference_interfaces = reference_contact
        .as_ref()
        .map(|contact| contact.physical_contacts().frozen_interfaces())
        .or(fixed_interfaces);
    let approximate_interfaces = approximate_contact
        .as_ref()
        .map(|contact| contact.frozen_interfaces())
        .or(fixed_interfaces);
    let (response, reference) = RobinResponse::for_steady_enclosure(
        cx,
        reference_problem,
        fixed_interfaces,
        reference_contact.as_ref(),
        linear,
        reference_temperature,
        regions,
        enclosure.surfaces(),
        reference_radiosity,
        residual_tolerance_w,
        energy_tolerance_w,
    )?;
    let approximate = equilibrium::endpoint(
        cx,
        approximate_problem,
        approximate_interfaces,
        approximate_temperature,
        residual_tolerance_w,
        approximate_radiosity.clone(),
        &approximate_radiosity,
        &|_| Ok(approximate_boundary.clone()),
        Some(enclosure.surfaces()),
    )?;
    for value in [
        approximate.energy.source_w,
        approximate.energy.neumann_out_w,
        approximate.energy.robin_out_w,
        approximate.energy.dirichlet_in_w,
        approximate.energy.closure_w,
        approximate.energy.scale_w,
    ] {
        checked(value)?;
    }
    let radiosity_tangent = enclosure.linearize(cx, &reference_means, heat_tolerance_w)?;
    let mut flux_derivatives = vec![0.0; p * p];
    let mut flux_weights = vec![0.0; p];
    let mut radiosity_iterations = 0usize;
    for i in 0..p {
        poll(cx, i)?;
        flux_weights.fill(0.0);
        flux_weights[i] = 1.0;
        let row = radiosity_tangent.pullback(
            cx,
            &flux_weights,
            linear.tolerance,
            linear.max_iterations,
        )?;
        radiosity_iterations = radiosity_iterations
            .checked_add(row.iterations)
            .ok_or_else(|| invalid("enclosure goal radiosity iteration count overflow"))?;
        // The shared radiosity owner checks every returned physical transpose
        // residual. These rows include all reflected paths, also Robin/Robin.
        for (j, &value) in row.temperatures.iter().enumerate() {
            flux_derivatives[i * p + j] = checked(value)?;
        }
    }
    let names: Vec<_> = enclosure
        .surfaces()
        .iter()
        .map(|surface| surface.name())
        .collect();
    let flux_feedback = MeanFluxFeedback::new(
        cx,
        problem.mesh,
        &reference_boundary,
        &names,
        &flux_derivatives,
        max_feedback_entries
            - robin_entries.ok_or_else(|| invalid("enclosure Robin factor allowance overflow"))?,
    )?;
    let dual = response.pullback_endpoint_with_mean_and_flux_feedback(
        cx,
        &h_slopes[..m],
        &local_reference_slopes[..m],
        reference_feedback,
        Some(&flux_feedback),
        full_nodal_weights,
        max_feedback_entries,
    )?;
    let primal_relative_residual = physical_relative_residual(
        cx,
        reference_problem,
        reference_interfaces,
        reference_temperature,
        enclosure.surfaces(),
        &reference,
    )?;
    let mut nodal_contributions = vec![0.0; n];
    let mut signed_residual_change = 0.0;
    if reference.free_residual_w.len() != approximate.free_residual_w.len() {
        return Err(invalid(
            "enclosure goal fields have different physical free-node layouts",
        ));
    }
    for (i, (&(vertex, a), &(other, b))) in reference
        .free_residual_w
        .iter()
        .zip(&approximate.free_residual_w)
        .enumerate()
    {
        if i % 256 == 0 {
            poll(cx, i)?;
        }
        if vertex != other {
            return Err(invalid("enclosure goal free-node identities differ"));
        }
        let contribution = checked(dual.nodal_load[vertex] * checked(a - b)?)?;
        nodal_contributions[vertex] = contribution;
        signed_residual_change = checked(signed_residual_change + contribution)?;
    }
    let mut signed_goal_change = 0.0;
    for vertex in 0..n {
        if vertex % 256 == 0 {
            poll(cx, vertex)?;
        }
        let difference = checked(reference_temperature[vertex] - approximate_temperature[vertex])?;
        signed_goal_change =
            checked(signed_goal_change + checked(full_nodal_weights[vertex] * difference)?)?;
    }
    let linearization_remainder = checked(signed_goal_change - signed_residual_change)?;
    let dual_iterations = dual
        .iterations
        .checked_add(radiosity_iterations)
        .ok_or_else(|| invalid("enclosure goal total dual iteration count overflow"))?;
    poll(cx, dual_iterations)?;
    Ok(DiscreteGoalComparison {
        signed_residual_change,
        linearization_remainder,
        signed_goal_change,
        nodal_contributions,
        primal_relative_residual,
        dual_relative_residual: dual.relative_residual,
        dual_iterations,
        uses_nonlinear_jacobian: response.uses_nonlinear_jacobian(),
    })
}

fn checked_radiosity(
    cx: &Cx<'_>,
    enclosure: &Enclosure<'_>,
    means: &[f64],
    heat_tolerance_w: f64,
) -> Result<RadiosityReport, ConductionError> {
    let report = enclosure.solve(cx, means)?;
    let areas = match enclosure {
        Enclosure::Closed(model) => model.view_factors().areas_m2(),
        Enclosure::Open(model, _) => model.view_factors().areas_m2(),
    };
    let largest_area = areas.iter().copied().fold(0.0_f64, f64::max);
    for (i, &value) in report
        .radiosity_w_m2
        .iter()
        .chain(&report.irradiation_w_m2)
        .chain(&report.net_outward_flux_w_m2)
        .chain(&report.net_outward_heat_w)
        .enumerate()
    {
        if i % 256 == 0 {
            poll(cx, i)?;
        }
        checked(value)?;
    }
    checked(report.enclosure_energy_scale_w)?;
    if checked(report.linear_residual_max_w_m2 * largest_area)? > heat_tolerance_w
        || checked(report.enclosure_energy_closure_w)?.abs() > heat_tolerance_w
    {
        return Err(invalid(
            "enclosure goal radiosity equations or closed energy sum exceed the original heat allowance",
        ));
    }
    Ok(report)
}

/// Normalize the already-owned physical residual by its eliminated physical
/// load, including original fixed-temperature lifts and direct radiative flux.
/// This is diagnostic only; admission above uses the unscaled watt targets.
fn physical_relative_residual(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    temperature: &[f64],
    surfaces: &[super::RadiationSurface],
    endpoint: &equilibrium::Endpoint,
) -> Result<f64, ConductionError> {
    let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
    let system = assemble_operator_scaled_with_interfaces(
        cx,
        problem.mesh,
        problem.boundary,
        problem.material,
        problem.source,
        temperature,
        None,
        interfaces,
        problem.element_materials,
    )?;
    let mut rhs = dofs.gather(&system.load);
    for (i, &vertex) in dofs.free().iter().enumerate() {
        if i % 256 == 0 {
            poll(cx, i)?;
        }
        let (columns, values) = system.operator.row(vertex);
        for (&column, &value) in columns.iter().zip(values) {
            if dofs.slot_of(column).is_none() {
                rhs[i] = checked(rhs[i] - checked(value * dofs.prescribed()[column])?)?;
            }
        }
    }
    for (surface, &flux) in surfaces
        .iter()
        .zip(&endpoint.report.nonlinear_radiosity.net_outward_flux_w_m2)
    {
        for (i, &slot) in surface.face_slots().iter().enumerate() {
            if i % 256 == 0 {
                poll(cx, i)?;
            }
            let face = &problem.mesh.boundary()[slot];
            let nodal_heat = checked(face.area * (flux / 3.0))?;
            for &vertex in &face.vertices {
                if let Some(row) = dofs.slot_of(vertex as usize) {
                    rhs[row] = checked(rhs[row] - nodal_heat)?;
                }
            }
        }
    }
    let load_norm = checked(fs_solver::norm2(&rhs))?;
    checked(endpoint.report.residual_w / load_norm.max(f64::MIN_POSITIVE))
}
