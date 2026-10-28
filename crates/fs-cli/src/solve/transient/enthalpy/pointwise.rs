//! Pointwise radiation uses the native enthalpy Newton solve and endpoint audit.
//! No frozen-radiation iteration or separate caloric history is introduced.
use fs_conduction::{
    PointwiseAmbientRadiation, PointwiseEnclosureRadiation, PointwiseEnclosureTolerance,
};
use fs_conduction::radiation::GrayDiffuseEnclosure;
use fs_conduction::transient::enthalpy::EnthalpyEndpointEvaluation;
use super::*;

enum Owner<'a> {
    Ambient(PointwiseAmbientRadiation<'a>),
    Enclosure(PointwiseEnclosureRadiation<'a>),
}

impl<'a> Owner<'a> {
    fn new(
        cx: &Cx<'_>, problem: ConductionProblem<'a>,
        lowering: &'a radiation::LoweredRadiation,
        enclosure: Option<&'a GrayDiffuseEnclosure>,
    ) -> Result<Self, SolveRefusal> {
        if !lowering.is_pointwise() || lowering.is_enclosure() != enclosure.is_some() {
            return Err(bad("pointwise radiation lost its exact declared spatial model"));
        }
        match enclosure {
            Some(model) => PointwiseEnclosureRadiation::new(
                cx, problem.mesh, problem.boundary, model,
                PointwiseEnclosureTolerance {
                    absolute_heat_w: lowering.config.balance_tolerance_w,
                    relative_heat: lowering.config.balance_relative_tolerance,
                },
            ).map(Self::Enclosure).map_err(lower),
            None => PointwiseAmbientRadiation::new(
                cx, problem.mesh, problem.boundary, &lowering.patches,
            ).map(Self::Ambient).map_err(lower),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self, cx: &Cx<'_>, engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        old: &[f64], candidate: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    ) -> Result<EnthalpyEndpointEvaluation, SolveRefusal> {
        match self {
            Self::Ambient(owner) => engine.evaluate_endpoint_with_pointwise_ambient_radiation(
                cx, problem, interfaces, old, candidate, dt_s, config, owner,
            ),
            Self::Enclosure(owner) => engine.evaluate_endpoint_with_pointwise_enclosure_radiation(
                cx, problem, interfaces, old, candidate, dt_s, config, owner,
            ),
        }.map_err(failed)
    }

    #[allow(clippy::too_many_arguments)]
    fn advance(
        &self, cx: &Cx<'_>, engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        old: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepSolution, SolveRefusal> {
        match self {
            Self::Ambient(owner) => engine.advance_with_pointwise_ambient_radiation(
                cx, problem, interfaces, old, dt_s, config, owner,
            ),
            Self::Enclosure(owner) => engine.advance_with_pointwise_enclosure_radiation(
                cx, problem, interfaces, old, dt_s, config, owner,
            ),
        }.map_err(failed)
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    cx: &Cx<'_>, engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    old: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    lowering: &radiation::LoweredRadiation, enclosure: Option<&GrayDiffuseEnclosure>,
    deadline: Option<(Instant, f64)>,
) -> Result<radiative::Step, SolveRefusal> {
    poll(cx, deadline)?;
    // Reuse the prepared face geometry/radiosity factorization for this whole
    // Newton solve and both physical audits. An outer convection trial binds its
    // own boundary, while all attempts retain this exact old enthalpy.
    let owner = Owner::new(cx, problem, lowering, enclosure)?;
    let initial = engine.prescribed_endpoint_initial_guess(cx, problem.boundary, old)
        .map_err(|error| match error {
            HeterogeneousEnthalpyError::Enthalpy(inner) => failed(inner),
            other => bad(format!("pointwise initial bath state refused: {other}")),
        })?;
    let initial = owner.evaluate(cx, engine, problem, interfaces, old, &initial, dt_s, config)?;
    let threshold = finite(config.newton.absolute_tolerance.max(
        config.newton.relative_tolerance * initial.residual_norm_j.max(f64::MIN_POSITIVE),
    ))?;
    let endpoint = owner.advance(cx, engine, problem, interfaces, old, dt_s, config)?;
    let physical = owner.evaluate(
        cx, engine, problem, interfaces, old, &endpoint.specific_enthalpy_j_kg, dt_s, config,
    )?;
    if physical.residual_norm_j > threshold || physical.energy_residual_j.abs() > config.energy_tolerance_j {
        return Err(bad("pointwise radiation failed its original physical residual or energy allowance"));
    }
    let krylov_iterations = endpoint.newton.history.iter().try_fold(0usize, |count, row| {
        count.checked_add(row.linear_iterations)
            .ok_or_else(|| bad("pointwise Newton Krylov work count overflows"))
    })?;
    let exchange = radiation::pointwise_receipt(
        endpoint.pointwise_radiation.as_ref(), endpoint.pointwise_enclosure_radiation.as_ref(),
        lowering, physical.convective_out_w,
    )?;
    let evidence = format!(
        "{{\"model\":\"pointwise-p1\",\"solid_solves\":1,\"solid_iterations\":{},\"krylov_iterations\":{krylov_iterations},\"physical_residual_j\":{},\"physical_residual_threshold_j\":{},\"physical_energy_residual_j\":{},\"physical_dirichlet_in_w\":{},\"convective_out_w\":{},\"prescribed_neumann_out_w\":{},\"nonlinear_radiative_out_w\":{},\"exchange\":{exchange}}}",
        endpoint.newton.iterations, number(physical.residual_norm_j)?, number(threshold)?,
        number(physical.energy_residual_j)?, number(physical.dirichlet_in_w)?,
        number(physical.convective_out_w)?, number(physical.neumann_out_w)?,
        number(physical.radiation_out_w)?,
    );
    let net_input_w = finite(physical.source_w + physical.dirichlet_in_w
        - physical.neumann_out_w - physical.convective_out_w - physical.radiation_out_w)?;
    poll(cx, deadline)?;
    Ok(radiative::Step {
        updates: endpoint.newton.iterations,
        endpoint, ambient: None, enclosure: None, evidence,
        solid_solves: 1, krylov_iterations, net_input_w,
        residual_norm_j: physical.residual_norm_j,
        residual_tolerance_j: threshold,
        energy_residual_j: physical.energy_residual_j,
        neumann_out_w: physical.neumann_out_w,
        dirichlet_in_w: physical.dirichlet_in_w,
    })
}

#[allow(clippy::too_many_arguments)]
pub(in crate::solve::transient) fn evaluate(
    cx: &Cx<'_>, engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    old: &[f64], candidate: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    lowering: &radiation::LoweredRadiation, enclosure: Option<&GrayDiffuseEnclosure>,
) -> Result<EnthalpyEndpointEvaluation, SolveRefusal> {
    Owner::new(cx, problem, lowering, enclosure)?
        .evaluate(cx, engine, problem, interfaces, old, candidate, dt_s, config)
}
