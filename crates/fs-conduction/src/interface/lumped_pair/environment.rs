//! Finite contact and exposed-surface reservoirs on the same caloric endpoints.
//!
//! The native finite-contact solve owns the single signed internal transfer.
//! The native lumped endpoint owns convection and radiation. Bounded fixed-point
//! trials compose those owners without advancing either accepted reservoir.

use super::*;
use crate::lumped::{LumpedThermalEnvironment, endpoint::LumpedEnthalpyEndpoint};

/// Original contact receipt and independently evaluated environmental heat.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactEnvironmentPairReport<R> {
    /// Unmodified contact solve. Its supplied heat includes the held boundary
    /// heat used in the final trial, in addition to each original internal source.
    pub contact: R,
    pub exposed_area_m2: [f64; 2],
    /// Signed contact transfer, debited from A and credited to B exactly once.
    pub heat_a_to_b_j: f64,
    /// Actual endpoint boundary powers integrated over the physical interval.
    pub convection_heat_j: [f64; 2],
    pub radiation_heat_j: [f64; 2],
    /// Per-body storage minus its original internal, contact and boundary heat.
    pub energy_residual_j: [f64; 2],
    /// Combined storage minus original internal and environmental heat.
    pub total_energy_residual_j: f64,
    /// Conservative sum of the worst contact and exposed-environment bounds.
    pub maximum_biot: f64,
    /// Outer boundary-feedback trials; material-contact work is reported inside.
    pub iterations: usize,
}

/// Advance two finite bodies while their uncovered surfaces exchange heat.
///
/// Contact area is removed from BOTH original surface areas before convection
/// and radiation are applied. The already resolved interface law is held over
/// the interval. Each environment retains its original material transport,
/// reservoir temperatures and source admission. Internal heat, finite contact
/// and boundary heat enter one backward-Euler endpoint balance; neither body
/// is first advanced against an isolated bath.
///
/// Each bounded trial reuses `advance_contact_pair` with held boundary heat,
/// then rechecks the original environmental powers at both proposed endpoints.
/// Both local and combined actual energy residuals must pass. The two returned
/// proposals must be published together. Failure or cancellation changes none
/// of the inputs. Convergence is not guaranteed; refine the interval or use the
/// existing spatial/network owner when this uniform-body reduction fails.
///
/// The native whole-surface environment admission remains conservative even
/// with partial exposure. Its exposed contribution is added to the contact
/// Biot bound; separate passing bounds cannot spend the same 0.1 allowance twice.
/// Fixed geometry, negligible interface storage and the original caloric work
/// convention remain assumptions. Algebraic closure does not bound time error.
#[allow(clippy::too_many_arguments)]
pub fn advance_contact_pair_with_environments(
    cx: &Cx<'_>,
    inputs: [UniformEnthalpyStepInput<'_>; 2],
    environments: [&LumpedThermalEnvironment; 2],
    interface: &InterfaceResistance,
    area_m2: f64,
    energy_tolerance_j: f64,
    iteration: ContactIterationConfig,
) -> Result<
    [UniformEnthalpyStep<ContactEnvironmentPairReport<ContactPairReport>>; 2],
    ConductionError,
> {
    advance_with_environments(
        cx,
        inputs,
        environments,
        area_m2,
        energy_tolerance_j,
        iteration,
        |trial| {
            let pair = advance_contact_pair(
                cx,
                trial,
                [&environments[0].transport, &environments[1].transport],
                interface,
                area_m2,
                energy_tolerance_j,
            )?;
            let biot = pair[0]
                .report
                .contact
                .as_ref()
                .map_or(0.0, |c| c.maximum_biot);
            Ok((pair, biot))
        },
    )
}

impl MaterialContactLaw<'_> {
    /// Source-dependent contact and independent exposed-surface baths together.
    ///
    /// The original contact owner pins its selected claim and resolves it at
    /// both actual temperatures/phases on every outer boundary trial. Source
    /// context and both thermal histories remain unchanged until joint acceptance.
    /// `contact_iteration` bounds each inner material solve; `boundary_iteration`
    /// bounds their outer count. Their product bounds the total inner solves.
    /// Coverage, actual heat checks and uniform-body limits are the same as
    /// [`advance_contact_pair_with_environments`].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_environments(
        &self,
        cx: &Cx<'_>,
        inputs: [UniformEnthalpyStepInput<'_>; 2],
        environments: [&LumpedThermalEnvironment; 2],
        area_m2: f64,
        energy_tolerance_j: f64,
        contact_iteration: ContactIterationConfig,
        boundary_iteration: ContactIterationConfig,
    ) -> Result<
        [UniformEnthalpyStep<ContactEnvironmentPairReport<MaterialContactPairReport>>; 2],
        ConductionError,
    > {
        advance_with_environments(
            cx,
            inputs,
            environments,
            area_m2,
            energy_tolerance_j,
            boundary_iteration,
            |trial| {
                let pair = self.advance(
                    cx,
                    trial,
                    [&environments[0].transport, &environments[1].transport],
                    area_m2,
                    energy_tolerance_j,
                    contact_iteration,
                )?;
                let biot = pair[0].report.endpoint_maximum_biot;
                Ok((pair, biot))
            },
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_with_environments<R: Clone>(
    cx: &Cx<'_>,
    inputs: [UniformEnthalpyStepInput<'_>; 2],
    environments: [&LumpedThermalEnvironment; 2],
    area_m2: f64,
    energy_tolerance_j: f64,
    iteration: ContactIterationConfig,
    mut contact: impl FnMut(
        [UniformEnthalpyStepInput<'_>; 2],
    ) -> Result<([UniformEnthalpyStep<R>; 2], f64), ConductionError>,
) -> Result<[UniformEnthalpyStep<ContactEnvironmentPairReport<R>>; 2], ConductionError> {
    poll(cx)?;
    admit_inputs(inputs, area_m2, energy_tolerance_j)?;
    if iteration.max_iterations == 0
        || !positive(iteration.relaxation)
        || iteration.relaxation > 1.0
        || !positive(iteration.residual_tolerance_j)
    {
        return Err(invalid(
            "boundary feedback needs bounded iterations, relaxation in (0,1] and a positive joule allowance",
        ));
    }
    let endpoints = [
        environments[0].prepare_endpoint(inputs[0])?,
        environments[1].prepare_endpoint(inputs[1])?,
    ];
    let exposed_area_m2 = inputs.map(|input| input.surface_area_m2 - area_m2);
    let fractions = [0, 1].map(|side| exposed_area_m2[side] / inputs[side].surface_area_m2);
    let environment_biot = (fractions[0] * endpoints[0].maximum_biot())
        .max(fractions[1] * endpoints[1].maximum_biot());
    let local_tolerances =
        environments.map(|e| energy_tolerance_j.min(e.maximum_thermal_residual_j));
    let mut trial = inputs.map(|input| input.initial);
    for index in 0..iteration.max_iterations {
        poll(cx)?;
        let held = boundary_heat(&endpoints, environments, trial, fractions)?;
        let mut proposed_inputs = inputs;
        for side in 0..2 {
            proposed_inputs[side].internal_heat_j += held[side][0] + held[side][1];
        }
        let (pair, contact_biot) = contact(proposed_inputs)?;
        let states = [pair[0].state, pair[1].state];
        let actual = boundary_heat(&endpoints, environments, states, fractions)?;
        let boundary = actual.map(|heat| heat[0] + heat[1]);
        let external = [0, 1].map(|side| pair[side].external_heat_j + boundary[side]);
        let residuals = [0, 1].map(|side| {
            inputs[side].mass_kg
                * (states[side].specific_enthalpy_j_kg()
                    - inputs[side].initial.specific_enthalpy_j_kg())
                - inputs[side].internal_heat_j
                - external[side]
        });
        let total = residuals[0] + residuals[1];
        let maximum_biot = contact_biot + environment_biot;
        if !maximum_biot.is_finite() || maximum_biot > LUMPED_BIOT_CEILING {
            return Err(invalid(
                "combined contact and environmental Biot exceeds the uniform-body limit; use spatial thermal transport",
            ));
        }
        if !external
            .into_iter()
            .chain(residuals)
            .chain([total])
            .all(f64::is_finite)
        {
            return Err(invalid(
                "combined contact and environmental heat is not representable",
            ));
        }
        let boundary_residual = (0..2)
            .map(|side| (boundary[side] - (held[side][0] + held[side][1])).abs())
            .fold(0.0_f64, f64::max);
        if boundary_residual <= iteration.residual_tolerance_j
            && (0..2).all(|side| residuals[side].abs() <= local_tolerances[side])
            && total.abs() <= energy_tolerance_j
        {
            poll(cx)?;
            let report = ContactEnvironmentPairReport {
                contact: pair[0].report.clone(),
                exposed_area_m2,
                heat_a_to_b_j: pair[1].external_heat_j,
                convection_heat_j: actual.map(|heat| heat[0]),
                radiation_heat_j: actual.map(|heat| heat[1]),
                energy_residual_j: residuals,
                total_energy_residual_j: total,
                maximum_biot,
                iterations: index + 1,
            };
            return Ok([0, 1].map(|side| UniformEnthalpyStep {
                state: states[side],
                external_heat_j: external[side],
                energy_residual_tolerance_j: local_tolerances[side],
                report: report.clone(),
            }));
        }
        for side in 0..2 {
            let h = (1.0 - iteration.relaxation) * trial[side].specific_enthalpy_j_kg()
                + iteration.relaxation * states[side].specific_enthalpy_j_kg();
            trial[side] = inputs[side]
                .curve
                .state_at_specific_enthalpy(h)
                .map_err(|error| invalid(&error.to_string()))?;
        }
    }
    Err(invalid(
        "boundary feedback exhausted its bounded iterations before closing both actual heat balances",
    ))
}

fn boundary_heat(
    endpoints: &[LumpedEnthalpyEndpoint<'_>; 2],
    environments: [&LumpedThermalEnvironment; 2],
    states: [EquilibriumPhaseState; 2],
    fractions: [f64; 2],
) -> Result<[[f64; 2]; 2], ConductionError> {
    let mut heat = [[0.0; 2]; 2];
    for side in 0..2 {
        // Use the exact accepted/trial h. A round trip through a joule increment
        // could move it to a neighboring represented caloric state.
        let sample = endpoints[side].evaluate_specific_enthalpy(
            states[side].specific_enthalpy_j_kg(),
            environments[side].temperature_k,
            0.0,
        )?;
        heat[side] = [
            fractions[side] * sample.convection_heat_j,
            fractions[side] * sample.radiation_heat_j,
        ];
    }
    Ok(heat)
}
