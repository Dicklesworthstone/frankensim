//! Original conducting surfaces on one common finite gas and distinct bodies.
//!
//! The native multiple-body owner supplies the only gas inversion, heat solve,
//! temporal controller and physical account. Each original material surface
//! supplies its own nonlinear coefficient at that shared gas endpoint. No
//! sequential reservoir updates or body-to-body heat redistribution occur.

use super::*;
use fs_couple::finite_gas_thermal::bodies::{
    FiniteGasBodiesEnvironment, FiniteGasBodiesError, FiniteGasBody, MAX_BODIES,
    window::{
        FiniteGasBodiesWindowCandidate, FiniteGasBodiesWindowError, FiniteGasBodiesWindowSolve,
    },
};

/// One physical body and its original conducting material/film surface.
///
/// Bodies may share immutable material sources while retaining distinct IDs,
/// caloric states, geometry and prescribed heat. The surface must bind this
/// body's exact caloric object. Each positive film covers the complete original
/// body area; an explicitly zero-area film insulates it.
#[derive(Clone, Copy, Debug)]
pub struct ConjugateGasBody<'b, 'a, 'm, F: ConjugateGasFilm> {
    pub body: FiniteGasBody<'b>,
    pub surface: &'b ConjugateGasSurface<'a, 'm, F>,
}

type WindowError = FiniteGasBodiesWindowError<WallConvectionError<ConjugateGasError>>;

fn source_error(id: u64, error: ConjugateGasError) -> WindowError {
    match error {
        ConjugateGasError::Thermal(DepositionError::Cancelled) => WindowError::Cancelled,
        error => WindowError::Step(FiniteGasBodiesError::Body {
            id,
            error: ConvectionError::Convection(WallConvectionError::Source(error)),
        }),
    }
}

/// Advance several material bodies and their conducting surfaces against ONE gas.
///
/// Each original surface is refreshed during every native coarse, fine and
/// retried interval. Its layer-work audit uses the actual interval duration and
/// one body share of the authored source-work tolerance, so the combined layer
/// allowance is independent of body count. Initial and final states of EVERY
/// proposal must remain supported by their original caloric/conductivity source,
/// even for insulation and rejected coarse steps.
///
/// Native IDs, masses, geometry, source energies and the shared gas inventory
/// remain fixed for this window. The controller alone apportions prescribed
/// window heats and retains accepted fine trajectories. Source failures and
/// cancellation return no partially advanced gas/body collection. Sources must
/// be repeatable; effects outside the physical owners cannot be rolled back.
///
/// The original pure-phase, quasi-steady normal-wall, complete exposed-area,
/// disabled body-node radiation and lumped-core/Biot requirements still apply.
/// This adds no lateral gas transport, layer storage or new temporal-order claim.
pub fn adaptive_isochoric_bodies_candidate<F, G>(
    cx: &Cx<'_>,
    environment: &FiniteGasBodiesEnvironment<'_, G>,
    bindings: &[ConjugateGasBody<'_, '_, '_, F>],
    controls: FiniteGasBodiesWindowSolve,
    wall: WallConvectionSolve,
) -> Result<
    FiniteGasBodiesWindowCandidate<
        WallConvectionSample<ConjugateGasConvectionSample<F::Sample, Temperature>>,
        G,
    >,
    FiniteGasBodiesWindowError<WallConvectionError<ConjugateGasError>>,
>
where
    F: ConjugateGasFilm,
    G: IsochoricGasStorage<State = F::State>,
{
    cx.checkpoint().map_err(|_| WindowError::Cancelled)?;
    // Bound the only preparation vectors before reserving; native history and
    // numerical work remain governed by the original shared window controls.
    if environment.maximum_bodies == 0
        || environment.maximum_bodies > MAX_BODIES
        || bindings.is_empty()
        || bindings.len() > environment.maximum_bodies
    {
        return Err(WindowError::Invalid(
            "bounded nonempty conducting body set required",
        ));
    }
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_unstable_by_key(|binding| binding.body.id);
    if ordered
        .windows(2)
        .any(|pair| pair[0].body.id == pair[1].body.id)
    {
        return Err(WindowError::Invalid(
            "conducting surfaces require unique physical body IDs",
        ));
    }
    for binding in &ordered {
        cx.checkpoint().map_err(|_| WindowError::Cancelled)?;
        binding
            .surface
            .admit_isochoric(
                environment.gas.state(),
                binding.body.environment,
                binding.body.input,
            )
            .map_err(|error| source_error(binding.body.id, error))?;
    }
    let bodies: Vec<_> = ordered.iter().map(|binding| binding.body).collect();
    let source_tolerance_j = environment.exchange.source_tolerance_j / bodies.len() as f64;
    let surface = |id| {
        ordered
            .binary_search_by_key(&id, |binding| binding.body.id)
            .map(|index| ordered[index].surface)
            .map_err(|_| {
                ConjugateGasError::Invalid("native body has no original conducting surface")
            })
    };
    environment
        .advance_adaptive_with_wall_convection_checked(
            cx,
            &bodies,
            controls,
            wall,
            |body, gas, temperature| {
                surface(body.id)?.resolve_isochoric_coefficient(
                    cx,
                    gas,
                    temperature,
                    body.input,
                    source_tolerance_j,
                )
            },
            |body, state| {
                let surface = surface(body.id)?;
                surface.path.admit_input(body.input)?;
                surface.path.admit_state(state)
            },
        )
        .map_err(|error| match error {
            WindowError::Step(FiniteGasBodiesError::Body {
                error:
                    ConvectionError::Convection(WallConvectionError::Source(
                        ConjugateGasError::Thermal(DepositionError::Cancelled),
                    )),
                ..
            }) => WindowError::Cancelled,
            other => other,
        })
}
