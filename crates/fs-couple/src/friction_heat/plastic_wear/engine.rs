//! Shared native caloric iteration for the two concrete sliding/wear owners.
//! Mechanical payloads remain their original receipts; this is no new law.

use std::cell::Cell;

use fs_solver::NewtonKrylovState;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};

use super::*;
use super::super::{FrictionHeatInput, FrictionHeatStep};

pub(super) struct MechanicalWork {
    pub contact_work_residual_j: f64,
    pub frictional_work_j: f64,
    pub temperature_parameter_work_j: f64,
    pub wear_parameter_work_j: f64,
    pub normal_plastic_work_j: f64,
    pub normal_work_j: f64,
}

/// Audit each original normal continuation and actual tangential work separately.
pub(super) fn mechanical_work<E>(
    temperature: &FiniteGapPlasticContinuationReceipt,
    wear: &[&FiniteGapPlasticContinuationReceipt],
    cells: &[PlasticCellSliding],
    frictional_work_j: f64,
    tolerance_j: f64,
    cx: &Cx<'_>,
) -> Result<MechanicalWork, PlasticSlidingHeatError<E>> {
    let mut maximum = 0.0_f64;
    for normal in std::iter::once(temperature).chain(wear.iter().copied()) {
        let residual = finite(
            normal.elastic_energy_change_j + normal.response.plastic_work_j
                - normal.normal_work_j
                - normal.parameter_work_j,
        )?;
        maximum = maximum.max(residual.abs());
    }
    let mut work = 0.0;
    let mut correction = 0.0;
    for (i, cell) in cells.iter().enumerate() {
        if i % 64 == 0 {
            checkpoint(cx)?;
        }
        let response = &cell.response;
        let supplied = finite(
            response.traction_n[0] * response.slip_increment_m[0]
                + response.traction_n[1] * response.slip_increment_m[1],
        )?;
        let residual = finite(response.stored_change_j + response.dissipated_j - supplied)?;
        maximum = maximum.max(residual.abs());
        let next = work + response.dissipated_j;
        correction += if work.abs() >= response.dissipated_j.abs() {
            (work - next) + response.dissipated_j
        } else {
            (response.dissipated_j - next) + work
        };
        work = next;
    }
    maximum = maximum.max(finite(work + correction - frictional_work_j)?.abs());
    if maximum > tolerance_j {
        return Err(PlasticSlidingHeatError::Invalid(
            "native contact work exceeds its independent energy tolerance",
        ));
    }
    let mut wear_parameter_work_j = 0.0;
    let mut normal_plastic_work_j = temperature.response.plastic_work_j;
    let mut normal_work_j = temperature.normal_work_j;
    for normal in wear {
        wear_parameter_work_j = finite(wear_parameter_work_j + normal.parameter_work_j)?;
        normal_plastic_work_j = finite(normal_plastic_work_j + normal.response.plastic_work_j)?;
        normal_work_j = finite(normal_work_j + normal.normal_work_j)?;
    }
    Ok(MechanicalWork {
        contact_work_residual_j: maximum,
        frictional_work_j,
        temperature_parameter_work_j: temperature.parameter_work_j,
        wear_parameter_work_j,
        normal_plastic_work_j,
        normal_work_j,
    })
}

pub(super) struct PreparedMechanics<S, M> {
    pub source: S,
    pub mechanics: M,
    pub work: MechanicalWork,
    pub removed_volume_m3: [f64; 2],
    pub removed_volume_fraction: [f64; 2],
    pub refreshed_contact_area_m2: f64,
    pub endpoint_contact_area_m2: f64,
    pub normal_plastic_heat_j: [f64; 2],
    pub nonthermal_plastic_work_j: f64,
}

struct Staged<P, S, M> {
    prepared: PreparedMechanics<S, M>,
    payload: P,
    conductivity_w_m_k: [f64; 2],
    contact_resistance: Option<PlasticThermalResistanceSample>,
}

pub(super) struct Solved<'t, P, S, M> {
    pub source: S,
    pub thermal: DepositionSimulation<'t>,
    pub successor_route: Option<FrictionHeatRoute>,
    pub report: PlasticSlidingHeatReport<P, M>,
}

fn numerical<E>(
    error: FiniteDifferenceError<PlasticSlidingHeatError<E>>,
) -> PlasticSlidingHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => PlasticSlidingHeatError::Invalid(message),
    }
}

/// Original finite-difference/Newton and finite-cell/interface residual, shared
/// without extrapolating any native state or assuming a new heat partition.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn solve<'t, P, E, S, M>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    route: &FrictionHeatRoute,
    previous_contact_area_m2: f64,
    varying_route: bool,
    controls: PlasticSlidingHeatControls,
    energy_guess_j: &[f64],
    interface_temperature_guess_k: f64,
    prepare: impl Fn(&DepositionSimulation<'t>, f64, f64)
        -> Result<PreparedMechanics<S, M>, PlasticSlidingHeatError<E>>,
    sample_resistance: impl Fn(&PreparedMechanics<S, M>, f64)
        -> Result<Option<PlasticThermalResistanceSample>, PlasticSlidingHeatError<E>>,
    resolve: impl Fn(&DepositionSimulation<'t>, f64)
        -> Result<PlasticSlidingHeatSource<P>, E>,
    cx: &Cx<'_>,
) -> Result<Solved<'t, P, S, M>, PlasticSlidingHeatError<E>> {
    let bound_area = route.contact_area_m2();
    if (previous_contact_area_m2 - bound_area).abs()
        > 64.0 * f64::EPSILON * previous_contact_area_m2.max(bound_area)
    {
        return Err(PlasticSlidingHeatError::Invalid(
            "accepted native pressure area differs from its thermal contact",
        ));
    }
    controls.newton.validate()?;
    let cells = thermal.cells().len();
    let dimension = cells.checked_add(1).ok_or(PlasticSlidingHeatError::Invalid(
        "thermal dimension overflow",
    ))?;
    if controls.maximum_cells == 0
        || cells > controls.maximum_cells
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || !interface_temperature_guess_k.is_finite()
        || interface_temperature_guess_k <= 0.0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.contact_work_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(PlasticSlidingHeatError::Invalid(
            "finite bounded thermal unknowns, iteration controls and independent energy tolerances required",
        ));
    }
    let check = || checkpoint(cx);
    let step = FrictionHeatStep::new(thermal, target_s, route, check)?;
    let duration_s = step.duration_s();
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        check()?;
        if trials.get() >= controls.maximum_trials {
            return Err(PlasticSlidingHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(PlasticSlidingHeatError::Invalid(
                "complete finite energy/interface-temperature iterate required",
            ));
        }
        let interface_temperature = finite(
            interface_temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k,
        )?;
        let prepare_input = |endpoint: &DepositionSimulation<'t>, temperature| {
            let prepared = prepare(endpoint, temperature, duration_s)?;
            Ok::<_, PlasticSlidingHeatError<E>>((prepared.endpoint_contact_area_m2, prepared))
        };
        let heat_input = |endpoint: &DepositionSimulation<'t>,
                          temperature,
                          prepared: PreparedMechanics<S, M>| {
            let PlasticSlidingHeatSource {
                conductivity_w_m_k,
                transport,
            } = resolve(endpoint, temperature).map_err(PlasticSlidingHeatError::Source)?;
            let contact_resistance = sample_resistance(&prepared, temperature)?;
            let resistance = contact_resistance.as_ref().map_or(
                route.area_resistance_m2_k_w(),
                |sample| sample.area_resistance_m2_k_w,
            );
            let StepCoefficients {
                mut sources_w,
                conductances_w_k,
                payload,
            } = transport;
            for (side, cell) in route.cells().into_iter().enumerate() {
                let heat_j = prepared.normal_plastic_heat_j[side];
                if heat_j > 0.0 {
                    let power_w = finite(heat_j / duration_s)?;
                    if power_w == 0.0 {
                        return Err(DepositionError::Arithmetic.into());
                    }
                    let source_w = sources_w.entry(cell).or_insert(0.0);
                    *source_w = finite(*source_w + power_w)?;
                }
            }
            Ok::<_, PlasticSlidingHeatError<E>>((
                FrictionHeatInput {
                    // The interface receives friction once; declared normal
                    // plastic heat enters bodies. Parameter work is not heat.
                    dissipated_j: prepared.work.frictional_work_j,
                    conductivity_w_m_k,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            prepared,
                            payload,
                            conductivity_w_m_k,
                            contact_resistance,
                        },
                    },
                },
                resistance,
            ))
        };
        if varying_route {
            step.evaluate_with_area_and_resistance(
                &x[..cells],
                interface_temperature,
                prepare_input,
                heat_input,
                false,
                check,
            )
            .map(|area_trial| {
                let (trial, route) = area_trial.into_trial_and_route();
                (trial, Some(route))
            })
        } else {
            step.evaluate(
                &x[..cells],
                interface_temperature,
                |endpoint, temperature| {
                    let (_, prepared) = prepare_input(endpoint, temperature)?;
                    heat_input(endpoint, temperature, prepared).map(|(input, _)| input)
                },
                check,
            )
            .map(|trial| (trial, None))
        }
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x).map(|(trial, _)| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.push(0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        cx.checkpoint().is_err() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    check()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(PlasticSlidingHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(PlasticSlidingHeatError::NotConverged(result.report));
    }
    let (trial, successor_route) = evaluate(&solver.x)?;
    let (thermal, staged, heat) = trial.accept(controls.heat, check)?;
    let Staged {
        prepared,
        payload,
        conductivity_w_m_k,
        contact_resistance,
    } = staged;
    let PreparedMechanics {
        source,
        mechanics,
        work,
        removed_volume_m3,
        removed_volume_fraction,
        refreshed_contact_area_m2,
        endpoint_contact_area_m2,
        normal_plastic_heat_j,
        nonthermal_plastic_work_j,
    } = prepared;
    let other_heat_j = finite(heat.external_heat_j - normal_plastic_heat_j[0] - normal_plastic_heat_j[1])?;
    check()?;
    Ok(Solved {
        source,
        thermal,
        successor_route,
        report: PlasticSlidingHeatReport {
            mechanics,
            heat,
            source: payload,
            conductivity_w_m_k,
            contact_resistance,
            previous_contact_area_m2,
            refreshed_contact_area_m2,
            endpoint_contact_area_m2,
            nonlinear: result.report,
            trials: trials.get(),
            contact_work_residual_j: work.contact_work_residual_j,
            frictional_work_j: work.frictional_work_j,
            temperature_parameter_work_j: work.temperature_parameter_work_j,
            wear_parameter_work_j: work.wear_parameter_work_j,
            normal_plastic_work_j: work.normal_plastic_work_j,
            normal_plastic_heat_j,
            nonthermal_plastic_work_j,
            other_heat_j,
            normal_work_j: work.normal_work_j,
            removed_volume_m3,
            removed_volume_fraction,
        },
    })
}
