//! Exact held-boundary exchange on the original Cp, Cv or measured storage law.
//!
//! The native cell owns its energy coordinate, inventory and boundary. Ambient
//! convection uses only its uncovered area; an optional support reservoir keeps
//! its separate heat account. Both baths are held, so their combined linear
//! power has one conductance-weighted temperature and one exact caloric solve.
//! This is an isolated uniform-body proposal, not a second timeline integrator.

use crate::ConductionError;
use crate::lumped::{BiotGate, lumped_error};
use fs_exec::Cx;
use fs_material::thermal::transport::conduction::{CaloricConductionLayer, CaloricHeatExchange};
use fs_thermal::enthalpy::DepositedCell;

/// Held geometry and numerical acceptance for one isolated cell.
///
/// The characteristic length must describe the body's internal conduction
/// distance. No shape, free expansion or contact geometry is inferred from Cp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedCaloricExchangeConfig {
    pub duration_s: f64,
    pub characteristic_length_m: f64,
    /// Absolute tolerance for stored-energy and separate-bath rounding [J].
    /// This admits representation error, not time-discretization error.
    pub energy_tolerance_j: f64,
}

/// One source-bound proposal; neither the input cell nor a timeline is modified.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedCaloricExchange {
    pub cell_id: u64,
    pub caloric_identity: [u8; 32],
    pub duration_s: f64,
    pub maximum_biot: f64,
    pub total_conductance_w_k: f64,
    /// A reduction of the two held linear baths, not a new material state.
    pub effective_reservoir_temperature_k: f64,
    pub caloric: CaloricHeatExchange,
    pub heat_into_body_j: f64,
    pub stored_energy_change_j: f64,
    pub heat_residual_j: f64,
    /// Signed heat from the original uncovered ambient film into the body.
    pub ambient_heat_j: f64,
    /// Signed heat from the original support reservoir into the body.
    pub reservoir_heat_j: f64,
    pub heat_partition_residual_j: f64,
    /// d Q / d e_initial [kg] in the original h, u or sensible coordinate.
    pub heat_initial_energy_derivative_kg: f64,
    /// d Q / d T_ambient [J/K], with both conductances and inventory held.
    pub heat_ambient_temperature_derivative_j_k: f64,
    /// d Q / d T_support [J/K]; zero for an absent or insulated reservoir.
    pub heat_reservoir_temperature_derivative_j_k: f64,
    /// d Q / d duration [W], retaining the unrounded final driving difference.
    pub heat_duration_derivative_w: f64,
}

/// Evolve the original isolated cell under its held ambient/support boundary.
///
/// Original Cp stores h, fixed-density Cv stores u, and measured capacities
/// retain their declared sensible coordinate and reference amounts. No Cp/Cv
/// conversion, resampling, frozen capacity or timestep loop is introduced.
/// Source support and native inventory admission apply to both endpoints.
///
/// The conservative Biot bound uses the largest active surface coefficient and
/// the minimum of EVERY original conductivity knot. The caller supplies the
/// internal characteristic length and retains the original geometry/work
/// convention. This proposal does not advance expansion, contacts with other
/// finite cells, radiation, deposition or any other concurrent physical owner.
/// A nonzero applied source refuses; it is not silently omitted.
///
/// Each bath keeps its own signed heat, including steady heat flowing between
/// unequal baths when the body stores no net energy. Absolute energy and bath
/// partition residuals must satisfy the caller's declared tolerance before a
/// proposal is returned. A combined bath that rounds nonzero original forcing
/// to equilibrium refuses instead of silently erasing heat. Zero duration or
/// zero total conductance preserves the
/// native state, without excusing invalid geometry or source context.
pub fn solve_lumped_caloric_exchange(
    cx: &Cx<'_>,
    cell: &DepositedCell,
    layer: &CaloricConductionLayer<'_, '_>,
    gate: BiotGate,
    config: LumpedCaloricExchangeConfig,
) -> Result<LumpedCaloricExchange, ConductionError> {
    checkpoint(cx)?;
    let boundary = boundary_conductances(cell, config)?;
    let law = layer.caloric_law();
    law.validate_cell(cell).map_err(|error| invalid(format!(
        "original caloric inventory refuses the initial cell: {error:?}"
    )))?;
    let initial_state = layer.state_at_specific_energy(cell.specific_h_j_kg)
        .map_err(|error| invalid(format!("initial caloric/conductivity support: {error:?}")))?;
    let caloric_identity = law.caloric_identity().ok_or_else(|| invalid(
        "the original caloric law has no identity".into()
    ))?;

    let mut minimum_k = f64::INFINITY;
    for (index, knot) in layer.conductivity_curve().knots().iter().enumerate() {
        if index % 64 == 0 {
            checkpoint(cx)?;
        }
        minimum_k = minimum_k.min(knot[1]);
    }
    let maximum_biot =
        (boundary.maximum_coefficient_w_m2_k / minimum_k) * config.characteristic_length_m;
    if !maximum_biot.is_finite()
        || maximum_biot > gate.ceiling()
        || (boundary.maximum_coefficient_w_m2_k > 0.0 && maximum_biot == 0.0)
    {
        return Err(invalid(format!(
            "source-complete Biot {maximum_biot} exceeds ceiling {} or is not representable; use spatial caloric transport",
            gate.ceiling()
        )));
    }
    let conductance = boundary.ambient_w_k + boundary.reservoir_w_k;
    if !conductance.is_finite() {
        return Err(invalid("combined boundary conductance is not representable".into()));
    }
    let (ambient_weight, reservoir_weight, effective_temperature) = if conductance == 0.0 {
        (0.0, 0.0, cell.boundary.ambient_k)
    } else if boundary.reservoir_w_k == 0.0 {
        (1.0, 0.0, cell.boundary.ambient_k)
    } else if boundary.ambient_w_k == 0.0 {
        (0.0, 1.0, boundary.reservoir_temperature_k)
    } else {
        let ambient_weight = boundary.ambient_w_k / conductance;
        let reservoir_weight = boundary.reservoir_w_k / conductance;
        // Interpolate from the colder bath. Equal temperatures stay bitwise
        // equal even when the two rounded conductance fractions do not sum
        // exactly to one; two equal baths must not manufacture heat.
        let effective = if cell.boundary.ambient_k <= boundary.reservoir_temperature_k {
            reservoir_weight.mul_add(
                boundary.reservoir_temperature_k - cell.boundary.ambient_k,
                cell.boundary.ambient_k,
            )
        } else {
            ambient_weight.mul_add(
                cell.boundary.ambient_k - boundary.reservoir_temperature_k,
                boundary.reservoir_temperature_k,
            )
        };
        if ambient_weight == 0.0 || reservoir_weight == 0.0
            || !effective.is_finite() || effective <= 0.0
        {
            return Err(invalid("combined reservoir temperature is not representable".into()));
        }
        (ambient_weight, reservoir_weight, effective)
    };
    // Reducing two baths to one absolute temperature can round a nonzero
    // original forcing to exact equilibrium. Native energy may still resolve
    // that heat even though T cannot. Refuse instead of returning a false
    // zero-heat proposal with a misleading zero residual.
    if config.duration_s > 0.0
        && boundary.ambient_w_k > 0.0
        && boundary.reservoir_w_k > 0.0
        && effective_temperature == initial_state.temperature_k
    {
        let ambient_power = boundary.ambient_w_k
            * (cell.boundary.ambient_k - initial_state.temperature_k);
        let reservoir_power = boundary.reservoir_w_k
            * (boundary.reservoir_temperature_k - initial_state.temperature_k);
        let ambient_roundoff = boundary.ambient_w_k.mul_add(
            cell.boundary.ambient_k - initial_state.temperature_k, -ambient_power,
        );
        let reservoir_roundoff = boundary.reservoir_w_k.mul_add(
            boundary.reservoir_temperature_k - initial_state.temperature_k, -reservoir_power,
        );
        if !ambient_power.is_finite() || !reservoir_power.is_finite()
            || ambient_power + reservoir_power != 0.0
            || ambient_roundoff + reservoir_roundoff != 0.0
            || (ambient_power == 0.0 && cell.boundary.ambient_k != initial_state.temperature_k)
            || (reservoir_power == 0.0
                && boundary.reservoir_temperature_k != initial_state.temperature_k)
        {
            return Err(invalid(
                "combined bath rounds nonzero original boundary forcing to equilibrium".into()
            ));
        }
    }
    let exposure = (conductance / cell.mass_kg) * config.duration_s;
    if !exposure.is_finite()
        || (conductance > 0.0 && config.duration_s > 0.0 && exposure == 0.0)
    {
        return Err(invalid("specific boundary exposure is not representable".into()));
    }
    let caloric = layer.linear_heat_exchange(
        cell.specific_h_j_kg,
        effective_temperature,
        exposure,
    ).map_err(|error| invalid(format!("original caloric exchange: {error:?}")))?;

    let mut endpoint = cell.clone();
    endpoint.specific_h_j_kg = caloric.final_specific_energy_j_kg;
    law.validate_cell(&endpoint).map_err(|error| invalid(format!(
        "original caloric inventory refuses the proposed cell: {error:?}"
    )))?;
    layer.state_at_specific_energy(endpoint.specific_h_j_kg)
        .map_err(|error| invalid(format!("final caloric/conductivity support: {error:?}")))?;

    let heat = cell.mass_kg * caloric.specific_heat_j_kg;
    let stored = cell.mass_kg * caloric.represented_energy_change_j_kg;
    let residual = cell.mass_kg * caloric.heat_residual_j_kg;
    if !heat.is_finite() || !stored.is_finite() || !residual.is_finite()
        || (caloric.specific_heat_j_kg != 0.0 && heat == 0.0)
        || (caloric.represented_energy_change_j_kg != 0.0 && stored == 0.0)
        || residual.abs() > config.energy_tolerance_j
        || (heat - stored).abs() > config.energy_tolerance_j
    {
        return Err(invalid(format!(
            "stored-energy residual {residual} J exceeds the admitted representation or {} J tolerance",
            config.energy_tolerance_j
        )));
    }

    // Separate the net stored heat from the bath-to-bath through-flow. This
    // avoids subtracting a rounded mean temperature in weak exchanges.
    let (ambient_heat, reservoir_heat) = if boundary.ambient_w_k == 0.0 {
        (0.0, heat)
    } else if boundary.reservoir_w_k == 0.0 {
        (heat, 0.0)
    } else {
        let through_heat = ((boundary.ambient_w_k * reservoir_weight) * config.duration_s)
            * (cell.boundary.ambient_k - boundary.reservoir_temperature_k);
        let ambient_heat = ambient_weight.mul_add(heat, through_heat);
        let reservoir_heat = reservoir_weight.mul_add(heat, -through_heat);
        (ambient_heat, reservoir_heat)
    };
    let partition_residual = (ambient_heat + reservoir_heat) - heat;
    if !ambient_heat.is_finite() || !reservoir_heat.is_finite()
        || !partition_residual.is_finite()
        || partition_residual.abs() > config.energy_tolerance_j
    {
        return Err(invalid(format!(
            "separate-bath heat residual {partition_residual} J exceeds {} J tolerance",
            config.energy_tolerance_j
        )));
    }

    let inlet_derivative = -cell.mass_kg * caloric.inlet_energy_loss;
    let ambient_derivative = cell.mass_kg
        * (caloric.reservoir_energy_derivative_j_kg_k * ambient_weight);
    let reservoir_derivative = cell.mass_kg
        * (caloric.reservoir_energy_derivative_j_kg_k * reservoir_weight);
    let duration_derivative = if config.duration_s == 0.0 {
        // At zero time the state is unchanged, but the duration derivative is
        // still the ORIGINAL nonzero boundary power, even if a combined bath
        // temperature rounds to the body's temperature.
        boundary.ambient_w_k.mul_add(
            cell.boundary.ambient_k - initial_state.temperature_k,
            boundary.reservoir_w_k
                * (boundary.reservoir_temperature_k - initial_state.temperature_k),
        )
    } else {
        conductance * caloric.remaining_temperature_difference_k
    };
    if [inlet_derivative, ambient_derivative, reservoir_derivative, duration_derivative]
        .iter().any(|value| !value.is_finite())
    {
        return Err(invalid("boundary heat Jacobian is not representable".into()));
    }
    checkpoint(cx)?;
    Ok(LumpedCaloricExchange {
        cell_id: cell.id.0,
        caloric_identity,
        duration_s: config.duration_s,
        maximum_biot,
        total_conductance_w_k: conductance,
        effective_reservoir_temperature_k: effective_temperature,
        caloric,
        heat_into_body_j: heat,
        stored_energy_change_j: stored,
        heat_residual_j: residual,
        ambient_heat_j: ambient_heat,
        reservoir_heat_j: reservoir_heat,
        heat_partition_residual_j: partition_residual,
        heat_initial_energy_derivative_kg: inlet_derivative,
        heat_ambient_temperature_derivative_j_k: ambient_derivative,
        heat_reservoir_temperature_derivative_j_k: reservoir_derivative,
        heat_duration_derivative_w: duration_derivative,
    })
}

struct BoundaryConductances {
    ambient_w_k: f64,
    reservoir_w_k: f64,
    reservoir_temperature_k: f64,
    maximum_coefficient_w_m2_k: f64,
}

fn boundary_conductances(
    cell: &DepositedCell,
    config: LumpedCaloricExchangeConfig,
) -> Result<BoundaryConductances, ConductionError> {
    let positive = |value: f64| value.is_finite() && value > 0.0;
    let nonnegative = |value: f64| value.is_finite() && value >= 0.0;
    if !positive(cell.mass_kg)
        || !cell.specific_h_j_kg.is_finite()
        || !positive(cell.geometry.volume_m3)
        || !nonnegative(cell.geometry.surface_area_m2)
        || cell.geometry.start_m.iter().chain(&cell.geometry.end_m)
            .any(|coordinate| !coordinate.is_finite())
        || !positive(cell.boundary.ambient_k)
        || !nonnegative(cell.boundary.htc_w_m2_k)
        || !nonnegative(config.duration_s)
        || !positive(config.characteristic_length_m)
        || !nonnegative(config.energy_tolerance_j)
    {
        return Err(invalid("invalid cell, held boundary, duration, length or tolerance".into()));
    }
    if !cell.boundary.source_w.is_finite() || cell.boundary.source_w != 0.0 {
        return Err(invalid("exact reservoir exchange requires zero applied cell source".into()));
    }
    let (reservoir_area, reservoir_w_k, reservoir_temperature_k, reservoir_coefficient) =
        if let Some(reservoir) = cell.boundary.reservoir {
            if !nonnegative(reservoir.area_m2)
                || reservoir.area_m2 > cell.geometry.surface_area_m2
                || !nonnegative(reservoir.conductance_w_k)
                || !positive(reservoir.temperature_k)
                || (reservoir.area_m2 == 0.0 && reservoir.conductance_w_k > 0.0)
            {
                return Err(invalid("invalid original support reservoir area or conductance".into()));
            }
            let coefficient = if reservoir.conductance_w_k == 0.0 {
                0.0
            } else {
                reservoir.conductance_w_k / reservoir.area_m2
            };
            if !coefficient.is_finite()
                || (reservoir.conductance_w_k > 0.0 && coefficient == 0.0)
            {
                return Err(invalid("support surface coefficient is not representable".into()));
            }
            (reservoir.area_m2, reservoir.conductance_w_k, reservoir.temperature_k, coefficient)
        } else {
            (0.0, 0.0, cell.boundary.ambient_k, 0.0)
        };
    let free_area = cell.geometry.surface_area_m2 - reservoir_area;
    let ambient_w_k = cell.boundary.htc_w_m2_k * free_area;
    if !ambient_w_k.is_finite()
        || (cell.boundary.htc_w_m2_k > 0.0 && free_area > 0.0 && ambient_w_k == 0.0)
    {
        return Err(invalid("ambient surface conductance is not representable".into()));
    }
    let ambient_coefficient = if ambient_w_k == 0.0 {
        0.0
    } else {
        cell.boundary.htc_w_m2_k
    };
    Ok(BoundaryConductances {
        ambient_w_k,
        reservoir_w_k,
        reservoir_temperature_k,
        maximum_coefficient_w_m2_k: ambient_coefficient.max(reservoir_coefficient),
    })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "lumped-caloric-exchange",
        at: 0,
    })
}

fn invalid(detail: String) -> ConductionError {
    lumped_error("source-caloric bath", detail)
}
