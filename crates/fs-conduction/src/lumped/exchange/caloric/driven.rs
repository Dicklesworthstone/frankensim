//! Held applied power and original linear baths on one native caloric body.
//! The shared material kernel owns the temperature-dependent energy integral.

use super::{
    LumpedCaloricExchangeConfig, PreparedCaloricBoundary, checkpoint, invalid,
    prepare_caloric_exchange, stored_heat_account,
};
use crate::ConductionError;
use crate::lumped::BiotGate;
use fs_exec::Cx;
use fs_material::thermal::transport::conduction::{
    CaloricConductionLayer, CaloricDrivenHeatExchange,
};
use fs_math::eft::{two_prod, two_sum};
use fs_thermal::enthalpy::DepositedCell;

/// Exact held-boundary proposal including signed applied cell heat.
///
/// Geometry and all boundary powers/coefficients are held over the window.
/// The original material determines h, u or declared sensible storage; no
/// thermodynamic meaning is inferred from the sign of the applied power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedCaloricDrivenExchange {
    pub cell_id: u64,
    pub caloric_identity: [u8; 32],
    pub duration_s: f64,
    pub maximum_biot: f64,
    pub total_conductance_w_k: f64,
    /// Sum of the original bath exposures G_i*dt/m [J/(kg K)].
    pub exposure_j_kg_k: f64,
    pub caloric: CaloricDrivenHeatExchange,
    /// Net heat from both baths and the applied source into the body.
    pub heat_into_body_j: f64,
    pub stored_energy_change_j: f64,
    pub heat_residual_j: f64,
    pub ambient_heat_j: f64,
    pub reservoir_heat_j: f64,
    /// Signed original source_w * duration, counted exactly once.
    pub source_heat_j: f64,
    /// Ambient + support + source minus modeled net heat, compensated [J].
    pub heat_partition_residual_j: f64,
    /// Ambient + support + source minus represented native storage [J].
    pub energy_balance_residual_j: f64,
    /// d Q / d original initial specific energy [kg].
    pub heat_initial_energy_derivative_kg: f64,
    pub heat_ambient_temperature_derivative_j_k: f64,
    pub heat_reservoir_temperature_derivative_j_k: f64,
    /// d Q / d applied source power [s], at fixed initial state and baths.
    pub heat_source_power_derivative_s: f64,
    /// d Q / d duration [W], including the original source at zero time.
    pub heat_duration_derivative_w: f64,
}

/// Exchange with held ambient/support baths AND a constant signed heat source.
///
/// The exact equation is C_specific(T)*dT/ds = Q0 - X*(T-T_initial),
/// for normalized time s in [0,1]. Q0 is the accurately accumulated initial
/// driving heat over the window, and X is the sum of the original exposures.
/// It never constructs an absolute equivalent-bath temperature. A negative
/// mathematical equilibrium therefore does not assert negative physical
/// kelvin; the actual trajectory must stay inside original material support.
///
/// Zero conductance is an exact native-energy impulse, with the original
/// capacity-weighted time-mean temperature. It uses no artificial tiny film.
/// Weak signed forcing survives an unchanged rounded absolute temperature.
/// Source and each bath keep separate heat accounts and analytic responses.
///
/// The original inventory, conductivity/Biot and endpoint admission are shared
/// with the passive entry point. This returns a proposal without advancing a
/// timeline or any concurrent source, expansion, radiation or finite-body
/// contact owner. Unrepresentable integrated forcing or heat accounts refuse.
pub fn solve_lumped_caloric_driven_exchange(
    cx: &Cx<'_>,
    cell: &DepositedCell,
    layer: &CaloricConductionLayer<'_, '_>,
    gate: BiotGate,
    config: LumpedCaloricExchangeConfig,
) -> Result<LumpedCaloricDrivenExchange, ConductionError> {
    try_solve_driven_exchange(cx, cell, layer, gate, config)?.ok_or_else(|| {
        invalid("integrated initial forcing or exposure is not representable".into())
    })
}

// None is only an arithmetic-range refusal during drive preparation. The
// passive compatibility entry may retain its older finite-exposure route for
// a saturated window whose initial driving heat would overflow. Source-domain
// and inventory failures never trigger that fallback.
pub(super) fn try_solve_driven_exchange(
    cx: &Cx<'_>,
    cell: &DepositedCell,
    layer: &CaloricConductionLayer<'_, '_>,
    gate: BiotGate,
    config: LumpedCaloricExchangeConfig,
) -> Result<Option<LumpedCaloricDrivenExchange>, ConductionError> {
    let prepared = prepare_caloric_exchange(cx, cell, layer, gate, config)?;
    let Some(drive) = DrivenInputs::new(cell, &prepared, config.duration_s) else {
        return Ok(None);
    };
    let caloric = layer.driven_heat_exchange(
        cell.specific_h_j_kg,
        drive.initial_j_kg,
        drive.exposure_j_kg_k,
    ).map_err(|error| invalid(format!("original driven caloric exchange: {error:?}")))?;
    let [heat, stored, residual] = stored_heat_account(
        cell, layer,
        caloric.final_specific_energy_j_kg,
        caloric.specific_heat_j_kg,
        caloric.represented_energy_change_j_kg,
        caloric.heat_residual_j_kg,
        config.energy_tolerance_j,
    )?;

    // Retain relative mean travel. Subtracting two rounded absolute
    // temperatures would erase boundary loss in weak source-driven heating.
    let [ambient_specific, reservoir_specific] = drive.bath_heats(
        caloric.specific_heat_j_kg,
        caloric.mean_temperature_change_k,
        cell.boundary.ambient_k - prepared.boundary.reservoir_temperature_k,
    );
    let ambient_heat = cell.mass_kg * ambient_specific;
    let reservoir_heat = cell.mass_kg * reservoir_specific;
    let source_heat = cell.boundary.source_w * config.duration_s;
    if !ambient_heat.is_finite() || !reservoir_heat.is_finite() || !source_heat.is_finite()
        || (ambient_specific != 0.0 && ambient_heat == 0.0)
        || (reservoir_specific != 0.0 && reservoir_heat == 0.0)
        || (cell.boundary.source_w != 0.0 && config.duration_s > 0.0 && source_heat == 0.0)
    {
        return Err(invalid("individual source or bath heat is not representable".into()));
    }
    let partition = compensated_sum(&[ambient_heat, reservoir_heat, source_heat, -heat])
        .ok_or_else(|| invalid("source/bath heat summation is not representable".into()))?;
    let balance = compensated_sum(&[ambient_heat, reservoir_heat, source_heat, -stored])
        .ok_or_else(|| invalid("native-energy heat summation is not representable".into()))?;
    if partition.abs() > config.energy_tolerance_j
        || balance.abs() > config.energy_tolerance_j
    {
        return Err(invalid(format!(
            "source/bath partition {partition} J or native balance {balance} J exceeds {} J tolerance",
            config.energy_tolerance_j
        )));
    }

    let source_gain = caloric.source_heat_gain;
    let inlet_derivative = -cell.mass_kg * caloric.inlet_energy_loss;
    let ambient_derivative = cell.mass_kg * (drive.ambient_exposure * source_gain);
    let reservoir_derivative = cell.mass_kg * (drive.reservoir_exposure * source_gain);
    let source_derivative = config.duration_s * source_gain;
    // A relative endpoint drive may underflow in J/kg before conversion to
    // watts. Use the original physical power times the retained attenuation
    // whenever that power is representable, including at zero duration.
    let duration_derivative = if let Some(power) = original_initial_power(cell, &prepared) {
        power * caloric.inlet_energy_gain
    } else if config.duration_s > 0.0 {
        product_ratio(caloric.remaining_drive_j_kg, cell.mass_kg, config.duration_s)
            .ok_or_else(|| invalid("endpoint source/bath power is not representable".into()))?
    } else {
        return Err(invalid("initial source/bath power is not representable".into()));
    };
    if [inlet_derivative, ambient_derivative, reservoir_derivative,
        source_derivative, duration_derivative]
        .iter().any(|value| !value.is_finite())
    {
        return Err(invalid("driven boundary heat Jacobian is not representable".into()));
    }
    checkpoint(cx)?;
    Ok(Some(LumpedCaloricDrivenExchange {
        cell_id: cell.id.0,
        caloric_identity: prepared.caloric_identity,
        duration_s: config.duration_s,
        maximum_biot: prepared.maximum_biot,
        total_conductance_w_k: prepared.total_conductance_w_k,
        exposure_j_kg_k: drive.exposure_j_kg_k,
        caloric,
        heat_into_body_j: heat,
        stored_energy_change_j: stored,
        heat_residual_j: residual,
        ambient_heat_j: ambient_heat,
        reservoir_heat_j: reservoir_heat,
        source_heat_j: source_heat,
        heat_partition_residual_j: partition,
        energy_balance_residual_j: balance,
        heat_initial_energy_derivative_kg: inlet_derivative,
        heat_ambient_temperature_derivative_j_k: ambient_derivative,
        heat_reservoir_temperature_derivative_j_k: reservoir_derivative,
        heat_source_power_derivative_s: source_derivative,
        heat_duration_derivative_w: duration_derivative,
    }))
}

struct DrivenInputs {
    ambient_exposure: f64,
    reservoir_exposure: f64,
    exposure_j_kg_k: f64,
    ambient_initial: (f64, f64),
    reservoir_initial: (f64, f64),
    initial_j_kg: f64,
    source_j_kg: f64,
}

impl DrivenInputs {
    fn new(cell: &DepositedCell, prepared: &PreparedCaloricBoundary, duration: f64) -> Option<Self> {
        let boundary = prepared.boundary;
        let ambient_exposure = product_ratio(boundary.ambient_w_k, duration, cell.mass_kg)?;
        let reservoir_exposure = product_ratio(boundary.reservoir_w_k, duration, cell.mass_kg)?;
        let source = product_ratio(cell.boundary.source_w, duration, cell.mass_kg)?;
        let exposure_j_kg_k = ambient_exposure + reservoir_exposure;
        if !exposure_j_kg_k.is_finite() {
            return None;
        }
        let ambient_gap = cell.boundary.ambient_k - prepared.initial_temperature_k;
        let reservoir_gap = boundary.reservoir_temperature_k - prepared.initial_temperature_k;
        let ambient_initial = two_prod(ambient_exposure, ambient_gap);
        let reservoir_initial = two_prod(reservoir_exposure, reservoir_gap);
        if (ambient_exposure != 0.0 && ambient_gap != 0.0 && ambient_initial.0 == 0.0)
            || (reservoir_exposure != 0.0 && reservoir_gap != 0.0 && reservoir_initial.0 == 0.0)
        {
            return None;
        }
        let initial_j_kg = compensated_sum(&[
            ambient_initial.0, reservoir_initial.0, source,
            ambient_initial.1, reservoir_initial.1,
        ])?;
        Some(Self {
            ambient_exposure,
            reservoir_exposure,
            exposure_j_kg_k,
            ambient_initial,
            reservoir_initial,
            initial_j_kg,
            source_j_kg: source,
        })
    }

    fn bath_heats(&self, heat: f64, mean_travel: f64, bath_difference: f64) -> [f64; 2] {
        let exposures = [self.ambient_exposure, self.reservoir_exposure];
        let initials = [self.ambient_initial, self.reservoir_initial];
        let mut accounts = [0.0; 2];
        let total = self.exposure_j_kg_k;
        let weights = if total == 0.0 { [0.0; 2] } else {
            exposures.map(|value| value / total)
        };
        // Use the larger fraction with the smaller exposure so the parallel
        // through-flow coefficient does not lose a weak but representable bath.
        let cross = if exposures[0] <= exposures[1] {
            exposures[0] * weights[1]
        } else {
            exposures[1] * weights[0]
        };
        let through = cross * bath_difference;
        let bath_heat = heat - self.source_j_kg;
        for side in 0..2 {
            let feedback = exposures[side] * mean_travel;
            let direct = (-exposures[side]).mul_add(mean_travel, initials[side].0)
                + initials[side].1;
            // Relative mean travel is best for weak source heating: subtracting
            // source heat from net heat would lose the small bath loss.
            let direct_scale = initials[side].0.abs() + initials[side].1.abs()
                + feedback.abs();
            let signed_through = if side == 0 { through } else { -through };
            let balanced = weights[side].mul_add(bath_heat, signed_through);
            // At saturation the mean may round to the bath, while finite stored
            // heat survives. Net bath heat plus through-flow retains that heat.
            // Choose the algebra with the smaller rounding-amplification scale.
            let balanced_scale = weights[side] * (heat.abs() + self.source_j_kg.abs())
                + through.abs();
            accounts[side] = if total > 0.0 && balanced.is_finite()
                && (exposures[side] == 0.0 || weights[side] > 0.0)
                && balanced_scale < direct_scale
            {
                balanced
            } else {
                direct
            };
        }
        accounts
    }
}

fn original_initial_power(
    cell: &DepositedCell,
    prepared: &PreparedCaloricBoundary,
) -> Option<f64> {
    let ambient = two_prod(
        prepared.boundary.ambient_w_k,
        cell.boundary.ambient_k - prepared.initial_temperature_k,
    );
    let reservoir = two_prod(
        prepared.boundary.reservoir_w_k,
        prepared.boundary.reservoir_temperature_k - prepared.initial_temperature_k,
    );
    compensated_sum(&[
        ambient.0, reservoir.0, cell.boundary.source_w, ambient.1, reservoir.1,
    ])
}

// Three possible multiplication orders avoid introducing an overflowing
// intermediate when the finite physical product/ratio itself is representable.
// A nonzero drive or exposure must not silently become zero.
fn product_ratio(a: f64, b: f64, denominator: f64) -> Option<f64> {
    if a == 0.0 || b == 0.0 {
        return Some(0.0);
    }
    for value in [(a / denominator) * b, a * (b / denominator), (a * b) / denominator] {
        if value.is_finite() && value != 0.0 {
            return Some(value);
        }
    }
    None
}

fn compensated_sum(values: &[f64]) -> Option<f64> {
    let mut sum = 0.0;
    let mut correction = 0.0;
    for &value in values {
        if !value.is_finite() {
            return None;
        }
        let (next, error) = two_sum(sum, value);
        if !next.is_finite() {
            return None;
        }
        sum = next;
        correction += error;
    }
    let result = sum + correction;
    result.is_finite().then_some(result)
}
