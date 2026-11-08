//! Pressure-selected gas flow with source-density thermal acceleration.
//!
//! Every trial reuses the native stationary species, NASA-9 enthalpy, convection
//! and hydraulic solve. Only the inlet mass rate changes. A constant-area path
//! without side injections has dp = sum(dp_loss) + beta*(mdot/A)*(u_out - u_in).
//! Thermodynamic pressure remains the source's common reference; the mechanical
//! pressure envelope and omitted kinetic/dissipative powers are explicitly bounded.

use fs_exec::Cx;
use fs_material::gas_field::steady::{MAX_STEADY_GAS_CELLS, SteadyAxialGas, SteadyGasSolve};
use fs_qty::Pressure;

use super::{
    GasCellTransport, GasTransportFeedback, MaterialSteadyGasError, MaterialSteadyGasSolution,
    check, positive, solve_material_gas_transport,
};
use crate::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::PressureDrivenBudget;

/// Original pressure controls, declared velocity profile and shared work/energy bounds.
#[derive(Clone, Copy, Debug)]
pub struct GasPressureBudget {
    /// Positive flow bracket, residual tolerances, full-solve count and omitted
    /// hydraulic-power allowance. Every evaluated trial must be supported.
    pub pressure: PressureDrivenBudget,
    /// Constant momentum-flux factor beta = mean(u²)/mean(u)², at least one.
    pub momentum_flux_factor: f64,
    /// Constant kinetic-energy-flux factor alpha = mean(u³)/mean(u)³.
    /// Unidirectional profiles require alpha >= beta² (up to roundoff).
    /// Plug flow uses (1, 1); developed circular laminar flow uses (4/3, 2).
    /// These are explicit model inputs, never inferred from a heat correlation.
    pub kinetic_energy_flux_factor: f64,
    /// Total actual native thermal substeps over ALL flow/feedback trials.
    pub maximum_substeps: usize,
    /// Bound on the sum of absolute regional kinetic-energy flux changes [W].
    /// These are omitted from the native sensible/chemical enthalpy equation.
    pub maximum_neglected_kinetic_power_w: f64,
}

/// One pressure-closed native solution, without a second publication solve.
#[derive(Debug)]
pub struct PressureDrivenGasSolution {
    /// Complete original thermal/species field and local material responses.
    pub gas: MaterialSteadyGasSolution,
    pub mass_flow_kg_per_s: f64,
    pub driving_pressure_pa: f64,
    /// Retained cross-sectional momentum and kinetic-energy profile factors.
    pub momentum_flux_factor: f64,
    pub kinetic_energy_flux_factor: f64,
    /// Signed source-density momentum contribution beta*G*(u_out-u_in) [Pa].
    pub acceleration_pressure_pa: f64,
    /// sum(dp_loss) + acceleration - drive [Pa].
    pub pressure_residual_pa: f64,
    /// sum(dp_loss) + max_cell |beta*G*(u_cell-u_in)| [Pa].
    /// Bounds mechanical pressure departures without cancellation.
    pub common_pressure_bound_pa: f64,
    /// Signed total outlet-minus-inlet kinetic-energy flow [W].
    pub kinetic_power_change_w: f64,
    /// Sum of absolute regional kinetic-energy flow changes [W].
    pub neglected_kinetic_power_bound_w: f64,
    /// Full native material-transport solves, including bracket endpoints.
    pub flow_evaluations: usize,
    /// Actual native thermal substeps across every pressure/feedback trial.
    pub substeps: usize,
}

/// Refusal never publishes a partial gas field or changes the supplied problem.
#[derive(Debug)]
pub enum GasPressureError {
    Invalid(&'static str),
    Transport(MaterialSteadyGasError),
    NotBracketed {
        lower_residual_pa: f64,
        upper_residual_pa: f64,
    },
    NotConverged {
        evaluations: usize,
        last_residual_pa: f64,
    },
    SubstepBudget {
        performed: usize,
        maximum: usize,
    },
    PressureVariation {
        bound_pa: f64,
        allowance_pa: f64,
    },
    NeglectedHydraulicPower {
        power_w: f64,
        allowance_w: f64,
    },
    NeglectedKineticPower {
        bound_w: f64,
        allowance_w: f64,
    },
}
impl core::fmt::Display for GasPressureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pressure-driven material gas: {self:?}")
    }
}
impl std::error::Error for GasPressureError {}
impl From<MaterialSteadyGasError> for GasPressureError {
    fn from(error: MaterialSteadyGasError) -> Self {
        Self::Transport(error)
    }
}

/// Select positive through-flow from an imposed static-pressure difference.
///
/// Geometry must have one constant flow area and no side injections; every
/// trial recomputes local density, viscosity, convection and enthalpy transport.
/// Friction and thermal acceleration both enter momentum closure with the same
/// declared, constant cross-sectional velocity-profile shape. All original
/// caloric/source, Reynolds and low-Mach limits remain active. The source EOS
/// stays at its common pressure: this is a bounded low-pressure-drop model,
/// not a compressible pressure field, choke/valve model or full energy solve.
///
/// The caller's bracket selects a candidate; heating can make the response
/// nonmonotone, so no uniqueness or dynamical-stability claim is made. A failed
/// trial is propagated, never skipped or replaced by a different source law.
/// The template inlet flow is replaced; its thermal guesses remain unchanged.
/// Omitted hydraulic and regional kinetic powers must fit their own allowances.
///
/// # Errors
/// Invalid geometry/controls, unsupported trials, cancellation, exhausted shared
/// work, excessive pressure/power, missing sign bracket or pressure nonconvergence.
pub fn solve_pressure_driven_material_gas(
    cx: &Cx<'_>,
    problem: &SteadyAxialGas,
    transport: &[GasCellTransport],
    drive: Pressure,
    controls: SteadyGasSolve,
    feedback: GasTransportFeedback,
    budget: GasPressureBudget,
) -> Result<PressureDrivenGasSolution, GasPressureError> {
    check(cx)?;
    let pressure = budget.pressure;
    let (mut low, mut high) = pressure.mass_flow_bracket_kg_per_s;
    let nonnegative = |value: f64| value.is_finite() && value >= 0.0;
    let beta_squared = budget.momentum_flux_factor * budget.momentum_flux_factor;
    let alpha = budget.kinetic_energy_flux_factor;
    if !positive(drive.value())
        || !positive(low)
        || !positive(high)
        || low >= high
        || !positive(pressure.pressure_absolute_tolerance_pa)
        || !(0.0..1.0).contains(&pressure.pressure_relative_tolerance)
        || !(2..=128).contains(&pressure.max_flow_evaluations)
        || !nonnegative(pressure.maximum_neglected_hydraulic_power_w)
        || !nonnegative(budget.maximum_neglected_kinetic_power_w)
        || !beta_squared.is_finite()
        || budget.momentum_flux_factor < 1.0
        || !alpha.is_finite()
        || alpha < 1.0
        || (alpha < beta_squared && beta_squared - alpha > 64.0 * f64::EPSILON * beta_squared)
        || budget.maximum_substeps == 0
        || feedback.maximum_substeps == 0
        || feedback.maximum_substeps > 1_000_000
        || !positive(feedback.maximum_pressure_drop_fraction)
        || feedback.maximum_pressure_drop_fraction > 0.1
        || !positive(problem.static_pressure_pa)
    {
        return Err(GasPressureError::Invalid(
            "invalid pressure bracket, velocity profile, work or omission budget",
        ));
    }
    if !problem.injections.is_empty()
        || transport.is_empty()
        || transport.len() > MAX_STEADY_GAS_CELLS
        || transport.len() != problem.cells.len()
    {
        return Err(GasPressureError::Invalid(
            "pressure selection needs a nonempty axial path without side injections",
        ));
    }
    let area = transport[0].resistance.flow_area().value();
    if transport.iter().any(|law| {
        (law.resistance.flow_area().value() - area).abs()
            > 64.0 * f64::EPSILON * area.max(law.resistance.flow_area().value())
    }) {
        return Err(GasPressureError::Invalid(
            "thermal acceleration requires one constant flow area",
        ));
    }
    let pressure_allowance = problem.static_pressure_pa * feedback.maximum_pressure_drop_fraction;
    if !positive(pressure_allowance) || drive.value() > pressure_allowance {
        return Err(GasPressureError::PressureVariation {
            bound_pa: drive.value(),
            allowance_pa: pressure_allowance,
        });
    }
    let mut trial = problem.clone();
    let mut evaluations = 0;
    let mut total_substeps = 0;
    let mut probe = |mass: f64| -> Result<PressureDrivenGasSolution, GasPressureError> {
        check(cx)?;
        if total_substeps >= budget.maximum_substeps {
            return Err(GasPressureError::SubstepBudget {
                performed: total_substeps,
                maximum: budget.maximum_substeps,
            });
        }
        trial.inlet_mass_flow_kg_s = mass;
        let gas = solve_material_gas_transport(
            cx,
            &trial,
            transport,
            controls,
            GasTransportFeedback {
                maximum_substeps: feedback
                    .maximum_substeps
                    .min(budget.maximum_substeps - total_substeps),
                ..feedback
            },
        )?;
        evaluations += 1;
        total_substeps += gas.substeps;
        let inlet = gas
            .thermal
            .field()
            .state_at(0.0)
            .map_err(MaterialSteadyGasError::from)?;
        let inlet_speed = (mass / inlet.state().density) / area;
        let mass_flux = mass / area;
        let mut previous_speed = inlet_speed;
        let mut acceleration_pressure_pa = 0.0;
        let mut maximum_acceleration: f64 = 0.0;
        let mut kinetic_power_change_w = 0.0;
        let mut neglected_kinetic_power_bound_w = 0.0;
        for cell in &gas.cells {
            check(cx)?;
            let speed = (mass / cell.gas.state().density) / area;
            acceleration_pressure_pa =
                budget.momentum_flux_factor * mass_flux * (speed - inlet_speed);
            maximum_acceleration = maximum_acceleration.max(acceleration_pressure_pa.abs());
            let kinetic =
                (0.5 * mass) * (speed - previous_speed) * (speed + previous_speed) * alpha;
            kinetic_power_change_w += kinetic;
            neglected_kinetic_power_bound_w += kinetic.abs();
            previous_speed = speed;
        }
        let common_pressure_bound_pa = gas.total_pressure_drop_pa + maximum_acceleration;
        let pressure_residual_pa =
            gas.total_pressure_drop_pa + acceleration_pressure_pa - drive.value();
        if ![
            inlet_speed,
            mass_flux,
            acceleration_pressure_pa,
            common_pressure_bound_pa,
            pressure_residual_pa,
            kinetic_power_change_w,
            neglected_kinetic_power_bound_w,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(GasPressureError::Invalid(
                "gas momentum or kinetic-energy response is not representable",
            ));
        }
        if common_pressure_bound_pa > pressure_allowance {
            return Err(GasPressureError::PressureVariation {
                bound_pa: common_pressure_bound_pa,
                allowance_pa: pressure_allowance,
            });
        }
        if gas.hydraulic_dissipation_w > pressure.maximum_neglected_hydraulic_power_w {
            return Err(GasPressureError::NeglectedHydraulicPower {
                power_w: gas.hydraulic_dissipation_w,
                allowance_w: pressure.maximum_neglected_hydraulic_power_w,
            });
        }
        if neglected_kinetic_power_bound_w > budget.maximum_neglected_kinetic_power_w {
            return Err(GasPressureError::NeglectedKineticPower {
                bound_w: neglected_kinetic_power_bound_w,
                allowance_w: budget.maximum_neglected_kinetic_power_w,
            });
        }
        check(cx)?;
        Ok(PressureDrivenGasSolution {
            gas,
            mass_flow_kg_per_s: mass,
            driving_pressure_pa: drive.value(),
            momentum_flux_factor: budget.momentum_flux_factor,
            kinetic_energy_flux_factor: alpha,
            acceleration_pressure_pa,
            pressure_residual_pa,
            common_pressure_bound_pa,
            kinetic_power_change_w,
            neglected_kinetic_power_bound_w,
            flow_evaluations: evaluations,
            substeps: total_substeps,
        })
    };
    let closed = |value: &PressureDrivenGasSolution| {
        let predicted = value.gas.total_pressure_drop_pa + value.acceleration_pressure_pa;
        value.pressure_residual_pa.abs()
            <= pressure.pressure_absolute_tolerance_pa.max(
                pressure.pressure_relative_tolerance * drive.value().abs().max(predicted.abs()),
            )
    };
    let a = probe(low)?;
    if closed(&a) {
        return Ok(a);
    }
    let mut low_residual = a.pressure_residual_pa;
    drop(a);
    let b = probe(high)?;
    if closed(&b) {
        return Ok(b);
    }
    let mut last_residual = b.pressure_residual_pa;
    if low_residual.is_sign_positive() == last_residual.is_sign_positive() {
        return Err(GasPressureError::NotBracketed {
            lower_residual_pa: low_residual,
            upper_residual_pa: last_residual,
        });
    }
    drop(b);
    let mut performed = 2;
    for _ in 2..pressure.max_flow_evaluations {
        let mass = 0.5 * low + 0.5 * high;
        if mass == low || mass == high {
            break;
        }
        let value = probe(mass)?;
        performed = value.flow_evaluations;
        if closed(&value) {
            return Ok(value);
        }
        last_residual = value.pressure_residual_pa;
        if last_residual.is_sign_positive() == low_residual.is_sign_positive() {
            low = mass;
            low_residual = last_residual;
        } else {
            high = mass;
        }
    }
    Err(GasPressureError::NotConverged {
        evaluations: performed,
        last_residual_pa: last_residual,
    })
}
