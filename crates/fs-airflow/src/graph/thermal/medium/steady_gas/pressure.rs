//! Pressure-selected gas flow with source-density thermal acceleration.
//!
//! Every trial reuses the native stationary species, NASA-9 enthalpy, convection
//! and hydraulic solve. Only the inlet mass rate changes. A constant-area path
//! without side injections has dp = sum(dp_loss) + beta*(mdot/A)*(u_out - u_in).
//! Thermodynamic pressure remains the source's common reference; the mechanical
//! pressure envelope and omitted mechanical-energy powers are explicitly bounded.
//! Opt-in gravity reuses the same source-density head and conservative work owner.

/// Source-consistent flow splitting and enthalpy mixing between pressure ports.
pub mod parallel;

use fs_exec::Cx;
use fs_material::gas_field::steady::{MAX_STEADY_GAS_CELLS, SteadyAxialGas, SteadyGasSolve};
use fs_qty::Pressure;

use crate::graph::thermal::medium::series::hydraulic::pressure::{
    HydrostaticPressure, SeriesGravity, SeriesPressureError,
};

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

/// Original gas-flow controls plus the explicit omission bound for gravitational energy.
#[derive(Clone, Copy, Debug)]
pub struct GasGravityBudget {
    pub flow: GasPressureBudget,
    /// Sum of absolute regional potential-energy flow changes [W], not just
    /// the net endpoint change. Finite and nonnegative; no energy becomes heat.
    pub maximum_neglected_potential_power_w: f64,
}

/// Original geometry, source-density head and separate conservative gravity work.
#[derive(Clone, Debug, PartialEq)]
pub struct GasGravityState {
    pub geometry: SeriesGravity,
    pub hydrostatic: HydrostaticPressure,
    /// -mdot * g·(outlet-inlet) [W]; exactly zero on a closed geometric path.
    pub potential_power_change_w: f64,
    /// Sum of absolute regional potential-energy flow changes [W].
    pub neglected_potential_power_bound_w: f64,
    /// Sum of absolute regional hydrostatic pressure gains [Pa].
    pub hydrostatic_pressure_bound_pa: f64,
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
    /// sum(dp_loss) + acceleration - drive - hydrostatic gain [Pa].
    pub pressure_residual_pa: f64,
    /// sum(dp_loss) + max_cell |beta*G*(u_cell-u_in)|, plus any absolute
    /// regional hydrostatic gains [Pa]. No head cancellation hides a departure.
    pub common_pressure_bound_pa: f64,
    /// Signed total outlet-minus-inlet kinetic-energy flow [W].
    pub kinetic_power_change_w: f64,
    /// Sum of absolute regional kinetic-energy flow changes [W].
    pub neglected_kinetic_power_bound_w: f64,
    /// Source-density gravity response when explicitly requested.
    pub gravity: Option<GasGravityState>,
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
    Gravity(SeriesPressureError),
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
    NeglectedPotentialPower {
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
    solve_pressure_inner(
        cx, problem, transport, drive, controls, feedback, budget, None,
    )
}

/// Select the same native gas flow with explicit source-density hydrostatic head.
///
/// The imposed pressure may be signed or zero. Each connected gravity segment
/// maps to one native cell, using its accepted well-mixed bulk density. Its
/// endpoint displacement cannot exceed that cell's physical axial path length.
/// Geometry and uniform acceleration are authored; no ambient density or axis
/// is inferred. For a chimney, supply the outside column's pressure difference
/// explicitly, for example from original material states and `SeriesGravity`.
///
/// Pressure closure is loss + acceleration = imposed pressure + sum(rho*g·dr).
/// Gravity adds no heat: the native enthalpy solve omits potential-energy flow,
/// whose absolute regional changes must fit the declared watt allowance.
/// A closed geometric path has zero net gravity work but can have buoyant head;
/// it retains the original prescribed inlet reservoir, not a periodic thermal
/// closure or a closed-inventory thermosiphon. Refine cells to resolve density.
///
/// # Errors
/// Original flow/source refusals, incompatible geometry, excessive pressure or
/// omitted potential power, missing positive-flow bracket or cancellation.
pub fn solve_pressure_gravity_driven_material_gas(
    cx: &Cx<'_>,
    problem: &SteadyAxialGas,
    transport: &[GasCellTransport],
    drive: Pressure,
    gravity: &SeriesGravity,
    controls: SteadyGasSolve,
    feedback: GasTransportFeedback,
    budget: GasGravityBudget,
) -> Result<PressureDrivenGasSolution, GasPressureError> {
    solve_pressure_inner(
        cx,
        problem,
        transport,
        drive,
        controls,
        feedback,
        budget.flow,
        Some((gravity, budget.maximum_neglected_potential_power_w)),
    )
}

fn solve_pressure_inner(
    cx: &Cx<'_>,
    problem: &SteadyAxialGas,
    transport: &[GasCellTransport],
    drive: Pressure,
    controls: SteadyGasSolve,
    feedback: GasTransportFeedback,
    budget: GasPressureBudget,
    gravity: Option<(&SeriesGravity, f64)>,
) -> Result<PressureDrivenGasSolution, GasPressureError> {
    check(cx)?;
    let pressure = budget.pressure;
    let (mut low, mut high) = pressure.mass_flow_bracket_kg_per_s;
    let nonnegative = |value: f64| value.is_finite() && value >= 0.0;
    let beta_squared = budget.momentum_flux_factor * budget.momentum_flux_factor;
    let alpha = budget.kinetic_energy_flux_factor;
    if !drive.value().is_finite()
        || (gravity.is_none() && !positive(drive.value()))
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
    if let Some((path, allowance)) = gravity {
        if path.section_count() != problem.cells.len() || !nonnegative(allowance) {
            return Err(GasPressureError::Invalid(
                "one gravity segment per cell and finite nonnegative potential-power allowance required",
            ));
        }
        let mut previous = 0.0;
        for (cell, pair) in problem.cells.iter().zip(path.vertices_m().windows(2)) {
            check(cx)?;
            let length = cell.end_m - previous;
            let chord_squared: f64 = (0..3)
                .map(|axis| {
                    let ratio = (pair[1][axis] - pair[0][axis]) / length;
                    ratio * ratio
                })
                .sum();
            if !positive(length)
                || !chord_squared.is_finite()
                || chord_squared > 1.0 + 128.0 * f64::EPSILON
            {
                return Err(GasPressureError::Invalid(
                    "gravity endpoint displacement exceeds the native axial cell length",
                ));
            }
            previous = cell.end_m;
        }
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
    if !positive(pressure_allowance) || drive.value().abs() > pressure_allowance {
        return Err(GasPressureError::PressureVariation {
            bound_pa: drive.value().abs(),
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
        let gravity_state = gravity
            .map(|(path, allowance)| observe_gravity(cx, path, &gas, mass, allowance))
            .transpose()?;
        let mut common_pressure_bound_pa = gas.total_pressure_drop_pa + maximum_acceleration;
        let mut effective_drive = drive.value();
        if let Some(state) = &gravity_state {
            common_pressure_bound_pa += state.hydrostatic_pressure_bound_pa;
            effective_drive += state.hydrostatic.pressure_gain_pa();
        }
        let pressure_residual_pa =
            gas.total_pressure_drop_pa + acceleration_pressure_pa - effective_drive;
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
            gravity: gravity_state,
            flow_evaluations: evaluations,
            substeps: total_substeps,
        })
    };
    let closed = |value: &PressureDrivenGasSolution| {
        let predicted = value.gas.total_pressure_drop_pa + value.acceleration_pressure_pa;
        let effective_drive = drive.value()
            + value
                .gravity
                .as_ref()
                .map_or(0.0, |state| state.hydrostatic.pressure_gain_pa());
        value.pressure_residual_pa.abs()
            <= pressure.pressure_absolute_tolerance_pa.max(
                pressure.pressure_relative_tolerance * effective_drive.abs().max(predicted.abs()),
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

fn observe_gravity(
    cx: &Cx<'_>,
    path: &SeriesGravity,
    gas: &MaterialSteadyGasSolution,
    mass: f64,
    maximum_neglected_potential_power_w: f64,
) -> Result<GasGravityState, GasPressureError> {
    let hydrostatic = path
        .evaluate_densities(
            gas.cells.iter().map(|cell| cell.gas.state().density),
            mass,
            || cx.checkpoint().map_err(|_| SeriesPressureError::Cancelled),
        )
        .map_err(GasPressureError::Gravity)?;
    let potential_power_change_w = -hydrostatic.gravitational_power_w();
    let mut hydrostatic_pressure_bound_pa = 0.0;
    let mut neglected_potential_power_bound_w = 0.0;
    for (&head, &work) in hydrostatic
        .section_pressure_gains_pa()
        .iter()
        .zip(path.section_specific_work_j_kg())
    {
        check(cx)?;
        hydrostatic_pressure_bound_pa += head.abs();
        neglected_potential_power_bound_w += (mass * work).abs();
    }
    if !hydrostatic_pressure_bound_pa.is_finite() || !neglected_potential_power_bound_w.is_finite()
    {
        return Err(GasPressureError::Invalid(
            "regional gravity pressure or potential power overflows",
        ));
    }
    // Endpoint work and the stable total head can differ from rounded section
    // sums by an ulp; the reported bounds must also cover those retained totals.
    hydrostatic_pressure_bound_pa =
        hydrostatic_pressure_bound_pa.max(hydrostatic.pressure_gain_pa().abs());
    neglected_potential_power_bound_w =
        neglected_potential_power_bound_w.max(potential_power_change_w.abs());
    if neglected_potential_power_bound_w > maximum_neglected_potential_power_w {
        return Err(GasPressureError::NeglectedPotentialPower {
            bound_w: neglected_potential_power_bound_w,
            allowance_w: maximum_neglected_potential_power_w,
        });
    }
    check(cx)?;
    Ok(GasGravityState {
        geometry: path.clone(),
        hydrostatic,
        potential_power_change_w,
        neglected_potential_power_bound_w,
        hydrostatic_pressure_bound_pa,
    })
}
