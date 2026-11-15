//! `friction_interface`: a dry-friction law of `fs-tribo` between two scene
//! ends, on the kernel's load line.
//!
//! The law is one rung of [`FrictionLaw`]: Coulomb, linear in speed,
//! Stribeck, hyperbolic or shifted-power weakening, authored in the scene,
//! or the measured rosin fit `fs-tribo` cites by name. The component is
//! [`FrictionInterface`], which owns the stick, slip, arrest and break-away
//! rows, the stop location and the heat; nothing here evaluates friction or
//! integrates anything.

use fs_ports::component::{Component, PortSpec};
use fs_system::{NodeId, SystemBuilder};
use fs_tribo::tangential::component::{FrictionInterface, FrictionState, InterfaceLoad};
use fs_tribo::{DryTangentialState, FrictionLaw, TangentialSolveConfig};

use super::{J, Refusal, fields, number, number_or};

/// What the kernel holds for a friction interface.
pub(super) type State = FrictionState<DryTangentialState>;

/// The measured rosin fit, named as `fs-tribo` names it.
const ROSIN: &str = "rosin_woodhouse_galluzzo_2025_eq13";

/// Every rung offered here is independent of temperature, so the state an
/// interface starts from is the same at any positive temperature. The
/// interface temperature of each step reaches the law through `heat`.
const START_K: f64 = 1.0;

/// Default budgets of the law's root search and of the interface, as the
/// tests of `fs-tribo` and `fs-couple` declare them.
const MAX_ITERATIONS: f64 = 64.0;
const FORCE_ABSOLUTE_TOLERANCE_N: f64 = 1.0e-12;
const FORCE_RELATIVE_TOLERANCE: f64 = 1.0e-12;
const IMPULSE_TOLERANCE_N_S: f64 = 1.0e-10;
const WORK_TOLERANCE_J: f64 = 1.0e-12;
/// The default ceiling of the traction search, in static capacities.
const FORCE_CEILING_IN_STATIC_CAPACITIES: f64 = 1.0e3;

/// One rung of the dry-friction ladder: the rosin fit by name, or an object
/// with one key naming the rung and holding its coefficients.
fn law(value: &J, what: &str) -> Result<FrictionLaw, String> {
    let rungs = "coulomb, velocity_dependent, stribeck, hyperbolic or shifted_power";
    if let Some(name) = value.as_str() {
        return if name == ROSIN {
            Ok(FrictionLaw::rosin_woodhouse_galluzzo_2025_eq13())
        } else {
            Err(format!("{what}: the named law is {ROSIN:?}, not {name:?}"))
        };
    }
    let Some([(rung, parameters)]) = value.as_object() else {
        return Err(format!(
            "{what}: a law is {ROSIN:?} or an object with one key, {rungs}"
        ));
    };
    let inner = format!("{what}.{rung}");
    let law = match rung.as_str() {
        "coulomb" => {
            fields(parameters, &["static_mu", "kinetic_mu"], &inner)?;
            FrictionLaw::Coulomb {
                static_mu: number(parameters, "static_mu", &inner)?,
                kinetic_mu: number(parameters, "kinetic_mu", &inner)?,
            }
        }
        "velocity_dependent" => {
            fields(
                parameters,
                &["static_mu", "mu_zero", "slope_per_speed_s_m"],
                &inner,
            )?;
            FrictionLaw::VelocityDependent {
                static_mu: number(parameters, "static_mu", &inner)?,
                mu_zero: number(parameters, "mu_zero", &inner)?,
                slope_per_speed: number(parameters, "slope_per_speed_s_m", &inner)?,
            }
        }
        "stribeck" => {
            fields(
                parameters,
                &[
                    "static_mu",
                    "kinetic_mu",
                    "characteristic_speed_m_s",
                    "viscous_per_speed_s_m",
                ],
                &inner,
            )?;
            FrictionLaw::Stribeck {
                static_mu: number(parameters, "static_mu", &inner)?,
                kinetic_mu: number(parameters, "kinetic_mu", &inner)?,
                characteristic_speed: number(parameters, "characteristic_speed_m_s", &inner)?,
                viscous_per_speed: number_or(parameters, "viscous_per_speed_s_m", 0.0, &inner)?,
            }
        }
        "hyperbolic" => {
            fields(
                parameters,
                &["static_mu", "kinetic_mu", "characteristic_speed_m_s"],
                &inner,
            )?;
            FrictionLaw::Hyperbolic {
                static_mu: number(parameters, "static_mu", &inner)?,
                kinetic_mu: number(parameters, "kinetic_mu", &inner)?,
                characteristic_speed: number(parameters, "characteristic_speed_m_s", &inner)?,
            }
        }
        "shifted_power" => {
            fields(
                parameters,
                &[
                    "static_mu",
                    "kinetic_mu",
                    "characteristic_speed_m_s",
                    "exponent",
                ],
                &inner,
            )?;
            FrictionLaw::ShiftedPower {
                static_mu: number(parameters, "static_mu", &inner)?,
                kinetic_mu: number(parameters, "kinetic_mu", &inner)?,
                characteristic_speed: number(parameters, "characteristic_speed_m_s", &inner)?,
                exponent: number(parameters, "exponent", &inner)?,
            }
        }
        other => {
            return Err(format!(
                "{what}: unknown law {other:?}; one of {ROSIN:?}, {rungs}"
            ));
        }
    };
    law.validate()
        .map_err(|error| format!("{inner}: {error}"))?;
    Ok(law)
}

/// The budgets of the law's root search and of the interface, each with
/// its default when the scene leaves it out.
fn load(
    definition: &J,
    normal_n: f64,
    law: &FrictionLaw,
    what: &str,
) -> Result<InterfaceLoad, String> {
    let inner = format!("{what}.solve");
    let solve = definition.get("solve").unwrap_or(&J::Null);
    if !matches!(solve, J::Null) {
        fields(
            solve,
            &[
                "max_iterations",
                "maximum_force_n",
                "force_absolute_tolerance_n",
                "force_relative_tolerance",
                "impulse_tolerance_n_s",
                "work_tolerance_j",
            ],
            &inner,
        )?;
    }
    let static_capacity_n = law
        .static_coefficient()
        .map_err(|error| format!("{what}: {error}"))?
        * normal_n;
    let iterations = number_or(solve, "max_iterations", MAX_ITERATIONS, &inner)?;
    if !(1.0..=128.0).contains(&iterations) || iterations.fract() != 0.0 {
        return Err(format!(
            "{inner}: max_iterations is a whole number from 1 to 128"
        ));
    }
    Ok(InterfaceLoad {
        normal_n,
        nominal_pressure_pa: None,
        solve: TangentialSolveConfig {
            max_iterations: iterations as usize,
            maximum_force_n: number_or(
                solve,
                "maximum_force_n",
                (FORCE_CEILING_IN_STATIC_CAPACITIES * static_capacity_n).max(1.0),
                &inner,
            )?,
            force_absolute_tolerance_n: number_or(
                solve,
                "force_absolute_tolerance_n",
                FORCE_ABSOLUTE_TOLERANCE_N,
                &inner,
            )?,
            force_relative_tolerance: number_or(
                solve,
                "force_relative_tolerance",
                FORCE_RELATIVE_TOLERANCE,
                &inner,
            )?,
        },
        impulse_tolerance_n_s: number_or(
            solve,
            "impulse_tolerance_n_s",
            IMPULSE_TOLERANCE_N_S,
            &inner,
        )?,
        work_tolerance_j: number_or(solve, "work_tolerance_j", WORK_TOLERANCE_J, &inner)?,
    })
}

/// A `friction_interface` element: the interface on the kernel's load line,
/// starting at the slip the scene declares in the movers' sign.
pub(super) fn interface(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(
        definition,
        &["kind", "normal_n", "law", "slip_m_s", "solve"],
        what,
    )?;
    let normal_n = number(definition, "normal_n", what)?;
    let law = law(
        definition
            .get("law")
            .ok_or_else(|| format!("{what}: missing law"))?,
        &format!("{what}.law"),
    )?;
    let load = load(definition, normal_n, &law, what)?;
    let slip = number_or(definition, "slip_m_s", 0.0, what)?;
    if !slip.is_finite() {
        return Err(format!("{what}.slip_m_s must be finite").into());
    }
    let component = FrictionInterface::new(move |_kelvin: f64| Ok(law.clone()), load)
        .map_err(|error| format!("{what}: {error}"))?;
    // Scene velocities use the mover's sign. Each end receives minus its
    // mover's flow, so the owner's slip has the opposite sign.
    let state = component
        .start(START_K, -slip)
        .map_err(|error| format!("{what}: {error}"))?;
    let ports = component.ports().to_vec();
    let node = builder
        .add_on_line(name, component, state)
        .map_err(|error| Refusal::Admission(error.to_string()))?;
    Ok((node, ports))
}

/// The interface's slip in the movers' sign, and whether it is stuck.
pub(super) fn values(state: Option<&State>) -> Vec<(&'static str, f64)> {
    // The owner's slip is flow at a less flow at b in the ports' signs.
    vec![
        (
            "slip_m_s",
            state.map_or(f64::NAN, |state| -state.slip_m_s()),
        ),
        (
            "stuck",
            state.map_or(f64::NAN, |state| if state.at_rest() { 1.0 } else { 0.0 }),
        ),
    ]
}
