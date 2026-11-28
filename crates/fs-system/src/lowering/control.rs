//! Controllers from the description: a declared PID law, lowered onto the
//! kernel's held [`SampledController`]. It reads a signal a lumped body
//! publishes and commands an ideal effort source, or a heater, on the
//! actuator's port.

use fs_ports::PortKind;
use fs_ports::component::{Carrier, PortRole};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{ControllerDecl, ControllerLaw};
use fs_scenario::entity::{EntityCatalog, EntityId};

use super::{
    FREQUENCY_DIMS, LoweringError, LoweringGraph, POWER_DIMS, PartLoweringContext,
    ProvenanceRecorder, SourceLoweringReceipt, TEMPERATURE_DIMS, TIME_DIMS,
};
use crate::elements::{
    ControlledSource, ControllerState, HeatCapacity, Heater, Mass, SampledController,
    SignalLimiter, SourceAccount,
};
use crate::{NodeId, PortRef};

/// The signals a lumped rigid body publishes when a controller reads it:
/// position and velocity, or angle and angular velocity for a rotation.
const BODY_SIGNALS: [&str; 2] = ["position", "velocity"];

/// The signal a lumped thermal body publishes when a controller reads it.
const THERMAL_SIGNAL: &str = "temperature";

/// Whether a declared controller reads one of `signals` on the part.
fn sensed(ctx: &PartLoweringContext<'_>, signals: &[&str]) -> bool {
    let input = ctx.draft.input();
    let mut sensors = input.controllers.iter().flat_map(|c| &c.sensor_ports);
    sensors.any(|(e, port)| *e == ctx.part.entity && signals.contains(&port.as_str()))
}

/// The part's mass, publishing its signals when a declared controller reads
/// them, and unchanged otherwise.
pub(super) fn publish_if_sensed(ctx: &PartLoweringContext<'_>, mass: Mass) -> Mass {
    if sensed(ctx, &BODY_SIGNALS) {
        mass.publishing()
    } else {
        mass
    }
}

/// The part's heat capacity, publishing its temperature when a declared
/// controller reads it, and unchanged otherwise.
pub(super) fn publish_heat_if_sensed(
    ctx: &PartLoweringContext<'_>,
    cap: HeatCapacity,
) -> HeatCapacity {
    if sensed(ctx, &[THERMAL_SIGNAL]) {
        cap.publishing()
    } else {
        cap
    }
}

/// Expose on `entity` the signals `node` publishes, if it publishes any.
pub(super) fn expose_signals(graph: &mut LoweringGraph, entity: EntityId, node: NodeId) {
    for name in BODY_SIGNALS.into_iter().chain([THERMAL_SIGNAL]) {
        if let Ok(port) = graph.port(node, name) {
            graph.expose_port(entity, name, port);
        }
    }
}

/// `num / den * s^time`, exponent by exponent.
fn ratio(num: Dims, den: Dims, time: i8) -> Dims {
    Dims(core::array::from_fn(|i| {
        num.0[i] - den.0[i] + time * TIME_DIMS.0[i]
    }))
}

fn refuse(scope: &str, detail: String) -> LoweringError {
    LoweringError::Invalid {
        item: scope.to_owned(),
        detail,
    }
}

fn missing(scope: &str, entity: EntityId, port: &str) -> LoweringError {
    LoweringError::MissingPort {
        item: scope.to_owned(),
        entity,
        port: port.to_owned(),
    }
}

/// The PID law's `(kp, ki, kd, setpoint)`; any other law refuses by name.
fn pid(ctrl: &ControllerDecl, scope: &str) -> Result<[QtyAny; 4], LoweringError> {
    let other = match &ctrl.law {
        ControllerLaw::Pid {
            kp,
            ki,
            kd,
            setpoint,
        } => return Ok([*kp, *ki, *kd, *setpoint]),
        ControllerLaw::LookupTable { .. } => "lookup-table",
        ControllerLaw::CardBacked { .. } => "card-backed",
    };
    Err(refuse(
        scope,
        format!("a {other} law does not lower yet; declare a PID law"),
    ))
}

/// The one signal a PID law reads, and its dimensions: a lumped body's
/// position or velocity in its power port's coordinates, or its temperature.
fn sensor_signal(
    ctrl: &ControllerDecl,
    scope: &str,
    catalog: &EntityCatalog,
    graph: &LoweringGraph,
) -> Result<(PortRef, Dims), LoweringError> {
    let [(sensed, sensor)] = ctrl.sensor_ports.as_slice() else {
        let n = ctrl.sensor_ports.len();
        return Err(refuse(
            scope,
            format!("a PID law reads one sensor, not {n}"),
        ));
    };
    if !BODY_SIGNALS.contains(&sensor.as_str()) && sensor != THERMAL_SIGNAL {
        return Err(refuse(
            scope,
            format!(
                "sensor '{sensor}': a lumped body publishes position, velocity and \
                 temperature, and no other sensor signal lowers yet"
            ),
        ));
    }
    let port = graph
        .resolve_entity_port(catalog, *sensed, sensor)
        .ok_or_else(|| missing(scope, *sensed, sensor))?;
    let body = PortRef {
        node: port.node,
        port: 0,
    };
    let dims = match graph.port_meta.get(&body) {
        Some(&(Carrier::Power(kind), ..)) if sensor == "position" => {
            ratio(kind.canonical_flow_dimensions(), Dims::NONE, 1)
        }
        Some(&(Carrier::Power(kind), ..)) => kind.canonical_flow_dimensions(),
        Some(&(Carrier::Heat, ..)) => TEMPERATURE_DIMS,
        _ => {
            let why = format!("sensor '{sensor}' is not on a lumped body");
            return Err(refuse(scope, why));
        }
    };
    Ok((port, dims))
}

/// What a controller drives: an ideal effort source in a power domain, or a
/// heater on a thermal node.
#[derive(Clone, Copy)]
enum Drive {
    Effort(PortKind),
    Heat,
}

/// The actuator's port, one a commanded source can drive, and the drive.
fn actuator(
    ctrl: &ControllerDecl,
    scope: &str,
    catalog: &EntityCatalog,
    graph: &LoweringGraph,
) -> Result<(PortRef, Drive), LoweringError> {
    let (entity, port) = (ctrl.actuator_port.0, ctrl.actuator_port.1.as_str());
    let target = graph
        .resolve_entity_port(catalog, entity, port)
        .ok_or_else(|| missing(scope, entity, port))?;
    match graph.port_meta.get(&target) {
        Some(&(Carrier::Power(kind), PortRole::EffortIn, 1)) => Ok((target, Drive::Effort(kind))),
        Some(&(Carrier::Heat, PortRole::FlowIn, _)) => Ok((target, Drive::Heat)),
        _ => Err(refuse(
            scope,
            format!(
                "actuator '{port}' takes no commanded effort or heat; name a lumped body's \
                 power or thermal port"
            ),
        )),
    }
}

/// Lower one declared controller: a PID law on one published sensor signal,
/// sampled at its declared rate and held between samples, commanding an
/// ideal effort source on the actuator's power port, or a heater behind a
/// limiter at 0 on its thermal port. Its receipt row joins the sources as
/// `controller:<name>`.
///
/// At sample `k` of period `t`, with `e[-1] = 0` and `S[-1] = 0`:
/// `e[k] = setpoint - y[k]`, `S[k] = S[k-1] + e[k]` and
/// `u[k] = kp e[k] + ki t S[k] + kd (e[k] - e[k-1]) / t`.
pub(super) fn lower_controller(
    ctrl: &ControllerDecl,
    catalog: &EntityCatalog,
    graph: &mut LoweringGraph,
) -> Result<SourceLoweringReceipt, LoweringError> {
    let scope = format!("controller:{}", ctrl.name);
    let [kp, ki, kd, setpoint] = pid(ctrl, &scope)?;
    let (sensor_port, sensor) = sensor_signal(ctrl, &scope, catalog, graph)?;
    let (target, drive) = actuator(ctrl, &scope, catalog, graph)?;
    let effort = match drive {
        Drive::Effort(kind) => kind.canonical_effort_dimensions(),
        Drive::Heat => POWER_DIMS,
    };

    let mut recorder = ProvenanceRecorder::default();
    for (label, q, want) in [
        ("kp", kp, ratio(effort, sensor, 0)),
        ("ki", ki, ratio(effort, sensor, -1)),
        ("kd", kd, ratio(effort, sensor, 1)),
        ("setpoint", setpoint, sensor),
    ] {
        if q.dims != want {
            let (got, need) = (q.dims.0, want.0);
            let why =
                format!("{label} has dimensions {got:?}; this sensor and actuator need {need:?}");
            return Err(refuse(&scope, why));
        }
        recorder.record_declared_field(label, q.value, q.dims, &scope, format!("law.pid.{label}"));
    }
    let rate = ctrl.sample_rate.value;
    recorder.record_declared_field("sample_rate", rate, FREQUENCY_DIMS, &scope, "sample_rate");

    // State x = (b, S[k-1], e[k-1]) with a bias b = 1, so the setpoint enters
    // through A and C: e[k] = er b + eu y[k].
    let (r, t) = (setpoint.value, 1.0 / rate);
    let (er, eu) = (r, -1.0);
    let (gi, gd) = (ki.value * t, kd.value / t);
    let gain = kp.value + gi + gd;
    let element = |source| LoweringError::Element {
        item: scope.clone(),
        source,
    };
    let law = SampledController::new(
        t,
        3,
        1,
        1,
        vec![1.0, 0.0, 0.0, er, 1.0, 0.0, er, 0.0, 0.0],
        vec![0.0, eu, eu],
        vec![gain * er, gi, -gd],
        vec![gain * eu],
    )
    .map_err(element)?;
    let memory = ControllerState::new(&law, vec![1.0, 0.0, 0.0]).map_err(element)?;

    let law_name = format!("controller.{}", ctrl.name);
    let law = graph.add_with_tangent(&law_name, None, law, memory)?;
    let (input, output) = (graph.port(law, "input0")?, graph.port(law, "output")?);
    graph.connect_signal(sensor_port, input)?;
    let created_nodes = wire(graph, &scope, law_name, output, target, drive)?;

    let (_, parameter_sources) = recorder.finish();
    let (sensed, signal) = &ctrl.sensor_ports[0];
    let (entity, port) = &ctrl.actuator_port;
    Ok(SourceLoweringReceipt {
        source_name: scope,
        target_description: format!("{}.{port} from {}.{signal}", entity.token(), sensed.token()),
        created_nodes,
        parameter_sources,
    })
}

/// Join the law's `output` to its drive on `target`; returns the new nodes'
/// names, the law's first. A heater does not cool, so its command passes a
/// limiter at `[0, f64::MAX]` first.
fn wire(
    graph: &mut LoweringGraph,
    scope: &str,
    law_name: String,
    output: PortRef,
    target: PortRef,
    drive: Drive,
) -> Result<Vec<String>, LoweringError> {
    let drive_name = format!("{law_name}.drive");
    let account = SourceAccount::default();
    let mut created = vec![law_name];
    match drive {
        Drive::Effort(kind) => {
            let source = ControlledSource::effort(kind);
            let node = graph.add_with_tangent_and_residual(&drive_name, None, source, account)?;
            graph.connect_signal(output, graph.port(node, "command")?)?;
            graph.attach_common_flow(target, graph.port(node, "p")?, 1.0);
        }
        Drive::Heat => {
            let limit_name = format!("{}.limit", created[0]);
            let limiter = SignalLimiter::new(vec![0.0], vec![f64::MAX]).map_err(|source| {
                LoweringError::Element {
                    item: scope.to_owned(),
                    source,
                }
            })?;
            let limit = graph.add_with_tangent(&limit_name, None, limiter, ())?;
            graph.connect_signal(output, graph.port(limit, "input")?)?;
            let heater = Heater::commanded();
            let node = graph.add_with_tangent_and_residual(&drive_name, None, heater, account)?;
            graph.connect_signal(graph.port(limit, "output")?, graph.port(node, "command")?)?;
            graph.attach_common_effort(target, graph.port(node, "heat")?, 1.0);
            created.push(limit_name);
        }
    }
    created.push(drive_name);
    Ok(created)
}
