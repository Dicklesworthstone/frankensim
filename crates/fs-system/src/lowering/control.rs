//! Controllers from the description: a declared PID law, lowered onto the
//! kernel's held [`SampledController`]. It reads a signal a lumped body
//! publishes and commands an ideal effort source on the actuator's port.

use fs_ports::PortKind;
use fs_ports::component::{Carrier, PortRole};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{ControllerDecl, ControllerLaw};
use fs_scenario::entity::{EntityCatalog, EntityId};

use super::{
    FREQUENCY_DIMS, LoweringError, LoweringGraph, PartLoweringContext, ProvenanceRecorder,
    SourceLoweringReceipt, TIME_DIMS,
};
use crate::elements::{ControlledSource, ControllerState, Mass, SampledController, SourceAccount};
use crate::{NodeId, PortRef};

/// The signals a lumped rigid body publishes when a controller reads it:
/// position and velocity, or angle and angular velocity for a rotation.
const BODY_SIGNALS: [&str; 2] = ["position", "velocity"];

/// The part's mass, publishing its signals when a declared controller reads
/// them, and unchanged otherwise.
pub(super) fn publish_if_sensed(ctx: &PartLoweringContext<'_>, mass: Mass) -> Mass {
    let input = ctx.draft.input();
    let mut sensors = input.controllers.iter().flat_map(|c| &c.sensor_ports);
    if sensors.any(|(e, port)| *e == ctx.part.entity && BODY_SIGNALS.contains(&port.as_str())) {
        mass.publishing()
    } else {
        mass
    }
}

/// Expose on `entity` the signals `node` publishes, if it publishes any.
pub(super) fn expose_signals(graph: &mut LoweringGraph, entity: EntityId, node: NodeId) {
    for name in BODY_SIGNALS {
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
/// position or velocity, in its power port's coordinates.
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
    if !BODY_SIGNALS.contains(&sensor.as_str()) {
        return Err(refuse(
            scope,
            format!(
                "sensor '{sensor}': a lumped rigid body publishes position and velocity, \
                 and no other sensor signal lowers yet"
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
    let Some(&(Carrier::Power(kind), ..)) = graph.port_meta.get(&body) else {
        return Err(refuse(
            scope,
            format!("sensor '{sensor}' is not on a power-port body"),
        ));
    };
    let flow = kind.canonical_flow_dimensions();
    let dims = if sensor == "position" {
        ratio(flow, Dims::NONE, 1)
    } else {
        flow
    };
    Ok((port, dims))
}

/// The actuator's power port, an effort input a commanded source can drive,
/// and its kind.
fn actuator(
    ctrl: &ControllerDecl,
    scope: &str,
    catalog: &EntityCatalog,
    graph: &LoweringGraph,
) -> Result<(PortRef, PortKind), LoweringError> {
    let (entity, port) = (ctrl.actuator_port.0, ctrl.actuator_port.1.as_str());
    let target = graph
        .resolve_entity_port(catalog, entity, port)
        .ok_or_else(|| missing(scope, entity, port))?;
    match graph.port_meta.get(&target) {
        Some(&(Carrier::Power(kind), PortRole::EffortIn, 1)) => Ok((target, kind)),
        Some(&(Carrier::Heat, ..)) => Err(refuse(
            scope,
            format!("actuator '{port}' is a heat port, and no commanded heater exists yet"),
        )),
        _ => Err(refuse(
            scope,
            format!("actuator '{port}' takes no commanded effort; name a lumped body's power port"),
        )),
    }
}

/// Lower one declared controller: a PID law on one published sensor signal,
/// sampled at its declared rate and held between samples, commanding an
/// ideal effort source on the actuator's power port. Its receipt row joins
/// the sources as `controller:<name>`.
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
    let (target, kind) = actuator(ctrl, &scope, catalog, graph)?;
    let effort = kind.canonical_effort_dimensions();

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
    let drive_name = format!("{law_name}.drive");
    let law = graph.add_with_tangent(&law_name, None, law, memory)?;
    let drive = ControlledSource::effort(kind);
    let drive =
        graph.add_with_tangent_and_residual(&drive_name, None, drive, SourceAccount::default())?;
    let (input, output) = (graph.port(law, "input0")?, graph.port(law, "output")?);
    let (command, p) = (graph.port(drive, "command")?, graph.port(drive, "p")?);
    graph.connect_signal(sensor_port, input)?;
    graph.connect_signal(output, command)?;
    graph.attach_common_flow(target, p, 1.0);

    let (_, parameter_sources) = recorder.finish();
    let (sensed, signal) = &ctrl.sensor_ports[0];
    let (entity, port) = &ctrl.actuator_port;
    Ok(SourceLoweringReceipt {
        source_name: scope,
        target_description: format!("{}.{port} from {}.{signal}", entity.token(), sensed.token()),
        created_nodes: vec![law_name, drive_name],
        parameter_sources,
    })
}
