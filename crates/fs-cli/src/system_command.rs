//! `frankensim system`: a device described as data, run on the composition
//! kernel.
//!
//! The scene names elements from [`fs_system::elements`] and the connections
//! between their ports; this module lowers it to a [`SystemBuilder`], steps
//! the sealed [`System`] for the requested duration, and reports the final
//! states, sampled states and the energy and entropy ledger. Nothing here
//! integrates anything: every number comes from the kernel and its elements.

use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::Read as _;
use std::path::Path;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::{Carrier, PortRole, PortSpec, PortTangent};
use fs_system::elements::{
    Conductance, Damper, Force, Gyrator, HeatCapacity, Heater, Mass, Motion, SourceAccount,
    Spring, Waveform,
};
use fs_system::{
    NodeId, PortRef, RouteTemperature, SolveConfig, System, SystemBuilder, SystemError, heat_sink,
};

use crate::json_read::JsonValue as J;
use crate::{CommandOutput, exit};

const SCHEMA: &str = "frankensim.system.v0";
const RESULT_SCHEMA: &str = "frankensim.system.result.v0";
const DIAGNOSTIC_SCHEMA: &str = "frankensim.system.diagnostic.v0";
const MAX_SCENE_BYTES: u64 = 1024 * 1024;
// A looped block is solved with a dense matrix over its interface inputs:
// about 120 MB per step at 2,400 elements, growing with the square.
const MAX_ELEMENTS: usize = 2_000;
const MAX_STEPS: f64 = 10_000_000.0;
const MAX_SAMPLES: usize = 100_000;
const NO_CLAIM: &str = "ideal lumped elements with constant parameters integrated by the implicit midpoint rule; no material cards, contact, friction, phase change or step-size error control; the ledger states what the discrete system did, not what a physical device would do";

const HELP: &str = "Usage: frankensim [--json] system <scene.json>

Run a lumped mechanical, rotational, electrical and thermal device described
as data on the composition kernel, and report its final state, sampled states
and the energy and entropy ledger of the run.

Scene schema: frankensim.system.v0. Top-level fields:
  duration_s        simulated time
  step_s            largest step; the kernel may take shorter ones
  sample_every_s    optional; record every element's state at the end of the
                    first step at or after each multiple of this interval
  solver            optional {relative_tolerance, absolute_tolerance, max_iterations}
  elements          object: name -> element
  connections       array of connections

Elements (ports in brackets; SI units in every field name):
  mass           mass_kg, position_m = 0, velocity_m_s = 0            [p]
  spring         stiffness_n_m, cubic_n_m3 = 0, extension_m = 0       [a, b]
  damper         damping_n_s_m                                        [a, b, heat]
  force          force_n (waveform)                                   [p]
  heat_capacity  capacity_j_k, temperature_k                          [heat]
  conductance    conductance_w_k                                      [a, b]
  heater         power_w (waveform, never negative)                   [heat]
The same four laws about a shaft, with angle, angular velocity and torque:
  inertia        inertia_kg_m2, angle_rad = 0, angular_velocity_rad_s = 0  [p]
  torsion_spring stiffness_n_m_rad, cubic_n_m_rad3 = 0, twist_rad = 0  [a, b]
  rotary_damper  damping_n_m_s_rad                                    [a, b, heat]
  torque         torque_n_m (waveform)                                [p]
and in a current loop, with charge, current and voltage (an inductor is the
loop's node; a capacitor and a resistor sit in series in it):
  inductor       inductance_h, charge_c = 0, current_a = 0            [p]
  capacitor      capacitance_f, charge_c = 0                          [a, b]
  resistor       resistance_ohm                                       [a, b, heat]
  voltage        voltage_v (waveform)                                 [p]
Between two domains, without loss (an ideal DC machine with a = electrical,
b = rotation and ratio in N m/A: torque ratio * i on the rotor, back-EMF
ratio * w against the loop):
  gyrator        ratio, a, b (translation | rotation | electrical)    [a, b]
A waveform is a number, {\"sine\": {amplitude, frequency_hz, phase_rad = 0,
offset = 0}} or {\"table\": [[time_s, value], ...]} (linear, ends held).

Connections, ports written element.port:
  {\"join\": [ports]}       the ports meet at one node: a mass (inertia,
                          inductor) and the ends that share its velocity
                          (angular velocity, current), or a heat capacity and
                          the elements at its temperature
  {\"route\": port, \"to\": [[port, fraction], ...],
   \"rest\": {\"sink\": name, \"temperature_k\": T},
   \"temperature\": \"sink-equivalent\" | {\"fixed_k\": T}}
                          send the heat leaving a port by fraction
  {\"prescribe\": port, <quantity>: waveform}
                          hold a port's input from outside: velocity_m_s for a
                          spring, damper or force end (0 anchors it), force_n
                          for a mass; angular_velocity_rad_s and torque_n_m
                          about a shaft; current_a and voltage_v in a loop (0
                          current grounds a series element's far end);
                          temperature_k for a damper, resistor, conductance or
                          heater, heat_rate_w for a heat capacity

Each element's port a is the end at mass A and b the end at mass B; the
extension of a spring is x_B - x_A less its unloaded length; a capacitor's
charge_c is the same difference, so one charged by a current entering at a
reads negative. Every port must be connected. See examples/system/.
";

pub fn run(args: &[OsString], json: bool) -> CommandOutput {
    if args.len() == 1 && matches!(args[0].to_str(), Some("--help" | "-h")) {
        return CommandOutput {
            exit_code: exit::SUCCESS,
            stdout: HELP.into(),
            stderr: String::new(),
        };
    }
    if args.len() != 1 {
        return refusal(exit::USAGE, "system-usage", HELP, json);
    }
    let path = Path::new(&args[0]);
    let mut device = match load(path) {
        Ok(device) => device,
        Err(Refusal::Input(message)) => return refusal(exit::INPUT, "system-input", &message, json),
        Err(Refusal::Admission(message)) => {
            return refusal(exit::REFUSED, "system-admission", &message, json);
        }
    };
    match device.run() {
        Ok(()) => CommandOutput {
            exit_code: exit::SUCCESS,
            stdout: if json {
                device.render_json(path)
            } else {
                device.render_text(path)
            },
            stderr: String::new(),
        },
        Err(RunError::Step(error)) => {
            let code = if matches!(error, SystemError::Cancelled) {
                exit::CANCELLED
            } else {
                exit::REFUSED
            };
            let message = format!(
                "step {} at t = {} s: {error}",
                device.system.steps() + 1,
                device.system.time_s()
            );
            refusal(code, "system-step", &message, json)
        }
        Err(RunError::TooManySteps) => {
            let message = format!(
                "the kernel took {} steps and reached t = {} s of {} s; components keep \
                 shortening the step",
                device.system.steps(),
                device.system.time_s(),
                device.duration_s
            );
            refusal(exit::REFUSED, "system-steps", &message, json)
        }
    }
}

enum RunError {
    Step(SystemError),
    TooManySteps,
}

enum Refusal {
    /// The file could not be read or does not follow the schema.
    Input(String),
    /// The kernel refused the description.
    Admission(String),
}

impl From<String> for Refusal {
    fn from(message: String) -> Self {
        Self::Input(message)
    }
}

impl From<&str> for Refusal {
    fn from(message: &str) -> Self {
        Self::Input(message.to_owned())
    }
}

/// A power domain of the lumped elements.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Domain {
    Translation,
    Rotation,
    Electrical,
}

/// Field names of one domain, SI units in each.
struct Names {
    /// Parameter of the node element (mass, inertia, inductor).
    inertance: &'static str,
    /// The node's displacement and flow.
    displacement: &'static str,
    flow: &'static str,
    /// Parameters and state of the two-port storage element.
    stiffness: &'static str,
    cubic: &'static str,
    extension: &'static str,
    /// Parameter of the two-port dissipative element.
    resistance: &'static str,
    /// The effort: what a source prescribes and a dissipative element returns.
    effort: &'static str,
}

impl Domain {
    fn parse(name: &str) -> Option<Self> {
        Some(match name {
            "translation" => Self::Translation,
            "rotation" => Self::Rotation,
            "electrical" => Self::Electrical,
            _ => return None,
        })
    }

    fn port_kind(self) -> PortKind {
        match self {
            Self::Translation => PortKind::MechanicalForceVelocity,
            Self::Rotation => PortKind::RotationalTorqueAngularVelocity,
            Self::Electrical => PortKind::ElectricalVoltageCurrent,
        }
    }

    fn of(kind: PortKind) -> Option<Self> {
        [Self::Translation, Self::Rotation, Self::Electrical]
            .into_iter()
            .find(|domain| domain.port_kind() == kind)
    }

    fn names(self) -> Names {
        match self {
            Self::Translation => Names {
                inertance: "mass_kg",
                displacement: "position_m",
                flow: "velocity_m_s",
                stiffness: "stiffness_n_m",
                cubic: "cubic_n_m3",
                extension: "extension_m",
                resistance: "damping_n_s_m",
                effort: "force_n",
            },
            Self::Rotation => Names {
                inertance: "inertia_kg_m2",
                displacement: "angle_rad",
                flow: "angular_velocity_rad_s",
                stiffness: "stiffness_n_m_rad",
                cubic: "cubic_n_m_rad3",
                extension: "twist_rad",
                resistance: "damping_n_m_s_rad",
                effort: "torque_n_m",
            },
            // A capacitor is declared by its capacitance; its stiffness is
            // the reciprocal and it has no cubic term.
            Self::Electrical => Names {
                inertance: "inductance_h",
                displacement: "charge_c",
                flow: "current_a",
                stiffness: "capacitance_f",
                cubic: "",
                extension: "charge_c",
                resistance: "resistance_ohm",
                effort: "voltage_v",
            },
        }
    }
}

/// What an element is: one of four laws in a power domain, a gyrator
/// between two, or a thermal element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// The node of a domain: momentum over flow (mass, inertia, inductor).
    Inertance(Domain),
    /// Effort from displacement (spring, torsion spring, capacitor).
    Compliance(Domain),
    /// Effort from flow, with a heat port (damper, rotary damper, resistor).
    Resistance(Domain),
    /// A prescribed effort (force, torque, voltage).
    Source(Domain),
    Gyrator,
    HeatCapacity,
    Conductance,
    Heater,
}

impl Kind {
    const TABLE: [(&'static str, Self); 16] = [
        ("mass", Self::Inertance(Domain::Translation)),
        ("spring", Self::Compliance(Domain::Translation)),
        ("damper", Self::Resistance(Domain::Translation)),
        ("force", Self::Source(Domain::Translation)),
        ("inertia", Self::Inertance(Domain::Rotation)),
        ("torsion_spring", Self::Compliance(Domain::Rotation)),
        ("rotary_damper", Self::Resistance(Domain::Rotation)),
        ("torque", Self::Source(Domain::Rotation)),
        ("inductor", Self::Inertance(Domain::Electrical)),
        ("capacitor", Self::Compliance(Domain::Electrical)),
        ("resistor", Self::Resistance(Domain::Electrical)),
        ("voltage", Self::Source(Domain::Electrical)),
        ("gyrator", Self::Gyrator),
        ("heat_capacity", Self::HeatCapacity),
        ("conductance", Self::Conductance),
        ("heater", Self::Heater),
    ];

    fn parse(name: &str) -> Option<Self> {
        Self::TABLE
            .iter()
            .find(|(known, _)| *known == name)
            .map(|(_, kind)| *kind)
    }

    fn name(self) -> &'static str {
        Self::TABLE
            .iter()
            .find(|(_, kind)| *kind == self)
            .map_or("element", |(name, _)| name)
    }
}

struct Element {
    name: String,
    kind: Kind,
    node: NodeId,
    ports: Vec<PortSpec>,
    /// Initial temperature of a heat capacity, used as the starting value of
    /// the temperature inputs joined to it.
    temperature_k: Option<f64>,
}

struct Device {
    system: System,
    elements: Vec<Element>,
    connections: usize,
    duration_s: f64,
    step_s: f64,
    sample_every_s: Option<f64>,
    samples: Vec<(f64, Vec<(String, f64)>)>,
    iterations: u64,
    shortenings: u64,
}

fn read(path: &Path) -> Result<String, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let metadata = file.metadata().map_err(|e| e.to_string())?;
    if !metadata.is_file() || metadata.len() > MAX_SCENE_BYTES {
        return Err(format!(
            "{}: expected a regular file of at most {MAX_SCENE_BYTES} bytes",
            path.display()
        ));
    }
    let mut bytes = Vec::new();
    file.take(MAX_SCENE_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() as u64 > MAX_SCENE_BYTES {
        return Err("input grew beyond its byte limit".into());
    }
    String::from_utf8(bytes).map_err(|e| e.to_string())
}

fn fields(value: &J, allowed: &[&str], what: &str) -> Result<(), String> {
    let object = value
        .as_object()
        .ok_or_else(|| format!("{what}: expected an object"))?;
    match object
        .iter()
        .find(|(key, _)| !allowed.contains(&key.as_str()))
    {
        Some((key, _)) => Err(format!("{what}: unsupported field {key:?}")),
        None => Ok(()),
    }
}

fn number(value: &J, key: &str, what: &str) -> Result<f64, String> {
    value
        .get(key)
        .ok_or_else(|| format!("{what}: missing {key}"))?
        .as_f64()
        .ok_or_else(|| format!("{what}: {key} must be a number"))
}

fn number_or(value: &J, key: &str, default: f64, what: &str) -> Result<f64, String> {
    match value.get(key) {
        None => Ok(default),
        Some(field) => field
            .as_f64()
            .ok_or_else(|| format!("{what}: {key} must be a number")),
    }
}

fn positive(value: &J, key: &str, what: &str) -> Result<f64, String> {
    let number = number(value, key, what)?;
    if number.is_finite() && number > 0.0 {
        Ok(number)
    } else {
        Err(format!("{what}: {key} must be positive and finite"))
    }
}

fn waveform(value: &J, what: &str) -> Result<Waveform, String> {
    if let Some(constant) = value.as_f64() {
        return Waveform::constant(constant).map_err(|e| format!("{what}: {e}"));
    }
    let object = value
        .as_object()
        .ok_or_else(|| format!("{what}: a waveform is a number, a sine or a table"))?;
    match object {
        [(key, sine)] if key == "sine" => {
            let inner = format!("{what}.sine");
            fields(sine, &["amplitude", "frequency_hz", "phase_rad", "offset"], &inner)?;
            Waveform::sine(
                number(sine, "amplitude", &inner)?,
                number(sine, "frequency_hz", &inner)?,
                number_or(sine, "phase_rad", 0.0, &inner)?,
                number_or(sine, "offset", 0.0, &inner)?,
            )
            .map_err(|e| format!("{inner}: {e}"))
        }
        [(key, table)] if key == "table" => {
            let rows = table
                .as_array()
                .ok_or_else(|| format!("{what}.table: expected an array of [time_s, value]"))?;
            let knots = rows
                .iter()
                .map(|row| match row.as_array() {
                    Some([t, v]) => t.as_f64().zip(v.as_f64()),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| format!("{what}.table: every row is [time_s, value]"))?;
            Waveform::table(knots).map_err(|e| format!("{what}: {e}"))
        }
        _ => Err(format!(
            "{what}: a waveform object has exactly one key, sine or table"
        )),
    }
}

fn add<C>(
    builder: &mut SystemBuilder,
    name: &str,
    element: C,
    state: C::State,
) -> Result<(NodeId, Vec<PortSpec>), Refusal>
where
    C: PortTangent + 'static,
    C::State: 'static,
{
    let ports = element.ports().to_vec();
    let node = builder
        .add_with_tangent(name, element, state)
        .map_err(|e| Refusal::Admission(e.to_string()))?;
    Ok((node, ports))
}

/// A spring, torsion spring or capacitor. A capacitor is declared by its
/// capacitance and has no cubic term.
fn compliance(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    domain: Domain,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    let names = domain.names();
    let (stiffness, cubic) = if domain == Domain::Electrical {
        fields(definition, &["kind", names.stiffness, names.extension], what)?;
        (1.0 / positive(definition, names.stiffness, what)?, 0.0)
    } else {
        fields(definition, &["kind", names.stiffness, names.cubic, names.extension], what)?;
        (
            number(definition, names.stiffness, what)?,
            number_or(definition, names.cubic, 0.0, what)?,
        )
    };
    let extension = number_or(definition, names.extension, 0.0, what)?;
    if !extension.is_finite() {
        return Err(format!("{what}: {} must be finite", names.extension).into());
    }
    let spring = Spring::in_domain(domain.port_kind(), stiffness, cubic).map_err(|_| {
        format!("{what}: {} must be positive and the cubic term finite", names.stiffness)
    })?;
    add(builder, name, spring, extension)
}

fn gyrator(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(definition, &["kind", "ratio", "a", "b"], what)?;
    let side = |key: &str| {
        definition
            .str_field(key)
            .and_then(Domain::parse)
            .ok_or_else(|| format!("{what}: {key} is one of translation, rotation, electrical"))
    };
    let (a, b) = (side("a")?, side("b")?);
    let gyrator = Gyrator::new(a.port_kind(), b.port_kind(), number(definition, "ratio", what)?)
        .map_err(|e| Refusal::Input(format!("{what}: {e}")))?;
    add(builder, name, gyrator, ())
}

fn element(builder: &mut SystemBuilder, name: &str, definition: &J) -> Result<Element, Refusal> {
    let what = format!("elements.{name}");
    let kind_name = definition
        .str_field("kind")
        .ok_or_else(|| format!("{what}: missing kind"))?;
    let kind = Kind::parse(kind_name).ok_or_else(|| format!("{what}: unknown kind {kind_name:?}"))?;
    let invalid = |e: fs_system::elements::ElementError| Refusal::Input(format!("{what}: {e}"));
    let mut temperature_k = None;
    let (node, ports) = match kind {
        Kind::Inertance(domain) => {
            let names = domain.names();
            fields(definition, &["kind", names.inertance, names.displacement, names.flow], &what)?;
            let state = Motion {
                position_m: number_or(definition, names.displacement, 0.0, &what)?,
                velocity_m_s: number_or(definition, names.flow, 0.0, &what)?,
            };
            if !(state.position_m.is_finite() && state.velocity_m_s.is_finite()) {
                return Err(format!(
                    "{what}: {} and {} must be finite",
                    names.displacement, names.flow
                )
                .into());
            }
            let inertance = number(definition, names.inertance, &what)?;
            let node = Mass::in_domain(domain.port_kind(), inertance)
                .map_err(|_| format!("{what}: {} must be positive", names.inertance))?;
            add(builder, name, node, state)?
        }
        Kind::Compliance(domain) => compliance(builder, name, definition, domain, &what)?,
        Kind::Resistance(domain) => {
            let names = domain.names();
            fields(definition, &["kind", names.resistance], &what)?;
            let resistance = number(definition, names.resistance, &what)?;
            let damper = Damper::in_domain(domain.port_kind(), resistance)
                .map_err(|_| format!("{what}: {} must be positive", names.resistance))?;
            add(builder, name, damper, ())?
        }
        Kind::Source(domain) => {
            let names = domain.names();
            fields(definition, &["kind", names.effort], &what)?;
            let wave = definition
                .get(names.effort)
                .ok_or_else(|| format!("{what}: missing {}", names.effort))?;
            let source = Force::in_domain(
                domain.port_kind(),
                waveform(wave, &format!("{what}.{}", names.effort))?,
            );
            add(builder, name, source, SourceAccount::default())?
        }
        Kind::Gyrator => gyrator(builder, name, definition, &what)?,
        Kind::HeatCapacity => {
            fields(definition, &["kind", "capacity_j_k", "temperature_k"], &what)?;
            let start = positive(definition, "temperature_k", &what)?;
            temperature_k = Some(start);
            let body = HeatCapacity::new(number(definition, "capacity_j_k", &what)?).map_err(invalid)?;
            add(builder, name, body, start)?
        }
        Kind::Conductance => {
            fields(definition, &["kind", "conductance_w_k"], &what)?;
            let wall =
                Conductance::new(number(definition, "conductance_w_k", &what)?).map_err(invalid)?;
            add(builder, name, wall, ())?
        }
        Kind::Heater => {
            fields(definition, &["kind", "power_w"], &what)?;
            let wave = definition
                .get("power_w")
                .ok_or_else(|| format!("{what}: missing power_w"))?;
            let heater = Heater::new(waveform(wave, &format!("{what}.power_w"))?).map_err(invalid)?;
            add(builder, name, heater, SourceAccount::default())?
        }
    };
    Ok(Element {
        name: name.to_owned(),
        kind,
        node,
        ports,
        temperature_k,
    })
}

/// A port named in the scene, resolved.
#[derive(Clone, Copy)]
struct Resolved {
    element: usize,
    port: PortRef,
    carrier: Carrier,
    role: PortRole,
}

fn resolve(
    builder: &SystemBuilder,
    elements: &[Element],
    text: &str,
) -> Result<Resolved, Refusal> {
    let (name, port) = text
        .rsplit_once('.')
        .ok_or_else(|| format!("port {text:?}: write it as element.port"))?;
    let element = elements
        .iter()
        .position(|element| element.name == name)
        .ok_or_else(|| format!("port {text:?}: no element named {name:?}"))?;
    let spec = elements[element]
        .ports
        .iter()
        .find(|spec| spec.id().as_str() == port)
        .ok_or_else(|| {
            let known = elements[element]
                .ports
                .iter()
                .map(|spec| spec.id().as_str())
                .collect::<Vec<_>>()
                .join(", ");
            format!(
                "port {text:?}: a {} has ports {known}",
                elements[element].kind.name()
            )
        })?;
    let port = builder
        .port(elements[element].node, port)
        .map_err(|e| Refusal::Admission(e.to_string()))?;
    Ok(Resolved {
        element,
        port,
        carrier: spec.carrier(),
        role: spec.role(),
    })
}

fn port_text(value: &J, what: &str) -> Result<String, String> {
    value
        .as_str()
        .map(str::to_owned)
        .ok_or_else(|| format!("{what}: a port is a string element.port"))
}

/// Starting values for temperature inputs inside loops.
type Guesses = Vec<(PortRef, f64)>;

fn join(
    builder: &mut SystemBuilder,
    elements: &[Element],
    ports: &[Resolved],
    guesses: &mut Guesses,
    what: &str,
) -> Result<(), Refusal> {
    if ports.len() < 2 {
        return Err(format!("{what}: a join needs at least two ports").into());
    }
    // The node is the one mass or heat capacity in the join: the mass sets
    // the velocity the ends it carries move with, the heat capacity the
    // temperature of the elements that sit at it.
    let is_node = |port: &Resolved| {
        matches!(
            (port.carrier, port.role),
            (Carrier::Power(_), PortRole::EffortIn) | (Carrier::Heat, PortRole::FlowIn)
        )
    };
    let nodes = ports.iter().filter(|port| is_node(port)).count();
    if nodes != 1 {
        return Err(format!(
            "{what}: a join has exactly one node port (a mass, inertia, inductor or heat \
             capacity) and the ports that meet at it; this one has {nodes}"
        )
        .into());
    }
    let setter = ports.iter().position(is_node).unwrap_or_default();
    let setter_role = ports[setter].role;
    let members: Vec<PortRef> = ports
        .iter()
        .enumerate()
        .filter(|&(index, _)| index != setter)
        .map(|(_, p)| p.port)
        .collect();
    let result = if setter_role == PortRole::EffortIn {
        builder.common_flow(ports[setter].port, &members)
    } else {
        builder.common_effort(ports[setter].port, &members)
    };
    result.map_err(|e| Refusal::Admission(format!("{what}: {e}")))?;
    if let (Carrier::Heat, Some(kelvin)) =
        (ports[setter].carrier, elements[ports[setter].element].temperature_k)
    {
        guesses.extend(members.iter().map(|&member| (member, kelvin)));
    }
    Ok(())
}

fn route(
    builder: &mut SystemBuilder,
    elements: &[Element],
    connection: &J,
    guesses: &mut Guesses,
    what: &str,
) -> Result<(), Refusal> {
    fields(connection, &["route", "to", "rest", "temperature"], what)?;
    let source = resolve(builder, elements, &port_text(connection.get("route").unwrap_or(&J::Null), what)?)?;
    let mut sinks = Vec::new();
    // Receivers' starting temperatures, for the source's starting value.
    let mut weights = Vec::new();
    if let Some(to) = connection.get("to") {
        let rows = to
            .as_array()
            .ok_or_else(|| format!("{what}: to is an array of [port, fraction]"))?;
        for row in rows {
            let (port, fraction) = match row.as_array() {
                Some([port, fraction]) => (port_text(port, what)?, fraction.as_f64()),
                _ => (String::new(), None),
            };
            let fraction = fraction.ok_or_else(|| format!("{what}: to is an array of [port, fraction]"))?;
            let sink = resolve(builder, elements, &port)?;
            if let Some(kelvin) = elements[sink.element].temperature_k {
                weights.push((fraction, kelvin));
            }
            sinks.push((sink.port, fraction));
        }
    }
    let rest = match connection.get("rest") {
        None => None,
        Some(rest) => {
            let inner = format!("{what}.rest");
            fields(rest, &["sink", "temperature_k"], &inner)?;
            let name = rest
                .str_field("sink")
                .ok_or_else(|| format!("{inner}: missing sink"))?;
            let kelvin = positive(rest, "temperature_k", &inner)?;
            let allocated: f64 = sinks.iter().map(|(_, fraction)| fraction).sum();
            weights.push(((1.0 - allocated).max(0.0), kelvin));
            Some(heat_sink(name, kelvin).map_err(|e| format!("{inner}: {e:?}"))?)
        }
    };
    let temperature = match connection.get("temperature") {
        None => RouteTemperature::SinkEquivalent,
        Some(J::Str(text)) if text == "sink-equivalent" => RouteTemperature::SinkEquivalent,
        Some(fixed) => {
            let inner = format!("{what}.temperature");
            fields(fixed, &["fixed_k"], &inner)
                .map_err(|_| format!("{inner}: \"sink-equivalent\" or {{\"fixed_k\": T}}"))?;
            RouteTemperature::Fixed(positive(fixed, "fixed_k", &inner)?)
        }
    };
    builder
        .heat_route(source.port, &sinks, rest, temperature)
        .map_err(|e| Refusal::Admission(format!("{what}: {e}")))?;
    // The source receives sum(f) / sum(f / T) of the receivers.
    let (sum, inverse) = weights
        .iter()
        .filter(|(fraction, _)| *fraction > 0.0)
        .fold((0.0, 0.0), |(s, i), (f, t)| (s + f, i + f / t));
    let start = match temperature {
        RouteTemperature::Fixed(kelvin) => Some(kelvin),
        RouteTemperature::SinkEquivalent if inverse > 0.0 => Some(sum / inverse),
        RouteTemperature::SinkEquivalent => None,
    };
    if let Some(kelvin) = start {
        guesses.push((source.port, kelvin));
    }
    Ok(())
}

fn prescribe(
    builder: &mut SystemBuilder,
    elements: &[Element],
    connection: &J,
    what: &str,
) -> Result<(), Refusal> {
    let target = resolve(builder, elements, &port_text(connection.get("prescribe").unwrap_or(&J::Null), what)?)?;
    let domain = match target.carrier {
        Carrier::Power(kind) => Domain::of(kind),
        _ => None,
    };
    let quantity = match (target.carrier, target.role, domain) {
        (Carrier::Heat, PortRole::EffortIn, _) => "temperature_k",
        (Carrier::Heat, PortRole::FlowIn, _) => "heat_rate_w",
        (Carrier::Power(_), PortRole::EffortIn, Some(domain)) => domain.names().effort,
        (Carrier::Power(_), PortRole::FlowIn, Some(domain)) => domain.names().flow,
        _ => return Err(format!("{what}: this port cannot be prescribed").into()),
    };
    fields(connection, &["prescribe", quantity], what)
        .map_err(|_| format!("{what}: this port takes exactly one quantity, {quantity}"))?;
    let wave = waveform(
        connection
            .get(quantity)
            .ok_or_else(|| format!("{what}: missing {quantity}"))?,
        &format!("{what}.{quantity}"),
    )?;
    if quantity == "temperature_k" && !(wave.lowest() > 0.0) {
        return Err(format!("{what}: temperature_k must stay positive").into());
    }
    // A port carried by a node receives minus the node's flow, so a port
    // prescribed to move at v (or to carry the current i) receives -v.
    let prescribed_flow =
        matches!((target.carrier, target.role), (Carrier::Power(_), PortRole::FlowIn));
    let sign = if prescribed_flow { -1.0 } else { 1.0 };
    builder
        .boundary(
            target.port,
            Box::new(move |t0, dt, values| values[0] = sign * wave.mean(t0, dt)),
        )
        .map_err(|e| Refusal::Admission(format!("{what}: {e}")))?;
    Ok(())
}

fn solver(definition: &J) -> Result<SolveConfig, String> {
    let what = "solver";
    fields(
        definition,
        &["relative_tolerance", "absolute_tolerance", "max_iterations"],
        what,
    )?;
    let mut config = SolveConfig::default();
    config.relative_tolerance =
        number_or(definition, "relative_tolerance", config.relative_tolerance, what)?;
    config.absolute_tolerance =
        number_or(definition, "absolute_tolerance", config.absolute_tolerance, what)?;
    let iterations = number_or(
        definition,
        "max_iterations",
        f64::from(config.max_iterations),
        what,
    )?;
    let admitted = |value: f64| value.is_finite() && value >= 0.0;
    if !admitted(config.relative_tolerance) || !admitted(config.absolute_tolerance) {
        return Err("solver: tolerances must be finite and not negative".into());
    }
    if !(1.0..=1000.0).contains(&iterations) || iterations.fract() != 0.0 {
        return Err("solver: max_iterations must be a whole number from 1 to 1000".into());
    }
    config.max_iterations = iterations as u32;
    Ok(config)
}

fn load(path: &Path) -> Result<Device, Refusal> {
    let text = read(path)?;
    let root = J::parse(&text).map_err(|e| e.to_string())?;
    fields(
        &root,
        &[
            "schema",
            "duration_s",
            "step_s",
            "sample_every_s",
            "solver",
            "elements",
            "connections",
        ],
        "scene",
    )?;
    if root.str_field("schema") != Some(SCHEMA) {
        return Err(format!("scene: schema must be {SCHEMA:?}").into());
    }
    let duration_s = positive(&root, "duration_s", "scene")?;
    let step_s = positive(&root, "step_s", "scene")?;
    if duration_s / step_s > MAX_STEPS {
        return Err(format!("scene: more than {MAX_STEPS} steps of step_s").into());
    }
    let sample_every_s = if root.get("sample_every_s").is_some() {
        let every = positive(&root, "sample_every_s", "scene")?;
        if duration_s / every > MAX_SAMPLES as f64 {
            return Err(format!("scene: more than {MAX_SAMPLES} samples").into());
        }
        Some(every)
    } else {
        None
    };
    let definitions = root
        .get("elements")
        .and_then(J::as_object)
        .ok_or("scene: elements must be an object of name -> element")?;
    if definitions.is_empty() || definitions.len() > MAX_ELEMENTS {
        return Err(format!("scene: from 1 to {MAX_ELEMENTS} elements").into());
    }
    let mut builder = SystemBuilder::new();
    let mut elements = Vec::with_capacity(definitions.len());
    for (name, definition) in definitions {
        elements.push(element(&mut builder, name, definition)?);
    }
    let connections = root
        .get("connections")
        .and_then(J::as_array)
        .ok_or("scene: connections must be an array")?;
    let mut guesses = Guesses::new();
    for (index, connection) in connections.iter().enumerate() {
        let what = format!("connections[{index}]");
        if let Some(list) = connection.get("join") {
            fields(connection, &["join"], &what)?;
            let ports = list
                .as_array()
                .ok_or_else(|| format!("{what}: join is an array of ports"))?
                .iter()
                .map(|port| resolve(&builder, &elements, &port_text(port, &what)?))
                .collect::<Result<Vec<_>, _>>()?;
            join(&mut builder, &elements, &ports, &mut guesses, &what)?;
        } else if connection.get("route").is_some() {
            route(&mut builder, &elements, connection, &mut guesses, &what)?;
        } else if connection.get("prescribe").is_some() {
            prescribe(&mut builder, &elements, connection, &what)?;
        } else {
            return Err(format!("{what}: expected join, route or prescribe").into());
        }
    }
    let mut system = builder
        .seal()
        .map_err(|e| Refusal::Admission(e.to_string()))?;
    if let Some(definition) = root.get("solver") {
        system.set_config(solver(definition)?);
    }
    for (port, kelvin) in guesses {
        system.guess(port, &[kelvin]);
    }
    Ok(Device {
        system,
        elements,
        connections: connections.len(),
        duration_s,
        step_s,
        sample_every_s,
        samples: Vec::new(),
        iterations: 0,
        shortenings: 0,
    })
}

impl Device {
    fn run(&mut self) -> Result<(), RunError> {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let key = StreamKey {
                seed: 0,
                kernel_id: 0x5953,
                tile: 0,
                iteration: 0,
            };
            let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
            self.sample();
            let mut next_sample = self.sample_every_s.unwrap_or(f64::INFINITY);
            // A step shorter than this would only round the end time.
            let slack = 1.0e-12 * self.duration_s;
            while self.duration_s - self.system.time_s() > slack {
                if self.system.steps() as f64 >= 2.0 * MAX_STEPS {
                    return Err(RunError::TooManySteps);
                }
                let remaining = self.duration_s - self.system.time_s();
                let report = self
                    .system
                    .step(self.step_s.min(remaining), &cx)
                    .map_err(RunError::Step)?;
                self.iterations += u64::from(report.iterations);
                self.shortenings += u64::from(report.shortenings);
                if report.time_s >= next_sample - slack {
                    self.sample();
                    let every = self.sample_every_s.unwrap_or(f64::INFINITY);
                    while next_sample <= report.time_s + slack {
                        next_sample += every;
                    }
                }
            }
            if self.sample_every_s.is_some()
                && self.samples.last().is_none_or(|(t, _)| *t < self.system.time_s())
            {
                self.sample();
            }
            Ok(())
        })
    }

    fn sample(&mut self) {
        if self.sample_every_s.is_none() || self.samples.len() > MAX_SAMPLES {
            return;
        }
        let values = self
            .elements
            .iter()
            .flat_map(|element| {
                self.values(element)
                    .into_iter()
                    .map(|(key, value)| (format!("{}.{key}", element.name), value))
            })
            .collect();
        self.samples.push((self.system.time_s(), values));
    }

    /// The state of an element, and for stateless ones the step-mean values
    /// of the last step.
    fn values(&self, element: &Element) -> Vec<(&'static str, f64)> {
        let system = &self.system;
        // Before the first step there is no last step: those values are null.
        let output = |port: &str| match system.port(element.node, port) {
            Some(port) if system.steps() > 0 => system.last_output(port)[0],
            _ => f64::NAN,
        };
        match element.kind {
            Kind::Inertance(domain) => {
                let names = domain.names();
                let motion = system.state::<Motion>(element.node).copied().unwrap_or(Motion {
                    position_m: f64::NAN,
                    velocity_m_s: f64::NAN,
                });
                vec![(names.displacement, motion.position_m), (names.flow, motion.velocity_m_s)]
            }
            Kind::Compliance(domain) => vec![(
                domain.names().extension,
                system.state::<f64>(element.node).copied().unwrap_or(f64::NAN),
            )],
            Kind::HeatCapacity => vec![(
                "temperature_k",
                system.state::<f64>(element.node).copied().unwrap_or(f64::NAN),
            )],
            Kind::Source(_) | Kind::Heater => vec![(
                "energy_in_j",
                system
                    .state::<SourceAccount>(element.node)
                    .map_or(f64::NAN, |account| account.energy_in_j),
            )],
            Kind::Resistance(domain) => {
                vec![(domain.names().effort, output("a")), ("heat_w", -output("heat"))]
            }
            Kind::Gyrator => vec![("effort_a", output("a")), ("effort_b", output("b"))],
            Kind::Conductance => vec![("heat_rate_w", output("a"))],
        }
    }

    fn render_text(&self, path: &Path) -> String {
        let system = &self.system;
        let ledger = system.ledger();
        let (blocks, loops) = system.plan_shape();
        let mut out = String::new();
        let _ = writeln!(out, "frankensim system {}", path.display());
        let _ = writeln!(
            out,
            "{} elements, {} connections; plan: {blocks} block(s), {loops} solved by interface iteration",
            self.elements.len(),
            self.connections
        );
        let _ = writeln!(
            out,
            "t = {} s after {} steps ({} interface iterations, {} shortened)\n",
            system.time_s(),
            system.steps(),
            self.iterations,
            self.shortenings
        );
        let _ = writeln!(out, "{:<16} {:<14} state", "element", "kind");
        for element in &self.elements {
            let state = self
                .values(element)
                .iter()
                .map(|(key, value)| format!("{key} = {value:.9e}"))
                .collect::<Vec<_>>()
                .join(", ");
            let _ = writeln!(out, "{:<16} {:<14} {state}", element.name, element.kind.name());
        }
        let _ = writeln!(
            out,
            "\n{:<16} {:>16} {:>16} {:>16} {:>16}",
            "ledger [J]", "port energy in", "stored change", "dissipated", "entropy [J/K]"
        );
        for (element, account) in self.elements.iter().zip(&ledger.nodes) {
            let _ = writeln!(
                out,
                "{:<16} {:>16.6e} {:>16.6e} {:>16.6e} {:>16.6e}",
                element.name,
                account.port_energy_in_j.value(),
                account.stored_change_j.value(),
                account.dissipated_j.value(),
                account.entropy_generated_j_per_k.value()
            );
        }
        for boundary in &ledger.boundaries {
            let _ = writeln!(
                out,
                "boundary {}: energy in {:.6e} J, entropy in {:.6e} J/K",
                boundary.port,
                boundary.energy_in_j.value(),
                boundary.entropy_in_j_per_k.value()
            );
        }
        for sink in &ledger.sinks {
            let _ = writeln!(out, "sink {}: heat {:.6e} J", sink.name, sink.heat_j.value());
        }
        let _ = writeln!(
            out,
            "\nfirst law: stored change {:.6e} J = boundary energy in {:.6e} J - sink heat {:.6e} J \
             - numerical loss {:.6e} J + connection defect {:.3e} J + residual {:.3e} J",
            ledger.stored_change_j(),
            ledger.boundary_energy_in_j(),
            ledger.sink_heat_j(),
            ledger.numerical_loss_j(),
            ledger.connection_defect_j(),
            ledger.energy_residual_j()
        );
        let _ = writeln!(
            out,
            "second law: entropy generated {:.6e} J/K",
            ledger.entropy_generated_j_per_k()
        );
        if !self.samples.is_empty() {
            let _ = writeln!(out, "\n{} samples (use --json to read them)", self.samples.len());
        }
        let _ = writeln!(out, "\nNo claim: {NO_CLAIM}.");
        out
    }

    fn render_json(&self, path: &Path) -> String {
        let system = &self.system;
        let ledger = system.ledger();
        let (blocks, loops) = system.plan_shape();
        let mut out = String::new();
        let _ = write!(
            out,
            "{{\"schema\":{},\"scene\":{},\"time_s\":{},\"steps\":{},\"iterations\":{},\"shortenings\":{},\"plan\":{{\"blocks\":{blocks},\"loops\":{loops}}},\"elements\":[",
            quote(RESULT_SCHEMA),
            quote(&path.display().to_string()),
            num(system.time_s()),
            system.steps(),
            self.iterations,
            self.shortenings
        );
        for (index, (element, account)) in self.elements.iter().zip(&ledger.nodes).enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"name\":{},\"kind\":{},\"state\":{},\"ledger\":{{\"port_energy_in_j\":{},\"stored_change_j\":{},\"dissipated_j\":{},\"numerical_loss_j\":{},\"entropy_generated_j_per_k\":{},\"balance_defect_j\":{}}}}}",
                quote(&element.name),
                quote(element.kind.name()),
                object(self.values(element).iter().map(|(k, v)| ((*k).to_owned(), *v))),
                num(account.port_energy_in_j.value()),
                num(account.stored_change_j.value()),
                num(account.dissipated_j.value()),
                num(account.numerical_loss_j.value()),
                num(account.entropy_generated_j_per_k.value()),
                num(account.balance_defect_j.value())
            );
        }
        out.push_str("],\"boundaries\":[");
        for (index, boundary) in ledger.boundaries.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"port\":{},\"energy_in_j\":{},\"entropy_in_j_per_k\":{}}}",
                quote(&boundary.port),
                num(boundary.energy_in_j.value()),
                num(boundary.entropy_in_j_per_k.value())
            );
        }
        out.push_str("],\"sinks\":[");
        for (index, sink) in ledger.sinks.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"name\":{},\"heat_j\":{}}}",
                quote(&sink.name),
                num(sink.heat_j.value())
            );
        }
        let _ = write!(
            out,
            "],\"totals\":{{\"stored_change_j\":{},\"boundary_energy_in_j\":{},\"sink_heat_j\":{},\"numerical_loss_j\":{},\"energy_residual_j\":{},\"connection_defect_j\":{},\"entropy_generated_j_per_k\":{}}},\"samples\":[",
            num(ledger.stored_change_j()),
            num(ledger.boundary_energy_in_j()),
            num(ledger.sink_heat_j()),
            num(ledger.numerical_loss_j()),
            num(ledger.energy_residual_j()),
            num(ledger.connection_defect_j()),
            num(ledger.entropy_generated_j_per_k())
        );
        for (index, (time_s, values)) in self.samples.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"time_s\":{},\"values\":{}}}",
                num(*time_s),
                object(values.iter().map(|(k, v)| (k.clone(), *v)))
            );
        }
        let _ = writeln!(out, "],\"no_claim\":{}}}", quote(NO_CLAIM));
        out
    }
}

fn num(value: f64) -> String {
    if value.is_finite() {
        format!("{value:?}")
    } else {
        "null".into()
    }
}

fn object(entries: impl Iterator<Item = (String, f64)>) -> String {
    let body = entries
        .map(|(key, value)| format!("{}:{}", quote(&key), num(value)))
        .collect::<Vec<_>>()
        .join(",");
    format!("{{{body}}}")
}

fn quote(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c < '\u{20}' => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn refusal(exit_code: u8, code: &str, message: &str, json: bool) -> CommandOutput {
    CommandOutput {
        exit_code,
        stdout: String::new(),
        stderr: if json {
            format!(
                "{{\"schema\":{},\"code\":{},\"message\":{}}}\n",
                quote(DIAGNOSTIC_SCHEMA),
                quote(code),
                quote(message)
            )
        } else {
            format!("{code}: {message}\n")
        },
    }
}
