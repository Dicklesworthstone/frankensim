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

use std::rc::Rc;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_mbd::component::network::{DampingLaw, ViscousNetwork};
use fs_mbd::component::relaxation::{DriveFrame, RelaxationNetwork, RelaxationState, ShiftLaw};
use fs_mbd::spatial_loads::internal_impulse::InternalWrenchImpulseControls;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::driven::DrivenViscousNetworkLimits;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::relaxation::{
    RelaxationNetworkLimits, RelaxationNetworkLink, SpatialRelaxationBranch,
};
use fs_mbd::spatial_loads::internal_impulse::viscous::network::{
    ViscousNetworkLimits, ViscousNetworkLink,
};
use fs_mbd::spatial_loads::internal_motion::viscous::network::relaxation::RelaxationNetworkMotionControls;
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
};
use fs_mbd::spatial_loads::internal_motion::viscous::network::ViscousNetworkMotionControls;
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_ports::PortKind;
use fs_ports::component::{Carrier, Component, PortRole, PortSpec, PortTangent};
use fs_system::constitutive::{
    CapacitorChargeCurve, InductorFluxCurve, MAX_CAPACITOR_KNOTS, MAX_INDUCTOR_KNOTS,
    NonlinearCapacitor, NonlinearInductor,
};
use fs_system::elements::{
    Admittance, Conductance, Damper, FlowSource, Force, Gyrator, HeatCapacity, Heater, Mass,
    Motion, SourceAccount, Spring, Storage, Waveform,
};
use fs_system::{
    AdaptiveConfig, AdaptiveError, AdaptiveReport, Checkpoint, ErrorNorm, ErrorTolerance,
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
const NO_CLAIM: &str = "lumped elements with constant parameters or authored nonlinear electrical storage curves, using work-conjugate discrete-gradient stepping; no material-card resolution, contact, friction, phase change or step-size error control; the ledger states what the discrete system did, not what a physical device would do";
const NETWORK_NO_CLAIM: &str = "lumped elements with constant parameters or authored nonlinear electrical storage curves, using work-conjugate discrete-gradient stepping, and rigid_network or viscoelastic_network elements stepped by the fs-mbd network owners (their connector sources are exact at fixed geometry and their composition is second order for moving geometry); bodies have their principal axes along their own axes and their centre of mass at their origin; connector damping and branch stiffness are constant and diagonal, and a branch is a rate law in its first body's axes, not a spring with a rest geometry; no material-card resolution, contact, friction, phase change or step-size error control; the ledger states what the discrete system did, not what a physical device would do";
/// Bodies and connectors one `rigid_network` may hold: its owner sweeps the
/// connectors of the whole bank on every step.
const MAX_NETWORK_BODIES: usize = 64;
const MAX_NETWORK_CONNECTORS: usize = 256;
const ADAPTIVE_NO_CLAIM: &str = "step-doubling estimates control the declared state coordinates and source-work increments under the selected order assumption; they are local numerical estimates, not global error bounds or physical validation; dissipated heat and interface balances retain the kernel's separate ledger audits; no material-card resolution, contact, friction or phase change";

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
  adaptive          optional {order = 1, min_step_s = 1e-12, max_attempts = 32,
                    tolerances: {element.state_field: {absolute, relative}, ...}}
  elements          object: name -> element
  connections       array of connections

Adaptive stepping requires an explicit tolerance for every dynamic state
coordinate and source energy_in_j account. Absolute tolerances are positive
in that field's SI units; relative tolerances are nonnegative. Source-work
tolerances apply to each trial interval's energy increment. Derived endpoint
voltage/current and stateless last-step means are not independent accounts.
The default order is 1; choose order 2 only for a smooth second-order model.
The refined two-half-step state is retained; rejected trials leave no state
or ledger entries. step_s remains the largest interval.
The final externally capped interval may be shorter than min_step_s; an
accuracy rejection still cannot reduce that final interval further.

Elements (ports in brackets; SI units in every field name):
  mass           mass_kg, position_m = 0, velocity_m_s = 0            [p]
  spring         stiffness_n_m, cubic_n_m3 = 0, extension_m = 0       [a, b]
  damper         damping_n_s_m                                        [a, b, heat]
  force          force_n (waveform)                                   [p]
  heat_capacity  capacity_j_k, temperature_k                          [heat]
                 (and a signal output, temperature: its temperature at the
                 start of the step; it may stay unconnected)
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
For a parallel circuit, the capacitor sets the node voltage and the branches
share that voltage (current is positive into each component):
  node_capacitor capacitance_f, charge_c = 0                          [p]
  resistor_branch resistance_ohm                                      [a, b, heat]
  current_source current_a (waveform; negative supplies the node)     [p]
Source-owned reversible electrical curves, with both polarities and an
interior (0,0) knot; both coordinates must increase strictly:
  nonlinear_capacitor charge_voltage_curve = [[charge_c, voltage_v], ...],
                      charge_c = 0                                  [p]
  nonlinear_inductor flux_current_curve = [[flux_wb, current_a], ...],
                     flux_wb = 0                                    [p]
The capacitor receives current and the inductor receives voltage. Curves
refuse extrapolation; reported voltage/current are endpoint values, while
the kernel exchanges work-conjugate means over each step.
Between two domains, without loss (an ideal DC machine with a = electrical,
b = rotation and ratio in N m/A: torque ratio * i on the rotor, back-EMF
ratio * w against the loop):
  gyrator        ratio, a, b (translation | rotation | electrical)    [a, b]
Rigid bodies in three dimensions joined by viscous connectors, stepped by the
fs-mbd network owner. Vectors are [x, y, z] in world axes unless said
otherwise; a body's principal axes are its own axes and its centre of mass
is its origin:
  rigid_network  bodies = [{mass_kg, principal_inertia_kg_m2, position_m,
                   orientation = [1, 0, 0, 0] (w, x, y, z), velocity_m_s = 0,
                   angular_velocity_rad_s = 0 (in the body's axes)}, ...]
                 connectors = [{bodies = [i, j], attachments_m = [point on i,
                   point on j] (each in its body's axes; default the centres
                   of mass), damping_n_s_m = 0, rotary_damping_n_m_s_rad = 0
                   (three each, along body i's axes)}, ...]
                 gravity_m_s2 = 0, torque_ports = [] (body indices),
                 balance_tolerance = 1e-9 (the owner's allowance on a step's
                   momentum, work and speed balances, in SI units)
                 [force<i> per body, heat<e> per connector, torque<i> for
                  each body named in torque_ports]
  At most 64 bodies and 256 connectors. force<i> takes the world force on the
  centre of mass and torque<i> the torque in the body's axes; prescribe both
  with three waveforms. heat<e> is routed like a damper's heat.
  viscoelastic_network  the same fields, stepped by the fs-mbd relaxation
                 owner (which solves bodies and branches as one dense system
                 and refuses one too large for it at the first step); a
                 connector may also have
                   branches = [{stiffness_n_m = 0, rotary_stiffness_n_m_rad
                     = 0 (three each), relaxation_time_s (a Maxwell branch;
                     without it the branch never relaxes), extension_m = 0,
                     twist_rad = 0 (the branch's starting elastic stretch)}]
                   drive = \"body\" | \"world\" (an actuator: the axes its
                     held force and couple are in)
                 [as rigid_network, and drive<e> and couple<e> for each
                  driven connector: the force [N] and couple [N m] held on
                  its second attachment, three waveforms each]
  Either network may follow temperature:
                 shift = {\"arrhenius\": {activation_energy_j_mol,
                   reference_k, valid_k = [lowest, highest]}} or
                   {\"wlf\": {c1, c2_k, reference_k, valid_k}}
                 the factor tau(T) / tau_ref of the fs-material owner; it
                 multiplies every connector's damping, and in a
                 viscoelastic_network every relaxation time as well. Each
                 connector then has a signal input temperature<e> [K], read
                 once a step; a temperature outside valid_k refuses the step.
A waveform is a number, {\"sine\": {amplitude, frequency_hz, phase_rad = 0,
offset = 0}} or {\"table\": [[time_s, value], ...]} (linear, ends held).

Connections, ports written element.port:
  {\"join\": [ports]}       the ports meet at one node: a mass (inertia,
                          inductor) and the ends that share its velocity
                          (angular velocity, current), or a heat capacity and
                          the elements at its temperature
  {\"common_effort\": [ports]}
                          one effort-setting port (node_capacitor,
                          nonlinear_capacitor or voltage)
                          and branches sharing its voltage (inductor,
                          resistor_branch, current_source)
  {\"transformer\": [effort_in_port, flow_in_port], \"ratio\": r}
                          lossless signed flow ratio r; effort returns with
                          the same factor; endpoints may be in different domains
  A member of a join or common_effort may be written
  {\"port\": element.port, \"ratio\": r}: it then sits behind its own
  transformer of ratio r, so one node carries several loads at different
  ratios and in different domains (a pinion's inertia with a motor torque
  and a rack's spring, damper and load at the pinion radius). Power ports
  only; the node port takes no ratio.
  {\"route\": port, \"to\": [[port, fraction], ...],
   \"rest\": {\"sink\": name, \"temperature_k\": T},
   \"temperature\": \"sink-equivalent\" | {\"fixed_k\": T}}
                          send the heat leaving a port by fraction
  {\"signal\": port, \"to\": [ports]}
                          copy a signal output (heat_capacity.temperature) to
                          signal inputs (a network's temperature<e>); it
                          carries no energy
  {\"prescribe\": port, <quantity>: waveform}
                          hold a port's input from outside: velocity_m_s for a
                          spring, damper or force end (0 anchors it), force_n
                          for a mass; angular_velocity_rad_s and torque_n_m
                          about a shaft; current_a and voltage_v in a loop (0
                          current grounds a series element's far end); voltage_v
                          grounds a resistor_branch end; current_a into a
                          node_capacitor or nonlinear_capacitor increases its charge;
                          temperature_k for a damper, resistor, conductance or
                          heater, heat_rate_w for a heat capacity; for a
                          network of rigid bodies, force_n or torque_n_m as
                          an array of three waveforms ([0, 0, 0] leaves a
                          body free or an actuator idle), and temperature_k
                          for a temperature<e> input held from outside

Each element's port a is the end at mass A and b the end at mass B; the
extension of a spring is x_B - x_A less its unloaded length; a series
capacitor's charge_c is the same difference, so one charged by a current
entering at a reads negative. A node_capacitor instead stores entering charge. Every port must be connected. See examples/system/.
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
        Err(RunError::Adaptive(error)) => {
            let code = if matches!(error, AdaptiveError::Step(SystemError::Cancelled)) {
                exit::CANCELLED
            } else {
                exit::REFUSED
            };
            let message = format!(
                "step {} at t = {} s: {error}",
                device.system.steps() + 1,
                device.system.time_s()
            );
            refusal(code, "system-adaptive", &message, json)
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
    Adaptive(AdaptiveError),
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

/// A lumped power or thermal element, or a lossless domain converter.
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
    /// Charge storage setting the voltage of a parallel electrical node.
    NodeCapacitor,
    /// Current from a voltage difference, with a heat port.
    ResistorBranch,
    /// A prescribed current, positive into the source.
    CurrentSource,
    /// Charge storage governed by a source-owned V(q) curve.
    NonlinearCapacitor,
    /// Flux storage governed by a source-owned I(flux) curve.
    NonlinearInductor,
    Gyrator,
    /// Rigid bodies and viscous connectors, stepped by their fs-mbd owner.
    RigidNetwork,
    /// The same with Maxwell branches and actuators on the connectors.
    ViscoelasticNetwork,
    HeatCapacity,
    Conductance,
    Heater,
}

impl Kind {
    const TABLE: [(&'static str, Self); 23] = [
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
        ("node_capacitor", Self::NodeCapacitor),
        ("resistor_branch", Self::ResistorBranch),
        ("current_source", Self::CurrentSource),
        ("nonlinear_capacitor", Self::NonlinearCapacitor),
        ("nonlinear_inductor", Self::NonlinearInductor),
        ("gyrator", Self::Gyrator),
        ("rigid_network", Self::RigidNetwork),
        ("viscoelastic_network", Self::ViscoelasticNetwork),
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
    /// Convert a node capacitor's endpoint charge to its endpoint voltage.
    capacitance_f: Option<f64>,
    /// Original laws for endpoint probes; prepared data are shared on clone.
    charge_voltage_curve: Option<CapacitorChargeCurve>,
    flux_current_curve: Option<InductorFluxCurve>,
    /// The bodies of a network, for reading velocities from momenta, and for
    /// each of its connectors the number of its Maxwell branches.
    bodies: Vec<MassProperties>,
    branches: Vec<usize>,
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
    adaptive: Option<AdaptiveSettings>,
}

#[derive(Clone, Copy)]
enum Coordinate {
    Position,
    Velocity,
    Scalar,
    SourceWork,
}

struct ErrorAccount {
    node: NodeId,
    coordinate: Coordinate,
    tolerance: ErrorTolerance,
}

impl ErrorAccount {
    fn read(&self, checkpoint: &Checkpoint) -> f64 {
        match self.coordinate {
            Coordinate::Position => checkpoint.state::<Motion>(self.node)
                .map_or(f64::NAN, |state| state.position_m),
            Coordinate::Velocity => checkpoint.state::<Motion>(self.node)
                .map_or(f64::NAN, |state| state.velocity_m_s),
            Coordinate::Scalar => checkpoint.state::<f64>(self.node).copied().unwrap_or(f64::NAN),
            Coordinate::SourceWork => checkpoint.state::<SourceAccount>(self.node)
                .map_or(f64::NAN, |state| state.energy_in_j),
        }
    }
}

struct AdaptiveSettings {
    config: AdaptiveConfig,
    accounts: Vec<ErrorAccount>,
    attempts: u64,
    trial_steps: u64,
    last_error_norm: f64,
    max_error_norm: f64,
    suggested_dt_s: f64,
}

impl AdaptiveSettings {
    fn metric(&self, initial: &Checkpoint, coarse: &Checkpoint, refined: &Checkpoint) -> f64 {
        let mut norm = ErrorNorm::new();
        for account in &self.accounts {
            let (start, full, half) = (account.read(initial), account.read(coarse), account.read(refined));
            if matches!(account.coordinate, Coordinate::SourceWork) {
                // Work is controlled over this interval, so a large cumulative
                // source account cannot loosen a later interval's tolerance.
                norm.include(0.0, full - start, half - start, account.tolerance);
            } else {
                norm.include(start, full, half, account.tolerance);
            }
        }
        norm.value()
    }

    fn record(&mut self, report: &AdaptiveReport) {
        self.attempts += u64::from(report.attempts);
        self.trial_steps += report.trial_steps;
        self.last_error_norm = report.error_norm;
        self.max_error_norm = self.max_error_norm.max(report.error_norm);
        self.suggested_dt_s = report.suggested_dt_s;
    }
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

fn storage_knots(
    definition: &J,
    key: &str,
    coordinates: &str,
    maximum: usize,
    what: &str,
) -> Result<Vec<(f64, f64)>, String> {
    let rows = definition
        .get(key)
        .and_then(J::as_array)
        .ok_or_else(|| format!("{what}.{key}: expected an array of {coordinates}"))?;
    if !(3..=maximum).contains(&rows.len()) {
        return Err(format!("{what}.{key}: expected 3..={maximum} knots"));
    }
    rows.iter()
        .map(|row| match row.as_array() {
            Some([coordinate, conjugate]) => coordinate.as_f64().zip(conjugate.as_f64()),
            _ => None,
        })
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| format!("{what}.{key}: every knot is {coordinates}"))
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

/// Three finite numbers.
fn point(value: &J, what: &str) -> Result<[f64; 3], String> {
    let numbers = match value.as_array() {
        Some([x, y, z]) => x.as_f64().zip(y.as_f64()).zip(z.as_f64()),
        _ => None,
    };
    match numbers {
        Some(((x, y), z)) if x.is_finite() && y.is_finite() && z.is_finite() => Ok([x, y, z]),
        _ => Err(format!("{what}: expected three finite numbers")),
    }
}

/// The three numbers of field `key`, or `default` where the field is absent.
fn triple(
    value: &J,
    key: &str,
    default: Option<[f64; 3]>,
    what: &str,
) -> Result<[f64; 3], String> {
    match value.get(key) {
        Some(field) => point(field, &format!("{what}.{key}")),
        None => default.ok_or_else(|| format!("{what}: missing {key}")),
    }
}

fn vec3(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

/// A whole number from 0 to `count - 1`.
fn index_below(value: &J, count: usize) -> Option<usize> {
    value
        .as_f64()
        .filter(|v| v.fract() == 0.0 && *v >= 0.0 && *v < count as f64)
        .map(|v| v as usize)
}

/// One body of a rigid network: its mass properties and its starting state.
fn network_body(row: &J, what: &str) -> Result<(MassProperties, RigidBodyState), Refusal> {
    fields(
        row,
        &[
            "mass_kg",
            "principal_inertia_kg_m2",
            "position_m",
            "orientation",
            "velocity_m_s",
            "angular_velocity_rad_s",
        ],
        what,
    )?;
    let refused = |e: fs_mbd::DynamicsError| Refusal::Input(format!("{what}: {e}"));
    let properties = MassProperties::new(
        positive(row, "mass_kg", what)?,
        Vec3::ZERO,
        vec3(triple(row, "principal_inertia_kg_m2", None, what)?),
    )
    .map_err(refused)?;
    let orientation = match row.get("orientation") {
        None => UnitQuaternion::IDENTITY,
        Some(field) => {
            let numbers = match field.as_array() {
                Some([w, x, y, z]) => w.as_f64().zip(x.as_f64()).zip(y.as_f64().zip(z.as_f64())),
                _ => None,
            };
            let ((w, x), (y, z)) = numbers
                .ok_or_else(|| format!("{what}: orientation is a quaternion [w, x, y, z]"))?;
            UnitQuaternion::new(w, x, y, z).map_err(refused)?
        }
    };
    let pose =
        Pose::new(vec3(triple(row, "position_m", None, what)?), orientation).map_err(refused)?;
    let velocity = triple(row, "velocity_m_s", Some([0.0; 3]), what)?;
    let mass = properties.mass();
    let momentum = Vec3::new(mass * velocity[0], mass * velocity[1], mass * velocity[2]);
    let spin = properties
        .angular_momentum_body_checked(vec3(triple(
            row,
            "angular_velocity_rad_s",
            Some([0.0; 3]),
            what,
        )?))
        .map_err(refused)?;
    let state = RigidBodyState::new(pose, momentum, spin).map_err(refused)?;
    Ok((properties, state))
}

/// One viscous connector of a network of `bodies` bodies. `extra` names the
/// further fields of the row that the caller reads itself.
fn network_connector(
    row: &J,
    bodies: usize,
    extra: &[&str],
    what: &str,
) -> Result<ViscousNetworkLink, Refusal> {
    let mut allowed = vec!["bodies", "attachments_m", "damping_n_s_m", "rotary_damping_n_m_s_rad"];
    allowed.extend_from_slice(extra);
    fields(row, &allowed, what)?;
    let ends = match row.get("bodies").and_then(J::as_array) {
        Some([a, b]) => index_below(a, bodies).zip(index_below(b, bodies)),
        _ => None,
    };
    let (first, second) = ends
        .filter(|(a, b)| a != b)
        .ok_or_else(|| format!("{what}: bodies is two different indices into bodies"))?;
    let arms_body_m = match row.get("attachments_m") {
        None => [Vec3::ZERO; 2],
        Some(field) => match field.as_array() {
            Some([a, b]) => {
                let inner = format!("{what}.attachments_m");
                [vec3(point(a, &inner)?), vec3(point(b, &inner)?)]
            }
            _ => return Err(format!("{what}: attachments_m is two points").into()),
        },
    };
    let linear = triple(row, "damping_n_s_m", Some([0.0; 3]), what)?;
    let rotary = triple(row, "rotary_damping_n_m_s_rad", Some([0.0; 3]), what)?;
    let mut coefficients = [[0.0; 6]; 6];
    for axis in 0..3 {
        coefficients[axis][axis] = linear[axis];
        coefficients[3 + axis][3 + axis] = rotary[axis];
    }
    let damping = PassiveBushingMatrix::try_new(
        coefficients,
        BushingCoordinateScales {
            translation_m: 1.0,
            rotation_rad: 1.0,
        },
    )
    .map_err(|e| format!("{what}: damping is not a passive form: {e:?}"))?;
    Ok(ViscousNetworkLink {
        bodies: [first, second],
        arms_body_m,
        damping,
    })
}

/// A scene's temperature shift `a_T = tau(T) / tau_ref`, evaluated by the
/// fs-material owner, offered to the two fs-mbd network Components: it
/// scales a viscous network's damping forms, and a viscoelastic network's
/// Kelvin viscosity and relaxation times together.
#[derive(Debug)]
struct SceneShift(TemperatureShift);

impl SceneShift {
    fn at(&self, temperature_k: f64) -> Result<f64, String> {
        self.0
            .evaluate(temperature_k)
            .map(|state| state.factor())
            .map_err(|e| e.to_string())
    }
}

impl ShiftLaw for SceneShift {
    fn factor(&self, _connector: usize, temperature_k: f64) -> Result<f64, String> {
        self.at(temperature_k)
    }
}

impl DampingLaw for SceneShift {
    fn at_temperature(
        &self,
        _connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> Result<PassiveBushingMatrix, String> {
        reference
            .scaled(self.at(temperature_k)?)
            .map_err(|e| format!("{e:?}"))
    }
}

/// The `shift` field of a network: an Arrhenius or WLF calibration with its
/// reference temperature and the band it is valid in.
fn scene_shift(definition: &J, what: &str) -> Result<Option<TemperatureShift>, Refusal> {
    let Some(shift) = definition.get("shift") else {
        return Ok(None);
    };
    let what = format!("{what}.shift");
    let (model, body) = match shift.as_object() {
        Some([(key, body)]) if key == "arrhenius" => {
            fields(body, &["activation_energy_j_mol", "reference_k", "valid_k"], &what)?;
            let activation_energy_j_mol = number(body, "activation_energy_j_mol", &what)?;
            (TemperatureShiftModel::Arrhenius { activation_energy_j_mol }, body)
        }
        Some([(key, body)]) if key == "wlf" => {
            fields(body, &["c1", "c2_k", "reference_k", "valid_k"], &what)?;
            let (c1, c2_k) = (number(body, "c1", &what)?, number(body, "c2_k", &what)?);
            (TemperatureShiftModel::WilliamsLandelFerry { c1, c2_k }, body)
        }
        _ => {
            return Err(format!("{what}: expected {{\"arrhenius\": ...}} or {{\"wlf\": ...}}").into());
        }
    };
    let valid = match body.get("valid_k").and_then(J::as_array) {
        Some([lowest, highest]) => lowest.as_f64().zip(highest.as_f64()),
        _ => None,
    };
    let (lowest, highest) =
        valid.ok_or_else(|| format!("{what}: valid_k is [lowest, highest] in kelvin"))?;
    TemperatureShift::new(model, number(body, "reference_k", &what)?, [lowest, highest])
        .map(Some)
        .map_err(|e| Refusal::Input(format!("{what}: {e}")))
}

/// What the two network elements share in a scene: gravity, the bodies and
/// their starting states, the connector rows, the owner's balance controls,
/// the bodies that get a torque port and the temperature shift, if any.
struct NetworkBase<'a> {
    gravity: Gravity,
    properties: Vec<MassProperties>,
    states: Vec<RigidBodyState>,
    connectors: &'a [J],
    source: ViscousConnectorControls,
    torques: Vec<usize>,
    shift: Option<TemperatureShift>,
}

fn network_base<'a>(definition: &'a J, what: &str) -> Result<NetworkBase<'a>, Refusal> {
    fields(
        definition,
        &[
            "kind",
            "bodies",
            "connectors",
            "gravity_m_s2",
            "torque_ports",
            "balance_tolerance",
            "shift",
        ],
        what,
    )?;
    let shift = scene_shift(definition, what)?;
    let gravity = Gravity::new(vec3(triple(definition, "gravity_m_s2", Some([0.0; 3]), what)?))
        .map_err(|e| format!("{what}.gravity_m_s2: {e}"))?;
    let rows = definition
        .get("bodies")
        .and_then(J::as_array)
        .filter(|rows| (1..=MAX_NETWORK_BODIES).contains(&rows.len()))
        .ok_or_else(|| format!("{what}: bodies is an array of 1 to {MAX_NETWORK_BODIES} bodies"))?;
    let mut properties = Vec::with_capacity(rows.len());
    let mut states = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        let (body, state) = network_body(row, &format!("{what}.bodies[{index}]"))?;
        properties.push(body);
        states.push(state);
    }
    let connectors = definition
        .get("connectors")
        .and_then(J::as_array)
        .filter(|rows| rows.len() <= MAX_NETWORK_CONNECTORS)
        .ok_or_else(|| {
            format!("{what}: connectors is an array of at most {MAX_NETWORK_CONNECTORS} connectors")
        })?;
    let tolerance = match definition.get("balance_tolerance") {
        None => 1.0e-9,
        Some(_) => positive(definition, "balance_tolerance", what)?,
    };
    let source = ViscousConnectorControls {
        impulse: InternalWrenchImpulseControls {
            angular_origin_world_m: Vec3::ZERO,
            linear_tolerance_ns: tolerance,
            angular_tolerance_nms: tolerance,
            work_tolerance_j: tolerance,
        },
        linear_speed_tolerance_m_s: tolerance,
        angular_speed_tolerance_rad_s: tolerance,
    };
    let mut torques = Vec::new();
    if let Some(list) = definition.get("torque_ports") {
        let rows = list
            .as_array()
            .ok_or_else(|| format!("{what}: torque_ports is an array of body indices"))?;
        for row in rows {
            torques.push(
                index_below(row, properties.len())
                    .ok_or_else(|| format!("{what}: torque_ports is an array of body indices"))?,
            );
        }
    }
    Ok(NetworkBase {
        gravity,
        properties,
        states,
        connectors,
        source,
        torques,
        shift,
    })
}

/// A network element lowered into the builder: its node and ports, its
/// bodies, and for each connector the number of its Maxwell branches.
struct NetworkParts {
    node: NodeId,
    ports: Vec<PortSpec>,
    bodies: Vec<MassProperties>,
    branches: Vec<usize>,
}

/// A `rigid_network`: rigid bodies joined by viscous connectors, lowered to
/// the fs-mbd network Component.
fn rigid_network(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<NetworkParts, Refusal> {
    let base = network_base(definition, what)?;
    let links = base
        .connectors
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let inner = format!("{what}.connectors[{index}]");
            network_connector(row, base.properties.len(), &[], &inner)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let branches = vec![0; links.len()];
    let controls = ViscousNetworkMotionControls {
        source: base.source,
        source_substeps: 1,
        limits: ViscousNetworkLimits {
            maximum_bodies: base.properties.len(),
            maximum_connectors: links.len(),
            maximum_pair_relaxations: 1_000_000,
            maximum_live_bytes: 64 << 20,
        },
    };
    // With a shift, every connector's damping form follows its own
    // `temperature<e>` signal.
    let built = match base.shift {
        None => ViscousNetwork::new(base.properties.clone(), links, base.gravity, controls),
        Some(shift) => ViscousNetwork::temperature_dependent(
            base.properties.clone(),
            links,
            base.gravity,
            controls,
            Rc::new(SceneShift(shift)),
        ),
    };
    let mut network = built.map_err(|e| Refusal::Input(format!("{what}: {e}")))?;
    for &body in &base.torques {
        network = network
            .torqued(body)
            .map_err(|e| Refusal::Input(format!("{what}.torque_ports: {e}")))?;
    }
    let ports = network.ports().to_vec();
    let node = builder
        .add(name, network, base.states)
        .map_err(|e| Refusal::Admission(e.to_string()))?;
    Ok(NetworkParts {
        node,
        ports,
        bodies: base.properties,
        branches,
    })
}

/// The Maxwell branches of one connector and their starting memory. A branch
/// of stiffness `k` on an axis has the factor `sqrt(k)` there, and an
/// elastic extension `x` of it is the memory amplitude `sqrt(k) x`, whose
/// stored energy is `k x^2 / 2`.
fn network_branches(
    row: &J,
    what: &str,
) -> Result<(Vec<SpatialRelaxationBranch>, Vec<[f64; 6]>), Refusal> {
    let Some(list) = row.get("branches") else {
        return Ok((Vec::new(), Vec::new()));
    };
    let rows = list
        .as_array()
        .ok_or_else(|| format!("{what}: branches is an array of branches"))?;
    let mut branches = Vec::with_capacity(rows.len());
    let mut memory = Vec::with_capacity(rows.len());
    for (index, branch) in rows.iter().enumerate() {
        let inner = format!("{what}.branches[{index}]");
        fields(
            branch,
            &[
                "stiffness_n_m",
                "rotary_stiffness_n_m_rad",
                "relaxation_time_s",
                "extension_m",
                "twist_rad",
            ],
            &inner,
        )?;
        let linear = triple(branch, "stiffness_n_m", Some([0.0; 3]), &inner)?;
        let rotary = triple(branch, "rotary_stiffness_n_m_rad", Some([0.0; 3]), &inner)?;
        let extension = triple(branch, "extension_m", Some([0.0; 3]), &inner)?;
        let twist = triple(branch, "twist_rad", Some([0.0; 3]), &inner)?;
        let mut factors = [[0.0; 6]; 6];
        let mut amplitude = [0.0; 6];
        for axis in 0..3 {
            for (slot, stiffness, stretch) in [
                (axis, linear[axis], extension[axis]),
                (3 + axis, rotary[axis], twist[axis]),
            ] {
                if stiffness < 0.0 {
                    return Err(format!("{inner}: a stiffness is never negative").into());
                }
                if stiffness == 0.0 && stretch != 0.0 {
                    return Err(format!(
                        "{inner}: an extension or twist needs a stiffness on the same axis"
                    )
                    .into());
                }
                factors[slot][slot] = stiffness.sqrt();
                amplitude[slot] = stiffness.sqrt() * stretch;
            }
        }
        let built = match branch.get("relaxation_time_s") {
            Some(_) => SpatialRelaxationBranch::from_gram_rows(
                factors,
                positive(branch, "relaxation_time_s", &inner)?,
            ),
            None => SpatialRelaxationBranch::from_equilibrium_gram_rows(factors),
        }
        .map_err(|e| format!("{inner}: {e:?}"))?;
        branches.push(built);
        memory.push(amplitude);
    }
    Ok((branches, memory))
}

/// A `viscoelastic_network`: rigid bodies joined by connectors with viscous
/// damping, Maxwell branches and, on request, an actuator, lowered to the
/// fs-mbd relaxation-network Component.
fn viscoelastic_network(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<NetworkParts, Refusal> {
    let base = network_base(definition, what)?;
    let mut links = Vec::with_capacity(base.connectors.len());
    let mut memory = Vec::with_capacity(base.connectors.len());
    let mut drives = Vec::new();
    for (index, row) in base.connectors.iter().enumerate() {
        let inner = format!("{what}.connectors[{index}]");
        let viscous =
            network_connector(row, base.properties.len(), &["branches", "drive"], &inner)?;
        let (branches, amplitudes) = network_branches(row, &inner)?;
        match row.get("drive").map(J::as_str) {
            None => {}
            Some(Some("body")) => drives.push((index, DriveFrame::MaterialFirst)),
            Some(Some("world")) => drives.push((index, DriveFrame::World)),
            Some(_) => return Err(format!("{inner}: drive is \"body\" or \"world\"").into()),
        }
        links.push(RelaxationNetworkLink { viscous, branches });
        memory.push(amplitudes);
    }
    let branches: Vec<usize> = memory.iter().map(Vec::len).collect();
    let controls = RelaxationNetworkMotionControls {
        source: base.source,
        limits: RelaxationNetworkLimits {
            source: DrivenViscousNetworkLimits {
                maximum_bodies: base.properties.len(),
                maximum_connectors: links.len(),
                maximum_live_bytes: 128_000_000,
                maximum_preparation_terms: 10_000_000_000,
            },
            maximum_branches: branches.iter().sum(),
        },
    };
    let refused = |e: fs_mbd::component::relaxation::RelaxationNetworkBuildError| {
        Refusal::Input(format!("{what}: {e}"))
    };
    let mut network =
        RelaxationNetwork::new(base.properties.clone(), links, base.gravity, controls)
            .map_err(refused)?;
    for (connector, frame) in drives {
        network = network.driven(connector, frame).map_err(refused)?;
    }
    // With a shift, every connector's Kelvin viscosity and relaxation
    // times follow its own `temperature<e>` signal.
    if let Some(shift) = base.shift {
        network = network.shifted(Rc::new(SceneShift(shift))).map_err(refused)?;
    }
    for &body in &base.torques {
        network = network.torqued(body).map_err(refused)?;
    }
    let ports = network.ports().to_vec();
    let state = RelaxationState {
        bodies: base.states,
        memory,
    };
    let node = builder
        .add(name, network, state)
        .map_err(|e| Refusal::Admission(e.to_string()))?;
    Ok(NetworkParts {
        node,
        ports,
        bodies: base.properties,
        branches,
    })
}

fn element(builder: &mut SystemBuilder, name: &str, definition: &J) -> Result<Element, Refusal> {
    let what = format!("elements.{name}");
    let kind_name = definition
        .str_field("kind")
        .ok_or_else(|| format!("{what}: missing kind"))?;
    let kind = Kind::parse(kind_name).ok_or_else(|| format!("{what}: unknown kind {kind_name:?}"))?;
    let invalid = |e: fs_system::elements::ElementError| Refusal::Input(format!("{what}: {e}"));
    let mut temperature_k = None;
    let mut capacitance_f = None;
    let mut charge_voltage_curve = None;
    let mut flux_current_curve = None;
    let mut bodies = Vec::new();
    let mut branches = Vec::new();
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
        Kind::NodeCapacitor => {
            fields(definition, &["kind", "capacitance_f", "charge_c"], &what)?;
            let capacity = positive(definition, "capacitance_f", &what)?;
            let charge = number_or(definition, "charge_c", 0.0, &what)?;
            if !charge.is_finite() {
                return Err(format!("{what}: charge_c must be finite").into());
            }
            let storage = Storage::in_domain(Domain::Electrical.port_kind(), capacity)
                .map_err(invalid)?;
            capacitance_f = Some(capacity);
            add(builder, name, storage, charge)?
        }
        Kind::ResistorBranch => {
            fields(definition, &["kind", "resistance_ohm"], &what)?;
            let resistance = positive(definition, "resistance_ohm", &what)?;
            let branch = Admittance::in_domain(
                Domain::Electrical.port_kind(),
                1.0 / resistance,
            )
            .map_err(invalid)?;
            add(builder, name, branch, ())?
        }
        Kind::CurrentSource => {
            fields(definition, &["kind", "current_a"], &what)?;
            let wave = definition
                .get("current_a")
                .ok_or_else(|| format!("{what}: missing current_a"))?;
            let source = FlowSource::in_domain(
                Domain::Electrical.port_kind(),
                waveform(wave, &format!("{what}.current_a"))?,
            );
            add(builder, name, source, SourceAccount::default())?
        }
        Kind::NonlinearCapacitor => {
            fields(definition, &["kind", "charge_voltage_curve", "charge_c"], &what)?;
            let curve = CapacitorChargeCurve::try_new(storage_knots(
                definition, "charge_voltage_curve", "[charge_c, voltage_v]",
                MAX_CAPACITOR_KNOTS, &what,
            )?)
            .map_err(|error| format!("{what}.charge_voltage_curve: {error:?}"))?;
            let charge = number_or(definition, "charge_c", 0.0, &what)?;
            curve.voltage_v(charge)
                .map_err(|error| format!("{what}.charge_c: {error:?}"))?;
            charge_voltage_curve = Some(curve.clone());
            add(builder, name, NonlinearCapacitor::new(curve), charge)?
        }
        Kind::NonlinearInductor => {
            fields(definition, &["kind", "flux_current_curve", "flux_wb"], &what)?;
            let curve = InductorFluxCurve::try_new(storage_knots(
                definition, "flux_current_curve", "[flux_wb, current_a]",
                MAX_INDUCTOR_KNOTS, &what,
            )?)
            .map_err(|error| format!("{what}.flux_current_curve: {error:?}"))?;
            let flux = number_or(definition, "flux_wb", 0.0, &what)?;
            curve.current_a(flux)
                .map_err(|error| format!("{what}.flux_wb: {error:?}"))?;
            flux_current_curve = Some(curve.clone());
            add(builder, name, NonlinearInductor::new(curve), flux)?
        }
        Kind::Gyrator => gyrator(builder, name, definition, &what)?,
        Kind::RigidNetwork | Kind::ViscoelasticNetwork => {
            let parts = if kind == Kind::RigidNetwork {
                rigid_network(builder, name, definition, &what)?
            } else {
                viscoelastic_network(builder, name, definition, &what)?
            };
            bodies = parts.bodies;
            branches = parts.branches;
            (parts.node, parts.ports)
        }
        Kind::HeatCapacity => {
            fields(definition, &["kind", "capacity_j_k", "temperature_k"], &what)?;
            let start = positive(definition, "temperature_k", &what)?;
            temperature_k = Some(start);
            // Every heat capacity of a scene has the signal output
            // `temperature`; unconnected, it changes nothing.
            let body = HeatCapacity::new(number(definition, "capacity_j_k", &what)?)
                .map_err(invalid)?
                .publishing();
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
        capacitance_f,
        charge_voltage_curve,
        flux_current_curve,
        bodies,
        branches,
    })
}

/// A port named in the scene, resolved.
#[derive(Clone, Copy)]
struct Resolved {
    element: usize,
    port: PortRef,
    carrier: Carrier,
    role: PortRole,
    /// Coefficients the port carries each way.
    width: usize,
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
        width: spec.width(),
    })
}

fn port_text(value: &J, what: &str) -> Result<String, String> {
    value
        .as_str()
        .map(str::to_owned)
        .ok_or_else(|| format!("{what}: a port is a string element.port"))
}

/// One entry of a junction: `"element.port"`, or
/// `{"port": "element.port", "ratio": r}` for a member that sits behind an
/// ideal transformer of that ratio.
fn member(value: &J, what: &str) -> Result<(String, Option<f64>), String> {
    if value.as_object().is_some() {
        fields(value, &["port", "ratio"], what)?;
        let port = value
            .get("port")
            .ok_or_else(|| format!("{what}: a member with a ratio names its port"))?;
        Ok((port_text(port, what)?, Some(number(value, "ratio", what)?)))
    } else {
        Ok((port_text(value, what)?, None))
    }
}

/// Starting values for temperature inputs inside loops.
type Guesses = Vec<(PortRef, f64)>;

fn join(
    builder: &mut SystemBuilder,
    elements: &[Element],
    ports: &[(Resolved, Option<f64>)],
    guesses: &mut Guesses,
    common_effort: bool,
    what: &str,
) -> Result<(), Refusal> {
    let ratios: Vec<Option<f64>> = ports.iter().map(|(_, ratio)| *ratio).collect();
    let ports: Vec<&Resolved> = ports.iter().map(|(port, _)| port).collect();
    if ports.len() < 2 {
        return Err(format!("{what}: a join needs at least two ports").into());
    }
    // A legacy join uses a flow-setting power node or a temperature-setting
    // heat node. The explicit common_effort form also admits a power node
    // that sets effort, such as a capacitor setting the voltage of branches.
    let is_node = |port: &Resolved| {
        if common_effort {
            port.role == PortRole::FlowIn
                && matches!(port.carrier, Carrier::Power(_) | Carrier::Heat)
        } else {
            matches!(
                (port.carrier, port.role),
                (Carrier::Power(_), PortRole::EffortIn) | (Carrier::Heat, PortRole::FlowIn)
            )
        }
    };
    let nodes = ports.iter().filter(|port| is_node(port)).count();
    if nodes != 1 {
        if common_effort {
            return Err(format!(
                "{what}: common_effort has exactly one effort-setting port \
                 (such as node_capacitor.p) and the branches sharing its effort; \
                 this one has {nodes}"
            )
            .into());
        }
        return Err(format!(
            "{what}: a join has exactly one node port (a mass, inertia, inductor or heat \
             capacity) and the ports that meet at it; this one has {nodes}"
        )
        .into());
    }
    let setter = ports.iter().position(|port| is_node(port)).unwrap_or_default();
    let setter_role = ports[setter].role;
    if ratios[setter].is_some() {
        return Err(format!(
            "{what}: the node port sets what the others share; a ratio belongs to a member"
        )
        .into());
    }
    let members: Vec<PortRef> = ports
        .iter()
        .enumerate()
        .filter(|&(index, _)| index != setter)
        .map(|(_, p)| p.port)
        .collect();
    let result = if ratios.iter().any(Option::is_some) {
        // Members behind ratios: power ports of any domain, each through
        // its own ideal transformer.
        let geared: Vec<(PortRef, f64)> = ports
            .iter()
            .zip(&ratios)
            .enumerate()
            .filter(|&(index, _)| index != setter)
            .map(|(_, (p, ratio))| (p.port, ratio.unwrap_or(1.0)))
            .collect();
        if setter_role == PortRole::EffortIn {
            builder.common_flow_ratios(ports[setter].port, &geared)
        } else {
            builder.common_effort_ratios(ports[setter].port, &geared)
        }
    } else if setter_role == PortRole::EffortIn {
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

fn transformer(
    builder: &mut SystemBuilder,
    elements: &[Element],
    connection: &J,
    what: &str,
) -> Result<(), Refusal> {
    fields(connection, &["transformer", "ratio"], what)?;
    let Some([a, b]) = connection.get("transformer").and_then(J::as_array) else {
        return Err(format!("{what}: transformer is [effort_in_port, flow_in_port]").into());
    };
    let a = resolve(builder, elements, &port_text(a, what)?)?;
    let b = resolve(builder, elements, &port_text(b, what)?)?;
    builder
        .transformer(a.port, b.port, number(connection, "ratio", what)?)
        .map_err(|e| Refusal::Admission(format!("{what}: {e}")))
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

/// `{"signal": output, "to": [inputs]}`: copy a signal output (a heat
/// capacity's start-of-step temperature) to signal inputs (the temperatures
/// a network's shift is read at).
fn signal(
    builder: &mut SystemBuilder,
    elements: &[Element],
    connection: &J,
    guesses: &mut Guesses,
    what: &str,
) -> Result<(), Refusal> {
    fields(connection, &["signal", "to"], what)?;
    let source = resolve(
        builder,
        elements,
        &port_text(connection.get("signal").unwrap_or(&J::Null), what)?,
    )?;
    let targets = connection
        .get("to")
        .and_then(J::as_array)
        .filter(|rows| !rows.is_empty())
        .ok_or_else(|| format!("{what}: to is an array of the ports that read the signal"))?;
    for target in targets {
        let target = resolve(builder, elements, &port_text(target, what)?)?;
        builder
            .signal(source.port, target.port)
            .map_err(|e| Refusal::Admission(format!("{what}: {e}")))?;
        // The kernel solves a signal with the rest of its block and starts
        // from the publisher's temperature.
        if let Some(kelvin) = elements[source.element].temperature_k {
            guesses.push((target.port, kelvin));
        }
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
        // The only signal inputs of a scene are the temperatures a network's
        // shift is read at.
        (Carrier::Heat, PortRole::EffortIn, _) | (Carrier::Signal, PortRole::SignalIn, _) => {
            "temperature_k"
        }
        (Carrier::Heat, PortRole::FlowIn, _) => "heat_rate_w",
        (Carrier::Power(_), PortRole::EffortIn, Some(domain)) => domain.names().effort,
        (Carrier::Power(_), PortRole::FlowIn, Some(domain)) => domain.names().flow,
        _ => return Err(format!("{what}: this port cannot be prescribed").into()),
    };
    fields(connection, &["prescribe", quantity], what)
        .map_err(|_| format!("{what}: this port takes exactly one quantity, {quantity}"))?;
    let given = connection
        .get(quantity)
        .ok_or_else(|| format!("{what}: missing {quantity}"))?;
    // A scalar port takes one waveform; a port of several components (a
    // rigid network's force or torque) takes one per component.
    let waves = if target.width == 1 {
        vec![waveform(given, &format!("{what}.{quantity}"))?]
    } else {
        given
            .as_array()
            .filter(|rows| rows.len() == target.width)
            .ok_or_else(|| {
                format!(
                    "{what}: {quantity} for this port is an array of {} waveforms",
                    target.width
                )
            })?
            .iter()
            .enumerate()
            .map(|(index, row)| waveform(row, &format!("{what}.{quantity}[{index}]")))
            .collect::<Result<Vec<_>, _>>()?
    };
    if quantity == "temperature_k" && !waves.iter().all(|wave| wave.lowest() > 0.0) {
        return Err(format!("{what}: temperature_k must stay positive").into());
    }
    // A series element end carried by a node receives minus the node's
    // flow. A node capacitor instead receives current into its own stored
    // charge, so its prescribed current keeps the component's positive sign.
    let prescribed_flow =
        matches!((target.carrier, target.role), (Carrier::Power(_), PortRole::FlowIn))
            && !matches!(elements[target.element].kind, Kind::NodeCapacitor | Kind::NonlinearCapacitor);
    let sign = if prescribed_flow { -1.0 } else { 1.0 };
    builder
        .boundary(
            target.port,
            Box::new(move |t0, dt, values| {
                for (value, wave) in values.iter_mut().zip(&waves) {
                    *value = sign * wave.mean(t0, dt);
                }
            }),
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

fn adaptive(definition: &J, elements: &[Element], step_s: f64) -> Result<AdaptiveSettings, String> {
    fields(definition, &["order", "min_step_s", "max_attempts", "tolerances"], "adaptive")?;
    let defaults = AdaptiveConfig::default();
    let order = number_or(definition, "order", f64::from(defaults.order), "adaptive")?;
    let attempts = number_or(definition, "max_attempts", f64::from(defaults.max_attempts), "adaptive")?;
    let minimum = number_or(definition, "min_step_s", defaults.min_step_s, "adaptive")?;
    if !(1.0..=32.0).contains(&order) || order.fract() != 0.0 {
        return Err("adaptive.order must be a whole number from 1 to 32".into());
    }
    if !(1.0..=f64::from(u32::MAX)).contains(&attempts) || attempts.fract() != 0.0 {
        return Err("adaptive.max_attempts must be a positive whole number fitting u32".into());
    }
    if !(minimum.is_finite() && minimum > 0.0 && minimum <= step_s) {
        return Err("adaptive.min_step_s must be positive, finite and no larger than step_s".into());
    }
    let tolerances = definition.get("tolerances").and_then(J::as_object)
        .ok_or("adaptive.tolerances must name every dynamic state coordinate and source-work account")?;
    let mut expected = Vec::new();
    for element in elements {
        let coordinates = match element.kind {
            Kind::Inertance(domain) => vec![
                (domain.names().displacement, Coordinate::Position),
                (domain.names().flow, Coordinate::Velocity),
            ],
            Kind::Compliance(domain) => vec![(domain.names().extension, Coordinate::Scalar)],
            Kind::NodeCapacitor | Kind::NonlinearCapacitor => vec![("charge_c", Coordinate::Scalar)],
            Kind::NonlinearInductor => vec![("flux_wb", Coordinate::Scalar)],
            Kind::HeatCapacity => vec![("temperature_k", Coordinate::Scalar)],
            Kind::Source(_) | Kind::CurrentSource | Kind::Heater => vec![("energy_in_j", Coordinate::SourceWork)],
            Kind::Resistance(_) | Kind::ResistorBranch | Kind::Gyrator | Kind::Conductance => Vec::new(),
            Kind::RigidNetwork | Kind::ViscoelasticNetwork => {
                return Err(format!(
                    "adaptive: {} is a {}, whose bodies have no error account here; \
                     run it with a fixed step_s",
                    element.name,
                    element.kind.name()
                ));
            }
        };
        for (field, coordinate) in coordinates {
            expected.push((format!("{}.{field}", element.name), element.node, coordinate));
        }
    }
    if expected.is_empty() {
        return Err("adaptive stepping needs at least one dynamic state or source-work account".into());
    }
    for (name, _) in tolerances {
        if !expected.iter().any(|(known, _, _)| known == name) {
            return Err(format!("adaptive.tolerances.{name}: not a dynamic state coordinate or source-work account"));
        }
    }
    let mut accounts = Vec::with_capacity(expected.len());
    for (name, node, coordinate) in expected {
        let what = format!("adaptive.tolerances.{name}");
        let value = tolerances.iter().find(|(key, _)| key == &name).map(|(_, value)| value)
            .ok_or_else(|| format!("{what}: missing explicit tolerance"))?;
        fields(value, &["absolute", "relative"], &what)?;
        let absolute = positive(value, "absolute", &what)?;
        let relative = number(value, "relative", &what)?;
        if !(relative.is_finite() && relative >= 0.0) {
            return Err(format!("{what}.relative must be finite and nonnegative"));
        }
        accounts.push(ErrorAccount { node, coordinate, tolerance: ErrorTolerance { absolute, relative } });
    }
    Ok(AdaptiveSettings {
        config: AdaptiveConfig { order: order as u32, min_step_s: minimum, max_attempts: attempts as u32 },
        accounts,
        attempts: 0,
        trial_steps: 0,
        last_error_norm: f64::NAN,
        max_error_norm: 0.0,
        suggested_dt_s: step_s,
    })
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
            "adaptive",
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
        if let Some(list) = connection
            .get("join")
            .or_else(|| connection.get("common_effort"))
        {
            let common_effort = connection.get("common_effort").is_some();
            let operation = if common_effort { "common_effort" } else { "join" };
            fields(connection, &[operation], &what)?;
            let ports = list
                .as_array()
                .ok_or_else(|| format!("{what}: {operation} is an array of ports"))?
                .iter()
                .map(|entry| {
                    let (text, ratio) = member(entry, &what)?;
                    Ok((resolve(&builder, &elements, &text)?, ratio))
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            join(&mut builder, &elements, &ports, &mut guesses, common_effort, &what)?;
        } else if connection.get("transformer").is_some() {
            transformer(&mut builder, &elements, connection, &what)?;
        } else if connection.get("route").is_some() {
            route(&mut builder, &elements, connection, &mut guesses, &what)?;
        } else if connection.get("signal").is_some() {
            signal(&mut builder, &elements, connection, &mut guesses, &what)?;
        } else if connection.get("prescribe").is_some() {
            prescribe(&mut builder, &elements, connection, &what)?;
        } else {
            return Err(format!(
                "{what}: expected join, common_effort, transformer, route, signal or prescribe"
            )
            .into());
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
    let adaptive = root.get("adaptive")
        .map(|definition| adaptive(definition, &elements, step_s)).transpose()?;
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
        adaptive,
    })
}

impl Device {
    fn no_claim(&self) -> &'static str {
        if self.adaptive.is_some() {
            ADAPTIVE_NO_CLAIM
        } else if self.elements.iter().any(|element| !element.bodies.is_empty()) {
            NETWORK_NO_CLAIM
        } else {
            NO_CLAIM
        }
    }

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
                let time_s = if let Some(adaptive) = &mut self.adaptive {
                    let proposed = adaptive.suggested_dt_s.min(self.step_s).min(remaining);
                    let mut config = adaptive.config;
                    // A caller-imposed final horizon is not an accuracy retry.
                    // Its own retries must still respect this smaller floor.
                    config.min_step_s = config.min_step_s.min(remaining);
                    let report = self.system.step_adaptive(
                        proposed, config,
                        |initial, coarse, refined| adaptive.metric(initial, coarse, refined),
                        &cx,
                    ).map_err(RunError::Adaptive)?;
                    self.iterations += report.iterations;
                    self.shortenings += report.shortenings;
                    adaptive.record(&report);
                    report.time_s
                } else {
                    let report = self.system.step(self.step_s.min(remaining), &cx)
                        .map_err(RunError::Step)?;
                    self.iterations += u64::from(report.iterations);
                    self.shortenings += u64::from(report.shortenings);
                    report.time_s
                };
                if time_s >= next_sample - slack {
                    self.sample();
                    let every = self.sample_every_s.unwrap_or(f64::INFINITY);
                    while next_sample <= time_s + slack {
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
    fn values(&self, element: &Element) -> Vec<(String, f64)> {
        if matches!(element.kind, Kind::RigidNetwork | Kind::ViscoelasticNetwork) {
            return self.network_values(element);
        }
        self.lumped_values(element)
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect()
    }

    /// Every body of a network (position and velocity of its centre of mass
    /// in world axes, orientation, angular velocity in its own axes), the
    /// mean heat rate of every connector over the last step, and the energy
    /// stored in every Maxwell branch.
    fn network_values(&self, element: &Element) -> Vec<(String, f64)> {
        let system = &self.system;
        let (states, memory) = if element.kind == Kind::ViscoelasticNetwork {
            let state = system.state::<RelaxationState>(element.node);
            (
                state.map(|state| state.bodies.as_slice()),
                state.map(|state| state.memory.as_slice()),
            )
        } else {
            let bank = system.state::<Vec<RigidBodyState>>(element.node);
            (bank.map(Vec::as_slice), None)
        };
        let unknown = Vec3::new(f64::NAN, f64::NAN, f64::NAN);
        let mut values = Vec::with_capacity(13 * element.bodies.len() + element.branches.len());
        for (index, properties) in element.bodies.iter().enumerate() {
            let state = states.and_then(|states| states.get(index));
            let position = state.map_or(unknown, |state| state.pose().position_world());
            let velocity = state
                .and_then(|state| state.center_of_mass_velocity_world(*properties).ok())
                .unwrap_or(unknown);
            let spin = state
                .and_then(|state| {
                    properties
                        .angular_velocity_body_checked(state.angular_momentum_body())
                        .ok()
                })
                .unwrap_or(unknown);
            for (field, unit, vector) in [
                ("position", "m", position),
                ("velocity", "m_s", velocity),
                ("angular_velocity", "rad_s", spin),
            ] {
                for (axis, component) in [("x", vector.x), ("y", vector.y), ("z", vector.z)] {
                    values.push((format!("body{index}.{field}_{axis}_{unit}"), component));
                }
            }
            let orientation = state.map_or([f64::NAN; 4], |state| {
                state.pose().orientation().components()
            });
            for (axis, component) in ["w", "x", "y", "z"].into_iter().zip(orientation) {
                values.push((format!("body{index}.orientation_{axis}"), component));
            }
        }
        for (connector, &branches) in element.branches.iter().enumerate() {
            let heat_w = match system.port(element.node, &format!("heat{connector}")) {
                Some(port) if system.steps() > 0 => -system.last_output(port)[0],
                _ => f64::NAN,
            };
            values.push((format!("connector{connector}.heat_w"), heat_w));
            for branch in 0..branches {
                let energy_j = memory
                    .and_then(|memory| memory.get(connector)?.get(branch))
                    .map_or(f64::NAN, |amplitude| {
                        0.5 * amplitude.iter().map(|a| a * a).sum::<f64>()
                    });
                values.push((format!("connector{connector}.branch{branch}.energy_j"), energy_j));
            }
        }
        values
    }

    fn lumped_values(&self, element: &Element) -> Vec<(&'static str, f64)> {
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
            Kind::NodeCapacitor => {
                let charge = system.state::<f64>(element.node).copied().unwrap_or(f64::NAN);
                vec![
                    ("charge_c", charge),
                    ("voltage_v", charge / element.capacitance_f.unwrap_or(f64::NAN)),
                ]
            }
            Kind::NonlinearCapacitor => {
                let charge = system.state::<f64>(element.node).copied().unwrap_or(f64::NAN);
                let voltage = element.charge_voltage_curve.as_ref()
                    .and_then(|curve| curve.voltage_v(charge).ok()).unwrap_or(f64::NAN);
                vec![("charge_c", charge), ("voltage_v", voltage)]
            }
            Kind::NonlinearInductor => {
                let flux = system.state::<f64>(element.node).copied().unwrap_or(f64::NAN);
                let current = element.flux_current_curve.as_ref()
                    .and_then(|curve| curve.current_a(flux).ok()).unwrap_or(f64::NAN);
                vec![("flux_wb", flux), ("current_a", current)]
            }
            Kind::HeatCapacity => vec![(
                "temperature_k",
                system.state::<f64>(element.node).copied().unwrap_or(f64::NAN),
            )],
            Kind::Source(_) | Kind::CurrentSource | Kind::Heater => vec![(
                "energy_in_j",
                system
                    .state::<SourceAccount>(element.node)
                    .map_or(f64::NAN, |account| account.energy_in_j),
            )],
            Kind::Resistance(domain) => {
                vec![(domain.names().effort, output("a")), ("heat_w", -output("heat"))]
            }
            Kind::ResistorBranch => {
                vec![("current_a", output("a")), ("heat_w", -output("heat"))]
            }
            Kind::Gyrator => vec![("effort_a", output("a")), ("effort_b", output("b"))],
            Kind::Conductance => vec![("heat_rate_w", output("a"))],
            // Read by `network_values`.
            Kind::RigidNetwork | Kind::ViscoelasticNetwork => Vec::new(),
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
        if let Some(adaptive) = &self.adaptive {
            let _ = writeln!(
                out,
                "adaptive order {}: {} attempts, {} trial steps; last normalized error {:.6e}, maximum accepted {:.6e}; suggested next interval {:.6e} s\n",
                adaptive.config.order, adaptive.attempts, adaptive.trial_steps,
                adaptive.last_error_norm, adaptive.max_error_norm, adaptive.suggested_dt_s
            );
        }
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
        let _ = writeln!(out, "\nNo claim: {}.", self.no_claim());
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
                object(self.values(element).into_iter()),
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
        out.push(']');
        if let Some(adaptive) = &self.adaptive {
            let _ = write!(
                out,
                ",\"adaptive\":{{\"order\":{},\"attempts\":{},\"trial_steps\":{},\"last_error_norm\":{},\"max_error_norm\":{},\"suggested_dt_s\":{}}}",
                adaptive.config.order, adaptive.attempts, adaptive.trial_steps,
                num(adaptive.last_error_norm), num(adaptive.max_error_norm), num(adaptive.suggested_dt_s)
            );
        }
        let _ = writeln!(out, ",\"no_claim\":{}}}", quote(self.no_claim()));
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
