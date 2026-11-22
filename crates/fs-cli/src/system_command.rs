//! `frankensim system`: a device described as data, run on the composition
//! kernel.
//!
//! The scene names elements from [`fs_system::elements`] and the connections
//! between their ports; this module lowers it to a [`SystemBuilder`], steps
//! the sealed [`System`] for the requested duration, and reports the final
//! states, sampled states and the energy, entropy and stream-inventory ledger.
//! Nothing here integrates anything: every number comes from the kernel and its elements.

mod adaptive_metric;
mod conduction;
mod contact;
mod dirac;
mod friction;
mod gas;
mod hydraulic;
mod materials;
mod modal;
mod phs;
mod snapshot;
mod solid;
mod timeline;

use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::Read as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use std::rc::Rc;

use fs_conduction::transient::enthalpy::material::MaterialEnthalpyTransport;
use fs_control::component::{SampledSpeedController, SampledSpeedState};
use fs_control::speed::{ShaftPowerLimits, SpeedControllerParameters};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_gas::component::ClosedGasParcel;
use fs_material::conductor::component::{OhmicResistor, Terminals};
use fs_material::conductor::graph::MaterialOhmicNode;
use fs_material::fluid_transport::laminar_pipe::{LaminarPipeComponent, MaterialLaminarPipe};
use fs_material::phase::component::UniformEnthalpyBody;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
};
use fs_material::thermal::conduction::component::LayerConductor;
use fs_material::thermal::conduction::{
    ConductivityCurve, LayerConductionConfig, LayerGeometry, MAX_CONDUCTIVITY_KNOTS,
};
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
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
};
use fs_mbd::spatial_loads::internal_motion::viscous::network::ViscousNetworkMotionControls;
use fs_mbd::spatial_loads::internal_motion::viscous::network::relaxation::RelaxationNetworkMotionControls;
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_ports::PortKind;
use fs_ports::component::{
    Carrier, Component, PortRole, PortSpec, PortTangent, STREAM_CARRIED, StreamQuantity,
};
use fs_qty::chemistry::SpeciesId;
use fs_system::constitutive::{
    CapacitorChargeCurve, InductorFluxCurve, MAX_CAPACITOR_KNOTS, MAX_INDUCTOR_KNOTS,
    NonlinearCapacitor, NonlinearInductor,
};
use fs_system::elements::{
    Admittance, AffineSignal, Conductance, ControlledSource, ControllerState, Damper, FlowSource,
    Force, GasContents, GasVolume, Gyrator, HeatCapacity, HeatNode, Heater, IdealDiode, IdealSwitch, Mass,
    ModulatedTransformer, Motion, Resistance, Restriction, SampledController, SignalClock,
    SignalSource, Slider, SliderEnd, Slip, SourceAccount, Spring, Storage, Waveform,
};
use fs_system::ledger::{Account, CarriedAccount, Ledger};
use fs_system::{
    AdaptiveConfig, AdaptiveError, AdaptiveReport, Checkpoint, ErrorNorm, ErrorTolerance, NodeId,
    PortRef, RouteTemperature, SolveConfig, System, SystemBuilder, SystemError, heat_sink,
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
/// Quantities one stream of a scene may carry besides mass and enthalpy.
const MAX_CARRIED: usize = 16;
/// Bodies one `heat_node` may join.
const MAX_NODE_BODIES: usize = 64;
/// Coordinates in one dense, positive-definite resistance law.
const MAX_RESISTANCE_COORDINATES: usize = 64;
/// Each input, output or memory dimension of an authored signal law.
const MAX_SIGNAL_COORDINATES: usize = 64;
/// Authored caloric support retained by one uniform enthalpy body.
const MAX_ENTHALPY_KNOTS: usize = 4096;
const ADAPTIVE_NO_CLAIM: &str = "step-doubling estimates control the declared state coordinates and source-work increments under the selected order assumption; they are local numerical estimates, not global error bounds or physical validation; dissipated heat and interface balances retain the kernel's separate ledger audits; no material-card resolution, contact, friction or phase change";

const HELP: &str = "Usage: frankensim [--json] system <scene.json>
       [--until <time_s> [--whole-steps]] [--checkpoint <new-file>]
       [--resume <file> --resume-root <hex>]

Run a lumped mechanical, rotational, electrical, hydraulic, gas and thermal
device described as data on the composition kernel, and report its final state,
sampled states and the energy, entropy and measured stream-inventory ledger of
the run.

--until stops at an absolute simulated time within the scene's duration; zero
can save the admitted initial state. --checkpoint saves the accepted kernel
state, sampled history and adaptive driver state to a new file, and reports
its resume_root. Existing files are never overwritten. To continue in another
process, supply --resume together with that independently retained 64-digit
hex root. The scene contents and executable must match the saved run. A scene
whose components lack native state codecs refuses checkpoint creation.

--until cuts the step it falls in, and that cut is a step boundary a run with
no pause does not have. With --whole-steps, which needs --until, no step is
cut: the run stops at the end of the first step that ends at or after that
time, and reports the time it stopped at. A run saved there and resumed takes
the steps a run with no pause takes.

For example:
  frankensim --json system scene.json --until 0.5 --checkpoint prefix.fssnap
  frankensim --json system scene.json --resume prefix.fssnap --resume-root <root>

Scene schema: frankensim.system.v0. Top-level fields:
  duration_s        simulated time
  step_s            largest step; the kernel may take shorter ones
  sample_every_s    optional; record every element's state at the end of the
                    first step at or after each multiple of this interval
  solver            optional {relative_tolerance, absolute_tolerance, max_iterations}
  adaptive          optional {order = 1, min_step_s = 1e-12, max_attempts = 32,
                    tolerances: {element.state_field: {absolute, relative}, ...}}
  materials         optional object: name -> thermal material-card declaration
  elements          object: name -> element
  connections       array of connections

Adaptive stepping requires an explicit tolerance for every dynamic state
coordinate and source energy_in_j account. Absolute tolerances are positive
in that field's SI units; relative tolerances are nonnegative. Source-work
tolerances apply to each trial interval's energy increment. Derived endpoint
voltage/current and stateless last-step means are not independent accounts.
A gas_volume needs mass_kg, energy_j, volume_m3 and carried<i>_kg for each
carried quantity; a closed gas_cavity needs volume_m3 and temperature_k.
A streaming gas_cavity additionally needs mass_kg, energy_j and
species<i>_mass_kg for every species in its declared stream catalog. Its
material-transfer step is first order, so adaptive.order must be 1.
Every network body needs body<i>.position_m, body<i>.velocity_m_s,
body<i>.angular_velocity_rad_s and body<i>.orientation_rad under its element
name. These control complete vector differences and rotation angle, not
individual quaternion coefficients. Each viscoelastic branch also needs
connector<i>.branch<j>.amplitude_sqrt_j for its complete signed memory vector.
Vector relative scales use initial/refined lengths; orientation uses angular
motion during the interval. All body and branch accounts are mandatory.
The default order is 1; choose order 2 only for a smooth second-order model.
Affine maps, signal sources, modulated restrictions, published inertance
observations and controlled commands supplied by signals or point-sampled
value admit order 1 only. SI-named controlled-source prescriptions retain
their interval-mean rules. sampled_controller and
speed_controller require fixed-step operation; their sample events still cap
steps exactly, but no adaptive discrete-memory metric is supplied. ideal_diode
and ideal_switch also require fixed stepping; no switching error metric is supplied.
The refined two-half-step state is retained; rejected trials leave no state
or ledger entries. step_s remains the largest interval.
The final externally capped interval may be shorter than min_step_s; an
accuracy rejection still cannot reduce that final interval further.

Thermal material declarations prepare original laws from normalized card packs:
  pack               path to a .fsmcdpk file, relative to the scene
  phase              solid | liquid (explicit pure phase)
  temperature_axis   the source's temperature coordinate name
  temperature_kind   legacy | dimensional | absolute-temperature
  state              fixed source context: [{name, value, kind}, ...], or [];
                     value is a unit-bearing string; omit the temperature axis
  temperature_min_k, temperature_max_k  admitted source support
  reference_specific_enthalpy_j_kg     enthalpy at temperature_min_k
  heat_capacity_property, density_property  original dimensional property names
  conductivity_property  optional original k(T) property name
  conductivity_kind  dimensional (default) | thermal-conductivity
  resistivity_property  optional original electrical resistivity property name
                     (ohm m), for card-backed resistors
  selection          single-claim-only (default) | prefer-observation-backed
The temperature convention and fixed context must match the source. Preparation
retains nominal Cp(T), density(T), and optional k(T) and resistivity(T); no
initial-value freezing or extrapolation. Reports retain content identities, claims and source receipts;
source uncertainty is retained but not propagated. At most 64 declarations,
4 MiB per pack and 64 MiB total pack bytes. No phase transition is inferred.

Elements (ports in brackets; physical fields use SI units):
  mass           mass_kg, position_m = 0, velocity_m_s = 0            [p]
                 publish = false; true adds scalar position and velocity
                 signals from accepted state at the start of the step.
                 publish_state is an equivalent alias; if both are supplied
                 they must agree. Both default to false.
                 The same option applies to inertia and inductor, with
                 angle/angular velocity or charge/current in their SI units.
  spring         stiffness_n_m, cubic_n_m3 = 0, extension_m = 0       [a, b]
  damper         damping_n_s_m                                        [a, b, heat]
  resistance_matrix domain, matrix (square, 1 to 64 rows)              [u, heat]
                 domain is translation, rotation, electrical or hydraulic;
                 matrix is symmetric positive definite, in SI effort/flow units
                 (N s/m, N m s/rad, ohm or Pa s/m^3). u receives a flow vector and
                 returns matrix times that vector; heat must be connected.
  force          force_n (waveform)                                   [p]
  affine_signal  matrix (m by n), bias (m numbers)      [input0..input<n-1>, output]
                 y = bias + matrix * u; n scalar inputs, one m-wide output.
  signal_source  waveforms (array of waveforms)                       [output]
                 One output coefficient per waveform, evaluated at the start
                 of each step. The source clock and definitions are replayable;
                 table knots cap steps. Signals carry no physical energy.
  sampled_controller sample_period_s, a, b, c, d, initial
                                                     [input0..input<m-1>, output]
                 Discrete law y = C*x + D*u, x_next = A*x + B*u at t = k*period,
                 starting at t = 0; output is held between samples. initial
                 gives n state coordinates; A is n by n, B n by m, C p by n,
                 D p by m. n may be zero (A and B [], C rows []); m,p >= 1.
                 Use fixed step_s; sample events cap larger physical steps.
  controlled_effort domain                                           [p, command]
  controlled_flow domain                                             [p, command]
                 domain is translation, rotation, electrical or hydraulic.
                 command is a scalar held effort or flow in that domain's SI unit.
                 Prescribe force_n, torque_n_m, voltage_v or pressure_pa for effort;
                 velocity_m_s, angular_velocity_rad_s, current_a or flow_m3_s for flow.
                 Those SI-named prescriptions supply the interval mean;
                 value instead supplies a signed start-of-step point sample.
                 Flow is positive INTO the source: negative current or liquid
                 flow supplies a node_capacitor or accumulator. Reports the declared
                 SI field, command and signed energy_in_j, negative when
                 energy is supplied.
  Signal matrices and vectors contain finite numbers, with at most 64 inputs,
  outputs or memory coordinates. Signal gains, units and initial memory are
  authored; no limiting, power-supply model or control stability is inferred.
  Reports expose output0, output1, ...; sampled controllers also expose x0,
  x1, ..., sample_count and time_s. signal_source exposes accepted time_s.
  speed_controller sample_period_s, proportional_gain_nm_s_per_rad,
                 integral_gain_nm_per_rad, anti_windup_time_s,
                 positive_torque_limit_nm, negative_torque_limit_nm,
                 motoring_power_w, regeneration_power_w,
                 integral_torque_nm = 0
                 [reference, feedforward, rotor_speed, stator_speed, torque]
                 fs-control's sampled PI policy for rotor-minus-stator speed;
                 torque drives a rotational controlled_effort.command. Its
                 four inputs are scalar signals: reference in relative rad/s,
                 feedforward in N m, both measured shaft speeds in rad/s.
                 Prescribe them with relative_speed_rad_s, torque_n_m and
                 angular_velocity_rad_s respectively, or connect signals.
                 Prescribed signals are sampled at the interval's start.
                 The kernel stops at each sample; torque is held across plant
                 substeps. Reports integral_torque_nm, held_torque_nm,
                 sample_phase_s and sample_count. Fixed step_s only: adaptive
                 refuses this discrete policy. No physical actuator, sensor
                 dynamics, current allocation or stability claim is inferred.
  velocity_source velocity_m_s (waveform)                             [p, span]
                 (prescribed motion: span publishes its start/end velocities;
                 join p to a slider end and signal span to that end's span input)
  slider         normal_n, static_coefficient, kinetic_coefficient,
                 slip_m_s = 0, prescribed_ends = []                  [a, b, heat]
                 (ideal dry friction; 0 < kinetic <= static; slip_m_s is
                 velocity at a minus velocity at b, matching the connected
                 bodies' initial velocities; zero starts stuck. Each named
                 prescribed end, \"a\" or \"b\", adds a signal input span_a or
                 span_b. Tell every accelerating prescribed end its span.
                 Reports slip_m_s, stuck (1 or 0), force_n on a, and heat_w;
                 requires fixed step_s and the kernel's load-line solve.)
  friction_interface normal_n, law, slip_m_s = 0, solve = {}         [a, b, heat]
                 (an fs-tribo dry-friction law between a and b. law is
                 \"rosin_woodhouse_galluzzo_2025_eq13\" (the measured rosin
                 fit fs-tribo cites) or one object key with its coefficients:
                   coulomb {static_mu, kinetic_mu}
                   velocity_dependent {static_mu, mu_zero, slope_per_speed_s_m}
                   stribeck {static_mu, kinetic_mu, characteristic_speed_m_s,
                     viscous_per_speed_s_m = 0}
                   hyperbolic {static_mu, kinetic_mu, characteristic_speed_m_s}
                   shifted_power {static_mu, kinetic_mu,
                     characteristic_speed_m_s, exponent}
                 slip_m_s as for slider. solve holds budgets, defaults:
                 max_iterations 64, force_absolute_tolerance_n 1e-12,
                 force_relative_tolerance 1e-12, maximum_force_n 1000 static
                 capacities (at least 1 N), impulse_tolerance_n_s 1e-10,
                 work_tolerance_j 1e-12. A moving surface is a prescribed
                 end or a velocity_source joined to one. A prescribed end
                 is held over each step: no span input, so an accelerating
                 one advances in jumps. reads_ends = [] names
                 the ends, \"a\" or \"b\", joined to a neighbour that tells
                 the end of its flow (a modal_structure attachment); the
                 kernel wires their end inputs a.end and b.end. Reports
                 slip_m_s, stuck, force_n on a and heat_w; requires fixed
                 step_s.)
  sliding_interface normal_n, law, slip_rate_m_s, solve = {} [heat, traction]
                 (the same laws driven at a prescribed slip rate
                 [x, y] m/s along the interface's two tangent axes, under
                 the constant normal_n; a traction signal output [x, y] N.
                 solve as for friction_interface, less
                 impulse_tolerance_n_s. Join heat to a heat capacity or heat
                 node, as a heater's. The prescribed motion's work is the
                 interface's own account, as a prescribed force's is.
                 Reports work_j, traction_x_n, traction_y_n and heat_w;
                 requires fixed step_s.)
  hertz_patch    normal_n, radius_m, static_mu, kinetic_mu,
                 motion_m_s = [0, 0, 0], solve = {}         [wrench, heat]
                 (fs-tribo's joint slip and spin friction of one circular
                 Hertz patch under a constant normal load and radius, on
                 the load line. wrench (width 3) receives the motion of
                 side a against side b, [v_x, v_y, radius_m w] m/s along
                 the patch's two tangent axes and about its normal, and
                 returns minus the wrench on side a, [F_x, F_y, M /
                 radius_m] N; connect it with a power_map of members from
                 the bodies' force and torque ports: rows t1, t2 and 0 for
                 a force port, r x t1, r x t2 and radius_m n for a torque
                 port with contact arm r, all negated for side a.
                 motion_m_s is the motion that map gives at the start; a
                 nonzero one starts sliding. solve holds the owner's
                 budgets, defaults: maximum_static_iterations 96,
                 maximum_iterations 256, maximum_evaluations 512,
                 force_tolerance_n 2e-9, moment_tolerance_nm 2e-10,
                 displacement_tolerance_m 2e-13, energy_tolerance_j 2e-11,
                 relative_matrix_tolerance 1e-13, maximum_force_n 1000
                 static capacities (at least 1 N), maximum_spin_moment_nm
                 that times radius_m (at least 1 N m). A step in which the
                 slip or the spin would reverse is refused. Join heat as a
                 heater's. Reports slip_x_m_s, slip_y_m_s, spin_speed_m_s,
                 owner_steps, force_x_n, force_y_n, spin_moment_nm and
                 heat_w; requires fixed step_s.)
  modal_structure sample_rate_hz, modes | network, attachments
                                                     [<attachments>, heat]
                 (fs-nlmodal's linear structure in mass-normalized damped
                 modes. modes: [{frequency_hz, damping_ratio,
                 displacement = 0, velocity = 0}], in modal coordinates.
                 attachments: {<id>: {domain, weights}}: a power port per
                 attachment, one row of weights per port row, one weight per
                 mode (mode shapes at the point, 1/sqrt(kg) for a force).
                 Or network: {stiffness, mass, rayleigh = [0, 0]}: a
                 lumped network's symmetric K and positive definite M over
                 its 1 to 256 degrees of freedom (N/m and kg for
                 translation), whose modes fs-nlmodal solves; rayleigh
                 [a, b] damps it by C = a M + b K, so mode j has the ratio
                 (a + b w_j^2) / (2 w_j). Each attachment row then weights
                 the degrees of freedom (one weight per degree of freedom),
                 and start = {displacement, velocity}, one entry per degree
                 of freedom (default at rest), is carried to the modes. A K
                 that does not hold every degree of freedom (a free body)
                 is refused.
                 sample_rate_hz is the owner's clock: modes must lie below
                 0.9 of its Nyquist frequency; steps of any length are exact
                 for the effort held over them. Each attachment also tells
                 the end of its flow, <id>.end. Reports every mode's
                 displacement and velocity, each attachment's flow at the
                 end of the last step, and heat_w; requires fixed step_s.)
  contact_interface stiffness, exponent, gap_m = 0, closure_m = 0, solve,
                 reads_ends = []                                     [a, b]
                 (fs-dcontact's one-point normal contact on the load line:
                 the reaction R = stiffness [x - gap_m]_+^exponent pushes a
                 and b apart, where the closure x is the displacement of a
                 towards b (x' = flow at a less flow at b) and starts at
                 closure_m. exponent >= 1 (1.5 is Hertz); stiffness in
                 N/m^exponent. Each step's reaction is the secant of the
                 law's potential over the step's move, so a lossless
                 contact gives back the work it takes. solve holds the
                 search budgets: maximum_force_n and maximum_penetration_m
                 (required ceilings), max_iterations 100 (1 to 128),
                 force_absolute_tolerance_n 1e-9, force_relative_tolerance
                 1e-12. reads_ends as for friction_interface. Reports
                 closure_m and reaction_n; requires fixed step_s.)
  heat_capacity  capacity_j_k, temperature_k                          [heat]
                 (and a signal output, temperature: its temperature at the
                 start of the step; it may stay unconnected)
  enthalpy_body  mass_kg, specific_enthalpy_j_kg,                     [heat]
                 phase = \"solid\" | \"liquid\" | \"solid-liquid\",
                 enthalpy_curve = [[h_j_kg, temperature_k,
                   liquid_mass_fraction, bulk_density_kg_m3], ...]
                 (2 to 4096 knots; h strictly increases. Solid-liquid charts
                 start fully solid and end fully liquid; temperature and
                 liquid fraction never decrease. Equal temperatures retain
                 latent heat. Single-phase charts have strictly rising T
                 and the declared phase at every knot. The original material
                 owner supplies sensible/latent storage and entropy-conjugate
                 heat-port temperature; endpoint probes report h, T, liquid
                 fraction and density. Initial h is required: T alone cannot
                 locate a state on a latent plateau. No extrapolation, phase
                 kinetics or mechanical expansion; authored charts confer no
                 material-pack provenance.)
                 Alternatively, material = \"name\" replaces phase and
                 enthalpy_curve with the declared source calorics.
  conductance    conductance_w_k                                      [a, b]
  layer_conductor conductivity_a, conductivity_b                     [a, b]
                 (each 2 to 4096 [temperature_k, k_w_m_k] knots),
                 geometry = {kind: planar, distances_m: [a, b], area_m2}
                   or {kind: cylindrical, radii_m: [a, b],
                       contact_radius_m, length_m}
                   or {kind: spherical, radii_m: [a, b], contact_radius_m},
                 contact_resistance_m2_k_w = 0,
                 local_solver = {max_iterations = 80,
                   absolute_heat_tolerance_w = 1e-9,
                   relative_heat_tolerance = 1e-12}
                 (positive piecewise-linear k(T), increasing temperature
                 knots; the original layer owner integrates the charts and
                 solves the contact jump. Radii must straddle the contact.
                 Stateless source-free layers at fixed geometry; no internal
                 layer storage or source extrapolation; authored charts confer
                 no material provenance.
                 heat_rate_w is positive from a to b.)
                 Alternatively, materials = [\"a\", \"b\"] replaces both
                 conductivity charts; both named sources must prepare k(T).
  heat_node      conductances_w_k (one per body), temperature_k (optional)
                                                [source, body0, body1, ...]
                 (a plane of no heat capacity between bodies: join each
                 body<k> to a heat capacity and a heater to source, or
                 prescribe heat_rate_w on source. The heat divides by where
                 the bodies' temperatures put the node; no fraction is
                 given. temperature_k is not a state: with a heater on
                 source it is where the solve for the node starts.)
  thermal_timeline laws, cells, links = [], events = [], ports = [],
                 environment = false, publish = [], max_step_s,
                 energy_tolerance_j = 1e-9
                         [cell<id> per port, environment, temperature<id>]
                 (enthalpy cells stepped by the fs-thermal timeline owner's
                 explicit rule. laws names 1 to 64 charts, each {phase,
                 enthalpy_curve} as for enthalpy_body or {material: \"name\"}.
                 A cell is {id, law, mass_kg, specific_enthalpy_j_kg,
                 volume_m3, surface_area_m2, start_m, end_m, boundary}; a
                 boundary is {ambient_k, htc_w_m2_k = 0, source_w = 0,
                 reservoir = {area_m2, conductance_w_k, temperature_k}},
                 and a cell without one is insulated. A link is {a, b,
                 conductance_w_k, area_m2, id = its index}, made at the
                 start; area_m2 > 0 is the surface it covers on each cell,
                 which no film then reaches. Events are {time_s, contact: link with its id},
                 {time_s, separation: link id} or {time_s, boundary: {cell,
                 ...}}, in time order; every cell is born at the start.
                 ports lists the cells with a heat port cell<id>: it takes
                 the heat rate into the cell and returns the cell's mean
                 temperature over the step. environment, required when a
                 boundary exchanges heat, returns the mean rate heat entered
                 through the boundaries; prescribe its temperature, which
                 the owner does not read. publish lists the cells with a
                 signal output temperature<id>: the temperature at the start
                 of the step. Steps are cut to max_step_s, to the owner's
                 stability bound and to event times. Reports each cell's
                 temperature_k and specific_enthalpy_j_kg, and with an
                 environment port environment_heat_w. Fixed step_s only. Its
                 state is saved by --checkpoint and resumed.)
  conduction_body box_m [Lx, Ly, Lz], cells [nx, ny, nz], conductivity_w_m_k
                 (or material, newton), volumetric_heat_capacity_j_m3_k
                 (or reference_temperature_k), temperature_k, ports,
                 traces, solve                          [one per port face]
                 (fs-conduction's backward-Euler finite-element body on the
                 box [0, Lx] x [0, Ly] x [0, Lz], nx x ny x nz cells of six
                 tetrahedra, at most 64 a side and 4096 in all, starting at
                 a uniform temperature_k. ports names box faces, x_min,
                 x_max, y_min, y_max, z_min, z_max: each takes the heat rate
                 into the body and returns its face's area-mean end
                 temperature. traces names faces whose port carries one
                 value for every vertex of the face: a heat load in, the
                 end temperature out; a trace_contact joins two of them. A
                 face is under ports or under traces; the faces not named
                 are insulated. solve gives energy_tolerance_j (the owner's
                 energy gate on each step, required) and the Krylov
                 tolerance 1e-12, max_iterations 20000 and restart 60.
                 In place of conductivity_w_m_k, material names a scene
                 material with a conductivity_property: the body takes the
                 k(T) curve its card gives over the material's declared
                 span, which must hold temperature_k; its ports and traces
                 are any body's.
                 A k that varies steps by the owner's Newton solve; newton
                 gives max_iterations 32, residual_rtol 1e-10 and
                 residual_atol_j 1e-10. A step that would leave the span is
                 refused. The heat capacity stays the declared constant,
                 unless the body omits volumetric_heat_capacity_j_m3_k: then
                 its heat capacity and density are the material card's too.
                 fs-conduction's enthalpy owner stores specific enthalpy on
                 the chart of the card's Cp(T) over the span, with each
                 node's mass the card's density at reference_temperature_k
                 (required) times its volume, fixed; every step is the
                 owner's Newton solve under newton, the tangent's solves
                 take solve's Krylov policy, and a trace_contact refuses
                 the body. Reports each port's <face>.temperature_k, each
                 trace's <face>.min_temperature_k and
                 <face>.max_temperature_k, the field's min_temperature_k and
                 max_temperature_k, and for a body of the card's heat
                 capacity reference_density_kg_m3. Fixed step_s only. Its
                 field is saved by --checkpoint and resumed.)
  trace_contact  between [\"body.face\", \"body.face\"], offset_m [x, y, z],
                 resistance_m2_k_w, source (or pack, state, selection,
                 follows), mapped                                     [a, b]
                 (fs-conduction's contact between two traces of two
                 conduction bodies, carried by the kernel's conductance
                 over interface points: heat crosses at each point of the
                 faces, and the contact generates entropy there. offset_m
                 is the second body's origin in the first's frame; the
                 owner refuses traces that do not lie in one plane.
                 resistance_m2_k_w is one declared area-specific resistance
                 for the whole patch and source says where it comes from;
                 both are reported. In their place pack names a compiled
                 interface pack beside the scene and state the point its
                 card is read at, once, as for a solid; selection is
                 single-claim-only unless given. The report then carries
                 the pack, card, claim and receipt identities. follows
                 names the card coordinate that takes the mean of the two
                 traces' temperatures, or two coordinates, the first
                 trace's and the second's: the card is then asked at every
                 interface point of every step with state as its reference
                 point, and a state it has no answer for refuses the step.
                 A contact has one of the two sources. Without mapped the
                 two faces must have the same points and triangles. mapped
                 gives the budgets of the owner's common refinement for
                 faces that do not match: plane_tolerance_m,
                 coverage_relative_tolerance, max_pair_tests,
                 max_overlap_triangles, all four. Join a to the first trace
                 and b to the second; any other join of them is refused.
                 Reports heat_rate_w from a to b, conductance_w_k,
                 resistance_m2_k_w and interface_points, and for a contact
                 that follows its card min_resistance_m2_k_w and
                 max_resistance_m2_k_w over the interface points at the
                 last step's end.)
  heater         power_w (waveform, never negative)                   [heat]
Gas in volumes joined by restrictions (a perfect gas with constant heat
capacities; all joined volumes must hold the same gas):
  gas_volume     gas_constant_j_kg_k, heat_capacity_ratio, volume_m3,
                 pressure_pa, temperature_k, carried = [] (the mass fraction
                 of each quantity carried with the gas), wall = false
                                                      [stream, heat, p if wall]
  restriction    conductance_kg_s_pa, carried = 0 (how many quantities)
                                                                      [a, b]
  modulated_restriction
                 the same fields, plus a dimensionless opening signal
                                                             [a, b, opening]
  A modulated restriction multiplies conductance by opening in [0, 1]:
  zero shuts the valve; one restores its full conductance. Connect a signal
  to opening, or prescribe opening as a waveform held at the step's start.
  An out-of-range opening is refused by the native proposal at runtime.
  Source-backed restrictions of either kind use stream_basis instead of
  carried (below). Modulated restrictions admit adaptive order 1 only.
  A restriction passes conductance (p_a - p_b) [kg/s] with the enthalpy and
  carried quantities of the side the gas leaves. Join each end to a volume's
  stream port; several restrictions may meet at one volume. heat takes the
  heat rate into the gas and returns its temperature, like a heat capacity.
  With wall = true, p is the pressure on a moving wall: join it to a mass
  as {\"port\": volume.p, \"ratio\": area_m2}, positive for a wall whose
  motion in the mass's positive direction enlarges the volume.
A gas of cited thermochemistry (fs-gas: frozen chemistry, ideal
gas, the calorics of its NASA-9 packs):
  gas_cavity     gas (a record \"thermochemical T_K p_Pa mole|mass N
                 species fraction ...\"), sources = [{species, nasa9,
                 transport}, ...] (canonical pack paths relative to the
                 scene, one triple per species), volume_m3,
                 temperature_min_k, temperature_max_k (the bracket of the
                 owner's caloric solve), energy_tolerance_j_kg = 1e-7,
                 max_iterations = 128, insulated = false,
                 stream_mass_tolerance_kg (optional; enables material flow)
                                        [p, heat unless insulated, stream if enabled]
  p receives the compression flow [m3/s], positive when the cavity
  shrinks, and returns the step's work pressure [Pa]: prescribe flow_m3_s
  on it, or join it to a mass as {\"port\": cavity.p, \"ratio\": area_m2},
  positive for a piston whose motion in the mass's positive direction
  enlarges the cavity. heat is a heat capacity's. Each step
  resolves the packs at its end state, transport claims included, so the
  gas stays inside every claim's validity (the seed-v1 transport cards
  hold from 75 to 101.325 kPa).
  stream_mass_tolerance_kg must be positive and finite; it is the native
  inventory-transfer tolerance, not a relaxation of the kernel's mass audit.
  It cannot be combined with insulated: prescribe zero heat_rate_w on heat
  for an adiabatic streaming cavity. Its stream carries mass, enthalpy and
  each catalog species' mass. A restriction between these cavities names
  stream_basis: \"cavity.stream\" instead of carried, copying that port's
  source identity and coordinate order. Both ends still need ordinary joins;
  mismatched source catalogs or energy references refuse at connection.
  The reference works in any element order. Streaming state reports include
  energy_j in the native species datum, species<i>_mass_kg and
  mass_fraction<i> in the declared catalog order;
  stream_species in JSON names the species, including depleted zero slots.
The same four laws about a shaft, with angle, angular velocity and torque:
  inertia        inertia_kg_m2, angle_rad = 0, angular_velocity_rad_s = 0  [p]
  torsion_spring stiffness_n_m_rad, cubic_n_m_rad3 = 0, twist_rad = 0  [a, b]
  rotary_damper  damping_n_m_s_rad                                    [a, b, heat]
  torque         torque_n_m (waveform)                                [p]
and in a current loop, with charge, current and voltage (an inductor is the
loop's node; a capacitor and a resistor sit in series in it):
  inductor       inductance_h, charge_c = 0, current_a = 0            [p]
  capacitor      capacitance_f, charge_c = 0                          [a, b]
  resistor       resistance_ohm, or material, length_m, area_m2       [a, b, heat]
  voltage        voltage_v (waveform)                                 [p]
For a parallel circuit, the capacitor sets the node voltage and the branches
share that voltage (current is positive into each component):
  node_capacitor capacitance_f, charge_c = 0                          [p]
  resistor_branch resistance_ohm, or material, length_m, area_m2      [a, b, heat]
  current_source current_a (waveform; negative supplies the node)     [p]
  ideal_diode    no parameters; i >= 0, v <= 0, i*v = 0              [p]
  ideal_switch   no parameters; closed signal is exactly 0 or 1      [p, closed]
An ideal device's p receives current and returns voltage, positive power into
the device. Join a resistor_branch.a to a voltage.p, node_capacitor.p and
ideal device.p for a series RC circuit; prescribe 0 V at the resistor's b.
A closed switch has zero voltage; an open switch has zero current. Connect
a signal_source to closed or prescribe value there, using held_table for
binary schedules. Reports expose current_a, voltage_v and the switch's closed
command from the last step (null initially). Fixed stepping only; these ideal
relations contain no semiconductor loss, recovery, arcing or impulse model.
A resistor or resistor_branch that names a scene material whose declaration
names a resistivity_property is that material's uniform conductor of length
length_m and cross-section area_m2: in every trial fs-material's Ohmic node
reads the card's rho(T) length_m / area_m2 at the temperature heat receives,
and all of the Joule heat leaves through heat. A temperature outside the
material's declared span refuses the step. Its report adds temperature_k, the
temperature heat received in the last step, and resistance_ohm there. Fixed
geometry: no thermal expansion, skin effect or contact resistance.
The same parallel forms carry a liquid, with volume, flow and pressure:
  accumulator    capacitance_m3_pa, volume_m3 = 0                     [p]
                 (volume_m3 is the volume taken in; pressure_pa = volume_m3 /
                 capacitance_m3_pa)
  hydraulic_resistor resistance_pa_s_m3                              [a, b, heat]
  flow_source    flow_m3_s (waveform; negative supplies the node)     [p]
  pressure       pressure_pa (waveform)                               [p]
  carreau_pipe   viscosity_zero_pa_s, viscosity_infinite_pa_s,
                 time_constant_s, power_index (0 < n <= 1), transition (a > 0),
                 rate_support_per_s = [lowest, highest] (lowest 0: the core of
                 the pipe is at rest), shift (as a network's, below; required),
                 shift_mode = \"vertical\" | \"vertical-and-rate\",
                 radius_m, length_m, maximum_wall_rate_per_s       [a, b, heat]
                 (fully developed laminar flow of a liquid whose viscosity at
                 the shift's reference temperature is mu_inf + (mu_0 - mu_inf)
                 (1 + (time_constant_s rate)^a)^((n - 1)/a), solved over the
                 radius by the fs-material owner. heat takes the liquid's one
                 temperature, which multiplies the viscosity by the shift factor
                 (and with vertical-and-rate the rate too), and returns all the
                 viscous heat. Reports flow_m3_s from a to b, heat_w and the
                 temperature_k read in the last step. A wall rate above
                 maximum_wall_rate_per_s or outside the support, or a
                 temperature outside valid_k, refuses the step. No inertia,
                 compressibility, entrance loss or turbulence.)
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
  gyrator        ratio, a, b                                        [a, b]
  modulated_transformer a, b                                        [a, b, ratio]
                 a and b name translation, rotation, electrical or hydraulic.
                 a receives flow and b receives effort; the scalar ratio
                 signal sets e_a = ratio * e_b, f_b = -ratio * f_a.
                 Prescribe ratio as a waveform or connect a scalar signal.
                 Signed and zero ratios are valid. Reports the last-step
                 ratio, effort_a and flow_b in the declared domains' SI units.
                 This is an ideal held-ratio power map, not crank/cam geometry
                 or a displacement constraint; no ratio-actuator work is modeled.
An authored lossless N-port relation in one scalar power domain:
  dirac_junction domain, ports = [name, ...], incidence (n x n) [each port]
                 domain is translation, rotation, electrical or hydraulic;
                 name 2 to 64 distinct ports, in matrix row/column order.
                 Every port receives effort and returns flow = incidence * effort.
                 incidence must be finite and exactly skew, with zero diagonal;
                 singular and zero matrices are valid. Entries carry the declared
                 domain's SI flow/effort units. Reports each port's last-step
                 effort and flow as <port>.<SI quantity>, null before a step.
                 The junction has no state, heat outlet or adaptive error account.
A linear port-Hamiltonian system written as matrices, stepped by the exact
held-input owner of fs-system:
  port_hamiltonian energy, structure, dissipation (each n x n, 1 <= n <= 64),
                 input (n x m), ports = [{name, domain, role}, ...] (m of
                 them), state = [0, ...] (n numbers)   [each port, heat]
                 (dx/dt = (structure - dissipation) energy x + input u, with
                 stored energy x' energy x / 2 in J. energy and dissipation
                 are symmetric positive semidefinite, structure is skew. The
                 domain is translation, rotation, electrical or hydraulic.
                 An effort_in port receives an effort and returns the step mean of its row
                 of input' energy x as a flow; a flow_in port receives a flow
                 and returns an effort. heat is routed like a damper's and
                 carries the work dissipation absorbs. Reports x0, x1, ...,
                 each port's returned mean and heat_w; adaptive tolerances
                 name <element>.x<i>. Exact for held inputs only. Its state
                 is saved by --checkpoint and resumed.)
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
  A body may give solid = {pack (a material card pack, relative to the
  scene), state = [{name, value, kind}, ...] (the card's coordinates),
  selection = \"single-claim-only\", shape} in place of mass_kg and
  principal_inertia_kg_m2. The shape is centred on the body's origin along
  its axes: \"sphere\" (radius_m), \"cuboid\" (sides_m, three), \"ellipsoid\"
  (semiaxes_m, three) or \"cylinder\" (radius_m, length_m along z,
  inner_radius_m = 0 for a tube). Its mass and principal inertia are the
  shape's volume moments times the density the card gives at that state.
  With temperature_k, the body is freely expanded about its origin from the
  state's temperature to that one by the card's linear expansion
  coefficient; its mass is kept. The state then needs a T or exactly one
  absolute-temperature coordinate, and everything else in it stays fixed.
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
offset = 0}}, {\"table\": [[time_s, value], ...]} (linear, ends held), or
{\"held_table\": [[time_s, value], ...]} (constant between knots; a knot's
new value applies at that time). Both tables hold their first/last values
outside their span and require nonempty, finite, strictly increasing times.
Times and values are bounded to magnitude 1e150. Held tables schedule each
knot after the first and integrate each constant segment exactly.

Connections, ports written element.port:
  {\"join\": [ports]}       the ports meet at one node: a mass (inertia,
                          inductor) and the ends that share its velocity
                          (angular velocity, current), or a resistor_branch
                          end setting series current, or a heat capacity and
                          the elements at its temperature
  {\"common_effort\": [ports]}
                          one effort-setting port (node_capacitor,
                          nonlinear_capacitor, voltage, accumulator or
                          pressure) and branches sharing its voltage or
                          pressure (inductor, resistor_branch,
                          current_source, hydraulic_resistor, flow_source,
                          carreau_pipe)
  {\"transformer\": [effort_in_port, flow_in_port], \"ratio\": r}
                          lossless signed flow ratio r; effort returns with
                          the same factor; endpoints may be in different domains
  {\"power_map\": [effort_in_port, flow_in_port], \"matrix\": [[...], ...]}
                          a constant rectangular power map P: one row per
                          coefficient of the second port, one column per
                          coefficient of the first. Flow passes as -P f;
                          effort returns by the metric adjoint (P transpose
                          for ordinary ports). For a 3-D body's force port
                          and a scalar spring: matrix = [[nx, ny, nz]].
                          Entries carry the SI conversion and orientation;
                          heat, stream and signal ports are not power maps.
  {\"power_map\": [[{\"port\": effort_in_port, \"matrix\": [[...], ...]},
                  ...], flow_in_port]}
                          several effort-in ports onto one flow-in port, each
                          with its own map P_k: the flow-in port receives
                          -sum P_k f_k, and each member the adjoint of its
                          own map (a hertz_patch's wrench from a body's
                          force and torque ports)
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
                          copy a signal output (heat_capacity.temperature or
                          velocity_source.span, a published mass observation,
                          or a control element's output) to matching signal inputs
                          (a network's temperature<e>, slider.span_a/span_b,
                          modulated_transformer.ratio,
                          modulated_restriction.opening, input<i> or command);
                          it carries no energy
  {\"prescribe\": port, <quantity>: waveform}
                          hold a port's input from outside: velocity_m_s for a
                          spring, damper or force end (0 anchors it), force_n
                          for a mass; angular_velocity_rad_s and torque_n_m
                          about a shaft; current_a and voltage_v in a loop (0
                          current grounds a series element's far end); voltage_v
                          grounds a resistor_branch end; current_a into a
                          node_capacitor or nonlinear_capacitor increases its charge;
                          pressure_pa holds a hydraulic_resistor or carreau_pipe
                          end; flow_m3_s into an accumulator increases its volume;
                          flow_m3_s on a gas_cavity's p compresses it;
                          temperature_k for a damper, resistor, conductance,
                          heater, heat node body or carreau_pipe heat port,
                          heat_rate_w for a heat
                          capacity or a heat node's source; for a
                          network of rigid bodies, force_n or torque_n_m as
                          an array of three waveforms ([0, 0, 0] leaves a
                          body free or an actuator idle), and temperature_k
                          for a temperature<e> input held from outside; for a
                          gas volume, heat_rate_w on heat, and stream on its
                          stream port as an array of waveforms (mass rate
                          [kg/s], enthalpy rate [W], then each carried rate
                          [kg/s], all into the volume; zeros close it); for a
                          restriction's end held from outside, contents
                          (pressure [Pa], specific enthalpy [J/kg], then each
                          carried mass fraction); ratio for a
                          modulated_transformer's ratio signal; opening for
                          modulated_restriction.opening; value for
                          affine/sampled_controller inputs, ideal_switch.closed
                          or controlled-source commands (the latter also accept
                          their SI field name).
  A prescribed value, valve opening and speed_controller input are point
  samples at the step's start, matching signal_source. SI-named controlled-source commands
  and physical prescriptions use interval means.
  On slider.span_a/span_b, prescribe velocity_m_s as the same waveform
  moving that end; its two endpoint values are supplied, not its step mean.
  Prescribed table knots are scheduled automatically, including vector and
  stream inputs and slider endpoint spans. The kernel ends a step at the next
  knot; smooth waveforms still need a suitable step size or adaptive tolerances.

For a two-ended spring or damper, a is the end at mass A and b the end at
mass B. A spring's extension is x_B - x_A less its unloaded length; a series
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
    let options = match SceneOptions::parse(args) {
        Ok(options) => options,
        Err(message) => return refusal(exit::USAGE, "system-usage", &message, json),
    };
    let path = options.scene;
    let text = match read(path) {
        Ok(text) => text,
        Err(message) => return refusal(exit::INPUT, "system-input", &message, json),
    };
    let mut scene = match SystemScene::from_text(&text, path) {
        Ok(scene) => scene,
        Err(error) => return refusal(error.exit_code, error.code, &error.message, json),
    };
    if let Some((resume, root)) = options.resume.zip(options.resume_root) {
        let bytes = match read_checkpoint(resume) {
            Ok(bytes) => bytes,
            Err(message) => {
                return refusal(exit::INPUT, "system-checkpoint-file", &message, json);
            }
        };
        if let Err(error) = scene.resume_checkpoint(&bytes, root) {
            return refusal(error.exit_code, error.code, &error.message, json);
        }
    }
    let until_s = options.until_s.unwrap_or(scene.duration_s());
    if until_s < scene.time_s() || until_s > scene.duration_s() {
        return refusal(
            exit::INPUT,
            "system-until",
            &format!(
                "--until must be between the accepted time {} s and scene duration {} s",
                scene.time_s(),
                scene.duration_s()
            ),
            json,
        );
    }
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let key = StreamKey {
            seed: 0,
            kernel_id: 0x5953,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        // With --whole-steps every interval has the ceiling a run with no
        // pause gives it, the scene's horizon, so no step is cut at `until_s`.
        let ceiling_s = if options.whole_steps {
            scene.duration_s()
        } else {
            until_s
        };
        while !scene.is_finished() && scene.time_s() < until_s {
            scene.advance(ceiling_s - scene.time_s(), &cx)?;
        }
        Ok::<(), SystemSceneError>(())
    });
    if let Err(error) = result {
        return refusal(error.exit_code, error.code, &error.message, json);
    }
    let saved_root = if let Some(destination) = options.checkpoint {
        let (bytes, root) = match scene.checkpoint() {
            Ok(checkpoint) => checkpoint,
            Err(error) => return refusal(error.exit_code, error.code, &error.message, json),
        };
        if let Err(message) = write_checkpoint(destination, &bytes) {
            return refusal(exit::INPUT, "system-checkpoint-file", &message, json);
        }
        Some(fs_blake3::ContentHash(root).to_hex())
    } else {
        None
    };
    let saved = options.checkpoint.zip(saved_root.as_deref());
    let stdout = if json {
        scene.device.render_json(path, saved)
    } else {
        let mut report = scene.render_text();
        if let Some((destination, root)) = saved {
            let _ = writeln!(
                report,
                "\nCheckpoint: {}\nResume root: {root}",
                destination.display()
            );
        }
        report
    };
    CommandOutput {
        exit_code: exit::SUCCESS,
        stdout,
        stderr: String::new(),
    }
}

struct SceneOptions<'a> {
    scene: &'a Path,
    until_s: Option<f64>,
    whole_steps: bool,
    checkpoint: Option<&'a Path>,
    resume: Option<&'a Path>,
    resume_root: Option<[u8; 32]>,
}

impl<'a> SceneOptions<'a> {
    fn parse(args: &'a [OsString]) -> Result<Self, String> {
        let Some(scene) = args.first() else {
            return Err(HELP.into());
        };
        if scene.to_str().is_some_and(|text| text.starts_with('-')) {
            return Err(HELP.into());
        }
        let mut options = Self {
            scene: Path::new(scene),
            until_s: None,
            whole_steps: false,
            checkpoint: None,
            resume: None,
            resume_root: None,
        };
        let mut rest = args[1..].iter();
        while let Some(name) = rest.next() {
            // The one option with no value.
            if name.to_str() == Some("--whole-steps") {
                if options.whole_steps {
                    return Err(format!("--whole-steps is given twice\n{HELP}"));
                }
                options.whole_steps = true;
                continue;
            }
            let Some(value) = rest.next() else {
                return Err(format!("every system option requires a value\n{HELP}"));
            };
            match name.to_str() {
                Some("--until") if options.until_s.is_none() => {
                    let value = value
                        .to_str()
                        .and_then(|value| value.parse::<f64>().ok())
                        .filter(|value| value.is_finite() && *value >= 0.0)
                        .ok_or("--until must be a finite, nonnegative time in seconds")?;
                    options.until_s = Some(value);
                }
                Some("--checkpoint") if options.checkpoint.is_none() => {
                    options.checkpoint = Some(Path::new(value));
                }
                Some("--resume") if options.resume.is_none() => {
                    options.resume = Some(Path::new(value));
                }
                Some("--resume-root") if options.resume_root.is_none() => {
                    options.resume_root = Some(
                        value
                            .to_str()
                            .and_then(fs_blake3::ContentHash::from_hex)
                            .ok_or("--resume-root must contain exactly 64 hexadecimal digits")?
                            .0,
                    );
                }
                _ => {
                    return Err(format!(
                        "unknown or repeated system option {name:?}\n{HELP}"
                    ));
                }
            }
        }
        if options.resume.is_some() != options.resume_root.is_some() {
            return Err("--resume and --resume-root must be supplied together".into());
        }
        if options.whole_steps && options.until_s.is_none() {
            return Err("--whole-steps requires --until".into());
        }
        Ok(options)
    }
}

fn read_checkpoint(path: &Path) -> Result<Vec<u8>, String> {
    let read = || -> std::io::Result<Vec<u8>> {
        let file = std::fs::File::open(path)?;
        let mut bytes = Vec::new();
        file.take(snapshot::MAX_CHECKPOINT_BYTES as u64 + 1)
            .read_to_end(&mut bytes)?;
        Ok(bytes)
    };
    let bytes = read().map_err(|error| format!("checkpoint {}: {error}", path.display()))?;
    if bytes.len() > snapshot::MAX_CHECKPOINT_BYTES {
        return Err(format!(
            "checkpoint {} exceeds {} bytes",
            path.display(),
            snapshot::MAX_CHECKPOINT_BYTES
        ));
    }
    Ok(bytes)
}

fn write_checkpoint(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let write = || -> std::io::Result<()> {
        // A failed write cannot truncate a previous accepted checkpoint or
        // the scene itself. The caller receives a root only after sync succeeds.
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        file.write_all(bytes)?;
        file.sync_all()?;
        #[cfg(unix)]
        {
            let parent = path
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            std::fs::File::open(parent)?.sync_all()?;
        }
        Ok(())
    };
    write().map_err(|error| format!("checkpoint {}: {error}", path.display()))
}

/// A scene admission, probe or execution refusal, using the CLI's diagnostics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SystemSceneError {
    /// Stable process class, including [`exit::CANCELLED`] for cancellation.
    pub exit_code: u8,
    /// Stable diagnostic code; no output is printed by the library API.
    pub code: &'static str,
    /// Named input or execution detail in SI units where applicable.
    pub message: String,
}

impl std::fmt::Display for SystemSceneError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}

impl std::error::Error for SystemSceneError {}

impl From<Refusal> for SystemSceneError {
    fn from(error: Refusal) -> Self {
        let (exit_code, code, message) = match error {
            Refusal::Input(message) => (exit::INPUT, "system-input", message),
            Refusal::Admission(message) => (exit::REFUSED, "system-admission", message),
        };
        Self {
            exit_code,
            code,
            message,
        }
    }
}

/// An admitted native scene driven by the same kernel path as `system`.
///
/// Construction consumes a bounded snapshot of the scene's UTF-8 contents;
/// the scene file is never reopened. `source_path` is retained as the source
/// label and base location for relative material packs, which are read once
/// during admission. Checkpoints retain the accepted kernel generation,
/// adaptive driver memory and report samples, and bind the exact source and
/// executable. This v0 runtime is not the canonical assembly schema or a
/// browser binding.
/// The host supplies the cancellation gate and arena through `Cx`; this
/// facade adds no wall-clock or memory governor beyond the existing owners.
pub struct SystemScene {
    device: Device,
    source_path: PathBuf,
    source_identity: [u8; 32],
}

impl SystemScene {
    /// Admit at most 1 MiB of UTF-8 scene data without reading its source path.
    pub fn from_bytes(
        bytes: &[u8],
        source_path: impl AsRef<Path>,
    ) -> Result<Self, SystemSceneError> {
        if bytes.len() as u64 > MAX_SCENE_BYTES {
            return Err(Refusal::Input(format!("scene: at most {MAX_SCENE_BYTES} bytes")).into());
        }
        let text = std::str::from_utf8(bytes).map_err(|error| Refusal::Input(error.to_string()))?;
        Self::from_text(text, source_path)
    }

    /// Admit a retained text snapshot with the same bounds and rules as the CLI.
    pub fn from_text(text: &str, source_path: impl AsRef<Path>) -> Result<Self, SystemSceneError> {
        if text.len() as u64 > MAX_SCENE_BYTES {
            return Err(Refusal::Input(format!("scene: at most {MAX_SCENE_BYTES} bytes")).into());
        }
        Ok(Self {
            device: load_text(text, source_path.as_ref())?,
            source_path: source_path.as_ref().to_owned(),
            source_identity: *fs_blake3::hash_domain(
                "org.frankensim.fs-cli.system-scene.source.v1",
                text.as_bytes(),
            ).as_bytes(),
        })
    }

    /// The caller's source label and location; admission does not canonicalize it.
    #[must_use]
    pub fn source_path(&self) -> &Path {
        &self.source_path
    }

    /// Accepted physical time, never the time of an unpublished trial.
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.device.system.time_s()
    }

    /// The scene's declared final horizon.
    #[must_use]
    pub fn duration_s(&self) -> f64 {
        self.device.duration_s
    }

    /// Accepted kernel steps, including refinements accepted by the adaptive driver.
    #[must_use]
    pub fn steps(&self) -> u64 {
        self.device.system.steps()
    }

    /// Whether the existing final-horizon rounding allowance has been reached.
    #[must_use]
    pub fn is_finished(&self) -> bool {
        self.device.is_finished()
    }

    /// Read the accepted energy, entropy and stream-inventory accounts.
    #[must_use]
    pub fn ledger(&self) -> &Ledger {
        self.device.system.ledger()
    }

    /// Publish at most one driver interval, capped by the caller, scene and horizon.
    ///
    /// Returns true after acceptance and false when already finished. The
    /// caller's ceiling must be positive and finite. Adaptive retries and
    /// component events may shorten it further. A refusal leaves the accepted
    /// state, ledger, samples and driver statistics unchanged.
    pub fn advance(&mut self, max_interval_s: f64, cx: &Cx<'_>) -> Result<bool, SystemSceneError> {
        if !max_interval_s.is_finite() || max_interval_s <= 0.0 {
            return Err(SystemSceneError {
                exit_code: exit::INPUT,
                code: "system-interval",
                message: "maximum interval [s] must be positive and finite".into(),
            });
        }
        self.device
            .advance(max_interval_s, cx)
            .map_err(|error| self.run_error(error))
    }

    /// Continue from the accepted state to the declared horizon using the caller's context.
    pub fn run(&mut self, cx: &Cx<'_>) -> Result<(), SystemSceneError> {
        while self.advance(self.device.step_s, cx)? {}
        Ok(())
    }

    /// Read one report coordinate, e.g. `store.liquid_mass_fraction`.
    ///
    /// Unknown or ambiguous names refuse. `None` is the report's `null`:
    /// a stateless last-step output has no value before the first step.
    pub fn probe(&self, path: &str) -> Result<Option<f64>, SystemSceneError> {
        let refused = |reason| SystemSceneError {
            exit_code: exit::INPUT,
            code: "system-probe",
            message: format!("probe {path:?}: {reason}"),
        };
        let mut found = None;
        for element in &self.device.elements {
            let Some(field) = path
                .strip_prefix(element.name.as_str())
                .and_then(|tail| tail.strip_prefix('.'))
            else {
                continue;
            };
            if let Some((_, value)) = self
                .device
                .values(element)
                .into_iter()
                .find(|(name, _)| name == field)
            {
                if found.is_some() {
                    return Err(refused("ambiguous report coordinate"));
                }
                found = Some(value.is_finite().then_some(value));
            }
        }
        found.ok_or_else(|| refused("unknown report coordinate"))
    }

    /// Project the existing JSON report, including accepted samples and ledger accounts.
    #[must_use]
    pub fn render_json(&self) -> String {
        self.device.render_json(&self.source_path, None)
    }

    /// Project the existing human-readable report from the same accepted state.
    #[must_use]
    pub fn render_text(&self) -> String {
        self.device.render_text(&self.source_path)
    }

    fn run_error(&self, error: RunError) -> SystemSceneError {
        let (exit_code, code, message) = match error {
            RunError::Step(error) => (
                if matches!(error, SystemError::Cancelled) {
                    exit::CANCELLED
                } else {
                    exit::REFUSED
                },
                "system-step",
                format!(
                    "step {} at t = {} s: {error}",
                    self.steps() + 1,
                    self.time_s()
                ),
            ),
            RunError::Adaptive(error) => (
                if matches!(error, AdaptiveError::Step(SystemError::Cancelled)) {
                    exit::CANCELLED
                } else {
                    exit::REFUSED
                },
                "system-adaptive",
                format!(
                    "step {} at t = {} s: {error}",
                    self.steps() + 1,
                    self.time_s()
                ),
            ),
            RunError::TooManySteps => (
                exit::REFUSED,
                "system-steps",
                format!(
                    "the kernel took {} steps and reached t = {} s of {} s; components keep shortening the step",
                    self.steps(),
                    self.time_s(),
                    self.duration_s(),
                ),
            ),
        };
        SystemSceneError {
            exit_code,
            code,
            message,
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
    /// Pressure and volume flow of an incompressible liquid.
    Hydraulic,
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
    /// The power domain named by a generic source, converter or matrix law.
    fn parse(name: &str) -> Option<Self> {
        Some(match name {
            "translation" => Self::Translation,
            "rotation" => Self::Rotation,
            "electrical" => Self::Electrical,
            "hydraulic" => Self::Hydraulic,
            _ => return None,
        })
    }

    fn port_kind(self) -> PortKind {
        match self {
            Self::Translation => PortKind::MechanicalForceVelocity,
            Self::Rotation => PortKind::RotationalTorqueAngularVelocity,
            Self::Electrical => PortKind::ElectricalVoltageCurrent,
            Self::Hydraulic => PortKind::FluidPressureFlux,
        }
    }

    fn of(kind: PortKind) -> Option<Self> {
        [
            Self::Translation,
            Self::Rotation,
            Self::Electrical,
            Self::Hydraulic,
        ]
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
            // An accumulator is declared by its capacitance, the volume it
            // takes in per pascal; it holds that volume times its pressure.
            Self::Hydraulic => Names {
                inertance: "inertance_pa_s2_m3",
                displacement: "volume_m3",
                flow: "flow_m3_s",
                stiffness: "capacitance_m3_pa",
                cubic: "",
                extension: "volume_m3",
                resistance: "resistance_pa_s_m3",
                effort: "pressure_pa",
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
    /// A positive-definite resistance over several coordinates of one domain.
    ResistanceMatrix,
    /// A prescribed effort (force, torque, voltage).
    Source(Domain),
    /// An authored dense affine map over scalar signal inputs.
    AffineSignal,
    /// Replayable point-sampled signal waveforms on an accepted clock.
    SignalSource,
    /// An authored discrete recurrence with a sampled, held output.
    SampledController,
    /// A sampled speed policy with owner-defined PI and command limits.
    SpeedController,
    /// An ideal effort source commanded through a scalar signal.
    ControlledEffort,
    /// An ideal flow source, with flow positive into the source.
    ControlledFlow,
    /// Prescribed translation, publishing the velocities at the step's ends.
    VelocitySource,
    /// Constant-coefficient dry friction on the kernel's load line.
    Slider,
    /// A dry-friction law of `fs-tribo` on the kernel's load line.
    FrictionInterface,
    /// The same laws under a prescribed slip rate, `fs-tribo`'s sliding
    /// interface, with its heat on a heat port.
    SlidingInterface,
    /// The joint slip and spin friction of one circular Hertz patch,
    /// `fs-tribo`'s, on the kernel's load line, with its heat on a heat port.
    HertzPatch,
    /// A linear structure in damped modes, `fs-nlmodal`'s, acted on through
    /// its attachments.
    ModalStructure,
    /// A one-point normal power-law contact between two ends,
    /// `fs-dcontact`'s, on the kernel's load line.
    ContactInterface,
    /// Storage setting the effort of a parallel node: charge and voltage
    /// (node capacitor), volume and pressure (accumulator).
    NodeCapacitor(Domain),
    /// Flow from an effort difference, with a heat port.
    ResistorBranch(Domain),
    /// A prescribed flow, positive into the source.
    CurrentSource(Domain),
    /// An ideal one-way electrical relation on the kernel's load line.
    IdealDiode,
    /// An ideal electrical open/short relation with a binary signal input.
    IdealSwitch,
    /// The fs-material laminar pipe of a Carreau-Yasuda liquid, whose
    /// viscosity follows the temperature at its heat port.
    CarreauPipe,
    /// Charge storage governed by a source-owned V(q) curve.
    NonlinearCapacitor,
    /// Flux storage governed by a source-owned I(flux) curve.
    NonlinearInductor,
    Gyrator,
    /// A lossless power map whose scalar ratio comes from a signal.
    ModulatedTransformer,
    /// An authored lossless scalar N-port law in one power domain.
    DiracJunction,
    /// Authored linear port-Hamiltonian matrices, stepped exactly for held
    /// inputs by the fs-system owner, with a heat port.
    PortHamiltonian,
    /// Rigid bodies and viscous connectors, stepped by their fs-mbd owner.
    RigidNetwork,
    /// The same with Maxwell branches and actuators on the connectors.
    ViscoelasticNetwork,
    HeatCapacity,
    /// Fixed material mass on an authored sensible/latent enthalpy chart.
    EnthalpyBody,
    Conductance,
    /// Original temperature-dependent layer integrals and contact jump.
    LayerConductor,
    /// A thermal node of no capacity joined to bodies by conductances, with
    /// a heat source at the node: the plane a heat is generated on.
    HeatNode,
    /// Enthalpy cells, their links and boundaries, stepped by the fs-thermal
    /// timeline owner.
    ThermalTimeline,
    /// A box meshed in tetrahedra, `fs-conduction`'s backward-Euler
    /// finite-element body, with heat ports on its named faces.
    ConductionBody,
    /// `fs-conduction`'s contact between two traces of two conduction
    /// bodies, carried by the kernel's conductance over interface points.
    TraceContact,
    /// A well-mixed volume of perfect gas: the node streams meet at.
    GasVolume,
    /// A linear flow restriction between two gas volumes.
    Restriction,
    /// The same native stream law with a held dimensionless opening signal.
    ModulatedRestriction,
    Heater,
    /// A closed, uniform gas of NASA-9 thermochemistry from cited packs,
    /// `fs-gas`'s cavity, with a pressure-work port and a heat port.
    GasCavity,
}

impl Kind {
    const TABLE: [(&'static str, Self); 57] = [
        ("mass", Self::Inertance(Domain::Translation)),
        ("spring", Self::Compliance(Domain::Translation)),
        ("damper", Self::Resistance(Domain::Translation)),
        ("resistance_matrix", Self::ResistanceMatrix),
        ("force", Self::Source(Domain::Translation)),
        ("affine_signal", Self::AffineSignal),
        ("signal_source", Self::SignalSource),
        ("sampled_controller", Self::SampledController),
        ("speed_controller", Self::SpeedController),
        ("controlled_effort", Self::ControlledEffort),
        ("controlled_flow", Self::ControlledFlow),
        ("velocity_source", Self::VelocitySource),
        ("slider", Self::Slider),
        ("friction_interface", Self::FrictionInterface),
        ("sliding_interface", Self::SlidingInterface),
        ("hertz_patch", Self::HertzPatch),
        ("modal_structure", Self::ModalStructure),
        ("contact_interface", Self::ContactInterface),
        ("inertia", Self::Inertance(Domain::Rotation)),
        ("torsion_spring", Self::Compliance(Domain::Rotation)),
        ("rotary_damper", Self::Resistance(Domain::Rotation)),
        ("torque", Self::Source(Domain::Rotation)),
        ("inductor", Self::Inertance(Domain::Electrical)),
        ("capacitor", Self::Compliance(Domain::Electrical)),
        ("resistor", Self::Resistance(Domain::Electrical)),
        ("voltage", Self::Source(Domain::Electrical)),
        ("node_capacitor", Self::NodeCapacitor(Domain::Electrical)),
        ("resistor_branch", Self::ResistorBranch(Domain::Electrical)),
        ("current_source", Self::CurrentSource(Domain::Electrical)),
        ("ideal_diode", Self::IdealDiode),
        ("ideal_switch", Self::IdealSwitch),
        ("accumulator", Self::NodeCapacitor(Domain::Hydraulic)),
        (
            "hydraulic_resistor",
            Self::ResistorBranch(Domain::Hydraulic),
        ),
        ("flow_source", Self::CurrentSource(Domain::Hydraulic)),
        ("pressure", Self::Source(Domain::Hydraulic)),
        ("carreau_pipe", Self::CarreauPipe),
        ("nonlinear_capacitor", Self::NonlinearCapacitor),
        ("nonlinear_inductor", Self::NonlinearInductor),
        ("gyrator", Self::Gyrator),
        ("modulated_transformer", Self::ModulatedTransformer),
        ("dirac_junction", Self::DiracJunction),
        ("port_hamiltonian", Self::PortHamiltonian),
        ("rigid_network", Self::RigidNetwork),
        ("viscoelastic_network", Self::ViscoelasticNetwork),
        ("heat_capacity", Self::HeatCapacity),
        ("enthalpy_body", Self::EnthalpyBody),
        ("conductance", Self::Conductance),
        ("layer_conductor", Self::LayerConductor),
        ("heat_node", Self::HeatNode),
        ("thermal_timeline", Self::ThermalTimeline),
        ("conduction_body", Self::ConductionBody),
        ("trace_contact", Self::TraceContact),
        ("gas_volume", Self::GasVolume),
        ("restriction", Self::Restriction),
        ("modulated_restriction", Self::ModulatedRestriction),
        ("heater", Self::Heater),
        ("gas_cavity", Self::GasCavity),
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
    /// A controlled source driven by a signal or prescribed through `value`
    /// uses a held sample; SI-named prescriptions retain interval-mean semantics.
    point_sampled_command: bool,
    /// Initial temperature of a heat capacity, used as the starting value of
    /// the temperature inputs joined to it.
    temperature_k: Option<f64>,
    /// The same for single ports of an element whose ports start at
    /// different temperatures (a thermal timeline's cells); a port named here
    /// takes this value in place of `temperature_k`.
    port_temperatures_k: Vec<(String, f64)>,
    /// A node capacitor's capacitance [F] or an accumulator's [m³/Pa]: its
    /// endpoint charge or volume over this is its endpoint voltage or pressure.
    node_capacity: Option<f64>,
    /// Original laws for endpoint probes; prepared data are shared on clone.
    charge_voltage_curve: Option<CapacitorChargeCurve>,
    flux_current_curve: Option<InductorFluxCurve>,
    /// The bodies of a network, for reading velocities from momenta, and for
    /// each of its connectors the number of its Maxwell branches.
    bodies: Vec<MassProperties>,
    branches: Vec<usize>,
    /// The gas of a gas volume, for reading pressure and temperature from
    /// what it holds.
    gas: Option<GasVolume>,
    /// Named prepared sources actually selected by this component.
    material_sources: Vec<String>,
    /// Each network body declared as a solid: its index and report object.
    solids: Vec<(usize, String)>,
    /// The state coordinates of a port-Hamiltonian element, reported and
    /// error-controlled as x0, x1, ...
    coordinates: usize,
    /// A Hertz patch's radius [m], which scales its port's third coefficient:
    /// its moment is the radius times that coefficient's effort.
    radius_m: Option<f64>,
    /// What the scene declared for a trace contact: the two traces it was
    /// computed for and its resistance with its source.
    contact: Option<conduction::Declared>,
    /// A conduction body whose state is specific enthalpy: the paired
    /// material whose chart reads it as temperature.
    enthalpy: Option<&'static MaterialEnthalpyTransport>,
    /// A modal structure's angular frequencies [rad/s], one per mode, for the
    /// energy norm of its error account.
    modal: Option<&'static [f64]>,
    /// A card-backed resistor's node, for reporting its resistance at the
    /// temperature its heat port read.
    ohmic: Option<MaterialOhmicNode>,
    /// A card-backed pipe, for reporting its viscosity at the temperature
    /// its heat port read.
    laminar: Option<MaterialLaminarPipe>,
}

/// A controlled source's scalar command uses the SI unit of what its physical
/// port returns. Signal connections carry those authored numbers unchanged.
fn controlled_quantity(element: &Element) -> Option<&'static str> {
    let flow = match element.kind {
        Kind::ControlledEffort => false,
        Kind::ControlledFlow => true,
        _ => return None,
    };
    let domain = element.ports.iter().find_map(|port| match port.carrier() {
        Carrier::Power(kind) => Domain::of(kind),
        _ => None,
    })?;
    Some(if flow {
        domain.names().flow
    } else {
        domain.names().effort
    })
}

struct Device {
    system: System,
    elements: Vec<Element>,
    materials: materials::Materials,
    connections: usize,
    duration_s: f64,
    step_s: f64,
    sample_every_s: Option<f64>,
    samples: Vec<(f64, Vec<(String, f64)>)>,
    next_sample_s: f64,
    iterations: u64,
    shortenings: u64,
    adaptive: Option<AdaptiveSettings>,
}

#[derive(Clone)]
enum Coordinate {
    Position,
    Velocity,
    Scalar,
    Enthalpy,
    SourceWork,
    Network(adaptive_metric::NetworkCoordinate),
    /// One coordinate of a port-Hamiltonian element's state vector.
    State(usize),
    /// What a gas volume holds.
    Held(Held),
    /// A gas cavity's volume, and the temperature of its gas.
    CavityVolume,
    CavityTemperature,
    /// A streaming cavity's native extensive inventory and caloric energy.
    CavityMass,
    CavityEnergy,
    /// Stable stream identity, independent of the active mixture's ordering.
    CavitySpecies(SpeciesId),
}

/// A coordinate of what a gas volume holds: its mass, its internal energy,
/// its volume, or the mass of one carried quantity.
#[derive(Clone, Copy)]
enum Held {
    Mass,
    Energy,
    Volume,
    Carried(usize),
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
            Coordinate::Enthalpy => checkpoint.state::<EquilibriumPhaseState>(self.node)
                .map_or(f64::NAN, |state| state.specific_enthalpy_j_kg()),
            Coordinate::SourceWork => checkpoint.state::<SourceAccount>(self.node)
                .map_or(f64::NAN, |state| state.energy_in_j),
            Coordinate::Network(_) => f64::NAN,
            Coordinate::State(index) => checkpoint
                .state::<Vec<f64>>(self.node)
                .and_then(|state| state.get(index))
                .copied()
                .unwrap_or(f64::NAN),
            Coordinate::Held(held) => checkpoint
                .state::<GasContents>(self.node)
                .and_then(|gas| match held {
                    Held::Mass => Some(gas.mass_kg),
                    Held::Energy => Some(gas.energy_j),
                    Held::Volume => Some(gas.volume_m3),
                    Held::Carried(index) => gas.carried_kg.get(index).copied(),
                })
                .unwrap_or(f64::NAN),
            Coordinate::CavityVolume => checkpoint
                .state::<ClosedGasParcel>(self.node)
                .map_or(f64::NAN, ClosedGasParcel::volume_m3),
            Coordinate::CavityTemperature => checkpoint
                .state::<ClosedGasParcel>(self.node)
                .map_or(f64::NAN, |gas| gas.state().state().temperature),
            Coordinate::CavityMass => checkpoint
                .state::<ClosedGasParcel>(self.node)
                .map_or(f64::NAN, ClosedGasParcel::mass_kg),
            Coordinate::CavityEnergy => checkpoint
                .state::<ClosedGasParcel>(self.node)
                .map_or(f64::NAN, gas::energy_j),
            Coordinate::CavitySpecies(ref species) => checkpoint
                .state::<ClosedGasParcel>(self.node)
                .and_then(|gas| gas.species_mass_kg(species).ok())
                .unwrap_or(f64::NAN),
        }
    }
}

/// A gas volume's error accounts: its mass, its internal energy and its
/// volume, and the mass of each quantity its stream carries.
fn gas_volume_accounts(element: &Element) -> Vec<(String, NodeId, Coordinate)> {
    let carried: usize = element.ports.iter().map(PortSpec::carried_quantities).sum();
    let held = [
        ("mass_kg".to_owned(), Held::Mass),
        ("energy_j".to_owned(), Held::Energy),
        ("volume_m3".to_owned(), Held::Volume),
    ];
    let carried = (0..carried).map(|index| (format!("carried{index}_kg"), Held::Carried(index)));
    held.into_iter()
        .chain(carried)
        .map(|(field, held)| {
            let name = format!("{}.{field}", element.name);
            (name, element.node, Coordinate::Held(held))
        })
        .collect()
}

/// Source-backed gas keeps its native temperature and volume coordinates.
/// Material transfer additionally controls every extensive inventory and
/// the caloric energy, using the stream catalog even after species depletion.
fn gas_cavity_accounts(element: &Element) -> Vec<(String, NodeId, Coordinate)> {
    let mut coordinates = vec![
        ("volume_m3".to_owned(), Coordinate::CavityVolume),
        ("temperature_k".to_owned(), Coordinate::CavityTemperature),
    ];
    if let Some(basis) = gas::stream_basis(&element.ports) {
        coordinates.push(("mass_kg".to_owned(), Coordinate::CavityMass));
        coordinates.push(("energy_j".to_owned(), Coordinate::CavityEnergy));
        for (index, quantity) in basis.quantities().iter().enumerate() {
            if let StreamQuantity::SpeciesMass(species) = quantity {
                coordinates.push((
                    format!("species{index}_mass_kg"),
                    Coordinate::CavitySpecies(species.clone()),
                ));
            }
        }
    }
    coordinates
        .into_iter()
        .map(|(field, coordinate)| {
            (format!("{}.{field}", element.name), element.node, coordinate)
        })
        .collect()
}

/// The refusal of an element whose kind has no error account, with the
/// reason for its kind.
fn no_error_account(element: &Element) -> String {
    let why = match element.kind {
        Kind::SlidingInterface => "its law's state is not given as numbers",
        Kind::ThermalTimeline => "its cells are born and joined at events",
        Kind::SampledController => {
            "its discrete memory and sample events have no step-doubling error metric here"
        }
        Kind::SpeedController => {
            "its sample count, phase and held command form a discrete policy, not a continuous error coordinate"
        }
        Kind::IdealDiode | Kind::IdealSwitch => {
            "its switching relation has no event-resolved step-doubling error metric here"
        }
        // A slider, a friction interface, a Hertz patch, a contact interface.
        _ => {
            "its force jumps where it sticks, slips or strikes, and an estimate by step \
             doubling does not fall with the step there"
        }
    };
    format!(
        "adaptive: {} is a {}, whose state has no error account here ({why}); \
         run it with a fixed step_s",
        element.name,
        element.kind.name()
    )
}

/// The one error account over the whole state of a conduction body or of a
/// modal structure, with its name: the body's field, in kelvin or in the
/// specific enthalpy an enthalpy body holds, or the structure's energy norm.
fn whole_state_account(element: &Element) -> (&'static str, Coordinate) {
    let (name, account) = match (element.modal, element.enthalpy) {
        (Some(omegas), _) => (
            "energy_norm_sqrt_j",
            adaptive_metric::NetworkCoordinate::Modes { omegas },
        ),
        (None, Some(_)) => (
            "specific_enthalpy_j_kg",
            adaptive_metric::NetworkCoordinate::Field,
        ),
        (None, None) => ("temperature_k", adaptive_metric::NetworkCoordinate::Field),
    };
    (name, Coordinate::Network(account))
}

/// The supported scene order for kinds that bound `adaptive.order`. Physical
/// owners' step orders are measured by halving a fixed step. A
/// conduction body steps by backward Euler. A modal structure is exact for a
/// held effort and holds one that varies within a step, which leaves the
/// coupled step at order two. A gas volume steps by the trapezoidal rule, and
/// a closed gas cavity measures order two. The native streaming cavity's
/// explicit material transfer and subsequent heat/work solve are first order.
/// Point-sampled commands and accepted-state observations conservatively limit
/// control scenes to order one.
/// None: the order is the scene's to establish.
fn step_order(element: &Element) -> Option<u32> {
    if element.point_sampled_command {
        return Some(1);
    }
    match element.kind {
        Kind::ConductionBody
        | Kind::AffineSignal
        | Kind::SignalSource
        | Kind::ModulatedRestriction => Some(1),
        Kind::Inertance(_)
            if element
                .ports
                .iter()
                .any(|port| port.carrier() == Carrier::Signal) =>
        {
            Some(1)
        }
        Kind::GasCavity => Some(if gas::stream_basis(&element.ports).is_some() {
            1
        } else {
            2
        }),
        Kind::ModalStructure | Kind::GasVolume => Some(2),
        _ => None,
    }
}

/// The step-doubling driver trusts the declared order when it grows the step,
/// so a scene may not declare an order above that of a kind it accounts.
fn declared_order_fits(order: f64, elements: &[Element]) -> Result<(), String> {
    for element in elements {
        if let Some(own) = step_order(element).filter(|own| order > f64::from(*own)) {
            if matches!(
                element.kind,
                Kind::AffineSignal
                    | Kind::SignalSource
                    | Kind::ControlledEffort
                    | Kind::ControlledFlow
                    | Kind::ModulatedRestriction
                    | Kind::Inertance(_)
            ) {
                return Err(format!(
                    "adaptive.order: {} uses held start-of-step control signals; the supported order is 1, but the scene declares order {order}",
                    element.name,
                ));
            }
            return Err(format!(
                "adaptive.order: {} is a {}, whose step is of order {own}; the scene declares order {order}",
                element.name,
                element.kind.name()
            ));
        }
    }
    Ok(())
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
            if let Coordinate::Network(coordinate) = account.coordinate {
                coordinate.include(
                    &mut norm,
                    account.node,
                    initial,
                    coarse,
                    refined,
                    account.tolerance,
                );
                continue;
            }
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
        .ok_or_else(|| format!("{what}: a waveform is a number, sine, table or held_table"))?;
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
        [(key, table)] if key == "table" || key == "held_table" => {
            let rows = table
                .as_array()
                .ok_or_else(|| format!("{what}.{key}: expected an array of [time_s, value]"))?;
            let knots = rows
                .iter()
                .map(|row| match row.as_array() {
                    Some([t, v]) => t.as_f64().zip(v.as_f64()),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| format!("{what}.{key}: every row is [time_s, value]"))?;
            let waveform = if key == "held_table" {
                Waveform::held_table(knots)
            } else {
                Waveform::table(knots)
            };
            waveform.map_err(|e| format!("{what}: {e}"))
        }
        _ => Err(format!(
            "{what}: a waveform object has exactly one key, sine, table or held_table"
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

/// Lower original authored or source-resolved calorics without freezing Cp.
fn enthalpy_body(
    definition: &J,
    materials: &materials::Materials,
    what: &str,
) -> Result<(UniformEnthalpyBody, EquilibriumPhaseState), String> {
    if definition.get("material").is_some() {
        fields(
            definition,
            &["kind", "mass_kg", "specific_enthalpy_j_kg", "material"],
            what,
        )?;
        let name = definition
            .str_field("material")
            .ok_or_else(|| format!("{what}.material: expected a material name"))?;
        let source = materials
            .get(name)
            .ok_or_else(|| format!("{what}.material: unknown material {name:?}"))?;
        let initial = source
            .caloric
            .state_at_specific_enthalpy(number(definition, "specific_enthalpy_j_kg", what)?)
            .map_err(|error| format!("{what}.specific_enthalpy_j_kg: {error}"))?;
        let body = UniformEnthalpyBody::new(
            source.caloric.clone(),
            positive(definition, "mass_kg", what)?,
        )
        .map_err(|error| format!("{what}: {error}"))?;
        return Ok((body, initial));
    }
    fields(
        definition,
        &[
            "kind",
            "mass_kg",
            "specific_enthalpy_j_kg",
            "phase",
            "enthalpy_curve",
        ],
        what,
    )?;
    let curve = authored_enthalpy_curve(definition, what)?;
    let initial = curve
        .state_at_specific_enthalpy(number(definition, "specific_enthalpy_j_kg", what)?)
        .map_err(|error| format!("{what}.specific_enthalpy_j_kg: {error}"))?;
    let body = UniformEnthalpyBody::new(curve, positive(definition, "mass_kg", what)?)
        .map_err(|error| format!("{what}: {error}"))?;
    Ok((body, initial))
}

/// An authored sensible and latent chart: `phase` and `enthalpy_curve`.
fn authored_enthalpy_curve(
    definition: &J,
    what: &str,
) -> Result<EquilibriumEnthalpyPhaseCurve, String> {
    let phase = definition
        .str_field("phase")
        .ok_or_else(|| format!("{what}.phase: expected solid, liquid or solid-liquid"))?;
    if !matches!(phase, "solid" | "liquid" | "solid-liquid") {
        return Err(format!(
            "{what}.phase: expected solid, liquid or solid-liquid"
        ));
    }
    let rows = definition
        .get("enthalpy_curve")
        .and_then(J::as_array)
        .filter(|rows| (2..=MAX_ENTHALPY_KNOTS).contains(&rows.len()))
        .ok_or_else(|| format!("{what}.enthalpy_curve: expected 2..={MAX_ENTHALPY_KNOTS} knots"))?;
    let mut knots = Vec::with_capacity(rows.len());
    // This identifies the authored law, not a verified material card. Use
    // physical values in a fixed order, independent of JSON whitespace.
    let mut authored = phase.as_bytes().to_vec();
    authored.push(0);
    for (index, row) in rows.iter().enumerate() {
        let Some([h, t, liquid, density]) = row.as_array() else {
            return Err(format!(
                "{what}.enthalpy_curve[{index}]: expected [h_j_kg, temperature_k, \
                 liquid_mass_fraction, bulk_density_kg_m3]"
            ));
        };
        let mut values = [0.0; 4];
        for (slot, value) in values.iter_mut().zip([h, t, liquid, density]) {
            *slot = value
                .as_f64()
                .filter(|value| value.is_finite())
                .ok_or_else(|| {
                    format!("{what}.enthalpy_curve[{index}]: four finite numbers required")
                })?;
            // The sign of an authored zero is not a distinct material.
            if *slot == 0.0 {
                *slot = 0.0;
            }
            authored.extend_from_slice(&slot.to_le_bytes());
        }
        knots.push(EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: values[0],
            temperature_k: values[1],
            liquid_mass_fraction: values[2],
            bulk_density_kg_m3: values[3],
        });
    }
    let identity =
        fs_blake3::hash_domain("frankensim.system.authored_enthalpy_curve.v0", &authored);
    match phase {
        "solid" => EquilibriumEnthalpyPhaseCurve::try_single_phase(
            identity,
            SolidLiquidPhase::Solid,
            knots,
        ),
        "liquid" => EquilibriumEnthalpyPhaseCurve::try_single_phase(
            identity,
            SolidLiquidPhase::Liquid,
            knots,
        ),
        _ => EquilibriumEnthalpyPhaseCurve::try_new(identity, knots),
    }
    .map_err(|error| format!("{what}.enthalpy_curve: {error}"))
}

fn conductivity_curve(definition: &J, key: &str, what: &str) -> Result<ConductivityCurve, String> {
    let rows = definition
        .get(key)
        .and_then(J::as_array)
        .filter(|rows| (2..=MAX_CONDUCTIVITY_KNOTS).contains(&rows.len()))
        .ok_or_else(|| format!("{what}.{key}: expected 2..={MAX_CONDUCTIVITY_KNOTS} knots"))?;
    let knots = rows
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let values = match row.as_array() {
                Some([temperature, conductivity]) => {
                    temperature.as_f64().zip(conductivity.as_f64())
                }
                _ => None,
            };
            match values {
                Some((temperature, conductivity))
                    if temperature.is_finite() && conductivity.is_finite() =>
                {
                    Ok([temperature, conductivity])
                }
                _ => Err(format!(
                    "{what}.{key}[{index}]: expected finite [temperature_k, k_w_m_k]"
                )),
            }
        })
        .collect::<Result<Vec<_>, String>>()?;
    ConductivityCurve::new(knots).map_err(|error| format!("{what}.{key}: {error}"))
}

/// Preserve the physical geometry and original k(T) integrals at the ports.
fn layer_conductor(
    definition: &J,
    materials: &materials::Materials,
    what: &str,
) -> Result<LayerConductor, String> {
    fields(
        definition,
        &[
            "kind",
            "conductivity_a",
            "conductivity_b",
            "materials",
            "geometry",
            "contact_resistance_m2_k_w",
            "local_solver",
        ],
        what,
    )?;
    let curves = if let Some(names) = definition.get("materials") {
        if definition.get("conductivity_a").is_some() || definition.get("conductivity_b").is_some()
        {
            return Err(format!(
                "{what}.materials: cannot combine material sources with authored conductivity tables"
            ));
        }
        let Some([a, b]) = names.as_array() else {
            return Err(format!("{what}.materials: expected two material names"));
        };
        let curve = |entry: &J, index| {
            let name = entry
                .as_str()
                .ok_or_else(|| format!("{what}.materials[{index}]: expected a material name"))?;
            let source = materials
                .get(name)
                .ok_or_else(|| format!("{what}.materials[{index}]: unknown material {name:?}"))?;
            source.conductivity.clone().ok_or_else(|| {
                format!("{what}.materials[{index}]: material {name:?} has no prepared conductivity")
            })
        };
        [curve(a, 0)?, curve(b, 1)?]
    } else {
        [
            conductivity_curve(definition, "conductivity_a", what)?,
            conductivity_curve(definition, "conductivity_b", what)?,
        ]
    };
    let geometry_name = format!("{what}.geometry");
    let geometry = definition
        .get("geometry")
        .ok_or_else(|| format!("{geometry_name}: required"))?;
    let pair = |key: &str| {
        let values = match geometry.get(key).and_then(J::as_array) {
            Some([a, b]) => a.as_f64().zip(b.as_f64()),
            _ => None,
        };
        match values {
            Some((a, b)) if a.is_finite() && b.is_finite() => Ok([a, b]),
            _ => Err(format!(
                "{geometry_name}.{key}: expected two finite numbers"
            )),
        }
    };
    let geometry = match geometry.str_field("kind") {
        Some("planar") => {
            fields(
                geometry,
                &["kind", "distances_m", "area_m2"],
                &geometry_name,
            )?;
            LayerGeometry::planar(
                pair("distances_m")?,
                positive(geometry, "area_m2", &geometry_name)?,
            )
        }
        Some("cylindrical") => {
            fields(
                geometry,
                &["kind", "radii_m", "contact_radius_m", "length_m"],
                &geometry_name,
            )?;
            LayerGeometry::cylindrical(
                pair("radii_m")?,
                positive(geometry, "contact_radius_m", &geometry_name)?,
                positive(geometry, "length_m", &geometry_name)?,
            )
        }
        Some("spherical") => {
            fields(
                geometry,
                &["kind", "radii_m", "contact_radius_m"],
                &geometry_name,
            )?;
            LayerGeometry::spherical(
                pair("radii_m")?,
                positive(geometry, "contact_radius_m", &geometry_name)?,
            )
        }
        _ => {
            return Err(format!(
                "{geometry_name}.kind: expected planar, cylindrical or spherical"
            ));
        }
    }
    .map_err(|error| format!("{geometry_name}: {error}"))?;
    let resistance = number_or(definition, "contact_resistance_m2_k_w", 0.0, what)?;
    if !resistance.is_finite() || resistance < 0.0 {
        return Err(format!(
            "{what}.contact_resistance_m2_k_w: must be finite and not negative"
        ));
    }
    let mut config = LayerConductionConfig {
        max_iterations: 80,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12,
    };
    if let Some(local) = definition.get("local_solver") {
        let local_name = format!("{what}.local_solver");
        fields(
            local,
            &[
                "max_iterations",
                "absolute_heat_tolerance_w",
                "relative_heat_tolerance",
            ],
            &local_name,
        )?;
        let iterations = number_or(local, "max_iterations", 80.0, &local_name)?;
        if !(1.0..=1000.0).contains(&iterations) || iterations.fract() != 0.0 {
            return Err(format!(
                "{local_name}.max_iterations: must be a whole number from 1 to 1000"
            ));
        }
        config.max_iterations = iterations as usize;
        config.absolute_heat_tolerance_w = number_or(
            local,
            "absolute_heat_tolerance_w",
            config.absolute_heat_tolerance_w,
            &local_name,
        )?;
        config.relative_heat_tolerance = number_or(
            local,
            "relative_heat_tolerance",
            config.relative_heat_tolerance,
            &local_name,
        )?;
        for (key, value) in [
            (
                "absolute_heat_tolerance_w",
                config.absolute_heat_tolerance_w,
            ),
            ("relative_heat_tolerance", config.relative_heat_tolerance),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(format!(
                    "{local_name}.{key}: must be finite and not negative"
                ));
            }
        }
    }
    LayerConductor::new(curves, geometry, resistance, config)
        .map_err(|error| format!("{what}: {error}"))
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

/// A finite matrix in the row-major order used by the kernel's owners.
fn coefficient_matrix(
    value: Option<&J>,
    row_count: usize,
    column_count: usize,
    what: &str,
) -> Result<Vec<f64>, String> {
    let rows = value
        .and_then(J::as_array)
        .filter(|rows| rows.len() == row_count)
        .ok_or_else(|| format!("{what}: matrix row count must be {row_count}"))?;
    let entries = row_count
        .checked_mul(column_count)
        .ok_or_else(|| format!("{what}: matrix dimensions are too large"))?;
    let mut matrix = Vec::with_capacity(entries);
    for (row_index, row) in rows.iter().enumerate() {
        let values = row
            .as_array()
            .filter(|values| values.len() == column_count)
            .ok_or_else(|| {
                format!("{what}: matrix[{row_index}] column count must be {column_count}")
            })?;
        for (column, value) in values.iter().enumerate() {
            let coefficient = value
                .as_f64()
                .filter(|value| value.is_finite())
                .ok_or_else(|| {
                    format!("{what}: matrix[{row_index}][{column}] must be a finite number")
                })?;
            matrix.push(coefficient);
        }
    }
    Ok(matrix)
}

fn resistance_matrix(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(definition, &["kind", "domain", "matrix"], what)?;
    let domain = definition
        .str_field("domain")
        .and_then(Domain::parse)
        .ok_or_else(|| format!("{what}: domain is one of translation, rotation, electrical, hydraulic"))?;
    let rows = definition
        .get("matrix")
        .and_then(J::as_array)
        .filter(|rows| (1..=MAX_RESISTANCE_COORDINATES).contains(&rows.len()))
        .ok_or_else(|| {
            format!("{what}: matrix must have 1 to {MAX_RESISTANCE_COORDINATES} square rows")
        })?;
    let matrix = coefficient_matrix(definition.get("matrix"), rows.len(), rows.len(), what)?;
    let resistance = Resistance::new(domain.port_kind(), matrix)
        .map_err(|error| Refusal::Input(format!("{what}: {error}")))?;
    add(builder, name, resistance, ())
}

/// The nonempty rectangular shape that admits a bounded signal map.
fn signal_shape(definition: &J, key: &str, what: &str) -> Result<(usize, usize), String> {
    let rows = definition
        .get(key)
        .and_then(J::as_array)
        .filter(|rows| (1..=MAX_SIGNAL_COORDINATES).contains(&rows.len()))
        .ok_or_else(|| {
            format!("{what}.{key}: expected 1 to {MAX_SIGNAL_COORDINATES} matrix rows")
        })?;
    let columns = rows[0]
        .as_array()
        .filter(|row| (1..=MAX_SIGNAL_COORDINATES).contains(&row.len()))
        .ok_or_else(|| {
            format!("{what}.{key}: expected 1 to {MAX_SIGNAL_COORDINATES} matrix columns")
        })?
        .len();
    Ok((rows.len(), columns))
}

fn signal_coefficients(value: Option<&J>, count: usize, what: &str) -> Result<Vec<f64>, String> {
    value
        .and_then(J::as_array)
        .filter(|values| values.len() == count)
        .ok_or_else(|| format!("{what}: expected {count} coefficients"))?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            value
                .as_f64()
                .filter(|value| value.is_finite())
                .ok_or_else(|| format!("{what}[{index}]: expected a finite number"))
        })
        .collect()
}

fn affine_signal(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(definition, &["kind", "matrix", "bias"], what)?;
    let (outputs, inputs) = signal_shape(definition, "matrix", what)?;
    let matrix = coefficient_matrix(definition.get("matrix"), outputs, inputs, what)?;
    let bias = signal_coefficients(definition.get("bias"), outputs, &format!("{what}.bias"))?;
    let signal = AffineSignal::new(inputs, outputs, matrix, bias)
        .map_err(|error| format!("{what}: {error}"))?;
    add(builder, name, signal, ())
}

fn signal_source(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(definition, &["kind", "waveforms"], what)?;
    let values = definition
        .get("waveforms")
        .and_then(J::as_array)
        .filter(|values| (1..=MAX_SIGNAL_COORDINATES).contains(&values.len()))
        .ok_or_else(|| {
            format!("{what}.waveforms: expected 1 to {MAX_SIGNAL_COORDINATES} waveforms")
        })?;
    let waves = values
        .iter()
        .enumerate()
        .map(|(index, value)| waveform(value, &format!("{what}.waveforms[{index}]")))
        .collect::<Result<Vec<_>, _>>()?;
    let source = SignalSource::new(waves).map_err(|error| format!("{what}: {error}"))?;
    add(builder, name, source, SignalClock::default())
}

fn sampled_controller(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(
        definition,
        &["kind", "sample_period_s", "a", "b", "c", "d", "initial"],
        what,
    )?;
    let memory = definition
        .get("initial")
        .and_then(J::as_array)
        .filter(|values| values.len() <= MAX_SIGNAL_COORDINATES)
        .ok_or_else(|| {
            format!("{what}.initial: expected 0 to {MAX_SIGNAL_COORDINATES} coordinates")
        })?
        .len();
    let initial = signal_coefficients(
        definition.get("initial"),
        memory,
        &format!("{what}.initial"),
    )?;
    let (outputs, inputs) = signal_shape(definition, "d", what)?;
    let matrix = |key: &str, rows, columns| {
        coefficient_matrix(definition.get(key), rows, columns, &format!("{what}.{key}"))
    };
    let controller = SampledController::new(
        positive(definition, "sample_period_s", what)?,
        memory,
        inputs,
        outputs,
        matrix("a", memory, memory)?,
        matrix("b", memory, inputs)?,
        matrix("c", outputs, memory)?,
        matrix("d", outputs, inputs)?,
    )
    .map_err(|error| format!("{what}: {error}"))?;
    let state = ControllerState::new(&controller, initial)
        .map_err(|error| format!("{what}.initial: {error}"))?;
    add(builder, name, controller, state)
}

fn controlled_source(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    kind: Kind,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(definition, &["kind", "domain"], what)?;
    let domain = definition
        .str_field("domain")
        .and_then(Domain::parse)
        .ok_or_else(|| format!("{what}: domain is one of translation, rotation, electrical, hydraulic"))?;
    let source = if kind == Kind::ControlledEffort {
        ControlledSource::effort(domain.port_kind())
    } else {
        ControlledSource::flow(domain.port_kind())
    };
    add(builder, name, source, SourceAccount::default())
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
            .ok_or_else(|| format!("{what}: {key} is one of translation, rotation, electrical, hydraulic"))
    };
    let (a, b) = (side("a")?, side("b")?);
    let gyrator = Gyrator::new(a.port_kind(), b.port_kind(), number(definition, "ratio", what)?)
        .map_err(|e| Refusal::Input(format!("{what}: {e}")))?;
    add(builder, name, gyrator, ())
}

fn modulated_transformer(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(definition, &["kind", "a", "b"], what)?;
    let side = |key: &str| {
        definition
            .str_field(key)
            .and_then(Domain::parse)
            .ok_or_else(|| format!("{what}: {key} is one of translation, rotation, electrical, hydraulic"))
    };
    let (a, b) = (side("a")?, side("b")?);
    add(
        builder,
        name,
        ModulatedTransformer::new(a.port_kind(), b.port_kind()),
        (),
    )
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

/// One body of a rigid network: its mass properties, its starting state and,
/// for a body declared as a `solid`, the report of what its mass came from.
fn network_body(
    row: &J,
    scene: &Path,
    what: &str,
) -> Result<(MassProperties, RigidBodyState, Option<String>), Refusal> {
    fields(
        row,
        &[
            "mass_kg",
            "principal_inertia_kg_m2",
            "solid",
            "position_m",
            "orientation",
            "velocity_m_s",
            "angular_velocity_rad_s",
        ],
        what,
    )?;
    let refused = |e: fs_mbd::DynamicsError| Refusal::Input(format!("{what}: {e}"));
    let (mass_kg, principal, report) = match row.get("solid") {
        Some(_) if row.get("mass_kg").is_some() || row.get("principal_inertia_kg_m2").is_some() => {
            return Err(format!(
                "{what}: a solid's mass and inertia come from its shape and card; \
                 mass_kg and principal_inertia_kg_m2 are not given with it"
            )
            .into());
        }
        Some(declared) => {
            let body = solid::solid(declared, scene, &format!("{what}.solid"))?;
            (
                body.mass_kg,
                body.principal_inertia_kg_m2,
                Some(body.report),
            )
        }
        None => (
            positive(row, "mass_kg", what)?,
            triple(row, "principal_inertia_kg_m2", None, what)?,
            None,
        ),
    };
    let properties = MassProperties::new(mass_kg, Vec3::ZERO, vec3(principal)).map_err(refused)?;
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
    Ok((properties, state, report))
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

    /// Bind the complete original calibration, including its refusal domain.
    /// Sampling the law would lose identity at the reference temperature,
    /// where every admitted calibration has the same exact factor of one.
    fn identity(&self) -> [u8; 32] {
        let mut identity =
            fs_blake3::DomainHasher::new("org.frankensim.fs-cli.system-scene.temperature-shift.v1");
        match self.0.model() {
            TemperatureShiftModel::Arrhenius { activation_energy_j_mol } => {
                identity.update(&[0]);
                identity.update(&activation_energy_j_mol.to_bits().to_le_bytes());
            }
            TemperatureShiftModel::WilliamsLandelFerry { c1, c2_k } => {
                identity.update(&[1]);
                identity.update(&c1.to_bits().to_le_bytes());
                identity.update(&c2_k.to_bits().to_le_bytes());
            }
        }
        let [lowest, highest] = self.0.validity_k();
        for value in [self.0.reference_temperature_k(), lowest, highest] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        *identity.finalize().as_bytes()
    }
}

impl ShiftLaw for SceneShift {
    fn factor(&self, _connector: usize, temperature_k: f64) -> Result<f64, String> {
        self.at(temperature_k)
    }

    fn checkpoint_identity(&self) -> Option<[u8; 32]> {
        Some(self.identity())
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

    fn checkpoint_identity(&self) -> Option<[u8; 32]> {
        Some(self.identity())
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
    /// Each body declared as a solid: its index and report object.
    solids: Vec<(usize, String)>,
}

fn network_base<'a>(
    definition: &'a J,
    scene: &Path,
    what: &str,
) -> Result<NetworkBase<'a>, Refusal> {
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
    let mut solids = Vec::new();
    for (index, row) in rows.iter().enumerate() {
        let (body, state, solid) = network_body(row, scene, &format!("{what}.bodies[{index}]"))?;
        properties.push(body);
        states.push(state);
        if let Some(report) = solid {
            solids.push((index, format!("{{\"body\":{index},{report}}}")));
        }
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
        solids,
    })
}

/// A network element lowered into the builder: its node and ports, its
/// bodies, and for each connector the number of its Maxwell branches.
struct NetworkParts {
    node: NodeId,
    ports: Vec<PortSpec>,
    bodies: Vec<MassProperties>,
    branches: Vec<usize>,
    solids: Vec<(usize, String)>,
}

/// A `rigid_network`: rigid bodies joined by viscous connectors, lowered to
/// the fs-mbd network Component.
fn rigid_network(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    scene: &Path,
    what: &str,
) -> Result<NetworkParts, Refusal> {
    let base = network_base(definition, scene, what)?;
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
        solids: base.solids,
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
    scene: &Path,
    what: &str,
) -> Result<NetworkParts, Refusal> {
    let base = network_base(definition, scene, what)?;
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
        solids: base.solids,
    })
}

fn speed_controller(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "sample_period_s",
            "proportional_gain_nm_s_per_rad",
            "integral_gain_nm_per_rad",
            "anti_windup_time_s",
            "positive_torque_limit_nm",
            "negative_torque_limit_nm",
            "motoring_power_w",
            "regeneration_power_w",
            "integral_torque_nm",
        ],
        what,
    )?;
    let parameters = SpeedControllerParameters {
        sample_period_s: number(definition, "sample_period_s", what)?,
        proportional_gain_nm_s_per_rad: number(definition, "proportional_gain_nm_s_per_rad", what)?,
        integral_gain_nm_per_rad: number(definition, "integral_gain_nm_per_rad", what)?,
        anti_windup_time_s: number(definition, "anti_windup_time_s", what)?,
    };
    let limits = ShaftPowerLimits {
        positive_torque_limit_nm: number(definition, "positive_torque_limit_nm", what)?,
        negative_torque_limit_nm: number(definition, "negative_torque_limit_nm", what)?,
        motoring_power_w: number(definition, "motoring_power_w", what)?,
        regeneration_power_w: number(definition, "regeneration_power_w", what)?,
    };
    let invalid = |error: fs_control::ControlError| Refusal::Input(format!("{what}: {error}"));
    let controller = SampledSpeedController::new(parameters, limits).map_err(invalid)?;
    let initial = SampledSpeedState::new(number_or(definition, "integral_torque_nm", 0.0, what)?)
        .map_err(invalid)?;
    add(builder, name, controller, initial)
}

fn element(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    materials: &materials::Materials,
    prepared: &mut conduction::Prepared,
    scene: &Path,
    existing: &[Element],
) -> Result<Element, Refusal> {
    let what = format!("elements.{name}");
    let kind_name = definition
        .str_field("kind")
        .ok_or_else(|| format!("{what}: missing kind"))?;
    let kind = Kind::parse(kind_name).ok_or_else(|| format!("{what}: unknown kind {kind_name:?}"))?;
    let invalid = |e: fs_system::elements::ElementError| Refusal::Input(format!("{what}: {e}"));
    let mut temperature_k = None;
    let mut port_temperatures_k = Vec::new();
    let mut node_capacity = None;
    let mut charge_voltage_curve = None;
    let mut flux_current_curve = None;
    let mut bodies = Vec::new();
    let mut branches = Vec::new();
    let mut gas = None;
    let mut material_sources = Vec::new();
    let mut solids = Vec::new();
    let mut coordinates = 0;
    let mut radius_m = None;
    let mut contact = None;
    let mut enthalpy = None;
    let mut modal = None;
    let mut ohmic = None;
    let mut laminar = None;
    let (node, ports) = match kind {
        Kind::Inertance(domain) => {
            let names = domain.names();
            fields(
                definition,
                &["kind", names.inertance, names.displacement, names.flow, "publish", "publish_state"],
                &what,
            )?;
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
            let flag = |key: &str| -> Result<Option<bool>, String> {
                match definition.get(key) {
                    None => Ok(None),
                    Some(J::Bool(value)) => Ok(Some(*value)),
                    _ => Err(format!("{what}.{key}: expected a boolean")),
                }
            };
            let publish = flag("publish")?;
            let publish_state = flag("publish_state")?;
            if matches!((publish, publish_state), (Some(a), Some(b)) if a != b) {
                return Err(format!("{what}: publish and publish_state must agree when both are supplied").into());
            }
            let node = if publish.or(publish_state).unwrap_or(false) { node.publishing() } else { node };
            add(builder, name, node, state)?
        }
        Kind::Compliance(domain) => compliance(builder, name, definition, domain, &what)?,
        Kind::Resistance(Domain::Electrical) if definition.get("material").is_some() => {
            let (resistor, material) =
                ohmic_resistor(definition, materials, Terminals::Current, &what)?;
            ohmic = Some(resistor.node().clone());
            material_sources.push(material);
            add(builder, name, resistor, ())?
        }
        Kind::Resistance(domain) => {
            let names = domain.names();
            fields(definition, &["kind", names.resistance], &what)?;
            let resistance = number(definition, names.resistance, &what)?;
            let damper = Damper::in_domain(domain.port_kind(), resistance)
                .map_err(|_| format!("{what}: {} must be positive", names.resistance))?;
            add(builder, name, damper, ())?
        }
        Kind::ResistanceMatrix => resistance_matrix(builder, name, definition, &what)?,
        Kind::AffineSignal => affine_signal(builder, name, definition, &what)?,
        Kind::SignalSource => signal_source(builder, name, definition, &what)?,
        Kind::SampledController => sampled_controller(builder, name, definition, &what)?,
        Kind::SpeedController => speed_controller(builder, name, definition, &what)?,
        Kind::ControlledEffort | Kind::ControlledFlow => {
            controlled_source(builder, name, definition, kind, &what)?
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
        Kind::VelocitySource => {
            fields(definition, &["kind", "velocity_m_s"], &what)?;
            let wave = definition
                .get("velocity_m_s")
                .ok_or_else(|| format!("{what}: missing velocity_m_s"))?;
            let source = FlowSource::new(waveform(wave, &format!("{what}.velocity_m_s"))?)
                .publishing();
            add(builder, name, source, SourceAccount::default())?
        }
        Kind::Slider => slider(builder, name, definition, &what)?,
        Kind::FrictionInterface => friction::interface(builder, name, definition, &what)?,
        Kind::SlidingInterface => friction::sliding(builder, name, definition, &what)?,
        Kind::HertzPatch => {
            let (node, ports, radius) = friction::hertz_patch(builder, name, definition, &what)?;
            radius_m = Some(radius);
            (node, ports)
        }
        Kind::ModalStructure => modal::structure(builder, name, definition, &what, &mut modal)?,
        Kind::ContactInterface => contact::interface(builder, name, definition, &what)?,
        Kind::NodeCapacitor(domain) => {
            let names = domain.names();
            fields(
                definition,
                &["kind", names.stiffness, names.extension],
                &what,
            )?;
            let capacity = positive(definition, names.stiffness, &what)?;
            let content = number_or(definition, names.extension, 0.0, &what)?;
            if !content.is_finite() {
                return Err(format!("{what}: {} must be finite", names.extension).into());
            }
            let storage = Storage::in_domain(domain.port_kind(), capacity).map_err(invalid)?;
            node_capacity = Some(capacity);
            add(builder, name, storage, content)?
        }
        Kind::IdealDiode => {
            fields(definition, &["kind"], &what)?;
            let diode = IdealDiode::new();
            let ports = diode.ports().to_vec();
            let node = builder
                .add_on_line(name, diode, ())
                .map_err(|error| Refusal::Admission(error.to_string()))?;
            (node, ports)
        }
        Kind::IdealSwitch => {
            fields(definition, &["kind"], &what)?;
            let switch = IdealSwitch::new();
            let ports = switch.ports().to_vec();
            let node = builder
                .add_on_line(name, switch, ())
                .map_err(|error| Refusal::Admission(error.to_string()))?;
            (node, ports)
        }
        Kind::ResistorBranch(Domain::Electrical) if definition.get("material").is_some() => {
            let (resistor, material) =
                ohmic_resistor(definition, materials, Terminals::Voltage, &what)?;
            ohmic = Some(resistor.node().clone());
            material_sources.push(material);
            add(builder, name, resistor, ())?
        }
        Kind::ResistorBranch(Domain::Hydraulic) if definition.get("material").is_some() => {
            let (pipe, material) = laminar_pipe(definition, materials, &what)?;
            laminar = Some(pipe.pipe().clone());
            material_sources.push(material);
            add(builder, name, pipe, ())?
        }
        Kind::ResistorBranch(domain) => {
            let names = domain.names();
            fields(definition, &["kind", names.resistance], &what)?;
            let resistance = positive(definition, names.resistance, &what)?;
            let branch =
                Admittance::in_domain(domain.port_kind(), 1.0 / resistance).map_err(invalid)?;
            add(builder, name, branch, ())?
        }
        Kind::CurrentSource(domain) => {
            let names = domain.names();
            fields(definition, &["kind", names.flow], &what)?;
            let wave = definition
                .get(names.flow)
                .ok_or_else(|| format!("{what}: missing {}", names.flow))?;
            let source = FlowSource::in_domain(
                domain.port_kind(),
                waveform(wave, &format!("{what}.{}", names.flow))?,
            );
            add(builder, name, source, SourceAccount::default())?
        }
        Kind::CarreauPipe => {
            let pipe = hydraulic::carreau_pipe(definition, &what)?;
            add(builder, name, pipe, ())?
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
        Kind::ModulatedTransformer => modulated_transformer(builder, name, definition, &what)?,
        Kind::DiracJunction => add(builder, name, dirac::junction(definition, name, &what)?, ())?,
        Kind::PortHamiltonian => {
            let (owner, state) = phs::port_hamiltonian(definition, &what)?;
            coordinates = state.len();
            add(builder, name, owner, state)?
        }
        Kind::RigidNetwork | Kind::ViscoelasticNetwork => {
            let parts = if kind == Kind::RigidNetwork {
                rigid_network(builder, name, definition, scene, &what)?
            } else {
                viscoelastic_network(builder, name, definition, scene, &what)?
            };
            bodies = parts.bodies;
            branches = parts.branches;
            solids = parts.solids;
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
        Kind::EnthalpyBody => {
            let (body, state) = enthalpy_body(definition, materials, &what)?;
            if let Some(name) = definition.str_field("material") {
                material_sources.push(name.to_owned());
            }
            temperature_k = Some(state.temperature_k());
            add(builder, name, body, state)?
        }
        Kind::Conductance => {
            fields(definition, &["kind", "conductance_w_k"], &what)?;
            let wall =
                Conductance::new(number(definition, "conductance_w_k", &what)?).map_err(invalid)?;
            add(builder, name, wall, ())?
        }
        Kind::LayerConductor => {
            let conductor = layer_conductor(definition, materials, &what)?;
            if let Some(names) = definition.get("materials").and_then(J::as_array) {
                material_sources.extend(names.iter().filter_map(J::as_str).map(str::to_owned));
            }
            add(builder, name, conductor, ())?
        }
        Kind::HeatNode => {
            fields(definition, &["kind", "conductances_w_k", "temperature_k"], &what)?;
            let rows = definition
                .get("conductances_w_k")
                .and_then(J::as_array)
                .filter(|rows| (1..=MAX_NODE_BODIES).contains(&rows.len()))
                .ok_or_else(|| {
                    format!(
                        "{what}: conductances_w_k is an array of 1 to {MAX_NODE_BODIES} \
                         conductances, one per body"
                    )
                })?;
            let conductances = rows
                .iter()
                .map(|row| {
                    row.as_f64()
                        .ok_or_else(|| format!("{what}: conductances_w_k must be numbers"))
                })
                .collect::<Result<Vec<f64>, String>>()?;
            // The node stores nothing. A temperature given here is only
            // where the solve starts for what is joined to its source.
            if definition.get("temperature_k").is_some() {
                temperature_k = Some(positive(definition, "temperature_k", &what)?);
            }
            add(builder, name, HeatNode::new(&conductances).map_err(invalid)?, ())?
        }
        Kind::ThermalTimeline => {
            let parts = timeline::timeline(builder, name, definition, materials, &what)?;
            port_temperatures_k = parts.port_temperatures_k;
            material_sources = parts.material_sources;
            (parts.node, parts.ports)
        }
        Kind::ConductionBody => {
            let added = conduction::body(builder, name, definition, materials, prepared, &what)?;
            temperature_k = Some(added.temperature_k);
            material_sources.extend(added.material);
            enthalpy = added.enthalpy;
            (added.node, added.ports)
        }
        Kind::TraceContact => {
            let (node, ports, declared) = conduction::wall(builder, name, prepared, &what)?;
            contact = Some(declared);
            (node, ports)
        }
        Kind::GasVolume => {
            let (volume, contents) = gas_volume(definition, &what)?;
            temperature_k = Some(volume.temperature_k(&contents));
            gas = Some(volume.clone());
            add(builder, name, volume, contents)?
        }
        Kind::Restriction | Kind::ModulatedRestriction => {
            let law = restriction(definition, &what, builder, existing)?;
            let law = if kind == Kind::ModulatedRestriction {
                law.modulated()
            } else {
                law
            };
            add(builder, name, law, ())?
        }
        Kind::Heater => {
            fields(definition, &["kind", "power_w"], &what)?;
            let wave = definition
                .get("power_w")
                .ok_or_else(|| format!("{what}: missing power_w"))?;
            let heater = Heater::new(waveform(wave, &format!("{what}.power_w"))?).map_err(invalid)?;
            add(builder, name, heater, SourceAccount::default())?
        }
        Kind::GasCavity => {
            let (added, start_k) = gas::cavity(builder, name, definition, scene, &what)?;
            temperature_k = Some(start_k);
            added
        }
    };
    Ok(Element {
        name: name.to_owned(),
        kind,
        node,
        ports,
        point_sampled_command: false,
        temperature_k,
        port_temperatures_k,
        node_capacity,
        charge_voltage_curve,
        flux_current_curve,
        bodies,
        branches,
        gas,
        material_sources,
        solids,
        coordinates,
        radius_m,
        contact,
        enthalpy,
        modal,
        ohmic,
        laminar,
    })
}

/// A resistor whose resistance is a scene material's `rho(T) L / A`, which
/// fs-material's Ohmic node reads again at the temperature its heat port
/// receives, and the material's name.
fn ohmic_resistor(
    definition: &J,
    materials: &materials::Materials,
    terminals: Terminals,
    what: &str,
) -> Result<(OhmicResistor, String), Refusal> {
    fields(
        definition,
        &["kind", "material", "length_m", "area_m2"],
        what,
    )?;
    let name = definition
        .str_field("material")
        .ok_or_else(|| format!("{what}.material: expected a material name"))?;
    let material = materials
        .get(name)
        .ok_or_else(|| format!("{what}.material: unknown material {name:?}"))?;
    let length_m = positive(definition, "length_m", what)?;
    let area_m2 = positive(definition, "area_m2", what)?;
    let node = material
        .ohmic_node(length_m, area_m2)
        .map_err(|error| format!("{what}.material: {error}"))?;
    let resistor =
        OhmicResistor::new(node, terminals).map_err(|error| format!("{what}.material: {error}"))?;
    Ok((resistor, name.to_owned()))
}

/// A hydraulic resistor that is a laminar circular pipe of a scene
/// material's liquid, whose viscosity fs-material's pipe reads again at the
/// temperature its heat port receives, and the material's name.
fn laminar_pipe(
    definition: &J,
    materials: &materials::Materials,
    what: &str,
) -> Result<(LaminarPipeComponent, String), Refusal> {
    fields(
        definition,
        &["kind", "material", "radius_m", "length_m"],
        what,
    )?;
    let name = definition
        .str_field("material")
        .ok_or_else(|| format!("{what}.material: expected a material name"))?;
    let material = materials
        .get(name)
        .ok_or_else(|| format!("{what}.material: unknown material {name:?}"))?;
    let radius_m = positive(definition, "radius_m", what)?;
    let length_m = positive(definition, "length_m", what)?;
    let pipe = material
        .laminar_pipe(radius_m, length_m)
        .map_err(|error| format!("{what}.material: {error}"))?;
    Ok((LaminarPipeComponent::new(pipe), name.to_owned()))
}

/// A card-backed pipe's temperature, as its heat port read it in the last
/// step, and the card's viscosity there. A hydraulic resistor of constant
/// resistance has neither, and its report is unchanged.
fn laminar_values(system: &System, element: &Element) -> Vec<(&'static str, f64)> {
    let Some(pipe) = &element.laminar else {
        return Vec::new();
    };
    let temperature = match system.port(element.node, "heat") {
        Some(port) if system.steps() > 0 => system.last_input(port)[0],
        _ => f64::NAN,
    };
    let viscosity = pipe
        .viscosity_and_slope(temperature)
        .map_or(f64::NAN, |(viscosity, _)| viscosity);
    vec![
        ("temperature_k", temperature),
        ("viscosity_pa_s", viscosity),
    ]
}

/// A card-backed resistor's temperature, as its heat port read it in the
/// last step, and the card's resistance there. A resistor of constant
/// resistance has neither, and its report is unchanged.
fn ohmic_values(system: &System, element: &Element) -> Vec<(&'static str, f64)> {
    let Some(node) = &element.ohmic else {
        return Vec::new();
    };
    let temperature = match system.port(element.node, "heat") {
        Some(port) if system.steps() > 0 => system.last_input(port)[0],
        _ => f64::NAN,
    };
    let resistance = node
        .resolve_at_temperature_k(temperature)
        .map_or(f64::NAN, |conductor| conductor.resistance_ohm());
    vec![
        ("temperature_k", temperature),
        ("resistance_ohm", resistance),
    ]
}

fn slider(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "normal_n",
            "static_coefficient",
            "kinetic_coefficient",
            "slip_m_s",
            "prescribed_ends",
        ],
        what,
    )?;
    let mut law = Slider::new(
        number(definition, "normal_n", what)?,
        number(definition, "static_coefficient", what)?,
        number(definition, "kinetic_coefficient", what)?,
    )
    .map_err(|error| format!("{what}: {error}"))?;
    if let Some(ends) = definition.get("prescribed_ends") {
        let ends = ends
            .as_array()
            .filter(|ends| ends.len() <= 2)
            .ok_or_else(|| format!("{what}.prescribed_ends: an array containing a and/or b"))?;
        let mut seen = [false; 2];
        for end in ends {
            let (index, end) = match end.as_str() {
                Some("a") => (0, SliderEnd::A),
                Some("b") => (1, SliderEnd::B),
                _ => return Err(format!("{what}.prescribed_ends: each end must be a or b").into()),
            };
            if seen[index] {
                return Err(format!("{what}.prescribed_ends: an end appears twice").into());
            }
            seen[index] = true;
            law = law.told(end);
        }
    }
    let slip = number_or(definition, "slip_m_s", 0.0, what)?;
    if !slip.is_finite() {
        return Err(format!("{what}.slip_m_s must be finite").into());
    }
    // Scene velocities use the mover's sign. Each slider end receives
    // minus its mover's flow, so the owner's slip has the opposite sign.
    let state = if slip == 0.0 {
        Slip::stuck()
    } else {
        Slip::sliding(-slip)
    };
    let ports = law.ports().to_vec();
    let node = builder
        .add_on_line(name, law, state)
        .map_err(|error| Refusal::Admission(error.to_string()))?;
    Ok((node, ports))
}

/// A `restriction` between gas volumes.
fn restriction(
    definition: &J,
    what: &str,
    builder: &SystemBuilder,
    elements: &[Element],
) -> Result<Restriction, Refusal> {
    fields(
        definition,
        &["kind", "conductance_kg_s_pa", "carried", "stream_basis"],
        what,
    )?;
    let conductance = number(definition, "conductance_kg_s_pa", what)?;
    if definition.get("stream_basis").is_some() {
        if definition.get("carried").is_some() {
            return Err(format!(
                "{what}: stream_basis determines the carried coordinates; omit carried"
            )
            .into());
        }
        let reference = definition
            .str_field("stream_basis")
            .ok_or_else(|| format!("{what}.stream_basis: expected a port name, cavity.stream"))?;
        let target = resolve(builder, elements, reference)?;
        let owner = &elements[target.element];
        let spec = owner
            .ports
            .iter()
            .find(|port| port.carrier() == Carrier::Stream)
            .filter(|port| {
                owner.kind == Kind::GasCavity
                    && target.carrier == Carrier::Stream
                    && port.stream_basis().is_some()
            })
            .ok_or_else(|| {
                format!(
                    "{what}.stream_basis: {reference:?} must name the stream port of a \
                     gas_cavity with stream_mass_tolerance_kg"
                )
            })?;
        return Restriction::for_stream(conductance, spec)
            .map_err(|error| Refusal::Input(format!("{what}: {error}")));
    }
    let carried = number_or(definition, "carried", 0.0, what)?;
    if !(0.0..=MAX_CARRIED as f64).contains(&carried) || carried.fract() != 0.0 {
        return Err(format!(
            "{what}: carried is the number of carried quantities, 0 to {MAX_CARRIED}"
        )
        .into());
    }
    Restriction::new(
        conductance,
        carried as usize,
    )
    .map_err(|e| Refusal::Input(format!("{what}: {e}")))
}

/// A `gas_volume`: the gas, and what it holds at the start.
fn gas_volume(definition: &J, what: &str) -> Result<(GasVolume, GasContents), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "gas_constant_j_kg_k",
            "heat_capacity_ratio",
            "volume_m3",
            "pressure_pa",
            "temperature_k",
            "carried",
            "wall",
        ],
        what,
    )?;
    let invalid = |e: fs_system::elements::ElementError| Refusal::Input(format!("{what}: {e}"));
    let fractions = match definition.get("carried") {
        None => Vec::new(),
        Some(list) => list
            .as_array()
            .filter(|rows| rows.len() <= MAX_CARRIED)
            .and_then(|rows| rows.iter().map(J::as_f64).collect::<Option<Vec<f64>>>())
            .ok_or_else(|| {
                format!(
                    "{what}: carried is an array of at most {MAX_CARRIED} mass fractions, \
                     one per carried quantity"
                )
            })?,
    };
    let wall = match definition.get("wall") {
        None => false,
        Some(J::Bool(wall)) => *wall,
        Some(_) => return Err(format!("{what}: wall is true or false").into()),
    };
    let volume = GasVolume::new(
        number(definition, "gas_constant_j_kg_k", what)?,
        number(definition, "heat_capacity_ratio", what)?,
        fractions.len(),
    )
    .map_err(invalid)?;
    let volume = if wall { volume.with_wall() } else { volume };
    let contents = volume
        .contents(
            number(definition, "pressure_pa", what)?,
            number(definition, "temperature_k", what)?,
            number(definition, "volume_m3", what)?,
            &fractions,
        )
        .map_err(invalid)?;
    Ok((volume, contents))
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
    /// The temperature the port starts at, if its element states one.
    temperature_k: Option<f64>,
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
    let temperature_k = elements[element]
        .port_temperatures_k
        .iter()
        .find(|(id, _)| id == port)
        .map(|&(_, kelvin)| kelvin)
        .or(elements[element].temperature_k);
    let port = builder
        .port(elements[element].node, port)
        .map_err(|e| Refusal::Admission(e.to_string()))?;
    Ok(Resolved {
        element,
        port,
        carrier: spec.carrier(),
        role: spec.role(),
        width: spec.width(),
        temperature_k,
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

/// Starting values for inputs inside loops, in each port's coefficient order.
type Guesses = Vec<(PortRef, Vec<f64>)>;

/// A trace contact was computed for two traces. Its port `a` is joined to
/// the first and its port `b` to the second, each alone: any other trace,
/// even one of the same width, has other points.
fn contact_wiring(
    builder: &SystemBuilder,
    elements: &[Element],
    ports: &[(Resolved, Option<f64>)],
    what: &str,
) -> Result<(), Refusal> {
    for (port, _) in ports {
        let element = &elements[port.element];
        let Some(declared) = &element.contact else {
            continue;
        };
        let end = usize::from(builder.port(element.node, "a").ok() != Some(port.port));
        let trace = resolve(builder, elements, &declared.between[end])?.port;
        if ports.len() != 2 || !ports.iter().any(|(other, _)| other.port == trace) {
            return Err(format!(
                "{what}: {}.{} is the contact's side at {}; join it to that trace alone",
                element.name,
                ["a", "b"][end],
                declared.between[end]
            )
            .into());
        }
    }
    Ok(())
}

fn join(
    builder: &mut SystemBuilder,
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
            // A gas volume is the node its streams meet at, as it is for
            // its heat.
            matches!(
                (port.carrier, port.role),
                (Carrier::Power(_), PortRole::EffortIn)
                    | (Carrier::Heat | Carrier::Stream, PortRole::FlowIn)
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
            "{what}: a join has exactly one node port (a mass, inertia, inductor, heat \
             capacity, a heat node's source, or a gas volume's stream or heat) and the \
             ports that meet at it; this one has {nodes}"
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
    if let (Carrier::Heat, Some(kelvin)) = (ports[setter].carrier, ports[setter].temperature_k) {
        guesses.extend(
            members
                .iter()
                .map(|&member| (member, vec![kelvin; ports[setter].width])),
        );
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

/// A rectangular virtual-work map, with the adjoint owned by the kernel.
fn power_map(
    builder: &mut SystemBuilder,
    elements: &[Element],
    connection: &J,
    what: &str,
) -> Result<(), Refusal> {
    fields(connection, &["power_map", "matrix"], what)?;
    let Some([a, b]) = connection.get("power_map").and_then(J::as_array) else {
        return Err(format!(
            "{what}: power_map is [effort_in_port, flow_in_port], or [[members], flow_in_port]"
        )
        .into());
    };
    // Several effort-in ports onto one flow-in port: each member has its own
    // matrix, and the flow-in port receives minus the sum of their maps.
    if let Some(members) = a.as_array() {
        if connection.get("matrix").is_some() {
            return Err(format!(
                "{what}: a power_map of members has a matrix per member, not one for the map"
            )
            .into());
        }
        let b = resolve(builder, elements, &port_text(b, what)?)?;
        let mut maps = Vec::with_capacity(members.len());
        for (index, member) in members.iter().enumerate() {
            let inner = format!("{what}.power_map[0][{index}]");
            fields(member, &["port", "matrix"], &inner)?;
            let port = member
                .get("port")
                .ok_or_else(|| format!("{inner}: missing port"))?;
            let a = resolve(builder, elements, &port_text(port, &inner)?)?;
            let map = coefficient_matrix(member.get("matrix"), b.width, a.width, &inner)?;
            maps.push((a.port, map));
        }
        return builder
            .power_map_members(maps, b.port)
            .map_err(|error| Refusal::Admission(format!("{what}: {error}")));
    }
    let a = resolve(builder, elements, &port_text(a, what)?)?;
    let b = resolve(builder, elements, &port_text(b, what)?)?;
    let map = coefficient_matrix(connection.get("matrix"), b.width, a.width, what)?;
    builder
        .power_map(a.port, b.port, map)
        .map_err(|error| Refusal::Admission(format!("{what}: {error}")))
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
            if let Some(kelvin) = sink.temperature_k {
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
        guesses.push((source.port, vec![kelvin]));
    }
    Ok(())
}

/// `{"signal": output, "to": [inputs], "matrix": [[...]]}`: copy a signal,
/// or apply the same finite rectangular map to it at each destination.
/// Rows follow the receiving port and columns follow the source port.
fn signal(
    builder: &mut SystemBuilder,
    elements: &mut [Element],
    connection: &J,
    guesses: &mut Guesses,
    what: &str,
) -> Result<(), Refusal> {
    fields(connection, &["signal", "to", "matrix"], what)?;
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
        let map = connection
            .get("matrix")
            .map(|value| coefficient_matrix(Some(value), target.width, source.width, what))
            .transpose()?;
        // Keep temperature-informed starts consistent with the connection.
        // A rectangular map may give each receiving coefficient a different
        // value, so the old uniform broadcast is only the copy case.
        let initial = source.temperature_k.map(|kelvin| {
            if let Some(map) = &map {
                map.chunks_exact(source.width)
                    .map(|row| row.iter().map(|&coefficient| coefficient * kelvin).sum())
                    .collect::<Vec<f64>>()
            } else {
                vec![kelvin; target.width]
            }
        });
        if initial
            .as_ref()
            .is_some_and(|values| values.iter().any(|value| !value.is_finite()))
        {
            return Err(format!("{what}: initial mapped signal must be finite").into());
        }
        let wired = if let Some(map) = map {
            builder.signal_map(source.port, target.port, map)
        } else {
            builder.signal(source.port, target.port)
        };
        wired.map_err(|e| Refusal::Admission(format!("{what}: {e}")))?;
        if matches!(
            elements[target.element].kind,
            Kind::ControlledEffort | Kind::ControlledFlow
        ) {
            elements[target.element].point_sampled_command = true;
        }
        if let Some(initial) = initial {
            guesses.push((target.port, initial));
        }
    }
    Ok(())
}

fn prescribe(
    builder: &mut SystemBuilder,
    elements: &mut [Element],
    connection: &J,
    what: &str,
) -> Result<(), Refusal> {
    let target = resolve(
        builder,
        elements,
        &port_text(connection.get("prescribe").unwrap_or(&J::Null), what)?,
    )?;
    if target.carrier == Carrier::Signal
        && target.role == PortRole::SignalIn
        && elements[target.element].kind == Kind::Slider
    {
        fields(connection, &["prescribe", "velocity_m_s"], what)?;
        let wave = connection
            .get("velocity_m_s")
            .ok_or_else(|| format!("{what}: a slider span takes velocity_m_s"))?;
        let wave = waveform(wave, &format!("{what}.velocity_m_s"))?;
        let events = wave.event_times_s().collect();
        return builder
            .scheduled_boundary(
                target.port,
                Box::new(move |start, dt, values| {
                    values[0] = wave.value_at(start);
                    values[1] = wave.value_at(start + dt);
                }),
                events,
            )
            .map_err(|error| Refusal::Admission(format!("{what}: {error}")));
    }
    let domain = match target.carrier {
        Carrier::Power(kind) => Domain::of(kind),
        _ => None,
    };
    let element = &elements[target.element];
    let controlled = controlled_quantity(element);
    let value_command = controlled.is_some() && connection.get("value").is_some();
    let point_sample = target.carrier == Carrier::Signal
        && (value_command
            || matches!(
                element.kind,
                Kind::AffineSignal
                    | Kind::SampledController
                    | Kind::SpeedController
                    | Kind::IdealSwitch
                    | Kind::ModulatedRestriction
            ));
    let signal_quantity = if element.kind == Kind::ModulatedRestriction {
        Some("opening")
    } else if element.kind == Kind::SpeedController {
        Some(
            if builder.port(element.node, "reference").ok() == Some(target.port) {
                "relative_speed_rad_s"
            } else if builder.port(element.node, "feedforward").ok() == Some(target.port) {
                "torque_n_m"
            } else {
                "angular_velocity_rad_s"
            },
        )
    } else if value_command
        || matches!(
            element.kind,
            Kind::AffineSignal | Kind::SampledController | Kind::IdealSwitch
        )
    {
        Some("value")
    } else {
        controlled
    };
    let quantity = match (target.carrier, target.role, domain, signal_quantity) {
        (Carrier::Signal, PortRole::SignalIn, _, Some(quantity)) => quantity,
        (Carrier::Signal, PortRole::SignalIn, _, _)
            if element.kind == Kind::ModulatedTransformer =>
        {
            "ratio"
        }
        // Control signals, spans, openings and ratios were handled above; remaining
        // signals are temperatures at which a network's shift is read.
        (Carrier::Heat, PortRole::EffortIn, _, _) | (Carrier::Signal, PortRole::SignalIn, _, _) => {
            "temperature_k"
        }
        (Carrier::Heat, PortRole::FlowIn, _, _) => "heat_rate_w",
        (Carrier::Power(_), PortRole::EffortIn, Some(domain), _) => domain.names().effort,
        (Carrier::Power(_), PortRole::FlowIn, Some(domain), _) => domain.names().flow,
        // What crosses into a gas volume: mass rate [kg/s], enthalpy rate
        // [W] and the rate of each carried quantity [kg/s].
        (Carrier::Stream, PortRole::FlowIn, _, _) => "stream",
        // What a restriction's end is held at: pressure [Pa], specific
        // enthalpy [J/kg] and the mass fraction of each carried quantity.
        (Carrier::Stream, PortRole::EffortIn, _, _) => "contents",
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
    let mut events: Vec<f64> = waves.iter().flat_map(|wave| wave.event_times_s()).collect();
    events.sort_by(f64::total_cmp);
    // A vector/stream port uses the union of all coordinates' knots. Treat
    // signed zeros as the same physical time as the builder's ordering does.
    events.dedup_by(|a, b| (*a).partial_cmp(&*b) == Some(std::cmp::Ordering::Equal));
    // A series element end carried by a node receives minus the node's
    // flow. Node capacitors and ideal electrical devices instead use current
    // positive into their own p port, so their prescriptions keep that sign;
    // a gas cavity likewise keeps its compression flow's.
    let prescribed_flow = matches!(
        (target.carrier, target.role),
        (Carrier::Power(_), PortRole::FlowIn)
    ) && !matches!(
        elements[target.element].kind,
        Kind::NodeCapacitor(_)
            | Kind::NonlinearCapacitor
            | Kind::IdealDiode
            | Kind::IdealSwitch
            | Kind::GasCavity
    );
    let sign = if prescribed_flow { -1.0 } else { 1.0 };
    builder
        .scheduled_boundary(
            target.port,
            Box::new(move |t0, dt, values| {
                for (value, wave) in values.iter_mut().zip(&waves) {
                    *value = sign
                        * if point_sample {
                            wave.value_at(t0)
                        } else {
                            wave.mean(t0, dt)
                        };
                }
            }),
            events,
        )
        .map_err(|e| Refusal::Admission(format!("{what}: {e}")))?;
    if value_command && target.carrier == Carrier::Signal {
        elements[target.element].point_sampled_command = true;
    }
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
            Kind::Compliance(domain) | Kind::NodeCapacitor(domain) => {
                vec![(domain.names().extension, Coordinate::Scalar)]
            }
            Kind::NonlinearCapacitor => vec![("charge_c", Coordinate::Scalar)],
            Kind::NonlinearInductor => vec![("flux_wb", Coordinate::Scalar)],
            Kind::HeatCapacity => vec![("temperature_k", Coordinate::Scalar)],
            Kind::EnthalpyBody => vec![("specific_enthalpy_j_kg", Coordinate::Enthalpy)],
            Kind::Source(_)
            | Kind::VelocitySource
            | Kind::CurrentSource(_)
            | Kind::ControlledEffort
            | Kind::ControlledFlow
            | Kind::Heater => {
                vec![("energy_in_j", Coordinate::SourceWork)]
            }
            Kind::Resistance(_)
            | Kind::AffineSignal
            | Kind::SignalSource
            | Kind::ResistanceMatrix
            | Kind::ResistorBranch(_)
            | Kind::CarreauPipe
            | Kind::Gyrator
            | Kind::ModulatedTransformer
            | Kind::DiracJunction
            | Kind::Conductance
            | Kind::TraceContact
            | Kind::LayerConductor
            | Kind::HeatNode
            | Kind::Restriction
            | Kind::ModulatedRestriction => Vec::new(),
            Kind::PortHamiltonian => {
                for index in 0..element.coordinates {
                    expected.push((
                        format!("{}.x{index}", element.name),
                        element.node,
                        Coordinate::State(index),
                    ));
                }
                Vec::new()
            }
            Kind::RigidNetwork | Kind::ViscoelasticNetwork => {
                for (body, &properties) in element.bodies.iter().enumerate() {
                    for (field, quantity) in adaptive_metric::BodyQuantity::ACCOUNTS {
                        expected.push((
                            format!("{}.body{body}.{field}", element.name),
                            element.node,
                            Coordinate::Network(adaptive_metric::NetworkCoordinate::Body {
                                body,
                                properties,
                                quantity,
                            }),
                        ));
                    }
                }
                for (connector, &count) in element.branches.iter().enumerate() {
                    for branch in 0..count {
                        expected.push((
                            format!(
                                "{}.connector{connector}.branch{branch}.amplitude_sqrt_j",
                                element.name
                            ),
                            element.node,
                            Coordinate::Network(adaptive_metric::NetworkCoordinate::Branch {
                                connector,
                                branch,
                            }),
                        ));
                    }
                }
                Vec::new()
            }
            Kind::GasVolume => {
                expected.extend(gas_volume_accounts(element));
                Vec::new()
            }
            Kind::GasCavity => {
                expected.extend(gas_cavity_accounts(element));
                Vec::new()
            }
            Kind::ConductionBody | Kind::ModalStructure => vec![whole_state_account(element)],
            Kind::Slider
            | Kind::FrictionInterface
            | Kind::SlidingInterface
            | Kind::HertzPatch
            | Kind::ContactInterface
            | Kind::SampledController
            | Kind::SpeedController
            | Kind::IdealDiode
            | Kind::IdealSwitch
            | Kind::ThermalTimeline => return Err(no_error_account(element)),
        };
        for (field, coordinate) in coordinates {
            expected.push((format!("{}.{field}", element.name), element.node, coordinate));
        }
    }
    declared_order_fits(order, elements)?;
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

fn load_text(text: &str, source_path: &Path) -> Result<Device, Refusal> {
    let root = J::parse(text).map_err(|e| e.to_string())?;
    fields(
        &root,
        &[
            "schema",
            "duration_s",
            "step_s",
            "sample_every_s",
            "solver",
            "adaptive",
            "materials",
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
    let materials = materials::load(&root, source_path)?;
    // A trace contact is computed from two bodies, so those are built first.
    let mut prepared = conduction::prepare(definitions, &materials, source_path)?;
    let mut elements = Vec::with_capacity(definitions.len());
    // Source-identified restrictions copy a cavity's native port declaration.
    // Build those after the storages, regardless of the scene's element order.
    for deferred in [false, true] {
        for (name, definition) in definitions {
            let needs_basis = matches!(
                definition.str_field("kind"),
                Some("restriction") | Some("modulated_restriction")
            ) && definition.get("stream_basis").is_some();
            if needs_basis != deferred {
                continue;
            }
            let added = element(
                &mut builder,
                name,
                definition,
                &materials,
                &mut prepared,
                source_path,
                &elements,
            )?;
            elements.push(added);
        }
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
            contact_wiring(&builder, &elements, &ports, &what)?;
            join(&mut builder, &ports, &mut guesses, common_effort, &what)?;
        } else if connection.get("transformer").is_some() {
            transformer(&mut builder, &elements, connection, &what)?;
        } else if connection.get("power_map").is_some() {
            power_map(&mut builder, &elements, connection, &what)?;
        } else if connection.get("route").is_some() {
            route(&mut builder, &elements, connection, &mut guesses, &what)?;
        } else if connection.get("signal").is_some() {
            signal(&mut builder, &mut elements, connection, &mut guesses, &what)?;
        } else if connection.get("prescribe").is_some() {
            prescribe(&mut builder, &mut elements, connection, &what)?;
        } else {
            return Err(format!(
                "{what}: expected join, common_effort, transformer, power_map, route, \
                 signal or prescribe"
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
    for (port, values) in guesses {
        system.guess(port, &values);
    }
    let adaptive = root.get("adaptive")
        .map(|definition| adaptive(definition, &elements, step_s)).transpose()?;
    let mut device = Device {
        system,
        elements,
        materials,
        connections: connections.len(),
        duration_s,
        step_s,
        sample_every_s,
        samples: Vec::new(),
        next_sample_s: sample_every_s.unwrap_or(f64::INFINITY),
        iterations: 0,
        shortenings: 0,
        adaptive,
    };
    // The initial sample belongs to admission, so a cancelled first advance
    // does not mutate the accepted report or duplicate its initial row.
    device.sample();
    Ok(device)
}

impl Device {
    #[allow(clippy::too_many_lines)] // one narrowing sentence per element kind that has one
    fn no_claim(&self) -> std::borrow::Cow<'static, str> {
        let base = if self.adaptive.is_some() {
            ADAPTIVE_NO_CLAIM
        } else if self
            .elements
            .iter()
            .any(|element| !element.bodies.is_empty())
        {
            NETWORK_NO_CLAIM
        } else {
            NO_CLAIM
        };
        let mut description = std::borrow::Cow::Borrowed(base);
        if self
            .elements
            .iter()
            .any(|element| matches!(element.kind, Kind::IdealDiode | Kind::IdealSwitch))
        {
            description = std::borrow::Cow::Owned(format!(
                "{description}; ideal electrical devices enforce their open/short relation on \
                 each representative step with binary held switch commands; no semiconductor \
                 losses, recovery, breakdown, arcing, impulses or interior zero-crossing model"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::Slider)
        {
            // Keep existing scene reports unchanged; only a scene with a
            // slider replaces the former blanket friction exclusion.
            let base = description.replace("contact, friction,", "geometric contact, impact,");
            description = std::borrow::Cow::Owned(format!(
                "{base}; dry friction uses constant static/kinetic coefficients and normal force, \
                 with first-order switching steps and explicitly declared initial slip"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::FrictionInterface)
        {
            let base = description.replace("contact, friction,", "geometric contact, impact,");
            description = std::borrow::Cow::Owned(format!(
                "{base}; friction interfaces step an fs-tribo dry law (authored coefficients or \
                 the cited rosin fit, independent of temperature) under a constant normal force, \
                 with every dissipated joule sent to their heat port; a prescribed end is held \
                 over each step, so an accelerating one advances in jumps; no contact geometry, \
                 normal dynamics, wear or third body"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::SlidingInterface)
        {
            let base = description.replace("contact, friction,", "geometric contact, impact,");
            description = std::borrow::Cow::Owned(format!(
                "{base}; sliding interfaces step the same fs-tribo dry laws at a constant \
                 prescribed slip rate and normal force, with every dissipated joule sent to their \
                 heat port and the prescribed motion's work kept as the interface's own account; \
                 no contact geometry, normal dynamics, wear or third body"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::HertzPatch)
        {
            let base = description.replace("contact, friction,", "geometric contact, impact,");
            description = std::borrow::Cow::Owned(format!(
                "{base}; Hertz patches step fs-tribo's joint slip and spin friction of one \
                 circular patch with constant static and kinetic coefficients, normal load and \
                 radius, through constant maps of their bodies' ports fixed at admission, with \
                 every dissipated joule sent to their heat port; a step in which the slip or the \
                 spin would reverse is refused, not shortened; no normal contact law, rolling, \
                 partial slip, wear or third body"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::ModalStructure)
        {
            description = std::borrow::Cow::Owned(format!(
                "{description}; modal structures step fs-nlmodal's damped modes exactly under \
                 the effort held over each step, from authored frequencies, damping ratios and \
                 mode-shape weights, or from the modes fs-nlmodal solves from an authored \
                 stiffness and mass under Rayleigh damping; no mode solve from geometry or \
                 material cards"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::ContactInterface)
        {
            // The scene has a contact, and a strike through it is an impact:
            // neither stays in the list of what the scene has none of.
            let base = description
                .replace("contact, friction,", "friction,")
                .replace("geometric contact, impact,", "geometric contact,");
            description = std::borrow::Cow::Owned(format!(
                "{base}; contact interfaces step fs-dcontact's one-point normal power law \
                 with an authored stiffness, exponent and gap, its reaction the secant of the \
                 law's potential over each step's move, settled on the load line by bracketed \
                 bisection within the scene's budgets; no contact loss, no events at touch, no \
                 contact geometry"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::ConductionBody)
        {
            // A scene with a trace contact does have contact traces.
            let contacts = self
                .elements
                .iter()
                .any(|element| element.kind == Kind::TraceContact);
            let without = if contacts {
                "mesh or radiation"
            } else {
                "mesh, radiation or contact trace"
            };
            // A body that names a material takes the material's k(T); one
            // that declares no heat capacity takes the card's Cp(T) and
            // density too, and otherwise the heat capacity stays the scene's
            // constant.
            let material_body = self.elements.iter().any(|element| {
                element.kind == Kind::ConductionBody && !element.material_sources.is_empty()
            });
            let card_capacity = self
                .elements
                .iter()
                .any(|element| element.enthalpy.is_some());
            let capacity = if card_capacity {
                "volumetric heat capacity or a material card's Cp(T) and density, whose specific \
                 enthalpy fs-conduction's enthalpy owner stores on the card's chart, with masses \
                 fixed by the card's density at the body's reference temperature and no thermal \
                 expansion"
            } else {
                "volumetric heat capacity"
            };
            let (conductivity, dependent) = if material_body {
                (
                    "one declared conductivity or a scene material's k(T) curve, stepped by \
                     fs-conduction's Newton endpoint and refused off the curve's span, and one \
                     declared",
                    if card_capacity {
                        ""
                    } else {
                        "temperature-dependent heat capacity (a body takes only a material's \
                         k(T)), "
                    },
                )
            } else {
                (
                    "one declared conductivity and",
                    "temperature-dependent conductivity, ",
                )
            };
            description = std::borrow::Cow::Owned(format!(
                "{description}; conduction bodies step fs-conduction's backward-Euler \
                 finite-element field on an authored box with {conductivity} {capacity}, each \
                 port a whole box face taking a uniform flux and returning the face's area-mean \
                 end temperature, the other faces insulated; the entropy generated inside a body \
                 is reported as zero, so the entropy totals are a lower bound; no internal \
                 source, {dependent}imported {without}"
            ));
            if contacts {
                // A card's value is measured, but it is read once unless
                // the contact follows its card.
                let any = |asked: fn(&conduction::Declared) -> bool| {
                    self.elements
                        .iter()
                        .any(|element| element.contact.as_ref().is_some_and(asked))
                };
                let carded = any(conduction::Declared::is_card_backed);
                let resistance = if any(conduction::Declared::follows_its_card) {
                    "a number declared in the scene or an interface card's value: read once at \
                     the declared state, or, for a contact that names the card coordinates its \
                     two traces' temperatures take, asked at every interface point of every \
                     step at those temperatures, with the step refused where the card has no \
                     answer; it follows no pressure field, the card's uncertainty is not \
                     propagated, and nothing checks that the two bodies are the card's two \
                     surfaces"
                } else if carded {
                    "a number declared in the scene or an interface card's value read once at \
                     the declared state: it follows neither the faces' temperatures nor any \
                     pressure field, the card's uncertainty is not propagated, and nothing \
                     checks that the two bodies are the card's two surfaces"
                } else {
                    "declared in the scene and not measured, with no dependence on temperature \
                     or pressure"
                };
                // A trace contact is a thermal contact, so the list of what
                // the scene has none of names the mechanical one (unless a
                // contact interface took contact off the list); a card
                // behind the contact is a card resolution. No comma after
                // "friction": the adaptive sentence goes on "or phase change".
                let mut base =
                    description.replace("contact, friction", "mechanical contact, friction");
                if carded {
                    base = base.replace(
                        "no material-card resolution,",
                        "no material-card resolution except a trace contact's interface card,",
                    );
                }
                description = std::borrow::Cow::Owned(format!(
                    "{base}; a trace carries one heat load and one end temperature for \
                     every vertex of its face, and a trace contact carries heat between two \
                     traces at fs-conduction's interface points through one area-specific \
                     resistance for the whole patch, {resistance}; a contact generates entropy \
                     where the heat crosses, the ledger books the nodal sum at its ports, and \
                     no order is claimed for the gap between the two"
                ));
            }
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::EnthalpyBody)
        {
            let base = description
                .replace(", phase change", "")
                .replace(" or phase change", "");
            description = std::borrow::Cow::Owned(format!(
                "{base}; enthalpy bodies use authored bounded sensible and latent storage \
                 with entropy-conjugate heat ports at fixed mass and negligible pressure-volume \
                 work; no material-pack provenance, phase kinetics or mechanical expansion"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::ThermalTimeline)
        {
            let base = description
                .replace(", phase change", "")
                .replace(" or phase change", "");
            description = std::borrow::Cow::Owned(format!(
                "{base}; thermal timelines are stepped by the fs-thermal owner's explicit rule \
                 under its own stability bound, through authored link conductances, films, \
                 reservoirs and sources, on charts that confer no material-pack provenance unless \
                 a declared source supplies them; the entropy generated inside a timeline is \
                 reported as zero, so the entropy totals are a lower bound; every cell is born at \
                 the start; no implicit conduction, radiation, contact geometry or expansion"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::LayerConductor)
        {
            let base = description.replace(
                "no material-card resolution, contact,",
                "no material-card resolution, geometric contact,",
            );
            description = std::borrow::Cow::Owned(format!(
                "{base}; layer conductors integrate authored bounded k(T) charts with fixed \
                 planar or radial geometry and areal contact resistance; no verified \
                 material-pack provenance, internal layer storage or source extrapolation"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::SpeedController)
        {
            description = std::borrow::Cow::Owned(format!(
                "{description}; speed_controller executes fs-control's sampled PI policy \
                 with authored torque and power limits and holds each command until the \
                 next sample; ideal sampled measurements and torque actuation, with no \
                 machine current allocation, sensor dynamics or stability claim"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| matches!(element.kind, Kind::ControlledEffort | Kind::ControlledFlow))
        {
            description = std::borrow::Cow::Owned(format!(
                "{description}; controlled sources exchange signed physical work while their \
                 scalar commands carry no energy and are held over each step in the \
                 declared SI unit; no saturation, power-supply or controller stability \
                 model is inferred"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::ModulatedTransformer)
        {
            let base = description.replace(
                "constant parameters or",
                "constant parameters, signal-modulated lossless ratios or",
            );
            description = std::borrow::Cow::Owned(format!(
                "{base}; transformer ratios are held over each step in declared SI conversion \
                 units; no crank/cam geometry, displacement constraint, ratio-actuator work \
                 or control stability is inferred"
            ));
        }
        if self.elements.iter().any(|element| {
            matches!(
                element.kind,
                Kind::AffineSignal
                    | Kind::SignalSource
                    | Kind::SampledController
            )
        }) {
            description = std::borrow::Cow::Owned(format!(
                "{description}; control gains, units and memory are authored; measurements and \
                 reference waveforms are sampled from accepted start-of-step state, and sampled \
                 controllers update only at their declared clock events; ideal actuator commands \
                 are held and their signed physical work is included in the ledger; no saturation, \
                 power-supply dynamics or control stability is inferred"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::CarreauPipe)
        {
            let base = description.replace(
                "constant parameters or",
                "constant parameters, authored Carreau-Yasuda pipe laws or",
            );
            description = std::borrow::Cow::Owned(format!(
                "{base}; a carreau_pipe is the fs-material quasi-steady, fully developed \
                 laminar pipe at one uniform temperature, with its authored temperature shift; \
                 no fluid inertia, compressibility, entrance loss, transition or turbulence, \
                 and no measured-material claim"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::PortHamiltonian)
        {
            let base = description.replace(
                "constant parameters or",
                "constant parameters, authored linear port-Hamiltonian matrices or",
            );
            description = std::borrow::Cow::Owned(format!(
                "{base}; a port_hamiltonian is stepped by the fs-system exact owner, whose map \
                 is exact only for inputs held over each step: inputs from other elements are \
                 the kernel's step means, so a coupled scene is not exact; quadratic storage \
                 only, no nonlinear storage or switching"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| !element.material_sources.is_empty())
        {
            let base = description
                .replace("no material-card resolution except a trace contact's interface card, ", "source resolution is limited to declared pure-phase thermal laws and trace contacts' interface cards; no ")
                .replace("no material-card resolution, ", "source resolution is limited to declared pure-phase thermal laws; no ")
                .replace("enthalpy bodies use authored bounded", "enthalpy bodies use authored or source-resolved bounded")
                .replace("layer conductors integrate authored bounded", "layer conductors integrate authored or source-resolved bounded")
                .replace("no material-pack provenance, phase kinetics", "authored charts carry no material-pack provenance; no phase kinetics")
                .replace("no verified material-pack provenance, internal layer storage", "authored charts carry no verified material-pack provenance; no internal layer storage");
            description = std::borrow::Cow::Owned(format!(
                "{base}; source-backed thermal components retain original pack, card, \
                 caloric and claim identities and use nominal Cp(T), density(T) and k(T) \
                 on the declared fixed context; source uncertainties are retained but \
                 not propagated, and source resolution does not imply experimental qualification"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.kind == Kind::GasCavity)
        {
            let base = description
                .replace(
                    "no material-card resolution, ",
                    "no material-card resolution besides gas cavities' packs, ",
                )
                .replace(
                    "limited to declared pure-phase thermal laws; no ",
                    "limited to declared pure-phase thermal laws and gas cavities' packs; no ",
                );
            let claim = if self.elements.iter().any(|element| {
                element.kind == Kind::GasCavity && gas::stream_basis(&element.ports).is_some()
            }) {
                "gas_cavity calorics and pressure use the original NASA-9 and species packs; \
                 streaming cavities use fs-gas's native inventory and energy transfer, with \
                 first-order transfer and independently audited total and species mass; restrictions \
                 use the donor's enthalpy and composition, and exact source/enthalpy identity \
                 is required at every stream connection; transport packs bound admitted states \
                 but do not supply the authored restriction conductance; stream entropy is not \
                 carried or audited, and no reactions, condensation, real-gas equation of \
                 state, internal waves or wall heat capacity is modeled"
            } else {
                "a gas_cavity holds a closed, uniform, frozen ideal-gas mixture whose \
                 calorics and pressure are its NASA-9 and species packs', stepped by the fs-gas \
                 discrete-gradient heat and work solve; its transport packs are resolved at \
                 every step end and bound the admitted states though no transport is used; no \
                 reactions, condensation, real-gas equation of state, internal waves, leakage \
                 or wall heat capacity"
            };
            description = std::borrow::Cow::Owned(format!("{base}; {claim}"));
        }
        if self
            .elements
            .iter()
            .any(|element| !element.solids.is_empty())
        {
            let base = description
                .replace(
                    "no material-card resolution besides ",
                    "no material-card resolution besides solid bodies' densities and ",
                )
                .replace(
                    "no material-card resolution, ",
                    "no material-card resolution besides solid bodies' densities, ",
                )
                .replace(
                    "limited to declared pure-phase thermal laws",
                    "limited to declared pure-phase thermal laws, solid bodies' densities",
                );
            let expanded = self
                .elements
                .iter()
                .flat_map(|element| &element.solids)
                .any(|(_, report)| solid::expanded(report));
            description = std::borrow::Cow::Owned(if expanded {
                format!(
                    "{base}; a solid body's mass and principal inertia are the closed-form \
                     volume moments of its declared shape, centred on its origin along its axes, \
                     times the nominal density its card gives at the declared state: a uniform, \
                     rigid body of those dimensions. A solid given a temperature_k is that body \
                     freely and uniformly expanded about its origin by its card's linear \
                     expansion coefficient integrated from the state's temperature, its mass \
                     kept: no thermal stress, gradient or constraint, no tolerance or density \
                     field, and neither the density's nor the coefficient's uncertainty \
                     propagated"
                )
            } else {
                format!(
                    "{base}; a solid body's mass and principal inertia are the closed-form \
                     volume moments of its declared shape, centred on its origin along its axes, \
                     times the nominal density its card gives at the declared state: a uniform, \
                     rigid body of those dimensions, with no thermal expansion, tolerance or \
                     density field, and the density's uncertainty not propagated"
                )
            });
        }
        if self.elements.iter().any(|element| element.ohmic.is_some()) {
            // Last, so that every list the earlier sentences rewrote is
            // there to extend. A card-backed resistor names a material, so
            // the source-resolution sentence above is present.
            let base = description.replacen(
                "constant parameters",
                "constant parameters, card-backed resistances",
                1,
            );
            let base = if base.contains("pure-phase thermal laws; no ") {
                base.replace(
                    "pure-phase thermal laws; no ",
                    "pure-phase thermal laws and resistors' resistivity; no ",
                )
            } else {
                base.replace(
                    "pure-phase thermal laws",
                    "pure-phase thermal laws, resistors' resistivity",
                )
            };
            description = std::borrow::Cow::Owned(format!(
                "{base}; a card-backed resistor is fs-material's Ohmic node of fixed length and \
                 cross-section, its resistance the card's nominal rho(T) L / A at the \
                 temperature its heat port receives and refused off the declared span: DC \
                 conduction with no thermal expansion, skin effect, contact resistance or heat \
                 storage of its own, and the resistivity's uncertainty not propagated"
            ));
        }
        if self
            .elements
            .iter()
            .any(|element| element.laminar.is_some())
        {
            // After the resistors' sentence, for the same reason. A
            // card-backed pipe names a material too.
            let base = description.replacen(
                "constant parameters",
                "constant parameters, card-backed pipes",
                1,
            );
            let base = if base.contains("pure-phase thermal laws; no ") {
                base.replace(
                    "pure-phase thermal laws; no ",
                    "pure-phase thermal laws and pipes' viscosity; no ",
                )
            } else {
                base.replace(
                    "pure-phase thermal laws",
                    "pure-phase thermal laws, pipes' viscosity",
                )
            };
            description = std::borrow::Cow::Owned(format!(
                "{base}; a card-backed pipe is fs-material's laminar circular pipe of fixed \
                 radius and length, its conductance pi r^4 / (8 mu L) with the card's nominal \
                 viscosity at the temperature its heat port receives and refused off the \
                 declared span: fully developed, incompressible Hagen-Poiseuille flow, laminar \
                 by assumption with no Reynolds number checked, all of its viscous heat at that \
                 one temperature, no entrance length, wall slip, pressure-dependent viscosity, \
                 heat storage or heat carried by the flow, and the viscosity's uncertainty not \
                 propagated"
            ));
        }
        description
    }

    fn is_finished(&self) -> bool {
        self.duration_s - self.system.time_s() <= 1.0e-12 * self.duration_s
    }

    fn advance(&mut self, max_interval_s: f64, cx: &Cx<'_>) -> Result<bool, RunError> {
        if self.is_finished() {
            return Ok(false);
        }
        cx.checkpoint()
            .map_err(|_| RunError::Step(SystemError::Cancelled))?;
        if self.system.steps() as f64 >= 2.0 * MAX_STEPS {
            return Err(RunError::TooManySteps);
        }
        let remaining = (self.duration_s - self.system.time_s()).min(max_interval_s);
        let time_s = if let Some(adaptive) = &mut self.adaptive {
            let mut config = adaptive.config;
            // A caller-imposed horizon is not an accuracy retry. Its own
            // retries must still respect this smaller floor.
            config.min_step_s = config.min_step_s.min(remaining);
            // A preceding short caller horizon may have produced a suggestion
            // below the configured floor. Restore this call's effective floor
            // before proposing another interval.
            let proposed = adaptive
                .suggested_dt_s
                .max(config.min_step_s)
                .min(self.step_s)
                .min(remaining);
            let report = self
                .system
                .step_adaptive(
                    proposed,
                    config,
                    |initial, coarse, refined| adaptive.metric(initial, coarse, refined),
                    cx,
                )
                .map_err(RunError::Adaptive)?;
            self.iterations += report.iterations;
            self.shortenings += report.shortenings;
            adaptive.record(&report);
            report.time_s
        } else {
            let report = self
                .system
                .step(self.step_s.min(remaining), cx)
                .map_err(RunError::Step)?;
            self.iterations += u64::from(report.iterations);
            self.shortenings += u64::from(report.shortenings);
            report.time_s
        };
        let slack = 1.0e-12 * self.duration_s;
        if time_s >= self.next_sample_s - slack {
            self.sample();
            let every = self.sample_every_s.unwrap_or(f64::INFINITY);
            while self.next_sample_s <= time_s + slack {
                self.next_sample_s += every;
            }
        }
        if self.is_finished()
            && self.sample_every_s.is_some()
            && self.samples.last().is_none_or(|(t, _)| *t < time_s)
        {
            self.sample();
        }
        Ok(true)
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

    /// Accepted state and the last step's outputs (means for power ports,
    /// point samples or held commands for control signals).
    fn values(&self, element: &Element) -> Vec<(String, f64)> {
        if matches!(
            element.kind,
            Kind::AffineSignal | Kind::SignalSource | Kind::SampledController
        ) {
            return self.signal_values(element);
        }
        if matches!(element.kind, Kind::RigidNetwork | Kind::ViscoelasticNetwork) {
            return self.network_values(element);
        }
        if element.kind == Kind::GasVolume {
            return self.gas_values(element);
        }
        if element.kind == Kind::GasCavity {
            return gas::values(&self.system, element.node, &element.ports);
        }
        if element.kind == Kind::ModulatedRestriction {
            return self.modulated_restriction_values(element);
        }
        if element.kind == Kind::ModalStructure {
            return modal::values(&self.system, element.node, &element.ports);
        }
        if element.kind == Kind::ConductionBody {
            return conduction::values(
                &self.system,
                element.node,
                &element.ports,
                element.enthalpy,
            );
        }
        if let Some(declared) = &element.contact {
            return conduction::wall_values(&self.system, element.node, declared);
        }
        if element.kind == Kind::ResistanceMatrix {
            return self.resistance_values(element);
        }
        if element.kind == Kind::PortHamiltonian {
            return self.phs_values(element);
        }
        if element.kind == Kind::DiracJunction {
            return self.dirac_values(element);
        }
        if element.kind == Kind::ThermalTimeline {
            let system = &self.system;
            let environment = match system.port(element.node, "environment") {
                Some(port) if system.steps() > 0 => Some(system.last_output(port)[0]),
                Some(_) => Some(f64::NAN),
                None => None,
            };
            return timeline::values(system.state::<timeline::State>(element.node), environment);
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

    /// Control outputs, followed by any accepted clock and discrete memory.
    fn signal_values(&self, element: &Element) -> Vec<(String, f64)> {
        let width = element
            .ports
            .iter()
            .find(|port| port.id().as_str() == "output")
            .map_or(0, PortSpec::output_width);
        let output = self
            .system
            .port(element.node, "output")
            .filter(|_| self.system.steps() > 0)
            .map(|port| self.system.last_output(port));
        let mut values: Vec<_> = (0..width)
            .map(|index| {
                (
                    format!("output{index}"),
                    output
                        .and_then(|values| values.get(index))
                        .copied()
                        .unwrap_or(f64::NAN),
                )
            })
            .collect();
        if let Some(clock) = self.system.state::<SignalClock>(element.node) {
            values.push(("time_s".to_owned(), clock.time_s()));
        }
        if let Some(state) = self.system.state::<ControllerState>(element.node) {
            values.extend(
                state
                    .coordinates()
                    .iter()
                    .enumerate()
                    .map(|(index, &value)| (format!("x{index}"), value)),
            );
            values.push(("sample_count".to_owned(), state.samples() as f64));
            values.push(("time_s".to_owned(), state.time_s()));
        }
        values
    }

    /// A resistance has no state: these are its last-step mean efforts and heat.
    fn resistance_values(&self, element: &Element) -> Vec<(String, f64)> {
        let Some(port) = element.ports.iter().find(|port| port.id().as_str() == "u") else {
            return Vec::new();
        };
        let Carrier::Power(kind) = port.carrier() else {
            return Vec::new();
        };
        let Some(domain) = Domain::of(kind) else {
            return Vec::new();
        };
        let output = |name: &str| {
            self.system
                .port(element.node, name)
                .filter(|_| self.system.steps() > 0)
                .map(|port| self.system.last_output(port))
        };
        let efforts = output("u");
        let mut values: Vec<_> = (0..port.width())
            .map(|row| {
                (
                    format!("u{row}.{}", domain.names().effort),
                    efforts
                        .and_then(|values| values.get(row))
                        .copied()
                        .unwrap_or(f64::NAN),
                )
            })
            .collect();
        let heat_w = output("heat")
            .and_then(|values| values.first())
            .map_or(f64::NAN, |value| -value);
        values.push(("heat_w".to_owned(), heat_w));
        values
    }

    /// The last step's conjugate effort and flow at every named Dirac port.
    fn dirac_values(&self, element: &Element) -> Vec<(String, f64)> {
        let mut values = Vec::with_capacity(2 * element.ports.len());
        for port in &element.ports {
            let Carrier::Power(kind) = port.carrier() else {
                continue;
            };
            let Some(domain) = Domain::of(kind) else {
                continue;
            };
            let name = port.id().as_str();
            let accepted = self
                .system
                .port(element.node, name)
                .filter(|_| self.system.steps() > 0);
            let effort = accepted
                .and_then(|port| self.system.last_input(port).first().copied())
                .unwrap_or(f64::NAN);
            let flow = accepted
                .and_then(|port| self.system.last_output(port).first().copied())
                .unwrap_or(f64::NAN);
            values.push((format!("{name}.{}", domain.names().effort), effort));
            values.push((format!("{name}.{}", domain.names().flow), flow));
        }
        values
    }

    /// A port-Hamiltonian element: every state coordinate, then for each
    /// power port the mean of what it returned over the last step (a flow
    /// for an effort_in port, an effort for a flow_in port), and the heat.
    fn phs_values(&self, element: &Element) -> Vec<(String, f64)> {
        let state = self.system.state::<Vec<f64>>(element.node);
        let mut values: Vec<_> = (0..element.coordinates)
            .map(|index| {
                let coordinate = state.and_then(|state| state.get(index)).copied();
                (format!("x{index}"), coordinate.unwrap_or(f64::NAN))
            })
            .collect();
        let output = |name: &str| {
            self.system
                .port(element.node, name)
                .filter(|_| self.system.steps() > 0)
                .and_then(|port| self.system.last_output(port).first().copied())
                .unwrap_or(f64::NAN)
        };
        for port in &element.ports {
            let Carrier::Power(kind) = port.carrier() else {
                continue;
            };
            let Some(domain) = Domain::of(kind) else {
                continue;
            };
            let names = domain.names();
            let returned = if port.role() == PortRole::EffortIn {
                names.flow
            } else {
                names.effort
            };
            let name = port.id().as_str();
            values.push((format!("{name}.{returned}"), output(name)));
        }
        values.push(("heat_w".to_owned(), -output("heat")));
        values
    }

    #[allow(clippy::too_many_lines)] // one report row per element kind
    fn lumped_values(&self, element: &Element) -> Vec<(&'static str, f64)> {
        let system = &self.system;
        // Before the first step there is no last step: those values are null.
        let output = |port: &str| match system.port(element.node, port) {
            Some(port) if system.steps() > 0 => system.last_output(port)[0],
            _ => f64::NAN,
        };
        let input = |port: &str| match system.port(element.node, port) {
            Some(port) if system.steps() > 0 => system.last_input(port)[0],
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
            Kind::NodeCapacitor(domain) => {
                let names = domain.names();
                let content = system
                    .state::<f64>(element.node)
                    .copied()
                    .unwrap_or(f64::NAN);
                vec![
                    (names.extension, content),
                    (
                        names.effort,
                        content / element.node_capacity.unwrap_or(f64::NAN),
                    ),
                ]
            }
            Kind::IdealDiode => vec![("current_a", input("p")), ("voltage_v", output("p"))],
            Kind::IdealSwitch => vec![
                ("current_a", input("p")),
                ("voltage_v", output("p")),
                ("closed", input("closed")),
            ],
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
            Kind::EnthalpyBody => {
                let state = system.state::<EquilibriumPhaseState>(element.node);
                vec![
                    ("specific_enthalpy_j_kg", state.map_or(f64::NAN, |s| s.specific_enthalpy_j_kg())),
                    ("temperature_k", state.map_or(f64::NAN, |s| s.temperature_k())),
                    ("liquid_mass_fraction", state.map_or(f64::NAN, |s| s.liquid_mass_fraction())),
                    ("bulk_density_kg_m3", state.map_or(f64::NAN, |s| s.bulk_density_kg_m3())),
                ]
            }
            Kind::Source(_) | Kind::VelocitySource | Kind::CurrentSource(_) | Kind::Heater => {
                vec![(
                    "energy_in_j",
                    system
                        .state::<SourceAccount>(element.node)
                        .map_or(f64::NAN, |account| account.energy_in_j),
                )]
            }
            Kind::ControlledEffort | Kind::ControlledFlow => {
                let mut values = vec![(
                    "energy_in_j",
                    system.state::<SourceAccount>(element.node)
                        .map_or(f64::NAN, |account| account.energy_in_j),
                )];
                if let Some(quantity) = controlled_quantity(element) {
                    values.push((quantity, output("p")));
                }
                values.push(("command", output("p")));
                values
            }
            Kind::SpeedController => {
                let state = system.state::<SampledSpeedState>(element.node);
                vec![
                    (
                        "integral_torque_nm",
                        state.map_or(f64::NAN, SampledSpeedState::integral_torque_nm),
                    ),
                    (
                        "held_torque_nm",
                        state.map_or(f64::NAN, SampledSpeedState::held_torque_nm),
                    ),
                    (
                        "sample_phase_s",
                        state.map_or(f64::NAN, SampledSpeedState::phase_s),
                    ),
                    (
                        "sample_count",
                        state.map_or(f64::NAN, |state| state.samples() as f64),
                    ),
                ]
            }
            Kind::Resistance(domain) => {
                let mut values = vec![
                    (domain.names().effort, output("a")),
                    ("heat_w", -output("heat")),
                ];
                values.extend(ohmic_values(system, element));
                values
            }
            Kind::Slider => {
                let state = system.state::<Slip>(element.node);
                vec![
                    ("slip_m_s", state.map_or(f64::NAN, |state| -state.slip_m_s())),
                    ("stuck", state.map_or(f64::NAN, |state| if state.is_stuck() { 1.0 } else { 0.0 })),
                    ("force_n", output("a")),
                    ("heat_w", -output("heat")),
                ]
            }
            Kind::FrictionInterface => {
                let mut values = friction::values(system.state::<friction::State>(element.node));
                values.extend([("force_n", output("a")), ("heat_w", -output("heat"))]);
                values
            }
            // The last step's traction entering the interface, both axes.
            Kind::SlidingInterface => {
                let traction = match system.port(element.node, "traction") {
                    Some(port) if system.steps() > 0 => system.last_output(port),
                    _ => &[f64::NAN, f64::NAN][..],
                };
                let mut values =
                    friction::sliding_values(system.state::<friction::Sliding>(element.node));
                values.extend([
                    ("traction_x_n", traction[0]),
                    ("traction_y_n", traction[1]),
                    ("heat_w", -output("heat")),
                ]);
                values
            }
            Kind::ContactInterface => contact::values(system.state::<contact::State>(element.node)),
            // The wrench the last step's patch acted with on side a: its port
            // returns minus it, the third coefficient over the radius.
            Kind::HertzPatch => {
                let wrench = match system.port(element.node, "wrench") {
                    Some(port) if system.steps() > 0 => system.last_output(port),
                    _ => &[f64::NAN; 3][..],
                };
                let radius_m = element.radius_m.unwrap_or(f64::NAN);
                let mut values =
                    friction::hertz_values(system.state::<friction::Patch>(element.node));
                values.extend([
                    ("force_x_n", -wrench[0]),
                    ("force_y_n", -wrench[1]),
                    ("spin_moment_nm", -wrench[2] * radius_m),
                    ("heat_w", -output("heat")),
                ]);
                values
            }
            Kind::ResistorBranch(domain) => {
                let mut values = vec![
                    (domain.names().flow, output("a")),
                    ("heat_w", -output("heat")),
                ];
                values.extend(ohmic_values(system, element));
                values.extend(laminar_values(system, element));
                values
            }
            // The last step's flow from a to b, its heat, and the
            // temperature the liquid's viscosity was read at.
            Kind::CarreauPipe => {
                let temperature = match system.port(element.node, "heat") {
                    Some(port) if system.steps() > 0 => system.last_input(port)[0],
                    _ => f64::NAN,
                };
                vec![
                    ("flow_m3_s", output("a")),
                    ("heat_w", -output("heat")),
                    ("temperature_k", temperature),
                ]
            }
            Kind::Gyrator => vec![("effort_a", output("a")), ("effort_b", output("b"))],
            Kind::ModulatedTransformer => {
                let ratio = match system.port(element.node, "ratio") {
                    Some(port) if system.steps() > 0 => system.last_input(port)[0],
                    _ => f64::NAN,
                };
                vec![("ratio", ratio), ("effort_a", output("a")), ("flow_b", output("b"))]
            }
            Kind::Conductance | Kind::LayerConductor => vec![("heat_rate_w", output("a"))],
            // The temperature the node's balance gave it in the last step.
            Kind::HeatNode => vec![("temperature_k", output("source"))],
            // What crossed from `a` to `b` during the last step.
            Kind::Restriction => {
                let crossing = match system.port(element.node, "a") {
                    Some(port) if system.steps() > 0 => system.last_output(port),
                    _ => &[f64::NAN, f64::NAN][..],
                };
                vec![("mass_rate_kg_s", crossing[0]), ("enthalpy_rate_w", crossing[1])]
            }
            // Read by their vector-valued report functions.
            Kind::RigidNetwork
            | Kind::AffineSignal
            | Kind::SignalSource
            | Kind::SampledController
            | Kind::ViscoelasticNetwork
            | Kind::GasVolume
            | Kind::ModulatedRestriction
            | Kind::ResistanceMatrix
            | Kind::DiracJunction
            | Kind::PortHamiltonian
            | Kind::ModalStructure
            | Kind::ThermalTimeline
            | Kind::ConductionBody
            | Kind::TraceContact
            | Kind::GasCavity => Vec::new(),
        }
    }

    /// The opening actually used and native transport from a to b.
    /// Stream coordinates retain their source catalog order.
    fn modulated_restriction_values(&self, element: &Element) -> Vec<(String, f64)> {
        let crossing = self
            .system
            .port(element.node, "a")
            .filter(|_| self.system.steps() > 0)
            .map(|port| self.system.last_output(port));
        let rate = |index: usize| {
            crossing
                .and_then(|values| values.get(index).copied())
                .unwrap_or(f64::NAN)
        };
        let opening = self
            .system
            .port(element.node, "opening")
            .filter(|_| self.system.steps() > 0)
            .and_then(|port| self.system.last_input(port).first().copied())
            .unwrap_or(f64::NAN);
        let mut values = vec![
            ("mass_rate_kg_s".to_owned(), rate(0)),
            ("enthalpy_rate_w".to_owned(), rate(1)),
            ("opening".to_owned(), opening),
        ];
        let carried = element
            .ports
            .iter()
            .find(|port| port.id().as_str() == "a")
            .map_or(0, PortSpec::carried_quantities);
        for index in 0..carried {
            values.push((
                format!("carried{index}_rate_kg_s"),
                rate(STREAM_CARRIED + index),
            ));
        }
        values
    }

    /// What a gas volume holds: pressure, temperature, mass, volume and the
    /// mass fraction of each carried quantity.
    fn gas_values(&self, element: &Element) -> Vec<(String, f64)> {
        let held = self.system.state::<GasContents>(element.node);
        let (Some(held), Some(gas)) = (held, element.gas.as_ref()) else {
            return Vec::new();
        };
        let mut values = vec![
            ("pressure_pa".to_owned(), gas.pressure_pa(held)),
            ("temperature_k".to_owned(), gas.temperature_k(held)),
            ("mass_kg".to_owned(), held.mass_kg),
            ("volume_m3".to_owned(), held.volume_m3),
        ];
        for (index, carried) in held.carried_kg.iter().enumerate() {
            values.push((format!("carried{index}"), carried / held.mass_kg));
        }
        values
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
        stream_ledger_text(&mut out, ledger);
        if !self.samples.is_empty() {
            let _ = writeln!(out, "\n{} samples (use --json to read them)", self.samples.len());
        }
        if !self.materials.is_empty() {
            let _ = writeln!(
                out,
                "\nMaterial sources (nominal values; uncertainty not propagated):"
            );
            for (name, material) in &self.materials {
                let _ = writeln!(
                    out,
                    "{name}: caloric {}, {} source queries",
                    material.caloric.identity().to_hex(),
                    material.receipts().count()
                );
            }
            for element in &self.elements {
                if !element.material_sources.is_empty() {
                    let _ = writeln!(
                        out,
                        "{} <- {}",
                        element.name,
                        element.material_sources.join(", ")
                    );
                }
            }
        }
        out.push_str(&solid::text_report(&self.elements));
        out.push_str(&conduction::text_report(&self.elements));
        let _ = writeln!(out, "\nNo claim: {}.", self.no_claim());
        out
    }

    fn render_json(&self, path: &Path, checkpoint: Option<(&Path, &str)>) -> String {
        let system = &self.system;
        let ledger = system.ledger();
        let (blocks, loops) = system.plan_shape();
        let mut out = String::new();
        let _ = write!(
            out,
            "{{\"schema\":{},\"scene\":{},\"time_s\":{},\"steps\":{},\"iterations\":{},\"shortenings\":{},\"ledger_steps\":{},\"plan\":{{\"blocks\":{blocks},\"loops\":{loops}}},\"elements\":[",
            quote(RESULT_SCHEMA),
            quote(&path.display().to_string()),
            num(system.time_s()),
            system.steps(),
            self.iterations,
            self.shortenings,
            ledger.steps
        );
        for (index, (element, account)) in self.elements.iter().zip(&ledger.nodes).enumerate() {
            if index > 0 {
                out.push(',');
            }
            let mass = account.mass.as_ref().map_or_else(
                || "null".into(),
                |mass| {
                    inventory_json(
                        mass.port_mass_in_kg.value(),
                        mass.stored_change_kg.value(),
                        mass.balance_defect_kg.value(),
                        mass.audited_steps,
                    )
                },
            );
            let _ = write!(
                out,
                "{{\"name\":{},\"kind\":{},\"state\":{},\"ledger\":{{\"port_energy_in_j\":{},\"stored_change_j\":{},\"dissipated_j\":{},\"numerical_loss_j\":{},\"entropy_generated_j_per_k\":{},\"balance_defect_j\":{},\"has_streams\":{},\"mass\":{},\"carried_kg\":{}}}",
                quote(&element.name),
                quote(element.kind.name()),
                object(self.values(element).into_iter()),
                num(account.port_energy_in_j.value()),
                num(account.stored_change_j.value()),
                num(account.dissipated_j.value()),
                num(account.numerical_loss_j.value()),
                num(account.entropy_generated_j_per_k.value()),
                num(account.balance_defect_j.value()),
                account.has_streams,
                mass,
                carried_accounts_json(&account.carried)
            );
            // A trace contact's resistance is the scene's declaration or a
            // card's value: the report repeats it, with its source or with
            // the card's identities, beside the contact's rows.
            if let Some(declared) = &element.contact {
                out.push_str(&declared.json_member());
            }
            if element.kind == Kind::GasCavity && gas::stream_basis(&element.ports).is_some() {
                let species = gas::species_names(&element.ports)
                    .into_iter()
                    .map(quote)
                    .collect::<Vec<_>>()
                    .join(",");
                let _ = write!(out, ",\"stream_species\":[{species}]");
            }
            out.push('}');
        }
        out.push_str("],\"boundaries\":[");
        for (index, boundary) in ledger.boundaries.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"port\":{},\"energy_in_j\":{},\"entropy_in_j_per_k\":{},\"mass_in_kg\":{},\"carried_in_kg\":{}}}",
                quote(&boundary.port),
                num(boundary.energy_in_j.value()),
                num(boundary.entropy_in_j_per_k.value()),
                num(boundary.mass_in_kg.value()),
                sums_json(&boundary.carried_in)
            );
        }
        out.push_str("],\"connections\":[");
        for (index, connection) in ledger.connections.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"kernel_index\":{index},\"energy_defect_j\":{},\"mass_defect_kg\":{},\"carried_defect_kg\":{}}}",
                num(connection.energy_defect_j.value()),
                num(connection.mass_defect_kg.value()),
                sums_json(&connection.carried_defect)
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
        let carried_residuals = (0..carried_count(ledger))
            .map(|index| optional_num(ledger.carried_residual(index)))
            .collect::<Vec<_>>()
            .join(",");
        let _ = write!(
            out,
            "],\"totals\":{{\"stored_change_j\":{},\"boundary_energy_in_j\":{},\"sink_heat_j\":{},\"numerical_loss_j\":{},\"energy_residual_j\":{},\"connection_defect_j\":{},\"entropy_generated_j_per_k\":{},\"mass_residual_kg\":{},\"carried_residual_kg\":[{carried_residuals}]}},\"samples\":[",
            num(ledger.stored_change_j()),
            num(ledger.boundary_energy_in_j()),
            num(ledger.sink_heat_j()),
            num(ledger.numerical_loss_j()),
            num(ledger.energy_residual_j()),
            num(ledger.connection_defect_j()),
            num(ledger.entropy_generated_j_per_k()),
            optional_num(ledger.mass_residual_kg())
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
                adaptive.config.order,
                adaptive.attempts,
                adaptive.trial_steps,
                num(adaptive.last_error_norm),
                num(adaptive.max_error_norm),
                num(adaptive.suggested_dt_s)
            );
        }
        if !self.materials.is_empty() {
            let reports = self
                .materials
                .values()
                .map(|source| source.report_json())
                .collect::<Vec<_>>()
                .join(",");
            let bindings = self
                .elements
                .iter()
                .filter(|element| !element.material_sources.is_empty())
                .map(|element| {
                    format!(
                        "{}:[{}]",
                        quote(&element.name),
                        element
                            .material_sources
                            .iter()
                            .map(|name| quote(name))
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            let _ = write!(
                out,
                ",\"materials\":[{reports}],\"material_bindings\":{{{bindings}}}"
            );
        }
        out.push_str(&solid::json_report(&self.elements));
        if let Some((destination, root)) = checkpoint {
            let _ = write!(
                out,
                ",\"checkpoint\":{{\"path\":{},\"resume_root\":{}}}",
                quote(&destination.display().to_string()),
                quote(root)
            );
        }
        let _ = writeln!(out, ",\"no_claim\":{}}}", quote(&self.no_claim()));
        out
    }
}

// Scene stream coordinates are gas mass fractions, so their inventories
// and transports are kg. The kernel's generic carried coordinates do not
// imply these units or attach a species identity to an ordinal.
fn inventory_json(port_in: f64, stored_change: f64, defect: f64, audited_steps: u64) -> String {
    format!(
        "{{\"port_in_kg\":{},\"stored_change_kg\":{},\"balance_defect_kg\":{},\"audited_steps\":{audited_steps}}}",
        num(port_in),
        num(stored_change),
        num(defect)
    )
}

fn carried_accounts_json(accounts: &[Option<CarriedAccount>]) -> String {
    let entries = accounts
        .iter()
        .map(|account| {
            account.as_ref().map_or_else(
                || "null".into(),
                |account| {
                    inventory_json(
                        account.port_in.value(),
                        account.stored_change.value(),
                        account.balance_defect.value(),
                        account.audited_steps,
                    )
                },
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{entries}]")
}

fn sums_json(sums: &[Account]) -> String {
    format!(
        "[{}]",
        sums.iter()
            .map(|sum| num(sum.value()))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn carried_count(ledger: &Ledger) -> usize {
    ledger
        .nodes
        .iter()
        .map(|node| node.carried.len())
        .max()
        .unwrap_or(0)
}

fn optional_num(value: Option<f64>) -> String {
    value.map_or_else(|| "null".into(), num)
}

fn stream_ledger_text(out: &mut String, ledger: &Ledger) {
    if !ledger.nodes.iter().any(|node| node.has_streams) {
        return;
    }
    let _ = writeln!(
        out,
        "\nstream ledger [kg]: measured/accepted steps; carried ordinals are gas constituent masses"
    );
    for node in ledger.nodes.iter().filter(|node| node.has_streams) {
        if let Some(mass) = &node.mass {
            let _ = writeln!(
                out,
                "{} mass: port in {:.6e}, stored change {:.6e}, defect {:.3e}; {}/{} steps",
                node.name,
                mass.port_mass_in_kg.value(),
                mass.stored_change_kg.value(),
                mass.balance_defect_kg.value(),
                mass.audited_steps,
                ledger.steps
            );
        } else {
            let _ = writeln!(
                out,
                "{} mass: unavailable; 0/{} steps",
                node.name, ledger.steps
            );
        }
        for (index, account) in node.carried.iter().enumerate() {
            if let Some(account) = account {
                let _ = writeln!(
                    out,
                    "{} carried[{index}]: port in {:.6e}, stored change {:.6e}, defect {:.3e}; {}/{} steps",
                    node.name,
                    account.port_in.value(),
                    account.stored_change.value(),
                    account.balance_defect.value(),
                    account.audited_steps,
                    ledger.steps
                );
            } else {
                let _ = writeln!(
                    out,
                    "{} carried[{index}]: unavailable; 0/{} steps",
                    node.name, ledger.steps
                );
            }
        }
    }
    for boundary in &ledger.boundaries {
        let _ = writeln!(
            out,
            "boundary {}: mass in {:.6e} kg, carried in {} kg",
            boundary.port,
            boundary.mass_in_kg.value(),
            sums_json(&boundary.carried_in)
        );
    }
    for (index, connection) in ledger.connections.iter().enumerate() {
        let _ = writeln!(
            out,
            "kernel connection {index}: mass defect {:.3e} kg, carried defects {} kg",
            connection.mass_defect_kg.value(),
            sums_json(&connection.carried_defect)
        );
    }
    let residual = |value: Option<f64>| {
        value.map_or_else(
            || "unavailable (incomplete measurement coverage)".into(),
            |value| format!("{value:.3e} kg"),
        )
    };
    let _ = writeln!(
        out,
        "stream mass residual: {}",
        residual(ledger.mass_residual_kg())
    );
    for index in 0..carried_count(ledger) {
        let _ = writeln!(
            out,
            "carried[{index}] residual: {}",
            residual(ledger.carried_residual(index))
        );
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

#[cfg(test)]
mod tests {
    use super::{J, Kind};

    /// Claims, declared before the first run: every kind name this command
    /// lowers is the `kind` of an element of some scene in `examples/system`.
    /// The checkpoint tests save and resume every scene there, so a kind
    /// with an example is a kind whose save and resume a test shows. A kind
    /// added to the table fails here, by name, until a scene uses it.
    #[test]
    fn every_kind_is_the_kind_of_an_element_of_some_example_scene() {
        let directory = concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/system");
        let mut used = std::collections::BTreeSet::new();
        let mut scenes = 0;
        for entry in std::fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().is_none_or(|kind| kind != "json") {
                continue;
            }
            let root = J::parse(&std::fs::read_to_string(&path).unwrap()).unwrap();
            let elements = root.get("elements").and_then(J::as_object);
            for (_, definition) in elements.expect("a scene's elements") {
                used.insert(definition.str_field("kind").unwrap().to_owned());
            }
            scenes += 1;
        }
        let unused: Vec<&str> = Kind::TABLE
            .iter()
            .map(|(name, _)| *name)
            .filter(|name| !used.contains(*name))
            .collect();
        assert!(
            unused.is_empty(),
            "no element of the {scenes} scenes in examples/system has kind {unused:?}: \
             add a scene that uses each, so that a test saves and resumes it"
        );
    }
}
