//! Unfrozen canonical assembly model draft (bead `frankensim-rc-1109-root-3gc12.11.1`).
//!
//! Extends the frozen `scenario-assembly-v1` vocabulary without modifying v1
//! wire bytes or `scenario.ir`. Covers:
//! - Parts (`RigidBody`, `FlexibleBody`, `SpringBeam`, plus all v1 part kinds)
//!   with explicit per-part physics sets and fidelity policies (rungs + error budget);
//! - Kinematic joints (`Fixed`, `Revolute`, `Prismatic`, `Helical`) and gear-mesh
//!   relations;
//! - Card-bound interfaces (power ports, mechanical/thermal/electrical contacts,
//!   glue) with surface-feature pairs, explicit state, and third-body lubricants;
//! - Sources with typed time histories on ports or parameters, signal-port
//!   controllers, and an explicit time plan (horizon, step policy, sample rate);
//! - Extended environment (applied magnetic field, relative humidity, gas composition);
//! - Probes (microphone receiver model, thermocouple, timing, strain gauge, field);
//! - Parametric mechanical/primitive generators (`Gear`, `SpiralSpring`, `Primitive`)
//!   alongside `BoreProfile`, `OutlineWithHoles`, and `ArchingSurface`.
//!
//! Structural validation refuses invalid names, units, identities, bounds, or
//! references by name. It performs no solver lowering, meshing, or material
//! lookup (those belong to the lowering registry in bead `.11.2`).

#![deny(unsafe_code)]

use std::collections::BTreeSet;
use std::fmt;

use fs_blake3::{ContentHash, hash_bytes};
use fs_qty::{Dims, QtyAny};

use crate::assembly::{
    DimensionBound, GeometrySource, MaterialRef, ParameterKey, ParameterPath, ReferenceState,
    UnknownProperty,
};
use crate::entity::{Entity, EntityCatalog, EntityId, EntityKind};
use crate::{Environment, FrameId, Scenario};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);
const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const ACCEL: Dims = Dims([1, 0, -2, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const MAGNETIC_FIELD_TESLA: Dims = Dims([0, 1, -2, 0, -1, 0]);

/// Extended physical part kinds covering mechanical, acoustic, thermal, fluid,
/// and coupled multi-domain devices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DraftPartKind {
    /// Vibrating string / wire.
    String,
    /// Thin shell / plate.
    Shell,
    /// Continuum solid (thermal, elastic, coupled).
    Solid,
    /// Acoustic or fluid duct.
    Duct,
    /// Enclosed fluid or acoustic cavity.
    Cavity,
    /// Valve or tone-hole pattern.
    ValvePattern,
    /// Prescribed / driven body.
    DrivenBody,
    /// 6-DOF rigid body.
    RigidBody,
    /// Flexible / modal / continuum-deformable body.
    FlexibleBody,
    /// Discrete spring, beam, reed, mainspring, or hairspring.
    SpringBeam,
}

impl DraftPartKind {
    /// Canonical stable label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::String => "string",
            Self::Shell => "shell",
            Self::Solid => "solid",
            Self::Duct => "duct",
            Self::Cavity => "cavity",
            Self::ValvePattern => "valve_pattern",
            Self::DrivenBody => "driven_body",
            Self::RigidBody => "rigid_body",
            Self::FlexibleBody => "flexible_body",
            Self::SpringBeam => "spring_beam",
        }
    }
}

/// Active physics domain owned by a part.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PhysicsDomain {
    /// Rigid or flexible body mechanics.
    Mechanics,
    /// Linear or nonlinear acoustics.
    Acoustics,
    /// Heat transfer / thermodynamics.
    Thermal,
    /// Compressible or incompressible fluid flow.
    Fluid,
    /// Electromagnetics / electromechanical fields.
    Electromagnetics,
    /// Contact mechanics, friction, lubrication, and wear.
    Tribology,
    /// Species transport / thermochemistry.
    Chemistry,
}

impl PhysicsDomain {
    /// Canonical stable label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Mechanics => "mechanics",
            Self::Acoustics => "acoustics",
            Self::Thermal => "thermal",
            Self::Fluid => "fluid",
            Self::Electromagnetics => "electromagnetics",
            Self::Tribology => "tribology",
            Self::Chemistry => "chemistry",
        }
    }
}

/// Discretisation / model fidelity rung available for a part.
#[derive(Debug, Clone, PartialEq)]
pub enum FidelityRung {
    /// Lumped-parameter ODE/DAE model.
    Lumped,
    /// Truncated modal basis.
    Modal {
        /// Number of retained modes (positive).
        mode_count: u32,
    },
    /// Reduced-order state-space basis.
    ReducedOrder {
        /// Basis dimension (positive).
        basis_dim: u32,
    },
    /// 1D finite-volume / finite-cell discretisation.
    OneDField {
        /// Number of spatial cells (positive).
        cells: u32,
    },
    /// Finite-element discretisation.
    FiniteElement {
        /// Polynomial order (`1..=6`).
        order: u8,
        /// Target element size [m] (positive finite).
        h_target: QtyAny,
    },
    /// Lattice-Boltzmann discretisation.
    LatticeBoltzmann {
        /// Lattice spacing [m] (positive finite).
        dx: QtyAny,
    },
    /// Boundary-element discretisation.
    BoundaryElement {
        /// Number of surface panels (positive).
        panels: u32,
    },
}

/// Explicit per-part accuracy and compute budget.
#[derive(Debug, Clone, PartialEq)]
pub struct FidelityErrorBudget {
    /// Target relative error tolerance in `(0, 1)`.
    pub relative_tol: f64,
    /// Target absolute error tolerance in the part's primary quantity units (>= 0).
    pub absolute_tol: QtyAny,
    /// Optional work-unit cap for this part.
    pub max_work_units: Option<u64>,
}

/// Per-part fidelity ladder and error budget.
#[derive(Debug, Clone, PartialEq)]
pub struct FidelityPolicy {
    /// Ordered hierarchy of admissible fidelity rungs (at least one).
    pub rungs: Vec<FidelityRung>,
    /// Accuracy and work budget governing rung selection/refinement.
    pub error_budget: FidelityErrorBudget,
}

/// Parametric involute gear descriptor (beside bore/outline/arching generators).
#[derive(Debug, Clone, PartialEq)]
pub struct GearGenerator {
    /// Number of teeth (`>= 3`).
    pub teeth: u32,
    /// Metric module [m] (positive finite).
    pub module: QtyAny,
    /// Axial face width [m] (positive finite).
    pub face_width: QtyAny,
    /// Pressure angle [rad] in `(0, pi/4]`.
    pub pressure_angle: QtyAny,
    /// Dimensionless profile shift coefficient `x` in `[-1.0, 1.0]`.
    pub profile_shift: f64,
    /// Central bore radius [m] (`>= 0`, strictly below root radius).
    pub bore_radius: QtyAny,
}

impl GearGenerator {
    /// Pitch radius `m * z / 2` [m].
    #[must_use]
    pub fn pitch_radius_m(&self) -> f64 {
        0.5 * self.module.value * f64::from(self.teeth)
    }

    /// Dedendum / root radius `(z/2 - 1.25 + x) * m` [m].
    #[must_use]
    pub fn root_radius_m(&self) -> f64 {
        (0.5 * f64::from(self.teeth) - 1.25 + self.profile_shift) * self.module.value
    }

    /// Convert to the v1 `GeometrySource::Parametric` parameter table.
    #[must_use]
    pub fn to_geometry_source(&self) -> GeometrySource {
        GeometrySource::Parametric {
            generator: "Gear.draft.v1".to_owned(),
            params: vec![
                (
                    ParameterKey::new("teeth", Dims::NONE),
                    QtyAny::dimensionless(f64::from(self.teeth)),
                ),
                (ParameterKey::new("module", LENGTH), self.module),
                (ParameterKey::new("face_width", LENGTH), self.face_width),
                (
                    ParameterKey::new("pressure_angle", Dims::NONE),
                    self.pressure_angle,
                ),
                (
                    ParameterKey::new("profile_shift", Dims::NONE),
                    QtyAny::dimensionless(self.profile_shift),
                ),
                (ParameterKey::new("bore_radius", LENGTH), self.bore_radius),
            ],
        }
    }
}

/// Parametric Archimedean spiral spring descriptor (mainspring / hairspring).
#[derive(Debug, Clone, PartialEq)]
pub struct SpiralSpringGenerator {
    /// Inner arbor radius [m] (positive finite).
    pub inner_radius: QtyAny,
    /// Outer barrel/collet radius [m] (strictly greater than `inner_radius`).
    pub outer_radius: QtyAny,
    /// Number of active coils (positive finite).
    pub turns: f64,
    /// Radial blade thickness [m] (positive finite, must fit inside radial span).
    pub blade_thickness: QtyAny,
    /// Axial blade height [m] (positive finite).
    pub blade_height: QtyAny,
}

impl SpiralSpringGenerator {
    /// Convert to the v1 `GeometrySource::Parametric` parameter table.
    #[must_use]
    pub fn to_geometry_source(&self) -> GeometrySource {
        GeometrySource::Parametric {
            generator: "SpiralSpring.draft.v1".to_owned(),
            params: vec![
                (ParameterKey::new("inner_radius", LENGTH), self.inner_radius),
                (ParameterKey::new("outer_radius", LENGTH), self.outer_radius),
                (
                    ParameterKey::new("turns", Dims::NONE),
                    QtyAny::dimensionless(self.turns),
                ),
                (
                    ParameterKey::new("blade_thickness", LENGTH),
                    self.blade_thickness,
                ),
                (ParameterKey::new("blade_height", LENGTH), self.blade_height),
            ],
        }
    }
}

/// Parametric solid primitive shapes.
#[derive(Debug, Clone, PartialEq)]
pub enum PrimitiveShape {
    /// Axis-aligned rectangular block `[lx, ly, lz]` [m].
    Box {
        /// X extent [m].
        lx: QtyAny,
        /// Y extent [m].
        ly: QtyAny,
        /// Z extent [m].
        lz: QtyAny,
    },
    /// Right circular cylinder [m].
    Cylinder {
        /// Radius [m].
        radius: QtyAny,
        /// Axial height [m].
        height: QtyAny,
    },
    /// Sphere [m].
    Sphere {
        /// Radius [m].
        radius: QtyAny,
    },
    /// Hollow circular tube [m].
    Tube {
        /// Inner radius [m].
        inner_radius: QtyAny,
        /// Outer radius [m] (strictly greater than `inner_radius`).
        outer_radius: QtyAny,
        /// Axial length [m].
        length: QtyAny,
    },
}

impl PrimitiveShape {
    /// Convert to the v1 `GeometrySource::Parametric` parameter table.
    #[must_use]
    pub fn to_geometry_source(&self) -> GeometrySource {
        match self {
            Self::Box { lx, ly, lz } => GeometrySource::Parametric {
                generator: "Primitive.Box.draft.v1".to_owned(),
                params: vec![
                    (ParameterKey::new("lx", LENGTH), *lx),
                    (ParameterKey::new("ly", LENGTH), *ly),
                    (ParameterKey::new("lz", LENGTH), *lz),
                ],
            },
            Self::Cylinder { radius, height } => GeometrySource::Parametric {
                generator: "Primitive.Cylinder.draft.v1".to_owned(),
                params: vec![
                    (ParameterKey::new("radius", LENGTH), *radius),
                    (ParameterKey::new("height", LENGTH), *height),
                ],
            },
            Self::Sphere { radius } => GeometrySource::Parametric {
                generator: "Primitive.Sphere.draft.v1".to_owned(),
                params: vec![(ParameterKey::new("radius", LENGTH), *radius)],
            },
            Self::Tube {
                inner_radius,
                outer_radius,
                length,
            } => GeometrySource::Parametric {
                generator: "Primitive.Tube.draft.v1".to_owned(),
                params: vec![
                    (ParameterKey::new("inner_radius", LENGTH), *inner_radius),
                    (ParameterKey::new("outer_radius", LENGTH), *outer_radius),
                    (ParameterKey::new("length", LENGTH), *length),
                ],
            },
        }
    }
}

/// Geometry source for a draft part: either any v1 `GeometrySource` (imported
/// chart or bore/outline/arching parametric table) or one of the mechanical /
/// primitive parametric generators.
#[derive(Debug, Clone, PartialEq)]
pub enum DraftGeometrySource {
    /// Existing v1 parametric or imported geometry source.
    Standard(GeometrySource),
    /// Parametric involute gear.
    Gear(GearGenerator),
    /// Parametric Archimedean spiral spring.
    SpiralSpring(SpiralSpringGenerator),
    /// Parametric geometric primitive.
    Primitive(PrimitiveShape),
}

/// One part in the unfrozen canonical assembly draft.
#[derive(Debug, Clone, PartialEq)]
pub struct DraftPartDecl {
    /// Human-readable unique part name within the draft.
    pub name: String,
    /// Persistent catalog entity identity (`EntityKind::Part`).
    pub entity: EntityId,
    /// Physical part class.
    pub kind: DraftPartKind,
    /// Geometry definition (standard or mechanical/primitive generator).
    pub geometry: DraftGeometrySource,
    /// Material card references with explicit query points and unknowns.
    pub materials: Vec<MaterialRef>,
    /// Active physics domains owned by this part (non-empty).
    pub physics: BTreeSet<PhysicsDomain>,
    /// Fidelity ladder and error budget for this part.
    pub fidelity: FidelityPolicy,
    /// Coordinate/dimensional uncertainty bounds.
    pub uncertainty: Vec<DimensionBound>,
    /// Declared initial state. A constructor that needs one refuses a part
    /// without it; none is assumed.
    pub initial: Option<DraftInitialState>,
}

/// Where a part starts, how it moves and how warm it is at the start of the
/// time plan, in the assembly frame.
#[derive(Debug, Clone, PartialEq)]
pub struct DraftInitialState {
    /// Position of the part's body origin [m].
    pub position: [QtyAny; 3],
    /// Orientation of the part's body frame, a unit quaternion `[w, x, y, z]`.
    pub orientation: [f64; 4],
    /// Linear velocity of the body origin [m/s].
    pub linear_velocity: [QtyAny; 3],
    /// Angular velocity [rad/s].
    pub angular_velocity: [QtyAny; 3],
    /// Uniform initial temperature [K]; `None` means the environment's ambient
    /// temperature.
    pub temperature: Option<QtyAny>,
}

/// Kinematic joint classes matching `fs-mbd::articulated::JointKind`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum JointKind {
    /// Zero relative degrees of freedom.
    Fixed,
    /// Single rotational degree of freedom about `axis`.
    Revolute,
    /// Single translational degree of freedom along `axis`.
    Prismatic,
    /// Coupled rotation and axial translation along `axis` with `pitch`.
    Helical,
}

/// Optional travel limits on a 1-DOF joint.
#[derive(Debug, Clone, PartialEq)]
pub struct JointLimits {
    /// Lower travel limit (`Dims::NONE` [rad] for revolute/helical, `LENGTH` [m] for prismatic).
    pub lower: QtyAny,
    /// Upper travel limit (must satisfy `lower <= upper`).
    pub upper: QtyAny,
}

/// Kinematic joint joining two parts (or a part to ground when `parent` is the root assembly).
#[derive(Debug, Clone, PartialEq)]
pub struct JointDecl {
    /// Unique joint name.
    pub name: String,
    /// Parent entity (assembly root for ground, or a part).
    pub parent: EntityId,
    /// Child part entity (must differ from `parent`).
    pub child: EntityId,
    /// Frame in which `axis` is expressed.
    pub frame: FrameId,
    /// Joint class (`Fixed`, `Revolute`, `Prismatic`, `Helical`).
    pub kind: JointKind,
    /// Joint axis vector in `frame` (must be finite and nonzero for non-fixed joints).
    pub axis: [f64; 3],
    /// Lead / pitch [m/rad] (required positive/nonzero finite for `Helical`, forbidden otherwise).
    pub pitch: Option<QtyAny>,
    /// Optional joint travel limits (forbidden on `Fixed`).
    pub limits: Option<JointLimits>,
    /// Optional viscous joint damping parameter (any finite nonnegative quantity).
    pub damping: Option<QtyAny>,
    /// Optional elastic joint stiffness parameter (any finite nonnegative quantity).
    pub stiffness: Option<QtyAny>,
}

/// Gear-mesh kinematic and compliance relation between two gears.
#[derive(Debug, Clone, PartialEq)]
pub struct GearMeshDecl {
    /// Unique gear-mesh name.
    pub name: String,
    /// Driving gear part entity.
    pub driver: EntityId,
    /// Driven gear part entity (must differ from `driver`).
    pub driven: EntityId,
    /// Driving gear tooth count (`>= 3`).
    pub teeth_driver: u32,
    /// Driven gear tooth count (`>= 3`).
    pub teeth_driven: u32,
    /// Operating pressure angle [rad] in `(0, pi/4]`.
    pub pressure_angle: QtyAny,
    /// Operating center distance [m] (positive finite).
    pub center_distance: QtyAny,
    /// Optional circumferential backlash [m] (nonnegative finite).
    pub backlash: Option<QtyAny>,
    /// Content identity of the mesh contact/efficiency card.
    pub mesh_card: ContentHash,
}

/// Third-body intervening medium (lubricant, thermal interface material, oxide layer).
#[derive(Debug, Clone, PartialEq)]
pub struct ThirdBodyDecl {
    /// Material card identity for the third body.
    pub material_card: ContentHash,
    /// Nominal film / layer thickness [m] (positive finite).
    pub film_thickness: QtyAny,
    /// Thermodynamic / shear query point for the third-body card.
    pub query_point: Vec<(ParameterKey, QtyAny)>,
}

/// Card-bound interface laws across mechanical, thermal, electrical, and power-port boundaries,
/// and card-less heat routes and shear films.
#[derive(Debug, Clone, PartialEq)]
pub enum DraftInterfaceKind {
    /// Typed power-conjugate port schema.
    PowerPort {
        /// Schema content hash.
        schema: ContentHash,
    },
    /// Mechanical / tribological contact with optional third-body lubricant and initial state.
    MechanicalContact {
        /// Contact/friction law card hash.
        law_card: ContentHash,
        /// Optional lubricant / third-body film.
        third_body: Option<ThirdBodyDecl>,
        /// Explicit interface state variables (e.g. slip, wear depth, preload).
        state: Vec<(ParameterKey, QtyAny)>,
    },
    /// Bonded / adhesive interface.
    Glue {
        /// Adhesive law card hash.
        law_card: ContentHash,
    },
    /// Thermal interface conductance / radiation / contact resistance.
    ThermalContact {
        /// Thermal interface card hash.
        law_card: ContentHash,
        /// Optional thermal interface material (TIM) layer.
        third_body: Option<ThirdBodyDecl>,
        /// Explicit interface state variables (e.g. contact pressure, surface roughness).
        state: Vec<(ParameterKey, QtyAny)>,
    },
    /// Electrical contact resistance / constriction / switching interface.
    ElectricalContact {
        /// Electrical contact card hash.
        law_card: ContentHash,
        /// Optional oxide / contact-film third body.
        third_body: Option<ThirdBodyDecl>,
        /// Explicit interface state variables (e.g. film thickness, contact force).
        state: Vec<(ParameterKey, QtyAny)>,
    },
    /// Heat bookkeeping, not a material interface, so no card: a declared
    /// share of the heat dissipated at `feature_a` delivered to `feature_b`.
    HeatRoute {
        /// Explicit route state (e.g. the dimensionless `fraction`).
        state: Vec<(ParameterKey, QtyAny)>,
    },
    /// A Newtonian shear film; no card, since its one measured property is
    /// the viscosity on its third body's card.
    ShearFilm {
        /// The fluid, its film thickness and its query point.
        third_body: ThirdBodyDecl,
        /// Explicit film state (the wetted `area`).
        state: Vec<(ParameterKey, QtyAny)>,
    },
}

/// Interface binding two surface/region features through a card-backed law.
#[derive(Debug, Clone, PartialEq)]
pub struct DraftInterfaceDecl {
    /// Unique interface name.
    pub name: String,
    /// Catalog interface entity identity (`EntityKind::Interface`).
    pub entity: EntityId,
    /// First surface/region/part feature identity.
    pub feature_a: EntityId,
    /// Second surface/region/part feature identity (must differ from `feature_a`).
    pub feature_b: EntityId,
    /// Interface law kind, card identity, state, and optional third body.
    pub kind: DraftInterfaceKind,
}

/// Target of an external source excitation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SourceTarget {
    /// Named port on a catalog entity.
    Port {
        /// Owning part or interface entity.
        entity: EntityId,
        /// Stable port name.
        port: String,
    },
    /// Dimensioned parameter path on an entity.
    Parameter(ParameterPath),
}

/// Prescribed time history for a source excitation.
#[derive(Debug, Clone, PartialEq)]
pub enum TimeHistory {
    /// Constant dimensioned value.
    Constant(QtyAny),
    /// Sinusoidal excitation `offset + amplitude * sin(2*pi*frequency*t + phase)`.
    Harmonic {
        /// Amplitude (defines output dimensions).
        amplitude: QtyAny,
        /// Frequency [Hz] (`FREQUENCY`, positive finite).
        frequency: QtyAny,
        /// Phase [rad] (`Dims::NONE`, finite).
        phase: QtyAny,
        /// DC offset (same dimensions as `amplitude`).
        offset: QtyAny,
    },
    /// Piecewise-linear time table with strictly increasing time knots [s].
    PiecewiseLinear {
        /// Ordered `(time [s], value)` pairs (at least 2 knots).
        knots: Vec<(QtyAny, QtyAny)>,
    },
    /// Step transition at `t_switch` [s].
    Step {
        /// Switching time [s] (`TIME`, finite >= 0).
        t_switch: QtyAny,
        /// Value before `t_switch`.
        before: QtyAny,
        /// Value at and after `t_switch` (same dims as `before`).
        after: QtyAny,
    },
}

/// External source acting on a port or parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceDecl {
    /// Unique source name.
    pub name: String,
    /// Target port or parameter path.
    pub target: SourceTarget,
    /// Prescribed time history.
    pub history: TimeHistory,
}

/// Feedback control law acting on signal ports.
#[derive(Debug, Clone, PartialEq)]
pub enum ControllerLaw {
    /// Proportional-integral-derivative feedback controller.
    Pid {
        /// Proportional gain.
        kp: QtyAny,
        /// Integral gain.
        ki: QtyAny,
        /// Derivative gain.
        kd: QtyAny,
        /// Target setpoint.
        setpoint: QtyAny,
    },
    /// 1D lookup schedule `(input, output)` with strictly increasing inputs.
    LookupTable {
        /// Ordered breakpoints (at least 2).
        breakpoints: Vec<(QtyAny, QtyAny)>,
    },
    /// Card-backed control law identified by content hash.
    CardBacked {
        /// Control law card identity.
        law_card: ContentHash,
        /// Dimensioned controller parameters.
        params: Vec<(ParameterKey, QtyAny)>,
    },
}

/// Closed-loop controller wired between signal ports.
#[derive(Debug, Clone, PartialEq)]
pub struct ControllerDecl {
    /// Unique controller name.
    pub name: String,
    /// Sensor signal inputs `(entity, port_name)` (at least one).
    pub sensor_ports: Vec<(EntityId, String)>,
    /// Actuator signal output `(entity, port_name)`.
    pub actuator_port: (EntityId, String),
    /// Discrete controller update rate [Hz] (`FREQUENCY`, positive finite).
    pub sample_rate: QtyAny,
    /// Control law specification.
    pub law: ControllerLaw,
}

/// Time-stepping policy for a coupled simulation run.
#[derive(Debug, Clone, PartialEq)]
pub enum StepPolicy {
    /// Fixed step size [s] (`TIME`, positive finite, `<= horizon`).
    Fixed {
        /// Step duration [s].
        dt: QtyAny,
    },
    /// Adaptive error-controlled stepping.
    Adaptive {
        /// Initial step [s].
        dt_initial: QtyAny,
        /// Minimum allowed step [s].
        dt_min: QtyAny,
        /// Maximum allowed step [s].
        dt_max: QtyAny,
        /// Relative error tolerance in `(0, 1)`.
        relative_tol: f64,
        /// Declared time order (`>= 1`).
        order: u8,
    },
}

/// How the thermal timeline integrates in time. The mechanical networks
/// have one step each and run beside any scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TimeScheme {
    /// The timeline's explicit step.
    #[default]
    Explicit,
    /// Backward Euler, on the owner's implicit conduction graph.
    BackwardEuler,
    /// Two-stage singly diagonally implicit Runge-Kutta, on the same graph.
    Sdirk2,
}

/// Simulation time horizon, stepping policy, and default output sample rate.
#[derive(Debug, Clone, PartialEq)]
pub struct TimePlan {
    /// Total simulated duration [s] (`TIME`, positive finite).
    pub horizon: QtyAny,
    /// Time-step control policy.
    pub step_policy: StepPolicy,
    /// Default output sampling frequency [Hz] (`FREQUENCY`, positive finite).
    pub default_sample_rate: QtyAny,
    /// The thermal timeline's integration scheme.
    pub scheme: TimeScheme,
}

/// One species mole fraction in the ambient gas mixture.
#[derive(Debug, Clone, PartialEq)]
pub struct GasComponent {
    /// Canonical species name (e.g. `"N2"`, `"O2"`, `"H2O"`).
    pub species: String,
    /// Mole fraction in `(0, 1]`.
    pub mole_fraction: f64,
}

/// Extended ambient environment adding applied magnetic field, humidity, and
/// explicit gas composition atop [`Environment`].
#[derive(Debug, Clone, PartialEq)]
pub struct DraftEnvironment {
    /// Base gravity, ambient temperature, and ambient pressure.
    pub base: Environment,
    /// Uniform applied magnetic flux density `[Bx, By, Bz]` [T].
    pub magnetic_field: [QtyAny; 3],
    /// Relative humidity in `[0.0, 1.0]`.
    pub relative_humidity: f64,
    /// Ambient gas mole fractions (when non-empty, must sum to `1.0` within `1e-9`).
    pub gas_composition: Vec<GasComponent>,
}

impl DraftEnvironment {
    /// Standard Earth lab environment with zero magnetic field, 45% RH, and dry air composition.
    #[must_use]
    pub fn earth_lab_air() -> Self {
        Self {
            base: Environment::earth_lab(),
            magnetic_field: [
                QtyAny::new(0.0, MAGNETIC_FIELD_TESLA),
                QtyAny::new(0.0, MAGNETIC_FIELD_TESLA),
                QtyAny::new(0.0, MAGNETIC_FIELD_TESLA),
            ],
            relative_humidity: 0.45,
            gas_composition: vec![
                GasComponent {
                    species: "N2".to_owned(),
                    mole_fraction: 0.78,
                },
                GasComponent {
                    species: "O2".to_owned(),
                    mole_fraction: 0.21,
                },
                GasComponent {
                    species: "Ar".to_owned(),
                    mole_fraction: 0.01,
                },
            ],
        }
    }
}

/// Directional pickup pattern for an acoustic microphone receiver.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MicrophonePattern {
    /// Omnidirectional pressure pickup.
    Omnidirectional,
    /// Cardioid first-order directional pickup.
    Cardioid,
    /// Bidirectional figure-eight pickup.
    FigureEight,
    /// Supercardioid pickup.
    Supercardioid,
}

/// Crossing edge polarity for a timing probe (e.g. balance-wheel beat counter).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TimingEdge {
    /// Rising crossing.
    Rising,
    /// Falling crossing.
    Falling,
    /// Both rising and falling crossings (e.g. watch beat ticks).
    Both,
}

/// Receiver / transducer model for a physical probe.
#[derive(Debug, Clone, PartialEq)]
pub enum ProbeKind {
    /// Acoustic microphone receiver with position, aim vector, and polar pattern.
    Microphone {
        /// Receiver position [m].
        position: [QtyAny; 3],
        /// Aim unit/direction vector (must be finite and nonzero).
        aim_direction: [f64; 3],
        /// Polar response pattern.
        pattern: MicrophonePattern,
    },
    /// Contact or immersed thermocouple temperature probe.
    Thermocouple {
        /// Instrumented part or surface entity.
        location: EntityId,
        /// Mounting position [m] in `frame`.
        position: [QtyAny; 3],
        /// Optional thermal lag time constant [s] (`TIME`, positive finite).
        time_constant: Option<QtyAny>,
    },
    /// Event-timing probe (escapement beat rate, impact interval, zero crossing).
    Timing {
        /// Observed part or joint entity.
        target: EntityId,
        /// Threshold quantity triggering a timestamped event.
        threshold: QtyAny,
        /// Triggering crossing polarity.
        edge: TimingEdge,
    },
    /// Surface strain gauge along a declared axis.
    StrainGauge {
        /// Instrumented part or surface entity.
        target: EntityId,
        /// Gauge mounting position [m] in `frame`.
        position: [QtyAny; 3],
        /// Sensitive strain direction (finite nonzero).
        gauge_direction: [f64; 3],
    },
    /// Generic dimensioned field probe.
    FieldQuantity {
        /// Target part/region entity.
        target: EntityId,
        /// Observation position [m] in `frame`.
        position: [QtyAny; 3],
        /// Observed physical quantity and dimensions.
        quantity: ParameterKey,
    },
}

/// Declared probe with receiver model, frame, and explicit sample rate.
#[derive(Debug, Clone, PartialEq)]
pub struct ProbeDecl {
    /// Unique probe name.
    pub name: String,
    /// Coordinate frame for positions/directions.
    pub frame: FrameId,
    /// Sampling rate [Hz] (`FREQUENCY`, positive finite).
    pub sample_rate: QtyAny,
    /// Receiver model and target.
    pub kind: ProbeKind,
}

/// Complete untrusted authoring input for an unfrozen canonical assembly draft.
#[derive(Debug, Clone, PartialEq)]
pub struct DraftAssemblyInput {
    /// Device / assembly name.
    pub name: String,
    /// Root assembly entity in `catalog`.
    pub root: EntityId,
    /// Persistent entity catalog.
    pub catalog: EntityCatalog,
    /// Explicit geometric reference/load state (never defaulted).
    pub reference_state: Option<ReferenceState>,
    /// Extended ambient environment.
    pub environment: DraftEnvironment,
    /// Simulation horizon, step policy, and default sample rate.
    pub time_plan: TimePlan,
    /// Physical parts with physics domains and fidelity policies.
    pub parts: Vec<DraftPartDecl>,
    /// Kinematic joints.
    pub joints: Vec<JointDecl>,
    /// Gear-mesh relations.
    pub gear_meshes: Vec<GearMeshDecl>,
    /// Card-bound interfaces across domains.
    pub interfaces: Vec<DraftInterfaceDecl>,
    /// External sources on ports or parameters.
    pub sources: Vec<SourceDecl>,
    /// Closed-loop signal-port controllers.
    pub controllers: Vec<ControllerDecl>,
    /// Physical probes and receiver models.
    pub probes: Vec<ProbeDecl>,
}

/// Structured named refusal from [`DraftAssembly::admit`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DraftAssemblyError {
    /// Named element or section that failed validation (e.g. `"joint:pivot"`).
    pub item: String,
    /// Human-readable reason for refusal.
    pub reason: String,
}

impl DraftAssemblyError {
    fn new(item: impl Into<String>, reason: impl Into<String>) -> Self {
        Self {
            item: item.into(),
            reason: reason.into(),
        }
    }
}

impl fmt::Display for DraftAssemblyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "draft assembly refusal at {}: {}",
            self.item, self.reason
        )
    }
}

impl std::error::Error for DraftAssemblyError {}

/// Structurally admitted unfrozen canonical assembly draft bound to a [`Scenario`].
#[derive(Debug, Clone, PartialEq)]
pub struct DraftAssembly {
    scenario_ir_hash: ContentHash,
    input: DraftAssemblyInput,
}

fn check_name(name: &str, scope: &str) -> Result<(), DraftAssemblyError> {
    if name.trim().is_empty() || name.len() > 256 || name.chars().any(char::is_control) {
        return Err(DraftAssemblyError::new(
            scope,
            "name must be non-empty, <= 256 bytes, and free of control characters",
        ));
    }
    Ok(())
}

fn check_initial_state(initial: &DraftInitialState, scope: &str) -> Result<(), DraftAssemblyError> {
    let scope = format!("{scope}.initial");
    for (field, values, dims) in [
        ("position", &initial.position, LENGTH),
        ("linear_velocity", &initial.linear_velocity, VELOCITY),
        ("angular_velocity", &initial.angular_velocity, FREQUENCY),
    ] {
        for (axis, value) in values.iter().enumerate() {
            let item = format!("{scope}.{field}[{axis}]");
            check_qty(*value, dims, &item, false, false)?;
        }
    }
    let q = initial.orientation;
    let norm_sq = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if !q.iter().all(|c| c.is_finite()) || (norm_sq - 1.0).abs() > 1.0e-12 {
        return Err(DraftAssemblyError::new(
            format!("{scope}.orientation"),
            "orientation must be a finite unit quaternion [w, x, y, z]",
        ));
    }
    if let Some(temperature) = initial.temperature {
        let item = format!("{scope}.temperature");
        check_qty(temperature, TEMPERATURE, &item, true, false)?;
    }
    Ok(())
}

fn check_qty(
    q: QtyAny,
    expected_dims: Dims,
    item: &str,
    positive_only: bool,
    nonnegative_only: bool,
) -> Result<(), DraftAssemblyError> {
    if !q.value.is_finite() {
        return Err(DraftAssemblyError::new(item, "quantity must be finite"));
    }
    if q.dims != expected_dims {
        return Err(DraftAssemblyError::new(
            item,
            format!(
                "dimension mismatch: expected {:?}, got {:?}",
                expected_dims.0, q.dims.0
            ),
        ));
    }
    if positive_only && !(q.value > 0.0) {
        return Err(DraftAssemblyError::new(
            item,
            format!("quantity must be strictly positive, got {}", q.value),
        ));
    }
    if nonnegative_only && !(q.value >= 0.0) {
        return Err(DraftAssemblyError::new(
            item,
            format!("quantity must be nonnegative, got {}", q.value),
        ));
    }
    Ok(())
}

fn check_vec3_nonzero(v: [f64; 3], item: &str) -> Result<(), DraftAssemblyError> {
    if v.iter().any(|c| !c.is_finite()) {
        return Err(DraftAssemblyError::new(
            item,
            "direction vector must be finite",
        ));
    }
    let norm2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if !(norm2.is_finite() && norm2 > 0.0) {
        return Err(DraftAssemblyError::new(
            item,
            "direction vector must be nonzero",
        ));
    }
    Ok(())
}

fn check_hash(h: ContentHash, item: &str) -> Result<(), DraftAssemblyError> {
    if h.0 == [0_u8; 32] {
        return Err(DraftAssemblyError::new(
            item,
            "content hash identity must be nonzero",
        ));
    }
    Ok(())
}

fn check_entity(
    catalog: &EntityCatalog,
    id: EntityId,
    expected: Option<EntityKind>,
    root: EntityId,
    item: &str,
) -> Result<(), DraftAssemblyError> {
    if !catalog.is_active(id) {
        return Err(DraftAssemblyError::new(
            item,
            format!("entity {} is not active in catalog", id.token()),
        ));
    }
    if let Some(kind) = expected
        && id.kind() != kind
    {
        return Err(DraftAssemblyError::new(
            item,
            format!(
                "entity {} has kind {:?}, expected {:?}",
                id.token(),
                id.kind(),
                kind
            ),
        ));
    }
    let mut cur = Some(id);
    for _ in 0..=catalog.entities().len() {
        if cur == Some(root) {
            return Ok(());
        }
        cur = cur.and_then(|e| catalog.get(e)).and_then(Entity::parent);
        if cur.is_none() {
            break;
        }
    }
    Err(DraftAssemblyError::new(
        item,
        format!("entity {} is not contained under root assembly", id.token()),
    ))
}

fn check_frame(scenario: &Scenario, frame: FrameId, item: &str) -> Result<(), DraftAssemblyError> {
    if frame.0 != 0 && !scenario.frames.frames.iter().any(|f| f.id == frame) {
        return Err(DraftAssemblyError::new(
            item,
            format!("frame {} is not declared in scenario", frame.0),
        ));
    }
    Ok(())
}

fn check_params(params: &[(ParameterKey, QtyAny)], item: &str) -> Result<(), DraftAssemblyError> {
    let mut seen = BTreeSet::new();
    for (k, q) in params {
        check_name(&k.name, item)?;
        if !seen.insert(&k.name) {
            return Err(DraftAssemblyError::new(
                item,
                format!("duplicate parameter '{}'", k.name),
            ));
        }
        check_qty(*q, k.dims, &format!("{item}.{}", k.name), false, false)?;
    }
    Ok(())
}

fn check_third_body(tb: &ThirdBodyDecl, item: &str) -> Result<(), DraftAssemblyError> {
    check_hash(
        tb.material_card,
        &format!("{item}.third_body.material_card"),
    )?;
    check_qty(
        tb.film_thickness,
        LENGTH,
        &format!("{item}.third_body.film_thickness"),
        true,
        false,
    )?;
    check_params(&tb.query_point, &format!("{item}.third_body.query_point"))
}

#[allow(clippy::too_many_lines)]
fn check_geometry(
    geom: &DraftGeometrySource,
    part_entity: EntityId,
    catalog: &EntityCatalog,
    item: &str,
) -> Result<(), DraftAssemblyError> {
    match geom {
        DraftGeometrySource::Standard(src) => match src {
            GeometrySource::Parametric { generator, params } => {
                check_name(generator, &format!("{item}.geometry.generator"))?;
                check_params(params, &format!("{item}.geometry.params"))?;
                super::generators::validate_source(src, 65_536).map_err(|e| {
                    DraftAssemblyError::new(format!("{item}.geometry.generator"), e.to_string())
                })?;
            }
            GeometrySource::Imported {
                fingerprint,
                ledger_row,
                ..
            } => {
                check_hash(fingerprint.hash(), &format!("{item}.geometry.fingerprint"))?;
                check_name(ledger_row, &format!("{item}.geometry.ledger_row"))?;
                if catalog.get(part_entity).and_then(Entity::fingerprint) != Some(*fingerprint) {
                    return Err(DraftAssemblyError::new(
                        item,
                        "imported geometry fingerprint does not match catalog entity",
                    ));
                }
            }
        },
        DraftGeometrySource::Gear(g) => {
            if g.teeth < 3 {
                return Err(DraftAssemblyError::new(
                    format!("{item}.gear.teeth"),
                    format!("gear requires at least 3 teeth, got {}", g.teeth),
                ));
            }
            check_qty(
                g.module,
                LENGTH,
                &format!("{item}.gear.module"),
                true,
                false,
            )?;
            check_qty(
                g.face_width,
                LENGTH,
                &format!("{item}.gear.face_width"),
                true,
                false,
            )?;
            check_qty(
                g.pressure_angle,
                Dims::NONE,
                &format!("{item}.gear.pressure_angle"),
                true,
                false,
            )?;
            if g.pressure_angle.value > core::f64::consts::FRAC_PI_4 {
                return Err(DraftAssemblyError::new(
                    format!("{item}.gear.pressure_angle"),
                    "pressure angle must be <= pi/4 rad",
                ));
            }
            if !g.profile_shift.is_finite() || !(-1.0..=1.0).contains(&g.profile_shift) {
                return Err(DraftAssemblyError::new(
                    format!("{item}.gear.profile_shift"),
                    "profile shift must be finite in [-1, 1]",
                ));
            }
            check_qty(
                g.bore_radius,
                LENGTH,
                &format!("{item}.gear.bore_radius"),
                false,
                true,
            )?;
            let root_r = g.root_radius_m();
            if !(root_r.is_finite() && root_r > 0.0 && g.bore_radius.value < root_r) {
                return Err(DraftAssemblyError::new(
                    format!("{item}.gear.bore_radius"),
                    format!(
                        "bore radius {} m must be strictly smaller than root radius {} m",
                        g.bore_radius.value, root_r
                    ),
                ));
            }
        }
        DraftGeometrySource::SpiralSpring(s) => {
            check_qty(
                s.inner_radius,
                LENGTH,
                &format!("{item}.spiral_spring.inner_radius"),
                true,
                false,
            )?;
            check_qty(
                s.outer_radius,
                LENGTH,
                &format!("{item}.spiral_spring.outer_radius"),
                true,
                false,
            )?;
            if s.outer_radius.value <= s.inner_radius.value {
                return Err(DraftAssemblyError::new(
                    format!("{item}.spiral_spring.outer_radius"),
                    "outer_radius must exceed inner_radius",
                ));
            }
            if !(s.turns.is_finite() && s.turns > 0.0) {
                return Err(DraftAssemblyError::new(
                    format!("{item}.spiral_spring.turns"),
                    "turns must be finite and positive",
                ));
            }
            check_qty(
                s.blade_thickness,
                LENGTH,
                &format!("{item}.spiral_spring.blade_thickness"),
                true,
                false,
            )?;
            check_qty(
                s.blade_height,
                LENGTH,
                &format!("{item}.spiral_spring.blade_height"),
                true,
                false,
            )?;
            let radial_span = s.outer_radius.value - s.inner_radius.value;
            let stacked = s.turns * s.blade_thickness.value;
            if !(stacked < radial_span) {
                return Err(DraftAssemblyError::new(
                    format!("{item}.spiral_spring.blade_thickness"),
                    format!(
                        "coils overlap: turns * blade_thickness ({stacked} m) >= radial span ({radial_span} m)"
                    ),
                ));
            }
        }
        DraftGeometrySource::Primitive(p) => match p {
            PrimitiveShape::Box { lx, ly, lz } => {
                check_qty(*lx, LENGTH, &format!("{item}.primitive.lx"), true, false)?;
                check_qty(*ly, LENGTH, &format!("{item}.primitive.ly"), true, false)?;
                check_qty(*lz, LENGTH, &format!("{item}.primitive.lz"), true, false)?;
            }
            PrimitiveShape::Cylinder { radius, height } => {
                check_qty(
                    *radius,
                    LENGTH,
                    &format!("{item}.primitive.radius"),
                    true,
                    false,
                )?;
                check_qty(
                    *height,
                    LENGTH,
                    &format!("{item}.primitive.height"),
                    true,
                    false,
                )?;
            }
            PrimitiveShape::Sphere { radius } => {
                check_qty(
                    *radius,
                    LENGTH,
                    &format!("{item}.primitive.radius"),
                    true,
                    false,
                )?;
            }
            PrimitiveShape::Tube {
                inner_radius,
                outer_radius,
                length,
            } => {
                check_qty(
                    *inner_radius,
                    LENGTH,
                    &format!("{item}.primitive.inner_radius"),
                    true,
                    false,
                )?;
                check_qty(
                    *outer_radius,
                    LENGTH,
                    &format!("{item}.primitive.outer_radius"),
                    true,
                    false,
                )?;
                if outer_radius.value <= inner_radius.value {
                    return Err(DraftAssemblyError::new(
                        format!("{item}.primitive.outer_radius"),
                        "tube outer_radius must exceed inner_radius",
                    ));
                }
                check_qty(
                    *length,
                    LENGTH,
                    &format!("{item}.primitive.length"),
                    true,
                    false,
                )?;
            }
        },
    }
    Ok(())
}

impl DraftAssembly {
    /// Structurally admit an unfrozen canonical assembly draft bound to `scenario`.
    ///
    /// # Errors
    /// Returns [`DraftAssemblyError`] naming the exact offending element and reason.
    #[allow(clippy::too_many_lines)]
    pub fn admit(
        input: DraftAssemblyInput,
        scenario: &Scenario,
    ) -> Result<Self, DraftAssemblyError> {
        check_name(&input.name, "assembly.name")?;
        if input.reference_state.is_none() {
            return Err(DraftAssemblyError::new(
                "assembly.reference_state",
                "reference_state must be explicitly declared (never defaulted)",
            ));
        }
        check_entity(
            &input.catalog,
            input.root,
            Some(EntityKind::Assembly),
            input.root,
            "assembly.root",
        )?;
        if !input.catalog.verify_receipts() {
            return Err(DraftAssemblyError::new(
                "assembly.catalog",
                "entity catalog receipt chain failed verification",
            ));
        }
        let findings = input.catalog.validate();
        if !findings.is_empty() {
            return Err(DraftAssemblyError::new(
                "assembly.catalog",
                format!("catalog validation findings: {findings:?}"),
            ));
        }
        for placement in input.catalog.placements() {
            check_frame(scenario, placement.frame(), "assembly.catalog.placements")?;
        }

        // 1. Environment validation (both base gravity/T/p and extended fields).
        for (i, g) in input.environment.base.gravity.iter().enumerate() {
            check_qty(
                *g,
                ACCEL,
                &format!("environment.base.gravity[{i}]"),
                false,
                false,
            )?;
        }
        check_qty(
            input.environment.base.ambient_temperature,
            TEMPERATURE,
            "environment.base.ambient_temperature",
            true,
            false,
        )?;
        check_qty(
            input.environment.base.ambient_pressure,
            PRESSURE,
            "environment.base.ambient_pressure",
            false,
            true,
        )?;
        for (i, b) in input.environment.magnetic_field.iter().enumerate() {
            check_qty(
                *b,
                MAGNETIC_FIELD_TESLA,
                &format!("environment.magnetic_field[{i}]"),
                false,
                false,
            )?;
        }
        if !input.environment.relative_humidity.is_finite()
            || !(0.0..=1.0).contains(&input.environment.relative_humidity)
        {
            return Err(DraftAssemblyError::new(
                "environment.relative_humidity",
                format!(
                    "relative humidity must be in [0, 1], got {}",
                    input.environment.relative_humidity
                ),
            ));
        }
        if !input.environment.gas_composition.is_empty() {
            let mut species_seen = BTreeSet::new();
            let mut sum = 0.0_f64;
            for comp in &input.environment.gas_composition {
                check_name(&comp.species, "environment.gas_composition")?;
                if !species_seen.insert(&comp.species) {
                    return Err(DraftAssemblyError::new(
                        format!("environment.gas_composition:{}", comp.species),
                        "duplicate gas species",
                    ));
                }
                if !(comp.mole_fraction.is_finite()
                    && comp.mole_fraction > 0.0
                    && comp.mole_fraction <= 1.0)
                {
                    return Err(DraftAssemblyError::new(
                        format!("environment.gas_composition:{}", comp.species),
                        format!(
                            "mole_fraction must be in (0, 1], got {}",
                            comp.mole_fraction
                        ),
                    ));
                }
                sum += comp.mole_fraction;
            }
            if (sum - 1.0).abs() > 1.0e-9 {
                return Err(DraftAssemblyError::new(
                    "environment.gas_composition",
                    format!("gas mole fractions must sum to 1.0, got {sum}"),
                ));
            }
        }

        // 2. Time plan validation.
        check_qty(
            input.time_plan.horizon,
            TIME,
            "time_plan.horizon",
            true,
            false,
        )?;
        check_qty(
            input.time_plan.default_sample_rate,
            FREQUENCY,
            "time_plan.default_sample_rate",
            true,
            false,
        )?;
        match &input.time_plan.step_policy {
            StepPolicy::Fixed { dt } => {
                check_qty(*dt, TIME, "time_plan.step_policy.dt", true, false)?;
                if dt.value > input.time_plan.horizon.value {
                    return Err(DraftAssemblyError::new(
                        "time_plan.step_policy.dt",
                        "fixed dt must not exceed horizon",
                    ));
                }
            }
            StepPolicy::Adaptive {
                dt_initial,
                dt_min,
                dt_max,
                relative_tol,
                order,
            } => {
                check_qty(*dt_min, TIME, "time_plan.step_policy.dt_min", true, false)?;
                check_qty(
                    *dt_initial,
                    TIME,
                    "time_plan.step_policy.dt_initial",
                    true,
                    false,
                )?;
                check_qty(*dt_max, TIME, "time_plan.step_policy.dt_max", true, false)?;
                if !(dt_min.value <= dt_initial.value
                    && dt_initial.value <= dt_max.value
                    && dt_max.value <= input.time_plan.horizon.value)
                {
                    return Err(DraftAssemblyError::new(
                        "time_plan.step_policy",
                        "adaptive step bounds must satisfy dt_min <= dt_initial <= dt_max <= horizon",
                    ));
                }
                if !(relative_tol.is_finite() && *relative_tol > 0.0 && *relative_tol < 1.0) {
                    return Err(DraftAssemblyError::new(
                        "time_plan.step_policy.relative_tol",
                        "relative_tol must be finite in (0, 1)",
                    ));
                }
                if *order == 0 {
                    return Err(DraftAssemblyError::new(
                        "time_plan.step_policy.order",
                        "adaptive order must be >= 1",
                    ));
                }
            }
        }

        // 3. Parts validation.
        if input.parts.is_empty() {
            return Err(DraftAssemblyError::new(
                "assembly.parts",
                "assembly must declare at least one part",
            ));
        }
        let mut part_names = BTreeSet::new();
        let mut part_entities = BTreeSet::new();
        for part in &input.parts {
            let scope = format!("part:{}", part.name);
            check_name(&part.name, &scope)?;
            if !part_names.insert(&part.name) {
                return Err(DraftAssemblyError::new(&scope, "duplicate part name"));
            }
            check_entity(
                &input.catalog,
                part.entity,
                Some(EntityKind::Part),
                input.root,
                &scope,
            )?;
            if !part_entities.insert(part.entity) {
                return Err(DraftAssemblyError::new(&scope, "duplicate part entity"));
            }
            if part.physics.is_empty() {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "part must declare at least one PhysicsDomain",
                ));
            }
            check_geometry(&part.geometry, part.entity, &input.catalog, &scope)?;
            if let Some(initial) = &part.initial {
                check_initial_state(initial, &scope)?;
            }
            if part.fidelity.rungs.is_empty() {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "fidelity policy must declare at least one rung",
                ));
            }
            for (r_idx, rung) in part.fidelity.rungs.iter().enumerate() {
                let r_scope = format!("{scope}.fidelity.rungs[{r_idx}]");
                match rung {
                    FidelityRung::Lumped => {}
                    FidelityRung::Modal { mode_count } if *mode_count > 0 => {}
                    FidelityRung::ReducedOrder { basis_dim } if *basis_dim > 0 => {}
                    FidelityRung::OneDField { cells } if *cells > 0 => {}
                    FidelityRung::BoundaryElement { panels } if *panels > 0 => {}
                    FidelityRung::FiniteElement { order, h_target } => {
                        if !(1..=6).contains(order) {
                            return Err(DraftAssemblyError::new(
                                &r_scope,
                                "FEM order must be in 1..=6",
                            ));
                        }
                        check_qty(*h_target, LENGTH, &r_scope, true, false)?;
                    }
                    FidelityRung::LatticeBoltzmann { dx } => {
                        check_qty(*dx, LENGTH, &r_scope, true, false)?;
                    }
                    _ => {
                        return Err(DraftAssemblyError::new(
                            &r_scope,
                            "fidelity rung resolution parameter must be positive",
                        ));
                    }
                }
            }
            let eb = &part.fidelity.error_budget;
            if !(eb.relative_tol.is_finite() && eb.relative_tol > 0.0 && eb.relative_tol < 1.0) {
                return Err(DraftAssemblyError::new(
                    format!("{scope}.fidelity.error_budget.relative_tol"),
                    "relative_tol must be finite in (0, 1)",
                ));
            }
            if !eb.absolute_tol.value.is_finite() || eb.absolute_tol.value < 0.0 {
                return Err(DraftAssemblyError::new(
                    format!("{scope}.fidelity.error_budget.absolute_tol"),
                    "absolute_tol must be finite and nonnegative",
                ));
            }
            for (u_idx, b) in part.uncertainty.iter().enumerate() {
                let u_scope = format!("{scope}.uncertainty[{u_idx}]");
                check_name(&b.property.name, &u_scope)?;
                check_qty(
                    b.lower,
                    b.property.dims,
                    &format!("{u_scope}.lower"),
                    false,
                    false,
                )?;
                check_qty(
                    b.upper,
                    b.property.dims,
                    &format!("{u_scope}.upper"),
                    false,
                    false,
                )?;
                if b.lower.value > b.upper.value {
                    return Err(DraftAssemblyError::new(
                        &u_scope,
                        "uncertainty bound lower must be <= upper",
                    ));
                }
            }
            for (m_idx, mat) in part.materials.iter().enumerate() {
                let m_scope = format!("{scope}.materials[{m_idx}]");
                check_hash(mat.card, &m_scope)?;
                check_frame(scenario, mat.frame, &m_scope)?;
                check_params(&mat.query_point, &m_scope)?;
                let mut unknown_names = BTreeSet::new();
                for UnknownProperty::Unknown { property, bound } in &mat.unknown {
                    check_name(&property.name, &m_scope)?;
                    if !unknown_names.insert(&property.name) {
                        return Err(DraftAssemblyError::new(
                            &m_scope,
                            format!("duplicate unknown property '{}'", property.name),
                        ));
                    }
                    if mat.query_point.iter().any(|(k, _)| k.name == property.name) {
                        return Err(DraftAssemblyError::new(
                            &m_scope,
                            format!(
                                "unknown property '{}' cannot also carry a query value",
                                property.name
                            ),
                        ));
                    }
                    if let Some(b) = bound {
                        if b.property != *property {
                            return Err(DraftAssemblyError::new(
                                &m_scope,
                                format!(
                                    "unknown property bound does not match '{}'",
                                    property.name
                                ),
                            ));
                        }
                        check_qty(
                            b.lower,
                            b.property.dims,
                            &format!("{m_scope}.bound.lower"),
                            false,
                            false,
                        )?;
                        check_qty(
                            b.upper,
                            b.property.dims,
                            &format!("{m_scope}.bound.upper"),
                            false,
                            false,
                        )?;
                        if b.lower.value > b.upper.value {
                            return Err(DraftAssemblyError::new(
                                &m_scope,
                                "unknown property bound lower must be <= upper",
                            ));
                        }
                    }
                }
            }
        }

        // 4. Joints validation.
        let mut joint_names = BTreeSet::new();
        for joint in &input.joints {
            let scope = format!("joint:{}", joint.name);
            check_name(&joint.name, &scope)?;
            if !joint_names.insert(&joint.name) {
                return Err(DraftAssemblyError::new(&scope, "duplicate joint name"));
            }
            check_entity(&input.catalog, joint.parent, None, input.root, &scope)?;
            if joint.parent != input.root && !part_entities.contains(&joint.parent) {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "joint parent must be the root assembly or a declared part in this assembly",
                ));
            }
            check_entity(
                &input.catalog,
                joint.child,
                Some(EntityKind::Part),
                input.root,
                &scope,
            )?;
            if !part_entities.contains(&joint.child) {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "joint child must be a declared part in this assembly",
                ));
            }
            if joint.parent == joint.child {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "joint parent and child must be distinct entities",
                ));
            }
            check_frame(scenario, joint.frame, &scope)?;
            match joint.kind {
                JointKind::Fixed => {
                    if joint.limits.is_some() {
                        return Err(DraftAssemblyError::new(
                            &scope,
                            "fixed joint cannot declare travel limits",
                        ));
                    }
                    if joint.pitch.is_some() {
                        return Err(DraftAssemblyError::new(
                            &scope,
                            "fixed joint cannot declare helical pitch",
                        ));
                    }
                }
                JointKind::Revolute | JointKind::Prismatic | JointKind::Helical => {
                    check_vec3_nonzero(joint.axis, &format!("{scope}.axis"))?;
                    if joint.kind == JointKind::Helical {
                        let Some(pitch) = joint.pitch else {
                            return Err(DraftAssemblyError::new(
                                &scope,
                                "helical joint requires pitch [m/rad]",
                            ));
                        };
                        check_qty(pitch, LENGTH, &format!("{scope}.pitch"), false, false)?;
                        if pitch.value == 0.0 {
                            return Err(DraftAssemblyError::new(
                                &scope,
                                "helical pitch must be nonzero",
                            ));
                        }
                    } else if joint.pitch.is_some() {
                        return Err(DraftAssemblyError::new(
                            &scope,
                            "pitch is valid only for Helical joints",
                        ));
                    }
                    if let Some(lim) = &joint.limits {
                        let expected_dims = match joint.kind {
                            JointKind::Prismatic => LENGTH,
                            JointKind::Revolute | JointKind::Helical => Dims::NONE,
                            JointKind::Fixed => unreachable!(),
                        };
                        check_qty(
                            lim.lower,
                            expected_dims,
                            &format!("{scope}.limits.lower"),
                            false,
                            false,
                        )?;
                        check_qty(
                            lim.upper,
                            expected_dims,
                            &format!("{scope}.limits.upper"),
                            false,
                            false,
                        )?;
                        if lim.lower.value > lim.upper.value {
                            return Err(DraftAssemblyError::new(
                                &scope,
                                "joint limits lower must be <= upper",
                            ));
                        }
                    }
                }
            }
            if let Some(d) = joint.damping
                && (!d.value.is_finite() || d.value < 0.0)
            {
                return Err(DraftAssemblyError::new(
                    format!("{scope}.damping"),
                    "joint damping must be finite and nonnegative",
                ));
            }
            if let Some(k) = joint.stiffness
                && (!k.value.is_finite() || k.value < 0.0)
            {
                return Err(DraftAssemblyError::new(
                    format!("{scope}.stiffness"),
                    "joint stiffness must be finite and nonnegative",
                ));
            }
        }

        // 5. Gear meshes validation.
        let mut mesh_names = BTreeSet::new();
        for mesh in &input.gear_meshes {
            let scope = format!("gear_mesh:{}", mesh.name);
            check_name(&mesh.name, &scope)?;
            if !mesh_names.insert(&mesh.name) {
                return Err(DraftAssemblyError::new(&scope, "duplicate gear_mesh name"));
            }
            check_entity(
                &input.catalog,
                mesh.driver,
                Some(EntityKind::Part),
                input.root,
                &scope,
            )?;
            check_entity(
                &input.catalog,
                mesh.driven,
                Some(EntityKind::Part),
                input.root,
                &scope,
            )?;
            if !part_entities.contains(&mesh.driver) || !part_entities.contains(&mesh.driven) {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "gear mesh driver and driven must be declared parts in this assembly",
                ));
            }
            if mesh.driver == mesh.driven {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "gear mesh driver and driven must be distinct parts",
                ));
            }
            if mesh.teeth_driver < 3 || mesh.teeth_driven < 3 {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "gear mesh requires >= 3 teeth on both driver and driven gears",
                ));
            }
            check_qty(
                mesh.pressure_angle,
                Dims::NONE,
                &format!("{scope}.pressure_angle"),
                true,
                false,
            )?;
            if mesh.pressure_angle.value > core::f64::consts::FRAC_PI_4 {
                return Err(DraftAssemblyError::new(
                    format!("{scope}.pressure_angle"),
                    "pressure angle must be <= pi/4 rad",
                ));
            }
            check_qty(
                mesh.center_distance,
                LENGTH,
                &format!("{scope}.center_distance"),
                true,
                false,
            )?;
            if let Some(b) = mesh.backlash {
                check_qty(b, LENGTH, &format!("{scope}.backlash"), false, true)?;
            }
            check_hash(mesh.mesh_card, &format!("{scope}.mesh_card"))?;
        }

        // 6. Interfaces validation.
        let mut iface_names = BTreeSet::new();
        let mut iface_entities = BTreeSet::new();
        for iface in &input.interfaces {
            let scope = format!("interface:{}", iface.name);
            check_name(&iface.name, &scope)?;
            if !iface_names.insert(&iface.name) {
                return Err(DraftAssemblyError::new(&scope, "duplicate interface name"));
            }
            check_entity(
                &input.catalog,
                iface.entity,
                Some(EntityKind::Interface),
                input.root,
                &scope,
            )?;
            if !iface_entities.insert(iface.entity) {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "duplicate interface entity",
                ));
            }
            check_entity(&input.catalog, iface.feature_a, None, input.root, &scope)?;
            check_entity(&input.catalog, iface.feature_b, None, input.root, &scope)?;
            if iface.feature_a == iface.feature_b {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "interface feature_a and feature_b must be distinct",
                ));
            }
            if let Some(pair) = input.catalog.get(iface.entity).and_then(Entity::pair) {
                let same_pair = (pair.from() == iface.feature_a && pair.to() == iface.feature_b)
                    || (pair.from() == iface.feature_b && pair.to() == iface.feature_a);
                if !same_pair {
                    return Err(DraftAssemblyError::new(
                        &scope,
                        "interface feature_a/feature_b do not match catalog interface pair endpoints",
                    ));
                }
            }
            match &iface.kind {
                DraftInterfaceKind::PowerPort { schema } => {
                    check_hash(*schema, &format!("{scope}.schema"))?;
                }
                DraftInterfaceKind::Glue { law_card } => {
                    check_hash(*law_card, &format!("{scope}.law_card"))?;
                }
                DraftInterfaceKind::MechanicalContact {
                    law_card,
                    third_body,
                    state,
                }
                | DraftInterfaceKind::ThermalContact {
                    law_card,
                    third_body,
                    state,
                }
                | DraftInterfaceKind::ElectricalContact {
                    law_card,
                    third_body,
                    state,
                } => {
                    check_hash(*law_card, &format!("{scope}.law_card"))?;
                    if let Some(tb) = third_body {
                        check_third_body(tb, &scope)?;
                    }
                    check_params(state, &format!("{scope}.state"))?;
                }
                DraftInterfaceKind::HeatRoute { state } => {
                    check_params(state, &format!("{scope}.state"))?;
                }
                DraftInterfaceKind::ShearFilm { third_body, state } => {
                    check_third_body(third_body, &scope)?;
                    check_params(state, &format!("{scope}.state"))?;
                }
            }
        }

        // 7. Sources validation.
        let mut src_names = BTreeSet::new();
        for src in &input.sources {
            let scope = format!("source:{}", src.name);
            check_name(&src.name, &scope)?;
            if !src_names.insert(&src.name) {
                return Err(DraftAssemblyError::new(&scope, "duplicate source name"));
            }
            let expected_value_dims = match &src.target {
                SourceTarget::Port { entity, port } => {
                    check_entity(&input.catalog, *entity, None, input.root, &scope)?;
                    check_name(port, &format!("{scope}.port"))?;
                    None
                }
                SourceTarget::Parameter(path) => {
                    check_entity(&input.catalog, path.entity, None, input.root, &scope)?;
                    check_name(&path.parameter.name, &format!("{scope}.parameter"))?;
                    Some(path.parameter.dims)
                }
            };
            let history_dims = match &src.history {
                TimeHistory::Constant(v) => {
                    if !v.value.is_finite() {
                        return Err(DraftAssemblyError::new(
                            &scope,
                            "constant source must be finite",
                        ));
                    }
                    v.dims
                }
                TimeHistory::Harmonic {
                    amplitude,
                    frequency,
                    phase,
                    offset,
                } => {
                    if !amplitude.value.is_finite() {
                        return Err(DraftAssemblyError::new(&scope, "amplitude must be finite"));
                    }
                    check_qty(
                        *frequency,
                        FREQUENCY,
                        &format!("{scope}.frequency"),
                        true,
                        false,
                    )?;
                    check_qty(*phase, Dims::NONE, &format!("{scope}.phase"), false, false)?;
                    check_qty(
                        *offset,
                        amplitude.dims,
                        &format!("{scope}.offset"),
                        false,
                        false,
                    )?;
                    amplitude.dims
                }
                TimeHistory::PiecewiseLinear { knots } => {
                    if knots.len() < 2 {
                        return Err(DraftAssemblyError::new(
                            &scope,
                            "piecewise-linear time history requires at least 2 knots",
                        ));
                    }
                    let val_dims = knots[0].1.dims;
                    for (idx, (t, v)) in knots.iter().enumerate() {
                        check_qty(*t, TIME, &format!("{scope}.knots[{idx}].t"), false, false)?;
                        check_qty(
                            *v,
                            val_dims,
                            &format!("{scope}.knots[{idx}].v"),
                            false,
                            false,
                        )?;
                        if idx > 0 && !(t.value > knots[idx - 1].0.value) {
                            return Err(DraftAssemblyError::new(
                                &scope,
                                "piecewise-linear knot times must be strictly increasing",
                            ));
                        }
                    }
                    val_dims
                }
                TimeHistory::Step {
                    t_switch,
                    before,
                    after,
                } => {
                    check_qty(*t_switch, TIME, &format!("{scope}.t_switch"), false, true)?;
                    if !before.value.is_finite() {
                        return Err(DraftAssemblyError::new(
                            &scope,
                            "before value must be finite",
                        ));
                    }
                    check_qty(*after, before.dims, &format!("{scope}.after"), false, false)?;
                    before.dims
                }
            };
            if let Some(expected) = expected_value_dims
                && history_dims != expected
            {
                return Err(DraftAssemblyError::new(
                    &scope,
                    format!(
                        "source history dimensions {:?} do not match target parameter dimensions {:?}",
                        history_dims.0, expected.0
                    ),
                ));
            }
        }

        // 8. Controllers validation.
        let mut ctrl_names = BTreeSet::new();
        for ctrl in &input.controllers {
            let scope = format!("controller:{}", ctrl.name);
            check_name(&ctrl.name, &scope)?;
            if !ctrl_names.insert(&ctrl.name) {
                return Err(DraftAssemblyError::new(&scope, "duplicate controller name"));
            }
            if ctrl.sensor_ports.is_empty() {
                return Err(DraftAssemblyError::new(
                    &scope,
                    "controller requires at least one sensor port",
                ));
            }
            for (entity, port) in &ctrl.sensor_ports {
                check_entity(&input.catalog, *entity, None, input.root, &scope)?;
                check_name(port, &format!("{scope}.sensor_port"))?;
            }
            check_entity(
                &input.catalog,
                ctrl.actuator_port.0,
                None,
                input.root,
                &scope,
            )?;
            check_name(&ctrl.actuator_port.1, &format!("{scope}.actuator_port"))?;
            check_qty(
                ctrl.sample_rate,
                FREQUENCY,
                &format!("{scope}.sample_rate"),
                true,
                false,
            )?;
            match &ctrl.law {
                ControllerLaw::Pid {
                    kp,
                    ki,
                    kd,
                    setpoint,
                } => {
                    for (lbl, q) in [("kp", kp), ("ki", ki), ("kd", kd), ("setpoint", setpoint)] {
                        if !q.value.is_finite() {
                            return Err(DraftAssemblyError::new(
                                format!("{scope}.{lbl}"),
                                "PID parameter must be finite",
                            ));
                        }
                    }
                }
                ControllerLaw::LookupTable { breakpoints } => {
                    if breakpoints.len() < 2 {
                        return Err(DraftAssemblyError::new(
                            &scope,
                            "lookup table requires at least 2 breakpoints",
                        ));
                    }
                    let in_dims = breakpoints[0].0.dims;
                    let out_dims = breakpoints[0].1.dims;
                    for (idx, (x, y)) in breakpoints.iter().enumerate() {
                        check_qty(
                            *x,
                            in_dims,
                            &format!("{scope}.breakpoints[{idx}].x"),
                            false,
                            false,
                        )?;
                        check_qty(
                            *y,
                            out_dims,
                            &format!("{scope}.breakpoints[{idx}].y"),
                            false,
                            false,
                        )?;
                        if idx > 0 && !(x.value > breakpoints[idx - 1].0.value) {
                            return Err(DraftAssemblyError::new(
                                &scope,
                                "lookup breakpoints must be strictly increasing in input",
                            ));
                        }
                    }
                }
                ControllerLaw::CardBacked { law_card, params } => {
                    check_hash(*law_card, &format!("{scope}.law_card"))?;
                    check_params(params, &format!("{scope}.params"))?;
                }
            }
        }

        // 9. Probes validation.
        let mut probe_names = BTreeSet::new();
        for probe in &input.probes {
            let scope = format!("probe:{}", probe.name);
            check_name(&probe.name, &scope)?;
            if !probe_names.insert(&probe.name) {
                return Err(DraftAssemblyError::new(&scope, "duplicate probe name"));
            }
            check_frame(scenario, probe.frame, &scope)?;
            check_qty(
                probe.sample_rate,
                FREQUENCY,
                &format!("{scope}.sample_rate"),
                true,
                false,
            )?;
            match &probe.kind {
                ProbeKind::Microphone {
                    position,
                    aim_direction,
                    ..
                } => {
                    for (i, p) in position.iter().enumerate() {
                        check_qty(*p, LENGTH, &format!("{scope}.position[{i}]"), false, false)?;
                    }
                    check_vec3_nonzero(*aim_direction, &format!("{scope}.aim_direction"))?;
                }
                ProbeKind::Thermocouple {
                    location,
                    position,
                    time_constant,
                } => {
                    check_entity(&input.catalog, *location, None, input.root, &scope)?;
                    for (i, p) in position.iter().enumerate() {
                        check_qty(*p, LENGTH, &format!("{scope}.position[{i}]"), false, false)?;
                    }
                    if let Some(tau) = time_constant {
                        check_qty(*tau, TIME, &format!("{scope}.time_constant"), true, false)?;
                    }
                }
                ProbeKind::Timing {
                    target, threshold, ..
                } => {
                    check_entity(&input.catalog, *target, None, input.root, &scope)?;
                    if !threshold.value.is_finite() {
                        return Err(DraftAssemblyError::new(
                            format!("{scope}.threshold"),
                            "timing threshold must be finite",
                        ));
                    }
                }
                ProbeKind::StrainGauge {
                    target,
                    position,
                    gauge_direction,
                } => {
                    check_entity(&input.catalog, *target, None, input.root, &scope)?;
                    for (i, p) in position.iter().enumerate() {
                        check_qty(*p, LENGTH, &format!("{scope}.position[{i}]"), false, false)?;
                    }
                    check_vec3_nonzero(*gauge_direction, &format!("{scope}.gauge_direction"))?;
                }
                ProbeKind::FieldQuantity {
                    target,
                    position,
                    quantity,
                } => {
                    check_entity(&input.catalog, *target, None, input.root, &scope)?;
                    for (i, p) in position.iter().enumerate() {
                        check_qty(*p, LENGTH, &format!("{scope}.position[{i}]"), false, false)?;
                    }
                    check_name(&quantity.name, &format!("{scope}.quantity"))?;
                }
            }
        }

        let scenario_ir_hash = hash_bytes(crate::ir::write_ir(scenario).as_bytes());
        Ok(Self {
            scenario_ir_hash,
            input,
        })
    }

    /// Borrow the admitted draft input.
    #[must_use]
    pub const fn input(&self) -> &DraftAssemblyInput {
        &self.input
    }

    /// Exact canonical hash of the bound [`Scenario`].
    #[must_use]
    pub const fn scenario_ir_hash(&self) -> ContentHash {
        self.scenario_ir_hash
    }

    /// Clone and re-admit the in-memory draft representation against `scenario`,
    /// proving structural round-trip invariance before schema freeze.
    ///
    /// # Errors
    /// Propagates any structural or scenario-binding refusal.
    pub fn round_trip(&self, scenario: &Scenario) -> Result<Self, DraftAssemblyError> {
        let re_admitted = Self::admit(self.input.clone(), scenario)?;
        if re_admitted.scenario_ir_hash != self.scenario_ir_hash {
            return Err(DraftAssemblyError::new(
                "assembly.scenario",
                "bound scenario IR hash changed during round-trip",
            ));
        }
        Ok(re_admitted)
    }

    /// Deterministic content digest of this admitted draft and its bound scenario.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let canonical = format!("{:?}|{:02x?}", self.input, self.scenario_ir_hash.0);
        hash_bytes(canonical.as_bytes())
    }
}
