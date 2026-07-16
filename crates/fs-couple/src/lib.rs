//! fs-couple — multiphysics composition through port-Hamiltonian Dirac
//! structures. Layer: L3.
//!
//! Ad-hoc FSI staggering suffers added-mass instabilities and energy drift.
//! The implemented DIRAC interconnection is LOSSLESS BY CONSTRUCTION:
//! power-conjugate [`Port`]s use equal effort and opposite flow, so their net
//! interface power is exactly zero. [`EnergyAudit`] records caller-supplied
//! interface balances as a G0 bug alarm. Neither invariant alone proves that
//! the coupled components, discretizations, transfers, iterations, time
//! integrators, sources, or a finite accounting window are passive.
//!
//! [`PortSchema`] v2 is the dependency-light, versioned description carried by
//! new coupling relations. It makes identity, dimensions, shape, coordinates,
//! clock, power pairing, and conservation roles explicit. The four primitive
//! relation descriptors distinguish conservative junctions, storage,
//! dissipation, and sources/reservoirs instead of smuggling all four claims
//! into a lossless topology.
//!
//! For the hard, strongly-coupled cases, [`AitkenRelaxation`] gives dynamic
//! interface relaxation: on the classic ADDED-MASS-INSTABILITY fixture (a light
//! structure in a dense fluid) naive staggering diverges, while Aitken-relaxed
//! coupling converges — demonstrated by [`iterate_fixed_relaxation`] vs
//! [`iterate_aitken`]. Deterministic; depends only on the neutral `fs-iface`
//! vocabulary and `fs-qty`'s six-base dimension vector.

use core::num::NonZeroUsize;

use fs_iface::SpaceType;
use fs_qty::{Dims, Force, Power, Pressure, Temperature, Velocity, VolumetricFlowRate};

/// Current public port-schema version.
pub const PORT_SCHEMA_VERSION: u16 = 2;

/// The physical type of a power-conjugate port (its effort/flow pair).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortKind {
    /// Mechanical: force (effort) × velocity (flow).
    MechanicalForceVelocity,
    /// Fluid: pressure (effort) × volumetric flux (flow).
    FluidPressureFlux,
    /// Thermal: temperature (effort) × entropy flow.
    ThermalTemperatureEntropy,
}

impl PortKind {
    /// The effort dimension of the three scalar seed port kinds.
    #[must_use]
    pub const fn scalar_effort_dimensions(self) -> Dims {
        match self {
            Self::MechanicalForceVelocity => Force::DIMS,
            Self::FluidPressureFlux => Pressure::DIMS,
            Self::ThermalTemperatureEntropy => Temperature::DIMS,
        }
    }

    /// The flow dimension of the three scalar seed port kinds.
    #[must_use]
    pub const fn scalar_flow_dimensions(self) -> Dims {
        match self {
            Self::MechanicalForceVelocity => Velocity::DIMS,
            Self::FluidPressureFlux => VolumetricFlowRate::DIMS,
            // Entropy flow is W/K.
            Self::ThermalTemperatureEntropy => Dims([2, 1, -3, -1, 0, 0]),
        }
    }

    /// Migrate one legacy scalar kind into an explicit v2 schema.
    ///
    /// Identity, coordinates, and time remain caller supplied: the migration
    /// must not invent any of the Five Explicits.
    ///
    /// # Errors
    ///
    /// Returns a structured schema error if the supplied metadata does not
    /// form an admissible scalar power pairing.
    pub fn scalar_seed_schema(
        self,
        id: StableId,
        coordinates: CoordinateBinding,
        timestamp: PortTimestamp,
    ) -> Result<PortSchema, CoupleError> {
        PortSchema::try_new(
            id,
            self,
            self.scalar_effort_dimensions(),
            self.scalar_flow_dimensions(),
            PortValueShape::Scalar,
            coordinates,
            PowerPairing::ScalarProduct,
            timestamp,
            [ConservationRole::Energy],
        )
    }
}

/// A validated, stable machine identifier used by port and relation schemas.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StableId(String);

impl StableId {
    /// Validate a stable identifier.
    ///
    /// The admitted alphabet is intentionally transport-safe and canonical:
    /// ASCII alphanumerics plus `-`, `_`, `.`, `:`, and `/`. The first byte
    /// must be alphanumeric.
    ///
    /// # Errors
    ///
    /// [`CoupleError::InvalidStableId`] for an empty or non-canonical value.
    pub fn new(value: impl Into<String>) -> Result<Self, CoupleError> {
        let value = value.into();
        let mut chars = value.chars();
        let valid_first = chars.next().is_some_and(|c| c.is_ascii_alphanumeric());
        let valid_tail =
            chars.all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.' | ':' | '/'));
        if !valid_first || !valid_tail {
            return Err(CoupleError::InvalidStableId { value });
        }
        Ok(Self(value))
    }

    /// Borrow the canonical identifier.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// The shape paired by a port's effort and flow coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortValueShape {
    /// One scalar effort and one scalar flow.
    Scalar,
    /// A finite vector with a statically non-zero component count.
    Vector(NonZeroUsize),
    /// A finite tensor in a declared basis.
    Tensor {
        /// Row count.
        rows: NonZeroUsize,
        /// Column count.
        columns: NonZeroUsize,
    },
    /// A field trace in a neutral function-space role.
    Field {
        /// Component count at each field point.
        components: NonZeroUsize,
        /// FEEC/interface function-space role.
        space: SpaceType,
    },
}

impl PortValueShape {
    /// Construct a non-empty vector shape.
    ///
    /// # Errors
    /// [`CoupleError::EmptyPortShape`] when `components == 0`.
    pub fn vector(components: usize) -> Result<Self, CoupleError> {
        NonZeroUsize::new(components)
            .map(Self::Vector)
            .ok_or(CoupleError::EmptyPortShape)
    }

    /// Construct a non-empty tensor shape.
    ///
    /// # Errors
    /// [`CoupleError::EmptyPortShape`] when either extent is zero.
    pub fn tensor(rows: usize, columns: usize) -> Result<Self, CoupleError> {
        let rows = NonZeroUsize::new(rows).ok_or(CoupleError::EmptyPortShape)?;
        let columns = NonZeroUsize::new(columns).ok_or(CoupleError::EmptyPortShape)?;
        Ok(Self::Tensor { rows, columns })
    }

    /// Construct a non-empty field shape.
    ///
    /// # Errors
    /// [`CoupleError::EmptyPortShape`] when `components == 0`.
    pub fn field(components: usize, space: SpaceType) -> Result<Self, CoupleError> {
        let components = NonZeroUsize::new(components).ok_or(CoupleError::EmptyPortShape)?;
        Ok(Self::Field { components, space })
    }
}

/// The positive-coordinate convention of a port.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortOrientation {
    /// Positive flow leaves the component that owns the port.
    OutwardFromOwner,
    /// Positive values follow the declared frame/basis orientation.
    AlongFrame,
    /// Positive values oppose the declared frame/basis orientation.
    AgainstFrame,
}

impl PortOrientation {
    fn composes_with(self, other: Self) -> bool {
        // PR-1's executable scalar relation is proven only in the standard
        // component-owned convention: both flows are positive outward, hence
        // their algebraic values sum to zero. Common-frame orientations need
        // an explicit public pullback before they may be interconnected.
        matches!(
            (self, other),
            (Self::OutwardFromOwner, Self::OutwardFromOwner)
        )
    }
}

/// Basis, reference frame, and sign convention for a port value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoordinateBinding {
    basis: StableId,
    frame: StableId,
    orientation: PortOrientation,
}

impl CoordinateBinding {
    /// Bind a port to an explicit basis, frame, and orientation.
    #[must_use]
    pub fn new(basis: StableId, frame: StableId, orientation: PortOrientation) -> Self {
        Self {
            basis,
            frame,
            orientation,
        }
    }

    /// Declared basis identifier.
    #[must_use]
    pub fn basis(&self) -> &StableId {
        &self.basis
    }

    /// Declared frame identifier.
    #[must_use]
    pub fn frame(&self) -> &StableId {
        &self.frame
    }

    /// Declared positive orientation.
    #[must_use]
    pub fn orientation(&self) -> PortOrientation {
        self.orientation
    }
}

/// A deterministic port timestamp in a named logical clock domain.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortTimestamp {
    clock: StableId,
    tick: u64,
}

impl PortTimestamp {
    /// A timestamp. The clock defines the unit and epoch of `tick`.
    #[must_use]
    pub fn new(clock: StableId, tick: u64) -> Self {
        Self { clock, tick }
    }

    /// Clock-domain identifier.
    #[must_use]
    pub fn clock(&self) -> &StableId {
        &self.clock
    }

    /// Logical clock tick.
    #[must_use]
    pub fn tick(&self) -> u64 {
        self.tick
    }
}

/// How effort and flow coordinates are contracted into instantaneous power.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PowerPairing {
    /// Scalar multiplication.
    ScalarProduct,
    /// Euclidean component-wise dot product.
    EuclideanDot,
    /// A declared field duality/integral pairing with its integration-measure
    /// dimensions (for example area for a boundary traction/velocity pair).
    FieldDuality {
        /// Dimensions contributed by the integration measure.
        measure_dimensions: Dims,
    },
}

/// Conserved or audited quantities transported by a port.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ConservationRole {
    /// Energy/power exchange.
    Energy,
    /// Total mass.
    Mass,
    /// Amount of substance or species/element amount.
    Amount,
    /// Linear momentum.
    LinearMomentum,
    /// Angular momentum.
    AngularMomentum,
    /// Entropy.
    Entropy,
    /// Electric charge.
    ElectricCharge,
}

/// Versioned schema for one typed effort/flow port.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortSchema {
    version: u16,
    id: StableId,
    kind: PortKind,
    effort_dimensions: Dims,
    flow_dimensions: Dims,
    shape: PortValueShape,
    coordinates: CoordinateBinding,
    power_pairing: PowerPairing,
    timestamp: PortTimestamp,
    conservation_roles: Vec<ConservationRole>,
}

impl PortSchema {
    /// Construct and structurally admit a v2 port schema.
    ///
    /// This PR-1 admission proves only that the shape/pairing agree and that
    /// effort × flow (plus the declared measure for field duality) has
    /// dimensions of power. Port-kind-specific semantic dimension vocabularies
    /// are the PR-2 gate.
    ///
    /// # Errors
    ///
    /// Returns a structured error for dimension overflow, a non-power
    /// dimension product, an incompatible shape/pairing, or omission of the
    /// mandatory energy conservation role.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: StableId,
        kind: PortKind,
        effort_dimensions: Dims,
        flow_dimensions: Dims,
        shape: PortValueShape,
        coordinates: CoordinateBinding,
        power_pairing: PowerPairing,
        timestamp: PortTimestamp,
        conservation_roles: impl IntoIterator<Item = ConservationRole>,
    ) -> Result<Self, CoupleError> {
        let pairing_matches_shape = matches!(
            (shape, power_pairing),
            (PortValueShape::Scalar, PowerPairing::ScalarProduct)
                | (PortValueShape::Vector(_), PowerPairing::EuclideanDot)
                | (PortValueShape::Tensor { .. }, PowerPairing::EuclideanDot)
                | (
                    PortValueShape::Field { .. },
                    PowerPairing::FieldDuality { .. }
                )
        );
        if !pairing_matches_shape {
            return Err(CoupleError::PairingShapeMismatch {
                shape,
                pairing: power_pairing,
            });
        }

        let pointwise_product = effort_dimensions.checked_plus(flow_dimensions).ok_or(
            CoupleError::PortDimensionOverflow {
                effort: effort_dimensions,
                flow: flow_dimensions,
            },
        )?;
        let product = match power_pairing {
            PowerPairing::FieldDuality { measure_dimensions } => pointwise_product
                .checked_plus(measure_dimensions)
                .ok_or(CoupleError::PortMeasureDimensionOverflow {
                    pointwise_product,
                    measure: measure_dimensions,
                })?,
            PowerPairing::ScalarProduct | PowerPairing::EuclideanDot => pointwise_product,
        };
        if product != Power::DIMS {
            return Err(CoupleError::PortPowerDimensionMismatch {
                effort: effort_dimensions,
                flow: flow_dimensions,
                product,
            });
        }

        let mut conservation_roles: Vec<_> = conservation_roles.into_iter().collect();
        conservation_roles.sort_unstable();
        conservation_roles.dedup();
        if !conservation_roles.contains(&ConservationRole::Energy) {
            return Err(CoupleError::MissingEnergyConservationRole);
        }

        Ok(Self {
            version: PORT_SCHEMA_VERSION,
            id,
            kind,
            effort_dimensions,
            flow_dimensions,
            shape,
            coordinates,
            power_pairing,
            timestamp,
            conservation_roles,
        })
    }

    /// Schema version.
    #[must_use]
    pub fn version(&self) -> u16 {
        self.version
    }

    /// Stable port identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// Physical port vocabulary entry.
    #[must_use]
    pub fn kind(&self) -> PortKind {
        self.kind
    }

    /// Effort dimensions.
    #[must_use]
    pub fn effort_dimensions(&self) -> Dims {
        self.effort_dimensions
    }

    /// Flow dimensions.
    #[must_use]
    pub fn flow_dimensions(&self) -> Dims {
        self.flow_dimensions
    }

    /// Value/field shape.
    #[must_use]
    pub fn shape(&self) -> PortValueShape {
        self.shape
    }

    /// Coordinate binding.
    #[must_use]
    pub fn coordinates(&self) -> &CoordinateBinding {
        &self.coordinates
    }

    /// Power contraction.
    #[must_use]
    pub fn power_pairing(&self) -> PowerPairing {
        self.power_pairing
    }

    /// Clock/timestamp binding.
    #[must_use]
    pub fn timestamp(&self) -> &PortTimestamp {
        &self.timestamp
    }

    /// Canonically sorted, duplicate-free conservation roles.
    #[must_use]
    pub fn conservation_roles(&self) -> &[ConservationRole] {
        &self.conservation_roles
    }

    fn first_conjugacy_mismatch(&self, other: &Self) -> Option<&'static str> {
        if self.id == other.id {
            return Some("stable_id");
        }
        if self.kind != other.kind {
            return Some("kind");
        }
        if self.effort_dimensions != other.effort_dimensions {
            return Some("effort_dimensions");
        }
        if self.flow_dimensions != other.flow_dimensions {
            return Some("flow_dimensions");
        }
        if self.shape != other.shape {
            return Some("shape");
        }
        if self.coordinates.basis != other.coordinates.basis {
            return Some("basis");
        }
        if self.coordinates.frame != other.coordinates.frame {
            return Some("frame");
        }
        if !self
            .coordinates
            .orientation
            .composes_with(other.coordinates.orientation)
        {
            return Some("orientation");
        }
        if self.power_pairing != other.power_pairing {
            return Some("power_pairing");
        }
        if self.timestamp != other.timestamp {
            return Some("clock_timestamp");
        }
        if self.conservation_roles != other.conservation_roles {
            return Some("conservation_roles");
        }
        None
    }
}

/// A power port: an effort/flow pair. `power = effort × flow`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Port {
    /// The effort variable (force / pressure / temperature).
    pub effort: f64,
    /// The flow variable (velocity / flux / entropy flow).
    pub flow: f64,
    /// The physical type.
    pub kind: PortKind,
}

impl Port {
    /// A port.
    #[must_use]
    pub fn new(effort: f64, flow: f64, kind: PortKind) -> Port {
        Port { effort, flow, kind }
    }

    /// The instantaneous power `effort × flow`.
    #[must_use]
    pub fn power(&self) -> f64 {
        self.effort * self.flow
    }

    /// Are two ports POWER-CONJUGATE (composable) — the same physical
    /// effort/flow type? (The composition-time type discipline.)
    #[must_use]
    pub fn conjugate_to(&self, other: &Port) -> bool {
        self.kind == other.kind
    }
}

/// A scalar runtime value bound to an explicit [`PortSchema`].
#[derive(Debug, Clone, PartialEq)]
pub struct SchemaPort {
    schema: PortSchema,
    effort: f64,
    flow: f64,
}

impl SchemaPort {
    fn new(schema: PortSchema, effort: f64, flow: f64) -> Self {
        Self {
            schema,
            effort,
            flow,
        }
    }

    /// Explicit schema carried by this value.
    #[must_use]
    pub fn schema(&self) -> &PortSchema {
        &self.schema
    }

    /// Scalar effort value in coherent SI units.
    #[must_use]
    pub fn effort(&self) -> f64 {
        self.effort
    }

    /// Scalar flow value in coherent SI units.
    #[must_use]
    pub fn flow(&self) -> f64 {
        self.flow
    }

    /// Instantaneous scalar power.
    #[must_use]
    pub fn power(&self) -> f64 {
        self.effort * self.flow
    }
}

/// A schema-bound two-port Dirac interconnection.
#[derive(Debug, Clone, PartialEq)]
pub struct SchemaInterconnection {
    /// Side A, with caller-supplied flow sign.
    pub port_a: SchemaPort,
    /// Side B, with the balancing flow sign.
    pub port_b: SchemaPort,
    /// Net scalar interface power; zero for finite admitted inputs.
    pub interface_power: f64,
}

/// A conservative Dirac/Stokes–Dirac junction descriptor.
///
/// PR-1 implements the exact scalar two-port seed. Multi-port matrix/field
/// relations remain a later operator lane, but the schema admission is already
/// general over scalar, vector, tensor, and field shapes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConservativeJunction {
    id: StableId,
    port_a: PortSchema,
    port_b: PortSchema,
}

impl ConservativeJunction {
    /// Admit a lossless two-port relation between conjugate schemas.
    ///
    /// # Errors
    /// [`CoupleError::IncompatiblePortSchemas`] localizes the first metadata
    /// mismatch, including reused port IDs and any orientation other than the
    /// PR-1 scalar seed's outward-from-owner convention on both sides.
    pub fn new(id: StableId, port_a: PortSchema, port_b: PortSchema) -> Result<Self, CoupleError> {
        if id == port_a.id || id == port_b.id {
            return Err(CoupleError::DuplicateIdentity {
                id: id.as_str().to_string(),
            });
        }
        if let Some(field) = port_a.first_conjugacy_mismatch(&port_b) {
            return Err(CoupleError::IncompatiblePortSchemas {
                a: port_a.id.as_str().to_string(),
                b: port_b.id.as_str().to_string(),
                field,
            });
        }
        Ok(Self { id, port_a, port_b })
    }

    /// Junction identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// The admitted pair, in deterministic A/B order.
    #[must_use]
    pub fn ports(&self) -> (&PortSchema, &PortSchema) {
        (&self.port_a, &self.port_b)
    }

    /// Evaluate the migrated scalar seed: shared effort, opposite flow.
    ///
    /// # Errors
    /// Refuses a non-scalar schema or non-finite runtime input. This keeps the
    /// schema-bound path fail-closed; the legacy raw [`Port`] path remains for
    /// migration comparison only.
    pub fn interconnect_scalar(
        &self,
        effort: f64,
        flow: f64,
    ) -> Result<SchemaInterconnection, CoupleError> {
        if self.port_a.shape != PortValueShape::Scalar {
            return Err(CoupleError::ScalarOperationRequiresScalarPort {
                id: self.port_a.id.as_str().to_string(),
                shape: self.port_a.shape,
            });
        }
        if !effort.is_finite() {
            return Err(CoupleError::NonFinitePortValue { field: "effort" });
        }
        if !flow.is_finite() {
            return Err(CoupleError::NonFinitePortValue { field: "flow" });
        }
        let side_power = effort * flow;
        if !side_power.is_finite() {
            return Err(CoupleError::NonFinitePortValue { field: "power" });
        }
        let port_a = SchemaPort::new(self.port_a.clone(), effort, flow);
        let port_b = SchemaPort::new(self.port_b.clone(), effort, -flow);
        Ok(SchemaInterconnection {
            interface_power: side_power + -side_power,
            port_a,
            port_b,
        })
    }

    fn require_added_mass_fixture_schema(&self) -> Result<(), CoupleError> {
        if self.port_a.kind != PortKind::MechanicalForceVelocity {
            return Err(CoupleError::AddedMassFixtureRequiresMechanicalPort {
                kind: self.port_a.kind,
            });
        }
        if self.port_a.shape != PortValueShape::Scalar {
            return Err(CoupleError::ScalarOperationRequiresScalarPort {
                id: self.port_a.id.as_str().to_string(),
                shape: self.port_a.shape,
            });
        }
        Ok(())
    }

    /// Run the legacy fixed-relaxation added-mass fixture through this
    /// schema-bound mechanical junction.
    ///
    /// This is a migration bridge, not a general FSI operator. It lets retained
    /// goldens prove that PortSchema v2 preserves the original scalar fixture
    /// bit-for-bit while downstream domain adapters migrate.
    ///
    /// # Errors
    /// Refuses a non-mechanical or non-scalar junction and propagates
    /// [`CoupleError::NonFinitePortValue`] if an iteration residual cannot be
    /// represented by the schema-bound scalar exchange.
    pub fn iterate_added_mass_fixed(
        &self,
        mu: f64,
        c: f64,
        x0: f64,
        omega: f64,
        max_steps: usize,
        tol: f64,
    ) -> Result<FsiResult, CoupleError> {
        self.require_added_mass_fixture_schema()?;
        let mut x = x0;
        for step in 1..=max_steps {
            let raw_residual = interface_map(mu, c, x) - x;
            // Unit coherent-SI effort is a migration witness: the numerical
            // residual is round-tripped through the typed Dirac pair without
            // changing its bits, and the balancing side is exercised at every
            // iteration. This does not turn the scalar map into a physical FSI
            // discretization.
            let exchange = self.interconnect_scalar(1.0, raw_residual)?;
            let residual = exchange.port_a.flow();
            x += omega * residual;
            if !x.is_finite() || x.abs() > BLOWUP {
                return Ok(FsiResult {
                    converged: false,
                    steps: step,
                    solution: x,
                    final_residual: f64::INFINITY,
                });
            }
            if residual.abs() < tol {
                return Ok(FsiResult {
                    converged: true,
                    steps: step,
                    solution: x,
                    final_residual: residual.abs(),
                });
            }
        }
        let raw_residual = interface_map(mu, c, x) - x;
        let residual = self
            .interconnect_scalar(1.0, raw_residual)?
            .port_a
            .flow()
            .abs();
        Ok(FsiResult {
            converged: residual < tol,
            steps: max_steps,
            solution: x,
            final_residual: residual,
        })
    }

    /// Run the legacy Aitken added-mass fixture through this schema-bound
    /// mechanical junction.
    ///
    /// # Errors
    /// Refuses a non-mechanical or non-scalar junction and propagates
    /// [`CoupleError::NonFinitePortValue`] if an iteration residual cannot be
    /// represented by the schema-bound scalar exchange.
    #[allow(clippy::too_many_arguments)]
    pub fn iterate_added_mass_aitken(
        &self,
        mu: f64,
        c: f64,
        x0: f64,
        omega_init: f64,
        omega_max: f64,
        max_steps: usize,
        tol: f64,
    ) -> Result<FsiResult, CoupleError> {
        self.require_added_mass_fixture_schema()?;
        let mut x = x0;
        let mut aitken = AitkenRelaxation::new(omega_init, omega_max);
        for step in 1..=max_steps {
            let raw_residual = interface_map(mu, c, x) - x;
            let exchange = self.interconnect_scalar(1.0, raw_residual)?;
            let residual = exchange.port_a.flow();
            if residual.abs() < tol {
                return Ok(FsiResult {
                    converged: true,
                    steps: step,
                    solution: x,
                    final_residual: residual.abs(),
                });
            }
            let omega = aitken.next_omega(residual);
            x += omega * residual;
            if !x.is_finite() || x.abs() > BLOWUP {
                return Ok(FsiResult {
                    converged: false,
                    steps: step,
                    solution: x,
                    final_residual: f64::INFINITY,
                });
            }
        }
        let raw_residual = interface_map(mu, c, x) - x;
        let residual = self
            .interconnect_scalar(1.0, raw_residual)?
            .port_a
            .flow()
            .abs();
        Ok(FsiResult {
            converged: residual < tol,
            steps: max_steps,
            solution: x,
            final_residual: residual,
        })
    }
}

/// Thermodynamic potential represented by a storage element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StoragePotential {
    /// Hamiltonian energy state.
    Hamiltonian,
    /// Free-energy state under an explicit thermodynamic chart.
    FreeEnergy,
}

/// A typed storage primitive.
///
/// The state and constitutive-gradient operator are durable references, not
/// hidden executable closures; domain crates implement the public operator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StorageElement {
    id: StableId,
    port: PortSchema,
    potential: StoragePotential,
    state_schema: StableId,
    state_dimension: NonZeroUsize,
    constitutive_gradient: StableId,
}

impl StorageElement {
    /// Construct a storage descriptor with an explicit state and gradient.
    ///
    /// # Errors
    /// Refuses identity aliasing between the relation and its port.
    pub fn new(
        id: StableId,
        port: PortSchema,
        potential: StoragePotential,
        state_schema: StableId,
        state_dimension: NonZeroUsize,
        constitutive_gradient: StableId,
    ) -> Result<Self, CoupleError> {
        reject_relation_port_alias(&id, &port)?;
        Ok(Self {
            id,
            port,
            potential,
            state_schema,
            state_dimension,
            constitutive_gradient,
        })
    }

    /// Stable relation identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// Exposed power port.
    #[must_use]
    pub fn port(&self) -> &PortSchema {
        &self.port
    }

    /// Hamiltonian or free-energy chart.
    #[must_use]
    pub fn potential(&self) -> StoragePotential {
        self.potential
    }

    /// Durable state-schema identifier.
    #[must_use]
    pub fn state_schema(&self) -> &StableId {
        &self.state_schema
    }

    /// Number of state coordinates consumed by the gradient operator.
    #[must_use]
    pub fn state_dimension(&self) -> NonZeroUsize {
        self.state_dimension
    }

    /// Durable constitutive-gradient operator identifier.
    #[must_use]
    pub fn constitutive_gradient(&self) -> &StableId {
        &self.constitutive_gradient
    }
}

/// Physical family of a dissipative constitutive relation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DissipationLaw {
    /// Electrical or generalized resistance.
    Resistive,
    /// Dry or rate-dependent friction.
    Frictional,
    /// Viscous momentum loss.
    Viscous,
    /// Thermal conduction.
    Conductive,
    /// Plastic flow.
    Plastic,
}

/// Evidence required before a dissipative relation may make a sign claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DissipationEvidence {
    /// A durable monotonicity proof/receipt.
    Monotonicity(StableId),
    /// A durable nonnegative-production proof/receipt.
    NonnegativeProduction(StableId),
}

/// A typed dissipative primitive with mandatory evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DissipativeRelation {
    id: StableId,
    port: PortSchema,
    law: DissipationLaw,
    constitutive_operator: StableId,
    evidence: DissipationEvidence,
}

impl DissipativeRelation {
    /// Construct an evidence-bound dissipative relation.
    ///
    /// # Errors
    /// Refuses identity aliasing between the relation and its port.
    pub fn new(
        id: StableId,
        port: PortSchema,
        law: DissipationLaw,
        constitutive_operator: StableId,
        evidence: DissipationEvidence,
    ) -> Result<Self, CoupleError> {
        reject_relation_port_alias(&id, &port)?;
        Ok(Self {
            id,
            port,
            law,
            constitutive_operator,
            evidence,
        })
    }

    /// Stable relation identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// Exposed power port.
    #[must_use]
    pub fn port(&self) -> &PortSchema {
        &self.port
    }

    /// Constitutive loss family.
    #[must_use]
    pub fn law(&self) -> DissipationLaw {
        self.law
    }

    /// Public constitutive-operator identifier.
    #[must_use]
    pub fn constitutive_operator(&self) -> &StableId {
        &self.constitutive_operator
    }

    /// Proof/receipt that licenses the dissipation sign claim.
    #[must_use]
    pub fn evidence(&self) -> &DissipationEvidence {
        &self.evidence
    }
}

/// What crosses a source/reservoir accounting boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceClass {
    /// Prescribed effort, such as voltage or temperature.
    PrescribedEffort,
    /// Prescribed flow, such as current or heat input.
    PrescribedFlow,
    /// Environment/fuel/body treated as a reservoir exchange.
    Reservoir,
}

/// How a boundary contribution appears in a closed accounting window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryTreatment {
    /// The source term is explicitly included inside the audited model.
    IncludedSourceTerm,
    /// The exchange crosses to an explicitly external reservoir.
    ExternalReservoirExchange,
}

/// Explicit boundary that prevents a source from disappearing from an audit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AccountingBoundary {
    id: StableId,
    coordinates: CoordinateBinding,
    treatment: BoundaryTreatment,
}

impl AccountingBoundary {
    /// Declare a signed accounting boundary.
    #[must_use]
    pub fn new(id: StableId, coordinates: CoordinateBinding, treatment: BoundaryTreatment) -> Self {
        Self {
            id,
            coordinates,
            treatment,
        }
    }

    /// Boundary identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// Basis, frame, and positive contribution convention.
    #[must_use]
    pub fn coordinates(&self) -> &CoordinateBinding {
        &self.coordinates
    }

    /// Whether this is an included source or external exchange.
    #[must_use]
    pub fn treatment(&self) -> BoundaryTreatment {
        self.treatment
    }
}

/// A typed source or reservoir primitive with an explicit audit boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceOrReservoir {
    id: StableId,
    port: PortSchema,
    class: SourceClass,
    boundary: AccountingBoundary,
}

impl SourceOrReservoir {
    /// Construct a source/reservoir descriptor.
    ///
    /// # Errors
    /// Refuses identity aliasing among the relation, port, and boundary, or an
    /// [`CoupleError::AccountingBoundaryMismatch`] between their coordinate
    /// bindings.
    pub fn new(
        id: StableId,
        port: PortSchema,
        class: SourceClass,
        boundary: AccountingBoundary,
    ) -> Result<Self, CoupleError> {
        reject_relation_port_alias(&id, &port)?;
        if id == boundary.id || port.id == boundary.id {
            return Err(CoupleError::DuplicateIdentity {
                id: boundary.id.as_str().to_string(),
            });
        }
        if port.coordinates != boundary.coordinates {
            return Err(CoupleError::AccountingBoundaryMismatch {
                port: port.id.as_str().to_string(),
                boundary: boundary.id.as_str().to_string(),
            });
        }
        Ok(Self {
            id,
            port,
            class,
            boundary,
        })
    }

    /// Stable relation identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// Exposed power port.
    #[must_use]
    pub fn port(&self) -> &PortSchema {
        &self.port
    }

    /// Prescribed-effort, prescribed-flow, or reservoir class.
    #[must_use]
    pub fn class(&self) -> SourceClass {
        self.class
    }

    /// Explicit signed accounting boundary.
    #[must_use]
    pub fn boundary(&self) -> &AccountingBoundary {
        &self.boundary
    }
}

/// Closed vocabulary of the four port-thermodynamic relation primitives.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PortPrimitive {
    /// Lossless topology only.
    ConservativeJunction(ConservativeJunction),
    /// Stored-energy relation.
    StorageElement(StorageElement),
    /// Evidence-bound dissipative relation.
    DissipativeRelation(DissipativeRelation),
    /// Explicit source/reservoir boundary.
    SourceOrReservoir(SourceOrReservoir),
}

fn reject_relation_port_alias(id: &StableId, port: &PortSchema) -> Result<(), CoupleError> {
    if id == &port.id {
        return Err(CoupleError::DuplicateIdentity {
            id: id.as_str().to_string(),
        });
    }
    Ok(())
}

/// A structured coupling failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoupleError {
    /// The ports are not power-conjugate (mismatched physical types).
    IncompatiblePorts {
        /// The first port's kind.
        a: PortKind,
        /// The second port's kind.
        b: PortKind,
    },
    /// A stable identifier was empty or outside the canonical alphabet.
    InvalidStableId {
        /// Rejected value.
        value: String,
    },
    /// A vector, tensor, or field shape had a zero extent.
    EmptyPortShape,
    /// Effort/flow exponent addition overflowed the six-base vector.
    PortDimensionOverflow {
        /// Effort dimensions.
        effort: Dims,
        /// Flow dimensions.
        flow: Dims,
    },
    /// A field pairing's pointwise product plus measure dimensions overflowed.
    PortMeasureDimensionOverflow {
        /// Checked effort × flow dimensions before integration.
        pointwise_product: Dims,
        /// Integration-measure dimensions.
        measure: Dims,
    },
    /// Effort × flow did not have watt dimensions.
    PortPowerDimensionMismatch {
        /// Effort dimensions.
        effort: Dims,
        /// Flow dimensions.
        flow: Dims,
        /// Checked product dimensions.
        product: Dims,
    },
    /// The declared contraction cannot consume the declared shape.
    PairingShapeMismatch {
        /// Value/field shape.
        shape: PortValueShape,
        /// Requested contraction.
        pairing: PowerPairing,
    },
    /// Every power port must declare energy as a conservation role.
    MissingEnergyConservationRole,
    /// Two port schemas disagree on a localized conjugacy field.
    IncompatiblePortSchemas {
        /// First port ID.
        a: String,
        /// Second port ID.
        b: String,
        /// First incompatible schema field.
        field: &'static str,
    },
    /// Stable identities aliased within one relation.
    DuplicateIdentity {
        /// Reused identity.
        id: String,
    },
    /// The scalar seed evaluator was called on a non-scalar schema.
    ScalarOperationRequiresScalarPort {
        /// Port ID.
        id: String,
        /// Actual shape.
        shape: PortValueShape,
    },
    /// A schema-bound runtime value was non-finite.
    NonFinitePortValue {
        /// Offending input field.
        field: &'static str,
    },
    /// The retained scalar added-mass fixture is mechanical-only.
    AddedMassFixtureRequiresMechanicalPort {
        /// Actual port kind.
        kind: PortKind,
    },
    /// A source boundary used a different basis/frame/orientation than its port.
    AccountingBoundaryMismatch {
        /// Source port ID.
        port: String,
        /// Accounting-boundary ID.
        boundary: String,
    },
}

/// A Dirac interconnection of two ports: shared effort, opposite flow — so the
/// interface power `e·f + e·(−f) = 0` EXACTLY (power-conserving by construction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interconnection {
    /// The first (side A) port.
    pub port_a: Port,
    /// The second (side B) port.
    pub port_b: Port,
    /// The net interface power (`0` by construction).
    pub interface_power: f64,
}

/// Interconnect two subsystems at a shared effort and flow through a Dirac
/// structure (effort continuity + flow balance).
///
/// # Errors
/// [`CoupleError::IncompatiblePorts`] if the ports are not power-conjugate.
pub fn interconnect(
    kind_a: PortKind,
    kind_b: PortKind,
    effort: f64,
    flow: f64,
) -> Result<Interconnection, CoupleError> {
    if kind_a != kind_b {
        return Err(CoupleError::IncompatiblePorts {
            a: kind_a,
            b: kind_b,
        });
    }
    let port_a = Port::new(effort, flow, kind_a);
    let port_b = Port::new(effort, -flow, kind_b);
    Ok(Interconnection {
        interface_power: port_a.power() + port_b.power(),
        port_a,
        port_b,
    })
}

/// The net interface power of a set of ports (`Σ effort·flow`) — `0` for a
/// power-conserving interconnection.
#[must_use]
pub fn interface_power(ports: &[Port]) -> f64 {
    ports.iter().map(Port::power).sum()
}

/// A caller-fed interface-balance audit.
///
/// The legacy [`EnergyAudit::is_passive`] name checks only whether every
/// recorded scalar interface imbalance stays within tolerance. It is not a
/// whole-system passivity certificate.
#[derive(Debug, Clone, Default)]
pub struct EnergyAudit {
    balances: Vec<f64>,
}

impl EnergyAudit {
    /// A fresh audit.
    #[must_use]
    pub fn new() -> EnergyAudit {
        EnergyAudit {
            balances: Vec::new(),
        }
    }

    /// Record one exchange's net interface power.
    pub fn record(&mut self, net_interface_power: f64) {
        self.balances.push(net_interface_power);
    }

    /// The worst interface power generation seen (the bug-alarm metric).
    ///
    /// A recorded NaN interface power means the coupling numerically broke
    /// down — the single worst thing this audit exists to catch. `f64::max`
    /// SILENTLY DROPS NaN (`f64::max(0.0, NaN) == 0.0`), so a plain fold would
    /// report zero imbalance and let the legacy `is_passive` predicate return
    /// true for a blown-up coupling. Poison instead: any NaN balance makes the
    /// metric NaN, and `NaN <= tol` is false, so the audit fails closed.
    /// (`±∞` already survives `f64::max` and alarms correctly.)
    #[must_use]
    pub fn max_generation(&self) -> f64 {
        if self.balances.iter().any(|b| b.is_nan()) {
            return f64::NAN;
        }
        self.balances.iter().map(|b| b.abs()).fold(0.0, f64::max)
    }

    /// Is every recorded interface-power imbalance within `tol`?
    ///
    /// This legacy name does not establish component or closed-window
    /// passivity; callers must audit those obligations separately.
    #[must_use]
    pub fn is_passive(&self, tol: f64) -> bool {
        self.max_generation() <= tol
    }
}

/// Scalar Aitken (Δ²) dynamic relaxation for the strongly-coupled interface
/// fixed point.
#[derive(Debug, Clone)]
pub struct AitkenRelaxation {
    omega: f64,
    omega_max: f64,
    prev_residual: Option<f64>,
}

impl AitkenRelaxation {
    /// A relaxer with an initial ω and a magnitude cap.
    #[must_use]
    pub fn new(omega_init: f64, omega_max: f64) -> AitkenRelaxation {
        AitkenRelaxation {
            omega: omega_init,
            omega_max,
            prev_residual: None,
        }
    }

    /// The Aitken relaxation factor for the current residual:
    /// `ωₖ = −ωₖ₋₁ · rₖ₋₁ / (rₖ − rₖ₋₁)` (scalar), magnitude-capped.
    pub fn next_omega(&mut self, residual: f64) -> f64 {
        if let Some(prev) = self.prev_residual {
            let dr = residual - prev;
            if dr.abs() > 1e-14 {
                let w = -self.omega * prev / dr;
                self.omega = w.clamp(-self.omega_max, self.omega_max);
            }
        }
        self.prev_residual = Some(residual);
        self.omega
    }
}

/// The result of an FSI interface fixed-point iteration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FsiResult {
    /// Did it converge (residual below `tol`) without blowing up?
    pub converged: bool,
    /// Iterations taken (or `max_steps` if it did not converge).
    pub steps: usize,
    /// The final interface value.
    pub solution: f64,
    /// The final residual magnitude.
    pub final_residual: f64,
}

// The classic linearized added-mass interface map: H(x) = −μ·x + c, where μ is
// the added-mass ratio (fluid added mass / structural mass). Naive staggering
// (ω = 1) converges only for μ < 1; a dense fluid on a light structure (μ ≥ 1)
// diverges.
fn interface_map(mu: f64, c: f64, x: f64) -> f64 {
    -mu * x + c
}

const BLOWUP: f64 = 1e12;

/// Iterate the added-mass interface fixed point with FIXED under-relaxation
/// `omega`. Diverges (naive staggering, `omega = 1`) when the added-mass ratio
/// `mu >= 1`.
#[must_use]
pub fn iterate_fixed_relaxation(
    mu: f64,
    c: f64,
    x0: f64,
    omega: f64,
    max_steps: usize,
    tol: f64,
) -> FsiResult {
    let mut x = x0;
    for step in 1..=max_steps {
        let r = interface_map(mu, c, x) - x;
        x += omega * r;
        if !x.is_finite() || x.abs() > BLOWUP {
            return FsiResult {
                converged: false,
                steps: step,
                solution: x,
                final_residual: f64::INFINITY,
            };
        }
        if r.abs() < tol {
            return FsiResult {
                converged: true,
                steps: step,
                solution: x,
                final_residual: r.abs(),
            };
        }
    }
    let r = (interface_map(mu, c, x) - x).abs();
    FsiResult {
        converged: r < tol,
        steps: max_steps,
        solution: x,
        final_residual: r,
    }
}

/// Iterate the same interface fixed point with AITKEN dynamic relaxation, which
/// stabilizes and accelerates it even for `mu >= 1` (the added-mass fix).
#[must_use]
pub fn iterate_aitken(
    mu: f64,
    c: f64,
    x0: f64,
    omega_init: f64,
    omega_max: f64,
    max_steps: usize,
    tol: f64,
) -> FsiResult {
    let mut x = x0;
    let mut aitken = AitkenRelaxation::new(omega_init, omega_max);
    for step in 1..=max_steps {
        let r = interface_map(mu, c, x) - x;
        if r.abs() < tol {
            return FsiResult {
                converged: true,
                steps: step,
                solution: x,
                final_residual: r.abs(),
            };
        }
        let omega = aitken.next_omega(r);
        x += omega * r;
        if !x.is_finite() || x.abs() > BLOWUP {
            return FsiResult {
                converged: false,
                steps: step,
                solution: x,
                final_residual: f64::INFINITY,
            };
        }
    }
    let r = (interface_map(mu, c, x) - x).abs();
    FsiResult {
        converged: r < tol,
        steps: max_steps,
        solution: x,
        final_residual: r,
    }
}
