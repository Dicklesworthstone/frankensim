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
//! [`StreamPort`] is the distinct multi-conserved-flux bundle for transported
//! mass, constituent amount, momentum, energy, and entropy. Its single
//! [`StreamEnergyChart`] is structural, and alternate thermodynamic charts
//! require context-bound exact pressure-work or Euler/Legendre crosswalks.
//! [`AccountingWindowSpec`] and [`WindowAuditReport`] separately audit
//! pre-integrated, window-bound storage and port contributions without
//! pretending that a logical clock tick is a physical duration.
//!
//! For the hard, strongly-coupled cases, [`AitkenRelaxation`] gives dynamic
//! interface relaxation: on the classic ADDED-MASS-INSTABILITY fixture (a light
//! structure in a dense fluid) naive staggering diverges, while Aitken-relaxed
//! coupling converges — demonstrated by [`iterate_fixed_relaxation`] vs
//! [`iterate_aitken`]. Deterministic; depends only on the neutral `fs-iface`
//! vocabulary and `fs-qty`'s six-base dimension vector.

use core::num::NonZeroUsize;

use fs_iface::SpaceType;
use fs_qty::chemistry::{ElementId, SpeciesId};
use fs_qty::{
    Amount, Density, Dims, ElectricCharge, Energy, Force, Mass, MassFlowRate, Power, Pressure, Qty,
    Temperature, Velocity, VolumetricFlowRate,
};

/// Current public port-schema version.
pub const PORT_SCHEMA_VERSION: u16 = 2;

/// Current public stream-port version.
pub const STREAM_PORT_VERSION: u16 = 1;

/// Current closed accounting-window schema version.
pub const ACCOUNTING_WINDOW_VERSION: u16 = 1;

/// Moles per second.
pub type AmountFlowRate = Qty<0, 0, -1, 0, 0, 1>;

/// Watts per kelvin.
pub type EntropyFlowRate = Qty<2, 1, -3, -1, 0, 0>;

/// Joules per kelvin.
pub type Entropy = Qty<2, 1, -2, -1, 0, 0>;

/// Joules per kilogram.
pub type SpecificEnergy = Qty<2, 0, -2, 0, 0, 0>;

/// The physical type of a power-conjugate port (its effort/flow pair).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortKind {
    /// Mechanical: force (effort) × velocity (flow).
    MechanicalForceVelocity,
    /// Fluid: pressure (effort) × volumetric flux (flow).
    FluidPressureFlux,
    /// Thermal: temperature (effort) × entropy flow.
    ThermalTemperatureEntropy,
    /// Rotational: torque × angular velocity.
    RotationalTorqueAngularVelocity,
    /// Electrical: voltage × current.
    ElectricalVoltageCurrent,
    /// Magnetic: magnetomotive force × magnetic-flux rate.
    MagneticMmfFluxRate,
    /// Chemical: electrochemical potential × amount flow.
    ChemicalPotentialAmountFlow,
}

impl PortKind {
    /// Whether this kind belongs to the retained raw scalar migration oracle.
    const fn is_legacy_scalar_seed(self) -> bool {
        matches!(
            self,
            Self::MechanicalForceVelocity
                | Self::FluidPressureFlux
                | Self::ThermalTemperatureEntropy
        )
    }

    /// Canonical generalized effort dimensions for this physical kind.
    #[must_use]
    pub const fn canonical_effort_dimensions(self) -> Dims {
        match self {
            Self::MechanicalForceVelocity => Force::DIMS,
            Self::FluidPressureFlux => Pressure::DIMS,
            Self::ThermalTemperatureEntropy => Temperature::DIMS,
            // Torque is N·m. Its dimensions match energy, while its semantic
            // kind remains distinct.
            Self::RotationalTorqueAngularVelocity => Dims([2, 1, -2, 0, 0, 0]),
            // Voltage is W/A.
            Self::ElectricalVoltageCurrent => Dims([2, 1, -3, 0, -1, 0]),
            // Magnetomotive force is ampere-turn; turn is dimensionless.
            Self::MagneticMmfFluxRate => Dims([0, 0, 0, 0, 1, 0]),
            // Electrochemical potential is J/mol.
            Self::ChemicalPotentialAmountFlow => Dims([2, 1, -2, 0, 0, -1]),
        }
    }

    /// Canonical generalized flow dimensions for this physical kind.
    #[must_use]
    pub const fn canonical_flow_dimensions(self) -> Dims {
        match self {
            Self::MechanicalForceVelocity => Velocity::DIMS,
            Self::FluidPressureFlux => VolumetricFlowRate::DIMS,
            // Entropy flow is W/K.
            Self::ThermalTemperatureEntropy => Dims([2, 1, -3, -1, 0, 0]),
            // Angular velocity is rad/s; radian is dimensionless.
            Self::RotationalTorqueAngularVelocity => Dims([0, 0, -1, 0, 0, 0]),
            Self::ElectricalVoltageCurrent => Dims([0, 0, 0, 0, 1, 0]),
            // Magnetic-flux rate is Wb/s.
            Self::MagneticMmfFluxRate => Dims([2, 1, -3, 0, -1, 0]),
            // Amount flow is mol/s.
            Self::ChemicalPotentialAmountFlow => Dims([0, 0, -1, 0, 0, 1]),
        }
    }

    /// Minimum conservation roles admitted at this staged schema version.
    ///
    /// Energy is universal. PR-2 schema-only kinds whose flow is itself the
    /// rate of a named conserved quantity must also declare that quantity;
    /// dimensional compatibility alone is not enough. The original three
    /// scalar seeds retain their PR-1 Energy-only role vectors as migration
    /// goldens; PR-4 owns stronger closed-window transport audits.
    #[must_use]
    pub fn required_conservation_roles(self) -> &'static [ConservationRole] {
        match self {
            Self::MechanicalForceVelocity
            | Self::FluidPressureFlux
            | Self::ThermalTemperatureEntropy
            | Self::MagneticMmfFluxRate => &[ConservationRole::Energy],
            Self::RotationalTorqueAngularVelocity => {
                &[ConservationRole::Energy, ConservationRole::AngularMomentum]
            }
            Self::ElectricalVoltageCurrent => {
                &[ConservationRole::Energy, ConservationRole::ElectricCharge]
            }
            Self::ChemicalPotentialAmountFlow => {
                &[ConservationRole::Energy, ConservationRole::Amount]
            }
        }
    }

    /// Construct one canonical scalar member of this port kind.
    ///
    /// Identity, coordinates, and time remain caller supplied: schema
    /// construction must not invent any of the Five Explicits. Chemical ports
    /// use this for a single species coordinate; multi-species ports use
    /// [`PortValueShape::Vector`] with the same per-component dimensions.
    /// Chemical affinity/reaction-extent pairs are a distinct vocabulary and
    /// are not represented by this kind.
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
            self.canonical_effort_dimensions(),
            self.canonical_flow_dimensions(),
            PortValueShape::Scalar,
            coordinates,
            PowerPairing::ScalarProduct,
            timestamp,
            self.required_conservation_roles().iter().copied(),
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
    /// A field pairing with separate neutral function-space roles for effort
    /// and flow.
    Field {
        /// Component count at each field point.
        components: NonZeroUsize,
        /// FEEC/interface function-space role of the effort coordinate.
        effort_space: SpaceType,
        /// FEEC/interface function-space role of the flow coordinate.
        flow_space: SpaceType,
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
    pub fn field(
        components: usize,
        effort_space: SpaceType,
        flow_space: SpaceType,
    ) -> Result<Self, CoupleError> {
        let components = NonZeroUsize::new(components).ok_or(CoupleError::EmptyPortShape)?;
        Ok(Self::Field {
            components,
            effort_space,
            flow_space,
        })
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

/// Which pointwise field variable receives the integration measure when it is
/// compared with the generalized dimensions of its [`PortKind`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldMeasureSide {
    /// The effort is a density, such as traction; effort × measure is force.
    Effort,
    /// The flow is a density, such as entropy flux; flow × measure is entropy
    /// flow.
    Flow,
}

/// Effort or flow, used by localized dimension diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortVariable {
    /// Effort coordinate.
    Effort,
    /// Flow coordinate.
    Flow,
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
        /// Which pointwise variable is promoted to generalized dimensions by
        /// that measure.
        measure_side: FieldMeasureSide,
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
    /// Admission proves that the shape/pairing agree, effort × flow (plus the
    /// declared measure for field duality) has dimensions of power, and the
    /// measure-adjusted generalized coordinates match the declared port kind.
    /// It does not prove the downstream constitutive law or adapter.
    ///
    /// # Errors
    ///
    /// Returns a structured error for dimension overflow, a non-power
    /// dimension product, kind-specific dimension mismatch, an incompatible
    /// shape/pairing, or omission of a kind-required conservation role.
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

        admit_port_dimensions(kind, effort_dimensions, flow_dimensions, power_pairing)?;

        let mut conservation_roles: Vec<_> = conservation_roles.into_iter().collect();
        conservation_roles.sort_unstable();
        conservation_roles.dedup();
        if !conservation_roles.contains(&ConservationRole::Energy) {
            return Err(CoupleError::MissingEnergyConservationRole);
        }
        for &role in kind.required_conservation_roles() {
            if !conservation_roles.contains(&role) {
                return Err(CoupleError::MissingPortKindConservationRole { kind, role });
            }
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

fn admit_port_dimensions(
    kind: PortKind,
    effort_dimensions: Dims,
    flow_dimensions: Dims,
    power_pairing: PowerPairing,
) -> Result<(), CoupleError> {
    admit_power_dimensions(effort_dimensions, flow_dimensions, power_pairing)?;
    let (generalized_effort, generalized_flow) =
        generalized_port_dimensions(effort_dimensions, flow_dimensions, power_pairing)?;
    let expected_effort = kind.canonical_effort_dimensions();
    if generalized_effort != expected_effort {
        return Err(CoupleError::PortKindDimensionMismatch {
            kind,
            side: PortVariable::Effort,
            expected: expected_effort,
            actual: generalized_effort,
        });
    }
    let expected_flow = kind.canonical_flow_dimensions();
    if generalized_flow != expected_flow {
        return Err(CoupleError::PortKindDimensionMismatch {
            kind,
            side: PortVariable::Flow,
            expected: expected_flow,
            actual: generalized_flow,
        });
    }
    Ok(())
}

fn admit_power_dimensions(
    effort_dimensions: Dims,
    flow_dimensions: Dims,
    power_pairing: PowerPairing,
) -> Result<(), CoupleError> {
    let pointwise_product = effort_dimensions.checked_plus(flow_dimensions).ok_or(
        CoupleError::PortDimensionOverflow {
            effort: effort_dimensions,
            flow: flow_dimensions,
        },
    )?;
    let product = match power_pairing {
        PowerPairing::FieldDuality {
            measure_dimensions, ..
        } => pointwise_product.checked_plus(measure_dimensions).ok_or(
            CoupleError::PortMeasureDimensionOverflow {
                pointwise_product,
                measure: measure_dimensions,
            },
        )?,
        PowerPairing::ScalarProduct | PowerPairing::EuclideanDot => pointwise_product,
    };
    if product != Power::DIMS {
        return Err(CoupleError::PortPowerDimensionMismatch {
            effort: effort_dimensions,
            flow: flow_dimensions,
            product,
        });
    }
    Ok(())
}

fn generalized_port_dimensions(
    effort_dimensions: Dims,
    flow_dimensions: Dims,
    power_pairing: PowerPairing,
) -> Result<(Dims, Dims), CoupleError> {
    match power_pairing {
        PowerPairing::FieldDuality {
            measure_dimensions,
            measure_side: FieldMeasureSide::Effort,
        } => Ok((
            effort_dimensions.checked_plus(measure_dimensions).ok_or(
                CoupleError::PortMeasureApplicationOverflow {
                    side: PortVariable::Effort,
                    dimensions: effort_dimensions,
                    measure: measure_dimensions,
                },
            )?,
            flow_dimensions,
        )),
        PowerPairing::FieldDuality {
            measure_dimensions,
            measure_side: FieldMeasureSide::Flow,
        } => Ok((
            effort_dimensions,
            flow_dimensions.checked_plus(measure_dimensions).ok_or(
                CoupleError::PortMeasureApplicationOverflow {
                    side: PortVariable::Flow,
                    dimensions: flow_dimensions,
                    measure: measure_dimensions,
                },
            )?,
        )),
        PowerPairing::ScalarProduct | PowerPairing::EuclideanDot => {
            Ok((effort_dimensions, flow_dimensions))
        }
    }
}

/// Canonical constituent axis entry carried by a [`StreamPort`].
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StreamConstituentId {
    /// A chemical species amount rate.
    Species(SpeciesId),
    /// An elemental amount rate.
    Element(ElementId),
}

impl StreamConstituentId {
    fn diagnostic_label(&self) -> String {
        match self {
            Self::Species(id) => format!("species:{id}"),
            Self::Element(id) => format!("element:{id}"),
        }
    }
}

/// One signed species or element amount rate in a stream bundle.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamConstituentFlow {
    id: StreamConstituentId,
    amount_flow: AmountFlowRate,
}

impl StreamConstituentFlow {
    /// Construct one finite constituent amount rate.
    ///
    /// # Errors
    ///
    /// [`CoupleError::NonFiniteStreamValue`] when `amount_flow` is not finite.
    pub fn try_new(
        id: StreamConstituentId,
        amount_flow: AmountFlowRate,
    ) -> Result<Self, CoupleError> {
        ensure_finite_stream_value("constituent_amount_flow", None, amount_flow.value())?;
        Ok(Self { id, amount_flow })
    }

    /// Species or element identifier.
    #[must_use]
    pub fn id(&self) -> &StreamConstituentId {
        &self.id
    }

    /// Signed amount rate in moles per second.
    #[must_use]
    pub fn amount_flow(&self) -> AmountFlowRate {
        self.amount_flow
    }
}

/// Normative pressure/deviatoric split used by stream-energy charts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamStressWorkConvention {
    /// Cauchy tension is positive with `sigma = -p I + tau`. Integrated
    /// outward boundary power is `p Q - integral_A((tau n) dot u) dA`.
    /// Accordingly, [`DeviatoricStressWork`] stores the already integrated,
    /// signed `-integral_A((tau n) dot u) dA` contribution.
    CauchyTensionPositiveOutwardPower,
}

/// Context to which an energy chart and every exact crosswalk are bound.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StreamChartBinding {
    port_id: StableId,
    state_schema: StableId,
    constituent_basis: StableId,
    constituent_axis: Vec<StreamConstituentId>,
    chemical_reference_state: StableId,
    coordinates: CoordinateBinding,
    timestamp: PortTimestamp,
    gravity_datum: StableId,
    stress_convention: StreamStressWorkConvention,
}

impl StreamChartBinding {
    /// Bind a chart to one boundary, state, explicit constituent axis,
    /// coordinate convention, clock, gravity datum, and stress convention.
    ///
    /// # Errors
    ///
    /// Refuses an empty axis or a duplicate species/element coordinate.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        port_id: StableId,
        state_schema: StableId,
        constituent_basis: StableId,
        constituent_axis: impl IntoIterator<Item = StreamConstituentId>,
        chemical_reference_state: StableId,
        coordinates: CoordinateBinding,
        timestamp: PortTimestamp,
        gravity_datum: StableId,
        stress_convention: StreamStressWorkConvention,
    ) -> Result<Self, CoupleError> {
        let mut constituent_axis: Vec<_> = constituent_axis.into_iter().collect();
        canonicalize_constituent_axis(&mut constituent_axis)?;
        Ok(Self {
            port_id,
            state_schema,
            constituent_basis,
            constituent_axis,
            chemical_reference_state,
            coordinates,
            timestamp,
            gravity_datum,
            stress_convention,
        })
    }

    /// Stable stream-port identifier.
    #[must_use]
    pub fn port_id(&self) -> &StableId {
        &self.port_id
    }

    /// State schema used by the chart.
    #[must_use]
    pub fn state_schema(&self) -> &StableId {
        &self.state_schema
    }

    /// Ordered species/element basis artifact.
    #[must_use]
    pub fn constituent_basis(&self) -> &StableId {
        &self.constituent_basis
    }

    /// Canonically ordered species/element coordinates covered by the basis.
    #[must_use]
    pub fn constituent_axis(&self) -> &[StreamConstituentId] {
        &self.constituent_axis
    }

    /// Chemical reference state shared by every energy term.
    #[must_use]
    pub fn chemical_reference_state(&self) -> &StableId {
        &self.chemical_reference_state
    }

    /// Coordinate, frame, and orientation binding.
    #[must_use]
    pub fn coordinates(&self) -> &CoordinateBinding {
        &self.coordinates
    }

    /// Logical clock and tick.
    #[must_use]
    pub fn timestamp(&self) -> &PortTimestamp {
        &self.timestamp
    }

    /// Gravity datum used by the `g z` term.
    #[must_use]
    pub fn gravity_datum(&self) -> &StableId {
        &self.gravity_datum
    }

    /// Cauchy-stress sign/splitting convention.
    #[must_use]
    pub fn stress_convention(&self) -> StreamStressWorkConvention {
        self.stress_convention
    }

    fn first_mismatch(&self, other: &Self) -> Option<&'static str> {
        if self.port_id != other.port_id {
            return Some("port_id");
        }
        if self.state_schema != other.state_schema {
            return Some("state_schema");
        }
        if self.constituent_basis != other.constituent_basis {
            return Some("constituent_basis");
        }
        if self.constituent_axis != other.constituent_axis {
            return Some("constituent_axis");
        }
        if self.chemical_reference_state != other.chemical_reference_state {
            return Some("chemical_reference_state");
        }
        if self.coordinates != other.coordinates {
            return Some("coordinates");
        }
        if self.timestamp != other.timestamp {
            return Some("clock_timestamp");
        }
        if self.gravity_datum != other.gravity_datum {
            return Some("gravity_datum");
        }
        if self.stress_convention != other.stress_convention {
            return Some("stress_convention");
        }
        None
    }
}

/// Exact thermodynamic identity named by a retained proof reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamIdentity {
    /// `h = e + p/rho` and its pressure-work rate form.
    PressureWork,
    /// Euler extensivity plus the enthalpy/Helmholtz/Gibbs transforms.
    EulerLegendre,
    /// Partition proving that a state potential excludes an explicit chemical
    /// power contribution.
    ChemicalEnergyPartition,
}

/// Durable reference to externally retained exact-identity evidence.
///
/// This crate checks identity kind, complete context binding, and the numeric
/// equality where the chart exposes all terms. It does not execute the named
/// verifier or validate the referenced artifact contents.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactIdentityProofRef {
    identity: StreamIdentity,
    receipt_id: StableId,
    verifier_id: StableId,
    statement_digest: StableId,
    binding: Box<StreamChartBinding>,
}

impl ExactIdentityProofRef {
    /// Name one retained exact-identity proof artifact.
    #[must_use]
    pub fn new(
        identity: StreamIdentity,
        receipt_id: StableId,
        verifier_id: StableId,
        statement_digest: StableId,
        binding: StreamChartBinding,
    ) -> Self {
        Self {
            identity,
            receipt_id,
            verifier_id,
            statement_digest,
            binding: Box::new(binding),
        }
    }

    /// Identity established by the referenced artifact.
    #[must_use]
    pub fn identity(&self) -> StreamIdentity {
        self.identity
    }

    /// Context bound by the proof statement.
    #[must_use]
    pub fn binding(&self) -> &StreamChartBinding {
        &self.binding
    }

    /// Retained proof receipt identifier.
    #[must_use]
    pub fn receipt_id(&self) -> &StableId {
        &self.receipt_id
    }

    /// Verifier implementation that produced the receipt.
    #[must_use]
    pub fn verifier_id(&self) -> &StableId {
        &self.verifier_id
    }

    /// Digest of the exact identity statement.
    #[must_use]
    pub fn statement_digest(&self) -> &StableId {
        &self.statement_digest
    }
}

/// Unadmitted chemical-energy ownership declaration.
#[derive(Debug, Clone, PartialEq)]
pub enum ChemicalEnergyInput {
    /// The selected state potential already contains chemical energy.
    IncludedInStatePotential {
        /// Chemical reference state used by that potential.
        reference_state: StableId,
    },
    /// The state potential excludes chemistry and a proved species-potential
    /// contribution is added exactly once.
    ExplicitSpeciesPotentials {
        /// Chemical reference state used by the explicit contribution.
        reference_state: StableId,
        /// Signed `sum(mu_i * n_dot_i)` contribution in watts.
        power_rate: Power,
        /// Proof that the selected state potential excludes this contribution.
        partition_proof: ExactIdentityProofRef,
    },
    /// Invalid request retained so admission can issue a structured
    /// double-counting refusal.
    IncludedAndExplicitSpeciesPotentials {
        /// Chemical reference state shared by the two conflicting terms.
        reference_state: StableId,
        /// Duplicate explicit chemical contribution.
        power_rate: Power,
        /// Claimed partition evidence, which cannot legalize double ownership.
        partition_proof: ExactIdentityProofRef,
    },
}

/// Exclusive admitted ownership of chemical energy in one chart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChemicalEnergyMode {
    /// Chemical energy is embedded in the state potential.
    IncludedInStatePotential,
    /// Chemical energy is a separately proved species-potential contribution.
    ExplicitSpeciesPotentials,
}

/// Admitted, exactly-once chemical-energy accounting.
#[derive(Debug, Clone, PartialEq)]
pub struct ChemicalEnergyAccounting {
    binding: StreamChartBinding,
    mode: ChemicalEnergyMode,
    explicit_power_rate: Power,
    partition_proof: Option<ExactIdentityProofRef>,
}

impl ChemicalEnergyAccounting {
    /// Admit one exclusive chemical-energy ownership declaration.
    ///
    /// # Errors
    ///
    /// Refuses reference/proof binding mismatches, non-finite explicit power,
    /// or a declaration that embeds and separately adds chemical energy.
    pub fn try_new(
        binding: &StreamChartBinding,
        input: ChemicalEnergyInput,
    ) -> Result<Self, CoupleError> {
        match input {
            ChemicalEnergyInput::IncludedInStatePotential { reference_state } => {
                require_stable_id_match(
                    "chemical_reference_state",
                    &binding.chemical_reference_state,
                    &reference_state,
                )?;
                Ok(Self {
                    binding: binding.clone(),
                    mode: ChemicalEnergyMode::IncludedInStatePotential,
                    explicit_power_rate: Power::new(0.0),
                    partition_proof: None,
                })
            }
            ChemicalEnergyInput::ExplicitSpeciesPotentials {
                reference_state,
                power_rate,
                partition_proof,
            } => {
                require_species_constituent_axis(binding)?;
                require_stable_id_match(
                    "chemical_reference_state",
                    &binding.chemical_reference_state,
                    &reference_state,
                )?;
                require_identity_proof(
                    binding,
                    &partition_proof,
                    StreamIdentity::ChemicalEnergyPartition,
                )?;
                ensure_finite_stream_value("explicit_chemical_power", None, power_rate.value())?;
                Ok(Self {
                    binding: binding.clone(),
                    mode: ChemicalEnergyMode::ExplicitSpeciesPotentials,
                    explicit_power_rate: power_rate,
                    partition_proof: Some(partition_proof),
                })
            }
            ChemicalEnergyInput::IncludedAndExplicitSpeciesPotentials { .. } => {
                Err(CoupleError::DoubleCountedChemicalEnergy)
            }
        }
    }

    /// Admitted ownership mode.
    #[must_use]
    pub fn mode(&self) -> ChemicalEnergyMode {
        self.mode
    }

    /// Complete chart context bound by this ownership decision.
    #[must_use]
    pub fn binding(&self) -> &StreamChartBinding {
        &self.binding
    }

    /// Separate chemical power contribution, or zero when embedded.
    #[must_use]
    pub fn explicit_power_rate(&self) -> Power {
        self.explicit_power_rate
    }

    /// Exact partition proof when chemical power is explicit.
    #[must_use]
    pub fn partition_proof(&self) -> Option<&ExactIdentityProofRef> {
        self.partition_proof.as_ref()
    }
}

/// Velocity and gravitational terms shared by stream energy charts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamKinematics {
    velocity: [Velocity; 3],
    specific_gravitational_energy: SpecificEnergy,
}

impl StreamKinematics {
    /// Admit finite velocity components and a finite `g z` term.
    ///
    /// # Errors
    ///
    /// [`CoupleError::NonFiniteStreamValue`] for any non-finite input.
    pub fn try_new(
        velocity: [Velocity; 3],
        specific_gravitational_energy: SpecificEnergy,
    ) -> Result<Self, CoupleError> {
        for (index, component) in velocity.into_iter().enumerate() {
            ensure_finite_stream_value("velocity", Some(index), component.value())?;
        }
        ensure_finite_stream_value(
            "specific_gravitational_energy",
            None,
            specific_gravitational_energy.value(),
        )?;
        let kinematics = Self {
            velocity,
            specific_gravitational_energy,
        };
        ensure_finite_stream_value(
            "specific_transport_energy",
            None,
            kinematics.specific_transport_energy(),
        )?;
        Ok(kinematics)
    }

    fn specific_transport_energy(self) -> f64 {
        let speed_squared = self
            .velocity
            .into_iter()
            .map(|component| component.value() * component.value())
            .sum::<f64>();
        0.5 * speed_squared + self.specific_gravitational_energy.value()
    }
}

/// Signed deviatoric-stress work with explicit operator/evidence ownership.
#[derive(Debug, Clone, PartialEq)]
pub struct DeviatoricStressWork {
    binding: StreamChartBinding,
    power_rate: Power,
    operator_id: StableId,
    evidence_id: StableId,
}

impl DeviatoricStressWork {
    /// Declare one finite deviatoric-stress work rate bound to a stream chart.
    ///
    /// # Errors
    ///
    /// [`CoupleError::NonFiniteStreamValue`] for a non-finite rate.
    pub fn try_new(
        binding: &StreamChartBinding,
        power_rate: Power,
        operator_id: StableId,
        evidence_id: StableId,
    ) -> Result<Self, CoupleError> {
        ensure_finite_stream_value("deviatoric_stress_work", None, power_rate.value())?;
        Ok(Self {
            binding: binding.clone(),
            power_rate,
            operator_id,
            evidence_id,
        })
    }

    /// Complete boundary/frame/stress context covered by the evidence.
    #[must_use]
    pub fn binding(&self) -> &StreamChartBinding {
        &self.binding
    }

    /// Signed work rate in watts.
    #[must_use]
    pub fn power_rate(&self) -> Power {
        self.power_rate
    }

    /// Operator that produced the signed stress-work contribution.
    #[must_use]
    pub fn operator_id(&self) -> &StableId {
        &self.operator_id
    }

    /// Retained evidence for that operator contribution.
    #[must_use]
    pub fn evidence_id(&self) -> &StableId {
        &self.evidence_id
    }
}

/// Canonical moving-stream enthalpy chart.
///
/// Its signed boundary contribution is
/// `mass_flow * (h + |velocity|^2 / 2 + g z) + deviatoric_work`, plus an
/// explicit chemical term only when a partition proof says `h` excludes it.
#[derive(Debug, Clone, PartialEq)]
pub struct MovingStreamEnthalpyChart {
    binding: StreamChartBinding,
    specific_enthalpy: SpecificEnergy,
    kinematics: StreamKinematics,
    deviatoric_work: DeviatoricStressWork,
    chemical_energy: ChemicalEnergyAccounting,
}

impl MovingStreamEnthalpyChart {
    /// Admit a canonical enthalpy chart.
    ///
    /// # Errors
    ///
    /// Refuses non-finite enthalpy or chemical ownership bound to a different
    /// stream context.
    pub fn try_new(
        binding: StreamChartBinding,
        specific_enthalpy: SpecificEnergy,
        kinematics: StreamKinematics,
        deviatoric_work: DeviatoricStressWork,
        chemical_energy: ChemicalEnergyAccounting,
    ) -> Result<Self, CoupleError> {
        ensure_finite_stream_value("specific_enthalpy", None, specific_enthalpy.value())?;
        require_stream_binding(&binding, deviatoric_work.binding())?;
        require_stream_binding(&binding, chemical_energy.binding())?;
        Ok(Self {
            binding,
            specific_enthalpy,
            kinematics,
            deviatoric_work,
            chemical_energy,
        })
    }

    fn energy_rate(&self, mass_flow: MassFlowRate) -> Power {
        Power::new(
            mass_flow.value()
                * (self.specific_enthalpy.value() + self.kinematics.specific_transport_energy())
                + self.deviatoric_work.power_rate.value()
                + self.chemical_energy.explicit_power_rate.value(),
        )
    }
}

/// Exact pressure-work crosswalk between enthalpy and internal-energy charts.
#[derive(Debug, Clone, PartialEq)]
pub struct PressureWorkCrosswalk {
    proof: ExactIdentityProofRef,
    mass_flow: MassFlowRate,
    density: Density,
    specific_enthalpy: SpecificEnergy,
    specific_internal_energy: SpecificEnergy,
    pressure: Pressure,
    volume_flow: VolumetricFlowRate,
    pressure_work_rate: Power,
}

impl PressureWorkCrosswalk {
    /// Recompute and admit the exact identities `h = e + p/rho`,
    /// `volume_flow = mass_flow/rho`, and
    /// `mass_flow * (h - e) = p * volume_flow`.
    ///
    /// # Errors
    ///
    /// Refuses a wrong proof kind, non-finite/non-positive state, or any
    /// identity that is not bit-exact under the pinned evaluation order.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        proof: ExactIdentityProofRef,
        mass_flow: MassFlowRate,
        density: Density,
        specific_enthalpy: SpecificEnergy,
        specific_internal_energy: SpecificEnergy,
        pressure: Pressure,
        volume_flow: VolumetricFlowRate,
    ) -> Result<Self, CoupleError> {
        if proof.identity != StreamIdentity::PressureWork {
            return Err(CoupleError::WrongStreamIdentityProof {
                expected: StreamIdentity::PressureWork,
                actual: proof.identity,
            });
        }
        for (field, value) in [
            ("crosswalk_mass_flow", mass_flow.value()),
            ("crosswalk_density", density.value()),
            ("crosswalk_specific_enthalpy", specific_enthalpy.value()),
            (
                "crosswalk_specific_internal_energy",
                specific_internal_energy.value(),
            ),
            ("crosswalk_pressure", pressure.value()),
            ("crosswalk_volume_flow", volume_flow.value()),
        ] {
            ensure_finite_stream_value(field, None, value)?;
        }
        if density.value() <= 0.0 {
            return Err(CoupleError::NonPositiveStreamDensity);
        }

        require_exact_stream_identity(
            StreamIdentity::PressureWork,
            "specific_enthalpy",
            specific_internal_energy.value() + pressure.value() / density.value(),
            specific_enthalpy.value(),
        )?;
        require_exact_stream_identity(
            StreamIdentity::PressureWork,
            "volume_flow",
            mass_flow.value() / density.value(),
            volume_flow.value(),
        )?;
        let enthalpy_gap_rate =
            mass_flow.value() * (specific_enthalpy.value() - specific_internal_energy.value());
        let pressure_work_rate = pressure.value() * volume_flow.value();
        require_exact_stream_identity(
            StreamIdentity::PressureWork,
            "pressure_work_rate",
            enthalpy_gap_rate,
            pressure_work_rate,
        )?;

        Ok(Self {
            proof,
            mass_flow,
            density,
            specific_enthalpy,
            specific_internal_energy,
            pressure,
            volume_flow,
            pressure_work_rate: Power::new(pressure_work_rate),
        })
    }

    /// Context-bound proof reference.
    #[must_use]
    pub fn proof(&self) -> &ExactIdentityProofRef {
        &self.proof
    }

    /// Signed mass rate used by the exact rate identity.
    #[must_use]
    pub fn mass_flow(&self) -> MassFlowRate {
        self.mass_flow
    }

    /// Specific internal energy in joules per kilogram.
    #[must_use]
    pub fn specific_internal_energy(&self) -> SpecificEnergy {
        self.specific_internal_energy
    }

    /// Exact pressure-work contribution in watts.
    #[must_use]
    pub fn pressure_work_rate(&self) -> Power {
        self.pressure_work_rate
    }

    /// Exact state terms retained by the crosswalk.
    #[must_use]
    pub fn state_terms(&self) -> (Density, SpecificEnergy, Pressure, VolumetricFlowRate) {
        (
            self.density,
            self.specific_enthalpy,
            self.pressure,
            self.volume_flow,
        )
    }
}

/// Internal-energy chart with explicit pressure plus deviatoric Cauchy work.
#[derive(Debug, Clone, PartialEq)]
pub struct InternalEnergyCauchyWorkChart {
    binding: StreamChartBinding,
    pressure_crosswalk: PressureWorkCrosswalk,
    kinematics: StreamKinematics,
    deviatoric_work: DeviatoricStressWork,
    chemical_energy: ChemicalEnergyAccounting,
}

impl InternalEnergyCauchyWorkChart {
    /// Admit an internal-energy chart only when its pressure-work and chemical
    /// ownership evidence bind the exact same stream context.
    ///
    /// # Errors
    ///
    /// [`CoupleError::StreamChartBindingMismatch`] for any context mismatch.
    pub fn try_new(
        binding: StreamChartBinding,
        pressure_crosswalk: PressureWorkCrosswalk,
        kinematics: StreamKinematics,
        deviatoric_work: DeviatoricStressWork,
        chemical_energy: ChemicalEnergyAccounting,
    ) -> Result<Self, CoupleError> {
        require_stream_binding(&binding, pressure_crosswalk.proof.binding())?;
        require_stream_binding(&binding, deviatoric_work.binding())?;
        require_stream_binding(&binding, chemical_energy.binding())?;
        Ok(Self {
            binding,
            pressure_crosswalk,
            kinematics,
            deviatoric_work,
            chemical_energy,
        })
    }

    fn energy_rate(&self) -> Power {
        let mass_flow = self.pressure_crosswalk.mass_flow.value();
        Power::new(
            mass_flow
                * (self.pressure_crosswalk.specific_internal_energy.value()
                    + self.kinematics.specific_transport_energy())
                + self.pressure_crosswalk.pressure_work_rate.value()
                + self.deviatoric_work.power_rate.value()
                + self.chemical_energy.explicit_power_rate.value(),
        )
    }
}

/// Specific thermodynamic potential selected from an exact Euler/Legendre
/// crosswalk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConjugatePotentialKind {
    /// Internal energy.
    InternalEnergy,
    /// Enthalpy.
    Enthalpy,
    /// Helmholtz free energy.
    Helmholtz,
    /// Gibbs free energy.
    Gibbs,
}

/// Exact mixture Euler identity and its three canonical Legendre transforms.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerLegendreCrosswalk {
    proof: ExactIdentityProofRef,
    internal_energy: SpecificEnergy,
    temperature_entropy_term: SpecificEnergy,
    pressure_volume_term: SpecificEnergy,
    chemical_potential_term: SpecificEnergy,
    enthalpy: SpecificEnergy,
    helmholtz: SpecificEnergy,
    gibbs: SpecificEnergy,
}

impl EulerLegendreCrosswalk {
    /// Recompute the specific-energy identities
    /// `u = T*s - p*v + sum(mu_i*n_i)`, `h = u + p*v`,
    /// `f = u - T*s`, and `g = h - T*s` exactly.
    ///
    /// # Errors
    ///
    /// Refuses the wrong proof kind, non-finite values, or any non-exact
    /// identity under the pinned evaluation order.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        proof: ExactIdentityProofRef,
        internal_energy: SpecificEnergy,
        temperature_entropy_term: SpecificEnergy,
        pressure_volume_term: SpecificEnergy,
        chemical_potential_term: SpecificEnergy,
        enthalpy: SpecificEnergy,
        helmholtz: SpecificEnergy,
        gibbs: SpecificEnergy,
    ) -> Result<Self, CoupleError> {
        if proof.identity != StreamIdentity::EulerLegendre {
            return Err(CoupleError::WrongStreamIdentityProof {
                expected: StreamIdentity::EulerLegendre,
                actual: proof.identity,
            });
        }
        for (field, value) in [
            ("euler_internal_energy", internal_energy.value()),
            (
                "euler_temperature_entropy_term",
                temperature_entropy_term.value(),
            ),
            ("euler_pressure_volume_term", pressure_volume_term.value()),
            (
                "euler_chemical_potential_term",
                chemical_potential_term.value(),
            ),
            ("euler_enthalpy", enthalpy.value()),
            ("euler_helmholtz", helmholtz.value()),
            ("euler_gibbs", gibbs.value()),
        ] {
            ensure_finite_stream_value(field, None, value)?;
        }
        require_exact_stream_identity(
            StreamIdentity::EulerLegendre,
            "mixture_euler",
            (temperature_entropy_term.value() - pressure_volume_term.value())
                + chemical_potential_term.value(),
            internal_energy.value(),
        )?;
        require_exact_stream_identity(
            StreamIdentity::EulerLegendre,
            "enthalpy_legendre",
            internal_energy.value() + pressure_volume_term.value(),
            enthalpy.value(),
        )?;
        require_exact_stream_identity(
            StreamIdentity::EulerLegendre,
            "helmholtz_legendre",
            internal_energy.value() - temperature_entropy_term.value(),
            helmholtz.value(),
        )?;
        require_exact_stream_identity(
            StreamIdentity::EulerLegendre,
            "gibbs_legendre",
            enthalpy.value() - temperature_entropy_term.value(),
            gibbs.value(),
        )?;
        Ok(Self {
            proof,
            internal_energy,
            temperature_entropy_term,
            pressure_volume_term,
            chemical_potential_term,
            enthalpy,
            helmholtz,
            gibbs,
        })
    }

    /// Context-bound proof reference.
    #[must_use]
    pub fn proof(&self) -> &ExactIdentityProofRef {
        &self.proof
    }

    /// Select one exactly crosswalked specific potential.
    #[must_use]
    pub fn potential(&self, kind: ConjugatePotentialKind) -> SpecificEnergy {
        match kind {
            ConjugatePotentialKind::InternalEnergy => self.internal_energy,
            ConjugatePotentialKind::Enthalpy => self.enthalpy,
            ConjugatePotentialKind::Helmholtz => self.helmholtz,
            ConjugatePotentialKind::Gibbs => self.gibbs,
        }
    }

    /// Reconstruct canonical specific enthalpy from the selected coordinate.
    ///
    /// A Legendre potential is a chart coordinate, not itself the transported
    /// moving-stream energy. The retained exact dual-pair terms translate each
    /// selection back to enthalpy before the stream rate is formed.
    #[must_use]
    pub fn enthalpy_from(&self, kind: ConjugatePotentialKind) -> SpecificEnergy {
        let value = match kind {
            ConjugatePotentialKind::InternalEnergy => {
                self.internal_energy.value() + self.pressure_volume_term.value()
            }
            ConjugatePotentialKind::Enthalpy => self.enthalpy.value(),
            ConjugatePotentialKind::Helmholtz => {
                (self.helmholtz.value() + self.temperature_entropy_term.value())
                    + self.pressure_volume_term.value()
            }
            ConjugatePotentialKind::Gibbs => {
                self.gibbs.value() + self.temperature_entropy_term.value()
            }
        };
        SpecificEnergy::new(value)
    }

    /// Euler dual-pair terms retained by the exact crosswalk.
    #[must_use]
    pub fn euler_terms(&self) -> (SpecificEnergy, SpecificEnergy, SpecificEnergy) {
        (
            self.temperature_entropy_term,
            self.pressure_volume_term,
            self.chemical_potential_term,
        )
    }
}

/// Energy chart using one potential from an exact Euler/Legendre family.
#[derive(Debug, Clone, PartialEq)]
pub struct ConjugatePotentialChart {
    binding: StreamChartBinding,
    crosswalk: EulerLegendreCrosswalk,
    selected_potential: ConjugatePotentialKind,
    kinematics: StreamKinematics,
    deviatoric_work: DeviatoricStressWork,
    chemical_energy: ChemicalEnergyAccounting,
}

impl ConjugatePotentialChart {
    /// Admit a conjugate-potential chart bound to one exact Euler/Legendre
    /// receipt. Its chemical term is already inside the Euler identity.
    ///
    /// # Errors
    ///
    /// Refuses context mismatch or a separate chemical-power term that would
    /// double count the Euler chemical-potential contribution.
    pub fn try_new(
        binding: StreamChartBinding,
        crosswalk: EulerLegendreCrosswalk,
        selected_potential: ConjugatePotentialKind,
        kinematics: StreamKinematics,
        deviatoric_work: DeviatoricStressWork,
        chemical_energy: ChemicalEnergyAccounting,
    ) -> Result<Self, CoupleError> {
        require_stream_binding(&binding, crosswalk.proof.binding())?;
        require_stream_binding(&binding, deviatoric_work.binding())?;
        require_stream_binding(&binding, chemical_energy.binding())?;
        require_species_constituent_axis(&binding)?;
        if chemical_energy.mode != ChemicalEnergyMode::IncludedInStatePotential {
            return Err(CoupleError::DoubleCountedChemicalEnergy);
        }
        require_exact_stream_identity(
            StreamIdentity::EulerLegendre,
            "selected_potential_to_enthalpy",
            crosswalk.enthalpy_from(selected_potential).value(),
            crosswalk.enthalpy.value(),
        )?;
        Ok(Self {
            binding,
            crosswalk,
            selected_potential,
            kinematics,
            deviatoric_work,
            chemical_energy,
        })
    }

    fn energy_rate(&self, mass_flow: MassFlowRate) -> Power {
        Power::new(
            mass_flow.value()
                * (self
                    .crosswalk
                    .enthalpy_from(self.selected_potential)
                    .value()
                    + self.kinematics.specific_transport_energy())
                + self.deviatoric_work.power_rate.value()
                + self.chemical_energy.explicit_power_rate.value(),
        )
    }
}

/// Exactly one selected energy-accounting chart for a stream boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum StreamEnergyChart {
    /// Canonical moving-stream enthalpy chart.
    MovingStreamEnthalpy(Box<MovingStreamEnthalpyChart>),
    /// Internal energy plus exact pressure and deviatoric Cauchy work.
    InternalEnergyCauchyWork(Box<InternalEnergyCauchyWorkChart>),
    /// Exact Euler/Legendre conjugate-potential chart.
    ConjugatePotential(Box<ConjugatePotentialChart>),
}

impl StreamEnergyChart {
    /// Complete context bound by the selected chart.
    #[must_use]
    pub fn binding(&self) -> &StreamChartBinding {
        match self {
            Self::MovingStreamEnthalpy(chart) => &chart.binding,
            Self::InternalEnergyCauchyWork(chart) => &chart.binding,
            Self::ConjugatePotential(chart) => &chart.binding,
        }
    }

    fn energy_rate(&self, mass_flow: MassFlowRate) -> Result<Power, CoupleError> {
        match self {
            Self::MovingStreamEnthalpy(chart) => Ok(chart.energy_rate(mass_flow)),
            Self::InternalEnergyCauchyWork(chart) => {
                require_exact_stream_identity(
                    StreamIdentity::PressureWork,
                    "stream_mass_flow_binding",
                    chart.pressure_crosswalk.mass_flow.value(),
                    mass_flow.value(),
                )?;
                Ok(chart.energy_rate())
            }
            Self::ConjugatePotential(chart) => Ok(chart.energy_rate(mass_flow)),
        }
    }
}

/// One admitted stream boundary carrying all coupled conserved rates together.
///
/// Unlike a scalar effort/flow [`PortSchema`], a stream is a bundled transport
/// object: mass, constituent amount, three momentum components, energy, and
/// entropy share one identity, basis/frame/orientation, clock, and energy chart.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamPort {
    version: u16,
    binding: StreamChartBinding,
    mass_flow: MassFlowRate,
    constituent_flows: Vec<StreamConstituentFlow>,
    momentum_flow: [Force; 3],
    energy_flow: Power,
    entropy_flow: EntropyFlowRate,
    energy_chart: StreamEnergyChart,
}

impl StreamPort {
    /// Admit a complete stream bundle and verify its selected chart reproduces
    /// the declared energy rate exactly.
    ///
    /// # Errors
    ///
    /// Refuses missing/duplicate constituents, non-finite rates, non-outward
    /// orientation, chart context mismatch, or an exact energy-rate mismatch.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        binding: StreamChartBinding,
        mass_flow: MassFlowRate,
        constituent_flows: impl IntoIterator<Item = StreamConstituentFlow>,
        momentum_flow: [Force; 3],
        energy_flow: Power,
        entropy_flow: EntropyFlowRate,
        energy_chart: StreamEnergyChart,
    ) -> Result<Self, CoupleError> {
        if binding.coordinates.orientation != PortOrientation::OutwardFromOwner {
            return Err(CoupleError::StreamPortRequiresOutwardOrientation {
                actual: binding.coordinates.orientation,
            });
        }
        require_stream_binding(&binding, energy_chart.binding())?;
        ensure_finite_stream_value("mass_flow", None, mass_flow.value())?;
        for (index, component) in momentum_flow.into_iter().enumerate() {
            ensure_finite_stream_value("momentum_flow", Some(index), component.value())?;
        }
        ensure_finite_stream_value("energy_flow", None, energy_flow.value())?;
        ensure_finite_stream_value("entropy_flow", None, entropy_flow.value())?;

        let mut constituent_flows: Vec<_> = constituent_flows.into_iter().collect();
        if constituent_flows.is_empty() {
            return Err(CoupleError::EmptyStreamConstituents);
        }
        constituent_flows.sort_by(|left, right| left.id.cmp(&right.id));
        for pair in constituent_flows.windows(2) {
            if pair[0].id == pair[1].id {
                return Err(CoupleError::DuplicateStreamConstituent {
                    id: pair[0].id.diagnostic_label(),
                });
            }
        }
        let actual_axis: Vec<_> = constituent_flows
            .iter()
            .map(|flow| flow.id.clone())
            .collect();
        if actual_axis != binding.constituent_axis {
            return Err(CoupleError::StreamConstituentAxisMismatch {
                expected: binding
                    .constituent_axis
                    .iter()
                    .map(StreamConstituentId::diagnostic_label)
                    .collect(),
                actual: actual_axis
                    .iter()
                    .map(StreamConstituentId::diagnostic_label)
                    .collect(),
            });
        }

        let accounted_energy = energy_chart.energy_rate(mass_flow)?;
        ensure_finite_stream_value("accounted_energy_flow", None, accounted_energy.value())?;
        if accounted_energy.value().to_bits() != energy_flow.value().to_bits() {
            return Err(CoupleError::StreamEnergyFlowMismatch {
                accounted_bits: accounted_energy.value().to_bits(),
                declared_bits: energy_flow.value().to_bits(),
            });
        }

        Ok(Self {
            version: STREAM_PORT_VERSION,
            binding,
            mass_flow,
            constituent_flows,
            momentum_flow,
            energy_flow,
            entropy_flow,
            energy_chart,
        })
    }

    /// Stream schema version.
    #[must_use]
    pub fn version(&self) -> u16 {
        self.version
    }

    /// Stable stream identifier and complete context.
    #[must_use]
    pub fn binding(&self) -> &StreamChartBinding {
        &self.binding
    }

    /// Signed mass rate in kilograms per second.
    #[must_use]
    pub fn mass_flow(&self) -> MassFlowRate {
        self.mass_flow
    }

    /// Canonically ordered species/element amount rates.
    #[must_use]
    pub fn constituent_flows(&self) -> &[StreamConstituentFlow] {
        &self.constituent_flows
    }

    /// Signed momentum-rate vector in newtons.
    #[must_use]
    pub fn momentum_flow(&self) -> [Force; 3] {
        self.momentum_flow
    }

    /// Signed energy rate in watts.
    #[must_use]
    pub fn energy_flow(&self) -> Power {
        self.energy_flow
    }

    /// Signed entropy rate in watts per kelvin.
    #[must_use]
    pub fn entropy_flow(&self) -> EntropyFlowRate {
        self.entropy_flow
    }

    /// Structurally unique selected energy chart.
    #[must_use]
    pub fn energy_chart(&self) -> &StreamEnergyChart {
        &self.energy_chart
    }

    /// Fixed bundled conservation roles in canonical enum order.
    #[must_use]
    pub fn conservation_roles(&self) -> &'static [ConservationRole] {
        &[
            ConservationRole::Energy,
            ConservationRole::Mass,
            ConservationRole::Amount,
            ConservationRole::LinearMomentum,
            ConservationRole::Entropy,
        ]
    }
}

fn ensure_finite_stream_value(
    field: &'static str,
    index: Option<usize>,
    value: f64,
) -> Result<(), CoupleError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(CoupleError::NonFiniteStreamValue { field, index })
    }
}

fn canonicalize_constituent_axis(
    constituent_axis: &mut [StreamConstituentId],
) -> Result<(), CoupleError> {
    if constituent_axis.is_empty() {
        return Err(CoupleError::EmptyStreamConstituents);
    }
    constituent_axis.sort();
    for pair in constituent_axis.windows(2) {
        if pair[0] == pair[1] {
            return Err(CoupleError::DuplicateStreamConstituent {
                id: pair[0].diagnostic_label(),
            });
        }
    }
    Ok(())
}

fn require_species_constituent_axis(binding: &StreamChartBinding) -> Result<(), CoupleError> {
    binding
        .constituent_axis
        .iter()
        .find(|constituent| !matches!(constituent, StreamConstituentId::Species(_)))
        .map_or(Ok(()), |constituent| {
            Err(CoupleError::StreamChemicalEnergyRequiresSpeciesAxis {
                constituent: constituent.diagnostic_label(),
            })
        })
}

fn require_stream_binding(
    expected: &StreamChartBinding,
    actual: &StreamChartBinding,
) -> Result<(), CoupleError> {
    expected.first_mismatch(actual).map_or(Ok(()), |field| {
        Err(CoupleError::StreamChartBindingMismatch { field })
    })
}

fn require_identity_proof(
    binding: &StreamChartBinding,
    proof: &ExactIdentityProofRef,
    expected: StreamIdentity,
) -> Result<(), CoupleError> {
    if proof.identity != expected {
        return Err(CoupleError::WrongStreamIdentityProof {
            expected,
            actual: proof.identity,
        });
    }
    require_stream_binding(binding, proof.binding())
}

fn require_stable_id_match(
    field: &'static str,
    expected: &StableId,
    actual: &StableId,
) -> Result<(), CoupleError> {
    if expected == actual {
        Ok(())
    } else {
        Err(CoupleError::StreamChartBindingMismatch { field })
    }
}

fn require_exact_stream_identity(
    identity: StreamIdentity,
    check: &'static str,
    expected: f64,
    actual: f64,
) -> Result<(), CoupleError> {
    ensure_finite_stream_value(check, None, expected)?;
    ensure_finite_stream_value(check, None, actual)?;
    if expected.to_bits() == actual.to_bits() {
        Ok(())
    } else {
        Err(CoupleError::StreamIdentityMismatch { identity, check })
    }
}

/// A raw scalar power-coordinate container retained for migration.
///
/// Construction performs no schema admission. New rotational, electrical,
/// magnetic, and chemical kinds must use [`PortSchema`] before composition;
/// [`Port::conjugate_to`] and [`interconnect`] deliberately refuse them here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Port {
    /// The raw effort coordinate.
    pub effort: f64,
    /// The raw flow coordinate.
    pub flow: f64,
    /// The declared physical kind; not schema-admitted by this container.
    pub kind: PortKind,
}

impl Port {
    /// Construct an unadmitted raw scalar container.
    #[must_use]
    pub fn new(effort: f64, flow: f64, kind: PortKind) -> Self {
        Self { effort, flow, kind }
    }

    /// Raw effort coordinate.
    #[must_use]
    pub fn effort(&self) -> f64 {
        self.effort
    }

    /// Raw flow coordinate.
    #[must_use]
    pub fn flow(&self) -> f64 {
        self.flow
    }

    /// Declared raw physical kind.
    #[must_use]
    pub fn kind(&self) -> PortKind {
        self.kind
    }

    /// The instantaneous power `effort × flow`.
    #[must_use]
    pub fn power(&self) -> f64 {
        self.effort * self.flow
    }

    /// Whether two raw ports are composable by the retained three-kind scalar
    /// migration oracle.
    #[must_use]
    pub fn conjugate_to(&self, other: &Port) -> bool {
        self.kind.is_legacy_scalar_seed() && self.kind == other.kind
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

/// Closed logical interval over which already-integrated balances are audited.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AccountingWindowInterval {
    id: StableId,
    start: PortTimestamp,
    end: PortTimestamp,
}

impl AccountingWindowInterval {
    /// Admit one non-empty interval in a single named logical clock.
    ///
    /// Logical tick differences are not interpreted as seconds. Boundary
    /// contributions require separate integration evidence.
    ///
    /// # Errors
    ///
    /// Refuses different clock domains or `end.tick <= start.tick`.
    pub fn try_new(
        id: StableId,
        start: PortTimestamp,
        end: PortTimestamp,
    ) -> Result<Self, CoupleError> {
        if start.clock != end.clock {
            return Err(CoupleError::AccountingWindowContextMismatch {
                field: "interval_clock",
            });
        }
        if end.tick <= start.tick {
            return Err(CoupleError::InvalidAccountingWindowInterval {
                start_tick: start.tick,
                end_tick: end.tick,
            });
        }
        Ok(Self { id, start, end })
    }

    /// Stable window identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// Inclusive initial logical timestamp.
    #[must_use]
    pub fn start(&self) -> &PortTimestamp {
        &self.start
    }

    /// Inclusive final logical timestamp.
    #[must_use]
    pub fn end(&self) -> &PortTimestamp {
        &self.end
    }

    fn contains(&self, timestamp: &PortTimestamp) -> bool {
        timestamp.clock == self.start.clock
            && (self.start.tick..=self.end.tick).contains(&timestamp.tick)
    }
}

/// Semantic role and subject bound by one retained window-evidence statement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowEvidenceRole {
    /// Evidence that one exact declared boundary/source manifest is complete.
    ManifestClosure {
        /// Canonically ordered complete row subjects covered by the statement.
        rows: Vec<WindowManifestRowSubject>,
    },
    /// Evidence binding one manifest row to its local port.
    ManifestBinding {
        /// Expected contribution row.
        contribution_id: StableId,
        /// Local port named by the row.
        local_port_id: StableId,
    },
    /// Evidence for one exact stored-state endpoint.
    Inventory {
        /// Initial or final state covered by the statement.
        endpoint: WindowEndpoint,
        /// Common stored/transfer energy-reference contract.
        energy_reference: StableId,
    },
    /// Evidence for physical-time integration and sign normalization.
    IntegratedTransfer {
        /// Filled contribution row.
        contribution_id: StableId,
        /// Local port whose values were integrated.
        local_port_id: StableId,
        /// Common stored/transfer energy-reference contract.
        energy_reference: StableId,
    },
    /// Boundary-temperature evidence for one contribution.
    BoundaryTemperature {
        /// Contribution whose heat/radiation entropy uses the temperature.
        contribution_id: StableId,
    },
    /// Species/direct-element projection evidence for one exact window axis.
    ElementProjection {
        /// Canonically ordered element axis named by the statement.
        elements: Vec<ElementId>,
    },
    /// Direct-coulomb, species-charge, or neutrality projection evidence.
    ChargeProjection {
        /// Exact charge policy and basis named by the statement.
        subject: WindowChargeEvidenceSubject,
    },
}

impl WindowEvidenceRole {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            Self::ManifestClosure { .. } => "manifest_closure",
            Self::ManifestBinding { .. } => "manifest_binding",
            Self::Inventory { .. } => "inventory",
            Self::IntegratedTransfer { .. } => "integrated_transfer",
            Self::BoundaryTemperature { .. } => "boundary_temperature",
            Self::ElementProjection { .. } => "element_projection",
            Self::ChargeProjection { .. } => "charge_projection",
        }
    }

    /// Bind closure evidence to the complete canonical manifest structure.
    #[must_use]
    pub fn manifest_closure(manifest: &[WindowManifestEntry]) -> Self {
        let mut rows: Vec<_> = manifest
            .iter()
            .map(WindowManifestRowSubject::from_entry)
            .collect();
        rows.sort_by(|left, right| left.contribution_id.cmp(&right.contribution_id));
        Self::ManifestClosure { rows }
    }
}

/// Exact charge-policy subject named by a retained evidence statement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowChargeEvidenceSubject {
    /// Direct coulomb values under this basis/sign convention.
    DirectCoulomb {
        /// Charge basis and sign convention.
        basis: StableId,
    },
    /// Species amounts projected to coulombs under this complete basis.
    SpeciesProjection {
        /// Species-charge and molar-charge convention.
        basis: StableId,
    },
    /// Exact zero charge under a retained neutrality claim.
    ProvenNeutral,
}

/// Retained reference to evidence covering one exact accounting interval.
///
/// `fs-couple` checks the interval binding but does not execute the named
/// verifier or validate the referenced artifact contents.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WindowEvidenceRef {
    receipt_id: StableId,
    verifier_id: StableId,
    statement_digest: StableId,
    role: WindowEvidenceRole,
    interval: Box<AccountingWindowInterval>,
}

impl WindowEvidenceRef {
    /// Name one retained interval-bound evidence artifact.
    #[must_use]
    pub fn new(
        receipt_id: StableId,
        verifier_id: StableId,
        statement_digest: StableId,
        role: WindowEvidenceRole,
        interval: AccountingWindowInterval,
    ) -> Self {
        Self {
            receipt_id,
            verifier_id,
            statement_digest,
            role,
            interval: Box::new(interval),
        }
    }

    /// Evidence receipt identifier.
    #[must_use]
    pub fn receipt_id(&self) -> &StableId {
        &self.receipt_id
    }

    /// Verifier implementation named by the receipt.
    #[must_use]
    pub fn verifier_id(&self) -> &StableId {
        &self.verifier_id
    }

    /// Digest of the retained evidence statement.
    #[must_use]
    pub fn statement_digest(&self) -> &StableId {
        &self.statement_digest
    }

    /// Exact semantic role and subject named by the statement.
    #[must_use]
    pub fn role(&self) -> &WindowEvidenceRole {
        &self.role
    }

    /// Exact logical interval covered by the statement.
    #[must_use]
    pub fn interval(&self) -> &AccountingWindowInterval {
        self.interval.as_ref()
    }
}

/// One canonically identified elemental amount in an inventory or transfer.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowElementAmount {
    element: ElementId,
    amount: Amount,
}

impl WindowElementAmount {
    /// Construct one finite elemental amount.
    ///
    /// # Errors
    ///
    /// [`CoupleError::NonFiniteAccountingWindowValue`] for a non-finite amount.
    pub fn try_new(element: ElementId, amount: Amount) -> Result<Self, CoupleError> {
        ensure_finite_window_value("element_amount", Some(element.as_str()), amount.value())?;
        Ok(Self { element, amount })
    }

    /// Element identifier.
    #[must_use]
    pub fn element(&self) -> &ElementId {
        &self.element
    }

    /// Signed amount in moles.
    #[must_use]
    pub fn amount(&self) -> Amount {
        self.amount
    }
}

/// Integrated extensive values used by one side of a window balance.
///
/// Boundary transfers use the single convention "positive into the audited
/// system". Inventory snapshots use ordinary stored values.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowBalance {
    energy: Energy,
    mass: Mass,
    element_amounts: Vec<WindowElementAmount>,
    charge: ElectricCharge,
    entropy: Entropy,
}

impl WindowBalance {
    /// Admit finite values and canonicalize the element axis.
    ///
    /// # Errors
    ///
    /// Refuses non-finite values or duplicate element identifiers.
    pub fn try_new(
        energy: Energy,
        mass: Mass,
        element_amounts: impl IntoIterator<Item = WindowElementAmount>,
        charge: ElectricCharge,
        entropy: Entropy,
    ) -> Result<Self, CoupleError> {
        ensure_finite_window_value("energy", None, energy.value())?;
        ensure_finite_window_value("mass", None, mass.value())?;
        ensure_finite_window_value("charge", None, charge.value())?;
        ensure_finite_window_value("entropy", None, entropy.value())?;
        let mut element_amounts: Vec<_> = element_amounts.into_iter().collect();
        element_amounts.sort_by(|left, right| left.element.cmp(&right.element));
        for pair in element_amounts.windows(2) {
            if pair[0].element == pair[1].element {
                return Err(CoupleError::DuplicateAccountingWindowElement {
                    element: pair[0].element.as_str().to_string(),
                });
            }
        }
        Ok(Self {
            energy,
            mass,
            element_amounts,
            charge,
            entropy,
        })
    }

    /// Stored or integrated energy in joules.
    #[must_use]
    pub fn energy(&self) -> Energy {
        self.energy
    }

    /// Stored or integrated mass in kilograms.
    #[must_use]
    pub fn mass(&self) -> Mass {
        self.mass
    }

    /// Canonically ordered elemental amounts.
    #[must_use]
    pub fn element_amounts(&self) -> &[WindowElementAmount] {
        &self.element_amounts
    }

    /// Stored or integrated electric charge in coulombs.
    #[must_use]
    pub fn charge(&self) -> ElectricCharge {
        self.charge
    }

    /// Stored or integrated entropy in joules per kelvin.
    #[must_use]
    pub fn entropy(&self) -> Entropy {
        self.entropy
    }
}

/// Boundary-temperature evidence used to interpret heat/radiation entropy.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryTemperatureReference {
    /// No heat- or radiation-mediated entropy is present in this contribution.
    NotApplicable {
        /// Retained rationale for the no-claim boundary.
        rationale: StableId,
    },
    /// One strictly positive boundary temperature was used for integration.
    Constant {
        /// Absolute temperature in kelvin.
        temperature: Temperature,
        /// Evidence for the temperature and its applicability to this row.
        evidence: WindowEvidenceRef,
    },
    /// A caller-integrated temperature profile is retained by identity.
    Profile {
        /// Temperature-profile artifact identifier.
        profile: StableId,
        /// Evidence for the profile and its applicability to this row.
        evidence: WindowEvidenceRef,
    },
}

impl BoundaryTemperatureReference {
    /// Declare that no boundary-temperature interpretation is applicable.
    #[must_use]
    pub fn not_applicable(rationale: StableId) -> Self {
        Self::NotApplicable { rationale }
    }

    /// Admit one finite, strictly positive constant boundary temperature.
    ///
    /// # Errors
    ///
    /// Returns [`CoupleError::NonPositiveBoundaryTemperature`] for NaN,
    /// infinity, zero, or a negative absolute temperature.
    pub fn try_constant(
        temperature: Temperature,
        evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(CoupleError::NonPositiveBoundaryTemperature);
        }
        Ok(Self::Constant {
            temperature,
            evidence,
        })
    }

    /// Retain a caller-integrated boundary-temperature profile by identity.
    #[must_use]
    pub fn profile(profile: StableId, evidence: WindowEvidenceRef) -> Self {
        Self::Profile { profile, evidence }
    }

    fn evidence(&self) -> Option<&WindowEvidenceRef> {
        match self {
            Self::NotApplicable { .. } => None,
            Self::Constant { evidence, .. } | Self::Profile { evidence, .. } => Some(evidence),
        }
    }
}

/// Explicit integrated boundary-entropy mechanisms, all positive into system.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryEntropyBreakdown {
    advected: Entropy,
    diffusive_chemical: Entropy,
    heat: Entropy,
    radiation: Entropy,
    temperature_reference: BoundaryTemperatureReference,
}

impl BoundaryEntropyBreakdown {
    /// Admit all four mandatory signed entropy-transfer slots.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite term or a non-finite fixed-order total.
    pub fn try_new(
        advected: Entropy,
        diffusive_chemical: Entropy,
        heat: Entropy,
        radiation: Entropy,
        temperature_reference: BoundaryTemperatureReference,
    ) -> Result<Self, CoupleError> {
        for (field, value) in [
            ("entropy_advected", advected.value()),
            ("entropy_diffusive_chemical", diffusive_chemical.value()),
            ("entropy_heat", heat.value()),
            ("entropy_radiation", radiation.value()),
        ] {
            ensure_finite_window_value(field, None, value)?;
        }
        let breakdown = Self {
            advected,
            diffusive_chemical,
            heat,
            radiation,
            temperature_reference,
        };
        if (heat.value().abs() > 0.0 || radiation.value().abs() > 0.0)
            && matches!(
                &breakdown.temperature_reference,
                BoundaryTemperatureReference::NotApplicable { .. }
            )
        {
            return Err(CoupleError::AccountingWindowTemperatureReferenceRequired);
        }
        ensure_finite_window_value("entropy_boundary_total", None, breakdown.total().value())?;
        Ok(breakdown)
    }

    /// Advected entropy transfer.
    #[must_use]
    pub fn advected(&self) -> Entropy {
        self.advected
    }

    /// Diffusive and chemical-potential entropy transfer.
    #[must_use]
    pub fn diffusive_chemical(&self) -> Entropy {
        self.diffusive_chemical
    }

    /// Heat-mediated entropy transfer under explicit boundary temperatures.
    #[must_use]
    pub fn heat(&self) -> Entropy {
        self.heat
    }

    /// Radiation-mediated entropy transfer.
    #[must_use]
    pub fn radiation(&self) -> Entropy {
        self.radiation
    }

    /// Boundary-temperature evidence used for heat/radiation interpretation.
    #[must_use]
    pub fn temperature_reference(&self) -> &BoundaryTemperatureReference {
        &self.temperature_reference
    }

    /// Fixed-order total boundary entropy transfer.
    #[must_use]
    pub fn total(&self) -> Entropy {
        Entropy::new(
            ((self.advected.value() + self.diffusive_chemical.value()) + self.heat.value())
                + self.radiation.value(),
        )
    }
}

/// Whether and how elemental conservation is claimed for a window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowElementSchema {
    /// No elemental balance is applicable; the rationale is retained.
    NotApplicable {
        /// Evidence/rationale identifier explaining the non-chemical window.
        rationale: StableId,
    },
    /// Exact named element axis is audited after an external projection.
    Audited {
        /// Canonically ordered element identifiers.
        elements: Vec<ElementId>,
        /// Evidence binding species/direct-element projection to this axis.
        projection_evidence: WindowEvidenceRef,
    },
}

impl WindowElementSchema {
    /// Declare a non-chemical window explicitly.
    #[must_use]
    pub fn not_applicable(rationale: StableId) -> Self {
        Self::NotApplicable { rationale }
    }

    /// Admit a non-empty duplicate-free elemental audit axis.
    ///
    /// # Errors
    ///
    /// Refuses an empty or duplicate element axis.
    pub fn try_audited(
        elements: impl IntoIterator<Item = ElementId>,
        projection_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        let mut elements: Vec<_> = elements.into_iter().collect();
        if elements.is_empty() {
            return Err(CoupleError::EmptyAccountingWindowElementAxis);
        }
        elements.sort();
        for pair in elements.windows(2) {
            if pair[0] == pair[1] {
                return Err(CoupleError::DuplicateAccountingWindowElement {
                    element: pair[0].as_str().to_string(),
                });
            }
        }
        require_window_evidence_role(
            &projection_evidence,
            &WindowEvidenceRole::ElementProjection {
                elements: elements.clone(),
            },
        )?;
        Ok(Self::Audited {
            elements,
            projection_evidence,
        })
    }

    /// Element axis, empty only for an explicitly non-chemical window.
    #[must_use]
    pub fn elements(&self) -> &[ElementId] {
        match self {
            Self::NotApplicable { .. } => &[],
            Self::Audited { elements, .. } => elements,
        }
    }

    /// Projection receipt for an audited elemental axis.
    #[must_use]
    pub fn projection_evidence(&self) -> Option<&WindowEvidenceRef> {
        match self {
            Self::NotApplicable { .. } => None,
            Self::Audited {
                projection_evidence,
                ..
            } => Some(projection_evidence),
        }
    }
}

/// Evidence policy governing the window's direct coulomb balance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowChargeSchema {
    /// Charge is supplied directly in coulombs under a named basis/convention.
    DirectCoulomb {
        /// Charge basis and sign convention.
        basis: StableId,
        /// Evidence binding the direct values to the declared convention.
        projection_evidence: WindowEvidenceRef,
    },
    /// Species amounts were projected through an explicit charge/Faraday map.
    SpeciesProjection {
        /// Species-charge basis and molar-charge convention.
        basis: StableId,
        /// Evidence for the species-to-coulomb projection.
        projection_evidence: WindowEvidenceRef,
    },
    /// Every charge value must be exactly zero under retained neutrality proof.
    ProvenNeutral {
        /// Evidence for the bound neutrality claim.
        projection_evidence: WindowEvidenceRef,
    },
}

impl WindowChargeSchema {
    /// Evidence governing the charge values admitted by this schema.
    #[must_use]
    pub fn projection_evidence(&self) -> &WindowEvidenceRef {
        match self {
            Self::DirectCoulomb {
                projection_evidence,
                ..
            }
            | Self::SpeciesProjection {
                projection_evidence,
                ..
            }
            | Self::ProvenNeutral {
                projection_evidence,
            } => projection_evidence,
        }
    }

    fn evidence_role(&self) -> WindowEvidenceRole {
        let subject = match self {
            Self::DirectCoulomb { basis, .. } => WindowChargeEvidenceSubject::DirectCoulomb {
                basis: basis.clone(),
            },
            Self::SpeciesProjection { basis, .. } => {
                WindowChargeEvidenceSubject::SpeciesProjection {
                    basis: basis.clone(),
                }
            }
            Self::ProvenNeutral { .. } => WindowChargeEvidenceSubject::ProvenNeutral,
        };
        WindowEvidenceRole::ChargeProjection { subject }
    }
}

/// Typed acceptance tolerances for a physical accounting-window audit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowAuditTolerances {
    energy: Energy,
    mass: Mass,
    amount: Amount,
    charge: ElectricCharge,
    entropy: Entropy,
}

impl WindowAuditTolerances {
    /// Admit finite nonnegative tolerances for every audited balance.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite or negative tolerance.
    pub fn try_new(
        energy: Energy,
        mass: Mass,
        amount: Amount,
        charge: ElectricCharge,
        entropy: Entropy,
    ) -> Result<Self, CoupleError> {
        for (field, value) in [
            ("energy_tolerance", energy.value()),
            ("mass_tolerance", mass.value()),
            ("amount_tolerance", amount.value()),
            ("charge_tolerance", charge.value()),
            ("entropy_tolerance", entropy.value()),
        ] {
            ensure_nonnegative_window_value(field, value)?;
        }
        Ok(Self {
            energy,
            mass,
            amount,
            charge,
            entropy,
        })
    }

    /// First-law tolerance in joules.
    #[must_use]
    pub fn energy(self) -> Energy {
        self.energy
    }

    /// Mass-balance tolerance in kilograms.
    #[must_use]
    pub fn mass(self) -> Mass {
        self.mass
    }

    /// Per-element tolerance in moles.
    #[must_use]
    pub fn amount(self) -> Amount {
        self.amount
    }

    /// Charge-balance tolerance in coulombs.
    #[must_use]
    pub fn charge(self) -> ElectricCharge {
        self.charge
    }

    /// One-sided entropy-generation tolerance in joules per kelvin.
    #[must_use]
    pub fn entropy(self) -> Entropy {
        self.entropy
    }
}

/// Explicit environment for the optional Gouy-Stodola exergy diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub struct ExergyEnvironment {
    id: StableId,
    temperature: Temperature,
}

impl ExergyEnvironment {
    /// Admit a finite, strictly positive reference temperature.
    ///
    /// # Errors
    ///
    /// [`CoupleError::NonPositiveExergyTemperature`] otherwise.
    pub fn try_new(id: StableId, temperature: Temperature) -> Result<Self, CoupleError> {
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(CoupleError::NonPositiveExergyTemperature);
        }
        Ok(Self { id, temperature })
    }

    /// Environment identifier.
    #[must_use]
    pub fn id(&self) -> &StableId {
        &self.id
    }

    /// Reference temperature.
    #[must_use]
    pub fn temperature(&self) -> Temperature {
        self.temperature
    }
}

/// Topological role of one expected window contribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowExchangeRole {
    /// One side of an internal two-port transfer; both sides share an exchange.
    InternalTransfer,
    /// Exchange with an explicitly external reservoir.
    ExternalExchange,
    /// Source term included inside the audited model.
    IncludedSource,
}

/// Physical vocabulary carried by a window-manifest port.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowPortKind {
    /// Ordinary power-conjugate port.
    Power(PortKind),
    /// Multi-balance stream bundle.
    Stream,
}

/// Complete structural subject of one manifest row, excluding its evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WindowManifestRowSubject {
    contribution_id: StableId,
    exchange_id: StableId,
    local_port_id: StableId,
    counterparty_id: StableId,
    role: WindowExchangeRole,
    port_kind: WindowPortKind,
    stream_binding: Option<Box<StreamChartBinding>>,
    timestamp: PortTimestamp,
    boundary: Option<AccountingBoundary>,
}

impl WindowManifestRowSubject {
    fn from_entry(entry: &WindowManifestEntry) -> Self {
        Self {
            contribution_id: entry.contribution_id.clone(),
            exchange_id: entry.exchange_id.clone(),
            local_port_id: entry.local_port_id.clone(),
            counterparty_id: entry.counterparty_id.clone(),
            role: entry.role,
            port_kind: entry.port_kind,
            stream_binding: entry.stream_binding.clone(),
            timestamp: entry.timestamp.clone(),
            boundary: entry.boundary.clone(),
        }
    }
}

/// One expected contribution in a complete accounting-window manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WindowManifestEntry {
    contribution_id: StableId,
    exchange_id: StableId,
    local_port_id: StableId,
    counterparty_id: StableId,
    role: WindowExchangeRole,
    port_kind: WindowPortKind,
    stream_binding: Option<Box<StreamChartBinding>>,
    timestamp: PortTimestamp,
    boundary: Option<AccountingBoundary>,
    binding_evidence: WindowEvidenceRef,
}

impl WindowManifestEntry {
    /// Build both sides of one admitted internal power-port pair.
    ///
    /// # Errors
    ///
    /// Refuses schemas that are not conjugate under the existing junction
    /// admission rules or evidence bound to a different window.
    pub fn internal_power_pair(
        exchange_id: StableId,
        contribution_a_id: StableId,
        port_a: &PortSchema,
        contribution_b_id: StableId,
        port_b: &PortSchema,
        binding_evidence_a: WindowEvidenceRef,
        binding_evidence_b: WindowEvidenceRef,
    ) -> Result<[Self; 2], CoupleError> {
        if let Some(field) = port_a.first_conjugacy_mismatch(port_b) {
            return Err(CoupleError::IncompatiblePortSchemas {
                a: port_a.id.as_str().to_string(),
                b: port_b.id.as_str().to_string(),
                field,
            });
        }
        let first = Self::internal(
            contribution_a_id,
            exchange_id.clone(),
            port_a.id.clone(),
            port_b.id.clone(),
            WindowPortKind::Power(port_a.kind),
            None,
            port_a.timestamp.clone(),
            binding_evidence_a,
        )?;
        let second = Self::internal(
            contribution_b_id,
            exchange_id,
            port_b.id.clone(),
            port_a.id.clone(),
            WindowPortKind::Power(port_b.kind),
            None,
            port_b.timestamp.clone(),
            binding_evidence_b,
        )?;
        Ok([first, second])
    }

    /// Build both sides of one context-compatible internal stream pair.
    ///
    /// # Errors
    ///
    /// Refuses reused IDs or mismatched state, constituent, frame, clock,
    /// gravity, or stress-work context.
    pub fn internal_stream_pair(
        exchange_id: StableId,
        contribution_a_id: StableId,
        port_a: &StreamPort,
        contribution_b_id: StableId,
        port_b: &StreamPort,
        binding_evidence_a: WindowEvidenceRef,
        binding_evidence_b: WindowEvidenceRef,
    ) -> Result<[Self; 2], CoupleError> {
        require_internal_stream_pair(port_a.binding(), port_b.binding())?;
        let first = Self::internal(
            contribution_a_id,
            exchange_id.clone(),
            port_a.binding.port_id.clone(),
            port_b.binding.port_id.clone(),
            WindowPortKind::Stream,
            Some(Box::new(port_a.binding.clone())),
            port_a.binding.timestamp.clone(),
            binding_evidence_a,
        )?;
        let second = Self::internal(
            contribution_b_id,
            exchange_id,
            port_b.binding.port_id.clone(),
            port_a.binding.port_id.clone(),
            WindowPortKind::Stream,
            Some(Box::new(port_b.binding.clone())),
            port_b.binding.timestamp.clone(),
            binding_evidence_b,
        )?;
        Ok([first, second])
    }

    /// Bind an ordinary port to an external accounting boundary.
    ///
    /// # Errors
    ///
    /// Refuses the wrong boundary treatment, coordinate mismatch, non-outward
    /// orientation, or identity aliasing.
    #[allow(clippy::too_many_arguments)]
    pub fn external_power_port(
        contribution_id: StableId,
        exchange_id: StableId,
        counterparty_id: StableId,
        port: &PortSchema,
        boundary: AccountingBoundary,
        binding_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        Self::external(
            contribution_id,
            exchange_id,
            port.id.clone(),
            counterparty_id,
            WindowPortKind::Power(port.kind),
            None,
            port.coordinates.clone(),
            port.timestamp.clone(),
            boundary,
            binding_evidence,
        )
    }

    /// Bind a stream port to an external accounting boundary.
    ///
    /// # Errors
    ///
    /// Refuses the wrong boundary treatment, coordinate mismatch, non-outward
    /// orientation, or identity aliasing.
    #[allow(clippy::too_many_arguments)]
    pub fn external_stream_port(
        contribution_id: StableId,
        exchange_id: StableId,
        counterparty_id: StableId,
        port: &StreamPort,
        boundary: AccountingBoundary,
        binding_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        Self::external(
            contribution_id,
            exchange_id,
            port.binding.port_id.clone(),
            counterparty_id,
            WindowPortKind::Stream,
            Some(Box::new(port.binding.clone())),
            port.binding.coordinates.clone(),
            port.binding.timestamp.clone(),
            boundary,
            binding_evidence,
        )
    }

    /// Bind an admitted included source and its explicit accounting boundary.
    ///
    /// # Errors
    ///
    /// Refuses a boundary not declared as an included source term.
    pub fn included_source(
        contribution_id: StableId,
        source: &SourceOrReservoir,
        binding_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        if source.boundary.treatment != BoundaryTreatment::IncludedSourceTerm {
            return Err(CoupleError::AccountingWindowBoundaryTreatmentMismatch {
                boundary: source.boundary.id.as_str().to_string(),
                expected: BoundaryTreatment::IncludedSourceTerm,
                actual: source.boundary.treatment,
            });
        }
        Self::from_parts(
            contribution_id,
            source.id.clone(),
            source.port.id.clone(),
            source.boundary.id.clone(),
            WindowExchangeRole::IncludedSource,
            WindowPortKind::Power(source.port.kind),
            None,
            source.port.timestamp.clone(),
            Some(source.boundary.clone()),
            binding_evidence,
        )
    }

    /// Bind an admitted external source/reservoir exchange.
    ///
    /// # Errors
    ///
    /// Refuses a boundary not declared as an outward external reservoir exchange.
    pub fn external_source(
        contribution_id: StableId,
        source: &SourceOrReservoir,
        binding_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        if source.boundary.treatment != BoundaryTreatment::ExternalReservoirExchange {
            return Err(CoupleError::AccountingWindowBoundaryTreatmentMismatch {
                boundary: source.boundary.id.as_str().to_string(),
                expected: BoundaryTreatment::ExternalReservoirExchange,
                actual: source.boundary.treatment,
            });
        }
        if source.boundary.coordinates.orientation != PortOrientation::OutwardFromOwner {
            return Err(CoupleError::AccountingWindowBoundaryRequiresOutward {
                boundary: source.boundary.id.as_str().to_string(),
            });
        }
        Self::from_parts(
            contribution_id,
            source.id.clone(),
            source.port.id.clone(),
            source.boundary.id.clone(),
            WindowExchangeRole::ExternalExchange,
            WindowPortKind::Power(source.port.kind),
            None,
            source.port.timestamp.clone(),
            Some(source.boundary.clone()),
            binding_evidence,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn internal(
        contribution_id: StableId,
        exchange_id: StableId,
        local_port_id: StableId,
        counterparty_id: StableId,
        port_kind: WindowPortKind,
        stream_binding: Option<Box<StreamChartBinding>>,
        timestamp: PortTimestamp,
        binding_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        Self::from_parts(
            contribution_id,
            exchange_id,
            local_port_id,
            counterparty_id,
            WindowExchangeRole::InternalTransfer,
            port_kind,
            stream_binding,
            timestamp,
            None,
            binding_evidence,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn external(
        contribution_id: StableId,
        exchange_id: StableId,
        local_port_id: StableId,
        counterparty_id: StableId,
        port_kind: WindowPortKind,
        stream_binding: Option<Box<StreamChartBinding>>,
        coordinates: CoordinateBinding,
        timestamp: PortTimestamp,
        boundary: AccountingBoundary,
        binding_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        if boundary.treatment != BoundaryTreatment::ExternalReservoirExchange {
            return Err(CoupleError::AccountingWindowBoundaryTreatmentMismatch {
                boundary: boundary.id.as_str().to_string(),
                expected: BoundaryTreatment::ExternalReservoirExchange,
                actual: boundary.treatment,
            });
        }
        if coordinates != boundary.coordinates {
            return Err(CoupleError::AccountingBoundaryMismatch {
                port: local_port_id.as_str().to_string(),
                boundary: boundary.id.as_str().to_string(),
            });
        }
        if coordinates.orientation != PortOrientation::OutwardFromOwner {
            return Err(CoupleError::AccountingWindowBoundaryRequiresOutward {
                boundary: boundary.id.as_str().to_string(),
            });
        }
        Self::from_parts(
            contribution_id,
            exchange_id,
            local_port_id,
            counterparty_id,
            WindowExchangeRole::ExternalExchange,
            port_kind,
            stream_binding,
            timestamp,
            Some(boundary),
            binding_evidence,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_parts(
        contribution_id: StableId,
        exchange_id: StableId,
        local_port_id: StableId,
        counterparty_id: StableId,
        role: WindowExchangeRole,
        port_kind: WindowPortKind,
        stream_binding: Option<Box<StreamChartBinding>>,
        timestamp: PortTimestamp,
        boundary: Option<AccountingBoundary>,
        binding_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        for (left, right) in [
            (&contribution_id, &exchange_id),
            (&contribution_id, &local_port_id),
            (&contribution_id, &counterparty_id),
            (&exchange_id, &local_port_id),
            (&exchange_id, &counterparty_id),
            (&local_port_id, &counterparty_id),
        ] {
            if left == right {
                return Err(CoupleError::DuplicateIdentity {
                    id: left.as_str().to_string(),
                });
            }
        }
        if let Some(boundary) = &boundary {
            for id in [&contribution_id, &exchange_id, &local_port_id] {
                if id == &boundary.id {
                    return Err(CoupleError::DuplicateIdentity {
                        id: boundary.id.as_str().to_string(),
                    });
                }
            }
        }
        if matches!(port_kind, WindowPortKind::Stream) != stream_binding.is_some() {
            return Err(CoupleError::AccountingWindowContextMismatch {
                field: "manifest_stream_binding",
            });
        }
        require_window_evidence_role(
            &binding_evidence,
            &WindowEvidenceRole::ManifestBinding {
                contribution_id: contribution_id.clone(),
                local_port_id: local_port_id.clone(),
            },
        )?;
        Ok(Self {
            contribution_id,
            exchange_id,
            local_port_id,
            counterparty_id,
            role,
            port_kind,
            stream_binding,
            timestamp,
            boundary,
            binding_evidence,
        })
    }

    /// Expected contribution identifier.
    #[must_use]
    pub fn contribution_id(&self) -> &StableId {
        &self.contribution_id
    }

    /// Pair/relation/exchange group used for localization.
    #[must_use]
    pub fn exchange_id(&self) -> &StableId {
        &self.exchange_id
    }

    /// Port whose signed integrated values occupy this manifest row.
    #[must_use]
    pub fn local_port_id(&self) -> &StableId {
        &self.local_port_id
    }

    /// Opposite port, reservoir, source boundary, or environment identity.
    #[must_use]
    pub fn counterparty_id(&self) -> &StableId {
        &self.counterparty_id
    }

    /// Internal, external, or included-source role.
    #[must_use]
    pub fn role(&self) -> WindowExchangeRole {
        self.role
    }

    /// Physical port vocabulary.
    #[must_use]
    pub fn port_kind(&self) -> WindowPortKind {
        self.port_kind
    }

    /// Complete thermodynamic context retained for a stream row.
    #[must_use]
    pub fn stream_binding(&self) -> Option<&StreamChartBinding> {
        self.stream_binding.as_deref()
    }

    /// Port schema timestamp retained for clock/window admission.
    #[must_use]
    pub fn timestamp(&self) -> &PortTimestamp {
        &self.timestamp
    }

    /// Explicit accounting boundary for external/source entries.
    #[must_use]
    pub fn boundary(&self) -> Option<&AccountingBoundary> {
        self.boundary.as_ref()
    }

    /// Evidence binding this contribution to its local port; closure evidence
    /// separately binds the complete row including counterparty and topology.
    #[must_use]
    pub fn binding_evidence(&self) -> &WindowEvidenceRef {
        &self.binding_evidence
    }
}

/// Complete declared boundary/source set and acceptance policy for one window.
#[derive(Debug, Clone, PartialEq)]
pub struct AccountingWindowSpec {
    version: u16,
    interval: AccountingWindowInterval,
    energy_reference: StableId,
    element_schema: WindowElementSchema,
    charge_schema: WindowChargeSchema,
    entropy_convention: StableId,
    manifest: Vec<WindowManifestEntry>,
    tolerances: WindowAuditTolerances,
    environment: Option<ExergyEnvironment>,
    closure_evidence: WindowEvidenceRef,
}

impl AccountingWindowSpec {
    /// Admit a complete, canonically ordered expected-contribution manifest.
    ///
    /// # Errors
    ///
    /// Refuses duplicate contribution or local-port IDs, out-of-window
    /// timestamps, malformed internal pairs, or foreign evidence. An empty
    /// manifest is a valid isolated window only under the explicit closure
    /// evidence supplied here.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        interval: AccountingWindowInterval,
        energy_reference: StableId,
        element_schema: WindowElementSchema,
        charge_schema: WindowChargeSchema,
        entropy_convention: StableId,
        manifest: impl IntoIterator<Item = WindowManifestEntry>,
        tolerances: WindowAuditTolerances,
        environment: Option<ExergyEnvironment>,
        closure_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        require_window_interval(&interval, closure_evidence.interval())?;
        if let Some(projection_evidence) = element_schema.projection_evidence() {
            require_window_interval(&interval, projection_evidence.interval())?;
            require_window_evidence_role(
                projection_evidence,
                &WindowEvidenceRole::ElementProjection {
                    elements: element_schema.elements().to_vec(),
                },
            )?;
        }
        require_window_interval(&interval, charge_schema.projection_evidence().interval())?;
        require_window_evidence_role(
            charge_schema.projection_evidence(),
            &charge_schema.evidence_role(),
        )?;
        let mut manifest: Vec<_> = manifest.into_iter().collect();
        manifest.sort_by(|left, right| left.contribution_id.cmp(&right.contribution_id));
        require_window_evidence_role(
            &closure_evidence,
            &WindowEvidenceRole::manifest_closure(&manifest),
        )?;
        for pair in manifest.windows(2) {
            if pair[0].contribution_id == pair[1].contribution_id {
                return Err(CoupleError::DuplicateAccountingWindowContribution {
                    contribution: pair[0].contribution_id.as_str().to_string(),
                });
            }
        }
        let mut local_ports: Vec<_> = manifest.iter().map(|entry| &entry.local_port_id).collect();
        local_ports.sort();
        for pair in local_ports.windows(2) {
            if pair[0] == pair[1] {
                return Err(CoupleError::DuplicateAccountingWindowPort {
                    port: pair[0].as_str().to_string(),
                });
            }
        }
        for entry in &manifest {
            if !interval.contains(&entry.timestamp) {
                return Err(CoupleError::AccountingWindowContextMismatch {
                    field: "manifest_timestamp",
                });
            }
            require_window_interval(&interval, entry.binding_evidence.interval())?;
        }
        let requires_elements = manifest.iter().any(|entry| {
            matches!(
                entry.port_kind,
                WindowPortKind::Stream
                    | WindowPortKind::Power(PortKind::ChemicalPotentialAmountFlow)
            )
        });
        if requires_elements && matches!(&element_schema, WindowElementSchema::NotApplicable { .. })
        {
            return Err(CoupleError::AccountingWindowElementSchemaRequired);
        }
        validate_window_manifest_groups(&manifest)?;
        validate_window_stream_contexts(&manifest)?;
        Ok(Self {
            version: ACCOUNTING_WINDOW_VERSION,
            interval,
            energy_reference,
            element_schema,
            charge_schema,
            entropy_convention,
            manifest,
            tolerances,
            environment,
            closure_evidence,
        })
    }

    /// Accounting-window schema version.
    #[must_use]
    pub fn version(&self) -> u16 {
        self.version
    }

    /// Exact logical interval.
    #[must_use]
    pub fn interval(&self) -> &AccountingWindowInterval {
        &self.interval
    }

    /// Common state/stream energy-reference contract for this audit.
    #[must_use]
    pub fn energy_reference(&self) -> &StableId {
        &self.energy_reference
    }

    /// Element projection/applicability schema.
    #[must_use]
    pub fn element_schema(&self) -> &WindowElementSchema {
        &self.element_schema
    }

    /// Direct-coulomb/species-projection/neutrality charge policy.
    #[must_use]
    pub fn charge_schema(&self) -> &WindowChargeSchema {
        &self.charge_schema
    }

    /// Nonequilibrium convention shared by every entropy ledger row.
    #[must_use]
    pub fn entropy_convention(&self) -> &StableId {
        &self.entropy_convention
    }

    /// Canonically ordered expected contribution rows.
    #[must_use]
    pub fn manifest(&self) -> &[WindowManifestEntry] {
        &self.manifest
    }

    /// Typed physical acceptance tolerances.
    #[must_use]
    pub fn tolerances(&self) -> WindowAuditTolerances {
        self.tolerances
    }

    /// Optional environment for `T0 * S_gen`.
    #[must_use]
    pub fn environment(&self) -> Option<&ExergyEnvironment> {
        self.environment.as_ref()
    }

    /// Evidence claiming that the declared manifest is complete.
    #[must_use]
    pub fn closure_evidence(&self) -> &WindowEvidenceRef {
        &self.closure_evidence
    }

    fn manifest_entry(&self, contribution_id: &StableId) -> Option<&WindowManifestEntry> {
        self.manifest
            .binary_search_by(|entry| entry.contribution_id.cmp(contribution_id))
            .ok()
            .map(|index| &self.manifest[index])
    }
}

/// Initial or final endpoint of an accounting window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowEndpoint {
    /// Stored state at `interval.start`.
    Initial,
    /// Stored state at `interval.end`.
    Final,
}

/// Evidence-bound stored extensive values at one exact window endpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowInventorySnapshot {
    interval: AccountingWindowInterval,
    energy_reference: StableId,
    element_schema: WindowElementSchema,
    charge_schema: WindowChargeSchema,
    entropy_convention: StableId,
    endpoint: WindowEndpoint,
    balance: WindowBalance,
    state_evidence: WindowEvidenceRef,
}

impl WindowInventorySnapshot {
    /// Admit one endpoint inventory against a window's exact element schema.
    ///
    /// # Errors
    ///
    /// Refuses foreign evidence, an element-axis mismatch, or negative stored
    /// mass/element amounts. Energy and entropy reference values may be signed.
    pub fn try_new(
        spec: &AccountingWindowSpec,
        endpoint: WindowEndpoint,
        balance: WindowBalance,
        state_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        require_window_interval(&spec.interval, state_evidence.interval())?;
        require_window_evidence_role(
            &state_evidence,
            &WindowEvidenceRole::Inventory {
                endpoint,
                energy_reference: spec.energy_reference.clone(),
            },
        )?;
        require_window_element_axis(&spec.element_schema, &balance, "inventory_element_axis")?;
        require_window_charge_admission(&spec.charge_schema, &balance, "inventory_charge")?;
        if balance.mass.value() < 0.0 {
            return Err(CoupleError::NegativeAccountingWindowInventory { field: "mass" });
        }
        if balance
            .element_amounts
            .iter()
            .any(|entry| entry.amount.value() < 0.0)
        {
            return Err(CoupleError::NegativeAccountingWindowInventory {
                field: "element_amount",
            });
        }
        Ok(Self {
            interval: spec.interval.clone(),
            energy_reference: spec.energy_reference.clone(),
            element_schema: spec.element_schema.clone(),
            charge_schema: spec.charge_schema.clone(),
            entropy_convention: spec.entropy_convention.clone(),
            endpoint,
            balance,
            state_evidence,
        })
    }

    /// Initial or final endpoint.
    #[must_use]
    pub fn endpoint(&self) -> WindowEndpoint {
        self.endpoint
    }

    /// Stored extensive values.
    #[must_use]
    pub fn balance(&self) -> &WindowBalance {
        &self.balance
    }

    /// Evidence covering this state snapshot.
    #[must_use]
    pub fn state_evidence(&self) -> &WindowEvidenceRef {
        &self.state_evidence
    }
}

/// Pre-integrated transfer values, normalized positive into the audited system.
#[derive(Debug, Clone, PartialEq)]
pub struct IntegratedWindowTransfer {
    energy_into: Energy,
    mass_into: Mass,
    element_amounts_into: Vec<WindowElementAmount>,
    charge_into: ElectricCharge,
    entropy_into: BoundaryEntropyBreakdown,
    entropy_convention: StableId,
}

impl IntegratedWindowTransfer {
    /// Admit one finite, explicit integrated transfer decomposition.
    ///
    /// The caller must integrate in physical time and normalize an
    /// owner-outward stream by negating its integrated rates. The retained
    /// integration evidence carries that obligation.
    ///
    /// # Errors
    ///
    /// Refuses duplicate element IDs or a non-finite fixed-order total.
    pub fn try_new(
        energy_into: Energy,
        mass_into: Mass,
        element_amounts_into: impl IntoIterator<Item = WindowElementAmount>,
        charge_into: ElectricCharge,
        entropy_into: BoundaryEntropyBreakdown,
        entropy_convention: StableId,
    ) -> Result<Self, CoupleError> {
        let balance = WindowBalance::try_new(
            energy_into,
            mass_into,
            element_amounts_into,
            charge_into,
            entropy_into.total(),
        )?;
        Ok(Self {
            energy_into: balance.energy,
            mass_into: balance.mass,
            element_amounts_into: balance.element_amounts,
            charge_into: balance.charge,
            entropy_into,
            entropy_convention,
        })
    }

    /// Integrated energy transfer, positive into system.
    #[must_use]
    pub fn energy_into(&self) -> Energy {
        self.energy_into
    }

    /// Integrated mass transfer, positive into system.
    #[must_use]
    pub fn mass_into(&self) -> Mass {
        self.mass_into
    }

    /// Integrated elemental transfers, positive into system.
    #[must_use]
    pub fn element_amounts_into(&self) -> &[WindowElementAmount] {
        &self.element_amounts_into
    }

    /// Integrated charge transfer, positive into system.
    #[must_use]
    pub fn charge_into(&self) -> ElectricCharge {
        self.charge_into
    }

    /// Explicit advected/diffusive/heat/radiation entropy transfer.
    #[must_use]
    pub fn entropy_into(&self) -> &BoundaryEntropyBreakdown {
        &self.entropy_into
    }

    /// Nonequilibrium convention governing the four entropy slots.
    #[must_use]
    pub fn entropy_convention(&self) -> &StableId {
        &self.entropy_convention
    }

    fn as_balance(&self) -> WindowBalance {
        WindowBalance {
            energy: self.energy_into,
            mass: self.mass_into,
            element_amounts: self.element_amounts_into.clone(),
            charge: self.charge_into,
            entropy: self.entropy_into.total(),
        }
    }
}

/// One complete manifest row filled with pre-integrated physical values.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowBoundaryContribution {
    manifest_entry: WindowManifestEntry,
    energy_reference: StableId,
    element_schema: WindowElementSchema,
    charge_schema: WindowChargeSchema,
    transfer: IntegratedWindowTransfer,
    integration_evidence: WindowEvidenceRef,
}

impl WindowBoundaryContribution {
    /// Fill one expected manifest row.
    ///
    /// # Errors
    ///
    /// Refuses an undeclared ID, element-axis mismatch, or evidence bound to a
    /// different interval.
    pub fn try_new(
        spec: &AccountingWindowSpec,
        contribution_id: &StableId,
        transfer: IntegratedWindowTransfer,
        integration_evidence: WindowEvidenceRef,
    ) -> Result<Self, CoupleError> {
        let manifest_entry = spec.manifest_entry(contribution_id).ok_or_else(|| {
            CoupleError::AccountingWindowContributionSetMismatch {
                missing: Vec::new(),
                unexpected: vec![contribution_id.as_str().to_string()],
            }
        })?;
        require_window_interval(&spec.interval, integration_evidence.interval())?;
        require_window_evidence_role(
            &integration_evidence,
            &WindowEvidenceRole::IntegratedTransfer {
                contribution_id: manifest_entry.contribution_id.clone(),
                local_port_id: manifest_entry.local_port_id.clone(),
                energy_reference: spec.energy_reference.clone(),
            },
        )?;
        if let Some(temperature_evidence) = transfer.entropy_into.temperature_reference.evidence() {
            require_window_interval(&spec.interval, temperature_evidence.interval())?;
            require_window_evidence_role(
                temperature_evidence,
                &WindowEvidenceRole::BoundaryTemperature {
                    contribution_id: manifest_entry.contribution_id.clone(),
                },
            )?;
        }
        require_window_element_axis(
            &spec.element_schema,
            &transfer.as_balance(),
            "transfer_element_axis",
        )?;
        require_window_charge_admission(
            &spec.charge_schema,
            &transfer.as_balance(),
            "transfer_charge",
        )?;
        validate_window_transfer_shape(manifest_entry, &transfer)?;
        Ok(Self {
            manifest_entry: manifest_entry.clone(),
            energy_reference: spec.energy_reference.clone(),
            element_schema: spec.element_schema.clone(),
            charge_schema: spec.charge_schema.clone(),
            transfer,
            integration_evidence,
        })
    }

    /// Complete expected manifest row.
    #[must_use]
    pub fn manifest_entry(&self) -> &WindowManifestEntry {
        &self.manifest_entry
    }

    /// Signed values normalized positive into system.
    #[must_use]
    pub fn transfer(&self) -> &IntegratedWindowTransfer {
        &self.transfer
    }

    /// Physical-time integration and sign-normalization evidence.
    #[must_use]
    pub fn integration_evidence(&self) -> &WindowEvidenceRef {
        &self.integration_evidence
    }
}

/// Conservation/inequality axis named by a physical window violation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowBalanceLaw {
    /// First-law energy closure.
    FirstLaw,
    /// Total mass closure.
    Mass,
    /// One elemental-amount closure.
    Element(ElementId),
    /// Electric-charge closure.
    Charge,
    /// Nonnegative entropy generation.
    EntropyGeneration,
    /// Pair-local energy closure for an internal transfer.
    InternalFirstLaw,
    /// Pair-local mass closure for an internal transfer.
    InternalMass,
    /// Pair-local elemental closure for an internal transfer.
    InternalElement(ElementId),
    /// Pair-local charge closure for an internal transfer.
    InternalCharge,
    /// Pair-local entropy-transfer closure for an internal transfer.
    InternalEntropyTransfer,
}

/// Structured physical failure retaining all nonzero involved ledger rows.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowAuditViolation {
    law: WindowBalanceLaw,
    residual_si: f64,
    tolerance_si: f64,
    contribution_ids: Vec<StableId>,
    port_ids: Vec<StableId>,
    exchange_ids: Vec<StableId>,
}

impl WindowAuditViolation {
    /// Failed balance or one-sided inequality.
    #[must_use]
    pub fn law(&self) -> &WindowBalanceLaw {
        &self.law
    }

    /// Signed residual in the coherent SI unit implied by [`Self::law`].
    #[must_use]
    pub fn residual_si(&self) -> f64 {
        self.residual_si
    }

    /// Typed policy tolerance expressed in the same coherent SI unit.
    #[must_use]
    pub fn tolerance_si(&self) -> f64 {
        self.tolerance_si
    }

    /// Canonically ordered nonzero contribution rows involved in this axis.
    #[must_use]
    pub fn contribution_ids(&self) -> &[StableId] {
        &self.contribution_ids
    }

    /// Canonically ordered local ports involved in this axis.
    #[must_use]
    pub fn port_ids(&self) -> &[StableId] {
        &self.port_ids
    }

    /// Canonically ordered pair/relation/exchange groups involved in this axis.
    #[must_use]
    pub fn exchange_ids(&self) -> &[StableId] {
        &self.exchange_ids
    }
}

/// Per-pair/relation energy and entropy ledger for localization.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowExchangeLedger {
    exchange_id: StableId,
    role: WindowExchangeRole,
    contribution_ids: Vec<StableId>,
    port_ids: Vec<StableId>,
    energy_into: Energy,
    mass_into: Mass,
    element_amounts_into: Vec<WindowElementAmount>,
    charge_into: ElectricCharge,
    entropy_into: Entropy,
}

impl WindowExchangeLedger {
    /// Pair/relation/exchange identifier.
    #[must_use]
    pub fn exchange_id(&self) -> &StableId {
        &self.exchange_id
    }

    /// Internal transfer, external exchange, or included source.
    #[must_use]
    pub fn role(&self) -> WindowExchangeRole {
        self.role
    }

    /// Canonically ordered contribution rows in this group.
    #[must_use]
    pub fn contribution_ids(&self) -> &[StableId] {
        &self.contribution_ids
    }

    /// Canonically ordered local ports in this group.
    #[must_use]
    pub fn port_ids(&self) -> &[StableId] {
        &self.port_ids
    }

    /// Net integrated energy into the audited system from this group.
    #[must_use]
    pub fn energy_into(&self) -> Energy {
        self.energy_into
    }

    /// Net integrated mass into the audited system from this group.
    #[must_use]
    pub fn mass_into(&self) -> Mass {
        self.mass_into
    }

    /// Net integrated elemental amounts into the system from this group.
    #[must_use]
    pub fn element_amounts_into(&self) -> &[WindowElementAmount] {
        &self.element_amounts_into
    }

    /// Net integrated charge into the audited system from this group.
    #[must_use]
    pub fn charge_into(&self) -> ElectricCharge {
        self.charge_into
    }

    /// Net integrated entropy into the audited system from this group.
    #[must_use]
    pub fn entropy_into(&self) -> Entropy {
        self.entropy_into
    }
}

/// Typed residuals for all PR-4 closed-window laws.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowResiduals {
    first_law: Energy,
    mass: Mass,
    elements: Vec<WindowElementAmount>,
    charge: ElectricCharge,
    entropy_generation: Entropy,
}

impl WindowResiduals {
    /// `delta stored energy - total energy into system`.
    #[must_use]
    pub fn first_law(&self) -> Energy {
        self.first_law
    }

    /// `delta stored mass - total mass into system`.
    #[must_use]
    pub fn mass(&self) -> Mass {
        self.mass
    }

    /// Per-element `delta stored - total into system` residuals.
    #[must_use]
    pub fn elements(&self) -> &[WindowElementAmount] {
        &self.elements
    }

    /// `delta stored charge - total charge into system`.
    #[must_use]
    pub fn charge(&self) -> ElectricCharge {
        self.charge
    }

    /// `delta stored entropy - total boundary/source entropy into system`.
    #[must_use]
    pub fn entropy_generation(&self) -> Entropy {
        self.entropy_generation
    }
}

/// Optional Gouy-Stodola diagnostic, not a complete availability balance.
#[derive(Debug, Clone, PartialEq)]
pub struct ExergyLedger {
    environment: ExergyEnvironment,
    entropy_generation: Entropy,
    destruction: Energy,
}

impl ExergyLedger {
    /// Explicit reference environment.
    #[must_use]
    pub fn environment(&self) -> &ExergyEnvironment {
        &self.environment
    }

    /// Audited entropy generation used in the product.
    #[must_use]
    pub fn entropy_generation(&self) -> Entropy {
        self.entropy_generation
    }

    /// `T0 * S_gen` in joules; may be slightly negative within tolerance.
    #[must_use]
    pub fn destruction(&self) -> Energy {
        self.destruction
    }
}

/// Retained green or red result of one structurally complete physical audit.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowAuditReport {
    spec: AccountingWindowSpec,
    initial: WindowInventorySnapshot,
    final_snapshot: WindowInventorySnapshot,
    contributions: Vec<WindowBoundaryContribution>,
    exchanges: Vec<WindowExchangeLedger>,
    residuals: WindowResiduals,
    violations: Vec<WindowAuditViolation>,
    exergy: Option<ExergyLedger>,
}

impl WindowAuditReport {
    /// Audit a complete set of pre-integrated window contributions.
    ///
    /// A physically red window returns `Ok(Self)` with retained violations.
    /// `Err` is reserved for malformed, incomplete, foreign, duplicate, or
    /// non-finite evidence/input.
    pub fn audit(
        spec: AccountingWindowSpec,
        initial: WindowInventorySnapshot,
        final_snapshot: WindowInventorySnapshot,
        contributions: impl IntoIterator<Item = WindowBoundaryContribution>,
    ) -> Result<Self, CoupleError> {
        require_window_snapshot(&spec, &initial, WindowEndpoint::Initial)?;
        require_window_snapshot(&spec, &final_snapshot, WindowEndpoint::Final)?;
        let mut contributions: Vec<_> = contributions.into_iter().collect();
        contributions.sort_by(|left, right| {
            left.manifest_entry
                .contribution_id
                .cmp(&right.manifest_entry.contribution_id)
        });
        for pair in contributions.windows(2) {
            if pair[0].manifest_entry.contribution_id == pair[1].manifest_entry.contribution_id {
                return Err(CoupleError::DuplicateAccountingWindowContribution {
                    contribution: pair[0].manifest_entry.contribution_id.as_str().to_string(),
                });
            }
        }
        require_complete_window_contributions(&spec, &contributions)?;
        let residuals = compute_window_residuals(&spec, &initial, &final_snapshot, &contributions)?;
        let exchanges = build_window_exchange_ledger(&contributions)?;
        let violations = build_window_violations(&spec, &residuals, &contributions, &exchanges);
        let exergy = spec
            .environment
            .clone()
            .map(|environment| {
                let destruction =
                    environment.temperature.value() * residuals.entropy_generation.value();
                ensure_finite_window_value("exergy_destruction", None, destruction)?;
                Ok(ExergyLedger {
                    environment,
                    entropy_generation: residuals.entropy_generation,
                    destruction: Energy::new(destruction),
                })
            })
            .transpose()?;

        Ok(Self {
            spec,
            initial,
            final_snapshot,
            contributions,
            exchanges,
            residuals,
            violations,
            exergy,
        })
    }

    /// Window specification and completeness evidence.
    #[must_use]
    pub fn spec(&self) -> &AccountingWindowSpec {
        &self.spec
    }

    /// Initial retained inventory.
    #[must_use]
    pub fn initial(&self) -> &WindowInventorySnapshot {
        &self.initial
    }

    /// Final retained inventory.
    #[must_use]
    pub fn final_snapshot(&self) -> &WindowInventorySnapshot {
        &self.final_snapshot
    }

    /// Complete canonical per-port ledger with raw signed contributions.
    #[must_use]
    pub fn contributions(&self) -> &[WindowBoundaryContribution] {
        &self.contributions
    }

    /// Canonical pair/relation/exchange localization ledger.
    #[must_use]
    pub fn exchanges(&self) -> &[WindowExchangeLedger] {
        &self.exchanges
    }

    /// All computed physical residuals, including raw entropy generation.
    #[must_use]
    pub fn residuals(&self) -> &WindowResiduals {
        &self.residuals
    }

    /// All physical failures; empty means the declared balances are green.
    #[must_use]
    pub fn violations(&self) -> &[WindowAuditViolation] {
        &self.violations
    }

    /// Whether every equality and the one-sided entropy inequality passes.
    #[must_use]
    pub fn is_green(&self) -> bool {
        self.violations.is_empty()
    }

    /// Optional narrow `T0 * S_gen` diagnostic.
    #[must_use]
    pub fn exergy(&self) -> Option<&ExergyLedger> {
        self.exergy.as_ref()
    }
}

fn ensure_finite_window_value(
    field: &'static str,
    context: Option<&str>,
    value: f64,
) -> Result<(), CoupleError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(CoupleError::NonFiniteAccountingWindowValue {
            field,
            context: context.map(str::to_string),
        })
    }
}

fn ensure_nonnegative_window_value(field: &'static str, value: f64) -> Result<(), CoupleError> {
    ensure_finite_window_value(field, None, value)?;
    if value < 0.0 {
        Err(CoupleError::NegativeAccountingWindowTolerance { field })
    } else {
        Ok(())
    }
}

fn require_window_evidence_role(
    evidence: &WindowEvidenceRef,
    expected: &WindowEvidenceRole,
) -> Result<(), CoupleError> {
    if &evidence.role == expected {
        Ok(())
    } else if evidence.role.diagnostic_label() == expected.diagnostic_label() {
        Err(CoupleError::AccountingWindowContextMismatch {
            field: "evidence_subject",
        })
    } else {
        Err(CoupleError::AccountingWindowEvidenceRoleMismatch {
            expected: expected.diagnostic_label(),
            actual: evidence.role.diagnostic_label(),
        })
    }
}

fn validate_window_manifest_groups(manifest: &[WindowManifestEntry]) -> Result<(), CoupleError> {
    let mut exchange_ids: Vec<_> = manifest
        .iter()
        .map(|entry| entry.exchange_id.clone())
        .collect();
    exchange_ids.sort();
    exchange_ids.dedup();
    for exchange_id in exchange_ids {
        let group: Vec<_> = manifest
            .iter()
            .filter(|entry| entry.exchange_id == exchange_id)
            .collect();
        let role = group[0].role;
        if group.iter().any(|entry| entry.role != role) {
            return Err(CoupleError::MalformedAccountingWindowExchange {
                exchange: exchange_id.as_str().to_string(),
                field: "mixed_roles",
            });
        }
        if role != WindowExchangeRole::InternalTransfer {
            continue;
        }
        if group.len() != 2 {
            return Err(CoupleError::MalformedAccountingWindowExchange {
                exchange: exchange_id.as_str().to_string(),
                field: "internal_cardinality",
            });
        }
        let [first, second] = [group[0], group[1]];
        if first.local_port_id != second.counterparty_id
            || first.counterparty_id != second.local_port_id
        {
            return Err(CoupleError::MalformedAccountingWindowExchange {
                exchange: exchange_id.as_str().to_string(),
                field: "internal_reciprocity",
            });
        }
        if first.port_kind != second.port_kind {
            return Err(CoupleError::MalformedAccountingWindowExchange {
                exchange: exchange_id.as_str().to_string(),
                field: "internal_port_kind",
            });
        }
        if first.timestamp != second.timestamp {
            return Err(CoupleError::MalformedAccountingWindowExchange {
                exchange: exchange_id.as_str().to_string(),
                field: "internal_timestamp",
            });
        }
        if first.boundary.is_some() || second.boundary.is_some() {
            return Err(CoupleError::MalformedAccountingWindowExchange {
                exchange: exchange_id.as_str().to_string(),
                field: "internal_boundary",
            });
        }
    }
    Ok(())
}

fn validate_window_stream_contexts(manifest: &[WindowManifestEntry]) -> Result<(), CoupleError> {
    let mut bindings = manifest
        .iter()
        .filter_map(|entry| entry.stream_binding.as_deref());
    let Some(reference) = bindings.next() else {
        return Ok(());
    };
    for binding in bindings {
        let mismatch = if reference.state_schema != binding.state_schema {
            Some("stream_state_schema")
        } else if reference.constituent_basis != binding.constituent_basis {
            Some("stream_constituent_basis")
        } else if reference.constituent_axis != binding.constituent_axis {
            Some("stream_constituent_axis")
        } else if reference.chemical_reference_state != binding.chemical_reference_state {
            Some("stream_chemical_reference_state")
        } else if reference.coordinates != binding.coordinates {
            Some("stream_coordinates")
        } else if reference.timestamp.clock != binding.timestamp.clock {
            Some("stream_clock")
        } else if reference.gravity_datum != binding.gravity_datum {
            Some("stream_gravity_datum")
        } else if reference.stress_convention != binding.stress_convention {
            Some("stream_stress_convention")
        } else {
            None
        };
        if let Some(field) = mismatch {
            return Err(CoupleError::AccountingWindowContextMismatch { field });
        }
    }
    Ok(())
}

fn require_window_interval(
    expected: &AccountingWindowInterval,
    actual: &AccountingWindowInterval,
) -> Result<(), CoupleError> {
    if expected == actual {
        Ok(())
    } else {
        Err(CoupleError::AccountingWindowContextMismatch { field: "interval" })
    }
}

fn require_internal_stream_pair(
    first: &StreamChartBinding,
    second: &StreamChartBinding,
) -> Result<(), CoupleError> {
    if first.port_id == second.port_id {
        return Err(CoupleError::DuplicateIdentity {
            id: first.port_id.as_str().to_string(),
        });
    }
    let mismatch = if first.state_schema != second.state_schema {
        Some("state_schema")
    } else if first.constituent_basis != second.constituent_basis {
        Some("constituent_basis")
    } else if first.constituent_axis != second.constituent_axis {
        Some("constituent_axis")
    } else if first.chemical_reference_state != second.chemical_reference_state {
        Some("chemical_reference_state")
    } else if first.coordinates != second.coordinates {
        Some("coordinates")
    } else if first.timestamp != second.timestamp {
        Some("clock_timestamp")
    } else if first.gravity_datum != second.gravity_datum {
        Some("gravity_datum")
    } else if first.stress_convention != second.stress_convention {
        Some("stress_convention")
    } else {
        None
    };
    mismatch.map_or(Ok(()), |field| {
        Err(CoupleError::IncompatibleStreamPorts {
            a: first.port_id.as_str().to_string(),
            b: second.port_id.as_str().to_string(),
            field,
        })
    })
}

fn require_window_element_axis(
    schema: &WindowElementSchema,
    balance: &WindowBalance,
    context: &'static str,
) -> Result<(), CoupleError> {
    let expected = schema.elements();
    let actual: Vec<_> = balance
        .element_amounts
        .iter()
        .map(|entry| &entry.element)
        .collect();
    if expected.iter().eq(actual.iter().copied()) {
        Ok(())
    } else {
        Err(CoupleError::AccountingWindowElementAxisMismatch {
            context,
            expected: expected
                .iter()
                .map(|element| element.as_str().to_string())
                .collect(),
            actual: actual
                .iter()
                .map(|element| element.as_str().to_string())
                .collect(),
        })
    }
}

fn require_window_charge_admission(
    schema: &WindowChargeSchema,
    balance: &WindowBalance,
    context: &'static str,
) -> Result<(), CoupleError> {
    if matches!(schema, WindowChargeSchema::ProvenNeutral { .. })
        && balance.charge.value().abs() > 0.0
    {
        Err(CoupleError::AccountingWindowChargeSchemaMismatch { context })
    } else {
        Ok(())
    }
}

fn reject_window_transfer_field(
    entry: &WindowManifestEntry,
    field: &'static str,
    nonzero: bool,
) -> Result<(), CoupleError> {
    if nonzero {
        Err(CoupleError::AccountingWindowPortTransferMismatch {
            port: entry.local_port_id.as_str().to_string(),
            field,
        })
    } else {
        Ok(())
    }
}

fn validate_window_transfer_shape(
    entry: &WindowManifestEntry,
    transfer: &IntegratedWindowTransfer,
) -> Result<(), CoupleError> {
    let WindowPortKind::Power(kind) = entry.port_kind else {
        return Ok(());
    };
    reject_window_transfer_field(
        entry,
        "mass_on_power_port",
        transfer.mass_into.value().abs() > 0.0,
    )?;
    if kind != PortKind::ChemicalPotentialAmountFlow {
        reject_window_transfer_field(
            entry,
            "element_amount_on_nonchemical_power_port",
            transfer
                .element_amounts_into
                .iter()
                .any(|element| element.amount.value().abs() > 0.0),
        )?;
    }
    if !matches!(
        kind,
        PortKind::ElectricalVoltageCurrent | PortKind::ChemicalPotentialAmountFlow
    ) {
        reject_window_transfer_field(
            entry,
            "charge_on_noncharge_power_port",
            transfer.charge_into.value().abs() > 0.0,
        )?;
    }
    if kind != PortKind::ThermalTemperatureEntropy {
        return reject_window_transfer_field(
            entry,
            "entropy_on_nonthermal_power_port",
            [
                transfer.entropy_into.advected.value(),
                transfer.entropy_into.diffusive_chemical.value(),
                transfer.entropy_into.heat.value(),
                transfer.entropy_into.radiation.value(),
            ]
            .into_iter()
            .any(|value| value.abs() > 0.0),
        );
    }

    reject_window_transfer_field(
        entry,
        "advected_entropy_on_thermal_power_port",
        transfer.entropy_into.advected.value().abs() > 0.0,
    )?;
    reject_window_transfer_field(
        entry,
        "diffusive_entropy_on_thermal_power_port",
        transfer.entropy_into.diffusive_chemical.value().abs() > 0.0,
    )?;
    match &transfer.entropy_into.temperature_reference {
        BoundaryTemperatureReference::NotApplicable { .. } => reject_window_transfer_field(
            entry,
            "thermal_temperature_reference",
            transfer.energy_into.value().abs() > 0.0
                || transfer.entropy_into.heat.value().abs() > 0.0
                || transfer.entropy_into.radiation.value().abs() > 0.0,
        ),
        BoundaryTemperatureReference::Constant { temperature, .. } => {
            if transfer.entropy_into.radiation.value().abs() > 0.0 {
                return Err(CoupleError::AccountingWindowPortTransferMismatch {
                    port: entry.local_port_id.as_str().to_string(),
                    field: "constant_temperature_radiation_requires_profile",
                });
            }
            let expected_energy = temperature.value() * transfer.entropy_into.heat.value();
            ensure_finite_window_value(
                "constant_temperature_thermal_energy",
                Some(entry.contribution_id.as_str()),
                expected_energy,
            )?;
            if expected_energy.to_bits() == transfer.energy_into.value().to_bits() {
                Ok(())
            } else {
                Err(CoupleError::AccountingWindowPortTransferMismatch {
                    port: entry.local_port_id.as_str().to_string(),
                    field: "constant_temperature_energy_entropy_identity",
                })
            }
        }
        BoundaryTemperatureReference::Profile { .. } => Ok(()),
    }
}

fn require_window_snapshot(
    spec: &AccountingWindowSpec,
    snapshot: &WindowInventorySnapshot,
    endpoint: WindowEndpoint,
) -> Result<(), CoupleError> {
    require_window_interval(&spec.interval, &snapshot.interval)?;
    require_window_interval(&spec.interval, snapshot.state_evidence.interval())?;
    if snapshot.energy_reference != spec.energy_reference {
        return Err(CoupleError::AccountingWindowContextMismatch {
            field: "snapshot_energy_reference",
        });
    }
    if snapshot.element_schema != spec.element_schema {
        return Err(CoupleError::AccountingWindowContextMismatch {
            field: "snapshot_element_schema",
        });
    }
    if snapshot.charge_schema != spec.charge_schema {
        return Err(CoupleError::AccountingWindowContextMismatch {
            field: "snapshot_charge_schema",
        });
    }
    if snapshot.entropy_convention != spec.entropy_convention {
        return Err(CoupleError::AccountingWindowContextMismatch {
            field: "snapshot_entropy_convention",
        });
    }
    require_window_evidence_role(
        &snapshot.state_evidence,
        &WindowEvidenceRole::Inventory {
            endpoint,
            energy_reference: spec.energy_reference.clone(),
        },
    )?;
    if snapshot.endpoint != endpoint {
        return Err(CoupleError::AccountingWindowEndpointMismatch {
            expected: endpoint,
            actual: snapshot.endpoint,
        });
    }
    require_window_element_axis(
        &spec.element_schema,
        &snapshot.balance,
        "inventory_element_axis",
    )?;
    require_window_charge_admission(&spec.charge_schema, &snapshot.balance, "inventory_charge")
}

fn require_complete_window_contributions(
    spec: &AccountingWindowSpec,
    contributions: &[WindowBoundaryContribution],
) -> Result<(), CoupleError> {
    let expected: Vec<_> = spec
        .manifest
        .iter()
        .map(|entry| entry.contribution_id.clone())
        .collect();
    let actual: Vec<_> = contributions
        .iter()
        .map(|entry| entry.manifest_entry.contribution_id.clone())
        .collect();
    let missing: Vec<_> = expected
        .iter()
        .filter(|id| actual.binary_search(id).is_err())
        .map(|id| id.as_str().to_string())
        .collect();
    let unexpected: Vec<_> = actual
        .iter()
        .filter(|id| expected.binary_search(id).is_err())
        .map(|id| id.as_str().to_string())
        .collect();
    if !missing.is_empty() || !unexpected.is_empty() {
        return Err(CoupleError::AccountingWindowContributionSetMismatch {
            missing,
            unexpected,
        });
    }
    for (manifest, contribution) in spec.manifest.iter().zip(contributions) {
        if manifest != &contribution.manifest_entry {
            return Err(CoupleError::AccountingWindowContextMismatch {
                field: "manifest_entry",
            });
        }
        if contribution.energy_reference != spec.energy_reference {
            return Err(CoupleError::AccountingWindowContextMismatch {
                field: "contribution_energy_reference",
            });
        }
        if contribution.element_schema != spec.element_schema {
            return Err(CoupleError::AccountingWindowContextMismatch {
                field: "contribution_element_schema",
            });
        }
        if contribution.charge_schema != spec.charge_schema {
            return Err(CoupleError::AccountingWindowContextMismatch {
                field: "contribution_charge_schema",
            });
        }
        require_window_interval(&spec.interval, contribution.integration_evidence.interval())?;
        require_window_evidence_role(
            &contribution.integration_evidence,
            &WindowEvidenceRole::IntegratedTransfer {
                contribution_id: manifest.contribution_id.clone(),
                local_port_id: manifest.local_port_id.clone(),
                energy_reference: spec.energy_reference.clone(),
            },
        )?;
        if let Some(temperature_evidence) = contribution
            .transfer
            .entropy_into
            .temperature_reference
            .evidence()
        {
            require_window_interval(&spec.interval, temperature_evidence.interval())?;
            require_window_evidence_role(
                temperature_evidence,
                &WindowEvidenceRole::BoundaryTemperature {
                    contribution_id: manifest.contribution_id.clone(),
                },
            )?;
        }
        require_window_element_axis(
            &spec.element_schema,
            &contribution.transfer.as_balance(),
            "transfer_element_axis",
        )?;
        require_window_charge_admission(
            &spec.charge_schema,
            &contribution.transfer.as_balance(),
            "transfer_charge",
        )?;
        validate_window_transfer_shape(manifest, &contribution.transfer)?;
        if contribution.transfer.entropy_convention != spec.entropy_convention {
            return Err(CoupleError::AccountingWindowContextMismatch {
                field: "entropy_convention",
            });
        }
    }
    Ok(())
}

fn fixed_order_window_sum(
    field: &'static str,
    values: impl IntoIterator<Item = f64>,
) -> Result<f64, CoupleError> {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        ensure_finite_window_value(field, None, value)?;
        let next = sum + value;
        ensure_finite_window_value(field, None, next)?;
        let adjustment = if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        correction += adjustment;
        ensure_finite_window_value(field, None, correction)?;
        sum = next;
    }
    let total = sum + correction;
    ensure_finite_window_value(field, None, total)?;
    Ok(total)
}

fn window_residual(
    field: &'static str,
    initial: f64,
    final_value: f64,
    inputs: impl IntoIterator<Item = f64>,
) -> Result<f64, CoupleError> {
    let change = final_value - initial;
    ensure_finite_window_value(field, None, change)?;
    let total_input = fixed_order_window_sum(field, inputs)?;
    let residual = change - total_input;
    ensure_finite_window_value(field, None, residual)?;
    Ok(residual)
}

fn compute_window_residuals(
    spec: &AccountingWindowSpec,
    initial: &WindowInventorySnapshot,
    final_snapshot: &WindowInventorySnapshot,
    contributions: &[WindowBoundaryContribution],
) -> Result<WindowResiduals, CoupleError> {
    let first_law = window_residual(
        "first_law_residual",
        initial.balance.energy.value(),
        final_snapshot.balance.energy.value(),
        contributions
            .iter()
            .map(|entry| entry.transfer.energy_into.value()),
    )?;
    let mass = window_residual(
        "mass_residual",
        initial.balance.mass.value(),
        final_snapshot.balance.mass.value(),
        contributions
            .iter()
            .map(|entry| entry.transfer.mass_into.value()),
    )?;
    let charge = window_residual(
        "charge_residual",
        initial.balance.charge.value(),
        final_snapshot.balance.charge.value(),
        contributions
            .iter()
            .map(|entry| entry.transfer.charge_into.value()),
    )?;
    let entropy_generation = window_residual(
        "entropy_generation",
        initial.balance.entropy.value(),
        final_snapshot.balance.entropy.value(),
        contributions
            .iter()
            .map(|entry| entry.transfer.entropy_into.total().value()),
    )?;
    let mut elements = Vec::with_capacity(spec.element_schema.elements().len());
    for (index, element) in spec.element_schema.elements().iter().enumerate() {
        let residual = window_residual(
            "element_residual",
            initial.balance.element_amounts[index].amount.value(),
            final_snapshot.balance.element_amounts[index].amount.value(),
            contributions
                .iter()
                .map(|entry| entry.transfer.element_amounts_into[index].amount.value()),
        )?;
        elements.push(WindowElementAmount {
            element: element.clone(),
            amount: Amount::new(residual),
        });
    }
    Ok(WindowResiduals {
        first_law: Energy::new(first_law),
        mass: Mass::new(mass),
        elements,
        charge: ElectricCharge::new(charge),
        entropy_generation: Entropy::new(entropy_generation),
    })
}

fn build_window_exchange_ledger(
    contributions: &[WindowBoundaryContribution],
) -> Result<Vec<WindowExchangeLedger>, CoupleError> {
    let mut exchange_ids: Vec<_> = contributions
        .iter()
        .map(|entry| entry.manifest_entry.exchange_id.clone())
        .collect();
    exchange_ids.sort();
    exchange_ids.dedup();
    let mut ledgers = Vec::with_capacity(exchange_ids.len());
    for exchange_id in exchange_ids {
        let group: Vec<_> = contributions
            .iter()
            .filter(|entry| entry.manifest_entry.exchange_id == exchange_id)
            .collect();
        let energy_into = fixed_order_window_sum(
            "exchange_energy_into",
            group.iter().map(|entry| entry.transfer.energy_into.value()),
        )?;
        let mass_into = fixed_order_window_sum(
            "exchange_mass_into",
            group.iter().map(|entry| entry.transfer.mass_into.value()),
        )?;
        let charge_into = fixed_order_window_sum(
            "exchange_charge_into",
            group.iter().map(|entry| entry.transfer.charge_into.value()),
        )?;
        let entropy_into = fixed_order_window_sum(
            "exchange_entropy_into",
            group
                .iter()
                .map(|entry| entry.transfer.entropy_into.total().value()),
        )?;
        let mut element_amounts_into = Vec::new();
        if let Some(first) = group.first() {
            element_amounts_into.reserve(first.transfer.element_amounts_into.len());
            for (index, element) in first.transfer.element_amounts_into.iter().enumerate() {
                let amount = fixed_order_window_sum(
                    "exchange_element_amount_into",
                    group
                        .iter()
                        .map(|entry| entry.transfer.element_amounts_into[index].amount.value()),
                )?;
                element_amounts_into.push(WindowElementAmount {
                    element: element.element.clone(),
                    amount: Amount::new(amount),
                });
            }
        }
        let contribution_ids = group
            .iter()
            .map(|entry| entry.manifest_entry.contribution_id.clone())
            .collect();
        let mut port_ids: Vec<_> = group
            .iter()
            .map(|entry| entry.manifest_entry.local_port_id.clone())
            .collect();
        port_ids.sort();
        port_ids.dedup();
        ledgers.push(WindowExchangeLedger {
            exchange_id,
            role: group[0].manifest_entry.role,
            contribution_ids,
            port_ids,
            energy_into: Energy::new(energy_into),
            mass_into: Mass::new(mass_into),
            element_amounts_into,
            charge_into: ElectricCharge::new(charge_into),
            entropy_into: Entropy::new(entropy_into),
        });
    }
    Ok(ledgers)
}

fn transfer_value_for_law(transfer: &IntegratedWindowTransfer, law: &WindowBalanceLaw) -> f64 {
    match law {
        WindowBalanceLaw::FirstLaw | WindowBalanceLaw::InternalFirstLaw => {
            transfer.energy_into.value()
        }
        WindowBalanceLaw::Mass | WindowBalanceLaw::InternalMass => transfer.mass_into.value(),
        WindowBalanceLaw::Element(element) | WindowBalanceLaw::InternalElement(element) => transfer
            .element_amounts_into
            .iter()
            .find(|entry| &entry.element == element)
            .map_or(0.0, |entry| entry.amount.value()),
        WindowBalanceLaw::Charge | WindowBalanceLaw::InternalCharge => transfer.charge_into.value(),
        WindowBalanceLaw::EntropyGeneration | WindowBalanceLaw::InternalEntropyTransfer => {
            transfer.entropy_into.total().value()
        }
    }
}

fn window_violation(
    law: WindowBalanceLaw,
    residual_si: f64,
    tolerance_si: f64,
    contributions: &[WindowBoundaryContribution],
) -> WindowAuditViolation {
    let involved: Vec<_> = contributions
        .iter()
        .filter(|entry| transfer_value_for_law(&entry.transfer, &law).abs() > 0.0)
        .collect();
    let contribution_ids = involved
        .iter()
        .map(|entry| entry.manifest_entry.contribution_id.clone())
        .collect();
    let mut port_ids: Vec<_> = involved
        .iter()
        .map(|entry| entry.manifest_entry.local_port_id.clone())
        .collect();
    port_ids.sort();
    port_ids.dedup();
    let mut exchange_ids: Vec<_> = involved
        .iter()
        .map(|entry| entry.manifest_entry.exchange_id.clone())
        .collect();
    exchange_ids.sort();
    exchange_ids.dedup();
    WindowAuditViolation {
        law,
        residual_si,
        tolerance_si,
        contribution_ids,
        port_ids,
        exchange_ids,
    }
}

fn internal_exchange_violation(
    law: WindowBalanceLaw,
    residual_si: f64,
    tolerance_si: f64,
    exchange: &WindowExchangeLedger,
) -> WindowAuditViolation {
    WindowAuditViolation {
        law,
        residual_si,
        tolerance_si,
        contribution_ids: exchange.contribution_ids.clone(),
        port_ids: exchange.port_ids.clone(),
        exchange_ids: vec![exchange.exchange_id.clone()],
    }
}

fn build_window_violations(
    spec: &AccountingWindowSpec,
    residuals: &WindowResiduals,
    contributions: &[WindowBoundaryContribution],
    exchanges: &[WindowExchangeLedger],
) -> Vec<WindowAuditViolation> {
    let mut violations = Vec::new();
    for (law, residual, tolerance) in [
        (
            WindowBalanceLaw::FirstLaw,
            residuals.first_law.value(),
            spec.tolerances.energy.value(),
        ),
        (
            WindowBalanceLaw::Mass,
            residuals.mass.value(),
            spec.tolerances.mass.value(),
        ),
        (
            WindowBalanceLaw::Charge,
            residuals.charge.value(),
            spec.tolerances.charge.value(),
        ),
    ] {
        if residual.abs() > tolerance {
            violations.push(window_violation(law, residual, tolerance, contributions));
        }
    }
    for element in &residuals.elements {
        if element.amount.value().abs() > spec.tolerances.amount.value() {
            violations.push(window_violation(
                WindowBalanceLaw::Element(element.element.clone()),
                element.amount.value(),
                spec.tolerances.amount.value(),
                contributions,
            ));
        }
    }
    if residuals.entropy_generation.value() < -spec.tolerances.entropy.value() {
        violations.push(window_violation(
            WindowBalanceLaw::EntropyGeneration,
            residuals.entropy_generation.value(),
            spec.tolerances.entropy.value(),
            contributions,
        ));
    }
    for exchange in exchanges
        .iter()
        .filter(|exchange| exchange.role == WindowExchangeRole::InternalTransfer)
    {
        for (law, residual, tolerance) in [
            (
                WindowBalanceLaw::InternalFirstLaw,
                exchange.energy_into.value(),
                spec.tolerances.energy.value(),
            ),
            (
                WindowBalanceLaw::InternalMass,
                exchange.mass_into.value(),
                spec.tolerances.mass.value(),
            ),
            (
                WindowBalanceLaw::InternalCharge,
                exchange.charge_into.value(),
                spec.tolerances.charge.value(),
            ),
            (
                WindowBalanceLaw::InternalEntropyTransfer,
                exchange.entropy_into.value(),
                spec.tolerances.entropy.value(),
            ),
        ] {
            if residual.abs() > tolerance {
                violations.push(internal_exchange_violation(
                    law, residual, tolerance, exchange,
                ));
            }
        }
        for element in &exchange.element_amounts_into {
            if element.amount.value().abs() > spec.tolerances.amount.value() {
                violations.push(internal_exchange_violation(
                    WindowBalanceLaw::InternalElement(element.element.clone()),
                    element.amount.value(),
                    spec.tolerances.amount.value(),
                    exchange,
                ));
            }
        }
    }
    violations
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
    /// Applying a field measure to its declared density side overflowed.
    PortMeasureApplicationOverflow {
        /// Density side receiving the measure.
        side: PortVariable,
        /// Pointwise dimensions on that side.
        dimensions: Dims,
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
    /// A watt-dimensional pair used dimensions belonging to another port kind.
    PortKindDimensionMismatch {
        /// Declared physical kind.
        kind: PortKind,
        /// Effort or flow mismatch.
        side: PortVariable,
        /// Canonical generalized dimensions for the kind.
        expected: Dims,
        /// Declared generalized dimensions after field measure application.
        actual: Dims,
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
    /// A physical port kind omitted the conserved quantity implied by its
    /// effort/flow semantics.
    MissingPortKindConservationRole {
        /// Declared physical kind.
        kind: PortKind,
        /// Required but absent role.
        role: ConservationRole,
    },
    /// The raw scalar interconnection oracle received a v2 schema-only kind.
    LegacyInterconnectRequiresSeedKind {
        /// Rejected physical kind.
        kind: PortKind,
    },
    /// One scalar or vector component of a stream bundle was not finite.
    NonFiniteStreamValue {
        /// Offending semantic field.
        field: &'static str,
        /// Component index for vector-valued fields.
        index: Option<usize>,
    },
    /// A stream bundle omitted its species/element amount axis.
    EmptyStreamConstituents,
    /// A stream declared the same constituent more than once.
    DuplicateStreamConstituent {
        /// Repeated species or element identifier.
        id: String,
    },
    /// The stream rates did not exactly match the proof-bound constituent axis.
    StreamConstituentAxisMismatch {
        /// Canonical species/element axis bound by the chart proofs.
        expected: Vec<String>,
        /// Canonical species/element axis carried by the stream rates.
        actual: Vec<String>,
    },
    /// Species-potential accounting was requested on an axis containing an
    /// element coordinate.
    StreamChemicalEnergyRequiresSpeciesAxis {
        /// First non-species coordinate on the canonical axis.
        constituent: String,
    },
    /// Stream flow signs are admitted only in the owner-outward convention.
    StreamPortRequiresOutwardOrientation {
        /// Rejected orientation.
        actual: PortOrientation,
    },
    /// A chart, proof, or chemical partition was bound to a different stream
    /// context.
    StreamChartBindingMismatch {
        /// First mismatching context field.
        field: &'static str,
    },
    /// A proof reference named the wrong thermodynamic identity.
    WrongStreamIdentityProof {
        /// Identity required by the admission path.
        expected: StreamIdentity,
        /// Identity named by the supplied proof.
        actual: StreamIdentity,
    },
    /// A pressure-work or Euler/Legendre equality was not bit-exact.
    StreamIdentityMismatch {
        /// Identity being checked.
        identity: StreamIdentity,
        /// Localized equality within that identity.
        check: &'static str,
    },
    /// A pressure-work crosswalk used zero or negative density.
    NonPositiveStreamDensity,
    /// Chemical energy was owned by both the state potential and a separate
    /// species-potential contribution.
    DoubleCountedChemicalEnergy,
    /// The declared stream energy rate differed from the selected chart.
    StreamEnergyFlowMismatch {
        /// Bit pattern recomputed from the admitted chart.
        accounted_bits: u64,
        /// Bit pattern declared by the stream bundle.
        declared_bits: u64,
    },
    /// Initial and final logical timestamps used different clocks or evidence
    /// named a different complete window context.
    AccountingWindowContextMismatch {
        /// First mismatching context field.
        field: &'static str,
    },
    /// An accounting window had no positive logical extent.
    InvalidAccountingWindowInterval {
        /// Initial logical tick.
        start_tick: u64,
        /// Final logical tick.
        end_tick: u64,
    },
    /// A stored, integrated, derived, or tolerance value was not finite.
    NonFiniteAccountingWindowValue {
        /// Offending semantic field.
        field: &'static str,
        /// Optional element/contribution context.
        context: Option<String>,
    },
    /// An elemental inventory/transfer axis repeated an identifier.
    DuplicateAccountingWindowElement {
        /// Repeated element identifier.
        element: String,
    },
    /// An audited elemental schema omitted its axis.
    EmptyAccountingWindowElementAxis,
    /// A stream or chemical port was admitted without an elemental schema.
    AccountingWindowElementSchemaRequired,
    /// Retained evidence named the wrong semantic role or subject.
    AccountingWindowEvidenceRoleMismatch {
        /// Required evidence role.
        expected: &'static str,
        /// Supplied evidence role.
        actual: &'static str,
    },
    /// An exchange group had mixed roles, incomplete topology, or bad reciprocity.
    MalformedAccountingWindowExchange {
        /// Exchange identifier.
        exchange: String,
        /// First malformed structural field.
        field: &'static str,
    },
    /// A proven-neutral charge schema received a nonzero coulomb value.
    AccountingWindowChargeSchemaMismatch {
        /// Inventory or transfer context.
        context: &'static str,
    },
    /// A power-port row carried a balance axis inconsistent with its port kind.
    AccountingWindowPortTransferMismatch {
        /// Local port identifier.
        port: String,
        /// Incompatible balance field or exact thermal identity.
        field: &'static str,
    },
    /// A physical audit tolerance was negative.
    NegativeAccountingWindowTolerance {
        /// Offending tolerance field.
        field: &'static str,
    },
    /// The optional exergy environment used a non-finite or non-positive `T0`.
    NonPositiveExergyTemperature,
    /// A constant boundary temperature was non-finite or non-positive.
    NonPositiveBoundaryTemperature,
    /// Heat/radiation entropy was supplied without boundary-temperature evidence.
    AccountingWindowTemperatureReferenceRequired,
    /// Two stream sides disagreed on a localized transfer-context field.
    IncompatibleStreamPorts {
        /// First stream-port ID.
        a: String,
        /// Second stream-port ID.
        b: String,
        /// First incompatible context field.
        field: &'static str,
    },
    /// A window manifest assigned the wrong source/external treatment.
    AccountingWindowBoundaryTreatmentMismatch {
        /// Boundary identifier.
        boundary: String,
        /// Required treatment for the constructor.
        expected: BoundaryTreatment,
        /// Actual boundary treatment.
        actual: BoundaryTreatment,
    },
    /// An external accounting boundary lacked the owner-outward convention.
    AccountingWindowBoundaryRequiresOutward {
        /// Boundary identifier.
        boundary: String,
    },
    /// A manifest or filled ledger repeated one contribution ID.
    DuplicateAccountingWindowContribution {
        /// Repeated contribution identifier.
        contribution: String,
    },
    /// A manifest assigned more than one row to the same local port.
    DuplicateAccountingWindowPort {
        /// Repeated local port identifier.
        port: String,
    },
    /// A stored or transferred elemental axis differed from the declared axis.
    AccountingWindowElementAxisMismatch {
        /// Inventory or transfer context.
        context: &'static str,
        /// Declared canonical element axis.
        expected: Vec<String>,
        /// Actual canonical element axis.
        actual: Vec<String>,
    },
    /// A stored mass or elemental amount was negative.
    NegativeAccountingWindowInventory {
        /// Offending stored field.
        field: &'static str,
    },
    /// Filled contribution IDs did not exactly cover the declared manifest.
    AccountingWindowContributionSetMismatch {
        /// Declared rows without values.
        missing: Vec<String>,
        /// Supplied rows absent from the manifest.
        unexpected: Vec<String>,
    },
    /// Initial/final inventory snapshots were supplied in the wrong slots.
    AccountingWindowEndpointMismatch {
        /// Required endpoint.
        expected: WindowEndpoint,
        /// Actual endpoint.
        actual: WindowEndpoint,
    },
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
/// [`CoupleError::IncompatiblePorts`] if the ports are not power-conjugate, or
/// [`CoupleError::LegacyInterconnectRequiresSeedKind`] if a v2 schema-only kind is
/// passed to this migration oracle.
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
    if !kind_a.is_legacy_scalar_seed() {
        return Err(CoupleError::LegacyInterconnectRequiresSeedKind { kind: kind_a });
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
