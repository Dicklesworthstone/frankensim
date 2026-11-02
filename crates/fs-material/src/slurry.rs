//! Source-bound mineral/liquid flow accounting with additive constituent volumes.
//!
//! This deliberately small reference model uses `volume_flow = mass_flow / density`
//! for each solid constituent and one liquid carrier. There is no entrained gas,
//! volume contraction, reaction, dissolution, settling, pore-volume correction,
//! rheology, liberation, or thermal equilibration model. A mineral's reference
//! density is not promoted to a measured bulk-ore density or a plant calibration.
//! Source query conditions and uncertainty stay inspectable; uncertainty is not
//! silently propagated into a certified slurry-density interval.
//!
//! Streams contain nonnegative kg/s rates. Components are sorted by identity,
//! and sums use a fixed numeric order with compensation. Mixing preserves each
//! constituent; the base splitter is nonselective, while [`separation`] accepts
//! explicit component and carrier partitions. Density and composition are undefined
//! for a zero-flow stream. Arithmetic that overflows or erases a positive phase
//! through underflow is refused. All operations have explicit small size bounds.
//!
//! The separate [`rheology`] module evaluates source-qualified steady shear and
//! laminar pipe flow; its parameters and temperature require explicit declarations.

use fs_blake3::ContentHash;
use fs_matdb::{
    MatDbError, MaterialAnswer, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy,
};
use fs_qty::{Density, Dimensionless, Dims, MassFlowRate, QuantitySpec, VolumetricFlowRate};
use std::collections::BTreeMap;

/// Prescribed selective component/carrier partitions with conserved outlet streams.
#[path = "slurry/separation.rs"]
pub mod separation;

/// Source-qualified Herschel–Bulkley/Bingham rheology and laminar pipe hydraulics.
#[path = "slurry/rheology.rs"]
pub mod rheology;

/// Maximum distinct dry constituents admitted by a stream.
pub const MAX_SOLID_COMPONENTS: usize = 128;
/// Maximum inlet streams admitted by one mixing operation.
pub const MAX_MIXING_STREAMS: usize = 128;
const MAX_COMPONENT_ID_BYTES: usize = 256;
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);

/// Admission or arithmetic refusal; no silently repaired physical inputs.
#[derive(Debug, Clone, PartialEq)]
pub enum SlurryError {
    /// Original material query or receipt-verification refusal.
    Material(MatDbError),
    /// Invalid input, unsupported operation, or unrepresentable arithmetic.
    Invalid(&'static str),
    /// One component identity was used with different source/query evidence.
    SourceMismatch {
        /// Caller-defined identity that cannot be merged.
        component: String,
    },
    /// Addition alone cannot reach the requested solids mass fraction.
    WaterRemovalRequired {
        /// Required liquid removal, in kg/s; no removal has been performed.
        removal: MassFlowRate,
    },
    /// A zero-flow stream has no composition to dilute.
    UndefinedComposition,
}

impl core::fmt::Display for SlurryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "slurry material source: {error}"),
            Self::Invalid(what) => write!(f, "slurry: {what}"),
            Self::SourceMismatch { component } => {
                write!(f, "slurry: incompatible density sources for {component:?}")
            }
            Self::WaterRemovalRequired { removal } => write!(
                f,
                "slurry: target requires removing {} kg/s of liquid",
                removal.value()
            ),
            Self::UndefinedComposition => write!(f, "slurry: zero flow has no composition"),
        }
    }
}
impl std::error::Error for SlurryError {}
impl From<MatDbError> for SlurryError {
    fn from(value: MatDbError) -> Self {
        Self::Material(value)
    }
}

type Result<T> = core::result::Result<T, SlurryError>;

fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(SlurryError::Invalid(message))
    }
}

fn positive(value: f64, message: &'static str) -> Result<f64> {
    require(value.is_finite() && value > 0.0, message)?;
    Ok(value)
}

/// Immutable density sample with its real pack, query, evidence and usage receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct DensitySource {
    pack_id: String,
    pack_hash: ContentHash,
    query: QueryPoint,
    answer: MaterialAnswer,
}

impl DensitySource {
    /// Query and verify a positive density from an admitted normalized pack.
    ///
    /// The property must be `density`, `density-reported-range-minimum`, or
    /// `density-reported-range-maximum`. A source range never becomes a midpoint
    /// or a probability distribution: the caller must select the named endpoint.
    pub fn resolve(
        pack: &NormalizedPack,
        property: &str,
        query: &QueryPoint,
        policy: SelectionPolicy,
    ) -> Result<Self> {
        require(
            matches!(
                property,
                "density" | "density-reported-range-minimum" | "density-reported-range-maximum"
            ),
            "select an explicit density property or reported-range endpoint",
        )?;
        let answer = pack.claims().query_typed(
            &PropertyKey::with_quantity(property, QuantitySpec::dimensional(DENSITY)),
            query,
            policy,
        )?;
        pack.claims().verify_receipt(&answer.receipt)?;
        positive(
            answer.evidence.value.value,
            "source density must be finite and positive",
        )?;
        Ok(Self {
            pack_id: pack.pack_id().to_owned(),
            pack_hash: pack.content_hash(),
            query: query.clone(),
            answer,
        })
    }

    /// Selected density, in kg/m³, without invented uncertainty propagation.
    #[must_use]
    pub fn density(&self) -> Density {
        Density::new(self.answer.evidence.value.value)
    }

    /// Source-manifest pack identifier.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }

    /// Content identity of the complete normalized pack.
    #[must_use]
    pub const fn pack_hash(&self) -> ContentHash {
        self.pack_hash
    }

    /// Original typed conditions; unknown reference temperature stays unknown.
    #[must_use]
    pub const fn query(&self) -> &QueryPoint {
        &self.query
    }

    /// Original evidence, uncertainty and verified property-usage receipt.
    #[must_use]
    pub const fn answer(&self) -> &MaterialAnswer {
        &self.answer
    }
}

/// One identified constituent's flow and selected density evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct SlurryComponent {
    id: String,
    mass_flow: MassFlowRate,
    density_source: DensitySource,
}

impl SlurryComponent {
    /// Assign a stable nonblank identity and finite nonnegative flow in kg/s.
    ///
    /// Identity is supplied by the caller, not inferred from a mineral name.
    /// A component's source remains unchanged by subsequent flow operations.
    pub fn new(
        id: impl Into<String>,
        mass_flow: MassFlowRate,
        density_source: DensitySource,
    ) -> Result<Self> {
        let id = id.into();
        require(
            !id.trim().is_empty() && id.len() <= MAX_COMPONENT_ID_BYTES,
            "component identity must contain 1 to 256 bytes and cannot be blank",
        )?;
        require(
            mass_flow.value().is_finite() && mass_flow.value() >= 0.0,
            "component mass flow must be finite and nonnegative",
        )?;
        let component = Self {
            id,
            mass_flow: MassFlowRate::new(if mass_flow.value() == 0.0 {
                0.0
            } else {
                mass_flow.value()
            }),
            density_source,
        };
        component.volume_flow()?;
        Ok(component)
    }

    /// Caller-defined material identity.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Constituent mass rate in kg/s.
    #[must_use]
    pub const fn mass_flow(&self) -> MassFlowRate {
        self.mass_flow
    }

    /// Selected source density and its complete original context.
    #[must_use]
    pub const fn density_source(&self) -> &DensitySource {
        &self.density_source
    }

    fn volume_flow(&self) -> Result<f64> {
        if self.mass_flow.value() == 0.0 {
            Ok(0.0)
        } else {
            positive(
                self.mass_flow.value() / self.density_source.density().value(),
                "constituent volume flow is unrepresentable",
            )
        }
    }

    fn with_flow(&self, mass_flow: f64) -> Result<Self> {
        Self::new(
            self.id.clone(),
            MassFlowRate::new(mass_flow),
            self.density_source.clone(),
        )
    }
}

/// Dry-solid mass divided by total slurry mass; distinct from volume fraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidsMassFraction(Dimensionless);

impl SolidsMassFraction {
    /// Admit a finite fraction in the closed interval [0, 1].
    pub fn new(value: Dimensionless) -> Result<Self> {
        require(
            value.value().is_finite() && (0.0..=1.0).contains(&value.value()),
            "solids mass fraction must be finite and between zero and one",
        )?;
        Ok(Self(Dimensionless::new(value.value().abs())))
    }

    /// Unitless numeric fraction, not percentage.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Typed dimensionless value with the mass-basis wrapper removed explicitly.
    #[must_use]
    pub const fn quantity(self) -> Dimensionless {
        self.0
    }
}

/// Solid constituent volume divided by total additive volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidsVolumeFraction(Dimensionless);

impl SolidsVolumeFraction {
    /// Unitless numeric fraction, not percentage or a solids mass fraction.
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0.value()
    }

    /// Typed dimensionless value with the volume-basis wrapper removed explicitly.
    #[must_use]
    pub const fn quantity(self) -> Dimensionless {
        self.0
    }
}

/// Aggregate reference properties; all extensive rates are finite and nonnegative.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlurryState {
    /// Sum of all solid constituent mass rates, in kg/s.
    pub dry_solids_mass_flow: MassFlowRate,
    /// Liquid carrier mass rate, in kg/s.
    pub liquid_mass_flow: MassFlowRate,
    /// Combined solids and liquid mass rate, in kg/s.
    pub total_mass_flow: MassFlowRate,
    /// Sum of solid constituent `mass_flow / density`, in m³/s.
    pub solid_volume_flow: VolumetricFlowRate,
    /// Liquid carrier volume rate, in m³/s.
    pub liquid_volume_flow: VolumetricFlowRate,
    /// Additive solid and liquid volume rate, in m³/s.
    pub total_volume_flow: VolumetricFlowRate,
    /// Total mass/volume ratio, in kg/m³; `None` for zero total flow.
    pub density: Option<Density>,
    /// Dry-solid mass fraction; `None` for zero total flow.
    pub solids_mass_fraction: Option<SolidsMassFraction>,
    /// Additive solid volume fraction; `None` for zero total flow.
    pub solids_volume_fraction: Option<SolidsVolumeFraction>,
}

/// One finite stream under the explicit additive-volume, no-gas assumption.
#[derive(Debug, Clone, PartialEq)]
pub struct SlurryStream {
    solids: Vec<SlurryComponent>,
    liquid: SlurryComponent,
    state: SlurryState,
}

/// Successful addition of the stream's original liquid carrier.
#[derive(Debug, Clone, PartialEq)]
pub struct SlurryDilution {
    /// New stream with conserved solid component flows and retained sources.
    pub stream: SlurryStream,
    /// Nonnegative amount of liquid added, in kg/s.
    pub added_liquid_flow: MassFlowRate,
}

impl SlurryStream {
    /// Construct an explicitly additive-volume stream, without gas or contraction.
    ///
    /// Up to 128 uniquely identified solids are sorted into canonical ID order.
    /// The liquid must have a distinct identity. A zero-flow stream is admitted
    /// with `None` for density and composition; no pure-water value is invented.
    pub fn new_additive_volume(
        mut solids: Vec<SlurryComponent>,
        liquid: SlurryComponent,
    ) -> Result<Self> {
        require(
            solids.len() <= MAX_SOLID_COMPONENTS,
            "too many solid components",
        )?;
        solids.sort_by(|left, right| left.id.cmp(&right.id));
        require(
            solids.windows(2).all(|pair| pair[0].id != pair[1].id)
                && solids.iter().all(|solid| solid.id != liquid.id),
            "solid and liquid component identities must be distinct",
        )?;
        let solid_mass = sum(solids.iter().map(|solid| solid.mass_flow.value()))?;
        let liquid_mass = liquid.mass_flow.value();
        let total_mass = sum([solid_mass, liquid_mass])?;
        let solid_volume = sum(solids
            .iter()
            .map(SlurryComponent::volume_flow)
            .collect::<Result<Vec<_>>>()?)?;
        let liquid_volume = liquid.volume_flow()?;
        let total_volume = sum([solid_volume, liquid_volume])?;
        let (density, mass_fraction, volume_fraction) = if total_mass == 0.0 {
            (None, None, None)
        } else {
            (
                Some(Density::new(positive(
                    total_mass / total_volume,
                    "slurry density is unrepresentable",
                )?)),
                Some(SolidsMassFraction(Dimensionless::new(fraction(
                    solid_mass,
                    liquid_mass,
                    total_mass,
                )?))),
                Some(SolidsVolumeFraction(Dimensionless::new(fraction(
                    solid_volume,
                    liquid_volume,
                    total_volume,
                )?))),
            )
        };
        Ok(Self {
            solids,
            liquid,
            state: SlurryState {
                dry_solids_mass_flow: MassFlowRate::new(solid_mass),
                liquid_mass_flow: MassFlowRate::new(liquid_mass),
                total_mass_flow: MassFlowRate::new(total_mass),
                solid_volume_flow: VolumetricFlowRate::new(solid_volume),
                liquid_volume_flow: VolumetricFlowRate::new(liquid_volume),
                total_volume_flow: VolumetricFlowRate::new(total_volume),
                density,
                solids_mass_fraction: mass_fraction,
                solids_volume_fraction: volume_fraction,
            },
        })
    }

    /// Canonically ordered solid components, including any zero-flow constituents.
    #[must_use]
    pub fn solids(&self) -> &[SlurryComponent] {
        &self.solids
    }

    /// Liquid carrier with its original source density and query context.
    #[must_use]
    pub const fn liquid(&self) -> &SlurryComponent {
        &self.liquid
    }

    /// Evaluated finite rates and optional intensive properties.
    #[must_use]
    pub const fn state(&self) -> SlurryState {
        self.state
    }

    /// Add only the original liquid until the requested solids mass fraction.
    ///
    /// The result retains all source evidence. A target above the current mass
    /// fraction reports the required liquid removal without performing it.
    /// Positive solid flow cannot reach exactly zero fraction with finite water.
    pub fn dilute_to_mass_fraction(&self, target: SolidsMassFraction) -> Result<SlurryDilution> {
        let current = self
            .state
            .solids_mass_fraction
            .ok_or(SlurryError::UndefinedComposition)?;
        if target == current {
            return Ok(SlurryDilution {
                stream: self.clone(),
                added_liquid_flow: MassFlowRate::new(0.0),
            });
        }
        let solid = self.state.dry_solids_mass_flow.value();
        require(solid > 0.0, "liquid addition cannot create a solid phase")?;
        require(
            target.value() > 0.0,
            "positive solids cannot reach zero mass fraction with finite liquid",
        )?;
        // Multiplication first avoids an avoidable ratio overflow for very small
        // solid flows. A zero result for a positive product is refused below.
        let numerator = solid * (1.0 - target.value());
        let required_liquid = numerator / target.value();
        require(
            required_liquid.is_finite()
                && (target.value() == 1.0 || (numerator > 0.0 && required_liquid > 0.0)),
            "dilution liquid flow is unrepresentable",
        )?;
        let delta = required_liquid - self.liquid.mass_flow.value();
        if delta < 0.0 {
            return Err(SlurryError::WaterRemovalRequired {
                removal: MassFlowRate::new(-delta),
            });
        }
        positive(delta, "requested dilution change is unrepresentable")?;
        Ok(SlurryDilution {
            stream: Self::new_additive_volume(
                self.solids.clone(),
                self.liquid.with_flow(required_liquid)?,
            )?,
            added_liquid_flow: MassFlowRate::new(delta),
        })
    }

    /// Divide every constituent by the same first-outlet fraction in [0, 1].
    ///
    /// This is a nonselective splitter, not separation or recovery. The second
    /// outlet uses the constituent remainder, preserving mass to floating-point
    /// rounding. An endpoint produces one valid zero-flow stream.
    pub fn split(&self, first_fraction: Dimensionless) -> Result<(Self, Self)> {
        let share = SolidsMassFraction::new(first_fraction)?.value();
        let split_component = |component: &SlurryComponent| -> Result<_> {
            let original = component.mass_flow.value();
            let first = original * share;
            let second = original - first;
            require(
                original == 0.0 || share == 0.0 || share == 1.0 || (first > 0.0 && second > 0.0),
                "positive split outlet is unrepresentable",
            )?;
            Ok((component.with_flow(first)?, component.with_flow(second)?))
        };
        let (first_solids, second_solids): (Vec<_>, Vec<_>) = self
            .solids
            .iter()
            .map(split_component)
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .unzip();
        let (first_liquid, second_liquid) = split_component(&self.liquid)?;
        Ok((
            Self::new_additive_volume(first_solids, first_liquid)?,
            Self::new_additive_volume(second_solids, second_liquid)?,
        ))
    }

    /// Conservatively combine 1 to 128 compatible streams in canonical order.
    ///
    /// All carriers must have the same identity and density source/query. Solids
    /// with the same identity must also have identical source/query evidence.
    /// Distinct solid identities stay distinct. Input order cannot change a sum;
    /// this does not perform temperature equilibration or reconcile source data.
    pub fn mix(streams: &[&Self]) -> Result<Self> {
        require(
            !streams.is_empty() && streams.len() <= MAX_MIXING_STREAMS,
            "mixing requires between one and 128 streams",
        )?;
        let carrier = &streams[0].liquid;
        let mut constituents: BTreeMap<&str, (&SlurryComponent, Vec<f64>)> = BTreeMap::new();
        let mut liquid_rates = Vec::with_capacity(streams.len());
        for stream in streams {
            require(
                stream.liquid.id == carrier.id,
                "mixed streams must use the same liquid carrier identity",
            )?;
            same_source(carrier, &stream.liquid)?;
            liquid_rates.push(stream.liquid.mass_flow.value());
            for solid in &stream.solids {
                if let Some((original, rates)) = constituents.get_mut(solid.id.as_str()) {
                    same_source(original, solid)?;
                    rates.push(solid.mass_flow.value());
                } else {
                    require(
                        constituents.len() < MAX_SOLID_COMPONENTS,
                        "mixed stream has too many solid components",
                    )?;
                    constituents.insert(&solid.id, (solid, vec![solid.mass_flow.value()]));
                }
            }
        }
        let solids = constituents
            .into_values()
            .map(|(component, rates)| component.with_flow(sum(rates)?))
            .collect::<Result<Vec<_>>>()?;
        Self::new_additive_volume(solids, carrier.with_flow(sum(liquid_rates)?)?)
    }
}

fn same_source(left: &SlurryComponent, right: &SlurryComponent) -> Result<()> {
    if left.density_source == right.density_source {
        Ok(())
    } else {
        Err(SlurryError::SourceMismatch {
            component: left.id.clone(),
        })
    }
}

// Nonnegative sums, sorted by magnitude, with Neumaier compensation. Inputs are
// bounded by the stream/component limits before reaching this private helper.
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut values: Vec<_> = values.into_iter().collect();
    values.sort_by(f64::total_cmp);
    let mut total = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        let next = total + value;
        require(next.is_finite(), "combined flow overflows")?;
        correction += if total >= value {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    let result = total + correction;
    require(
        result.is_finite() && result >= 0.0,
        "combined flow is unrepresentable",
    )?;
    Ok(result)
}

fn fraction(solid: f64, liquid: f64, total: f64) -> Result<f64> {
    if solid == 0.0 {
        return Ok(0.0);
    }
    if liquid == 0.0 {
        return Ok(1.0);
    }
    let result = solid / total;
    require(
        result > 0.0 && result < 1.0,
        "positive phase fraction is unrepresentable",
    )?;
    Ok(result)
}

#[cfg(test)]
#[path = "slurry/tests.rs"]
mod tests;
