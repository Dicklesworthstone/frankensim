//! Original triplet-pair encounter coefficients and singlet branching.
//!
//! The source convention is `R = k N_t^2 / V` PAIR EVENTS per second, with
//! `k` in m^3/s and populations in individual-emitter counts. Each pair removes
//! two triplets. A coefficient in `dn_t/dt = -k_loss n_t^2` must be explicitly
//! converted to `k_loss/2` before becoming this claim; molar coefficients need
//! their own explicit conversion. No diffusion coefficient, spin factor or
//! material-name estimate is substituted for either original source role.
//!
//! The physical owner supplies its active volume and excitation energies.
//! This snapshot does not evolve populations or decide whether a productive
//! pair can supply the singlet energy. Those belong to the existing thermal
//! annihilation law and its finite-population coupling. Temperature refresh
//! preserves the original card, selected claims and all other query conditions.

use super::{
    ContentHash, Dims, DomainHasher, MaterialCard, MaterialPropertySelection, PhosphorescenceError,
    QuantitySpec, QueryPoint, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermal::luminescence::transfer::annihilation::TripletAnnihilation;

/// SI m^3/s pair-event coefficient, with `R = k N_t^2 / V`.
pub const TRIPLET_ANNIHILATION_PAIR_RATE_PROPERTY: &str = "luminescence_triplet_pair_event_rate";
/// Probability that one annihilated triplet pair produces one singlet.
pub const TRIPLET_ANNIHILATION_SINGLET_YIELD_PROPERTY: &str =
    "luminescence_triplet_pair_singlet_yield";
/// Cubic metres per second, not a first-order removal rate or molar coefficient.
pub const TRIPLET_ANNIHILATION_PAIR_RATE_DIMS: Dims = Dims([3, 0, -1, 0, 0, 0]);
const TEMPERATURE_QUANTITY: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::AbsoluteTemperature,
    ValueForm::Static,
));

/// Immutable source bundle at an exact physical condition. Source identities
/// retain attribution and support; they do not establish empirical accuracy.
#[derive(Clone, Debug, PartialEq)]
pub struct TripletAnnihilationStatePoint {
    resolved: ResolvedMaterialStatePoint,
    temperature_axis: Option<String>,
    temperature_k: Option<f64>,
    pair_rate_m3_s: f64,
    singlet_yield: f64,
    identity: ContentHash,
}

impl TripletAnnihilationStatePoint {
    /// Resolve both explicit roles. A zero encounter rate or zero productive
    /// yield is valid; missing and ambiguous claims are never defaulted.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        temperature_axis: Option<&str>,
    ) -> Result<Self, PhosphorescenceError> {
        let requirements = [
            ScalarPropertyRequirement::try_new(
                TRIPLET_ANNIHILATION_PAIR_RATE_PROPERTY,
                TRIPLET_ANNIHILATION_PAIR_RATE_DIMS,
                ScalarAdmissibility::NonNegative,
            )?,
            ScalarPropertyRequirement::try_new(
                TRIPLET_ANNIHILATION_SINGLET_YIELD_PROPERTY,
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )?,
        ];
        Self::try_from_material(
            resolve_material_state_point(card, point, &requirements, selection)?,
            temperature_axis,
        )
    }

    /// Adopt a complete source bundle, rechecking units, branching probability
    /// and scalar support even when its original requirements were permissive.
    /// Only the named absolute-temperature coordinate is refreshed implicitly.
    /// Other supported external conditions stay fixed; field/history, frequency
    /// and unowned mechanical response refuse.
    pub fn try_from_material(
        resolved: ResolvedMaterialStatePoint,
        temperature_axis: Option<&str>,
    ) -> Result<Self, PhosphorescenceError> {
        let temperature_k = temperature_axis
            .map(|axis| {
                let value = resolved.state_coordinate(axis);
                if resolved.axis_quantities().get(axis) != Some(&TEMPERATURE_QUANTITY)
                    || !value.is_some_and(|t| t.is_finite() && t > 0.0)
                {
                    return Err(invalid(
                        "owned temperature requires positive explicitly typed kelvin",
                    ));
                }
                Ok(value.expect("admitted temperature coordinate"))
            })
            .transpose()?;
        let coordinate = [(temperature_axis.unwrap_or(""), TEMPERATURE_QUANTITY)];
        let owned = &coordinate[..usize::from(temperature_axis.is_some())];
        let parameter = |name: &'static str, dims| {
            let error = || PhosphorescenceError::Property {
                name,
                reason: "requires the original nonnegative scalar quantity at owned temperature and fixed external conditions",
            };
            let property = resolved.property(name).ok_or_else(error)?;
            let requirement = property.requirement();
            let value = property.value_si();
            if requirement.quantity() != QuantitySpec::dimensional(dims)
                || requirement.hardness_test().is_some()
                || requirement.elastic_component().is_some()
                || requirement.strain_component().is_some()
                || requirement.stress_component().is_some()
                || !value.is_finite()
                || value < 0.0
                || !property.has_only_owned_coordinates_and_fixed_external_conditions(owned)
            {
                return Err(error());
            }
            Ok(value)
        };
        let pair_rate_m3_s = parameter(
            TRIPLET_ANNIHILATION_PAIR_RATE_PROPERTY,
            TRIPLET_ANNIHILATION_PAIR_RATE_DIMS,
        )?;
        let singlet_yield = parameter(TRIPLET_ANNIHILATION_SINGLET_YIELD_PROPERTY, Dims::NONE)?;
        if singlet_yield > 1.0 {
            return Err(invalid("singlet probability per pair must lie in [0, 1]"));
        }
        let mut hash = DomainHasher::new("org.frankensim.material.triplet-pair-annihilation.v1");
        hash.update(resolved.identity().as_bytes());
        hash.update(&[u8::from(temperature_axis.is_some())]);
        if let Some(axis) = temperature_axis {
            hash.update(&(axis.len() as u64).to_le_bytes());
            hash.update(axis.as_bytes());
        }
        Ok(Self {
            resolved,
            temperature_axis: temperature_axis.map(str::to_owned),
            temperature_k,
            pair_rate_m3_s,
            singlet_yield,
            identity: hash.finalize(),
        })
    }

    /// Complete original card, query, requirements and selected claim receipts.
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }
    /// Identity of the complete source snapshot and its temperature ownership.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Original m^3/s pair-event coefficient, not the twice-larger loss coefficient.
    pub const fn pair_rate_m3_s(&self) -> f64 {
        self.pair_rate_m3_s
    }
    /// Original probability of one singlet per removed triplet pair.
    pub const fn singlet_yield(&self) -> f64 {
        self.singlet_yield
    }
    /// Explicitly owned absolute-temperature coordinate, when supplied.
    pub fn temperature_axis(&self) -> Option<&str> {
        self.temperature_axis.as_deref()
    }
    /// Actual temperature at the owned coordinate, in kelvin.
    pub const fn temperature_k(&self) -> Option<f64> {
        self.temperature_k
    }

    /// Lower into the existing physical law with the caller's explicit active
    /// volume. No source query or population evolution occurs here.
    pub fn law(&self, active_volume_m3: f64) -> Result<TripletAnnihilation, PhosphorescenceError> {
        Ok(TripletAnnihilation::try_new(
            self.pair_rate_m3_s,
            active_volume_m3,
            self.singlet_yield,
        )?)
    }

    /// Refresh the complete point with the same original card and selected claims.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, PhosphorescenceError> {
        Self::try_from_material(
            self.resolved.at_query_point(card, point)?,
            self.temperature_axis.as_deref(),
        )
    }

    /// Refresh only the owned temperature. `Some(T)` is required exactly when
    /// an axis was bound; a caller-supplied temperature is never ignored.
    pub fn at_temperature(
        &self,
        card: &MaterialCard,
        temperature_k: Option<f64>,
    ) -> Result<Self, PhosphorescenceError> {
        if temperature_k.is_some() != self.temperature_axis.is_some()
            || temperature_k.is_some_and(|t| !t.is_finite() || t <= 0.0)
        {
            return Err(invalid(
                "positive temperature must match the explicitly owned axis",
            ));
        }
        let point = if let (Some(axis), Some(value)) = (&self.temperature_axis, temperature_k) {
            self.resolved
                .query()
                .with(axis.as_str(), value)
                .map_err(|_| invalid("invalid annihilation temperature coordinate"))?
        } else {
            self.resolved.query().clone()
        };
        self.at_query_point(card, &point)
    }
}

fn invalid(reason: &'static str) -> PhosphorescenceError {
    PhosphorescenceError::Property {
        name: "triplet-pair annihilation",
        reason,
    }
}
