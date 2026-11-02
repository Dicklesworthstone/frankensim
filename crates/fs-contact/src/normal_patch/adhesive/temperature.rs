//! Original-source temperature changes of normal adhesive equilibrium.
//!
//! Both bulk elastic states and the interface work of adhesion are refreshed
//! independently. The original radius, attraction range, approximation and
//! applicability envelope remain fixed. This is an endpoint parameter change:
//! it neither supplies thermal expansion nor guarantees support or stability
//! along an unspecified temperature path. Losing the selected equilibrium
//! refuses; branch jumps and their energy require the existing event owners.
//!
//! The signed contact-potential change is retained separately from mechanical
//! work by a held-load actuator. Their difference is parameter work under this
//! selected control, NOT friction heat, entropy production or a caloric law.

use fs_blake3::ContentHash;
use fs_matdb::MaterialCard;
use fs_material::state_point::{IsotropicSolidStatePoint, MaterialStatePointError};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::fixture::JkrLoadingError;
use super::{JkrBranch, JkrError, JkrResponse, JkrState, SourceJkrSphereContact};
use crate::interface_binding::temperature::bind_isotropic_elastic_interface_at_temperatures;
use crate::interface_binding::{
    BoundNormalInterfaceState, InterfaceBindingError, bind_normal_interface_state,
};

/// Independently specified positive absolute temperatures, in kelvin.
///
/// Each value updates the original canonical `T` coordinate in its own complete
/// source query. Neither bulk temperature is replaced by the interface value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTemperatures {
    /// Temperature of the original ordered surface-A solid.
    pub surface_a_k: f64,
    /// Temperature of the original ordered surface-B solid.
    pub surface_b_k: f64,
    /// Temperature of the original ordered interface claims.
    pub interface_k: f64,
}

/// The physical normal control retained during a temperature change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AdhesiveThermalControl {
    /// Signed indentation is held; the normal actuator does no work.
    HeldIndentation,
    /// The stable contacting equilibrium follows one held applied load.
    HeldLoad {
        /// Constant signed compressive-positive actuator load, N.
        force_n: f64,
    },
}

/// Immutable refreshed owner, endpoint responses and signed work accounting.
///
/// The old owner is identified without cloning its cards a second time. The
/// response and refreshed model are read-only so consumers cannot import
/// edited force, area, potential or source identity into a refresh result.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceAdhesiveThermalRefresh<M, R> {
    original_model_identity: ContentHash,
    source: M,
    before: R,
    response: R,
    control: AdhesiveThermalControl,
    potential_change_j: f64,
    mechanical_work_j: f64,
    parameter_work_j: f64,
}

impl<M, R> SourceAdhesiveThermalRefresh<M, R> {
    /// Exact original model, including original source, geometry and envelope.
    #[must_use]
    pub const fn original_model_identity(&self) -> ContentHash {
        self.original_model_identity
    }
    /// Same normal-law family rebuilt from the refreshed original sources.
    #[must_use]
    pub const fn source(&self) -> &M {
        &self.source
    }
    /// Move the refreshed owner into an accepted native state without cloning
    /// its retained material and interface source bundles.
    #[must_use]
    pub fn into_source(self) -> M {
        self.source
    }
    /// Canonical original equilibrium under the selected control.
    #[must_use]
    pub const fn before(&self) -> &R {
        &self.before
    }
    /// Canonical refreshed equilibrium under the same control.
    #[must_use]
    pub const fn response(&self) -> &R {
        &self.response
    }
    /// Explicit held indentation or load, never inferred from an endpoint force.
    #[must_use]
    pub const fn control(&self) -> AdhesiveThermalControl {
        self.control
    }
    /// Actual signed change of elastic-plus-surface potential, J.
    #[must_use]
    pub const fn potential_change_j(&self) -> f64 {
        self.potential_change_j
    }
    /// Normal actuator work: zero at held indentation, F times delta(delta) at
    /// held load. This does not include tangential work or frictional heat.
    #[must_use]
    pub const fn mechanical_work_j(&self) -> f64 {
        self.mechanical_work_j
    }
    /// Potential change minus normal actuator work, J; may have either sign.
    /// A coupled thermal owner must supply its own thermodynamic interpretation.
    #[must_use]
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }

    fn new(
        original_model_identity: ContentHash,
        source: M,
        before: R,
        response: R,
        control: AdhesiveThermalControl,
        potential_change_j: f64,
        mechanical_work_j: f64,
    ) -> Result<Self, AdhesiveTemperatureError> {
        let parameter_work_j = potential_change_j - mechanical_work_j;
        if [potential_change_j, mechanical_work_j, parameter_work_j]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(AdhesiveTemperatureError::Invalid(
                "unrepresentable normal thermal parameter work",
            ));
        }
        Ok(Self {
            original_model_identity,
            source,
            before,
            response,
            control,
            potential_change_j,
            mechanical_work_j,
            parameter_work_j,
        })
    }
}

/// Original JKR branch and refreshed JKR equilibrium with complete sources.
pub type SourceJkrThermalRefresh =
    SourceAdhesiveThermalRefresh<SourceJkrSphereContact, JkrResponse>;

/// Atomic source, normal approximation, branch or work-accounting refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveTemperatureError {
    /// Exact original material/interface query refused.
    Material(MaterialStatePointError),
    /// Ordered bulk/interface source binding refused.
    Binding(InterfaceBindingError),
    /// Original JKR approximation or displacement equilibrium refused.
    Jkr(JkrError),
    /// Original JKR force-controlled equilibrium refused.
    JkrLoading(JkrLoadingError),
    /// Unsupported temperature ownership or nonfinite arithmetic.
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for AdhesiveTemperatureError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<InterfaceBindingError> for AdhesiveTemperatureError {
    fn from(error: InterfaceBindingError) -> Self {
        Self::Binding(error)
    }
}
impl From<JkrError> for AdhesiveTemperatureError {
    fn from(error: JkrError) -> Self {
        Self::Jkr(error)
    }
}
impl From<JkrLoadingError> for AdhesiveTemperatureError {
    fn from(error: JkrLoadingError) -> Self {
        Self::JkrLoading(error)
    }
}
impl core::fmt::Display for AdhesiveTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "normal adhesive temperature source: {error}"),
            Self::Binding(error) => write!(f, "normal adhesive temperature binding: {error}"),
            Self::Jkr(error) => write!(f, "normal adhesive temperature: {error}"),
            Self::JkrLoading(error) => write!(f, "normal adhesive temperature: {error}"),
            Self::Invalid(reason) => write!(f, "normal adhesive temperature: {reason}"),
        }
    }
}
impl std::error::Error for AdhesiveTemperatureError {}

fn temperature_point(
    original: fs_matdb::QueryPoint,
    temperature_k: f64,
) -> Result<fs_matdb::QueryPoint, AdhesiveTemperatureError> {
    let quantity = original.axis_quantities().get("T").copied();
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if !temperature_k.is_finite()
        || temperature_k <= 0.0
        || !original
            .axes()
            .get("T")
            .is_some_and(|t| t.is_finite() && *t > 0.0)
        || quantity
            .is_some_and(|q| q != absolute && q != QuantitySpec::dimensional(Temperature::DIMS))
    {
        return Err(AdhesiveTemperatureError::Invalid(
            "normal refresh requires the original absolute-temperature T coordinate",
        ));
    }
    original.with("T", temperature_k).map_err(|source| {
        MaterialStatePointError::Query {
            property: "normal adhesive temperature:T".into(),
            source,
        }
        .into()
    })
}

fn bulk_at_temperature(
    original: &IsotropicSolidStatePoint,
    card: &MaterialCard,
    temperature_k: f64,
) -> Result<IsotropicSolidStatePoint, AdhesiveTemperatureError> {
    if card.content_hash() != original.resolved().card_identity() {
        return Err(AdhesiveTemperatureError::Invalid(
            "normal refresh requires the original ordered bulk cards",
        ));
    }
    let point = temperature_point(original.resolved().query(), temperature_k)?;
    if original.resolved().state_coordinate("T") == Some(temperature_k) {
        // Preserve the original receipt policy and identity for a no-op.
        Ok(original.clone())
    } else {
        Ok(original.at_query_point(card, &point)?)
    }
}

impl BoundNormalInterfaceState {
    /// Refresh both original bulk states and interface claims at independently
    /// specified temperatures, retaining exact selected claims and source roles.
    ///
    /// `material_cards` must be the original ordered immutable cards, including
    /// on a no-op temperature query. Only canonical absolute `T` changes; every
    /// other source coordinate and quantity schema is retained. A second source
    /// temperature, rate, strain or history coordinate refuses. The original
    /// isotropic parameterization, yield admission and coupled properties are
    /// preserved by the material owners. No fallback claim is selected.
    pub fn at_temperatures(
        &self,
        material_cards: [&MaterialCard; 2],
        temperatures: AdhesiveTemperatures,
    ) -> Result<Self, AdhesiveTemperatureError> {
        let elastic = self.elastic();
        if elastic
            .surface_a_state()
            .resolved()
            .properties()
            .iter()
            .chain(elastic.surface_b_state().resolved().properties())
            .chain(self.interface_state().properties())
            .any(|p| !p.has_only_fixed_external_conditions(Some("T"), false))
        {
            return Err(AdhesiveTemperatureError::Invalid(
                "normal refresh requires one owned temperature and fixed external conditions",
            ));
        }
        let a = bulk_at_temperature(
            elastic.surface_a_state(),
            material_cards[0],
            temperatures.surface_a_k,
        )?;
        let b = bulk_at_temperature(
            elastic.surface_b_state(),
            material_cards[1],
            temperatures.surface_b_k,
        )?;
        let original = self.interface_state();
        let point = temperature_point(original.query(), temperatures.interface_k)?;
        let interface = if original.state_coordinate("T") == Some(temperatures.interface_k) {
            original.clone()
        } else {
            original.at_query_point(elastic.interface().card(), &point)?
        };
        let refreshed =
            bind_isotropic_elastic_interface_at_temperatures(elastic.interface(), &a, &b, &point)?;
        Ok(bind_normal_interface_state(&refreshed, &interface)?)
    }
}

impl SourceJkrSphereContact {
    /// Rebuild the same JKR approximation from independently refreshed original
    /// sources. Radius, attraction range and every applicability limit stay fixed.
    /// This owner-only operation does not transfer accepted branch history.
    pub fn at_temperatures(
        &self,
        material_cards: [&MaterialCard; 2],
        temperatures: AdhesiveTemperatures,
    ) -> Result<Self, AdhesiveTemperatureError> {
        Ok(Self::new(
            self.normal_state()
                .at_temperatures(material_cards, temperatures)?,
            self.effective_radius_m(),
            self.interaction_range_m(),
            self.envelope(),
        )?)
    }

    /// Retain an original attached or detached branch at its exact signed
    /// indentation while the three temperatures change. A thermal pull-off
    /// cannot silently detach, nucleate, reuse an old radius or discard energy.
    /// An absent new attached equilibrium returns the original JKR refusal.
    pub fn refresh_at_held_indentation(
        &self,
        previous: JkrState,
        material_cards: [&MaterialCard; 2],
        temperatures: AdhesiveTemperatures,
    ) -> Result<SourceJkrThermalRefresh, AdhesiveTemperatureError> {
        if previous.model_identity != self.identity() {
            return Err(JkrError::StateMismatch.into());
        }
        let delta = previous.indentation_m;
        let before = match previous.branch {
            JkrBranch::Attached => self.evaluate_attached(delta)?,
            JkrBranch::Detached => {
                self.detached_state(delta)?;
                self.detached_response(delta)
            }
        };
        let source = self.at_temperatures(material_cards, temperatures)?;
        let response = match previous.branch {
            JkrBranch::Attached => source.evaluate_attached(delta)?,
            JkrBranch::Detached => {
                source.detached_state(delta)?;
                source.detached_response(delta)
            }
        };
        let change = response.reversible_energy_j - before.reversible_energy_j;
        SourceAdhesiveThermalRefresh::new(
            self.identity(),
            source,
            before,
            response,
            AdhesiveThermalControl::HeldIndentation,
            change,
            0.0,
        )
    }

    /// Refresh the stable force-controlled contacting branch at one held load.
    /// The result separates actual potential change from F times delta(delta).
    /// This explicit equilibrium query neither attaches a detached history nor
    /// claims the intervening thermal path; loads below either pull-off refuse.
    pub fn refresh_at_held_load(
        &self,
        force_n: f64,
        material_cards: [&MaterialCard; 2],
        temperatures: AdhesiveTemperatures,
    ) -> Result<SourceJkrThermalRefresh, AdhesiveTemperatureError> {
        let before = self.evaluate_force_controlled(force_n)?;
        let source = self.at_temperatures(material_cards, temperatures)?;
        let response = source.evaluate_force_controlled(force_n)?;
        let change = response.reversible_energy_j - before.reversible_energy_j;
        let work = force_n * (response.indentation_m - before.indentation_m);
        SourceAdhesiveThermalRefresh::new(
            self.identity(),
            source,
            before,
            response,
            AdhesiveThermalControl::HeldLoad { force_n },
            change,
            work,
        )
    }
}
