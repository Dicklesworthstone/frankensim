//! Bind immutable material/interface identities into executable contact inputs.
//!
//! `fs-matdb` owns ordered surface identity and evidence-bearing property
//! claims; `fs-tribo` owns dependency-light dry constitutive leaves. This
//! module is the L3 bridge between them. It selects only an explicitly offered,
//! versioned model card and never infers a law from material names or upgrades
//! claim authority. A successful identity bind alone is not evidence that the
//! card contains every property a chosen law needs.

use core::fmt;
use std::collections::BTreeMap;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, InterfaceSystemCard, MaterialStateId};
use fs_material::state_point::{
    IsotropicSolidStatePoint, MaterialStatePointError, ResolvedInterfaceStatePoint,
};
use fs_qty::Dims;
use fs_tribo::{
    InputAuthority, InterfaceMedium, InterfaceSystemRef, surface_excitation::UniformSurfaceTrace,
};

/// Canonical ordered-interface property consumed by adhesive normal laws.
pub const ADHESION_ENERGY_PROPERTY: &str = "adhesion-energy";
/// SI dimensions of adhesion energy per area, J/m2 = kg/s2.
pub const ADHESION_ENERGY_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);

/// Stateless nonadhesive Hertz model-card law identifier.
pub const NORMAL_HERTZ_LAW_ID: &str = "fs-contact.normal-hertz";
/// Sphere-only Hunt--Crossley model-card law identifier.
pub const NORMAL_HUNT_CROSSLEY_SPHERE_LAW_ID: &str = "fs-contact.normal-hunt-crossley-sphere";
/// Executable model-card law version understood by this bridge.
pub const NORMAL_CONTACT_LAW_VERSION: u32 = 1;
/// Stateless normal-law state schema understood by this bridge.
pub const NORMAL_CONTACT_STATE_SCHEMA_VERSION: u32 = 1;

const CHARACTERISTIC_RATE_PARAMETER: &str = "characteristic-rate";
const MAX_PATCH_TO_RADIUS_PARAMETER: &str = "max-patch-to-radius";
const MAX_STRAIN_PARAMETER: &str = "max-strain";
const MAX_PATCH_TO_DEPTH_PARAMETER: &str = "max-patch-to-depth";
const MAX_PATCH_TO_LAYER_PARAMETER: &str = "max-patch-to-layer";
const MAX_PRESSURE_TO_YIELD_PARAMETER: &str = "max-pressure-to-yield";
const MAX_RATE_RATIO_PARAMETER: &str = "max-rate-ratio";
const HUNT_CROSSLEY_DISSIPATION_PARAMETER: &str = "dissipation";
const SPEED_DIMS: Dims = Dims([1, 0, -1, 0, 0, 0]);
const INVERSE_SPEED_DIMS: Dims = Dims([-1, 0, 1, 0, 0, 0]);

/// One ordered dry interface identity reconstructed from its immutable card.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundDryInterfaceSystem {
    card: InterfaceSystemCard,
    card_identity: ContentHash,
    surface_materials: [MaterialStateId; 2],
    texture_frame_ids: [String; 2],
    interface: InterfaceSystemRef,
}

/// Ordered isotropic elastic properties admitted against one dry interface.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundIsotropicElasticInterface {
    interface: BoundDryInterfaceSystem,
    state_point: Vec<(String, f64)>,
    surface_a_state: IsotropicSolidStatePoint,
    surface_b_state: IsotropicSolidStatePoint,
    reduced_modulus_pa: f64,
    limiting_yield_stress_pa: f64,
    identity: ContentHash,
}

/// Elastic bulk states plus the interface-state datum needed by the normal
/// contact ladder.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundNormalInterfaceState {
    elastic: BoundIsotropicElasticInterface,
    interface_state: ResolvedInterfaceStatePoint,
    adhesion_energy_j_per_m2: f64,
    identity: ContentHash,
}

/// Executable normal-response family selected by an immutable interface card.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BoundNormalContactLaw {
    /// Nonadhesive Hertz response; local geometry selects the analytic rung.
    ElasticHertz,
    /// Passive Hunt--Crossley augmentation, currently admitted only for a
    /// sphere/plane patch by the generic normal-law implementation.
    HuntCrossleySphere {
        /// Card-resolved dissipation coefficient [s/m].
        dissipation_s_per_m: f64,
    },
}

/// Explicit selection when an interface card carries normal-response models.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NormalContactModelSelection {
    /// Exactly one supported normal-response card must exist.
    SingleSupported,
    /// Select one exact immutable constitutive-model-card identity.
    Pinned(ContentHash),
}

/// Normal interface state bound to one executable constitutive-model card.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundNormalContactModel {
    normal_state: BoundNormalInterfaceState,
    model_card: ConstitutiveModelCard,
    law: BoundNormalContactLaw,
    characteristic_rate_m_per_s: f64,
    limits: crate::normal_patch::ApplicabilityLimits,
    identity: ContentHash,
}

impl BoundNormalContactModel {
    /// Complete ordered bulk/interface state consumed by this model.
    #[must_use]
    pub const fn normal_state(&self) -> &BoundNormalInterfaceState {
        &self.normal_state
    }

    /// Complete immutable model card, including validity and provenance.
    #[must_use]
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.model_card
    }

    /// Executable law family and any card-resolved coefficient.
    #[must_use]
    pub const fn law(&self) -> BoundNormalContactLaw {
        self.law
    }

    /// Card-resolved characteristic rate [m/s].
    #[must_use]
    pub const fn characteristic_rate_m_per_s(&self) -> f64 {
        self.characteristic_rate_m_per_s
    }

    /// Card-resolved applicability limits, including its temperature domain.
    #[must_use]
    pub const fn limits(&self) -> crate::normal_patch::ApplicabilityLimits {
        self.limits
    }

    /// Identity binding the ordered physical state and complete model card.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

impl BoundNormalInterfaceState {
    /// Ordered bulk elastic state and dry interface identity.
    #[must_use]
    pub const fn elastic(&self) -> &BoundIsotropicElasticInterface {
        &self.elastic
    }

    /// Complete resolved interface-property state and usage receipts.
    #[must_use]
    pub const fn interface_state(&self) -> &ResolvedInterfaceStatePoint {
        &self.interface_state
    }

    /// Resolved adhesion energy [J/m2]. Zero is an explicit card datum.
    #[must_use]
    pub const fn adhesion_energy_j_per_m2(&self) -> f64 {
        self.adhesion_energy_j_per_m2
    }

    /// Identity binding bulk states, interface state, and derived input.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

impl BoundIsotropicElasticInterface {
    /// Complete ordered dry interface identity.
    #[must_use]
    pub const fn interface(&self) -> &BoundDryInterfaceSystem {
        &self.interface
    }

    /// Complete common physical state point of both bulk materials.
    #[must_use]
    pub fn state_point(&self) -> &[(String, f64)] {
        &self.state_point
    }

    /// Look up one state coordinate by exact axis name.
    #[must_use]
    pub fn state_coordinate(&self, axis: &str) -> Option<f64> {
        self.state_point
            .binary_search_by(|(name, _)| name.as_str().cmp(axis))
            .ok()
            .map(|index| self.state_point[index].1)
    }

    /// Surface-A resolved material-state identity.
    #[must_use]
    pub const fn surface_a_state_identity(&self) -> ContentHash {
        self.surface_a_state.resolved().identity()
    }

    /// Surface-B resolved material-state identity.
    #[must_use]
    pub const fn surface_b_state_identity(&self) -> ContentHash {
        self.surface_b_state.resolved().identity()
    }

    /// Complete surface-A material state and property-use receipts.
    #[must_use]
    pub const fn surface_a_state(&self) -> &IsotropicSolidStatePoint {
        &self.surface_a_state
    }

    /// Complete surface-B material state and property-use receipts.
    #[must_use]
    pub const fn surface_b_state(&self) -> &IsotropicSolidStatePoint {
        &self.surface_b_state
    }

    /// Two-half-space Hertz reduced modulus [Pa].
    #[must_use]
    pub const fn reduced_modulus_pa(&self) -> f64 {
        self.reduced_modulus_pa
    }

    /// Smaller of the two state-point yield stresses [Pa].
    #[must_use]
    pub const fn limiting_yield_stress_pa(&self) -> f64 {
        self.limiting_yield_stress_pa
    }

    /// Identity binding ordered interface and both resolved material states.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

impl BoundDryInterfaceSystem {
    /// Complete immutable ordered interface card, including claims and models.
    #[must_use]
    pub const fn card(&self) -> &InterfaceSystemCard {
        &self.card
    }

    /// Content identity of the complete ordered fs-matdb interface card.
    #[must_use]
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    /// Bulk material-state identities, in ordered interface roles.
    #[must_use]
    pub const fn surface_materials(&self) -> &[MaterialStateId; 2] {
        &self.surface_materials
    }

    /// Opaque texture-frame identities, in ordered interface roles.
    #[must_use]
    pub fn texture_frame_ids(&self) -> [&str; 2] {
        [&self.texture_frame_ids[0], &self.texture_frame_ids[1]]
    }

    /// Dependency-light interface identity consumed by fs-contact/fs-tribo laws.
    #[must_use]
    pub const fn interface(&self) -> &InterfaceSystemRef {
        &self.interface
    }

    /// Verify that two sampled traces are exactly the card's ordered textures.
    pub fn admit_surface_traces(
        &self,
        surface_a: &UniformSurfaceTrace,
        surface_b: &UniformSurfaceTrace,
    ) -> Result<(), InterfaceBindingError> {
        for (role, expected, observed) in [
            (
                "surface_a",
                self.texture_frame_ids[0].as_str(),
                surface_a.texture_frame_id(),
            ),
            (
                "surface_b",
                self.texture_frame_ids[1].as_str(),
                surface_b.texture_frame_id(),
            ),
        ] {
            if expected != observed {
                return Err(InterfaceBindingError::TextureFrameMismatch {
                    role,
                    expected: expected.to_owned(),
                    observed: observed.to_owned(),
                });
            }
        }
        Ok(())
    }
}

/// Bind one complete ordered fs-matdb card to the dry-contact identity seam.
///
/// The card must declare the exact medium token `dry`; gas/liquid films are
/// different mechanism owners and cannot be re-labeled as dry tribology. The
/// caller supplies an authority ceiling because card identity alone does not
/// adjudicate the authority of every property claim it contains.
pub fn bind_dry_interface_system_card(
    card: &InterfaceSystemCard,
    authority: InputAuthority,
) -> Result<BoundDryInterfaceSystem, InterfaceBindingError> {
    if card.medium() != "dry" {
        return Err(InterfaceBindingError::UnsupportedMedium {
            observed: card.medium().to_owned(),
        });
    }
    let card_identity = card.content_hash();
    let interface = InterfaceSystemRef::new(
        format!("fs-matdb/interface-system-card/{card_identity}"),
        card.history(),
        format!("fs-matdb/interface-system-card/{card_identity}"),
        authority,
        InterfaceMedium::Dry,
    )
    .map_err(|error| InterfaceBindingError::Tribo(error.to_string()))?;
    Ok(BoundDryInterfaceSystem {
        card: card.clone(),
        card_identity,
        surface_materials: [
            card.surface_a().material.clone(),
            card.surface_b().material.clone(),
        ],
        texture_frame_ids: [
            card.surface_a().texture_frame.clone(),
            card.surface_b().texture_frame.clone(),
        ],
        interface,
    })
}

/// Bind two resolved isotropic solid states to their exact ordered interface.
///
/// Both material-state ids and all state-point coordinates must agree with
/// the interface roles and with each other. This prevents a contact law from
/// combining, for example, room-temperature disc elasticity with hot-base
/// elasticity or a reversed material pair. The result derives only elastic
/// modulus and a conservative yield applicability scalar; it does not select
/// a normal law, friction model, adhesion model, or validity limit.
pub fn bind_isotropic_elastic_interface(
    interface: &BoundDryInterfaceSystem,
    surface_a: &IsotropicSolidStatePoint,
    surface_b: &IsotropicSolidStatePoint,
) -> Result<BoundIsotropicElasticInterface, InterfaceBindingError> {
    for (role, expected, observed) in [
        (
            "surface_a",
            &interface.surface_materials[0],
            surface_a.resolved().material(),
        ),
        (
            "surface_b",
            &interface.surface_materials[1],
            surface_b.resolved().material(),
        ),
    ] {
        if expected != observed {
            return Err(InterfaceBindingError::MaterialStateMismatch {
                role,
                expected: expected.clone(),
                observed: observed.clone(),
            });
        }
    }
    if !state_points_exact_eq(
        surface_a.resolved().query_point(),
        surface_b.resolved().query_point(),
    ) {
        return Err(InterfaceBindingError::StatePointMismatch);
    }
    let reduced = surface_a
        .reduced_modulus_against(surface_b)
        .map_err(InterfaceBindingError::MaterialState)?;
    let limiting_yield_stress_pa = surface_a.yield_stress_pa().min(surface_b.yield_stress_pa());
    if !limiting_yield_stress_pa.is_finite() || limiting_yield_stress_pa <= 0.0 {
        return Err(InterfaceBindingError::InvalidDerived {
            quantity: "limiting_yield_stress_pa",
        });
    }
    let mut hasher = DomainHasher::new("org.frankensim.fs-contact.isotropic-interface.v1");
    hasher.update(interface.card_identity().as_bytes());
    hasher.update(reduced.surface_a_state_identity.as_bytes());
    hasher.update(reduced.surface_b_state_identity.as_bytes());
    hasher.update(&reduced.value_pa.to_bits().to_le_bytes());
    hasher.update(&limiting_yield_stress_pa.to_bits().to_le_bytes());
    Ok(BoundIsotropicElasticInterface {
        interface: interface.clone(),
        state_point: surface_a.resolved().query_point().to_vec(),
        surface_a_state: surface_a.clone(),
        surface_b_state: surface_b.clone(),
        reduced_modulus_pa: reduced.value_pa,
        limiting_yield_stress_pa,
        identity: hasher.finalize(),
    })
}

/// Bind the normal-contact interface datum to the exact ordered elastic state.
///
/// The interface property bundle must originate from the same immutable card
/// and the same complete state point as both bulk material bundles. The
/// current normal-law ladder requires an explicit nonnegative
/// [`ADHESION_ENERGY_PROPERTY`]; absence refuses rather than silently assuming
/// a nonadhesive interface.
pub fn bind_normal_interface_state(
    elastic: &BoundIsotropicElasticInterface,
    interface_state: &ResolvedInterfaceStatePoint,
) -> Result<BoundNormalInterfaceState, InterfaceBindingError> {
    if interface_state.card_identity() != elastic.interface.card_identity {
        return Err(InterfaceBindingError::InterfaceStateCardMismatch {
            expected: elastic.interface.card_identity,
            observed: interface_state.card_identity(),
        });
    }
    if interface_state.surface_materials() != elastic.interface.surface_materials() {
        return Err(InterfaceBindingError::InterfaceStateSurfaceMismatch);
    }
    if !state_points_exact_eq(interface_state.query_point(), elastic.state_point()) {
        return Err(InterfaceBindingError::StatePointMismatch);
    }
    let adhesion = interface_state.property(ADHESION_ENERGY_PROPERTY).ok_or(
        InterfaceBindingError::MissingInterfaceProperty {
            property: ADHESION_ENERGY_PROPERTY,
        },
    )?;
    if adhesion.requirement().dims() != ADHESION_ENERGY_DIMS {
        return Err(InterfaceBindingError::InterfacePropertyDimensionMismatch {
            property: ADHESION_ENERGY_PROPERTY,
            expected: ADHESION_ENERGY_DIMS,
            observed: adhesion.requirement().dims(),
        });
    }
    let adhesion_energy_j_per_m2 = adhesion.value_si();
    if !adhesion_energy_j_per_m2.is_finite() || adhesion_energy_j_per_m2 < 0.0 {
        return Err(InterfaceBindingError::InvalidDerived {
            quantity: "adhesion_energy_j_per_m2",
        });
    }
    let mut hasher = DomainHasher::new("org.frankensim.fs-contact.normal-interface-state.v1");
    hasher.update(elastic.identity.as_bytes());
    hasher.update(interface_state.identity().as_bytes());
    hasher.update(&adhesion_energy_j_per_m2.to_bits().to_le_bytes());
    Ok(BoundNormalInterfaceState {
        elastic: elastic.clone(),
        interface_state: interface_state.clone(),
        adhesion_energy_j_per_m2,
        identity: hasher.finalize(),
    })
}

/// Bind one interface-card constitutive model to the exact normal state.
///
/// Model choice, damping, rate scale, ratio limits, and temperature validity
/// come from the immutable ordered-interface card. No material name or absent
/// parameter implies a preset. Geometry-owned half-space and layer extents
/// remain assembly inputs because they are not material properties.
pub fn bind_normal_contact_model(
    normal_state: &BoundNormalInterfaceState,
    selection: NormalContactModelSelection,
) -> Result<BoundNormalContactModel, InterfaceBindingError> {
    let models = normal_state.elastic.interface.card.models();
    let supported: Vec<&ConstitutiveModelCard> = models
        .iter()
        .filter(|model| is_supported_normal_law(&model.law.0))
        .collect();
    let model = match selection {
        NormalContactModelSelection::SingleSupported => match supported.as_slice() {
            [model] => *model,
            [] => return Err(InterfaceBindingError::MissingNormalContactModel),
            many => {
                return Err(InterfaceBindingError::AmbiguousNormalContactModel {
                    candidates: many.iter().map(|model| model.content_hash()).collect(),
                });
            }
        },
        NormalContactModelSelection::Pinned(expected) => models
            .iter()
            .find(|model| model.content_hash() == expected)
            .ok_or(InterfaceBindingError::PinnedNormalContactModelMissing { expected })?,
    };
    if !is_supported_normal_law(&model.law.0) {
        return Err(InterfaceBindingError::UnsupportedNormalContactLaw {
            law: model.law.0.clone(),
            version: model.law_version,
        });
    }
    if model.law_version != NORMAL_CONTACT_LAW_VERSION
        || model.state_schema_version != NORMAL_CONTACT_STATE_SCHEMA_VERSION
        || model.initial_state != InitialStatePolicy::ZeroInternalState
    {
        return Err(InterfaceBindingError::UnsupportedNormalContactLaw {
            law: model.law.0.clone(),
            version: model.law_version,
        });
    }

    let point: BTreeMap<String, f64> = normal_state
        .interface_state
        .query_point()
        .iter()
        .cloned()
        .collect();
    if !model.validity.contains(&point) {
        return Err(InterfaceBindingError::NormalContactModelOutsideValidity {
            model: model.content_hash(),
        });
    }
    let (min_temperature_k, max_temperature_k) = model
        .validity
        .bound("T")
        .filter(|(minimum, maximum)| {
            minimum.is_finite() && maximum.is_finite() && *minimum > 0.0 && minimum <= maximum
        })
        .ok_or(InterfaceBindingError::MissingNormalContactTemperatureValidity)?;

    let characteristic_rate_m_per_s =
        required_model_parameter(model, CHARACTERISTIC_RATE_PARAMETER, SPEED_DIMS, true)?;
    let limits = crate::normal_patch::ApplicabilityLimits {
        max_patch_to_radius: required_model_parameter(
            model,
            MAX_PATCH_TO_RADIUS_PARAMETER,
            Dims::NONE,
            true,
        )?,
        max_strain: required_model_parameter(model, MAX_STRAIN_PARAMETER, Dims::NONE, true)?,
        max_patch_to_depth: required_model_parameter(
            model,
            MAX_PATCH_TO_DEPTH_PARAMETER,
            Dims::NONE,
            true,
        )?,
        max_patch_to_layer: required_model_parameter(
            model,
            MAX_PATCH_TO_LAYER_PARAMETER,
            Dims::NONE,
            true,
        )?,
        max_pressure_to_yield: required_model_parameter(
            model,
            MAX_PRESSURE_TO_YIELD_PARAMETER,
            Dims::NONE,
            true,
        )?,
        max_rate_ratio: required_model_parameter(
            model,
            MAX_RATE_RATIO_PARAMETER,
            Dims::NONE,
            true,
        )?,
        min_temperature_k,
        max_temperature_k,
    };
    let law = match model.law.0.as_str() {
        NORMAL_HERTZ_LAW_ID => BoundNormalContactLaw::ElasticHertz,
        NORMAL_HUNT_CROSSLEY_SPHERE_LAW_ID => BoundNormalContactLaw::HuntCrossleySphere {
            dissipation_s_per_m: required_model_parameter(
                model,
                HUNT_CROSSLEY_DISSIPATION_PARAMETER,
                INVERSE_SPEED_DIMS,
                false,
            )?,
        },
        _ => unreachable!("supported normal law filter and dispatch must agree"),
    };
    validate_parameter_roster(model, law)?;

    let mut hasher = DomainHasher::new("org.frankensim.fs-contact.normal-model-binding.v1");
    hasher.update(normal_state.identity.as_bytes());
    hasher.update(model.content_hash().as_bytes());
    Ok(BoundNormalContactModel {
        normal_state: normal_state.clone(),
        model_card: model.clone(),
        law,
        characteristic_rate_m_per_s,
        limits,
        identity: hasher.finalize(),
    })
}

fn is_supported_normal_law(law: &str) -> bool {
    matches!(
        law,
        NORMAL_HERTZ_LAW_ID | NORMAL_HUNT_CROSSLEY_SPHERE_LAW_ID
    )
}

fn required_model_parameter(
    model: &ConstitutiveModelCard,
    name: &'static str,
    expected_dims: Dims,
    strictly_positive: bool,
) -> Result<f64, InterfaceBindingError> {
    let parameter = model
        .parameters
        .get(name)
        .ok_or(InterfaceBindingError::MissingNormalContactModelParameter { parameter: name })?;
    if parameter.dims != expected_dims {
        return Err(
            InterfaceBindingError::NormalContactModelParameterDimensionMismatch {
                parameter: name,
                expected: expected_dims,
                observed: parameter.dims,
            },
        );
    }
    if !parameter.value.is_finite()
        || (strictly_positive && parameter.value <= 0.0)
        || (!strictly_positive && parameter.value < 0.0)
    {
        return Err(InterfaceBindingError::InvalidNormalContactModelParameter { parameter: name });
    }
    Ok(parameter.value)
}

fn validate_parameter_roster(
    model: &ConstitutiveModelCard,
    law: BoundNormalContactLaw,
) -> Result<(), InterfaceBindingError> {
    let expected = match law {
        BoundNormalContactLaw::ElasticHertz => &COMMON_NORMAL_PARAMETERS[..],
        BoundNormalContactLaw::HuntCrossleySphere { .. } => &HUNT_NORMAL_PARAMETERS[..],
    };
    if model.parameters.len() != expected.len()
        || !model
            .parameters
            .keys()
            .all(|parameter| expected.contains(&parameter.as_str()))
    {
        return Err(InterfaceBindingError::NormalContactModelParameterRosterMismatch);
    }
    Ok(())
}

const COMMON_NORMAL_PARAMETERS: [&str; 7] = [
    CHARACTERISTIC_RATE_PARAMETER,
    MAX_PATCH_TO_RADIUS_PARAMETER,
    MAX_STRAIN_PARAMETER,
    MAX_PATCH_TO_DEPTH_PARAMETER,
    MAX_PATCH_TO_LAYER_PARAMETER,
    MAX_PRESSURE_TO_YIELD_PARAMETER,
    MAX_RATE_RATIO_PARAMETER,
];
const HUNT_NORMAL_PARAMETERS: [&str; 8] = [
    CHARACTERISTIC_RATE_PARAMETER,
    MAX_PATCH_TO_RADIUS_PARAMETER,
    MAX_STRAIN_PARAMETER,
    MAX_PATCH_TO_DEPTH_PARAMETER,
    MAX_PATCH_TO_LAYER_PARAMETER,
    MAX_PRESSURE_TO_YIELD_PARAMETER,
    MAX_RATE_RATIO_PARAMETER,
    HUNT_CROSSLEY_DISSIPATION_PARAMETER,
];

fn state_points_exact_eq(left: &[(String, f64)], right: &[(String, f64)]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|((left_name, left_value), (right_name, right_value))| {
                left_name == right_name && left_value.to_bits() == right_value.to_bits()
            })
}

/// Structured refusal from the material-card/contact identity bridge.
#[derive(Clone, Debug, PartialEq)]
pub enum InterfaceBindingError {
    /// The selected contact mechanism is not the dry fs-tribo lane.
    UnsupportedMedium {
        /// Exact medium token retained by the ordered card.
        observed: String,
    },
    /// A sampled trace belongs to a different surface texture frame.
    TextureFrameMismatch {
        /// Ordered surface role.
        role: &'static str,
        /// Card-declared texture frame.
        expected: String,
        /// Trace-declared texture frame.
        observed: String,
    },
    /// A resolved material state was supplied in the wrong ordered role.
    MaterialStateMismatch {
        /// Ordered surface role.
        role: &'static str,
        /// Interface-card material state.
        expected: MaterialStateId,
        /// Offered resolved material state.
        observed: MaterialStateId,
    },
    /// The two bulk properties were evaluated at different physical points.
    StatePointMismatch,
    /// Interface properties came from a different immutable ordered card.
    InterfaceStateCardMismatch {
        /// Card bound to the elastic interface.
        expected: ContentHash,
        /// Card bound to the property state.
        observed: ContentHash,
    },
    /// Interface-state ordered material roles differed from the elastic bind.
    InterfaceStateSurfaceMismatch,
    /// No executable normal-response model card was supplied.
    MissingNormalContactModel,
    /// More than one supported normal-response model survived selection.
    AmbiguousNormalContactModel {
        /// Candidate immutable model identities.
        candidates: Vec<ContentHash>,
    },
    /// A pinned immutable model identity is absent from the interface card.
    PinnedNormalContactModelMissing {
        /// Requested complete model-card identity.
        expected: ContentHash,
    },
    /// The law id, law version, state schema, or initial-state policy is not executable here.
    UnsupportedNormalContactLaw {
        /// Offered stable law id.
        law: String,
        /// Offered law semantic version.
        version: u32,
    },
    /// The exact interface state lies outside the model-card validity box.
    NormalContactModelOutsideValidity {
        /// Complete model-card identity.
        model: ContentHash,
    },
    /// A finite positive `T` validity range is required by the normal law.
    MissingNormalContactTemperatureValidity,
    /// A required dimensioned law parameter was absent.
    MissingNormalContactModelParameter {
        /// Stable parameter name.
        parameter: &'static str,
    },
    /// A law parameter used the wrong SI dimensions.
    NormalContactModelParameterDimensionMismatch {
        /// Stable parameter name.
        parameter: &'static str,
        /// Executable-law dimensions.
        expected: Dims,
        /// Card-declared dimensions.
        observed: Dims,
    },
    /// A required law parameter was non-finite or outside its numerical domain.
    InvalidNormalContactModelParameter {
        /// Stable parameter name.
        parameter: &'static str,
    },
    /// Missing or unknown parameters made the executable mapping ambiguous.
    NormalContactModelParameterRosterMismatch,
    /// A consumer-required ordered-interface property was absent.
    MissingInterfaceProperty {
        /// Exact missing property key.
        property: &'static str,
    },
    /// A resolved interface property was admitted under the wrong dimensions.
    InterfacePropertyDimensionMismatch {
        /// Exact property key.
        property: &'static str,
        /// Consumer-required dimensions.
        expected: Dims,
        /// Dimensions used to resolve the property bundle.
        observed: Dims,
    },
    /// A material state-point operation refused.
    MaterialState(MaterialStatePointError),
    /// An ordered derived interface scalar was invalid.
    InvalidDerived {
        /// Stable derived quantity name.
        quantity: &'static str,
    },
    /// The dependency-light fs-tribo identity refused construction.
    Tribo(String),
}

impl fmt::Display for InterfaceBindingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for InterfaceBindingError {}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId,
        LawParameter, MaterialCard, PropertyClaim, PropertyKey, PropertyValue, Provenance,
        QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_interface_state_point, resolve_isotropic_solid_state_point,
    };
    use fs_qty::{Density, Dims, Pressure};
    use fs_tribo::surface_excitation::SurfaceTraceBoundary;

    use super::*;

    fn material(chemistry: &str) -> MaterialStateId {
        MaterialStateId {
            chemistry: chemistry.to_owned(),
            phase: "solid".to_owned(),
            process: "declared-test-state".to_owned(),
            revision: 0,
        }
    }

    fn normal_model(law: &str, dissipation_s_per_m: Option<f64>) -> ConstitutiveModelCard {
        let mut parameters = BTreeMap::new();
        for (name, dims, value) in [
            (CHARACTERISTIC_RATE_PARAMETER, SPEED_DIMS, 1.0),
            (MAX_PATCH_TO_RADIUS_PARAMETER, Dims::NONE, 0.1),
            (MAX_STRAIN_PARAMETER, Dims::NONE, 0.01),
            (MAX_PATCH_TO_DEPTH_PARAMETER, Dims::NONE, 0.1),
            (MAX_PATCH_TO_LAYER_PARAMETER, Dims::NONE, 0.1),
            (MAX_PRESSURE_TO_YIELD_PARAMETER, Dims::NONE, 0.2),
            (MAX_RATE_RATIO_PARAMETER, Dims::NONE, 0.1),
        ] {
            parameters.insert(name.to_owned(), LawParameter { value, dims });
        }
        if let Some(value) = dissipation_s_per_m {
            parameters.insert(
                HUNT_CROSSLEY_DISSIPATION_PARAMETER.to_owned(),
                LawParameter {
                    value,
                    dims: INVERSE_SPEED_DIMS,
                },
            );
        }
        ConstitutiveModelCard {
            law: LawId(law.to_owned()),
            law_version: NORMAL_CONTACT_LAW_VERSION,
            parameters,
            state_schema_version: NORMAL_CONTACT_STATE_SCHEMA_VERSION,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            sources: Vec::new(),
            provenance: Provenance {
                source: format!("synthetic {law} card"),
                license: "CC0-1.0".to_owned(),
                artifact: None,
            },
        }
    }

    fn card_with_models(medium: &str, models: Vec<ConstitutiveModelCard>) -> InterfaceSystemCard {
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(ADHESION_ENERGY_PROPERTY, ADHESION_ENERGY_DIMS),
                value: PropertyValue::Scalar {
                    value: 0.0,
                    dims: ADHESION_ENERGY_DIMS,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::HalfWidth {
                    half_width: 0.0,
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic explicitly nonadhesive interface".to_owned(),
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
            })
            .expect("adhesion claim");
        InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: material("copper-c110"),
                texture_frame: "disc-edge/profile-17".to_owned(),
            },
            SurfaceSpec {
                material: material("soda-lime-glass"),
                texture_frame: "base-track/profile-4".to_owned(),
            },
            SystemContext {
                medium: medium.to_owned(),
                third_body: None,
                environment: "air-293K-40pctRH".to_owned(),
                history: "cleaned-and-run-in-100-cycles".to_owned(),
            },
            claims,
            models,
        )
        .expect("complete ordered interface identity")
    }

    fn card(medium: &str) -> InterfaceSystemCard {
        card_with_models(medium, vec![normal_model(NORMAL_HERTZ_LAW_ID, None)])
    }

    fn trace(frame: &str) -> UniformSurfaceTrace {
        UniformSurfaceTrace::new(
            frame,
            format!("profilometer/{frame}"),
            InputAuthority::SyntheticFixture,
            1.0e-6,
            vec![0.0; 8],
            SurfaceTraceBoundary::Periodic,
        )
        .expect("synthetic trace")
    }

    fn solid_card(
        chemistry: &str,
        density: f64,
        young: f64,
        poisson: f64,
        yield_stress: f64,
    ) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, density),
            ("young_modulus", Pressure::DIMS, young),
            ("poisson_ratio", Dims::NONE, poisson),
            ("yield_stress", Pressure::DIMS, yield_stress),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                    uncertainty: UncertaintyModel::RelativeHalfWidth {
                        fraction: 0.01,
                        confidence: 0.95,
                    },
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: format!("synthetic {chemistry} {name}"),
                        license: "CC0-1.0".to_owned(),
                        artifact: None,
                    },
                })
                .expect("synthetic claim");
        }
        MaterialCard::assemble(material(chemistry), claims, Vec::new()).expect("solid card")
    }

    fn state(card: &MaterialCard, temperature_k: f64) -> IsotropicSolidStatePoint {
        resolve_isotropic_solid_state_point(
            card,
            &QueryPoint::new()
                .with("T", temperature_k)
                .expect("state point"),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("resolved solid")
    }

    fn interface_state(
        card: &InterfaceSystemCard,
        temperature_k: f64,
    ) -> ResolvedInterfaceStatePoint {
        resolve_interface_state_point(
            card,
            &QueryPoint::new()
                .with("T", temperature_k)
                .expect("state point"),
            &[ScalarPropertyRequirement::try_new(
                ADHESION_ENERGY_PROPERTY,
                ADHESION_ENERGY_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .expect("normal interface requirement")],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("resolved interface")
    }

    #[test]
    fn g0_card_binding_retains_order_materials_textures_and_identity() {
        let card = card("dry");
        let bound = bind_dry_interface_system_card(&card, InputAuthority::CallerDeclared)
            .expect("dry card binds");
        assert_eq!(bound.card_identity(), card.content_hash());
        assert_eq!(bound.surface_materials()[0].chemistry, "copper-c110");
        assert_eq!(bound.surface_materials()[1].chemistry, "soda-lime-glass");
        assert_eq!(
            bound.texture_frame_ids(),
            ["disc-edge/profile-17", "base-track/profile-4"]
        );
        assert_eq!(
            bound.interface().ordered_system_id(),
            format!("fs-matdb/interface-system-card/{}", card.content_hash())
        );
        bound
            .admit_surface_traces(
                &trace("disc-edge/profile-17"),
                &trace("base-track/profile-4"),
            )
            .expect("ordered traces match their card roles");
    }

    #[test]
    fn g0_surface_reversal_and_non_dry_medium_refuse() {
        let bound = bind_dry_interface_system_card(&card("dry"), InputAuthority::CallerDeclared)
            .expect("dry card binds");
        assert!(matches!(
            bound.admit_surface_traces(
                &trace("base-track/profile-4"),
                &trace("disc-edge/profile-17"),
            ),
            Err(InterfaceBindingError::TextureFrameMismatch {
                role: "surface_a",
                ..
            })
        ));
        assert!(matches!(
            bind_dry_interface_system_card(&card("oil-film"), InputAuthority::CallerDeclared),
            Err(InterfaceBindingError::UnsupportedMedium { .. })
        ));
    }

    #[test]
    fn g0_isotropic_binding_derives_properties_and_refuses_role_or_state_rebinding() {
        let interface_card = card("dry");
        let dry = bind_dry_interface_system_card(&interface_card, InputAuthority::CallerDeclared)
            .expect("dry interface");
        let copper = solid_card("copper-c110", 8960.0, 117.0e9, 0.34, 70.0e6);
        let glass = solid_card("soda-lime-glass", 2500.0, 72.0e9, 0.22, 1.0e9);
        let copper_room = state(&copper, 293.15);
        let glass_room = state(&glass, 293.15);
        let bound = bind_isotropic_elastic_interface(&dry, &copper_room, &glass_room)
            .expect("ordered room-temperature states bind");
        let normal = bind_normal_interface_state(&bound, &interface_state(&interface_card, 293.15))
            .expect("normal interface state binds");
        let model =
            bind_normal_contact_model(&normal, NormalContactModelSelection::SingleSupported)
                .expect("single Hertz model binds");
        assert_eq!(normal.adhesion_energy_j_per_m2(), 0.0);
        assert_eq!(model.law(), BoundNormalContactLaw::ElasticHertz);
        assert_eq!(model.characteristic_rate_m_per_s(), 1.0);
        assert_eq!(model.limits().min_temperature_k, 280.0);
        assert_eq!(normal.elastic().identity(), bound.identity());
        assert!(bound.reduced_modulus_pa() > 0.0);
        assert_eq!(bound.limiting_yield_stress_pa(), 70.0e6);
        assert_eq!(
            bound.surface_a_state_identity(),
            copper_room.resolved().identity()
        );
        assert!(matches!(
            bind_isotropic_elastic_interface(&dry, &glass_room, &copper_room),
            Err(InterfaceBindingError::MaterialStateMismatch {
                role: "surface_a",
                ..
            })
        ));
        let glass_hot = state(&glass, 300.0);
        assert!(matches!(
            bind_isotropic_elastic_interface(&dry, &copper_room, &glass_hot),
            Err(InterfaceBindingError::StatePointMismatch)
        ));

        assert!(matches!(
            bind_normal_interface_state(&bound, &interface_state(&interface_card, 300.0)),
            Err(InterfaceBindingError::StatePointMismatch)
        ));
    }

    #[test]
    fn g0_normal_model_binding_is_card_driven_and_refuses_missing_or_ambiguous_models() {
        let copper = solid_card("copper-c110", 8960.0, 117.0e9, 0.34, 70.0e6);
        let glass = solid_card("soda-lime-glass", 2500.0, 72.0e9, 0.22, 1.0e9);
        let copper_room = state(&copper, 293.15);
        let glass_room = state(&glass, 293.15);

        let hunt_card = card_with_models(
            "dry",
            vec![normal_model(
                NORMAL_HUNT_CROSSLEY_SPHERE_LAW_ID,
                Some(0.075),
            )],
        );
        let hunt_dry = bind_dry_interface_system_card(&hunt_card, InputAuthority::CallerDeclared)
            .expect("hunt interface");
        let hunt_elastic = bind_isotropic_elastic_interface(&hunt_dry, &copper_room, &glass_room)
            .expect("hunt elastic state");
        let hunt_normal =
            bind_normal_interface_state(&hunt_elastic, &interface_state(&hunt_card, 293.15))
                .expect("hunt normal state");
        let hunt =
            bind_normal_contact_model(&hunt_normal, NormalContactModelSelection::SingleSupported)
                .expect("hunt model");
        assert_eq!(
            hunt.law(),
            BoundNormalContactLaw::HuntCrossleySphere {
                dissipation_s_per_m: 0.075,
            }
        );

        let missing_card = card_with_models("dry", Vec::new());
        let missing_dry =
            bind_dry_interface_system_card(&missing_card, InputAuthority::CallerDeclared)
                .expect("model-free dry interface");
        let missing_elastic =
            bind_isotropic_elastic_interface(&missing_dry, &copper_room, &glass_room)
                .expect("model-free elastic state");
        let missing_normal =
            bind_normal_interface_state(&missing_elastic, &interface_state(&missing_card, 293.15))
                .expect("model-free normal state");
        assert_eq!(
            bind_normal_contact_model(
                &missing_normal,
                NormalContactModelSelection::SingleSupported,
            ),
            Err(InterfaceBindingError::MissingNormalContactModel)
        );

        let ambiguous_card = card_with_models(
            "dry",
            vec![
                normal_model(NORMAL_HERTZ_LAW_ID, None),
                normal_model(NORMAL_HUNT_CROSSLEY_SPHERE_LAW_ID, Some(0.075)),
            ],
        );
        let ambiguous_dry =
            bind_dry_interface_system_card(&ambiguous_card, InputAuthority::CallerDeclared)
                .expect("ambiguous dry interface");
        let ambiguous_elastic =
            bind_isotropic_elastic_interface(&ambiguous_dry, &copper_room, &glass_room)
                .expect("ambiguous elastic state");
        let ambiguous_normal = bind_normal_interface_state(
            &ambiguous_elastic,
            &interface_state(&ambiguous_card, 293.15),
        )
        .expect("ambiguous normal state");
        assert!(matches!(
            bind_normal_contact_model(
                &ambiguous_normal,
                NormalContactModelSelection::SingleSupported,
            ),
            Err(InterfaceBindingError::AmbiguousNormalContactModel { .. })
        ));
    }
}
