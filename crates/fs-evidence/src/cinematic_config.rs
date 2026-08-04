//! Versioned cinematic composition and minimally partitioned identities.
//!
//! Referenced components are already typed/versioned artifacts. Locator hints
//! remain operational metadata and never substitute for content identities.

use core::fmt;
use fs_blake3::{ContentHash, hash_domain};

/// Exact schema version accepted by [`CinematicConfig`].
pub const CINEMATIC_CONFIG_SCHEMA_VERSION: u16 = 1;
const MAX_LOCATOR_BYTES: usize = 1024;
const MAX_NAMESPACE_BYTES: usize = 96;
const MAX_ASSETS_PER_CLASS: usize = 1024;
const TRAJECTORY_DOMAIN: &str = "org.frankensim.cinematic-config.trajectory.v1";
const IMAGE_DOMAIN: &str = "org.frankensim.cinematic-config.image.v1";
const AUDIO_DOMAIN: &str = "org.frankensim.cinematic-config.audio.v1";
const MUX_DOMAIN: &str = "org.frankensim.cinematic-config.mux.v1";
const COMPOSITION_DOMAIN: &str = "org.frankensim.cinematic-config.composition.v1";
const ASSET_BYTES_DOMAIN: &str = "org.frankensim.cinematic-config.asset-bytes.v1";

/// Semantic role prevents cross-wiring equal-shaped references.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum CinematicComponentRole {
    /// Accepted Euler trajectory artifact.
    Trajectory = 1,
    /// Master frame/sample/cut timeline.
    Timeline,
    /// Camera path and optics.
    Camera,
    /// Scene charts/meshes.
    SceneGeometry,
    /// Scene instance transforms and mapping.
    InstanceMapping,
    /// Sampler, integrator, depth, and shutter configuration.
    Renderer,
    /// AOV, denoise, color, and display pipeline.
    ImagePipeline,
    /// Simulation-derived acoustic excitation.
    AudioExcitation,
    /// Modal or other declared sound synthesis model.
    SoundModel,
    /// Listener/microphone geometry.
    Microphone,
    /// Room/spatialization configuration.
    Room,
    /// Render-specific resource profile.
    RenderBudgetProfile,
    /// Audio-specific resource profile.
    AudioBudgetProfile,
}

/// Typed, versioned content reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CinematicComponentRef {
    role: CinematicComponentRole,
    identity: ContentHash,
    version: u32,
}

impl CinematicComponentRef {
    /// Bind a semantic role to a nonzero versioned content identity.
    pub fn try_new(
        role: CinematicComponentRole,
        identity: ContentHash,
        version: u32,
    ) -> Result<Self, CinematicConfigError> {
        check_hash(identity)?;
        if version == 0 {
            return Err(CinematicConfigError::InvalidComponentVersion(role));
        }
        Ok(Self {
            role,
            identity,
            version,
        })
    }

    #[must_use]
    /// Semantic component role.
    pub const fn role(self) -> CinematicComponentRole {
        self.role
    }
    #[must_use]
    /// Content identity.
    pub const fn identity(self) -> ContentHash {
        self.identity
    }
    #[must_use]
    /// Component schema/configuration version.
    pub const fn version(self) -> u32 {
        self.version
    }
}

/// Interpretation bound together with external bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum CinematicAssetInterpretation {
    /// Wavelength-dependent reflectance data.
    SpectralReflectance = 1,
    /// Wavelength-dependent emission data.
    SpectralEmission,
    /// Linear scene-referred texture samples.
    LinearSceneTexture,
    /// Display-referred encoded texture samples.
    DisplayReferredTexture,
    /// Geometry whose coordinates are metres.
    GeometryMeters,
    /// SI acoustic impulse-response samples.
    AcousticImpulseResponseSi,
}

/// Asset bytes, interpretation, version, and relocatable lookup hint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicAssetBinding {
    content_identity: ContentHash,
    interpretation: CinematicAssetInterpretation,
    version: u32,
    locator_hint: String,
}

impl CinematicAssetBinding {
    fn from_identity(
        content_identity: ContentHash,
        interpretation: CinematicAssetInterpretation,
        version: u32,
        locator_hint: String,
    ) -> Result<Self, CinematicConfigError> {
        check_hash(content_identity)?;
        if version == 0 {
            return Err(CinematicConfigError::InvalidAssetVersion);
        }
        validate_locator(&locator_hint)?;
        Ok(Self {
            content_identity,
            interpretation,
            version,
            locator_hint,
        })
    }

    /// Bind the supplied bytes directly, avoiding path-based identity.
    pub fn from_bytes(
        bytes: &[u8],
        interpretation: CinematicAssetInterpretation,
        version: u32,
        locator_hint: String,
    ) -> Result<Self, CinematicConfigError> {
        Self::from_identity(
            hash_domain(ASSET_BYTES_DOMAIN, bytes),
            interpretation,
            version,
            locator_hint,
        )
    }

    /// Relocate an admitted asset without changing its content identity.
    pub fn with_locator_hint(&self, locator_hint: String) -> Result<Self, CinematicConfigError> {
        Self::from_identity(
            self.content_identity,
            self.interpretation,
            self.version,
            locator_hint,
        )
    }

    /// Refuse stale or substituted bytes at asset load time.
    pub fn verify_bytes(&self, bytes: &[u8]) -> Result<(), CinematicConfigError> {
        if hash_domain(ASSET_BYTES_DOMAIN, bytes) == self.content_identity {
            Ok(())
        } else {
            Err(CinematicConfigError::AssetContentMismatch)
        }
    }

    #[must_use]
    /// Identity of the bound bytes.
    pub const fn content_identity(&self) -> ContentHash {
        self.content_identity
    }
    #[must_use]
    /// Declared interpretation of those bytes.
    pub const fn interpretation(&self) -> CinematicAssetInterpretation {
        self.interpretation
    }
    #[must_use]
    /// Non-authoritative lookup hint.
    pub fn locator_hint(&self) -> &str {
        &self.locator_hint
    }

    fn key(&self) -> (ContentHash, CinematicAssetInterpretation, u32) {
        (self.content_identity, self.interpretation, self.version)
    }
}

/// Composition-level unit basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicConfigUnits {
    /// SI mechanics plus radians and 48 kHz sample-frame clocks.
    SiMetersKilogramsSecondsRadians = 1,
}

/// Explicit admitted capability bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CinematicCapabilities(u32);

impl CinematicCapabilities {
    /// Image rendering capability.
    pub const RENDER: u32 = 1 << 0;
    /// Sound synthesis capability.
    pub const AUDIO: u32 = 1 << 1;
    /// Quarantined out-of-process mux capability.
    pub const QUARANTINED_MUX: u32 = 1 << 2;
    const KNOWN: u32 = Self::RENDER | Self::AUDIO | Self::QUARANTINED_MUX;

    /// Admit a nonempty subset of the known capability bits.
    pub fn try_new(bits: u32) -> Result<Self, CinematicConfigError> {
        if bits == 0 || bits & !Self::KNOWN != 0 {
            return Err(CinematicConfigError::InvalidCapabilities(bits));
        }
        Ok(Self(bits))
    }
    #[must_use]
    /// Canonical capability bitset.
    pub const fn bits(self) -> u32 {
        self.0
    }
    const fn contains(self, bits: u32) -> bool {
        self.0 & bits == bits
    }
}

/// Supported non-authoritative delivery codecs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicMuxCodec {
    /// AV1 video and Opus audio in Matroska.
    Av1OpusMatroska = 1,
    /// H.265 video and PCM audio in QuickTime.
    H265PcmQuickTime,
}

/// Optional quarantined out-of-process media assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CinematicMuxRequest {
    /// Do not request a muxed derivative.
    None,
    /// Request a separately distributed quarantined adapter.
    QuarantinedAdapter {
        /// Content identity of the adapter binary/receipt.
        adapter_identity: ContentHash,
        /// Nonzero adapter version.
        adapter_version: u32,
        /// Requested delivery codec.
        codec: CinematicMuxCodec,
    },
}

/// Stable logical artifact namespace plus relocatable physical hint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicArtifactRoot {
    logical_namespace: String,
    locator_hint: String,
}

impl CinematicArtifactRoot {
    /// Validate a logical namespace and its current lookup hint.
    pub fn try_new(
        logical_namespace: String,
        locator_hint: String,
    ) -> Result<Self, CinematicConfigError> {
        validate_namespace(&logical_namespace)?;
        validate_locator(&locator_hint)?;
        Ok(Self {
            logical_namespace,
            locator_hint,
        })
    }
    #[must_use]
    /// Identity-bearing logical namespace.
    pub fn logical_namespace(&self) -> &str {
        &self.logical_namespace
    }
    #[must_use]
    /// Non-authoritative current storage hint.
    pub fn locator_hint(&self) -> &str {
        &self.locator_hint
    }
}

/// Complete caller input; optional fields make missing Five Explicits refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicConfigInput {
    /// Exact configuration schema version.
    pub schema_version: u16,
    /// Mandatory units explicit.
    pub units: Option<CinematicConfigUnits>,
    /// Mandatory seed explicit; zero remains a valid declared seed.
    pub seed: Option<u64>,
    /// Mandatory capability explicit.
    pub capabilities: Option<CinematicCapabilities>,
    /// Mandatory render resource budget identity.
    pub render_budget_profile: Option<CinematicComponentRef>,
    /// Mandatory audio resource budget identity.
    pub audio_budget_profile: Option<CinematicComponentRef>,
    /// Source trajectory.
    pub trajectory: CinematicComponentRef,
    /// Shared master timeline.
    pub timeline: CinematicComponentRef,
    /// Image-only camera configuration.
    pub camera: CinematicComponentRef,
    /// Renderable geometry.
    pub scene_geometry: CinematicComponentRef,
    /// Geometry instance mapping.
    pub instance_mapping: CinematicComponentRef,
    /// Render integrator configuration.
    pub renderer: CinematicComponentRef,
    /// Image finishing configuration.
    pub image_pipeline: CinematicComponentRef,
    /// Audio excitation derived from the trajectory.
    pub audio_excitation: CinematicComponentRef,
    /// Sound synthesis configuration.
    pub sound_model: CinematicComponentRef,
    /// Microphone/listener configuration.
    pub microphone: CinematicComponentRef,
    /// Room/spatialization configuration.
    pub room: CinematicComponentRef,
    /// Order-insensitive material assets.
    pub material_assets: Vec<CinematicAssetBinding>,
    /// Order-insensitive lighting assets.
    pub light_assets: Vec<CinematicAssetBinding>,
    /// Explicit environment; there is no hidden default.
    pub environment_asset: CinematicAssetBinding,
    /// Logical output namespace and relocatable hint.
    pub artifact_root: CinematicArtifactRoot,
    /// Optional external delivery derivative.
    pub mux_request: CinematicMuxRequest,
}

/// Validated configuration and partition identities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicConfig {
    input: CinematicConfigInput,
    trajectory_identity: ContentHash,
    image_identity: ContentHash,
    audio_identity: ContentHash,
    mux_identity: ContentHash,
    composition_identity: ContentHash,
}

impl CinematicConfig {
    /// Validate/canonicalize all fields and derive partition identities.
    pub fn try_new(mut input: CinematicConfigInput) -> Result<Self, CinematicConfigError> {
        if input.schema_version != CINEMATIC_CONFIG_SCHEMA_VERSION {
            return Err(CinematicConfigError::UnsupportedSchemaVersion(
                input.schema_version,
            ));
        }
        let units = input.units.ok_or(CinematicConfigError::MissingUnits)?;
        let seed = input.seed.ok_or(CinematicConfigError::MissingSeed)?;
        let capabilities = input
            .capabilities
            .ok_or(CinematicConfigError::MissingCapabilities)?;
        let render_budget = input
            .render_budget_profile
            .ok_or(CinematicConfigError::MissingBudget)?;
        let audio_budget = input
            .audio_budget_profile
            .ok_or(CinematicConfigError::MissingBudget)?;
        if !capabilities.contains(CinematicCapabilities::RENDER | CinematicCapabilities::AUDIO) {
            return Err(CinematicConfigError::MissingRequiredCapability);
        }
        if matches!(
            input.mux_request,
            CinematicMuxRequest::QuarantinedAdapter { .. }
        ) && !capabilities.contains(CinematicCapabilities::QUARANTINED_MUX)
        {
            return Err(CinematicConfigError::MissingMuxCapability);
        }
        validate_mux(input.mux_request)?;
        if input.material_assets.is_empty() {
            return Err(CinematicConfigError::MissingMaterialAssets);
        }
        if input.light_assets.is_empty() {
            return Err(CinematicConfigError::MissingLightAssets);
        }
        if input.material_assets.len() > MAX_ASSETS_PER_CLASS
            || input.light_assets.len() > MAX_ASSETS_PER_CLASS
        {
            return Err(CinematicConfigError::TooManyAssets);
        }
        for (component, expected) in [
            (render_budget, CinematicComponentRole::RenderBudgetProfile),
            (audio_budget, CinematicComponentRole::AudioBudgetProfile),
            (input.trajectory, CinematicComponentRole::Trajectory),
            (input.timeline, CinematicComponentRole::Timeline),
            (input.camera, CinematicComponentRole::Camera),
            (input.scene_geometry, CinematicComponentRole::SceneGeometry),
            (
                input.instance_mapping,
                CinematicComponentRole::InstanceMapping,
            ),
            (input.renderer, CinematicComponentRole::Renderer),
            (input.image_pipeline, CinematicComponentRole::ImagePipeline),
            (
                input.audio_excitation,
                CinematicComponentRole::AudioExcitation,
            ),
            (input.sound_model, CinematicComponentRole::SoundModel),
            (input.microphone, CinematicComponentRole::Microphone),
            (input.room, CinematicComponentRole::Room),
        ] {
            if component.role != expected {
                return Err(CinematicConfigError::ComponentRoleMismatch {
                    expected,
                    got: component.role,
                });
            }
        }
        canonicalize_assets(&mut input.material_assets)?;
        canonicalize_assets(&mut input.light_assets)?;

        let trajectory_identity = hash_components(TRAJECTORY_DOMAIN, &[input.trajectory]);
        let image_identity = image_hash(&input, render_budget, units, seed);
        let audio_identity = audio_hash(&input, audio_budget, units, seed);
        let mux_identity = mux_hash(&input, image_identity, audio_identity);
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&input.schema_version.to_le_bytes());
        for identity in [
            trajectory_identity,
            image_identity,
            audio_identity,
            mux_identity,
        ] {
            push_hash(&mut bytes, identity);
        }
        push_string(&mut bytes, &input.artifact_root.logical_namespace);
        let composition_identity = hash_domain(COMPOSITION_DOMAIN, &bytes);
        Ok(Self {
            input,
            trajectory_identity,
            image_identity,
            audio_identity,
            mux_identity,
            composition_identity,
        })
    }

    #[must_use]
    /// Identity of source trajectory work only.
    pub const fn trajectory_identity(&self) -> ContentHash {
        self.trajectory_identity
    }
    #[must_use]
    /// Identity of all inputs needed to produce image masters.
    pub const fn image_identity(&self) -> ContentHash {
        self.image_identity
    }
    #[must_use]
    /// Identity of all inputs needed to produce audio masters.
    pub const fn audio_identity(&self) -> ContentHash {
        self.audio_identity
    }
    #[must_use]
    /// Identity of image/audio/timeline and mux request.
    pub const fn mux_identity(&self) -> ContentHash {
        self.mux_identity
    }
    #[must_use]
    /// Identity of the complete composition and logical output namespace.
    pub const fn composition_identity(&self) -> ContentHash {
        self.composition_identity
    }

    /// Canonical compact composition preimage. Component detail remains in
    /// the referenced content-addressed artifacts; locator hints are omitted.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.input.schema_version.to_le_bytes());
        for identity in [
            self.trajectory_identity,
            self.image_identity,
            self.audio_identity,
            self.mux_identity,
        ] {
            push_hash(&mut bytes, identity);
        }
        push_string(&mut bytes, &self.input.artifact_root.logical_namespace);
        bytes
    }
    #[must_use]
    /// Canonicalized admitted input; locator hints remain available here.
    pub const fn input(&self) -> &CinematicConfigInput {
        &self.input
    }
}

fn image_hash(
    input: &CinematicConfigInput,
    budget: CinematicComponentRef,
    units: CinematicConfigUnits,
    seed: u64,
) -> ContentHash {
    let mut bytes = Vec::new();
    for component in [
        input.trajectory,
        input.timeline,
        input.camera,
        input.scene_geometry,
        input.instance_mapping,
        input.renderer,
        input.image_pipeline,
        budget,
    ] {
        push_component(&mut bytes, component);
    }
    push_explicits(&mut bytes, units, seed, CinematicCapabilities::RENDER);
    push_assets(&mut bytes, &input.material_assets);
    push_assets(&mut bytes, &input.light_assets);
    push_asset(&mut bytes, &input.environment_asset);
    hash_domain(IMAGE_DOMAIN, &bytes)
}

fn audio_hash(
    input: &CinematicConfigInput,
    budget: CinematicComponentRef,
    units: CinematicConfigUnits,
    seed: u64,
) -> ContentHash {
    let mut bytes = Vec::new();
    for component in [
        input.trajectory,
        input.timeline,
        input.audio_excitation,
        input.sound_model,
        input.microphone,
        input.room,
        budget,
    ] {
        push_component(&mut bytes, component);
    }
    push_explicits(&mut bytes, units, seed, CinematicCapabilities::AUDIO);
    hash_domain(AUDIO_DOMAIN, &bytes)
}

fn mux_hash(input: &CinematicConfigInput, image: ContentHash, audio: ContentHash) -> ContentHash {
    let mut bytes = Vec::new();
    push_hash(&mut bytes, image);
    push_hash(&mut bytes, audio);
    push_component(&mut bytes, input.timeline);
    match input.mux_request {
        CinematicMuxRequest::None => bytes.push(0),
        CinematicMuxRequest::QuarantinedAdapter {
            adapter_identity,
            adapter_version,
            codec,
        } => {
            bytes.push(1);
            push_hash(&mut bytes, adapter_identity);
            bytes.extend_from_slice(&adapter_version.to_le_bytes());
            bytes.push(codec as u8);
        }
    }
    hash_domain(MUX_DOMAIN, &bytes)
}

fn hash_components(domain: &str, components: &[CinematicComponentRef]) -> ContentHash {
    let mut bytes = Vec::new();
    for component in components {
        push_component(&mut bytes, *component);
    }
    hash_domain(domain, &bytes)
}

fn push_component(bytes: &mut Vec<u8>, value: CinematicComponentRef) {
    bytes.push(value.role as u8);
    push_hash(bytes, value.identity);
    bytes.extend_from_slice(&value.version.to_le_bytes());
}

fn push_explicits(
    bytes: &mut Vec<u8>,
    units: CinematicConfigUnits,
    seed: u64,
    relevant_capability: u32,
) {
    bytes.push(units as u8);
    bytes.extend_from_slice(&seed.to_le_bytes());
    bytes.extend_from_slice(&relevant_capability.to_le_bytes());
}

fn push_assets(bytes: &mut Vec<u8>, assets: &[CinematicAssetBinding]) {
    bytes.extend_from_slice(
        &u16::try_from(assets.len())
            .unwrap_or(u16::MAX)
            .to_le_bytes(),
    );
    for asset in assets {
        push_asset(bytes, asset);
    }
}

fn push_asset(bytes: &mut Vec<u8>, asset: &CinematicAssetBinding) {
    push_hash(bytes, asset.content_identity);
    bytes.push(asset.interpretation as u8);
    bytes.extend_from_slice(&asset.version.to_le_bytes());
}

fn push_hash(bytes: &mut Vec<u8>, value: ContentHash) {
    bytes.extend_from_slice(value.as_bytes());
}
fn push_string(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&u32::try_from(value.len()).unwrap_or(u32::MAX).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

fn canonicalize_assets(assets: &mut [CinematicAssetBinding]) -> Result<(), CinematicConfigError> {
    assets.sort_unstable_by_key(CinematicAssetBinding::key);
    if assets.windows(2).any(|pair| pair[0].key() == pair[1].key()) {
        return Err(CinematicConfigError::DuplicateAsset);
    }
    Ok(())
}

fn validate_mux(request: CinematicMuxRequest) -> Result<(), CinematicConfigError> {
    if let CinematicMuxRequest::QuarantinedAdapter {
        adapter_identity,
        adapter_version,
        ..
    } = request
    {
        check_hash(adapter_identity)?;
        if adapter_version == 0 {
            return Err(CinematicConfigError::InvalidMuxAdapterVersion);
        }
    }
    Ok(())
}

fn validate_locator(value: &str) -> Result<(), CinematicConfigError> {
    if value.is_empty() || value.len() > MAX_LOCATOR_BYTES || value.chars().any(char::is_control) {
        Err(CinematicConfigError::InvalidLocator)
    } else {
        Ok(())
    }
}

fn validate_namespace(value: &str) -> Result<(), CinematicConfigError> {
    let bytes = value.as_bytes();
    if bytes.is_empty()
        || bytes.len() > MAX_NAMESPACE_BYTES
        || !bytes[0].is_ascii_lowercase()
        || !bytes.iter().all(|byte| {
            byte.is_ascii_lowercase()
                || byte.is_ascii_digit()
                || matches!(byte, b'-' | b'_' | b'.' | b'/')
        })
    {
        Err(CinematicConfigError::InvalidArtifactNamespace)
    } else {
        Ok(())
    }
}

fn check_hash(value: ContentHash) -> Result<(), CinematicConfigError> {
    if value.as_bytes().iter().all(|byte| *byte == 0) {
        Err(CinematicConfigError::MissingContentIdentity)
    } else {
        Ok(())
    }
}

/// Stable configuration-admission refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicConfigError {
    /// Schema version is not supported.
    UnsupportedSchemaVersion(u16),
    /// Units explicit is absent.
    MissingUnits,
    /// Seed explicit is absent.
    MissingSeed,
    /// Capability explicit is absent.
    MissingCapabilities,
    /// Render or audio budget explicit is absent.
    MissingBudget,
    /// Render and audio capabilities are not both admitted.
    MissingRequiredCapability,
    /// A mux was requested without quarantined-adapter capability.
    MissingMuxCapability,
    /// Capability bits are empty or unknown.
    InvalidCapabilities(u32),
    /// A required identity was all zero.
    MissingContentIdentity,
    /// A component version was zero.
    InvalidComponentVersion(CinematicComponentRole),
    /// A reference was wired into the wrong semantic port.
    ComponentRoleMismatch {
        /// Port role.
        expected: CinematicComponentRole,
        /// Supplied role.
        got: CinematicComponentRole,
    },
    /// An asset version was zero.
    InvalidAssetVersion,
    /// A locator was empty, excessive, or contained control text.
    InvalidLocator,
    /// Logical artifact namespace was malformed.
    InvalidArtifactNamespace,
    /// Identical asset bindings appeared twice.
    DuplicateAsset,
    /// The composition omitted every material asset.
    MissingMaterialAssets,
    /// The composition omitted every light asset.
    MissingLightAssets,
    /// One asset class exceeded its bounded canonical count.
    TooManyAssets,
    /// Loaded bytes disagreed with their admitted content identity.
    AssetContentMismatch,
    /// Mux adapter version was zero.
    InvalidMuxAdapterVersion,
}

impl CinematicConfigError {
    /// Stable machine diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::UnsupportedSchemaVersion(_) => "cinematic-config-unsupported-schema",
            Self::MissingUnits => "cinematic-config-missing-units",
            Self::MissingSeed => "cinematic-config-missing-seed",
            Self::MissingCapabilities => "cinematic-config-missing-capabilities",
            Self::MissingBudget => "cinematic-config-missing-budget",
            Self::MissingRequiredCapability => "cinematic-config-missing-render-audio-capability",
            Self::MissingMuxCapability => "cinematic-config-missing-mux-capability",
            Self::InvalidCapabilities(_) => "cinematic-config-invalid-capabilities",
            Self::MissingContentIdentity => "cinematic-config-missing-content-identity",
            Self::InvalidComponentVersion(_) => "cinematic-config-invalid-component-version",
            Self::ComponentRoleMismatch { .. } => "cinematic-config-component-role-mismatch",
            Self::InvalidAssetVersion => "cinematic-config-invalid-asset-version",
            Self::InvalidLocator => "cinematic-config-invalid-locator",
            Self::InvalidArtifactNamespace => "cinematic-config-invalid-artifact-namespace",
            Self::DuplicateAsset => "cinematic-config-duplicate-asset",
            Self::MissingMaterialAssets => "cinematic-config-missing-material-assets",
            Self::MissingLightAssets => "cinematic-config-missing-light-assets",
            Self::TooManyAssets => "cinematic-config-too-many-assets",
            Self::AssetContentMismatch => "cinematic-config-asset-content-mismatch",
            Self::InvalidMuxAdapterVersion => "cinematic-config-invalid-mux-version",
        }
    }
}

impl fmt::Display for CinematicConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {self:?}", self.code())
    }
}
impl std::error::Error for CinematicConfigError {}
