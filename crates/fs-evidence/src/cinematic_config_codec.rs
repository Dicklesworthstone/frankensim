//! Strict, bounded transport for complete cinematic composition documents.
//!
//! [`crate::cinematic_config::CinematicConfig::canonical_bytes`] is an
//! identity preimage, not a reconstructible user file. This module owns the
//! complete v1 document instead: every component reference, named quality
//! profile, interpreted external asset, output namespace, and optional mux
//! request is retained. Asset identities are still not trusted from text;
//! admission resolves the declared locator, hashes the returned bytes in
//! cancellation-bounded tiles using the same asset-identity domain, and
//! compares the expected identity.

use core::fmt;
use std::collections::BTreeMap;

use fs_blake3::{ContentHash, DomainHasher};

use crate::{
    cinematic_budget::{
        CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION, CinematicBudgetError, CinematicQualityProfile,
        CinematicQualityTier,
    },
    cinematic_config::{
        ASSET_BYTES_DOMAIN, CINEMATIC_CONFIG_SCHEMA_VERSION, CinematicArtifactRoot,
        CinematicAssetBinding, CinematicAssetInterpretation, CinematicCapabilities,
        CinematicComponentRef, CinematicComponentRole, CinematicConfig, CinematicConfigError,
        CinematicConfigInput, CinematicConfigUnits, CinematicMuxCodec, CinematicMuxRequest,
    },
};

/// Exact reconstructible document schema accepted by v1 readers.
pub const CINEMATIC_CONFIG_DOCUMENT_SCHEMA: &str = "frankensim.cinematic-config-document.v1";
/// Maximum complete UTF-8 document bytes.
pub const MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES: usize = 1024 * 1024;
/// Maximum physical lines, including comments and blanks.
pub const MAX_CINEMATIC_CONFIG_DOCUMENT_LINES: usize = 4_096;
/// Maximum declarations in either order-insensitive asset class.
pub const MAX_CINEMATIC_CONFIG_DOCUMENT_ASSETS_PER_CLASS: usize = 1_024;
/// Default maximum bytes returned for one resolved external asset.
pub const MAX_CINEMATIC_RESOLVED_ASSET_BYTES: usize = 256 * 1024 * 1024;
/// Default maximum bytes returned across all resolved external assets.
pub const MAX_CINEMATIC_RESOLVED_ASSET_TOTAL_BYTES: usize = 1024 * 1024 * 1024;
/// Maximum bytes hashed between caller cancellation checkpoints.
pub const CINEMATIC_ASSET_HASH_TILE_BYTES: usize = 64 * 1_024;

const MAX_LOCATOR_BYTES: usize = 1_024;

/// Location-independent class used in resolver diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CinematicDocumentAssetClass {
    /// Material/reflectance asset.
    Material,
    /// Light/emission asset.
    Light,
    /// Unique environment asset.
    Environment,
}

impl CinematicDocumentAssetClass {
    /// Stable field-path stem.
    #[must_use]
    pub const fn field_name(self) -> &'static str {
        match self {
            Self::Material => "material_asset",
            Self::Light => "light_asset",
            Self::Environment => "environment_asset",
        }
    }
}

/// Bounded external-input failure returned by a caller-supplied resolver.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CinematicAssetAccessError {
    /// Cancellation was requested while loading the asset.
    Cancelled,
    /// The locator could not be opened or read.
    Unavailable,
    /// The asset exceeded the resolver's admitted byte ceiling.
    TooLarge,
    /// The resolver could not reserve bounded storage.
    Capacity,
}

/// Explicit post-resolution byte envelope enforced before asset hashing.
///
/// The resolver owns its I/O allocation strategy. This budget prevents bytes
/// returned by that caller from entering configuration admission or hashing
/// outside a declared per-asset and aggregate bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct CinematicAssetAdmissionBudget {
    /// Maximum returned bytes for one declaration.
    pub max_asset_bytes: usize,
    /// Maximum returned bytes across material, light, and environment assets.
    pub max_total_asset_bytes: usize,
}

impl CinematicAssetAdmissionBudget {
    /// Default CLI-scale envelope.
    pub const DEFAULT: Self = Self {
        max_asset_bytes: MAX_CINEMATIC_RESOLVED_ASSET_BYTES,
        max_total_asset_bytes: MAX_CINEMATIC_RESOLVED_ASSET_TOTAL_BYTES,
    };
}

impl Default for CinematicAssetAdmissionBudget {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Expected asset identity plus its interpretation and relocatable hint.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicAssetDeclaration {
    expected_identity: ContentHash,
    interpretation: CinematicAssetInterpretation,
    version: u32,
    locator_hint: String,
}

impl CinematicAssetDeclaration {
    /// Expected domain-separated identity of the resolved bytes.
    #[must_use]
    pub const fn expected_identity(&self) -> ContentHash {
        self.expected_identity
    }

    /// Declared interpretation of those bytes.
    #[must_use]
    pub const fn interpretation(&self) -> CinematicAssetInterpretation {
        self.interpretation
    }

    /// Nonzero interpretation/schema version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// Relocatable lookup hint; it is not part of composition identity.
    #[must_use]
    pub fn locator_hint(&self) -> &str {
        &self.locator_hint
    }

    fn key(&self) -> (ContentHash, CinematicAssetInterpretation, u32) {
        (self.expected_identity, self.interpretation, self.version)
    }
}

/// Complete decoded document whose assets have not yet been resolved.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicConfigDocument {
    quality_profile: CinematicQualityTier,
    units: CinematicConfigUnits,
    seed: u64,
    capabilities: CinematicCapabilities,
    render_budget_profile: CinematicComponentRef,
    audio_budget_profile: CinematicComponentRef,
    trajectory: CinematicComponentRef,
    timeline: CinematicComponentRef,
    camera: CinematicComponentRef,
    scene_geometry: CinematicComponentRef,
    instance_mapping: CinematicComponentRef,
    renderer: CinematicComponentRef,
    image_pipeline: CinematicComponentRef,
    audio_excitation: CinematicComponentRef,
    sound_model: CinematicComponentRef,
    microphone: CinematicComponentRef,
    room: CinematicComponentRef,
    material_assets: Vec<CinematicAssetDeclaration>,
    light_assets: Vec<CinematicAssetDeclaration>,
    environment_asset: CinematicAssetDeclaration,
    artifact_root: CinematicArtifactRoot,
    mux_request: CinematicMuxRequest,
}

impl CinematicConfigDocument {
    /// Decode one complete UTF-8 document under hard byte and line ceilings.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, CinematicConfigDocumentError> {
        if bytes.len() > MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES {
            return Err(CinematicConfigDocumentError::DocumentTooLarge {
                bytes: bytes.len(),
                maximum: MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES,
            });
        }
        let source =
            core::str::from_utf8(bytes).map_err(|_| CinematicConfigDocumentError::InvalidUtf8)?;
        Self::from_str(source)
    }

    /// Decode one complete document. Comments begin with `#`; all other
    /// nonblank lines are strict `key=value` records.
    #[allow(clippy::too_many_lines)]
    pub fn from_str(source: &str) -> Result<Self, CinematicConfigDocumentError> {
        if source.len() > MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES {
            return Err(CinematicConfigDocumentError::DocumentTooLarge {
                bytes: source.len(),
                maximum: MAX_CINEMATIC_CONFIG_DOCUMENT_BYTES,
            });
        }
        let mut unique = BTreeMap::<String, String>::new();
        let mut materials = Vec::new();
        let mut lights = Vec::new();
        let mut line_count = 0usize;
        for (zero_line, raw) in source.lines().enumerate() {
            line_count = line_count.saturating_add(1);
            if line_count > MAX_CINEMATIC_CONFIG_DOCUMENT_LINES {
                return Err(CinematicConfigDocumentError::TooManyLines);
            }
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let line_number = zero_line + 1;
            let (key, value) = line
                .split_once('=')
                .ok_or(CinematicConfigDocumentError::MalformedLine { line: line_number })?;
            let key = key.trim();
            let value = value.trim();
            if key.is_empty() || value.is_empty() {
                return Err(CinematicConfigDocumentError::MalformedLine { line: line_number });
            }
            match key {
                "material_asset" => {
                    push_asset_declaration(&mut materials, value, "material_asset")?;
                }
                "light_asset" => {
                    push_asset_declaration(&mut lights, value, "light_asset")?;
                }
                "schema"
                | "quality_profile"
                | "units"
                | "seed"
                | "capabilities"
                | "render_budget_profile"
                | "audio_budget_profile"
                | "trajectory"
                | "timeline"
                | "camera"
                | "scene_geometry"
                | "instance_mapping"
                | "renderer"
                | "image_pipeline"
                | "audio_excitation"
                | "sound_model"
                | "microphone"
                | "room"
                | "environment_asset"
                | "artifact_namespace"
                | "artifact_root"
                | "mux" => {
                    if unique.insert(key.to_owned(), value.to_owned()).is_some() {
                        return Err(CinematicConfigDocumentError::DuplicateField {
                            field: key.to_owned(),
                        });
                    }
                }
                _ => {
                    return Err(CinematicConfigDocumentError::UnknownField {
                        field: key.to_owned(),
                    });
                }
            }
        }
        if source.as_bytes().last().is_some_and(|byte| *byte == b'\r') {
            return Err(CinematicConfigDocumentError::InvalidLineEnding);
        }

        let schema = take_required(&mut unique, "schema")?;
        if schema != CINEMATIC_CONFIG_DOCUMENT_SCHEMA {
            return Err(CinematicConfigDocumentError::UnsupportedSchema);
        }
        let quality_profile =
            parse_quality_profile(&take_required(&mut unique, "quality_profile")?)?;
        let units = parse_units(&take_required(&mut unique, "units")?)?;
        let seed = parse_u64(&take_required(&mut unique, "seed")?, "seed")?;
        let capabilities = parse_capabilities(&take_required(&mut unique, "capabilities")?)?;
        let render_budget_profile = parse_component(
            &take_required(&mut unique, "render_budget_profile")?,
            "render_budget_profile",
            CinematicComponentRole::RenderBudgetProfile,
        )?;
        let audio_budget_profile = parse_component(
            &take_required(&mut unique, "audio_budget_profile")?,
            "audio_budget_profile",
            CinematicComponentRole::AudioBudgetProfile,
        )?;
        let trajectory = parse_component(
            &take_required(&mut unique, "trajectory")?,
            "trajectory",
            CinematicComponentRole::Trajectory,
        )?;
        let timeline = parse_component(
            &take_required(&mut unique, "timeline")?,
            "timeline",
            CinematicComponentRole::Timeline,
        )?;
        let camera = parse_component(
            &take_required(&mut unique, "camera")?,
            "camera",
            CinematicComponentRole::Camera,
        )?;
        let scene_geometry = parse_component(
            &take_required(&mut unique, "scene_geometry")?,
            "scene_geometry",
            CinematicComponentRole::SceneGeometry,
        )?;
        let instance_mapping = parse_component(
            &take_required(&mut unique, "instance_mapping")?,
            "instance_mapping",
            CinematicComponentRole::InstanceMapping,
        )?;
        let renderer = parse_component(
            &take_required(&mut unique, "renderer")?,
            "renderer",
            CinematicComponentRole::Renderer,
        )?;
        let image_pipeline = parse_component(
            &take_required(&mut unique, "image_pipeline")?,
            "image_pipeline",
            CinematicComponentRole::ImagePipeline,
        )?;
        let audio_excitation = parse_component(
            &take_required(&mut unique, "audio_excitation")?,
            "audio_excitation",
            CinematicComponentRole::AudioExcitation,
        )?;
        let sound_model = parse_component(
            &take_required(&mut unique, "sound_model")?,
            "sound_model",
            CinematicComponentRole::SoundModel,
        )?;
        let microphone = parse_component(
            &take_required(&mut unique, "microphone")?,
            "microphone",
            CinematicComponentRole::Microphone,
        )?;
        let room = parse_component(
            &take_required(&mut unique, "room")?,
            "room",
            CinematicComponentRole::Room,
        )?;
        let environment_asset = parse_asset_declaration(
            &take_required(&mut unique, "environment_asset")?,
            "environment_asset",
        )?;
        let artifact_namespace = take_required(&mut unique, "artifact_namespace")?;
        let artifact_locator = take_required(&mut unique, "artifact_root")?;
        let artifact_root = CinematicArtifactRoot::try_new(artifact_namespace, artifact_locator)
            .map_err(|error| match error {
                CinematicConfigError::InvalidArtifactNamespace => {
                    CinematicConfigDocumentError::InvalidField {
                        field: "artifact_namespace",
                    }
                }
                CinematicConfigError::InvalidLocator => {
                    CinematicConfigDocumentError::InvalidField {
                        field: "artifact_root",
                    }
                }
                other => CinematicConfigDocumentError::Config(other),
            })?;
        let mux_request = parse_mux(&take_required(&mut unique, "mux")?)?;
        if matches!(mux_request, CinematicMuxRequest::QuarantinedAdapter { .. })
            && capabilities.bits() & CinematicCapabilities::QUARANTINED_MUX == 0
        {
            return Err(CinematicConfigDocumentError::InvalidField {
                field: "capabilities",
            });
        }
        debug_assert!(unique.is_empty());

        canonicalize_declarations(&mut materials, CinematicDocumentAssetClass::Material)?;
        canonicalize_declarations(&mut lights, CinematicDocumentAssetClass::Light)?;
        if materials.is_empty() {
            return Err(CinematicConfigDocumentError::MissingField {
                field: "material_asset".to_owned(),
            });
        }
        if lights.is_empty() {
            return Err(CinematicConfigDocumentError::MissingField {
                field: "light_asset".to_owned(),
            });
        }

        let document = Self {
            quality_profile,
            units,
            seed,
            capabilities,
            render_budget_profile,
            audio_budget_profile,
            trajectory,
            timeline,
            camera,
            scene_geometry,
            instance_mapping,
            renderer,
            image_pipeline,
            audio_excitation,
            sound_model,
            microphone,
            room,
            material_assets: materials,
            light_assets: lights,
            environment_asset,
            artifact_root,
            mux_request,
        };
        document.quality_profile()?;
        Ok(document)
    }

    /// Frozen named quality profile selected by the complete document.
    pub fn quality_profile(&self) -> Result<CinematicQualityProfile, CinematicConfigDocumentError> {
        let profile = CinematicQualityProfile::canonical(self.quality_profile)
            .map_err(CinematicConfigDocumentError::Budget)?;
        let expected = profile.identity();
        let expected_version = u32::from(CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION);
        for (field, reference) in [
            ("render_budget_profile", self.render_budget_profile),
            ("audio_budget_profile", self.audio_budget_profile),
        ] {
            if reference.identity() != expected || reference.version() != expected_version {
                return Err(CinematicConfigDocumentError::BudgetProfileMismatch { field });
            }
        }
        Ok(profile)
    }

    /// Exact caller-declared (not authenticated) capability bits.
    #[must_use]
    pub const fn capabilities(&self) -> CinematicCapabilities {
        self.capabilities
    }

    /// Expected trajectory artifact reference.
    #[must_use]
    pub const fn trajectory(&self) -> CinematicComponentRef {
        self.trajectory
    }

    /// Optional quarantined mux request.
    #[must_use]
    pub const fn mux_request(&self) -> CinematicMuxRequest {
        self.mux_request
    }

    /// Identity-bearing logical output namespace. The physical locator is not
    /// exposed here so diagnostic callers need not accidentally leak it.
    #[must_use]
    pub fn artifact_namespace(&self) -> &str {
        self.artifact_root.logical_namespace()
    }

    /// Current physical output hint. Callers must apply their own filesystem
    /// policy before writing; this string carries no authority.
    #[must_use]
    pub fn artifact_locator_hint(&self) -> &str {
        self.artifact_root.locator_hint()
    }

    /// Resolve every external asset, verify expected byte identities, and
    /// re-enter the authoritative [`CinematicConfig::try_new`] constructor.
    ///
    /// The resolver owns cancellation while acquiring bytes. After return,
    /// `checkpoint` is called before every 64 KiB identity-hash tile and once
    /// after the final (including empty) tile sequence. A checkpoint refusal
    /// publishes no admitted configuration.
    pub fn admit_with_asset_resolver<F, C>(
        &self,
        budget: CinematicAssetAdmissionBudget,
        mut checkpoint: C,
        mut resolve: F,
    ) -> Result<CinematicConfig, CinematicConfigDocumentError>
    where
        F: FnMut(
            CinematicDocumentAssetClass,
            usize,
            &CinematicAssetDeclaration,
        ) -> Result<Vec<u8>, CinematicAssetAccessError>,
        C: FnMut() -> Result<(), CinematicAssetAccessError>,
    {
        if budget.max_asset_bytes == 0
            || budget.max_total_asset_bytes == 0
            || budget.max_asset_bytes > budget.max_total_asset_bytes
        {
            return Err(CinematicConfigDocumentError::InvalidAssetBudget);
        }
        self.quality_profile()?;
        let mut resolved_bytes = 0usize;
        let material_assets = resolve_assets(
            CinematicDocumentAssetClass::Material,
            &self.material_assets,
            budget,
            &mut resolved_bytes,
            &mut checkpoint,
            &mut resolve,
        )?;
        let light_assets = resolve_assets(
            CinematicDocumentAssetClass::Light,
            &self.light_assets,
            budget,
            &mut resolved_bytes,
            &mut checkpoint,
            &mut resolve,
        )?;
        let environment_asset = resolve_asset(
            CinematicDocumentAssetClass::Environment,
            0,
            &self.environment_asset,
            budget,
            &mut resolved_bytes,
            &mut checkpoint,
            &mut resolve,
        )?;
        CinematicConfig::try_new(CinematicConfigInput {
            schema_version: CINEMATIC_CONFIG_SCHEMA_VERSION,
            units: Some(self.units),
            seed: Some(self.seed),
            capabilities: Some(self.capabilities),
            render_budget_profile: Some(self.render_budget_profile),
            audio_budget_profile: Some(self.audio_budget_profile),
            trajectory: self.trajectory,
            timeline: self.timeline,
            camera: self.camera,
            scene_geometry: self.scene_geometry,
            instance_mapping: self.instance_mapping,
            renderer: self.renderer,
            image_pipeline: self.image_pipeline,
            audio_excitation: self.audio_excitation,
            sound_model: self.sound_model,
            microphone: self.microphone,
            room: self.room,
            material_assets,
            light_assets,
            environment_asset,
            artifact_root: self.artifact_root.clone(),
            mux_request: self.mux_request,
        })
        .map_err(CinematicConfigDocumentError::Config)
    }

    /// Canonical complete document text. Declaration order and insignificant
    /// input whitespace/comments do not affect this representation. The final
    /// record has no trailing newline, keeping canonicalization closed at the
    /// exact public byte ceiling.
    #[must_use]
    pub fn to_canonical_text(&self) -> String {
        let mut out = String::new();
        push_line(&mut out, "schema", CINEMATIC_CONFIG_DOCUMENT_SCHEMA);
        push_line(
            &mut out,
            "quality_profile",
            quality_profile_name(self.quality_profile),
        );
        push_line(&mut out, "units", "si-m-kg-s-rad");
        push_line(&mut out, "seed", &self.seed.to_string());
        push_line(
            &mut out,
            "capabilities",
            capabilities_name(self.capabilities),
        );
        for (field, component) in [
            ("render_budget_profile", self.render_budget_profile),
            ("audio_budget_profile", self.audio_budget_profile),
            ("trajectory", self.trajectory),
            ("timeline", self.timeline),
            ("camera", self.camera),
            ("scene_geometry", self.scene_geometry),
            ("instance_mapping", self.instance_mapping),
            ("renderer", self.renderer),
            ("image_pipeline", self.image_pipeline),
            ("audio_excitation", self.audio_excitation),
            ("sound_model", self.sound_model),
            ("microphone", self.microphone),
            ("room", self.room),
        ] {
            push_line(&mut out, field, &component_text(component));
        }
        for declaration in &self.material_assets {
            push_line(&mut out, "material_asset", &asset_text(declaration));
        }
        for declaration in &self.light_assets {
            push_line(&mut out, "light_asset", &asset_text(declaration));
        }
        push_line(
            &mut out,
            "environment_asset",
            &asset_text(&self.environment_asset),
        );
        push_line(
            &mut out,
            "artifact_namespace",
            self.artifact_root.logical_namespace(),
        );
        push_line(&mut out, "artifact_root", self.artifact_root.locator_hint());
        push_line(&mut out, "mux", &mux_text(self.mux_request));
        let trailing = out.pop();
        debug_assert_eq!(trailing, Some('\n'));
        out
    }
}

fn resolve_assets<F, C>(
    class: CinematicDocumentAssetClass,
    declarations: &[CinematicAssetDeclaration],
    budget: CinematicAssetAdmissionBudget,
    resolved_bytes: &mut usize,
    checkpoint: &mut C,
    resolve: &mut F,
) -> Result<Vec<CinematicAssetBinding>, CinematicConfigDocumentError>
where
    F: FnMut(
        CinematicDocumentAssetClass,
        usize,
        &CinematicAssetDeclaration,
    ) -> Result<Vec<u8>, CinematicAssetAccessError>,
    C: FnMut() -> Result<(), CinematicAssetAccessError>,
{
    let mut assets = Vec::new();
    assets
        .try_reserve_exact(declarations.len())
        .map_err(|_| CinematicConfigDocumentError::Capacity)?;
    for (index, declaration) in declarations.iter().enumerate() {
        assets.push(resolve_asset(
            class,
            index,
            declaration,
            budget,
            resolved_bytes,
            checkpoint,
            resolve,
        )?);
    }
    Ok(assets)
}

fn resolve_asset<F, C>(
    class: CinematicDocumentAssetClass,
    index: usize,
    declaration: &CinematicAssetDeclaration,
    budget: CinematicAssetAdmissionBudget,
    resolved_bytes: &mut usize,
    checkpoint: &mut C,
    resolve: &mut F,
) -> Result<CinematicAssetBinding, CinematicConfigDocumentError>
where
    F: FnMut(
        CinematicDocumentAssetClass,
        usize,
        &CinematicAssetDeclaration,
    ) -> Result<Vec<u8>, CinematicAssetAccessError>,
    C: FnMut() -> Result<(), CinematicAssetAccessError>,
{
    let bytes = resolve(class, index, declaration)
        .map_err(|kind| CinematicConfigDocumentError::AssetAccess { class, index, kind })?;
    let next_total = resolved_bytes.checked_add(bytes.len()).ok_or(
        CinematicConfigDocumentError::AssetAccess {
            class,
            index,
            kind: CinematicAssetAccessError::TooLarge,
        },
    )?;
    if bytes.len() > budget.max_asset_bytes || next_total > budget.max_total_asset_bytes {
        return Err(CinematicConfigDocumentError::AssetAccess {
            class,
            index,
            kind: CinematicAssetAccessError::TooLarge,
        });
    }
    *resolved_bytes = next_total;
    let mut hasher = DomainHasher::new(ASSET_BYTES_DOMAIN);
    for tile in bytes.chunks(CINEMATIC_ASSET_HASH_TILE_BYTES) {
        checkpoint().map_err(|kind| CinematicConfigDocumentError::AssetAccess {
            class,
            index,
            kind,
        })?;
        hasher.update(tile);
    }
    checkpoint().map_err(|kind| CinematicConfigDocumentError::AssetAccess {
        class,
        index,
        kind,
    })?;
    let binding = CinematicAssetBinding::from_identity(
        hasher.finalize(),
        declaration.interpretation,
        declaration.version,
        declaration.locator_hint.clone(),
    )
    .map_err(CinematicConfigDocumentError::Config)?;
    if binding.content_identity() != declaration.expected_identity {
        return Err(CinematicConfigDocumentError::AssetIdentityMismatch { class, index });
    }
    Ok(binding)
}

fn take_required(
    fields: &mut BTreeMap<String, String>,
    field: &'static str,
) -> Result<String, CinematicConfigDocumentError> {
    fields
        .remove(field)
        .ok_or_else(|| CinematicConfigDocumentError::MissingField {
            field: field.to_owned(),
        })
}

fn push_asset_declaration(
    output: &mut Vec<CinematicAssetDeclaration>,
    value: &str,
    field: &'static str,
) -> Result<(), CinematicConfigDocumentError> {
    if output.len() >= MAX_CINEMATIC_CONFIG_DOCUMENT_ASSETS_PER_CLASS {
        return Err(CinematicConfigDocumentError::TooManyAssets { field });
    }
    output.push(parse_asset_declaration(value, field)?);
    Ok(())
}

fn parse_component(
    value: &str,
    field: &'static str,
    role: CinematicComponentRole,
) -> Result<CinematicComponentRef, CinematicConfigDocumentError> {
    let (version, identity) = value
        .split_once(':')
        .ok_or(CinematicConfigDocumentError::InvalidField { field })?;
    if identity.contains(':') {
        return Err(CinematicConfigDocumentError::InvalidField { field });
    }
    let version = parse_u32(version, field)?;
    let identity = parse_hash(identity, field)?;
    CinematicComponentRef::try_new(role, identity, version)
        .map_err(CinematicConfigDocumentError::Config)
}

fn parse_asset_declaration(
    value: &str,
    field: &'static str,
) -> Result<CinematicAssetDeclaration, CinematicConfigDocumentError> {
    let mut parts = value.splitn(4, ':');
    let interpretation = parts
        .next()
        .and_then(parse_interpretation)
        .ok_or(CinematicConfigDocumentError::InvalidField { field })?;
    let version = parse_u32(
        parts
            .next()
            .ok_or(CinematicConfigDocumentError::InvalidField { field })?,
        field,
    )?;
    let expected_identity = parse_hash(
        parts
            .next()
            .ok_or(CinematicConfigDocumentError::InvalidField { field })?,
        field,
    )?;
    let locator_hint = parts
        .next()
        .ok_or(CinematicConfigDocumentError::InvalidField { field })?;
    if locator_hint.is_empty()
        || locator_hint.len() > MAX_LOCATOR_BYTES
        || locator_hint.chars().any(char::is_control)
    {
        return Err(CinematicConfigDocumentError::InvalidField { field });
    }
    if version == 0 || expected_identity.as_bytes().iter().all(|byte| *byte == 0) {
        return Err(CinematicConfigDocumentError::InvalidField { field });
    }
    Ok(CinematicAssetDeclaration {
        expected_identity,
        interpretation,
        version,
        locator_hint: locator_hint.to_owned(),
    })
}

fn canonicalize_declarations(
    declarations: &mut [CinematicAssetDeclaration],
    class: CinematicDocumentAssetClass,
) -> Result<(), CinematicConfigDocumentError> {
    declarations.sort_unstable_by_key(CinematicAssetDeclaration::key);
    if declarations
        .windows(2)
        .any(|pair| pair[0].key() == pair[1].key())
    {
        return Err(CinematicConfigDocumentError::DuplicateAsset { class });
    }
    Ok(())
}

fn parse_quality_profile(
    value: &str,
) -> Result<CinematicQualityTier, CinematicConfigDocumentError> {
    match value {
        "storyboard-smoke" => Ok(CinematicQualityTier::StoryboardSmoke),
        "daily-1080p" => Ok(CinematicQualityTier::Daily1080p),
        "qualification-4k-frame" => Ok(CinematicQualityTier::Qualification4kFrame),
        "final-4k" => Ok(CinematicQualityTier::Final4k),
        _ => Err(CinematicConfigDocumentError::InvalidField {
            field: "quality_profile",
        }),
    }
}

const fn quality_profile_name(value: CinematicQualityTier) -> &'static str {
    match value {
        CinematicQualityTier::StoryboardSmoke => "storyboard-smoke",
        CinematicQualityTier::Daily1080p => "daily-1080p",
        CinematicQualityTier::Qualification4kFrame => "qualification-4k-frame",
        CinematicQualityTier::Final4k => "final-4k",
    }
}

fn parse_units(value: &str) -> Result<CinematicConfigUnits, CinematicConfigDocumentError> {
    match value {
        "si-m-kg-s-rad" => Ok(CinematicConfigUnits::SiMetersKilogramsSecondsRadians),
        _ => Err(CinematicConfigDocumentError::InvalidField { field: "units" }),
    }
}

fn parse_capabilities(value: &str) -> Result<CinematicCapabilities, CinematicConfigDocumentError> {
    let mut bits = 0u32;
    for capability in value.split(',') {
        let bit = match capability {
            "render" => CinematicCapabilities::RENDER,
            "audio" => CinematicCapabilities::AUDIO,
            "quarantined-mux" => CinematicCapabilities::QUARANTINED_MUX,
            _ => {
                return Err(CinematicConfigDocumentError::InvalidField {
                    field: "capabilities",
                });
            }
        };
        if bits & bit != 0 {
            return Err(CinematicConfigDocumentError::InvalidField {
                field: "capabilities",
            });
        }
        bits |= bit;
    }
    if bits & (CinematicCapabilities::RENDER | CinematicCapabilities::AUDIO)
        != CinematicCapabilities::RENDER | CinematicCapabilities::AUDIO
    {
        return Err(CinematicConfigDocumentError::InvalidField {
            field: "capabilities",
        });
    }
    CinematicCapabilities::try_new(bits).map_err(CinematicConfigDocumentError::Config)
}

fn capabilities_name(capabilities: CinematicCapabilities) -> &'static str {
    match capabilities.bits() {
        bits if bits == CinematicCapabilities::RENDER | CinematicCapabilities::AUDIO => {
            "render,audio"
        }
        bits if bits
            == CinematicCapabilities::RENDER
                | CinematicCapabilities::AUDIO
                | CinematicCapabilities::QUARANTINED_MUX =>
        {
            "render,audio,quarantined-mux"
        }
        _ => "invalid",
    }
}

fn parse_mux(value: &str) -> Result<CinematicMuxRequest, CinematicConfigDocumentError> {
    if value == "none" {
        return Ok(CinematicMuxRequest::None);
    }
    let mut parts = value.split(':');
    let codec = match parts.next() {
        Some("av1-opus-matroska") => CinematicMuxCodec::Av1OpusMatroska,
        Some("h265-pcm-quicktime") => CinematicMuxCodec::H265PcmQuickTime,
        _ => return Err(CinematicConfigDocumentError::InvalidField { field: "mux" }),
    };
    let version = parse_u32(
        parts
            .next()
            .ok_or(CinematicConfigDocumentError::InvalidField { field: "mux" })?,
        "mux",
    )?;
    let identity = parse_hash(
        parts
            .next()
            .ok_or(CinematicConfigDocumentError::InvalidField { field: "mux" })?,
        "mux",
    )?;
    if parts.next().is_some() {
        return Err(CinematicConfigDocumentError::InvalidField { field: "mux" });
    }
    Ok(CinematicMuxRequest::QuarantinedAdapter {
        adapter_identity: identity,
        adapter_version: version,
        codec,
    })
}

fn mux_text(value: CinematicMuxRequest) -> String {
    match value {
        CinematicMuxRequest::None => "none".to_owned(),
        CinematicMuxRequest::QuarantinedAdapter {
            adapter_identity,
            adapter_version,
            codec,
        } => format!(
            "{}:{adapter_version}:{}",
            match codec {
                CinematicMuxCodec::Av1OpusMatroska => "av1-opus-matroska",
                CinematicMuxCodec::H265PcmQuickTime => "h265-pcm-quicktime",
            },
            adapter_identity.to_hex(),
        ),
    }
}

fn parse_interpretation(value: &str) -> Option<CinematicAssetInterpretation> {
    match value {
        "spectral-reflectance" => Some(CinematicAssetInterpretation::SpectralReflectance),
        "spectral-emission" => Some(CinematicAssetInterpretation::SpectralEmission),
        "linear-scene-texture" => Some(CinematicAssetInterpretation::LinearSceneTexture),
        "display-referred-texture" => Some(CinematicAssetInterpretation::DisplayReferredTexture),
        "geometry-meters" => Some(CinematicAssetInterpretation::GeometryMeters),
        "acoustic-impulse-response-si" => {
            Some(CinematicAssetInterpretation::AcousticImpulseResponseSi)
        }
        _ => None,
    }
}

const fn interpretation_name(value: CinematicAssetInterpretation) -> &'static str {
    match value {
        CinematicAssetInterpretation::SpectralReflectance => "spectral-reflectance",
        CinematicAssetInterpretation::SpectralEmission => "spectral-emission",
        CinematicAssetInterpretation::LinearSceneTexture => "linear-scene-texture",
        CinematicAssetInterpretation::DisplayReferredTexture => "display-referred-texture",
        CinematicAssetInterpretation::GeometryMeters => "geometry-meters",
        CinematicAssetInterpretation::AcousticImpulseResponseSi => "acoustic-impulse-response-si",
    }
}

fn component_text(value: CinematicComponentRef) -> String {
    format!("{}:{}", value.version(), value.identity().to_hex())
}

fn asset_text(value: &CinematicAssetDeclaration) -> String {
    format!(
        "{}:{}:{}:{}",
        interpretation_name(value.interpretation),
        value.version,
        value.expected_identity.to_hex(),
        value.locator_hint,
    )
}

fn push_line(output: &mut String, key: &str, value: &str) {
    output.push_str(key);
    output.push('=');
    output.push_str(value);
    output.push('\n');
}

fn parse_hash(
    value: &str,
    field: &'static str,
) -> Result<ContentHash, CinematicConfigDocumentError> {
    let hash =
        ContentHash::from_hex(value).ok_or(CinematicConfigDocumentError::InvalidField { field })?;
    if hash.as_bytes().iter().all(|byte| *byte == 0) {
        return Err(CinematicConfigDocumentError::InvalidField { field });
    }
    Ok(hash)
}

fn parse_u32(value: &str, field: &'static str) -> Result<u32, CinematicConfigDocumentError> {
    if value.is_empty()
        || !value.bytes().all(|byte| byte.is_ascii_digit())
        || (value.len() > 1 && value.starts_with('0'))
    {
        return Err(CinematicConfigDocumentError::InvalidField { field });
    }
    value
        .parse::<u32>()
        .ok()
        .filter(|value| *value > 0)
        .ok_or(CinematicConfigDocumentError::InvalidField { field })
}

fn parse_u64(value: &str, field: &'static str) -> Result<u64, CinematicConfigDocumentError> {
    if value.is_empty()
        || !value.bytes().all(|byte| byte.is_ascii_digit())
        || (value.len() > 1 && value.starts_with('0'))
    {
        return Err(CinematicConfigDocumentError::InvalidField { field });
    }
    value
        .parse::<u64>()
        .map_err(|_| CinematicConfigDocumentError::InvalidField { field })
}

/// Stable refusal from complete-document decoding or asset-backed admission.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CinematicConfigDocumentError {
    /// Input exceeded the hard config-document byte cap.
    DocumentTooLarge {
        /// Observed bytes.
        bytes: usize,
        /// Hard maximum.
        maximum: usize,
    },
    /// Input was not UTF-8.
    InvalidUtf8,
    /// Physical line cap was exceeded.
    TooManyLines,
    /// A nonblank, noncomment line was not `key=value`.
    MalformedLine {
        /// One-based line number; line contents are intentionally omitted.
        line: usize,
    },
    /// A final bare carriage return made line handling ambiguous.
    InvalidLineEnding,
    /// Schema was absent from the supported closed set.
    UnsupportedSchema,
    /// A field outside the closed grammar appeared.
    UnknownField {
        /// Field name only; values are never retained in the error.
        field: String,
    },
    /// A singleton field appeared more than once.
    DuplicateField {
        /// Duplicate key.
        field: String,
    },
    /// A required field or repeated asset class was absent.
    MissingField {
        /// Missing key.
        field: String,
    },
    /// A known field had malformed or unsupported content.
    InvalidField {
        /// Stable field key.
        field: &'static str,
    },
    /// A repeated asset class exceeded its hard count cap.
    TooManyAssets {
        /// Repeated key.
        field: &'static str,
    },
    /// Two declarations had one semantic asset key.
    DuplicateAsset {
        /// Asset class.
        class: CinematicDocumentAssetClass,
    },
    /// Named profile and one profile reference did not bind the same bytes.
    BudgetProfileMismatch {
        /// Refusing reference field.
        field: &'static str,
    },
    /// Caller could not resolve bounded bytes for one declaration.
    AssetAccess {
        /// Asset class.
        class: CinematicDocumentAssetClass,
        /// Canonical declaration index.
        index: usize,
        /// Bounded access class.
        kind: CinematicAssetAccessError,
    },
    /// Resolved bytes disagreed with the document's expected identity.
    AssetIdentityMismatch {
        /// Asset class.
        class: CinematicDocumentAssetClass,
        /// Canonical declaration index.
        index: usize,
    },
    /// A bounded output allocation failed.
    Capacity,
    /// The caller supplied a zero or internally inconsistent asset envelope.
    InvalidAssetBudget,
    /// Underlying quality-profile admission refused.
    Budget(CinematicBudgetError),
    /// Underlying composition-schema admission refused.
    Config(CinematicConfigError),
}

impl CinematicConfigDocumentError {
    /// Stable machine-readable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::DocumentTooLarge { .. } => "cinematic-document-too-large",
            Self::InvalidUtf8 => "cinematic-document-invalid-utf8",
            Self::TooManyLines => "cinematic-document-too-many-lines",
            Self::MalformedLine { .. } => "cinematic-document-malformed-line",
            Self::InvalidLineEnding => "cinematic-document-invalid-line-ending",
            Self::UnsupportedSchema => "cinematic-document-unsupported-schema",
            Self::UnknownField { .. } => "cinematic-document-unknown-field",
            Self::DuplicateField { .. } => "cinematic-document-duplicate-field",
            Self::MissingField { .. } => "cinematic-document-missing-field",
            Self::InvalidField { .. } => "cinematic-document-invalid-field",
            Self::TooManyAssets { .. } => "cinematic-document-too-many-assets",
            Self::DuplicateAsset { .. } => "cinematic-document-duplicate-asset",
            Self::BudgetProfileMismatch { .. } => "cinematic-document-budget-profile-mismatch",
            Self::AssetAccess { .. } => "cinematic-document-asset-access",
            Self::AssetIdentityMismatch { .. } => "cinematic-document-asset-identity-mismatch",
            Self::Capacity => "cinematic-document-capacity",
            Self::InvalidAssetBudget => "cinematic-document-invalid-asset-budget",
            Self::Budget(error) => error.code(),
            Self::Config(error) => error.code(),
        }
    }

    /// Stable logical field path without locator or source-value disclosure.
    #[must_use]
    pub fn field_path(&self) -> String {
        match self {
            Self::UnknownField { .. } => "config.<unknown>".to_owned(),
            Self::DuplicateField { field } | Self::MissingField { field } => {
                format!("config.{field}")
            }
            Self::InvalidField { field }
            | Self::TooManyAssets { field }
            | Self::BudgetProfileMismatch { field } => format!("config.{field}"),
            Self::DuplicateAsset { class } => {
                format!("config.{}", class.field_name())
            }
            Self::AssetAccess { class, index, .. }
            | Self::AssetIdentityMismatch { class, index } => {
                format!("config.{}[{index}]", class.field_name())
            }
            Self::MalformedLine { line } => format!("config.line[{line}]"),
            Self::UnsupportedSchema => "config.schema".to_owned(),
            Self::Budget(_) => "config.quality_profile".to_owned(),
            Self::Config(_) => "config".to_owned(),
            Self::InvalidAssetBudget => "config.assets".to_owned(),
            Self::DocumentTooLarge { .. }
            | Self::InvalidUtf8
            | Self::TooManyLines
            | Self::InvalidLineEnding
            | Self::Capacity => "config".to_owned(),
        }
    }
}

impl fmt::Display for CinematicConfigDocumentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{} at {}", self.code(), self.field_path())
    }
}

impl std::error::Error for CinematicConfigDocumentError {}
