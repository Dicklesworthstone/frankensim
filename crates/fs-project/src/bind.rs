//! Validate-time resolution of material and interface card bindings
//! (bead f85xj.6.4): every declared region resolves to a matdb
//! `MaterialCard` and every declared interface to an
//! `InterfaceSystemCard`, with the envelope-vs-validity check performed
//! BEFORE any solve — the product's first line of defense against silent
//! extrapolation.
//!
//! The resolver never invents authority. Every property value comes out
//! of `fs-matdb`'s receipted query path (so refusals are matdb's own
//! typed refusals mapped to project violations with fixes), every query
//! receipt is retained for the run ledger, `Unstated` uncertainty is
//! surfaced up front as an advisory (it caps downstream evidence at
//! Estimated), and a card with coexisting conflicting claims resolves
//! only through the explicit `claim` pin recorded in the project file —
//! there is no auto-pick path.
//!
//! Coverage logic: a claim's `ValidityDomain` bounds are per-axis
//! intervals, so a domain that contains both endpoints of the admitted
//! range contains the whole range. The resolver therefore queries each
//! property at BOTH endpoints and additionally requires the SAME claim
//! to be selected at both — two different claims covering one endpoint
//! each would leave the interior ambiguous or stitched, which is a
//! refusal, not a resolution.

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimId, ClaimSet, InterfaceSystemCard, MatDbError, MaterialAnswer, MaterialCard,
    PropertyUsageReceipt, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::Dims;
use fs_regime::RegimeAuditCard;
use fs_scenario::Violation;

use crate::spec::{EntityDecl, ProjectSpec, dims};

/// Canonical matdb property name for bulk thermal conductivity, aligned
/// with `fs-conduction`'s consumption (W/(m·K)).
pub const THERMAL_CONDUCTIVITY_PROPERTY: &str = "thermal-conductivity";

/// SI dimensions of thermal conductivity, `[m, kg, s, K, A, mol]`.
pub const THERMAL_CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);

/// Canonical matdb property name for the interface area-specific thermal
/// contact resistance `R''` (m²·K/W), aligned with `fs-conduction`.
pub const CONTACT_RESISTANCE_PROPERTY: &str = "area-specific-thermal-contact-resistance";

/// SI dimensions of `R''`, `[m, kg, s, K, A, mol]`.
pub const CONTACT_RESISTANCE_DIMS: Dims = Dims([0, -1, 3, 1, 0, 0]);

/// Validity/query axis name for temperature, aligned with matdb cards.
pub const TEMPERATURE_AXIS: &str = "T";

/// One property the product requires a binding's card to answer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RequiredProperty {
    /// Property name as spelled on matdb claims.
    pub property: String,
    /// Required SI dimensions of the answered sample.
    pub dims: Dims,
}

/// What the resolution demands from every binding: which properties each
/// material and interface card must answer, over which axis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BindingRequirements {
    /// The validity axis queried (v1 resolves over temperature only).
    pub temperature_axis: String,
    /// Properties every region's material card must answer.
    pub material_properties: Vec<RequiredProperty>,
    /// Properties every interface's system card must answer.
    pub interface_properties: Vec<RequiredProperty>,
}

impl BindingRequirements {
    /// The steady thermal-conduction requirement set of the ratified
    /// cooling vertical: bulk conductivity per region, area-specific
    /// contact resistance per interface.
    #[must_use]
    pub fn thermal_steady_v1() -> BindingRequirements {
        BindingRequirements {
            temperature_axis: TEMPERATURE_AXIS.to_string(),
            material_properties: vec![RequiredProperty {
                property: THERMAL_CONDUCTIVITY_PROPERTY.to_string(),
                dims: THERMAL_CONDUCTIVITY_DIMS,
            }],
            interface_properties: vec![RequiredProperty {
                property: CONTACT_RESISTANCE_PROPERTY.to_string(),
                dims: CONTACT_RESISTANCE_DIMS,
            }],
        }
    }
}

/// Caller-supplied card store keyed by content hash. The library computes
/// each key from the card's own `content_hash()`, so a key can never lie
/// about its card; authenticity of the collection itself (WHO supplied
/// these cards) is the caller's trust channel, not this type's claim.
#[derive(Debug, Default)]
pub struct CardLibrary {
    materials: BTreeMap<String, MaterialCard>,
    interfaces: BTreeMap<String, InterfaceSystemCard>,
}

impl CardLibrary {
    /// An empty library.
    #[must_use]
    pub fn new() -> CardLibrary {
        CardLibrary::default()
    }

    /// Store a material card under its own content hash; returns the hex
    /// key a project binding must reference.
    pub fn insert_material(&mut self, card: MaterialCard) -> String {
        let key = card.content_hash().to_hex();
        self.materials.insert(key.clone(), card);
        key
    }

    /// Store an interface-system card under its own content hash;
    /// returns the hex key a project binding must reference.
    pub fn insert_interface(&mut self, card: InterfaceSystemCard) -> String {
        let key = card.content_hash().to_hex();
        self.interfaces.insert(key.clone(), card);
        key
    }

    /// Look up a material card by its 64-hex content hash.
    #[must_use]
    pub fn material(&self, hex: &str) -> Option<&MaterialCard> {
        self.materials.get(hex)
    }

    /// Look up an interface-system card by its 64-hex content hash.
    #[must_use]
    pub fn interface(&self, hex: &str) -> Option<&InterfaceSystemCard> {
        self.interfaces.get(hex)
    }
}

/// One retained matdb usage receipt, ledger-ready: the receipt, its
/// canonical bytes, its content hash, and a deterministic human context.
#[derive(Debug, Clone, PartialEq)]
pub struct RetainedReceipt {
    /// Deterministic context line (target, property, query point).
    pub context: String,
    /// The receipt itself, replayable via `ClaimSet::verify_receipt`.
    pub receipt: PropertyUsageReceipt,
    /// Hex content hash of the canonical receipt bytes.
    pub receipt_hash: String,
    /// The canonical portable bytes.
    pub bytes: Vec<u8>,
}

/// One property resolved through the receipted query path at both
/// endpoints of the binding's admitted range.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedProperty {
    /// Property name.
    pub property: String,
    /// Value at the low endpoint of the admitted range.
    pub value_lo: f64,
    /// Value at the high endpoint of the admitted range.
    pub value_hi: f64,
    /// Answered sample dimensions.
    pub dims: Dims,
    /// The claim's stated uncertainty model.
    pub uncertainty: UncertaintyModel,
    /// True when the uncertainty is `Unstated`: downstream evidence for
    /// anything consuming this property is capped at Estimated.
    pub unstated_uncertainty: bool,
    /// Hex content hash of the selected claim (same at both endpoints).
    pub selected_claim: String,
    /// Owner-neutral product-audit projection bound to the exact project card,
    /// selected claim, portable query schema, and claim validity domain.
    pub regime_card: RegimeAuditCard,
    /// The selected claim's provenance source citation.
    pub provenance_source: String,
    /// The selected claim's provenance license.
    pub provenance_license: String,
    /// Retained receipt for the low-endpoint query.
    pub receipt_lo: RetainedReceipt,
    /// Retained receipt for the high-endpoint query.
    pub receipt_hi: RetainedReceipt,
}

/// What a resolved binding is bound to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BindingTarget {
    /// A declared region entity.
    Region(String),
    /// A declared interface entity.
    Interface(String),
}

impl BindingTarget {
    fn describe(&self) -> String {
        match self {
            BindingTarget::Region(name) => format!("region `{name}`"),
            BindingTarget::Interface(name) => format!("interface `{name}`"),
        }
    }
}

/// One fully resolved binding row of the material resolution table.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedBinding {
    /// The declared entity this binding resolves.
    pub target: BindingTarget,
    /// Hex content hash of the bound card.
    pub card: String,
    /// The card's own identity rendering (material state or interface
    /// system summary).
    pub card_identity: String,
    /// Source channel declared in the project file (recorded, not
    /// authenticated by this resolver).
    pub declared_source: String,
    /// Low endpoint (K) of the range the properties were resolved over.
    pub range_lo: f64,
    /// High endpoint (K) of the range the properties were resolved over.
    pub range_hi: f64,
    /// The explicit claim pin, when the project file recorded one.
    pub pinned_claim: Option<String>,
    /// Resolved properties in requirement order.
    pub properties: Vec<ResolvedProperty>,
}

/// A non-refusing finding surfaced up front (validation still passes).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Advisory {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// What was found.
    pub what: String,
    /// Why it matters downstream.
    pub note: String,
}

/// The complete material resolution: violations (empty means every
/// binding resolved), advisories, and the per-binding resolution table
/// with retained receipts.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MaterialResolution {
    /// Refusals, each with an actionable fix.
    pub violations: Vec<Violation>,
    /// Up-front warnings that do not refuse (e.g. Unstated uncertainty).
    pub advisories: Vec<Advisory>,
    /// Resolution table rows, in project declaration order (materials
    /// before interfaces).
    pub bindings: Vec<ResolvedBinding>,
}

impl MaterialResolution {
    /// True when every binding resolved with no refusal.
    #[must_use]
    pub fn admissible(&self) -> bool {
        self.violations.is_empty()
    }

    /// Every retained receipt, in table order (low endpoint before high),
    /// ready for ledger retention.
    pub fn receipts(&self) -> impl Iterator<Item = &RetainedReceipt> {
        self.bindings.iter().flat_map(|binding| {
            binding
                .properties
                .iter()
                .flat_map(|property| [&property.receipt_lo, &property.receipt_hi])
        })
    }

    /// Exact matdb claim cards consumed by this resolution, sorted and
    /// deduplicated by their content-bound audit identity.
    ///
    /// These projections carry no invented discrepancy, ambition,
    /// calibration, or validation authority. Their name binds the project card
    /// and selected property claim hashes; their version binds the portable
    /// query-receipt schema; their validity is cloned from that immutable
    /// selected claim.
    #[must_use]
    pub fn regime_audit_cards(&self) -> Vec<RegimeAuditCard> {
        let mut cards = BTreeMap::new();
        for binding in &self.bindings {
            for property in &binding.properties {
                cards
                    .entry(property.regime_card.name.clone())
                    .or_insert_with(|| property.regime_card.clone());
            }
        }
        cards.into_values().collect()
    }

    /// Deterministic text table: one line per resolved property with the
    /// full region -> card -> property -> uncertainty -> receipts chain,
    /// for CLI logs and the e2e battery.
    #[must_use]
    pub fn render_table(&self) -> String {
        use core::fmt::Write as _;

        let mut out = String::new();
        for binding in &self.bindings {
            for property in &binding.properties {
                let pin = binding.pinned_claim.as_deref().unwrap_or("-");
                let _ = writeln!(
                    out,
                    "{} | card {} ({}) | source {} | range [{}, {}] K | pin {} | {} = [{}, {}] {} | uncertainty {} | claim {} from {} ({}) | regime-card {}@{} | receipts {},{}",
                    binding.target.describe(),
                    binding.card,
                    binding.card_identity,
                    binding.declared_source,
                    binding.range_lo,
                    binding.range_hi,
                    pin,
                    property.property,
                    property.value_lo,
                    property.value_hi,
                    property.dims.unit_string(),
                    render_uncertainty(&property.uncertainty),
                    property.selected_claim,
                    property.provenance_source,
                    property.provenance_license,
                    property.regime_card.name,
                    property.regime_card.version,
                    property.receipt_lo.receipt_hash,
                    property.receipt_hi.receipt_hash,
                );
            }
        }
        out
    }
}

fn render_uncertainty(uncertainty: &UncertaintyModel) -> String {
    match uncertainty {
        UncertaintyModel::Unstated => {
            "UNSTATED (downstream evidence capped at Estimated)".to_string()
        }
        UncertaintyModel::HalfWidth {
            half_width,
            confidence,
        } => format!("±{half_width} @ confidence {confidence}"),
        UncertaintyModel::RelativeHalfWidth {
            fraction,
            confidence,
        } => format!("±{fraction} (relative) @ confidence {confidence}"),
    }
}

fn violation(code: &'static str, what: impl Into<String>, fix: impl Into<String>) -> Violation {
    Violation {
        code,
        what: what.into(),
        fix: fix.into(),
    }
}

/// The per-target operating range the material data must cover: from the
/// envelope's low ambient up to the hotter of the high ambient and every
/// thermal limit declared on the target region(s). A region with a 398 K
/// junction limit needs data valid to 398 K to certify compliance —
/// the solve is not allowed to discover mid-run that its material table
/// ends below the requirement it is judging.
struct RequiredRange {
    lo: f64,
    hi: f64,
}

/// Resolve every material and interface binding in `spec` against the
/// caller's `library` under `requirements`. Structural validation
/// ([`ProjectSpec::validate`]) is a separate, prior concern: this
/// resolver re-guards only what it must to stay total, and names its own
/// preconditions instead of panicking on garbage.
#[must_use]
pub fn resolve_bindings(
    spec: &ProjectSpec,
    library: &CardLibrary,
    requirements: &BindingRequirements,
) -> MaterialResolution {
    let mut resolution = MaterialResolution::default();

    let (Some(envelope), Some(assembly), Some(materials), Some(interface_cards)) = (
        &spec.envelope,
        &spec.assembly,
        &spec.materials,
        &spec.interface_cards,
    ) else {
        resolution.violations.push(violation(
            "project-binding-preconditions",
            "material resolution needs the `envelope`, `assembly`, `materials`, and `interface-cards` sections",
            "declare the missing sections; run structural validation first to see every omission",
        ));
        return resolution;
    };

    let ambient_lo = envelope.ambient_lo;
    let ambient_hi = envelope.ambient_hi;
    if ambient_lo.dims != dims::TEMPERATURE
        || ambient_hi.dims != dims::TEMPERATURE
        || !ambient_lo.value.is_finite()
        || !ambient_hi.value.is_finite()
        || ambient_lo.value > ambient_hi.value
    {
        resolution.violations.push(violation(
            "project-binding-preconditions",
            "the operating envelope's ambient range is not a finite, ordered temperature range",
            "fix the envelope section first; the envelope-vs-domain check is meaningless without it",
        ));
        return resolution;
    }

    // Entity view: region names, and interface names with their sides.
    let mut region_names: BTreeSet<&str> = BTreeSet::new();
    let mut interface_sides: BTreeMap<&str, (&str, &str)> = BTreeMap::new();
    for decl in assembly {
        match decl {
            EntityDecl::Region { name, .. } => {
                region_names.insert(name.as_str());
            }
            EntityDecl::Interface { name, from, to, .. } => {
                interface_sides.insert(name.as_str(), (from.as_str(), to.as_str()));
            }
            EntityDecl::Assembly { .. } | EntityDecl::Part { .. } => {}
        }
    }

    // Per-region required ceiling: high ambient raised by every thermal
    // limit declared on the region.
    let mut required_hi: BTreeMap<&str, f64> = region_names
        .iter()
        .map(|name| (*name, ambient_hi.value))
        .collect();
    if let Some(limits) = &spec.requirements {
        for limit in limits {
            if limit.limit.dims == dims::TEMPERATURE
                && limit.limit.value.is_finite()
                && let Some(ceiling) = required_hi.get_mut(limit.region.as_str())
            {
                *ceiling = ceiling.max(limit.limit.value);
            }
        }
    }

    resolve_material_bindings(
        &mut resolution,
        materials,
        library,
        requirements,
        &region_names,
        &required_hi,
        ambient_lo.value,
    );
    resolve_interface_bindings(
        &mut resolution,
        interface_cards,
        library,
        requirements,
        &interface_sides,
        &required_hi,
        ambient_lo.value,
        ambient_hi.value,
    );
    resolution
}

fn resolve_material_bindings(
    resolution: &mut MaterialResolution,
    materials: &[crate::spec::MaterialBinding],
    library: &CardLibrary,
    requirements: &BindingRequirements,
    region_names: &BTreeSet<&str>,
    required_hi: &BTreeMap<&str, f64>,
    ambient_lo: f64,
) {
    let mut bound: BTreeSet<&str> = BTreeSet::new();
    for binding in materials {
        let region = binding.region.as_str();
        if !region_names.contains(region) {
            resolution.violations.push(violation(
                "project-material-target-kind",
                format!(
                    "material binding targets `{region}`, which is not a declared region entity"
                ),
                "bind materials to region entities; parts and assemblies carry no material of their own",
            ));
            continue;
        }
        if !bound.insert(region) {
            resolution.violations.push(violation(
                "project-material-binding-duplicate",
                format!("region `{region}` is bound to more than one material card"),
                "bind each region exactly once; a region has one manufactured state",
            ));
            continue;
        }
        let target = BindingTarget::Region(binding.region.clone());
        let Some(card) = library.material(&binding.card) else {
            resolution.violations.push(violation(
                "project-material-card-unknown",
                format!(
                    "{} references card {}, which is not in the supplied card library",
                    target.describe(),
                    binding.card
                ),
                "load the pack that carries this card, or correct the card hash in the binding",
            ));
            continue;
        };
        let card_identity = card.id().to_string();
        if binding.state != card_identity {
            resolution.violations.push(violation(
                "project-material-state-mismatch",
                format!(
                    "{} declares state `{}` but card {} is `{card_identity}`",
                    target.describe(),
                    binding.state,
                    binding.card
                ),
                format!(
                    "spell the manufactured state exactly as the card does: `{card_identity}` (a material is a manufactured state, never a name)"
                ),
            ));
            continue;
        }
        let range = material_range(resolution, binding, required_hi, ambient_lo);
        let Some(range) = range else {
            continue;
        };
        let pin = parse_pin(resolution, &target, binding.claim.as_deref());
        let Ok(pin) = pin else {
            continue;
        };
        resolve_card_properties(
            resolution,
            target,
            card.claims(),
            card_identity,
            binding.card.clone(),
            binding.source.clone(),
            &range,
            pin,
            &requirements.temperature_axis,
            &requirements.material_properties,
        );
    }
    for region in region_names {
        if !bound.contains(region) {
            resolution.violations.push(violation(
                "project-material-unbound-region",
                format!("region `{region}` has no material binding"),
                "every region needs exactly one material card with state, admitted range, and source",
            ));
        }
    }
}

/// The admitted range a material binding declares, checked to be a
/// finite ordered temperature range that covers the region's required
/// operating range.
fn material_range(
    resolution: &mut MaterialResolution,
    binding: &crate::spec::MaterialBinding,
    required_hi: &BTreeMap<&str, f64>,
    ambient_lo: f64,
) -> Option<RequiredRange> {
    let lo = binding.temp_lo;
    let hi = binding.temp_hi;
    if lo.dims != dims::TEMPERATURE
        || hi.dims != dims::TEMPERATURE
        || !lo.value.is_finite()
        || !hi.value.is_finite()
        || lo.value > hi.value
    {
        resolution.violations.push(violation(
            "project-binding-range-invalid",
            format!(
                "the admitted temperature range for `{}` is not a finite, ordered temperature range",
                binding.region
            ),
            "state the admitted range as two ordered finite temperatures",
        ));
        return None;
    }
    let needed_hi = required_hi
        .get(binding.region.as_str())
        .copied()
        .unwrap_or(f64::NEG_INFINITY);
    if lo.value > ambient_lo || hi.value < needed_hi {
        resolution.violations.push(violation(
            "project-material-envelope-uncovered",
            format!(
                "region `{}` must be covered over [{ambient_lo}, {needed_hi}] K (envelope ambient range raised to its thermal limits) but the binding admits only [{}, {}] K",
                binding.region, lo.value, hi.value
            ),
            "widen the admitted range only if the card's data actually covers it; otherwise choose a card qualified for this envelope or change the design envelope",
        ));
        return None;
    }
    Some(RequiredRange {
        lo: lo.value,
        hi: hi.value,
    })
}

#[allow(clippy::too_many_arguments)] // one call site; the argument list IS the binding row
fn resolve_interface_bindings(
    resolution: &mut MaterialResolution,
    interface_cards: &[crate::spec::InterfaceCardBinding],
    library: &CardLibrary,
    requirements: &BindingRequirements,
    interface_sides: &BTreeMap<&str, (&str, &str)>,
    required_hi: &BTreeMap<&str, f64>,
    ambient_lo: f64,
    ambient_hi: f64,
) {
    let mut bound: BTreeSet<&str> = BTreeSet::new();
    for binding in interface_cards {
        let name = binding.interface.as_str();
        let Some((from, to)) = interface_sides.get(name) else {
            resolution.violations.push(violation(
                "project-interface-target-kind",
                format!(
                    "interface card binding targets `{name}`, which is not a declared interface entity"
                ),
                "bind interface-system cards to interface entities declared in the assembly",
            ));
            continue;
        };
        if !bound.insert(name) {
            resolution.violations.push(violation(
                "project-interface-binding-duplicate",
                format!("interface `{name}` is bound to more than one interface-system card"),
                "bind each interface exactly once",
            ));
            continue;
        }
        let target = BindingTarget::Interface(binding.interface.clone());
        let Some(card) = library.interface(&binding.card) else {
            resolution.violations.push(violation(
                "project-interface-card-unknown",
                format!(
                    "{} references card {}, which is not in the supplied card library",
                    target.describe(),
                    binding.card
                ),
                "load the pack that carries this card, or correct the card hash in the binding",
            ));
            continue;
        };
        // An interface must be valid wherever either side can operate.
        let hi = required_hi
            .get(from)
            .copied()
            .unwrap_or(ambient_hi)
            .max(required_hi.get(to).copied().unwrap_or(ambient_hi));
        let range = RequiredRange { lo: ambient_lo, hi };
        let pin = parse_pin(resolution, &target, binding.claim.as_deref());
        let Ok(pin) = pin else {
            continue;
        };
        let card_identity = format!(
            "{} ~ {} [{}]",
            card.surface_a().material,
            card.surface_b().material,
            card.medium()
        );
        resolve_card_properties(
            resolution,
            target,
            card.claims(),
            card_identity,
            binding.card.clone(),
            binding.source.clone(),
            &range,
            pin,
            &requirements.temperature_axis,
            &requirements.interface_properties,
        );
    }
    for name in interface_sides.keys() {
        if !bound.contains(name) {
            resolution.violations.push(violation(
                "project-interface-unbound",
                format!("interface `{name}` has no interface-system card binding"),
                "every declared interface needs exactly one TIM/contact system card",
            ));
        }
    }
}

/// Parse the optional explicit claim pin. `Err(())` means a violation was
/// already recorded.
fn parse_pin(
    resolution: &mut MaterialResolution,
    target: &BindingTarget,
    pin: Option<&str>,
) -> Result<Option<ClaimId>, ()> {
    let Some(text) = pin else {
        return Ok(None);
    };
    let Some(hash) = ContentHash::from_hex(text) else {
        resolution.violations.push(violation(
            "project-binding-pin-malformed",
            format!(
                "{} pins claim `{text}`, which is not a 64-hex claim content hash",
                target.describe()
            ),
            "pin the exact claim by its full content hash as listed in the conflict refusal",
        ));
        return Err(());
    };
    Ok(Some(ClaimId(hash)))
}

#[allow(clippy::too_many_arguments)] // two call sites; the argument list IS the binding row
fn resolve_card_properties(
    resolution: &mut MaterialResolution,
    target: BindingTarget,
    claims: &ClaimSet,
    card_identity: String,
    card_hex: String,
    declared_source: String,
    range: &RequiredRange,
    pin: Option<ClaimId>,
    axis: &str,
    properties: &[RequiredProperty],
) {
    let mut resolved = Vec::new();
    let mut clean = true;
    for property in properties {
        match resolve_property(
            resolution, &target, claims, &card_hex, range, pin, axis, property,
        ) {
            Some(row) => resolved.push(row),
            None => clean = false,
        }
    }
    if clean {
        resolution.bindings.push(ResolvedBinding {
            target,
            card: card_hex,
            card_identity,
            declared_source,
            range_lo: range.lo,
            range_hi: range.hi,
            pinned_claim: pin.map(|claim| claim.0.to_hex()),
            properties: resolved,
        });
    }
}

fn resolve_property(
    resolution: &mut MaterialResolution,
    target: &BindingTarget,
    claims: &ClaimSet,
    card_hex: &str,
    range: &RequiredRange,
    pin: Option<ClaimId>,
    axis: &str,
    property: &RequiredProperty,
) -> Option<ResolvedProperty> {
    let low = query_endpoint(resolution, target, claims, axis, range.lo, pin, property)?;
    let high = query_endpoint(resolution, target, claims, axis, range.hi, pin, property)?;

    if low.receipt.selected != high.receipt.selected {
        resolution.violations.push(violation(
            "project-binding-domain-split",
            format!(
                "{}: no single `{}` claim covers [{}, {}] K — claim {} answers the low end and claim {} the high end",
                target.describe(),
                property.property,
                range.lo,
                range.hi,
                low.receipt.selected.0.to_hex(),
                high.receipt.selected.0.to_hex()
            ),
            "bind a card revision whose claim covers the whole admitted range, or pin one claim that does; stitching two claims across the range is not a resolution",
        ));
        return None;
    }
    let sample = &low.evidence.value;
    if sample.dims != property.dims {
        resolution.violations.push(violation(
            "project-binding-property-dims",
            format!(
                "{}: `{}` answered with dimensions {:?}, expected {:?}",
                target.describe(),
                property.property,
                sample.dims,
                property.dims
            ),
            "the card's claim carries the wrong dimensions for this property; fix the card data, not the consumer",
        ));
        return None;
    }
    if sample.uncertainty == UncertaintyModel::Unstated {
        resolution.advisories.push(Advisory {
            code: "binding-uncertainty-unstated",
            what: format!(
                "{}: `{}` from claim {} has UNSTATED uncertainty",
                target.describe(),
                property.property,
                low.receipt.selected.0.to_hex()
            ),
            note: "every downstream result consuming this property is capped at Estimated; state a half-width on the claim to lift the cap".to_string(),
        });
    }

    let selected = low.receipt.selected;
    let Some((_, selected_claim)) = claims
        .claims_for(&property.property)
        .into_iter()
        .find(|(id, _)| *id == selected)
    else {
        resolution.violations.push(violation(
            "project-binding-selected-claim-missing",
            format!(
                "{}: matdb selected claim {} for `{}`, but the immutable claim set no longer exposes it",
                target.describe(),
                selected.0.to_hex(),
                property.property,
            ),
            "refuse the resolution and investigate claim-set/query identity corruption",
        ));
        return None;
    };
    let selected_hex = selected.0.to_hex();
    let regime_card = RegimeAuditCard::new(
        format!("fs-matdb:claim:{card_hex}:{selected_hex}"),
        format!("property-usage-receipt-v{}", low.receipt.schema_version),
        selected_claim.validity.clone(),
    );

    let receipt_lo = retain_receipt(resolution, target, &low, range.lo, property)?;
    let receipt_hi = retain_receipt(resolution, target, &high, range.hi, property)?;
    Some(ResolvedProperty {
        property: property.property.clone(),
        value_lo: low.evidence.value.value,
        value_hi: high.evidence.value.value,
        dims: sample.dims,
        uncertainty: sample.uncertainty.clone(),
        unstated_uncertainty: sample.uncertainty == UncertaintyModel::Unstated,
        selected_claim: selected_hex,
        regime_card,
        provenance_source: selected_claim.provenance.source.clone(),
        provenance_license: selected_claim.provenance.license.clone(),
        receipt_lo,
        receipt_hi,
    })
}

fn query_endpoint(
    resolution: &mut MaterialResolution,
    target: &BindingTarget,
    claims: &ClaimSet,
    axis: &str,
    at: f64,
    pin: Option<ClaimId>,
    property: &RequiredProperty,
) -> Option<MaterialAnswer> {
    let point = match QueryPoint::new().with(axis, at) {
        Ok(point) => point,
        Err(error) => {
            resolution.violations.push(violation(
                "project-binding-query",
                format!(
                    "{}: cannot form the query point {axis} = {at}: {error}",
                    target.describe()
                ),
                "the admitted range must be finite; fix the binding's range",
            ));
            return None;
        }
    };
    let answer = match pin {
        Some(pinned) => claims.query_pinned(&property.property, &point, pinned),
        None => claims.query(&property.property, &point, SelectionPolicy::SingleClaimOnly),
    };
    match answer {
        Ok(answer) => Some(answer),
        Err(error) => {
            resolution
                .violations
                .push(query_violation(target, property, axis, at, &error));
            None
        }
    }
}

/// Map a matdb refusal to a project violation with a product-level fix.
/// The refusal text itself travels in `what`, so the user sees matdb's
/// own diagnosis plus what to do about it in the project file.
fn query_violation(
    target: &BindingTarget,
    property: &RequiredProperty,
    axis: &str,
    at: f64,
    error: &MatDbError,
) -> Violation {
    let context = format!(
        "{}: `{}` at {axis} = {at} K: {error}",
        target.describe(),
        property.property
    );
    match error {
        MatDbError::UnknownProperty { .. } => violation(
            "project-binding-property-missing",
            context,
            format!(
                "the bound card carries no `{}` claims; bind a card that states this property",
                property.property
            ),
        ),
        MatDbError::NoClaimInDomain { .. } | MatDbError::OutsideKnotSpan { .. } => violation(
            "project-binding-domain-uncovered",
            context,
            "the card's validity does not cover the admitted range: narrow the binding's admitted range only if the envelope allows it, or bind a card qualified for these temperatures — the resolver never extrapolates",
        ),
        MatDbError::AmbiguousSelection { candidates, .. } => violation(
            "project-binding-claims-conflict",
            context,
            format!(
                "the card carries {} coexisting claims for this property and the system never auto-picks: record an explicit `:claim <hex>` pin in the binding, choosing among: {}",
                candidates.len(),
                candidates
                    .iter()
                    .map(|claim| claim.0.to_hex())
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        ),
        MatDbError::PinnedClaimUnknown { .. } => violation(
            "project-binding-pin-unknown",
            context,
            "the pinned claim is not on the bound card: the pin is stale or belongs to another card revision; re-run without the pin to list the candidates",
        ),
        MatDbError::PinnedClaimOutOfDomain { .. } => violation(
            "project-binding-pin-domain",
            context,
            "the pinned claim does not cover the admitted range; pin a claim valid over the whole range — a pin never bypasses the extrapolation refusal",
        ),
        MatDbError::MissingQueryAxis { .. } => violation(
            "project-binding-axis",
            context,
            "the card's claims depend on an axis this resolution does not query (v1 queries temperature only); bind a card parameterized by temperature alone",
        ),
        _ => violation(
            "project-binding-query",
            context,
            "the matdb query refused; the upstream refusal names the exact cause",
        ),
    }
}

fn retain_receipt(
    resolution: &mut MaterialResolution,
    target: &BindingTarget,
    answer: &MaterialAnswer,
    at: f64,
    property: &RequiredProperty,
) -> Option<RetainedReceipt> {
    let bytes = match answer.receipt.to_bytes() {
        Ok(bytes) => bytes,
        Err(error) => {
            resolution.violations.push(violation(
                "project-binding-receipt",
                format!(
                    "{}: the `{}` usage receipt is not portable: {error:?}",
                    target.describe(),
                    property.property
                ),
                "a non-portable receipt cannot be retained for replay; this is an upstream card/claim defect",
            ));
            return None;
        }
    };
    let receipt_hash = answer.receipt.content_hash().to_hex();
    Some(RetainedReceipt {
        context: format!("{} {} @ {at} K", target.describe(), property.property),
        receipt: answer.receipt.clone(),
        receipt_hash,
        bytes,
    })
}
