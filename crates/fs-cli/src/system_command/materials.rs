//! Material-card preparation for the existing generic thermal components and
//! card-backed resistors.
//!
//! This cold boundary resolves original source support once. Thermal trials
//! use the material owners' analytical caloric and conductivity laws, and a
//! resistor's trials read the resistivity source line its Ohmic node retains;
//! they do not query a database, reconstruct source receipts, or freeze Cp, k
//! or the resistivity at the initial temperature.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::io::Read as _;
use std::path::Path;
use std::sync::Arc;

use fs_blake3::hash_bytes;
use fs_conduction::transient::enthalpy::material::{
    ExactMaterialEnthalpySource, MaterialEnthalpyTransport,
};
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialCard, NormalizedMaterialCardPack, PropertyKey,
    PropertyLineSlice, PropertyUsageReceipt, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::conductor::graph::MaterialOhmicNode;
use fs_material::conductor::{
    ElectricalPropertyBasis, ResolvedConductorState, resolve_conductor_state_with_property,
};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve,
    SolidLiquidPhase,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_material::thermal::conduction::ConductivityCurve;
use fs_material::thermal::conduction::source::ResolvedConductionLayer;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};

use crate::cards::{MAX_CARD_PACK_BYTES, MAX_CARD_PACK_SOURCE_BYTES};
use crate::material_state_input::{quantity, state_point};

use super::{J, fields, num, number, positive, quote};

// Match the existing file-driven electrothermal command's source budgets.
const MAX_MATERIALS: usize = 64;
const MAX_TOTAL_PACK_BYTES: u64 = 64 * 1024 * 1024;

pub(super) type Materials = BTreeMap<String, PreparedMaterial>;

/// Owned original laws, with the actual source receipts retained alongside.
pub(super) struct PreparedMaterial {
    pub caloric: EquilibriumEnthalpyPhaseCurve,
    pub conductivity: Option<ConductivityCurve>,
    conductivity_source: Option<Arc<PropertyLineSlice>>,
    report: String,
    source: CardSource,
}

/// The card and the selection the laws were prepared from: its property
/// keys, the declared context and span, the phase and the enthalpy datum.
struct CardSource {
    card: MaterialCard,
    heat_capacity: PropertyKey,
    density: PropertyKey,
    conductivity: Option<PropertyKey>,
    /// The card's resistivity resolved at the span's lower end with the
    /// context unchanged, if the declaration names one.
    resistivity: Option<ResolvedConductorState>,
    axis: String,
    /// The declared context without the temperature, and the temperature's
    /// declared quantity, if any.
    context: QueryPoint,
    temperature_kind: Option<QuantitySpec>,
    lower: QueryPoint,
    upper: QueryPoint,
    /// The declared temperature span [K].
    span_k: [f64; 2],
    selection: SelectionPolicy,
    phase: SolidLiquidPhase,
    lower_specific_enthalpy_j_kg: f64,
}

/// The declared context with its temperature axis at `temperature_k`, with
/// the declared quantity if there is one.
fn at_temperature(
    context: &QueryPoint,
    axis: &str,
    kind: Option<QuantitySpec>,
    temperature_k: f64,
) -> Result<QueryPoint, MatDbError> {
    match kind {
        Some(q) => context.clone().with_quantity(axis, q, temperature_k),
        None => context.clone().with(axis, temperature_k),
    }
}

impl PreparedMaterial {
    /// fs-conduction's paired storage and transport of this material,
    /// prepared by the owner from the same card, keys, context and span:
    /// the chart integrated from the card's Cp(T) and density, the card's
    /// k(T) as isotropic conduction, and the density read once at
    /// `reference_k` with the context unchanged, which fixes the masses.
    ///
    /// # Errors
    /// A material whose declaration names no conductivity, and every
    /// refusal of the owner's (a reference temperature off the span among
    /// them), as text.
    pub fn enthalpy_transport(
        &self,
        reference_k: f64,
    ) -> Result<MaterialEnthalpyTransport, String> {
        let source = &self.source;
        let conductivity = source
            .conductivity
            .as_ref()
            .ok_or("the declaration names no conductivity_property")?;
        let reference = at_temperature(
            &source.context,
            &source.axis,
            source.temperature_kind,
            reference_k,
        )
        .map_err(|error| error.to_string())?;
        MaterialEnthalpyTransport::from_isotropic_card_exact(
            &source.card,
            ExactMaterialEnthalpySource {
                caloric: HeatCapacitySourceQuery {
                    heat_capacity: &source.heat_capacity,
                    density: &source.density,
                    temperature_axis: &source.axis,
                    lower: &source.lower,
                    upper: &source.upper,
                    selection: source.selection,
                },
                conductivity,
                reference: &reference,
                phase: source.phase,
                lower_specific_enthalpy_j_kg: source.lower_specific_enthalpy_j_kg,
            },
        )
        .map_err(|error| error.to_string())
    }

    /// fs-material's Ohmic node of this material's resistivity over the
    /// declared span, for a uniform conductor of `length_m` and `area_m2`:
    /// the card's `rho(T) L / A`, read again at each temperature asked, with
    /// every other coordinate at the declared context.
    ///
    /// # Errors
    /// A material whose declaration names no resistivity, and every refusal
    /// of the node's (a geometry that is not positive and finite among them),
    /// as text.
    pub fn ohmic_node(&self, length_m: f64, area_m2: f64) -> Result<MaterialOhmicNode, String> {
        let source = &self.source;
        let original = source
            .resistivity
            .clone()
            .ok_or("the material's declaration names no resistivity_property")?;
        MaterialOhmicNode::try_new(
            &source.card,
            original,
            &source.axis,
            source.span_k,
            length_m,
            area_m2,
        )
        .map_err(|error| error.to_string())
    }

    /// The actual queries performed while preparing the original source laws.
    /// These are not queries at each runtime trial, which uses the compiled laws.
    pub fn receipts(&self) -> impl Iterator<Item = &PropertyUsageReceipt> {
        self.caloric
            .source_receipts()
            .iter()
            .chain(
                self.conductivity_source
                    .iter()
                    .flat_map(|line| line.knots.iter().map(|(_, answer)| &answer.receipt)),
            )
            .chain(
                self.source
                    .resistivity
                    .iter()
                    .map(|state| &state.source_property().answer().receipt),
            )
    }

    pub fn report_json(&self) -> &str {
        &self.report
    }
}

fn text<'a>(definition: &'a J, key: &str, what: &str) -> Result<&'a str, String> {
    definition
        .str_field(key)
        .filter(|value| !value.trim().is_empty())
        .ok_or_else(|| format!("{what}.{key}: a nonempty string is required"))
}

/// Read one pack file of at most `limit` bytes; `field` names the scene field
/// that gave its path. Gas cavities read their source packs through this too.
pub(super) fn read_pack(path: &Path, limit: u64, field: &str) -> Result<Vec<u8>, String> {
    let file = std::fs::File::open(path)
        .map_err(|error| format!("{field}: {}: {error}", path.display()))?;
    let metadata = file
        .metadata()
        .map_err(|error| format!("{field}: {error}"))?;
    if !metadata.is_file() || metadata.len() > limit {
        return Err(format!(
            "{field}: expected a regular file of at most {limit} bytes"
        ));
    }
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("{field}: {error}"))?;
    if bytes.len() as u64 > limit {
        return Err(format!("{field}: input grew beyond its byte limit"));
    }
    Ok(bytes)
}

/// Load the optional named material map before lowering any elements.
/// Pack paths are relative to the scene; path spellings confer no identity.
pub(super) fn load(root: &J, scene_path: &Path) -> Result<Materials, String> {
    let Some(declarations) = root.get("materials") else {
        return Ok(Materials::new());
    };
    let declarations = declarations
        .as_object()
        .filter(|items| items.len() <= MAX_MATERIALS)
        .ok_or_else(|| {
            format!("materials: expected an object of at most {MAX_MATERIALS} entries")
        })?;
    let mut materials = Materials::new();
    let mut total_bytes = 0;
    for (name, declaration) in declarations {
        let what = format!("materials.{name}");
        if name.trim().is_empty() || materials.contains_key(name) {
            return Err(format!(
                "{what}: material names must be nonempty and unique"
            ));
        }
        fields(
            declaration,
            &[
                "pack",
                "phase",
                "temperature_axis",
                "temperature_kind",
                "state",
                "temperature_min_k",
                "temperature_max_k",
                "reference_specific_enthalpy_j_kg",
                "heat_capacity_property",
                "density_property",
                "conductivity_property",
                "conductivity_kind",
                "resistivity_property",
                "selection",
            ],
            &what,
        )?;
        let pack_path = text(declaration, "pack", &what)?;
        if pack_path.len() > MAX_CARD_PACK_SOURCE_BYTES {
            return Err(format!(
                "{what}.pack: path exceeds {MAX_CARD_PACK_SOURCE_BYTES} bytes"
            ));
        }
        let pack_path = scene_path
            .parent()
            .unwrap_or(Path::new("."))
            .join(pack_path);
        let bytes = read_pack(
            &pack_path,
            MAX_CARD_PACK_BYTES.min(MAX_TOTAL_PACK_BYTES - total_bytes),
            &format!("{what}.pack"),
        )?;
        total_bytes += bytes.len() as u64;
        let pack = NormalizedMaterialCardPack::from_bytes(&bytes)
            .map_err(|error| format!("{what}.pack: {error}"))?;
        let mut material = prepare(declaration, &pack, &what)?;
        material.report = report(name, declaration, &pack, &bytes, &material);
        materials.insert(name.clone(), material);
    }
    Ok(materials)
}

/// The declaration's `conductivity_property` under its `conductivity_kind`
/// (dimensional unless named), or none when it names no conductivity.
fn conductivity_key(declaration: &J, what: &str) -> Result<Option<PropertyKey>, String> {
    if declaration.get("conductivity_property").is_none() {
        if declaration.get("conductivity_kind").is_some() {
            return Err(format!(
                "{what}.conductivity_kind: requires conductivity_property"
            ));
        }
        return Ok(None);
    }
    let kind = if declaration.get("conductivity_kind").is_some() {
        text(declaration, "conductivity_kind", what)?
    } else {
        "dimensional"
    };
    let conductivity_quantity = quantity(kind, THERMAL_CONDUCTIVITY_DIMS)
        .map_err(|error| format!("{what}.conductivity_kind: {error}"))?;
    let static_conductivity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    let Some(conductivity_quantity) = conductivity_quantity.filter(|q| {
        *q == QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS) || *q == static_conductivity
    }) else {
        return Err(format!(
            "{what}.conductivity_kind: expected dimensional or thermal-conductivity"
        ));
    };
    Ok(Some(PropertyKey::with_quantity(
        text(declaration, "conductivity_property", what)?,
        conductivity_quantity,
    )))
}

fn prepare(
    declaration: &J,
    pack: &NormalizedMaterialCardPack,
    what: &str,
) -> Result<PreparedMaterial, String> {
    let phase = match text(declaration, "phase", what)? {
        "solid" => SolidLiquidPhase::Solid,
        "liquid" => SolidLiquidPhase::Liquid,
        _ => {
            return Err(format!(
                "{what}.phase: sourced calorics require solid or liquid"
            ));
        }
    };
    let axis = text(declaration, "temperature_axis", what)?;
    let temperature_kind = quantity(
        text(declaration, "temperature_kind", what)?,
        Temperature::DIMS,
    )
    .map_err(|error| format!("{what}.temperature_kind: {error}"))?;
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if temperature_kind
        .is_some_and(|q| q != absolute && q != QuantitySpec::dimensional(Temperature::DIMS))
    {
        return Err(format!(
            "{what}.temperature_kind: an absolute temperature convention is required"
        ));
    }
    let context = state_point(declaration).map_err(|error| format!("{what}.state: {error}"))?;
    if context.axes().contains_key(axis) {
        return Err(format!(
            "{what}.state: the temperature axis {axis:?} is declared separately"
        ));
    }
    let lower_k = positive(declaration, "temperature_min_k", what)?;
    let upper_k = positive(declaration, "temperature_max_k", what)?;
    if lower_k >= upper_k {
        return Err(format!(
            "{what}: temperature_min_k must be less than temperature_max_k"
        ));
    }
    let point = |temperature| -> Result<QueryPoint, String> {
        at_temperature(&context, axis, temperature_kind, temperature)
            .map_err(|error| format!("{what}.state: {error}"))
    };
    let lower = point(lower_k)?;
    let upper = point(upper_k)?;
    let reference_h = number(declaration, "reference_specific_enthalpy_j_kg", what)?;
    if !reference_h.is_finite() {
        return Err(format!(
            "{what}.reference_specific_enthalpy_j_kg: must be finite"
        ));
    }
    let selection = match declaration.get("selection") {
        None => SelectionPolicy::SingleClaimOnly,
        Some(value) => SelectionPolicy::from_tag(
            value
                .as_str()
                .ok_or_else(|| format!("{what}.selection: expected a policy name"))?,
        )
        .map_err(|error| format!("{what}.selection: {error}"))?,
    };
    let heat_capacity = PropertyKey::new(
        text(declaration, "heat_capacity_property", what)?,
        SPECIFIC_HEAT_CAPACITY_DIMS,
    );
    let density = PropertyKey::new(text(declaration, "density_property", what)?, Density::DIMS);
    let resolved = ResolvedHeatCapacityCurve::try_new_exact(
        pack.card(),
        HeatCapacitySourceQuery {
            heat_capacity: &heat_capacity,
            density: &density,
            temperature_axis: axis,
            lower: &lower,
            upper: &upper,
            selection,
        },
        phase,
        reference_h,
    )
    .map_err(|error| format!("{what}.calorics: {error}"))?;
    let material_point = ThermalMaterialPoint {
        material: pack.card(),
        context: &lower,
        temperature_axis: axis,
        liquid_fraction_axis: None,
    };
    // A continuous source line is not by itself a static material law. Keep
    // the existing owner's exclusions for hidden frequency, independently
    // varying temperatures and undeclared dynamical source coordinates.
    let requirements = [&heat_capacity, &density]
        .into_iter()
        .map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
        })
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| format!("{what}.calorics: {error}"))?;
    let material_selection = match selection {
        SelectionPolicy::SingleClaimOnly => MaterialPropertySelection::SingleClaimOnly,
        SelectionPolicy::PreferObservationBacked => {
            MaterialPropertySelection::PreferObservationBacked
        }
    };
    let properties = resolve_material_state_point(
        pack.card(),
        &lower,
        &requirements,
        material_selection.clone(),
    )
    .map_err(|error| format!("{what}.calorics: {error}"))?;
    for property in properties.properties() {
        material_point
            .admit_static(property)
            .map_err(|error| format!("{what}.calorics: {error}"))?;
    }
    let conductivity_key = conductivity_key(declaration, what)?;
    let (conductivity, conductivity_source) = match &conductivity_key {
        Some(key) => {
            let layer = ResolvedConductionLayer::new(
                resolved.curve(),
                material_point,
                key,
                ClaimSelection::Policy(selection),
                phase,
                [lower_k, upper_k],
                || Ok(()),
            )
            .map_err(|error| format!("{what}.conductivity_property: {error}"))?;
            (
                Some(layer.conductivity_curve().clone()),
                Some(Arc::clone(layer.source_line())),
            )
        }
        None => (None, None),
    };
    // A resistivity is resolved where the calorics were, at the span's lower
    // end; a resistor's node binds its geometry and reads the same source
    // line across the span.
    let resistivity = if declaration.get("resistivity_property").is_some() {
        let basis = ElectricalPropertyBasis::Resistivity;
        let key = PropertyKey::new(
            text(declaration, "resistivity_property", what)?,
            basis.dimensions(),
        );
        Some(
            resolve_conductor_state_with_property(
                pack.card(),
                &lower,
                material_selection,
                basis,
                &key,
            )
            .map_err(|error| format!("{what}.resistivity_property: {error}"))?,
        )
    } else {
        None
    };
    Ok(PreparedMaterial {
        caloric: resolved.curve().clone(),
        conductivity,
        conductivity_source,
        report: String::new(),
        source: CardSource {
            card: pack.card().clone(),
            heat_capacity,
            density,
            conductivity: conductivity_key,
            resistivity,
            axis: axis.to_owned(),
            context,
            temperature_kind,
            lower,
            upper,
            span_k: [lower_k, upper_k],
            selection,
            phase,
            lower_specific_enthalpy_j_kg: reference_h,
        },
    })
}

fn uncertainty_json(uncertainty: &UncertaintyModel) -> String {
    match *uncertainty {
        UncertaintyModel::Unstated => "{\"kind\":\"unstated\"}".into(),
        UncertaintyModel::HalfWidth {
            half_width,
            confidence,
        } => format!(
            "{{\"kind\":\"absolute-half-width\",\"half_width\":{},\"confidence\":{}}}",
            num(half_width),
            num(confidence),
        ),
        UncertaintyModel::RelativeHalfWidth {
            fraction,
            confidence,
        } => format!(
            "{{\"kind\":\"relative-half-width\",\"fraction\":{},\"confidence\":{}}}",
            num(fraction),
            num(confidence),
        ),
    }
}

fn report(
    name: &str,
    declaration: &J,
    pack: &NormalizedMaterialCardPack,
    bytes: &[u8],
    material: &PreparedMaterial,
) -> String {
    let mut receipts = String::from("[");
    let mut claims = BTreeSet::new();
    for (index, receipt) in material.receipts().enumerate() {
        if index > 0 {
            receipts.push(',');
        }
        let coordinates = receipt
            .query_point
            .iter()
            .map(|(axis, value)| format!("{}:{}", quote(axis), num(*value)))
            .collect::<Vec<_>>()
            .join(",");
        let quantities = receipt
            .axis_quantities
            .iter()
            .map(|(axis, quantity)| format!("{}:{:?}", quote(axis), quantity.canonical_bytes()))
            .collect::<Vec<_>>()
            .join(",");
        let _ = write!(
            receipts,
            "{{\"receipt_hash\":{},\"property\":{},\"selected_claim\":{},\"policy\":{},\"observation_backed\":{},\"query_point\":{{{coordinates}}},\"axis_quantities\":{{{quantities}}}}}",
            quote(&receipt.content_hash().to_hex()),
            quote(&receipt.property),
            quote(&receipt.selected.0.to_hex()),
            quote(receipt.policy),
            receipt.observation_backed,
        );
        claims.insert(receipt.selected);
    }
    receipts.push(']');
    let claims = claims.into_iter().filter_map(|id| {
        let claim = pack.card().claims().claim(id)?;
        Some(format!(
            "{{\"claim_hash\":{},\"property\":{},\"source\":{},\"license\":{},\"uncertainty\":{}}}",
            quote(&id.0.to_hex()), quote(claim.key.name()), quote(&claim.provenance.source),
            quote(&claim.provenance.license), uncertainty_json(&claim.uncertainty),
        ))
    }).collect::<Vec<_>>().join(",");
    let bounds = material.caloric.knots();
    format!(
        "{{\"name\":{},\"pack_hash\":{},\"pack_bytes_hash\":{},\"card_hash\":{},\"card_state\":{},\"caloric_hash\":{},\"phase\":{},\"temperature_axis\":{},\"temperature_kind\":{},\"temperature_min_k\":{},\"temperature_max_k\":{},\"reference_specific_enthalpy_j_kg\":{},\"has_conductivity\":{},\"uncertainty_propagated\":false,\"claims\":[{claims}],\"source_receipts\":{receipts}}}",
        quote(name),
        quote(&pack.content_hash().to_hex()),
        quote(&hash_bytes(bytes).to_hex()),
        quote(&pack.card().content_hash().to_hex()),
        quote(&pack.card().id().to_string()),
        quote(&material.caloric.identity().to_hex()),
        quote(declaration.str_field("phase").unwrap_or("")),
        quote(declaration.str_field("temperature_axis").unwrap_or("")),
        quote(declaration.str_field("temperature_kind").unwrap_or("")),
        num(bounds[0].temperature_k),
        num(bounds[bounds.len() - 1].temperature_k),
        num(bounds[0].specific_enthalpy_j_kg),
        material.conductivity.is_some(),
    )
}
