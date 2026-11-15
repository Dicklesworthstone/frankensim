//! A rigid body declared as a solid: a shape with its dimensions and a
//! material card. The body's mass and principal inertia are the shape's
//! closed-form volume moments (fs-mbd-material `VolumeMoments`) times the
//! density that the card gives at the declared state. They are compiled by
//! the same owner (fs-mbd-material `resolve_body`) that keeps the density's
//! usage receipt. Nothing typed in the scene is a mass or an inertia.
//!
//! A solid with a `temperature_k` is that body freely expanded from the
//! state's temperature: fs-material integrates the card's linear expansion
//! coefficient along the path, and fs-mbd-material `expand_body` scales the
//! shape about the body origin while keeping its mass.

use std::fmt::Write as _;
use std::path::Path;

use fs_blake3::hash_bytes;
use fs_ga::Vec3;
use fs_matdb::{MaterialCard, NormalizedMaterialCardPack, QueryPoint, SelectionPolicy};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, integrate_isotropic_thermal_expansion,
    thermal_path_axis,
};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd_material::thermal::{ThermalExpansionBinding, expand_body};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};

use crate::cards::{MAX_CARD_PACK_BYTES, MAX_CARD_PACK_SOURCE_BYTES};
use crate::material_state_input::state_point;

use super::materials::read_pack;
use super::{Element, J, fields, num, number_or, positive, quote, triple};

const SHAPES: &str = "sphere, cuboid, ellipsoid or cylinder";

/// One body's mass properties as compiled from its shape and card, and the
/// report of what they came from.
pub(super) struct Solid {
    pub mass_kg: f64,
    /// About the body's own axes, which are the shape's principal axes.
    pub principal_inertia_kg_m2: [f64; 3],
    /// The members of a JSON object: shape, card and receipt identities,
    /// density, volume, mass and principal inertia, and for a body expanded
    /// to a `temperature_k` its temperatures, strain, expanded volume and
    /// expansion receipt.
    pub report: String,
}

fn shape_geometry(row: &J, shape: &str, what: &str) -> Result<VolumeMoments, String> {
    let dimensions: &[&str] = match shape {
        "sphere" => &["radius_m"],
        "cuboid" => &["sides_m"],
        "ellipsoid" => &["semiaxes_m"],
        "cylinder" => &["radius_m", "length_m", "inner_radius_m"],
        _ => return Err(format!("{what}.shape: expected {SHAPES}")),
    };
    for key in [
        "radius_m",
        "inner_radius_m",
        "length_m",
        "sides_m",
        "semiaxes_m",
    ] {
        if row.get(key).is_some() && !dimensions.contains(&key) {
            return Err(format!("{what}.{key}: not a dimension of a {shape}"));
        }
    }
    let vector = |key| triple(row, key, None, what).map(|[x, y, z]| Vec3::new(x, y, z));
    match shape {
        "sphere" => VolumeMoments::sphere(positive(row, "radius_m", what)?),
        "cuboid" => VolumeMoments::cuboid(vector("sides_m")?),
        "ellipsoid" => VolumeMoments::ellipsoid(vector("semiaxes_m")?),
        _ => VolumeMoments::hollow_cylinder(
            number_or(row, "inner_radius_m", 0.0, what)?,
            positive(row, "radius_m", what)?,
            positive(row, "length_m", what)?,
        ),
    }
    .map_err(|error| format!("{what}: {error}"))
}

/// The report member that opens what an expanded solid's report adds.
const EXPANDED: &str = ",\"temperature_k\":";

/// Whether a solid's report is that of a body expanded to a `temperature_k`.
pub(super) fn expanded(report: &str) -> bool {
    report.contains(EXPANDED)
}

/// A solid at `temperature_k`: the body freely expanded from the state's
/// temperature about its origin, and the members its report gains. `None`
/// when the solid gives no temperature.
fn expansion(
    row: &J,
    body: &ResolvedBody,
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    what: &str,
) -> Result<Option<(SpatialInertia, String)>, String> {
    if row.get("temperature_k").is_none() {
        return Ok(None);
    }
    let kelvin = positive(row, "temperature_k", what)?;
    let axis = thermal_path_axis(point).ok_or_else(|| {
        format!("{what}.temperature_k: the state has no absolute-temperature axis to expand along")
    })?;
    let hot = point
        .clone()
        .with(axis, kelvin)
        .map_err(|error| format!("{what}.temperature_k: {error}"))?;
    let path = integrate_isotropic_thermal_expansion(card, point, &hot, selection)
        .map_err(|error| format!("{what}.temperature_k: {error}"))?;
    let expanded = expand_body(
        body,
        &[ThermalExpansionBinding {
            key: "solid",
            expansion: &path,
            fixed_point_m: Vec3::new(0.0, 0.0, 0.0),
        }],
        || false,
    )
    .map_err(|error| format!("{what}.temperature_k: {error}"))?;
    let members = format!(
        "{EXPANDED}{},\"reference_temperature_k\":{},\"free_linear_strain\":{},\"expanded_volume_m3\":{},\"expansion_hash\":{},\"expansion_claim\":{}",
        num(kelvin),
        num(point.axes()[axis]),
        num(path.free_linear_strain()),
        num(expanded.regions()[0].geometry().volume_m3()),
        quote(&path.identity().to_hex()),
        quote(&path.selected_claim().0.to_hex()),
    );
    Ok(Some((expanded.inertia(), members)))
}

/// Compile the body of `row.solid`. Pack paths are relative to the scene.
pub(super) fn solid(row: &J, scene: &Path, what: &str) -> Result<Solid, String> {
    fields(
        row,
        &[
            "pack",
            "state",
            "selection",
            "temperature_k",
            "shape",
            "radius_m",
            "inner_radius_m",
            "length_m",
            "sides_m",
            "semiaxes_m",
        ],
        what,
    )?;
    let shape = row
        .str_field("shape")
        .ok_or_else(|| format!("{what}.shape: expected {SHAPES}"))?;
    let geometry = shape_geometry(row, shape, what)?;
    let pack_path = row
        .str_field("pack")
        .filter(|path| !path.trim().is_empty() && path.len() <= MAX_CARD_PACK_SOURCE_BYTES)
        .ok_or_else(|| format!("{what}.pack: a card pack path is required"))?;
    let pack_path = scene.parent().unwrap_or(Path::new(".")).join(pack_path);
    let bytes = read_pack(&pack_path, MAX_CARD_PACK_BYTES, &format!("{what}.pack"))?;
    let pack = NormalizedMaterialCardPack::from_bytes(&bytes)
        .map_err(|error| format!("{what}.pack: {error}"))?;
    let point = state_point(row).map_err(|error| format!("{what}.state: {error}"))?;
    let selection = match row.get("selection") {
        None => MaterialPropertySelection::SingleClaimOnly,
        Some(tag) => {
            match tag
                .as_str()
                .map(SelectionPolicy::from_tag)
                .ok_or_else(|| format!("{what}.selection: expected a policy name"))?
                .map_err(|error| format!("{what}.selection: {error}"))?
            {
                SelectionPolicy::SingleClaimOnly => MaterialPropertySelection::SingleClaimOnly,
                SelectionPolicy::PreferObservationBacked => {
                    MaterialPropertySelection::PreferObservationBacked
                }
            }
        }
    };
    let body = resolve_body(
        &[RegionBinding {
            key: "solid",
            geometry,
            material: pack.card(),
            point: &point,
            selection: selection.clone(),
        }],
        1,
        || false,
    )
    .map_err(|error| format!("{what}: {error}"))?;
    let hot = expansion(row, &body, pack.card(), &point, selection, what)?;
    let inertia = hot
        .as_ref()
        .map_or_else(|| body.inertia(), |(inertia, _)| *inertia);
    let tensor = inertia.inertia_com().m;
    // Every shape here is centred on the body origin with its principal axes
    // along the body's axes; anything else is not this body.
    if [1, 2, 3, 5, 6, 7].iter().any(|&i| tensor[i] != 0.0)
        || inertia.center_of_mass() != Vec3::new(0.0, 0.0, 0.0)
    {
        return Err(format!(
            "{what}: the shape is not centred on its principal axes"
        ));
    }
    let density = body.regions()[0]
        .material()
        .property(DENSITY_PROPERTY)
        .ok_or_else(|| format!("{what}: the card resolved no density"))?;
    let receipt = &density.answer().receipt;
    let principal = [tensor[0], tensor[4], tensor[8]];
    let mut report = format!(
        "\"shape\":{},\"pack_hash\":{},\"pack_bytes_hash\":{},\"card_hash\":{},\"card_state\":{},\"density_kg_m3\":{},\"density_receipt_hash\":{},\"selected_claim\":{},\"volume_m3\":{},\"mass_kg\":{},\"principal_inertia_kg_m2\":[{},{},{}]",
        quote(shape),
        quote(&pack.content_hash().to_hex()),
        quote(&hash_bytes(&bytes).to_hex()),
        quote(&pack.card().content_hash().to_hex()),
        quote(&pack.card().id().to_string()),
        num(density.value_si()),
        quote(&receipt.content_hash().to_hex()),
        quote(&receipt.selected.0.to_hex()),
        num(geometry.volume_m3()),
        num(inertia.mass()),
        num(principal[0]),
        num(principal[1]),
        num(principal[2]),
    );
    if let Some((_, members)) = &hot {
        report.push_str(members);
    }
    Ok(Solid {
        mass_kg: inertia.mass(),
        principal_inertia_kg_m2: principal,
        report,
    })
}

/// The text report's lines for the bodies declared as solids: each one's
/// mass and principal inertia. Empty when a scene has none.
pub(super) fn text_report(elements: &[Element]) -> String {
    let mut out = String::new();
    if elements.iter().all(|element| element.solids.is_empty()) {
        return out;
    }
    out.push_str(
        "\nSolid bodies (mass and inertia from a shape and the density of its card; \
         --json gives the receipts):\n",
    );
    for element in elements {
        for &(index, _) in &element.solids {
            let body = element.bodies[index];
            let inertia = body.principal_inertia_body();
            let _ = writeln!(
                out,
                "{}.body{index}: {} kg, principal inertia [{}, {}, {}] kg m^2",
                element.name,
                num(body.mass()),
                num(inertia.x),
                num(inertia.y),
                num(inertia.z)
            );
        }
    }
    out
}

/// The JSON report's `solids` member, with its leading comma: per element,
/// the report object of each body declared as a solid. Empty when a scene
/// has none.
pub(super) fn json_report(elements: &[Element]) -> String {
    let solids = elements
        .iter()
        .filter(|element| !element.solids.is_empty())
        .map(|element| {
            let reports = element.solids.iter().map(|(_, report)| report.as_str());
            format!(
                "{}:[{}]",
                quote(&element.name),
                reports.collect::<Vec<_>>().join(",")
            )
        })
        .collect::<Vec<_>>();
    if solids.is_empty() {
        String::new()
    } else {
        format!(",\"solids\":{{{}}}", solids.join(","))
    }
}
