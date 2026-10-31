//! Explicit material-data ingress for the existing bounded transient experiment.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use super::{
    HIGH_K, HeatCapacitySourceQuery, LOW_K, MaterialCard, MaterialEnthalpySource,
    MaterialEnthalpyTransport, MaterialStateId, NormalizedMaterialCardPack, PropertyKey,
    QuantitySpec, QueryPoint, Result, SelectionPolicy, SolidLiquidPhase, read_pack,
};

pub const USAGE: &str = "material_transient <nasa316|aa2024|nasapeek|source> <pack-or-descriptor> <thickness_m> <cells> <duration_s> <steps> [<profile2|source> <pack-or-descriptor2> <thickness2_m> <contact_R_m2_K_per_W>]
The source descriptor starts with frankensim.material-transient.source.v1, then:
pack = path/to/material.fsmatpk
phase = solid
temperature-axis = temperature
density = density-property-name
heat-capacity = specific-heat-property-name
conductivity = scalar-property-name
axis source-condition-name = value_in_SI
For orthotropic transport use three conductivity names and add
principal-axes = e1x e1y e1z e2x e2y e2z e3x e3y e3z
Rows are finite orthonormal principal directions in the mesh frame. Phase is
explicitly solid or liquid; neither selects or models a phase transition.
All six named fields are required. Unknown/duplicate fields or axes refuse.
Pack paths are relative to the descriptor. Axis types and property units come
from the original normalized pack; all numerical coordinates use coherent SI.
No source condition, missing property or extrapolated interval is inferred.
Limits: 64 KiB/descriptor, 128 lines, 64 axes, 16 MiB/pack. The experiment retains
300 K initial/reference and 360/300 K reservoirs, a 20x20 mm section, 1000 W/(m2 K)
end films and insulated sides; every source must support the full 300–360 K path.";

const HEADER: &str = "frankensim.material-transient.source.v1";
const MAX_BYTES: usize = 64 * 1024;

#[derive(Debug)]
struct SourceInput {
    pack: PathBuf,
    phase: SolidLiquidPhase,
    temperature_axis: String,
    density: String,
    heat_capacity: String,
    conductivity: Vec<String>,
    principal_axes: Option<[[f64; 3]; 3]>,
    axes: BTreeMap<String, f64>,
}

fn number(text: &str) -> Result<f64> {
    let value: f64 = text.parse()?;
    if !value.is_finite() {
        return Err("source descriptor numbers must be finite".into());
    }
    Ok(value)
}

fn parse(text: &str, base: &Path) -> Result<SourceInput> {
    if text.len() > MAX_BYTES || text.lines().count() > 128 {
        return Err("source descriptor exceeds 64 KiB or 128 lines".into());
    }
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'));
    if lines.next() != Some(HEADER) {
        return Err(format!("source descriptor must start with {HEADER}").into());
    }
    let mut fields = BTreeMap::new();
    let mut axes = BTreeMap::new();
    for line in lines {
        let (name, value) = line
            .split_once('=')
            .ok_or("source descriptor needs name = value")?;
        let (name, value) = (name.trim(), value.trim());
        if value.is_empty() {
            return Err(format!("empty source descriptor field {name}").into());
        }
        if let Some(axis) = name.strip_prefix("axis ").map(str::trim) {
            if axis.is_empty()
                || axes.len() == 64
                || axes.insert(axis.to_owned(), number(value)?).is_some()
            {
                return Err("source axes need at most 64 distinct, nonempty names".into());
            }
        } else if !matches!(
            name,
            "pack"
                | "phase"
                | "temperature-axis"
                | "density"
                | "heat-capacity"
                | "conductivity"
                | "principal-axes"
        ) {
            return Err(format!("unknown source descriptor field {name}").into());
        } else if fields.insert(name, value).is_some() {
            return Err(format!("duplicate source descriptor field {name}").into());
        }
    }
    let required = |name| {
        fields
            .get(name)
            .copied()
            .ok_or_else(|| format!("missing source descriptor field {name}"))
    };
    let phase = match required("phase")? {
        "solid" => SolidLiquidPhase::Solid,
        "liquid" => SolidLiquidPhase::Liquid,
        _ => return Err("source phase must explicitly be solid or liquid".into()),
    };
    let temperature_axis = required("temperature-axis")?.to_owned();
    if axes.contains_key(&temperature_axis) {
        return Err("the transient experiment owns its temperature axis".into());
    }
    let conductivity: Vec<_> = required("conductivity")?
        .split_whitespace()
        .map(str::to_owned)
        .collect();
    let principal_axes = fields
        .get("principal-axes")
        .map(|text| -> Result<_> {
            let components: Vec<_> = text.split_whitespace().map(number).collect::<Result<_>>()?;
            if components.len() != 9 {
                return Err("principal-axes requires nine row-major components".into());
            }
            Ok(std::array::from_fn(|i| {
                std::array::from_fn(|j| components[3 * i + j])
            }))
        })
        .transpose()?;
    if !matches!(
        (conductivity.len(), principal_axes.is_some()),
        (1, false) | (3, true)
    ) {
        return Err(
            "select one conductivity without a frame, or three conductivities with principal-axes"
                .into(),
        );
    }
    Ok(SourceInput {
        pack: base.join(required("pack")?),
        phase,
        temperature_axis,
        density: required("density")?.into(),
        heat_capacity: required("heat-capacity")?.into(),
        conductivity,
        principal_axes,
        axes,
    })
}

// Reuse the original source keys and coordinate descriptors. A name alone
// cannot silently select between dimensionally/semantically different keys.
fn source_line(
    card: &MaterialCard,
    input: &SourceInput,
) -> Result<(Vec<PropertyKey>, [QueryPoint; 2])> {
    let mut keys = Vec::new();
    let mut types: BTreeMap<String, Option<QuantitySpec>> = BTreeMap::new();
    for name in [&input.density, &input.heat_capacity]
        .into_iter()
        .chain(&input.conductivity)
    {
        let claims = card.claims().claims_for(name);
        let (_, first) = claims
            .first()
            .ok_or_else(|| format!("missing material property {name}"))?;
        if claims.iter().any(|(_, claim)| claim.key != first.key) {
            return Err(
                format!("material property {name} has ambiguous complete source keys").into(),
            );
        }
        keys.push(first.key.clone());
        for (_, claim) in claims {
            for axis in claim.validity.bounds().keys() {
                let quantity = claim.validity.axis_quantities().get(axis).copied();
                if types
                    .insert(axis.clone(), quantity)
                    .is_some_and(|previous| previous != quantity)
                {
                    return Err(format!("source axis {axis} has conflicting quantity types").into());
                }
            }
        }
    }
    let point = |temperature| -> Result<QueryPoint> {
        let mut point = QueryPoint::new();
        for (name, value) in input
            .axes
            .iter()
            .map(|(name, &value)| (name.as_str(), value))
            .chain(std::iter::once((
                input.temperature_axis.as_str(),
                temperature,
            )))
        {
            let quantity = types.get(name).ok_or_else(|| {
                format!("axis {name} is not declared by the selected material properties")
            })?;
            point = match quantity {
                Some(quantity) => point.with_quantity(name, *quantity, value)?,
                None => point.with(name, value)?,
            };
        }
        Ok(point)
    };
    Ok((keys, [point(LOW_K)?, point(HIGH_K)?]))
}

fn resolve(card: &MaterialCard, input: &SourceInput) -> Result<MaterialEnthalpyTransport> {
    let (keys, [lower, upper]) = source_line(card, input)?;
    let source = MaterialEnthalpySource {
        caloric: HeatCapacitySourceQuery {
            density: &keys[0],
            heat_capacity: &keys[1],
            temperature_axis: &input.temperature_axis,
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        conductivity: &keys[2],
        reference: &lower,
        phase: input.phase,
        lower_specific_enthalpy_j_kg: 0.0,
        maximum_interpolation_error_j_kg: 0.01,
    };
    Ok(match input.principal_axes {
        None => MaterialEnthalpyTransport::from_isotropic_card(card, source)?,
        Some(axes) => MaterialEnthalpyTransport::from_orthotropic_card(
            card,
            MaterialEnthalpySource {
                conductivity: [&keys[2], &keys[3], &keys[4]],
                caloric: source.caloric,
                reference: source.reference,
                phase: source.phase,
                lower_specific_enthalpy_j_kg: source.lower_specific_enthalpy_j_kg,
                maximum_interpolation_error_j_kg: source.maximum_interpolation_error_j_kg,
            },
            axes,
        )?,
    })
}

pub fn load(path: &Path) -> Result<(NormalizedMaterialCardPack, MaterialEnthalpyTransport)> {
    let mut text = String::new();
    std::fs::File::open(path)?
        .take(MAX_BYTES as u64 + 1)
        .read_to_string(&mut text)?;
    let input = parse(&text, path.parent().unwrap_or_else(|| Path::new(".")))?;
    let pack = read_pack(&input.pack)?;
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: match input.phase {
                SolidLiquidPhase::Solid => "solid",
                _ => "liquid",
            }
            .into(),
            process:
                "source-selected stationary thermal specimen; original source conditions retained"
                    .into(),
            revision: 0,
        },
        pack,
    )?;
    let material = resolve(pack.card(), &input)?;
    Ok((pack, material))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, PropertyClaim, PropertyValue, Provenance, UncertaintyModel,
    };
    use fs_qty::{Density, Dims, Pressure, Temperature};

    const INPUT: &str = "frankensim.material-transient.source.v1\npack = material.fsmatpk\nphase = solid\ntemperature-axis = temp\ndensity = rho\nheat-capacity = cp\nconductivity = k\naxis pressure = 101325\n";

    #[test]
    fn parses_relative_pack_and_explicit_scalar_or_principal_sources() {
        let input = parse(INPUT, Path::new("case")).unwrap();
        assert_eq!(input.pack, Path::new("case/material.fsmatpk"));
        assert_eq!(input.axes["pressure"], 101325.0);
        assert!(input.principal_axes.is_none());
        let directional = INPUT.replace(
            "conductivity = k",
            "conductivity = k1 k2 k3\nprincipal-axes = 0 1 0 -1 0 0 0 0 1",
        );
        let input = parse(&directional, Path::new("case")).unwrap();
        assert_eq!(input.conductivity, ["k1", "k2", "k3"]);
        assert_eq!(input.principal_axes.unwrap()[0], [0.0, 1.0, 0.0]);
        parse(include_str!("nasa316.source"), Path::new(".")).unwrap();
    }

    #[test]
    fn refuses_ambiguous_incomplete_and_oversized_descriptions() {
        for tail in [
            "\npack = other",
            "\naxis pressure = 2",
            "\nunknown = 1",
            "\naxis temp = 325",
            "\nprincipal-axes = 1 0 0 0 1 0 0 0 1",
            "\naxis bad = NaN",
        ] {
            assert!(parse(&format!("{INPUT}{tail}"), Path::new(".")).is_err());
        }
        assert!(
            parse(&INPUT.replace("density = rho\n", ""), Path::new("."))
                .unwrap_err()
                .to_string()
                .contains("missing source descriptor field density")
        );
        assert!(
            parse(
                &INPUT.replace("conductivity = k", "conductivity = a b c"),
                Path::new(".")
            )
            .is_err()
        );
        assert!(
            parse(
                &INPUT.replace("phase = solid", "phase = melting"),
                Path::new(".")
            )
            .is_err()
        );
        assert!(parse(&"x".repeat(MAX_BYTES + 1), Path::new(".")).is_err());
    }

    fn card() -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("rho", Density::DIMS, 1000.0),
            ("cp", Dims([2, 0, -2, -1, 0, 0]), 1000.0),
            ("k", Dims([1, 1, -3, -1, 0, 0]), 2.0),
            ("k-longitudinal", Dims([1, 1, -3, -1, 0, 0]), 20.0),
            ("k-through", Dims([1, 1, -3, -1, 0, 0]), 1.0),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Curve {
                        abscissa: "temp".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (360.0, value)],
                        dims,
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "temp",
                            QuantitySpec::dimensional(Temperature::DIMS),
                            300.0,
                            360.0,
                        )
                        .with_quantity(
                            "pressure",
                            QuantitySpec::dimensional(Pressure::DIMS),
                            100000.0,
                            102000.0,
                        ),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: "synthetic numerical fixture, not measured data".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "fixture".into(),
                phase: "solid".into(),
                process: "numerical test".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn original_source_quantities_and_conditions_reach_the_shared_thermal_owner() {
        let card = card();
        let input = parse(INPUT, Path::new(".")).unwrap();
        let material = resolve(&card, &input).unwrap();
        assert_eq!(material.reference_density_kg_m3(), 1000.0);
        assert_eq!(
            material.reference_query().axis_quantities()["pressure"],
            QuantitySpec::dimensional(Pressure::DIMS)
        );
        for receipt in material.conductivity().receipts() {
            card.claims().verify_receipt(receipt).unwrap();
        }
        let turned = INPUT.replace(
            "conductivity = k",
            "conductivity = k-longitudinal k k-through\nprincipal-axes = 0 1 0 -1 0 0 0 0 1",
        );
        let directional = resolve(&card, &parse(&turned, Path::new(".")).unwrap()).unwrap();
        assert_eq!(
            directional.conductivity().tensor_at(330.0).unwrap(),
            [[2.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(
            directional.caloric().curve().identity(),
            material.caloric().curve().identity()
        );
        for bad in [
            INPUT.replace("axis pressure = 101325\n", ""),
            INPUT.replace("101325", "200000"),
            INPUT.replace("conductivity = k", "conductivity = missing"),
            format!("{INPUT}axis not-a-source-axis = 0\n"),
        ] {
            assert!(resolve(&card, &parse(&bad, Path::new(".")).unwrap()).is_err());
        }
    }

    #[test]
    fn descriptor_file_and_original_pack_decode_preserve_data_selected_material() {
        let directory = std::env::temp_dir().join(format!(
            "frankensim-material-transient-source-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
        ));
        std::fs::create_dir(&directory).unwrap();
        let pack = super::super::NormalizedPack::new(
            "synthetic-transient-source-fixture",
            "native source-input test",
            fs_blake3::hash_domain(
                "org.frankensim.tests.material-transient.source.v1",
                b"synthetic numerical fixture",
            ),
            "CC0-1.0; synthetic numerical test data",
            card().claims().clone(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        std::fs::write(directory.join("material.fsmatpk"), pack.to_bytes()).unwrap();
        for (name, description, expected) in [
            ("scalar.source", INPUT.to_owned(), [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]),
            ("turned.source", INPUT.replace("conductivity = k",
                "conductivity = k-longitudinal k k-through\nprincipal-axes = 0 1 0 -1 0 0 0 0 1"),
                [[2.0,0.0,0.0],[0.0,20.0,0.0],[0.0,0.0,1.0]]),
        ] {
            let path = directory.join(name);
            std::fs::write(&path, description).unwrap();
            let (loaded, material) = load(&path).unwrap();
            assert_eq!(material.card_identity(), loaded.card().content_hash());
            assert_eq!(material.reference_density_kg_m3(), 1000.0);
            assert_eq!(material.conductivity().tensor_at(330.0).unwrap(), expected);
            for receipt in material.caloric().receipts().iter().chain(material.conductivity().receipts()) {
                loaded.card().claims().verify_receipt(receipt).unwrap();
            }
        }
    }
}
