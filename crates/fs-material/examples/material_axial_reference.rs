//! Run a source-bound, monotonic axial MODULUS REFERENCE.
//!
//! All material values come from the supplied canonical source pack. The
//! caller supplies geometry, loading direction, source conditions, durations
//! and a reference strain ceiling. The zero-origin linear reference and its
//! work are authored calculations, not measured stress curves, recoverable
//! energy, actual relaxation, yield/failure predictions, or 3-D elasticity.
//!
//! Arguments:
//! ```text
//! material_axial_reference PACK tension|compression DENSITY_KEY MODULUS_KEY LENGTH_M AREA_M2 STRAIN_CEILING
//!   [--axis NAME VALUE_SI]... [--speed-axis NAME] [--step DISPLACEMENT_M DURATION_S]...
//! ```
//! At most 64 source axes and 64 strictly monotonic steps are accepted. Axis
//! quantity types come from the selected source claims. A source speed axis
//! must be explicitly owned and its m/s value must match each imposed ramp.
//! With unknown source test speed, the ramp duration remains an authored
//! scenario choice and makes no source-rate or relaxation claim.
//!
//! PEEK tension, using the new ensinger-tecapeek-natural pack:
//! ```text
//! material_axial_reference peek.matdb-pack tension density tensile-modulus 0.05 0.0001 0.001
//!   --axis source-grade-tecapeek-natural 1 --axis source-temperature-known 0
//!   --axis source-moisture-condition-known 0 --axis source-pressure-known 0
//!   --axis test-speed 1.6666666666666667e-5 --speed-axis test-speed
//!   --step 0.00005 3
//! ```
//! The authored 50 mm / 100 mm2 specimen gives nominal force 420 N, mass
//! 0.00655 kg and model loading work 0.0105 J from the source's 4.2 GPa and
//! 1310 kg/m3. The source's test specimen geometry is not inferred.
//!
//! ROHACELL 51 IG-F compression, with unknown test rate preserved:
//! ```text
//! material_axial_reference rohacell.matdb-pack compression density compression-modulus 0.05 0.0001 0.001
//!   --axis source-grade-rohacell-51-igf-2022 1 --axis source-temperature-known 0
//!   --axis source-moisture-condition-known 0 --axis source-test-rate-known 0
//!   --axis full-sheet-at-least-10mm 1 --step -0.00005 1
//! ```
//! This reference gives -4.3 N, 0.00026 kg and 0.0001075 J from 43 MPa and
//! 52 kg/m3. Its one-second duration is authored; foam recovery is unknown.
//! Other Ensinger grades use their own grade marker and exact source protocol.
//! PA66 additionally requires immediately-after-machining=1.

#![deny(unsafe_code)]

use fs_matdb::{
    MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, QueryPoint,
};
use fs_material::axial_reference::{
    AxialLoadingMode, AxialReferenceRequest, ResolvedAxialModulusReference,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::QuantitySpec;
use std::collections::BTreeMap;
use std::error::Error;
use std::io::Read;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const USAGE: &str = "material_axial_reference PACK tension|compression DENSITY_KEY MODULUS_KEY LENGTH_M AREA_M2 STRAIN_CEILING [--axis NAME VALUE_SI]... [--speed-axis NAME] [--step DISPLACEMENT_M DURATION_S]...";
fn invalid(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::new(std::io::ErrorKind::InvalidInput, message.into()).into()
}
fn number(text: &str) -> Result<f64> {
    let value = text.parse::<f64>()?;
    if !value.is_finite() {
        return Err(invalid("all numerical arguments must be finite"));
    }
    Ok(value)
}
fn requirement(card: &MaterialCard, name: &str) -> Result<ScalarPropertyRequirement> {
    let claims = card.claims().claims_for(name);
    let (_, first) = claims
        .first()
        .ok_or_else(|| invalid(format!("missing source property {name}")))?;
    if claims.iter().any(|(_, claim)| claim.key != first.key) {
        return Err(invalid(format!(
            "property {name} has ambiguous complete source keys"
        )));
    }
    Ok(ScalarPropertyRequirement::try_with_key(
        &first.key,
        ScalarAdmissibility::StrictlyPositive,
    )?)
}
fn source_point(
    card: &MaterialCard,
    names: &[&str],
    coordinates: &BTreeMap<String, f64>,
) -> Result<QueryPoint> {
    let mut types: BTreeMap<String, Option<QuantitySpec>> = BTreeMap::new();
    for name in names {
        for (_, claim) in card.claims().claims_for(name) {
            for axis in claim.validity.bounds().keys() {
                let quantity = claim.validity.axis_quantities().get(axis).copied();
                if let Some(previous) = types.insert(axis.clone(), quantity) {
                    if previous != quantity {
                        return Err(invalid(format!(
                            "source axis {axis} has conflicting quantity types"
                        )));
                    }
                }
            }
        }
    }
    let mut point = QueryPoint::new();
    for (name, value) in coordinates {
        let quantity = types.get(name).ok_or_else(|| {
            invalid(format!(
                "axis {name} is not declared by the selected source properties"
            ))
        })?;
        point = match quantity {
            Some(quantity) => point.with_quantity(name, *quantity, *value)?,
            None => point.with(name, *value)?,
        };
    }
    Ok(point)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.as_slice() == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() < 7 || args.len() > 512 || args.iter().any(|s| s.len() > 4096) {
        return Err(invalid(USAGE));
    }
    let mode = match args[1].as_str() {
        "tension" => AxialLoadingMode::Tension,
        "compression" => AxialLoadingMode::Compression,
        _ => return Err(invalid("select tension or compression explicitly")),
    };
    let (length, area, limit) = (number(&args[4])?, number(&args[5])?, number(&args[6])?);
    let mut axes = BTreeMap::new();
    let mut speed_axis = None;
    let mut steps = Vec::new();
    let mut i = 7;
    while i < args.len() {
        match args[i].as_str() {
            "--axis" if i + 2 < args.len() => {
                if axes
                    .insert(args[i + 1].clone(), number(&args[i + 2])?)
                    .is_some()
                {
                    return Err(invalid("duplicate source axis"));
                }
                i += 3;
            }
            "--speed-axis" if i + 1 < args.len() => {
                if speed_axis.replace(args[i + 1].clone()).is_some() {
                    return Err(invalid("duplicate source-speed owner"));
                }
                i += 2;
            }
            "--step" if i + 2 < args.len() => {
                steps.push((number(&args[i + 1])?, number(&args[i + 2])?));
                i += 3;
            }
            _ => {
                return Err(invalid(format!(
                    "unknown or incomplete argument: {}",
                    args[i]
                )));
            }
        }
        if axes.len() > 64 || steps.len() > 64 {
            return Err(invalid(
                "at most 64 axes and 64 loading steps are supported",
            ));
        }
    }
    if steps.is_empty() {
        return Err(invalid("at least one explicit loading step is required"));
    }
    let mut bytes = Vec::new();
    std::fs::File::open(&args[0])?
        .take(16 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() > 16 * 1024 * 1024 {
        return Err(invalid("source pack exceeds 16 MiB"));
    }
    let pack = NormalizedPack::from_bytes(&bytes)?;
    let pack_id = pack.pack_id().to_owned();
    let pack_identity = pack.content_hash();
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack_id.clone(),
            phase: "solid".into(),
            process: "explicit reported-modulus axial reference; source context retained".into(),
            revision: 0,
        },
        pack,
    )?;
    let point = source_point(pack.card(), &[&args[2], &args[3]], &axes)?;
    let request = AxialReferenceRequest {
        density: requirement(pack.card(), &args[2])?,
        modulus: requirement(pack.card(), &args[3])?,
        mode,
        maximum_absolute_strain: limit,
        crosshead_speed_axis: speed_axis,
    };
    let material = ResolvedAxialModulusReference::resolve(
        pack.card(),
        &point,
        &request,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let mut coupon = material.coupon(length, area)?;
    let mut responses = Vec::with_capacity(steps.len());
    for (displacement, duration) in steps {
        responses.push(coupon.advance_displacement(displacement, duration)?);
    }
    // The whole requested path is checked before any successful result is emitted.
    println!("frankensim.axial-modulus-reference.v1");
    println!(
        "scope\tauthored zero-origin linear reference; no unloading, relaxation, recoverable-energy or failure claim"
    );
    println!("pack\t{pack_id}");
    println!("pack-identity\t{pack_identity}");
    println!("card-identity\t{}", material.source().card_identity());
    println!("source-state-identity\t{}", material.source().identity());
    println!("reference-law-identity\t{}", material.identity());
    println!("loading-mode\t{mode:?}");
    println!("reference-length-m\t{length:.17e}");
    println!("reference-area-m2\t{area:.17e}");
    println!("authored-strain-ceiling\t{limit:.17e}");
    println!(
        "reference-stiffness-n-per-m\t{:.17e}",
        coupon.stiffness_n_m()
    );
    println!("mass-kg\t{:.17e}", coupon.mass_kg());
    println!(
        "source-crosshead-speed-axis\t{}",
        request
            .crosshead_speed_axis
            .as_deref()
            .unwrap_or("unstated; ramp duration authored")
    );
    for (name, value) in material.source().query_point() {
        println!(
            "source-axis\t{name}\t{value:.17e}\t{:?}",
            material.source().axis_quantities().get(name)
        );
    }
    for property in material.source().properties() {
        println!(
            "source-property\t{}\t{:.17e}\t{}\t{:?}",
            property.requirement().name(),
            property.value_si(),
            property.answer().receipt.selected.0,
            property.answer().evidence.value.uncertainty
        );
    }
    println!(
        "step\telapsed-s\tdisplacement-m\tengineering-strain\treference-stress-pa\treference-force-n\tmodel-loading-work-j\tincremental-model-loading-work-j\tresponse-identity"
    );
    for (index, response) in responses.iter().enumerate() {
        println!(
            "{}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{}",
            index + 1,
            response.elapsed_s,
            response.displacement_m,
            response.engineering_strain,
            response.reference_stress_pa,
            response.reference_force_n,
            response.model_loading_work_j,
            response.incremental_model_loading_work_j,
            response.identity
        );
    }
    Ok(())
}
