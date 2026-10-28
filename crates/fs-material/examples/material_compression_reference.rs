//! Run a bounded compression reference from original stress/strain observations.
//!
//! Arguments:
//! ```text
//! material_compression_reference PACK DENSITY_KEY STRESS_KEY STRAIN_AXIS LENGTH_M AREA_M2
//!   [--axis NAME VALUE_SI]... [--speed-axis NAME] [--interpolate-observations]
//!   [--step SHORTENING_M DURATION_S]...
//! ```
//! Shortening, compressive strain, stress and force are positive magnitudes.
//! The initial state is the FIRST supported source strain. Loading work starts
//! there: the program never invents the missing path from an unloaded origin.
//! A known source crosshead speed must match every ramp. With unknown source
//! speed, durations remain authored choices and establish no rate dependence.
//!
//! By default the source must permit interpolation. The explicit
//! `--interpolate-observations` flag authors a piecewise-linear reference
//! between discrete observations without changing their original receipts.
//! This reference supplies no unloading, recovery, creep, hysteresis, stored
//! energy, heat generation or failure prediction. Geometry is caller supplied;
//! source specimen geometry and lateral contraction are not inferred.
//!
//! For the glass-filled PEEK pack, select `compressive-stress` against
//! `compressive-strain`, its grade/unknown-condition markers and the exact
//! `test-speed`. The 1%, 2% and 5% source observations remain distinct from
//! intermediate values supplied only by the chosen reference model.

#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, QueryPoint, SelectionPolicy,
};
use fs_material::compression_reference::{
    CompressionReferenceInterpolation, CompressionReferenceRequest, ResolvedCompressionReference,
};
use fs_qty::QuantitySpec;
use std::collections::BTreeMap;
use std::error::Error;
use std::io::Read;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const USAGE: &str = "material_compression_reference PACK DENSITY_KEY STRESS_KEY STRAIN_AXIS LENGTH_M AREA_M2 [--axis NAME VALUE_SI]... [--speed-axis NAME] [--interpolate-observations] [--step SHORTENING_M DURATION_S]...";
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
fn property_key(card: &MaterialCard, name: &str) -> Result<PropertyKey> {
    let claims = card.claims().claims_for(name);
    let (_, first) = claims
        .first()
        .ok_or_else(|| invalid(format!("missing source property {name}")))?;
    if claims.iter().any(|(_, claim)| claim.key != first.key) {
        return Err(invalid(format!(
            "property {name} has ambiguous complete source keys"
        )));
    }
    Ok(first.key.clone())
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
    if args.len() < 6 || args.len() > 512 || args.iter().any(|s| s.len() > 4096) {
        return Err(invalid(USAGE));
    }
    let (length, area) = (number(&args[4])?, number(&args[5])?);
    let mut axes = BTreeMap::new();
    let mut speed_axis = None;
    let mut interpolate_observations = false;
    let mut steps = Vec::new();
    let mut i = 6;
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
            "--interpolate-observations" => {
                if interpolate_observations {
                    return Err(invalid("duplicate interpolation choice"));
                }
                interpolate_observations = true;
                i += 1;
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
            process: "bounded compression reference; original source context retained".into(),
            revision: 0,
        },
        pack,
    )?;
    let point = source_point(pack.card(), &[&args[1], &args[2]], &axes)?;
    let request = CompressionReferenceRequest {
        density: property_key(pack.card(), &args[1])?,
        stress: property_key(pack.card(), &args[2])?,
        strain_axis: args[3].clone(),
        crosshead_speed_axis: speed_axis,
        interpolation: if interpolate_observations {
            CompressionReferenceInterpolation::PiecewiseLinearBetweenObservations
        } else {
            CompressionReferenceInterpolation::SourceInterpolationOnly
        },
    };
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let material =
        ResolvedCompressionReference::resolve(pack.card(), &point, &request, selection, selection)?;
    let mut coupon = material.coupon(length, area)?;
    let mut responses = Vec::with_capacity(steps.len());
    for (shortening, duration) in steps {
        responses.push(coupon.advance_shortening(shortening, duration)?);
    }
    // Validate the whole requested path before emitting any successful result.
    println!("frankensim.compression-reference.v1");
    println!(
        "scope\tbounded loading reference; no unloading, creep, recovery, stored-energy or failure claim"
    );
    println!("interpolation-choice\t{:?}", request.interpolation);
    println!("pack\t{pack_id}");
    println!("pack-identity\t{pack_identity}");
    println!("card-identity\t{}", material.card_identity());
    println!("reference-law-identity\t{}", material.identity());
    println!("reference-length-m\t{length:.17e}");
    println!("reference-area-m2\t{area:.17e}");
    let [lower, upper] = material.strain_domain();
    println!("source-strain-domain\t{lower:.17e}\t{upper:.17e}");
    println!(
        "initial-shortening-m\t{:.17e}",
        coupon.initial_shortening_m()
    );
    println!("work-reference\tzero loading work at the first supported source strain");
    println!("mass-kg\t{:.17e}", coupon.mass_kg());
    println!(
        "source-crosshead-speed-axis\t{}",
        request
            .crosshead_speed_axis
            .as_deref()
            .unwrap_or("unstated; ramp duration authored")
    );
    for (name, value) in point.axes() {
        println!(
            "source-axis\t{name}\t{value:.17e}\t{:?}",
            point.axis_quantities().get(name)
        );
    }
    let density = material.density_answer();
    println!(
        "source-density\t{:.17e}\t{}\t{:?}",
        density.evidence.value.value,
        density.receipt.selected.0,
        density.evidence.value.uncertainty
    );
    for (strain, answer) in material.observations() {
        println!(
            "source-stress-point\t{strain:.17e}\t{:.17e}\t{}\t{:?}\t{:?}",
            answer.evidence.value.value,
            answer.receipt.selected.0,
            answer.receipt.decision,
            answer.evidence.value.uncertainty
        );
    }
    println!(
        "step\telapsed-s\tshortening-m\tcompressive-strain\tnominal-stress-pa\tforce-n\ttangent-n-per-m\tmodel-loading-work-j\tincremental-model-loading-work-j\tresponse-identity"
    );
    for (index, response) in responses.iter().enumerate() {
        println!(
            "{}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{}",
            index + 1,
            response.elapsed_s,
            response.shortening_m,
            response.compressive_strain,
            response.nominal_stress_pa,
            response.force_n,
            response.tangent_n_m,
            response.model_loading_work_j,
            response.incremental_model_loading_work_j,
            response.identity
        );
    }
    Ok(())
}
