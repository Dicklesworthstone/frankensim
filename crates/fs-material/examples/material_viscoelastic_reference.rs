//! Run a source-bound incremental Prony specimen under mixed strain/force control.
//!
//! ```text
//! material_viscoelastic_reference PACK LENGTH_M AREA_M2 MIN_INCREMENTAL_STRAIN MAX_INCREMENTAL_STRAIN MAX_ABS_STRAIN_RATE_S MAX_ELAPSED_S
//!   --time-domain-reference [--axis NAME VALUE_SI]...
//!   [--source-length-axis NAME] [--source-area-axis NAME]
//!   [--strain VALUE DURATION_S | --force VALUE_N DURATION_S]...
//! ```
//! Endpoints are deviations from a fixed source operating state. They are not
//! additions to the previous endpoint. The initially relaxed incremental state
//! retains its branch history across positive-duration ramps, holds and control
//! switches. Geometry and history bounds are authored explicitly.
//!
//! Time-domain use of fitted spectra is a model reference, not independent
//! experimental validation. This small incremental uniaxial model supplies no
//! finite-strain crushing, permanent set, damage or temperature evolution.

#![deny(unsafe_code)]

use fs_matdb::{
    MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, QueryPoint,
    SelectionPolicy,
};
use fs_material::visco::source::{
    PronyReferenceLimits, PronyReferenceRequest, ResolvedPronyReference,
};
use fs_qty::QuantitySpec;
use std::{collections::BTreeMap, error::Error, io::Read};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const USAGE: &str = "material_viscoelastic_reference PACK LENGTH_M AREA_M2 MIN_INCREMENTAL_STRAIN MAX_INCREMENTAL_STRAIN MAX_ABS_STRAIN_RATE_S MAX_ELAPSED_S --time-domain-reference [--axis NAME VALUE_SI]... [--source-length-axis NAME] [--source-area-axis NAME] [--strain VALUE DURATION_S | --force VALUE_N DURATION_S]...";

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

fn spectral_property(name: &str) -> bool {
    matches!(
        name,
        "density" | "prony-branch-count" | "prony-equilibrium-modulus"
    ) || ["prony-modulus-", "prony-relaxation-time-"]
        .iter()
        .any(|prefix| {
            name.strip_prefix(*prefix)
                .is_some_and(|suffix| suffix.parse::<usize>().is_ok_and(|branch| branch > 0))
        })
}

fn source_point(card: &MaterialCard, coordinates: &BTreeMap<String, f64>) -> Result<QueryPoint> {
    let mut types: BTreeMap<String, Option<QuantitySpec>> = BTreeMap::new();
    for (_, claim) in card.claims().claims_ordered() {
        if !spectral_property(claim.key.name()) {
            continue;
        }
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
    let mut point = QueryPoint::new();
    for (name, value) in coordinates {
        point = match types.get(name).ok_or_else(|| {
            invalid(format!(
                "axis {name} is not declared by the source coefficient properties"
            ))
        })? {
            Some(quantity) => point.with_quantity(name, *quantity, *value)?,
            None => point.with(name, *value)?,
        };
    }
    Ok(point)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(1025).collect();
    if args.as_slice() == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() < 7 || args.len() > 1024 || args.iter().any(|s| s.len() > 4096) {
        return Err(invalid(USAGE));
    }
    let (length, area) = (number(&args[1])?, number(&args[2])?);
    let mut request = PronyReferenceRequest {
        limits: PronyReferenceLimits {
            minimum_incremental_strain: number(&args[3])?,
            maximum_incremental_strain: number(&args[4])?,
            maximum_absolute_strain_rate_s: number(&args[5])?,
            maximum_elapsed_s: number(&args[6])?,
        },
        original_length_axis: None,
        original_area_axis: None,
    };
    let mut axes = BTreeMap::new();
    let mut time_domain_reference = false;
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
            "--source-length-axis" | "--source-area-axis" if i + 1 < args.len() => {
                let axis = if args[i] == "--source-length-axis" {
                    &mut request.original_length_axis
                } else {
                    &mut request.original_area_axis
                };
                if axis.replace(args[i + 1].clone()).is_some() {
                    return Err(invalid("duplicate source-geometry axis"));
                }
                i += 2;
            }
            "--time-domain-reference" => {
                if time_domain_reference {
                    return Err(invalid("duplicate time-domain reference choice"));
                }
                time_domain_reference = true;
                i += 1;
            }
            "--strain" | "--force" if i + 2 < args.len() => {
                let duration = number(&args[i + 2])?;
                if duration <= 0.0 {
                    return Err(invalid(
                        "each ramp or hold needs a positive duration; zero-time jumps are unsupported",
                    ));
                }
                steps.push((
                    args[i].trim_start_matches("--").to_owned(),
                    number(&args[i + 1])?,
                    duration,
                ));
                i += 3;
            }
            _ => {
                return Err(invalid(format!(
                    "unknown or incomplete argument: {}",
                    args[i]
                )));
            }
        }
        if axes.len() > 64 || steps.len() > 256 {
            return Err(invalid(
                "at most 64 source axes and 256 history steps are supported",
            ));
        }
    }
    if !time_domain_reference {
        return Err(invalid(
            "--time-domain-reference is required; fitted spectra do not establish independent transient validation",
        ));
    }
    if steps.is_empty() {
        return Err(invalid(
            "at least one --strain or --force history step is required",
        ));
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
            process: "source-fitted incremental Prony time-domain reference".into(),
            revision: 0,
        },
        pack,
    )?;
    let point = source_point(pack.card(), &axes)?;
    let material = ResolvedPronyReference::resolve(
        pack.card(),
        &point,
        &request,
        SelectionPolicy::SingleClaimOnly,
    )?;
    // An axis attached only to an unselected competing claim cannot make the
    // selected coefficients appear to support that condition.
    for axis in point.axes().keys() {
        if !material.properties().iter().any(|(_, answer)| {
            pack.card()
                .claims()
                .claim(answer.receipt.selected)
                .is_some_and(|claim| claim.validity.bounds().contains_key(axis))
        }) {
            return Err(invalid(format!(
                "source axis {axis} is not owned by any selected coefficient claim"
            )));
        }
    }
    let mut coupon = material.coupon(length, area)?;
    let mut responses = Vec::with_capacity(steps.len());
    for (control, target, duration) in steps {
        let response = if control == "strain" {
            coupon.advance_strain(target, duration)?
        } else {
            coupon.advance_force(target, duration)?
        };
        responses.push((control, target, duration, response));
    }
    // Validate the complete requested path before emitting a successful result.
    println!("frankensim.viscoelastic-reference.v1");
    println!(
        "scope\tsource-fitted incremental linear model; no independent transient validation, finite-strain crushing, permanent set, damage or thermal-evolution claim"
    );
    println!(
        "reference-state\trelaxed incremental state about the fixed source operating point; unknown absolute preload stress and energy are excluded"
    );
    println!("time-domain-reference\texplicitly-authored");
    println!("pack\t{pack_id}");
    println!("pack-identity\t{pack_identity}");
    println!("card-identity\t{}", material.card_identity());
    println!("reference-law-identity\t{}", material.identity());
    println!("original-length-m\t{length:.17e}");
    println!("original-area-m2\t{area:.17e}");
    match coupon.mass_kg() {
        Some(mass) => println!("mass-kg\t{mass:.17e}"),
        None => println!("mass-kg\tunstated"),
    }
    println!(
        "authored-strain-limits\t{:.17e}\t{:.17e}",
        request.limits.minimum_incremental_strain, request.limits.maximum_incremental_strain
    );
    println!(
        "authored-max-absolute-strain-rate-per-s\t{:.17e}",
        request.limits.maximum_absolute_strain_rate_s
    );
    println!(
        "authored-max-elapsed-s\t{:.17e}",
        request.limits.maximum_elapsed_s
    );
    println!(
        "source-length-axis\t{}",
        request
            .original_length_axis
            .as_deref()
            .unwrap_or("unbound; geometry authored")
    );
    println!(
        "source-area-axis\t{}",
        request
            .original_area_axis
            .as_deref()
            .unwrap_or("unbound; geometry authored")
    );
    for (name, value) in point.axes() {
        println!(
            "source-axis\t{name}\t{value:.17e}\t{:?}",
            point.axis_quantities().get(name)
        );
    }
    for (name, answer) in material.properties() {
        println!(
            "source-property\t{name}\t{:.17e}\t{}\t{}\t{:?}\t{:?}",
            answer.evidence.value.value,
            answer.receipt.selected.0,
            answer.receipt.content_hash(),
            answer.receipt.decision,
            answer.evidence.value.uncertainty
        );
    }
    println!(
        "equilibrium-modulus-pa\t{:.17e}",
        material.equilibrium_modulus_pa()
    );
    for (index, (modulus, seconds)) in material.branches().iter().enumerate() {
        println!(
            "maxwell-branch\t{}\t{modulus:.17e}\t{seconds:.17e}",
            index + 1
        );
    }
    println!(
        "step\tcontrol\ttarget\tduration-s\telapsed-s\tincremental-strain\tincremental-displacement-m\tincremental-stress-pa\tincremental-force-n\tstored-energy-j\tmodel-work-j\tdissipated-energy-j\tincremental-model-work-j\tincremental-dissipated-energy-j\tenergy-balance-residual-j\tresponse-identity"
    );
    for (index, (control, target, duration, response)) in responses.iter().enumerate() {
        println!(
            "{}\t{control}\t{target:.17e}\t{duration:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}\t{}",
            index + 1,
            response.elapsed_s,
            response.incremental_strain,
            response.incremental_displacement_m,
            response.incremental_stress_pa,
            response.incremental_force_n,
            response.stored_energy_j,
            response.model_work_j,
            response.dissipated_energy_j,
            response.incremental_model_work_j,
            response.incremental_dissipated_energy_j,
            response.energy_balance_residual_j,
            response.identity
        );
    }
    Ok(())
}
