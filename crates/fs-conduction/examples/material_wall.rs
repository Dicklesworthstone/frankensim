//! A source-backed steady wall assembled from the committed CFAST reference packs.
//!
//! First compile the three manifests with xtask matdb-pack, then pass their
//! output files in gypsum / glass-fiber insulation / common-brick order.
//! The source numbers are nominal CFAST model inputs with unknown test conditions.
//! The film/contact resistances below are declared scenario assumptions.
//!
//! A general source-bound nonlinear mode also accepts explicit layer packs,
//! geometry, reservoirs, films, contacts and source axes. Run `--help` for syntax.
//! Axis values use SI units; their types come from the source declarations.
//! No engineering-model selector is supplied implicitly. Every layer must
//! support the complete reservoir interval; no extrapolation is attempted.

#![deny(unsafe_code)]

use std::collections::BTreeMap;
use std::error::Error;
use std::ffi::OsString;
use std::io::Read;
use std::path::Path;

use fs_conduction::layered::{LayeredWall, NonlinearWallOptions, PlanarThermalLayer};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, QueryPoint, SelectionPolicy,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const MATERIALS: [(&str, f64, &str); 3] = [
    (
        "gypsum-board-5-8-cfast",
        0.016,
        "source-product-and-formulation-known",
    ),
    (
        "glass-fiber-insulation-cfast",
        0.088,
        "source-product-and-install-density-known",
    ),
    (
        "brick-common-cfast-sp1041",
        0.076,
        "source-material-processing-known",
    ),
];

fn source_point(temperature_k: f64, source_condition: &str) -> QueryPoint {
    let mut point = QueryPoint::new();
    for axis in [
        "source-test-temperature-known",
        "source-moisture-condition-known",
        "source-measurement-method-known",
        source_condition,
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 0.0)
            .expect("literal dimensionless source-context coordinate");
    }
    point
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature_k,
        )
        .expect("literal positive absolute temperature")
}

fn load(path: &Path, slug: &str) -> Result<NormalizedMaterialCardPack, Box<dyn Error>> {
    let pack = read_pack(path)?;
    if pack.pack_id() != slug {
        return Err(format!("expected {slug}, received {}", pack.pack_id()).into());
    }
    Ok(NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process: "CFAST nominal layer-input reference; unspecified test conditions".into(),
            revision: 0,
        },
        pack,
    )?)
}

fn cfast(args: &[OsString]) -> Result<(), Box<dyn Error>> {
    if args.len() != 3 {
        return Err(
            "usage: material_wall <gypsum.fsmatpk> <glass-fiber.fsmatpk> <brick.fsmatpk>".into(),
        );
    }
    let property = PropertyKey::with_quantity(
        "thermal-conductivity",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    );
    let mut layers = Vec::new();
    for (argument, (slug, thickness_m, condition)) in args.iter().zip(MATERIALS) {
        let pack = load(Path::new(argument), slug)?;
        let layer = PlanarThermalLayer::from_card(
            slug,
            thickness_m,
            pack.card(),
            &property,
            "temperature",
            &[
                source_point(273.15, condition),
                source_point(293.15, condition),
            ],
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )?;
        println!(
            "material={slug} card={} thickness_m={thickness_m}",
            pack.card().content_hash()
        );
        layers.push(layer);
    }
    // The values here are explicit example assumptions, not source claims.
    let wall = LayeredWall::new(layers, vec![0.01, 0.002], [0.13, 0.04])?;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let solution = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 51,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        wall.solve(&cx, [293.15, 273.15], 10.0)
    })?;
    println!("Nominal steady planar model; CFAST source test temperature/moisture remain unknown.");
    println!(
        "R_m2_K_per_W={} U_W_per_m2_K={} flux_W_per_m2={} heat_W={}",
        solution.total_resistance_m2_k_per_w,
        solution.thermal_transmittance_w_per_m2_k,
        solution.heat_flux_w_per_m2,
        solution.heat_rate_w
    );
    for (layer, result) in solution
        .wall
        .layers()
        .iter()
        .zip(&solution.layer_temperatures)
    {
        println!(
            "layer={} left_K={} right_K={} Fourier_residual_W_per_m2={} source_receipts={}",
            layer.name(),
            result.left_k,
            result.right_k,
            result.flux_residual_w_per_m2,
            layer.conductivity().receipts().len()
        );
    }
    println!(
        "maximum_Fourier_residual_W_per_m2={} right_boundary_residual_K={}",
        solution.maximum_flux_residual_w_per_m2, solution.right_boundary_residual_k
    );
    Ok(())
}

type CliResult<T> = Result<T, Box<dyn Error>>;
const USAGE: &str = "material_wall <gypsum.fsmatpk> <glass-fiber.fsmatpk> <brick.fsmatpk>\nmaterial_wall nonlinear --left-k K --right-k K --area-m2 A [--films Rleft Rright] --layer PACK THICKNESS_M [--axis NAME=VALUE_SI]... [--contact R --layer PACK THICKNESS_M ...]\nAxes apply to the preceding layer. Contacts and films use m2 K/W. Omitted films prescribe face temperatures. Supply each contact explicitly, including zero. Limits: 16 layers, 64 axes/layer, 512 arguments, 128 KiB argument bytes, 16 MiB/pack. Reservoir temperatures must differ.";

fn invalid(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::new(std::io::ErrorKind::InvalidInput, message.into()).into()
}

fn number(text: &str) -> CliResult<f64> {
    let value = text.parse::<f64>()?;
    if !value.is_finite() {
        return Err(invalid("all numerical arguments must be finite"));
    }
    Ok(value)
}

fn read_pack(path: &Path) -> CliResult<NormalizedPack> {
    const LIMIT: u64 = 16 * 1024 * 1024;
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take(LIMIT + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 > LIMIT {
        return Err(invalid("source pack exceeds 16 MiB"));
    }
    Ok(NormalizedPack::from_bytes(&bytes)?)
}

struct LayerInput {
    path: String,
    thickness_m: f64,
    axes: BTreeMap<String, f64>,
}

fn source_line(
    card: &MaterialCard,
    coordinates: &BTreeMap<String, f64>,
    temperatures: [f64; 2],
) -> CliResult<(PropertyKey, [QueryPoint; 2])> {
    let claims = card.claims().claims_for("thermal-conductivity");
    let (_, first) = claims
        .first()
        .ok_or_else(|| invalid("missing thermal-conductivity"))?;
    if claims.iter().any(|(_, claim)| claim.key != first.key) {
        return Err(invalid(
            "thermal-conductivity has ambiguous complete source keys",
        ));
    }
    let mut types: BTreeMap<String, Option<QuantitySpec>> = BTreeMap::new();
    for (_, claim) in &claims {
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
    let point = |temperature| -> CliResult<QueryPoint> {
        let mut point = QueryPoint::new();
        for (name, value) in coordinates
            .iter()
            .map(|(k, v)| (k.as_str(), *v))
            .chain(std::iter::once(("temperature", temperature)))
        {
            let quantity = types.get(name).ok_or_else(|| {
                invalid(format!(
                    "axis {name} is not declared by thermal-conductivity"
                ))
            })?;
            point = match quantity {
                Some(quantity) => point.with_quantity(name, *quantity, value)?,
                None => point.with(name, value)?,
            };
        }
        Ok(point)
    };
    Ok((
        first.key.clone(),
        [point(temperatures[0])?, point(temperatures[1])?],
    ))
}

fn nonlinear(args: &[String]) -> CliResult<()> {
    let (mut left, mut right, mut area, mut films) = (None, None, None, None);
    let mut inputs: Vec<LayerInput> = Vec::new();
    let mut contacts = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--left-k" | "--right-k" | "--area-m2" if i + 1 < args.len() => {
                let slot = match args[i].as_str() {
                    "--left-k" => &mut left,
                    "--right-k" => &mut right,
                    _ => &mut area,
                };
                if slot.replace(number(&args[i + 1])?).is_some() {
                    return Err(invalid("duplicate boundary or area argument"));
                }
                i += 2;
            }
            "--films" if i + 2 < args.len() => {
                if films
                    .replace([number(&args[i + 1])?, number(&args[i + 2])?])
                    .is_some()
                {
                    return Err(invalid("duplicate films"));
                }
                i += 3;
            }
            "--layer" if i + 2 < args.len() => {
                if inputs.len() == 16 || contacts.len() != inputs.len() {
                    return Err(invalid(
                        "at most 16 layers; declare a contact before each later layer",
                    ));
                }
                inputs.push(LayerInput {
                    path: args[i + 1].clone(),
                    thickness_m: number(&args[i + 2])?,
                    axes: BTreeMap::new(),
                });
                i += 3;
            }
            "--axis" if i + 1 < args.len() && contacts.len() + 1 == inputs.len() => {
                let (name, value) = args[i + 1]
                    .split_once('=')
                    .ok_or_else(|| invalid("axis syntax is NAME=VALUE_SI"))?;
                if name.is_empty() || name == "temperature" {
                    return Err(invalid(
                        "source axis needs a name; temperature is owned by the reservoirs",
                    ));
                }
                let axes = &mut inputs.last_mut().expect("layer count checked").axes;
                if axes.len() == 64 || axes.insert(name.into(), number(value)?).is_some() {
                    return Err(invalid("at most 64 distinct source axes per layer"));
                }
                i += 2;
            }
            "--contact" if i + 1 < args.len() && contacts.len() + 1 == inputs.len() => {
                contacts.push(number(&args[i + 1])?);
                i += 2;
            }
            _ => {
                return Err(invalid(format!(
                    "unknown, incomplete or misplaced argument: {}",
                    args[i]
                )));
            }
        }
    }
    let (left, right, area) = match (left, right, area) {
        (Some(left), Some(right), Some(area))
            if left > 0.0 && right > 0.0 && area > 0.0 && left != right =>
        {
            (left, right, area)
        }
        _ => {
            return Err(invalid(
                "declare positive, distinct --left-k and --right-k, and positive --area-m2",
            ));
        }
    };
    if inputs.is_empty() || contacts.len() + 1 != inputs.len() {
        return Err(invalid(
            "require at least one layer and exactly one explicit contact per adjacent pair",
        ));
    }
    let mut layers = Vec::new();
    let mut identities = Vec::new();
    for (index, input) in inputs.iter().enumerate() {
        let pack = read_pack(Path::new(&input.path))?;
        let pack_id = pack.pack_id().to_owned();
        identities.push((pack_id.clone(), pack.content_hash()));
        let pack = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: pack_id,
                phase: "solid".into(),
                process: "steady planar conductivity reference; explicit source context".into(),
                revision: 0,
            },
            pack,
        )?;
        let (property, points) =
            source_line(pack.card(), &input.axes, [left.min(right), left.max(right)])?;
        layers.push(PlanarThermalLayer::from_card(
            format!("layer-{}", index + 1),
            input.thickness_m,
            pack.card(),
            &property,
            "temperature",
            &points,
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )?);
    }
    let wall = LayeredWall::new(layers, contacts, films.unwrap_or([0.0; 2]))?;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let solution = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 51,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        wall.solve_temperature_dependent(
            &cx,
            [left, right],
            area,
            NonlinearWallOptions {
                max_iterations: 128,
                temperature_tolerance_k: 1.0e-8,
                flux_absolute_tolerance_w_per_m2: 1.0e-8,
                flux_relative_tolerance: 1.0e-10,
            },
        )
    })?;
    // Emit a result only after every source query and the complete solve succeed.
    println!("frankensim.nonlinear-material-wall.v1");
    println!(
        "scope=nominal steady planar reference; residuals are numerical diagnostics, not material uncertainty"
    );
    println!(
        "left_K={left:.17e} right_K={right:.17e} area_m2={area:.17e} films_m2_K_per_W={:?} contacts_m2_K_per_W={:?}",
        wall.surface_resistances_m2_k_per_w(),
        wall.contact_resistances_m2_k_per_w()
    );
    println!(
        "flux_W_per_m2={:.17e} heat_W={:.17e} effective_R_m2_K_per_W={:.17e} effective_U_W_per_m2_K={:.17e}",
        solution.heat_flux_w_per_m2,
        solution.heat_rate_w,
        solution.effective_resistance_m2_k_per_w,
        solution.effective_transmittance_w_per_m2_k
    );
    for ((layer, result), (pack_id, pack_hash)) in wall
        .layers()
        .iter()
        .zip(&solution.layer_temperatures)
        .zip(&identities)
    {
        println!(
            "layer={} pack={} pack_identity={} card_identity={} thickness_m={:.17e} left_K={:.17e} right_K={:.17e} Fourier_residual_W_per_m2={:.17e}",
            layer.name(),
            pack_id,
            pack_hash,
            layer.material_card().expect("card-bound layer"),
            layer.thickness_m(),
            result.left_k,
            result.right_k,
            result.flux_residual_w_per_m2
        );
        for ((temperature, conductivity), receipt) in layer
            .conductivity()
            .knots()
            .iter()
            .zip(layer.conductivity().receipts())
        {
            println!(
                "source_receipt layer={} temperature_K={temperature:.17e} conductivity_W_per_m_K={conductivity:.17e} identity={} claim={} policy={} decision={:?} point={:?} axis_quantities={:?} source_hashes={:?}",
                layer.name(),
                receipt.content_hash(),
                receipt.selected.0,
                receipt.policy,
                receipt.decision,
                receipt.query_point,
                receipt.axis_quantities,
                receipt.source_hashes
            );
        }
    }
    println!(
        "maximum_Fourier_residual_W_per_m2={:.17e} cold_boundary_residual_K={:.17e} iterations={} flux_bracket_W_per_m2={:?} numerical_limits={:?}",
        solution.maximum_flux_residual_w_per_m2,
        solution.cold_boundary_residual_k,
        solution.iterations,
        solution.flux_bracket_w_per_m2,
        solution.options
    );
    Ok(())
}

fn main() -> CliResult<()> {
    let args: Vec<_> = std::env::args_os().skip(1).take(513).collect();
    if args.len() > 512
        || args.iter().any(|s| s.len() > 4096)
        || args.iter().map(|s| s.len()).sum::<usize>() > 128 * 1024
    {
        return Err(invalid("argument size limit exceeded"));
    }
    if args.len() == 1 && args[0] == "--help" {
        println!("{USAGE}");
        return Ok(());
    }
    if args.first().is_some_and(|s| s == "nonlinear") {
        let args = args[1..]
            .iter()
            .map(|s| {
                s.to_str()
                    .map(str::to_owned)
                    .ok_or_else(|| invalid("nonlinear arguments must be UTF-8"))
            })
            .collect::<CliResult<Vec<_>>>()?;
        nonlinear(&args)
    } else {
        cfast(&args)
    }
}
