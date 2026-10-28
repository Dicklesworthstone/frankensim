//! A source-backed steady wall assembled from the committed CFAST reference packs.
//!
//! First compile the three manifests with xtask matdb-pack, then pass their
//! output files in gypsum / glass-fiber insulation / common-brick order.
//! The source numbers are nominal CFAST model inputs with unknown test conditions.
//! The film/contact resistances below are declared scenario assumptions.

#![deny(unsafe_code)]

use std::error::Error;
use std::path::Path;

use fs_conduction::layered::{LayeredWall, PlanarThermalLayer};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, PropertyKey,
    QueryPoint, SelectionPolicy,
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
    let pack = NormalizedPack::from_bytes(&std::fs::read(path)?)?;
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

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args_os().skip(1).collect::<Vec<_>>();
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
