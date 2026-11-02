//! G0/G3: actual mineral/water packs feed Bond throughput and slurry accounting.
//! The blend, work index and dilution target are authored scenarios, not ore tests.
#![deny(unsafe_code)]

use std::{path::Path, process::Command};

use fs_matdb::{MatDbError, NormalizedPack, Provenance, QueryPoint, SelectionPolicy};
use fs_material::{
    comminution::{BondGrindingCard, BondSizeDomain, BondTestContext, BondWorkIndex},
    slurry::{DensitySource, SlurryComponent, SlurryError, SlurryStream, SolidsMassFraction},
};
use fs_qty::{
    Dimensionless, Dims, Length, MassFlowRate, Power, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn compile(folder: &str) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "frankensim-mining-slurry-{}-{nonce}-{folder}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{folder}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{folder}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    pack
}

fn mineral_point() -> QueryPoint {
    [
        ("mineral-reference-model", 1.0),
        ("source-test-temperature-known", 0.0),
        ("source-test-pressure-known", 0.0),
        ("bulk-ore-qualified", 0.0),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |point, (axis, value)| {
        point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

fn water_point(kelvin: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap()
        .with_quantity("phase-liquid", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
}

fn source(pack: &NormalizedPack, property: &str, point: &QueryPoint) -> DensitySource {
    let source =
        DensitySource::resolve(pack, property, point, SelectionPolicy::SingleClaimOnly).unwrap();
    assert_eq!(source.pack_id(), pack.pack_id());
    assert_eq!(source.pack_hash(), pack.content_hash());
    assert_eq!(source.query(), point);
    assert_eq!(source.answer().receipt.property, property);
    assert!(source.answer().receipt.observation_backed);
    pack.claims()
        .verify_receipt(&source.answer().receipt)
        .unwrap();
    source
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g3_bond_feed_dilution_and_split_recombine_retain_real_sources() {
    let feed_size = Length::new(0.0025);
    let product_size = Length::new(0.0001);
    let grinding = BondGrindingCard::new(
        BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap(),
        BondTestContext {
            specimen: "authored mineral blend; not a measured ore".into(),
            method: "authored Bond ball-mill scenario with prepared feed".into(),
            closing_screen: Length::new(150e-6),
            provenance: Provenance {
                source: "authored integration fixture, not a site calibration".into(),
                license: "project test fixture".into(),
                artifact: None,
            },
        },
        BondSizeDomain::new([feed_size; 2], [product_size; 2]).unwrap(),
    )
    .unwrap();
    close(
        grinding
            .specific_energy(feed_size, product_size)
            .unwrap()
            .value(),
        43_200.0,
    );
    let dry = grinding
        .dry_throughput(
            feed_size,
            product_size,
            Power::new(12e6),
            Dimensionless::new(0.9),
        )
        .unwrap();
    // 12 MW * 0.9 / 43200 J/kg = 250 kg/s, independently of dilution.
    close(dry.value(), 250.0);
    let mut solids = Vec::new();
    for (mineral, fraction, density) in [
        ("magnetite", 0.5, 5175.0),
        ("pyrite", 0.2, 5018.0),
        ("dolomite", 0.3, 2860.0),
    ] {
        let pack = compile(&format!("{mineral}-handbook-mineralogy-reference"));
        let selected = source(&pack, "density", &mineral_point());
        close(selected.density().value(), density);
        solids.push(
            SlurryComponent::new(mineral, MassFlowRate::new(dry.value() * fraction), selected)
                .unwrap(),
        );
    }
    let water_pack = compile("water-liquid-iapws-sr6-08");
    assert_eq!(water_pack.pack_id(), "water-liquid-iapws-sr6-08-0.1mpa");
    let water = source(&water_pack, "density", &water_point(298.15));
    close(water.density().value(), 997.047013399765);
    let initial = SlurryStream::new_additive_volume(
        solids,
        SlurryComponent::new("process-water", MassFlowRate::new(50.0), water).unwrap(),
    )
    .unwrap();
    let diluted = initial
        .dilute_to_mass_fraction(SolidsMassFraction::new(Dimensionless::new(0.65)).unwrap())
        .unwrap();
    assert_eq!(diluted.stream.solids(), initial.solids());
    assert_eq!(
        diluted.stream.liquid().density_source(),
        initial.liquid().density_source()
    );

    // Independent Decimal50 reference using the cited source densities:
    // Vs=125/5175+50/5018+75/2860; mw=250*(1-.65)/.65; Vw=mw/rhow.
    let state = diluted.stream.state();
    close(diluted.added_liquid_flow.value(), 84.61538461538462);
    close(state.dry_solids_mass_flow.value(), 250.0);
    close(state.liquid_mass_flow.value(), 134.61538461538462);
    close(state.total_mass_flow.value(), 384.61538461538464);
    close(state.solid_volume_flow.value(), 0.06034249473087049);
    close(state.liquid_volume_flow.value(), 0.13501407938264463);
    close(state.total_volume_flow.value(), 0.19535657411351512);
    close(state.density.unwrap().value(), 1968.786493931336);
    close(state.solids_mass_fraction.unwrap().value(), 0.65);
    close(
        state.solids_volume_fraction.unwrap().value(),
        0.3088838704542776,
    );

    let (first, second) = diluted.stream.split(Dimensionless::new(0.37)).unwrap();
    for (outlet, share) in [(&first, 0.37), (&second, 0.63)] {
        close(
            outlet.state().total_mass_flow.value(),
            state.total_mass_flow.value() * share,
        );
        close(
            outlet.state().total_volume_flow.value(),
            state.total_volume_flow.value() * share,
        );
        close(
            outlet.state().density.unwrap().value(),
            state.density.unwrap().value(),
        );
        close(outlet.state().solids_mass_fraction.unwrap().value(), 0.65);
    }
    let restored = SlurryStream::mix(&[&first, &second]).unwrap();
    assert_eq!(restored, SlurryStream::mix(&[&second, &first]).unwrap());
    for (original, restored) in diluted.stream.solids().iter().zip(restored.solids()) {
        assert_eq!(original.id(), restored.id());
        close(original.mass_flow().value(), restored.mass_flow().value());
        assert_eq!(original.density_source(), restored.density_source());
    }
    close(
        restored.liquid().mass_flow().value(),
        state.liquid_mass_flow.value(),
    );
    assert_eq!(
        restored.liquid().density_source(),
        initial.liquid().density_source()
    );
    close(
        restored.state().density.unwrap().value(),
        state.density.unwrap().value(),
    );
}

#[test]
fn g3_range_endpoints_and_liquid_conditions_never_become_silent_defaults() {
    let mineral = compile("chalcopyrite-handbook-mineralogy-reference");
    assert!(matches!(
        DensitySource::resolve(
            &mineral,
            "density",
            &mineral_point(),
            SelectionPolicy::SingleClaimOnly
        ),
        Err(SlurryError::Material(MatDbError::UnknownProperty { .. }))
    ));
    let water_pack = compile("water-liquid-iapws-sr6-08");
    let water = source(&water_pack, "density", &water_point(298.15));
    let mut alternatives = Vec::new();
    for (property, expected) in [
        ("density-reported-range-minimum", 4100.0),
        ("density-reported-range-maximum", 4300.0),
    ] {
        let selected = source(&mineral, property, &mineral_point());
        close(selected.density().value(), expected);
        let stream = SlurryStream::new_additive_volume(
            vec![SlurryComponent::new("chalcopyrite", MassFlowRate::new(1.0), selected).unwrap()],
            SlurryComponent::new("process-water", MassFlowRate::new(1.0), water.clone()).unwrap(),
        )
        .unwrap();
        close(stream.state().solid_volume_flow.value(), 1.0 / expected);
        alternatives.push(stream);
    }
    assert!(matches!(
        SlurryStream::mix(&[&alternatives[0], &alternatives[1]]),
        Err(SlurryError::SourceMismatch { component }) if component == "chalcopyrite"
    ));
    for wrong in [
        QueryPoint::new(),
        water_point(370.0),
        water_point(298.15).with("pressure", 101_325.0).unwrap(),
        water_point(298.15).with("phase-liquid", 0.0).unwrap(),
    ] {
        assert!(matches!(
            DensitySource::resolve(
                &water_pack,
                "density",
                &wrong,
                SelectionPolicy::SingleClaimOnly
            ),
            Err(SlurryError::Material(_))
        ));
    }
    // A second supported temperature is still a different source query; mixing
    // cannot silently average density or claim thermal equilibration.
    let warmer = source(&water_pack, "density", &water_point(303.15));
    let warm_stream = SlurryStream::new_additive_volume(
        alternatives[0].solids().to_vec(),
        SlurryComponent::new("process-water", MassFlowRate::new(1.0), warmer).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        SlurryStream::mix(&[&alternatives[0], &warm_stream]),
        Err(SlurryError::SourceMismatch { component }) if component == "process-water"
    ));
}
