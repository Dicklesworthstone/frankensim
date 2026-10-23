//! Render-path material ingress (bead `frankensim-irp-v5-instrument-realism-pexfp.32`,
//! DONE-WHEN 3) on the plate-valve and flat-plate performance files: a section is either a
//! pack-backed material record, resolved through `fs_couple::material_source`
//! from a material-card pack compiled at the named state, or bare numbers with
//! an explicit `fixture LABEL`. Bare numbers without a label refuse, and every
//! section's source reaches the performance's section provenance (the render
//! sidecar's `section_sources`).

#[path = "support/material_packs.rs"]
mod material_packs;

use fs_couple::bernoulli_aperture::performance::file::PlateValvePerformance;
use fs_exec::CancelGate;
use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use material_packs::ALUMINUM_2024_T3;

const INPUT: &str = include_str!("../examples/plate-valve.performance");
const LABELLED: &str =
    "section isotropic 0.0003 900 4000000000 0.3 fixture authored-illustrative-blade-not-measured";

/// The 2024-T3 seed as a material-card pack compiled at state `as-tabulated`.
fn aluminum_card_pack() -> NormalizedMaterialCardPack {
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "aluminum-2024-t3".into(),
            phase: "solid".into(),
            process: "as-tabulated".into(),
            revision: 0,
        },
        NormalizedPack::from_bytes(ALUMINUM_2024_T3).expect("aluminum pack"),
    )
    .expect("card pack")
}

fn with_section(row: &str) -> String {
    assert_eq!(
        INPUT.matches(LABELLED).count(),
        1,
        "the example carries one labelled section"
    );
    INPUT.replace(LABELLED, row)
}

fn load_with_aluminum(text: &str) -> Result<PlateValvePerformance, String> {
    PlateValvePerformance::from_bytes_with_all_sources(
        text.as_bytes(),
        64,
        &CancelGate::new(),
        None,
        |path| {
            if path == "aluminum.fsmatpk" {
                Ok(aluminum_card_pack())
            } else {
                Err(format!("no pack {path}"))
            }
        },
    )
    .map_err(|e| e.to_string())
}

#[test]
fn bare_number_without_fixture_label_refuses() {
    let labelled = PlateValvePerformance::from_bytes(INPUT.as_bytes(), 64, &CancelGate::new())
        .expect("labelled");
    assert_eq!(
        labelled.section_provenance(),
        ["fixture authored-illustrative-blade-not-measured"]
    );
    for row in [
        "section isotropic 0.0003 900 4000000000 0.3",
        "section isotropic 0.0003 900 4000000000 0.3 fixture",
        "section isotropic 0.0003 900 4000000000 0.3 label authored",
        "section isotropic 0.0003 900 4000000000 0.3 fixture a b",
    ] {
        let error = match PlateValvePerformance::from_bytes(
            with_section(row).as_bytes(),
            64,
            &CancelGate::new(),
        ) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("{row}: bare numbers without a fixture label must refuse"),
        };
        assert!(error.contains("fixture LABEL"), "{row}: {error}");
    }
    // A pack-backed section needs a loader: the default constructor has none.
    let pack_row = "section material isotropic 0.00012 aluminum.fsmatpk as-tabulated T=300";
    let error = match PlateValvePerformance::from_bytes(
        with_section(pack_row).as_bytes(),
        64,
        &CancelGate::new(),
    ) {
        Err(e) => e.to_string(),
        Ok(_) => panic!("a pack record without a loader must refuse"),
    };
    assert!(error.contains("material-card loader"), "{error}");
}

#[test]
fn pack_record_matches_direct_resolution_bitwise() {
    let pack = aluminum_card_pack();
    // Resolve the section's three properties directly, at the record's point.
    let temperature = pack
        .card()
        .claims()
        .claims_for("density")
        .first()
        .and_then(|(_, claim)| claim.validity.axis_quantities().get("T").copied())
        .unwrap_or(QuantitySpec::dimensional(fs_qty::Temperature::DIMS));
    let point = QueryPoint::new()
        .with_quantity("T", temperature, 300.0)
        .expect("point");
    let requirement = |name: &str, dims| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite)
            .expect("requirement")
    };
    let state = resolve_material_state_point(
        pack.card(),
        &point,
        &[
            requirement("density", Density::DIMS),
            requirement("young_modulus", Pressure::DIMS),
            requirement("poisson_ratio", Dims::NONE),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("direct resolution");
    let value = |name: &str| state.property(name).expect("resolved").value_si();
    let (rho, e, nu) = (
        value("density"),
        value("young_modulus"),
        value("poisson_ratio"),
    );
    // The example's Maxwell region must match its section's equilibrium law,
    // so the aluminum files carry the resolved modulus and Poisson ratio there.
    let aluminum = |row: &str| {
        let region = "region 0 4000000000 0.3 ";
        let text = with_section(row);
        assert_eq!(
            text.matches(region).count(),
            1,
            "the example's one Maxwell region"
        );
        text.replace(region, &format!("region 0 {e:?} {nu:?} "))
    };
    // The same thin aluminum blade, once by record and once by the directly
    // resolved numbers as a labelled fixture (Debug formatting round-trips).
    let by_record = load_with_aluminum(&aluminum(
        "section material isotropic 0.00012 aluminum.fsmatpk as-tabulated T=300",
    ))
    .expect("pack-backed section");
    let by_numbers = load_with_aluminum(&aluminum(&format!(
        "section isotropic 0.00012 {rho:?} {e:?} {nu:?} fixture direct-resolution"
    )))
    .expect("direct numbers");
    let reduction = |p: &PlateValvePerformance| {
        let plate = p
            .renderer()
            .system()
            .aperture()
            .plate_reduction()
            .expect("plate");
        [
            plate.mass_kg(),
            plate.stiffness_n_m(),
            plate.pressure_area_m2(),
        ]
        .map(f64::to_bits)
    };
    assert_eq!(
        reduction(&by_record),
        reduction(&by_numbers),
        "the record resolves to the same plate, bit for bit"
    );
    let provenance = &by_record.section_provenance()[0];
    assert!(
        provenance.starts_with("pack aluminum.fsmatpk ")
            && provenance.contains(&pack.content_hash().to_hex())
            && provenance.ends_with("state as-tabulated at T=300.0"),
        "{provenance}"
    );
    assert_eq!(
        by_numbers.section_provenance(),
        ["fixture direct-resolution"]
    );
    // The state must be named exactly, and a property the pack lacks refuses by name.
    for (row, needle) in [
        (
            "section material isotropic 0.00012 aluminum.fsmatpk water-soaked T=300",
            "not the named state 'water-soaked'",
        ),
        (
            // The resolver names the first missing property in canonical order.
            "section material orthotropic12 0.00012 0 aluminum.fsmatpk as-tabulated T=300",
            "UnknownProperty { property: \"poisson_ratio_12\" }",
        ),
        (
            "section material isotropic 0.00012 missing.fsmatpk as-tabulated T=300",
            "no pack missing.fsmatpk",
        ),
        (
            "section material isotropic 0.00012 ../aluminum.fsmatpk as-tabulated T=300",
            "relative",
        ),
    ] {
        let error = match load_with_aluminum(&aluminum(row)) {
            Err(e) => e,
            Ok(_) => panic!("{row} must refuse"),
        };
        assert!(error.contains(needle), "{row}: {error}");
    }
}

const PLATE_MESH: &str = include_str!("../examples/plate-mesh.performance");
const PLATE_REGION: &str = "section 0.004 600.0 11000000000.0 1400000000.0 0.3 800000000.0 0.35 fixture authored-illustrative-plate-not-measured";

fn plate_mesh_with(row: &str) -> String {
    assert_eq!(
        PLATE_MESH.matches(PLATE_REGION).count(),
        1,
        "the example's second region"
    );
    PLATE_MESH.replace(PLATE_REGION, row)
}

fn plate_with_aluminum(
    text: &str,
) -> Result<fs_couple::render::plate::file::PlatePerformance, String> {
    fs_couple::render::plate::file::PlatePerformance::from_bytes_with_all_sources(
        text.as_bytes(),
        64,
        None,
        |path| {
            if path == "aluminum.fsmatpk" {
                Ok(aluminum_card_pack())
            } else {
                Err(format!("no pack {path}"))
            }
        },
    )
    .map_err(|e| e.to_string())
}

/// The flat-plate grammar takes the same record and label: an unlabelled region
/// refuses, the labelled example keeps its labels, and a pack-backed region
/// renders with its pack, hash, state and point in the section provenance.
#[test]
fn plate_mesh_sections_take_the_same_record_or_fixture_label() {
    use fs_couple::pcm_wav::observation::PressureRenderer;
    let labelled = plate_with_aluminum(PLATE_MESH).expect("labelled example");
    assert_eq!(
        labelled.section_provenance(),
        ["fixture authored-illustrative-plate-not-measured"; 2]
    );
    let bare =
        plate_mesh_with("section 0.004 600.0 11000000000.0 1400000000.0 0.3 800000000.0 0.35");
    let error = plate_with_aluminum(&bare)
        .err()
        .expect("an unlabelled region refuses");
    assert!(error.contains("fixture LABEL"), "{error}");
    let packed = plate_with_aluminum(&plate_mesh_with(
        "section material isotropic 0.0018 aluminum.fsmatpk as-tabulated T=300",
    ))
    .expect("pack-backed region");
    let provenance = &packed.section_provenance()[1];
    assert!(
        provenance.starts_with("pack aluminum.fsmatpk ")
            && provenance.contains(&aluminum_card_pack().content_hash().to_hex())
            && provenance.ends_with("state as-tabulated at T=300.0"),
        "{provenance}"
    );
    assert!(packed.info().retained_modes > 0);
    let mut renderer = packed.into_renderer();
    let mut block = vec![0.0; 64];
    renderer.block(&mut block).expect("renders");
    assert!(block.iter().all(|x| x.is_finite()));
    let wrong =
        plate_mesh_with("section material isotropic 0.0018 aluminum.fsmatpk water-soaked T=300");
    let error = plate_with_aluminum(&wrong)
        .err()
        .expect("a wrong state refuses");
    assert!(
        error.contains("not the named state 'water-soaked'"),
        "{error}"
    );
    let default = fs_couple::render::plate::file::PlatePerformance::from_bytes(
        plate_mesh_with("section material isotropic 0.0018 aluminum.fsmatpk as-tabulated T=300")
            .as_bytes(),
        64,
    );
    assert!(
        default
            .err()
            .expect("no loader")
            .to_string()
            .contains("material-card loader")
    );
}
