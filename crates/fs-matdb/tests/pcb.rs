//! G0/G3 battery for the evidence-carrying PCB laminate first rung.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, CopperCoverage, InterpolationPolicy, MaterialCard, MaterialStateId,
    PCB_THERMAL_CONDUCTIVITY_DIMS, PcbConductivityDatum, PcbHomogenizationError, PcbLayer,
    PcbPrincipalFrame, PcbScaleSeparation, PcbStackup, PcbViaCorrection, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};

const PROPERTY: &str = "thermal_conductivity";

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} differs from {expected:.17e} by more than {tolerance:.3e}"
    );
}

fn material_card(
    chemistry: &str,
    process: &str,
    conductivity: f64,
    uncertainty: UncertaintyModel,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(PROPERTY, PCB_THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: conductivity,
                dims: PCB_THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 400.0),
            uncertainty,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("{chemistry} thermal-data fixture"),
                license: "test-only".to_string(),
                artifact: None,
            },
        })
        .expect("conductivity claim");
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "solid".to_string(),
            process: process.to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("material card")
}

fn datum(card: &MaterialCard) -> PcbConductivityDatum {
    let point = QueryPoint::new().with("T", 300.0).expect("query point");
    PcbConductivityDatum::from_card(card, PROPERTY, &point, SelectionPolicy::SingleClaimOnly)
        .expect("conductivity datum")
}

fn coverage(source: &str, nominal: f64, lower: f64, upper: f64) -> CopperCoverage {
    CopperCoverage::new(
        source,
        nominal,
        lower,
        upper,
        Provenance {
            source: format!("ODB++ coverage export {source}"),
            license: "test-only".to_string(),
            artifact: None,
        },
    )
    .expect("coverage")
}

fn exact_cards() -> (MaterialCard, MaterialCard) {
    (
        material_card(
            "C11000",
            "rolled-foil",
            400.0,
            UncertaintyModel::HalfWidth {
                half_width: 0.0,
                confidence: 0.95,
            },
        ),
        material_card(
            "FR4",
            "cured-laminate",
            0.25,
            UncertaintyModel::HalfWidth {
                half_width: 0.0,
                confidence: 0.95,
            },
        ),
    )
}

fn reference_stackup() -> PcbStackup {
    let (copper, fr4) = exact_cards();
    let plane = PcbLayer::new(
        "L1-plane",
        0.2e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/L1", 1.0, 1.0, 1.0),
    )
    .expect("plane");
    let core = PcbLayer::new(
        "core",
        0.8e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/core", 0.0, 0.0, 0.0),
    )
    .expect("core");
    PcbStackup::new(
        "reference-two-layer",
        vec![plane, core],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(25.0e-6, 0.05).expect("separation"),
    )
    .expect("stackup")
}

#[test]
fn hand_calculated_parallel_series_and_structural_bounds_match() {
    let result = reference_stackup().homogenize().expect("homogenize");
    let expected_parallel = 0.2 * 400.0 + 0.8 * 0.25;
    let expected_series = 1.0 / (0.2 / 400.0 + 0.8 / 0.25);
    close(
        result.principal().nominal_w_mk[0],
        expected_parallel,
        1.0e-13,
    );
    close(
        result.principal().nominal_w_mk[1],
        expected_parallel,
        1.0e-13,
    );
    close(result.principal().nominal_w_mk[2], expected_series, 1.0e-15);
    close(
        result.structural_bounds().reuss_w_mk,
        expected_series,
        1.0e-15,
    );
    close(
        result.structural_bounds().voigt_w_mk,
        expected_parallel,
        1.0e-13,
    );
    assert!(
        result.structural_bounds().reuss_w_mk <= result.principal().nominal_w_mk[2]
            && result.principal().nominal_w_mk[2] <= result.principal().nominal_w_mk[0]
            && result.principal().nominal_w_mk[0] <= result.structural_bounds().voigt_w_mk
    );
    assert_eq!(result.via_correction(), PcbViaCorrection::NotModeled);
    assert!(result.material_uncertainty_complete());
    assert_eq!(result.material_uses().len(), 4);
    assert_eq!(result.coverage_influences().len(), 2);
    println!(
        "{{\"suite\":\"fs-matdb-pcb\",\"case\":\"hand-stackup\",\"status\":\"pass\",\
         \"k_in_plane\":{},\"k_through\":{},\"reuss\":{},\"voigt\":{}}}",
        result.principal().nominal_w_mk[0],
        result.principal().nominal_w_mk[2],
        result.structural_bounds().reuss_w_mk,
        result.structural_bounds().voigt_w_mk
    );
}

#[test]
fn coverage_bounds_propagate_as_one_shared_directional_source() {
    let (copper, fr4) = exact_cards();
    let layer = PcbLayer::new(
        "signal-layer",
        1.0e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/signal", 0.50, 0.40, 0.60),
    )
    .expect("layer");
    let result = PcbStackup::new(
        "bounded-coverage",
        vec![layer],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(10.0e-6, 0.02).expect("separation"),
    )
    .expect("stackup")
    .homogenize()
    .expect("homogenize");

    let influence = &result.coverage_influences()[0];
    assert_eq!(influence.source_id, "coverage/signal");
    assert_eq!(
        influence.principal_at_lower_w_mk[0].to_bits(),
        influence.principal_at_lower_w_mk[1].to_bits()
    );
    assert_eq!(
        influence.principal_at_upper_w_mk[0].to_bits(),
        influence.principal_at_upper_w_mk[1].to_bits()
    );
    assert!(
        influence.principal_at_lower_w_mk[0] < result.principal().nominal_w_mk[0]
            && result.principal().nominal_w_mk[0] < influence.principal_at_upper_w_mk[0]
    );
    assert_eq!(
        influence.principal_at_lower_w_mk[0].to_bits(),
        influence.principal_at_lower_w_mk[2].to_bits(),
        "one physical layer has no series/parallel anisotropy"
    );
    assert_eq!(
        result.principal().lower_w_mk.map(f64::to_bits),
        influence.principal_at_lower_w_mk.map(f64::to_bits)
    );
    assert_eq!(
        result.principal().upper_w_mk.map(f64::to_bits),
        influence.principal_at_upper_w_mk.map(f64::to_bits)
    );
}

#[test]
fn stated_material_bands_and_coverage_boxes_enclose_every_corner() {
    let copper = material_card(
        "C11000",
        "rolled-foil",
        400.0,
        UncertaintyModel::RelativeHalfWidth {
            fraction: 0.05,
            confidence: 0.95,
        },
    );
    let fr4 = material_card(
        "FR4",
        "cured-laminate",
        0.25,
        UncertaintyModel::RelativeHalfWidth {
            fraction: 0.20,
            confidence: 0.95,
        },
    );
    let layer = PcbLayer::new(
        "bounded-layer",
        1.0e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/bounded", 0.5, 0.4, 0.6),
    )
    .expect("layer");
    let result = PcbStackup::new(
        "bounded-materials",
        vec![layer],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(10.0e-6, 0.02).expect("separation"),
    )
    .expect("stackup")
    .homogenize()
    .expect("homogenize");
    assert!(result.material_uncertainty_complete());

    for coverage in [0.4_f64, 0.6] {
        for copper in [380.0, 420.0] {
            for matrix in [0.20, 0.30] {
                let corner = coverage.mul_add(copper, (1.0 - coverage) * matrix);
                assert!(
                    result.principal().lower_w_mk[0] <= corner
                        && corner <= result.principal().upper_w_mk[0]
                );
            }
        }
    }
}

#[test]
fn unstated_material_uncertainty_stays_explicitly_incomplete() {
    let copper = material_card("C11000", "rolled-foil", 400.0, UncertaintyModel::Unstated);
    let fr4 = material_card(
        "FR4",
        "cured-laminate",
        0.25,
        UncertaintyModel::HalfWidth {
            half_width: 0.0,
            confidence: 0.95,
        },
    );
    let layer = PcbLayer::new(
        "layer",
        1.0e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/unstated", 0.5, 0.4, 0.6),
    )
    .expect("layer");
    let result = PcbStackup::new(
        "unstated-material",
        vec![layer],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(10.0e-6, 0.02).expect("separation"),
    )
    .expect("stackup")
    .homogenize()
    .expect("homogenize");
    assert!(!result.material_uncertainty_complete());
}

#[test]
fn zero_coverage_and_single_layer_degenerate_to_the_matrix() {
    let (copper, fr4) = exact_cards();
    let layer = PcbLayer::new(
        "dielectric-only",
        1.0e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/zero", 0.0, 0.0, 0.0),
    )
    .expect("layer");
    let result = PcbStackup::new(
        "matrix-only",
        vec![layer],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(10.0e-6, 0.02).expect("separation"),
    )
    .expect("stackup")
    .homogenize()
    .expect("homogenize");
    assert_eq!(
        result.principal().nominal_w_mk.map(f64::to_bits),
        [0.25_f64; 3].map(f64::to_bits)
    );
    assert_eq!(result.tensor_w_mk()[0][0].to_bits(), 0.25_f64.to_bits());
    assert_eq!(result.tensor_w_mk()[1][1].to_bits(), 0.25_f64.to_bits());
    assert_eq!(result.tensor_w_mk()[2][2].to_bits(), 0.25_f64.to_bits());
}

#[test]
fn frame_rotation_preserves_symmetry_and_moves_principal_directions() {
    let original = reference_stackup().homogenize().expect("original");
    let rotated_frame =
        PcbPrincipalFrame::new([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
            .expect("right-handed rotation");
    let layers = reference_stackup().layers().to_vec();
    let rotated = PcbStackup::new(
        "rotated",
        layers,
        rotated_frame,
        PcbScaleSeparation::new(25.0e-6, 0.05).expect("separation"),
    )
    .expect("stackup")
    .homogenize()
    .expect("homogenize");
    assert_eq!(
        rotated.principal().nominal_w_mk.map(f64::to_bits),
        original.principal().nominal_w_mk.map(f64::to_bits)
    );
    close(
        rotated.tensor_w_mk()[0][0],
        original.principal().nominal_w_mk[2],
        1.0e-15,
    );
    close(
        rotated.tensor_w_mk()[1][1],
        original.principal().nominal_w_mk[1],
        1.0e-13,
    );
    close(
        rotated.tensor_w_mk()[2][2],
        original.principal().nominal_w_mk[0],
        1.0e-13,
    );
    for row in 0..3 {
        for column in 0..3 {
            assert_eq!(
                rotated.tensor_w_mk()[row][column].to_bits(),
                rotated.tensor_w_mk()[column][row].to_bits()
            );
        }
    }
}

#[test]
fn scale_separation_and_malformed_inputs_refuse_typed() {
    let (copper, fr4) = exact_cards();
    let layer = PcbLayer::new(
        "too-coarse",
        1.0e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/coarse", 0.5, 0.4, 0.6),
    )
    .expect("layer");
    let error = PcbStackup::new(
        "outside-domain",
        vec![layer],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(0.2e-3, 0.1).expect("rule"),
    )
    .expect("stackup")
    .homogenize()
    .expect_err("feature/thickness ratio must refuse");
    assert!(matches!(
        error,
        PcbHomogenizationError::ScaleSeparation { .. }
    ));
    assert!(matches!(
        CopperCoverage::new(
            "bad",
            0.5,
            -0.1,
            0.6,
            Provenance {
                source: "fixture".to_string(),
                license: "test-only".to_string(),
                artifact: None,
            },
        ),
        Err(PcbHomogenizationError::InvalidField {
            field: "copper-coverage",
            ..
        })
    ));
    assert!(PcbPrincipalFrame::new([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]).is_err());
}

#[test]
fn stack_order_is_identity_bearing_even_when_the_series_value_is_equal() {
    let stackup = reference_stackup();
    let forward = stackup.homogenize().expect("forward");
    let mut reversed_layers = stackup.layers().to_vec();
    reversed_layers.reverse();
    let reverse = PcbStackup::new(
        "reference-two-layer",
        reversed_layers,
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(25.0e-6, 0.05).expect("separation"),
    )
    .expect("reverse stack")
    .homogenize()
    .expect("reverse");
    assert_eq!(
        forward.principal().nominal_w_mk.map(f64::to_bits),
        reverse.principal().nominal_w_mk.map(f64::to_bits)
    );
    assert_ne!(
        forward.identity(),
        reverse.identity(),
        "physical stack order must remain provenance even when this rung is permutation-invariant"
    );
}
