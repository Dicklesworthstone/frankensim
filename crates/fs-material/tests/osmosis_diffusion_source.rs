//! Material-backed liquid kinetics through the unchanged shared physical owner.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId,
    LawParameter, MaterialCard, MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::{
    osmosis::{
        AqueousSolution,
        diffusion::{
            SolutionDiffusionError, SolutionDiffusionSolve,
            source::{
                BINARY_SOLUTION_DIFFUSIVITY_DIMS as DIMS, BinarySolutionDiffusionSource,
                MaterialSolutionMaxwellStefan, SolutionDiffusionSourceAxes,
            },
        },
        nonideal::VirialSolutionModel,
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{
    MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{SpeciesId, transport::maxwell_stefan::MaxwellStefanMolarSolveV1};
use std::{collections::BTreeMap, sync::Arc};

const T: f64 = 298.15;
const VBAR: f64 = 18e-6;
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn species(name: &str) -> (SpeciesId, MolarMass) {
    (
        id(name),
        MolarMass::new(if name == "water" { 0.018 } else { 0.06 }),
    )
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured binary mobility data".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn model(n: usize) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, VBAR, 1.0, vec![200.0; n], vec![0.0; n * n]).unwrap())
}
fn solution(model: &Arc<VirialSolutionModel>, c: &[f64]) -> AqueousSolution {
    AqueousSolution::nonideal(model.clone(), 1e-3, c.iter().map(|v| v * 1e-3).collect()).unwrap()
}
fn curve(name: &str, axis: &str, range: [f64; 2], scale: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, DIMS),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: typed(QuantityKind::AmountConcentration).dims(),
            knots: range
                .into_iter()
                .map(|c| (c, scale * (1.0 + c / 100.0)))
                .collect(),
            dims: DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", typed(QuantityKind::AbsoluteTemperature), 290.0, 310.0)
            .with_quantity(
                axis,
                typed(QuantityKind::AmountConcentration),
                range[0],
                range[1],
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: provenance(),
    }
}
fn card(claims: Vec<PropertyClaim>) -> (Arc<MaterialCard>, Vec<ClaimId>) {
    let mut set = ClaimSet::new();
    let ids = claims
        .into_iter()
        .map(|claim| set.insert_claim(claim).unwrap())
        .collect();
    (
        Arc::new(
            MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "explicit solution pair".into(),
                    phase: "liquid".into(),
                    process: "manufactured reference".into(),
                    revision: 0,
                },
                set,
                vec![],
            )
            .unwrap(),
        ),
        ids,
    )
}
fn axes(names: &[&str]) -> SolutionDiffusionSourceAxes {
    SolutionDiffusionSourceAxes {
        temperature: "T".into(),
        concentrations: core::iter::once("water")
            .chain(names.iter().copied())
            .map(|name| (id(name), format!("c_{name}")))
            .collect(),
    }
}
fn binding(
    pair: [&str; 2],
    material: &Arc<MaterialCard>,
    property: &str,
) -> BinarySolutionDiffusionSource {
    BinarySolutionDiffusionSource {
        species: pair.map(id),
        material: material.clone(),
        property: PropertyKey::new(property, DIMS),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn controls() -> SolutionDiffusionSolve {
    SolutionDiffusionSolve {
        friction: MaxwellStefanMolarSolveV1::default(),
        entropy_tolerance_w_m3_k: 1e-10,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn g1_actual_concentration_changes_pair_flux_and_returns_original_receipts() {
    let model = model(1);
    let (material, ids) = card(vec![curve("D", "c_A", [10.0, 150.0], 1e-9)]);
    let law = MaterialSolutionMaxwellStefan::new(
        model.clone(),
        species("water"),
        vec![species("A")],
        vec![binding(["water", "A"], &material, "D")],
        axes(&["A"]),
        &solution(&model, &[40.0]),
        || false,
    )
    .unwrap();
    let mut flux = [0.0];
    let report = law
        .flux_into(
            &solution(&model, &[100.0]),
            &solution(&model, &[20.0]),
            0.005,
            controls(),
            &mut law.workspace().unwrap(),
            &mut flux,
            || false,
        )
        .unwrap();
    let expected = 1.6e-9 * 60.0 * (1.0 + 60.0 * VBAR) * 5.0_f64.ln() / 0.005;
    close(flux[0], expected, 2e-18);
    assert_eq!(law.sources()[0].selected_claim(), ids[0]);
    let pair = &report.sources.pairs()[0];
    for (receipt, expected) in [
        (&pair.left, 2e-9),
        (&pair.face, 1.6e-9),
        (&pair.right, 1.2e-9),
    ] {
        let property = &receipt.properties()[0];
        close(property.value_si(), expected, 1e-24);
        assert_eq!(property.answer().receipt.selected, ids[0]);
        assert_eq!(receipt.card_identity(), material.content_hash());
    }
    let mut frozen = [0.0];
    law.reference_law()
        .flux_into(
            &solution(&model, &[100.0]),
            &solution(&model, &[20.0]),
            0.005,
            controls(),
            &mut law.workspace().unwrap(),
            &mut frozen,
            || false,
        )
        .unwrap();
    close(flux[0] / frozen[0], 8.0 / 7.0, 1e-14);
    assert!(report.constitutive.chemical_entropy_production_w_m3_k > 0.0);
}

#[test]
fn g4_pinned_source_refuses_unsupported_endpoints_despite_an_admissible_face() {
    let model = model(1);
    let (material, ids) = card(vec![
        curve("D", "c_A", [40.0, 80.0], 1e-9),
        curve("D", "c_A", [10.0, 150.0], 2e-9),
    ]);
    let mut source = binding(["water", "A"], &material, "D");
    source.selection = MaterialPropertySelection::PinnedByProperty(vec![("D".into(), ids[0])]);
    let law = MaterialSolutionMaxwellStefan::new(
        model.clone(),
        species("water"),
        vec![species("A")],
        vec![source],
        axes(&["A"]),
        &solution(&model, &[60.0]),
        || false,
    )
    .unwrap();
    let mut output = [123.0];
    let error = law
        .flux_into(
            &solution(&model, &[100.0]),
            &solution(&model, &[20.0]),
            0.005,
            controls(),
            &mut law.workspace().unwrap(),
            &mut output,
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, SolutionDiffusionError::Material(_)));
    assert_eq!(output, [123.0]);
    assert_eq!(law.sources()[0].selected_claim(), ids[0]);
    assert_eq!(
        law.flux_into(
            &solution(&model, &[70.0]),
            &solution(&model, &[50.0]),
            0.005,
            controls(),
            &mut law.workspace().unwrap(),
            &mut output,
            || true
        ),
        Err(SolutionDiffusionError::Cancelled)
    );
    assert_eq!(output, [123.0]);
}

#[test]
fn g3_complete_pair_and_coordinate_permutations_preserve_original_species_order() {
    let model = model(2);
    let (material, _) = card(vec![
        curve("Dz", "c_Z", [1.0, 150.0], 1e-9),
        curve("Da", "c_A", [1.0, 150.0], 3e-9),
        curve("Dza", "c_A", [1.0, 150.0], 5e-10),
    ]);
    let run = |permuted| {
        let names = if permuted { ["A", "Z"] } else { ["Z", "A"] };
        let (left, right) = if permuted {
            ([100.0, 20.0], [10.0, 21.0])
        } else {
            ([20.0, 100.0], [21.0, 10.0])
        };
        let mut pairs = vec![
            binding(["Z", "water"], &material, "Dz"),
            binding(["water", "A"], &material, "Da"),
            binding(["Z", "A"], &material, "Dza"),
        ];
        let mut coordinates = axes(&names);
        if permuted {
            pairs.reverse();
            coordinates.concentrations.reverse();
        }
        let law = MaterialSolutionMaxwellStefan::new(
            model.clone(),
            species("water"),
            names.map(species).to_vec(),
            pairs,
            coordinates,
            &solution(&model, &left),
            || false,
        )
        .unwrap();
        let mut flux = [0.0; 2];
        law.flux_into(
            &solution(&model, &left),
            &solution(&model, &right),
            0.01,
            controls(),
            &mut law.workspace().unwrap(),
            &mut flux,
            || false,
        )
        .unwrap();
        flux
    };
    let original = run(false);
    let permuted = run(true);
    close(original[0], permuted[1], 2e-18);
    close(original[1], permuted[0], 2e-18);
}

#[test]
fn g4_actual_joint_surface_support_is_required_beyond_the_validity_box() {
    let model = model(1);
    let mut claim = curve("D", "c_A", [10.0, 150.0], 1e-9);
    claim.value = PropertyValue::Surface2D {
        axes: ["T".into(), "c_A".into()],
        axis_dims: [
            typed(QuantityKind::AbsoluteTemperature).dims(),
            typed(QuantityKind::AmountConcentration).dims(),
        ],
        samples: vec![
            ([290.0, 10.0], 1e-9),
            ([310.0, 10.0], 2e-9),
            ([310.0, 150.0], 3e-9),
        ],
        triangles: vec![[0, 1, 2]],
        dims: DIMS,
    };
    let (card, _) = card(vec![claim]);
    // At T=298.15 the supported triangle ends at c=67.05, though the declared
    // rectangular validity permits c=150. A c=60 face cannot qualify c=100.
    let law = MaterialSolutionMaxwellStefan::new(
        model.clone(),
        species("water"),
        vec![species("A")],
        vec![binding(["water", "A"], &card, "D")],
        axes(&["A"]),
        &solution(&model, &[60.0]),
        || false,
    )
    .unwrap();
    let mut output = [456.0];
    assert!(
        law.flux_into(
            &solution(&model, &[100.0]),
            &solution(&model, &[20.0]),
            0.005,
            controls(),
            &mut law.workspace().unwrap(),
            &mut output,
            || false
        )
        .is_err()
    );
    assert_eq!(output, [456.0]);
}

#[test]
fn g4_no_missing_pair_or_unowned_pressure_default_is_admitted() {
    let model = model(1);
    let reference = solution(&model, &[60.0]);
    let (material, _) = card(vec![curve("D", "c_A", [10.0, 150.0], 1e-9)]);
    assert!(
        MaterialSolutionMaxwellStefan::new(
            model.clone(),
            species("water"),
            vec![species("A")],
            vec![],
            axes(&["A"]),
            &reference,
            || false
        )
        .is_err()
    );
    let mut pressure = curve("D", "c_A", [10.0, 150.0], 1e-9);
    pressure.validity = pressure.validity.with_quantity(
        "pressure",
        typed(QuantityKind::Pressure),
        90000.0,
        110000.0,
    );
    let (pressure, _) = card(vec![pressure]);
    assert!(
        MaterialSolutionMaxwellStefan::new(
            model.clone(),
            species("water"),
            vec![species("A")],
            vec![binding(["water", "A"], &pressure, "D")],
            axes(&["A"]),
            &reference,
            || false
        )
        .is_err()
    );
    let mut wrong = binding(["water", "A"], &material, "D");
    wrong.property = PropertyKey::new("D", typed(QuantityKind::Pressure).dims());
    assert!(
        MaterialSolutionMaxwellStefan::new(
            model,
            species("water"),
            vec![species("A")],
            vec![wrong],
            axes(&["A"]),
            &reference,
            || false
        )
        .is_err()
    );
}

#[test]
fn g4_retained_source_charge_includes_unselected_data_observations_and_models() {
    let base = curve("D", "c_A", [10.0, 150.0], 1e-9);
    let (small, _) = card(vec![base.clone()]);
    let mut claims = ClaimSet::new();
    claims.insert_claim(base).unwrap();
    let mut unused = curve("unused", "c_A", [10.0, 150.0], 1e-9);
    if let PropertyValue::Curve { knots, .. } = &mut unused.value {
        *knots = (0..1000).map(|i| (10.0 + 0.1 * i as f64, 1e-9)).collect();
    }
    claims.insert_claim(unused).unwrap();
    claims
        .register_observation(ObservationDataset {
            specimen: "unselected specimen".repeat(100),
            method: "manufactured method".into(),
            artifact: small.content_hash(),
            caveats: "full retained source".into(),
            provenance: provenance(),
        })
        .unwrap();
    let model_card = ConstitutiveModelCard {
        law: LawId("unselected calibration".into()),
        law_version: 1,
        parameters: BTreeMap::from([(
            "unselected_parameter".into(),
            LawParameter {
                value: 1.0,
                dims: DIMS,
            },
        )]),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained(),
        sources: vec![small.content_hash()],
        provenance: provenance(),
    };
    let large =
        Arc::new(MaterialCard::assemble(small.id().clone(), claims, vec![model_card]).unwrap());
    assert!(
        large.estimated_retained_payload_bytes() > small.estimated_retained_payload_bytes() + 16000
    );
    let model = model(1);
    let make = |card: &Arc<MaterialCard>| {
        MaterialSolutionMaxwellStefan::new(
            model.clone(),
            species("water"),
            vec![species("A")],
            vec![binding(["water", "A"], card, "D")],
            axes(&["A"]),
            &solution(&model, &[60.0]),
            || false,
        )
        .unwrap()
    };
    assert_eq!(
        make(&large).estimated_heap_payload_bytes() - make(&small).estimated_heap_payload_bytes(),
        large.estimated_retained_payload_bytes() - small.estimated_retained_payload_bytes()
    );
}
