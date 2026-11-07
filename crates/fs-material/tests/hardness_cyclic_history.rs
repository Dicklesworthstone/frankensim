//! Synthetic original-source checks on native cyclic histories; no calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::cyclic_history::{ResolvedCyclicHardness, SourceBoundHardnessLaw};
use fs_material::hardness::history_field::{AcceptedHardnessField, HardnessFieldError};
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::plastic::chaboche::source::{
    BackstressProperties, CHABOCHE_BACKSTRESS_MODULUS_PROPERTY, CHABOCHE_DYNAMIC_RECOVERY_PROPERTY,
    ResolvedChabocheBranches, resolve_chaboche_branches, resolve_chaboche_j2,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY,
};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY, resolve_voce_j2,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_material::{SmallStrainLaw, Voigt};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec};

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const GATES: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-9,
    yield_pa: 1e-9,
};
const TEMPERATURE: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::AbsoluteTemperature,
    ValueForm::Static,
));

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", TEMPERATURE, temperature)
        .unwrap()
}

fn fixture(equal_branches: bool) -> (MaterialCard, PropertyKey) {
    let provenance = Provenance {
        source: "synthetic cyclic-hardness fixture; no measured data".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let environment = ValidityDomain::unconstrained().with_quantity("T", TEMPERATURE, 280.0, 380.0);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.3),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 30.0),
        (KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 20.0),
        (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, 0.6),
        (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, 80.0),
        (CHABOCHE_BACKSTRESS_MODULUS_PROPERTY, Pressure::DIMS, 300.0),
        (CHABOCHE_DYNAMIC_RECOVERY_PROPERTY, Dims::NONE, 250.0),
        ("C1", Pressure::DIMS, 300.0),
        ("gamma1", Dims::NONE, 250.0),
        (
            "C2",
            Pressure::DIMS,
            if equal_branches { 300.0 } else { 80.0 },
        ),
        (
            "gamma2",
            Dims::NONE,
            if equal_branches { 250.0 } else { 30.0 },
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: environment.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic cyclic specimen at its accumulated plastic strain".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"cyclic hardness"),
            caveats: "numerical regression only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic spherical indenter, radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic projected-area protocol revision 1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: ALPHA.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 100.0), (0.01, 300.0), (0.08, 150.0)],
                dims: Pressure::DIMS,
            },
            validity: environment.with_quantity(
                ALPHA,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                0.08,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    (
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic cyclic solid".into(),
                phase: "solid".into(),
                process: "synthetic independent cyclic and indentation tests".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        key,
    )
}

fn branches(card: &MaterialCard, temperature: f64) -> ResolvedChabocheBranches<2> {
    resolve_chaboche_branches(
        card,
        &point(temperature),
        SELECT,
        BASIS,
        0.05,
        [
            BackstressProperties::new("C1", "gamma1"),
            BackstressProperties::new("C2", "gamma2"),
        ],
    )
    .unwrap()
}

fn check_endpoint<L: SourceBoundHardnessLaw>(
    model: &ResolvedCyclicHardness<L>,
    card: &MaterialCard,
    key: &PropertyKey,
    strain: &Voigt,
    state: &L::State,
    alpha: f64,
) -> f64 {
    let hardness = model
        .hardness_pa_at_accepted(model.plastic(), strain, state, GATES.yield_pa)
        .unwrap();
    let expected = if alpha <= 0.01 {
        100.0 + 20_000.0 * alpha
    } else {
        300.0 - 150.0 * (alpha - 0.01) / 0.07
    };
    assert!((hardness - expected).abs() < 1e-10);
    let typed = model
        .resolve_hardness_at_accepted(card, model.plastic(), strain, state, GATES.yield_pa)
        .unwrap();
    assert!((typed.hardness_pa() - hardness).abs() < 1e-10);
    assert_eq!(typed.key(), key);
    assert_eq!(typed.resolved().state_coordinate(ALPHA), Some(alpha));
    assert_eq!(typed.resolved().state_coordinate("T"), Some(300.0));
    assert_eq!(
        typed.property().answer().receipt.selected,
        model.hardness_source().property().answer().receipt.selected
    );
    hardness
}

#[test]
fn voce_loading_unloading_reversal_observes_saturation_without_linear_projection() {
    let (card, key) = fixture(false);
    let plastic = resolve_voce_j2(&card, &point(300.0), SELECT, BASIS, 0.05).unwrap();
    let model = ResolvedCyclicHardness::new(plastic, &card, key.clone(), SELECT, 0.08).unwrap();
    let mut state = model.plastic().initial_state();
    let mut strain = [0.0; 6];
    let mut previous_hardness = 100.0;
    for (step, axial) in [0.008, 0.0078, -0.008, 0.008].into_iter().enumerate() {
        let target = [axial, 0.0, 0.0, 0.0, 0.0, 0.0];
        let next = model
            .plastic()
            .law()
            .advance_with_energy(&strain, &state, &target, GATES)
            .unwrap();
        let before = next.clone();
        let hardness = check_endpoint(&model, &card, &key, &target, &next.state, next.state.alpha);
        assert_eq!(next, before);
        if step == 0 {
            assert!(
                model
                    .plastic()
                    .backbone()
                    .law()
                    .validate_committed_history(&target, &next.state, GATES.yield_pa)
                    .is_err(),
                "the nonlinear accepted state must not be projected onto bare J2"
            );
        }
        if step == 1 {
            assert_eq!(next.state.alpha, state.alpha);
            assert_eq!(hardness, previous_hardness);
        }
        assert!(next.state.alpha >= state.alpha);
        strain = target;
        state = next.state;
        previous_hardness = hardness;
    }
    assert!(state.alpha > 0.01);
}

#[test]
fn single_and_multiple_chaboche_branches_keep_native_nonproportional_recovery() {
    let (card, key) = fixture(false);
    let single = resolve_chaboche_j2(&card, &point(300.0), SELECT, BASIS, 0.05).unwrap();
    let single = ResolvedCyclicHardness::new(single, &card, key.clone(), SELECT, 0.08).unwrap();
    let target = [0.007, -0.002, -0.001, 0.0015, 0.0, 0.0];
    let next = single
        .plastic()
        .law()
        .advance_with_energy(&[0.0; 6], &single.plastic().initial_state(), &target, GATES)
        .unwrap();
    check_endpoint(
        &single,
        &card,
        &key,
        &target,
        &next.state,
        next.state.plastic.alpha,
    );
    assert!(
        next.state.backstress_pa[0]
            .iter()
            .any(|value| *value != 0.0)
    );

    let model =
        ResolvedCyclicHardness::new(branches(&card, 300.0), &card, key.clone(), SELECT, 0.08)
            .unwrap();
    let mut state = model.plastic().initial_state();
    let mut strain = [0.0; 6];
    let mut recovery = 0.0;
    for target in [
        target,
        [0.003, 0.006, -0.002, -0.004, 0.001, 0.002],
        [-0.006, 0.001, 0.002, 0.004, -0.002, 0.001],
    ] {
        let next = model
            .plastic()
            .law()
            .advance_with_energy(&strain, &state, &target, GATES)
            .unwrap();
        let before = next.clone();
        check_endpoint(
            &model,
            &card,
            &key,
            &target,
            &next.state,
            next.state.plastic.alpha,
        );
        assert_eq!(next, before);
        assert!(next.state.plastic.alpha > state.plastic.alpha);
        assert_ne!(next.state.backstress_pa[0], next.state.backstress_pa[1]);
        recovery += next.recovery_dissipation_j_m3;
        state = next.state;
        strain = target;
    }
    assert!(recovery > 0.0);
}

#[test]
fn branch_source_order_is_binding_even_when_numeric_branches_are_equal() {
    let (card, key) = fixture(true);
    let model =
        ResolvedCyclicHardness::new(branches(&card, 300.0), &card, key, SELECT, 0.08).unwrap();
    let reordered = resolve_chaboche_branches(
        &card,
        &point(300.0),
        SELECT,
        BASIS,
        0.05,
        [
            BackstressProperties::new("C2", "gamma2"),
            BackstressProperties::new("C1", "gamma1"),
        ],
    )
    .unwrap();
    assert_eq!(model.plastic().law(), reordered.law());
    assert!(
        model
            .hardness_pa_at_accepted(
                &reordered,
                &[0.0; 6],
                &reordered.initial_state(),
                GATES.yield_pa,
            )
            .is_err()
    );
    let changed_temperature = branches(&card, 301.0);
    assert_eq!(model.plastic().law(), changed_temperature.law());
    assert!(
        model
            .hardness_pa_at_accepted(
                &changed_temperature,
                &[0.0; 6],
                &changed_temperature.initial_state(),
                GATES.yield_pa,
            )
            .is_err()
    );
}

#[test]
fn complete_state_yield_and_prepared_support_refuse_before_endpoint_export() {
    let (card, key) = fixture(false);
    let native = branches(&card, 300.0);
    let strain = [0.008, 0.0, 0.0, 0.001, 0.0, 0.0];
    let accepted = native
        .law()
        .response(&strain, &native.initial_state())
        .unwrap()
        .state;
    let model =
        ResolvedCyclicHardness::new(native.clone(), &card, key.clone(), SELECT, 0.08).unwrap();
    let before = accepted.clone();
    let mut corrupted = accepted.clone();
    corrupted.backstress_pa[1][0] += 0.1;
    assert!(
        model
            .hardness_pa_at_accepted(&native, &strain, &corrupted, GATES.yield_pa)
            .is_err()
    );
    corrupted = accepted.clone();
    corrupted.plastic.alpha = f64::NAN;
    assert!(
        model
            .hardness_pa_at_accepted(&native, &strain, &corrupted, GATES.yield_pa)
            .is_err()
    );
    let virgin = native.initial_state();
    native.law().validate_history(&strain, &virgin).unwrap();
    assert!(
        model
            .hardness_pa_at_accepted(&native, &strain, &virgin, GATES.yield_pa)
            .is_err(),
        "incoming-history admissibility alone does not certify an accepted cyclic yield state"
    );
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .hardness_pa_at_accepted(&native, &strain, &accepted, tolerance)
                .is_err()
        );
    }
    let bounded = ResolvedCyclicHardness::new(
        native.clone(),
        &card,
        key,
        SELECT,
        accepted.plastic.alpha * 0.5,
    )
    .unwrap();
    assert!(
        bounded
            .hardness_pa_at_accepted(&native, &strain, &accepted, GATES.yield_pa)
            .is_err()
    );
    assert!(
        bounded
            .resolve_hardness_at_accepted(&card, &native, &strain, &accepted, GATES.yield_pa)
            .is_err()
    );
    assert_eq!(accepted, before);
}

#[test]
fn compact_cyclic_field_binds_native_family_and_ordered_equal_branch_sources() {
    let (card, key) = fixture(true);
    let native = branches(&card, 300.0);
    let model =
        ResolvedCyclicHardness::new(native.clone(), &card, key.clone(), SELECT, 0.08).unwrap();
    let strains = [
        [0.007, -0.002, -0.001, 0.0015, 0.0, 0.0],
        [0.003, 0.006, -0.002, -0.004, 0.001, 0.002],
    ];
    let first = native
        .law()
        .response(&strains[0], &native.initial_state())
        .unwrap()
        .state;
    let second = native.law().response(&strains[1], &first).unwrap().state;
    let states = [first, second];
    let field = AcceptedHardnessField::from_cyclic(
        &model,
        (0..2).map(|point| (&native, &strains[point], &states[point])),
        GATES.yield_pa,
    )
    .unwrap();
    assert_eq!(field.source_identity(), model.identity());
    assert_eq!(field.mechanical_source(), native.material());
    assert_eq!(field.hardness_source(), model.hardness_source());
    for point in 0..2 {
        assert_eq!(
            field.sample(point).unwrap().alpha(),
            states[point].plastic.alpha
        );
        assert_eq!(
            field.hardness_pa(point),
            Some(
                model
                    .hardness_pa_at_accepted(
                        &native,
                        &strains[point],
                        &states[point],
                        GATES.yield_pa,
                    )
                    .unwrap()
            )
        );
    }
    let rebuilt =
        ResolvedCyclicHardness::new(native.clone(), &card, key.clone(), SELECT, 0.08).unwrap();
    assert_eq!(model.identity(), rebuilt.identity());
    let narrower =
        ResolvedCyclicHardness::new(native.clone(), &card, key.clone(), SELECT, 0.04).unwrap();
    assert_ne!(model.identity(), narrower.identity());
    let reordered = resolve_chaboche_branches(
        &card,
        &point(300.0),
        SELECT,
        BASIS,
        0.05,
        [
            BackstressProperties::new("C2", "gamma2"),
            BackstressProperties::new("C1", "gamma1"),
        ],
    )
    .unwrap();
    let reordered_model =
        ResolvedCyclicHardness::new(reordered.clone(), &card, key.clone(), SELECT, 0.08).unwrap();
    assert_eq!(native.law(), reordered.law());
    assert_ne!(model.identity(), reordered_model.identity());
    assert!(matches!(
        AcceptedHardnessField::from_cyclic(
            &model,
            (0..2).map(|point| (
                if point == 0 { &native } else { &reordered },
                &strains[point],
                &states[point],
            )),
            GATES.yield_pa,
        ),
        Err(HardnessFieldError::Point { point: 1, .. })
    ));

    // Same full source bundle and prepared hardness, different actual native
    // family. The field identity must retain which complete history was checked.
    let voce =
        ResolvedCyclicHardness::new(native.backbone().clone(), &card, key, SELECT, 0.08).unwrap();
    assert_eq!(voce.plastic().material(), native.material());
    assert_ne!(voce.identity(), model.identity());
    let virgin = voce.plastic().initial_state();
    let voce_field = AcceptedHardnessField::from_cyclic(
        &voce,
        std::iter::once((voce.plastic(), &[0.0; 6], &virgin)),
        GATES.yield_pa,
    )
    .unwrap();
    assert_eq!(voce_field.source_identity(), voce.identity());
    assert_eq!(voce_field.sample(0).unwrap().alpha(), 0.0);
    assert_eq!(voce_field.hardness_pa(0), Some(100.0));
}
