//! Native Hex8 histories feed local indentation hardness without point averaging.
//! Synthetic coefficients isolate the mechanics/property seam, not calibration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::cyclic_history::ResolvedCyclicHardness;
use fs_material::hardness::history_field::HardnessFieldError;
use fs_material::hardness::plastic_history::{
    ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA, ResolvedJ2Hardness,
};
use fs_material::plastic::chaboche::source::{BackstressProperties, resolve_chaboche_branches};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY, resolve_j2_plasticity,
};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec};
use fs_solid::plastic3::hex::hardness::{HexHardnessError, HexHardnessField};
use fs_solid::plastic3::hex::{HexAssemblyBudget, SmallStrainHexProblem, SmallStrainHexSolution};
use fs_solid::{Formulation, IncrementSettings, SmallStrainLawField};
use fs_tribo::{
    ArchardLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, LocalHardnessArchardSample,
    WearState,
};

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;
const TOL: f64 = 1e-7;
const E: f64 = 1200.0;
const NU: f64 = 0.3;
const HI: f64 = 30.0;
const SLOPE: f64 = 20_000.0;
const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
];
const CELLS: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cancelled_cx(false, f)
}

fn with_cancelled_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate, arena,
            StreamKey { seed: 7, kernel_id: 31, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        ))
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} != {expected}");
}

fn card(region: usize) -> (MaterialCard, PropertyKey) {
    let provenance = Provenance {
        source: format!("synthetic Hex8 hardness region {region}; no measured claim"),
        license: "CC0-1.0".into(), artifact: None,
    };
    let environment = ValidityDomain::unconstrained().with("T", 300.0, 300.0);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims::NONE, NU),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0 + region as f64),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HI),
        (KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 20.0),
        (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, 0.6),
        (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, 80.0),
        ("C1", Pressure::DIMS, 300.0), ("gamma1", Dims::NONE, 250.0),
        ("C2", Pressure::DIMS, 80.0), ("gamma2", Dims::NONE, 30.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: environment.clone(), uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: provenance.clone(),
        }).unwrap();
    }
    let observation = claims.register_observation(ObservationDataset {
        specimen: format!("synthetic Hex8 specimen {region}"),
        method: "synthetic force/projected-area indentation".into(),
        artifact: fs_blake3::hash_domain("fixture", &[region as u8]),
        caveats: "numerical regression only".into(), provenance: provenance.clone(),
    }).unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(HardnessTestContext::new(
            "synthetic spherical indenter, radius 1 mm",
            vec![HardnessLoadStep::new(
                QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
            ).unwrap()],
            "synthetic projected-area protocol revision 1", observation,
        ).unwrap()).unwrap();
    let base = 100.0 + 400.0 * region as f64;
    claims.insert_claim(PropertyClaim {
        key: key.clone(),
        value: PropertyValue::Curve {
            abscissa: ALPHA.into(), abscissa_dims: Dims::NONE,
            knots: vec![(0.0, base), (0.08, base + 0.08 * SLOPE)], dims: Pressure::DIMS,
        },
        validity: environment.with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.08),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![observation], provenance,
    }).unwrap();
    (MaterialCard::assemble(MaterialStateId {
        chemistry: format!("synthetic Hex8 region {region}"), phase: "solid".into(),
        process: "synthetic constitutive and indentation tests".into(), revision: 0,
    }, claims, vec![]).unwrap(), key)
}

fn source(card: &MaterialCard, key: &PropertyKey) -> ResolvedJ2Hardness {
    let native = resolve_j2_plasticity(
        card, &QueryPoint::new().with("T", 300.0).unwrap(), SELECT,
        HardeningPropertyBasis::PlasticStrainModulus, 0.05,
    ).unwrap();
    ResolvedJ2Hardness::new(native, card, key.clone(), SELECT, 0.08).unwrap()
}

fn prescribed(scale: f64) -> Vec<(usize, f64)> {
    // Trilinear native geometry represents this bilinear displacement exactly.
    // Strain and accumulated plastic history differ within a single cell.
    NODES.iter().flat_map(|x| [scale * x[0] * (1.0 + 0.4 * x[1]), 0.0, 0.0])
        .enumerate().collect()
}

fn problem<'a, L: SmallStrainLaw>(
    laws: SmallStrainLawField<'a, L>,
    prescribed: &'a [(usize, f64)],
    formulation: Formulation,
) -> SmallStrainHexProblem<'a, L> {
    SmallStrainHexProblem {
        nodes_m: &NODES, hexahedra: &CELLS, laws, formulation,
        strain_limit: 0.03, rotation_limit_rad: 0.03, prescribed_m: prescribed,
        nodal_forces_n: &[], load_path: &[1.0], budget: HexAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn native_hex_quadrature_hardness_tracks_local_work_hardening_and_reversal() {
    let cards = [card(0), card(1)];
    let sources = [source(&cards[0].0, &cards[0].1), source(&cards[1].0, &cards[1].1)];
    let region = [0, 1, 0, 1, 1, 0, 1, 0];
    let laws = region.map(|r| sources[r].plastic().clone());
    let native_field = [laws.as_slice()];
    let source_field = [region.map(|r| &sources[r])];
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let mut accepted: Option<SmallStrainHexSolution<_>> = None;
        let mut previous = Vec::new();
        for (step, scale) in [0.006, 0.00599, -0.006].into_iter().enumerate() {
            let prescribed = prescribed(scale);
            let p = problem(SmallStrainLawField::PerQuadrature(&native_field), &prescribed, formulation);
            let next = with_cx(|cx| match &accepted {
                Some(old) => p.solve_from(old.displacement_m.as_flattened(), &old.quadrature_states, cx),
                None => p.solve(cx),
            }).unwrap();
            let before = next.clone();
            let hardness = with_cx(|cx| p.observe_hardness(
                &next, HexHardnessField::PerQuadrature(&source_field), TOL, cx,
            )).unwrap();
            assert_eq!(hardness.len(), 8);
            for q in 0..8 {
                if step == 0 {
                    // Independent radial-return oracle, including tensor shear.
                    // B-bar changes only trace, leaving this deviatoric norm intact.
                    let gauss = |i: usize| {
                        let sign = if i == 0 { -1.0 } else { 1.0 };
                        0.5 * (1.0 + sign / 3.0_f64.sqrt())
                    };
                    let (x, y) = (gauss(q / 4), gauss((q / 2) % 2));
                    let exx = scale * (1.0 + 0.4 * y);
                    let exy = 0.2 * scale * x;
                    let mu = E / (2.0 * (1.0 + NU));
                    let trial = 2.0 * mu * (exx * exx + 3.0 * exy * exy).sqrt();
                    let alpha = (trial - 1.0 - region[q] as f64).max(0.0) / (3.0 * mu + HI);
                    near(next.quadrature_states[q].alpha, alpha, 1e-12);
                }
                near(hardness[q], 100.0 + 400.0 * region[q] as f64 + SLOPE * next.quadrature_states[q].alpha, 1e-9);
                if step == 1 { near(hardness[q], previous[q], 1e-9); }
                if step == 2 { assert!(hardness[q] > previous[q]); }
            }
            assert_eq!(next, before, "observation must not update any native state");
            assert!(next.quadrature_states.windows(2).any(|pair| (pair[0].alpha - pair[1].alpha).abs() > 1e-4));
            // No fixed-source export can relabel these different original materials.
            assert!(matches!(with_cx(|cx| p.export_hardness_field(&next, &sources[0], TOL, cx)),
                Err(HexHardnessError::Field(HardnessFieldError::Point { point: 1, .. }))));
            previous = hardness;
            accepted = Some(next);
        }
    }
}

#[test]
fn compact_hex_export_preserves_local_hardness_for_wear_and_rejects_late_mismatch() {
    let (card, key) = card(0);
    let source = source(&card, &key);
    let prescribed = prescribed(0.008);
    let mut p = problem(SmallStrainLawField::Uniform(source.plastic()), &prescribed, Formulation::BBar);
    let solution = with_cx(|cx| p.solve(cx)).unwrap();
    let original = solution.clone();
    let field = with_cx(|cx| p.export_hardness_field(&solution, &source, TOL, cx)).unwrap();
    let per_cell = [&source];
    let h = with_cx(|cx| p.observe_hardness(&solution, HexHardnessField::PerElement(&per_cell), TOL, cx)).unwrap();
    assert_eq!(field.len(), 8);
    assert_eq!(field.source_identity(), source.identity());
    assert_eq!(field.hardness_source().key(), &key);
    for q in 0..8 {
        assert_eq!(field.hardness_pa(q), Some(h[q]));
        assert_eq!(field.sample(q).unwrap().alpha(), solution.quadrature_states[q].alpha);
    }
    let interface = InterfaceSystemRef::new(
        "fixture/Hex8->counterface", "fixture/local-pressure-slip", "fixture/independent-H",
        InputAuthority::SyntheticFixture, InterfaceMedium::Dry,
    ).unwrap();
    let wear = ArchardLaw::advance_local_hardness_pressure_weighted_sliding(
        &interface, &mut WearState::default(), (0..8).map(|q| LocalHardnessArchardSample {
            wear_coefficient: 1e-5, pressure_weighted_sliding_j: 1.0 + q as f64,
            hardness_pa: field.hardness_pa(q).unwrap(),
        }),
    ).unwrap().increment_m3;
    let expected = (0..8).map(|q| 1e-5 * (1.0 + q as f64) / h[q]).sum::<f64>();
    near(wear, expected, 1e-18);
    let averaged = 1e-5 * 36.0 / (h.iter().sum::<f64>() / 8.0);
    assert!((wear - averaged).abs() > 1e-4 * wear);

    let limited = ResolvedJ2Hardness::new(
        source.plastic().clone(), &card, key, SELECT, solution.quadrature_states[7].alpha * 0.5,
    ).unwrap();
    let mut selected = [[&source; 8]];
    selected[0][7] = &limited;
    assert!(matches!(with_cx(|cx| p.observe_hardness(&solution, HexHardnessField::PerQuadrature(&selected), TOL, cx)),
        Err(HexHardnessError::J2 { element: 0, point: 7, .. })));
    let mut malformed = solution.clone();
    malformed.element_offsets[1] = 7;
    assert!(matches!(with_cx(|cx| p.export_hardness_field(&malformed, &source, TOL, cx)), Err(HexHardnessError::Invalid(_))));
    p.budget.maximum_quadrature_points = 7;
    assert!(matches!(with_cx(|cx| p.export_hardness_field(&solution, &source, TOL, cx)), Err(HexHardnessError::Invalid(_))));
    p.budget.maximum_quadrature_points = 8;
    assert!(matches!(with_cancelled_cx(true, |cx| p.export_hardness_field(&solution, &source, TOL, cx)),
        Err(HexHardnessError::Cancelled)));
    assert!(matches!(with_cancelled_cx(true, |cx| p.observe_hardness(&solution, HexHardnessField::Uniform(&source), TOL, cx)),
        Err(HexHardnessError::Cancelled)));
    assert!(matches!(with_cx(|cx| p.observe_hardness(&solution, HexHardnessField::PerQuadrature(&[]), TOL, cx)),
        Err(HexHardnessError::Invalid(_))));
    assert_eq!(with_cx(|cx| p.export_hardness_field(&solution, &source, TOL, cx)).unwrap(), field);
    assert_eq!(solution, original);
}

#[test]
fn cyclic_hex_hardness_observes_complete_independent_branches_at_every_point() {
    let (card, key) = card(0);
    let law = resolve_chaboche_branches(
        &card, &QueryPoint::new().with("T", 300.0).unwrap(), SELECT,
        HardeningPropertyBasis::MixedPlasticStrainModuli, 0.05,
        [BackstressProperties::new("C1", "gamma1"), BackstressProperties::new("C2", "gamma2")],
    ).unwrap();
    let source = ResolvedCyclicHardness::new(law, &card, key, SELECT, 0.08).unwrap();
    let mut accepted: Option<SmallStrainHexSolution<_>> = None;
    for scale in [0.006, -0.006, 0.004] {
        let prescribed = prescribed(scale);
        let p = problem(SmallStrainLawField::Uniform(source.plastic()), &prescribed, Formulation::BBar);
        let solution = with_cx(|cx| match &accepted {
            Some(old) => p.solve_from(old.displacement_m.as_flattened(), &old.quadrature_states, cx),
            None => p.solve(cx),
        }).unwrap();
        let before = solution.clone();
        let h = with_cx(|cx| p.observe_cyclic_hardness(&solution, HexHardnessField::Uniform(&source), TOL, cx)).unwrap();
        let field = with_cx(|cx| p.export_cyclic_hardness_field(&solution, &source, TOL, cx)).unwrap();
        for q in 0..8 {
            let state = &solution.quadrature_states[q];
            assert_ne!(state.backstress_pa[0], state.backstress_pa[1]);
            near(h[q], 100.0 + SLOPE * state.plastic.alpha, 1e-9);
            assert_eq!(field.hardness_pa(q), Some(h[q]));
            assert_eq!(field.sample(q).unwrap().alpha(), state.plastic.alpha);
        }
        let mut invalid = solution.clone();
        invalid.quadrature_states[7].backstress_pa[1][0] += 1.0;
        assert!(matches!(with_cx(|cx| p.observe_cyclic_hardness(&invalid, HexHardnessField::Uniform(&source), TOL, cx)),
            Err(HexHardnessError::Cyclic { element: 0, point: 7, .. })));
        assert!(matches!(with_cx(|cx| p.export_cyclic_hardness_field(&invalid, &source, TOL, cx)),
            Err(HexHardnessError::Field(HardnessFieldError::Point { point: 7, .. }))));
        assert_eq!(solution, before);
        accepted = Some(solution);
    }
}
