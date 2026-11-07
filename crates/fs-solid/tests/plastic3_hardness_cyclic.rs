//! Native cyclic equilibrium and local hardness on two independent solid coupons.
//! Synthetic source cards and uniform physical tractions isolate the coupling;
//! these are numerical regressions, not measured material qualification.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::cyclic_history::ResolvedCyclicHardness;
use fs_material::hardness::history_field::HardnessFieldError;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::plastic::chaboche::source::{
    BackstressProperties, ResolvedChabocheBranches, resolve_chaboche_branches,
};
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
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec};
use fs_solid::plastic3::hardness::{
    CyclicHardnessField, CyclicHardnessObservationError, HardnessFieldExportError,
};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const YIELD_TOLERANCE_PA: f64 = 1e-7;
const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [2.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [2.0, 1.0, 0.0],
    [2.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [4, 5, 6, 7]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 9,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn material(region: usize) -> (MaterialCard, PropertyKey) {
    let provenance = Provenance {
        source: format!("synthetic cyclic coupon {region}; no measured data"),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let environment = ValidityDomain::unconstrained().with("T", 300.0, 300.0);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.3),
        (YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0),
        (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 30.0),
        (KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 20.0),
        (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, 0.6),
        (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, 80.0),
        ("C1", Pressure::DIMS, 300.0),
        ("gamma1", Dims::NONE, 250.0),
        ("C2", Pressure::DIMS, 80.0 + 40.0 * region as f64),
        ("gamma2", Dims::NONE, 30.0),
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
            specimen: format!("synthetic cyclic coupon {region}"),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", &[region as u8]),
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
    let base = [100.0, 500.0][region];
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: ALPHA.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, base), (0.08, base + 1600.0)],
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
                chemistry: format!("synthetic cyclic region {region}"),
                phase: "solid".into(),
                process: "synthetic cyclic and indentation tests".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        key,
    )
}

fn source(
    card: &MaterialCard,
    key: &PropertyKey,
) -> ResolvedCyclicHardness<ResolvedChabocheBranches<2>> {
    let native = resolve_chaboche_branches(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        SELECT,
        BASIS,
        0.05,
        [
            BackstressProperties::new("C1", "gamma1"),
            BackstressProperties::new("C2", "gamma2"),
        ],
    )
    .unwrap();
    ResolvedCyclicHardness::new(native, card, key.clone(), SELECT, 0.08).unwrap()
}

fn displacement(strain: Voigt) -> [f64; 12] {
    [
        0.0,
        0.0,
        0.0,
        strain[0],
        0.0,
        0.0,
        2.0 * strain[3],
        strain[1],
        0.0,
        2.0 * strain[5],
        2.0 * strain[4],
        strain[2],
    ]
}

fn supports(count: usize) -> Vec<(usize, f64)> {
    (0..count)
        .flat_map(|element| [0, 1, 2, 4, 5, 8].map(|dof| (12 * element + dof, 0.0)))
        .collect()
}

fn tractions(stresses: &[Voigt]) -> Vec<(usize, f64)> {
    // Actual consistent boundary forces for unit Cartesian gradients and V=1/6.
    stresses
        .iter()
        .enumerate()
        .flat_map(|(element, s)| {
            [
                (3, s[0]),
                (6, s[3]),
                (7, s[1]),
                (9, s[5]),
                (10, s[4]),
                (11, s[2]),
            ]
            .map(|(dof, force)| (12 * element + dof, force / 6.0))
        })
        .collect()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, allowance {tolerance}"
    );
}

#[test]
fn heterogeneous_native_cyclic_equilibrium_preserves_branches_and_local_hardness_atomically() {
    let cards = [material(0), material(1)];
    let sources = [
        source(&cards[0].0, &cards[0].1),
        source(&cards[1].0, &cards[1].1),
    ];
    let laws = [sources[0].plastic().clone(), sources[1].plastic().clone()];
    let source_refs = [&sources[0], &sources[1]];
    let supports = supports(2);
    let initial_strains = [
        [0.007, -0.002, -0.001, 0.0015, 0.0, 0.0],
        [0.003, 0.006, -0.002, -0.004, 0.001, 0.002],
    ];
    let mut states: Vec<_> = laws
        .iter()
        .zip(initial_strains)
        .map(|(law, strain)| {
            law.law()
                .response(&strain, &law.initial_state())
                .unwrap()
                .state
        })
        .collect();
    let mut u: Vec<_> = initial_strains.into_iter().flat_map(displacement).collect();
    for targets in [
        [
            [0.003, 0.006, -0.002, -0.004, 0.001, 0.002],
            [-0.006, 0.001, 0.002, 0.004, -0.002, 0.001],
        ],
        [
            [-0.006, 0.001, 0.002, 0.004, -0.002, 0.001],
            [0.007, -0.002, -0.001, -0.0015, 0.0, 0.0],
        ],
    ] {
        let expected: Vec<_> = (0..2)
            .map(|element| {
                laws[element]
                    .law()
                    .response(&targets[element], &states[element])
                    .unwrap()
            })
            .collect();
        let loads = tractions(
            &expected
                .iter()
                .map(|step| step.stress_pa)
                .collect::<Vec<_>>(),
        );
        let problem = SmallStrainTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m: &supports,
            nodal_forces_n: &loads,
            load_path: &[1.0],
            budget: TetAssemblyBudget::default(),
            settings: IncrementSettings::default(),
        };
        let mechanical = with_cx(|cx| problem.solve_from(&u, &states, cx)).unwrap();
        assert!(mechanical.increments[0].residual_history.len() > 1);
        assert_eq!(mechanical.increments[0].indefinite_fallbacks, 0);
        let before = mechanical.clone();
        let hardness = with_cx(|cx| {
            problem.observe_cyclic_hardness(
                &mechanical,
                CyclicHardnessField::PerElement(&source_refs),
                YIELD_TOLERANCE_PA,
                cx,
            )
        })
        .unwrap();
        for element in 0..2 {
            let state = &mechanical.element_states[element];
            assert!(state.plastic.alpha > states[element].plastic.alpha);
            close(
                state.plastic.alpha,
                expected[element].state.plastic.alpha,
                2e-9,
            );
            close(
                hardness[element],
                [100.0, 500.0][element] + 20_000.0 * state.plastic.alpha,
                1e-9,
            );
            for component in 0..6 {
                close(
                    mechanical.element_strain[element][component],
                    targets[element][component],
                    2e-9,
                );
                close(
                    mechanical.element_stress[element][component],
                    expected[element].stress_pa[component],
                    2e-8,
                );
                for branch in 0..2 {
                    close(
                        state.backstress_pa[branch][component],
                        expected[element].state.backstress_pa[branch][component],
                        2e-8,
                    );
                }
            }
            assert_ne!(state.backstress_pa[0], state.backstress_pa[1]);
            let typed = sources[element]
                .resolve_hardness_at_accepted(
                    &cards[element].0,
                    &laws[element],
                    &mechanical.element_strain[element],
                    state,
                    YIELD_TOLERANCE_PA,
                )
                .unwrap();
            assert_eq!(typed.key(), &cards[element].1);
            close(typed.hardness_pa(), hardness[element], 1e-9);
        }
        assert_ne!(hardness[0], hardness[1]);
        assert!(matches!(
            with_cx(|cx| problem.export_cyclic_hardness_field(
                &mechanical,
                &sources[0],
                YIELD_TOLERANCE_PA,
                cx,
            )),
            Err(HardnessFieldExportError::Field(HardnessFieldError::Point {
                point: 1,
                ..
            }))
        ));
        let limited = ResolvedCyclicHardness::new(
            laws[1].clone(),
            &cards[1].0,
            cards[1].1.clone(),
            SELECT,
            mechanical.element_states[1].plastic.alpha * 0.5,
        )
        .unwrap();
        for bad_sources in [[&sources[0], &limited], [&sources[0], &sources[0]]] {
            assert!(matches!(
                with_cx(|cx| problem.observe_cyclic_hardness(
                    &mechanical,
                    CyclicHardnessField::PerElement(&bad_sources),
                    YIELD_TOLERANCE_PA,
                    cx,
                )),
                Err(CyclicHardnessObservationError::Element { element: 1, .. })
            ));
        }
        assert!(matches!(
            with_cx(|cx| problem.observe_cyclic_hardness(
                &mechanical,
                CyclicHardnessField::PerElement(&source_refs[..1]),
                YIELD_TOLERANCE_PA,
                cx,
            )),
            Err(CyclicHardnessObservationError::Invalid(_))
        ));
        let mut bad_history = mechanical.clone();
        bad_history.element_states[1].backstress_pa[1][0] += 1.0;
        assert!(matches!(
            with_cx(|cx| problem.observe_cyclic_hardness(
                &bad_history,
                CyclicHardnessField::PerElement(&source_refs),
                YIELD_TOLERANCE_PA,
                cx,
            )),
            Err(CyclicHardnessObservationError::Element { element: 1, .. })
        ));
        assert_eq!(mechanical, before);
        assert_eq!(
            with_cx(|cx| problem.observe_cyclic_hardness(
                &mechanical,
                CyclicHardnessField::PerElement(&source_refs),
                YIELD_TOLERANCE_PA,
                cx,
            ))
            .unwrap(),
            hardness
        );
        u = mechanical
            .displacement_m
            .iter()
            .flatten()
            .copied()
            .collect();
        states = mechanical.element_states;
    }
}

#[test]
fn uniform_voce_solid_observes_its_native_saturated_endpoint() {
    let (card, key) = material(0);
    let law = resolve_voce_j2(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        SELECT,
        BASIS,
        0.05,
    )
    .unwrap();
    let source = ResolvedCyclicHardness::new(law, &card, key, SELECT, 0.08).unwrap();
    let target = [0.008, 0.0, 0.0, 0.001, 0.0, 0.0];
    let old = source.plastic().initial_state();
    let stress = source.plastic().stress(&target, &old);
    let loads = tractions(&[stress]);
    let supports = supports(1);
    let problem = SmallStrainTetProblem {
        nodes_m: &NODES[..4],
        tetrahedra: &TETS[..1],
        laws: SmallStrainLawField::Uniform(source.plastic()),
        prescribed_m: &supports,
        nodal_forces_n: &loads,
        load_path: &[1.0],
        budget: TetAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    };
    let solution = with_cx(|cx| problem.solve(cx)).unwrap();
    let field = with_cx(|cx| {
        problem.observe_cyclic_hardness(
            &solution,
            CyclicHardnessField::Uniform(&source),
            YIELD_TOLERANCE_PA,
            cx,
        )
    })
    .unwrap();
    let exported = with_cx(|cx| {
        problem.export_cyclic_hardness_field(&solution, &source, YIELD_TOLERANCE_PA, cx)
    })
    .unwrap();
    assert_eq!(exported.source_identity(), source.identity());
    assert_eq!(exported.hardness_pa(0), Some(field[0]));
    assert_eq!(
        exported.sample(0).unwrap().alpha(),
        solution.element_states[0].alpha
    );
    for component in 0..6 {
        close(
            solution.element_strain[0][component],
            target[component],
            2e-9,
        );
    }
    assert!(solution.element_states[0].alpha > 0.0);
    close(
        field[0],
        100.0 + 20_000.0 * solution.element_states[0].alpha,
        1e-9,
    );
    assert!(
        source
            .plastic()
            .backbone()
            .law()
            .validate_committed_history(
                &solution.element_strain[0],
                &solution.element_states[0],
                YIELD_TOLERANCE_PA,
            )
            .is_err()
    );
}

#[test]
fn native_chaboche_export_keeps_independent_spatial_histories_under_one_fixed_source() {
    let (card, key) = material(0);
    let source = source(&card, &key);
    let native = source.plastic();
    let initial = [
        [0.007, -0.002, -0.001, 0.0015, 0.0, 0.0],
        [0.003, 0.006, -0.002, -0.004, 0.001, 0.002],
    ];
    let targets = [
        [0.003, 0.006, -0.002, -0.004, 0.001, 0.002],
        [-0.006, 0.001, 0.002, 0.004, -0.002, 0.001],
    ];
    let old: Vec<_> = initial
        .iter()
        .map(|strain| {
            native
                .law()
                .response(strain, &native.initial_state())
                .unwrap()
                .state
        })
        .collect();
    let expected: Vec<_> = (0..2)
        .map(|point| native.law().response(&targets[point], &old[point]).unwrap())
        .collect();
    let loads = tractions(
        &expected
            .iter()
            .map(|point| point.stress_pa)
            .collect::<Vec<_>>(),
    );
    let supports = supports(2);
    let problem = SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(native),
        prescribed_m: &supports,
        nodal_forces_n: &loads,
        load_path: &[1.0],
        budget: TetAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    };
    let u: Vec<_> = initial.into_iter().flat_map(displacement).collect();
    let mechanical = with_cx(|cx| problem.solve_from(&u, &old, cx)).unwrap();
    let before = mechanical.clone();
    let field = with_cx(|cx| {
        problem.export_cyclic_hardness_field(&mechanical, &source, YIELD_TOLERANCE_PA, cx)
    })
    .unwrap();
    assert_eq!(field.len(), 2);
    assert_ne!(field.hardness_pa(0), field.hardness_pa(1));
    assert_eq!(field.source_identity(), source.identity());
    assert_eq!(field.mechanical_source(), native.material());
    assert_eq!(field.hardness_source().key(), &key);
    for point in 0..2 {
        let alpha = mechanical.element_states[point].plastic.alpha;
        assert_eq!(field.sample(point).unwrap().alpha(), alpha);
        close(alpha, expected[point].state.plastic.alpha, 2e-9);
        close(
            field.hardness_pa(point).unwrap(),
            100.0 + 20_000.0 * alpha,
            1e-9,
        );
        for branch in 0..2 {
            for component in 0..6 {
                close(
                    mechanical.element_states[point].backstress_pa[branch][component],
                    expected[point].state.backstress_pa[branch][component],
                    2e-8,
                );
            }
        }
    }
    let mut corrupt = mechanical.clone();
    corrupt.element_states[1].backstress_pa[1][0] += 1.0;
    assert!(matches!(
        with_cx(|cx| problem.export_cyclic_hardness_field(
            &corrupt,
            &source,
            YIELD_TOLERANCE_PA,
            cx,
        )),
        Err(HardnessFieldExportError::Field(HardnessFieldError::Point {
            point: 1,
            ..
        }))
    ));
    assert_eq!(mechanical, before);
}
