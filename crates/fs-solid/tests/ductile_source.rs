//! G1/G3: original material-card ductility changes native tet acceptance.
//! Synthetic coefficients test the integration, not a calibrated real alloy.
#![deny(unsafe_code)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy,
    LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::ductile::source::{
    DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY, LODE_PARAMETER_AXIS,
    JOHNSON_COOK_INITIATION_LAW_ID, PLASTIC_STRAIN_RATE_AXIS, STRESS_TRIAXIALITY_AXIS,
    resolve_ductile_locus, resolve_ductile_stress_locus, resolve_johnson_cook_initiation,
};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ResolvedJ2Plasticity, resolve_j2_plasticity,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::ductile::{
    DuctileAdvanceError, DuctileHistory, DuctileLawField, DuctileSolidError, DuctileStepConditions,
    DuctileTemperatures,
};
use fs_solid::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget};

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn fixed_point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), 300.0)
        .unwrap()
}
fn calibration_point() -> QueryPoint {
    fixed_point()
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
        )
        .unwrap()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            0.01,
        )
        .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic native ductility regression".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn material(fracture_strain: f64) -> MaterialCard {
    let fixed = ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 300.0);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, 1000.0),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("yield_stress", Pressure::DIMS, 1.0),
        ("isotropic_hardening_modulus", Pressure::DIMS, 100.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: fixed.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: STRESS_TRIAXIALITY_AXIS.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(-1.0, fracture_strain), (1.0, fracture_strain)],
                dims: Dims::NONE,
            },
            validity: fixed
                .with_quantity(
                    STRESS_TRIAXIALITY_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    -1.0,
                    1.0,
                )
                .with_quantity(
                    PLASTIC_STRAIN_RATE_AXIS,
                    QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                    1e-8,
                    1.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::TabulatedOnly,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic equal-stiffness equal-yield alloy".into(),
            phase: "solid".into(),
            process: "synthetic coupon".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(card: &MaterialCard) -> ResolvedJ2Plasticity {
    resolve_j2_plasticity(
        card,
        &fixed_point(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::PlasticStrainModulus,
        0.05,
    )
    .unwrap()
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 5,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const TARGET: [f64; 1] = [0.006];
fn problem<'a>(
    law: &'a ResolvedJ2Plasticity,
    prescribed: &'a [(usize, f64)],
) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: prescribed,
        nodal_forces_n: &[],
        load_path: &TARGET,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn sourced_ductility_changes_acceptance_and_source_swaps_refuse_transactionally() {
    let high_card = material(0.02);
    let low_card = material(0.003);
    let high_native = resolve(&high_card);
    let low_native = resolve(&low_card);
    assert_eq!(
        high_native.law(),
        low_native.law(),
        "only fracture calibration differs"
    );
    let request = DuctileSourceRequest {
        temperature_axis: "T".into(),
        maximum_points: 2,
    };
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let high =
        resolve_ductile_locus(&high_card, &calibration_point(), &request, selection).unwrap();
    let low = resolve_ductile_locus(&low_card, &calibration_point(), &request, selection).unwrap();
    // Prescribe axial displacement; two lateral DOFs remain free for actual equilibrium.
    let mut prescribed = Vec::new();
    for (node, x) in NODES.iter().enumerate() {
        prescribed.push((3 * node, x[0]));
        if x[1] == 0.0 {
            prescribed.push((3 * node + 1, 0.0));
        }
        if x[2] == 0.0 {
            prescribed.push((3 * node + 2, 0.0));
        }
    }
    let high_problem = problem(&high_native, &prescribed);
    let low_problem = problem(&low_native, &prescribed);
    let conditions = DuctileStepConditions {
        duration_s: 1.0,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    };
    with_cx(|cx| {
        let high_history = DuctileHistory::from_source_bound_j2_tet(
            &high_problem,
            DuctileLawField::Uniform(&high),
            cx,
        )
        .unwrap();
        assert!(matches!(
            DuctileHistory::from_source_bound_j2_tet(
                &high_problem,
                DuctileLawField::Uniform(&low),
                cx
            ),
            Err(DuctileSolidError::Input(_))
        ));
        let virgin = [high_native.initial_state()];
        let accepted = high_history
            .advance_source_bound_j2_tet(&high_problem, &[0.0; 12], &virgin, conditions, cx)
            .unwrap();
        // Independent monotonic uniaxial equations: alpha=(E*epsilon-sigma_y)/(E+H).
        let expected_alpha = (1000.0 * TARGET[0] - 1.0) / 1100.0;
        assert!((accepted.mechanics.element_states[0].alpha - expected_alpha).abs() < 1e-10);
        assert!(
            (accepted.mechanics.element_stress[0][0] - (1.0 + 100.0 * expected_alpha)).abs() < 1e-8
        );
        assert!(
            (accepted.ductility.points()[0].utilization() - expected_alpha / 0.02).abs() < 1e-8
        );
        let low_history = DuctileHistory::from_source_bound_j2_tet(
            &low_problem,
            DuctileLawField::Uniform(&low),
            cx,
        )
        .unwrap();
        match low_history.advance_source_bound_j2_tet(
            &low_problem,
            &[0.0; 12],
            &[low_native.initial_state()],
            conditions,
            cx,
        ) {
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Initiated {
                point,
                candidate,
            })) => {
                assert_eq!(point, 0);
                assert!(
                    (candidate.state.accumulated_plastic_strain() - expected_alpha).abs() < 1e-10
                );
                assert!(candidate.state.utilization() > 1.0);
            }
            other => panic!("expected initiation from the lower sourced ductility: {other:?}"),
        }
        assert_eq!(low_history.points()[0].utilization(), 0.0);
        let displacement: Vec<f64> = accepted
            .mechanics
            .displacement_m
            .iter()
            .flatten()
            .copied()
            .collect();
        let saved = accepted.ductility.points().to_vec();
        // Even an already accepted run must reject a later card swap before solving.
        assert!(matches!(
            accepted.ductility.advance_source_bound_j2_tet(
                &low_problem,
                &displacement,
                &accepted.mechanics.element_states,
                conditions,
                cx
            ),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))
        ));
        assert_eq!(accepted.ductility.points(), saved.as_slice());
        assert_eq!(high_history.points()[0].utilization(), 0.0);
    });
}

#[test]
fn original_lode_surface_reaches_the_assembled_tet_history() {
    let base = material(0.02);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut row = original.clone();
        if row.key.name() == FRACTURE_PLASTIC_STRAIN_PROPERTY {
            row.value = PropertyValue::Surface2D {
                axes: [STRESS_TRIAXIALITY_AXIS.into(), LODE_PARAMETER_AXIS.into()],
                axis_dims: [Dims::NONE; 2],
                samples: [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]
                    .into_iter()
                    .map(|at| (at, 0.02 + 0.01 * at[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [0, 2, 3]],
                dims: Dims::NONE,
            };
            row.validity = row.validity.with_quantity(
                LODE_PARAMETER_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                -1.0,
                1.0,
            );
            row.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(row).unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let native = resolve(&card);
    let locus = resolve_ductile_stress_locus(
        &card,
        &calibration_point()
            .with_quantity(
                LODE_PARAMETER_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                -1.0,
            )
            .unwrap(),
        &DuctileSourceRequest {
            temperature_axis: "T".into(),
            maximum_points: 4,
        },
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap();
    // The two affine gradients have equal trace and equivalent strain. Their
    // native J2 accumulated plastic strain agrees, while xi is 0 versus +1.
    let a = 2.0 / 3.0_f64.sqrt();
    let shear: Vec<_> = NODES
        .iter()
        .enumerate()
        .flat_map(|(node, p)| {
            [
                (3 * node, 2.0 * p[1]),
                (3 * node + 1, 0.0),
                (3 * node + 2, 0.0),
            ]
        })
        .collect();
    let meridian: Vec<_> = NODES
        .iter()
        .enumerate()
        .flat_map(|(node, p)| {
            [
                (3 * node, a * p[0]),
                (3 * node + 1, -0.5 * a * p[1]),
                (3 * node + 2, -0.5 * a * p[2]),
            ]
        })
        .collect();
    with_cx(|cx| {
        let conditions = DuctileStepConditions {
            duration_s: 1.0,
            temperatures_k: DuctileTemperatures::Uniform(300.0),
        };
        let run = |supports: &[(usize, f64)]| {
            let specimen = problem(&native, supports);
            let history = DuctileHistory::from_source_bound_j2_tet(
                &specimen,
                DuctileLawField::Uniform(&locus),
                cx,
            )
            .unwrap();
            history
                .advance_source_bound_j2_tet(
                    &specimen,
                    &[0.0; 12],
                    &[native.initial_state()],
                    conditions,
                    cx,
                )
                .unwrap()
        };
        let shear = run(&shear);
        let meridian = run(&meridian);
        let shear_alpha = shear.mechanics.element_states[0].alpha;
        let meridian_alpha = meridian.mechanics.element_states[0].alpha;
        assert!((shear_alpha - meridian_alpha).abs() < 1e-12);
        let mu = 1000.0 / 2.6;
        let expected_alpha = (2.0 * mu * 3.0_f64.sqrt() * TARGET[0] - 1.0) / (3.0 * mu + 100.0);
        assert!((shear_alpha - expected_alpha).abs() < 1e-12);
        assert!((shear.ductility.points()[0].utilization() - expected_alpha / 0.02).abs() < 1e-10);
        assert!(
            (meridian.ductility.points()[0].utilization() - expected_alpha / 0.03).abs() < 1e-10
        );
    });
}

fn johnson_cook_card() -> MaterialCard {
    let base = material(0.02);
    let thermal = ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 600.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        if original.key.name() != FRACTURE_PLASTIC_STRAIN_PROPERTY {
            let mut claim = original.clone();
            // Deliberately constant synthetic mechanical coefficients isolate
            // the temperature/rate effects in the initiation calibration.
            claim.validity = thermal.clone();
            claims.insert_claim(claim).unwrap();
        }
    }
    let model = ConstitutiveModelCard {
        law: LawId(JOHNSON_COOK_INITIATION_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("d1", 0.006, Dims::NONE),
            ("d2", 0.0, Dims::NONE),
            ("d3", 0.0, Dims::NONE),
            ("d4", 0.04, Dims::NONE),
            ("d5", -0.5, Dims::NONE),
            ("reference_strain_rate_s_inv", 1.0, Dims([0, 0, -1, 0, 0, 0])),
            ("transition_temperature_k", 300.0, Temperature::DIMS),
            ("melting_temperature_k", 900.0, Temperature::DIMS),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: thermal
            .with_quantity(STRESS_TRIAXIALITY_AXIS, QuantitySpec::dimensional(Dims::NONE), -1.0, 1.0)
            .with_quantity(
                PLASTIC_STRAIN_RATE_AXIS,
                QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                1e-8, 1.0,
            ),
        sources: vec![fs_blake3::hash_bytes(b"synthetic native Johnson-Cook calibration")],
        provenance: provenance(),
    };
    MaterialCard::assemble(base.id().clone(), claims, vec![model]).unwrap()
}

fn resolve_at(card: &MaterialCard, temperature_k: f64) -> ResolvedJ2Plasticity {
    resolve_j2_plasticity(
        card,
        &QueryPoint::new().with_quantity("T", temperature(), temperature_k).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::PlasticStrainModulus,
        0.05,
    ).unwrap()
}

fn axial_supports() -> Vec<(usize, f64)> {
    NODES.iter().enumerate().flat_map(|(node, x)| {
        let mut supports = vec![(3 * node, x[0])];
        if x[1] == 0.0 { supports.push((3 * node + 1, 0.0)); }
        if x[2] == 0.0 { supports.push((3 * node + 2, 0.0)); }
        supports
    }).collect()
}

#[test]
fn sourced_johnson_cook_changes_native_acceptance_with_rate_and_actual_temperature() {
    let card = johnson_cook_card();
    let locus = resolve_johnson_cook_initiation(
        &card, card.models()[0].content_hash(), &calibration_point(), "T",
    ).unwrap();
    let cold = resolve_at(&card, 300.0);
    let hot = resolve_at(&card, 600.0);
    assert_eq!(cold.law(), hot.law());
    let supports = axial_supports();
    let alpha = (1000.0 * TARGET[0] - 1.0) / 1100.0;
    with_cx(|cx| {
        let specimen = problem(&cold, &supports);
        let history = DuctileHistory::from_source_bound_j2_tet(
            &specimen, DuctileLawField::Uniform(&locus), cx,
        ).unwrap();
        let cold_slow = history.advance_source_bound_j2_tet(
            &specimen, &[0.0; 12], &[cold.initial_state()],
            DuctileStepConditions {
                duration_s: 1.0, temperatures_k: DuctileTemperatures::Uniform(300.0),
            }, cx,
        ).unwrap();
        let cold_fast = history.advance_source_bound_j2_tet(
            &specimen, &[0.0; 12], &[cold.initial_state()],
            DuctileStepConditions {
                duration_s: 0.1, temperatures_k: DuctileTemperatures::Uniform(300.0),
            }, cx,
        ).unwrap();
        for (accepted, dt) in [(&cold_slow, 1.0), (&cold_fast, 0.1)] {
            let expected = alpha / (0.006 * (1.0 + 0.04 * (alpha / dt).ln()));
            assert!((accepted.mechanics.element_states[0].alpha - alpha).abs() < 1e-10);
            assert!((accepted.ductility.points()[0].utilization() - expected).abs() < 1e-8);
        }
        assert_eq!(cold_slow.mechanics, cold_fast.mechanics);
        assert!(cold_fast.ductility.points()[0].utilization() < cold_slow.ductility.points()[0].utilization());
        let hot_specimen = problem(&hot, &supports);
        let refused = history.advance_source_bound_j2_tet(
            &hot_specimen, &[0.0; 12], &[hot.initial_state()],
            DuctileStepConditions {
                duration_s: 1.0, temperatures_k: DuctileTemperatures::Uniform(600.0),
            }, cx,
        ).unwrap_err();
        match refused {
            DuctileAdvanceError::Ductility(DuctileSolidError::Initiated { point, candidate }) => {
                assert_eq!(point, 0);
                let expected = alpha / (0.006 * (1.0 + 0.04 * alpha.ln()) * 0.75);
                assert!((candidate.state.utilization() - expected).abs() < 1e-8);
                assert!(candidate.state.initiated());
            }
            other => panic!("expected physical-temperature initiation, got {other:?}"),
        }
        // Merely changing the criterion temperature must not silently heat a
        // mechanical law frozen at 300 K, even though the fit supports 600 K.
        assert!(matches!(history.advance_source_bound_j2_tet(
            &specimen, &[0.0; 12], &[cold.initial_state()],
            DuctileStepConditions {
                duration_s: 1.0, temperatures_k: DuctileTemperatures::Uniform(600.0),
            }, cx,
        ), Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))));
        assert_eq!(history.points()[0].utilization(), 0.0);
    });
}

#[test]
fn sourced_coupled_candidate_binds_its_own_temperature_and_refuses_unbound_authority() {
    let card = johnson_cook_card();
    let locus = resolve_johnson_cook_initiation(
        &card, card.models()[0].content_hash(), &calibration_point(), "T",
    ).unwrap();
    let cold = resolve_at(&card, 300.0);
    let hot = resolve_at(&card, 450.0);
    let supports = axial_supports();
    with_cx(|cx| {
        let cold_specimen = problem(&cold, &supports);
        // Authored numerical coefficients remain usable through ordinary
        // entry points; this specifically sourced entry cannot label them so.
        assert!(matches!(DuctileHistory::from_source_bound_j2_tet(
            &cold_specimen, DuctileLawField::Uniform(locus.law()), cx,
        ), Err(DuctileSolidError::Input(_))));
        let history = DuctileHistory::from_source_bound_j2_tet(
            &cold_specimen, DuctileLawField::Uniform(&locus), cx,
        ).unwrap();
        let old = [cold.initial_state()];
        let candidate = problem(&hot, &supports).solve_from(&[0.0; 12], &old, cx).unwrap();
        let conditions = DuctileStepConditions {
            duration_s: 0.1, temperatures_k: DuctileTemperatures::Uniform(450.0),
        };
        let gated = history.check_material_tet_candidate(
            &old, core::slice::from_ref(&cold), &candidate, core::slice::from_ref(&hot),
            conditions, cx,
        ).unwrap();
        let alpha = candidate.element_states[0].alpha;
        let expected = alpha / (0.006 * (1.0 + 0.04 * (alpha / 0.1).ln()) * 0.875);
        assert!((gated.points()[0].utilization() - expected).abs() < 1e-12);
        assert!(history.check_material_tet_candidate(
            &old, core::slice::from_ref(&cold), &candidate, core::slice::from_ref(&cold),
            conditions, cx,
        ).is_err());
        assert!(history.validate_material_history(
            &old, core::slice::from_ref(&cold), conditions, cx,
        ).is_err());
        assert_eq!(history.points()[0].utilization(), 0.0);
        assert_eq!(old, [cold.initial_state()]);
    });
}
