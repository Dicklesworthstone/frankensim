//! Actual nodal temperatures, original material sources and native mechanics.
//! Synthetic analytic laws; these are numerical tests, not material validation.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solid::linear3::hex::{HexConnectivity, HexElasticProblem, hex20_point, hex27_point};
use fs_solid::linear3::source::{
    NodalThermalTemperatures, QuadraticThermalSource, QuadraticThermalSources,
    ResolvedQuadraticThermalStrains, SolidMaterialError, ThermalExpansionModel,
    ThermalExpansionSource,
};
use fs_solid::linear3::tet10::{QuadraticTetProblem, tet10_point};
use fs_solid::linear3::{
    QuadraticThermalResponse, QuadraticThermalStrainField, TetAssemblyBudget, TetElasticError,
    TetElasticMaterial, TetMaterialField, TetStaticSolveConfig, TetThermalLoad,
    TetThermalStrainField,
};
use tet_support::{block, material, promote, with_cx};

const FRAME: ContentHash = ContentHash([184; 32]);
const ROTATION: [[f64; 3]; 3] = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn context(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
        .with("p", 125_000.0)
        .unwrap()
}
fn card(isotropic: bool, slope: f64, upper: f64) -> MaterialCard {
    let coefficients = if isotropic {
        vec![(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, 1e-5)]
    } else {
        ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip([1e-5, 0.0, 0.0])
            .collect()
    };
    let mut claims = ClaimSet::new();
    for (name, alpha) in coefficients {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, INVERSE_TEMPERATURE_DIMS),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [290.0, 350.0, upper]
                        .map(|t| (t, alpha + slope * (t - 300.0)))
                        .to_vec(),
                    dims: INVERSE_TEMPERATURE_DIMS,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 290.0, upper)
                    .with("p", 1e5, 2e5),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic nodal-temperature expansion law".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal solid".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn reference(
    card: &MaterialCard,
    model: ThermalExpansionModel,
    law: &TetElasticMaterial,
) -> ThermalExpansionSource {
    model
        .resolve(
            card,
            &context(300.0),
            &context(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            law,
            FRAME,
            0.01,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

enum Problem<'a> {
    Tet(QuadraticTetProblem<'a>),
    Hex(HexElasticProblem<'a>),
}
impl Problem<'_> {
    fn nodes(&self) -> &[[f64; 3]] {
        match self {
            Self::Tet(p) => p.nodes_m,
            Self::Hex(p) => p.nodes_m,
        }
    }
    fn resolve(
        &self,
        temperatures: NodalThermalTemperatures<'_>,
        sources: QuadraticThermalSources<'_>,
        cx: &Cx<'_>,
    ) -> Result<ResolvedQuadraticThermalStrains, SolidMaterialError> {
        match self {
            Self::Tet(p) => p.resolve_nodal_thermal_strains(temperatures, sources, FRAME, 0.01, cx),
            Self::Hex(p) => p.resolve_nodal_thermal_strains(temperatures, sources, FRAME, 0.01, cx),
        }
    }
    fn response(
        &self,
        field: &ResolvedQuadraticThermalStrains,
        cx: &Cx<'_>,
    ) -> QuadraticThermalResponse {
        match self {
            Self::Tet(p) => {
                p.solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx)
            }
            Self::Hex(p) => {
                p.solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx)
            }
        }
        .unwrap()
    }
    fn load<'t>(
        &self,
        field: impl Into<QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> TetThermalLoad {
        let field = field.into();
        match self {
            Self::Tet(p) => p.assemble_thermal_load(field, cx),
            Self::Hex(p) => p.assemble_thermal_load(field, cx),
        }
        .unwrap()
    }
    fn positions(&self) -> Vec<[f64; 3]> {
        match self {
            Self::Tet(p) => {
                let rule = duffy_quadrature(p.quadrature.stiffness_points_per_dir);
                p.elements
                    .iter()
                    .enumerate()
                    .flat_map(|(element, nodes)| {
                        rule.iter()
                            .enumerate()
                            .map(move |(point, &(coordinate, _))| {
                                tet10_point(
                                    nodes.map(|n| p.nodes_m[n]),
                                    coordinate,
                                    1e-12,
                                    element,
                                    point,
                                )
                                .unwrap()
                                .position_m
                            })
                    })
                    .collect()
            }
            Self::Hex(p) => {
                let rule = p.quadrature.points().unwrap();
                match p.elements {
                    HexConnectivity::Hex20(elements) => elements
                        .iter()
                        .enumerate()
                        .flat_map(|(element, nodes)| {
                            rule.iter()
                                .enumerate()
                                .map(move |(point, &(coordinate, _))| {
                                    hex20_point(
                                        nodes.map(|n| p.nodes_m[n]),
                                        coordinate,
                                        1e-12,
                                        element,
                                        point,
                                    )
                                    .unwrap()
                                    .position_m
                                })
                        })
                        .collect(),
                    HexConnectivity::Hex27(elements) => elements
                        .iter()
                        .enumerate()
                        .flat_map(|(element, nodes)| {
                            rule.iter()
                                .enumerate()
                                .map(move |(point, &(coordinate, _))| {
                                    hex27_point(
                                        nodes.map(|n| p.nodes_m[n]),
                                        coordinate,
                                        1e-12,
                                        element,
                                        point,
                                    )
                                    .unwrap()
                                    .position_m
                                })
                        })
                        .collect(),
                }
            }
        }
    }
}
fn fixed_face(nodes: &[[f64; 3]], rotated: bool) -> Vec<usize> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[usize::from(rotated)] == 0.0)
        .flat_map(|(node, _)| [3 * node, 3 * node + 1, 3 * node + 2])
        .collect()
}
fn each_problem(
    law: &TetElasticMaterial,
    rotated: bool,
    cx: &Cx<'_>,
    mut check: impl FnMut(Problem<'_>),
) {
    let rotate = |nodes: &mut [[f64; 3]]| {
        if rotated {
            for p in nodes {
                *p = [-p[1], p[0], p[2]];
            }
        }
    };
    let (mut nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
    rotate(&mut nodes);
    let fixed = fixed_face(&nodes, rotated);
    check(Problem::Tet(QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &elements,
        materials: TetMaterialField::Uniform(law),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: tet_support::Q,
    }));
    for family in hex_support::FAMILIES {
        let mut mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
        rotate(&mut mesh.nodes_m);
        let fixed = fixed_face(&mesh.nodes_m, rotated);
        check(Problem::Hex(hex_support::problem(&mesh, law, &fixed)));
    }
}

#[test]
fn g1_temperature_is_interpolated_before_nonconstant_expansion_is_integrated() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let card = card(true, 1e-7, 450.0);
        let reference = reference(&card, ThermalExpansionModel::Isotropic, &law);
        let source = QuadraticThermalSource {
            card: &card,
            reference_expansion: reference.reference(),
            model: ThermalExpansionModel::Isotropic,
        };
        each_problem(&law, false, cx, |problem| {
            let original: Vec<_> = problem.nodes().iter().map(|p| 310.0 + 5.0 * p[2]).collect();
            let current: Vec<_> = problem
                .nodes()
                .iter()
                .map(|p| 340.0 + 40.0 * p[2])
                .collect();
            let resolved = problem
                .resolve(
                    NodalThermalTemperatures {
                        reference_k: &original,
                        current_k: &current,
                    },
                    QuadraticThermalSources::Uniform(&source),
                    cx,
                )
                .unwrap();
            assert_eq!(resolved.points_per_element(), 27);
            assert_eq!(resolved.strains().len(), resolved.sources().len());
            for ((strain, receipt), position) in resolved
                .strains()
                .iter()
                .zip(resolved.sources())
                .zip(problem.positions())
            {
                let from = 310.0 + 5.0 * position[2];
                let to = 340.0 + 40.0 * position[2];
                let primitive = |t: f64| 1e-5 * (t - 300.0) + 0.5e-7 * (t - 300.0).powi(2);
                let expected = primitive(to) - primitive(from);
                close(strain.reference_temperature_k(), from, 2e-12);
                close(strain.temperature_k(), to, 2e-12);
                for axis in 0..3 {
                    close(strain.free_strain_mandel()[axis], expected, 1e-15);
                }
                assert!((expected - (1e-5 + 1e-7 * (to - 300.0)) * (to - from)).abs() > 1e-5);
                assert_eq!(receipt.current().card_identity(), card.content_hash());
                assert_eq!(receipt.current().query().axes()["p"], 125_000.0);
                assert_eq!(
                    receipt.current().query().axis_quantities()["T"],
                    temperature_quantity()
                );
                assert_eq!(
                    receipt.current().properties()[0].answer().receipt.selected,
                    reference.reference().properties()[0]
                        .answer()
                        .receipt
                        .selected
                );
            }
        });
    });
}

#[test]
fn g1_g3_nodal_gradient_and_rotated_source_axes_produce_actual_compatible_bending() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let card = card(false, 0.0, 450.0);
        for rotated in [false, true] {
            let model = ThermalExpansionModel::OrientedOrthotropic {
                principal_to_mesh: if rotated { ROTATION } else { IDENTITY },
                principal_frame: ContentHash([196; 32]),
            };
            let reference = reference(&card, model, &law);
            let source = QuadraticThermalSource {
                card: &card,
                reference_expansion: reference.reference(),
                model,
            };
            each_problem(&law, rotated, cx, |problem| {
                let original = vec![300.0; problem.nodes().len()];
                let current: Vec<_> = problem
                    .nodes()
                    .iter()
                    .map(|p| 300.0 + 100.0 * p[2])
                    .collect();
                let resolved = problem
                    .resolve(
                        NodalThermalTemperatures {
                            reference_k: &original,
                            current_k: &current,
                        },
                        QuadraticThermalSources::Uniform(&source),
                        cx,
                    )
                    .unwrap();
                let response = problem.response(&resolved, cx);
                for (p, displacement) in problem.nodes().iter().zip(response.displacement_m()) {
                    let x = p[usize::from(rotated)];
                    let axial = 0.001 * x * p[2];
                    let expected = if rotated {
                        [0.0, axial, -0.0005 * x * x]
                    } else {
                        [axial, 0.0, -0.0005 * x * x]
                    };
                    for (actual, expected) in displacement.iter().zip(expected) {
                        close(*actual, expected, 2e-11);
                    }
                }
                assert!(response.elastic_energy_j() < 1e-14);
                assert!(response.elements().iter().all(|element| {
                    element
                        .mean_stress_mandel_pa()
                        .iter()
                        .all(|stress| stress.abs() < 1e-5)
                }));
                assert!(
                    resolved
                        .sources()
                        .iter()
                        .all(|source| source.transform().is_some())
                );
            });
        }
    });
}

#[test]
fn g3_uniform_temperature_matches_the_original_source_integrated_load_exactly() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let card = card(true, 1e-7, 450.0);
        let model = ThermalExpansionModel::Isotropic;
        let reference = reference(&card, model, &law);
        let source = QuadraticThermalSource {
            card: &card,
            reference_expansion: reference.reference(),
            model,
        };
        let selected = reference.reference().properties()[0]
            .answer()
            .receipt
            .selected;
        let uniform = model
            .resolve(
                &card,
                &context(300.0),
                &context(340.0),
                MaterialPropertySelection::PinnedByProperty(vec![(
                    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.into(),
                    selected,
                )]),
                &law,
                FRAME,
                0.01,
            )
            .unwrap();
        each_problem(&law, false, cx, |problem| {
            let original = vec![300.0; problem.nodes().len()];
            let current = vec![340.0; problem.nodes().len()];
            let resolved = problem
                .resolve(
                    NodalThermalTemperatures {
                        reference_k: &original,
                        current_k: &current,
                    },
                    QuadraticThermalSources::Uniform(&source),
                    cx,
                )
                .unwrap();
            assert!(resolved.sources().iter().all(|source| source == &uniform));
            let old = problem.load(TetThermalStrainField::Uniform(uniform.thermal_strain()), cx);
            let current = problem.load(&resolved, cx);
            assert_eq!(current.full_equivalent_force_n, old.full_equivalent_force_n);
            assert_eq!(current.reduced_force_n, old.reduced_force_n);
        });
    });
}

#[test]
fn g0_source_mismatch_counts_and_temperature_domain_refuse_without_mutation() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let card = card(true, 1e-7, 450.0);
        let other = self::card(true, 2e-7, 450.0);
        let reference = reference(&card, ThermalExpansionModel::Isotropic, &law);
        let source = QuadraticThermalSource {
            card: &card,
            reference_expansion: reference.reference(),
            model: ThermalExpansionModel::Isotropic,
        };
        each_problem(&law, false, cx, |problem| {
            let original = vec![300.0; problem.nodes().len()];
            let current = vec![340.0; problem.nodes().len()];
            let values = NodalThermalTemperatures {
                reference_k: &original,
                current_k: &current,
            };
            let before = problem
                .resolve(values, QuadraticThermalSources::Uniform(&source), cx)
                .unwrap();
            let bad_source = QuadraticThermalSource {
                card: &other,
                ..source
            };
            assert!(matches!(
                problem.resolve(values, QuadraticThermalSources::Uniform(&bad_source), cx),
                Err(SolidMaterialError::Solid(
                    TetElasticError::InvalidThermalStrain { .. }
                ))
            ));
            assert!(
                problem
                    .resolve(values, QuadraticThermalSources::PerElement(&[]), cx)
                    .is_err()
            );
            assert!(
                problem
                    .resolve(
                        NodalThermalTemperatures {
                            reference_k: &original[..original.len() - 1],
                            current_k: &current,
                        },
                        QuadraticThermalSources::Uniform(&source),
                        cx
                    )
                    .is_err()
            );
            for bad in [0.0, -1.0, f64::NAN, 500.0] {
                let bad_current = vec![bad; current.len()];
                assert!(
                    problem
                        .resolve(
                            NodalThermalTemperatures {
                                reference_k: &original,
                                current_k: &bad_current,
                            },
                            QuadraticThermalSources::Uniform(&source),
                            cx
                        )
                        .is_err()
                );
            }
            let after = problem
                .resolve(values, QuadraticThermalSources::Uniform(&source), cx)
                .unwrap();
            assert_eq!(before.strains(), after.strains());
            assert_eq!(before.sources(), after.sources());
        });
    });
}

#[test]
fn g0_positive_in_domain_nodes_do_not_hide_high_order_temperature_overshoot() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let card = card(true, 0.0, 400.0);
        let reference = reference(&card, ThermalExpansionModel::Isotropic, &law);
        let source = QuadraticThermalSource {
            card: &card,
            reference_expansion: reference.reference(),
            model: ThermalExpansionModel::Isotropic,
        };
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let problem = hex_support::problem(&mesh, &law, &[]);
            let original = vec![300.0; mesh.nodes_m.len()];
            // Q2 interpolation of these valid nodal values reaches 412.5 K.
            let current: Vec<_> = mesh
                .nodes_m
                .iter()
                .map(|p| 300.0 + 300.0 * p[2] - 200.0 * p[2] * p[2])
                .collect();
            assert!(current.iter().all(|t| (300.0..=400.0).contains(t)));
            assert!(matches!(
                problem.resolve_nodal_thermal_strains(
                    NodalThermalTemperatures {
                        reference_k: &original,
                        current_k: &current,
                    },
                    QuadraticThermalSources::Uniform(&source),
                    FRAME,
                    0.01,
                    cx
                ),
                Err(SolidMaterialError::Material(_))
            ));
        }
    });
}
