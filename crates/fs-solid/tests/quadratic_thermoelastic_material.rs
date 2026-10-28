//! Source-backed C(T) and expansion on the original quadratic solid owners.
//! Independent polynomial work/energy oracles on a unit cube; synthetic data.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solid::linear3::hex::{HexElasticProblem, HexQuadrature};
use fs_solid::linear3::source::*;
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature};
use fs_solid::linear3::{
    QuadraticPrestressedAssembly, QuadraticThermalResponse, QuadraticThermalStrainField,
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial, TetMaterialField,
    TetStaticSolveConfig, TetThermalLoad,
};
use tet_support::{block, promote, quadratic_form, with_cx};

const FRAME: ContentHash = ContentHash([184; 32]);
const YOUNG: f64 = 2e6;
const SLOPE: f64 = -1e4;
const B: f64 = 1e-3;
fn quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn context(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", quantity(), t).unwrap()
}
fn expansion_model() -> ThermalExpansionModel {
    ThermalExpansionModel::OrientedOrthotropic {
        principal_to_mesh: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        principal_frame: ContentHash([197; 32]),
    }
}
fn fixture(
    elastic_upper: f64,
) -> (
    MaterialCard,
    ResolvedMaterialStatePoint,
    ThermalExpansionSource,
    TetElasticMaterial,
) {
    let mut claims = ClaimSet::new();
    let mut properties = vec![
        (
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            YOUNG,
            SLOPE,
            elastic_upper,
        ),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0, 0.0, elastic_upper),
        // Density has NO support at any heated quadrature point. It belongs
        // only to the original reference inventory and must not be re-queried.
        (DENSITY_PROPERTY, Density::DIMS, 1.0, 0.0, 300.0),
    ];
    properties.extend(
        ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip([1e-5, 0.0, 0.0])
            .map(|(name, value)| (name, INVERSE_TEMPERATURE_DIMS, value, 0.0, 450.0)),
    );
    for (name, dims, value, slope, upper) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (290.0, value - 10.0 * slope),
                        (upper, value + (upper - 300.0) * slope),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    290.0,
                    upper,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic local thermoelastic polynomial fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic elastic solid".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let elastic = ThermalStiffnessModel::Isotropic
        .resolve_reference(
            &card,
            &context(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    let law = ThermalStiffnessModel::Isotropic
        .material(&elastic, 1.0, elastic.identity(), 0.01)
        .unwrap();
    let expansion = expansion_model()
        .resolve(
            &card,
            &context(300.0),
            &context(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            &law,
            FRAME,
            0.01,
        )
        .unwrap();
    (card, elastic, expansion, law)
}
fn source<'a>(
    card: &'a MaterialCard,
    elastic: &'a ResolvedMaterialStatePoint,
    expansion: &'a ThermalExpansionSource,
) -> QuadraticThermoelasticSource<'a> {
    QuadraticThermoelasticSource {
        expansion: QuadraticThermalSource {
            card,
            reference_expansion: expansion.reference(),
            model: expansion_model(),
        },
        reference_elastic: elastic,
        stiffness_model: ThermalStiffnessModel::Isotropic,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance={tolerance:e}"
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
        source: &QuadraticThermoelasticSource<'_>,
        cx: &Cx<'_>,
    ) -> Result<ResolvedQuadraticThermalStrains, SolidMaterialError> {
        let original = vec![300.0; self.nodes().len()];
        let current: Vec<_> = self.nodes().iter().map(|p| 300.0 + 100.0 * p[2]).collect();
        let t = NodalThermalTemperatures {
            reference_k: &original,
            current_k: &current,
        };
        let sources = QuadraticThermoelasticSources::Uniform(source);
        match self {
            Self::Tet(p) => p.resolve_nodal_thermoelasticity(t, sources, FRAME, 0.01, cx),
            Self::Hex(p) => p.resolve_nodal_thermoelasticity(t, sources, FRAME, 0.01, cx),
        }
    }
    fn assembly(
        &self,
        field: Option<&ResolvedQuadraticThermalStrains>,
        cx: &Cx<'_>,
    ) -> TetElasticAssembly {
        match (self, field) {
            (Self::Tet(p), Some(f)) => p.assemble_thermoelastic(f, cx),
            (Self::Hex(p), Some(f)) => p.assemble_thermoelastic(f, cx),
            (Self::Tet(p), None) => p.assemble(cx),
            (Self::Hex(p), None) => p.assemble(cx),
        }
        .unwrap()
    }
    fn load(
        &self,
        field: QuadraticThermalStrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<TetThermalLoad, TetElasticError> {
        match self {
            Self::Tet(p) => p.assemble_thermal_load(field, cx),
            Self::Hex(p) => p.assemble_thermal_load(field, cx),
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
    fn tangent(
        &self,
        field: &ResolvedQuadraticThermalStrains,
        forces: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> QuadraticPrestressedAssembly {
        match self {
            Self::Tet(p) => p.solve_thermal_tangent_with_nodal_forces(
                field,
                forces,
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            ),
            Self::Hex(p) => p.solve_thermal_tangent_with_nodal_forces(
                field,
                forces,
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            ),
        }
        .unwrap()
    }
}
fn each_problem(
    law: &TetElasticMaterial,
    fixed_axis: usize,
    cx: &Cx<'_>,
    mut check: impl FnMut(Problem<'_>),
) {
    let fixed = |nodes: &[[f64; 3]]| {
        nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[fixed_axis] == 0.0)
            .flat_map(|(node, _)| [3 * node, 3 * node + 1, 3 * node + 2])
            .collect::<Vec<_>>()
    };
    let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
    let constraints = fixed(&nodes);
    check(Problem::Tet(QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &elements,
        materials: TetMaterialField::Uniform(law),
        fixed_dofs: &constraints,
        budget: TetAssemblyBudget::standard(),
        // Variable C times gradient prestress reaches degree four. The caller
        // deliberately raises the original stiffness rule; no hidden averaging.
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 4,
            mass_points_per_dir: 4,
        },
    }));
    for family in hex_support::FAMILIES {
        let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
        let constraints = fixed(&mesh.nodes_m);
        let mut problem = hex_support::problem(&mesh, law, &constraints);
        problem.quadrature = HexQuadrature { points_per_dir: 4 };
        check(Problem::Hex(problem));
    }
}
fn mode(nodes: &[[f64; 3]], free: &[usize]) -> Vec<f64> {
    free.iter()
        .map(|&dof| {
            if dof % 3 == 0 {
                nodes[dof / 3][0] * nodes[dof / 3][2]
            } else {
                0.0
            }
        })
        .collect()
}
fn coupled<'a>(
    field: &'a ResolvedQuadraticThermalStrains,
    materials: &'a [TetElasticMaterial],
) -> QuadraticThermalStrainField<'a> {
    QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
        states: field.strains(),
        materials,
        points_per_element: field.points_per_element(),
    }
}

#[test]
fn g1_local_elasticity_controls_stiffness_thermal_work_energy_and_prestress_at_fixed_mass() {
    with_cx(|cx| {
        let (card, elastic, expansion, law) = fixture(450.0);
        let source = source(&card, &elastic, &expansion);
        each_problem(&law, 2, cx, |problem| {
            let field = problem.resolve(&source, cx).unwrap();
            let hot = problem.assembly(Some(&field), cx);
            let cold = problem.assembly(None, cx);
            let mode = mode(problem.nodes(), &hot.free_dofs);
            // E(z)=E0(1-z/2), epsilon_free_xx=b*z, trial u=(xz,0,0).
            // Integrals over the unit cube, independent of FE quadrature:
            // u'Ku = E0 integral (1-z/2)(z^2+1/6) dz = E0/3.
            close(quadratic_form(&hot.stiffness, &mode), YOUNG / 3.0, 1e-7);
            close(quadratic_form(&hot.mass, &mode), 1.0 / 9.0, 1e-13);
            close(hot.total_mass_kg, 1.0, 1e-13);
            assert_eq!(hot.total_mass_kg, cold.total_mass_kg);
            for row in 0..hot.free_dofs.len() {
                assert_eq!(hot.mass.row(row), cold.mass.row(row));
            }
            assert!(
                field
                    .elastic_sources()
                    .unwrap()
                    .iter()
                    .all(|s| s.property(DENSITY_PROPERTY).is_none())
            );
            let load = problem.load(field.field(), cx).unwrap();
            let thermal_work: f64 = mode
                .iter()
                .zip(&load.reduced_force_n)
                .map(|(u, f)| u * f)
                .sum();
            close(thermal_work, YOUNG * B * 5.0 / 24.0, 1e-10);
            // Equal opposing DEAD nodal forces hold u=0 exactly. Their actual
            // constitutive stress must also drive the geometric stiffness.
            let forces: Vec<_> = load
                .full_equivalent_force_n
                .iter()
                .enumerate()
                .map(|(dof, f)| (dof, -f))
                .collect();
            let tangent = problem.tangent(&field, &forces, cx);
            assert!(
                tangent
                    .response()
                    .displacement_m()
                    .iter()
                    .flatten()
                    .all(|u| *u == 0.0)
            );
            close(
                tangent.response().elastic_energy_j(),
                YOUNG * B * B * 5.0 / 48.0,
                1e-11,
            );
            close(
                quadratic_form(tangent.elastic_stiffness(), &mode),
                YOUNG / 3.0,
                1e-7,
            );
            close(
                quadratic_form(tangent.geometric_stiffness(), &mode),
                -YOUNG * B * 3.0 / 20.0,
                1e-8,
            );
        });
    });
}

#[test]
fn g1_actual_equilibrium_uses_the_same_softened_operator_as_thermal_loading() {
    with_cx(|cx| {
        let (card, elastic, expansion, law) = fixture(450.0);
        let source = source(&card, &elastic, &expansion);
        each_problem(&law, 0, cx, |problem| {
            let field = problem.resolve(&source, cx).unwrap();
            let response = problem.response(&field, cx);
            // This compatible quadratic expansion is stress free for EVERY
            // positive local E(z). A cold K with a hot thermal load fails it.
            for (p, u) in problem.nodes().iter().zip(response.displacement_m()) {
                for (actual, expected) in
                    u.iter().zip([B * p[0] * p[2], 0.0, -0.5 * B * p[0] * p[0]])
                {
                    close(*actual, expected, 5e-11);
                }
            }
            assert!(response.elastic_energy_j() < 1e-14);
            assert!(
                response
                    .elements()
                    .iter()
                    .all(|e| e.mean_stress_mandel_pa().iter().all(|s| s.abs() < 1e-5))
            );
        });
    });
}

#[test]
fn g0_local_material_counts_density_and_late_identity_mismatch_refuse_and_coefficients_bind_identity()
 {
    with_cx(|cx| {
        let (card, elastic, expansion, law) = fixture(450.0);
        let source = source(&card, &elastic, &expansion);
        each_problem(&law, 2, cx, |problem| {
            let field = problem.resolve(&source, cx).unwrap();
            let original = problem.load(field.field(), cx).unwrap();
            let mut materials = field.elastic_materials().unwrap().to_vec();
            assert!(matches!(
                problem.load(coupled(&field, &materials[..materials.len() - 1]), cx),
                Err(TetElasticError::MaterialCountMismatch { .. })
            ));
            let last = materials.len() - 1;
            materials[last].density_kg_m3 = 2.0;
            assert!(matches!(
                problem.load(coupled(&field, &materials), cx),
                Err(TetElasticError::InvalidMaterial { .. })
            ));
            materials[last] = field.elastic_materials().unwrap()[last].clone();
            materials[last].material_state_identity = [198; 32];
            assert!(matches!(
                problem.load(coupled(&field, &materials), cx),
                Err(TetElasticError::ThermalElasticStateMismatch { .. })
            ));
            // Deliberately retain an input identity while changing only the
            // last local tensor: load identity must bind actual coefficient bits.
            materials[last] = TetElasticMaterial::try_new(
                1.0,
                3e6,
                0.0,
                field.elastic_materials().unwrap()[last].material_state_identity,
            )
            .unwrap();
            assert_ne!(
                problem
                    .load(coupled(&field, &materials), cx)
                    .unwrap()
                    .identity,
                original.identity
            );
            assert_eq!(problem.load(field.field(), cx).unwrap(), original);
        });
    });
}

#[test]
fn g0_elastic_source_domain_is_enforced_independently_of_expansion_support() {
    with_cx(|cx| {
        let (card, elastic, expansion, law) = fixture(350.0);
        let source = source(&card, &elastic, &expansion);
        each_problem(&law, 2, cx, |problem| {
            assert!(matches!(
                problem.resolve(&source, cx),
                Err(SolidMaterialError::Material(_))
            ));
        });
    });
}
