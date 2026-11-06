//! G1 directional tensile specimens; G3 frame/pointwise covariance; G4 refusals.
//! All data are synthetic analytic fixtures, not a measured failure calibration.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;
#[path = "orthotropic_tensile_strength/tsai_wu.rs"]
mod tsai_wu;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ORTHOTROPIC_POISSON_RATIO_PROPERTIES,
    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    resolve_orthotropic_elastic_state_point,
};
use fs_material::strength::OrthotropicTensileLimit;
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::{
    QuadraticThermalResponse, QuadraticThermalStrainField, TetAssemblyBudget,
    TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig, TetThermalStrainField,
    TetThermalStrainState, TetThermalStressReport,
};
use fs_solid::strength::TensileAssessmentError;
use fs_solid::strength::orthotropic::{
    OrthotropicTensileMaterialBinding as Binding, OrthotropicTensileMaterialField as Field,
};

const LIMITS: [&str; 3] = ["tensile_axis_1", "tensile_axis_2", "tensile_axis_3"];
const E: [f64; 3] = [2e6, 1e6, 0.5e6];
const FRAME: ContentHash = ContentHash([221; 32]);
const ORIENTATION: ContentHash = ContentHash([222; 32]);
const IDENTITY: [[f64; 3]; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TET: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const HEX: [[f64; 3]; 8] = [
    [0., 0., 0.],
    [1., 0., 0.],
    [1., 1., 0.],
    [0., 1., 0.],
    [0., 0., 1.],
    [1., 0., 1.],
    [1., 1., 1.],
    [0., 1., 1.],
];

fn with_cx<R>(cancelled: bool, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 221,
                kernel_id: 221,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn card(tag: u32) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut insert = |name: &str, dims: Dims, low: f64, high: f64| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300., low), (500., high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300., 500.),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic directional specimen".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    insert(DENSITY_PROPERTY, Density::DIMS, 1000., 1000.);
    for axis in 0..3 {
        insert(
            ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[axis],
            Pressure::DIMS,
            E[axis],
            E[axis],
        );
        insert(
            ORTHOTROPIC_POISSON_RATIO_PROPERTIES[axis],
            Dims::NONE,
            0.,
            0.,
        );
        insert(
            ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[axis],
            Pressure::DIMS,
            0.4e6,
            0.4e6,
        );
        let limit = [1000., 100., 50.][axis];
        insert(LIMITS[axis], Pressure::DIMS, limit, 0.5 * limit);
        insert(
            ["compression_1", "compression_2", "compression_3"][axis],
            Pressure::DIMS,
            1.5 * limit,
            0.75 * limit,
        );
        let shear = [100., 40., 25.][axis];
        insert(
            ["strength_12", "strength_23", "strength_31"][axis],
            Pressure::DIMS,
            shear,
            0.5 * shear,
        );
        let interaction = [0.2, -0.3, 0.1][axis];
        insert(
            ["f12", "f23", "f31"][axis],
            Dims::NONE,
            interaction,
            interaction,
        );
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic orthotropic specimen".into(),
            phase: "solid".into(),
            process: format!("analytic test {tag}"),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn binding(card: &MaterialCard, temperature: f64, rotation: [[f64; 3]; 3]) -> Binding {
    let point = QueryPoint::new().with("T", temperature).unwrap();
    let state = resolve_orthotropic_elastic_state_point(
        card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        0.01,
    )
    .unwrap();
    let limit = OrthotropicTensileLimit::resolve(
        card,
        &point,
        LIMITS,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    Binding::try_from_resolved_orthotropic_state(&state, limit, rotation, ORIENTATION).unwrap()
}

fn thermal(binding: &Binding, material_strain: [f64; 6]) -> TetThermalStrainState {
    let spatial = fs_material::tensor::rotate(&material_strain, &binding.principal_to_world());
    TetThermalStrainState::try_new(
        binding.limit().resolved().state_coordinate("T").unwrap(),
        300.,
        fs_material::tensor::to_mandel(&spatial),
        binding.material().material_state_identity,
        ContentHash([223; 32]),
        0.01,
    )
    .unwrap()
}

fn tet_report(binding: &Binding, strain: [f64; 6], cx: &Cx<'_>) -> TetThermalStressReport {
    let fixed: Vec<_> = (0..12).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &TET,
        tetrahedra: &[[0, 1, 2, 3]],
        materials: TetMaterialField::Uniform(binding.material()),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
    };
    let state = thermal(binding, strain);
    let field = TetThermalStrainField::Uniform(&state);
    let solved = problem
        .solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
        .unwrap();
    problem
        .recover_thermal_stress(&solved, field, FRAME, 0.01, cx)
        .unwrap()
}

#[test]
fn g1_g3_directional_tensile_demand_uses_material_axes_under_rotation() {
    with_cx(false, |cx| {
        let source = card(0);
        for rotation in [IDENTITY, fs_material::tensor::rotation([0., 0., 1.], 0.63)] {
            let binding = binding(&source, 300., rotation);
            // Both directions carry 200 Pa: weak axis 2 fails despite axis 1's
            // strength being ten times larger. An isotropic max eigenvalue
            // compared against the longitudinal strength would miss this.
            let response = tet_report(&binding, [-0.0001, -0.0002, 0., 0., 0., 0.], cx);
            let report = response
                .assess_orthotropic_tensile_strength(Field::Uniform(&binding), 1, cx)
                .unwrap();
            let assessment = report.governing_sample().unwrap().assessment();
            assert_eq!(assessment.governing_axis(), 1);
            for (actual, expected) in assessment
                .tensile_utilizations()
                .into_iter()
                .zip([0.2, 2., 0.])
            {
                assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
            }
            assert!(report.reaches_limit());
            assert_eq!(report.response_identity(), response.identity());
            assert_eq!(
                report,
                response
                    .assess_orthotropic_tensile_strength(
                        Field::PerElement(&[binding.clone()]),
                        1,
                        cx
                    )
                    .unwrap()
            );
            for strain in [[0.001; 6], [0., 0., 0., 0.001, 0., 0.]] {
                let compressed = tet_report(&binding, strain, cx);
                let assessment = compressed
                    .assess_orthotropic_tensile_strength(Field::Uniform(&binding), 1, cx)
                    .unwrap();
                assert!(
                    assessment
                        .governing_sample()
                        .unwrap()
                        .assessment()
                        .utilization()
                        < 1e-12
                );
            }
        }
    });
}

fn check_bending(response: &QuadraticThermalResponse, binding: &Binding, peak: f64, cx: &Cx<'_>) {
    let report = response
        .assess_orthotropic_tensile_strength(Field::Uniform(binding), 64, cx)
        .unwrap();
    let governing = report.governing_sample().unwrap();
    assert_eq!(governing.assessment().governing_axis(), 1);
    assert!((governing.assessment().utilization() - peak / 100.).abs() < 1e-10);
    assert!(governing.quadrature_point().is_some());
    assert!(report.reaches_limit());
    assert!(
        response.elements()[0]
            .mean_stress_mandel_pa()
            .iter()
            .all(|v| v.abs() < 1e-9)
    );
    // Independently integrate 1/2 E_2 [0.0004(z-1/2)]^2 on the unit cube.
    assert!((response.elastic_energy_j() - E[1] * 0.0004_f64.powi(2) / 24.).abs() < 1e-12);
}

#[test]
fn g1_q1_q2_directional_bending_peaks_survive_zero_mean() {
    with_cx(false, |cx| {
        let binding = binding(
            &card(0),
            300.,
            fs_material::tensor::rotation([0., 0., 1.], 0.63),
        );
        let gauss = [-1. / 3.0_f64.sqrt(), 1. / 3.0_f64.sqrt()];
        let states: Vec<_> = (0..8)
            .map(|p| thermal(&binding, [0., 0.0002 * gauss[p % 2], 0., 0., 0., 0.]))
            .collect();
        let fixed: Vec<_> = (0..24).collect();
        let p = Hex8ElasticProblem {
            nodes_m: &HEX,
            hexahedra: &[[0, 1, 2, 3, 4, 5, 6, 7]],
            materials: TetMaterialField::Uniform(binding.material()),
            fixed_dofs: &fixed,
            formulation: Formulation::Standard,
            budget: TetAssemblyBudget::standard(),
            stiffness_points_per_dir: 2,
        };
        let response = p
            .solve_thermal_response(
                QuadraticThermalStrainField::AtQuadrature {
                    states: &states,
                    points_per_element: 8,
                },
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            )
            .unwrap();
        check_bending(&response, &binding, E[1] * 0.0002 / 3.0_f64.sqrt(), cx);
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.; 3], family, cx);
            let fixed: Vec<_> = (0..mesh.nodes_m.len() * 3).collect();
            let p = hex_support::problem(&mesh, binding.material(), &fixed);
            let states: Vec<_> = p
                .quadrature
                .points()
                .unwrap()
                .iter()
                .map(|(r, _)| thermal(&binding, [0., 0.0002 * r[2], 0., 0., 0., 0.]))
                .collect();
            let response = p
                .solve_thermal_response(
                    QuadraticThermalStrainField::AtQuadrature {
                        states: &states,
                        points_per_element: 27,
                    },
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.01,
                    cx,
                )
                .unwrap();
            check_bending(
                &response,
                &binding,
                E[1] * 0.0002 * (3.0_f64 / 5.).sqrt(),
                cx,
            );
        }
    });
}

#[test]
fn g3_original_pointwise_orientation_and_temperature_control_strength() {
    with_cx(false, |cx| {
        let source = card(0);
        let cold = binding(&source, 300., IDENTITY);
        let hot = binding(
            &source,
            500.,
            fs_material::tensor::rotation([0., 0., 1.], 0.63),
        );
        let bindings: Vec<_> = (0..8)
            .map(|p| if p < 4 { cold.clone() } else { hot.clone() })
            .collect();
        let materials: Vec<_> = bindings.iter().map(|b| b.material().clone()).collect();
        let states: Vec<_> = bindings
            .iter()
            .map(|b| thermal(b, [0., -0.00008, 0., 0., 0., 0.]))
            .collect();
        let fixed: Vec<_> = (0..24).collect();
        let p = Hex8ElasticProblem {
            nodes_m: &HEX,
            hexahedra: &[[0, 1, 2, 3, 4, 5, 6, 7]],
            materials: TetMaterialField::Uniform(cold.material()),
            fixed_dofs: &fixed,
            formulation: Formulation::Standard,
            budget: TetAssemblyBudget::standard(),
            stiffness_points_per_dir: 2,
        };
        let response = p
            .solve_thermal_response(
                QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
                    states: &states,
                    materials: &materials,
                    points_per_element: 8,
                },
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            )
            .unwrap();
        let report = response
            .assess_orthotropic_tensile_strength(Field::PerPoint(&bindings), 8, cx)
            .unwrap();
        for (i, sample) in report.samples().iter().enumerate() {
            assert!(
                (sample.assessment().utilization() - if i < 4 { 0.8 } else { 1.6 }).abs() < 1e-12
            );
        }
        assert_eq!(
            report.governing_sample().unwrap().quadrature_point(),
            Some(4)
        );
        assert_eq!(
            response.assess_orthotropic_tensile_strength(Field::Uniform(&cold), 8, cx),
            Err(TensileAssessmentError::SourceMismatch {
                element: 0,
                point: Some(4)
            })
        );
    });
}

#[test]
fn g4_source_orientation_budget_count_and_cancellation_refuse_atomically() {
    let source = card(0);
    let original = binding(&source, 300., IDENTITY);
    let rotated = binding(
        &source,
        300.,
        fs_material::tensor::rotation([0., 0., 1.], 0.63),
    );
    let hot = binding(&source, 500., IDENTITY);
    let other = binding(&card(1), 300., IDENTITY);
    let state = resolve_orthotropic_elastic_state_point(
        &source,
        &QueryPoint::new().with("T", 300.).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.01,
    )
    .unwrap();
    for limit in [hot.limit(), other.limit()] {
        assert!(
            Binding::try_from_resolved_orthotropic_state(
                &state,
                limit.clone(),
                IDENTITY,
                ORIENTATION
            )
            .is_err()
        );
    }
    for rotation in [
        [[1., 1., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., -1.]],
    ] {
        assert!(
            Binding::try_from_resolved_orthotropic_state(
                &state,
                original.limit().clone(),
                rotation,
                ORIENTATION
            )
            .is_err()
        );
    }
    let response = with_cx(false, |cx| {
        tet_report(&original, [0., -0.0002, 0., 0., 0., 0.], cx)
    });
    let saved = response.clone();
    with_cx(false, |cx| {
        assert_eq!(
            response.assess_orthotropic_tensile_strength(Field::Uniform(&original), 0, cx),
            Err(TensileAssessmentError::Budget)
        );
        assert!(matches!(
            response.assess_orthotropic_tensile_strength(Field::PerElement(&[]), 1, cx),
            Err(TensileAssessmentError::Invalid(_))
        ));
        for mismatched in [&rotated, &hot, &other] {
            assert_eq!(
                response.assess_orthotropic_tensile_strength(Field::Uniform(mismatched), 1, cx),
                Err(TensileAssessmentError::SourceMismatch {
                    element: 0,
                    point: None
                })
            );
        }
    });
    with_cx(true, |cx| {
        assert_eq!(
            response.assess_orthotropic_tensile_strength(Field::Uniform(&original), 1, cx),
            Err(TensileAssessmentError::Cancelled)
        )
    });
    assert_eq!(response, saved);
}
