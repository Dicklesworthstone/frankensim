//! G1 manufactured tensile loading, G3 quadrature/temperature response and G4 refusal.
//! Synthetic material data test implementation, not a measured failure population.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY,
    YOUNG_MODULUS_PROPERTY, resolve_isotropic_elastic_state_point,
};
use fs_material::strength::IsotropicTensileLimit;
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::{
    QuadraticThermalResponse, QuadraticThermalStrainField, TetAssemblyBudget,
    TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig, TetThermalStrainField,
    TetThermalStrainState, TetThermalStressReport,
};
use fs_solid::strength::{
    TensileAssessmentError, TensileMaterialBinding, TensileMaterialField,
};

const LIMIT: &str = "synthetic_tensile_limit";
const FRAME: ContentHash = ContentHash([217; 32]);
const E: f64 = 2e6;
const HEX_NODES: [[f64; 3]; 8] = [
    [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
    [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.],
];
const HEX: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const TET_NODES: [[f64; 3]; 4] = [
    [0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.],
];

fn with_cx<R>(cancelled: bool, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| run(&Cx::new(
        &gate, arena,
        StreamKey { seed: 217, kernel_id: 217, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic,
    )))
}

fn card(tag: u32) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, low, high) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0, 1.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E, E),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0, 0.0),
        (LIMIT, Pressure::DIMS, 1000.0, 500.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS, knots: vec![(300.0, low), (500.0, high)], dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance { source: "synthetic tensile test".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic isotropic solid".into(),
        phase: "solid".into(), process: format!("analytic test {tag}"), revision: 0 }, claims, vec![]).unwrap()
}

fn binding(card: &MaterialCard, temperature: f64) -> TensileMaterialBinding {
    let query = QueryPoint::new().with("T", temperature).unwrap();
    let elastic = resolve_isotropic_elastic_state_point(
        card, &query, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let limit = IsotropicTensileLimit::resolve(
        card, &query, LIMIT, MaterialPropertySelection::SingleClaimOnly).unwrap();
    TensileMaterialBinding::from_isotropic_elastic(&elastic, limit).unwrap()
}

fn thermal(binding: &TensileMaterialBinding, strain: [f64; 6]) -> TetThermalStrainState {
    let t = binding.limit().resolved().state_coordinate("T").unwrap();
    TetThermalStrainState::try_new(t, 300.0, strain,
        binding.material().material_state_identity, ContentHash([218; 32]), 0.01).unwrap()
}

fn tet_report(binding: &TensileMaterialBinding, strain: [f64; 6], cx: &Cx<'_>)
    -> TetThermalStressReport {
    let fixed: Vec<_> = (0..12).collect();
    let p = TetLinearElasticProblem { nodes_m: &TET_NODES, tetrahedra: &[[0, 1, 2, 3]],
        materials: TetMaterialField::Uniform(binding.material()), fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard() };
    let state = thermal(binding, strain);
    let field = TetThermalStrainField::Uniform(&state);
    let solved = p.solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx).unwrap();
    p.recover_thermal_stress(&solved, field, FRAME, 0.01, cx).unwrap()
}

#[test]
fn g1_tet_tension_compression_hydrostatic_and_mandel_shear_use_actual_stress() {
    with_cx(false, |cx| {
        let material = binding(&card(0), 300.0);
        for (strain, expected) in [
            ([-0.001, 0., 0., 0., 0., 0.], 2000.0),
            ([0.001, 0.001, 0.001, 0., 0., 0.], 0.0),
            ([-0.001, -0.001, -0.001, 0., 0., 0.], 2000.0),
            ([0., 0., 0., -0.001 * std::f64::consts::SQRT_2, 0., 0.], 2000.0),
        ] {
            let stress = tet_report(&material, strain, cx);
            let report = stress.assess_tensile_strength(TensileMaterialField::Uniform(&material), 1, cx).unwrap();
            let worst = report.governing_sample().unwrap();
            assert_eq!(report.response_identity(), stress.identity());
            assert_eq!(worst.quadrature_point(), None);
            assert!((worst.assessment().tensile_stress_pa() - expected).abs() < 1e-8);
            assert!((worst.assessment().utilization() - expected / 1000.0).abs() < 1e-10);
            assert_eq!(report.reaches_limit(), expected >= 1000.0);
            assert_eq!(report, stress.assess_tensile_strength(
                TensileMaterialField::Uniform(&material), 1, cx).unwrap());
        }
    });
}

fn check_bending(response: &QuadraticThermalResponse, material: &TensileMaterialBinding,
    expected_peak: f64, cx: &Cx<'_>) {
    let element = &response.elements()[0];
    assert!(element.mean_stress_mandel_pa().iter().all(|v| v.abs() < 1e-8));
    let report = response.assess_tensile_strength(TensileMaterialField::Uniform(material), 64, cx).unwrap();
    let sample = report.governing_sample().unwrap();
    assert!(sample.quadrature_point().is_some());
    assert!((sample.assessment().tensile_stress_pa() - expected_peak).abs() < 1e-8);
    assert!(report.reaches_limit());
    // Independently integrate E/2 [0.002(z-1/2)]^2 on the unit cube.
    assert!((response.elastic_energy_j() - E * 0.002_f64.powi(2) / 24.0).abs() < 1e-10);
    for i in 0..6 {
        let integrated: f64 = element.stress_points().iter()
            .map(|p| p.integration_weight_m3() * p.stress_mandel_pa()[i]).sum();
        assert!((integrated / element.quadrature_volume_m3() - element.mean_stress_mandel_pa()[i]).abs() < 1e-9);
    }
    for axis in 0..3 {
        assert!(response.nodal_internal_force_n().iter().map(|f| f[axis]).sum::<f64>().abs() < 1e-8);
    }
}

#[test]
fn g1_q1_q2_bending_tensile_hotspots_survive_zero_element_mean() {
    with_cx(false, |cx| {
        let material = binding(&card(0), 300.0);
        // The clamped unit cube has sigma_xx=-E*0.002(z-1/2), with zero mean.
        let gauss = [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
        let states: Vec<_> = (0..8).map(|point| thermal(&material,
            [0.001 * gauss[point % 2], 0., 0., 0., 0., 0.])).collect();
        let fixed: Vec<_> = (0..24).collect();
        let p = Hex8ElasticProblem { nodes_m: &HEX_NODES, hexahedra: &HEX,
            materials: TetMaterialField::Uniform(material.material()), fixed_dofs: &fixed,
            formulation: Formulation::Standard, budget: TetAssemblyBudget::standard(),
            stiffness_points_per_dir: 2 };
        let response = p.solve_thermal_response(QuadraticThermalStrainField::AtQuadrature {
            states: &states, points_per_element: 8 }, TetStaticSolveConfig::default(), FRAME, 0.01, cx).unwrap();
        check_bending(&response, &material, E * 0.001 / 3.0_f64.sqrt(), cx);
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.; 3], family, cx);
            let fixed: Vec<_> = (0..mesh.nodes_m.len() * 3).collect();
            let p = hex_support::problem(&mesh, material.material(), &fixed);
            let rule = p.quadrature.points().unwrap();
            let states: Vec<_> = rule.iter().map(|(r, _)| thermal(&material,
                [0.001 * r[2], 0., 0., 0., 0., 0.])).collect();
            let response = p.solve_thermal_response(QuadraticThermalStrainField::AtQuadrature {
                states: &states, points_per_element: 27 }, TetStaticSolveConfig::default(), FRAME, 0.01, cx).unwrap();
            check_bending(&response, &material, E * 0.001 * (3.0_f64 / 5.0).sqrt(), cx);
        }
    });
}

#[test]
fn g3_local_temperature_limits_follow_original_quadrature_materials() {
    with_cx(false, |cx| {
        let card = card(0);
        let cold = binding(&card, 300.0);
        let hot = binding(&card, 500.0);
        let bindings: Vec<_> = (0..8).map(|p| if p < 4 { cold.clone() } else { hot.clone() }).collect();
        let materials: Vec<_> = bindings.iter().map(|b| b.material().clone()).collect();
        let states: Vec<_> = bindings.iter().map(|b| thermal(b, [-0.0004, 0., 0., 0., 0., 0.])).collect();
        let fixed: Vec<_> = (0..24).collect();
        let p = Hex8ElasticProblem { nodes_m: &HEX_NODES, hexahedra: &HEX,
            materials: TetMaterialField::Uniform(cold.material()), fixed_dofs: &fixed,
            formulation: Formulation::Standard, budget: TetAssemblyBudget::standard(),
            stiffness_points_per_dir: 2 };
        let response = p.solve_thermal_response(QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
            states: &states, materials: &materials, points_per_element: 8 },
            TetStaticSolveConfig::default(), FRAME, 0.01, cx).unwrap();
        let report = response.assess_tensile_strength(TensileMaterialField::PerPoint(&bindings), 8, cx).unwrap();
        for (i, sample) in report.samples().iter().enumerate() {
            assert!((sample.assessment().utilization() - if i < 4 { 0.8 } else { 1.6 }).abs() < 1e-10);
        }
        assert_eq!(report.governing_sample().unwrap().quadrature_point(), Some(4));
        assert_eq!(response.assess_tensile_strength(TensileMaterialField::Uniform(&cold), 8, cx),
            Err(TensileAssessmentError::SourceMismatch { element: 0, point: Some(4) }));
    });
}

#[test]
fn g4_source_context_counts_budgets_and_cancellation_refuse_without_mutation() {
    let source = card(0);
    let material = binding(&source, 300.0);
    let hot = binding(&source, 500.0);
    let other = binding(&card(1), 300.0);
    let elastic = resolve_isotropic_elastic_state_point(&source,
        &QueryPoint::new().with("T", 300.0).unwrap(), MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert!(TensileMaterialBinding::from_isotropic_elastic(&elastic, hot.limit().clone()).is_err());
    assert!(TensileMaterialBinding::from_isotropic_elastic(&elastic, other.limit().clone()).is_err());
    let stress = with_cx(false, |cx| tet_report(&material, [-0.001, 0., 0., 0., 0., 0.], cx));
    let saved = stress.clone();
    with_cx(false, |cx| {
        assert_eq!(stress.assess_tensile_strength(TensileMaterialField::Uniform(&material), 0, cx),
            Err(TensileAssessmentError::Budget));
        assert!(matches!(stress.assess_tensile_strength(TensileMaterialField::PerElement(&[]), 1, cx),
            Err(TensileAssessmentError::Invalid(_))));
        assert_eq!(stress.assess_tensile_strength(TensileMaterialField::Uniform(&hot), 1, cx),
            Err(TensileAssessmentError::SourceMismatch { element: 0, point: None }));
    });
    with_cx(true, |cx| assert_eq!(stress.assess_tensile_strength(
        TensileMaterialField::Uniform(&material), 1, cx), Err(TensileAssessmentError::Cancelled)));
    assert_eq!(stress, saved);
}
