//! Actual sourced elasticity and S-N fatigue share one original material state.
//! Synthetic curves verify binding and arithmetic, not component service life.
#![deny(unsafe_code)]

#[path = "tet10/support.rs"]
#[allow(dead_code)]
mod support;

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::source::{
    SnLifeMeasure, SnSourceRequest, SnStressMeasure, resolve_log_log_sn,
};
use fs_material::fatigue::{MeanStressCorrection, RainflowLimits, SnCurve};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use fs_solid::fatigue::{
    TetFatigueConfig, TetFatigueError, TetFatigueGauge, TetFatigueHistoryConfig, TetFatiguePath,
};
use fs_solid::linear3::source::SolidStiffnessModel;
use fs_solid::linear3::{
    TetElasticMaterial, TetLinearElasticProblem, TetMaterialField, TetThermalStrainField,
    TetThermalStrainState,
};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const FRAME: ContentHash = ContentHash([119; 32]);
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Fixture {
    card: MaterialCard,
    fit_point: QueryPoint,
    bundle: Arc<ResolvedMaterialStatePoint>,
    material: TetElasticMaterial,
    thermal: TetThermalStrainState,
    curve: Arc<SnCurve>,
}

fn request() -> SnSourceRequest {
    SnSourceRequest {
        stress_measure: SnStressMeasure::Amplitude,
        life_measure: SnLifeMeasure::Cycles,
        temperature_axis: "absolute-temperature".into(),
        stress_ratio_axis: "calibration-R".into(),
        surface_and_geometry: "synthetic nominal-stress coupon".into(),
        environment: "synthetic fixed humidity".into(),
        load_basis: "fully reversed normal stress amplitude, full cycles".into(),
        maximum_points: 2,
    }
}

fn fixture(revision: bool) -> Fixture {
    let request = request();
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let dimensionless = QuantitySpec::dimensional(Dims::NONE);
    let physical_domain = ValidityDomain::unconstrained()
        .with_quantity(&request.temperature_axis, temperature, 290.0, 330.0)
        .with_quantity("humidity", dimensionless, 0.0, 1.0);
    let physical_point = QueryPoint::new()
        .with_quantity(&request.temperature_axis, temperature, 300.0)
        .unwrap()
        .with_quantity("humidity", dimensionless, 0.2)
        .unwrap();
    let life = request.life_measure.axis_name();
    let fit_domain = physical_domain
        .clone()
        .with_quantity(&request.stress_ratio_axis, dimensionless, -1.0, -1.0)
        .with_quantity(life, dimensionless, 100.0, 1000.0);
    let fit_point = physical_point
        .clone()
        .with_quantity(&request.stress_ratio_axis, dimensionless, -1.0)
        .unwrap()
        .with_quantity(life, dimensionless, 100.0)
        .unwrap();
    let mut claims = ClaimSet::new();
    let mut insert = |key, value, validity, interpolation| {
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic native S-N source fixture; not measured".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for (name, dims, value) in [
        ("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, 200e9),
        ("poisson_ratio", Dims::NONE, 0.0),
    ] {
        insert(
            PropertyKey::new(name, dims),
            PropertyValue::Scalar { value, dims },
            physical_domain.clone(),
            InterpolationPolicy::ConstantWithinValidity,
        );
    }
    insert(
        request.stress_measure.property_key(),
        PropertyValue::Curve {
            abscissa: life.into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(100.0, 10e6), (1000.0, 1e6)],
            dims: Pressure::DIMS,
        },
        fit_domain,
        InterpolationPolicy::TabulatedOnly,
    );
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic native fatigue alloy".into(),
            phase: "solid".into(),
            process: "synthetic coupon process".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    if revision {
        card = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    }
    let bundle = Arc::new(
        SolidStiffnessModel::Isotropic
            .resolve_reference(
                &card,
                &physical_point,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
    );
    let material = SolidStiffnessModel::Isotropic
        .material(
            &bundle,
            bundle.property("density").unwrap().value_si(),
            bundle.identity(),
            0.01,
        )
        .unwrap();
    let thermal = TetThermalStrainState::try_new(
        300.0,
        300.0,
        [0.0; 6],
        bundle.identity().0,
        ContentHash([120; 32]),
        0.01,
    )
    .unwrap();
    let curve = Arc::new(resolve_log_log_sn(&card, &fit_point, &request, SINGLE).unwrap());
    Fixture {
        card,
        fit_point,
        bundle,
        material,
        thermal,
        curve,
    }
}

fn gauge(source: &Fixture) -> TetFatigueGauge {
    TetFatigueGauge {
        element: 0,
        normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: source.bundle.identity().0,
        material_source: Some(source.bundle.clone()),
        curve: source.curve.clone(),
        mean_stress: MeanStressCorrection::Uncorrected,
    }
}

fn config() -> TetFatigueConfig {
    TetFatigueConfig {
        solver: Default::default(),
        maximum_absolute_mandel_strain: 0.01,
        history: RainflowLimits {
            maximum_samples: 40,
            maximum_residue: 16,
        },
        maximum_gauges: 1,
        maximum_total_residue: 16,
    }
}

fn history() -> TetFatigueHistoryConfig {
    let config = config();
    TetFatigueHistoryConfig {
        maximum_absolute_mandel_strain: config.maximum_absolute_mandel_strain,
        history: config.history,
        maximum_gauges: config.maximum_gauges,
        maximum_total_residue: config.maximum_total_residue,
    }
}

#[test]
fn original_sn_and_elastic_bundle_drive_actual_static_cycles_without_changing_mechanics() {
    let source = fixture(false);
    let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&source.material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    support::with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let thermal = TetThermalStrainField::Uniform(&source.thermal);
        let mut path =
            TetFatiguePath::new(&problem, thermal, FRAME, vec![gauge(&source)], config(), cx)
                .unwrap();
        assert!(std::ptr::eq(
            path.material_source(0).unwrap(),
            source.bundle.as_ref()
        ));
        assert_eq!(source.bundle.state_coordinate("T"), None);
        assert_eq!(
            source.bundle.state_coordinate("absolute-temperature"),
            Some(300.0)
        );
        assert_eq!(source.bundle.state_coordinate("calibration-R"), None);
        assert_eq!(source.bundle.state_coordinate("cycles-to-failure"), None);
        let retained = path.curve(0).unwrap().source().unwrap();
        assert_eq!(retained.card_identity(), source.card.content_hash());
        assert_eq!(retained.query(), &source.fit_point);
        for (index, stress) in [2e6, -2e6, 2e6].into_iter().enumerate() {
            // Unit tetrahedron: V=1/6 and only x displacement at node1 is free.
            let load = geometry
                .assemble(&[], VolumeLoad3::None, &[(3, stress / 6.0)], cx)
                .unwrap();
            let original = load
                .solve(&problem, thermal, Default::default(), FRAME, 0.01, cx)
                .unwrap();
            let step = path.advance(index as f64, &load, cx).unwrap();
            assert_eq!(step.mechanics.stress, original.stress);
            assert_eq!(step.mechanics.displacement, original.displacement);
            assert!((step.signed_stress_pa[0] - stress).abs() < 0.01);
        }
        let result = path.finish(cx).unwrap();
        assert_eq!(result[0].samples, 3);
        assert!((result[0].miner_damage - 0.002).abs() < 1e-10); // N=1e9/(2MPa)=500.
        assert_eq!(
            path.curve(0).unwrap().source().unwrap(),
            source.curve.source().unwrap()
        );
    });
}

#[test]
fn static_and_inertial_owners_refuse_foreign_sources_and_preserve_bound_checkpoint_evidence() {
    let source = fixture(false);
    let foreign = fixture(true);
    let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&source.material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    support::with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, problem.budget, 4, cx).unwrap();
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let thermal = TetThermalStrainField::Uniform(&source.thermal);
        let mut initial = ZERO;
        initial[1][0] = 2e6 / 200e9;
        let native = load
            .begin_transient(
                &problem,
                Some(thermal),
                0.0,
                &initial,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let bound = gauge(&source);
        let mut missing = bound.clone();
        missing.material_source = None;
        let mut wrong_bundle = bound.clone();
        wrong_bundle.material_source = Some(foreign.bundle);
        let mut wrong_card = bound.clone();
        wrong_card.curve = foreign.curve;
        let different_humidity = source
            .fit_point
            .clone()
            .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.8)
            .unwrap();
        let mut wrong_context = bound.clone();
        wrong_context.curve = Arc::new(
            resolve_log_log_sn(&source.card, &different_humidity, &request(), SINGLE).unwrap(),
        );
        for bad in [missing, wrong_bundle, wrong_card, wrong_context] {
            assert!(matches!(
                TetFatiguePath::new(&problem, thermal, FRAME, vec![bad.clone()], config(), cx),
                Err(TetFatigueError::SourceMismatch { .. })
            ));
            assert!(matches!(
                native.clone().with_fatigue(FRAME, vec![bad], history(), cx),
                Err(TetFatigueError::SourceMismatch { .. })
            ));
        }
        let mut original = native.clone();
        let mut observed = native
            .with_fatigue(FRAME, vec![bound], history(), cx)
            .unwrap();
        assert!(std::ptr::eq(
            observed.material_source(0).unwrap(),
            source.bundle.as_ref()
        ));
        let h = 2.0 * (std::f64::consts::PI / 32.0).tan() / (10.0_f64 * 200e9 / 7800.0).sqrt();
        for _ in 0..16 {
            let mechanical = original.advance(&load, h, cx).unwrap();
            assert_eq!(
                observed.advance(&load, h, cx).unwrap().mechanics,
                mechanical
            );
        }
        let mut checkpoint = observed.clone();
        for _ in 0..16 {
            let mechanical = original.advance(&load, h, cx).unwrap();
            assert_eq!(
                observed.advance(&load, h, cx).unwrap().mechanics,
                mechanical
            );
            checkpoint.advance(&load, h, cx).unwrap();
        }
        assert_eq!(observed.mechanics().state(), original.state());
        assert_eq!(checkpoint.mechanics().state(), observed.mechanics().state());
        assert_eq!(checkpoint.summaries(), observed.summaries());
        let result = observed.finish(cx).unwrap();
        assert_eq!(checkpoint.finish(cx).unwrap(), result);
        assert!((result[0].miner_damage - 0.002).abs() < 1e-10);
        assert!(std::ptr::eq(
            checkpoint.material_source(0).unwrap(),
            source.bundle.as_ref()
        ));
        assert_eq!(
            checkpoint.curve(0).unwrap().source().unwrap(),
            source.curve.source().unwrap()
        );
    });
}
