//! Synthetic native rod cycles verify mechanics/observation coupling, not service life.
use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyLimits, StressClosureTolerance,
};
use fs_material::fatigue::source::{
    SnLifeMeasure, SnSourceRequest, SnStressMeasure, resolve_log_log_sn,
};
use fs_material::fatigue::{
    FatigueScope, MeanStressCorrection, NormalStressDirection, RainflowLimits, RainflowMiner,
    SnCurve, SnPoint,
};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
    LayerIsotropicParameters,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters, MaterialPropertySelection,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::rod::material_field::fatigue::{
    RodFatigueConfig, RodFatigueError, RodFatigueGauge, RodFatiguePoint, RodFindleyGauge,
    SpatialRodFatiguePath,
};
use fs_solid::rod::material_field::reference::RodReferenceGeometry;
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};
use fs_solid::rod::{TipLoad, WorldTipLoad};

const E: f64 = 1e6;
const R: f64 = 0.05;
const TOLERANCE: f64 = 1e-9;

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-8),
        "actual={actual}, expected={expected}"
    );
}
fn scope(temperature: f64) -> FatigueScope {
    FatigueScope {
        source_reference: "synthetic rod S-N fixture; not measured".into(),
        material: "synthetic rod alloy".into(),
        process: "synthetic fatigue coupon".into(),
        surface_and_geometry: "unnotched synthetic coupon".into(),
        environment: "synthetic fixed dry condition".into(),
        load_basis: "signed normal-stress amplitude, full cycles".into(),
        temperature_k: temperature,
    }
}
fn curve(temperature: f64, maximum_amplitude: f64) -> Arc<SnCurve> {
    Arc::new(
        SnCurve::new(
            // An S-N scope names the material "chemistry; phase", as an S-N source
            // writes it; a Findley scope names the chemistry.
            FatigueScope {
                material: "synthetic rod alloy; solid".into(),
                ..scope(temperature)
            },
            vec![
                SnPoint {
                    amplitude_pa: 1.0,
                    cycles_to_failure: 1e8,
                },
                SnPoint {
                    amplitude_pa: maximum_amplitude,
                    cycles_to_failure: (1e4 / maximum_amplitude).powi(2),
                },
            ],
            Some(0.1),
        )
        .unwrap(),
    )
}
fn source_request() -> SnSourceRequest {
    SnSourceRequest {
        stress_measure: SnStressMeasure::Amplitude,
        life_measure: SnLifeMeasure::Cycles,
        temperature_axis: "T".into(),
        stress_ratio_axis: "R".into(),
        surface_and_geometry: "unnotched synthetic coupon".into(),
        environment: "synthetic fixed dry condition".into(),
        load_basis: "signed normal-stress amplitude, full cycles".into(),
        maximum_points: 2,
    }
}
fn point(temperature: f64, quantity: QuantitySpec) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity, temperature)
        .unwrap()
}
fn card(young: f64, alpha: f64, revision: u32, quantity: QuantitySpec) -> MaterialCard {
    let domain = ValidityDomain::unconstrained().with_quantity("T", quantity, 290.0, 410.0);
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS),
        ("shear_modulus", 0.4 * young, Pressure::DIMS),
        (
            "linear_thermal_expansion_coefficient",
            alpha,
            INVERSE_TEMPERATURE_DIMS,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: domain.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let request = source_request();
    let life = request.life_measure.axis_name();
    claims
        .insert_claim(PropertyClaim {
            key: request.stress_measure.property_key(),
            value: PropertyValue::Curve {
                abscissa: life.into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(1.0, 1e4), (1e8, 1.0)],
                dims: Pressure::DIMS,
            },
            validity: domain
                .with_quantity("R", QuantitySpec::dimensional(Dims::NONE), -1.0, -1.0)
                .with_quantity(life, QuantitySpec::dimensional(Dims::NONE), 1.0, 1e8),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::TabulatedOnly,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    // A lineage starts at revision 0; each later revision supersedes the one before.
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic rod alloy".into(),
            phase: "solid".into(),
            process: "synthetic fatigue coupon".into(),
            revision: 0,
        },
        claims.clone(),
        vec![],
    )
    .unwrap();
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, claims.clone(), vec![]).unwrap();
    }
    card
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic native rod fixture".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn material(
    card: &MaterialCard,
    temperature: f64,
    quantity: QuantitySpec,
) -> ResolvedMaterialStatePoint {
    let requirements = [
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("shear_modulus", Pressure::DIMS),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    resolve_material_state_point(
        card,
        &point(temperature, quantity),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn section(materials: Vec<ResolvedMaterialStatePoint>) -> ConcentricShearSection {
    let count = materials.len();
    ConcentricElasticSection::try_new(
        0.0,
        materials
            .into_iter()
            .enumerate()
            .map(|(i, material)| ConcentricMaterialLayer {
                outer_radius_m: R * (i + 1) as f64 / count as f64,
                material,
            })
            .collect(),
    )
    .unwrap()
    .with_shear_modulus("shear_modulus", 0.8)
    .unwrap()
}
fn native(youngs: &[f64]) -> SpatialMaterialRod {
    let quantity = QuantitySpec::dimensional(Temperature::DIMS);
    let section = section(
        youngs
            .iter()
            .map(|&young| material(&card(young, 0.0, 0, quantity), 300.0, quantity))
            .collect(),
    );
    SpatialMaterialRod::try_straight(0.5, vec![section; 2]).unwrap()
}
fn isotropic_native(youngs: &[f64]) -> SpatialMaterialRod {
    let original = native(youngs);
    let layers = original.sections()[0].elastic().layers().to_vec();
    let bindings = layers
        .iter()
        .map(|layer| {
            let requirement = |name| layer.material.property(name).unwrap().requirement().clone();
            LayerIsotropicParameters {
                density: requirement("density"),
                parameters: IsotropicElasticParameters::YoungShear {
                    young_modulus: requirement("young_modulus"),
                    shear_modulus: requirement("shear_modulus"),
                },
            }
        })
        .collect::<Vec<_>>();
    let section =
        ConcentricElasticSection::try_new_with_isotropic_parameters(0.0, layers, &bindings)
            .unwrap()
            .with_isotropic_shear(0.8)
            .unwrap();
    SpatialMaterialRod::try_straight(0.5, vec![section; 2]).unwrap()
}
fn config() -> RodFatigueConfig {
    RodFatigueConfig {
        history: RainflowLimits {
            maximum_samples: 256,
            maximum_residue: 32,
        },
        maximum_gauges: 8,
        maximum_total_residue: 256,
        maximum_total_findley_samples: 2048,
        maximum_absolute_axial_strain: 0.02,
        maximum_engineering_shear: 0.02,
        maximum_unassessed_transverse_force_n: 1e-7,
    }
}
fn gauge(layer: usize, coordinates: [f64; 2], direction: [f64; 3]) -> RodFatigueGauge {
    RodFatigueGauge {
        point: RodFatiguePoint {
            segment: 0,
            layer,
            coordinates_yz_m: coordinates,
        },
        normal_in_material_frame: NormalStressDirection::new(direction).unwrap(),
        temperature_axis: "T".into(),
        curve: curve(300.0, 1e4),
        mean_stress: MeanStressCorrection::Uncorrected,
    }
}
fn findley(point: RodFatiguePoint) -> RodFindleyGauge {
    RodFindleyGauge {
        point,
        temperature_axis: "T".into(),
        calibration: Arc::new(FindleyCalibration::new(scope(300.0), 0.2, 2000.0).unwrap()),
        normals_in_material_frame: vec![[1.0, 0.0, 0.0]],
        limits: FindleyLimits {
            maximum_samples: 256,
            maximum_planes: 4,
            maximum_projection_work: 1024,
            maximum_circle_work: 100_000,
        },
        closure: StressClosureTolerance {
            absolute_stress_pa: 1e-4,
            relative: 1e-8,
        },
    }
}
fn axial(force: f64) -> TipLoad {
    TipLoad {
        force: [force, 0.0, 0.0],
        moment: [0.0; 3],
    }
}
fn torque(torque: f64) -> WorldTipLoad {
    WorldTipLoad {
        force: [0.0; 3],
        moment: [torque, 0.0, 0.0],
    }
}

#[test]
fn layered_axial_cycles_preview_and_partitioned_continuation_match_hand_count() {
    let gate = CancelGate::new();
    let rod = native(&[E, 2.0 * E]);
    let force = rod.sections()[0].elastic().axial_stiffness_n() * 0.001;
    let mut path = SpatialRodFatiguePath::new(
        rod,
        vec![
            gauge(0, [0.02, 0.0], [1.0, 0.0, 0.0]),
            gauge(1, [0.04, 0.0], [1.0, 0.0, 0.0]),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    assert_eq!(path.sample_count(), 1);
    path.advance(&[axial(force)], TOLERANCE, &gate).unwrap();
    let committed = path.summaries();
    assert!(
        committed
            .iter()
            .all(|s| s.miner_damage == 0.0 && !s.finished)
    );
    let preview = path.preview_fatigue(&gate).unwrap();
    near(
        preview[0].miner_damage,
        0.5 * (500.0_f64 / 1e4).powi(2),
        2e-7,
    );
    near(
        preview[1].miner_damage,
        0.5 * (1000.0_f64 / 1e4).powi(2),
        2e-7,
    );
    assert_eq!(path.summaries(), committed);
    near(
        path.observation().normal_samples()[0].stress_pa[0],
        1000.0,
        1e-8,
    );
    near(
        path.observation().normal_samples()[1].stress_pa[0],
        2000.0,
        1e-8,
    );
    let mut split = path.clone();
    path.advance(&[axial(-force), axial(force), axial(0.0)], TOLERANCE, &gate)
        .unwrap();
    for target in [axial(-force), axial(force), axial(0.0)] {
        split.advance(&[target], TOLERANCE, &gate).unwrap();
    }
    assert_eq!(path.rod().positions_m(), split.rod().positions_m());
    assert_eq!(path.summaries(), split.summaries());
    let expected = path.preview_fatigue(&gate).unwrap();
    let done = path.finish(&gate).unwrap();
    assert_eq!(done, expected);
    assert_eq!(done, split.finish(&gate).unwrap());
    near(done[0].miner_damage, 0.0125, 2e-7);
    near(done[1].miner_damage, 0.05, 2e-7);
    assert_eq!(path.sample_count(), 5);
    assert!(matches!(
        path.advance(&[axial(0.0)], TOLERANCE, &gate),
        Err(RodFatigueError::Finished)
    ));
}

#[test]
fn torsion_drives_signed_oblique_normal_life_and_original_findley_shear() {
    let gate = CancelGate::new();
    let rod = native(&[E]);
    let tau = 1000.0;
    let t = tau / (0.4 * E * R) * rod.sections()[0].torsional_stiffness_n_m2();
    let oblique = gauge(0, [R, 0.0], [1.0, 0.0, 1.0]);
    let mut path = SpatialRodFatiguePath::new(
        rod,
        vec![gauge(0, [R, 0.0], [1.0, 0.0, 0.0]), oblique.clone()],
        vec![findley(oblique.point)],
        config(),
        &gate,
    )
    .unwrap();
    path.advance_world(&[torque(t)], TOLERANCE, &gate).unwrap();
    let sample = path.observation().normal_samples()[1];
    near(sample.stress_pa[5], tau, 2e-7);
    near(
        oblique
            .normal_in_material_frame
            .project(&sample.stress_pa)
            .unwrap(),
        tau,
        2e-7,
    );
    path.advance_world(&[torque(-t)], TOLERANCE, &gate).unwrap();
    near(
        path.observation().normal_samples()[1].stress_pa[5],
        -tau,
        2e-7,
    );
    path.advance_world_and_close_findley_cycle(&[torque(0.0)], TOLERANCE, &gate)
        .unwrap();
    let result = path.observation().findley_results()[0].unwrap();
    near(result.shear_amplitude_pa, tau, 2e-7);
    near(result.usage_factor, 0.5, 2e-7);
    assert_eq!(result.samples, 4);
    let before = path.summaries();
    path.restart_findley_cycle(&gate).unwrap();
    assert_eq!(path.sample_count(), 4);
    assert_eq!(path.summaries(), before);
    path.advance_world_and_close_findley_cycle(
        &[torque(t), torque(-t), torque(0.0)],
        TOLERANCE,
        &gate,
    )
    .unwrap();
    let done = path.finish(&gate).unwrap();
    assert_eq!(done[0].miner_damage, 0.0);
    assert!(done[1].miner_damage > 0.0);
    assert_eq!(path.observation().findley_results()[0].unwrap().samples, 4);
}

#[test]
fn transverse_force_drives_native_fixed_point_findley_and_normal_life() {
    let gate = CancelGate::new();
    let rod = isotropic_native(&[E]);
    let oblique = gauge(0, [0.0; 2], [1.0, 1.0, 0.0]);
    let mut limits = config();
    limits.maximum_unassessed_transverse_force_n = 0.0;
    let load = |force| WorldTipLoad {
        force: [0.0, force, 0.0],
        moment: [0.0; 3],
    };
    let mut omitted =
        SpatialRodFatiguePath::new(rod.clone(), vec![oblique.clone()], vec![], limits, &gate)
            .unwrap();
    assert!(
        omitted
            .advance_world(&[load(0.5)], TOLERANCE, &gate)
            .is_err()
    );
    assert_eq!(omitted.sample_count(), 1);
    let mut path = SpatialRodFatiguePath::new_with_transverse_shear(
        rod,
        vec![oblique.clone()],
        vec![findley(oblique.point)],
        limits,
        &gate,
    )
    .unwrap();
    let mut stresses = Vec::new();
    for force in [0.5, -0.5] {
        path.advance_world(&[load(force)], TOLERANCE, &gate)
            .unwrap();
        let sample = path.observation().normal_samples()[0];
        let (native_force, _) = path.rod().segment_resultants(0).unwrap();
        let area = core::f64::consts::PI * R * R;
        // E/G=2.5 gives nu=1/4 and a centerline multiplier 7/5. No Q/A
        // uniform stress or the mechanical kappa=0.8 enters this reference.
        near(
            sample.stress_pa[3],
            (7.0 / 5.0) * native_force[1] / area,
            2e-12,
        );
        assert_eq!(sample.unassessed_transverse_force_n, [0.0; 2]);
        assert_eq!(
            path.observation().findley_samples()[0].stress_pa,
            sample.stress_pa
        );
        stresses.push(sample.stress_pa);
    }
    path.advance_world_and_close_findley_cycle(&[load(0.0)], TOLERANCE, &gate)
        .unwrap();
    let result = path.observation().findley_results()[0].unwrap();
    let amplitude = 0.5 * (stresses[0][3] - stresses[1][3]).abs();
    near(result.shear_amplitude_pa, amplitude, 2e-7);
    near(
        result.maximum_normal_stress_pa,
        stresses[0][0].max(stresses[1][0]).max(0.0),
        2e-6,
    );
    assert!(result.usage_factor > 0.0);
    assert!(path.finish(&gate).unwrap()[0].miner_damage > 0.0);
}

#[test]
fn transverse_recovery_refuses_foreign_family_and_combined_overstrain_atomically() {
    let gate = CancelGate::new();
    let observation = gauge(0, [0.0; 2], [1.0, 1.0, 0.0]);
    for rod in [native(&[E]), isotropic_native(&[E, 2.0 * E])] {
        assert!(
            SpatialRodFatiguePath::new_with_transverse_shear(
                rod,
                vec![observation.clone()],
                vec![],
                config(),
                &gate,
            )
            .is_err()
        );
    }
    let mut limits = config();
    limits.maximum_engineering_shear = 1e-4;
    let rod = isotropic_native(&[E]);
    let torsion = rod.sections()[0].torsional_stiffness_n_m2() * 6e-5 / R;
    let mut path = SpatialRodFatiguePath::new_with_transverse_shear(
        rod,
        vec![observation],
        vec![],
        limits,
        &gate,
    )
    .unwrap();
    let before = path.clone();
    assert!(
        path.advance_world(
            &[WorldTipLoad {
                force: [0.0, 0.15, 0.0],
                moment: [torsion, 0.0, 0.0],
            }],
            TOLERANCE,
            &gate
        )
        .is_err()
    );
    assert_eq!(path.rod().positions_m(), before.rod().positions_m());
    assert_eq!(path.rod().frames(), before.rod().frames());
    assert_eq!(path.summaries(), before.summaries());
    assert_eq!(
        path.observation().normal_samples(),
        before.observation().normal_samples()
    );
    assert_eq!(path.sample_count(), 1);
}

#[test]
fn bonded_layer_shear_drives_the_original_native_fatigue_history() {
    let gate = CancelGate::new();
    let rod = isotropic_native(&[E, 2.0 * E]);
    let center = gauge(0, [0.0; 2], [1.0, 1.0, 0.0]);
    let mut limits = config();
    limits.maximum_unassessed_transverse_force_n = 0.0;
    let mut path = SpatialRodFatiguePath::new_with_layered_transverse_shear(
        rod.clone(),
        vec![center.clone()],
        vec![findley(center.point)],
        limits,
        &gate,
    )
    .unwrap();
    let t = 0.5_f64;
    let b = -(3.5 / 4.0) * (1.0 - t * t) * t * t / (3.0 - t * t);
    let a_core = b - 3.5 / 4.0 + b / (t * t);
    let factor = -2.0 * a_core / (1.25 * (2.0 - t.powi(4)));
    let load = |force| WorldTipLoad {
        force: [0.0, force, 0.0],
        moment: [0.0; 3],
    };
    for force in [0.5, -0.5] {
        path.advance_world(&[load(force)], TOLERANCE, &gate)
            .unwrap();
        let (resultant, _) = path.rod().segment_resultants(0).unwrap();
        let sample = path.observation().normal_samples()[0];
        near(
            sample.stress_pa[3],
            factor * resultant[1] / (core::f64::consts::PI * R * R),
            2e-11,
        );
        near(
            sample.stress_pa[5],
            factor * resultant[2] / (core::f64::consts::PI * R * R),
            2e-11,
        );
        assert_eq!(sample.unassessed_transverse_force_n, [0.0; 2]);
        assert_eq!(
            sample.stress_pa,
            path.observation().findley_samples()[0].stress_pa
        );
    }
    path.advance_world_and_close_findley_cycle(&[load(0.0)], TOLERANCE, &gate)
        .unwrap();
    assert!(
        path.observation().findley_results()[0]
            .unwrap()
            .usage_factor
            > 0.0
    );
    assert!(path.finish(&gate).unwrap()[0].miner_damage > 0.0);
    assert_eq!(path.rod().sections(), rod.sections());
}

#[test]
fn pure_bending_preserves_opposite_material_point_histories() {
    let gate = CancelGate::new();
    let rod = native(&[E]);
    let moment = rod.sections()[0].elastic().bending_stiffness_n_m2() * 0.02;
    let target = |m| WorldTipLoad {
        force: [0.0; 3],
        moment: [0.0, 0.0, m],
    };
    let mut path = SpatialRodFatiguePath::new(
        rod,
        vec![
            gauge(0, [R, 0.0], [1.0, 0.0, 0.0]),
            gauge(0, [-R, 0.0], [1.0, 0.0, 0.0]),
        ],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    path.advance_world(&[target(moment)], TOLERANCE, &gate)
        .unwrap();
    near(
        path.observation().normal_samples()[0].stress_pa[0],
        -1000.0,
        2e-7,
    );
    near(
        path.observation().normal_samples()[1].stress_pa[0],
        1000.0,
        2e-7,
    );
    path.advance_world(
        &[target(-moment), target(moment), target(0.0)],
        TOLERANCE,
        &gate,
    )
    .unwrap();
    let done = path.finish(&gate).unwrap();
    near(done[0].miner_damage, 0.0125, 2e-6);
    near(done[1].miner_damage, 0.0125, 2e-6);
}

#[test]
fn fixed_thermal_mismatch_supplies_actual_layer_mean_stress() {
    let gate = CancelGate::new();
    let quantity = QuantitySpec::dimensional(Temperature::DIMS);
    let cards = [card(E, 1e-5, 0, quantity), card(E, 3e-5, 0, quantity)];
    let section = section(cards.iter().map(|c| material(c, 400.0, quantity)).collect());
    let paths = cards
        .iter()
        .map(|c| {
            integrate_isotropic_thermal_expansion(
                c,
                &point(300.0, quantity),
                &point(400.0, quantity),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
        })
        .collect();
    let thermal =
        ConcentricThermoelasticSection::try_new(section.elastic().clone(), paths).unwrap();
    let force = section.elastic().axial_stiffness_n() * 0.0005;
    let mut rod = SpatialMaterialRod::try_straight(0.5, vec![section])
        .unwrap()
        .with_thermal_sections(vec![thermal], 0.01, &gate)
        .unwrap();
    rod.solve_static(&TipLoad::default(), 1, TOLERANCE, &gate)
        .unwrap();
    let mut gauges = vec![
        gauge(0, [0.02, 0.0], [1.0, 0.0, 0.0]),
        gauge(1, [0.04, 0.0], [1.0, 0.0, 0.0]),
    ];
    let correction = MeanStressCorrection::GoodmanTensionOnly {
        ultimate_strength_pa: 1e4,
    };
    for gauge in &mut gauges {
        gauge.curve = curve(400.0, 1e4);
        gauge.mean_stress = correction;
    }
    let mut path = SpatialRodFatiguePath::new(rod, gauges, vec![], config(), &gate).unwrap();
    for (sample, expected) in path
        .observation()
        .normal_samples()
        .iter()
        .zip([1500.0, -500.0])
    {
        near(sample.stress_pa[0], expected, 2e-8);
    }
    assert_eq!(path.sample_count(), 1);
    let targets = [axial(force), axial(-force), axial(force), axial(0.0)];
    path.advance(&targets, TOLERANCE, &gate).unwrap();
    let done = path.finish(&gate).unwrap();
    for (index, mean) in [1500.0, -500.0].into_iter().enumerate() {
        let mut reference =
            RainflowMiner::new(curve(400.0, 1e4), correction, config().history).unwrap();
        reference
            .push_block(
                &[mean, mean + 500.0, mean - 500.0, mean + 500.0, mean],
                &mut || false,
            )
            .unwrap();
        near(
            done[index].miner_damage,
            reference.finish().unwrap().miner_damage,
            2e-7,
        );
    }
    assert!(done[0].miner_damage > done[1].miner_damage);
}

#[test]
fn late_gauge_and_findley_closure_refusals_restore_native_pose_and_all_histories() {
    let gate = CancelGate::new();
    let rod = native(&[E]);
    let force = rod.sections()[0].elastic().axial_stiffness_n() * 0.001;
    let first = gauge(0, [R, 0.0], [1.0, 0.0, 0.0]);
    let mut second = first.clone();
    second.curve = curve(300.0, 600.0);
    let mut path =
        SpatialRodFatiguePath::new(rod.clone(), vec![first, second], vec![], config(), &gate)
            .unwrap();
    let before = path.clone();
    assert!(matches!(
        path.advance(
            &[axial(force), axial(-force), axial(force), axial(-force)],
            TOLERANCE,
            &gate
        ),
        Err(RodFatigueError::Gauge { gauge: 1, .. })
    ));
    assert_eq!(path.rod().positions_m(), before.rod().positions_m());
    assert_eq!(path.rod().frames(), before.rod().frames());
    assert_eq!(path.summaries(), before.summaries());
    assert_eq!(
        path.observation().normal_samples(),
        before.observation().normal_samples()
    );
    assert_eq!(path.sample_count(), 1);

    let g = gauge(0, [R, 0.0], [1.0, 0.0, 1.0]);
    let t = rod.sections()[0].torsional_stiffness_n_m2() * 0.05;
    let mut path = SpatialRodFatiguePath::new(
        rod,
        vec![g.clone()],
        vec![findley(g.point)],
        config(),
        &gate,
    )
    .unwrap();
    let mut uninterrupted = path.clone();
    assert!(matches!(
        path.advance_world_and_close_findley_cycle(
            &[torque(t), torque(-t), torque(0.5 * t)],
            TOLERANCE,
            &gate
        ),
        Err(RodFatigueError::Findley { .. })
    ));
    assert_eq!(path.sample_count(), 1);
    assert_eq!(path.rod().frames(), uninterrupted.rod().frames());
    assert_eq!(path.summaries(), uninterrupted.summaries());
    for run in [&mut path, &mut uninterrupted] {
        run.advance_world_and_close_findley_cycle(
            &[torque(t), torque(-t), torque(0.0)],
            TOLERANCE,
            &gate,
        )
        .unwrap();
    }
    assert_eq!(
        path.finish(&gate).unwrap(),
        uninterrupted.finish(&gate).unwrap()
    );
    assert_eq!(
        path.observation().findley_results(),
        uninterrupted.observation().findley_results()
    );
}

#[test]
fn source_identity_temperature_quantity_and_annulus_cannot_be_substituted() {
    let gate = CancelGate::new();
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let original = card(E, 0.0, 0, temperature);
    let request = source_request();
    // A tabulated-only S-N source is selected at one of its own knots (1 and 1e8 cycles).
    let query = point(300.0, temperature)
        .with_quantity("R", QuantitySpec::dimensional(Dims::NONE), -1.0)
        .unwrap()
        .with_quantity(
            request.life_measure.axis_name(),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    let sourced = |card: &MaterialCard| {
        Arc::new(
            resolve_log_log_sn(
                card,
                &query,
                &request,
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .unwrap(),
        )
    };
    let rod = SpatialMaterialRod::try_straight(
        0.5,
        vec![section(vec![material(&original, 300.0, temperature)])],
    )
    .unwrap();
    let mut correct = gauge(0, [R, 0.0], [1.0, 0.0, 0.0]);
    correct.curve = sourced(&original);
    let accepted =
        SpatialRodFatiguePath::new(rod.clone(), vec![correct.clone()], vec![], config(), &gate)
            .unwrap();
    assert_eq!(
        accepted
            .observation()
            .material_source(0)
            .unwrap()
            .card_identity(),
        original.content_hash()
    );
    let mut wrong = correct.clone();
    wrong.curve = sourced(&card(E, 0.0, 1, temperature));
    assert!(matches!(
        SpatialRodFatiguePath::new(rod.clone(), vec![wrong], vec![], config(), &gate),
        Err(RodFatigueError::SourceMismatch { .. })
    ));
    wrong = correct.clone();
    wrong.temperature_axis = "absolute-temperature".into();
    assert!(matches!(
        SpatialRodFatiguePath::new(rod.clone(), vec![wrong], vec![], config(), &gate),
        Err(RodFatigueError::SourceMismatch { .. })
    ));
    wrong = correct;
    wrong.point.coordinates_yz_m = [1.1 * R, 0.0];
    assert!(SpatialRodFatiguePath::new(rod, vec![wrong], vec![], config(), &gate).is_err());

    // "T" typed as an untyped pressure is a sampled stress, which a rod layer refuses
    // as a frozen elastic coordinate before any gauge exists.
    let untyped = QuantitySpec::dimensional(Pressure::DIMS);
    let untyped_card = card(E, 0.0, 0, untyped);
    assert!(matches!(
        SpatialMaterialRod::try_straight(
            0.5,
            vec![section(vec![material(&untyped_card, 300.0, untyped)])],
        ),
        Err(SpatialRodError::Solid(
            fs_solid::SolidError::MaterialRefused { .. }
        ))
    ));
    // A static typed pressure is a fixed external condition the rod admits; the
    // gauge refuses it as a temperature coordinate.
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let wrong_card = card(E, 0.0, 0, pressure);
    let wrong_rod = SpatialMaterialRod::try_straight(
        0.5,
        vec![section(vec![material(&wrong_card, 300.0, pressure)])],
    )
    .unwrap();
    assert!(matches!(
        SpatialRodFatiguePath::new(
            wrong_rod,
            vec![gauge(0, [R, 0.0], [1.0, 0.0, 0.0])],
            vec![],
            config(),
            &gate
        ),
        Err(RodFatigueError::SourceMismatch { .. })
    ));
}

#[test]
fn native_transverse_force_and_outer_fiber_domains_require_explicit_admission() {
    let gate = CancelGate::new();
    let rod = native(&[E]);
    let mut path = SpatialRodFatiguePath::new(
        rod.clone(),
        vec![gauge(0, [R, 0.0], [1.0, 0.0, 1.0])],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    assert!(matches!(
        path.advance(
            &[TipLoad {
                force: [0.0, 0.1, 0.0],
                moment: [0.0; 3]
            }],
            TOLERANCE,
            &gate
        ),
        Err(RodFatigueError::Invalid(_))
    ));
    assert_eq!(path.sample_count(), 1);
    assert_eq!(path.rod().positions_m(), rod.positions_m());
    let mut declared = config();
    declared.maximum_unassessed_transverse_force_n = 1.0;
    let mut admitted = SpatialRodFatiguePath::new(
        rod.clone(),
        vec![gauge(0, [R, 0.0], [1.0, 0.0, 0.0])],
        vec![],
        declared,
        &gate,
    )
    .unwrap();
    admitted
        .advance(
            &[TipLoad {
                force: [0.0, 0.1, 0.0],
                moment: [0.0; 3],
            }],
            TOLERANCE,
            &gate,
        )
        .unwrap();
    assert!(
        admitted.observation().normal_samples()[0].unassessed_transverse_force_n[0].abs() > 0.09
    );

    let mut bounded = config();
    bounded.maximum_engineering_shear = 1e-5;
    let t = rod.sections()[0].torsional_stiffness_n_m2() * 0.05;
    let mut axis = SpatialRodFatiguePath::new(
        rod,
        vec![gauge(0, [0.0; 2], [1.0, 0.0, 0.0])],
        vec![],
        bounded,
        &gate,
    )
    .unwrap();
    assert!(axis.advance_world(&[torque(t)], TOLERANCE, &gate).is_err());
    assert_eq!(axis.sample_count(), 1); // A neutral-axis gauge cannot hide surface overstrain.
}

#[test]
fn intrinsic_twist_cancellation_and_global_sample_budget_preserve_accepted_state() {
    let gate = CancelGate::new();
    let section = native(&[E]).sections()[0].clone();
    let frames = (0..=2)
        .map(|i| So3::exp(So3Tangent::new(Vec3::new(0.2 * i as f64, 0.0, 0.0))).unwrap())
        .collect();
    let reference = RodReferenceGeometry::try_new(
        0.5,
        vec![[0.0; 3], [0.25, 0.0, 0.0], [0.5, 0.0, 0.0]],
        frames,
    )
    .unwrap();
    let rod = SpatialMaterialRod::try_from_reference(reference, vec![section; 2]).unwrap();
    let mut small = config();
    small.history.maximum_samples = 2;
    let mut path = SpatialRodFatiguePath::new(
        rod.clone(),
        vec![gauge(0, [R, 0.0], [1.0, 0.0, 1.0])],
        vec![],
        small,
        &gate,
    )
    .unwrap();
    assert_eq!(path.observation().normal_samples()[0].stress_pa, [0.0; 6]);
    assert_eq!(path.preview_fatigue(&gate).unwrap()[0].miner_damage, 0.0);
    assert!(matches!(
        path.advance_world(&[torque(0.0), torque(0.0)], TOLERANCE, &gate),
        Err(RodFatigueError::Budget)
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        path.advance_world(&[torque(0.0)], TOLERANCE, &cancelled),
        Err(RodFatigueError::Cancelled)
    ));
    assert!(matches!(
        path.preview_fatigue(&cancelled),
        Err(RodFatigueError::Cancelled)
    ));
    assert_eq!(path.rod().frames(), rod.frames());
    assert_eq!(path.sample_count(), 1);
    path.advance_world(&[torque(0.0)], TOLERANCE, &gate)
        .unwrap();
    assert_eq!(path.finish(&gate).unwrap()[0].miner_damage, 0.0);
}
