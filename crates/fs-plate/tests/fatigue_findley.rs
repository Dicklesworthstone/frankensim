//! Native DKT tensor recovery and multiaxial Findley block semantics.
//! Synthetic fits check implementation; they are not measured material data.
#[path = "fatigue/fixture.rs"]
#[allow(dead_code)]
mod fixture;
use fixture::*;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyError, FindleyLimits, StressClosureTolerance,
    source::{FindleyProperty, FindleySourceRequest, resolve_findley},
};
use fs_material::state_point::{
    IsotropicElasticStatePoint, MaterialPropertySelection, resolve_isotropic_elastic_state_point,
};
use fs_plate::fatigue::{
    FatiguePlate, PlateFatigueError, PlateFatigueSection, PlateFindleyBudget, PlateFindleyGauge,
};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use std::sync::Arc;

const K: f64 = 0.25;
const LIMIT: f64 = 5e6;
fn limits(samples: u64) -> FindleyLimits {
    FindleyLimits {
        maximum_samples: samples,
        maximum_planes: 8,
        maximum_projection_work: samples * 8,
        maximum_circle_work: 100_000,
    }
}
fn budget(gauges: &[PlateFindleyGauge]) -> PlateFindleyBudget {
    PlateFindleyBudget {
        maximum_total_samples: gauges.iter().map(|g| g.limits.maximum_samples).sum(),
        maximum_total_planes: gauges.iter().map(|g| g.normals_in_mesh_frame.len()).sum(),
        maximum_total_projection_work: gauges
            .iter()
            .map(|g| g.limits.maximum_projection_work)
            .sum(),
        maximum_total_circle_work: gauges.iter().map(|g| g.limits.maximum_circle_work).sum(),
    }
}
fn options(tension: f64) -> AssemblyOptions {
    AssemblyOptions {
        pretension: tension,
        support: EdgeSupport::Clamped,
    }
}
fn field(mesh: &PlateMesh, k: [f64; 3]) -> Vec<f64> {
    mesh.nodes
        .iter()
        .flat_map(|&(x, y)| {
            [
                -0.5 * (k[0] * x * x + k[1] * y * y + k[2] * x * y),
                -k[0] * x - 0.5 * k[2] * y,
                -k[1] * y - 0.5 * k[2] * x,
            ]
        })
        .collect()
}
fn findley_gauge(
    f: &Fixture,
    element: usize,
    depth_m: f64,
    normals: Vec<[f64; 3]>,
    samples: u64,
) -> PlateFindleyGauge {
    let material = f.material.resolved();
    let calibration = FindleyCalibration::new(
        FatigueScope {
            source_reference: "synthetic plate Findley fit".into(),
            material: material.material().chemistry.clone(),
            process: material.material().process.clone(),
            surface_and_geometry: "synthetic smooth homogeneous plate".into(),
            environment: "synthetic fixed atmosphere".into(),
            load_basis: "representative closed elastic block".into(),
            temperature_k: 300.0,
        },
        K,
        LIMIT,
    )
    .unwrap();
    PlateFindleyGauge {
        element,
        surface_point: 0,
        depth_m,
        temperature_axis: "T".into(),
        calibration: Arc::new(calibration),
        normals_in_mesh_frame: normals,
        limits: limits(samples),
        closure: StressClosureTolerance {
            absolute_stress_pa: 1e-7,
            relative: 1e-12,
        },
    }
}

#[test]
fn native_twist_pretension_and_two_windows_keep_usage_separate_from_miner() {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.01, 1, 1);
    let mean = 0.5e6;
    let amplitude = 2e6;
    let z = H / 2.0;
    let shear_modulus = E / (2.0 * (1.0 + NU));
    let gauges = vec![
        findley_gauge(&f, 0, z, vec![[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], 8),
        findley_gauge(&f, 0, -z, vec![[1.0, 0.0, 0.0]], 8),
    ];
    let plate = FatiguePlate::assemble_with_findley(
        &mesh,
        vec![section(&f); 2],
        &[],
        &options(mean * H),
        vec![gauge(&f, 0, z, [1.0, 1.0])],
        gauges.clone(),
        config(100),
        budget(&gauges),
        &mut || false,
    )
    .unwrap();
    let original = FatiguePlate::assemble(
        &mesh,
        vec![section(&f); 2],
        &[],
        &options(mean * H),
        vec![gauge(&f, 0, z, [1.0, 1.0])],
        config(100),
        &mut || false,
    )
    .unwrap();
    assert_eq!(plate.model().k.to_dense(), original.model().k.to_dense());
    assert_eq!(plate.model().m.to_dense(), original.model().m.to_dense());
    let positive = field(&mesh, [0.0, 0.0, amplitude / (shear_modulus * z)]);
    let negative: Vec<_> = positive.iter().map(|v| -*v).collect();
    let mut observed = plate.observe(&negative, &mut || false).unwrap();
    for (i, tensor) in observed.findley_stresses_pa().iter().enumerate() {
        close(tensor[0], mean, 1e-12);
        close(tensor[1], mean, 1e-12);
        close(
            tensor[3],
            if i == 0 { -amplitude } else { amplitude },
            1e-12,
        );
        assert_eq!([tensor[2], tensor[4], tensor[5]], [0.0; 3]);
    }
    observed.push(&positive, &mut || false).unwrap();
    observed.push(&negative, &mut || false).unwrap();
    let preview = observed.preview_findley(&mut || false).unwrap();
    assert!(observed.findley_results().iter().all(Option::is_none));
    for result in &preview {
        close(result.shear_amplitude_pa, amplitude, 1e-12);
        close(result.maximum_normal_stress_pa, mean, 1e-12);
        close(result.usage_factor, (amplitude + K * mean) / LIMIT, 1e-12);
    }
    let normal = observed.summaries();
    assert_eq!(observed.finish_findley(&mut || false).unwrap(), preview);
    assert!(!observed.is_finished());
    assert!(matches!(
        observed.push(&positive, &mut || false),
        Err(PlateFatigueError::FindleyGauge {
            error: FindleyError::Finished,
            ..
        })
    ));
    assert_eq!(observed.summaries(), normal);
    assert!(observed.findley_results().iter().all(Option::is_some));
    observed.restart_findley(&mut || false).unwrap();
    assert_eq!(observed.summaries(), normal);
    assert_eq!(observed.findley_sample_count(0), Some(1));
    assert!(observed.findley_results().iter().all(Option::is_none));
    observed.push(&positive, &mut || false).unwrap();
    observed.push(&negative, &mut || false).unwrap();
    let summaries = observed.finish(&mut || false).unwrap();
    close(
        summaries[0].miner_damage,
        2.0 * amplitude / LIFE_STRESS,
        1e-11,
    );
    for result in observed.findley_results().into_iter().flatten() {
        close(result.usage_factor, (amplitude + K * mean) / LIMIT, 1e-12);
    }
    assert_eq!(observed.finish(&mut || false).unwrap(), summaries);
}

#[test]
fn rotating_shear_path_on_a_tilted_plane_recovers_its_two_dimensional_radius() {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.01, 1, 1);
    let amplitude = 1e6;
    let z = H / 2.0;
    let g = E / (2.0 * (1.0 + NU));
    let gauges = vec![findley_gauge(&f, 0, z, vec![[1.0, 0.0, 1.0]], 5)];
    let plate = FatiguePlate::assemble_with_findley(
        &mesh,
        vec![section(&f); 2],
        &[],
        &options(0.0),
        vec![],
        gauges.clone(),
        config(100),
        budget(&gauges),
        &mut || false,
    )
    .unwrap();
    // sigma_xx=2A cos(theta), tau_xy=sqrt(2) A sin(theta), sigma_yy=0.
    // On n=(1,0,1)/sqrt(2), shear traction traces a circle of radius A
    // and maximum normal stress is A. Invert independent plane-stress compliance.
    let states: Vec<_> = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0)]
        .into_iter()
        .map(|(c, s)| {
            field(
                &mesh,
                [
                    2.0 * amplitude * c / (E * z),
                    -2.0 * NU * amplitude * c / (E * z),
                    2.0_f64.sqrt() * amplitude * s / (g * z),
                ],
            )
        })
        .collect();
    let mut observed = plate.observe(&states[0], &mut || false).unwrap();
    for state in &states[1..] {
        observed.push(state, &mut || false).unwrap();
    }
    assert!(observed.finish(&mut || false).unwrap().is_empty());
    let result = observed.findley_results()[0].unwrap();
    close(result.shear_amplitude_pa, amplitude, 1e-12);
    close(result.maximum_normal_stress_pa, amplitude, 1e-12);
    close(result.findley_parameter_pa, amplitude * (1.0 + K), 1e-12);
    assert_eq!(result.samples, 5);
    assert_eq!(result.evaluated_planes, 1);
}

#[test]
fn late_budget_cancellation_and_open_gauge_finish_restore_both_models() {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.01, 1, 1);
    let make = |second_samples| {
        let gauges = vec![
            findley_gauge(&f, 0, H / 2.0, vec![[1.0, 0.0, 0.0]], 8),
            findley_gauge(&f, 1, H / 2.0, vec![[1.0, 0.0, 0.0]], second_samples),
        ];
        FatiguePlate::assemble_with_findley(
            &mesh,
            vec![section(&f); 2],
            &[],
            &options(0.0),
            vec![gauge(&f, 0, H / 2.0, [1.0, 1.0])],
            gauges.clone(),
            config(100),
            budget(&gauges),
            &mut || false,
        )
        .unwrap()
    };
    let plate = make(3);
    let positive = field(&mesh, [0.0, 0.0, 0.03]);
    let negative: Vec<_> = positive.iter().map(|v| -*v).collect();
    let mut observed = plate.observe(&negative, &mut || false).unwrap();
    observed.push(&positive, &mut || false).unwrap();
    let mut trial = observed.clone();
    let mut polls = 0;
    trial
        .push(&negative, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let old = observed.summaries();
    let tensors = observed.findley_stresses_pa().to_vec();
    let mut seen = 0;
    assert_eq!(
        observed.push(&negative, &mut || {
            seen += 1;
            seen == polls
        }),
        Err(PlateFatigueError::Cancelled),
    );
    assert_eq!(observed.summaries(), old);
    assert_eq!(observed.findley_stresses_pa(), tensors.as_slice());
    assert_eq!(observed.findley_sample_count(0), Some(2));
    assert_eq!(observed.findley_sample_count(1), Some(2));
    observed.push(&negative, &mut || false).unwrap();
    let old = observed.summaries();
    assert!(matches!(
        observed.push(&positive, &mut || false),
        Err(PlateFatigueError::FindleyGauge {
            gauge: 1,
            error: FindleyError::Budget
        })
    ));
    assert_eq!(observed.summaries(), old);
    assert_eq!(observed.findley_sample_count(0), Some(3));
    assert_eq!(observed.findley_sample_count(1), Some(3));

    let plate = make(8);
    let mut open = plate.observe(&negative, &mut || false).unwrap();
    open.push(&positive, &mut || false).unwrap();
    let mut almost = negative.clone();
    // Node 2 belongs only to triangle 1. The first gauge closes exactly,
    // while the second gauge retains a genuine full-tensor closure defect.
    almost[3 * 2] += 1e-7;
    open.push(&almost, &mut || false).unwrap();
    let normal = open.summaries();
    assert!(matches!(
        open.finish(&mut || false),
        Err(PlateFatigueError::FindleyGauge {
            gauge: 1,
            error: FindleyError::NotClosed { .. }
        })
    ));
    assert!(!open.is_finished());
    assert_eq!(open.summaries(), normal);
    assert!(open.findley_results().iter().all(Option::is_none));
    open.push(&negative, &mut || false).unwrap();
    open.finish(&mut || false).unwrap();
    assert!(open.findley_results().iter().all(Option::is_some));
}

fn sourced(
    youngs: f64,
    humidity: f64,
) -> (Arc<IsotropicElasticStatePoint>, Arc<FindleyCalibration>) {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let domain = ValidityDomain::unconstrained()
        .with_quantity("T", absolute, 290.0, 310.0)
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
    let point = QueryPoint::new()
        .with_quantity("T", absolute, 300.0)
        .unwrap()
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), humidity)
        .unwrap();
    let mut rows = ClaimSet::new();
    let mut values = vec![
        (PropertyKey::new("density", Density::DIMS), 7800.0),
        (PropertyKey::new("young_modulus", Pressure::DIMS), youngs),
        (PropertyKey::new("poisson_ratio", Dims::NONE), NU),
    ];
    values.extend(
        FindleyProperty::ALL
            .into_iter()
            .map(|p| (p.property_key(), [K, LIMIT][p.index()])),
    );
    for (key, value) in values {
        rows.insert_claim(PropertyClaim {
            value: PropertyValue::Scalar {
                value,
                dims: key.dims(),
            },
            key,
            validity: domain.clone(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic common plate and Findley fit".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic plate alloy".into(),
            phase: "solid".into(),
            process: "synthetic fixed process".into(),
            revision: 0,
        },
        rows,
        vec![],
    )
    .unwrap();
    let material = Arc::new(
        resolve_isotropic_elastic_state_point(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    );
    let calibration = Arc::new(
        resolve_findley(
            &card,
            &point,
            &FindleySourceRequest {
                temperature_axis: "T".into(),
                surface_and_geometry: "synthetic smooth plate".into(),
                environment: "synthetic fixed humidity".into(),
                load_basis: "synthetic common periodic fit".into(),
            },
            [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 2],
        )
        .unwrap(),
    );
    (material, calibration)
}

#[test]
fn original_findley_source_context_and_aggregate_admission_are_enforced() {
    let f = fixture(E, 0.2);
    let (material, calibration) = sourced(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.01, 1, 1);
    let mut g = findley_gauge(&f, 0, H / 2.0, vec![[1.0, 0.0, 0.0]], 5);
    g.calibration = calibration.clone();
    let sections = vec![
        PlateFatigueSection {
            material: material.clone(),
            thickness_m: H
        };
        2
    ];
    let original = budget(std::slice::from_ref(&g));
    let make = |sections: Vec<PlateFatigueSection>, cap| {
        FatiguePlate::assemble_with_findley(
            &mesh,
            sections,
            &[],
            &options(0.0),
            vec![],
            vec![g.clone()],
            config(20),
            cap,
            &mut || false,
        )
    };
    let plate = make(sections.clone(), original).unwrap();
    assert_eq!(
        plate
            .findley_gauge(0)
            .unwrap()
            .calibration
            .source()
            .unwrap()
            .card_identity(),
        material.resolved().card_identity(),
    );
    let mut too_small = original;
    too_small.maximum_total_samples -= 1;
    assert!(matches!(
        make(sections, too_small),
        Err(PlateFatigueError::Budget)
    ));
    for (youngs, humidity) in [(E, 0.7), (0.5 * E, 0.2)] {
        let (other, _) = sourced(youngs, humidity);
        assert!(matches!(
            make(
                vec![
                    PlateFatigueSection {
                        material: other,
                        thickness_m: H
                    };
                    2
                ],
                original
            ),
            Err(PlateFatigueError::FindleySourceMismatch { gauge: 0, .. }),
        ));
    }
}
