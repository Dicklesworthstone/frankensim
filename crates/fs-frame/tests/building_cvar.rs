//! G1/G3: material catalog mass selection on native building histories.
//! Material cards are authored numerical fixtures, not measured calibration.

use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_frame::cvar::building::{
    BuildingCvarControls, BuildingCvarError, BuildingMotionSource, building_cvar_from_motions,
    building_cvar_mass_min,
};
use fs_frame::{BuildingFrame, BuildingSpec, RcMember};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fiber::MenegottoPintoSteel;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims, Pressure, QtyAny};
use fs_scenario::ensemble::{RealizationBudget, SpectrumModel, StochasticEnsemble};
use fs_scenario::signal::{Interp, TimeSignal};
use fs_solid::fiber::material::{
    ConcreteFiberProperties, FiberMaterial, RcSection, RcSectionGeometry,
};
use fs_solid::frame2d::FrameMassFormulation;

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn concrete(density: f64) -> FiberMaterial {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", density, Density::DIMS),
        ("fcc", 30e6, Pressure::DIMS),
        ("eps_cc", 0.002, Dims::NONE),
        ("ec", 30e9, Pressure::DIMS),
        ("eps_cu", 0.005, Dims::NONE),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.01,
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored building risk fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "numerical concrete".into(),
            phase: "solid".into(),
            process: "authored fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    FiberMaterial::resolve_concrete(
        &card,
        &QueryPoint::new(),
        &ConcreteFiberProperties {
            density: "density".into(),
            compressive_strength: "fcc".into(),
            peak_strain: "eps_cc".into(),
            initial_modulus: "ec".into(),
            ultimate_strain: "eps_cu".into(),
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn section(member: RcMember, material: &FiberMaterial) -> RcSection {
    RcSection::new(
        RcSectionGeometry {
            depth_m: member.depth,
            width_m: member.width,
            cover_m: 0.04,
            concrete_layers: 16,
            steel_area_per_face_m2: member.steel_area,
            steel_centroid_y_m: 0.5 * member.depth - 0.06,
        },
        material.clone(),
        material.clone(),
        FiberMaterial::steel(
            MenegottoPintoSteel::new(200e9, 450e6, 0.01).unwrap(),
            7800.0,
        )
        .unwrap(),
    )
    .unwrap()
}

fn spec(scale: f64) -> BuildingSpec {
    BuildingSpec {
        stories: 1,
        bays: 1,
        floor_mass: 20_000.0,
        sections_per_member: 3,
        scale,
        ..BuildingSpec::default()
    }
}

fn building(spec: BuildingSpec, density: f64) -> BuildingFrame {
    let material = concrete(density);
    BuildingFrame::with_sections(
        spec,
        section(spec.column, &material),
        section(spec.beam, &material),
        FrameMassFormulation::ConsistentReference,
    )
    .unwrap()
}

fn ensemble() -> StochasticEnsemble {
    StochasticEnsemble {
        name: "material-catalog-risk".into(),
        seed: 29,
        members: 3,
        duration: QtyAny::new(0.16, Dims([0, 0, 1, 0, 0, 0])),
        dt: QtyAny::new(0.02, Dims([0, 0, 1, 0, 0, 0])),
        model: SpectrumModel::KanaiTajimi {
            s0: 1e-4,
            omega_g: QtyAny::new(12.5, Dims([0, 0, -1, 0, 0, 0])),
            zeta_g: 0.6,
        },
    }
}

fn controls() -> BuildingCvarControls {
    BuildingCvarControls {
        beta: 0.5,
        drift_limit: 1.0,
        max_candidates: 3,
        max_members: 3,
        realization_budget: RealizationBudget {
            max_samples: 8,
            max_work: 128,
        },
    }
}

fn physical_mass(spec: BuildingSpec, density: f64) -> f64 {
    let line = |m: RcMember| {
        density * (m.depth * m.width - 2.0 * m.steel_area)
            + 7800.0 * 2.0 * m.steel_area
    };
    spec.stories as f64
        * (spec.floor_mass
            + (spec.bays + 1) as f64 * spec.story_height * line(spec.column) * spec.scale
            + spec.bays as f64 * spec.bay_width * line(spec.beam))
}

#[test]
fn minimum_actual_mass_preserves_sources_and_exact_native_fractional_risk() {
    let ensemble = ensemble();
    let light = building(spec(1.4), 1000.0);
    let catalog = vec![building(spec(0.7), 4000.0), light.clone(), light];
    let identity = catalog[1].material_sections().unwrap()[0]
        .materials()[0].source().unwrap().identity();
    let first_motion = ensemble.realize(0).unwrap();
    let before = catalog[1].clone().run(&first_motion.values, ensemble.dt.value).unwrap();
    let design = with_cx(false, |cx| {
        building_cvar_mass_min(&ensemble, &catalog, controls(), cx).unwrap()
    });
    assert_eq!(design.selected_index, 1);
    assert!(catalog[1].spec().scale > catalog[0].spec().scale);
    assert!(std::ptr::eq(design.selected, &catalog[1]));
    assert!(std::ptr::eq(design.catalog, catalog.as_slice()));
    match &design.motions {
        BuildingMotionSource::Ensemble(retained) => {
            assert!(std::ptr::eq(*retained, &ensemble));
        }
        _ => panic!("ensemble source was replaced"),
    }
    assert_eq!(design.evaluations.len(), catalog.len());
    for (index, (model, evaluation)) in catalog.iter().zip(&design.evaluations).enumerate() {
        let density = if index == 0 { 4000.0 } else { 1000.0 };
        let mass = physical_mass(*model.spec(), density);
        assert!((evaluation.mass_kg - mass).abs() <= 1e-12 * mass);
        let direct: Vec<f64> = (0..ensemble.members)
            .map(|member| {
                let motion = ensemble.realize(member).unwrap();
                model.clone().run(&motion.values, ensemble.dt.value).unwrap().peak_drift_ratio
            })
            .collect();
        assert!(direct.iter().all(|loss| loss.is_finite() && *loss > 0.0));
        assert_eq!(evaluation.peak_drifts, direct);
        assert_eq!(evaluation.tail, fs_robust::empirical_cvar(&direct, 0.5).unwrap());
        let mut sorted = direct;
        sorted.sort_by(f64::total_cmp);
        let expected = (sorted[2] + 0.5 * sorted[1]) / 1.5;
        assert!((evaluation.tail.cvar() - expected).abs() <= 1e-12 * expected);
        assert_eq!(evaluation.tail.boundary_rank(), 2);
        assert_eq!(evaluation.tail.boundary_weight(), 0.5);
    }
    assert!(design.evaluations[1].mass_kg < design.evaluations[0].mass_kg);
    assert_eq!(design.evaluations[1].mass_kg, design.evaluations[2].mass_kg);
    assert_eq!(
        identity,
        design.selected.material_sections().unwrap()[0]
            .materials()[0].source().unwrap().identity()
    );
    assert_eq!(before, catalog[1].clone().run(&first_motion.values, ensemble.dt.value).unwrap());
}

#[test]
fn no_feasible_catalog_returns_the_native_best_risk_without_mutating_models() {
    let ensemble = ensemble();
    let catalog = vec![building(spec(1.0), 2400.0)];
    let direct: Vec<f64> = (0..ensemble.members)
        .map(|member| {
            let motion = ensemble.realize(member).unwrap();
            catalog[0].clone().run(&motion.values, ensemble.dt.value).unwrap().peak_drift_ratio
        })
        .collect();
    let expected = fs_robust::empirical_cvar(&direct, 0.5).unwrap().cvar();
    assert!(expected > 0.0);
    with_cx(false, |cx| {
        match building_cvar_mass_min(
            &ensemble,
            &catalog,
            BuildingCvarControls { drift_limit: 0.0, ..controls() },
            cx,
        ) {
            Err(BuildingCvarError::NoFeasibleCandidate { limit, best_cvar }) => {
                assert_eq!(limit, 0.0);
                assert_eq!(best_cvar, expected);
            }
            _ => panic!("expected a completed infeasible-catalog result"),
        }
    });
    let motion = ensemble.realize(0).unwrap();
    assert_eq!(
        catalog[0].clone().run(&motion.values, ensemble.dt.value).unwrap().peak_drift_ratio,
        direct[0],
    );
}

#[test]
fn cancellation_caps_and_realization_budgets_refuse() {
    let ensemble = ensemble();
    let catalog = vec![building(spec(1.0), 2400.0)];
    with_cx(true, |cx| {
        assert!(matches!(
            building_cvar_mass_min(&ensemble, &catalog, controls(), cx),
            Err(BuildingCvarError::Cancelled)
        ));
    });
    with_cx(false, |cx| {
        for cap in [
            BuildingCvarControls { max_candidates: 0, ..controls() },
            BuildingCvarControls { max_members: 2, ..controls() },
        ] {
            assert!(building_cvar_mass_min(&ensemble, &catalog, cap, cx).is_err());
        }
        for budget in [
            RealizationBudget { max_samples: 7, max_work: 128 },
            RealizationBudget { max_samples: 8, max_work: 1 },
        ] {
            assert!(matches!(
                building_cvar_mass_min(
                    &ensemble,
                    &catalog,
                    BuildingCvarControls { realization_budget: budget, ..controls() },
                    cx,
                ),
                Err(BuildingCvarError::Realization { member: 0, .. })
            ));
        }
    });
}

#[test]
fn missing_materials_and_incompatible_comparisons_are_named() {
    let ensemble = ensemble();
    let source = building(spec(1.0), 2400.0);
    let legacy = vec![BuildingFrame::new(spec(1.0)).unwrap()];
    let incompatible = vec![
        source,
        building(BuildingSpec { floor_mass: 21_000.0, ..spec(1.0) }, 2400.0),
    ];
    with_cx(false, |cx| {
        assert!(matches!(
            building_cvar_mass_min(&ensemble, &legacy, controls(), cx),
            Err(BuildingCvarError::MissingMaterialSections { candidate: 0 })
        ));
        assert!(matches!(
            building_cvar_mass_min(&ensemble, &incompatible, controls(), cx),
            Err(BuildingCvarError::IncompatibleCandidate { candidate: 1 })
        ));
        assert!(building_cvar_mass_min(&ensemble, &[], controls(), cx).is_err());
        assert!(building_cvar_mass_min(
            &ensemble,
            &incompatible[..1],
            BuildingCvarControls { beta: f64::NAN, ..controls() },
            cx,
        ).is_err());
    });
}

fn recorded_motion(dt: f64, values: &[f64]) -> TimeSignal {
    TimeSignal::Table {
        times: (0..values.len()).map(|i| i as f64 * dt).collect(),
        values: values.to_vec(),
        dims: Dims([1, 0, -2, 0, 0, 0]),
        interp: Interp::Linear,
    }
}

#[test]
fn realized_ensemble_and_retained_record_suite_have_identical_evaluations() {
    let ensemble = ensemble();
    let catalog = vec![building(spec(1.0), 2400.0), building(spec(1.3), 1800.0)];
    let records: Vec<_> = (0..ensemble.members)
        .map(|member| {
            let motion = ensemble.realize(member).unwrap();
            TimeSignal::Table {
                times: motion.times,
                values: motion.values,
                dims: Dims([1, 0, -2, 0, 0, 0]),
                interp: Interp::Linear,
            }
        })
        .collect();
    let stochastic = with_cx(false, |cx| {
        building_cvar_mass_min(&ensemble, &catalog, controls(), cx).unwrap()
    });
    let recorded = with_cx(false, |cx| {
        building_cvar_from_motions(
            BuildingMotionSource::Records(&records),
            &catalog,
            controls(),
            cx,
        )
        .unwrap()
    });
    assert_eq!(recorded.selected_index, stochastic.selected_index);
    assert_eq!(recorded.evaluations.len(), stochastic.evaluations.len());
    for (a, b) in recorded.evaluations.iter().zip(&stochastic.evaluations) {
        assert_eq!(a.mass_kg, b.mass_kg);
        assert_eq!(a.peak_drifts, b.peak_drifts);
        assert_eq!(a.tail, b.tail);
    }
    match &recorded.motions {
        BuildingMotionSource::Records(retained) => {
            assert!(std::ptr::eq(*retained, records.as_slice()));
        }
        _ => panic!("recorded source was replaced"),
    }
    match &stochastic.motions {
        BuildingMotionSource::Ensemble(retained) => {
            assert!(std::ptr::eq(*retained, &ensemble));
        }
        _ => panic!("ensemble source was replaced"),
    }
}

#[test]
fn records_preserve_initial_samples_individual_steps_and_equal_record_weight() {
    let catalog = vec![building(spec(1.0), 2400.0)];
    let records = vec![
        recorded_motion(0.01, &[0.3, 0.0, 0.0, 0.0]),
        recorded_motion(0.02, &[0.2, -0.1, 0.05, 0.0, -0.03, 0.01]),
        recorded_motion(0.04, &[0.12, 0.0, -0.04, 0.01, 0.0, -0.02, 0.0, 0.01]),
    ];
    let report = with_cx(false, |cx| {
        building_cvar_from_motions(
            BuildingMotionSource::Records(&records),
            &catalog,
            BuildingCvarControls {
                realization_budget: RealizationBudget { max_samples: 8, max_work: 0 },
                ..controls()
            },
            cx,
        )
        .unwrap()
    });
    let direct: Vec<f64> = records.iter()
        .map(|record| {
            let TimeSignal::Table { times, values, .. } = record else {
                unreachable!()
            };
            catalog[0].clone().run(values, times[1]).unwrap().peak_drift_ratio
        })
        .collect();
    assert!(direct[0] > 0.0);
    // Omitting the nonzero t=0 sample would turn the first motion into rest.
    assert!(catalog[0].clone().run(&[0.0; 3], 0.01).unwrap().peak_drift_ratio < 1e-10);
    assert_eq!(report.evaluations[0].peak_drifts, direct);
    assert_eq!(
        report.evaluations[0].tail,
        fs_robust::empirical_cvar(&direct, 0.5).unwrap(),
    );
    let mut ordered = direct;
    ordered.sort_by(f64::total_cmp);
    let expected = (ordered[2] + 0.5 * ordered[1]) / 1.5;
    assert!((report.evaluations[0].tail.cvar() - expected).abs() <= 1e-12 * expected);
    assert_eq!(report.evaluations[0].tail.boundary_weight(), 0.5);
}

#[test]
fn all_record_shapes_are_preflighted_before_any_candidate_history() {
    let catalog = vec![building(spec(1.0), 2400.0)];
    let valid = recorded_motion(0.02, &[0.1; 8]);
    let acceleration = Dims([1, 0, -2, 0, 0, 0]);
    for defect in 0..9 {
        let invalid = if defect == 8 {
            TimeSignal::Constant(QtyAny::new(0.0, acceleration))
        } else {
            let mut record = valid.clone();
            let TimeSignal::Table { times, values, dims, interp } = &mut record else {
                unreachable!()
            };
            match defect {
                0 => *dims = Dims::NONE,
                1 => *interp = Interp::Hold,
                2 => times[1] = 0.03,
                3 => times.iter_mut().for_each(|time| *time += 1.0),
                4 => values[0] = f64::NAN,
                5 => { times.truncate(1); values.truncate(1); }
                6 => { times.push(0.16); values.push(0.0); }
                7 => { values.pop(); }
                _ => unreachable!(),
            }
            record
        };
        // These finite first-record samples would overflow native forcing.
        // The malformed later record must be refused before any history solve.
        let suite = vec![recorded_motion(0.02, &[f64::MAX; 2]), invalid];
        with_cx(false, |cx| {
            assert!(matches!(
                building_cvar_from_motions(
                    BuildingMotionSource::Records(&suite),
                    &catalog,
                    controls(),
                    cx,
                ),
                Err(BuildingCvarError::Record { member: 1, .. })
            ), "record defect {defect} did not fail suite preflight");
        });
    }
    with_cx(false, |cx| {
        assert!(building_cvar_from_motions(
            BuildingMotionSource::Records(&[]),
            &catalog,
            controls(),
            cx,
        ).is_err());
    });
}
