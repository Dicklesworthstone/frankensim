//! Actual P1 load/stress solves feeding the empirical fatigue owner.
//! Every S-N curve below is SYNTHETIC, not a measured material life claim.
#[path = "tet10/support.rs"]
mod support;

use std::sync::Arc;
use support::{block, with_cx};
use fs_blake3::ContentHash;
use fs_material::fatigue::{FatigueError, FatigueScope, MeanStressCorrection, RainflowLimits, SnCurve, SnPoint};
use fs_solid::fatigue::{TetFatigueConfig, TetFatigueError, TetFatigueGauge, TetFatiguePath};
use fs_solid::linear3::{TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState};
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};

fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(7800.0, 200e9, 0.25, [7; 32]).unwrap()
}
fn thermal() -> TetThermalStrainState {
    TetThermalStrainState::try_new(293.15, 293.15, [0.0; 6], [7; 32], ContentHash([8; 32]), 0.01).unwrap()
}
fn curve(maximum: f64, temperature: f64) -> Arc<SnCurve> {
    Arc::new(SnCurve::new(FatigueScope {
        source_reference: "synthetic N(a)=1e9/a; not measured".into(),
        material: "synthetic fixture".into(), process: "fixture".into(),
        surface_and_geometry: "smooth axial fixture".into(), environment: "fixture".into(),
        load_basis: "fully reversed axial stress".into(), temperature_k: temperature,
    }, vec![SnPoint { amplitude_pa: 1e6, cycles_to_failure: 1000.0 },
        SnPoint { amplitude_pa: maximum, cycles_to_failure: 1e9 / maximum }], None).unwrap())
}
fn gauge(element: usize) -> TetFatigueGauge {
    TetFatigueGauge { element, normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: [7; 32], curve: curve(10e6, 293.15),
        mean_stress: MeanStressCorrection::Uncorrected }
}
fn config() -> TetFatigueConfig {
    TetFatigueConfig { solver: Default::default(), maximum_absolute_mandel_strain: 0.01,
        history: RainflowLimits { maximum_samples: 100, maximum_residue: 32 },
        maximum_gauges: 8, maximum_total_residue: 256 }
}
fn rollers(nodes: &[[f64; 3]]) -> Vec<usize> {
    nodes.iter().enumerate().flat_map(|(i, p)|
        (0..3).filter_map(move |c| (p[c] == 0.0).then_some(3*i+c))).collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:e} != {expected:e}");
}

#[test]
fn native_tension_compression_retains_signed_cycles_and_original_mechanics() {
    for refinement in [1, 2] {
        let mesh = block([refinement; 3], [1.0; 3]);
        let material = material(); let thermal = thermal(); let fixed = rollers(&mesh.nodes);
        let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
            materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget: Default::default() };
        with_cx(|cx| {
            let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 1000, cx).unwrap();
            let right: Vec<_> = geometry.faces().iter().enumerate()
                .filter_map(|(i, face)| (face.outward_normal[0] > 0.99).then_some(i)).collect();
            let mut path = TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
                ContentHash([9; 32]), vec![gauge(0), gauge(mesh.tets.len()-1)], config(), cx).unwrap();
            for (sample, stress) in [2e6, -2e6, 2e6].into_iter().enumerate() {
                let load = geometry.assemble(&[SurfaceLoad3::Traction { faces: &right,
                    pascals: [stress, 0.0, 0.0] }], VolumeLoad3::None, &[], cx).unwrap();
                let step = path.advance(sample as f64, &load, cx).unwrap();
                for signed in step.signed_stress_pa { close(signed, stress, 0.01); }
                for (point, displacement) in mesh.nodes.iter().zip(step.mechanics.displacement.displacement_m()) {
                    close(displacement[0], stress / 200e9 * point[0], 1e-12);
                }
                let original = load.solve(&problem, TetThermalStrainField::Uniform(&thermal),
                    Default::default(), ContentHash([9; 32]), 0.01, cx).unwrap();
                assert_eq!(step.mechanics.stress, original.stress);
                assert_eq!(step.mechanics.displacement, original.displacement);
            }
            let result = path.finish(cx).unwrap();
            for item in &result {
                assert_eq!(item.samples, 3);
                assert_eq!(item.full_cycles * 2 + item.half_cycles, 2);
                close(item.miner_damage, 0.002, 1e-10);
            }
            assert_eq!(path.finish(cx).unwrap(), result);
            assert_eq!(path.curve(0).unwrap().scope().source_reference, "synthetic N(a)=1e9/a; not measured");
        });
    }
}

#[test]
fn actual_positive_mean_stress_uses_the_named_goodman_policy() {
    let mesh = block([1; 3], [1.0; 3]);
    let material = material(); let thermal = thermal(); let fixed = rollers(&mesh.nodes);
    let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget: Default::default() };
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 100, cx).unwrap();
        let right: Vec<_> = geometry.faces().iter().enumerate()
            .filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
        let mut corrected = gauge(0);
        corrected.mean_stress = MeanStressCorrection::GoodmanTensionOnly { ultimate_strength_pa: 20e6 };
        let mut path = TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
            ContentHash([9; 32]), vec![gauge(0), corrected], config(), cx).unwrap();
        for (i, stress) in [2e6, 8e6, 2e6].into_iter().enumerate() {
            let load = geometry.assemble(&[SurfaceLoad3::Traction { faces: &right,
                pascals: [stress, 0.0, 0.0] }], VolumeLoad3::None, &[], cx).unwrap();
            path.advance(i as f64, &load, cx).unwrap();
        }
        let results = path.finish(cx).unwrap();
        close(results[0].miner_damage, 0.003, 1e-10);
        close(results[1].miner_damage, 0.004, 1e-10);
        close(results[1].maximum_corrected_amplitude_pa, 4e6, 0.02);
    });
}

#[test]
fn a_later_gauge_refusal_rolls_back_the_entire_sample_and_finalization() {
    let mesh = block([1; 3], [1.0; 3]);
    let material = material(); let thermal = thermal(); let fixed = rollers(&mesh.nodes);
    let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget: Default::default() };
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 100, cx).unwrap();
        let all: Vec<_> = (0..geometry.faces().len()).collect();
        let mut narrow = gauge(0); narrow.curve = curve(1.5e6, 293.15);
        let mut path = TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
            ContentHash([9; 32]), vec![gauge(0), narrow], config(), cx).unwrap();
        for (i, pressure) in [2e6, -2e6, 2e6].into_iter().enumerate() {
            let load = geometry.assemble(&[SurfaceLoad3::Pressure { faces: &all, pascals: pressure }],
                VolumeLoad3::None, &[], cx).unwrap();
            path.advance(i as f64, &load, cx).unwrap();
        }
        let before = path.summaries();
        let load = geometry.assemble(&[SurfaceLoad3::Pressure { faces: &all, pascals: -2e6 }],
            VolumeLoad3::None, &[], cx).unwrap();
        assert!(matches!(path.advance(3.0, &load, cx), Err(TetFatigueError::Gauge {
            gauge: 1, error: FatigueError::AmplitudeOutsideCurve { .. } })));
        assert_eq!(path.summaries(), before);
        assert_eq!(path.last_time_s(), Some(2.0));
        assert!(matches!(path.finish(cx), Err(TetFatigueError::Gauge { gauge: 1, .. })));
        assert_eq!(path.summaries(), before);
        assert!(path.summaries().iter().all(|s| !s.finished));
    });
}

#[test]
fn source_temperature_and_work_budgets_cannot_be_silently_rebound() {
    let mesh = block([1; 3], [1.0; 3]);
    let material = material(); let thermal = thermal(); let fixed = rollers(&mesh.nodes);
    let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget: Default::default() };
    with_cx(|cx| {
        let mut foreign = gauge(0); foreign.material_state_identity = [88; 32];
        let mut hot = gauge(0); hot.curve = curve(10e6, 350.0);
        for bad in [foreign, hot] {
            assert!(matches!(TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
                ContentHash([9; 32]), vec![bad], config(), cx), Err(TetFatigueError::SourceMismatch { .. })));
        }
        let mut bounded = config(); bounded.maximum_total_residue = 31;
        assert!(matches!(TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
            ContentHash([9; 32]), vec![gauge(0)], bounded, cx), Err(TetFatigueError::Budget)));
        assert!(TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
            ContentHash([9; 32]), vec![gauge(mesh.tets.len())], config(), cx).is_err());
    });
}

#[test]
fn cancelled_or_foreign_mesh_steps_preserve_the_last_accepted_time() {
    let mesh = block([1; 3], [1.0; 3]);
    let foreign = block([1; 3], [2.0; 3]);
    let material = material(); let thermal = thermal(); let fixed = rollers(&mesh.nodes);
    let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget: Default::default() };
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 100, cx).unwrap();
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut path = TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
            ContentHash([9; 32]), vec![gauge(0)], config(), cx).unwrap();
        path.advance(0.0, &load, cx).unwrap();
        let before = path.summaries();
        assert!(path.advance(0.0, &load, cx).is_err());
        let other = TetLoadGeometry3::new(&foreign.nodes, &foreign.tets, problem.budget, 100, cx).unwrap();
        let other_load = other.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        assert!(matches!(path.advance(1.0, &other_load, cx), Err(TetFatigueError::Load(_))));
        let gate = fs_exec::CancelGate::new_clock_free(); gate.request();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let stopped = fs_exec::Cx::new(&gate, arena, fs_exec::StreamKey {
                seed: 58, kernel_id: 58, tile: 0, iteration: 0,
            }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
            assert!(matches!(path.advance(1.0, &load, &stopped), Err(TetFatigueError::Cancelled)));
            assert!(matches!(path.finish(&stopped), Err(TetFatigueError::Cancelled)));
        });
        assert_eq!(path.summaries(), before);
        assert_eq!(path.last_time_s(), Some(0.0));
        path.advance(1.0, &load, cx).unwrap();
        assert_eq!(path.finish(cx).unwrap()[0].miner_damage, 0.0);
        assert!(matches!(path.advance(2.0, &load, cx), Err(TetFatigueError::Finished)));
    });
}
