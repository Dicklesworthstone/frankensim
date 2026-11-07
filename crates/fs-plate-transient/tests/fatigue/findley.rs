//! Native trajectory integration of representative multiaxial fatigue windows.
use super::*;
use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyLimits, StressClosureTolerance,
};
use fs_plate::fatigue::{PlateFindleyBudget, PlateFindleyGauge};
use std::sync::Arc;

const K_FINDLEY: f64 = 0.3;
const LIMIT_PA: f64 = 1e6;

fn shear_plate(samples: [u64; 2]) -> FatiguePlate {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.02, 2, 2);
    let element = mesh.tris.iter().position(|t| t.contains(&4)).unwrap();
    let state = f.material.resolved().material();
    let calibration = Arc::new(
        FindleyCalibration::new(
            FatigueScope {
                source_reference: "synthetic Findley vibration fixture; not measured".into(),
                material: state.chemistry.clone(),
                process: state.process.clone(),
                surface_and_geometry: "synthetic homogeneous smooth plate".into(),
                environment: "synthetic fixed environment".into(),
                load_basis: "resolved closed elastic stress block".into(),
                temperature_k: 300.0,
            },
            K_FINDLEY,
            LIMIT_PA,
        )
        .unwrap(),
    );
    let findley = samples
        .into_iter()
        .enumerate()
        .map(|(i, count)| PlateFindleyGauge {
            element,
            surface_point: 0,
            depth_m: if i == 0 { H / 2.0 } else { -H / 2.0 },
            temperature_axis: "T".into(),
            calibration: calibration.clone(),
            normals_in_mesh_frame: vec![[1.0, 0.0, 1.0]],
            limits: FindleyLimits {
                maximum_samples: count,
                maximum_planes: 1,
                maximum_projection_work: count,
                maximum_circle_work: 100_000,
            },
            closure: StressClosureTolerance {
                absolute_stress_pa: 1e-5,
                relative: 1e-8,
            },
        })
        .collect();
    FatiguePlate::assemble_with_findley(
        &mesh,
        vec![section(&f); mesh.tris.len()],
        &PlateMesh::rectangle_boundary(2, 2),
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
        vec![gauge(&f, element, H / 2.0, [1.0, 0.0])],
        findley,
        config(1000),
        PlateFindleyBudget {
            maximum_total_samples: samples.iter().sum(),
            maximum_total_planes: 2,
            maximum_total_projection_work: samples.iter().sum(),
            maximum_total_circle_work: 200_000,
        },
        &mut || false,
    )
    .unwrap()
}

fn oscillator(plate: &FatiguePlate) -> (Vec<f64>, Vec<f64>, f64, f64) {
    let p = plate.model();
    assert_eq!(p.free, 3);
    let w = p.dof_map[12].unwrap();
    let omega = (p.k.to_dense()[3 * w + w] / p.m.to_dense()[3 * w + w]).sqrt();
    let h = 2.0 * (std::f64::consts::PI / 64.0).tan() / omega;
    let mut q = vec![0.0; 3];
    q[w] = 1e-7;
    (q, vec![0.0; 3], omega, h)
}

#[test]
fn g1_two_native_periods_keep_miner_continuous_across_findley_windows() {
    let plate = shear_plate([65, 65]);
    let (q, v, _, h) = oscillator(&plate);
    let stepper = method(3, h);
    let mut run = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.0),
        stepper,
        initial(&q, &v),
        &mut || false,
    )
    .unwrap();
    let initial_tensor = run.observation().findley_stresses_pa()[0];
    // For n=(1,0,1)/sqrt(2), plane stress gives sigma_n=sigma_xx/2.
    // A proportional signed cycle has a straight shear-traction path with
    // radius sqrt(sigma_xx^2/4 + tau_xy^2/2), regardless of sample ordering.
    let shear = (initial_tensor[0].powi(2) / 4.0 + initial_tensor[3].powi(2) / 2.0).sqrt();
    let maximum_normal = initial_tensor[0].abs() / 2.0;
    let expected_usage = (shear + K_FINDLEY * maximum_normal) / LIMIT_PA;
    assert!(expected_usage > 0.0);
    let normal_amplitude = run.observation().stresses_pa()[0].abs();
    let mut native = run.state().clone();
    let mut completed_first = None;
    for period in 0..2 {
        for i in 1..=64 {
            run.advance(&mut || false).unwrap();
            stepper
                .step(&mut native, &dynamics(&plate, 0.0), &v)
                .unwrap();
            native.history.clear();
            assert_eq!(run.state(), &native);
            if i == 32 {
                let before = run.summaries();
                assert!(run.preview_findley(&mut || false).is_err());
                assert_eq!(run.summaries(), before);
                assert!(run.findley_results().iter().all(Option::is_none));
            }
        }
        let motion = run.state().clone();
        let miner = run.summaries();
        let preview = run.preview_findley(&mut || false).unwrap();
        assert_eq!(run.state(), &motion);
        assert_eq!(run.summaries(), miner);
        assert!(run.findley_results().iter().all(Option::is_none));
        let result = run.finish_findley(&mut || false).unwrap();
        assert_eq!(preview, result);
        assert_eq!(run.state(), &motion);
        assert_eq!(run.summaries(), miner);
        for gauge in &result {
            assert_eq!(gauge.samples, 65);
            close(gauge.shear_amplitude_pa, shear, 5e-8);
            close(gauge.maximum_normal_stress_pa, maximum_normal, 5e-8);
            close(gauge.usage_factor, expected_usage, 5e-8);
        }
        assert_eq!(run.finish_findley(&mut || false).unwrap(), result);
        if period == 0 {
            completed_first = Some(result);
            assert!(run.advance(&mut || false).is_err());
            assert_eq!(run.state(), &motion);
            assert_eq!(run.summaries(), miner);
            assert!(run.restart_findley(&mut || true).is_err());
            assert_eq!(
                run.findley_results()[0],
                Some(completed_first.as_ref().unwrap()[0])
            );
            run.restart_findley(&mut || false).unwrap();
            assert_eq!(run.state(), &motion);
            assert_eq!(run.summaries(), miner);
            assert_eq!(run.observation().findley_sample_count(0), Some(1));
            assert!(run.findley_results().iter().all(Option::is_none));
        }
    }
    let result = run.finish(&mut || false).unwrap();
    close(
        result[0].miner_damage,
        2.0 * normal_amplitude / LIFE_STRESS,
        5e-8,
    );
    assert!(completed_first.is_some());
    assert!(run.findley_results().iter().all(Option::is_some));
}

#[test]
fn g4_late_findley_gauge_budget_or_cancellation_rolls_back_native_motion() {
    let plate = shear_plate([6, 3]);
    let (q, v, _, h) = oscillator(&plate);
    let mut run = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.0),
        method(3, h),
        initial(&q, &v),
        &mut || false,
    )
    .unwrap();
    let old = run.state().clone();
    let miner = run.summaries();
    let tensor = run.observation().findley_stresses_pa().to_vec();
    let mut trial = run.clone();
    let mut calls = 0;
    trial
        .advance(&mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(
        run.advance(&mut || {
            seen += 1;
            seen == calls
        })
        .is_err()
    );
    assert_eq!(run.state(), &old);
    assert_eq!(run.summaries(), miner);
    assert_eq!(run.observation().findley_stresses_pa(), tensor.as_slice());
    for index in 0..2 {
        assert_eq!(run.observation().findley_sample_count(index), Some(1));
    }
    for _ in 0..2 {
        run.advance(&mut || false).unwrap();
    }
    let old = run.state().clone();
    let miner = run.summaries();
    let tensor = run.observation().findley_stresses_pa().to_vec();
    // Gauge zero can append; gauge one cannot. Neither its accepted history
    // nor the motion or ordinary rainflow counters may retain that late suffix.
    assert!(run.advance(&mut || false).is_err());
    assert_eq!(run.state(), &old);
    assert_eq!(run.summaries(), miner);
    assert_eq!(run.observation().findley_stresses_pa(), tensor.as_slice());
    for index in 0..2 {
        assert_eq!(run.observation().findley_sample_count(index), Some(3));
    }
}

#[test]
fn g3_damped_motion_is_not_silently_declared_a_periodic_stress_cycle() {
    let plate = shear_plate([65, 65]);
    let (q, v, omega, h) = oscillator(&plate);
    let mut run = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.1 * omega),
        method(3, h),
        initial(&q, &v),
        &mut || false,
    )
    .unwrap();
    for _ in 0..64 {
        run.advance(&mut || false).unwrap();
    }
    let motion = run.state().clone();
    let miner = run.summaries();
    assert!(run.finish_findley(&mut || false).is_err());
    assert_eq!(run.state(), &motion);
    assert_eq!(run.summaries(), miner);
    assert!(run.findley_results().iter().all(Option::is_none));
    // A caller may explicitly finish conditioning and begin a new window;
    // the complete nonperiodic history still belongs to the Miner observer.
    run.restart_findley(&mut || false).unwrap();
    assert_eq!(run.state(), &motion);
    assert_eq!(run.summaries(), miner);
    assert_eq!(run.observation().findley_sample_count(0), Some(1));
}
