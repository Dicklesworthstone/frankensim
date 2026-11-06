//! Cyclic plastic frame deformation feeding the shared bounded strain-life owner.
//! This authored calibration is an implementation demonstration, not measured life.
//! Run: cargo run -p fs-solid --example frame_fatigue -- --demo 20
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::elastic::IsotropicElastic;
use fs_material::fatigue::strain_life::{CmbParameters, StrainLifeCurve, StrainLifeMeanStress};
use fs_material::fatigue::{FatigueScope, RainflowLimits};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::reduction::UniaxialJ2;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::fatigue::{FrameFatigueConfig, FrameFatiguePath, FrameFiberFatigueGauge};
use fs_solid::{Frame2d, Geometry, SectionModel};

fn run(cycles: usize, cx: &Cx<'_>) -> Result<(), Box<dyn std::error::Error>> {
    let elastic = IsotropicElastic::new(200e9, 0.3, 0.04)?;
    let law = FiberLaw::J2(UniaxialJ2::new(J2Plasticity::new_mixed(
        elastic, 200e6, 2e9, 18e9,
    )?)?);
    let mut frame = Frame2d::new(Geometry::Linear);
    let left = frame.add_node(0.0, 0.0);
    let right = frame.add_node(1.0, 0.0);
    frame.fix(left, [true; 3]);
    frame.fix(right, [false, true, true]);
    frame.add_element(left, right, 3, 0.0, &|| {
        SectionModel::Fiber(Section::new(
            [-0.1, 0.1]
                .map(|y| Fiber {
                    y,
                    area: 1e-3,
                    law: law.clone(),
                })
                .to_vec(),
        ))
    })?;
    let curve = Arc::new(StrainLifeCurve::new(
        FatigueScope {
            source_reference: "authored synthetic demonstration; not measured".into(),
            material: "synthetic cyclic metal".into(),
            process: "synthetic fixed condition".into(),
            surface_and_geometry: "axial coupon".into(),
            environment: "synthetic dry atmosphere".into(),
            load_basis: "total mechanical axial strain amplitude".into(),
            temperature_k: 300.0,
        },
        CmbParameters {
            youngs_modulus_pa: 200e9,
            fatigue_strength_coefficient_pa: 1e9,
            fatigue_strength_exponent: -0.5,
            fatigue_ductility_coefficient: 0.005,
            fatigue_ductility_exponent: -0.5,
        },
        1.0,
        1e12,
    )?);
    let gauge = FrameFiberFatigueGauge {
        element: 0,
        section: 0,
        fiber: 0,
        material_source: None,
        curve,
        mean_stress: StrainLifeMeanStress::Morrow,
    };
    let mut run = FrameFatiguePath::new(
        frame,
        vec![gauge],
        FrameFatigueConfig {
            history: RainflowLimits {
                maximum_samples: (2 * cycles + 2) as u64,
                maximum_residue: 64,
            },
            maximum_gauges: 1,
            maximum_total_residue: 64,
            maximum_absolute_axial_strain: 0.04,
        },
        cx,
    )?;
    let pattern = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    run.displacement_control(&pattern, 3, &[0.003], cx)?;
    for _ in 0..cycles {
        // Native plastic state and unresolved rainflow reversals survive blocks.
        run.displacement_control(&pattern, 3, &[-0.003, 0.003], cx)?;
    }
    let sample = run.samples()[0];
    let committed = run.summaries();
    let provisional = run.preview_finish(cx)?;
    assert_eq!(run.summaries(), committed);
    assert!(!run.is_finished());
    let completed = run.finish(cx)?;
    assert_eq!(provisional, completed);
    println!("synthetic cyclic frame; no service-life qualification");
    println!("cycles={cycles} accepted_samples={}", run.sample_count());
    println!(
        "axial_strain={:.6e} committed_plastic_stress_pa={:.6e}",
        sample.mechanical_strain, sample.normal_stress_pa
    );
    println!(
        "committed_damage={:.8e} damage_if_stopped_now={:.8e}",
        committed.fibers[0].miner_damage, completed.fibers[0].miner_damage
    );
    println!(
        "counted_full_cycles={} counted_half_cycles={}",
        completed.fibers[0].full_cycles, completed.fibers[0].half_cycles
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) != Some("--demo") || args.len() > 2 {
        return Err(
            "usage: frame_fatigue --demo [cycles: 1..5000]; synthetic calibration only".into(),
        );
    }
    let cycles = args
        .get(1)
        .map(|v| v.parse::<usize>())
        .transpose()?
        .unwrap_or(20);
    if !(1..=5000).contains(&cycles) {
        return Err("cycles must be in 1..=5000".into());
    }
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(cycles, &cx)
    })
}
