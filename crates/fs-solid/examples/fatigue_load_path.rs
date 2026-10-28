//! cargo run -p fs-solid --example fatigue_load_path -- --demo 20
//! SYNTHETIC arithmetic demonstration only; no measured material life claim.
//! Pressure here is signed normal traction: negative values pull outward.

use std::sync::Arc;
use fs_blake3::ContentHash;
use fs_material::fatigue::{FatigueScope, MeanStressCorrection, RainflowLimits, SnCurve, SnPoint};
use fs_solid::fatigue::{TetFatigueConfig, TetFatigueGauge, TetFatiguePath};
use fs_solid::linear3::{TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState};
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 || args[0] != "--demo" {
        return Err("usage: fatigue_load_path --demo CYCLES (1..=5000); synthetic data only".into());
    }
    let cycles: u32 = args[1].parse()?;
    if cycles == 0 || cycles > 5000 { return Err("CYCLES must lie in 1..=5000".into()); }
    let nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let tets = [[0, 1, 2, 3]];
    let material = TetElasticMaterial::try_new(7800.0, 200e9, 0.25, [7; 32])?;
    let thermal = TetThermalStrainState::try_new(293.15, 293.15, [0.0; 6],
        [7; 32], ContentHash([8; 32]), 0.01)?;
    let problem = TetLinearElasticProblem { nodes_m: &nodes, tetrahedra: &tets,
        materials: TetMaterialField::Uniform(&material), fixed_dofs: &[0, 1, 2, 4, 5, 8],
        budget: Default::default() };
    let curve = Arc::new(SnCurve::new(FatigueScope {
        source_reference: "SYNTHETIC N(amplitude)=1e9/amplitude; not measured".into(),
        material: "demonstration only".into(), process: "synthetic".into(),
        surface_and_geometry: "synthetic fixed-plane observation".into(),
        environment: "synthetic".into(), load_basis: "fully reversed signed normal stress".into(),
        temperature_k: 293.15,
    }, vec![SnPoint { amplitude_pa: 1e6, cycles_to_failure: 1000.0 },
        SnPoint { amplitude_pa: 10e6, cycles_to_failure: 100.0 }], None)?);
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = fs_exec::Cx::new(&gate, arena, fs_exec::StreamKey {
            seed: 58, kernel_id: 58, tile: 0, iteration: 0,
        }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        let geometry = TetLoadGeometry3::new(&nodes, &tets, problem.budget, 4, &cx)?;
        let faces: Vec<_> = (0..geometry.faces().len()).collect();
        let mut path = TetFatiguePath::new(&problem, TetThermalStrainField::Uniform(&thermal),
            ContentHash([9; 32]), vec![TetFatigueGauge {
                element: 0, normal_in_mesh_frame: [1.0, 0.0, 0.0], material_state_identity: [7; 32],
                curve, mean_stress: MeanStressCorrection::Uncorrected,
            }], TetFatigueConfig { solver: Default::default(), maximum_absolute_mandel_strain: 0.01,
                history: RainflowLimits { maximum_samples: u64::from(2*cycles+1), maximum_residue: 16 },
                maximum_gauges: 1, maximum_total_residue: 16 }, &cx)?;
        for sample in 0..=2*cycles {
            let pressure = if sample % 2 == 0 { 2e6 } else { -2e6 };
            let load = geometry.assemble(&[SurfaceLoad3::Pressure { faces: &faces, pascals: pressure }],
                VolumeLoad3::None, &[], &cx)?;
            path.advance(f64::from(sample), &load, &cx)?;
        }
        let summary = path.finish(&cx)?[0];
        println!("SYNTHETIC demonstration, not validated component life");
        println!("samples={} full_cycles={} half_cycles={}", summary.samples, summary.full_cycles, summary.half_cycles);
        println!("miner_damage={:.12} analytic_expected={:.12}", summary.miner_damage, f64::from(cycles)*0.002);
        Ok(())
    })
}
