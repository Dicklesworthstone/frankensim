//! cargo run -p fs-solid --example fatigue_crack_growth -- --demo 1000000
//!
//! SYNTHETIC Paris-law demonstration; coefficients are not measured material data.
//! Three native solves resolve one tensile 2 -> 8 -> 2 MPa cycle. The caller then
//! explicitly declares identical quasistatic cycles with fixed material, thermal
//! state, supports and uncracked mechanics. Constant Y=1 is a declared geometry
//! approximation: a denotes crack half-length; growth does not change stiffness.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGrowthLimits, ParisDomain, ParisLaw, ParisParameters,
};
use fs_solid::fatigue::crack_growth::{
    CrackStressTolerance, TetCrackGrowthConfig, TetCrackGrowthGauge, TetCrackGrowthPath,
};
use fs_solid::linear3::{
    TetElasticMaterial, TetLinearElasticProblem, TetMaterialField, TetThermalStrainField,
    TetThermalStrainState,
};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const INITIAL_LENGTH_M: f64 = 1e-3;
const REFERENCE_DELTA_K: f64 = 1e6;
const REFERENCE_RATE: f64 = 1e-10;
const MINIMUM_STRESS_PA: f64 = 2e6;
const MAXIMUM_STRESS_PA: f64 = 8e6;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 || args[0] != "--demo" {
        return Err(
            "usage: fatigue_crack_growth --demo CYCLES (1..=1000000); synthetic data only".into(),
        );
    }
    let cycles: u64 = args[1].parse()?;
    if !(1..=1_000_000).contains(&cycles) {
        return Err("CYCLES must lie in 1..=1000000".into());
    }
    let nodes = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let tetrahedra = [[0, 1, 2, 3]];
    // Only node 1's x displacement is free; nu=0 gives a uniaxial stress state.
    let fixed_dofs: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
    let material = TetElasticMaterial::try_new(7800.0, 200e9, 0.0, [104; 32])?;
    let thermal = TetThermalStrainState::try_new(
        300.0,
        300.0,
        [0.0; 6],
        [104; 32],
        ContentHash([105; 32]),
        0.01,
    )?;
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &tetrahedra,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed_dofs,
        budget: Default::default(),
    };
    let law = Arc::new(ParisLaw::new(
        FatigueScope {
            source_reference: "SYNTHETIC m=2 arithmetic fixture; not measured".into(),
            material: "synthetic elastic demonstration".into(),
            process: "synthetic".into(),
            surface_and_geometry: "declared constant Y=1; a is crack half-length".into(),
            environment: "synthetic fixed conditions".into(),
            load_basis: "tensile Mode-I full stress-intensity range".into(),
            temperature_k: 300.0,
        },
        ParisParameters {
            reference_delta_k_pa_sqrt_m: REFERENCE_DELTA_K,
            reference_rate_m_per_cycle: REFERENCE_RATE,
            exponent: 2.0,
        },
        ParisDomain {
            minimum_delta_k_pa_sqrt_m: 1.0,
            maximum_delta_k_pa_sqrt_m: 1e8,
            minimum_load_ratio: 0.0,
            maximum_load_ratio: 0.95,
            maximum_k_pa_sqrt_m: 1e9,
        },
    )?);
    let allowance = CrackStressTolerance {
        absolute_stress_pa: 1e-5,
        relative: 1e-9,
    };
    let gauge = TetCrackGrowthGauge {
        element: 0,
        normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: [104; 32],
        material_source: None,
        law,
        geometry: ConstantGeometryCrack {
            geometry_factor: 1.0,
            minimum_length_m: 1e-5,
            maximum_length_m: 0.01,
        },
        initial_length_m: INITIAL_LENGTH_M,
        limits: CrackGrowthLimits {
            maximum_blocks: 1,
            maximum_cycles: cycles as f64,
        },
        closure: allowance,
        mode_i: allowance,
    };
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 104, kernel_id: 104, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        let geometry = TetLoadGeometry3::new(&nodes, &tetrahedra, problem.budget, 4, &cx)?;
        let mut path = TetCrackGrowthPath::new(
            &problem, TetThermalStrainField::Uniform(&thermal), ContentHash([106; 32]),
            vec![gauge], Default::default(),
            TetCrackGrowthConfig {
                maximum_absolute_mandel_strain: 0.01,
                maximum_gauges: 1,
                maximum_samples: 3,
            }, &cx,
        )?;
        for (sample, stress) in [MINIMUM_STRESS_PA, MAXIMUM_STRESS_PA, MINIMUM_STRESS_PA]
            .into_iter().enumerate()
        {
            // Unit tetrahedron: V=1/6 m^3 and stiffness E/6, so F=stress/6 N.
            let load = geometry.assemble(&[], VolumeLoad3::None, &[(3, stress / 6.0)], &cx)?;
            path.advance(sample as f64, &load, &cx)?;
        }
        let steps = path.close_repeated_cycle(cycles, &cx)?;
        let step = &steps[0];
        path.finish(&cx)?;
        // Independent m=2 solution: da/dN = rate * (Y*DeltaStress/Kref)^2 * pi*a.
        let normalized_range = (MAXIMUM_STRESS_PA - MINIMUM_STRESS_PA) / REFERENCE_DELTA_K;
        let expected = INITIAL_LENGTH_M
            * (REFERENCE_RATE * normalized_range.powi(2) * std::f64::consts::PI * cycles as f64).exp();
        println!("SYNTHETIC demonstration, not validated component life");
        println!("Declared identical quasistatic cycles; uncracked mechanics and Y=1 stay fixed.");
        println!("native_samples={} repeated_full_cycles={cycles}", path.observed_samples());
        println!("normal_stress_min_pa={MINIMUM_STRESS_PA:.6e} normal_stress_max_pa={MAXIMUM_STRESS_PA:.6e}");
        println!("initial_crack_length_m={:.12e} final_crack_length_m={:.12e}",
            step.before.crack_length_m, step.after.crack_length_m);
        println!("initial_delta_k_pa_sqrt_m={:.12e} final_delta_k_pa_sqrt_m={:.12e}",
            step.initial_delta_k_pa_sqrt_m, step.final_delta_k_pa_sqrt_m);
        println!("analytic_m2_length_m={expected:.12e} relative_error={:.6e}",
            (step.after.crack_length_m - expected).abs() / expected);
        Ok(())
    })
}
