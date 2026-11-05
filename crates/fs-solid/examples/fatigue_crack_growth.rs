//! cargo run -p fs-solid --example fatigue_crack_growth -- --demo 1000000
//! cargo run -p fs-solid --example fatigue_crack_growth -- --demo 1000000 --varying-geometry
//! cargo run -p fs-solid --example fatigue_crack_growth -- --demo 1000000 --finite-width
//!
//! SYNTHETIC Paris-law demonstration; coefficients are not measured material data.
//! Three native solves resolve one tensile 2 -> 8 -> 2 MPa cycle. The caller then
//! explicitly declares identical quasistatic cycles with fixed material, thermal
//! state, supports and uncracked mechanics. Constant Y=1 and the explicitly
//! synthetic Y(a)=sqrt(a/a0) table retain independent analytic reference solutions.
//! The finite-width mode declares a long elastic plate of full width 0.02 m under
//! remote uniform gross-section tension. Adaptive secant-helper sampling requests
//! a 1e-5 relative-Y interpolation criterion with a 1025-point cap. The native
//! uncracked coupon supplies nominal stress; plate/crack geometry is declared
//! separately. Here a is physical half-crack length, full crack length is 2a,
//! and growth leaves stiffness fixed. The interpolation criterion is analytic;
//! floating-point intervals and material lifetime are not certified.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGeometry, CrackGeometryPoint, CrackGeometryTableContext,
    CrackGrowthLimits, FiniteWidthCentreCrack, ParisDomain, ParisLaw, ParisParameters,
    TabulatedCrackGeometry,
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
const FULL_PLATE_WIDTH_M: f64 = 0.02;
const FINITE_WIDTH_RELATIVE_Y_TOLERANCE: f64 = 1e-5;
const FINITE_WIDTH_MAXIMUM_POINTS: usize = 1025;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(2..=3).contains(&args.len())
        || args[0] != "--demo"
        || args
            .get(2)
            .is_some_and(|flag| flag != "--varying-geometry" && flag != "--finite-width")
    {
        return Err(
            "usage: fatigue_crack_growth --demo CYCLES [--varying-geometry | --finite-width] (1..=1000000); synthetic data only".into(),
        );
    }
    let varying_geometry = args.get(2).is_some_and(|flag| flag == "--varying-geometry");
    let finite_width = args.get(2).is_some_and(|flag| flag == "--finite-width");
    let finite_plate = if finite_width {
        Some(FiniteWidthCentreCrack::new(FULL_PLATE_WIDTH_M)?)
    } else {
        None
    };
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
            surface_and_geometry: if finite_width {
                "centred through crack in a long elastic plate; full W=0.02 m; physical half-length a".into()
            } else if varying_geometry {
                "synthetic Y(a)=sqrt(a/a0) table; a is crack half-length".into()
            } else {
                "declared constant Y=1; a is crack half-length".into()
            },
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
    let crack_geometry = if let Some(plate) = finite_plate {
        let minimum = 0.5 * INITIAL_LENGTH_M;
        let maximum = plate.maximum_half_length_m();
        // The bounded helper retains these exact endpoints and chooses the
        // interior knots from its analytic relative-factor interpolation bound.
        CrackGeometry::Tabulated(Arc::new(plate.tabulate_to_tolerance(
            minimum,
            maximum,
            FINITE_WIDTH_RELATIVE_Y_TOLERANCE,
            FINITE_WIDTH_MAXIMUM_POINTS,
        )?))
    } else if varying_geometry {
        CrackGeometry::Tabulated(Arc::new(TabulatedCrackGeometry::new(
            CrackGeometryTableContext {
                source_reference: "SYNTHETIC power profile; not a calibrated component geometry"
                    .into(),
                crack_length_convention: "a is crack half-length in metres; a0=0.001 m".into(),
                nominal_stress_basis:
                    "uncracked uniform tensile stress on the fixed x-normal plane".into(),
            },
            [0.5_f64, 1.0, 1.005, 2.0, 10.0]
                .into_iter()
                .map(|ratio| CrackGeometryPoint {
                    crack_length_m: ratio * INITIAL_LENGTH_M,
                    geometry_factor: ratio.sqrt(),
                })
                .collect(),
            5,
        )?))
    } else {
        ConstantGeometryCrack {
            geometry_factor: 1.0,
            minimum_length_m: 1e-5,
            maximum_length_m: 0.01,
        }
        .into()
    };
    let gauge = TetCrackGrowthGauge {
        element: 0,
        normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: [104; 32],
        material_source: None,
        law,
        geometry: crack_geometry,
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
        // Independent m=2 integrals: constant Y gives da/dN=k*a, while
        // Y(a)=sqrt(a/a0) gives da/dN=k*a^2/a0.
        let normalized_range = (MAXIMUM_STRESS_PA - MINIMUM_STRESS_PA) / REFERENCE_DELTA_K;
        let exposure = REFERENCE_RATE * normalized_range.powi(2) * std::f64::consts::PI * cycles as f64;
        let expected = if finite_width {
            None
        } else if varying_geometry {
            Some(INITIAL_LENGTH_M / (1.0 - exposure))
        } else {
            Some(INITIAL_LENGTH_M * exposure.exp())
        };
        println!("SYNTHETIC demonstration, not validated component life");
        println!("Declared identical quasistatic cycles; uncracked mechanics stay fixed.");
        println!("geometry={}", if finite_width {
            "finite-width centred through crack; elastic secant factor in a log-log table"
        } else if varying_geometry {
            "declared log-log table Y(a)=sqrt(a/a0), evaluated as the crack grows"
        } else {
            "declared constant Y=1"
        });
        println!("native_samples={} repeated_full_cycles={cycles}", path.observed_samples());
        println!("normal_stress_min_pa={MINIMUM_STRESS_PA:.6e} normal_stress_max_pa={MAXIMUM_STRESS_PA:.6e}");
        println!("initial_crack_length_m={:.12e} final_crack_length_m={:.12e}",
            step.before.crack_length_m, step.after.crack_length_m);
        println!("initial_delta_k_pa_sqrt_m={:.12e} final_delta_k_pa_sqrt_m={:.12e}",
            step.initial_delta_k_pa_sqrt_m, step.final_delta_k_pa_sqrt_m);
        if let Some(expected) = expected {
            println!("analytic_m2_length_m={expected:.12e} relative_error={:.6e}",
                (step.after.crack_length_m - expected).abs() / expected);
        }
        if let Some(plate) = finite_plate {
            println!("remote_uniform_gross_section_tension; physical half-length a; full crack length=2a");
            println!("full_plate_width_m={:.6e} initial_full_crack_length_m={:.12e} final_full_crack_length_m={:.12e}",
                plate.plate_width_m(), 2.0 * step.before.crack_length_m, 2.0 * step.after.crack_length_m);
            println!("maximum_supported_half_length_m={:.6e} final_full_crack_to_width_ratio={:.6e} supported_ratio_max=0.8",
                plate.maximum_half_length_m(), 2.0 * step.after.crack_length_m / plate.plate_width_m());
            if let Some(CrackGeometry::Tabulated(table)) = path.geometry(0) {
                let (minimum, maximum) = table.length_bounds();
                println!("geometry_table_points={} selection=adaptive relative_y_tolerance={FINITE_WIDTH_RELATIVE_Y_TOLERANCE:.6e} maximum_points={FINITE_WIDTH_MAXIMUM_POINTS} half_length_min_m={minimum:.6e} half_length_max_m={maximum:.6e}",
                    table.points().len());
                println!("initial_tabulated_y={:.12e} final_tabulated_y={:.12e} final_direct_secant_y={:.12e}",
                    table.factor_at(step.before.crack_length_m)?,
                    table.factor_at(step.after.crack_length_m)?,
                    plate.factor_at(step.after.crack_length_m)?);
                // The secant factor increases across this supported interval.
                // Its table extrema bound da/dN=k*Y(a)^2*a, giving independent
                // exponential growth bounds, not an interpolation-error estimate.
                let minimum_y = table.points()[0].geometry_factor;
                let maximum_y = table.points()[table.points().len() - 1].geometry_factor;
                let lower = INITIAL_LENGTH_M * (exposure * minimum_y.powi(2)).exp();
                let upper = INITIAL_LENGTH_M * (exposure * maximum_y.powi(2)).exp();
                println!("bounded_m2_length_min_m={lower:.12e} bounded_m2_length_max_m={upper:.12e} within_growth_bounds={}",
                    (lower..=upper).contains(&step.after.crack_length_m));
                println!("crack_length_increased={} within_geometry_table={}",
                    step.after.crack_length_m > step.before.crack_length_m,
                    (minimum..=maximum).contains(&step.after.crack_length_m));
                println!("Requested tolerance controls analytic relative-Y interpolation; floating-point intervals and material life are not certified.");
            }
        }
        Ok(())
    })
}
