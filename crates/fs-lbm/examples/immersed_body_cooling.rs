//! A finite heated block cooled by native D3Q19 flow and D3Q7 heat transport.
//!
//! Run: cargo run -p fs-lbm --example immersed_body_cooling -- 400 0.025
//!
//! All coordinates are lattice units. theta=(T-T_inlet)/delta_T, the fluid
//! cell heat capacity is one, and the block's TOTAL capacity is 40. The
//! transverse directions are periodic; this is a repeated-cell demonstration,
//! not an isolated-body correlation or a grid-converged heat-transfer claim.
//! The block is internally isothermal with a declared surface resistance.

use std::error::Error;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::{BoundarySpec3, CollisionModel3, FaceBoundary3};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3};
use fs_lbm::thermal::boundary::{FiniteThermalWall, ThermalWall};
use fs_lbm::thermal::open::immersed::ImmersedThermalBody3;
use fs_lbm::thermal::open::{ThermalOpenBoundaries3, ThermalOpenFace3, ThermalOpenLimits3};

fn simulate(cx: &Cx<'_>, steps: usize, speed: f64) -> Result<(), Box<dyn Error>> {
    let dims = [24, 10, 10];
    let cells = dims.into_iter().product();
    let is_solid =
        |[x, y, z]: [usize; 3]| (8..10).contains(&x) && (4..6).contains(&y) && (4..6).contains(&z);
    let grid = Grid3::<D3q19, 19>::try_new_equilibrium(
        dims,
        CollisionModel3::Bgk { tau: 0.9 },
        [0.0; 3],
        [false, true, true],
        GridBudget3 {
            maximum_cells: cells,
            maximum_mach: 0.2,
        },
        is_solid,
        |_| (1.0, [speed, 0.0, 0.0]),
    )?;
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[0] = FaceBoundary3::Velocity {
        velocity: [speed, 0.0, 0.0],
    };
    faces[1] = FaceBoundary3::Pressure { density: 1.0 };
    let mut flow = OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces))?;
    let mut thermal_faces = [ThermalOpenFace3::Periodic; 6];
    thermal_faces[0] = ThermalOpenFace3::Reservoir { temperature: 0.0 };
    // Explicit outlet temperature also defines any local backflow.
    thermal_faces[1] = ThermalOpenFace3::Reservoir { temperature: 0.0 };
    let mut thermal = flow.prepare_thermal_field(
        cx,
        &vec![0.0; cells],
        0.85,
        ThermalOpenBoundaries3 {
            faces: thermal_faces,
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: cells,
            maximum_storage_bytes: 4 * 1024 * 1024,
            temperature_range: [-0.1, 2.0],
            balance_tolerance: 1e-10,
        },
    )?;
    let body_cells: Vec<_> = (0..cells)
        .filter(|&cell| {
            is_solid([
                cell % dims[0],
                cell / dims[0] % dims[1],
                cell / (dims[0] * dims[1]),
            ])
        })
        .collect();
    thermal.set_immersed_thermal_bodies(
        &[ImmersedThermalBody3 {
            solid_cells: &body_cells,
            solid: FiniteThermalWall {
                temperature: 1.0,
                heat_capacity: 40.0,
                resistance: 0.5,
                heat_per_step: 0.02,
            },
        }],
        128 * 1024,
    )?;
    let heat = vec![0.0; cells];
    let mut worst_defect = 0.0_f64;
    println!("# lattice units; total block capacity=40; source=0.02 until halfway, then zero");
    println!("step,body_theta,maximum_fluid_theta,body_heat_into_fluid,combined_defect");
    for step in 1..=steps {
        if step == steps / 2 + 1 {
            thermal.set_immersed_body_heat(0, 0.0)?;
        }
        let receipt =
            flow.step_serial_with_thermal(&mut thermal, cx, None, &heat, None, [0.0; 3])?;
        worst_defect = worst_defect.max(receipt.heat.combined_defect.abs());
        if step == 1 || step == steps || step % 20 == 0 {
            let body = thermal.immersed_body(0).ok_or("missing finite block")?;
            let exchange = body.last_exchange.ok_or("missing block heat balance")?;
            let maximum = (0..cells)
                .filter_map(|cell| thermal.temperature(cell))
                .fold(f64::NEG_INFINITY, f64::max);
            println!(
                "{step},{:.9},{maximum:.9},{:.9e},{:.3e}",
                body.solid.temperature, exchange.heat_in, receipt.heat.combined_defect
            );
        }
    }
    println!("# maximum absolute combined defect={worst_defect:.3e}");
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    if arguments.len() > 2 {
        return Err("usage: immersed_body_cooling [STEPS=400] [INLET_SPEED=0.025]".into());
    }
    let steps: usize = arguments.first().map_or(Ok(400), |s| s.parse())?;
    let speed: f64 = arguments.get(1).map_or(Ok(0.025), |s| s.parse())?;
    if !(2..=10_000).contains(&steps) || !speed.is_finite() || !(0.0..=0.05).contains(&speed) {
        return Err("need 2..10000 steps and an inlet speed in [0,0.05] lattice units".into());
    }
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 228,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        simulate(&cx, steps, speed)
    })
}
