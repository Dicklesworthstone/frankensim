//! A finite heated block cooled by native D3Q19 flow and D3Q7 heat transport.
//!
//! Run: cargo run -p fs-lbm --example immersed_body_cooling -- 400 0.025 resolved
//! Use `lumped` as the last argument to impose an isothermal solid interior.
//!
//! All coordinates are lattice units. theta=(T-T_inlet)/delta_T, the fluid
//! cell heat capacity is one, and the block's TOTAL capacity is 40. The
//! transverse directions are periodic; this is a repeated-cell demonstration,
//! not an isolated-body correlation or a grid-converged heat-transfer claim.
//! The default resolves all 27 solid voxels with conductivity 0.2 and a
//! surface-contact resistance of 0.5. Only the interior core receives heat.

use std::error::Error;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::{BoundarySpec3, CollisionModel3, FaceBoundary3};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3};
use fs_lbm::thermal::boundary::{FiniteThermalWall, ThermalWall};
use fs_lbm::thermal::open::immersed::ImmersedThermalBody3;
use fs_lbm::thermal::open::immersed::voxel::ImmersedThermalVoxel3;
use fs_lbm::thermal::open::{ThermalOpenBoundaries3, ThermalOpenFace3, ThermalOpenLimits3};

fn simulate(cx: &Cx<'_>, steps: usize, speed: f64, resolved: bool) -> Result<(), Box<dyn Error>> {
    let dims = [24, 10, 10];
    let cells = dims.into_iter().product();
    let is_solid =
        |[x, y, z]: [usize; 3]| (8..11).contains(&x) && (4..7).contains(&y) && (4..7).contains(&z);
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
    let (body_count, heated_node) = if resolved {
        let core = 9 + dims[0] * 5 + dims[0] * dims[1] * 5;
        let voxels: Vec<_> = body_cells
            .iter()
            .map(|&cell| ImmersedThermalVoxel3 {
                cell,
                temperature: 1.0,
                heat_capacity: 40.0 / body_cells.len() as f64,
                conductivity: 0.2,
                contact_resistance: 0.5,
                heat_per_step: if cell == core { 0.02 } else { 0.0 },
            })
            .collect();
        let mapping = thermal.set_immersed_thermal_voxels(&voxels, Default::default())?;
        (
            mapping.cells().len(),
            mapping.node_for_cell(core).ok_or("missing heated core")?,
        )
    } else {
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
        (1, 0)
    };
    let heat = vec![0.0; cells];
    let mut worst_defect = 0.0_f64;
    println!("# lattice units; total block capacity=40; source=0.02 until halfway, then zero");
    println!("# resolved solid={resolved}; thermal nodes={body_count}");
    println!(
        "step,mean_body_theta,min_body_theta,max_body_theta,maximum_fluid_theta,body_heat_into_fluid,combined_defect"
    );
    for step in 1..=steps {
        if step == steps / 2 + 1 {
            thermal.set_immersed_body_heat(heated_node, 0.0)?;
        }
        let receipt =
            flow.step_serial_with_thermal(&mut thermal, cx, None, &heat, None, [0.0; 3])?;
        worst_defect = worst_defect.max(receipt.heat.combined_defect.abs());
        if step == 1 || step == steps || step % 20 == 0 {
            let mut body_heat = 0.0;
            let mut mean = 0.0;
            let mut minimum = f64::INFINITY;
            let mut maximum_body = f64::NEG_INFINITY;
            for node in 0..body_count {
                let body = thermal.immersed_body(node).ok_or("missing solid node")?;
                body_heat += body
                    .last_exchange
                    .ok_or("missing solid heat balance")?
                    .heat_in;
                mean += body.solid.temperature * (body.solid.heat_capacity / 40.0);
                minimum = minimum.min(body.solid.temperature);
                maximum_body = maximum_body.max(body.solid.temperature);
            }
            let maximum = (0..cells)
                .filter_map(|cell| thermal.temperature(cell))
                .fold(f64::NEG_INFINITY, f64::max);
            println!(
                "{step},{mean:.9},{minimum:.9},{maximum_body:.9},{maximum:.9},{body_heat:.9e},{:.3e}",
                receipt.heat.combined_defect
            );
        }
    }
    println!("# maximum absolute combined defect={worst_defect:.3e}");
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    if arguments.len() > 3 {
        return Err(
            "usage: immersed_body_cooling [STEPS=400] [INLET_SPEED=0.025] [resolved|lumped]".into(),
        );
    }
    let steps: usize = arguments.first().map_or(Ok(400), |s| s.parse())?;
    let speed: f64 = arguments.get(1).map_or(Ok(0.025), |s| s.parse())?;
    let resolved = match arguments.get(2).map_or("resolved", String::as_str) {
        "resolved" => true,
        "lumped" => false,
        _ => return Err("solid model must be resolved or lumped".into()),
    };
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
        simulate(&cx, steps, speed, resolved)
    })
}
