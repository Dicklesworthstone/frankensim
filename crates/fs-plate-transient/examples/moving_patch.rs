//! Moving pressure footprint with bending couples and load-history adjoints.
//!
//! cargo run -p fs-plate-transient --example moving_patch
//!
//! Uses the original assembled plate, integrator, observations and discrete
//! adjoint. All geometry/material values are illustrative SI inputs; no contact
//! or experimentally validated damping model is inferred.
#[path = "moving_plate.rs"]
#[allow(dead_code)]
mod moving;

use fs_plate::loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint};
use fs_plate_transient::{MovingPlateLoad, MovingPlateLoadBudget, PlateLoadKeyframe};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let (mesh, plate) = moving::build_plate()?;
    let footprint = [
        [-0.025, -0.02],
        [0.025, -0.02],
        [-0.025, 0.02],
        [0.025, 0.02],
    ]
    .map(|point| PlatePatchPoint {
        point,
        area_m2: 0.0005,
    });
    let knots = [
        PlateLoadKeyframe {
            time_s: 0.0,
            center_m: [0.12, 0.27],
            load: [1.0, 0.002, -0.003],
        },
        PlateLoadKeyframe {
            time_s: moving::STEP * 6.0,
            center_m: [0.5, 0.32],
            load: [1.7, -0.004, 0.002],
        },
        PlateLoadKeyframe {
            time_s: moving::STEP * moving::STEPS as f64,
            center_m: [0.88, 0.27],
            load: [0.6, 0.003, 0.001],
        },
    ];
    let load = MovingPlateLoad::new(
        &mesh,
        &plate,
        &footprint,
        &knots,
        MovingPlateLoadBudget {
            max_keyframes: knots.len(),
            patch: PlatePatchBudget {
                max_points: footprint.len(),
                location: PlateLoadBudget {
                    max_nodes: mesh.nodes.len(),
                    max_triangles: mesh.tris.len(),
                },
            },
        },
    )?;
    let run = moving::evaluate_load(
        &mesh,
        &plate,
        moving::dynamics_parameters([1.1, 0.9, 3.0, 2e-4, 0.0]),
        &load,
    )?;
    println!("time_s,sensor_deflection_m");
    for (time, deflection) in &run.samples {
        println!("{time:.6},{deflection:.12e}");
    }
    println!(
        "sampled_loss={:.12e}; checkpoints={}; forward_replays={}",
        run.loss, run.checkpoints, run.replays
    );
    for (knot, gradient) in run.gradient[4..].chunks_exact(3).enumerate() {
        println!(
            "knot={knot}; d_loss/d_Fz={:.12e}; d_loss/d_Mx={:.12e}; d_loss/d_My={:.12e}",
            gradient[0], gradient[1], gradient[2]
        );
    }
    Ok(())
}
