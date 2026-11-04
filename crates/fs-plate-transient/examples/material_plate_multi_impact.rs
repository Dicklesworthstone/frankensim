//! Two material-backed bodies interacting through one spatial DKT plate.
//!
//! `cargo run -p fs-plate-transient --example material_plate_multi_impact`
//!
//! The first two scenarios change only body 0's initial velocity. Body 1 has the
//! same initial motion, source, mass and footprint in both; its changed response
//! measures coupling through the shared plate. The third adds smooth 5 ms body
//! force pulses. Applied work uses the integrator's stage-power quadrature.
//! All cards are explicit synthetic
//! fixtures, not experimental evidence. The actual 8 mm plate thickness also
//! supplies both sources' half-space depth and layer thickness. Fixed footprints
//! transfer total forces and sample plate motion, not changing Hertz pressure
//! distributions. Energy is physical storage, with no exact step-balance claim.

#[path = "material_plate_impact.rs"]
#[allow(dead_code)]
mod single;

use fs_plate::{
    PlateMesh, PlateModel,
    loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint, PlatePatchStencil},
};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateImpact,
    PlateImpactBudget, PlateImpactSet, PlateImpactSetBudget, PlateLoad,
};
use fs_time::galpha::{
    SecondOrderState,
    initialization::{InitialSolveConfig, initial_workspace_components, second_order_acceleration},
    second_order_adjoint::trajectory::StructuralTrajectoryModel,
};
use std::error::Error;

struct BodyDrive {
    amplitudes_n: [f64; 2],
}
impl BodyDrive {
    fn profile(time: f64) -> Result<f64, String> {
        if !time.is_finite() || time < 0.0 {
            return Err("body-drive time must be finite and nonnegative".into());
        }
        let s = time / (100.0 * single::STEP_S);
        Ok(if s <= 1.0 {
            16.0 * s * s * (1.0 - s) * (1.0 - s)
        } else {
            0.0
        })
    }
}
impl PlateLoad for BodyDrive {
    fn parameter_count(&self) -> usize {
        2
    }
    fn forcing(&self, time: f64, out: &mut [f64]) -> Result<(), String> {
        if out.len() != 2 || self.amplitudes_n.iter().any(|x| !x.is_finite()) {
            return Err("two finite body-drive amplitudes required".into());
        }
        let shape = Self::profile(time)?;
        for (f, amplitude) in out.iter_mut().zip(self.amplitudes_n) {
            *f = shape * amplitude;
        }
        Ok(())
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        if seed.len() != 2 || out.len() != 2 || seed.iter().any(|x| !x.is_finite()) {
            return Err("two finite body-drive seeds required".into());
        }
        let shape = Self::profile(time)?;
        for (bar, seed) in out.iter_mut().zip(seed) {
            *bar = shape * seed;
        }
        Ok(())
    }
}

fn patch(
    mesh: &PlateMesh,
    plate: &PlateModel,
    center_x: f64,
) -> Result<PlatePatchStencil, Box<dyn Error>> {
    let points = [
        [-1.5e-4, -1.5e-4],
        [1.5e-4, -1.5e-4],
        [-1.5e-4, 1.5e-4],
        [1.5e-4, 1.5e-4],
    ]
    .map(|offset| PlatePatchPoint {
        point: [center_x + offset[0], 0.15 + offset[1]],
        area_m2: 9.0e-8,
    });
    Ok(PlatePatchStencil::locate(
        mesh,
        plate,
        &points,
        PlatePatchBudget {
            max_points: points.len(),
            location: PlateLoadBudget {
                max_nodes: mesh.nodes.len(),
                max_triangles: mesh.tris.len(),
            },
        },
    )?)
}

struct Scenario {
    source_ids: [String; 2],
    plate_state_id: String,
    peak_force_n: [f64; 2],
    peak_plate_displacement_m: [f64; 2],
    final_body_displacement_m: [f64; 2],
    final_body_velocity_m_per_s: [f64; 2],
    initial_energy_j: f64,
    final_energy_j: f64,
    drive_amplitudes_n: [f64; 2],
    sampled_body_drive_work_j: f64,
}

fn run(first_speed: f64, drive: Option<&BodyDrive>) -> Result<Scenario, Box<dyn Error>> {
    let first_source = single::synthetic_source("synthetic shared plate", 1.0e7, 900.0, 0.4)?;
    let second_source = single::synthetic_source("synthetic shared plate", 1.0e7, 900.0, 0.8)?;
    let first_bulk = first_source
        .model()
        .normal_state()
        .elastic()
        .surface_b_state_identity();
    let second_bulk = second_source
        .model()
        .normal_state()
        .elastic()
        .surface_b_state_identity();
    if first_bulk != second_bulk {
        return Err("both contacts must retain the same plate material state".into());
    }
    let source_ids = [
        first_source.model().identity().to_string(),
        second_source.model().identity().to_string(),
    ];
    let (mesh, plate) = single::build_plate(&first_source)?;
    let no_load = NoPlateLoad;
    // One physical plate and one applied-load owner. Copies borrow these exact
    // same owners; the set assembles their shared inertia/stiffness only once.
    let dynamics = PlateDynamics::new(
        &plate,
        PlateDynamicsParameters {
            stiffness_scale: 1.0,
            mass_scale: 1.0,
            mass_damping_per_s: 0.0,
            stiffness_damping_s: 0.0,
        },
        &no_load,
        PlateDynamicsBudget {
            max_dofs: plate.free,
            max_nonzeros: plate.k.nnz() + plate.m.nnz(),
            max_full_dofs: plate.dof_map.len(),
            max_load_parameters: 0,
        },
        &mut || false,
    )?;
    let body_budget = PlateImpactBudget {
        max_free_dofs: plate.free + 1,
        max_contact_points: 1,
    };
    let first = PlateImpact::from_source(
        dynamics,
        patch(&mesh, &plate, 0.17)?,
        0.001,
        first_source,
        single::GAP_M,
        body_budget,
    )?;
    let second = PlateImpact::from_source(
        dynamics,
        patch(&mesh, &plate, 0.23)?,
        0.0015,
        second_source,
        single::GAP_M,
        body_budget,
    )?;
    let impact = PlateImpactSet::new(
        vec![first, second],
        PlateImpactSetBudget {
            max_bodies: 2,
            max_free_dofs: plate.free + 2,
            max_contact_points: 2,
            max_parameters: 6,
        },
    )?;
    let impact = if let Some(drive) = drive {
        impact.with_body_forcing(drive)?
    } else {
        impact
    };
    let n = impact.dimension();
    let q = vec![0.0; n];
    let mut v = vec![0.0; n];
    v[plate.free] = first_speed;
    v[plate.free + 1] = 0.005;
    let mut forcing = vec![0.0; n];
    impact.forcing(0.0, &mut forcing)?;
    let config = InitialSolveConfig {
        restart: 2,
        max_cycles: 3,
        tolerance: 1.0e-12,
    };
    let a = second_order_acceleration(
        &impact,
        &q,
        &v,
        &forcing,
        &impact.mass_preconditioner(),
        config,
        initial_workspace_components(n, 0, config).ok_or("initial workspace overflow")?,
        &mut || false,
    )?;
    let mut state = SecondOrderState::new(0.0, &q, &v, &a.value);
    let integrator = single::method(n, single::STEP_S);
    let initial_energy_j = impact.probe(&state.q, &state.v)?.total_energy_j;
    let mut peak_force_n = [0.0_f64; 2];
    let mut peak_plate_displacement_m = [0.0_f64; 2];
    let mut sampled_body_drive_work_j = 0.0;
    for _ in 0..100 {
        let previous_velocity = [state.v[plate.free], state.v[plate.free + 1]];
        impact.forcing(integrator.forcing_time(state.t)?, &mut forcing)?;
        integrator.step(&mut state, &impact, &forcing)?;
        // single::method uses rho=1, so its actual stage velocity is the
        // endpoint average. This measures applied work, not total energy error.
        for i in 0..2 {
            sampled_body_drive_work_j += single::STEP_S
                * forcing[plate.free + i]
                * 0.5
                * (previous_velocity[i] + state.v[plate.free + i]);
        }
        let sample = impact.probe(&state.q, &state.v)?;
        for (i, body) in sample.bodies.iter().enumerate() {
            peak_force_n[i] = peak_force_n[i].max(body.normal_force_n);
            peak_plate_displacement_m[i] =
                peak_plate_displacement_m[i].max(body.patch.displacement.abs());
        }
    }
    let final_sample = impact.probe(&state.q, &state.v)?;
    Ok(Scenario {
        source_ids,
        plate_state_id: first_bulk.to_string(),
        peak_force_n,
        peak_plate_displacement_m,
        final_body_displacement_m: std::array::from_fn(|i| {
            final_sample.bodies[i].indenter_displacement_m
        }),
        final_body_velocity_m_per_s: std::array::from_fn(|i| {
            final_sample.bodies[i].indenter_velocity_m_per_s
        }),
        initial_energy_j,
        final_energy_j: final_sample.total_energy_j,
        drive_amplitudes_n: drive.map_or([0.0; 2], |input| input.amplitudes_n),
        sampled_body_drive_work_j,
    })
}

fn write_summary(name: &str, first_speed: f64, result: &Scenario) {
    println!(
        "{name},{},{},{},{first_speed:.8e},5.00000000e-3,{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
        result.source_ids[0],
        result.source_ids[1],
        result.plate_state_id,
        result.drive_amplitudes_n[0],
        result.drive_amplitudes_n[1],
        result.peak_force_n[0],
        result.peak_force_n[1],
        result.peak_plate_displacement_m[0],
        result.peak_plate_displacement_m[1],
        result.final_body_displacement_m[0],
        result.final_body_displacement_m[1],
        result.final_body_velocity_m_per_s[0],
        result.final_body_velocity_m_per_s[1],
        result.initial_energy_j,
        result.final_energy_j,
        result.sampled_body_drive_work_j
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    let baseline = run(0.0, None)?;
    let simultaneous = run(0.008, None)?;
    let drive = BodyDrive {
        amplitudes_n: [1.0e-3, 0.5e-3],
    };
    let controlled = run(0.008, Some(&drive))?;
    println!(
        "scenario,source_0_id,source_1_id,plate_state_id,initial_body_0_velocity_m_per_s,initial_body_1_velocity_m_per_s,drive_0_amplitude_n,drive_1_amplitude_n,peak_force_0_n,peak_force_1_n,peak_plate_displacement_0_m,peak_plate_displacement_1_m,final_body_0_displacement_m,final_body_1_displacement_m,final_body_0_velocity_m_per_s,final_body_1_velocity_m_per_s,initial_total_energy_j,final_total_energy_j,sampled_body_drive_work_j"
    );
    write_summary("first_body_stationary", 0.0, &baseline);
    write_summary("both_bodies_moving", 0.008, &simultaneous);
    write_summary("both_bodies_driven", 0.008, &controlled);
    eprintln!(
        "body_1_final_velocity_change_m_per_s={:.12e}; body_1_peak_force_change_n={:.12e}",
        simultaneous.final_body_velocity_m_per_s[1] - baseline.final_body_velocity_m_per_s[1],
        simultaneous.peak_force_n[1] - baseline.peak_force_n[1]
    );
    eprintln!(
        "drive_body_1_final_velocity_change_m_per_s={:.12e}; sampled_body_drive_work_j={:.12e}",
        controlled.final_body_velocity_m_per_s[1] - simultaneous.final_body_velocity_m_per_s[1],
        controlled.sampled_body_drive_work_j
    );
    Ok(())
}
