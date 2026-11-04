//! Two-way impact of a finite-mass body on a spatial clamped DKT plate.
//!
//! Run: cargo run -p fs-plate-transient --example plate_impact
//! Standard output is CSV; the short run summary goes to standard error.
//!
//! All geometry, elastic properties and contact coefficients below are explicit
//! illustrative SI inputs, not a material-calibrated impact model. The finite
//! footprint distributes the resultant force and observes plate motion through
//! the same area-weighted stencil. Its area does not multiply the already-total
//! contact stiffness. Positive body travel closes the gap; positive plate travel
//! opens it. The body mass, every free plate DOF, and the existing Hunt--Crossley
//! loss participate in the same implicit generalized-alpha step.
//!
//! Reported power is instantaneous physical dissipation. Reported energy is
//! physical kinetic plus plate/contact potential energy. Generalized-alpha also
//! has numerical dissipation; endpoint powers are not an exact step energy budget.

use fs_dcontact::Obstacle;
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateMesh, PlateSection, assemble,
    loading::{PlateLoadBudget, PlatePatchBudget, PlatePatchPoint, PlatePatchStencil},
};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters, PlateImpact,
    PlateImpactBudget, PlateImpactProbe,
};
use fs_solver::NewtonKrylovConfig;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderState,
    initialization::{InitialSolveConfig, initial_workspace_components, second_order_acceleration},
    second_order_adjoint::trajectory::StructuralTrajectoryModel,
};
use std::error::Error;

const STEP_S: f64 = 5.0e-5;
const STEPS: usize = 200;
const RHO_INF: f64 = 0.8;

fn write_sample(time_s: f64, probe: &PlateImpactProbe) {
    println!(
        "{time_s:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
        probe.indenter_displacement_m,
        probe.indenter_velocity_m_per_s,
        probe.patch.displacement,
        probe.patch.velocity,
        probe.maximum_penetration_m,
        probe.normal_force_n,
        probe.total_energy_j,
        probe.contact_dissipated_power_w,
        probe.plate_dissipated_power_w,
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    let mesh = PlateMesh::rectangle(1.0, 0.6, 4, 4);
    let section = PlateSection::isotropic(3.0e9, 0.28, 0.008, 600.0)?;
    let plate = assemble(
        &mesh,
        &section,
        &PlateMesh::rectangle_boundary(4, 4),
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
    )?;
    // Four equal-area samples of a 5 cm by 5 cm fixed footprint.
    let points = [
        PlatePatchPoint {
            point: [0.4875, 0.2875],
            area_m2: 0.000625,
        },
        PlatePatchPoint {
            point: [0.5125, 0.2875],
            area_m2: 0.000625,
        },
        PlatePatchPoint {
            point: [0.4875, 0.3125],
            area_m2: 0.000625,
        },
        PlatePatchPoint {
            point: [0.5125, 0.3125],
            area_m2: 0.000625,
        },
    ];
    let patch = PlatePatchStencil::locate(
        &mesh,
        &plate,
        &points,
        PlatePatchBudget {
            max_points: points.len(),
            location: PlateLoadBudget {
                max_nodes: mesh.nodes.len(),
                max_triangles: mesh.tris.len(),
            },
        },
    )?;
    let no_load = NoPlateLoad;
    let dynamics = PlateDynamics::new(
        &plate,
        PlateDynamicsParameters {
            stiffness_scale: 1.0,
            mass_scale: 1.0,
            mass_damping_per_s: 0.1,
            stiffness_damping_s: 1.0e-5,
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
    // One total normal contact law, with a 50 micrometre initial clearance.
    // K is in N / m^(3/2); chi is in s/m. The original law supplies unilateral
    // activation and the non-attractive unloading cap.
    let law = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![5.0e-5],
        vec![1.0],
        2.0e6,
        1.5,
        "illustrative total Hertz-type footprint; not material calibrated".into(),
    )?
    .with_internal_loss(1.0)?;
    let impact = PlateImpact::new(
        dynamics,
        patch,
        0.02,
        law,
        PlateImpactBudget {
            max_free_dofs: plate.free + 1,
            max_contact_points: 1,
        },
    )?;
    let n = impact.dimension();
    let body = impact.plate_dofs();
    let q = vec![0.0; n];
    let mut v = vec![0.0; n];
    v[body] = 0.2; // 20 g body approaches the plate at 0.2 m/s.
    let mut forcing = vec![0.0; n];
    impact.forcing(0.0, &mut forcing)?;
    let initial_config = InitialSolveConfig {
        restart: 2,
        max_cycles: 3,
        tolerance: 1.0e-12,
    };
    let acceleration = second_order_acceleration(
        &impact,
        &q,
        &v,
        &forcing,
        &impact.mass_preconditioner(),
        initial_config,
        initial_workspace_components(n, 0, initial_config).ok_or("initial workspace overflow")?,
        &mut || false,
    )?;
    let mut state = SecondOrderState::new(0.0, &q, &v, &acceleration.value);
    let method = OperatorGeneralizedAlpha::new(
        n,
        STEP_S,
        RHO_INF,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1.0e-10,
                relative_tolerance: 1.0e-10,
                linear_restart: n,
                max_linear_cycles: 4,
                forcing_minimum: 1.0e-12,
                forcing_maximum: 1.0e-10,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 16,
        },
    );
    let initial = impact.probe(&state.q, &state.v)?;
    let mut maximum_penetration: f64 = 0.0;
    let mut peak_force: f64 = 0.0;
    let mut rebound_observed = false;
    println!(
        "time_s,body_displacement_m,body_velocity_m_per_s,plate_displacement_m,plate_velocity_m_per_s,penetration_m,normal_force_n,total_energy_j,contact_dissipated_power_w,plate_dissipated_power_w"
    );
    write_sample(state.t, &initial);
    for _ in 0..STEPS {
        impact.forcing(method.forcing_time(state.t)?, &mut forcing)?;
        method.step(&mut state, &impact, &forcing)?;
        let probe = impact.probe(&state.q, &state.v)?;
        maximum_penetration = maximum_penetration.max(probe.maximum_penetration_m);
        peak_force = peak_force.max(probe.normal_force_n);
        rebound_observed |= probe.indenter_velocity_m_per_s < 0.0;
        write_sample(state.t, &probe);
    }
    let final_probe = impact.probe(&state.q, &state.v)?;
    eprintln!(
        "plate_dofs={}, steps={}, initial_energy_j={:.8e}, final_energy_j={:.8e}, maximum_penetration_m={maximum_penetration:.8e}, peak_force_n={peak_force:.8e}, rebound_observed={rebound_observed}",
        body, state.steps, initial.total_energy_j, final_probe.total_energy_j,
    );
    Ok(())
}
