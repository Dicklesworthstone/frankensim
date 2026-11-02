//! Source-scaled fluid preparation, held loading and physical momentum transfer.
use super::*;
use fs_couple::material_immersed::dry::{ResolvedImmersedBodies, ResolvedImmersedBody};
use fs_couple::material_immersed::settling::SubmergedSphere;
use fs_lbm::d3q19::{BoundaryGrid3, DryBodyCouplingConfig3};

fn fields(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let [nx, ny, nz] = grid.dimensions();
    let mut values = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                values.push(grid.populations(x, y, z));
            }
        }
    }
    values
}
fn physical_momentum_si(model: &RetainedFluidSpheres<'_, ResolvedLiquidState, 2>) -> [f64; 3] {
    let mut total = [0.0; 3];
    // Independent dimensional conversion: cell mass rho dx^3 times dx/dt.
    let unit = model.scaling().transport().density_kg_m3() * DX.powi(4) / DT;
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let f = model.grid().populations(x, y, z);
                let force = model.grid().force_density(x, y, z);
                let mut p = force.map(|f| 0.5 * f);
                for (value, (ex, ey, ez)) in f.into_iter().zip(E3) {
                    for (a, e) in [ex, ey, ez].into_iter().enumerate() {
                        p[a] += value * f64::from(e);
                    }
                }
                for a in 0..3 {
                    total[a] += unit * p[a];
                }
            }
        }
    }
    for i in 0..2 {
        let state = model.state_si(i).unwrap();
        for a in 0..3 {
            total[a] += model.sources()[i].material().inertia().mass() * state.velocity_m_s[a];
        }
    }
    total
}

#[test]
fn sourced_fluid_load_closes_si_momentum_without_rebind_impulse_or_reset() {
    let density = 3.0;
    let fluid = fluid(density, 3e-4);
    let warmer = super::fluid(density, 6e-4);
    let source = body(6.0);
    let mut model = assembly(&fluid, &source);
    let gate = CancelGate::new();
    let velocity = [0.0002, -0.0001, 0.00015];
    model.initialize_uniform_flow_si(velocity, &gate).unwrap();
    let initial_bodies = [model.state_si(0).unwrap(), model.state_si(1).unwrap()];
    for a in 0..3 {
        close(
            model.grid().velocity(2, 3, 4)[a] * DX / DT,
            velocity[a],
            2e-17,
        );
    }
    let load = [1e-4, -2e-4, 1.5e-4];
    let before = physical_momentum_si(&model);
    model.rebind_uniform_fluid_force_si(load, &gate).unwrap();
    for a in 0..3 {
        close(physical_momentum_si(&model)[a], before[a], 2e-19);
        close(
            model.grid().force_density(2, 3, 4)[a],
            load[a] * DT * DT / (density * DX),
            2e-18,
        );
        close(
            model.uniform_fluid_force_density_n_m3().unwrap()[a],
            load[a],
            2e-18,
        );
    }
    assert_eq!(
        [model.state_si(0).unwrap(), model.state_si(1).unwrap()],
        initial_bodies
    );
    assert_eq!(model.accepted_steps(), 0);
    let accepted = snapshot(&model);
    model.rebind_uniform_fluid_force_si(load, &gate).unwrap();
    assert_eq!(snapshot(&model), accepted);
    model
        .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
        .unwrap();
    let after = physical_momentum_si(&model);
    for a in 0..3 {
        close(
            after[a] - before[a],
            load[a] * (8.0 * DX).powi(3) * DT,
            2e-16,
        );
    }
    assert_eq!(model.accepted_steps(), 1);
    let accepted = snapshot(&model);
    assert!(model.initialize_uniform_flow_si([0.0; 3], &gate).is_err());
    assert_eq!(snapshot(&model), accepted);
    model.rebind_fluid(&warmer, &gate).unwrap();
    assert_eq!(snapshot(&model), accepted);
    for a in 0..3 {
        close(
            model.uniform_fluid_force_density_n_m3().unwrap()[a],
            load[a],
            2e-18,
        );
    }
    let before_switch = physical_momentum_si(&model);
    model
        .rebind_uniform_fluid_force_si([0.0; 3], &gate)
        .unwrap();
    for a in 0..3 {
        close(physical_momentum_si(&model)[a], before_switch[a], 2e-19);
    }
    assert_eq!(model.accepted_steps(), 1);
}

#[test]
fn dry_single_assembly_and_plane_wrappers_retain_initial_body_and_shared_clock() {
    let fluid = fluid(1.0, 1e-4);
    let source = body(1.0);
    let gate = CancelGate::new();
    let initial = SphereStateSi {
        center_m: [0.04; 3],
        velocity_m_s: [0.0; 3],
        ..initial()[0]
    };
    let mut single =
        ResolvedImmersedBody::sphere(&fluid, domain(), &source, RADIUS, initial, [0.0; 3], &gate)
            .unwrap();
    let load = [1e-5, -2e-5, 3e-5];
    let flow = [0.0001, -0.0002, 0.0001];
    single.rebind_uniform_fluid_force_si(load, &gate).unwrap();
    single.initialize_uniform_flow_si(flow, &gate).unwrap();
    assert_eq!(single.state_si().unwrap(), initial);
    for a in 0..3 {
        close(single.grid().velocity(3, 4, 2)[a] * DX / DT, flow[a], 2e-17);
        close(
            single.uniform_fluid_force_density_n_m3().unwrap()[a],
            load[a],
            2e-18,
        );
    }
    let mut assembly = ResolvedImmersedBodies::spheres(
        &fluid,
        domain(),
        [SubmergedSphere::new(&source, RADIUS).unwrap()],
        [initial],
        [0.0; 3],
        film(),
        &gate,
    )
    .unwrap();
    assembly.initialize_uniform_flow_si(flow, &gate).unwrap();
    assembly.rebind_uniform_fluid_force_si(load, &gate).unwrap();
    let before = fields(assembly.grid());
    let mut planes = assembly.with_prescribed_planes(&[], &gate).unwrap();
    planes.rebind_uniform_fluid_force_si(load, &gate).unwrap();
    assert_eq!(fields(planes.assembly().grid()), before);
    assert_eq!(planes.assembly().state_si(0).unwrap(), initial);
    for a in 0..3 {
        close(
            planes.uniform_fluid_force_density_n_m3().unwrap()[a],
            load[a],
            2e-18,
        );
    }
    planes
        .step_forces_si(
            [[0.0; 6]],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(planes.assembly().accepted_steps(), 1);
    let accepted = fields(planes.assembly().grid());
    assert!(planes.initialize_uniform_flow_si(flow, &gate).is_err());
    assert_eq!(fields(planes.assembly().grid()), accepted);
}

#[test]
fn invalid_speed_force_and_cancel_preserve_material_flow_and_body_history() {
    let fluid = fluid(1.0, 1e-4);
    let source = body(2.0);
    let mut model = assembly(&fluid, &source);
    let before = snapshot(&model);
    let bodies = [model.state_si(0).unwrap(), model.state_si(1).unwrap()];
    let gate = CancelGate::new();
    for velocity in [[0.002, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
        assert!(model.initialize_uniform_flow_si(velocity, &gate).is_err());
        assert_eq!(snapshot(&model), before);
    }
    assert!(
        model
            .rebind_uniform_fluid_force_si([f64::INFINITY, 0.0, 0.0], &gate)
            .is_err()
    );
    assert_eq!(snapshot(&model), before);
    gate.request();
    assert!(
        model
            .initialize_uniform_flow_si([0.0001; 3], &gate)
            .is_err()
    );
    assert!(
        model
            .rebind_uniform_fluid_force_si([1e-5; 3], &gate)
            .is_err()
    );
    assert_eq!(snapshot(&model), before);
    assert_eq!(
        [model.state_si(0).unwrap(), model.state_si(1).unwrap()],
        bodies
    );
    assert_eq!(model.uniform_fluid_force_density_n_m3().unwrap(), [0.0; 3]);
    assert_eq!(model.accepted_steps(), 0);
}
