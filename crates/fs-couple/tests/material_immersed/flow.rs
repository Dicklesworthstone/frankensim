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

fn spatial_loads(scale: f64) -> Vec<[f64; 3]> {
    (0..512)
        .map(|index| {
            let x = index % 8;
            let y = (index / 8) % 8;
            let z = index / 64;
            [
                scale * 1e-5 * (1 + x + 2 * z) as f64,
                scale * -2e-5 * (1 + y) as f64,
                scale * 3e-5 * (1 + z + y) as f64,
            ]
        })
        .collect()
}

// Independent SI volume-force and fixed-origin torque integration, with no
// material conversion helper or native forcing receipt used by this oracle.
fn imposed_impulse_si(field: &[[f64; 3]], uniform: [f64; 3]) -> [f64; 6] {
    let mut total = [0.0; 6];
    for (index, force) in field.iter().enumerate() {
        let arm = [index % 8, (index / 8) % 8, index / 64].map(|i| (i as f64 + 0.5) * DX);
        let p: [f64; 3] = core::array::from_fn(|a| (force[a] + uniform[a]) * DX.powi(3) * DT);
        let l = [
            arm[1] * p[2] - arm[2] * p[1],
            arm[2] * p[0] - arm[0] * p[2],
            arm[0] * p[1] - arm[1] * p[0],
        ];
        for a in 0..3 {
            total[a] += p[a];
            total[a + 3] += l[a];
        }
    }
    total
}

// Full physical second moments distinguish a conservative force reframe from
// an equilibrium reset after viscous stresses have developed near the bodies.
fn physical_moments(grid: &BoundaryGrid3) -> Vec<(f64, [f64; 3], [[f64; 3]; 3])> {
    let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    (0..512)
        .map(|index| {
            let [x, y, z] = [index % 8, (index / 8) % 8, index / 64];
            let f = grid.populations(x, y, z);
            let force = grid.force_density(x, y, z);
            let rho = f.iter().sum::<f64>();
            let p: [f64; 3] = core::array::from_fn(|a| {
                f.iter().zip(directions).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]
            });
            let u = p.map(|p| p / rho);
            let stress = core::array::from_fn(|a| {
                core::array::from_fn(|b| {
                    f.iter()
                        .zip(directions)
                        .map(|(f, e)| f * e[a] * e[b])
                        .sum::<f64>()
                        + 0.5 * (u[a] * force[b] + force[a] * u[b])
                })
            });
            (rho, p, stress)
        })
        .collect()
}

#[test]
fn spatial_si_histories_keep_stress_and_supply_real_coupled_force_and_torque() {
    let fluid = fluid(3.0, 3e-4);
    let source = body(6.0);
    let mut model = assembly(&fluid, &source);
    let gate = CancelGate::new_clock_free();
    let uniform = [1e-4, -2e-4, 1.5e-4];
    model
        .initialize_uniform_flow_si([0.0002, -0.0001, 0.00015], &gate)
        .unwrap();
    model.rebind_uniform_fluid_force_si(uniform, &gate).unwrap();
    model
        .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
        .unwrap();
    for scale in [1.0, -0.7, 0.0] {
        let field = spatial_loads(scale);
        let before = physical_momentum_si(&model);
        let stress_before = physical_moments(model.grid());
        let bodies = [model.state_si(0).unwrap(), model.state_si(1).unwrap()];
        let steps = model.accepted_steps();
        model.rebind_fluid_force_field_si(&field, &gate).unwrap();
        assert_eq!(model.accepted_steps(), steps);
        assert_eq!(
            [model.state_si(0).unwrap(), model.state_si(1).unwrap()],
            bodies
        );
        for (index, (actual, old)) in physical_moments(model.grid())
            .iter()
            .zip(stress_before)
            .enumerate()
        {
            let [x, y, z] = [index % 8, (index / 8) % 8, index / 64];
            close(actual.0, old.0, 4e-15);
            for a in 0..3 {
                close(actual.1[a], old.1[a], 4e-15);
                close(
                    model.grid().force_density(x, y, z)[a],
                    (uniform[a] + field[index][a]) * DT * DT / (3.0 * DX),
                    2e-18,
                );
                for b in 0..3 {
                    close(actual.2[a][b], old.2[a][b], 4e-15);
                }
            }
        }
        for a in 0..3 {
            close(physical_momentum_si(&model)[a], before[a], 3e-19);
            close(
                model.uniform_fluid_force_density_n_m3().unwrap()[a],
                uniform[a],
                2e-18,
            );
        }
        let accepted = snapshot(&model);
        model.rebind_fluid_force_field_si(&field, &gate).unwrap();
        assert_eq!(snapshot(&model), accepted);
        let report = model
            .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
            .unwrap();
        let after = physical_momentum_si(&model);
        let imposed = imposed_impulse_si(&field, uniform);
        for a in 0..3 {
            close(after[a] - before[a], imposed[a], 2e-16);
            close(
                model
                    .impulse_ns(report.coupling.forcing.impulse[a])
                    .unwrap(),
                imposed[a],
                2e-16,
            );
            close(
                model
                    .angular_impulse_nms(report.coupling.forcing.angular_impulse[a])
                    .unwrap(),
                imposed[a + 3],
                2e-17,
            );
        }
    }
}

#[test]
fn malformed_cancelled_and_late_invalid_si_fields_preserve_history_and_retry() {
    let fluid = fluid(3.0, 3e-4);
    let source = body(6.0);
    let mut model = assembly(&fluid, &source);
    let mut control = assembly(&fluid, &source);
    let gate = CancelGate::new_clock_free();
    for model in [&mut model, &mut control] {
        model
            .rebind_uniform_fluid_force_si([1e-5; 3], &gate)
            .unwrap();
        model
            .rebind_fluid_force_field_si(&spatial_loads(1.0), &gate)
            .unwrap();
        model
            .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
            .unwrap();
    }
    let accepted = snapshot(&model);
    let force_before: Vec<_> = (0..512)
        .map(|index| {
            model
                .grid()
                .force_density(index % 8, (index / 8) % 8, index / 64)
        })
        .collect();
    let bodies = [model.state_si(0).unwrap(), model.state_si(1).unwrap()];
    let mut invalid = spatial_loads(-0.7);
    assert!(matches!(
        model.rebind_fluid_force_field_si(&invalid[..511], &gate),
        Err(MaterialImmersedError::Invalid(_))
    ));
    invalid[511][2] = f64::NAN;
    assert!(matches!(
        model.rebind_fluid_force_field_si(&invalid, &gate),
        Err(MaterialImmersedError::Invalid(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        model
            .rebind_fluid_force_field_si(&spatial_loads(-0.7), &cancelled)
            .is_err()
    );
    assert_eq!(snapshot(&model), accepted);
    for (index, force) in force_before.iter().enumerate() {
        assert_eq!(
            &model
                .grid()
                .force_density(index % 8, (index / 8) % 8, index / 64),
            force
        );
    }
    assert_eq!(
        [model.state_si(0).unwrap(), model.state_si(1).unwrap()],
        bodies
    );
    assert_eq!(model.accepted_steps(), 1);
    for model in [&mut model, &mut control] {
        model
            .rebind_fluid_force_field_si(&spatial_loads(-0.7), &gate)
            .unwrap();
        model
            .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
            .unwrap();
    }
    assert_eq!(snapshot(&model), snapshot(&control));
    assert_eq!(
        [model.state_si(0).unwrap(), model.state_si(1).unwrap()],
        [control.state_si(0).unwrap(), control.state_si(1).unwrap()]
    );
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
    let field = spatial_loads(1.0);
    let flow = [0.0001, -0.0002, 0.0001];
    single.rebind_uniform_fluid_force_si(load, &gate).unwrap();
    single.rebind_fluid_force_field_si(&field, &gate).unwrap();
    single.initialize_uniform_flow_si(flow, &gate).unwrap();
    assert_eq!(single.state_si().unwrap(), initial);
    for a in 0..3 {
        close(single.grid().velocity(3, 4, 2)[a] * DX / DT, flow[a], 2e-17);
        close(
            single.uniform_fluid_force_density_n_m3().unwrap()[a],
            load[a],
            2e-18,
        );
        close(
            single.grid().force_density(3, 4, 2)[a],
            (load[a] + field[(2 * 8 + 4) * 8 + 3][a]) * DT * DT / DX,
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
    assembly.rebind_fluid_force_field_si(&field, &gate).unwrap();
    let before = fields(assembly.grid());
    let mut planes = assembly.with_prescribed_planes(&[], &gate).unwrap();
    planes.rebind_uniform_fluid_force_si(load, &gate).unwrap();
    planes.rebind_fluid_force_field_si(&field, &gate).unwrap();
    assert_eq!(fields(planes.assembly().grid()), before);
    assert_eq!(planes.assembly().state_si(0).unwrap(), initial);
    for a in 0..3 {
        close(
            planes.uniform_fluid_force_density_n_m3().unwrap()[a],
            load[a],
            2e-18,
        );
    }
    let report = planes
        .step_forces_si(
            [[0.0; 6]],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let imposed = imposed_impulse_si(&field, load);
    for a in 0..6 {
        let full = report.motion.represented_fluid_forcing_impulse[a]
            + report.motion.interior_forcing_impulses[0][a];
        let si = if a < 3 {
            planes.assembly().impulse_ns(full).unwrap()
        } else {
            planes.assembly().angular_impulse_nms(full).unwrap()
        };
        close(si, imposed[a], 2e-16);
        close(report.motion.combined_momentum_residual[a], 0.0, 3e-10);
    }
    assert!(
        report.motion.interior_forcing_impulses[0][0] > 0.0,
        "occupied fluid receives the spatial load before the dry correction subtracts it once"
    );
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
