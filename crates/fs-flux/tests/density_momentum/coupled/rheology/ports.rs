//! Nonzero-shear pressure responses and finite RCR history at the same endpoint.

use super::*;
use fs_flux::material::{
    MaterialStokesWindkesselHistory as History, MaterialStokesWindkesselPort as Port,
};
use fs_flux::ns::boundary::conservative::coupled::rheology::CoupledDensityMomentumRheologyStep;
use fs_flux::ns::boundary::conservative::coupled::windkessel::DensityMomentumWindkesselError as RcrError;

fn port_curve(cell: &Cell<1>) -> CarreauYasuda {
    CarreauYasuda::new(
        0.6 + 0.02 * cell.density_kg_m3 + 0.002 * cell.specific[0],
        0.05,
        4.0,
        0.45,
        2.0,
        [0.0, 100.0],
    )
    .unwrap()
}

fn sides(mesh: &TriMesh) -> (Vec<usize>, Vec<usize>) {
    let side = |x| {
        (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e) && mesh.edges[e].mid[0] == x)
            .collect()
    };
    (side(0.0), side(1.0))
}

fn volume_rate(system: &FluxSystem<'_>, port: &[usize], velocity: &[f64]) -> f64 {
    port.iter()
        .map(|&e| system.mesh.edges[e].len * velocity[2 * e])
        .sum()
}

fn check_endpoint(
    system: &FluxSystem<'_>,
    accepted: &CoupledDensityMomentumRheologyStep<1>,
    old: &[Cell<1>],
    inflow: &[Option<Cell<1>>],
    dt: f64,
    control: Controls,
    samples: &[(Cell<1>, f64, f64)],
) {
    assert!(accepted.step.iterations > 1);
    assert!(accepted.step.nonlinear_residual <= nonlinear().relative_tolerance);
    assert!(accepted.viscous_heating_w > 0.0);
    let step = &accepted.step.step;
    assert!(step.cells.iter().zip(old).any(|(a, b)| a != b));
    assert!(accepted.cells.iter().any(|c| c.strain_rate_per_s > 0.05));
    for (cell, response) in accepted.cells.iter().enumerate() {
        assert_eq!(samples[cell].0, step.cells[cell]);
        assert_eq!(samples[cell].1, response.strain_rate_per_s);
        assert_eq!(samples[cell].2, response.dynamic_viscosity_pa_s);
        let rate = system
            .try_cell_kinematics(&step.flow.field.x, cell)
            .unwrap()
            .shear_rate();
        near(response.strain_rate_per_s, rate, 2e-13);
        near(
            response.dynamic_viscosity_pa_s,
            port_curve(&step.cells[cell]).viscosity_pa_s(rate).unwrap(),
            2e-13,
        );
    }
    let plan = system
        .prepare_density_transport(
            &step.flow.field.x[..system.n_u],
            control.depth_m,
            control.relative_volume_divergence_tolerance,
            control.maximum_transport_storage_bytes,
            || false,
        )
        .unwrap();
    let mut scratch = plan
        .workspace(control.maximum_transport_storage_bytes)
        .unwrap();
    let mut transported = vec![Cell::default(); old.len()];
    let report = plan
        .advance_into(
            old,
            inflow,
            dt,
            control.maximum_courant,
            &mut transported,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(step.cells, transported);
    assert_eq!(step.transport, report);
    near(report.mass_balance_residual_kg, 0.0, 2e-12);
    near(report.carried_balance_residual[0], 0.0, 2e-11);
}

#[test]
fn controlled_sheared_channel_resolves_density_rheology_and_port_pressure_together() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let (west, east) = sides(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            if west.contains(&edge) || east.contains(&edge) {
                Boundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                Boundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let original: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| {
            if p[0] < 0.5 {
                Cell {
                    density_kg_m3: 2.0,
                    specific: [3.0],
                }
            } else {
                Cell {
                    density_kg_m3: 8.0,
                    specific: [19.0],
                }
            }
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| {
            mesh.is_boundary(e).then_some(Cell {
                density_kg_m3: 6.0,
                specific: [9.0],
            })
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    // Annotated: on macOS an extra `IntoIterator` impl in the dependency
    // graph leaves `.zip(velocity)` below without a unique type otherwise.
    let mut invariant: Option<(Vec<f64>, f64, f64)> = None;
    for reference in [1.0, 13.0] {
        let samples = RefCell::new(vec![(Cell::default(), 0.0, 0.0); original.len()]);
        let force_calls = Counter::new(0);
        let source_calls = Counter::new(0);
        let source = |cell, point, state: &Cell<1>, rate| {
            assert_eq!(point, mesh.centroids[cell]);
            source_calls.set(source_calls.get() + 1);
            let mu = port_curve(state).viscosity_pa_s(rate)?;
            samples.borrow_mut()[cell] = (*state, rate, mu);
            Ok(mu)
        };
        let control = controls(reference, 1.7);
        let solved = system
            .density_momentum_coupled_bdf1_rheology_flow_control_step(
                numerics(),
                &source,
                &|_, _, _| {
                    force_calls.set(force_calls.get() + 1);
                    [0.0; 2]
                },
                &boundaries,
                &[&west],
                &[-0.15],
                2e-11,
                &previous,
                &original,
                &inflow,
                0.06,
                control,
                nonlinear(),
                || false,
            )
            .unwrap();
        check_endpoint(
            &system,
            &solved.step,
            &original,
            &inflow,
            0.06,
            control,
            &samples.borrow(),
        );
        // One source query per cell/assembly, with no source query in native
        // pressure responses or reporting; force uses six native volume nodes.
        assert_eq!(force_calls.get(), 6 * source_calls.get());
        let step = &solved.step.step.step;
        near(solved.outward_volume_rates_m2_s[0], -0.15, 2e-11);
        near(
            volume_rate(&system, &west, &step.flow.field.x),
            -0.15,
            2e-11,
        );
        near(volume_rate(&system, &east, &step.flow.field.x), 0.15, 2e-11);
        assert!(solved.maximum_flow_residual_m2_s <= 2e-11);
        let reaction = step
            .flow
            .reactions
            .resultant(&west, [0.0; 2], &[None])
            .unwrap();
        near(
            reaction.imposed.force[0] * reference,
            -solved.pressure_offsets_pa[0],
            2e-8,
        );
        let plateau = system
            .density_momentum_coupled_bdf1_flow_control_step(
                numerics(),
                &|_, _, state| port_curve(state).viscosity_pa_s(0.0).unwrap(),
                &|_, _, _| [0.0; 2],
                &boundaries,
                &[&west],
                &[-0.15],
                2e-11,
                &previous,
                &original,
                &inflow,
                0.06,
                control,
                nonlinear(),
                || false,
            )
            .unwrap();
        assert!(solved.pressure_offsets_pa[0] + 1e-4 < plateau.pressure_offsets_pa[0]);
        let values = (
            step.flow.field.x[..system.n_u].to_vec(),
            solved.pressure_offsets_pa[0],
            solved.step.viscous_heating_w,
        );
        if let Some((velocity, pressure, heat)) = &invariant {
            for (&a, &b) in values.0.iter().zip(velocity) {
                near(a, b, 2e-10);
            }
            near(values.1, *pressure, 2e-8);
            near(values.2, *heat, 2e-10);
        } else {
            invariant = Some(values);
        }
    }
    assert!(previous.iter().all(|v| *v == 0.0));
}

#[test]
fn sheared_rcr_flow_preserves_original_storage_laws_across_two_density_updates() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 2));
    let system = FluxSystem::new(&mesh);
    let (west, east) = sides(&mesh);
    let ports = [&west[..], &east[..]];
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| {
            if west.contains(&edge) || east.contains(&edge) {
                Boundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                Boundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let initial: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| {
            if p[0] < 0.5 {
                Cell {
                    density_kg_m3: 2.0,
                    specific: [3.0],
                }
            } else {
                Cell {
                    density_kg_m3: 8.0,
                    specific: [19.0],
                }
            }
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|e| {
            mesh.is_boundary(e).then_some(Cell {
                density_kg_m3: 6.0,
                specific: [9.0],
            })
        })
        .collect();
    let loads = [
        Port {
            reservoir_pressure_pa: 12.0,
            proximal_resistance_pa_s_per_m2: 1.0,
            distal_resistance_pa_s_per_m2: 20.0,
            compliance_m2_per_pa: 0.03,
            inertance_pa_s2_per_m2: 0.2,
        },
        Port {
            reservoir_pressure_pa: 0.0,
            proximal_resistance_pa_s_per_m2: 2.0,
            distal_resistance_pa_s_per_m2: 30.0,
            compliance_m2_per_pa: 0.02,
            inertance_pa_s2_per_m2: 0.3,
        },
    ];
    let mut invariant = Vec::new();
    for reference in [1.0, 13.0] {
        let mut old_cells = initial.clone();
        let mut old_velocity = vec![0.0; system.n_u];
        let mut old_x = [1.0, 0.2];
        for (interval, dt) in [0.08, 0.05].into_iter().enumerate() {
            let saved = (old_cells.clone(), old_velocity.clone(), old_x);
            let samples = RefCell::new(vec![(Cell::default(), 0.0, 0.0); old_cells.len()]);
            let source = |cell, _, state: &Cell<1>, rate| {
                let mu = port_curve(state).viscosity_pa_s(rate)?;
                samples.borrow_mut()[cell] = (*state, rate, mu);
                Ok(mu)
            };
            let control = controls(reference, 1.7);
            let solved = system
                .density_momentum_coupled_bdf1_rheology_windkessel_step(
                    numerics(),
                    &source,
                    &|_, _, _| [0.0; 2],
                    &boundaries,
                    &ports,
                    &loads,
                    History {
                        velocity: &old_velocity,
                        compliance_pressures_pa: &old_x,
                    },
                    &old_cells,
                    &inflow,
                    dt,
                    control,
                    nonlinear(),
                    2e-9,
                    2e-11,
                    || false,
                )
                .unwrap();
            check_endpoint(
                &system,
                &solved.flow,
                &old_cells,
                &inflow,
                dt,
                control,
                &samples.borrow(),
            );
            assert_eq!((old_cells.clone(), old_velocity.clone(), old_x), saved);
            let step = &solved.flow.step.step;
            assert!(solved.outward_volume_rates_m2_s[1] > 0.0);
            for i in 0..2 {
                let q = volume_rate(&system, ports[i], &step.flow.field.x);
                let old_q = volume_rate(&system, ports[i], &old_velocity);
                let load = loads[i];
                let balance = solved.balances[i];
                let x = (load.compliance_m2_per_pa * old_x[i] / dt + q)
                    / (load.compliance_m2_per_pa / dt + 1.0 / load.distal_resistance_pa_s_per_m2);
                let pressure = load.reservoir_pressure_pa
                    + load.proximal_resistance_pa_s_per_m2 * q
                    + load.inertance_pa_s2_per_m2 * (q - old_q) / dt
                    + x;
                near(solved.outward_volume_rates_m2_s[i], q, 2e-11);
                near(balance.compliance_pressure_pa, x, 2e-9);
                near(solved.pressure_offsets_pa[i], pressure, 2e-8);
                near(
                    balance.compliance_energy_j_per_m,
                    0.5 * load.compliance_m2_per_pa * x * x,
                    2e-10,
                );
                near(
                    balance.inertance_energy_j_per_m,
                    0.5 * load.inertance_pa_s2_per_m2 * q * q,
                    2e-10,
                );
                let defect = (load.compliance_m2_per_pa * (x - old_x[i]).powi(2)
                    + load.inertance_pa_s2_per_m2 * (q - old_q).powi(2))
                    / (2.0 * dt);
                near(balance.bdf_energy_defect_w_per_m, defect, 2e-10);
                near(balance.port_power_w_per_m, pressure * q, 2e-8);
                near(balance.power_balance_residual_w_per_m, 0.0, 2e-9);
                near(balance.pressure_law_residual_pa, 0.0, 2e-8);
                near(balance.compliance_flow_residual_m2_s, 0.0, 2e-11);
            }
            if interval == 0 {
                let plateau = system
                    .density_momentum_coupled_bdf1_windkessel_step(
                        numerics(),
                        &|_, _, state| port_curve(state).viscosity_pa_s(0.0).unwrap(),
                        &|_, _, _| [0.0; 2],
                        &boundaries,
                        &ports,
                        &loads,
                        History {
                            velocity: &old_velocity,
                            compliance_pressures_pa: &old_x,
                        },
                        &old_cells,
                        &inflow,
                        dt,
                        control,
                        nonlinear(),
                        2e-9,
                        2e-11,
                        || false,
                    )
                    .unwrap();
                assert!(
                    solved.outward_volume_rates_m2_s[1]
                        > plateau.outward_volume_rates_m2_s[1] + 1e-5
                );
            }
            let values = (
                step.flow.field.x[..system.n_u].to_vec(),
                solved.pressure_offsets_pa.clone(),
                solved.flow.viscous_heating_w,
            );
            if reference == 1.0 {
                invariant.push(values);
            } else {
                for (&a, &b) in values.0.iter().zip(&invariant[interval].0) {
                    near(a, b, 2e-10);
                }
                for (&a, &b) in values.1.iter().zip(&invariant[interval].1) {
                    near(a, b, 2e-8);
                }
                near(values.2, invariant[interval].2, 2e-10);
            }
            old_x = [
                solved.balances[0].compliance_pressure_pa,
                solved.balances[1].compliance_pressure_pa,
            ];
            old_velocity = step.flow.field.x[..system.n_u].to_vec();
            old_cells = step.cells.clone();
        }
    }
    // A late source-domain refusal must survive both native hydraulic reduction
    // and physical result conversion, without advancing the authored history.
    let history = vec![0.0; system.n_u];
    let old_x = [1.0, 0.2];
    let narrow = CarreauYasuda::new(0.8, 0.05, 4.0, 0.45, 2.0, [0.0, 0.01]).unwrap();
    let calls = Counter::new(0);
    let error = system
        .density_momentum_coupled_bdf1_rheology_windkessel_step(
            numerics(),
            &|_, _, _, rate| {
                calls.set(calls.get() + 1);
                narrow.viscosity_pa_s(rate)
            },
            &|_, _, _| [0.0; 2],
            &boundaries,
            &ports,
            &loads,
            History {
                velocity: &history,
                compliance_pressures_pa: &old_x,
            },
            &initial,
            &inflow,
            0.08,
            controls(13.0, 1.7),
            nonlinear(),
            2e-9,
            2e-11,
            || false,
        )
        .err()
        .unwrap();
    assert!(calls.get() > initial.len());
    assert!(
        matches!(error, RcrError::Momentum(Error::Material { error, .. }) if error == narrow.viscosity_pa_s(1.0).unwrap_err())
    );
    assert!(history.iter().all(|v| *v == 0.0));
    assert_eq!(old_x, [1.0, 0.2]);
}
