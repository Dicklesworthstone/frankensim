//! Both implicit live-flow entries retain the selected thermal face scheme.

use super::*;
use fs_lbm::conjugate::ViscousHeatControls;
use fs_lbm::conjugate::enthalpy::advection::feedback::implicit::{
    ImplicitFlowControls, ImplicitUnsteadyFlowControls,
};

fn controls() -> ImplicitFlowControls {
    ImplicitFlowControls {
        max_unknowns: 4,
        max_newton_iterations: 32,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 4e-7,
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 2e-8,
            relative_tolerance: 0.0,
            linear_restart: 4,
            ..Default::default()
        },
    }
}

fn heating() -> ViscousHeatControls {
    ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-7,
    }
}

#[test]
fn implicit_pressure_flow_uses_reconstructed_endpoint_h_with_source_mu_and_physical_work() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Pressure {
        pressure_pa: 100020.0,
    };
    faces[1] = FvBoundary::Pressure {
        pressure_pa: 100000.0,
    };
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::wall();
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-13;
    flow.pressure_tolerance = 1e-12;
    flow.momentum_tolerance = 1e-10;
    let card = feedback_card(100.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old_h = [20.0, 40.0, 80.0, 140.0];
    let old = profile(&domain, 1.0, &caloric, old_h);
    let plan = plan(&old);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow, (), &gate, 2048)
        .unwrap()
        .with_pressure_work_heating(heating())
        .unwrap()
        .with_reconstruction(Arc::clone(&plan))
        .unwrap();
    let dt = 0.2;
    let (next, report) = observer
        .step_implicit(&old, dt, controls(), &[0.0; 4], || Ok(()))
        .unwrap();
    let h: [f64; 4] = std::array::from_fn(|c| next.cells()[&CellId(c as u64)].specific_h_j_kg);
    let q = report.endpoint.flow.field.flux_x(2, 0, 0);
    assert!(q > 0.0);
    let face = traces(h, true, 0.0);
    let power = report.endpoint.viscous_power.as_ref().unwrap();
    for c in 0..4 {
        let residual = h[c] - old_h[c] - dt * (power.cell_heat_w[c] + q * (face[c] - face[c + 1]));
        assert!(residual.abs() < 1e-7, "cell {c} residual {residual}");
        near(
            report.endpoint.samples[c]
                .transport
                .dynamic_viscosity_pa_s(),
            viscosity(100.0, h[c]),
        );
    }
    // The nonzero correction would fail these local balances if an implicit
    // adapter silently discarded the selected scheme and used donor h.
    assert!((face[2] - h[1]).abs() > 1.0);
    near(report.net_stream_enthalpy_j, -dt * q * h[3]);
    near(
        next.energy_j() - old.energy_j(),
        report.mechanical.generated_heat_j + report.net_stream_enthalpy_j,
    );
    assert!(report.combined_residual_j.abs() < controls().total_energy_tolerance_j);
    assert!(std::ptr::eq(
        observer.reconstruction().unwrap(),
        plan.as_ref()
    ));
    assert!(observer.applied().is_none());
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn implicit_retained_momentum_uses_reconstructed_endpoint_traces_and_can_resume() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-13;
    flow.pressure_tolerance = 1e-12;
    flow.momentum_tolerance = 1e-10;
    let card = card(claims(1.0, 100.0, false));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old_h = [200.0, 100.0, 100.0, 0.0];
    let original = profile(&domain, 1.0, &caloric, old_h);
    let cycle = [0usize, 2, 3, 1];
    // An authored closed centreline: each segment's midpoint is its actual
    // voxel centre. The original plan follows arclength around these bends.
    let vertices = [
        [1.0, 0.0, 0.5],
        [0.0, 1.0, 0.5],
        [1.0, 2.0, 0.5],
        [2.0, 1.0, 0.5],
    ];
    let events = cycle
        .into_iter()
        .enumerate()
        .map(|(i, c)| {
            let mut cell = original.cells()[&CellId(c as u64)].clone();
            cell.geometry.start_m = vertices[i];
            cell.geometry.end_m = vertices[(i + 1) % 4];
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    let old = DepositionSimulation::new(
        events,
        &[(MAT, caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 64,
            max_events: 4,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0;
    let stencils: Vec<_> = (0..4)
        .map(|i| AxialCellStencil {
            previous: CellId(cycle[(i + 3) % 4] as u64),
            cell: CellId(cycle[i] as u64),
            next: CellId(cycle[(i + 1) % 4] as u64),
        })
        .collect();
    let plan =
        Arc::new(LimitedEnthalpyReconstruction::try_new(&old, &stencils, 4, || Ok(())).unwrap());
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let force = vec![
        [-100.0, 100.0, 0.0],
        [-100.0, -100.0, 0.0],
        [100.0, 100.0, 0.0],
        [100.0, -100.0, 0.0],
    ];
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow, (), &gate, 2048)
        .unwrap()
        .with_unsteady_body_force_heating(force, heating())
        .unwrap()
        .with_reconstruction(Arc::clone(&plan))
        .unwrap();
    let mut current = old.clone();
    let mut retained = observer.clone();
    let mut old_q = 0.0;
    for target in [0.5, 1.0] {
        let dt = target - current.time_s();
        let q = (400.0 + 4.0 * old_q / dt) / (4.0 / dt + 2400.0);
        let (next, accepted, report) = retained
            .step_unsteady_implicit(
                &current,
                target,
                ImplicitUnsteadyFlowControls {
                    thermal: controls(),
                    rheology: None,
                },
                &[1.0; 4],
                || Ok(()),
            )
            .unwrap();
        near(accepted.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q);
        let h: [f64; 4] = cycle.map(|c| next.cells()[&CellId(c as u64)].specific_h_j_kg);
        let mut trace = h;
        for i in 0..4 {
            let left = h[i] - h[(i + 3) % 4];
            let right = h[(i + 1) % 4] - h[i];
            if left * right > 0.0 {
                trace[i] += 0.5 * left.signum() * left.abs().min(right.abs());
            }
        }
        assert!(
            trace
                .iter()
                .zip(h)
                .any(|(trace, h)| (trace - h).abs() > 1.0)
        );
        for (i, c) in cycle.into_iter().enumerate() {
            let stored = h[i] - current.cells()[&CellId(c as u64)].specific_h_j_kg;
            let watts = 600.0 * q * q + q * (trace[(i + 3) % 4] - trace[i]);
            assert!((stored - dt * watts).abs() < 1e-7);
        }
        near(report.deposited_heat_j, dt * 2400.0 * q * q);
        near(
            report.thermal.stored_energy_change_j,
            report.deposited_heat_j,
        );
        near(
            report.mechanical.initial_kinetic_energy_j,
            2.0 * old_q * old_q,
        );
        near(report.mechanical.final_kinetic_energy_j, 2.0 * q * q);
        near(
            report.mechanical.backward_euler_loss_j,
            2.0 * (q - old_q).powi(2),
        );
        assert!(report.combined_residual_j.abs() < controls().total_energy_tolerance_j);
        assert!(std::ptr::eq(
            accepted.reconstruction().unwrap(),
            plan.as_ref()
        ));
        current = next;
        retained = accepted;
        old_q = q;
    }
    assert!(observer.applied().is_none());
    near(observer.unsteady_flow().unwrap().field.flux_x(1, 1, 0), 0.0);
    assert_eq!(old.time_s(), 0.0);
}
