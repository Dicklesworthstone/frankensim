//! The same end model must survive actual solid diffusion, conductivity
//! reweighting and motion without relabeling its numerical loss as Joule heat.
use super::*;
use fs_em::{
    planar::eddy::{moving::MovingEddyControls, nonlinear::moving::NonlinearMovingEddyControls},
    winding::{
        end::EndWindingMode,
        nonlinear::{NonlinearWindingCandidate, NonlinearWindingControls},
    },
};

fn mode_driven<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
    bias: f64,
) -> NonlinearWindingCircuit<'m> {
    NonlinearWindingCircuit::new_with_solid_regions_and_end_modes(
        initial(
            problem,
            sources(card, None),
            &[4, 5, 6, 7],
            &boundary(problem, bias),
        ),
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 777.0 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles(),
            end_inductance_h: vec![0.0],
        },
        vec![port()],
        vec![EndWindingMode {
            inductance_h: 0.05,
            turns: vec![2.0],
        }],
        vec![],
        bindings(),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap()
}

fn moving_controls() -> NonlinearWindingControls<NonlinearMovingEddyControls> {
    let base = controls();
    NonlinearWindingControls {
        field: NonlinearMovingEddyControls {
            moving: MovingEddyControls {
                relative_tolerance: 1e-12,
                absolute_residual_a_m: 1e-13,
                restart: 8,
                max_cycles: 50,
                max_krylov_scalars: 10_000,
                max_cell_magnetic_peclet: 2.0,
                energy_tolerance_j: 1e-9,
                current_balance_tolerance_a: 1e-9,
            },
            absolute_residual_tolerance_a_m: 1e-11,
            relative_residual_tolerance: 1e-11,
            max_newton_iterations: 30,
            max_backtracks: 30,
        },
        maximum_field_trials: base.maximum_field_trials,
        maximum_newton_iterations: base.maximum_newton_iterations,
        difference_step_a: base.difference_step_a,
        terminal_voltage_tolerance_v: base.terminal_voltage_tolerance_v,
        power_tolerance_w: base.power_tolerance_w,
        energy_tolerance_j: base.energy_tolerance_j,
        newton: base.newton,
    }
}

fn same<F>(a: &NonlinearWindingCandidate<'_, F>, b: &NonlinearWindingCandidate<'_, F>) {
    near(a.circuit().currents_a()[0], b.circuit().currents_a()[0]);
    near(
        a.circuit().flux_linkages_wb()[0],
        b.circuit().flux_linkages_wb()[0],
    );
    near(
        a.circuit().resistance_ohm()[0],
        b.circuit().resistance_ohm()[0],
    );
    for (&x, &y) in a
        .circuit()
        .field()
        .potential_wb_m()
        .iter()
        .zip(b.circuit().field().potential_wb_m())
    {
        near(x, y);
    }
    for (&x, &y) in a.circuit().profiles_per_m2()[0]
        .iter()
        .zip(&b.circuit().profiles_per_m2()[0])
    {
        near(x, y);
    }
    let (a, b) = (a.step(), b.step());
    near(a.terminal_voltage_v[0], b.terminal_voltage_v[0]);
    near(a.inductive_voltage_v[0], b.inductive_voltage_v[0]);
    near(a.motional_voltage_v[0], b.motional_voltage_v[0]);
    near(a.stored_before_j, b.stored_before_j);
    near(a.stored_after_j, b.stored_after_j);
    near(a.stored_change_j, b.stored_change_j);
    near(a.heat_j, b.heat_j);
    near(a.numerical_dissipation_j, b.numerical_dissipation_j);
    near(a.mechanical_work_j, b.mechanical_work_j);
    near(a.supplied_work_j, b.supplied_work_j);
    assert!(a.cells.is_empty() && b.cells.is_empty());
    assert!(a.energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(b.energy_residual_j.abs() <= controls().energy_tolerance_j);
}

#[test]
fn end_modes_survive_solid_diffusion_reweighted_sources_and_nonzero_motion() {
    let card = common::source(true, SIGMA);
    let problem = common::material::problem(&card, mesh());
    let bias = 0.15;
    let mut diagonal = driven(&problem, &card, bias, 0.2);
    let original = mode_driven(&problem, &card, bias);
    let mut modal = original.clone();
    assert_eq!(modal.end_inductance_matrix_h(), [0.2]);
    let mut old_current = 0.0;
    for (target, imposed, hot_triangle) in [(0.05, 3.0, None), (0.1, 0.8, Some(0))] {
        let a = diagonal
            .step_to_candidate(
                target,
                sources(&card, hot_triangle),
                vec![],
                &[imposed],
                controls(),
                &mut || false,
            )
            .unwrap();
        let b = modal
            .step_to_candidate(
                target,
                sources(&card, hot_triangle),
                vec![],
                &[imposed],
                controls(),
                &mut || false,
            )
            .unwrap();
        same(&a, &b);
        near(
            b.step().numerical_dissipation_j - b.step().field.numerical_dissipation_j,
            0.1 * (imposed - old_current).powi(2),
        );
        assert_eq!(b.circuit().end_inductance(), original.end_inductance());
        if hot_triangle.is_some() {
            assert_ne!(b.circuit().profiles_per_m2(), original.profiles_per_m2());
        }
        old_current = imposed;
        diagonal = a.into_parts().0;
        modal = b.into_parts().0;
    }
    // Uniform endpoint conductivity admits nonzero Eulerian through-flow in
    // each solid region. The old diffusion cochain is retained across refresh.
    let hot = common::material::conductor(&card, 500.0, ElectricalPropertyBasis::Resistivity);
    let refreshed: Vec<_> = sources(&card, None)
        .into_iter()
        .map(|mut s| {
            s.conductor = hot.clone();
            s
        })
        .collect();
    let velocity: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .map(|p| {
            [
                if p[0] <= 1.0 {
                    0.04
                } else if p[0] >= 2.0 {
                    -0.04
                } else {
                    0.0
                },
                0.0,
            ]
        })
        .collect();
    let imposed = -0.3;
    let a = diagonal
        .step_with_eulerian_velocity_and_sources_candidate(
            0.2,
            &velocity,
            refreshed.clone(),
            vec![],
            &[imposed],
            moving_controls(),
            &mut || false,
        )
        .unwrap();
    let b = modal
        .step_with_eulerian_velocity_and_sources_candidate(
            0.2,
            &velocity,
            refreshed,
            vec![],
            &[imposed],
            moving_controls(),
            &mut || false,
        )
        .unwrap();
    same(&a, &b);
    near(
        b.step().numerical_dissipation_j - b.step().field.numerical_dissipation_j,
        0.1 * (imposed - old_current).powi(2),
    );
    assert!(b.step().heat_j > 0.0);
    assert!(b.step().motional_voltage_v[0].abs() > 1e-4);
    assert_eq!(b.circuit().end_inductance(), original.end_inductance());
    assert_eq!(original.currents_a(), [0.0]);
    assert_eq!(original.time_s(), 0.0);
}
