//! Actual native circuits must respond to spatial H(B) and conductor data.
#[path = "nonlinear_winding_circuit/common.rs"]
mod common;
#[path = "nonlinear_winding_circuit/thermal.rs"]
mod thermal;

use common::{b_from_current, conductors, controls, density_energy, h, linkage, near};
use fs_em::{
    planar::eddy::{EddyConductor, nonlinear::PlanarNonlinearEddy},
    winding::nonlinear::{
        NonlinearWindingCircuit, NonlinearWindingError, NonlinearWindingInitial, NonlinearWindings,
    },
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::{Branch, CircuitGraph};

fn bisect(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low) <= 0.0 && residual(high) >= 0.0);
    for _ in 0..80 {
        let mid = f64::midpoint(low, high);
        if residual(mid) > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    f64::midpoint(low, high)
}

#[test]
fn source_h_b_and_resistivity_change_voltage_driven_current_field_heat_and_continuation() {
    let depth = 0.3;
    let dt = 0.2;
    let graph_r = 0.3;
    let mut material_outcomes = Vec::new();
    for saturating in [false, true] {
        for sigma in [18.0, 36.0] {
            for sign in [1.0, -1.0] {
                let card = common::source(saturating, sigma);
                let problem = common::problem(&card, depth);
                let original = common::single(&problem, &card, graph_r, 0.0, 0.0);
                let rw = 72.0 * depth / sigma;
                near(original.resistance_ohm()[0], rw);
                let voltage = 9.0 * sign;
                let magnitude = bisect(0.0, if saturating { 13.0 } else { 4.0 }, |i| {
                    (rw + graph_r) * i + linkage(i, depth, saturating) / dt - 9.0
                });
                let expected = sign * magnitude;
                let candidate = original
                    .step_to_candidate(
                        dt,
                        vec![],
                        conductors(&card, 300.0),
                        &[voltage],
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                let step = candidate.step();
                let next = candidate.circuit();
                near(next.currents_a()[0], expected);
                let b = b_from_current(expected, saturating);
                near(next.field().potential_wb_m()[6], b / 2.0);
                near(next.field().potential_wb_m()[7], -b / 2.0);
                near(next.flux_linkages_wb()[0], 2.0 * depth * b);
                near(step.terminal_voltage_v[0], voltage - graph_r * expected);
                let physical_heat = (rw + graph_r) * expected * expected * dt;
                let storage = 2.0 * depth * density_energy(b, saturating);
                let numerical = expected * linkage(expected, depth, saturating) - storage;
                near(step.stored_after_j, storage);
                near(step.stored_before_j, 0.0);
                near(step.heat_j, physical_heat);
                near(
                    step.cells.iter().map(|c| c.heat_j).sum(),
                    rw * expected * expected * dt,
                );
                near(
                    step.graph_resistors[0].heat_j,
                    graph_r * expected * expected * dt,
                );
                near(step.numerical_dissipation_j, numerical);
                assert!(numerical > 0.1);
                near(step.supplied_work_j, voltage * expected * dt);
                // Private terminal work is NOT external supply. The remaining
                // native circuit stores no energy and dissipates only graph R.
                near(step.native.supplied, step.graph_resistors[0].heat_j);
                near(step.native.dissipated, step.graph_resistors[0].heat_j);
                assert!((step.native.supplied - step.supplied_work_j).abs() > 0.1);
                assert!(step.field.cells.is_empty());
                near(step.field.heat_j, 0.0);
                assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
                assert!(
                    step.terminal_residual_v[0].abs() < controls().terminal_voltage_tolerance_v
                );
                assert!(step.field_trials > 1);
                for cell in &step.cells {
                    near(cell.electric_field_v_m, cell.current_density_a_m2 / sigma);
                    near(
                        cell.power_w,
                        cell.current_density_a_m2.powi(2) * cell.volume_m3 / sigma,
                    );
                    assert_eq!(
                        cell.source.conductor.material().card_identity(),
                        card.content_hash()
                    );
                }
                let dt2 = 0.1;
                let next_magnitude = bisect(magnitude, if saturating { 13.0 } else { 4.0 }, |i| {
                    (rw + graph_r) * i
                        + (linkage(i, depth, saturating) - linkage(magnitude, depth, saturating))
                            / dt2
                        - 9.0
                });
                let second = next
                    .step_to_candidate(
                        dt + dt2,
                        vec![],
                        conductors(&card, 300.0),
                        &[voltage],
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                near(second.circuit().currents_a()[0], sign * next_magnitude);
                near(second.step().stored_before_j, step.stored_after_j);
                near(second.circuit().time_s(), dt + dt2);
                assert!(core::ptr::eq(second.circuit().field().problem(), &problem));
                near(original.time_s(), 0.0);
                near(original.currents_a()[0], 0.0);
                assert!(original.field().potential_wb_m().iter().all(|v| *v == 0.0));
                if sigma == 18.0 && sign == 1.0 {
                    material_outcomes.push((expected, storage, physical_heat));
                }
            }
        }
    }
    near(material_outcomes[0].0, 2.0);
    near(material_outcomes[1].0, 3.0);
    assert!(material_outcomes[1].1 < material_outcomes[0].1);
    assert!(material_outcomes[1].2 > 2.0 * material_outcomes[0].2);
}

#[test]
fn two_spatial_windings_and_native_capacitor_match_an_independent_mutual_circuit_oracle() {
    let depth = 0.3;
    let dt = 0.2;
    let vs = 9.0;
    let rp = 0.3;
    let rs = 0.7;
    let c = 0.4;
    let q0 = 0.8;
    let ends = [0.1, 0.2];
    let rw = 2.4;
    for saturating in [false, true] {
        let card = common::source(saturating, 18.0);
        let problem = common::problem(&card, depth);
        let profiles = common::profiles(2);
        let original = NonlinearWindingCircuit::new(
            common::initial_field(&problem, &profiles, &[0.0, 0.0]),
            CircuitGraph {
                node_count: 5,
                branches: vec![
                    (1, 0, Branch::VoltageSource { port: 0 }),
                    (1, 2, Branch::Inductor { henries: 777.0 }),
                    (2, 0, Branch::Resistor { ohms: rp }),
                    (3, 0, Branch::Inductor { henries: 999.0 }),
                    (3, 4, Branch::Capacitor { farads: c }),
                    (4, 0, Branch::Resistor { ohms: rs }),
                ],
                transformers: vec![],
            },
            NonlinearWindings {
                branch_indices: vec![1, 3],
                profiles_per_m2: profiles,
                end_inductance_h: ends.to_vec(),
            },
            conductors(&card, 300.0),
            NonlinearWindingInitial {
                currents_a: vec![0.0, 0.0],
                remaining_fluxes_wb: vec![],
                charges_c: vec![q0],
                external_inputs: vec![0.0],
            },
            common::config(),
            &mut || false,
        )
        .unwrap();
        // Exact coarse-mesh magnetic relation H(2a)=I1+I2 and native midpoint
        // charge equation q1=q0-dt*I2 reduce the coupled problem to one scalar.
        let a1 = rw + rp + ends[0] / dt;
        let a2 = rw + rs + dt / (2.0 * c) + ends[1] / dt;
        let total = bisect(-4.0, 4.0, |s| {
            s - vs / a1 - q0 / (c * a2) + linkage(s, depth, saturating) / dt * (1.0 / a1 + 1.0 / a2)
        });
        let g = linkage(total, depth, saturating) / dt;
        let currents = [(vs - g) / a1, (q0 / c - g) / a2];
        near(currents[0] + currents[1], total);
        assert!(currents[0] > 0.0 && currents[1] < 0.0);
        let q1 = q0 - dt * currents[1];
        let candidate = original
            .step_to_candidate(
                dt,
                vec![],
                conductors(&card, 300.0),
                &[vs],
                controls(),
                &mut || false,
            )
            .unwrap();
        let next = candidate.circuit();
        let step = candidate.step();
        for i in 0..2 {
            near(next.currents_a()[i], currents[i]);
            near(
                next.flux_linkages_wb()[i],
                linkage(total, depth, saturating) + ends[i] * currents[i],
            );
            near(next.resistance_ohm()[i], rw);
        }
        near(next.state()[next.dae().charge_index[0]], q1);
        near(step.terminal_voltage_v[0], vs - rp * currents[0]);
        near(
            step.terminal_voltage_v[1],
            (q0 + q1) / (2.0 * c) - rs * currents[1],
        );
        let b = b_from_current(total, saturating);
        let magnetic = 2.0 * depth * density_energy(b, saturating);
        let end_energy = 0.5 * ends[0] * currents[0].powi(2) + 0.5 * ends[1] * currents[1].powi(2);
        let capacitor_change = (q1 * q1 - q0 * q0) / (2.0 * c);
        let physical = dt * ((rw + rp) * currents[0].powi(2) + (rw + rs) * currents[1].powi(2));
        let numerical = total * linkage(total, depth, saturating) - magnetic + end_energy;
        near(step.native.delta_h, capacitor_change);
        near(step.stored_before_j, q0 * q0 / (2.0 * c));
        near(
            step.stored_after_j,
            magnetic + end_energy + q1 * q1 / (2.0 * c),
        );
        near(step.heat_j, physical);
        near(step.numerical_dissipation_j, numerical);
        near(step.supplied_work_j, vs * currents[0] * dt);
        near(
            magnetic + end_energy + capacitor_change + physical + numerical,
            vs * currents[0] * dt,
        );
        near(next.field().potential_wb_m()[6], b / 2.0);
        near(next.field().potential_wb_m()[7], -b / 2.0);
        assert!(
            step.cells
                .iter()
                .any(|cell| cell.winding == 1 && cell.heat_j > 0.0)
        );
        assert_eq!(
            step.graph_resistors
                .iter()
                .map(|r| r.branch)
                .collect::<Vec<_>>(),
            vec![2, 5]
        );
    }
}

#[test]
fn winding_and_resolved_core_eddy_heat_are_distinct_and_match_the_same_scalar_field() {
    let depth = 0.3;
    let dt = 0.2;
    let sigma = 18.0;
    let card = common::source(true, sigma);
    let problem = common::problem(&card, depth);
    let profiles = vec![common::profiles(2).remove(0)];
    let zeros = problem
        .solve(
            &[0.0; 8],
            &common::boundary(),
            &[0.0; 8],
            common::material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let material = common::material::conductor(&card, 300.0, ElectricalPropertyBasis::Resistivity);
    let core: Vec<_> = [2, 3, 6, 7]
        .into_iter()
        .map(|triangle| EddyConductor {
            triangle,
            region: if triangle < 4 { 7 } else { 8 },
            conductor: material.clone(),
        })
        .collect();
    let field = PlanarNonlinearEddy::new_with_impressed_current(
        &zeros,
        core.clone(),
        &[0, 1, 4, 5],
        &common::boundary(),
        &[0.0; 8],
        0.0,
        common::material::controls(),
        &mut || false,
    )
    .unwrap();
    let winding: Vec<_> = conductors(&card, 300.0)
        .into_iter()
        .filter(|c| [0, 1, 4, 5].contains(&c.triangle))
        .collect();
    let original = NonlinearWindingCircuit::new(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 1.0 }),
                (2, 0, Branch::Resistor { ohms: 0.3 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        winding.clone(),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap();
    // Each isolated core half has C_aa=sigma*d/36. Its weak field equation
    // is I=H(2a)+sigma*delta(a)/(72dt), independently of the terminal equation.
    let a = bisect(0.0, 2.0, |a| {
        4.0 * depth * a / dt + (2.4 + 0.3) * (h(2.0 * a, true) + sigma * a / (72.0 * dt)) - 9.0
    });
    let i = h(2.0 * a, true) + sigma * a / (72.0 * dt);
    let qcore = sigma * depth * a * a / (18.0 * dt);
    let candidate = original
        .step_to_candidate(dt, core, winding, &[9.0], controls(), &mut || false)
        .unwrap();
    near(candidate.circuit().currents_a()[0], i);
    near(candidate.circuit().field().potential_wb_m()[6], a);
    near(candidate.step().field.heat_j, qcore);
    near(
        candidate.step().cells.iter().map(|c| c.heat_j).sum(),
        2.4 * i * i * dt,
    );
    near(candidate.step().heat_j, qcore + 2.7 * i * i * dt);
    assert!(qcore > 0.0 && candidate.step().field.numerical_dissipation_j > 0.0);
    for core in &candidate.step().field.cells {
        assert!(
            !candidate
                .step()
                .cells
                .iter()
                .any(|w| w.source.triangle == core.source.triangle)
        );
    }
    assert!(candidate.step().energy_residual_j.abs() < controls().energy_tolerance_j);
}

#[test]
fn source_refresh_refusals_work_budget_and_late_cancellation_preserve_both_histories() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let original = common::single(&problem, &card, 0.3, 0.0, 0.0);
    let old_state = original.state().to_vec();
    let old_field = original.field().potential_wb_m().to_vec();
    let warm = original
        .step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 400.0),
            &[9.0],
            controls(),
            &mut || false,
        )
        .unwrap();
    let expected = bisect(0.0, 13.0, |i| {
        (1.8 + 0.3) * i + linkage(i, 0.3, true) / 0.2 - 9.0
    });
    near(warm.circuit().currents_a()[0], expected);
    near(warm.circuit().resistance_ohm()[0], 1.8);
    let cold = original
        .step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[9.0],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(warm.circuit().currents_a()[0] < cold.circuit().currents_a()[0]);

    let other = common::source(false, 18.0);
    assert!(
        original
            .step_to_candidate(
                0.2,
                vec![],
                conductors(&other, 300.0),
                &[9.0],
                controls(),
                &mut || false
            )
            .is_err()
    );
    let mut wrong_basis = conductors(&card, 300.0);
    wrong_basis[0].conductor =
        common::material::conductor(&card, 300.0, ElectricalPropertyBasis::Conductivity);
    assert!(
        original
            .step_to_candidate(0.2, vec![], wrong_basis, &[9.0], controls(), &mut || false)
            .is_err()
    );
    let mut limited = controls();
    limited.maximum_field_trials = 1;
    assert!(matches!(
        original.step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[9.0],
            limited,
            &mut || false
        ),
        Err(NonlinearWindingError::TrialBudget)
    ));
    let mut no_field = controls();
    no_field.field.nonlinear.max_newton_iterations = 0;
    assert!(
        original
            .step_to_candidate(
                0.2,
                vec![],
                conductors(&card, 300.0),
                &[9.0],
                no_field,
                &mut || false
            )
            .is_err()
    );

    let mut checkpoints = 0;
    let success = original
        .step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[9.0],
            controls(),
            &mut || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    let cancelled = original.step_to_candidate(
        0.2,
        vec![],
        conductors(&card, 300.0),
        &[9.0],
        controls(),
        &mut || {
            calls += 1;
            calls >= checkpoints
        },
    );
    assert!(matches!(cancelled, Err(NonlinearWindingError::Cancelled)));
    assert_eq!(original.state(), old_state);
    assert_eq!(original.field().potential_wb_m(), old_field);
    near(original.time_s(), 0.0);
    near(original.currents_a()[0], 0.0);
    let retry = original
        .step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[9.0],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(retry.circuit().state(), success.circuit().state());
    assert_eq!(
        retry.circuit().field().potential_wb_m(),
        success.circuit().field().potential_wb_m()
    );
    assert_eq!(retry.circuit().currents_a(), success.circuit().currents_a());
}
