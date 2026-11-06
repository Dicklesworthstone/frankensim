//! End-region mutual flux must enter native nonlinear terminal equations,
//! capacitor storage and physical loss accounting through the same factors.
use super::{bisect, common};
use common::{b_from_current, conductors, controls, density_energy, linkage, near};
use fs_em::winding::{
    end::EndWindingMode,
    nonlinear::{
        NonlinearWindingCircuit, NonlinearWindingInitial, NonlinearWindingPort, NonlinearWindings,
        resistors::GraphResistorMaterial,
    },
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::{Branch, CircuitGraph};

#[test]
fn mutual_ends_and_native_capacitor_follow_independent_nonlinear_loop_equations() {
    let (depth, dt, vs, rp, rs, c, q0, rw) = (0.3, 0.2, 9.0, 0.3, 0.7, 0.4, 0.8, 2.4);
    // Diagonal [.1,.2], plus two independent physical common-flux paths.
    let ends = [0.33, -0.26, -0.26, 0.82];
    for saturating in [false, true] {
        let card = common::source(saturating, 18.0);
        let problem = common::problem(&card, depth);
        let profiles = common::profiles(2);
        let original = NonlinearWindingCircuit::new_with_end_modes(
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
                end_inductance_h: vec![0.1, 0.2],
            },
            vec![NonlinearWindingPort::CurrentInput; 2],
            vec![
                EndWindingMode {
                    inductance_h: 0.15,
                    turns: vec![1.0, -2.0],
                },
                EndWindingMode {
                    inductance_h: 0.02,
                    turns: vec![2.0, 1.0],
                },
            ],
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
        for (&actual, expected) in original.end_inductance_matrix_h().iter().zip(ends) {
            near(actual, expected);
        }
        // H(B)=I1+I2 on this exact coarse-mesh mode. The native capacitor
        // remains trapezoidal; only the magnetic field/end law uses BE.
        let a = rw + rp + ends[0] / dt;
        let d = rw + rs + dt / (2.0 * c) + ends[3] / dt;
        let cross = ends[1] / dt;
        let determinant = a * d - cross * cross;
        let currents_at = |s| {
            let g = linkage(s, depth, saturating) / dt;
            [
                (d * (vs - g) - cross * (q0 / c - g)) / determinant,
                (a * (q0 / c - g) - cross * (vs - g)) / determinant,
            ]
        };
        let total = bisect(-4.0, 4.0, |s| {
            let i = currents_at(s);
            s - i[0] - i[1]
        });
        let i = currents_at(total);
        let q1 = q0 - dt * i[1];
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
        for w in 0..2 {
            near(next.currents_a()[w], i[w]);
            near(
                next.flux_linkages_wb()[w],
                linkage(total, depth, saturating) + ends[2 * w] * i[0] + ends[2 * w + 1] * i[1],
            );
        }
        near(next.state()[next.dae().charge_index[0]], q1);
        near(step.terminal_voltage_v[0], vs - rp * i[0]);
        near(
            step.terminal_voltage_v[1],
            (q0 + q1) / (2.0 * c) - rs * i[1],
        );
        let magnetic = 2.0 * depth * density_energy(b_from_current(total, saturating), saturating);
        let end_energy = 0.05 * i[0] * i[0]
            + 0.1 * i[1] * i[1]
            + 0.075 * (i[0] - 2.0 * i[1]).powi(2)
            + 0.01 * (2.0 * i[0] + i[1]).powi(2);
        let capacitor_change = (q1 * q1 - q0 * q0) / (2.0 * c);
        let heat = dt * ((rw + rp) * i[0] * i[0] + (rw + rs) * i[1] * i[1]);
        let numerical = total * linkage(total, depth, saturating) - magnetic + end_energy;
        near(step.stored_before_j, q0 * q0 / (2.0 * c));
        near(
            step.stored_after_j,
            magnetic + end_energy + q1 * q1 / (2.0 * c),
        );
        near(
            step.stored_change_j,
            magnetic + end_energy + capacitor_change,
        );
        near(step.native.delta_h, capacitor_change);
        near(step.heat_j, heat);
        near(step.numerical_dissipation_j, numerical);
        near(step.supplied_work_j, vs * i[0] * dt);
        near(step.stored_change_j + heat + numerical, vs * i[0] * dt);
        near(
            next.field().potential_wb_m()[6],
            b_from_current(total, saturating) / 2.0,
        );
        assert_eq!(next.end_inductance(), original.end_inductance());
        assert_eq!(original.currents_a(), [0.0, 0.0]);
        assert_eq!(original.time_s(), 0.0);
        assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
}

#[test]
fn reversed_end_coupling_retains_live_history_through_hot_resistor_and_winding_refresh() {
    let (depth, dt, secondary_v) = (0.3, 0.2, 1.0);
    let ends = [0.25, -0.3, -0.3, 0.8];
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, depth);
    let mut profiles = common::profiles(2);
    for p in &mut profiles[1] {
        *p = -*p;
    }
    let original = NonlinearWindingCircuit::new_with_end_modes(
        common::initial_field(&problem, &profiles, &[0.0, 0.0]),
        CircuitGraph {
            node_count: 4,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 777.0 }),
                (2, 0, Branch::VoltageSource { port: 1 }),
                (3, 2, Branch::Inductor { henries: 999.0 }),
                (3, 0, Branch::Resistor { ohms: 0.7 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1, 3],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.1, 0.2],
        },
        vec![
            NonlinearWindingPort::VoltageInput {
                initial_voltage_v: 0.0,
                residual_scale_ohm: 2.0,
                current_tolerance_a: 1e-8,
            },
            NonlinearWindingPort::CurrentInput,
        ],
        vec![EndWindingMode {
            inductance_h: 0.15,
            turns: vec![1.0, -2.0],
        }],
        conductors(&card, 300.0),
        NonlinearWindingInitial {
            currents_a: vec![0.0, 0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0, 0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap();
    let end_energy =
        |i: [f64; 2]| 0.05 * i[0] * i[0] + 0.1 * i[1] * i[1] + 0.075 * (i[0] - 2.0 * i[1]).powi(2);
    let mut accepted = original.clone();
    let mut old = [0.0, 0.0];
    for (temperature, primary_i) in [(300.0, 3.0), (400.0, 2.0)] {
        let resistor =
            common::material::conductor(&card, temperature, ElectricalPropertyBasis::Resistivity)
                .into_uniform_conductor(0.7 * 18.0, 1.0)
                .unwrap();
        let refreshed = accepted
            .with_graph_resistors(
                vec![GraphResistorMaterial {
                    branch: 4,
                    conductor: resistor,
                }],
                &mut || false,
            )
            .unwrap();
        assert_eq!(refreshed.end_inductance(), original.end_inductance());
        assert_eq!(refreshed.state(), accepted.state());
        assert_eq!(
            refreshed.field().potential_wb_m(),
            accepted.field().potential_wb_m()
        );
        assert_eq!(refreshed.currents_a(), accepted.currents_a());
        // The source's resistivity is affine between 300 K and 500 K.
        let ratio = 1.0 + (temperature - 300.0) / 200.0;
        let (rw, rs) = (2.4 * ratio, 0.7 * ratio);
        let old_total = old[0] - old[1];
        let g = |s| (linkage(s, depth, true) - linkage(old_total, depth, true)) / dt;
        let secondary = |s| {
            (g(s) - secondary_v - ends[1] * (primary_i - old[0]) / dt + ends[3] * old[1] / dt)
                / (rw + rs + ends[3] / dt)
        };
        let total = bisect(-4.0, 5.0, |s| s - primary_i + secondary(s));
        let i = [primary_i, secondary(total)];
        let di = [i[0] - old[0], i[1] - old[1]];
        let primary_v = g(total) + (ends[0] * di[0] + ends[1] * di[1]) / dt + rw * i[0];
        let candidate = refreshed
            .step_to_candidate(
                accepted.time_s() + dt,
                vec![],
                conductors(&card, temperature),
                &[primary_i, secondary_v],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        let magnetic = 2.0 * depth * density_energy(b_from_current(total, true), true);
        let magnetic_old = 2.0 * depth * density_energy(b_from_current(old_total, true), true);
        let numerical = total * g(total) * dt - (magnetic - magnetic_old) + end_energy(di);
        let heat = dt * (rw * i[0] * i[0] + (rw + rs) * i[1] * i[1]);
        let work = dt * (i[0] * primary_v - secondary_v * i[1]);
        for w in 0..2 {
            near(candidate.circuit().currents_a()[w], i[w]);
            near(candidate.circuit().resistance_ohm()[w], rw);
            near(
                candidate.circuit().flux_linkages_wb()[w],
                [1.0, -1.0][w] * linkage(total, depth, true)
                    + ends[2 * w] * i[0]
                    + ends[2 * w + 1] * i[1],
            );
        }
        near(step.terminal_voltage_v[0], primary_v);
        near(step.terminal_voltage_v[1], -secondary_v - rs * i[1]);
        near(step.stored_before_j, magnetic_old + end_energy(old));
        near(step.stored_after_j, magnetic + end_energy(i));
        near(
            step.stored_change_j,
            magnetic - magnetic_old + end_energy(i) - end_energy(old),
        );
        near(step.numerical_dissipation_j, numerical);
        near(step.heat_j, heat);
        near(step.graph_resistors[0].heat_j, dt * rs * i[1] * i[1]);
        near(step.supplied_work_j, work);
        near(step.stored_change_j + heat + numerical, work);
        assert_eq!(
            candidate.circuit().end_inductance(),
            original.end_inductance()
        );
        assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
        old = i;
        accepted = candidate.into_parts().0;
    }
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.currents_a(), [0.0, 0.0]);
}
