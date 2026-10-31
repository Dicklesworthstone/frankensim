//! Native current drive must act through actual nonlinear linkage and storage.
//! All material laws and exact coarse-mesh oracles here are manufactured.
use super::{bisect, common};
use common::{b_from_current, conductors, controls, density_energy, linkage, near};
use fs_em::{
    planar::{
        eddy::{EddyConductor, nonlinear::PlanarNonlinearEddy},
        nonlinear::PlanarNonlinearMagnetostatics,
    },
    winding::nonlinear::{
        NonlinearWindingCircuit, NonlinearWindingError, NonlinearWindingInitial,
        NonlinearWindingPort, NonlinearWindings, resistors::GraphResistorMaterial,
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::{Branch, CircuitGraph};

fn voltage_port(initial_voltage_v: f64, residual_scale_ohm: f64) -> NonlinearWindingPort {
    NonlinearWindingPort::VoltageInput {
        initial_voltage_v,
        residual_scale_ohm,
        current_tolerance_a: 1e-8,
    }
}

fn driven<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
    initial_current: f64,
    initial_source: f64,
    end: f64,
    port: NonlinearWindingPort,
) -> Result<NonlinearWindingCircuit<'m>, NonlinearWindingError> {
    let profiles = common::profiles(1);
    NonlinearWindingCircuit::new_with_ports(
        common::initial_field(problem, &profiles, &[initial_current]),
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
            profiles_per_m2: profiles,
            end_inductance_h: vec![end],
        },
        vec![port],
        conductors(card, 300.0),
        NonlinearWindingInitial {
            currents_a: vec![initial_current],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![initial_source],
        },
        common::config(),
        &mut || false,
    )
}

#[test]
fn ideal_current_drive_resolves_signed_nonlinear_voltage_end_storage_and_history() {
    let depth = 0.3;
    let end = 0.2;
    let rw = 1.2;
    for saturating in [false, true] {
        let card = common::source(saturating, 18.0);
        let problem = common::problem(&card, depth);
        for sign in [-1.0, 1.0] {
            // Changing only the numerical residual scale must not manufacture
            // a resistance or alter any terminal voltage, heat or stored energy.
            for scale in [0.2, 5.0] {
                let i0 = sign * 0.4;
                let original =
                    driven(&problem, &card, i0, i0, end, voltage_port(rw * i0, scale)).unwrap();
                assert!(original.dae().flux_index.is_empty());
                assert!(original.dae().charge_index.is_empty());
                let mut accepted = original.clone();
                for (target, current) in [(0.2, sign * 3.0), (0.3, sign * 1.5)] {
                    let dt = target - accepted.time_s();
                    let previous = accepted.currents_a()[0];
                    let before_b = b_from_current(previous, saturating);
                    let after_b = b_from_current(current, saturating);
                    let field_before = 2.0 * depth * density_energy(before_b, saturating);
                    let field_after = 2.0 * depth * density_energy(after_b, saturating);
                    let delta_linkage =
                        linkage(current, depth, saturating) - linkage(previous, depth, saturating);
                    let voltage = (delta_linkage + end * (current - previous)) / dt + rw * current;
                    let numerical = current * delta_linkage - (field_after - field_before)
                        + 0.5 * end * (current - previous).powi(2);
                    let candidate = accepted
                        .step_to_candidate(
                            target,
                            vec![],
                            conductors(&card, 300.0),
                            &[current],
                            controls(),
                            &mut || false,
                        )
                        .unwrap();
                    let step = candidate.step();
                    near(candidate.circuit().currents_a()[0], current);
                    near(
                        candidate.circuit().field().potential_wb_m()[6],
                        after_b / 2.0,
                    );
                    near(step.terminal_voltage_v[0], voltage);
                    near(step.native.y[0], voltage);
                    near(step.native.y[1], -current);
                    near(step.native.supplied, 0.0);
                    near(step.native.dissipated, 0.0);
                    near(
                        step.stored_before_j,
                        field_before + 0.5 * end * previous.powi(2),
                    );
                    near(
                        step.stored_after_j,
                        field_after + 0.5 * end * current.powi(2),
                    );
                    near(step.heat_j, dt * rw * current.powi(2));
                    near(step.numerical_dissipation_j, numerical);
                    near(step.supplied_work_j, dt * current * voltage);
                    near(
                        step.flux_linkages_wb[0][0],
                        linkage(previous, depth, saturating) + end * previous,
                    );
                    near(
                        step.flux_linkages_wb[1][0],
                        linkage(current, depth, saturating) + end * current,
                    );
                    assert!(step.terminal_current_residual_a[0].abs() <= 1e-8);
                    assert!(
                        step.terminal_residual_v[0].abs()
                            <= controls().terminal_voltage_tolerance_v
                    );
                    assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
                    accepted = candidate.into_parts().0;
                }
                near(original.currents_a()[0], i0);
                near(original.time_s(), 0.0);
                near(
                    original.flux_linkages_wb()[0],
                    linkage(i0, depth, saturating) + end * i0,
                );
            }
        }
    }
}

#[test]
fn smooth_current_ramp_converges_to_physical_work_and_removes_be_loss_under_refinement() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let mut losses = Vec::new();
    let mut errors = Vec::new();
    for count in [6, 12] {
        let mut accepted = driven(&problem, &card, 0.0, 0.0, 0.2, voltage_port(0.0, 1.0)).unwrap();
        let mut supplied = 0.0;
        let mut heat = 0.0;
        let mut numerical = 0.0;
        for index in 1..=count {
            // I(t)=3*t over one second. Both partitions put a node at the
            // manufactured saturation corner I=1, so each interval has exact
            // piecewise-affine linkage and an independent BE-loss oracle.
            let target = f64::from(index) / f64::from(count);
            let current = 3.0 * target;
            let candidate = accepted
                .step_to_candidate(
                    target,
                    vec![],
                    conductors(&card, 300.0),
                    &[current],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            supplied += candidate.step().supplied_work_j;
            heat += candidate.step().heat_j;
            numerical += candidate.step().numerical_dissipation_j;
            accepted = candidate.into_parts().0;
        }
        let n = f64::from(count);
        near(
            heat,
            1.2 * 9.0 * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n),
        );
        near(numerical, 2.25 / n);
        // Exact continuous ramp: final magnetic+end energy is 1.8 J, and
        // integral(R*I(t)^2 dt)=3.6 J. Temporal BE loss is not physical heat.
        near(supplied, 1.8 + heat + numerical);
        losses.push(numerical);
        errors.push((supplied - 5.4).abs());
    }
    near(losses[0], 2.0 * losses[1]);
    assert!(errors[1] < 0.51 * errors[0]);
}

#[test]
fn voltage_port_preserves_native_capacitor_charge_and_unbound_inductor_flux() {
    let depth = 0.3;
    let dt = 0.2;
    let c = 0.4;
    let l = 0.9;
    let source = 30.0;
    let rw = 1.2;
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, depth);
    let profiles = common::profiles(1);
    let original = NonlinearWindingCircuit::new_with_ports(
        common::initial_field(&problem, &profiles, &[0.0]),
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 777.0 }),
                (1, 0, Branch::Capacitor { farads: c }),
                (1, 0, Branch::Inductor { henries: l }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        vec![voltage_port(0.0, 1.0)],
        conductors(&card, 300.0),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![0.0],
            charges_c: vec![0.0],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap();
    // Native midpoint storage gives q1=2*C*V, phi1=dt*V. Kirchhoff therefore
    // requires I+(2*C/dt+dt/(2*L))*V=source, independently of the field solve.
    let required_voltage = |i| linkage(i, depth, true) / dt + rw * i;
    let current = bisect(0.0, 13.0, |i| {
        i + (2.0 * c / dt + dt / (2.0 * l)) * required_voltage(i) - source
    });
    let voltage = required_voltage(current);
    let candidate = original
        .step_to_candidate(
            dt,
            vec![],
            conductors(&card, 300.0),
            &[source],
            controls(),
            &mut || false,
        )
        .unwrap();
    let next = candidate.circuit();
    let step = candidate.step();
    let q = 2.0 * c * voltage;
    let phi = dt * voltage;
    let native_storage = q * q / (2.0 * c) + phi * phi / (2.0 * l);
    let field_storage = 2.0 * depth * density_energy(b_from_current(current, true), true);
    near(next.currents_a()[0], current);
    near(step.terminal_voltage_v[0], voltage);
    near(next.state()[next.dae().charge_index[0]], q);
    near(next.state()[next.dae().flux_index[0]], phi);
    near(step.native.delta_h, native_storage);
    near(step.stored_before_j, 0.0);
    near(step.stored_after_j, native_storage + field_storage);
    near(step.heat_j, dt * rw * current.powi(2));
    near(
        step.numerical_dissipation_j,
        current * linkage(current, depth, true) - field_storage,
    );
    near(step.supplied_work_j, dt * source * voltage);
    assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
    // Reinitialization against the trial voltage would contradict this zero
    // initial charge; the finite native solve must evolve it instead.
    assert_eq!(original.state()[original.dae().charge_index[0]], 0.0);
    assert_eq!(original.state()[original.dae().flux_index[0]], 0.0);
}

#[test]
fn mixed_causality_reversed_mutual_winding_and_graph_material_refresh_are_power_reciprocal() {
    let depth = 0.3;
    let dt = 0.2;
    let current_primary = 3.0;
    let secondary_source = 1.0;
    let rw = 2.4;
    let ends = [0.1, 0.2];
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, depth);
    let mut profiles = common::profiles(2);
    for p in &mut profiles[1] {
        *p = -*p;
    }
    let original = NonlinearWindingCircuit::new_with_ports(
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
            end_inductance_h: ends.to_vec(),
        },
        vec![voltage_port(0.0, 2.0), NonlinearWindingPort::CurrentInput],
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
    for temperature in [300.0, 400.0] {
        let resistor =
            common::material::conductor(&card, temperature, ElectricalPropertyBasis::Resistivity)
                .into_uniform_conductor(0.7 * 18.0, 1.0)
                .unwrap();
        let rs = resistor.resistance_ohm();
        let refreshed = original
            .with_graph_resistors(
                vec![GraphResistorMaterial {
                    branch: 4,
                    conductor: resistor,
                }],
                &mut || false,
            )
            .unwrap();
        assert_eq!(refreshed.state(), original.state());
        assert_eq!(refreshed.terminal_ports(), original.terminal_ports());
        let a2 = rw + rs + ends[1] / dt;
        // Reversed spatial turns give H(B)=I1-I2 and lambda2=-lambda1_field.
        // Reversed second terminals give V2=-Vs-Rs*I2. Both signs must survive
        // the independent native power receipt, not just the Newton residual.
        let total = bisect(-1.0, 5.0, |s| {
            s - current_primary + (linkage(s, depth, true) / dt - secondary_source) / a2
        });
        let g = linkage(total, depth, true) / dt;
        let current_secondary = (g - secondary_source) / a2;
        let primary_voltage = g + (rw + ends[0] / dt) * current_primary;
        let secondary_voltage = -secondary_source - rs * current_secondary;
        let candidate = refreshed
            .step_to_candidate(
                dt,
                vec![],
                conductors(&card, 300.0),
                &[current_primary, secondary_source],
                controls(),
                &mut || false,
            )
            .unwrap();
        let next = candidate.circuit();
        let step = candidate.step();
        let magnetic = 2.0 * depth * density_energy(b_from_current(total, true), true);
        let end_energy =
            0.5 * (ends[0] * current_primary.powi(2) + ends[1] * current_secondary.powi(2));
        let heat = dt * (rw * current_primary.powi(2) + (rw + rs) * current_secondary.powi(2));
        let work = dt * (current_primary * primary_voltage - secondary_source * current_secondary);
        near(next.currents_a()[0], current_primary);
        near(next.currents_a()[1], current_secondary);
        near(
            next.flux_linkages_wb()[0],
            linkage(total, depth, true) + ends[0] * current_primary,
        );
        near(
            next.flux_linkages_wb()[1],
            -linkage(total, depth, true) + ends[1] * current_secondary,
        );
        near(step.terminal_voltage_v[0], primary_voltage);
        near(step.terminal_voltage_v[1], secondary_voltage);
        near(step.native.y[0], primary_voltage);
        near(step.native.y[1], -current_secondary);
        near(step.native.y[2], -current_primary);
        near(step.native.y[3], secondary_voltage);
        near(step.stored_after_j, magnetic + end_energy);
        near(step.heat_j, heat);
        near(step.supplied_work_j, work);
        near(
            step.numerical_dissipation_j,
            total * linkage(total, depth, true) - magnetic + end_energy,
        );
        near(step.native.dissipated, dt * rs * current_secondary.powi(2));
        near(
            step.graph_resistors[0].heat_j,
            dt * rs * current_secondary.powi(2),
        );
        assert_eq!(step.graph_materials.len(), 1);
        assert!(step.terminal_current_residual_a[0].abs() < 1e-8);
        assert_eq!(step.terminal_current_residual_a[1], 0.0);
        assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
    }
}

#[test]
fn current_drive_resolves_nonlinear_core_eddy_shielding_and_separate_joule_heat() {
    let depth = 0.3;
    let dt = 0.2;
    let sigma = 18.0;
    let current = 3.0;
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
    let original = NonlinearWindingCircuit::new_with_ports(
        field,
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
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        vec![voltage_port(0.0, 1.0)],
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
    let a = bisect(0.0, 2.0, |a| {
        common::h(2.0 * a, true) + sigma * a / (72.0 * dt) - current
    });
    let core_heat = sigma * depth * a * a / (18.0 * dt);
    let winding_heat = 2.4 * current * current * dt;
    let voltage = 4.0 * depth * a / dt + 2.4 * current;
    let candidate = original
        .step_to_candidate(dt, core, winding, &[current], controls(), &mut || false)
        .unwrap();
    near(candidate.circuit().currents_a()[0], current);
    near(candidate.circuit().field().potential_wb_m()[6], a);
    near(candidate.step().terminal_voltage_v[0], voltage);
    near(candidate.step().field.heat_j, core_heat);
    near(
        candidate.step().cells.iter().map(|cell| cell.heat_j).sum(),
        winding_heat,
    );
    near(candidate.step().heat_j, core_heat + winding_heat);
    near(candidate.step().supplied_work_j, dt * current * voltage);
    assert!(2.0 * a < b_from_current(current, true));
    assert!(candidate.step().numerical_dissipation_j > 0.0);
    assert!(candidate.step().energy_residual_j.abs() < controls().energy_tolerance_j);
}

#[test]
fn initial_current_checks_physical_residual_budget_and_cancellation_are_transactional() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    assert!(matches!(
        driven(&problem, &card, 0.0, 1.0, 0.0, voltage_port(0.0, 1.0)),
        Err(NonlinearWindingError::Invalid(
            "native initial winding current disagrees with accepted field history"
        ))
    ));
    for port in [
        voltage_port(f64::NAN, 1.0),
        voltage_port(0.0, 0.0),
        voltage_port(0.0, f64::INFINITY),
        NonlinearWindingPort::VoltageInput {
            initial_voltage_v: 0.0,
            residual_scale_ohm: 1.0,
            current_tolerance_a: 0.0,
        },
    ] {
        assert!(driven(&problem, &card, 0.0, 0.0, 0.0, port).is_err());
    }
    let original = driven(&problem, &card, 0.0, 0.0, 0.0, voltage_port(0.0, 1.0)).unwrap();
    let old_state = original.state().to_vec();
    let old_field = original.field().potential_wb_m().to_vec();
    let mut loose = controls();
    loose.newton.absolute_tolerance = 10.0;
    assert!(matches!(
        original.step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[0.001],
            loose,
            &mut || false
        ),
        Err(NonlinearWindingError::Balance(
            "native terminal current and actual field current do not close"
        ))
    ));
    let mut limited = controls();
    limited.maximum_field_trials = 1;
    assert!(matches!(
        original.step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[3.0],
            limited,
            &mut || false
        ),
        Err(NonlinearWindingError::TrialBudget)
    ));
    let mut checkpoints = 0;
    let success = original
        .step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[3.0],
            controls(),
            &mut || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        original.step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[3.0],
            controls(),
            &mut || {
                calls += 1;
                calls >= checkpoints
            }
        ),
        Err(NonlinearWindingError::Cancelled)
    ));
    assert_eq!(original.state(), old_state);
    assert_eq!(original.field().potential_wb_m(), old_field);
    near(original.currents_a()[0], 0.0);
    near(original.flux_linkages_wb()[0], 0.0);
    near(original.time_s(), 0.0);
    let retry = original
        .step_to_candidate(
            0.2,
            vec![],
            conductors(&card, 300.0),
            &[3.0],
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
