//! Actuator supply/recovery, positive heat and original finite caloric feedback.
use super::*;
use fs_couple::rigid_damping_heat::driven::{
    DrivenDampingImplicitConfig, DrivenDampingThermalConfig, HeldInternalEffort,
    advance_candidate as advance_driven, advance_implicit_candidate as advance_driven_implicit,
    advance_temperature_candidate as advance_driven_temperature,
    advance_temperature_implicit_candidate as advance_driven_temperature_implicit,
};

fn pair_input<'a>(
    bodies: &'a Bodies,
    routes: &'a [BTreeMap<CellId, f64>],
) -> RigidDampingThermalInput<'a> {
    let mut input = bodies.input(routes);
    input.states = &bodies.states[..2];
    input.properties = &bodies.properties[..2];
    input.links = &bodies.links[..1];
    input.loads.wrenches = &bodies.loads[..2];
    input
}
fn controls() -> DrivenDampingThermalConfig {
    let baseline = config();
    DrivenDampingThermalConfig {
        mechanical: baseline.mechanical.source,
        thermal: baseline.thermal,
        maximum_heat_bindings: baseline.maximum_heat_bindings,
        energy_tolerance_j: baseline.energy_tolerance_j,
    }
}
fn effort(force: f64) -> HeldInternalEffort {
    HeldInternalEffort::MaterialFirst([force, 0.0, 0.0, 0.0, 0.0, 0.0])
}
fn scalar_reference(c: f64, force: f64, time: f64) -> (f64, f64, f64) {
    // Two original 2-kg bodies give relative mobility 1/kg and reduced
    // mass 1 kg. u = v_B-v_A obeys u'=F-c*u, u(0)=-2 m/s.
    let terminal = force / c;
    let transient = -2.0 - terminal;
    let first = -(-c * time).exp_m1() / c;
    let second = -(-2.0 * c * time).exp_m1() / (2.0 * c);
    let displacement = terminal * time + transient * first;
    let heat = c
        * (terminal * terminal * time
            + 2.0 * terminal * transient * first
            + transient * transient * second);
    (
        terminal + transient * (-c * time).exp(),
        force * displacement,
        heat,
    )
}

#[test]
fn implicit_driven_pair_keeps_signed_work_separate_from_heat_and_feedback() {
    let bodies = Bodies::new([0.7, 0.0]);
    let calorics = linear_calorics();
    let initial = stiff_thermal(&calorics, 8);
    let graph = stiff_graph(&calorics, &initial);
    let routes = [BTreeMap::from([(CellId(0), 0.25), (CellId(1), 0.5)])];
    let input = pair_input(&bodies, &routes);
    let observer = Observe::default();
    for scheme in [
        ImplicitThermalScheme::BackwardEuler,
        ImplicitThermalScheme::Sdirk2,
    ] {
        let common = implicit_config(scheme);
        let implicit = DrivenDampingImplicitConfig {
            mechanical: controls().mechanical,
            thermal: common.thermal,
            scheme,
            maximum_heat_bindings: common.maximum_heat_bindings,
            energy_tolerance_j: common.energy_tolerance_j,
        };
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for force in [-1.4, 1.4] {
            let (out, native) = advance_driven_implicit(
                input,
                effort(force),
                &initial,
                &observer,
                0.5,
                &graph,
                implicit,
                &mut workspace,
                || false,
            )
            .unwrap();
            let (velocity, work, heat) = scalar_reference(0.7, force, 0.5);
            let actual = out.states()[1]
                .center_of_mass_velocity_world(bodies.properties[1])
                .unwrap()
                .x
                - out.states()[0]
                    .center_of_mass_velocity_world(bodies.properties[0])
                    .unwrap()
                    .x;
            near(actual, velocity, 3e-12);
            near(out.report().actuator_work_j, work, 3e-12);
            near(out.report().dissipated_energy_j, heat, 3e-12);
            near(out.report().deposited_heat_j, 0.75 * heat, 3e-12);
            near(out.report().unallocated_energy_j, 0.25 * heat, 3e-12);
            near(native.energy.stored_energy_change_j, 0.75 * heat, 1e-9);
            near(
                out.report().mechanical_energy_change_j,
                0.5 * (velocity * velocity - 4.0),
                3e-12,
            );
            near(out.report().combined_energy_residual_j, 0.0, 1e-8);
            assert_eq!(out.report().actuator_work_j.is_sign_negative(), force > 0.0);
            assert!(heat > 0.0);
            assert_eq!(native.energy.substeps, 2);
            assert_eq!(native.energy, out.report().thermal);
            assert_eq!(out.observer().intervals.len(), 2);
        }

        let law = LinearTemperatureDamping {
            beta_per_k: 0.5,
            maximum_k: 303.0,
        };
        let (feedback, native) = advance_driven_temperature_implicit(
            input,
            effort(0.4),
            &initial,
            &observer,
            &law,
            &[CellId(0)],
            0.5,
            &graph,
            implicit,
            &mut workspace,
            || false,
        )
        .unwrap();
        let temperature = feedback.thermal().temperature_k(CellId(0)).unwrap();
        assert!(temperature > 300.0);
        assert_eq!(feedback.temperatures_k(), [temperature]);
        near(
            feedback.endpoint_links()[0].damping.native_coefficients()[0][0],
            0.7 * (1.0 + 0.5 * (temperature - 300.0)),
            2e-14,
        );
        assert!(feedback.report().actuator_work_j < 0.0);
        assert!(feedback.report().deposited_heat_j > 0.0);
        near(feedback.report().combined_energy_residual_j, 0.0, 1e-8);
        assert_eq!(native.energy, feedback.report().thermal);
    }
    assert_eq!(initial.time_s(), 0.0);
    assert!(
        initial
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
    assert_eq!(observer, Observe::default());
}

#[test]
fn g1_driven_pair_deposits_only_exact_damping_heat_with_supply_and_recovery() {
    let bodies = Bodies::new([0.7, 0.0]);
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let routes = [BTreeMap::from([(CellId(0), 0.25), (CellId(1), 0.5)])];
    for (force, time) in [(-1.4, 0.5), (1.4, 0.2)] {
        let candidate = advance_driven(
            pair_input(&bodies, &routes),
            effort(force),
            &initial,
            &(),
            time,
            controls(),
            || false,
        )
        .unwrap();
        let (velocity, work, heat) = scalar_reference(0.7, force, time);
        let actual = candidate.states()[1]
            .center_of_mass_velocity_world(bodies.properties[1])
            .unwrap()
            .x
            - candidate.states()[0]
                .center_of_mass_velocity_world(bodies.properties[0])
                .unwrap()
                .x;
        near(actual, velocity, 2e-12);
        near(candidate.report().actuator_work_j, work, 2e-12);
        near(candidate.report().dissipated_energy_j, heat, 2e-12);
        near(candidate.report().mechanical_work_j, 0.0, 0.0);
        near(
            candidate.report().mechanical_energy_change_j,
            0.5 * (velocity * velocity - 4.0),
            3e-12,
        );
        near(candidate.report().deposited_heat_j, 0.75 * heat, 2e-12);
        near(candidate.report().unallocated_energy_j, 0.25 * heat, 1e-12);
        near(candidate.report().external_thermal_heat_j, 0.0, 1e-13);
        near(candidate.report().combined_energy_residual_j, 0.0, 5e-12);
        assert!(heat > 0.0);
        if force > 0.0 {
            assert!(candidate.report().actuator_work_j < 0.0);
        }
        for (i, fraction) in [(0, 0.25), (1, 0.5)] {
            let id = CellId(i as u64);
            let specific = fraction * heat / initial.cells()[&id].mass_kg;
            near(
                candidate.thermal().cells()[&id].specific_h_j_kg,
                specific,
                2e-12,
            );
            near(
                candidate.thermal().temperature_k(id).unwrap(),
                calorics[i].temperature_k(specific).unwrap(),
                2e-12,
            );
        }
    }
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g0_drive_without_damping_supplies_work_without_thermal_deposition() {
    let bodies = Bodies::new([0.0, 0.0]);
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let time = 0.5;
    let force = 0.4;
    let candidate = advance_driven(
        pair_input(&bodies, &routes),
        effort(force),
        &initial,
        &(),
        time,
        controls(),
        || false,
    )
    .unwrap();
    let velocity = -2.0 + force * time;
    near(
        candidate.report().actuator_work_j,
        force * (-2.0 * time + 0.5 * force * time * time),
        2e-12,
    );
    near(
        candidate.report().mechanical_energy_change_j,
        0.5 * (velocity * velocity - 4.0),
        2e-12,
    );
    assert_eq!(candidate.report().dissipated_energy_j, 0.0);
    assert_eq!(candidate.report().deposited_heat_j, 0.0);
    assert_eq!(candidate.thermal().cells(), initial.cells());
}

#[test]
fn g0_native_events_and_observer_sources_share_the_driven_heat_clock() {
    let bodies = Bodies::new([0.7, 0.0]);
    let calorics = laws();
    let initial = thermal(&calorics, true, true);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let observer = Observe {
        supplied: true,
        ..Default::default()
    };
    let mut config = controls();
    config.thermal.max_dt_s = 0.6;
    let candidate = advance_driven(
        pair_input(&bodies, &routes),
        effort(0.35),
        &initial,
        &observer,
        1.0,
        config,
        || false,
    )
    .unwrap();
    assert_eq!(candidate.observer().events, [0.35]);
    assert_eq!(
        candidate.observer().intervals,
        [(0.0, 0.35, 0.35), (0.35, 0.35, 0.7), (0.7, 1.0 - 0.7, 1.0)]
    );
    let (_, work, heat) = scalar_reference(0.7, 0.35, 1.0);
    near(candidate.report().actuator_work_j, work, 3e-12);
    near(candidate.report().dissipated_energy_j, heat, 3e-12);
    near(
        candidate.report().external_thermal_heat_j,
        2.0 * 0.35 + 4.0 * 0.65 + 1.0 + 0.25,
        3e-12,
    );
    near(candidate.report().combined_energy_residual_j, 0.0, 5e-12);
    assert_eq!(candidate.report().pair_relaxations, 6);
    assert!(observer.intervals.is_empty());
}

fn original_feedback_reference(time: f64, force: f64, beta: f64) -> [f64; 2] {
    // Independent continuous original mechanical/caloric ODE. All heat goes
    // to cell0, whose original m*Cp = 1 J/K, hence T-300 = Q.
    let derivative = |state: [f64; 2]| {
        let c = 0.7 * (1.0 + beta * state[1]);
        [force - c * state[0], c * state[0] * state[0]]
    };
    let add = |a: [f64; 2], b: [f64; 2], h: f64| [a[0] + h * b[0], a[1] + h * b[1]];
    let mut state = [-2.0, 0.0];
    let steps = 32768;
    let h = time / steps as f64;
    for _ in 0..steps {
        let a = derivative(state);
        let b = derivative(add(state, a, 0.5 * h));
        let c = derivative(add(state, b, 0.5 * h));
        let d = derivative(add(state, c, h));
        for i in 0..2 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}
#[test]
fn g1_driven_temperature_feedback_refines_to_original_nonlinear_ode() {
    let bodies = Bodies::new([0.7, 0.0]);
    let calorics = linear_calorics();
    let initial = thermal(&calorics, false, true);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let law = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    let force = 0.6;
    let target = 0.5;
    let exact = original_feedback_reference(target, force, law.beta_per_k);
    let mut errors = Vec::new();
    for steps in [16, 32, 64, 128] {
        let mut config = controls();
        config.thermal.max_dt_s = target / steps as f64;
        let candidate = advance_driven_temperature(
            pair_input(&bodies, &routes),
            effort(force),
            &initial,
            &(),
            &law,
            &[CellId(0)],
            target,
            config,
            || false,
        )
        .unwrap();
        let velocity = candidate.states()[1]
            .center_of_mass_velocity_world(bodies.properties[1])
            .unwrap()
            .x
            - candidate.states()[0]
                .center_of_mass_velocity_world(bodies.properties[0])
                .unwrap()
                .x;
        let heat = candidate.report().dissipated_energy_j;
        let temperature = candidate.thermal().temperature_k(CellId(0)).unwrap();
        near(temperature - 300.0, heat, 3e-12);
        near(
            0.5 * (velocity * velocity - 4.0) + heat - candidate.report().actuator_work_j,
            0.0,
            1e-11,
        );
        near(candidate.report().combined_energy_residual_j, 0.0, 1e-11);
        assert_eq!(candidate.temperatures_k(), &[temperature]);
        near(
            candidate.endpoint_links()[0].damping.native_coefficients()[0][0],
            0.7 * (1.0 + law.beta_per_k * (temperature - 300.0)),
            0.0,
        );
        errors.push((velocity - exact[0]).abs() + (heat - exact[1]).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            (1.8..2.2).contains(&(pair[0] / pair[1])),
            "first-order feedback: {errors:?}"
        );
    }
}

#[test]
fn g4_driven_source_observer_and_cancellation_refusals_are_transactional() {
    let bodies = Bodies::new([0.7, 0.0]);
    let calorics = linear_calorics();
    let initial = thermal(&calorics, false, true);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let observer = Observe {
        fail_at: Some(2),
        ..Default::default()
    };
    assert!(matches!(
        advance_driven(
            pair_input(&bodies, &routes),
            effort(0.6),
            &initial,
            &observer,
            0.25,
            controls(),
            || false,
        ),
        Err(RigidDampingThermalError::Observer(_))
    ));
    let law = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 300.01,
    };
    assert!(matches!(
        advance_driven_temperature(
            pair_input(&bodies, &routes),
            effort(0.6),
            &initial,
            &(),
            &law,
            &[CellId(0)],
            0.125,
            controls(),
            || false,
        ),
        Err(RigidDampingThermalError::Source { connector: 0, .. })
    ));
    let mut polls = 0;
    let expected = advance_driven(
        pair_input(&bodies, &routes),
        effort(0.6),
        &initial,
        &(),
        0.25,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut repeated = 0;
    assert!(
        advance_driven(
            pair_input(&bodies, &routes),
            effort(0.6),
            &initial,
            &(),
            0.25,
            controls(),
            || {
                repeated += 1;
                repeated >= polls - 1
            },
        )
        .is_err()
    );
    let retry = advance_driven(
        pair_input(&bodies, &routes),
        effort(0.6),
        &initial,
        &(),
        0.25,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(retry.states(), expected.states());
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.intervals.is_empty());
}

#[test]
fn g0_driven_pair_requires_original_canonical_topology_and_temperature_coordinate() {
    let bodies = Bodies::new([0.7, 0.0]);
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    assert!(
        advance_driven(
            bodies.input(&routes),
            effort(0.6),
            &initial,
            &(),
            0.125,
            controls(),
            || false
        )
        .is_err()
    );
    let law = LinearTemperatureDamping {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    assert!(
        advance_driven_temperature(
            pair_input(&bodies, &routes),
            effort(0.6),
            &initial,
            &(),
            &law,
            &[],
            0.125,
            controls(),
            || false,
        )
        .is_err()
    );
    let mut link = bodies.links[0].clone();
    link.bodies = [1, 0];
    let links = [link];
    let mut input = pair_input(&bodies, &routes);
    input.links = &links;
    assert!(
        advance_driven(input, effort(0.6), &initial, &(), 0.125, controls(), || {
            false
        })
        .is_err()
    );
}
