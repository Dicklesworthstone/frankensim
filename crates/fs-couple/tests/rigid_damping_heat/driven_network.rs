//! Shared-body cycles, counteracting actuators and finite caloric feedback.
use super::*;
use fs_couple::rigid_damping_heat::driven::network::{
    DrivenNetworkDampingImplicitConfig, DrivenNetworkDampingThermalConfig, HeldInternalEffort,
    advance_candidate as advance_graph, advance_implicit_candidate as advance_graph_implicit,
    advance_temperature_candidate as advance_graph_temperature,
    advance_temperature_implicit_candidate as advance_graph_temperature_implicit,
};
use fs_mbd::spatial_loads::{
    internal_impulse::viscous::network::driven::DrivenViscousNetworkLimits,
    internal_motion::viscous::network::driven::DrivenViscousNetworkMotionControls,
};

fn graph_controls() -> DrivenNetworkDampingThermalConfig {
    let base = config();
    DrivenNetworkDampingThermalConfig {
        mechanical: DrivenViscousNetworkMotionControls {
            source: base.mechanical.source,
            limits: DrivenViscousNetworkLimits {
                maximum_bodies: 4,
                maximum_connectors: 4,
                maximum_live_bytes: 128_000_000,
                maximum_preparation_terms: 5_000_000_000,
            },
        },
        thermal: base.thermal,
        maximum_heat_bindings: 8,
        energy_tolerance_j: base.energy_tolerance_j,
    }
}
fn held(force: f64) -> HeldInternalEffort {
    HeldInternalEffort::MaterialFirst([force, 0.0, 0.0, 0.0, 0.0, 0.0])
}
fn cycle_bodies() -> Bodies {
    let mut bodies = Bodies::new([0.7, 0.7]);
    for i in 0..3 {
        bodies.states[i] = RigidBodyState::new(
            bodies.states[i].pose(),
            Vec3::new(2.0 * (1.0 - i as f64), 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap();
    }
    bodies
}
fn cycle_links() -> [ViscousNetworkLink; 3] {
    [axial([0, 1], 0.7), axial([1, 2], 0.7), axial([0, 2], 0.7)]
}
fn cycle_routes() -> [BTreeMap<CellId, f64>; 3] {
    [
        BTreeMap::from([(CellId(0), 0.5), (CellId(1), 0.5)]),
        BTreeMap::from([(CellId(1), 1.0)]),
        BTreeMap::from([(CellId(0), 0.25)]),
    ]
}
fn graph_input<'a>(
    bodies: &'a Bodies,
    links: &'a [ViscousNetworkLink],
    routes: &'a [BTreeMap<CellId, f64>],
) -> RigidDampingThermalInput<'a> {
    RigidDampingThermalInput {
        links,
        ..bodies.input(routes)
    }
}
fn cycle_reference(force: f64, time: f64) -> (f64, f64, f64) {
    // Three original 2-kg masses with velocities [-u,0,u]. Equal
    // c=0.7 links (01,12,02), efforts (F,F,2F), give u'=1.5*(F-c*u).
    // u(0)=-1. Per-link work/heat are (1,1,4) times these integrals.
    let c = 0.7;
    let rate = 1.5 * c;
    let terminal = force / c;
    let tail = -1.0 - terminal;
    let first = -(-rate * time).exp_m1() / rate;
    let second = -(-2.0 * rate * time).exp_m1() / (2.0 * rate);
    (
        terminal + tail * (-rate * time).exp(),
        force * (terminal * time + tail * first),
        c * (terminal * terminal * time + 2.0 * terminal * tail * first + tail * tail * second),
    )
}

#[test]
fn implicit_driven_cycle_preserves_original_edge_accounts_and_temperature_bindings() {
    let bodies = cycle_bodies();
    let links = cycle_links();
    let routes = cycle_routes();
    let input = graph_input(&bodies, &links, &routes);
    let calorics = linear_calorics();
    let initial = stiff_thermal(&calorics, 8);
    let graph = stiff_graph(&calorics, &initial);
    let observer = Observe::default();
    for scheme in [
        ImplicitThermalScheme::BackwardEuler,
        ImplicitThermalScheme::Sdirk2,
    ] {
        let common = implicit_config(scheme);
        let implicit = DrivenNetworkDampingImplicitConfig {
            mechanical: graph_controls().mechanical,
            thermal: common.thermal,
            scheme,
            maximum_heat_bindings: common.maximum_heat_bindings,
            energy_tolerance_j: common.energy_tolerance_j,
        };
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        for force in [-0.35, 0.4] {
            let efforts = [held(force), held(force), held(2.0 * force)];
            let (out, native) = advance_graph_implicit(
                input,
                &efforts,
                &initial,
                &observer,
                0.5,
                &graph,
                implicit,
                &mut workspace,
                || false,
            )
            .unwrap();
            let (u, work, heat) = cycle_reference(force, 0.5);
            for (i, expected) in [(0, -u), (1, 0.0), (2, u), (3, -0.5)] {
                near(
                    out.states()[i]
                        .center_of_mass_velocity_world(bodies.properties[i])
                        .unwrap()
                        .x,
                    expected,
                    3e-11,
                );
            }
            for (i, factor) in [1.0, 1.0, 4.0].into_iter().enumerate() {
                near(out.report().connector_heat_j[i], factor * heat, 4e-11);
                near(
                    out.report().connector_actuator_work_j[i],
                    factor * work,
                    4e-11,
                );
                assert!(out.report().connector_heat_j[i] > 0.0);
            }
            near(out.report().dissipated_energy_j, 6.0 * heat, 5e-11);
            near(out.report().actuator_work_j, 6.0 * work, 5e-11);
            near(out.report().deposited_heat_j, 3.0 * heat, 5e-11);
            near(out.report().unallocated_energy_j, 3.0 * heat, 5e-11);
            near(native.energy.stored_energy_change_j, 3.0 * heat, 1e-9);
            near(out.report().combined_energy_residual_j, 0.0, 1e-8);
            assert_eq!(work.is_sign_negative(), force > 0.0);
            assert_eq!(native.energy.substeps, 2);
            assert_eq!(native.energy, out.report().thermal);
        }

        let law = LinearTemperatureDamping {
            beta_per_k: 0.5,
            maximum_k: 303.0,
        };
        let cells = [CellId(0), CellId(1), CellId(0)];
        let (out, native) = advance_graph_temperature_implicit(
            input,
            &[held(0.4), held(0.4), held(0.8)],
            &initial,
            &observer,
            &law,
            &cells,
            0.5,
            &graph,
            implicit,
            &mut workspace,
            || false,
        )
        .unwrap();
        for (i, cell) in cells.into_iter().enumerate() {
            let temperature = out.thermal().temperature_k(cell).unwrap();
            assert!(temperature > 300.0);
            assert_eq!(out.temperatures_k()[i], temperature);
            near(
                out.endpoint_links()[i].damping.native_coefficients()[0][0],
                0.7 * (1.0 + 0.5 * (temperature - 300.0)),
                2e-14,
            );
        }
        assert!(out.report().actuator_work_j < 0.0);
        assert!(out.report().deposited_heat_j > 0.0);
        near(out.report().combined_energy_residual_j, 0.0, 1e-8);
        assert_eq!(native.energy, out.report().thermal);
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
fn g1_driven_cycle_keeps_each_original_heat_and_work_with_finite_calorics() {
    let bodies = cycle_bodies();
    let links = cycle_links();
    let routes = cycle_routes();
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    for force in [-0.35, 0.4] {
        let efforts = [held(force), held(force), held(2.0 * force)];
        let target = 0.4;
        let result = advance_graph(
            graph_input(&bodies, &links, &routes),
            &efforts,
            &initial,
            &(),
            target,
            graph_controls(),
            || false,
        )
        .unwrap();
        let (u, work, heat) = cycle_reference(force, target);
        for (i, expected) in [(0, -u), (1, 0.0), (2, u), (3, -0.5)] {
            near(
                result.states()[i]
                    .center_of_mass_velocity_world(bodies.properties[i])
                    .unwrap()
                    .x,
                expected,
                2e-11,
            );
        }
        for (i, factor) in [1.0, 1.0, 4.0].into_iter().enumerate() {
            near(result.report().connector_heat_j[i], factor * heat, 2e-11);
            near(
                result.report().connector_actuator_work_j[i],
                factor * work,
                2e-11,
            );
            assert!(result.report().connector_heat_j[i] > 0.0);
        }
        near(result.report().dissipated_energy_j, 6.0 * heat, 4e-11);
        near(result.report().actuator_work_j, 6.0 * work, 4e-11);
        near(
            result.report().mechanical_energy_change_j,
            2.0 * (u * u - 1.0),
            4e-11,
        );
        near(result.report().deposited_heat_j, 3.0 * heat, 3e-11);
        near(result.report().unallocated_energy_j, 3.0 * heat, 3e-11);
        near(result.report().combined_energy_residual_j, 0.0, 5e-11);
        assert_eq!(result.report().pair_relaxations, 0);
        assert_eq!(result.report().mechanical_work_j, 0.0);
        for i in 0..2 {
            let id = CellId(i as u64);
            let specific = 1.5 * heat / initial.cells()[&id].mass_kg;
            near(
                result.thermal().cells()[&id].specific_h_j_kg,
                specific,
                2e-11,
            );
            near(
                result.thermal().temperature_k(id).unwrap(),
                calorics[i].temperature_k(specific).unwrap(),
                2e-11,
            );
        }
        assert_eq!(work.is_sign_negative(), force > 0.0);
    }
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_parallel_actuator_recovery_cancels_supply_without_erasing_damper_heat() {
    let bodies = Bodies::new([0.7, 0.0]);
    let links = [axial([0, 1], 0.7), axial([0, 1], 0.0)];
    let routes = [
        BTreeMap::from([(CellId(0), 1.0)]),
        BTreeMap::from([(CellId(1), 1.0)]),
    ];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let mut input = graph_input(&bodies, &links, &routes);
    input.states = &bodies.states[..2];
    input.properties = &bodies.properties[..2];
    input.loads.wrenches = &bodies.loads[..2];
    let target = 0.5;
    let result = advance_graph(
        input,
        &[held(0.4), held(-0.4)],
        &initial,
        &(),
        target,
        graph_controls(),
        || false,
    )
    .unwrap();
    // Opposite original efforts cancel in momentum; their separate powers
    // remain nonzero. Only the first link dissipates c*u².
    let displacement = -2.0 * -(-0.7 * target).exp_m1() / 0.7;
    let heat = 2.0 * -(-1.4 * target).exp_m1();
    near(
        result.report().connector_actuator_work_j[0],
        0.4 * displacement,
        2e-11,
    );
    near(
        result.report().connector_actuator_work_j[1],
        -0.4 * displacement,
        2e-11,
    );
    near(result.report().actuator_work_j, 0.0, 2e-11);
    near(result.report().connector_heat_j[0], heat, 2e-11);
    assert_eq!(result.report().connector_heat_j[1], 0.0);
    assert_eq!(result.thermal().cells()[&CellId(1)].specific_h_j_kg, 0.0);
    near(
        result.report().mechanical_energy_change_j + heat,
        0.0,
        4e-11,
    );
    near(result.report().deposited_heat_j, heat, 2e-11);
}

fn feedback_reference(force: f64, beta: f64, target: f64) -> [f64; 2] {
    // Independent original ODE, with all six weighted loss contributions
    // deposited into the original cell0 capacity of 1 J/K.
    let f = |s: [f64; 2]| {
        let c = 0.7 * (1.0 + beta * s[1]);
        [1.5 * (force - c * s[0]), 6.0 * c * s[0] * s[0]]
    };
    let add = |a: [f64; 2], b: [f64; 2], h: f64| [a[0] + h * b[0], a[1] + h * b[1]];
    let mut s = [-1.0, 0.0];
    let h = target / 32768.0;
    for _ in 0..32768 {
        let a = f(s);
        let b = f(add(s, a, 0.5 * h));
        let c = f(add(s, b, 0.5 * h));
        let d = f(add(s, c, h));
        for i in 0..2 {
            s[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    s
}
#[test]
fn g1_driven_cycle_temperature_feedback_refines_to_original_coupled_ode() {
    let bodies = cycle_bodies();
    let links = cycle_links();
    let routes = core::array::from_fn::<_, 3, _>(|_| BTreeMap::from([(CellId(0), 1.0)]));
    let calorics = linear_calorics();
    let initial = thermal(&calorics, false, true);
    let law = LinearTemperatureDamping {
        beta_per_k: 0.5,
        maximum_k: 303.0,
    };
    let force = 0.2;
    let target = 0.3;
    let exact = feedback_reference(force, law.beta_per_k, target);
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let mut controls = graph_controls();
        controls.thermal.max_dt_s = target / steps as f64;
        let result = advance_graph_temperature(
            graph_input(&bodies, &links, &routes),
            &[held(force), held(force), held(2.0 * force)],
            &initial,
            &(),
            &law,
            &[CellId(0); 3],
            target,
            controls,
            || false,
        )
        .unwrap();
        let u = result.states()[2]
            .center_of_mass_velocity_world(bodies.properties[2])
            .unwrap()
            .x;
        let heat = result.report().dissipated_energy_j;
        let temperature = result.thermal().temperature_k(CellId(0)).unwrap();
        near(temperature - 300.0, heat, 3e-11);
        near(
            2.0 * (u * u - 1.0) + heat - result.report().actuator_work_j,
            0.0,
            8e-11,
        );
        near(result.report().combined_energy_residual_j, 0.0, 8e-11);
        assert_eq!(result.temperatures_k(), &[temperature; 3]);
        for link in result.endpoint_links() {
            near(
                link.damping.native_coefficients()[0][0],
                0.7 * (1.0 + 0.5 * (temperature - 300.0)),
                0.0,
            );
        }
        errors.push((u - exact[0]).abs() + (heat - exact[1]).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            (1.8..2.2).contains(&(pair[0] / pair[1])),
            "first-order feedback: {errors:?}"
        );
    }
}

#[test]
fn g4_native_graph_heat_obeys_events_and_refuses_atomically() {
    let bodies = cycle_bodies();
    let links = cycle_links();
    let routes = cycle_routes();
    let calorics = laws();
    let initial = thermal(&calorics, true, true);
    let efforts = [held(0.2), held(0.2), held(0.4)];
    let observer = Observe {
        supplied: true,
        ..Default::default()
    };
    let mut controls = graph_controls();
    controls.thermal.max_dt_s = 0.6;
    let mut polls = 0;
    let expected = advance_graph(
        graph_input(&bodies, &links, &routes),
        &efforts,
        &initial,
        &observer,
        1.0,
        controls,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(expected.observer().events, [0.35]);
    assert_eq!(
        expected.observer().intervals,
        [(0.0, 0.35, 0.35), (0.35, 0.35, 0.7), (0.7, 1.0 - 0.7, 1.0)]
    );
    let (_, work, heat) = cycle_reference(0.2, 1.0);
    near(expected.report().actuator_work_j, 6.0 * work, 5e-11);
    near(expected.report().dissipated_energy_j, 6.0 * heat, 5e-11);
    near(
        expected.report().external_thermal_heat_j,
        2.0 * 0.35 + 4.0 * 0.65 + 1.0 + 0.25,
        5e-11,
    );
    let mut repeated = 0;
    assert!(matches!(
        advance_graph(
            graph_input(&bodies, &links, &routes),
            &efforts,
            &initial,
            &observer,
            1.0,
            controls,
            || {
                repeated += 1;
                repeated >= polls - 1
            },
        ),
        Err(RigidDampingThermalError::Cancelled)
    ));
    let failing = Observe {
        fail_at: Some(2),
        ..Default::default()
    };
    assert!(matches!(
        advance_graph(
            graph_input(&bodies, &links, &routes),
            &efforts,
            &initial,
            &failing,
            1.0,
            controls,
            || false
        ),
        Err(RigidDampingThermalError::Observer(_))
    ));
    assert!(
        advance_graph(
            graph_input(&bodies, &links, &routes),
            &efforts[..2],
            &initial,
            &observer,
            1.0,
            controls,
            || false
        )
        .is_err()
    );
    let limited = DrivenNetworkDampingThermalConfig {
        mechanical: DrivenViscousNetworkMotionControls {
            limits: DrivenViscousNetworkLimits {
                maximum_bodies: 3,
                ..controls.mechanical.limits
            },
            ..controls.mechanical
        },
        ..controls
    };
    assert!(
        advance_graph(
            graph_input(&bodies, &links, &routes),
            &efforts,
            &initial,
            &observer,
            1.0,
            limited,
            || false
        )
        .is_err()
    );
    let cold = thermal(&calorics, false, true);
    let law = LinearTemperatureDamping {
        beta_per_k: 0.5,
        maximum_k: 300.001,
    };
    assert!(matches!(
        advance_graph_temperature(
            graph_input(&bodies, &links, &routes),
            &efforts,
            &cold,
            &(),
            &law,
            &[CellId(0); 3],
            0.125,
            graph_controls(),
            || false,
        ),
        Err(RigidDampingThermalError::Source { .. })
    ));
    let retry = advance_graph(
        graph_input(&bodies, &links, &routes),
        &efforts,
        &initial,
        &observer,
        1.0,
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(retry.states(), expected.states());
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.intervals.is_empty() && failing.intervals.is_empty());
}
