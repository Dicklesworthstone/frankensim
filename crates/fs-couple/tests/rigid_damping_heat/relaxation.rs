//! Original Maxwell storage, finite heat and accepted history continuation.

#[path = "relaxation/temperature.rs"]
mod temperature;
use super::*;
use fs_couple::rigid_damping_heat::{
    driven::network::{DrivenNetworkDampingThermalConfig, advance_candidate as advance_driven},
    relaxation::{
        HeldInternalEffort, RelaxationDampingImplicitConfig, RelaxationDampingThermalConfig,
        RelaxationDampingThermalInput, RelaxationNetworkLimits, RelaxationNetworkMotionControls,
        SpatialRelaxationBranch, advance_candidate as advance_memory,
        advance_implicit_candidate as advance_memory_implicit,
    },
};
use fs_mbd::spatial_loads::{
    internal_impulse::viscous::network::driven::DrivenViscousNetworkLimits,
    internal_motion::viscous::network::driven::DrivenViscousNetworkMotionControls,
};

fn memory_controls() -> RelaxationDampingThermalConfig {
    let original = config();
    RelaxationDampingThermalConfig {
        mechanical: RelaxationNetworkMotionControls {
            source: original.mechanical.source,
            limits: RelaxationNetworkLimits {
                source: DrivenViscousNetworkLimits {
                    maximum_bodies: 2,
                    maximum_connectors: 1,
                    maximum_live_bytes: 128_000_000,
                    maximum_preparation_terms: 10_000_000_000,
                },
                maximum_branches: 4,
            },
        },
        thermal: original.thermal,
        maximum_heat_bindings: 2,
        energy_tolerance_j: original.energy_tolerance_j,
    }
}
fn memory_implicit_controls(scheme: ImplicitThermalScheme) -> RelaxationDampingImplicitConfig {
    let common = implicit_config(scheme);
    RelaxationDampingImplicitConfig {
        mechanical: memory_controls().mechanical,
        thermal: common.thermal,
        scheme,
        maximum_heat_bindings: 2,
        energy_tolerance_j: common.energy_tolerance_j,
    }
}
struct MemoryPair {
    bodies: Bodies,
    links: [ViscousNetworkLink; 1],
    branches: [Vec<SpatialRelaxationBranch>; 1],
    memory: [Vec<[f64; 6]>; 1],
    efforts: [HeldInternalEffort; 1],
}
impl MemoryPair {
    fn new(c: f64, root_k: f64, tau: f64, q: f64, u: f64, force: f64) -> Self {
        let mut bodies = Bodies::new([c, 0.0]);
        for (i, speed) in [-0.5 * u, 0.5 * u].into_iter().enumerate() {
            bodies.states[i] = RigidBodyState::new(
                bodies.states[i].pose(),
                Vec3::new(2.0 * speed, 0.0, 0.0),
                Vec3::ZERO,
            )
            .unwrap();
        }
        let mut rows = [[0.0; 6]; 6];
        rows[0][0] = root_k;
        Self {
            bodies,
            links: [axial([0, 1], c)],
            branches: [vec![
                SpatialRelaxationBranch::from_gram_rows(rows, tau).unwrap(),
            ]],
            memory: [vec![[q, 0.0, 0.0, 0.0, 0.0, 0.0]]],
            efforts: [HeldInternalEffort::MaterialFirst([
                force, 0.0, 0.0, 0.0, 0.0, 0.0,
            ])],
        }
    }
    fn input<'a>(
        &'a self,
        routes: &'a [BTreeMap<CellId, f64>],
    ) -> RelaxationDampingThermalInput<'a> {
        RelaxationDampingThermalInput {
            mechanical: RigidDampingThermalInput {
                states: &self.bodies.states[..2],
                properties: &self.bodies.properties[..2],
                links: &self.links,
                loads: HeldExternalLoads {
                    wrenches: &self.bodies.loads[..2],
                    gravity: self.bodies.gravity,
                },
                ..self.bodies.input(routes)
            },
            branches: &self.branches,
            memory: &self.memory,
            efforts: &self.efforts,
        }
    }
}

// Independent original two-body equations. Both masses are 2 kg, so the
// relative mobility is 1/kg and relative kinetic storage is u²/2.
// q is the original Maxwell energy amplitude: H_memory=q²/2.
fn reference(c: f64, root_k: f64, tau: f64, force: f64, u: f64, q: f64, t: f64) -> [f64; 4] {
    let f = |s: [f64; 4]| {
        [
            force - c * s[0] - root_k * s[1],
            root_k * s[0] - s[1] / tau,
            c * s[0] * s[0] + s[1] * s[1] / tau,
            force * s[0],
        ]
    };
    let add = |a: [f64; 4], b: [f64; 4], h: f64| core::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [u, q, 0.0, 0.0];
    let h = t / 16384.0;
    for _ in 0..16384 {
        let a = f(state);
        let b = f(add(state, a, h * 0.5));
        let c = f(add(state, b, h * 0.5));
        let d = f(add(state, c, h));
        for i in 0..4 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}
fn relative(states: &[RigidBodyState], bodies: &Bodies) -> f64 {
    states[1]
        .center_of_mass_velocity_world(bodies.properties[1])
        .unwrap()
        .x
        - states[0]
            .center_of_mass_velocity_world(bodies.properties[0])
            .unwrap()
            .x
}

#[test]
fn implicit_maxwell_retries_and_continuation_preserve_original_memory_and_work() {
    let fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, 0.7);
    let routes = [BTreeMap::from([(CellId(0), 0.7), (CellId(1), 0.2)])];
    let input = fixture.input(&routes);
    let calorics = laws();
    let initial = stiff_thermal(&calorics, 128);
    let saved_cells = initial.cells().clone();
    let saved_memory = fixture.memory.clone();
    let graph = stiff_graph(&calorics, &initial);
    let observer = Observe::default();
    let expected = reference(0.4, 1.3, 0.2, 0.7, -2.0, 0.6, 0.25);
    for scheme in [
        ImplicitThermalScheme::BackwardEuler,
        ImplicitThermalScheme::Sdirk2,
    ] {
        let mut controls = memory_implicit_controls(scheme);
        controls.thermal.maximum_newton_iterations = 2;
        controls.thermal.cell_tolerance_j = 1e-8;
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (whole, native) = advance_memory_implicit(
            input,
            &initial,
            &observer,
            0.25,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        assert!(native.rejected_intervals > 0);
        assert_eq!(
            whole.observer().hooks[8],
            native.energy.substeps + native.rejected_intervals
        );
        assert_eq!(whole.observer().intervals.len(), native.energy.substeps);
        assert_eq!(native.energy, whole.report().thermal);
        near(
            relative(whole.states(), &fixture.bodies),
            expected[0],
            8e-11,
        );
        near(whole.memory()[0][0][0], expected[1], 8e-11);
        near(whole.report().dissipated_energy_j, expected[2], 8e-11);
        near(whole.report().actuator_work_j, expected[3], 8e-11);
        near(
            whole.report().memory_energy_change_j,
            0.5 * (expected[1] * expected[1] - 0.6 * 0.6),
            8e-11,
        );
        near(whole.report().deposited_heat_j, 0.9 * expected[2], 8e-11);
        near(
            whole.report().unallocated_energy_j,
            0.1 * expected[2],
            8e-11,
        );
        near(whole.report().combined_energy_residual_j, 0.0, 1e-8);

        // Split at an actually accepted endpoint. Rejected larger intervals
        // cannot alter the memory that either continuation starts from.
        let split = whole.observer().intervals[0].2;
        let (first, first_native) = advance_memory_implicit(
            input,
            &initial,
            &observer,
            split,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        let continued = RelaxationDampingThermalInput {
            mechanical: RigidDampingThermalInput {
                time_s: first.report().time_s,
                states: first.states(),
                ..input.mechanical
            },
            memory: first.memory(),
            ..input
        };
        let (second, second_native) = advance_memory_implicit(
            continued,
            first.thermal(),
            first.observer(),
            0.25,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        assert_eq!(second.states(), whole.states());
        assert_eq!(second.memory(), whole.memory());
        assert_eq!(second.thermal().cells(), whole.thermal().cells());
        assert_eq!(
            first_native.energy.substeps + second_native.energy.substeps,
            native.energy.substeps
        );
        for (a, b, all) in [
            (
                first.report().dissipated_energy_j,
                second.report().dissipated_energy_j,
                whole.report().dissipated_energy_j,
            ),
            (
                first.report().actuator_work_j,
                second.report().actuator_work_j,
                whole.report().actuator_work_j,
            ),
            (
                first.report().memory_energy_change_j,
                second.report().memory_energy_change_j,
                whole.report().memory_energy_change_j,
            ),
        ] {
            near(a + b, all, 8e-11);
        }
        let failing = Observe {
            fail_at: Some(2),
            ..Default::default()
        };
        assert!(matches!(
            advance_memory_implicit(
                input,
                &initial,
                &failing,
                0.25,
                &graph,
                controls,
                &mut workspace,
                || false
            ),
            Err(RigidDampingThermalError::Observer(
                DepositionError::InvalidInput("synthetic observer refusal")
            ))
        ));
        assert!(failing.intervals.is_empty());
    }
    assert_eq!(fixture.memory, saved_memory);
    assert_eq!(initial.cells(), &saved_cells);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(observer, Observe::default());
}

#[test]
fn g1_preloaded_memory_accelerates_bodies_and_heats_original_finite_cells() {
    let fixture = MemoryPair::new(0.1, 2.0, 0.3, 2.0, 0.0, 0.0);
    let routes = [BTreeMap::from([(CellId(0), 0.3), (CellId(1), 0.4)])];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let target = 0.4;
    let expected = reference(0.1, 2.0, 0.3, 0.0, 0.0, 2.0, target);
    let out = advance_memory(
        fixture.input(&routes),
        &initial,
        &(),
        target,
        memory_controls(),
        || false,
    )
    .unwrap();
    near(relative(out.states(), &fixture.bodies), expected[0], 3e-11);
    near(out.memory()[0][0][0], expected[1], 3e-11);
    near(out.report().connector_heat_j[0], expected[2], 3e-11);
    near(out.report().dissipated_energy_j, expected[2], 3e-11);
    near(
        out.report().memory_energy_change_j,
        0.5 * (expected[1] * expected[1] - 4.0),
        3e-11,
    );
    near(out.report().mechanical_energy_change_j, -expected[2], 5e-11);
    near(out.report().deposited_heat_j, 0.7 * expected[2], 3e-11);
    near(out.report().unallocated_energy_j, 0.3 * expected[2], 3e-11);
    near(out.report().combined_energy_residual_j, 0.0, 5e-11);
    assert!(out.report().memory_energy_change_j < 0.0);
    assert!(relative(out.states(), &fixture.bodies).abs() > 0.1);
    assert_eq!(out.report().connector_actuator_work_j, [0.0]);
    assert_eq!(out.report().pair_relaxations, 0);
    for (id, fraction) in [(CellId(0), 0.3), (CellId(1), 0.4)] {
        let h = fraction * expected[2] / initial.cells()[&id].mass_kg;
        near(out.thermal().cells()[&id].specific_h_j_kg, h, 3e-11);
        near(
            out.thermal().temperature_k(id).unwrap(),
            calorics[id.0 as usize].temperature_k(h).unwrap(),
            3e-11,
        );
    }
    assert_eq!(fixture.memory[0][0][0], 2.0);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_driven_memory_continues_with_the_same_original_caloric_and_elastic_history() {
    let fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, 0.7);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let input = fixture.input(&routes);
    let expected = advance_memory(input, &initial, &(), 0.75, memory_controls(), || false).unwrap();
    let first = advance_memory(input, &initial, &(), 0.375, memory_controls(), || false).unwrap();
    let continued = RelaxationDampingThermalInput {
        mechanical: RigidDampingThermalInput {
            time_s: first.report().time_s,
            states: first.states(),
            ..input.mechanical
        },
        memory: first.memory(),
        ..input
    };
    let second = advance_memory(
        continued,
        first.thermal(),
        first.observer(),
        0.75,
        memory_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(second.states(), expected.states());
    assert_eq!(second.memory(), expected.memory());
    assert_eq!(second.thermal().cells(), expected.thermal().cells());
    near(
        first.report().dissipated_energy_j + second.report().dissipated_energy_j,
        expected.report().dissipated_energy_j,
        3e-11,
    );
    near(
        first.report().actuator_work_j + second.report().actuator_work_j,
        expected.report().actuator_work_j,
        3e-11,
    );
    near(
        first.report().memory_energy_change_j + second.report().memory_energy_change_j,
        expected.report().memory_energy_change_j,
        3e-11,
    );
    let original = reference(0.4, 1.3, 0.2, 0.7, -2.0, 0.6, 0.75);
    near(
        relative(second.states(), &fixture.bodies),
        original[0],
        3e-11,
    );
    near(second.memory()[0][0][0], original[1], 3e-11);
    near(expected.report().dissipated_energy_j, original[2], 3e-11);
    near(expected.report().actuator_work_j, original[3], 3e-11);
    near(expected.report().combined_energy_residual_j, 0.0, 5e-11);
}

#[test]
fn g3_empty_memory_matches_existing_driven_thermal_graph_without_elastic_storage() {
    let mut fixture = MemoryPair::new(0.7, 1.0, 0.2, 0.0, -2.0, 0.4);
    fixture.branches[0].clear();
    fixture.memory[0].clear();
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let controls = memory_controls();
    let input = fixture.input(&routes);
    let expected = advance_driven(
        input.mechanical,
        input.efforts,
        &initial,
        &(),
        0.5,
        DrivenNetworkDampingThermalConfig {
            mechanical: DrivenViscousNetworkMotionControls {
                source: controls.mechanical.source,
                limits: controls.mechanical.limits.source,
            },
            thermal: controls.thermal,
            maximum_heat_bindings: controls.maximum_heat_bindings,
            energy_tolerance_j: controls.energy_tolerance_j,
        },
        || false,
    )
    .unwrap();
    let out = advance_memory(input, &initial, &(), 0.5, controls, || false).unwrap();
    near(
        relative(out.states(), &fixture.bodies),
        relative(expected.states(), &fixture.bodies),
        3e-11,
    );
    near(
        out.report().dissipated_energy_j,
        expected.report().dissipated_energy_j,
        3e-11,
    );
    near(
        out.report().actuator_work_j,
        expected.report().actuator_work_j,
        3e-11,
    );
    near(
        out.thermal().temperature_k(CellId(0)).unwrap(),
        expected.thermal().temperature_k(CellId(0)).unwrap(),
        3e-11,
    );
    assert_eq!(out.report().memory_energy_change_j, 0.0);
    assert_eq!(expected.report().memory_energy_change_j, 0.0);
    assert!(out.memory()[0].is_empty());
}

#[test]
fn g4_memory_and_finite_heat_refuse_and_retry_as_one_candidate() {
    let fixture = MemoryPair::new(0.1, 2.0, 0.3, 2.0, 0.0, 0.2);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let observer = Observe::default();
    let input = fixture.input(&routes);
    let original_states = fixture.bodies.states;
    let original_memory = fixture.memory.clone();
    let mut polls = 0;
    let expected = advance_memory(input, &initial, &observer, 0.5, memory_controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        advance_memory(input, &initial, &observer, 0.5, memory_controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(RigidDampingThermalError::Cancelled)
    ));
    let failing = Observe {
        fail_at: Some(2),
        ..Default::default()
    };
    assert!(matches!(
        advance_memory(input, &initial, &failing, 0.5, memory_controls(), || false),
        Err(RigidDampingThermalError::Observer(_))
    ));
    let mut limited = memory_controls();
    limited.mechanical.limits.maximum_branches = 0;
    assert!(advance_memory(input, &initial, &observer, 0.5, limited, || false).is_err());
    let absent = [Vec::new()];
    assert!(
        advance_memory(
            RelaxationDampingThermalInput {
                memory: &absent,
                ..input
            },
            &initial,
            &observer,
            0.5,
            memory_controls(),
            || false,
        )
        .is_err()
    );
    let [mut first_law, second_law] = laws();
    first_law.maximum = 1e-8;
    let cold_laws = [first_law, second_law];
    let limited_thermal = thermal(&cold_laws, false, true);
    assert!(
        advance_memory(
            input,
            &limited_thermal,
            &observer,
            0.5,
            memory_controls(),
            || false,
        )
        .is_err()
    );
    assert_eq!(limited_thermal.time_s(), 0.0);
    let retry =
        advance_memory(input, &initial, &observer, 0.5, memory_controls(), || false).unwrap();
    assert_eq!(retry.states(), expected.states());
    assert_eq!(retry.memory(), expected.memory());
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(fixture.bodies.states, original_states);
    assert_eq!(fixture.memory, original_memory);
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.intervals.is_empty() && failing.intervals.is_empty());
}
