//! Actual finite-cell temperature drives original Kelvin and Maxwell kinetics.
use super::*;
use fs_couple::rigid_damping_heat::{
    driven::network::advance_temperature_candidate as advance_viscous_temperature,
    relaxation::temperature::{
        CalibratedRelaxationDamping, RelaxationTemperatureError,
        advance_candidate as advance_shifted_memory,
        advance_implicit_candidate as advance_shifted_memory_implicit,
    },
    temperature::shift::CalibratedKelvinDamping,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const GAS_CONSTANT: f64 = 8.31446261815324;

pub(crate) fn arrhenius(reference_k: f64, high_k: f64) -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 20_000.0 * GAS_CONSTANT,
        },
        reference_k,
        [280.0, high_k],
    )
    .unwrap()
}

#[test]
fn implicit_calibrated_maxwell_heat_updates_kinetics_and_preserves_continuation() {
    let fixture = MemoryPair::new(0.4, 1.3, 0.2, 2.0, -2.0, 0.7);
    let routes = [BTreeMap::from([(CellId(0), 0.7), (CellId(1), 0.2)])];
    let input = fixture.input(&routes);
    let calorics = laws();
    let initial = stiff_thermal(&calorics, 8);
    let saved_cells = initial.cells().clone();
    let saved_memory = fixture.memory.clone();
    let graph = stiff_graph(&calorics, &initial);
    let shifts = [arrhenius(300.0, 360.0)];
    let binding =
        CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &shifts).unwrap();
    let observer = Observe::default();
    for scheme in [
        ImplicitThermalScheme::BackwardEuler,
        ImplicitThermalScheme::Sdirk2,
    ] {
        let controls = memory_implicit_controls(scheme);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (whole, native) = advance_shifted_memory_implicit(
            input,
            &initial,
            &observer,
            &binding,
            &[CellId(0)],
            0.5,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        let (first, _) = advance_shifted_memory_implicit(
            input,
            &initial,
            &observer,
            &binding,
            &[CellId(0)],
            0.25,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        assert_eq!(native.energy.substeps, 2);
        assert_eq!(native.rejected_intervals, 0);
        assert_eq!(native.energy, whole.report().thermal);

        // Independent physical u/q/heat/work equations on the same accepted
        // held-temperature intervals. The first hot cell changes viscosity
        // AND the original Maxwell time, while its elastic form stays fixed.
        let first_reference = reference(0.4, 1.3, 0.2, 0.7, -2.0, 2.0, 0.25);
        let middle_t = first.thermal().temperature_k(CellId(0)).unwrap();
        let middle_factor = (20_000.0 * (1.0 / middle_t - 1.0 / 300.0)).exp();
        let second_reference = reference(
            0.4 * middle_factor,
            1.3,
            0.2 * middle_factor,
            0.7,
            first_reference[0],
            first_reference[1],
            0.25,
        );
        near(
            relative(whole.states(), &fixture.bodies),
            second_reference[0],
            1e-10,
        );
        near(whole.memory()[0][0][0], second_reference[1], 1e-10);
        near(
            whole.report().dissipated_energy_j,
            first_reference[2] + second_reference[2],
            1e-10,
        );
        near(
            whole.report().actuator_work_j,
            first_reference[3] + second_reference[3],
            1e-10,
        );
        near(
            whole.report().memory_energy_change_j,
            0.5 * (second_reference[1] * second_reference[1] - 4.0),
            1e-10,
        );
        assert!(middle_factor < 1.0);
        assert!(whole.law_evaluations() > 1);
        let temperature = whole.thermal().temperature_k(CellId(0)).unwrap();
        assert_eq!(whole.temperatures_k(), [temperature]);
        near(
            whole.shift_factors()[0],
            (20_000.0 * (1.0 / temperature - 1.0 / 300.0)).exp(),
            2e-14,
        );
        near(
            whole.report().unallocated_energy_j,
            0.1 * whole.report().dissipated_energy_j,
            1e-10,
        );
        near(whole.report().combined_energy_residual_j, 0.0, 1e-8);

        let continued = RelaxationDampingThermalInput {
            mechanical: RigidDampingThermalInput {
                time_s: first.report().time_s,
                states: first.states(),
                ..input.mechanical
            },
            memory: first.memory(),
            ..input
        };
        let (second, _) = advance_shifted_memory_implicit(
            continued,
            first.thermal(),
            first.observer(),
            &binding,
            &[CellId(0)],
            0.5,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        assert_eq!(second.states(), whole.states());
        assert_eq!(second.memory(), whole.memory());
        assert_eq!(second.thermal().cells(), whole.thermal().cells());
        assert_eq!(second.shift_factors(), whole.shift_factors());
        near(
            first.report().dissipated_energy_j + second.report().dissipated_energy_j,
            whole.report().dissipated_energy_j,
            1e-10,
        );

        let narrow_shifts = [arrhenius(300.0, 300.001)];
        let narrow =
            CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &narrow_shifts)
                .unwrap();
        assert!(matches!(
            advance_shifted_memory_implicit(
                input,
                &initial,
                &observer,
                &narrow,
                &[CellId(0)],
                0.25,
                &graph,
                controls,
                &mut workspace,
                || false
            ),
            Err(RigidDampingThermalError::Source {
                error: RelaxationTemperatureError::Shift(_),
                ..
            })
        ));
    }
    assert_eq!(fixture.memory, saved_memory);
    assert_eq!(fixture.branches[0][0].relaxation_time_s(), Some(0.2));
    assert_eq!(initial.cells(), &saved_cells);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(observer, Observe::default());
}

#[test]
fn g1_calibrated_isothermal_kinetics_match_independent_shifted_body_history_and_heat() {
    let fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, 0.7);
    // Heat leaves the mechanical source, so the explicit temperature cell stays
    // at 300 K. It remains a finite original caloric coordinate, not an inferred
    // interface average or an artificial fixed-temperature boundary.
    let routes = [BTreeMap::new()];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let wlf = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 0.8,
            c2_k: 60.0,
        },
        290.0,
        [280.0, 360.0],
    )
    .unwrap();
    let sources = [
        (
            arrhenius(290.0, 360.0),
            (20_000.0_f64 * (1.0 / 300.0 - 1.0 / 290.0)).exp(),
        ),
        (wlf, 10.0_f64.powf(-0.8 * 10.0 / 70.0)),
    ];
    for (shift, factor) in sources {
        let shifts = [shift];
        let binding =
            CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &shifts)
                .unwrap();
        let out = advance_shifted_memory(
            fixture.input(&routes),
            &initial,
            &(),
            &binding,
            &[CellId(0)],
            0.5,
            memory_controls(),
            || false,
        )
        .unwrap();
        let expected = reference(0.4 * factor, 1.3, 0.2 * factor, 0.7, -2.0, 0.6, 0.5);
        near(relative(out.states(), &fixture.bodies), expected[0], 5e-11);
        near(out.memory()[0][0][0], expected[1], 5e-11);
        near(out.report().dissipated_energy_j, expected[2], 5e-11);
        near(out.report().actuator_work_j, expected[3], 5e-11);
        near(out.report().unallocated_energy_j, expected[2], 5e-11);
        near(out.report().deposited_heat_j, 0.0, 0.0);
        near(out.shift_factors()[0], factor, 4e-15);
        assert_eq!(out.temperatures_k(), [300.0]);
        assert_eq!(out.law_evaluations(), 1);
        assert_eq!(fixture.branches[0][0].relaxation_time_s(), Some(0.2));
        assert_eq!(fixture.memory[0][0][0], 0.6);
        near(out.report().combined_energy_residual_j, 0.0, 1e-10);
    }
}

// Independent physical ODE: relative kinetic storage u²/2, original Maxwell
// storage q²/2, finite nonlinear cell E=m*(a*dT+b*dT²/2), and signed work.
// No production shift evaluation, dense source, loss observer or caloric inverse
// is used by this reference.
fn coupled_reference(t: f64) -> [f64; 4] {
    let rhs = |s: [f64; 4]| {
        let e = s[2] / 2.0;
        let temperature = 300.0 + 2.0 * e / (0.5 + (0.25 + 0.4 * e).sqrt());
        let factor = (20_000.0 * (1.0 / temperature - 1.0 / 300.0)).exp();
        let viscosity = 0.4 * factor;
        let tau = 0.2 * factor;
        [
            0.7 - viscosity * s[0] - 1.3 * s[1],
            1.3 * s[0] - s[1] / tau,
            viscosity * s[0] * s[0] + s[1] * s[1] / tau,
            0.7 * s[0],
        ]
    };
    let add = |a: [f64; 4], b: [f64; 4], h: f64| core::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [-2.0, 2.0, 0.0, 0.0];
    let h = t / 16384.0;
    for _ in 0..16384 {
        let a = rhs(state);
        let b = rhs(add(state, a, 0.5 * h));
        let c = rhs(add(state, b, 0.5 * h));
        let d = rhs(add(state, c, h));
        for i in 0..4 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn g1_self_heating_relaxation_refines_to_independent_coupled_physical_ode() {
    let fixture = MemoryPair::new(0.4, 1.3, 0.2, 2.0, -2.0, 0.7);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let shifts = [arrhenius(300.0, 360.0)];
    let binding =
        CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &shifts).unwrap();
    let expected = coupled_reference(0.5);
    let e = expected[2] / 2.0;
    let expected_t = 300.0 + 2.0 * e / (0.5 + (0.25 + 0.4 * e).sqrt());
    let mut errors = [0.0; 3];
    for (rung, steps) in [16, 32, 64].into_iter().enumerate() {
        let mut controls = memory_controls();
        controls.thermal.max_dt_s = 0.5 / steps as f64;
        let out = advance_shifted_memory(
            fixture.input(&routes),
            &initial,
            &(),
            &binding,
            &[CellId(0)],
            0.5,
            controls,
            || false,
        )
        .unwrap();
        let temperature = out.thermal().temperature_k(CellId(0)).unwrap();
        errors[rung] = (relative(out.states(), &fixture.bodies) - expected[0]).abs()
            + (out.memory()[0][0][0] - expected[1]).abs()
            + (out.report().dissipated_energy_j - expected[2]).abs()
            + (out.report().actuator_work_j - expected[3]).abs()
            + (temperature - expected_t).abs();
        assert!(temperature > 300.5);
        assert!(out.shift_factors()[0] < 0.9);
        assert!(out.law_evaluations() > 1);
        near(out.temperatures_k()[0], temperature, 0.0);
        near(
            out.shift_factors()[0],
            (20_000.0 * (1.0 / temperature - 1.0 / 300.0)).exp(),
            2e-14,
        );
        near(out.report().combined_energy_residual_j, 0.0, 2e-10);
        near(
            out.report().mechanical_energy_change_j + out.report().dissipated_energy_j
                - out.report().actuator_work_j,
            0.0,
            2e-10,
        );
        near(
            out.thermal().cells()[&CellId(0)].specific_h_j_kg * 2.0,
            out.report().dissipated_energy_j,
            2e-10,
        );
    }
    for pair in errors.windows(2) {
        assert!((1.8..2.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn g1_original_calibration_and_memory_continue_across_native_thermal_events() {
    let fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, 0.7);
    let routes = [BTreeMap::from([(CellId(0), 0.7), (CellId(1), 0.2)])];
    let calorics = laws();
    // Births are processed on the original timeline; an external heat change
    // at .35 must be respected by both the one-call and continuation paths.
    let initial = thermal(&calorics, true, false);
    let shifts = [arrhenius(300.0, 360.0)];
    let binding =
        CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &shifts).unwrap();
    let input = fixture.input(&routes);
    let mut controls = memory_controls();
    controls.maximum_heat_bindings = 2;
    let one = advance_shifted_memory(
        input,
        &initial,
        &(),
        &binding,
        &[CellId(0)],
        0.75,
        controls,
        || false,
    )
    .unwrap();
    // The existing event splits [.25,.375] at .35. Its following full interval
    // ends at .475, so continuing there preserves the exact same native grid.
    let first = advance_shifted_memory(
        input,
        &initial,
        &(),
        &binding,
        &[CellId(0)],
        0.475,
        controls,
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
    let second = advance_shifted_memory(
        continued,
        first.thermal(),
        first.observer(),
        &binding,
        &[CellId(0)],
        0.75,
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(second.states(), one.states());
    assert_eq!(second.memory(), one.memory());
    assert_eq!(second.thermal().cells(), one.thermal().cells());
    assert_eq!(second.temperatures_k(), one.temperatures_k());
    assert_eq!(second.shift_factors(), one.shift_factors());
    for (a, b, whole) in [
        (
            first.report().dissipated_energy_j,
            second.report().dissipated_energy_j,
            one.report().dissipated_energy_j,
        ),
        (
            first.report().actuator_work_j,
            second.report().actuator_work_j,
            one.report().actuator_work_j,
        ),
        (
            first.report().memory_energy_change_j,
            second.report().memory_energy_change_j,
            one.report().memory_energy_change_j,
        ),
        (
            first.report().unallocated_energy_j,
            second.report().unallocated_energy_j,
            one.report().unallocated_energy_j,
        ),
    ] {
        near(a + b, whole, 2e-10);
    }
    near(
        one.report().unallocated_energy_j,
        0.1 * one.report().dissipated_energy_j,
        2e-10,
    );
    assert!(initial.cells().is_empty());
}

#[test]
fn g3_empty_histories_match_the_existing_calibrated_viscous_feedback() {
    let mut fixture = MemoryPair::new(0.7, 1.0, 0.2, 0.0, -2.0, 0.4);
    fixture.branches[0].clear();
    fixture.memory[0].clear();
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let shifts = [arrhenius(300.0, 360.0)];
    let binding =
        CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &shifts).unwrap();
    let viscous = CalibratedKelvinDamping::try_new(&fixture.links, &shifts).unwrap();
    let input = fixture.input(&routes);
    let control = memory_controls();
    let one = advance_shifted_memory(
        input,
        &initial,
        &(),
        &binding,
        &[CellId(0)],
        0.5,
        control,
        || false,
    )
    .unwrap();
    let expected = advance_viscous_temperature(
        input.mechanical,
        input.efforts,
        &initial,
        &(),
        &viscous,
        &[CellId(0)],
        0.5,
        DrivenNetworkDampingThermalConfig {
            mechanical: DrivenViscousNetworkMotionControls {
                source: control.mechanical.source,
                limits: control.mechanical.limits.source,
            },
            thermal: control.thermal,
            maximum_heat_bindings: control.maximum_heat_bindings,
            energy_tolerance_j: control.energy_tolerance_j,
        },
        || false,
    )
    .unwrap();
    near(
        relative(one.states(), &fixture.bodies),
        relative(expected.states(), &fixture.bodies),
        4e-11,
    );
    near(
        one.report().dissipated_energy_j,
        expected.report().dissipated_energy_j,
        4e-11,
    );
    near(
        one.report().actuator_work_j,
        expected.report().actuator_work_j,
        4e-11,
    );
    near(one.temperatures_k()[0], expected.temperatures_k()[0], 4e-11);
    assert_eq!(one.report().memory_energy_change_j, 0.0);
    assert!(one.memory()[0].is_empty());
}

#[test]
fn g4_shift_domain_reference_and_late_refusals_preserve_all_original_histories() {
    let fixture = MemoryPair::new(0.1, 2.0, 0.3, 2.0, 0.0, 0.2);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let shifts = [arrhenius(300.0, 360.0)];
    let binding =
        CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &shifts).unwrap();
    let observer = Observe::default();
    let input = fixture.input(&routes);
    let original_memory = fixture.memory.clone();
    let mut polls = 0;
    let expected = advance_shifted_memory(
        input,
        &initial,
        &observer,
        &binding,
        &[CellId(0)],
        0.5,
        memory_controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        advance_shifted_memory(
            input,
            &initial,
            &observer,
            &binding,
            &[CellId(0)],
            0.5,
            memory_controls(),
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(RigidDampingThermalError::Cancelled)
    ));
    let limited = [arrhenius(300.0, 300.001)];
    let limited_binding =
        CalibratedRelaxationDamping::try_new(&fixture.links, &fixture.branches, &limited).unwrap();
    // Exactly one native interval: the final endpoint refresh is the first
    // place this actual caloric temperature leaves the calibrated support.
    assert!(matches!(
        advance_shifted_memory(
            input,
            &initial,
            &observer,
            &limited_binding,
            &[CellId(0)],
            0.125,
            memory_controls(),
            || false,
        ),
        Err(RigidDampingThermalError::Source {
            error: RelaxationTemperatureError::Shift(_),
            ..
        })
    ));
    let changed = [vec![
        SpatialRelaxationBranch::from_gram_rows(*fixture.branches[0][0].gram_rows(), 0.6).unwrap(),
    ]];
    assert!(matches!(
        advance_shifted_memory(
            RelaxationDampingThermalInput {
                branches: &changed,
                ..input
            },
            &initial,
            &observer,
            &binding,
            &[CellId(0)],
            0.5,
            memory_controls(),
            || false,
        ),
        Err(RigidDampingThermalError::Source {
            error: RelaxationTemperatureError::Reference { connector: 0 },
            ..
        })
    ));
    assert!(
        advance_shifted_memory(
            input,
            &initial,
            &observer,
            &binding,
            &[CellId(99)],
            0.5,
            memory_controls(),
            || false,
        )
        .is_err()
    );
    let failing = Observe {
        fail_at: Some(2),
        ..Default::default()
    };
    assert!(matches!(
        advance_shifted_memory(
            input,
            &initial,
            &failing,
            &binding,
            &[CellId(0)],
            0.5,
            memory_controls(),
            || false,
        ),
        Err(RigidDampingThermalError::Observer(_))
    ));
    let retry = advance_shifted_memory(
        input,
        &initial,
        &observer,
        &binding,
        &[CellId(0)],
        0.5,
        memory_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(retry.states(), expected.states());
    assert_eq!(retry.memory(), expected.memory());
    assert_eq!(retry.thermal().cells(), expected.thermal().cells());
    assert_eq!(retry.shift_factors(), expected.shift_factors());
    assert_eq!(fixture.memory, original_memory);
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.intervals.is_empty() && failing.intervals.is_empty());
}
