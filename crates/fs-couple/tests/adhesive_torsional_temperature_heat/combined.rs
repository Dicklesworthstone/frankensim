//! One thermally relaxed disk carries normal, vector shear and twist together.
//! The parent supplies original material cards and the independently soluble
//! two-cell thermal system; no separate mode radius or fitted heat law is used.

use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::combined::motion::JkrCombinedMotionBudget;
use fs_couple::friction_heat::adhesive::combined_temperature::*;
use fs_mbd::tangential_friction::adhesive::jkr_combined::*;

const PRELOAD: [f64; 3] = [3e-5, -2e-5, 0.004];

fn combined_controls() -> RigidJkrCombinedControls {
    let old = native_controls();
    RigidJkrCombinedControls {
        motion: JkrCombinedMotionBudget {
            maximum_step_s: old.motion.maximum_step_s,
            maximum_coordinate_increment: [0.001, 0.001, 0.001, 0.1],
            maximum_rate: [100.0; 4],
            maximum_effort: [1e4; 4],
            maximum_iterations: 64,
            maximum_evaluations: 8192,
            maximum_backtracks: 24,
            coordinate_absolute_tolerance: [1e-16; 4],
            coordinate_relative_tolerance: 1e-12,
            rate_absolute_tolerance: [1e-13; 4],
            rate_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-13,
            energy_relative_tolerance: 1e-10,
        },
        maximum_bodies: old.maximum_bodies,
        geometry_tolerance_m: old.geometry_tolerance_m,
        displacement_tolerance_m: old.displacement_tolerance_m,
        twist_tolerance_rad: old.twist_tolerance_rad,
        velocity_tolerance_m_s: old.velocity_tolerance_m_s,
        angular_velocity_tolerance_rad_s: old.angular_velocity_tolerance_rad_s,
        linear_momentum_tolerance_ns: old.linear_momentum_tolerance_ns,
        angular_momentum_tolerance_nms: old.angular_momentum_tolerance_nms,
        energy_tolerance_j: old.energy_tolerance_j,
    }
}

impl Pair {
    fn combined<'a>(&self, source: &'a SourceJkrMixedMode) -> PreparedRigidJkrCombinedPair<'a> {
        PreparedRigidJkrCombinedPair::new(
            &self.bodies,
            &self.properties,
            RigidJkrCombinedContact {
                binding: RigidContactBinding {
                    body_a: 0,
                    body_b: 1,
                    arm_a_body: Vec3::new(0.0, 0.0, 0.5),
                    arm_b_body: Vec3::new(0.0, 0.0, -0.5),
                },
                frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
                reference_approach_m: 0.0,
                source,
            },
            combined_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    }

    fn hold_combined(
        &mut self,
        source: &SourceJkrMixedMode,
        previous: Option<RigidJkrCombinedState>,
        preload: [f64; 3],
    ) -> RigidJkrCombinedState {
        let native = self.combined(source);
        let gate = CancelGate::new_clock_free();
        let state = previous.unwrap_or_else(|| native.initialize(preload, &gate).unwrap());
        let effort = state.motion_state().contact().force;
        let force = Vec3::new(effort[1], effort[2], effort[0]);
        let couple = Vec3::new(0.0, 0.0, effort[3]);
        let arms = native.reaction_arms_body();
        // Holding the full native pair requires the original point-force arm
        // torque as well as the free torsional couple, in each body's frame.
        let loads = core::array::from_fn(|i| {
            let sign = if i == 0 { 1.0 } else { -1.0 };
            let force_world = force.scale(sign);
            let orientation = self.bodies[i].pose().orientation();
            Wrench {
                force_world,
                torque_body: arms[i]
                    .cross(orientation.rotate_world_to_body(force_world))
                    .add(orientation.rotate_world_to_body(couple.scale(sign))),
            }
        });
        let step = native
            .advance_with_wrenches(state, loads, DT, &gate)
            .unwrap();
        self.bodies = step.bodies_after;
        step.state_after
    }
}

// Independently differentiated full potential at R=1 and delta=0:
// E*a^3 + 8*Gbar*a*phi^2 + 4*Gstar*|u|^2/a = 2*pi*W.
// The larger root is the stable opening branch. These tests stay well away
// from its fold; the sign and curvature checks make that bracket explicit.
fn independent_combined(t: AdhesiveTemperatures, preload: [f64; 3]) -> (f64, [f64; 4]) {
    let ea = 1e6 * (1.0 - 0.01 * (t.surface_a_k - 300.0));
    let eb = 2e6 * (1.0 - 0.01 * (t.surface_b_k - 300.0));
    let e = ((1.0 - 0.25_f64.powi(2)) * (1.0 / ea + 1.0 / eb)).recip();
    let gbar = (2.5 / ea + 2.5 / eb).recip();
    let gstar = (1.75 * (2.5 / ea + 2.5 / eb)).recip();
    let w = 1.0 - 0.03 * (t.interface_k - 300.0);
    let shear2 = preload[0] * preload[0] + preload[1] * preload[1];
    let phi2 = preload[2] * preload[2];
    let unsheared = independent_source(t, preload[2]).0;
    let residual =
        |a: f64| e * a.powi(3) + 8.0 * gbar * a * phi2 + 4.0 * gstar * shear2 / a - 2.0 * PI * w;
    let a = if shear2 == 0.0 {
        unsheared
    } else {
        let mut low = 0.5 * unsheared;
        let mut high = unsheared;
        assert!(residual(low) < 0.0 && residual(high) > 0.0);
        for _ in 0..100 {
            let mid = 0.5 * (low + high);
            if residual(mid) > 0.0 {
                high = mid;
            } else {
                low = mid;
            }
        }
        0.5 * (low + high)
    };
    assert!(4.0 * e * a.powi(3) + 16.0 * gbar * a * phi2 - 2.0 * PI * w > 0.0);
    (
        a,
        [
            e * a.powi(5) / 5.0,
            4.0 * gstar * a * shear2,
            (8.0 / 3.0) * gbar * a.powi(3) * phi2,
            -PI * w * a * a,
        ],
    )
}

fn independent_combined_feedback(q: [f64; 2]) -> (AdhesiveTemperatures, f64, [f64; 4]) {
    let room = AdhesiveTemperatures {
        surface_a_k: 300.0,
        surface_b_k: 300.0,
        interface_k: 300.0,
    };
    let original = independent_combined(room, PRELOAD).0;
    let residual =
        |a: f64| a - independent_combined(temperatures_at_area(PI * a * a, q), PRELOAD).0;
    let mut low = 0.2 * original;
    let mut high = 2.0 * original;
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if residual(mid) > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    let a = 0.5 * (low + high);
    let temperatures = temperatures_at_area(PI * a * a, q);
    (
        temperatures,
        a,
        independent_combined(temperatures, PRELOAD).1,
    )
}

#[test]
fn one_shared_thermal_radius_updates_all_four_stores_and_retains_native_history() {
    let data = common::Fixture::new();
    let source = source(&data);
    let caloric = curves(&source);
    let mut pair = Pair::new();
    let state = pair.hold_combined(&source, None, PRELOAD);
    let native = pair.combined(&source);
    let before = state.motion_state().contact();
    assert!(state.coordinates()[0].abs() < 1e-13);
    assert!(before.tangential_elastic_energy_j > 0.0 && before.torsional_elastic_energy_j > 0.0);
    let thermal = finite_thermal(&caloric, before.contact_area_m2);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let cells = thermal.cells().clone();
    let q = [0.4, -0.08];
    let candidate = propose_jkr_combined_temperature_heat(
        &native,
        state,
        &thermal,
        binding(&data, &caloric, &route),
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| resolve(q, t),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let report = candidate.report();
    let (expected_t, expected_a, storage) = independent_combined_feedback(q);
    close(
        report.temperatures.surface_a_k,
        expected_t.surface_a_k,
        2e-7,
    );
    close(
        report.temperatures.surface_b_k,
        expected_t.surface_b_k,
        2e-7,
    );
    close(
        report.temperatures.interface_k,
        expected_t.interface_k,
        2e-7,
    );
    close(
        candidate.route().contact_area_m2(),
        PI * expected_a * expected_a,
        2e-12,
    );
    let after = report.temperature.after;
    close(after.contact_radius_m, expected_a, 2e-11);
    for (observed, expected) in [
        after.normal_elastic_energy_j,
        after.tangential_elastic_energy_j,
        after.torsional_elastic_energy_j,
        after.surface_energy_j,
    ]
    .into_iter()
    .zip(storage)
    {
        close(observed, expected, 1e-11);
    }
    close(
        report.temperature.parameter_work_j,
        storage.iter().sum::<f64>() - before.reversible_energy_j,
        2e-11,
    );
    for ((change, old), expected) in [
        report.temperature.normal_elastic_change_j,
        report.temperature.tangential_elastic_change_j,
        report.temperature.torsional_elastic_change_j,
        report.temperature.surface_energy_change_j,
    ]
    .into_iter()
    .zip([
        before.normal_elastic_energy_j,
        before.tangential_elastic_energy_j,
        before.torsional_elastic_energy_j,
        before.surface_energy_j,
    ])
    .zip(storage)
    {
        close(change, expected - old, 2e-11);
    }
    assert!(report.temperature.tangential_elastic_change_j.abs() > 1e-9);
    assert!(report.temperature.torsional_elastic_change_j.abs() > 1e-9);
    assert!(after.contact_area_m2 < before.contact_area_m2);
    // Omitting the elastic shear would relax a different disk, even at the
    // same three temperatures. Its area cannot price this loaded contact.
    let shear_free_radius = independent_source(report.temperatures, PRELOAD[2]).0;
    assert!((shear_free_radius - after.contact_radius_m).abs() > 1e-6);
    assert_eq!(
        report.parameter_work_source,
        "prescribed source-coefficient controller"
    );
    close(
        report.heat.thermal.stored_energy_change_j,
        q.iter().sum(),
        2e-9,
    );
    assert_eq!(report.heat.dissipated_j, 0.0);
    assert_eq!(report.heat.entropy.unwrap().friction_generation_j_k, 0.0);
    close(report.combined_energy_residual_j, 0.0, 4e-9);
    assert_eq!(
        candidate.state().coordinates().map(f64::to_bits),
        state.coordinates().map(f64::to_bits)
    );
    assert_eq!(
        candidate.state().rates().map(f64::to_bits),
        state.rates().map(f64::to_bits)
    );
    assert_eq!(
        candidate.state().pending_com_travel_m(),
        state.pending_com_travel_m()
    );
    assert_eq!(
        candidate.state().pending_angular_displacement_body_rad(),
        state.pending_angular_displacement_body_rad()
    );
    assert_eq!(candidate.state().steps(), state.steps());
    assert_eq!(candidate.state().time_s(), state.time_s());
    assert_eq!(candidate.thermal().time_s(), state.time_s());
    assert_eq!(
        report.temperatures.surface_a_k,
        candidate.thermal().temperature_k(CellId(0)).unwrap()
    );
    assert_eq!(
        report.temperatures.surface_b_k,
        candidate.thermal().temperature_k(CellId(1)).unwrap()
    );
    assert_eq!(candidate.bodies(), &pair.bodies);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(report.source, report.temperatures);
    let patch = candidate.thermal().contact_patch(LINK).unwrap();
    close(patch.link.area_m2, candidate.route().contact_area_m2(), 0.0);
    close(
        patch.link.conductance_w_k,
        series_conductance_w_k(patch.link.area_m2, [HALF_LENGTH; 2], K).unwrap(),
        1e-12,
    );
    let (successor, bodies, next_state, _, _, _) = candidate.into_parts();
    pair.bodies = bodies;
    let later = pair.hold_combined(&successor, Some(next_state), PRELOAD);
    assert_eq!(later.steps(), state.steps() + 1);
    assert_eq!(later.time_s(), 2.0 * DT);
    for (observed, expected) in later.coordinates()[1..].iter().zip(PRELOAD) {
        close(*observed, expected, 1e-12);
    }
}

#[test]
fn zero_shear_recovers_the_original_native_torsional_thermal_owner() {
    let data = common::Fixture::new();
    let source = source(&data);
    let caloric = curves(&source);
    let mut combined = Pair::new();
    let state = combined.hold_combined(&source, None, [0.0, 0.0, PRELOAD[2]]);
    let native = combined.combined(&source);
    let thermal = finite_thermal(&caloric, state.motion_state().contact().contact_area_m2);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let gate = CancelGate::new_clock_free();
    let actual = propose_jkr_combined_temperature_heat(
        &native,
        state,
        &thermal,
        binding(&data, &caloric, &route),
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| resolve([0.4, -0.08], t),
        &gate,
    )
    .unwrap();
    let mut original = Pair::new();
    let initial = original
        .prepared(&source)
        .initialize(PRELOAD[2], &gate)
        .unwrap();
    let old_state = original.hold_phase(&source, Some(initial));
    let old_native = original.prepared(&source);
    let old_thermal = finite_thermal(&caloric, old_state.motion_state().contact().contact_area_m2);
    let old_route = FrictionHeatRoute::from_link(&old_thermal, LINK, K).unwrap();
    let expected = propose_jkr_torsional_temperature_heat(
        &old_native,
        old_state,
        &old_thermal,
        binding(&data, &caloric, &old_route),
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| resolve([0.4, -0.08], t),
        &gate,
    )
    .unwrap();
    let a = actual.report();
    let b = expected.report();
    close(a.temperatures.surface_a_k, b.temperatures.surface_a_k, 2e-7);
    close(a.temperatures.surface_b_k, b.temperatures.surface_b_k, 2e-7);
    close(a.temperatures.interface_k, b.temperatures.interface_k, 2e-7);
    close(
        a.temperature.after.contact_radius_m,
        b.temperature.after.contact_radius_m,
        2e-11,
    );
    close(
        a.temperature.parameter_work_j,
        b.temperature.parameter_work_j,
        2e-11,
    );
    assert_eq!(a.temperature.after.tangential_elastic_energy_j, 0.0);
    assert_eq!(a.temperature.tangential_elastic_change_j, 0.0);
    close(
        a.heat.thermal.stored_energy_change_j,
        b.heat.thermal.stored_energy_change_j,
        2e-9,
    );
    for id in [CellId(0), CellId(1)] {
        close(
            actual.thermal().temperature_k(id).unwrap(),
            expected.thermal().temperature_k(id).unwrap(),
            2e-7,
        );
    }
}

#[test]
fn combined_receiver_domain_work_budget_and_final_cancellation_preserve_retry() {
    let data = common::Fixture::new();
    let source = source(&data);
    let caloric = curves(&source);
    let mut pair = Pair::new();
    let initial = pair
        .combined(&source)
        .initialize(PRELOAD, &CancelGate::new_clock_free())
        .unwrap();
    let state = pair.hold_combined(&source, Some(initial), PRELOAD);
    let native = pair.combined(&source);
    let thermal = finite_thermal(&caloric, state.motion_state().contact().contact_area_m2);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let cells = thermal.cells().clone();
    let bodies = pair.bodies;
    let original = binding(&data, &caloric, &route);
    let invoke = |b, s, limits| {
        propose_jkr_combined_temperature_heat(
            &native,
            s,
            &thermal,
            b,
            limits,
            &[0.0; 2],
            300.0,
            |_, t, _| resolve([0.4, -0.08], t),
            &CancelGate::new_clock_free(),
        )
    };
    let mut wrong = original;
    wrong.caloric.swap(0, 1);
    assert!(invoke(wrong, state, controls()).is_err());
    wrong = original;
    wrong.material_cards.swap(0, 1);
    assert!(invoke(wrong, state, controls()).is_err());
    wrong = original;
    wrong.body_ids.swap(0, 1);
    assert!(invoke(wrong, state, controls()).is_err());
    wrong = original;
    wrong.parameter_work_source = " ";
    assert!(invoke(wrong, state, controls()).is_err());
    assert!(invoke(original, initial, controls()).is_err());
    let mut limited = controls();
    limited.maximum_absolute_parameter_work_j = 1e-15;
    assert!(invoke(original, state, limited).is_err());
    limited = controls();
    limited.maximum_fractional_area_change = 0.0;
    assert!(invoke(original, state, limited).is_err());
    limited = controls();
    limited.thermal.maximum_trials = 1;
    assert!(invoke(original, state, limited).is_err());
    let stale = thermal
        .apply_heat_j(
            &BTreeMap::from([(CellId(0), CAPACITY[0])]),
            1e-12,
            1e-12,
            || Ok(()),
        )
        .unwrap()
        .0;
    let attempted = Cell::new(false);
    let refusal = propose_jkr_combined_temperature_heat(
        &native,
        state,
        &stale,
        original,
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| {
            attempted.set(true);
            resolve([0.4, -0.08], t)
        },
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(
        refusal,
        Err(JkrCombinedTemperatureHeatError::Invalid(
            "incoming bulk source temperatures differ from their accepted caloric cells"
        ))
    ));
    assert!(!attempted.get());
    let calls = Cell::new(0usize);
    let pilot = propose_jkr_combined_temperature_heat(
        &native,
        state,
        &thermal,
        original,
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| {
            calls.set(calls.get() + 1);
            resolve([0.4, -0.08], t)
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let final_call = calls.get();
    assert_eq!(final_call, pilot.report().thermal_trials);
    assert!(final_call > 2);
    calls.set(0);
    let gate = CancelGate::new_clock_free();
    let cancelled = propose_jkr_combined_temperature_heat(
        &native,
        state,
        &thermal,
        original,
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| {
            calls.set(calls.get() + 1);
            if calls.get() == final_call {
                gate.request();
            }
            resolve([0.4, -0.08], t)
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls.get(), final_call);
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(pair.bodies, bodies);
    assert_eq!(native.contact().source, &source);
    let retry = invoke(original, state, controls()).unwrap();
    assert_eq!(retry.state(), pilot.state());
    assert_eq!(retry.source(), pilot.source());
    assert_eq!(retry.report().temperature, pilot.report().temperature);
}
