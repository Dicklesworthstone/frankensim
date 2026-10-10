//! G1/G3: actual-binary adaptive h history, accepted heat and bounded refusal.
//! Synthetic charts exercise discrete behavior; no continuum error bound.
use super::*;

fn policy(tolerance_h: f64) -> J {
    J::parse(&format!(
        r#"{{"absolute_tolerance_k":1000000,"absolute_specific_enthalpy_tolerance_j_kg":{tolerance_h},"relative_tolerance":0,"minimum_trial_step_s":0.000001,"max_trials":4000}}"#
    )).unwrap()
}

fn request(mixed: bool, tolerance_h: f64) -> J {
    let mut input = J::parse(FIXTURE).unwrap();
    remove(&mut input, "radiation");
    let schedule = member(&mut input, "transient");
    put(schedule, "max_steps", number(4000.0));
    put(schedule, "adaptive", policy(tolerance_h));
    if mixed {
        let material = member(schedule, "enthalpy");
        remove(material, "initial_specific_enthalpy_j_kg");
        put(material, "initial_specific_enthalpies_j_kg",
            J::Array([500.0, 2000.0, 2100.0, 2200.0].into_iter().map(number).collect()));
    }
    input
}

fn history(result: &J) -> &[J] {
    result.get("transient").unwrap().get("history").unwrap().as_array().unwrap()
}

fn assert_accepted_accounting(result: &J) {
    let run = result.get("transient").unwrap();
    let adaptive = run.get("adaptive").unwrap();
    assert_eq!(adaptive.str_field("method"), Some("backward-euler-enthalpy-step-doubling"));
    let accepted = n(adaptive, "accepted_trials");
    near(n(run, "steps"), 2.0 * accepted, 0.0);
    near(n(adaptive, "trials") - n(adaptive, "rejected_trials"), accepted, 0.0);
    assert!(n(adaptive, "largest_accepted_error_ratio") <= 1.0);
    let rows = history(result);
    assert_eq!(rows.len(), n(run, "steps") as usize + 1);
    let (mut time, mut input, mut air, mut storage, mut accepted_work) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (index, row) in rows[1..].iter().enumerate() {
        let endpoint = n(row, "time_s");
        let dt = n(row, "dt_s");
        assert!(endpoint > time && dt > 0.0);
        near(dt, endpoint - time, 1e-14);
        assert!(!(time < 0.4 && endpoint > 0.4), "trial crossed a workload event");
        if index % 2 == 0 {
            assert!(matches!(row.get("estimated_local_error_ratio"), Some(J::Null)));
        } else {
            assert!(n(row, "estimated_local_error_ratio") <= 1.0);
        }
        input += dt * n(row, "source_w");
        air += dt * n(row, "air_heat_gain_w");
        storage += n(row, "stored_energy_change_j");
        accepted_work += n(row, "coupling_iterations");
        time = endpoint;
    }
    near(time, 0.8, 1e-14);
    near(input, 400.0, 1e-7);
    near(n(run, "input_energy_j"), input, 1e-7);
    near(n(run, "air_energy_gain_j"), air, 1e-7);
    near(n(run, "stored_energy_change_j"), storage, 1e-7);
    near(storage - input + air, n(run, "energy_residual_j"), 1e-7);
    assert!(n(run, "forward_solid_solves") > accepted_work,
        "coarse and rejected solid solves must remain charged");
}

#[test]
fn adaptive_plateau_matches_analytic_energy_and_carries_h_between_cycles() {
    let input = request(false, 0.01);
    let result = run(&input);
    assert_accepted_accounting(&result);
    let trajectory = result.get("transient").unwrap();
    near(n(trajectory.get("adaptive").unwrap(), "rejected_trials"), 0.0, 0.0);
    for (actual, expected) in values(&result, "solid_specific_enthalpies_j_kg")
        .iter().zip(enthalpies(0.8, None)) {
        near(*actual, expected, 2e-6);
    }
    for row in history(&result) {
        let expected = enthalpies(n(row, "time_s"), None);
        near(n(row, "minimum_specific_enthalpy_j_kg"),
            expected.into_iter().fold(f64::INFINITY, f64::min), 2e-6);
        near(n(row, "maximum_specific_enthalpy_j_kg"),
            expected.into_iter().fold(f64::NEG_INFINITY, f64::max), 2e-6);
    }

    let mut repeated = input.clone();
    put(member(&mut repeated, "transient"), "repeat",
        J::parse(r#"{"cycles":3,"max_total_steps":100}"#).unwrap());
    let repeated = run(&repeated);
    let cycles = repeated.get("repeated_cycles").unwrap();
    near(n(cycles, "cycles_completed"), 3.0, 0.0);
    near(n(cycles, "total_accepted_steps"), 24.0, 0.0);
    near(n(cycles, "input_energy_j"), 1200.0, 1e-7);
    near(n(cycles, "air_energy_gain_j"), 2.4 * heat(None).1, 1e-6);
    let mass = RHO / 24.0;
    let once = enthalpies(0.8, None);
    let actual = values(&repeated, "solid_specific_enthalpies_j_kg");
    for (&actual, once) in actual.iter().zip(once) {
        near(actual, 2000.0 + 3.0 * (once - 2000.0), 6e-6);
        assert!((1000.0..3000.0).contains(&actual));
    }
    near(n(cycles, "stored_energy_change_j"),
        actual.iter().map(|h| mass * (h - 2000.0)).sum(), 1e-6);
    // The last-cycle report remains local even though the cycle owner sums all work.
    assert_accepted_accounting(&repeated);
    assert!(n(cycles, "total_solid_solves")
        > n(repeated.get("transient").unwrap(), "total_solid_solves"));
}

#[test]
fn latent_h_discrepancy_refines_and_accepted_pair_replays_as_fixed_history() {
    // The sensible first node changes the heat flux into three latent nodes.
    // A huge temperature tolerance isolates the required full-field h criterion.
    let loose = run(&request(true, 1_000_000.0));
    let input = request(true, 0.002);
    let result = run(&input);
    assert_accepted_accounting(&result);
    let trajectory = result.get("transient").unwrap();
    assert!(n(trajectory.get("adaptive").unwrap(), "rejected_trials") > 0.0);
    assert!(n(trajectory, "steps") > n(loose.get("transient").unwrap(), "steps"));
    for temperature in &values(&result, "solid_temperatures_k")[1..] {
        near(*temperature, 350.0, 1e-8);
    }

    // Re-solve exactly the accepted half-step durations, with one fixed step
    // per interval. Any leaked coarse/rejected h changes this final field.
    let intervals = history(&result)[1..].iter().map(|row| {
        J::Object(vec![
            ("duration_s".into(), number(n(row, "dt_s"))),
            ("steps".into(), number(1.0)),
            ("power_scale".into(), number(if n(row, "interval") == 0.0 { 1.0 } else { 0.0 })),
        ])
    }).collect();
    let mut replay = input.clone();
    let schedule = member(&mut replay, "transient");
    remove(schedule, "adaptive");
    put(schedule, "max_step_s", number(1.0));
    put(schedule, "intervals", J::Array(intervals));
    let replay = run(&replay);
    for field in ["solid_specific_enthalpies_j_kg", "solid_temperatures_k", "solid_liquid_mass_fractions"] {
        let expected = values(&replay, field);
        for (actual, expected) in values(&result, field).iter().zip(expected) {
            near(*actual, expected, 2e-6);
        }
    }
    for field in ["input_energy_j", "stored_energy_change_j", "air_energy_gain_j"] {
        near(n(trajectory, field), n(replay.get("transient").unwrap(), field), 2e-6);
    }
}

#[test]
fn adaptive_enthalpy_refuses_exhausted_trials_resolution_and_endpoint_budgets() {
    let input = request(true, 0.002);
    let mut missing = input.clone();
    remove(member(member(&mut missing, "transient"), "adaptive"),
        "absolute_specific_enthalpy_tolerance_j_kg");
    assert!(refuses(&missing).contains("absolute_specific_enthalpy_tolerance_j_kg"));

    let mut trial_budget = input.clone();
    put(member(member(&mut trial_budget, "transient"), "adaptive"), "max_trials", number(1.0));
    assert!(refuses(&trial_budget).contains("trial budget"));

    let mut minimum = input.clone();
    put(member(member(&mut minimum, "transient"), "adaptive"),
        "minimum_trial_step_s", number(0.2));
    assert!(refuses(&minimum).contains("minimum"));

    let mut endpoints = input.clone();
    put(member(&mut endpoints, "transient"), "max_steps", number(8.0));
    assert!(refuses(&endpoints).contains("endpoint budget"));

    let mut repeated = input.clone();
    put(member(&mut repeated, "transient"), "repeat",
        J::parse(r#"{"cycles":3,"max_total_steps":24}"#).unwrap());
    assert!(refuses(&repeated).contains("endpoint budget"));

    let mut adjoint = input;
    put(member(&mut adjoint, "transient"), "adjoint",
        J::parse(r#"{"qoi":"final","max_checkpoint_bytes":1048576}"#).unwrap());
    assert!(refuses(&adjoint).contains("fixed timesteps"));
}
