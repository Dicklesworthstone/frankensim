use super::*;
use crate::brass_loop::nonlinear::{lip_face, tests::specimen};
use fs_duct::nonlinear::bore::{BoreSpec, NonlinearBore};
use fs_duct::{Duct, Segment};

fn bore(gas: &FrozenGas, dt: f64) -> NonlinearBore {
    NonlinearBore::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.006, length: 0.20 }] },
        gas.clone(),
        293.15,
        101_325.0,
        None,
        BoreSpec {
            maximum_cell_m: 0.01,
            inlet_volume_m3: 1.5e-6,
            step_s: dt,
            cfl: 0.45,
            wall: None,
        },
    ).unwrap()
}

fn stationary(gas: &FrozenGas, pressure: f64, temperature: f64) -> Primitive {
    gas.primitive(
        gas.conserved(pressure / (gas.gas_constant() * temperature), 0.0, pressure).unwrap(),
        || false,
    ).unwrap()
}

fn spec(volume: f64, feed: f64, pressure: f64, temperature: f64) -> FiniteMouthSpec {
    FiniteMouthSpec {
        volume_m3: volume,
        feed_area_m2: feed,
        initial_pressure_pa: pressure,
        initial_temperature_k: temperature,
    }
}

// This advances the actual nonlinear bore and its actual lip-face law. The
// compact volume sees the same native PortStage values as the production lip.
fn advance(
    gas: &FrozenGas,
    field: &mut NonlinearBore,
    supply: &mut FiniteMouth,
    mouth: &mut Primitive,
    opening: f64,
    dt: f64,
) -> FiniteMouthReport {
    advance_with(gas, field, supply, mouth, opening, dt, 0.0).0
}

// The same step with a series offset `series_pa` handed to the trial; also returns the mouth state each stage's lip face was offered.
fn advance_with(
    gas: &FrozenGas,
    field: &mut NonlinearBore,
    supply: &mut FiniteMouth,
    mouth: &mut Primitive,
    opening: f64,
    dt: f64,
    series_pa: f64,
) -> (FiniteMouthReport, Vec<Primitive>) {
    let before_field = field.solver().totals();
    let before_mouth = supply.state(*mouth);
    let inlet = field.inlet_area_m2();
    let mut trial = supply.trial(*mouth, dt, inlet, series_pa);
    let mut offered = Vec::new();
    let step = field.step(|stage| {
        let flux = trial.flux(gas, stage, |p| {
            offered.push(p);
            Ok(lip_face(gas, stage.left.1, p, opening, inlet)?.flux)
        })?;
        Ok(Some(flux))
    }).unwrap();
    *mouth = supply.commit(trial, &step.report).unwrap();
    let report = *supply.report().unwrap();
    let after_field = field.solver().totals();
    // The only external exchanges are the held feed and the radiating end.
    // There are no wall sources in this fixture.
    let mass_error = after_field.mass + report.after.mass_kg
        - before_field.mass - before_mouth.mass_kg
        - report.feed_mass_kg - step.report.right_boundary_change.mass;
    let energy_error = after_field.energy + report.after.energy_j
        - before_field.energy - before_mouth.energy_j
        - report.feed_energy_j - step.report.right_boundary_change.energy;
    assert!(mass_error.abs() <= 2e-12 * (before_field.mass + before_mouth.mass_kg));
    assert!(energy_error.abs() <= 2e-12 * (before_field.energy + before_mouth.energy_j).abs());
    assert!(report.mass_residual_kg.abs() <= 2e-12 * before_mouth.mass_kg);
    assert!(report.energy_residual_j.abs() <= 2e-12 * before_mouth.energy_j.abs());
    assert_eq!(report.bore_mass_kg, step.report.left_boundary_change.mass);
    assert_eq!(report.bore_energy_j, step.report.left_boundary_change.energy);
    assert!(report.maximum_cfl <= 0.45);
    (report, offered)
}

fn sealed_discharge(dt: f64) -> f64 {
    let gas = specimen();
    let mut field = bore(&gas, dt);
    let (mut supply, mut mouth) =
        FiniteMouth::new(&gas, spec(2e-6, 0.0, 150_000.0, 293.15), field.room()).unwrap();
    let initial = supply.state(mouth);
    for _ in 0..(160e-6 / dt).round() as usize {
        let report = advance(&gas, &mut field, &mut supply, &mut mouth, 1e-5, dt);
        assert_eq!(report.feed_mass_kg, 0.0);
        assert_eq!(report.feed_energy_j, 0.0);
        assert!(report.bore_mass_kg > 0.0 && report.bore_energy_j > 0.0);
    }
    let final_state = supply.state(mouth);
    assert!(final_state.pressure_pa < initial.pressure_pa - 100.0);
    assert!(final_state.temperature_k < initial.temperature_k);
    assert!(final_state.mass_kg < initial.mass_kg);
    // Independent adiabatic tank law for this manufactured Cp/R=3.5+.001*T:
    // ln(m/m0)=2.5 ln(T/T0)+.001*(T-T0). The remaining gas is isentropic.
    let log_mass = (final_state.mass_kg / initial.mass_kg).ln();
    let (mut lo, mut hi) = (200.0_f64, initial.temperature_k);
    for _ in 0..80 {
        let t = lo.midpoint(hi);
        let residual = 2.5 * (t / initial.temperature_k).ln()
            + 0.001 * (t - initial.temperature_k) - log_mass;
        if residual > 0.0 { hi = t; } else { lo = t; }
    }
    (final_state.temperature_k - lo.midpoint(hi)).abs()
}

#[test]
fn sealed_supply_pressure_and_cooling_converge_to_the_variable_cp_tank_law() {
    let coarse = sealed_discharge(8e-6);
    let fine = sealed_discharge(4e-6);
    assert!(fine <= 0.4 * coarse + 2e-9, "coarse={coarse}, fine={fine}");
    assert!(fine < 5e-4, "temperature error={fine} K");
}

#[test]
fn bore_backflow_fills_and_heats_the_finite_supply_with_the_same_energy() {
    let gas = specimen();
    let dt = 2e-6;
    let mut field = bore(&gas, dt);
    let (mut supply, mut mouth) =
        FiniteMouth::new(&gas, spec(1e-5, 0.0, 95_000.0, 250.0), field.room()).unwrap();
    let before = supply.state(mouth);
    for _ in 0..12 {
        let report = advance(&gas, &mut field, &mut supply, &mut mouth, 5e-6, dt);
        assert!(report.bore_mass_kg < 0.0 && report.bore_energy_j < 0.0);
    }
    let after = supply.state(mouth);
    assert!(after.mass_kg > before.mass_kg);
    assert!(after.energy_j > before.energy_j);
    assert!(after.pressure_pa > before.pressure_pa);
    assert!(after.temperature_k > before.temperature_k);
}

#[test]
fn feed_reversal_uses_the_donors_actual_enthalpy_and_changes_inventory() {
    let gas = specimen();
    let dt = 2e-6;
    let mut field = bore(&gas, dt);
    let feed = stationary(&gas, 140_000.0, 293.15);
    let (mut supply, mut mouth) =
        FiniteMouth::new(&gas, spec(1e-5, 2e-6, 101_325.0, 293.15), feed).unwrap();
    let before = supply.state(mouth);
    for _ in 0..12 {
        let report = advance(&gas, &mut field, &mut supply, &mut mouth, 0.0, dt);
        assert!(report.feed_mass_kg > 0.0 && report.feed_energy_j > 0.0);
        let h = gas.stagnation_enthalpy(feed.temperature).unwrap();
        assert!((report.feed_energy_j / report.feed_mass_kg - h).abs() < 1e-10 * h);
    }
    let filled = supply.state(mouth);
    assert!(filled.pressure_pa > before.pressure_pa);
    assert!(filled.temperature_k > before.temperature_k);
    supply.set_feed(stationary(&gas, 90_000.0, 293.15));
    for _ in 0..12 {
        let report = advance(&gas, &mut field, &mut supply, &mut mouth, 0.0, dt);
        assert!(report.feed_mass_kg < 0.0 && report.feed_energy_j < 0.0);
        let delivered = report.feed_energy_j / report.feed_mass_kg;
        let lo = gas.stagnation_enthalpy(report.after.temperature_k).unwrap();
        let hi = gas.stagnation_enthalpy(report.before.temperature_k).unwrap();
        assert!(delivered >= lo - 1e-9 * lo && delivered <= hi + 1e-9 * hi);
    }
    assert!(supply.state(mouth).mass_kg < filled.mass_kg);
}

#[test]
fn repeated_predictors_and_refused_correctors_do_not_publish_storage() {
    let gas = specimen();
    let dt = 2e-6;
    let mut field = bore(&gas, dt);
    let (mut supply, mut mouth) =
        FiniteMouth::new(&gas, spec(2e-5, 1e-6, 106_325.0, 293.15), field.room()).unwrap();
    let before = supply.state(mouth);
    let cells = field.solver().states().to_vec();
    let clock = field.solver().time();
    let inlet = field.inlet_area_m2();
    let mut trial = supply.trial(mouth, dt, inlet, 0.0);
    let refused = field.step(|stage| {
        let flux = trial.flux(&gas, stage, |p| {
            Ok(lip_face(&gas, stage.left.1, p, 4e-6, inlet)?.flux)
        })?;
        if stage.stage == 0 {
            let replay = trial.flux(&gas, stage, |p| {
                Ok(lip_face(&gas, stage.left.1, p, 4e-6, inlet)?.flux)
            })?;
            assert_eq!(replay, flux);
        } else {
            return Err(Error::Invalid("test refusal after the complete compact-volume corrector"));
        }
        Ok(Some(flux))
    });
    assert!(refused.is_err());
    assert_eq!(field.solver().states(), cells);
    assert_eq!(field.solver().time(), clock);
    assert_eq!(supply.state(mouth), before);
    assert!(supply.report().is_none());
    let report = advance(&gas, &mut field, &mut supply, &mut mouth, 4e-6, dt);
    assert!(report.bore_mass_kg > 0.0);
}

#[test]
fn source_admission_and_compact_volume_cfl_refuse_before_publication() {
    let gas = specimen();
    let dt = 2e-6;
    let mut field = bore(&gas, dt);
    for bad in [
        spec(0.0, 0.0, 105_325.0, 293.15),
        spec(1e-5, -1.0, 105_325.0, 293.15),
        spec(1e-5, 0.0, f64::NAN, 293.15),
        spec(1e-5, 0.0, 105_325.0, 1200.0),
    ] {
        assert!(FiniteMouth::new(&gas, bad, field.room()).is_err());
    }
    let (supply, mouth) =
        FiniteMouth::new(&gas, spec(1e-12, 0.0, 105_325.0, 293.15), field.room()).unwrap();
    let before = supply.state(mouth);
    let cells = field.solver().states().to_vec();
    let inlet = field.inlet_area_m2();
    let mut trial = supply.trial(mouth, dt, inlet, 0.0);
    let refused = field.step(|stage| {
        let flux = trial.flux(&gas, stage, |p| {
            Ok(lip_face(&gas, stage.left.1, p, 4e-6, inlet)?.flux)
        })?;
        Ok(Some(flux))
    }).unwrap_err();
    assert!(format!("{refused}").contains("CflExceeded"));
    assert_eq!(field.solver().states(), cells);
    assert_eq!(supply.state(mouth), before);
    assert!(supply.report().is_none());
}

// ---- the series offset of the noise rows (the closure of their work on a finite mouth) ----

// Two identical runs, one with no offset and one with `series_pa`, side by side for `steps` steps.
fn paired_runs(volume: f64, feed: f64, pressure: f64, opening: f64, series_pa: f64, steps: usize)
    -> (Vec<(FiniteMouthReport, Vec<Primitive>)>, Vec<(FiniteMouthReport, Vec<Primitive>)>) {
    let gas = specimen();
    let dt = 2e-6;
    let run = |series: f64| {
        let mut field = bore(&gas, dt);
        let feed_state = stationary(&gas, 140_000.0, 293.15);
        let (mut supply, mut mouth) =
            FiniteMouth::new(&gas, spec(volume, feed, pressure, 293.15), feed_state).unwrap();
        (0..steps).map(|_| advance_with(&gas, &mut field, &mut supply, &mut mouth, opening, dt, series)).collect::<Vec<_>>()
    };
    (run(0.0), run(series_pa))
}

#[test]
fn without_an_offset_the_lip_is_offered_the_stored_gas_itself_and_no_work_is_reported() {
    let gas = specimen();
    let dt = 2e-6;
    let mut field = bore(&gas, dt);
    let (mut supply, mut mouth) =
        FiniteMouth::new(&gas, spec(2e-5, 1e-6, 106_325.0, 293.15), field.room()).unwrap();
    for _ in 0..20 {
        let before = mouth;
        let (report, offered) = advance_with(&gas, &mut field, &mut supply, &mut mouth, 4e-6, dt, 0.0);
        assert_eq!(offered.len(), 2, "one offered state per stage");
        // Stage 0 is the stored gas bit for bit; the offset path is not entered.
        assert_eq!(offered[0], before);
        assert_eq!(report.series_work_j.to_bits(), 0.0_f64.to_bits());
    }
}

#[test]
fn the_offset_reaches_both_stages_of_the_lip_face_at_the_stored_temperature() {
    let (plain, offset) = paired_runs(2e-5, 0.0, 106_325.0, 4e-6, 300.0, 12);
    for (a, b) in plain.iter().zip(&offset) {
        assert_eq!(a.1.len(), 2);
        assert_eq!(b.1.len(), 2);
        for stage in 0..2 {
            // Within the 12 steps the two runs' stored gas differs by far less than the offset, so each offered stage is the
            // plain run's stage plus the offset, at the stored temperature.
            let rise = b.1[stage].pressure - a.1[stage].pressure;
            assert!((rise - 300.0).abs() < 2.0, "stage {stage} was offered {rise} Pa more, not 300");
            assert!((b.1[stage].temperature - a.1[stage].temperature).abs() < 1e-3 * a.1[stage].temperature);
        }
    }
}

#[test]
fn the_feed_nozzle_and_the_inventory_never_see_the_offset() {
    // A feed above the stored gas: if the offset leaked into the feed's back pressure, the 5 kPa offset on a 38.7 kPa drop would
    // change the first step's feed mass by several percent. The two runs differ only through the (tiny) stored-gas change.
    let (plain, offset) = paired_runs(2e-5, 2e-6, 101_325.0, 2e-6, 5_000.0, 6);
    for (a, b) in plain.iter().zip(&offset) {
        let (fa, fb) = (a.0.feed_mass_kg, b.0.feed_mass_kg);
        assert!(fa > 0.0 && (fb - fa).abs() <= 2e-3 * fa, "feed mass {fa} against {fb}");
    }
}

#[test]
fn a_constant_offset_closes_the_audits_and_reports_the_work_it_would_do() {
    let (_, steps) = {
        let gas = specimen();
        let dt = 2e-6;
        let mut field = bore(&gas, dt);
        let (mut supply, mut mouth) =
            FiniteMouth::new(&gas, spec(2e-5, 0.0, 106_325.0, 293.15), field.room()).unwrap();
        let mut steps = Vec::new();
        for _ in 0..40 {
            // The harness audits mass and energy of the bore plus the stored gas to 2e-12 on every step.
            steps.push(advance_with(&gas, &mut field, &mut supply, &mut mouth, 4e-6, dt, 400.0));
        }
        (supply.state(mouth), steps)
    };
    let gas = specimen();
    let specific = gas.gas_constant() * 293.15;
    let (mut work, mut expected) = (0.0, 0.0);
    for (report, offered) in &steps {
        assert!(report.bore_mass_kg > 0.0);
        work += report.series_work_j;
        // Volume delivered to the bore at the density the lip was offered, independently from the bore's own mass receipt.
        let density = 0.5 * (offered[0].density + offered[1].density);
        expected += 400.0 * report.bore_mass_kg / density;
        assert!(report.energy_residual_j.abs() <= 2e-12 * report.before.energy_j.abs());
        assert!(report.mass_residual_kg.abs() <= 2e-12 * report.before.mass_kg);
    }
    assert!(work > 0.0 && (work - expected).abs() <= 2e-3 * expected, "work {work} J against {expected} J (RT {specific})");
}

#[test]
fn a_negative_offset_and_a_backflow_phase_keep_the_inventory_exact_and_the_work_signed() {
    let gas = specimen();
    let dt = 2e-6;
    // A stored gas below the room: with an opening the bore sucks from the room (backflow) unless the offset overcomes it.
    let mut field = bore(&gas, dt);
    let (mut supply, mut mouth) =
        FiniteMouth::new(&gas, spec(1e-5, 0.0, 95_000.0, 250.0), field.room()).unwrap();
    for _ in 0..8 {
        let (report, _) = advance_with(&gas, &mut field, &mut supply, &mut mouth, 5e-6, dt, 1_000.0);
        assert!(report.bore_mass_kg < 0.0, "the offset is too small to turn the flow round");
        assert!(report.series_work_j < 0.0, "a positive offset against a reversed flow does negative work");
    }
    let mut field = bore(&gas, dt);
    let (mut supply, mut mouth) =
        FiniteMouth::new(&gas, spec(2e-5, 0.0, 106_325.0, 293.15), field.room()).unwrap();
    for _ in 0..8 {
        let (report, _) = advance_with(&gas, &mut field, &mut supply, &mut mouth, 4e-6, dt, -300.0);
        assert!(report.bore_mass_kg > 0.0);
        assert!(report.series_work_j < 0.0, "a negative offset on a forward flow does negative work");
    }
}

#[test]
fn an_offset_that_leaves_no_positive_pressure_is_refused_before_publication() {
    let gas = specimen();
    let dt = 2e-6;
    let mut field = bore(&gas, dt);
    let (supply, mouth) =
        FiniteMouth::new(&gas, spec(2e-5, 0.0, 106_325.0, 293.15), field.room()).unwrap();
    let before = supply.state(mouth);
    let cells = field.solver().states().to_vec();
    let inlet = field.inlet_area_m2();
    let mut trial = supply.trial(mouth, dt, inlet, -200_000.0);
    let refused = field.step(|stage| {
        let flux = trial.flux(&gas, stage, |p| Ok(lip_face(&gas, stage.left.1, p, 4e-6, inlet)?.flux))?;
        Ok(Some(flux))
    });
    assert!(refused.is_err());
    assert_eq!(field.solver().states(), cells);
    assert_eq!(supply.state(mouth), before);
    assert!(supply.report().is_none());
}
