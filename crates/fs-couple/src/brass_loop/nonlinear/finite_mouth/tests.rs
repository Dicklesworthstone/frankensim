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
    let before_field = field.solver().totals();
    let before_mouth = supply.state(*mouth);
    let inlet = field.inlet_area_m2();
    let mut trial = supply.trial(*mouth, dt, inlet);
    let step = field.step(|stage| {
        let flux = trial.flux(gas, stage, |p| {
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
    report
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
    let mut trial = supply.trial(mouth, dt, inlet);
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
    let mut trial = supply.trial(mouth, dt, inlet);
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
