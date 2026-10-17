use super::*;
use crate::cavity_phs::piston::tests::{gas, controls, drive};
use std::convert::Infallible;

fn feedback_controls() -> HeatExchangeSolve {
    HeatExchangeSolve { heat_bracket_j: [-1.0, 1.0], source_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-7, maximum_iterations: 128 }
}
fn piston() -> GasPiston { GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap() }
fn law(port: ThermalPort, gas: &ResolvedThermochemicalGasState)
    -> Result<ConductanceSample<(u64, f64)>, Infallible>
{
    Ok(ConductanceSample { conductance_w_k: gas.state().thermal_conductivity * 20.0,
        source: (port.id, gas.state().temperature) })
}

#[test]
fn one_contact_reuses_single_boundary_and_empty_reuses_isolated_motion() {
    let old = piston();
    let port = ThermalPort { id: 7, temperature_k: 500.0 };
    let many = old.thermal_ports_candidate(drive(0.02), controls(), feedback_controls(),
        &[port], 1, law, || false).unwrap();
    let single = old.thermal_boundary_candidate(drive(0.02), controls(), feedback_controls(),
        port.temperature_k, |gas| law(port, gas), || false).unwrap();
    assert_eq!(many.piston(), single.piston());
    assert_eq!(many.step().piston, single.step().piston);
    assert_eq!(many.step().evaluated_heat_j, single.step().evaluated_heat_j);
    assert_eq!(many.step().exchange_trials, single.step().exchange_trials);
    assert_eq!(&many.source()[0].source, single.source());
    let empty = old.thermal_ports_candidate(drive(0.02), controls(), feedback_controls(),
        &[], 0, law, || false).unwrap();
    let mut isolated = old.clone();
    let report = isolated.step(drive(0.02), controls(), || false).unwrap();
    assert_eq!(empty.piston(), &isolated);
    assert_eq!(empty.step().piston, report);
    assert!(empty.source().is_empty());
}

#[test]
fn equal_opposite_heat_survives_zero_net_and_every_source_is_sampled() {
    let old = piston();
    let ports = [ThermalPort { id: 9, temperature_k: 500.0 },
        ThermalPort { id: 3, temperature_k: 300.0 }];
    let solve = |ports: &[ThermalPort]| old.thermal_ports_candidate(drive(0.0), controls(), feedback_controls(),
        ports, 2, |p, _| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 2.0, source: p.id }),
        || false).unwrap();
    let result = solve(&ports);
    assert_eq!(result.piston(), &old);
    assert_eq!(result.step().exchanged_heat_j, 0.0);
    assert_eq!(result.source()[0].port.id, 3);
    // 2 W/K across the -100 K and +100 K port-to-gas differences, one step.
    let dt = drive(0.0).duration_s;
    assert!((result.source()[1].heat_j - 2.0 * 100.0 * dt).abs() <= 1e-15);
    assert_eq!(result.source()[0].heat_j, -result.source()[1].heat_j);
    assert_eq!(result, solve(&[ports[1], ports[0]]));
    let missing = old.thermal_ports_candidate(drive(0.0), controls(), feedback_controls(), &ports, 2,
        |_, _| Err::<ConductanceSample<()>, _>("missing conductivity"), || false);
    assert!(matches!(missing, Err(ThermalPistonError::Interface("missing conductivity"))));
}

#[test]
fn material_heat_split_uses_one_actual_midpoint_and_changes_motion() {
    let old = piston();
    let ports = [ThermalPort { id: 1, temperature_k: 500.0 }, ThermalPort { id: 2, temperature_k: 350.0 }];
    let run = |factor| old.thermal_ports_candidate(drive(0.02), controls(), feedback_controls(), &ports, 2,
        |p, gas| {
            let mut sample = law(p, gas).unwrap();
            if p.id == 1 { sample.conductance_w_k *= factor; }
            Ok::<_, Infallible>(sample)
        }, || false).unwrap();
    let a = run(1.0); let b = run(3.0);
    assert!(b.piston().velocity_m_s() > a.piston().velocity_m_s());
    assert!(b.piston().gas().state().state().sound_speed > a.piston().gas().state().state().sound_speed);
    for result in [&a, &b] {
        let midpoint = result.step().piston.midpoint.state().temperature;
        for source in result.source() {
            assert_eq!(source.source.1, midpoint);
            assert_eq!(source.heat_j, drive(0.0).duration_s * (source.conductance_w_k * (source.port.temperature_k - midpoint)));
        }
        assert!(result.source()[0].heat_j > 0.0 && result.source()[1].heat_j < 0.0);
        assert!(result.step().source_residual_j.abs() <= feedback_controls().source_tolerance_j);
        assert!(result.step().total_energy_residual_j.abs() <= feedback_controls().energy_tolerance_j);
    }
}

#[test]
fn bad_contacts_budget_and_late_cancellation_return_no_endpoint() {
    let old = piston(); let saved = old.clone();
    let port = ThermalPort { id: 7, temperature_k: 500.0 };
    assert!(old.thermal_ports_candidate(drive(0.0), controls(), feedback_controls(), &[port], 0, law, || false).is_err());
    assert!(old.thermal_ports_candidate(drive(0.0), controls(), feedback_controls(), &[port,port], 2, law, || false).is_err());
    for temperature_k in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(old.thermal_ports_candidate(drive(0.0), controls(), feedback_controls(),
            &[ThermalPort { temperature_k, ..port }], 1, law, || false).is_err());
    }
    let mut count = 0;
    old.thermal_ports_candidate(drive(0.0), controls(), feedback_controls(), &[port], 1, law,
        || { count += 1; false }).unwrap();
    let mut actual = 0;
    assert!(matches!(old.thermal_ports_candidate(drive(0.0), controls(), feedback_controls(), &[port], 1, law,
        || { actual += 1; actual == count }), Err(ThermalPistonError::Cancelled)));
    assert_eq!(old, saved);
}
