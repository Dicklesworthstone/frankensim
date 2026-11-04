use super::*;
use fs_material::healing::kinetics::{MobilityIntegration, ShiftedMobility};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn kinetics(upper_k: f64) -> ShiftedMobility {
    let shift = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 2.0,
            c2_k: 100.0,
        },
        400.0,
        [330.0, upper_k],
    )
    .unwrap();
    ShiftedMobility::new(
        "test-mobility".into(),
        shift,
        0.25,
        Some(350.0),
        MobilityIntegration {
            absolute_dose_tolerance: 1e-12,
            relative_dose_tolerance: 1e-10,
            maximum_evaluations: 4097,
            maximum_depth: 24,
        },
    )
    .unwrap()
}

fn heated_schedule(temperature: f64) -> Vec<TimedEvent> {
    let mut events = schedule(temperature, 0.0);
    for id in [CellId(1), CellId(2)] {
        events.push(event(0.0, Event::Boundary(id, boundary(0.0, 400.0))));
    }
    events
}

#[test]
fn shifted_kinetics_uses_the_actual_thermal_ramp_and_develops_the_original_bond() {
    let law = kinetics(600.0);
    let initial = simulation(heated_schedule(350.0), &law, &[7], 1000);
    // Each original cell has C=2 J/K and receives 400 W. Their identical
    // temperatures rise linearly from 350 to 450 K in 0.5 s. The native
    // observer supplies that path to kinetics without an intermediate table.
    let expected = 0.183_067_380_238_716_9;
    for dt in [0.5, 0.125] {
        let (result, audit, _) = initial.advance_to(0.5, config(dt), || Ok(())).unwrap();
        close(audit.residual_j, 0.0, 1e-9);
        close(
            result.thermal().temperature_k(CellId(1)).unwrap(),
            450.0,
            1e-12,
        );
        close(
            result.thermal().temperature_k(CellId(2)).unwrap(),
            450.0,
            1e-12,
        );
        close(result.records()[&7].history().dose(), expected, 2e-11);
        close(result.records()[&7].history().contact_time_s(), 0.5, 1e-14);
        let frozen = result.freeze_bond(7, &map()).unwrap();
        let cohesive = frozen.bond.law.unwrap();
        close(cohesive.strength_pa(), 2e4 * expected, 5e-7);
        close(cohesive.fracture_energy_j_m2(), 10.0 * expected, 3e-10);
        // Exercise the actual native cohesive response created by that thermal
        // history, including peak force and the work to separate its real area.
        let peak = cohesive
            .bonded_point()
            .trial([cohesive.onset_separation_m(), 0.0, 0.0], 1e-10)
            .unwrap();
        close(
            frozen.area_m2 * peak.traction_pa[0],
            frozen.area_m2 * 2e4 * expected,
            1e-10,
        );
        let failed = peak
            .point
            .trial([cohesive.failure_separation_m(), 0.0, 0.0], 1e-10)
            .unwrap();
        assert!(failed.point.failed());
        assert_eq!(failed.traction_pa, [0.0; 3]);
        close(
            frozen.area_m2 * (peak.energy.work + failed.energy.work),
            frozen.area_m2 * 10.0 * expected,
            7e-14,
        );
        close(
            frozen.area_m2 * failed.point.dissipated_j_m2(),
            frozen.area_m2 * 10.0 * expected,
            7e-14,
        );
        let reclosed = failed.point.trial([-1e-5, 0.0, 0.0], 1e-10).unwrap();
        assert!(reclosed.traction_pa[0] < 0.0);
        assert!(reclosed.point.failed());
        assert_eq!(
            reclosed.point.dissipated_j_m2(),
            failed.point.dissipated_j_m2()
        );
        let reopened = reclosed
            .point
            .trial([cohesive.onset_separation_m(), 0.0, 0.0], 1e-10)
            .unwrap();
        assert_eq!(reopened.traction_pa, [0.0; 3]);
        assert!(reopened.point.failed());
    }
    assert!(initial.records().is_empty());
}

#[test]
fn shifted_kinetics_domain_refusal_preserves_both_original_owners() {
    let law = kinetics(450.0);
    let initial = simulation(heated_schedule(440.0), &law, &[7], 1000);
    let (accepted, _, _) = initial.advance_to(0.0, config(0.1), || Ok(())).unwrap();
    // The caloric law admits 460 K; only the original kinetic calibration ends
    // at 450 K. Its refusal must also discard the proposed thermal endpoint.
    assert!(matches!(
        accepted.advance_to(0.1, config(0.1), || Ok(())),
        Err(BondingError::Healing {
            error: fs_material::healing::HealingError::Temperature(_),
            ..
        })
    ));
    assert_eq!(accepted.thermal().time_s(), 0.0);
    assert_eq!(accepted.thermal().temperature_k(CellId(1)).unwrap(), 440.0);
    assert_eq!(accepted.records()[&7].history().dose(), 0.0);
    assert_eq!(accepted.records()[&7].history().contact_time_s(), 0.0);
}
