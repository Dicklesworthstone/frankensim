//! Accepted source histories retain native moments without changing the burn.

use super::*;
use fs_thermochem::fusion::neutrons::{FusionNeutronEmissionV1, FusionNeutronIntervalV1};

pub(in crate::fusion) fn check_history(
    history: Option<&[FusionNeutronIntervalV1]>,
    accepted: usize,
    total: FusionNeutronEmissionV1,
) {
    let history = history.expect("requested source history");
    assert_eq!(history.len(), accepted);
    let mut amounts = [solve::Sum::default(); 2];
    let mut energies = [solve::Sum::default(); 2];
    let mut end = Time::new(0.0);
    for interval in history {
        assert_eq!(interval.start(), end);
        assert!(interval.end().value() > interval.start().value());
        assert_eq!(
            interval.emission().duration().value(),
            interval.end().value() - interval.start().value(),
        );
        for (r, channel) in interval.emission().channels().iter().enumerate() {
            assert_eq!(channel.reaction(), total.channels()[r].reaction());
            amounts[r].add(channel.amount().value());
            energies[r].add(channel.kinetic_energy().value());
        }
        end = interval.end();
    }
    assert_eq!(end, total.duration());
    for r in 0..2 {
        close(
            amounts[r].value,
            total.channels()[r].amount().value(),
            3e-14,
            0.0,
        );
        close(
            energies[r].value,
            total.channels()[r].kinetic_energy().value(),
            3e-14,
            0.0,
        );
    }
}

#[test]
fn g3_dt_history_follows_changing_source_excludes_retries_and_preserves_physics() {
    let (plasma, initial) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let original = initial.clone();
    let transport = DtFusionTransport {
        energy_confinement: Some(
            FusionEnergyConfinement::new(
                Time::new(0.02),
                Time::new(0.03),
                Temperature::new(2.0 * KELVIN_PER_KEV),
                Temperature::new(KELVIN_PER_KEV),
            )
            .unwrap(),
        ),
        ..DtFusionTransport::default()
    };
    let mut c = controls();
    c.relative_tolerance = 1e-7;
    c.initial_step_s = 0.02;
    c.maximum_step_s = 0.02;
    assert!(!c.record_neutron_history);
    let plain = advance_dt_fusion_with_transport(
        &initial,
        &plasma,
        heating(0.8),
        &transport,
        Time::new(0.02),
        c,
        || false,
    )
    .unwrap();
    assert!(plain.neutron_history().is_none());
    c.record_neutron_history = true;
    let mut polls = 0;
    let recorded = advance_dt_fusion_with_transport(
        &initial,
        &plasma,
        heating(0.8),
        &transport,
        Time::new(0.02),
        c,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(recorded.state(), plain.state());
    assert_eq!(recorded.balance(), plain.balance());
    assert_eq!(recorded.neutron_emission(), plain.neutron_emission());
    assert!(recorded.balance().rejected_steps > 0);
    check_history(
        recorded.neutron_history(),
        recorded.balance().accepted_steps,
        recorded.neutron_emission().unwrap(),
    );
    let intervals = recorded.neutron_history().unwrap();
    let rate = |interval: &FusionNeutronIntervalV1| {
        interval.emission().total_particle_count() / interval.emission().duration().value()
    };
    assert!(rate(intervals.last().unwrap()) < 0.8 * rate(&intervals[0]));
    // Refuse on the final cancellation boundary, after the accepted history
    // has been built. Neither history nor a partial physical result escapes.
    let mut replay_polls = 0;
    assert_eq!(
        advance_dt_fusion_with_transport(
            &initial,
            &plasma,
            heating(0.8),
            &transport,
            Time::new(0.02),
            c,
            || {
                replay_polls += 1;
                replay_polls == polls
            },
        ),
        Err(DtFusionError::Cancelled)
    );
    assert_eq!(initial, original);
}

#[test]
fn g0_history_distinguishes_disabled_empty_and_zero_emission_intervals() {
    let (plasma, initial) = source::plasma([0.0, 0.0, 5e19], 1.0, [10.0, 5.0]);
    let mut c = controls();
    let disabled =
        advance_dt_fusion(&initial, &plasma, heating(0.8), Time::new(0.0), c, || false).unwrap();
    assert!(disabled.neutron_history().is_none());
    c.record_neutron_history = true;
    for duration in [0.0, 0.01] {
        let result = advance_dt_fusion(
            &initial,
            &plasma,
            heating(0.8),
            Time::new(duration),
            c,
            || false,
        )
        .unwrap();
        check_history(
            result.neutron_history(),
            result.balance().accepted_steps,
            result.neutron_emission().unwrap(),
        );
        assert_eq!(result.state().amounts(), initial.amounts());
        close(
            result.state().heavy_temperature().value(),
            initial.heavy_temperature().value(),
            3e-12,
            1e-5,
        );
        close(
            result.state().electron_temperature().value(),
            initial.electron_temperature().value(),
            3e-12,
            1e-5,
        );
        assert_eq!(
            result.neutron_history().unwrap().is_empty(),
            duration == 0.0
        );
        assert!(
            result
                .neutron_history()
                .unwrap()
                .iter()
                .all(|interval| interval.emission().total_particle_count() == 0.0)
        );
    }
    c.maximum_steps = 1;
    assert_eq!(
        advance_dt_fusion(&initial, &plasma, heating(0.8), Time::new(0.01), c, || {
            false
        },),
        Err(DtFusionError::StepBudget)
    );
}
