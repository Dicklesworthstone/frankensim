//! Both neutron channels follow accepted burn intervals and original knots.

use super::*;
use crate::fusion::tests::history::check_history;
use crate::fusion::volume::FusionVolumeKnot;
use fs_qty::Volume;

fn geometry(knots: &[(f64, f64)]) -> FusionVolumeHistory {
    FusionVolumeHistory::new(
        knots
            .iter()
            .map(|&(time, volume)| FusionVolumeKnot {
                time: Time::new(time),
                volume: Volume::new(volume),
            })
            .collect(),
        || false,
    )
    .unwrap()
}

#[test]
fn g3_network_history_keeps_both_channels_and_every_volume_boundary() {
    let (plasma, initial) = fixture::plasma(MIXED_DENSITY, 1.0, [15.0, 4.0]);
    let volume = geometry(&[(0.0, 1.0), (0.003, 0.9), (0.01, 1.05)]);
    let transport = FusionNetworkTransport::default();
    let mut c = controls();
    c.relative_tolerance = 1e-6;
    let plain = advance_fusion_network_with_volume(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        &volume,
        c,
        || false,
    )
    .unwrap();
    c.record_neutron_history = true;
    let recorded = advance_fusion_network_with_volume(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        &volume,
        c,
        || false,
    )
    .unwrap();
    assert!(plain.neutron_history().is_none());
    assert_eq!(recorded.state(), plain.state());
    assert_eq!(recorded.balance(), plain.balance());
    let emission = recorded.neutron_emission().unwrap();
    assert!(
        emission
            .channels()
            .iter()
            .all(|channel| channel.particle_count() > 0.0)
    );
    check_history(
        recorded.neutron_history(),
        recorded.balance().accepted_steps,
        emission,
    );
    for knot in &volume.knots()[1..] {
        assert!(
            recorded
                .neutron_history()
                .unwrap()
                .iter()
                .any(|interval| interval.end() == knot.time)
        );
    }
}

#[test]
fn g1_history_retains_a_late_source_smaller_than_cumulative_roundoff() {
    let (plasma, initial) = fixture::plasma([0.0, 5e19, 5e19, 0.0, 0.0], 1.0, [10.0, 5.0]);
    // An early hot burn is followed by adiabatic cooling/dilution and a short
    // authored final hold. The last neutron contribution remains representable
    // even though adding it to the full-run total cannot change that total.
    let volume = geometry(&[
        (0.0, 1.0),
        (0.001, 1.0),
        (0.002, 64.0),
        (0.002_000_000_000_001, 64.0),
    ]);
    let mut c = controls();
    c.relative_tolerance = 1e-6;
    c.record_neutron_history = true;
    let result = advance_fusion_network_with_volume(
        &initial,
        &plasma,
        FusionNetworkChannels::new([true, false, false, false]),
        heating(),
        &FusionNetworkTransport::default(),
        &volume,
        c,
        || false,
    )
    .unwrap();
    let total = result.neutron_emission().unwrap();
    check_history(
        result.neutron_history(),
        result.balance().accepted_steps,
        total,
    );
    let last = result.neutron_history().unwrap().last().unwrap();
    assert_eq!(last.start(), Time::new(0.002));
    let amount = last.emission().total_amount().value();
    let energy = last.emission().total_kinetic_energy().value();
    assert!(amount > 0.0 && energy > 0.0);
    assert_eq!(
        total.total_amount().value() + amount,
        total.total_amount().value()
    );
    assert_eq!(
        total.total_kinetic_energy().value() + energy,
        total.total_kinetic_energy().value()
    );
}
