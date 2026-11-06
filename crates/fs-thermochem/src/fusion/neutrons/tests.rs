use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g0_channel_identity_number_and_energy_survive_the_handoff() {
    let amounts = [Amount::new(2e-6), Amount::new(3e-6)];
    let per_particle = [14.1e6 * 1.602_176_634e-19, 2.45e6 * 1.602_176_634e-19];
    let energies = core::array::from_fn(|i| {
        Energy::new(amounts[i].value() * 6.022_140_76e23 * per_particle[i])
    });
    let emission = FusionNeutronEmissionV1::new(Time::new(0.25), amounts, energies).unwrap();
    assert_eq!(emission.duration(), Time::new(0.25));
    let reactions = [
        FusionReactionV1::DeuteriumTritium,
        FusionReactionV1::DeuteriumDeuteriumNeutron,
    ];
    for (i, channel) in emission.channels().iter().enumerate() {
        assert_eq!(channel.reaction(), reactions[i]);
        assert_eq!(channel.amount(), amounts[i]);
        assert_eq!(channel.kinetic_energy(), energies[i]);
        close(
            channel.particle_count(),
            amounts[i].value() * 6.022_140_76e23,
        );
        close(
            channel.mean_particle_energy().unwrap().value(),
            per_particle[i],
        );
    }
    close(emission.total_amount().value(), 5e-6);
    close(emission.total_particle_count(), 5e-6 * 6.022_140_76e23);
    assert_eq!(emission.total_kinetic_energy(), energies[0] + energies[1]);
}

#[test]
fn g0_empty_channels_have_no_mean_or_invented_particles() {
    for duration in [0.0, 3.0] {
        let emission = FusionNeutronEmissionV1::new(
            Time::new(duration),
            [Amount::new(0.0); 2],
            [Energy::new(0.0); 2],
        )
        .unwrap();
        assert_eq!(emission.total_amount().value(), 0.0);
        assert_eq!(emission.total_kinetic_energy().value(), 0.0);
        assert_eq!(emission.total_particle_count(), 0.0);
        for channel in emission.channels() {
            assert_eq!(channel.particle_count(), 0.0);
            assert_eq!(channel.mean_particle_energy(), None);
        }
    }
    let only_dd = FusionNeutronEmissionV1::new(
        Time::new(1.0),
        [Amount::new(0.0), Amount::new(1e-6)],
        [Energy::new(0.0), Energy::new(4.0)],
    )
    .unwrap();
    assert_eq!(only_dd.channels()[0].mean_particle_energy(), None);
    assert!(only_dd.channels()[1].mean_particle_energy().is_some());
}

#[test]
fn g0_invalid_duration_population_and_energy_are_refused() {
    use NeutronEmissionErrorV1 as E;
    for duration in [-1.0, f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
        assert_eq!(
            FusionNeutronEmissionV1::new(
                Time::new(duration),
                [Amount::new(0.0); 2],
                [Energy::new(0.0); 2],
            ),
            Err(E::InvalidDuration)
        );
    }
    assert_eq!(
        FusionNeutronEmissionV1::new(
            Time::new(0.0),
            [Amount::new(1.0), Amount::new(0.0)],
            [Energy::new(1.0), Energy::new(0.0)],
        ),
        Err(E::InvalidDuration)
    );
    for (i, reaction) in [
        FusionReactionV1::DeuteriumTritium,
        FusionReactionV1::DeuteriumDeuteriumNeutron,
    ]
    .into_iter()
    .enumerate()
    {
        for invalid in [-1.0, f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
            let mut amounts = [Amount::new(1.0); 2];
            amounts[i] = Amount::new(invalid);
            assert_eq!(
                FusionNeutronEmissionV1::new(Time::new(1.0), amounts, [Energy::new(1.0); 2]),
                Err(E::InvalidAmount(reaction))
            );
            let mut energies = [Energy::new(1.0); 2];
            energies[i] = Energy::new(invalid);
            assert_eq!(
                FusionNeutronEmissionV1::new(Time::new(1.0), [Amount::new(1.0); 2], energies),
                Err(E::InvalidEnergy(reaction))
            );
        }
        for (amount, energy) in [(0.0, 1.0), (1.0, 0.0)] {
            let mut amounts = [Amount::new(1.0); 2];
            let mut energies = [Energy::new(1.0); 2];
            amounts[i] = Amount::new(amount);
            energies[i] = Energy::new(energy);
            assert_eq!(
                FusionNeutronEmissionV1::new(Time::new(1.0), amounts, energies),
                Err(E::InconsistentChannel(reaction))
            );
        }
    }
}

#[test]
fn g0_unrepresentable_counts_means_and_combined_accounts_refuse() {
    let smallest = f64::from_bits(1);
    for (amount, energy) in [(f64::MAX, 1.0), (smallest, f64::MAX), (1.0, smallest)] {
        assert!(matches!(
            FusionNeutronEmissionV1::new(
                Time::new(1.0),
                [Amount::new(amount), Amount::new(0.0)],
                [Energy::new(energy), Energy::new(0.0)],
            ),
            Err(NeutronEmissionErrorV1::Arithmetic(_))
        ));
    }
    for (amounts, energies) in [
        ([1.0; 2], [0.75 * f64::MAX; 2]),
        ([0.75 * (f64::MAX / 6.022_140_76e23); 2], [1.0; 2]),
    ] {
        assert!(matches!(
            FusionNeutronEmissionV1::new(
                Time::new(1.0),
                amounts.map(Amount::new),
                energies.map(Energy::new),
            ),
            Err(NeutronEmissionErrorV1::Arithmetic(_))
        ));
    }
    // Tiny, representable inventories remain valid; subnormal amount is not
    // itself an underflow. The mean still lies within finite positive range.
    let tiny = FusionNeutronEmissionV1::new(
        Time::new(1.0),
        [Amount::new(smallest), Amount::new(0.0)],
        [Energy::new(smallest), Energy::new(0.0)],
    )
    .unwrap();
    assert!(tiny.channels()[0].particle_count() > 0.0);
    assert!(tiny.channels()[0].mean_particle_energy().unwrap().value() > 0.0);
}

#[test]
fn g3_scaling_integrated_accounts_retains_each_channel_mean() {
    let original = FusionNeutronEmissionV1::new(
        Time::new(2.0),
        [Amount::new(1e-6), Amount::new(4e-6)],
        [Energy::new(8.0), Energy::new(3.0)],
    )
    .unwrap();
    for factor in [1e-12, 3.0, 1e12] {
        let scaled = FusionNeutronEmissionV1::new(
            Time::new(7.0),
            original
                .channels()
                .map(|c| Amount::new(c.amount().value() * factor)),
            original
                .channels()
                .map(|c| Energy::new(c.kinetic_energy().value() * factor)),
        )
        .unwrap();
        for (a, b) in original.channels().iter().zip(scaled.channels()) {
            close(
                a.mean_particle_energy().unwrap().value(),
                b.mean_particle_energy().unwrap().value(),
            );
        }
        close(
            scaled.total_amount().value(),
            original.total_amount().value() * factor,
        );
        close(
            scaled.total_kinetic_energy().value(),
            original.total_kinetic_energy().value() * factor,
        );
    }
}

#[test]
fn g0_interval_preserves_original_clock_and_complete_channel_accounts() {
    let start = Time::new(0.2);
    let end = Time::new(0.6);
    let emission = FusionNeutronEmissionV1::new(
        Time::new(end.value() - start.value()),
        [Amount::new(2e-6), Amount::new(3e-6)],
        [Energy::new(5.0), Energy::new(7.0)],
    )
    .unwrap();
    let interval = FusionNeutronIntervalV1::new(start, end, emission).unwrap();
    let copy = interval;
    assert_eq!(interval, copy);
    assert_eq!(interval.start(), start);
    assert_eq!(interval.end(), end);
    assert_eq!(interval.emission(), &emission);
}

#[test]
fn g0_interval_retains_empty_emission_only_on_a_positive_matching_window() {
    let emission =
        FusionNeutronEmissionV1::new(Time::new(0.5), [Amount::new(0.0); 2], [Energy::new(0.0); 2])
            .unwrap();
    let interval = FusionNeutronIntervalV1::new(Time::new(2.0), Time::new(2.5), emission).unwrap();
    assert_eq!(interval.emission(), &emission);
    assert_eq!(interval.emission().total_particle_count(), 0.0);
    let zero_duration =
        FusionNeutronEmissionV1::new(Time::new(0.0), [Amount::new(0.0); 2], [Energy::new(0.0); 2])
            .unwrap();
    assert_eq!(
        FusionNeutronIntervalV1::new(Time::new(2.0), Time::new(2.5), zero_duration),
        Err(NeutronEmissionErrorV1::InconsistentIntervalDuration)
    );
}

#[test]
fn g3_interval_rejects_invalid_clock_windows_without_retiming_emission() {
    let emission = FusionNeutronEmissionV1::new(
        Time::new(0.5),
        [Amount::new(1e-6), Amount::new(0.0)],
        [Energy::new(2.0), Energy::new(0.0)],
    )
    .unwrap();
    for (start, end) in [
        (-1.0, -0.5),
        (-0.5, 0.0),
        (1.0, 1.0),
        (1.0, 0.5),
        (f64::NAN, 0.5),
        (0.0, f64::NAN),
        (f64::INFINITY, f64::INFINITY),
        (0.0, f64::INFINITY),
        (f64::NEG_INFINITY, 0.5),
    ] {
        assert_eq!(
            FusionNeutronIntervalV1::new(Time::new(start), Time::new(end), emission),
            Err(NeutronEmissionErrorV1::InvalidInterval)
        );
    }
    assert_eq!(
        FusionNeutronIntervalV1::new(Time::new(1.0), Time::new(2.0), emission),
        Err(NeutronEmissionErrorV1::InconsistentIntervalDuration)
    );
}

#[test]
fn g3_interval_requires_exact_duration_at_ordinary_and_large_clock_offsets() {
    for (start, end, authored_duration) in [
        (0.2, 0.6, 0.4),
        (2.0_f64.powi(52), 2.0_f64.powi(52) + 0.75, 0.75),
    ] {
        assert_ne!(end - start, authored_duration);
        let emission = FusionNeutronEmissionV1::new(
            Time::new(authored_duration),
            [Amount::new(1e-6), Amount::new(0.0)],
            [Energy::new(2.0), Energy::new(0.0)],
        )
        .unwrap();
        assert_eq!(
            FusionNeutronIntervalV1::new(Time::new(start), Time::new(end), emission),
            Err(NeutronEmissionErrorV1::InconsistentIntervalDuration)
        );
        let matching = FusionNeutronEmissionV1::new(
            Time::new(end - start),
            emission.channels().map(|channel| channel.amount()),
            emission.channels().map(|channel| channel.kinetic_energy()),
        )
        .unwrap();
        let interval =
            FusionNeutronIntervalV1::new(Time::new(start), Time::new(end), matching).unwrap();
        assert_eq!(interval.emission(), &matching);
    }
}
