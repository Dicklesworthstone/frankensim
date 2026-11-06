use super::*;

#[test]
fn g1_each_donor_volume_sets_its_own_well_mixed_residence_half_life() {
    let volumes = [Volume::new(2.0), Volume::new(8.0)];
    let flow = VolumetricFlowRate::new(1.0);
    let forward = NuclideTransfer::from_volumetric_flow(0, 1, &volumes, flow, 3).unwrap();
    let reverse = NuclideTransfer::from_volumetric_flow(1, 0, &volumes, flow, 3).unwrap();
    let extraction = NuclideOutflow::from_volumetric_flow(1, &volumes, flow, 3).unwrap();
    assert_eq!(forward.rates_per_s, [0.5; 3]);
    assert_eq!(reverse.rates_per_s, [0.125; 3]);
    assert_eq!(extraction.rates_per_s, reverse.rates_per_s);
    // A fixed 1 m3/s flow flushes half the tracer in V*ln(2) seconds.
    for (rate, half_time) in [(forward.rates_per_s[0], 2.0), (reverse.rates_per_s[0], 8.0)] {
        let remaining =
            fs_math::decay::remaining_fraction(rate * half_time * core::f64::consts::LN_2);
        assert!((remaining - 0.5).abs() < 1e-15);
    }
}

#[test]
fn g3_scaling_carrier_flow_and_volumes_preserves_the_same_tracer_operator() {
    let original = [Volume::new(2.0), Volume::new(8.0)];
    let scaled = [Volume::new(16.0), Volume::new(64.0)];
    let transfer =
        NuclideTransfer::from_volumetric_flow(0, 1, &original, VolumetricFlowRate::new(1.0), 2)
            .unwrap();
    let larger =
        NuclideTransfer::from_volumetric_flow(0, 1, &scaled, VolumetricFlowRate::new(8.0), 2)
            .unwrap();
    assert_eq!(transfer, larger);
    let outflow =
        NuclideOutflow::from_volumetric_flow(1, &original, VolumetricFlowRate::new(1.0), 2)
            .unwrap();
    let larger =
        NuclideOutflow::from_volumetric_flow(1, &scaled, VolumetricFlowRate::new(8.0), 2).unwrap();
    assert_eq!(outflow, larger);
    // No intermediate product overflows when the finite quotient is one.
    let extremes = [Volume::new(f64::MAX), Volume::new(1.0)];
    let extreme = NuclideTransfer::from_volumetric_flow(
        0,
        1,
        &extremes,
        VolumetricFlowRate::new(f64::MAX),
        2,
    )
    .unwrap();
    assert_eq!(extreme.rates_per_s, [1.0; 2]);
}

#[test]
fn g0_zero_and_resolved_small_rates_are_preserved_without_clamping_or_repair() {
    let volumes = [Volume::new(1.0), Volume::new(2.0)];
    let closed = NuclideTransfer::from_volumetric_flow(
        0,
        1,
        &volumes,
        VolumetricFlowRate::new(0.0),
        MAX_DECAY_NODES,
    )
    .unwrap();
    assert_eq!(closed.from_compartment, 0);
    assert_eq!(closed.to_compartment, 1);
    assert_eq!(closed.rates_per_s, vec![0.0; MAX_DECAY_NODES]);
    let stopped =
        NuclideOutflow::from_volumetric_flow(1, &volumes, VolumetricFlowRate::new(0.0), 1).unwrap();
    assert_eq!(stopped.rates_per_s, [0.0]);

    let smallest = f64::from_bits(1);
    let tiny =
        NuclideTransfer::from_volumetric_flow(0, 1, &volumes, VolumetricFlowRate::new(smallest), 1)
            .unwrap();
    assert_eq!(tiny.rates_per_s, [smallest]);
    assert_eq!(
        NuclideOutflow::from_volumetric_flow(1, &volumes, VolumetricFlowRate::new(smallest), 1),
        Err(FissionError::Arithmetic)
    );
    let tiny_volume = [Volume::new(0.5)];
    assert_eq!(
        NuclideOutflow::from_volumetric_flow(0, &tiny_volume, VolumetricFlowRate::new(f64::MAX), 1),
        Err(FissionError::Arithmetic)
    );
}

#[test]
fn g0_both_referenced_volumes_and_directed_path_dimensions_are_required() {
    let valid = [Volume::new(1.0), Volume::new(2.0)];
    for (from, to) in [(0, 0), (2, 0), (0, 2)] {
        assert!(
            NuclideTransfer::from_volumetric_flow(
                from,
                to,
                &valid,
                VolumetricFlowRate::new(1.0),
                2
            )
            .is_err()
        );
    }
    assert!(
        NuclideOutflow::from_volumetric_flow(2, &valid, VolumetricFlowRate::new(1.0), 2).is_err()
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        for volumes in [[Volume::new(bad), valid[1]], [valid[0], Volume::new(bad)]] {
            assert!(
                NuclideTransfer::from_volumetric_flow(
                    0,
                    1,
                    &volumes,
                    VolumetricFlowRate::new(0.0),
                    2
                )
                .is_err()
            );
        }
        assert!(
            NuclideOutflow::from_volumetric_flow(
                0,
                &[Volume::new(bad)],
                VolumetricFlowRate::new(0.0),
                1
            )
            .is_err()
        );
    }
    for flow in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            NuclideTransfer::from_volumetric_flow(0, 1, &valid, VolumetricFlowRate::new(flow), 2)
                .is_err()
        );
        assert!(
            NuclideOutflow::from_volumetric_flow(0, &valid, VolumetricFlowRate::new(flow), 2)
                .is_err()
        );
    }
    for nuclides in [0, MAX_DECAY_NODES + 1] {
        assert!(
            NuclideTransfer::from_volumetric_flow(
                0,
                1,
                &valid,
                VolumetricFlowRate::new(1.0),
                nuclides
            )
            .is_err()
        );
        assert!(
            NuclideOutflow::from_volumetric_flow(0, &valid, VolumetricFlowRate::new(1.0), nuclides)
                .is_err()
        );
    }
}
