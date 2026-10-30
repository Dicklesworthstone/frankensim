//! Actual shared receiver and caller-owned accepted port flows, not a new solver.
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::pcm_wav::baffled::gas::{GasPropagation, GasPropagationLimits};
use fs_couple::pcm_wav::baffled::openings::{
    CircularOpening,
    duct::{DuctOpening, DuctPressure, DuctRadiationSpec, DuctReconstruction, DuctView},
};
use fs_couple::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_material::gas::{GasSpec, GasState};

fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101_325.0).unwrap()
}
fn graph() -> TubeNetworkSpec {
    TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.8 }],
        sections: vec![TubeSection { nodes: [0,1], length_m: 0.2, radius_m: 0.006, max_length_error_m: 0.001 }],
        sound_speed_m_s: gas(293.15).sound_speed, max_wave_memory_bytes: 1 << 20 }
}
fn spec() -> DuctRadiationSpec {
    DuctRadiationSpec {
        openings: vec![DuctOpening { network: 0, node: 1, center_m: [-0.025, 0.0] },
            DuctOpening { network: 1, node: 1, center_m: [0.025, 0.0] }],
        exterior: Some(gas(310.0)),
        receiver: CircularOutletReceiver { position_m: [0.1, 0.0, 0.4],
            radial_rings: 4, angular_points: 16, maximum_frequency_hz: 1000.0 },
        reconstruction: DuctReconstruction::Differentiating { mechanics_rate_hz: 192_000, output_rate_hz: 48_000 },
    }
}
fn receiver(s: DuctRadiationSpec) -> Result<DuctPressure, String> {
    let g = graph();
    DuctPressure::new(&[DuctView { graph: &g, gases: &[gas(293.15)] },
        DuctView { graph: &g, gases: &[gas(320.0)] }], s)
}
fn flows(n: usize) -> [f64; 2] {
    [1e-5 * (n as f64 * 0.037).sin(), -2e-5 * (n as f64 * 0.023).cos()]
}

#[test]
fn shared_packets_equal_complete_direct_traces_and_preserve_pending_copies() {
    let s = spec();
    let mut actual = receiver(s.clone()).unwrap();
    let openings: Vec<_> = s.openings.iter().map(|p| CircularOpening {
        center_m: p.center_m, radius_m: 0.006,
    }).collect();
    let exterior = s.exterior.unwrap();
    let mut reference = BaffledPressure::circular_openings_multirate_differentiating(
        &openings, 192_000, 48_000, s.receiver,
        RayleighMedium { density: exterior.density, sound_speed: exterior.sound_speed }).unwrap();
    let mut copy: Option<DuctPressure> = None;
    let mut trace = [0.0; 8];
    for n in 0..2400 {
        // The invalid second port may alter scratch but not a pending frame or FIR.
        let before = (actual.mechanical_samples(), actual.pending_frames(), actual.samples_rendered());
        assert!(actual.push_frame(|network, _| if network == 0 { Some(123.0) } else { None }).is_err());
        assert!(actual.push_frame(|_, _| Some(f64::NAN)).is_err());
        assert_eq!(before, (actual.mechanical_samples(), actual.pending_frames(), actual.samples_rendered()));
        let q = flows(n);
        trace[2*(n%4)..2*(n%4)+2].copy_from_slice(&q);
        let p = actual.push_frame(|network, node| { assert_eq!(node,1); Some(q[network]) }).unwrap();
        if let Some(c) = &mut copy {
            assert_eq!(p, c.push_frame(|network, _| Some(q[network])).unwrap());
        } else if n == 201 { copy = Some(actual.clone()); assert_eq!(actual.pending_frames(),2); }
        if n%4 == 3 { assert_eq!(p, Some(reference.step_trace(&trace).unwrap())); }
        else { assert_eq!(p,None); }
    }
    assert_eq!(actual.source(), &s);
    assert_eq!(*actual.exterior_gas(), exterior);
    assert_eq!(actual.mechanical_samples(),2400);
    assert_eq!(actual.samples_rendered(),600);
    assert_eq!(actual.pending_frames(),0);
}

#[test]
fn legacy_local_gas_and_single_opening_remain_exact() {
    let mut s = spec(); s.openings.truncate(1); s.openings[0].center_m = [0.0;2];
    s.exterior = None;
    s.reconstruction = DuctReconstruction::FiniteDifference { sample_rate_hz: 48_000 };
    let local = gas(293.15);
    let mut a = receiver(s.clone()).unwrap();
    let mut b = BaffledPressure::circular_outlet(0.006,48_000,s.receiver,
        RayleighMedium { density: local.density, sound_speed: local.sound_speed }).unwrap();
    assert_eq!(*a.exterior_gas(),local);
    for n in 0..600 {
        let q = flows(n)[0];
        assert_eq!(a.push_frame(|_,_|Some(q)).unwrap(), Some(b.step(&[q]).unwrap()));
    }
}

#[test]
fn incomplete_gas_address_and_environment_data_never_fall_back() {
    let mut s = spec(); s.exterior = None; assert!(receiver(s).is_err());
    let mut s = spec(); s.openings[1].network = 0; assert!(receiver(s).is_err());
    let mut s = spec(); s.openings[1].network = 2; assert!(receiver(s).is_err());
    let mut s = spec(); s.openings[1].node = 0; assert!(receiver(s).is_err());
    let mut s = spec(); s.openings[1].center_m = s.openings[0].center_m; assert!(receiver(s).is_err());
    let mut s = spec(); s.exterior.as_mut().unwrap().density *= 2.0; assert!(receiver(s).is_err());
    let mut s = spec(); s.exterior = Some(GasState::try_new(&GasSpec::dry_air_ussa1976(),310.0,90_000.0).unwrap());
    assert!(receiver(s).is_err());
    let mut s = spec(); s.reconstruction = DuctReconstruction::Differentiating {
        mechanics_rate_hz: 190_000, output_rate_hz: 48_000 };
    assert!(receiver(s).is_err());
    let mut s = spec(); s.openings.truncate(1);
    assert!(DuctPressure::new(&[DuctView { graph: &graph(), gases: &[] }],s).is_err());
}

#[test]
fn material_paths_retain_original_ports_exterior_and_multirate_packets() {
    let s = spec();
    let law = GasPropagation::classical(GasPropagationLimits::surface_reference());
    let mut actual = receiver(s.clone())
        .unwrap()
        .with_observer_gas(&law, &mut || false)
        .unwrap();
    let exterior = s.exterior.unwrap();
    let report = actual.observer_gas().unwrap();
    assert_eq!(report.gas, exterior);
    assert_eq!(report.surface_filter_indices.len(), 2 * 4 * 16);
    assert!(
        report.filters.len() > 1,
        "different surface distances must remain distinct"
    );
    let shortest = report
        .filters
        .iter()
        .map(|r| r.range_m)
        .fold(f64::INFINITY, f64::min);
    let longest = report.filters.iter().map(|r| r.range_m).fold(0.0, f64::max);
    assert!(longest - shortest > 0.005);
    let openings: Vec<_> = s
        .openings
        .iter()
        .map(|p| CircularOpening {
            center_m: p.center_m,
            radius_m: 0.006,
        })
        .collect();
    let (mut reference, _) = BaffledPressure::circular_openings_multirate_gas_differentiating(
        &openings,
        192_000,
        48_000,
        s.receiver,
        &exterior,
        &law,
        &mut || false,
    )
    .unwrap();
    let mut trace = [0.0; 8];
    let mut copy = None::<DuctPressure>;
    for n in 0..2400 {
        let q = flows(n);
        trace[2 * (n % 4)..2 * (n % 4) + 2].copy_from_slice(&q);
        let pressure = actual
            .push_frame(|network, node| {
                assert_eq!(node, 1);
                Some(q[network])
            })
            .unwrap();
        if let Some(c) = &mut copy {
            assert_eq!(
                pressure,
                c.push_frame(|network, _| Some(q[network])).unwrap()
            );
        } else if n == 201 {
            copy = Some(actual.clone());
        }
        if n % 4 == 3 {
            assert_eq!(pressure, Some(reference.step_trace(&trace).unwrap()));
        } else {
            assert_eq!(pressure, None);
        }
    }
    assert_eq!(actual.source(), &s);
    assert_eq!(actual.mechanical_samples(), 2400);
    assert!(actual.with_observer_gas(&law, &mut || false).is_err());
}

#[test]
fn material_observation_refuses_late_repeated_cancelled_or_implicit_reconstruction() {
    let law = GasPropagation::classical(GasPropagationLimits::surface_reference());
    let mut warm = receiver(spec()).unwrap();
    warm.push_frame(|_, _| Some(0.0)).unwrap();
    assert!(warm.with_observer_gas(&law, &mut || false).is_err());
    let cold = receiver(spec()).unwrap();
    assert!(cold.clone().with_observer_gas(&law, &mut || true).is_err());
    assert_eq!(cold.mechanical_samples(), 0);
    let selected = cold.with_observer_gas(&law, &mut || false).unwrap();
    assert!(selected.with_observer_gas(&law, &mut || false).is_err());
    let mut s = spec();
    s.reconstruction = DuctReconstruction::FiniteDifference {
        sample_rate_hz: 48_000,
    };
    assert!(
        receiver(s)
            .unwrap()
            .with_observer_gas(&law, &mut || false)
            .is_err()
    );
    let mut bounded = law;
    bounded.limits.max_unique_paths = 1;
    assert!(
        receiver(spec())
            .unwrap()
            .with_observer_gas(&bounded, &mut || false)
            .is_err()
    );
}
