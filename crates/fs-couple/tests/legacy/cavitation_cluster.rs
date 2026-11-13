#![cfg(feature = "cavitation")]

use fs_aeroac::curle3d::moving::{MovingCompactLimits, MovingCompactObserver};
use fs_aeroac::monopole::MonopoleLimits;
use fs_couple::cavitation::cluster::{MaterialBubbleCluster, MaterialBubbleClusterError};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_flux::cavitation::{
    BubbleCluster, BubbleClusterError, BubbleClusterLimits, BubbleClusterNumerics, BubbleGas,
    BubbleLimits, BubbleNumerics, BubbleState, ClusterBubble, RayleighPlessetBubble,
};
use fs_material::cavitation::CavitationState;
use fs_material::fluid_acoustics::FluidAcoustics;

const RADIUS: f64 = 50e-6;
const SEPARATION: f64 = 1e-3;
const AMBIENT: f64 = 4000.0;
const RECEIVER: [f64; 3] = [-0.03, 0.0, 0.0];

fn liquid() -> CavitationState {
    CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap()
}

fn native(source: CavitationState) -> BubbleCluster {
    let bubble = RayleighPlessetBubble::try_new(
        source,
        BubbleGas::equilibrium(source, RADIUS, AMBIENT, 1.4).unwrap(),
        BubbleLimits {
            minimum_radius_m: 0.5 * RADIUS,
            maximum_radius_m: 2.0 * RADIUS,
            maximum_wall_mach_number: 0.01,
            maximum_vapor_density_ratio: 0.01,
        },
    )
    .unwrap();
    BubbleCluster::try_new(
        vec![
            ClusterBubble {
                bubble,
                center_m: [0.0; 3],
            },
            ClusterBubble {
                bubble,
                center_m: [SEPARATION, 0.0, 0.0],
            },
        ],
        vec![
            BubbleState::try_new(0.0, 0.97 * RADIUS, 0.01).unwrap(),
            BubbleState::try_new(0.0, 1.02 * RADIUS, -0.02).unwrap(),
        ],
        0.0,
        AMBIENT,
        BubbleClusterLimits {
            maximum_bubbles: 16,
            maximum_radius_to_separation: 0.1,
            maximum_inertia_row_sum: 0.5,
            maximum_liquid_pressure_ratio: 0.01,
            maximum_source_angular_frequency_rad_s: core::f64::consts::TAU * 20_000.0,
            maximum_retardation_phase_rad: 0.1,
        },
        AMBIENT,
    )
    .unwrap()
}

fn numerics() -> BubbleClusterNumerics {
    BubbleClusterNumerics {
        radial: BubbleNumerics {
            relative_tolerance: 1e-11,
            radius_absolute_tolerance_m: 1e-18,
            velocity_absolute_tolerance_m_s: 1e-12,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-14,
            maximum_step_s: 1e-7,
            maximum_attempts: 20_000,
        },
        absolute_energy_tolerance_j: 1e-18,
    }
}

fn observer_limits() -> MonopoleLimits {
    MonopoleLimits {
        maximum_ka: 0.1,
        minimum_distance_radii: 10.0,
        maximum_pressure_ratio: 0.01,
    }
}

#[test]
fn sources_retain_the_joint_accelerations_centers_and_original_material() {
    let source = liquid();
    let owner = MaterialBubbleCluster::new(&source, native(source)).unwrap();
    assert!(std::ptr::eq(owner.source(), &source));
    assert_eq!(owner.material(), source);
    assert!((owner.source_bandwidth_hz() - 20_000.0).abs() < 1e-10);
    let before = owner.cluster().states().to_vec();
    let joint = owner.cluster().observe(AMBIENT).unwrap();
    let sources = owner.compact_source_samples(AMBIENT).unwrap();
    let motions = owner.source_motions().unwrap();
    assert_eq!(sources.len(), joint.bubbles.len());

    for (index, bubble) in joint.bubbles.iter().enumerate() {
        let member = owner.cluster().members()[index];
        let isolated = member.bubble.observe(&bubble.state, AMBIENT).unwrap();
        assert_ne!(
            bubble.wall_acceleration_m_s2,
            isolated.wall_acceleration_m_s2
        );
        assert_eq!(sources[index].source_time_s, joint.time_s);
        assert_eq!(
            sources[index].strength.volume_flow_m3_s,
            bubble.volume_rate_m3_s
        );
        assert_eq!(
            sources[index].strength.volume_flow_rate_m3_s2,
            bubble.volume_acceleration_m3_s2,
        );
        let radius = bubble.state.radius_m();
        let velocity = bubble.state.wall_velocity_m_s();
        let independent = 4.0
            * core::f64::consts::PI
            * radius
            * (radius * bubble.wall_acceleration_m_s2 + 2.0 * velocity * velocity);
        assert!(
            (sources[index].strength.volume_flow_rate_m3_s2 - independent).abs()
                < 1e-14 * independent.abs() + 1e-18
        );
        assert_eq!(sources[index].strength.force_on_fluid_n, [0.0; 3]);
        assert_eq!(sources[index].strength.force_rate_on_fluid_n_s, [0.0; 3]);
        assert_eq!(motions[index].reference_position_m(), member.center_m);
        assert_eq!(motions[index].velocity_m_s(), [0.0; 3]);
    }
    assert_eq!(owner.cluster().states(), before.as_slice());

    // Equal acoustic rho/c does not excuse a different viscous radial equation.
    let other = CavitationState::try_new(300.0, 1000.0, 0.002, 0.072, 2000.0, 1500.0).unwrap();
    assert!(matches!(
        MaterialBubbleCluster::new(&other, native(source)),
        Err(MaterialBubbleClusterError::MaterialMismatch),
    ));
}

#[test]
fn dropping_native_proposal_preserves_all_states_and_commit_replays_direct_advance() {
    let mut owner = native(liquid());
    let before = owner.states().to_vec();
    let expected = native(liquid())
        .advance(10e-6, &numerics(), &|_| AMBIENT)
        .unwrap();
    {
        let proposal = owner
            .propose_controlled(10e-6, &numerics(), &|_| AMBIENT, &mut || false)
            .unwrap();
        assert_eq!(proposal.step(), &expected);
    }
    assert_eq!(owner.time_s(), 0.0);
    assert_eq!(owner.states(), before.as_slice());
    let accepted = owner
        .propose_controlled(10e-6, &numerics(), &|_| AMBIENT, &mut || false)
        .unwrap()
        .commit();
    assert_eq!(accepted, expected);
    assert_eq!(owner.time_s(), 10e-6);
    for (state, bubble) in owner.states().iter().zip(&expected.observation.bubbles) {
        assert_eq!(*state, bubble.state);
    }
}

#[test]
fn observed_step_preserves_native_dynamics_and_separate_retarded_paths() {
    let source = liquid();
    let mut owner = MaterialBubbleCluster::new(&source, native(source)).unwrap();
    let observer = owner
        .monopole_observer(RECEIVER, observer_limits())
        .unwrap();
    assert_eq!(observer.receiver_position_m(), RECEIVER);
    let expected = native(source)
        .advance(20e-6, &numerics(), &|_| AMBIENT)
        .unwrap();
    let step = owner
        .advance_observed(20e-6, &numerics(), &|_| AMBIENT, &observer, &mut || false)
        .unwrap();
    assert_eq!(step.cluster, expected);
    assert_eq!(step.sources, owner.compact_source_samples(AMBIENT).unwrap());
    for (index, sample) in step.acoustic.iter().enumerate() {
        let distance = 0.03 + index as f64 * SEPARATION;
        let expected_pressure = source.liquid_density_kg_m3()
            * expected.observation.bubbles[index].volume_acceleration_m3_s2
            / (4.0 * core::f64::consts::PI * distance);
        assert!(
            (sample.pressure_pa - expected_pressure).abs()
                < 2e-14 * expected_pressure.abs() + 1e-16
        );
        assert_eq!(sample.source_time_s, 20e-6);
        assert!(
            (sample.arrival_time_s - (20e-6 + distance / source.liquid_sound_speed_m_s())).abs()
                < 2e-19
        );
    }
    assert_ne!(
        step.acoustic[0].arrival_time_s,
        step.acoustic[1].arrival_time_s
    );
}

#[test]
fn observer_refusal_and_final_cancellation_roll_back_all_members_then_retry_exactly() {
    let source = liquid();
    let mut owner = MaterialBubbleCluster::new(&source, native(source)).unwrap();
    let before = owner.cluster().states().to_vec();
    let strict = owner
        .monopole_observer(
            RECEIVER,
            MonopoleLimits {
                maximum_pressure_ratio: 1e-20,
                ..observer_limits()
            },
        )
        .unwrap();
    assert!(matches!(
        owner.advance_observed(20e-6, &numerics(), &|_| AMBIENT, &strict, &mut || false,),
        Err(MaterialBubbleClusterError::Acoustic(_)),
    ));
    assert_eq!(owner.cluster().time_s(), 0.0);
    assert_eq!(owner.cluster().states(), before.as_slice());

    let observer = owner
        .monopole_observer(RECEIVER, observer_limits())
        .unwrap();
    let mut reference = MaterialBubbleCluster::new(&source, native(source)).unwrap();
    let mut polls = 0;
    let expected = reference
        .advance_observed(20e-6, &numerics(), &|_| AMBIENT, &observer, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled_polls = 0;
    assert!(matches!(
        owner.advance_observed(20e-6, &numerics(), &|_| AMBIENT, &observer, &mut || {
            cancelled_polls += 1;
            cancelled_polls == polls
        },),
        Err(MaterialBubbleClusterError::Dynamics(
            BubbleClusterError::Cancelled
        )),
    ));
    assert_eq!(cancelled_polls, polls);
    assert_eq!(owner.cluster().time_s(), 0.0);
    assert_eq!(owner.cluster().states(), before.as_slice());
    let retried = owner
        .advance_observed(20e-6, &numerics(), &|_| AMBIENT, &observer, &mut || false)
        .unwrap();
    assert_eq!(retried, expected);
}

#[test]
fn geometry_and_medium_mismatch_refuse_before_evaluating_the_drive() {
    let source = liquid();
    let owner = MaterialBubbleCluster::new(&source, native(source)).unwrap();
    assert!(
        owner
            .monopole_observer([0.0; 3], observer_limits())
            .is_err()
    );
    assert!(
        owner
            .monopole_observer([0.0009, 0.0, 0.0], observer_limits())
            .is_err()
    );
    let observer = owner
        .monopole_observer(RECEIVER, observer_limits())
        .unwrap();

    let mut members = owner.cluster().members().to_vec();
    members.swap(0, 1);
    let changed = BubbleCluster::try_new(
        members,
        owner.cluster().states().to_vec(),
        0.0,
        AMBIENT,
        owner.cluster().limits(),
        AMBIENT,
    )
    .unwrap();
    let mut other = MaterialBubbleCluster::new(&source, changed).unwrap();
    assert!(matches!(
        other.advance_observed(
            10e-6,
            &numerics(),
            &|_| panic!("mismatched observer must refuse before drive evaluation"),
            &observer,
            &mut || false,
        ),
        Err(MaterialBubbleClusterError::ObserverMismatch),
    ));
    assert_eq!(other.cluster().time_s(), 0.0);

    let other_source =
        CavitationState::try_new(300.0, 990.0, 0.001, 0.072, 2000.0, 1500.0).unwrap();
    let mut other = MaterialBubbleCluster::new(&other_source, native(other_source)).unwrap();
    assert!(matches!(
        other.advance_observed(
            10e-6,
            &numerics(),
            &|_| panic!("mismatched medium must refuse before drive evaluation"),
            &observer,
            &mut || false,
        ),
        Err(MaterialBubbleClusterError::ObserverMismatch),
    ));
}

fn rendered_history(source_step_s: f64) -> Vec<f64> {
    let source = liquid();
    let mut owner = MaterialBubbleCluster::new(&source, native(source)).unwrap();
    let mut histories: Vec<_> = owner
        .compact_source_samples(AMBIENT)
        .unwrap()
        .into_iter()
        .map(|sample| vec![sample])
        .collect();
    let endpoints = (80e-6 / source_step_s).round() as usize;
    for index in 1..=endpoints {
        let step = owner
            .advance_compact_sources(
                index as f64 * source_step_s,
                &numerics(),
                &|_| AMBIENT,
                &mut || false,
            )
            .unwrap();
        for (history, sample) in histories.iter_mut().zip(step.sources) {
            history.push(sample);
        }
    }
    let states = owner.cluster().states().to_vec();
    let tracks: Vec<_> = histories
        .iter()
        .zip(owner.source_motions().unwrap())
        .map(|(samples, motion)| CompactSourceHistory {
            motion,
            samples,
            source_bandwidth_hz: owner.source_bandwidth_hz(),
        })
        .collect();
    let observer = MovingCompactObserver::new(
        source.acoustic_state().unwrap(),
        RECEIVER,
        MovingCompactLimits {
            maximum_mach: 0.0,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let mut renderer = MovingCompactRenderer::new(
        observer,
        &tracks,
        MovingPressureConfig {
            first_receiver_time_s: 0.031 / 1500.0 + 1e-6,
            sample_rate_hz: 500_000,
            samples: 32,
            max_block: 32,
            receiver_bandwidth_hz: 20_000.0,
            maximum_time_residual_s: 1e-15,
        },
        MovingPressureBudget {
            maximum_sources: 2,
            maximum_source_samples: 200,
            maximum_output_samples: 32,
            maximum_payload_bytes: 100_000,
            maximum_source_receiver_pairs: 64,
        },
        &mut || false,
    )
    .unwrap();
    let mut pressure = vec![0.0; 32];
    renderer.block(&mut pressure).unwrap();
    assert_eq!(owner.cluster().states(), states.as_slice());
    pressure
}

#[test]
fn coherent_microphone_history_converges_to_independently_advanced_retarded_sources() {
    let coarse = rendered_history(2e-6);
    let fine = rendered_history(1e-6);
    // Independently integrate each observation time; no history interpolation,
    // single shared delay, summed source acceleration or renderer in the oracle.
    let mut retarded_clusters = [native(liquid()), native(liquid())];
    let mut oracle = Vec::new();
    for index in 0..32 {
        let receiver_time = 0.031 / 1500.0 + 1e-6 + index as f64 / 500_000.0;
        let mut pressure = 0.0;
        for (member, cluster) in retarded_clusters.iter_mut().enumerate() {
            let distance = 0.03 + member as f64 * SEPARATION;
            let step = cluster
                .advance(receiver_time - distance / 1500.0, &numerics(), &|_| AMBIENT)
                .unwrap();
            pressure += 1000.0 * step.observation.bubbles[member].volume_acceleration_m3_s2
                / (4.0 * core::f64::consts::PI * distance);
        }
        oracle.push(pressure);
    }
    let scale = oracle.iter().copied().map(f64::abs).fold(0.0, f64::max);
    assert!(scale > 0.0);
    let error = |samples: &[f64]| {
        samples
            .iter()
            .zip(&oracle)
            .map(|(sample, expected)| (sample - expected).abs())
            .fold(0.0, f64::max)
    };
    let coarse_error = error(&coarse);
    let fine_error = error(&fine);
    assert!(
        fine_error < 0.4 * coarse_error,
        "{fine_error} vs {coarse_error}"
    );
    assert!(fine_error < 1e-3 * scale, "{fine_error} at scale {scale}");
}
