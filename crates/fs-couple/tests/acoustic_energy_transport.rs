//! G1/G3 energy transport uses coherent physical fields and the original fluid.

use fs_aeroac::curle3d::energy::AcousticEnergyTransport;
use fs_aeroac::curle3d::moving::flow::{UniformFlowFrame, UniformFlowReceiver};
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource, ReflectedTrajectoryBudget, TimedSourceKinematics,
    TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    TimedCompactSource, UniformSourceMotion,
};
use fs_aeroac::curle3d::quadrupole::recorded::QuadrupoleSampleClock;
use fs_aeroac::curle3d::quadrupole::transient::history::{
    PreparedQuadrupoleHistory, QuadrupoleHistoryBudget, QuadrupoleHistoryEnds,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::observation::moving::field::{CompactFieldSample, MovingMonopoleRenderer};
use fs_couple::pcm_wav::observation::moving::force::{ForceHistoryQuadrature, MovingForceRenderer};
use fs_couple::pcm_wav::observation::moving::{
    CompactMotion, CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget,
    MovingPressureConfig,
};
use fs_couple::pcm_wav::observation::quadrupole::QuadrupolePressureBudget;
use fs_couple::pcm_wav::observation::quadrupole::moving::MovingQuadrupoleSourceWindow;
use fs_couple::pcm_wav::observation::quadrupole::moving::field::{
    MovingQuadrupoleFieldBudget, MovingQuadrupoleFieldRenderer,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(3.2, 210.0).unwrap()
}
fn limits() -> MovingCompactLimits {
    MovingCompactLimits {
        maximum_mach: 0.8,
        minimum_distance_m: 0.01,
    }
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.1,
        sample_rate_hz: 1000,
        samples: 3,
        max_block: 3,
        receiver_bandwidth_hz: 20.0,
        maximum_time_residual_s: 1e-12,
    }
}
fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 2,
        maximum_source_samples: 8,
        maximum_output_samples: 3,
        maximum_payload_bytes: 131_072,
        maximum_source_receiver_pairs: 12,
    }
}
fn quadrature() -> ForceHistoryQuadrature {
    ForceHistoryQuadrature {
        absolute_tolerance_n_s_m3: 1e-13,
        maximum_subdivision_depth: 20,
        maximum_evaluations: 100_000,
    }
}
fn history(scale: f64) -> [TimedCompactSource; 2] {
    [0.0, 0.3].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: CompactSourceStrength {
            volume_flow_m3_s: scale * (1.0 + source_time_s) * 1e-4,
            volume_flow_rate_m3_s2: scale * 1e-4,
            ..Default::default()
        },
    })
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs() + 1e-24,
        "actual={actual:e}, expected={expected:e}"
    );
}
fn check(report: &AcousticEnergyTransport, field: CompactFieldSample, flow: [f64; 3]) {
    let rho = medium().density_kg_m3();
    let c = medium().sound_speed_m_s();
    let energy = 0.5
        * (field.pressure_pa.powi(2) / (rho * c * c)
            + rho
                * field
                    .particle_velocity_m_s
                    .iter()
                    .map(|v| v * v)
                    .sum::<f64>());
    close(report.energy_density_j_m3, energy);
    for axis in 0..3 {
        let work = field.pressure_pa * field.particle_velocity_m_s[axis];
        close(report.pressure_work_flux_w_m2[axis], work);
        close(report.energy_flux_w_m2[axis], work + flow[axis] * energy);
    }
}

#[test]
fn fixed_lab_flow_uses_original_material_and_physical_fields_without_advancing_audio() {
    let flow = [40.0, 0.0, 0.0];
    let receiver = UniformFlowReceiver::new(medium(), flow, [3.0, 0.0, 0.0], limits()).unwrap();
    let samples = history(1.0);
    let sources = [CompactSourceHistory {
        motion: receiver.bind_source([0.0; 3]).unwrap(),
        samples: &samples,
        source_bandwidth_hz: 1.0,
    }];
    let original =
        MovingCompactRenderer::with_uniform_flow(&sources, config(), budget(), &mut || false)
            .unwrap();
    assert_eq!(original.mean_flow_velocity_m_s(), flow);
    let mic = FirstOrder::new(0.0, [0.0, 1.0, 0.0]).unwrap();
    let mut renderer =
        MovingMonopoleRenderer::new(receiver, &sources, mic, config(), budget(), &mut || false)
            .unwrap();
    let mut fields = [CompactFieldSample::default(); 3];
    renderer.field_block(&mut fields).unwrap();
    assert_eq!(renderer.mean_flow_velocity_m_s(), flow);
    for (index, field) in fields.into_iter().enumerate() {
        assert_eq!(
            field.receiver_time_s,
            renderer.receiver_time_at(index as u64).unwrap()
        );
        // The transverse figure-eight pickup is null despite nonzero energy.
        assert_eq!(field.microphone_pressure_pa, 0.0);
        assert!(field.pressure_pa.abs() > 0.0);
        let report = renderer.energy_transport(&field).unwrap();
        check(&report, field, flow);
        assert!(report.energy_density_j_m3 > 0.0);
        assert!(report.energy_flux_w_m2[0] > report.pressure_work_flux_w_m2[0]);
        assert_eq!(renderer.samples_rendered(), 3);
    }
    assert_eq!(original.samples_rendered(), 0);
    assert_eq!(renderer.config(), config());
}

#[test]
fn moving_sources_interfere_before_energy_without_becoming_a_mean_flow() {
    let positive = history(1.0);
    let negative = history(-1.0);
    let motion = UniformSourceMotion::new(0.0, [0.0; 3], [12.0, -3.0, 0.0]).unwrap();
    let a = CompactSourceHistory {
        motion,
        samples: &positive,
        source_bandwidth_hz: 1.0,
    };
    let b = CompactSourceHistory {
        samples: &negative,
        ..a
    };
    let observer = MovingCompactObserver::new(medium(), [3.0, 1.0, 0.5], limits()).unwrap();
    let mut energies = [0.0; 3];
    for (case, sources) in [&[a][..], &[a, a][..], &[a, b][..]].into_iter().enumerate() {
        let mut renderer = MovingMonopoleRenderer::new(
            observer,
            sources,
            FirstOrder::default(),
            config(),
            budget(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(renderer.mean_flow_velocity_m_s(), [0.0; 3]);
        let mut field = [CompactFieldSample::default()];
        renderer.field_block(&mut field).unwrap();
        let report = renderer.energy_transport(&field[0]).unwrap();
        check(&report, field[0], [0.0; 3]);
        assert_eq!(report.energy_flux_w_m2, report.pressure_work_flux_w_m2);
        assert_eq!(renderer.samples_rendered(), 1);
        energies[case] = report.energy_density_j_m3;
    }
    assert!(energies[0] > 0.0);
    close(energies[1], 4.0 * energies[0]);
    assert_eq!(energies[2], 0.0);
}

fn path(
    position: [f64; 3],
    velocity: [f64; 3],
    frame: UniformFlowFrame,
) -> PreparedKinematicTrajectory {
    let knots = [0.0, 0.3].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: SourceKinematics {
            position_m: core::array::from_fn(|i| position[i] + velocity[i] * source_time_s),
            velocity_m_s: velocity,
            acceleration_m_s2: [0.0; 3],
        },
    });
    PreparedKinematicTrajectory::in_uniform_flow(
        medium(),
        &knots,
        frame,
        0.8,
        TrajectoryBudget {
            maximum_segments: 2,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 8192,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

#[test]
fn moving_receivers_and_stress_images_retain_fluid_flow_not_microphone_translation() {
    let position = [1.0, 2.0, 1.0];
    let reflection_budget = ReflectedTrajectoryBudget {
        maximum_payload_bytes: 131_072,
        maximum_bound_nodes: 8192,
        maximum_subdivision_depth: 20,
    };
    let pair_budget = MovingReceiverBudget {
        maximum_segment_pairs: 16,
        maximum_bound_nodes: 8192,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    };
    let pair_config = MovingReceiverConfig {
        receiver_time_interval_s: [0.1, 0.11],
        minimum_distance_m: 0.01,
        maximum_time_residual_s: 1e-12,
    };
    for flow in [[0.0; 3], [30.0, 7.0, 0.0]] {
        let frame = UniformFlowFrame::new(flow, 0.0, [0.0; 3]).unwrap();
        let source_path = path(position, [-7.0, 1.0, 0.0], frame);
        let receiver = path([4.0, 0.0, 2.0], [11.0, -3.0, 0.0], frame);
        let pair = PreparedMovingReceiverSource::new(
            &source_path,
            &receiver,
            pair_config,
            pair_budget,
            &mut || false,
        )
        .unwrap();
        assert_eq!((&pair).mean_flow_velocity_m_s(), flow);
        let samples = [0.0_f64, 0.3].map(|source_time_s| TimedCompactSource {
            source_time_s,
            strength: CompactSourceStrength {
                force_on_fluid_n: [source_time_s.powi(2), 0.0, 0.0],
                force_rate_on_fluid_n_s: [2.0 * source_time_s, 0.0, 0.0],
                ..Default::default()
            },
        });
        let tracks = [CompactSourceHistory {
            motion: &pair,
            samples: &samples,
            source_bandwidth_hz: 1.0,
        }];
        let mut force = MovingForceRenderer::new_with_zero_force_prehistory(
            &receiver,
            &tracks,
            FirstOrder::default(),
            config(),
            budget(),
            quadrature(),
            &mut || false,
        )
        .unwrap();
        let mut field = [CompactFieldSample::default()];
        force.field_block(&mut field).unwrap();
        assert_eq!(force.mean_flow_velocity_m_s(), flow);
        let report = force.energy_transport(&field[0]).unwrap();
        check(&report, field[0], flow);
        assert!(report.energy_density_j_m3 > 0.0);
        assert_eq!(force.samples_rendered(), 1);

        let plane =
            IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], IdealPlaneBoundary::Rigid).unwrap();
        let walls = MovingReflectionObserver::in_uniform_flow(
            &receiver,
            &[plane],
            2,
            reflection_budget,
            &mut || false,
        )
        .unwrap();
        let reflected = PreparedReflectedMovingReceiverSource::new(
            &source_path,
            &walls,
            pair_config,
            reflection_budget,
            pair_budget,
            &mut || false,
        )
        .unwrap();
        assert_eq!((&reflected).mean_flow_velocity_m_s(), flow);
        let moments = [[0.0; 6], [0.09, -0.018, 0.0, 0.027, 0.0, 0.0]];
        let history = PreparedQuadrupoleHistory::new(
            &[position],
            QuadrupoleSampleClock::new(0.0, 0.3, 0.0).unwrap(),
            &[0.0, 0.3],
            &moments,
            &[QuadrupoleHistoryEnds {
                initial_rate_n_m_s: [0.0; 6],
                final_rate_n_m_s: [0.6, -0.12, 0.0, 0.18, 0.0, 0.0],
                initial_integral_n_m_s: [0.0; 6],
            }],
            QuadrupoleHistoryBudget {
                maximum_sources: 1,
                maximum_samples_per_source: 2,
                maximum_workspace_bytes: 16_384,
            },
            || false,
        )
        .unwrap();
        let sources = [MovingQuadrupoleSourceWindow {
            motion: &reflected,
            source_time_interval_s: [0.0, 0.3],
            source_bandwidth_hz: 1.0,
        }];
        let mut stress = MovingQuadrupoleFieldRenderer::new_with_zero_stress_prehistory(
            &walls,
            &sources,
            &history,
            FirstOrder::default(),
            config(),
            MovingQuadrupoleFieldBudget {
                pressure: QuadrupolePressureBudget {
                    maximum_sources: 1,
                    maximum_output_samples: 3,
                    maximum_source_path_samples: 6,
                    maximum_scratch_bytes: 4096,
                },
                maximum_history_bytes: 16_384,
                quadrature: quadrature(),
            },
            &mut || false,
        )
        .unwrap();
        stress.field_block(&mut field).unwrap();
        assert_eq!(stress.mean_flow_velocity_m_s(), flow);
        let report = stress.energy_transport(&field[0]).unwrap();
        check(&report, field[0], flow);
        assert!(report.energy_density_j_m3 > 0.0);
        assert_eq!(stress.samples_rendered(), 1);
    }
}
