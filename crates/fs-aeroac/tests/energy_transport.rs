//! Energy follows the coherent physical field and original mean-fluid frame.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::CompactFieldStatistics;
use fs_aeroac::curle3d::energy::{instantaneous, time_averaged};
use fs_aeroac::curle3d::moving::flow::{UniformFlowCompactObserver, UniformFlowFrame};
use fs_aeroac::curle3d::moving::trajectory::{
    MovingMultipoleHistory, MovingMultipoleSource, MovingReceiverBudget, MovingReceiverConfig,
    MovingSourceArrayBudget, PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedMovingSourceArray, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{CompactSourceStrength, MovingCompactLimits, SourceKinematics};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.25, 20.0).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn upstream_downstream_and_crossflow_waves_transport_at_their_group_velocity() {
    let fluid = medium();
    let wind = [7.0, -2.0, 1.0];
    let pressure = 3.0;
    for n in [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]] {
        let velocity = n.map(|v| v * pressure / fluid.characteristic_impedance_pa_s_m());
        let result = instantaneous(fluid, wind, pressure, velocity).unwrap();
        let energy = pressure * pressure / fluid.adiabatic_bulk_modulus_pa();
        near(result.energy_density_j_m3, energy, 4e-16);
        for i in 0..3 {
            near(
                result.energy_flux_w_m2[i],
                (wind[i] + 20.0 * n[i]) * energy,
                5e-16,
            );
        }
    }
    // Transverse/vortical motion carries kinetic energy even at a pressure null.
    let shear = instantaneous(fluid, wind, 0.0, [0.0, 2.0, 0.0]).unwrap();
    assert_eq!(shear.pressure_work_flux_w_m2, [0.0; 3]);
    assert_eq!(shear.energy_density_j_m3, 2.5);
    assert_eq!(shear.energy_flux_w_m2, wind.map(|v| 2.5 * v));
    let rest = instantaneous(fluid, [0.0; 3], 3.0, [-0.0, 0.5, -0.7]).unwrap();
    assert_eq!(
        rest.energy_flux_w_m2.map(f64::to_bits),
        rest.pressure_work_flux_w_m2.map(f64::to_bits)
    );
}

#[test]
fn averaging_reuses_original_second_moments_without_another_half_factor() {
    let p = C64::new(3.0, -1.0);
    let u = [C64::new(0.2, 0.1), C64::new(-0.1, 0.3), C64::new(0.1, -0.2)];
    let work = u.map(|v| (p * v.conj()).scale(0.5));
    let mut moments = CompactFieldStatistics {
        mean_square_pressure_pa2: 0.5 * p.norm_sq(),
        velocity_covariance_m2_s2: core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    C64::from_re(0.5 * u[i].norm_sq())
                } else {
                    (u[i] * u[j].conj()).scale(0.5)
                }
            })
        }),
        active_intensity_w_m2: work.map(|v| v.re),
        reactive_intensity_w_m2: work.map(|v| v.im),
        scope: "manufactured coherent spectral line",
    };
    let wind = [4.0, -2.0, 1.0];
    let actual = time_averaged(medium(), wind, &moments).unwrap();
    let mut energy = 0.0;
    let mut flux = [0.0; 3];
    for phase in [
        C64::ONE,
        C64::new(0.0, 1.0),
        C64::from_re(-1.0),
        C64::new(0.0, -1.0),
    ] {
        let sample =
            instantaneous(medium(), wind, (p * phase).re, u.map(|v| (v * phase).re)).unwrap();
        energy += sample.energy_density_j_m3 * 0.25;
        for i in 0..3 {
            flux[i] += sample.energy_flux_w_m2[i] * 0.25;
        }
    }
    near(actual.energy_density_j_m3, energy, 1e-15);
    for i in 0..3 {
        near(actual.energy_flux_w_m2[i], flux[i], 1e-15);
    }
    moments.velocity_covariance_m2_s2[1][1] = C64::from_re(-1.0);
    assert!(time_averaged(medium(), wind, &moments).is_err());
}

fn monopole_energy(
    position: [f64; 3],
    time: f64,
) -> fs_aeroac::curle3d::energy::AcousticEnergyTransport {
    let wind = [4.0, -2.0, 1.0];
    let path = UniformFlowCompactObserver::new(
        medium(),
        wind,
        position,
        [0.0; 3],
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let prepared = path.prepare_receiver(time, 1e-12).unwrap();
    let t = prepared.source_time_s();
    let field = prepared
        .pressure_in_fluid_frame()
        .evaluate_monopole_field(0.03 * (11.0 * t).sin(), 0.33 * (11.0 * t).cos())
        .unwrap();
    instantaneous(
        medium(),
        wind,
        field.sample.pressure_pa,
        field.particle_velocity_m_s,
    )
    .unwrap()
}

fn derivative(h: f64, mut f: impl FnMut(f64) -> f64) -> f64 {
    (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
}

#[test]
fn native_convected_near_fields_close_source_free_energy_continuity() {
    let position = [1.2, -0.7, 0.9];
    let time = 0.37;
    let mut previous = f64::INFINITY;
    for h in [0.04, 0.02, 0.01] {
        let dt = derivative(h / 20.0, |t| {
            monopole_energy(position, time + t).energy_density_j_m3
        });
        let mut div = 0.0;
        let mut pressure_work_div = 0.0;
        for i in 0..3 {
            div += derivative(h, |d| {
                let mut x = position;
                x[i] += d;
                monopole_energy(x, time).energy_flux_w_m2[i]
            });
            pressure_work_div += derivative(h, |d| {
                let mut x = position;
                x[i] += d;
                monopole_energy(x, time).pressure_work_flux_w_m2[i]
            });
        }
        let error = (dt + div).abs() / dt.abs().max(div.abs());
        assert!(error < 0.12 * previous, "h={h}, continuity error={error:e}");
        assert!((dt + pressure_work_div).abs() > 0.01 * dt.abs());
        previous = error;
    }
    assert!(previous < 1e-8);
}

fn trajectory(
    position: [f64; 3],
    velocity: [f64; 3],
    frame: Option<UniformFlowFrame>,
) -> PreparedKinematicTrajectory {
    let knots = [-1.0, 1.0].map(|t| TimedSourceKinematics {
        source_time_s: t,
        kinematics: SourceKinematics {
            position_m: core::array::from_fn(|i| position[i] + velocity[i] * t),
            velocity_m_s: velocity,
            acceleration_m_s2: [0.0; 3],
        },
    });
    let budget = TrajectoryBudget {
        maximum_segments: 1,
        maximum_payload_bytes: 65536,
        maximum_bound_nodes: 16384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    };
    match frame {
        Some(f) => PreparedKinematicTrajectory::in_uniform_flow(
            medium(),
            &knots,
            f,
            0.8,
            budget,
            &mut || false,
        ),
        None => PreparedKinematicTrajectory::new(medium(), &knots, 0.8, budget, &mut || false),
    }
    .unwrap()
}

#[test]
fn moving_array_binds_original_mean_flow_preserves_coherence_and_cancellation() {
    for frame in [
        None,
        Some(UniformFlowFrame::new([3.0, -1.0, 0.5], 0.0, [0.0; 3]).unwrap()),
    ] {
        let source = trajectory([0.0; 3], [0.6, 0.2, -0.1], frame);
        let receiver = trajectory([3.0, 2.0, 1.0], [-2.0, 0.3, 0.1], frame);
        let pair = PreparedMovingReceiverSource::new(
            &source,
            &receiver,
            MovingReceiverConfig {
                receiver_time_interval_s: [0.0, 0.3],
                minimum_distance_m: 0.1,
                maximum_time_residual_s: 1e-12,
            },
            MovingReceiverBudget {
                maximum_segment_pairs: 64,
                maximum_bound_nodes: 16384,
                maximum_subdivision_depth: 20,
                maximum_root_iterations: 96,
            },
            &mut || false,
        )
        .unwrap();
        let sources = [pair, pair];
        let array = PreparedMovingSourceArray::new(
            &sources,
            MovingSourceArrayBudget {
                maximum_sources: 2,
                maximum_source_paths: 2,
                maximum_retarded_iterations: 192,
            },
            &mut || false,
        )
        .unwrap();
        let mut strength = |_, t: f64| {
            Ok(CompactSourceStrength {
                volume_flow_m3_s: 0.02 * (4.0 * t).cos(),
                volume_flow_rate_m3_s2: -0.08 * (4.0 * t).sin(),
                ..CompactSourceStrength::default()
            })
        };
        let field = array
            .force_field_at(
                0.2,
                &mut strength,
                &mut |_, _, _, _| Ok([0.0; 3]),
                &mut || false,
            )
            .unwrap();
        let result = array
            .force_energy_at(
                0.2,
                &mut strength,
                &mut |_, _, _, _| Ok([0.0; 3]),
                &mut || false,
            )
            .unwrap();
        assert_eq!(result.field, field);
        let expected = instantaneous(
            medium(),
            frame.map_or([0.0; 3], |f| f.velocity_m_s()),
            field.pressure_pa,
            field.particle_velocity_m_s,
        )
        .unwrap();
        assert_eq!(result.transport, expected);
        let mut polls = 0;
        let bright = array
            .multipole_energy_at(
                0.2,
                &mut |i, t| {
                    Ok(MovingMultipoleSource {
                        compact: strength(i, t)?,
                        ..MovingMultipoleSource::default()
                    })
                },
                &mut |_, _, _, _| Ok(MovingMultipoleHistory::default()),
                &mut || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        assert_eq!(bright.field, field);
        let cancelled = array.multipole_energy_at(
            0.2,
            &mut |i, t| {
                Ok(MovingMultipoleSource {
                    compact: strength(i, t)?,
                    ..MovingMultipoleSource::default()
                })
            },
            &mut |_, _, _, _| Ok(MovingMultipoleHistory::default()),
            &mut || {
                polls -= 1;
                polls == 0
            },
        );
        assert!(matches!(cancelled, Err(AeroacError::Cancelled { .. })));
        let dark = array
            .multipole_energy_at(
                0.2,
                &mut |i, t| {
                    let mut q = strength(i, t)?;
                    if i == 1 {
                        q.volume_flow_m3_s = -q.volume_flow_m3_s;
                        q.volume_flow_rate_m3_s2 = -q.volume_flow_rate_m3_s2;
                    }
                    Ok(MovingMultipoleSource {
                        compact: q,
                        ..MovingMultipoleSource::default()
                    })
                },
                &mut |_, _, _, _| Ok(MovingMultipoleHistory::default()),
                &mut || false,
            )
            .unwrap();
        assert_eq!(dark.transport.energy_density_j_m3, 0.0);
        assert_eq!(dark.transport.energy_flux_w_m2, [0.0; 3]);
    }
}

#[test]
fn scaled_energy_avoids_intermediate_overflow_and_refuses_nonfinite_transport() {
    let large_modulus = FluidAcousticState::try_new(1e100, 1e100).unwrap();
    near(
        instantaneous(large_modulus, [0.0; 3], 1e250, [0.0; 3])
            .unwrap()
            .energy_density_j_m3,
        5e199,
        1e-15,
    );
    let low_density = FluidAcousticState::try_new(1e-200, 1e100).unwrap();
    near(
        instantaneous(low_density, [0.0; 3], 0.0, [1e200, 0.0, 0.0])
            .unwrap()
            .energy_density_j_m3,
        5e199,
        1e-15,
    );
    assert!(instantaneous(medium(), [f64::NAN, 0.0, 0.0], 1.0, [0.0; 3]).is_err());
    assert!(instantaneous(medium(), [0.0; 3], 1.0, [f64::INFINITY, 0.0, 0.0]).is_err());
    assert!(instantaneous(large_modulus, [1e200, 0.0, 0.0], 1e250, [0.0; 3]).is_err());
    // Reduction adds no Mach restriction beyond the original propagation owner.
    assert!(instantaneous(medium(), [1000.0, 0.0, 0.0], 1.0, [0.0; 3]).is_ok());
}
