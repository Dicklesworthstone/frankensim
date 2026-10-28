//! Real boundary/SGS stepping with spatially resolved inlet traces.

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundaryProfileError3, BoundarySpec3, Face3, FaceBoundary3, InletFlux3,
};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_math::det;

fn grid(axis: usize, upper: bool) -> BoundaryGrid3 {
    let mut faces = [FaceBoundary3::Periodic; 6];
    let mut velocity = [0.0; 3];
    velocity[axis] = if upper { -0.04 } else { 0.04 };
    faces[2 * axis + usize::from(upper)] = FaceBoundary3::Velocity { velocity };
    faces[2 * axis + usize::from(!upper)] = FaceBoundary3::Pressure { density: 1.0 };
    BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::new(faces))
}

fn bits(grid: &BoundaryGrid3) -> Vec<u64> {
    let mut values = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                values.extend(grid.populations(x, y, z).map(f64::to_bits));
            }
        }
    }
    values
}

#[test]
fn resolved_profiles_drive_all_six_face_orientations_with_controlled_bulk_flow() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for upper in [false, true] {
            let face = Face3::ALL[2 * axis + usize::from(upper)];
            let sign = f64::from(face.normal()[axis]);
            let mut domain = grid(axis, upper);
            let mut calls = 0;
            let report = domain
                .set_velocity_profile(face, InletFlux3::MeanNormal(-0.04), &gate, |p| {
                    calls += 1;
                    let mut u = [0.0; 3];
                    u[axis] = -sign
                        * (0.047
                            + 0.003
                                * det::sin(
                                    core::f64::consts::TAU * p[(axis + 1) % 3] as f64 / 8.0,
                                ));
                    u[(axis + 1) % 3] = 0.002 * p[(axis + 2) % 3] as f64;
                    Ok(u)
                })
                .unwrap();
            assert_eq!(calls, 64);
            assert_eq!(report.sampled_cells, 64);
            assert!((sign * report.prescribed_mean_velocity_lu[axis] + 0.04).abs() < 1e-15);
            assert!((report.normal_velocity_correction_lu - 0.007).abs() < 1e-15);
            let profile = domain.velocity_profile(face).unwrap().to_vec();
            domain.step();
            let mut index = 0;
            for z in 0..8 {
                for y in 0..8 {
                    for x in 0..8 {
                        let p = [x, y, z];
                        if p[axis] != if upper { 7 } else { 0 } {
                            continue;
                        }
                        let actual = domain.velocity(x, y, z);
                        for a in 0..3 {
                            assert!(
                                (actual[a] - profile[index][a]).abs() < 2e-15,
                                "{face:?} {p:?}: {actual:?}"
                            );
                        }
                        index += 1;
                    }
                }
            }
            assert_eq!(index, 64);
        }
    }
}

#[test]
fn uniform_profile_and_clear_preserve_the_original_boundary_trajectory() {
    let gate = CancelGate::new_clock_free();
    let mut ordinary = grid(0, false);
    let mut profiled = grid(0, false);
    let report = profiled
        .set_velocity_profile(Face3::XMin, InletFlux3::MeanNormal(-0.04), &gate, |_| {
            Ok([0.04, 0.0, 0.0])
        })
        .unwrap();
    assert_eq!(report.normal_velocity_correction_lu, 0.0);
    for _ in 0..3 {
        ordinary.step();
        profiled.step();
        assert_eq!(bits(&ordinary), bits(&profiled));
    }
    profiled.clear_velocity_profile(Face3::XMin);
    assert!(profiled.velocity_profile(Face3::XMin).is_none());
    ordinary.step();
    profiled.step();
    assert_eq!(bits(&ordinary), bits(&profiled));
}

fn inflow(p: [usize; 3], tick: usize) -> Result<[f64; 3], BoundaryProfileError3> {
    let phase = core::f64::consts::TAU * (p[1] as f64 + p[2] as f64 - 0.04 * tick as f64) / 8.0;
    Ok([
        0.04 + 0.004 * det::sin(phase),
        -0.003 * det::sin(phase),
        -0.001 * det::sin(phase),
    ])
}

#[test]
fn evolving_profiles_feed_the_same_serial_and_pooled_sgs_step() {
    let pool = TilePool::for_host(2, 0x1A1E7);
    let gate = CancelGate::new_clock_free();
    let mut serial = grid(0, false);
    let mut pooled = grid(0, false);
    let mut serial_work = serial.sgs_workspace().unwrap();
    let mut pooled_work = pooled.sgs_workspace().unwrap();
    let config = SgsConfig {
        model: SgsModel::Wale { cw: 0.5 },
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    };
    let mut addresses = Vec::new();
    let mut active_sgs = false;
    for tick in 1..=5 {
        for domain in [&mut serial, &mut pooled] {
            domain
                .set_velocity_profile(Face3::XMin, InletFlux3::MeanNormal(-0.04), &gate, |p| {
                    inflow(p, tick)
                })
                .unwrap();
        }
        addresses.push(serial.velocity_profile(Face3::XMin).unwrap().as_ptr());
        let a = serial
            .step_with_sgs(&gate, config, &mut serial_work, None)
            .unwrap();
        let b = pooled
            .step_pooled_with_sgs(&pool, &gate, config, &mut pooled_work, None)
            .unwrap();
        active_sgs |= a.sgs.maximum_eddy_viscosity_lu > 0.0;
        assert_eq!(a, b);
        assert_eq!(bits(&serial), bits(&pooled));
    }
    assert!(
        active_sgs,
        "spatial boundary gradients must reach the closure"
    );
    assert_eq!(addresses[0], addresses[2]);
    assert_eq!(addresses[1], addresses[3]);
}

#[test]
fn refused_or_cancelled_update_keeps_the_entire_previous_trace() {
    let gate = CancelGate::new_clock_free();
    let mut domain = grid(0, false);
    domain
        .set_velocity_profile(Face3::XMin, InletFlux3::Unconstrained, &gate, |p| {
            inflow(p, 0)
        })
        .unwrap();
    let before = domain.velocity_profile(Face3::XMin).unwrap().to_vec();
    let populations = bits(&domain);
    let mut calls = 0;
    let error = domain
        .set_velocity_profile(Face3::XMin, InletFlux3::Unconstrained, &gate, |_| {
            calls += 1;
            if calls == 20 {
                Err(BoundaryProfileError3::Source("field refused".into()))
            } else {
                Ok([0.03, 0.0, 0.0])
            }
        })
        .unwrap_err();
    assert!(matches!(error, BoundaryProfileError3::Source(_)));
    assert_eq!(domain.velocity_profile(Face3::XMin).unwrap(), before);
    let error = domain
        .set_velocity_profile(Face3::XMin, InletFlux3::MeanNormal(-0.17), &gate, |_| {
            Ok([0.04, 0.05, 0.0])
        })
        .unwrap_err();
    assert!(matches!(
        error,
        BoundaryProfileError3::InvalidVelocity { .. }
    ));
    assert_eq!(domain.velocity_profile(Face3::XMin).unwrap(), before);
    let cancelled = CancelGate::new_clock_free();
    let error = domain
        .set_velocity_profile(Face3::XMin, InletFlux3::Unconstrained, &cancelled, |_| {
            cancelled.request();
            Ok([0.02, 0.0, 0.0])
        })
        .unwrap_err();
    assert_eq!(error, BoundaryProfileError3::Cancelled);
    assert_eq!(domain.velocity_profile(Face3::XMin).unwrap(), before);
    assert_eq!(bits(&domain), populations);
    assert_eq!(
        domain
            .set_velocity_profile(Face3::XMax, InletFlux3::Unconstrained, &gate, |_| panic!(
                "pressure face must refuse before sampling"
            ))
            .unwrap_err(),
        BoundaryProfileError3::NotVelocityFace
    );
    let mut walled = BoundaryGrid3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        BoundarySpec3::velocity_pressure_duct_z([0.0, 0.0, 0.04], 1.0),
    );
    assert_eq!(
        walled
            .set_velocity_profile(Face3::ZMin, InletFlux3::MeanNormal(-0.04), &gate, |_| Ok([
                0.0, 0.0, 0.04
            ]))
            .unwrap_err(),
        BoundaryProfileError3::MixedWallSeam
    );
}
