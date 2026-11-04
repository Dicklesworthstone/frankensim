//! Held pressure traces on the original forced regularized boundary.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryForcingError3, BoundaryGrid3, BoundaryProfileError3, BoundarySpec3,
    Face3, FaceBoundary3, Q3,
};

fn grid(axis: usize, upper: bool, walled: bool) -> BoundaryGrid3 {
    let lateral = if walled { FaceBoundary3::stationary_wall() } else { FaceBoundary3::Periodic };
    let mut faces = [lateral; 6];
    faces[2 * axis + usize::from(upper)] = FaceBoundary3::Pressure { density: 1.002 };
    faces[2 * axis + usize::from(!upper)] = FaceBoundary3::Pressure { density: 1.0 };
    BoundaryGrid3::new(
        8, 8, 8, 0.83, [1e-6, -2e-6, 3e-6], BoundarySpec3::new(faces),
    )
}

fn bits(grid: &BoundaryGrid3) -> Vec<[u64; Q3]> {
    (0..512)
        .filter(|&i| !grid.is_solid(i % 8, (i / 8) % 8, i / 64))
        .map(|i| grid.populations(i % 8, (i / 8) % 8, i / 64).map(f64::to_bits))
        .collect()
}

fn density(cell: [usize; 3], axis: usize) -> f64 {
    1.002 + 0.0002 * cell[(axis + 1) % 3] as f64
        + 0.0001 * cell[(axis + 2) % 3] as f64
}

#[test]
fn uniform_pressure_profile_and_clear_keep_all_six_native_trajectories() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for upper in [false, true] {
            for walled in [false, true] {
                let face = Face3::ALL[2 * axis + usize::from(upper)];
                let mut ordinary = grid(axis, upper, walled);
                let mut profiled = grid(axis, upper, walled);
                let before = bits(&profiled);
                let report = profiled.set_pressure_profile(face, &gate, |_| Ok(1.002)).unwrap();
                assert_eq!(report.sampled_cells, 64);
                assert_eq!(report.minimum_density_lu, 1.002);
                assert_eq!(report.maximum_density_lu, 1.002);
                assert_eq!(bits(&profiled), before);
                for _ in 0..3 {
                    ordinary.step();
                    profiled.step();
                    assert_eq!(bits(&ordinary), bits(&profiled));
                }
                profiled.clear_pressure_profile(face);
                assert!(profiled.pressure_profile(face).is_none());
                ordinary.step();
                profiled.step();
                assert_eq!(bits(&ordinary), bits(&profiled));
            }
        }
    }
}

#[test]
fn spatial_pressure_profiles_prescribe_density_on_all_six_faces() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for upper in [false, true] {
            let face = Face3::ALL[2 * axis + usize::from(upper)];
            let mut domain = grid(axis, upper, false);
            let mut sampled = Vec::new();
            let report = domain.set_pressure_profile(face, &gate, |cell| {
                sampled.push(cell);
                Ok(density(cell, axis))
            }).unwrap();
            assert_eq!(sampled.len(), 64);
            assert_eq!(report.minimum_density_lu, 1.002);
            assert_eq!(report.maximum_density_lu, density([7; 3], axis));
            let held = domain.pressure_profile(face).unwrap().to_vec();
            domain.step();
            for (cell, expected) in sampled.into_iter().zip(held) {
                assert_eq!(cell[axis], if upper { 7 } else { 0 });
                assert_eq!(expected, density(cell, axis));
                let actual = domain.density(cell[0], cell[1], cell[2]);
                assert!((actual - expected).abs() < 3e-15, "{face:?} {cell:?}: {actual} != {expected}");
            }
        }
    }
}

#[test]
fn pressure_profile_refusal_reuses_buffers_and_preserves_live_target_precedence() {
    let gate = CancelGate::new_clock_free();
    let mut domain = grid(2, false, false);
    let face = Face3::ZMin;
    let before = bits(&domain);
    let mut addresses = Vec::new();
    for offset in [0.0, 0.001, 0.002] {
        domain.set_pressure_profile(face, &gate, |cell| Ok(density(cell, 2) + offset)).unwrap();
        addresses.push(domain.pressure_profile(face).unwrap().as_ptr());
    }
    assert_eq!(addresses[0], addresses[2]);
    assert_ne!(addresses[0], addresses[1]);
    let held = domain.pressure_profile(face).unwrap().to_vec();
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            domain.set_pressure_profile(face, &gate, |cell| {
                Ok(if cell[0] == 7 && cell[1] == 7 { bad } else { 1.001 })
            }),
            Err(BoundaryProfileError3::InvalidDensity { cell: [7, 7, 0] }),
        );
        assert_eq!(domain.pressure_profile(face).unwrap(), held);
    }
    let reason = BoundaryProfileError3::Source("reservoir field unavailable".into());
    assert_eq!(
        domain.set_pressure_profile(face, &gate, |_| Err(reason.clone())),
        Err(reason),
    );
    let cancelled = CancelGate::new_clock_free();
    assert_eq!(
        domain.set_pressure_profile(face, &cancelled, |_| {
            cancelled.request();
            Ok(1.001)
        }),
        Err(BoundaryProfileError3::Cancelled),
    );
    assert_eq!(domain.pressure_profile(face).unwrap(), held);
    let old = domain.boundaries();
    let mut faces = Face3::ALL.map(|face| old.face(face));
    faces[4] = FaceBoundary3::Pressure { density: 1.001 };
    assert!(matches!(
        domain.rebind_open_boundaries(BoundarySpec3::new(faces), &gate),
        Err(BoundaryForcingError3::Invalid(_)),
    ));
    assert_eq!(domain.boundaries(), old);
    faces[4] = old.face(face);
    faces[5] = FaceBoundary3::Pressure { density: 1.001 };
    domain.rebind_open_boundaries(BoundarySpec3::new(faces), &gate).unwrap();
    assert_eq!(domain.pressure_profile(face).unwrap(), held);
    assert_eq!(bits(&domain), before);
    // A failed scratch proposal must not leak into the next complete profile.
    domain.set_pressure_profile(face, &gate, |cell| Ok(density(cell, 2))).unwrap();
    assert_eq!(domain.pressure_profile(face).unwrap().as_ptr(), addresses[1]);
    let mut fresh = grid(2, false, false);
    fresh.rebind_open_boundaries(BoundarySpec3::new(faces), &gate).unwrap();
    fresh.set_pressure_profile(face, &gate, |cell| Ok(density(cell, 2))).unwrap();
    domain.step();
    fresh.step();
    assert_eq!(bits(&domain), bits(&fresh));
    let accepted = bits(&domain);
    domain.clear_pressure_profile(face);
    faces[4] = FaceBoundary3::Pressure { density: 1.001 };
    domain.rebind_open_boundaries(BoundarySpec3::new(faces), &gate).unwrap();
    assert_eq!(bits(&domain), accepted);
}

#[test]
fn pressure_profiles_retain_original_face_and_solid_admission() {
    let gate = CancelGate::new_clock_free();
    let mut wrong = BoundaryGrid3::new(
        8, 8, 8, 0.83, [0.0; 3],
        BoundarySpec3::velocity_pressure_duct_z([0.0, 0.0, 0.01], 1.0),
    );
    assert_eq!(
        wrong.set_pressure_profile(Face3::ZMin, &gate, |_| panic!("wrong face sampled")),
        Err(BoundaryProfileError3::NotPressureFace),
    );
    let mut blocked = grid(2, false, false);
    blocked.voxelize_sdf(|p| if p == [3.5, 3.5, 0.5] { -1.0 } else { 1.0 });
    let before = bits(&blocked);
    let report = blocked.set_pressure_profile(Face3::ZMin, &gate, |cell| {
        assert_ne!(cell, [3, 3, 0], "solid nodes must never be sampled");
        Ok(1.001)
    }).unwrap();
    assert_eq!(report.sampled_cells, 63);
    assert_eq!(blocked.pressure_profile(Face3::ZMin).unwrap().len(), 64);
    assert_eq!(bits(&blocked), before);
    let mut empty = grid(2, false, false);
    empty.voxelize_sdf(|p| p[2] - 1.0);
    assert_eq!(
        empty.set_pressure_profile(Face3::ZMin, &gate, |_| panic!("solid face sampled")),
        Err(BoundaryProfileError3::EmptyFace),
    );
    assert!(empty.pressure_profile(Face3::ZMin).is_none());
}
