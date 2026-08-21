//! E4.7b battery (bead wf-root-guzez.5.19): the V-06a exactness
//! doctrine crossed onto the fs-vpm hybrid wake, plus the V-10-side
//! conversion symmetry. Per-item oracles throughout: on-plane normal
//! cancellation probe by probe; mirror∘coarsen == coarsen∘mirror
//! BITWISE per row per station; mirror∘aggregate == aggregate∘mirror
//! BITWISE per cell moment; the ledger exclusion proven by EXECUTING
//! the forbidden materialization and watching the per-station impulse
//! ledger corrupt; stable image identities; plane admission at the
//! exact clearance boundary (cap AND cap+1); determinism golden.
//! Repro: cargo test -p fs-vpm --test images3d_battery

use fs_vpm::coarsen3d::{coarsen_oldest, station_impulse};
use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode};
use fs_vpm::filament3d::{FilamentWake, WakeRateCertificate};
use fs_vpm::images3d::{
    ImagePlane, ImageSystem, grounded_hybrid_velocity, image_identity, unlawful_materialize_images,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-vpm-images3d\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;

/// Lifting line at z = 1 m over the plane z = 0; wake convects
/// downstream and slightly UP so every node stays strictly above.
fn line_at(z: f64) -> Vec<[f64; 3]> {
    (0..=8).map(|i| [0.0, i as f64 - 4.0, z]).collect()
}

fn flying_wake(rows: usize) -> FilamentWake {
    let cert = WakeRateCertificate {
        shed_hz: 120.0,
        n_stations: 8,
        max_rows: rows + 4,
    };
    let mut wake = FilamentWake::new(cert, line_at(1.0)).unwrap();
    let g: Vec<f64> = (0..8)
        .map(|s| {
            let y = (s as f64 - 3.5) / 4.0;
            5.0 * (1.0 - y * y).max(0.0)
        })
        .collect();
    for _ in 0..rows {
        wake.shed(&g, [-13.0 * DT, 0.0, 0.004]).unwrap();
    }
    wake
}

fn mag(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

const PLANE: ImagePlane = ImagePlane { z0_m: 0.0 };

#[test]
fn v06a_on_plane_normal_velocity_cancels() {
    let mut near = flying_wake(64);
    let far = FarField::aggregate(&mut near, 32, 8, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let images = ImageSystem::build(&near, &far, PLANE).unwrap();
    // Probe the plane itself, spanning under the line and the wake.
    let mut worst = 0.0f64;
    for i in 0..5 {
        for j in 0..3 {
            let p = [-6.0 + 3.0 * i as f64, -3.0 + 3.0 * j as f64, 0.0];
            let v = grounded_hybrid_velocity(&near, &far, &images, p);
            let rel = v[2].abs() / mag(v).max(1e-300);
            assert!(rel < 1e-9, "normal residual at ({},{}) = {rel}", p[0], p[1]);
            worst = worst.max(rel);
        }
    }
    jlog("v06a-cross", &format!("\"worst_normal_rel\":{worst}"));
}

#[test]
fn conversion_symmetry_mirror_commutes_with_coarsen() {
    // coarsen(mirror(w)) vs mirror(coarsen(w)) — with the plane at
    // z = 0 the reflection is EXACT in floats, so the commutation is
    // asserted BITWISE per row per station (never a totals-only sum).
    let base = flying_wake(16);
    let mut path_a = base.clone(); // coarsen then mirror (via images).
    coarsen_oldest(&mut path_a, 4).unwrap();
    let far_a = FarField::aggregate(&mut path_a, 4, 2, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let img_a = ImageSystem::build(&path_a, &far_a, PLANE).unwrap();

    // Mirror first: build a mirrored TWIN of the base wake by the
    // lawful constructor, then coarsen/aggregate the twin.
    let empty_far = FarField {
        mode: WakeCoreEvolutionMode::Frozen,
        cells: Vec::new(),
    };
    let img_base = ImageSystem::build(&base, &empty_far, PLANE).unwrap();
    let mut path_b = img_base.near_image_wake().clone();
    coarsen_oldest(&mut path_b, 4).unwrap();
    let far_b = FarField::aggregate(&mut path_b, 4, 2, WakeCoreEvolutionMode::Frozen, DT).unwrap();

    // Row-level commutation: mirrored-then-coarsened rows must equal
    // the coarsened-then-mirrored rows bit for bit.
    let a_img_rows = &img_a.near_image_wake().rows;
    assert_eq!(a_img_rows.len(), path_b.rows.len());
    for (r, (ra, rb)) in a_img_rows.iter().zip(path_b.rows.iter()).enumerate() {
        for s in 0..8 {
            assert_eq!(
                ra.gamma[s].to_bits(),
                rb.gamma[s].to_bits(),
                "gamma row {r} station {s}"
            );
        }
        for (n, (na, nb)) in ra.nodes.iter().zip(rb.nodes.iter()).enumerate() {
            for c in 0..3 {
                assert_eq!(na[c].to_bits(), nb[c].to_bits(), "node row {r}.{n}.{c}");
            }
        }
    }
    // Cell-level commutation: mirror(aggregate) vs aggregate(mirror),
    // every moment bitwise.
    let a_cells = img_a.far_image_cells();
    assert_eq!(a_cells.len(), far_b.cells.len());
    for (ci, (ca, cb)) in a_cells.iter().zip(far_b.cells.iter()).enumerate() {
        for c in 0..3 {
            assert_eq!(
                ca.centroid[c].to_bits(),
                cb.centroid[c].to_bits(),
                "cell {ci} centroid {c}"
            );
            assert_eq!(
                ca.w_sum[c].to_bits(),
                cb.w_sum[c].to_bits(),
                "cell {ci} w_sum {c}"
            );
        }
        assert_eq!(ca.w_abs_sum.to_bits(), cb.w_abs_sum.to_bits());
        assert_eq!(ca.radius_m.to_bits(), cb.radius_m.to_bits());
        for j in 0..3 {
            for l in 0..3 {
                assert_eq!(ca.dip[j][l].to_bits(), cb.dip[j][l].to_bits());
                for m in 0..3 {
                    assert_eq!(
                        ca.quad[j][l][m].to_bits(),
                        cb.quad[j][l][m].to_bits(),
                        "cell {ci} quad {j}{l}{m}"
                    );
                }
            }
        }
    }
    jlog(
        "conversion-symmetry",
        &format!(
            "\"rows_bitwise\":{},\"cells_bitwise\":{}",
            a_img_rows.len(),
            a_cells.len()
        ),
    );
}

#[test]
fn ledger_exclusion_proven_by_forbidden_materialization() {
    let mut near = flying_wake(32);
    let far = FarField::aggregate(&mut near, 8, 4, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let impulse_before: Vec<f64> = (0..8).map(|s| station_impulse(&near, s)).collect();
    let digest_before = near.digest();
    // LAWFUL path: build + evaluate images. The physical ledgers and
    // digest must be BIT-IDENTICAL afterward — images never enter.
    let images = ImageSystem::build(&near, &far, PLANE).unwrap();
    let _ = grounded_hybrid_velocity(&near, &far, &images, [0.3, 0.0, 0.9]);
    for s in 0..8 {
        assert_eq!(
            station_impulse(&near, s).to_bits(),
            impulse_before[s].to_bits(),
            "station {s} ledger touched by lawful images"
        );
    }
    assert_eq!(near.digest(), digest_before, "digest touched");
    // FORBIDDEN materialization (executed): the same images appended
    // as physical rows corrupt the per-station impulse visibly —
    // proving the exclusion is real AND the oracle is sensitive.
    let mut corrupted = near.clone();
    unlawful_materialize_images(&mut corrupted, PLANE);
    let mut worst = 0.0f64;
    for s in 0..8 {
        let after = station_impulse(&corrupted, s);
        let rel = ((after - impulse_before[s]) / impulse_before[s].abs().max(1e-12)).abs();
        worst = worst.max(rel);
    }
    assert!(
        worst > 0.2,
        "materialized images must visibly corrupt the impulse ledger: {worst}"
    );
    jlog(
        "ledger-exclusion",
        &format!("\"lawful_bitwise\":true,\"forbidden_worst_rel\":{worst}"),
    );
}

#[test]
fn identity_stability_and_plane_admission() {
    let mut near = flying_wake(16);
    let far = FarField::aggregate(&mut near, 4, 2, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    // Stability: two builds, same identity; a different plane and a
    // different source both move it.
    let a = ImageSystem::build(&near, &far, PLANE).unwrap();
    let b = ImageSystem::build(&near, &far, PLANE).unwrap();
    assert_eq!(a.identity, b.identity, "identity stable");
    let other = ImageSystem::build(&near, &far, ImagePlane { z0_m: -1.0 }).unwrap();
    assert_ne!(a.identity, other.identity, "plane enters identity");
    assert_ne!(
        a.identity,
        image_identity("someone-else", PLANE),
        "source enters identity"
    );
    // Admission at the EXACT clearance boundary. The lowest node of
    // this wake is the lifting line at z = 1.0: a plane one ulp below
    // admits; AT the node refuses; NaN refuses.
    let below = 1.0 - f64::EPSILON;
    assert!(ImagePlane { z0_m: below }.admit(&near).is_ok());
    for bad in [1.0, 1.5, f64::NAN] {
        assert_eq!(
            ImagePlane { z0_m: bad }.admit(&near).unwrap_err().code,
            "image-plane-invalid",
            "plane {bad}"
        );
    }
    // The lawful constructor routes through admission.
    assert_eq!(
        ImageSystem::build(&near, &far, ImagePlane { z0_m: 1.0 })
            .unwrap_err()
            .code,
        "image-plane-invalid"
    );
    jlog("identity", &format!("\"id\":\"{}\"", a.identity));
}

#[test]
fn ground_effect_direction_and_determinism_golden() {
    let mut near = flying_wake(64);
    let far = FarField::aggregate(&mut near, 32, 8, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let images = ImageSystem::build(&near, &far, PLANE).unwrap();
    // Physics: at the lifting line the image system OPPOSES the wake
    // downwash (ground effect reduces induced downwash) — direction
    // asserted, magnitude REPORTED.
    let probe = [0.3, 0.0, 1.05];
    let free = {
        let f = far.eval(probe);
        let n = near.induced_velocity(probe);
        [n[0] + f[0], n[1] + f[1], n[2] + f[2]]
    };
    let grounded = grounded_hybrid_velocity(&near, &far, &images, probe);
    assert!(
        grounded[2].abs() < free[2].abs(),
        "ground effect must reduce downwash: {} vs {}",
        grounded[2],
        free[2]
    );
    // Determinism golden over grounded velocities at fixed probes.
    let run = || {
        let mut b = Vec::new();
        for i in 0..4 {
            let p = [
                0.5 * i as f64 - 1.0,
                0.7 * i as f64 - 1.0,
                0.6 + 0.2 * i as f64,
            ];
            let v = grounded_hybrid_velocity(&near, &far, &images, p);
            for c in v {
                b.extend_from_slice(&c.to_bits().to_le_bytes());
            }
        }
        fs_blake3_hex(&b)
    };
    let a = run();
    assert_eq!(a, run(), "bit-identical twice");
    jlog(
        "golden",
        &format!(
            "\"digest\":\"{a}\",\"free_w\":{},\"grounded_w\":{}",
            free[2], grounded[2]
        ),
    );
    assert_eq!(
        a, "cebf414b1ba1b5086b71afb372ab0b3f8bebf39f056e066d34f334b3d827f503",
        "grounded-hybrid golden moved — determinism regression or an \
         intentional image/multipole change requiring the golden-bump protocol"
    );
}

fn fs_blake3_hex(b: &[u8]) -> String {
    fs_blake3::hash_domain("org.frankensim.fs-vpm.images3d.battery.v1", b).to_hex()
}
