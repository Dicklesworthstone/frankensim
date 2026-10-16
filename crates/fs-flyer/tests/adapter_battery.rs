//! E3.2-iii adapter battery (bead wf-root-guzez.4.2.3): bit-exact
//! round-trips through fs-geom and fs-mbd, the double-cover law (canonical
//! path REFUSES the negative representative; normalizing path resigns it),
//! rotation-semantics agreement between the spine's trajectory and fs-mbd's
//! rotate_body_to_world, and an integrated spine→adapter→spine round-trip.
//! Repro: cargo test -p fs-flyer --test adapter_battery

use fs_flyer::adapter::{
    geom_to_pos, mbd_to_quat, mbd_to_state, mbd_to_vec, pos_to_geom, quat_to_mbd_canonical,
    quat_to_mbd_normalizing, state_to_mbd, vec_to_mbd,
};
use fs_flyer::spine::{Loads, RigidBody, SixDofState, advance};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-adapter\",\"case\":\"{case}\",{payload}}}");
}

fn awkward() -> [f64; 3] {
    // Values with no short binary representation, incl. a negative zero.
    [0.1 + 0.2, -0.0, 1.0e-308]
}

#[test]
fn vector_and_point_round_trips_are_bit_exact() {
    let v = awkward();
    let via_geom = geom_to_pos(pos_to_geom(v));
    let via_mbd = mbd_to_vec(vec_to_mbd(v));
    for i in 0..3 {
        assert_eq!(
            v[i].to_bits(),
            via_geom[i].to_bits(),
            "geom round-trip bit {i}"
        );
        assert_eq!(
            v[i].to_bits(),
            via_mbd[i].to_bits(),
            "mbd round-trip bit {i}"
        );
    }
    // Negative zero survives (bit_exact means the SIGN of zero too).
    assert_eq!((-0.0f64).to_bits(), via_geom[1].to_bits());
    jlog(
        "round-trips",
        "\"bit_exact\":true,\"includes\":\"-0.0, subnormal-adjacent\"",
    );
}

#[test]
fn double_cover_law_two_spellings_no_middle_ground() {
    // A canonical unit quaternion (positive w) round-trips bit-exactly.
    let q = {
        // Build an exactly-normalized quaternion via fs-mbd itself.
        fs_mbd::UnitQuaternion::new(0.9, 0.1, -0.3, 0.2)
            .unwrap()
            .components()
    };
    let rt = mbd_to_quat(quat_to_mbd_canonical(q).unwrap());
    for i in 0..4 {
        assert_eq!(
            q[i].to_bits(),
            rt[i].to_bits(),
            "canonical round-trip bit {i}"
        );
    }
    // The NEGATIVE representative: canonical path REFUSES (typed, with the
    // normalizing repair ranked first)...
    let neg = [-q[0], -q[1], -q[2], -q[3]];
    let refusal = quat_to_mbd_canonical(neg).unwrap_err();
    assert_eq!(refusal.code, "quaternion-not-canonical");
    assert!(refusal.ranked_repairs[0].contains("normalizing"));
    // ...while the normalizing path resigns it to the SAME rotation.
    let resigned = quat_to_mbd_normalizing(neg).unwrap().components();
    for i in 0..4 {
        assert!(
            (resigned[i] - q[i]).abs() < 1e-15,
            "resign must recover q at {i}"
        );
    }
    // Degenerate quaternions refuse on both paths.
    assert_eq!(
        quat_to_mbd_normalizing([0.0; 4]).unwrap_err().code,
        "quaternion-invalid"
    );
    assert_eq!(
        quat_to_mbd_canonical([f64::NAN, 0.0, 0.0, 0.0])
            .unwrap_err()
            .code,
        "quaternion-not-canonical"
    );
    jlog(
        "double-cover",
        "\"canonical_refuses_negative\":true,\"normalizing_resigns\":true",
    );
}

#[test]
fn rotation_semantics_agree_between_spine_and_mbd() {
    // Spin the spine 90° about +z (yaw right), then rotate body-x into the
    // world through the ADAPTED quaternion: it must land on world +y
    // (frd/NED: nose right of north = east = +y).
    let body = RigidBody {
        mass_kg: 1.0,
        inertia_kgm2: [1.0, 1.0, 1.0],
    };
    let start = SixDofState {
        pos_m: [0.0; 3],
        vel_mps: [0.0; 3],
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: [0.0, 0.0, core::f64::consts::FRAC_PI_2], // 90°/s yaw
    };
    let dt = 1.0 / 120.0;
    let (end, _) = advance(&body, &start, 0.0, dt, 120, |_, _| Loads {
        force_n: [0.0; 3],
        moment_nm: [0.0; 3],
    })
    .unwrap();
    let q = quat_to_mbd_normalizing(end.quat).unwrap();
    let world_x = q.rotate_body_to_world(vec_to_mbd([1.0, 0.0, 0.0]));
    assert!(world_x.x.abs() < 1e-9, "x residue {}", world_x.x);
    assert!(
        (world_x.y - 1.0).abs() < 1e-9,
        "body-x must land on +y, got {}",
        world_x.y
    );
    assert!(world_x.z.abs() < 1e-12);
    // And the inverse direction agrees.
    let back = q.rotate_world_to_body(world_x);
    assert!((back.x - 1.0).abs() < 1e-12 && back.y.abs() < 1e-9);
    jlog(
        "rotation-semantics",
        &format!("\"world_x\":[{},{},{}]", world_x.x, world_x.y, world_x.z),
    );
}

#[test]
fn full_state_round_trip_through_the_seam() {
    // Run a coupled trajectory, cross the seam both ways, and require the
    // rebuilt state to be BIT-IDENTICAL (the replay-path guarantee).
    let body = RigidBody {
        mass_kg: 340.17,
        inertia_kgm2: [1787.0, 367.4, 1820.9],
    };
    let start = SixDofState {
        pos_m: [10.0, -3.0, -2.0],
        vel_mps: [13.9, 0.0, 0.0],
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: [0.01, -0.02, 0.03],
    };
    let (end, digests) = advance(&body, &start, 0.0, 1.0 / 120.0, 240, |t, s| Loads {
        force_n: [5.0 * (2.0 * t).sin(), 1.0, 3336.0 - 30.0 * s.vel_mps[2]],
        moment_nm: [2.0, 8.0 * t.cos(), -1.0],
    })
    .unwrap();
    let (p, v, q, w) = state_to_mbd(&end).expect("a live trajectory stays canonical here");
    let rebuilt = mbd_to_state(p, v, q, w);
    assert_eq!(
        fs_flyer::spine::tick_digest(999, &end),
        fs_flyer::spine::tick_digest(999, &rebuilt),
        "seam round-trip must preserve the state digest bit-exactly"
    );
    jlog(
        "state-round-trip",
        &format!(
            "\"ticks\":240,\"final_digest\":\"{}\"",
            digests.last().unwrap()
        ),
    );
}
