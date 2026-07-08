//! Battery for the vortex particle method (fs-vpm). Each test is an analytic
//! vortex fixture: a single vortex induces a tangential Gamma/(2 pi r) field and
//! stays put, a counter-rotating pair self-propels at Gamma/(2 pi d), a
//! co-rotating pair conserves total circulation and its centroid.

use fs_vpm::{VortexParticle, induced_velocity, simulate, total_circulation, vorticity_centroid};
use std::f64::consts::{PI, TAU};

fn mag(v: [f64; 2]) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}

#[test]
fn a_single_vortex_induces_the_analytic_tangential_field() {
    // Gamma = 2 pi at the origin -> |u| = Gamma/(2 pi r) = 1/r, purely tangential.
    let v = [VortexParticle::new([0.0, 0.0], TAU)];
    let u1 = induced_velocity(&v, [1.0, 0.0], 0.0);
    assert!((mag(u1) - 1.0).abs() < 1e-12); // 1/r at r=1
    assert!(u1[0].abs() < 1e-12 && (u1[1] - 1.0).abs() < 1e-12); // tangential (+y)
    let u2 = induced_velocity(&v, [2.0, 0.0], 0.0);
    assert!((mag(u2) - 0.5).abs() < 1e-12); // 1/r at r=2
    // tangential everywhere: u . r = 0.
    let p = [1.3, -0.7];
    let u = induced_velocity(&v, p, 0.0);
    assert!((u[0] * p[0] + u[1] * p[1]).abs() < 1e-12);
}

#[test]
fn a_single_vortex_does_not_move_itself() {
    let start = [VortexParticle::new([0.3, -0.2], TAU)];
    let end = simulate(&start, 0.01, 100, 0.0);
    assert!((end[0].pos[0] - 0.3).abs() < 1e-12);
    assert!((end[0].pos[1] - (-0.2)).abs() < 1e-12);
}

#[test]
fn a_counter_rotating_pair_self_propels_at_the_analytic_speed() {
    // +Gamma at (0, 0.5), -Gamma at (0, -0.5): d = 1, speed Gamma/(2 pi d) = 1 in +x.
    let start = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    let t = 0.4;
    let end = simulate(&start, 0.01, 40, 0.0);
    // translated +x by speed*t = 0.4, vertical positions unchanged.
    assert!((end[0].pos[0] - t).abs() < 1e-6, "x = {}", end[0].pos[0]);
    assert!((end[1].pos[0] - t).abs() < 1e-6);
    assert!((end[0].pos[1] - 0.5).abs() < 1e-9 && (end[1].pos[1] + 0.5).abs() < 1e-9);
    // separation preserved (rigid translation).
    let sep = mag([end[0].pos[0] - end[1].pos[0], end[0].pos[1] - end[1].pos[1]]);
    assert!((sep - 1.0).abs() < 1e-9);
}

#[test]
fn a_co_rotating_pair_conserves_circulation_and_centroid() {
    // two +Gamma vortices orbit their shared centroid at the origin.
    let start = [
        VortexParticle::new([0.5, 0.0], TAU),
        VortexParticle::new([-0.5, 0.0], TAU),
    ];
    let c0 = total_circulation(&start);
    let end = simulate(&start, 0.005, 200, 1e-6);
    // total circulation is invariant.
    assert!((total_circulation(&end) - c0).abs() < 1e-12);
    assert!((c0 - 2.0 * TAU).abs() < 1e-12);
    // the centroid of vorticity stays at the origin.
    let centroid = vorticity_centroid(&end).unwrap();
    assert!(mag(centroid) < 1e-3, "centroid drifted to {centroid:?}");
    // they actually moved (orbiting), so this is not a trivial fixed point.
    assert!(mag([end[0].pos[0] - 0.5, end[0].pos[1]]) > 0.05);
}

#[test]
fn a_symmetric_pair_has_no_defined_centroid() {
    // total circulation zero -> centroid undefined (guarded).
    let pair = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    assert!(total_circulation(&pair).abs() < 1e-12);
    assert!(vorticity_centroid(&pair).is_none());
}

#[test]
fn the_desingularized_core_bounds_the_velocity() {
    // with a finite core, the self-field is finite even at the particle.
    let v = [VortexParticle::new([0.0, 0.0], TAU)];
    let u = induced_velocity(&v, [0.0, 0.0], 0.1);
    assert!(u[0].abs() < 1e-12 && u[1].abs() < 1e-12); // exactly at the particle
    // just off-core the speed is bounded by the desingularized kernel.
    let near = induced_velocity(&v, [0.01, 0.0], 0.1);
    assert!(mag(near) < 1.0 / (2.0 * PI) * 10.0); // far below the singular 1/r
}

#[test]
fn the_method_is_deterministic() {
    let start = [
        VortexParticle::new([0.5, 0.0], TAU),
        VortexParticle::new([-0.5, 0.0], TAU),
    ];
    let a = simulate(&start, 0.005, 50, 1e-6);
    let b = simulate(&start, 0.005, 50, 1e-6);
    assert_eq!(a[0].pos[0].to_bits(), b[0].pos[0].to_bits());
}
