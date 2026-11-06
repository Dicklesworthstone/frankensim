//! The orifice-terminal admittance of the first open holes (`open_hole_port_admittance`): one port is the
//! existing one-port network impedance, the matrix is reciprocal and passive, and its size is bounded.
use crate::exterior::input_impedance_holes_and_ports;
use crate::{Duct, HoleState, LossModel, Segment, Termination, WallModel};
use fs_material::gas::GasState;
use fs_math::c64::C64;

fn gas() -> GasState {
    GasState::try_new_moist_air(293.15, 101_325.0, 0.5).unwrap()
}

fn hole(radius: f64, chimney: f64, state: HoleState) -> Segment {
    Segment::ToneHole { hole_radius: radius, chimney_height: chimney, bore_radius: 0.00745, state }
}

fn cylinder(length: f64) -> Segment {
    Segment::Cylinder { radius: 0.00745, length }
}

/// Three open holes of different size, one closed hole between them, and a bell load.
fn duct() -> Duct {
    Duct {
        segments: vec![
            cylinder(0.10),
            hole(0.00135, 0.0127, HoleState::Open),
            cylinder(0.03),
            hole(0.00235, 0.0069, HoleState::Closed),
            cylinder(0.02),
            hole(0.00235, 0.0069, HoleState::Open),
            cylinder(0.02),
            hole(0.0020, 0.0050, HoleState::Open),
            cylinder(0.35),
        ],
    }
}

fn run(duct: &Duct, f_hz: f64, ports: usize) -> (Vec<crate::OpenHoleFlow>, Vec<C64>) {
    let state = gas();
    let (_, flows, y) = input_impedance_holes_and_ports(
        duct, &state, &state, 2.0 * core::f64::consts::PI * f_hz, LossModel::Bessel,
        Termination::FlangedOpen, WallModel::Rigid, ports,
    )
    .unwrap();
    (flows, y)
}

#[test]
fn one_port_admittance_is_the_inverse_of_the_one_port_network_impedance() {
    // Z_net of `open_hole_port_reflection` is Z0 (1 + S) / (1 - S); the ladder must reproduce it for every
    // open hole taken as the only port when the others are linear.
    let duct = duct();
    for f in [150.0, 300.0, 700.0, 1500.0, 3000.0] {
        let (flows, y) = run(&duct, f, 1);
        assert_eq!(y.len(), 1);
        let s = flows[0].shunt_port_reflection;
        let z0 = flows[0].shunt_port_reference_pa_s_m3;
        let z_net = (C64::ONE + s) / (C64::ONE - s);
        let expected = (C64::from_re(1.0) / z_net).scale(1.0 / z0);
        let err = (y[0] - expected).abs() / expected.abs();
        assert!(err < 1.0e-9, "{f} Hz: ladder {:?} vs one-port {:?}, relative {err:e}", y[0], expected);
    }
}

#[test]
fn port_admittance_is_reciprocal_and_passive() {
    let duct = duct();
    for f in [150.0, 400.0, 900.0, 1800.0, 2600.0, 4000.0] {
        let (flows, y) = run(&duct, f, 3);
        assert_eq!(flows.len(), 3);
        assert_eq!(y.len(), 9);
        let scale = y.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for i in 0..3 {
            for j in 0..3 {
                let d = (y[i * 3 + j] - y[j * 3 + i]).abs();
                assert!(d <= 1.0e-9 * scale, "{f} Hz: Y[{i}][{j}] != Y[{j}][{i}] by {d:e}");
            }
        }
        // x^H Y x has a nonnegative real part for every x: the network only absorbs.
        for k in 0..64 {
            let x: Vec<C64> = (0..3)
                .map(|i| C64::new(((k * 7 + i * 3) as f64).sin(), ((k * 5 + i * 11) as f64).cos()))
                .collect();
            let mut q = C64::ZERO;
            for i in 0..3 {
                for j in 0..3 {
                    q = q + C64::new(x[i].re, -x[i].im) * y[i * 3 + j] * x[j];
                }
            }
            assert!(q.re >= -1.0e-9 * scale, "{f} Hz, probe {k}: Re(x^H Y x) = {:e}", q.re);
        }
    }
}

#[test]
fn the_port_count_is_bounded_and_limited_by_the_open_holes() {
    let duct = duct();
    let (_, y) = run(&duct, 700.0, 2);
    assert_eq!(y.len(), 4);
    let (_, y0) = run(&duct, 700.0, 0);
    assert!(y0.is_empty());
    // The duct has three open holes: a request for more than the cap refuses, never truncates silently.
    let state = gas();
    let refused = input_impedance_holes_and_ports(
        &duct, &state, &state, 2.0 * core::f64::consts::PI * 700.0, LossModel::Bessel,
        Termination::FlangedOpen, WallModel::Rigid, crate::compose::MAX_ORIFICE_PORTS + 1,
    );
    assert!(refused.is_err());
}

#[test]
fn the_leading_block_does_not_change_with_the_port_count() {
    // A port that is not asked for stays a linear open hole, which is exactly a port with zero drop.
    let duct = duct();
    for f in [200.0, 800.0, 2000.0] {
        let (_, y1) = run(&duct, f, 1);
        let (_, y2) = run(&duct, f, 2);
        let (_, y3) = run(&duct, f, 3);
        assert!((y1[0] - y2[0]).abs() <= 1.0e-9 * y1[0].abs(), "{f} Hz: Y11 differs between 1 and 2 ports");
        assert!((y1[0] - y3[0]).abs() <= 1.0e-9 * y1[0].abs(), "{f} Hz: Y11 differs between 1 and 3 ports");
        for (i, j) in [(0, 1), (1, 0), (1, 1)] {
            let a = y2[i * 2 + j];
            let b = y3[i * 3 + j];
            assert!((a - b).abs() <= 1.0e-9 * a.abs().max(1.0e-12), "{f} Hz: block ({i},{j}) differs between 2 and 3 ports");
        }
    }
}

fn pair_duct(before: f64, after: f64) -> Duct {
    Duct {
        segments: vec![
            cylinder(before),
            hole(0.00235, 0.0069, HoleState::Open),
            cylinder(0.0002),
            hole(0.00235, 0.0069, HoleState::Open),
            cylinder(after),
        ],
    }
}

#[test]
fn the_antisymmetric_drive_of_a_coincident_pair_never_reaches_the_bore() {
    // Equal and opposite drops (+e, -e) on two identical holes at one node push a flow around the two chimneys and
    // back through the node, so Y11 - Y12 is the chimney's own admittance and cannot depend on the bore that carries
    // the pair. The symmetric combination Y11 + Y12 does depend on it; a wrong sign of the off-diagonal term would
    // swap the two.
    for f in [150.0, 300.0, 600.0] {
        let (_, ya) = run(&pair_duct(0.12, 0.30), f, 2);
        let (_, yb) = run(&pair_duct(0.20, 0.55), f, 2);
        let local_a = ya[0] - ya[1];
        let local_b = yb[0] - yb[1];
        assert!((local_a - local_b).abs() <= 1.0e-2 * local_a.abs(), "{f} Hz: Y11-Y12 {local_a:?} vs {local_b:?}");
        let shared_a = ya[0] + ya[1];
        let shared_b = yb[0] + yb[1];
        assert!((shared_a - shared_b).abs() > 5.0 * (local_a - local_b).abs(), "{f} Hz: Y11+Y12 should feel the bore much more than Y11-Y12 does");
    }
}
