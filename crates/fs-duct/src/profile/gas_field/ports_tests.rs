//! The orifice-terminal admittance on a continuous gas field (`response_with_wall_holes_and_ports`): a uniform field is the
//! homogeneous result bit for bit, on a graded field one port is the inverse of the one-port network impedance and the matrix
//! is reciprocal and passive, and under the open-hole interaction the port network is the COUPLED one (every other opening
//! linear, all radiating into one exterior with their mutual impedance): it is reciprocal and passive, it differs from the
//! independent network, and as the coupling vanishes it is the independent network. Gauntlet G0 (algebraic laws) on a graded
//! two-gas field.
use super::tests::field;
use super::*;

const RADIUS: f64 = 0.0073;
const OMEGAS_HZ: [f64; 5] = [150.0, 400.0, 900.0, 1800.0, 2600.0];

fn cylinder(length: f64) -> Segment {
    Segment::Cylinder { radius: RADIUS, length }
}
fn hole(radius: f64, chimney: f64) -> Segment {
    Segment::ToneHole { hole_radius: radius, chimney_height: chimney, bore_radius: RADIUS, state: crate::HoleState::Open }
}
/// Three open holes of different size on a 0.57 m bore. The field must cover the duct's length EXACTLY, so these
/// four lengths are chosen to sum to 0.57 bit for bit in the order the code adds them.
fn duct() -> Duct {
    Duct {
        segments: vec![
            cylinder(0.285),
            hole(0.0021, 0.0035),
            cylinder(0.0475),
            hole(0.0025, 0.0040),
            cylinder(0.0475),
            hole(0.0018, 0.0030),
            cylinder(0.19),
        ],
    }
}
fn graded_field() -> AxialGasField {
    field([285.0, 335.0], [0.8, 0.4], 101325.0)
}
fn graded() -> GasFieldDuct {
    GasFieldDuct::new(&duct(), &graded_field(), RiccatiSteps(32)).unwrap()
}
fn run(d: &GasFieldDuct, hz: f64, ports: usize) -> Result<(DuctResponse, Vec<OpenHoleFlow>, Vec<C64>), DuctError> {
    run_in(d, hz, ports, None)
}
fn run_in(d: &GasFieldDuct, hz: f64, ports: usize, exterior: Option<&GasState>)
    -> Result<(DuctResponse, Vec<OpenHoleFlow>, Vec<C64>), DuctError> {
    d.response_with_wall_holes_and_ports(
        core::f64::consts::TAU * hz, LossModel::Bessel, Termination::FlangedOpen, exterior, WallModel::Rigid, ports,
    )
}

#[test]
fn a_uniform_field_gives_the_homogeneous_ports_bit_for_bit() {
    let material = field([300.0, 300.0], [0.6, 0.6], 101325.0);
    let gas = *material.state_at(0.0).unwrap().state();
    let uniform = GasFieldDuct::new(&duct(), &material, RiccatiSteps(32)).unwrap();
    assert!(uniform.field.uniform_state().is_some());
    for hz in OMEGAS_HZ {
        let (_, _, y) = run(&uniform, hz, 3).unwrap();
        let (_, _, direct) = crate::exterior::input_impedance_holes_and_ports(
            &duct(), &gas, &gas, core::f64::consts::TAU * hz, LossModel::Bessel, Termination::FlangedOpen,
            WallModel::Rigid, 3,
        )
        .unwrap();
        assert_eq!(y.len(), 9);
        for (a, b) in y.iter().zip(&direct) {
            assert_eq!((a.re.to_bits(), a.im.to_bits()), (b.re.to_bits(), b.im.to_bits()), "{hz} Hz");
        }
    }
}

#[test]
fn on_a_graded_field_one_port_is_the_inverse_of_the_one_port_network_impedance() {
    let d = graded();
    for hz in OMEGAS_HZ {
        let (_, flows, y) = run(&d, hz, 1).unwrap();
        assert_eq!(y.len(), 1);
        let s = flows[0].shunt_port_reflection;
        let z0 = flows[0].shunt_port_reference_pa_s_m3;
        let z_net = (C64::ONE + s) / (C64::ONE - s);
        let expected = (C64::from_re(1.0) / z_net).scale(1.0 / z0);
        let err = (y[0] - expected).abs() / expected.abs();
        assert!(err < 1.0e-9, "{hz} Hz: {:?} vs {:?}, relative {err:e}", y[0], expected);
    }
}

#[test]
fn the_graded_three_port_admittance_is_reciprocal_passive_and_independent_of_the_port_count() {
    let d = graded();
    for hz in OMEGAS_HZ {
        let (_, flows, y3) = run(&d, hz, 3).unwrap();
        assert_eq!((flows.len(), y3.len()), (3, 9));
        let scale = y3.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for i in 0..3 {
            for j in 0..3 {
                assert!((y3[i * 3 + j] - y3[j * 3 + i]).abs() <= 1.0e-9 * scale, "{hz} Hz: Y[{i}][{j}] != Y[{j}][{i}]");
            }
        }
        for k in 0..64 {
            let x: Vec<C64> = (0..3)
                .map(|i| C64::new(((k * 7 + i * 3) as f64).sin(), ((k * 5 + i * 11) as f64).cos()))
                .collect();
            let mut q = C64::ZERO;
            for i in 0..3 {
                for j in 0..3 {
                    q = q + C64::new(x[i].re, -x[i].im) * y3[i * 3 + j] * x[j];
                }
            }
            assert!(q.re >= -1.0e-9 * scale, "{hz} Hz, probe {k}: Re(x^H Y x) = {:e}", q.re);
        }
        // The leading block does not depend on how many further holes are made ports.
        let (_, _, y1) = run(&d, hz, 1).unwrap();
        assert!((y1[0] - y3[0]).abs() <= 1.0e-9 * y1[0].abs(), "{hz} Hz: Y11 differs between 1 and 3 ports");
        // Zero ports is empty, more than the cap refuses.
        assert!(run(&d, hz, 0).unwrap().2.is_empty());
        assert!(run(&d, hz, crate::compose::MAX_ORIFICE_PORTS + 1).is_err());
    }
}

/// A card always declares its exterior; the coupling puts every opening into that one exterior, so every run in these tests
/// names it (with no exterior the plain path would use each hole's own gas, which is a different problem).
fn declared_exterior() -> GasState {
    let material = graded_field();
    *material.state_at(material.length_m()).unwrap().state()
}

fn coupled(epsilon: f64) -> GasFieldDuct {
    graded().with_open_hole_interaction(OpenHoleInteraction::new(epsilon).unwrap())
}

#[test]
fn the_interaction_changes_the_radiated_flows_and_the_port_network_with_them() {
    let exterior = declared_exterior();
    let plain = graded();
    let coupled = coupled(2.0);
    let (mut largest_flow_change, mut largest_port_change) = (0.0_f64, 0.0_f64);
    for hz in OMEGAS_HZ {
        let (_, f_plain, y_plain) = run_in(&plain, hz, 3, Some(&exterior)).unwrap();
        let (_, f_coupled, y_coupled) = run_in(&coupled, hz, 3, Some(&exterior)).unwrap();
        let scale = y_plain.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for (a, b) in y_plain.iter().zip(&y_coupled) {
            largest_port_change = largest_port_change.max((*a - *b).abs() / scale);
        }
        for (a, b) in f_plain.iter().zip(&f_coupled) {
            largest_flow_change = largest_flow_change
                .max((a.u_hole_over_u_in - b.u_hole_over_u_in).abs() / a.u_hole_over_u_in.abs());
        }
    }
    assert!(largest_flow_change > 1.0e-9,
        "the coupling must act on the radiated mouth flows, or this test proves nothing: largest relative change {largest_flow_change:e}");
    assert!(largest_port_change > 1.0e-6,
        "the port network must be the coupled one, not the independent one: largest relative change {largest_port_change:e}");
}

#[test]
fn a_vanishing_coupling_gives_the_independent_port_network() {
    // eps -> infinity switches the mutual impedance off: the coupled orifice-terminal network is then the independent
    // chain's, entry by entry (the two are built by different routes, so this ties them together).
    let exterior = declared_exterior();
    let plain = graded();
    let weak = coupled(1.0e15);
    for hz in OMEGAS_HZ {
        let (_, _, y_plain) = run_in(&plain, hz, 3, Some(&exterior)).unwrap();
        let (_, _, y_weak) = run_in(&weak, hz, 3, Some(&exterior)).unwrap();
        let scale = y_plain.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for (a, b) in y_plain.iter().zip(&y_weak) {
            assert!((*a - *b).abs() <= 1.0e-8 * scale, "{hz} Hz: {a:?} vs {b:?}");
        }
    }
}

#[test]
fn the_coupled_three_port_admittance_is_reciprocal_passive_and_independent_of_the_port_count() {
    let exterior = declared_exterior();
    let d = coupled(2.0);
    for hz in OMEGAS_HZ {
        let (_, flows, y3) = run_in(&d, hz, 3, Some(&exterior)).unwrap();
        assert_eq!((flows.len(), y3.len()), (3, 9));
        let scale = y3.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for i in 0..3 {
            for j in 0..3 {
                assert!((y3[i * 3 + j] - y3[j * 3 + i]).abs() <= 1.0e-8 * scale, "{hz} Hz: Y[{i}][{j}] != Y[{j}][{i}]");
            }
        }
        for k in 0..64 {
            let x: Vec<C64> = (0..3)
                .map(|i| C64::new(((k * 7 + i * 3) as f64).sin(), ((k * 5 + i * 11) as f64).cos()))
                .collect();
            let mut q = C64::ZERO;
            for i in 0..3 {
                for j in 0..3 {
                    q = q + C64::new(x[i].re, -x[i].im) * y3[i * 3 + j] * x[j];
                }
            }
            assert!(q.re >= -1.0e-8 * scale, "{hz} Hz, probe {k}: Re(x^H Y x) = {:e}", q.re);
        }
        // The leading block does not depend on how many further holes are made ports: the others stay linear openings.
        let (_, _, y1) = run_in(&d, hz, 1, Some(&exterior)).unwrap();
        assert!((y1[0] - y3[0]).abs() <= 1.0e-8 * y1[0].abs(), "{hz} Hz: Y11 differs between 1 and 3 ports");
        assert!(run_in(&d, hz, 0, Some(&exterior)).unwrap().2.is_empty());
        assert!(run_in(&d, hz, crate::compose::MAX_ORIFICE_PORTS + 1, Some(&exterior)).is_err());
    }
}
