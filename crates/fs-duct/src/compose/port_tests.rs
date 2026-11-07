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

/// Eight open holes of graded size and spacing on the same 7.45 mm bore (a fuller fingering than the three-hole duct).
fn eight_hole_duct() -> Duct {
    let mut segments = vec![cylinder(0.08)];
    for k in 0..8 {
        segments.push(hole(0.0013 + 0.0002 * k as f64, 0.004 + 0.001 * k as f64, HoleState::Open));
        segments.push(cylinder(0.015 + 0.002 * k as f64));
    }
    segments.push(cylinder(0.2));
    Duct { segments }
}

#[test]
fn eight_ports_are_reciprocal_passive_and_leave_the_leading_block_alone() {
    // The ladder is general in the port count (the cap is 16): the same laws as for three ports hold for eight.
    assert!(crate::compose::MAX_ORIFICE_PORTS >= 8);
    let duct = eight_hole_duct();
    for f in [200.0, 700.0, 1500.0, 3000.0] {
        let (flows, y8) = run(&duct, f, 8);
        assert_eq!((flows.len(), y8.len()), (8, 64));
        let scale = y8.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for i in 0..8 {
            for j in 0..8 {
                let d = (y8[i * 8 + j] - y8[j * 8 + i]).abs();
                assert!(d <= 1.0e-9 * scale, "{f} Hz: Y[{i}][{j}] != Y[{j}][{i}] by {d:e}");
            }
        }
        for k in 0..64 {
            let x: Vec<C64> = (0..8)
                .map(|i| C64::new(((k * 7 + i * 3) as f64).sin(), ((k * 5 + i * 11) as f64).cos()))
                .collect();
            let mut q = C64::ZERO;
            for i in 0..8 {
                for j in 0..8 {
                    q = q + C64::new(x[i].re, -x[i].im) * y8[i * 8 + j] * x[j];
                }
            }
            assert!(q.re >= -1.0e-9 * scale, "{f} Hz, probe {k}: Re(x^H Y x) = {:e}", q.re);
        }
        let (_, y4) = run(&duct, f, 4);
        for i in 0..4 {
            for j in 0..4 {
                let (a, b) = (y4[i * 4 + j], y8[i * 8 + j]);
                assert!((a - b).abs() <= 1.0e-9 * scale, "{f} Hz: block ({i},{j}) differs between 4 and 8 ports");
            }
        }
    }
}

/// Fourteen open holes of graded size and spacing on the same 7.45 mm bore: the clarinet's G4 fingering has 14.
fn fourteen_hole_duct() -> Duct {
    let mut segments = vec![cylinder(0.06)];
    for k in 0..14 {
        segments.push(hole(0.0012 + 0.00012 * k as f64, 0.004 + 0.0005 * k as f64, HoleState::Open));
        segments.push(cylinder(0.012 + 0.0015 * k as f64));
    }
    segments.push(cylinder(0.15));
    Duct { segments }
}

#[test]
fn fourteen_ports_are_reciprocal_passive_and_leave_the_leading_block_alone() {
    // Beyond the old cap of 12: the same laws hold for fourteen, and the leading 8 x 8 block is the eight-port network.
    assert!(crate::compose::MAX_ORIFICE_PORTS >= 14);
    let duct = fourteen_hole_duct();
    for f in [200.0, 700.0, 1500.0, 3000.0] {
        let (flows, y14) = run(&duct, f, 14);
        assert_eq!((flows.len(), y14.len()), (14, 196));
        let scale = y14.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for i in 0..14 {
            for j in 0..14 {
                let d = (y14[i * 14 + j] - y14[j * 14 + i]).abs();
                assert!(d <= 1.0e-9 * scale, "{f} Hz: Y[{i}][{j}] != Y[{j}][{i}] by {d:e}");
            }
        }
        for k in 0..64 {
            let x: Vec<C64> = (0..14)
                .map(|i| C64::new(((k * 7 + i * 3) as f64).sin(), ((k * 5 + i * 11) as f64).cos()))
                .collect();
            let mut q = C64::ZERO;
            for i in 0..14 {
                for j in 0..14 {
                    q = q + C64::new(x[i].re, -x[i].im) * y14[i * 14 + j] * x[j];
                }
            }
            assert!(q.re >= -1.0e-9 * scale, "{f} Hz, probe {k}: Re(x^H Y x) = {:e}", q.re);
        }
        let (_, y8) = run(&duct, f, 8);
        for i in 0..8 {
            for j in 0..8 {
                let (a, b) = (y8[i * 8 + j], y14[i * 14 + j]);
                assert!((a - b).abs() <= 1.0e-9 * scale, "{f} Hz: block ({i},{j}) differs between 8 and 14 ports");
            }
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

/// The dense implementation `open_hole_port_admittance` had before its nodal matrix was solved as the tridiagonal
/// matrix it is, frozen verbatim (n x n storage, `lu_complex`, one `solve` per port) as the oracle for the
/// bit-identity tests below.
mod frozen {
    use crate::DuctError;
    use crate::compose::{
        ComposeOptions, LadderBlock, MAX_ORIFICE_PORTS, PreparedChain, PreparedOpenHole, bad, bound, finite,
        pure_shunt_admittance,
    };
    use fs_la::eigen_complex::lu_complex;
    use fs_math::c64::C64;

    pub(super) fn dense_open_hole_port_admittance(
        chain: &PreparedChain,
        open_holes: &[PreparedOpenHole],
        ports: usize,
        z_load: Option<C64>,
        options: ComposeOptions,
    ) -> Result<Vec<C64>, DuctError> {
        if ports > MAX_ORIFICE_PORTS {
            return Err(bad("at most sixteen open holes can carry a coupled orifice"));
        }
        let m = ports.min(open_holes.len());
        if m == 0 {
            return Ok(Vec::new());
        }
        let holes = &open_holes[..m];
        let reference = chain.reference_impedance;
        let n = chain.blocks.len();
        let mut port_of = vec![None; n];
        for (j, hole) in holes.iter().enumerate() {
            if hole.block_index >= n || port_of[hole.block_index].is_some() {
                return Err(bad("coupled open holes need distinct blocks inside the chain"));
            }
            port_of[hole.block_index] = Some(j);
        }
        let mut ladder: Vec<LadderBlock> = Vec::with_capacity(n);
        let mut size = 1usize;
        for (i, block) in chain.blocks.iter().enumerate() {
            match port_of[i] {
                Some(j) => {
                    if holes[j].half_series_impedance != C64::ZERO {
                        size += 2;
                    }
                    ladder.push(LadderBlock::Port(j));
                }
                None => match pure_shunt_admittance(block) {
                    Some(g) => ladder.push(LadderBlock::Shunt(g)),
                    None => {
                        let (ym, condition) = block.admittance(options.maximum_condition)?;
                        bound(condition, options.maximum_condition)?;
                        size += 1;
                        ladder.push(LadderBlock::Series(ym));
                    }
                },
            }
        }
        let mut y = vec![C64::ZERO; size * size];
        let stamp = |y: &mut [C64], a: usize, b: usize, m2: [C64; 4]| {
            y[a * size + a] = y[a * size + a] + m2[0];
            y[a * size + b] = y[a * size + b] + m2[1];
            y[b * size + a] = y[b * size + a] + m2[2];
            y[b * size + b] = y[b * size + b] + m2[3];
        };
        let (mut left, mut next) = (0usize, 1usize);
        let mut internal = vec![0usize; m];
        for item in &ladder {
            match item {
                LadderBlock::Series(ym) => {
                    stamp(&mut y, left, next, *ym);
                    left = next;
                    next += 1;
                }
                LadderBlock::Shunt(g) => {
                    y[left * size + left] = y[left * size + left] + *g;
                }
                LadderBlock::Port(j) => {
                    let z = holes[*j].half_series_impedance;
                    if z == C64::ZERO {
                        internal[*j] = left;
                    } else {
                        let g = z.recip();
                        let series = [g, -g, -g, g];
                        stamp(&mut y, left, next, series);
                        stamp(&mut y, next, next + 1, series);
                        internal[*j] = next;
                        left = next + 1;
                        next += 2;
                    }
                }
            }
        }
        debug_assert_eq!(next, size);
        y[0] = y[0] + C64::from_re(1.0 / reference);
        match z_load {
            None => {}
            Some(z) if finite(z) && z != C64::ZERO => y[left * size + left] = y[left * size + left] + z.recip(),
            _ => return Err(DuctError::Singular),
        }
        let mut c_hat = Vec::with_capacity(m);
        for (j, hole) in holes.iter().enumerate() {
            let open = C64::from_re(hole.sigma) / hole.open_chimney_impedance;
            let [closed_p, closed_u] = hole.closed_return_coordinates;
            if closed_p == C64::ZERO && closed_u != C64::ZERO {
                return Err(DuctError::Singular);
            }
            let closed = if closed_u == C64::ZERO {
                C64::ZERO
            } else {
                closed_u / closed_p.scale(hole.shunt_port_reference_pa_s_m3)
            };
            if !finite(open) || !finite(closed) {
                return Err(DuctError::Singular);
            }
            let node = internal[j];
            y[node * size + node] = y[node * size + node] + open + closed;
            c_hat.push(open);
        }
        let lu = lu_complex(&y, size).map_err(|_| DuctError::Singular)?;
        let mut port = vec![C64::ZERO; m * m];
        for j in 0..m {
            let mut v = vec![C64::ZERO; size];
            v[internal[j]] = v[internal[j]] + c_hat[j];
            lu.solve(&mut v);
            for i in 0..m {
                let delta = if i == j { C64::ONE } else { C64::ZERO };
                let entry = -(c_hat[i] * (v[internal[i]] - delta));
                if !finite(entry) {
                    return Err(DuctError::Singular);
                }
                port[i * m + j] = entry;
            }
        }
        Ok(port)
    }
}

/// Bits with the sign of a zero removed: the one thing the tridiagonal graph may differ in.
fn bits(v: &[C64]) -> Vec<(u64, u64)> {
    v.iter().map(|c| ((c.re + 0.0).to_bits(), (c.im + 0.0).to_bits())).collect()
}

/// The admittance of one duct at one frequency from the tridiagonal solve and from the frozen dense oracle, on the
/// same prepared chain (both must agree on refusing as well).
fn both(duct: &Duct, f_hz: f64, ports: usize) -> (Result<Vec<C64>, crate::DuctError>, Result<Vec<C64>, crate::DuctError>) {
    use crate::compose::{ComposeOptions, duct_chain_in_exterior, exterior_options, open_hole_port_admittance};
    let state = gas();
    let omega = 2.0 * core::f64::consts::PI * f_hz;
    let options: ComposeOptions = exterior_options(duct, &state, &state);
    let (chain, open_holes) = duct_chain_in_exterior(
        duct, &state, &state, omega, LossModel::Bessel, WallModel::Rigid, options, true, || false,
    )
    .unwrap();
    let radius = duct.segments.last().unwrap().outlet_radius();
    let (z_load, _) = crate::termination_impedance(Termination::FlangedOpen, &state, radius, omega).unwrap();
    (
        open_hole_port_admittance(&chain, &open_holes, ports, z_load, options),
        frozen::dense_open_hole_port_admittance(&chain, &open_holes, ports, z_load, options),
    )
}

fn assert_same(duct: &Duct, f_hz: f64, ports: usize) {
    let (new, old) = both(duct, f_hz, ports);
    match (new, old) {
        (Ok(a), Ok(b)) => assert_eq!(bits(&a), bits(&b), "{f_hz} Hz, {ports} ports"),
        (Err(a), Err(b)) => assert_eq!(format!("{a:?}"), format!("{b:?}"), "{f_hz} Hz, {ports} ports"),
        (a, b) => panic!("{f_hz} Hz, {ports} ports: new {a:?} against frozen dense {b:?}"),
    }
}

#[test]
fn the_tridiagonal_solve_reproduces_the_dense_port_admittance_bit_for_bit() {
    for duct in [duct(), eight_hole_duct(), fourteen_hole_duct()] {
        for f in [90.0, 150.0, 300.0, 520.0, 700.0, 1500.0, 2600.0, 4000.0, 9000.0, 15000.0] {
            for ports in [1, 2, 3, 5, 8, 14] {
                assert_same(&duct, f, ports);
            }
        }
    }
}

/// A bore of the size of the A5 clarinet card: 113 segments, 24 tone holes of which 10 are open.
fn clarinet_like_duct() -> Duct {
    let mut segments = Vec::new();
    let mut holes = 0;
    for s in 0..89 {
        let r = |i: usize| 0.0072 + 0.0008 * i as f64 / 89.0;
        if s % 7 == 0 {
            segments.push(Segment::Cylinder { radius: r(s), length: 0.011 });
        } else {
            segments.push(Segment::Cone { inlet_radius: r(s), outlet_radius: r(s + 1), length: 0.011 });
        }
        if holes < 24 && (s % 4 == 3 || s >= 80) {
            let state = if holes % 5 < 2 { HoleState::Open } else { HoleState::Closed };
            segments.push(Segment::ToneHole {
                hole_radius: 0.0015 + 0.0001 * (holes % 5) as f64,
                chimney_height: 0.005 + 0.0007 * (holes % 4) as f64,
                bore_radius: r(s + 1),
                state,
            });
            holes += 1;
        }
    }
    assert_eq!((segments.len(), holes), (113, 24));
    Duct { segments }
}

/// The bins of the multiport sweep of a 0.9 s clarinet render (`48000 / 16384` Hz apart), every `stride`-th one.
fn sweep_bins(stride: usize) -> impl Iterator<Item = f64> {
    (1..8192usize).step_by(stride).map(|k| k as f64 * 48_000.0 / 16_384.0)
}

#[test]
fn the_clarinet_sized_ladder_matches_the_dense_one_on_sampled_sweep_bins() {
    let duct = clarinet_like_duct();
    for f in sweep_bins(997) {
        for ports in [4, 10] {
            assert_same(&duct, f, ports);
        }
    }
}

/// Every bin of the sweep, one run only (minutes): `cargo test --release -p fs-duct -- --ignored
/// clarinet_sized_ladder_matches_the_dense_one_on_every_sweep_bin`.
#[test]
#[ignore = "8191 bins of a 113-segment bore against the dense oracle; run on its own in release"]
fn clarinet_sized_ladder_matches_the_dense_one_on_every_sweep_bin() {
    let duct = clarinet_like_duct();
    for f in sweep_bins(1) {
        assert_same(&duct, f, 4);
    }
}
