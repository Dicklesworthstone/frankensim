//! The delayed free-field mutual term: (1) two outlets against the exact
//! j w rho exp(-j k d)/(4 pi d) and its mutual power rho w^2 sinc(kd)/(4 pi c);
//! (2) a linear two-outlet network against a frequency-domain solve with the
//! same matrix; (3) energy from rest; refusals of the cold admission.
use super::*;
use crate::bernoulli_aperture::radiation::{BaffledRadiationLoad, RadiationFlange};
use core::f64::consts::{PI, TAU};
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::RelaxationTerm;
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment};

const RATE: f64 = 288_000.0;
const DT: f64 = 1.0 / RATE;
const RHO: f64 = 1.2041;
const C: f64 = 343.23;
const BAND: f64 = 2800.0;
const HOLE: f64 = 0.004;
const BORE: f64 = 0.009;
const D: f64 = 0.0301;

fn gate() -> fs_exec::CancelGate { fs_exec::CancelGate::new_clock_free() }
fn compact(flange: RadiationFlange) -> BaffledRadiationLoad {
    BaffledRadiationLoad::with_flange(HOLE, RHO, C, DT, BAND, flange, &gate()).unwrap()
}
fn loads() -> [RelaxationImpedanceSpec; 2] {
    [compact(RadiationFlange::Baffled).load(), compact(RadiationFlange::Unflanged).load()]
}
fn two(d: f64) -> Result<OutletMutual, FlueError> {
    OutletMutual::admit(&[2, 3], &[[0.0, 0.0], [d, 0.0]], &loads(), RHO, C, DT, BAND, &gate())
}
fn what(result: Result<OutletMutual, FlueError>) -> &'static str {
    match result { Err(FlueError::Invalid { what }) => what, other => panic!("expected a refusal, got {other:?}") }
}
fn zero() -> C64 { C64::new(0.0, 0.0) }
fn rel(a: C64, b: C64) -> f64 { (a - b).abs() / b.abs() }

/// x = Re(X exp(j w k dt)) over whole periods: X from one DFT bin, k absolute.
fn bin(x: &[f64], omega: f64, start: usize) -> C64 {
    let mut sum = zero();
    for (i, v) in x.iter().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let (s, c) = fs_math::det::sin_cos(omega * (start + i) as f64 * DT);
        sum = sum + C64::new(v * c, -v * s);
    }
    #[allow(clippy::cast_precision_loss)]
    let n = x.len() as f64;
    sum.scale(2.0 / n)
}

/// Bilinear-realized impedance of a one-term compact load, exp(+j w t).
fn self_z(load: &RelaxationImpedanceSpec, omega: f64) -> C64 {
    let term = load.terms()[0];
    let s = C64::new(0.0, 2.0 * fs_math::det::tan(0.5 * omega * DT) / DT);
    (s.scale(term.resistance_pa_s_m3)) / (s + C64::from_re(term.rate_per_s))
}

#[test]
fn two_outlet_kernel_realizes_the_free_field_term_and_its_mutual_power() {
    let mutual = two(D).unwrap();
    let record = mutual.record().clone();
    assert_eq!(record.min_delay_samples, 24, "{record:?}");
    assert!(record.max_band_error < 0.01, "{record:?}");
    assert!(record.min_pivot_ratio > 0.0 && record.grid_points >= 28_799, "{record:?}");
    let pair = mutual.pairs()[0];
    for hz in [500.0, 1125.0, 2250.0] {
        let omega = TAU * hz;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let period = (RATE / hz).round() as usize;
        let (start, steps) = (4096, 4096 + 32 * period);
        let wave = |k: usize| {
            #[allow(clippy::cast_precision_loss)]
            let t = k as f64 * DT;
            fs_math::det::cos(omega * t)
        };
        // (a) only the second outlet moves: the first sees M q_b.
        // (b) the mirror: the second sees the same sequence, bit for bit.
        // (c) both move in phase: the mean mutual power is Re M.
        let mut runs = [mutual.clone(), mutual.clone(), mutual.clone()];
        let mut seen = [vec![], vec![], vec![]];
        let mut power = 0.0;
        let mut out = [0.0; 2];
        for k in 0..steps {
            for (r, run) in runs.iter_mut().enumerate() {
                run.exterior_pressures(&mut out);
                let q = wave(k);
                let flows = [[0.0, q], [q, 0.0], [q, q]][r];
                if k >= start {
                    match r {
                        0 => seen[0].push(out[0]),
                        1 => seen[1].push(out[1]),
                        _ => power += out[0] * flows[0] + out[1] * flows[1],
                    }
                }
                if r == 2 && k >= start { seen[2].push(q); }
                run.push_flows(&flows);
            }
        }
        assert_eq!(seen[0], seen[1], "{hz} Hz: reciprocity");
        let h = bin(&seen[0], omega, start) / bin(&seen[2], omega, start);
        let realized = pair.realized(omega, DT);
        let exact = pair.exact(omega, RHO, C);
        assert!(rel(h, realized) < 1e-9, "{hz} Hz: {h:?} vs realized {realized:?}");
        assert!(rel(h, exact) <= record.max_band_error + 1e-9, "{hz} Hz: {h:?} vs exact {exact:?}");
        // Mean of sum p q over whole periods: Re M (each outlet 1/2 Re M).
        #[allow(clippy::cast_precision_loss)]
        let mean = power / seen[2].len() as f64;
        let kd = omega * D / C;
        let sinc = fs_math::det::sin(kd) / kd;
        let free_field = RHO * omega * omega * sinc / (4.0 * PI * C);
        assert!((mean / realized.re - 1.0).abs() < 1e-9, "{hz} Hz: {mean} vs {}", realized.re);
        assert!((mean / free_field - 1.0).abs() < 0.05, "{hz} Hz: {mean} vs {free_field}");
        assert!((exact.re / free_field - 1.0).abs() < 1e-12);
    }
}

fn network() -> WaveguideNetwork {
    let [a, b] = loads();
    let nodes = [NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Relaxation { load: a }, NetworkNode::Relaxation { load: b }];
    let zb = RHO * C / (PI * BORE * BORE);
    let zh = RHO * C / (PI * HOLE * HOLE);
    let segs = [
        NetworkSegment { nodes: [0, 1], one_way_samples: 37, impedance_pa_s_m3: zb },
        NetworkSegment { nodes: [1, 2], one_way_samples: 13, impedance_pa_s_m3: zh },
        NetworkSegment { nodes: [1, 3], one_way_samples: 21, impedance_pa_s_m3: zh },
    ];
    let mut net = WaveguideNetwork::new(&nodes, &segs, DT, 1 << 22).unwrap();
    net.attach_terminal_exterior_pressures(&[2, 3]).unwrap();
    net
}

fn solve(mut a: Vec<Vec<C64>>, mut b: Vec<C64>) -> Vec<C64> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs())).unwrap();
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let factor = a[row][col] / a[col][col];
            for k in col..n { let t = a[col][k]; a[row][k] = a[row][k] - factor * t; }
            let t = b[col];
            b[row] = b[row] - factor * t;
        }
    }
    let mut x = vec![zero(); n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for k in row + 1..n { s = s - a[row][k] * x[k]; }
        x[row] = s / a[row][row];
    }
    x
}

/// Steady outlet flows [U_A, U_B] for inlet wave X = 1 with mutual term `m`.
/// Unknowns: waves leaving the junction into lines 1, 2 (o10, o20), waves
/// leaving the outlets (o11, o21), and the outlet flows.
#[allow(clippy::many_single_char_names)]
fn frequency_domain(omega: f64, m: C64) -> [C64; 2] {
    let delay = |n: f64| { let (s, c) = fs_math::det::sin_cos(omega * n * DT); C64::new(c, -s) };
    let (z0, z1, z2) = (delay(37.0), delay(13.0), delay(21.0));
    let zb = RHO * C / (PI * BORE * BORE);
    let zh = RHO * C / (PI * HOLE * HOLE);
    let (y0, yh) = (1.0 / zb, 1.0 / zh);
    let (w0, wh) = (y0 / (y0 + 2.0 * yh), yh / (y0 + 2.0 * yh));
    let [la, lb] = loads();
    let (sa, sb) = (self_z(&la, omega), self_z(&lb, omega));
    let r = |x: f64| C64::from_re(x);
    let o = zero();
    let a = vec![
        vec![r(1.0), o, z1.scale(1.0 - 2.0 * wh), z2.scale(-2.0 * wh), o, o],
        vec![o, r(1.0), z1.scale(-2.0 * wh), z2.scale(1.0 - 2.0 * wh), o, o],
        vec![-z1, o, r(1.0), o, r(zh), o],
        vec![o, -z2, o, r(1.0), o, r(zh)],
        vec![z1.scale(2.0), o, o, o, -(r(zh) + sa), -m],
        vec![o, z2.scale(2.0), o, o, -m, -(r(zh) + sb)],
    ];
    let src = z0.scale(2.0 * w0);
    let x = solve(a, vec![src, src, o, o, o, o]);
    [x[4], x[5]]
}

#[test]
fn a_linear_two_outlet_network_matches_its_frequency_domain_solve() {
    for hz in [500.0, 1125.0, 2250.0] {
        let omega = TAU * hz;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let period = (RATE / hz).round() as usize;
        let (start, steps) = (8192, 8192 + 32 * period);
        let mut net = network();
        let mut mutual = OutletMutual::admit(&[2, 3], &[[0.0, 0.0], [D, 0.0]], &loads(), RHO, C, DT, BAND, &gate()).unwrap();
        let (mut ua, mut ub) = (vec![], vec![]);
        for k in 0..steps {
            mutual.apply(&mut net).unwrap();
            #[allow(clippy::cast_precision_loss)]
            let x = fs_math::det::cos(omega * k as f64 * DT);
            let frame = net.step(x).unwrap();
            mutual.record_step(&net);
            assert!(frame.balance_residual_j().abs()
                <= 1e-11 * (frame.stored_energy_j + frame.inlet_work_j.abs() + 1e-30), "{hz} Hz step {k}");
            if k >= start {
                ua.push(net.node_frame(2).unwrap().net_flow_into_node_m3_s);
                ub.push(net.node_frame(3).unwrap().net_flow_into_node_m3_s);
            }
        }
        let simulated = [bin(&ua, omega, start), bin(&ub, omega, start)];
        let pair = mutual.pairs()[0];
        let coupled = frequency_domain(omega, pair.realized(omega, DT));
        let alone = frequency_domain(omega, zero());
        for i in 0..2 {
            assert!(rel(simulated[i], coupled[i]) < 1e-6, "{hz} Hz outlet {i}: {:?} vs {:?}", simulated[i], coupled[i]);
            // Planted negative: the mutual term is not negligible here.
            assert!(rel(alone[i], coupled[i]) > 1e-3, "{hz} Hz outlet {i}: the coupling vanished");
        }
    }
}

#[test]
fn from_rest_the_exterior_never_returns_more_than_it_took() {
    let mut net = network();
    let mut mutual = OutletMutual::admit(&[2, 3], &[[0.0, 0.0], [D, 0.0]], &loads(), RHO, C, DT, BAND, &gate()).unwrap();
    let (mut supplied, mut taken, mut peak) = (0.0f64, 0.0f64, 0.0f64);
    for k in 0..40_000 {
        mutual.apply(&mut net).unwrap();
        #[allow(clippy::cast_precision_loss)]
        let x = if k < 400 { 40.0 * fs_math::det::sin(0.11 * k as f64) * fs_math::det::sin(PI * k as f64 / 400.0) } else { 0.0 };
        let frame = net.step(x).unwrap();
        mutual.record_step(&net);
        supplied += frame.inlet_work_j;
        for node in [2, 3] {
            let o = net.node_frame(node).unwrap();
            taken += o.pressure_pa * o.net_flow_into_node_m3_s * DT;
        }
        peak = peak.max(frame.stored_energy_j);
        let scale = peak + supplied.abs() + 1e-30;
        assert!(frame.balance_residual_j().abs() <= 1e-11 * scale, "step {k}: {frame:?}");
        assert!(taken >= -1e-11 * scale, "step {k}: the outlets returned {taken} J net");
        // Waves hold what the inlet supplied less what the outlets took.
        assert!(frame.wave_stored_energy_j <= supplied + 1e-11 * scale, "step {k}: wave energy above supplied");
    }
    assert!(peak > 0.0);
    assert!(net.stored_energy_j() < 1e-9 * peak, "{} of {peak}", net.stored_energy_j());
}

#[test]
fn the_cold_admission_refuses_what_it_cannot_certify() {
    assert!(what(two(1.5 * C * DT)).contains("explicit"));
    assert!(what(two(0.0)).contains("distinct centers"));
    assert!(what(OutletMutual::admit(&[2], &[[0.0; 2]], &loads()[..1], RHO, C, DT, BAND, &gate())).contains("two to 64"));
    assert!(what(OutletMutual::admit(&[2, 2], &[[0.0; 2], [D, 0.0]], &loads(), RHO, C, DT, BAND, &gate())).contains("distinct terminals"));
    assert!(what(OutletMutual::admit(&[2, 3], &[[0.0; 2], [D, 0.0]], &loads(), RHO, C, DT, 0.5 / DT, &gate())).contains("Nyquist"));
    // A band where the central difference is far from d/dt fails the 5 % gate.
    assert!(what(OutletMutual::admit(&[2, 3], &[[0.0; 2], [D, 0.0]], &loads(), RHO, C, DT, 60_000.0, &gate())).contains("5 %"));
    // Self loads too weak for the mutual term at high frequency: no Cholesky factor.
    let base = SeriesImpedanceSpec { resistance_pa_s_m3: 0.0, inertance_pa_s2_m3: 0.0, compliance_m3_pa: None };
    let weak = RelaxationImpedanceSpec::new(base, &[RelaxationTerm { resistance_pa_s_m3: 1.0e3, rate_per_s: 1.0e5 }]).unwrap();
    assert!(what(OutletMutual::admit(&[2, 3], &[[0.0; 2], [D, 0.0]], &[weak, weak], RHO, C, DT, BAND, &gate())).contains("not passive"));
    let [a, b] = loads();
    let exact = a.with_sampling(RelaxationSampling::ExactHeldFlow);
    assert!(what(OutletMutual::admit(&[2, 3], &[[0.0; 2], [D, 0.0]], &[exact, b], RHO, C, DT, BAND, &gate())).contains("bilinear"));
}
