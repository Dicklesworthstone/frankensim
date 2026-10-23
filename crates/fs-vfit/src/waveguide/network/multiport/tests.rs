use super::*;
use crate::mimo::check_mimo_passivity;
use crate::model::PoleTerm;
use crate::waveguide::network::NetworkSegment;
use fs_math::c64::C64;

const DT: f64 = 1.0 / 48_000.0;

fn graph() -> WaveguideNetwork {
    WaveguideNetwork::new(
        &[
            NetworkNode::Inlet,
            NetworkNode::Junction,
            NetworkNode::Termination { reflection: -1.0 },
            NetworkNode::Termination { reflection: -1.0 },
        ],
        &[
            NetworkSegment {
                nodes: [0, 1],
                one_way_samples: 1,
                impedance_pa_s_m3: 2.0,
            },
            NetworkSegment {
                nodes: [1, 2],
                one_way_samples: 1,
                impedance_pa_s_m3: 3.0,
            },
            NetworkSegment {
                nodes: [1, 3],
                one_way_samples: 7,
                impedance_pa_s_m3: 4.0,
            },
        ],
        DT,
        1 << 20,
    )
    .unwrap()
}
fn model() -> MimoModel {
    MimoModel {
        ports: 2,
        poles: vec![],
        residues: vec![],
        d: vec![2.0, 0.5, 0.5, 3.0],
        e: vec![0.0; 4],
        reciprocal: true,
    }
}
fn spec(model: MimoModel) -> MimoTerminalSpec {
    MimoTerminalSpec {
        nodes: vec![2, 3],
        model,
        band_rad_s: [1.0, core::f64::consts::PI / DT],
        sampled_passivity_points: 129,
        max_states: 32,
        max_memory_bytes: 1 << 20,
    }
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-12 * (a.abs() + b.abs() + 1e-20),
        "{a} != {b}"
    );
}
fn snapshot(graph: &WaveguideNetwork) -> Vec<u64> {
    let mut out = vec![graph.accepted_steps];
    for line in &graph.lines {
        out.push(line.head as u64);
        out.extend(line.waves.iter().chain(&line.energy).map(|x| x.to_bits()));
    }
    for z in graph.mimo_pole_states().unwrap() {
        out.extend([z.re.to_bits(), z.im.to_bits()]);
    }
    for f in graph.mimo_port_frames().unwrap() {
        out.extend([
            f.pressure_pa.to_bits(),
            f.flow_m3_s.to_bits(),
            f.work_j.to_bits(),
        ]);
    }
    out
}

#[test]
fn mutual_pressure_reacts_before_the_other_bore_branch_arrives() {
    let mut net = graph().with_mimo_terminals(spec(model())).unwrap();
    net.step(1.0).unwrap();
    net.step(0.0).unwrap();
    let a = net.arriving(net.ports[net.offsets[2]]);
    assert_eq!(net.arriving(net.ports[net.offsets[3]]), 0.0);
    assert!(a > 0.0);
    let f = net.step(0.0).unwrap();
    // Independent 2x2 inverse of diag(3,4)+[[2,.5],[.5,3]].
    let u = [2.0 * a * 7.0 / 34.75, -a / 34.75];
    let p = [2.0 * a - 3.0 * u[0], -4.0 * u[1]];
    for i in 0..2 {
        let port = net.mimo_port_frames().unwrap()[i];
        near(port.flow_m3_s, u[i]);
        near(port.pressure_pa, p[i]);
        near(port.work_j, p[i] * u[i] * DT);
    }
    assert!(net.node_frame(3).unwrap().pressure_pa > 0.0);
    assert!(net.mimo_port_frames().unwrap()[1].work_j < 0.0);
    assert!(f.coupled_terminal_work_j > 0.0);
    assert_eq!(f.terminal_loss_j, 0.0);
    assert!(f.balance_residual_j().abs() < 1e-18);
    assert!(net.set_terminal_reflection(2, 0.0).is_err());
}

#[test]
fn rational_states_match_the_original_zoh_and_the_wave_work_window() {
    let mut m = model();
    m.poles = vec![
        PoleTerm::Real {
            pole: -400.0,
            residue: 0.0,
        },
        PoleTerm::Pair {
            pole: C64::new(-900.0, 2000.0),
            residue: C64::ZERO,
        },
    ];
    m.residues = vec![
        vec![
            C64::from_re(200.0),
            C64::from_re(50.0),
            C64::from_re(50.0),
            C64::from_re(180.0),
        ],
        vec![
            C64::new(20.0, 10.0),
            C64::new(2.0, 1.0),
            C64::new(2.0, 1.0),
            C64::new(15.0, 5.0),
        ],
    ];
    let mut reference = m.zoh(DT).unwrap();
    let mut net = graph().with_mimo_terminals(spec(m)).unwrap();
    let mut total_work = 0.0;
    let mut total_inlet = 0.0;
    for sample in 0..256 {
        let old = snapshot(&net);
        let input = if sample < 37 {
            f64::from(sample % 7) * 0.1 - 0.3
        } else {
            0.0
        };
        let preview = net.preview_step(input).unwrap();
        assert_eq!(snapshot(&net), old);
        let frame = net.step(input).unwrap();
        assert_eq!(frame, preview);
        let ports = net.mimo_port_frames().unwrap();
        let expected = reference.step_mean(&[ports[0].flow_m3_s, ports[1].flow_m3_s]);
        for i in 0..2 {
            near(expected[i], ports[i].pressure_pa);
        }
        assert_eq!(reference.pole_states(), net.mimo_pole_states().unwrap());
        total_work += frame.coupled_terminal_work_j;
        total_inlet += frame.inlet_work_j;
        assert!(frame.balance_residual_j().abs() < 1e-17);
    }
    near(net.stored_energy_j() + total_work, total_inlet);
    let before = snapshot(&net);
    assert!(net.step(f64::MAX).is_err());
    assert_eq!(snapshot(&net), before);
    let mut other = graph();
    assert!(net.step_pair(&mut other, [1.0, f64::NAN]).is_err());
    assert_eq!(snapshot(&net), before);
    assert!(net.install_mimo_terminals(spec(model())).is_err());
}

#[test]
fn sampled_response_is_the_actual_start_sample_transfer() {
    let mut m = model();
    m.poles = vec![PoleTerm::Real {
        pole: -400.0,
        residue: 0.0,
    }];
    m.residues = vec![vec![
        C64::from_re(200.0),
        C64::from_re(50.0),
        C64::from_re(50.0),
        C64::from_re(180.0),
    ]];
    let discrete = m.zoh(DT).unwrap();
    let omega = 2.0 * core::f64::consts::PI * 173.0;
    let (s, c) = fs_math::det::sin_cos(omega * DT);
    let z = C64::new(c, s);
    let decay = fs_math::det::exp(-400.0 * DT);
    let response = discrete.eval_discrete_iw(omega);
    let mean_response = discrete.eval_mean_discrete_iw(omega);
    let phi1 = (1.0 - decay) / (400.0 * DT);
    for i in 0..4 {
        let expected = C64::from_re(m.d[i])
            + (z - C64::from_re(decay))
                .recip()
                .scale(m.residues[0][i].re * (1.0 - decay) / 400.0);
        near(response[i].re, expected.re);
        near(response[i].im, expected.im);
        let mean_expected = C64::from_re(m.d[i] + m.residues[0][i].re * (1.0 - phi1) / 400.0)
            + (expected - C64::from_re(m.d[i])).scale(phi1);
        near(mean_response[i].re, mean_expected.re);
        near(mean_response[i].im, mean_expected.im);
    }
}

#[test]
fn exact_mean_exchange_admits_a_passive_model_rejected_by_start_sampling() {
    let mut m = model();
    m.d = vec![1.0, 0.0, 0.0, 1.0];
    m.poles = vec![PoleTerm::Real {
        pole: -1.0e6,
        residue: 0.0,
    }];
    m.residues = vec![vec![
        C64::from_re(1.0e8),
        C64::ZERO,
        C64::ZERO,
        C64::from_re(1.0e8),
    ]];
    assert!(
        check_mimo_passivity(&m, (1.0, core::f64::consts::PI / DT))
            .unwrap()
            .passive
    );
    assert!(!m.zoh(DT).unwrap().sampled_passivity(129).passive);
    assert!(m.zoh(DT).unwrap().mean_sampled_passivity(129).passive);
    let mut net = graph();
    assert!(net.install_mimo_terminals(spec(m)).is_ok());
    assert!(net.mimo_terminal_spec().is_some());
    assert_eq!(net.stored_energy_j(), 0.0);
}

#[test]
fn exact_mean_pressure_closes_analytic_one_pole_work_at_small_and_stiff_steps() {
    // G1: x'=-a*x+u, p=D*u+R*x. E=R*x^2/2, loss=D*u^2+a*R*x^2.
    // No fitted energy accumulator supplies any term in this independent audit.
    let a = 400.0_f64;
    let r = 500.0_f64;
    let mut model = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Real {
            pole: -a,
            residue: 0.0,
        }],
        residues: vec![vec![C64::from_re(r)]],
        d: vec![0.0],
        e: vec![0.0],
        reciprocal: true,
    };
    for (dt, direct) in [(1.0e-10, 0.0), (0.002, 0.25), (0.05, 0.0)] {
        model.d[0] = direct;
        let mut exact = model.zoh(dt).unwrap();
        let mut start_output = model.zoh(dt).unwrap();
        let mut x = 0.0_f64;
        for u in [1.0, -0.4, 0.0, 0.8, 0.0] {
            let h = a * dt;
            let exponential = (-h).exp();
            let integral_exp = -(-h).exp_m1() / a;
            let next = exponential * x + u * integral_exp;
            let equilibrium = u / a;
            let offset = x - equilibrium;
            let (mean_x, integral_x2) = if h < 1.0e-5 {
                // Cancellation-free small-h Taylor mean. Simpson integrates the
                // analytic trajectory's limiting quadratic x(t)^2; its error
                // here is O(h^2) relative and below the stated tolerance.
                let mean = x + (u - a * x) * dt * (0.5 - h / 6.0 + h * h / 24.0);
                let half = (-0.5 * h).exp() * x - u * (-0.5 * h).exp_m1() / a;
                (mean, dt * (x * x + 4.0 * half * half + next * next) / 6.0)
            } else {
                // Closed-form integral of (equilibrium + offset*exp(-a*t))^2.
                (
                    equilibrium + offset * integral_exp / dt,
                    equilibrium * equilibrium * dt
                        + 2.0 * equilibrium * offset * integral_exp
                        + offset * offset * -(-2.0 * h).exp_m1() / (2.0 * a),
                )
            };
            let pressure = exact.step_mean(&[u])[0];
            near(pressure, direct * u + r * mean_x);
            near(start_output.step(&[u])[0], direct * u + r * x);
            assert_eq!(exact.pole_states(), start_output.pole_states());
            near(exact.pole_states()[0].re, next);
            let before = 0.5 * r * x * x;
            let after = 0.5 * r * next * next;
            let loss = direct * u * u * dt + a * r * integral_x2;
            let work = pressure * u * dt;
            let scale = before.abs() + after.abs() + loss.abs() + work.abs() + f64::MIN_POSITIVE;
            assert!(
                (after - before + loss - work).abs() <= 5.0e-12 * scale,
                "dt={dt}, x={x}, u={u}, delta={}, loss={loss}, work={work}",
                after - before
            );
            x = next;
        }
    }
}

#[test]
fn unsupported_physics_and_false_source_shapes_refuse() {
    let mut bad = spec(model());
    bad.band_rad_s[1] *= 1.0 - 1e-9;
    assert!(graph().with_mimo_terminals(bad).is_err());
    let mut bad = spec(model());
    bad.model.e[0] = 1.0;
    assert!(graph().with_mimo_terminals(bad).is_err());
    let mut bad = spec(model());
    bad.nodes[1] = bad.nodes[0];
    assert!(graph().with_mimo_terminals(bad).is_err());
    let mut bad = spec(model());
    bad.model.d.pop();
    assert!(graph().with_mimo_terminals(bad).is_err());
    let mut bad = spec(model());
    bad.max_memory_bytes = 0;
    assert!(graph().with_mimo_terminals(bad).is_err());
    let mut net = graph();
    net.step(0.0).unwrap();
    assert!(net.install_mimo_terminals(spec(model())).is_err());
}

#[test]
fn network_and_field_consumers_share_the_same_load_admission() {
    let source = spec(model());
    let limits = MimoLoadLimits {
        band_rad_s: source.band_rad_s,
        sampled_passivity_points: source.sampled_passivity_points,
        max_states: source.max_states,
        max_memory_bytes: source.max_memory_bytes,
    };
    let prepared = PreparedMimoLoad::new(&source.model, DT, limits).unwrap();
    let network = graph().with_mimo_terminals(source).unwrap();
    let (continuous, sampled) = network.mimo_passivity().unwrap();
    assert_eq!(continuous.class, prepared.reports().0.class);
    assert_eq!(continuous.worst, prepared.reports().0.worst);
    assert_eq!(sampled.class, prepared.reports().1.class);
    assert_eq!(sampled.worst, prepared.reports().1.worst);
    assert_eq!(network.mimo_pole_states().unwrap(), prepared.pole_states());
    for bad_limits in [
        MimoLoadLimits { max_memory_bytes: 0, ..limits },
        MimoLoadLimits { sampled_passivity_points: 16, ..limits },
        MimoLoadLimits { band_rad_s: [1.0, 0.9 * limits.band_rad_s[1]], ..limits },
    ] {
        let mut source = spec(model());
        source.band_rad_s = bad_limits.band_rad_s;
        source.sampled_passivity_points = bad_limits.sampled_passivity_points;
        source.max_states = bad_limits.max_states;
        source.max_memory_bytes = bad_limits.max_memory_bytes;
        let direct = PreparedMimoLoad::new(&source.model, DT, bad_limits).unwrap_err();
        let network = graph().with_mimo_terminals(source).err().unwrap();
        assert_eq!(network.0, direct.0);
    }
}
