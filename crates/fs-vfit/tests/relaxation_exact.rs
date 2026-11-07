//! Physical branch trajectories and whole-network work for exact memory.

use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{
    ComplianceRelaxationTerm, RelaxationImpedance, RelaxationImpedanceSpec, RelaxationSampling,
    RelaxationTerm,
};
use fs_vfit::waveguide::network::admittance::{
    AdmittanceTerm, RelaxationAdmittance, RelaxationAdmittanceSpec,
};
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};

fn base(r: f64, l: f64, c: Option<f64>) -> SeriesImpedanceSpec {
    SeriesImpedanceSpec {
        resistance_pa_s_m3: r,
        inertance_pa_s2_m3: l,
        compliance_m3_pa: c,
    }
}
fn spec() -> RelaxationImpedanceSpec {
    RelaxationImpedanceSpec::new(
        base(0.4, 0.3, Some(0.7)),
        &[RelaxationTerm {
            resistance_pa_s_m3: 0.6,
            rate_per_s: 2.0,
        }],
    )
    .unwrap()
    .with_compliance_terms(&[
        ComplianceRelaxationTerm {
            resistance_pa_s_m3: 0.8,
            rate_per_s: 3.0,
        },
        ComplianceRelaxationTerm {
            resistance_pa_s_m3: 0.2,
            rate_per_s: 40.0,
        },
    ])
    .unwrap()
    .with_sampling(RelaxationSampling::ExactHeldFlow)
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * scale.max(f64::MIN_POSITIVE),
        "{actual:.17e} != {expected:.17e}, scale {scale:.17e}"
    );
}

// The independent analytic interval-mean pressure needed to impose a held flow.
fn incident_for_flow(load: &RelaxationImpedance, z: f64, dt: f64, flow: f64) -> f64 {
    let spec = load.spec();
    let base = spec.base();
    let state = load.base_state();
    let mut pressure = base.resistance_pa_s_m3 * flow
        + 2.0 * base.inertance_pa_s2_m3 / dt * (flow - state.inertive_flow_m3_s)
        + state.compliance_pressure_pa
        + base.compliance_m3_pa.map_or(0.0, |c| 0.5 * dt / c * flow);
    for (&q, term) in load.branch_flows_m3_s().iter().zip(spec.terms()) {
        let r = dt * term.rate_per_s;
        let mean = match spec.sampling() {
            RelaxationSampling::ExactHeldFlow => -(-r).exp_m1() / r,
            RelaxationSampling::BilinearMidpoint => 1.0 / (1.0 + 0.5 * r),
        };
        pressure += term.resistance_pa_s_m3 * mean * (flow - q);
    }
    for (&p, term) in load
        .branch_pressures_pa()
        .iter()
        .zip(spec.compliance_terms())
    {
        let r = dt * term.rate_per_s;
        let (old, driven) = match spec.sampling() {
            RelaxationSampling::ExactHeldFlow => {
                let b = -(-r).exp_m1() / r;
                (b, 1.0 - b)
            }
            RelaxationSampling::BilinearMidpoint => {
                (1.0 / (1.0 + 0.5 * r), 0.5 * r / (1.0 + 0.5 * r))
            }
        };
        pressure += old * p + driven * term.resistance_pa_s_m3 * flow;
    }
    0.5 * (pressure + z * flow)
}

#[test]
fn exact_mixed_branches_share_the_solved_flow_and_integrate_physical_heat() {
    let (dt, z) = (0.02, 1.0);
    let spec = spec();
    let mut load = RelaxationImpedance::new(spec, z, dt).unwrap();
    for n in 0..180 {
        let old_base = load.base_state();
        let old_q = load.branch_flows_m3_s().to_vec();
        let old_p = load.branch_pressures_pa().to_vec();
        let before = load.stored_energy_j();
        let frame = load
            .step(if n < 120 {
                f64::from(n % 13) - 6.0
            } else {
                0.0
            })
            .unwrap();
        let u = frame.port.flow_m3_s;
        let mut pressure = 0.4 * u
            + 0.3 * (frame.port.state.inertive_flow_m3_s - old_base.inertive_flow_m3_s) / dt
            + 0.5 * (frame.port.state.compliance_pressure_pa + old_base.compliance_pressure_pa);
        let mut energy = 0.15 * frame.port.state.inertive_flow_m3_s.powi(2)
            + 0.35 * frame.port.state.compliance_pressure_pa.powi(2);
        let mut heat = 0.4 * u * u * dt;
        for (i, term) in spec.terms().iter().enumerate() {
            let r = term.rate_per_s * dt;
            let decay = (-r).exp();
            let next = old_q[i] + (-(-r).exp_m1()) * (u - old_q[i]);
            close(
                frame.branch_flows_m3_s()[i],
                next,
                next.abs() + old_q[i].abs() + u.abs(),
            );
            pressure += term.resistance_pa_s_m3 * (-(-r).exp_m1() / r) * (u - old_q[i]);
            let l = term.resistance_pa_s_m3 / term.rate_per_s;
            energy += 0.5 * l * next * next;
            heat += 0.5 * l * (u - old_q[i]).powi(2) * (1.0 - decay * decay);
        }
        for (i, term) in spec.compliance_terms().iter().enumerate() {
            let r = term.rate_per_s * dt;
            let target = term.resistance_pa_s_m3 * u;
            let next = target + (old_p[i] - target) * (-r).exp();
            close(
                frame.branch_pressures_pa()[i],
                next,
                next.abs() + old_p[i].abs() + target.abs(),
            );
            let b = -(-r).exp_m1() / r;
            pressure += b * old_p[i] + (1.0 - b) * target;
            energy += 0.5 * next * next / (term.resistance_pa_s_m3 * term.rate_per_s);
            // Independent quadrature of the analytic resistor pressure, not
            // the Prony heat formula or a work-minus-storage subtraction.
            let mut integral = 0.0;
            let intervals = 400;
            for k in 0..=intervals {
                let t = dt * f64::from(k) / f64::from(intervals);
                let p = target + (old_p[i] - target) * (-term.rate_per_s * t).exp();
                let weight = if k == 0 || k == intervals {
                    1.0
                } else if k % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                integral += weight * p * p / term.resistance_pa_s_m3;
            }
            heat += integral * dt / (3.0 * f64::from(intervals));
        }
        close(
            frame.port.pressure_pa,
            pressure,
            frame.port.pressure_pa.abs() + pressure.abs() + 1.0,
        );
        close(frame.port.stored_energy_j, energy, energy);
        close(frame.port.dissipated_energy_j, heat, heat);
        close(
            frame.port.balance_residual_j(),
            0.0,
            before + energy + heat + frame.port.supplied_work_j.abs(),
        );
    }
}

#[test]
fn fast_rc_memory_relaxes_without_midpoint_ringing_and_weak_heat_is_retained() {
    let spec = RelaxationImpedanceSpec::new(base(0.0, 0.0, None), &[])
        .unwrap()
        .with_compliance_terms(&[ComplianceRelaxationTerm {
            resistance_pa_s_m3: 2.0,
            rate_per_s: 1e6,
        }])
        .unwrap();
    let dt = 0.01;
    let mut exact = RelaxationImpedance::new(
        spec.with_sampling(RelaxationSampling::ExactHeldFlow),
        1.0,
        dt,
    )
    .unwrap();
    let mut midpoint = RelaxationImpedance::new(spec, 1.0, dt).unwrap();
    for load in [&mut exact, &mut midpoint] {
        load.step(incident_for_flow(load, 1.0, dt, 1.0)).unwrap();
    }
    let before = exact.stored_energy_j();
    let frame = exact.step(incident_for_flow(&exact, 1.0, dt, 0.0)).unwrap();
    close(frame.port.flow_m3_s, 0.0, 1.0);
    assert!(frame.branch_pressures_pa()[0].abs() < 1e-12);
    close(frame.port.dissipated_energy_j, before, before);
    close(frame.port.balance_residual_j(), 0.0, before);
    let old = midpoint.branch_pressures_pa()[0];
    let ring = midpoint
        .step(incident_for_flow(&midpoint, 1.0, dt, 0.0))
        .unwrap();
    assert!(ring.branch_pressures_pa()[0] < -0.99 * old);

    let weak = RelaxationImpedanceSpec::new(base(0.0, 0.0, None), &[])
        .unwrap()
        .with_compliance_terms(&[ComplianceRelaxationTerm {
            resistance_pa_s_m3: 1.0,
            rate_per_s: 1e-12,
        }])
        .unwrap()
        .with_sampling(RelaxationSampling::ExactHeldFlow);
    let frame = RelaxationImpedance::new(weak, 1.0, 1.0)
        .unwrap()
        .step(1.0)
        .unwrap();
    let expected = frame.port.flow_m3_s.powi(2) * 1e-24 / 3.0;
    assert!(frame.port.dissipated_energy_j > 0.0);
    close(frame.port.dissipated_energy_j, expected, expected);
}

#[test]
fn held_flow_subdivision_preserves_both_memories_and_their_independent_heat() {
    let spec = RelaxationImpedanceSpec::new(
        base(0.2, 0.0, None),
        &[RelaxationTerm {
            resistance_pa_s_m3: 0.7,
            rate_per_s: 7000.0,
        }],
    )
    .unwrap()
    .with_compliance_terms(&[ComplianceRelaxationTerm {
        resistance_pa_s_m3: 0.3,
        rate_per_s: 30.0,
    }])
    .unwrap()
    .with_sampling(RelaxationSampling::ExactHeldFlow);
    let dt = 0.2;
    let mut coarse = RelaxationImpedance::new(spec, 1.0, dt).unwrap();
    let mut fine = RelaxationImpedance::new(spec, 1.0, dt / 100.0).unwrap();
    for u in [0.4, -0.2, 0.0, 0.0] {
        let frame = coarse.step(incident_for_flow(&coarse, 1.0, dt, u)).unwrap();
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..100 {
            let small = fine
                .step(incident_for_flow(&fine, 1.0, dt / 100.0, u))
                .unwrap();
            heat += small.port.dissipated_energy_j;
            work += small.port.supplied_work_j;
        }
        close(
            coarse.branch_flows_m3_s()[0],
            fine.branch_flows_m3_s()[0],
            1.0,
        );
        close(
            coarse.branch_pressures_pa()[0],
            fine.branch_pressures_pa()[0],
            1.0,
        );
        close(
            frame.port.dissipated_energy_j,
            heat,
            heat + frame.port.stored_energy_j,
        );
        close(
            frame.port.supplied_work_j,
            work,
            heat + frame.port.stored_energy_j + work.abs(),
        );
    }
}

#[test]
fn defaults_empty_branches_and_refused_trials_keep_their_existing_semantics() {
    let mid = spec().with_sampling(RelaxationSampling::BilinearMidpoint);
    let original = RelaxationImpedanceSpec::new(mid.base(), mid.terms())
        .unwrap()
        .with_compliance_terms(mid.compliance_terms())
        .unwrap();
    assert_eq!(original, mid);
    let mut a = RelaxationImpedance::new(original, 1.0, 0.02).unwrap();
    let mut b = RelaxationImpedance::new(mid, 1.0, 0.02).unwrap();
    for n in 0..100 {
        assert_eq!(
            a.step(f64::from(n % 7)).unwrap(),
            b.step(f64::from(n % 7)).unwrap()
        );
    }
    let empty = RelaxationImpedanceSpec::new(mid.base(), &[]).unwrap();
    let mut a = RelaxationImpedance::new(empty, 1.0, 0.02).unwrap();
    let mut b = RelaxationImpedance::new(
        empty.with_sampling(RelaxationSampling::ExactHeldFlow),
        1.0,
        0.02,
    )
    .unwrap();
    for n in 0..100 {
        assert_eq!(
            a.step(f64::from(n % 7)).unwrap(),
            b.step(f64::from(n % 7)).unwrap()
        );
    }
    let mut exact = RelaxationImpedance::new(spec(), 1.0, 0.02).unwrap();
    exact.step(0.7).unwrap();
    let mut untouched = exact;
    for invalid in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(exact.step(invalid).is_err());
    }
    assert_eq!(exact.step(-0.3).unwrap(), untouched.step(-0.3).unwrap());
}

#[test]
fn exact_series_rc_admittance_uses_the_same_pressure_conjugate_physics() {
    let (g, rate, z, dt) = (0.7, 5000.0, 1.3, 0.02);
    let spec = RelaxationAdmittanceSpec::new(
        0.2,
        0.0,
        &[AdmittanceTerm {
            conductance_m3_pa_s: g,
            rate_per_s: rate,
        }],
    )
    .unwrap()
    .with_sampling(RelaxationSampling::ExactHeldFlow);
    assert_eq!(spec.sampling(), RelaxationSampling::ExactHeldFlow);
    let mut load = RelaxationAdmittance::new(spec, z, dt).unwrap();
    for n in 0..200 {
        let old = load.branch_pressures_pa()[0];
        let before = load.stored_energy_j();
        let frame = load
            .step(if n < 100 { f64::from(n % 9) - 4.0 } else { 0.0 })
            .unwrap();
        let expected = frame.pressure_pa + (old - frame.pressure_pa) * (-rate * dt).exp();
        close(
            load.branch_pressures_pa()[0],
            expected,
            old.abs() + frame.pressure_pa.abs(),
        );
        let heat = 0.2 * frame.pressure_pa.powi(2) * dt
            + 0.5 * g / rate * (old - frame.pressure_pa).powi(2) * (-(-2.0 * rate * dt).exp_m1());
        close(frame.dissipated_energy_j, heat, heat);
        close(
            frame.balance_residual_j(),
            0.0,
            before + frame.stored_energy_j + heat + frame.supplied_work_j.abs(),
        );
    }
}

#[test]
fn terminal_shunt_and_series_network_loads_use_exact_memory_in_one_work_balance() {
    for kind in 0..3 {
        let load = spec();
        let (nodes, segments, loaded) = if kind == 0 {
            (
                vec![NetworkNode::Inlet, NetworkNode::Relaxation { load }],
                vec![NetworkSegment {
                    nodes: [0, 1],
                    one_way_samples: 3,
                    impedance_pa_s_m3: 1.0,
                }],
                1,
            )
        } else {
            (
                vec![
                    NetworkNode::Inlet,
                    if kind == 1 {
                        NetworkNode::Shunt { load }
                    } else {
                        NetworkNode::Series { load }
                    },
                    NetworkNode::Termination { reflection: -0.6 },
                ],
                vec![
                    NetworkSegment {
                        nodes: [0, 1],
                        one_way_samples: 3,
                        impedance_pa_s_m3: 1.0,
                    },
                    NetworkSegment {
                        nodes: [1, 2],
                        one_way_samples: 5,
                        impedance_pa_s_m3: 2.0,
                    },
                ],
                1,
            )
        };
        let dt = 0.02;
        let mut network = WaveguideNetwork::new(&nodes, &segments, dt, 1 << 20).unwrap();
        let mut dissipated = 0.0;
        for n in 0..900 {
            let before = network.stored_energy_j();
            let previous = network.load_compliance_pressures(loaded).unwrap().to_vec();
            let frame = network
                .step(if n < 250 {
                    f64::from(n % 17) - 8.0
                } else {
                    0.0
                })
                .unwrap();
            let node = network.node_frame(loaded).unwrap();
            let u = if kind == 0 {
                node.net_flow_into_node_m3_s
            } else {
                node.load_flow_m3_s
            };
            for ((&old, &next), term) in previous
                .iter()
                .zip(network.load_compliance_pressures(loaded).unwrap())
                .zip(load.compliance_terms())
            {
                let target = term.resistance_pa_s_m3 * u;
                let expected = target + (old - target) * (-term.rate_per_s * dt).exp();
                close(next, expected, next.abs() + old.abs() + target.abs());
            }
            dissipated += frame.terminal_loss_j + frame.interior_loss_j;
            close(
                frame.balance_residual_j(),
                0.0,
                before
                    + frame.stored_energy_j
                    + frame.inlet_work_j.abs()
                    + frame.terminal_loss_j
                    + frame.interior_loss_j,
            );
        }
        assert!(dissipated > 0.0);
    }
}
