use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, ApertureTerminal, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::loss::fit_boundary_loss;
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, ApertureNetworkFrame, NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn cavity() -> HelmholtzLoadSpec {
    HelmholtzLoadSpec { volume_m3: 1e-4, neck_radius_m: 0.003,
        effective_neck_length_m: 0.02, resistance_pa_s_m3: 2e5 }
}
fn excess(w: f64, scale: f64) -> f64 {
    // Independent synthetic target; no measured material authority is claimed.
    scale * (2e5 * w * w / (w * w + 500.0_f64.powi(2))
        + 3e5 * w * w / (w * w + 6000.0_f64.powi(2)))
}
fn terminal(scale: f64) -> NetworkNode {
    let base = cavity().impedance(1.2, 343.0).unwrap();
    if scale == 0.0 { return NetworkNode::Impedance { load: base }; }
    let train = [500.0, 6000.0].map(|w| (w, excess(w, scale)));
    let checks = [1000.0, 2000.0, 4000.0].map(|w| (w, excess(w, scale)));
    fit_boundary_loss(base, &train, &checks, 1e-7, 1e-10).unwrap().termination()
}
fn model(scale: f64, budget: u64) -> ApertureNetwork {
    let cell = 343.0 * 1e-5;
    let sections = [(0, 1, 8, 0.007), (1, 2, 12, 0.009), (1, 3, 5, 0.003)]
        .map(|(a, b, n, radius)| TubeSection { nodes: [a, b], length_m: f64::from(n) * cell,
            radius_m: radius, max_length_error_m: 1e-12 }).to_vec();
    let spec = TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.8 }, terminal(scale)],
        sections, sound_speed_m_s: 343.0, max_wave_memory_bytes: 1 << 20 };
    let mechanics = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35, density_kg_m3: 1.2,
        impedance_pa_s_m3: spec.inlet_impedance(1.2).unwrap(), time_step_s: 1e-5, max_steps: budget,
    };
    let obstacle = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic contact fixture".into()).unwrap().with_internal_loss(5.0).unwrap();
    ApertureNetwork::new(DynamicAperture::new(mechanics, ApertureState {
        opening_m: 4e-4, opening_velocity_m_s: 0.0 }, obstacle).unwrap(), spec).unwrap()
}
fn frame_bits(f: ApertureNetworkFrame) -> Vec<u64> {
    [f.aperture.state.opening_m, f.aperture.state.opening_velocity_m_s,
     f.aperture.outgoing_pressure_pa, f.aperture.bore_pressure_pa,
     f.stored_energy_j, f.storage_change_j, f.dissipated_energy_j,
     f.network.wave_stored_energy_j, f.network.load_stored_energy_j,
     f.upstream_work_j, f.body_work_j].map(f64::to_bits).to_vec()
}

#[test]
fn fit_retains_independent_check_error_and_rejects_inaccurate_positive_fallback() {
    let base = cavity().impedance(1.2, 343.0).unwrap();
    let train = [500.0, 6000.0].map(|w| (w, excess(w, 1.0)));
    let checks = [1000.0, 2000.0, 4000.0].map(|w| (w, excess(w, 1.0)));
    let fit = fit_boundary_loss(base, &train, &checks, 1e-7, 1e-10).unwrap();
    assert!(fit.max_checked_error_pa_s_m3 < 1e-7);
    assert!((fit.load.terms()[0].resistance_pa_s_m3 - 2e5).abs() < 1e-6);
    assert!((fit.load.terms()[1].resistance_pa_s_m3 - 3e5).abs() < 1e-6);
    let mut corrupt = checks;
    corrupt[1].1 *= 1.2;
    assert!(fit_boundary_loss(base, &train, &corrupt, 1e-7, 1e-10).is_err());
    // Positive Foster terms cannot realize a decreasing resistance curve.
    assert!(fit_boundary_loss(base, &[(500.0, 1e6), (6000.0, 1e3)], &[(1000.0, 5e5)], 1.0, 0.001).is_err());
}

#[test]
fn changed_loss_changes_actual_valve_motion_only_after_the_return_path() {
    let (mut a, mut b) = (model(0.5, 1024), model(2.0, 1024));
    let (mut dp, mut dy) = (0.0_f64, 0.0_f64);
    let drive = TubeDrive { upstream_pressure_pa: 800.0, body_flow_m3_s: 0.0 };
    for n in 0..1024 {
        let x = a.step(drive).unwrap();
        let y = b.step(drive).unwrap();
        if n < 26 {
            assert_eq!(x.aperture.state.opening_m.to_bits(), y.aperture.state.opening_m.to_bits());
            assert_eq!(x.aperture.bore_pressure_pa.to_bits(), y.aperture.bore_pressure_pa.to_bits());
        }
        dp = dp.max((x.aperture.bore_pressure_pa - y.aperture.bore_pressure_pa).abs());
        dy = dy.max((x.aperture.state.opening_m - y.aperture.state.opening_m).abs());
    }
    assert!(dp > 0.1 && dy > 1e-9, "loss must affect physical feedback, not just an output filter");
}

#[test]
fn full_balance_includes_loss_memory_and_unforced_network_decays() {
    let mut model = model(1.0, 1024);
    let mut missing_storage_error = 0.0_f64;
    let mut before_load = 0.0;
    for n in 0..1024 {
        let before = model.stored_energy_j();
        let f = model.step(TubeDrive { upstream_pressure_pa: if n < 256 { 1200.0 } else { 0.0 },
            body_flow_m3_s: if n < 400 { 2e-7 * (0.1 * f64::from(n)).sin() } else { 0.0 } }).unwrap();
        let scale = (before + f.stored_energy_j + f.dissipated_energy_j
            + f.upstream_work_j.abs() + f.body_work_j.abs()).max(f64::MIN_POSITIVE);
        assert!(f.balance_residual_j().abs() < 3e-10 * scale);
        assert!(f.dissipated_energy_j >= 0.0);
        if n >= 400 { assert!(f.stored_energy_j <= before + 3e-10 * scale); }
        let delta = f.network.load_stored_energy_j - before_load;
        missing_storage_error = missing_storage_error.max((f.balance_residual_j() - delta).abs());
        before_load = f.network.load_stored_energy_j;
    }
    assert!(missing_storage_error > 1e-12);
}

#[test]
fn cancellation_budget_and_refusal_retain_wave_valve_and_loss_histories() {
    let inputs: Vec<_> = (0..512).map(|n| TubeDrive {
        upstream_pressure_pa: if n < 256 { 800.0 } else { 0.0 },
        body_flow_m3_s: 2e-7 * (0.1 * f64::from(n)).sin(),
    }).collect();
    let mut a = model(1.0, 512);
    let expected: Vec<_> = inputs.iter().map(|d| frame_bits(a.step(*d).unwrap())).collect();
    let mut b = model(1.0, 137);
    let mut actual: Vec<_> = inputs[..100].iter().map(|d| frame_bits(b.step(*d).unwrap())).collect();
    let gate = CancelGate::new_clock_free();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let mut out = vec![ApertureNetworkFrame::default(); 412];
    let p = b.advance_block(&inputs[100..], &mut out, &cancelled).unwrap();
    assert_eq!((p.completed, p.terminal), (0, ApertureTerminal::Cancelled));
    assert!(out.iter().all(|f| f.aperture.step == 0));
    let p = b.advance_block(&inputs[100..], &mut out, &gate).unwrap();
    assert_eq!((p.completed, p.terminal), (37, ApertureTerminal::BudgetExhausted));
    actual.extend(out[..37].iter().map(|f| frame_bits(*f)));
    assert!(out[37..].iter().all(|f| f.aperture.step == 0));
    b.extend_step_budget(512).unwrap();
    actual.extend(inputs[137..].iter().map(|d| frame_bits(b.step(*d).unwrap())));
    assert_eq!(actual, expected);
    b.extend_step_budget(600).unwrap();
    let state = b.aperture().state();
    let node = *b.node_frame(3).unwrap();
    let energy = b.stored_energy_j();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(b.step(TubeDrive { upstream_pressure_pa: 800.0, body_flow_m3_s: bad }).is_err());
    }
    assert!(b.set_terminal_reflection(3, 0.0).is_err());
    assert_eq!(b.aperture().state(), state);
    assert_eq!(*b.node_frame(3).unwrap(), node);
    assert_eq!(b.stored_energy_j().to_bits(), energy.to_bits());
}

#[test]
fn table_admission_refuses_invalid_scale_reused_checks_and_nonfinite_values() {
    let base = cavity().impedance(1.2, 343.0).unwrap();
    let train = [(500.0, excess(500.0, 1.0)), (6000.0, excess(6000.0, 1.0))];
    assert!(fit_boundary_loss(base, &train, &train[..1], 1.0, 0.01).is_err());
    assert!(fit_boundary_loss(base, &train, &[(100.0, 1.0)], 1.0, 0.01).is_err());
    assert!(fit_boundary_loss(base, &[(1e200, 1.0), (2e200, 1.0)], &[(1.5e200, 1.0)], 1.0, 0.01).is_err());
    assert!(fit_boundary_loss(base, &train, &[(1000.0, f64::NAN)], 1.0, 0.01).is_err());
    assert!(fit_boundary_loss(base, &train, &[(1000.0, 1.0)], f64::INFINITY, 0.01).is_err());
}
