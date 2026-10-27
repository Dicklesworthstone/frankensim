use super::*;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::flue_loop::network::FlueNetworkSpec;
use crate::flue_loop::{FlueControl, FlueGeometry, JetStepCoefficients, RetardedVelocity};
use fs_material::gas::GasSpec;

fn gas(temperature: f64, pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap()
}
fn model(size: [f64; 2], cells: [usize; 2], gas: GasState) -> RectangularMouthLoad {
    RectangularMouthLoad::new(
        size,
        cells,
        gas,
        1.0 / 96_000.0,
        1000.0,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn actual_rectangle_refines_toward_the_analytic_uniform_aperture_inertance() {
    let g = gas(293.15, 101_325.0);
    let (a, b) = (0.010_f64, 0.004_f64);
    // Independent integral of (a-|x|)(b-|y|)/sqrt(x*x+y*y) over the
    // difference rectangle, not a second copy of the cell/self-disc rule.
    let integral = 2.0 * (a * a * b * (b / a).asinh() + a * b * b * (a / b).asinh())
        + (2.0 / 3.0) * (a.powi(3) + b.powi(3) - a.hypot(b).powi(3));
    let exact = g.density * integral / (core::f64::consts::TAU * (a * b).powi(2));
    let mut previous = f64::INFINITY;
    // Coarse self cells cannot certify themselves against the same biased rule.
    assert!(
        RectangularMouthLoad::new(
            [a, b],
            [4, 4],
            g,
            1.0 / 96_000.0,
            1000.0,
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    for n in [8, 12, 16] {
        let load = model([a, b], [n, n], g);
        let spec = load.load();
        let term = spec.terms()[0];
        let inertia = term.resistance_pa_s_m3 / term.rate_per_s;
        assert!((inertia / exact - 1.0).abs() < 1e-13);
        let error = load.static_quadrature_relative_error();
        assert!(
            error < previous,
            "rectangle quadrature must converge: {n}, {error}"
        );
        previous = error;
        let coefficient = term.resistance_pa_s_m3 / term.rate_per_s.powi(2);
        assert!(
            (coefficient / (g.density / (core::f64::consts::TAU * g.sound_speed)) - 1.0).abs()
                < 1e-13
        );
        assert!(load.max_complex_relative_error() <= 0.05);
        assert!(load.max_resistance_relative_error() <= 0.05);
    }
    assert!(previous < 0.02);
    let squareish = model([a, b], [16, 16], g).load().terms()[0];
    let slender = model([2.0 * a, 0.5 * b], [32, 4], g).load().terms()[0];
    assert!(
        (squareish.resistance_pa_s_m3 / squareish.rate_per_s
            - slender.resistance_pa_s_m3 / slender.rate_per_s)
            .abs()
            > 0.1 * exact,
        "equal area is not equal radiation inertance; retain the physical aspect ratio"
    );
}

#[test]
fn exterior_density_and_sound_speed_recompute_the_passive_coefficients() {
    let base = model([0.01, 0.004], [8, 4], gas(293.15, 101_325.0));
    let dense = model([0.01, 0.004], [8, 4], gas(293.15, 202_650.0));
    let warm = model([0.01, 0.004], [8, 4], gas(323.15, 101_325.0));
    let a = base.load().terms()[0];
    let b = dense.load().terms()[0];
    let c = warm.load().terms()[0];
    assert!((b.resistance_pa_s_m3 / a.resistance_pa_s_m3 - 2.0).abs() < 1e-12);
    assert!((b.rate_per_s / a.rate_per_s - 1.0).abs() < 1e-12);
    assert!(
        (c.rate_per_s / a.rate_per_s - warm.exterior().sound_speed / base.exterior().sound_speed)
            .abs()
            < 1e-12
    );
    assert_ne!(
        base.impedance_at(800.0, true).unwrap(),
        warm.impedance_at(800.0, true).unwrap()
    );
}

#[test]
fn shared_memory_work_is_storage_plus_escaped_energy_including_release() {
    let g = gas(293.15, 101_325.0);
    let mut load = MouthRadiation::new(
        model([0.01, 0.004], [8, 4], g),
        g.density * g.sound_speed / (core::f64::consts::PI * 0.006 * 0.006),
    )
    .unwrap();
    let (mut work, mut escaped, mut returned_work) = (0.0, 0.0, false);
    for i in 0..512 {
        let q = if i < 64 {
            1e-5
        } else if i < 80 {
            1e-6
        } else {
            0.0
        };
        let (direct, history) = load.affine_pressure();
        let old_flow = load.state.branch_flows_m3_s()[0];
        let old_energy = load.state.stored_energy_j();
        let (next, frame) = load.candidate(q).unwrap();
        let scale = frame.port.supplied_work_j.abs()
            + frame.port.storage_change_j.abs()
            + frame.port.dissipated_energy_j;
        // The balance is exact in real arithmetic, but storage_change is a
        // difference of two O(E) energies and the port pressure is 2a - Z*U
        // after an O(eps*U) flow solve. Neither rounding shrinks with this
        // step's work, so a converged step is roundoff of E and of the
        // characteristic-wave power (Z + direct)*U^2*dt; each enters at most
        // ~7 eps of this floor, and 64 eps matches candidate()'s flow gate.
        let floor = 64.0
            * f64::EPSILON
            * (old_energy
                + frame.port.stored_energy_j
                + (load.port_impedance + direct)
                    * q.abs().max(old_flow.abs()).powi(2)
                    * load.model.time_step_s);
        let tolerance = 1e-12 * scale.max(1e-30) + floor;
        assert!(frame.port.balance_residual_j().abs() <= tolerance);
        assert!(frame.port.dissipated_energy_j >= 0.0);
        let expected_work = (direct * q + history) * q * load.model.time_step_s;
        assert!((frame.port.supplied_work_j - expected_work).abs() <= tolerance);
        returned_work |= frame.port.supplied_work_j < 0.0;
        work += frame.port.supplied_work_j;
        escaped += frame.port.dissipated_energy_j;
        load = next;
    }
    assert!(
        returned_work,
        "near-field energy must return through the port, not become fake loss"
    );
    assert!(escaped > 0.0);
    assert!((work - escaped - load.state.stored_energy_j()).abs() < 1e-12 * escaped);
}

fn voice() -> FlueVoice {
    let g = gas(293.15, 101_325.0);
    let card = fs_aeroac::jetcard::mint_tonal_interim_card().unwrap();
    let geometry = FlueGeometry {
        flue_width_m: 0.01,
        flue_height_m: 0.001,
        cut_up_m: 0.004,
        labium_offset_m: 0.0,
        jet_angle_rad: 0.0,
        bore_area_m2: core::f64::consts::PI * 0.006 * 0.006,
        mouth_loss_rel: 0.0,
    };
    let source = FlueNetworkSpec {
        graph: TubeNetworkSpec {
            nodes: vec![
                NetworkNode::Inlet,
                NetworkNode::Termination { reflection: -0.9 },
            ],
            sections: vec![TubeSection {
                nodes: [0, 1],
                length_m: 0.3,
                radius_m: 0.006,
                max_length_error_m: 0.005,
            }],
            sound_speed_m_s: g.sound_speed,
            max_wave_memory_bytes: 1_000_000,
        },
        gases: vec![g],
        viscothermal: vec![],
        minimum_jet_speed_m_s: 1.0,
        maximum_jet_mach: 0.3,
        max_jet_history_samples: 4096,
    };
    FlueVoice::from_network(
        Some(&card),
        90.0,
        1e-6,
        geometry,
        source,
        96_000,
        3,
        0,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn same_sample_jet_feedback_sees_the_retained_radiation_pressure() {
    let mut v = voice();
    v.attach_mouth_radiation(model([0.01, 0.004], [8, 4], *v.gas()), [0.01, 0.004])
        .unwrap();
    assert_eq!((v.inertance, v.resistance), (0.0, 0.0));
    v.mouth_radiation = Some(v.mouth_radiation.unwrap().candidate(1e-5).unwrap().0);
    let old = v.mouth_radiation.unwrap().state.stored_energy_j();
    let coefficients = JetStepCoefficients {
        q_max: 2e-4,
        half_width_m: 0.0005,
        offset_m: 0.0,
        band_alpha: 0.1,
        hp_alpha: 0.7,
        leak_rate: 10.0,
        eta_override: None,
    };
    let delayed = RetardedVelocity {
        past_m_s: 0.1,
        current_weight: 0.4,
    };
    let boundary = crate::wind_line::WindBoundary {
        direct_reflection: 0.25,
        history_pressure_pa: 0.02,
    };
    let (step, trial, loss) = v
        .solve_jet_mouth(boundary, delayed, 1e-4, coefficients, None)
        .unwrap();
    let (_, memory) = v
        .mouth_radiation
        .unwrap()
        .candidate(step.flow_m3_s)
        .unwrap();
    let residual = trial.source_pa - step.pressure_pa - memory.port.pressure_pa;
    assert!(residual.abs() < 1e-11 * trial.source_pa.abs());
    assert_eq!(loss, 0.0);
    assert_eq!(
        v.mouth_radiation.unwrap().state.stored_energy_j(),
        old,
        "root trials must not advance memory"
    );
    let accepted_jet = v
        .jet_trial(
            delayed.past_m_s + delayed.current_weight * step.flow_m3_s / v.geometry.mouth_area_m2(),
            1e-4,
            coefficients,
        )
        .unwrap();
    assert!((trial.eta_m - accepted_jet.eta_m).abs() < 1e-12 * trial.eta_m.abs());
    v.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    v.step_block(&mut [0.0; 512]).unwrap();
    for d in v.diagnostics() {
        let scale = d.source_work_j.abs()
            + d.pipe_work_j.abs()
            + d.stored_delta_j.abs()
            + d.radiation_loss_j;
        assert!(d.ledger_defect_j.abs() < 1e-10 * scale);
        assert!(d.radiation_loss_j >= 0.0);
        assert_eq!(
            d.numerical_dissipation_j, 0.0,
            "the replaced backward-Euler inertia must not remain"
        );
    }
}

#[test]
fn invalid_bands_sources_and_double_loading_refuse() {
    let g = gas(293.15, 101_325.0);
    for (size, cells, dt, band) in [
        ([0.1, 0.1], [4, 4], 1.0 / 96_000.0, 1000.0),
        ([0.01, 0.004], [4, 4], 1.0 / 1000.0, 1000.0),
        ([0.01, 0.004], [33, 33], 1.0 / 96_000.0, 1000.0),
        ([0.01, 0.004], [0, 4], 1.0 / 96_000.0, 1000.0),
    ] {
        assert!(
            RectangularMouthLoad::new(size, cells, g, dt, band, &CancelGate::new_clock_free())
                .is_err()
        );
    }
    for load in [
        model([0.01, 0.004], [8, 4], gas(293.15, 202_650.0)),
        model([0.008, 0.005], [8, 4], g),
    ] {
        assert!(voice().attach_mouth_radiation(load, [0.01, 0.004]).is_err());
    }
    let load = model([0.01, 0.004], [8, 4], g);
    let mut v = voice();
    v.geometry.mouth_loss_rel = 0.01;
    assert!(v.attach_mouth_radiation(load, [0.01, 0.004]).is_err());
    let mut v = voice();
    v.attach_mouth_radiation(load, [0.01, 0.004]).unwrap();
    assert!(v.attach_mouth_radiation(load, [0.01, 0.004]).is_err());
}
