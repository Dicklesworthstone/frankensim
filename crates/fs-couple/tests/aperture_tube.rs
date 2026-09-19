//! Coupled physical consumer checks. All constitutive data here are synthetic.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{
    ApertureDrive, ApertureState, ApertureTerminal, DynamicAperture, DynamicApertureSpec,
};
use fs_couple::bernoulli_aperture::tube::{ApertureTube, TubeDrive, TubeFrame, UniformTubeSpec};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

const DT: f64 = 1e-5;
fn tube_spec(delay: usize, reflection: f64) -> UniformTubeSpec {
    UniformTubeSpec {
        length_m: delay as f64 * (343.0 * DT), radius_m: 0.007,
        sound_speed_m_s: 343.0, terminal_reflection: reflection,
        max_length_error_m: 0.0, max_wave_memory_bytes: 1 << 20,
    }
}
fn aperture(spec: UniformTubeSpec, budget: u64, state: ApertureState) -> DynamicAperture {
    DynamicAperture::new(DynamicApertureSpec {
        aperture: BernoulliAperture {
            rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0,
        },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35,
        density_kg_m3: 1.2, impedance_pa_s_m3: spec.characteristic_impedance(1.2).unwrap(),
        time_step_s: DT, max_steps: budget,
    }, state, Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic tube/aperture fixture".into()).unwrap().with_internal_loss(5.0).unwrap()).unwrap()
}
fn model(delay: usize, reflection: f64, budget: u64) -> ApertureTube {
    let spec = tube_spec(delay, reflection);
    ApertureTube::new(aperture(spec, budget, ApertureState {
        opening_m: 4e-4, opening_velocity_m_s: 0.0,
    }), spec).unwrap()
}
fn frame_bits(f: TubeFrame) -> Vec<u64> {
    let a = f.aperture;
    let w = f.waveguide;
    let mut bits = vec![a.step];
    bits.extend([
        a.time_s, a.state.opening_m, a.state.opening_velocity_m_s,
        a.midpoint_opening_m, a.outgoing_pressure_pa, a.bore_pressure_pa,
        a.jet_flow_m3_s, a.swept_flow_m3_s, a.bore_flow_m3_s, a.flow_residual_m3_s,
        a.stored_energy_j, a.storage_change_j, a.dissipated_energy_j, a.pressure_work_j,
        w.incoming_pressure_pa, w.inlet_pressure_pa, w.inlet_flow_m3_s,
        w.terminal_pressure_pa, w.terminal_flow_m3_s, w.stored_energy_j,
        w.storage_change_j, w.inlet_work_j, w.terminal_loss_j,
        f.stored_energy_j, f.storage_change_j, f.dissipated_energy_j,
        f.upstream_work_j, f.body_work_j,
    ].map(f64::to_bits));
    bits
}

#[test]
fn geometry_sets_finite_transit_and_reflected_feedback_not_an_authored_tone() {
    let render = |delay| {
        let mut m = model(delay, -0.8, 128);
        (0..128).map(|_| m.step(TubeDrive {
            upstream_pressure_pa: 800.0, body_flow_m3_s: 0.0,
        }).unwrap()).collect::<Vec<_>>()
    };
    let short = render(8);
    let long = render(12);
    assert!(short[..8].iter().all(|f| f.waveguide.terminal_pressure_pa == 0.0));
    assert!(short[8].waveguide.terminal_pressure_pa.abs() > 1.0);
    assert!(short[..16].iter().all(|f| f.waveguide.incoming_pressure_pa == 0.0));
    assert!(short[16].waveguide.incoming_pressure_pa.abs() > 1.0);
    for i in 0..8 { assert_eq!(frame_bits(short[i]), frame_bits(long[i])); }
    assert!(short[16..].iter().zip(&long[16..])
        .any(|(a,b)| (a.aperture.bore_pressure_pa - b.aperture.bore_pressure_pa).abs() > 1.0));
    assert!(short[16..].iter().zip(&long[16..])
        .any(|(a,b)| (a.aperture.state.opening_m - b.aperture.state.opening_m).abs() > 1e-9));
}

#[test]
fn total_energy_closes_with_separate_upstream_and_body_supplies() {
    let mut missing_body_defect = 0.0_f64;
    for reflection in [-1.0, -0.8, 0.0, 0.5, 1.0] {
        for initial in [
            ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 },
            ApertureState { opening_m: -1e-4, opening_velocity_m_s: -0.2 },
        ] {
            let spec = tube_spec(8, reflection);
            let mut m = ApertureTube::new(aperture(spec, 512, initial), spec).unwrap();
            for i in 0..512 {
                let before = m.stored_energy_j();
                let drive = TubeDrive {
                    upstream_pressure_pa: if i < 128 { 1200.0 } else { 0.0 },
                    body_flow_m3_s: if i < 256 { 2e-7 * (f64::from(i) * 0.11).sin() } else { 0.0 },
                };
                let f = m.step(drive).unwrap();
                let scale = before + f.stored_energy_j + f.dissipated_energy_j
                    + f.upstream_work_j.abs() + f.body_work_j.abs();
                let tolerance = 3e-10 * scale.max(f64::MIN_POSITIVE);
                assert!(f.balance_residual_j().abs() <= tolerance);
                assert!((f.stored_energy_j - before - f.storage_change_j).abs() <= tolerance);
                assert_eq!(m.stored_energy_j().to_bits(), f.stored_energy_j.to_bits());
                assert_eq!(f.aperture.bore_pressure_pa.to_bits(), f.waveguide.inlet_pressure_pa.to_bits());
                assert_eq!(f.aperture.bore_flow_m3_s.to_bits(), f.waveguide.inlet_flow_m3_s.to_bits());
                assert!(f.dissipated_energy_j >= 0.0);
                // Direct source work, not the sum of the two reported residuals.
                let body_work = f.aperture.bore_pressure_pa * (drive.body_flow_m3_s * DT);
                assert_eq!(body_work.to_bits(), f.body_work_j.to_bits());
                missing_body_defect = missing_body_defect.max(
                    (f.storage_change_j + f.dissipated_energy_j - f.upstream_work_j).abs());
            }
        }
    }
    assert!(missing_body_defect > 1e-12, "omitting body-source work must fail the physical balance");
}

#[test]
fn unforced_release_dissipates_the_declared_initial_energy() {
    let spec = tube_spec(8, -0.8);
    let mut m = ApertureTube::new(aperture(spec, 800, ApertureState {
        opening_m: 4.4e-4, opening_velocity_m_s: 0.2,
    }), spec).unwrap();
    let initial = m.stored_energy_j();
    let mut energy = initial;
    for _ in 0..800 {
        let f = m.step(TubeDrive::default()).unwrap();
        assert!(f.stored_energy_j <= energy + 3e-10 * initial);
        assert_eq!(f.upstream_work_j, 0.0);
        assert_eq!(f.body_work_j, 0.0);
        energy = f.stored_energy_j;
    }
    assert!(energy >= 0.0 && energy < 0.1 * initial);
}

#[test]
fn cancelled_and_budgeted_resume_retains_the_whole_wave_and_mechanical_state() {
    let inputs: Vec<_> = (0..128).map(|i| TubeDrive {
        upstream_pressure_pa: if i < 50 { 800.0 } else { 0.0 },
        body_flow_m3_s: if i % 3 == 0 { -1e-7 } else { 2e-7 },
    }).collect();
    let gate = CancelGate::new_clock_free();
    let mut reference = model(8, -0.8, 128);
    let expected: Vec<_> = inputs.iter().map(|d| reference.step(*d).unwrap()).collect();
    let mut resumed = model(8, -0.8, 12);
    let sentinel = TubeFrame { stored_energy_j: -1.0, ..TubeFrame::default() };
    let mut actual = vec![sentinel; 128];
    resumed.advance_block(&inputs[..8], &mut actual[..8], &gate).unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let pause = resumed.advance_block(&inputs[8..], &mut actual[8..], &cancelled).unwrap();
    assert_eq!(pause.completed, 0);
    assert_eq!(pause.terminal, ApertureTerminal::Cancelled);
    assert!(actual[8..].iter().all(|f| *f == sentinel));
    let exhaustion = resumed.advance_block(&inputs[8..], &mut actual[8..], &gate).unwrap();
    assert_eq!(exhaustion.completed, 4);
    assert_eq!(exhaustion.terminal, ApertureTerminal::BudgetExhausted);
    assert_eq!(resumed.aperture().accepted_steps(), 12);
    assert!(actual[12..].iter().all(|f| *f == sentinel));
    resumed.extend_step_budget(128).unwrap();
    resumed.advance_block(&inputs[12..], &mut actual[12..], &gate).unwrap();
    for (a,b) in actual.into_iter().zip(expected) { assert_eq!(frame_bits(a), frame_bits(b)); }
}

#[test]
fn refused_input_does_not_shift_waves_or_poison_the_next_valid_step() {
    let mut a = model(8, -0.8, 128);
    let mut b = model(8, -0.8, 128);
    let valid = TubeDrive { upstream_pressure_pa: 800.0, body_flow_m3_s: 0.0 };
    for _ in 0..32 { a.step(valid).unwrap(); b.step(valid).unwrap(); }
    let energy = a.stored_energy_j();
    let state = a.aperture().state();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(a.step(TubeDrive { body_flow_m3_s: bad, ..valid }).is_err());
        assert_eq!(a.aperture().state(), state);
        assert_eq!(a.aperture().accepted_steps(), 32);
        assert_eq!(a.stored_energy_j().to_bits(), energy.to_bits());
    }
    for _ in 0..32 { assert_eq!(frame_bits(a.step(valid).unwrap()), frame_bits(b.step(valid).unwrap())); }
}

#[test]
fn geometry_approximation_memory_and_clock_are_admitted_explicitly() {
    let state = ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 };
    let mut spec = tube_spec(8, -0.8);
    spec.length_m = 8.25 * (343.0 * DT);
    assert!(ApertureTube::new(aperture(spec, 4, state), spec).is_err());
    spec.max_length_error_m = 0.3 * (343.0 * DT);
    let m = ApertureTube::new(aperture(spec, 4, state), spec).unwrap();
    assert_eq!(m.one_way_samples(), 8);
    assert_eq!(m.spec().length_m.to_bits(), spec.length_m.to_bits());
    assert_ne!(m.represented_length_m().to_bits(), spec.length_m.to_bits());
    assert!((m.represented_length_m() - spec.length_m).abs() <= spec.max_length_error_m);
    spec.max_wave_memory_bytes = 1;
    assert!(ApertureTube::new(aperture(spec, 4, state), spec).is_err());
    let spec = tube_spec(8, -0.8);
    let a = aperture(spec, 4, state);
    let mut incompatible = spec;
    incompatible.radius_m *= 2.0;
    assert!(ApertureTube::new(a, incompatible).is_err());
    let mut advanced = aperture(spec, 4, state);
    advanced.step(ApertureDrive::default()).unwrap();
    assert!(ApertureTube::new(advanced, spec).is_err());
}
