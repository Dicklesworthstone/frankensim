use crate::driving_point::C64;
use super::*;
use fs_duct::{Duct, Segment};
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_qty::semantic::{Composition, CompositionBasis};

fn gas(co2: bool) -> GasState {
    let (species, fractions) = if co2 {
        (vec![GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::CarbonDioxide], vec![0.4, 0.2, 0.4])
    } else { (vec![GasSpecies::Nitrogen, GasSpecies::Oxygen], vec![0.7, 0.3]) };
    resolve_gas_mixture(&species,
        &Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(), 293.15, 100000.0).unwrap()
}
fn location() -> CircularOutletReceiver {
    CircularOutletReceiver { position_m: [0.0, 0.0, 1.0], radial_rings: 2,
        angular_points: 8, maximum_frequency_hz: 2000.0 }
}
fn duct(a: GasState, b: GasState) -> RegionalDuct {
    RegionalDuct::new(&Duct { segments: vec![
        Segment::Cylinder { radius: 0.0015, length: 0.2 },
        Segment::Cone { inlet_radius: 0.0015, outlet_radius: 0.0012, length: 0.3 },
    ] }, &[a, b]).unwrap()
}

fn duct_with_open_and_vented_holes(air: GasState) -> RegionalDuct {
    let hole = |state| Segment::ToneHole {
        hole_radius: 0.0006, chimney_height: 0.003, bore_radius: 0.0015, state,
    };
    RegionalDuct::new(&Duct { segments: vec![
        Segment::Cylinder { radius: 0.0015, length: 0.12 },
        Segment::ToneHole { hole_radius: 0.0005, chimney_height: 0.003,
            bore_radius: 0.0015, state: HoleState::Open },
        Segment::Cylinder { radius: 0.0015, length: 0.17 },
        hole(HoleState::Closed),
        hole(HoleState::Vent(0.0)),
        hole(HoleState::Vent(0.25)),
        Segment::Cylinder { radius: 0.0015, length: 0.21 },
    ] }, &[air; 7]).unwrap()
}

#[test]
fn explicit_common_mouth_speeds_preserve_observers_and_scaled_shunt_bits() {
    let air = gas(false);
    let d = duct_with_open_and_vented_holes(air);
    let mut receiver = location();
    receiver.maximum_frequency_hz = 700.0;
    let (mut common_bell, mut common_holes, common_shunt) = d
        .baffled_outlet_and_open_holes_and_shunt_delayed_at_with_port_scale(
            8000, 256, receiver, WallModel::Rigid, air.sound_speed, 0.5,
        ).unwrap();
    let (mut explicit_bell, mut explicit_holes, explicit_shunt) = d
        .baffled_outlet_and_open_holes_and_shunt_with_delay_speeds(
            8000, 256, receiver, WallModel::Rigid, air.sound_speed,
            &[air.sound_speed; 2], 0.5,
        ).unwrap();
    assert_eq!(common_holes.len(), 2);
    assert_eq!(explicit_holes.len(), 2);
    assert_eq!(common_bell.info().fft_samples, explicit_bell.info().fft_samples);
    let mut common_shunt = common_shunt.unwrap();
    let mut explicit_shunt = explicit_shunt.unwrap();
    assert_eq!(common_shunt.r2_pa_s_m3.to_bits(), explicit_shunt.r2_pa_s_m3.to_bits());
    for sample in 0..1024 {
        let pressure = 0.1 * (0.13 * sample as f64).sin();
        assert_eq!(common_bell.step(pressure).unwrap().to_bits(),
            explicit_bell.step(pressure).unwrap().to_bits());
        for (common, explicit) in common_holes.iter_mut().zip(&mut explicit_holes) {
            assert_eq!(common.step(pressure).unwrap().to_bits(), explicit.step(pressure).unwrap().to_bits());
        }
        assert_eq!(common_shunt.step(pressure).unwrap().to_bits(),
            explicit_shunt.step(pressure).unwrap().to_bits());
        assert_eq!(common_shunt.s12_history().to_bits(), explicit_shunt.s12_history().to_bits());
    }
}

#[test]
fn each_open_or_vented_mouth_uses_its_own_original_rayleigh_speed_and_site() {
    let air = gas(false);
    let d = duct_with_open_and_vented_holes(air);
    let sites = d.open_hole_sites();
    assert_eq!(sites.len(), 2);
    assert!((d.total_axial_length_m() - 0.5).abs() < 1e-14);
    let mut receiver = location();
    receiver.position_m = [0.05, 0.1, 0.8];
    receiver.maximum_frequency_hz = 700.0;
    let bell_speed = air.sound_speed * 1.1;
    let speeds = [air.sound_speed * 0.8, air.sound_speed * 1.3];
    let (bell, holes, _) = d.baffled_outlet_and_open_holes_and_shunt_with_delay_speeds(
        8000, 256, receiver, WallModel::Rigid, bell_speed, &speeds, 1.0,
    ).unwrap();
    let bell_rayleigh = BaffledPressure::circular_outlet(0.0015, 8000, receiver,
        RayleighMedium { density: air.density, sound_speed: bell_speed }).unwrap();
    assert_eq!(bell.info().exterior_delay_samples, bell_rayleigh.delay_samples);
    for (hole, ((axial_position, radius), speed)) in
        holes.iter().zip([(0.12, 0.0005), (0.29, 0.0003)].into_iter().zip(speeds))
    {
        let mut expected_receiver = receiver;
        expected_receiver.position_m[0] += 0.5 - axial_position;
        assert!((hole.info().receiver.position_m[0] - expected_receiver.position_m[0]).abs() < 1e-14);
        assert_eq!(hole.info().radius_m, radius);
        let mut original = BaffledPressure::circular_outlet(radius, 8000, expected_receiver,
            RayleighMedium { density: air.density, sound_speed: speed }).unwrap();
        assert_eq!(hole.info().exterior_delay_samples, original.delay_samples);
        let mut radiation = hole.radiation.clone();
        let mut previous = 0.0;
        let mut energy = 0.0;
        for sample in 0..400 {
            let flow = 1e-5 * (0.17 * sample as f64).sin();
            let expected = original.step(&[flow]).unwrap();
            let actual = radiation.push((flow - previous) * 8000.0).unwrap();
            previous = flow;
            assert!((actual - expected).abs() < 1e-12);
            energy += expected * expected;
        }
        assert!(energy > 0.0);
    }
}

#[test]
fn mouth_speeds_admit_exact_open_count_and_finite_positive_values_before_observation() {
    let air = gas(false);
    let d = duct_with_open_and_vented_holes(air);
    // Deliberately invalid observation clock: speed admission must happen first.
    let refused = |bell_speed, hole_speeds: &[f64]| d
        .baffled_outlet_and_open_holes_and_shunt_with_delay_speeds(
            0, 0, location(), WallModel::Rigid, bell_speed, hole_speeds, 1.0,
        ).unwrap_err();
    for speeds in [&[][..], &[air.sound_speed][..], &[air.sound_speed; 3][..]] {
        assert!(matches!(refused(air.sound_speed, speeds),
            OutletError::Transfer(DrivingPointError::Invalid {
                what: "one propagation speed is required per open or vented hole"
            })));
    }
    for speed in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        for (bell_speed, speeds) in [(speed, [air.sound_speed; 2]),
            (air.sound_speed, [air.sound_speed, speed])]
        {
            assert!(matches!(refused(bell_speed, &speeds),
                OutletError::Transfer(DrivingPointError::Invalid {
                    what: "mouth propagation speeds must be finite and positive"
                })));
        }
    }
    let no_holes = duct(air, air);
    assert!(matches!(no_holes.baffled_outlet_and_open_holes_and_shunt_with_delay_speeds(
        0, 0, location(), WallModel::Rigid, air.sound_speed, &[air.sound_speed], 1.0,
    ).unwrap_err(), OutletError::Transfer(DrivingPointError::Invalid {
        what: "one propagation speed is required per open or vented hole"
    })));
}

#[test]
fn a_matched_delay_becomes_the_exact_delayed_physical_difference_not_a_wrapped_echo() {
    let rate = 48000;
    let delay = 40;
    let zc = 2.0e6;
    let time = delay as f64 / f64::from(rate);
    let (residual, actual_delay) = acceleration_residual(rate, 256, time, |omega| {
        // Independent matched-line solution: Zin=Zc; Qin=p+/Zc;
        // Qout(t)=Qin(t-tau). No use of the implementation's transfer formula.
        Ok(C64::new((omega*time).cos(), (omega*time).sin()).scale(1.0 / zc))
    }).unwrap();
    assert_eq!(actual_delay, delay);
    // A whole-sample arrival has no pre-arrival spread: nothing replays.
    let latency = outlet_pre_arrival_latency(256);
    let band = (100.0, 4800.0);
    let replay = replay_gain(&residual, latency, rate, band).0;
    let replay_without_latency = replay_gain(&residual, 0, rate, band).0;
    // FFT roundoff only: 1e-12 sits 9 decades below the replay gate.
    assert_eq!(latency, 64);
    assert!(
        replay <= 1e-12 && replay_without_latency <= 1e-12,
        "{replay:e} {replay_without_latency:e}"
    );
    let kernel = realize_with_latency(delay, &residual, latency);
    assert_eq!(kernel.len(), 256 + delay);
    // The window's last `latency` taps sit between the exact zeros and the
    // arrival; a whole-sample arrival leaves only FFT roundoff there. Without
    // them the kernel is the exact delayed difference, checked as before.
    let realized = kernel;
    assert!(
        realized[delay..delay + latency]
            .iter()
            .all(|value| value.abs() < 1e-12)
    );
    let kernel: Vec<f64> = realized[..delay]
        .iter()
        .chain(&realized[delay + latency..])
        .copied()
        .collect();
    for (n, value) in kernel.iter().enumerate() {
        let expected = if n == delay { f64::from(rate)/zc }
            else if n == delay+1 { -f64::from(rate)/zc } else { 0.0 };
        assert!((value-expected).abs() < 1e-12, "tap {n}: {value} versus {expected}");
        if n < delay { assert_eq!(*value, 0.0); }
    }
    // The realized kernel: the common observer latency moves the exact
    // difference `latency` later.
    let (kernel, delay) = (realized, delay + latency);
    let mut filter = DelayedFilter::from_impulse_response(1.0/f64::from(rate), kernel).unwrap();
    let mut values = Vec::new();
    for n in 0..700 {
        let input = (0.073*n as f64).sin() + 0.2*(0.19*n as f64).cos();
        values.push(input);
        let current = if n>=delay {values[n-delay]} else {0.0};
        let previous = if n>delay {values[n-delay-1]} else {0.0};
        let expected = (current-previous)*f64::from(rate)/zc;
        assert!((filter.push(input).unwrap()-expected).abs() < 1e-12);
    }
}

#[test]
fn extracted_spatial_filter_reproduces_the_shared_rayleigh_receiver_without_differentiating_twice() {
    let medium = RayleighMedium { density: 1.3, sound_speed: 320.0 };
    let mut receiver = location(); receiver.position_m = [0.2, 0.1, 1.0];
    let mut original = BaffledPressure::circular_outlet(0.01, 48000, receiver, medium).unwrap();
    let mut filter = rayleigh_acceleration_filter(
        BaffledPressure::circular_outlet(0.01, 48000, receiver, medium).unwrap(), 48000).unwrap();
    let mut previous = 0.0;
    for n in 0..1000 {
        let volume_flow = 1e-5*(0.061*n as f64).sin();
        let acceleration = (volume_flow-previous)*48000.0;
        previous = volume_flow;
        let expected = original.step(&[volume_flow]).unwrap();
        let actual = filter.push(acceleration).unwrap();
        assert!((actual-expected).abs() < 1e-12, "sample {n}");
    }
}

#[test]
fn regional_outlet_uses_terminal_geometry_and_gas_with_a_causal_material_dependent_response() {
    let air = gas(false); let enriched = gas(true);
    let baseline = duct(air, air);
    let regional = duct(air, enriched);
    let mut a = baseline.baffled_outlet(48000, 1024, location(), None).unwrap();
    let mut b = regional.baffled_outlet(48000, 1024, location(), None).unwrap();
    assert_eq!(b.info().gas, enriched);
    assert_eq!(b.info().radius_m, 0.0012);
    assert!(b.info().axial_delay_samples > a.info().axial_delay_samples);
    assert!(b.info().exterior_delay_samples.0 > a.info().exterior_delay_samples.0);
    let first = b.info().axial_delay_samples + b.info().exterior_delay_samples.0;
    let mut energy = 0.0; let mut change = 0.0;
    for n in 0..2000 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let pa = a.step(input).unwrap(); let pb = b.step(input).unwrap();
        assert!(pa.is_finite() && pb.is_finite());
        if n < first { assert_eq!(pb, 0.0); }
        energy += pb*pb; change += (pa-pb)*(pa-pb);
    }
    assert!(energy > 1e-12 && change > 1e-4*energy);
}

#[test]
fn construction_refusals_and_invalid_input_preserve_the_complete_healthy_cascade() {
    let d = duct(gas(false), gas(true));
    // The next power of two above the line ladder's ceiling refuses.
    for n in [0, 255, 257, 2 * MAX_OUTLET_FFT_SAMPLES, usize::MAX] {
        assert!(d.baffled_outlet(48000, n, location(), None).is_err());
    }
    assert!(d.baffled_outlet(0, 1024, location(), None).is_err());
    let mut receiver = location(); receiver.position_m[2] = -1.0;
    assert!(d.baffled_outlet(48000, 1024, receiver, None).is_err());
    let mut a = d.baffled_outlet(48000, 1024, location(), None).unwrap();
    let mut b = a.clone();
    for n in 0..400 {
        assert!(a.step(f64::NAN).is_err());
        let input = (n as f64 * 0.3).sin();
        assert_eq!(a.step(input).unwrap(), b.step(input).unwrap());
    }
}

#[test]
fn a_second_stage_overflow_poison_cannot_resume_with_inconsistent_histories() {
    let mut observer = duct(gas(false), gas(false)).baffled_outlet(48000, 256, location(), None).unwrap();
    observer.acceleration = DelayedFilter::from_impulse_response(1.0/48000.0,
        vec![f64::MAX, 0.0, 0.0, 0.0]).unwrap();
    observer.radiation = DelayedFilter::from_impulse_response(1.0/48000.0,
        vec![2.0, 0.0, 0.0, 0.0]).unwrap();
    assert!(observer.step(1.0).is_err());
    let history = observer.acceleration.history();
    assert!(observer.step(0.0).is_err());
    assert_eq!(observer.acceleration.history(), history);
}

#[test]
fn open_hole_shunt_line_is_port_scale_invariant_passive_at_small_amplitude_and_uniform_across_holes() {
    let air = gas(false);
    let register_duct = RegionalDuct::new(
        &Duct {
            segments: vec![
                Segment::Cylinder { radius: 0.00745, length: 0.15 },
                Segment::ToneHole {
                    hole_radius: 0.00135,
                    chimney_height: 0.0127,
                    bore_radius: 0.00745,
                    state: fs_duct::HoleState::Open,
                },
                Segment::Cylinder { radius: 0.00745, length: 0.11 },
            ],
        },
        &[air, air, air],
    )
    .unwrap();
    let (_, _, shunt_1) = register_duct
        .baffled_outlet_and_open_holes_and_shunt_delayed_at_with_port_scale(
            48000,
            1024,
            location(),
            WallModel::Rigid,
            air.sound_speed,
            1.0,
        )
        .unwrap();
    let (_, _, shunt_half) = register_duct
        .baffled_outlet_and_open_holes_and_shunt_delayed_at_with_port_scale(
            48000,
            1024,
            location(),
            WallModel::Rigid,
            air.sound_speed,
            0.5,
        )
        .unwrap();
    let (_, _, shunt_double_r2) = register_duct
        .baffled_outlet_and_open_holes_and_shunt_delayed_at_with_port_scale(
            48000,
            1024,
            location(),
            WallModel::Rigid,
            air.sound_speed,
            2.0,
        )
        .unwrap();
    let mut shunt_1 = shunt_1.expect("open hole should build a WDF orifice shunt line");
    let mut shunt_half = shunt_half.expect("0.5x R2 open hole should build a WDF shunt line");
    let mut shunt_double_r2 =
        shunt_double_r2.expect("2.0x R2 open hole should build a WDF shunt line");
    assert_eq!(shunt_1.u_crit_m3_s, 0.0);
    assert!((shunt_half.r2_pa_s_m3 - 0.5 * shunt_1.r2_pa_s_m3).abs() < 1.0e-6);
    assert!((shunt_double_r2.r2_pa_s_m3 - 2.0 * shunt_1.r2_pa_s_m3).abs() < 1.0e-6);

    // Gate 2: R2 port-resistance scaling (0.5x, 1.0x, 2.0x) must leave the solved
    // orifice flow U_h(t) and reflected reed wave delta p^-(t) invariant to < 1e-9.
    let mut shunt_double_drive = shunt_1.clone();
    let mut shunt_soft = shunt_1.clone();
    let mut shunt_soft_half = shunt_1.clone();
    let mut shunt_no_vortex = shunt_1.without_vortex_loss();
    let mut shunt_no_vortex_half = shunt_1.without_vortex_loss();

    let mut max_dp_1: f64 = 0.0;
    let mut max_dp_diff_half: f64 = 0.0;
    let mut max_dp_diff_double: f64 = 0.0;
    let mut max_u_1: f64 = 0.0;
    let mut max_u_diff_half: f64 = 0.0;
    let mut max_u_diff_double: f64 = 0.0;
    let mut energy_single = 0.0;
    let mut energy_double = 0.0;
    let mut max_soft_u: f64 = 0.0;
    let mut max_soft_diff: f64 = 0.0;
    let mut max_soft_half_u: f64 = 0.0;
    let mut max_soft_half_diff: f64 = 0.0;

    for n in 0..512 {
        let phase = 2.0 * std::f64::consts::PI * 900.0 * n as f64 / 48000.0;
        let p_in = 1500.0 * phase.sin();
        let u1 = shunt_1.step(p_in).unwrap();
        let u_half = shunt_half.step(p_in).unwrap();
        let u_dbl_r2 = shunt_double_r2.step(p_in).unwrap();
        let u_dbl_drive = shunt_double_drive.step(2.0 * p_in).unwrap();
        assert!(u1.is_finite() && u_half.is_finite() && u_dbl_r2.is_finite() && u_dbl_drive.is_finite());

        let dp1 = shunt_1.s12_history();
        let dp_half = shunt_half.s12_history();
        let dp_dbl_r2 = shunt_double_r2.s12_history();
        let dp_dbl_drive = shunt_double_drive.s12_history();

        max_dp_1 = max_dp_1.max(dp1.abs());
        max_dp_diff_half = max_dp_diff_half.max((dp_half - dp1).abs());
        max_dp_diff_double = max_dp_diff_double.max((dp_dbl_r2 - dp1).abs());
        max_u_1 = max_u_1.max(u1.abs());
        max_u_diff_half = max_u_diff_half.max((u_half - u1).abs());
        max_u_diff_double = max_u_diff_double.max((u_dbl_r2 - u1).abs());

        energy_single += dp1 * dp1;
        energy_double += dp_dbl_drive * dp_dbl_drive;

        // Gate 3: small-amplitude quadratic deviation from vortex-free model is < 1e-6 at 1 mPa
        // and scales linearly with amplitude (halving p_soft halves rel within 1%).
        let p_soft = 1.0e-3 * phase.sin();
        let u_soft = shunt_soft.step(p_soft).unwrap();
        let u_lin = shunt_no_vortex.step(p_soft).unwrap();
        assert_eq!(shunt_no_vortex.s12_history(), 0.0);
        max_soft_u = max_soft_u.max(u_lin.abs());
        max_soft_diff = max_soft_diff.max((u_soft - u_lin).abs());

        let p_soft_half = 0.5 * p_soft;
        let u_soft_half = shunt_soft_half.step(p_soft_half).unwrap();
        let u_lin_half = shunt_no_vortex_half.step(p_soft_half).unwrap();
        max_soft_half_u = max_soft_half_u.max(u_lin_half.abs());
        max_soft_half_diff = max_soft_half_diff.max((u_soft_half - u_lin_half).abs());
    }

    assert!(
        max_dp_diff_half / max_dp_1 < 1.0e-9 && max_dp_diff_double / max_dp_1 < 1.0e-9,
        "R2 port scaling changed reflected reed wave: rel_half={}, rel_double={}",
        max_dp_diff_half / max_dp_1,
        max_dp_diff_double / max_dp_1
    );
    assert!(
        max_u_diff_half / max_u_1 < 1.0e-9 && max_u_diff_double / max_u_1 < 1.0e-9,
        "R2 port scaling changed orifice flow: rel_half={}, rel_double={}",
        max_u_diff_half / max_u_1,
        max_u_diff_double / max_u_1
    );
    let rel_soft = max_soft_diff / max_soft_u;
    let rel_soft_half = max_soft_half_diff / max_soft_half_u;
    assert!(
        rel_soft < 1.0e-6 && (rel_soft_half / (0.5 * rel_soft) - 1.0).abs() < 0.01,
        "small-amplitude orifice flow deviation failed linear-scaling gate: rel_1mPa={rel_soft}, rel_0.5mPa={rel_soft_half}"
    );
    assert!(
        energy_double > 4.0 * energy_single,
        "quadratic orifice loss should grow superlinearly with amplitude: single={energy_single}, double={energy_double}"
    );

    // Gate 4: one law for every open hole with K proportional to b^{-4}.
    let tonehole_duct = RegionalDuct::new(
        &Duct {
            segments: vec![
                Segment::Cylinder { radius: 0.00745, length: 0.10 },
                Segment::ToneHole {
                    hole_radius: 0.00235,
                    chimney_height: 0.0069,
                    bore_radius: 0.00745,
                    state: fs_duct::HoleState::Open,
                },
                Segment::Cylinder { radius: 0.00745, length: 0.20 },
            ],
        },
        &[air, air, air],
    )
    .unwrap();
    let (_, _, tonehole_shunt) = tonehole_duct
        .baffled_outlet_and_open_holes_and_shunt_delayed_at(
            48000,
            1024,
            location(),
            WallModel::Rigid,
            air.sound_speed,
        )
        .unwrap();
    let tonehole_shunt =
        tonehole_shunt.expect("2.35 mm open tone hole must obey the same orifice law");
    let expected_k_ratio = (0.00235_f64 / 0.00135_f64).powi(4);
    let actual_k_ratio = shunt_1.k_quad_pa_s2_m6 / tonehole_shunt.k_quad_pa_s2_m6;
    assert!(
        (actual_k_ratio - expected_k_ratio).abs() / expected_k_ratio < 1.0e-12,
        "actual_k_ratio={actual_k_ratio}, expected_k_ratio={expected_k_ratio}"
    );
}

// Replay gate at the receiver (MB 16347): planted mouths on a 1024-point grid at
// 48 kHz (bins of 46.875 Hz), latency 256, so the replayed taps are [512, 768).
const GATE_N: usize = 1024;
const GATE_RATE: u32 = 48_000;
const GATE_BAND: (f64, f64) = (100.0, 20_000.0);

/// `direct` at the arrival (tap 0) and `late` at tap 600, inside the replayed span.
fn planted(direct: f64, late: f64) -> Vec<f64> {
    let mut residual = vec![0.0; GATE_N];
    residual[0] = direct;
    residual[600] = late;
    residual
}

fn db(ratio: f64) -> f64 {
    20.0 * ratio.log10()
}

/// Bin `k` of one tap `value` at sample `at`, summed directly (no FFT).
fn bin_of(value: f64, at: usize, k: usize) -> C64 {
    let phase = -core::f64::consts::TAU * (k * at) as f64 / GATE_N as f64;
    C64::new(value * phase.cos(), value * phase.sin())
}

#[test]
fn receiver_gate_refuses_a_loud_mouths_replay_at_minus_50_db_re_the_summed_direct_sound() {
    let latency = outlet_pre_arrival_latency(GATE_N);
    assert_eq!(latency, 256);
    let unit: &[f64] = &[1.0];
    // Direct sound 1 + 0.01 = 1.01 at every bin; mouth 0 replays -50 dB re that.
    let loud = planted(1.0, 1.01 * 10.0_f64.powf(-50.0 / 20.0));
    let quiet = planted(0.01, 0.0);
    let mouths = [
        ReceiverMouth {
            residual: &loud,
            axial_delay_samples: 0,
            rayleigh: Some(unit),
        },
        ReceiverMouth {
            residual: &quiet,
            axial_delay_samples: 0,
            rayleigh: Some(unit),
        },
    ];
    let worst = receiver_replay(&mouths, latency, GATE_RATE, GATE_BAND).unwrap();
    assert!((db(worst.gain) + 50.0).abs() < 1e-6, "{worst:?}");
    assert!(worst.direct_db.abs() < 1e-6, "{worst:?}");
    assert!((worst.late_db + 50.0).abs() < 1e-6, "{worst:?}");
    assert_eq!(worst.mouth, 0);
    assert!(worst.mouth_share_db.abs() < 1e-6, "{worst:?}");
    let refusal = admit_replay(&mouths, latency, GATE_RATE, GATE_BAND).unwrap_err();
    let text = refusal.to_string();
    assert!(
        matches!(
            &refusal,
            OutletError::Replay {
                mouth: 0,
                fft_samples: GATE_N,
                ..
            }
        ),
        "{text}"
    );
    assert!(
        text.contains("-50.0 dB re the receiver's direct sound")
            && text.contains("mouth 0 (0 = bell) carries 0.0 dB"),
        "{text}"
    );
}

#[test]
fn receiver_gate_refuses_replays_small_against_each_mouth_but_not_against_their_cancelling_sum() {
    let latency = outlet_pre_arrival_latency(GATE_N);
    let unit: &[f64] = &[1.0];
    // Mouth 1 arrives 4 samples later, inverted, at 0.99: the direct paths cancel
    // to 0.01 (-40 dB re mouth 0) at bin 256 = 12 kHz, the band's only bin where
    // 4 k is a multiple of 1024. Each replay is -80 dB re its own mouth's direct
    // path, and the two arrive in phase there: 1.99e-4 / 0.01, -34.0 dB.
    let a = planted(1.0, 1.0e-4);
    let b = planted(-0.99, 0.99e-4);
    let mouths = [
        ReceiverMouth {
            residual: &a,
            axial_delay_samples: 0,
            rayleigh: Some(unit),
        },
        ReceiverMouth {
            residual: &b,
            axial_delay_samples: 4,
            rayleigh: Some(unit),
        },
    ];
    let worst = receiver_replay(&mouths, latency, GATE_RATE, GATE_BAND).unwrap();
    assert_eq!(worst.frequency_hz, 12_000.0, "{worst:?}");
    assert!(
        (db(worst.gain) - db(1.99e-4 / 0.01)).abs() < 1e-6,
        "{worst:?}"
    );
    assert!(matches!(
        admit_replay(&mouths, latency, GATE_RATE, GATE_BAND),
        Err(OutletError::Replay { .. })
    ));
    // Each mouth against its own direct path passes (the previous gate) ...
    for residual in [&a, &b] {
        let (own, _) = replay_gain(residual, latency, GATE_RATE, GATE_BAND);
        assert!((db(own) + 80.0).abs() < 1e-6, "{own:e}");
    }
    // ... and so does the loudest mouth's direct path as the reference: at that
    // bin, summed directly (no FFT), the replayed sum is -74 dB re mouth 0's.
    let k = 256;
    let direct = [bin_of(1.0, 0, k), bin_of(-0.99, 4, k)];
    let late = bin_of(1.0e-4, 600, k) + bin_of(0.99e-4, 604, k);
    let loudest = direct.iter().map(|d| d.abs()).fold(0.0, f64::max);
    let loudest_reference = db(late.abs() / loudest);
    assert!(
        (loudest_reference + 74.0).abs() < 0.1,
        "{loudest_reference:.2} dB"
    );
    // The gate's FFT sum equals the direct sum at that bin.
    let total = (direct[0] + direct[1]).abs();
    assert!(
        (late.abs() / total / worst.gain - 1.0).abs() < 1e-9,
        "{worst:?}"
    );
    // The 4-sample offset may sit in the axial delay or in the Rayleigh taps.
    let shifted: &[f64] = &[0.0, 0.0, 0.0, 0.0, 1.0];
    let moved = [
        mouths[0],
        ReceiverMouth {
            residual: &b,
            axial_delay_samples: 0,
            rayleigh: Some(shifted),
        },
    ];
    assert_eq!(
        receiver_replay(&moved, latency, GATE_RATE, GATE_BAND).unwrap(),
        worst
    );
}

#[test]
fn receiver_gate_passes_a_bell_100_db_under_a_hole_with_its_own_replay_at_minus_45_db() {
    let latency = outlet_pre_arrival_latency(GATE_N);
    let unit: &[f64] = &[1.0];
    let bell = planted(1.0e-5, 1.0e-5 * 10.0_f64.powf(-45.0 / 20.0));
    let hole = planted(1.0, 0.0);
    // Alone, the bell refuses at -45 dB re its own direct path (the previous gate) ...
    let alone = admit_replay(
        &[ReceiverMouth::alone(&bell)],
        latency,
        GATE_RATE,
        GATE_BAND,
    );
    assert!(
        matches!(&alone, Err(OutletError::Replay { gain, .. }) if (db(*gain) + 45.0).abs() < 1e-6),
        "{alone:?}"
    );
    // ... while at the receiver its replay lies 145 dB under the summed direct sound.
    let mouths = [
        ReceiverMouth {
            residual: &bell,
            axial_delay_samples: 0,
            rayleigh: Some(unit),
        },
        ReceiverMouth {
            residual: &hole,
            axial_delay_samples: 0,
            rayleigh: Some(unit),
        },
    ];
    let (gain, _) = admit_replay(&mouths, latency, GATE_RATE, GATE_BAND).unwrap();
    assert!((db(gain) + 145.0).abs() < 1e-3, "{:.4} dB", db(gain));
}
