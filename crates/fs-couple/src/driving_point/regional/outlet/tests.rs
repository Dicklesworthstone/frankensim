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

#[test]
fn a_matched_delay_becomes_the_exact_delayed_physical_difference_not_a_wrapped_echo() {
    let rate = 48000;
    let delay = 40;
    let zc = 2.0e6;
    let time = delay as f64 / f64::from(rate);
    let (kernel, actual_delay) = acceleration_kernel(rate, 256, time, |omega| {
        // Independent matched-line solution: Zin=Zc; Qin=p+/Zc;
        // Qout(t)=Qin(t-tau). No use of the implementation's transfer formula.
        Ok(C64::new((omega*time).cos(), (omega*time).sin()).scale(1.0 / zc))
    }).unwrap();
    assert_eq!(actual_delay, delay);
    assert_eq!(kernel.len(), 256 + delay);
    for (n, value) in kernel.iter().enumerate() {
        let expected = if n == delay { f64::from(rate)/zc }
            else if n == delay+1 { -f64::from(rate)/zc } else { 0.0 };
        assert!((value-expected).abs() < 1e-12, "tap {n}: {value} versus {expected}");
        if n < delay { assert_eq!(*value, 0.0); }
    }
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
    let mut shunt_no_vortex = shunt_1.without_vortex_loss();

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

        // Gate 3: small-amplitude limit matches the vortex-free model to rounding.
        let p_soft = 1.0e-3 * phase.sin();
        let u_soft = shunt_soft.step(p_soft).unwrap();
        let u_lin = shunt_no_vortex.step(p_soft).unwrap();
        assert_eq!(shunt_no_vortex.s12_history(), 0.0);
        max_soft_u = max_soft_u.max(u_lin.abs());
        max_soft_diff = max_soft_diff.max((u_soft - u_lin).abs());
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
    assert!(
        max_soft_diff / max_soft_u < 1.0e-9,
        "small-amplitude orifice flow diverged from vortex-free model: rel={}",
        max_soft_diff / max_soft_u
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
