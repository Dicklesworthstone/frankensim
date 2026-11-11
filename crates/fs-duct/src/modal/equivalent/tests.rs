use super::*;
use crate::{LossModel, input_impedance_load};
use fs_material::gas::GasSpec;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

// A Bessel horn r(s) = r_in (u / (u - s))^flare through r_out at s = length, as cones on the curve.
fn bessel_cones(r_in: f64, r_out: f64, length: f64, flare: f64, pieces: usize) -> Vec<Segment> {
    let q = det::pow(r_out / r_in, 1.0 / flare);
    let u = q * length / (q - 1.0);
    let step = length / pieces as f64;
    let mut previous = r_in;
    (1..=pieces)
        .map(|k| {
            let radius = if k == pieces {
                r_out
            } else {
                r_in * det::pow(u / (u - step * k as f64), flare)
            };
            let cone = Segment::Cone {
                inlet_radius: previous,
                outlet_radius: radius,
                length: step,
            };
            previous = radius;
            cone
        })
        .collect()
}

// The tenor-trombone slide, neck, tuning slide, sleeve and four Bessel bell horns of the Braden (2006) template
// the brass cards use (without the mouthpiece and leadpipe): two radius steps and a 64-cone bell to 108.6 mm.
fn trombone_like() -> Duct {
    let mut segments = vec![
        Segment::Cylinder {
            radius: 0.00695,
            length: 1.39439,
        },
        Segment::Cylinder {
            radius: 0.00725,
            length: 0.26910,
        },
        Segment::Cone {
            inlet_radius: 0.00755,
            outlet_radius: 0.00990,
            length: 0.2649,
        },
        Segment::Cylinder {
            radius: 0.00990,
            length: 0.06732,
        },
        Segment::Cone {
            inlet_radius: 0.01053,
            outlet_radius: 0.01105,
            length: 0.05336,
        },
    ];
    for (r_in, r_out, length, flare) in [
        (0.01105, 0.01507, 0.15042, 1.260),
        (0.01507, 0.02228, 0.15042, 0.894),
        (0.02228, 0.04118, 0.09685, 0.494),
        (0.04118, 0.10860, 0.05593, 1.110),
    ] {
        segments.extend(bessel_cones(r_in, r_out, length, flare, 16));
    }
    Duct { segments }
}

// A cylinder with a Bessel bell to 40 mm: a mouth whose flanged radiation limit (ka = 24, 33 kHz) is reachable with a short kernel step.
fn bell_bore() -> Duct {
    let mut segments = vec![Segment::Cylinder {
        radius: 0.0065,
        length: 0.5,
    }];
    segments.extend(bessel_cones(0.0065, 0.04, 0.15, 1.0, 16));
    Duct { segments }
}

fn z_inf(gas: &GasState, radius: f64) -> f64 {
    gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius)
}

#[test]
fn with_one_mode_the_equivalent_load_is_the_baffled_piston_at_the_actual_mouth() {
    let gas = air();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0065,
                length: 0.5,
            },
            Segment::Cone {
                inlet_radius: 0.0065,
                outlet_radius: 0.05,
                length: 0.2,
            },
        ],
    };
    let mouth = PreparedModalDuct::new(&duct, 1, 1)
        .unwrap()
        .mouth_radius_m();
    for f in [50.0, 200.0, 800.0, 2000.0, 5000.0] {
        let omega = core::f64::consts::TAU * f;
        let z = equivalent_plane_mouth_impedance(&duct, &gas, &gas, omega, LossModel::Bessel, 1, 1)
            .unwrap();
        let (r, x) =
            fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, mouth, omega).unwrap();
        let exact = C64::new(r, x);
        let error = (z - exact).abs() / exact.abs();
        assert!(error <= 1e-10, "{f} Hz: relative error {error:e}");
    }
}

#[test]
fn the_equivalent_load_of_the_trombone_like_bore_has_no_negative_resistance() {
    // Bound fixed before the first run: Re Z_eq >= -1e-3 Z_inf from 40 to 3000 Hz (the smallest value measured was +3.1e-3 Z_inf at 40 Hz).
    let gas = air();
    let duct = trombone_like();
    let mouth = PreparedModalDuct::new(&duct, 3, 1)
        .unwrap()
        .mouth_radius_m();
    let zi = z_inf(&gas, mouth);
    let mut worst = (f64::INFINITY, 0.0);
    let mut f = 40.0;
    while f <= 3000.0 {
        let z = equivalent_plane_mouth_impedance(
            &duct,
            &gas,
            &gas,
            core::f64::consts::TAU * f,
            LossModel::Bessel,
            3,
            1,
        )
        .unwrap();
        assert!(z.re.is_finite() && z.im.is_finite(), "{f} Hz");
        if z.re / zi < worst.0 {
            worst = (z.re / zi, f);
        }
        f += 20.0;
    }
    assert!(
        worst.0 >= -1e-3,
        "Re Z_eq / Z_inf = {:e} at {} Hz",
        worst.0,
        worst.1
    );
}

#[test]
fn on_a_uniform_cylinder_the_plane_transfer_matrix_with_the_equivalent_load_is_the_multimodal_bore()
{
    // Inside a uniform tube the plane mode decouples from the higher modes, so the plane transfer matrix
    // terminated with the plane-mode equivalent load is exactly the multimodal input impedance.
    let gas = air();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.05,
            length: 0.4,
        }],
    };
    for f in [100.0, 450.0, 777.0, 1500.0, 2500.0, 3000.0] {
        let omega = core::f64::consts::TAU * f;
        let modal = mm_input_impedance_flanged(&duct, &gas, &gas, omega, LossModel::Bessel, 3, 1)
            .unwrap()
            .plane_impedance;
        let z = equivalent_plane_mouth_impedance(&duct, &gas, &gas, omega, LossModel::Bessel, 3, 1)
            .unwrap();
        let plane = input_impedance_load(&duct, &gas, omega, LossModel::Bessel, z, None)
            .unwrap()
            .impedance;
        let error = (plane - modal).abs() / modal.abs();
        assert!(error <= 1e-6, "{f} Hz: relative error {error:e}");
    }
}

#[test]
fn the_causal_kernel_reproduces_the_equivalent_load_between_its_bins() {
    // Bounds fixed before the first run: |Z_kernel - Z_eq| <= 1e-3 Z_inf at every off-grid midpoint up to 0.9 of the claimed band (4 kHz, bins 10 Hz),
    // the weights sum to one once the dropped tail's sum is added back (1e-9), and the dropped energy is within the declared tolerance (1e-10).
    // Measured: 3.5e-5 at 3.5 kHz, sum -3.4e-10, dropped energy 9.2e-11.
    const OFF_GRID_BOUND: f64 = 1e-3;
    const SUM_BOUND: f64 = 1e-9;
    let gas = air();
    let duct = bell_bore();
    let spec = MouthKernelSpec {
        step_s: 1.25e-5,
        bin_hz: 10.0,
        band_hz: 4000.0,
        tail_tolerance: 1e-10,
    };
    let kernel = plane_mouth_kernel(&duct, &gas, &gas, LossModel::Bessel, 3, 1, spec).unwrap();
    let sum: f64 = kernel.weights.iter().sum::<f64>() + kernel.dropped_sum;
    assert!((sum - 1.0).abs() <= SUM_BOUND, "sum - 1 = {:e}", sum - 1.0);
    assert!(kernel.dropped_tail_energy <= spec.tail_tolerance);
    assert!(
        kernel.evaluated_hz
            >= 0.9 * 24.0 * gas.sound_speed / (core::f64::consts::TAU * kernel.mouth_radius_m)
    );
    let top = (spec.band_hz / kernel.bin_hz).floor() as usize;
    let mut worst = (0.0_f64, 0.0_f64);
    for k in (1..=(0.9 * top as f64) as usize).step_by(9) {
        let f = (k as f64 + 0.5) * kernel.bin_hz;
        let omega = core::f64::consts::TAU * f;
        let exact =
            equivalent_plane_mouth_impedance(&duct, &gas, &gas, omega, LossModel::Bessel, 3, 1)
                .unwrap();
        let error = (kernel.impedance_at(omega) - exact).abs() / kernel.z_inf;
        if error > worst.0 {
            worst = (error, f);
        }
    }
    assert!(
        worst.0 <= OFF_GRID_BOUND,
        "worst {:e} at {} Hz",
        worst.0,
        worst.1
    );
}

#[test]
fn a_kernel_step_too_coarse_to_reach_the_radiation_limit_or_a_claim_above_a_third_of_the_evaluation_is_refused()
 {
    let gas = air();
    let duct = bell_bore();
    let spec = |step_s: f64, band_hz: f64| MouthKernelSpec {
        step_s,
        bin_hz: 10.0,
        band_hz,
        tail_tolerance: 1e-10,
    };
    // 25 us has a Nyquist frequency of 20 kHz: the evaluation reaches 19.6 kHz of the 33 kHz limit (ka 14 of 24), where the omitted resistance
    // still leaks 1.4e-3 Z_inf into the reactance at 3.5 kHz (measured before the refusal existed), so that step is refused.
    assert!(
        plane_mouth_kernel(
            &duct,
            &gas,
            &gas,
            LossModel::Bessel,
            3,
            1,
            spec(2.5e-5, 4000.0)
        )
        .is_err()
    );
    // A claim of 12 kHz is above a third of the 33 kHz evaluation.
    assert!(
        plane_mouth_kernel(
            &duct,
            &gas,
            &gas,
            LossModel::Bessel,
            3,
            1,
            spec(1.25e-5, 12_000.0)
        )
        .is_err()
    );
    assert!(
        plane_mouth_kernel(
            &duct,
            &gas,
            &gas,
            LossModel::Bessel,
            3,
            1,
            spec(1.25e-5, 0.0)
        )
        .is_err()
    );
    assert!(
        plane_mouth_kernel(
            &duct,
            &gas,
            &gas,
            LossModel::Bessel,
            0,
            1,
            spec(1.25e-5, 4000.0)
        )
        .is_err()
    );
}

// Rule 6 (MB 24783, B2): equal gases realize, bit for bit, the kernel the builder realized before its reference impedance became the BORE gas's.
// The digests below were read on the unchanged builder (tip e667e944f7, hz3) from the same two inputs, then the builder was changed.
fn kernel_digest(kernel: &MouthKernel) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    let mut eat = |x: f64| {
        for byte in x.to_bits().to_le_bytes() {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    for x in [kernel.step_s, kernel.z_inf, kernel.mouth_radius_m, kernel.dropped_tail_energy, kernel.dropped_sum, kernel.bin_hz, kernel.evaluated_hz] {
        eat(x);
    }
    eat(kernel.weights.len() as f64);
    for &w in &kernel.weights {
        eat(w);
    }
    hash
}

const PRE_CHANGE_DIGEST_BELL_BESSEL: u64 = 0xdc68_0437_c392_7380; // 181 taps
const PRE_CHANGE_DIGEST_TROMBONE_LOSSLESS: u64 = 0xe87b_7179_13be_4757; // 525 taps

#[test]
fn equal_gases_realize_the_kernel_they_always_did() {
    let gas = air();
    let spec = MouthKernelSpec { step_s: 1.25e-5, bin_hz: 10.0, band_hz: 4000.0, tail_tolerance: 1e-10 };
    let bell = plane_mouth_kernel(&bell_bore(), &gas, &gas, LossModel::Bessel, 3, 1, spec).unwrap();
    let trombone = plane_mouth_kernel(&trombone_like(), &gas, &gas, LossModel::Lossless, 3, 1, spec).unwrap();
    let (a, b) = (kernel_digest(&bell), kernel_digest(&trombone));
    println!("DIGESTS bell/bessel {a:#018x} trombone/lossless {b:#018x} (taps {} and {})", bell.weights.len(), trombone.weights.len());
    assert_eq!(a, PRE_CHANGE_DIGEST_BELL_BESSEL, "bell/bessel digest {a:#018x}");
    assert_eq!(b, PRE_CHANGE_DIGEST_TROMBONE_LOSSLESS, "trombone/lossless digest {b:#018x}");
}

// B2 (MB 24783): two gases. The bore holds the bore gas, the exterior is the room; the load is the two-gas flanged law and the kernel's reference impedance is the
// bore gas's. 301 K dry air against 293.15 K dry air has the `rho c` ratio of the 27 C saturated breath state against the room (0.987, a 1.3 % mismatch).
fn warm_bore_gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 301.0, 101_325.0).unwrap()
}

#[test]
fn a_kernels_reference_impedance_is_the_bore_gas_and_its_radiation_limit_the_exterior() {
    let (bore, room) = (warm_bore_gas(), air());
    let spec = MouthKernelSpec { step_s: 1.25e-5, bin_hz: 10.0, band_hz: 4000.0, tail_tolerance: 1e-10 };
    let kernel = plane_mouth_kernel(&bell_bore(), &bore, &room, LossModel::Lossless, 3, 1, spec).unwrap();
    assert_eq!(kernel.z_inf, z_inf(&bore, kernel.mouth_radius_m));
    let ratio = kernel.z_inf / z_inf(&room, kernel.mouth_radius_m);
    println!("B2: rho c bore / room {ratio:.5}, reflection {:.5}", (1.0 - ratio) / (1.0 + ratio));
    assert!((ratio - 1.0).abs() > 5e-3, "the two gases must differ in rho c: {ratio}");
    assert!(kernel.evaluated_hz >= 0.9 * 24.0 * room.sound_speed / (core::f64::consts::TAU * kernel.mouth_radius_m));
}

#[test]
fn on_a_uniform_cylinder_the_two_gas_equivalent_load_closes_the_bore_gas_to_the_two_gas_multimodal_bore() {
    // Bound fixed before the first run: the one-gas test's 1e-6. The plane transfer matrix of the bore gas, closed by the equivalent load of the bore gas against the
    // exterior, is the multimodal bore of mm_input_impedance_flanged(duct, bore gas, exterior).
    let (bore, room) = (warm_bore_gas(), air());
    let duct = Duct { segments: vec![Segment::Cylinder { radius: 0.05, length: 0.4 }] };
    for f in [100.0, 450.0, 777.0, 1500.0, 2500.0, 3000.0] {
        let omega = core::f64::consts::TAU * f;
        let modal = mm_input_impedance_flanged(&duct, &bore, &room, omega, LossModel::Bessel, 3, 1).unwrap().plane_impedance;
        let z = equivalent_plane_mouth_impedance(&duct, &bore, &room, omega, LossModel::Bessel, 3, 1).unwrap();
        let plane = input_impedance_load(&duct, &bore, omega, LossModel::Bessel, z, None).unwrap().impedance;
        let error = (plane - modal).abs() / modal.abs();
        assert!(error <= 1e-6, "{f} Hz: relative error {error:e}");
    }
}

#[test]
fn at_low_ka_the_two_gas_equivalent_load_of_a_wide_straight_pipe_is_the_two_gas_flanged_pipe_law() {
    // Bound fixed before the first run: 1e-6. A uniform pipe of radius 0.05 m and length 0.4 m has no mode conversion, and its first higher mode (cutoff ka 3.83) decays
    // by exp(-76 x 0.4) ~ 6e-14 over the pipe at these frequencies, so the plane-mode load with the higher modes eliminated is `flanged_pipe_impedance_with_exterior`
    // (the exterior owns the radiation, the bore gas the eliminated pipe modes). The one-gas law at the exterior is the control that the gases matter.
    // Measured (hz3, run duct1): the relative difference is 0 to 3.5e-16 at 60, 150, 400 and 900 Hz, and the one-gas law at the exterior differs from the two-gas law by 3.5e-4 to 5.2e-4.
    let (bore, room) = (warm_bore_gas(), air());
    let radius = 0.05;
    let duct = Duct { segments: vec![Segment::Cylinder { radius, length: 0.4 }] };
    for f in [60.0, 150.0, 400.0, 900.0] {
        let omega = core::f64::consts::TAU * f;
        let z = equivalent_plane_mouth_impedance(&duct, &bore, &room, omega, LossModel::Lossless, 3, 1).unwrap();
        let law = crate::flanged_pipe_impedance_with_exterior(&bore, &room, radius, omega, 3).unwrap();
        let one_gas = crate::flanged_pipe_impedance(&room, radius, omega, 3).unwrap();
        let error = (z - law).abs() / law.abs();
        let gases_matter = (law - one_gas).abs() / law.abs();
        println!("B2: {f} Hz two-gas law: relative difference {error:.3e}; the one-gas law at the exterior differs by {gases_matter:.3e}");
        assert!(error <= 1e-6, "{f} Hz: relative error {error:e}");
        assert!(gases_matter > 1e-4, "{f} Hz: the bore gas does not reach the law");
    }
}

#[test]
fn the_causal_kernel_of_two_gases_reproduces_the_two_gas_load_between_its_bins() {
    // Bound fixed before the first run: 2e-3 Z_inf at every off-grid midpoint up to 0.9 of the claimed band (the one-gas test's bound is 1e-3, measured 3.5e-5). The
    // extra is the resistance the realized load omits above the evaluated top, where it is the bore gas's Z_inf and the radiation resistance tends to the exterior's: a
    // step of delta = 1.3 % leaks (delta / pi) ln((f_top + f) / (f_top - f)) into the reactance, about 9e-4 at 3.6 kHz of a 32.7 kHz top.
    // Measured (hz3, tip e667e944f7 export, debug build): 3.98e-4 Z_inf at 3525 Hz, weights sum - 1 = -3.6e-10, dropped energy 9.6e-11.
    const OFF_GRID_BOUND: f64 = 2e-3;
    // Bessel loss as in the one-gas test (the Lossless kernel of this bore sums to one within 7.8e-9, which is not what is asked here).
    let (bore, room) = (warm_bore_gas(), air());
    let duct = bell_bore();
    let spec = MouthKernelSpec { step_s: 1.25e-5, bin_hz: 10.0, band_hz: 4000.0, tail_tolerance: 1e-10 };
    let kernel = plane_mouth_kernel(&duct, &bore, &room, LossModel::Bessel, 3, 1, spec).unwrap();
    let sum: f64 = kernel.weights.iter().sum::<f64>() + kernel.dropped_sum;
    println!("B2: weights sum - 1 = {:e}, dropped tail energy {:e}", sum - 1.0, kernel.dropped_tail_energy);
    assert!((sum - 1.0).abs() <= 1e-9, "sum - 1 = {:e}", sum - 1.0);
    assert!(kernel.dropped_tail_energy <= spec.tail_tolerance);
    let top = (spec.band_hz / kernel.bin_hz).floor() as usize;
    let mut worst = (0.0_f64, 0.0_f64);
    for k in (1..=(0.9 * top as f64) as usize).step_by(9) {
        let f = (k as f64 + 0.5) * kernel.bin_hz;
        let omega = core::f64::consts::TAU * f;
        let exact = equivalent_plane_mouth_impedance(&duct, &bore, &room, omega, LossModel::Bessel, 3, 1).unwrap();
        let error = (kernel.impedance_at(omega) - exact).abs() / kernel.z_inf;
        if error > worst.0 {
            worst = (error, f);
        }
    }
    println!("B2: worst off-grid error {:.3e} Z_inf at {} Hz (bound {OFF_GRID_BOUND:e})", worst.0, worst.1);
    assert!(worst.0 <= OFF_GRID_BOUND, "worst {:e} at {} Hz", worst.0, worst.1);
}
