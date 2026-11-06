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
