use super::*;
use crate::{Segment, segment_matrix};
use fs_material::gas::GasSpec;
use fs_mms::{LadderSide, OrderGate, RefinementLadder};

const TAU: f64 = core::f64::consts::TAU;
fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn closed() -> C64 {
    C64::new(f64::INFINITY, 0.0)
}

// Independent exponential-horn solution: set y=exp(m*x/2)*p in
// Webster's equation. Then y''+(k^2-m^2/4)y=0. This real form also
// treats imaginary beta and the removable singularity at cutoff.
fn exponential_exact(gas: &GasState, r0: f64, m: f64, length: f64, omega: f64) -> [C64; 4] {
    let k = omega / gas.sound_speed;
    let b2 = k * k - m * m / 4.0;
    let (c, s) = if b2 > 0.0 {
        let b = det::sqrt(b2);
        (det::cos(b * length), det::sin(b * length) / b)
    } else if b2 < 0.0 {
        let b = det::sqrt(-b2);
        let e = det::exp(b * length);
        ((e + 1.0 / e) / 2.0, (e - 1.0 / e) / (2.0 * b))
    } else {
        (1.0, length)
    };
    let area = core::f64::consts::PI * r0 * r0;
    let e = det::exp(m * length / 2.0);
    [
        C64::from_re(e * (c - m * s / 2.0)),
        C64::new(0.0, -omega * gas.density / area / e * s),
        C64::new(0.0, -e * area * k * k / (omega * gas.density) * s),
        C64::from_re((c + m * s / 2.0) / e),
    ]
}
fn load_matrix(m: [C64; 4], load: C64) -> C64 {
    if load.re == f64::INFINITY {
        m[0] / m[2]
    } else {
        (m[0] * load + m[1]) / (m[2] * load + m[3])
    }
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}
fn matrix_error(a: [C64; 4], b: [C64; 4], zref: f64) -> f64 {
    let weights = [1.0, 1.0 / zref, zref, 1.0];
    let numerator = (0..4)
        .map(|i| (a[i] - b[i]).abs() * weights[i])
        .fold(0.0, f64::max);
    let denominator = (0..4).map(|i| b[i].abs() * weights[i]).fold(0.0, f64::max);
    numerator / denominator
}
fn order(case: &str, ns: &[usize], errors: Vec<f64>, p: f64) {
    let ladder =
        RefinementLadder::new(ns.iter().map(|&n| 1.0 / n as f64).collect(), errors).unwrap();
    let verdict = OrderGate { theoretical: p }
        .check(case, LadderSide::Primal, &ladder)
        .unwrap();
    eprintln!("{}", verdict.json_line(true));
}

#[test]
fn riccati_exponential_horn_matches_closed_form() {
    let start = std::time::Instant::now();
    let gas = air();
    let profile = RadiusProfile::exponential(0.005, 2.0, 0.5).unwrap();
    for i in 0..100 {
        // The declared band includes both sides of c*m/(4*pi).
        let hz = 5.0 + 345.0 * f64::from(i) / 99.0;
        let omega = TAU * hz;
        let exact = exponential_exact(&gas, 0.005, 2.0, 0.5, omega);
        // Coordinator ruling2858: exact closed-form two-ports also
        // satisfy1e-12; the impedance sweep retains its separate1e-9 gate.
        let transfer = profile_matrix(
            &profile,
            &gas,
            omega,
            LossModel::Lossless,
            RiccatiSteps(4096),
        )
        .unwrap();
        let zref = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.005 * 0.005);
        let transfer_error = matrix_error(transfer.matrix, exact, zref);
        assert!(
            transfer_error < 1e-12,
            "exponential2-port hz={hz}: {transfer_error:e}"
        );
        for load in [C64::ZERO, closed()] {
            let got = profile_impedance(
                &profile,
                &gas,
                omega,
                LossModel::Lossless,
                load,
                RiccatiSteps(4096),
            )
            .unwrap();
            let want = load_matrix(exact, load);
            let error = relative(got.z_in, want);
            assert!(error < 1e-9, "hz={hz} load={load:?}: rel={error:e}");
            assert_eq!(got.steps_used, 4096);
            eprintln!(
                "{{\"case\":\"riccati_exponential_horn_matches_closed_form\",\"inputs\":\"analytic S=S0 exp(2x), r0=5mm L=.5m\",\"hz\":{hz},\"closed\":{},\"observed\":{error:e},\"oracle\":\"transformed constant-coefficient Webster equation\",\"tol\":1e-9,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
                load.re == f64::INFINITY,
                start.elapsed().as_micros()
            );
        }
    }
    // Cutoff itself is a removable beta=0 limit, not a refusal.
    let omega = gas.sound_speed;
    let got = profile_impedance(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        closed(),
        RiccatiSteps(4096),
    )
    .unwrap();
    assert!(
        relative(
            got.z_in,
            load_matrix(exponential_exact(&gas, 0.005, 2.0, 0.5, omega), closed())
        ) < 1e-9
    );
}

fn power_basis(gas: &GasState, r0: f64, x0: f64, x: f64, omega: f64) -> [C64; 4] {
    let xx = x0 + x;
    let k = omega / gas.sound_speed;
    let z = k * xx;
    let (s, c) = (det::sin(z), det::cos(z));
    let j0 = s / z;
    let y0 = -c / z;
    let j1 = s / (z * z) - c / z;
    let y1 = -c / (z * z) - s / z;
    let radius = r0 * (xx / x0) * (xx / x0);
    let flow = core::f64::consts::PI * radius * radius / (omega * gas.density);
    [
        C64::from_re(j1 / xx),
        C64::from_re(y1 / xx),
        C64::new(0.0, -flow * (k * j0 - 3.0 * j1 / xx) / xx),
        C64::new(0.0, -flow * (k * y0 - 3.0 * y1 / xx) / xx),
    ]
}
#[test]
fn riccati_power_horn_matches_spherical_bessel_form() {
    let start = std::time::Instant::now();
    let gas = air();
    let (r0, x0, length) = (0.005, 0.2, 0.3);
    let profile = RadiusProfile::power(r0, x0, 2.0, length).unwrap();
    for i in 0..100 {
        let omega = TAU * (20.0 + 430.0 * f64::from(i) / 99.0);
        let e0 = power_basis(&gas, r0, x0, 0.0, omega);
        let el = power_basis(&gas, r0, x0, length, omega);
        let det = el[0] * el[3] - el[1] * el[2];
        let inv = [el[3] / det, -el[1] / det, -el[2] / det, el[0] / det];
        let exact = mul2(e0, inv);
        for load in [C64::ZERO, closed()] {
            let got = profile_impedance(
                &profile,
                &gas,
                omega,
                LossModel::Lossless,
                load,
                RiccatiSteps(4096),
            )
            .unwrap();
            let error = relative(got.z_in, load_matrix(exact, load));
            assert!(error < 1e-9, "omega={omega}: rel={error:e}");
            eprintln!(
                "{{\"case\":\"riccati_power_horn_matches_spherical_bessel_form\",\"inputs\":\"analytic r=.005((.2+x)/.2)^2 L=.3m\",\"omega\":{omega},\"observed\":{error:e},\"oracle\":\"j1(kX)/X,y1(kX)/X elementary solutions\",\"tol\":1e-9,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
                start.elapsed().as_micros()
            );
        }
    }
}

#[test]
fn riccati_cylinder_and_cone_reduce_to_tmm() {
    let start = std::time::Instant::now();
    let gas = air();
    let omega = TAU * 270.0;
    let r0 = 0.005;
    let length = 0.3;
    let zref = gas.density * gas.sound_speed / (core::f64::consts::PI * r0 * r0);
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        let profile = RadiusProfile::power(r0, 0.3, 0.0, length).unwrap();
        let wave = segment_wave(&gas, r0, omega, loss).unwrap();
        let exact = segment_matrix(&Segment::Cylinder { radius: r0, length }, &wave).unwrap();
        let got = profile_matrix(&profile, &gas, omega, loss, RiccatiSteps(4096)).unwrap();
        assert!(matrix_error(got.matrix, exact, zref) < 1e-12);
    }
    // Only a LOSSLESS whole cone has constant specific line parameters.
    // A finite lossy cone freezes them at one radius; its continuous
    // counterpart is verified by the separate slice refinement test.
    let profile = RadiusProfile::power(r0, 0.3, 1.0, length).unwrap();
    let wave = segment_wave(&gas, 0.0075, omega, LossModel::Lossless).unwrap();
    let exact = segment_matrix(
        &Segment::Cone {
            inlet_radius: r0,
            outlet_radius: 0.01,
            length,
        },
        &wave,
    )
    .unwrap();
    let got = profile_matrix(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        RiccatiSteps(4096),
    )
    .unwrap();
    assert!(matrix_error(got.matrix, exact, zref) < 1e-12);
    // A live falsifier for substituting a single frozen Bessel cone
    // for continuous propagation. These are different variable-coefficient
    // models; the slice ladder, rather than equality, links their limits.
    let wave = segment_wave(&gas, 0.0075, omega, LossModel::Bessel).unwrap();
    let frozen = segment_matrix(
        &Segment::Cone {
            inlet_radius: r0,
            outlet_radius: 0.01,
            length,
        },
        &wave,
    )
    .unwrap();
    let continuous =
        profile_matrix(&profile, &gas, omega, LossModel::Bessel, RiccatiSteps(4096)).unwrap();
    let discrepancy = matrix_error(continuous.matrix, frozen, zref);
    assert!(
        discrepancy > 1e-8,
        "frozen-radius model falsifier was not live"
    );
    eprintln!(
        "{{\"case\":\"finite_loss_cone_model_form\",\"inputs\":\"analytic 5..10mm cone L.3m at270Hz\",\"observed\":{discrepancy:e},\"oracle\":\"continuous local ZK vs single frozen-radius spherical cone\",\"tol\":\">1e-8 falsifier liveness\",\"verdict\":\"pass\",\"elapsed_us\":{}}}",
        start.elapsed().as_micros()
    );
}

#[test]
fn riccati_step_ladder_observed_order() {
    let start = std::time::Instant::now();
    let gas = air();
    let profile = RadiusProfile::exponential(0.005, 2.0, 0.5).unwrap();
    let omega = TAU * 1500.0;
    let exact = load_matrix(exponential_exact(&gas, 0.005, 2.0, 0.5, omega), closed());
    let ns = [25, 50, 100, 200];
    let mut errors = Vec::new();
    for n in ns {
        let got = profile_impedance(
            &profile,
            &gas,
            omega,
            LossModel::Lossless,
            closed(),
            RiccatiSteps(n),
        )
        .unwrap();
        let error = relative(got.z_in, exact);
        errors.push(error);
        eprintln!(
            "{{\"case\":\"riccati_step_ladder_observed_order\",\"inputs\":\"analytic exponential horn m2 L.5 at1500Hz\",\"axis\":\"Magnus steps\",\"rung\":{n},\"observed\":{error:e},\"oracle\":\"closed-form impedance\",\"tol\":\"order4+/-.2\",\"verdict\":\"pass\",\"elapsed_us\":{}}}",
            start.elapsed().as_micros()
        );
    }
    assert!(errors[0] > 1e-6, "liveness: {errors:?}");
    order("riccati_step_ladder_observed_order", &ns, errors, 4.0);
}

#[test]
fn spline_profile_converges_to_analytic() {
    let start = std::time::Instant::now();
    let gas = air();
    let omega = TAU * 200.0;
    let exact = load_matrix(exponential_exact(&gas, 0.005, 2.0, 0.5, omega), closed());
    let ns = [8, 16, 32, 64];
    let mut errors = Vec::new();
    for spans in ns {
        let stations: Vec<_> = (0..=spans)
            .map(|i| {
                let x = 0.5 * i as f64 / spans as f64;
                (x, 0.005 * det::exp(x))
            })
            .collect();
        let endpoints = [0.005, 0.005 * det::exp(0.5)];
        let profile = RadiusProfile::spline(stations.clone(), Some(endpoints)).unwrap();
        assert_eq!(profile.stations(), Some(stations.as_slice()));
        for &(x, r) in &stations {
            assert!((profile.radius(x).unwrap() - r).abs() < 1e-17);
        }
        let got = profile_impedance(
            &profile,
            &gas,
            omega,
            LossModel::Lossless,
            closed(),
            RiccatiSteps(8192),
        )
        .unwrap();
        let error = relative(got.z_in, exact);
        errors.push(error);
        eprintln!(
            "{{\"case\":\"spline_profile_converges_to_analytic\",\"inputs\":\"analytic exp station values and exact endpoint derivatives; clamped C2 cubic\",\"axis\":\"station spans\",\"rung\":{spans},\"observed\":{error:e},\"oracle\":\"analytic exponential impedance\",\"tol\":\"order4+/-.2\",\"verdict\":\"pass\",\"elapsed_us\":{}}}",
            start.elapsed().as_micros()
        );
    }
    order("spline_profile_converges_to_analytic", &ns, errors, 4.0);
    let natural =
        RadiusProfile::spline(vec![(1.0, 0.005), (1.2, 0.006), (1.4, 0.007)], None).unwrap();
    assert!((natural.radius(0.1).unwrap() - 0.0055).abs() < 1e-17);
}

#[test]
fn cone_slice_ladder_converges() {
    let start = std::time::Instant::now();
    let gas = air();
    let omega = TAU * 110.0;
    let (rin, rout, length) = (0.004, 0.020, 0.6);
    let profile = RadiusProfile::power(rin, length * rin / (rout - rin), 1.0, length).unwrap();
    let exact =
        profile_matrix(&profile, &gas, omega, LossModel::Bessel, RiccatiSteps(8192)).unwrap();
    let zref = gas.density * gas.sound_speed / (core::f64::consts::PI * rin * rin);
    let ns = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];
    let mut errors = Vec::new();
    let mut values = Vec::new();
    for n in ns {
        let t = cone_matrix(
            rin,
            rout,
            length,
            &gas,
            omega,
            LossModel::Bessel,
            ConeSlices::Fixed(n),
        )
        .unwrap();
        assert_eq!(t.steps_used, n);
        let error = matrix_error(t.matrix, exact.matrix, zref);
        errors.push(error);
        values.push(t.matrix);
        eprintln!(
            "{{\"case\":\"cone_slice_ladder_converges\",\"inputs\":\"analytic linear cone4..20mm L.6m air293.15K at110Hz; shared exact ZK\",\"axis\":\"spherical slices\",\"rung\":{n},\"observed\":{error:e},\"oracle\":\"8192-step continuous Magnus plane wave\",\"tol\":\"order2+/-.2; finest pair<1e-8\",\"verdict\":\"pass\",\"elapsed_us\":{}}}",
            start.elapsed().as_micros()
        );
    }
    assert_ne!(
        values[3], values[4],
        "32->64 must change the evaluated port"
    );
    assert!(
        errors[errors.len() - 2..].iter().all(|&e| e < 1e-8),
        "{errors:?}"
    );
    // Ruling2858: the literal32/64 pair cannot satisfy1e-8 for a finite
    // lossy cone. Preserve every initial rung and refine through4096;
    // the actual finest pair must meet that SAME band. Use the declared64+
    // asymptotic range for the shared strict order gate.
    order(
        "cone_slice_ladder_converges",
        &ns[4..],
        errors[4..].to_vec(),
        2.0,
    );
    let auto = ConeSlices::Auto.count(rin, rout).unwrap();
    assert!(auto > 12);
}

#[test]
fn cone_slice_midpoint_uses_shared_zk_bitwise() {
    // Ruling2858 replaces frozen-cone equality by lossless closed-form
    // equality, the lossy ladder, AND this bitwise material-owner check.
    let gas = air();
    let (rin, rout, length, omega, n) = (0.004, 0.020, 0.6, TAU * 110.0, 64_usize);
    let mut expected = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    for i in 0..n {
        let t0 = i as f64 / n as f64;
        let t1 = (i + 1) as f64 / n as f64;
        let ra = rin + (rout - rin) * t0;
        let rb = rin + (rout - rin) * t1;
        // Observe the material seam the production cone loop calls.
        let actual = crate::cone_midpoint_wave(&gas, ra, rb, omega, LossModel::Bessel).unwrap();
        let owner = segment_wave(&gas, f64::midpoint(ra, rb), omega, LossModel::Bessel).unwrap();
        for (a, b) in [
            (actual.wavenumber, owner.wavenumber),
            (
                actual.characteristic_impedance,
                owner.characteristic_impedance,
            ),
            (actual.specific_impedance, owner.specific_impedance),
            (
                C64::from_re(actual.shear_number),
                C64::from_re(owner.shear_number),
            ),
        ] {
            assert_eq!(a.re.to_bits(), b.re.to_bits(), "slice{i}");
            assert_eq!(a.im.to_bits(), b.im.to_bits(), "slice{i}");
        }
        let s = segment_matrix(
            &Segment::Cone {
                inlet_radius: ra,
                outlet_radius: rb,
                length: length / n as f64,
            },
            &owner,
        )
        .unwrap();
        expected = mul2(expected, s);
    }
    let actual = cone_matrix(
        rin,
        rout,
        length,
        &gas,
        omega,
        LossModel::Bessel,
        ConeSlices::Fixed(n),
    )
    .unwrap();
    for (a, b) in actual.matrix.into_iter().zip(expected) {
        assert_eq!(a.re.to_bits(), b.re.to_bits());
        assert_eq!(a.im.to_bits(), b.im.to_bits());
    }
}

#[test]
fn profile_refusals() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.005, 0.0, 0.1).unwrap();
    assert!(
        profile_matrix(
            &profile,
            &gas,
            TAU * 100.0,
            LossModel::Lossless,
            RiccatiSteps(MAX_RICCATI_STEPS)
        )
        .is_ok()
    );
    for n in [0, MAX_RICCATI_STEPS + 1] {
        let e = profile_matrix(
            &profile,
            &gas,
            TAU * 100.0,
            LossModel::Lossless,
            RiccatiSteps(n),
        )
        .unwrap_err();
        assert!(e.to_string().starts_with("FS-DUCT-"));
    }
    assert!(
        cone_matrix(
            0.005,
            0.005,
            0.1,
            &gas,
            TAU * 100.0,
            LossModel::Lossless,
            ConeSlices::Fixed(MAX_CONE_SLICES)
        )
        .is_ok()
    );
    for n in [0, MAX_CONE_SLICES + 1] {
        assert!(ConeSlices::Fixed(n).count(0.005, 0.01).is_err());
    }
    assert!(RadiusProfile::spline(vec![(0.0, 0.005); MAX_PROFILE_STATIONS + 1], None).is_err());
    let cap_stations = (0..MAX_PROFILE_STATIONS)
        .map(|i| (i as f64, 0.005))
        .collect();
    assert!(RadiusProfile::spline(cap_stations, None).is_ok());
    for rows in [
        vec![(1.0, 0.005), (0.0, 0.006)],
        vec![(0.0, 0.0), (1.0, 0.005)],
        vec![(0.0, 0.005), (f64::NAN, 0.006)],
    ] {
        assert!(
            RadiusProfile::spline(rows, None)
                .unwrap_err()
                .to_string()
                .starts_with("FS-DUCT-")
        );
    }
    assert!(RadiusProfile::spline(vec![(0.0, 0.005), (1.0, 0.005)], Some([-1.0, 1.0])).is_err());
    assert!(RadiusProfile::exponential(0.0, 1.0, 1.0).is_err());
    assert!(RadiusProfile::exponential(1.0, f64::INFINITY, 1.0).is_err());
    assert!(RadiusProfile::power(0.005, 0.0, 2.0, 1.0).is_err());
    assert!(profile.radius(-f64::MIN_POSITIVE).is_err());
    assert!(profile.radius(f64::NAN).is_err());
    assert!(profile_matrix(&profile, &gas, 0.0, LossModel::Lossless, RiccatiSteps(10)).is_err());
    assert!(
        profile_impedance(
            &profile,
            &gas,
            TAU * 100.0,
            LossModel::Bessel,
            C64::new(f64::NAN, 0.0),
            RiccatiSteps(10)
        )
        .is_err()
    );
    let mut polls = 0;
    let out = profile_matrix_controlled(
        &profile,
        &gas,
        TAU * 100.0,
        LossModel::Bessel,
        RiccatiSteps(20),
        &mut || {
            polls += 1;
            polls == 4
        },
    );
    assert!(matches!(
        out,
        Err(DuctError::Cancelled { completed_steps: 3 })
    ));
    assert_eq!(polls, 4);
}

#[test]
fn profile_bitwise_repeat() {
    let gas = air();
    let profile = RadiusProfile::exponential(0.004, 3.0, 0.3).unwrap();
    let a = profile_impedance(
        &profile,
        &gas,
        TAU * 200.0,
        LossModel::Bessel,
        C64::ZERO,
        RiccatiSteps(512),
    )
    .unwrap();
    let b = profile_impedance(
        &profile,
        &gas,
        TAU * 200.0,
        LossModel::Bessel,
        C64::ZERO,
        RiccatiSteps(512),
    )
    .unwrap();
    for (x, y) in [
        (a.z_in, b.z_in),
        (a.p_out_over_u_in, b.p_out_over_u_in),
        (a.u_out_over_u_in, b.u_out_over_u_in),
    ] {
        assert_eq!(x.re.to_bits(), y.re.to_bits());
        assert_eq!(x.im.to_bits(), y.im.to_bits());
    }
    assert_eq!(a.min_shear_number.to_bits(), b.min_shear_number.to_bits());
}
