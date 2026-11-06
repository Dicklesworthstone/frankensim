//! Calibration curves of the Smith & Woodhouse (2000) rosin law,
//! `fs_tribo::rosin_thermal`: the shifted power law of a `FrictionLaw`
//! (Woodhouse & Galluzzo 2025, eq. 13) beside the paper's Eq. (2). Each test
//! prints JSON lines. The rosin and the rig's perspex are Smith & Woodhouse
//! (2000) Table 1; the rig geometry is their section 5.3, as in Galluzzo,
//! Woodhouse & Mansour (2017) Table 2, "steady sliding".

use fs_math::det;
use fs_tribo::rosin_thermal::{
    RosinCalibration, SmithWoodhouse, SmithWoodhouseParams, SteadySlidingCurve, ThermalSolid,
};
use fs_tribo::{
    FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStepInput, TriboError,
};

/// Rosin layer density [kg/m^3], Table 1 (the pack
/// `rosin-smith-woodhouse-jmps2000`).
const ROSIN_DENSITY: f64 = 1010.0;
/// Rosin layer specific heat [J/kg/K], Table 1.
const ROSIN_SPECIFIC_HEAT: f64 = 250.0;
/// Perspex, Table 1: the rod and the wedge of the steady-sliding rig.
const PERSPEX: ThermalSolid = ThermalSolid {
    conductivity_w_per_m_k: 0.2,
    density_kg_per_m3: 1200.0,
    specific_heat_j_per_kg_k: 1460.0,
};

/// Eq. (13) of Woodhouse & Galluzzo (2025), written out independently.
fn eq13(v: f64) -> f64 {
    0.2 * det::pow(v + 0.011, -0.4)
}

fn params(cal: RosinCalibration, step_s: f64, history_capacity: usize) -> SmithWoodhouseParams {
    SmithWoodhouseParams {
        rosin_density_kg_per_m3: ROSIN_DENSITY,
        rosin_specific_heat_j_per_kg_k: ROSIN_SPECIFIC_HEAT,
        steady_fit: cal.curve,
        calibration: cal.rig,
        calibration_speed_max_mps: cal.speed_max_mps,
        contact: cal.rig,
        step_s,
        history_capacity,
    }
}

fn step(normal_n: f64, slip_m: f64, dt_s: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n,
        slip_increment_m: [slip_m, 0.0],
        dt_s,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 1.0e4,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-14,
        },
    }
}

fn speeds(n: u32, lo: f64, hi: f64) -> impl Iterator<Item = f64> {
    (0..=n).map(move |i| {
        if i == n {
            hi
        } else {
            lo * det::exp(f64::from(i) / f64::from(n) * det::ln(hi / lo))
        }
    })
}

#[test]
fn g1_eq13_curve_is_the_friction_law_bit_for_bit() {
    let law = FrictionLaw::rosin_woodhouse_galluzzo_2025_eq13();
    let curve = SteadySlidingCurve::from_friction_law(&law).unwrap();
    let mut checked = 0;
    for v in speeds(240, 1e-6, 2.0) {
        let own = curve.coefficient(v);
        let shared = law.kinetic_coefficient(v).unwrap();
        assert_eq!(own.to_bits(), shared.to_bits(), "v {v}: {own} vs {shared}");
        assert_eq!(curve.coefficient(-v).to_bits(), own.to_bits(), "odd in v");
        checked += 1;
    }
    let zero = curve.coefficient(0.0);
    assert_eq!(zero.to_bits(), law.static_coefficient().unwrap().to_bits());
    // Hyperbolic maps to exponent 1; that law has its own quotient form, so
    // the two agree to rounding: a few roundings in each fraction, so 8 eps.
    let cello = FrictionLaw::Hyperbolic {
        static_mu: 0.8,
        kinetic_mu: 0.4,
        characteristic_speed: 0.2,
    };
    let hyper = SteadySlidingCurve::from_friction_law(&cello).unwrap();
    let mut worst = 0.0_f64;
    for v in speeds(240, 1e-6, 2.0) {
        let shared = cello.kinetic_coefficient(v).unwrap();
        worst = worst.max((hyper.coefficient(v) - shared).abs() / shared);
    }
    println!(
        "{{\"law\":\"rosin-thermal\",\"test\":\"curve_from_friction_law\",\"eq13_bitwise_speeds\":{checked},\"hyperbolic_rel_err\":{worst:e},\"tol\":{:e}}}",
        8.0 * f64::EPSILON
    );
    assert!(
        worst <= 8.0 * f64::EPSILON,
        "hyperbolic as exponent 1: {worst:e}"
    );
}

#[test]
fn g1_eq13_calibration_reproduces_eq13_through_eq20() {
    // Eq. (20) inverted along eq. (13) defines k_y(T). Along the calibration
    // curve the law must return T from an independent Eq. (20), and
    // A k_y(T)/N must return eq. (13). The rig conducts into its perspex.
    let cal = RosinCalibration::woodhouse_galluzzo_2025_eq13(PERSPEX).unwrap();
    assert_eq!(cal.rig.radius_m.to_bits(), 500e-6_f64.to_bits());
    assert_eq!(cal.rig.reference_normal_n.to_bits(), 3.0_f64.to_bits());
    assert_eq!(cal.rig.layer_thickness_m.to_bits(), 10e-6_f64.to_bits());
    assert_eq!(cal.rig.stationary, PERSPEX);
    assert_eq!(cal.rig.moving, PERSPEX);
    let law = SmithWoodhouse::new(params(cal, 1e-4, 1)).unwrap();
    let (radius, layer, load) = (500e-6, 10e-6, 3.0);
    let area = core::f64::consts::PI * radius * radius;
    let k_rig = PERSPEX.conductivity_w_per_m_k;
    let d = k_rig / (PERSPEX.density_kg_per_m3 * PERSPEX.specific_heat_j_per_kg_k);
    let (mut worst_t, mut worst_mu, mut last_t) = (0.0_f64, 0.0_f64, 0.0_f64);
    for v in speeds(60, 1e-4, 0.3) {
        // Eq. (20): mu N v = v a delta rho_r c_r T + T A K_2 sqrt(3v/(8 a D_2)).
        let oracle_t = eq13(v) * load * v
            / (v * radius * layer * ROSIN_DENSITY * ROSIN_SPECIFIC_HEAT
                + area * k_rig * det::sqrt(3.0 * v / (8.0 * radius * d)));
        let t = law.calibration_steady_temperature_k(v).unwrap();
        assert!(
            t > last_t,
            "T_ss increases with speed: {t} after {last_t} at {v}"
        );
        last_t = t;
        let mu = area * law.yield_stress_pa(t).unwrap() / load;
        worst_t = worst_t.max((t - oracle_t).abs() / oracle_t);
        worst_mu = worst_mu.max((mu - eq13(v)).abs() / eq13(v));
    }
    let static_mu = area * law.yield_stress_pa(0.0).unwrap() / load;
    println!(
        "{{\"law\":\"rosin-thermal\",\"test\":\"eq13_calibration\",\"oracle\":\"eq20,eq13\",\"t_rel_err\":{worst_t:e},\"mu_rel_err\":{worst_mu:e},\"tol\":1e-12,\"static_mu\":{static_mu},\"temperature_limit_k\":{},\"ky0_pa\":{},\"ky_at_limit_pa\":{}}}",
        law.temperature_limit_k(),
        law.yield_stress_pa(0.0).unwrap(),
        law.yield_stress_pa(law.temperature_limit_k()).unwrap()
    );
    assert!(worst_t < 1e-12, "Eq. (20) temperature {worst_t:e}");
    assert!(worst_mu < 1e-12, "eq. (13) through k_y(T) {worst_mu:e}");
    assert!(
        (static_mu - 0.2 * det::pow(0.011, -0.4)).abs() < 1e-14,
        "k_y(0) A/N {static_mu}"
    );
}

#[test]
fn g1_eq13_authors_extended_range_is_opt_in() {
    // The default keeps the measured 0.3 m/s: its limit is the one the row
    // printed before this constructor existed (39.06036861910653 K), and k_y
    // just above it refuses. The authors-extended range (Woodhouse & Galluzzo
    // 2025 Fig. 7, to 5 m/s) has the same curve and rig, admits that
    // temperature, follows an independent Eq. (20) and eq. (13) out to 5 m/s,
    // and refuses above its own limit.
    let measured = RosinCalibration::woodhouse_galluzzo_2025_eq13(PERSPEX).unwrap();
    let extended =
        RosinCalibration::woodhouse_galluzzo_2025_eq13_authors_extended(PERSPEX).unwrap();
    assert_eq!(measured.speed_max_mps.to_bits(), 0.3_f64.to_bits());
    assert_eq!(extended.speed_max_mps.to_bits(), 5.0_f64.to_bits());
    assert_eq!(
        (extended.curve, extended.rig),
        (measured.curve, measured.rig)
    );
    let base = SmithWoodhouse::new(params(measured, 1e-4, 1)).unwrap();
    let wide = SmithWoodhouse::new(params(extended, 1e-4, 1)).unwrap();
    assert_eq!(
        base.temperature_limit_k().to_bits(),
        39.060_368_619_106_53_f64.to_bits()
    );
    let refused = |law: &SmithWoodhouse, t: f64| {
        matches!(
            law.yield_stress_pa(t),
            Err(TriboError::OutsideApplicability {
                field: "smith_woodhouse_temperature_rise",
                maximum,
                ..
            }) if maximum.to_bits() == law.temperature_limit_k().to_bits()
        )
    };
    let just_above = 39.1;
    assert!(refused(&base, just_above), "the default refuses 39.1 K");
    let admitted = wide.yield_stress_pa(just_above).unwrap();
    // Below the measured limit both give the same k_y(T); their bisections
    // start from different brackets, so to rounding.
    let mut worst_same = 0.0_f64;
    for t in [0.0, 1.0, 10.0, 30.0, 39.0] {
        let (a, b) = (
            base.yield_stress_pa(t).unwrap(),
            wide.yield_stress_pa(t).unwrap(),
        );
        worst_same = worst_same.max((a - b).abs() / a);
    }
    let (radius, layer, load) = (500e-6, 10e-6, 3.0);
    let area = core::f64::consts::PI * radius * radius;
    let k_rig = PERSPEX.conductivity_w_per_m_k;
    let d = k_rig / (PERSPEX.density_kg_per_m3 * PERSPEX.specific_heat_j_per_kg_k);
    let oracle = |v: f64| {
        eq13(v) * load * v
            / (v * radius * layer * ROSIN_DENSITY * ROSIN_SPECIFIC_HEAT
                + area * k_rig * det::sqrt(3.0 * v / (8.0 * radius * d)))
    };
    let (mut worst_t, mut worst_mu, mut last_t) = (0.0_f64, 0.0_f64, 0.0_f64);
    for v in speeds(60, 0.3, 5.0) {
        let t = wide.calibration_steady_temperature_k(v).unwrap();
        assert!(
            t > last_t,
            "T_ss increases with speed: {t} after {last_t} at {v}"
        );
        last_t = t;
        let mu = area * wide.yield_stress_pa(t).unwrap() / load;
        worst_t = worst_t.max((t - oracle(v)).abs() / oracle(v));
        worst_mu = worst_mu.max((mu - eq13(v)).abs() / eq13(v));
    }
    let limit = wide.temperature_limit_k();
    assert!(
        refused(&wide, limit + 1e-9),
        "the extended range refuses above its limit"
    );
    assert!(base.calibration_steady_temperature_k(0.31).is_err());
    // No speed range on this rig reaches the about 53 K of the authors' Fig. 7:
    // the oracle's T_ss peaks near 6.5 m/s, a few hundredths of a kelvin above
    // the 5 m/s limit.
    let (mut ceiling, mut at) = (0.0_f64, 0.0_f64);
    for v in speeds(4000, 0.3, 100.0) {
        if oracle(v) > ceiling {
            (ceiling, at) = (oracle(v), v);
        }
    }
    println!(
        "{{\"law\":\"rosin-thermal\",\"test\":\"eq13_authors_extended\",\"oracle\":\"eq20,eq13\",\"measured_limit_k\":{},\"extended_limit_k\":{limit},\"k_y_at_39_1_k_pa\":{admitted},\"same_below_rel_err\":{worst_same:e},\"t_rel_err\":{worst_t:e},\"mu_rel_err\":{worst_mu:e},\"tol\":1e-12,\"oracle_ceiling_k\":{ceiling},\"ceiling_speed_mps\":{at}}}",
        base.temperature_limit_k()
    );
    assert!(
        worst_same < 1e-12,
        "k_y(T) below the measured limit {worst_same:e}"
    );
    assert!(worst_t < 1e-12, "Eq. (20) temperature {worst_t:e}");
    assert!(worst_mu < 1e-12, "eq. (13) through k_y(T) {worst_mu:e}");
    assert!((limit - oracle(5.0)).abs() < 1e-12 * limit, "limit {limit}");
    assert!(limit > 45.2 && limit < 45.3, "limit {limit}");
    assert!(
        ceiling >= limit && ceiling < 45.4 && at > 5.0 && at < 8.0,
        "ceiling {ceiling} at {at}"
    );
}

#[test]
fn g1_eq13_steady_sliding_approaches_eq13() {
    // A port held at constant speed on the calibration rig returns eq. (13)
    // in the limit. The stationary half-space keeps drawing K T/sqrt(pi D t),
    // so mu(t) = mu_13(v) + c t^{-1/2} + O(t^{-1}); the oracles are those of
    // the Eq. (2) test in benchmark_laws.rs.
    let cal = RosinCalibration::woodhouse_galluzzo_2025_eq13(PERSPEX).unwrap();
    let (radius, load) = (500e-6, 3.0);
    for v in [0.05, 0.1, 0.2] {
        let transit = 8.0 * radius / (3.0 * core::f64::consts::PI * v);
        let h = transit / 16.0;
        let marks = [0.25, 1.0, 4.0];
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let ends: Vec<usize> = marks.iter().map(|t| (t / h).round() as usize).collect();
        let law = SmithWoodhouse::new(params(cal, h, ends[2])).unwrap();
        let mut s = law.initial_state();
        let mut mu = Vec::with_capacity(ends[2]);
        for _ in 0..ends[2] {
            let next = law.advance(&s, &step(load, v * h, h)).unwrap();
            mu.push(next.traction_n[0] / load);
            s = next.state;
        }
        let oracle = eq13(v);
        #[allow(clippy::cast_precision_loss)]
        let times: Vec<f64> = ends.iter().map(|&n| n as f64 * h).collect();
        let dev: Vec<f64> = ends.iter().map(|&n| mu[n - 1] - oracle).collect();
        let (r1, r2) = (det::sqrt(times[1]), det::sqrt(times[2]));
        let extrapolated = (r2 * mu[ends[2] - 1] - r1 * mu[ends[1] - 1]) / (r2 - r1);
        let ratio = dev[2] / dev[1];
        let next_order = det::sqrt(transit / times[1]);
        let residual = (extrapolated - oracle).abs();
        println!(
            "{{\"law\":\"rosin-thermal\",\"test\":\"eq13_steady_sliding\",\"v_mps\":{v},\"h_s\":{h:e},\"deviation\":{dev:?},\"ratio_4s_1s\":{ratio},\"ratio_band\":{next_order},\"richardson_residual\":{residual:e},\"residual_bound\":{},\"temperature_rise_k\":{}}}",
            0.1 * dev[2].abs(),
            s.temperature_rise_k()
        );
        assert!(
            dev[0] > dev[1] && dev[1] > dev[2] && dev[2] > 0.0,
            "approach from above: {dev:?}"
        );
        assert!((ratio - 0.5).abs() <= next_order, "t^-1/2 law: {ratio}");
        assert!(
            residual <= 0.1 * dev[2].abs(),
            "Richardson limit {extrapolated} vs eq. (13) {oracle}"
        );
    }
}

#[test]
fn g0_calibration_curves_refuse_by_name() {
    let unsupported = |law: FrictionLaw| {
        matches!(
            SteadySlidingCurve::from_friction_law(&law),
            Err(TriboError::InvalidInput {
                field: "smith_woodhouse_calibration_curve_unsupported_law"
            })
        )
    };
    assert!(unsupported(FrictionLaw::Coulomb {
        static_mu: 0.5,
        kinetic_mu: 0.4
    }));
    assert!(unsupported(FrictionLaw::Stribeck {
        static_mu: 0.8,
        kinetic_mu: 0.4,
        characteristic_speed: 0.1,
        viscous_per_speed: 0.0,
    }));
    assert!(unsupported(FrictionLaw::VelocityDependent {
        static_mu: 0.5,
        mu_zero: 0.5,
        slope_per_speed: -0.1,
    }));
    // An invalid law refuses through its own validation.
    assert!(
        SteadySlidingCurve::from_friction_law(&FrictionLaw::ShiftedPower {
            static_mu: 0.5,
            kinetic_mu: 0.6,
            characteristic_speed: 0.01,
            exponent: 0.4,
        })
        .is_err()
    );
    // A curve that decays faster than |v|^{-1/2} makes T_ss fall at speed:
    // the certificate refuses it.
    let mut cal = RosinCalibration::woodhouse_galluzzo_2025_eq13(PERSPEX).unwrap();
    cal.curve = SteadySlidingCurve::from_friction_law(&FrictionLaw::ShiftedPower {
        static_mu: 1.2,
        kinetic_mu: 0.0,
        characteristic_speed: 0.011,
        exponent: 0.9,
    })
    .unwrap();
    assert!(matches!(
        SmithWoodhouse::new(params(cal, 1e-4, 1)),
        Err(TriboError::InvalidInput {
            field: "smith_woodhouse_calibration_temperature_not_certified_increasing"
        })
    ));
    // The named calibration refuses a rig solid it cannot use.
    let mut cold = PERSPEX;
    cold.conductivity_w_per_m_k = 0.0;
    assert!(RosinCalibration::woodhouse_galluzzo_2025_eq13(cold).is_err());
    // A curve with kinetic above static refuses inside the params too.
    let mut p = params(
        RosinCalibration::woodhouse_galluzzo_2025_eq13(PERSPEX).unwrap(),
        1e-4,
        1,
    );
    p.steady_fit = SteadySlidingCurve::ShiftedPower {
        static_mu: 0.5,
        kinetic_mu: 0.6,
        characteristic_speed_mps: 0.01,
        exponent: 0.4,
    };
    assert!(matches!(
        SmithWoodhouse::new(p),
        Err(TriboError::InvalidInput {
            field: "smith_woodhouse_curve_kinetic_mu_above_static"
        })
    ));
    println!("{{\"law\":\"rosin-thermal\",\"test\":\"calibration_refusals\",\"pass\":true}}");
}
