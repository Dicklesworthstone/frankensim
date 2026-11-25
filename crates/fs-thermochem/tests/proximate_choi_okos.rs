//! G0/G2 ASHRAE Ch. 19 / Choi–Okos (1986) proximate thermophysical oracles,
//! Kirchhoff-consistent freezing enthalpy conservation, JANAF-298.15K datum
//! alignment, and planted-negative boundary tests.

use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Dims, Pressure, Temperature};
use fs_thermochem::MassSpecificEnergyQuantityV1;
use fs_thermochem::caloric::proximate::{
    CHOI_OKOS_LAW_ID_V1, CHOI_OKOS_LAW_VERSION_V1, ChoiOkosEvaluator, ChoiOkosTable,
    ConductivityMixingModel, DEFAULT_MUSHY_ZONE_WIDTH_K, H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG,
    ICE_FUSION_ENTHALPY_0C_J_KG, ICE_POINT_K, JANAF_FORMATION_DATUM_ID_V1,
    JANAF_REFERENCE_PRESSURE_PA, JANAF_REFERENCE_TEMPERATURE_K, MAXIMUM_ADMITTED_TEMPERATURE_K,
    MAXIMUM_VALID_TEMPERATURE_K, MINIMUM_ADMITTED_TEMPERATURE_K, MINIMUM_VALID_TEMPERATURE_K,
    ProximateComposition, ProximateErrorV1, ProximateEvalConfig, ProximateFractions,
    ValidityStatus, classify_temperature_validity, estimate_freezing_point_raoult,
    water_enthalpy_janaf_from_if97, water_if97_janaf_offset,
    water_if97_reference_enthalpy_298k_1bar, water_state_enthalpy_janaf,
};
use fs_thermochem::water::if97::{If97Region, If97SaturationPair, If97State};

fn celsius(t_c: f64) -> Temperature {
    Temperature::new(ICE_POINT_K + t_c)
}

fn lean_beef() -> ProximateComposition {
    // ASHRAE Ch. 19 lean beef (~74% water, 21% protein, 4% fat, 1% ash), T_f = -1.7 C.
    ProximateComposition::new(
        ProximateFractions {
            water: 0.74,
            protein: 0.21,
            fat: 0.04,
            carbohydrate: 0.0,
            fiber: 0.0,
            ash: 0.01,
        },
        0.084,
        ICE_POINT_K - 1.7,
    )
    .expect("lean beef composition should admit")
}

fn chicken_breast() -> ProximateComposition {
    // ASHRAE Ch. 19 lean chicken breast (74% water, 23% protein, 1.5% fat, 1.5% ash), T_f = -2.2 C.
    ProximateComposition::new(
        ProximateFractions {
            water: 0.74,
            protein: 0.23,
            fat: 0.015,
            carbohydrate: 0.0,
            fiber: 0.0,
            ash: 0.015,
        },
        0.092,
        ICE_POINT_K - 2.2,
    )
    .expect("chicken breast composition should admit")
}

fn russet_potato() -> ProximateComposition {
    // ASHRAE Ch. 19 raw potato (79% water, 2% protein, 0.1% fat, 16% carb, 1.8% fiber, 1.1% ash), T_f = -1.1 C.
    ProximateComposition::new(
        ProximateFractions {
            water: 0.79,
            protein: 0.02,
            fat: 0.001,
            carbohydrate: 0.16,
            fiber: 0.018,
            ash: 0.011,
        },
        0.008,
        ICE_POINT_K - 1.1,
    )
    .expect("russet potato composition should admit")
}

fn butter() -> ProximateComposition {
    // Unsalted butter (~81.1% fat, 17.5% water, 0.85% protein, 0.06% carb, 0.49% ash), T_f = -0.5 C.
    ProximateComposition::new(
        ProximateFractions {
            water: 0.175,
            protein: 0.0085,
            fat: 0.811,
            carbohydrate: 0.0006,
            fiber: 0.0,
            ash: 0.0049,
        },
        0.0034,
        ICE_POINT_K - 0.5,
    )
    .expect("butter composition should admit")
}

fn white_bread() -> ProximateComposition {
    // ASHRAE white bread crumb (~37% water, 9% protein, 3% fat, 47% carb, 2% fiber, 2% ash, 55% gas porosity).
    ProximateComposition::new(
        ProximateFractions {
            water: 0.37,
            protein: 0.09,
            fat: 0.03,
            carbohydrate: 0.47,
            fiber: 0.02,
            ash: 0.02,
        },
        0.036,
        ICE_POINT_K - 2.5,
    )
    .expect("white bread composition should admit")
    .with_porosity_gas(0.55)
    .expect("bread porosity should admit")
}

fn vanilla_ice_cream() -> ProximateComposition {
    // Vanilla ice cream (~61% water, 3.5% protein, 11% fat, 23.6% carb, 0.9% ash, T_f = -5.6 C, 40% air overrun).
    ProximateComposition::new(
        ProximateFractions {
            water: 0.61,
            protein: 0.035,
            fat: 0.11,
            carbohydrate: 0.236,
            fiber: 0.0,
            ash: 0.009,
        },
        0.014,
        ICE_POINT_K - 5.6,
    )
    .expect("ice cream composition should admit")
    .with_porosity_gas(0.40)
    .expect("ice cream porosity should admit")
}

/// 5-point Gauss–Legendre quadrature over `[a, b]` subdivided into `panels` uniform panels.
fn integrate_gl5<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, panels: usize) -> f64 {
    const NODES: [f64; 5] = [
        -0.906_179_845_938_664,
        -0.538_469_310_105_683_1,
        0.0,
        0.538_469_310_105_683_1,
        0.906_179_845_938_664,
    ];
    const WEIGHTS: [f64; 5] = [
        0.236_926_885_056_189_1,
        0.478_628_670_499_366_5,
        0.568_888_888_888_888_9,
        0.478_628_670_499_366_5,
        0.236_926_885_056_189_1,
    ];
    let h = (b - a) / (panels as f64);
    let mut sum = 0.0;
    for i in 0..panels {
        let x0 = a + (i as f64) * h;
        let mid = x0 + 0.5 * h;
        let half = 0.5 * h;
        let mut panel_sum = 0.0;
        for (&xi, &wi) in NODES.iter().zip(WEIGHTS.iter()) {
            panel_sum += wi * f(mid + half * xi);
        }
        sum += half * panel_sum;
    }
    sum
}

/// Integrate `c_p_apparent(T)` piecewise across all nonsmooth / transition knots
/// (`T_f - delta_T`, `T_f`, and `273.15 K`) so high-order quadrature resolves
/// both the mushy transition and the `0 C` Choi–Okos liquid-water polynomial seam.
fn integrate_cp_apparent_piecewise(
    eval: &ChoiOkosEvaluator,
    comp: &ProximateComposition,
    t1_k: f64,
    t2_k: f64,
) -> f64 {
    let mut knots = vec![t1_k, t2_k];
    for candidate in [
        comp.initial_freezing_point_k() - eval.config().mushy_zone_width_k,
        comp.initial_freezing_point_k(),
        ICE_POINT_K,
    ] {
        if candidate > t1_k && candidate < t2_k {
            knots.push(candidate);
        }
    }
    knots.sort_by(f64::total_cmp);
    knots.dedup();
    let mut total = 0.0;
    for seg in knots.windows(2) {
        let panels = if (seg[1] - seg[0]) <= 0.2 { 64 } else { 128 };
        total += integrate_gl5(
            |t| {
                eval.evaluate(comp, Temperature::new(t))
                    .expect("interior quadrature point should admit")
                    .cp_apparent()
                    .value()
            },
            seg[0],
            seg[1],
            panels,
        );
    }
    total
}

fn assert_within_rel(
    food: &str,
    t_c: f64,
    prop: &str,
    calc: f64,
    ashrae_ref: f64,
    max_rel_err: f64,
) {
    let rel_err = (calc - ashrae_ref).abs() / ashrae_ref.abs();
    assert!(
        rel_err <= max_rel_err,
        "ASHRAE oracle mismatch: food={food}, T_C={t_c}, prop={prop}, calc={calc:.6}, ashrae_ref={ashrae_ref:.6}, rel_err={rel_err:.4e} (limit {max_rel_err:.4e})"
    );
}

#[test]
fn ashrae_ch19_golden_oracle_six_reference_foods_and_pure_water() {
    let eval = ChoiOkosEvaluator::ashrae_default();

    // 1. Lean beef (~74% water, unfrozen at 20 C & 95 C, frozen at -30 C).
    let beef = lean_beef();
    let beef_20 = eval.evaluate(&beef, celsius(20.0)).expect("beef 20 C");
    assert_eq!(beef_20.validity(), ValidityStatus::Valid);
    assert_eq!(beef_20.ice_mass_fraction().to_bits(), 0.0_f64.to_bits());
    assert!((beef.solids() - 0.26).abs() < 1e-12);
    assert!(beef_20.thermal_diffusivity_sensible_m2_s() > 1.0e-7);
    assert_eq!(
        beef_20.thermal_diffusivity_sensible_m2_s().to_bits(),
        beef_20.thermal_diffusivity_apparent_m2_s().to_bits()
    );
    assert_within_rel(
        "lean_beef",
        20.0,
        "rho",
        beef_20.density().value(),
        1053.0,
        0.03,
    );
    assert_within_rel(
        "lean_beef",
        20.0,
        "cp",
        beef_20.cp_apparent().value(),
        3609.0,
        0.03,
    );
    assert_within_rel(
        "lean_beef",
        20.0,
        "k_me",
        beef_20.k_maxwell_eucken().value(),
        0.49,
        0.08,
    );

    let beef_95 = eval.evaluate(&beef, celsius(95.0)).expect("beef 95 C");
    assert_eq!(beef_95.validity(), ValidityStatus::Valid);
    assert_within_rel(
        "lean_beef",
        95.0,
        "cp",
        beef_95.cp_apparent().value(),
        3660.0,
        0.03,
    );
    assert!(beef_95.k_effective().value() > beef_20.k_effective().value());

    let beef_m30 = eval.evaluate(&beef, celsius(-30.0)).expect("beef -30 C");
    assert_eq!(beef_m30.validity(), ValidityStatus::Valid);
    let expected_xi_beef_m30 = (0.74 - 0.084) * (1.0 - (-1.7) / (-30.0));
    assert_within_rel(
        "lean_beef",
        -30.0,
        "X_i",
        beef_m30.ice_mass_fraction(),
        expected_xi_beef_m30,
        1e-12,
    );
    // Ice conductivity (~2.5 W/m K) makes frozen beef > 3x more conductive than unfrozen beef
    // (ASHRAE Ch. 19 lean beef at -30 C parallel/isotropic ~1.57 W/(m K)).
    assert_within_rel(
        "lean_beef",
        -30.0,
        "k_me",
        beef_m30.k_maxwell_eucken().value(),
        1.57,
        0.08,
    );
    assert!(beef_m30.density().value() < beef_20.density().value());

    // 2. Chicken breast (74% water).
    let chicken = chicken_breast();
    let ch_20 = eval
        .evaluate(&chicken, celsius(20.0))
        .expect("chicken 20 C");
    assert_within_rel(
        "chicken_breast",
        20.0,
        "rho",
        ch_20.density().value(),
        1066.0,
        0.03,
    );
    assert_within_rel(
        "chicken_breast",
        20.0,
        "cp",
        ch_20.cp_apparent().value(),
        3605.0,
        0.03,
    );
    assert_within_rel(
        "chicken_breast",
        20.0,
        "k_me",
        ch_20.k_maxwell_eucken().value(),
        0.495,
        0.08,
    );

    let ch_m20 = eval
        .evaluate(&chicken, celsius(-20.0))
        .expect("chicken -20 C");
    assert!(ch_m20.ice_mass_fraction() > 0.55);
    assert!(ch_m20.k_effective().value() > 2.0 * ch_20.k_effective().value());

    // 3. Russet potato (79% water).
    let pot = russet_potato();
    let pot_20 = eval.evaluate(&pot, celsius(20.0)).expect("potato 20 C");
    assert_within_rel(
        "russet_potato",
        20.0,
        "rho",
        pot_20.density().value(),
        1076.0,
        0.03,
    );
    assert_within_rel(
        "russet_potato",
        20.0,
        "cp",
        pot_20.cp_apparent().value(),
        3643.0,
        0.03,
    );
    assert_within_rel(
        "russet_potato",
        20.0,
        "k_me",
        pot_20.k_maxwell_eucken().value(),
        0.535,
        0.08,
    );
    let pot_m30 = eval.evaluate(&pot, celsius(-30.0)).expect("potato -30 C");
    assert!(pot_m30.ice_mass_fraction() > 0.72);
    assert!(pot_m30.k_effective().value() > 1.5);

    // 4. Butter (~81% fat).
    let but = butter();
    let but_20 = eval.evaluate(&but, celsius(20.0)).expect("butter 20 C");
    assert_within_rel("butter", 20.0, "rho", but_20.density().value(), 932.0, 0.03);
    assert_within_rel(
        "butter",
        20.0,
        "cp",
        but_20.cp_apparent().value(),
        2386.0,
        0.03,
    );
    assert_within_rel(
        "butter",
        20.0,
        "k_me",
        but_20.k_maxwell_eucken().value(),
        0.214,
        0.08,
    );
    // Fat dominates butter conductivity so butter k is much lower than lean beef k.
    assert!(but_20.k_effective().value() < 0.5 * beef_20.k_effective().value());

    // 5. White bread (porous crumb, phi_g = 0.55).
    let bread = white_bread();
    let br_20 = eval.evaluate(&bread, celsius(20.0)).expect("bread 20 C");
    assert_within_rel(
        "white_bread",
        20.0,
        "rho_bulk",
        br_20.density().value(),
        570.0,
        0.03,
    );
    assert_within_rel(
        "white_bread",
        20.0,
        "cp",
        br_20.cp_apparent().value(),
        2595.0,
        0.03,
    );
    assert_within_rel(
        "white_bread",
        20.0,
        "k_me",
        br_20.k_maxwell_eucken().value(),
        0.16,
        0.08,
    );

    // 6. Vanilla ice cream (depressed T_f = -5.6 C, 40% air overrun).
    let ic = vanilla_ice_cream();
    let ic_5 = eval.evaluate(&ic, celsius(5.0)).expect("ice cream 5 C");
    assert_eq!(ic_5.ice_mass_fraction().to_bits(), 0.0_f64.to_bits());
    assert_within_rel(
        "ice_cream",
        5.0,
        "cp",
        ic_5.cp_apparent().value(),
        3212.0,
        0.03,
    );
    assert_within_rel(
        "ice_cream",
        5.0,
        "k_me",
        ic_5.k_maxwell_eucken().value(),
        0.24,
        0.08,
    );
    // At -3 C (above T_f = -5.6 C), ice cream is still strictly unfrozen!
    let ic_m3 = eval.evaluate(&ic, celsius(-3.0)).expect("ice cream -3 C");
    assert_eq!(ic_m3.ice_mass_fraction().to_bits(), 0.0_f64.to_bits());
    // At -20 C, ice fraction follows exact Chen-Heldman: (0.61 - 0.014) * (1 - 5.6 / 20.0).
    let ic_m20 = eval.evaluate(&ic, celsius(-20.0)).expect("ice cream -20 C");
    assert_within_rel(
        "ice_cream",
        -20.0,
        "X_i",
        ic_m20.ice_mass_fraction(),
        (0.61 - 0.014) * (1.0 - 5.6 / 20.0),
        1e-12,
    );
    assert!(ic_m20.k_effective().value() > 1.8 * ic_5.k_effective().value());
}

#[test]
fn pure_water_and_near_pure_water_mushy_zone_never_nan_or_blow_up() {
    let eval = ChoiOkosEvaluator::ashrae_default();
    let water = ProximateComposition::pure_water();

    // Sweep across frozen, mushy boundary, exact 0 C, and hot liquid water.
    for t_c in [
        -30.0,
        -1.0,
        -DEFAULT_MUSHY_ZONE_WIDTH_K,
        -0.5 * DEFAULT_MUSHY_ZONE_WIDTH_K,
        -1e-6,
        -1e-12,
        0.0,
        1e-12,
        20.0,
        95.0,
    ] {
        let res = eval
            .evaluate(&water, celsius(t_c))
            .expect("pure water should evaluate cleanly without NaN");
        assert!(res.ice_mass_fraction().is_finite(), "X_i finite at {t_c} C");
        assert!(
            res.d_ice_fraction_d_t().is_finite(),
            "dX_i/dT finite at {t_c} C"
        );
        assert!(res.d_ice_fraction_d_t() <= 0.0, "dX_i/dT <= 0 at {t_c} C");
        assert!(res.cp_apparent().value().is_finite() && res.cp_apparent().value() > 0.0);
        assert!(res.specific_enthalpy_0c_datum().value().is_finite());
        assert!(res.specific_enthalpy_janaf().value().is_finite());
        if t_c >= 0.0 {
            assert_eq!(
                res.ice_mass_fraction().to_bits(),
                0.0_f64.to_bits(),
                "strictly zero ice at T >= T_f ({t_c} C)"
            );
            assert_eq!(
                res.d_ice_fraction_d_t().to_bits(),
                0.0_f64.to_bits(),
                "strictly zero dX_i/dT at T >= T_f ({t_c} C)"
            );
        }
        if t_c <= -DEFAULT_MUSHY_ZONE_WIDTH_K {
            assert_eq!(
                res.ice_mass_fraction().to_bits(),
                1.0_f64.to_bits(),
                "100% ice at T <= T_f - delta_T ({t_c} C)"
            );
            assert_eq!(
                res.d_ice_fraction_d_t().to_bits(),
                0.0_f64.to_bits(),
                "zero dX_i/dT in pure solid ice ({t_c} C)"
            );
        }
    }

    // Pure-water exact limits: solid ice at -30 C and -1 C, mushy midpoint at -0.025 C, liquid at 0/20/95 C.
    let w_m30 = eval.evaluate(&water, celsius(-30.0)).unwrap();
    assert_within_rel(
        "pure_water",
        -30.0,
        "rho_ice",
        w_m30.density().value(),
        920.8113,
        1e-6,
    );
    assert_within_rel(
        "pure_water",
        -30.0,
        "cp_ice",
        w_m30.cp_apparent().value(),
        1879.993,
        1e-6,
    );
    assert_within_rel(
        "pure_water",
        -30.0,
        "k_ice",
        w_m30.k_effective().value(),
        2.498_453,
        1e-6,
    );

    let t_mid = celsius(-0.5 * DEFAULT_MUSHY_ZONE_WIDTH_K);
    let w_mid = eval.evaluate(&water, t_mid).unwrap();
    // Exact Hermite value at the binary64-represented temperature (ICE_POINT_K = 273.15 is not
    // dyadic, so 273.15 - 273.125 = 0.025 - 2.27e-14, shifting s by 4.55e-13 and f by 6.82e-13).
    let s_rep = (ICE_POINT_K - t_mid.value()) / DEFAULT_MUSHY_ZONE_WIDTH_K;
    let f_rep = s_rep * s_rep * (3.0 - 2.0 * s_rep);
    assert!((w_mid.ice_mass_fraction() - f_rep).abs() < 1e-15);
    assert!((w_mid.ice_mass_fraction() - 0.5).abs() < 1e-12);
    assert!((w_mid.d_ice_fraction_d_t() - (-1.5 / DEFAULT_MUSHY_ZONE_WIDTH_K)).abs() < 1e-12);

    let w_0 = eval.evaluate(&water, celsius(0.0)).unwrap();
    assert_within_rel(
        "pure_water",
        0.0,
        "cp_w",
        w_0.cp_apparent().value(),
        4176.2,
        1e-12,
    );
    assert_within_rel(
        "pure_water",
        0.0,
        "k_w",
        w_0.k_effective().value(),
        0.57109,
        1e-12,
    );
    assert_eq!(
        w_0.specific_enthalpy_0c_datum().value().to_bits(),
        0.0_f64.to_bits()
    );

    let w_20 = eval.evaluate(&water, celsius(20.0)).unwrap();
    assert_within_rel(
        "pure_water",
        20.0,
        "rho_w",
        w_20.density().value(),
        995.74,
        1e-4,
    );
    assert_within_rel(
        "pure_water",
        20.0,
        "cp_w",
        w_20.cp_apparent().value(),
        4176.57,
        1e-4,
    );
    let w_95 = eval.evaluate(&water, celsius(95.0)).unwrap();
    assert_within_rel(
        "pure_water",
        95.0,
        "cp_w",
        w_95.cp_apparent().value(),
        4216.96,
        1e-4,
    );

    // Near-pure water (0 < 273.15 - T_f < delta_T): T_f = -0.01 C, delta_T = 0.05 K.
    let near_pure = ProximateComposition::pure_water()
        .with_freezing_point(celsius(-0.01))
        .unwrap();
    let mut prev_xi = 2.0;
    for step in 0..=200 {
        let t_c = -0.10 + 0.12 * (step as f64) / 200.0;
        let (xi, dxi) = eval
            .ice_fraction_and_derivative(&near_pure, celsius(t_c))
            .unwrap();
        assert!(xi.is_finite() && (0.0..=1.0).contains(&xi));
        assert!(dxi.is_finite() && dxi <= 0.0);
        assert!(
            xi <= prev_xi + 1e-15,
            "X_i must be non-increasing at t_c={t_c}"
        );
        prev_xi = xi;
    }
    // Check C^1 continuity of X_i and dX_i/dT at the lower mushy seam t_m = -0.06 C.
    let eps = 1e-9;
    let (xi_lo, dxi_lo) = eval
        .ice_fraction_and_derivative(&near_pure, celsius(-0.06 - eps))
        .unwrap();
    let (xi_hi, dxi_hi) = eval
        .ice_fraction_and_derivative(&near_pure, celsius(-0.06 + eps))
        .unwrap();
    assert!((xi_lo - xi_hi).abs() < 1e-6);
    assert!((dxi_lo - dxi_hi).abs() < 1e-4);
}

#[test]
fn kirchhoff_enthalpy_integral_and_analytical_cp_apparent_closure() {
    let eval = ChoiOkosEvaluator::ashrae_default();

    // 1. Verify Kirchhoff law: L_f(T) drops ~33% from 0 C to -40 C (NOT constant 333.6 kJ/kg!).
    let table = eval.table();
    let lf_0 = table.latent_heat_fusion(0.0);
    let lf_m40 = table.latent_heat_fusion(-40.0);
    assert_eq!(lf_0.to_bits(), ICE_FUSION_ENTHALPY_0C_J_KG.to_bits());
    // Also verify QuadraticCelsiusPoly::derivative matches central difference.
    let poly_fd = (table.protein.cp.eval(40.0 + 1e-5) - table.protein.cp.eval(40.0 - 1e-5)) / 2e-5;
    assert!((poly_fd - table.protein.cp.derivative(40.0)).abs() < 1e-7);
    assert!(
        lf_m40 < 0.75 * lf_0 && lf_m40 > 0.60 * lf_0,
        "Kirchhoff L_f(-40 C)={lf_m40:.1} J/kg should drop ~33% below L_f(0 C)={lf_0:.1} J/kg"
    );
    // Verify dL_f/dT == c_{p,w}(T) - c_{p,ice}(T) via central finite difference at -20 C.
    let dt = 1e-5;
    let dlf_fd =
        (table.latent_heat_fusion(-20.0 + dt) - table.latent_heat_fusion(-20.0 - dt)) / (2.0 * dt);
    let dlf_exact = table.water_cp(-20.0) - table.ice_cp.eval(-20.0);
    assert!(
        (dlf_fd - dlf_exact).abs() / dlf_exact.abs() < 1e-8,
        "Kirchhoff derivative mismatch: fd={dlf_fd}, exact={dlf_exact}"
    );

    // 2. Verify analytical c_p_apparent(T) == dh/dT against central finite differences
    //    AND integral_{T1}^{T2} c_p_apparent(T) dT == h(T2) - h(T1) to < 1e-10 across T_f.
    let foods = [
        ("lean_beef", lean_beef(), -40.0, 95.0),
        ("russet_potato", russet_potato(), -30.0, 95.0),
        ("butter", butter(), -30.0, 95.0),
        ("vanilla_ice_cream", vanilla_ice_cream(), -30.0, 20.0),
        (
            "pure_water",
            ProximateComposition::pure_water(),
            -30.0,
            95.0,
        ),
    ];

    for (name, comp, t1_c, t2_c) in foods {
        let t1_k = ICE_POINT_K + t1_c;
        let t2_k = ICE_POINT_K + t2_c;
        let h1 = eval.evaluate(&comp, Temperature::new(t1_k)).unwrap();
        let h2 = eval.evaluate(&comp, Temperature::new(t2_k)).unwrap();

        let dh_0c =
            h2.specific_enthalpy_0c_datum().value() - h1.specific_enthalpy_0c_datum().value();
        let dh_janaf = h2.specific_enthalpy_janaf().value() - h1.specific_enthalpy_janaf().value();
        let integral_cp = integrate_cp_apparent_piecewise(&eval, &comp, t1_k, t2_k);

        let rel_err = (integral_cp - dh_0c).abs() / dh_0c.abs();
        assert!(
            rel_err < 1e-10,
            "Kirchhoff integral mismatch for {name} [{t1_c}..{t2_c} C]: integral={integral_cp:.6}, dh={dh_0c:.6}, rel_err={rel_err:.3e}"
        );
        // JANAF shift is strictly temperature-independent across freezing!
        let datum_diff = (dh_janaf - dh_0c).abs() / dh_0c.abs();
        assert!(
            datum_diff < 1e-9,
            "JANAF datum shift leaked temperature dependence for {name}: diff={datum_diff:.3e}"
        );

        // Pointwise central finite-difference verification of cp_apparent = dh/dT
        // across frozen, mushy-interior, and unfrozen points (avoiding the isolated 0 C seam).
        let tf_c = comp.initial_freezing_point_k() - ICE_POINT_K;
        let probes = [
            -25.0,
            tf_c - 0.5 * DEFAULT_MUSHY_ZONE_WIDTH_K,
            tf_c - 0.5,
            20.0,
            80.0,
        ];
        for probe_c in probes {
            if (probe_c).abs() < 1e-3 {
                continue;
            }
            let h_step = 1e-5;
            let p_minus = eval.evaluate(&comp, celsius(probe_c - h_step)).unwrap();
            let p_mid = eval.evaluate(&comp, celsius(probe_c)).unwrap();
            let p_plus = eval.evaluate(&comp, celsius(probe_c + h_step)).unwrap();
            let fd_cp = (p_plus.specific_enthalpy_0c_datum().value()
                - p_minus.specific_enthalpy_0c_datum().value())
                / (2.0 * h_step);
            let exact_cp = p_mid.cp_apparent().value();
            let cp_rel = (fd_cp - exact_cp).abs() / exact_cp.abs();
            assert!(
                cp_rel < 1e-6,
                "analytical cp_apparent mismatch for {name} at {probe_c:.4} C: fd={fd_cp:.6}, exact={exact_cp:.6}, rel={cp_rel:.3e}"
            );
        }
    }

    // ASHRAE -40 C reference enthalpy accessor check: zero at -40 C, matches h(20 C) - h(-40 C).
    let beef = lean_beef();
    assert_eq!(
        eval.enthalpy_ashrae_minus40c(&beef, celsius(-40.0))
            .unwrap()
            .value()
            .to_bits(),
        0.0_f64.to_bits()
    );
    let h_beef_20_ashrae = eval
        .enthalpy_ashrae_minus40c(&beef, celsius(20.0))
        .unwrap()
        .value();
    // ASHRAE Ch. 19 lean beef enthalpy from -40 C to +20 C is ~360-390 kJ/kg (sensible + latent).
    assert!(
        h_beef_20_ashrae > 330_000.0 && h_beef_20_ashrae < 420_000.0,
        "lean beef enthalpy (-40 C to +20 C) out of physical range: {h_beef_20_ashrae:.1} J/kg"
    );
}

#[test]
fn janaf_298k_formation_datum_aligns_choi_okos_and_if97_streams() {
    assert_eq!(JANAF_FORMATION_DATUM_ID_V1, "janaf-298.15K-formation");
    // 1. Liquid H2O formation enthalpy constant check: -15.866 MJ/kg.
    assert!(
        (H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG - (-15.866e6)).abs() < 1_000.0,
        "expected H2O(l) h_f0(298.15 K) ~= -15.866 MJ/kg, got {H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG:.2}"
    );

    // 2. Pure water Choi-Okos at 298.15 K equals exact H2O(l) formation enthalpy.
    let eval = ChoiOkosEvaluator::ashrae_default();
    let water = ProximateComposition::pure_water();
    let w_298 = eval
        .evaluate(&water, Temperature::new(JANAF_REFERENCE_TEMPERATURE_K))
        .unwrap();
    assert!(
        (w_298.specific_enthalpy_janaf().value() - H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG).abs()
            < 1e-6
    );

    // 3. IF97 liquid water at (298.15 K, 1 bar) shifted via water_state_enthalpy_janaf
    //    matches H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG identically!
    let if97_ref = water_if97_reference_enthalpy_298k_1bar().unwrap();
    assert!(
        if97_ref.value() > 104_000.0 && if97_ref.value() < 106_000.0,
        "IF97 h(298.15 K, 1 bar) expected ~104.9 kJ/kg, got {:.2}",
        if97_ref.value()
    );
    let if97_298 = If97State::at_pressure_temperature(
        Pressure::new(JANAF_REFERENCE_PRESSURE_PA),
        Temperature::new(JANAF_REFERENCE_TEMPERATURE_K),
        If97Region::Region1,
    )
    .unwrap();
    let h_if97_298_janaf = water_state_enthalpy_janaf(&if97_298).unwrap().value();
    assert!(
        (h_if97_298_janaf - H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG).abs() < 1e-8,
        "IF97 JANAF shift at 298.15 K mismatch: {h_if97_298_janaf} vs {H2O_LIQUID_FORMATION_ENTHALPY_298K_J_KG}"
    );

    // 4. Cross-check IF97 saturated steam at 298.15 K shifted to JANAF datum against
    //    NIST-JANAF H2O(g) formation enthalpy (-241.826 kJ/mol / 0.018015268 = -13.42339 MJ/kg).
    let sat_298 =
        If97SaturationPair::at_temperature(Temperature::new(JANAF_REFERENCE_TEMPERATURE_K))
            .unwrap();
    let h_steam_janaf = water_state_enthalpy_janaf(sat_298.vapor()).unwrap().value();
    let nist_h2o_gas_hf0 = -241_826.0 / 0.018_015_268;
    let steam_rel = (h_steam_janaf - nist_h2o_gas_hf0).abs() / nist_h2o_gas_hf0.abs();
    assert!(
        steam_rel < 5e-4,
        "IF97-JANAF vapor enthalpy {h_steam_janaf:.1} vs NIST H2O(g) {nist_h2o_gas_hf0:.1} (rel {steam_rel:.4e})"
    );

    // 5. Cross-check Choi-Okos pure water vs IF97 liquid water at 1 bar from 5 C to 95 C.
    for t_c in [5.0, 25.0, 60.0, 95.0] {
        let t = celsius(t_c);
        let co_h = eval
            .evaluate(&water, t)
            .unwrap()
            .specific_enthalpy_janaf()
            .value();
        let if97_s = If97State::at_pressure_temperature(
            Pressure::new(JANAF_REFERENCE_PRESSURE_PA),
            t,
            If97Region::Region1,
        )
        .unwrap();
        let if97_h = water_state_enthalpy_janaf(&if97_s).unwrap().value();
        assert!(
            (co_h - if97_h).abs() < 1_000.0,
            "Choi-Okos vs IF97 JANAF liquid water at {t_c} C diverged: co={co_h:.1}, if97={if97_h:.1}"
        );
    }

    // 6. Food stream offset scales linearly with water mass fraction X_w.
    let beef = lean_beef();
    let b_298 = eval
        .evaluate(&beef, Temperature::new(JANAF_REFERENCE_TEMPERATURE_K))
        .unwrap();
    let expected_shift = beef.water() * eval.table().water_janaf_datum_offset_j_kg();
    assert!(
        ((b_298.specific_enthalpy_janaf().value() - b_298.specific_enthalpy_0c_datum().value())
            - expected_shift)
            .abs()
            < 1e-6
    );
    assert!(water_if97_janaf_offset().unwrap().value() < -15.9e6);
}

#[test]
fn conductivity_models_respect_wiener_bounds_and_pure_limits() {
    let base_eval = ChoiOkosEvaluator::ashrae_default();
    let foods = [
        lean_beef(),
        chicken_breast(),
        russet_potato(),
        butter(),
        white_bread(),
        vanilla_ice_cream(),
        ProximateComposition::pure_water(),
    ];

    for comp in &foods {
        for t_c in [-30.0, -5.0, 20.0, 95.0] {
            let r = base_eval.evaluate(comp, celsius(t_c)).unwrap();
            let k_ser = r.k_series().value();
            let k_par = r.k_parallel().value();
            let k_me = r.k_maxwell_eucken().value();
            let k_kop = r.k_kopelman().value();
            assert!(
                k_ser <= k_par + 1e-14,
                "Wiener order violated: k_ser={k_ser} > k_par={k_par}"
            );
            assert!(
                k_me >= k_ser - 1e-13 && k_me <= k_par + 1e-13,
                "Maxwell-Eucken outside Wiener bounds: k_ser={k_ser}, k_me={k_me}, k_par={k_par}"
            );
            assert!(
                k_kop >= k_ser - 1e-13 && k_kop <= k_par + 1e-13,
                "Kopelman outside Wiener bounds: k_ser={k_ser}, k_kop={k_kop}, k_par={k_par}"
            );

            // Check Geometric and Krischer sweeps across anisotropy_factor in [0, 1].
            let geom_eval = base_eval
                .clone()
                .with_config(ProximateEvalConfig {
                    mushy_zone_width_k: DEFAULT_MUSHY_ZONE_WIDTH_K,
                    conductivity_model: ConductivityMixingModel::Geometric,
                })
                .unwrap();
            let k_geom = geom_eval
                .evaluate(comp, celsius(t_c))
                .unwrap()
                .k_effective()
                .value();
            assert!(k_geom >= k_ser - 1e-13 && k_geom <= k_par + 1e-13);

            for (f, expect_par, expect_ser) in [
                (0.0, true, false),
                (0.25, false, false),
                (0.5, false, false),
                (0.75, false, false),
                (1.0, false, true),
            ] {
                let kr_eval = base_eval
                    .clone()
                    .with_config(ProximateEvalConfig {
                        mushy_zone_width_k: DEFAULT_MUSHY_ZONE_WIDTH_K,
                        conductivity_model: ConductivityMixingModel::Krischer {
                            anisotropy_factor: f,
                        },
                    })
                    .unwrap();
                let k_kr = kr_eval
                    .evaluate(comp, celsius(t_c))
                    .unwrap()
                    .k_effective()
                    .value();
                assert!(k_kr >= k_ser - 1e-13 && k_kr <= k_par + 1e-13);
                if expect_par {
                    assert!((k_kr - k_par).abs() < 1e-14);
                }
                if expect_ser {
                    assert!((k_kr - k_ser).abs() < 1e-14);
                }
            }
        }
    }

    // Pure component limit: all conductivity mixing models collapse identically.
    let pure_fat = ProximateComposition::new(
        ProximateFractions {
            water: 0.0,
            protein: 0.0,
            fat: 1.0,
            carbohydrate: 0.0,
            fiber: 0.0,
            ash: 0.0,
        },
        0.0,
        ICE_POINT_K,
    )
    .unwrap();
    let rf = base_eval.evaluate(&pure_fat, celsius(100.0)).unwrap();
    assert!((rf.k_parallel().value() - rf.k_series().value()).abs() < 1e-15);
    assert!((rf.k_maxwell_eucken().value() - rf.k_parallel().value()).abs() < 1e-15);
    assert!((rf.k_kopelman().value() - rf.k_parallel().value()).abs() < 1e-15);
    // Verify ASHRAE erratum-corrected fat thermal conductivity at 100 C is positive (~0.151 W/m K).
    assert!((rf.k_effective().value() - 0.151_331_1).abs() < 1e-5);
}

#[test]
fn validity_warnings_and_hard_range_refusals_and_card_roundtrip() {
    let eval = ChoiOkosEvaluator::ashrae_default();
    let beef = lean_beef();

    // Calibration bounds [-40 C, +150 C] inclusive -> Valid.
    assert_eq!(
        eval.evaluate(&beef, Temperature::new(MINIMUM_VALID_TEMPERATURE_K))
            .unwrap()
            .validity(),
        ValidityStatus::Valid
    );
    assert_eq!(
        eval.evaluate(&beef, Temperature::new(MAXIMUM_VALID_TEMPERATURE_K))
            .unwrap()
            .validity(),
        ValidityStatus::Valid
    );

    // Below calibration [-50 C, -40 C) -> BelowCalibrationWarning.
    assert_eq!(
        eval.evaluate(&beef, Temperature::new(MINIMUM_ADMITTED_TEMPERATURE_K))
            .unwrap()
            .validity(),
        ValidityStatus::BelowCalibrationWarning
    );

    // Above 150 C crust desiccation (+150.01 C .. +250 C) -> AboveDesiccationWarning.
    assert_eq!(
        eval.evaluate(&beef, celsius(150.01)).unwrap().validity(),
        ValidityStatus::AboveDesiccationWarning
    );
    assert_eq!(
        eval.evaluate(&beef, Temperature::new(MAXIMUM_ADMITTED_TEMPERATURE_K))
            .unwrap()
            .validity(),
        ValidityStatus::AboveDesiccationWarning
    );

    // Outside [-50 C, +250 C] -> typed TemperatureOutOfRange refusal.
    assert_eq!(
        classify_temperature_validity(celsius(-50.01)),
        ValidityStatus::RefusedOutOfRange
    );
    assert_eq!(
        classify_temperature_validity(celsius(250.01)),
        ValidityStatus::RefusedOutOfRange
    );
    assert!(matches!(
        eval.evaluate(&beef, celsius(-50.01)),
        Err(ProximateErrorV1::TemperatureOutOfRange { .. })
    ));
    assert!(matches!(
        eval.evaluate(&beef, celsius(250.01)),
        Err(ProximateErrorV1::TemperatureOutOfRange { .. })
    ));

    // Card roundtrip via fs-matdb ConstitutiveModelCard.
    let ashrae = ChoiOkosTable::ashrae_1986();
    let card = ConstitutiveModelCard {
        law: LawId(CHOI_OKOS_LAW_ID_V1.to_string()),
        law_version: CHOI_OKOS_LAW_VERSION_V1,
        parameters: ashrae.to_card_parameters(),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with(
            "T",
            MINIMUM_VALID_TEMPERATURE_K,
            MAXIMUM_VALID_TEMPERATURE_K,
        ),
        sources: Vec::new(),
        provenance: Provenance {
            source: ashrae.citation.clone(),
            license: "CC0-1.0".to_string(),
            artifact: None,
        },
    };
    let loaded = ChoiOkosTable::from_card(&card).expect("valid Choi-Okos card should load");
    let loaded_eval = ChoiOkosEvaluator::new(loaded, ProximateEvalConfig::default()).unwrap();
    let r_builtin = eval.evaluate(&beef, celsius(20.0)).unwrap();
    let r_card = loaded_eval.evaluate(&beef, celsius(20.0)).unwrap();
    assert_eq!(r_builtin, r_card);

    // Planted negative card mutations refuse with typed MaterialCard error.
    let mut bad_law = card.clone();
    bad_law.law = LawId("wrong-law".to_string());
    assert!(matches!(
        ChoiOkosTable::from_card(&bad_law),
        Err(ProximateErrorV1::MaterialCard { .. })
    ));

    let mut bad_dims = card.clone();
    bad_dims.parameters.insert(
        "fusion_enthalpy_0c".to_string(),
        LawParameter {
            value: ICE_FUSION_ENTHALPY_0C_J_KG,
            dims: Dims::NONE,
        },
    );
    assert!(matches!(
        ChoiOkosTable::from_card(&bad_dims),
        Err(ProximateErrorV1::MaterialCard { .. })
    ));

    let mut extra_param = card.clone();
    extra_param.parameters.insert(
        "rogue_tuning_knob".to_string(),
        LawParameter {
            value: 1.0,
            dims: Dims::NONE,
        },
    );
    assert!(matches!(
        ChoiOkosTable::from_card(&extra_param),
        Err(ProximateErrorV1::MaterialCard { .. })
    ));

    let mut extra_axis = card;
    extra_axis.validity = extra_axis.validity.with("P", 10_000.0, 200_000.0);
    assert!(matches!(
        ChoiOkosTable::from_card(&extra_axis),
        Err(ProximateErrorV1::MaterialCard { .. })
    ));
}

#[test]
fn planted_negative_composition_inputs_refuse() {
    let valid_fracs = ProximateFractions {
        water: 0.74,
        protein: 0.21,
        fat: 0.04,
        carbohydrate: 0.0,
        fiber: 0.0,
        ash: 0.01,
    };

    // 1. Mass fractions != 1 (outside 1e-4 tolerance).
    assert!(matches!(
        ProximateComposition::new(
            ProximateFractions {
                water: 0.70,
                ..valid_fracs
            },
            0.05,
            ICE_POINT_K - 1.5
        ),
        Err(ProximateErrorV1::CompositionSumOutOfTolerance { .. })
    ));
    assert!(matches!(
        ProximateComposition::new(
            ProximateFractions {
                water: 0.7402,
                ..valid_fracs
            },
            0.05,
            ICE_POINT_K - 1.5
        ),
        Err(ProximateErrorV1::CompositionSumOutOfTolerance { .. })
    ));

    // 2. Negative fraction.
    assert!(matches!(
        ProximateComposition::new(
            ProximateFractions {
                water: 0.75,
                protein: 0.21,
                fat: 0.04,
                carbohydrate: -0.01,
                fiber: 0.0,
                ash: 0.01,
            },
            0.05,
            ICE_POINT_K - 1.5
        ),
        Err(ProximateErrorV1::InvalidFraction {
            field: "carbohydrate",
            ..
        })
    ));

    // 3. NaN / Inf inputs in fractions, x_b, T_f, porosity.
    assert!(matches!(
        ProximateComposition::new(
            ProximateFractions {
                water: f64::NAN,
                ..valid_fracs
            },
            0.0,
            ICE_POINT_K
        ),
        Err(ProximateErrorV1::InvalidFraction { field: "water", .. })
    ));
    assert!(matches!(
        ProximateComposition::new(valid_fracs, f64::NAN, ICE_POINT_K - 1.0),
        Err(ProximateErrorV1::InvalidFraction { field: "x_b", .. })
    ));
    assert!(matches!(
        ProximateComposition::new(valid_fracs, 0.05, f64::NAN),
        Err(ProximateErrorV1::InvalidFreezingPoint { .. })
    ));
    assert!(matches!(
        lean_beef().with_porosity_gas(f64::NAN),
        Err(ProximateErrorV1::InvalidPorosity { .. })
    ));
    assert!(matches!(
        lean_beef().with_porosity_gas(1.0),
        Err(ProximateErrorV1::InvalidPorosity { .. })
    ));
    assert!(matches!(
        lean_beef().with_porosity_gas(-0.01),
        Err(ProximateErrorV1::InvalidPorosity { .. })
    ));

    // 4. T_f > 273.15 K or T_f <= 0 K.
    assert!(matches!(
        ProximateComposition::new(valid_fracs, 0.05, 273.16),
        Err(ProximateErrorV1::InvalidFreezingPoint { .. })
    ));
    assert!(matches!(
        ProximateComposition::new(valid_fracs, 0.05, 0.0),
        Err(ProximateErrorV1::InvalidFreezingPoint { .. })
    ));

    // 5. x_b > water (both constructor and with_bound_water).
    assert!(matches!(
        ProximateComposition::new(valid_fracs, 0.75, ICE_POINT_K - 1.5),
        Err(ProximateErrorV1::BoundWaterExceedsWater { .. })
    ));
    assert!(matches!(
        lean_beef().with_bound_water(0.75),
        Err(ProximateErrorV1::BoundWaterExceedsWater { .. })
    ));
    let updated_xb = lean_beef().with_bound_water(0.05).unwrap();
    assert_eq!(updated_xb.bound_water().to_bits(), 0.05_f64.to_bits());
}

#[test]
fn planted_negative_evaluation_config_raoult_and_determinism() {
    let eval = ChoiOkosEvaluator::ashrae_default();

    // 1. T <= 0 K, NaN, Inf, and non-finite IF97 enthalpy.
    for bad_t in [0.0, -10.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            eval.evaluate(&lean_beef(), Temperature::new(bad_t)),
            Err(ProximateErrorV1::InvalidTemperature { .. })
        ));
    }
    assert!(matches!(
        water_enthalpy_janaf_from_if97(MassSpecificEnergyQuantityV1::new(f64::NAN)),
        Err(ProximateErrorV1::InvalidIf97Reference { .. })
    ));

    // 2. Invalid mushy-zone width and invalid Krischer anisotropy factor.
    for bad_dt in [0.0, -0.05, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            eval.clone().with_config(ProximateEvalConfig {
                mushy_zone_width_k: bad_dt,
                ..ProximateEvalConfig::default()
            }),
            Err(ProximateErrorV1::InvalidMushyZoneWidth { .. })
        ));
    }
    for bad_f in [-0.01, 1.01, f64::NAN] {
        assert!(matches!(
            eval.clone().with_config(ProximateEvalConfig {
                mushy_zone_width_k: DEFAULT_MUSHY_ZONE_WIDTH_K,
                conductivity_model: ConductivityMixingModel::Krischer {
                    anisotropy_factor: bad_f,
                },
            }),
            Err(ProximateErrorV1::InvalidConductivityParameter { .. })
        ));
    }

    // 3. Raoult freezing-point estimator + composition helper sanity and refusals.
    let tf_pure = estimate_freezing_point_raoult(1.0, 0.0, 0.0, 0.18).unwrap();
    assert_eq!(tf_pure.value().to_bits(), ICE_POINT_K.to_bits());
    let tf_sugary = estimate_freezing_point_raoult(0.65, 0.02, 0.25, 0.342).unwrap();
    assert!(tf_sugary.value() < ICE_POINT_K - 1.0 && tf_sugary.value() > ICE_POINT_K - 5.0);
    assert!(matches!(
        estimate_freezing_point_raoult(0.65, 0.02, 0.25, 0.0),
        Err(ProximateErrorV1::InvalidFraction { .. })
    ));
    let auto_beef = ProximateComposition::from_fractions_with_raoult_freezing(
        ProximateFractions {
            water: 0.74,
            protein: 0.21,
            fat: 0.04,
            carbohydrate: 0.0,
            fiber: 0.0,
            ash: 0.01,
        },
        0.4,
        0.030,
    )
    .expect("Raoult helper should admit valid fractions");
    assert!((auto_beef.bound_water() - 0.084).abs() < 1e-12);
    assert!(
        auto_beef.initial_freezing_point_k() < ICE_POINT_K
            && auto_beef.initial_freezing_point_k() > ICE_POINT_K - 3.0
    );

    // 4. Bit-identical deterministic repeatability.
    let r1 = eval.evaluate(&lean_beef(), celsius(-12.5)).unwrap();
    let r2 = eval.evaluate(&lean_beef(), celsius(-12.5)).unwrap();
    assert_eq!(
        r1.density().value().to_bits(),
        r2.density().value().to_bits()
    );
    assert_eq!(
        r1.cp_apparent().value().to_bits(),
        r2.cp_apparent().value().to_bits()
    );
    assert_eq!(
        r1.specific_enthalpy_janaf().value().to_bits(),
        r2.specific_enthalpy_janaf().value().to_bits()
    );
    assert_eq!(
        r1.k_effective().value().to_bits(),
        r2.k_effective().value().to_bits()
    );
}
