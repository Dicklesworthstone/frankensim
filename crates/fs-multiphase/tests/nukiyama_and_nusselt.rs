//! G0/G1/G2/G3 verification oracles and planted-negative boundary tests for `fs-multiphase`
//! (`KM2.1`: full Nukiyama boiling curve with Schluender/Thome binary suppression +
//! Lienhard-Dhir CHF + `C^1` transition + Bromley/Leidenfrost + Nusselt/Colburn-Hougen condensation).

use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_math::det::{pow, sin};
use fs_multiphase::{
    BoilingCondensationInterfaceV1, BoilingCurveV1, BoilingErrorV1, BoilingFluidKindV1,
    BoilingRegimeV1, BoilingSurfaceCardV1, CondensationErrorV1, CondensationFluidPropertiesV1,
    CondensationGeometryV1, CondensationModelV1, CondensationRegimeV1,
    LIENHARD_DHIR_FLAT_PLATE_C_CHF_V1, LeidenfrostModeV1, NonCondensableGasV1,
    ROHSENOW_SURFACE_LAW_ID_V1, STANDARD_GRAVITY_M_S2, SurfacePhaseChangeModeV1,
    ZUBER_INFINITE_PLATE_C_CHF_V1,
};
use fs_qty::{Dims, Pressure, Temperature};
use std::collections::BTreeMap;

#[test]
fn g2_incropera_1atm_water_nukiyama_curve_csf_shifts_and_denver_altitude() {
    let p_1atm = Pressure::new(101_325.0);
    let ss304 = BoilingSurfaceCardV1::stainless_steel_304_v1();
    let cu_scored = BoilingSurfaceCardV1::scored_copper_v1();
    let cast_iron = BoilingSurfaceCardV1::seasoned_cast_iron_v1();

    let curve_ss = BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, ss304, p_1atm).unwrap();
    let curve_cu =
        BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, cu_scored, p_1atm).unwrap();
    let curve_ci =
        BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, cast_iron, p_1atm).unwrap();
    let curve_zuber = BoilingCurveV1::prepare(
        BoilingFluidKindV1::PureWater,
        ss304.with_chf_constant(ZUBER_INFINITE_PLATE_C_CHF_V1),
        p_1atm,
    )
    .unwrap();

    // 1. Lienhard-Dhir C_CHF = 0.149 reproduces Incropera Eq. 10.6 (~1.261 MW/m^2 +- 5%),
    //    while Zuber 0.131 reproduces ~1.109 MW/m^2.
    let b_ss = curve_ss.regime_boundaries();
    assert!(
        (b_ss.q_chf_w_m2 - 1.261e6).abs() / 1.261e6 < 0.05,
        "1 atm water Lienhard-Dhir CHF expected ~1.261 MW/m^2, got {:.4} MW/m^2",
        b_ss.q_chf_w_m2 / 1e6
    );
    let b_zuber = curve_zuber.regime_boundaries();
    assert!(
        (b_zuber.q_chf_w_m2 - 1.109e6).abs() / 1.109e6 < 0.05,
        "1 atm water Zuber CHF expected ~1.109 MW/m^2, got {:.4} MW/m^2",
        b_zuber.q_chf_w_m2 / 1e6
    );

    // 2. Saturated pool Leidenfrost minimum lies in 195-220 C on SS304 / cast iron,
    //    while subcooled impinging dynamic skillet Leidenfrost shifts upward into 250-335 C.
    let t_min_ss_c = b_ss.t_leidenfrost_k - 273.15;
    let t_min_ci_c = curve_ci.regime_boundaries().t_leidenfrost_k - 273.15;
    assert!(
        (194.0..=222.0).contains(&t_min_ss_c),
        "SS304 saturated pool T_min expected ~195-220 C, got {t_min_ss_c:.2} C"
    );
    assert!(
        (194.0..=222.0).contains(&t_min_ci_c),
        "cast iron saturated pool T_min expected ~195-220 C, got {t_min_ci_c:.2} C"
    );

    let dynamic_skillet = cast_iron.with_leidenfrost_mode(LeidenfrostModeV1::SubcooledImpinging {
        subcooling_k: 75.0,
        dynamic_multiplier: 1.05,
    });
    let curve_dyn =
        BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, dynamic_skillet, p_1atm).unwrap();
    let t_min_dyn_c = curve_dyn.regime_boundaries().t_leidenfrost_k - 273.15;
    assert!(
        (250.0..=335.0).contains(&t_min_dyn_c) && t_min_dyn_c > t_min_ci_c + 45.0,
        "dynamic subcooled skillet Leidenfrost expected 250-335 C, got {t_min_dyn_c:.2} C"
    );

    // 3. Rohsenow C_sf shift: scored copper (0.0068) yields ~(0.013/0.0068)^3 ≈ 7x higher
    //    nucleate heat flux than SS304 (0.013) at the same wall superheat (dT_e = 8 K).
    let t_sat = curve_ss.saturation_temperature().value();
    let q_ss = curve_ss
        .evaluate(Temperature::new(t_sat + 8.0))
        .unwrap()
        .heat_flux_w_m2;
    let q_cu = curve_cu
        .evaluate(Temperature::new(t_sat + 8.0))
        .unwrap()
        .heat_flux_w_m2;
    let expected_ratio = pow(ss304.c_sf / cu_scored.c_sf, 3.0);
    assert!(
        ((q_cu / q_ss) - expected_ratio).abs() / expected_ratio < 0.02,
        "Cu/SS304 Rohsenow ratio expected {expected_ratio:.3}, got {:.3}",
        q_cu / q_ss
    );

    // 4. JANAF-298.15K-formation datum consistency: h_v_janaf - h_l_janaf == h_fg identically.
    let ev_nuc = curve_ss.evaluate(Temperature::new(t_sat + 10.0)).unwrap();
    let h_fg = curve_ss.fluid_properties().latent_heat_j_kg;
    assert!(
        ((ev_nuc.vapor_enthalpy_janaf_j_kg - ev_nuc.liquid_enthalpy_janaf_j_kg) - h_fg).abs()
            < 1e-6,
        "JANAF enthalpy gap must equal IAPWS latent heat h_fg"
    );
    assert!(
        ev_nuc.liquid_enthalpy_janaf_j_kg < -15.0e6,
        "liquid water enthalpy must carry negative JANAF formation datum (~-15.5 MJ/kg)"
    );

    // 5. High-altitude Denver (5,280 ft -> P = 83.5 kPa): T_sat drops to ~94.6 C and CHF drops.
    let curve_denver = BoilingCurveV1::prepare(
        BoilingFluidKindV1::PureWater,
        ss304,
        Pressure::new(83_500.0),
    )
    .unwrap();
    let t_denver_c = curve_denver.saturation_temperature().value() - 273.15;
    assert!(
        (t_denver_c - 94.64).abs() < 0.25,
        "Denver 83.5 kPa T_sat expected ~94.64 C, got {t_denver_c:.3} C"
    );
    assert!(
        curve_denver.regime_boundaries().q_chf_w_m2 < b_ss.q_chf_w_m2,
        "lower ambient pressure at Denver altitude must reduce CHF"
    );
}

#[test]
fn g2_wine_deglazing_schlunder_thome_suppression_and_azeotrope_zero_suppression() {
    let p_1atm = Pressure::new(101_325.0);
    let ss304 = BoilingSurfaceCardV1::stainless_steel_304_v1();

    // 13% ABV table wine corresponds to ~4.7 mol% ethanol (x_EtOH = 0.047).
    let wine = BoilingCurveV1::prepare(
        BoilingFluidKindV1::EthanolWaterNrtl {
            ethanol_mole_fraction: 0.047,
        },
        ss304,
        p_1atm,
    )
    .unwrap();
    let dt_bp_wine = wine.fluid_properties().boiling_range_k;
    assert!(
        dt_bp_wine > 5.0,
        "wine NRTL boiling range Delta T_bp must be > 5 K, got {dt_bp_wine:.3} K"
    );

    // Evaluate well inside nucleate boiling (dT_e = 8 K and 10 K): Schluender/Thome-Shakir must
    // suppress nucleate HTC by 30-60% relative to ideal pure-fluid Rohsenow.
    let t_bub = wine.saturation_temperature().value();
    for &dt_test in &[8.0, 10.0] {
        let ev_wine = wine.evaluate(Temperature::new(t_bub + dt_test)).unwrap();
        assert_eq!(ev_wine.regime, BoilingRegimeV1::NucleateBoiling);
        let suppression_pct = (1.0 - ev_wine.mixture_suppression_factor) * 100.0;
        assert!(
            (30.0..=60.0).contains(&suppression_pct),
            "wine deglazing Schluender-Thome suppression at dT={dt_test} K expected 30-60%, got {suppression_pct:.2}% \
             (T_bub={t_bub:.2} K, Delta T_bp={dt_bp_wine:.2} K, q_mix={:.1} kW/m^2, q_id={:.1} kW/m^2)",
            ev_wine.heat_flux_w_m2 / 1e3,
            ev_wine.ideal_rohsenow_flux_w_m2 / 1e3
        );
    }

    // Ethanol lowers surface tension, shifting saturated Leidenfrost downward vs pure water.
    let pure_water = BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, ss304, p_1atm).unwrap();
    assert!(
        wine.fluid_properties().surface_tension_n_m
            < pure_water.fluid_properties().surface_tension_n_m,
        "wine ethanol must reduce surface tension below pure water"
    );
    assert!(
        wine.regime_boundaries().t_leidenfrost_k < pure_water.regime_boundaries().t_leidenfrost_k,
        "wine Leidenfrost temperature must shift downward below pure water"
    );

    // Near the ethanol-water azeotrope (x_EtOH = 0.897), Delta T_bp -> 0 and suppression vanishes.
    let azeotrope = BoilingCurveV1::prepare(
        BoilingFluidKindV1::EthanolWaterNrtl {
            ethanol_mole_fraction: 0.897,
        },
        ss304,
        p_1atm,
    )
    .unwrap();
    assert!(
        azeotrope.fluid_properties().boiling_range_k < 0.05,
        "azeotrope boiling range must vanish, got {:.4} K",
        azeotrope.fluid_properties().boiling_range_k
    );
    let ev_az = azeotrope
        .evaluate(Temperature::new(
            azeotrope.saturation_temperature().value() + 10.0,
        ))
        .unwrap();
    assert!(
        ev_az.mixture_suppression_factor > 0.99,
        "azeotrope mixture suppression factor must approach 1.0, got {:.5}",
        ev_az.mixture_suppression_factor
    );
}

#[test]
fn g0_nukiyama_c1_continuity_and_exact_tangent_across_all_five_regimes() {
    let curve = BoilingCurveV1::prepare(
        BoilingFluidKindV1::PureWater,
        BoilingSurfaceCardV1::stainless_steel_304_v1(),
        Pressure::new(101_325.0),
    )
    .unwrap();
    let b = curve.regime_boundaries();

    // Check C^1 continuity right across every regime boundary.
    // At T_sat, q_nc = C_nc * dT^(4/3) has exact lim_{dT->0+} q = 0 and dq/dT = 0.
    let eps_sat = 1e-8;
    let sat_left = curve
        .evaluate(Temperature::new(b.t_sat_k - eps_sat))
        .unwrap();
    let sat_right = curve
        .evaluate(Temperature::new(b.t_sat_k + eps_sat))
        .unwrap();
    assert_eq!(sat_left.heat_flux_w_m2, 0.0);
    assert_eq!(sat_left.dq_dtw_w_m2_k, 0.0);
    assert!(sat_right.heat_flux_w_m2 < 1e-6 && sat_right.dq_dtw_w_m2_k < 2.0);

    for &knot in &[b.t_onb_k, b.t_chf_k, b.t_leidenfrost_k] {
        let eps = 1e-5;
        let left = curve.evaluate(Temperature::new(knot - eps)).unwrap();
        let right = curve.evaluate(Temperature::new(knot + eps)).unwrap();
        let q_jump = (right.heat_flux_w_m2 - left.heat_flux_w_m2).abs();
        let dq_jump = (right.dq_dtw_w_m2_k - left.dq_dtw_w_m2_k).abs();
        assert!(
            q_jump <= 1e-3 * right.heat_flux_w_m2.max(100.0) + 2.0,
            "heat flux discontinuity at knot {knot:.3} K: left={}, right={}",
            left.heat_flux_w_m2,
            right.heat_flux_w_m2
        );
        // Slope must be continuous (at CHF and Leidenfrost both sides approach 0).
        if (knot - b.t_chf_k).abs() < 1e-9 || (knot - b.t_leidenfrost_k).abs() < 1e-9 {
            assert!(
                left.dq_dtw_w_m2_k.abs() < 50.0 && right.dq_dtw_w_m2_k.abs() < 50.0,
                "zero-slope extremum violated at {knot:.3} K: left={}, right={}",
                left.dq_dtw_w_m2_k,
                right.dq_dtw_w_m2_k
            );
        } else {
            assert!(
                dq_jump / right.dq_dtw_w_m2_k.abs().max(1.0) < 1e-2,
                "slope jump at knot {knot:.3} K: left={}, right={}",
                left.dq_dtw_w_m2_k,
                right.dq_dtw_w_m2_k
            );
        }
    }

    // Verify exact dq''/dT_w vs central finite difference (< 1e-8) at 14 interior points spanning all 5 regimes.
    let sweep = [
        b.t_sat_k - 5.0,
        b.t_sat_k + 0.5,
        b.t_sat_k + 2.0,
        b.t_onb_k,
        b.t_onb_k + 4.0,
        b.t_chf_k - 3.0,
        b.t_chf_k - 0.5,
        b.t_chf_k + 5.0,
        0.5 * (b.t_chf_k + b.t_leidenfrost_k),
        b.t_leidenfrost_k - 5.0,
        b.t_leidenfrost_k + 3.0,
        b.t_leidenfrost_k + 12.0,
        b.t_leidenfrost_k + 60.0,
        b.t_leidenfrost_k + 180.0,
    ];
    for &tw in &sweep {
        let ev = curve.evaluate(Temperature::new(tw)).unwrap();
        let h = 1e-4;
        let qp = curve
            .evaluate(Temperature::new(tw + h))
            .unwrap()
            .heat_flux_w_m2;
        let qm = curve
            .evaluate(Temperature::new(tw - h))
            .unwrap()
            .heat_flux_w_m2;
        let fd = (qp - qm) / (2.0 * h);
        let rel = (ev.dq_dtw_w_m2_k - fd).abs() / fd.abs().max(1.0);
        assert!(
            rel < 1e-8,
            "FD tangent failed (T_w={tw:.3} K, T_sat={:.3} K, regime={:?}, q_flux={:.2}, dq_dT={}, fd={fd}, rel={rel:e})",
            b.t_sat_k,
            ev.regime,
            ev.heat_flux_w_m2,
            ev.dq_dtw_w_m2_k
        );
    }
}

#[test]
fn g1_coupled_pan_quench_sweeps_leidenfrost_transition_chf_nucleate_with_exact_conservation() {
    // Lumped skillet wall (C_wall = 2000 J/K, wetted area A = 0.04 m^2) starting in film boiling
    // at 280 C (553.15 K) quenched by 1 atm saturated water down to 105 C (378.15 K).
    // Integrates C_wall * dT_w/dt = -A * q''(T_w) via implicit midpoint with safeguarded damped
    // Newton + curve.recommend_step_cap(), verifying zero 2-cycling and < 1e-11 energy defect.
    let curve = BoilingCurveV1::prepare(
        BoilingFluidKindV1::PureWater,
        BoilingSurfaceCardV1::stainless_steel_304_v1(),
        Pressure::new(101_325.0),
    )
    .unwrap();

    let c_wall = 2_000.0_f64;
    let area = 0.04_f64;
    let t_init = 553.15_f64;
    let mut t_w = t_init;
    let mut q_cumulative_j = 0.0_f64;
    let mut saw_film = false;
    let mut saw_transition = false;
    let mut saw_nucleate = false;
    let mut total_steps = 0_usize;

    while t_w > 378.15 && total_steps < 2_000 {
        let mut dt = 0.25_f64;
        loop {
            // Solve implicit midpoint: R(T1) = C*(T1 - T0) + dt*A*q''(0.5*(T0+T1)) = 0
            let mut t1 = t_w;
            let mut converged = false;
            for _iter in 0..40 {
                let t_mid = 0.5 * (t_w + t1);
                let ev = curve.evaluate(Temperature::new(t_mid)).unwrap();
                let res = c_wall * (t1 - t_w) + dt * area * ev.heat_flux_w_m2;
                if res.abs() <= 1e-11 * c_wall {
                    converged = true;
                    break;
                }
                let jac = c_wall + 0.5 * dt * area * ev.dq_dtw_w_m2_k;
                if jac <= 0.1 * c_wall {
                    // Negative slope in transition boiling exceeded step stiffness: shorten dt.
                    break;
                }
                let step = (res / jac).clamp(-20.0, 20.0);
                t1 -= step;
            }
            if !converged {
                dt *= 0.5;
                assert!(
                    dt > 1e-6,
                    "damped Newton failed to converge even at tiny dt"
                );
                continue;
            }
            let advice = curve
                .recommend_step_cap(Temperature::new(t_w), Temperature::new(t1), dt, 6.0)
                .unwrap();
            if advice.shorten_recommended && dt > 0.01 {
                dt = advice.recommended_max_dt_s;
                continue;
            }
            let t_mid = 0.5 * (t_w + t1);
            let ev_mid = curve.evaluate(Temperature::new(t_mid)).unwrap();
            match ev_mid.regime {
                BoilingRegimeV1::FilmBoiling => saw_film = true,
                BoilingRegimeV1::TransitionBoiling => saw_transition = true,
                BoilingRegimeV1::NucleateBoiling => saw_nucleate = true,
                _ => {}
            }
            q_cumulative_j += dt * area * ev_mid.heat_flux_w_m2;
            t_w = t1;
            total_steps += 1;
            break;
        }
    }

    assert!(
        saw_film && saw_transition && saw_nucleate,
        "quench must traverse Film -> Transition -> Nucleate"
    );
    assert!(t_w <= 378.15, "skillet must quench below 105 C");
    let stored_loss = c_wall * (t_init - t_w);
    let rel_energy_defect = (stored_loss - q_cumulative_j).abs() / stored_loss;
    assert!(
        rel_energy_defect < 1e-11,
        "pan-quench discrete energy defect must be < 1e-11, got {rel_energy_defect:e}"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn g0_g2_nusselt_condensation_oracles_and_colburn_hougen_air_degradation() {
    // 1. Exact analytical vertical-plate, inclined-plate, and horizontal-cylinder laminar Nusselt.
    let props = CondensationFluidPropertiesV1 {
        liquid_density_kg_m3: 960.0,
        vapor_density_kg_m3: 0.60,
        liquid_viscosity_pa_s: 2.8e-4,
        liquid_conductivity_w_m_k: 0.68,
        liquid_cp_j_kg_k: 4_216.0,
        latent_heat_j_kg: 2.257e6,
        vapor_enthalpy_janaf_j_kg: -13.29e6,
        liquid_enthalpy_janaf_j_kg: -15.55e6,
    };
    let t_dew = 373.15;
    let tw = 368.15;
    let dt = t_dew - tw;
    let h_fg_p = props.latent_heat_j_kg + 0.68 * props.liquid_cp_j_kg_k * dt;

    let vert = CondensationModelV1::fixed_properties(
        Pressure::new(101_325.0),
        Temperature::new(t_dew),
        CondensationGeometryV1::VerticalPlate { height_m: 0.04 },
        props,
        false,
    )
    .unwrap();
    let ev_vert = vert.evaluate(Temperature::new(tw)).unwrap();
    let expected_vert = 0.943
        * pow(
            (STANDARD_GRAVITY_M_S2
                * props.liquid_density_kg_m3
                * (props.liquid_density_kg_m3 - props.vapor_density_kg_m3)
                * pow(props.liquid_conductivity_w_m_k, 3.0)
                * h_fg_p)
                / (props.liquid_viscosity_pa_s * dt * 0.04),
            0.25,
        );
    assert!(
        (ev_vert.film_htc_w_m2_k - expected_vert).abs() / expected_vert < 1e-10,
        "vertical plate Nusselt analytical oracle error: calc={}, ref={expected_vert}",
        ev_vert.film_htc_w_m2_k
    );

    // Inclined plate at 30 deg: scales by sin(pi/6)^(1/4) = 0.5^(1/4).
    let tilt = core::f64::consts::FRAC_PI_6;
    let incl = CondensationModelV1::fixed_properties(
        Pressure::new(101_325.0),
        Temperature::new(t_dew),
        CondensationGeometryV1::InclinedPlate {
            length_m: 0.04,
            tilt_from_horizontal_rad: tilt,
        },
        props,
        false,
    )
    .unwrap();
    let ev_incl = incl.evaluate(Temperature::new(tw)).unwrap();
    let expected_incl = expected_vert * pow(sin(tilt), 0.25);
    assert!((ev_incl.film_htc_w_m2_k - expected_incl).abs() / expected_incl < 1e-10);

    // Horizontal cylinder: exact ratio 0.729 / 0.943 for equal D = L.
    let cyl = CondensationModelV1::fixed_properties(
        Pressure::new(101_325.0),
        Temperature::new(t_dew),
        CondensationGeometryV1::HorizontalCylinder { diameter_m: 0.04 },
        props,
        false,
    )
    .unwrap();
    let ev_cyl = cyl.evaluate(Temperature::new(tw)).unwrap();
    let expected_cyl = expected_vert * (0.729 / 0.943);
    assert!((ev_cyl.film_htc_w_m2_k - expected_cyl).abs() / expected_cyl < 1e-10);

    // 2. Combi-oven 100% pure steam vs 30% steam (70% air) Colburn-Hougen degradation.
    let p_1atm = Pressure::new(101_325.0);
    let geom = CondensationGeometryV1::VerticalPlate { height_m: 0.15 };
    let tw_food = Temperature::new(313.15); // 40 C dough/protein surface

    let pure_steam = CondensationModelV1::pure_steam(p_1atm, geom).unwrap();
    let ev_100 = pure_steam.evaluate(tw_food).unwrap();
    assert!(
        (3_000.0..=15_000.0).contains(&ev_100.effective_htc_w_m2_k),
        "100% steam condensation HTC expected 3,000-15,000 W/(m^2 K), got {:.1}",
        ev_100.effective_htc_w_m2_k
    );

    let mut prev_q = ev_100.heat_flux_w_m2;
    let mut prev_h = ev_100.effective_htc_w_m2_k;
    for w_air in [0.02, 0.10, 0.30, 0.50, 0.70] {
        let gas = NonCondensableGasV1::steam_air(w_air, Temperature::new(373.15), 25.0).unwrap();
        let model = CondensationModelV1::steam_air(p_1atm, geom, gas).unwrap();
        let ev = model.evaluate(tw_food).unwrap();
        assert!(
            ev.heat_flux_w_m2 < prev_q && ev.effective_htc_w_m2_k < prev_h,
            "air fraction {w_air} must monotonically degrade condensation flux and HTC"
        );
        assert!(
            ev.interface_temperature.value() < ev.dew_temperature.value(),
            "non-condensable air must depress T_i below T_dew"
        );
        prev_q = ev.heat_flux_w_m2;
        prev_h = ev.effective_htc_w_m2_k;
    }
    // 70% air (30% steam combi mode) cuts HTC dramatically (> 5x below 100% pure steam).
    assert!(
        prev_h < 0.2 * ev_100.effective_htc_w_m2_k,
        "30% steam combi mode HTC ({prev_h:.1}) must drop >5x below 100% pure steam ({:.1})",
        ev_100.effective_htc_w_m2_k
    );

    // 3. Handoff guard: T_w >= T_dew produces zero condensation.
    let dry_wall = pure_steam.evaluate(Temperature::new(375.0)).unwrap();
    assert_eq!(dry_wall.regime, CondensationRegimeV1::NoCondensation);
    assert_eq!(dry_wall.heat_flux_w_m2, 0.0);

    // 4. Unified BoilingCondensationInterfaceV1 single-source ownership + matdb card load.
    let mut params = BTreeMap::new();
    params.insert(
        "c_sf".into(),
        LawParameter {
            value: 0.0130,
            dims: Dims::NONE,
        },
    );
    params.insert(
        "c_chf".into(),
        LawParameter {
            value: LIENHARD_DHIR_FLAT_PLATE_C_CHF_V1,
            dims: Dims::NONE,
        },
    );
    params.insert(
        "contact_angle_rad".into(),
        LawParameter {
            value: 1.15,
            dims: Dims::NONE,
        },
    );
    params.insert(
        "emissivity".into(),
        LawParameter {
            value: 0.20,
            dims: Dims::NONE,
        },
    );
    params.insert(
        "characteristic_length_m".into(),
        LawParameter {
            value: 0.25,
            dims: Dims([1, 0, 0, 0, 0, 0]),
        },
    );
    params.insert(
        "wall_k".into(),
        LawParameter {
            value: 16.2,
            dims: Dims([1, 1, -3, -1, 0, 0]),
        },
    );
    params.insert(
        "wall_rho".into(),
        LawParameter {
            value: 7_900.0,
            dims: Dims([-3, 1, 0, 0, 0, 0]),
        },
    );
    params.insert(
        "wall_cp".into(),
        LawParameter {
            value: 500.0,
            dims: Dims([2, 0, -2, -1, 0, 0]),
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(ROHSENOW_SURFACE_LAW_ID_V1.into()),
        law_version: 1,
        parameters: params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 273.15, 900.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Incropera Table 10.1 SS304".into(),
            license: "test".into(),
            artifact: None,
        },
    };
    let loaded_surf = BoilingSurfaceCardV1::from_card(&card).unwrap();
    let boiling =
        BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, loaded_surf, p_1atm).unwrap();
    let gas_50 = NonCondensableGasV1::steam_air(0.50, Temperature::new(373.15), 25.0).unwrap();
    let cond_50 = CondensationModelV1::steam_air(p_1atm, geom, gas_50).unwrap();
    let iface = BoilingCondensationInterfaceV1::new(boiling, cond_50).unwrap();

    // T_w > T_sat -> Boiling (> 0 flux)
    let ev_hot = iface.evaluate(Temperature::new(385.0)).unwrap();
    assert!(matches!(ev_hot.mode, SurfacePhaseChangeModeV1::Boiling(_)));
    assert!(ev_hot.signed_heat_flux_w_m2 > 0.0);

    // T_dew < T_w < T_sat -> NoPhaseChange (zero flux, handed off to fs-convection)
    let t_mid = 0.5
        * (iface.condensation().dew_temperature().value()
            + iface.boiling().saturation_temperature().value());
    let ev_mid = iface.evaluate(Temperature::new(t_mid)).unwrap();
    assert_eq!(ev_mid.mode, SurfacePhaseChangeModeV1::NoPhaseChange);
    assert_eq!(ev_mid.signed_heat_flux_w_m2, 0.0);

    // T_w < T_dew -> Condensation (< 0 wall-to-fluid flux)
    let ev_cold = iface.evaluate(Temperature::new(320.0)).unwrap();
    assert!(matches!(
        ev_cold.mode,
        SurfacePhaseChangeModeV1::Condensation(_)
    ));
    assert!(ev_cold.signed_heat_flux_w_m2 < 0.0);
    assert!(ev_cold.dq_dtw_w_m2_k > 0.0);
}

#[test]
#[allow(clippy::too_many_lines)]
fn planted_negatives_return_typed_errors_with_zero_panics_or_nans() {
    let ss304 = BoilingSurfaceCardV1::stainless_steel_304_v1();
    let p_1atm = Pressure::new(101_325.0);

    // 1. Round-2 Physics #3: pure triglyceride cooking oils MUST refuse Nukiyama boiling.
    for (name, smoke_k, flash_k) in [
        ("canola_oil", 493.15, 588.15),
        ("peanut_oil", 505.15, 588.15),
        ("extra_virgin_olive_oil", 463.15, 583.15),
    ] {
        assert!(matches!(
            BoilingCurveV1::prepare(
                BoilingFluidKindV1::TriglycerideOil {
                    name,
                    smoke_point_k: smoke_k,
                    flash_point_k: flash_k,
                },
                ss304,
                p_1atm,
            ),
            Err(BoilingErrorV1::NonBoilingLipidPhase { .. })
        ));
    }

    // 2. Invalid pressures (<= 0, NaN, inf) -> typed InvalidPressure.
    for bad_p in [0.0, -101_325.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let p = Pressure::new(bad_p);
        assert!(matches!(
            BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, ss304, p),
            Err(BoilingErrorV1::InvalidPressure { .. })
        ));
        assert!(matches!(
            CondensationModelV1::pure_steam(
                p,
                CondensationGeometryV1::VerticalPlate { height_m: 0.1 }
            ),
            Err(CondensationErrorV1::InvalidPressure { .. })
        ));
    }

    // 3. Invalid C_sf <= 0, NaN, bad emissivity, bad contact angle -> InvalidParameter.
    let mut bad_csf = ss304;
    bad_csf.c_sf = 0.0;
    assert!(matches!(
        BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, bad_csf, p_1atm),
        Err(BoilingErrorV1::InvalidParameter(_))
    ));
    bad_csf.c_sf = -0.01;
    assert!(matches!(
        BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, bad_csf, p_1atm),
        Err(BoilingErrorV1::InvalidParameter(_))
    ));
    bad_csf.c_sf = f64::NAN;
    assert!(matches!(
        BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, bad_csf, p_1atm),
        Err(BoilingErrorV1::InvalidParameter(_))
    ));

    // 4. Invalid ethanol mole fraction -> InvalidComposition.
    for bad_x in [-0.05, 1.05, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            BoilingCurveV1::prepare(
                BoilingFluidKindV1::EthanolWaterNrtl {
                    ethanol_mole_fraction: bad_x,
                },
                ss304,
                p_1atm,
            ),
            Err(BoilingErrorV1::InvalidComposition(_))
        ));
    }

    // 5. Invalid wall temperatures (T <= 0, NaN, inf) -> InvalidTemperature.
    let curve = BoilingCurveV1::prepare(BoilingFluidKindV1::PureWater, ss304, p_1atm).unwrap();
    let cond = CondensationModelV1::pure_steam(
        p_1atm,
        CondensationGeometryV1::VerticalPlate { height_m: 0.1 },
    )
    .unwrap();
    for bad_t in [0.0, -50.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let t = Temperature::new(bad_t);
        assert!(matches!(
            curve.evaluate(t),
            Err(BoilingErrorV1::InvalidTemperature { .. })
        ));
        assert!(matches!(
            cond.evaluate(t),
            Err(CondensationErrorV1::InvalidTemperature { .. })
        ));
    }

    // 6. Invalid air mass fraction outside [0, 1] or NaN -> InvalidNonCondensableFraction.
    for bad_w in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            NonCondensableGasV1::steam_air(bad_w, Temperature::new(373.15), 25.0),
            Err(CondensationErrorV1::InvalidNonCondensableFraction { .. })
        ));
    }

    // 7. Invalid condensation geometry -> InvalidGeometry.
    for bad_geom in [
        CondensationGeometryV1::VerticalPlate { height_m: 0.0 },
        CondensationGeometryV1::HorizontalCylinder { diameter_m: -0.1 },
        CondensationGeometryV1::InclinedPlate {
            length_m: 0.1,
            tilt_from_horizontal_rad: 0.0,
        },
        CondensationGeometryV1::InclinedPlate {
            length_m: 0.1,
            tilt_from_horizontal_rad: 2.0,
        },
    ] {
        assert!(matches!(
            CondensationModelV1::pure_steam(p_1atm, bad_geom),
            Err(CondensationErrorV1::InvalidGeometry(_))
        ));
    }
}
