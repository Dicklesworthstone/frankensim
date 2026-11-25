//! G0/G2/G3 integration oracles and planted-negative boundary tests for KM0.1:
//! asymmetric NRTL/UNIQUAC excess-Gibbs VLE and GAB/Oswin moisture sorption isotherms.

use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_math::det::ln;
use fs_qty::{Dims, Pressure, Temperature};
use fs_thermochem::caloric::sorption::{
    GAB_ISOTHERM_LAW_ID_V1, GabIsothermV1, OSWIN_ISOTHERM_LAW_ID_V1, OswinIsothermV1,
    SorptionErrorV1, SorptionIsothermV1, WATER_MOLAR_MASS_KG_PER_MOL_V1,
};
use fs_thermochem::solution::nrtl::{NRTL_LAW_ID_V1, NrtlModelV1};
use fs_thermochem::solution::uniquac::{UNIQUAC_LAW_ID_V1, UniquacComponentV1, UniquacModelV1};
use fs_thermochem::solution::{
    ExcessGibbsErrorV1, FlashPhaseRegimeV1, ModifiedRaoultVleV1, VaporPressureLawV1,
    VleSolverOptionsV1,
};
use fs_thermochem::water::SaturatedWaterState;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

/// Deterministic splitmix64 PRNG for simplex sampling without external dependencies.
struct DetRng(u64);
impl DetRng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    fn next_unit(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0 / ((1_u64 << 53) as f64))
    }
    fn simplex(&mut self, n: usize) -> Vec<f64> {
        let mut raw: Vec<f64> = (0..n).map(|_| -ln((self.next_unit()).max(1e-12))).collect();
        let sum: f64 = raw.iter().sum();
        for v in &mut raw {
            *v /= sum;
        }
        raw
    }
}

#[test]
fn g0_nrtl_and_uniquac_gibbs_duhem_identity_holds_across_simplex_and_temperatures() {
    // 4-component asymmetric polar mixture for NRTL and UNIQUAC.
    let species = vec![
        SpeciesId::new("C2H5OH").unwrap(),
        SpeciesId::new("H2O").unwrap(),
        SpeciesId::new("CH3COOH").unwrap(),
        SpeciesId::new("C4H8O2").unwrap(),
    ];
    let n = species.len();
    let nrtl_a = vec![
        0.0, -0.8009, 0.22, -0.35, 3.4578, 0.0, -0.45, 1.15, -0.18, 0.62, 0.0, 0.14, 0.48, -0.92,
        -0.27, 0.0,
    ];
    let nrtl_b = vec![
        0.0, 246.18, -110.0, 175.0, -586.0809, 0.0, 320.0, -290.0, 165.0, -210.0, 0.0, 95.0,
        -135.0, 410.0, 125.0, 0.0,
    ];
    let nrtl_alpha = vec![
        0.0, 0.30, 0.25, 0.32, 0.30, 0.0, 0.35, 0.28, 0.25, 0.35, 0.0, 0.20, 0.32, 0.28, 0.20, 0.0,
    ];
    let nrtl = NrtlModelV1::new(
        species.clone(),
        nrtl_a,
        nrtl_b,
        nrtl_alpha,
        Temperature::new(273.0),
        Temperature::new(450.0),
    )
    .unwrap();

    let uq_comps = vec![
        UniquacComponentV1::modified(species[0].clone(), 2.1055, 1.972, 0.92).unwrap(),
        UniquacComponentV1::modified(species[1].clone(), 0.92, 1.40, 1.00).unwrap(),
        UniquacComponentV1::standard(species[2].clone(), 2.2024, 2.072).unwrap(),
        UniquacComponentV1::standard(species[3].clone(), 3.4786, 3.116).unwrap(),
    ];
    let uq_a = vec![
        0.0, -1.20, 0.15, -0.40, 1.45, 0.0, -0.30, 0.85, -0.12, 0.38, 0.0, 0.10, 0.28, -0.75,
        -0.19, 0.0,
    ];
    let uq_b = vec![
        0.0, 420.0, -85.0, 130.0, -510.0, 0.0, 210.0, -260.0, 95.0, -140.0, 0.0, 70.0, -105.0,
        310.0, 88.0, 0.0,
    ];
    let uniquac = UniquacModelV1::new(
        uq_comps,
        uq_a,
        uq_b,
        Temperature::new(273.0),
        Temperature::new(450.0),
    )
    .unwrap();

    let mut rng = DetRng(0x5eed_cafe_beef_0001);
    let temps = [273.15, 298.15, 333.15, 351.30, 373.15, 420.0, 450.0];

    for &t_k in &temps {
        let t = Temperature::new(t_k);
        // Include pure vertices, binary edges, and 16 random interior simplex points.
        let mut points: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.5, 0.5, 0.0, 0.0],
            vec![0.25, 0.25, 0.25, 0.25],
            vec![1e-8, 0.99999998, 1e-8, 0.0],
        ];
        for _ in 0..16 {
            points.push(rng.simplex(n));
        }

        for x in &points {
            let nrtl_ev = nrtl.evaluate(t, x).unwrap();
            let uq_ev = uniquac.evaluate(t, x).unwrap();

            for (label, ev) in [("NRTL", &nrtl_ev), ("UNIQUAC", &uq_ev)] {
                // 1. Euler identity: sum_i x_i ln(gamma_i) == g^E/(RT).
                let euler: f64 = ev
                    .mole_fractions()
                    .iter()
                    .zip(ev.ln_activity_coefficients())
                    .map(|(xi, lg)| xi * lg)
                    .sum();
                let euler_res = (euler - ev.reduced_excess_gibbs()).abs();
                assert!(
                    euler_res < 1e-12,
                    "{label} Euler mismatch at T={t_k}, x={x:?}: residual={euler_res:e}"
                );

                // 2. Gibbs-Helmholtz Euler identity: sum_i x_i h_i^E == h^E.
                let h_euler: f64 = ev
                    .mole_fractions()
                    .iter()
                    .zip(ev.partial_excess_enthalpy_j_per_mol())
                    .map(|(xi, hi)| xi * hi)
                    .sum();
                let h_res = (h_euler - ev.excess_enthalpy_j_per_mol()).abs();
                assert!(
                    h_res < 1e-9,
                    "{label} enthalpy Euler mismatch at T={t_k}, x={x:?}: residual={h_res:e}"
                );

                // 3. Column Gibbs-Duhem residuals: |sum_i x_i d(ln gamma_i)/dx_j| < 1e-12.
                for (j, col_res) in ev.gibbs_duhem_column_residuals().into_iter().enumerate() {
                    assert!(
                        col_res.abs() < 1e-12,
                        "{label} Gibbs-Duhem column {j} failed at T={t_k}, x={x:?}: residual={col_res:e}"
                    );
                }

                // 4. Directional simplex Gibbs-Duhem residual < 1e-12.
                let dir = [0.6, -0.25, -0.20, -0.15];
                let gd_dir = ev.gibbs_duhem_residual(&dir).unwrap();
                assert!(
                    gd_dir.abs() < 1e-12,
                    "{label} directional Gibbs-Duhem failed at T={t_k}, x={x:?}: residual={gd_dir:e}"
                );

                // 5. Jacobian symmetry d(ln gamma_i)/dx_j == d(ln gamma_j)/dx_i.
                for i in 0..n {
                    for j in (i + 1)..n {
                        let a_ij = ev.d_ln_gamma_dx_entry(i, j).unwrap();
                        let a_ji = ev.d_ln_gamma_dx_entry(j, i).unwrap();
                        assert!(
                            (a_ij - a_ji).abs() < 1e-12,
                            "{label} Jacobian asymmetry ({i},{j}) at T={t_k}: {a_ij} vs {a_ji}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g2_ethanol_water_1atm_azeotrope_txy_table_and_denver_altitude_drop() {
    let vle = NrtlModelV1::ethanol_water_vle_v1();
    let opts = VleSolverOptionsV1::default();
    let p_1atm = Pressure::new(101_325.0);

    // 1. Locate exact 1 atm azeotrope where y_EtOH(x) == x on [0.80, 0.99].
    let mut lo = 0.80_f64;
    let mut hi = 0.99_f64;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let bp = vle.bubble_point_t(p_1atm, &[mid, 1.0 - mid], opts).unwrap();
        if bp.vapor_y[0] > mid {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let x_az = 0.5 * (lo + hi);
    let bp_az = vle
        .bubble_point_t(p_1atm, &[x_az, 1.0 - x_az], opts)
        .unwrap();
    let t_az = bp_az.temperature.value();

    assert!(
        (t_az - 351.30).abs() <= 0.15,
        "1 atm EtOH-H2O azeotrope T out of tolerance: T_calc={t_az:.4} K, T_ref=351.30 K, residual={:.4} K",
        t_az - 351.30
    );
    assert!(
        (x_az - 0.894).abs() <= 0.005,
        "1 atm EtOH-H2O azeotrope x out of tolerance: x_calc={x_az:.5}, x_ref=0.894, residual={:.5}",
        x_az - 0.894
    );
    assert!(
        (bp_az.vapor_y[0] - x_az).abs() < 1e-9,
        "azeotrope must satisfy y_EtOH == x_EtOH: y={}, x={x_az}",
        bp_az.vapor_y[0]
    );

    // Minimum-boiling azeotrope must boil strictly below BOTH pure ethanol and pure water.
    let bp_pure_etoh = vle.bubble_point_t(p_1atm, &[1.0, 0.0], opts).unwrap();
    let bp_pure_h2o = vle.bubble_point_t(p_1atm, &[0.0, 1.0], opts).unwrap();
    assert!(
        t_az < bp_pure_etoh.temperature.value() && t_az < bp_pure_h2o.temperature.value(),
        "azeotrope T={t_az} K must lie below pure EtOH ({} K) and pure H2O ({} K)",
        bp_pure_etoh.temperature.value(),
        bp_pure_h2o.temperature.value()
    );

    // 2. Compare 1 atm T-x-y table against DECHEMA/NIST published ethanol-water literature points.
    let lit_table: &[(f64, f64, f64)] = &[
        // (x_EtOH, T_ref_K, y_EtOH_ref)
        (0.000, 373.15, 0.000),
        (0.050, 363.65, 0.332),
        (0.100, 359.65, 0.442),
        (0.200, 356.25, 0.531),
        (0.300, 354.85, 0.576),
        (0.500, 352.95, 0.654),
        (0.700, 351.85, 0.753),
        (0.894, 351.30, 0.894),
        (1.000, 351.45, 1.000),
    ];
    for &(x_etoh, t_ref, y_ref) in lit_table {
        let bp = vle
            .bubble_point_t(p_1atm, &[x_etoh, 1.0 - x_etoh], opts)
            .unwrap();
        let t_err = (bp.temperature.value() - t_ref).abs();
        let y_err = (bp.vapor_y[0] - y_ref).abs();
        assert!(
            t_err <= 0.60,
            "EtOH-H2O 1 atm T mismatch: (T_calc={:.2} K, T_ref={t_ref:.2} K, P=101325 Pa, x={x_etoh}, y_calc={:.4}, y_ref={y_ref:.4}, residual={t_err:.3} K)",
            bp.temperature.value(),
            bp.vapor_y[0]
        );
        assert!(
            y_err <= 0.015,
            "EtOH-H2O 1 atm y mismatch: (T={:.2} K, P=101325 Pa, x={x_etoh}, y_calc={:.4}, y_ref={y_ref:.4}, residual={y_err:.4})",
            bp.temperature.value(),
            bp.vapor_y[0]
        );
    }

    // 3. Boiling range Delta T_bp = T_dew - T_bub for wine deglazing (~5 mol% EtOH) vs azeotrope.
    let (wine_bub, wine_dew, wine_dt_bp) =
        vle.boiling_range_k(p_1atm, &[0.05, 0.95], opts).unwrap();
    assert!(
        wine_dt_bp > 5.0 && wine_dew.temperature.value() > wine_bub.temperature.value(),
        "5 mol% wine must have substantial boiling range Delta T_bp > 5 K, got {wine_dt_bp:.3} K"
    );
    let (_, _, az_dt_bp) = vle
        .boiling_range_k(p_1atm, &[x_az, 1.0 - x_az], opts)
        .unwrap();
    assert!(
        az_dt_bp.abs() < 1e-5,
        "azeotrope boiling range Delta T_bp must collapse to 0, got {az_dt_bp:e} K"
    );

    // 4. High-altitude Denver (5,280 ft -> P = 83.5 kPa): boiling point drops cleanly.
    let p_denver = Pressure::new(83_500.0);
    let denver_h2o = vle.bubble_point_t(p_denver, &[0.0, 1.0], opts).unwrap();
    let denver_az = vle
        .bubble_point_t(p_denver, &[x_az, 1.0 - x_az], opts)
        .unwrap();
    let denver_h2o_c = denver_h2o.temperature.value() - 273.15;
    assert!(
        (denver_h2o_c - 94.64).abs() < 0.25,
        "Denver 83.5 kPa water boiling point expected ~94.64 C, got {denver_h2o_c:.3} C"
    );
    assert!(
        denver_az.temperature.value() < t_az - 4.0,
        "Denver azeotrope must drop > 4 K below 1 atm: denver={} K vs 1atm={t_az} K",
        denver_az.temperature.value()
    );

    // 5. Isothermal PT flash checks mass balance, phase equilibrium, and single-phase regimes.
    let flash_2ph = vle
        .flash_pt(Temperature::new(358.0), p_1atm, &[0.20, 0.80], opts)
        .unwrap();
    assert_eq!(flash_2ph.phase_regime, FlashPhaseRegimeV1::TwoPhase);
    assert!(flash_2ph.vapor_fraction > 0.0 && flash_2ph.vapor_fraction < 1.0);
    assert!(flash_2ph.relative_residual < 1e-10);
    for i in 0..2 {
        let mb = (1.0 - flash_2ph.vapor_fraction) * flash_2ph.liquid_x[i]
            + flash_2ph.vapor_fraction * flash_2ph.vapor_y[i];
        let z_i = [0.20, 0.80][i];
        assert!(
            (mb - z_i).abs() < 1e-11,
            "flash mass balance component {i}: got {mb}, expected {z_i}"
        );
    }
    let flash_sub = vle
        .flash_pt(Temperature::new(320.0), p_1atm, &[0.20, 0.80], opts)
        .unwrap();
    assert_eq!(flash_sub.phase_regime, FlashPhaseRegimeV1::SubcooledLiquid);
    assert_eq!(flash_sub.vapor_fraction, 0.0);

    let flash_sup = vle
        .flash_pt(Temperature::new(385.0), p_1atm, &[0.20, 0.80], opts)
        .unwrap();
    assert_eq!(flash_sup.phase_regime, FlashPhaseRegimeV1::SuperheatedVapor);
    assert_eq!(flash_sup.vapor_fraction, 1.0);
}

#[test]
fn g2_uniquac_ethanol_water_and_wagner_vapor_pressure_oracle() {
    let vle = UniquacModelV1::ethanol_water_vle_v1();
    let opts = VleSolverOptionsV1::default();
    let p_1atm = Pressure::new(101_325.0);

    let mut lo = 0.80_f64;
    let mut hi = 0.99_f64;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let bp = vle.bubble_point_t(p_1atm, &[mid, 1.0 - mid], opts).unwrap();
        if bp.vapor_y[0] > mid {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let x_az = 0.5 * (lo + hi);
    let bp_az = vle
        .bubble_point_t(p_1atm, &[x_az, 1.0 - x_az], opts)
        .unwrap();
    assert!(
        (bp_az.temperature.value() - 351.30).abs() < 0.25,
        "UNIQUAC EtOH-H2O azeotrope T={:.3} K",
        bp_az.temperature.value()
    );
    assert!(
        (0.85..=0.92).contains(&x_az),
        "UNIQUAC EtOH-H2O azeotrope x={x_az:.4} out of physical window"
    );

    // Wagner 2.5-5 propane check at normal boiling point 231.02 K (Poling/Prausnitz/O'Connell 5th ed., Appendix A).
    let wagner_c3h8 = VaporPressureLawV1::wagner_25_5(
        Temperature::new(369.83),
        Pressure::new(4.248e6),
        -6.76368,
        1.55481,
        -1.5872,
        -2.024,
        Temperature::new(150.0),
        Temperature::new(369.83),
    )
    .unwrap();
    let pt_c3 = wagner_c3h8.evaluate(Temperature::new(231.02)).unwrap();
    assert!(
        (pt_c3.pressure_pa - 101_325.0).abs() / 101_325.0 < 2e-3,
        "Wagner 2.5-5 propane NBP P_sat={:.1} Pa vs 101325.0 Pa",
        pt_c3.pressure_pa
    );
    let dt = 1e-3;
    let p_plus = wagner_c3h8.evaluate(Temperature::new(231.02 + dt)).unwrap();
    let p_minus = wagner_c3h8.evaluate(Temperature::new(231.02 - dt)).unwrap();
    let fd_dp = (p_plus.pressure_pa - p_minus.pressure_pa) / (2.0 * dt);
    assert!(
        (pt_c3.dp_dt_pa_per_k - fd_dp).abs() / fd_dp < 1e-8,
        "Wagner 2.5-5 dP/dT analytic={} vs fd={fd_dp}",
        pt_c3.dp_dt_pa_per_k
    );

    // Wagner 3-6 water check against IAPWS at 373.15 K (Poling/Prausnitz/O'Connell 5th ed., Appendix A).
    let wagner_h2o = VaporPressureLawV1::wagner_3_6(
        Temperature::new(647.14),
        Pressure::new(22.064e6),
        -7.77224,
        1.45684,
        -2.71942,
        -1.41336,
        Temperature::new(273.16),
        Temperature::new(647.14),
    )
    .unwrap();
    let pt = wagner_h2o.evaluate(Temperature::new(373.15)).unwrap();
    let iapws = VaporPressureLawV1::iapws_water()
        .evaluate(Temperature::new(373.15))
        .unwrap();
    assert!(
        (pt.pressure_pa - iapws.pressure_pa).abs() / iapws.pressure_pa < 1e-3,
        "Wagner 3-6 water P_sat={:.1} Pa vs IAPWS={:.1} Pa",
        pt.pressure_pa,
        iapws.pressure_pa
    );
    let h2o_p = wagner_h2o.evaluate(Temperature::new(373.15 + dt)).unwrap();
    let h2o_m = wagner_h2o.evaluate(Temperature::new(373.15 - dt)).unwrap();
    let h2o_fd = (h2o_p.pressure_pa - h2o_m.pressure_pa) / (2.0 * dt);
    assert!(
        (pt.dp_dt_pa_per_k - h2o_fd).abs() / h2o_fd < 1e-8,
        "Wagner 3-6 water dP/dT analytic={} vs fd={h2o_fd}",
        pt.dp_dt_pa_per_k
    );

    let wagner_36 = VaporPressureLawV1::wagner_3_6(
        Temperature::new(513.92),
        Pressure::new(6.148e6),
        -8.51838,
        0.34163,
        -5.73683,
        8.32581,
        Temperature::new(273.15),
        Temperature::new(513.92),
    )
    .unwrap();
    let e_pt = wagner_36.evaluate(Temperature::new(351.45)).unwrap();
    let e_p = wagner_36.evaluate(Temperature::new(351.45 + dt)).unwrap();
    let e_m = wagner_36.evaluate(Temperature::new(351.45 - dt)).unwrap();
    let e_fd = (e_p.pressure_pa - e_m.pressure_pa) / (2.0 * dt);
    assert!(
        (e_pt.pressure_pa - 101_325.0).abs() / 101_325.0 < 0.02,
        "Wagner 3-6 ethanol normal boiling pressure={:.1} Pa",
        e_pt.pressure_pa
    );
    assert!(
        (e_pt.dp_dt_pa_per_k - e_fd).abs() / e_fd < 1e-8,
        "Wagner 3-6 dP/dT analytic={} vs fd={e_fd}",
        e_pt.dp_dt_pa_per_k
    );
}

#[test]
fn g2_sucrose_syrup_and_nacl_brine_boiling_point_elevation() {
    let opts = VleSolverOptionsV1::default();
    let p_1atm = Pressure::new(101_325.0);
    let mw_h2o = 18.01528;

    // 1. Water-Sucrose syrup BPE vs Sugar Technologists Manual (Bubnik et al.).
    let sucrose_vle = NrtlModelV1::water_sucrose_vle_v1();
    let t_pure_h2o = sucrose_vle
        .bubble_point_t(p_1atm, &[1.0, 0.0], opts)
        .unwrap()
        .temperature
        .value();
    let mw_sucrose = 342.2965;
    let sucrose_refs: &[(f64, f64, f64)] = &[
        // (wt_frac, bpe_ref_K, tol_K)
        (0.20, 0.35, 0.15),
        (0.40, 1.15, 0.20),
        (0.50, 1.85, 0.25),
        (0.60, 3.10, 0.30),
        (0.70, 5.50, 0.45),
    ];
    let mut prev_bpe = 0.0;
    for &(w, bpe_ref, tol) in sucrose_refs {
        let n_s = w / mw_sucrose;
        let n_w = (1.0 - w) / mw_h2o;
        let x_w = n_w / (n_w + n_s);
        let bp = sucrose_vle
            .bubble_point_t(p_1atm, &[x_w, 1.0 - x_w], opts)
            .unwrap();
        let bpe = bp.temperature.value() - t_pure_h2o;
        assert!(
            bpe > prev_bpe,
            "sucrose BPE must increase monotonically with concentration"
        );
        assert!(
            (bpe - bpe_ref).abs() <= tol,
            "sucrose {w:.0}% BPE mismatch: calc={bpe:.3} K, ref={bpe_ref:.3} K"
        );
        // Vapor must be 100% water (sucrose is non-volatile).
        assert!((bp.vapor_y[0] - 1.0).abs() < 1e-14 && bp.vapor_y[1] == 0.0);
        prev_bpe = bpe;
    }

    // 2. Water-NaCl brine BPE on dissociated-ion basis (M_ion = 29.22 g/mol).
    let nacl_vle = NrtlModelV1::water_nacl_vle_v1();
    let mw_nacl_ion = 29.22;
    let mut prev_brine_bpe = 0.0;
    for w in [0.05, 0.10, 0.15, 0.20, 0.26] {
        let n_s = w / mw_nacl_ion;
        let n_w = (1.0 - w) / mw_h2o;
        let x_w = n_w / (n_w + n_s);
        let bp = nacl_vle
            .bubble_point_t(p_1atm, &[x_w, 1.0 - x_w], opts)
            .unwrap();
        let bpe = bp.temperature.value() - t_pure_h2o;
        assert!(
            bpe > prev_brine_bpe,
            "brine BPE must increase monotonically"
        );
        assert!((bp.vapor_y[0] - 1.0).abs() < 1e-14 && bp.vapor_y[1] == 0.0);
        prev_brine_bpe = bpe;
    }
    assert!(
        (7.0..=10.5).contains(&prev_brine_bpe),
        "saturated 26 wt% NaCl brine BPE expected ~8.6 K, got {prev_brine_bpe:.2} K"
    );

    // 3. Two-phase isothermal flash of boiling 50 wt% sucrose syrup concentrates liquid and vents pure steam.
    let flash_syrup = sucrose_vle
        .flash_pt(Temperature::new(376.5), p_1atm, &[0.95, 0.05], opts)
        .unwrap();
    assert_eq!(flash_syrup.phase_regime, FlashPhaseRegimeV1::TwoPhase);
    assert!(flash_syrup.vapor_fraction > 0.0 && flash_syrup.vapor_fraction < 1.0);
    assert!((flash_syrup.vapor_y[0] - 1.0).abs() < 1e-12 && flash_syrup.vapor_y[1] == 0.0);
    assert!(
        flash_syrup.liquid_x[1] > 0.05,
        "boiling off steam must concentrate sucrose"
    );
}

#[test]
fn g0_g2_gab_and_oswin_food_sorption_roundtrip_and_isosteric_heat_oracles() {
    let gab_foods: &[(&str, GabIsothermV1)] = &[
        ("potato_starch", GabIsothermV1::potato_starch_v1()),
        ("bread_crumb", GabIsothermV1::bread_crumb_v1()),
        ("lean_beef", GabIsothermV1::lean_beef_v1()),
        ("potato_tissue", GabIsothermV1::potato_tissue_v1()),
    ];
    let oswin_foods: &[(&str, OswinIsothermV1)] = &[
        ("bread_oswin", OswinIsothermV1::bread_oswin_v1()),
        ("meat_oswin", OswinIsothermV1::meat_oswin_v1()),
        ("starch_oswin", OswinIsothermV1::starch_oswin_v1()),
    ];

    let t25 = Temperature::new(298.15);
    let t80 = Temperature::new(353.15);
    let h_fg_25 = SaturatedWaterState::at_temperature(t25)
        .unwrap()
        .latent_enthalpy()
        .value();

    // GAB food sorption checks at 25 C and 80 C.
    for &(name, model) in gab_foods {
        let mut prev_x25 = 0.0;
        let mut prev_q_low = f64::INFINITY;
        for &aw in &[1e-4, 0.05, 0.15, 0.30, 0.50, 0.70, 0.90, 0.995] {
            let ev25 = model.evaluate_from_activity(aw, t25).unwrap();
            let ev80 = model.evaluate_from_activity(aw, t80).unwrap();
            let x25 = ev25.moisture_dry_basis();
            let x80 = ev80.moisture_dry_basis();

            // Monotonic increase with a_w, and moisture decreases as T rises at fixed a_w.
            assert!(x25 > prev_x25, "{name} X_eq must increase with a_w");
            assert!(
                x80 < x25,
                "{name} X_eq at 80 C ({x80:.4}) must lie below 25 C ({x25:.4}) at a_w={aw}"
            );
            prev_x25 = x25;

            // Roundtrip |a_w(X_eq(a_w, T), T) - a_w| < 1e-12.
            let aw_rt25 = model.water_activity(x25, t25).unwrap();
            let aw_rt80 = model.water_activity(x80, t80).unwrap();
            assert!(
                (aw_rt25 - aw).abs() < 1e-12,
                "{name} GAB 25 C roundtrip error at aw={aw}: got {aw_rt25}"
            );
            assert!(
                (aw_rt80 - aw).abs() < 1e-12,
                "{name} GAB 80 C roundtrip error at aw={aw}: got {aw_rt80}"
            );

            // Net isosteric heat q_st^net > 0, monotonically decreasing with moisture,
            // and total Q_st > h_fg at low moisture.
            let q_net = ev25.net_isosteric_heat_j_per_mol();
            assert!(
                q_net > 0.0,
                "{name} net q_st must be positive at aw={aw}, got {q_net}"
            );
            assert!(
                q_net <= prev_q_low + 1e-9,
                "{name} net q_st must decrease as moisture rises: {q_net} > {prev_q_low}"
            );
            prev_q_low = q_net;

            if aw <= 0.2 {
                let q_tot = model.total_isosteric_heat_j_per_kg(x25, t25).unwrap();
                assert!(
                    q_tot > h_fg_25 + 1e5,
                    "{name} low-moisture Q_st ({q_tot:.0} J/kg) must exceed h_fg ({h_fg_25:.0} J/kg)"
                );
            }

            // Analytical Clausius-Clapeyron q_st matches central finite difference < 1e-8.
            let dt = 1e-3;
            let aw_p = model
                .water_activity(x25, Temperature::new(t25.value() + dt))
                .unwrap();
            let aw_m = model
                .water_activity(x25, Temperature::new(t25.value() - dt))
                .unwrap();
            let q_fd = R * t25.value() * t25.value() * (ln(aw_p) - ln(aw_m)) / (2.0 * dt);
            assert!(
                (q_net - q_fd).abs() / q_net.abs().max(1.0) < 1e-8,
                "{name} GAB q_st FD check at aw={aw}: analytic={q_net}, fd={q_fd}"
            );
        }

        // Zero-moisture continuous limit.
        let dry = model.evaluate_from_moisture(0.0, t25).unwrap();
        assert_eq!(dry.water_activity(), 0.0);
        assert!(dry.daw_dx() > 0.0);
        assert!(dry.net_isosteric_heat_j_per_mol() > 0.0);
    }

    // Oswin food sorption checks (both linear and Arrhenius variants + unified dispatch enum).
    for &(name, model) in oswin_foods {
        let unified = SorptionIsothermV1::Oswin(model);
        assert_eq!(unified.water_activity(0.0, t25).unwrap(), 0.0);
        for &aw in &[1e-4, 0.05, 0.25, 0.55, 0.85, 0.995] {
            let ev25 = unified.evaluate_from_activity(aw, t25).unwrap();
            let x25 = ev25.moisture_dry_basis();
            let aw_rt = unified.water_activity(x25, t25).unwrap();
            assert!(
                (aw_rt - aw).abs() < 1e-12,
                "{name} Oswin roundtrip error at aw={aw}: got {aw_rt}"
            );
            assert!(
                ev25.net_isosteric_heat_j_per_mol() > 0.0,
                "{name} Oswin net q_st must be positive at aw={aw}"
            );
            let q_tot = unified.total_isosteric_heat_j_per_kg(x25, t25).unwrap();
            let expected_tot =
                h_fg_25 + ev25.net_isosteric_heat_j_per_mol() / WATER_MOLAR_MASS_KG_PER_MOL_V1;
            assert!((q_tot - expected_tot).abs() < 1e-6);

            let dt = 1e-3;
            let aw_p = model
                .water_activity(x25, Temperature::new(t25.value() + dt))
                .unwrap();
            let aw_m = model
                .water_activity(x25, Temperature::new(t25.value() - dt))
                .unwrap();
            let q_fd = R * t25.value() * t25.value() * (ln(aw_p) - ln(aw_m)) / (2.0 * dt);
            assert!(
                (ev25.net_isosteric_heat_j_per_mol() - q_fd).abs()
                    / ev25.net_isosteric_heat_j_per_mol().abs().max(1.0)
                    < 1e-8,
                "{name} Oswin q_st FD check at aw={aw}: analytic={}, fd={q_fd}",
                ev25.net_isosteric_heat_j_per_mol()
            );
        }
    }
}

#[test]
fn matdb_card_roundtrip_and_bitwise_determinism() {
    let prov = Provenance {
        source: "KM0.1 test fixture".into(),
        license: "test".into(),
        artifact: None,
    };

    // 1. NRTL from_card matches direct ethanol_water_aspen_v1 bit-for-bit.
    let mut nrtl_params = BTreeMap::new();
    nrtl_params.insert(
        "a_0_1".into(),
        LawParameter {
            value: -0.8009,
            dims: Dims::NONE,
        },
    );
    nrtl_params.insert(
        "a_1_0".into(),
        LawParameter {
            value: 3.4578,
            dims: Dims::NONE,
        },
    );
    nrtl_params.insert(
        "b_0_1".into(),
        LawParameter {
            value: 246.18,
            dims: Temperature::DIMS,
        },
    );
    nrtl_params.insert(
        "b_1_0".into(),
        LawParameter {
            value: -586.0809,
            dims: Temperature::DIMS,
        },
    );
    nrtl_params.insert(
        "alpha_0_1".into(),
        LawParameter {
            value: 0.3,
            dims: Dims::NONE,
        },
    );
    let nrtl_card = ConstitutiveModelCard {
        law: LawId(NRTL_LAW_ID_V1.into()),
        law_version: 1,
        parameters: nrtl_params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 273.15, 423.15),
        sources: Vec::new(),
        provenance: prov.clone(),
    };
    let species = vec![
        SpeciesId::new("C2H5OH").unwrap(),
        SpeciesId::new("H2O").unwrap(),
    ];
    let nrtl_loaded = NrtlModelV1::from_card(&nrtl_card, species.clone()).unwrap();
    let nrtl_direct = NrtlModelV1::ethanol_water_aspen_v1();
    let ev_a = nrtl_loaded
        .evaluate(Temperature::new(351.30), &[0.42, 0.58])
        .unwrap();
    let ev_b = nrtl_direct
        .evaluate(Temperature::new(351.30), &[0.42, 0.58])
        .unwrap();
    assert_eq!(
        ev_a.reduced_excess_gibbs().to_bits(),
        ev_b.reduced_excess_gibbs().to_bits()
    );
    assert_eq!(
        ev_a.ln_activity_coefficients()[0].to_bits(),
        ev_b.ln_activity_coefficients()[0].to_bits()
    );

    // 2. UNIQUAC from_card matches direct ethanol_water_aspen_v1 bit-for-bit.
    let mut uq_params = BTreeMap::new();
    uq_params.insert(
        "a_0_1".into(),
        LawParameter {
            value: -2.4936,
            dims: Dims::NONE,
        },
    );
    uq_params.insert(
        "a_1_0".into(),
        LawParameter {
            value: 2.0046,
            dims: Dims::NONE,
        },
    );
    uq_params.insert(
        "b_0_1".into(),
        LawParameter {
            value: 756.9477,
            dims: Temperature::DIMS,
        },
    );
    uq_params.insert(
        "b_1_0".into(),
        LawParameter {
            value: -728.9705,
            dims: Temperature::DIMS,
        },
    );
    let uq_card = ConstitutiveModelCard {
        law: LawId(UNIQUAC_LAW_ID_V1.into()),
        law_version: 1,
        parameters: uq_params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 273.15, 423.15),
        sources: Vec::new(),
        provenance: prov.clone(),
    };
    let uq_comps = vec![
        UniquacComponentV1::standard(species[0].clone(), 2.1055, 1.972).unwrap(),
        UniquacComponentV1::standard(species[1].clone(), 0.92, 1.40).unwrap(),
    ];
    let uq_loaded = UniquacModelV1::from_card(&uq_card, uq_comps).unwrap();
    let uq_direct = UniquacModelV1::ethanol_water_aspen_v1();
    let uev_a = uq_loaded
        .evaluate(Temperature::new(351.30), &[0.42, 0.58])
        .unwrap();
    let uev_b = uq_direct
        .evaluate(Temperature::new(351.30), &[0.42, 0.58])
        .unwrap();
    assert_eq!(
        uev_a.reduced_excess_gibbs().to_bits(),
        uev_b.reduced_excess_gibbs().to_bits()
    );

    // 3. GAB and Oswin from_card bit-exact matches.
    let mut gab_params = BTreeMap::new();
    gab_params.insert(
        "x_m".into(),
        LawParameter {
            value: 0.095,
            dims: Dims::NONE,
        },
    );
    gab_params.insert(
        "c0".into(),
        LawParameter {
            value: 0.015,
            dims: Dims::NONE,
        },
    );
    gab_params.insert(
        "delta_h_c".into(),
        LawParameter {
            value: 17_500.0,
            dims: Dims([2, 1, -2, 0, 0, -1]),
        },
    );
    gab_params.insert(
        "k0".into(),
        LawParameter {
            value: 0.60,
            dims: Dims::NONE,
        },
    );
    gab_params.insert(
        "delta_h_k".into(),
        LawParameter {
            value: 650.0,
            dims: Dims([2, 1, -2, 0, 0, -1]),
        },
    );
    let gab_card = ConstitutiveModelCard {
        law: LawId(GAB_ISOTHERM_LAW_ID_V1.into()),
        law_version: 1,
        parameters: gab_params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 273.15, 423.15),
        sources: Vec::new(),
        provenance: prov.clone(),
    };
    let gab_loaded = GabIsothermV1::from_card(&gab_card).unwrap();
    let gab_direct = GabIsothermV1::potato_starch_v1();
    let g_a = gab_loaded
        .evaluate_from_activity(0.55, Temperature::new(310.0))
        .unwrap();
    let g_b = gab_direct
        .evaluate_from_activity(0.55, Temperature::new(310.0))
        .unwrap();
    assert_eq!(
        g_a.moisture_dry_basis().to_bits(),
        g_b.moisture_dry_basis().to_bits()
    );
    assert_eq!(
        g_a.net_isosteric_heat_j_per_mol().to_bits(),
        g_b.net_isosteric_heat_j_per_mol().to_bits()
    );

    let mut oswin_params = BTreeMap::new();
    oswin_params.insert(
        "k_ref".into(),
        LawParameter {
            value: 0.105,
            dims: Dims::NONE,
        },
    );
    oswin_params.insert(
        "dk_dt".into(),
        LawParameter {
            value: -4.5e-4,
            dims: Dims([0, 0, 0, -1, 0, 0]),
        },
    );
    oswin_params.insert(
        "t_ref".into(),
        LawParameter {
            value: 298.15,
            dims: Temperature::DIMS,
        },
    );
    oswin_params.insert(
        "n".into(),
        LawParameter {
            value: 0.43,
            dims: Dims::NONE,
        },
    );
    let oswin_card = ConstitutiveModelCard {
        law: LawId(OSWIN_ISOTHERM_LAW_ID_V1.into()),
        law_version: 1,
        parameters: oswin_params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 273.15, 423.15),
        sources: Vec::new(),
        provenance: prov,
    };
    let oswin_loaded = OswinIsothermV1::from_card(&oswin_card).unwrap();
    let oswin_direct = OswinIsothermV1::bread_oswin_v1();
    let o_a = oswin_loaded
        .evaluate_from_activity(0.65, Temperature::new(310.0))
        .unwrap();
    let o_b = oswin_direct
        .evaluate_from_activity(0.65, Temperature::new(310.0))
        .unwrap();
    assert_eq!(
        o_a.moisture_dry_basis().to_bits(),
        o_b.moisture_dry_basis().to_bits()
    );
}

#[test]
fn planted_negatives_return_typed_errors_with_zero_panics_or_nans() {
    let nrtl = NrtlModelV1::ethanol_water_aspen_v1();
    let uniquac = UniquacModelV1::ethanol_water_aspen_v1();
    let vle = NrtlModelV1::ethanol_water_vle_v1();
    let t_ok = Temperature::new(350.0);
    let opts = VleSolverOptionsV1::default();

    // Invalid mole fractions: negative, > 1, sum != 1, wrong length, NaN, inf.
    let bad_compositions: &[&[f64]] = &[
        &[-0.1, 1.1],
        &[1.2, -0.2],
        &[0.3, 0.3],
        &[0.6, 0.6],
        &[0.5],
        &[0.5, 0.3, 0.2],
        &[f64::NAN, 1.0],
        &[f64::INFINITY, 0.0],
        &[f64::NEG_INFINITY, 1.0],
    ];
    for &bad_x in bad_compositions {
        assert!(
            matches!(
                nrtl.evaluate(t_ok, bad_x),
                Err(ExcessGibbsErrorV1::InvalidComposition(_))
            ),
            "NRTL must refuse composition {bad_x:?}"
        );
        assert!(
            matches!(
                uniquac.evaluate(t_ok, bad_x),
                Err(ExcessGibbsErrorV1::InvalidComposition(_))
            ),
            "UNIQUAC must refuse composition {bad_x:?}"
        );
        assert!(
            matches!(
                vle.bubble_point_t(Pressure::new(101_325.0), bad_x, opts),
                Err(ExcessGibbsErrorV1::InvalidComposition(_))
            ),
            "VLE bubble_point_t must refuse composition {bad_x:?}"
        );
    }

    // Invalid temperatures: T <= 0, NaN, inf, out of domain.
    for bad_t in [
        0.0,
        -25.0,
        200.0,
        500.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        let t = Temperature::new(bad_t);
        assert!(matches!(
            nrtl.evaluate(t, &[0.5, 0.5]),
            Err(ExcessGibbsErrorV1::TemperatureOutOfDomain { .. })
        ));
        assert!(matches!(
            uniquac.evaluate(t, &[0.5, 0.5]),
            Err(ExcessGibbsErrorV1::TemperatureOutOfDomain { .. })
        ));
    }

    // Invalid pressures: P <= 0, NaN, inf.
    for bad_p in [0.0, -101_325.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let p = Pressure::new(bad_p);
        assert!(matches!(
            vle.k_values(t_ok, p, &[0.5, 0.5]),
            Err(ExcessGibbsErrorV1::InvalidPressure { .. })
        ));
        assert!(matches!(
            vle.bubble_point_t(p, &[0.5, 0.5], opts),
            Err(ExcessGibbsErrorV1::InvalidPressure { .. })
        ));
        assert!(matches!(
            vle.flash_pt(t_ok, p, &[0.5, 0.5], opts),
            Err(ExcessGibbsErrorV1::InvalidPressure { .. })
        ));
    }

    // Invalid non-volatile liquid cannot boil, and non-volatile vapor cannot condense in dew solve.
    let sucrose_vle = NrtlModelV1::water_sucrose_vle_v1();
    assert!(matches!(
        sucrose_vle.bubble_point_t(Pressure::new(101_325.0), &[0.0, 1.0], opts),
        Err(ExcessGibbsErrorV1::InvalidComposition(_))
    ));
    assert!(matches!(
        sucrose_vle.dew_point_t(Pressure::new(101_325.0), &[0.5, 0.5], opts),
        Err(ExcessGibbsErrorV1::InvalidComposition(_))
    ));

    // Invalid NRTL/UNIQUAC model construction: duplicate species, asymmetric alpha, bad r/q.
    let sp_dup = vec![
        SpeciesId::new("H2O").unwrap(),
        SpeciesId::new("H2O").unwrap(),
    ];
    assert!(matches!(
        NrtlModelV1::new(
            sp_dup,
            vec![0.0; 4],
            vec![0.0; 4],
            vec![0.0; 4],
            Temperature::new(280.0),
            Temperature::new(380.0),
        ),
        Err(ExcessGibbsErrorV1::InvalidModel(_))
    ));
    assert!(matches!(
        NrtlModelV1::new(
            vec![SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            vec![0.0; 4],
            vec![0.0; 4],
            vec![0.0, 0.3, 0.4, 0.0],
            Temperature::new(280.0),
            Temperature::new(380.0),
        ),
        Err(ExcessGibbsErrorV1::InvalidModel(_))
    ));
    assert!(matches!(
        UniquacComponentV1::standard(SpeciesId::new("A").unwrap(), 0.0, 1.0),
        Err(ExcessGibbsErrorV1::InvalidModel(_))
    ));
    assert!(matches!(
        ModifiedRaoultVleV1::new(nrtl.clone(), vec![VaporPressureLawV1::iapws_water()]),
        Err(ExcessGibbsErrorV1::InvalidModel(_))
    ));

    // Sorption planted negatives: a_w <= 0, a_w >= 1, T <= 0 / out of domain, X_w < 0,
    // moisture above GAB saturation ceiling, and singular GAB denominator 1 - K*a_w <= 0.
    let gab = GabIsothermV1::potato_starch_v1();
    let oswin = OswinIsothermV1::bread_oswin_v1();
    let t_sorp = Temperature::new(298.15);

    for bad_aw in [
        0.0,
        -0.1,
        1.0,
        1.05,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert!(matches!(
            gab.equilibrium_moisture(bad_aw, t_sorp),
            Err(SorptionErrorV1::WaterActivityOutOfDomain { .. })
        ));
        assert!(matches!(
            oswin.equilibrium_moisture(bad_aw, t_sorp),
            Err(SorptionErrorV1::WaterActivityOutOfDomain { .. })
        ));
    }

    for bad_t in [0.0, -10.0, 200.0, 500.0, f64::NAN, f64::INFINITY] {
        let t = Temperature::new(bad_t);
        assert!(matches!(
            gab.equilibrium_moisture(0.5, t),
            Err(SorptionErrorV1::TemperatureOutOfDomain { .. })
        ));
        assert!(matches!(
            oswin.equilibrium_moisture(0.5, t),
            Err(SorptionErrorV1::TemperatureOutOfDomain { .. })
        ));
    }

    for bad_x in [-1e-6, -0.5, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            gab.water_activity(bad_x, t_sorp),
            Err(SorptionErrorV1::MoistureOutOfDomain { .. })
        ));
        assert!(matches!(
            oswin.water_activity(bad_x, t_sorp),
            Err(SorptionErrorV1::MoistureOutOfDomain { .. })
        ));
    }

    // GAB with K(T) < 1 has a finite saturation moisture ceiling X_w(1^-, T): exceeding it refuses cleanly.
    assert!(matches!(
        gab.water_activity(10.0, t_sorp),
        Err(SorptionErrorV1::MoistureOutOfDomain { .. })
    ));

    // Singular GAB denominator 1 - K(T)*a_w <= 0 (K = 1.25, a_w = 0.80 => 1 - K*a_w = 0).
    let singular_gab = GabIsothermV1::new(
        0.09,
        12.0,
        0.0,
        1.25,
        0.0,
        Temperature::new(273.15),
        Temperature::new(400.0),
    )
    .unwrap();
    assert!(matches!(
        singular_gab.equilibrium_moisture(0.80, t_sorp),
        Err(SorptionErrorV1::SingularGabDenominator { .. })
    ));
    assert!(matches!(
        singular_gab.equilibrium_moisture(0.95, t_sorp),
        Err(SorptionErrorV1::SingularGabDenominator { .. })
    ));

    // Dimensional mismatch on matdb cards must be refused rather than silently misinterpreting units.
    let mut bad_nrtl_params = BTreeMap::new();
    bad_nrtl_params.insert(
        "a_0_1".into(),
        LawParameter {
            value: -0.8,
            dims: Dims::NONE,
        },
    );
    bad_nrtl_params.insert(
        "a_1_0".into(),
        LawParameter {
            value: 3.4,
            dims: Dims::NONE,
        },
    );
    bad_nrtl_params.insert(
        "b_0_1".into(),
        LawParameter {
            value: 246.0,
            dims: Dims::NONE,
        },
    ); // wrong dims!
    bad_nrtl_params.insert(
        "b_1_0".into(),
        LawParameter {
            value: -586.0,
            dims: Temperature::DIMS,
        },
    );
    bad_nrtl_params.insert(
        "alpha_0_1".into(),
        LawParameter {
            value: 0.3,
            dims: Dims::NONE,
        },
    );
    let bad_nrtl_card = ConstitutiveModelCard {
        law: LawId(NRTL_LAW_ID_V1.into()),
        law_version: 1,
        parameters: bad_nrtl_params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 273.15, 423.15),
        sources: Vec::new(),
        provenance: Provenance {
            source: "planted negative".into(),
            license: "test".into(),
            artifact: None,
        },
    };
    assert!(matches!(
        NrtlModelV1::from_card(
            &bad_nrtl_card,
            vec![
                SpeciesId::new("C2H5OH").unwrap(),
                SpeciesId::new("H2O").unwrap()
            ]
        ),
        Err(ExcessGibbsErrorV1::InvalidModel(_))
    ));
}
