//! Integration, physical-oracle, conservation, and conformance tests for
//! [`fs_multiphase::MultiphaseVesselComponent`] (`KM2.2`).

use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_multiphase::{
    BoilingCurveV1, BoilingFluidKindV1, BoilingSurfaceCardV1, COOKING_OIL_LAW_ID_V1,
    LIQUID_STREAM_WIDTH, MultiphaseVesselComponent, MultiphaseVesselConfig, MultiphaseVesselState,
    OilCardV1, OpenVentSpecV1, PROBE_FFA_PCT, PROBE_FOAM_H_M, PROBE_FOND_ADHERED_KG,
    PROBE_M_LIQ_KG, PROBE_SAFETY_OR_BOILOVER_FLAG, PROBE_T_COLDZONE_K, PROBE_T_LIQ_K,
    PROBE_T_SMOKE_K, PROBE_TPM_PCT, PROBE_VISC_PA_S, PROBE_X_ETOH, PROBE_X_SOLIDS,
    VAPOR_STREAM_WIDTH, VesselErrorV1, VesselInitialStateV1,
};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{
    Component, ContractViolation, Proposal, STREAM_CARRIED, STREAM_ENTHALPY_RATE, STREAM_MASS_RATE,
    audit_carried_reactions, audit_carried_step_with_production, audit_mass_step, audit_step,
};
use fs_qty::{Dims, Pressure, Temperature};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 22,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn audit_full_step(
    comp: &MultiphaseVesselComponent,
    old: &MultiphaseVesselState,
    new: &MultiphaseVesselState,
    dt_s: f64,
    inputs: &[f64],
    outputs: &[f64],
    receipt: &fs_ports::component::Receipt,
) {
    let ports = comp.ports();
    let layout = comp.layout();

    let e_old = comp.stored_energy_j(old).expect("e_old");
    let e_new = comp.stored_energy_j(new).expect("e_new");
    let step_audit = audit_step(ports, layout, dt_s, inputs, outputs, receipt, e_old, e_new)
        .expect("First-Law audit_step must pass");
    let e_scale = e_old.abs().max(e_new.abs()).max(1.0);
    assert!(
        step_audit.balance_defect_j.abs() <= 1.0e-11 * e_scale,
        "energy defect {} exceeded 1e-11 relative to {}",
        step_audit.balance_defect_j,
        e_scale
    );

    let m_old = comp.stored_mass_kg(old).expect("m_old");
    let m_new = comp.stored_mass_kg(new).expect("m_new");
    let mass_audit = audit_mass_step(ports, layout, dt_s, inputs, outputs, m_old, m_new)
        .expect("mass conservation audit_mass_step must pass")
        .expect("vessel exposes stored mass");
    let m_scale = m_old.unwrap_or(0.0).max(m_new.unwrap_or(0.0)).max(1.0);
    assert!(
        mass_audit.balance_defect_kg.abs() <= 1.0e-11 * m_scale,
        "mass defect {} exceeded 1e-11",
        mass_audit.balance_defect_kg
    );

    let mut productions = Vec::with_capacity(layout.carried_len());
    for carried in 0..layout.carried_len() {
        let c_old = comp.stored_carried_quantity(old, carried).expect("c_old");
        let c_new = comp.stored_carried_quantity(new, carried).expect("c_new");
        let prod = comp
            .carried_production(old, new, dt_s, carried)
            .expect("prod");
        audit_carried_step_with_production(
            ports, layout, dt_s, inputs, outputs, carried, c_old, c_new, prod,
        )
        .unwrap_or_else(|e| panic!("carried slot {carried} failed conservation audit: {e}"));
        productions.push(prod);
    }

    audit_carried_reactions(layout, &productions, comp.reaction_basis())
        .expect("carried reactions must close elemental & mass stoichiometry");
}

#[test]
fn fs_ports_conformance_and_durable_checkpoint_round_trip() {
    let comp =
        MultiphaseVesselComponent::new(MultiphaseVesselConfig::saucepan_default()).expect("vessel");
    let state = comp
        .initial_state(VesselInitialStateV1::wine_deglaze(0.60, 0.11, 0.025, 348.0))
        .expect("state");
    let stale = comp
        .initial_state(VesselInitialStateV1::water(0.25, 310.0))
        .expect("stale");

    // Nonzero wall heat + wine pour-in on liquid_stream + active stirring.
    let mut pour = MultiphaseVesselComponent::liquid_stream_for_pure_species(0, 0.015, 300.0, None);
    let etoh = MultiphaseVesselComponent::liquid_stream_for_pure_species(1, 0.002, 300.0, None);
    pour[STREAM_MASS_RATE] += etoh[STREAM_MASS_RATE];
    pour[STREAM_ENTHALPY_RATE] += etoh[STREAM_ENTHALPY_RATE];
    pour[STREAM_CARRIED + 1] += etoh[STREAM_CARRIED + 1];

    let inputs = comp
        .build_inputs(
            850.0,
            &[-120.0],
            pour,
            [0.0; VAPOR_STREAM_WIDTH],
            90.0,
            101_325.0,
        )
        .expect("inputs");

    with_cx(|cx| {
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: 2.0,
            inputs: &inputs,
        };
        let observed = check_component(&comp, &case, cx).expect("check_component must pass");
        assert!(matches!(observed, Observed::Stepped { .. }));

        // Durable StateIdentityCharterV2 codec roundtrip & deterministic identity_words.
        let charter = MultiphaseVesselComponent::STATE_CHARTER.expect("charter declared");
        assert_eq!(charter.owner, "fs-multiphase::MultiphaseVessel");
        assert_eq!(charter.state_family, "stirred-vle-v1");
        assert_eq!(charter.codec_version, 1);

        let mut next = stale.clone();
        let mut out = vec![0.0; comp.layout().output_len()];
        let prop = comp
            .propose(&state, 2.0, &inputs, &mut next, &mut out, cx)
            .expect("propose");
        assert!(matches!(prop, Proposal::Step(_)));

        let mut words = Vec::new();
        assert!(comp.encode_state(&next, &mut words));
        let decoded = comp.decode_state(&words).expect("decode_state");
        assert_eq!(decoded, next);
    });
}

#[test]
fn emergent_usda_wine_simmer_ethanol_retention_and_fond_dissolution_conservation() {
    // Drive a saucepan containing 1.0 kg red wine (~13% ABV => 10.5 wt% EtOH) + 30 g adhered fond
    // against a simmering SS304 pan wall held at T_wall = 101.5 C using KM2.1 BoilingCurveV1
    // (Rohsenow + KM0.1 NRTL Schluender-Thome) and open-pot donor-rule venting.
    // ZERO hardcoded decay tables: USDA retention (~40-50% at 15m, ~25-40% at 30m, ~15-30% at 60m)
    // emerges purely from KM0.1 NRTL VLE + T_bub rise + KM2.1 boiling flux.
    let mut cfg = MultiphaseVesselConfig::saucepan_default();
    cfg.base_area_m2 = 0.0314;
    cfg.free_surface_area_m2 = 0.0314;
    let comp = MultiphaseVesselComponent::new(cfg).expect("saucepan");

    let mut state = comp
        .initial_state(VesselInitialStateV1::wine_deglaze(1.0, 0.105, 0.030, 364.2))
        .expect("wine state");
    let initial_etoh_kg = state.m_etoh_liq_kg;
    let initial_total_fond_kg = state.m_fond_adhered_kg + state.m_fond_dissolved_kg;

    let surface = BoilingSurfaceCardV1::stainless_steel_304_v1();
    let p_atm = 101_325.0;
    let t_wall_k = 273.15 + 101.5;
    let dt_s = 10.0;
    let mut out = vec![0.0; comp.layout().output_len()];

    let mut ret_15m = None;
    let mut ret_30m = None;
    let mut ret_60m = None;

    with_cx(|cx| {
        for step in 1..=360 {
            let x_etoh = state.ethanol_mole_fraction();
            let curve = BoilingCurveV1::prepare(
                BoilingFluidKindV1::EthanolWaterNrtl {
                    ethanol_mole_fraction: x_etoh,
                },
                surface,
                Pressure::new(p_atm),
            )
            .expect("KM2.1 boiling curve");
            let eval = curve
                .evaluate(Temperature::new(t_wall_k))
                .expect("evaluate boiling");
            // Gentle simmer active nucleate floor fraction (~18% of pan base).
            let q_wall = 0.18 * comp.config().base_area_m2 * eval.heat_flux_w_m2;

            let inputs = comp
                .open_vent_inputs(
                    &state,
                    dt_s,
                    q_wall,
                    &[0.0],
                    [0.0; LIQUID_STREAM_WIDTH],
                    [0.0; VAPOR_STREAM_WIDTH],
                    30.0,
                    p_atm,
                    OpenVentSpecV1::open_pot(),
                )
                .expect("vent inputs");

            let mut next = state.clone();
            let prop = comp
                .propose(&state, dt_s, &inputs, &mut next, &mut out, cx)
                .expect("step");
            let Proposal::Step(receipt) = prop else {
                panic!("unexpected shorten at step {step}");
            };
            audit_full_step(&comp, &state, &next, dt_s, &inputs, &out, &receipt);

            // Exact internal fond ownership closure: adhered + dissolved is constant to < 1e-12.
            let fond_sum = next.m_fond_adhered_kg + next.m_fond_dissolved_kg;
            assert!(
                (fond_sum - initial_total_fond_kg).abs() < 1.0e-12,
                "fond conservation violated: {fond_sum} vs {initial_total_fond_kg}"
            );

            state = next;
            let ret_pct = 100.0 * (state.m_etoh_liq_kg / initial_etoh_kg);
            if step == 90 {
                ret_15m = Some(ret_pct);
            } else if step == 180 {
                ret_30m = Some(ret_pct);
            } else if step == 360 {
                ret_60m = Some(ret_pct);
            }
        }
    });

    let r15 = ret_15m.expect("15m");
    let r30 = ret_30m.expect("30m");
    let r60 = ret_60m.expect("60m");

    assert!(
        (25.0..=65.0).contains(&r15),
        "15 min ethanol retention out of USDA simmer VLE band: {r15:.1}%"
    );
    assert!(
        (15.0..=50.0).contains(&r30) && r30 < r15,
        "30 min ethanol retention out of band: {r30:.1}% (15m={r15:.1}%)"
    );
    assert!(
        (5.0..=38.0).contains(&r60) && r60 < r30,
        "60 min ethanol retention out of band: {r60:.1}% (30m={r30:.1}%)"
    );

    assert!(
        state.m_fond_dissolved_kg > 0.85 * initial_total_fond_kg,
        "expected most fond dissolved, dissolved={} of {}",
        state.m_fond_dissolved_kg,
        initial_total_fond_kg
    );
    let probes = &out[comp.layout().output_range(comp.probes_port())];
    assert!((probes[PROBE_FOND_ADHERED_KG] - state.m_fond_adhered_kg).abs() < 1.0e-15);
    assert!(probes[PROBE_X_ETOH] < 0.025);
}

#[test]
fn bikerman_pasta_pot_boilover_triggers_and_collapses_on_power_drop() {
    let comp =
        MultiphaseVesselComponent::new(MultiphaseVesselConfig::stockpot_default()).expect("pot");
    let mut state = comp
        .initial_state(VesselInitialStateV1::pasta_pot(7.5, 0.080, 0.075, 372.5))
        .expect("pasta pot");
    let mut out = vec![0.0; comp.layout().output_len()];
    let p_atm = 101_325.0;
    let dt_s = 2.0;

    let mut boil_over_seen = false;
    let mut peak_foam_m = 0.0_f64;

    with_cx(|cx| {
        for _ in 0..35 {
            let inputs = comp
                .open_vent_inputs(
                    &state,
                    dt_s,
                    3_500.0,
                    &[0.0],
                    [0.0; LIQUID_STREAM_WIDTH],
                    [0.0; VAPOR_STREAM_WIDTH],
                    0.0,
                    p_atm,
                    OpenVentSpecV1::open_pot(),
                )
                .expect("inputs");
            let mut next = state.clone();
            let Proposal::Step(receipt) = comp
                .propose(&state, dt_s, &inputs, &mut next, &mut out, cx)
                .expect("propose")
            else {
                panic!("unexpected shorten during pasta boil");
            };
            audit_full_step(&comp, &state, &next, dt_s, &inputs, &out, &receipt);
            let probes = &out[comp.layout().output_range(comp.probes_port())];
            peak_foam_m = peak_foam_m.max(probes[PROBE_FOAM_H_M]);
            if (probes[PROBE_SAFETY_OR_BOILOVER_FLAG] - 1.0).abs() < 1.0e-12 {
                boil_over_seen = true;
            }
            state = next;
        }
        assert!(
            boil_over_seen,
            "expected Bikerman boil_over_flag == 1.0 under 3.5 kW starchy boil (peak foam={peak_foam_m:.3} m)"
        );

        // Drop burner power to gentle 150 W simmer -> foam collapses and boil_over clears.
        for _ in 0..25 {
            let inputs = comp
                .open_vent_inputs(
                    &state,
                    dt_s,
                    150.0,
                    &[0.0],
                    [0.0; LIQUID_STREAM_WIDTH],
                    [0.0; VAPOR_STREAM_WIDTH],
                    0.0,
                    p_atm,
                    OpenVentSpecV1::open_pot(),
                )
                .expect("inputs");
            let mut next = state.clone();
            let Proposal::Step(receipt) = comp
                .propose(&state, dt_s, &inputs, &mut next, &mut out, cx)
                .expect("propose")
            else {
                panic!("unexpected shorten");
            };
            audit_full_step(&comp, &state, &next, dt_s, &inputs, &out, &receipt);
            state = next;
        }
    });

    let probes = &out[comp.layout().output_range(comp.probes_port())];
    assert!(
        probes[PROBE_FOAM_H_M] < 0.25 * peak_foam_m,
        "foam must collapse on power drop: now {} vs peak {}",
        probes[PROBE_FOAM_H_M],
        peak_foam_m
    );
    assert!(
        probes[PROBE_SAFETY_OR_BOILOVER_FLAG].abs() < 1.0e-12,
        "boil_over flag must clear after power drop, got {}",
        probes[PROBE_SAFETY_OR_BOILOVER_FLAG]
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn butter_emulsion_flashes_water_near_100c_then_browns_milk_solids_unlike_dry_ghee_and_canola() {
    let run_pan = |oil: OilCardV1| -> (Vec<f64>, MultiphaseVesselState) {
        let comp = MultiphaseVesselComponent::new(MultiphaseVesselConfig::butter_skillet_default(
            oil.clone(),
        ))
        .expect("pan");
        let mut state = comp
            .initial_state(VesselInitialStateV1::butter_pan(&oil, 0.100, 320.0))
            .expect("state");
        let mut out = vec![0.0; comp.layout().output_len()];
        let mut temps_c = Vec::new();
        with_cx(|cx| {
            let mut t: f64 = 0.0;
            let t_end: f64 = 150.0;
            let mut dt: f64 = 1.0;
            while t < t_end - 1.0e-9 {
                let step_dt = dt.min(t_end - t);
                let inputs = comp
                    .open_vent_inputs(
                        &state,
                        step_dt,
                        420.0,
                        &[0.0],
                        [0.0; LIQUID_STREAM_WIDTH],
                        [0.0; VAPOR_STREAM_WIDTH],
                        20.0,
                        101_325.0,
                        OpenVentSpecV1::open_pot(),
                    )
                    .expect("inputs");
                let mut next = state.clone();
                match comp
                    .propose(&state, step_dt, &inputs, &mut next, &mut out, cx)
                    .expect("propose")
                {
                    Proposal::Step(receipt) => {
                        audit_full_step(&comp, &state, &next, step_dt, &inputs, &out, &receipt);
                        t += step_dt;
                        state = next;
                        temps_c.push(state.t_liq_k - 273.15);
                        dt = (step_dt * 1.5).min(1.0);
                    }
                    Proposal::Shorten { dt_s } => {
                        assert!(dt_s > 0.0 && dt_s < step_dt);
                        dt = dt_s;
                    }
                }
            }
        });
        (temps_c, state)
    };

    let (t_unsalted, st_unsalted) = run_pan(OilCardV1::butter_unsalted());
    let (t_salted, st_salted) = run_pan(OilCardV1::butter_salted());
    let (t_ghee, st_ghee) = run_pan(OilCardV1::clarified_butter_ghee());
    let (_t_canola, st_canola) = run_pan(OilCardV1::canola_refined());

    let unsalted_plateau_steps = t_unsalted
        .iter()
        .filter(|&&tc| (99.0..=102.0).contains(&tc))
        .count();
    let salted_plateau_steps = t_salted
        .iter()
        .filter(|&&tc| (99.5..=104.5).contains(&tc))
        .count();
    let ghee_plateau_steps = t_ghee
        .iter()
        .filter(|&&tc| (99.0..=102.0).contains(&tc))
        .count();

    assert!(
        unsalted_plateau_steps >= 60,
        "unsalted butter must plateau near 100 C while flashing 17.5% water, got {unsalted_plateau_steps} steps"
    );
    assert!(
        salted_plateau_steps >= 50,
        "salted butter must plateau near 101-103 C while flashing 16% salted water, got {salted_plateau_steps} steps"
    );
    assert!(
        ghee_plateau_steps <= 5,
        "clarified butter/ghee has negligible water and must not plateau near 100 C, got {ghee_plateau_steps}"
    );

    assert!(
        st_unsalted.t_liq_k > 393.15,
        "unsalted butter must climb past 120 C after flash, got {} K",
        st_unsalted.t_liq_k
    );
    assert!(
        st_unsalted.milk_solid_browning > 0.05,
        "milk solids must brown in beurre noisette regime, got {}",
        st_unsalted.milk_solid_browning
    );
    assert!(
        st_salted.t_liq_k > 393.15 && st_salted.milk_solid_browning > 0.05,
        "salted butter must also enter milk-solid browning after water flash"
    );
    assert!(
        st_canola.t_liq_k > st_unsalted.t_liq_k + 40.0
            && st_ghee.t_liq_k > st_unsalted.t_liq_k + 40.0,
        "dry canola/ghee heat much faster than 17% water butter emulsion"
    );
    assert!(st_canola.milk_solid_browning.abs() < 1.0e-15);
}

#[test]
#[allow(clippy::too_many_lines)]
fn deep_fryer_multi_day_aging_cold_zone_stoichiometry_and_safety_flags() {
    let canola = OilCardV1::canola_refined();
    let comp =
        MultiphaseVesselComponent::new(MultiphaseVesselConfig::deep_fryer_default(canola.clone()))
            .expect("fryer");
    let mut state = comp
        .initial_state(VesselInitialStateV1::oil_fryer(
            &canola, 11.5, 4.5, 450.15, 328.15,
        ))
        .expect("fryer state");
    let fresh_visc = comp.effective_viscosity_pa_s(&state, 0.0);
    let fresh_smoke = comp.smoke_point_k(&state);
    let mut out = vec![0.0; comp.layout().output_len()];

    let dt_s = 60.0;
    let total_steps: usize = 40 * 60;
    let mut t_cz_after_4h = None;

    with_cx(|cx| {
        for step in 1..=total_steps {
            let food_water = MultiphaseVesselComponent::liquid_stream_for_pure_species(
                0,
                0.002,
                373.15,
                Some(&canola),
            );
            let q_thermostat = 8_825.0 + 350.0 * (450.15 - state.t_liq_k);
            let inputs = comp
                .open_vent_inputs(
                    &state,
                    dt_s,
                    q_thermostat,
                    &[-2_000.0, -2_000.0],
                    food_water,
                    [0.0; VAPOR_STREAM_WIDTH],
                    0.0,
                    101_325.0,
                    OpenVentSpecV1::open_fryer(),
                )
                .expect("fryer inputs");

            let mut next = state.clone();
            let Proposal::Step(receipt) = comp
                .propose(&state, dt_s, &inputs, &mut next, &mut out, cx)
                .expect("fryer step")
            else {
                panic!("unexpected shorten at fryer step {step}");
            };

            if step <= 5 || step.is_multiple_of(240) || step == total_steps {
                audit_full_step(&comp, &state, &next, dt_s, &inputs, &out, &receipt);
            }

            if step == 240 {
                t_cz_after_4h = Some((next.t_liq_k, next.t_coldzone_k));
            }
            state = next;
        }
    });

    let (t_act_4h, t_cz_4h) = t_cz_after_4h.expect("4h");
    let cz_drop_k = t_act_4h - t_cz_4h;
    assert!(
        (75.0..=125.0).contains(&cz_drop_k),
        "cold-zone drop after 4h should be ~83-111 K (150-200 F), got {cz_drop_k:.1} K"
    );

    let probes = &out[comp.layout().output_range(comp.probes_port())];
    let aged_ffa_pct = probes[PROBE_FFA_PCT];
    let aged_tpm_pct = probes[PROBE_TPM_PCT];
    let aged_smoke_k = probes[PROBE_T_SMOKE_K];
    let aged_visc = probes[PROBE_VISC_PA_S];

    assert!(
        (1.2..=5.0).contains(&aged_ffa_pct),
        "40h fryer FFA% out of expected range: {aged_ffa_pct:.2}%"
    );
    assert!(
        (20.0..=55.0).contains(&aged_tpm_pct),
        "40h fryer TPM% out of expected range: {aged_tpm_pct:.2}%"
    );
    assert!(
        fresh_smoke - aged_smoke_k >= 35.0,
        "expected >= 35 K smoke-point drop, got {:.1} K",
        fresh_smoke - aged_smoke_k
    );
    assert!(
        aged_visc / fresh_visc > 1.6,
        "expected significant TPM viscosity rise, ratio = {:.2}x",
        aged_visc / fresh_visc
    );
    assert!(
        state.m_glycerol_kg > 0.0,
        "hydrolysis must produce glycerol"
    );
    assert!(
        probes[PROBE_T_COLDZONE_K] < probes[PROBE_T_LIQ_K],
        "cold zone probe must remain cooler than active vat"
    );

    let h_still = comp.effective_immersed_htc_w_m2_k(&state, 0.0, 0.0);
    let h_bubbling = comp.effective_immersed_htc_w_m2_k(&state, 0.0, 0.025);
    assert!(
        h_bubbling > 2.5 * h_still,
        "fryer steam bubbling must boost immersed HTC (>2.5x), got {:.2}x",
        h_bubbling / h_still
    );

    let mut hot_flash = state.clone();
    hot_flash.t_liq_k = canola.t_flash_k + 2.0;
    let flags_flash = comp.safety_flags(&hot_flash, 0.0, 0.0);
    assert!(
        flags_flash.smoking
            && flags_flash.flashpoint_exceeded
            && !flags_flash.autoignition_exceeded
    );

    let mut hot_auto = state;
    hot_auto.t_liq_k = canola.t_autoignition_k + 5.0;
    let flags_auto = comp.safety_flags(&hot_auto, 0.0, 0.0);
    assert!(flags_auto.flashpoint_exceeded && flags_auto.autoignition_exceeded);
}

#[test]
fn planted_negative_hydrolysis_stoichiometry_fault_is_refused_by_audit_carried_reactions() {
    let comp = MultiphaseVesselComponent::new(MultiphaseVesselConfig::deep_fryer_default(
        OilCardV1::canola_refined(),
    ))
    .expect("fryer");
    let layout = comp.layout();

    let mut corrupted = vec![Some(0.0); layout.carried_len()];
    corrupted[0] = Some(-3.0 * 1.0e-5 * fs_multiphase::MOLAR_MASS_H2O_KG_MOL);
    corrupted[2] = Some(-1.0e-5 * fs_multiphase::MOLAR_MASS_TRIGLYCERIDE_KG_MOL);
    corrupted[3] = Some(3.0 * 1.0e-5 * fs_multiphase::MOLAR_MASS_FFA_KG_MOL);
    corrupted[4] = Some(0.0); // Missing glycerol breaks C, H, O and total mass closure.

    let err = audit_carried_reactions(layout, &corrupted, comp.reaction_basis())
        .expect_err("corrupted oil hydrolysis stoichiometry must be refused");
    assert!(
        matches!(
            err,
            ContractViolation::ElementalImbalance { .. }
                | ContractViolation::ReactionMassImbalance { .. }
        ),
        "expected elemental or mass imbalance violation, got {err:?}"
    );
}

#[test]
fn dryout_to_zero_stress_test_emits_shorten_and_zero_nans_or_negative_mass() {
    let mut cfg = MultiphaseVesselConfig::skillet_default();
    // Thin interior skillet thermal skin (120 J/K) keeps temperature finite even at m_liq == 0.
    cfg.dry_vessel_heat_capacity_j_k = 120.0;
    let comp = MultiphaseVesselComponent::new(cfg).expect("skillet");
    let mut state = comp
        .initial_state(VesselInitialStateV1::water(0.025, 372.0))
        .expect("state");
    let mut out = vec![0.0; comp.layout().output_len()];

    let mut shorten_count = 0usize;
    let mut dt = 5.0;
    let mut steps = 0usize;

    with_cx(|cx| {
        while state.m_h2o_liq_kg > 5.0e-5 && steps < 400 {
            steps += 1;
            let f_wet = comp.wetting_fraction(&state);
            let q_wall =
                1_800.0 * f_wet.max(0.05) * ((460.0 - state.t_liq_k) / 87.0).clamp(0.0, 1.2);
            let inputs = comp
                .open_vent_inputs(
                    &state,
                    dt,
                    q_wall,
                    &[0.0],
                    [0.0; LIQUID_STREAM_WIDTH],
                    [0.0; VAPOR_STREAM_WIDTH],
                    0.0,
                    101_325.0,
                    OpenVentSpecV1::open_pot(),
                )
                .expect("inputs");
            let mut next = state.clone();
            match comp
                .propose(&state, dt, &inputs, &mut next, &mut out, cx)
                .expect("propose must not panic or NaN on dryout")
            {
                Proposal::Shorten { dt_s } => {
                    assert!(dt_s > 0.0 && dt_s < dt);
                    shorten_count += 1;
                    dt = dt_s;
                }
                Proposal::Step(receipt) => {
                    audit_full_step(&comp, &state, &next, dt, &inputs, &out, &receipt);
                    assert!(next.m_h2o_liq_kg >= 0.0 && next.m_h2o_liq_kg.is_finite());
                    assert!(next.t_liq_k.is_finite() && next.t_liq_k > 0.0 && next.t_liq_k < 500.0);
                    assert!(out.iter().all(|v| v.is_finite()));
                    state = next;
                }
            }
        }
    });

    assert!(
        shorten_count >= 1,
        "expected Proposal::Shorten on dryout approach"
    );
    assert!(
        state.m_h2o_liq_kg < 1.0e-4,
        "expected water to boil down near zero, remaining = {}",
        state.m_h2o_liq_kg
    );
    assert!(comp.wetting_fraction(&state) < 0.05);
}

#[test]
fn altitude_denver_boiling_drop_pressure_cooker_elevation_and_sauce_reduction_rheology() {
    let comp =
        MultiphaseVesselComponent::new(MultiphaseVesselConfig::saucepan_default()).expect("pot");
    let water = comp
        .initial_state(VesselInitialStateV1::water(2.0, 350.0))
        .expect("water");

    let (t_sea, _, _) = comp
        .bubble_point_and_vapor_composition(&water, 101_325.0)
        .expect("sea");
    let (t_denver, _, _) = comp
        .bubble_point_and_vapor_composition(&water, 83_400.0)
        .expect("denver");
    let (t_pressure_12psig, _, _) = comp
        .bubble_point_and_vapor_composition(&water, 101_325.0 + 82_737.0)
        .expect("12 psig");

    assert!(((t_sea - 273.15) - 99.97).abs() < 0.15);
    assert!(
        ((t_denver - 273.15) - 94.6).abs() < 0.35,
        "Denver 83.4 kPa T_sat must be ~94.6 C, got {:.2} C",
        t_denver - 273.15
    );
    assert!(
        ((t_pressure_12psig - 273.15) - 117.7).abs() < 0.6,
        "12 psig pressure cooker T_sat must be ~117.7 C, got {:.2} C",
        t_pressure_12psig - 273.15
    );

    let dilute_sauce = comp
        .initial_state(VesselInitialStateV1::sauce(
            1.0, 0.04, 0.02, 0.01, 0.02, 370.0,
        ))
        .expect("dilute");
    let reduced_sauce = comp
        .initial_state(VesselInitialStateV1::sauce(
            0.22, 0.04, 0.02, 0.01, 0.02, 370.0,
        ))
        .expect("reduced");

    let (_, _, bpe_dil) = comp
        .bubble_point_and_vapor_composition(&dilute_sauce, 101_325.0)
        .expect("bpe_dil");
    let (_, _, bpe_red) = comp
        .bubble_point_and_vapor_composition(&reduced_sauce, 101_325.0)
        .expect("bpe_red");
    assert!(bpe_red > 3.0 * bpe_dil, "reduced sauce must elevate BPE");

    let mu_dil = comp.effective_viscosity_pa_s(&dilute_sauce, 0.0);
    let mu_red_still = comp.effective_viscosity_pa_s(&reduced_sauce, 0.0);
    let mu_red_stirred = comp.effective_viscosity_pa_s(&reduced_sauce, 120.0);
    assert!(mu_red_still > 5.0 * mu_dil, "reduction must thicken sauce");
    assert!(
        mu_red_stirred < mu_red_still,
        "Carreau shear-thinning must reduce viscosity under stirring"
    );

    let flags_still = comp.safety_flags(&reduced_sauce, 650.0, 0.0);
    let flags_stirred = comp.safety_flags(&reduced_sauce, 650.0, 150.0);
    assert!(
        flags_still.scorch_risk && !flags_stirred.scorch_risk,
        "stopping stirring on reduced thick sauce must raise scorch risk while stirring prevents it"
    );

    with_cx(|cx| {
        let inputs = comp
            .build_inputs(
                300.0,
                &[0.0],
                [0.0; LIQUID_STREAM_WIDTH],
                [0.0; VAPOR_STREAM_WIDTH],
                60.0,
                83_400.0,
            )
            .expect("inputs");
        let mut next = reduced_sauce.clone();
        let mut out = vec![0.0; comp.layout().output_len()];
        let Proposal::Step(receipt) = comp
            .propose(&reduced_sauce, 1.0, &inputs, &mut next, &mut out, cx)
            .expect("step")
        else {
            panic!("unexpected shorten");
        };
        audit_full_step(&comp, &reduced_sauce, &next, 1.0, &inputs, &out, &receipt);
        let probes = &out[comp.layout().output_range(comp.probes_port())];
        assert!(probes[PROBE_X_SOLIDS] > 0.20);
        assert!(probes[PROBE_M_LIQ_KG] > 0.25);
        assert!(probes[PROBE_T_LIQ_K] > 365.0);
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn oil_card_loads_from_fs_matdb_and_refuses_planted_invalid_cards_and_inputs() {
    let mut params = BTreeMap::new();
    let ins = |map: &mut BTreeMap<String, LawParameter>, k: &str, v: f64, d: [i8; 6]| {
        map.insert(
            k.to_owned(),
            LawParameter {
                value: v,
                dims: Dims(d),
            },
        );
    };
    ins(
        &mut params,
        "lipid_mass_fraction",
        0.999,
        [0, 0, 0, 0, 0, 0],
    );
    ins(
        &mut params,
        "water_mass_fraction",
        0.001,
        [0, 0, 0, 0, 0, 0],
    );
    ins(
        &mut params,
        "milk_solids_mass_fraction",
        0.0,
        [0, 0, 0, 0, 0, 0],
    );
    ins(
        &mut params,
        "salt_nacl_mass_fraction",
        0.0,
        [0, 0, 0, 0, 0, 0],
    );
    ins(&mut params, "t_ref", 298.15, [0, 0, 0, 1, 0, 0]);
    ins(&mut params, "rho_ref", 914.0, [-3, 1, 0, 0, 0, 0]);
    ins(&mut params, "drho_dt", -0.66, [-3, 1, 0, -1, 0, 0]);
    ins(&mut params, "cp_ref", 1985.0, [2, 0, -2, -1, 0, 0]);
    ins(&mut params, "dcp_dt", 2.25, [2, 0, -2, -2, 0, 0]);
    ins(&mut params, "k_ref", 0.168, [1, 1, -3, -1, 0, 0]);
    ins(&mut params, "dk_dt", -0.00021, [1, 1, -3, -2, 0, 0]);
    ins(&mut params, "mu_0", 2.4e-6, [-1, 1, -1, 0, 0, 0]);
    ins(&mut params, "e_visc", 24900.0, [2, 1, -2, 0, 0, -1]);
    ins(&mut params, "sigma_ref", 0.032, [0, 1, -2, 0, 0, 0]);
    ins(&mut params, "dsigma_dt", -6.7e-5, [0, 1, -2, -1, 0, 0]);
    ins(&mut params, "t_smoke_fresh", 477.15, [0, 0, 0, 1, 0, 0]);
    ins(&mut params, "t_flash", 588.15, [0, 0, 0, 1, 0, 0]);
    ins(&mut params, "t_autoignition", 651.15, [0, 0, 0, 1, 0, 0]);
    ins(&mut params, "k_hyd", 110000.0, [0, 0, -1, 0, 0, 0]);
    ins(&mut params, "ea_hyd", 84000.0, [2, 1, -2, 0, 0, -1]);
    ins(&mut params, "k_ox", 360000.0, [0, 0, -1, 0, 0, 0]);
    ins(&mut params, "ea_ox", 88000.0, [2, 1, -2, 0, 0, -1]);
    ins(&mut params, "k_poly", 200000.0, [0, 0, -1, 0, 0, 0]);
    ins(&mut params, "ea_poly", 94000.0, [2, 1, -2, 0, 0, -1]);
    ins(&mut params, "beta_ffa", -1250.0, [0, 0, 0, 1, 0, 0]);
    ins(&mut params, "beta_tpm", -120.0, [0, 0, 0, 1, 0, 0]);
    ins(&mut params, "gamma_tpm", 2.6, [0, 0, 0, 0, 0, 0]);
    ins(&mut params, "ffa_fresh", 0.0005, [0, 0, 0, 0, 0, 0]);
    ins(&mut params, "tpm_discard_limit", 0.25, [0, 0, 0, 0, 0, 0]);

    let card = ConstitutiveModelCard {
        law: LawId(COOKING_OIL_LAW_ID_V1.to_owned()),
        law_version: 1,
        parameters: params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 273.15, 588.15),
        sources: Vec::new(),
        provenance: Provenance {
            source: "data/matdb/kitchen-v1/oils/oil-canola-refined".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    let loaded = OilCardV1::from_cooking_oil_card(&card).expect("valid canola card");
    assert!((loaded.t_smoke_fresh_k - 477.15).abs() < 1.0e-12);
    assert!((loaded.beta_smoke_ffa_k - 1250.0).abs() < 1.0e-12);

    // Planted negative 1: wrong law id refused.
    let mut bad_law = card;
    bad_law.law = LawId("wrong-oil-law".to_owned());
    assert!(matches!(
        OilCardV1::from_cooking_oil_card(&bad_law),
        Err(VesselErrorV1::InvalidCard { .. })
    ));

    // Planted negative 2: inconsistent stream mass flow refused.
    let comp =
        MultiphaseVesselComponent::new(MultiphaseVesselConfig::saucepan_default()).expect("v");
    let state = comp
        .initial_state(VesselInitialStateV1::water(1.0, 350.0))
        .expect("s");
    let mut bad_liq = [0.0; LIQUID_STREAM_WIDTH];
    bad_liq[STREAM_MASS_RATE] = 0.05;
    let bad_inputs = comp
        .build_inputs(
            0.0,
            &[0.0],
            bad_liq,
            [0.0; VAPOR_STREAM_WIDTH],
            0.0,
            101_325.0,
        )
        .expect("inputs");
    with_cx(|cx| {
        let mut next = state.clone();
        let mut out = vec![0.0; comp.layout().output_len()];
        assert!(matches!(
            comp.propose(&state, 1.0, &bad_inputs, &mut next, &mut out, cx),
            Err(VesselErrorV1::InconsistentStreamMassFlow { .. })
        ));

        // Planted negative 3: negative stir_rpm refused.
        let neg_stir = comp
            .build_inputs(
                0.0,
                &[0.0],
                [0.0; LIQUID_STREAM_WIDTH],
                [0.0; VAPOR_STREAM_WIDTH],
                -10.0,
                101_325.0,
            )
            .expect("inputs");
        assert!(matches!(
            comp.propose(&state, 1.0, &neg_stir, &mut next, &mut out, cx),
            Err(VesselErrorV1::InvalidScalar { .. })
        ));
    });
}
