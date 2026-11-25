//! G0/G2/G3 verification oracles, `fs-ports` conformance checks, frosting/defrost
//! transients, condenser-fouling diagnostics, and air-cooled ice harvester tests
//! for `fs_multiphase::cycle` (`KM2.3`).

use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_multiphase::{
    CompressorStateV1, CondenserEvaporatorHXComponent, CondenserEvaporatorHXConfig, CycleErrorV1,
    CycleOperatingPointV1, ExpansionThrottleComponent, ExpansionThrottleConfig,
    ExpansionThrottleStateV1, IceHarvesterComponent, IceHarvesterConfig, PR_REFRIGERANT_LAW_ID_V1,
    RefrigerantModelV1, RefrigerationCompressorComponent, RefrigerationCompressorConfig,
};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{
    Component, Proposal, RefrigerantSpecies, audit_carried_step, audit_mass_step, audit_step,
};
use fs_qty::Dims;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 23,
                kernel_id: 333,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn verify_component_step<C>(
    comp: &C,
    old: &C::State,
    new: &C::State,
    dt_s: f64,
    inputs: &[f64],
    outputs: &[f64],
    receipt: &fs_ports::component::Receipt,
) where
    C: Component,
{
    let ports = comp.ports();
    let layout = fs_ports::component::PortLayout::new(ports).expect("layout");
    let e_old = comp.stored_energy_j(old).expect("e_old");
    let e_new = comp.stored_energy_j(new).expect("e_new");
    let audit = audit_step(ports, &layout, dt_s, inputs, outputs, receipt, e_old, e_new)
        .expect("First-Law audit_step must pass");
    let e_scale = e_old.abs().max(e_new.abs()).max(1.0);
    assert!(
        audit.balance_defect_j.abs() <= 1.0e-11 * e_scale + 1.0e-9,
        "energy balance defect {} exceeded tolerance (scale {})",
        audit.balance_defect_j,
        e_scale
    );

    let m_old = comp.stored_mass_kg(old).expect("m_old");
    let m_new = comp.stored_mass_kg(new).expect("m_new");
    audit_mass_step(ports, &layout, dt_s, inputs, outputs, m_old, m_new)
        .expect("mass conservation audit_mass_step must pass");

    for c in 0..layout.carried_len() {
        let c_old = comp.stored_carried_quantity(old, c).expect("c_old");
        let c_new = comp.stored_carried_quantity(new, c).expect("c_new");
        audit_carried_step(ports, &layout, dt_s, inputs, outputs, c, c_old, c_new)
            .unwrap_or_else(|e| panic!("carried coordinate {c} failed conservation: {e}"));
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn g0_first_and_second_law_across_all_five_refrigerants_including_r448a_glide_and_r744_transcritical()
 {
    // 1. Subcritical walk-in cooler operating point across R-290, R-134a, R-404A, R-448A.
    let subcrit_spec = CycleOperatingPointV1 {
        t_evap_k: 268.15, // -5 C evaporator
        suction_superheat_k: 6.0,
        t_ambient_k: 303.15, // 30 C kitchen ambient
        ua_cond_clean_w_k: 260.0,
        condenser_fouling_factor: 0.0,
        subcooling_k: 4.0,
        swept_volume_m3_s: 1.1e-3,
        eta_isentropic: 0.72,
        eta_motor: 0.88,
        shell_loss_fraction: 0.85,
        clearance_ratio: 0.045,
        supercritical_pressure_pa: None,
    };

    for species in [
        RefrigerantSpecies::R290,
        RefrigerantSpecies::R134a,
        RefrigerantSpecies::R404A,
        RefrigerantSpecies::R448A,
    ] {
        let model = RefrigerantModelV1::for_species(species).expect("model");
        let sum = model.solve_closed_cycle(subcrit_spec).expect("cycle");
        assert!(
            !sum.is_supercritical_gas_cooler,
            "{species:?} should operate subcritically at 30 C"
        );
        assert!(
            sum.first_law_relative_defect < 1.0e-11,
            "{species:?} First Law defect {}",
            sum.first_law_relative_defect
        );
        assert!(
            sum.cop_cooling > 1.5 && sum.cop_cooling < sum.cop_carnot,
            "{species:?} COP={:.3} must be positive and strictly < Carnot={:.3}",
            sum.cop_cooling,
            sum.cop_carnot
        );
        assert!(
            sum.internal_entropy_gen_w_k > 0.0,
            "{species:?} must generate positive internal entropy"
        );
        assert!(
            (sum.throttle.outlet_enthalpy_j_kg - sum.throttle.inlet.enthalpy_j_kg).abs() < 1.0e-12,
            "{species:?} TXV must be exactly isenthalpic"
        );
        if species == RefrigerantSpecies::R448A {
            assert!(
                (5.0..=7.2).contains(&sum.evap_envelope.glide_k),
                "R-448A evaporator glide expected ~5.5-6.7 K, got {:.2} K",
                sum.evap_envelope.glide_k
            );
        }
    }

    // 2. R-744 (CO2): subcritical (15 C cool ambient), supercritical 95 F (35 C kitchen),
    //    and supercritical 115 F (46.1 C Phoenix ambient) with declared 18% near-critical band!
    let r744 = RefrigerantModelV1::r744().expect("r744");
    let r744_sub = r744
        .solve_closed_cycle(CycleOperatingPointV1 {
            t_evap_k: 263.15,    // -10 C
            t_ambient_k: 288.15, // 15 C cool ambient -> subcritical condenser
            ua_cond_clean_w_k: 420.0,
            subcooling_k: 3.0,
            swept_volume_m3_s: 2.2e-4,
            ..subcrit_spec
        })
        .expect("R-744 subcritical");
    assert!(!r744_sub.is_supercritical_gas_cooler);
    assert!(r744_sub.first_law_relative_defect < 1.0e-11);
    assert!(r744_sub.cop_cooling > 1.2 && r744_sub.cop_cooling < r744_sub.cop_carnot);

    let r744_95f = r744
        .solve_closed_cycle(CycleOperatingPointV1 {
            t_evap_k: 263.15,    // -10 C
            t_ambient_k: 308.15, // 95 F (35 C) kitchen -> supercritical gas-cooler (> 30.98 C)
            swept_volume_m3_s: 2.2e-4,
            supercritical_pressure_pa: Some(9_500_000.0), // 95 bar
            ..subcrit_spec
        })
        .expect("R-744 95F supercritical");
    assert!(r744_95f.is_supercritical_gas_cooler);
    assert!(r744_95f.near_critical_warning);
    assert!((r744_95f.near_critical_density_error_band - 0.18).abs() < 1.0e-12);
    assert!(r744_95f.first_law_relative_defect < 1.0e-11);
    assert!(r744_95f.cop_cooling > 0.8 && r744_95f.cop_cooling < r744_95f.cop_carnot);
    assert!(r744_95f.internal_entropy_gen_w_k > 0.0);

    let r744_115f = r744
        .solve_closed_cycle(CycleOperatingPointV1 {
            t_evap_k: 263.15,    // -10 C
            t_ambient_k: 319.26, // 115 F (46.1 C) Phoenix ambient
            swept_volume_m3_s: 2.2e-4,
            supercritical_pressure_pa: Some(11_500_000.0), // 115 bar
            ..subcrit_spec
        })
        .expect("R-744 115F Phoenix supercritical");
    assert!(r744_115f.is_supercritical_gas_cooler);
    assert!(r744_115f.first_law_relative_defect < 1.0e-11);
    assert!(
        r744_115f.cop_cooling > 0.5
            && r744_115f.cop_cooling < r744_95f.cop_cooling
            && r744_115f.cop_cooling < r744_115f.cop_carnot
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn fs_ports_conformance_and_durable_codec_roundtrips_for_all_four_components() {
    with_cx(|cx| {
        // 1. RefrigerationCompressorComponent
        let comp = RefrigerationCompressorComponent::new(
            RefrigerationCompressorConfig::commercial_default(RefrigerantSpecies::R290),
        )
        .expect("comp");
        let p_suc = comp.fluid().bubble_pressure_at_t(268.15).expect("ps");
        let p_dis = comp.fluid().bubble_pressure_at_t(315.15).expect("pd");
        let suc_pt = comp
            .fluid()
            .state_pt(274.15, p_suc, PengRobinsonPhaseV1::Vapor)
            .expect("suc");
        let c_state = CompressorStateV1::new(310.0).expect("s");
        let c_stale = CompressorStateV1::new(340.0).expect("stale");
        let c_in = comp.build_inputs(
            230.0,
            p_suc,
            suc_pt.enthalpy_j_kg,
            p_dis,
            suc_pt.enthalpy_j_kg + 80_000.0,
            298.15,
            0.0,
        );
        let obs = check_component(
            &comp,
            &Case {
                state: &c_state,
                stale: &c_stale,
                dt_s: 5.0,
                inputs: &c_in,
            },
            cx,
        )
        .expect("compressor check_component");
        assert!(matches!(obs, Observed::Stepped { .. }));
        let mut words = Vec::new();
        assert!(comp.encode_state(&c_state, &mut words));
        assert_eq!(comp.decode_state(&words), Some(c_state));
        assert_eq!(
            RefrigerationCompressorComponent::STATE_CHARTER
                .expect("charter")
                .owner,
            "fs-multiphase::RefrigerationCompressor"
        );

        // 2. CondenserEvaporatorHXComponent (both CondenserGasCooler and FrostingEvaporator modes)
        let cond = CondenserEvaporatorHXComponent::new(
            CondenserEvaporatorHXConfig::condenser_default(RefrigerantSpecies::R744),
        )
        .expect("cond");
        let r744_gc = RefrigerantModelV1::r744()
            .expect("r744")
            .state_pt(380.0, 9_500_000.0, PengRobinsonPhaseV1::Vapor)
            .expect("gc");
        let h_state = cond.initial_state(325.0).expect("hs");
        let h_stale = cond.initial_state(310.0).expect("hst");
        let h_in = cond
            .build_inputs(
                &h_state,
                5.0,
                9_500_000.0,
                r744_gc.enthalpy_j_kg,
                0.025,
                308.15,
                0.0,
                0.0,
                0.0,
                0.0,
            )
            .expect("h_in");
        assert!(matches!(
            check_component(
                &cond,
                &Case {
                    state: &h_state,
                    stale: &h_stale,
                    dt_s: 5.0,
                    inputs: &h_in,
                },
                cx,
            )
            .expect("condenser check_component"),
            Observed::Stepped { .. }
        ));
        words.clear();
        assert!(cond.encode_state(&h_state, &mut words));
        assert_eq!(cond.decode_state(&words), Some(h_state));

        let evap = CondenserEvaporatorHXComponent::new(
            CondenserEvaporatorHXConfig::evaporator_default(RefrigerantSpecies::R448A),
        )
        .expect("evap");
        let e_state = evap.initial_state(265.0).expect("es");
        let e_stale = evap.initial_state(270.0).expect("est");
        let e_in = evap
            .build_inputs(
                &e_state, 5.0, 450_000.0, 240_000.0, 0.015, 276.0, 0.0, 1.0, 2.0e-4, 0.0,
            )
            .expect("e_in");
        assert!(matches!(
            check_component(
                &evap,
                &Case {
                    state: &e_state,
                    stale: &e_stale,
                    dt_s: 5.0,
                    inputs: &e_in,
                },
                cx,
            )
            .expect("evaporator check_component"),
            Observed::Stepped { .. }
        ));

        // 3. ExpansionThrottleComponent
        let txv = ExpansionThrottleComponent::new(ExpansionThrottleConfig::txv_default(
            RefrigerantSpecies::R134a,
        ))
        .expect("txv");
        let r134a = RefrigerantModelV1::r134a().expect("r134a");
        let p_hi = r134a.bubble_pressure_at_t(315.15).expect("phi");
        let p_lo = r134a.bubble_pressure_at_t(268.15).expect("plo");
        let liq_in = r134a
            .state_pt(310.15, p_hi, PengRobinsonPhaseV1::Liquid)
            .expect("liq");
        let t_state = ExpansionThrottleStateV1::new(0.45).expect("ts");
        let t_stale = ExpansionThrottleStateV1::new(0.80).expect("tst");
        let t_in = txv.build_inputs(p_hi, liq_in.enthalpy_j_kg, p_lo, 8.5, 310.15);
        assert!(matches!(
            check_component(
                &txv,
                &Case {
                    state: &t_state,
                    stale: &t_stale,
                    dt_s: 2.0,
                    inputs: &t_in,
                },
                cx,
            )
            .expect("txv check_component"),
            Observed::Stepped { .. }
        ));
        words.clear();
        assert!(txv.encode_state(&t_state, &mut words));
        assert_eq!(txv.decode_state(&words), Some(t_state));

        // 4. IceHarvesterComponent
        let ice = IceHarvesterComponent::new(IceHarvesterConfig::commercial_air_cooled_default())
            .expect("ice");
        let i_state = ice.initial_state(5.0, 12.0, 283.15).expect("is");
        let i_stale = ice.initial_state(20.0, 8.0, 290.0).expect("ist");
        let i_in = ice.build_inputs(208.0, 294.26, 0.0055, 283.15, 0.0);
        assert!(matches!(
            check_component(
                &ice,
                &Case {
                    state: &i_state,
                    stale: &i_stale,
                    dt_s: 30.0,
                    inputs: &i_in,
                },
                cx,
            )
            .expect("ice harvester check_component"),
            Observed::Stepped { .. }
        ));
        words.clear();
        assert!(ice.encode_state(&i_state, &mut words));
        assert_eq!(ice.decode_state(&words), Some(i_state));
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn evaporator_frosting_degrades_ua_and_defrost_drains_meltwater_restoring_capacity() {
    let evap = CondenserEvaporatorHXComponent::new(
        CondenserEvaporatorHXConfig::evaporator_default(RefrigerantSpecies::R448A),
    )
    .expect("evap");
    let mut state = evap.initial_state(264.0).expect("init");
    let ua_clean = evap.effective_ua_w_k(&state, 1.0);
    let mut out = vec![0.0; evap.layout().output_len()];
    let dt_s = 60.0;

    with_cx(|cx| {
        // Phase 1: 4 hours of humid walk-in door traffic (0.25 g/s frost deposition onto coil).
        let mut prev_ua = ua_clean;
        for _ in 0..240 {
            let inputs = evap
                .build_inputs(
                    &state, dt_s, 430_000.0, 245_000.0, 0.040, 275.15, // 2 C walk-in box air
                    0.0,    // defrost off
                    1.0,    // full fan airflow
                    2.5e-4, // humid moisture frosting rate [kg/s]
                    0.0,
                )
                .expect("inputs");
            let mut next = state.clone();
            let Proposal::Step(receipt) = evap
                .propose(&state, dt_s, &inputs, &mut next, &mut out, cx)
                .expect("frost step")
            else {
                panic!("unexpected shorten");
            };
            verify_component_step(&evap, &state, &next, dt_s, &inputs, &out, &receipt);
            let ua_now = evap.effective_ua_w_k(&next, 1.0);
            assert!(
                ua_now <= prev_ua + 1.0e-12,
                "UA must degrade monotonically during frosting: {ua_now} vs {prev_ua}"
            );
            prev_ua = ua_now;
            state = next;
        }

        let frosted_mass = state.m_frost_kg;
        let ua_frosted = evap.effective_ua_w_k(&state, 1.0);
        assert!(
            frosted_mass > 3.0,
            "expected >3 kg frost accumulated, got {frosted_mass}"
        );
        assert!(
            ua_frosted < 0.70 * ua_clean,
            "frost layer must materially degrade UA_evap (<70% clean): {ua_frosted:.1} vs {ua_clean:.1}"
        );

        // Phase 2: 20-minute 2.2 kW electric defrost cycle draining liquid meltwater OUT via drain_stream!
        let mut saw_box_warming_spike = false;
        for _ in 0..20 {
            let inputs = evap
                .build_inputs(
                    &state, dt_s, 430_000.0, 245_000.0, 0.0, // compressor off during defrost
                    275.15, 2_200.0, // electric defrost heater
                    1.0, 0.0, 1.0, // drain 100% of pan meltwater out via drain_stream!
                )
                .expect("defrost inputs");
            let mut next = state.clone();
            let Proposal::Step(receipt) = evap
                .propose(&state, dt_s, &inputs, &mut next, &mut out, cx)
                .expect("defrost step")
            else {
                panic!("unexpected shorten");
            };
            verify_component_step(&evap, &state, &next, dt_s, &inputs, &out, &receipt);
            let probes = &out[evap.layout().output_range(evap.probes_port())];
            if probes[9] > 10.0 {
                saw_box_warming_spike = true;
            }
            state = next;
        }

        // Flush final meltwater pan step so all melted water exits drain_stream.
        let flush_in = evap
            .build_inputs(
                &state, dt_s, 430_000.0, 245_000.0, 0.0, 275.15, 0.0, 1.0, 0.0, 1.0,
            )
            .expect("flush");
        let mut flushed = state.clone();
        let Proposal::Step(receipt) = evap
            .propose(&state, dt_s, &flush_in, &mut flushed, &mut out, cx)
            .expect("flush step")
        else {
            panic!("unexpected shorten");
        };
        verify_component_step(&evap, &state, &flushed, dt_s, &flush_in, &out, &receipt);
        state = flushed;

        let ua_restored = evap.effective_ua_w_k(&state, 1.0);
        assert!(
            state.m_frost_kg < 0.05,
            "defrost must clear frost layer, remaining = {} kg",
            state.m_frost_kg
        );
        assert!(
            (state.cumulative_drain_kg - frosted_mass).abs() < 0.08,
            "drained meltwater {:.3} kg must match melted frost {:.3} kg (never re-evaporated!)",
            state.cumulative_drain_kg,
            frosted_mass
        );
        assert!(
            ua_restored > 0.98 * ua_clean,
            "defrost must restore UA_evap to clean baseline: {ua_restored:.1} vs {ua_clean:.1}"
        );
        assert!(
            saw_box_warming_spike,
            "defrost must produce measurable temporary box-air warming spike"
        );
    });
}

#[test]
fn planted_25pct_condenser_fouling_raises_compressor_power_current_and_head_pressure() {
    // 1. Cycle-level verification across R-404A.
    let r404a = RefrigerantModelV1::r404a().expect("r404a");
    let clean_spec = CycleOperatingPointV1 {
        t_evap_k: 263.15,
        suction_superheat_k: 6.0,
        t_ambient_k: 298.15,
        ua_cond_clean_w_k: 220.0,
        condenser_fouling_factor: 0.0,
        subcooling_k: 4.0,
        swept_volume_m3_s: 1.15e-3,
        eta_isentropic: 0.72,
        eta_motor: 0.88,
        shell_loss_fraction: 0.85,
        clearance_ratio: 0.045,
        supercritical_pressure_pa: None,
    };
    let fouled_spec = CycleOperatingPointV1 {
        condenser_fouling_factor: 0.25, // +25% condenser fouling
        ..clean_spec
    };
    let clean = r404a.solve_closed_cycle(clean_spec).expect("clean");
    let fouled = r404a.solve_closed_cycle(fouled_spec).expect("fouled");

    assert!(
        fouled.p_high_pa > clean.p_high_pa * 1.04,
        "+25% fouling must raise head pressure noticeably: {:.2} vs {:.2} bar",
        fouled.p_high_pa / 1.0e5,
        clean.p_high_pa / 1.0e5
    );
    assert!(
        fouled.w_comp_w > clean.w_comp_w * 1.02,
        "+25% fouling must raise compressor power: {:.1} vs {:.1} W",
        fouled.w_comp_w,
        clean.w_comp_w
    );
    assert!(fouled.cop_cooling < clean.cop_cooling);

    // 2. Component-level compressor current & discharge pressure rise under +25% fouling.
    let comp = RefrigerationCompressorComponent::new(
        RefrigerationCompressorConfig::commercial_default(RefrigerantSpecies::R404A),
    )
    .expect("comp");
    let state = CompressorStateV1::new(315.0).expect("s");
    let in_clean = comp.build_inputs(
        208.0,
        clean.p_evap_pa,
        clean.compression.suction.enthalpy_j_kg,
        clean.p_high_pa,
        clean.compression.discharge.enthalpy_j_kg,
        298.15,
        0.0,
    );
    let in_foul = comp.build_inputs(
        208.0,
        clean.p_evap_pa,
        clean.compression.suction.enthalpy_j_kg,
        clean.p_high_pa,
        clean.compression.discharge.enthalpy_j_kg,
        298.15,
        0.25,
    );
    with_cx(|cx| {
        let mut s_c = state.clone();
        let mut o_c = vec![0.0; comp.layout().output_len()];
        let mut s_f = state.clone();
        let mut o_f = vec![0.0; comp.layout().output_len()];
        comp.propose(&state, 10.0, &in_clean, &mut s_c, &mut o_c, cx)
            .expect("clean propose");
        comp.propose(&state, 10.0, &in_foul, &mut s_f, &mut o_f, cx)
            .expect("foul propose");
        let i_clean = o_c[comp.layout().output_range(0)][0];
        let i_foul = o_f[comp.layout().output_range(0)][0];
        let pr_clean = o_c[comp.layout().output_range(5)][4];
        let pr_foul = o_f[comp.layout().output_range(5)][4];
        assert!(
            i_foul > i_clean * 1.03,
            "fouling must elevate compressor current: {i_foul:.3} A vs {i_clean:.3} A"
        );
        assert!(pr_foul > pr_clean * 1.05);
    });
}

#[test]
fn air_cooled_ice_harvester_reproduces_25_to_40pct_ambient_derating_waste_heat_and_bin_conservation()
 {
    let harvester = IceHarvesterComponent::new(IceHarvesterConfig::commercial_air_cooled_default())
        .expect("ice");

    // AHRI-810 standard vs hot kitchen ambient:
    // 70 F (294.26 K) air / 50 F (283.15 K) water  vs  95 F (308.15 K) air / 70 F (294.26 K) water.
    let run_one_hour = |t_room_k: f64, t_water_k: f64| -> (f64, f64, f64) {
        let mut state = harvester.initial_state(0.0, 25.0, t_water_k).expect("init");
        let mut out = vec![0.0; harvester.layout().output_len()];
        let dt_s = 60.0;
        let mut total_room_w_sum = 0.0;
        let feed_rate = 0.0055; // steady potable water feed [kg/s]
        let initial_mass = state.total_mass_kg();
        with_cx(|cx| {
            for _ in 0..60 {
                let inputs = harvester.build_inputs(208.0, t_room_k, feed_rate, t_water_k, 0.0);
                let mut next = state.clone();
                let Proposal::Step(receipt) = harvester
                    .propose(&state, dt_s, &inputs, &mut next, &mut out, cx)
                    .expect("ice step")
                else {
                    panic!("unexpected shorten");
                };
                verify_component_step(&harvester, &state, &next, dt_s, &inputs, &out, &receipt);
                let q_room_w = -out[harvester.layout().output_range(1)][0];
                total_room_w_sum += q_room_w;
                state = next;
            }
        });
        // Exact mass conservation over 1 full hour: m_final - m_initial == 3600 * feed_rate!
        let expected_mass = initial_mass + 3_600.0 * feed_rate;
        assert!(
            (state.total_mass_kg() - expected_mass).abs() < 1.0e-11 * expected_mass,
            "1-hour ice harvester mass conservation violated"
        );
        let avg_room_heat_w = total_room_w_sum / 60.0;
        let probes = &out[harvester.layout().output_range(5)];
        (state.ice_bin_mass_kg, avg_room_heat_w, probes[5])
    };

    let (ice_70f, q_room_70f, cop_70f) = run_one_hour(294.26, 283.15);
    let (ice_95f, q_room_95f, cop_95f) = run_one_hour(308.15, 294.26);

    let derating_pct = 100.0 * (ice_70f - ice_95f) / ice_70f;
    assert!(
        (25.0..=40.0).contains(&derating_pct),
        "70 F vs 95 F ice harvest derating must lie in 25-40%, got {derating_pct:.1}% ({ice_70f:.2} -> {ice_95f:.2} kg/h)"
    );
    assert!(
        (1_500.0..=3_500.0).contains(&q_room_70f) && (1_500.0..=3_500.0).contains(&q_room_95f),
        "continuous condenser waste heat to kitchen air must lie in 1.5-3.5 kW: 70F={q_room_70f:.0} W, 95F={q_room_95f:.0} W"
    );
    assert!(
        cop_95f < cop_70f && cop_95f > 1.0,
        "hotter kitchen ambient must lower COP: 95F={cop_95f:.2} vs 70F={cop_70f:.2}"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn refrigerant_card_loads_from_fs_matdb_and_planted_negatives_refuse_cleanly() {
    // 1. Load R-290 directly from a kitchen-v1 peng-robinson-refrigerant-v1 card.
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
    ins(&mut params, "component_count", 1.0, [0, 0, 0, 0, 0, 0]);
    ins(
        &mut params,
        "mixture_critical_temperature",
        369.89,
        [0, 0, 0, 1, 0, 0],
    );
    ins(
        &mut params,
        "mixture_critical_pressure",
        4_251_200.0,
        [-1, 1, -2, 0, 0, 0],
    );
    ins(
        &mut params,
        "subcritical_accurate_t_max",
        351.4,
        [0, 0, 0, 1, 0, 0],
    );
    ins(
        &mut params,
        "near_critical_density_error_band",
        0.025,
        [0, 0, 0, 0, 0, 0],
    );
    ins(&mut params, "comp_0_mole_fraction", 1.0, [0, 0, 0, 0, 0, 0]);
    ins(
        &mut params,
        "comp_0_critical_temperature",
        369.89,
        [0, 0, 0, 1, 0, 0],
    );
    ins(
        &mut params,
        "comp_0_critical_pressure",
        4_251_200.0,
        [-1, 1, -2, 0, 0, 0],
    );
    ins(
        &mut params,
        "comp_0_acentric_factor",
        0.1521,
        [0, 0, 0, 0, 0, 0],
    );
    ins(
        &mut params,
        "comp_0_molar_mass",
        0.044_095_6,
        [0, 1, 0, 0, 0, -1],
    );
    ins(
        &mut params,
        "comp_0_peneloux_c",
        -4.739_282e-6,
        [3, 0, 0, 0, 0, -1],
    );

    let card = ConstitutiveModelCard {
        law: LawId(PR_REFRIGERANT_LAW_ID_V1.to_owned()),
        law_version: 1,
        parameters: params,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 210.0, 351.4),
        sources: Vec::new(),
        provenance: Provenance {
            source: "data/matdb/kitchen-v1/refrigerants/r290-propane".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    let loaded =
        RefrigerantModelV1::from_card(RefrigerantSpecies::R290, &card, 1_670.0).expect("load r290");
    assert!((loaded.critical_temperature_k() - 369.89).abs() < 1.0e-12);

    // 2. Planted negatives: negative pressure, negative temperature, superheat inversion,
    //    discharge < suction inversion, negative fouling factor -> typed errors, zero panics/NaNs!
    assert!(matches!(
        loaded.state_pt(300.0, -100_000.0, PengRobinsonPhaseV1::Vapor),
        Err(CycleErrorV1::InvalidPressure { .. })
    ));
    assert!(matches!(
        loaded.state_pt(-10.0, 500_000.0, PengRobinsonPhaseV1::Vapor),
        Err(CycleErrorV1::InvalidTemperature { .. })
    ));
    assert!(matches!(
        loaded.isentropic_compression(275.0, 1_000_000.0, 500_000.0, 0.72),
        Err(CycleErrorV1::ThermodynamicInversion(_))
    ));
    assert!(matches!(
        loaded.isenthalpic_throttle(310.0, 500_000.0, 1_000_000.0),
        Err(CycleErrorV1::ThermodynamicInversion(_))
    ));

    let comp = RefrigerationCompressorComponent::new(
        RefrigerationCompressorConfig::commercial_default(RefrigerantSpecies::R290),
    )
    .expect("comp");
    let c_state = CompressorStateV1::new(310.0).expect("s");
    with_cx(|cx| {
        let mut next = c_state.clone();
        let mut out = vec![0.0; comp.layout().output_len()];
        // Negative fouling factor refused.
        let bad_foul = comp.build_inputs(
            230.0,
            400_000.0,
            600_000.0,
            1_400_000.0,
            650_000.0,
            298.15,
            -0.1,
        );
        assert!(matches!(
            comp.propose(&c_state, 1.0, &bad_foul, &mut next, &mut out, cx),
            Err(CycleErrorV1::InvalidScalar { .. })
        ));
        // Negative superheat on TXV refused.
        let txv = ExpansionThrottleComponent::new(ExpansionThrottleConfig::txv_default(
            RefrigerantSpecies::R290,
        ))
        .expect("txv");
        let t_state = ExpansionThrottleStateV1::new(0.5).expect("ts");
        let mut t_next = t_state.clone();
        let mut t_out = vec![0.0; txv.layout().output_len()];
        let bad_sh = txv.build_inputs(1_400_000.0, 300_000.0, 400_000.0, -2.0, 310.0);
        assert!(matches!(
            txv.propose(&t_state, 1.0, &bad_sh, &mut t_next, &mut t_out, cx),
            Err(CycleErrorV1::ThermodynamicInversion(_))
        ));
    });
}
