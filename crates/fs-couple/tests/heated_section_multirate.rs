//! Equivalence and 100x multirate acceptance tests for the migrated
//! `examples/heated_section_sound.rs` composition (`apps/observatory`,
//! bead `frankensim-rc-1109-root-3gc12.5.9`).

#[allow(dead_code, clippy::all, clippy::pedantic)]
#[path = "../examples/heated_section_sound.rs"]
mod example;

use fs_conduction::transient::deposition::{Boundary, CellId, DepositionBudget, StepConfig};
use fs_couple::material_assembly::{
    AcousticMaterialBindings,
    concentric::thermal::heating::compile_heated_concentric_material_assembly,
};
use fs_evidence::Color;
use fs_material::conductor::OhmicDrive;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::section::thermal::AxialConstraint;
use fs_material::section_temperature::{
    ConcentricSectionTemperatureLaw,
    heating::{HeatedConcentricSection, SectionThermalInterface},
};
use fs_nlmodal::component::ModalAcousticState;
use fs_ports::component::Component;

struct PairingReference {
    temperatures_k: [f64; 2],
    electrical_energy_j: f64,
    stored_heat_j: f64,
    hot_tension_n: f64,
    layer_stresses_pa: Vec<f64>,
    hot_pa: Vec<f64>,
}

fn run_legacy_pairing(
    core: &fs_matdb::MaterialCard,
    sleeve: &fs_matdb::MaterialCard,
    charts: &[fs_material::phase::ResolvedHeatCapacityCurve; 2],
    reference: &fs_material::section::ConcentricElasticSection,
) -> PairingReference {
    let law =
        ConcentricSectionTemperatureLaw::try_new(reference.clone(), &[core, sleeve], &["T", "T"])
            .unwrap();
    let two_pi_l = 2.0 * core::f64::consts::PI * example::LENGTH;
    let core_node = example::RADII[0] / 2.0_f64.sqrt();
    let sleeve_node =
        ((example::RADII[0] * example::RADII[0] + example::RADII[1] * example::RADII[1]) / 2.0)
            .sqrt();
    let interface = SectionThermalInterface {
        area_m2: two_pi_l * example::RADII[0],
        resistance_geometry_per_m: [
            (example::RADII[0] / core_node).ln() / two_pi_l,
            (sleeve_node / example::RADII[0]).ln() / two_pi_l,
        ],
        contact_resistance_k_w: 0.0,
    };
    let boundary = Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    };
    let pairing_cold = HeatedConcentricSection::try_new(
        law,
        &[&charts[0], &charts[1]],
        example::LENGTH,
        &[boundary; 2],
        vec![interface],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 200_000,
        },
    )
    .unwrap()
    .with_axial_constraint(AxialConstraint::FixedStrain(1.0e-4), 0.01)
    .unwrap();
    let (pairing_hot, pairing_rep) = pairing_cold
        .advance_to(
            0.1,
            OhmicDrive::Voltage(0.2),
            StepConfig {
                max_dt_s: 0.001,
                energy_tolerance_j: 1.0e-7,
            },
            || Ok(()),
        )
        .unwrap();
    let axial = pairing_hot.axial_state().unwrap().response();
    let hot_pa = compile_heated_concentric_material_assembly(
        &example::template(),
        &pairing_hot,
        &AcousticMaterialBindings::default(),
    )
    .unwrap()
    .realize()
    .unwrap()
    .pressure_pa;
    let mut layer_stresses_pa = Vec::with_capacity(axial.layers().len());
    for l in axial.layers() {
        layer_stresses_pa.push(l.stress_pa());
    }
    PairingReference {
        temperatures_k: [
            pairing_hot.thermal().temperature_k(CellId(0)).unwrap(),
            pairing_hot.thermal().temperature_k(CellId(1)).unwrap(),
        ],
        electrical_energy_j: pairing_rep.electrical_energy_j,
        stored_heat_j: pairing_rep.thermal.stored_energy_change_j,
        hot_tension_n: axial.force_n(),
        layer_stresses_pa,
        hot_pa,
    }
}

#[test]
fn migrated_example_matches_pairing_and_closes_ledger() {
    let core = example::demo_with_expansion(1.0, Some(1.0e-5)).unwrap();
    let sleeve = example::demo_with_expansion(2.0, Some(2.0e-5)).unwrap();
    let charts = [
        example::caloric(&core).unwrap(),
        example::caloric(&sleeve).unwrap(),
    ];
    let reference = example::section([&core, &sleeve], true).unwrap();
    let pairing = run_legacy_pairing(&core, &sleeve, &charts, &reference);

    // Migrated example (`synthesize`) running the 100x multirate `fs_system::System`.
    let support = Some(example::FixedSupport::try_new(1.0e-4, 0.01).unwrap());
    let migrated = example::synthesize([&core, &sleeve], 0.2, 0.1, support).unwrap();

    for i in 0..2 {
        assert!(
            (migrated.temperatures_k[i] - pairing.temperatures_k[i]).abs() < 1.0e-4,
            "layer {i} T mismatch: migrated={} pairing={}",
            migrated.temperatures_k[i],
            pairing.temperatures_k[i]
        );
    }
    assert!(
        (migrated.electrical_energy_j - pairing.electrical_energy_j).abs()
            / pairing.electrical_energy_j
            < 1.0e-4
    );
    assert!(
        (migrated.stored_heat_j - pairing.stored_heat_j).abs() / pairing.stored_heat_j < 1.0e-4
    );
    assert!(migrated.thermal_residual_j.abs() < 1.0e-12);
    assert!(
        (migrated.hot_tension_n - pairing.hot_tension_n).abs() / pairing.hot_tension_n < 1.0e-4
    );
    let stresses = migrated.hot_layer_stresses_pa.as_ref().unwrap();
    for (m, p) in stresses.iter().zip(&pairing.layer_stresses_pa) {
        assert!((m - p).abs() / p.abs().max(1.0) < 1.0e-4);
    }
    let peak = pairing
        .hot_pa
        .iter()
        .map(|p| p.abs())
        .fold(0.0_f64, f64::max);
    let max_diff = pairing
        .hot_pa
        .iter()
        .zip(&migrated.hot_pa)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(peak > 0.0 && max_diff / peak < 1.0e-4);

    // Inspect the exact window ledger directly on the migrated multirate rig.
    let mut rig = example::build_multirate_section(
        [&core, &sleeve],
        &charts,
        &reference,
        migrated.cold_tension_n,
        0.2,
        0.0,
        0.0,
    )
    .unwrap();
    example::step_multirate_section(&mut rig, 0.1, example::ACOUSTIC_SUBSTEPS, 16).unwrap();
    let ledger = rig.sys.ledger();
    assert!(ledger.energy_residual_j().abs() < 1.0e-12);
    assert!(ledger.connection_defect_j().abs() < 1.0e-10);
    assert!(ledger.numerical_loss_j().abs() < 1.0e-12);
    assert!(ledger.entropy_generated_j_per_k() > 0.0);
    match ledger.outcome().energy {
        Color::Verified { .. } | Color::Validated { .. } => {}
        Color::Estimated { dispersion, .. } => assert!(dispersion < 1.0e-10),
    }
    assert!(matches!(ledger.outcome().entropy, Color::Verified { .. }));
}

#[test]
fn multirate_100x_damped_string_conserves_modal_heat_into_thermal_layers() {
    let core = example::demo_with_expansion(1.0, None).unwrap();
    let sleeve = example::demo_with_expansion(2.0, None).unwrap();
    let charts = [
        example::caloric(&core).unwrap(),
        example::caloric(&sleeve).unwrap(),
    ];
    let reference = example::section([&core, &sleeve], false).unwrap();

    // 100x multirate run: 4 macro steps of 6.25 ms, 100 fast 16 kHz substeps each.
    let mut mr = example::build_multirate_section(
        [&core, &sleeve],
        &charts,
        &reference,
        10.0,
        0.2,
        2.0e-2,
        0.4,
    )
    .unwrap();
    let init_modal_e = mr.modal_owner.stored_energy_j(&mr.init_modes).unwrap();
    example::step_multirate_section(&mut mr, 0.025, example::ACOUSTIC_SUBSTEPS, 4).unwrap();

    // Single-rate 1x reference over the same 25 ms (400 steps of 62.5 us).
    let mut sr = example::build_multirate_section(
        [&core, &sleeve],
        &charts,
        &reference,
        10.0,
        0.2,
        2.0e-2,
        0.4,
    )
    .unwrap();
    example::step_multirate_section(&mut sr, 0.025, 1, 400).unwrap();

    let ledger = mr.sys.ledger();
    let elec_j = ledger.nodes[mr.core_r.index()].dissipated_j.value()
        + ledger.nodes[mr.sleeve_r.index()].dissipated_j.value();
    let modal_diss_j = ledger.nodes[mr.string_modes.index()].dissipated_j.value();
    let stored_heat_j = ledger.nodes[mr.core_body.index()].stored_change_j.value()
        + ledger.nodes[mr.sleeve_body.index()].stored_change_j.value();
    let final_modal_state = mr
        .sys
        .state::<Vec<ModalAcousticState>>(mr.string_modes)
        .unwrap();
    let final_modal_e = mr.modal_owner.stored_energy_j(final_modal_state).unwrap();

    // Nonzero fast modal damping heat crosses the 100x rate boundary into the slow
    // enthalpy bodies and closes the whole-system energy account to roundoff.
    assert!(modal_diss_j > 1.0e-3);
    assert!((modal_diss_j - (init_modal_e - final_modal_e)).abs() < 1.0e-12);
    assert!((stored_heat_j - (elec_j + modal_diss_j)).abs() < 1.0e-12);
    assert!(ledger.energy_residual_j().abs() < 1.0e-12);
    assert!(ledger.connection_defect_j().abs() < 1.0e-12);
    assert!(ledger.numerical_loss_j().abs() < 1.0e-12);
    assert!(matches!(ledger.outcome().entropy, Color::Verified { .. }));

    // 100x multirate agrees with the 1x single-rate reference.
    let sr_ledger = sr.sys.ledger();
    let sr_stored_j = sr_ledger.nodes[sr.core_body.index()]
        .stored_change_j
        .value()
        + sr_ledger.nodes[sr.sleeve_body.index()]
            .stored_change_j
            .value();
    let sr_modal_diss = sr_ledger.nodes[sr.string_modes.index()]
        .dissipated_j
        .value();
    assert!((modal_diss_j - sr_modal_diss).abs() / sr_modal_diss < 1.0e-10);
    assert!((stored_heat_j - sr_stored_j).abs() / sr_stored_j < 1.0e-5);
    for (mr_node, sr_node) in [
        (mr.core_body, sr.core_body),
        (mr.sleeve_body, sr.sleeve_body),
    ] {
        let t_mr = mr
            .sys
            .state::<EquilibriumPhaseState>(mr_node)
            .unwrap()
            .temperature_k();
        let t_sr = sr
            .sys
            .state::<EquilibriumPhaseState>(sr_node)
            .unwrap()
            .temperature_k();
        assert!(
            (t_mr - t_sr).abs() < 2.0e-4,
            "100x vs 1x T mismatch: {t_mr} vs {t_sr}"
        );
    }
}
