//! Original HEM flow and enthalpy conversion, with independent PR references.
#![cfg(feature = "material")]

#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;
use common::{MIXTURE_MASS, model, source};
use fs_gas::material::phase_nozzle::{
    MixtureNozzleControl as Control, MixtureNozzleError as Error, MixtureNozzleFlow as Flow,
    MixtureNozzleStage, equilibrium_mixture_nozzle,
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::ThermochemicalFlashError;
use fs_qty::{Area, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}",
        actual - expected
    );
}
fn controls() -> Control {
    Control::new(
        [Temperature::new(180.0), Temperature::new(240.0)],
        Pressure::new(1.5e6),
    )
}
fn nozzle(
    original: &ThermochemicalPengRobinsonModel,
    back: f64,
    area: f64,
    control: Control,
    cancelled: impl FnMut() -> bool,
) -> Result<Flow, Error> {
    equilibrium_mixture_nozzle(
        original,
        Temperature::new(230.0),
        Pressure::new(3e6),
        Pressure::new(back),
        Area::new(area),
        control,
        cancelled,
    )
}

#[test]
fn g0_prepared_entropy_paths_retain_the_original_inlet_across_pressure_queries() {
    let original = source();
    let control = controls().isentropic;
    let path = original
        .prepare_phase_entropy_path(Temperature::new(230.0), Pressure::new(3e6), control, || {
            false
        })
        .unwrap();
    assert_eq!(path.preparation_pt_evaluations(), 1);
    for pressure in [3e6, 2.5e6, 1.8e6, 2.5e6] {
        let a = path
            .solve_increment(Pressure::new(pressure), 0.0, 64, || false)
            .unwrap();
        let b = original
            .phase_equilibrium_pressure_entropy_increment(
                Temperature::new(230.0),
                Pressure::new(3e6),
                Pressure::new(pressure),
                0.0,
                control,
                || false,
            )
            .unwrap();
        assert_eq!(a.state(), b.state());
        assert_eq!(a.enthalpy_change_j_kg(), b.enthalpy_change_j_kg());
        assert_eq!(a.accepted_increment_j_kg_k(), b.accepted_increment_j_kg_k());
        assert_eq!(a.pt_evaluations() + 1, b.pt_evaluations());
        assert_eq!(
            a.eos_evaluation_upper_bound() + path.preparation_eos_evaluation_upper_bound(),
            b.eos_evaluation_upper_bound()
        );
    }
    let rest = path
        .solve_increment(Pressure::new(3e6), 0.0, 0, || false)
        .unwrap();
    assert_eq!(rest.pt_evaluations(), 0);
    assert_eq!(rest.enthalpy_change_j_kg(), Some(0.0));
    assert!(matches!(
        path.solve_increment(Pressure::new(2.5e6), 0.0, 0, || false),
        Err(ThermochemicalFlashError::WorkLimit {
            flash_evaluations: 0
        })
    ));
    assert!(matches!(
        path.solve_increment(Pressure::new(2.5e6), 0.0, 64, || true),
        Err(ThermochemicalFlashError::Cancelled)
    ));
    assert_eq!(
        path.solve_increment(Pressure::new(3e6), 0.0, 0, || false)
            .unwrap(),
        rest
    );
}

#[test]
fn g1_two_phase_sonic_flow_matches_independent_pr_isentrope_and_flux_maximum() {
    // Independent oracle solves PR fugacity and entropy simultaneously for T,x,y,
    // differentiates rho(p)|s at refined steps, and separately maximizes rho*v.
    // It uses original fixture inputs, not the production PS/acoustic solver.
    let original = source();
    let flow = nozzle(&original, 1e6, 1e-4, controls(), || false).unwrap();
    assert!(flow.choked);
    close(
        flow.outlet.state().pressure_pa(),
        1_859_777.623_720_862_2,
        2.0,
    );
    close(
        flow.outlet.state().temperature_k(),
        215.790_125_620_931_7,
        5e-5,
    );
    close(
        flow.outlet.state().density_kg_m3(),
        52.822_407_977_044_385,
        1e-4,
    );
    close(flow.velocity_m_s, 181.367_503_109_847_04, 3e-4);
    close(flow.mass_flux_kg_m2_s, 9580.268_243_046_206, 0.01);
    close(flow.mass_flow.value(), 0.958_026_824_304_620_6, 1e-6);
    close(
        flow.inlet.vapor_mass_fraction().unwrap(),
        0.408_352_234_633_990_06,
        2e-8,
    );
    close(
        flow.outlet.state().vapor_mass_fraction().unwrap(),
        0.434_611_145_162_811,
        2e-7,
    );
    close(
        flow.outlet.enthalpy_change_j_kg().unwrap(),
        -16447.085_592_150_186,
        0.06,
    );
    assert!(
        flow.outlet.state().vapor_mass_fraction().unwrap()
            > flow.inlet.vapor_mass_fraction().unwrap()
    );
    assert!(flow.outlet.residual_j_kg_k().abs() <= controls().isentropic.entropy_tolerance_j_kg_k);
    assert!(flow.stagnation_enthalpy_residual_j_kg.abs() < 1e-9);
    close(flow.mach_number, 1.0, 2e-8);
    let acoustic = flow.acoustics.as_ref().unwrap();
    assert_eq!(acoustic.source(), flow.outlet.state());
    let bracket = flow.sonic_pressure_bracket_pa.unwrap();
    assert!(
        bracket[0] <= flow.outlet.state().pressure_pa()
            && flow.outlet.state().pressure_pa() <= bracket[1]
    );
    assert!(
        bracket[1] - bracket[0]
            <= controls().absolute_pressure_tolerance_pa
                + controls().relative_pressure_tolerance * 3e6
    );
}

#[test]
fn g1_unchoked_exit_uses_native_outlet_density_and_equilibrium_sound_speed() {
    let original = source();
    let flow = nozzle(&original, 2.5e6, 1e-4, controls(), || false).unwrap();
    assert!(!flow.choked);
    assert_eq!(flow.outlet.state().pressure_pa(), 2.5e6);
    assert_eq!(flow.sonic_pressure_bracket_pa, None);
    assert_eq!(flow.pressure_evaluations, 2);
    close(
        flow.outlet.state().temperature_k(),
        224.379_060_551_629_8,
        2e-6,
    );
    close(flow.velocity_m_s, 111.027_693_113_095_44, 2e-5);
    close(
        flow.acoustics.as_ref().unwrap().state().sound_speed_m_s(),
        174.569_468_996_362_99,
        2e-3,
    );
    close(flow.mass_flow.value(), 0.810_932_674_987_624_3, 2e-7);
    assert!(flow.mach_number < 1.0);
    assert!(flow.sonic_residual_m2_s2.unwrap() > 0.0);
    assert!(flow.stagnation_enthalpy_residual_j_kg.abs() < 1e-9);
    // An upstream-density mass flux would be about 23% too high.
    assert!(
        (flow.mass_flow.value() - 1e-4 * flow.inlet.density_kg_m3() * flow.velocity_m_s).abs()
            > 0.18
    );
}

#[test]
fn g1_dilute_nozzle_recovers_analytic_choking_and_finite_pressure_pr_reference() {
    let original = source();
    let control = Control::new(
        [Temperature::new(280.0), Temperature::new(410.0)],
        Pressure::new(40.0),
    );
    let solve = |back| {
        equilibrium_mixture_nozzle(
            &original,
            Temperature::new(400.0),
            Pressure::new(100.0),
            Pressure::new(back),
            Area::new(1.0),
            control,
            || false,
        )
        .unwrap()
    };
    let flow = solve(20.0);
    assert!(flow.choked);
    let gamma = 4.0 / 3.0;
    let rs = R / MIXTURE_MASS;
    let t_star = 2.0 * 400.0 / (gamma + 1.0);
    let p_star = 100.0 * (2.0_f64 / (gamma + 1.0)).powf(gamma / (gamma - 1.0));
    let speed = (gamma * rs * t_star).sqrt();
    let flux = p_star * speed / (rs * t_star);
    // The actual finite-100-Pa PR model has small nonzero departures.
    close(flow.outlet.state().temperature_k(), t_star, 1.2e-4);
    close(flow.outlet.state().pressure_pa(), p_star, 5e-6);
    close(flow.velocity_m_s, speed, 5e-4);
    close(flow.mass_flux_kg_m2_s, flux, 2.5e-7);
    close(flow.mass_flow.value(), 0.177_255_034_989_548_48, 2e-8);
    let unchoked = solve(80.0);
    assert!(!unchoked.choked);
    close(
        unchoked.outlet.state().temperature_k(),
        378.296_603_041_566_75,
        3e-6,
    );
    close(unchoked.velocity_m_s, 250.227_764_845_557_4, 2e-5);
    close(unchoked.mass_flux_kg_m2_s, 0.146_736_752_202_682_4, 2e-8);
}

#[test]
fn g3_choked_flow_is_back_pressure_independent_and_area_scales_only_mass_rate() {
    let original = source();
    let a = nozzle(&original, 1e6, 1e-4, controls(), || false).unwrap();
    for (back, area) in [(0.0, 1e-4), (0.7e6, 2e-4)] {
        let b = nozzle(&original, back, area, controls(), || false).unwrap();
        assert_eq!(a.outlet, b.outlet);
        assert_eq!(a.acoustics, b.acoustics);
        assert_eq!(a.velocity_m_s, b.velocity_m_s);
        assert_eq!(a.mass_flux_kg_m2_s, b.mass_flux_kg_m2_s);
        assert_eq!(a.pt_evaluations, b.pt_evaluations);
        assert_eq!(b.mass_flow.value(), a.mass_flow.value() * (area / 1e-4));
        assert!(b.choked);
    }
}

#[test]
fn g3_species_enthalpy_and_entropy_datums_do_not_change_nozzle_flow() {
    let original = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let a = nozzle(&original, 1e6, 1e-4, controls(), || false).unwrap();
    let b = nozzle(&shifted, 1e6, 1e-4, controls(), || false).unwrap();
    assert_eq!(a.mass_flow, b.mass_flow);
    assert_eq!(a.velocity_m_s, b.velocity_m_s);
    assert_eq!(
        a.outlet.state().temperature_k(),
        b.outlet.state().temperature_k()
    );
    assert_eq!(
        a.outlet.state().pressure_pa(),
        b.outlet.state().pressure_pa()
    );
    assert_eq!(a.outlet.state().phase(), b.outlet.state().phase());
    assert_eq!(
        a.outlet.enthalpy_change_j_kg(),
        b.outlet.enthalpy_change_j_kg()
    );
    assert_eq!(a.outlet.residual_j_kg_k(), b.outlet.residual_j_kg_k());
    assert_eq!(a.pt_evaluations, b.pt_evaluations);
}

#[test]
fn g0_exact_rest_and_closed_area_need_one_preparation_and_no_acoustic_work() {
    let original = source();
    let mut control = controls();
    control.maximum_pt_evaluations = 1;
    control.maximum_pressure_evaluations = 1;
    for (back, area) in [(3e6, 1e-4), (0.0, 0.0)] {
        let flow = nozzle(&original, back, area, control, || false).unwrap();
        assert_eq!(flow.mass_flow.value(), 0.0);
        assert_eq!(flow.velocity_m_s, 0.0);
        assert_eq!(flow.outlet.enthalpy_change_j_kg(), Some(0.0));
        assert_eq!(flow.outlet.state(), &flow.inlet);
        assert_eq!(flow.outlet.pt_evaluations(), 0);
        assert_eq!(flow.pt_evaluations, 1);
        assert_eq!(flow.pressure_evaluations, 1);
        assert_eq!(flow.acoustics, None);
        assert!(!flow.choked);
    }
}

#[test]
fn g0_source_domains_phase_kinks_and_insufficient_pressure_support_refuse() {
    let original = source();
    let mut control = controls();
    control.minimum_pressure_pa = 2.5e6;
    assert!(matches!(nozzle(&original, 1e6, 1e-4, control, || false),
        Err(Error::OutsidePressureBracket { minimum_pressure_pa, .. }) if minimum_pressure_pa == 2.5e6));
    control = controls();
    control.isentropic.temperature_bracket = [Temperature::new(225.0), Temperature::new(240.0)];
    assert!(matches!(
        nozzle(&original, 1e6, 1e-4, control, || false),
        Err(Error::Equilibrium {
            stage: MixtureNozzleStage::Isentropic,
            ..
        })
    ));
    let crossing = Control::new(
        [Temperature::new(170.0), Temperature::new(280.0)],
        Pressure::new(1e6),
    );
    // Independent dew crossing is at 2.414422846 MPa. The first 2.25 MPa
    // scan state has two phases and remains subsonic (a²-v²≈32246 m²/s²).
    assert!(matches!(
        equilibrium_mixture_nozzle(&original, Temperature::new(260.0), Pressure::new(3e6),
            Pressure::new(1e6), Area::new(1e-4), crossing, || false),
        Err(Error::PhaseBoundary { pressure_pa, initial_phases: 1, trial_phases: 2 }) if pressure_pa == 2.25e6
    ));
    control = controls();
    control.acoustics.phase.flash.log_k_tolerance *= 2.0;
    assert!(matches!(
        nozzle(&original, 1e6, 1e-4, control, || false),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        nozzle(&original, 3.1e6, 1e-4, controls(), || false),
        Err(Error::Invalid(_))
    ));
}

#[test]
fn g0_shared_work_budget_includes_preparation_all_pressure_trials_and_acoustics() {
    let original = source();
    let accepted = nozzle(&original, 1e6, 1e-4, controls(), || false).unwrap();
    let mut control = controls();
    control.maximum_pt_evaluations = accepted.pt_evaluations;
    control.maximum_pressure_evaluations = accepted.pressure_evaluations;
    assert_eq!(
        nozzle(&original, 1e6, 1e-4, control, || false).unwrap(),
        accepted
    );
    control.maximum_pt_evaluations -= 1;
    assert!(matches!(
        nozzle(&original, 1e6, 1e-4, control, || false),
        Err(Error::WorkLimit { .. })
            | Err(Error::Equilibrium {
                source: ThermochemicalFlashError::WorkLimit { .. },
                ..
            })
    ));
    control.maximum_pt_evaluations = accepted.pt_evaluations;
    control.maximum_pressure_evaluations -= 1;
    assert!(matches!(nozzle(&original, 1e6, 1e-4, control, || false),
        Err(Error::WorkLimit { pressure_evaluations, .. })
        if pressure_evaluations == control.maximum_pressure_evaluations));
}

#[test]
fn g4_final_cancellation_publishes_no_flow_and_replay_is_exact() {
    let original = source();
    let mut polls = 0;
    let accepted = nozzle(&original, 2.5e6, 1e-4, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut count = 0;
    assert_eq!(
        nozzle(&original, 2.5e6, 1e-4, controls(), || {
            count += 1;
            count == polls
        })
        .unwrap_err(),
        Error::Cancelled
    );
    assert_eq!(count, polls);
    assert_eq!(
        nozzle(&original, 2.5e6, 1e-4, controls(), || false).unwrap(),
        accepted
    );
}
