//! Original-law derivatives and explicitly equilibrated acoustic consumers.
#[path = "common/nonideal_equilibrium.rs"]
mod common;
use common::{MOLAR_MASS, controls, model};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics, FrozenHomogeneousMixture};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::equilibrium::{
    PureFluidEquilibriumPhase as Phase, PureFluidEquilibriumState,
};
use fs_qty::Temperature;

fn sample(
    source: &ThermochemicalPengRobinsonModel,
    temperature: f64,
    density: f64,
) -> PureFluidEquilibriumState {
    source
        .evaluate_pure_equilibrium_density_temperature(
            density,
            Temperature::new(temperature),
            controls(),
            || false,
        )
        .unwrap()
}
fn at_quality(source: &ThermochemicalPengRobinsonModel, quality: f64) -> PureFluidEquilibriumState {
    let saturation = source
        .pure_saturation_at_temperature(
            Temperature::new(150.0),
            controls().pressure_bracket,
            controls().saturation_options,
        )
        .unwrap();
    let phase = saturation.saturation().at_vapor_fraction(quality).unwrap();
    sample(source, 150.0, phase.mass_density_kg_m3())
}
fn close(value: f64, expected: f64, tolerance: f64) {
    assert!(
        (value - expected).abs() <= tolerance * expected.abs(),
        "{value} != {expected}"
    );
}
fn entropy_neighbor(
    source: &ThermochemicalPengRobinsonModel,
    density: f64,
    entropy: f64,
) -> PureFluidEquilibriumState {
    // Independent entropy constraint: no acoustic derivative or response Cv
    // enters this bounded bisection through full original equilibrium states.
    let mut lo = 130.0;
    let mut hi = 180.0;
    for _ in 0..55 {
        let state = sample(source, lo + 0.5 * (hi - lo), density);
        if state.specific_entropy_j_kg_k() < entropy {
            lo = state.temperature_k();
        } else {
            hi = state.temperature_k();
        }
    }
    let state = sample(source, lo + 0.5 * (hi - lo), density);
    assert!((state.specific_entropy_j_kg_k() - entropy).abs() < 1e-7);
    state
}

#[test]
fn g1_equilibrium_sound_matches_independent_entropy_preserving_compression() {
    let source = model(0.0, [100.0, 250.0]);
    // Independent NumPy PR roots + Maxwell-area brentq, then an isentropic
    // temperature solve using full analytical entropy, not the new Jacobian.
    for (quality, oracle_c) in [
        (0.01, 29.335_537_597_650_53),
        (0.5, 178.869_227_661_574_48),
        (0.9, 252.240_558_185_836_8),
    ] {
        let center = at_quality(&source, quality);
        let acoustic = center.equilibrium_acoustics().unwrap();
        close(acoustic.state().sound_speed_m_s(), oracle_c, 2e-6);
        let rho = center.density_kg_m3();
        let delta = rho * 1e-4;
        let below = entropy_neighbor(&source, rho - delta, center.specific_entropy_j_kg_k());
        let above = entropy_neighbor(&source, rho + delta, center.specific_entropy_j_kg_k());
        let numerical_c2 = (above.pressure_pa() - below.pressure_pa()) / (2.0 * delta);
        close(
            acoustic.derivatives().pressure_density_at_entropy_m2_s2,
            numerical_c2,
            3e-6,
        );
        close(
            acoustic
                .derivatives()
                .temperature_density_at_entropy_k_m3_kg,
            (above.temperature_k() - below.temperature_k()) / (2.0 * delta),
            3e-6,
        );
    }
}

#[test]
fn g1_pressure_entropy_and_conserved_energy_jacobians_match_original_law_differences() {
    let source = model(0.0, [100.0, 250.0]);
    let center = at_quality(&source, 0.5);
    let response = center.equilibrium_acoustics().unwrap();
    let jacobian = response.derivatives();
    let rho = center.density_kg_m3();
    let dt = 0.001;
    let lower_t = sample(&source, 150.0 - dt, rho);
    let upper_t = sample(&source, 150.0 + dt, rho);
    close(
        jacobian.pressure_temperature_at_density_pa_k,
        (upper_t.pressure_pa() - lower_t.pressure_pa()) / (2.0 * dt),
        3e-7,
    );
    close(
        jacobian.entropy_temperature_at_density_j_kg_k2,
        (upper_t.specific_entropy_j_kg_k() - lower_t.specific_entropy_j_kg_k()) / (2.0 * dt),
        3e-7,
    );
    close(
        jacobian.isochoric_heat_capacity_j_kg_k,
        (upper_t.specific_internal_energy_j_kg() - lower_t.specific_internal_energy_j_kg())
            / (2.0 * dt),
        3e-7,
    );
    let dr = rho * 1e-4;
    let lower_rho = sample(&source, 150.0, rho - dr);
    let upper_rho = sample(&source, 150.0, rho + dr);
    assert_eq!(jacobian.pressure_density_at_temperature_m2_s2, 0.0);
    assert_eq!(lower_rho.pressure_pa(), upper_rho.pressure_pa());
    close(
        jacobian.entropy_density_at_temperature_j_m3_kg2_k,
        (upper_rho.specific_entropy_j_kg_k() - lower_rho.specific_entropy_j_kg_k()) / (2.0 * dr),
        3e-7,
    );
    let energy = center.specific_internal_energy_j_kg();
    let low_u = source
        .reconstruct_pure_density_energy(rho, energy - 10.0, controls(), || false)
        .unwrap();
    let high_u = source
        .reconstruct_pure_density_energy(rho, energy + 10.0, controls(), || false)
        .unwrap();
    close(
        jacobian.pressure_energy_at_density_kg_m3,
        (high_u.state().pressure_pa() - low_u.state().pressure_pa()) / 20.0,
        3e-6,
    );
    let low_rho_u = source
        .reconstruct_pure_density_energy(rho - dr, energy, controls(), || false)
        .unwrap();
    let high_rho_u = source
        .reconstruct_pure_density_energy(rho + dr, energy, controls(), || false)
        .unwrap();
    close(
        jacobian.pressure_density_at_energy_m2_s2,
        (high_rho_u.state().pressure_pa() - low_rho_u.state().pressure_pa()) / (2.0 * dr),
        3e-6,
    );
    close(
        jacobian.pressure_density_at_entropy_m2_s2,
        jacobian.pressure_density_at_energy_m2_s2
            + (center.pressure_pa() / rho) * (jacobian.pressure_energy_at_density_kg_m3 / rho),
        3e-15,
    );
}

#[test]
fn g3_equilibrated_response_is_distinct_from_frozen_wood_and_preserves_source_datums() {
    let source = model(0.0, [100.0, 250.0]);
    for quality in [0.01, 0.1, 0.5, 0.9] {
        let center = at_quality(&source, quality);
        let response = center.equilibrium_acoustics().unwrap();
        assert!(core::ptr::eq(response.source(), &center));
        let Phase::Coexistence(mixture) = center.phase() else {
            panic!("expected coexistence");
        };
        let ideal_cv = center.reference().properties().cv().value();
        let phase = |root: &fs_thermochem::eos::peng_robinson::PengRobinsonPhaseRootV1| {
            let mechanics = root.caloric_mechanics(ideal_cv, MOLAR_MASS).unwrap();
            FluidAcousticState::try_new(root.mass_density_kg_per_m3(), mechanics.sound_speed_m_s)
                .unwrap()
        };
        let frozen = FrozenHomogeneousMixture::from_mass_fraction(
            phase(mixture.saturation().liquid()),
            phase(mixture.saturation().vapor()),
            quality,
        )
        .unwrap();
        assert!(response.state().sound_speed_m_s() < frozen.state().sound_speed_m_s());
        let consumer: &dyn FluidAcoustics = &response;
        assert_eq!(consumer.acoustic_state().unwrap(), response.state());
        assert_eq!(
            consumer.acoustic_state().unwrap().density_kg_m3(),
            center.density_kg_m3()
        );
        let shifted = at_quality(&model(1e20, [100.0, 250.0]), quality);
        let shifted_response = shifted.equilibrium_acoustics().unwrap();
        assert_eq!(shifted_response.state(), response.state());
        assert_eq!(shifted_response.derivatives(), response.derivatives());
    }
}

#[test]
fn g0_original_homogeneous_mechanics_and_nondifferentiable_phase_endpoints_remain_explicit() {
    let source = model(0.0, [100.0, 250.0]);
    for volume in [0.00005, 0.0005] {
        let state = sample(&source, 170.0, MOLAR_MASS / volume);
        let root = match state.phase() {
            Phase::Liquid(root) | Phase::Vapor(root) => root,
            Phase::Coexistence(_) => panic!("expected original homogeneous branch"),
        };
        let original = root
            .caloric_mechanics(state.reference().properties().cv().value(), MOLAR_MASS)
            .unwrap();
        let response = state.equilibrium_acoustics().unwrap();
        assert_eq!(response.state().sound_speed_m_s(), original.sound_speed_m_s);
        assert!(response.derivatives().pressure_density_at_temperature_m2_s2 > 0.0);
        close(
            response.derivatives().pressure_density_at_entropy_m2_s2,
            response.derivatives().pressure_density_at_energy_m2_s2
                + state.pressure_pa() / state.density_kg_m3()
                    * (response.derivatives().pressure_energy_at_density_kg_m3
                        / state.density_kg_m3()),
            3e-15,
        );
    }
    let saturation = source
        .pure_saturation_at_temperature(
            Temperature::new(150.0),
            controls().pressure_bracket,
            controls().saturation_options,
        )
        .unwrap();
    for endpoint in [0.0, 1.0] {
        assert!(
            saturation
                .saturation()
                .at_vapor_fraction(endpoint)
                .unwrap()
                .equilibrium_response(saturation.reference().properties().cv().value())
                .is_err()
        );
    }
}
