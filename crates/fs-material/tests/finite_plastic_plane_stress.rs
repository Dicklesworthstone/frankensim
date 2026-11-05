//! Independent sheet-forming references and the derivative of transverse
//! equilibrium. Authored constitutive fixtures are not material calibration.

use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::plane_stress::FinitePlaneStressResponse;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalError,
    FiniteJ2ThermalOverstress, FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2State};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const H: f64 = 8e6;
const LIMITS: FiniteJ2Limits = FiniteJ2Limits {
    minimum_jacobian: 0.05,
    maximum_elastic_log_strain: 0.8,
    maximum_plastic_log_strain: 3.0,
    maximum_accumulated_plastic_strain: 3.0,
};

fn parameters() -> FiniteJ2Parameters {
    FiniteJ2Parameters {
        shear_modulus_pa: MU,
        bulk_modulus_pa: K,
        yield_stress_pa: Y,
        hardening_modulus_pa: H,
    }
}

fn law() -> FiniteJ2 {
    FiniteJ2::new(parameters(), LIMITS).unwrap()
}

fn settings() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 1e-13,
        ..PlaneStressSettings::default()
    }
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}, difference={:.8e}",
        (actual - expected).abs()
    );
}

fn traction_free(response: &FinitePlaneStressResponse) {
    assert!(response.thickness_stretch() > 0.0);
    assert!(response.relative_transverse_residual <= settings().relative_tolerance);
    for component in [2, 5, 8] {
        near(response.material.first_piola_pa[component], 0.0, 3e-4, 0.0);
    }
}

fn metric_determinant(state: &FiniteJ2State) -> f64 {
    let [xx, yy, zz, xy, yz, zx] = state.inverse_plastic_cauchy_green;
    xx * yy * zz + 2.0 * xy * yz * zx - xx * yz * yz - yy * zx * zx - zz * xy * xy
}

// Under equibiaxial loading tau = diag(s*q,s*q,0). Thus the axial
// elastic Hencky strain is s*q*C, while plastic axial strain is s*gamma/2.
// This scalar reference uses no spectral decomposition or local Newton solve.
fn equibiaxial_reference(stretch: f64, yield_stress: f64, radius_slope: f64) -> (f64, f64, f64) {
    let axial = stretch.ln();
    let compliance = 1.0 / (6.0 * MU) + 2.0 / (9.0 * K);
    let gamma = (axial.abs() - yield_stress * compliance) / (0.5 + radius_slope * compliance);
    assert!(gamma > 0.0);
    let q = yield_stress + radius_slope * gamma;
    let thickness = (-2.0 * axial + axial.signum() * 2.0 * q / (3.0 * K)).exp();
    (gamma, q, thickness)
}

#[test]
fn independent_biaxial_tension_and_compression_retain_isochoric_flow_and_elastic_volume() {
    for hardening in [0.0, H] {
        let material = FiniteJ2::new(
            FiniteJ2Parameters {
                hardening_modulus_pa: hardening,
                ..parameters()
            },
            LIMITS,
        )
        .unwrap();
        for stretch in [0.78_f64, 1.23] {
            let old = material.initial_state();
            let response = material
                .plane_stress_response(&[stretch, 0.0, 0.0, stretch], &old, [0.0, 0.0, 1.0], settings())
                .unwrap();
            let (gamma, q, thickness) = equibiaxial_reference(stretch, Y, hardening);
            let sign = stretch.ln().signum();
            traction_free(&response);
            near(response.thickness_stretch(), thickness, 2e-12, 2e-12);
            near(response.material.plastic_increment, gamma, 2e-12, 2e-12);
            near(response.material.state.alpha, gamma, 2e-12, 2e-12);
            for component in [0, 3] {
                near(response.piola_pa[component], sign * q / stretch, 3e-4, 2e-11);
            }
            near(response.piola_pa[1], 0.0, 1e-8, 0.0);
            near(response.piola_pa[2], 0.0, 1e-8, 0.0);
            let expected_metric = [(-sign * gamma).exp(), (-sign * gamma).exp(), (2.0 * sign * gamma).exp(), 0.0, 0.0, 0.0];
            for (actual, expected) in response.material.state.inverse_plastic_cauchy_green.iter().zip(expected_metric) {
                near(*actual, expected, 2e-11, 2e-11);
            }
            near(metric_determinant(&response.material.state), 1.0, 3e-12, 0.0);
            let expected_j = (sign * 2.0 * q / (3.0 * K)).exp();
            near(response.material.jacobian, expected_j, 3e-12, 2e-12);
            assert!((response.material.jacobian - 1.0).abs() > 0.005);
            assert_eq!(response.thickness_stretch() > 1.0, stretch < 1.0);
            let stored = q * q / (6.0 * MU) + 2.0 * q * q / (9.0 * K) + 0.5 * hardening * gamma * gamma;
            near(response.material.stored_energy_j_m3, stored, 1e-4, 3e-11);
            near(response.material.plastic_dissipation_j_m3, Y * gamma, 1e-5, 2e-11);
            near(response.material.incremental_potential_j_m3, stored + Y * gamma, 1e-4, 3e-11);
            assert_eq!(old, material.initial_state());
        }
    }
}

#[test]
fn unloaded_sheet_keeps_permanent_thickness_and_hardening_without_additional_flow() {
    let material = law();
    for stretch in [0.82_f64, 1.19] {
        let loaded = material
            .plane_stress_response(&[stretch, 0.0, 0.0, stretch], &material.initial_state(), [0.0, 0.0, 1.0], settings())
            .unwrap();
        let old = loaded.material.state.clone();
        let sign = stretch.ln().signum();
        let residual_stretch = (0.5 * sign * old.alpha).exp();
        let unloaded = material
            .plane_stress_response(&[residual_stretch, 0.0, 0.0, residual_stretch], &old, loaded.transverse_column(), settings())
            .unwrap();
        traction_free(&unloaded);
        for stress in unloaded.material.first_piola_pa {
            near(stress, 0.0, 3e-4, 0.0);
        }
        near(unloaded.thickness_stretch(), (-sign * old.alpha).exp(), 2e-12, 2e-12);
        near(unloaded.material.jacobian, 1.0, 3e-12, 0.0);
        assert_eq!(unloaded.material.state, old);
        assert_eq!(unloaded.material.plastic_increment, 0.0);
        assert_eq!(unloaded.material.plastic_dissipation_j_m3, 0.0);
        near(unloaded.material.stored_energy_j_m3, 0.5 * H * old.alpha * old.alpha, 1e-7, 1e-12);
        assert!((unloaded.thickness_stretch() - 1.0).abs() > 0.1);
        assert_eq!(old, loaded.material.state);
    }
}

#[test]
fn retained_transverse_shear_is_relaxed_in_both_directions() {
    let material = law();
    let (xz, yz) = (0.23, -0.17);
    // Cp^-1 = Fp^-1 Fp^-T for the isochoric triangular residual gradient
    // Fp = [1,0,xz; 0,1,yz; 0,0,1]. The positive accumulation admits this
    // metric within the native reachability bound.
    let old = FiniteJ2State {
        inverse_plastic_cauchy_green: [1.0 + xz * xz, 1.0 + yz * yz, 1.0, xz * yz, -yz, -xz],
        alpha: 0.5,
        creep_alpha: 0.0,
    };
    material.validate_history(&old).unwrap();
    let preserved = old.clone();
    let response = material
        .plane_stress_response(&[1.0, 0.0, 0.0, 1.0], &old, [0.0, 0.0, 1.0], settings())
        .unwrap();
    traction_free(&response);
    for (actual, expected) in response.transverse_column().into_iter().zip([xz, yz, 1.0]) {
        near(actual, expected, 3e-12, 3e-12);
    }
    for stress in response.material.first_piola_pa {
        near(stress, 0.0, 3e-4, 0.0);
    }
    assert_eq!(response.material.state, preserved);
    assert_eq!(response.material.plastic_increment, 0.0);
    assert_eq!(old, preserved);
}

#[test]
fn condensed_tangent_and_potential_gradient_use_frozen_history_on_every_perturbation() {
    let material = law();
    let initial = material.initial_state();
    let retained = material
        .response(&[1.0, 0.12, 0.07, 0.0, 0.98, -0.04, 0.0, 0.0, 1.0 / 0.98], &initial)
        .unwrap()
        .state;
    for (in_plane, old) in [([1.001, 0.002, -0.001, 0.999], &initial), ([1.08, 0.16, -0.035, 0.97], &retained)] {
        let preserved = old.clone();
        let response = material.plane_stress_response(&in_plane, old, [0.0, 0.0, 1.0], settings()).unwrap();
        traction_free(&response);
        let mut squared_error = 0.0;
        let mut squared_reference = 0.0;
        for column in 0..4 {
            let h = 2e-6;
            let mut plus = in_plane;
            let mut minus = in_plane;
            plus[column] += h;
            minus[column] -= h;
            let positive = material.plane_stress_response(&plus, old, response.transverse_column(), settings()).unwrap();
            let negative = material.plane_stress_response(&minus, old, response.transverse_column(), settings()).unwrap();
            near((positive.material.incremental_potential_j_m3 - negative.material.incremental_potential_j_m3) / (2.0 * h), response.piola_pa[column], 0.2, 3e-6);
            for row in 0..4 {
                let numerical = (positive.piola_pa[row] - negative.piola_pa[row]) / (2.0 * h);
                let analytic = response.tangent_pa[row][column];
                near(analytic, numerical, 5.0, 4e-6);
                near(analytic, response.tangent_pa[column][row], 1e-5, 3e-12);
                squared_error += (analytic - numerical).powi(2);
                squared_reference += numerical * numerical;
            }
        }
        assert!((squared_error / squared_reference).sqrt() < 6e-7);
        assert_eq!(old, &preserved);
    }
}

#[test]
fn rate_creep_and_combined_sheets_keep_one_physical_interval_and_native_accounts() {
    let eta = 60e6;
    let stretch: f64 = 1.14;
    let in_plane = [stretch, 0.0, 0.0, stretch];
    let mut previous_rate_increment = 0.0;
    for duration in [0.03, 0.3] {
        let kinetics = FiniteOverstressParameters { viscosity_pa_s: eta, stress_scale_pa: Y, exponent: 1.0 };
        let rate = FiniteOverstressJ2::new(law(), kinetics, duration).unwrap();
        let old = rate.initial_state();
        let response = rate.plane_stress_response(&in_plane, &old, [0.0, 0.0, 1.0], settings()).unwrap();
        let (gamma, q, thickness) = equibiaxial_reference(stretch, Y, H + eta / duration);
        traction_free(&response);
        near(response.material.plastic_increment, gamma, 3e-12, 3e-11);
        near(response.thickness_stretch(), thickness, 3e-12, 3e-12);
        near(response.material.equivalent_kirchhoff_pa, q, 1e-3, 3e-11);
        near(response.material.plastic_dissipation_j_m3, (Y + eta / duration * gamma) * gamma, 1e-4, 3e-11);
        assert!(gamma > previous_rate_increment);
        previous_rate_increment = gamma;
        assert_eq!(response.material, rate.response(&response.deformation_gradient, &old).unwrap());
        assert_eq!(old, rate.initial_state());

        let creep = FiniteNortonCreep::new(MU, K, LIMITS, kinetics, duration).unwrap();
        let old = creep.initial_state();
        let response = creep.plane_stress_response(&in_plane, &old, [0.0, 0.0, 1.0], settings()).unwrap();
        let (gamma, q, thickness) = equibiaxial_reference(stretch, 0.0, eta / duration);
        traction_free(&response);
        near(response.material.plastic_increment, gamma, 3e-12, 3e-11);
        near(response.thickness_stretch(), thickness, 3e-12, 3e-12);
        near(response.material.plastic_dissipation_j_m3, q * gamma, 1e-4, 3e-11);
        near(response.material.incremental_potential_j_m3, response.material.stored_energy_j_m3 + 0.5 * q * gamma, 1e-4, 3e-11);
        assert_eq!(response.material, creep.response(&response.deformation_gradient, &old).unwrap());
        assert_eq!(old, creep.initial_state());

        let combined = FiniteNortonJ2::new(law(), kinetics, duration, 3.0).unwrap();
        let old = combined.initial_state();
        let response = combined.plane_stress_response(&in_plane, &old, [0.0, 0.0, 1.0], settings()).unwrap();
        let compliance = 1.0 / (6.0 * MU) + 2.0 / (9.0 * K);
        let q = (stretch.ln() + Y / (2.0 * H)) / (compliance + 1.0 / (2.0 * H) + duration / (2.0 * eta));
        let plastic = (q - Y) / H;
        let creep = q * duration / eta;
        assert!(plastic > 0.0 && creep > 0.0);
        traction_free(&response);
        near(response.material.plastic_increment, plastic, 3e-12, 3e-11);
        near(response.material.creep_increment, creep, 3e-12, 3e-11);
        near(response.material.state.alpha, plastic, 3e-12, 3e-11);
        near(response.material.state.creep_alpha, creep, 3e-12, 3e-11);
        near(response.material.plastic_dissipation_j_m3, Y * plastic, 1e-4, 3e-11);
        near(response.material.creep_dissipation_j_m3, q * creep, 1e-4, 3e-11);
        near(response.material.incremental_potential_j_m3, response.material.stored_energy_j_m3 + Y * plastic + 0.5 * q * creep, 1e-4, 3e-11);
        assert_eq!(response.material, combined.response(&response.deformation_gradient, &old).unwrap());
        assert_eq!(old, combined.initial_state());
    }
}

#[test]
fn hot_sheet_uses_total_thermal_volume_softened_properties_and_selected_rate_law() {
    let parameters = FiniteJ2ThermalParameters {
        reference_temperature_k: 300.0,
        temperature_range_k: [280.0, 500.0],
        reference: parameters(),
        slopes: FiniteJ2TemperatureSlopes { shear_modulus_pa_k: -20_000.0, bulk_modulus_pa_k: -30_000.0, yield_stress_pa_k: -4_000.0, hardening_modulus_pa_k: -10_000.0 },
        logarithmic_expansion_per_k: 2e-4,
        logarithmic_expansion_slope_per_k2: 1e-7,
    };
    let kinetics = FiniteJ2ThermalOverstress { reference_viscosity_pa_s: 60e6, stress_scale_pa: Y, exponent: 1.0, activation_temperature_k: 1800.0 };
    let family = AffineThermalFiniteJ2::new(parameters, LIMITS).unwrap().with_overstress(kinetics).unwrap();
    let (temperature, duration) = (450.0, 0.08);
    let frozen = family.at_temperature_for_step(temperature, duration).unwrap();
    let theta = (2e-4_f64 * 150.0 + 0.5 * 1e-7 * 150.0_f64.powi(2)).exp();
    let old = family.initial_state();
    let free = frozen.plane_stress_response(&[theta, 0.0, 0.0, theta], &old, [0.0, 0.0, theta], settings()).unwrap();
    traction_free(&free);
    near(free.thickness_stretch(), theta, 2e-14, 2e-14);
    near(free.material.jacobian, theta.powi(3), 2e-14, 2e-14);
    assert_eq!(free.material.state, old);

    let native = FiniteOverstressJ2::new(
        FiniteJ2::new(FiniteJ2Parameters { shear_modulus_pa: MU - 20_000.0 * 150.0, bulk_modulus_pa: K - 30_000.0 * 150.0, yield_stress_pa: Y - 4_000.0 * 150.0, hardening_modulus_pa: H - 10_000.0 * 150.0 }, LIMITS).unwrap(),
        FiniteOverstressParameters { viscosity_pa_s: 60e6 * (1800.0_f64 * (1.0 / temperature - 1.0 / 300.0)).exp(), stress_scale_pa: Y, exponent: 1.0 },
        duration,
    ).unwrap();
    let mechanical = [1.1, 0.08, 0.0, 0.96];
    let expected = native.plane_stress_response(&mechanical, &old, [0.0, 0.0, 1.0], settings()).unwrap();
    let response = frozen.plane_stress_response(&mechanical.map(|value| theta * value), &old, [0.0, 0.0, theta], settings()).unwrap();
    traction_free(&response);
    for i in 0..9 {
        near(response.deformation_gradient[i], theta * expected.deformation_gradient[i], 3e-12, 3e-12);
    }
    for row in 0..4 {
        near(response.piola_pa[row], expected.piola_pa[row] / theta, 1e-3, 3e-11);
        for column in 0..4 {
            near(response.tangent_pa[row][column], expected.tangent_pa[row][column] / theta.powi(2), 1e-3, 3e-11);
        }
    }
    near(response.material.state.alpha, expected.material.state.alpha, 3e-12, 3e-11);
    near(response.material.plastic_dissipation_j_m3, expected.material.plastic_dissipation_j_m3, 1e-4, 3e-11);
    near(response.material.stored_energy_j_m3, expected.material.stored_energy_j_m3, 1e-4, 3e-11);
    assert_eq!(response.material, frozen.response(&response.deformation_gradient, &old).unwrap());
    assert_eq!(old, family.initial_state());
    let inspection = family.at_temperature_for_inspection(temperature).unwrap();
    assert!(matches!(inspection.plane_stress_response(&mechanical, &old, [0.0, 0.0, theta], settings()), Err(FiniteJ2ThermalError::DurationRequired)));
}

#[test]
fn invalid_geometry_exhausted_iteration_and_material_limits_refuse_without_changing_history() {
    let material = law();
    let old = material.initial_state();
    let preserved = old.clone();
    for (in_plane, seed) in [([1.0, 0.0, 0.0, -1.0], [0.0, 0.0, 1.0]), ([1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0]), ([f64::NAN, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0])] {
        assert!(material.plane_stress_response(&in_plane, &old, seed, settings()).is_err());
        assert_eq!(old, preserved);
    }
    assert!(material.plane_stress_response(&[1.23, 0.0, 0.0, 1.23], &old, [0.0, 0.0, 1.0], PlaneStressSettings { max_iterations: 1, ..settings() }).is_err());
    let limited = FiniteJ2::new(parameters(), FiniteJ2Limits { maximum_accumulated_plastic_strain: 0.02, ..LIMITS }).unwrap();
    assert!(limited.plane_stress_response(&[1.23, 0.0, 0.0, 1.23], &old, [0.0, 0.0, 1.0], settings()).is_err());
    assert_eq!(old, preserved);
}
