//! Independent Stokes geometry/response checks and exact passive source limits.
//! These authored fixtures qualify the mathematical reduction, not a helium
//! operating envelope or measured material-property interpolation.

use fs_material::superfluid::{
    MAX_NORMAL_VISCOUS_MODES, NormalViscousFactors, NormalViscousSpectrum, TwoFluidProperties,
    TwoFluidState,
};

fn state(viscosity: f64) -> TwoFluidState {
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 12.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: 0.0,
        normal_viscosity_pa_s: viscosity,
    })
    .unwrap()
}

fn spectrum(circle: bool, count: usize) -> NormalViscousSpectrum {
    if circle {
        NormalViscousSpectrum::circular_capillary(state(1.0), 1.0, count).unwrap()
    } else {
        NormalViscousSpectrum::parallel_plates(state(1.0), 1.0, count).unwrap()
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn every_prefix_retains_normal_inertia_and_exact_poiseuille_mobility() {
    for circle in [false, true] {
        let dc = if circle { 1.0 / 8.0 } else { 1.0 / 12.0 };
        for count in [1, 2, 4, 8, 16, 32] {
            let model = spectrum(circle, count);
            assert_eq!(model.resolved_mode_count(), count);
            assert_eq!(model.modes().len(), count + 1);
            assert!(model.modes().len() <= MAX_NORMAL_VISCOUS_MODES);
            assert_eq!(model.material(), &state(1.0));
            assert_eq!(model.modes().iter().map(|m| m.weight).sum::<f64>(), 1.0);
            assert!(
                model
                    .modes()
                    .iter()
                    .all(|m| m.weight > 0.0 && m.decay_rate_s_inv > 0.0)
            );
            let measured_dc = model
                .modes()
                .iter()
                .map(|m| m.weight / m.decay_rate_s_inv)
                .sum::<f64>();
            close(measured_dc, dc, 4e-16);
            close(model.dc_mobility_s().unwrap(), dc, 4e-16);
            close(model.wall().permeability_m2(), dc, 4e-16);
            let tail = model.truncation();
            assert!(tail.modeled_tail_decay_rate_s_inv >= tail.first_omitted_decay_rate_s_inv);
            close(
                tail.omitted_weight / tail.modeled_tail_decay_rate_s_inv,
                tail.omitted_dc_mobility_fraction * dc,
                1e-17,
            );
        }
    }
}

#[test]
fn circular_roots_satisfy_the_independent_j0_integral() {
    let model = spectrum(true, 32);
    // J0(x) = (1/pi) integral_0^pi cos(x*sin(theta)) dtheta.
    // Periodic smooth integrand; this quadrature is independent of the
    // high-precision power series used to produce the root constants.
    for root in model.modes()[..32]
        .iter()
        .map(|m| m.decay_rate_s_inv.sqrt())
        .chain([model.truncation().first_omitted_decay_rate_s_inv.sqrt()])
    {
        let samples = 4096;
        let integral: f64 = (0..samples)
            .map(|i| {
                let theta = core::f64::consts::PI * (i as f64 + 0.5) / (samples as f64);
                (root * theta.sin()).cos()
            })
            .sum::<f64>()
            / (samples as f64);
        close(integral, 0.0, 4e-15);
    }
}

#[test]
fn harmonic_response_converges_to_closed_form_womersley_oracles_inside_tail_bound() {
    // 120-digit independent power-series evaluations of:
    // circle: [1 - 2 I1(q)/(q I0(q))]/(i*omega),
    // slit:   [1 - tanh(q/2)/(q/2)]/(i*omega), q=sqrt(i*omega).
    // Entries below are divided by exact DC mobility; nu=length=1.
    let circle_oracle = [
        (1.0, 0.9721847306857618, -0.16186236141036338),
        (10.0, 0.2795497501957999, -0.4257706707249414),
        (100.0, 0.010499846222102805, -0.06866996253216768),
    ];
    let slit_oracle = [
        (1.0, 0.9899837629725949, -0.09898528459806184),
        (10.0, 0.500610428266516, -0.4941524395021566),
        (100.0, 0.016929819602278855, -0.1030293605464417),
    ];
    for (circle, oracle) in [(true, circle_oracle), (false, slit_oracle)] {
        for (omega, re, im) in oracle {
            let mut previous_error = f64::INFINITY;
            for count in [1, 4, 16, 32] {
                let model = spectrum(circle, count);
                let dc = model.dc_mobility_s().unwrap();
                let observed = model.harmonic_mobility_s(omega).unwrap();
                assert!(observed.re >= 0.0 && observed.im < 0.0);
                let error =
                    ((observed.re / dc - re).powi(2) + (observed.im / dc - im).powi(2)).sqrt();
                let bound = model.harmonic_mobility_error_bound_s(omega).unwrap() / dc;
                assert!(
                    error <= bound + 3e-14,
                    "circle={circle}, N={count}, omega={omega}, error={error}, bound={bound}"
                );
                assert!(error < previous_error);
                previous_error = error;
            }
            assert!(previous_error < 2e-8);
        }
        let model = spectrum(circle, 32);
        let dc = model.harmonic_mobility_s(0.0).unwrap();
        close(dc.re, model.dc_mobility_s().unwrap(), 0.0);
        assert_eq!(dc.im, 0.0);
        assert_eq!(model.harmonic_mobility_error_bound_s(0.0).unwrap(), 0.0);
        let omega = 1e15;
        let high = model.harmonic_mobility_s(omega).unwrap();
        close(-omega * high.im, 1.0, 4e-15);
        assert!(omega * high.re < 1e-10);
    }
}

#[test]
fn profile_storage_converges_to_parabolic_energy_and_startup_upper_approximation() {
    for (circle, exact_shape_factor) in [(true, 4.0 / 3.0), (false, 6.0 / 5.0)] {
        let mut previous_error = f64::INFINITY;
        for count in [1, 4, 16, 32] {
            let model = spectrum(circle, count);
            let dc = model.dc_mobility_s().unwrap();
            let shape = model
                .modes()
                .iter()
                .map(|m| m.weight / (m.decay_rate_s_inv * dc).powi(2))
                .sum::<f64>();
            assert!(shape > 1.0 && shape < exact_shape_factor);
            let error = exact_shape_factor - shape;
            assert!(error < previous_error);
            previous_error = error;
        }
        assert!(previous_error < 1e-8);
    }
    for time in [1e-4, 0.01, 0.5] {
        // Independent planar startup series using exact full DC minus the
        // exponentially decaying remainder. The uncomputed remainder after
        // 512 modes is negligible at the smallest tested time.
        let exact = 1.0 / 12.0
            - (0..512)
                .map(|m| {
                    let pole = ((2 * m + 1) as f64 * core::f64::consts::PI).powi(2);
                    8.0 / pole.powi(2) * (-pole * time).exp()
                })
                .sum::<f64>();
        let mut previous = f64::INFINITY;
        for count in [1, 4, 16, 32] {
            let model = spectrum(false, count);
            let approximation = model
                .modes()
                .iter()
                .map(|m| m.weight / m.decay_rate_s_inv * (-(-m.decay_rate_s_inv * time).exp_m1()))
                .sum::<f64>();
            assert!(approximation >= exact - 1e-16);
            assert!(approximation <= previous + 1e-16);
            previous = approximation;
        }
    }
}

#[test]
fn prepared_sources_preserve_exact_work_heat_and_wall_impulse() {
    for s in [0.0, 1e-12, 0.01, 0.124, 0.125, 1.0, 7.0, 100.0, 1e12] {
        let factors = NormalViscousFactors::try_new(s, 1.0).unwrap();
        for (initial, forcing_increment) in [(0.7, -0.9), (-0.7, 0.9), (0.0, 1.3), (1.2, 0.0)] {
            let end = factors.decay * initial + factors.phi1 * forcing_increment;
            let mean = factors.phi1 * initial + factors.phi2 * forcing_increment;
            let heat = (factors.heat_initial_root
                * (initial + factors.heat_shift * forcing_increment))
                .powi(2)
                + (factors.heat_forced_root * forcing_increment).powi(2);
            let wall_impulse = -(factors.loss * initial + factors.drag_forced * forcing_increment);
            assert!(heat >= 0.0);
            close(
                0.5 * (end * end - initial * initial),
                forcing_increment * mean - heat,
                3e-15,
            );
            close(end - initial, forcing_increment + wall_impulse, 8e-16);
            if s > 0.0 && s <= 7.0 {
                // Independent Simpson integral of rate*z(t)^2 (4096 panels),
                // not the completion-of-square formula used by production.
                let panels = 4096;
                let mut integral = 0.0;
                for i in 0..=panels {
                    let t = i as f64 / panels as f64;
                    let z = initial * (-s * t).exp() + forcing_increment * (-(-s * t).exp_m1()) / s;
                    let weight = if i == 0 || i == panels {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    integral += weight * s * z * z;
                }
                close(heat, integral / (3 * panels) as f64, 3e-12);
            }
        }
    }
    let model = spectrum(true, 16);
    let prepared = model.prepare(0.07).unwrap();
    assert_eq!(prepared.duration_s(), 0.07);
    for (source, mode) in prepared.modes().iter().zip(model.modes()) {
        assert_eq!(source.weight, mode.weight);
        assert_eq!(
            source.factors,
            NormalViscousFactors::try_new(mode.decay_rate_s_inv, 0.07).unwrap()
        );
    }
    // Roots preserve representable heat coefficients when their squares would
    // underflow. A finite rate*duration may also overflow without NaN factors.
    let tiny = NormalViscousFactors::try_new(f64::from_bits(1), 1.0).unwrap();
    assert!(tiny.heat_forced_root > 0.0);
    let product_underflow = NormalViscousFactors::try_new(1e-300, 1e-300).unwrap();
    assert!(product_underflow.heat_initial_root > 0.0 && product_underflow.heat_forced_root > 0.0);
    let stiff = NormalViscousFactors::try_new(1e300, 1e300).unwrap();
    assert_eq!(stiff.decay, 0.0);
    assert_eq!(stiff.loss, 1.0);
    assert!(stiff.heat_forced_root > 0.0);
    close(stiff.heat_forced_root / 1e-300, 1.0, 4e-15);
}

#[test]
fn zero_viscosity_is_inviscid_and_invalid_domains_are_explicit() {
    let model = NormalViscousSpectrum::parallel_plates(state(0.0), 1.0, 4).unwrap();
    assert!(model.modes().iter().all(|m| m.decay_rate_s_inv == 0.0));
    assert!(model.dc_mobility_s().is_err());
    assert!(model.harmonic_mobility_s(0.0).is_err());
    let response = model.harmonic_mobility_s(2.0).unwrap();
    assert_eq!(response.re, 0.0);
    close(response.im, -0.5, 2e-16);
    assert_eq!(model.harmonic_mobility_error_bound_s(2.0).unwrap(), 0.0);
    for mode in model.prepare(1e200).unwrap().modes() {
        assert_eq!(mode.factors.decay, 1.0);
        assert_eq!(mode.factors.phi1, 1.0);
        assert_eq!(mode.factors.phi2, 0.5);
        assert_eq!(mode.factors.heat_initial_root, 0.0);
        assert_eq!(mode.factors.heat_forced_root, 0.0);
    }
    for count in [0, 33, usize::MAX] {
        assert!(NormalViscousSpectrum::circular_capillary(state(1.0), 1.0, count).is_err());
        assert!(NormalViscousSpectrum::parallel_plates(state(1.0), 1.0, count).is_err());
    }
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(NormalViscousSpectrum::circular_capillary(state(1.0), invalid, 1).is_err());
        assert!(NormalViscousSpectrum::parallel_plates(state(1.0), invalid, 1).is_err());
        assert!(model.prepare(invalid).is_err());
        assert!(model.harmonic_mobility_s(invalid).is_err());
        assert!(model.harmonic_mobility_error_bound_s(invalid).is_err());
        assert!(NormalViscousFactors::try_new(invalid, 1.0).is_err());
    }
    assert!(NormalViscousSpectrum::parallel_plates(state(1.0), 1e-154, 32).is_err());
    assert_eq!(NormalViscousFactors::try_new(3.0, 0.0).unwrap().phi2, 0.5);
}
