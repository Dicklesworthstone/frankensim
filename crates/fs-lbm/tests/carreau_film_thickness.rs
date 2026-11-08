//! G1/G3: a signed film flux can have two physical thicknesses. The caller
//! selects a monotone interval; original constitutive stress, work and heat
//! remain authoritative throughout the inverse. All parameters are synthetic.

use fs_lbm::film::carreau::{
    CarreauFilmControls, CarreauFilmError, CarreauShearedFilm, CarreauShearedFilmFlowSpec,
    CarreauShearedFilmSpec,
};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const MU: f64 = 0.125;
const TURN: f64 = 1.0 / 1024.0;

fn law(parameters: [f64; 5]) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            parameters[0],
            parameters[1],
            parameters[2],
            parameters[3],
            parameters[4],
            [0.0, 1e6],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}

fn constant_law() -> ThermalCarreauYasuda {
    law([MU, 0.03125, 0.4, 1.0, 2.0])
}

fn controls() -> CarreauFilmControls {
    let mut c = CarreauFilmControls::default();
    c.slit.local.relative_constitutive_tolerance = 1e-12;
    c.slit.local.relative_quadrature_tolerance = 1e-10;
    c.slit.maximum_wall_stress_iterations = 1;
    c.relative_flow_tolerance = 1e-10;
    c
}

fn opposed(flow: f64, interval: [f64; 2]) -> CarreauShearedFilmFlowSpec {
    CarreauShearedFilmFlowSpec {
        density_kg_m3: 1024.0,
        tangential_acceleration_m_s2: 1.0,
        surface_shear_pa: -1.0,
        flow_per_width_m2_s: flow,
        thickness_interval_m: interval,
    }
}

fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(scale),
        "{actual:e} != {expected:e}, scale={scale:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn opposing_forcing_recovers_both_closed_form_thicknesses_and_their_signs() {
    // With x=h/TURN, q/(TURN²/MU)=x³/3-x²/2. For q=-TURN²/(12 MU),
    // the two positive roots are x=1/2 and x=(1+sqrt(3))/2.
    let magnitude = TURN * TURN / (12.0 * MU);
    for orientation in [-1.0, 1.0] {
        for (interval, expected_h, branch_sign) in [
            ([0.25 * TURN, TURN], 0.5 * TURN, -1.0),
            ([TURN, 1.5 * TURN], 0.5 * (1.0 + 3.0_f64.sqrt()) * TURN, 1.0),
        ] {
            let mut spec = opposed(-orientation * magnitude, interval);
            spec.tangential_acceleration_m_s2 *= orientation;
            spec.surface_shear_pa *= orientation;
            let solution =
                CarreauShearedFilm::at_flow(constant_law(), 300.0, spec, 1000.0, controls())
                    .unwrap();
            let film = solution.film();
            near(film.thickness_m(), expected_h, TURN, 2e-9);
            near(
                film.flow_per_width_m2_s(),
                spec.flow_per_width_m2_s,
                magnitude,
                2e-9,
            );
            assert_eq!(solution.spec(), spec);
            assert!(solution.profile_evaluations() <= controls().maximum_thickness_iterations);
            let wall_rate = (spec.surface_shear_pa
                + spec.density_kg_m3 * spec.tangential_acceleration_m_s2 * expected_h)
                / MU;
            near(
                solution.flow_thickness_derivative_m_s(),
                expected_h * wall_rate,
                TURN / MU,
                2e-9,
            );
            assert_eq!(
                solution.flow_thickness_derivative_m_s().signum(),
                orientation * branch_sign
            );
            near(
                solution.thickness_flow_derivative_s_m().unwrap()
                    * solution.flow_thickness_derivative_m_s(),
                1.0,
                1.0,
                2e-14,
            );
            assert!(film.power().viscous_dissipation_w_m2 > 0.0);
            near(
                film.power().body_work_w_m2 + film.power().surface_work_w_m2,
                film.power().viscous_dissipation_w_m2,
                film.power().viscous_dissipation_w_m2,
                2e-12,
            );
        }
    }
}

#[test]
fn zero_flux_keeps_counterflow_and_turning_point_has_no_finite_inverse_tangent() {
    let solution = CarreauShearedFilm::at_flow(
        constant_law(),
        300.0,
        opposed(0.0, [TURN, 2.0 * TURN]),
        1000.0,
        controls(),
    )
    .unwrap();
    let film = solution.film();
    near(film.thickness_m(), 1.5 * TURN, TURN, 2e-9);
    near(film.flow_per_width_m2_s(), 0.0, TURN * TURN / MU, 2e-9);
    assert!(film.surface_velocity_m_s() < 0.0);
    assert!(film.sample(0.25 * film.thickness_m()).unwrap().velocity_m_s > 0.0);
    assert!(film.power().viscous_dissipation_w_m2 > 0.0);
    near(
        film.power().viscous_dissipation_w_m2,
        film.power().surface_work_w_m2,
        film.power().viscous_dissipation_w_m2,
        2e-9,
    );

    let at_turn = CarreauShearedFilm::new(
        constant_law(),
        300.0,
        CarreauShearedFilmSpec {
            density_kg_m3: 1024.0,
            thickness_m: TURN,
            tangential_acceleration_m_s2: 1.0,
            surface_shear_pa: -1.0,
        },
        1000.0,
        controls().slit,
    )
    .unwrap();
    near(
        at_turn.flow_per_width_m2_s(),
        -TURN * TURN / (6.0 * MU),
        TURN * TURN / MU,
        2e-14,
    );
    for interval in [[0.5 * TURN, TURN], [TURN, 2.0 * TURN]] {
        let solution = CarreauShearedFilm::at_flow(
            constant_law(),
            300.0,
            opposed(at_turn.flow_per_width_m2_s(), interval),
            1000.0,
            controls(),
        )
        .unwrap();
        assert_eq!(solution.film().thickness_m(), TURN);
        assert_eq!(solution.profile_evaluations(), 2);
        assert_eq!(solution.flow_thickness_derivative_m_s(), 0.0);
        assert!(solution.thickness_flow_derivative_s_m().is_none());
    }
}

#[test]
fn pure_traction_inverse_retains_shifted_yasuda_material_and_analytic_sensitivity() {
    let source = law([0.08, 0.006, 0.4, 0.45, 0.7]);
    let (temperature, rate, h) = (340.0, 30.0, 0.003);
    let mu = source
        .evaluate(temperature, rate)
        .unwrap()
        .dynamic_viscosity_pa_s();
    for sign in [-1.0, 1.0] {
        let spec = CarreauShearedFilmFlowSpec {
            density_kg_m3: 950.0,
            tangential_acceleration_m_s2: 0.0,
            surface_shear_pa: sign * mu * rate,
            flow_per_width_m2_s: sign * 0.5 * h * h * rate,
            thickness_interval_m: [0.001, 0.006],
        };
        let solved =
            CarreauShearedFilm::at_flow(source, temperature, spec, 1000.0, controls()).unwrap();
        near(solved.film().thickness_m(), h, h, 2e-9);
        near(
            solved.film().surface_velocity_m_s(),
            sign * h * rate,
            h * rate,
            2e-9,
        );
        near(
            solved.flow_thickness_derivative_m_s(),
            sign * h * rate,
            h * rate,
            2e-9,
        );
        near(
            solved.thickness_flow_derivative_s_m().unwrap(),
            sign / (h * rate),
            1.0 / (h * rate),
            2e-9,
        );
        assert_eq!(solved.film().law(), &source);
        assert_eq!(solved.film().temperature_k(), temperature);
        assert_eq!(solved.film().power().body_work_w_m2, 0.0);
    }
}

// Independent stress inversion for mu=0.006+0.074[1+(0.4 gamma)²]^(-0.275).
// Platform powf and bisection do not call the production material inverse.
fn reference_rate(stress: f64) -> f64 {
    let (mut low, mut high) = (0.0, stress.abs() / 0.006);
    for _ in 0..80 {
        let rate = 0.5 * (low + high);
        let mu = 0.006 + 0.074 * (1.0 + (0.4 * rate).powi(2)).powf(-0.275);
        if mu * rate > stress.abs() {
            high = rate;
        } else {
            low = rate;
        }
    }
    (0.5 * (low + high)).copysign(stress)
}

fn reference_integrals(h: f64, force: f64, traction: f64) -> [f64; 3] {
    let panels = 2048;
    let dz = h / panels as f64;
    let mut sums = [0.0; 3];
    for i in 0..=panels {
        let z = i as f64 * dz;
        let stress = traction + force * z;
        let rate = reference_rate(stress);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in sums.iter_mut().zip([rate, z * rate, stress * rate]) {
            *sum += weight * value;
        }
    }
    sums.map(|sum| sum * dz / 3.0)
}

#[test]
fn nonlinear_thickness_profile_and_heat_match_independent_momentum_integrals() {
    let (expected_h, force, traction) = (0.0012, 1900.0, -0.8);
    let exact = reference_integrals(expected_h, force, traction);
    let source = law([0.08, 0.006, 0.4, 0.45, 2.0]);
    let spec = CarreauShearedFilmFlowSpec {
        density_kg_m3: 950.0,
        tangential_acceleration_m_s2: 2.0,
        surface_shear_pa: traction,
        flow_per_width_m2_s: exact[1],
        thickness_interval_m: [0.0005, 0.002],
    };
    let solved = CarreauShearedFilm::at_flow(source, 300.0, spec, 1000.0, controls()).unwrap();
    let film = solved.film();
    near(film.thickness_m(), expected_h, expected_h, 3e-7);
    near(film.surface_velocity_m_s(), exact[0], exact[0].abs(), 3e-7);
    near(film.flow_per_width_m2_s(), exact[1], exact[1].abs(), 3e-7);
    near(
        film.power().viscous_dissipation_w_m2,
        exact[2],
        exact[2],
        3e-7,
    );
    let derivative = expected_h * reference_rate(traction + force * expected_h);
    near(
        solved.flow_thickness_derivative_m_s(),
        derivative,
        derivative.abs(),
        3e-7,
    );
    near(solved.flow_residual_m2_s(), 0.0, exact[1].abs(), 3e-9);
    assert_eq!(
        film.power().body_work_w_m2,
        force * film.flow_per_width_m2_s()
    );
    assert_eq!(film.law(), &source);
}

#[test]
fn ambiguous_branch_unreachable_flux_source_limits_and_exhausted_work_refuse() {
    let target = -TURN * TURN / (12.0 * MU);
    let spec = opposed(target, [0.25 * TURN, TURN]);
    assert!(matches!(
        CarreauShearedFilm::at_flow(
            constant_law(),
            300.0,
            opposed(target, [0.5 * TURN, 1.5 * TURN]),
            1000.0,
            controls()
        ),
        Err(CarreauFilmError::NonMonotoneThicknessInterval { .. })
    ));
    assert!(matches!(
        CarreauShearedFilm::at_flow(
            constant_law(),
            300.0,
            CarreauShearedFilmFlowSpec {
                flow_per_width_m2_s: -target,
                ..spec
            },
            1000.0,
            controls()
        ),
        Err(CarreauFilmError::FlowOutsideThicknessInterval { .. })
    ));
    for flow in [0.0, 1e-6] {
        assert!(
            CarreauShearedFilm::at_flow(
                constant_law(),
                300.0,
                CarreauShearedFilmFlowSpec {
                    flow_per_width_m2_s: flow,
                    tangential_acceleration_m_s2: 0.0,
                    surface_shear_pa: 0.0,
                    ..spec
                },
                1000.0,
                controls()
            )
            .is_err()
        );
    }
    // The root could lie below the rate ceiling, but the declared upper bound
    // does not. The inverse must refuse the supplied domain, not clip it.
    let unsupported = CarreauShearedFilmFlowSpec {
        surface_shear_pa: 0.0,
        flow_per_width_m2_s: 1e-6,
        thickness_interval_m: [1.0 / 4096.0, 1.0 / 256.0],
        ..spec
    };
    assert!(matches!(
        CarreauShearedFilm::at_flow(constant_law(), 300.0, unsupported, 10.0, controls()),
        Err(CarreauFilmError::Material(_))
    ));
    assert!(matches!(
        CarreauShearedFilm::at_flow(constant_law(), 270.0, spec, 1000.0, controls()),
        Err(CarreauFilmError::Material(_))
    ));
    let mut bounded = controls();
    bounded.maximum_thickness_iterations = 2;
    assert!(matches!(
        CarreauShearedFilm::at_flow(constant_law(), 300.0, spec, 1000.0, bounded),
        Err(CarreauFilmError::UnresolvedFlow { iterations: 2, .. })
    ));
    let mut checks = 0;
    let result = CarreauShearedFilm::at_flow_with_cancel(
        constant_law(),
        300.0,
        spec,
        1000.0,
        controls(),
        || {
            checks += 1;
            checks == 20
        },
    );
    assert!(matches!(result, Err(CarreauFilmError::Cancelled)));
    assert_eq!(checks, 20);
}

#[test]
fn admitted_tolerance_retains_actual_flux_and_heat_instead_of_substituting_the_request() {
    let h = 1.0 / 2048.0;
    let fixed = CarreauShearedFilm::new(
        constant_law(),
        300.0,
        CarreauShearedFilmSpec {
            density_kg_m3: 1024.0,
            thickness_m: h,
            tangential_acceleration_m_s2: 1.0,
            surface_shear_pa: 1.0,
        },
        1000.0,
        controls().slit,
    )
    .unwrap();
    let spec = CarreauShearedFilmFlowSpec {
        density_kg_m3: 1024.0,
        tangential_acceleration_m_s2: 1.0,
        surface_shear_pa: 1.0,
        flow_per_width_m2_s: 1.001 * fixed.flow_per_width_m2_s(),
        thickness_interval_m: [h, 1.0 / 512.0],
    };
    let mut budget = controls();
    budget.relative_flow_tolerance = 0.01;
    budget.maximum_thickness_iterations = 2;
    let solved = CarreauShearedFilm::at_flow(constant_law(), 300.0, spec, 1000.0, budget).unwrap();
    assert_eq!(solved.film(), &fixed);
    assert_eq!(solved.profile_evaluations(), 2);
    assert!(solved.flow_residual_m2_s() < 0.0);
    assert_eq!(
        solved.film().power().body_work_w_m2,
        1024.0 * fixed.flow_per_width_m2_s()
    );
    assert_ne!(
        solved.film().power().body_work_w_m2,
        1024.0 * spec.flow_per_width_m2_s
    );
}
