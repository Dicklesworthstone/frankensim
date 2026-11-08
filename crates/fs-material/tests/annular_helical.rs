//! G0/G1/G3: one invariant viscosity for coupled cylindrical shear.
use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::annular::helical::{
    CarreauHelicalProfile, HelicalCarreauControls,
};
use fs_material::fluid_transport::rheology::annular::{
    CarreauAnnularControls, LaminarCarreauAnnulus,
};
use fs_material::fluid_transport::rheology::pipe::CarreauPipeControls;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const A: f64 = 0.01;
const B: f64 = 0.04;
const LENGTH: f64 = 0.5;

fn law(constant: bool, mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    let p = if constant {
        [3.0, 3.0, 0.6, 0.5, 1.0]
    } else {
        [12.0, 0.0, 0.6, 0.5, 1.0]
    };
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(p[0], p[1], p[2], p[3], p[4], [0.0, 1e9]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        mode,
    )
}
fn owner(constant: bool) -> LaminarCarreauAnnulus {
    LaminarCarreauAnnulus::new(
        law(constant, ViscosityShiftMode::VerticalAndRate),
        300.0,
        A,
        B,
        LENGTH,
        1e5,
    )
    .unwrap()
}
fn controls() -> HelicalCarreauControls {
    HelicalCarreauControls {
        local: CarreauPipeControls {
            relative_constitutive_tolerance: 1e-13,
            relative_quadrature_tolerance: 2e-11,
            ..CarreauPipeControls::default()
        },
        relative_velocity_tolerance: 1e-10,
        ..HelicalCarreauControls::default()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1e-15),
        "actual {actual:e}, expected {expected:e}"
    );
}
fn simpson(f: impl Fn(f64) -> f64, lower: f64, upper: f64) -> f64 {
    let panels = 8192;
    let mut sum = 0.0;
    for i in 0..=panels {
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * f(lower + (upper - lower) * f64::from(i) / f64::from(panels));
    }
    (upper - lower) * sum / (3.0 * f64::from(panels))
}
// Independent algebraic inverse for mu=12/sqrt(1+0.6*gamma).
fn inverse(stress: f64) -> f64 {
    let z = stress * 0.6 / 12.0;
    (z * z + z * (z * z + 4.0).sqrt()) / 1.2
}
fn oracle(r: f64, inner: [f64; 2], pressure_stress: f64) -> ([f64; 2], [f64; 2]) {
    let g = 2.0 * pressure_stress / (B - A);
    let tau = [
        A / r * inner[0] - g * (r * r - A * A) / (2.0 * r),
        (A / r).powi(2) * inner[1],
    ];
    let norm = tau[0].hypot(tau[1]);
    let gamma = if norm == 0.0 {
        [0.0; 2]
    } else {
        tau.map(|v| v / norm * inverse(norm))
    };
    (tau, gamma)
}
fn generalized_forces(profile: &CarreauHelicalProfile) -> [f64; 2] {
    let f = profile.flow();
    [
        0.5 * (f.axial_wall_forces_n[1] - f.axial_wall_forces_n[0]),
        f.wall_torques_n_m[1] / A,
    ]
}

#[test]
fn g1_newtonian_helical_flow_recovers_independent_axial_and_rotational_solutions() {
    let owner = owner(true);
    let mu = 3.0;
    let log = (B / A).ln();
    let area = core::f64::consts::PI * (B * B - A * A);
    let axial_resistance = core::f64::consts::TAU * mu * LENGTH / log;
    let torque_resistance =
        4.0 * core::f64::consts::PI * mu * LENGTH * A * A * B * B / (B * B - A * A);
    let conductance = core::f64::consts::PI / (8.0 * mu * LENGTH)
        * (B.powi(4) - A.powi(4) - (B * B - A * A).powi(2) / log);
    for pressure in [-2000.0, 0.0, 2000.0] {
        let axial = [-0.2, 0.8];
        let angular = [11.0, -7.0];
        let du = axial[1] - axial[0];
        let dw = angular[1] - angular[0];
        let profile = owner
            .helical_at_pressure_and_walls(pressure, axial, angular, controls())
            .unwrap();
        let f = profile.flow();
        let c = mu * du / log + pressure * (B * B - A * A) / (4.0 * LENGTH * log);
        let flow = area * axial[0]
            + core::f64::consts::PI * du * (B * B - (B * B - A * A) / (2.0 * log))
            + conductance * pressure;
        let heat = conductance * pressure * pressure
            + axial_resistance * du * du
            + torque_resistance * dw * dw;
        near(f.volumetric_flow_m3_s, flow, 1e-9);
        near(f.dissipated_power_w, heat, 2e-9);
        near(f.wall_torques_n_m[1], torque_resistance * dw, 1e-9);
        near(f.wall_resistance_n_s_m[0][0], axial_resistance, 1e-9);
        near(
            f.wall_resistance_n_s_m[1][1],
            torque_resistance / (A * A),
            1e-9,
        );
        assert_eq!(f.wall_resistance_n_s_m[0][1], 0.0);
        assert!(f.work_residual_w.abs() <= 2e-9 * heat);
        near(
            f.axial_wall_forces_n.iter().sum::<f64>(),
            -pressure * area,
            1e-12,
        );
        for fraction in [0.0, 0.17, 0.5, 0.83, 1.0] {
            let r = A + fraction * (B - A);
            let s = profile.sample_gap_fraction(fraction).unwrap();
            let u =
                axial[0] + c / mu * (r / A).ln() - pressure * (r * r - A * A) / (4.0 * mu * LENGTH);
            let omega = angular[0] + dw * (1.0 - A * A / (r * r)) / (1.0 - A * A / (B * B));
            let gamma = [
                (c / r - pressure * r / (2.0 * LENGTH)) / mu,
                2.0 * dw * A * A * B * B / ((B * B - A * A) * r * r),
            ];
            assert!((s.axial_velocity_m_s - u).abs() < 2e-9);
            assert!((s.angular_velocity_rad_s - omega).abs() < 2e-8);
            near(s.shear_rate_per_s, gamma[0].hypot(gamma[1]), 1e-9);
            near(
                s.dissipation_w_m3,
                mu * (gamma[0] * gamma[0] + gamma[1] * gamma[1]),
                2e-9,
            );
        }
    }
}

#[test]
fn g1_nonlinear_helical_flow_recovers_independent_vector_stress_inverse_and_volume_heat() {
    let owner = owner(false);
    let inner = [60.0, -30.0];
    let p = 20.0;
    let pressure = 2.0 * LENGTH * p / (B - A);
    let du = simpson(|r| oracle(r, inner, p).1[0], A, B);
    let dw = simpson(|r| oracle(r, inner, p).1[1] / r, A, B);
    let axial = [0.2, 0.2 + du];
    let angular = [1.0, 1.0 + dw];
    let profile = owner
        .helical_at_pressure_and_walls(pressure, axial, angular, controls())
        .unwrap();
    let f = profile.flow();
    let flow = core::f64::consts::PI * (B * B - A * A) * axial[0]
        + core::f64::consts::PI * simpson(|r| (B * B - r * r) * oracle(r, inner, p).1[0], A, B);
    let heat = core::f64::consts::TAU
        * LENGTH
        * simpson(
            |r| {
                let (tau, gamma) = oracle(r, inner, p);
                r * (tau[0] * gamma[0] + tau[1] * gamma[1])
            },
            A,
            B,
        );
    near(-f.wall_tractions_pa[0][0], inner[0], 2e-8);
    near(-f.wall_tractions_pa[0][1], inner[1], 2e-8);
    near(f.volumetric_flow_m3_s, flow, 2e-8);
    near(f.dissipated_power_w, heat, 3e-8);
    assert!(f.work_residual_w.abs() < 3e-9 * heat);
    assert!(f.newton_iterations > 0 && f.newton_iterations < controls().maximum_newton_iterations);
    for fraction in [0.0, 0.19, 0.53, 0.88, 1.0] {
        let r = A + fraction * (B - A);
        let sample = profile.sample_gap_fraction(fraction).unwrap();
        let (tau, gamma) = oracle(r, inner, p);
        for i in 0..2 {
            near(sample.radial_shear_stress_pa[i], tau[i], 2e-8);
            near(sample.shear_rate_components_per_s[i], gamma[i], 3e-8);
        }
        let expected_u = axial[0] + simpson(|r| oracle(r, inner, p).1[0], A, r);
        let expected_omega = angular[0] + simpson(|r| oracle(r, inner, p).1[1] / r, A, r);
        assert!((sample.axial_velocity_m_s - expected_u).abs() < 1e-9);
        assert!((sample.angular_velocity_rad_s - expected_omega).abs() < 1e-7);
        let material = owner
            .law()
            .evaluate(300.0, gamma[0].hypot(gamma[1]))
            .unwrap();
        near(
            sample.dynamic_viscosity_pa_s,
            material.dynamic_viscosity_pa_s(),
            3e-8,
        );
    }
    // The coupled invariant changes torque relative to a rotation-only profile.
    let separate = owner
        .rotation_at_wall_speeds(angular, CarreauAnnularControls::default())
        .unwrap();
    assert!((separate.rotation().wall_torques_n_m[1] / f.wall_torques_n_m[1] - 1.0).abs() > 1e-3);
}

#[test]
fn g1_coupled_wall_resistance_is_reciprocal_positive_and_matches_independent_differences() {
    let owner = owner(false);
    let pressure = 500.0;
    let axial = [-0.1, 0.3];
    let angular = [-4.0, 8.0];
    let profile = owner
        .helical_at_pressure_and_walls(pressure, axial, angular, controls())
        .unwrap();
    let r = profile.flow().wall_resistance_n_s_m;
    assert_eq!(r[0][1], r[1][0]);
    assert!(r[0][0] > 0.0 && r[1][1] > 0.0 && r[0][0] * r[1][1] > r[0][1] * r[1][0]);
    assert!(r[0][1].abs() > 1e-5);
    let h = 1e-5;
    for column in 0..2 {
        let mut lo_u = axial;
        let mut hi_u = axial;
        let mut lo_w = angular;
        let mut hi_w = angular;
        if column == 0 {
            lo_u[0] += 0.5 * h;
            lo_u[1] -= 0.5 * h;
            hi_u[0] -= 0.5 * h;
            hi_u[1] += 0.5 * h;
        } else {
            lo_w[0] += 0.5 * h / A;
            lo_w[1] -= 0.5 * h / A;
            hi_w[0] -= 0.5 * h / A;
            hi_w[1] += 0.5 * h / A;
        }
        let low = generalized_forces(
            &owner
                .helical_at_pressure_and_walls(pressure, lo_u, lo_w, controls())
                .unwrap(),
        );
        let high = generalized_forces(
            &owner
                .helical_at_pressure_and_walls(pressure, hi_u, hi_w, controls())
                .unwrap(),
        );
        for row in 0..2 {
            near(r[row][column], (high[row] - low[row]) / (2.0 * h), 3e-6);
        }
    }
}

#[test]
fn g3_pure_limits_rigid_motion_reversal_and_common_motion_preserve_physical_heat() {
    let owner = owner(false);
    let pure = CarreauAnnularControls {
        local: controls().local,
        ..CarreauAnnularControls::default()
    };
    let axial = owner
        .profile_at_pressure_and_walls(800.0, [-0.1, 0.2], pure)
        .unwrap();
    let helical = owner
        .helical_at_pressure_and_walls(800.0, [-0.1, 0.2], [7.0, 7.0], controls())
        .unwrap();
    near(
        helical.flow().volumetric_flow_m3_s,
        axial.flow().volumetric_flow_m3_s,
        2e-8,
    );
    near(
        helical.flow().dissipated_power_w,
        axial.flow().dissipated_power_w,
        2e-8,
    );
    assert_eq!(helical.flow().wall_torques_n_m, [-0.0, 0.0]);
    let rotation = owner.rotation_at_wall_speeds([-2.0, 8.0], pure).unwrap();
    let helical = owner
        .helical_at_pressure_and_walls(0.0, [0.4, 0.4], [-2.0, 8.0], controls())
        .unwrap();
    near(
        helical.flow().wall_torques_n_m[1],
        rotation.rotation().wall_torques_n_m[1],
        2e-8,
    );
    near(
        helical.flow().dissipated_power_w,
        rotation.rotation().dissipated_power_w,
        2e-8,
    );
    let rigid = owner
        .helical_at_pressure_and_walls(0.0, [3.0, 3.0], [123.0, 123.0], controls())
        .unwrap();
    assert_eq!(rigid.flow().dissipated_power_w, 0.0);
    for fraction in [0.0, 0.2, 0.7, 1.0] {
        let s = rigid.sample_gap_fraction(fraction).unwrap();
        assert_eq!(s.shear_rate_per_s, 0.0);
        assert_eq!(s.axial_velocity_m_s, 3.0);
        assert_eq!(s.angular_velocity_rad_s, 123.0);
    }
    let direct = owner
        .helical_at_pressure_and_walls(80.0, [0.0, 1.0], [0.0, 10.0], controls())
        .unwrap();
    let shifted = owner
        .helical_at_pressure_and_walls(80.0, [1e6, 1e6 + 1.0], [1e9, 1e9 + 10.0], controls())
        .unwrap();
    let reverse = owner
        .helical_at_pressure_and_walls(-80.0, [0.0, -1.0], [0.0, -10.0], controls())
        .unwrap();
    assert_eq!(
        direct.flow().dissipated_power_w,
        shifted.flow().dissipated_power_w
    );
    assert_eq!(
        direct.flow().mechanical_power_w,
        shifted.flow().mechanical_power_w
    );
    near(
        reverse.flow().dissipated_power_w,
        direct.flow().dissipated_power_w,
        1e-12,
    );
    near(
        reverse.flow().wall_torques_n_m[1],
        -direct.flow().wall_torques_n_m[1],
        1e-12,
    );
}

#[test]
fn g3_both_thermal_shift_modes_retain_the_original_coupled_constitutive_identity() {
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let law = law(false, mode);
        let cold = LaminarCarreauAnnulus::new(law, 300.0, A, B, LENGTH, 1e5).unwrap();
        let warm = LaminarCarreauAnnulus::new(law, 330.0, A, B, LENGTH, 1e5).unwrap();
        let factor = law.temperature_shift().evaluate(330.0).unwrap().factor();
        let cold = cold
            .helical_at_pressure_and_walls(500.0, [0.0, 0.2], [0.0, 9.0], controls())
            .unwrap();
        let (pressure, velocity_scale, stress_scale) = match mode {
            ViscosityShiftMode::Vertical => (factor * 500.0, 1.0, factor),
            ViscosityShiftMode::VerticalAndRate => (500.0, factor.recip(), 1.0),
        };
        let warm = warm
            .helical_at_pressure_and_walls(
                pressure,
                [0.0, 0.2 * velocity_scale],
                [0.0, 9.0 * velocity_scale],
                controls(),
            )
            .unwrap();
        near(
            warm.flow().wall_torques_n_m[1],
            stress_scale * cold.flow().wall_torques_n_m[1],
            2e-8,
        );
        near(
            warm.flow().volumetric_flow_m3_s,
            velocity_scale * cold.flow().volumetric_flow_m3_s,
            2e-8,
        );
        near(
            warm.flow().dissipated_power_w,
            stress_scale * velocity_scale * cold.flow().dissipated_power_w,
            3e-8,
        );
        for i in 0..2 {
            for j in 0..2 {
                near(
                    warm.flow().wall_resistance_n_s_m[i][j],
                    factor * cold.flow().wall_resistance_n_s_m[i][j],
                    3e-8,
                );
            }
        }
    }
}

#[test]
fn g1_thin_gap_uses_vector_shear_and_g0_budgets_refuse_transactionally() {
    let a = 1.0;
    let b = a + 1e-8;
    let gap = b - a;
    let owner = LaminarCarreauAnnulus::new(
        law(false, ViscosityShiftMode::VerticalAndRate),
        300.0,
        a,
        b,
        LENGTH,
        100.0,
    )
    .unwrap();
    let profile = owner
        .helical_at_pressure_and_walls(0.0, [0.0, 2.0 * gap], [0.0, 3.0 * gap / a], controls())
        .unwrap();
    let s = profile.sample_gap_fraction(0.5).unwrap();
    near(s.shear_rate_components_per_s[0], 2.0, 1e-7);
    near(s.shear_rate_components_per_s[1], 3.0, 1e-7);
    near(s.shear_rate_per_s, 13.0_f64.sqrt(), 1e-7);
    let material = owner.law().evaluate(300.0, 13.0_f64.sqrt()).unwrap();
    near(
        s.dynamic_viscosity_pa_s,
        material.dynamic_viscosity_pa_s(),
        1e-7,
    );
    let mu = material.dynamic_viscosity_pa_s();
    let tangent = material.stress_tangent_pa_s();
    let coefficient = core::f64::consts::TAU * LENGTH * a / gap;
    let expected_cross = coefficient * (tangent - mu) * 6.0 / 13.0;
    near(
        profile.flow().wall_resistance_n_s_m[0][1],
        expected_cross,
        2e-7,
    );

    let owner = super_owner();
    let good = || owner.helical_at_pressure_and_walls(500.0, [0.0, 0.2], [0.0, 9.0], controls());
    let original = good().unwrap();
    let mut budget = controls();
    budget.maximum_constitutive_queries = 2;
    assert!(
        owner
            .helical_at_pressure_and_walls(500.0, [0.0, 0.2], [0.0, 9.0], budget)
            .is_err()
    );
    assert!(
        owner
            .helical_at_pressure_and_walls(f64::NAN, [0.0; 2], [0.0; 2], controls())
            .is_err()
    );
    assert!(
        owner
            .helical_at_pressure_and_walls(1e20, [0.0; 2], [0.0; 2], controls())
            .is_err()
    );
    let mut calls = 0;
    let interrupted = owner.helical_with_check(500.0, [0.0, 0.2], [0.0, 9.0], controls(), || {
        calls += 1;
        if calls == 7 {
            Err(MaterialError::State {
                what: "interrupt helical trial".into(),
            })
        } else {
            Ok(())
        }
    });
    assert!(interrupted.is_err());
    assert_eq!(calls, 7);
    assert_eq!(good().unwrap(), original);
    assert!(original.sample_gap_fraction(-0.1).is_err());
    assert!(
        original
            .sample_gap_fraction_with_check(0.5, || Err(MaterialError::State {
                what: "interrupt sample".into()
            }))
            .is_err()
    );
}
fn super_owner() -> LaminarCarreauAnnulus {
    owner(false)
}
