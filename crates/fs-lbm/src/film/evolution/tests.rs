//! Independent conservative/energetic and constitutive references for the
//! locally parallel film model. Synthetic coefficients are not calibration.

use super::*;
use fs_material::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn geometry(boundary: ThinFilmBoundary, normal: f64) -> ThinFilmGeometry {
    ThinFilmGeometry {
        spacing_m: 0.5,
        normal_acceleration_m_s2: normal,
        boundary,
    }
}
fn material(sigma: f64) -> FilmMaterial {
    FilmMaterial {
        density_kg_m3: 2.0,
        kinematic_viscosity_m2_s: 0.25,
        surface_tension_n_m: sigma,
    }
}
fn limits() -> ThinFilmLimits {
    ThinFilmLimits {
        minimum_thickness_m: 0.005,
        maximum_thickness_m: 0.2,
        characteristic_length_m: 4.0,
        maximum_height_to_length: 0.06,
        maximum_slope: 0.25,
        maximum_reduced_reynolds: 0.05,
    }
}
fn controls() -> ThinFilmControls {
    ThinFilmControls {
        maximum_newton_iterations: 48,
        maximum_endpoint_evaluations: 4096,
        maximum_profile_evaluations: 1_000_000,
        height_scale_m: 0.1,
        finite_difference_step: 1e-5,
        cell_tolerance_m: 2e-12,
        volume_tolerance_m2: 1e-11,
        energy_tolerance_j_per_m: 1e-11,
        ..ThinFilmControls::default()
    }
}
fn source(exponent: f64) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(0.5, 0.05, 8.0, exponent, 2.0, [0.0, 1e4]).unwrap(),
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
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
// Independent quadratic energy; closed ends have no fictitious wrap link.
fn energy(h: &[f64], rho: f64, normal: f64, sigma: f64, dx: f64, periodic: bool) -> f64 {
    let mut value = 0.5 * rho * normal * dx * h.iter().map(|v| v * v).sum::<f64>();
    for i in 0..h.len() - 1 {
        value += 0.5 * sigma * (h[i + 1] - h[i]).powi(2) / dx;
    }
    if periodic {
        value += 0.5 * sigma * (h[0] - h[h.len() - 1]).powi(2) / dx;
    }
    value
}

#[test]
fn periodic_gravity_and_capillary_decay_match_the_independent_discrete_fourier_rate() {
    let (n, h0, amplitude, dx) = (12, 0.1, 1e-5, 0.5);
    let shape: Vec<_> = (0..n)
        .map(|i| (2.0 * std::f64::consts::PI * i as f64 / n as f64).cos())
        .collect();
    let heights: Vec<_> = shape.iter().map(|&v| h0 + amplitude * v).collect();
    // Separate gravity, pure capillarity, and their combined linear rate.
    for (normal, sigma) in [(1.5, 0.0), (0.0, 0.04), (1.5, 0.04)] {
        let m = material(sigma);
        let k2 = 4.0 * (std::f64::consts::PI / n as f64).sin().powi(2) / (dx * dx);
        let mu = m.density_kg_m3 * m.kinematic_viscosity_m2_s;
        let rate = h0 * h0 * h0 / (3.0 * mu) * (m.density_kg_m3 * normal * k2 + sigma * k2 * k2);
        let dt = 0.4 / rate;
        let mut line = ThinFilmLine::newtonian(
            m,
            geometry(ThinFilmBoundary::Periodic, normal),
            heights.clone(),
            0.0,
        )
        .unwrap();
        assert_eq!(line.face_count(), n);
        let step = line
            .advance(
                dt,
                0.0,
                &vec![0.0; line.face_count()],
                limits(),
                controls(),
                &heights,
                || false,
            )
            .unwrap();
        let measured = 2.0 / n as f64
            * line
                .heights_m()
                .iter()
                .zip(&shape)
                .map(|(&h, &c)| (h - h0) * c)
                .sum::<f64>();
        let expected = amplitude / (1.0 + rate * dt);
        near(measured, expected, 3e-4 * amplitude);
        near(
            dx * line
                .heights_m()
                .iter()
                .zip(&heights)
                .map(|(&a, &b)| a - b)
                .sum::<f64>(),
            0.0,
            1e-11,
        );
        assert!(step.viscous_heat_j_per_m > 0.0);
        assert!(step.stored_energy_change_j_per_m < 0.0);
        assert!(step.backward_euler_energy_defect_j_per_m >= 0.0);
        near(step.energy_residual_j_per_m, 0.0, 1e-11);
        assert!(step.maximum_cell_residual_m <= controls().cell_tolerance_m);
    }
}

#[test]
fn closed_passive_leveling_conserves_volume_and_separates_be_defect_from_heat() {
    let m = material(0.04);
    let g = geometry(ThinFilmBoundary::Closed, 1.5);
    let initial = vec![0.10, 0.12, 0.14, 0.09, 0.11, 0.08];
    let mut line = ThinFilmLine::newtonian(m, g, initial.clone(), 0.0).unwrap();
    for target in [20.0, 40.0] {
        let old = line.heights_m().to_vec();
        let step = line
            .advance(
                target,
                0.0,
                &vec![0.0; line.face_count()],
                limits(),
                controls(),
                &old,
                || false,
            )
            .unwrap();
        let next = line.heights_m();
        let delta: Vec<_> = next.iter().zip(&old).map(|(&a, &b)| a - b).collect();
        let change =
            energy(next, 2.0, 1.5, 0.04, 0.5, false) - energy(&old, 2.0, 1.5, 0.04, 0.5, false);
        let defect = energy(&delta, 2.0, 1.5, 0.04, 0.5, false);
        near(0.5 * delta.iter().sum::<f64>(), 0.0, 1e-11);
        near(
            step.volume_change_m2,
            0.5 * delta.iter().sum::<f64>(),
            2e-14,
        );
        near(step.stored_energy_change_j_per_m, change, 3e-14);
        near(step.backward_euler_energy_defect_j_per_m, defect, 3e-14);
        assert!(change < 0.0 && defect > 0.0 && step.viscous_heat_j_per_m > 0.0);
        near(change + defect + step.viscous_heat_j_per_m, 0.0, 1e-11);
        near(step.body_work_j_per_m, 0.0, 0.0);
        near(step.surface_work_j_per_m, 0.0, 0.0);
    }
    near(line.accepted_time_s(), 40.0, 0.0);
    near(line.heights_m().iter().sum(), initial.iter().sum(), 2e-11);
}

// No slit/film solver enters this reference: bisection of the complete original
// constitutive evaluator plus composite Simpson integration in surface depth.
fn original_profile(
    law: ThermalCarreauYasuda,
    temperature: f64,
    h: f64,
    force: f64,
    tau: f64,
) -> [f64; 3] {
    let n = 512;
    let mut result = [0.0; 3];
    for i in 0..=n {
        let z = h * i as f64 / n as f64;
        let stress = tau + force * z;
        let (mut lo, mut hi) = (0.0, 20.0);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if mid == lo || mid == hi {
                break;
            }
            let value = mid
                * law
                    .evaluate(temperature, mid)
                    .unwrap()
                    .dynamic_viscosity_pa_s();
            if value < stress.abs() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let rate = (0.5 * (lo + hi)).copysign(stress);
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in result.iter_mut().zip([z * rate, rate, stress * rate]) {
            *sum += w * value * h / (3.0 * n as f64);
        }
    }
    result
}

#[test]
fn original_carreau_sources_recover_newtonian_motion_and_resolve_true_thermal_shear_thinning() {
    let heights = vec![0.10, 0.11, 0.12, 0.09, 0.08, 0.10];
    let g = geometry(ThinFilmBoundary::Periodic, 1.5);
    for temperature in [300.0, 330.0] {
        let law = source(1.0);
        let mu = law
            .evaluate(temperature, 0.0)
            .unwrap()
            .dynamic_viscosity_pa_s();
        let mut linear = ThinFilmLine::newtonian(
            FilmMaterial {
                kinematic_viscosity_m2_s: mu / 2.0,
                ..material(0.04)
            },
            g,
            heights.clone(),
            0.0,
        )
        .unwrap();
        let mut original = ThinFilmLine::carreau(
            law,
            &vec![temperature; heights.len()],
            2.0,
            0.04,
            20.0,
            g,
            heights.clone(),
            0.0,
        )
        .unwrap();
        let shear = vec![0.03; linear.face_count()];
        let a = linear
            .advance(0.4, 0.3, &shear, limits(), controls(), &heights, || false)
            .unwrap();
        let b = original
            .advance(0.4, 0.3, &shear, limits(), controls(), &heights, || false)
            .unwrap();
        for (&x, &y) in linear.heights_m().iter().zip(original.heights_m()) {
            near(x, y, 8e-11);
        }
        for (x, y) in a.faces.iter().zip(&b.faces) {
            near(x.flow_per_width_m2_s, y.flow_per_width_m2_s, 2e-10);
            near(x.viscous_dissipation_w_m2, y.viscous_dissipation_w_m2, 2e-9);
        }
    }
    let mut fluxes = Vec::new();
    for temperature in [300.0, 330.0] {
        let law = source(0.4);
        let flat = vec![0.1; 4];
        let mut line = ThinFilmLine::carreau(
            law,
            &[temperature; 4],
            2.0,
            0.0,
            20.0,
            geometry(ThinFilmBoundary::Periodic, 0.0),
            flat.clone(),
            0.0,
        )
        .unwrap();
        let report = line
            .advance(0.1, 2.0, &[0.1; 4], limits(), controls(), &flat, || false)
            .unwrap();
        let reference = original_profile(law, temperature, 0.1, 4.0, 0.1);
        for face in &report.faces {
            near(face.thickness_m, 0.1, 2e-12);
            near(face.flow_per_width_m2_s, reference[0], 2e-9 * reference[0]);
            near(face.surface_velocity_m_s, reference[1], 2e-9 * reference[1]);
            near(
                face.viscous_dissipation_w_m2,
                reference[2],
                2e-9 * reference[2],
            );
        }
        near(4.0 * reference[0] + 0.1 * reference[1], reference[2], 2e-12);
        fluxes.push(report.faces[0].flow_per_width_m2_s);
    }
    assert!(fluxes[1] > 1.05 * fluxes[0]);
}

#[test]
fn opposing_body_and_surface_forces_select_the_interior_godunov_extremum_under_reversal() {
    let hstar = 0.02;
    let heights = vec![0.01, 0.03, 0.01, 0.03];
    let mu = 0.5;
    for sign in [-1.0, 1.0] {
        let force = 2.0 * sign;
        let tau = -force * hstar;
        let mut line = ThinFilmLine::newtonian(
            material(0.0),
            geometry(ThinFilmBoundary::Periodic, 0.0),
            heights.clone(),
            0.0,
        )
        .unwrap();
        let report = line
            .advance(
                0.1,
                force / 2.0,
                &[tau; 4],
                limits(),
                controls(),
                &heights,
                || false,
            )
            .unwrap();
        let expected_q = -force * hstar * hstar * hstar / (6.0 * mu);
        let expected_u = -force * hstar * hstar / (2.0 * mu);
        let mut extrema = 0;
        for face in &report.faces {
            if (face.thickness_m - hstar).abs() < 1e-12 {
                extrema += 1;
                near(face.flow_per_width_m2_s, expected_q, 2e-14);
                near(face.surface_velocity_m_s, expected_u, 2e-14);
                near(face.wall_shear_pa, 0.0, 2e-14);
                near(
                    face.viscous_dissipation_w_m2,
                    force * expected_q + tau * expected_u,
                    2e-14,
                );
            }
        }
        // Both increasing and decreasing Riemann states exist; reversing
        // force switches min to max and changes the active pair of faces.
        assert_eq!(extrema, 2);
        near(report.volume_change_m2, 0.0, 1e-11);
        near(report.energy_residual_j_per_m, 0.0, 1e-11);
    }
}

#[test]
fn invalid_sources_work_limits_and_final_cancellation_keep_geometry_history_atomic() {
    let initial = vec![0.10, 0.12, 0.09, 0.11];
    let g = geometry(ThinFilmBoundary::Periodic, 1.5);
    let make = || ThinFilmLine::newtonian(material(0.04), g, initial.clone(), 0.0).unwrap();
    assert!(
        ThinFilmLine::carreau(
            source(0.4),
            &[300.0, 300.0, 300.0, 361.0],
            2.0,
            0.04,
            20.0,
            g,
            initial.clone(),
            0.0
        )
        .is_err()
    );
    let mut victim = make();
    assert!(
        victim
            .advance(
                0.4,
                0.3,
                &[0.03; 4],
                limits(),
                controls(),
                &[f64::NAN; 4],
                || false
            )
            .is_err()
    );
    assert!(
        victim
            .advance(
                0.4,
                0.3,
                &[0.03; 4],
                limits(),
                ThinFilmControls {
                    maximum_endpoint_evaluations: 1,
                    ..controls()
                },
                &initial,
                || false
            )
            .is_err()
    );
    assert_eq!(victim.heights_m(), initial);
    assert_eq!(victim.accepted_time_s(), 0.0);
    let mut limited = ThinFilmLine::carreau(
        source(0.4),
        &[300.0; 4],
        2.0,
        0.04,
        0.01,
        g,
        initial.clone(),
        0.0,
    )
    .unwrap();
    assert!(
        limited
            .advance(0.4, 2.0, &[0.1; 4], limits(), controls(), &initial, || {
                false
            })
            .is_err()
    );
    assert_eq!(limited.heights_m(), initial);
    assert_eq!(limited.accepted_time_s(), 0.0);
    let mut reference = make();
    let mut polls = 0;
    let expected = reference
        .advance(0.4, 0.3, &[0.03; 4], limits(), controls(), &initial, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(
        victim
            .advance(0.4, 0.3, &[0.03; 4], limits(), controls(), &initial, || {
                calls += 1;
                calls == polls
            })
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(victim.heights_m(), initial);
    assert_eq!(victim.accepted_time_s(), 0.0);
    let accepted = victim
        .advance(0.4, 0.3, &[0.03; 4], limits(), controls(), &initial, || {
            false
        })
        .unwrap();
    assert_eq!(victim.heights_m(), reference.heights_m());
    assert_eq!(victim.accepted_time_s(), reference.accepted_time_s());
    assert_eq!(
        accepted.viscous_heat_j_per_m.to_bits(),
        expected.viscous_heat_j_per_m.to_bits()
    );
    assert_eq!(
        accepted.energy_residual_j_per_m.to_bits(),
        expected.energy_residual_j_per_m.to_bits()
    );
}
