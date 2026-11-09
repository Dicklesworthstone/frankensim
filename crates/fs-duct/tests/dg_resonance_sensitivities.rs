//! G1/G3: geometric resonance/decay and physical-mode tangents, not fitted pitch.
use fs_duct::dg::continuous::poles::{DgResonancePole, PoleOptions, PoleResponseOptions};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgSidePort, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::WallPin;

const WORK: usize = 1_000_000_000;
fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn cylinder(length: f64, elements: usize, radius: f64, outlet: DgBoundary, lossy: bool) -> DgDuct {
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        &air(),
        DgMesh { elements, order: 2 },
        DgBoundary::Port(PortId(7)),
        outlet,
    )
    .unwrap();
    if lossy {
        duct.with_viscothermal(4, [80.0, 1200.0], 0.2).unwrap()
    } else {
        duct
    }
}
fn options(omega: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: omega,
        step_tolerance: 1e-12,
        residual_tolerance: 1e-11,
        max_iterations: 30,
        max_dimension: 1024,
        max_dense_bytes: 128 * 1024 * 1024,
        max_work: 2_000_000_000,
    }
}
fn response_options(omega: f64) -> PoleResponseOptions {
    PoleResponseOptions {
        spectral_scale: omega,
        pole_step_tolerance: 1e-10,
        pole_residual_tolerance: 1e-10,
        linear_residual_tolerance: 1e-11,
        max_dimension: 1024,
        max_vectors: 1,
        max_dense_bytes: 128 * 1024 * 1024,
        max_work: 2_000_000_000,
    }
}
fn pole(duct: &DgDuct, seed: C64) -> DgResonancePole {
    let state = duct
        .harmonic_response(
            0.98 * seed.im.abs(),
            &[(PortId(7), C64::ONE)],
            WORK,
            1e-10,
            &mut || false,
        )
        .unwrap()
        .energy_coordinates;
    duct.resonance_pole(seed, &state, &options(seed.abs()), || false)
        .unwrap()
}
fn near(actual: C64, expected: C64, relative: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < relative * scale.max(expected.abs()),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn lossless_length_derivative_converges_and_transports_physical_fields_and_ports() {
    let length = 0.5;
    let omega = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let mut errors = Vec::new();
    for elements in [2, 4, 8] {
        let duct = cylinder(length, elements, 0.01, DgBoundary::IdealOpen, false)
            .with_side_ports(&[DgSidePort {
                id: PortId(11),
                downstream_element: elements / 2,
            }])
            .unwrap();
        let original = pole(&duct, C64::new(0.0, -omega));
        let prepared = duct
            .prepare_resonance_poles(omega, WORK as u64, || false)
            .unwrap();
        let result = prepared
            .log_length_sensitivity(&original, &response_options(omega), || false)
            .unwrap();
        near(
            result.laplace_pole_derivative_per_s,
            -original.laplace_pole_per_s,
            1e-10,
            omega,
        );
        errors.push((result.laplace_pole_derivative_per_s.im / omega - 1.0).abs());
        let gauge = original.energy_coordinates[result.normalization_pivot];
        for (&derivative, &field) in result
            .pressure_derivative_pa
            .iter()
            .zip(&original.pressure_pa)
        {
            near(
                derivative,
                (field / gauge).scale(-0.5),
                1e-9,
                original.pressure_pa[0].abs() / gauge.abs(),
            );
        }
        for (&derivative, &field) in result
            .volume_flow_derivative_m3_s
            .iter()
            .zip(&original.volume_flow_m3_s)
        {
            near(derivative, (field / gauge).scale(-0.5), 1e-9, 1.0);
        }
        for (derivative, field) in result.port_pressure_derivatives.iter().zip(&original.ports) {
            assert_eq!(derivative.id, field.id);
            near(
                derivative.pressure_pa,
                (field.pressure_pa / gauge).scale(-0.5),
                1e-9,
                1.0,
            );
        }
        assert!(
            result
                .energy_coordinate_derivative
                .iter()
                .all(|x| x.abs() < 1e-9)
        );
        assert!(result.oscillation_frequency_derivative_hz().unwrap() < 0.0);
        assert!(
            result.original_relative_residual < 1e-10
                && result.derivative_relative_residual < 1e-11
        );
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 1e-6, "{errors:?}");
}

// Scalar telegraph dispersion of the retained ORIGINAL local wall law:
// F(s,L)=p(s)q(s)+(pi*c/(2L))². Therefore ds/dln(L)=2*omega0²/F_s,
// which generally changes decay as well as pitch and is not simply -s.
fn telegraph(model: &fs_phs::ZkDiffusive, s: C64) -> (C64, C64) {
    let (mut p, mut q) = (s, s + C64::from_re(model.series_drag_per_s()));
    let (mut dp, mut dq) = (C64::ONE, C64::ONE);
    for (terms, value, derivative) in [
        (model.shunt_terms(), &mut p, &mut dp),
        (model.series_terms(), &mut q, &mut dq),
    ] {
        for term in terms {
            let denominator = s + C64::from_re(term.rate_per_s);
            *value = *value + (s / denominator).scale(term.gain_per_s);
            *derivative = *derivative
                + C64::from_re(term.gain_per_s * term.rate_per_s) / (denominator * denominator);
        }
    }
    (p * q, dp * q + p * dq)
}

#[test]
fn dispersive_length_tangent_converges_and_matches_rebuilt_geometry() {
    let length = 0.3;
    let omega = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let duct = cylinder(length, elements, 0.001, DgBoundary::IdealOpen, true);
        let model = duct.viscothermal_model().unwrap();
        let mut expected = C64::new(-10.0, -omega);
        for _ in 0..40 {
            let (value, derivative) = telegraph(model, expected);
            let step = (value + C64::from_re(omega * omega)) / derivative;
            expected = expected - step;
            if step.abs() < 1e-13 * omega {
                break;
            }
        }
        let expected_derivative = C64::from_re(2.0 * omega * omega) / telegraph(model, expected).1;
        let original = pole(&duct, expected);
        let prepared = duct
            .prepare_resonance_poles(omega, WORK as u64, || false)
            .unwrap();
        let result = prepared
            .log_length_sensitivity(&original, &response_options(omega), || false)
            .unwrap();
        assert!(result.reference_energy_coordinates.len() > 2 * duct.field_dofs());
        errors.push((result.laplace_pole_derivative_per_s - expected_derivative).abs() / omega);
        assert!(
            (result.laplace_pole_derivative_per_s + original.laplace_pole_per_s).abs()
                > 1e-3 * omega
        );
        if elements == 4 {
            let h = 1e-5_f64;
            let lower = cylinder(
                length * (-h).exp(),
                elements,
                0.001,
                DgBoundary::IdealOpen,
                true,
            );
            let upper = cylinder(
                length * h.exp(),
                elements,
                0.001,
                DgBoundary::IdealOpen,
                true,
            );
            let lower = pole(&lower, original.laplace_pole_per_s);
            let upper = pole(&upper, original.laplace_pole_per_s);
            let finite = (upper.laplace_pole_per_s - lower.laplace_pole_per_s).scale(0.5 / h);
            near(result.laplace_pole_derivative_per_s, finite, 2e-6, omega);
            let pivot = result.normalization_pivot;
            let low_gauge = lower.energy_coordinates[pivot];
            let high_gauge = upper.energy_coordinates[pivot];
            for (i, &tangent) in result.energy_coordinate_derivative.iter().enumerate() {
                let finite = (upper.energy_coordinates[i] / high_gauge
                    - lower.energy_coordinates[i] / low_gauge)
                    .scale(0.5 / h);
                near(tangent, finite, 2e-5, 1.0);
            }
            let pscale = original
                .pressure_pa
                .iter()
                .map(|p| p.abs())
                .fold(0.0_f64, f64::max)
                / original.energy_coordinates[pivot].abs();
            for (i, &tangent) in result.pressure_derivative_pa.iter().enumerate() {
                let finite = (upper.pressure_pa[i] / high_gauge - lower.pressure_pa[i] / low_gauge)
                    .scale(0.5 / h);
                near(tangent, finite, 2e-5, pscale);
            }
        }
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 1e-5, "{errors:?}");
}

#[test]
fn length_changes_boundary_decay_but_not_a_local_gas_loaded_wall_resonance() {
    let gas = air();
    let (radius, length) = (0.01, 0.3);
    let z = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let duct = cylinder(
        length,
        4,
        radius,
        DgBoundary::resistive(3.0 * z).unwrap(),
        false,
    );
    let omega = core::f64::consts::PI * gas.sound_speed / length;
    let original = pole(
        &duct,
        C64::new(
            -gas.sound_speed * core::f64::consts::LN_2 / (2.0 * length),
            -omega,
        ),
    );
    let result = duct
        .prepare_resonance_poles(omega, WORK as u64, || false)
        .unwrap()
        .log_length_sensitivity(&original, &response_options(omega), || false)
        .unwrap();
    near(
        result.laplace_pole_derivative_per_s,
        -original.laplace_pole_per_s,
        1e-9,
        omega,
    );
    assert!(original.amplitude_decay_rate_per_s > 0.0);
    assert!(result.amplitude_decay_derivative_per_s() < 0.0);

    // Uniform sealed-cylinder radial vibration has local gas stiffness B*2/r,
    // independent of axial length. Rescaling every generator entry would fail.
    let pin = WallPin {
        surface_density: 0.2,
        stiffness_per_area: 1e7,
        resistance: 300.0,
    };
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        &gas,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_compliant_wall(pin)
    .unwrap();
    let beta = 2.0 * gas.density * gas.sound_speed * gas.sound_speed / radius;
    let decay = pin.resistance / (2.0 * pin.surface_density);
    let omega = ((pin.stiffness_per_area + beta) / pin.surface_density - decay * decay).sqrt();
    let s = C64::new(-decay, -omega);
    let coupling = (beta / pin.surface_density).sqrt();
    let wall_omega = (pin.stiffness_per_area / pin.surface_density).sqrt();
    let m = duct.field_dofs();
    let mut seed: Vec<_> = duct
        .to_energy_coords(&vec![1.0; m], &vec![0.0; m])
        .into_iter()
        .map(C64::from_re)
        .collect();
    for i in 0..m {
        seed[2 * m + i] = seed[i].scale(-wall_omega / coupling);
        seed[3 * m + i] = -s * seed[i].scale(1.0 / coupling);
    }
    let original = duct
        .resonance_pole(s, &seed, &options(omega), || false)
        .unwrap();
    let result = duct
        .prepare_resonance_poles(omega, WORK as u64, || false)
        .unwrap()
        .log_length_sensitivity(&original, &response_options(omega), || false)
        .unwrap();
    assert!(result.laplace_pole_derivative_per_s.abs() < 1e-9 * omega);
}

#[test]
fn sensitivity_is_gauge_invariant_and_refuses_foreign_modes_limits_and_cancellation() {
    let duct = cylinder(0.5, 2, 0.01, DgBoundary::IdealOpen, false);
    let omega = core::f64::consts::PI * air().sound_speed;
    let original = pole(&duct, C64::new(0.0, -omega));
    let prepared = duct
        .prepare_resonance_poles(omega, WORK as u64, || false)
        .unwrap();
    let limits = response_options(omega);
    let reference = prepared
        .log_length_sensitivity(&original, &limits, || false)
        .unwrap();
    let mut changed = original.clone();
    for x in &mut changed.energy_coordinates {
        *x = *x * C64::new(1e100, -3e100);
    }
    changed.pressure_pa.fill(C64::new(f64::NAN, 0.0));
    changed.relative_residual = f64::NAN;
    let result = prepared
        .log_length_sensitivity(&changed, &limits, || false)
        .unwrap();
    near(
        result.laplace_pole_derivative_per_s,
        reference.laplace_pole_derivative_per_s,
        1e-11,
        omega,
    );
    changed.energy_coordinates[0] = C64::ZERO;
    assert!(
        prepared
            .log_length_sensitivity(&changed, &limits, || false)
            .is_err()
    );
    assert!(
        prepared
            .log_length_sensitivity(&original, &limits, || true)
            .is_err()
    );
    assert!(
        prepared
            .log_length_sensitivity(
                &original,
                &PoleResponseOptions {
                    max_work: 1,
                    ..limits
                },
                || false
            )
            .is_err()
    );
    assert!(
        prepared
            .log_length_sensitivity(
                &original,
                &PoleResponseOptions {
                    max_dense_bytes: 1,
                    ..limits
                },
                || false
            )
            .is_err()
    );
    assert_eq!(
        prepared
            .log_length_sensitivity(&original, &limits, || false)
            .unwrap(),
        reference
    );
}
