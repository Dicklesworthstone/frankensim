//! G1 exact loaded mechanics/descriptor jets; G3 boundary normalization and
//! source-phase invariance; G0 refusal and immutable-factor replay.
use fs_math::c64::C64;
use fs_phs::harmonic::scattering::{
    ScatteringBudget, ScatteringComplexJet, ScatteringFrequencyParameterResponse, ScatteringModel,
    ScatteringPowerJet,
};
use fs_phs::harmonic::{
    HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection as Direction,
};
use fs_phs::{DescriptorPortHamiltonian, PortHamiltonian, QuadraticStorage};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 16,
        max_dense_bytes: 4_000_000,
        max_dense_work: 8_000_000,
        relative_residual_tolerance: 2e-11,
    }
}
fn wave_budget() -> ScatteringBudget {
    ScatteringBudget {
        maximum_ports: 4,
        maximum_payload_bytes: 1_000_000,
        maximum_work: 2_000_000,
    }
}
fn loaded(
    m: f64,
    k: f64,
    c: f64,
    g: f64,
    z: f64,
    algebraic: Option<f64>,
    wave_budget: ScatteringBudget,
) -> ScatteringModel {
    let storage = Box::new(QuadraticStorage::new(vec![k, 0., 0., 1. / m], 2).unwrap());
    let physical = if let Some(r) = algebraic {
        // q'=v; p'=-kq-cv+lambda+g*u; 0=-v-r*lambda.
        // The algebraic reaction supplies physical loss |v|²/(2r).
        let source = DescriptorPortHamiltonian::new(
            3,
            2,
            1,
            vec![0., 1., 0., -1., 0., 1., 0., -1., 0.],
            vec![0., 0., 0., 0., c, 0., 0., 0., r],
            vec![0., g, 0.],
            storage,
        )
        .unwrap();
        HarmonicModel::from_descriptor(&source, budget(), || false).unwrap()
    } else {
        let source = PortHamiltonian::new(
            2,
            1,
            vec![0., 1., -1., 0.],
            vec![0., 0., 0., c],
            vec![0., g],
            storage,
        )
        .unwrap();
        HarmonicModel::from_ode(&source, budget(), || false).unwrap()
    };
    ScatteringModel::new(physical, &[(0, z)], wave_budget, || false).unwrap()
}
fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * expected.abs().max(1e-4),
        "{actual:?} != {expected:?}"
    );
}
fn real(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn values(jet: ScatteringComplexJet) -> [C64; 5] {
    [
        jet.value,
        jet.frequency,
        jet.second_frequency,
        jet.parameter,
        jet.mixed,
    ]
}
fn powers(jet: ScatteringPowerJet) -> [f64; 5] {
    [
        jet.value,
        jet.frequency,
        jet.second_frequency,
        jet.parameter,
        jet.mixed,
    ]
}
fn compare(actual: [C64; 5], expected: [C64; 5]) {
    for (a, b) in actual.into_iter().zip(expected) {
        near(a, b);
    }
}
fn scale_parameter(a: [C64; 5], value: f64, derivative: f64) -> [C64; 5] {
    [
        a[0].scale(value),
        a[1].scale(value),
        a[2].scale(value),
        a[3].scale(value) + a[0].scale(derivative),
        a[4].scale(value) + a[1].scale(derivative),
    ]
}
fn norm(a: [C64; 5]) -> [f64; 5] {
    [
        a[0].norm_sq(),
        2. * (a[0].conj() * a[1]).re,
        2. * (a[1].norm_sq() + (a[0].conj() * a[2]).re),
        2. * (a[0].conj() * a[3]).re,
        2. * ((a[1].conj() * a[3]).re + (a[0].conj() * a[4]).re),
    ]
}
fn weighted_norm(a: [C64; 5], coefficient: f64, derivative: f64) -> [f64; 5] {
    let n = norm(a);
    [
        coefficient * n[0],
        coefficient * n[1],
        coefficient * n[2],
        coefficient * n[3] + derivative * n[0],
        coefficient * n[4] + derivative * n[1],
    ]
}
fn source(omega: f64) -> [C64; 5] {
    let (delay, delay_theta) = (0.23, 0.07);
    let a = C64::new(0.7, -0.2) * C64::new((omega * delay).cos(), (omega * delay).sin());
    let parameter = a * (C64::new(0.13, -0.04) + C64::new(0., omega * delay_theta));
    [
        a,
        a * C64::new(0., delay),
        a.scale(-delay * delay),
        parameter,
        parameter * C64::new(0., delay) + a * C64::new(0., delay_theta),
    ]
}
fn physical_fields(
    v: &ScatteringFrequencyParameterResponse,
    index: usize,
    effort: bool,
) -> [C64; 5] {
    let p = &v.physical;
    let c = &p.frequency_response;
    let f = &c.frequency_response;
    let t = &p.parameter_response;
    if effort {
        [
            f.response.effort[index],
            f.effort_derivative[index],
            c.effort_second_derivative[index],
            t.effort_derivative[index],
            p.effort_mixed_derivative[index],
        ]
    } else {
        [
            f.response.state[index],
            f.state_derivative[index],
            c.state_second_derivative[index],
            t.state_derivative[index],
            p.state_mixed_derivative[index],
        ]
    }
}

#[test]
fn exact_full_source_mass_loss_coupling_reference_and_algebraic_reaction_jets() {
    let (m, k, g, z, dm, dk, dc, dg, dz, dr) =
        (1.7, 6.8, 0.85, 0.7, 0.2, -0.4, 0.03, 0.1, 0.08, 0.11);
    for algebraic in [None, Some(1.4_f64)] {
        for c in [0.0, 0.12] {
            let owner = loaded(m, k, c, g, z, algebraic, wave_budget());
            let n = if algebraic.is_some() { 3 } else { 2 };
            let mut dq = vec![0.; n * n];
            dq[0] = dk;
            dq[n + 1] = -dm / (m * m);
            let mut dloss = vec![0.; n * n];
            dloss[n + 1] = dc;
            if algebraic.is_some() {
                dloss[2 * n + 2] = dr;
            }
            let mut dport = vec![0.; n];
            dport[1] = dg;
            let direction = Direction {
                storage: Some(&dq),
                dissipation: Some(&dloss),
                ports: Some(&dport),
                ..Default::default()
            };
            for omega in [0.25, 2.0, 2.2, 4.0] {
                let a = source(omega);
                let response = owner
                    .factor(omega, || false)
                    .unwrap()
                    .response_with_frequency_parameter_derivative(
                        &[a[0]],
                        direction,
                        Some(&[a[1]]),
                        Some(&[a[2]]),
                        Some(&[a[3]]),
                        Some(&[a[4]]),
                        Some(&[dz]),
                        || false,
                    )
                    .unwrap();
                let effective_loss = c + g * g / z + algebraic.map_or(0., |r| 1. / r);
                let loss_theta = dc + 2. * g * dg / z
                    - g * g * dz / (z * z)
                    - algebraic.map_or(0., |r| dr / (r * r));
                let d = C64::new(k - m * omega * omega, -omega * effective_loss);
                let dw = C64::new(-2. * m * omega, -effective_loss);
                let dww = C64::from_re(-2. * m);
                let dt = C64::new(dk - dm * omega * omega, -omega * loss_theta);
                let dwt = C64::new(-2. * dm * omega, -loss_theta);
                let force = scale_parameter(a, 2. * g / z, 2. * dg / z - 2. * g * dz / (z * z));
                let q = force[0] / d;
                let qw = (force[1] - dw * q) / d;
                let qt = (force[3] - dt * q) / d;
                let qww = (force[2] - dww * q - (dw * qw).scale(2.)) / d;
                let qwt = (force[4] - dwt * q - dw * qt - dt * qw) / d;
                let q = [q, qw, qww, qt, qwt];
                let s = C64::new(0., -omega);
                let sw = C64::new(0., -1.);
                let v = [
                    s * q[0],
                    s * q[1] + sw * q[0],
                    s * q[2] + (sw * q[1]).scale(2.),
                    s * q[3],
                    s * q[4] + sw * q[3],
                ];
                compare(physical_fields(&response, 0, false), q);
                compare(
                    physical_fields(&response, 1, false),
                    scale_parameter(v, m, dm),
                );
                compare(
                    physical_fields(&response, 0, true),
                    scale_parameter(q, k, dk),
                );
                compare(physical_fields(&response, 1, true), v);
                if let Some(r) = algebraic {
                    let lambda = scale_parameter(v, -1. / r, dr / (r * r));
                    compare(physical_fields(&response, 2, false), lambda);
                    compare(physical_fields(&response, 2, true), lambda);
                }
                let y = scale_parameter(v, g, dg);
                let p = &response.physical;
                let curvature = &p.frequency_response;
                let frequency = &curvature.frequency_response;
                compare(
                    [
                        frequency.response.output[0],
                        frequency.output_derivative[0],
                        curvature.output_second_derivative[0],
                        p.parameter_response.output_derivative[0],
                        p.output_mixed_derivative[0],
                    ],
                    y,
                );
                let outgoing = std::array::from_fn(|j| y[j] - a[j]);
                let numerator = std::array::from_fn(|j| a[j].scale(2.) - y[j]);
                compare(
                    values(response.physical_inputs[0]),
                    scale_parameter(numerator, 1. / z, -dz / (z * z)),
                );
                compare(values(response.waves[0].incident), a);
                compare(values(response.waves[0].outgoing), outgoing);
                let incident_power = weighted_norm(a, 0.5 / z, -0.5 * dz / (z * z));
                let outgoing_power = weighted_norm(outgoing, 0.5 / z, -0.5 * dz / (z * z));
                let potential = weighted_norm(q, 0.25 * k, 0.25 * dk);
                let kinetic = weighted_norm(v, 0.25 * m, 0.25 * dm);
                let heat = weighted_norm(
                    v,
                    0.5 * (c + algebraic.map_or(0., |r| 1. / r)),
                    0.5 * (dc - algebraic.map_or(0., |r| dr / (r * r))),
                );
                let energy = [
                    frequency.response.mean_stored_energy_j,
                    frequency.mean_stored_energy_derivative,
                    curvature.mean_stored_energy_second_derivative,
                    p.parameter_response.mean_stored_energy_derivative,
                    p.mean_stored_energy_mixed_derivative,
                ];
                let physical_heat = [
                    frequency.response.dissipated_power_w,
                    frequency.dissipated_power_derivative,
                    curvature.dissipated_power_second_derivative,
                    p.parameter_response.dissipated_power_derivative,
                    p.dissipated_power_mixed_derivative,
                ];
                let physical_balance = [
                    frequency.response.power_residual_w,
                    frequency.power_residual_derivative,
                    curvature.power_residual_second_derivative,
                    p.parameter_response.power_residual_derivative,
                    p.power_residual_mixed_derivative,
                ];
                for order in 0..5 {
                    real(
                        powers(response.waves[0].incident_power)[order],
                        incident_power[order],
                    );
                    real(
                        powers(response.waves[0].outgoing_power)[order],
                        outgoing_power[order],
                    );
                    real(energy[order], potential[order] + kinetic[order]);
                    real(physical_heat[order], heat[order]);
                    real(powers(response.waves[0].physical_power)[order], heat[order]);
                    real(powers(response.waves[0].power_residual)[order], 0.);
                    real(powers(response.wave_power_residual)[order], 0.);
                    real(physical_balance[order], 0.);
                }
            }
        }
    }
}

#[test]
fn exact_zero_reflection_has_finite_curvature_and_mixed_wave_jets() {
    let owner = loaded(1., 4., 0.5, 1., 2., None, wave_budget());
    let factor = owner.factor(2., || false).unwrap();
    let damping = [0., 0., 0., 0.1];
    let run = |delay: f64| {
        let a = C64::new((2. * delay).cos(), (2. * delay).sin());
        factor
            .response_with_frequency_parameter_derivative(
                &[a],
                Direction {
                    dissipation: Some(&damping),
                    ..Default::default()
                },
                Some(&[a * C64::new(0., delay)]),
                Some(&[a.scale(-delay * delay)]),
                None,
                None,
                None,
                || false,
            )
            .unwrap()
    };
    let fixed = run(0.);
    let delayed = run(0.3);
    assert!(fixed.waves[0].outgoing.value.abs() < 1e-14);
    assert!(fixed.waves[0].outgoing.frequency.abs() > 0.1);
    assert!(fixed.waves[0].outgoing.parameter.abs() > 0.01);
    for (a, b) in powers(fixed.waves[0].outgoing_power)
        .into_iter()
        .zip(powers(delayed.waves[0].outgoing_power))
    {
        real(a, b);
    }
    assert_eq!(fixed, run(0.));
}

#[test]
fn full_jet_refusal_and_cancellation_leave_closed_factor_reusable() {
    let owner = loaded(1., 4., 0.2, 1., 0.4, None, wave_budget());
    let factor = owner.factor(2., || false).unwrap();
    let input = [C64::new(0.3, -0.2)];
    let expected = factor.response(&input, || false).unwrap();
    assert!(
        factor
            .response_with_frequency_parameter_derivative(
                &input,
                Direction::default(),
                Some(&[]),
                None,
                None,
                None,
                None,
                || false,
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        factor.response_with_frequency_parameter_derivative(
            &input,
            Direction::default(),
            None,
            None,
            None,
            None,
            None,
            || {
                polls += 1;
                polls > 20
            },
        ),
        Err(HarmonicError::Cancelled)
    ));
    assert_eq!(expected, factor.response(&input, || false).unwrap());
    let mut small = wave_budget();
    small.maximum_payload_bytes = 3000;
    let limited = loaded(1., 4., 0.2, 1., 0.4, None, small);
    let factor = limited.factor(2., || false).unwrap();
    assert!(
        factor
            .response_with_frequency_parameter_derivative(
                &input,
                Direction::default(),
                None,
                None,
                None,
                None,
                None,
                || false,
            )
            .is_err()
    );
    assert!(factor.response(&input, || false).is_ok());
}

#[test]
fn maxwell_memory_and_ordinary_source_work_retain_curvature_and_mixed_jets() {
    let point = |theta: f64| {
        let source = PortHamiltonian::new(
            2,
            2,
            vec![0., 1., -1., 0.],
            vec![0., 0., 0., 0.2],
            vec![0., 0., 0.6 + 0.02 * theta, 0.9 - 0.05 * theta],
            Box::new(QuadraticStorage::new(vec![110., 0., 0., 2.5], 2).unwrap()),
        )
        .unwrap()
        .with_relaxation_branches(vec![fs_phs::RelaxationBranch {
            projection: vec![1., 0.],
            stiffness: 300. + 20. * theta,
            relaxation_time_s: 0.025 + 0.003 * theta,
        }])
        .unwrap();
        ScatteringModel::new(
            HarmonicModel::from_ode(&source, budget(), || false).unwrap(),
            &[(1, 0.3 + 0.04 * theta)],
            wave_budget(),
            || false,
        )
        .unwrap()
    };
    let drives = |omega: f64, theta: f64| {
        let base = [C64::new(0.1, 0.2), C64::new(0.7, -0.3)];
        let delay = [0.17, 0.08];
        let dt = [0.002, -0.004];
        let gain = [0.11, 0.05];
        let a: [C64; 2] = std::array::from_fn(|p| {
            let angle = omega * (delay[p] + theta * dt[p]);
            base[p] * C64::new(angle.cos(), angle.sin()).scale(1. + theta * gain[p])
        });
        let aw: [C64; 2] = std::array::from_fn(|p| a[p] * C64::new(0., delay[p] + theta * dt[p]));
        (a, aw)
    };
    let omega = 30.;
    let (a, aw) = drives(omega, 0.);
    let delays = [0.17, 0.08];
    let ddelay = [0.002, -0.004];
    let gains = [0.11, 0.05];
    let aww: [C64; 2] = std::array::from_fn(|p| a[p].scale(-delays[p] * delays[p]));
    let at: [C64; 2] = std::array::from_fn(|p| a[p] * C64::new(gains[p], omega * ddelay[p]));
    let awt: [C64; 2] =
        std::array::from_fn(|p| at[p] * C64::new(0., delays[p]) + a[p] * C64::new(0., ddelay[p]));
    let cross = -10. / 300.0_f64.sqrt();
    let dq = [20., 0., cross, 0., 0., 0., cross, 0., 0.];
    let dr = [0., 0., 0., 0., 0., 0., 0., 0., -0.003 / (0.025 * 0.025)];
    let dg = [0., 0., 0.02, -0.05, 0., 0.];
    let owner = point(0.);
    let result = owner
        .factor(omega, || false)
        .unwrap()
        .response_with_frequency_parameter_derivative(
            &a,
            Direction {
                storage: Some(&dq),
                dissipation: Some(&dr),
                ports: Some(&dg),
                ..Default::default()
            },
            Some(&aw),
            Some(&aww),
            Some(&at),
            Some(&awt),
            Some(&[0.04]),
            || false,
        )
        .unwrap();
    assert_eq!(result.waves[0].original_port, 1);
    assert_eq!(
        values(result.physical_inputs[0]),
        [a[0], aw[0], aww[0], at[0], awt[0]]
    );
    let sample = |w: f64, theta: f64| {
        let (input, first) = drives(w, theta);
        point(theta)
            .factor(w, || false)
            .unwrap()
            .response_with_frequency_derivative(&input, Some(&first), || false)
            .unwrap()
    };
    let ht = 1e-5;
    let hw = 1e-4;
    let minus_parameter = sample(omega, -ht);
    let plus_parameter = sample(omega, ht);
    let minus_frequency = sample(omega - hw, 0.);
    let plus_frequency = sample(omega + hw, 0.);
    let compare_complex = |actual: C64, lo: C64, hi: C64, h: f64| {
        let expected = (hi - lo).scale(0.5 / h);
        let rounding = 128. * f64::EPSILON * (lo.abs() + hi.abs()) / h;
        assert!(
            (actual - expected).abs() <= 3e-6 * expected.abs() + rounding + 1e-12,
            "{actual:?} != {expected:?}"
        );
    };
    let compare_real = |actual: f64, lo: f64, hi: f64, h: f64| {
        let expected = (hi - lo) * 0.5 / h;
        let rounding = 128. * f64::EPSILON * (lo.abs() + hi.abs()) / h;
        assert!(
            (actual - expected).abs() <= 3e-6 * expected.abs() + rounding + 1e-12,
            "{actual:e} != {expected:e}"
        );
    };
    let physical = &result.physical;
    let curvature = &physical.frequency_response;
    for (actual, lo, hi, h) in [
        (
            &physical.state_mixed_derivative,
            &minus_parameter.state_derivative,
            &plus_parameter.state_derivative,
            ht,
        ),
        (
            &physical.effort_mixed_derivative,
            &minus_parameter.effort_derivative,
            &plus_parameter.effort_derivative,
            ht,
        ),
        (
            &physical.output_mixed_derivative,
            &minus_parameter.output_derivative,
            &plus_parameter.output_derivative,
            ht,
        ),
        (
            &curvature.state_second_derivative,
            &minus_frequency.state_derivative,
            &plus_frequency.state_derivative,
            hw,
        ),
        (
            &curvature.effort_second_derivative,
            &minus_frequency.effort_derivative,
            &plus_frequency.effort_derivative,
            hw,
        ),
        (
            &curvature.output_second_derivative,
            &minus_frequency.output_derivative,
            &plus_frequency.output_derivative,
            hw,
        ),
    ] {
        for ((&actual, &lo), &hi) in actual.iter().zip(lo).zip(hi) {
            compare_complex(actual, lo, hi, h);
        }
    }
    for p in 0..2 {
        compare_complex(
            result.physical_inputs[p].mixed,
            minus_parameter.physical_input_derivative[p],
            plus_parameter.physical_input_derivative[p],
            ht,
        );
    }
    compare_complex(
        result.waves[0].outgoing.mixed,
        minus_parameter.wave_derivatives[0].outgoing,
        plus_parameter.wave_derivatives[0].outgoing,
        ht,
    );
    for (actual, lo, hi, h) in [
        (
            physical.mean_stored_energy_mixed_derivative,
            minus_parameter.mean_stored_energy_derivative,
            plus_parameter.mean_stored_energy_derivative,
            ht,
        ),
        (
            physical.dissipated_power_mixed_derivative,
            minus_parameter.dissipated_power_derivative,
            plus_parameter.dissipated_power_derivative,
            ht,
        ),
        (
            physical.supplied_power_mixed_derivative,
            minus_parameter.supplied_power_derivative,
            plus_parameter.supplied_power_derivative,
            ht,
        ),
        (
            curvature.mean_stored_energy_second_derivative,
            minus_frequency.mean_stored_energy_derivative,
            plus_frequency.mean_stored_energy_derivative,
            hw,
        ),
        (
            curvature.dissipated_power_second_derivative,
            minus_frequency.dissipated_power_derivative,
            plus_frequency.dissipated_power_derivative,
            hw,
        ),
        (
            result.waves[0].outgoing_power.mixed,
            minus_parameter.wave_derivatives[0].outgoing_power_w,
            plus_parameter.wave_derivatives[0].outgoing_power_w,
            ht,
        ),
    ] {
        compare_real(actual, lo, hi, h);
    }
    for residual in powers(result.wave_power_residual) {
        real(residual, 0.);
    }
}
