use super::*;

fn tensor(matrix: [[C64; 3]; 3]) -> ElasticSurfaceCondition {
    ElasticSurfaceCondition::TensorImpedance(ElasticSurfaceImpedanceTensor::new(matrix).unwrap())
}

fn diagonal(values: [C64; 3]) -> [[C64; 3]; 3] {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { C64::ZERO }))
}

#[test]
fn g3_diagonal_tensor_recovers_native_axisymmetric_reflection() {
    let state = material(1.0);
    let normal = C64::new(2e6, 3e6);
    let tangent = C64::new(7e6, -2e6);
    let matrix = diagonal([tangent, tangent, normal]);
    for mode in [ElasticMode::P, ElasticMode::Sv] {
        for cosine in [1.0, 0.9, 0.5] {
            let scalar = IsotropicElasticSurface::new(&state, support(normal, tangent))
                .unwrap()
                .prepare_incidence(mode, cosine)
                .unwrap();
            let coupled = IsotropicElasticSurface::new(&state, tensor(matrix))
                .unwrap()
                .prepare_incidence(mode, cosine)
                .unwrap();
            let input = core::array::from_fn(|i| {
                if scalar.incoming_mode(MODES[i]).is_some() {
                    C64::new(0.001 * (i + 1) as f64, -0.002)
                } else {
                    C64::ZERO
                }
            });
            let a = scalar.drive(input).unwrap();
            let b = coupled.drive(input).unwrap();
            assert_eq!(a.material_state_identity, b.material_state_identity);
            for i in 0..3 {
                close(
                    a.reflected_velocity_m_s[i],
                    b.reflected_velocity_m_s[i],
                    2e-17,
                );
                close(
                    a.boundary_velocity_m_s[i],
                    b.boundary_velocity_m_s[i],
                    2e-17,
                );
                close(a.boundary_traction_pa[i], b.boundary_traction_pa[i], 2e-9);
            }
            assert!(
                (a.absorbed_normal_power_w_m2 - b.absorbed_normal_power_w_m2).abs()
                    < 2e-12 * a.total_incoming_power_w_m2
            );
        }
    }
}

// The two tangent eigenloads are Zs and 3 Zs, rotated by 45 degrees. Their
// scalar reflection eigenvalues are 0 and -1/2, independently of our solve.
#[test]
fn g1_rotated_tangential_support_converts_sh_to_sv() {
    let state = material(1.0);
    let zs = 5e6;
    let matrix = [
        [C64::from_re(2.0 * zs), C64::from_re(-zs), C64::ZERO],
        [C64::from_re(-zs), C64::from_re(2.0 * zs), C64::ZERO],
        [C64::ZERO, C64::ZERO, C64::from_re(1e7)],
    ];
    let prepared = IsotropicElasticSurface::new(&state, tensor(matrix))
        .unwrap()
        .prepare_incidence(ElasticMode::Sh, 1.0)
        .unwrap();
    let coefficients = prepared.velocity_coefficients(ElasticMode::Sh).unwrap();
    close(coefficients[0], C64::ZERO, 1e-14);
    close(coefficients[1], C64::from_re(0.25), 1e-14);
    close(coefficients[2], C64::from_re(-0.25), 1e-14);
    let response = prepared
        .drive([C64::ZERO, C64::ZERO, C64::from_re(0.001)])
        .unwrap();
    assert!(
        (response.absorbed_normal_power_w_m2 / response.total_incoming_power_w_m2 - 0.875).abs()
            < 1e-14
    );
}

// In power-normalized Cartesian coordinates, Z'=I+lambda*(ex ez^T+ez ex^T).
// Its exact reflection is (lambda² I - 2 lambda K)/(4-lambda²) in the xz block.
#[test]
fn g1_normal_tangent_coupling_converts_incident_p_into_shear() {
    let state = material(1.0);
    let lambda = 0.5;
    let zs: f64 = 5e6;
    let zp: f64 = 1e7;
    let coupling = lambda * (zs * zp).sqrt();
    let matrix = [
        [C64::from_re(zs), C64::ZERO, C64::from_re(coupling)],
        [C64::ZERO, C64::from_re(zs), C64::ZERO],
        [C64::from_re(coupling), C64::ZERO, C64::from_re(zp)],
    ];
    let prepared = IsotropicElasticSurface::new(&state, tensor(matrix))
        .unwrap()
        .prepare_incidence(ElasticMode::P, 1.0)
        .unwrap();
    let coefficients = prepared.velocity_coefficients(ElasticMode::P).unwrap();
    let denominator = 4.0 - lambda * lambda;
    close(
        coefficients[0],
        C64::from_re(-lambda * lambda / denominator),
        2e-14,
    );
    close(
        coefficients[1],
        C64::from_re(-2.0 * lambda / denominator * (zp / zs).sqrt()),
        2e-14,
    );
    close(coefficients[2], C64::ZERO, 1e-14);
}

#[test]
fn g0_noncommuting_resistance_reactance_preserve_coherent_work_and_evanescence() {
    let state = material(1.0);
    let resistance = [
        [2e6, 0.3e6, 0.2e6],
        [0.3e6, 1e6, -0.1e6],
        [0.2e6, -0.1e6, 4e6],
    ];
    let reactance = [[4e6, 0.8e6, -2e6], [0.8e6, -3e6, 1e6], [-2e6, 1e6, 2e6]];
    for loss_scale in [0.0, 1e-12, 1.0] {
        let matrix = core::array::from_fn(|i| {
            core::array::from_fn(|j| C64::new(loss_scale * resistance[i][j], reactance[i][j]))
        });
        let law = ElasticSurfaceImpedanceTensor::new(matrix).unwrap();
        assert_eq!(law.matrix_pa_s_m(), &matrix);
        for (mode, cosine) in [(ElasticMode::P, 0.7), (ElasticMode::Sv, 0.5)] {
            let prepared =
                IsotropicElasticSurface::new(&state, ElasticSurfaceCondition::TensorImpedance(law))
                    .unwrap()
                    .prepare_incidence(mode, cosine)
                    .unwrap();
            let input = core::array::from_fn(|i| {
                if prepared.incoming_mode(MODES[i]).is_some() {
                    C64::new(0.001 * (i + 1) as f64, -0.0007 * (3 - i) as f64)
                } else {
                    C64::ZERO
                }
            });
            let response = prepared.drive(input).unwrap();
            let power = response.total_incoming_power_w_m2;
            assert!(
                (power - response.total_reflected_power_w_m2 - response.absorbed_normal_power_w_m2)
                    .abs()
                    < 3e-12 * power
            );
            assert!(
                (response.boundary_normal_flux_w_m2 - response.absorbed_normal_power_w_m2).abs()
                    < 3e-12 * power
            );
            let v = response.boundary_velocity_m_s;
            let mut independent_loss = 0.0;
            for i in 0..3 {
                let zv = (0..3).fold(C64::ZERO, |sum, j| sum + matrix[i][j] * v[j]);
                close(response.boundary_traction_pa[i], -zv, 2e-9);
                independent_loss += 0.5
                    * (0..3)
                        .map(|j| loss_scale * resistance[i][j] * (v[i].conj() * v[j]).re)
                        .sum::<f64>();
            }
            if loss_scale == 0.0 {
                assert_eq!(response.absorbed_normal_power_w_m2, 0.0);
            } else {
                assert!(response.absorbed_normal_power_w_m2 > 0.0);
                assert!(
                    (response.absorbed_normal_power_w_m2 / independent_loss - 1.0).abs() < 2e-13
                );
            }
            if mode == ElasticMode::Sv {
                assert_eq!(response.reflected_normal_power_w_m2[0], 0.0);
                assert!(response.reflected_velocity_m_s[0].abs() > 1e-5);
            }
        }
    }
}

#[test]
fn g0_rank_deficient_resistance_keeps_the_undamped_coherent_direction() {
    let state = material(1.0);
    let r = C64::from_re(5e6);
    let matrix = [[r, r, C64::ZERO], [r, r, C64::ZERO], [C64::ZERO; 3]];
    let prepared = IsotropicElasticSurface::new(&state, tensor(matrix))
        .unwrap()
        .prepare_incidence(ElasticMode::Sv, 1.0)
        .unwrap();
    let response = prepared
        .drive([C64::ZERO, C64::from_re(0.001), C64::from_re(-0.001)])
        .unwrap();
    assert!(response.absorbed_normal_power_w_m2 < 1e-25 * response.total_incoming_power_w_m2);
    close(
        response.reflected_velocity_m_s[1],
        C64::from_re(0.001),
        1e-17,
    );
    close(
        response.reflected_velocity_m_s[2],
        C64::from_re(-0.001),
        1e-17,
    );
    assert!(
        prepared
            .drive([C64::ZERO, C64::from_re(0.001), C64::from_re(0.001)])
            .unwrap()
            .absorbed_normal_power_w_m2
            > 0.0
    );
}

#[test]
fn g1_augmented_tensor_coordinates_retain_stiff_wall_motion() {
    let state = material(1.0);
    let z = 1e30;
    let prepared = IsotropicElasticSurface::new(&state, tensor(diagonal([C64::from_re(z); 3])))
        .unwrap()
        .prepare_incidence(ElasticMode::P, 1.0)
        .unwrap();
    let response = prepared
        .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
        .unwrap();
    let velocity = 0.002 * (1e7 / z) / (1.0 + 1e7 / z);
    let loss = 0.5 * z * velocity * velocity;
    assert!((response.boundary_velocity_m_s[2].re / velocity - 1.0).abs() < 2e-14);
    assert!((response.absorbed_normal_power_w_m2 / loss - 1.0).abs() < 2e-14);
}

#[test]
fn g0_g4_tensor_passivity_source_domain_and_cancelled_preparation_refuse() {
    for mut matrix in [
        diagonal([C64::ONE; 3]),
        diagonal([C64::from_re(1e30), C64::ONE, C64::ONE]),
    ] {
        matrix[1][2] = C64::from_re(2.0);
        matrix[2][1] = C64::from_re(2.0);
        assert!(ElasticSurfaceImpedanceTensor::new(matrix).is_err());
    }
    let mut matrix = diagonal([C64::ONE; 3]);
    matrix[0][1] = C64::new(0.0, 1.0);
    assert!(ElasticSurfaceImpedanceTensor::new(matrix).is_err());
    matrix[0][1] = C64::new(f64::NAN, 0.0);
    assert!(ElasticSurfaceImpedanceTensor::new(matrix).is_err());
    let state = material(1.0);
    let owner =
        IsotropicElasticSurface::new(&state, tensor(diagonal([C64::from_re(5e6); 3]))).unwrap();
    let mut polls = 0;
    assert!(
        owner
            .prepare_incidence_cancellable(ElasticMode::P, 0.8, &mut || {
                polls += 1;
                polls == 10
            })
            .is_err()
    );
    assert_eq!(polls, 10);
    let prepared = owner.prepare_incidence(ElasticMode::P, 0.8).unwrap();
    assert!(
        prepared
            .drive([C64::from_re(1e6), C64::ZERO, C64::ZERO])
            .is_err()
    );
    assert_eq!(
        prepared
            .drive([C64::ZERO; 3])
            .unwrap()
            .absorbed_normal_power_w_m2,
        0.0
    );
}

#[test]
fn g0_original_gram_channels_preserve_null_and_subroundoff_directional_loss() {
    let gram = [[3e15, 0.0, 0.0], [4e15, 0.0, 0.0], [0.0; 3]];
    let null_velocity = [C64::from_re(4.0), C64::from_re(-3.0), C64::ZERO];
    let no_loss = ElasticSurfaceImpedanceTensor::from_resistance_gram(gram, [[0.0; 3]; 3]).unwrap();
    assert_eq!(no_loss.absorbed_power(null_velocity).unwrap(), 0.0);

    // This original weak channel is too small to change any f64 Cartesian
    // coefficient beside the stiff channel, but its local work is resolved.
    let mut weak = gram;
    weak[0][1] = 4e-15;
    weak[1][1] = -3e-15;
    let weak_loss =
        ElasticSurfaceImpedanceTensor::from_resistance_gram(weak, [[0.0; 3]; 3]).unwrap();
    assert_eq!(weak_loss.matrix_pa_s_m(), no_loss.matrix_pa_s_m());
    let expected = 0.5 * (25e-15_f64).powi(2);
    assert!((weak_loss.absorbed_power(null_velocity).unwrap() / expected - 1.0).abs() < 2e-15);

    let mut nonreciprocal = [[0.0; 3]; 3];
    nonreciprocal[0][1] = 1.0;
    assert!(ElasticSurfaceImpedanceTensor::from_resistance_gram(gram, nonreciprocal).is_err());
    assert!(
        ElasticSurfaceImpedanceTensor::from_resistance_gram([[f64::INFINITY; 3]; 3], [[0.0; 3]; 3])
            .is_err()
    );
    assert!(
        ElasticSurfaceImpedanceTensor::from_resistance_gram([[1e300; 3]; 3], [[0.0; 3]; 3])
            .is_err()
    );
    assert!(
        ElasticSurfaceImpedanceTensor::from_resistance_gram([[1e-300; 3]; 3], [[0.0; 3]; 3])
            .is_err()
    );
}
