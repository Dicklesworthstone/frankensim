//! Independent manufactured Maxwell/Debye references; no measured-material claim.
#[path = "relaxation/thermal.rs"]
mod thermal;
#[path = "relaxation/transient.rs"]
mod transient;
#[path = "relaxation/equilibrium.rs"]
mod equilibrium;

use super::*;
use fs_material::piezoelectric::tensor::relaxation::{
    DirectedDebyeProperties, RelaxingTensorPiezoelectric,
};
use fs_material::visco_tensor::source::TensorBranchProperties;
use fs_math::c64::C64;
use fs_qty::Time;

const MECHANICAL_TIME_S: f64 = 0.002;
const DIELECTRIC_TIME_S: f64 = 0.01;
const PERMITTIVITY_INCREMENT_F_M: f64 = 3e-9;

fn diagonal_branch() -> C {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i == j {
                [1e10, 2e10, 3e10, 1.5e10, 2.5e10, 4e10][i]
            } else {
                0.0
            }
        })
    })
}

fn fixture(
    branch: C,
    direction: [f64; 3],
    mut alter: impl FnMut(&mut PropertyClaim),
) -> (
    MaterialCard,
    TensorBranchProperties,
    DirectedDebyeProperties,
) {
    let (c, e, p) = coefficients();
    let original = source(c, e, p, |_| {});
    let mechanical = TensorBranchProperties {
        stiffness_upper: core::array::from_fn(|i| format!("piezo_memory_c_{i:02}")),
        relaxation_time: "piezo_memory_tau".into(),
    };
    let dielectric = DirectedDebyeProperties {
        permittivity_increment: "piezo_polarization_increment".into(),
        relaxation_time: "piezo_polarization_tau".into(),
        direction: core::array::from_fn(|i| format!("piezo_polarization_direction_{i}")),
    };
    let mut values = vec![
        (
            mechanical.relaxation_time.clone(),
            Time::DIMS,
            MECHANICAL_TIME_S,
        ),
        (
            dielectric.relaxation_time.clone(),
            Time::DIMS,
            DIELECTRIC_TIME_S,
        ),
        (
            dielectric.permittivity_increment.clone(),
            ELECTRIC_PERMITTIVITY_DIMS,
            PERMITTIVITY_INCREMENT_F_M,
        ),
    ];
    let mut entry = 0;
    for i in 0..6 {
        for j in i..6 {
            values.push((
                mechanical.stiffness_upper[entry].clone(),
                Pressure::DIMS,
                branch[i][j],
            ));
            entry += 1;
        }
    }
    for (name, value) in dielectric.direction.iter().zip(direction) {
        values.push((name.clone(), Dims::NONE, value));
    }
    let mut claims = original
        .claims()
        .claims_ordered()
        .map(|(_, claim)| claim.clone())
        .collect::<Vec<_>>();
    for (name, dims, value) in values {
        claims.push(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored anisotropic piezoelectric relaxation fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        });
    }
    let mut set = ClaimSet::new();
    for mut claim in claims {
        alter(&mut claim);
        set.insert_claim(claim).unwrap();
    }
    (
        MaterialCard::assemble(original.id().clone(), set, vec![]).unwrap(),
        mechanical,
        dielectric,
    )
}

fn base_law(card: &MaterialCard) -> TensorPiezoelectric {
    TensorPiezoelectric::resolve(
        card,
        &point().with("omega", 100.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap()
}

fn bound(
    card: &MaterialCard,
    mechanical: &TensorBranchProperties,
    dielectric: &DirectedDebyeProperties,
) -> RelaxingTensorPiezoelectric {
    RelaxingTensorPiezoelectric::resolve(
        base_law(card),
        card,
        vec![mechanical.clone()],
        vec![dielectric.clone()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn complex_near(actual: C64, expected: C64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:?} != {expected:?}"
    );
}

fn complex_apply<const R: usize, const N: usize>(
    matrix: &[[f64; N]; R],
    vector: &[C64; N],
) -> [C64; R] {
    core::array::from_fn(|i| {
        matrix[i]
            .iter()
            .zip(vector)
            .fold(C64::ZERO, |sum, (&a, &z)| sum + z.scale(a))
    })
}

fn quadratic<const N: usize>(matrix: &[[f64; N]; N], vector: &[C64; N]) -> f64 {
    let product = complex_apply(matrix, vector);
    vector
        .iter()
        .zip(product)
        .map(|(&a, b)| (a.conj() * b).re)
        .sum()
}

#[test]
fn g1_complex_coefficients_match_maxwell_debye_reference_and_frequency_limits() {
    let branch = diagonal_branch();
    let (card, mechanical, dielectric) = fixture(branch, [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let (c, _, epsilon) = coefficients();
    let dc = material.frequency_coefficients(0.0).unwrap();
    assert_eq!(dc.stiffness_increment_mandel_pa, [[C64::ZERO; 6]; 6]);
    assert_eq!(
        dc.permittivity_increment_f_m[2][2],
        C64::from_re(PERMITTIVITY_INCREMENT_F_M)
    );
    for omega in [1e-6, 1e-2, 100.0, 500.0, 1e6, 1e10] {
        let response = material.frequency_coefficients(omega).unwrap();
        let x = omega * MECHANICAL_TIME_S;
        let y = omega * DIELECTRIC_TIME_S;
        // Independent scalar solutions of tau*z_dot + z = strain/field,
        // expressed with the public exp(-i*omega*t) phasor convention.
        let maxwell = C64::new(x * x / (1.0 + x * x), -x / (1.0 + x * x));
        let debye = C64::new(1.0 / (1.0 + y * y), y / (1.0 + y * y));
        for i in 0..6 {
            for j in 0..6 {
                let increment = maxwell.scale(branch[i][j]);
                complex_near(
                    response.stiffness_increment_mandel_pa[i][j],
                    increment,
                    1e-9,
                    3e-13,
                );
                complex_near(
                    response.stiffness_mandel_pa[i][j],
                    C64::from_re(c[i][j]) + increment,
                    1e-5,
                    3e-13,
                );
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                let increment = if i == 2 && j == 2 {
                    debye.scale(PERMITTIVITY_INCREMENT_F_M)
                } else {
                    C64::ZERO
                };
                complex_near(
                    response.permittivity_increment_f_m[i][j],
                    increment,
                    1e-22,
                    3e-13,
                );
                complex_near(
                    response.permittivity_f_m[i][j],
                    C64::from_re(epsilon[i][j]) + increment,
                    1e-22,
                    3e-13,
                );
            }
        }
        assert!(response.stiffness_increment_mandel_pa[0][0].im < 0.0);
        assert!(response.permittivity_increment_f_m[2][2].im > 0.0);
    }
    let low = material.frequency_coefficients(1e-6).unwrap();
    let high = material.frequency_coefficients(1e10).unwrap();
    assert!(low.stiffness_increment_mandel_pa[0][0].abs() < 3e-9 * branch[0][0]);
    near(
        high.stiffness_increment_mandel_pa[0][0].re,
        branch[0][0],
        1e-4,
        1e-12,
    );
    near(
        low.permittivity_increment_f_m[2][2].re,
        PERMITTIVITY_INCREMENT_F_M,
        1e-22,
        1e-12,
    );
    assert!(high.permittivity_increment_f_m[2][2].abs() < 2e-8 * PERMITTIVITY_INCREMENT_F_M);
}

#[test]
fn g1_branch_storage_heat_and_constitutive_mean_work_are_independent_accounts() {
    let branch = diagonal_branch();
    let direction = [0.6, 0.8, 0.0];
    let (card, mechanical, dielectric) = fixture(branch, direction, |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let strain = [
        C64::new(1e-5, -2e-5),
        C64::new(-2e-5, 0.5e-5),
        C64::new(3e-5, 1e-5),
        C64::new(4e-6, -3e-6),
        C64::new(-5e-6, 2e-6),
        C64::new(2e-6, 7e-6),
    ];
    let field = [
        C64::new(1000.0, 300.0),
        C64::new(-2000.0, 900.0),
        C64::new(3000.0, -700.0),
    ];
    let omega = 500.0;
    let response = material.harmonic_response(omega, strain, field).unwrap();
    let x = omega * MECHANICAL_TIME_S;
    let y = omega * DIELECTRIC_TIME_S;
    let projected_field = direction
        .iter()
        .zip(field)
        .fold(C64::ZERO, |sum, (&a, z)| sum + z.scale(a));
    let mechanical_energy = 0.25 * quadratic(&branch, &strain) * x * x / (1.0 + x * x);
    let dielectric_energy =
        0.25 * PERMITTIVITY_INCREMENT_F_M * projected_field.norm_sq() / (1.0 + y * y);
    let mechanical_heat = 2.0 * mechanical_energy / MECHANICAL_TIME_S;
    let dielectric_heat = 2.0 * dielectric_energy * omega * omega * DIELECTRIC_TIME_S;
    for (actual, expected) in [
        (
            response.average_mechanical_memory_energy_j_m3,
            mechanical_energy,
        ),
        (
            response.average_dielectric_memory_energy_j_m3,
            dielectric_energy,
        ),
        (response.average_mechanical_heat_w_m3, mechanical_heat),
        (response.average_dielectric_heat_w_m3, dielectric_heat),
    ] {
        assert!(actual > 0.0);
        near(actual, expected, 1e-12, 5e-13);
    }
    let (c, _, epsilon) = coefficients();
    near(
        response.average_stored_energy_j_m3,
        0.25 * (quadratic(&c, &strain) + quadratic(&epsilon, &field))
            + mechanical_energy
            + dielectric_energy,
        1e-10,
        5e-13,
    );
    let derivative = |z: C64| C64::new(omega * z.im, -omega * z.re);
    let mechanical_work = response
        .stress_mandel_pa
        .iter()
        .zip(strain)
        .map(|(&stress, s)| 0.5 * (stress * derivative(s).conj()).re)
        .sum::<f64>();
    let electrical_work = field
        .iter()
        .zip(response.displacement_c_m2)
        .map(|(&e, d)| 0.5 * (e * derivative(d).conj()).re)
        .sum::<f64>();
    near(
        mechanical_work + electrical_work,
        mechanical_heat + dielectric_heat,
        1e-8,
        2e-12,
    );
    // Resolve each physical branch separately: neither can manufacture the
    // other's heat, and their extra stores add without duplicating the base.
    for (mechanical_only, branches, poles) in [
        (true, vec![mechanical.clone()], vec![]),
        (false, vec![], vec![dielectric.clone()]),
    ] {
        let one = RelaxingTensorPiezoelectric::resolve(
            base_law(&card),
            &card,
            branches,
            poles,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .harmonic_response(omega, strain, field)
        .unwrap();
        near(
            one.average_mechanical_heat_w_m3,
            if mechanical_only {
                mechanical_heat
            } else {
                0.0
            },
            1e-12,
            5e-13,
        );
        near(
            one.average_dielectric_heat_w_m3,
            if mechanical_only {
                0.0
            } else {
                dielectric_heat
            },
            1e-12,
            5e-13,
        );
        near(
            one.average_stored_energy_j_m3,
            response.average_stored_energy_j_m3
                - if mechanical_only {
                    dielectric_energy
                } else {
                    mechanical_energy
                },
            1e-10,
            5e-13,
        );
    }
}

#[test]
fn g1_rank_deficient_memory_and_directed_polarization_preserve_null_spaces() {
    let v = [1.0, -1.0, 0.0, 2.0, 0.0, 0.0];
    let branch = core::array::from_fn(|i| core::array::from_fn(|j| 1e10 * v[i] * v[j]));
    let direction = [0.6, 0.8, 0.0];
    let (card, mechanical, dielectric) = fixture(branch, direction, |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let omega = 500.0;
    let null_strain = [
        C64::new(2e-5, -1e-5),
        C64::new(2e-5, -1e-5),
        C64::new(3e-5, 4e-6),
        C64::ZERO,
        C64::new(-1e-5, 2e-5),
        C64::ZERO,
    ];
    let null_field = [C64::ZERO, C64::ZERO, C64::new(2000.0, -500.0)];
    let null = material
        .harmonic_response(omega, null_strain, null_field)
        .unwrap();
    for actual in null.additional_stress_mandel_pa {
        complex_near(actual, C64::ZERO, 1e-9, 0.0);
    }
    for actual in null.polarization_c_m2 {
        complex_near(actual, C64::ZERO, 1e-22, 0.0);
    }
    for value in [
        null.average_mechanical_memory_energy_j_m3,
        null.average_dielectric_memory_energy_j_m3,
        null.average_mechanical_heat_w_m3,
        null.average_dielectric_heat_w_m3,
    ] {
        near(value, 0.0, 1e-20, 0.0);
    }
    let strain = v.map(|a| C64::new(1e-5, -0.5e-5).scale(a));
    let field = direction.map(|a| C64::new(1000.0, 300.0).scale(a));
    let active = material.harmonic_response(omega, strain, field).unwrap();
    let expected_stress = complex_apply(&branch, &strain);
    let x = omega * MECHANICAL_TIME_S;
    let y = omega * DIELECTRIC_TIME_S;
    let maxwell = C64::new(x * x / (1.0 + x * x), -x / (1.0 + x * x));
    let polarization = C64::new(1000.0, 300.0)
        * C64::new(1.0 / (1.0 + y * y), y / (1.0 + y * y)).scale(PERMITTIVITY_INCREMENT_F_M);
    for (actual, expected) in active
        .additional_stress_mandel_pa
        .into_iter()
        .zip(expected_stress)
    {
        complex_near(actual, expected * maxwell, 1e-8, 5e-13);
    }
    for (actual, component) in active.polarization_c_m2.into_iter().zip(direction) {
        complex_near(actual, polarization.scale(component), 1e-21, 5e-13);
    }
    assert!(active.average_mechanical_heat_w_m3 > 0.0);
    assert!(active.average_dielectric_heat_w_m3 > 0.0);
}

#[test]
fn g3_corotation_preserves_complex_fields_storage_heat_and_original_sources() {
    let v = [1.0, -1.0, 0.0, 2.0, 0.0, 0.0];
    let branch = core::array::from_fn(|i| core::array::from_fn(|j| 1e10 * v[i] * v[j]));
    let (card, mechanical, dielectric) = fixture(branch, [0.6, 0.8, 0.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
    let m = mandel_rotation(q);
    let target = ContentHash([29; 32]);
    let rotated = material.rotated(target, q).unwrap();
    let rebound = RelaxingTensorPiezoelectric::resolve(
        material.base().rotated(target, q).unwrap(),
        &card,
        vec![mechanical],
        vec![dielectric],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let strain = [
        C64::new(1e-5, 2e-5),
        C64::new(-2e-5, 1e-5),
        C64::new(3e-5, -2e-5),
        C64::new(4e-6, 7e-6),
        C64::new(-5e-6, 1e-6),
        C64::new(2e-6, -3e-6),
    ];
    let field = [
        C64::new(1000.0, -500.0),
        C64::new(-2000.0, 700.0),
        C64::new(3000.0, 900.0),
    ];
    let original = material.harmonic_response(170.0, strain, field).unwrap();
    for law in [&rotated, &rebound] {
        let moved = law
            .harmonic_response(170.0, complex_apply(&m, &strain), complex_apply(&q, &field))
            .unwrap();
        for (actual, expected) in moved
            .stress_mandel_pa
            .into_iter()
            .zip(complex_apply(&m, &original.stress_mandel_pa))
        {
            complex_near(actual, expected, 1e-7, 5e-13);
        }
        for (actual, expected) in moved
            .displacement_c_m2
            .into_iter()
            .zip(complex_apply(&q, &original.displacement_c_m2))
        {
            complex_near(actual, expected, 1e-17, 5e-13);
        }
        for (actual, expected) in [
            (
                moved.average_stored_energy_j_m3,
                original.average_stored_energy_j_m3,
            ),
            (
                moved.average_mechanical_memory_energy_j_m3,
                original.average_mechanical_memory_energy_j_m3,
            ),
            (
                moved.average_dielectric_memory_energy_j_m3,
                original.average_dielectric_memory_energy_j_m3,
            ),
            (
                moved.average_mechanical_heat_w_m3,
                original.average_mechanical_heat_w_m3,
            ),
            (
                moved.average_dielectric_heat_w_m3,
                original.average_dielectric_heat_w_m3,
            ),
        ] {
            near(actual, expected, 1e-9, 5e-13);
        }
        assert_eq!(law.material(), material.material());
        assert_eq!(law.base().material(), material.base().material());
        assert_eq!(law.base().source_frame_identity(), FRAME);
        assert_eq!(law.base().frame_identity(), target);
        assert_eq!(law.frequency_band_rad_s(), material.frequency_band_rad_s());
    }
    assert_eq!(material.material().card_identity(), card.content_hash());
}

#[test]
fn g0_source_identity_units_domain_direction_frames_and_frequency_band_are_enforced() {
    let branch = diagonal_branch();
    let (card, mechanical, dielectric) = fixture(branch, [0.0, 0.0, 1.0], |_| {});
    let material = bound(&card, &mechanical, &dielectric);
    let resolve = |card: &MaterialCard,
                   mechanical: &TensorBranchProperties,
                   dielectric: &DirectedDebyeProperties| {
        RelaxingTensorPiezoelectric::resolve(
            base_law(card),
            card,
            vec![mechanical.clone()],
            vec![dielectric.clone()],
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    let (other, _, _) = fixture(branch, [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_memory_tau" {
            claim.value = PropertyValue::Scalar {
                value: 0.003,
                dims: Time::DIMS,
            };
        }
    });
    assert!(
        RelaxingTensorPiezoelectric::resolve(
            base_law(&card),
            &other,
            vec![mechanical.clone()],
            vec![dielectric.clone()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    let (wrong_units, m, d) = fixture(branch, [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_memory_tau" {
            claim.key = PropertyKey::new("piezo_memory_tau", Dims::NONE);
            claim.value = PropertyValue::Scalar {
                value: 0.002,
                dims: Dims::NONE,
            };
        }
    });
    assert!(resolve(&wrong_units, &m, &d).is_err());
    let (missing, m, d) = fixture(branch, [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_polarization_direction_2" {
            claim.key = PropertyKey::new("unrelated_direction", Dims::NONE);
        }
    });
    assert!(resolve(&missing, &m, &d).is_err());
    let (outside_domain, m, d) = fixture(branch, [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_memory_tau" {
            claim.validity = ValidityDomain::unconstrained().with("T", 310.0, 400.0);
        }
    });
    assert!(resolve(&outside_domain, &m, &d).is_err());
    let (sampled_spectrum, m, d) = fixture(branch, [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_polarization_tau" {
            claim.value = PropertyValue::Curve {
                abscissa: "omega".into(),
                abscissa_dims: fs_qty::Frequency::DIMS,
                knots: vec![(10.0, 0.01), (1000.0, 0.02)],
                dims: Time::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            claim.validity = claim.validity.clone().with("omega", 10.0, 1000.0);
        }
    });
    assert!(resolve(&sampled_spectrum, &m, &d).is_err());
    for name in [
        "piezo_memory_tau",
        "piezo_polarization_tau",
        "piezo_polarization_increment",
    ] {
        let (zero_parameter, m, d) = fixture(branch, [0.0, 0.0, 1.0], |claim| {
            if claim.key.name() == name {
                let PropertyValue::Scalar { value, .. } = &mut claim.value else {
                    unreachable!()
                };
                *value = 0.0;
            }
        });
        assert!(resolve(&zero_parameter, &m, &d).is_err());
    }
    for direction in [[0.0; 3], [0.0, 0.0, 2.0]] {
        let (bad_direction, m, d) = fixture(branch, direction, |_| {});
        assert!(resolve(&bad_direction, &m, &d).is_err());
    }
    let mut indefinite = branch;
    indefinite[0][0] = -1e10;
    let (unstable, m, d) = fixture(indefinite, [0.0, 0.0, 1.0], |_| {});
    assert!(resolve(&unstable, &m, &d).is_err());
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(material.rotated(ContentHash([0; 32]), identity).is_err());
    assert!(
        material
            .rotated(FRAME, rotation([0.0, 0.0, 1.0], 0.4))
            .is_err()
    );
    assert!(
        material
            .rotated(
                ContentHash([30; 32]),
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            )
            .is_err()
    );
    let (limited, m, d) = fixture(branch, [0.0, 0.0, 1.0], |claim| {
        if claim.key.name() == "piezo_memory_tau" {
            claim.validity = claim.validity.clone().with("omega", 10.0, 1000.0);
        }
        if claim.key.name() == "piezo_polarization_increment" {
            claim.validity = claim.validity.clone().with("omega", 50.0, 500.0);
        }
    });
    let limited = bound(&limited, &m, &d);
    assert_eq!(limited.frequency_band_rad_s(), [50.0, 500.0]);
    for omega in [50.0, 500.0] {
        assert!(limited.frequency_coefficients(omega).is_ok());
    }
    for omega in [-1.0, 0.0, 49.0, 501.0, f64::NAN, f64::INFINITY] {
        assert!(limited.frequency_coefficients(omega).is_err());
        assert!(
            limited
                .harmonic_response(omega, [C64::ZERO; 6], [C64::ZERO; 3])
                .is_err()
        );
    }
    assert!(
        limited
            .harmonic_response(100.0, [C64::new(f64::NAN, 0.0); 6], [C64::ZERO; 3])
            .is_err()
    );
}
