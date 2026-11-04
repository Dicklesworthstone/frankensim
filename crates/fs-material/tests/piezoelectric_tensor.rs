//! G0/G1/G3 synthetic tensor and source checks; no measured-material claim.
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::tensor::{STRESS_CHARGE_DIMS, TensorPiezoelectric};
use fs_material::piezoelectric::thermal::ThermalLongitudinalPiezoelectric;
use fs_material::piezoelectric::{COMPLIANCE_DIMS, PiezoelectricError, STRAIN_CHARGE_DIMS};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::tensor::{mandel_basis, mandel_rotation, rotation};
use fs_qty::{Density, Dims, Pressure};

const FRAME: ContentHash = ContentHash([7; 32]);
type C = [[f64; 6]; 6];
type E = [[f64; 6]; 3];
type P = [[f64; 3]; 3];

fn coefficients() -> (C, E, P) {
    let mut c = [[0.0; 6]; 6];
    for (i, diagonal) in [12e10, 11e10, 10e10, 5e10, 6e10, 7e10]
        .into_iter()
        .enumerate()
    {
        c[i][i] = diagonal;
    }
    for (i, j, value) in [(0, 1, 3e10), (0, 2, 2e10), (1, 2, 2.5e10), (0, 3, 2e9)] {
        c[i][j] = value;
        c[j][i] = value;
    }
    (
        c,
        [
            [2.0, 3.0, -1.0, 8.0, 0.5, 2.0],
            [1.0, -4.0, 2.0, 1.0, 5.0, -0.5],
            [-6.0, -6.0, 14.0, 0.25, 1.0, 3.0],
        ],
        [
            [2e-8, 1e-9, 2e-10],
            [1e-9, 1.7e-8, -4e-10],
            [2e-10, -4e-10, 1.1e-8],
        ],
    )
}

fn source(c: C, e: E, epsilon: P, mut alter: impl FnMut(&mut PropertyClaim)) -> MaterialCard {
    let mut values: Vec<(String, Dims, f64)> = vec![("density".into(), Density::DIMS, 7500.0)];
    for i in 0..6 {
        for j in 0..6 {
            values.push((
                format!("piezoelectric_stiffness_e_mandel_{}{}", i + 1, j + 1),
                Pressure::DIMS,
                c[i][j],
            ));
        }
    }
    for i in 0..3 {
        for j in 0..6 {
            values.push((
                format!("piezoelectric_stress_charge_mandel_{}{}", i + 1, j + 1),
                STRESS_CHARGE_DIMS,
                e[i][j],
            ));
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            values.push((
                format!("piezoelectric_permittivity_s_{}{}", i + 1, j + 1),
                ELECTRIC_PERMITTIVITY_DIMS,
                epsilon[i][j],
            ));
        }
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in values {
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored anisotropic piezoelectric tensor fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        };
        alter(&mut claim);
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reciprocal solid".into(),
            phase: "solid".into(),
            process: "authored three-dimensional Mandel stress-charge coefficients".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
}
fn law(c: C, e: E, p: P) -> TensorPiezoelectric {
    TensorPiezoelectric::resolve(
        &source(c, e, p, |_| {}),
        &point(),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap()
}
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn apply<const R: usize, const C: usize>(matrix: &[[f64; C]; R], x: &[f64; C]) -> [f64; R] {
    core::array::from_fn(|i| matrix[i].iter().zip(x).map(|(a, b)| a * b).sum())
}

fn strain_charge_source(s: C, d: E, epsilon_t: P) -> MaterialCard {
    source(s, d, epsilon_t, |claim| {
        for (old, new, dims) in [
            (
                "piezoelectric_stiffness_e_mandel_",
                "piezoelectric_compliance_e_mandel_",
                COMPLIANCE_DIMS,
            ),
            (
                "piezoelectric_stress_charge_mandel_",
                "piezoelectric_strain_charge_mandel_",
                STRAIN_CHARGE_DIMS,
            ),
            (
                "piezoelectric_permittivity_s_",
                "piezoelectric_permittivity_t_",
                ELECTRIC_PERMITTIVITY_DIMS,
            ),
        ] {
            if let Some(suffix) = claim.key.name().strip_prefix(old) {
                let name = format!("{new}{suffix}");
                let PropertyValue::Scalar { value, .. } = claim.value else {
                    unreachable!()
                };
                claim.key = PropertyKey::new(name, dims);
                claim.value = PropertyValue::Scalar { value, dims };
                break;
            }
        }
    })
}

fn strain_charge_coefficients() -> (C, E, P) {
    let mut s: C = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i == j {
                [2e-11, 3e-11, 2.5e-11, 4e-11, 5e-11, 6e-11][i]
            } else {
                0.0
            }
        })
    });
    s[0][1] = -0.5e-11;
    s[1][0] = s[0][1];
    let d = [
        [2e-11, -1e-11, 1e-11, 4e-11, 0.0, 1e-10],
        [-3e-11, 1e-10, -4e-11, 0.0, 2e-10, 0.0],
        [-1e-10, -1e-10, 3e-10, 0.0, 0.0, 0.0],
    ];
    let epsilon_t = [
        [3e-8, 1e-9, 0.0],
        [1e-9, 2.5e-8, -2e-10],
        [0.0, -2e-10, 2e-8],
    ];
    (s, d, epsilon_t)
}

#[test]
fn g1_strain_charge_datasheet_form_retains_sources_and_reciprocal_response() {
    let (s, d, epsilon_t) = strain_charge_coefficients();
    let card = strain_charge_source(s, d, epsilon_t);
    let material = TensorPiezoelectric::resolve_strain_charge(
        &card,
        &point(),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    assert_eq!(
        TensorPiezoelectric::strain_charge_requirements()
            .unwrap()
            .len(),
        64
    );
    assert_eq!(material.material().card_identity(), card.content_hash());
    assert_eq!(material.material().properties().len(), 64);
    assert!(
        material
            .material()
            .property("piezoelectric_compliance_e_mandel_12")
            .is_some()
    );
    assert!(
        material
            .material()
            .property("piezoelectric_stiffness_e_mandel_12")
            .is_none()
    );
    let stress = [1e6, -2e6, 0.5e6, 1e5, -2e5, 3e5];
    let field = [2000.0, -1000.0, 4000.0];
    let mechanical = apply(&s, &stress);
    let strain =
        core::array::from_fn(|i| mechanical[i] + (0..3).map(|j| d[j][i] * field[j]).sum::<f64>());
    let direct = apply(&d, &stress);
    let dielectric = apply(&epsilon_t, &field);
    let displacement: [f64; 3] = core::array::from_fn(|i| direct[i] + dielectric[i]);
    let response = material.field_response(strain, field).unwrap();
    for i in 0..6 {
        near(response.stress_mandel_pa[i], stress[i], 5e-8, 1e-13);
    }
    for i in 0..3 {
        near(response.displacement_c_m2[i], displacement[i], 1e-17, 1e-13);
    }
    let charge = material.charge_response(strain, displacement).unwrap();
    for i in 0..3 {
        near(charge.electric_field_v_m[i], field[i], 1e-8, 1e-13);
    }
    let stress_form = law(
        *material.stiffness_mandel_pa(),
        *material.coupling_c_m2(),
        *material.permittivity_f_m(),
    );
    near(
        response.stored_energy_j_m3,
        stress_form
            .field_response(strain, field)
            .unwrap()
            .stored_energy_j_m3,
        1e-11,
        1e-13,
    );
    let energy = 0.5
        * (strain.iter().zip(stress).map(|(s, t)| s * t).sum::<f64>()
            + field
                .iter()
                .zip(displacement)
                .map(|(e, d)| e * d)
                .sum::<f64>());
    near(response.stored_energy_j_m3, energy, 1e-11, 1e-13);
    let mut unstable = d;
    unstable[0][0] = 1e-8;
    // Both diagonal source blocks remain SPD; the COUPLED Schur complement fails.
    assert!(
        TensorPiezoelectric::resolve_strain_charge(
            &strain_charge_source(s, unstable, epsilon_t),
            &point(),
            MaterialPropertySelection::SingleClaimOnly,
            FRAME
        )
        .is_err()
    );
}

#[test]
fn g1_full_tensor_direct_converse_and_charge_energy_gradients() {
    let (c, e, p) = coefficients();
    let material = law(c, e, p);
    let strain = [1.1e-5, -2.2e-5, 3.3e-5, 4e-6, -5e-6, 2e-6];
    let field = [1000.0, -2000.0, 3000.0];
    let response = material.field_response(strain, field).unwrap();
    let elastic = apply(&c, &strain);
    let polarization = apply(&e, &strain);
    let dielectric = apply(&p, &field);
    for i in 0..6 {
        let stress = elastic[i] - (0..3).map(|j| e[j][i] * field[j]).sum::<f64>();
        near(response.stress_mandel_pa[i], stress, 1e-9, 2e-15);
    }
    for i in 0..3 {
        near(
            response.displacement_c_m2[i],
            polarization[i] + dielectric[i],
            1e-18,
            2e-15,
        );
    }
    let reconstructed = 0.5
        * (strain
            .iter()
            .zip(response.stress_mandel_pa)
            .map(|(s, t)| s * t)
            .sum::<f64>()
            + field
                .iter()
                .zip(response.displacement_c_m2)
                .map(|(e, d)| e * d)
                .sum::<f64>());
    near(response.stored_energy_j_m3, reconstructed, 1e-12, 2e-15);
    let charge = material
        .charge_response(strain, response.displacement_c_m2)
        .unwrap();
    for i in 0..3 {
        near(charge.electric_field_v_m[i], field[i], 1e-10, 2e-14);
    }
    for i in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[i] += 1e-8;
        minus[i] -= 1e-8;
        let derivative = (material
            .charge_response(plus, response.displacement_c_m2)
            .unwrap()
            .stored_energy_j_m3
            - material
                .charge_response(minus, response.displacement_c_m2)
                .unwrap()
                .stored_energy_j_m3)
            / 2e-8;
        near(derivative, response.stress_mandel_pa[i], 2e-5, 1e-10);
    }
    for i in 0..3 {
        let mut plus = response.displacement_c_m2;
        let mut minus = plus;
        plus[i] += 1e-7;
        minus[i] -= 1e-7;
        let derivative = (material
            .charge_response(strain, plus)
            .unwrap()
            .stored_energy_j_m3
            - material
                .charge_response(strain, minus)
                .unwrap()
                .stored_energy_j_m3)
            / 2e-7;
        near(derivative, field[i], 2e-6, 1e-10);
    }
}

#[test]
fn g1_open_short_and_signed_longitudinal_limits() {
    let c: C = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1e11 } else { 0.0 }));
    let p: P = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 2e-8 } else { 0.0 }));
    for coupling in [-8.0, 0.0, 8.0] {
        let mut e = [[0.0; 6]; 3];
        e[2][2] = coupling;
        let material = law(c, e, p);
        let longitudinal = material.longitudinal([0.0, 0.0, 1.0]).unwrap();
        near(
            longitudinal.compliance_at_fixed_field_pa_inv(),
            1e-11,
            1e-26,
            1e-15,
        );
        near(
            longitudinal.strain_charge_coefficient_m_v(),
            coupling / 1e11,
            1e-26,
            1e-15,
        );
        near(
            longitudinal.axially_clamped_permittivity_f_m(),
            2e-8,
            1e-23,
            1e-15,
        );
        near(
            longitudinal.open_circuit_modulus_pa(),
            1e11 + coupling * coupling / 2e-8,
            1e-5,
            1e-15,
        );
        let mut strain = [0.0; 6];
        strain[2] = 1e-5;
        let short = material.field_response(strain, [0.0; 3]).unwrap();
        let open = material.charge_response(strain, [0.0; 3]).unwrap();
        near(short.stress_mandel_pa[2], 1e6, 1e-9, 1e-15);
        near(
            open.stress_mandel_pa[2],
            (1e11 + coupling * coupling / 2e-8) * 1e-5,
            1e-8,
            1e-15,
        );
        near(
            material.open_circuit_stiffness_mandel_pa()[2][2],
            1e11 + coupling * coupling / 2e-8,
            1e-5,
            1e-15,
        );
        near(
            open.electric_field_v_m[2],
            -coupling * 1e-5 / 2e-8,
            1e-9,
            1e-15,
        );
        strain[2] = coupling * 1000.0 / 1e11;
        let free = material.field_response(strain, [0.0, 0.0, 1000.0]).unwrap();
        near(free.stress_mandel_pa[2], 0.0, 1e-9, 0.0);
        near(
            free.displacement_c_m2[2],
            (2e-8 + coupling * coupling / 1e11) * 1000.0,
            1e-18,
            1e-15,
        );
    }
}

#[test]
fn g1_longitudinal_reduction_recovers_anisotropic_strain_charge_work() {
    let (s, d, epsilon_t) = strain_charge_coefficients();
    let card = strain_charge_source(s, d, epsilon_t);
    let tensor = TensorPiezoelectric::resolve_strain_charge(
        &card,
        &point(),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap();
    let n = [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0];
    let reduced = tensor.longitudinal(n).unwrap();
    let basis = mandel_basis();
    let r: [f64; 6] = core::array::from_fn(|a| {
        (0..3)
            .map(|i| (0..3).map(|j| basis[a][i][j] * n[i] * n[j]).sum::<f64>())
            .sum()
    });
    let scalar = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let compliance = scalar(&r, &apply(&s, &r));
    let direct = scalar(&n, &apply(&d, &r));
    let dielectric = scalar(&n, &apply(&epsilon_t, &n));
    near(
        reduced.compliance_at_fixed_field_pa_inv(),
        compliance,
        1e-25,
        1e-13,
    );
    near(
        reduced.strain_charge_coefficient_m_v(),
        direct,
        1e-24,
        1e-13,
    );
    near(
        reduced.permittivity_at_fixed_stress_f_m(),
        dielectric,
        1e-22,
        1e-13,
    );
    assert_eq!(reduced.material(), tensor.material());
    assert_eq!(reduced.material().card_identity(), card.content_hash());
    assert_eq!(reduced.material().properties().len(), 64);
    assert!(
        reduced
            .material()
            .property("piezoelectric_compliance_33_at_constant_field")
            .is_none()
    );
    let receipt = reduced.tensor_reduction().unwrap();
    assert_eq!(receipt.tensor_identity, tensor.identity());
    assert_eq!(receipt.frame_identity, FRAME);
    for i in 0..3 {
        near(receipt.direction[i], n[i], 2e-16, 0.0);
    }
    // The original full strain-charge source is the independent reference. It
    // includes every lateral/shear compliance and coupling, not just C33/e33.
    for (axial_stress, axial_field) in [(2e6, 0.0), (0.0, -4000.0), (2e6, -4000.0)] {
        let stress = r.map(|x| x * axial_stress);
        let field = n.map(|x| x * axial_field);
        let mechanical = apply(&s, &stress);
        let strain: [f64; 6] = core::array::from_fn(|a| {
            mechanical[a] + (0..3).map(|i| d[i][a] * field[i]).sum::<f64>()
        });
        let charge = apply(&d, &stress);
        let dielectric = apply(&epsilon_t, &field);
        let displacement: [f64; 3] = core::array::from_fn(|i| charge[i] + dielectric[i]);
        let axial_strain = scalar(&r, &strain);
        let axial_displacement = scalar(&n, &displacement);
        let full = tensor.field_response(strain, field).unwrap();
        let response = reduced.response(axial_strain, axial_displacement).unwrap();
        near(response.stress_pa, axial_stress, 5e-8, 1e-13);
        near(response.electric_field_v_m, axial_field, 2e-9, 1e-13);
        near(
            response.stored_energy_j_m3,
            full.stored_energy_j_m3,
            1e-11,
            1e-13,
        );
        near(
            response.stored_energy_j_m3,
            0.5 * (axial_stress * axial_strain + axial_field * axial_displacement),
            1e-11,
            1e-13,
        );
    }
}

#[test]
fn g3_longitudinal_tensor_reduction_preserves_orientation_and_relaxation() {
    let (c, e, p) = coefficients();
    let tensor = law(c, e, p);
    let axial = tensor.longitudinal([0.0, 0.0, 1.0]).unwrap();
    assert!((axial.compliance_at_fixed_field_pa_inv() * c[2][2] - 1.0).abs() > 0.01);
    assert!((axial.strain_charge_coefficient_m_v() * c[2][2] - e[2][2]).abs() > 0.1);
    let direction = [1.0, 2.0, 2.0];
    let reduced = tensor.longitudinal(direction).unwrap();
    assert_eq!(reduced, tensor.longitudinal([10.0, 20.0, 20.0]).unwrap());
    let reversed = tensor.longitudinal(direction.map(|x| -x)).unwrap();
    near(
        reversed.compliance_at_fixed_field_pa_inv(),
        reduced.compliance_at_fixed_field_pa_inv(),
        0.0,
        1e-14,
    );
    near(
        reversed.permittivity_at_fixed_stress_f_m(),
        reduced.permittivity_at_fixed_stress_f_m(),
        0.0,
        1e-14,
    );
    near(
        reversed.strain_charge_coefficient_m_v(),
        -reduced.strain_charge_coefficient_m_v(),
        0.0,
        1e-14,
    );
    assert_ne!(reversed.identity(), reduced.identity());
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
    let rotated = tensor.rotated(ContentHash([9; 32]), q).unwrap();
    let moved = rotated.longitudinal(apply(&q, &direction)).unwrap();
    near(
        moved.compliance_at_fixed_field_pa_inv(),
        reduced.compliance_at_fixed_field_pa_inv(),
        1e-25,
        1e-13,
    );
    near(
        moved.strain_charge_coefficient_m_v(),
        reduced.strain_charge_coefficient_m_v(),
        1e-24,
        1e-13,
    );
    near(
        moved.permittivity_at_fixed_stress_f_m(),
        reduced.permittivity_at_fixed_stress_f_m(),
        1e-22,
        1e-13,
    );
    assert_eq!(moved.material(), reduced.material());
    assert_eq!(
        moved.tensor_reduction().unwrap().tensor_identity,
        rotated.identity()
    );
    assert_eq!(
        moved.tensor_reduction().unwrap().frame_identity,
        ContentHash([9; 32])
    );
    assert_ne!(moved.identity(), reduced.identity());
}

#[test]
fn g0_longitudinal_tensor_direction_and_thermal_boundaries_are_explicit() {
    let (c, e, p) = coefficients();
    let card = source(c, e, p, |_| {});
    let tensor = law(c, e, p);
    for direction in [[0.0; 3], [f64::NAN, 0.0, 1.0], [1.0, f64::INFINITY, 0.0]] {
        assert!(tensor.longitudinal(direction).is_err());
    }
    let reduced = tensor.longitudinal([1.0, 0.0, 0.0]).unwrap();
    for direction in [
        [f64::MAX, 0.0, 0.0],
        [f64::MIN_POSITIVE / 4.0, 0.0, 0.0],
        [1.0, -0.0, 0.0],
    ] {
        assert_eq!(tensor.longitudinal(direction).unwrap(), reduced);
    }
    assert!(matches!(
        ThermalLongitudinalPiezoelectric::new(&card, &reduced, "T", [250.0, 400.0]),
        Err(PiezoelectricError::Invalid(
            "tensor-reduced piezoelectric material requires a tensor-aware thermal law"
        ))
    ));
}

#[test]
fn g3_simultaneous_rotation_preserves_work_and_all_tensor_ranks() {
    let (c, e, p) = coefficients();
    let material = law(c, e, p);
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
    let m = mandel_rotation(q);
    let target = ContentHash([8; 32]);
    let rotated = material.rotated(target, q).unwrap();
    let strain = [1e-5, -2e-5, 3e-5, 4e-6, -5e-6, 2e-6];
    let field = [1000.0, -2000.0, 3000.0];
    let base = material.field_response(strain, field).unwrap();
    let moved = rotated
        .field_response(apply(&m, &strain), apply(&q, &field))
        .unwrap();
    let stress = apply(&m, &base.stress_mandel_pa);
    let displacement = apply(&q, &base.displacement_c_m2);
    for i in 0..6 {
        near(moved.stress_mandel_pa[i], stress[i], 3e-8, 3e-14);
    }
    for i in 0..3 {
        near(moved.displacement_c_m2[i], displacement[i], 2e-18, 3e-14);
    }
    near(
        moved.stored_energy_j_m3,
        base.stored_energy_j_m3,
        2e-12,
        3e-14,
    );
    let charge = rotated
        .charge_response(apply(&m, &strain), displacement)
        .unwrap();
    for i in 0..3 {
        near(
            charge.electric_field_v_m[i],
            moved.electric_field_v_m[i],
            1e-9,
            3e-14,
        );
    }
    assert_eq!(rotated.material(), material.material());
    assert_eq!(rotated.source_frame_identity(), FRAME);
    assert_eq!(rotated.frame_identity(), target);
    assert_ne!(rotated.identity(), material.identity());
    assert_eq!(
        rotated.identity(),
        material.rotated(target, q).unwrap().identity()
    );
}

#[test]
fn g0_source_units_support_completeness_and_stability_refuse() {
    let (c, e, p) = coefficients();
    let card = source(c, e, p, |_| {});
    let resolve = |card: &MaterialCard| {
        TensorPiezoelectric::resolve(
            card,
            &point(),
            MaterialPropertySelection::SingleClaimOnly,
            FRAME,
        )
    };
    let material = resolve(&card).unwrap();
    assert_eq!(material.material().card_identity(), card.content_hash());
    assert_eq!(material.material().properties().len(), 64);
    assert_eq!(material.density_kg_m3(), 7500.0);
    assert_eq!(TensorPiezoelectric::requirements().unwrap().len(), 64);
    assert!(
        TensorPiezoelectric::resolve(
            &card,
            &point().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            FRAME
        )
        .is_err()
    );
    assert!(
        TensorPiezoelectric::resolve(
            &card,
            &point(),
            MaterialPropertySelection::SingleClaimOnly,
            ContentHash([0; 32])
        )
        .is_err()
    );
    let missing = source(c, e, p, |claim| {
        if claim.key.name() == "piezoelectric_stress_charge_mandel_36" {
            claim.key = PropertyKey::new("unrelated_coupling", STRESS_CHARGE_DIMS);
        }
    });
    assert!(resolve(&missing).is_err());
    let wrong_units = source(c, e, p, |claim| {
        if claim.key.name() == "piezoelectric_permittivity_s_11" {
            claim.key = PropertyKey::new(claim.key.name(), Dims::NONE);
            claim.value = PropertyValue::Scalar {
                value: 2e-8,
                dims: Dims::NONE,
            };
        }
    });
    assert!(resolve(&wrong_units).is_err());
    let dynamic = source(c, e, p, |claim| {
        claim.validity = claim.validity.clone().with("electric_field", -1.0, 1.0);
    });
    assert!(
        TensorPiezoelectric::resolve(
            &dynamic,
            &point().with("electric_field", 0.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            FRAME
        )
        .is_err()
    );
    let mut asymmetric = c;
    asymmetric[1][0] += 1.0;
    assert!(resolve(&source(asymmetric, e, p, |_| {})).is_err());
    let mut unstable = c;
    unstable[0][1] = 2e11;
    unstable[1][0] = 2e11;
    assert!(resolve(&source(unstable, e, p, |_| {})).is_err());
    let mut unstable_epsilon = p;
    unstable_epsilon[0][1] = 1e-7;
    unstable_epsilon[1][0] = 1e-7;
    assert!(resolve(&source(c, e, unstable_epsilon, |_| {})).is_err());
    assert!(material.field_response([f64::NAN; 6], [0.0; 3]).is_err());
    assert!(material.charge_response([0.0; 6], [f64::MAX; 3]).is_err());
    assert!(
        material
            .rotated(
                ContentHash([8; 32]),
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            )
            .is_err()
    );
    assert!(
        material
            .rotated(
                ContentHash([8; 32]),
                [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            )
            .is_err()
    );
    assert!(
        material
            .rotated(FRAME, rotation([0.0, 0.0, 1.0], 0.1))
            .is_err()
    );
}
