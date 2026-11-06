//! Full-source thermal reduction, independently checked in strain-charge form.
#[path = "thermal/stress_charge.rs"]
mod stress_charge;

use super::*;
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn thermal_point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", temperature(), t).unwrap()
}

fn thermal_source(mut alter: impl FnMut(&mut PropertyClaim)) -> MaterialCard {
    let (s, d, p) = strain_charge_coefficients();
    let original = strain_charge_source(s, d, p);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        let PropertyValue::Scalar { value, dims } = claim.value else {
            unreachable!()
        };
        let slope = slope(claim.key.name(), value);
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, value - 20.0 * slope), (340.0, value + 40.0 * slope)],
            dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", temperature(), 280.0, 340.0);
        alter(&mut claim);
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn slope(name: &str, value: f64) -> f64 {
    if name.starts_with("piezoelectric_compliance_e_mandel") {
        0.001 * value
    } else if name.starts_with("piezoelectric_strain_charge_mandel") {
        -0.002 * value
    } else if name.starts_with("piezoelectric_permittivity_t") {
        0.0015 * value
    } else {
        0.0
    }
}

fn tensor(card: &MaterialCard) -> TensorPiezoelectric {
    TensorPiezoelectric::resolve_strain_charge(
        card,
        &thermal_point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )
    .unwrap()
}

#[test]
fn g1_full_tensor_temperature_reduction_preserves_cross_derivatives_and_potentials() {
    let card = thermal_source(|_| {});
    let tensor = tensor(&card);
    let n = [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0];
    let law = ThermalLongitudinalPiezoelectric::from_tensor(
        &card, &tensor, n, "T", [280.0, 340.0],
    )
    .unwrap();
    assert_eq!(law.original(), &tensor.longitudinal(n).unwrap());
    let (s, d, p) = strain_charge_coefficients();
    let basis = mandel_basis();
    let r: [f64; 6] = core::array::from_fn(|a| {
        (0..3)
            .flat_map(|i| (0..3).map(move |j| basis[a][i][j] * n[i] * n[j]))
            .sum()
    });
    let dot = |a: &[f64], b: &[f64]| {
        a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>()
    };
    let base = [
        dot(&r, &apply(&s, &r)),
        dot(&n, &apply(&d, &r)),
        dot(&n, &apply(&p, &n)),
    ];
    let derivatives = [base[0] * 0.001, base[1] * -0.002, base[2] * 0.0015];
    let (strain, charge, volume) = (8e-6, 2.5e-4, 8e-6);
    for t in [280.0, 300.0, 325.0, 340.0] {
        let at = law.at_temperature(t).unwrap();
        law.at_material(at.material()).unwrap();
        assert_eq!(at.material().material().card_identity(), card.content_hash());
        assert_eq!(at.material().material().properties().len(), 64);
        let [s, d, p] = core::array::from_fn(|i| base[i] + (t - 300.0) * derivatives[i]);
        near(at.material().compliance_at_fixed_field_pa_inv(), s, 1e-26, 2e-13);
        near(at.material().strain_charge_coefficient_m_v(), d, 1e-25, 2e-13);
        near(at.material().permittivity_at_fixed_stress_f_m(), p, 1e-22, 2e-13);
        let determinant = s * p - d * d;
        let stress = (p * strain - d * charge) / determinant;
        let field = (s * charge - d * strain) / determinant;
        let derivative = [
            derivatives[0] * stress + derivatives[1] * field,
            derivatives[1] * stress + derivatives[2] * field,
        ];
        let free = 0.5 * volume * (strain * stress + charge * field);
        let entropy = 0.5 * volume * (stress * derivative[0] + field * derivative[1]);
        let heat_capacity = -t * volume
            * (p * derivative[0] * derivative[0]
                - 2.0 * d * derivative[0] * derivative[1]
                + s * derivative[1] * derivative[1])
            / determinant;
        let state = at.state(strain, charge, volume).unwrap();
        near(state.free_energy_j, free, 1e-18, 2e-13);
        near(state.entropy_j_k, entropy, 1e-21, 2e-13);
        near(state.internal_energy_j, free + t * entropy, 1e-18, 2e-13);
        near(state.heat_capacity_j_k, heat_capacity, 1e-22, 2e-13);
        let lower_bound = law
            .minimum_heat_capacity_j_k(strain, charge, volume)
            .unwrap();
        assert!(state.heat_capacity_j_k >= lower_bound);
    }
    // Independent temperature differentiation at fixed actual strain/charge.
    let state = |t| {
        law.at_temperature(t)
            .unwrap()
            .state(strain, charge, volume)
            .unwrap()
    };
    let h = 0.01;
    near(
        state(315.0).entropy_j_k,
        -(state(315.0 + h).free_energy_j - state(315.0 - h).free_energy_j) / (2.0 * h),
        1e-20,
        2e-8,
    );
    near(
        state(315.0).heat_capacity_j_k,
        (state(315.0 + h).internal_energy_j - state(315.0 - h).internal_energy_j) / (2.0 * h),
        1e-20,
        2e-8,
    );
}

#[test]
fn g1_full_tensor_thermal_signed_coupling_can_cross_zero() {
    let card = thermal_source(|claim| {
        if claim.key.name() == "piezoelectric_strain_charge_mandel_33" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                *knots = vec![(280.0, -4e-11), (340.0, 8e-11)];
            }
        }
    });
    let tensor = tensor(&card);
    let law = ThermalLongitudinalPiezoelectric::from_tensor(
        &card, &tensor, [0.0, 0.0, 1.0], "T", [280.0, 340.0],
    )
    .unwrap();
    for t in [280.0, 300.0, 320.0, 340.0] {
        let at = law.at_temperature(t).unwrap();
        near(
            at.material().strain_charge_coefficient_m_v(),
            (t - 300.0) * 2e-12,
            1e-25,
            2e-13,
        );
        law.at_material(at.material()).unwrap();
    }
    let state = |t| {
        law.at_temperature(t)
            .unwrap()
            .state(2e-5, 3e-4, 1e-5)
            .unwrap()
    };
    near(
        state(300.0).entropy_j_k,
        -(state(300.01).free_energy_j - state(299.99).free_energy_j) / 0.02,
        1e-20,
        2e-8,
    );
}

#[test]
fn g3_tensor_thermal_projection_corotates_and_retains_oriented_electrodes() {
    let card = thermal_source(|_| {});
    let tensor = tensor(&card);
    let n = [1.0, 2.0, -2.0];
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
    let rotated = tensor.rotated(ContentHash([9; 32]), q).unwrap();
    let base = ThermalLongitudinalPiezoelectric::from_tensor(
        &card, &tensor, n, "T", [280.0, 340.0],
    )
    .unwrap();
    let moved = ThermalLongitudinalPiezoelectric::from_tensor(
        &card, &rotated, apply(&q, &n), "T", [280.0, 340.0],
    )
    .unwrap();
    let reversed = ThermalLongitudinalPiezoelectric::from_tensor(
        &card, &tensor, n.map(|v| -v), "T", [280.0, 340.0],
    )
    .unwrap();
    for t in [280.0, 300.0, 320.0, 340.0] {
        let a = base.at_temperature(t).unwrap();
        let b = moved.at_temperature(t).unwrap();
        let c = reversed.at_temperature(t).unwrap();
        let x = a.state(1e-5, 2e-4, 1e-5).unwrap();
        for y in [
            b.state(1e-5, 2e-4, 1e-5).unwrap(),
            c.state(1e-5, -2e-4, 1e-5).unwrap(),
        ] {
            near(x.free_energy_j, y.free_energy_j, 1e-18, 2e-12);
            near(x.entropy_j_k, y.entropy_j_k, 1e-21, 2e-12);
            near(x.heat_capacity_j_k, y.heat_capacity_j_k, 1e-22, 2e-12);
        }
        moved.at_material(b.material()).unwrap();
        assert!(base.at_material(b.material()).is_err());
        assert!(base.at_material(c.material()).is_err());
    }
}

#[test]
fn g0_tensor_thermal_admits_all_source_directions_and_refuses_unsupported_paths() {
    let card = thermal_source(|_| {});
    let valid = tensor(&card);
    let bind = |card: &MaterialCard, tensor: &TensorPiezoelectric| {
        ThermalLongitudinalPiezoelectric::from_tensor(
            card, tensor, [0.0, 0.0, 1.0], "T", [280.0, 340.0],
        )
    };
    let thermal_law = bind(&card, &valid).unwrap();
    assert!(thermal_law.at_temperature(279.0).is_err());
    assert!(thermal_law.at_temperature(f64::NAN).is_err());
    let unsupported = thermal_source(|claim| {
        if claim.key.name() == "piezoelectric_compliance_e_mandel_44" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                knots.insert(1, (300.0, 4.1e-11));
            }
        }
    });
    assert!(bind(&unsupported, &tensor(&unsupported)).is_err());
    // The axial projection remains stable, while an unexcited transverse
    // coupling makes the COMPLETE source block unstable at the endpoints.
    let unstable = thermal_source(|claim| {
        if claim.key.name() == "piezoelectric_strain_charge_mandel_14" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                *knots = vec![(280.0, -1.92e-9), (340.0, 3.96e-9)];
            }
        }
    });
    assert!(bind(&unstable, &tensor(&unstable)).is_err());
    // Affine C/e/epsilonS data generally do not define affine s/d/epsilonT.
    // Freezing that distinction would lose the actual thermal derivatives.
    let (c, e, p) = coefficients();
    let stress_card = source(c, e, p, |_| {});
    let stress_tensor = law(c, e, p);
    assert!(
        ThermalLongitudinalPiezoelectric::from_tensor(
            &stress_card, &stress_tensor, [0.0, 0.0, 1.0], "T", [280.0, 340.0],
        )
        .is_err()
    );
}
