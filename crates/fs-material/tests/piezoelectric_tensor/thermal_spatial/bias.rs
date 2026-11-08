//! Authored anisotropic thermo-piezoelectric data, not material qualification.
use super::*;

const ALPHA: [f64; 6] = [1e-5, -2e-6, 8e-6, 3e-6, -1e-6, 2e-6];
const PYRO: [f64; 3] = [2e-4, -1e-4, 3e-4];

fn bias_card(
    base: (C, E, P),
    slopes: (C, E, P),
    alpha: [f64; 6],
    pyro: [f64; 3],
    mut alter: impl FnMut(&mut PropertyClaim) -> bool,
) -> MaterialCard {
    let original = thermal_card(base, slopes, |_| {});
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    for (i, requirement) in ThermalTensorPiezoelectric::thermal_bias_requirements()
        .unwrap()
        .into_iter()
        .enumerate()
    {
        let dims = requirement.dims();
        let mut claim = PropertyClaim {
            key: PropertyKey::new(requirement.name(), dims),
            value: PropertyValue::Scalar {
                value: if i < 6 { alpha[i] } else { pyro[i - 6] },
                dims,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                temperature(),
                280.0,
                340.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "authored Mandel expansion and clamped pyroelectric coefficients".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        };
        if alter(&mut claim) {
            claims.insert_claim(claim).unwrap();
        }
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn constant_card() -> MaterialCard {
    bias_card(
        coefficients(),
        ([[0.0; 6]; 6], [[0.0; 6]; 3], [[0.0; 3]; 3]),
        ALPHA,
        PYRO,
        |_| true,
    )
}

fn biased(card: &MaterialCard) -> ThermalTensorPiezoelectric {
    bind(card)
        .with_thermal_bias(300.0, MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}

#[test]
fn g1_free_clamped_and_open_responses_distinguish_primary_and_secondary_pyroelectricity() {
    let card = constant_card();
    let law = biased(&card);
    let at = law.at_temperature(320.0).unwrap();
    let alpha_dt = ALPHA.map(|a| 20.0 * a);
    let pyro_dt = PYRO.map(|p| 20.0 * p);
    assert_eq!(at.thermal_strain_mandel(), alpha_dt);
    assert_eq!(at.clamped_polarization_c_m2(), pyro_dt);
    let (c, e, epsilon) = coefficients();
    let free = at.field_response(alpha_dt, [0.0; 3]).unwrap();
    assert_eq!(free.stress_mandel_pa, [0.0; 6]);
    assert_eq!(free.stored_energy_j_m3, 0.0);
    let secondary = apply(&e, &alpha_dt);
    for i in 0..3 {
        near(
            free.displacement_c_m2[i],
            pyro_dt[i] + secondary[i],
            1e-18,
            1e-13,
        );
    }
    let clamped = at.field_response([0.0; 6], [0.0; 3]).unwrap();
    assert_eq!(clamped.displacement_c_m2, pyro_dt);
    let stress = apply(&c, &alpha_dt).map(|value| -value);
    for i in 0..6 {
        near(clamped.stress_mandel_pa[i], stress[i], 1e-8, 1e-13);
        near(at.stress_offset_mandel_pa()[i], stress[i], 1e-8, 1e-13);
    }
    near(
        clamped.stored_energy_j_m3,
        0.5 * dot(&alpha_dt, &apply(&c, &alpha_dt)),
        1e-10,
        1e-13,
    );
    let open = at.charge_response([0.0; 6], [0.0; 3]).unwrap();
    let field = solve_dielectric(epsilon, pyro_dt.map(|p| -p));
    for i in 0..3 {
        near(open.electric_field_v_m[i], field[i], 1e-8, 1e-13);
    }
    for i in 0..6 {
        let expected = stress[i] - (0..3).map(|j| e[j][i] * field[j]).sum::<f64>();
        near(open.stress_mandel_pa[i], expected, 1e-8, 1e-13);
    }
    assert_eq!(open.displacement_c_m2, [0.0; 3]);
    let source = at.thermal_bias_source().unwrap();
    assert_eq!(source.properties().len(), 9);
    assert_eq!(source.card_identity(), card.content_hash());
    assert_eq!(source.query(), point_at(320.0));
    assert_eq!(at.thermal_bias_reference_temperature_k(), Some(300.0));
    assert_eq!(at.material().material().properties().len(), 64);
    law.at_material(at.material()).unwrap();
    let different_reference = bind(&card)
        .with_thermal_bias(310.0, MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .at_temperature(320.0)
        .unwrap();
    assert_ne!(
        at.thermal_bias_identity(),
        different_reference.thermal_bias_identity()
    );
    let reference = law.at_temperature(300.0).unwrap();
    let strain = [2e-5; 6];
    let field = [1e3, -2e3, 3e3];
    assert_eq!(
        reference.field_response(strain, field).unwrap(),
        reference.material().field_response(strain, field).unwrap()
    );
    let unbiased = bind(&card).at_temperature(320.0).unwrap();
    assert!(unbiased.thermal_bias_source().is_none());
    assert_eq!(unbiased.thermal_strain_mandel(), [0.0; 6]);
    assert_eq!(
        unbiased.field_response(strain, field).unwrap(),
        unbiased.material().field_response(strain, field).unwrap()
    );
}

fn independent_free(t: f64, s: [f64; 6], d: [f64; 3], volume: f64) -> f64 {
    let base = coefficients();
    let slopes = dense_slopes(base.0, base.1);
    let c = at(base.0, slopes.0, t);
    let e = at(base.1, slopes.1, t);
    let epsilon = at(base.2, slopes.2, t);
    let x = core::array::from_fn(|i| s[i] - ALPHA[i] * (t - 300.0));
    let piezo = apply(&e, &s);
    let q = core::array::from_fn(|i| d[i] - piezo[i] - PYRO[i] * (t - 300.0));
    0.5 * volume * (dot(&x, &apply(&c, &x)) + dot(&q, &solve_dielectric(epsilon, q)))
}

#[test]
fn g1_g3_biased_potential_gradients_entropy_capacity_and_corotation_share_one_law() {
    let base = coefficients();
    let card = bias_card(base, dense_slopes(base.0, base.1), ALPHA, PYRO, |_| true);
    let law = biased(&card);
    let q = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.73);
    let rotated = resolve_tensor(&card)
        .rotated(ContentHash([9; 32]), q)
        .unwrap();
    let moved = ThermalTensorPiezoelectric::new(&card, &rotated, "T", [280.0, 340.0])
        .unwrap()
        .with_thermal_bias(300.0, MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let s = [2e-4, -1e-4, 3e-4, 8e-5, -3e-5, 2e-5];
    let d = [3e-3, -1e-3, 2e-3];
    let volume = 1e-5;
    let lower = law
        .minimum_heat_capacity_j_k(norm(&s), norm(&d), volume)
        .unwrap();
    near(
        lower,
        moved
            .minimum_heat_capacity_j_k(norm(&s), norm(&d), volume)
            .unwrap(),
        1e-15,
        3e-12,
    );
    for t in [280.0, 300.0, 315.0, 340.0] {
        let current = law.at_temperature(t).unwrap();
        let mapped = moved.at_temperature(t).unwrap();
        let state = current.state(s, d, volume).unwrap();
        let rotated_state = mapped
            .state(apply(&mandel_rotation(q), &s), apply(&q, &d), volume)
            .unwrap();
        near(
            state.free_energy_j,
            independent_free(t, s, d, volume),
            1e-14,
            3e-12,
        );
        near(
            state.free_energy_j,
            rotated_state.free_energy_j,
            1e-14,
            3e-12,
        );
        near(state.entropy_j_k, rotated_state.entropy_j_k, 1e-15, 3e-12);
        near(
            state.internal_energy_j,
            rotated_state.internal_energy_j,
            1e-12,
            3e-12,
        );
        near(
            state.heat_capacity_j_k,
            rotated_state.heat_capacity_j_k,
            1e-15,
            3e-12,
        );
        assert!(state.heat_capacity_j_k >= lower);
        let response = current.charge_response(s, d).unwrap();
        let mapped_response = mapped
            .charge_response(apply(&mandel_rotation(q), &s), apply(&q, &d))
            .unwrap();
        for (actual, expected) in mapped_response
            .stress_mandel_pa
            .iter()
            .zip(apply(&mandel_rotation(q), &response.stress_mandel_pa))
        {
            near(*actual, expected, 1e-7, 3e-12);
        }
        for (actual, expected) in mapped_response
            .electric_field_v_m
            .iter()
            .zip(apply(&q, &response.electric_field_v_m))
        {
            near(*actual, expected, 1e-8, 3e-12);
        }
    }
    let t = 315.0;
    let current = law.at_temperature(t).unwrap();
    let state = current.state(s, d, volume).unwrap();
    let response = current.charge_response(s, d).unwrap();
    let h = 0.02;
    let free = |t| independent_free(t, s, d, volume);
    near(
        state.entropy_j_k,
        -(free(t + h) - free(t - h)) / (2.0 * h),
        1e-12,
        3e-6,
    );
    near(
        state.heat_capacity_j_k,
        -t * (free(t + h) - 2.0 * free(t) + free(t - h)) / (h * h),
        1e-10,
        3e-6,
    );
    near(
        state.internal_energy_j,
        free(t) + t * state.entropy_j_k,
        1e-12,
        3e-12,
    );
    let h = 1e-8;
    for i in 0..6 {
        let mut plus = s;
        let mut minus = s;
        plus[i] += h;
        minus[i] -= h;
        let gradient = (independent_free(t, plus, d, volume)
            - independent_free(t, minus, d, volume))
            / (2.0 * h * volume);
        near(response.stress_mandel_pa[i], gradient, 1e-2, 3e-8);
    }
    for i in 0..3 {
        let mut plus = d;
        let mut minus = d;
        plus[i] += h;
        minus[i] -= h;
        let gradient = (independent_free(t, s, plus, volume)
            - independent_free(t, s, minus, volume))
            / (2.0 * h * volume);
        near(response.electric_field_v_m[i], gradient, 1e-4, 3e-8);
    }
}

#[test]
fn g1_directional_capacity_preserves_the_stress_free_caloric_reference() {
    let (mut c, _, _) = coefficients();
    for (i, row) in c.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            if i != j {
                *value = 0.0;
            }
        }
    }
    let e = [[0.0; 6], [0.0; 6], [0.0, 0.0, 12.0, 0.0, 0.0, 0.0]];
    let epsilon = [[2e-8, 0.0, 0.0], [0.0, 1.5e-8, 0.0], [0.0, 0.0, 1e-8]];
    let alpha = [0.0, 0.0, 1e-5, 0.0, 0.0, 0.0];
    let pyro = [0.0, 0.0, 1e-4];
    let card = bias_card(
        (c, e, epsilon),
        ([[0.0; 6]; 6], [[0.0; 6]; 3], [[0.0; 3]; 3]),
        alpha,
        pyro,
        |_| true,
    );
    let law = biased(&card);
    let volume = 1e-9;
    // Czz*alpha_zz² + pz²/epsilon_zz = 10+1 J/(m³ K²).
    let lower = law.minimum_heat_capacity_j_k(1e-3, 0.1, volume).unwrap();
    assert!(lower <= -340.0 * volume * 11.0);
    near(lower, -340.0 * volume * 11.0, 0.0, 1e-10);
    assert!(7500.0 * volume + lower > 0.0);
    for t in [280.0, 300.0, 320.0, 340.0] {
        let at = law.at_temperature(t).unwrap();
        let clamped = at.state([0.0; 6], [0.0; 3], volume).unwrap();
        near(
            clamped.free_energy_j,
            0.5 * volume * 11.0 * (t - 300.0).powi(2),
            1e-18,
            1e-13,
        );
        near(
            clamped.entropy_j_k,
            -volume * 11.0 * (t - 300.0),
            1e-20,
            1e-13,
        );
        near(
            clamped.internal_energy_j,
            -0.5 * volume * 11.0 * (t * t - 300.0_f64.powi(2)),
            1e-17,
            1e-13,
        );
        near(clamped.heat_capacity_j_k, -t * volume * 11.0, 1e-20, 1e-13);
        let free = at
            .field_response(at.thermal_strain_mandel(), [0.0; 3])
            .unwrap();
        let reference = at
            .state(at.thermal_strain_mandel(), free.displacement_c_m2, volume)
            .unwrap();
        near(reference.free_energy_j, 0.0, 1e-30, 0.0);
        near(reference.entropy_j_k, 0.0, 1e-20, 0.0);
        near(reference.internal_energy_j, 0.0, 1e-17, 0.0);
    }
}

#[test]
fn g0_bias_sources_refuse_missing_units_nonconstant_paths_and_implicit_reference_changes() {
    let card = constant_card();
    for reference in [279.0, 341.0, f64::NAN] {
        assert!(
            bind(&card)
                .with_thermal_bias(reference, MaterialPropertySelection::SingleClaimOnly)
                .is_err()
        );
    }
    assert!(
        biased(&card)
            .with_thermal_bias(300.0, MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    for case in 0..4 {
        let card = bias_card(
            coefficients(),
            ([[0.0; 6]; 6], [[0.0; 6]; 3], [[0.0; 3]; 3]),
            ALPHA,
            PYRO,
            |claim| {
                if claim.key.name() != "piezoelectric_pyroelectric_clamped_3" {
                    return true;
                }
                match case {
                    0 => return false,
                    1 => {
                        claim.key = PropertyKey::new(
                            "piezoelectric_pyroelectric_total_3",
                            claim.key.dims(),
                        );
                    }
                    2 => {
                        let dims = fs_material::state_point::INVERSE_TEMPERATURE_DIMS;
                        claim.key = PropertyKey::new(claim.key.name(), dims);
                        claim.value = PropertyValue::Scalar {
                            value: PYRO[2],
                            dims,
                        };
                    }
                    _ => {
                        claim.value = PropertyValue::Curve {
                            abscissa: "T".into(),
                            abscissa_dims: Temperature::DIMS,
                            knots: vec![(280.0, 0.8 * PYRO[2]), (340.0, 1.4 * PYRO[2])],
                            dims: claim.key.dims(),
                        };
                        claim.interpolation = InterpolationPolicy::LinearInside;
                    }
                }
                true
            },
        );
        assert!(
            bind(&card)
                .with_thermal_bias(300.0, MaterialPropertySelection::SingleClaimOnly)
                .is_err()
        );
    }
}
