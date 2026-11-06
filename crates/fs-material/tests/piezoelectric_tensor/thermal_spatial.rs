//! G0/G1/G3 full fixed-(S,D) thermal source tests. Synthetic coefficients;
//! closed independent potentials and tensor elimination are the references.
use super::*;
use fs_material::piezoelectric::thermal::ThermalTensorPiezoelectric;
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ))
}
fn point_at(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", temperature(), t).unwrap()
}
fn thermal_card(
    base: (C, E, P),
    slopes: (C, E, P),
    mut alter: impl FnMut(&mut PropertyClaim),
) -> MaterialCard {
    let (cp, ep, pp) = slopes;
    source(base.0, base.1, base.2, |claim| {
        let PropertyValue::Scalar { value, dims } = claim.value else { unreachable!() };
        let name = claim.key.name();
        let slope = if name == "density" {
            0.0
        } else {
            let suffix = name.rsplit_once('_').unwrap().1.as_bytes();
            let (i, j) = ((suffix[0]-b'1') as usize, (suffix[1]-b'1') as usize);
            if name.starts_with("piezoelectric_stiffness_e_mandel") { cp[i][j] }
            else if name.starts_with("piezoelectric_stress_charge_mandel") { ep[i][j] }
            else { pp[i][j] }
        };
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, value-20.0*slope), (340.0, value+40.0*slope)], dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", temperature(), 280.0, 340.0);
        alter(claim);
    })
}
fn resolve_tensor(card: &MaterialCard) -> TensorPiezoelectric {
    TensorPiezoelectric::resolve(
        card, &point_at(300.0), MaterialPropertySelection::SingleClaimOnly, FRAME,
    ).unwrap()
}
fn bind(card: &MaterialCard) -> ThermalTensorPiezoelectric {
    ThermalTensorPiezoelectric::new(card, &resolve_tensor(card), "T", [280.0, 340.0]).unwrap()
}
fn at<const R: usize, const N: usize>(
    matrix: [[f64; N]; R], slope: [[f64; N]; R], t: f64,
) -> [[f64; N]; R] {
    core::array::from_fn(|i| core::array::from_fn(|j| {
        matrix[i][j]+(t-300.0)*slope[i][j]
    }))
}
fn norm<const N: usize>(x: &[f64; N]) -> f64 {
    x.iter().fold(0.0_f64, |norm, value| norm.hypot(*value))
}
fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a*b).sum()
}

#[test]
fn g1_full_tensor_thermal_potentials_match_independent_diagonal_rational_energy() {
    let mut c = [[0.0; 6]; 6];
    let mut cp = [[0.0; 6]; 6];
    for (i, row) in c.iter_mut().enumerate() {
        row[i] = [12e10, 11e10, 10e10, 5e10, 6e10, 7e10][i];
        cp[i][i] = [2e8, -1e8, 3e8, -2e7, 4e7, 1e7][i];
    }
    let mut e = [[0.0; 6]; 3];
    let mut ep = [[0.0; 6]; 3];
    for (i, j, value, slope) in [(0, 0, 2.0, 0.02), (1, 3, -3.0, -0.015), (2, 2, 14.0, 0.04)] {
        e[i][j] = value;
        ep[i][j] = slope;
    }
    let p: P = core::array::from_fn(|i| core::array::from_fn(|j| {
        if i == j { [2e-8, 1.7e-8, 1.1e-8][i] } else { 0.0 }
    }));
    let pp: P = core::array::from_fn(|i| core::array::from_fn(|j| {
        if i == j { [-2e-11, 1e-11, 2e-11][i] } else { 0.0 }
    }));
    let card = thermal_card((c, e, p), (cp, ep, pp), |_| {});
    let law = bind(&card);
    let strain = [2e-5, -1e-5, 3e-5, 8e-6, -3e-6, 2e-6];
    let displacement = [3e-4, -1e-4, 2e-4];
    let volume = 1e-5;
    let lower = law.minimum_heat_capacity_j_k(norm(&strain), norm(&displacement), volume).unwrap();
    for t in [280.0, 300.0, 325.0, 340.0] {
        let current_c = at(c, cp, t);
        let current_e = at(e, ep, t);
        let current_p = at(p, pp, t);
        let current = law.at_temperature(t).unwrap();
        law.at_material(current.material()).unwrap();
        assert_eq!(current.material().material().properties().len(), 64);
        assert_eq!(current.material().material().card_identity(), card.content_hash());
        assert_eq!(current.material().source_frame_identity(), FRAME);
        let state = current.state(strain, displacement, volume).unwrap();
        let mut free = 0.5*volume*(0..6).map(|i| current_c[i][i]*strain[i]*strain[i]).sum::<f64>();
        let mut entropy = -0.5*volume*(0..6).map(|i| cp[i][i]*strain[i]*strain[i]).sum::<f64>();
        let mut capacity = 0.0;
        for (i, j) in [(0, 0), (1, 3), (2, 2)] {
            // Independently differentiate (D_i-e_ij S_j)^2/(2 epsilon_i).
            // q' has the opposite sign to e', testing the mixed thermal term.
            let q = displacement[i]-current_e[i][j]*strain[j];
            let qp = -ep[i][j]*strain[j];
            let epsilon = current_p[i][i];
            let epsilon_prime = pp[i][i];
            free += 0.5*volume*q*q/epsilon;
            entropy -= volume*(q*qp/epsilon-0.5*q*q*epsilon_prime/(epsilon*epsilon));
            capacity -= t*volume/epsilon*(qp-q*epsilon_prime/epsilon).powi(2);
        }
        near(state.free_energy_j, free, 1e-18, 3e-12);
        near(state.entropy_j_k, entropy, 1e-21, 3e-12);
        near(state.internal_energy_j, free+t*entropy, 1e-18, 3e-12);
        near(state.heat_capacity_j_k, capacity, 1e-21, 3e-12);
        assert!(state.heat_capacity_j_k <= 0.0);
        assert!(state.heat_capacity_j_k >= lower);
    }
    assert_eq!(law.at_temperature(300.0).unwrap().material().identity(), law.original().identity());
    assert_eq!(law.card().content_hash(), card.content_hash());
    assert_eq!(law.temperature_axis(), "T");
    assert_eq!(law.temperature_range_k(), [280.0, 340.0]);
}

#[test]
fn g1_affine_full_stiffness_contributes_entropy_but_exactly_zero_fixed_strain_capacity() {
    let (c, e, p) = coefficients();
    let cp: C = core::array::from_fn(|i| core::array::from_fn(|j| 0.001*c[i][j]));
    let card = thermal_card((c, e, p), (cp, [[0.0; 6]; 3], [[0.0; 3]; 3]), |_| {});
    let law = bind(&card);
    let strain = [2e-5, -1e-5, 3e-5, 8e-6, -3e-6, 2e-6];
    let displacement = [3e-4, -1e-4, 2e-4];
    let volume = 1e-5;
    let entropy = -0.5*volume*dot(&strain, &apply(&cp, &strain));
    for t in [280.0, 300.0, 340.0] {
        let state = law.at_temperature(t).unwrap().state(strain, displacement, volume).unwrap();
        near(state.entropy_j_k, entropy, 1e-21, 3e-12);
        assert_eq!(state.heat_capacity_j_k, 0.0);
        near(state.internal_energy_j, state.free_energy_j+t*entropy, 1e-18, 3e-12);
    }
    assert_eq!(law.minimum_heat_capacity_j_k(norm(&strain), norm(&displacement), volume).unwrap(), 0.0);
    let zero = law.at_temperature(300.0).unwrap().state([0.0; 6], [0.0; 3], volume).unwrap();
    assert_eq!((zero.free_energy_j, zero.entropy_j_k, zero.internal_energy_j, zero.heat_capacity_j_k),
        (0.0, 0.0, 0.0, 0.0));
}

fn dense_slopes(c: C, e: E) -> (C, E, P) {
    let mut cp: C = core::array::from_fn(|i| core::array::from_fn(|j| {
        0.0002*(i+j+1) as f64*c[i][j]
    }));
    cp[0][4] = 2e7;
    cp[4][0] = cp[0][4];
    cp[2][5] = -4e7;
    cp[5][2] = cp[2][5];
    let ep: E = core::array::from_fn(|i| core::array::from_fn(|j| {
        0.001*(2+i) as f64*e[i][j]+0.003*(j as f64-i as f64)
    }));
    let pp = [[2e-11, 5e-12, 2e-12], [5e-12, -3e-11, -3e-12], [2e-12, -3e-12, 1e-11]];
    (cp, ep, pp)
}

// Independent three-dimensional partial-pivot elimination, not the native
// material's scaled Cholesky factor or its thermal state implementation.
fn solve_dielectric(mut a: P, mut b: [f64; 3]) -> [f64; 3] {
    for j in 0..3 {
        let pivot = (j..3).max_by(|&i, &k| a[i][j].abs().total_cmp(&a[k][j].abs())).unwrap();
        a.swap(j, pivot);
        b.swap(j, pivot);
        for i in j+1..3 {
            let ratio = a[i][j]/a[j][j];
            for k in j..3 {
                let pivot_value = a[j][k];
                a[i][k] -= ratio*pivot_value;
            }
            let pivot_value = b[j];
            b[i] -= ratio*pivot_value;
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        x[i] = (b[i]-(i+1..3).map(|j| a[i][j]*x[j]).sum::<f64>())/a[i][i];
    }
    x
}
fn independent_free(t: f64, strain: [f64; 6], displacement: [f64; 3], volume: f64) -> f64 {
    let (c, e, p) = coefficients();
    let (cp, ep, pp) = dense_slopes(c, e);
    let c = at(c, cp, t);
    let e = at(e, ep, t);
    let p = at(p, pp, t);
    let polarization = apply(&e, &strain);
    let charge = core::array::from_fn(|i| displacement[i]-polarization[i]);
    let field = solve_dielectric(p, charge);
    0.5*volume*(dot(&strain, &apply(&c, &strain))+dot(&charge, &field))
}

#[test]
fn g1_g3_full_thermal_source_matches_independent_tensor_derivatives_and_corotation() {
    let base = coefficients();
    let slopes = dense_slopes(base.0, base.1);
    let card = thermal_card(base, slopes, |_| {});
    let tensor = resolve_tensor(&card);
    let law = bind(&card);
    let q = rotation([1.0/3.0_f64.sqrt(); 3], 0.73);
    let rotated = tensor.rotated(ContentHash([9; 32]), q).unwrap();
    let moved = ThermalTensorPiezoelectric::new(&card, &rotated, "T", [280.0, 340.0]).unwrap();
    let strain = [2e-5, -1e-5, 3e-5, 8e-6, -3e-6, 2e-6];
    let displacement = [3e-4, -1e-4, 2e-4];
    let volume = 1e-5;
    let strain_norm = norm(&strain);
    let displacement_norm = norm(&displacement);
    let lower = law.minimum_heat_capacity_j_k(strain_norm, displacement_norm, volume).unwrap();
    let moved_lower = moved.minimum_heat_capacity_j_k(strain_norm, displacement_norm, volume).unwrap();
    near(lower, moved_lower, 1e-18, 3e-12);
    near(law.minimum_heat_capacity_j_k(2.0*strain_norm, 2.0*displacement_norm, volume).unwrap(),
        4.0*lower, 1e-18, 3e-12);
    near(law.minimum_heat_capacity_j_k(strain_norm, displacement_norm, 2.0*volume).unwrap(),
        2.0*lower, 1e-18, 3e-12);
    for t in [280.0, 295.0, 300.0, 320.0, 340.0] {
        let current = law.at_temperature(t).unwrap();
        let mapped = moved.at_temperature(t).unwrap();
        law.at_material(current.material()).unwrap();
        moved.at_material(mapped.material()).unwrap();
        assert!(law.at_material(mapped.material()).is_err());
        for (s_sign, d_sign) in [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)] {
            let s = strain.map(|x| x*s_sign);
            let d = displacement.map(|x| x*d_sign);
            let state = current.state(s, d, volume).unwrap();
            let rotated_state = mapped.state(apply(&mandel_rotation(q), &s), apply(&q, &d), volume).unwrap();
            near(state.free_energy_j, independent_free(t, s, d, volume), 1e-18, 3e-12);
            near(state.free_energy_j, rotated_state.free_energy_j, 1e-18, 3e-12);
            near(state.entropy_j_k, rotated_state.entropy_j_k, 1e-21, 3e-12);
            near(state.internal_energy_j, rotated_state.internal_energy_j, 1e-18, 3e-12);
            near(state.heat_capacity_j_k, rotated_state.heat_capacity_j_k, 1e-21, 3e-12);
            assert!(state.heat_capacity_j_k <= 0.0);
            assert!(state.heat_capacity_j_k >= lower);
        }
    }
    let t = 315.0;
    let h = 0.05;
    let state = law.at_temperature(t).unwrap().state(strain, displacement, volume).unwrap();
    let free = |t| independent_free(t, strain, displacement, volume);
    near(state.entropy_j_k, -(free(t+h)-free(t-h))/(2.0*h), 1e-20, 3e-7);
    near(state.heat_capacity_j_k, -t*(free(t+h)-2.0*free(t)+free(t-h))/(h*h), 1e-16, 3e-6);
    near(state.internal_energy_j, free(t)+t*state.entropy_j_k, 1e-18, 3e-12);
}

#[test]
fn g0_full_thermal_tensor_retains_sources_and_refuses_unsupported_or_unstable_paths() {
    let base = coefficients();
    let slopes = dense_slopes(base.0, base.1);
    let card = thermal_card(base, slopes, |_| {});
    let law = bind(&card);
    assert!(law.at_temperature(279.0).is_err());
    assert!(law.at_temperature(f64::NAN).is_err());
    for (s, d, v) in [(-1.0, 0.0, 1.0), (0.0, -1.0, 1.0), (0.0, 0.0, 0.0),
        (f64::NAN, 0.0, 1.0), (0.0, f64::INFINITY, 1.0)] {
        assert!(law.minimum_heat_capacity_j_k(s, d, v).is_err());
    }
    let current = law.at_temperature(300.0).unwrap();
    assert!(current.state([f64::NAN; 6], [0.0; 3], 1.0).is_err());
    assert!(current.state([0.0; 6], [0.0; 3], -1.0).is_err());
    let foreign = thermal_card(base, slopes, |claim| {
        if claim.key.name() == "density" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                for (_, value) in knots { *value = 7600.0; }
            }
        }
    });
    assert!(law.at_material(&resolve_tensor(&foreign)).is_err());
    let nonlinear = thermal_card(base, slopes, |claim| {
        if claim.key.name() == "piezoelectric_stress_charge_mandel_31" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                knots.insert(1, (300.0, -5.99));
            }
        }
    });
    assert!(ThermalTensorPiezoelectric::new(
        &nonlinear, &resolve_tensor(&nonlinear), "T", [280.0, 340.0],
    ).is_err());
    let unstable = thermal_card(base, slopes, |claim| {
        if ["piezoelectric_permittivity_s_12", "piezoelectric_permittivity_s_21"]
            .contains(&claim.key.name()) {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                *knots = vec![(280.0, -4e-8), (340.0, 8e-8)];
            }
        }
    });
    assert!(ThermalTensorPiezoelectric::new(
        &unstable, &resolve_tensor(&unstable), "T", [280.0, 340.0],
    ).is_err());
    let changing_density = thermal_card(base, slopes, |claim| {
        if claim.key.name() == "density" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value { knots[1].1 += 100.0; }
        }
    });
    assert!(ThermalTensorPiezoelectric::new(
        &changing_density, &resolve_tensor(&changing_density), "T", [280.0, 340.0],
    ).is_err());
    // Constant source coefficients still require explicit absolute-temperature
    // semantics; an untyped legacy coordinate cannot acquire kelvin by guessing.
    let constant = source(base.0, base.1, base.2, |_| {});
    let untyped = TensorPiezoelectric::resolve(&constant, &point(),
        MaterialPropertySelection::SingleClaimOnly, FRAME).unwrap();
    assert!(ThermalTensorPiezoelectric::new(
        &constant, &untyped, "T", [280.0, 340.0],
    ).is_err());
}

#[test]
fn g0_full_thermal_admits_original_and_independently_resolved_same_claim_policies() {
    let base = coefficients();
    let card = thermal_card(base, dense_slopes(base.0, base.1), |_| {});
    let original = resolve_tensor(&card);
    let law = ThermalTensorPiezoelectric::new(&card, &original, "T", [280.0, 340.0]).unwrap();
    let original_at = law.at_material(&original).unwrap();
    assert_eq!(original_at.material().identity(), original.identity());
    assert_eq!(law.at_temperature(300.0).unwrap().material().identity(), original.identity());
    let pins = original.material().properties().iter().map(|property| {
        (property.requirement().name().to_owned(), property.answer().receipt.selected)
    }).collect::<Vec<_>>();
    let strain = [2e-5, -1e-5, 3e-5, 8e-6, -3e-6, 2e-6];
    let displacement = [3e-4, -1e-4, 2e-4];
    for temperature in [300.0, 320.0] {
        let expected = law.at_temperature(temperature).unwrap()
            .state(strain, displacement, 1e-5).unwrap();
        for policy in [
            MaterialPropertySelection::SingleClaimOnly,
            MaterialPropertySelection::PinnedByProperty(pins.clone()),
        ] {
            let independently_resolved = TensorPiezoelectric::resolve(
                &card, &point_at(temperature), policy.clone(), FRAME,
            ).unwrap();
            if temperature == 300.0 && matches!(policy, MaterialPropertySelection::PinnedByProperty(_)) {
                // This fixture must exercise distinct receipt-policy hashes.
                assert_ne!(independently_resolved.material().identity(), original.material().identity());
            }
            let admitted = law.at_material(&independently_resolved).unwrap();
            assert_eq!(admitted.material().identity(), independently_resolved.identity());
            let actual = admitted.state(strain, displacement, 1e-5).unwrap();
            assert_eq!(actual.free_energy_j, expected.free_energy_j);
            assert_eq!(actual.entropy_j_k, expected.entropy_j_k);
            assert_eq!(actual.internal_energy_j, expected.internal_energy_j);
            assert_eq!(actual.heat_capacity_j_k, expected.heat_capacity_j_k);
        }
    }
}
