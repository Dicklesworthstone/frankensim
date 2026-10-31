//! G1/G3 layerwise thermal-force, curvature, energy and frame oracles.
//! Synthetic numerical fixtures; no experimental material validation is claimed.

use fs_material::laminate::thermal::{
    LaminaThermalStrain, LaminateThermalError, PreparedLaminateThermal,
};
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use std::sync::Arc;

fn layer(e: f64, h: f64, angle: f64) -> LaminaLayer {
    LaminaLayer {
        thickness_m: h,
        orientation_rad: angle,
        density_kg_m3: 2.0,
        material: None,
        law: LayerLaw::PlaneStress {
            stiffness_pa: [e, 0.0, 0.0, 0.0, e, 0.0, 0.0, 0.0, 0.5 * e],
            transverse_shear_pa: [0.5 * e; 2],
        },
    }
}

fn prepare(layers: Vec<LaminaLayer>, profiles: &[LaminaThermalStrain]) -> PreparedLaminateThermal {
    let section = Arc::new(LaminateSection::try_new(layers, ReferenceSurface::Midsurface).unwrap());
    let prepared =
        PreparedLaminateThermal::try_new(Arc::clone(&section), profiles, || false).unwrap();
    assert!(core::ptr::eq(prepared.section(), section.as_ref()));
    prepared
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tol={tolerance:e}"
    );
}

#[test]
fn g1_bimetal_free_curvature_retains_incompatible_layer_energy_and_stress() {
    let (e, h) = (12.0, 1.0);
    for sign in [-1.0, 1.0] {
        let a = sign * 0.002;
        let b = sign * 0.006;
        let delta = b - a;
        let state = prepare(
            vec![layer(e, 0.5 * h, 0.0); 2],
            &[
                LaminaThermalStrain::uniform([a, a, 0.0]).unwrap(),
                LaminaThermalStrain::uniform([b, b, 0.0]).unwrap(),
            ],
        );
        near(state.membrane_force_n_m()[0], e * h * (a + b) / 2.0, 1e-16);
        near(state.bending_moment_n()[0], e * h * h * delta / 8.0, 1e-17);
        let membrane = [(a + b) / 2.0, (a + b) / 2.0, 0.0];
        let curvature = [1.5 * delta / h, 1.5 * delta / h, 0.0];
        let response = state.evaluate(membrane, curvature, || false).unwrap();
        for n in response.membrane_force_n_m {
            near(n, 0.0, 1e-16);
        }
        for m in response.bending_moment_n {
            near(m, 0.0, 1e-17);
        }
        // Independent integral of two linear residual strains: E*h*delta²/32
        // per normal direction. Net force/moment vanish, but storage does not.
        near(
            response.strain_energy_j_m2,
            e * h * delta * delta / 16.0,
            1e-19,
        );
        assert!(response.strain_energy_j_m2 > 0.0);
        let lower = state.layer_stress(0, membrane, curvature).unwrap();
        let upper = state.layer_stress(1, membrane, curvature).unwrap();
        near(lower.bottom_pa[0], -e * delta / 4.0, 1e-16);
        near(lower.top_pa[0], e * delta / 2.0, 1e-16);
        near(upper.bottom_pa[0], -e * delta / 2.0, 1e-16);
        near(upper.top_pa[0], e * delta / 4.0, 1e-16);
        near(state.section().i0_kg_m2(), 2.0 * h, 0.0);
    }
}

#[test]
fn g1_unequal_layer_stiffness_and_thickness_use_actual_abd_coupling() {
    let (e0, e1, h0, h1, a0, a1) = (8.0, 20.0, 0.4, 0.6, 0.001, 0.007);
    let state = prepare(
        vec![layer(e0, h0, 0.0), layer(e1, h1, 0.0)],
        &[
            LaminaThermalStrain::uniform([a0, 0.0, 0.0]).unwrap(),
            LaminaThermalStrain::uniform([a1, 0.0, 0.0]).unwrap(),
        ],
    );
    let spans = [(e0, -0.5_f64, -0.1_f64, a0), (e1, -0.1_f64, 0.5_f64, a1)];
    let mut a = 0.0;
    let mut b = 0.0;
    let mut d = 0.0;
    let mut n = 0.0;
    let mut m = 0.0;
    for (e, z0, z1, free) in spans {
        let q0 = z1 - z0;
        let q1 = (z1 * z1 - z0 * z0) / 2.0;
        let q2 = (z1.powi(3) - z0.powi(3)) / 3.0;
        a += e * q0;
        b += e * q1;
        d += e * q2;
        n += e * free * q0;
        m += e * free * q1;
    }
    let mid = (d * n - b * m) / (a * d - b * b);
    let curvature = (a * m - b * n) / (a * d - b * b);
    let response = state
        .evaluate([mid, 0.0, 0.0], [curvature, 0.0, 0.0], || false)
        .unwrap();
    near(response.membrane_force_n_m[0], 0.0, 2e-16);
    near(response.bending_moment_n[0], 0.0, 2e-17);
    let energy = spans
        .into_iter()
        .map(|(e, z0, z1, free)| {
            let r = mid - free;
            0.5 * e
                * (r * r * (z1 - z0)
                    + r * curvature * (z1 * z1 - z0 * z0)
                    + curvature * curvature * (z1.powi(3) - z0.powi(3)) / 3.0)
        })
        .sum::<f64>();
    near(response.strain_energy_j_m2, energy, 1e-18);
}

#[test]
fn g1_affine_thermal_profile_has_zero_stress_at_compatible_strain() {
    let mid = [0.003, -0.002, 0.001];
    let gradient = [0.01, -0.03, 0.02];
    let h = 0.08;
    let bottom = core::array::from_fn(|j| mid[j] - h * gradient[j] / 2.0);
    let top = core::array::from_fn(|j| mid[j] + h * gradient[j] / 2.0);
    let state = prepare(
        vec![layer(25.0, h, 0.0)],
        &[LaminaThermalStrain::from_faces(bottom, top).unwrap()],
    );
    let response = state.evaluate(mid, gradient, || false).unwrap();
    near(response.strain_energy_j_m2, 0.0, 1e-35);
    for n in response.membrane_force_n_m {
        near(n, 0.0, 1e-17);
    }
    for m in response.bending_moment_n {
        near(m, 0.0, 1e-18);
    }
    let stress = state.layer_stress(0, mid, gradient).unwrap();
    for s in stress.bottom_pa.into_iter().chain(stress.top_pa) {
        near(s, 0.0, 1e-16);
    }
}

#[test]
fn g3_reference_surface_shift_preserves_energy_and_physical_stress() {
    let layers = vec![layer(8.0, 0.4, 0.0), layer(20.0, 0.6, 0.3)];
    let profiles = [
        LaminaThermalStrain::from_faces([0.001, -0.002, 0.0], [0.004, 0.001, 0.001]).unwrap(),
        LaminaThermalStrain::from_faces([0.005, 0.002, -0.001], [0.007, 0.003, 0.0]).unwrap(),
    ];
    let a = prepare(layers.clone(), &profiles);
    let shift = -0.19;
    let section =
        Arc::new(LaminateSection::try_new(layers, ReferenceSurface::Offset(0.5 + shift)).unwrap());
    let b = PreparedLaminateThermal::try_new(section, &profiles, || false).unwrap();
    let mid = [0.003, 0.002, -0.001];
    let gradient = [0.01, -0.005, 0.002];
    let shifted = core::array::from_fn(|j| mid[j] + shift * gradient[j]);
    let ra = a.evaluate(mid, gradient, || false).unwrap();
    let rb = b.evaluate(shifted, gradient, || false).unwrap();
    near(ra.strain_energy_j_m2, rb.strain_energy_j_m2, 3e-19);
    for j in 0..3 {
        near(ra.membrane_force_n_m[j], rb.membrane_force_n_m[j], 3e-16);
        near(
            rb.bending_moment_n[j],
            ra.bending_moment_n[j] - shift * ra.membrane_force_n_m[j],
            5e-17,
        );
        near(
            b.bending_moment_n()[j],
            a.bending_moment_n()[j] - shift * a.membrane_force_n_m()[j],
            5e-17,
        );
    }
    for i in 0..2 {
        let sa = a.layer_stress(i, mid, gradient).unwrap();
        let sb = b.layer_stress(i, shifted, gradient).unwrap();
        for j in 0..3 {
            near(sa.bottom_pa[j], sb.bottom_pa[j], 3e-16);
            near(sa.top_pa[j], sb.top_pa[j], 3e-16);
        }
    }
}

#[test]
fn g3_rotated_orthotropic_layer_preserves_engineering_shear_and_work() {
    let h = 0.02;
    let mut ply = layer(1.0, h, core::f64::consts::FRAC_PI_4);
    ply.law = LayerLaw::PlaneStress {
        stiffness_pa: [90.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 12.0],
        transverse_shear_pa: [12.0, 15.0],
    };
    let state = prepare(
        vec![ply],
        &[LaminaThermalStrain::uniform([0.002, -0.001, 0.0008]).unwrap()],
    );
    // Exact 45-degree tensor rotation, including gamma=2 epsilon12.
    let free = [0.0001, 0.0009, 0.003];
    let stress = [0.075 - 0.0096, 0.075 + 0.0096, 0.105];
    let restrained = state.evaluate([0.0; 3], [0.0; 3], || false).unwrap();
    for j in 0..3 {
        near(state.membrane_force_n_m()[j], h * stress[j], 2e-18);
        near(restrained.membrane_force_n_m[j], -h * stress[j], 2e-18);
    }
    near(
        restrained.strain_energy_j_m2,
        0.5 * h * (90.0 * 0.002_f64.powi(2) + 30.0 * 0.001_f64.powi(2) + 12.0 * 0.0008_f64.powi(2)),
        2e-20,
    );
    near(
        state
            .evaluate(free, [0.0; 3], || false)
            .unwrap()
            .strain_energy_j_m2,
        0.0,
        1e-34,
    );
    // Work derivative in every membrane and curvature component.
    let membrane = [0.004, -0.003, 0.001];
    let gradient = [0.007, -0.003, 0.004];
    let response = state.evaluate(membrane, gradient, || false).unwrap();
    for dof in 0..6 {
        let mut plus = [membrane, gradient];
        let mut minus = plus;
        let step = 1e-5;
        plus[dof / 3][dof % 3] += step;
        minus[dof / 3][dof % 3] -= step;
        let derivative = (state
            .evaluate(plus[0], plus[1], || false)
            .unwrap()
            .strain_energy_j_m2
            - state
                .evaluate(minus[0], minus[1], || false)
                .unwrap()
                .strain_energy_j_m2)
            / (2.0 * step);
        let expected = if dof < 3 {
            response.membrane_force_n_m[dof]
        } else {
            response.bending_moment_n[dof - 3]
        };
        near(derivative, expected, 5e-15);
    }
}

#[test]
fn g0_bad_profiles_viscoelastic_history_and_cancellation_refuse() {
    assert!(LaminaThermalStrain::uniform([f64::NAN, 0.0, 0.0]).is_err());
    assert!(LaminaThermalStrain::from_faces([0.0; 3], [0.0, f64::INFINITY, 0.0]).is_err());
    let section = Arc::new(
        LaminateSection::try_new(vec![layer(20.0, 0.2, 0.0)], ReferenceSurface::Midsurface)
            .unwrap(),
    );
    let profile = LaminaThermalStrain::uniform([0.001, 0.002, 0.0]).unwrap();
    assert!(matches!(
        PreparedLaminateThermal::try_new(Arc::clone(&section), &[], || false),
        Err(LaminateThermalError::LayerCount { .. })
    ));
    assert!(matches!(
        PreparedLaminateThermal::try_new(Arc::clone(&section), &[profile], || true),
        Err(LaminateThermalError::Cancelled)
    ));
    let state = PreparedLaminateThermal::try_new(section, &[profile], || false).unwrap();
    assert!(matches!(
        state.evaluate([0.0; 3], [0.0; 3], || true),
        Err(LaminateThermalError::Cancelled)
    ));
    assert!(state.evaluate([0.0; 3], [f64::NAN; 3], || false).is_err());
    assert!(matches!(
        state.layer_stress(1, [0.0; 3], [0.0; 3]),
        Err(LaminateThermalError::LayerIndex { .. })
    ));
    near(
        state
            .evaluate([0.0; 3], [0.0; 3], || false)
            .unwrap()
            .strain_energy_j_m2,
        state.restrained_energy_j_m2(),
        0.0,
    );
    let diagonal =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 10.0 } else { 0.0 }));
    let prony =
        fs_material::visco::tensor::TensorProny::new(diagonal, None, vec![], (1.0, 100.0)).unwrap();
    let mut ply = layer(20.0, 0.2, 0.0);
    ply.law = LayerLaw::Prony(prony);
    let section =
        Arc::new(LaminateSection::try_new(vec![ply], ReferenceSurface::Midsurface).unwrap());
    assert!(matches!(
        PreparedLaminateThermal::try_new(section, &[profile], || false),
        Err(LaminateThermalError::ViscoelasticLayer { layer: 0 })
    ));
}

#[test]
fn g1_integrated_material_sources_drive_distinct_layers_without_endpoint_cte_shortcut() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::state_point::{
        INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        MaterialPropertySelection, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
        integrate_isotropic_thermal_expansion, integrate_orthotropic_thermal_expansion,
    };
    use fs_qty::Temperature;

    let mut claims = ClaimSet::new();
    for (name, scale) in [
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, 1.0),
        (ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[0], 2.0),
        (
            ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[1],
            -0.5,
        ),
        (ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[2], 3.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, INVERSE_TEMPERATURE_DIMS),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (300.0, scale * 1e-5),
                        (350.0, scale * 4e-5),
                        (400.0, scale * 2e-5),
                    ],
                    dims: INVERSE_TEMPERATURE_DIMS,
                },
                interpolation: InterpolationPolicy::LinearInside,
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic layered thermal integration oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test laminate".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let at = |t| QueryPoint::new().with("T", t).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let iso =
        integrate_isotropic_thermal_expansion(&card, &at(300.0), &at(400.0), selection.clone())
            .unwrap();
    let ortho =
        integrate_orthotropic_thermal_expansion(&card, &at(300.0), &at(400.0), selection.clone())
            .unwrap();
    let state = prepare(
        vec![layer(12.0, 0.5, 0.0); 2],
        &[
            LaminaThermalStrain::from_isotropic(&iso).unwrap(),
            LaminaThermalStrain::from_orthotropic(&ortho).unwrap(),
        ],
    );
    let integral = 0.00275; // two exact trapezoids, not alpha(400) * 100 = 0.002
    near(state.membrane_force_n_m()[0], 6.0 * 3.0 * integral, 2e-17);
    near(state.membrane_force_n_m()[1], 6.0 * 0.5 * integral, 2e-17);
    near(state.bending_moment_n()[0], 1.5 * integral, 2e-18);
    near(state.bending_moment_n()[1], -1.5 * 1.5 * integral, 2e-18);
    assert!(
        integrate_isotropic_thermal_expansion(&card, &at(300.0), &at(410.0), selection).is_err()
    );
}
