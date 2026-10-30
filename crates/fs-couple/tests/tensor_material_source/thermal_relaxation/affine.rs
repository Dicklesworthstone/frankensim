//! Independent Helmholtz/entropy/internal-energy checks, including real cooling.
use super::*;
use fs_couple::tensor_material::solid::thermal::relaxation::AffineThermalBudget;
use fs_material::visco_tensor::source::thermal::AffineTensorThermal;

fn branch_slope() -> [[f64; 6]; 6] {
    let mut m = [[0.0; 6]; 6];
    for (i, v) in [0.8, -0.1, 0.4, 0.3, -0.2, 0.5].into_iter().enumerate() {
        m[i][i] = v;
    }
    m[0][1] = 0.05;
    m[1][0] = 0.05;
    m[0][2] = 0.2;
    m[2][0] = 0.2;
    m
}
fn claims(with_branch: bool) -> (Vec<PropertyClaim>, TensorPronyProperties) {
    let (mut claims, properties) = fixture(1000.0, 10.0, if with_branch { 0.4 } else { 0.0 });
    claims.push(scalar(
        "specific_heat_capacity",
        Dims([2, 0, -2, -1, 0, 0]),
        CP,
    ));
    let mut slopes = BTreeMap::from([("young_modulus".to_owned(), -2.0)]);
    if let Some(branch) = properties.branches.first() {
        let matrix = branch_slope();
        let mut entry = 0;
        for row in 0..6 {
            for column in row..6 {
                slopes.insert(branch.stiffness_upper[entry].clone(), matrix[row][column]);
                entry += 1;
            }
        }
    }
    for claim in &mut claims {
        let Some(&slope) = slopes.get(claim.key.name()).filter(|s| **s != 0.0) else {
            continue;
        };
        let (value, dims) = match claim.value {
            PropertyValue::Scalar { value, dims } => (value, dims),
            _ => unreachable!(),
        };
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: [280.0, 300.0, 310.0, 330.0]
                .map(|t| (t, value + slope * (t - 300.0)))
                .to_vec(),
            dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
    }
    (claims, properties)
}
fn budget() -> AffineThermalBudget {
    AffineThermalBudget {
        max_absolute_coefficient_work_j: 0.1,
        relative_memory_tolerance: 1e-11,
    }
}

// Physical strain coordinates: x, xdot, Maxwell spring strain, accumulated
// dashpot heat and force work. The reference owns no matrix factor or FE code.
fn trajectory(mut y: [f64; 5], t: f64, dt: f64, branch: bool) -> [f64; 5] {
    let e = 1000.0 - 2.0 * (t - 300.0);
    let b = if branch {
        280.0 + 0.4 * (t - 300.0)
    } else {
        0.0
    };
    let rhs = |z: [f64; 5]| {
        [
            z[1],
            (0.02 - (e * z[0] + b * z[2]) / 6.0) / (10.0 / 60.0),
            z[1] - z[2] / 0.02,
            b * z[2] * z[2] / (6.0 * 0.02),
            0.02 * z[1],
        ]
    };
    let h = dt / 2048.0;
    for _ in 0..2048 {
        let a = rhs(y);
        let b = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| y[i] + h * c[i]));
        y = core::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}
fn entropy(y: [f64; 5], branch: bool) -> f64 {
    -(-2.0 * y[0] * y[0] + if branch { 0.4 * y[2] * y[2] } else { 0.0 }) / 12.0
}
fn potential(y: [f64; 5], t: f64, branch: bool) -> f64 {
    ((1000.0 - 2.0 * (t - 300.0)) * y[0] * y[0]
        + if branch {
            (280.0 + 0.4 * (t - 300.0)) * y[2] * y[2]
        } else {
            0.0
        })
        / 12.0
}
fn endpoint(y: [f64; 5], t0: f64, dt: f64, branch: bool) -> ([f64; 5], f64, f64) {
    let mut lo = 280.0;
    let mut hi = 330.0;
    for _ in 0..48 {
        let t = 0.5 * (lo + hi);
        let next = trajectory(y, t, dt, branch);
        let reversible = t * (entropy(y, branch) - entropy(next, branch));
        let residual = (10.0 / 6.0) * CP * (t - t0) - (next[3] - y[3] + reversible);
        if residual > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = 0.5 * (lo + hi);
    let next = trajectory(y, t, dt, branch);
    (next, t, t * (entropy(y, branch) - entropy(next, branch)))
}

#[test]
fn affine_tensors_close_physical_internal_energy_and_preserve_all_native_accounts() {
    with_cx(|cx, gate| {
        let (claims, properties) = claims(true);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let initial = body
            .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
            .unwrap();
        let original = initial.material_history().state().to_vec();
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let mut thermal = heat(&[&curve], &[CP], 0.0);
        let mut model: Option<RelaxingSolid<'_>> = None;
        let mut y = [0.006, 0.08, 0.006, 0.0, 0.0];
        let mut temperature = 300.0;
        let mut coefficient = 0.0;
        for target in [0.005, 0.01] {
            let history = model
                .as_ref()
                .map(RelaxingSolid::history)
                .unwrap_or_else(|| initial.clone());
            let (next, t, rev) = endpoint(y, temperature, 0.005, true);
            let result = history
                .advance_affine_thermal_candidate(
                    &thermal,
                    target,
                    &FORCE,
                    &cells,
                    controls(),
                    budget(),
                    &[0.00005],
                    no_extra,
                    gate,
                )
                .unwrap();
            let actual = result.model().history();
            let report = result.report();
            let thermo = report.thermodynamics.as_ref().unwrap();
            near(actual.material_history().displacement_m()[0], next[0], 2e-9);
            near(actual.material_history().velocity_m_s()[0], next[1], 2e-8);
            near(
                actual.material_history().dissipated_heat_j(),
                next[3],
                2e-11,
            );
            near(result.thermal().temperature_k(CellId(0)).unwrap(), t, 2e-7);
            near(thermo.before[0].entropy_j_k, entropy(y, true), 2e-13);
            near(thermo.after[0].entropy_j_k, entropy(next, true), 2e-13);
            near(
                thermo.after[0].free_energy_j,
                potential(next, t, true),
                2e-11,
            );
            near(
                thermo.after[0].internal_energy_j,
                potential(next, t, true) + t * entropy(next, true),
                2e-11,
            );
            near(thermo.reversible_heat_j[0], rev, 3e-11);
            assert!(thermo.reversible_heat_j[0] < 0.0);
            assert!(report.mechanical.element_heat_j[0] > 0.0);
            near(
                report.deposited_material_heat_j,
                next[3] - y[3] + rev,
                3e-11,
            );
            let p = -(t - temperature) * entropy(y, true);
            coefficient += p;
            near(thermo.coefficient_work_j[0], p, 2e-11);
            near(
                actual.material_history().thermal_coefficient_work_j(),
                coefficient,
                3e-11,
            );
            assert_eq!(actual.material_history().parameter_work_j(), 0.0);
            assert!(actual.material_history().balance_residual_j().abs() < 2e-10);
            let kinetic = 0.5 * (10.0 / 60.0) * (next[1] * next[1] - y[1] * y[1]);
            let physical = kinetic + thermo.after[0].internal_energy_j
                - thermo.before[0].internal_energy_j
                + result.thermal().energy_j()
                - thermal.energy_j()
                - (next[4] - y[4]);
            near(physical, 0.0, 3e-10);
            assert!(report.total_energy_residual_j.abs() < 2e-10);
            assert_eq!(
                result.model().thermodynamics().unwrap(),
                thermo.after.as_slice()
            );
            let (new_model, new_thermal, _) = result.into_parts();
            model = Some(new_model);
            thermal = new_thermal;
            y = next;
            temperature = t;
        }
        assert_eq!(
            model.as_ref().unwrap().history().material_history().steps(),
            2
        );
        assert_eq!(initial.material_history().state(), original);
        assert_eq!(initial.material_history().thermal_coefficient_work_j(), 0.0);
    });
}

#[test]
fn a_lossless_affine_elastic_body_can_cool_its_finite_caloric_receiver() {
    with_cx(|cx, gate| {
        let (claims, properties) = claims(false);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let history = body
            .history(&[0.006], &[0.08], InitialMemory::Relaxed, gate)
            .unwrap();
        let thermal = heat(&[&curve], &[CP], 0.0);
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let (expected, t, rev) = endpoint([0.006, 0.08, 0.0, 0.0, 0.0], 300.0, 0.005, false);
        let next = history
            .advance_affine_thermal_candidate(
                &thermal,
                0.005,
                &FORCE,
                &cells,
                controls(),
                budget(),
                &[-0.0002],
                no_extra,
                gate,
            )
            .unwrap();
        assert!(t < 300.0);
        assert!(rev < 0.0);
        assert_eq!(next.report().mechanical.dynamics.dissipated, 0.0);
        assert!(next.report().deposited_material_heat_j < 0.0);
        near(next.thermal().temperature_k(CellId(0)).unwrap(), t, 2e-7);
        near(
            next.model().history().material_history().displacement_m()[0],
            expected[0],
            2e-9,
        );
        near(
            next.report()
                .thermodynamics
                .as_ref()
                .unwrap()
                .reversible_heat_j[0],
            rev,
            3e-11,
        );
        near(next.report().external_thermal_heat_j, 0.0, 1e-12);
        assert!(next.report().total_energy_residual_j.abs() < 2e-10);
    });
}

#[test]
fn full_anisotropic_entropy_uses_physical_mandel_strain_and_the_same_rotation() {
    let (claims, properties) = claims(true);
    let card = card(claims);
    let source = bind(&card, properties);
    let law = AffineTensorThermal::try_new(&source, &card, "T", [280.0, 330.0]).unwrap();
    let e = [0.002, -0.003, 0.001, 0.004, -0.001, 0.002];
    let mandel =
        core::array::from_fn::<_, 6, _>(|i| if i < 3 { e[i] } else { 2.0_f64.sqrt() * e[i] });
    let slope = branch_slope();
    let dot = |m: [[f64; 6]; 6]| {
        (0..6)
            .flat_map(|i| (0..6).map(move |j| mandel[i] * m[i][j] * mandel[j]))
            .sum::<f64>()
    };
    let expected_entropy =
        -0.5 * 0.2 * (-2.0 * mandel.iter().map(|x| x * x).sum::<f64>() + dot(slope));
    let mut samples = Vec::new();
    for (rotation, strain) in [
        (IDENTITY, e),
        (TURN, [e[1], e[0], e[2], -e[3], e[5], -e[4]]),
    ] {
        let material = source.rotated_law(rotation).unwrap();
        let columns: Vec<_> = (0..6)
            .map(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
            .collect();
        let projected = material
            .project_strain(&columns, 0.2, projection(), || false)
            .unwrap();
        let memory: Vec<_> = projected
            .branches()
            .iter()
            .flat_map(|branch| {
                branch
                    .rows()
                    .iter()
                    .map(|row| row.iter().zip(strain).map(|(a, e)| a * e).sum::<f64>())
            })
            .collect();
        let sample = projected
            .thermal_state(
                &law.in_frame(rotation).unwrap(),
                &strain,
                &memory,
                300.0,
                1e-11,
                || false,
            )
            .unwrap();
        near(sample.entropy_j_k, expected_entropy, 1e-13);
        near(
            sample.internal_energy_j,
            sample.free_energy_j + 300.0 * expected_entropy,
            1e-11,
        );
        samples.push(sample);
    }
    near(samples[0].free_energy_j, samples[1].free_energy_j, 1e-12);
    near(
        samples[0].internal_energy_j,
        samples[1].internal_energy_j,
        1e-12,
    );
}

#[test]
fn orthotropic_equilibrium_uses_one_affine_scale_and_refuses_nonproportional_moduli() {
    use fs_material::state_point::{
        ORTHOTROPIC_POISSON_RATIO_PROPERTIES, ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES,
        ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    };
    let build = |nonproportional| {
        let mut claims = vec![
            scalar("density", Density::DIMS, 10.0),
            scalar("strain_limit", Dims::NONE, 0.1),
        ];
        for name in ORTHOTROPIC_POISSON_RATIO_PROPERTIES {
            claims.push(scalar(name, Dims::NONE, 0.0));
        }
        for (index, (name, value)) in ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES
            .into_iter()
            .chain(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES)
            .zip([1200.0, 700.0, 400.0, 250.0, 150.0, 200.0])
            .enumerate()
        {
            let mut claim = scalar(name, Pressure::DIMS, value);
            let slope = if nonproportional && index == 5 {
                -0.001
            } else {
                -0.002
            };
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: [280.0, 300.0, 330.0]
                    .map(|t| (t, value * (1.0 + slope * (t - 300.0))))
                    .to_vec(),
                dims: Pressure::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            claims.push(claim);
        }
        let card = card(claims);
        let source = bind(
            &card,
            TensorPronyProperties {
                elasticity: TensorElasticity::Orthotropic,
                strain_limit: "strain_limit".into(),
                branches: vec![],
            },
        );
        (card, source)
    };
    let (card, source) = build(false);
    let law = AffineTensorThermal::try_new(&source, &card, "T", [280.0, 330.0]).unwrap();
    let strain = [0.002, -0.003, 0.001, 0.004, -0.001, 0.002];
    let columns: Vec<_> = (0..6)
        .map(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
        .collect();
    let projected = source
        .law()
        .project_strain(&columns, 0.2, projection(), || false)
        .unwrap();
    let state = projected
        .thermal_state(
            &law.in_frame(IDENTITY).unwrap(),
            &strain,
            &[],
            300.0,
            1e-11,
            || false,
        )
        .unwrap();
    // Physical tensor shear stores 2*G*epsilon_ij^2 per shear direction.
    let psi = 0.1
        * (1200.0 * strain[0] * strain[0]
            + 700.0 * strain[1] * strain[1]
            + 400.0 * strain[2] * strain[2]
            + 4.0
                * (250.0 * strain[3] * strain[3]
                    + 150.0 * strain[4] * strain[4]
                    + 200.0 * strain[5] * strain[5]));
    near(state.free_energy_j, psi, 1e-13);
    near(state.entropy_j_k, 0.002 * psi, 1e-13);
    near(state.internal_energy_j, 1.6 * psi, 1e-13);
    let (card, source) = build(true);
    assert!(AffineTensorThermal::try_new(&source, &card, "T", [280.0, 330.0]).is_err());
}

#[test]
fn nonaffine_sources_and_late_failed_trials_leave_every_accepted_account_unchanged() {
    with_cx(|cx, gate| {
        for key in ["young_modulus", "branch_c00", "density", "poisson_ratio"] {
            let (mut claims, properties) = claims(true);
            let claim = claims.iter_mut().find(|c| c.key.name() == key).unwrap();
            match &mut claim.value {
                PropertyValue::Curve { knots, .. } => knots[2].1 += 1.0,
                PropertyValue::Scalar { value, dims } => {
                    claim.value = PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![
                            (280.0, *value),
                            (300.0, *value),
                            (330.0, *value + if key == "density" { 1.0 } else { 0.05 }),
                        ],
                        dims: *dims,
                    };
                    claim.interpolation = InterpolationPolicy::LinearInside;
                }
                _ => unreachable!(),
            }
            let card = card(claims);
            let source = bind(&card, properties);
            assert!(AffineTensorThermal::try_new(&source, &card, "T", [280.0, 330.0]).is_err());
        }
        let (claims, properties) = claims(true);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let history = body
            .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
            .unwrap();
        let thermal = heat(&[&curve], &[CP], 0.0);
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let state = history.material_history().state().to_vec();
        let energy = thermal.energy_j();
        let tight = AffineThermalBudget {
            max_absolute_coefficient_work_j: 0.0,
            ..budget()
        };
        assert!(
            history
                .advance_affine_thermal_candidate(
                    &thermal,
                    0.005,
                    &FORCE,
                    &cells,
                    controls(),
                    tight,
                    &[0.00005],
                    no_extra,
                    gate
                )
                .is_err()
        );
        let late = CancelGate::new_clock_free();
        let failed = history.advance_affine_thermal_candidate(
            &thermal,
            0.005,
            &FORCE,
            &cells,
            controls(),
            budget(),
            &[0.00005],
            |t| {
                late.request();
                no_extra(t)
            },
            &late,
        );
        assert!(matches!(
            failed,
            Err(fs_couple::circuit_enthalpy::CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(history.material_history().state(), state);
        assert_eq!(thermal.energy_j(), energy);
        assert_eq!(history.material_history().thermal_coefficient_work_j(), 0.0);
        assert_eq!(history.material_history().parameter_work_j(), 0.0);
        let first = history
            .advance_affine_thermal_candidate(
                &thermal,
                0.005,
                &FORCE,
                &cells,
                controls(),
                budget(),
                &[0.00005],
                no_extra,
                gate,
            )
            .unwrap();
        let retry = history
            .advance_affine_thermal_candidate(
                &thermal,
                0.005,
                &FORCE,
                &cells,
                controls(),
                budget(),
                &[0.00005],
                no_extra,
                gate,
            )
            .unwrap();
        assert_eq!(
            first.model().history().material_history().state(),
            retry.model().history().material_history().state()
        );
        assert_eq!(first.thermal().energy_j(), retry.thermal().energy_j());
    });
}
