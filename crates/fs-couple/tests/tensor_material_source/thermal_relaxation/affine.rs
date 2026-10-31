//! Independent Helmholtz/entropy/internal-energy checks, including real cooling.
use super::*;
use fs_couple::tensor_material::solid::thermal::relaxation::AffineThermalBudget;
use fs_material::visco_tensor::source::thermal::{AffineTensorThermal, TensorThermalLaw};

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
fn trajectory(y: [f64; 5], t: f64, dt: f64, branch: bool) -> [f64; 5] {
    let e = 1000.0 - 2.0 * (t - 300.0);
    let b = if branch {
        280.0 + 0.4 * (t - 300.0)
    } else {
        0.0
    };
    trajectory_with_stiffness(y, dt, e, b)
}
fn trajectory_with_stiffness(mut y: [f64; 5], dt: f64, e: f64, b: f64) -> [f64; 5] {
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

const RATIONAL_CP: f64 = 0.001;
const RATIONAL_HEATER_W: f64 = 2.0;

fn rational_claims(cp: f64, with_branch: bool) -> (Vec<PropertyClaim>, TensorPronyProperties) {
    let (mut claims, mut properties) = claims(with_branch);
    properties.elasticity = TensorElasticity::IsotropicYoungShear;
    for claim in &mut claims {
        if let PropertyValue::Scalar { value, .. } = &mut claim.value {
            match claim.key.name() {
                "specific_heat_capacity" => *value = cp,
                "strain_limit" => *value = 0.01,
                _ => {}
            }
        }
    }
    let mut shear = scalar("shear_modulus", Pressure::DIMS, 500.0);
    shear.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: fs_qty::Temperature::DIMS,
        knots: [280.0, 300.0, 310.0, 330.0]
            .map(|t| (t, 500.0 + 0.5 * (t - 300.0)))
            .to_vec(),
        dims: Pressure::DIMS,
    };
    shear.interpolation = InterpolationPolicy::LinearInside;
    claims.push(shear);
    (claims, properties)
}

// Independent scalar Young/shear constitutive algebra. No source derivatives,
// projected factors or finite differences enter the reference thermodynamics.
fn rational_bulk(t: f64) -> [f64; 3] {
    let young = 1000.0 - 2.0 * (t - 300.0);
    let shear = 500.0 + 0.5 * (t - 300.0);
    let gap = 3.0 * shear - young;
    let cross = -2.0 * shear - 0.5 * young;
    [
        young * shear / (3.0 * gap),
        ((-2.0 * shear + 0.5 * young) * gap - 3.5 * young * shear) / (3.0 * gap * gap),
        2.0 * cross * cross / (gap * gap * gap),
    ]
}
fn rational_axial(t: f64) -> [f64; 3] {
    let [bulk, slope, curvature] = rational_bulk(t);
    [
        bulk + (4.0 / 3.0) * (500.0 + 0.5 * (t - 300.0)),
        slope + 2.0 / 3.0,
        curvature,
    ]
}
fn rational_state(y: [f64; 5], t: f64) -> [f64; 3] {
    let [modulus, slope, curvature] = rational_axial(t);
    let branch = 280.0 + 0.4 * (t - 300.0);
    [
        (modulus * y[0].powi(2) + branch * y[2].powi(2)) / 12.0,
        -(slope * y[0].powi(2) + 0.4 * y[2].powi(2)) / 12.0,
        -t * curvature * y[0].powi(2) / 12.0,
    ]
}
fn rational_endpoint(y: [f64; 5], t0: f64, dt: f64) -> ([f64; 5], f64, f64, f64) {
    let old = rational_state(y, t0);
    let at = |t: f64| {
        let next =
            trajectory_with_stiffness(y, dt, rational_axial(t)[0], 280.0 + 0.4 * (t - 300.0));
        let coefficient = rational_state(y, t)[0] - old[0];
        let reversible = -coefficient - t * rational_state(next, t)[1] + t0 * old[1];
        (next, coefficient, reversible)
    };
    let (mut low, mut high) = (280.0, 330.0);
    for _ in 0..48 {
        let t = 0.5 * (low + high);
        let (next, _, reversible) = at(t);
        let residual = (10.0 / 6.0) * RATIONAL_CP * (t - t0)
            - (next[3] - y[3] + reversible + RATIONAL_HEATER_W * dt);
        if residual > 0.0 {
            high = t;
        } else {
            low = t;
        }
    }
    let t = 0.5 * (low + high);
    let (next, coefficient, reversible) = at(t);
    (next, t, coefficient, reversible)
}

#[test]
fn rational_young_shear_finite_heat_refreshes_entropy_and_capacity_while_retaining_memory() {
    with_cx(|cx, gate| {
        let (claims, properties) = rational_claims(RATIONAL_CP, true);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let initial = body
            .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
            .unwrap();
        let accepted_state = initial.material_history().state().to_vec();
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let mut thermal = heat(&[&curve], &[RATIONAL_CP], 0.0);
        // The compatibility entry point must keep its narrower promise.
        assert!(
            initial
                .advance_affine_thermal_candidate(
                    &thermal,
                    0.005,
                    &FORCE,
                    &cells,
                    controls(),
                    budget(),
                    &[0.01],
                    no_extra,
                    gate,
                )
                .is_err()
        );
        let mut model: Option<RelaxingSolid<'_>> = None;
        let mut y = [0.006, 0.08, 0.006, 0.0, 0.0];
        let mut temperature = 300.0;
        let mut total_coefficient = 0.0;
        for target in [0.005, 0.01] {
            let history = model
                .as_ref()
                .map(RelaxingSolid::history)
                .unwrap_or_else(|| initial.clone());
            let (next, t, coefficient, reversible) = rational_endpoint(y, temperature, 0.005);
            let result = history
                .advance_material_thermal_candidate(
                    &thermal,
                    target,
                    &FORCE,
                    &cells,
                    controls(),
                    budget(),
                    &[0.01],
                    |_| {
                        Ok::<_, Infallible>(StepCoefficients {
                            sources_w: BTreeMap::from([(CellId(0), RATIONAL_HEATER_W)]),
                            conductances_w_k: BTreeMap::new(),
                            payload: (),
                        })
                    },
                    gate,
                )
                .unwrap();
            let actual = result.model().history();
            let native = actual.material_history();
            let report = result.report();
            let exchange = report.thermodynamics.as_ref().unwrap();
            let expected = rational_state(next, t);
            near(result.thermal().temperature_k(CellId(0)).unwrap(), t, 2e-7);
            near(native.displacement_m()[0], next[0], 2e-9);
            near(native.velocity_m_s()[0], next[1], 2e-8);
            near(native.dissipated_heat_j(), next[3], 2e-11);
            near(native.supplied_work_j(), next[4], 2e-11);
            near(
                exchange.before[0].entropy_j_k,
                rational_state(y, temperature)[1],
                2e-13,
            );
            near(exchange.after[0].free_energy_j, expected[0], 2e-11);
            near(exchange.after[0].entropy_j_k, expected[1], 2e-13);
            near(exchange.after[0].heat_capacity_j_k, expected[2], 2e-13);
            assert!(exchange.after[0].heat_capacity_j_k < 0.0);
            near(
                exchange.after[0].internal_energy_j,
                expected[0] + t * expected[1],
                2e-11,
            );
            near(exchange.coefficient_work_j[0], coefficient, 2e-11);
            near(exchange.reversible_heat_j[0], reversible, 3e-11);
            total_coefficient += coefficient;
            near(
                native.thermal_coefficient_work_j(),
                total_coefficient,
                3e-11,
            );
            assert_eq!(native.parameter_work_j(), 0.0);
            let physical = 0.5 * (10.0 / 60.0) * (next[1].powi(2) - y[1].powi(2))
                + exchange.after[0].internal_energy_j
                - exchange.before[0].internal_energy_j
                + result.thermal().energy_j()
                - thermal.energy_j()
                - (next[4] - y[4])
                - RATIONAL_HEATER_W * 0.005;
            near(physical, 0.0, 3e-10);
            near(report.total_energy_residual_j, 0.0, 2e-10);
            assert_eq!(
                result.model().thermodynamics().unwrap(),
                exchange.after.as_slice()
            );
            let (next_model, next_thermal, _) = result.into_parts();
            model = Some(next_model);
            thermal = next_thermal;
            temperature = t;
            y = next;
        }
        assert_eq!(
            model.as_ref().unwrap().history().material_history().steps(),
            2
        );
        assert_eq!(initial.material_history().state(), accepted_state);
        assert_eq!(initial.material_history().thermal_coefficient_work_j(), 0.0);
    });
}

#[test]
fn rational_held_volumetric_strain_closes_finite_heat_against_exact_internal_energy() {
    with_cx(|cx, gate| {
        let (claims, properties) = rational_claims(RATIONAL_CP, false);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = ViscoElements::new(
            SolidMeshRef::Tet4 {
                nodes_m: &NODES,
                elements: &[[0, 1, 2, 3]],
            },
            TetViscoField::Uniform(&TetViscoMaterial::from_source(&source, IDENTITY).unwrap()),
            &[0, 1, 2, 4, 5, 6, 8, 9, 10],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 128,
            },
            cx,
        )
        .unwrap();
        let body = TensorSolidBody::from_material_sources(
            elements
                .project_basis(&[vec![1.0, 1.0, 1.0]], controls().projection, cx)
                .unwrap(),
            vec![SolidMaterialSource {
                material: source,
                material_to_mesh: IDENTITY,
            }],
            dynamics(),
            gate,
        )
        .unwrap();
        let dilation: f64 = 0.004;
        let target_temperature = 310.0;
        let target_s = 0.005;
        let [bulk, slope, curvature] = rational_bulk(target_temperature);
        let old = rational_bulk(300.0);
        // Uniform epsilon_ii=x has psi=(9/2)*V*K*x^2. At the desired
        // temperature the fixed load balances every free nodal direction,
        // yielding a manufactured held-strain endpoint of the dynamic solve.
        let scale = 0.75 * dilation.powi(2);
        let initial_u = scale * (old[0] - 300.0 * old[1]);
        let final_u = scale * (bulk - target_temperature * slope);
        let supplied =
            (10.0 / 6.0) * RATIONAL_CP * (target_temperature - 300.0) + final_u - initial_u;
        let force = 0.5 * bulk * dilation;
        let forces = [
            [0.0; 3],
            [force, 0.0, 0.0],
            [0.0, force, 0.0],
            [0.0, 0.0, force],
        ];
        let history = body
            .history(&[dilation], &[0.0], InitialMemory::Relaxed, gate)
            .unwrap();
        let thermal = heat(&[&curve], &[RATIONAL_CP], 0.0);
        let context = query();
        let result = history
            .advance_material_thermal_candidate(
                &thermal,
                target_s,
                &forces,
                &[cell(0, &card, &curve, &context)],
                controls(),
                budget(),
                &[0.012],
                |_| {
                    Ok::<_, Infallible>(StepCoefficients {
                        sources_w: BTreeMap::from([(CellId(0), supplied / target_s)]),
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    })
                },
                gate,
            )
            .unwrap();
        let accepted = result.model().history();
        let exchange = result.report().thermodynamics.as_ref().unwrap();
        near(
            result.thermal().temperature_k(CellId(0)).unwrap(),
            target_temperature,
            2e-7,
        );
        near(
            accepted.material_history().displacement_m()[0],
            dilation,
            2e-10,
        );
        near(accepted.material_history().velocity_m_s()[0], 0.0, 2e-9);
        near(exchange.after[0].free_energy_j, scale * bulk, 2e-11);
        near(exchange.after[0].entropy_j_k, -scale * slope, 2e-13);
        near(
            exchange.after[0].heat_capacity_j_k,
            -target_temperature * scale * curvature,
            2e-13,
        );
        near(exchange.after[0].internal_energy_j, final_u, 2e-11);
        near(
            result.thermal().energy_j() - thermal.energy_j() + final_u - initial_u,
            supplied,
            3e-10,
        );
        assert_eq!(accepted.material_history().dissipated_heat_j(), 0.0);
        near(accepted.material_history().supplied_work_j(), 0.0, 2e-11);
        near(result.report().total_energy_residual_j, 0.0, 2e-10);
    });
}

#[test]
fn rational_capacity_refusal_covers_the_whole_strain_domain_and_preserves_inputs() {
    with_cx(|cx, gate| {
        let (claims, properties) = rational_claims(CP, false);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let derivative = TensorThermalLaw::try_new(&source, &card, "T", [280.0, 330.0])
            .unwrap()
            .in_frame(IDENTITY)
            .unwrap()
            .at_temperature(300.0)
            .unwrap();
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let bound = body.projection().material_points()[0]
            .thermal_capacity_bound_j_k2(&derivative, 300.0, [0.0; 6])
            .unwrap();
        let caloric_capacity = (10.0 / 6.0) * CP;
        assert!(
            caloric_capacity - 300.0 * rational_axial(300.0)[2] * 0.006_f64.powi(2) / 12.0 > 0.0
        );
        assert!(caloric_capacity < 330.0 * bound);
        let history = body
            .history(&[0.006], &[0.08], InitialMemory::Relaxed, gate)
            .unwrap();
        let before = history.material_history().state().to_vec();
        let thermal = heat(&[&curve], &[CP], 0.0);
        let energy = thermal.energy_j();
        let context = query();
        let calls = std::cell::Cell::new(0);
        let failed = history.advance_material_thermal_candidate(
            &thermal,
            0.005,
            &FORCE,
            &[cell(0, &card, &curve, &context)],
            controls(),
            budget(),
            &[0.00005],
            |t| {
                calls.set(calls.get() + 1);
                no_extra(t)
            },
            gate,
        );
        assert!(matches!(failed,
            Err(fs_couple::circuit_enthalpy::CircuitEnthalpyError::Invalid(message))
                if message.contains("positive total fixed-strain heat capacity")));
        assert_eq!(calls.get(), 0);
        assert_eq!(history.material_history().state(), before);
        assert_eq!(history.material_history().steps(), 0);
        assert_eq!(history.material_history().thermal_coefficient_work_j(), 0.0);
        assert_eq!(history.material_history().parameter_work_j(), 0.0);
        assert_eq!(thermal.energy_j(), energy);
    });
}
