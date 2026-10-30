//! G1/G3/G4: actual source heating, reciprocal inertia and retained Maxwell strain.
mod affine;

use super::*;
use fs_blake3::ContentHash;
use fs_couple::tensor_material::solid::thermal::relaxation::{
    AffineThermalBudget, ThermalExpansionCell,
};
use fs_solid::linear3::source::ThermalExpansionModel;

const ALPHA: f64 = 1e-5;
const HEATER_W: f64 = 0.1;

fn claims(branch: bool, alpha: f64) -> (Vec<PropertyClaim>, TensorPronyProperties) {
    let (mut claims, properties) = fixture(1000.0, 10.0, if branch { 0.4 } else { 0.0 });
    claims.push(scalar(
        "specific_heat_capacity",
        Dims([2, 0, -2, -1, 0, 0]),
        CP,
    ));
    claims.push(scalar(
        "linear_thermal_expansion_coefficient",
        Dims([0, 0, 0, -1, 0, 0]),
        alpha,
    ));
    (claims, properties)
}
fn expansion(reference: &QueryPoint) -> ThermalExpansionCell<'_> {
    ThermalExpansionCell {
        reference,
        model: ThermalExpansionModel::Isotropic,
        selection: MaterialPropertySelection::SingleClaimOnly,
        mesh_frame: ContentHash([91; 32]),
    }
}
fn budget() -> AffineThermalBudget {
    AffineThermalBudget {
        max_absolute_coefficient_work_j: 0.1,
        relative_memory_tolerance: 1e-11,
    }
}
fn heating(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::from([(CellId(0), HEATER_W)]),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}

// Independent physical coordinates [u,v,q_x,q_y,q_z,dashpot heat,force work].
// x/y motion is constrained, but their thermal stress, memory and energy remain.
fn psi(y: [f64; 7], t: f64) -> f64 {
    let free = ALPHA * (t - 300.0);
    (1000.0 * (2.0 * free * free + (y[0] - free).powi(2))
        + 400.0 * y[2] * y[2]
        + 40.0 * y[3] * y[3]
        + 40.0 * y[2] * y[3]
        + 280.0 * y[4] * y[4])
        / 12.0
}
fn entropy(y: [f64; 7], t: f64) -> f64 {
    ALPHA
        * (1000.0 * (y[0] - 3.0 * ALPHA * (t - 300.0)) + 420.0 * y[2] + 60.0 * y[3] + 280.0 * y[4])
        / 6.0
}
fn trajectory(old: [f64; 7], t0: f64, t: f64, dt: f64) -> ([f64; 7], f64) {
    let mut y = old;
    for q in &mut y[2..5] {
        *q -= ALPHA * (t - t0);
    }
    let coefficient = psi(y, t) - psi(old, t0);
    let rhs = |z: [f64; 7]| {
        [
            z[1],
            (0.02 - (1000.0 * (z[0] - ALPHA * (t - 300.0)) + 280.0 * z[4]) / 6.0) / (10.0 / 60.0),
            -z[2] / 0.02,
            -z[3] / 0.02,
            z[1] - z[4] / 0.02,
            (400.0 * z[2] * z[2] + 40.0 * z[3] * z[3] + 40.0 * z[2] * z[3] + 280.0 * z[4] * z[4])
                / (6.0 * 0.02),
            0.02 * z[1],
        ]
    };
    let h = dt / 512.0;
    for _ in 0..512 {
        let a = rhs(y);
        let b = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| y[i] + h * c[i]));
        y = core::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    (y, coefficient)
}
fn endpoint(y: [f64; 7], t0: f64, dt: f64) -> ([f64; 7], f64, f64, f64) {
    let mut lo = 280.0;
    let mut hi = 330.0;
    for _ in 0..48 {
        let t = 0.5 * (lo + hi);
        let (next, coefficient) = trajectory(y, t0, t, dt);
        let reversible = -coefficient - t * entropy(next, t) + t0 * entropy(y, t0);
        let residual = (10.0 / 6.0) * CP * (t - t0) - (next[5] - y[5] + reversible + HEATER_W * dt);
        if residual > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = 0.5 * (lo + hi);
    let (next, coefficient) = trajectory(y, t0, t, dt);
    (
        next,
        t,
        coefficient,
        -coefficient - t * entropy(next, t) + t0 * entropy(y, t0),
    )
}

#[test]
fn thermal_expansion_preserves_native_memory_and_closes_physical_energy_on_continuing_motion() {
    with_cx(|cx, gate| {
        let (claims, properties) = claims(true, ALPHA);
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
        let expansion = [expansion(&context)];
        let mut thermal = heat(&[&curve], &[CP], 0.0);
        let initial_caloric = thermal.energy_j();
        let mut model: Option<RelaxingSolid<'_>> = None;
        let mut y = [0.006, 0.08, 0.0, 0.0, 0.006, 0.0, 0.0];
        let mut temperature = 300.0;
        let mut coefficient = 0.0;
        let initial_physical =
            0.5 * (10.0 / 60.0) * y[1] * y[1] + psi(y, 300.0) + 300.0 * entropy(y, 300.0);
        for target in [0.005, 0.01] {
            let history = model
                .as_ref()
                .map(RelaxingSolid::history)
                .unwrap_or_else(|| initial.clone());
            let (next, t, coeff, rev) = endpoint(y, temperature, 0.005);
            let result = history
                .advance_expanding_thermal_candidate(
                    &thermal,
                    target,
                    &FORCE,
                    &cells,
                    &expansion,
                    controls(),
                    budget(),
                    &[0.0005],
                    heating,
                    gate,
                )
                .unwrap();
            let report = result.report();
            let accepted = result.model().history();
            let native = accepted.material_history();
            near(native.displacement_m()[0], next[0], 2e-10);
            near(native.velocity_m_s()[0], next[1], 2e-9);
            near(native.dissipated_heat_j(), next[5], 1e-10);
            near(native.supplied_work_j(), next[6], 1e-11);
            let actual_t = curve
                .curve()
                .state_at_specific_enthalpy(result.thermal().cells()[&CellId(0)].specific_h_j_kg)
                .unwrap()
                .temperature_k();
            near(actual_t, t, 2e-6);
            assert!(actual_t > temperature);
            coefficient += coeff;
            near(native.thermal_coefficient_work_j(), coefficient, 1e-10);
            assert_eq!(native.parameter_work_j(), 0.0);
            near(*native.state().last().unwrap(), 1.0, 2e-13);
            assert_eq!(native.steps(), if target == 0.005 { 1 } else { 2 });
            let exchange = report.thermodynamics.as_ref().unwrap();
            near(exchange.reversible_heat_j[0], rev, 1e-10);
            near(exchange.after[0].free_energy_j, psi(next, t), 1e-10);
            near(exchange.after[0].entropy_j_k, entropy(next, t), 1e-12);
            let point = &result.model().body().projection().material_points()[0];
            let columns: Vec<_> = (0..6)
                .map(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
                .collect();
            let factors = point
                .material()
                .project_strain(&columns, 1.0 / 6.0, projection(), || false)
                .unwrap();
            let q = [next[2], next[3], next[4], 0.0, 0.0, 0.0];
            for (row, z) in factors.branches()[0]
                .rows()
                .iter()
                .zip(&native.state()[2..8])
            {
                near(*z, row.iter().zip(q).map(|(a, b)| a * b).sum(), 2e-10);
            }
            near(native.balance_residual_j(), 0.0, 3e-10);
            let physical = 0.5 * (10.0 / 60.0) * native.velocity_m_s()[0].powi(2)
                + exchange.after[0].internal_energy_j;
            near(
                physical - initial_physical + result.thermal().energy_j()
                    - initial_caloric
                    - native.supplied_work_j()
                    - HEATER_W * target,
                0.0,
                5e-10,
            );
            near(report.total_energy_residual_j, 0.0, 3e-10);
            let (next_model, next_thermal, _) = result.into_parts();
            model = Some(next_model);
            thermal = next_thermal;
            y = next;
            temperature = t;
        }
        assert_eq!(initial.material_history().state(), original);
    });
}

#[test]
fn rotated_orthotropic_expansion_drives_real_motion_and_preserves_work_conjugate_entropy() {
    with_cx(|cx, gate| {
        use fs_material::state_point::ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES;
        let (mut claims, properties) = fixture(1000.0, 10.0, 0.0);
        claims.push(scalar(
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            CP,
        ));
        for (name, value) in ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip([ALPHA, 2.0 * ALPHA, 3.0 * ALPHA])
        {
            claims.push(scalar(name, Dims([0, 0, 0, -1, 0, 0]), value));
        }
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let history = body
            .history(&[0.0], &[0.0], InitialMemory::Relaxed, gate)
            .unwrap();
        let thermal = heat(&[&curve], &[CP], 0.0);
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let angle = core::f64::consts::FRAC_1_SQRT_2;
        let rotation = [[angle, -angle, 0.0], [angle, angle, 0.0], [0.0, 0.0, 1.0]];
        let mut endpoints = Vec::new();
        for frame in [IDENTITY, rotation] {
            let bindings = [ThermalExpansionCell {
                model: ThermalExpansionModel::OrientedOrthotropic {
                    principal_to_mesh: frame,
                    principal_frame: ContentHash([92; 32]),
                },
                ..expansion(&context)
            }];
            let result = history
                .advance_expanding_thermal_candidate(
                    &thermal,
                    0.005,
                    &[[0.0; 3]; 4],
                    &cells,
                    &bindings,
                    controls(),
                    budget(),
                    &[0.0005],
                    heating,
                    gate,
                )
                .unwrap();
            let t = curve
                .curve()
                .state_at_specific_enthalpy(result.thermal().cells()[&CellId(0)].specific_h_j_kg)
                .unwrap()
                .temperature_k();
            let accepted = result.model().history();
            let u = accepted.material_history().displacement_m()[0];
            assert!(u > 0.0);
            assert_eq!(accepted.material_history().supplied_work_j(), 0.0);
            let transform = fs_material::tensor::mandel_rotation(frame);
            let alpha: [f64; 6] = core::array::from_fn(|i| {
                transform[i][0] * ALPHA
                    + transform[i][1] * 2.0 * ALPHA
                    + transform[i][2] * 3.0 * ALPHA
            });
            let point = &result.model().body().projection().material_points()[0];
            for (free, a) in point.free_strain_mandel().iter().zip(alpha) {
                near(*free, a * (t - 300.0), 1e-14);
            }
            let elastic: [f64; 6] = core::array::from_fn(|i| {
                if i == 2 {
                    u - alpha[i] * (t - 300.0)
                } else {
                    -alpha[i] * (t - 300.0)
                }
            });
            let s = 1000.0 / 6.0 * alpha.iter().zip(elastic).map(|(a, e)| a * e).sum::<f64>();
            near(
                result.report().thermodynamics.as_ref().unwrap().after[0].entropy_j_k,
                s,
                1e-12,
            );
            near(result.report().total_energy_residual_j, 0.0, 3e-10);
            endpoints.push((t, u));
        }
        near(endpoints[0].0, endpoints[1].0, 2e-8);
        near(endpoints[0].1, endpoints[1].1, 1e-12);
    });
}

#[test]
fn zero_expansion_retains_the_original_trajectory_and_state_extent() {
    with_cx(|cx, gate| {
        let (claims, properties) = claims(true, 0.0);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let history = body
            .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
            .unwrap();
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let thermal = heat(&[&curve], &[CP], 0.0);
        let standard = history
            .advance_affine_thermal_candidate(
                &thermal,
                0.005,
                &FORCE,
                &cells,
                controls(),
                budget(),
                &[0.0005],
                heating,
                gate,
            )
            .unwrap();
        let expanded = history
            .advance_expanding_thermal_candidate(
                &thermal,
                0.005,
                &FORCE,
                &cells,
                &[expansion(&context)],
                controls(),
                budget(),
                &[0.0005],
                heating,
                gate,
            )
            .unwrap();
        assert_eq!(
            standard.model().history().material_history().state(),
            expanded.model().history().material_history().state()
        );
        assert_eq!(standard.thermal().energy_j(), expanded.thermal().energy_j());
        assert_eq!(
            standard.report().mechanical.element_heat_j,
            expanded.report().mechanical.element_heat_j
        );
    });
}

#[test]
fn unsupported_expansion_thermodynamics_and_late_cancellation_publish_nothing() {
    with_cx(|cx, gate| {
        for invalid_case in ["alpha", "stiffness", "capacity"] {
            let (mut claims, properties) = claims(
                true,
                if invalid_case == "capacity" {
                    1e-4
                } else {
                    ALPHA
                },
            );
            if invalid_case != "capacity" {
                let name = if invalid_case == "alpha" {
                    "linear_thermal_expansion_coefficient"
                } else {
                    "young_modulus"
                };
                let claim = claims.iter_mut().find(|c| c.key.name() == name).unwrap();
                let (v, dims) = match claim.value {
                    PropertyValue::Scalar { value, dims } => (value, dims),
                    _ => unreachable!(),
                };
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(280.0, v * 0.9), (300.0, v),
                        (330.0, if invalid_case == "stiffness" { v * 1.16 } else { v * 1.15 })],
                    dims,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
            let card = card(claims);
            let curve = caloric(&card);
            let source = bind(&card, properties);
            let elements = elements(&source, cx);
            let body = body(&elements, source, cx, gate);
            let history = body
                .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
                .unwrap();
            let before = history.material_history().state().to_vec();
            let thermal = heat(&[&curve], &[CP], 0.0);
            let energy = thermal.energy_j();
            let context = query();
            assert!(
                history
                    .advance_expanding_thermal_candidate(
                        &thermal,
                        0.005,
                        &FORCE,
                        &[cell(0, &card, &curve, &context)],
                        &[expansion(&context)],
                        controls(),
                        budget(),
                        &[0.0005],
                        heating,
                        gate
                    )
                    .is_err()
            );
            assert_eq!(history.material_history().state(), before);
            assert_eq!(thermal.energy_j(), energy);
        }
        let (claims, properties) = claims(true, ALPHA);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let history = body
            .history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate)
            .unwrap();
        let before = history.material_history().state().to_vec();
        let thermal = heat(&[&curve], &[CP], 0.0);
        let energy = thermal.energy_j();
        let context = query();
        let cells = [cell(0, &card, &curve, &context)];
        let bindings = [expansion(&context)];
        let late = CancelGate::new_clock_free();
        let failed = history.advance_expanding_thermal_candidate(
            &thermal,
            0.005,
            &FORCE,
            &cells,
            &bindings,
            controls(),
            budget(),
            &[0.0005],
            |t| {
                late.request();
                heating(t)
            },
            &late,
        );
        assert!(matches!(
            failed,
            Err(fs_couple::circuit_enthalpy::CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(history.material_history().state(), before);
        assert_eq!(thermal.energy_j(), energy);
        let result = history
            .advance_expanding_thermal_candidate(
                &thermal,
                0.005,
                &FORCE,
                &cells,
                &bindings,
                controls(),
                budget(),
                &[0.0005],
                heating,
                gate,
            )
            .unwrap();
        let next = result.model().history();
        let other_reference = query().with("T", 299.0).unwrap();
        assert!(
            next.advance_expanding_thermal_candidate(
                result.thermal(),
                0.01,
                &FORCE,
                &cells,
                &[expansion(&other_reference)],
                controls(),
                budget(),
                &[0.0005],
                heating,
                gate
            )
            .is_err()
        );
        assert!(
            next.advance_affine_thermal_candidate(
                result.thermal(),
                0.01,
                &FORCE,
                &cells,
                controls(),
                budget(),
                &[0.0005],
                heating,
                gate
            )
            .is_err()
        );
    });
}

#[test]
fn unchanged_singular_branch_keeps_all_represented_viscous_strain_during_expansion() {
    use fs_material::visco_tensor::{PronyBranch, TensorProny};
    let mut equilibrium = [[0.0; 6]; 6];
    for (i, row) in equilibrium.iter_mut().enumerate() {
        row[i] = 1000.0;
    }
    let mut factor = [[0.0; 6]; 6];
    factor[0][0] = 10.0;
    factor[0][2] = 20.0;
    let law = TensorProny::new_mandel(
        equilibrium,
        vec![PronyBranch::from_gram_rows(factor, 0.02).unwrap()],
        0.1,
    )
    .unwrap();
    let point = law
        .project_strain(
            &[[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]],
            0.25,
            projection(),
            || false,
        )
        .unwrap();
    let next = point
        .with_eigenstrain_mandel([0.001, 0.002, 0.003, 0.0, 0.0, 0.0])
        .unwrap();
    let memory = [0.07, 0.0, 0.0, 0.0, 0.0, 0.0];
    let changed = point
        .transfer_factor_memory(&next, &memory, 1e-11, || false)
        .unwrap();
    near(
        changed[0],
        0.07 - 0.5 * (10.0 * 0.001 + 20.0 * 0.003),
        1e-15,
    );
    let returned = next
        .transfer_factor_memory(&point, &changed, 1e-11, || false)
        .unwrap();
    for (a, b) in returned.iter().zip(memory) {
        near(*a, b, 1e-15);
    }
}
