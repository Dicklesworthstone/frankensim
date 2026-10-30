//! Combined source C(T)/alpha continuation and its conservative stability domain.
use super::*;
use fs_material::visco_tensor::source::thermal::AffineTensorThermal;

const SLOPE: f64 = 0.0001;

fn affine_claims(cp: f64) -> (Vec<PropertyClaim>, TensorPronyProperties) {
    let (mut claims, properties) = super::claims(true, ALPHA);
    for claim in &mut claims {
        if claim.key.name() == "specific_heat_capacity" {
            claim.value = PropertyValue::Scalar {
                value: cp, dims: Dims([2, 0, -2, -1, 0, 0]),
            };
        }
        if let PropertyValue::Scalar { value, dims } = claim.value {
            if dims == Pressure::DIMS {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(280.0, value * factor(280.0)), (330.0, value * factor(330.0))],
                    dims,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
        }
    }
    (claims, properties)
}
fn factor(t: f64) -> f64 { 1.0 + SLOPE * (t - 300.0) }

// Independent physical-strain potential; no projection/factor state is used.
fn psi(y: [f64; 7], t: f64) -> f64 { factor(t) * super::psi(y, t) }
fn entropy(y: [f64; 7], t: f64) -> f64 {
    factor(t) * super::entropy(y, t) - SLOPE * super::psi(y, t)
}
fn trajectory(old: [f64; 7], t0: f64, t: f64, dt: f64) -> ([f64; 7], f64) {
    let mut y = old;
    for q in &mut y[2..5] { *q -= ALPHA * (t - t0); }
    let coefficient = psi(y, t) - psi(old, t0);
    let f = factor(t);
    let rhs = |z: [f64; 7]| [
        z[1],
        (0.02 - f * (1000.0 * (z[0] - ALPHA * (t - 300.0)) + 280.0 * z[4]) / 6.0) / (10.0 / 60.0),
        -z[2] / 0.02,
        -z[3] / 0.02,
        z[1] - z[4] / 0.02,
        f * (400.0 * z[2] * z[2] + 40.0 * z[3] * z[3] + 40.0 * z[2] * z[3] + 280.0 * z[4] * z[4]) / (6.0 * 0.02),
        0.02 * z[1],
    ];
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
    let (mut lo, mut hi) = (280.0, 330.0);
    for _ in 0..48 {
        let t = 0.5 * (lo + hi);
        let (next, coefficient) = trajectory(y, t0, t, dt);
        let reversible = -coefficient - t * entropy(next, t) + t0 * entropy(y, t0);
        let residual = (10.0 / 6.0) * CP * (t - t0) - (next[5] - y[5] + reversible + HEATER_W * dt);
        if residual > 0.0 { hi = t; } else { lo = t; }
    }
    let t = 0.5 * (lo + hi);
    let (next, coefficient) = trajectory(y, t0, t, dt);
    (next, t, coefficient, -coefficient - t * entropy(next, t) + t0 * entropy(y, t0))
}
#[test]
fn affine_stiffness_and_expansion_retain_memory_and_close_continuing_energy() {
    with_cx(|cx, gate| {
        let (claims, properties) = affine_claims(CP);
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
fn affine_expansion_checks_full_domain_capacity_before_advancing_either_history() {
    with_cx(|cx, gate| {
        let cp = 0.00004;
        let (claims, properties) = affine_claims(cp);
        let card = card(claims);
        let curve = caloric(&card);
        let source = bind(&card, properties);
        let derivative = AffineTensorThermal::try_new(&source, &card, "T", [280.0, 330.0])
            .unwrap().in_frame(IDENTITY).unwrap();
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let point = &body.projection().material_points()[0];
        let alpha = [ALPHA, ALPHA, ALPHA, 0.0, 0.0, 0.0];
        let constant = point.expansion_capacity_coefficient_j_k2(alpha).unwrap();
        let combined = point.affine_expansion_capacity_bound_j_k2(&derivative, 300.0, alpha).unwrap();
        assert!((10.0 / 6.0) * cp > 330.0 * constant);
        assert!((10.0 / 6.0) * cp < 330.0 * combined);
        let history = body.history(&[0.006], &[0.08], InitialMemory::Instantaneous, gate).unwrap();
        let before = history.material_history().state().to_vec();
        let thermal = heat(&[&curve], &[cp], 0.0);
        let energy = thermal.energy_j();
        let context = query();
        assert!(history.advance_expanding_thermal_candidate(
            &thermal, 0.005, &FORCE, &[cell(0, &card, &curve, &context)],
            &[expansion(&context)], controls(), budget(), &[0.0005], heating, gate,
        ).is_err());
        assert_eq!(history.material_history().state(), before);
        assert_eq!(thermal.energy_j(), energy);
    });
}

#[test]
fn affine_expansion_refuses_memory_outside_the_capacity_envelope() {
    with_cx(|cx, gate| {
        let (claims, properties) = affine_claims(CP);
        let card = card(claims);
        let source = bind(&card, properties);
        let derivative = AffineTensorThermal::try_new(&source, &card, "T", [280.0, 330.0])
            .unwrap().in_frame(IDENTITY).unwrap();
        let elements = elements(&source, cx);
        let body = body(&elements, source, cx, gate);
        let columns: Vec<_> = (0..6)
            .map(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 })).collect();
        let point = body.projection().material_points()[0].material()
            .project_strain(&columns, 1.0 / 6.0, projection(), || false).unwrap();
        // Total and equilibrium elastic strain are zero; only retained spring
        // memory exceeds 2*L=0.2. It must not evade the coupled capacity bound.
        let q = [0.3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let memory: Vec<f64> = point.branches()[0].rows().iter()
            .map(|row| row.iter().zip(q).map(|(a, b)| a * b).sum()).collect();
        assert!(point.thermal_state_with_expansion(
            &derivative, &[0.0; 6], &memory, 300.0,
            [ALPHA, ALPHA, ALPHA, 0.0, 0.0, 0.0], 1e-11, || false,
        ).is_err());
        // Zero alpha retains the previously admitted affine-C thermodynamics.
        assert!(point.thermal_state(
            &derivative, &[0.0; 6], &memory, 300.0, 1e-11, || false,
        ).is_ok());
    });
}
