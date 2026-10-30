use super::*;

fn variable_cp_problem(beta: f64) -> Problem {
    let mut p = problem(1.0, 2.0, 1.0, 1.0, true);
    p.materials[0] = material_source(1.0, 2.0, 1.0, true, "liquid", |claim| {
        if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: [290.0, 300.0, 320.0, 350.0]
                    .map(|t| (t, 1.0 + beta * (t - 300.0)))
                    .to_vec(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
    });
    p.fluid = load_fluid(
        &input(1.0),
        &mut p.domain,
        &p.material_indices,
        &p.materials,
    )
    .unwrap();
    p
}

// Independent nonlinear three-control-volume equations with analytic
// h(theta) = theta + beta/2 * theta^2. This does not call the native chart,
// finite-volume assembly, material resolver, or nonlinear iteration.
fn analytic_step(old: [f64; 3], beta: f64, power: f64) -> ([f64; 3], f64, f64) {
    let h = |t: f64| t + 0.5 * beta * t * t;
    let mut t = old;
    let mut converged = false;
    for _ in 0..30 {
        let c = t.map(|v| 1.0 + beta * v);
        let enthalpy = t.map(h);
        let r = [
            10.0 * (enthalpy[0] - h(old[0])) + 3.0 * t[0] - t[1] + enthalpy[0],
            10.0 * (enthalpy[1] - h(old[1])) + 2.0 * t[1] - t[0] - t[2] + enthalpy[1]
                - enthalpy[0]
                - power,
            10.0 * (enthalpy[2] - h(old[2])) + t[2] - t[1] + enthalpy[2] - enthalpy[1],
        ];
        if r.iter().all(|v| v.abs() < 1e-13) {
            converged = true;
            break;
        }
        let mut d = [11.0 * c[0] + 3.0, 11.0 * c[1] + 2.0, 11.0 * c[2] + 1.0];
        let mut rhs = r.map(|v| -v);
        for i in 1..3 {
            let factor = (-1.0 - c[i - 1]) / d[i - 1];
            d[i] += factor;
            rhs[i] -= factor * rhs[i - 1];
        }
        let mut delta = [0.0; 3];
        delta[2] = rhs[2] / d[2];
        delta[1] = (rhs[1] + delta[2]) / d[1];
        delta[0] = (rhs[0] + delta[1]) / d[0];
        for i in 0..3 {
            t[i] += delta[i];
        }
    }
    assert!(converged);
    let stored = (0..3).map(|i| h(t[i]) - h(old[i])).sum::<f64>();
    let outflow = 0.1 * (2.0 * t[0] + h(t[2]));
    (t, stored, outflow)
}

#[test]
fn g1_g3_variable_source_cp_changes_storage_and_flow_through_drive_switch_off() {
    let mut peaks = Vec::new();
    for beta in [0.0, 0.08] {
        let p = variable_cp_problem(beta);
        let bytes = p.materials[0].pack.to_bytes();
        let flow = p.fluid.as_ref().unwrap().flow.clone();
        let r = request(request_text(), 3);
        let trajectory = march(&p, &r, &CancelGate::new()).unwrap();
        let output = render(&p, &trajectory, true);
        let value = J::parse(&output).unwrap();
        assert_eq!(
            value.str_field("fluid_advection"),
            Some("source-enthalpy-first-order-upwind")
        );
        let mut theta = [0.0; 3];
        for (i, row) in value
            .get("records")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
        {
            let power = if i < 4 { 4.5 } else { 0.0 };
            let (next, stored, outflow) = analytic_step(theta, beta, power);
            assert!((stored + outflow - 0.1 * power).abs() < 1e-13);
            // Source integration's bounded piecewise-linear chart is distinct
            // from these analytic Cp integrals; do not call it exact physics.
            assert!((row.f64_field("stored_energy_change_j").unwrap() - stored).abs() < 8e-5);
            assert!((row.f64_field("physical_outflow_j").unwrap() - outflow).abs() < 8e-5);
            assert!(row.f64_field("closure_j").unwrap().abs() < 2e-7);
            theta = next;
        }
        let actual = temperatures(&value);
        for i in 0..3 {
            assert!((actual[i] - 300.0 - theta[i]).abs() < 8e-5);
        }
        assert_eq!(p.materials[0].pack.to_bytes(), bytes);
        assert_eq!(p.fluid.as_ref().unwrap().flow, flow);
        assert_eq!(p.initial_temperature, vec![300.0; 3]);
        assert_eq!(
            output,
            render(&p, &march(&p, &r, &CancelGate::new()).unwrap(), true)
        );
        peaks.push(actual[1]);
    }
    assert!(peaks[1] < peaks[0] - 0.001);
}

#[test]
fn g0_variable_caloric_source_must_cover_each_inlet_and_accepted_cell_state() {
    let mut p = variable_cp_problem(0.08);
    let before = p.materials[0].pack.to_bytes();
    let tight = request_text().replace(
        "\"temperature_max_k\":350",
        "\"temperature_max_k\":300.00001",
    );
    assert!(march(&p, &request(&tight, 3), &CancelGate::new()).is_err());
    p.thermal.faces[1] = ThermalFace::Outflow {
        backflow_temperature: 351.0,
    };
    assert!(march(&p, &request(request_text(), 3), &CancelGate::new()).is_err());
    p.thermal.faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    assert!(march(&p, &request(request_text(), 3), &CancelGate::new()).is_ok());
    assert_eq!(p.materials[0].pack.to_bytes(), before);
    assert_eq!(p.initial_temperature, vec![300.0; 3]);
}
