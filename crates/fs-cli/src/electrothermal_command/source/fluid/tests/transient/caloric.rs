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

#[test]
fn g1_g3_steady_source_cp_uses_actual_enthalpy_flux_and_preserves_reference_flow() {
    let mut peaks = Vec::new();
    for beta in [0.0, 0.08] {
        let mut p = variable_cp_problem(beta);
        p.thermal.power_w = vec![0.0, 55.5, 0.0]; // Plus 4.5 W actual Joule heat.
        let bytes = p.materials[0].pack.to_bytes();
        let reference = p.fluid.as_ref().unwrap().reference.clone();
        let flow = p.fluid.as_ref().unwrap().flow.clone();
        assert_eq!(
            p.fluid.as_ref().unwrap().properties.specific_heat_j_kg_k,
            1.0
        );
        let result = p.solve(&CancelGate::new()).unwrap();
        // Independent steady balances with theta=T-300 and h=theta+beta/2*theta^2:
        // b=c=3*a+h(a), and total heat 2*a+h(b)=60 W. This scalar
        // relation uses neither source-chart inversion nor native assembly.
        let h = |theta: f64| theta + 0.5 * beta * theta * theta;
        let (mut lo, mut hi) = (0.0, 30.0);
        for _ in 0..64 {
            let a = 0.5 * lo + 0.5 * hi;
            if 2.0 * a + h(3.0 * a + h(a)) < 60.0 {
                lo = a;
            } else {
                hi = a;
            }
        }
        let a = 0.5 * lo + 0.5 * hi;
        let b = 3.0 * a + h(a);
        let expected = [300.0 + a, 300.0 + b, 300.0 + b];
        let energy = &result.coupled.thermal;
        for (actual, expected) in energy.temperature.iter().zip(expected) {
            assert!((actual - expected).abs() < 2e-8, "{actual} != {expected}");
        }
        assert_eq!(energy.report.scheme, ConvectionScheme::Upwind);
        assert!((result.coupled.electric.total_joule_w - 4.5).abs() < 1e-10);
        assert!((energy.report.balance.source_w - 60.0).abs() < 1e-10);
        assert!((energy.report.balance.advective_outflow_w - h(b)).abs() < 2e-8);
        assert!((energy.report.balance.boundary_outflow_w - 60.0).abs() < 2e-8);
        assert!(energy.report.balance.residual_w.abs() < 2e-8);
        assert_eq!(p.materials[0].pack.to_bytes(), bytes);
        assert_eq!(p.fluid.as_ref().unwrap().reference, reference);
        assert_eq!(p.fluid.as_ref().unwrap().flow, flow);
        assert_eq!(p.initial_temperature, vec![300.0; 3]);
        assert_eq!(p.solve(&CancelGate::new()).unwrap().coupled, result.coupled);
        peaks.push(energy.temperature[1]);
    }
    // Both sources have the SAME startup Cp; only the original Cp(T) differs.
    assert!(peaks[1] < peaks[0] - 10.0);
}

#[test]
fn g1_source_cp_minmod_closes_manufactured_steady_and_transient_cell_balances() {
    // Independent five-cell balance, with unit density, volume, k and flow.
    // h(theta)=theta+0.04*theta^2. The monotone profile has nonzero, unequal
    // limited face corrections; an upwind-only implementation cannot pass.
    let theta: [f64; 5] = [1.0, 2.0, 4.0, 7.0, 8.0];
    let h = theta.map(|t| t + 0.04 * t * t);
    let face_h = [h[0], h[1] + 0.56, h[2] + 1.24, h[3] + 0.8, h[4]];
    let mut steady_power = [0.0; 5];
    for cell in 0..5 {
        let conduction = if cell == 0 {
            2.0 * theta[0] + theta[0] - theta[1]
        } else if cell == 4 {
            theta[4] - theta[3]
        } else {
            2.0 * theta[cell] - theta[cell - 1] - theta[cell + 1]
        };
        let inflow = if cell == 0 { 0.0 } else { face_h[cell - 1] };
        steady_power[cell] = conduction + face_h[cell] - inflow;
    }
    assert!((steady_power.iter().sum::<f64>() - 12.56).abs() < 1e-13);

    for transient in [false, true] {
        let mut p = variable_cp_problem(0.08);
        p.domain = VoxelDomain::from_fn(5, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        p.material_indices = vec![0; 5];
        p.initial_temperature = vec![300.0; 5];
        p.fluid = load_fluid(&input(1.0), &mut p.domain, &p.material_indices, &p.materials).unwrap();
        p.electric.electrodes[1].lo = [4, 0, 0];
        p.electric.electrodes[1].hi = [5, 1, 1];
        p.electric.electrodes[1].drive = ElectrodeDrive::Voltage(0.0);
        p.energy.scheme = ConvectionScheme::EnthalpyMinmod;
        p.energy.tolerance = 1e-11;
        p.thermal.power_w = (0..5).map(|cell| {
            steady_power[cell] + if transient { 10.0 * h[cell] } else { 0.0 }
        }).collect();
        let original = p.materials[0].pack.to_bytes();
        let flow = p.fluid.as_ref().unwrap().flow.clone();
        if transient {
            let text = request_text().replace("\"steps\":8", "\"steps\":1")
                .replace("[[1,1,1],[5,0,1]]", "[[1,1,1]]");
            let r = request(&text, 5);
            let result = march(&p, &r, &CancelGate::new()).unwrap();
            let output = render(&p, &result, true);
            let value = J::parse(&output).unwrap();
            assert_eq!(value.str_field("fluid_advection"), Some("source-enthalpy-minmod"));
            for (actual, expected) in temperatures(&value).into_iter().zip(theta) {
                assert!((actual - 300.0 - expected).abs() < 2e-7);
            }
            let rows = value.get("records").unwrap().as_array().unwrap();
            assert_eq!(rows.len(), 1);
            assert!((rows[0].f64_field("stored_energy_change_j").unwrap() - 27.36).abs() < 2e-7);
            assert!((rows[0].f64_field("physical_outflow_j").unwrap() - 1.256).abs() < 2e-7);
            assert!((rows[0].f64_field("source_j").unwrap() - 28.616).abs() < 2e-7);
            assert!(rows[0].f64_field("closure_j").unwrap().abs() < 2e-7);
            assert_eq!(output, render(&p, &march(&p, &r, &CancelGate::new()).unwrap(), true));
        } else {
            let result = p.solve(&CancelGate::new()).unwrap();
            let energy = &result.coupled.thermal;
            assert_eq!(energy.report.scheme, ConvectionScheme::EnthalpyMinmod);
            for (actual, expected) in energy.temperature.iter().zip(theta) {
                assert!((actual - 300.0 - expected).abs() < 2e-7);
            }
            assert_eq!(result.coupled.electric.total_joule_w, 0.0);
            assert!((energy.report.balance.advective_outflow_w - 10.56).abs() < 2e-7);
            assert!((energy.report.balance.boundary_outflow_w - 12.56).abs() < 2e-7);
            assert!(energy.report.balance.residual_w.abs() < 2e-7);
            p.energy.scheme = ConvectionScheme::Upwind;
            let upwind = p.solve(&CancelGate::new()).unwrap();
            let difference = upwind.coupled.thermal.temperature.iter().zip(theta)
                .map(|(actual, expected)| (actual - 300.0 - expected).abs())
                .fold(0.0, f64::max);
            assert!(difference > 0.05);
        }
        assert_eq!(p.materials[0].pack.to_bytes(), original);
        assert_eq!(p.fluid.as_ref().unwrap().flow, flow);
        assert_eq!(p.initial_temperature, vec![300.0; 5]);
    }
}

#[test]
fn g0_steady_source_cp_requires_complete_support_and_constant_density() {
    for refusal in ["unbounded", "discrete", "density", "narrow"] {
        let mut p = variable_cp_problem(0.08);
        p.materials[0] = material_source(1.0, 2.0, 1.0, true, "liquid", |claim| {
            if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
                match refusal {
                    "unbounded" => claim.validity = ValidityDomain::unconstrained(),
                    "discrete" => {
                        claim.value = PropertyValue::Curve {
                            abscissa: "T".into(),
                            abscissa_dims: fs_qty::Temperature::DIMS,
                            knots: vec![(290.0, 0.5), (300.0, 1.0), (350.0, 5.0)],
                            dims: SPECIFIC_HEAT_CAPACITY_DIMS,
                        };
                        claim.interpolation = InterpolationPolicy::TabulatedOnly;
                    }
                    "narrow" => {
                        claim.validity = ValidityDomain::unconstrained().with("T", 299.0, 301.0)
                    }
                    _ => {}
                }
            }
            if refusal == "density" && claim.key.name() == DENSITY_PROPERTY {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(290.0, 1.0), (300.0, 1.0), (320.0, 0.9), (350.0, 1.0)],
                    dims: Density::DIMS,
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
        let original = p.materials[0].pack.to_bytes();
        assert!(
            p.solve(&CancelGate::new()).is_err(),
            "{refusal} must refuse"
        );
        assert_eq!(p.materials[0].pack.to_bytes(), original);
        assert_eq!(p.initial_temperature, vec![300.0; 3]);
    }
    let mut p = variable_cp_problem(0.08);
    p.thermal.faces[1] = ThermalFace::Outflow {
        backflow_temperature: 351.0,
    };
    assert!(p.solve(&CancelGate::new()).is_err());
    p.thermal.faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(p.solve(&gate), Err(ChtError::Cancelled)));
    assert!(p.solve(&CancelGate::new()).is_ok());
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
fn g1_transient_caloric_endpoint_requires_local_heat_with_loose_iteration_controls() {
    let mut p = variable_cp_problem(0.08);
    // Loose iteration-change and integrated-energy allowances must not admit
    // a stale Cp tangent that fails the requested local native heat tolerance.
    p.coupling.temperature_tolerance_k = 100.0;
    p.energy.tolerance = 1e-11;
    let text = request_text()
        .replace("\"energy_absolute_tolerance_j\":0.0000001", "\"energy_absolute_tolerance_j\":1000")
        .replace("\"energy_relative_tolerance\":0.000000001", "\"energy_relative_tolerance\":0");
    let r = request(&text, 3);
    let bytes = p.materials[0].pack.to_bytes();
    let trajectory = march(&p, &r, &CancelGate::new()).unwrap();
    let value = J::parse(&render(&p, &trajectory, true)).unwrap();
    let mut theta = [0.0; 3];
    for (step, row) in value.get("records").unwrap().as_array().unwrap().iter().enumerate() {
        let power = if step < 4 { 4.5 } else { 0.0 };
        let (next, stored, outflow) = analytic_step(theta, 0.08, power);
        assert!((row.f64_field("stored_energy_change_j").unwrap() - stored).abs() < 2e-7);
        assert!((row.f64_field("physical_outflow_j").unwrap() - outflow).abs() < 2e-7);
        assert!(row.f64_field("closure_j").unwrap().abs() < 2e-7);
        theta = next;
    }
    for (actual, expected) in temperatures(&value).into_iter().zip(theta) {
        assert!((actual - 300.0 - expected).abs() < 2e-7);
    }
    assert_eq!(p.materials[0].pack.to_bytes(), bytes);
    assert_eq!(p.initial_temperature, vec![300.0; 3]);
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
