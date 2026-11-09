//! Frozen signed thermal prestress reaches the original structural trajectory.
use super::*;

#[test]
fn native_and_hanging_prestress_preserve_reference_balance_mass_and_transposes() {
    let scale = 0.37;
    let mut cart = cartesian(0.1);
    cart.set_scales(&[scale]).unwrap();
    let mut tree = adaptive();
    tree.set_scales(&vec![scale; tree.cells()]).unwrap();
    assert!(tree.physical_nodes().len() > tree.nodes().len());
    for source in [
        CutElasticityRef3::Cartesian(&cart),
        CutElasticityRef3::Adaptive(&tree),
    ] {
        let (mass, prepared, nodes, native_n) = match source {
            CutElasticityRef3::Cartesian(op) => (
                op.prepare_material_mass(POINTS, poll).unwrap(),
                op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap(),
                op.nodes(),
                op.n(),
            ),
            CutElasticityRef3::Adaptive(op) => (
                op.prepare_material_mass(POINTS, poll).unwrap(),
                op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap(),
                op.nodes(),
                op.n(),
            ),
        };
        let occupancy = vec![0.61; mass.cells()];
        let gate = CancelGate::new_clock_free();
        let build = || {
            CutSecondOrderProblem3::new_with_eigenstrain(
                source, &mass, &occupancy, &prepared, 0.02, POINTS, native_n, &gate,
            )
            .unwrap()
        };
        let elastic = build();
        let direction = elastic
            .restrict_force(&nodes.iter().flat_map(|p| [0.0, p[0], 0.0]).collect::<Vec<_>>())
            .unwrap();
        let n = elastic.dimension();
        let seed: Vec<_> = (0..n).map(|i| 0.2 + 0.03 * (i % 11) as f64).collect();
        let mut elastic_action = vec![0.0; n];
        elastic.tangent_apply(&vec![0.0; n], &direction, &mut elastic_action);
        for sign in [-1.0, 0.0, 1.0] {
            let extension = EXPANSION * (1.0 + 0.5 * sign);
            let reference: Vec<_> = nodes
                .iter()
                .flat_map(|p| [extension * p[0], 0.0, 0.0])
                .collect();
            let model = build()
                .with_reference_prestress(&reference, POINTS, POINTS)
                .unwrap();
            let q0 = model.restrict_force(&reference).unwrap();
            let geometric = model.reference_prestress().unwrap();
            match source {
                CutElasticityRef3::Cartesian(op) => op.validate_prestress(geometric).unwrap(),
                CutElasticityRef3::Adaptive(op) => op.validate_prestress(geometric).unwrap(),
            }
            let mut geometric_reference = vec![0.0; native_n];
            geometric
                .apply_into(&reference, &mut geometric_reference, poll)
                .unwrap();
            if sign != 0.0 {
                // An uncentered Gq would incorrectly accelerate this reference.
                assert!(geometric_reference.iter().any(|v| v.abs() > 1e-8));
            }
            let mut expected = vec![0.0; n];
            let mut actual = expected.clone();
            elastic.internal_force(&q0, &mut expected);
            model.internal_force(&q0, &mut actual);
            assert_eq!(actual, expected);
            elastic.mass_apply(&direction, &mut expected);
            model.mass_apply(&direction, &mut actual);
            assert_eq!(actual, expected);

            let mut tangent = vec![0.0; n];
            let mut transpose = tangent.clone();
            model.tangent_apply(&q0, &direction, &mut tangent);
            model
                .tangent_transpose_apply(&q0, &seed, &mut transpose)
                .unwrap();
            close(dot(&seed, &tangent), dot(&direction, &transpose), 2e-12, 2e-12);
            // Independent continuum work for d=(0,x,0): integral sigma_xx.
            // The exact material xx modulus is 15 and retained volume is 0.73.
            close(
                dot(&direction, &tangent) - dot(&direction, &elastic_action),
                scale * 15.0 * (extension - EXPANSION) * 0.73,
                1e-10,
                2e-12,
            );
            let epsilon = 1e-7;
            let plus: Vec<_> = q0.iter().zip(&direction).map(|(q, d)| q + epsilon * d).collect();
            let minus: Vec<_> = q0.iter().zip(&direction).map(|(q, d)| q - epsilon * d).collect();
            model.internal_force(&plus, &mut actual);
            model.internal_force(&minus, &mut expected);
            for ((p, m), value) in actual.iter().zip(&expected).zip(&tangent) {
                close((p - m) / (2.0 * epsilon), *value, 2e-8, 2e-10);
            }
        }
    }
}

#[test]
fn prestressed_midpoint_trajectory_matches_original_dense_mass_and_signed_tangent() {
    let mut op = cartesian(0.1);
    op.set_scales(&[0.37]).unwrap();
    let mass = op.prepare_material_mass(POINTS, poll).unwrap();
    let occupancy = vec![0.61; op.cells()];
    let prepared = op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap();
    let gate = CancelGate::new_clock_free();
    let reference: Vec<_> = equilibrium(op.nodes()).iter().map(|q| 0.5 * q).collect();
    let model = CutSecondOrderProblem3::new_with_eigenstrain(
        CutElasticityRef3::Cartesian(&op), &mass, &occupancy, &prepared,
        0.02, POINTS, op.n(), &gate,
    )
    .unwrap()
    .with_reference_prestress(&reference, POINTS, POINTS)
    .unwrap();
    let n = model.dimension();
    assert!(n <= 24);
    let q0 = model.restrict_force(&reference).unwrap();
    let mut force = vec![0.0; op.n()];
    op.apply(&reference, &mut force);
    for (f, free) in force.iter_mut().zip(op.eigenstrain_load(&prepared, poll).unwrap()) {
        *f -= free;
    }
    let reduced_force = model.restrict_force(&force).unwrap();
    // Bounded independent dense midpoint oracle, using the original native
    // M/K/G actions rather than the dynamics adapter's tangent or recurrence.
    let mut m = vec![0.0; n * n];
    let mut k = m.clone();
    let mut unit = vec![0.0; op.n()];
    let mut action = unit.clone();
    let mut geometric = unit.clone();
    for (column, &dof) in model.free_dofs().iter().enumerate() {
        unit.fill(0.0);
        unit[dof] = 1.0;
        mass.apply_into(&occupancy, &unit, &mut action, poll).unwrap();
        for (row, &native) in model.free_dofs().iter().enumerate() {
            m[row * n + column] = action[native];
        }
        op.apply(&unit, &mut action);
        model.reference_prestress().unwrap()
            .apply_into(&unit, &mut geometric, poll).unwrap();
        for (row, &native) in model.free_dofs().iter().enumerate() {
            k[row * n + column] = action[native] + geometric[native];
        }
    }
    assert!((0..n).any(|i| (0..n).any(|j| i != j && m[i * n + j] != 0.0)));
    let multiply = |matrix: &[f64], x: &[f64]| -> Vec<f64> {
        matrix.chunks_exact(n).map(|row| dot(row, x)).collect()
    };
    let mass_factor = fs_la::factor::lu(&m, n).unwrap();
    let mut displacement: Vec<_> = (0..n).map(|i| 1e-5 * ((i % 5) as f64 - 2.0)).collect();
    let mut velocity: Vec<_> = (0..n).map(|i| 2e-5 * ((i % 7) as f64 - 3.0)).collect();
    let q: Vec<_> = q0.iter().zip(&displacement).map(|(q, d)| q + d).collect();
    let a = acceleration(&model, &q, &velocity, &reduced_force);
    let mut expected_a: Vec<_> = multiply(&k, &displacement).iter().map(|f| -f).collect();
    mass_factor.solve(&mut expected_a);
    for (a, expected) in a.iter().zip(&expected_a) {
        close(*a, *expected, 1e-10, 2e-12);
    }
    let mut state = SecondOrderState::new(0.0, &q, &velocity, &a);
    let dt = 0.04;
    let beta = 0.25 * dt * dt;
    let lhs: Vec<_> = m.iter().zip(&k).map(|(m, k)| m + beta * k).collect();
    let rhs_matrix: Vec<_> = m.iter().zip(&k).map(|(m, k)| m - beta * k).collect();
    let factor = fs_la::factor::lu(&lhs, n).unwrap();
    let stepper = method(n, dt, 1.0);
    for _ in 0..3 {
        let mut next = multiply(&rhs_matrix, &displacement);
        for (value, momentum) in next.iter_mut().zip(multiply(&m, &velocity)) {
            *value += dt * momentum;
        }
        factor.solve(&mut next);
        let next_velocity: Vec<_> = next.iter().zip(&displacement).zip(&velocity)
            .map(|((next, old), v)| 2.0 * (next - old) / dt - v).collect();
        let mut next_acceleration: Vec<_> = multiply(&k, &next).iter().map(|f| -f).collect();
        mass_factor.solve(&mut next_acceleration);
        model.advance(&stepper, &mut state, &force).unwrap();
        for ((actual, reference), expected) in state.q.iter().zip(&q0).zip(&next) {
            close(actual - reference, *expected, 1e-8, 2e-12);
        }
        for (actual, expected) in state.v.iter().zip(&next_velocity)
            .chain(state.a.iter().zip(&next_acceleration)) {
            close(*actual, *expected, 1e-8, 2e-11);
        }
        displacement = next;
        velocity = next_velocity;
    }
}

#[test]
fn reference_refusals_and_physical_endpoint_cancellation_publish_no_state() {
    let op = cartesian(1e-6);
    let mass = op.prepare_material_mass(POINTS, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let prepared = op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap();
    let gate = CancelGate::new_clock_free();
    let source = CutElasticityRef3::Cartesian(&op);
    let build = || CutSecondOrderProblem3::new_with_eigenstrain(
        source, &mass, &occupancy, &prepared, 0.02, POINTS, op.n(), &gate,
    ).unwrap();
    let reference: Vec<_> = op.nodes().iter()
        .flat_map(|p| [(EXPANSION - 0.25e-6) * p[0], 0.0, 0.0]).collect();
    assert!(CutSecondOrderProblem3::new(source, &mass, &occupancy, op.n(), &gate)
        .unwrap().with_reference_prestress(&reference, POINTS, POINTS).is_err());
    assert!(build().with_reference_prestress(&reference, 0, POINTS).is_err());
    assert!(build().with_reference_prestress(&reference, POINTS, 0).is_err());
    assert!(build().with_reference_prestress(&reference[..op.n() - 1], POINTS, POINTS).is_err());
    let mut invalid = reference.clone();
    invalid[0] = f64::NAN;
    assert!(build().with_reference_prestress(&invalid, POINTS, POINTS).is_err());
    let fixed = op.fixed().iter().position(|fixed| *fixed).unwrap();
    invalid = reference.clone();
    invalid[3 * fixed] = 1e-9;
    assert!(build().with_reference_prestress(&invalid, POINTS, POINTS).is_err());
    // The unexpanded state violates the original elastic source calibration.
    assert!(build().with_reference_prestress(&vec![0.0; op.n()], POINTS, POINTS).is_err());
    let model = build().with_reference_prestress(&reference, POINTS, POINTS).unwrap();
    let q0 = model.restrict_force(&reference).unwrap();
    let n = model.dimension();
    let zero = vec![0.0; n];
    let mut resistance = zero.clone();
    model.internal_force(&q0, &mut resistance);
    let force = model.expand_field(&resistance).unwrap();
    let dt = 1e-4;
    let velocity = model.restrict_force(&op.nodes().iter()
        .flat_map(|p| [1.5e-6 / dt * p[0], 0.0, 0.0]).collect::<Vec<_>>()).unwrap();
    let a = acceleration(&model, &q0, &velocity, &resistance);
    let mut state = SecondOrderState::new(0.0, &q0, &velocity, &a);
    let before = state.clone();
    // The midpoint remains inside calibration while the endpoint crosses it.
    assert!(matches!(model.advance(&method(n, dt, 1.0), &mut state, &force),
        Err(TimeSolveError::ModelRefused { phase: "accepted endpoint", .. })));
    assert_eq!(state, before);
    let a = acceleration(&model, &q0, &zero, &resistance);
    let mut state = SecondOrderState::new(0.0, &q0, &zero, &a);
    model.advance(&method(n, dt, 1.0), &mut state, &force).unwrap();
    let accepted = state.clone();
    let mut calls = 0;
    let result = method(n, dt, 1.0).step_controlled(&mut state, &model, &resistance, &mut || {
        calls += 1;
        if calls == 6 {
            gate.request();
        }
        gate.is_requested()
    });
    assert!(matches!(result, Err(TimeSolveError::Cancelled)));
    assert!(calls >= 6);
    assert_eq!(state, accepted);
    assert!(matches!(model.with_reference_prestress(&reference, POINTS, POINTS),
        Err(ElasticityError3::Cancelled)));
}
