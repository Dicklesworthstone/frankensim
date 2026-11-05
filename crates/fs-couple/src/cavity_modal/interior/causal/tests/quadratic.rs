use super::*;

#[test]
fn midside_only_inlet_and_opening_retain_original_force_work_and_volume_injection() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = constant(&geometry);
    let shape = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let inlet = model
        .quadratic_velocity_source(&[([0, 1, 2], shape)], &mut || false)
        .unwrap();
    let opening = model
        .quadratic_velocity_source(&[([0, 1, 3], shape)], &mut || false)
        .unwrap();
    let volume = model
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    let linear_feed = ExactVelocityPort {
        profile: vec![([0, 1, 2], [1.0; 3])],
    };
    let dt = 0.02;
    let mut coupled = CausalInteriorTime::new_with_boundary_sources(
        &model,
        dt,
        &walls(),
        &[linear_feed],
        &[&inlet, &opening],
        &[2],
        &[&volume],
        &[[1, 2, 3]],
        scalar_load(dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(coupled.driven_ports(), [0, 1]);
    assert_eq!(coupled.exterior_ports(), [2]);
    assert_eq!(coupled.interior().port_count(), 3);
    // g=4 lambda0 lambda1 on area 1/2: int(g)=1/6, while the work
    // projection onto the surviving phi0=1-x-y-z is int(phi0*g)=1/15.
    // Every corner value is zero, so a corner-only path would erase both
    // quadratic velocity columns and their mutual work completely.
    let b = 1.0 / 15.0;
    close(inlet.integrated_shape_m2(), 1.0 / 6.0, 3e-14);
    assert!(opening.has_nonzero_projection());
    let (mut q, mut p, mut x) = (0.0, 0.0, 0.0);
    for step in 0..32 {
        let linear = if step < 12 { 0.0 } else { 0.01 };
        let quadratic = if step < 12 { 0.02 } else { -0.005 };
        let injection = -0.0005;
        let source = B * linear + b * quadratic + 0.25 * injection;
        let (baseline, _) = oscillator(q, p, source, dt);
        let (unit, _) = oscillator(0.0, 0.0, b, dt);
        let phi = -(-A * dt).exp_m1() / (A * dt);
        let history = R * x * phi;
        let direct = D + R / A * (1.0 - phi);
        let velocity = (b * (baseline - q) / dt - history) / (b * unit / dt + direct);
        let (qn, pn) = oscillator(q, p, source - b * velocity, dt);
        let xn = (-A * dt).exp() * x - (-A * dt).exp_m1() / A * velocity;
        let result = coupled
            .step(&[linear, quadratic], &[injection], &mut || false)
            .unwrap();
        close(result.outward_velocity_m_s[0], velocity, 2e-10);
        close(
            result.exterior_mean_force_n[0],
            history + direct * velocity,
            2e-9,
        );
        close(result.interior_mean_force_n[0], b * (qn - q) / dt, 2e-9);
        close(coupled.state().pressure_pa()[0], pn, 3e-9);
        close(coupled.state().pressure_potential_pa_s()[0], qn, 2e-10);
        let CausalExteriorLoad::Rational(load) = coupled.exterior() else {
            panic!()
        };
        close(load.pole_states()[0].re, xn, 2e-11);
        let feed_work = (B * linear + b * quadratic) * (qn - q);
        let volume_work = 0.25 * injection * (qn - q);
        close(result.work.boundary_drive_work_j, feed_work, 2e-12);
        close(result.work.volume_source_work_j, volume_work, 2e-12);
        close(result.work.source_work_j, feed_work + volume_work, 2e-12);
        close(
            result.work.exterior_exchange_work_j,
            velocity * b * (qn - q),
            2e-12,
        );
        close(result.work.balance_residual_j(), 0.0, 2e-12);
        if step == 0 {
            assert!(result.outward_velocity_m_s[0] > 1e-7);
        }
        (q, p, x) = (qn, pn, xn);
    }
}

fn equivalent_fields(
    actual: &CausalInteriorTime<'_, '_, '_>,
    expected: &CausalInteriorTime<'_, '_, '_>,
) {
    for (a, b) in actual
        .state()
        .pressure_pa()
        .iter()
        .zip(expected.state().pressure_pa())
    {
        close(*a, *b, 3e-9);
    }
    for (a, b) in actual
        .state()
        .pressure_potential_pa_s()
        .iter()
        .zip(expected.state().pressure_potential_pa_s())
    {
        close(*a, *b, 3e-10);
    }
}

#[test]
fn prepared_linear_recovery_and_live_zero_net_midsides_keep_complete_mutual_history() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = constant(&geometry);
    let feed = ExactVelocityPort {
        profile: vec![([0, 1, 2], [1.0, 0.0, -0.4])],
    };
    let first = ExactVelocityPort {
        profile: vec![([0, 1, 3], [1.0; 3])],
    };
    let second = ExactVelocityPort {
        profile: vec![([0, 2, 3], [1.0, -1.0, 0.0])],
    };
    let uniform = model
        .quadratic_velocity_source(&[([0, 1, 3], [1.0; 6])], &mut || false)
        .unwrap();
    let linear = model
        .quadratic_velocity_source(
            &[([0, 2, 3], [1.0, -1.0, 0.0, 0.0, -0.5, 0.5])],
            &mut || false,
        )
        .unwrap();
    // Pure midside motion on the same opening: equal opposite edge shapes
    // integrate to zero while their original P2 pressure projection is live.
    let midsides = model
        .quadratic_velocity_source(&[([0, 2, 3], [0.0, 0.0, 0.0, 1.0, -1.0, 0.0])], &mut || {
            false
        })
        .unwrap();
    close(midsides.integrated_shape_m2(), 0.0, 3e-14);
    assert!(midsides.has_nonzero_projection());
    let volume = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    let dt = 0.01;
    let source = delayed_source(dt);
    let mut original = CausalInteriorTime::new(
        &model,
        dt,
        &[],
        &[feed.clone(), first.clone(), second.clone()],
        &[2, 1],
        &[&volume],
        &[],
        delayed_load(&source, dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut recovered = CausalInteriorTime::new_with_boundary_sources(
        &model,
        dt,
        &[],
        &[feed.clone()],
        &[&uniform, &linear],
        &[2, 1],
        &[&volume],
        &[],
        delayed_load(&source, dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut quadratic = CausalInteriorTime::new_with_boundary_sources(
        &model,
        dt,
        &[],
        &[feed],
        &[&uniform, &midsides],
        &[2, 1],
        &[&volume],
        &[],
        delayed_load(&source, dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut history: Vec<Vec<f64>> = Vec::new();
    let mut live = false;
    for step in 0..24 {
        let drive = 0.015 * (0.41 * f64::from(step)).cos();
        let injection = -0.0001;
        let expected = original
            .step(&[drive], &[injection], &mut || false)
            .unwrap();
        let actual = recovered
            .step(&[drive], &[injection], &mut || false)
            .unwrap();
        equivalent_fields(&recovered, &original);
        for (a, b) in actual
            .outward_velocity_m_s
            .iter()
            .zip(&expected.outward_velocity_m_s)
        {
            close(*a, *b, 3e-10);
        }
        close(
            actual.work.source_work_j,
            expected.work.source_work_j,
            3e-12,
        );
        close(
            actual.work.exterior_exchange_work_j,
            expected.work.exterior_exchange_work_j,
            3e-12,
        );
        let result = quadratic
            .step(&[drive], &[injection], &mut || false)
            .unwrap();
        history.push(result.outward_velocity_m_s.clone());
        live |= result.outward_velocity_m_s[0].abs() > 1e-7;
        assert_eq!(quadratic.state().pressure_pa().len(), 10);
        // Independent continuous step-response integrals confirm both original
        // mutual coordinates and fractional delays after the spatial change.
        for a in 0..2 {
            let mut expected_force = 0.0;
            for b in 0..2 {
                let entry = 2 * a + b;
                for k in 0..=step as usize {
                    let jump = history[k][b] - if k == 0 { 0.0 } else { history[k - 1][b] };
                    let start = (step as usize - k) as f64 * dt - source.delays[entry];
                    expected_force += jump
                        * (step_primitive(&source.model, entry, start + dt)
                            - step_primitive(&source.model, entry, start))
                        / dt;
                }
            }
            close(result.exterior_mean_force_n[a], expected_force, 3e-10);
            close(result.interior_mean_force_n[a], expected_force, 3e-10);
        }
        close(result.work.balance_residual_j(), 0.0, 3e-12);
    }
    assert!(live);
    let CausalExteriorLoad::Delayed(load) = quadratic.exterior() else {
        panic!()
    };
    assert_eq!(load.source_model(), &source);
    assert_eq!(load.reports().0.class, CertificateClass::GridOnly);
}

#[test]
fn prepared_source_identity_support_zero_projection_and_combined_caps_refuse() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = constant(&geometry);
    let foreign_model = constant(&geometry);
    let shape = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let opening = model
        .quadratic_velocity_source(&[([0, 1, 3], shape)], &mut || false)
        .unwrap();
    let zero_feed = model
        .quadratic_velocity_source(&[([3, 0, 1], [0.0; 6])], &mut || false)
        .unwrap();
    let foreign = foreign_model
        .quadratic_velocity_source(&[([0, 1, 3], shape)], &mut || false)
        .unwrap();
    let dt = 0.02;
    let prepare = |sources: &[&BoundaryVelocitySource<'_, '_, '_>],
                   selection: &[usize],
                   release: &[[u32; 3]]| {
        CausalInteriorTime::new_with_boundary_sources(
            &model,
            dt,
            &[],
            &[],
            sources,
            selection,
            &[],
            release,
            scalar_load(dt),
            limits(),
            &mut || false,
        )
    };
    assert!(prepare(&[&foreign], &[0], &[]).is_err());
    assert!(prepare(&[&zero_feed], &[0], &[]).is_err());
    assert!(prepare(&[&opening, &zero_feed], &[0], &[]).is_err());
    assert!(prepare(&[&opening], &[1], &[]).is_err());
    assert!(prepare(&[&opening], &[0], &[[3, 1, 0]]).is_err());
    assert!(
        CausalInteriorTime::new_with_boundary_sources(
            &model,
            dt,
            &[ResistiveBoundary {
                faces: vec![[1, 3, 0]],
                admittance_m_pa_s: 0.0
            }],
            &[],
            &[&opening],
            &[0],
            &[],
            &[],
            scalar_load(dt),
            limits(),
            &mut || false
        )
        .is_err()
    );
    // Count BOTH original source families before allocating any combined map.
    let feed = ExactVelocityPort {
        profile: vec![([0, 1, 2], [1.0; 3])],
    };
    assert!(
        CausalInteriorTime::new_with_boundary_sources(
            &model,
            dt,
            &[],
            &vec![feed.clone(); 128],
            &[&opening],
            &[128],
            &[],
            &[],
            scalar_load(dt),
            limits(),
            &mut || false
        )
        .is_err()
    );
    // Explicit empty prepared-source calls preserve every legacy operation.
    let linear_opening = ExactVelocityPort {
        profile: vec![([0, 1, 3], [1.0; 3])],
    };
    let ports = [feed, linear_opening];
    let mut legacy = CausalInteriorTime::new(
        &model,
        dt,
        &[],
        &ports,
        &[1],
        &[],
        &[],
        scalar_load(dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut empty = CausalInteriorTime::new_with_boundary_sources(
        &model,
        dt,
        &[],
        &ports,
        &[],
        &[1],
        &[],
        &[],
        scalar_load(dt),
        limits(),
        &mut || false,
    )
    .unwrap();
    for _ in 0..4 {
        let a = legacy.step(&[0.02], &[], &mut || false).unwrap();
        let b = empty.step(&[0.02], &[], &mut || false).unwrap();
        assert_eq!(legacy.state().pressure_pa(), empty.state().pressure_pa());
        assert_eq!(
            legacy.state().pressure_potential_pa_s(),
            empty.state().pressure_potential_pa_s()
        );
        assert_eq!(a.outward_velocity_m_s, b.outward_velocity_m_s);
        assert_eq!(
            a.work.source_work_j.to_bits(),
            b.work.source_work_j.to_bits()
        );
    }
}
