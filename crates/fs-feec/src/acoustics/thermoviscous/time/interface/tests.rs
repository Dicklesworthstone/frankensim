use super::*;

fn geometry(right: bool) -> (TetComplex, Vec<[f64; 3]>) {
    let (mesh, mut points) = crate::fixtures::single_tet();
    for p in &mut points {
        for x in p {
            *x *= 0.002;
        }
    }
    if right {
        points[1][0] *= -1.0;
        (TetComplex::from_tets(4, vec![[0, 2, 1, 3]]), points)
    } else {
        (mesh, points)
    }
}

fn model<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    right: bool,
    density: f64,
) -> ThermoviscousFe<'c> {
    let face = [0, 2, 3];
    let modes = if right {
        [[1.0, 3.0, 1.0], [-0.5, 0.5, -0.5]]
    } else {
        [[1.0; 3], [0.0, 1.0, 0.0]]
    };
    let mut ports: Vec<_> = modes
        .iter()
        .map(|&shape| BoundaryMode {
            traces: vec![(face, shape)],
        })
        .collect();
    ports.push(BoundaryMode::uniform(
        &mesh
            .faces
            .iter()
            .copied()
            .filter(|f| *f != face)
            .collect::<Vec<_>>(),
    ));
    ThermoviscousFe::new(
        mesh,
        points,
        ThermoviscousMedium {
            density_kg_m3: density,
            pressure_pa: 101325.0,
            temperature_k: 300.0,
            cp_j_kg_k: 1005.0,
            dynamic_viscosity_pa_s: 1.81e-5,
            thermal_conductivity_w_m_k: 0.0257,
            bulk_viscosity_pa_s: 1.0e-5,
        },
        &ThermoviscousBoundary {
            no_slip_faces: vec![],
            isothermal_faces: vec![],
            ports,
        },
        ThermoviscousLimits {
            max_dofs: MAX_LNS_DOFS,
            max_elements: MAX_LNS_ELEMENTS,
            max_dense_bytes: 16 * 1024 * 1024,
            max_element_aspect_ratio: 20.0,
        },
        &mut || false,
    )
    .unwrap()
}

fn prepared<'a, 'c>(model: &'a ThermoviscousFe<'c>, dt: f64) -> PreparedThermoviscousTime<'a, 'c> {
    model
        .prepare_time_with_caloric_reservoirs(
            dt,
            965.0,
            &[ThermoviscousCaloricReservoir {
                heat_capacity_j_kg_k: 40.0,
                relaxation_time_s: 8e-6,
            }],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap()
}

#[test]
fn g1_original_trace_basis_and_two_field_dae_reactions_match() {
    let (lm, lp) = geometry(false);
    let (rm, rp) = geometry(true);
    let left = model(&lm, &lp, false, 1.2);
    let right = model(&rm, &rp, true, 1.8);
    let l = prepared(&left, 1e-6);
    let r = prepared(&right, 1e-6);
    let pair = l
        .connect(
            &r,
            &[0, 1],
            &[0, 1],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    for (actual, expected) in pair.traction_map().iter().zip([1.0, -0.5, 2.0, 1.0]) {
        assert!((actual - expected).abs() < 1e-13);
    }
    assert_eq!(pair.external_pressure_modes(), &[(0, 2), (1, 2)]);
    assert_eq!(
        pair.system.state_dim(),
        l.system.state_dim() + r.system.state_dim()
    );
    let mut states = pair.zero_states();
    for (side, state) in states.iter_mut().enumerate() {
        for (i, p) in state.pressure_pa.iter_mut().enumerate() {
            *p = (i + side) as f64 * 0.1;
        }
        state.temperature_k[1] = 2e-4;
        state.molecular_temperature_k[0][0] = -1e-4;
    }
    let mut x = l.encode(&states[0]).unwrap();
    x.extend(r.encode(&states[1]).unwrap());
    let input = [0.7, -0.2];
    let lambda = pair.projection.reaction(&x, &input).unwrap();
    let effort = pair.system.effort(&x);
    let n = x.len();
    let (j, rmat, g) = pair.system.structure();
    let coupled: Vec<_> = (0..n)
        .map(|i| {
            (0..n)
                .map(|k| (j[i * n + k] - rmat[i * n + k]) * effort[k])
                .sum::<f64>()
                + (0..2).map(|p| g[i * 2 + p] * input[p]).sum::<f64>()
        })
        .collect();
    let mut start = 0;
    for (side, source) in [&l, &r].iter().enumerate() {
        let count = source.system.state_dim();
        let e = source.system.effort(&x[start..start + count]);
        let (j, r, g) = source.system.structure();
        let p = if side == 0 {
            [
                lambda[0] - 0.5 * lambda[1],
                2.0 * lambda[0] + lambda[1],
                input[0],
            ]
        } else {
            [lambda[0], lambda[1], input[1]]
        };
        for i in 0..count {
            let expected = (0..count)
                .map(|k| (j[i * count + k] - r[i * count + k]) * e[k])
                .sum::<f64>()
                + (0..3).map(|k| g[i * 3 + k] * p[k]).sum::<f64>();
            let scale = (0..count).map(|k| ((j[i * count + k] - r[i * count + k]) * e[k]).abs()).sum::<f64>()
                + (0..3).map(|k| (g[i * 3 + k] * p[k]).abs()).sum::<f64>();
            let scale = scale.max(coupled[start + i].abs()).max(1e-30);
            assert!((coupled[start + i] - expected).abs() < 2e-10 * scale);
        }
        start += count;
    }
    for row in pair.projection.constraint_rows().chunks(n) {
        let residual: f64 = row.iter().zip(&coupled).map(|(a, b)| a * b).sum();
        let bound: f64 = row.iter().zip(&coupled).map(|(a, b)| (a * b).abs()).sum();
        assert!(residual.abs() < 2e-10 * bound.max(1e-30));
    }
}

#[test]
fn g1_causal_transmission_retains_molecular_fields_and_independent_region_work() {
    let (lm, lp) = geometry(false);
    let (rm, rp) = geometry(true);
    let left = model(&lm, &lp, false, 1.2);
    let right = model(&rm, &rp, true, 1.8);
    let l = prepared(&left, 4e-6);
    let r = prepared(&right, 4e-6);
    let pair = l
        .connect(
            &r,
            &[0, 1],
            &[0, 1],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut states = pair.zero_states();
    states[0].temperature_k[1] = 2e-4;
    let original = states.clone();
    let frame = pair
        .step(
            [&states[0], &states[1]],
            &[1.0, 0.0],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    assert_eq!(states, original);
    assert!(frame.regions[1].work.final_energy_j > 0.0);
    assert!(frame.regions[0].work.molecular_relaxation_heat_j > 0.0);
    assert!(frame.regions[0].work.fluid_conductive_heat_j > 0.0);
    assert_eq!(
        frame.regions[1].state.molecular_temperature_k[0].len(),
        right.space.ndof
    );
    let scale = frame
        .regions
        .iter()
        .map(|r| r.work.initial_energy_j + r.work.final_energy_j)
        .sum::<f64>();
    assert!(frame.interface_work_j[0].abs() > 1e-20);
    assert!(frame.interface_work_j.iter().sum::<f64>().abs() < 1e-10 * scale);
    assert!(frame.balance_residual_j().abs() < 1e-10 * scale);
    for region in &frame.regions {
        assert!(region.work.balance_residual_j().abs() < 1e-10 * scale);
    }
    for a in 0..2 {
        let left: f64 = (0..2)
            .map(|i| pair.traction_map[i * 2 + a] * frame.regions[0].mean_inward_flow_m3_s[i])
            .sum();
        let right = frame.regions[1].mean_inward_flow_m3_s[a];
        assert!((left + right).abs() < 2e-10 * (left.abs() + right.abs()).max(1e-24));
    }
}

#[test]
fn g3_complete_coupled_trajectory_and_work_are_invariant_under_clock_subdivision() {
    let (lm, lp) = geometry(false);
    let (rm, rp) = geometry(true);
    let left = model(&lm, &lp, false, 1.2);
    let right = model(&rm, &rp, true, 1.8);
    let l = prepared(&left, 4e-6);
    let r = prepared(&right, 4e-6);
    let pair = l
        .connect(
            &r,
            &[0, 1],
            &[0, 1],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let lhalf = prepared(&left, 2e-6);
    let rhalf = prepared(&right, 2e-6);
    let half = lhalf
        .connect(
            &rhalf,
            &[0, 1],
            &[0, 1],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let states = pair.zero_states();
    let whole = pair
        .step(
            [&states[0], &states[1]],
            &[1.0, -0.2],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    let first = half
        .step(
            [&states[0], &states[1]],
            &[1.0, -0.2],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    let second = half
        .step(
            [&first.regions[0].state, &first.regions[1].state],
            &[1.0, -0.2],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    for (side, source) in [&l, &r].iter().enumerate() {
        let a = source.encode(&whole.regions[side].state).unwrap();
        let b = source.encode(&second.regions[side].state).unwrap();
        let error = a.iter().zip(&b).map(|(a, b)| (a - b).powi(2)).sum::<f64>();
        let norm = a.iter().map(|a| a * a).sum::<f64>();
        assert!(error < 1e-16 * norm);
        let reference = whole.interface_work_j[side];
        assert!(
            (reference - first.interface_work_j[side] - second.interface_work_j[side]).abs()
                < 1e-8 * reference.abs().max(1e-30)
        );
    }
}

#[test]
fn g0_wrong_geometry_shape_clock_capacity_and_initial_flow_refuse_without_updates() {
    let (lm, lp) = geometry(false);
    let (rm, rp) = geometry(true);
    let left = model(&lm, &lp, false, 1.2);
    let right = model(&rm, &rp, true, 1.8);
    let l = prepared(&left, 1e-6);
    let r = prepared(&right, 1e-6);
    let limits = ThermoviscousTimeLimits::default();
    assert!(
        l.connect(&l, &[0, 1], &[0, 1], limits, &mut || false)
            .is_err()
    );
    // Equal areas and nonzero overlap cannot substitute a different trace span.
    assert!(l.connect(&r, &[0], &[0], limits, &mut || false).is_err());
    assert!(
        l.connect(&r, &[0, 0], &[0, 1], limits, &mut || false)
            .is_err()
    );
    let rhalf = prepared(&right, 0.5e-6);
    assert!(
        l.connect(&rhalf, &[0, 1], &[0, 1], limits, &mut || false)
            .is_err()
    );
    assert!(
        l.connect(
            &r,
            &[0, 1],
            &[0, 1],
            ThermoviscousTimeLimits {
                max_states: 1,
                ..limits
            },
            &mut || false
        )
        .is_err()
    );
    let pair = l
        .connect(&r, &[0, 1], &[0, 1], limits, &mut || false)
        .unwrap();
    let mut states = pair.zero_states();
    states[0].velocity_m_s[0][0] = 0.1;
    let original = states.clone();
    assert!(
        pair.step(
            [&states[0], &states[1]],
            &[0.0, 0.0],
            [&[], &[]],
            &mut || false
        )
        .is_err()
    );
    assert_eq!(states, original);
    let states = pair.zero_states();
    let mut polls = 0;
    assert!(
        pair.step(
            [&states[0], &states[1]],
            &[1.0, 0.0],
            [&[], &[]],
            &mut || {
                polls += 1;
                polls == 2
            }
        )
        .is_err()
    );
    assert!(
        pair.step(
            [&states[0], &states[1]],
            &[1.0, 0.0],
            [&[], &[]],
            &mut || false
        )
        .is_ok()
    );
}
