use super::*;

fn geometry(right: bool, permuted: bool) -> (TetComplex, Vec<[f64; 3]>) {
    let (mesh, mut points) = crate::fixtures::single_tet();
    for point in &mut points {
        for value in point {
            *value *= 0.002;
        }
    }
    let tet = if right {
        points[1][0] *= -1.0;
        [0, 2, 1, 3]
    } else {
        mesh.tets[0]
    };
    if !permuted {
        return (TetComplex::from_tets(4, vec![tet]), points);
    }
    let old_at_new = [2, 0, 3, 1];
    let mut new_at_old = [0; 4];
    for (new, &old) in old_at_new.iter().enumerate() {
        new_at_old[old] = new as u32;
    }
    (
        TetComplex::from_tets(4, vec![tet.map(|v| new_at_old[v as usize])]),
        old_at_new.map(|i| points[i]).to_vec(),
    )
}

fn model<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    wall: bool,
) -> (ThermoviscousFe<'c>, [u32; 3]) {
    let cut = *mesh
        .faces
        .iter()
        .find(|face| face.iter().all(|&v| points[v as usize][0] == 0.0))
        .unwrap();
    let clamped: Vec<_> = mesh
        .faces
        .iter()
        .copied()
        .filter(|face| wall && face.iter().all(|&v| points[v as usize][2] == 0.0))
        .collect();
    let external: Vec<_> = mesh
        .faces
        .iter()
        .copied()
        .filter(|face| *face != cut && !clamped.contains(face))
        .collect();
    let right = points.iter().any(|p| p[0] < 0.0);
    let medium = ThermoviscousMedium {
        density_kg_m3: if right { 1.8 } else { 1.2 },
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        cp_j_kg_k: 1005.0,
        dynamic_viscosity_pa_s: if right { 5e-5 } else { 1.8e-5 },
        thermal_conductivity_w_m_k: if right { 0.04 } else { 0.025 },
        bulk_viscosity_pa_s: 1e-5,
    };
    let model = ThermoviscousFe::new(
        mesh,
        points,
        medium,
        &ThermoviscousBoundary {
            no_slip_faces: clamped.clone(),
            isothermal_faces: clamped,
            ports: vec![
                BoundaryMode::uniform(&[cut]),
                BoundaryMode::uniform(&external),
            ],
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    (model, cut)
}

fn coordinates(prepared: &PreparedThermoviscous<'_, '_>, field: &ThermoviscousField) -> Vec<C64> {
    let model = prepared.model;
    let mut result = vec![C64::ZERO; model.n];
    for (global, index) in model.velocity_map.iter().enumerate() {
        if let Some(i) = index {
            for c in 0..3 {
                result[c * model.velocity_count + i] = field.velocity_m_s[global][c];
            }
        }
    }
    result[model.pressure_start..model.temperature_start].copy_from_slice(&field.pressure_pa);
    for (global, index) in model.temperature_map.iter().enumerate() {
        if let Some(i) = index {
            result[model.temperature_start + i] =
                field.temperature_k[global].scale(1.0 / model.medium.temperature_k);
        }
    }
    result
}

fn assert_traces(
    pair: &PreparedConformingThermoviscous<'_, '_, '_>,
    response: &ConformingThermoviscousResponse,
) {
    let layout = pair.interface();
    let velocity_scale = response
        .regions
        .iter()
        .flat_map(|r| r.velocity_m_s.iter().flatten())
        .map(|v| v.abs())
        .fold(1e-30, f64::max);
    let thermal_scale = response
        .regions
        .iter()
        .flat_map(|r| &r.temperature_k)
        .map(|v| v.abs())
        .fold(1e-30, f64::max);
    for (&a, &b) in layout.trace_dofs[0].iter().zip(&layout.trace_dofs[1]) {
        for c in 0..3 {
            assert!(
                (response.regions[0].velocity_m_s[a][c] - response.regions[1].velocity_m_s[b][c])
                    .abs()
                    < 2e-8 * velocity_scale
            );
        }
        assert!(
            (response.regions[0].temperature_k[a] - response.regions[1].temperature_k[b]).abs()
                < 2e-8 * thermal_scale
        );
    }
    let scale = response.external_input_w.abs()
        + response
            .regions
            .iter()
            .map(|r| r.power.dissipated_w())
            .sum::<f64>();
    assert!(response.balance_residual_w().abs() <= 2e-9 * scale);
    for region in &response.regions {
        assert!(region.power.balance_residual_w().abs() < 2e-8 * scale);
    }
}

/// G1/G3: independent primal identification of physical P2 coefficients, not
/// another multiplier solve. Both material pressure spaces remain original.
#[test]
fn original_harmonic_fields_match_primal_identification_with_sources_and_constraints() {
    for (left_wall, right_wall) in [(false, false), (true, false), (true, true)] {
        let (lm, lp) = geometry(false, false);
        let (rm, rp) = geometry(true, true);
        let (left, lf) = model(&lm, &lp, left_wall);
        let (right, rf) = model(&rm, &rp, right_wall);
        let force = left
            .volume_force_source(
                left.volume_quadrature_order(),
                &mut |_, _| Ok::<_, ThermoviscousError>([0.0, 3e5, -2e5]),
                &mut || false,
            )
            .unwrap();
        let heat = left
            .volume_heat_source(
                left.volume_quadrature_order(),
                &mut |_, _| Ok::<_, ThermoviscousError>(8e5),
                &mut || false,
            )
            .unwrap();
        let l = left
            .prepare_with_heat_capacity(4000.0, C64::new(1000.0, -4.0), &mut || false)
            .unwrap();
        let r = right.prepare(4000.0, &mut || false).unwrap();
        let pair = l
            .connect_conforming(
                &r,
                &[lf],
                &[rf],
                ThermoviscousInterfaceLimits::default(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(pair.external_pressure_modes(), &[(0, 1), (1, 1)]);
        let drive = [C64::new(0.7, 0.2), C64::new(-0.2, 0.1)];
        let forces = [(&force, C64::new(0.8, -0.3))];
        let heats = [(&heat, C64::new(0.1, 0.2))];
        let response = pair
            .response_with_volume_sources(
                &drive,
                [&[]; 2],
                [&forces, &[]],
                [&heats, &[]],
                &mut || false,
            )
            .unwrap();
        assert_traces(&pair, &response);

        let prepared = [&l, &r];
        let starts = [0, left.n];
        let n = left.n + right.n;
        let mut parent: Vec<Option<usize>> = (0..n).map(Some).collect();
        for (&a, &b) in pair.interface().trace_dofs[0]
            .iter()
            .zip(&pair.interface().trace_dofs[1])
        {
            for c in 0..4 {
                let indices: [Option<usize>; 2] = std::array::from_fn(|s| {
                    let model = prepared[s].model;
                    let global = if s == 0 { a } else { b };
                    if c < 3 {
                        model.velocity_map[global].map(|i| starts[s] + c * model.velocity_count + i)
                    } else {
                        model.temperature_map[global]
                            .map(|i| starts[s] + model.temperature_start + i)
                    }
                });
                match indices {
                    [Some(a), Some(b)] => parent[b] = Some(a),
                    [Some(i), None] | [None, Some(i)] => parent[i] = None,
                    [None, None] => {}
                }
            }
        }
        let retained: Vec<_> = parent
            .iter()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let ids: Vec<_> = parent
            .iter()
            .map(|p| p.map(|i| retained.binary_search(&i).unwrap()))
            .collect();
        let k = retained.len();
        let mut a = vec![C64::ZERO; k * k];
        let mut rhs = vec![C64::ZERO; k];
        for side in 0..2 {
            let source = prepared[side];
            let model = source.model;
            let mut load: Vec<_> = model.loads[1]
                .iter()
                .map(|&v| drive[side].scale(v))
                .collect();
            heat_source::add_sources(
                model,
                &mut load,
                &[],
                if side == 0 { &forces } else { &[] },
                if side == 0 { &heats } else { &[] },
                &mut || false,
            )
            .unwrap();
            for i in 0..model.n {
                let Some(row) = ids[starts[side] + i] else {
                    continue;
                };
                rhs[row] = rhs[row] + load[i];
                for j in 0..model.n {
                    if let Some(col) = ids[starts[side] + j] {
                        a[row * k + col] = a[row * k + col]
                            + model.operator_entry(
                                i,
                                j,
                                std::f64::consts::TAU * source.frequency_hz,
                                source.cp,
                                &source.wall_relaxation,
                            );
                    }
                }
            }
        }
        // Independent direct primal LU; no Schur, interface loads, or cached factors.
        lu_complex_with_cancel(&a, k, &mut || false)
            .unwrap()
            .solve_with_cancel(&mut rhs, &mut || false)
            .unwrap();
        for side in 0..2 {
            let actual = coordinates(prepared[side], &response.regions[side]);
            let expected: Vec<_> = (0..actual.len())
                .map(|i| ids[starts[side] + i].map_or(C64::ZERO, |j| rhs[j]))
                .collect();
            let scale = expected.iter().map(|z| z.abs()).fold(1e-30, f64::max);
            let error = actual
                .iter()
                .zip(&expected)
                .map(|(&a, &b)| (a - b).abs())
                .fold(0.0, f64::max);
            assert!(error < 3e-8 * scale, "primal error={error} scale={scale}");
        }
    }
}

/// G1: a tangential volume force drives shear in BOTH material regions.
#[test]
fn tangential_excitation_transfers_shear_and_reuses_the_same_frequency_factors() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, false);
    let (left, lf) = model(&lm, &lp, false);
    let (right, rf) = model(&rm, &rp, false);
    let force = left
        .volume_force_source(
            left.volume_quadrature_order(),
            &mut |_, point| Ok::<_, ThermoviscousError>([0.0, 2e8 * point[0], 0.0]),
            &mut || false,
        )
        .unwrap();
    let l = left.prepare(2500.0, &mut || false).unwrap();
    let r = right.prepare(2500.0, &mut || false).unwrap();
    let pair = l
        .connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousInterfaceLimits::default(),
            &mut || false,
        )
        .unwrap();
    let one = [(&force, C64::new(1.0, 0.2))];
    let first = pair
        .response_with_volume_sources(
            &[C64::ZERO; 2],
            [&[]; 2],
            [&one, &[]],
            [&[]; 2],
            &mut || false,
        )
        .unwrap();
    assert_traces(&pair, &first);
    assert!(first.regions[1].power.viscous_w > 0.0);
    assert!(first.traction_pa.iter().any(|t| t[1].abs() > 1e-8));
    assert!(first.mechanical_interface_work_w[0] < 0.0);
    assert!(first.mechanical_interface_work_w[1] > 0.0);
    let phase = C64::new(0.0, 2.0);
    let two = [(&force, one[0].1 * phase)];
    let second = pair
        .response_with_volume_sources(
            &[C64::ZERO; 2],
            [&[]; 2],
            [&two, &[]],
            [&[]; 2],
            &mut || false,
        )
        .unwrap();
    for (a, b) in first.regions[1]
        .velocity_m_s
        .iter()
        .flatten()
        .zip(second.regions[1].velocity_m_s.iter().flatten())
    {
        assert!((*b - *a * phase).abs() <= 1e-8 * (a.abs() + 1e-20));
    }
    assert!(
        (second.external_input_w - 4.0 * first.external_input_w).abs()
            < 1e-8 * second.external_input_w
    );
}

/// G4: failed admissions/cancellation leave both prepared volume owners reusable.
#[test]
fn interface_refuses_mismatch_caps_and_cancellation_without_consuming_sources() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, true);
    let (left, lf) = model(&lm, &lp, false);
    let (right, rf) = model(&rm, &rp, false);
    let l = left.prepare(5000.0, &mut || false).unwrap();
    let r = right.prepare(5000.0, &mut || false).unwrap();
    let other = right.prepare(5001.0, &mut || false).unwrap();
    let limits = ThermoviscousInterfaceLimits::default();
    assert!(
        l.connect_conforming(&other, &[lf], &[rf], limits, &mut || false)
            .is_err()
    );
    assert!(
        l.connect_conforming(&l, &[lf], &[lf], limits, &mut || false)
            .is_err()
    );
    assert!(
        l.connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousInterfaceLimits {
                max_trace_dofs: 23,
                ..limits
            },
            &mut || false
        )
        .is_err()
    );
    assert!(
        l.connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousInterfaceLimits {
                max_dense_bytes: 1,
                ..limits
            },
            &mut || false
        )
        .is_err()
    );
    assert!(
        l.connect_conforming(&r, &[lf], &[rf], limits, &mut || true)
            .is_err()
    );
    let pair = l
        .connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousInterfaceLimits {
                max_trace_dofs: 24,
                ..limits
            },
            &mut || false,
        )
        .unwrap();
    assert!(
        pair.response(&[C64::ZERO], [&[]; 2], &mut || false)
            .is_err()
    );
    assert!(
        pair.response(&[C64::ZERO; 2], [&[]; 2], &mut || true)
            .is_err()
    );
    let response = pair
        .response(&[C64::ZERO; 2], [&[]; 2], &mut || false)
        .unwrap();
    assert_eq!(response.external_input_w, 0.0);
    assert_eq!(response.balance_residual_w(), 0.0);
    assert!(
        response
            .traction_pa
            .iter()
            .flatten()
            .all(|&x| x == C64::ZERO)
    );
}
