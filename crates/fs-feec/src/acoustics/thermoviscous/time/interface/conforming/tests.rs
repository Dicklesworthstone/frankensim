use super::*;

fn geometry(right: bool, permuted: bool) -> (TetComplex, Vec<[f64; 3]>) {
    let (mesh, points) = super::super::tests::geometry(right);
    if !permuted {
        return (mesh, points);
    }
    let old_at_new = [2, 0, 3, 1];
    let mut new_at_old = [0; 4];
    for (new, &old) in old_at_new.iter().enumerate() {
        new_at_old[old] = new as u32;
    }
    (
        TetComplex::from_tets(
            4,
            mesh.tets
                .iter()
                .map(|tet| tet.map(|v| new_at_old[v as usize]))
                .collect(),
        ),
        old_at_new.map(|old| points[old]).to_vec(),
    )
}

fn model<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    wall: bool,
    isothermal: bool,
) -> (ThermoviscousFe<'c>, [u32; 3]) {
    let cut = mesh
        .faces
        .iter()
        .copied()
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
    let result = ThermoviscousFe::new(
        mesh,
        points,
        ThermoviscousMedium {
            density_kg_m3: if points.iter().any(|p| p[0] < 0.0) {
                1.8
            } else {
                1.2
            },
            pressure_pa: 101325.0,
            temperature_k: 300.0,
            cp_j_kg_k: 1005.0,
            dynamic_viscosity_pa_s: 1.81e-5,
            thermal_conductivity_w_m_k: 0.0257,
            bulk_viscosity_pa_s: 1e-5,
        },
        &ThermoviscousBoundary {
            no_slip_faces: clamped.clone(),
            isothermal_faces: if isothermal { clamped } else { vec![] },
            ports: vec![
                BoundaryMode::uniform(&[cut]),
                BoundaryMode::uniform(&external),
            ],
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    (result, cut)
}

fn prepared<'a, 'c>(model: &'a ThermoviscousFe<'c>, dt: f64) -> PreparedThermoviscousTime<'a, 'c> {
    super::super::tests::prepared(model, dt)
}

fn assert_traces(
    pair: &PreparedConformingThermoviscousTimePair<'_, '_, '_>,
    states: [&ThermoviscousTimeState; 2],
) {
    let layout = pair.interface();
    for (&a, &b) in layout.trace_dofs[0].iter().zip(&layout.trace_dofs[1]) {
        for c in 0..3 {
            let scale = states
                .iter()
                .flat_map(|s| s.velocity_m_s.iter().map(|v| v[c].abs()))
                .fold(1e-24, f64::max);
            assert!(
                (states[0].velocity_m_s[a][c] - states[1].velocity_m_s[b][c]).abs() < 2e-9 * scale
            );
        }
        let scale = states
            .iter()
            .flat_map(|s| s.temperature_k.iter().map(|v| v.abs()))
            .fold(1e-24, f64::max);
        assert!((states[0].temperature_k[a] - states[1].temperature_k[b]).abs() < 2e-9 * scale);
    }
}

#[test]
fn complete_original_dae_matches_independent_primal_trace_elimination() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, true);
    let (left, lf) = model(&lm, &lp, false, false);
    let (right, rf) = model(&rm, &rp, false, false);
    let l = prepared(&left, 1e-6);
    let r = prepared(&right, 1e-6);
    let pair = l
        .connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let sizes = [l.system.state_dim(), r.system.state_dim()];
    let n = sizes.iter().sum::<usize>();
    assert_eq!(pair.system().state_dim(), n);
    assert_eq!(pair.inner.interface_count, 24);
    assert_eq!(pair.external_pressure_modes(), &[(0, 1), (1, 1)]);
    let layout = pair.interface();
    let mut parent: Vec<_> = (0..n).collect();
    let mut weights = vec![1.0; n];
    // Independent conforming trial space: identify the original physical
    // coefficients, retaining both original volume pressure/memory spaces.
    for (&a, &b) in layout.trace_dofs[0].iter().zip(&layout.trace_dofs[1]) {
        for c in 0..4 {
            let (ia, ib) = if c < 3 {
                (
                    c * left.velocity_count + left.velocity_map[a].unwrap(),
                    c * right.velocity_count + right.velocity_map[b].unwrap(),
                )
            } else {
                (
                    left.temperature_start + left.temperature_map[a].unwrap(),
                    right.temperature_start + right.temperature_map[b].unwrap(),
                )
            };
            parent[sizes[0] + ib] = ia;
            weights[sizes[0] + ib] = r.scale[ib] / l.scale[ia];
        }
    }
    let retained: Vec<_> = (0..n).filter(|&i| parent[i] == i).collect();
    let k = retained.len();
    let mut e = vec![0.0; n * k];
    for i in 0..n {
        e[i * k + retained.binary_search(&parent[i]).unwrap()] = weights[i];
    }
    let z: Vec<_> = (0..k).map(|i| ((17 * i + 3) % 13) as f64 * 1e-8).collect();
    let x: Vec<_> = (0..n)
        .map(|i| (0..k).map(|j| e[i * k + j] * z[j]).sum())
        .collect();
    pair.inner.projection.admit_state(&x).unwrap();
    let input = [0.7, -0.2];
    let mut q = vec![0.0; n * n];
    let mut force = Vec::new();
    let mut start = 0;
    for (side, source) in [&l, &r].into_iter().enumerate() {
        let size = sizes[side];
        let effort = source.system.effort(&x[start..start + size]);
        let (j, loss, g) = source.system.structure();
        let free: Vec<_> = (0..size)
            .map(|i| {
                (0..size)
                    .map(|jcol| (j[i * size + jcol] - loss[i * size + jcol]) * effort[jcol])
                    .sum::<f64>()
                    + g[i * 2 + 1] * input[side]
            })
            .collect();
        // Q * xdot = the ORIGINAL weak field residual, before any interface
        // projection or multiplier solve. Assemble E^T Q E and E^T residual.
        force.extend(source.system.effort(&free));
        for col in 0..size {
            let mut unit = vec![0.0; size];
            unit[col] = 1.0;
            for (row, value) in source.system.effort(&unit).into_iter().enumerate() {
                q[(start + row) * n + start + col] = value;
            }
        }
        start += size;
    }
    let mut qe = vec![0.0; n * k];
    for i in 0..n {
        for a in 0..k {
            qe[i * k + a] = (0..n).map(|j| q[i * n + j] * e[j * k + a]).sum();
        }
    }
    let mut mass = vec![0.0; k * k];
    let mut rhs = vec![0.0; k];
    for a in 0..k {
        rhs[a] = (0..n).map(|i| e[i * k + a] * force[i]).sum();
        for b in 0..k {
            mass[a * k + b] = (0..n).map(|i| e[i * k + a] * qe[i * k + b]).sum();
        }
    }
    cholesky(&mass, k).unwrap().solve(&mut rhs);
    let expected: Vec<_> = (0..n)
        .map(|i| (0..k).map(|a| e[i * k + a] * rhs[a]).sum::<f64>())
        .collect();
    let (j, loss, g) = pair.system().structure();
    let effort = pair.system().effort(&x);
    let actual: Vec<_> = (0..n)
        .map(|i| {
            (0..n)
                .map(|a| (j[i * n + a] - loss[i * n + a]) * effort[a])
                .sum::<f64>()
                + (0..2).map(|p| g[i * 2 + p] * input[p]).sum::<f64>()
        })
        .collect();
    let norm: f64 = expected.iter().map(|v| v * v).sum();
    let error: f64 = actual
        .iter()
        .zip(&expected)
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    assert!(
        error < 1e-20 * norm,
        "primal original-DAE relative squared error={}",
        error / norm
    );

    // Independently check heat-flux units: q [W/m2] is dual to theta=T'/T0,
    // with the actual surface mass. An extra T0 would pass a mere constraint
    // check but corrupt physical heat-flux amplitudes by a factor of 300.
    let surface = left.fe.boundary_mass(&[lf]).unwrap();
    let a = layout.trace_dofs[0][0];
    let row = left.temperature_start + left.temperature_map[a].unwrap();
    let (columns, values) = surface.row(a);
    let diagonal = values[columns.iter().position(|&col| col == a).unwrap()];
    let observed = pair.inner.projection.constraint_rows()[18 * n + row];
    assert!((observed - diagonal / l.scale[row]).abs() < 1e-11 * (diagonal / l.scale[row]));
}

#[test]
fn full_trace_transmits_shear_and_heat_with_separate_same_trajectory_work() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, true);
    let (left, lf) = model(&lm, &lp, false, false);
    let (right, rf) = model(&rm, &rp, false, false);
    let l = prepared(&left, 4e-6);
    let r = prepared(&right, 4e-6);
    let pair = l
        .connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut states = pair.zero_states();
    let apex = lp.iter().position(|p| p[0] > 0.0).unwrap();
    states[0].temperature_k[apex] = 0.02;
    states[0].velocity_m_s[apex][1] = 0.001;
    let original = states.clone();
    let frame = pair
        .step(
            [&states[0], &states[1]],
            &[0.0, 0.0],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    assert_eq!(states, original);
    assert_traces(&pair, [&frame.regions[0].state, &frame.regions[1].state]);
    assert!(frame.regions[1].work.final_energy_j > 0.0);
    assert!(frame.mean_traction_pa.iter().any(|t| t[1] != 0.0));
    assert!(frame.mean_inward_heat_flux_w_m2.iter().any(|q| *q != 0.0));
    let caloric = frame.caloric_heat_transfer_j.unwrap();
    assert!(caloric[0].abs() > 0.0);
    assert_eq!(caloric[0], -caloric[1]);
    let scale = frame
        .regions
        .iter()
        .map(|r| r.work.initial_energy_j + r.work.final_energy_j)
        .sum::<f64>();
    assert!(frame.mechanical_interface_work_j.iter().sum::<f64>().abs() < 1e-10 * scale);
    assert!(frame.conductive_interface_work_j.iter().sum::<f64>().abs() < 1e-10 * scale);
    assert!(frame.balance_residual_j().abs() < 1e-10 * scale);
    for side in 0..2 {
        assert_eq!(
            frame.regions[side].work.conductive_boundary_work_j,
            frame.conductive_interface_work_j[side]
        );
        assert!(frame.regions[side].work.balance_residual_j().abs() < 1e-10 * scale);
        assert_eq!(
            frame.regions[side].state.molecular_temperature_k[0].len(),
            10
        );
    }
}

#[test]
fn one_sided_stationary_edges_constrain_free_side_and_shared_null_traces_are_omitted() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, true);
    let (left, lf) = model(&lm, &lp, true, true);
    let (right, rf) = model(&rm, &rp, false, false);
    let (clamped_right, crf) = model(&rm, &rp, true, true);
    let l = prepared(&left, 1e-6);
    let r = prepared(&right, 1e-6);
    let cr = prepared(&clamped_right, 1e-6);
    let pair = l
        .connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(pair.interface().traction_indices.len(), 6);
    assert_eq!(pair.interface().heat_flux_indices.len(), 6);
    let states = pair.zero_states();
    let frame = pair
        .step(
            [&states[0], &states[1]],
            &[0.5, -0.2],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    assert!(frame.caloric_heat_transfer_j.is_some());
    assert_traces(&pair, [&frame.regions[0].state, &frame.regions[1].state]);
    for (&a, &b) in pair.interface().trace_dofs[0]
        .iter()
        .zip(&pair.interface().trace_dofs[1])
    {
        if left.velocity_map[a].is_none() {
            assert_eq!(frame.regions[0].state.velocity_m_s[a], [0.0; 3]);
            assert!(
                frame.regions[1].state.velocity_m_s[b]
                    .iter()
                    .all(|v| v.abs() < 1e-14)
            );
        }
        if left.temperature_map[a].is_none() {
            assert_eq!(frame.regions[0].state.temperature_k[a], 0.0);
            assert!(frame.regions[1].state.temperature_k[b].abs() < 1e-14);
        }
    }
    let both = l
        .connect_conforming(
            &cr,
            &[lf],
            &[crf],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(both.interface().traction_indices.len(), 3);
    assert_eq!(both.interface().heat_flux_indices.len(), 3);
    let initial = both.zero_states();
    let frame = both
        .step(
            [&initial[0], &initial[1]],
            &[0.5, -0.2],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    assert_eq!(frame.caloric_heat_transfer_j, None);
}

#[test]
fn incompatible_trace_geometry_clock_and_cancelled_candidates_never_update_sources() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, false);
    let (left, lf) = model(&lm, &lp, false, false);
    let (right, rf) = model(&rm, &rp, false, false);
    let l = prepared(&left, 1e-6);
    let r = prepared(&right, 1e-6);
    let half = prepared(&right, 0.5e-6);
    let limits = ThermoviscousTimeLimits::default();
    assert!(
        l.connect_conforming(&l, &[lf], &[lf], limits, &mut || false)
            .is_err()
    );
    assert!(
        l.connect_conforming(&r, &[lf, lf], &[rf], limits, &mut || false)
            .is_err()
    );
    assert!(
        l.connect_conforming(&half, &[lf], &[rf], limits, &mut || false)
            .is_err()
    );
    assert!(
        l.connect_conforming(
            &r,
            &[lf],
            &[rf],
            ThermoviscousTimeLimits {
                max_states: 1,
                ..limits
            },
            &mut || false
        )
        .is_err()
    );
    let pair = l
        .connect_conforming(&r, &[lf], &[rf], limits, &mut || false)
        .unwrap();
    let mut states = pair.zero_states();
    states[0].temperature_k[pair.interface().trace_dofs[0][0]] = 1e-4;
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
    let original = states.clone();
    let mut count = 0;
    assert!(
        pair.step(
            [&states[0], &states[1]],
            &[1.0, 0.0],
            [&[], &[]],
            &mut || {
                count += 1;
                count == 2
            }
        )
        .is_err()
    );
    assert_eq!(states, original);
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
