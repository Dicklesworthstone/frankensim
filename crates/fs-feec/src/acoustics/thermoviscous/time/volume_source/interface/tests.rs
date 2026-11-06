use super::*;
use fs_rep_mesh::TetComplex;

const VOLUME_M3: f64 = 0.002 * 0.002 * 0.002 / 6.0;

fn limits() -> ThermoviscousTimeLimits {
    let mut limits = ThermoviscousTimeLimits::default();
    // Full matching P2 vector/thermal traces add many independent work
    // observations. Admit their native cold work without dropping any field.
    limits.preparation.max_multiply_terms = 80_000_000_000;
    limits
}

fn geometry(right: bool, permuted: bool) -> (TetComplex, Vec<[f64; 3]>) {
    let (mesh, mut points) = crate::fixtures::single_tet();
    for point in &mut points {
        for coordinate in point {
            *coordinate *= 0.002;
        }
    }
    let mesh = if right {
        points[1][0] *= -1.0;
        TetComplex::from_tets(4, vec![[0, 2, 1, 3]])
    } else {
        mesh
    };
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
                .map(|tet| tet.map(|vertex| new_at_old[vertex as usize]))
                .collect(),
        ),
        old_at_new.map(|old| points[old]).to_vec(),
    )
}

fn model<'c>(
    mesh: &'c TetComplex,
    points: &[[f64; 3]],
    normal_basis: bool,
) -> (ThermoviscousFe<'c>, [u32; 3]) {
    let cut = mesh
        .faces
        .iter()
        .copied()
        .find(|face| face.iter().all(|&vertex| points[vertex as usize][0] == 0.0))
        .unwrap();
    let external: Vec<_> = mesh
        .faces
        .iter()
        .copied()
        .filter(|face| *face != cut)
        .collect();
    let right = points.iter().any(|point| point[0] < 0.0);
    let shapes = if normal_basis {
        if right {
            vec![[1.0, 3.0, 1.0], [-0.5, 0.5, -0.5]]
        } else {
            vec![[1.0; 3], [0.0, 1.0, 0.0]]
        }
    } else {
        vec![[1.0; 3]]
    };
    let mut ports: Vec<_> = shapes
        .into_iter()
        .map(|shape| BoundaryMode {
            traces: vec![(cut, shape)],
        })
        .collect();
    ports.push(BoundaryMode::uniform(&external));
    let model = ThermoviscousFe::new(
        mesh,
        points,
        ThermoviscousMedium {
            density_kg_m3: if right { 1.8 } else { 1.2 },
            pressure_pa: 101325.0,
            temperature_k: 300.0,
            cp_j_kg_k: 1005.0,
            dynamic_viscosity_pa_s: 1.81e-5,
            thermal_conductivity_w_m_k: 0.0257,
            bulk_viscosity_pa_s: 1e-5,
        },
        &ThermoviscousBoundary {
            no_slip_faces: vec![],
            isothermal_faces: vec![],
            ports,
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    (model, cut)
}

fn reservoir() -> [ThermoviscousCaloricReservoir; 1] {
    [ThermoviscousCaloricReservoir {
        heat_capacity_j_kg_k: 40.0,
        relaxation_time_s: 8e-6,
    }]
}

fn passive<'a, 'c>(
    model: &'a ThermoviscousFe<'c>,
    duration: f64,
) -> PreparedThermoviscousTime<'a, 'c> {
    model
        .prepare_time_with_caloric_reservoirs(duration, 965.0, &reservoir(), limits(), &mut || {
            false
        })
        .unwrap()
}

fn driven<'a, 'c>(
    model: &'a ThermoviscousFe<'c>,
    duration: f64,
    forces: &[&ThermoviscousBodyForceSource<'_, '_>],
    heat: &[&ThermoviscousHeatSource<'_, '_>],
) -> PreparedThermoviscousVolumeTime<'a, 'c> {
    model
        .prepare_time_with_volume_sources(
            duration,
            965.0,
            &reservoir(),
            None,
            None,
            forces,
            heat,
            limits(),
            &mut || false,
        )
        .unwrap()
}

fn close_field(left: &[f64], right: &[f64]) {
    assert_eq!(left.len(), right.len());
    let scale = left
        .iter()
        .chain(right)
        .fold(1e-24_f64, |s, v| s.max(v.abs()));
    let error = left
        .iter()
        .zip(right)
        .fold(0.0_f64, |e, (a, b)| e.max((a - b).abs()));
    assert!(
        error < 3e-9 * scale + 1e-25,
        "field error {error:e}, scale {scale:e}"
    );
}

fn close_state(left: &ThermoviscousTimeState, right: &ThermoviscousTimeState) {
    close_field(&left.pressure_pa, &right.pressure_pa);
    let lv: Vec<_> = left.velocity_m_s.iter().flatten().copied().collect();
    let rv: Vec<_> = right.velocity_m_s.iter().flatten().copied().collect();
    close_field(&lv, &rv);
    close_field(&left.temperature_k, &right.temperature_k);
    assert_eq!(
        left.molecular_temperature_k.len(),
        right.molecular_temperature_k.len()
    );
    for (a, b) in left
        .molecular_temperature_k
        .iter()
        .zip(&right.molecular_temperature_k)
    {
        close_field(a, b);
    }
    assert_eq!(
        left.conductive_wall_temperature_k.len(),
        right.conductive_wall_temperature_k.len()
    );
    for (a, b) in left
        .conductive_wall_temperature_k
        .iter()
        .zip(&right.conductive_wall_temperature_k)
    {
        close_field(a, b);
    }
    close_field(
        &left.wall_displacement_m_sqrt_kg,
        &right.wall_displacement_m_sqrt_kg,
    );
    close_field(
        &left.wall_velocity_m_sqrt_kg_per_s,
        &right.wall_velocity_m_sqrt_kg_per_s,
    );
}

fn work_values(work: ThermoviscousTimeWork) -> [f64; 11] {
    [
        work.delta_energy_j,
        work.pressure_work_j,
        work.structural_work_j,
        work.body_force_work_j,
        work.heat_source_work_j,
        work.conductive_boundary_work_j,
        work.fluid_viscous_heat_j,
        work.fluid_conductive_heat_j,
        work.structural_damping_heat_j,
        work.molecular_relaxation_heat_j,
        work.solid_conductive_heat_j,
    ]
}

fn scale(regions: &[ThermoviscousTimeStep; 2]) -> f64 {
    regions
        .iter()
        .map(|region| {
            region.work.initial_energy_j.abs()
                + region.work.final_energy_j.abs()
                + region.work.supplied_j().abs()
                + region.work.dissipated_j().abs()
        })
        .sum::<f64>()
        .max(1e-30)
}

fn close_work(actual: f64, expected: f64, energy: f64) {
    assert!(
        (actual - expected).abs() < 3e-9 * energy,
        "work {actual:e}, expected {expected:e}, scale {energy:e}"
    );
}

fn assert_full_trace(
    pair: &PreparedConformingThermoviscousVolumeTimePair<'_, '_, '_>,
    states: [&ThermoviscousTimeState; 2],
) {
    let layout = pair.interface();
    let velocity_scale = states
        .iter()
        .flat_map(|state| state.velocity_m_s.iter().flatten())
        .fold(1e-24_f64, |s, value| s.max(value.abs()));
    let temperature_scale = states
        .iter()
        .flat_map(|state| &state.temperature_k)
        .fold(1e-24_f64, |s, value| s.max(value.abs()));
    for (&left, &right) in layout.trace_dofs[0].iter().zip(&layout.trace_dofs[1]) {
        for axis in 0..3 {
            assert!(
                (states[0].velocity_m_s[left][axis] - states[1].velocity_m_s[right][axis]).abs()
                    < 3e-9 * velocity_scale
            );
        }
        assert!(
            (states[0].temperature_k[left] - states[1].temperature_k[right]).abs()
                < 3e-9 * temperature_scale
        );
    }
}

#[test]
fn conforming_volume_driven_and_passive_regions_share_full_fields_and_exact_window_work() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, true);
    let (left, lf) = model(&lm, &lp, false);
    let (right, rf) = model(&rm, &rp, false);
    let force = left
        .volume_force_source(
            4,
            &mut |_, _| Ok::<_, ThermoviscousError>([0.0, 1.0 / VOLUME_M3, 0.0]),
            &mut || false,
        )
        .unwrap();
    let heat = left
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME_M3), &mut || false)
        .unwrap();
    let duration = 4e-6;
    let l = driven(&left, duration, &[&force], &[&heat]);
    let r = passive(&right, duration);
    let pair = PreparedConformingThermoviscousVolumeTimePair::connect(
        &l,
        &r,
        &[lf],
        &[rf],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(pair.body_force_source_counts(), [1, 0]);
    assert_eq!(pair.heat_source_counts(), [1, 0]);
    assert_eq!(
        pair.system().state_dim(),
        l.system().state_dim() + r.system().state_dim()
    );
    assert_eq!(pair.external_pressure_modes(), &[(0, 1), (1, 1)]);
    let states = pair.zero_states();
    let original = states.clone();
    let force_drive = [5e-7];
    let heat_drive = [2e-4];
    let frame = pair
        .step(
            [&states[0], &states[1]],
            &[0.0, 0.0],
            [&[], &[]],
            [&force_drive, &[]],
            [&heat_drive, &[]],
            &mut || false,
        )
        .unwrap();
    assert_eq!(states, original);
    assert_full_trace(&pair, [&frame.regions[0].state, &frame.regions[1].state]);
    assert!(frame.regions[1].work.final_energy_j > 0.0);
    assert!(frame.regions[0].work.body_force_work_j.abs() > 0.0);
    assert!(frame.regions[0].work.heat_source_work_j.abs() > 0.0);
    assert_eq!(frame.regions[1].work.body_force_work_j, 0.0);
    assert_eq!(frame.regions[1].work.heat_source_work_j, 0.0);
    assert!(frame.mean_traction_pa.iter().any(|t| t[1].abs() > 0.0));
    assert!(
        frame
            .mean_inward_heat_flux_w_m2
            .iter()
            .any(|q| q.abs() > 0.0)
    );
    let energy = scale(&frame.regions);
    close_work(frame.mechanical_interface_work_j.iter().sum(), 0.0, energy);
    close_work(frame.conductive_interface_work_j.iter().sum(), 0.0, energy);
    close_work(frame.balance_residual_j(), 0.0, energy);
    let caloric = frame.caloric_heat_transfer_j.unwrap();
    assert_ne!(caloric[0], 0.0);
    assert_eq!(caloric[0], -caloric[1]);
    for region in &frame.regions {
        close_work(region.work.balance_residual_j(), 0.0, energy);
        assert!(region.mean_wall_velocity_m_sqrt_kg_per_s.is_empty());
        assert_eq!(region.state.molecular_temperature_k[0].len(), 10);
    }

    let lh = driven(&left, duration * 0.5, &[&force], &[&heat]);
    let rh = passive(&right, duration * 0.5);
    let half = PreparedConformingThermoviscousVolumeTimePair::connect(
        &lh,
        &rh,
        &[lf],
        &[rf],
        limits(),
        &mut || false,
    )
    .unwrap();
    let first = half
        .step(
            [&states[0], &states[1]],
            &[0.0, 0.0],
            [&[], &[]],
            [&force_drive, &[]],
            [&heat_drive, &[]],
            &mut || false,
        )
        .unwrap();
    let second = half
        .step(
            [&first.regions[0].state, &first.regions[1].state],
            &[0.0, 0.0],
            [&[], &[]],
            [&force_drive, &[]],
            [&heat_drive, &[]],
            &mut || false,
        )
        .unwrap();
    for side in 0..2 {
        close_state(&frame.regions[side].state, &second.regions[side].state);
        for ((whole, a), b) in work_values(frame.regions[side].work)
            .into_iter()
            .zip(work_values(first.regions[side].work))
            .zip(work_values(second.regions[side].work))
        {
            close_work(whole, a + b, energy);
        }
        close_work(
            frame.mechanical_interface_work_j[side],
            first.mechanical_interface_work_j[side] + second.mechanical_interface_work_j[side],
            energy,
        );
        close_work(
            frame.conductive_interface_work_j[side],
            first.conductive_interface_work_j[side] + second.conductive_interface_work_j[side],
            energy,
        );
    }

    // Source slots must never leak into the legacy pressure/wall observations.
    // Extra held-input columns may change numerical rounding; compare physical fields.
    let lb = passive(&left, duration);
    let legacy = lb
        .connect_conforming(&r, &[lf], &[rf], limits(), &mut || false)
        .unwrap();
    let mut seeded = pair.zero_states();
    let apex = lp.iter().position(|point| point[0] > 0.0).unwrap();
    seeded[0].temperature_k[apex] = 0.02;
    seeded[0].velocity_m_s[apex][1] = 0.001;
    let zero_source = pair
        .step(
            [&seeded[0], &seeded[1]],
            &[0.3, -0.1],
            [&[], &[]],
            [&[0.0], &[]],
            [&[0.0], &[]],
            &mut || false,
        )
        .unwrap();
    let passive_frame = legacy
        .step(
            [&seeded[0], &seeded[1]],
            &[0.3, -0.1],
            [&[], &[]],
            &mut || false,
        )
        .unwrap();
    let parity_scale = scale(&passive_frame.regions);
    for side in 0..2 {
        close_state(
            &zero_source.regions[side].state,
            &passive_frame.regions[side].state,
        );
        for (actual, expected) in work_values(zero_source.regions[side].work)
            .into_iter()
            .zip(work_values(passive_frame.regions[side].work))
        {
            close_work(actual, expected, parity_scale);
        }
        close_field(
            &zero_source.regions[side].mean_inward_flow_m3_s,
            &passive_frame.regions[side].mean_inward_flow_m3_s,
        );
        assert!(
            zero_source.regions[side]
                .mean_wall_velocity_m_sqrt_kg_per_s
                .is_empty()
        );
    }

    // Equal total count is insufficient: force and heat have distinct units/order.
    assert!(
        pair.step(
            [&states[0], &states[1]],
            &[0.0, 0.0],
            [&[], &[]],
            [&[], &[]],
            [&[0.0, 0.0], &[]],
            &mut || false,
        )
        .is_err()
    );
    assert!(
        pair.step(
            [&states[0], &states[1]],
            &[0.0, 0.0],
            [&[], &[]],
            [&[f64::NAN], &[]],
            [&heat_drive, &[]],
            &mut || false,
        )
        .is_err()
    );
    let mut polls = 0;
    let refused = pair.step(
        [&states[0], &states[1]],
        &[0.0, 0.0],
        [&[], &[]],
        [&force_drive, &[]],
        [&heat_drive, &[]],
        &mut || {
            polls += 1;
            polls == 3
        },
    );
    assert!(matches!(
        refused,
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    assert_eq!(polls, 3);
    assert_eq!(states, original);
    let retry = pair
        .step(
            [&states[0], &states[1]],
            &[0.0, 0.0],
            [&[], &[]],
            [&force_drive, &[]],
            [&heat_drive, &[]],
            &mut || false,
        )
        .unwrap();
    for side in 0..2 {
        assert_eq!(retry.regions[side].state, frame.regions[side].state);
    }
}

#[test]
fn normal_pair_preserves_original_region_and_profile_order_in_motion_and_work() {
    let (lm, lp) = geometry(false, false);
    let (rm, rp) = geometry(true, false);
    let (left, _) = model(&lm, &lp, true);
    let (right, _) = model(&rm, &rp, true);
    let lf0 = left
        .volume_force_source(
            4,
            &mut |_, _| Ok::<_, ThermoviscousError>([1.0 / VOLUME_M3, 0.0, 0.0]),
            &mut || false,
        )
        .unwrap();
    let lf1 = left
        .volume_force_source(
            4,
            &mut |_, point| {
                Ok::<_, ThermoviscousError>([0.0, (point[1] - 0.0005) / (0.002 * VOLUME_M3), 0.0])
            },
            &mut || false,
        )
        .unwrap();
    let lh = left
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME_M3), &mut || false)
        .unwrap();
    let rf = right
        .volume_force_source(
            4,
            &mut |_, _| Ok::<_, ThermoviscousError>([0.0, 0.0, 1.0 / VOLUME_M3]),
            &mut || false,
        )
        .unwrap();
    let rh0 = right
        .volume_heat_source(
            4,
            &mut |_, point| Ok((1.0 + point[2] / 0.002) / VOLUME_M3),
            &mut || false,
        )
        .unwrap();
    let rh1 = right
        .volume_heat_source(
            4,
            &mut |_, point| Ok((point[1] - 0.0005) / (0.002 * VOLUME_M3)),
            &mut || false,
        )
        .unwrap();
    let duration = 3e-6;
    let l = driven(&left, duration, &[&lf0, &lf1], &[&lh]);
    let r = driven(&right, duration, &[&rf], &[&rh0, &rh1]);
    let pair = PreparedThermoviscousVolumeTimePair::connect(
        &l,
        &r,
        &[0, 1],
        &[0, 1],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(pair.body_force_source_counts(), [2, 1]);
    assert_eq!(pair.heat_source_counts(), [1, 2]);
    assert_eq!(pair.external_pressure_modes(), &[(0, 2), (1, 2)]);
    for (actual, expected) in pair.traction_map().iter().zip([1.0, -0.5, 2.0, 1.0]) {
        assert!((actual - expected).abs() < 1e-13);
    }
    let states = pair.zero_states();
    let external = [0.3, -0.2];
    let force_left = [1e-6, -4e-7];
    let force_right = [-6e-7];
    let heat_left = [2e-4];
    let heat_right = [3e-4, -2e-4];
    let frame = pair
        .step(
            [&states[0], &states[1]],
            &external,
            [&[], &[]],
            [&force_left, &force_right],
            [&heat_left, &heat_right],
            &mut || false,
        )
        .unwrap();
    let energy = scale(&frame.regions);
    close_work(frame.interface_work_j.iter().sum(), 0.0, energy);
    close_work(frame.balance_residual_j(), 0.0, energy);
    for side in 0..2 {
        assert!(
            frame.regions[side]
                .mean_wall_velocity_m_sqrt_kg_per_s
                .is_empty()
        );
        assert_eq!(frame.regions[side].mean_inward_flow_m3_s.len(), 3);
        close_work(frame.regions[side].work.balance_residual_j(), 0.0, energy);
    }
    for a in 0..2 {
        let left_flow: f64 = (0..2)
            .map(|i| pair.traction_map()[i * 2 + a] * frame.regions[0].mean_inward_flow_m3_s[i])
            .sum();
        let right_flow = frame.regions[1].mean_inward_flow_m3_s[a];
        assert!(
            (left_flow + right_flow).abs() < 3e-9 * (left_flow.abs() + right_flow.abs()).max(1e-30)
        );
    }

    // Independently reconstruct every original regional differential row at
    // the accepted endpoint, including all source columns and interface reaction.
    let sources = [&l.fields, &r.fields];
    let mut x = sources[0].encode(&frame.regions[0].state).unwrap();
    x.extend(sources[1].encode(&frame.regions[1].state).unwrap());
    let input = [
        external[0],
        external[1],
        force_left[0],
        force_left[1],
        heat_left[0],
        force_right[0],
        heat_right[0],
        heat_right[1],
    ];
    let system = pair.system();
    let n = system.state_dim();
    let effort = system.effort(&x);
    let (j, loss, g) = system.structure();
    let coupled: Vec<_> = (0..n)
        .map(|row| {
            (0..n)
                .map(|column| (j[row * n + column] - loss[row * n + column]) * effort[column])
                .sum::<f64>()
                + (0..input.len())
                    .map(|p| g[row * input.len() + p] * input[p])
                    .sum::<f64>()
        })
        .collect();
    let lambda = &frame.endpoint_interface_traction_pa;
    let original_inputs = [
        vec![
            lambda[0] - 0.5 * lambda[1],
            2.0 * lambda[0] + lambda[1],
            external[0],
            force_left[0],
            force_left[1],
            heat_left[0],
        ],
        vec![
            lambda[0],
            lambda[1],
            external[1],
            force_right[0],
            heat_right[0],
            heat_right[1],
        ],
    ];
    let mut start = 0;
    for side in 0..2 {
        let source = sources[side].system();
        let size = source.state_dim();
        let ports = source.port_dim();
        let e = source.effort(&x[start..start + size]);
        let (j, loss, g) = source.structure();
        assert_eq!(ports, original_inputs[side].len());
        for row in 0..size {
            let drift: Vec<_> = (0..size)
                .map(|column| (j[row * size + column] - loss[row * size + column]) * e[column])
                .collect();
            let supply: Vec<_> = (0..ports)
                .map(|p| g[row * ports + p] * original_inputs[side][p])
                .collect();
            let expected = drift.iter().sum::<f64>() + supply.iter().sum::<f64>();
            let bound = drift
                .iter()
                .chain(&supply)
                .map(|v| v.abs())
                .sum::<f64>()
                .max(1e-30);
            assert!((coupled[start + row] - expected).abs() < 3e-9 * bound);
        }
        start += size;
    }

    // Permuting two mechanical shapes and two thermal shapes with their drives
    // preserves physical state and each source-kind subtotal without normalization.
    let swapped_l = driven(&left, duration, &[&lf1, &lf0], &[&lh]);
    let swapped_r = driven(&right, duration, &[&rf], &[&rh1, &rh0]);
    let swapped = PreparedThermoviscousVolumeTimePair::connect(
        &swapped_l,
        &swapped_r,
        &[0, 1],
        &[0, 1],
        limits(),
        &mut || false,
    )
    .unwrap();
    let permuted = swapped
        .step(
            [&states[0], &states[1]],
            &external,
            [&[], &[]],
            [&[force_left[1], force_left[0]], &force_right],
            [&heat_left, &[heat_right[1], heat_right[0]]],
            &mut || false,
        )
        .unwrap();
    for side in 0..2 {
        close_state(&frame.regions[side].state, &permuted.regions[side].state);
        for (actual, expected) in work_values(frame.regions[side].work)
            .into_iter()
            .zip(work_values(permuted.regions[side].work))
        {
            close_work(actual, expected, energy);
        }
        close_work(
            frame.interface_work_j[side],
            permuted.interface_work_j[side],
            energy,
        );
    }
}
