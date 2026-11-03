//! Original P2 mass, moving supports and retained mixed Maxwell memory.
use super::*;
use fs_solid::hyper3::mixed::visco::dynamics::SolidImpulseControls;
use fs_solver::CheckedCgConfig;

const HOLD: f64 = 0.007;

fn impulse_controls() -> SolidImpulseControls {
    SolidImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        maximum_workspace_components: 4096,
        linear_tolerance_kg_m_s: 1e-11,
        angular_tolerance_kg_m2_s: 1e-11,
        work_tolerance_j: 1e-11,
    }
}

// Independent exact integrals of the quadratic barycentric basis on an affine
// tetrahedron: corner/edge cross entries are negative, despite positive density.
fn mass(a: usize, b: usize) -> f64 {
    let coefficient = if a < 4 && b < 4 {
        if a == b { 6.0 } else { 1.0 }
    } else if a < 4 || b < 4 {
        let (corner, edge) = if a < 4 { (a, b - 4) } else { (b, a - 4) };
        if TET10_EDGES[edge].contains(&corner) {
            -4.0
        } else {
            -6.0
        }
    } else if a == b {
        32.0
    } else if TET10_EDGES[a - 4]
        .iter()
        .any(|v| TET10_EDGES[b - 4].contains(v))
    {
        16.0
    } else {
        8.0
    };
    DENSITY * VOLUME * coefficient / 420.0
}

fn relaxed<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a ViscoHyperelastic,
    cx: &Cx<'_>,
) -> MixedViscoTetSolid<'a> {
    let fixed: Vec<_> = nodes
        .iter()
        .enumerate()
        .flat_map(|(i, x)| {
            [
                0.002 * x[0] + 0.002 * x[1] * x[2],
                -0.001 * x[1],
                0.001 * x[2],
            ]
            .into_iter()
            .enumerate()
            .map(move |(c, u)| (3 * i + c, u))
        })
        .collect();
    let mut value = solid(
        nodes,
        elements,
        material,
        VolumetricTreatment::OriginalMaterial,
        cx,
    );
    for dt in [0.0, HOLD] {
        value = value
            .advance(
                dt,
                &fixed,
                &[],
                [0.0; 3],
                None,
                settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
    }
    assert!(value.quadrature_states().iter().all(|s| s.dissipated > 0.0));
    assert!(value.pressure_pa().iter().any(|p| p.abs() > 0.01));
    value
}

fn retained(before: &MixedViscoTetMotion<'_>, after: &MixedViscoTetMotion<'_>) {
    assert_eq!(
        after.solid().displacement_m(),
        before.solid().displacement_m()
    );
    assert_eq!(
        after.solid().pressure_nodes(),
        before.solid().pressure_nodes()
    );
    assert_eq!(after.solid().pressure_pa(), before.solid().pressure_pa());
    assert_eq!(
        after.solid().quadrature_states(),
        before.solid().quadrature_states()
    );
    assert_eq!(
        after.solid().first_piola_pa(),
        before.solid().first_piola_pa()
    );
    assert_eq!(
        after.solid().remainder_energy_j(),
        before.solid().remainder_energy_j()
    );
    assert_eq!(
        after.solid().bulk_energy_j(),
        before.solid().bulk_energy_j()
    );
    assert_eq!(after.solid().time_s(), before.solid().time_s());
    assert_eq!(
        after.reference_density_kg_m3(),
        before.reference_density_kg_m3()
    );
    assert_eq!(after.reference_mass_kg(), before.reference_mass_kg());
}

#[test]
fn edge_impulse_keeps_pressure_and_maxwell_history_with_signed_support_work() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let old_v = vec![[0.1, -0.02, 0.03]; nodes.len()];
        let motion = MixedViscoTetMotion::try_new(
            relaxed(&nodes, &elements, &material, cx),
            &[DENSITY],
            &old_v,
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let edge = TET10_EDGES
            .iter()
            .position(|e| e.contains(&2) && e.contains(&3))
            .unwrap();
        let local = 4 + edge;
        let node = elements[0][local];
        let free_dof = 3 * node;
        let held: Vec<_> = (0..3 * nodes.len()).filter(|&d| d != free_dof).collect();
        let delta_v = 0.02;
        let impulse = mass(local, local) * delta_v;
        let origin = [0.4, -0.3, 0.2];
        let event = motion
            .impulse_candidate(
                &held,
                &[(free_dof, 0.3 * impulse), (free_dof, 0.7 * impulse)],
                origin,
                impulse_controls(),
                cx,
            )
            .unwrap();
        retained(&motion, &event.motion);
        for (i, actual) in event.motion.velocity_mps().iter().enumerate() {
            for c in 0..3 {
                near(
                    actual[c],
                    old_v[i][c] + if 3 * i + c == free_dof { delta_v } else { 0.0 },
                    2e-12,
                );
            }
        }
        let mut support = vec![0.0; 3 * nodes.len()];
        for &(dof, value) in &event.impulse.support_impulses_kg_m_s {
            support[dof] = value;
        }
        let mut support_work = 0.0;
        for a in 0..10 {
            let node_a = elements[0][a];
            for c in 0..3 {
                let expected = if c == 0 && a != local {
                    mass(a, local) * delta_v
                } else {
                    0.0
                };
                near(support[3 * node_a + c], expected, 2e-13);
                support_work += expected * old_v[node_a][c];
            }
        }
        assert!(support[3 * elements[0][0]] < 0.0);
        assert!(support_work > 0.0);
        near(
            event.impulse.applied_work_j,
            impulse * (old_v[node][0] + 0.5 * delta_v),
            2e-13,
        );
        near(event.impulse.support_work_j, support_work, 2e-13);
        near(
            event.impulse.kinetic_change_j,
            event.impulse.applied_work_j + support_work,
            2e-13,
        );
        near(event.impulse.work_residual_j, 0.0, 2e-13);
        let x = std::array::from_fn::<_, 3, _>(|c| {
            nodes[node][c] + motion.solid().displacement_m()[node][c] - origin[c]
        });
        let expected_torque = [0.0, x[2] * impulse, -x[1] * impulse];
        for c in 0..3 {
            near(
                event.impulse.applied.angular_kg_m2_s[c],
                expected_torque[c],
                2e-13,
            );
            near(event.impulse.linear_residual_kg_m_s[c], 0.0, 2e-13);
            near(event.impulse.angular_residual_kg_m2_s[c], 0.0, 2e-13);
        }
        assert_eq!(motion.velocity_mps(), old_v);

        // The kicked quadratic shear coordinate remains free; the other DOFs
        // continue their accepted support velocities through the original solve.
        let dt = 0.001;
        let targets: Vec<_> = held
            .iter()
            .map(|&d| {
                (
                    d,
                    motion.solid().displacement_m()[d / 3][d % 3] + dt * old_v[d / 3][d % 3],
                )
            })
            .collect();
        let next = event
            .motion
            .advance(
                dt,
                &targets,
                &[],
                [0.0; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_ne!(
            next.motion.solid().displacement_m()[node][0],
            motion.solid().displacement_m()[node][0]
        );
        for (before, after) in motion
            .solid()
            .quadrature_states()
            .iter()
            .zip(next.motion.solid().quadrature_states())
        {
            let original = material
                .step_remainder_response(before, &after.deformation_gradient, dt)
                .unwrap();
            assert_eq!(*after, original.state);
        }
        near(next.motion.solid().time_s(), HOLD + dt, 1e-15);
        work_identity(&next);
    });
}

#[test]
fn signed_corner_impulses_produce_one_uniform_native_velocity_and_continue_freely() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            ),
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let velocity = [0.3, -0.2, 0.1];
        let impulses: Vec<_> = (0..10)
            .flat_map(|a| {
                let integrated_shape = if a < 4 { -VOLUME / 20.0 } else { VOLUME / 5.0 };
                let node = elements[0][a];
                (0..3).map(move |c| (3 * node + c, DENSITY * integrated_shape * velocity[c]))
            })
            .collect();
        let event = motion
            .impulse_candidate(&[], &impulses, [0.0; 3], impulse_controls(), cx)
            .unwrap();
        retained(&motion, &event.motion);
        for v in event.motion.velocity_mps() {
            for c in 0..3 {
                near(v[c], velocity[c], 2e-12);
            }
        }
        for c in 0..3 {
            near(
                event.impulse.after.linear_kg_m_s[c],
                DENSITY * VOLUME * velocity[c],
                2e-12,
            );
        }
        let torque = [
            0.25 * (velocity[2] - velocity[1]),
            0.25 * (velocity[0] - velocity[2]),
            0.25 * (velocity[1] - velocity[0]),
        ];
        for c in 0..3 {
            near(
                event.impulse.after.angular_kg_m2_s[c],
                DENSITY * VOLUME * torque[c],
                2e-12,
            );
        }
        near(
            event.impulse.kinetic_after_j,
            0.5 * DENSITY * VOLUME * velocity.iter().map(|v| v * v).sum::<f64>(),
            2e-12,
        );
        assert!(event.impulse.support_impulses_kg_m_s.is_empty());
        let dt = 0.01;
        let next = event
            .motion
            .advance(dt, &[], &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        for (u, v) in next
            .motion
            .solid()
            .displacement_m()
            .iter()
            .zip(next.motion.velocity_mps())
        {
            for c in 0..3 {
                near(u[c], dt * velocity[c], 2e-10);
                near(v[c], velocity[c], 2e-8);
            }
        }
        near(next.material_energy.dissipated_j, 0.0, 2e-12);
        work_identity(&next);
    });
}

#[test]
fn exact_pressure_and_failed_impulses_preserve_the_complete_accepted_motion() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let old_v = vec![[0.1, 0.0, 0.0]; nodes.len()];
        let motion = MixedViscoTetMotion::try_new(
            relaxed(&nodes, &elements, &material, cx),
            &[DENSITY],
            &old_v,
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let snapshot = || {
            (
                motion.solid().displacement_m().to_vec(),
                motion.solid().pressure_pa().to_vec(),
                motion.solid().quadrature_states().to_vec(),
                motion.solid().first_piola_pa().to_vec(),
                motion.velocity_mps().to_vec(),
                motion.solid().time_s(),
                motion.solid().remainder_energy_j(),
                motion.solid().bulk_energy_j(),
                motion.reference_density_kg_m3().to_vec(),
                motion.reference_mass_kg(),
            )
        };
        let before = snapshot();
        let impulse = [(0, 0.001)];
        let mut limited = impulse_controls();
        limited.maximum_workspace_components = 0;
        assert!(
            motion
                .impulse_candidate(&[], &impulse, [0.0; 3], limited, cx)
                .is_err()
        );
        assert!(
            motion
                .impulse_candidate(&[1, 1], &impulse, [0.0; 3], impulse_controls(), cx)
                .is_err()
        );
        assert!(
            motion
                .impulse_candidate(&[], &[(0, f64::NAN)], [0.0; 3], impulse_controls(), cx)
                .is_err()
        );
        let cancel = CancelGate::new_clock_free();
        cancel.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let stopped = Cx::new(
                &cancel,
                arena,
                StreamKey {
                    seed: 61,
                    kernel_id: 22,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                motion
                    .impulse_candidate(&[], &impulse, [0.0; 3], impulse_controls(), &stopped)
                    .is_err()
            );
        });
        assert_eq!(snapshot(), before);
        let first = motion
            .impulse_candidate(&[], &impulse, [0.0; 3], impulse_controls(), cx)
            .unwrap();
        let retry = motion
            .impulse_candidate(&[], &impulse, [0.0; 3], impulse_controls(), cx)
            .unwrap();
        assert_eq!(first.motion.velocity_mps(), retry.motion.velocity_mps());
        retained(&motion, &first.motion);

        let exact = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::ExactIncompressible,
                cx,
            ),
            &[DENSITY],
            &old_v,
            MASS_BYTES,
            cx,
        )
        .unwrap();
        assert!(matches!(
            exact.impulse_candidate(&[], &impulse, [0.0; 3], impulse_controls(), cx),
            Err(MixedViscoTetError::ExactIncompressibleImpulse)
        ));
        assert_eq!(exact.velocity_mps(), old_v);
        assert_eq!(exact.solid().time_s(), 0.0);
        assert!(exact.solid().pressure_pa().iter().all(|p| *p == 0.0));
    });
}

fn constrained_controls() -> fs_solid::hyper3::mixed::visco::dynamics::MixedImpulseControls {
    let mut balance = impulse_controls();
    balance.solve.tolerance = 1e-11;
    fs_solid::hyper3::mixed::visco::dynamics::MixedImpulseControls {
        balance,
        velocity_scale_m_s: 1.0,
        pressure_impulse_scale_pa_s: 1.0,
        work_scale_j: 1.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
        maximum_relative_volume_rate_per_s: 1e-10,
    }
}

#[test]
fn exact_pressure_kick_retains_relaxed_stretched_memory_and_native_continuation() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e300, true); // The unused bulk penalty must remain unused.
        let stretch: f64 = 1.02;
        let held: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(i, x)| {
                (0..3).filter_map(move |c| {
                    if c == 0 {
                        Some((3 * i, (stretch - 1.0) * x[0]))
                    } else {
                        (x[c] == 0.0).then_some((3 * i + c, 0.0))
                    }
                })
            })
            .collect();
        let mut accepted = solid(
            &nodes,
            &elements,
            &material,
            VolumetricTreatment::ExactIncompressible,
            cx,
        );
        for dt in [0.0, HOLD] {
            accepted = accepted
                .advance(
                    dt,
                    &held,
                    &[],
                    [0.0; 3],
                    None,
                    settings(),
                    controls().material,
                    cx,
                )
                .unwrap()
                .solid;
        }
        assert!(
            accepted
                .quadrature_states()
                .iter()
                .all(|s| s.dissipated > 0.0)
        );
        assert!(accepted.pressure_pa().iter().any(|p| p.abs() > 0.01));
        let old_v = [0.02, 0.01, -0.01];
        let delta = [0.03, -0.02, 0.01];
        let motion = MixedViscoTetMotion::try_new(
            accepted,
            &[DENSITY],
            &vec![old_v; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let gradients = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let cofactor = [1.0 / stretch, stretch.sqrt(), stretch.sqrt()];
        let pressure_impulse = 0.7;
        let mut applied = Vec::new();
        for (a, &node) in elements[0].iter().enumerate() {
            let row_mass = DENSITY * VOLUME * if a < 4 { -0.05 } else { 0.2 };
            for c in 0..3 {
                // Independent affine geometry: integral grad(N_corner)=0;
                // integral grad(N_edge)=V*(grad(L_i)+grad(L_j)).
                let volume = if a < 4 {
                    0.0
                } else {
                    let [i, j] = TET10_EDGES[a - 4];
                    VOLUME * cofactor[c] * (gradients[i][c] + gradients[j][c])
                };
                applied.push((
                    3 * node + c,
                    row_mass * delta[c] + pressure_impulse * volume,
                ));
            }
        }
        let event = motion
            .constrained_impulse_candidate(
                &[],
                &applied,
                None,
                [0.0; 3],
                constrained_controls(),
                cx,
            )
            .unwrap();
        retained(&motion, &event.motion);
        let inertia = event.impulse.kkt_inertia.unwrap();
        assert_eq!(inertia.positive, 30);
        assert_eq!(inertia.negative, 4);
        for v in event.motion.velocity_mps() {
            for c in 0..3 {
                near(v[c], old_v[c] + delta[c], 2e-9);
            }
        }
        for p in &event.impulse.pressure_impulses_pa_s {
            near(*p, pressure_impulse, 2e-8);
        }
        for rate in event
            .impulse
            .volume_rate_before_m3_s
            .iter()
            .chain(&event.impulse.volume_rate_after_m3_s)
        {
            near(
                *rate,
                0.0,
                constrained_controls().maximum_relative_volume_rate_per_s * VOLUME,
            );
        }
        for c in 0..3 {
            near(
                event.impulse.external.linear_kg_m_s[c],
                DENSITY * VOLUME * delta[c],
                2e-11,
            );
            near(event.impulse.pressure.linear_kg_m_s[c], 0.0, 2e-11);
            near(event.impulse.pressure.angular_kg_m2_s[c], 0.0, 2e-11);
        }
        let kinetic_change = 0.5
            * DENSITY
            * VOLUME
            * (0..3)
                .map(|c| (old_v[c] + delta[c]).powi(2) - old_v[c].powi(2))
                .sum::<f64>();
        near(event.impulse.external_work_j, kinetic_change, 2e-10);
        near(event.impulse.pressure_work_j, 0.0, 1e-11);
        near(
            event.impulse.momentum.kinetic_change_j,
            kinetic_change,
            2e-10,
        );

        let dt = 0.01;
        let next = event
            .motion
            .advance(dt, &[], &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        for (before, after) in motion
            .solid()
            .quadrature_states()
            .iter()
            .zip(next.motion.solid().quadrature_states())
        {
            let original = material
                .step_remainder_response(before, &after.deformation_gradient, dt)
                .unwrap();
            assert_eq!(*after, original.state);
        }
        near(next.motion.solid().time_s(), HOLD + dt, 1e-15);
        work_identity(&next);
    });
}

#[test]
fn constrained_entry_refuses_wrong_material_incompatible_rates_gauges_and_budgets() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let motion = |treatment, velocity: &[[f64; 3]]| {
            MixedViscoTetMotion::try_new(
                solid(&nodes, &elements, &material, treatment, cx),
                &[DENSITY],
                velocity,
                MASS_BYTES,
                cx,
            )
            .unwrap()
        };
        let zero = vec![[0.0; 3]; nodes.len()];
        let finite = motion(VolumetricTreatment::OriginalMaterial, &zero);
        assert!(matches!(
            finite.constrained_impulse_candidate(
                &[],
                &[],
                None,
                [0.0; 3],
                constrained_controls(),
                cx
            ),
            Err(MixedViscoTetError::FiniteCompliancePressureImpulse)
        ));
        let expansion: Vec<_> = nodes.iter().map(|x| x.map(|value| 0.02 * value)).collect();
        let incompatible = motion(VolumetricTreatment::ExactIncompressible, &expansion);
        let history = incompatible.solid().quadrature_states().to_vec();
        assert!(
            incompatible
                .constrained_impulse_candidate(&[], &[], None, [0.0; 3], constrained_controls(), cx)
                .is_err()
        );
        assert_eq!(incompatible.velocity_mps(), expansion);
        assert_eq!(incompatible.solid().quadrature_states(), history);
        let exact = motion(VolumetricTreatment::ExactIncompressible, &zero);
        let applied = [(0, 0.001)];
        assert!(
            exact
                .constrained_impulse_candidate(
                    &[],
                    &applied,
                    Some(elements[0][0]),
                    [0.0; 3],
                    constrained_controls(),
                    cx
                )
                .is_err()
        );
        let mut limited = constrained_controls();
        limited.maximum_factor_bytes = 1;
        assert!(
            exact
                .constrained_impulse_candidate(&[], &applied, None, [0.0; 3], limited, cx)
                .is_err()
        );
        assert_eq!(exact.velocity_mps(), zero);
        assert_eq!(exact.solid().time_s(), 0.0);
        assert!(exact.solid().pressure_pa().iter().all(|p| *p == 0.0));
    });
}

const EXACT_STRETCH: f64 = 1.04;

fn relaxed_exact<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a ViscoHyperelastic,
    cx: &Cx<'_>,
) -> MixedViscoTetSolid<'a> {
    let held: Vec<_> = nodes
        .iter()
        .enumerate()
        .flat_map(|(node, p)| {
            (0..3).filter_map(move |c| {
                if c == 0 {
                    Some((3 * node, (EXACT_STRETCH - 1.0) * p[0]))
                } else {
                    (p[c] == 0.0).then_some((3 * node + c, 0.0))
                }
            })
        })
        .collect();
    let mut value = solid(
        nodes,
        elements,
        material,
        VolumetricTreatment::ExactIncompressible,
        cx,
    );
    for dt in [0.0, HOLD] {
        value = value
            .advance(
                dt,
                &held,
                &[],
                [0.0; 3],
                None,
                settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
    }
    assert!(value.quadrature_states().iter().all(|s| s.dissipated > 0.0));
    assert!(value.pressure_pa().iter().any(|p| p.abs() > 0.01));
    value
}

#[test]
fn constrained_impulse_keeps_moving_support_work_and_a_true_pressure_gauge() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e300, true);
        let old_velocity = [0.1, 0.0, 0.0];
        let motion = MixedViscoTetMotion::try_new(
            relaxed_exact(&nodes, &elements, &material, cx),
            &[DENSITY],
            &vec![old_velocity; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        // Three tangent edge directions leave one true constant-pressure gauge.
        let free: Vec<_> = (0..3)
            .map(|c| {
                let edge = TET10_EDGES
                    .iter()
                    .position(|e| e.contains(&0) && e.contains(&(c + 1)))
                    .unwrap();
                3 * elements[0][4 + edge] + c
            })
            .collect();
        let held: Vec<_> = (0..3 * nodes.len()).filter(|d| !free.contains(d)).collect();
        let load = [(free[0], 0.1)];
        let event = motion
            .constrained_impulse_candidate(
                &held,
                &load,
                Some(elements[0][0]),
                [0.0; 3],
                constrained_controls(),
                cx,
            )
            .unwrap();
        retained(&motion, &event.motion);
        for v in event.motion.velocity_mps() {
            for c in 0..3 {
                near(v[c], old_velocity[c], 1e-10);
            }
        }
        // B[:,edge01.x]=V/(5*s)*[1,-1,0,0], with lambda0 fixed to zero.
        near(
            event.impulse.pressure_impulses_pa_s[1],
            -3.0 * EXACT_STRETCH,
            2e-8,
        );
        near(event.impulse.external_work_j, 0.01, 1e-11);
        near(event.impulse.momentum.support_work_j, -0.01, 1e-11);
        near(event.impulse.momentum.support.linear_kg_m_s[0], -0.1, 1e-11);
        near(event.impulse.momentum.kinetic_change_j, 0.0, 1e-11);
        near(event.impulse.pressure_work_j, 0.0, 1e-11);
        assert_eq!(motion.velocity_mps(), vec![old_velocity; nodes.len()]);
    });
}
