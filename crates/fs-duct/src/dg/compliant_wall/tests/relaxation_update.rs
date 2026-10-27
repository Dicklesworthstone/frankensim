//! Physical continuation and accounting of retained Maxwell internal variables.
use super::regional::base;
use super::*;
use fs_material::visco::GeneralizedMaxwell;

fn mapping(pairs: &[(usize, usize)]) -> Vec<DgWallRelaxationTransfer> {
    (0..2)
        .flat_map(|element| {
            pairs.iter().map(
                move |&(previous_arm, replacement_arm)| DgWallRelaxationTransfer {
                    element,
                    previous_arm,
                    replacement_arm,
                },
            )
        })
        .collect()
}

#[test]
fn relaxation_time_change_continues_the_independent_three_root_solution_without_work() {
    let gas = duct(DgBoundary::Closed);
    let beta = 2.0 * gas.gas().density * gas.gas().sound_speed.powi(2) / RADIUS;
    let omega: f64 = 10_000.0;
    let sigma = 1.0;
    let pin = WallPin {
        surface_density: sigma,
        stiffness_per_area: sigma * omega.powi(2) - beta,
        resistance: 0.0,
    };
    let arm = 10.0 * sigma * omega.powi(2);
    let old =
        GeneralizedMaxwell::new(pin.stiffness_per_area, vec![(arm, 1.0 / (3.0 * omega))]).unwrap();
    let new =
        GeneralizedMaxwell::new(pin.stiffness_per_area, vec![(arm, 1.0 / (6.0 * omega))]).unwrap();
    let mut dg = gas.with_compliant_wall_relaxation(pin, &old).unwrap();
    let old_propagator = dg.propagator(1e-6, 1e-13).unwrap();
    let m = dg.field_dofs();
    let mut state = dg.zero_state();
    let z0 = 10.0;
    let branch = &dg.compliant_wall.as_ref().unwrap().relaxation[0];
    let coordinates = dg.to_energy_coords(&vec![z0 / branch.pressure_scale; m], &vec![0.0; m]);
    state.tilde[branch.start..branch.start + m].copy_from_slice(&coordinates[..m]);
    state.initial_energy_j = state.stored_energy_j();
    let before = state.clone();
    let old_motion = dg.compliant_wall_state(&state).unwrap();
    let account = dg
        .rebind_compliant_wall_relaxation_regions(
            &mut state,
            &[pin; 2],
            &[Some(&new); 2],
            &mapping(&[(0, 0)]),
            0.0,
            &mut || false,
        )
        .unwrap();
    assert_eq!(account.parameter_work_j, 0.0);
    assert_eq!(state.tilde, before.tilde);
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    let new_motion = dg.compliant_wall_state(&state).unwrap();
    assert_eq!(
        new_motion.relaxation_pressure_pa,
        old_motion.relaxation_pressure_pa
    );
    close(
        new_motion.dissipated_power_w,
        2.0 * old_motion.dissipated_power_w,
        old_motion.dissipated_power_w,
    );
    let saved = state.clone();
    assert!(old_propagator.step(&mut state, &[]).is_err());
    assert_eq!(state, saved);

    // Sealed uniform gas gives p=-beta*w. The coupled characteristic is
    // (s+omega)(s+2*omega)(s+3*omega), with initial w=v=0 and z=z0.
    // This analytic solution includes the original gas/wall feedback.
    let propagator = dg.propagator(1e-6, 1e-13).unwrap();
    let mut previous_energy = state.stored_energy_j();
    for step in 1..=40 {
        let frame = propagator.step(&mut state, &[]).unwrap();
        propagator.sync_nodal(&mut state).unwrap();
        let t = step as f64 * propagator.dt();
        let [a, b, c] = [
            (-omega * t).exp(),
            (-2.0 * omega * t).exp(),
            (-3.0 * omega * t).exp(),
        ];
        let scale = -z0 / (2.0 * sigma * omega.powi(2));
        let w = scale * (a - 2.0 * b + c);
        let v = scale * omega * (-a + 4.0 * b - 3.0 * c);
        let z = z0 * (a - 5.0 * b + 5.0 * c);
        let actual = dg.compliant_wall_state(&state).unwrap();
        for node in 0..m {
            close(actual.displacement_m[node], w, z0 / (sigma * omega.powi(2)));
            close(actual.velocity_m_s[node], v, z0 / (sigma * omega));
            close(actual.relaxation_pressure_pa[0][node], z, z0);
            close(state.pressure()[node], -beta * w, z0);
        }
        assert!(state.stored_energy_j() <= previous_energy);
        previous_energy = state.stored_energy_j();
        assert_eq!(frame.ledger.cumulative_parameter_work_j, 0.0);
        close(
            state.stored_energy_j() + state.cumulative_dissipated_j(),
            before.initial_energy_j(),
            before.initial_energy_j(),
        );
    }
}

#[test]
fn reordered_regional_histories_preserve_extensions_momentum_and_physical_work() {
    let old_pin = wall(30.0);
    let old_law =
        GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001), (0.0, 0.3), (6e5, 0.0004)]).unwrap();
    let new_pins = [
        WallPin {
            surface_density: 0.3,
            stiffness_per_area: 1.6e7,
            resistance: 50.0,
        },
        WallPin {
            surface_density: 0.5,
            stiffness_per_area: 6e6,
            resistance: 70.0,
        },
    ];
    let new_laws = [
        GeneralizedMaxwell::new(1e6, vec![(1.2e6, 0.0002), (0.0, 0.3), (3e5, 0.002)]).unwrap(),
        GeneralizedMaxwell::new(2e6, vec![(4e5, 0.0007), (0.0, 0.3), (8e5, 0.001)]).unwrap(),
    ];
    for spatial in [false, true] {
        let gas = if spatial {
            base(true)
        } else {
            duct(DgBoundary::Port(PortId(0)))
                .with_viscothermal(8, [50.0, 2000.0], 0.01)
                .unwrap()
        };
        let mut dg = gas
            .with_compliant_wall_relaxation(old_pin, &old_law)
            .unwrap();
        let mut state = dg.project_initial_state(|x| 10.0 + x, |x| 1e-7 * (1.0 + x));
        for (i, value) in state.tilde[2 * dg.field_dofs()..].iter_mut().enumerate() {
            *value = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
        }
        state.initial_energy_j = state.stored_energy_j();
        let original = dg.clone();
        let before = state.clone();
        let motion = dg.compliant_wall_state(&state).unwrap();
        let prefix = dg.compliant_wall.as_ref().unwrap().start;
        let account = dg
            .rebind_compliant_wall_relaxation_regions(
                &mut state,
                &new_pins,
                &[Some(&new_laws[0]), Some(&new_laws[1])],
                &mapping(&[(0, 2), (2, 0)]),
                1.0,
                &mut || false,
            )
            .unwrap();
        let after = dg.compliant_wall_state(&state).unwrap();
        assert_eq!(&state.tilde[..prefix], &before.tilde[..prefix]);
        assert_eq!(dg.ports, original.ports);
        assert_eq!(dg.gas_field(), original.gas_field());
        assert_eq!(dg.viscothermal_model(), original.viscothermal_model());
        assert_eq!(motion.relaxation_pressure_pa.len(), 2);
        assert_eq!(after.relaxation_pressure_pa.len(), 4);
        let local = dg.mesh.order + 1;
        for node in 0..dg.field_dofs() {
            let e = node / local;
            close(
                after.displacement_m[node],
                motion.displacement_m[node],
                1e-4,
            );
            close(
                new_pins[e].surface_density * after.velocity_m_s[node],
                old_pin.surface_density * motion.velocity_m_s[node],
                1.0,
            );
            for (old_row, new_term) in [(0, 2), (1, 0)] {
                let old_term = [0, 2][old_row];
                let a = old_pin.stiffness_per_area / old_law.e_inf * old_law.terms[old_term].0;
                let b = new_pins[e].stiffness_per_area / new_laws[e].e_inf
                    * new_laws[e].terms[new_term].0;
                close(
                    after.relaxation_pressure_pa[2 * e + 1 - old_row][node] / b,
                    motion.relaxation_pressure_pa[old_row][node] / a,
                    1e-4,
                );
            }
        }
        // Independent physical geometric integration; no normalized-coordinate
        // energy formula or effective gas/wall mass is reused here.
        let (points, weights) = gauss_legendre_rule(8);
        let mut expected_work = 0.0;
        for (e, element) in dg.elements.iter().enumerate() {
            for (&point, &weight) in points.iter().zip(&weights) {
                let shape = dg.lgl.eval_basis(point);
                let interp = |values: &[f64]| -> f64 {
                    (0..local).map(|i| shape[i] * values[e * local + i]).sum()
                };
                let w = interp(&motion.displacement_m);
                let momentum = old_pin.surface_density * interp(&motion.velocity_m_s);
                let mut density =
                    (new_pins[e].stiffness_per_area - old_pin.stiffness_per_area) * w * w
                        + (1.0 / new_pins[e].surface_density - 1.0 / old_pin.surface_density)
                            * momentum
                            * momentum;
                for (row, old_term, new_term) in [(0, 0, 2), (1, 2, 0)] {
                    let a = old_pin.stiffness_per_area / old_law.e_inf * old_law.terms[old_term].0;
                    let b = new_pins[e].stiffness_per_area / new_laws[e].e_inf
                        * new_laws[e].terms[new_term].0;
                    let extension = interp(&motion.relaxation_pressure_pa[row]) / a;
                    density += (b - a) * extension * extension;
                }
                expected_work +=
                    0.5 * core::f64::consts::TAU * RADIUS * 0.5 * element.length * weight * density;
            }
        }
        close(
            account.parameter_work_j,
            expected_work,
            before.stored_energy_j(),
        );
        close(
            account.energy_after_j - account.energy_before_j,
            expected_work,
            before.stored_energy_j(),
        );
        assert_eq!(
            state.cumulative_parameter_work_j(),
            account.parameter_work_j
        );
        let reverse = dg
            .rebind_compliant_wall_relaxation_regions(
                &mut state,
                &[old_pin; 2],
                &[Some(&old_law); 2],
                &mapping(&[(2, 0), (0, 2)]),
                1.0,
                &mut || false,
            )
            .unwrap();
        assert_eq!(dg.operator_identity(), original.operator_identity());
        close(
            reverse.parameter_work_j,
            -account.parameter_work_j,
            before.stored_energy_j(),
        );
        for (&a, &b) in state.tilde.iter().zip(&before.tilde) {
            close(a, b, 1e-4);
        }
    }
}

#[test]
fn incomplete_duplicate_or_changed_support_and_budget_leave_nonzero_memories_untouched() {
    let law = GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001), (6e5, 0.0004)]).unwrap();
    let original = base(true)
        .with_compliant_wall_relaxation(wall(30.0), &law)
        .unwrap();
    let mut initial = original.zero_state();
    initial.tilde[2 * original.field_dofs()..].fill(1e-6);
    initial.initial_energy_j = initial.stored_energy_j();
    let complete = mapping(&[(0, 0), (1, 1)]);
    let mut duplicate = complete.clone();
    duplicate[1] = duplicate[0];
    let mut outside = complete.clone();
    outside[0].element = 2;
    let changed = GeneralizedMaxwell::new(1e6, vec![(4e5, 0.001), (6e5, 0.0004)]).unwrap();
    let dropped = GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001), (0.0, 0.0004)]).unwrap();
    for (spectrum, entries, budget) in [
        (&law, &complete[..3], 1.0),
        (&law, duplicate.as_slice(), 1.0),
        (&law, outside.as_slice(), 1.0),
        (&changed, complete.as_slice(), 0.0),
        (&dropped, complete.as_slice(), 1.0),
    ] {
        let mut dg = original.clone();
        let mut current = initial.clone();
        assert!(
            dg.rebind_compliant_wall_relaxation_regions(
                &mut current,
                &[wall(30.0); 2],
                &[Some(spectrum); 2],
                entries,
                budget,
                &mut || false
            )
            .is_err()
        );
        assert_eq!(current, initial);
        assert_eq!(dg.operator_identity(), original.operator_identity());
        assert!(
            dg.compliant_wall_state(&current)
                .unwrap()
                .relaxation_stored_energy_j
                > 0.0
        );
    }
}
