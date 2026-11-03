//! G1/G3: graded reference inertia in actual Maxwell memory and motion.
use super::*;

const NODAL_DENSITY: [f64; 4] = [30.0, 150.0, 30.0, 30.0];

fn prepared_solid<'a>(
    law: &'a ViscoHyperelastic,
    extension: f64,
    hold_s: f64,
    cx: &Cx<'_>,
) -> FiniteViscoTetSolid<'a> {
    let mut solid = FiniteViscoTetSolid::try_new(
        &NODES,
        &TETS,
        ViscoelasticLawField::Uniform(law),
        TetAssemblyBudget::default(),
        cx,
    )
    .unwrap();
    if extension != 0.0 {
        let fixed: Vec<_> = (0..12)
            .map(|dof| (dof, if dof == 3 { extension } else { 0.0 }))
            .collect();
        for dt in [0.0, hold_s] {
            solid = solid
                .advance(
                    dt,
                    &fixed,
                    &[],
                    [0.0; 3],
                    settings(),
                    controls().material,
                    cx,
                )
                .unwrap()
                .solid;
        }
    }
    solid
}

#[test]
fn graded_gravity_uses_exact_mass_shares_and_drives_the_native_free_body() {
    with_cx(|cx| {
        let law = material(true);
        let velocity0 = [0.2, -0.1, 0.05];
        let acceleration = [0.3, -0.2, 0.1];
        let dt = 0.02;
        let mut motion = FiniteViscoTetMotion::try_new_with_nodal_density(
            prepared_solid(&law, 0.0, 0.0, cx),
            &NODAL_DENSITY,
            &[velocity0; 4],
            cx,
        )
        .unwrap();
        close(motion.reference_mass_kg(), 10.0, 1e-13);
        for interval in 1..=3 {
            let loads = motion
                .reference_acceleration_load_n(&[acceleration; 4], cx)
                .unwrap();
            let mut nodal_load = [0.0; 12];
            for &(dof, force) in &loads {
                nodal_load[dof] += force;
            }
            // Independent integral of rho*N_i for rho=30+120*x. Equal cell
            // density would incorrectly assign 2.5 kg to every gravity row.
            for (node, mass) in [2.25, 3.25, 2.25, 2.25].into_iter().enumerate() {
                for c in 0..3 {
                    close(nodal_load[3 * node + c], mass * acceleration[c], 1e-14);
                }
            }
            let step = motion
                .advance(dt, &[], &loads, [0.0; 3], settings(), controls(), cx)
                .unwrap();
            let count = f64::from(interval);
            for node in 0..4 {
                for c in 0..3 {
                    close(
                        step.motion.velocity_mps()[node][c],
                        velocity0[c] + count * dt * acceleration[c],
                        2e-9,
                    );
                    close(
                        step.motion.solid().displacement_m()[node][c],
                        count * dt * velocity0[c]
                            + 0.5 * count * (count + 1.0) * dt * dt * acceleration[c],
                        2e-11,
                    );
                }
            }
            assert!(step.equilibrium.reactions_n.is_empty());
            assert_eq!(step.motion.reference_density_kg_m3(), &[60.0]);
            assert_eq!(
                step.motion.reference_nodal_density_kg_m3(),
                Some(NODAL_DENSITY.as_slice())
            );
            assert_eq!(step.motion.reference_mass_kg(), motion.reference_mass_kg());
            close(step.material_energy.stored_after_j, 0.0, 1e-12);
            close(step.material_energy.dissipated_j, 0.0, 1e-15);
            close(step.balance_error_j, 0.0, 2e-9);
            motion = step.motion;
        }
    });
}

#[test]
fn graded_oscillator_preserves_relaxed_memory_and_uses_exact_inertia() {
    with_cx(|cx| {
        let law = material(true);
        let amplitude = 1e-4;
        let hold = 0.007;
        let dt = 0.0025;
        let solid = prepared_solid(&law, amplitude, hold, cx);
        let old_states = solid.element_states().to_vec();
        let old_displacement = solid.displacement_m().to_vec();
        assert!(old_states[0].dissipated > 0.0);
        let mut motion = FiniteViscoTetMotion::try_new_with_nodal_density(
            solid,
            &NODAL_DENSITY,
            &[[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        assert_eq!(motion.solid().element_states(), old_states);
        assert_eq!(motion.solid().displacement_m(), old_displacement);
        assert_eq!(motion.solid().time_s(), hold);
        let strict = FiniteViscoMotionControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            motion.advance(dt, &supports(), &[], [0.0; 3], settings(), strict, cx),
            Err(FiniteViscoTetError::EnergyGate {
                quantity: "backward-Euler numerical kinetic loss",
                ..
            })
        ));
        assert_eq!(motion.solid().element_states(), old_states);
        assert_eq!(motion.solid().time_s(), hold);
        // Exact integral rho*N_1^2 = V*(sum(rho)+2*rho_1)/60 = 1.5 kg.
        // Replacing rho by its cell mean gives 1 kg and the wrong frequency.
        let mass = 1.5;
        let stiffness = 4.0 * MU * VOLUME;
        let branch_stiffness = 4.0 * BRANCH_MU * VOLUME / 3.0;
        let (mut u, mut v, mut q) = (
            amplitude,
            0.0,
            branch_stiffness * amplitude / (1.0 + hold / TAU),
        );
        let mut total_dissipation = old_states[0].dissipated * VOLUME;
        for _ in 0..20 {
            let ratio = 1.0 + dt / TAU;
            v = (mass * v - dt * (stiffness * u + q / ratio))
                / (mass + dt * dt * (stiffness + branch_stiffness / ratio));
            u += dt * v;
            q = (q + branch_stiffness * dt * v) / ratio;
            let previous = motion.solid().element_states()[0].clone();
            let step = motion
                .advance(dt, &supports(), &[], [0.0; 3], settings(), controls(), cx)
                .unwrap();
            close(
                step.motion.solid().displacement_m()[1][0],
                u,
                5e-4 * amplitude,
            );
            close(step.motion.velocity_mps()[1][0], v, 2e-6);
            let actual = &step.motion.solid().element_states()[0];
            let native = law
                .step_response(&previous, &actual.deformation_gradient, dt)
                .unwrap();
            assert_eq!(*actual, native.state);
            close(
                step.material_energy.dissipated_j,
                VOLUME * native.dissipated_pa,
                1e-15,
            );
            let actual_v = step.motion.velocity_mps()[1][0];
            close(
                step.inertia.kinetic_after_j,
                0.5 * mass * actual_v * actual_v,
                1e-14,
            );
            assert_eq!(step.motion.reference_mass_kg(), motion.reference_mass_kg());
            assert_eq!(
                step.motion.reference_nodal_density_kg_m3(),
                Some(NODAL_DENSITY.as_slice())
            );
            total_dissipation += step.material_energy.dissipated_j;
            motion = step.motion;
        }
        close(
            motion.solid().element_states()[0].dissipated * VOLUME,
            total_dissipation,
            1e-14,
        );
        close(motion.solid().time_s(), hold + 20.0 * dt, 1e-15);
    });
}

#[test]
fn uniform_nodal_density_retains_native_bits_and_invalid_fields_refuse() {
    with_cx(|cx| {
        let law = material(true);
        let velocity = [[0.2, -0.1, 0.05]; 4];
        let original = initial(&law, 0.0, velocity, cx);
        let bind = |density: &[f64], velocity: &[[f64; 3]]| {
            FiniteViscoTetMotion::try_new_with_nodal_density(
                prepared_solid(&law, 0.0, 0.0, cx),
                density,
                velocity,
                cx,
            )
        };
        let nodal = bind(&[DENSITY; 4], &velocity).unwrap();
        assert_eq!(original.reference_nodal_density_kg_m3(), None);
        assert_eq!(original.reference_mass_kg(), nodal.reference_mass_kg());
        let acceleration = [[0.3, -0.2, 0.1]; 4];
        let loads = original
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        assert_eq!(
            loads,
            nodal
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );
        let old = original
            .advance(0.02, &[], &loads, [0.0; 3], settings(), controls(), cx)
            .unwrap();
        let new = nodal
            .advance(0.02, &[], &loads, [0.0; 3], settings(), controls(), cx)
            .unwrap();
        assert_eq!(
            old.motion.solid().displacement_m(),
            new.motion.solid().displacement_m()
        );
        assert_eq!(old.motion.velocity_mps(), new.motion.velocity_mps());
        assert_eq!(
            old.motion.solid().element_states(),
            new.motion.solid().element_states()
        );
        assert_eq!(old.material_energy, new.material_energy);
        assert_eq!(old.inertia, new.inertia);
        for density in [
            &[1.0, 2.0, 3.0][..],
            &[1.0, 0.0, 2.0, 3.0],
            &[1.0, f64::NAN, 2.0, 3.0],
        ] {
            assert!(bind(density, &velocity).is_err());
        }
        assert!(bind(&NODAL_DENSITY, &velocity[..3]).is_err());
        assert!(bind(&NODAL_DENSITY, &[[f64::INFINITY; 3]; 4]).is_err());
        assert!(nodal.reference_acceleration_load_n(&[], cx).is_err());
        assert!(
            nodal
                .reference_acceleration_load_n(&[[f64::NAN; 3]; 4], cx)
                .is_err()
        );
        let solid = prepared_solid(&law, 0.0, 0.0, cx);
        let cancel = CancelGate::new_clock_free();
        cancel.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &cancel,
                arena,
                StreamKey {
                    seed: 42,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                FiniteViscoTetMotion::try_new_with_nodal_density(
                    solid,
                    &NODAL_DENSITY,
                    &velocity,
                    &cancelled,
                ),
                Err(FiniteViscoTetError::Native(
                    fs_solid::hyper3::HyperTetError::Cancelled
                ))
            ));
            assert!(matches!(
                nodal.reference_acceleration_load_n(&acceleration, &cancelled),
                Err(FiniteViscoTetError::Native(
                    fs_solid::hyper3::HyperTetError::Cancelled
                ))
            ));
        });
        assert_eq!(nodal.velocity_mps(), velocity);
        assert_eq!(
            nodal.reference_nodal_density_kg_m3(),
            Some([DENSITY; 4].as_slice())
        );
    });
}
