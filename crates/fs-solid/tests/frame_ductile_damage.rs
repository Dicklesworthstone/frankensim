//! G1/G3 native frame softening, local work and transactional physical clocks.
//! All coefficients are synthetic implementation controls, not alloy validation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::ductile::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileJ2Conditions, DuctileJ2Settings,
    DuctileWorkEvolution,
};
use fs_material::ductile::uniaxial::{UniaxialDuctileJ2, UniaxialDuctileJ2Law};
use fs_material::ductile::{
    DuctileDomain, JohnsonCookInitiation, JohnsonCookParameters, TabulatedDuctileLocus,
};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::UniaxialJ2;
use fs_solid::fiber::ductile::DuctileFiberLaw;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::ductile_damage::{
    FrameDuctileDamageConfig, FrameDuctileDamageError, FrameDuctileDamagePath,
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const PATTERN: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
const WEIGHTS: [f64; 3] = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];

// This explicitly graded synthetic fixture keeps h/Gp fixed across the three
// original Lobatto receivers to provide one uniform scalar assembly oracle.
// The material tests separately vary h at fixed Gp to check regularization.
fn receiving_model<L: DifferentiableDuctileInitiationLaw + Clone>(
    law: &UniaxialDuctileJ2Law<L>,
    scale: f64,
) -> UniaxialDuctileJ2Law<L> {
    let mut settings = law.model().settings();
    settings.characteristic_length_m *= scale;
    settings.post_initiation_plastic_work_j_m2 *= scale;
    UniaxialDuctileJ2::new(
        law.model().native().clone(),
        law.model().initiation_law().clone(),
        settings,
    )
    .unwrap()
    .with_evolution(law.model().evolution())
    .prepared(law.conditions())
    .unwrap()
}
fn receiving_law(law: &FiberLaw, scale: f64) -> FiberLaw {
    match law {
        FiberLaw::DuctileJ2(DuctileFiberLaw::Tabulated(law)) => {
            FiberLaw::DuctileJ2(receiving_model(law, scale).into())
        }
        FiberLaw::DuctileJ2(DuctileFiberLaw::JohnsonCook(law)) => {
            FiberLaw::DuctileJ2(receiving_model(law, scale).into())
        }
        other => other.clone(),
    }
}
fn context(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 74,
                kernel_id: 310,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx);
    });
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1.0),
        "{a:.16e} != {b:.16e}"
    );
}
fn native() -> UniaxialJ2 {
    UniaxialJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(2000.0, 0.3, 0.05).unwrap(),
            2.0,
            50.0,
            30.0,
        )
        .unwrap(),
    )
    .unwrap()
}
fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-1.0, 1.0],
        strain_rate_s_inv: [1e-10, 1e4],
        temperature_k: [300.0, 300.0],
    }
}
fn settings(length: f64, work: f64) -> DuctileJ2Settings {
    DuctileJ2Settings {
        characteristic_length_m: length,
        post_initiation_plastic_work_j_m2: work,
        maximum_damage_increment: 1.0,
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-9,
            yield_pa: 1e-8,
        },
    }
}
fn conditions() -> DuctileJ2Conditions {
    DuctileJ2Conditions {
        duration_s: 1.0,
        temperature_k: 300.0,
    }
}
fn model(capacity: f64, length: f64, work: f64) -> UniaxialDuctileJ2<TabulatedDuctileLocus> {
    UniaxialDuctileJ2::new(
        native(),
        TabulatedDuctileLocus::new(vec![(-1.0, capacity), (1.0, capacity)], domain()).unwrap(),
        settings(length, work),
    )
    .unwrap()
}
fn frame(law: FiberLaw, carrier: Option<FiberLaw>) -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    let a = frame.add_node(0.0, 0.0);
    let b = frame.add_node(1.0, 0.0);
    frame.fix(a, [true; 3]);
    frame.fix(b, [false, true, true]);
    frame.add_mass(b, 10.0);
    let next = std::cell::Cell::new(0usize);
    frame
        .add_element(a, b, 3, 0.0, &|| {
            let section = next.get();
            next.set(section + 1);
            let law = receiving_law(&law, 2.0 * WEIGHTS[section]);
            let fibers = match &carrier {
                Some(carrier) => vec![
                    Fiber {
                        y: -0.1,
                        area: 1.0,
                        law: carrier.clone(),
                    },
                    Fiber {
                        y: 0.0,
                        area: 1.0,
                        law: law.clone(),
                    },
                    Fiber {
                        y: 0.1,
                        area: 1.0,
                        law: carrier.clone(),
                    },
                ],
                None => [-0.1, 0.1]
                    .map(|y| Fiber {
                        y,
                        area: 1.0,
                        law: law.clone(),
                    })
                    .to_vec(),
            };
            SectionModel::Fiber(Section::new(fibers))
        })
        .unwrap();
    frame
}
fn config() -> FrameDuctileDamageConfig {
    FrameDuctileDamageConfig {
        maximum_points: 6,
        maximum_accepted_steps: 1000,
        maximum_absolute_axial_strain: 0.04,
    }
}
fn path(law: DuctileFiberLaw, cx: &Cx<'_>) -> FrameDuctileDamagePath {
    FrameDuctileDamagePath::new(frame(FiberLaw::DuctileJ2(law), None), 10.0, config(), cx).unwrap()
}

#[test]
fn native_displacement_control_follows_actual_softening_and_retains_local_work() {
    let material = model(0.001, 0.5, 0.01);
    let prepared = material.prepared(conditions()).unwrap();
    context(|_, cx| {
        let mut run = path(prepared.clone().into(), cx);
        let mut all_at_once = run.clone();
        let targets = [0.0005, 0.0015, 0.0025, 0.003, 0.0035, 0.004, 0.005, 0.006];
        let mut state = material.initial_state();
        let mut total = [0.0; 4];
        let mut forces = Vec::new();
        for strain in targets {
            let expected = prepared.response(strain, &state).unwrap();
            let block = run
                .displacement_control(&PATTERN, 3, &[strain], 1.0, cx)
                .unwrap();
            let force = run.frame().element_forces(0)[0];
            close(force, 2.0 * expected.stress_pa, 2e-8);
            close(block.mechanics[0].lambda, force, 2e-8);
            for point in run.points() {
                close(point.axial_strain, strain, 1e-11);
                close(point.damage, expected.state.damage(), 2e-8);
                assert_eq!(point.characteristic_length_m, WEIGHTS[point.section]);
                assert_eq!(point.reference_volume_m3, WEIGHTS[point.section]);
            }
            let account = block.energy;
            close(
                account.mechanical_work_j,
                2.0 * expected.energy.mechanical_work_j_m3,
                2e-8,
            );
            close(
                account.plastic_dissipation_j,
                2.0 * expected.energy.plastic_dissipation_j_m3,
                2e-8,
            );
            close(
                account.damage_release_j,
                2.0 * expected.energy.damage_release_j_m3,
                2e-8,
            );
            close(account.residual_j, 0.0, 1e-8);
            for (sum, increment) in total.iter_mut().zip([
                account.mechanical_work_j,
                account.plastic_dissipation_j,
                account.damage_release_j,
                account.numerical_dissipation_j,
            ]) {
                *sum += increment;
            }
            forces.push(force);
            state = expected.state;
        }
        assert!(
            forces[7] < forces[4],
            "physical post-initiation force must soften"
        );
        assert!(run.points().iter().all(|p| p.damage > 0.0));
        let block = all_at_once
            .displacement_control(&PATTERN, 3, &targets, 1.0, cx)
            .unwrap();
        assert_eq!(all_at_once.time_s(), run.time_s());
        assert_eq!(all_at_once.accepted_steps(), 8);
        for (a, b) in all_at_once.points().iter().zip(run.points()) {
            close(a.damage, b.damage, 1e-10);
        }
        for (actual, expected) in [
            block.energy.mechanical_work_j,
            block.energy.plastic_dissipation_j,
            block.energy.damage_release_j,
            block.energy.numerical_dissipation_j,
        ]
        .into_iter()
        .zip(total)
        {
            close(actual, expected, 2e-8);
        }
        let held = run
            .displacement_control(&PATTERN, 3, &[0.006], 0.25, cx)
            .unwrap();
        close(held.energy.mechanical_work_j, 0.0, 1e-10);
        close(held.energy.damage_release_j, 0.0, 1e-10);
        let damage = run.points()[0].damage;
        run.displacement_control(&PATTERN, 3, &[0.0055], 0.25, cx)
            .unwrap();
        assert!(run.points()[0].damage >= damage);
    });
}

#[test]
fn exact_local_failure_transfers_force_to_remaining_fibers_without_erasure() {
    let material =
        model(0.0002, 0.5, 0.0001).with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic);
    let prepared = material.prepared(conditions()).unwrap();
    let carrier = UniaxialJ2::new(
        J2Plasticity::new(
            IsotropicElastic::new(20000.0, 0.3, 0.05).unwrap(),
            2000.0,
            0.0,
        )
        .unwrap(),
    )
    .unwrap();
    context(|_, cx| {
        let specimen = frame(
            FiberLaw::DuctileJ2(prepared.clone().into()),
            Some(FiberLaw::J2(carrier.clone())),
        );
        let mut run = FrameDuctileDamagePath::new(specimen, 0.0, config(), cx).unwrap();
        let block = run
            .displacement_control(&PATTERN, 3, &[0.0005, 0.001, 0.002, 0.004], 0.5, cx)
            .unwrap();
        assert_eq!(run.points().len(), 3);
        assert!(
            run.points()
                .iter()
                .all(|p| p.damage == 1.0 && p.stress_pa == 0.0)
        );
        close(
            run.frame().element_forces(0)[0],
            2.0 * 20000.0 * 0.004,
            1e-8,
        );
        close(block.energy.residual_j, 0.0, 1e-8);
        let before = run.points().to_vec();
        run.displacement_control(&PATTERN, 3, &[0.003], 0.5, cx)
            .unwrap();
        assert_eq!(run.points().len(), before.len());
        assert!(
            run.points()
                .iter()
                .all(|p| p.damage == 1.0 && p.stress_pa == 0.0)
        );
    });
}

#[test]
fn physical_duration_changes_native_damage_response_and_clone_replay_is_exact() {
    let criterion = JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: 0.0015,
            d2: 0.0,
            d3: 0.0,
            d4: 0.08,
            d5: 0.0,
            reference_strain_rate_s_inv: 0.001,
            transition_temperature_k: 300.0,
            melting_temperature_k: 1000.0,
        },
        DuctileDomain {
            strain_rate_s_inv: [1e-7, 1e4],
            ..domain()
        },
    )
    .unwrap();
    let material = UniaxialDuctileJ2::new(native(), criterion, settings(0.5, 0.01)).unwrap();
    let prepared = material.prepared(conditions()).unwrap();
    context(|_, cx| {
        let mut slow = path(prepared.clone().into(), cx);
        let mut fast = slow.clone();
        let targets = [0.001, 0.002, 0.003, 0.004];
        slow.displacement_control(&PATTERN, 3, &targets, 1.0, cx)
            .unwrap();
        fast.displacement_control(&PATTERN, 3, &targets, 0.1, cx)
            .unwrap();
        assert!(slow.points()[0].damage > fast.points()[0].damage);
        assert!(slow.frame().element_forces(0)[0] < fast.frame().element_forces(0)[0]);
        let mut resumed = fast.clone();
        for run in [&mut fast, &mut resumed] {
            run.displacement_control(&PATTERN, 3, &[0.0045], 0.2, cx)
                .unwrap();
        }
        assert_eq!(fast.points(), resumed.points());
        assert_eq!(
            fast.frame().displacements(),
            resumed.frame().displacements()
        );
        assert_eq!(fast.time_s(), resumed.time_s());
    });
}

#[test]
fn native_motion_and_load_control_keep_the_original_intact_limit() {
    let material = model(1.0, 0.5, 0.01);
    let prepared = material.prepared(conditions()).unwrap();
    context(|_, cx| {
        let mut run = path(prepared.clone().into(), cx);
        let mut baseline = frame(FiberLaw::J2(native()), None);
        let load = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let expected = baseline.static_load(&load, 4).unwrap();
        let actual = run.static_load(&load, 4, 0.25, cx).unwrap();
        for (a, b) in actual.mechanics.iter().zip(&expected) {
            for (&u, &v) in a.u.iter().zip(&b.u) {
                close(u, v, 1e-11);
            }
        }
        let mut split = run.clone();
        let ground = [0.0, -0.3, -0.7, -0.7, -0.2, 0.0, 0.1];
        let damping = Rayleigh { a0: 0.02, a1: 0.0 };
        let expected = baseline
            .newmark_with_load(&ground, 0.002, (1.0, 0.0), damping, &load)
            .unwrap();
        let actual = run
            .newmark_with_load(&ground, 0.002, (1.0, 0.0), damping, &load, cx)
            .unwrap();
        assert_eq!(actual.mechanics.t, expected.t);
        for (a, b) in actual.mechanics.u.iter().zip(&expected.u) {
            for (&u, &v) in a.iter().zip(b) {
                close(u, v, 1e-10);
            }
        }
        for (a, b) in actual
            .mechanics
            .internal_work
            .iter()
            .zip(&expected.internal_work)
        {
            close(*a, *b, 1e-10);
        }
        split
            .newmark_with_load(&ground[..4], 0.002, (1.0, 0.0), damping, &load, cx)
            .unwrap();
        split
            .newmark_with_load(&ground[3..], 0.002, (1.0, 0.0), damping, &load, cx)
            .unwrap();
        assert_eq!(split.time_s(), run.time_s());
        assert_eq!(run.accepted_steps(), 10);
        for (a, b) in split.points().iter().zip(run.points()) {
            close(a.axial_strain, b.axial_strain, 1e-10);
        }
    });
}

#[test]
fn wrong_band_untimed_routes_late_refusal_budget_and_cancellation_preserve_state() {
    context(|gate, cx| {
        let wrong = model(0.001, 0.25, 0.01).prepared(conditions()).unwrap();
        assert!(
            FrameDuctileDamagePath::new(
                frame(FiberLaw::DuctileJ2(wrong.into()), None),
                0.0,
                config(),
                cx
            )
            .is_err()
        );
        let prepared = model(0.001, 0.5, 0.01).prepared(conditions()).unwrap();
        let mut original = frame(FiberLaw::DuctileJ2(prepared.clone().into()), None);
        assert!(original.static_load(&PATTERN, 1).is_err());
        assert!(
            original
                .displacement_control(&PATTERN, 3, &[0.001])
                .is_err()
        );
        assert!(
            original
                .newmark_with_load(
                    &[0.0, 0.0],
                    0.1,
                    (1.0, 0.0),
                    Rayleigh { a0: 0.0, a1: 0.0 },
                    &PATTERN
                )
                .is_err()
        );
        let mut run = FrameDuctileDamagePath::new(original, 0.0, config(), cx).unwrap();
        let before = run.clone();
        assert!(
            run.displacement_control(&PATTERN, 3, &[0.001, 0.08], 0.1, cx)
                .is_err()
        );
        assert_eq!(run.points(), before.points());
        assert_eq!(run.frame().displacements(), before.frame().displacements());
        assert_eq!(run.time_s(), before.time_s());
        assert_eq!(run.accepted_steps(), 0);
        let mut limits = config();
        limits.maximum_accepted_steps = 1;
        let mut bounded =
            FrameDuctileDamagePath::new(before.frame().clone(), 0.0, limits, cx).unwrap();
        assert!(matches!(
            bounded.displacement_control(&PATTERN, 3, &[0.001, 0.002], 0.1, cx),
            Err(FrameDuctileDamageError::Budget)
        ));
        gate.request();
        assert!(matches!(
            run.displacement_control(&PATTERN, 3, &[0.001], 0.1, cx),
            Err(FrameDuctileDamageError::Cancelled)
        ));
        assert_eq!(run.points(), before.points());
    });
}
