//! Native finite creep/recovery with actual free nodal motion and explicit time.
//! Manufactured Hencky/kinetic references verify the implementation, not an
//! alloy calibration, rupture model or convergence of a locking-prone P1 mesh.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::StrainHardeningCreep;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::creep::recovery::FiniteRecoveringNortonCreep;
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_material::plastic::finite::{
    FiniteInelasticMaterial, FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters,
};
use fs_solid::hyper3::plastic::dynamics::{FinitePlasticMotionControls, FinitePlasticTetMotion};
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetError, FinitePlasticTetSolid,
    FinitePlasticTetStep,
};
use fs_solid::{HyperTetError, HyperTetSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const VOLUME: f64 = 1.0 / 6.0;
const MU: f64 = 8000.0;
const BULK: f64 = 16000.0;
const ETA: f64 = 4000.0;
const DT: f64 = 0.2;
const PRIMARY: f64 = 0.003;
const RECOVERY: f64 = 2.0;
const PIOLA: f64 = 1000.0;
const COMPLIANCE: f64 = 1.0 / (9.0 * BULK) + 1.0 / (3.0 * MU);

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 189,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.25,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 2.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}

fn kinetics() -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: ETA,
        stress_scale_pa: PIOLA,
        exponent: 1.0,
    }
}

fn recovering(duration_s: f64) -> FiniteRecoveringNortonCreep {
    let native = FiniteNortonCreep::new(MU, BULK, limits(), kinetics(), duration_s)
        .unwrap()
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: PRIMARY,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap();
    FiniteRecoveringNortonCreep::new(native, RECOVERY).unwrap()
}

fn initial<'a, L: FiniteInelasticMaterial>(
    material: &'a L,
    cx: &Cx<'_>,
) -> FinitePlasticTetSolid<'a, L> {
    FinitePlasticTetSolid::try_new(
        &NODES,
        &TETS,
        FinitePlasticLawField::Uniform(material),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}

fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        ..HyperTetSettings::default()
    }
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-7,
        maximum_endpoint_work_defect_j: 25.0,
        yield_tolerance_pa: 1e-6,
    }
}

fn supports() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|d| ![3, 7, 11].contains(d))
        .map(|d| (d, 0.0))
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

#[derive(Clone, Copy, Default)]
struct Point {
    axial_log: f64,
    lateral_log: f64,
    creep: f64,
    resistance: f64,
    q: f64,
    dissipation: f64,
    flow_potential: f64,
}

// Independent free-sided Hencky + primary/recovery equations. At endpoint
// q=P*exp(a), a=p+q/E, and dc*(a_eff+r_old+dc)=a_eff*dt*q/eta.
// This oracle does not call the production constitutive or tangent evaluator.
fn loaded(old: Point) -> Point {
    let divisor = 1.0 + RECOVERY * DT;
    let effective_primary = PRIMARY * divisor;
    let increment = |q: f64| {
        let shifted = effective_primary + old.resistance;
        let drive = effective_primary * DT * q / ETA;
        2.0 * drive / (shifted + (shifted * shifted + 4.0 * drive).sqrt())
    };
    let residual = |q: f64| (q / PIOLA).ln() - old.creep - increment(q) - COMPLIANCE * q;
    let (mut lo, mut hi) = (PIOLA, 2.0 * PIOLA);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..100 {
        let q = 0.5 * (lo + hi);
        if residual(q) < 0.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let dc = increment(q);
    let creep = old.creep + dc;
    Point {
        axial_log: creep + COMPLIANCE * q,
        lateral_log: -0.5 * creep + q / (9.0 * BULK) - q / (6.0 * MU),
        creep,
        resistance: (old.resistance + dc) / divisor,
        q,
        dissipation: q * dc,
        // Integral of the inverse rate, distinct from endpoint heat q*dc.
        flow_potential: (ETA / DT)
            * (0.5 * (1.0 + old.resistance / effective_primary) * dc * dc
                + dc * dc * dc / (3.0 * effective_primary)),
    }
}

fn unloaded(old: Point) -> Point {
    Point {
        axial_log: old.creep,
        lateral_log: -0.5 * old.creep,
        creep: old.creep,
        resistance: old.resistance / (1.0 + RECOVERY * DT),
        ..Point::default()
    }
}

fn check(step: &FinitePlasticTetStep<'_, FiniteRecoveringNortonCreep>, old: Point, next: Point) {
    let state = &step.solid.element_states()[0];
    close(state.native.alpha, next.creep, 2e-9);
    close(state.resistance_strain, next.resistance, 2e-9);
    assert_eq!(state.native.creep_alpha, 0.0);
    for (actual, expected) in state.native.inverse_plastic_cauchy_green.iter().zip([
        (-2.0 * next.creep).exp(),
        next.creep.exp(),
        next.creep.exp(),
        0.0,
        0.0,
        0.0,
    ]) {
        close(*actual, expected, 5e-9);
    }
    for (axis, log) in [next.axial_log, next.lateral_log, next.lateral_log]
        .into_iter()
        .enumerate()
    {
        close(
            step.solid.displacement_m()[axis + 1][axis],
            log.exp() - 1.0,
            3e-9,
        );
    }
    close(
        step.solid.element_first_piola_pa()[0][0],
        next.q / next.axial_log.exp(),
        3e-6,
    );
    for component in 1..9 {
        close(step.solid.element_first_piola_pa()[0][component], 0.0, 3e-6);
    }
    close(
        step.energy.stored_before_j,
        0.5 * VOLUME * COMPLIANCE * old.q * old.q,
        1e-7,
    );
    close(
        step.energy.stored_after_j,
        0.5 * VOLUME * COMPLIANCE * next.q * next.q,
        1e-7,
    );
    close(
        step.energy.inelastic_dissipation_j(),
        VOLUME * next.dissipation,
        1e-7,
    );
    close(
        step.energy.flow_potential_j,
        VOLUME * next.flow_potential,
        1e-7,
    );
    close(
        step.energy.incremental_potential_j,
        step.energy.stored_after_j + VOLUME * next.flow_potential,
        1e-7,
    );
    close(step.energy.equilibrium_work_residual_j, 0.0, 1e-7);
    let work =
        VOLUME * next.q / next.axial_log.exp() * (next.axial_log.exp() - old.axial_log.exp());
    close(step.energy.dead_load_endpoint_work_j, work, 1e-7);
    close(step.energy.support_endpoint_work_j, 0.0, 1e-9);
}

#[test]
fn loaded_tetrahedron_recovers_resistance_during_dwell_and_creeps_faster_on_reload() {
    with_cx(|_, cx| {
        let material = recovering(DT);
        let supports = supports();
        let loads = [(3, VOLUME * PIOLA)];
        let mut solid = initial(&material, cx);
        assert_eq!(solid.duration_s(), Some(DT));
        let mut point = Point::default();
        for _ in 0..3 {
            let expected = loaded(point);
            let step = solid
                .advance(&supports, &loads, [0.0; 3], settings(), controls(), cx)
                .unwrap();
            check(&step, point, expected);
            assert!(step.energy.inelastic_dissipation_j() > step.energy.flow_potential_j);
            solid = step.solid;
            point = expected;
        }
        let accepted = solid.element_states().to_vec();
        let immediate = solid
            .advance(&supports, &loads, [0.0; 3], settings(), controls(), cx)
            .unwrap();
        check(&immediate, point, loaded(point));
        assert_eq!(solid.element_states(), accepted);
        let permanent_before_dwell = point.creep;
        for _ in 0..5 {
            let expected = unloaded(point);
            let step = solid
                .advance(&supports, &[], [0.0; 3], settings(), controls(), cx)
                .unwrap();
            check(&step, point, expected);
            solid = step.solid;
            point = expected;
        }
        close(
            solid.element_states()[0].native.alpha,
            permanent_before_dwell,
            2e-9,
        );
        assert!(solid.element_states()[0].resistance_strain < accepted[0].resistance_strain);
        let reloaded = solid
            .advance(&supports, &loads, [0.0; 3], settings(), controls(), cx)
            .unwrap();
        check(&reloaded, point, loaded(point));
        assert!(reloaded.solid.displacement_m()[1][0] > immediate.solid.displacement_m()[1][0]);
        assert!(
            reloaded.solid.element_states()[0].native.alpha
                > immediate.solid.element_states()[0].native.alpha
        );
    });
}

#[test]
fn stress_free_recovering_body_translates_with_original_mass_and_one_physical_clock() {
    with_cx(|gate, cx| {
        let material = recovering(DT);
        let supports = supports();
        let loaded = initial(&material, cx)
            .advance(
                &supports,
                &[(3, VOLUME * PIOLA)],
                [0.0; 3],
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let rested = loaded
            .solid
            .advance(&supports, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        let before = rested.solid.element_states()[0].clone();
        let displacement = rested.solid.displacement_m().to_vec();
        let velocity = [0.02, -0.01, 0.005];
        let motion =
            FinitePlasticTetMotion::try_new(rested.solid, &[6.0], &[velocity; 4], cx).unwrap();
        let controls = FinitePlasticMotionControls {
            material: controls(),
            maximum_numerical_kinetic_loss_j: 1e-8,
        };
        assert!(matches!(
            motion.advance(&[], &[], [0.0; 3], DT / 2.0, settings(), controls, cx),
            Err(FinitePlasticTetError::Native(
                HyperTetError::InvalidInput { .. }
            ))
        ));
        assert_eq!(motion.solid().element_states()[0], before);
        let step = motion
            .advance(&[], &[], [0.0; 3], DT, settings(), controls, cx)
            .unwrap();
        for node in 0..4 {
            for axis in 0..3 {
                close(
                    step.motion.solid().displacement_m()[node][axis],
                    displacement[node][axis] + DT * velocity[axis],
                    2e-8,
                );
                close(step.motion.velocity_mps()[node][axis], velocity[axis], 1e-7);
            }
        }
        let after = &step.motion.solid().element_states()[0];
        close(after.native.alpha, before.native.alpha, 2e-9);
        close(
            after.resistance_strain,
            before.resistance_strain / (1.0 + RECOVERY * DT),
            2e-9,
        );
        close(step.motion.reference_mass_kg(), 6.0 * VOLUME, 1e-12);
        close(
            step.inertia.kinetic_after_j,
            0.5 * 6.0 * VOLUME * velocity.iter().map(|v| v * v).sum::<f64>(),
            1e-9,
        );
        close(step.inertia.numerical_loss_j, 0.0, 1e-9);
        close(step.material_energy.inelastic_dissipation_j(), 0.0, 1e-8);
        close(step.balance_error_j, 0.0, 1e-7);
        gate.request();
        assert!(
            motion
                .advance(&[], &[], [0.0; 3], DT, settings(), controls, cx)
                .is_err()
        );
        assert_eq!(motion.solid().element_states()[0], before);
    });
}

#[test]
fn simultaneous_finite_creep_and_yielding_keep_separate_paths_heat_and_flow_potential() {
    with_cx(|_, cx| {
        let yield_stress = 300.0;
        let hardening = 8000.0;
        let native = FiniteJ2::new(
            FiniteJ2Parameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: BULK,
                yield_stress_pa: yield_stress,
                hardening_modulus_pa: hardening,
            },
            limits(),
        )
        .unwrap();
        let material = FiniteNortonJ2::new(native, kinetics(), DT, 3.0).unwrap();
        let axial: f64 = 0.15;
        let mut supports = supports();
        supports.push((3, axial.exp() - 1.0));
        let step = initial(&material, cx)
            .advance(&supports, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        // Independent free-sided n=1 endpoint on the active plastic surface.
        let q = (axial + yield_stress / hardening) / (COMPLIANCE + DT / ETA + 1.0 / hardening);
        let plastic = (q - yield_stress) / hardening;
        let creep = DT * q / ETA;
        let lateral = q / (9.0 * BULK) - q / (6.0 * MU) - 0.5 * (plastic + creep);
        let state = &step.solid.element_states()[0];
        close(state.alpha, plastic, 2e-9);
        close(state.creep_alpha, creep, 2e-9);
        close(step.solid.displacement_m()[2][1], lateral.exp() - 1.0, 3e-9);
        close(step.solid.displacement_m()[3][2], lateral.exp() - 1.0, 3e-9);
        close(
            step.solid.element_first_piola_pa()[0][0],
            q / axial.exp(),
            3e-6,
        );
        close(
            step.energy.plastic_dissipation_j,
            VOLUME * yield_stress * plastic,
            1e-7,
        );
        close(step.energy.creep_dissipation_j, VOLUME * q * creep, 1e-7);
        close(
            step.energy.flow_potential_j,
            VOLUME * (yield_stress * plastic + 0.5 * q * creep),
            1e-7,
        );
        close(
            step.energy.stored_after_j,
            VOLUME * (0.5 * COMPLIANCE * q * q + 0.5 * hardening * plastic * plastic),
            1e-7,
        );
        assert!(step.energy.inelastic_dissipation_j() > step.energy.flow_potential_j);
        close(step.energy.equilibrium_work_residual_j, 0.0, 1e-7);
    });
}

#[test]
fn one_body_cannot_mix_different_material_durations() {
    with_cx(|_, cx| {
        let nodes = [NODES[0], NODES[1], NODES[2], NODES[3], [-1.0, 0.0, 0.0]];
        let tets = [[0, 1, 2, 3], [0, 4, 3, 2]];
        let materials = [recovering(DT), recovering(2.0 * DT)];
        assert!(matches!(
            FinitePlasticTetSolid::try_new(
                &nodes,
                &tets,
                FinitePlasticLawField::PerElement(&materials),
                TetAssemblyBudget::standard(),
                cx,
            ),
            Err(FinitePlasticTetError::Native(
                HyperTetError::InvalidInput { .. }
            ))
        ));
    });
}
