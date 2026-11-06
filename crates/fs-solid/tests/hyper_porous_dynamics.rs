//! Finite crushing through the original porous return and consistent momentum.
//! Independent hydrostatic and P1 integral formulas check each physical account;
//! authored coefficients and density are software fixtures, not foam calibration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::porous::source::{FINITE_POROUS_LAW_ID, ResolvedFinitePorous};
use fs_qty::Dims;
use fs_solid::hyper3::plastic::porous::dynamics::{
    FinitePlasticMotionControls, FinitePorousTetMotion,
};
use fs_solid::hyper3::plastic::porous::{FinitePorousTetError, FinitePorousTetSolid};
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetError,
};
use fs_solid::{HyperTetError, HyperTetSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const V: f64 = 1.0 / 6.0;
const RHO: f64 = 6.0;
const NODAL_RHO: [f64; 4] = [3.0, 6.0, 9.0, 12.0];
const K: f64 = 160.0;
const Y: f64 = 2.0;
const H: f64 = 8.0;
const A: f64 = 1.2;

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 98,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate);
    });
}

fn source() -> ResolvedFinitePorous {
    const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
    let digest = fs_blake3::hash_bytes(b"authored porous motion fixture");
    let card = ConstitutiveModelCard {
        law: LawId(FINITE_POROUS_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("shear_modulus_pa", 80.0, PA),
            ("bulk_modulus_pa", K, PA),
            ("yield_stress_pa", Y, PA),
            ("hardening_modulus_pa", H, PA),
            ("pressure_sensitivity", A, Dims::NONE),
            ("minimum_jacobian", 0.2, Dims::NONE),
            ("maximum_elastic_log_strain", 0.5, Dims::NONE),
            ("maximum_plastic_log_strain", 2.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 3.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
        sources: vec![digest],
        provenance: Provenance {
            source: "authored porous motion fixture".into(),
            license: "test fixture".into(),
            artifact: Some(digest),
        },
    };
    ResolvedFinitePorous::try_from_model_card_at(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
    )
    .unwrap()
}

fn solid<'a>(
    law: &'a ResolvedFinitePorous,
    cx: &Cx<'_>,
) -> FinitePorousTetSolid<'a, ResolvedFinitePorous> {
    FinitePorousTetSolid::try_new(
        &NODES,
        &TETS,
        FinitePlasticLawField::Uniform(law),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}

fn motion<'a>(
    law: &'a ResolvedFinitePorous,
    graded: bool,
    velocity: &[[f64; 3]],
    cx: &Cx<'_>,
) -> FinitePorousTetMotion<'a, ResolvedFinitePorous> {
    if graded {
        FinitePorousTetMotion::try_new_with_nodal_density(solid(law, cx), &NODAL_RHO, velocity, cx)
    } else {
        FinitePorousTetMotion::try_new(solid(law, cx), &[RHO], velocity, cx)
    }
    .unwrap()
}

fn controls() -> FinitePlasticMotionControls {
    FinitePlasticMotionControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-7,
            maximum_endpoint_work_defect_j: 10.0,
            yield_tolerance_pa: 2e-7,
        },
        maximum_numerical_kinetic_loss_j: 10.0,
    }
}

fn settings() -> HyperTetSettings {
    HyperTetSettings {
        relative_tolerance: 1e-11,
        ..HyperTetSettings::default()
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:e}"
    );
}

// Exact integrals of rho*N_a*N_b on the original unit-coordinate tetrahedron.
fn mass_tile(graded: bool) -> [[f64; 4]; 4] {
    let rho = if graded { NODAL_RHO } else { [RHO; 4] };
    let sum = rho.iter().sum::<f64>();
    core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            V * if a == b {
                (sum + 2.0 * rho[a]) / 60.0
            } else {
                (sum + rho[a] + rho[b]) / 120.0
            }
        })
    })
}

#[test]
fn unrestrained_porous_motion_uses_original_graded_mass_and_source() {
    with_cx(|cx, _| {
        let law = source();
        let initial_v = [0.2, -0.1, 0.05];
        let dt = 0.02;
        for graded in [false, true] {
            for acceleration in [[0.0; 3], [1.0, -0.5, 0.25]] {
                let mut body = motion(&law, graded, &[initial_v; 4], cx);
                let mass = if graded { 7.5 } else { RHO } * V;
                let loads = body
                    .reference_acceleration_load_n(&[acceleration; 4], cx)
                    .unwrap();
                for node in 0..4 {
                    let row_mass = if graded {
                        V * (NODAL_RHO.iter().sum::<f64>() + NODAL_RHO[node]) / 20.0
                    } else {
                        mass / 4.0
                    };
                    for axis in 0..3 {
                        let force = loads
                            .iter()
                            .find(|(dof, _)| *dof == 3 * node + axis)
                            .map_or(0.0, |(_, f)| *f);
                        close(force, row_mass * acceleration[axis], 1e-14);
                    }
                }
                for interval in 1..=3 {
                    let n = f64::from(interval);
                    let v0: [f64; 3] =
                        core::array::from_fn(|i| initial_v[i] + (n - 1.0) * dt * acceleration[i]);
                    let v: [f64; 3] =
                        core::array::from_fn(|i| initial_v[i] + n * dt * acceleration[i]);
                    let u: [f64; 3] = core::array::from_fn(|i| {
                        n * dt * initial_v[i] + 0.5 * n * (n + 1.0) * dt * dt * acceleration[i]
                    });
                    let step = body
                        .advance(&[], &loads, [0.0; 3], dt, settings(), controls(), cx)
                        .unwrap();
                    assert!(step.equilibrium.reactions_n.is_empty());
                    for node in 0..4 {
                        for axis in 0..3 {
                            close(
                                step.motion.solid().displacement_m()[node][axis],
                                u[axis],
                                1e-10,
                            );
                            close(step.motion.velocity_mps()[node][axis], v[axis], 1e-8);
                        }
                    }
                    close(step.motion.reference_mass_kg(), mass, 1e-13);
                    close(
                        step.inertia.kinetic_before_j,
                        0.5 * mass * v0.iter().map(|a| a * a).sum::<f64>(),
                        1e-9,
                    );
                    close(
                        step.inertia.kinetic_after_j,
                        0.5 * mass * v.iter().map(|a| a * a).sum::<f64>(),
                        1e-9,
                    );
                    close(
                        step.inertia.numerical_loss_j,
                        0.5 * mass * dt * dt * acceleration.iter().map(|a| a * a).sum::<f64>(),
                        1e-9,
                    );
                    close(
                        step.material_energy.dead_load_endpoint_work_j,
                        mass * dt * (0..3).map(|i| acceleration[i] * v[i]).sum::<f64>(),
                        1e-9,
                    );
                    close(step.material_energy.stored_after_j, 0.0, 1e-12);
                    close(step.balance_error_j, 0.0, 1e-9);
                    assert_eq!(step.material_energy.plastic_dissipation_j, 0.0);
                    assert_eq!(
                        step.motion.solid().element_states(),
                        &[law.law().initial_state()]
                    );
                    close(step.motion.solid().element_jacobians()[0], 1.0, 1e-11);
                    assert_eq!(step.motion.solid().element_plastic_jacobians(), &[1.0]);
                    match step.motion.solid().material() {
                        FinitePlasticLawField::Uniform(retained) => {
                            assert!(std::ptr::eq(retained, &law));
                            assert_eq!(retained.context().axes()["T"], 300.0);
                        }
                        FinitePlasticLawField::PerElement(_) => panic!("source binding changed"),
                    }
                    assert_eq!(
                        step.motion.reference_nodal_density_kg_m3(),
                        graded.then_some(&NODAL_RHO[..])
                    );
                    body = step.motion;
                }
            }
        }
    });
}

#[test]
fn inertial_hydrostatic_crush_unload_reload_and_late_refusal_keep_all_accounts() {
    with_cx(|cx, gate| {
        let law = source();
        let dt = 0.1;
        let beta = A / (1.0 + A * A / 9.0).sqrt();
        let peak = 0.8_f64;
        let peak_alpha = (-3.0 * K * beta * peak.ln() - Y) / (H + K * beta * beta);
        let release = (-beta * peak_alpha / 3.0).exp();
        let held: Vec<_> = (0..12)
            .filter(|dof| ![3, 7, 11].contains(dof))
            .map(|dof| (dof, 0.0))
            .collect();
        for graded in [false, true] {
            let mut body = motion(&law, graded, &[[0.0; 3]; 4], cx);
            let mass = mass_tile(graded);
            let axial_mass = mass[1][1] + mass[2][2] + mass[3][3];
            let total_mass = mass.iter().flatten().sum::<f64>();
            let (mut old_u, mut old_v, mut old_alpha, mut log_jp, mut old_storage) =
                (0.0, 0.0, 0.0, 0.0, 0.0);
            for (interval, stretch) in [0.9, peak, release, (peak * release).sqrt()]
                .into_iter()
                .enumerate()
            {
                let trial_tau = K * (3.0 * stretch.ln() - log_jp);
                let delta_alpha =
                    (beta * trial_tau.abs() - Y - H * old_alpha).max(0.0) / (H + K * beta * beta);
                let alpha = old_alpha + delta_alpha;
                log_jp += beta * trial_tau.signum() * delta_alpha;
                let tau = K * (3.0 * stretch.ln() - log_jp);
                let piola = tau / stretch;
                let u = stretch - 1.0;
                let v = (u - old_u) / dt;
                let acceleration = (v - old_v) / dt;
                let forces: Vec<_> = [3, 7, 11]
                    .into_iter()
                    .enumerate()
                    .map(|(axis, dof)| (dof, V * piola + mass[axis + 1][axis + 1] * acceleration))
                    .collect();
                let before_states = body.solid().element_states().to_vec();
                let before_u = body.solid().displacement_m().to_vec();
                let before_v = body.velocity_mps().to_vec();
                let before_jp = body.solid().element_plastic_jacobians().to_vec();
                let step = body
                    .advance(&held, &forces, [0.0; 3], dt, settings(), controls(), cx)
                    .unwrap();
                assert_eq!(body.solid().element_states(), before_states);
                assert_eq!(body.solid().displacement_m(), before_u);
                assert_eq!(body.velocity_mps(), before_v);
                if interval == 2 {
                    let mut refuse = controls();
                    refuse.maximum_numerical_kinetic_loss_j = 0.0;
                    let error = body
                        .advance(&held, &forces, [0.0; 3], dt, settings(), refuse, cx)
                        .unwrap_err();
                    assert!(matches!(
                        error,
                        FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate {
                            quantity: "porous backward-Euler numerical kinetic loss",
                            ..
                        })
                    ));
                    assert_eq!(body.solid().element_states(), before_states);
                    assert_eq!(body.solid().displacement_m(), before_u);
                    assert_eq!(body.velocity_mps(), before_v);
                    assert_eq!(body.solid().element_plastic_jacobians(), before_jp);
                    let retry = body
                        .advance(&held, &forces, [0.0; 3], dt, settings(), controls(), cx)
                        .unwrap();
                    assert_eq!(retry.equilibrium, step.equilibrium);
                    assert_eq!(retry.material_energy, step.material_energy);
                    assert_eq!(retry.inertia, step.inertia);
                    assert_eq!(
                        retry.motion.solid().element_states(),
                        step.motion.solid().element_states()
                    );
                    assert_eq!(retry.motion.velocity_mps(), step.motion.velocity_mps());
                }
                let storage = V * (0.5 * tau * tau / K + 0.5 * H * alpha * alpha);
                let heat = V * Y * delta_alpha;
                let internal_work = 0.5 * piola * (u - old_u);
                let inertia_work = axial_mass * acceleration * (u - old_u);
                close(step.material_energy.stored_before_j, old_storage, 1e-8);
                close(step.material_energy.stored_after_j, storage, 1e-8);
                close(step.material_energy.plastic_dissipation_j, heat, 1e-8);
                close(
                    step.material_energy.internal_endpoint_work_j,
                    internal_work,
                    1e-8,
                );
                close(
                    step.material_energy.dead_load_endpoint_work_j,
                    internal_work + inertia_work,
                    1e-8,
                );
                close(step.material_energy.support_endpoint_work_j, 0.0, 0.0);
                close(
                    step.material_energy.endpoint_work_defect_j,
                    internal_work - storage + old_storage - heat,
                    1e-8,
                );
                close(
                    step.inertia.kinetic_before_j,
                    0.5 * axial_mass * old_v * old_v,
                    1e-8,
                );
                close(step.inertia.kinetic_after_j, 0.5 * axial_mass * v * v, 1e-8);
                close(
                    step.inertia.numerical_loss_j,
                    0.5 * axial_mass * (v - old_v).powi(2),
                    1e-8,
                );
                close(step.inertia.inertial_work_j, inertia_work, 1e-8);
                close(step.balance_error_j, 0.0, 1e-8);
                close(step.motion.solid().element_states()[0].alpha, alpha, 1e-9);
                close(
                    step.motion.solid().element_plastic_jacobians()[0],
                    log_jp.exp(),
                    1e-9,
                );
                close(
                    step.motion.solid().element_jacobians()[0],
                    stretch.powi(3),
                    1e-9,
                );
                close(
                    step.motion.current_density_kg_m3(cx).unwrap()[0] * stretch.powi(3) * V,
                    total_mass,
                    1e-8,
                );
                close(step.motion.reference_mass_kg(), total_mass, 1e-14);
                for (node, reference) in NODES.iter().enumerate() {
                    for axis in 0..3 {
                        close(
                            step.motion.solid().displacement_m()[node][axis],
                            u * reference[axis],
                            1e-9,
                        );
                        close(
                            step.motion.velocity_mps()[node][axis],
                            v * reference[axis],
                            1e-8,
                        );
                    }
                }
                for component in 0..9 {
                    close(
                        step.motion.solid().element_deformation()[0][component],
                        if component % 4 == 0 { stretch } else { 0.0 },
                        1e-9,
                    );
                    close(
                        step.motion.solid().element_first_piola_pa()[0][component],
                        if component % 4 == 0 { piola } else { 0.0 },
                        1e-8,
                    );
                }
                for &(dof, reaction) in &step.equilibrium.reactions_n {
                    let node = dof / 3;
                    let axis = dof % 3;
                    let internal = if node == 0 { -V * piola } else { 0.0 };
                    close(
                        reaction,
                        internal + mass[node][axis + 1] * acceleration,
                        1e-8,
                    );
                }
                if interval >= 2 {
                    assert_eq!(step.motion.solid().element_states(), before_states);
                    assert_eq!(step.material_energy.plastic_dissipation_j, 0.0);
                }
                body = step.motion;
                (old_u, old_v, old_alpha, old_storage) = (u, v, alpha, storage);
            }
            assert!(body.solid().element_plastic_jacobians()[0] < 0.55);
            for invalid_dt in [0.0, -1.0, f64::NAN] {
                assert!(
                    body.advance(&held, &[], [0.0; 3], invalid_dt, settings(), controls(), cx)
                        .is_err()
                );
            }
            if graded {
                let saved = body.solid().element_states().to_vec();
                let saved_velocity = body.velocity_mps().to_vec();
                gate.request();
                assert!(matches!(
                    body.advance(&held, &[], [0.0; 3], dt, settings(), controls(), cx),
                    Err(FinitePorousTetError::Solid(FinitePlasticTetError::Native(
                        HyperTetError::Cancelled
                    )))
                ));
                assert!(body.current_density_kg_m3(cx).is_err());
                assert_eq!(body.solid().element_states(), saved);
                assert_eq!(body.velocity_mps(), saved_velocity);
            }
        }
    });
}

#[test]
fn porous_motion_rejects_incomplete_mass_velocity_and_nonfinite_acceleration() {
    with_cx(|cx, _| {
        let law = source();
        for density in [&[][..], &[0.0][..], &[f64::INFINITY][..], &[RHO, RHO][..]] {
            assert!(
                FinitePorousTetMotion::try_new(solid(&law, cx), density, &[[0.0; 3]; 4], cx)
                    .is_err()
            );
        }
        assert!(
            FinitePorousTetMotion::try_new(solid(&law, cx), &[RHO], &[[0.0; 3]; 3], cx).is_err()
        );
        assert!(
            FinitePorousTetMotion::try_new_with_nodal_density(
                solid(&law, cx),
                &NODAL_RHO[..3],
                &[[0.0; 3]; 4],
                cx
            )
            .is_err()
        );
        let body = motion(&law, false, &[[0.0; 3]; 4], cx);
        assert!(
            body.reference_acceleration_load_n(&[[0.0; 3]; 3], cx)
                .is_err()
        );
        assert!(
            body.reference_acceleration_load_n(&[[f64::NAN; 3]; 4], cx)
                .is_err()
        );
    });
}
