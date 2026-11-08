//! G1/G3: original sheet inventory, free momentum, finite ductile flow and
//! separate material/kinetic work under the same physical interval.
#![deny(unsafe_code)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::ductile::evolution::DuctileJ2Settings;
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::FiniteDuctileJ2;
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::ductile::{DuctileStepConditions, DuctileTemperatures};
use fs_solid::hyper2d::plastic::SheetAssemblyBudget;
use fs_solid::hyper2d::plastic::ductile::dynamics::{
    FiniteDuctileSheetMotion, FiniteDuctileSheetMotionControls,
};
use fs_solid::hyper2d::plastic::ductile::{
    FiniteDuctileSheet, FiniteDuctileSheetError, SheetCrackBandGeometry,
};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::plastic3::SmallStrainTetError;
use fs_solid::{IncrementSettings, Mesh2};

const MU: f64 = 8000.0;
const K: f64 = 16000.0;
const Y: f64 = 20.0;
const H: f64 = 80.0;

fn with_cx(run: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 97,
                    kernel_id: 217,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn model(width: f64, viscosity: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    let native = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.6,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap();
    let locus = TabulatedDuctileLocus::new(
        vec![(-50.0, 0.005), (50.0, 0.005)],
        DuctileDomain {
            triaxiality: [-50.0, 50.0],
            strain_rate_s_inv: [1e-12, 1000.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    let law = FiniteDuctileJ2::new(
        native,
        locus,
        DuctileJ2Settings {
            characteristic_length_m: width,
            post_initiation_plastic_work_j_m2: 2.0,
            maximum_damage_increment: 1.0,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-8,
                yield_pa: 1e-7,
            },
        },
    )
    .unwrap();
    if viscosity == 0.0 {
        law
    } else {
        law.with_overstress(FiniteOverstressParameters {
            viscosity_pa_s: viscosity,
            stress_scale_pa: Y,
            exponent: 1.0,
        })
        .unwrap()
    }
}

fn bands(mesh: &Mesh2, cx: &Cx<'_>) -> Vec<SheetCrackBandGeometry> {
    (0..mesh.elems.len())
        .map(|e| SheetCrackBandGeometry::from_element(mesh, e, [1.0, 0.0], cx).unwrap())
        .collect()
}

fn sheet<'g, 'm>(
    mesh: &'g Mesh2,
    materials: &'m [&'m FiniteDuctileJ2<TabulatedDuctileLocus>],
    thickness: &'g [f64],
    bands: &'g [SheetCrackBandGeometry],
    cx: &Cx<'_>,
) -> FiniteDuctileSheet<'g, 'm, TabulatedDuctileLocus> {
    FiniteDuctileSheet::try_new(
        mesh,
        materials,
        thickness,
        bands,
        PlaneStressSettings {
            relative_tolerance: 1e-12,
            ..PlaneStressSettings::default()
        },
        SheetAssemblyBudget {
            maximum_nodes: 24,
            maximum_elements: 24,
            maximum_free_dofs: 48,
            minimum_scaled_jacobian: 1e-4,
        },
        cx,
    )
    .unwrap()
}

fn conditions(duration_s: f64) -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    }
}

fn controls() -> FiniteDuctileSheetMotionControls {
    FiniteDuctileSheetMotionControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-7,
            maximum_endpoint_work_defect_j: 10.0,
            yield_tolerance_pa: 2e-6,
        },
        maximum_numerical_kinetic_loss_j: 1000.0,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 + 2e-8 * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn centroid(mesh: &Mesh2, element: usize) -> [f64; 2] {
    let nodes = &mesh.elems[element];
    [0, 1].map(|i| nodes.iter().map(|&a| mesh.nodes[a][i]).sum::<f64>() / nodes.len() as f64)
}

#[test]
fn heterogeneous_thickness_and_density_keep_original_mass_moment_and_free_acceleration() {
    with_cx(|_, cx| {
        for mesh in [
            Mesh2::triangles(2.0, 1.0, 2, 1),
            Mesh2::quads(2.0, 1.0, 2, 1),
        ] {
            let bands = bands(&mesh, cx);
            let models: Vec<_> = bands.iter().map(|b| model(b.band_width_m(), 0.0)).collect();
            let materials: Vec<_> = models.iter().collect();
            let thickness: Vec<_> = (0..mesh.elems.len())
                .map(|e| 0.001 * (e + 1) as f64)
                .collect();
            let density: Vec<_> = (0..mesh.elems.len())
                .map(|e| 2700.0 + 500.0 * e as f64)
                .collect();
            let initial = sheet(&mesh, &materials, &thickness, &bands, cx);
            let velocity = vec![[0.2, 0.1]; mesh.nodes.len()];
            let motion =
                FiniteDuctileSheetMotion::try_new(initial, &density, &velocity, 1 << 20, cx)
                    .unwrap();
            let acceleration = [0.4, -0.3];
            let forces = motion
                .reference_acceleration_load_n(&vec![acceleration; mesh.nodes.len()], cx)
                .unwrap();
            let mut expected_mass = 0.0;
            let mut expected_moment = 0.0;
            for e in 0..mesh.elems.len() {
                let mass = bands[e].reference_area_m2() * thickness[e] * density[e];
                let center = centroid(&mesh, e);
                expected_mass += mass;
                expected_moment +=
                    mass * (center[0] * acceleration[1] - center[1] * acceleration[0]);
            }
            near(motion.reference_mass_kg(), expected_mass);
            let mut total = [0.0; 2];
            let mut moment = 0.0;
            for &(dof, force) in &forces {
                total[dof % 2] += force;
                moment += force
                    * if dof % 2 == 1 {
                        mesh.nodes[dof / 2][0]
                    } else {
                        -mesh.nodes[dof / 2][1]
                    };
            }
            for i in 0..2 {
                near(total[i], expected_mass * acceleration[i]);
            }
            near(moment, expected_moment);
            let dt = 0.07;
            let step = motion
                .advance(
                    &[],
                    &forces,
                    [0.0; 2],
                    conditions(dt),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let next = [0, 1].map(|i| velocity[0][i] + dt * acceleration[i]);
            for (u, v) in step
                .motion
                .sheet()
                .displacement_m()
                .iter()
                .zip(step.motion.velocity_mps())
            {
                for i in 0..2 {
                    near(u[i], dt * next[i]);
                    near(v[i], next[i]);
                }
            }
            near(
                step.inertia.kinetic_before_j,
                0.5 * expected_mass * (0.2_f64.powi(2) + 0.1_f64.powi(2)),
            );
            near(
                step.inertia.kinetic_after_j,
                0.5 * expected_mass * (next[0].powi(2) + next[1].powi(2)),
            );
            near(
                step.inertia.numerical_loss_j,
                0.5 * expected_mass * dt * dt * (acceleration[0].powi(2) + acceleration[1].powi(2)),
            );
            near(
                step.inertia.inertial_work_j,
                expected_mass * dt * (acceleration[0] * next[0] + acceleration[1] * next[1]),
            );
            near(step.balance_error_j, 0.0);
            assert_eq!(step.motion.reference_density_kg_m3(), density);
            assert_eq!(
                step.motion.reference_mass_kg().to_bits(),
                motion.reference_mass_kg().to_bits()
            );
            assert!(
                step.motion
                    .sheet()
                    .quadrature_states()
                    .iter()
                    .all(|s| s.damage() == 0.0 && s.native_state().alpha == 0.0)
            );
            assert_eq!(motion.velocity_mps(), velocity);
            assert_eq!(motion.sheet().physical_time_s(), 0.0);
        }
    });
}

#[test]
fn prescribed_rate_dependent_damage_and_held_relaxation_keep_material_and_kinetic_accounts_separate()
 {
    with_cx(|_, cx| {
        let dt = 0.2;
        let e = 0.03_f64;
        let viscosity = 600.0;
        for mesh in [
            Mesh2::triangles(1.0, 1.0, 1, 1),
            Mesh2::quads(1.0, 1.0, 1, 1),
        ] {
            let bands = bands(&mesh, cx);
            let models: Vec<_> = bands
                .iter()
                .map(|b| model(b.band_width_m(), viscosity))
                .collect();
            let materials: Vec<_> = models.iter().collect();
            let thickness: Vec<_> = (0..mesh.elems.len())
                .map(|e| 0.001 * (e + 1) as f64)
                .collect();
            let density = vec![7800.0; mesh.elems.len()];
            let initial = sheet(&mesh, &materials, &thickness, &bands, cx);
            let u: Vec<_> = mesh
                .nodes
                .iter()
                .map(|x| x.map(|v| e.exp_m1() * v))
                .collect();
            let prescribed: Vec<_> = u.as_flattened().iter().copied().enumerate().collect();
            let static_step = initial
                .advance(
                    &prescribed,
                    &[],
                    conditions(dt),
                    IncrementSettings::default(),
                    controls().material,
                    cx,
                )
                .unwrap();
            let motion = FiniteDuctileSheetMotion::try_new(
                initial,
                &density,
                &vec![[0.0; 2]; mesh.nodes.len()],
                1 << 20,
                cx,
            )
            .unwrap();
            let body = [3.0, -2.0];
            let step = motion
                .advance(
                    &prescribed,
                    &[],
                    body,
                    conditions(dt),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(step.point_energy, static_step.point_energy);
            assert_eq!(
                step.motion.sheet().quadrature_states(),
                static_step.sheet.quadrature_states()
            );
            let mass = motion.reference_mass_kg();
            // Exact unit-square moment, also equal on both constituent P1
            // triangles: mean(x²+y²)=2/3. Centroid-only mass would fail this.
            let kinetic = mass * (e.exp_m1() / dt).powi(2) / 3.0;
            near(step.inertia.kinetic_after_j, kinetic);
            near(step.inertia.numerical_loss_j, kinetic);
            near(step.inertia.inertial_work_j, 2.0 * kinetic);
            let mut body_work = 0.0;
            for element in 0..mesh.elems.len() {
                let center = centroid(&mesh, element);
                body_work += bands[element].reference_area_m2()
                    * thickness[element]
                    * e.exp_m1()
                    * (body[0] * center[0] + body[1] * center[1]);
            }
            near(step.material_energy.dead_load_endpoint_work_j, body_work);
            near(
                step.material_energy.support_endpoint_work_j + body_work,
                static_step.energy.support_endpoint_work_j + 2.0 * kinetic,
            );
            assert!(step.material_energy.plastic_dissipation_j > 0.0);
            assert!(
                step.motion
                    .sheet()
                    .quadrature_states()
                    .iter()
                    .all(|s| s.damage() > 0.0)
            );
            let c = 1.0 / (6.0 * MU) + 2.0 / (9.0 * K);
            let alpha = (e - c * Y) / (0.5 + c * (H + viscosity / dt));
            for state in step.motion.sheet().quadrature_states() {
                near(state.native_state().alpha, alpha);
            }
            let relaxed = step
                .motion
                .advance(
                    &prescribed,
                    &[],
                    [0.0; 2],
                    conditions(dt),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let next_alpha =
                (e - c * Y + c * viscosity * alpha / dt) / (0.5 + c * (H + viscosity / dt));
            for (old, new) in step
                .motion
                .sheet()
                .quadrature_states()
                .iter()
                .zip(relaxed.motion.sheet().quadrature_states())
            {
                near(new.native_state().alpha, next_alpha);
                assert!(new.damage() > old.damage());
            }
            near(relaxed.inertia.kinetic_change_j, -kinetic);
            near(relaxed.inertia.numerical_loss_j, kinetic);
            assert_eq!(relaxed.inertia.kinetic_after_j, 0.0);
            assert_eq!(relaxed.inertia.inertial_work_j, 0.0);
            assert!(relaxed.material_energy.plastic_dissipation_j > 0.0);
            assert!(relaxed.material_energy.damage_release_j > 0.0);
            near(relaxed.balance_error_j, 0.0);
            assert_eq!(
                relaxed.motion.reference_mass_kg().to_bits(),
                motion.reference_mass_kg().to_bits()
            );
            assert_eq!(relaxed.motion.sheet().physical_time_s(), 2.0 * dt);
            assert!(
                relaxed
                    .motion
                    .sheet()
                    .quadrature_thickness_m()
                    .iter()
                    .zip(&thickness)
                    .any(|(current, original)| current != original)
            );
        }
    });
}

#[test]
fn moving_free_deformation_dof_solves_with_damage_and_refusals_preserve_original_motion() {
    with_cx(|gate, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let bands = bands(&mesh, cx);
        let law = model(bands[0].band_width_m(), 0.0);
        let materials = [&law];
        let thickness = [0.001];
        let initial = sheet(&mesh, &materials, &thickness, &bands, cx);
        let zero = vec![[0.0; 2]; mesh.nodes.len()];
        assert!(
            FiniteDuctileSheetMotion::try_new(initial.clone(), &[], &zero, 1 << 20, cx).is_err()
        );
        assert!(
            FiniteDuctileSheetMotion::try_new(initial.clone(), &[-1.0], &zero, 1 << 20, cx)
                .is_err()
        );
        assert!(
            FiniteDuctileSheetMotion::try_new(initial.clone(), &[7800.0], &zero, 0, cx).is_err()
        );
        assert!(
            FiniteDuctileSheetMotion::try_new(initial.clone(), &[7800.0], &zero[..1], 1 << 20, cx)
                .is_err()
        );
        let motion =
            FiniteDuctileSheetMotion::try_new(initial, &[7800.0], &zero, 1 << 20, cx).unwrap();
        let dt = 0.002;
        let u: Vec<_> = mesh
            .nodes
            .iter()
            .map(|x| x.map(|v| 0.02_f64.exp_m1() * v))
            .collect();
        let free_node = mesh.nodes.iter().position(|p| *p == [1.0, 0.0]).unwrap();
        let free = 2 * free_node;
        let supports: Vec<_> = u
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .filter(|(d, _)| *d != free)
            .collect();
        let evaluation = motion
            .sheet()
            .evaluate_trial(&u, conditions(dt), true, cx)
            .unwrap();
        let acceleration: Vec<_> = u.iter().map(|v| v.map(|x| x / (dt * dt))).collect();
        let inertia = motion
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        let force = evaluation.internal_force_n[free_node][0]
            + inertia.iter().find(|(d, _)| *d == free).unwrap().1;
        let step = motion
            .advance(
                &supports,
                &[(free, force)],
                [0.0; 2],
                conditions(dt),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        near(
            step.motion.sheet().displacement_m()[free_node][0],
            u[free_node][0],
        );
        near(
            step.motion.velocity_mps()[free_node][0],
            u[free_node][0] / dt,
        );
        assert!(
            step.motion
                .sheet()
                .quadrature_states()
                .iter()
                .all(|s| s.damage() > 0.0)
        );
        assert!(step.material_energy.plastic_dissipation_j > 0.0);
        near(step.balance_error_j, 0.0);
        assert!(matches!(
            motion.advance(
                &supports,
                &[(free, force)],
                [0.0; 2],
                conditions(dt),
                IncrementSettings::default(),
                FiniteDuctileSheetMotionControls {
                    maximum_numerical_kinetic_loss_j: 0.0,
                    ..controls()
                },
                cx
            ),
            Err(FiniteDuctileSheetError::Energy {
                quantity: "backward-Euler numerical kinetic loss",
                ..
            })
        ));
        assert!(
            motion
                .advance(
                    &supports,
                    &[(free, force)],
                    [0.0; 2],
                    conditions(0.0),
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(motion.velocity_mps(), zero);
        assert_eq!(motion.sheet().displacement_m(), zero);
        assert!(
            motion
                .sheet()
                .quadrature_states()
                .iter()
                .all(|s| s.damage() == 0.0 && s.native_state().alpha == 0.0)
        );
        gate.request();
        assert!(matches!(
            motion.advance(
                &supports,
                &[(free, force)],
                [0.0; 2],
                conditions(dt),
                IncrementSettings::default(),
                controls(),
                cx
            ),
            Err(FiniteDuctileSheetError::Solid(
                SmallStrainTetError::Cancelled
            ))
        ));
        assert_eq!(motion.sheet().physical_time_s(), 0.0);
    });
}
