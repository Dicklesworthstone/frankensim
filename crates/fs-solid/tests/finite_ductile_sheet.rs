//! G1/G3: original P1/Q1 sheets through progressive finite ductile softening.
//! Synthetic coefficients verify mechanics and history; no alloy is qualified.
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
const ONSET: f64 = 0.005;
const GP: f64 = 2.0;
const THICKNESS: f64 = 0.001;
const DT: f64 = 0.2;
const C: f64 = 1.0 / (6.0 * MU) + 2.0 / (9.0 * K);

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 59,
                kernel_id: 197,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn model(width: f64, work: f64, damage_limit: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
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
        vec![(-2.0, ONSET), (2.0, ONSET)],
        DuctileDomain {
            triaxiality: [-2.0, 2.0],
            strain_rate_s_inv: [1e-12, 100.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    FiniteDuctileJ2::new(
        native,
        locus,
        DuctileJ2Settings {
            characteristic_length_m: width,
            post_initiation_plastic_work_j_m2: work,
            maximum_damage_increment: damage_limit,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-8,
                yield_pa: 1e-7,
            },
        },
    )
    .unwrap()
}

fn bands(mesh: &Mesh2, cx: &Cx<'_>) -> Vec<SheetCrackBandGeometry> {
    (0..mesh.elems.len())
        .map(|e| SheetCrackBandGeometry::from_element(mesh, e, [1.0, 0.0], cx).unwrap())
        .collect()
}

fn budget() -> SheetAssemblyBudget {
    SheetAssemblyBudget {
        maximum_nodes: 20,
        maximum_elements: 20,
        maximum_free_dofs: 40,
        minimum_scaled_jacobian: 1e-4,
    }
}

fn local() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 1e-12,
        ..PlaneStressSettings::default()
    }
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-8,
        maximum_endpoint_work_defect_j: 100.0,
        yield_tolerance_pa: 2e-6,
    }
}

fn conditions(duration_s: f64) -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    }
}

fn displacement(mesh: &Mesh2, e: f64, angle: f64) -> Vec<[f64; 2]> {
    let stretch = e.exp();
    let (s, c) = angle.sin_cos();
    mesh.nodes
        .iter()
        .map(|&[x, y]| [stretch * (c * x - s * y) - x, stretch * (s * x + c * y) - y])
        .collect()
}

fn prescribed(u: &[[f64; 2]]) -> Vec<(usize, f64)> {
    u.as_flattened().iter().copied().enumerate().collect()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 + 3e-8 * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

// Independent Hencky plane-stress relations in equibiaxial tension:
// e = C*q + alpha/2, q = Y + H*alpha + viscosity*(alpha-old_alpha)/dt.
fn equibiaxial_alpha(e: f64, old_alpha: f64, viscosity: f64, dt: f64) -> f64 {
    (e - C * Y + C * viscosity * old_alpha / dt) / (0.5 + C * H + C * viscosity / dt)
}

fn integrity(alpha: f64, width: f64, work: f64) -> f64 {
    let post = (alpha - ONSET).max(0.0);
    let resistance_integral = post * (Y + H * (ONSET + 0.5 * post));
    (-width * resistance_integral / work).exp()
}

#[test]
fn equibiaxial_softening_and_rotation_match_stress_thickness_and_physical_work() {
    with_cx(|_, cx| {
        let e = 0.02;
        let alpha = equibiaxial_alpha(e, 0.0, 0.0, DT);
        let q = Y + H * alpha;
        let stretch = e.exp();
        let transverse = (2.0 * q / (3.0 * K) - 2.0 * e).exp();
        let native_storage = C * q * q + 0.5 * H * alpha * alpha;
        for mesh in [
            Mesh2::triangles(1.0, 1.0, 1, 1),
            Mesh2::quads(1.0, 1.0, 1, 1),
        ] {
            let bands = bands(&mesh, cx);
            let models: Vec<_> = bands
                .iter()
                .map(|b| model(b.band_width_m(), GP, 1.0))
                .collect();
            let materials: Vec<_> = models.iter().collect();
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let initial = FiniteDuctileSheet::try_new(
                &mesh,
                &materials,
                &thickness,
                &bands,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let original = initial.quadrature_states().to_vec();
            let points_per_element = if mesh.elems[0].len() == 3 { 1 } else { 4 };
            assert_eq!(original.len(), points_per_element * mesh.elems.len());
            for angle in [0.0_f64, 0.61] {
                let supports = prescribed(&displacement(&mesh, e, angle));
                let step = initial
                    .advance(
                        &supports,
                        &[],
                        conditions(DT),
                        IncrementSettings::default(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                let (s, c) = angle.sin_cos();
                let mut storage = 0.0;
                let mut plastic = 0.0;
                let mut endpoint_work = 0.0;
                for (element, band) in bands.iter().enumerate() {
                    let g = integrity(alpha, band.band_width_m(), GP);
                    let volume = band.reference_area_m2() * THICKNESS;
                    near(step.sheet.reference_volume_m3(element).unwrap(), volume);
                    storage += volume * g * native_storage;
                    plastic += volume * g * Y * alpha;
                    endpoint_work += volume * 2.0 * g * q * (1.0 - c / stretch);
                    for state in &step.sheet.quadrature_states()
                        [element * points_per_element..(element + 1) * points_per_element]
                    {
                        near(state.native_state().alpha, alpha);
                        near(state.integrity(), g);
                        assert!(state.damage() > 0.1);
                        near(state.deformation()[8], transverse);
                        for (i, factor) in [(0, c), (1, -s), (3, s), (4, c)] {
                            near(state.first_piola_pa()[i], g * q * factor / stretch);
                        }
                        for i in [2, 5, 8] {
                            assert!(state.first_piola_pa()[i].abs() < 1e-6);
                        }
                        assert!((state.deformation()[8] - (-2.0 * e).exp()).abs() > 1e-5);
                    }
                }
                for current_thickness in step.sheet.quadrature_thickness_m() {
                    near(current_thickness, THICKNESS * transverse);
                }
                near(step.sheet.stored_energy_j(), storage);
                near(step.energy.stored_after_j, storage);
                near(step.energy.plastic_dissipation_j, plastic);
                near(step.energy.damage_release_j, 0.0);
                near(step.energy.internal_endpoint_work_j, endpoint_work);
                near(step.energy.support_endpoint_work_j, endpoint_work);
                near(step.energy.transverse_endpoint_work_j, 0.0);
                near(step.energy.equilibrium_work_residual_j, 0.0);
                near(step.energy.balance_error_j, 0.0);
                assert!(step.energy.absolute_endpoint_work_defect_j > 0.0);
                assert!(
                    step.energy.absolute_endpoint_work_defect_j + 1e-12
                        >= step.energy.endpoint_work_defect_j.abs()
                );
                assert_eq!(step.point_energy.len(), original.len());
                near(step.sheet.physical_time_s(), DT);
            }
            assert_eq!(initial.quadrature_states(), original);
            assert_eq!(initial.displacement_m(), &[[0.0; 2]; 4]);
            assert_eq!(initial.physical_time_s(), 0.0);
        }
    });
}

#[test]
fn dead_nodal_forces_reach_the_independently_calculated_damaged_equilibrium() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let bands = bands(&mesh, cx);
        // This synthetic branch still hardens after damage starts. A target
        // force on a descending branch would also admit an intact solution.
        let work = 20.0;
        let material = model(bands[0].band_width_m(), work, 1.0);
        let materials = [&material];
        let thickness = [THICKNESS];
        let initial = FiniteDuctileSheet::try_new(
            &mesh,
            &materials,
            &thickness,
            &bands,
            local(),
            budget(),
            cx,
        )
        .unwrap();
        let e = 0.02;
        let alpha = equibiaxial_alpha(e, 0.0, 0.0, DT);
        let g = integrity(alpha, bands[0].band_width_m(), work);
        let piola = g * (Y + H * alpha) / e.exp();
        let mut supports = Vec::new();
        let mut forces = Vec::new();
        for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
            for (component, coordinate) in [x, y].into_iter().enumerate() {
                if coordinate == 0.0 {
                    supports.push((2 * node + component, 0.0));
                } else {
                    forces.push((2 * node + component, 0.5 * THICKNESS * piola));
                }
            }
        }
        assert_eq!(forces.len(), 4);
        let step = initial
            .advance(
                &supports,
                &forces,
                conditions(DT),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        for (actual, expected) in step
            .sheet
            .displacement_m()
            .iter()
            .zip(displacement(&mesh, e, 0.0))
        {
            for i in 0..2 {
                near(actual[i], expected[i]);
            }
        }
        for state in step.sheet.quadrature_states() {
            near(state.native_state().alpha, alpha);
            near(state.integrity(), g);
            assert!(state.damage() > 0.01);
        }
        let evaluation = initial
            .evaluate_trial(step.sheet.displacement_m(), conditions(DT), false, cx)
            .unwrap();
        for &(dof, force) in &forces {
            near(evaluation.internal_force_n.as_flattened()[dof], force);
        }
        near(
            step.energy.dead_load_endpoint_work_j,
            2.0 * THICKNESS * piola * e.exp_m1(),
        );
        near(step.energy.support_endpoint_work_j, 0.0);
        near(step.energy.balance_error_j, 0.0);
        assert!(
            step.equilibrium
                .residual_history
                .last()
                .is_some_and(|r| *r <= 1e-10)
        );
        assert!(
            initial
                .quadrature_states()
                .iter()
                .all(|s| s.damage() == 0.0)
        );
    });
}

#[test]
fn complete_assembled_tangent_matches_frozen_history_force_differences_and_replays() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let bands = bands(&mesh, cx);
        let material = model(bands[0].band_width_m(), GP, 1.0);
        let materials = [&material];
        let thickness = [THICKNESS];
        let initial = FiniteDuctileSheet::try_new(
            &mesh,
            &materials,
            &thickness,
            &bands,
            local(),
            budget(),
            cx,
        )
        .unwrap();
        let sheet = initial
            .advance(
                &prescribed(&displacement(&mesh, 0.012, 0.0)),
                &[],
                conditions(DT),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .sheet;
        let old_states = sheet.quadrature_states().to_vec();
        let old_u = sheet.displacement_m().to_vec();
        let old_storage = sheet.stored_energy_j();
        assert!(old_states.iter().all(|state| state.damage() > 0.0));
        let trial: Vec<_> = mesh
            .nodes
            .iter()
            .map(|&[x, y]| {
                [
                    0.024_f64.exp_m1() * x + 0.004 * y + 0.001 * x * y,
                    0.017_f64.exp_m1() * y - 0.003 * x + 0.0005 * x * y,
                ]
            })
            .collect();
        let response = sheet
            .evaluate_trial(&trial, conditions(DT), true, cx)
            .unwrap();
        let replay = sheet
            .evaluate_trial(&trial, conditions(DT), true, cx)
            .unwrap();
        assert_eq!(response.internal_force_n, replay.internal_force_n);
        assert_eq!(response.tangent_n_m, replay.tangent_n_m);
        let tangent = response.tangent_n_m.unwrap();
        let dofs = 2 * mesh.nodes.len();
        assert_eq!((tangent.nrows(), tangent.ncols()), (dofs, dofs));
        let h = 1e-6;
        let mut error = 0.0_f64;
        let mut norm = 0.0_f64;
        for column in 0..dofs {
            let mut plus = trial.clone();
            let mut minus = trial.clone();
            plus[column / 2][column % 2] += h;
            minus[column / 2][column % 2] -= h;
            let plus = sheet
                .evaluate_trial(&plus, conditions(DT), false, cx)
                .unwrap();
            let minus = sheet
                .evaluate_trial(&minus, conditions(DT), false, cx)
                .unwrap();
            for row in 0..dofs {
                let numerical = (plus.internal_force_n.as_flattened()[row]
                    - minus.internal_force_n.as_flattened()[row])
                    / (2.0 * h);
                let analytical = tangent.get(row, column);
                error = error.hypot(numerical - analytical);
                norm = norm.hypot(analytical);
            }
        }
        assert!(
            norm > 0.0 && error / norm < 2e-5,
            "relative full tangent error {}",
            error / norm
        );
        let translated: Vec<_> = trial.iter().map(|u| [u[0] + 2.0, u[1] - 3.0]).collect();
        let shifted = sheet
            .evaluate_trial(&translated, conditions(DT), false, cx)
            .unwrap();
        for (&a, &b) in response
            .internal_force_n
            .as_flattened()
            .iter()
            .zip(shifted.internal_force_n.as_flattened())
        {
            near(a, b);
        }
        assert_eq!(sheet.quadrature_states(), old_states);
        assert_eq!(sheet.displacement_m(), old_u);
        assert_eq!(sheet.stored_energy_j(), old_storage);
        assert_eq!(sheet.physical_time_s(), DT);
    });
}

#[test]
fn rate_independent_holds_and_elastic_unloading_retain_damage_and_permanent_shape() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let bands = bands(&mesh, cx);
        let material = model(bands[0].band_width_m(), GP, 1.0);
        let materials = [&material];
        let thickness = [THICKNESS];
        let initial = FiniteDuctileSheet::try_new(
            &mesh,
            &materials,
            &thickness,
            &bands,
            local(),
            budget(),
            cx,
        )
        .unwrap();
        let supports = prescribed(&displacement(&mesh, 0.02, 0.0));
        let loaded = initial
            .advance(
                &supports,
                &[],
                conditions(DT),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .sheet;
        let retained = loaded.quadrature_states().to_vec();
        let held = loaded
            .advance(
                &supports,
                &[],
                conditions(0.4),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(held.sheet.quadrature_states(), retained);
        near(held.energy.plastic_dissipation_j, 0.0);
        near(held.energy.damage_release_j, 0.0);
        near(held.energy.internal_endpoint_work_j, 0.0);
        near(held.sheet.physical_time_s(), DT + 0.4);
        let alpha = equibiaxial_alpha(0.02, 0.0, 0.0, DT);
        let released = loaded
            .advance(
                &prescribed(&displacement(&mesh, 0.5 * alpha, 0.0)),
                &[],
                conditions(DT),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        for (old, state) in retained.iter().zip(released.sheet.quadrature_states()) {
            assert_eq!(state.native_state(), old.native_state());
            assert_eq!(state.damage(), old.damage());
            near(state.deformation()[8], (-alpha).exp());
            assert!(state.first_piola_pa().iter().all(|p| p.abs() < 1e-6));
        }
        let g = integrity(alpha, bands[0].band_width_m(), GP);
        near(
            released.sheet.stored_energy_j(),
            THICKNESS * g * 0.5 * H * alpha * alpha,
        );
        near(released.energy.plastic_dissipation_j, 0.0);
        near(released.energy.damage_release_j, 0.0);
        near(released.energy.balance_error_j, 0.0);
        assert_eq!(loaded.quadrature_states(), retained);
    });
}

#[test]
fn physical_duration_and_held_overstress_advance_plasticity_and_damage_exactly_once() {
    with_cx(|_, cx| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let bands = bands(&mesh, cx);
        let viscosity = 400.0;
        let material = model(bands[0].band_width_m(), GP, 1.0)
            .with_overstress(FiniteOverstressParameters {
                viscosity_pa_s: viscosity,
                stress_scale_pa: Y,
                exponent: 1.0,
            })
            .unwrap();
        let materials = [&material];
        let thickness = [THICKNESS];
        let initial = FiniteDuctileSheet::try_new(
            &mesh,
            &materials,
            &thickness,
            &bands,
            local(),
            budget(),
            cx,
        )
        .unwrap();
        let e = 0.02;
        let supports = prescribed(&displacement(&mesh, e, 0.0));
        let fast = initial
            .advance(
                &supports,
                &[],
                conditions(0.02),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let slow = initial
            .advance(
                &supports,
                &[],
                conditions(2.0),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        for (dt, step) in [(0.02, &fast), (2.0, &slow)] {
            let alpha = equibiaxial_alpha(e, 0.0, viscosity, dt);
            let q = Y + H * alpha + viscosity * alpha / dt;
            let g = integrity(alpha, bands[0].band_width_m(), GP);
            for state in step.sheet.quadrature_states() {
                near(state.native_state().alpha, alpha);
                near(state.integrity(), g);
                near(state.first_piola_pa()[0], g * q / e.exp());
                near(
                    state.deformation()[8],
                    (2.0 * q / (3.0 * K) - 2.0 * e).exp(),
                );
            }
            near(
                step.energy.plastic_dissipation_j,
                THICKNESS * g * (Y + viscosity * alpha / dt) * alpha,
            );
            near(step.energy.balance_error_j, 0.0);
            near(step.sheet.physical_time_s(), dt);
        }
        assert!(
            fast.sheet.quadrature_states()[0].damage() < slow.sheet.quadrature_states()[0].damage()
        );
        assert!(
            fast.sheet.quadrature_states()[0].first_piola_pa()[0]
                > slow.sheet.quadrature_states()[0].first_piola_pa()[0]
        );
        let retained = fast.sheet.quadrature_states().to_vec();
        let old_alpha = equibiaxial_alpha(e, 0.0, viscosity, 0.02);
        let alpha = equibiaxial_alpha(e, old_alpha, viscosity, DT);
        let q = Y + H * alpha + viscosity * (alpha - old_alpha) / DT;
        let g = integrity(alpha, bands[0].band_width_m(), GP);
        let held = fast
            .sheet
            .advance(
                &supports,
                &[],
                conditions(DT),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        for (old, state) in retained.iter().zip(held.sheet.quadrature_states()) {
            near(state.native_state().alpha, alpha);
            near(state.integrity(), g);
            near(state.first_piola_pa()[0], g * q / e.exp());
            assert!(state.damage() > old.damage());
            assert!(state.first_piola_pa()[0] < old.first_piola_pa()[0]);
        }
        assert!(held.energy.plastic_dissipation_j > 0.0);
        assert!(held.energy.damage_release_j > 0.0);
        near(held.energy.internal_endpoint_work_j, 0.0);
        near(held.energy.balance_error_j, 0.0);
        near(held.sheet.physical_time_s(), 0.02 + DT);
        assert_eq!(fast.sheet.quadrature_states(), retained);
        assert_eq!(fast.sheet.physical_time_s(), 0.02);
        assert!(
            initial
                .quadrature_states()
                .iter()
                .all(|s| s.native_state().alpha == 0.0)
        );
    });
}

#[test]
fn directional_band_geometry_scales_and_failed_endpoints_preserve_accepted_state() {
    with_cx(|gate, cx| {
        let mesh = Mesh2::quads(2.0, 1.0, 1, 1);
        let x_band = SheetCrackBandGeometry::from_element(&mesh, 0, [1.0, 0.0], cx).unwrap();
        let y_band = SheetCrackBandGeometry::from_element(&mesh, 0, [0.0, 1.0], cx).unwrap();
        let diagonal = SheetCrackBandGeometry::from_element(&mesh, 0, [1.0, 1.0], cx).unwrap();
        near(x_band.band_width_m(), 2.0);
        near(y_band.band_width_m(), 1.0);
        near(diagonal.band_width_m(), 2.0 * 2.0_f64.sqrt() / 3.0);
        for scale in [1.0, 2.0] {
            let mut scaled = Mesh2::quads(2.0 * scale, scale, 1, 1);
            for x in &mut scaled.nodes {
                x[0] += 3.0;
                x[1] -= 5.0;
            }
            let bands = bands(&scaled, cx);
            near(bands[0].band_width_m(), 2.0 * scale);
            let material = model(bands[0].band_width_m(), GP, 1.0);
            let materials = [&material];
            let thickness = [THICKNESS];
            let sheet = FiniteDuctileSheet::try_new(
                &scaled,
                &materials,
                &thickness,
                &bands,
                local(),
                budget(),
                cx,
            )
            .unwrap();
            let step = sheet
                .advance(
                    &prescribed(&displacement(&scaled, 0.012, 0.0)),
                    &[],
                    conditions(DT),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let alpha = equibiaxial_alpha(0.012, 0.0, 0.0, DT);
            for state in step.sheet.quadrature_states() {
                near(state.integrity(), integrity(alpha, 2.0 * scale, GP));
            }
            near(
                step.sheet.reference_volume_m3(0).unwrap(),
                2.0 * scale * scale * THICKNESS,
            );
        }
        assert!(SheetCrackBandGeometry::from_element(&mesh, 0, [0.0, 0.0], cx).is_err());
        let bands = vec![x_band];
        let thickness = [THICKNESS];
        let wrong = model(0.5 * bands[0].band_width_m(), GP, 1.0);
        assert!(
            FiniteDuctileSheet::try_new(
                &mesh,
                &[&wrong],
                &thickness,
                &bands,
                local(),
                budget(),
                cx,
            )
            .is_err()
        );
        let material = model(bands[0].band_width_m(), GP, 1.0);
        let materials = [&material];
        let sheet = FiniteDuctileSheet::try_new(
            &mesh,
            &materials,
            &thickness,
            &bands,
            local(),
            budget(),
            cx,
        )
        .unwrap()
        .advance(
            &prescribed(&displacement(&mesh, 0.008, 0.0)),
            &[],
            conditions(DT),
            IncrementSettings::default(),
            controls(),
            cx,
        )
        .unwrap()
        .sheet;
        let old_states = sheet.quadrature_states().to_vec();
        let old_u = sheet.displacement_m().to_vec();
        let old_storage = sheet.stored_energy_j();
        let supports = prescribed(&displacement(&mesh, 0.02, 0.0));
        let proposal = sheet
            .advance(
                &supports,
                &[],
                conditions(DT),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert!(proposal.energy.absolute_endpoint_work_defect_j > 0.0);
        let strict = FinitePlasticStepControls {
            maximum_endpoint_work_defect_j: 0.5 * proposal.energy.absolute_endpoint_work_defect_j,
            ..controls()
        };
        assert!(matches!(
            sheet.advance(
                &supports,
                &[],
                conditions(DT),
                IncrementSettings::default(),
                strict,
                cx,
            ),
            Err(FiniteDuctileSheetError::Energy { .. })
        ));
        assert!(
            sheet
                .advance(
                    &supports,
                    &[],
                    conditions(0.0),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        let inverted: Vec<_> = mesh.nodes.iter().map(|&[x, _]| [-2.0 * x, 0.0]).collect();
        assert!(
            sheet
                .evaluate_trial(&inverted, conditions(DT), true, cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            sheet.advance(
                &supports,
                &[],
                conditions(DT),
                IncrementSettings::default(),
                controls(),
                cx,
            ),
            Err(FiniteDuctileSheetError::Solid(
                SmallStrainTetError::Cancelled
            ))
        ));
        assert_eq!(sheet.quadrature_states(), old_states);
        assert_eq!(sheet.displacement_m(), old_u);
        assert_eq!(sheet.stored_energy_j(), old_storage);
        assert_eq!(sheet.physical_time_s(), DT);
    });
}
