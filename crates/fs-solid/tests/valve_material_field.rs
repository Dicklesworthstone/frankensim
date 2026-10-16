//! Synthetic regression specimens, not measurements of cane or an instrument.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::{TetElasticMaterial, TetMaterialField};
use fs_solid::reduce::material_field;
use fs_solid::reduce::{OrificeSpec, ValveCardRequest};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 77,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

struct Mesh {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    fixed: Vec<usize>,
    tip: Vec<usize>,
}

impl Mesh {
    fn new() -> Self {
        let (nx, ny, nz) = (2usize, 2usize, 2usize);
        let node = |i, j, k| (i * (ny + 1) + j) * (nz + 1) + k;
        let mut nodes = Vec::new();
        for i in 0..=nx {
            for j in 0..=ny {
                for k in 0..=nz {
                    nodes.push([
                        0.010 * i as f64 / nx as f64,
                        0.012 * j as f64 / ny as f64,
                        0.004 * k as f64 / nz as f64,
                    ]);
                }
            }
        }
        let mut tets = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let v = [
                        node(i, j, k),
                        node(i + 1, j, k),
                        node(i, j + 1, k),
                        node(i + 1, j + 1, k),
                        node(i, j, k + 1),
                        node(i + 1, j, k + 1),
                        node(i, j + 1, k + 1),
                        node(i + 1, j + 1, k + 1),
                    ];
                    for t in [
                        [0, 1, 3, 7],
                        [0, 3, 2, 7],
                        [0, 2, 6, 7],
                        [0, 6, 4, 7],
                        [0, 4, 5, 7],
                        [0, 5, 1, 7],
                    ] {
                        tets.push(t.map(|index| v[index]));
                    }
                }
            }
        }
        let mut fixed = Vec::new();
        let mut tip = Vec::new();
        for (n, p) in nodes.iter().enumerate() {
            if p[0] == 0.0 {
                fixed.extend([3 * n, 3 * n + 1, 3 * n + 2]);
            }
            if p[2] == 0.0 && p[0] >= 0.005 {
                tip.push(n);
            }
        }
        Self {
            nodes,
            tets,
            fixed,
            tip,
        }
    }

    fn request<'a>(
        &'a self,
        material: TetMaterialField<'a>,
    ) -> ValveCardRequest<'a, TetMaterialField<'a>> {
        ValveCardRequest {
            element_order: fs_solid::reduce::ValveElementOrder::Tet4,
            nodes_m: &self.nodes,
            tetrahedra: &self.tets,
            fixed_dofs: &self.fixed,
            material,
            orifice: OrificeSpec {
                face_nodes: &self.tip,
                opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0],
                opposing_plane_offset_m: -0.0004,
            },
            window_hz: (1.0, 50_000.0),
            retained_compliance_floor: 0.5,
            loss_factor: 0.1,
            source_id: "synthetic/material-field-regression/v1",
        }
    }
}

fn material(e: f64, rho: f64, identity: u8) -> TetElasticMaterial {
    TetElasticMaterial::try_new(rho, e, 0.30, [identity; 32]).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn uniform_and_explicit_identical_field_have_identical_cards() {
    with_cx(|cx| {
        let mesh = Mesh::new();
        let mat = material(1e6, 1050.0, 1);
        let repeated = vec![mat.clone(); mesh.tets.len()];
        let a = material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&mat)), cx)
            .unwrap();
        let b = material_field::reduce_valve(
            &mesh.request(TetMaterialField::PerElement(&repeated)),
            cx,
        )
        .unwrap();
        assert_eq!(a, b);
    });
}

#[test]
fn density_and_stiffness_swaps_reach_the_same_modal_pencil() {
    with_cx(|cx| {
        let mesh = Mesh::new();
        let base = material(1e6, 1050.0, 1);
        let dense = material(1e6, 2100.0, 2);
        let stiff = material(4e6, 1050.0, 3);
        let a = material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&base)), cx)
            .unwrap();
        let b = material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&dense)), cx)
            .unwrap();
        let c = material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&stiff)), cx)
            .unwrap();
        close(b.total_mass_kg, 2.0 * a.total_mass_kg, 1e-12);
        close(c.total_mass_kg, a.total_mass_kg, 1e-12);
        close(
            b.modes[0].frequency_hz,
            a.modes[0].frequency_hz / 2.0f64.sqrt(),
            1e-7,
        );
        close(c.modes[0].frequency_hz, 2.0 * a.modes[0].frequency_hz, 1e-7);
        close(b.modes[0].mass_kg, 2.0 * a.modes[0].mass_kg, 1e-6);
        close(
            c.modes[0].stiffness_n_m,
            4.0 * a.modes[0].stiffness_n_m,
            1e-6,
        );
        assert_ne!(a.material_identity, b.material_identity);
        assert_ne!(a.material_identity, c.material_identity);
    });
}

#[test]
fn nonuniform_material_assignment_preserves_actual_mass_and_provenance() {
    with_cx(|cx| {
        let mesh = Mesh::new();
        let light = material(1e6, 1000.0, 1);
        let heavy = material(1e6, 3000.0, 2);
        let mut field: Vec<_> = (0..mesh.tets.len())
            .map(|i| {
                if i < mesh.tets.len() / 2 {
                    light.clone()
                } else {
                    heavy.clone()
                }
            })
            .collect();
        let a =
            material_field::reduce_valve(&mesh.request(TetMaterialField::PerElement(&field)), cx)
                .unwrap();
        close(a.total_mass_kg, 0.010 * 0.012 * 0.004 * 2000.0, 1e-12);
        field.reverse();
        let b =
            material_field::reduce_valve(&mesh.request(TetMaterialField::PerElement(&field)), cx)
                .unwrap();
        close(a.total_mass_kg, b.total_mass_kg, 1e-12);
        assert_ne!(a.material_identity, b.material_identity);
        assert!((a.modes[0].frequency_hz / b.modes[0].frequency_hz - 1.0).abs() > 0.01);
    });
}

#[test]
fn full_anisotropic_tensor_is_not_replaced_by_an_effective_young_modulus() {
    with_cx(|cx| {
        let mesh = Mesh::new();
        let mut c = [[0.0; 6]; 6];
        for (i, row) in c.iter_mut().enumerate() {
            row[i] = 1e6;
        }
        c[0][0] = 4e6;
        let along = TetElasticMaterial::try_new_mandel(1050.0, c, [4; 32]).unwrap();
        c[0][0] = 1e6;
        c[1][1] = 4e6;
        let across = TetElasticMaterial::try_new_mandel(1050.0, c, [5; 32]).unwrap();
        let a = material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&along)), cx)
            .unwrap();
        let b = material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&across)), cx)
            .unwrap();
        close(a.total_mass_kg, b.total_mass_kg, 1e-12);
        assert!((a.modes[0].frequency_hz / b.modes[0].frequency_hz - 1.0).abs() > 0.01);
        assert_ne!(a.material_identity, b.material_identity);
    });
}

#[test]
fn pressure_port_keeps_sign_and_material_dependent_modes() {
    with_cx(|cx| {
        let mesh = Mesh::new();
        let mat = material(1e6, 1050.0, 1);
        let mut load = vec![[0.0; 3]; mesh.nodes.len()];
        for &n in &mesh.tip {
            load[n][2] = 1e-4 / mesh.tip.len() as f64;
        }
        let request = mesh.request(TetMaterialField::Uniform(&mat));
        let a = material_field::reduce_pressure_loaded_valve(&request, &load, cx).unwrap();
        for value in load.iter_mut().flatten() {
            *value = -*value;
        }
        let b = material_field::reduce_pressure_loaded_valve(&request, &load, cx).unwrap();
        assert_eq!(a.card(), b.card());
        for (&x, &y) in a.effective_areas_m2().iter().zip(b.effective_areas_m2()) {
            close(x, -y, 1e-12);
        }
        assert_ne!(a.identity(), b.identity());
        assert!(a.beating_reed(100.0, 0.01).is_ok());
        assert!(b.beating_reed(100.0, 0.01).is_err());
    });
}

#[test]
fn incomplete_or_mutated_material_fields_refuse_in_the_common_owner() {
    with_cx(|cx| {
        let mesh = Mesh::new();
        let mut mat = material(1e6, 1050.0, 1);
        assert!(
            material_field::reduce_valve(&mesh.request(TetMaterialField::PerElement(&[])), cx)
                .is_err()
        );
        mat.density_kg_m3 = f64::NAN;
        assert!(
            material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&mat)), cx)
                .is_err()
        );
        mat.density_kg_m3 = 1050.0;
        mat.material_state_identity = [0; 32];
        assert!(
            material_field::reduce_valve(&mesh.request(TetMaterialField::Uniform(&mat)), cx)
                .is_err()
        );
    });
}

#[test]
fn pressure_compliance_matches_direct_solve_and_checks_the_playing_image() {
    use fs_solid::linear3::{TetAssemblyBudget, TetLinearElasticProblem};
    use fs_solid::reduce::SingleModeValveAdmission;
    with_cx(|cx| {
        let mesh = Mesh::new();
        let mat = material(1e6, 1050.0, 1);
        let request = mesh.request(TetMaterialField::Uniform(&mat));
        let mut load = vec![[0.0; 3]; mesh.nodes.len()];
        for &n in &mesh.tip {
            load[n][2] = 1e-4 / mesh.tip.len() as f64;
        }
        let valve = material_field::reduce_pressure_loaded_valve(&request, &load, cx).unwrap();
        let assembly = TetLinearElasticProblem {
            nodes_m: &mesh.nodes,
            tetrahedra: &mesh.tets,
            fixed_dofs: &mesh.fixed,
            materials: TetMaterialField::Uniform(&mat),
            budget: TetAssemblyBudget::default(),
        }
        .assemble(cx)
        .unwrap();
        let symbolic =
            fs_sparse::SymbolicLdlt::analyze(&assembly.stiffness, fs_sparse::DirectOrdering::Amd)
                .unwrap();
        let factor = symbolic
            .factor(&assembly.stiffness, &fs_sparse::LdltOptions::default())
            .unwrap();
        let reduced_load: Vec<f64> = assembly
            .free_dofs
            .iter()
            .map(|&d| load[d / 3][d % 3])
            .collect();
        let displacement = factor.solve(&reduced_load);
        let volume: f64 = reduced_load
            .iter()
            .zip(&displacement)
            .map(|(a, b)| a * b)
            .sum();
        let opening: f64 = assembly
            .free_dofs
            .iter()
            .zip(&displacement)
            .filter(|(d, _)| **d % 3 == 2 && mesh.tip.contains(&(**d / 3)))
            .map(|(_, x)| *x)
            .sum::<f64>()
            / mesh.tip.len() as f64;
        let c = valve.pressure_compliance();
        close(c.exact_volume_m3_per_pa(), volume, 1e-12);
        close(c.exact_opening_m_per_pa(), opening, 1e-12);
        // This broad test window contains the full small synthetic pencil.
        close(c.retained_volume_m3_per_pa(), volume, 1e-6);
        close(c.retained_opening_m_per_pa(), opening, 1e-6);
        let admission = SingleModeValveAdmission {
            minimum_opening_compliance_fraction: 0.90,
            minimum_pressure_compliance_fraction: 0.90,
            maximum_static_opening_relative_error: 0.10,
        };
        let reed = valve.beating_reed_checked(100.0, 0.01, admission).unwrap();
        close(
            reed.stiffness_n_m * reed.rest_opening_m / reed.closing_pressure_pa,
            valve.effective_areas_m2()[0],
            1e-12,
        );
        let strict = SingleModeValveAdmission {
            minimum_opening_compliance_fraction: 0.999999,
            ..admission
        };
        assert!(valve.beating_reed_checked(100.0, 0.01, strict).is_err());
        assert!(
            valve
                .beating_reed_checked(f64::NAN, 0.01, admission)
                .is_err()
        );
        for x in load.iter_mut().flatten() {
            *x = -*x;
        }
        let reversed = material_field::reduce_pressure_loaded_valve(&request, &load, cx).unwrap();
        close(
            reversed.pressure_compliance().exact_volume_m3_per_pa(),
            volume,
            1e-12,
        );
        close(
            reversed.pressure_compliance().exact_opening_m_per_pa(),
            -opening,
            1e-12,
        );
        assert!(
            reversed
                .beating_reed_checked(100.0, 0.01, admission)
                .is_err()
        );
    });
}
