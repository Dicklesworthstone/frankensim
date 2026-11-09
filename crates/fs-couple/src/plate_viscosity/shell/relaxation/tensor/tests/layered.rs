//! Layer-resolved tensor damping on the unchanged native MITC mechanics.
use super::*;
use crate::render::schedule::force::coupled::ModalCouplingConfig;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};

fn layer_limits() -> TensorShellRelaxationLimits {
    TensorShellRelaxationLimits {
        max_material_points: 28,
        max_memory_states: 64,
        ..limits()
    }
}

fn layered_fixture() -> (MitcShellModel, Vec<ModePair>, Vec<TensorShellMaterial>) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    let mut layers = Vec::new();
    let mut moments = [0.0; 3];
    let mut membrane = [0.0; 9];
    let mut coupling = [0.0; 9];
    let mut bending = [0.0; 9];
    let mut transverse = 0.0;
    for (faces, young, density, tau, kappa, angle, loss_scale) in [
        ([-0.006, 0.003], 2e6, 1000.0, 0.02, 0.8, 0.37, 1.0),
        ([0.003, 0.015], 5e6, 1700.0, 0.13, 0.6, -0.52, 2.3),
    ] {
        let v = [18.0, -12.0, 20.0, 9.0, -14.0, 11.0];
        let viscosity = core::array::from_fn(|i| {
            core::array::from_fn(|j| loss_scale * (v[i] * v[j] + if i == j { 7000.0 } else { 0.0 }))
        });
        let original = source_parameters((0.0, 1e8), density, Some(viscosity), young, tau);
        layers.push(
            TensorShellMaterial::new(
                original,
                faces,
                fs_material::tensor::rotation([0.0, 0.0, 1.0], angle),
                kappa,
            )
            .unwrap(),
        );
        // Independently integrated isotropic equilibrium, before the arbitrary
        // full nonproportional Maxwell and Kelvin increments are introduced.
        let [lo, hi]: [f64; 2] = faces;
        let integrals = [
            hi - lo,
            (hi * hi - lo * lo) / 2.0,
            (hi.powi(3) - lo.powi(3)) / 3.0,
        ];
        let e = young / (1.0 - 0.25 * 0.25);
        let g = young / 2.5;
        let c = [e, 0.25 * e, 0.0, 0.25 * e, e, 0.0, 0.0, 0.0, g];
        for i in 0..9 {
            membrane[i] += integrals[0] * c[i];
            coupling[i] += integrals[1] * c[i];
            bending[i] += integrals[2] * c[i];
        }
        transverse += integrals[0] * kappa * g;
        for i in 0..3 {
            moments[i] += density * integrals[i];
        }
    }
    assert!(coupling[0] > 0.0 && moments[1] > 0.0);
    let section = ShearShellSection::new_coupled(
        SheetMatrix::from_coefficients(membrane).unwrap(),
        coupling,
        SheetMatrix::from_coefficients(bending).unwrap(),
        [transverse, 0.0, transverse],
        1.0,
        moments[0],
        moments[1],
        moments[2],
        0.0,
    )
    .unwrap();
    let mut supports: Vec<_> = (0..2)
        .map(|node| ShellConstraint {
            node,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        })
        .collect();
    supports.push(ShellConstraint {
        node: 2,
        mask: [false, true, false, true, false, true],
        frame: ConstraintFrame::Local([[0.8, 0.0, 0.6], [0.0, 1.0, 0.0], [-0.6, 0.0, 0.8]]),
    });
    let model = assemble_mitc3p(
        &mesh,
        &directors,
        &[section],
        &supports,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        1e-3,
        || false,
    )
    .unwrap();
    let n = model.pencil().free;
    let dense = |a: &fs_sparse::Csr| {
        (0..n)
            .flat_map(|i| (0..n).map(move |j| a.get(i, j)))
            .collect::<Vec<_>>()
    };
    let pairs =
        fs_modal::eigh_gen_dense(&dense(&model.pencil().k), &dense(&model.pencil().m), n).unwrap();
    (model, pairs, layers)
}

fn frequency_matrix(image: &TensorShellRelaxation, omega: f64) -> Vec<C64> {
    let n = image.mode_count();
    (0..n * n)
        .map(|index| {
            let i = index / n;
            let j = index % n;
            let equilibrium = image
                .points
                .iter()
                .map(|point| point.image.equilibrium_stiffness()[index])
                .sum();
            image.modal.rows.iter().zip(&image.modal.times).fold(
                C64::new(equilibrium, omega * image.viscosity.matrix_per_s()[index]),
                |sum, (row, tau)| {
                    sum + (C64::new(0.0, omega) / C64::new(tau.recip(), omega))
                        .scale(row[i] * row[j])
                },
            )
        })
        .collect()
}

#[test]
fn g1_layered_shell_matches_original_full_tensor_frequency_work() {
    let (model, pairs, layers) = layered_fixture();
    let modes: Vec<_> = pairs.iter().map(|p| p.phi.clone()).collect();
    let image = TensorShellRelaxation::from_layered_mitc_shell(
        &model,
        &[layers.clone()],
        &modes,
        layer_limits(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(image.points.len(), 28);
    assert_eq!(image.facets().len(), 2);
    assert_eq!(image.facet_materials(0).unwrap().len(), 2);
    assert!(image.facet_materials(1).is_none());
    assert_eq!(
        image.points.iter().map(|p| p.memory.len()).sum::<usize>(),
        196
    );
    assert!(image.points.iter().all(|p| p.facet() == 0));
    for (i, layer) in layers.iter().enumerate() {
        assert_eq!(
            image
                .points
                .iter()
                .filter(|p| p.material_index() == i)
                .count(),
            14
        );
        assert!(std::ptr::eq(image.facets()[i].source(), layer.source()));
    }
    let (mesh, directors, sections) = model.original_geometry();
    let observer =
        PreparedElementStrain::new(mesh, 0, directors, &sections[0], &mut || false).unwrap();
    let fields: Vec<_> = modes
        .iter()
        .map(|mode| {
            let full = model.expand(mode).unwrap();
            observer
                .real_field(&core::array::from_fn(|i| full[i]), &mut || false)
                .unwrap()
        })
        .collect();
    let n = modes.len();
    for omega in [0.0, 0.3, 7.0, 50.0, 5000.0, 1e6] {
        let mut expected = vec![C64::ZERO; n * n];
        for layer in &layers {
            // The original 3D complex law, constrained once as a whole, is an
            // independent frequency-domain reference for the memory image.
            let c = source_modulus(&layer.law, omega);
            let h = layer.faces_m[1] - layer.faces_m[0];
            let center = 0.5 * (layer.faces_m[0] + layer.faces_m[1]);
            for (lambda, weight) in integration_rule() {
                for z in [
                    center - h / (2.0 * 3.0_f64.sqrt()),
                    center + h / (2.0 * 3.0_f64.sqrt()),
                ] {
                    let columns: Vec<_> = fields
                        .iter()
                        .map(|field| {
                            let s = field.at_barycentric(lambda).unwrap();
                            [
                                s.membrane[0] + z * s.gradient_per_m[0],
                                s.membrane[1] + z * s.gradient_per_m[1],
                                s.membrane[2] + z * s.gradient_per_m[2],
                                layer.shear_correction.sqrt() * s.transverse_shear[1],
                                layer.shear_correction.sqrt() * s.transverse_shear[0],
                            ]
                        })
                        .collect();
                    let volume = mesh.facet(0).unwrap().area_m2 * weight * h * 0.5;
                    for i in 0..n {
                        for j in 0..n {
                            for a in 0..5 {
                                for b in 0..5 {
                                    expected[n * i + j] = expected[n * i + j]
                                        + c[a][b].scale(volume * columns[i][a] * columns[j][b]);
                                }
                            }
                        }
                    }
                }
            }
        }
        let actual = frequency_matrix(&image, omega);
        for i in 0..n {
            for j in 0..n {
                let scale = expected[n * i + i]
                    .re
                    .abs()
                    .max(expected[n * i + i].im.abs())
                    .sqrt()
                    * expected[n * j + j]
                        .re
                        .abs()
                        .max(expected[n * j + j].im.abs())
                        .sqrt();
                near(actual[n * i + j].re, expected[n * i + j].re, scale, 3e-9);
                near(actual[n * i + j].im, expected[n * i + j].im, scale, 3e-9);
            }
        }
    }
}

#[test]
fn g1_layered_sources_advance_through_the_original_coupled_heat_owner() {
    let (model, pairs, layers) = layered_fixture();
    let modes: Vec<_> = pairs.iter().map(|p| p.phi.clone()).collect();
    let gate = CancelGate::new_clock_free();
    let image = TensorShellRelaxation::from_layered_mitc_shell(
        &model,
        &[layers],
        &modes,
        layer_limits(),
        &gate,
    )
    .unwrap();
    let mut component = ModalAcousticTimeModel::try_new(
        192_000,
        pairs
            .iter()
            .map(|p| ModalAcousticMode {
                angular_frequency_rad_s: p.lambda.sqrt(),
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    component
        .restore_states(
            &(0..modes.len())
                .map(|i| ModalAcousticState {
                    displacement_m_sqrt_kg: (i as f64 + 1.0) * 1e-6,
                    velocity_m_sqrt_kg_per_s: (i as f64 - 1.0) * 0.001,
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    let connections = image.kelvin_connections(0);
    let count = connections.len();
    let base = CoupledModalSystem::new(
        vec![component],
        connections,
        ModalCouplingConfig {
            max_modes: 8,
            max_connections: 8,
            max_setup_terms: 100_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e6,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-8,
        },
        &gate,
    )
    .unwrap();
    let mut system = image
        .attach_to_with_kelvin_connections(base, 0, 0..count, InitialPlateMemory::Relaxed, &gate)
        .unwrap();
    let initial = system.total_energy_j().unwrap();
    let mut heat = 0.0;
    for _ in 0..32 {
        let frame = system.step(&vec![0.0; modes.len()]).unwrap();
        assert_eq!(frame.component_dissipation_j, 0.0);
        assert!(frame.connection_dissipation_j >= 0.0 && frame.material_dissipation_j >= 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        heat += frame.connection_dissipation_j + frame.material_dissipation_j;
        let q: Vec<_> = system.components()[0]
            .states()
            .iter()
            .map(|s| s.displacement_m_sqrt_kg)
            .collect();
        for i in 0..image.points.len() {
            image.point_strain(i, &q, system.material_memory()).unwrap();
        }
    }
    assert!(heat > 0.0);
    near(
        system.total_energy_j().unwrap() + heat,
        initial,
        initial,
        2e-8,
    );
}

#[test]
fn g3_homogeneous_entry_preserves_bits_and_identical_layer_partition_preserves_work() {
    let (model, pairs, material) = fixture(ShellMass::Consistent);
    let modes: Vec<_> = pairs.iter().map(|p| p.phi.clone()).collect();
    let gate = CancelGate::new_clock_free();
    let old = TensorShellRelaxation::from_mitc_shell(
        &model,
        &[material.clone()],
        &modes,
        limits(),
        &gate,
    )
    .unwrap();
    let single = TensorShellRelaxation::from_layered_mitc_shell(
        &model,
        &[vec![material.clone()]],
        &modes,
        limits(),
        &gate,
    )
    .unwrap();
    assert_eq!(old.modal.rows, single.modal.rows);
    assert_eq!(old.modal.times, single.modal.times);
    assert_eq!(
        old.viscosity.matrix_per_s(),
        single.viscosity.matrix_per_s()
    );
    let mut low = material.clone();
    low.faces_m[1] = 0.0017;
    let mut high = material;
    high.faces_m[0] = 0.0017;
    let split = TensorShellRelaxation::from_layered_mitc_shell(
        &model,
        &[vec![low, high]],
        &modes,
        layer_limits(),
        &gate,
    )
    .unwrap();
    for omega in [0.0, 1.0, 50.0, 5000.0] {
        let a = frequency_matrix(&old, omega);
        let b = frequency_matrix(&split, omega);
        let n = modes.len();
        for i in 0..n {
            for j in 0..n {
                let scale = a[n * i + i].re.abs().sqrt() * a[n * j + j].re.abs().sqrt();
                near(a[n * i + j].re, b[n * i + j].re, scale, 2e-10);
                near(a[n * i + j].im, b[n * i + j].im, scale, 2e-10);
            }
        }
    }
}

#[test]
fn g0_layered_shell_refuses_inconsistent_geometry_sources_and_bounds() {
    let (model, pairs, layers) = layered_fixture();
    let modes: Vec<_> = pairs.iter().map(|p| p.phi.clone()).collect();
    let gate = CancelGate::new_clock_free();
    let make = |stack: Vec<TensorShellMaterial>, budget| {
        TensorShellRelaxation::from_layered_mitc_shell(&model, &[stack], &modes, budget, &gate)
    };
    let accepted = make(layers.clone(), layer_limits()).unwrap();
    for delta in [-0.0001, 0.0001] {
        let mut gap = layers.clone();
        gap[1].faces_m[0] += delta;
        assert!(make(gap, layer_limits()).is_err());
    }
    let mut reversed = layers.clone();
    reversed.reverse();
    assert!(make(reversed, layer_limits()).is_err());
    assert!(make(vec![], layer_limits()).is_err());
    assert!(make(vec![layers[0].clone()], layer_limits()).is_err());
    let mut stiffness = layers.clone();
    stiffness[1].shear_correction *= 0.9;
    assert!(make(stiffness, layer_limits()).is_err());
    let mut shifted = layers.clone();
    for layer in &mut shifted {
        layer.faces_m = layer.faces_m.map(|z| z + 0.001);
    }
    assert!(make(shifted, layer_limits()).is_err());
    for budget in [
        TensorShellRelaxationLimits {
            max_material_points: 27,
            ..layer_limits()
        },
        TensorShellRelaxationLimits {
            max_memory_states: 0,
            ..layer_limits()
        },
        TensorShellRelaxationLimits {
            max_projection_terms: 1,
            ..layer_limits()
        },
    ] {
        assert!(make(layers.clone(), budget).is_err());
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        TensorShellRelaxation::from_layered_mitc_shell(
            &model,
            &[layers.clone()],
            &modes,
            layer_limits(),
            &cancelled,
        )
        .is_err()
    );
    let retry = make(layers, layer_limits()).unwrap();
    assert_eq!(accepted.modal.rows, retry.modal.rows);
    assert_eq!(accepted.modal.times, retry.modal.times);
}
