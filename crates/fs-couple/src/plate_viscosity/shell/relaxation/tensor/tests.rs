use super::super::tests::{phs, rk4};
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::sheet::SheetMatrix;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::visco_tensor::source::{
    TensorBranchProperties, TensorElasticity, TensorKelvinViscosityProperties,
    TensorPronyProperties,
};
use fs_math::c64::C64;
use fs_modal::ModePair;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_qty::{Density, Dims, Pressure};

mod kelvin;

fn source(band: (f64, f64), density: f64) -> Arc<ResolvedTensorMaterial> {
    source_with_kelvin(band, density, None)
}

fn source_with_kelvin(
    band: (f64, f64),
    density: f64,
    viscosity: Option<[[f64; 6]; 6]>,
) -> Arc<ResolvedTensorMaterial> {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims: Dims, value: f64| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290.0, 300.0)
                    .with("omega", band.0, band.1),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic full-tensor shear-shell law".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add("young_modulus", Pressure::DIMS, 2e6);
    add("poisson_ratio", Dims::NONE, 0.25);
    add("density", Density::DIMS, density);
    add("strain_limit", Dims::NONE, 0.1);
    add("tau", Dims([0, 0, 1, 0, 0, 0]), 0.02);
    // Exact integer coefficients make a strictly passive full source tensor.
    // In-plane and both transverse shears couple to the constrained zz history.
    let vector = [120.0, 70.0, -90.0, 45.0, 85.0, -65.0];
    let mut names = Vec::new();
    for i in 0..6 {
        for j in i..6 {
            let key = format!("branch-{i}-{j}");
            add(
                &key,
                Pressure::DIMS,
                vector[i] * vector[j] + if i == j { 5000.0 } else { 0.0 },
            );
            names.push(key);
        }
    }
    let kelvin = viscosity.map(|matrix| {
        let mut names = Vec::new();
        for i in 0..6 {
            for j in i..6 {
                let key = format!("viscosity-{i}-{j}");
                add(&key, Dims([-1, 1, -1, 0, 0, 0]), matrix[i][j]);
                names.push(key);
            }
        }
        TensorKelvinViscosityProperties {
            viscosity_upper: names.try_into().unwrap(),
        }
    });
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "MITC source-memory test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let resolved = ResolvedTensorMaterial::resolve(
        &card,
        &QueryPoint::new()
            .with("T", 295.0)
            .unwrap()
            .with("omega", band.0)
            .unwrap(),
        TensorPronyProperties {
            elasticity: TensorElasticity::Isotropic,
            strain_limit: "strain_limit".into(),
            branches: vec![TensorBranchProperties {
                stiffness_upper: names.try_into().unwrap(),
                relaxation_time: "tau".into(),
            }],
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    Arc::new(match kelvin {
        Some(properties) => resolved
            .with_kelvin_viscosity(
                &card,
                properties,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        None => resolved,
    })
}

fn facet(source: Arc<ResolvedTensorMaterial>, kappa: f64) -> TensorShellMaterial {
    TensorShellMaterial::new(
        source,
        [-0.005, 0.015],
        fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37),
        kappa,
    )
    .unwrap()
}

fn fixture(mass: ShellMass) -> (MitcShellModel, Vec<ModePair>, TensorShellMaterial) {
    fixture_with_source(mass, source((0.0, 1e8), 1000.0))
}

fn fixture_with_source(
    mass: ShellMass,
    source: Arc<ResolvedTensorMaterial>,
) -> (MitcShellModel, Vec<ModePair>, TensorShellMaterial) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    let material = facet(source, 0.8);
    let [h, i1, i2] = material.moments().unwrap();
    let c = 2e6 / (1.0 - 0.25 * 0.25);
    let sheet = |s| SheetMatrix::orthotropic(c * s, c * s, 0.25 * c * s, 8e5 * s).unwrap();
    let a = sheet(h);
    let section = ShearShellSection::new_coupled(
        a,
        [
            c * i1,
            0.25 * c * i1,
            0.0,
            0.25 * c * i1,
            c * i1,
            0.0,
            0.0,
            0.0,
            8e5 * i1,
        ],
        sheet(i2),
        [8e5 * h, 0.0, 8e5 * h],
        0.8,
        1000.0 * h,
        1000.0 * i1,
        1000.0 * i2,
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
        mass,
        MitcBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        1e-3,
        || false,
    )
    .unwrap();
    let n = model.pencil().free;
    let dense = |matrix: &fs_sparse::Csr| {
        (0..n)
            .flat_map(|i| (0..n).map(move |j| matrix.get(i, j)))
            .collect::<Vec<_>>()
    };
    let pairs =
        fs_modal::eigh_gen_dense(&dense(&model.pencil().k), &dense(&model.pencil().m), n).unwrap();
    (model, pairs, material)
}

fn limits() -> TensorShellRelaxationLimits {
    TensorShellRelaxationLimits {
        max_modes: 8,
        max_material_points: 32,
        max_memory_states: 18,
        max_projection_terms: 100_000_000,
    }
}

fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(1e-24),
        "{a:e} != {b:e}, scale={scale:e}"
    );
}

// Independent complex 3D Schur complement eliminates only sigma_zz. No use
// of the constrained point's energy matrix, factors or transformed rates.
fn source_modulus(law: &TensorProny, omega: f64) -> [[C64; 5]; 5] {
    let (storage, loss) = law.complex_modulus(omega).unwrap();
    let c: [[C64; 6]; 6] =
        core::array::from_fn(|i| core::array::from_fn(|j| C64::new(storage[i][j], loss[i][j])));
    let axes = [0, 1, 3, 4, 5];
    let weights = [
        1.0,
        1.0,
        core::f64::consts::FRAC_1_SQRT_2,
        core::f64::consts::FRAC_1_SQRT_2,
        core::f64::consts::FRAC_1_SQRT_2,
    ];
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (c[axes[i]][axes[j]] - c[axes[i]][2] * c[2][axes[j]] / c[2][2])
                .scale(weights[i] * weights[j])
        })
    })
}

#[test]
fn full_source_shell_keeps_shear_cross_work_and_complete_normal_history() {
    for mass in [ShellMass::Consistent, ShellMass::Lumped] {
        let (model, pairs, material) = fixture(mass);
        let modes: Vec<_> = pairs.iter().map(|m| m.phi.clone()).collect();
        let gate = CancelGate::new_clock_free();
        let image = TensorShellRelaxation::from_mitc_shell(
            &model,
            &[material.clone()],
            &modes,
            limits(),
            &gate,
        )
        .unwrap();
        assert_eq!(image.material_points().len(), 14);
        assert_eq!(
            image
                .material_points()
                .iter()
                .map(|p| p.memory_range().len())
                .sum::<usize>(),
            84
        );
        assert!(
            image.coupling_rows_per_s().len() <= 18,
            "equal-time original point histories must share exact driven coordinates"
        );
        assert!(std::ptr::eq(image.facets()[0].source(), material.source()));
        assert_eq!(image.band_rad_s(), (0.0, 1e8));
        assert!(
            image
                .relaxation_times_s()
                .iter()
                .any(|time| (*time - 0.02).abs() > 1e-8)
        );
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
        for omega in [0.0, 1.0, 50.0, 5000.0] {
            let c = source_modulus(&material.law, omega);
            let n = modes.len();
            let mut expected = vec![C64::ZERO; n * n];
            let h = 0.02;
            for (lambda, weight) in integration_rule() {
                for z in [
                    0.005 - h / 2.0 / 3.0_f64.sqrt(),
                    0.005 + h / 2.0 / 3.0_f64.sqrt(),
                ] {
                    let columns: Vec<[f64; 5]> = fields
                        .iter()
                        .map(|field| {
                            let s = field.at_barycentric(lambda).unwrap();
                            [
                                s.membrane[0] + z * s.gradient_per_m[0],
                                s.membrane[1] + z * s.gradient_per_m[1],
                                s.membrane[2] + z * s.gradient_per_m[2],
                                0.8_f64.sqrt() * s.transverse_shear[1],
                                0.8_f64.sqrt() * s.transverse_shear[0],
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
            for i in 0..n {
                for j in 0..n {
                    let equilibrium: f64 = image
                        .points
                        .iter()
                        .map(|point| point.image.equilibrium_stiffness()[n * i + j])
                        .sum();
                    let actual = image.modal.rows.iter().zip(&image.modal.times).fold(
                        C64::from_re(equilibrium),
                        |sum, (row, time)| {
                            sum + (C64::new(0.0, omega) / C64::new(time.recip(), omega))
                                .scale(row[i] * row[j])
                        },
                    );
                    let scale = (expected[n * i + i].re * expected[n * j + j].re).sqrt();
                    near(actual.re, expected[n * i + j].re, scale, 2e-9);
                    near(actual.im, expected[n * i + j].im, scale, 2e-9);
                }
            }
        }
        let q = [1e-5, -2e-5, 3e-5];
        let z = vec![0.0; image.coupling_rows_per_s().len()];
        let strain = image.point_strain(0, &q, &z).unwrap();
        assert!(strain[4].abs() + strain[5].abs() > 1e-9);
        let traction = image.points[0]
            .image
            .normal_stress_pa(&q, &image.point_memory(0, &z).unwrap())
            .unwrap();
        assert!(traction[1].abs() + traction[2].abs() > 1e-3);
        near(
            traction[0],
            0.0,
            traction[1].abs().max(traction[2].abs()),
            1e-10,
        );
        let memory: Vec<_> = z
            .iter()
            .enumerate()
            .map(|(i, _)| 1e-5 * (i as f64 - 2.0))
            .collect();
        let mut original_energy = 0.0;
        let mut original_power = 0.0;
        for (p, point) in image.points.iter().enumerate() {
            let local = image.point_memory(p, &memory).unwrap();
            for ((row, rate), z) in point
                .image
                .coupling_rows()
                .iter()
                .zip(point.image.relaxation_rates_per_s())
                .zip(local)
            {
                let residual = row.iter().zip(q).map(|(c, q)| c * q).sum::<f64>() - z;
                original_energy += 0.5 * residual * residual;
                original_power += rate * residual * residual;
            }
        }
        let mut shared_energy = 0.0;
        let mut shared_power = 0.0;
        for ((row, time), z) in image.modal.rows.iter().zip(&image.modal.times).zip(&memory) {
            let residual = row.iter().zip(q).map(|(c, q)| c * q).sum::<f64>() - z;
            shared_energy += 0.5 * residual * residual;
            shared_power += residual * residual / time;
        }
        near(shared_energy, original_energy, original_energy, 2e-10);
        near(shared_power, original_power, original_power, 2e-10);
    }
}

#[test]
fn source_shell_memory_drives_original_exact_runtime_and_independent_heat() {
    let (model, pairs, material) = fixture(ShellMass::Consistent);
    let modes: Vec<_> = pairs.iter().map(|m| m.phi.clone()).collect();
    let gate = CancelGate::new_clock_free();
    let image =
        TensorShellRelaxation::from_mitc_shell(&model, &[material], &modes, limits(), &gate)
            .unwrap();
    let omega2: Vec<_> = pairs.iter().map(|m| m.lambda).collect();
    let system = phs(&omega2, &image.modal, &gate);
    let states = system.state_dim();
    let dt = 1e-4;
    let prepared = PreparedExactZoh::new(
        &system,
        dt,
        LinearSplitBudget {
            max_modes: states,
            max_omega_dt: 0.25,
        },
    )
    .unwrap();
    let mut exact = vec![0.0; states];
    exact[0] = 1e-5;
    exact[3] = 0.003;
    let initial_energy = system.storage().hamiltonian(&exact);
    let mut reference = exact.clone();
    reference.extend([0.0, 0.0]);
    let force = [0.01, -0.02, 0.005];
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..24 {
        let step = prepared.step(&exact, &force).unwrap();
        assert!(step.dissipated >= 0.0);
        heat += step.dissipated;
        work += step.supplied;
        exact = step.x;
        for _ in 0..32 {
            rk4(&mut reference, &omega2, &image.modal, &force, dt / 32.0);
        }
    }
    assert!(heat > 0.0);
    near(
        system.storage().hamiltonian(&exact) + heat - work,
        initial_energy,
        initial_energy,
        2e-8,
    );
    near(heat, reference[states], heat, 2e-6);
    near(work, reference[states + 1], work.abs(), 2e-6);
    for i in 0..states {
        near(exact[i], reference[i], exact[i].abs().max(1e-8), 2e-6);
    }
    let q: Vec<_> = exact.iter().take(6).step_by(2).copied().collect();
    for point in 0..14 {
        image.point_strain(point, &q, &exact[6..]).unwrap();
    }
}

#[test]
fn source_shell_refuses_changed_mass_shear_band_basis_and_incomplete_history() {
    let (model, pairs, material) = fixture(ShellMass::Consistent);
    let modes: Vec<_> = pairs.iter().map(|m| m.phi.clone()).collect();
    let gate = CancelGate::new_clock_free();
    for incompatible in [
        facet(source((0.0, 1e8), 1100.0), 0.8),
        facet(source((0.0, 1e8), 1000.0), 0.7),
        facet(source((1.0, 1e8), 1000.0), 0.8),
        facet(source((0.0, 1.0), 1000.0), 0.8),
    ] {
        assert!(
            TensorShellRelaxation::from_mitc_shell(
                &model,
                &[incompatible],
                &modes,
                limits(),
                &gate
            )
            .is_err()
        );
    }
    for limit in [
        TensorShellRelaxationLimits {
            max_memory_states: 0,
            ..limits()
        },
        TensorShellRelaxationLimits {
            max_material_points: 13,
            ..limits()
        },
        TensorShellRelaxationLimits {
            max_projection_terms: 1,
            ..limits()
        },
    ] {
        assert!(
            TensorShellRelaxation::from_mitc_shell(
                &model,
                &[material.clone()],
                &modes,
                limit,
                &gate
            )
            .is_err()
        );
    }
    let mut bad = modes.clone();
    bad[0][0] *= 2.0;
    assert!(
        TensorShellRelaxation::from_mitc_shell(&model, &[material.clone()], &bad, limits(), &gate)
            .is_err()
    );
    let image = TensorShellRelaxation::from_mitc_shell(
        &model,
        &[material.clone()],
        &modes,
        limits(),
        &gate,
    )
    .unwrap();
    assert!(image.point_strain(0, &[0.0; 3], &[0.0; 84]).is_err());
    assert!(
        image
            .point_strain(0, &[1e6; 3], &vec![0.0; image.modal.rows.len()])
            .is_err()
    );
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(
        TensorShellRelaxation::from_mitc_shell(&model, &[material], &modes, limits(), &cancel)
            .is_err()
    );
}

#[test]
fn source_shell_retains_original_stiffness_and_refuses_unrelated_modal_mechanics() {
    use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
    let (model, pairs, material) = fixture(ShellMass::Consistent);
    let modes: Vec<_> = pairs.iter().map(|m| m.phi.clone()).collect();
    let gate = CancelGate::new_clock_free();
    let image = TensorShellRelaxation::from_mitc_shell(
        &model,
        &[material.clone()],
        &modes,
        limits(),
        &gate,
    )
    .unwrap();
    let component = |factor: f64, damping_ratio: f64| {
        ModalAcousticTimeModel::try_new(
            192_000,
            pairs
                .iter()
                .map(|p| ModalAcousticMode {
                    angular_frequency_rad_s: factor * p.lambda.sqrt(),
                    damping_ratio,
                    pressure_per_modal_velocity: C64::ZERO,
                })
                .collect(),
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap()
    };
    image.validate_component(&component(1.0, 0.0)).unwrap();
    assert!(image.validate_component(&component(1.1, 0.0)).is_err());
    assert!(image.validate_component(&component(1.0, 0.001)).is_err());
    for (&actual, expected) in image.equilibrium_omega_squared().iter().zip(&pairs) {
        near(actual, expected.lambda, expected.lambda, 1e-10);
    }
    let mut mixed = modes.clone();
    for i in 0..model.pencil().free {
        mixed[0][i] = (modes[0][i] + modes[1][i]) / 2.0_f64.sqrt();
        mixed[1][i] = (modes[0][i] - modes[1][i]) / 2.0_f64.sqrt();
    }
    // This change preserves the original mass metric but introduces elastic
    // cross work that independent modal oscillators cannot silently discard.
    validate_mass_basis(model.pencil(), &mixed, &gate).unwrap();
    assert!(
        TensorShellRelaxation::from_mitc_shell(&model, &[material], &mixed, limits(), &gate)
            .is_err()
    );
}
