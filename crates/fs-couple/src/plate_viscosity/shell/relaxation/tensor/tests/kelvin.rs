//! Combined source viscosity and normal-history loss on the original MITC basis.
use super::*;
use fs_phs::{PortHamiltonian, QuadraticStorage};

fn viscosity() -> [[f64; 6]; 6] {
    let v = [18.0, -12.0, 20.0, 9.0, -14.0, 11.0];
    core::array::from_fn(|i| {
        core::array::from_fn(|j| v[i] * v[j] + if i == j { 7000.0 } else { 0.0 })
    })
}

fn image_limits() -> TensorShellRelaxationLimits {
    TensorShellRelaxationLimits {
        max_memory_states: 24,
        ..limits()
    }
}

fn combined(
    mass: ShellMass,
) -> (
    MitcShellModel,
    Vec<ModePair>,
    TensorShellMaterial,
    TensorShellRelaxation,
) {
    let (model, pairs, material) = fixture_with_source(
        mass,
        source_with_kelvin((0.0, 1e8), 1000.0, Some(viscosity())),
    );
    let modes: Vec<_> = pairs.iter().map(|m| m.phi.clone()).collect();
    let image = TensorShellRelaxation::from_mitc_shell(
        &model,
        &[material.clone()],
        &modes,
        image_limits(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    (model, pairs, material, image)
}

#[test]
fn source_kelvin_shell_matches_original_3d_frequency_work_and_normal_history() {
    for mass in [ShellMass::Consistent, ShellMass::Lumped] {
        let (model, pairs, material, image) = combined(mass);
        assert!(material.source().kelvin_source().is_some());
        assert_eq!(
            image.points.iter().map(|p| p.memory.len()).sum::<usize>(),
            98
        );
        assert!(
            image
                .points
                .iter()
                .all(|p| p.image.point().normal_strain_coordinates() == 1)
        );
        assert!(image.viscosity.channel_count() > 0);
        let n = pairs.len();
        assert!(
            image
                .viscosity
                .matrix_per_s()
                .iter()
                .enumerate()
                .any(|(i, c)| i / n != i % n && c.abs() > 1e-4)
        );
        let (mesh, directors, sections) = model.original_geometry();
        let observer =
            PreparedElementStrain::new(mesh, 0, directors, &sections[0], &mut || false).unwrap();
        let fields: Vec<_> = pairs
            .iter()
            .map(|mode| {
                let full = model.expand(&mode.phi).unwrap();
                observer
                    .real_field(&core::array::from_fn(|i| full[i]), &mut || false)
                    .unwrap()
            })
            .collect();
        for omega in [0.0, 1.0, 50.0, 5000.0, 1e6] {
            // Independent original complex 6x6 law, with only sigma_zz
            // eliminated, integrated on the unmodified spatial trial fields.
            let c = source_modulus(&material.law, omega);
            let mut expected = vec![C64::ZERO; n * n];
            let h = material.faces_m[1] - material.faces_m[0];
            let center = material.faces_m[0] + 0.5 * h;
            for (lambda, weight) in integration_rule() {
                for z in [
                    center - h / 2.0 / 3.0_f64.sqrt(),
                    center + h / 2.0 / 3.0_f64.sqrt(),
                ] {
                    let columns: Vec<_> = fields
                        .iter()
                        .map(|field| {
                            let s = field.at_barycentric(lambda).unwrap();
                            [
                                s.membrane[0] + z * s.gradient_per_m[0],
                                s.membrane[1] + z * s.gradient_per_m[1],
                                s.membrane[2] + z * s.gradient_per_m[2],
                                material.shear_correction.sqrt() * s.transverse_shear[1],
                                material.shear_correction.sqrt() * s.transverse_shear[0],
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
                        .map(|p| p.image.equilibrium_stiffness()[n * i + j])
                        .sum();
                    let actual = image.modal.rows.iter().zip(&image.modal.times).fold(
                        C64::new(
                            equilibrium,
                            omega * image.viscosity.matrix_per_s()[n * i + j],
                        ),
                        |sum, (row, time)| {
                            sum + (C64::new(0.0, omega) / C64::new(time.recip(), omega))
                                .scale(row[i] * row[j])
                        },
                    );
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
                    near(actual.re, expected[n * i + j].re, scale, 3e-9);
                    near(actual.im, expected[n * i + j].im, scale, 3e-9);
                }
            }
        }
        let q = [1e-5, -2e-5, 3e-5];
        let v = [0.003, -0.005, 0.002];
        let history = vec![0.0; image.modal.rows.len()];
        for (i, point) in image.points.iter().enumerate() {
            let z = image.point_memory(i, &history).unwrap();
            assert!(point.image.original_memory(&z).is_err());
            assert!(point.image.normal_stress_pa(&q, &z).is_err());
            let stress = point
                .image
                .normal_stress_pa_with_velocity(&q, &v, &z)
                .unwrap();
            near(
                stress[0],
                0.0,
                stress[1].abs().max(stress[2].abs()).max(1.0),
                1e-9,
            );
            image.point_strain(i, &q, &history).unwrap();
        }
    }
}

#[test]
fn zero_source_kelvin_preserves_existing_maxwell_shell_bits() {
    let (model, pairs, plain) = fixture(ShellMass::Consistent);
    let zero = facet(
        source_with_kelvin((0.0, 1e8), 1000.0, Some([[0.0; 6]; 6])),
        0.8,
    );
    let modes: Vec<_> = pairs.iter().map(|m| m.phi.clone()).collect();
    let gate = CancelGate::new_clock_free();
    let make = |material| {
        TensorShellRelaxation::from_mitc_shell(&model, &[material], &modes, limits(), &gate)
            .unwrap()
    };
    let a = make(plain);
    let b = make(zero);
    assert!(a.kelvin_connections(0).is_empty());
    assert!(b.kelvin_connections(0).is_empty());
    assert_eq!(a.modal.rows, b.modal.rows);
    assert_eq!(
        a.modal
            .times
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>(),
        b.modal
            .times
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(a.equilibrium_omega_squared, b.equilibrium_omega_squared);
    for (a, b) in a.points.iter().zip(&b.points) {
        assert_eq!(a.image.coupling_rows(), b.image.coupling_rows());
        assert_eq!(
            a.image.relaxation_rates_per_s(),
            b.image.relaxation_rates_per_s()
        );
    }
}

#[test]
fn kelvin_shell_binding_checks_complete_physical_columns_and_coefficients() {
    let (_, _, _, image) = combined(ShellMass::Consistent);
    let connections = image.kelvin_connections(2);
    for (i, connection) in connections.iter().enumerate() {
        assert_eq!(connection.left.component, 2);
        assert_eq!(connection.right.component, 2);
        let mut column = vec![0.0; image.mode_count() + 3];
        for (out, (&a, &b)) in column[2..2 + image.mode_count()]
            .iter_mut()
            .zip(connection.left.shapes.iter().zip(&connection.right.shapes))
        {
            *out = a - b;
        }
        let validate = |c: &[f64], k, d, rest| image.validate_kelvin_channel(i, 2, c, k, d, rest);
        validate(&column, 0.0, connection.damping_n_s_m, 0.0).unwrap();
        assert!(validate(&column, 1.0, connection.damping_n_s_m, 0.0).is_err());
        assert!(validate(&column, 0.0, 0.5 * connection.damping_n_s_m, 0.0).is_err());
        assert!(validate(&column, 0.0, connection.damping_n_s_m, 1.0).is_err());
        column[0] = 1e-12;
        assert!(validate(&column, 0.0, connection.damping_n_s_m, 0.0).is_err());
        column[0] = 0.0;
        column[2] += 1e-12;
        assert!(validate(&column, 0.0, connection.damping_n_s_m, 0.0).is_err());
    }
}

#[test]
fn combined_shell_exact_motion_closes_separate_kelvin_and_memory_heat() {
    let (_, pairs, _, image) = combined(ShellMass::Consistent);
    let n = pairs.len();
    let d = 2 * n;
    let mut q = vec![0.0; d * d];
    let mut j = q.clone();
    let mut r = q.clone();
    let mut g = vec![0.0; d * n];
    for i in 0..n {
        q[d * (2 * i) + 2 * i] = pairs[i].lambda;
        q[d * (2 * i + 1) + 2 * i + 1] = 1.0;
        j[d * (2 * i) + 2 * i + 1] = 1.0;
        j[d * (2 * i + 1) + 2 * i] = -1.0;
        g[n * (2 * i + 1) + i] = 1.0;
        for k in 0..n {
            r[d * (2 * i + 1) + 2 * k + 1] = image.viscosity.matrix_per_s()[n * i + k];
        }
    }
    let gate = CancelGate::new_clock_free();
    let system = PortHamiltonian::new(
        d,
        n,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, d).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(image.branches(0, n, &gate).unwrap())
    .unwrap();
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
    let initial = system.storage().hamiltonian(&exact);
    let mut reference = exact.clone();
    reference.extend([0.0, 0.0, 0.0]); // independently integrated Kelvin heat, memory heat, input work
    let force = [0.01, -0.02, 0.005];
    let derivative = |x: &[f64]| {
        let mut out = vec![0.0; states + 3];
        for i in 0..n {
            out[2 * i] = x[2 * i + 1];
            out[2 * i + 1] = force[i] - pairs[i].lambda * x[2 * i];
            out[states + 2] += force[i] * x[2 * i + 1];
            for k in 0..n {
                let c = image.viscosity.matrix_per_s()[n * i + k];
                out[2 * i + 1] -= c * x[2 * k + 1];
                out[states] += x[2 * i + 1] * c * x[2 * k + 1];
            }
        }
        for (arm, (row, time)) in image.modal.rows.iter().zip(&image.modal.times).enumerate() {
            let extension = row
                .iter()
                .enumerate()
                .map(|(i, c)| c * x[2 * i])
                .sum::<f64>()
                - x[d + arm];
            out[d + arm] = extension / time;
            out[states + 1] += extension * extension / time;
            for i in 0..n {
                out[2 * i + 1] -= row[i] * extension;
            }
        }
        out
    };
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..24 {
        let step = prepared.step(&exact, &force).unwrap();
        assert!(step.dissipated >= 0.0);
        heat += step.dissipated;
        work += step.supplied;
        exact = step.x;
        for _ in 0..32 {
            let h = dt / 32.0;
            let stage = |rate: &[f64], scale: f64| {
                reference
                    .iter()
                    .zip(rate)
                    .map(|(x, k)| x + scale * k)
                    .collect::<Vec<_>>()
            };
            let k1 = derivative(&reference);
            let k2 = derivative(&stage(&k1, h / 2.0));
            let k3 = derivative(&stage(&k2, h / 2.0));
            let k4 = derivative(&stage(&k3, h));
            for i in 0..reference.len() {
                reference[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
        }
    }
    assert!(reference[states] > 0.0 && reference[states + 1] > 0.0);
    near(
        heat,
        reference[states] + reference[states + 1],
        initial,
        2e-8,
    );
    near(work, reference[states + 2], initial, 2e-8);
    near(
        system.storage().hamiltonian(&exact) - initial + heat - work,
        0.0,
        initial,
        1e-9,
    );
    for (&actual, &expected) in exact.iter().zip(&reference) {
        near(actual, expected, expected.abs().max(1e-8), 3e-7);
    }
}
