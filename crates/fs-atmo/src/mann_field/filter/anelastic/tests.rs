use super::*;
use crate::mann::MANN_TARGET_V1;
use crate::mann_field::fixed_step::FixedStepMannField;
use crate::mann_field::{HomogeneousMannField, MannQuadrature, MannQuadratureRule};

fn field() -> HomogeneousMannField {
    let model = MannSpectralModel::build_with_rule(
        17,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 12,
            min_kl: 0.05,
            max_kl: 3.0,
        },
        [2.0, -0.3, 0.0],
        MannQuadratureRule::BalancedEnergyV1,
    )
    .unwrap();
    let mut field = HomogeneousMannField::at_anchor(model, 39, vec![0.8; 12]).unwrap();
    field.advance_to(field.tick() + 3).unwrap();
    field
}

fn density(kappa: f64) -> ExponentialDensityProfile {
    ExponentialDensityProfile::new(1.2, 10.0, kappa, [0.0, 20.0]).unwrap()
}

fn prepare(
    field: &HomogeneousMannField,
    density: ExponentialDensityProfile,
    edges: [[f64; 3]; 3],
) -> MannAnelasticBoxMoments {
    MannAnelasticBoxMoments::prepare(field.model(), density, edges, 1_000_000, 1000, || false)
        .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} vs {b:e}"
    );
}

// Independent Gauss-Legendre spatial integration of the public point sampler.
// This does not evaluate any of the prepared exponential-product coefficients.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    let mut rule = Vec::new();
    for i in 0..n {
        let mut x = (core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 0.0;
        for _ in 0..30 {
            let (mut a, mut b) = (1.0, x);
            for degree in 2..=n {
                let c = ((2 * degree - 1) as f64 * x * b - (degree - 1) as f64 * a) / degree as f64;
                a = b;
                b = c;
            }
            derivative = n as f64 * (x * b - a) / (x * x - 1.0);
            let change = b / derivative;
            x -= change;
            if change.abs() < 2e-16 {
                break;
            }
        }
        rule.push((x, 1.0 / ((1.0 - x * x) * derivative * derivative)));
    }
    rule
}

fn integrate(
    field: &HomogeneousMannField,
    density: ExponentialDensityProfile,
    center: [f64; 3],
    edges: [[f64; 3]; 3],
    n: usize,
) -> MannFavreMoments {
    let rule = gauss(n);
    let mut values = Vec::new();
    let mut rho = 0.0;
    let mut mass = [0.0; 3];
    let mut momentum = [[0.0; 3]; 3];
    for &(a, wa) in &rule {
        for &(b, wb) in &rule {
            for &(c, wc) in &rule {
                let x = core::array::from_fn(|j| {
                    center[j] + a * edges[0][j] + b * edges[1][j] + c * edges[2][j]
                });
                let value = field.sample_anelastic(x, &density).unwrap();
                let weight = wa * wb * wc;
                rho += weight * value.density_kg_m3;
                for i in 0..3 {
                    mass[i] += weight * value.mass_flux_kg_m2_s[i];
                    for j in 0..3 {
                        momentum[i][j] += weight * value.mass_flux_kg_m2_s[i] * value.flow.u[j];
                    }
                }
                values.push((weight, value));
            }
        }
    }
    let mean = mass.map(|m| m / rho);
    let mut stress = [[0.0; 3]; 3];
    for (weight, value) in values {
        for i in 0..3 {
            for j in 0..3 {
                stress[i][j] += weight * value.density_kg_m3 / rho
                    * (value.flow.u[i] - mean[i])
                    * (value.flow.u[j] - mean[j]);
            }
        }
    }
    MannFavreMoments {
        mean_density_kg_m3: rho,
        mean_mass_flux_kg_m2_s: mass,
        favre_mean: FlowSample {
            u: mean,
            grad: [[0.0; 3]; 3],
        },
        favre_subfilter_stress_m2_s2: stress,
        convective_momentum_flux_pa: momentum,
        kinetic_energy_density_j_m3: 0.5 * (0..3).map(|i| momentum[i][i]).sum::<f64>(),
    }
}

#[test]
fn oblique_spatial_integrals_recover_mass_momentum_stress_and_energy() {
    let field = field();
    let center = [0.3, -0.6, 9.6];
    let edges = [[0.5, 0.1, 0.1], [0.1, 0.3, -0.12], [0.08, -0.05, 0.6]];
    for kappa in [0.0, 0.4, 3.0] {
        let density = density(kappa);
        let observation = prepare(&field, density, edges);
        let got = field
            .sample_anelastic_filtered_moments(center, &observation, || false)
            .unwrap();
        let coarse = integrate(&field, density, center, edges, 12);
        let reference = integrate(&field, density, center, edges, 20);
        close(got.mean_density_kg_m3, reference.mean_density_kg_m3, 2e-12);
        close(
            got.kinetic_energy_density_j_m3,
            reference.kinetic_energy_density_j_m3,
            3e-11,
        );
        for i in 0..3 {
            close(
                got.mean_mass_flux_kg_m2_s[i],
                reference.mean_mass_flux_kg_m2_s[i],
                3e-12,
            );
            close(got.favre_mean.u[i], reference.favre_mean.u[i], 3e-12);
            for j in 0..3 {
                close(
                    coarse.convective_momentum_flux_pa[i][j],
                    reference.convective_momentum_flux_pa[i][j],
                    2e-10,
                );
                close(
                    got.convective_momentum_flux_pa[i][j],
                    reference.convective_momentum_flux_pa[i][j],
                    3e-11,
                );
                close(
                    got.favre_subfilter_stress_m2_s2[i][j],
                    reference.favre_subfilter_stress_m2_s2[i][j],
                    3e-11,
                );
            }
        }
        close(
            got.favre_mean.grad[0][0] + got.favre_mean.grad[1][1] + got.favre_mean.grad[2][2],
            kappa * got.favre_mean.u[2],
            2e-13,
        );
        for axis in 0..3 {
            let mut a = center;
            let mut b = center;
            a[axis] += 1e-4;
            b[axis] -= 1e-4;
            let plus = integrate(&field, density, a, edges, 12);
            let minus = integrate(&field, density, b, edges, 12);
            for i in 0..3 {
                close(
                    got.favre_mean.grad[i][axis],
                    (plus.favre_mean.u[i] - minus.favre_mean.u[i]) / 2e-4,
                    2e-7,
                );
            }
        }
    }
}

#[test]
fn tiny_footprints_retain_the_positive_gradient_stress() {
    let field = field();
    let density = density(0.8);
    let center = [0.3, -0.6, 10.0];
    let gradient = field.sample_anelastic(center, &density).unwrap().flow.grad;
    for size in [1e-4, 1e-8, 1e-12] {
        let edges = [[0.5, 0.1, 0.1], [0.1, 0.3, -0.12], [0.08, -0.05, 0.6]]
            .map(|row| row.map(|v| v * size));
        let got = field
            .sample_anelastic_filtered_moments(center, &prepare(&field, density, edges), || false)
            .unwrap();
        let expected: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                edges
                    .iter()
                    .map(|edge| dot(gradient[i], *edge) * dot(gradient[j], *edge) / 3.0)
                    .sum()
            })
        });
        let scale = (0..3).map(|i| expected[i][i]).sum::<f64>();
        assert!(scale > 0.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (got.favre_subfilter_stress_m2_s2[i][j] - expected[i][j]).abs() < 3e-7 * scale
                );
            }
        }
        for direction in [[1.0, 0.0, 0.0], [0.3, -0.7, 0.4], [0.2, 0.8, -1.1]] {
            let variance = (0..3)
                .map(|i| direction[i] * dot(got.favre_subfilter_stress_m2_s2[i], direction))
                .sum::<f64>();
            assert!(variance >= -1e-12 * scale);
        }
    }
}

#[test]
fn constant_solenoidal_source_has_the_exact_density_induced_variance() {
    // A manufactured zero-frequency transverse source isolates the density
    // operator; it is not asserted to be a Mann spectral realization.
    let mut model = field().model().clone();
    for mode in &mut model.modes {
        mode.k = [0.0; 3];
        mode.polarization = [[0.0; 3]; 2];
    }
    model.modes[0].polarization[0] = [1.0, 0.0, 0.0];
    let field = HomogeneousMannField::at_anchor(model, 39, vec![0.8; 12]).unwrap();
    let center = [0.0, 0.0, 10.0];
    let edges = [[0.0; 3], [0.0; 3], [0.0, 0.0, 1.0]];
    let v = field.sample(center).unwrap().u[0];
    for a in [1e-10, 1e-4, 0.7, 8.0, 30.0, 700.0] {
        let density = ExponentialDensityProfile::new(1.2, 10.0, a, [9.0, 11.0]).unwrap();
        let got = field
            .sample_anelastic_filtered_moments(center, &prepare(&field, density, edges), || false)
            .unwrap();
        let inverse_mean = if a < 1e-3 {
            1.0 - a * a / 6.0 + 7.0 * a.powi(4) / 360.0
        } else {
            2.0 * a * (-a).exp() / (-(-2.0 * a).exp_m1())
        };
        let variance = if a < 1e-3 {
            a * a / 3.0 - a.powi(4) / 15.0 + 2.0 * a.powi(6) / 189.0
        } else {
            1.0 - inverse_mean * inverse_mean
        };
        let actual = got.favre_subfilter_stress_m2_s2[0][0];
        assert!(
            actual >= 0.0 && (actual / (v * v * variance) - 1.0).abs() < 2e-11,
            "a={a}: {actual:e}"
        );
        assert!((got.favre_mean.u[0] / (v * inverse_mean) - 1.0).abs() < 2e-12);
        close(
            got.convective_momentum_flux_pa[0][0] / got.mean_density_kg_m3,
            v * v,
            2e-12,
        );
        assert_eq!(got.mean_mass_flux_kg_m2_s, [1.2 * v, 0.0, 0.0]);
    }
}

#[test]
fn finite_volume_partition_and_closed_surface_preserve_physical_fluxes() {
    let field = field();
    let center = [0.3, -0.6, 10.0];
    let edges = [[0.4, 0.1, 0.05], [-0.04, 0.3, 0.1], [0.1, -0.1, 0.5]];
    let density = density(1.5);
    let whole = field
        .sample_anelastic_filtered_moments(center, &prepare(&field, density, edges), || false)
        .unwrap();
    let mut half_edges = edges;
    half_edges[0] = edges[0].map(|v| 0.5 * v);
    let half = prepare(&field, density, half_edges);
    let children = [-1.0, 1.0].map(|sign| {
        field
            .sample_anelastic_filtered_moments(
                core::array::from_fn(|i| center[i] + sign * 0.5 * edges[0][i]),
                &half,
                || false,
            )
            .unwrap()
    });
    close(
        whole.mean_density_kg_m3,
        0.5 * (children[0].mean_density_kg_m3 + children[1].mean_density_kg_m3),
        2e-13,
    );
    close(
        whole.kinetic_energy_density_j_m3,
        0.5 * (children[0].kinetic_energy_density_j_m3 + children[1].kinetic_energy_density_j_m3),
        2e-12,
    );
    for i in 0..3 {
        close(
            whole.mean_mass_flux_kg_m2_s[i],
            0.5 * (children[0].mean_mass_flux_kg_m2_s[i] + children[1].mean_mass_flux_kg_m2_s[i]),
            2e-13,
        );
        for j in 0..3 {
            close(
                whole.convective_momentum_flux_pa[i][j],
                0.5 * (children[0].convective_momentum_flux_pa[i][j]
                    + children[1].convective_momentum_flux_pa[i][j]),
                3e-12,
            );
            let combined = children
                .iter()
                .map(|child| {
                    0.5 * child.mean_density_kg_m3 / whole.mean_density_kg_m3
                        * (child.favre_subfilter_stress_m2_s2[i][j]
                            + (child.favre_mean.u[i] - whole.favre_mean.u[i])
                                * (child.favre_mean.u[j] - whole.favre_mean.u[j]))
                })
                .sum::<f64>();
            close(whole.favre_subfilter_stress_m2_s2[i][j], combined, 3e-12);
        }
    }
    let mut balance = 0.0;
    let mut scale = 0.0;
    for axis in 0..3 {
        let a = edges[(axis + 1) % 3];
        let b = edges[(axis + 2) % 3];
        let area = [
            4.0 * (a[1] * b[2] - a[2] * b[1]),
            4.0 * (a[2] * b[0] - a[0] * b[2]),
            4.0 * (a[0] * b[1] - a[1] * b[0]),
        ];
        let mut face_edges = edges;
        face_edges[axis] = [0.0; 3];
        let face = prepare(&field, density, face_edges);
        for sign in [-1.0, 1.0] {
            let point = core::array::from_fn(|i| center[i] + sign * edges[axis][i]);
            let mass = field
                .sample_anelastic_filtered_moments(point, &face, || false)
                .unwrap()
                .mean_mass_flux_kg_m2_s;
            let flux = sign * dot(area, mass);
            balance += flux;
            scale += flux.abs();
        }
    }
    assert!(balance.abs() < 2e-12 * scale);
    for kappa in [0.0, 0.2, 3.0] {
        let other = prepare(&field, super::tests::density(kappa), edges);
        assert_eq!(
            whole.mean_mass_flux_kg_m2_s,
            field
                .sample_anelastic_filtered_moments(center, &other, || false)
                .unwrap()
                .mean_mass_flux_kg_m2_s
        );
    }
}

#[test]
fn constant_density_and_points_preserve_original_arithmetic() {
    let field = field();
    let center = [0.3, -0.6, 10.0];
    for edges in [
        [[0.0; 3]; 3],
        [[0.2, 0.1, 0.0], [0.0; 3], [0.0; 3]],
        [[0.2, 0.1, 0.1], [0.0; 3], [0.0, 0.0, 0.3]],
    ] {
        let original =
            MannBoxMoments::prepare(field.model(), edges, 100_000, 1000, || false).unwrap();
        let new = prepare(&field, density(0.0), edges);
        assert_eq!(new.workspace_bytes(), original.workspace_bytes());
        let a = field
            .sample_filtered_moments(center, &original, || false)
            .unwrap();
        let b = field
            .sample_anelastic_filtered_moments(center, &new, || false)
            .unwrap();
        assert_eq!(a.mean, b.favre_mean);
        assert_eq!(a.subfilter_stress_m2_s2, b.favre_subfilter_stress_m2_s2);
    }
    let density = density(0.8);
    let point =
        MannAnelasticBoxMoments::prepare(field.model(), density, [[0.0; 3]; 3], 0, 0, || false)
            .unwrap();
    let original = field.sample_anelastic(center, &density).unwrap();
    let got = field
        .sample_anelastic_filtered_moments(center, &point, || false)
        .unwrap();
    assert_eq!(got.favre_mean, original.flow);
    assert_eq!(got.mean_mass_flux_kg_m2_s, original.mass_flux_kg_m2_s);
    assert_eq!(got.favre_subfilter_stress_m2_s2, [[0.0; 3]; 3]);
}

#[test]
fn source_binding_budgets_cancellation_and_both_clock_replays_are_preserved() {
    let field = field();
    let center = [0.3, -0.6, 10.0];
    let edges = [[0.2, 0.1, 0.1], [0.0; 3], [0.0, 0.0, 0.3]];
    let observation = prepare(&field, density(0.8), edges);
    let before = field.checkpoint();
    let expected = field
        .sample_anelastic_filtered_moments(center, &observation, || false)
        .unwrap();
    let bytes = 12 * 8 + 78 * 32;
    assert_eq!(observation.workspace_bytes(), bytes);
    assert_eq!(observation.mode_pairs(), 78);
    for (limit, pairs) in [(bytes - 1, 78), (bytes, 77)] {
        assert!(
            MannAnelasticBoxMoments::prepare(
                field.model(),
                density(0.8),
                edges,
                limit,
                pairs,
                || false
            )
            .is_err()
        );
    }
    assert!(
        MannAnelasticBoxMoments::prepare(field.model(), density(0.8), edges, bytes, 78, || false)
            .is_ok()
    );
    let mut calls = 0;
    field
        .sample_anelastic_filtered_moments(center, &observation, || {
            calls += 1;
            false
        })
        .unwrap();
    for at in [1, 5, calls] {
        let mut count = 0;
        assert!(
            field
                .sample_anelastic_filtered_moments(center, &observation, || {
                    count += 1;
                    count == at
                })
                .is_err()
        );
        assert_eq!(field.checkpoint(), before);
        assert_eq!(
            field
                .sample_anelastic_filtered_moments(center, &observation, || false)
                .unwrap(),
            expected
        );
    }
    let mut count = 0;
    assert!(
        MannAnelasticBoxMoments::prepare(field.model(), density(0.8), edges, bytes, 78, || {
            count += 1;
            count == 8
        })
        .is_err()
    );
    assert!(
        field
            .sample_anelastic_filtered_moments([0.0, 0.0, 0.2], &observation, || false)
            .is_err()
    );
    let old = MannSpectralModel::build(
        17,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 12,
            min_kl: 0.05,
            max_kl: 3.0,
        },
        [2.0, -0.3, 0.0],
    )
    .unwrap();
    let wrong =
        MannAnelasticBoxMoments::prepare(&old, density(0.8), edges, bytes, 78, || false).unwrap();
    assert_eq!(
        field
            .sample_anelastic_filtered_moments(center, &wrong, || false)
            .unwrap_err()
            .code,
        "mann-filter-source-mismatch"
    );
    let restored = HomogeneousMannField::resume(before.clone()).unwrap();
    assert_eq!(
        restored
            .sample_anelastic_filtered_moments(center, &observation, || false)
            .unwrap(),
        expected
    );
    let mut fixed = FixedStepMannField::from_snapshot(field, 191, 0.017, 88 * 12).unwrap();
    assert_eq!(
        fixed
            .sample_anelastic_filtered_moments(center, &observation, || false)
            .unwrap(),
        expected
    );
    fixed.advance_to(3, 4 * 12 * 3, || false).unwrap();
    let checkpoint = fixed.checkpoint();
    let mut restored = FixedStepMannField::resume(checkpoint.clone(), 88 * 12).unwrap();
    assert_eq!(
        fixed.sample_anelastic_filtered_moments(center, &observation, || false),
        restored.sample_anelastic_filtered_moments(center, &observation, || false)
    );
    assert!(
        fixed
            .sample_anelastic_filtered_moments(center, &observation, || true)
            .is_err()
    );
    assert_eq!(fixed.checkpoint(), checkpoint);
    fixed.advance_to(5, 4 * 12 * 2, || false).unwrap();
    restored.advance_to(5, 4 * 12 * 2, || false).unwrap();
    assert_eq!(fixed, restored);
}
