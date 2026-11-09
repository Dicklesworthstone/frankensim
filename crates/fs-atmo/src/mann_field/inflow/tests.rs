use super::*;
use crate::mann::MANN_TARGET_V1;
use crate::mann_field::{
    HomogeneousMannField, MannQuadrature, MannQuadratureRule, fixed_step::FixedStepMannField,
};

fn source() -> HomogeneousMannField {
    let model = MannSpectralModel::build_with_rule(
        43,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 8,
            min_kl: 0.05,
            max_kl: 3.0,
        },
        [7.8, -0.9, 0.2],
        MannQuadratureRule::BalancedEnergyV1,
    )
    .unwrap();
    let mut source = HomogeneousMannField::at_anchor(model, 173, vec![0.7; 8]).unwrap();
    source.advance_to(source.tick() + 7).unwrap();
    source
}

fn patch() -> MannInflowPatch {
    MannInflowPatch::new([0.3, -0.7, 4.0], [[0.2, 0.7, 0.1], [-0.1, 0.15, 0.65]]).unwrap()
}

fn frame(source: &HomogeneousMannField) -> MannInflowFrame {
    source
        .inflow_frame(
            patch(),
            MannInflowFrame::required_workspace_bytes(source.model().mode_count()).unwrap(),
            || false,
        )
        .unwrap()
}

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() < tol * (1.0 + a.abs().max(b.abs())),
        "{a:e} vs {b:e}"
    );
}

fn same_bits(a: FlowSample, b: FlowSample) {
    for i in 0..3 {
        assert_eq!(a.u[i].to_bits(), b.u[i].to_bits());
        for j in 0..3 {
            assert_eq!(a.grad[i][j].to_bits(), b.grad[i][j].to_bits());
        }
    }
}

// An independent Gauss-Legendre surface integral; no prepared face gain enters.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    let mut result = Vec::new();
    for i in 0..n {
        let mut x = (core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 0.0;
        for _ in 0..30 {
            let (mut a, mut b) = (1.0, x);
            for j in 2..=n {
                let c = ((2 * j - 1) as f64 * x * b - (j - 1) as f64 * a) / j as f64;
                a = b;
                b = c;
            }
            derivative = n as f64 * (x * b - a) / (x * x - 1.0);
            let delta = b / derivative;
            x -= delta;
            if delta.abs() < 2e-16 {
                break;
            }
        }
        result.push((x, 1.0 / ((1.0 - x * x) * derivative * derivative)));
    }
    result
}

fn surface_integral(
    patch: MannInflowPatch,
    mut sample: impl FnMut([f64; 3]) -> [f64; 3],
) -> [f64; 3] {
    let rule = gauss(16);
    let mut integral = [0.0; 3];
    for &(a, wa) in &rule {
        for &(b, wb) in &rule {
            let point = core::array::from_fn(|i| {
                patch.center_m[i] + a * patch.half_edges_m[0][i] + b * patch.half_edges_m[1][i]
            });
            let value = sample(point);
            for i in 0..3 {
                integral[i] += wa * wb * value[i];
            }
        }
    }
    integral
}

fn matrix_oracle(
    model: &MannSpectralModel,
    time: f64,
    amplitudes: [&[f64]; 4],
    frame: &mut MannInflowFrame,
) {
    let n = model.mode_count();
    let dofs = 4 * n;
    let patch = frame.patch;
    let points = [[0.2, -0.4, 3.8], [0.9, 0.2, 4.4]];
    let mut constraint = vec![0.0; dofs];
    let mut columns = [vec![[0.0; 3]; dofs], vec![[0.0; 3]; dofs]];
    let mut basis: [Vec<f64>; 4] = core::array::from_fn(|_| vec![0.0; n]);
    for j in 0..dofs {
        basis[j / n][j % n] = 1.0;
        let paths = core::array::from_fn(|p| basis[p].as_slice());
        constraint[j] = dot(
            surface_integral(patch, |x| {
                sample_modes(model, time, paths, x, None).unwrap().u
            }),
            patch.normal,
        );
        for p in 0..2 {
            columns[p][j] = sample_modes(model, time, paths, points[p], None).unwrap().u;
        }
        basis[j / n][j % n] = 0.0;
    }
    let variance = constraint.iter().map(|x| x * x).sum::<f64>();
    assert!(variance > 1e-5);
    let residual = (0..dofs)
        .map(|j| constraint[j] * amplitudes[j / n][j % n])
        .sum::<f64>();
    refresh(model, time, amplitudes, frame, || false).unwrap();
    for j in 0..dofs {
        close(
            frame.amplitudes[j / n][j % n],
            amplitudes[j / n][j % n] - constraint[j] * residual / variance,
            3e-12,
        );
    }
    let cross: [[f64; 3]; 2] = core::array::from_fn(|p| {
        core::array::from_fn(|i| (0..dofs).map(|j| columns[p][j][i] * constraint[j]).sum())
    });
    let mut measured = [[[[0.0; 3]; 3]; 2]; 2];
    for j in 0..dofs {
        basis[j / n][j % n] = 1.0;
        refresh(
            model,
            time,
            core::array::from_fn(|p| basis[p].as_slice()),
            frame,
            || false,
        )
        .unwrap();
        let observed = points.map(|x| frame.sample(x).unwrap().u);
        for p in 0..2 {
            for q in 0..2 {
                for a in 0..3 {
                    for b in 0..3 {
                        measured[p][q][a][b] += observed[p][a] * observed[q][b];
                    }
                }
            }
        }
        basis[j / n][j % n] = 0.0;
    }
    let mut removed = 0.0;
    for p in 0..2 {
        for q in 0..2 {
            for a in 0..3 {
                for b in 0..3 {
                    let original = (0..dofs)
                        .map(|j| columns[p][j][a] * columns[q][j][b])
                        .sum::<f64>();
                    let expected = original - cross[p][a] * cross[q][b] / variance;
                    close(measured[p][q][a][b], expected, 5e-12);
                    if p == q && a == b {
                        removed += original - measured[p][q][a][b];
                    }
                }
            }
        }
    }
    assert!(removed > 0.01); // The conditional covariance really changes.
    refresh(model, time, amplitudes, frame, || false).unwrap();
}

#[test]
fn full_conditional_covariance_matches_an_independent_surface_and_basis_matrix() {
    let source = source();
    let mut frame = frame(&source);
    matrix_oracle(
        source.model(),
        source.tick() as f64 / crate::TICK_HZ,
        core::array::from_fn(|i| source.paths[i].amplitudes()),
        &mut frame,
    );
    let mut fixed = FixedStepMannField::from_snapshot(source, 239, 0.017, 88 * 8).unwrap();
    fixed.advance_to(5, 4 * 8 * 5, || false).unwrap();
    let checkpoint = fixed.checkpoint();
    fixed.refresh_inflow_frame(&mut frame, || false).unwrap();
    matrix_oracle(
        fixed.model(),
        fixed.time_s(),
        core::array::from_fn(|i| checkpoint.amplitudes[i].as_slice()),
        &mut frame,
    );
    assert_eq!(frame.time_s(), fixed.time_s());
    assert_eq!(checkpoint, fixed.checkpoint());
}

#[test]
fn full_and_subdivided_inlets_preserve_bulk_mass_flux_and_turbulent_momentum() {
    let source = source();
    let frame = frame(&source);
    let patch = patch();
    let whole = MannBoxFilter::prepare(source.model(), patch.footprint(), 8 * 8, || false).unwrap();
    assert!(
        dot(
            source.sample_filtered(patch.center_m, &whole).unwrap().u,
            patch.normal
        )
        .abs()
            > 1e-5
    );
    let flow = frame.sample_filtered(patch.center_m, &whole).unwrap();
    assert!(dot(flow.u, patch.normal).abs() < 3e-13);
    let integrated = surface_integral(patch, |x| frame.sample(x).unwrap().u);
    for i in 0..3 {
        close(integrated[i], flow.u[i], 2e-12);
    }
    let mut tile_edges = patch.footprint();
    tile_edges[0] = tile_edges[0].map(|x| x / 4.0);
    tile_edges[1] = tile_edges[1].map(|x| x / 3.0);
    let tile = MannBoxFilter::prepare(source.model(), tile_edges, 8 * 8, || false).unwrap();
    let centers: Vec<[f64; 3]> = (0..4)
        .flat_map(|i| {
            (0..3).map(move |j| {
                core::array::from_fn(|k| {
                    patch.center_m[k]
                        + ((2 * i + 1) as f64 / 4.0 - 1.0) * patch.half_edges_m[0][k]
                        + ((2 * j + 1) as f64 / 3.0 - 1.0) * patch.half_edges_m[1][k]
                })
            })
        })
        .collect();
    let mut tiled = [0.0; 3];
    for point in &centers {
        let cell = frame.sample_filtered(*point, &tile).unwrap();
        for i in 0..3 {
            tiled[i] += cell.u[i] / 12.0;
        }
    }
    for i in 0..3 {
        close(tiled[i], flow.u[i], 2e-13);
    }
    assert!(dot(tiled, patch.normal).abs() < 3e-13);
    for kappa in [0.0, 0.4, 2.5] {
        let density = ExponentialDensityProfile::new(1.2, 4.0, kappa, [0.0, 10.0]).unwrap();
        let full = filter::anelastic::MannAnelasticBoxMoments::prepare(
            source.model(),
            density,
            patch.footprint(),
            100_000,
            100,
            || false,
        )
        .unwrap();
        let cells = filter::anelastic::MannAnelasticBoxMoments::prepare(
            source.model(),
            density,
            tile_edges,
            100_000,
            100,
            || false,
        )
        .unwrap();
        let expected = frame
            .sample_anelastic_filtered_moments(patch.center_m, &full, || false)
            .unwrap();
        assert!(dot(expected.mean_mass_flux_kg_m2_s, patch.normal).abs() < 5e-13);
        let integrated = surface_integral(patch, |x| {
            frame
                .sample_anelastic(x, &density)
                .unwrap()
                .mass_flux_kg_m2_s
        });
        assert!(dot(integrated, patch.normal).abs() < 3e-12);
        let mut mass = [0.0; 3];
        let mut momentum = [[0.0; 3]; 3];
        let mut energy = 0.0;
        for &point in &centers {
            let cell = frame
                .sample_anelastic_filtered_moments(point, &cells, || false)
                .unwrap();
            energy += cell.kinetic_energy_density_j_m3 / 12.0;
            for i in 0..3 {
                mass[i] += cell.mean_mass_flux_kg_m2_s[i] / 12.0;
                for j in 0..3 {
                    momentum[i][j] += cell.convective_momentum_flux_pa[i][j] / 12.0;
                }
            }
        }
        for i in 0..3 {
            close(mass[i], expected.mean_mass_flux_kg_m2_s[i], 2e-13);
            for j in 0..3 {
                close(
                    momentum[i][j],
                    expected.convective_momentum_flux_pa[i][j],
                    3e-12,
                );
            }
        }
        close(energy, expected.kinetic_energy_density_j_m3, 3e-12);
        assert!(
            (0..3)
                .map(|i| expected.favre_subfilter_stress_m2_s2[i][i])
                .sum::<f64>()
                > 0.001
        );
    }
}

#[test]
fn projected_local_gradients_and_closed_surface_flux_remain_solenoidal() {
    let source = source();
    let frame = frame(&source);
    for point in [[0.1, 0.2, 3.0], [0.8, -0.9, 4.3], [-1.2, 0.6, 5.7]] {
        let sample = frame.sample(point).unwrap();
        assert!((sample.grad[0][0] + sample.grad[1][1] + sample.grad[2][2]).abs() < 2e-13);
        for j in 0..3 {
            let mut a = point;
            let mut b = point;
            a[j] += 1e-5;
            b[j] -= 1e-5;
            let plus = frame.sample(a).unwrap().u;
            let minus = frame.sample(b).unwrap().u;
            for i in 0..3 {
                close(sample.grad[i][j], (plus[i] - minus[i]) / 2e-5, 5e-9);
            }
        }
    }
    let center = [0.5, 0.3, 4.5];
    let edges = [[0.6, 0.1, 0.2], [-0.1, 0.4, 0.1], [0.2, 0.1, 0.7]];
    let mut flux = 0.0;
    let mut scale = 0.0;
    for axis in 0..3 {
        let a = edges[(axis + 1) % 3];
        let b = edges[(axis + 2) % 3];
        let area = [
            4.0 * (a[1] * b[2] - a[2] * b[1]),
            4.0 * (a[2] * b[0] - a[0] * b[2]),
            4.0 * (a[0] * b[1] - a[1] * b[0]),
        ];
        let mut face = edges;
        face[axis] = [0.0; 3];
        let filter = MannBoxFilter::prepare(source.model(), face, 64, || false).unwrap();
        for sign in [-1.0, 1.0] {
            let point = core::array::from_fn(|i| center[i] + sign * edges[axis][i]);
            let value = sign * dot(area, frame.sample_filtered(point, &filter).unwrap().u);
            flux += value;
            scale += value.abs();
        }
    }
    assert!(flux.abs() < 2e-13 * scale);
}

#[test]
fn both_clocks_replay_and_cancelled_refresh_preserve_the_complete_accepted_frame() {
    let mut source = source();
    let before = source.checkpoint();
    let mut frame = frame(&source);
    assert_eq!(before, source.checkpoint());
    let saved = frame.clone();
    let point = [0.8, -0.9, 4.3];
    let sample = saved.sample(point).unwrap();
    source.advance_to(source.tick() + 9).unwrap();
    let accepted_source = source.checkpoint();
    let expected = source
        .inflow_frame(patch(), frame.workspace_bytes(), || false)
        .unwrap();
    let mut total = 0;
    source
        .refresh_inflow_frame(&mut frame, || {
            total += 1;
            false
        })
        .unwrap();
    for stop in [1, 5, total] {
        frame = saved.clone();
        let mut count = 0;
        assert!(
            source
                .refresh_inflow_frame(&mut frame, || {
                    count += 1;
                    count == stop
                })
                .is_err()
        );
        assert_eq!(frame, saved);
        same_bits(frame.sample(point).unwrap(), sample);
        assert_eq!(source.checkpoint(), accepted_source);
        source.refresh_inflow_frame(&mut frame, || false).unwrap();
        assert_eq!(frame, expected);
    }
    let restored = HomogeneousMannField::resume(accepted_source).unwrap();
    assert_eq!(
        restored
            .inflow_frame(patch(), frame.workspace_bytes(), || false)
            .unwrap(),
        frame
    );
    let mut fixed = FixedStepMannField::from_snapshot(source, 239, 0.017, 88 * 8).unwrap();
    assert_eq!(
        fixed
            .inflow_frame(patch(), frame.workspace_bytes(), || false)
            .unwrap(),
        frame
    );
    fixed.advance_to(3, 4 * 8 * 3, || false).unwrap();
    fixed.refresh_inflow_frame(&mut frame, || false).unwrap();
    let checkpoint = fixed.checkpoint();
    let saved = frame.clone();
    let mut resumed = FixedStepMannField::resume(checkpoint.clone(), 88 * 8).unwrap();
    assert_eq!(
        resumed
            .inflow_frame(patch(), frame.workspace_bytes(), || false)
            .unwrap(),
        frame
    );
    fixed.advance_to(5, 4 * 8 * 2, || false).unwrap();
    let mut count = 0;
    assert!(
        fixed
            .refresh_inflow_frame(&mut frame, || {
                count += 1;
                count == 7
            })
            .is_err()
    );
    assert_eq!(frame, saved);
    fixed.refresh_inflow_frame(&mut frame, || false).unwrap();
    resumed.advance_to(5, 4 * 8 * 2, || false).unwrap();
    assert_eq!(
        resumed
            .inflow_frame(patch(), frame.workspace_bytes(), || false)
            .unwrap(),
        frame
    );
    assert_eq!(fixed, resumed);
}

#[test]
fn storage_geometry_and_source_binding_are_explicit() {
    let source = source();
    let bytes = MannInflowFrame::required_workspace_bytes(8).unwrap();
    let frame = source.inflow_frame(patch(), bytes, || false).unwrap();
    assert_eq!(frame.workspace_bytes(), bytes);
    assert!(source.inflow_frame(patch(), bytes - 1, || false).is_err());
    assert!(MannInflowFrame::required_workspace_bytes(0).is_err());
    assert!(MannInflowFrame::required_workspace_bytes(MAX_MODES + 1).is_err());
    assert!(MannInflowFrame::required_workspace_bytes(MAX_MODES).is_ok());
    assert!(source.inflow_frame(patch(), bytes, || true).is_err());
    for edges in [
        [[0.0; 3], [1.0, 0.0, 0.0]],
        [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        [[1.0, f64::NAN, 0.0], [0.0, 1.0, 0.0]],
        [[f64::MAX, 0.0, 0.0], [0.0, f64::MAX, 0.0]],
    ] {
        assert!(MannInflowPatch::new([0.0; 3], edges).is_err());
    }
    assert!(MannInflowPatch::new([f64::INFINITY, 0.0, 0.0], patch().half_edges_m()).is_err());
    let model = MannSpectralModel::build(
        43,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 8,
            min_kl: 0.05,
            max_kl: 3.0,
        },
        [7.8, -0.9, 0.2],
    )
    .unwrap();
    let other = HomogeneousMannField::at_anchor(model, 173, vec![0.7; 8]).unwrap();
    let mut unchanged = frame.clone();
    assert_eq!(
        other
            .refresh_inflow_frame(&mut unchanged, || false)
            .unwrap_err()
            .code,
        "mann-filter-source-mismatch"
    );
    assert_eq!(unchanged, frame);
    let rectangular = MannInflowPatch::new([0.0; 3], [[0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]).unwrap();
    assert_eq!(rectangular.normal(), [1.0, 0.0, 0.0]);
    assert_eq!(rectangular.area_m2(), 24.0);
}

#[test]
fn an_identically_zero_flux_functional_is_an_exact_identity() {
    // Manufactured transverse geometry exercises the zero-variance constraint.
    let mut model = source().model().clone();
    for mode in &mut model.modes {
        mode.k = [0.0, 1.0, 0.0];
        mode.polarization = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let source = HomogeneousMannField::at_anchor(model, 173, vec![0.7; 8]).unwrap();
    let patch = MannInflowPatch::new([0.0; 3], [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]).unwrap();
    let frame = source.inflow_frame(patch, 10_000, || false).unwrap();
    assert!(!frame.active_constraint);
    for point in [[0.0; 3], [0.3, -0.7, 4.0], [2.0, 0.8, 7.0]] {
        same_bits(frame.sample(point).unwrap(), source.sample(point).unwrap());
    }
}
