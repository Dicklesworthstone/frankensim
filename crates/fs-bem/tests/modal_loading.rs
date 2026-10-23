//! M2.5-W4a: complete modal loading. Strict sphere tolerances are qualification
//! gates, not a claim that the existing collocation owner already meets them.
use fs_bem::helmholtz::rayleigh_modal::{self, AreaSample, PreparedSurface};
use fs_bem::helmholtz::{self, Formulation, IterativeOptions, Medium};
use fs_bem::modal_loading::*;
use fs_bem::panel3d::SpherePanels;
use fs_blake3::hash_domain;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};
use std::time::Instant;

fn air() -> Medium {
    Medium::air()
}
fn source() -> fs_blake3::ContentHash {
    hash_domain("analytic-fixture", "M2.5-W4a".as_bytes())
}
fn planar(cells: usize, modes: usize) -> Vec<AreaSample> {
    (0..cells)
        .map(|i| AreaSample {
            position_m: [0.02 * (i % 4) as f64, 0.025 * (i / 4) as f64],
            area_m2: 2e-4,
            mode_shape: (0..modes)
                .map(|m| {
                    if m == 0 {
                        1.0
                    } else {
                        ((i + 1) * (m + 1)) as f64 / 10.0 - 0.4
                    }
                })
                .collect(),
        })
        .collect()
}
fn load_planar(cells: &[AreaSample], omega: &[f64]) -> Result<ModalLoadingBake, ModalLoadingError> {
    modal_radiation_loading(
        SurfaceModeTraces::PlanarBaffled { samples: cells },
        LoadingRoute::PlanarBaffledRayleigh,
        omega,
        air(),
        source(),
        1_000_000,
        || false,
    )
}

/// G0: diagnostics belong to the original selected solve, including the
/// low-frequency CBIE route and the resonance-safe Burton--Miller route.
#[test]
fn bem_diagnostics_retain_the_selected_original_solve() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let medium = air();
    let shapes = vec![vec![1.0; surface.areas().len()]];
    let field = vec![C64::ONE; surface.areas().len()];
    let options = IterativeOptions {
        panel_rule: helmholtz::PanelRule::Triangle {
            max_pairs: surface.areas().len().pow(2),
        },
        ..IterativeOptions::default()
    };
    let omegas = [0.05 * medium.sound_speed, 2.5 * medium.sound_speed];
    let bake = modal_radiation_loading_routed(
        SurfaceModeTraces::ClosedBody {
            surface: &surface,
            shapes: &shapes,
        },
        LoadingRoute::ClosedBodyGeometryPolicy,
        &omegas,
        medium,
        source(),
        2,
        Some(&options),
        || false,
    )
    .unwrap();
    for (row, formulation) in bake
        .rows()
        .iter()
        .zip([Formulation::PlainCbie, Formulation::BurtonMiller])
    {
        let solution = helmholtz::solve_radiation_batch_with_rule(
            &surface,
            row.omega() / medium.sound_speed,
            medium,
            &[&field],
            formulation,
            options.panel_rule,
        )
        .unwrap();
        let diagnostics = row.bem_diagnostics().unwrap();
        assert_eq!(diagnostics.formulation, formulation);
        assert_eq!(
            diagnostics.minimum_panels_per_wavelength.to_bits(),
            solution[0].panels_per_wavelength.to_bits()
        );
        assert_eq!(
            diagnostics.maximum_condition_lower_bound.to_bits(),
            solution[0].condition_lower_bound.to_bits()
        );
    }
}

#[test]
fn modal_loading_caps() {
    let bem = LoadingRoute::ClosedBodyBem(Formulation::BurtonMiller);
    let ray = LoadingRoute::PlanarBaffledRayleigh;
    assert_eq!(admit_loading_work(8192, 256, 1, bem, 256 * 256).unwrap(), 1);
    assert!(matches!(
        admit_loading_work(8193, 1, 1, bem, 1),
        Err(ModalLoadingError::SolverRoute {
            route: "M2.5-W5c",
            ..
        })
    ));
    assert_eq!(admit_loading_work(20, 257, 1, bem, 257 * 257).unwrap(), 2);
    assert_eq!(source_batches(257).collect::<Vec<_>>(), [0..256, 256..257]);
    assert_eq!(admit_loading_work(1, 128, 1, ray, 128 * 128).unwrap(), 1);
    assert!(admit_loading_work(1, 129, 1, ray, 129 * 129).is_err());
    rayleigh_modal::admit_modal_work(u128::from(rayleigh_modal::MAX_MODAL_RAYLEIGH_WORK)).unwrap();
    assert!(
        rayleigh_modal::admit_modal_work(u128::from(rayleigh_modal::MAX_MODAL_RAYLEIGH_WORK) + 1)
            .is_err()
    );
    assert!(admit_loading_work(1, 2, 3, ray, 11).is_err());
    assert!(admit_loading_work(1, usize::MAX, 2, bem, usize::MAX).is_err());
    assert!(
        admit_loading_work(
            20,
            1,
            1,
            LoadingRoute::ClosedBodyBem(Formulation::PlainCbie),
            1
        )
        .is_err()
    );
    assert!(
        admit_loading_work(
            20,
            1,
            1,
            LoadingRoute::ClosedBodyBem(Formulation::BurtonMillerWrongAlphaSign),
            1
        )
        .is_err()
    );

    // Real Rayleigh preparation at the mode-count boundary, not only a planner.
    let cells = planar(1, 128);
    assert_eq!(load_planar(&cells, &[100.0]).unwrap().modes(), 128);
    assert!(load_planar(&planar(1, 129), &[100.0]).is_err());
    let cells = [AreaSample {
        position_m: [0.0; 2],
        area_m2: PI,
        mode_shape: vec![1.0],
    }];
    let omega6 = PI * air().sound_speed / 6.0;
    let prepared = PreparedSurface::new(&cells).unwrap();
    assert!(prepared.power_consistent_impedance(omega6, air()).is_ok());
    assert!(matches!(
        prepared.power_consistent_impedance(omega6 * (1.0 + 1e-9), air()),
        Err(helmholtz::HelmholtzError::TooCoarse { .. })
    ));

    let sphere = SpherePanels::icosphere(1.0, 0).unwrap();
    let size = sphere
        .areas()
        .iter()
        .copied()
        .fold(0.0_f64, f64::max)
        .sqrt();
    let k6 = TAU / (6.0 * size);
    let field = vec![C64::ONE; sphere.areas().len()];
    assert!(
        helmholtz::solve_radiation_batch(&sphere, k6, air(), &[&field], Formulation::BurtonMiller)
            .is_ok()
    );
    assert!(matches!(
        helmholtz::solve_radiation_batch(
            &sphere,
            k6 * (1.0 + 1e-9),
            air(),
            &[&field],
            Formulation::BurtonMiller
        ),
        Err(helmholtz::HelmholtzError::TooCoarse { .. })
    ));
}

#[test]
fn modal_loading_passivity() {
    let started = Instant::now();
    let cells = planar(12, 4);
    let bake = load_planar(&cells, &[50.0, 500.0, 2500.0]).unwrap();
    for row in bake.rows() {
        assert!(row.bem_diagnostics().is_none());
        let norm = row
            .impedance()
            .iter()
            .fold(0.0_f64, |s, z| s.hypot(z.abs()));
        assert!(row.minimum_resistance_eigenvalue() >= -1e-12 * norm);
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"modal-loading-passivity\",\"inputs\":\"analytic signed planar fields; physical area .0002m2/cell\",\"omega\":{},\"observed\":{},\"oracle\":\"minimum resistance eigenvalue >= -1e-12 Frobenius norm\",\"tol\":{},\"verdict\":\"pass\",\"elapsed_s\":{}}}",
            row.omega(),
            row.minimum_resistance_eigenvalue(),
            1e-12 * norm,
            started.elapsed().as_secs_f64()
        );
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(row.impedance()[i * 4 + j], row.impedance()[j * 4 + i]);
            }
        }
    }
    assert!(load_planar(&cells, &[500.0, 50.0]).is_err());
    assert!(load_planar(&cells, &[50.0, 50.0]).is_err());
    assert!(load_planar(&cells, &[f64::NAN]).is_err());
    let mut invalid = cells.clone();
    invalid[0].mode_shape.pop();
    assert!(load_planar(&invalid, &[50.0]).is_err());
    let mut invalid = cells.clone();
    invalid[0].mode_shape[0] = f64::INFINITY;
    assert!(load_planar(&invalid, &[50.0]).is_err());
    assert!(bake.at_omega(51.0).is_err());
    assert!(matches!(
        modal_radiation_loading(
            SurfaceModeTraces::PlanarBaffled { samples: &cells },
            LoadingRoute::PlanarBaffledRayleigh,
            &[50.0],
            air(),
            source(),
            100,
            || true
        ),
        Err(ModalLoadingError::Helmholtz(
            helmholtz::HelmholtzError::Cancelled
        ))
    ));
    assert!(
        modal_radiation_loading(
            SurfaceModeTraces::PlanarBaffled { samples: &cells },
            LoadingRoute::ClosedBodyBem(Formulation::BurtonMiller),
            &[50.0],
            air(),
            source(),
            100,
            || false
        )
        .is_err()
    );
}

#[test]
fn modal_loading_bitwise_repeat() {
    let cells = planar(8, 3);
    let a = load_planar(&cells, &[100.0, 900.0]).unwrap();
    let b = load_planar(&cells, &[100.0, 900.0]).unwrap();
    assert_eq!(a.declaration_hash(), b.declaration_hash());
    for (x, y) in a.rows().iter().zip(b.rows()) {
        assert_eq!(x.omega().to_bits(), y.omega().to_bits());
        assert_eq!(
            x.minimum_resistance_eigenvalue().to_bits(),
            y.minimum_resistance_eigenvalue().to_bits()
        );
        for (z, w) in x.impedance().iter().zip(y.impedance()) {
            assert_eq!(z.re.to_bits(), w.re.to_bits());
            assert_eq!(z.im.to_bits(), w.im.to_bits());
        }
    }
    let mut changed = cells.clone();
    changed[0].mode_shape[1] += 0.01;
    assert_ne!(
        a.declaration_hash(),
        load_planar(&changed, &[100.0, 900.0])
            .unwrap()
            .declaration_hash()
    );
    let changed_medium = modal_radiation_loading(
        SurfaceModeTraces::PlanarBaffled { samples: &cells },
        LoadingRoute::PlanarBaffledRayleigh,
        &[100.0, 900.0],
        Medium {
            density: 2.0 * air().density,
            ..air()
        },
        source(),
        1_000_000,
        || false,
    )
    .unwrap();
    assert_ne!(a.declaration_hash(), changed_medium.declaration_hash());
}

#[test]
fn bem_source_batches_preserve_all_257_fields_bitwise() {
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let mut shapes = vec![vec![0.0; surface.areas().len()]; 257];
    shapes[0].fill(1.0);
    shapes[256].fill(-0.5);
    let omega = air().sound_speed;
    let bake = modal_radiation_loading(
        SurfaceModeTraces::ClosedBody {
            surface: &surface,
            shapes: &shapes,
        },
        LoadingRoute::ClosedBodyBem(Formulation::BurtonMiller),
        &[omega],
        air(),
        source(),
        257 * 257,
        || false,
    )
    .unwrap();
    let mut expected = vec![C64::ZERO; 257 * 257];
    for range in source_batches(257) {
        let fields: Vec<Vec<C64>> = shapes[range.clone()]
            .iter()
            .map(|r| r.iter().copied().map(C64::from_re).collect())
            .collect();
        let refs: Vec<_> = fields.iter().map(Vec::as_slice).collect();
        let solutions = helmholtz::solve_radiation_batch(
            &surface,
            1.0,
            air(),
            &refs,
            Formulation::BurtonMiller,
        )
        .unwrap();
        for (s, solution) in range.zip(solutions) {
            for (r, shape) in shapes.iter().enumerate() {
                for ((&phi, &p), &area) in shape.iter().zip(&solution.pressure).zip(surface.areas())
                {
                    expected[r * 257 + s] = expected[r * 257 + s] + p.scale(phi * area);
                }
            }
        }
    }
    for (got, want) in bake.rows()[0].impedance().iter().zip(expected) {
        assert_eq!(got.re.to_bits(), want.re.to_bits());
        assert_eq!(got.im.to_bits(), want.im.to_bits());
    }
}

/// Independent double sum of the declared Rayleigh quadrature, including its
/// explicit equivalent-disc singular reactance and point-field Gram diagonal.
/// This checks kernel signs, signed cross-mode terms and both area factors;
/// continuous-cell discretization error is tested separately by the plate ladder.
#[test]
fn bem_and_rayleigh_routes_agree_on_baffled_plate() {
    let started = Instant::now();
    let cells = planar(12, 3);
    let omega = 1200.0;
    let medium = air();
    let k = omega / medium.sound_speed;
    let bake = load_planar(&cells, &[omega]).unwrap();
    let actual = bake.rows()[0].impedance();
    let mut expected = vec![C64::ZERO; 9];
    for (i, a) in cells.iter().enumerate() {
        for (j, b) in cells.iter().enumerate() {
            let r = (a.position_m[0] - b.position_m[0]).hypot(a.position_m[1] - b.position_m[1]);
            let kernel = if i == j {
                C64::new(
                    medium.density * omega * omega / (TAU * medium.sound_speed),
                    -medium.density * medium.sound_speed * (k * (a.area_m2 / PI).sqrt()).sin()
                        / a.area_m2,
                )
            } else {
                C64::new((k * r).sin(), -(k * r).cos()).scale(medium.density * omega / (TAU * r))
            };
            for row in 0..3 {
                for col in 0..3 {
                    expected[row * 3 + col] = expected[row * 3 + col]
                        + kernel
                            .scale(a.area_m2 * b.area_m2 * a.mode_shape[row] * b.mode_shape[col]);
                }
            }
        }
    }
    let scale = expected.iter().fold(0.0_f64, |s, z| s.hypot(z.abs()));
    let error = actual
        .iter()
        .zip(expected)
        .fold(0.0_f64, |s, (a, b)| s.hypot((*a - b).abs()))
        / scale;
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"modal-loading-rayleigh-double-integral\",\"inputs\":\"analytic planar signed fields; same declared point/disc quadrature\",\"observed\":{error},\"oracle\":\"independent Rayleigh double sum\",\"tol\":1e-6,\"verdict\":\"{}\",\"elapsed_s\":{}}}",
        if error <= 1e-6 { "pass" } else { "fail" },
        started.elapsed().as_secs_f64()
    );
    assert!(error <= 1e-6);
}

// Integral_0^L sin(m*pi*x/L) exp(-i*q*x) dx, derived by the two exponentials.
fn sine_transform(m: usize, length: f64, q: f64) -> C64 {
    let phase = q * length;
    let mp = m as f64 * PI;
    if (phase.abs() - mp).abs() < 1e-10 {
        return C64::new(0.0, -phase.signum() * length / 2.0);
    }
    let sign = if m % 2 == 0 { 1.0 } else { -1.0 };
    C64::new(1.0 - sign * phase.cos(), sign * phase.sin())
        .scale(length * mp / (mp * mp - phase * phase))
}
fn plate_efficiency(m: usize, n: usize, a: f64, b: f64, k: f64) -> f64 {
    let (polar, azimuth) = (128, 256);
    let mut integral = 0.0;
    for i in 0..=polar {
        let mu = i as f64 / polar as f64;
        let transverse = (1.0 - mu * mu).sqrt();
        let weight = if i == 0 || i == polar {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for j in 0..azimuth {
            let angle = TAU * (j as f64 + 0.5) / azimuth as f64;
            let transform = sine_transform(m, a, k * transverse * angle.cos())
                * sine_transform(n, b, k * transverse * angle.sin());
            integral += weight * transform.norm_sq();
        }
    }
    integral *= TAU / (3.0 * polar as f64 * azimuth as f64);
    k * k * integral / (PI * PI * a * b)
}

/// G2: sine-mode far-field integral, the physical construction used by Wallace,
/// JASA 51 (1972) 946-952, doi:10.1121/1.1912943. Original figure/equation digits
/// could not be verified from the available NTRS abstract (19720042393), so this
/// is NOT represented as a digitized-paper comparison. No fitted digits enter.
#[test]
fn baffled_plate_radiation_efficiency() {
    let started = Instant::now();
    let (a, b) = (0.4, 0.3);
    let modes = [(1, 1), (2, 1), (1, 2), (3, 1), (2, 2), (3, 2)];
    let k = 18.0;
    let omega = k * air().sound_speed;
    let oracle: Vec<_> = modes
        .iter()
        .map(|&(m, n)| plate_efficiency(m, n, a, b, k))
        .collect();
    let mut errors = [0.0_f64; 3];
    for (rung, (nx, ny)) in [(12, 9), (24, 18), (48, 36)].into_iter().enumerate() {
        let mut cells = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                let x = a * (i as f64 + 0.5) / nx as f64;
                let y = b * (j as f64 + 0.5) / ny as f64;
                cells.push(AreaSample {
                    position_m: [x, y],
                    area_m2: a * b / (nx * ny) as f64,
                    mode_shape: modes
                        .iter()
                        .map(|&(m, n)| {
                            (m as f64 * PI * x / a).sin() * (n as f64 * PI * y / b).sin()
                        })
                        .collect(),
                });
            }
        }
        let bake = load_planar(&cells, &[omega]).unwrap();
        let row = &bake.rows()[0];
        for (mode, &reference) in oracle.iter().enumerate() {
            let actual = row.impedance()[mode * 6 + mode].re
                / (air().density * air().sound_speed * a * b / 4.0);
            let error = (actual / reference - 1.0).abs();
            errors[rung] = errors[rung].max(error);
            println!(
                "{{\"suite\":\"fs-bem\",\"case\":\"modal-loading-baffled-plate\",\"inputs\":\"analytic sine plate .4x.3m; Wallace DOI 10.1121/1.1912943; independent angular integral\",\"axis\":\"panels\",\"rung\":{rung},\"panels\":{},\"mode\":{mode},\"observed\":{actual},\"oracle\":{reference},\"err\":{error},\"min_eig_re\":{},\"tol\":0.02,\"verdict\":\"{}\",\"elapsed_s\":{}}}",
                cells.len(),
                row.minimum_resistance_eigenvalue(),
                if error <= 0.02 { "pass" } else { "fail" },
                started.elapsed().as_secs_f64()
            );
        }
    }
    assert!(errors[2] <= 0.02, "first six modes: {errors:?}");
    assert!(
        errors[2] < errors[1] && errors[1] < errors[0],
        "panel axis must be live: {errors:?}"
    );
}

// Shared latitude vertices keep the retained triangle surface exactly closed.
// Unlike 320/1280/5120 icospheres, all THREE rungs resolve ka=10 at >=6 ppw.
fn resolved_sphere(latitude: usize) -> SpherePanels {
    let longitude = 2 * latitude;
    let mut vertices = vec![[0.0, 0.0, 1.0]];
    for j in 1..latitude {
        let theta = PI * j as f64 / latitude as f64;
        for i in 0..longitude {
            let phi = TAU * i as f64 / longitude as f64;
            vertices.push([
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ]);
        }
    }
    let south = vertices.len();
    vertices.push([0.0, 0.0, -1.0]);
    let ring = |j: usize, i: usize| 1 + j * longitude + i % longitude;
    let mut triangles = Vec::new();
    for i in 0..longitude {
        triangles.push([vertices[0], vertices[ring(0, i)], vertices[ring(0, i + 1)]]);
        for j in 0..latitude - 2 {
            let (a, b, c, d) = (
                ring(j, i),
                ring(j + 1, i),
                ring(j + 1, i + 1),
                ring(j, i + 1),
            );
            triangles.push([vertices[a], vertices[b], vertices[c]]);
            triangles.push([vertices[a], vertices[c], vertices[d]]);
        }
        triangles.push([
            vertices[ring(latitude - 2, i)],
            vertices[south],
            vertices[ring(latitude - 2, i + 1)],
        ]);
    }
    SpherePanels::from_triangles(triangles).unwrap()
}

/// The M2.5-W5 options a routed sphere row declares. Workers and the M2L row
/// cache are scheduling budgets only: every solve is bit-identical at any
/// count. Leaf capacity 32 keeps plain CBIE's leaves at 5.3 to 7.6 panel
/// radii on the 16128-81224 rungs, above CBIE_MIN_LEAF_RADII. Root kD
/// exceeds W5a's 10 only at ka = 10, which takes the directional route. The
/// translation and coefficient-slot budgets admit each row's reciprocity
/// twins at FMM orders 8 and 10 (9 and 11 on the directional route); budgets
/// change no bits.
fn routed_options(ka: f64) -> IterativeOptions {
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get().min(16));
    let directional = (ka > 5.0).then(|| {
        let mut d = fs_fmm::DirectionalOptions::for_tolerance(1e-3).unwrap();
        d.order = 7;
        d.low_freq_threshold = 0.5;
        d.leaf_capacity = 128;
        d.max_coefficient_slots = 1 << 28;
        // Declared: the order-7 base alone needs 8.0006e9 basis-transfer
        // work on the 16128-panel rung, above the 8e9 default (ts1,
        // 2026-10-23), and the finer rungs and twins need more.
        d.max_translation_work = 1_000_000_000_000;
        d
    });
    IterativeOptions {
        tol: 1e-9,
        fmm_order: 6,
        max_iters: 900,
        leaf_capacity: 32,
        directional,
        max_near_entries: 1 << 29,
        near_field_cache_bytes: 1 << 31,
        max_translation_work: 1_000_000_000_000,
        max_coefficient_slots: 1 << 27,
        workers,
        m2l_cache_bytes: 1 << 33,
        ..IterativeOptions::default()
    }
}

/// G2, sphere radius 1 m. Derivation: p_l=A*h_l^(1)(kr)*P_l and
/// v_l=k*p_l'/(i*w*rho) give z_l=i*rho*c*h_l/h_l'. For l=0 this
/// reduces to (-i*x)/(1-i*x). For translation, integral cos²(theta)dS=S/3;
/// h_1=-exp(ix)*(x+i)/x² gives S*rho*c/3 * [x^4-i*(2*x+x³)]/(4+x^4).
/// Verified source: Sutherland WR68-12, NASA 19690015869, Appendix B,
/// equations B-9 (printed p43), B-18 (p45), and Fig.8 (p26).
///
/// On the M2.5-W5 route. The dense rungs alone (latitudes 24/32/40, up to
/// 6240 panels) end at the flat-panel floor: finest errors 2.11e-3, 2.03e-3,
/// 1.85e-3 and 2.63e-3 at ka 0.05, pi, 5 and 10 (ts1, 8476279df). The ladder
/// now runs latitudes 40, 64, 100 and 143 (6240 dense, then 16128, 39600 and
/// 81224 iterative). Each row keeps its geometry-policy formulation on every
/// rung: plain CBIE below 0.8 of the Dirichlet bound, Burton-Miller above.
#[test]
fn sphere_modal_loading_closed_forms() {
    let started = Instant::now();
    let medium = air();
    let area = 4.0 * PI;
    let latitudes = [40, 64, 100, 143];
    // Every ka is measured before the gate, so one run reports the band.
    let mut unmet = Vec::new();
    for ka in [0.05_f64, 0.5, 1.0, PI, 5.0, 10.0] {
        let x2 = ka * ka;
        let x4 = x2 * x2;
        let oracle = [
            C64::new(x2, -ka).scale(area * medium.density * medium.sound_speed / (1.0 + x2)),
            C64::new(x4, -ka * (2.0 + x2))
                .scale(area * medium.density * medium.sound_speed / (3.0 * (4.0 + x4))),
        ];
        let options = routed_options(ka);
        let mut finest = None;
        let mut admitted = 0;
        let mut ladder: Vec<(usize, [f64; 2])> = Vec::new();
        for (rung, latitude) in latitudes.into_iter().enumerate() {
            let surface = resolved_sphere(latitude);
            let panels = surface.areas().len();
            let shapes = vec![
                vec![1.0; panels],
                surface.normals().iter().map(|n| n[2]).collect(),
            ];
            let solve_started = Instant::now();
            let result = helmholtz::with_dense_workers(options.workers, || {
                modal_radiation_loading_routed(
                    SurfaceModeTraces::ClosedBody {
                        surface: &surface,
                        shapes: &shapes,
                    },
                    LoadingRoute::ClosedBodyGeometryPolicy,
                    &[ka * medium.sound_speed],
                    medium,
                    source(),
                    4,
                    Some(&options),
                    || false,
                )
            });
            let route = if panels > helmholtz::MAX_DENSE_PANELS {
                "iterative"
            } else {
                "dense"
            };
            // A refused row's matrix is measured and logged, never admitted.
            let (z, min_eig, refusal) = match &result {
                Ok(bake) => (
                    bake.rows()[0].impedance().to_vec(),
                    bake.rows()[0].minimum_resistance_eigenvalue(),
                    None,
                ),
                Err(
                    ModalLoadingError::NonReciprocal { impedance, .. }
                    | ModalLoadingError::NonPassive { impedance, .. },
                ) => (
                    impedance.clone(),
                    f64::NAN,
                    result.as_ref().err().map(ToString::to_string),
                ),
                Err(error) => {
                    println!(
                        "{{\"suite\":\"fs-bem\",\"case\":\"modal-loading-sphere\",\"inputs\":\"analytic sphere a=1m\",\"axis\":\"panels\",\"rung\":{rung},\"ka\":{ka},\"panels\":{panels},\"route\":\"{route}\",\"err\":\"{}\",\"tol\":1e-3,\"verdict\":\"refused\",\"elapsed_s\":{}}}",
                        error,
                        started.elapsed().as_secs_f64()
                    );
                    continue;
                }
            };
            if refusal.is_none() {
                admitted += 1;
            }
            let mut errors = [0.0_f64; 2];
            for mode in 0..2 {
                let error = (z[mode * 2 + mode] - oracle[mode]).abs() / oracle[mode].abs();
                errors[mode] = error;
                // Panel size scales as panels^-1/2.
                let order = ladder.last().map_or(f64::NAN, |&(n, e)| {
                    (e[mode] / error).ln() / (panels as f64 / n as f64).sqrt().ln()
                });
                let verdict = match (&refusal, error <= 1e-3) {
                    (Some(_), _) => "refused",
                    (None, true) => "pass",
                    (None, false) => "fail",
                };
                println!(
                    "{{\"suite\":\"fs-bem\",\"case\":\"modal-loading-sphere\",\"inputs\":\"analytic sphere a=1m; NASA WR68-12 B-9/B-18\",\"axis\":\"panels\",\"rung\":{rung},\"ka\":{ka},\"panels\":{panels},\"route\":\"{route}\",\"mode\":{mode},\"observed\":[{},{}],\"oracle\":[{},{}],\"err\":{error},\"order\":{order},\"min_eig_re\":{min_eig},\"tol\":1e-3,\"verdict\":\"{verdict}\",\"refusal\":\"{}\",\"solve_s\":{},\"elapsed_s\":{}}}",
                    z[mode * 2 + mode].re,
                    z[mode * 2 + mode].im,
                    oracle[mode].re,
                    oracle[mode].im,
                    refusal.as_deref().unwrap_or(""),
                    solve_started.elapsed().as_secs_f64(),
                    started.elapsed().as_secs_f64()
                );
            }
            ladder.push((panels, errors));
            if rung == latitudes.len() - 1 && refusal.is_none() {
                finest = Some(errors[0].max(errors[1]));
            }
        }
        if !(admitted == latitudes.len() && finest.is_some_and(|e| e <= 1e-3)) {
            unmet.push((ka, finest));
        }
    }
    assert!(
        unmet.is_empty(),
        "strict sphere qualification unmet at (ka, finest): {unmet:?}; refine the numerical owner, never relax 1e-3"
    );
}

/// G5/G0: a routed bake at or below the dense cap is the dense bake, bitwise,
/// declaration identity included, in both policy formulations. A declared
/// triangle rule reaches that bake and its declaration. Above the cap
/// the unrouted owner names M2.5-W5c, while the routed one solves.
/// - A one-mode (pulsating) bake is reciprocal by construction: it is
///   admitted, worker-bitwise, and binds its numerical options (not its
///   workers) into the declaration.
/// - The two-mode bake is reciprocal to its measured error scale (the .105
///   G2 ruling): the defect is at most max(10 eps_gmres, eps_fmm), with
///   eps_fmm from the order 8 twin, and the receipt records every term. Its
///   pulsating entry is the one-mode bake's, bitwise. Dense rows record the
///   1e-12 bound and no scale.
#[test]
fn routed_modal_loading_is_the_dense_bake_below_the_cap() {
    let medium = air();
    let below = SpherePanels::icosphere(1.0, 2).unwrap();
    let shapes_of = |surface: &SpherePanels| {
        vec![
            vec![1.0; surface.areas().len()],
            surface.normals().iter().map(|n| n[2]).collect::<Vec<f64>>(),
        ]
    };
    let bake_modes =
        |surface: &SpherePanels, modes: usize, ka: f64, options: Option<&IterativeOptions>| {
            let mut shapes = shapes_of(surface);
            shapes.truncate(modes);
            modal_radiation_loading_routed(
                SurfaceModeTraces::ClosedBody {
                    surface,
                    shapes: &shapes,
                },
                LoadingRoute::ClosedBodyGeometryPolicy,
                &[ka * medium.sound_speed],
                medium,
                source(),
                4,
                options,
                || false,
            )
        };
    let bake = |surface: &SpherePanels, ka: f64, options: Option<&IterativeOptions>| {
        bake_modes(surface, 2, ka, options)
    };
    let bits = |b: &ModalLoadingBake| -> Vec<(u64, u64)> {
        b.rows()[0]
            .impedance()
            .iter()
            .map(|z| (z.re.to_bits(), z.im.to_bits()))
            .collect()
    };
    let options = routed_options(0.05);
    // ka 0.05 is plain CBIE and ka 3 Burton-Miller under the policy.
    for ka in [0.05, 3.0] {
        let dense = bake(&below, ka, None).unwrap();
        let routed = bake(&below, ka, Some(&options)).unwrap();
        assert_eq!(bits(&routed), bits(&dense));
        assert_eq!(routed.declaration_hash(), dense.declaration_hash());
        let row = &routed.rows()[0];
        assert_eq!(row.reciprocity_bound(), DENSE_RECIPROCITY_BOUND);
        assert_eq!(row.iterative_error_scale(), None);
    }
    // Below the cap the declared panel rule reaches the dense bake. A
    // Burton-Miller row (ka 3) is solve_radiation_batch_with_rule's, bitwise,
    // under its own declaration; a plain-CBIE row (ka 0.05) ignores the rule.
    for rule in [
        helmholtz::PanelRule::TriangleNear { max_pairs: 1 << 20 },
        helmholtz::PanelRule::Triangle { max_pairs: 1 << 20 },
    ] {
        let ruled = IterativeOptions {
            panel_rule: rule,
            ..options.clone()
        };
        let centroid = bake_modes(&below, 1, 3.0, None).unwrap();
        let bm = bake_modes(&below, 1, 3.0, Some(&ruled)).unwrap();
        let field: Vec<C64> = vec![C64::ONE; below.areas().len()];
        let k = 3.0 * medium.sound_speed / medium.sound_speed;
        let solved = helmholtz::solve_radiation_batch_with_rule(
            &below,
            k,
            medium,
            &[&field],
            Formulation::BurtonMiller,
            rule,
        )
        .unwrap();
        let mut force = C64::ZERO;
        for (&p, &area) in solved[0].pressure.iter().zip(below.areas()) {
            force = force + p.scale(area);
        }
        assert_eq!(bits(&bm), vec![(force.re.to_bits(), force.im.to_bits())]);
        assert_ne!(bits(&bm), bits(&centroid));
        assert_ne!(bm.declaration_hash(), centroid.declaration_hash());
        let cbie = bake_modes(&below, 1, 0.05, Some(&ruled)).unwrap();
        assert_eq!(
            bits(&cbie),
            bits(&bake_modes(&below, 1, 0.05, None).unwrap())
        );
    }
    let above = SpherePanels::icosphere(1.0, 5).unwrap();
    assert!(above.areas().len() > helmholtz::MAX_DENSE_PANELS);
    assert!(matches!(
        bake(&above, 0.05, None),
        Err(ModalLoadingError::SolverRoute {
            route: "M2.5-W5c",
            ..
        })
    ));
    // Past the cap the closed outward-surface admission is still live.
    let mut open = above.triangles().unwrap().to_vec();
    open.pop();
    let open = SpherePanels::from_triangles(open).unwrap();
    assert!(matches!(
        bake(&open, 0.05, Some(&options)),
        Err(ModalLoadingError::Helmholtz(_))
    ));
    let routed = bake_modes(&above, 1, 0.05, Some(&options)).unwrap();
    let diagnostics = routed.rows()[0].bem_diagnostics().unwrap();
    assert_eq!(diagnostics.formulation, Formulation::PlainCbie);
    assert_eq!(diagnostics.maximum_condition_lower_bound, 1.0);
    assert!(diagnostics.minimum_panels_per_wavelength >= 6.0);
    let serial = IterativeOptions {
        workers: 3,
        m2l_cache_bytes: 0,
        ..options.clone()
    };
    let repeat = bake_modes(&above, 1, 0.05, Some(&serial)).unwrap();
    assert_eq!(bits(&repeat), bits(&routed));
    assert_eq!(repeat.declaration_hash(), routed.declaration_hash());
    let looser = IterativeOptions {
        tol: 1e-7,
        ..options.clone()
    };
    assert_ne!(
        bake_modes(&above, 1, 0.05, Some(&looser))
            .unwrap()
            .declaration_hash(),
        routed.declaration_hash()
    );
    // A Burton-Miller row (ka 3) binds its panel rule: the triangle near rule
    // changes the numbers, so it is a different declaration. The smallest
    // lat-long sphere past the cap (8280 panels) and a loose GMRES tolerance
    // keep this cheap: the rule's effect is far above the solver error.
    let just_above = resolved_sphere(46);
    assert!(just_above.areas().len() > helmholtz::MAX_DENSE_PANELS);
    let centroid = IterativeOptions {
        tol: 1e-6,
        fmm_order: 4,
        ..options.clone()
    };
    let rule = IterativeOptions {
        panel_rule: helmholtz::PanelRule::TriangleNear { max_pairs: 1 << 24 },
        ..centroid.clone()
    };
    let centroid_bm = bake_modes(&just_above, 1, 3.0, Some(&centroid)).unwrap();
    let rule_bm = bake_modes(&just_above, 1, 3.0, Some(&rule)).unwrap();
    assert_ne!(rule_bm.declaration_hash(), centroid_bm.declaration_hash());
    assert_ne!(bits(&rule_bm), bits(&centroid_bm));
    let pair = bake(&above, 0.05, Some(&options)).unwrap();
    let row = &pair.rows()[0];
    let scale = row.iterative_error_scale().unwrap();
    // The exact bits, for comparison with the same bakes before the twins.
    let hex = |b: &ModalLoadingBake| {
        bits(b)
            .iter()
            .flat_map(|&(re, im)| [re, im])
            .map(|v| format!("{v:016x}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"routed-modal-loading-bits\",\"panels\":{},\"one_mode_z_bits\":\"{}\",\"two_mode_z_bits\":\"{}\"}}",
        above.areas().len(),
        hex(&routed),
        hex(&pair)
    );
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"routed-modal-loading-reciprocity\",\"inputs\":\"analytic icosphere-5, ka 0.05, plain CBIE\",\"panels\":{},\"gmres_tol\":{},\"fmm_orders\":{:?},\"eps_gmres\":{},\"eps_fmm\":{},\"eps_fmm_next\":{},\"observed\":{},\"bound\":{},\"verdict\":\"pass\",\"twins\":[{},{}]}}",
        above.areas().len(),
        options.tol,
        scale.fmm_orders,
        scale.gmres_residual,
        scale.fmm_error,
        scale.fmm_error_next,
        row.reciprocity_relative_defect(),
        row.reciprocity_bound(),
        scale.twins[0].to_json(),
        scale.twins[1].to_json()
    );
    assert_eq!(scale.fmm_orders, [6, 8, 10]);
    for (twin, order) in scale.twins.iter().zip([8, 10]) {
        assert_eq!(twin.fmm_order, order);
        assert_eq!(twin.max_translation_work, options.max_translation_work);
        assert_eq!(twin.max_coefficient_slots, options.max_coefficient_slots);
        assert_eq!(twin.applies, 2);
    }
    assert!(scale.gmres_residual <= options.tol && scale.estimate_is_resolved());
    assert_eq!(
        row.reciprocity_bound(),
        (GMRES_RECIPROCITY_FACTOR * scale.gmres_residual).max(scale.fmm_error)
    );
    assert!(row.reciprocity_relative_defect() <= row.reciprocity_bound());
    let z00 = row.impedance()[0];
    let single = routed.rows()[0].impedance()[0];
    assert_eq!(
        (z00.re.to_bits(), z00.im.to_bits()),
        (single.re.to_bits(), single.im.to_bits())
    );
}

/// G1: the low-ka translation (dipole) mode, which motivated the
/// geometry-policy route. Unit icosphere at 320, 1280 and 5120 panels,
/// ka = 0.05, below the sphere's plain-CBIE limit of 0.8 * pi * sqrt(3) / 2.
/// - The policy route equals a direct plain-CBIE batch solve bitwise, and
///   the Burton-Miller route bitwise above the limit.
/// - It converges at observed order >= 1.8 on both modes, and its finest
///   dipole resistance is within 1% of the closed form.
/// - The Burton-Miller route's errors are recorded beside it.
#[test]
fn geometry_policy_route_resolves_the_low_ka_dipole() {
    let medium = air();
    let area = 4.0 * PI;
    let ka = 0.05_f64;
    let (x2, x4) = (ka * ka, ka.powi(4));
    let oracle = [
        C64::new(x2, -ka).scale(area * medium.density * medium.sound_speed / (1.0 + x2)),
        C64::new(x4, -ka * (2.0 + x2))
            .scale(area * medium.density * medium.sound_speed / (3.0 * (4.0 + x4))),
    ];
    let shapes_of = |surface: &SpherePanels| {
        vec![
            vec![1.0; surface.areas().len()],
            surface.normals().iter().map(|n| n[2]).collect::<Vec<f64>>(),
        ]
    };
    let load = |surface: &SpherePanels, route, k: f64| {
        let shapes = shapes_of(surface);
        modal_radiation_loading(
            SurfaceModeTraces::ClosedBody {
                surface,
                shapes: &shapes,
            },
            route,
            &[k * medium.sound_speed],
            medium,
            source(),
            4,
            || false,
        )
        .unwrap()
        .rows()[0]
            .impedance()
            .to_vec()
    };
    let policy = LoadingRoute::ClosedBodyGeometryPolicy;
    let bm = LoadingRoute::ClosedBodyBem(Formulation::BurtonMiller);

    let coarse = SpherePanels::icosphere(1.0, 2).unwrap();
    let shapes = shapes_of(&coarse);
    let fields: Vec<Vec<C64>> = shapes
        .iter()
        .map(|s| s.iter().map(|&v| C64::from_re(v)).collect())
        .collect();
    let refs: Vec<_> = fields.iter().map(Vec::as_slice).collect();
    let solved =
        helmholtz::solve_radiation_batch(&coarse, ka, medium, &refs, Formulation::PlainCbie)
            .unwrap();
    let areas = coarse.areas();
    let direct: Vec<C64> = shapes
        .iter()
        .flat_map(|row| {
            solved.iter().map(move |solution| {
                row.iter()
                    .zip(&solution.pressure)
                    .zip(areas)
                    .fold(C64::ZERO, |s, ((&phi, &p), &a)| s + p.scale(phi * a))
            })
        })
        .collect();
    let routed = load(&coarse, policy, ka);
    assert!(
        routed
            .iter()
            .zip(&direct)
            .all(|(a, b)| a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits()),
        "below the limit the policy route must be the plain-CBIE solve"
    );
    let above = load(&coarse, policy, 3.0);
    assert!(
        above
            .iter()
            .zip(&load(&coarse, bm, 3.0))
            .all(|(a, b)| a.re.to_bits() == b.re.to_bits() && a.im.to_bits() == b.im.to_bits()),
        "above the limit the policy route must be the Burton-Miller solve"
    );

    let mut errors = [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]];
    let mut resistance = [0.0_f64; 2];
    for rung in 2..=4 {
        let surface = SpherePanels::icosphere(1.0, rung).unwrap();
        for (r, route) in [policy, bm].into_iter().enumerate() {
            let z = load(&surface, route, ka);
            for mode in 0..2 {
                let zm = z[mode * 2 + mode];
                let error = (zm - oracle[mode]).abs() / oracle[mode].abs();
                let resistance_error = (zm.re - oracle[mode].re).abs() / oracle[mode].re;
                errors[r][mode].push(error);
                if mode == 1 {
                    resistance[r] = resistance_error;
                }
                println!(
                    "{{\"suite\":\"fs-bem\",\"case\":\"policy-route-low-ka\",\"route\":\"{}\",\"panels\":{},\"ka\":{ka},\"mode\":{mode},\"err\":{error:e},\"resistance_err\":{resistance_error:e}}}",
                    if r == 0 {
                        "geometry-policy"
                    } else {
                        "burton-miller"
                    },
                    surface.areas().len()
                );
            }
        }
    }
    // Panel size halves per icosphere rung, so the observed order is log2
    // of the error ratio over the finest pair.
    let order = |e: &Vec<f64>| (e[1] / e[2]).log2();
    for mode in 0..2 {
        assert!(
            errors[0][mode].windows(2).all(|w| w[1] < w[0]) && order(&errors[0][mode]) >= 1.8,
            "policy route mode {mode} must converge at order >= 1.8: {:?}",
            errors[0][mode]
        );
    }
    assert!(
        resistance[0] <= 1e-2,
        "finest dipole resistance error {}",
        resistance[0]
    );
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"policy-route-low-ka-orders\",\"geometry_policy\":[{:.3},{:.3}],\"burton_miller\":[{:.3},{:.3}],\"dipole_resistance_err\":[{:e},{:e}]}}",
        order(&errors[0][0]),
        order(&errors[0][1]),
        order(&errors[1][0]),
        order(&errors[1][1]),
        resistance[0],
        resistance[1]
    );
}
