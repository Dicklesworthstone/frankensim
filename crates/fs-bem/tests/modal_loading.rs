//! M2.5-W4a: complete modal loading. Strict sphere tolerances are qualification
//! gates, not a claim that the existing collocation owner already meets them.
use fs_bem::helmholtz::rayleigh_modal::{self, AreaSample, PreparedSurface};
use fs_bem::helmholtz::{self, Formulation, Medium};
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

/// G2, sphere radius 1 m. Derivation: p_l=A*h_l^(1)(kr)*P_l and
/// v_l=k*p_l'/(i*w*rho) give z_l=i*rho*c*h_l/h_l'. For l=0 this
/// reduces to (-i*x)/(1-i*x). For translation, integral cos²(theta)dS=S/3;
/// h_1=-exp(ix)*(x+i)/x² gives S*rho*c/3 * [x^4-i*(2*x+x³)]/(4+x^4).
/// Verified source: Sutherland WR68-12, NASA 19690015869, Appendix B,
/// equations B-9 (printed p43), B-18 (p45), and Fig.8 (p26).
#[test]
fn sphere_modal_loading_closed_forms() {
    let started = Instant::now();
    let medium = air();
    let area = 4.0 * PI;
    for ka in [0.05_f64, 0.5, 1.0, PI, 5.0, 10.0] {
        let x2 = ka * ka;
        let x4 = x2 * x2;
        let oracle = [
            C64::new(x2, -ka).scale(area * medium.density * medium.sound_speed / (1.0 + x2)),
            C64::new(x4, -ka * (2.0 + x2))
                .scale(area * medium.density * medium.sound_speed / (3.0 * (4.0 + x4))),
        ];
        let mut finest = None;
        let mut admitted = 0;
        for (rung, latitude) in [24, 32, 40].into_iter().enumerate() {
            let surface = resolved_sphere(latitude);
            let shapes = vec![
                vec![1.0; surface.areas().len()],
                surface.normals().iter().map(|n| n[2]).collect(),
            ];
            let result = modal_radiation_loading(
                SurfaceModeTraces::ClosedBody {
                    surface: &surface,
                    shapes: &shapes,
                },
                LoadingRoute::ClosedBodyBem(Formulation::BurtonMiller),
                &[ka * medium.sound_speed],
                medium,
                source(),
                4,
                || false,
            );
            match result {
                Ok(bake) => {
                    admitted += 1;
                    let row = &bake.rows()[0];
                    let z = row.impedance();
                    let mut worst = 0.0_f64;
                    for mode in 0..2 {
                        let error = (z[mode * 2 + mode] - oracle[mode]).abs() / oracle[mode].abs();
                        worst = worst.max(error);
                        println!(
                            "{{\"suite\":\"fs-bem\",\"case\":\"modal-loading-sphere\",\"inputs\":\"analytic sphere a=1m; NASA WR68-12 B-9/B-18\",\"axis\":\"panels\",\"rung\":{rung},\"ka\":{ka},\"panels\":{},\"mode\":{mode},\"observed\":[{},{}],\"oracle\":[{},{}],\"err\":{error},\"min_eig_re\":{},\"tol\":1e-3,\"verdict\":\"{}\",\"elapsed_s\":{}}}",
                            surface.areas().len(),
                            z[mode * 2 + mode].re,
                            z[mode * 2 + mode].im,
                            oracle[mode].re,
                            oracle[mode].im,
                            row.minimum_resistance_eigenvalue(),
                            if error <= 1e-3 { "pass" } else { "fail" },
                            started.elapsed().as_secs_f64()
                        );
                    }
                    if rung == 2 {
                        finest = Some(worst);
                    }
                }
                Err(error) => println!(
                    "{{\"suite\":\"fs-bem\",\"case\":\"modal-loading-sphere\",\"inputs\":\"analytic sphere a=1m\",\"axis\":\"panels\",\"rung\":{rung},\"ka\":{ka},\"panels\":{},\"err\":\"{}\",\"tol\":1e-3,\"verdict\":\"refused\",\"elapsed_s\":{}}}",
                    surface.areas().len(),
                    error,
                    started.elapsed().as_secs_f64()
                ),
            }
        }
        assert!(
            admitted == 3 && finest.is_some_and(|e| e <= 1e-3),
            "strict sphere qualification unmet at ka={ka}: {finest:?}; refine the numerical owner, never relax 1e-3"
        );
    }
}
