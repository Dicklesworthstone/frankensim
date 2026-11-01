//! G1/G2/G3 multiple-fluid reflection against independent spherical waves.
//! Numerical material states are explicit fixtures, not measured materials.
use fs_bem::helmholtz::bc::{
    ForceDipole, IncidentField, IncidentSource, IncidentTraces, Monopole, PlaneWave,
};
use fs_bem::helmholtz::transmission::multiple::{Body, FluidBodies};
use fs_bem::helmholtz::transmission::{FluidInterface, Options, TransmissionError};
use fs_bem::helmholtz::{HelmholtzError, Medium, with_dense_workers};
use fs_bem::{near_field, panel3d::SpherePanels};
use fs_exec::CancelGate;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::c64::C64;
use std::sync::atomic::{AtomicUsize, Ordering};

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn air() -> Medium {
    Medium {
        density: 1.0,
        sound_speed: 1.0,
    }
}
fn plane() -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    })
}
fn exp(x: f64) -> C64 {
    C64::new(x.cos(), x.sin())
}
fn sphere(level: u32, center: [f64; 3]) -> SpherePanels {
    let original = SpherePanels::icosphere(1.0, level).unwrap();
    SpherePanels::from_triangles(
        original
            .triangles()
            .unwrap()
            .iter()
            .map(|t| t.map(|p| std::array::from_fn(|i| p[i] + center[i])))
            .collect(),
    )
    .unwrap()
}
fn norm(values: &[C64]) -> f64 {
    values
        .iter()
        .map(|z| z.re * z.re + z.im * z.im)
        .sum::<f64>()
        .sqrt()
}
fn error(actual: &[C64], expected: &[C64]) -> f64 {
    norm(
        &actual
            .iter()
            .zip(expected)
            .map(|(&a, &b)| a - b)
            .collect::<Vec<_>>(),
    ) / norm(expected)
}
fn zero_bits(values: &[C64]) {
    for z in values {
        assert_eq!((z.re.to_bits(), z.im.to_bits()), (0, 0));
    }
}
struct Counted {
    state: FluidAcousticState,
    queries: AtomicUsize,
}
impl FluidAcoustics for Counted {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.queries.fetch_add(1, Ordering::Relaxed);
        Ok(self.state)
    }
}

#[test]
fn transparent_bodies_preserve_coherent_volume_force_and_plane_sources_exactly() {
    let s0 = sphere(0, [0.0, 0.0, -1.6]);
    let s1 = sphere(0, [0.0, 0.0, 1.6]);
    let exterior = Counted {
        state: fluid(1.0, 1.0),
        queries: AtomicUsize::new(0),
    };
    let interiors = [
        Counted {
            state: exterior.state,
            queries: AtomicUsize::new(0),
        },
        Counted {
            state: exterior.state,
            queries: AtomicUsize::new(0),
        },
    ];
    let gate = CancelGate::new();
    let bodies = [
        Body {
            surface: &s0,
            material: &interiors[0],
        },
        Body {
            surface: &s1,
            material: &interiors[1],
        },
    ];
    let owner = FluidBodies::new(&bodies, &exterior, 0.6, 0.01, Options::default(), &gate).unwrap();
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.4, -0.2),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [0.0, 0.0, 4.0],
            volume_velocity_m3_s: C64::new(0.5, 0.3),
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [3.0, 0.0, 0.0],
            force_on_fluid_n: [C64::ONE, C64::ZERO, C64::new(0.2, -0.6)],
        }),
    ];
    let incident =
        IncidentTraces::from_sources(owner.surface(), 0.6, air(), &sources, 0.01).unwrap();
    for _ in 0..2 {
        let answer = owner.scatter(&incident, &gate).unwrap();
        zero_bits(answer.scattered_pressure_pa());
        zero_bits(answer.scattered_normal_velocity_m_s());
        assert_eq!(answer.total_pressure_pa(), incident.pressure_pa());
        assert_eq!(
            answer.total_normal_velocity_m_s(),
            incident.normal_velocity_m_s()
        );
    }
    let field = IncidentField::Sources(&sources);
    let solution = owner.scatter_field(field, 0.01, &gate).unwrap();
    let points = [[2.5, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let observed = solution
        .total_pressure_at(&points, 0.01, near_field::Options::default())
        .unwrap();
    for (&point, actual) in points.iter().zip(observed.pressure) {
        assert_eq!(
            actual,
            field
                .pressure_velocity_at(0.6, air(), point)
                .unwrap()
                .pressure_pa
        );
    }
    assert_eq!(exterior.queries.load(Ordering::Relaxed), 1);
    assert!(
        interiors
            .iter()
            .all(|x| x.queries.load(Ordering::Relaxed) == 1)
    );
    assert_eq!(owner.body_count(), 2);
    assert_eq!(owner.interior_state(1), Some(interiors[1].state));
    assert!(owner.body(2).is_none());
}

#[test]
fn one_body_agrees_with_original_and_a_transparent_neighbor_has_zero_reaction() {
    let s0 = sphere(1, [0.0, 0.0, 0.0]);
    let s1 = sphere(1, [0.0, 0.0, 3.0]);
    let (outside, inside) = (fluid(1.0, 1.0), fluid(1.8, 0.8));
    let gate = CancelGate::new();
    let original =
        FluidInterface::new(&s0, &outside, &inside, 0.6, Options::default(), &gate).unwrap();
    let bodies = [Body {
        surface: &s0,
        material: &inside,
    }];
    let one = FluidBodies::new(&bodies, &outside, 0.6, 0.01, Options::default(), &gate).unwrap();
    let mut pair = bodies.to_vec();
    pair.push(Body {
        surface: &s1,
        material: &outside,
    });
    let two = FluidBodies::new(&pair, &outside, 0.6, 0.01, Options::default(), &gate).unwrap();
    let a = original.scatter_field(plane(), 0.01, &gate).unwrap();
    let b = one.scatter_field(plane(), 0.01, &gate).unwrap();
    let c = two.scatter_field(plane(), 0.01, &gate).unwrap();
    let n = s0.areas().len();
    assert!(
        error(
            b.solution().outgoing_pressure_pa(),
            &a.solution().scattered().pressure
        ) < 1e-13
    );
    assert!(
        error(
            &c.solution().outgoing_pressure_pa()[..n],
            &a.solution().scattered().pressure
        ) < 1e-13
    );
    zero_bits(&c.solution().outgoing_pressure_pa()[n..]);
    zero_bits(&c.solution().outgoing_normal_velocity_m_s()[n..]);
    assert!(norm(&c.solution().scattered_pressure_pa()[n..]) > 0.01);
    // The transparent body's interior still carries the other body's field.
    let point = [0.0, 0.0, 3.0];
    let exterior = a
        .total_pressure_at(&[point], 0.01, near_field::Options::default())
        .unwrap()
        .pressure[0];
    let interior = c
        .solution()
        .interior_pressure_at(1, &[point], 0.01, near_field::Options::default())
        .unwrap()
        .pressure[0];
    assert!((interior - exterior).abs() / exterior.abs() < 0.015);
    assert!((exterior - exp(0.6 * 3.0)).abs() > 0.005);
}

// Independent spherical-wave reference. No panel kernels, original matrix
// assembly or fs-la solve are used. Translation coefficients are projected on
// a regular interior sphere by Gauss-Legendre integration. The retained-wave
// order is checked independently before comparison with the BEM mesh ladder.
fn bessel(x: f64, count: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut j = vec![0.0; count + 1];
    for (n, value) in j.iter_mut().enumerate() {
        let mut term = (1..=n).fold(1.0, |p, l| p * x / (2 * l + 1) as f64);
        let mut sum = term;
        for m in 1..80 {
            term *= -x * x / (2 * m * (2 * n + 2 * m + 1)) as f64;
            sum += term;
            if term.abs() < 1e-18 * sum.abs() {
                break;
            }
        }
        *value = sum;
    }
    let mut y = vec![0.0; count + 1];
    y[0] = -x.cos() / x;
    y[1] = -x.cos() / (x * x) - x.sin() / x;
    for n in 1..count {
        y[n + 1] = (2 * n + 1) as f64 * y[n] / x - y[n - 1];
    }
    let dj = (0..count).map(|n| n as f64 * j[n] / x - j[n + 1]).collect();
    let dy = (0..count).map(|n| n as f64 * y[n] / x - y[n + 1]).collect();
    (j, y, dj, dy)
}
fn legendre(x: f64, n: usize) -> Vec<f64> {
    let mut p = vec![1.0; n];
    if n > 1 {
        p[1] = x;
    }
    for l in 1..n - 1 {
        p[l + 1] = ((2 * l + 1) as f64 * x * p[l] - l as f64 * p[l - 1]) / (l + 1) as f64;
    }
    p
}
fn gauss(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            for _ in 0..20 {
                let p = legendre(x, n + 1);
                let dp = n as f64 * (x * p[n] - p[n - 1]) / (x * x - 1.0);
                let step = p[n] / dp;
                x -= step;
                if step.abs() < 2e-16 {
                    break;
                }
            }
            let p = legendre(x, n + 1);
            let dp = n as f64 * (x * p[n] - p[n - 1]) / (x * x - 1.0);
            (x, 2.0 / ((1.0 - x * x) * dp * dp))
        })
        .collect()
}
fn solve(mut a: Vec<C64>, mut b: Vec<C64>) -> Vec<C64> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        for j in 0..n {
            a.swap(k * n + j, p * n + j);
        }
        b.swap(k, p);
        for i in k + 1..n {
            let f = a[i * n + k] / a[k * n + k];
            for j in k + 1..n {
                a[i * n + j] = a[i * n + j] - f * a[k * n + j];
            }
            b[i] = b[i] - f * b[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            b[i] = b[i] - a[i * n + j] * b[j];
        }
        b[i] = b[i] / a[i * n + i];
    }
    b
}
fn sphere_pair(
    k: f64,
    centers: [f64; 2],
    fluids: [(f64, f64); 2],
    order: usize,
    coupled: bool,
) -> (Vec<C64>, Vec<C64>) {
    let n = 2 * order;
    let mut translation = vec![C64::ZERO; n * n];
    let mut scattering = vec![C64::ZERO; n];
    let mut transmission = vec![C64::ZERO; n];
    let mut incoming = vec![C64::ZERO; n];
    let quadrature = gauss(48);
    let (jo, yo, djo, dyo) = bessel(k, order);
    for body in 0..2 {
        let (rho, c) = fluids[body];
        let ki = k / c;
        let (ji, _, dji, _) = bessel(ki, order);
        let mut phase = exp(k * centers[body]);
        for l in 0..order {
            let index = body * order + l;
            let denominator = C64::new(
                k * ji[l] * djo[l] - ki * dji[l] * jo[l] / rho,
                k * ji[l] * dyo[l] - ki * dji[l] * yo[l] / rho,
            );
            scattering[index] =
                C64::from_re(ki * dji[l] * jo[l] / rho - k * ji[l] * djo[l]) / denominator;
            transmission[index] = C64::new(0.0, 1.0 / k) / denominator;
            incoming[index] = phase.scale((2 * l + 1) as f64);
            phase = phase * C64::new(0.0, 1.0);
        }
        if !coupled {
            continue;
        }
        let source = 1 - body;
        for &(mu, w) in &quadrature {
            let dz = centers[body] + mu - centers[source];
            let r = (1.0 - mu * mu + dz * dz).sqrt();
            let (j, y, _, _) = bessel(k * r, order);
            let ps = legendre(dz / r, order);
            let pt = legendre(mu, order);
            for l in 0..order {
                for m in 0..order {
                    let index = (body * order + l) * n + source * order + m;
                    translation[index] = translation[index]
                        + C64::new(j[m], y[m])
                            .scale(ps[m] * pt[l] * w * (2 * l + 1) as f64 / (2.0 * jo[l]));
                }
            }
        }
    }
    let mut matrix = vec![C64::ZERO; n * n];
    let mut drive = vec![C64::ZERO; n];
    for i in 0..n {
        for j in 0..n {
            matrix[i * n + j] = -scattering[i] * translation[i * n + j];
        }
        matrix[i * n + i] = C64::ONE;
        drive[i] = scattering[i] * incoming[i];
    }
    let outgoing = solve(matrix, drive);
    let far = [-1.0, 0.0, 1.0]
        .into_iter()
        .map(|mu| {
            let p = legendre(mu, order);
            let mut sum = C64::ZERO;
            for body in 0..2 {
                let mut phase = C64::new(0.0, -1.0) * exp(-k * mu * centers[body]);
                for l in 0..order {
                    sum = sum + phase * outgoing[body * order + l].scale(p[l] / k);
                    phase = phase * C64::new(0.0, -1.0);
                }
            }
            sum
        })
        .collect();
    let interior = (0..2)
        .map(|body| {
            let index = body * order;
            let mut excited = incoming[index];
            for j in 0..n {
                excited = excited + translation[index * n + j] * outgoing[j];
            }
            transmission[index] * excited
        })
        .collect();
    (far, interior)
}

#[test]
fn distinct_interacting_spheres_converge_to_independent_multiple_partial_waves() {
    let (k, centers, media) = (0.6, [-1.5, 1.5], [(1.6, 0.75), (0.7, 1.3)]);
    let expected = sphere_pair(k, centers, media, 7, true);
    let refined = sphere_pair(k, centers, media, 9, true);
    assert!(error(&expected.0, &refined.0) < 1e-7);
    assert!(error(&expected.1, &refined.1) < 1e-7);
    let independent = sphere_pair(k, centers, media, 7, false);
    assert!(
        error(&independent.0, &expected.0) > 0.01,
        "fixture must expose mutual reflection"
    );
    let exterior = fluid(1.0, 1.0);
    let interiors = [fluid(media[0].0, media[0].1), fluid(media[1].0, media[1].1)];
    let gate = CancelGate::new();
    let directions = [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let (mut last_far, mut last_interior) = (f64::INFINITY, f64::INFINITY);
    for level in 0..=2 {
        let surfaces = [
            sphere(level, [0.0, 0.0, centers[0]]),
            sphere(level, [0.0, 0.0, centers[1]]),
        ];
        let bodies = [
            Body {
                surface: &surfaces[0],
                material: &interiors[0],
            },
            Body {
                surface: &surfaces[1],
                material: &interiors[1],
            },
        ];
        let owner =
            FluidBodies::new(&bodies, &exterior, k, 0.01, Options::default(), &gate).unwrap();
        let field = owner.scatter_field(plane(), 0.01, &gate).unwrap();
        let far = field
            .solution()
            .far_field(&directions, 1_000_000, &gate)
            .unwrap();
        let inside: Vec<_> = (0..2)
            .map(|i| {
                field
                    .solution()
                    .interior_pressure_at(
                        i,
                        &[[0.0, 0.0, centers[i]]],
                        0.01,
                        near_field::Options::default(),
                    )
                    .unwrap()
                    .pressure[0]
            })
            .collect();
        let ef = error(&far, &expected.0);
        let ei = error(&inside, &expected.1);
        eprintln!("multiple spheres level={level}: far_error={ef:e}, interior_error={ei:e}");
        assert!(ef < 0.65 * last_far && ei < 0.65 * last_interior);
        last_far = ef;
        last_interior = ei;
    }
    // The inherited flat P0/disc owner has a 5% far-field conformance target
    // on this 320-panel-per-body rung; the full decreasing ladder is required.
    // This fixture does not provide a general spatial-error certificate.
    assert!(last_far < 0.05 && last_interior < 0.02);
}

#[test]
fn mutual_scattered_flux_retains_interference_and_coherent_drive_scaling() {
    let outside = fluid(1.0, 1.0);
    let a = fluid(1.6, 0.75);
    let b = fluid(0.7, 1.3);
    let gate = CancelGate::new();
    let surfaces = [sphere(1, [0.0, 0.0, -1.5]), sphere(1, [0.0, 0.0, 1.5])];
    let bodies = [
        Body {
            surface: &surfaces[0],
            material: &a,
        },
        Body {
            surface: &surfaces[1],
            material: &b,
        },
    ];
    let owner = FluidBodies::new(&bodies, &outside, 0.6, 0.01, Options::default(), &gate).unwrap();
    let original = owner.scatter_field(plane(), 0.01, &gate).unwrap();
    let scaled = owner
        .scatter_field(
            IncidentField::PlaneWave(PlaneWave {
                direction: [0.0, 0.0, 1.0],
                amplitude_pa: C64::new(0.0, 2.0),
            }),
            0.01,
            &gate,
        )
        .unwrap();
    let s = original.solution();
    let dot = |p: &[C64], v: &[C64]| {
        p.iter()
            .zip(v)
            .zip(owner.surface().areas())
            .fold(C64::ZERO, |sum, ((&p, &v), &a)| {
                sum + (p * v.conj()).scale(0.5 * a)
            })
    };
    let full = dot(s.scattered_pressure_pa(), s.scattered_normal_velocity_m_s());
    let isolated = dot(s.outgoing_pressure_pa(), s.outgoing_normal_velocity_m_s());
    assert!((full - s.scattered_work().surface_power).abs() < 1e-13 * full.abs());
    assert!(
        (full - isolated).abs() > 0.01 * full.abs(),
        "fixture must expose interference work"
    );
    assert!(
        (scaled.solution().scattered_work().surface_power - full.scale(4.0)).abs()
            < 1e-12 * full.abs()
    );
    let (lo, hi) = s.scattered_work_roundoff_interval();
    assert!(lo <= full.re && full.re <= hi);
    assert!(
        s.interface_work()
            .iter()
            .all(|(work, _)| work.surface_power.re.is_finite())
    );
}

#[test]
fn weak_contrast_array_resolves_forward_born_volume_without_incident_cancellation() {
    let exterior = fluid(1.0, 1.0);
    let gate = CancelGate::new();
    let k = 0.6;
    let surfaces = [sphere(1, [0.0, 0.0, -1.5]), sphere(1, [0.0, 0.0, 1.5])];
    let volume: f64 = surfaces
        .iter()
        .flat_map(|s| s.triangles().unwrap())
        .map(|t| {
            let [a, b, c] = *t;
            (a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum();
    for delta in [1e-6, 1e-10, 1e-12] {
        let interior = fluid(1.0, 1.0 / (1.0 + delta));
        let bodies = [
            Body {
                surface: &surfaces[0],
                material: &interior,
            },
            Body {
                surface: &surfaces[1],
                material: &interior,
            },
        ];
        let owner =
            FluidBodies::new(&bodies, &exterior, k, 0.01, Options::default(), &gate).unwrap();
        let field = owner.scatter_field(plane(), 0.01, &gate).unwrap();
        let actual = field
            .solution()
            .far_field(&[[0.0, 0.0, 1.0]], 1_000_000, &gate)
            .unwrap()[0];
        let ki = k / interior.sound_speed_m_s();
        let expected = C64::from_re((ki - k) * (ki + k) * volume / (4.0 * core::f64::consts::PI));
        assert!(
            (actual - expected).abs() / expected.abs() < 0.04,
            "delta={delta} actual={actual:?} expected={expected:?}"
        );
    }
}

#[test]
fn body_order_and_translation_preserve_the_complete_coherent_solution() {
    let outside = fluid(1.0, 1.0);
    let a = fluid(2.0, 0.8);
    let b = fluid(0.7, 1.2);
    let gate = CancelGate::new();
    let centers = [[0.0, 0.0, -1.5], [0.0, 0.0, 1.5]];
    let surfaces = centers.map(|p| sphere(0, p));
    let bodies = [
        Body {
            surface: &surfaces[0],
            material: &a,
        },
        Body {
            surface: &surfaces[1],
            material: &b,
        },
    ];
    let owner = FluidBodies::new(&bodies, &outside, 0.7, 0.01, Options::default(), &gate).unwrap();
    let reversed = FluidBodies::new(
        &[bodies[1], bodies[0]],
        &outside,
        0.7,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let d = [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let reference = owner.scatter_field(plane(), 0.01, &gate).unwrap();
    let other = reversed.scatter_field(plane(), 0.01, &gate).unwrap();
    let f = reference
        .solution()
        .far_field(&d, 1_000_000, &gate)
        .unwrap();
    let g = other.solution().far_field(&d, 1_000_000, &gate).unwrap();
    assert!(error(&f, &g) < 2e-12);
    let shift = [2.3, -0.8, 1.7];
    let translated = centers.map(|p| sphere(0, std::array::from_fn(|i| p[i] + shift[i])));
    let bodies = [
        Body {
            surface: &translated[0],
            material: &a,
        },
        Body {
            surface: &translated[1],
            material: &b,
        },
    ];
    let owner = FluidBodies::new(&bodies, &outside, 0.7, 0.01, Options::default(), &gate).unwrap();
    let shifted = owner.scatter_field(plane(), 0.01, &gate).unwrap();
    let actual = shifted.solution().far_field(&d, 1_000_000, &gate).unwrap();
    let expected: Vec<_> = d
        .iter()
        .zip(f)
        .map(|(d, f)| f * exp(0.7 * (shift[2] - (0..3).map(|i| d[i] * shift[i]).sum::<f64>())))
        .collect();
    assert!(error(&actual, &expected) < 2e-12);
}

#[test]
fn geometry_source_receiver_memory_and_cancellation_refusals_are_enforced() {
    let outside = fluid(1.0, 1.0);
    let inside = fluid(1.3, 0.8);
    let gate = CancelGate::new();
    let a = sphere(0, [0.0; 3]);
    let b = sphere(0, [0.0, 0.0, 3.0]);
    let overlap = sphere(0, [0.0, 0.0, 1.0]);
    let bodies = [
        Body {
            surface: &a,
            material: &inside,
        },
        Body {
            surface: &b,
            material: &inside,
        },
    ];
    assert!(FluidBodies::new(&[], &outside, 0.6, 0.01, Options::default(), &gate).is_err());
    let underflow = fluid(1e-200, 1e100);
    assert!(
        FluidBodies::new(
            &[Body {
                surface: &a,
                material: &underflow
            }],
            &underflow,
            1e-200,
            0.01,
            Options::default(),
            &gate,
        )
        .is_err()
    );
    assert!(
        FluidBodies::new(
            &[
                bodies[0],
                Body {
                    surface: &overlap,
                    material: &inside
                }
            ],
            &outside,
            0.6,
            0.01,
            Options::default(),
            &gate
        )
        .is_err()
    );
    assert!(FluidBodies::new(&bodies, &outside, 0.6, 2.0, Options::default(), &gate).is_err());
    assert!(
        FluidBodies::new(
            &bodies,
            &outside,
            0.6,
            0.01,
            Options {
                maximum_panels: 39,
                ..Options::default()
            },
            &gate
        )
        .is_err()
    );
    assert!(matches!(
        with_dense_workers(64, || FluidBodies::new(
            &bodies,
            &outside,
            0.6,
            0.01,
            Options {
                maximum_bytes: 1,
                ..Options::default()
            },
            &gate
        )),
        Err(TransmissionError::Helmholtz(
            HelmholtzError::CacheBudget { .. }
        ))
    ));
    let owner = FluidBodies::new(&bodies, &outside, 0.6, 0.01, Options::default(), &gate).unwrap();
    let singular = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        owner
            .scatter_field(IncidentField::Monopoles(&singular), 0.01, &gate)
            .is_err()
    );
    let field = owner.scatter_field(plane(), 0.01, &gate).unwrap();
    assert!(
        field
            .solution()
            .far_field(&[[0.0, 0.0, 1.0]], 1, &gate)
            .is_err()
    );
    assert!(
        field
            .solution()
            .far_field(&[[0.0, 0.0, 2.0]], 1_000_000, &gate)
            .is_err()
    );
    assert!(
        field
            .total_pressure_at(&[[0.0; 3]], 0.01, near_field::Options::default())
            .is_err()
    );
    assert!(
        field
            .solution()
            .interior_pressure_at(0, &[[0.0, 0.0, 3.0]], 0.01, near_field::Options::default())
            .is_err()
    );
    assert!(
        field
            .solution()
            .interior_pressure_at(2, &[[0.0; 3]], 0.01, near_field::Options::default())
            .is_err()
    );
    let wrong = IncidentTraces::from_field(owner.surface(), 0.7, air(), plane(), 0.01).unwrap();
    assert!(owner.scatter(&wrong, &gate).is_err());
    gate.request();
    assert!(matches!(
        field
            .solution()
            .far_field(&[[0.0, 0.0, 1.0]], 1_000_000, &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        owner.scatter_field(plane(), 0.01, &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        FluidBodies::new(&bodies, &outside, 0.6, 0.01, Options::default(), &gate),
        Err(TransmissionError::Helmholtz(HelmholtzError::Cancelled))
    ));
}
