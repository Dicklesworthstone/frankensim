//! Tests use original gas laws and independent spherical-wave continuity.
//! Synthetic relaxation reservoirs are numerical fixtures, not measured air.
use super::super::tests::{derivative, exp_i, outgoing, plane, regular, state};
use super::*;
use crate::helmholtz::bc::{ForceDipole, IncidentSource, Monopole, PlaneWave};
use crate::helmholtz::transmission::multiple as native;

fn sphere(level: u32, center: f64) -> SpherePanels {
    let original = SpherePanels::icosphere(1.0, level).unwrap();
    SpherePanels::from_triangles(
        original
            .triangles()
            .unwrap()
            .iter()
            .map(|t| t.map(|p| [p[0], p[1], p[2] + center]))
            .collect(),
    )
    .unwrap()
}
fn relative(a: &[C64], b: &[C64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let numerator: f64 = a.iter().zip(b).map(|(&a, &b)| (a - b).abs().powi(2)).sum();
    let denominator: f64 = b.iter().map(|b| b.abs().powi(2)).sum();
    (numerator / denominator.max(1e-60)).sqrt()
}
fn sample_error(a: IncidentSample, b: IncidentSample) -> f64 {
    ((a.pressure_pa - b.pressure_pa).abs() / b.pressure_pa.abs())
        .max(relative(&a.velocity_m_s, &b.velocity_m_s))
}
fn dot(p: &[C64], v: &[C64], areas: &[f64]) -> C64 {
    p.iter()
        .zip(v)
        .zip(areas)
        .fold(C64::ZERO, |sum, ((&p, &v), &a)| {
            sum + (p * v.conj()).scale(0.5 * a)
        })
}

#[test]
fn empty_reservoirs_recover_native_multiple_and_one_body_recovers_original_molecular_owner() {
    let gate = CancelGate::new();
    let outside = state(293.15, 101325.0, None);
    let a = state(353.15, 101325.0, None);
    let b = state(253.15, 101325.0, None);
    let surfaces = [sphere(0, -1.6), sphere(0, 1.6)];
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
    let molecular =
        MolecularFluidBodies::new(&bodies, &outside, 200.0, 0.01, Options::default(), &gate)
            .unwrap();
    let original = native::FluidBodies::new(
        &[
            native::Body {
                surface: &surfaces[0],
                material: a.gas(),
            },
            native::Body {
                surface: &surfaces[1],
                material: b.gas(),
            },
        ],
        outside.gas(),
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let actual = molecular.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let expected = original.scatter_field(plane(), 0.01, &gate).unwrap();
    for (actual, expected) in [
        (
            actual.total_pressure_pa(),
            expected.solution().total_pressure_pa(),
        ),
        (
            actual.total_normal_velocity_m_s(),
            expected.solution().total_normal_velocity_m_s(),
        ),
        (
            actual.outgoing_pressure_pa(),
            expected.solution().outgoing_pressure_pa(),
        ),
        (
            actual.outgoing_normal_velocity_m_s(),
            expected.solution().outgoing_normal_velocity_m_s(),
        ),
        (
            actual.scattered_pressure_pa(),
            expected.solution().scattered_pressure_pa(),
        ),
        (
            actual.scattered_normal_velocity_m_s(),
            expected.solution().scattered_normal_velocity_m_s(),
        ),
    ] {
        assert!(relative(actual, expected) < 2e-11);
    }
    assert!(std::ptr::eq(molecular.exterior_source(), &outside));
    assert!(std::ptr::eq(molecular.body(1).unwrap().material, &b));

    let outside = state(293.15, 101325.0, Some(0.1));
    let inside = state(353.15, 101325.0, Some(5.0));
    let body = [Body {
        surface: &surfaces[0],
        material: &inside,
    }];
    let multiple =
        MolecularFluidBodies::new(&body, &outside, 200.0, 0.01, Options::default(), &gate).unwrap();
    let single = MolecularFluidInterface::new(
        &surfaces[0],
        &outside,
        &inside,
        200.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let actual = multiple.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let expected = single.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    assert!(relative(actual.total_pressure_pa(), expected.total_pressure_pa()) < 2e-11);
    assert!(
        relative(
            actual.total_normal_velocity_m_s(),
            expected.total_normal_velocity_m_s()
        ) < 2e-11
    );
    for (point, interior) in [([0.0, 0.0, -1.3], true), ([0.2, 0.0, 2.0], false)] {
        let a = if interior {
            actual.interior_pressure_velocity_at(0, &[point], 0.01, 1000, &gate)
        } else {
            actual.total_pressure_velocity_at(&[point], 0.01, 1000, &gate)
        }
        .unwrap()[0];
        let b = if interior {
            expected.interior_pressure_velocity_at(&[point], 0.01, 1000, &gate)
        } else {
            expected.total_pressure_velocity_at(&[point], 0.01, 1000, &gate)
        }
        .unwrap()[0];
        assert!(sample_error(a, b) < 2e-11);
    }
}

#[test]
fn transparent_body_does_not_change_neighbors_and_still_transmits_their_full_field() {
    let gate = CancelGate::new();
    let outside = state(293.15, 101325.0, Some(0.1));
    let inside = state(353.15, 101325.0, Some(5.0));
    let surfaces = [sphere(0, -1.6), sphere(0, 1.6)];
    let owner = MolecularFluidBodies::new(
        &[
            Body {
                surface: &surfaces[0],
                material: &outside,
            },
            Body {
                surface: &surfaces[1],
                material: &inside,
            },
        ],
        &outside,
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let single = MolecularFluidInterface::new(
        &surfaces[1],
        &outside,
        &inside,
        200.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let all = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let active = single.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    assert!(
        all.outgoing_pressure_pa()[..20]
            .iter()
            .chain(&all.outgoing_normal_velocity_m_s()[..20])
            .all(|&v| v == C64::ZERO)
    );
    assert!(relative(&all.total_pressure_pa()[20..], active.total_pressure_pa()) < 2e-12);
    assert!(
        relative(
            &all.total_normal_velocity_m_s()[20..],
            active.total_normal_velocity_m_s()
        ) < 2e-12
    );
    let point = [0.0, 0.0, -1.3];
    let actual = all
        .interior_pressure_velocity_at(0, &[point], 0.01, 1000, &gate)
        .unwrap()[0];
    let expected = active
        .total_pressure_velocity_at(&[point], 0.01, 1000, &gate)
        .unwrap()[0];
    assert!(sample_error(actual, expected) < 2e-11);
    let incident = plane()
        .sample_wave(
            owner.wave,
            outside.gas().density * 200.0,
            None,
            point,
            &mut || false,
        )
        .unwrap();
    assert!((actual.pressure_pa - incident.pressure_pa).abs() > 1e-5);
    assert!(
        all.interior_pressure_velocity_at(0, &[point], 0.01, 39, &gate)
            .is_err()
    );

    let transparent = MolecularFluidBodies::new(
        &[
            Body {
                surface: &surfaces[0],
                material: &outside,
            },
            Body {
                surface: &surfaces[1],
                material: &outside,
            },
        ],
        &outside,
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.5, -0.2),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [2.5, 0.0, 0.0],
            volume_velocity_m3_s: C64::new(0.002, 0.001),
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [0.0, 2.5, 0.0],
            force_on_fluid_n: [C64::new(0.2, -0.1), C64::ZERO, C64::ONE],
        }),
    ];
    let incident = IncidentField::Sources(&sources);
    let field = transparent
        .scatter_field(incident, 0.01, 1000, &gate)
        .unwrap();
    assert!(
        field
            .scattered_pressure_pa()
            .iter()
            .chain(field.scattered_normal_velocity_m_s())
            .all(|&v| v == C64::ZERO)
    );
    for (point, body) in [
        ([0.1, 0.0, -1.3], Some(0)),
        ([0.1, 0.0, 1.4], Some(1)),
        ([0.0, 0.0, 3.0], None),
    ] {
        let expected = incident
            .sample_wave(
                transparent.wave,
                outside.gas().density * 200.0,
                None,
                point,
                &mut || false,
            )
            .unwrap();
        let actual = if let Some(body) = body {
            field.interior_pressure_velocity_at(body, &[point], 0.01, 1000, &gate)
        } else {
            field.total_pressure_velocity_at(&[point], 0.01, 1000, &gate)
        }
        .unwrap()[0];
        assert_eq!(actual.pressure_pa, expected.pressure_pa);
        assert_eq!(actual.velocity_m_s, expected.velocity_m_s);
    }
}

#[test]
fn coherent_molecular_work_keeps_cross_terms_and_reuses_the_coupled_factor() {
    let gate = CancelGate::new();
    let outside = state(293.15, 101325.0, Some(0.1));
    let a = state(220.0, 101325.0, Some(5.0));
    let b = state(480.0, 101325.0, Some(2.0));
    let surfaces = [sphere(1, -1.5), sphere(1, 1.5)];
    let owner = MolecularFluidBodies::new(
        &[
            Body {
                surface: &surfaces[0],
                material: &a,
            },
            Body {
                surface: &surfaces[1],
                material: &b,
            },
        ],
        &outside,
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let finite_source = IncidentSource::Monopole(Monopole {
        position_m: [2.5, 0.0, 0.0],
        volume_velocity_m3_s: C64::new(0.02, -0.03),
    });
    let p = IncidentSource::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    });
    let sources = [p, finite_source];
    let finite_sources = [finite_source];
    let first = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let second = owner
        .scatter_field(IncidentField::Sources(&finite_sources), 0.01, 1000, &gate)
        .unwrap();
    let sum = owner
        .scatter_field(IncidentField::Sources(&sources), 0.01, 1000, &gate)
        .unwrap();
    let p: Vec<_> = first
        .scattered_pressure_pa()
        .iter()
        .zip(second.scattered_pressure_pa())
        .map(|(&a, &b)| a + b)
        .collect();
    let v: Vec<_> = first
        .scattered_normal_velocity_m_s()
        .iter()
        .zip(second.scattered_normal_velocity_m_s())
        .map(|(&a, &b)| a + b)
        .collect();
    assert!(relative(sum.scattered_pressure_pa(), &p) < 2e-12);
    assert!(relative(sum.scattered_normal_velocity_m_s(), &v) < 2e-12);
    let work = dot(&p, &v, owner.surface().areas());
    let cross = dot(
        first.scattered_pressure_pa(),
        second.scattered_normal_velocity_m_s(),
        owner.surface().areas(),
    ) + dot(
        second.scattered_pressure_pa(),
        first.scattered_normal_velocity_m_s(),
        owner.surface().areas(),
    );
    assert!((work - sum.scattered_work().surface_power).abs() < 2e-12 * work.abs());
    assert!(
        (work
            - first.scattered_work().surface_power
            - second.scattered_work().surface_power
            - cross)
            .abs()
            < 2e-12 * work.abs()
    );
    assert!(
        cross.abs() > 1e-4 * work.abs(),
        "fixture must expose coherent source interference"
    );
    let local = dot(
        sum.outgoing_pressure_pa(),
        sum.outgoing_normal_velocity_m_s(),
        owner.surface().areas(),
    );
    assert!(
        (work - local).abs() > 1e-3 * work.abs(),
        "local outgoing densities are not complete physical work"
    );
    let (lo, hi) = sum.scattered_work_roundoff_interval();
    assert!(lo <= work.re && work.re <= hi);
    // Reordering body blocks changes neither fields nor the work for a body.
    let reversed = MolecularFluidBodies::new(
        &[
            Body {
                surface: &surfaces[1],
                material: &b,
            },
            Body {
                surface: &surfaces[0],
                material: &a,
            },
        ],
        &outside,
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let reversed = reversed
        .scatter_field(IncidentField::Sources(&sources), 0.01, 1000, &gate)
        .unwrap();
    let point = [[0.3, 0.2, 3.5]];
    assert!(
        sample_error(
            sum.total_pressure_velocity_at(&point, 0.01, 1000, &gate)
                .unwrap()[0],
            reversed
                .total_pressure_velocity_at(&point, 0.01, 1000, &gate)
                .unwrap()[0]
        ) < 2e-11
    );
    for body in 0..2 {
        let a = sum.interface_work()[body].0.surface_power;
        let b = reversed.interface_work()[1 - body].0.surface_power;
        assert!((a - b).abs() < 2e-11 * a.abs());
    }
}

// Continuum reference: project complex outgoing spherical waves onto each
// regular interior expansion, then solve the small mode system independently
// of all panel kernels, shared transmission matrices and the fs-la solver.
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

struct SpherePair {
    centers: [f64; 2],
    ko: C64,
    ki: [C64; 2],
    rho_o: f64,
    rho_i: [f64; 2],
    omega: f64,
    outgoing: Vec<C64>,
    transmitted: Vec<C64>,
    order: usize,
}
impl SpherePair {
    fn new(
        outside: &RelaxationState,
        inside: [&RelaxationState; 2],
        omega: f64,
        order: usize,
        coupled: bool,
    ) -> Self {
        let centers = [-1.5, 1.5];
        let ko = outside.relaxing_wavenumber(omega).unwrap();
        let ki = inside.map(|s| s.relaxing_wavenumber(omega).unwrap());
        let rho_o = outside.gas().density;
        let rho_i = inside.map(|s| s.gas().density);
        let n = 2 * order;
        let mut translation = vec![C64::ZERO; n * n];
        let mut scattering = vec![C64::ZERO; n];
        let mut transmission = vec![C64::ZERO; n];
        let mut incoming = vec![C64::ZERO; n];
        let quadrature = gauss(48);
        let (jo, ho) = (regular(ko, order + 1), outgoing(ko, order + 1));
        for body in 0..2 {
            let ji = regular(ki[body], order + 1);
            let eta = rho_i[body] / rho_o;
            let mut phase = exp_i(ko.scale(centers[body]));
            for l in 0..order {
                let index = body * order + l;
                let (djo, dho, dji) = (
                    derivative(&jo, ko, l),
                    derivative(&ho, ko, l),
                    derivative(&ji, ki[body], l),
                );
                let denominator = ko * ji[l] * dho - (ki[body] * dji * ho[l]).scale(1.0 / eta);
                scattering[index] =
                    ((ki[body] * dji * jo[l]).scale(1.0 / eta) - ko * ji[l] * djo) / denominator;
                transmission[index] = (C64::new(0.0, 1.0) / ko) / denominator;
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
                let h = outgoing(ko.scale(r), order + 1);
                let ps = legendre(dz / r, order);
                let pt = legendre(mu, order);
                for l in 0..order {
                    for m in 0..order {
                        let index = (body * order + l) * n + source * order + m;
                        translation[index] = translation[index]
                            + (h[m] / jo[l]).scale(ps[m] * pt[l] * w * (2 * l + 1) as f64 / 2.0);
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
        let transmitted = (0..n)
            .map(|i| {
                let mut exciting = incoming[i];
                for j in 0..n {
                    exciting = exciting + translation[i * n + j] * outgoing[j];
                }
                transmission[i] * exciting
            })
            .collect();
        Self {
            centers,
            ko,
            ki,
            rho_o,
            rho_i,
            omega,
            outgoing,
            transmitted,
            order,
        }
    }

    fn axis(&self, z: f64, interior: Option<usize>) -> IncidentSample {
        let (mut p, mut gradient) = (C64::ZERO, C64::ZERO);
        if interior.is_none() {
            p = exp_i(self.ko.scale(z));
            gradient = C64::new(0.0, 1.0) * self.ko * p;
        }
        for body in 0..2 {
            if interior.is_some_and(|b| b != body) {
                continue;
            }
            let dz = z - self.centers[body];
            let sign = dz.signum();
            let k = if interior.is_some() {
                self.ki[body]
            } else {
                self.ko
            };
            let argument = k.scale(dz.abs());
            let f = if interior.is_some() {
                regular(argument, self.order + 1)
            } else {
                outgoing(argument, self.order + 1)
            };
            let coefficients = if interior.is_some() {
                &self.transmitted
            } else {
                &self.outgoing
            };
            let mut parity = 1.0;
            for l in 0..self.order {
                let c = coefficients[body * self.order + l];
                p = p + (c * f[l]).scale(parity);
                gradient = gradient + (c * k * derivative(&f, argument, l)).scale(parity * sign);
                parity *= sign;
            }
        }
        let rho = interior.map_or(self.rho_o, |b| self.rho_i[b]);
        IncidentSample {
            pressure_pa: p,
            velocity_m_s: [
                C64::ZERO,
                C64::ZERO,
                gradient * C64::new(0.0, -1.0 / (rho * self.omega)),
            ],
        }
    }

    fn work(&self, body: usize) -> C64 {
        let k = self.ki[body];
        let j = regular(k, self.order + 1);
        let mut work = C64::ZERO;
        for l in 0..self.order {
            let coefficient = self.transmitted[body * self.order + l];
            let p = coefficient * j[l];
            let v = coefficient
                * k
                * derivative(&j, k, l)
                * C64::new(0.0, -1.0 / (self.rho_i[body] * self.omega));
            work = work + (p * v.conj()).scale(2.0 * core::f64::consts::PI / (2 * l + 1) as f64);
        }
        work
    }

    fn volume_heat(&self, body: usize, inside: &RelaxationState) -> f64 {
        let mut radial = 0.0;
        for i in 0..=256 {
            let r = i as f64 / 256.0;
            let j = regular(self.ki[body].scale(r), self.order);
            let norm: f64 = j
                .iter()
                .enumerate()
                .map(|(l, &j)| {
                    (self.transmitted[body * self.order + l] * j).abs().powi(2) / (2 * l + 1) as f64
                })
                .sum();
            let weight = if i == 0 || i == 256 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            radial += weight * r * r * norm;
        }
        let compressibility = inside
            .caloric_response(self.omega)
            .unwrap()
            .bulk_modulus_pa
            .recip();
        2.0 * core::f64::consts::PI * self.omega * compressibility.im * radial / (3.0 * 256.0)
    }
}

#[test]
fn relaxing_sphere_pair_converges_to_independent_all_orders_fields_and_heat() {
    let gate = CancelGate::new();
    let outside = state(293.15, 101325.0, Some(0.1));
    let interiors = [
        state(220.0, 101325.0, Some(5.0)),
        state(480.0, 101325.0, Some(2.0)),
    ];
    let omega = 200.0;
    let exact = SpherePair::new(&outside, [&interiors[0], &interiors[1]], omega, 9, true);
    let refined = SpherePair::new(&outside, [&interiors[0], &interiors[1]], omega, 11, true);
    let isolated = SpherePair::new(&outside, [&interiors[0], &interiors[1]], omega, 9, false);
    assert!(
        relative(&exact.outgoing, &isolated.outgoing) > 0.005,
        "fixture must detect mutual reflections"
    );
    for (z, body) in [(-1.15, Some(0)), (1.85, Some(1)), (-3.5, None), (3.5, None)] {
        assert!(sample_error(exact.axis(z, body), refined.axis(z, body)) < 1e-7);
    }
    for body in 0..2 {
        let work = exact.work(body);
        let heat = exact.volume_heat(body, &interiors[body]);
        assert!(work.re < 0.0 && heat > 0.0);
        assert!(
            (heat + work.re).abs() < 1e-8 * heat,
            "continuum interface work equals original molecular bulk heat"
        );
    }
    let (mut last_field, mut last_work) = (f64::INFINITY, f64::INFINITY);
    for level in 0..=2 {
        let surfaces = [sphere(level, -1.5), sphere(level, 1.5)];
        let owner = MolecularFluidBodies::new(
            &[
                Body {
                    surface: &surfaces[0],
                    material: &interiors[0],
                },
                Body {
                    surface: &surfaces[1],
                    material: &interiors[1],
                },
            ],
            &outside,
            omega,
            0.01,
            Options::default(),
            &gate,
        )
        .unwrap();
        let field = owner.scatter_field(plane(), 0.01, 10_000, &gate).unwrap();
        let mut field_error = 0.0_f64;
        for (z, body) in [(-1.15, Some(0)), (1.85, Some(1)), (-3.5, None), (3.5, None)] {
            let point = [[0.0, 0.0, z]];
            let observed = if let Some(body) = body {
                field.interior_pressure_velocity_at(body, &point, 0.01, 10_000, &gate)
            } else {
                field.total_pressure_velocity_at(&point, 0.01, 10_000, &gate)
            }
            .unwrap()[0];
            field_error = field_error.max(sample_error(observed, exact.axis(z, body)));
        }
        let mut work_error = 0.0_f64;
        for body in 0..2 {
            let expected = exact.work(body);
            let actual = field.interface_work()[body].0.surface_power;
            work_error = work_error.max((actual - expected).abs() / expected.abs());
        }
        eprintln!(
            "molecular multiple spheres level={level}: field_error={field_error:e}, work_error={work_error:e}"
        );
        assert!(field_error < last_field && work_error < last_work);
        last_field = field_error;
        last_work = work_error;
    }
    // Inherited flat P0/centroid discretization: a decreasing independent
    // three-mesh ladder, not a general continuum error certificate.
    assert!(last_field < 0.08 && last_work < 0.15);
}

#[test]
fn molecular_multiple_geometry_work_budgets_and_cancellation_are_enforced() {
    let gate = CancelGate::new();
    let outside = state(293.15, 101325.0, Some(0.1));
    let inside = state(353.15, 101325.0, Some(5.0));
    let mismatch = state(353.15, 102325.0, Some(5.0));
    let surfaces = [sphere(0, -1.5), sphere(0, 1.5), sphere(0, -1.0)];
    let body = |index, material| Body {
        surface: &surfaces[index],
        material,
    };
    let bodies = [body(0, &inside), body(1, &inside)];
    for options in [
        Options {
            maximum_bytes: 0,
            ..Options::default()
        },
        Options {
            maximum_panels: 39,
            ..Options::default()
        },
    ] {
        assert!(MolecularFluidBodies::new(&bodies, &outside, 200.0, 0.01, options, &gate).is_err());
    }
    assert!(
        MolecularFluidBodies::new(
            &[body(0, &inside), body(2, &inside)],
            &outside,
            200.0,
            0.01,
            Options::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        MolecularFluidBodies::new(
            &[body(0, &inside), body(1, &mismatch)],
            &outside,
            200.0,
            0.01,
            Options::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        MolecularFluidBodies::new(&bodies, &outside, 20_000.0, 0.01, Options::default(), &gate)
            .is_err()
    );
    let owner =
        MolecularFluidBodies::new(&bodies, &outside, 200.0, 0.01, Options::default(), &gate)
            .unwrap();
    assert!(owner.scatter_field(plane(), 0.01, 39, &gate).is_err());
    let sources = [Monopole {
        position_m: [0.0, 0.0, 1.5],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        owner
            .scatter_field(IncidentField::Monopoles(&sources), 0.01, 1000, &gate)
            .is_err()
    );
    let field = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    assert!(
        field
            .total_pressure_velocity_at(&[[0.0, 0.0, -1.5]], 0.01, 1000, &gate)
            .is_err()
    );
    assert!(
        field
            .interior_pressure_velocity_at(0, &[[0.0, 0.0, 1.5]], 0.01, 1000, &gate)
            .is_err()
    );
    assert!(
        field
            .interior_pressure_velocity_at(2, &[[0.0, 0.0, 1.5]], 0.01, 1000, &gate)
            .is_err()
    );
    assert!(
        field
            .total_pressure_velocity_at(&[[0.0, 0.0, 3.5]], 0.01, 39, &gate)
            .is_err()
    );
    gate.request();
    assert!(matches!(
        MolecularFluidBodies::new(&bodies, &outside, 200.0, 0.01, Options::default(), &gate),
        Err(TransmissionError::Helmholtz(HelmholtzError::Cancelled))
    ));
    assert!(matches!(
        owner.scatter_field(plane(), 0.01, 1000, &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        field.interior_pressure_velocity_at(0, &[[0.0, 0.0, -1.3]], 0.01, 1000, &gate),
        Err(HelmholtzError::Cancelled)
    ));
}
