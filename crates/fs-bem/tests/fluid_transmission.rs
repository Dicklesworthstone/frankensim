//! Penetrable-fluid scattering against independent spherical partial waves.
//! Numerical fluid states are explicit conformance fixtures, not measured data.
use fs_bem::helmholtz::bc::{IncidentField, IncidentSource, IncidentTraces, Monopole, PlaneWave};
use fs_bem::helmholtz::transmission::{FluidInterface, Options, TransmissionError};
use fs_bem::helmholtz::{HelmholtzError, Medium, far_field};
use fs_bem::near_field;
use fs_bem::panel3d::SpherePanels;
use fs_exec::CancelGate;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::c64::C64;
use std::sync::atomic::{AtomicUsize, Ordering};

fn fluid(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}
fn plane() -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    })
}
fn cexp(x: f64) -> C64 {
    C64::new(x.cos(), x.sin())
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-30)
}

// Independent partial-wave oracle for a unit sphere. j_n uses its entire
// power series, y_n its three-term recurrence; no BEM kernels or transfer
// matrices enter these continuity-derived coefficients.
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
fn coefficients(ko: f64, ki: f64, eta: f64) -> Vec<(C64, C64)> {
    let (jo, yo, djo, dyo) = bessel(ko, 18);
    let (ji, _, dji, _) = bessel(ki, 18);
    (0..18)
        .map(|n| {
            let denominator = C64::new(
                ko * ji[n] * djo[n] - ki * dji[n] * jo[n] / eta,
                ko * ji[n] * dyo[n] - ki * dji[n] * yo[n] / eta,
            );
            let numerator = ki * dji[n] * jo[n] / eta - ko * ji[n] * djo[n];
            // j_n h'_n-j'_n h_n=i/(k_o R)^2, with R=1.
            (
                C64::from_re(numerator) / denominator,
                C64::new(0.0, 1.0 / ko) / denominator,
            )
        })
        .collect()
}
fn sphere_far(ko: f64, ki: f64, eta: f64, cosine: f64) -> C64 {
    let mut p0 = 1.0;
    let mut p1 = cosine;
    let mut sum = C64::ZERO;
    for (n, (a, _)) in coefficients(ko, ki, eta).into_iter().enumerate() {
        let p = if n == 0 {
            1.0
        } else if n == 1 {
            cosine
        } else {
            let p = ((2 * n - 1) as f64 * cosine * p1 - (n - 1) as f64 * p0) / n as f64;
            p0 = p1;
            p1 = p;
            p
        };
        sum = sum + a * C64::new(0.0, -((2 * n + 1) as f64) * p / ko);
    }
    sum
}
fn sphere_exterior_axis(ko: f64, ki: f64, eta: f64, radius: f64) -> C64 {
    let (j, y, _, _) = bessel(ko * radius, 18);
    let mut phase = C64::ONE;
    let mut scattered = C64::ZERO;
    for (n, (a, _)) in coefficients(ko, ki, eta).into_iter().enumerate() {
        scattered = scattered + phase * a * C64::new(j[n], y[n]).scale((2 * n + 1) as f64);
        phase = phase * C64::new(0.0, 1.0);
    }
    cexp(ko * radius) + scattered
}

struct CountedFluid {
    state: FluidAcousticState,
    queries: AtomicUsize,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.queries.fetch_add(1, Ordering::Relaxed);
        Ok(self.state)
    }
}

#[test]
fn identical_fluids_preserve_coherent_traces_exactly_and_query_sources_once() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let exterior = CountedFluid {
        state: fluid(1.2, 340.0),
        queries: AtomicUsize::new(0),
    };
    let interior = CountedFluid {
        state: exterior.state,
        queries: AtomicUsize::new(0),
    };
    let gate = CancelGate::new();
    let interface = FluidInterface::new(
        &surface,
        &exterior,
        &interior,
        340.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.4, -0.7),
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [0.0, 0.0, 3.0],
            volume_velocity_m3_s: C64::new(0.002, -0.001),
        }),
    ];
    let medium = Medium {
        density: 1.2,
        sound_speed: 340.0,
    };
    let traces = IncidentTraces::from_sources(&surface, 1.0, medium, &sources, 0.01).unwrap();
    for _ in 0..2 {
        let solved = interface.scatter(&traces, &gate).unwrap();
        for z in solved
            .scattered()
            .pressure
            .iter()
            .chain(&solved.scattered().velocity)
        {
            assert_eq!(
                (z.re.to_bits(), z.im.to_bits()),
                (0.0f64.to_bits(), 0.0f64.to_bits())
            );
        }
        assert_eq!(solved.total_pressure_pa(), traces.pressure_pa());
        assert_eq!(
            solved.total_normal_velocity_m_s(),
            traces.normal_velocity_m_s()
        );
    }
    let solved = interface
        .scatter_field(IncidentField::Sources(&sources), 0.01, &gate)
        .unwrap();
    let point = [0.0, 0.0, 2.0];
    let observed = solved
        .total_pressure_at(&[point], 0.01, near_field::Options::default())
        .unwrap();
    let exact = IncidentField::Sources(&sources)
        .pressure_velocity_at(1.0, medium, point)
        .unwrap()
        .pressure_pa;
    assert_eq!(observed.pressure[0], exact);
    assert_eq!(exterior.queries.load(Ordering::Relaxed), 1);
    assert_eq!(interior.queries.load(Ordering::Relaxed), 1);
    assert_eq!(interface.angular_frequency_rad_s(), 340.0);
    assert_eq!(interface.interior_state(), interior.state);
    assert_eq!(interface.exterior_state(), exterior.state);
}

#[test]
fn density_and_speed_contrasts_converge_to_fluid_sphere_on_three_meshes() {
    let gate = CancelGate::new();
    let directions = [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    // All three meshes meet BOTH existing six-panel/wavelength admissions.
    for (rho_i, c_i, omega) in [(2.0, 1.0, 1.0), (1.0, 2.0 / 3.0, 1.0), (0.5, 0.5, 0.75)] {
        let exterior = fluid(1.0, 1.0);
        let interior = fluid(rho_i, c_i);
        let expected: Vec<_> = [-1.0, 0.0, 1.0]
            .into_iter()
            .map(|c| sphere_far(omega, omega / c_i, rho_i, c))
            .collect();
        let scale = expected.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
        let mut previous_error = f64::INFINITY;
        let mut previous_work = f64::INFINITY;
        for level in 0..=2 {
            let surface = SpherePanels::icosphere(1.0, level).unwrap();
            let interface = FluidInterface::new(
                &surface,
                &exterior,
                &interior,
                omega,
                Options::default(),
                &gate,
            )
            .unwrap();
            let field = interface.scatter_field(plane(), 0.01, &gate).unwrap();
            let solved = field.solution();
            let actual = far_field(
                &surface,
                solved.scattered(),
                solved.scattered().medium,
                &directions,
            );
            let error = actual
                .iter()
                .zip(&expected)
                .map(|(&a, &b)| (a - b).abs() / scale)
                .fold(0.0_f64, f64::max);
            let work =
                solved.interface_work.surface_power.re.abs() / solved.interface_work.apparent_power;
            assert!(
                error < previous_error,
                "rho={rho_i}, c={c_i}, mesh={level}, far error={error}"
            );
            assert!(
                work < previous_work,
                "rho={rho_i}, c={c_i}, mesh={level}, work defect={work}"
            );
            previous_error = error;
            previous_work = work;
            if level == 2 {
                assert!(
                    error < 0.055,
                    "resolved constant-panel sphere error {error}"
                );
                assert!(work < 0.01, "lossless interface work defect {work}");
            }
        }
    }
}

#[test]
fn actual_interior_and_total_exterior_receivers_converge_and_refuse_cross_wiring() {
    let gate = CancelGate::new();
    let exterior = fluid(1.0, 1.0);
    let interior = fluid(1.0, 2.0 / 3.0);
    let expected_inside = coefficients(1.0, 1.5, 1.0)[0].1;
    let expected_outside = sphere_exterior_axis(1.0, 1.5, 1.0, 2.0);
    let mut previous = (f64::INFINITY, f64::INFINITY);
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(1.0, level).unwrap();
        let interface = FluidInterface::new(
            &surface,
            &exterior,
            &interior,
            1.0,
            Options::default(),
            &gate,
        )
        .unwrap();
        let field = interface.scatter_field(plane(), 0.01, &gate).unwrap();
        let inside = field
            .solution()
            .interior_pressure_at(&[[0.0; 3]], 0.01, near_field::Options::default())
            .unwrap();
        let outside = field
            .total_pressure_at(&[[0.0, 0.0, 2.0]], 0.01, near_field::Options::default())
            .unwrap();
        let error = (
            relative(inside.pressure[0], expected_inside),
            relative(outside.pressure[0], expected_outside),
        );
        assert!(
            error.0 < previous.0 && error.1 < previous.1,
            "mesh={level}, error={error:?}"
        );
        previous = error;
        if level == 2 {
            assert!(error.0 < 0.04 && error.1 < 0.015, "{error:?}");
        }
        assert!(
            field
                .solution()
                .interior_pressure_at(&[[0.0, 0.0, 2.0]], 0.01, near_field::Options::default())
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
                .interior_pressure_at(
                    &[surface.centroids()[0]],
                    0.01,
                    near_field::Options::default()
                )
                .is_err()
        );
    }
}

#[test]
fn weak_speed_contrast_retains_the_forward_born_limit() {
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let gate = CancelGate::new();
    let exterior = fluid(1.0, 1.0);
    for delta in [1e-6, 1e-10, 1e-12] {
        let interior = fluid(1.0, 1.0 / (1.0 + delta));
        let interface = FluidInterface::new(
            &surface,
            &exterior,
            &interior,
            1.0,
            Options::default(),
            &gate,
        )
        .unwrap();
        let solved = interface.scatter_field(plane(), 0.01, &gate).unwrap();
        let far = far_field(
            &surface,
            solved.solution().scattered(),
            Medium {
                density: 1.0,
                sound_speed: 1.0,
            },
            &[[0.0, 0.0, 1.0]],
        )[0];
        let dk = interface.interior_wavenumber_rad_m() - 1.0;
        // First Born forward amplitude for unit volume radius:
        // (k_i²-k_o²)V/(4pi)=2*dk/3+O(dk²), independent of panel kernels.
        assert!(
            relative(far, C64::from_re(2.0 * dk / 3.0)) < 0.04,
            "delta={delta}, far={far:?}"
        );
    }
}

#[test]
fn cancellation_dense_caps_frequency_binding_and_regular_source_admission_are_enforced() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let exterior = fluid(1.0, 1.0);
    let interior = fluid(2.0, 1.0);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        FluidInterface::new(
            &surface,
            &exterior,
            &interior,
            1.0,
            Options::default(),
            &gate
        ),
        Err(TransmissionError::Helmholtz(HelmholtzError::Cancelled))
    ));
    let gate = CancelGate::new();
    let small = Options {
        maximum_panels: 20,
        ..Options::default()
    };
    assert!(matches!(
        FluidInterface::new(&surface, &exterior, &interior, 1.0, small, &gate),
        Err(TransmissionError::Helmholtz(HelmholtzError::WorkCap { .. }))
    ));
    let small = Options {
        maximum_bytes: 1,
        ..Options::default()
    };
    assert!(matches!(
        FluidInterface::new(&surface, &exterior, &interior, 1.0, small, &gate),
        Err(TransmissionError::Helmholtz(
            HelmholtzError::CacheBudget { .. }
        ))
    ));
    let fast = fluid(1.0, 0.01);
    assert!(matches!(
        FluidInterface::new(&surface, &exterior, &fast, 1.0, Options::default(), &gate),
        Err(TransmissionError::Helmholtz(
            HelmholtzError::TooCoarse { .. }
        ))
    ));
    let interface = FluidInterface::new(
        &surface,
        &exterior,
        &interior,
        1.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let sources = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        interface
            .scatter_field(IncidentField::Monopoles(&sources), 0.01, &gate)
            .is_err()
    );
    let source = [Monopole {
        position_m: surface.centroids()[0],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        interface
            .scatter_field(IncidentField::Monopoles(&source), 0.01, &gate)
            .is_err()
    );
    let medium = Medium {
        density: 1.0,
        sound_speed: 1.0,
    };
    let wrong = IncidentTraces::from_field(&surface, 0.5, medium, plane(), 0.01).unwrap();
    assert!(interface.scatter(&wrong, &gate).is_err());
    let correct = IncidentTraces::from_field(&surface, 1.0, medium, plane(), 0.01).unwrap();
    let mut disconnected = Vec::new();
    for shift in [-3.0, 3.0] {
        disconnected.extend(surface.triangles().unwrap().iter().map(|triangle| {
            triangle.map(|mut point| {
                point[0] += shift;
                point
            })
        }));
    }
    let disconnected = SpherePanels::from_triangles(disconnected).unwrap();
    assert!(
        FluidInterface::new(
            &disconnected,
            &exterior,
            &interior,
            1.0,
            Options::default(),
            &gate
        )
        .is_err()
    );
    gate.request();
    assert!(matches!(
        interface.scatter(&correct, &gate),
        Err(HelmholtzError::Cancelled)
    ));
}

#[test]
fn transmitting_body_retains_finite_source_phase_and_coherent_interior_response() {
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let exterior = CountedFluid {
        state: fluid(1.0, 1.0),
        queries: AtomicUsize::new(0),
    };
    let interior = CountedFluid {
        state: fluid(1.0, 2.0 / 3.0),
        queries: AtomicUsize::new(0),
    };
    let gate = CancelGate::new();
    let interface = FluidInterface::new(
        &surface,
        &exterior,
        &interior,
        1.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let monopole = IncidentSource::Monopole(Monopole {
        position_m: [0.4, -0.2, 3.0],
        volume_velocity_m3_s: C64::new(1.2, 0.3),
    });
    let mixed = [
        monopole,
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.02, -0.01),
        }),
    ];
    let single = [monopole];
    let monopole_transfer = coefficients(1.0, 1.5, 1.0)[0].1;
    for sources in [&single[..], &mixed[..]] {
        let incident = IncidentField::Sources(sources);
        let field = interface.scatter_field(incident, 0.01, &gate).unwrap();
        let free = incident
            .pressure_velocity_at(
                1.0,
                Medium {
                    density: 1.0,
                    sound_speed: 1.0,
                },
                [0.0; 3],
            )
            .unwrap()
            .pressure_pa;
        // At the centre only l=0 survives. Every regular incident field has
        // the same independent sphere coefficient multiplying its centre value.
        let expected = monopole_transfer * free;
        let actual = field
            .solution()
            .interior_pressure_at(&[[0.0; 3]], 0.01, near_field::Options::default())
            .unwrap()
            .pressure[0];
        assert!(
            relative(actual, expected) < 0.04,
            "actual={actual:?}, expected={expected:?}"
        );
        assert!(
            field
                .solution()
                .scattered()
                .pressure
                .iter()
                .any(|p| p.abs() > 1e-5)
        );
    }
    assert_eq!(exterior.queries.load(Ordering::Relaxed), 1);
    assert_eq!(interior.queries.load(Ordering::Relaxed), 1);
}

#[test]
fn parallel_lu_descriptor_workspace_is_included_before_dense_allocation() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let exterior = fluid(1.0, 1.0);
    let interior = fluid(2.0, 1.0);
    let gate = CancelGate::new();
    let m = 2 * surface.areas().len();
    let complex_bytes = std::mem::size_of::<C64>();
    let usize_bytes = std::mem::size_of::<usize>();
    let old_bound = (2 * m * m + 6 * m) * complex_bytes + m * usize_bytes;
    let workspace = fs_la::eigen_complex::lu_complex_blocked_workspace_bytes(m, 64).unwrap();
    let expected = 2 * m * m * complex_bytes + m * usize_bytes + m * complex_bytes + workspace;
    assert!(expected > old_bound);
    assert_eq!(
        fs_la::eigen_complex::lu_complex_blocked_workspace_bytes(64, 64),
        Some(0)
    );
    assert_eq!(
        fs_la::eigen_complex::lu_complex_blocked_workspace_bytes(usize::MAX, usize::MAX),
        None
    );
    let result = fs_bem::helmholtz::with_dense_workers(64, || {
        FluidInterface::new(
            &surface,
            &exterior,
            &interior,
            1.0,
            Options {
                maximum_bytes: old_bound,
                ..Options::default()
            },
            &gate,
        )
    });
    match result {
        Err(TransmissionError::Helmholtz(HelmholtzError::CacheBudget {
            required_bytes,
            maximum_bytes,
        })) => {
            assert_eq!(required_bytes, expected);
            assert_eq!(maximum_bytes, old_bound);
        }
        _ => panic!("64-worker descriptor storage must refuse the old serial-sized budget"),
    }
}
