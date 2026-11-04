//! G0 link conservation and G1 Fourier resistance checks on the original
//! D2Q5/D3Q7 collision and streaming owners, including nonequilibrium states.

use super::*;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g0_wall_admission_preserves_signed_scalar_coordinates_and_flux() {
    for wall in [
        ThermalWall::Temperature(-3.0),
        ThermalWall::Adiabatic,
        ThermalWall::HeatFlux(-0.2),
        ThermalWall::HeatFlux(0.0),
        ThermalWall::Convective {
            transfer: 0.0,
            ambient: -1.0,
        },
        ThermalWall::Convective {
            transfer: f64::MAX,
            ambient: f64::MAX,
        },
    ] {
        assert!(wall.validate().is_ok());
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for wall in [
            ThermalWall::Temperature(value),
            ThermalWall::HeatFlux(value),
            ThermalWall::Convective {
                transfer: value,
                ambient: 1.0,
            },
            ThermalWall::Convective {
                transfer: 0.0,
                ambient: value,
            },
        ] {
            assert!(wall.validate().is_err());
        }
    }
    assert!(
        ThermalWall::Convective {
            transfer: -0.01,
            ambient: 1.0
        }
        .validate()
        .is_err()
    );
}

#[test]
fn g0_wall_pairs_obey_the_prescribed_flux_and_implicit_film_law() {
    for weight in [1.0 / 6.0, 1.0 / 8.0] {
        for outgoing in [-0.03, 0.02, 0.4] {
            let (incoming, wall) = ThermalWall::Adiabatic.reflect(outgoing, weight);
            assert_eq!(incoming.to_bits(), outgoing.to_bits());
            near(wall, outgoing / weight, 2e-15);
            for q in [-0.01, 0.0, 0.03] {
                let (incoming, wall) = ThermalWall::HeatFlux(q).reflect(outgoing, weight);
                near(incoming - outgoing, q, 2e-16);
                near(weight * wall, outgoing + 0.5 * q, 2e-16);
            }
            for transfer in [0.0, 1e-8, 0.02, 0.25, 1.0, 4.0] {
                let ambient = 1.3;
                let (incoming, wall) =
                    ThermalWall::Convective { transfer, ambient }.reflect(outgoing, weight);
                near(incoming - outgoing, transfer * (ambient - wall), 4e-15);
                if transfer == 0.0 {
                    assert_eq!(incoming.to_bits(), outgoing.to_bits());
                }
            }
        }
    }
}

#[test]
fn g0_fixed_temperature_retains_original_bits_and_film_limits_stay_finite() {
    for weight in [1.0 / 6.0, 1.0 / 8.0] {
        for temperature in [-1.2, 0.0, 0.7, 1.0] {
            for outgoing in [-0.1, 0.01, 0.6] {
                let (actual, trace) =
                    ThermalWall::Temperature(temperature).reflect(outgoing, weight);
                let original = 2.0f64.mul_add(weight * temperature, -outgoing);
                assert_eq!(actual.to_bits(), original.to_bits());
                assert_eq!(trace.to_bits(), temperature.to_bits());
                let (strong, trace) = ThermalWall::Convective {
                    transfer: f64::MAX,
                    ambient: temperature,
                }
                .reflect(outgoing, weight);
                assert!(strong.is_finite() && trace.is_finite());
                near(strong, original, 2e-15);
                near(trace, temperature, 2e-15);
            }
        }
        // H=2w makes the reflected population exactly w*T_ambient. Its
        // convex combination is representable even if the naive difference
        // w*T_ambient-g_out overflows. The trace itself is unsupported here.
        let (incoming, _) = ThermalWall::Convective {
            transfer: 2.0 * weight,
            ambient: f64::MAX,
        }
        .reflect(-f64::MAX, weight);
        assert_eq!(incoming, weight * f64::MAX);
    }
}

fn storage<const Q: usize>(populations: &[[f64; Q]], plane: usize) -> f64 {
    populations[plane..populations.len() - plane]
        .iter()
        .flatten()
        .sum()
}

fn equilibrate<const D: usize, const Q: usize>(
    lattice: &AxisLattice<D, Q>,
    shape: [usize; D],
    tau: f64,
    walls: [ThermalWall; 2],
) -> (Vec<[f64; Q]>, [ThermalWallExchange; 2]) {
    let cells: usize = shape.iter().product();
    let plane = cells / shape[D - 1];
    let mut populations = vec![lattice.equilibrium(1.0, [0.0; D]); cells];
    let mut post = vec![0.0; Q * cells];
    let mut history = FluxHistory::default();
    let mut exchange = [ThermalWallExchange::default(); 2];
    let initial = storage(&populations, plane);
    let mut boundary_heat = 0.0;
    // The insulated/flux end has the slower quarter-wave diffusion mode.
    // At the smallest D3Q7 diffusivity, 10k steps leave a 3.5e-9 transient;
    // keep the analytic tolerance and let that physical mode decay instead.
    for _ in 0..20_000 {
        lattice.collide(
            &populations,
            &mut post,
            plane..cells - plane,
            tau,
            None,
            &mut history,
            |_| [0.0; D],
        );
        exchange = lattice.stream_with_walls(&mut populations, &post, shape, Some(walls));
        boundary_heat += exchange[0].heat_in + exchange[1].heat_in;
    }
    near(storage(&populations, plane) - initial, boundary_heat, 3e-10);
    (populations, exchange)
}

fn two_films<const D: usize, const Q: usize>(
    lattice: &AxisLattice<D, Q>,
    shape: [usize; D],
    tau: f64,
) {
    let n = shape[D - 1] - 2;
    let plane: usize = shape[..D - 1].iter().product();
    let alpha = lattice.cs2 * (tau - 0.5);
    let (hot, cold, bottom_h, top_h) = (1.4, 0.4, 0.03, 0.08);
    let flux = (hot - cold) / (1.0 / bottom_h + n as f64 / alpha + 1.0 / top_h);
    let bottom = hot - flux / bottom_h;
    let top = cold + flux / top_h;
    let (populations, exchange) = equilibrate(
        lattice,
        shape,
        tau,
        [
            ThermalWall::Convective {
                transfer: bottom_h,
                ambient: hot,
            },
            ThermalWall::Convective {
                transfer: top_h,
                ambient: cold,
            },
        ],
    );
    for layer in 1..=n {
        let expected = bottom - flux * (layer as f64 - 0.5) / alpha;
        for pop in &populations[plane * layer..plane * (layer + 1)] {
            near(pop.iter().sum(), expected, 2e-11);
        }
    }
    for (side, expected, heat) in [(0, bottom, flux), (1, top, -flux)] {
        near(exchange[side].heat_in, heat * plane as f64, 2e-12);
        near(exchange[side].mean_temperature, expected, 2e-11);
        near(exchange[side].min_temperature, expected, 2e-11);
        near(exchange[side].max_temperature, expected, 2e-11);
    }
}

#[test]
fn g1_two_robin_walls_recover_fourier_resistances_in_both_native_lattices() {
    for tau in [0.61, 1.0, 1.7] {
        two_films(&D2Q5, [3, 8], tau);
        two_films(&D3Q7, [3, 2, 8], tau);
    }
}

fn flux_wall<const D: usize, const Q: usize>(
    lattice: &AxisLattice<D, Q>,
    shape: [usize; D],
    tau: f64,
    flux: f64,
) {
    let n = shape[D - 1] - 2;
    let plane: usize = shape[..D - 1].iter().product();
    let alpha = lattice.cs2 * (tau - 0.5);
    let top = 0.8;
    let bottom = top + flux * n as f64 / alpha;
    let (populations, exchange) = equilibrate(
        lattice,
        shape,
        tau,
        [ThermalWall::HeatFlux(flux), ThermalWall::Temperature(top)],
    );
    for layer in 1..=n {
        let expected = bottom - flux * (layer as f64 - 0.5) / alpha;
        for pop in &populations[plane * layer..plane * (layer + 1)] {
            near(pop.iter().sum(), expected, 2e-11);
        }
    }
    near(exchange[0].mean_temperature, bottom, 2e-11);
    near(exchange[1].mean_temperature, top, 2e-15);
    near(exchange[0].heat_in, flux * plane as f64, 2e-13);
    near(exchange[1].heat_in, -flux * plane as f64, 2e-12);
}

#[test]
fn g1_prescribed_flux_has_no_spurious_relaxation_factor_or_wall_shift() {
    for tau in [0.61, 1.0, 1.7] {
        for flux in [-0.002, 0.004] {
            flux_wall(&D2Q5, [3, 8], tau, flux);
            flux_wall(&D3Q7, [3, 2, 8], tau, flux);
        }
    }
}

fn insulated<const D: usize, const Q: usize>(lattice: &AxisLattice<D, Q>, shape: [usize; D]) {
    let cells: usize = shape.iter().product();
    let plane = cells / shape[D - 1];
    let mut populations: Vec<[f64; Q]> = (0..cells)
        .map(|cell| {
            let mut pop = lattice.equilibrium(0.4 + 0.002 * cell as f64, [0.0; D]);
            // Finite non-equilibrium normal flux; insulation must conserve it too.
            pop[lattice.positive[D - 1]] += 0.003;
            pop[lattice.negative[D - 1]] -= 0.003;
            pop
        })
        .collect();
    let initial = storage(&populations, plane);
    let ghosts = [
        populations[..plane].to_vec(),
        populations[cells - plane..].to_vec(),
    ];
    let mut post = vec![0.0; Q * cells];
    let mut history = FluxHistory::default();
    for _ in 0..300 {
        lattice.collide(
            &populations,
            &mut post,
            plane..cells - plane,
            0.73,
            None,
            &mut history,
            |_| [0.0; D],
        );
        let report = lattice.stream_with_walls(
            &mut populations,
            &post,
            shape,
            Some([ThermalWall::Adiabatic; 2]),
        );
        assert_eq!(report[0].heat_in, 0.0);
        assert_eq!(report[1].heat_in, 0.0);
        for side in report {
            assert!(side.min_temperature <= side.mean_temperature + 2e-15);
            assert!(side.max_temperature >= side.mean_temperature - 2e-15);
        }
    }
    near(storage(&populations, plane), initial, 2e-11);
    assert_eq!(&populations[..plane], &ghosts[0]);
    assert_eq!(&populations[cells - plane..], &ghosts[1]);
}

#[test]
fn g0_insulated_nonequilibrium_transport_conserves_storage_and_ghosts() {
    insulated(&D2Q5, [5, 8]);
    insulated(&D3Q7, [3, 4, 8]);
}

#[test]
fn g0_trace_extrema_retain_nonfinite_refusals_after_later_finite_links() {
    let shape = [3, 4];
    let cells = 12;
    let mut populations = vec![[0.0; 5]; cells];
    let mut post = vec![0.1; 5 * cells];
    // Finite incoming/outgoing populations can still imply an unrepresentable
    // surface temperature. Later links must not hide it via f64::min/max.
    post[5 * 3 + D2Q5.negative[1]] = f64::MAX;
    let report = D2Q5.stream_with_walls(
        &mut populations,
        &post,
        shape,
        Some([ThermalWall::Adiabatic; 2]),
    );
    assert!(report[0].min_temperature.is_nan());
    assert!(report[0].max_temperature.is_nan());
    near(report[1].min_temperature, 0.6, 2e-15);
    near(report[1].max_temperature, 0.6, 2e-15);
}
