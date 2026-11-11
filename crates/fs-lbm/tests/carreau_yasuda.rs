//! G1/G3: the declared Yasuda crossover survives local collision and SI scaling.
//! All material numbers are synthetic declarations, not measured data.

use fs_lbm::rheology::{Rheology, RheologyField, TAU_CAP, TAU_FLOOR, channel_flow, update_tau};
use fs_lbm::{CS2, Cell, Grid, Q, equilibrium};
use fs_material::fluid_transport::rheology::CarreauYasuda;

// The D2Q9 lattice velocities and weights, in the order of `fs_lbm::equilibrium`. The crate keeps its own `E` and `W` crate-private (they have been since the
// crate's first commit), so this test carries them; `the_local_lattice_follows_the_crates_equilibrium_order` pins the order against the public `equilibrium`.
const E: [(i32, i32); Q] = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (-1, 1), (-1, -1), (1, -1)];
const W: [f64; Q] = [4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0];

#[test]
fn the_local_lattice_follows_the_crates_equilibrium_order() {
    for (rho, ux, uy) in [(1.0, 0.0, 0.0), (1.3, 0.04, -0.02), (0.8, -0.03, 0.05)] {
        let f = equilibrium(rho, ux, uy);
        for (q, e) in E.iter().enumerate() {
            let eu = f64::from(e.0) * ux + f64::from(e.1) * uy;
            let expected = W[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * (ux * ux + uy * uy));
            assert!((f[q] - expected).abs() <= 1e-14, "population {q}: {} against {expected}", f[q]);
        }
    }
}

fn solve_stress(law: Rheology, rate: f64, expected_tau: f64, previous: f64) -> (f64, f64) {
    let mut grid = Grid::uniform(2, 2, previous);
    let mut populations = equilibrium(1.0, 0.0, 0.0);
    for (q, e) in E.iter().enumerate() {
        populations[q] -= 3.0 * W[q] * expected_tau * rate * f64::from(e.0 * e.1);
    }
    grid.f.fill(populations);
    let stats = update_tau(&mut grid, law);
    assert_eq!((stats.floored, stats.capped), (0, 0));
    (grid.tau[0], stats.max_shear)
}

#[test]
fn yasuda_retains_original_material_curve_and_local_stress_root() {
    for a in [0.3, 0.7, 1.0, 2.0, 4.0] {
        let source = CarreauYasuda::new(0.3, 0.01, 40.0, 0.35, a, [0.0, 1e300]).unwrap();
        let law = Rheology::CarreauYasuda {
            nu0: 0.3,
            nu_inf: 0.01,
            lambda: 40.0,
            n: 0.35,
            a,
        };
        for rate in [0.0, 1e-20, 1e-6, 0.01, 1.0, 100.0, 1e10] {
            let expected = source.viscosity_pa_s(rate).unwrap();
            let actual = law.viscosity(rate);
            assert!(
                (actual / expected - 1.0).abs() < 2e-13,
                "a={a}, rate={rate}"
            );
            if !(1e-6..=0.1).contains(&rate) {
                continue;
            }
            let tau = 0.5 + expected / CS2;
            for previous in [TAU_FLOOR, 1.1, TAU_CAP] {
                let (solved, shear) = solve_stress(law, rate, tau, previous);
                assert!((solved - tau).abs() < 4e-12);
                assert!((shear - rate).abs() < 5e-13);
            }
        }
    }
}

#[test]
fn yasuda_keeps_constant_and_carreau_fast_paths_and_refuses_invalid_transition() {
    for (nu0, nu_inf, lambda, n) in [
        (0.2, 0.01, 40.0, 0.4),
        (0.1, 0.3, 1.0, 0.6),
        (0.2, 0.01, 4.0, 1.5),
    ] {
        let carreau = Rheology::Carreau {
            nu0,
            nu_inf,
            lambda,
            n,
        };
        let yasuda = Rheology::CarreauYasuda {
            nu0,
            nu_inf,
            lambda,
            n,
            a: 2.0,
        };
        for rate in [0.0, 1e-20, 1e-6, 0.01, 1.0, 100.0] {
            assert_eq!(
                yasuda.viscosity(rate).to_bits(),
                carreau.viscosity(rate).to_bits()
            );
            if rate <= 0.01 {
                let tau = 0.5 + carreau.viscosity(rate) / CS2;
                assert_eq!(
                    solve_stress(yasuda, rate, tau, 1.4),
                    solve_stress(carreau, rate, tau, 1.4)
                );
            }
        }
    }
    for (nu0, nu_inf, lambda, n) in [
        (0.2, 0.01, 0.0, 0.4),
        (0.2, 0.01, 40.0, 1.0),
        (0.2, 0.2, 40.0, 0.4),
    ] {
        let law = Rheology::CarreauYasuda {
            nu0,
            nu_inf,
            lambda,
            n,
            a: 0.7,
        };
        let tau = 0.5 + nu0 / CS2;
        assert_eq!(solve_stress(law, 0.01, tau, 2.7).0, tau);
    }
    for a in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let law = Rheology::CarreauYasuda {
            nu0: 0.2,
            nu_inf: 0.01,
            lambda: 40.0,
            n: 0.4,
            a,
        };
        assert!(RheologyField::new([1, 1], &[law]).is_err());
    }
}

#[test]
fn yasuda_extreme_transition_preserves_representable_viscosity() {
    for (nu0, nu_inf, lambda, n, a, rate) in [
        (0.3, 0.01, 1e300, 0.999, 0.7, 1e300),
        (0.3, 0.01, 1e300, 0.999, 2.0, 1e300),
        (0.3, 0.01, 1e-300, 0.35, 0.001, 1e-300),
        (1e300, 1e-250, 1e300, 0.35, 4.0, 1e300),
        (0.3, 0.01, 1.0, 0.35, 1e300, 10.0),
        (0.3, 0.01, 1.0, 0.35, 1e-300, 10.0),
    ] {
        let source = CarreauYasuda::new(nu0, nu_inf, lambda, n, a, [0.0, 1e300]).unwrap();
        let expected = source.viscosity_pa_s(rate).unwrap();
        let law = Rheology::CarreauYasuda {
            nu0,
            nu_inf,
            lambda,
            n,
            a,
        };
        let actual = law.viscosity(rate);
        assert!(actual.is_finite());
        assert!(
            (actual / expected - 1.0).abs() < 2e-12,
            "{law:?}, rate={rate}, {actual} != {expected}"
        );
    }
    // For a=1 and n=2 the curve is exactly nu0 + (nu0-nu_inf)*lambda*rate.
    // Form the reference without creating the overflowing shape.
    let thickening = Rheology::CarreauYasuda {
        nu0: 2e-200,
        nu_inf: 1e-200,
        lambda: 1e200,
        n: 2.0,
        a: 1.0,
    };
    let expected = 2e-200 + ((2e-200 - 1e-200) * 1e200) * 1e200;
    let actual = thickening.viscosity(1e200);
    assert!(actual.is_finite());
    assert!((actual / expected - 1.0).abs() < 2e-12);
}

#[test]
fn yasuda_physical_stress_survives_rescaling_and_viscosity_basis() {
    use fs_lbm::material::plan_scaling_from_material;
    use fs_lbm::rheology::physical::PhysicalRheology;
    use fs_material::fluid_transport::FluidTransportState;
    use fs_material::liquid_temperature::LiquidViscosityBasis;

    let source = CarreauYasuda::new(0.4, 0.02, 0.03, 0.4, 0.7, [0.0, 300.0]).unwrap();
    let density = 800.0;
    let material = FluidTransportState::try_new(density, 2000.0, 0.2, 0.08).unwrap();
    for refinement in [1.0, 2.0, 4.0] {
        let dx = 0.002 / refinement;
        let dt = 0.002 / (refinement * refinement);
        let scale = plan_scaling_from_material(&material, dx, dt, 0.01).unwrap();
        for basis in [
            LiquidViscosityBasis::Dynamic,
            LiquidViscosityBasis::Kinematic,
        ] {
            let divisor = if basis == LiquidViscosityBasis::Dynamic {
                1.0
            } else {
                density
            };
            let bound = scale
                .bind_rheology(
                    PhysicalRheology::CarreauYasuda {
                        viscosity_zero: 0.4 / divisor,
                        viscosity_infinite: 0.02 / divisor,
                        relaxation_time_s: 0.03,
                        n: 0.4,
                        a: 0.7,
                    },
                    basis,
                )
                .unwrap();
            assert!(core::ptr::eq(bound.scaling().material(), &material));
            for rate in [0.0, 0.01, 1.0, 20.0, 300.0] {
                let mu = source.viscosity_pa_s(rate).unwrap();
                let lattice_rate = rate * dt;
                let nu = bound.lattice_law().viscosity(lattice_rate);
                let stress = density * (dx / dt).powi(2) * nu * lattice_rate;
                assert!((stress - mu * rate).abs() < 5e-13 * (mu * rate).abs().max(1e-12));
            }
        }
    }
}

// Independent continuum momentum balance: nu(gamma)*gamma = |g*y|.
// Invert this monotone stress with bisection and integrate gamma from y to
// the stationary wall. Platform powf and Simpson integration are oracles,
// independent of the deterministic production curve and kinetic solver.
fn channel_reference(y: f64, half_height: f64, g: f64, a: f64) -> f64 {
    let shear = |coordinate: f64| {
        let stress = g * coordinate;
        let (mut lo, mut hi) = (0.0, stress / 0.03);
        for _ in 0..70 {
            let rate = 0.5 * (lo + hi);
            let nu = 0.03 + 0.12 * (1.0 + (5000.0 * rate).powf(a)).powf(-0.55 / a);
            if nu * rate > stress {
                hi = rate;
            } else {
                lo = rate;
            }
        }
        0.5 * (lo + hi)
    };
    let start = y.abs();
    let panels = 256;
    let dy = (half_height - start) / panels as f64;
    let mut sum = shear(start) + shear(half_height);
    for j in 1..panels {
        sum += if j % 2 == 0 { 2.0 } else { 4.0 } * shear(start + j as f64 * dy);
    }
    sum * dy / 3.0
}

#[test]
fn yasuda_channel_matches_continuum_in_both_native_dimensions() {
    let (height, g, a) = (25, 2e-6, 0.7);
    let law = Rheology::CarreauYasuda {
        nu0: 0.15,
        nu_inf: 0.03,
        lambda: 5000.0,
        n: 0.45,
        a,
    };
    let max_steps = 40_000;
    let (grid2, steps2, stats2) = channel_flow(2, height, g, law, max_steps);
    let (grid3, steps3, stats3) =
        fs_lbm::d3q19::plate_channel_flow3(2, 1, height, g, law, max_steps);
    assert!(steps2 < max_steps && steps3 < max_steps);
    assert_eq!(
        (stats2.floored, stats2.capped, stats3.floored, stats3.capped),
        (0, 0, 0, 0)
    );
    let center = (height - 1) as f64 / 2.0;
    let half_height = height as f64 / 2.0;
    let peak = channel_reference(0.0, half_height, g, a);
    // This broad transition differs materially from a silently substituted
    // fixed-a Carreau curve with otherwise identical material parameters.
    assert!((peak - channel_reference(0.0, half_height, g, 2.0)).abs() > 0.08 * peak);
    let mut errors = [0.0_f64; 2];
    for y in 0..height {
        let exact = channel_reference(y as f64 - center, half_height, g, a);
        let u2 = grid2.moments(grid2.idx(0, y + 1)).u[0];
        let u3 = grid3.moments3(grid3.idx(0, 0, y + 1)).1[0];
        errors[0] = errors[0].max((u2 - exact).abs());
        errors[1] = errors[1].max((u3 - exact).abs());
        assert!((u2 - u3).abs() < 1e-9);
    }
    assert!(
        errors.into_iter().all(|error| error / peak < 0.02),
        "relative errors={:?}",
        errors.map(|error| error / peak)
    );
    // `Grid::total_mass` sums the populations of the cells that are neither Gas nor Wall (core2.rs:1074, b835c8ff91,
    // 2026-07-08), and `channel_flow` makes its two wall rows Wall, so the 2-D mass is the fluid cells at density 1, not
    // `2 (height + 2)` (f782efa44f, 2026-11-08, counted the wall rows).
    let fluid_cells = grid2
        .flags
        .iter()
        .filter(|flag| !matches!(flag, Cell::Gas | Cell::Wall))
        .count();
    assert_eq!(fluid_cells, 2 * height);
    assert!((grid2.total_mass() - fluid_cells as f64).abs() < 1e-9);
    assert!((grid3.total_mass() - (2 * (height + 2)) as f64).abs() < 1e-9);
    println!(
        "Carreau-Yasuda a={a}: 2D={steps2}, 3D={steps3} steps, relative profile errors={:?}",
        errors.map(|error| error / peak)
    );
}
