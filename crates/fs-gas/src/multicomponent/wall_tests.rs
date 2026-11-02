//! Independent physical oracles for a reflecting wall at local mixture composition.
use crate::multicomponent::tests::{gas, source};
use crate::multicomponent::{
    Boundaries, Boundary, Conserved, Error, FaceFlux, MixtureEuler1d, MixtureGas, NormalCell,
    NormalPrimitive, SpatialCell, SpatialPrimitive, StepControl,
};
use fs_qty::{Length, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn caloric_gas(variable: bool, reference: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, if variable { 0.25 } else { 0.0 }, reference),
        source(
            "B",
            2.0 * R,
            if variable { 5.0 } else { 3.5 },
            if variable { 0.7 } else { 0.0 },
            -2.0 * reference,
        ),
    ])
}

fn reflecting_flux(gas: &MixtureGas, inside: &Conserved, normal: f64) -> ([f64; 4], FaceFlux) {
    let primitive = gas.primitive(inside.as_cell(), || false).unwrap();
    let mut flux = [f64::NAN; 4];
    let receipt = gas
        .wall_flux(inside.as_cell(), primitive, normal, &mut flux, &mut || {
            false
        })
        .unwrap();
    assert_eq!(flux[0], 0.0, "species A crossed a stationary wall");
    assert_eq!(flux[1], 0.0, "species B crossed a stationary wall");
    assert_eq!(
        flux[3], 0.0,
        "a stationary adiabatic wall transferred energy"
    );
    assert!(flux[2] > 0.0, "wall traction must not be tensile");
    assert!(receipt.signal_speed.is_finite() && receipt.signal_speed > 0.0);
    assert!(!receipt.used_hll_fallback);
    (flux, receipt)
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.15e}, expected={expected:.15e}, scale={scale:.9e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn moving_mixture_piston_closes_species_and_variable_caloric_work() {
    for densities in [[0.24_f64, 0.96], [0.96, 0.24], [0.0, 1.2], [1.2, 0.0]] {
        let rho1 = densities[0] + densities[1];
        let (a, b, r) = mixture_coefficients(densities);
        let (t1, speed) = (2.5_f64, 3.2_f64);
        let p1 = rho1 * r * t1;
        let mut expected = None;
        for reference in [0.0, -1e18, 1e18] {
            let gas = caloric_gas(true, reference);
            for normal in [-1.0, 1.0] {
                let wall_velocity = -normal * speed;
                let mut flux = [f64::NAN; 6];
                let report = gas
                    .moving_wall_spatial_flux_from_primitive(
                        SpatialPrimitive {
                            partial_densities: &densities,
                            normal_velocity: 0.0,
                            transverse_velocity: [3.0, -4.0],
                            pressure: p1,
                        },
                        normal,
                        wall_velocity,
                        &mut flux,
                        || false,
                    )
                    .unwrap();
                assert_eq!(&flux[..2], &[0.0; 2]);
                assert_eq!(&flux[3..5], &[0.0; 2]);
                assert_eq!(flux[5], flux[2] * wall_velocity);
                assert!(!report.used_hll_fallback);
                let p2 = flux[2];
                let w = (p2 - p1) / (rho1 * speed) - speed;
                let rho2 = rho1 * (speed + w) / w;
                let t2 = p2 / (rho2 * r);
                let delta_e = (a - r) * (t2 - t1) + 0.5 * b * (t2 * t2 - t1 * t1);
                // The exact translating-piston shock gives a uniform swept
                // slab. Its source-integrated thermal plus kinetic energy gain
                // must equal the laboratory pressure work; tangential energy
                // and formation references cancel by each species' mass law.
                let energy_rate = rho1 * (speed + w) * (delta_e + 0.5 * speed * speed);
                close(energy_rate, -normal * flux[5], energy_rate, 4e-9);
                for density in densities {
                    close(
                        density * (speed + w),
                        density * (rho2 / rho1) * w,
                        1.0,
                        4e-9,
                    );
                }
                let entropy = a * (t2 / t1).ln() + b * (t2 - t1) - r * (p2 / p1).ln();
                assert!(entropy > 0.0);
                let invariant = (p2.to_bits(), report.signal_speed.to_bits());
                if let Some(previous) = expected {
                    assert_eq!(invariant, previous);
                } else {
                    expected = Some(invariant);
                }
            }
        }
    }
}

#[test]
fn moving_mixture_wall_reduces_to_stationary_and_preserves_frame_covariance() {
    let gas = caloric_gas(true, 0.0);
    let densities = [0.25, 0.75];
    let planar = NormalPrimitive {
        partial_densities: &densities,
        normal_velocity: 0.5,
        transverse_velocity: 2.0,
        pressure: 1.5,
    };
    let spatial = SpatialPrimitive {
        partial_densities: &densities,
        normal_velocity: planar.normal_velocity,
        transverse_velocity: [planar.transverse_velocity, 0.0],
        pressure: planar.pressure,
    };
    let mut q2 = [0.0; 5];
    gas.normal_conserved_into(planar, &mut q2).unwrap();
    let cell2 = NormalCell {
        partial_densities: &q2[..2],
        normal_momentum: q2[2],
        transverse_momentum: q2[3],
        energy: q2[4],
    };
    let mut q3 = [0.0; 6];
    gas.spatial_conserved_into(spatial, &mut q3).unwrap();
    let cell3 = SpatialCell {
        partial_densities: &q3[..2],
        normal_momentum: q3[2],
        transverse_momentum: [q3[3], q3[4]],
        energy: q3[5],
    };
    for normal in [-1.0, 1.0] {
        for wall_velocity in [0.0, -0.25, 1.0] {
            let mut f2 = [0.0; 5];
            let mut f3 = [0.0; 6];
            let r2 = gas
                .moving_wall_normal_flux_from_primitive(
                    planar,
                    normal,
                    wall_velocity,
                    &mut f2,
                    || false,
                )
                .unwrap();
            let r3 = gas
                .moving_wall_spatial_flux_from_primitive(
                    spatial,
                    normal,
                    wall_velocity,
                    &mut f3,
                    || false,
                )
                .unwrap();
            assert_eq!(r2, r3);
            assert_eq!(f2[2], f3[2]);
            assert_eq!(f2[4], f3[5]);
            assert_eq!(&f3[..2], &[0.0; 2]);
            assert_eq!(&f3[3..5], &[0.0; 2]);
            let mut decoded2 = [0.0; 5];
            let mut decoded3 = [0.0; 6];
            let d2 = gas
                .moving_wall_normal_flux(cell2, normal, wall_velocity, &mut decoded2, || false)
                .unwrap();
            let d3 = gas
                .moving_wall_spatial_flux(cell3, normal, wall_velocity, &mut decoded3, || false)
                .unwrap();
            close(decoded2[2], f2[2], f2[2], 4e-9);
            close(decoded3[2], f3[2], f3[2], 4e-9);
            close(decoded2[4], f2[4], f2[2], 4e-9);
            close(decoded3[5], f3[5], f3[2], 4e-9);
            let boost = 4.0;
            let translated = SpatialPrimitive {
                normal_velocity: spatial.normal_velocity + boost,
                transverse_velocity: [5.0, -2.0],
                ..spatial
            };
            let mut boosted = [0.0; 6];
            let br = gas
                .moving_wall_spatial_flux_from_primitive(
                    translated,
                    normal,
                    wall_velocity + boost,
                    &mut boosted,
                    || false,
                )
                .unwrap();
            assert_eq!(br, r3);
            assert_eq!(boosted[2], f3[2]);
            close(boosted[5], f3[5] + boost * f3[2], f3[2], 2e-14);
            if wall_velocity == 0.0 {
                let mut stationary2 = [0.0; 5];
                let mut stationary3 = [0.0; 6];
                assert_eq!(
                    r2,
                    gas.stationary_wall_normal_flux_from_primitive(
                        planar,
                        normal,
                        &mut stationary2,
                        || false
                    )
                    .unwrap()
                );
                assert_eq!(
                    r3,
                    gas.stationary_wall_spatial_flux_from_primitive(
                        spatial,
                        normal,
                        &mut stationary3,
                        || false
                    )
                    .unwrap()
                );
                assert_eq!(stationary2, f2);
                assert_eq!(stationary3, f3);
                assert_eq!(
                    d2,
                    gas.stationary_wall_normal_flux(cell2, normal, &mut stationary2, || false)
                        .unwrap()
                );
                assert_eq!(
                    d3,
                    gas.stationary_wall_spatial_flux(cell3, normal, &mut stationary3, || false)
                        .unwrap()
                );
                assert_eq!(stationary2, decoded2);
                assert_eq!(stationary3, decoded3);
            }
        }
    }
}

#[test]
fn moving_mixture_wall_refuses_without_publishing_partial_work() {
    let gas = caloric_gas(false, 0.0);
    let densities = [0.7, 0.3];
    let input = SpatialPrimitive {
        partial_densities: &densities,
        normal_velocity: 0.0,
        transverse_velocity: [0.0; 2],
        pressure: 0.4,
    };
    let original = [19.0; 6];
    let mut output = original;
    for velocity in [-30.0, 4.0, f64::NAN, f64::INFINITY] {
        assert!(
            gas.moving_wall_spatial_flux_from_primitive(input, 1.0, velocity, &mut output, || {
                false
            })
            .is_err()
        );
        assert_eq!(output, original);
    }
    assert!(
        gas.moving_wall_spatial_flux_from_primitive(input, 1.0, -1.0, &mut output[..5], || false)
            .is_err()
    );
    assert_eq!(output, original);
    assert!(
        gas.moving_wall_spatial_flux_from_primitive(input, 0.0, -1.0, &mut output, || false)
            .is_err()
    );
    assert_eq!(output, original);
    let mut count = 0;
    gas.moving_wall_spatial_flux_from_primitive(input, 1.0, -1.0, &mut output, || {
        count += 1;
        false
    })
    .unwrap();
    for stop in [1, count / 2, count] {
        output = original;
        let mut polls = 0;
        assert_eq!(
            gas.moving_wall_spatial_flux_from_primitive(input, 1.0, -1.0, &mut output, || {
                polls += 1;
                polls == stop
            }),
            Err(Error::Cancelled)
        );
        assert_eq!(output, original);
    }
}

#[test]
fn arbitrary_constant_cp_compositions_match_analytic_shock_and_expansion_pressures() {
    let gas = caloric_gas(false, 0.0);
    let gamma = 1.4;
    // Both sources have Cp/R=3.5, so gamma is exactly 1.4 even though their
    // molar masses and the local mixture gas constants are different.
    for densities in [[0.7_f64, 0.3], [0.05, 1.75], [1.3, 0.0], [0.0, 0.8]] {
        let density = densities[0] + densities[1];
        let pressure = 1.2 * density;
        for velocity in [1e-6_f64, 0.2, 2.0, 6.0] {
            let q = gas.conserved(&densities, velocity, pressure).unwrap();
            let (flux, _) = reflecting_flux(&gas, &q, 1.0);
            // Exact gamma-law shock curve f(p*)=u, solved algebraically.
            let k = 0.5 * (gamma + 1.0) * density * velocity * velocity;
            let b = (gamma - 1.0) / (gamma + 1.0) * pressure;
            let exact = pressure + 0.5 * k + 0.5 * (k * k + 4.0 * k * (pressure + b)).sqrt();
            close(flux[2], exact, exact, 3e-9);
        }
        let pressure = 0.4 * density;
        let sound_speed = (gamma * pressure / density).sqrt();
        for speed in [0.01_f64, 0.5, 1.5, 3.0] {
            let q = gas.conserved(&densities, speed, pressure).unwrap();
            let (flux, _) = reflecting_flux(&gas, &q, -1.0);
            let factor = 1.0 - (gamma - 1.0) * speed / (2.0 * sound_speed);
            let exact = pressure * factor.powf(2.0 * gamma / (gamma - 1.0));
            close(flux[2], exact, pressure, 3e-9);
        }
    }
}

// For the declared two distinct source polynomials, the mass-specific law is
// Cp(T)=a+bT, R=Y_A+Y_B/2. This oracle directly integrates those polynomials;
// it does not call a frozen-gas or shared wall implementation.
fn mixture_coefficients(densities: [f64; 2]) -> (f64, f64, f64) {
    let rho = densities[0] + densities[1];
    let ya = densities[0] / rho;
    let yb = densities[1] / rho;
    (3.5 * ya + 2.5 * yb, 0.25 * ya + 0.35 * yb, ya + 0.5 * yb)
}

#[test]
fn different_variable_cp_species_satisfy_independent_wall_hugoniot_and_entropy() {
    let gas = caloric_gas(true, 0.0);
    for densities in [[0.24_f64, 0.96], [0.96, 0.24], [0.0, 1.2], [1.2, 0.0]] {
        let density = densities[0] + densities[1];
        let (a, b, r) = mixture_coefficients(densities);
        let t1 = 2.5;
        let pressure = density * r * t1;
        let velocity = 3.2;
        let q = gas.conserved(&densities, velocity, pressure).unwrap();
        let (flux, _) = reflecting_flux(&gas, &q, 1.0);
        let p2 = flux[2];
        // Infer shock speed, stopped density and temperature from independent
        // mass/momentum Rankine-Hugoniot equations, then test the energy jump.
        let shock_speed = (p2 - pressure) / (density * velocity) - velocity;
        assert!(shock_speed > 0.0);
        let rho2 = density * (velocity + shock_speed) / shock_speed;
        let t2 = p2 / (rho2 * r);
        let delta_e = (a - r) * (t2 - t1) + 0.5 * b * (t2 * t2 - t1 * t1);
        let energy_jump = 0.5 * (p2 + pressure) * (1.0 / density - 1.0 / rho2);
        close(delta_e, energy_jump, energy_jump.abs().max(1.0), 4e-9);
        let entropy_change = a * (t2 / t1).ln() + b * (t2 - t1) - r * (p2 / pressure).ln();
        assert!(entropy_change > 0.0 && t2 > t1 && rho2 > density);
        // A gamma frozen at the upstream temperature is an observably
        // different closure for these sources and this finite compression.
        let gamma1 = (a + b * t1) / (a - r + b * t1);
        let k = 0.5 * (gamma1 + 1.0) * density * velocity * velocity;
        let bp = (gamma1 - 1.0) / (gamma1 + 1.0) * pressure;
        let frozen_gamma = pressure + 0.5 * k + 0.5 * (k * k + 4.0 * k * (pressure + bp)).sqrt();
        assert!((p2 / frozen_gamma - 1.0).abs() > 1e-3);
    }
}

#[test]
fn local_mixture_rarefaction_matches_quadrature_and_preserves_reference_and_orientation() {
    let ordinary = caloric_gas(true, 0.0);
    for densities in [[0.2_f64, 0.8], [0.8, 0.2]] {
        let density = densities[0] + densities[1];
        let (a, b, r) = mixture_coefficients(densities);
        let t1 = 3.0;
        let speed = 0.8;
        let pressure = density * r * t1;
        let q = ordinary.conserved(&densities, -speed, pressure).unwrap();
        let (expected_flux, _) = reflecting_flux(&ordinary, &q, 1.0);
        let p2 = expected_flux[2];
        assert!(p2 < pressure);
        // Independently invert the analytic isentrope to obtain the wall T.
        let (mut lo, mut hi) = (0.001, t1);
        for _ in 0..80 {
            let t = 0.5 * (lo + hi);
            let entropy = a * (t / t1).ln() + b * (t - t1) - r * (p2 / pressure).ln();
            if entropy < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
        }
        let t2 = 0.5 * (lo + hi);
        // A separate fixed-grid Simpson integral evaluates int(Cp/c)dT.
        // Production uses bounded adaptive quadrature in sqrt(T) coordinates.
        let panels = 2000;
        let h = (t1 - t2) / panels as f64;
        let integrand = |t: f64| ((a + b * t) * (a - r + b * t) / (r * t)).sqrt();
        let mut integral = integrand(t2) + integrand(t1);
        for i in 1..panels {
            integral += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(t2 + i as f64 * h);
        }
        close(h * integral / 3.0, speed, speed, 4e-9);
        for toward_wall in [-speed, 2.0] {
            let q = ordinary
                .conserved(&densities, toward_wall, pressure)
                .unwrap();
            let (expected_flux, expected_receipt) = reflecting_flux(&ordinary, &q, 1.0);
            for reference in [-1e18, 0.0, 1e18] {
                let shifted = caloric_gas(true, reference);
                for normal in [-1.0, 1.0] {
                    let q = shifted
                        .conserved(&densities, normal * toward_wall, pressure)
                        .unwrap();
                    let (actual_flux, actual_receipt) = reflecting_flux(&shifted, &q, normal);
                    assert_eq!(actual_flux, expected_flux);
                    assert_eq!(
                        actual_receipt.signal_speed.to_bits(),
                        expected_receipt.signal_speed.to_bits()
                    );
                    assert_eq!(
                        actual_receipt.used_hll_fallback,
                        expected_receipt.used_hll_fallback
                    );
                }
            }
        }
    }
}

#[test]
fn retreating_wall_evolution_uses_positive_rarefaction_pressure_instead_of_tensile_hll() {
    let gas = caloric_gas(false, 0.0);
    let densities = [0.7, 0.3];
    let speed = 1.5;
    let pressure = 0.4;
    let state = gas.conserved(&densities, speed, pressure).unwrap();
    let mirrored = gas.conserved(&densities, -speed, pressure).unwrap();
    let mut ghost_flux = [0.0; 4];
    gas.hllc_flux(mirrored.as_cell(), state.as_cell(), &mut ghost_flux, || {
        false
    })
    .unwrap();
    // The previous wall path retained this momentum component and only zeroed
    // species and energy fluxes. At these data its HLL traction is negative.
    // Record it as a diagnostic; a future general flux improvement is allowed.
    eprintln!(
        "legacy mirrored mixture flux traction={:.9e}",
        ghost_flux[2]
    );
    let dx = 0.05;
    let mut flow = MixtureEuler1d::new(
        gas,
        Length::new(dx),
        0.7,
        vec![state.clone(); 20],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    let step = flow
        .step(Time::new(1e-8), StepControl::default(), || false)
        .unwrap();
    // A first-cell momentum update recovers the accepted wall traction at
    // dt -> 0 independently of the newly added direct wall API. The neighboring
    // interior flux is rho*u^2+p initially. RK evolution perturbs it by O(dt/dx).
    let inferred = (flow.cell(0).unwrap().momentum - state.momentum) * dx / step.elapsed.value()
        + speed * speed
        + pressure;
    let c = (1.4_f64 * pressure).sqrt();
    let exact = pressure * (1.0 - 0.2 * speed / c).powf(7.0);
    assert!(inferred > 0.0, "production wall still pulls: {inferred}");
    close(inferred, exact, pressure, 2e-6);
}

#[test]
fn reflecting_composition_transport_keeps_each_species_and_total_energy_inside() {
    let gas = caloric_gas(true, 0.0);
    let states: Vec<_> = (0..64)
        .map(|i| {
            gas.conserved(if i < 32 { &[0.8, 0.2] } else { &[0.1, 0.9] }, 1.2, 0.4)
                .unwrap()
        })
        .collect();
    let mut flow = MixtureEuler1d::new(
        gas,
        Length::new(1.0 / 64.0),
        0.7,
        states,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let initial = flow.totals();
    let mut impulse = 0.0;
    let end = 0.04;
    for _ in 0..10_000 {
        if flow.time().value() >= end {
            break;
        }
        let report = flow
            .step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(
            report.boundary_change.species_masses,
            [0.0; crate::multicomponent::MAX_SPECIES]
        );
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!(report.maximum_cfl <= 0.400_000_000_000_01);
        impulse += report.boundary_change.momentum + report.area_momentum_impulse;
        for i in 0..flow.cell_count() {
            let cell = flow.cell(i).unwrap();
            assert!(cell.partial_densities.iter().all(|rho| *rho >= 0.0));
            let primitive = flow.gas().primitive(cell, || false).unwrap();
            assert!(primitive.pressure > 0.0 && primitive.temperature > 0.0);
        }
    }
    assert_eq!(flow.time().value(), end);
    let final_state = flow.totals();
    for i in 0..flow.gas().species_count() {
        close(
            final_state.species_masses[i],
            initial.species_masses[i],
            1.0,
            2e-12,
        );
    }
    close(
        final_state.energy,
        initial.energy,
        initial.energy.abs().max(1.0),
        2e-11,
    );
    close(final_state.momentum - initial.momentum, impulse, 1.0, 2e-11);
    assert!(
        final_state.momentum < initial.momentum,
        "reflection did not decelerate the gas"
    );
}
