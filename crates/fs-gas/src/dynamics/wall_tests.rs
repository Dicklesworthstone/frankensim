use super::*;
use crate::dynamics::{CaloricDomain, NormalPrimitive, SpatialPrimitive};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};

fn gas(slope: f64, formation: f64) -> FrozenGas {
    let mut a = [0.0; 9];
    a[2] = 3.5;
    a[3] = slope;
    a[7] = formation;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 0.001, 50.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Manufactured thermally perfect wall reference".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("wall-reference").unwrap(),
    )
    .unwrap();
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let source = FrozenIdealGasMixtureModelV1::new(
        vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    FrozenGas::new(
        source,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn constant_cp_wall_compression_matches_analytic_shock_from_weak_to_strong() {
    let gas = gas(0.0, 0.0);
    for velocity in [1e-12_f64, 1e-5, 0.2, 2.0, 10.0] {
        let state = gas.conserved(1.0, velocity, 1.0).unwrap();
        let primitive = gas.primitive(state, || false).unwrap();
        let flux = gas.wall_flux(state, primitive, 1.0, &mut || false).unwrap();
        let k = 0.5 * 2.4 * velocity * velocity;
        let exact = 1.0 + 0.5 * k + 0.5 * (k * k + 4.0 * k * (1.0 + 1.0 / 6.0)).sqrt();
        close(flux.momentum, exact, 2e-10);
        assert_eq!(flux.mass, 0.0);
        assert_eq!(flux.energy, 0.0);
        assert!(!flux.used_hll_fallback);
    }
}

#[test]
fn constant_cp_retreat_matches_isentropic_rarefaction_and_never_pulls() {
    let gas = gas(0.0, 0.0);
    for speed in [0.01_f64, 0.5, 2.0, 3.0] {
        let state = gas.conserved(1.0, speed, 0.4).unwrap();
        let p = gas.primitive(state, || false).unwrap();
        let wall = gas.stopped_wall(p, -speed, &mut || false).unwrap();
        let factor = 1.0 - 0.2 * speed / (1.4_f64 * 0.4).sqrt();
        close(wall.pressure, 0.4 * factor.powf(7.0), 2e-10);
        close(wall.temperature, 0.4 * factor * factor, 2e-10);
        assert!(wall.pressure > 0.0 && wall.pressure < p.pressure);
        let f = gas.wall_flux(state, p, -1.0, &mut || false).unwrap();
        assert_eq!(f.mass, 0.0);
        assert_eq!(f.energy, 0.0);
        assert_eq!(f.momentum, wall.pressure);
    }
}

#[test]
fn variable_cp_wall_shock_closes_independent_mass_momentum_energy_and_entropy() {
    let gas = gas(0.4, -1e18);
    let q = 3.0;
    let state = gas.conserved(0.8, q, 1.6).unwrap();
    let before = gas.primitive(state, || false).unwrap();
    let wall = gas.stopped_wall(before, q, &mut || false).unwrap();
    let w = q / (wall.density / before.density - 1.0);
    let j1 = before.density * (q + w);
    let j2 = wall.density * w;
    close(j1, j2, 2e-11);
    close(
        before.pressure + j1 * (q + w),
        wall.pressure + j2 * w,
        2e-10,
    );
    // Cp/R = 3.5+0.4T, so this independently integrates the declared Cv/R.
    let delta_e = 2.5 * (wall.temperature - before.temperature)
        + 0.2 * (wall.temperature.powi(2) - before.temperature.powi(2));
    let hugoniot =
        0.5 * (wall.pressure + before.pressure) * (1.0 / before.density - 1.0 / wall.density);
    close(delta_e, hugoniot, 2e-10);
    let entropy = 3.5 * (wall.temperature / before.temperature).ln()
        + 0.4 * (wall.temperature - before.temperature)
        - (wall.pressure / before.pressure).ln();
    assert!(entropy > 0.0);
    let frozen_gamma = (3.5 + 0.4 * before.temperature) / (2.5 + 0.4 * before.temperature);
    let k = 0.5 * (frozen_gamma + 1.0) * before.density * q * q;
    let b = (frozen_gamma - 1.0) / (frozen_gamma + 1.0) * before.pressure;
    let frozen_pressure =
        before.pressure + 0.5 * k + 0.5 * (k * k + 4.0 * k * (before.pressure + b)).sqrt();
    assert!((wall.pressure / frozen_pressure - 1.0).abs() > 0.005);
}

#[test]
fn variable_cp_rarefaction_matches_independent_simpson_integral_and_entropy() {
    let gas = gas(0.4, 1e18);
    let state = gas.conserved(1.0, -1.0, 4.0).unwrap();
    let before = gas.primitive(state, || false).unwrap();
    let wall = gas.stopped_wall(before, -1.0, &mut || false).unwrap();
    let n = 2_000;
    let h = (before.temperature - wall.temperature) / n as f64;
    let f = |t: f64| ((3.5 + 0.4 * t) * (2.5 + 0.4 * t) / t).sqrt();
    let mut sum = f(wall.temperature) + f(before.temperature);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(wall.temperature + i as f64 * h);
    }
    close(h * sum / 3.0, 1.0, 2e-10);
    let entropy = 3.5 * (wall.temperature / before.temperature).ln()
        + 0.4 * (wall.temperature - before.temperature)
        - (wall.pressure / before.pressure).ln();
    assert!(entropy.abs() < 1e-10);
}

#[test]
fn wall_is_orientation_and_formation_reference_invariant() {
    let reference = gas(0.4, 0.0);
    for velocity in [-1.0_f64, 0.0, 2.0] {
        let state = reference.conserved(1.0, velocity, 3.0).unwrap();
        let p = reference.primitive(state, || false).unwrap();
        let expected = reference.wall_flux(state, p, 1.0, &mut || false).unwrap();
        let reverse = reference.conserved(1.0, -velocity, 3.0).unwrap();
        let reversed = reference
            .wall_flux(
                reverse,
                reference.primitive(reverse, || false).unwrap(),
                -1.0,
                &mut || false,
            )
            .unwrap();
        assert_eq!(expected, reversed);
        for formation in [-1e18, 1e18] {
            let shifted = gas(0.4, formation);
            let q = shifted.conserved(1.0, velocity, 3.0).unwrap();
            let f = shifted
                .wall_flux(q, shifted.primitive(q, || false).unwrap(), 1.0, &mut || {
                    false
                })
                .unwrap();
            assert_eq!(f, expected);
        }
    }
}

#[test]
fn wall_refuses_unsupported_hot_cold_and_cancelled_paths_without_extrapolation() {
    let gas = gas(0.0, 0.0);
    for (velocity, pressure, normal) in [(30.0, 1.0, 1.0), (4.0, 0.4, -1.0)] {
        let state = gas.conserved(1.0, velocity, pressure).unwrap();
        let primitive = gas.primitive(state, || false).unwrap();
        assert!(matches!(
            gas.wall_flux(state, primitive, normal, &mut || false),
            Err(Error::CaloricPath(_))
        ));
    }
    let state = gas.conserved(1.0, -1.0, 3.0).unwrap();
    let p = gas.primitive(state, || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        gas.wall_flux(state, p, 1.0, &mut || {
            polls += 1;
            polls >= 4
        }),
        Err(Error::Cancelled)
    );
}

#[test]
fn moving_piston_faces_match_analytic_shock_rarefaction_and_work() {
    let gas = gas(0.0, 0.0);
    let density = 0.8;
    let pressure = 1.6;
    let gamma = 1.4_f64;
    let sound = (gamma * pressure / density).sqrt();
    for speed in [-3.0_f64, -0.2, 0.0, 0.2, 1.0] {
        // Positive speed expands the domain; negative speed compresses it.
        let q = -speed;
        let exact = if q > 0.0 {
            let k = 0.5 * (gamma + 1.0) * density * q * q;
            let b = (gamma - 1.0) / (gamma + 1.0) * pressure;
            pressure + 0.5 * k + 0.5 * (k * k + 4.0 * k * (pressure + b)).sqrt()
        } else {
            pressure * (1.0 + (gamma - 1.0) * q / (2.0 * sound)).powf(7.0)
        };
        for normal in [-1.0, 1.0] {
            let wall_velocity = normal * speed;
            let input = NormalPrimitive {
                density,
                normal_velocity: 0.0,
                transverse_velocity: 2.0,
                pressure,
            };
            let face = gas
                .moving_wall_normal_flux_from_primitive(input, normal, wall_velocity, || false)
                .unwrap();
            close(face.normal_momentum, exact, 3e-10);
            assert_eq!(face.mass, 0.0);
            assert_eq!(face.transverse_momentum, 0.0);
            assert_eq!(face.energy, face.normal_momentum * wall_velocity);
            assert!(face.signal_speed + 8.0 * f64::EPSILON * sound >= sound + q.abs());
            assert!(!face.used_hll_fallback);
            let state = gas.normal_conserved(density, 0.0, 2.0, pressure).unwrap();
            let conservative = gas
                .moving_wall_normal_flux(state, normal, wall_velocity, || false)
                .unwrap();
            close(conservative.normal_momentum, face.normal_momentum, 3e-10);
            close(conservative.energy, face.energy, 3e-10);
            let spatial_input = SpatialPrimitive {
                density,
                normal_velocity: 0.0,
                transverse_velocity: [2.0, -3.0],
                pressure,
            };
            let spatial = gas
                .moving_wall_spatial_flux_from_primitive(
                    spatial_input,
                    normal,
                    wall_velocity,
                    || false,
                )
                .unwrap();
            assert_eq!(spatial.normal_momentum, face.normal_momentum);
            assert_eq!(spatial.energy, face.energy);
            assert_eq!(spatial.signal_speed, face.signal_speed);
            assert_eq!(spatial.transverse_momentum, [0.0; 2]);
            assert_eq!(spatial.mass, 0.0);
            let state3 = gas
                .spatial_conserved(density, 0.0, [2.0, -3.0], pressure)
                .unwrap();
            let conservative3 = gas
                .moving_wall_spatial_flux(state3, normal, wall_velocity, || false)
                .unwrap();
            close(
                conservative3.normal_momentum,
                spatial.normal_momentum,
                3e-10,
            );
            close(conservative3.energy, spatial.energy, 3e-10);
            if speed == 0.0 {
                assert_eq!(
                    face,
                    gas.stationary_wall_normal_flux_from_primitive(input, normal, || false)
                        .unwrap()
                );
                assert_eq!(
                    conservative,
                    gas.stationary_wall_normal_flux(state, normal, || false)
                        .unwrap()
                );
                assert_eq!(
                    spatial,
                    gas.stationary_wall_spatial_flux_from_primitive(spatial_input, normal, || {
                        false
                    })
                    .unwrap()
                );
                assert_eq!(
                    conservative3,
                    gas.stationary_wall_spatial_flux(state3, normal, || false)
                        .unwrap()
                );
            }
        }
    }
}

#[test]
fn moving_piston_work_equals_variable_caloric_shocked_slab_energy() {
    for reference in [-1e18, 0.0, 1e18] {
        let gas = gas(0.4, reference);
        let (rho1, t1, speed) = (0.8_f64, 2.0_f64, 3.0_f64);
        let p1 = rho1 * t1;
        let flux = gas
            .moving_wall_normal_flux_from_primitive(
                NormalPrimitive {
                    density: rho1,
                    normal_velocity: 0.0,
                    transverse_velocity: 7.0,
                    pressure: p1,
                },
                1.0,
                -speed,
                || false,
            )
            .unwrap();
        let p2 = flux.normal_momentum;
        let w = (p2 - p1) / (rho1 * speed) - speed;
        let rho2 = rho1 * (speed + w) / w;
        let t2 = p2 / rho2;
        // Independently integrate Cv(T)=2.5+0.4T and the piecewise-uniform
        // piston/shock solution: shocked length grows at w, upstream length
        // shrinks at speed+w. The unchanged transverse kinetic energy and
        // arbitrary formation reference cancel with the conserved mass.
        let delta_e = 2.5 * (t2 - t1) + 0.2 * (t2 * t2 - t1 * t1);
        let time = 0.01;
        let swept_mass = rho1 * (speed + w) * time;
        let energy_gain = swept_mass * (delta_e + 0.5 * speed * speed);
        close(rho2 * w * time, swept_mass, 3e-10);
        close(-rho2 * speed * w, p1 - p2, 3e-10);
        close(energy_gain, -flux.energy * time, 3e-10);
        assert!(energy_gain > 0.0);
        let entropy = 3.5 * (t2 / t1).ln() + 0.4 * (t2 - t1) - (p2 / p1).ln();
        assert!(entropy > 0.0);
    }
}

#[test]
fn moving_wall_flux_has_galilean_pressure_and_relative_wave_speeds() {
    let gas = gas(0.4, 1e18);
    let input = SpatialPrimitive {
        density: 1.0,
        normal_velocity: -0.25,
        transverse_velocity: [0.75, -1.5],
        pressure: 3.0,
    };
    for normal in [-1.0, 1.0] {
        let base = gas
            .moving_wall_spatial_flux_from_primitive(input, normal, 0.5, || false)
            .unwrap();
        for boost in [-8.0, 2.5, 30.0] {
            let translated = SpatialPrimitive {
                normal_velocity: input.normal_velocity + boost,
                transverse_velocity: [
                    input.transverse_velocity[0] + 3.0,
                    input.transverse_velocity[1] - 4.0,
                ],
                ..input
            };
            let flux = gas
                .moving_wall_spatial_flux_from_primitive(translated, normal, 0.5 + boost, || false)
                .unwrap();
            assert_eq!(flux.normal_momentum, base.normal_momentum);
            assert_eq!(flux.signal_speed, base.signal_speed);
            close(
                flux.energy,
                base.energy + boost * base.normal_momentum,
                2e-14,
            );
            assert_eq!(flux.mass, 0.0);
            assert_eq!(flux.transverse_momentum, [0.0; 2]);
        }
    }
}

#[test]
fn moving_wall_refuses_unsupported_motion_and_cancellation() {
    let gas = gas(0.0, 0.0);
    let input = NormalPrimitive {
        density: 1.0,
        normal_velocity: 0.0,
        transverse_velocity: 0.0,
        pressure: 0.4,
    };
    for speed in [-30.0, 4.0] {
        assert!(matches!(
            gas.moving_wall_normal_flux_from_primitive(input, 1.0, speed, || false),
            Err(Error::CaloricPath(_))
        ));
    }
    for speed in [f64::NAN, f64::INFINITY] {
        assert!(
            gas.moving_wall_normal_flux_from_primitive(input, 1.0, speed, || false)
                .is_err()
        );
    }
    let mut count = 0;
    gas.moving_wall_normal_flux_from_primitive(input, 1.0, -1.0, || {
        count += 1;
        false
    })
    .unwrap();
    for stop in [1, count / 2, count] {
        let mut polls = 0;
        assert_eq!(
            gas.moving_wall_normal_flux_from_primitive(input, 1.0, -1.0, || {
                polls += 1;
                polls == stop
            }),
            Err(Error::Cancelled)
        );
    }
}
