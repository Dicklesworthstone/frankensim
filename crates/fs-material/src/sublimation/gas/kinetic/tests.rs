use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, tests as fixture};

fn frame() -> MolecularSurfaceFrame {
    MolecularSurfaceFrame {
        surface_id: 1,
        outward_normal: [0.0, 0.0, 1.0],
        centroid_m: [0.2, 0.3, 0.4],
    }
}
fn reflection(catalog: &ThermochemicalGasModel, chi: f64) -> MaxwellReflection {
    MaxwellReflection::from_diffuse_fractions(
        catalog,
        &[
            (SpeciesId::new("A").unwrap(), chi),
            (fixture::species(), chi),
        ],
        [200.0, 900.0],
    )
    .unwrap()
}
fn stationary() -> MolecularInterfaceVelocity {
    MolecularInterfaceVelocity {
        gas_velocity_m_s: [0.0; 3],
        surface_velocity_m_s: [0.0; 3],
    }
}
#[allow(clippy::too_many_arguments)]
fn rate(
    reference: f64,
    ts: f64,
    tg: f64,
    pv: f64,
    beta: f64,
    chi: f64,
    velocity: MolecularInterfaceVelocity,
) -> KineticInterfaceRate {
    let catalog = fixture::catalog(reference);
    let binding = GasSurfaceBinding {
        accommodation_coefficient: beta,
        ..fixture::binding(reference)
    };
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let surface = GasChamberSurface {
        id: 1,
        law: &phase,
        binding: &binding,
        solid: fixture::solid(ts, chart),
        temperature_from_enthalpy: &temperature,
    };
    evaluate_molecular_interface(
        &surface,
        None,
        frame(),
        [0.0; 3],
        velocity,
        Some(&reflection(&catalog, chi)),
        &fixture::gas(&catalog, tg, pv),
        &catalog,
        fixture::controls().latent_compatibility_tolerance_j_kg,
        || false,
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * b.abs().max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn stationary_limit_matches_original_phase_beams_pressure_and_reflected_heat() {
    let catalog = fixture::catalog(0.0);
    for (ts, tg, pv) in [
        (400.0, 400.0, 1000.0),
        (350.0, 500.0, 4000.0),
        (700.0, 250.0, 0.0),
    ] {
        let chi = 0.4;
        let result = rate(0.0, ts, tg, pv, 0.1, chi, stationary());
        let original = molecular_mass_flux(
            &fixture::phase(),
            ts,
            0.1,
            VaporReservoir::Prescribed {
                partial_pressure_pa: pv,
                temperature_k: tg,
            },
        )
        .unwrap();
        for row in &result.species {
            let volatile = row.species == fixture::species();
            let p = if volatile { pv } else { 100000.0 - pv };
            let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / if volatile { 0.018015 } else { 0.02897 };
            let beta = if volatile { 0.1 } else { 0.0 };
            let je = if volatile {
                original.outgoing_kg_m2_s
            } else {
                0.0
            };
            let ji = p / (2.0 * core::f64::consts::PI * r * tg).sqrt();
            let cs = (core::f64::consts::PI * r * ts / 2.0).sqrt();
            let pressure =
                je * cs + 0.5 * p * (beta + (1.0 - beta) * (2.0 - chi + chi * (ts / tg).sqrt()));
            near(row.incident_mass_flux_kg_m2_s, ji);
            near(row.sticking_mass_flux_kg_m2_s, beta * ji);
            near(row.emitted_mass_flux_kg_m2_s, je);
            near(row.relative_momentum_flux_into_gas_pa[2], pressure);
            assert_eq!(&row.relative_momentum_flux_into_gas_pa[..2], &[0.0, 0.0]);
            let hs = catalog
                .species_specific_enthalpy_j_kg(&row.species, ts)
                .unwrap();
            let hg = catalog
                .species_specific_enthalpy_j_kg(&row.species, tg)
                .unwrap();
            near(
                row.laboratory_energy_flux_into_gas_w_m2,
                je * (hs - 0.5 * r * ts) - beta * ji * (hg - 0.5 * r * tg)
                    + (1.0 - beta) * chi * ji * ((hs - hg) - 0.5 * r * (ts - tg)),
            );
            assert_eq!(row.mechanical_power_into_gas_w_m2, 0.0);
        }
        near(
            result.mass_rate_into_gas_kg_s,
            original.net_kg_m2_s * result.exposed_solid_area_m2,
        );
        near(
            result.torque_on_condensed_material_nm[0],
            frame().centroid_m[1] * result.force_on_condensed_material_n[2],
        );
        near(
            result.torque_on_condensed_material_nm[1],
            -frame().centroid_m[0] * result.force_on_condensed_material_n[2],
        );
    }
}

// Independent finite-range Simpson quadrature of the shifted Gaussian. No
// production erfc, continued fraction or molecular-moment helper is called.
fn gaussian_integrals(a: f64) -> [f64; 3] {
    let count = 12000;
    let end = a.max(0.0) + 12.0;
    let h = end / count as f64;
    let mut integrals = [0.0; 3];
    for i in 0..=count {
        let y = h * i as f64;
        let weight = if i == 0 || i == count {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let density = (-0.5 * (y - a) * (y - a)).exp() / (2.0 * core::f64::consts::PI).sqrt();
        for (k, value) in integrals.iter_mut().enumerate() {
            *value += weight * y.powi(k as i32 + 1) * density;
        }
    }
    integrals.map(|value| value * h / 3.0)
}

#[test]
fn drifting_mixture_matches_independent_incident_mass_momentum_and_energy_integrals() {
    let catalog = fixture::catalog(0.0);
    for vn in [-600.0, 200.0, 800.0] {
        let ts = 450.0;
        let tg = 400.0;
        let chi = 0.6;
        let beta = 0.3;
        let velocity = MolecularInterfaceVelocity {
            gas_velocity_m_s: [30.0, -20.0, vn],
            ..stationary()
        };
        let result = rate(0.0, ts, tg, 2000.0, beta, chi, velocity);
        for row in &result.species {
            let volatile = row.species == fixture::species();
            let p = if volatile { 2000.0 } else { 98000.0 };
            let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / if volatile { 0.018015 } else { 0.02897 };
            let beta = if volatile { beta } else { 0.0 };
            let sigma = (r * tg).sqrt();
            let rho = p / (r * tg);
            let [m1, m2, m3] = gaussian_integrals(-vn / sigma);
            let ji = rho * sigma * m1;
            let normal_momentum = rho * sigma * sigma * m2;
            let hs = catalog
                .species_specific_enthalpy_j_kg(&row.species, ts)
                .unwrap();
            let hg = catalog
                .species_specific_enthalpy_j_kg(&row.species, tg)
                .unwrap();
            let incident_energy =
                rho * sigma * ((hg - 1.5 * r * tg + 0.5 * 1300.0) * m1 + 0.5 * r * tg * m3);
            let eout = hs - 0.5 * r * ts;
            let je = row.emitted_mass_flux_kg_m2_s;
            let cs = (core::f64::consts::PI * r * ts / 2.0).sqrt();
            let normal = je * cs
                + beta * normal_momentum
                + (1.0 - beta)
                    * (2.0 * (1.0 - chi) * normal_momentum + chi * (ji * cs + normal_momentum));
            near(row.incident_mass_flux_kg_m2_s, ji);
            near(row.relative_momentum_flux_into_gas_pa[2], normal);
            for (d, tangent) in [30.0, -20.0].into_iter().enumerate() {
                near(
                    row.relative_momentum_flux_into_gas_pa[d],
                    -(beta + (1.0 - beta) * chi) * ji * tangent,
                );
            }
            near(
                row.laboratory_energy_flux_into_gas_w_m2,
                je * eout - beta * incident_energy
                    + (1.0 - beta) * chi * (ji * eout - incident_energy),
            );
        }
    }
}

#[test]
fn specular_slip_has_no_drag_while_diffuse_slip_has_mechanical_work_and_crossing_energy() {
    let velocity = MolecularInterfaceVelocity {
        gas_velocity_m_s: [40.0, -20.0, 0.0],
        surface_velocity_m_s: [10.0, 5.0, 0.0],
    };
    let specular = rate(0.0, 400.0, 400.0, 1000.0, 0.0, 0.0, velocity);
    for (actual, expected) in specular
        .force_on_condensed_material_n
        .into_iter()
        .zip([0.0, 0.0, -10.0])
    {
        near(actual, expected);
    }
    assert_eq!(specular.laboratory_energy_rate_into_gas_w, 0.0);
    let diffuse = rate(0.0, 400.0, 400.0, 1000.0, 0.0, 1.0, velocity);
    for row in diffuse.species {
        let j = row.incident_mass_flux_kg_m2_s;
        near(row.relative_momentum_flux_into_gas_pa[0], -j * 30.0);
        near(row.relative_momentum_flux_into_gas_pa[1], j * 25.0);
        near(
            row.relative_energy_flux_beyond_gas_enthalpy_w_m2,
            -0.5 * j * (30.0 * 30.0 + 25.0 * 25.0),
        );
        near(
            row.mechanical_power_into_gas_w_m2,
            -j * (10.0 * 30.0 - 5.0 * 25.0),
        );
        assert_eq!(row.transported_surface_kinetic_energy_flux_w_m2, 0.0);
    }
}

#[test]
fn common_inertial_boost_preserves_force_and_transforms_transport_momentum_and_energy() {
    let velocity = MolecularInterfaceVelocity {
        gas_velocity_m_s: [32.0, -16.0, -64.0],
        surface_velocity_m_s: [8.0, 4.0, 16.0],
    };
    let boost = [100.0, -40.0, 8.0];
    let boosted = MolecularInterfaceVelocity {
        gas_velocity_m_s: core::array::from_fn(|d| velocity.gas_velocity_m_s[d] + boost[d]),
        surface_velocity_m_s: core::array::from_fn(|d| velocity.surface_velocity_m_s[d] + boost[d]),
    };
    let old = rate(0.0, 400.0, 500.0, 3000.0, 0.2, 0.7, velocity);
    let new = rate(0.0, 400.0, 500.0, 3000.0, 0.2, 0.7, boosted);
    assert_eq!(
        old.force_on_condensed_material_n,
        new.force_on_condensed_material_n
    );
    assert_eq!(
        old.torque_on_condensed_material_nm,
        new.torque_on_condensed_material_nm
    );
    for (a, b) in old.species.iter().zip(&new.species) {
        assert_eq!(
            a.relative_energy_flux_beyond_gas_enthalpy_w_m2,
            b.relative_energy_flux_beyond_gas_enthalpy_w_m2
        );
        let j = a.net_mass_flux_into_gas_kg_m2_s;
        for d in 0..3 {
            near(
                b.laboratory_momentum_flux_into_gas_pa[d],
                a.laboratory_momentum_flux_into_gas_pa[d] + boost[d] * j,
            );
        }
        let expected = a.laboratory_energy_flux_into_gas_w_m2
            + dot(boost, a.laboratory_momentum_flux_into_gas_pa).unwrap()
            + 0.5 * dot(boost, boost).unwrap() * j;
        near(b.laboratory_energy_flux_into_gas_w_m2, expected);
    }
}

#[test]
fn formation_reference_does_not_change_source_relative_energy_or_motion() {
    let velocity = MolecularInterfaceVelocity {
        gas_velocity_m_s: [10.0, 4.0, -30.0],
        surface_velocity_m_s: [5.0, -2.0, 8.0],
    };
    let old = rate(0.0, 400.0, 450.0, 2000.0, 0.1, 0.4, velocity);
    let shifted = rate(-1e10, 400.0, 450.0, 2000.0, 0.1, 0.4, velocity);
    assert_eq!(
        old.force_on_condensed_material_n,
        shifted.force_on_condensed_material_n
    );
    for (a, b) in old.species.iter().zip(&shifted.species) {
        near(
            a.relative_energy_flux_beyond_gas_enthalpy_w_m2,
            b.relative_energy_flux_beyond_gas_enthalpy_w_m2,
        );
        near(
            a.laboratory_energy_flux_beyond_gas_enthalpy_w_m2,
            b.laboratory_energy_flux_beyond_gas_enthalpy_w_m2,
        );
        assert_eq!(
            a.mechanical_power_into_gas_w_m2,
            b.mechanical_power_into_gas_w_m2
        );
    }
}

#[test]
fn sub_ulp_normal_drift_at_saturation_keeps_a_signed_phase_rate() {
    let zero = rate(0.0, 400.0, 400.0, 1000.0, 0.1, 0.0, stationary());
    let moving = rate(
        0.0,
        400.0,
        400.0,
        1000.0,
        0.1,
        0.0,
        MolecularInterfaceVelocity {
            gas_velocity_m_s: [0.0, 0.0, -1e-15],
            ..stationary()
        },
    );
    let a = zero
        .species
        .iter()
        .find(|row| row.species == fixture::species())
        .unwrap();
    let b = moving
        .species
        .iter()
        .find(|row| row.species == fixture::species())
        .unwrap();
    assert_eq!(a.net_mass_flux_into_gas_kg_m2_s, 0.0);
    assert_eq!(a.incident_mass_flux_kg_m2_s, b.incident_mass_flux_kg_m2_s);
    assert!(b.net_mass_flux_into_gas_kg_m2_s < 0.0);
    let expected =
        -0.1 * (1000.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015 * 400.0)) * 0.5e-15;
    assert!((b.net_mass_flux_into_gas_kg_m2_s / expected - 1.0).abs() < 1e-13);
}

#[test]
fn sharp_fusion_uses_solid_constituent_enthalpy_and_declared_exposed_area() {
    use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (hs, cp) = fixture::chart(&catalog, &binding);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        fs_blake3::ContentHash([91; 32]),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: hs - 200.0 * cp,
                temperature_k: 200.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 917.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: hs,
                temperature_k: 400.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 917.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: hs + 80000.0,
                temperature_k: 400.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    let source = SharpMeltingCaloric::new(&curve).unwrap();
    let temperature = |h| {
        curve.temperature_at_specific_enthalpy(h).map_err(|_| {
            fs_thermal::sublimation::SublimationError::Invalid("fixture caloric domain")
        })
    };
    let phase = fixture::phase();
    let surface = GasChamberSurface {
        id: 1,
        law: &phase,
        binding: &binding,
        solid: fs_thermal::sublimation::SolidSurfaceState {
            specific_enthalpy_j_kg: hs + 40000.0,
            ..fixture::solid(400.0, (hs, cp))
        },
        temperature_from_enthalpy: &temperature,
    };
    let gas = fixture::gas(&catalog, 400.0, 2000.0);
    let evaluate = |s: &GasChamberSurface<'_>, c| {
        evaluate_molecular_interface(
            s,
            c,
            frame(),
            [0.0; 3],
            stationary(),
            None,
            &gas,
            &catalog,
            1e-5,
            || false,
        )
    };
    assert!(evaluate(&surface, None).is_err());
    let mixed = evaluate(&surface, Some(source)).unwrap();
    let pure = evaluate(
        &GasChamberSurface {
            solid: fixture::solid(400.0, (hs, cp)),
            ..surface
        },
        None,
    )
    .unwrap();
    assert_eq!(mixed.species, pure.species);
    assert_eq!(
        mixed.force_on_condensed_material_n,
        pure.force_on_condensed_material_n
    );
    assert_eq!(mixed.exposed_solid_area_m2, surface.solid.area_m2);
    assert!(
        evaluate(
            &GasChamberSurface {
                solid: fs_thermal::sublimation::SolidSurfaceState {
                    specific_enthalpy_j_kg: hs + 80000.0,
                    ..surface.solid
                },
                ..surface
            },
            Some(source)
        )
        .is_err()
    );
}

#[test]
fn wrong_sources_frames_and_late_cancellation_leave_inputs_unchanged() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let temperature = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let surface = GasChamberSurface {
        id: 1,
        law: &phase,
        binding: &binding,
        solid: fixture::solid(400.0, chart),
        temperature_from_enthalpy: &temperature,
    };
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let saved = gas.clone();
    let evaluate = |f, c| {
        evaluate_molecular_interface(
            &surface,
            None,
            f,
            [0.0; 3],
            stationary(),
            None,
            &gas,
            c,
            1e-5,
            || false,
        )
    };
    assert!(
        evaluate(
            MolecularSurfaceFrame {
                outward_normal: [0.0, 0.0, 2.0],
                ..frame()
            },
            &catalog
        )
        .is_err()
    );
    assert!(
        evaluate(
            MolecularSurfaceFrame {
                surface_id: 2,
                ..frame()
            },
            &catalog
        )
        .is_err()
    );
    assert!(evaluate(frame(), &fixture::catalog(1e6)).is_err());
    let mut polls = 0;
    let result = evaluate_molecular_interface(
        &surface,
        None,
        frame(),
        [0.0; 3],
        stationary(),
        None,
        &gas,
        &catalog,
        1e-5,
        || {
            polls += 1;
            polls == 6
        },
    );
    assert!(result.is_err());
    assert_eq!(polls, 6);
    assert_eq!(gas, saved);
    assert_eq!(surface.solid, fixture::solid(400.0, chart));
}
