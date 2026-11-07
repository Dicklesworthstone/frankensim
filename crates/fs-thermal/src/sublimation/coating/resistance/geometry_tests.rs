//! G0/G1/G3 geometry, Fourier-circuit and original-physics handoff evidence.
use super::*;
use crate::sublimation::coating::resistance::{
    ResistiveSurfaceState, advance_resistive_surface_with_heat, heat_resistive_surface,
};
use crate::sublimation::{
    MOLAR_GAS_CONSTANT, SolidSurfaceState, SolverConfig, SublimationLaw, SublimationProperties,
    SurfaceKinetics, VaporReservoir,
};

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:e}, expected={expected:e}, relative={relative:e}"
    );
}

#[test]
fn planar_mass_controls_thickness_and_original_contact_series_resistance() {
    let layer = ConductiveCoatingGeometry::planar(1000.0, 2.0).unwrap();
    let result = layer.at_mass(1.0, 0.2, 0.01).unwrap();
    close(result.front.extent_m, 0.0005, 2e-15);
    assert_eq!(result.front.mass_kg, 1.0);
    assert_eq!(result.front.area_m2, 2.0);
    assert_eq!(result.contact_area_m2, 2.0);
    close(result.layer_resistance_k_w, 0.00125, 2e-15);
    close(result.contact_resistance_k_w, 0.005, 2e-15);
    close(result.conductance_w_k, 160.0, 2e-15);
    assert_eq!(layer.geometry().at_mass(1.0).unwrap(), result.front);
    assert_eq!(layer.geometry().occupied_volume_m3(1.0).unwrap(), 0.001);
    assert_eq!(layer.contact_area_m2(), result.contact_area_m2);
}

#[test]
fn thick_curved_shells_match_independent_fourier_solutions() {
    let density: f64 = 880.0;
    let inner: f64 = 0.2;
    let outer: f64 = 0.7;
    let length: f64 = 1.3;
    let k = 0.18;
    let contact = 0.04;
    // Deliberately use direct volume differences for these well-separated
    // radii; the production inversion uses factored shell identities.
    let spherical_mass = (4.0 * PI / 3.0) * (outer.powi(3) - inner.powi(3)) * density;
    let sphere = ConductiveCoatingGeometry::spherical(density, inner).unwrap();
    let spherical = sphere.at_mass(spherical_mass, k, contact).unwrap();
    let sphere_inner_area = 4.0 * PI * inner * inner;
    let spherical_resistance = (1.0 / inner - 1.0 / outer) / (4.0 * PI * k);
    close(spherical.front.extent_m, outer - inner, 2e-14);
    close(spherical.front.area_m2, 4.0 * PI * outer * outer, 2e-14);
    close(spherical.layer_resistance_k_w, spherical_resistance, 2e-14);
    close(
        spherical.contact_resistance_k_w,
        contact / sphere_inner_area,
        2e-14,
    );
    close(
        spherical.conductance_w_k,
        1.0 / (spherical_resistance + contact / sphere_inner_area),
        2e-14,
    );
    close(
        sphere
            .geometry()
            .occupied_volume_m3(spherical_mass)
            .unwrap(),
        (4.0 * PI / 3.0) * outer.powi(3),
        2e-14,
    );

    let cylindrical_mass = PI * length * (outer * outer - inner * inner) * density;
    let cylinder = ConductiveCoatingGeometry::cylindrical(density, inner, length).unwrap();
    let cylindrical = cylinder.at_mass(cylindrical_mass, k, contact).unwrap();
    let cylinder_inner_area = 2.0 * PI * inner * length;
    let cylindrical_resistance = (outer / inner).ln() / (2.0 * PI * k * length);
    close(cylindrical.front.extent_m, outer - inner, 2e-14);
    close(cylindrical.front.area_m2, 2.0 * PI * outer * length, 2e-14);
    close(
        cylindrical.layer_resistance_k_w,
        cylindrical_resistance,
        2e-14,
    );
    close(
        cylindrical.contact_resistance_k_w,
        contact / cylinder_inner_area,
        2e-14,
    );
    close(
        cylindrical.conductance_w_k,
        1.0 / (cylindrical_resistance + contact / cylinder_inner_area),
        2e-14,
    );
    close(
        cylinder
            .geometry()
            .occupied_volume_m3(cylindrical_mass)
            .unwrap(),
        PI * length * outer * outer,
        2e-14,
    );

    assert!(spherical.front.area_m2 > spherical.contact_area_m2);
    assert!(cylindrical.front.area_m2 > cylindrical.contact_area_m2);
}

#[test]
fn sub_ulp_shells_keep_positive_thickness_resistance_and_planar_limit() {
    let density = 1000.0;
    let radius = 1.0;
    let thickness = 1e-20;
    let conductivity = 0.25;
    assert_eq!(radius + thickness, radius);
    for geometry in [
        ConductiveCoatingGeometry::spherical(density, radius).unwrap(),
        ConductiveCoatingGeometry::cylindrical(density, radius, 0.3).unwrap(),
    ] {
        let area = geometry.contact_area_m2();
        let mass = density * area * thickness;
        let result = geometry.at_mass(mass, conductivity, 0.0).unwrap();
        close(result.front.extent_m, thickness, 3e-14);
        close(result.front.area_m2, area, 3e-14);
        close(
            result.layer_resistance_k_w,
            thickness / (conductivity * area),
            3e-14,
        );
        assert!(result.layer_resistance_k_w > 0.0);
        assert!(result.conductance_w_k.is_finite());
        let coated = geometry.at_mass(4.0 * mass, conductivity, 0.0).unwrap();
        close(
            coated.layer_resistance_k_w,
            4.0 * result.layer_resistance_k_w,
            3e-14,
        );
    }
}

#[test]
fn radial_series_subdivision_preserves_total_conduction_and_contact_area() {
    let density = 900.0;
    let inner: f64 = 0.01;
    let middle: f64 = 0.014;
    let outer: f64 = 0.03;
    let length = 0.4;
    let k = 0.6;
    let contact = 0.02;
    for spherical in [false, true] {
        let geometry = |radius| {
            if spherical {
                ConductiveCoatingGeometry::spherical(density, radius)
            } else {
                ConductiveCoatingGeometry::cylindrical(density, radius, length)
            }
            .unwrap()
        };
        let mass = |a: f64, b: f64| {
            density
                * if spherical {
                    (4.0 * PI / 3.0) * (b.powi(3) - a.powi(3))
                } else {
                    PI * length * (b * b - a * a)
                }
        };
        let whole = geometry(inner)
            .at_mass(mass(inner, outer), k, contact)
            .unwrap();
        let first = geometry(inner)
            .at_mass(mass(inner, middle), k, contact)
            .unwrap();
        let second = geometry(middle)
            .at_mass(mass(middle, outer), k, 0.0)
            .unwrap();
        close(
            whole.layer_resistance_k_w,
            first.layer_resistance_k_w + second.layer_resistance_k_w,
            3e-14,
        );
        close(
            whole.contact_resistance_k_w,
            first.contact_resistance_k_w,
            3e-14,
        );
        close(
            whole.conductance_w_k.recip(),
            first.conductance_w_k.recip() + second.conductance_w_k.recip(),
            3e-14,
        );
        assert!(whole.front.area_m2 > first.front.area_m2);
        assert!(whole.layer_resistance_k_w > first.layer_resistance_k_w);
    }
}

struct Law;
impl SublimationLaw for Law {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 400.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let pressure = 10.0 * (t / 230.0).sqrt();
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure / (2.0 * t),
            latent_heat_j_kg: 2.8e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}

fn controls() -> SolverConfig {
    SolverConfig {
        relative_tolerance: 1e-12,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 256,
    }
}

#[test]
fn mass_owned_circuit_drives_original_two_body_caloric_solver() {
    let geometry = ConductiveCoatingGeometry::planar(900.0, 0.01).unwrap();
    let initial = ResistiveSurfaceState {
        substrate_mass_kg: 1.0,
        substrate_specific_enthalpy_j_kg: 1000.0 * 330.0,
        deposit: SolidSurfaceState {
            mass_kg: 0.01,
            specific_enthalpy_j_kg: 2000.0 * 270.0,
            area_m2: 0.01,
        },
    };
    let circuit = geometry
        .at_mass(initial.deposit.mass_kg, 0.2, 0.002)
        .unwrap();
    let dt = 2.0;
    let step = heat_resistive_surface(
        &Law,
        initial,
        dt,
        circuit.conductance_w_k,
        0.0,
        0.0,
        controls(),
        |h| Ok(h / 1000.0),
        |h| Ok(h / 2000.0),
        || Ok(()),
    )
    .unwrap();
    // Independent solution of the 2x2 backward-Euler system.
    let support_capacity = initial.substrate_mass_kg * 1000.0;
    let deposit_capacity = initial.deposit.mass_kg * 2000.0;
    let contact = dt * circuit.conductance_w_k;
    let q = contact * (330.0 - 270.0)
        / (1.0 + contact * (1.0 / support_capacity + 1.0 / deposit_capacity));
    close(step.interface_heat_j, q, 1e-10);
    close(
        step.substrate_temperature_k,
        330.0 - q / support_capacity,
        1e-12,
    );
    close(
        step.deposit_temperature_k,
        270.0 + q / deposit_capacity,
        1e-12,
    );
    assert!(step.energy_residual_j.abs() < 1e-9);
    assert_eq!(step.state.deposit.mass_kg, initial.deposit.mass_kg);
}

#[test]
fn original_deposition_transaction_refreshes_actual_next_layer_resistance() {
    let geometry = ConductiveCoatingGeometry::planar(900.0, 0.001).unwrap();
    let mut state = ResistiveSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: 500.0 * 230.0,
        deposit: SolidSurfaceState {
            mass_kg: 0.001,
            specific_enthalpy_j_kg: 2000.0 * 230.0,
            area_m2: 0.001,
        },
    };
    let dt = 0.02;
    let flux = (10.0 - 30.0) * (0.018 / (2.0 * PI * MOLAR_GAS_CONSTANT * 230.0)).sqrt();
    for _ in 0..2 {
        let before = geometry.at_mass(state.deposit.mass_kg, 0.1, 0.003).unwrap();
        let result = advance_resistive_surface_with_heat(
            &Law,
            SurfaceKinetics {
                accommodation_coefficient: 1.0,
                mass_transfer_coefficient_m_s: None,
            },
            state,
            VaporReservoir::Prescribed {
                partial_pressure_pa: 30.0,
                temperature_k: 230.0,
            },
            dt,
            before.conductance_w_k,
            0.0,
            0.0,
            controls(),
            |h| Ok(h / 500.0),
            |h| Ok(h / 2000.0),
            || Ok(()),
        )
        .unwrap();
        let expected_mass = state.deposit.mass_kg - flux * before.front.area_m2 * dt;
        close(result.state.deposit.mass_kg, expected_mass, 3e-14);
        geometry
            .admit_density(
                Law.properties(result.phase.surface_temperature_k)
                    .unwrap()
                    .solid_density_kg_m3,
            )
            .unwrap();
        let after = geometry
            .at_mass(result.state.deposit.mass_kg, 0.1, 0.003)
            .unwrap();
        assert!(after.front.extent_m > before.front.extent_m);
        assert!(after.layer_resistance_k_w > before.layer_resistance_k_w);
        assert!(after.conductance_w_k < before.conductance_w_k);
        assert_eq!(after.contact_resistance_k_w, before.contact_resistance_k_w);
        assert!(result.mass_residual_kg.abs() < 1e-18);
        assert!(result.energy_residual_j.abs() < 1e-9);
        close(
            result.substrate_energy_change_j + result.deposit_energy_change_j,
            -result.phase.vapor_energy_j,
            1e-9,
        );
        state = result.state;
    }
}

#[test]
fn density_mismatch_bare_layers_and_unrepresentable_circuits_refuse() {
    let geometry = ConductiveCoatingGeometry::planar(900.0, 0.01).unwrap();
    assert!(geometry.admit_density(900.0).is_ok());
    assert!(geometry.admit_density(901.0).is_err());
    for (mass, k, contact) in [
        (0.0, 1.0, 0.0),
        (-1.0, 1.0, 0.0),
        (f64::NAN, 1.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, -1.0, 0.0),
        (1.0, f64::INFINITY, 0.0),
        (1.0, 1.0, -1.0),
        (1.0, 1.0, f64::INFINITY),
        (1.0, 1.0, f64::MAX),
        (1.0, f64::from_bits(1), 0.0),
    ] {
        assert!(
            geometry.at_mass(mass, k, contact).is_err(),
            "{mass}, {k}, {contact}"
        );
    }
    assert!(ConductiveCoatingGeometry::spherical(900.0, 0.0).is_err());
    assert!(ConductiveCoatingGeometry::cylindrical(900.0, 0.1, 0.0).is_err());
    assert!(ConductiveCoatingGeometry::spherical(900.0, 1e-200).is_err());
    assert!(ConductiveCoatingGeometry::cylindrical(900.0, 1e200, 1e200).is_err());
}
