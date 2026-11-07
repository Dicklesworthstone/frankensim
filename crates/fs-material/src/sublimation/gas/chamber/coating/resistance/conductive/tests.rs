use super::*;
use core::cell::Cell;
use core::f64::consts::PI;
use fs_thermal::sublimation::coating::resistance::{
    ResistiveSurfaceState, heat_resistive_surface_with_interface,
};
use fs_thermal::sublimation::{
    SolidSurfaceState, SolverConfig, SublimationLaw, SublimationProperties,
};

fn controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}

fn affine_curve() -> ConductivityCurve {
    // Manufactured original chart k(T) = 0.01*T with a crossed interior knot.
    ConductivityCurve::new(vec![[200.0, 2.0], [350.0, 3.5], [1000.0, 10.0]]).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual}, expected {expected}, tolerance {tolerance}"
    );
}

#[test]
fn nonlinear_coating_and_contact_recover_both_heat_directions_and_tangents() {
    let curve = affine_curve();
    let geometry = ConductiveCoatingGeometry::planar(1000.0, 0.02).unwrap();
    // m/(rho*A) = 0.04 m, C = A/d = 0.5 m, R_contact = 0.2 K/W.
    let interface =
        ConductiveCoatingInterface::new(geometry, 0.8, &curve, 0.004, controls()).unwrap();
    assert!(core::ptr::eq(interface.conductivity_curve(), &curve));
    close(interface.front().extent_m, 0.04, 1e-15);
    close(interface.contact_resistance_k_w(), 0.2, 1e-15);
    for (support, deposit, face, power) in
        [(580.0, 300.0, 500.0, 400.0), (220.0, 500.0, 300.0, -400.0)]
    {
        // Q = C*0.005*(T_face^2-T_deposit^2), T_support-T_face = R*Q.
        let response = interface.evaluate([support, deposit], || Ok(())).unwrap();
        close(response.power_to_reservoir_w, power, 1e-8);
        close(response.material_face_temperature_k, face, 1e-8);
        close(response.conductance_w_k, power / (support - deposit), 1e-10);
        close(response.heat_residual_w, 0.0, 1e-8);
        let denominator = 1.0 + 0.2 * 0.5 * 0.01 * face;
        close(
            response.endpoint_derivative_w_k[0],
            0.5 * 0.01 * face / denominator,
            1e-10,
        );
        close(
            response.endpoint_derivative_w_k[1],
            -0.5 * 0.01 * deposit / denominator,
            1e-10,
        );
        close(
            interface.heat_rate_w(support, deposit).unwrap(),
            power,
            1e-8,
        );
        // Circuit metadata can use the actual profile average without subtracting
        // a large contact resistance from a nearly equal total resistance.
        let average_k = curve
            .integral_w_m(deposit, response.material_face_temperature_k)
            .unwrap()
            / (response.material_face_temperature_k - deposit);
        let circuit = geometry.at_mass(0.8, average_k, 0.004).unwrap();
        close(circuit.conductance_w_k, response.conductance_w_k, 1e-10);
    }
}

#[test]
fn crossed_original_knots_change_flux_and_preserve_the_equilibrium_tangent() {
    let curve = ConductivityCurve::new(vec![
        [200.0, 1.0],
        [400.0, 4.0],
        [700.0, 2.0],
        [1000.0, 8.0],
    ])
    .unwrap();
    let geometry = ConductiveCoatingGeometry::planar(1000.0, 0.02).unwrap();
    let interface =
        ConductiveCoatingInterface::new(geometry, 0.8, &curve, 0.0, controls()).unwrap();
    // Exact trapezoids on [300,400], [400,700], [700,800] total 1525 W/m.
    // Freezing the deposit's k(300)=2.5 would instead give only 625 W.
    let response = interface.evaluate([800.0, 300.0], || Ok(())).unwrap();
    close(response.power_to_reservoir_w, 762.5, 1e-10);
    assert!((response.power_to_reservoir_w - 625.0).abs() > 100.0);
    assert_eq!(response.material_face_temperature_k, 800.0);
    assert_eq!(response.iterations, 0);
    close(interface.heat_rate_w(300.0, 800.0).unwrap(), -762.5, 1e-10);
    let equilibrium = interface.evaluate([400.0, 400.0], || Ok(())).unwrap();
    assert_eq!(equilibrium.power_to_reservoir_w, 0.0);
    close(equilibrium.conductance_w_k, 2.0, 1e-14);
    assert_eq!(equilibrium.endpoint_derivative_w_k, [2.0, -2.0]);
}

#[test]
fn curved_coatings_use_the_radial_fourier_metric_and_inner_contact_area() {
    let curve = affine_curve();
    let density = 1000.0;
    let inner = 0.02_f64;
    let outer = 0.03_f64;
    let length = 0.1;
    let cylinder = ConductiveCoatingGeometry::cylindrical(density, inner, length).unwrap();
    let sphere = ConductiveCoatingGeometry::spherical(density, inner).unwrap();
    for (geometry, mass, multiplier, contact_area, outer_area) in [
        (
            cylinder,
            PI * (outer * outer - inner * inner) * length * density,
            2.0 * PI * length / (outer / inner).ln(),
            2.0 * PI * inner * length,
            2.0 * PI * outer * length,
        ),
        (
            sphere,
            (4.0 * PI / 3.0) * (outer.powi(3) - inner.powi(3)) * density,
            4.0 * PI * inner * outer / (outer - inner),
            4.0 * PI * inner * inner,
            4.0 * PI * outer * outer,
        ),
    ] {
        let contact_areal = 1e-3;
        let contact_total = contact_areal / contact_area;
        let power = multiplier * 800.0;
        let support = 500.0 + contact_total * power;
        let interface =
            ConductiveCoatingInterface::new(geometry, mass, &curve, contact_areal, controls())
                .unwrap();
        let response = interface.evaluate([support, 300.0], || Ok(())).unwrap();
        close(response.power_to_reservoir_w, power, 2e-8);
        close(response.material_face_temperature_k, 500.0, 1e-8);
        close(interface.contact_resistance_k_w(), contact_total, 1e-14);
        close(interface.front().extent_m, outer - inner, 1e-15);
        close(interface.front().area_m2, outer_area, 1e-14);
        assert!(interface.front().area_m2 > interface.geometry().contact_area_m2());
    }
}

#[test]
fn sub_radius_ulp_shells_retain_positive_mass_and_finite_original_integral_flux() {
    let curve = affine_curve();
    let thickness = 1e-20;
    let radius = 1.0;
    let length = 0.5;
    assert_eq!(radius + thickness, radius);
    for (geometry, multiplier) in [
        (
            ConductiveCoatingGeometry::spherical(1000.0, radius).unwrap(),
            4.0 * PI * radius * radius / thickness,
        ),
        (
            ConductiveCoatingGeometry::cylindrical(1000.0, radius, length).unwrap(),
            2.0 * PI * radius * length / thickness,
        ),
    ] {
        let mass = geometry.geometry().at_extent(thickness).unwrap().mass_kg;
        let interface =
            ConductiveCoatingInterface::new(geometry, mass, &curve, 0.0, controls()).unwrap();
        assert!(interface.front().extent_m > 0.0);
        close(interface.front().extent_m / thickness, 1.0, 1e-12);
        let power = interface.heat_rate_w(400.0, 300.0).unwrap();
        // Exact material integral is 0.005*(400^2-300^2) = 350 W/m.
        close(power / (multiplier * 350.0), 1.0, 1e-12);
    }
}

struct HeatLaw;

impl SublimationLaw for HeatLaw {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 1000.0]
    }

    fn properties(&self, _: f64) -> Result<SublimationProperties, SublimationError> {
        Ok(SublimationProperties {
            saturation_pressure_pa: 1.0,
            saturation_pressure_derivative_pa_k: 0.0,
            latent_heat_j_kg: 1e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}

#[test]
fn original_two_temperature_owner_accepts_the_nonlinear_material_circuit_once() {
    let curve = affine_curve();
    let geometry = ConductiveCoatingGeometry::planar(1000.0, 1e-4).unwrap();
    let interface =
        ConductiveCoatingInterface::new(geometry, 1e-4, &curve, 1e-4, controls()).unwrap();
    // C=0.1 m and contact R=1 K/W. The manufactured final temperatures
    // (support=580, face=500, deposit=300 K) carry 80 W, hence 8 J in 0.1 s.
    // Independent constant cp charts store -8 J in the support and +8 J in
    // the deposit. Their positive masses and distinct caloric datums survive.
    let initial = ResistiveSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: 500.0 * 580.8 - 1e5,
        deposit: SolidSurfaceState {
            mass_kg: 1e-4,
            specific_enthalpy_j_kg: 1000.0 * 220.0 + 2e5,
            area_m2: 1e-4,
        },
    };
    let step = heat_resistive_surface_with_interface(
        &HeatLaw,
        initial,
        0.1,
        &interface,
        0.0,
        0.0,
        SolverConfig {
            max_iterations: 256,
            relative_tolerance: 1e-12,
            absolute_energy_tolerance_j: 1e-10,
        },
        |h| Ok((h + 1e5) / 500.0),
        |h| Ok((h - 2e5) / 1000.0),
        || Ok(()),
    )
    .unwrap();
    close(step.substrate_temperature_k, 580.0, 2e-8);
    close(step.deposit_temperature_k, 300.0, 2e-8);
    close(step.interface_heat_j, 8.0, 2e-9);
    close(step.substrate_energy_change_j, -8.0, 2e-9);
    close(step.deposit_energy_change_j, 8.0, 2e-9);
    close(step.substrate_energy_residual_j, 0.0, 2e-9);
    close(step.deposit_energy_residual_j, 0.0, 2e-9);
    close(step.energy_residual_j, 0.0, 2e-9);
    assert_eq!(step.state.substrate_mass_kg, initial.substrate_mass_kg);
    assert_eq!(step.state.deposit.mass_kg, initial.deposit.mass_kg);
}

#[test]
fn source_support_and_inner_cancellation_survive_the_interface_adapter() {
    let curve = ConductivityCurve::new(vec![[300.0, 3.0], [500.0, 5.0]]).unwrap();
    let geometry = ConductiveCoatingGeometry::planar(1000.0, 1e-4).unwrap();
    // The support may lie outside the deposit's source chart. Its physical
    // contact jump leaves the ACTUAL material profile between 300 and 400 K.
    let interface =
        ConductiveCoatingInterface::new(geometry, 1e-4, &curve, (300.0 / 35.0) * 1e-4, controls())
            .unwrap();
    let response = interface.evaluate([700.0, 300.0], || Ok(())).unwrap();
    close(response.power_to_reservoir_w, 35.0, 1e-8);
    close(response.material_face_temperature_k, 400.0, 1e-8);
    let perfect = ConductiveCoatingInterface::new(geometry, 1e-4, &curve, 0.0, controls()).unwrap();
    assert!(matches!(
        perfect.evaluate([700.0, 300.0], || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert!(ConductiveCoatingInterface::new(geometry, 0.0, &curve, 0.0, controls()).is_err());
    assert!(
        ConductiveCoatingInterface::new(
            geometry,
            1e-4,
            &curve,
            0.0,
            LayerConductionConfig {
                max_iterations: 0,
                ..controls()
            },
        )
        .is_err()
    );
    assert!(ReservoirLayerGeometry::from_geometric_multiplier(0.0, 1.0).is_err());
    assert!(ReservoirLayerGeometry::from_geometric_multiplier(1.0, f64::INFINITY).is_err());

    let mut direct_polls = 0;
    assert!(matches!(
        interface.evaluate([700.0, 300.0], || {
            direct_polls += 1;
            if direct_polls == 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(direct_polls, 2);
    // An immutable closure may borrow a shorter-lived original scope gate.
    // The trait still polls it inside the original layer solve.
    let polls = Cell::new(0);
    let checkpoint = || {
        polls.set(polls.get() + 1);
        if polls.get() == 2 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    let cancellable = interface.with_checkpoint(&checkpoint);
    assert_eq!(
        cancellable.heat_rate_w(700.0, 300.0),
        Err(SublimationError::Cancelled)
    );
    assert_eq!(polls.get(), 2);
}
