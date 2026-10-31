//! G1/G3: equilibrium ionization through the original FEEC current/heat owner.
mod common;
#[path = "../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;
use common::{config, ends, layered, limits, near};
use fs_em::conduction::{ConductionControls, ConductionError, PlanarConduction};
use fs_material::ionized_gas::IonizedGasConductivity;

fn controls() -> ConductionControls {
    ConductionControls {
        linear: config(),
        current_tolerance_a: 1e-6,
        power_tolerance_w: 1e-6,
    }
}

#[test]
fn g1_ionized_material_uses_original_field_current_joule_and_circuit_laws() {
    for t in [6000.0, 8000.0] {
        let gas = fixture::material(t, 100_000.0);
        let sigma = gas.conductivity_s_m();
        let problem =
            PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
                false
            })
            .unwrap();
        let field = problem
            .solve(&ends(2.0), controls(), &mut || false)
            .unwrap();
        // Unit length and width, depth=0.2: E=(-2,0), P=sigma E² volume.
        for (e, j) in field
            .electric_field_v_m()
            .iter()
            .zip(field.current_density_a_m2())
        {
            near(e[0], -2.0);
            near(e[1], 0.0);
            near(j[0], -2.0 * sigma);
            near(j[1], 0.0);
        }
        near(field.power_w(), 0.8 * sigma);
        assert!(field.power_balance_residual_w().abs() < 1e-6);
        let resistor = problem
            .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
            .unwrap();
        near(resistor.resistance_ohm(), 1.0 / (0.2 * sigma));
        let reversed = resistor
            .field_at_voltage(-2.0, 1e-6, 1e-6, &mut || false)
            .unwrap();
        near(reversed.power_w(), field.power_w());
    }
}

#[test]
fn g0_actual_pde_and_circuit_fields_cannot_cross_the_weak_drift_domain() {
    let gas = fixture::material(8000.0, 100_000.0);
    let excessive = 2.0 * gas.maximum_electric_field_v_m();
    let problem =
        PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
            false
        })
        .unwrap();
    assert!(matches!(
        problem.solve(&ends(excessive), controls(), &mut || false),
        Err(ConductionError::Invalid(_))
    ));
    let resistor = problem
        .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    assert!(matches!(
        resistor.field_at_voltage(excessive, 1e-6, 1e-6, &mut || false),
        Err(ConductionError::Invalid(_))
    ));
}

#[test]
fn g1_small_physical_drive_does_not_require_an_out_of_domain_unit_voltage_basis() {
    let chemistry = fixture::chemistry();
    let equilibrium = fixture::equilibrium(&chemistry, 8000.0, 100_000.0);
    let mut domain = fixture::limits();
    domain.maximum_drift_thermal_ratio = 1e-8;
    let gas = IonizedGasConductivity::new(
        &chemistry,
        equilibrium,
        fixture::rates(8000.0, 100_000.0, 1.0),
        domain,
        || false,
    )
    .unwrap();
    let sigma = gas.conductivity_s_m();
    let limit = gas.maximum_electric_field_v_m();
    assert!(limit < 0.001);
    let problem =
        PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![gas; 4], limits(), &mut || {
            false
        })
        .unwrap();
    assert!(
        problem
            .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
            .is_err()
    );
    let reference_v = 0.1 * limit;
    let resistor = problem
        .resistance_at_reference_voltage(&[2, 5], &[0, 3], reference_v, controls(), &mut || false)
        .unwrap();
    fixture::relative(resistor.reference_voltage_v(), reference_v, 1e-14);
    fixture::relative(
        resistor.reference_field().electric_field_v_m()[0][0],
        -reference_v,
        1e-10,
    );
    fixture::relative(resistor.resistance_ohm(), 1.0 / (0.2 * sigma), 1e-10);
    let actual_v = 0.5 * limit;
    let field = resistor
        .field_at_voltage(actual_v, 1e-8, 1e-10, &mut || false)
        .unwrap();
    fixture::relative(field.power_w(), 0.2 * sigma * actual_v * actual_v, 1e-10);
    fixture::relative(field.current_density_a_m2()[0][0], -sigma * actual_v, 1e-10);
    let current = 0.2 * sigma * actual_v;
    let driven = problem
        .solve_current_with_reference_voltage(
            &[2, 5],
            &[0, 3],
            current,
            reference_v,
            controls(),
            &mut || false,
        )
        .unwrap();
    fixture::relative(driven.power_w(), field.power_w(), 1e-10);
    assert!(
        resistor
            .field_at_voltage(2.0 * limit, 1e-8, 1e-10, &mut || false)
            .is_err()
    );
    assert!(
        problem
            .solve_current_with_reference_voltage(
                &[2, 5],
                &[0, 3],
                4.0 * current,
                reference_v,
                controls(),
                &mut || false,
            )
            .is_err()
    );
}
