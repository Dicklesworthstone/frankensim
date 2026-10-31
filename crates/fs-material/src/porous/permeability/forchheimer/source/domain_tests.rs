use super::*;

fn support() -> PorousGasSupport {
    PorousGasSupport {
        temperature_k: (250.0, 600.0),
        pressure_pa: (10000.0, 1e6),
        maximum_superficial_speed_m_s: 20.0,
    }
}

#[test]
fn original_source_coordinates_intersect_actual_temperature_pressure_and_speed() {
    let mut domain = support();
    for (axis, quantity, value, bounds) in [
        (
            "temperature",
            QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
            300.0,
            (290.0, 320.0),
        ),
        (
            "pressure",
            QuantitySpec::dimensional(fs_qty::Pressure::DIMS),
            150000.0,
            (120000.0, 200000.0),
        ),
        (
            SOURCE_SPEED_AXIS,
            QuantitySpec::dimensional(fs_qty::Velocity::DIMS),
            0.1,
            (0.0, 0.25),
        ),
    ] {
        admit_source_coordinate(&mut domain, axis, Some(quantity), Some(value), bounds).unwrap();
    }
    assert_eq!(domain.temperature_k, (290.0, 320.0));
    assert_eq!(domain.pressure_pa, (120000.0, 200000.0));
    assert_eq!(domain.maximum_superficial_speed_m_s, 0.25);
    admit_source_coordinate(
        &mut domain,
        SOURCE_SPEED_AXIS,
        Some(QuantitySpec::dimensional(fs_qty::Velocity::DIMS)),
        Some(0.1),
        (0.0, 2.0),
    )
    .unwrap();
    assert_eq!(domain.maximum_superficial_speed_m_s, 0.25);
}

#[test]
fn source_speed_requires_its_actual_basis_and_supported_zero_limit() {
    let velocity = Some(QuantitySpec::dimensional(fs_qty::Velocity::DIMS));
    for (quantity, bounds, value) in [
        (None, (0.0, 1.0), 0.1),
        (Some(QuantitySpec::dimensional(Dims::NONE)), (0.0, 1.0), 0.1),
        (velocity, (0.01, 1.0), 0.1),
        (velocity, (-1.0, 1.0), 0.1),
        (velocity, (0.0, 1.0), -0.1),
        (velocity, (0.0, f64::INFINITY), 0.1),
    ] {
        assert!(
            admit_source_coordinate(
                &mut support(),
                SOURCE_SPEED_AXIS,
                quantity,
                Some(value),
                bounds
            )
            .is_err()
        );
    }
    assert!(
        admit_source_coordinate(
            &mut support(),
            "pore_velocity",
            velocity,
            Some(0.1),
            (0.0, 1.0)
        )
        .is_err()
    );
}

#[test]
fn unowned_axes_and_disjoint_source_support_refuse() {
    assert!(
        admit_source_coordinate(
            &mut support(),
            "specimen_stress",
            Some(QuantitySpec::dimensional(fs_qty::Pressure::DIMS)),
            Some(150000.0),
            (120000.0, 200000.0)
        )
        .is_err()
    );
    assert!(
        admit_source_coordinate(
            &mut support(),
            "temperature",
            Some(QuantitySpec::dimensional(fs_qty::Temperature::DIMS)),
            Some(700.0),
            (650.0, 800.0)
        )
        .is_err()
    );
    assert!(
        admit_source_coordinate(
            &mut support(),
            "pressure",
            Some(QuantitySpec::dimensional(fs_qty::Pressure::DIMS)),
            Some(0.0),
            (0.0, 200000.0)
        )
        .is_err()
    );
}
