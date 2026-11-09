//! G0/G1/G3 source arithmetic, geometric mass continuity and Robin transport.
use fs_convection::ThermalConductivity;
use fs_convection::tube_bank::{TubeBankError, TubeBankGeometry, TubeBankLayout, evaluate};
use fs_qty::{Length, Temperature};

fn geometry(layout: TubeBankLayout, rows: u32) -> TubeBankGeometry {
    TubeBankGeometry::try_new(
        layout,
        rows,
        Length::new(0.01),
        Length::new(0.02),
        Length::new(0.02),
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_source_means_include_gap_reynolds_wall_prandtl_and_actual_layout() {
    // Independent 60-digit decimal evaluation of the two source formulas at
    // Re_max=20000, Pr=7, Pr_surface=3.5, four rows, equal pitches.
    // Equal pitches still represent two physically distinct arrangements.
    for (layout, expected) in [
        (TubeBankLayout::Inline, 298.36506027896674),
        (TubeBankLayout::Staggered, 284.16385885169537),
    ] {
        let result = evaluate(geometry(layout, 4), 10_000.0, 7.0, 3.5).unwrap();
        near(result.evidence().value, expected);
        assert_eq!(result.groups()["Re_max"], 20_000.0);
        assert_eq!(result.groups()["Pr_surface"], 3.5);
        assert!(result.evidence().model.in_domain);
        assert_eq!(result.card().discrepancy_rel, f64::INFINITY);
        assert!(result.source().identifier.contains("28978533677980"));
        let uncorrected = evaluate(geometry(layout, 4), 10_000.0, 7.0, 7.0).unwrap();
        near(
            result.evidence().value / uncorrected.evidence().value,
            2_f64.sqrt().sqrt(),
        );
    }
}

#[test]
fn g0_g3_gap_continuity_includes_two_staggered_passages_and_geometric_scaling() {
    near(
        geometry(TubeBankLayout::Inline, 4)
            .maximum_gap_speed(3.0)
            .unwrap(),
        6.0,
    );
    for scale in [0.1, 1.0, 100.0] {
        let bank = TubeBankGeometry::try_new(
            TubeBankLayout::Staggered,
            4,
            Length::new(0.01 * scale),
            Length::new(0.02 * scale),
            Length::new(0.01 * scale),
        )
        .unwrap();
        // Two diagonal gaps of sqrt(2)D-D carry the inlet streamtube 2D.
        near(bank.gap_velocity_ratio(), 1.0 / (2_f64.sqrt() - 1.0));
    }
    assert!(
        geometry(TubeBankLayout::Inline, 4)
            .maximum_gap_speed(0.0)
            .is_err()
    );
}

#[test]
fn g1_finite_rows_use_the_published_table_and_disclosed_interpolation() {
    for (layout, first_rows, first_factor) in [
        (TubeBankLayout::Inline, 1, 0.70),
        (TubeBankLayout::Staggered, 2, 0.76),
    ] {
        let full = evaluate(geometry(layout, 16), 5000.0, 0.71, 0.71)
            .unwrap()
            .evidence()
            .value;
        for (rows, factor) in [
            (first_rows, first_factor),
            (5, 0.93),
            (6, 0.945),
            (7, 0.96),
            (10, 0.98),
            (13, 0.99),
            (16, 1.0),
            (30, 1.0),
        ] {
            let nu = evaluate(geometry(layout, rows), 5000.0, 0.71, 0.71)
                .unwrap()
                .evidence()
                .value;
            near(nu / full, factor);
        }
    }
}

#[test]
fn g0_actual_gap_regime_geometry_and_source_property_limits_refuse() {
    for (layout, transverse, longitudinal, rows) in [
        (TubeBankLayout::Inline, 0.02, 0.02, 0),
        (TubeBankLayout::Staggered, 0.02, 0.02, 1),
        (TubeBankLayout::Inline, 0.01, 0.02, 4),
        (TubeBankLayout::Inline, 0.02, 0.01, 4),
        (TubeBankLayout::Staggered, 0.025, 0.01, 4),
        (TubeBankLayout::Staggered, f64::NAN, 0.02, 4),
    ] {
        assert!(
            TubeBankGeometry::try_new(
                layout,
                rows,
                Length::new(0.01),
                Length::new(transverse),
                Length::new(longitudinal)
            )
            .is_err()
        );
    }
    let bank = geometry(TubeBankLayout::Inline, 4);
    for (re, pr, wall_pr) in [
        (499.0, 0.71, 0.71),
        (100_001.0, 0.71, 0.71),
        (1000.0, 0.6, 0.71),
        (1000.0, 7.0, 0.6),
        (1000.0, 501.0, 7.0),
    ] {
        assert!(matches!(
            evaluate(bank, re, pr, wall_pr),
            Err(TubeBankError::OutOfDomain { .. })
        ));
    }
    for re in [500.0, 100_000.0] {
        assert!(evaluate(bank, re, 0.7, 500.0).is_ok());
    }
    assert!(evaluate(bank, 0.0, 0.71, 0.71).is_err());
}

#[test]
fn g3_typed_coefficient_and_robin_row_retain_heat_direction_and_sources() {
    let result = evaluate(geometry(TubeBankLayout::Inline, 4), 10_000.0, 7.0, 3.5).unwrap();
    let k = ThermalConductivity::new(0.5);
    let coefficient = result.heat_transfer_coefficient(k).unwrap();
    near(coefficient.value.value(), 50.0 * 298.36506027896674);
    for bulk in [300.0, 400.0] {
        let row = result.robin_boundary(k, Temperature::new(bulk)).unwrap();
        assert_eq!(row.coefficient(), &coefficient);
        let fs_conduction::ThermalBc::Robin { htc, t_ref } = row.boundary_condition() else {
            panic!("Robin required")
        };
        near(
            htc.at(0) * (350.0 - t_ref.at(0)),
            coefficient.value.value() * (350.0 - bulk),
        );
    }
    assert_ne!(
        coefficient.provenance,
        result
            .heat_transfer_coefficient(ThermalConductivity::new(1.0))
            .unwrap()
            .provenance
    );
    assert!(
        result
            .heat_transfer_coefficient(ThermalConductivity::new(0.0))
            .is_err()
    );
    assert!(result.robin_boundary(k, Temperature::new(0.0)).is_err());
}
