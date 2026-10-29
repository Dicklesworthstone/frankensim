use super::*;

fn column(pressure: f64, interface: f64, lower: f64, upper: f64) -> TwoLayerColumn {
    TwoLayerColumn::new(
        Length::new(0.0),
        Pressure::new(pressure),
        Length::new(interface),
        Density::new(lower),
        Density::new(upper),
    )
    .unwrap()
}
fn coefficient(value: f64) -> DischargeCoefficient {
    DischargeCoefficient::new(
        Dimensionless::new(value),
        SourceProvenance::new("Analytic orifice fixture", "analytic-vent-v1"),
        0.0,
        ToleranceBasis::Analytic,
    )
    .unwrap()
}
fn opening(bottom: f64, top: f64, width: f64, cd: f64) -> WallOpening {
    WallOpening::new(
        Length::new(bottom),
        Length::new(top),
        Length::new(width),
        coefficient(cd),
    )
    .unwrap()
}
fn solve(
    opening: &WallOpening,
    left: TwoLayerColumn,
    right: TwoLayerColumn,
    gravity: f64,
) -> VentFlow {
    wall_vent_flow(opening, left, right, Acceleration::new(gravity), || false).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0e-300),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn mass_balance(flow: &VentFlow) -> f64 {
    let r = flow.mass_rates();
    r.left.iter().chain(&r.right).map(|q| q.value()).sum()
}
fn gross(flow: &VentFlow) -> [[f64; 4]; 2] {
    let mut rates = [[0.0; 4]; 2];
    for s in flow.slabs() {
        let l = s.left_layer.index();
        let r = s.right_layer.index();
        if s.mass_rate.value() > 0.0 {
            rates[0][2 * l + r] += s.mass_rate.value();
        } else {
            rates[1][2 * r + l] -= s.mass_rate.value();
        }
    }
    rates
}

#[test]
fn equal_columns_and_explicit_closure_are_exactly_sealed() {
    let state = column(100_000.0, 1.25, 1.2, 0.6);
    let same = solve(&opening(0.0, 2.0, 1.0, 0.7), state, state, 9.81);
    assert!(same.slabs().is_empty());
    assert_eq!(same.mass_rates().left, [MassFlowRate::new(0.0); 2]);
    assert_eq!(same.mass_rates().right, [MassFlowRate::new(0.0); 2]);
    let unequal = column(100_050.0, 0.5, 2.0, 1.0);
    for vent in [opening(0.0, 2.0, 0.0, 0.7), opening(1.0, 1.0, 2.0, 0.7)] {
        let flow = solve(&vent, state, unequal, 9.81);
        assert_eq!(flow, same);
        let heat = flow
            .advect_enthalpy(
                [SpecificEnthalpy::new(1_000.0); 2],
                [SpecificEnthalpy::new(-300.0); 2],
            )
            .unwrap();
        assert_eq!(heat.left, [Power::new(0.0); 2]);
        assert_eq!(heat.right, [Power::new(0.0); 2]);
    }
}

#[test]
fn uniform_pressure_drop_reuses_original_orifice_normalization() {
    let discharge = DischargeCoefficient::sharp_edged().unwrap();
    let original =
        sharp_edged_orifice_loss("reference", Area::new(0.8), Density::new(1.2)).unwrap();
    assert_eq!(discharge.reference_loss().source, original.source);
    assert_eq!(
        discharge.reference_loss().uncertainty_rel,
        original.uncertainty_rel
    );
    assert_eq!(
        discharge.reference_loss().tolerance_basis,
        original.tolerance_basis
    );
    assert!(discharge.reference_loss().regime_audit_card().is_some());
    let vent = WallOpening::new(
        Length::new(0.0),
        Length::new(2.0),
        Length::new(0.4),
        discharge,
    )
    .unwrap();
    let left = column(100_020.0, 3.0, 1.2, 1.2);
    let right = column(100_000.0, 3.0, 1.2, 1.2);
    let flow = solve(&vent, left, right, 9.81);
    assert_eq!(flow.slabs().len(), 1);
    let expected = 1.2 * det::sqrt(20.0 / original.resistance.value());
    close(flow.slabs()[0].mass_rate.value(), expected, 2.0e-15);
    close(flow.mass_rates().left[0].value(), -expected, 2.0e-15);
    assert_eq!(mass_balance(&flow), 0.0);
    let reverse = solve(&vent, right, left, 9.81);
    assert_eq!(
        reverse.slabs()[0].mass_rate.value(),
        -flow.slabs()[0].mass_rate.value()
    );
}

#[test]
fn reversal_uses_actual_upstream_density_with_zero_gravity() {
    let vent = opening(0.0, 2.0, 0.4, 0.6);
    let left = column(100_020.0, 3.0, 1.2, 1.2);
    let right = column(100_000.0, 3.0, 0.8, 0.8);
    let flow = solve(&vent, left, right, 0.0);
    close(
        flow.slabs()[0].mass_rate.value(),
        0.6 * 0.8 * det::sqrt(2.0 * 1.2 * 20.0),
        2.0e-15,
    );
    assert_eq!(flow.slabs()[0].upstream_density, Density::new(1.2));
    let reversed_pressure = solve(
        &vent,
        column(100_000.0, 3.0, 1.2, 1.2),
        column(100_020.0, 3.0, 0.8, 0.8),
        0.0,
    );
    close(
        reversed_pressure.slabs()[0].mass_rate.value(),
        -0.6 * 0.8 * det::sqrt(2.0 * 0.8 * 20.0),
        2.0e-15,
    );
    assert_eq!(
        reversed_pressure.slabs()[0].upstream_density,
        Density::new(0.8)
    );
    assert_eq!(mass_balance(&reversed_pressure), 0.0);
}

#[test]
fn neutral_plane_counterflow_has_independent_analytic_integrals() {
    let vent = opening(0.0, 2.0, 0.5, 0.7);
    let left = column(99_990.0, 3.0, 1.0, 1.0);
    let right = column(100_000.0, 3.0, 2.0, 2.0);
    let flow = solve(&vent, left, right, 10.0);
    let slabs = flow.slabs();
    assert_eq!(slabs.len(), 2);
    assert_eq!(slabs[0].top, Length::new(1.0));
    assert_eq!(slabs[1].bottom, Length::new(1.0));
    assert_eq!(slabs[0].pressure_difference_top, Pressure::new(0.0));
    assert_eq!(slabs[1].pressure_difference_bottom, Pressure::new(0.0));
    close(
        slabs[0].mass_rate.value(),
        -(2.0 / 3.0) * 0.7 * 0.5 * det::sqrt(40.0),
        2.0e-15,
    );
    close(
        slabs[1].mass_rate.value(),
        (2.0 / 3.0) * 0.7 * 0.5 * det::sqrt(20.0),
        2.0e-15,
    );
    assert_eq!(mass_balance(&flow), 0.0);
    let reverse = solve(&vent, right, left, 10.0);
    for (a, b) in slabs.iter().zip(reverse.slabs()) {
        assert_eq!(a.bottom, b.bottom);
        assert_eq!(a.top, b.top);
        assert_eq!(a.mass_rate.value(), -b.mass_rate.value());
        assert_eq!(a.upstream_density, b.upstream_density);
    }
}

#[test]
fn both_layer_interfaces_and_three_neutral_planes_fill_the_six_slab_bound() {
    let left = column(100_002.5, 1.0, 2.0, 1.0);
    let right = column(100_000.0, 2.0, 1.5, 0.5);
    let flow = solve(&opening(0.0, 3.0, 0.8, 0.7), left, right, 10.0);
    assert_eq!(flow.slabs().len(), MAX_VENT_SLABS);
    let donors = [2.0, 1.5, 1.5, 1.0, 1.0, 0.5];
    let signs = [1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
    for (i, slab) in flow.slabs().iter().enumerate() {
        assert_eq!(slab.bottom, Length::new(0.5 * i as f64));
        assert_eq!(slab.top, Length::new(0.5 * (i + 1) as f64));
        assert_eq!(slab.height, Length::new(0.5));
        assert_eq!(
            slab.left_layer,
            if i < 2 { Layer::Lower } else { Layer::Upper }
        );
        assert_eq!(
            slab.right_layer,
            if i < 4 { Layer::Lower } else { Layer::Upper }
        );
        assert_eq!(slab.upstream_density, Density::new(donors[i]));
        let expected = signs[i] * (2.0 / 3.0) * 0.7 * 0.8 * 0.5 * det::sqrt(2.0 * donors[i] * 2.5);
        close(slab.mass_rate.value(), expected, 2.0e-15);
    }
    assert!(mass_balance(&flow).abs() < 4.0e-16);
    assert_eq!(flow, solve(&opening(0.0, 3.0, 0.8, 0.7), left, right, 10.0));
}

#[test]
fn arbitrary_vertical_refinement_preserves_each_directional_layer_transfer() {
    let left = column(100_002.5, 1.0, 2.0, 1.0);
    let right = column(100_000.0, 2.0, 1.5, 0.5);
    let full = gross(&solve(&opening(0.0, 3.0, 0.8, 0.7), left, right, 10.0));
    let heights = [0.0, 0.23, 0.88, 1.31, 2.17, 2.64, 3.0];
    let mut refined = [[0.0; 4]; 2];
    for pair in heights.windows(2) {
        let part = gross(&solve(
            &opening(pair[0], pair[1], 0.8, 0.7),
            left,
            right,
            10.0,
        ));
        for direction in 0..2 {
            for layer_pair in 0..4 {
                refined[direction][layer_pair] += part[direction][layer_pair];
            }
        }
    }
    for direction in 0..2 {
        for layer_pair in 0..4 {
            if full[direction][layer_pair] == 0.0 {
                assert_eq!(refined[direction][layer_pair], 0.0);
            } else {
                close(
                    refined[direction][layer_pair],
                    full[direction][layer_pair],
                    1.0e-14,
                );
            }
        }
    }
}

#[test]
fn zero_net_mass_counterflow_still_advects_nonzero_donor_enthalpy() {
    // rho_R/rho_L=8 makes the lower height 1 m balance the upper height 2 m:
    // sqrt(8)*1^(3/2) = sqrt(1)*2^(3/2), with a neutral plane at z=1 m.
    let flow = solve(
        &opening(0.0, 3.0, 1.0, 0.5),
        column(99_993.0, 4.0, 1.0, 1.0),
        column(100_000.0, 4.0, 8.0, 8.0),
        1.0,
    );
    let mass = (2.0 / 3.0) * 0.5 * det::sqrt(112.0);
    assert!(flow.mass_rates().left[0].value().abs() < 2.0e-15);
    let energy = flow
        .advect_enthalpy(
            [SpecificEnthalpy::new(1_000.0); 2],
            [SpecificEnthalpy::new(200.0); 2],
        )
        .unwrap();
    close(energy.left[0].value(), -800.0 * mass, 2.0e-15);
    close(energy.right[0].value(), 800.0 * mass, 2.0e-15);
    assert_eq!(energy.left[1].value(), 0.0);
    assert_eq!(energy.right[1].value(), 0.0);
    assert_eq!(energy.left[0].value(), -energy.right[0].value());
    let shifted = flow
        .advect_enthalpy(
            [SpecificEnthalpy::new(-9_000.0); 2],
            [SpecificEnthalpy::new(-9_800.0); 2],
        )
        .unwrap();
    close(shifted.left[0].value(), energy.left[0].value(), 2.0e-14);
}

#[test]
fn hydrostatic_reference_and_world_height_changes_preserve_exchange() {
    let left = column(100_002.5, 1.0, 2.0, 1.0);
    let right = column(100_000.0, 2.0, 1.5, 0.5);
    let original = solve(&opening(0.0, 3.0, 0.8, 0.7), left, right, 10.0);
    let relocated = |c: TwoLayerColumn| {
        TwoLayerColumn::new(
            Length::new(1_001.5),
            c.pressure_at(Length::new(1.5), Acceleration::new(10.0))
                .unwrap(),
            Length::new(c.interface_height().value() + 1_000.0),
            c.density(Layer::Lower),
            c.density(Layer::Upper),
        )
        .unwrap()
    };
    let moved = solve(
        &opening(1_000.0, 1_003.0, 0.8, 0.7),
        relocated(left),
        relocated(right),
        10.0,
    );
    assert_eq!(gross(&original), gross(&moved));
    assert_eq!(original.mass_rates(), moved.mass_rates());
}

#[test]
fn tiny_counterflow_near_a_global_coordinate_boundary_is_not_snapped_away() {
    let datum = 1_000_000_000.0;
    let make = |pressure, rho| {
        TwoLayerColumn::new(
            Length::new(datum),
            Pressure::new(pressure),
            Length::new(datum + 2.0),
            Density::new(rho),
            Density::new(rho),
        )
        .unwrap()
    };
    let left = make(100_000.0 - 1.0e-8, 1.0);
    let right = make(100_000.0, 2.0);
    let flow = solve(&opening(datum, datum + 1.0, 1.0, 0.7), left, right, 1.0);
    assert_eq!(flow.slabs().len(), 2);
    let tiny = &flow.slabs()[0];
    assert!(tiny.height.value() > 0.0);
    assert_eq!(tiny.bottom, tiny.top); // World-coordinate subtraction loses it.
    assert!(tiny.mass_rate.value() < 0.0);
    let dp = right.reference_pressure().value() - left.reference_pressure().value();
    close(
        tiny.mass_rate.value(),
        -(2.0 / 3.0) * 0.7 * dp * det::sqrt(4.0 * dp),
        2.0e-15,
    );
}

#[test]
fn invalid_domains_arithmetic_and_cancellation_refuse_without_partial_results() {
    assert!(
        TwoLayerColumn::new(
            Length::new(0.0),
            Pressure::new(1.0),
            Length::new(1.0),
            Density::new(0.0),
            Density::new(1.0)
        )
        .is_err()
    );
    for cd in [0.0, -0.1, 1.01, f64::NAN] {
        assert!(
            DischargeCoefficient::new(
                Dimensionless::new(cd),
                SourceProvenance::new("source", "id"),
                0.0,
                ToleranceBasis::Analytic
            )
            .is_err()
        );
    }
    assert!(
        DischargeCoefficient::new(
            Dimensionless::new(0.7),
            SourceProvenance::new("", "id"),
            0.0,
            ToleranceBasis::Analytic
        )
        .is_err()
    );
    assert!(
        DischargeCoefficient::new(
            Dimensionless::new(0.7),
            SourceProvenance::new("source", "id"),
            1.0,
            ToleranceBasis::Analytic
        )
        .is_err()
    );
    assert!(
        WallOpening::new(
            Length::new(1.0),
            Length::new(0.0),
            Length::new(1.0),
            coefficient(0.7)
        )
        .is_err()
    );
    assert!(
        WallOpening::new(
            Length::new(0.0),
            Length::new(1.0e-200),
            Length::new(1.0e-200),
            coefficient(0.7)
        )
        .is_err()
    );
    let vent = opening(0.0, 3.0, 0.8, 0.7);
    let left = column(100_002.5, 1.0, 2.0, 1.0);
    let right = column(100_000.0, 2.0, 1.5, 0.5);
    assert!(wall_vent_flow(&vent, left, right, Acceleration::new(-1.0), || false).is_err());
    assert!(
        wall_vent_flow(
            &vent,
            column(1.0, 1.0, 2.0, 1.0),
            right,
            Acceleration::new(10.0),
            || false
        )
        .is_err()
    );
    let expected = solve(&vent, left, right, 10.0);
    let mut calls = 0;
    wall_vent_flow(&vent, left, right, Acceleration::new(10.0), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, 4, calls] {
        let mut seen = 0;
        let error = wall_vent_flow(&vent, left, right, Acceleration::new(10.0), || {
            seen += 1;
            seen == stop
        })
        .unwrap_err();
        assert_eq!(error, VentError::Cancelled);
        assert_eq!(solve(&vent, left, right, 10.0), expected);
    }
    assert!(
        expected
            .advect_enthalpy(
                [SpecificEnthalpy::new(f64::NAN); 2],
                [SpecificEnthalpy::new(0.0); 2]
            )
            .is_err()
    );
    let large_flow = solve(&opening(0.0, 3.0, 2.4, 0.7), left, right, 10.0);
    assert!(
        large_flow
            .advect_enthalpy(
                [SpecificEnthalpy::new(f64::MAX); 2],
                [SpecificEnthalpy::new(0.0); 2]
            )
            .is_err()
    );
}
