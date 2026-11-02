//! Absolute mechanical pressure and shear share one source constitutive root.

use super::*;
use fs_qty::Pressure;

fn pressure_card(pressure_slope: f64, rate_slope: f64, half_width: f64) -> MaterialCard {
    let original = card(0.0, 1000.0, 0.001);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.validity =
            claim
                .validity
                .with("pressure", 100_000.0 - half_width, 100_000.0 + half_width);
        if claim.key.name() == "dynamic_viscosity" && (pressure_slope != 0.0 || rate_slope != 0.0) {
            let dims = claim.key.dims();
            claim.value = PropertyValue::Surface2D {
                axes: ["pressure".into(), "rate".into()],
                axis_dims: [Pressure::DIMS, Frequency::DIMS],
                samples: [
                    [100_000.0 - half_width, 0.0],
                    [100_000.0 + half_width, 0.0],
                    [100_000.0 - half_width, 1000.0],
                    [100_000.0 + half_width, 1000.0],
                ]
                .into_iter()
                .map(|p| {
                    (
                        p,
                        0.1 + pressure_slope * (p[0] - 100_000.0) + rate_slope * p[1],
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            claim.validity = claim.validity.with_quantity(
                "rate",
                QuantitySpec::dimensional(Frequency::DIMS),
                0.0,
                1000.0,
            );
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn binding(source: &MaterialCard, pressure: f64) -> IsothermalRheologyLaw<'_> {
    let original = bundle(source);
    let query = original.query().with("pressure", pressure).unwrap();
    IsothermalRheologyLaw::try_new_with_pressure(
        source,
        original.at_query_point(source, &query).unwrap(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        "T",
        "rate",
        "pressure",
    )
    .unwrap()
}

#[test]
fn link_loading_shares_source_pressure_and_shear_updates_and_refuses_atomically() {
    let gate = CancelGate::new_clock_free();
    let cells = pressure_slab();
    let source = pressure_card(0.01, 0.001, 5.0);
    let binding = binding(&source, 100_000.0);
    with_cx(&gate, |cx| {
        let mut cell = pressure_flow(&binding, &cells, 100_000.0, cx, &gate);
        let mut link = cell.clone();
        for _ in 0..3 {
            cell.step_with_gas_pressure(&gate, 0.05, |_| 100_000.5)
                .unwrap();
            link.step_with_gas_link_pressure(&gate, 0.05, |_| 100_000.5)
                .unwrap();
            assert_eq!(
                cell.lattice().grid().state_bits(),
                link.lattice().grid().state_bits()
            );
            for x in 2..6 {
                assert_eq!(cell.rheology_at([x, 0, 3]), link.rheology_at([x, 0, 3]));
                assert_eq!(cell.material_at([x, 0, 3]), link.material_at([x, 0, 3]));
            }
        }
        let bits = link.lattice().grid().state_bits();
        let mass = link.ledger_mass_kg().unwrap();
        let material = link.material_at([3, 0, 3]).unwrap().clone();
        let sample = link.rheology_at([3, 0, 3]).cloned().unwrap();
        for cancel in [false, true] {
            let gate = CancelGate::new_clock_free();
            let calls = std::cell::Cell::new(0);
            assert!(
                link.step_with_gas_link_pressure(&gate, 0.05, |_| {
                    calls.set(calls.get() + 1);
                    if calls.get() == 5 {
                        if cancel {
                            gate.request();
                        } else {
                            return 0.0;
                        }
                    }
                    100_000.5
                })
                .is_err()
            );
            assert!(calls.get() >= 5);
            assert_eq!(link.lattice().grid().state_bits(), bits);
            assert_eq!(link.ledger_mass_kg().unwrap(), mass);
            assert_eq!(link.material_at([3, 0, 3]), Some(&material));
            assert_eq!(link.rheology_at([3, 0, 3]), Some(&sample));
        }
    });
}

#[test]
fn g1_pressure_and_shear_use_actual_interior_state_and_one_constitutive_root() {
    let gate = CancelGate::new_clock_free();
    let cells = pressure_slab();
    let pressure_slope = 0.01;
    let rate_slope = 0.001;
    let source = pressure_card(pressure_slope, rate_slope, 5.0);
    let law = binding(&source, 100_000.0);
    with_cx(&gate, |cx| {
        let mut flow = pressure_flow(&law, &cells, 100_000.0, cx, &gate);
        for _ in 0..5 {
            flow.step_with_gas_pressure(&gate, 0.05, |p| {
                100_000.0 + if p[0] < 4.0 * DX { 0.75 } else { -0.75 }
            })
            .unwrap();
        }
        let coordinates = [3, 0, 3];
        let sample = flow.rheology_at(coordinates).unwrap();
        let native = flow.lattice().sample_at(3, 0, 3).unwrap();
        // Independent physical conversion: rho0=1000, dx/dt=1, cs²=1/3.
        let pressure = 100_000.0 + (native.density - 1.0) * 1000.0 / 3.0;
        let material = flow.material_at(coordinates).unwrap();
        let source_pressure = material
            .query_point()
            .iter()
            .find(|(axis, _)| axis == "pressure")
            .unwrap()
            .1;
        near(source_pressure, pressure, 2e-11);
        assert!((source_pressure - 100_000.0).abs() > 1e-4);
        assert!(sample.shear_rate_per_s > 1e-3);
        let viscosity =
            0.1 + pressure_slope * (pressure - 100_000.0) + rate_slope * sample.shear_rate_per_s;
        near(
            material.property("dynamic_viscosity").unwrap().value_si(),
            viscosity,
            2e-12,
        );
        near(sample.collision_tau, 0.5 + 3.0 * viscosity, 1e-11);
        // From the actual populations, r*(a + b*r) = stress_numerator/dt.
        // The stable quadratic expression is independent of the production root.
        let a = 0.8 + 3.0 * pressure_slope * (pressure - 100_000.0);
        let b = 3.0 * rate_slope;
        let rhs = numerator(&flow, [0.0; 3], coordinates) / DT;
        let expected_rate = 2.0 * rhs / (a + (a * a + 4.0 * b * rhs).sqrt());
        near(sample.shear_rate_per_s, expected_rate, 1e-8);
        assert_eq!(material.card_identity(), source.content_hash());
        assert_eq!(material.property("density").unwrap().value_si(), 1000.0);
        for property in material.properties() {
            assert_eq!(
                property.answer().receipt.selected,
                law.material()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    });
}

#[test]
fn g3_pressure_sensitive_material_changes_motion_while_constant_limit_is_exact() {
    let gate = CancelGate::new_clock_free();
    let cells = pressure_slab();
    let constant_source = pressure_card(0.0, 0.0, 5.0);
    let constant_pressure_law = binding(&constant_source, 100_000.0);
    let fixed_law = law(&constant_source);
    let source = pressure_card(0.01, 0.001, 5.0);
    let first_law = binding(&source, 100_000.0);
    let shifted_law = binding(&source, 100_002.0);
    with_cx(&gate, |cx| {
        let mut fixed = pressure_flow(&fixed_law, &cells, 100_000.0, cx, &gate);
        let mut constant = pressure_flow(&constant_pressure_law, &cells, 100_000.0, cx, &gate);
        let mut first = pressure_flow(&first_law, &cells, 100_000.0, cx, &gate);
        let mut shifted = pressure_flow(&shifted_law, &cells, 100_002.0, cx, &gate);
        near(first.rheology_at(CELL).unwrap().collision_tau, 0.8, 1e-12);
        near(
            shifted.rheology_at(CELL).unwrap().collision_tau,
            0.86,
            1e-12,
        );
        let load = |p: [f64; 3]| if p[0] < 4.0 * DX { 0.75 } else { -0.75 };
        for _ in 0..5 {
            fixed
                .step_with_gas_pressure(&gate, 0.05, |p| 100_000.0 + load(p))
                .unwrap();
            constant
                .step_with_gas_pressure(&gate, 0.05, |p| 100_000.0 + load(p))
                .unwrap();
            first
                .step_with_gas_pressure(&gate, 0.05, |p| 100_000.0 + load(p))
                .unwrap();
            shifted
                .step_with_gas_pressure(&gate, 0.05, |p| 100_002.0 + load(p))
                .unwrap();
            assert_eq!(
                fixed.lattice().grid().state_bits(),
                constant.lattice().grid().state_bits()
            );
            assert_eq!(
                fixed.ledger_mass_kg().unwrap(),
                constant.ledger_mass_kg().unwrap()
            );
            assert_eq!(fixed.rheology_at(CELL), constant.rheology_at(CELL));
        }
        assert_ne!(
            first.lattice().grid().state_bits(),
            shifted.lattice().grid().state_bits()
        );
        assert!(
            (mean_pressure_driven_velocity(&first) - mean_pressure_driven_velocity(&shifted)).abs()
                > 1e-12
        );
    });
}

#[test]
fn g4_pressure_source_refusal_and_cancellation_preserve_the_accepted_flow() {
    let gate = CancelGate::new_clock_free();
    let source = pressure_card(0.01, 0.001, 1e-5);
    let law = IsothermalRheologyLaw::try_new_with_pressure(
        &source,
        bundle(&source),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        "T",
        "rate",
        "pressure",
    )
    .unwrap();
    let cells = pressure_slab();
    with_cx(&gate, |cx| {
        assert!(
            MaterialFreeSurfaceRheology3::from_cells(
                &law,
                DX,
                DT,
                0.1,
                domain([0.0; 3]),
                &cells,
                solve(),
                cx,
                &gate,
            )
            .is_err()
        );
        assert!(
            MaterialFreeSurfaceRheology3::from_cells_with_pressure_reference(
                &law,
                DX,
                DT,
                0.1,
                domain([0.0; 3]),
                &cells,
                solve(),
                100_001.0,
                cx,
                &gate,
            )
            .is_err()
        );
        let mut flow = pressure_flow(&law, &cells, 100_000.0, cx, &gate);
        let bits = flow.lattice().grid().state_bits();
        let mass = flow.ledger_mass_kg().unwrap();
        let carry = flow.lattice().carry();
        let conversions = flow.lattice().conversions();
        let material = flow.material_at(CELL).unwrap().clone();
        let sample = *flow.rheology_at(CELL).unwrap();
        let unchanged = |flow: &MaterialFreeSurfaceRheology3<'_, '_>| {
            assert_eq!(flow.lattice().grid().state_bits(), bits);
            assert_eq!(flow.ledger_mass_kg().unwrap(), mass);
            assert_eq!(flow.lattice().carry(), carry);
            assert_eq!(flow.lattice().conversions(), conversions);
            assert_eq!(flow.material_at(CELL), Some(&material));
            assert_eq!(flow.rheology_at(CELL), Some(&sample));
        };
        let calls = std::cell::Cell::new(0);
        assert!(matches!(
            flow.step_with_gas_pressure(&gate, 0.05, |p| {
                calls.set(calls.get() + 1);
                100_000.0 + if p[0] < 4.0 * DX { 0.75 } else { -0.75 }
            }),
            Err(FreeSurfaceRheologyError::Material(_))
        ));
        assert!(
            calls.get() > 0,
            "source refusal follows actual pressure-loaded flow"
        );
        unchanged(&flow);
        let cancel = CancelGate::new_clock_free();
        assert!(
            flow.step_with_gas_pressure(&cancel, 0.05, |_| {
                cancel.request();
                100_000.0
            })
            .is_err()
        );
        unchanged(&flow);
        let mut reference = pressure_flow(&law, &cells, 100_000.0, cx, &gate);
        flow.step_with_gate(&gate).unwrap();
        reference.step_with_gate(&gate).unwrap();
        assert_eq!(
            flow.lattice().grid().state_bits(),
            reference.lattice().grid().state_bits()
        );
        assert_eq!(flow.material_at(CELL), reference.material_at(CELL));
    });
}

#[test]
fn g0_pressure_dependent_density_or_capillarity_refuses_the_incompressible_owner() {
    for property in ["density", "surface_tension"] {
        let original = pressure_card(0.01, 0.001, 5.0);
        let mut claims = ClaimSet::new();
        for (_, claim) in original.claims().claims_ordered() {
            let mut claim = claim.clone();
            if claim.key.name() == property {
                let value = if property == "density" { 1000.0 } else { 0.001 };
                claim.value = PropertyValue::Curve {
                    abscissa: "pressure".into(),
                    abscissa_dims: Pressure::DIMS,
                    knots: vec![(99_995.0, value), (100_005.0, value * 1.01)],
                    dims: claim.key.dims(),
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
            claims.insert_claim(claim).unwrap();
        }
        let source = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
        assert!(
            IsothermalRheologyLaw::try_new_with_pressure(
                &source,
                bundle(&source),
                LiquidViscosityBasis::Dynamic,
                LiquidCapillarity::Sourced,
                "T",
                "rate",
                "pressure",
            )
            .is_err()
        );
    }
}
