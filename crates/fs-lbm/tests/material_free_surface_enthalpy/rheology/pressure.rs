//! Actual liquid pressure joins the original thermal and constitutive timeline.

use super::*;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::Pressure;

const REFERENCE: f64 = 100_000.0;

// Existing source surfaces have two response axes. Exercise T/p and p/shear
// independently, while every runtime query retains all three coordinates.
fn pressure_card(second: &str, pressure_slope: f64, slope: f64, half_width: f64) -> MaterialCard {
    let original = card(0.1, 1000.0);
    let mut claims = ClaimSet::new();
    let (range, center, dims) = match second {
        "T" => ([300.0, 500.0], 400.0, Temperature::DIMS),
        "rate" => ([0.0, 1000.0], 0.0, Frequency::DIMS),
        _ => unreachable!(),
    };
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.validity =
            claim
                .validity
                .with("pressure", REFERENCE - half_width, REFERENCE + half_width);
        if claim.key.name() == "dynamic_viscosity" && (pressure_slope != 0.0 || slope != 0.0) {
            claim.value = PropertyValue::Surface2D {
                axes: ["pressure".into(), second.into()],
                axis_dims: [Pressure::DIMS, dims],
                samples: [
                    [REFERENCE - half_width, range[0]],
                    [REFERENCE + half_width, range[0]],
                    [REFERENCE - half_width, range[1]],
                    [REFERENCE + half_width, range[1]],
                ]
                .into_iter()
                .map(|p| {
                    (
                        p,
                        0.1 + pressure_slope * (p[0] - REFERENCE) + slope * (p[1] - center),
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: claim.key.dims(),
            };
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

fn with_pressure_liquid<R>(
    card: &MaterialCard,
    evolving_pressure: bool,
    exercise: impl FnOnce(
        Result<MaterialFreeSurfaceEnthalpy3<'_, '_>, FreeSurfaceEnthalpyError>,
        &CancelGate,
    ) -> R,
) -> R {
    let requirements = card
        .claims()
        .claims_ordered()
        .map(|(_, claim)| {
            ScalarPropertyRequirement::try_with_quantity(
                claim.key.name(),
                claim.key.quantity(),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let material = resolve_material_state_point(
        card,
        &query_point(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law =
        LiquidTemperatureLaw::try_from_material(card, material, "T", LiquidViscosityBasis::Dynamic)
            .unwrap();
    let source = if evolving_pressure {
        ResolvedLiquidSurfaceState::try_for_pressure_shear_rate_evolution(
            law.reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected,
            "T",
            "rate",
            "pressure",
        )
    } else {
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            law.reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            "T",
            "rate",
        )
    }
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut cells = initial(false);
    for cell in &mut cells {
        cell.velocity_m_s = [0.003 * (cell.coordinates[2] as f64 - 3.5), 0.0, 0.0];
    }
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .with_pressure_reference(REFERENCE)
            .unwrap();
        exercise(
            scaling.free_surface3_from_cells_with_rheology(
                domain(),
                &cells,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                solve(),
                |_| 200_000.0,
                cx,
                &gate,
            ),
            &gate,
        )
    })
}

fn pressure_load(position: [f64; 3]) -> f64 {
    REFERENCE + if position[0] < 0.0025 { 0.75 } else { -0.75 }
}

fn load(liquid: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>, gate: &CancelGate) {
    for _ in 0..3 {
        liquid
            .step_with_heat_and_gas_pressure(&[], gate, 0.05, pressure_load)
            .unwrap();
    }
}

#[test]
fn actual_pressure_and_shear_close_the_same_root_and_molecular_heat() {
    let source = pressure_card("rate", 0.01, 0.001, 5.0);
    with_pressure_liquid(&source, true, |result, gate| {
        let mut liquid = result.unwrap();
        load(&mut liquid, gate);
        let native = liquid.flow().lattice().sample_at(3, 3, 3).unwrap();
        // Independent units: rho0=1000, dx/dt=1 and cs²=1/3.
        let pressure = REFERENCE + (native.density - 1.0) * 1000.0 / 3.0;
        let material = liquid.material_at(CELL).unwrap();
        assert!((coordinate(material, "pressure") - pressure).abs() < 3e-11);
        assert!((pressure - REFERENCE).abs() > 1e-5);
        let a = 0.5 + 3.0 * (0.1 + 0.01 * (pressure - REFERENCE));
        let rhs = stress_numerator(&liquid, CELL, [0.0; 3]) / 0.001;
        let rate = 2.0 * rhs / (a + (a * a + 0.012 * rhs).sqrt());
        let mu = 0.1 + 0.01 * (pressure - REFERENCE) + 0.001 * rate;
        near(
            liquid.rheology_at(CELL).unwrap().shear_rate_per_s,
            rate,
            1e-8,
        );
        near(material.dynamic_viscosity_pa_s(), mu, 2e-12);
        near(
            liquid.rheology_at(CELL).unwrap().collision_tau,
            0.5 + 3.0 * mu,
            1e-11,
        );
        let before = liquid.total_enthalpy_j().unwrap();
        let report = liquid.step_with_viscous_heating(&[], gate).unwrap();
        let applied = report
            .applied
            .iter()
            .find(|s| s.coordinates == CELL)
            .unwrap();
        assert!((coordinate(&applied.material, "pressure") - pressure).abs() < 3e-11);
        near(
            applied.heat_j,
            mu * rate * rate * applied.liquid_volume_m3 * 0.001,
            2e-8,
        );
        assert!(report.flow.viscous_heat_j > 0.0);
        assert!(
            (liquid.total_enthalpy_j().unwrap() - before - report.flow.viscous_heat_j).abs() < 1e-8
        );
        let endpoint = liquid.flow().lattice().sample_at(3, 3, 3).unwrap();
        let expected = REFERENCE + (endpoint.density - 1.0) * 1000.0 / 3.0;
        assert!(
            (coordinate(liquid.material_at(CELL).unwrap(), "pressure") - expected).abs() < 3e-11
        );
    });
}

#[test]
fn temperature_pressure_feedback_changes_real_flow_and_constant_limit_keeps_bits() {
    let exercise = |source: &MaterialCard, pressure| {
        with_pressure_liquid(source, pressure, |result, gate| {
            let mut liquid = result.unwrap();
            load(&mut liquid, gate);
            let heat = [FreeSurfaceCellHeat {
                coordinates: CELL,
                energy_j: 0.01,
            }];
            let report = liquid.step_with_viscous_heating(&heat, gate).unwrap();
            assert!(
                (report.flow.enthalpy_after_j
                    - report.flow.enthalpy_before_j
                    - 0.01
                    - report.flow.viscous_heat_j)
                    .abs()
                    < 1e-8
            );
            let endpoint = liquid.material_at(CELL).unwrap();
            assert!(coordinate(endpoint, "T") > 400.1);
            // Retained claims, actual endpoint temperature and actual native pressure.
            let p = coordinate(endpoint, "pressure");
            let t = coordinate(endpoint, "T");
            let mu = endpoint.dynamic_viscosity_pa_s();
            let next = liquid.step_with_viscous_heating(&[], gate).unwrap();
            (
                liquid.flow().lattice().grid().state_bits(),
                liquid.total_enthalpy_j().unwrap().to_bits(),
                next.flow.viscous_heat_j.to_bits(),
                p,
                t,
                mu,
            )
        })
    };
    let constant = pressure_card("T", 0.0, 0.0, 5.0);
    let fixed = exercise(&constant, false);
    let evolving = exercise(&constant, true);
    assert_eq!(
        (&fixed.0, fixed.1, fixed.2),
        (&evolving.0, evolving.1, evolving.2)
    );
    let thermal = exercise(&pressure_card("T", 0.0, -0.0002, 5.0), true);
    let coupled = exercise(&pressure_card("T", 0.01, -0.0002, 5.0), true);
    near(
        coupled.5,
        0.1 + 0.01 * (coupled.3 - REFERENCE) - 0.0002 * (coupled.4 - 400.0),
        2e-12,
    );
    assert_ne!(thermal.0, coupled.0);
    assert_ne!(thermal.2, coupled.2);
}

#[test]
fn pressure_caloric_axes_and_late_pressure_failure_refuse_atomically() {
    let original = pressure_card("T", 0.01, -0.0002, 5.0);
    for (name, value) in [
        ("density", 1000.0),
        ("specific_heat_capacity", 2000.0),
        ("surface_tension", 0.001),
    ] {
        let mut claims = ClaimSet::new();
        for (_, claim) in original.claims().claims_ordered() {
            let mut claim = claim.clone();
            if claim.key.name() == name {
                // Equal sampled values still carry a pressure-response axis.
                claim.value = PropertyValue::Curve {
                    abscissa: "pressure".into(),
                    abscissa_dims: Pressure::DIMS,
                    knots: vec![(REFERENCE - 5.0, value), (REFERENCE + 5.0, value)],
                    dims: claim.key.dims(),
                };
            }
            claims.insert_claim(claim).unwrap();
        }
        let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
        with_pressure_liquid(&card, true, |result, _| {
            assert!(matches!(
                result,
                Err(FreeSurfaceEnthalpyError::Invalid(
                    "thermal pressure rheology requires pressure-independent density, caloric storage and surface tension"
                ))
            ))
        });
    }
    let narrow = pressure_card("T", 0.01, -0.0002, 1e-5);
    with_pressure_liquid(&narrow, true, |result, gate| {
        let mut liquid = result.unwrap();
        let before = (
            liquid.flow().lattice().grid().state_bits(),
            liquid.total_enthalpy_j().unwrap(),
            liquid.material_at(CELL).unwrap().clone(),
            *liquid.rheology_at(CELL).unwrap(),
            liquid.physical_time_s(),
        );
        let heat = [FreeSurfaceCellHeat {
            coordinates: CELL,
            energy_j: 0.001,
        }];
        assert!(matches!(
            liquid.step_with_heat_and_gas_pressure(&heat, gate, 0.05, pressure_load),
            Err(FreeSurfaceEnthalpyError::Material(_))
        ));
        let cancelled = CancelGate::new_clock_free();
        let calls = std::cell::Cell::new(0);
        assert!(
            liquid
                .step_with_heat_and_gas_pressure(&heat, &cancelled, 0.05, |p| {
                    calls.set(calls.get() + 1);
                    if calls.get() == 5 {
                        cancelled.request();
                    }
                    pressure_load(p)
                })
                .is_err()
        );
        assert!(calls.get() >= 5);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), before.0);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), before.1);
        assert_eq!(liquid.material_at(CELL), Some(&before.2));
        assert_eq!(liquid.rheology_at(CELL), Some(&before.3));
        assert_eq!(liquid.physical_time_s(), before.4);
    });
}
