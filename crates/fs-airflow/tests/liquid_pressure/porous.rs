//! Source permeability through the native heat/pressure/solid owners.

mod forchheimer;
mod wet;

use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::porous::{
    LiquidResistanceLaw, PorousLiquidResistance,
};
use fs_material::porous::permeability::{
    DarcyFlowConstraint, DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

fn support() -> PorousGasSupport {
    PorousGasSupport {
        temperature_k: (250.0, 450.0),
        pressure_pa: (90000.0, 1e9),
        maximum_superficial_speed_m_s: 1.0,
    }
}
fn skeleton(viscous_geometry: f64, domain: PorousGasSupport) -> DarcyGasBarrier {
    let dims = Dims([2, 0, 0, 0, 0, 0]);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("intrinsic_permeability", dims),
            value: PropertyValue::Scalar {
                value: 0.01 / (viscous_geometry * area().value()),
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", domain.temperature_k.0, domain.temperature_k.1)
                .with("pressure", domain.pressure_pa.0, domain.pressure_pa.1),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic original porous skeleton; not measurements".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic porous solid".into(),
            phase: "solid".into(),
            process: "analytical permeability test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let property = ScalarPropertyRequirement::try_new(
        "intrinsic_permeability",
        dims,
        ScalarAdmissibility::NonNegative,
    )
    .unwrap();
    let source = resolve_material_state_point(
        &card,
        &point(300.0),
        &[property],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let layers = [0.003, 0.007]
        .into_iter()
        .map(|thickness| {
            PorousGasLayer::from_material_isotropic(
                source.clone(),
                "intrinsic_permeability",
                thickness,
                None,
            )
            .unwrap()
        })
        .collect();
    DarcyGasBarrier::try_new(layers, [1.0, 0.0, 0.0], domain).unwrap()
}
fn porous(barrier: &DarcyGasBarrier, outlet: f64) -> PorousLiquidResistance<'_> {
    PorousLiquidResistance::try_new(
        barrier.continuum_geometry().unwrap(),
        area(),
        Length::new(0.01),
        (0.0, 2300.0),
        Pressure::new(outlet),
    )
    .unwrap()
}
fn mixed(
    barrier: &DarcyGasBarrier,
) -> [LiquidResistanceSection<'static, LiquidResistanceLaw<'_>>; 2] {
    [
        LiquidResistanceSection {
            region: "first",
            resistance: LiquidResistanceLaw::Porous(porous(barrier, 100000.0)),
        },
        LiquidResistanceSection {
            region: "last",
            resistance: LiquidResistanceLaw::Geometry(resistances()[1].resistance),
        },
    ]
}

#[test]
fn the_original_liquid_pressure_solver_observes_anisotropic_confinement() {
    // G1/G3: the same rotated skeleton and real liquid law, with only the
    // authored transverse boundary condition changed. The separate ordinary
    // hydraulic section remains in series with the porous section.
    let k0 = 0.01 / (5e9 * area().value());
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let tensor = IntrinsicPermeability::orthotropic(
        [k0, 9.0 * k0, 0.0],
        [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    let layer = PorousGasLayer::try_new(
        tensor,
        0.01,
        0.0,
        "synthetic rotated fully developed plug, no entrance-flow claim",
    )
    .unwrap();
    let original = DarcyGasBarrier::try_new(vec![layer], [1.0, 0.0, 0.0], support()).unwrap();
    for (constraint, effective_k_over_k0) in [
        (DarcyFlowConstraint::NormalPressureGradient, 5.0),
        (DarcyFlowConstraint::ConfinedNormalFlow, 1.8),
    ] {
        let barrier = original.with_flow_constraint(constraint).unwrap();
        with_path(
            false,
            1.0,
            LiquidPropertyTemperature::BulkMean,
            |cx, path| {
                let result = path
                    .solve_pressure_driven(
                        cx,
                        &[300.0, 300.0],
                        &mixed(&barrier),
                        Pressure::new(50.0),
                        local(),
                        pressure_budget(),
                    )
                    .unwrap();
                let geometry = 5e9 / effective_k_over_k0 + 5e9;
                let expected = 50.0 * 1000.0 / (0.002 * geometry);
                close(result.mass_flow_kg_per_s, expected);
                let observed = result.hydraulics.regions[0].porous.as_ref().unwrap();
                assert_eq!(observed.barrier.flow_constraint(), constraint);
                assert_eq!(observed.barrier.layers(), original.layers());
                close(result.hydraulics.total_pressure_drop_pa, 50.0);
                close(result.thermal.total_heat_rate_w, 0.0);
            },
        );
    }
}

#[test]
fn only_permeability_and_fluid_data_change_the_native_operating_point() {
    for permeability_scale in [1.0, 2.0] {
        let barrier = skeleton(5e9 / permeability_scale, support());
        for viscosity_scale in [1.0, 2.0] {
            with_path(
                false,
                viscosity_scale,
                LiquidPropertyTemperature::BulkMean,
                |cx, path| {
                    let result = path
                        .solve_pressure_driven(
                            cx,
                            &[300.0, 300.0],
                            &mixed(&barrier),
                            Pressure::new(50.0),
                            local(),
                            pressure_budget(),
                        )
                        .unwrap();
                    let geometry = 5e9 / permeability_scale + 5e9;
                    let expected = 50.0 * 1000.0 / (0.002 * viscosity_scale * geometry);
                    close(result.mass_flow_kg_per_s, expected);
                    let first = result.hydraulics.regions[0].porous.as_ref().unwrap();
                    assert_eq!(first.barrier, barrier);
                    close(first.inlet_pressure_pa, 100050.0);
                    close(
                        first.outlet_pressure_pa,
                        100000.0 + result.hydraulics.regions[1].loss.pressure_drop_pa,
                    );
                    assert!(result.hydraulics.regions[1].porous.is_none());
                    for region in &result.hydraulics.regions {
                        close(region.loss.volume_flow_m3_s, expected / 1000.0);
                        close(
                            region.loss.dissipated_power_w,
                            region.loss.pressure_drop_pa * expected / 1000.0,
                        );
                    }
                },
            );
        }
    }
}

fn hot_reference(mass: f64, porous_geometry: f64) -> (f64, f64) {
    let e = (-183.0 * 0.003 / (mass * 2000.0)).exp();
    let first = 400.0 + (300.0 - 400.0) * e;
    let last = 370.0 + (first - 370.0) * e;
    let mu_first = 0.005 - 0.5e-5 * (300.0 + first);
    let mu_last = 0.005 - 0.5e-5 * (first + last);
    (
        (porous_geometry * mu_first + 5e9 * mu_last) * mass / 1000.0,
        last,
    )
}
#[test]
fn original_permeability_closes_temperature_viscosity_mass_and_enthalpy_together() {
    with_path(
        true,
        1.0,
        LiquidPropertyTemperature::FilmMean,
        |cx, path| {
            for factor in [5e9, 2.5e9] {
                let barrier = skeleton(factor, support());
                let mut lo = 0.0002;
                let mut hi = 0.008;
                for _ in 0..80 {
                    let m = 0.5 * (lo + hi);
                    if hot_reference(m, factor).0 < 50.0 {
                        lo = m;
                    } else {
                        hi = m;
                    }
                }
                let expected_mass = 0.5 * (lo + hi);
                let result = path
                    .solve_pressure_driven(
                        cx,
                        &[400.0, 370.0],
                        &mixed(&barrier),
                        Pressure::new(50.0),
                        local(),
                        pressure_budget(),
                    )
                    .unwrap();
                close(result.mass_flow_kg_per_s, expected_mass);
                close(
                    result.thermal.regions[1].outlet.temperature_k(),
                    hot_reference(expected_mass, factor).1,
                );
                close(
                    result.thermal.total_heat_rate_w,
                    expected_mass * 2000.0 * (hot_reference(expected_mass, factor).1 - 300.0),
                );
                assert_eq!(
                    result.hydraulics.regions[0]
                        .porous
                        .as_ref()
                        .unwrap()
                        .barrier,
                    barrier
                );
            }
        },
    );
}

#[test]
fn original_permeability_dissipation_enters_the_same_fluid_enthalpy() {
    let barrier = skeleton(5e15, support());
    let law = porous(&barrier, 100000.0);
    let rows = [
        LiquidResistanceSection {
            region: "first",
            resistance: law,
        },
        LiquidResistanceSection {
            region: "last",
            resistance: law,
        },
    ];
    with_path(
        true,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| {
            let report = path
                .march_with_dissipation(cx, &[300.0, 300.0], &rows, local())
                .unwrap();
            let expected = heated_reference(0.001, true);
            close(report.thermal.regions[1].outlet.temperature_k(), expected.0);
            close(report.hydraulics.total_pressure_drop_pa, expected.1);
            close(report.hydraulics.total_dissipated_power_w, expected.2);
            close(
                report.enthalpy_gain_w,
                report.thermal.total_heat_rate_w + expected.2,
            );
            close(report.balance_residual_w, 0.0);
            let first = report.hydraulics.regions[0].porous.as_ref().unwrap();
            let last = report.hydraulics.regions[1].porous.as_ref().unwrap();
            close(first.inlet_pressure_pa, 100000.0 + expected.1);
            assert_eq!(first.outlet_pressure_pa, last.inlet_pressure_pa);
            assert_eq!(last.outlet_pressure_pa, 100000.0);
        },
    );
}

#[test]
fn source_limits_apply_to_discarded_flow_trials_not_only_the_root() {
    with_path(
        false,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| {
            // The real root is supported (50 Pa, speed 0.032 m/s), but the upper
            // bracket trial is not (160 Pa, speed 0.102 m/s). Never extrapolate it.
            for domain in [
                PorousGasSupport {
                    maximum_superficial_speed_m_s: 0.05,
                    ..support()
                },
                PorousGasSupport {
                    pressure_pa: (90000.0, 100080.0),
                    ..support()
                },
            ] {
                let barrier = skeleton(5e9, domain);
                assert!(matches!(
                    path.solve_pressure_driven(
                        cx,
                        &[300.0, 300.0],
                        &mixed(&barrier),
                        Pressure::new(50.0),
                        local(),
                        pressure_budget(),
                    ),
                    Err(PressureDrivenError::Resistance { .. })
                ));
            }
            let barrier = skeleton(
                5e9,
                PorousGasSupport {
                    temperature_k: (250.0, 310.0),
                    ..support()
                },
            );
            assert!(matches!(
                path.solve_pressure_driven(
                    cx,
                    &[400.0, 400.0],
                    &mixed(&barrier),
                    Pressure::new(50.0),
                    local(),
                    pressure_budget(),
                ),
                Err(PressureDrivenError::Resistance { .. })
            ));
            let barrier = skeleton(5e9, support());
            let incompatible = [
                LiquidResistanceSection {
                    region: "first",
                    resistance: porous(&barrier, 100000.0),
                },
                LiquidResistanceSection {
                    region: "last",
                    resistance: porous(&barrier, 100001.0),
                },
            ];
            assert!(matches!(
                path.solve_pressure_driven(
                    cx,
                    &[300.0, 300.0],
                    &incompatible,
                    Pressure::new(50.0),
                    local(),
                    pressure_budget(),
                ),
                Err(PressureDrivenError::Invalid(_))
            ));
        },
    );
}

#[test]
fn sealed_original_permeability_cannot_be_regularized_into_flow() {
    let layer = PorousGasLayer::try_new(
        fs_material::porous::permeability::IntrinsicPermeability::isotropic(0.0).unwrap(),
        0.01,
        0.0,
        "synthetic impermeable layer",
    )
    .unwrap();
    let barrier = DarcyGasBarrier::try_new(vec![layer], [1.0, 0.0, 0.0], support()).unwrap();
    assert!(
        PorousLiquidResistance::try_new(
            barrier.continuum_geometry().unwrap(),
            area(),
            Length::new(0.01),
            (0.0, 2300.0),
            Pressure::new(100000.0),
        )
        .is_err()
    );
}

#[test]
fn coupled_solid_callback_receives_the_permeability_selected_fluid_state() {
    use fs_airflow::conjugate::SolidRegionState;
    use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
    let barrier = skeleton(2.5e9, support());
    with_path(
        false,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| {
            let mut calls = 0;
            let result = path
                .solve_pressure_coupled_iqn(
                    cx,
                    &[300.0, 300.0],
                    &mixed(&barrier),
                    Pressure::new(50.0),
                    &solid_controls(),
                    |_, rows| {
                        calls += 1;
                        Ok(LiquidSolidResponse {
                            state: calls,
                            regions: rows
                                .iter()
                                .enumerate()
                                .map(|(i, row)| SolidRegionState {
                                    region: ["first", "last"][i].into(),
                                    area_m2: 0.003,
                                    mean_wall_temperature_k: 300.0,
                                    heat_rate_w: 0.0,
                                    mean_reference_temperature_k: Some(
                                        row.reference_temperature().value(),
                                    ),
                                })
                                .collect(),
                        })
                    },
                )
                .unwrap();
            close(result.operating_point.mass_flow_kg_per_s, 50.0 / 15000.0);
            assert_eq!(calls, result.total_solid_solves);
            assert_eq!(
                result.operating_point.hydraulics.regions[0]
                    .porous
                    .as_ref()
                    .unwrap()
                    .barrier,
                barrier
            );
        },
    );
}
