//! Independent surface balance, finite-gas energy/pressure, stiffness and rollback.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceEmissivity, FreeSurfaceInterfaceContact, FreeSurfaceInterfaceFace,
    FreeSurfaceInterfaceSolve,
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};

fn face(surroundings_k: f64) -> FreeSurfaceInterfaceFace {
    FreeSurfaceInterfaceFace {
        face: FreeSurfaceGasFace {
            liquid: EXPOSED,
            gas: [4, 3, 3],
            area_m2: 1e-6,
            distances_m: [0.0005; 2],
        },
        emissivity: 0.8.into(),
        surroundings_k,
    }
}

fn gas_controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 690.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

// Synthetic k_g(T) = 0.001 T. The original source law remains its owner;
// this test fixture only changes data, not the surface or thermal algorithm.
fn temperature_dependent_gas() -> ThermochemicalGasModel {
    gas_with_conductivity(PropertyValue::Curve {
        abscissa: "temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(200.0, 0.2), (700.0, 0.7)],
        dims: Dims([1, 1, -3, -1, 0, 0]),
    })
}

fn gas_with_conductivity(value: PropertyValue) -> ThermochemicalGasModel {
    let base = gas_source::model(0.0, 250000.0, 1.0);
    let mut transport = base.transport_materials().to_vec();
    for species in &mut transport {
        let mut claims = ClaimSet::new();
        for id in species.material.claims().observation_ids() {
            claims
                .register_observation(species.material.claims().observation(id).unwrap().clone())
                .unwrap();
        }
        for (_, claim) in species.material.claims().claims_ordered() {
            let mut claim = claim.clone();
            if claim.key.name() == "thermal_conductivity" {
                claim.value = value.clone();
            }
            claims.insert_claim(claim).unwrap();
        }
        species.material = MaterialCard::assemble(
            species.material.id().clone(),
            claims,
            species.material.models().to_vec(),
        )
        .unwrap();
    }
    ThermochemicalGasModel::new(base.thermodynamics().clone(), transport, None).unwrap()
}

#[test]
fn g1_unsupported_film_trial_backtracks_with_gas_outside_liquid_phase_range() {
    let liquid = card(0.02, 1000.0);
    // At the actual 100 kPa, this triangle only supports film T <= 320 K,
    // despite the original source's enclosing validity box extending to 700 K.
    let gas = gas_with_conductivity(PropertyValue::Surface2D {
        axes: ["temperature".into(), "pressure".into()],
        axis_dims: [Temperature::DIMS, fs_qty::Pressure::DIMS],
        samples: vec![
            ([200.0, 1000.0], 0.26),
            ([440.0, 1000.0], 0.26),
            ([200.0, 199000.0], 0.26),
        ],
        triangles: vec![[0, 1, 2]],
        dims: Dims([1, 1, -3, -1, 0, 0]),
    });
    let context = gas_source::point(250.0, 100000.0);
    assert!(gas.evaluate(&context).is_ok());
    assert!(gas.evaluate(&gas_source::point(325.0, 100000.0)).is_err());
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let surface = FreeSurfaceInterfaceFace {
            emissivity: 0.0.into(),
            ..face(400.0)
        };
        let report = thermal
            .step_with_surface_exchange(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &[FreeSurfaceInterfaceContact {
                    surface,
                    model: &gas,
                    context: &context,
                }],
                FreeSurfaceInterfaceSolve::default(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let endpoint = &report.surface.as_ref().unwrap().thermal_endpoint[0];
        // Initial Ts=400 has unsupported film T=325, and bulk gas 250 K is
        // outside the liquid's [300,600] K range. The physical interior root
        // is still supported and must not be refused by either trial choice.
        let exact = (0.5 * endpoint.liquid_temperature_k + 0.26 * 250.0) / (0.5 + 0.26);
        near(endpoint.surface_temperature_k, exact, 2e-11);
        assert!(endpoint.surface_temperature_k > 300.0);
        assert!(endpoint.gas.state().temperature < 320.0);
        assert!(report.surface.as_ref().unwrap().heat_from_gas_j < 0.0);
    });
}

// Independent bisection of actual W: no production conductance/root helper.
fn reference_surface(
    tl: f64,
    tg: f64,
    tr: f64,
    emissivity: f64,
    kl: f64,
    conductivity: impl Fn(f64) -> f64,
) -> (f64, f64, f64) {
    let (mut lo, mut hi) = (tl.min(tg).min(tr), tl.max(tg).max(tr));
    for _ in 0..80 {
        let ts = 0.5 * (lo + hi);
        let liquid = kl * 1e-6 / 0.0005 * (tl - ts);
        let gas = conductivity(0.5 * (ts + tg)) * 1e-6 / 0.0005 * (tg - ts);
        let radiant = emissivity * 5.670374419e-8 * 1e-6 * (tr.powi(4) - ts.powi(4));
        if liquid + gas + radiant > 0.0 {
            lo = ts;
        } else {
            hi = ts;
        }
    }
    let ts = 0.5 * (lo + hi);
    (
        ts,
        conductivity(0.5 * (ts + tg)) * 1e-6 / 0.0005 * (tg - ts),
        emissivity * 5.670374419e-8 * 1e-6 * (tr.powi(4) - ts.powi(4)),
    )
}

// Independent source-dependent nonlinear boundary, using fourth powers rather
// than the production secant/Schur equations.
fn reference_material_surface(tl: f64, tg: f64, tr: f64) -> (f64, f64, f64) {
    let (mut lo, mut hi) = (300.0_f64, 500.0_f64);
    for _ in 0..80 {
        let ts = 0.5 * (lo + hi);
        let eps = 0.2 + 0.003 * (ts - 300.0);
        let ql = 0.0005 * (tl - 200.0) * 1e-6 / 0.0005 * (tl - ts);
        let qg = 0.001 * 0.5 * (ts + tg) * 1e-6 / 0.0005 * (tg - ts);
        let qr = eps * 5.670374419e-8 * 1e-6 * (tr.powi(4) - ts.powi(4));
        if ql + qg + qr > 0.0 {
            lo = ts;
        } else {
            hi = ts;
        }
    }
    let ts = 0.5 * (lo + hi);
    (
        ts,
        0.001 * 0.5 * (ts + tg) * 1e-6 / 0.0005 * (tg - ts),
        (0.2 + 0.003 * (ts - 300.0)) * 5.670374419e-8 * 1e-6 * (tr.powi(4) - ts.powi(4)),
    )
}

#[test]
fn g1_original_surface_emissivity_drives_finite_gas_pressure_and_liquid_heat() {
    use fs_material::thermal::source::SURFACE_EMISSIVITY_PROPERTY;
    let liquid = radiation::emissive_card(1.0, 500.0);
    let law = radiation::emissive_law(&liquid);
    let chart = caloric(&liquid);
    with_exposed_law(&liquid, false, &chart, 200000.0, &law, |thermal, gate| {
        let mut gas = IsochoricGasParcel::new(
            temperature_dependent_gas(),
            gas_source::point(350.0, 100000.0),
            1e-6,
        )
        .unwrap();
        let gas_before = gas.state().specific_internal_energy_j_kg();
        let before = thermal.total_enthalpy_j().unwrap();
        let surface = FreeSurfaceInterfaceFace {
            emissivity: FreeSurfaceEmissivity::Material,
            ..face(800.0)
        };
        let report = thermal
            .step_with_finite_gas_surface_exchange(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &mut gas,
                &[surface],
                gas_controls(),
                FreeSurfaceInterfaceSolve::default(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let r = 8.314_462_618_153_24 / 0.028;
        let mass = 100000.0 * 1e-6 / (r * 350.0);
        let (mut tl, mut tg, mut qg_j, mut qr_j) = (400.0_f64, 350.0_f64, 0.0, 0.0);
        for _ in 0..4 {
            let (_, qg, qr) = reference_material_surface(tl, tg, 800.0);
            let heat_from_gas = 0.00025 * qg;
            tl += 0.00025 * (qg + qr) / 0.002;
            tg = (((2.0 + 0.001 * tg).powi(2) - 0.002 * heat_from_gas / (mass * r)).sqrt() - 2.0)
                / 0.001;
            qg_j += heat_from_gas;
            qr_j += 0.00025 * qr;
        }
        let exchange = report.surface.as_ref().unwrap();
        let endpoint = &exchange.thermal_endpoint[0];
        let (ts, qg, qr) = reference_material_surface(tl, tg, 800.0);
        assert_eq!(report.thermal.substeps, 4);
        near(endpoint.liquid_temperature_k, tl, 2e-11);
        near(endpoint.surface_temperature_k, ts, 2e-11);
        assert!(
            (ts - tl).abs() > 1.0,
            "surface must have its own physical T"
        );
        near(endpoint.emissivity, 0.2 + 0.003 * (ts - 300.0), 2e-11);
        let source = endpoint.emissivity_source.as_ref().unwrap();
        assert_eq!(
            source.answer().receipt.selected,
            law.reference_state()
                .material()
                .property(SURFACE_EMISSIVITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
        near(
            source
                .answer()
                .receipt
                .query_point
                .iter()
                .find(|(axis, _)| axis == "T")
                .unwrap()
                .1,
            ts,
            2e-11,
        );
        assert_eq!(
            source
                .answer()
                .receipt
                .query_point
                .iter()
                .find(|(axis, _)| axis == "pressure")
                .unwrap()
                .1,
            100000.0
        );
        near(endpoint.heat_from_gas_w, qg, 2e-8);
        near(endpoint.heat_from_radiation_w, qr, 2e-8);
        near(exchange.heat_from_gas_j, qg_j, 2e-8);
        near(exchange.heat_from_radiation_j, qr_j, 2e-8);
        near(gas.state().state().temperature, tg, 2e-11);
        near(gas.state().state().pressure, 100000.0 * tg / 350.0, 2e-11);
        let gas_change = gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - gas_before);
        near(
            thermal.total_enthalpy_j().unwrap() - before + gas_change,
            qr_j,
            2e-8,
        );
        assert!(
            report
                .finite_gas
                .unwrap()
                .combined_thermal_balance_defect_j
                .abs()
                < 1e-9
        );
    });
}

#[test]
fn g1_phase_conductivity_reaches_solved_surface_and_signed_energy_balance() {
    use fs_material::phase::conductivity::LiquidMassFractionConductivity;
    let liquid = conducting_card(0.0005);
    let gas = temperature_dependent_gas();
    let context = gas_source::point(350.0, 100000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        *thermal = thermal
            .clone()
            .with_phase_conductivity(
                LiquidMassFractionConductivity::declared(4.0, 3.0).unwrap(),
                gate,
            )
            .unwrap();
        let before = thermal.total_enthalpy_j().unwrap();
        let report = thermal
            .step_with_surface_exchange(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &[FreeSurfaceInterfaceContact {
                    surface: face(800.0),
                    model: &gas,
                    context: &context,
                }],
                FreeSurfaceInterfaceSolve::default(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let (mut t, mut gas_j, mut radiant_j) = (400.0, 0.0, 0.0);
        for _ in 0..4 {
            let (_, qg, qr) =
                reference_surface(t, 350.0, 800.0, 0.8, 3.0 * 0.0005 * (t - 200.0), |film| {
                    0.001 * film
                });
            t += 0.00025 * (qg + qr) / 0.002;
            gas_j += 0.00025 * qg;
            radiant_j += 0.00025 * qr;
        }
        let (ts, qg, qr) =
            reference_surface(t, 350.0, 800.0, 0.8, 3.0 * 0.0005 * (t - 200.0), |film| {
                0.001 * film
            });
        let exchange = report.surface.as_ref().unwrap();
        let endpoint = &exchange.thermal_endpoint[0];
        assert_eq!(report.thermal.substeps, 4);
        near(endpoint.liquid_temperature_k, t, 2e-12);
        near(endpoint.surface_temperature_k, ts, 2e-11);
        near(
            endpoint.liquid_conductivity.source_w_m_k,
            0.0005 * (t - 200.0),
            2e-12,
        );
        assert_eq!(endpoint.liquid_conductivity.phase_multiplier, 3.0);
        near(
            endpoint.liquid_conductance_w_k,
            3.0 * 0.0005 * (t - 200.0) * 1e-6 / 0.0005,
            2e-12,
        );
        near(endpoint.heat_from_gas_w, qg, 2e-8);
        near(endpoint.heat_from_radiation_w, qr, 2e-8);
        near(exchange.heat_from_gas_j, gas_j, 2e-8);
        near(exchange.heat_from_radiation_j, radiant_j, 2e-8);
        near(report.flow.external_heat_j, gas_j + radiant_j, 2e-8);
        near(
            thermal.total_enthalpy_j().unwrap(),
            before + gas_j + radiant_j,
            2e-12,
        );
        assert!(endpoint.balance_residual_w.abs() < 1e-10);
    });
}

#[test]
fn g1_shared_surface_matches_independent_nonlinear_balance_and_film_source() {
    let liquid = conducting_card(0.0005);
    let gas = temperature_dependent_gas();
    let context = gas_source::point(350.0, 100000.0);
    for emissivity in [0.0, 0.8] {
        with_exposed_liquid(&liquid, false, |thermal, gate| {
            let surface = FreeSurfaceInterfaceFace {
                emissivity: emissivity.into(),
                ..face(800.0)
            };
            let before = thermal.total_enthalpy_j().unwrap();
            let report = thermal
                .step_with_surface_exchange(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &[FreeSurfaceInterfaceContact {
                        surface,
                        model: &gas,
                        context: &context,
                    }],
                    FreeSurfaceInterfaceSolve::default(),
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            let exchange = report.surface.as_ref().unwrap();
            let (mut tl, mut gas_j, mut radiant_j) = (400.0, 0.0, 0.0);
            for _ in 0..4 {
                let (_, qg, qr) = reference_surface(
                    tl,
                    350.0,
                    800.0,
                    emissivity,
                    0.0005 * (tl - 200.0),
                    |film| 0.001 * film,
                );
                tl += 0.00025 * (qg + qr) / 0.002;
                gas_j += 0.00025 * qg;
                radiant_j += 0.00025 * qr;
            }
            let (ts, qg, qr) = reference_surface(
                tl,
                350.0,
                800.0,
                emissivity,
                0.0005 * (tl - 200.0),
                |film| 0.001 * film,
            );
            let endpoint = &exchange.thermal_endpoint[0];
            assert_eq!(report.thermal.substeps, 4);
            near(endpoint.liquid_temperature_k, tl, 2e-12);
            near(endpoint.surface_temperature_k, ts, 2e-11);
            assert!(
                tl - ts > 10.0,
                "finite liquid resistance must change surface T"
            );
            near(endpoint.gas.state().temperature, 0.5 * (ts + 350.0), 2e-11);
            near(
                endpoint.gas.state().thermal_conductivity,
                0.001 * 0.5 * (ts + 350.0),
                2e-11,
            );
            near(endpoint.heat_from_gas_w, qg, 2e-8);
            near(endpoint.heat_from_radiation_w, qr, 2e-8);
            near(endpoint.heat_into_liquid_w, qg + qr, 2e-8);
            assert!(endpoint.balance_residual_w.abs() < 1e-10);
            near(exchange.heat_from_gas_j, gas_j, 2e-8);
            near(exchange.heat_from_radiation_j, radiant_j, 2e-8);
            near(report.flow.external_heat_j, gas_j + radiant_j, 2e-8);
            near(report.thermal.external_heat_j, gas_j + radiant_j, 2e-8);
            near(
                thermal.total_enthalpy_j().unwrap(),
                before + gas_j + radiant_j,
                2e-12,
            );
            assert!(report.environment.is_none() && report.radiation.is_none());
        });
    }
}

#[test]
fn g1_finite_gas_receives_radiant_heat_through_surface_and_matches_nasa_storage() {
    let liquid = conducting_card(0.01);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut gas = IsochoricGasParcel::new(
            gas_source::model(100000.0, 250000.0, 100.0),
            gas_source::point(400.0, 100000.0),
            1e-8,
        )
        .unwrap();
        let gas_before = gas.state().specific_internal_energy_j_kg();
        let before = thermal.total_enthalpy_j().unwrap();
        let report = thermal
            .step_with_finite_gas_surface_exchange(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &mut gas,
                &[face(1000.0)],
                gas_controls(),
                FreeSurfaceInterfaceSolve::default(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.as_ref().unwrap();
        let finite = report.finite_gas.as_ref().unwrap();
        let r = 8.314_462_618_153_24 / 0.028;
        let mass = 100000.0 * 1e-8 / (r * 400.0);
        let (mut tl, mut tg, mut gas_j, mut radiant_j) = (400.0, 400.0, 0.0, 0.0);
        for _ in 0..4 {
            let pressure = 100000.0 * tg / 400.0;
            let kg = 2.6 * (0.5 + (pressure - 1000.0) / 249000.0);
            let (_, qg, qr) = reference_surface(tl, tg, 1000.0, 0.8, 0.01 * (tl - 200.0), |_| kg);
            let gas_heat = 0.00025 * qg;
            tl += 0.00025 * (qg + qr) / 0.002;
            // cv/R=2+0.001T, integrating NASA u before independent inversion.
            tg = (((2.0_f64 + 0.001 * tg).powi(2) - 0.002 * gas_heat / (mass * r)).sqrt() - 2.0)
                / 0.001;
            gas_j += gas_heat;
            radiant_j += 0.00025 * qr;
        }
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(finite.substeps, 4);
        assert!(
            exchange.heat_from_gas_j < 0.0,
            "radiation must heat the initially equal-T gas"
        );
        assert!(gas.state().state().temperature > 400.0);
        near(exchange.thermal_endpoint[0].liquid_temperature_k, tl, 2e-11);
        near(gas.state().state().temperature, tg, 2e-11);
        near(gas.state().state().pressure, 100000.0 * tg / 400.0, 2e-11);
        near(exchange.heat_from_gas_j, gas_j, 2e-8);
        near(exchange.heat_from_radiation_j, radiant_j, 2e-8);
        let gas_change = gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - gas_before);
        assert!((gas_change + exchange.heat_from_gas_j).abs() <= finite.gas_energy_tolerance_j);
        near(
            thermal.total_enthalpy_j().unwrap() - before + gas_change,
            exchange.heat_from_radiation_j,
            2e-8,
        );
        near(
            report.flow.external_heat_j,
            exchange.heat_from_gas_j + exchange.heat_from_radiation_j,
            2e-12,
        );
        assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g1_small_gas_capacity_limits_radiative_path_as_well_as_liquid_exchange() {
    let liquid = conducting_card(1e-6);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut gas = IsochoricGasParcel::new(
            gas_source::model(100000.0, 250000.0, 100.0),
            gas_source::point(400.0, 100000.0),
            1e-12,
        )
        .unwrap();
        let report = thermal
            .step_with_finite_gas_surface_exchange(
                &[],
                &[],
                FreeSurfaceConductionControls {
                    maximum_dt_s: 0.001,
                    ..controls()
                },
                &mut [],
                &[],
                &[],
                &mut gas,
                &[face(450.0)],
                gas_controls(),
                FreeSurfaceInterfaceSolve::default(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        assert!(
            report.thermal.substeps > 10,
            "finite gas/radiant decay must shorten the same thermal interval"
        );
        let surface = &report.surface.as_ref().unwrap().thermal_endpoint[0];
        assert!(surface.gas_decay_conductance_w_k > 10.0 * surface.liquid_decay_conductance_w_k);
        assert!(
            gas.state().state().temperature > 440.0 && gas.state().state().temperature <= 450.0
        );
        assert!((400.0..=450.0).contains(&surface.liquid_temperature_k));
        assert!(
            report
                .finite_gas
                .unwrap()
                .combined_thermal_balance_defect_j
                .abs()
                < 1e-9
        );
    });
}

#[test]
fn g3_hot_radiant_reservoir_outside_gas_domain_changes_native_heat_and_motion() {
    let liquid = conducting_card(1.0);
    let model = gas_source::model_with_temperature_limit(0.0, 250000.0, 100.0, 450.0);
    let context = gas_source::point(400.0, 100000.0);
    let run = |emissivity: f64| {
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let surface = FreeSurfaceInterfaceFace {
                emissivity: emissivity.into(),
                ..face(2000.0)
            };
            let report = thermal
                .step_with_surface_exchange(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &[FreeSurfaceInterfaceContact {
                        surface,
                        model: &model,
                        context: &context,
                    }],
                    FreeSurfaceInterfaceSolve::default(),
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            let exchange = report.surface.unwrap();
            assert!(exchange.thermal_endpoint[0].gas.state().temperature < 450.0);
            assert!(exchange.thermal_endpoint[0].surface_temperature_k < 450.0);
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal
                    .material_at(EXPOSED)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                exchange.heat_from_gas_j + exchange.heat_from_radiation_j,
            )
        })
    };
    let insulated = run(0.0);
    let heated = run(0.8);
    assert_eq!(insulated.2, 0.0);
    assert!(heated.2 > 0.0 && heated.1 < insulated.1);
    assert_ne!(heated.0, insulated.0);
}

#[test]
fn g4_surface_budget_and_late_unsupported_root_preserve_liquid_wall_and_gas() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut gas = IsochoricGasParcel::new(
            gas_source::model_with_temperature_limit(100000.0, 250000.0, 100.0, 410.0),
            gas_source::point(400.0, 100000.0),
            1e-7,
        )
        .unwrap();
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let bits = thermal.flow().lattice().grid().state_bits();
        let before = thermal.total_enthalpy_j().unwrap();
        let liquid_before = thermal.material_at(EXPOSED).unwrap().clone();
        let wall_before = bodies[0].specific_enthalpy_j_kg();
        let wall_source = bodies[0].material().clone();
        let gas_before = gas.state().clone();
        let gas_point = gas.point().clone();
        let bounded_gas = GasEnergySolve {
            maximum_temperature_k: 410.0,
            ..gas_controls()
        };
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        for failure in 0..3 {
            let heat = if failure == 1 {
                // The first 15 K substep is supported. Subsequent bulk liquid
                // heating needs a film above 410 K: a later root must refuse.
                vec![FreeSurfaceCellHeat {
                    coordinates: EXPOSED,
                    energy_j: 0.12,
                }]
            } else {
                Vec::new()
            };
            let solve = FreeSurfaceInterfaceSolve {
                maximum_iterations: if failure == 0 { 1 } else { 64 },
                ..FreeSurfaceInterfaceSolve::default()
            };
            assert!(
                thermal
                    .step_with_finite_gas_surface_exchange(
                        &[],
                        &heat,
                        controls(),
                        &mut bodies,
                        &[body::contact()],
                        &[],
                        &mut gas,
                        &[face(800.0)],
                        bounded_gas,
                        solve,
                        true,
                        if failure == 2 { &cancelled } else { gate },
                    )
                    .is_err()
            );
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
            assert_eq!(thermal.material_at(EXPOSED).unwrap(), &liquid_before);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
            assert_eq!(bodies[0].material(), &wall_source);
            assert_eq!(gas.state(), &gas_before);
            assert_eq!(gas.point(), &gas_point);
        }
        thermal
            .step_with_finite_gas_surface_exchange(
                &[],
                &[],
                controls(),
                &mut bodies,
                &[body::contact()],
                &[],
                &mut gas,
                &[face(800.0)],
                bounded_gas,
                FreeSurfaceInterfaceSolve::default(),
                true,
                gate,
            )
            .unwrap();
    });
}
