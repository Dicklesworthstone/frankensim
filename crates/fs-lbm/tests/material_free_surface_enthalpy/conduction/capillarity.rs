//! Constant sourced capillarity on the original thermal, rheology and VOF owners.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceGasFace;
use fs_lbm::material::free_surface::enthalpy::rheology::RheologySolve;
use fs_material::liquid::surface::{
    ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
};
use fs_material::state_point::{
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_qty::Frequency;

#[path = "../../nsf3_thermochemical/source.rs"]
mod gas_source;

#[path = "wetting.rs"]
mod wetting;

#[path = "electrical.rs"]
mod electrical;

const EXPOSED: [usize; 3] = [3, 3, 3];

fn source_card(sigma: f64, shear: bool, sigma_upper_k: f64) -> MaterialCard {
    let base = conducting_card(1.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        let name = claim.key.name();
        if name == SURFACE_TENSION_PROPERTY {
            claim.value = PropertyValue::Scalar { value: sigma, dims: SURFACE_TENSION_DIMS };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
            claim.validity = ValidityDomain::unconstrained()
                .with("T", 300.0, sigma_upper_k)
                .with("pressure", 100000.0, 100000.0);
            if sigma_upper_k < 500.0 {
                let mut fallback = claim.clone();
                fallback.value = PropertyValue::Scalar {
                    value: 2.0 * sigma, dims: SURFACE_TENSION_DIMS,
                };
                fallback.validity = ValidityDomain::unconstrained()
                    .with("T", sigma_upper_k + 0.000001, 500.0)
                    .with("pressure", 100000.0, 100000.0);
                claims.insert_claim(fallback).unwrap();
            }
        } else if name == "dynamic_viscosity" {
            let dims = claim.key.dims();
            if shear {
                claim.value = PropertyValue::Surface2D {
                    axes: ["T".into(), "rate".into()],
                    axis_dims: [Temperature::DIMS, Frequency::DIMS],
                    samples: [[300.0, 0.0], [500.0, 0.0], [300.0, 100.0], [500.0, 100.0]]
                        .into_iter()
                        .map(|p| (p, 0.1 - 0.0002 * (p[0] - 400.0) + 0.002 * p[1]))
                        .collect(),
                    triangles: vec![[0, 1, 2], [1, 3, 2]],
                    dims,
                };
                claim.validity = claim.validity.with_quantity(
                    "rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0, 100.0,
                );
            } else {
                claim.value = PropertyValue::Scalar { value: 0.1, dims };
                claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn bundle(card: &MaterialCard, shear: bool, t: f64) -> ResolvedMaterialStatePoint {
    let mut requirements = [
        ("density", Dims([-3, 1, 0, 0, 0, 0])),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0])),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS),
    ].into_iter().map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    }).collect::<Vec<_>>();
    requirements.push(ScalarPropertyRequirement::try_with_quantity(
        "thermal_conductivity",
        QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static)),
        ScalarAdmissibility::StrictlyPositive,
    ).unwrap());
    let query = if shear {
        point(t).with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0).unwrap()
    } else {
        point(t)
    };
    resolve_material_state_point(
        card, &query, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

fn liquid_source(card: &MaterialCard, shear: bool) -> ResolvedLiquidSurfaceState {
    let material = bundle(card, shear, 400.0);
    if shear {
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity(
            material, LiquidViscosityBasis::Dynamic, LiquidCapillarity::Sourced, "T", "rate",
        )
    } else {
        ResolvedLiquidSurfaceState::try_from_material(
            material, LiquidViscosityBasis::Dynamic, LiquidCapillarity::Sourced,
        )
    }.unwrap()
}

fn cells() -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut cells = Vec::new();
    for z in 1..=6 {
        for y in 1..=6 {
            for x in 0..=3 {
                cells.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    velocity_m_s: [0.0, 0.003 * (x as f64 - 1.5), 0.0],
                });
            }
        }
    }
    cells
}

fn with_capillary<R>(
    card: &MaterialCard,
    shear: bool,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(card, shear);
    let law = LiquidTemperatureLaw::try_from_material(
        card, source.material().clone(), "T", LiquidViscosityBasis::Dynamic,
    ).unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        let mut thermal = if shear {
            scaling.free_surface3_from_cells_with_rheology(
                domain(), &cells(), &law, &chart, FreeSurfaceEnthalpyControls::default(),
                RheologySolve {
                    maximum_iterations: 128,
                    relaxation: 0.5,
                    absolute_relaxation_tolerance: 1e-13,
                    relative_relaxation_tolerance: 1e-13,
                },
                |_| 200000.0, cx, &gate,
            )
        } else {
            scaling.free_surface3_from_cells(domain(), &cells(), cx).unwrap()
                .with_enthalpy(
                    &law, &chart, FreeSurfaceEnthalpyControls::default(), |_| 200000.0, &gate,
                )
        }.unwrap();
        run(&mut thermal, &gate)
    })
}

fn gas_controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 690.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

fn gas(conductivity: f64) -> IsochoricGasParcel {
    IsochoricGasParcel::new(
        gas_source::model(100000.0, 250000.0, conductivity),
        gas_source::point(450.0, 100000.0),
        1e-8,
    ).unwrap()
}

fn gas_face() -> FreeSurfaceGasFace {
    FreeSurfaceGasFace {
        liquid: EXPOSED, gas: [4, 3, 3], area_m2: 1e-6, distances_m: [0.0005; 2],
    }
}

#[test]
fn g1_constant_sigma_keeps_native_capillary_motion_and_uniform_enthalpy() {
    let mut outcomes = Vec::new();
    for sigma in [0.001, 0.004] {
        let card = source_card(sigma, false, 500.0);
        let source = liquid_source(&card, false);
        let gate = CancelGate::new_clock_free();
        let native = with_cx(&gate, |cx| {
            let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
            near(scaling.surface_tension_lattice(), sigma * 0.001_f64.powi(2)
                / (1000.0 * 0.001_f64.powi(3)), 2e-14);
            let mut flow = scaling.free_surface3_from_cells(domain(), &cells(), cx).unwrap();
            flow.step_with_gate(&gate).unwrap();
            flow.lattice().grid().state_bits()
        });
        outcomes.push(with_capillary(&card, false, |thermal, gate| {
            let before = thermal.total_enthalpy_j().unwrap();
            let report = thermal.step(gate).unwrap();
            let bits = thermal.flow().lattice().grid().state_bits();
            assert_eq!(bits, native, "thermal binding must retain native capillary dynamics");
            near(report.mass_after_kg, report.mass_before_kg, 2e-13);
            near(thermal.total_enthalpy_j().unwrap(), before, 2e-11);
            assert_eq!(report.external_heat_j, 0.0);
            assert_eq!(report.viscous_heat_j, 0.0);
            assert!(report.balance_defect_j.abs() < 1e-9);
            let state = thermal.material_at(EXPOSED).unwrap();
            near(temperature(state), 400.0, 2e-11);
            assert_eq!(state.material().properties().len(), 5);
            let retained = state.material().property(SURFACE_TENSION_PROPERTY).unwrap();
            assert_eq!(retained.answer().receipt.selected,
                source.material().property(SURFACE_TENSION_PROPERTY).unwrap().answer().receipt.selected);
            assert_eq!(retained.value_si(), sigma);
            bits
        }));
    }
    assert_ne!(outcomes[0], outcomes[1], "sigma data must change the actual native populations");
}

#[test]
fn g1_g3_capillarity_rheology_finite_gas_and_wall_share_the_native_heat_trajectory() {
    let run = |sigma, conductivity| {
        let card = source_card(sigma, true, 500.0);
        let solid = body::wall_card(0.5, 500.0);
        let chart = body::wall_chart(&solid);
        let context = point(450.0);
        with_capillary(&card, true, |thermal, gate| {
            let mut gas = gas(conductivity);
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
            let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j()
                + gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
            let report = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[body::contact()], &[],
                &mut gas, &[gas_face()], gas_controls(), false, gate,
            ).unwrap().conduction;
            let environment = report.environment.as_ref().unwrap();
            let wall = report.bodies.as_ref().unwrap();
            let finite = report.finite_gas.as_ref().unwrap();

            // Independent four-step reference. Geometry and storage are authored:
            // liquid 1e-6 kg * 2000 J/kg/K; wall 2e-6 kg * 500 J/kg/K.
            // NASA9 cv/R=2+0.001T gives a quadratic energy inversion. Preserve
            // the original gas formation datum: exchange changes u, not its zero.
            let r = 8.314_462_618_153_24 / 0.028;
            let mg = 100000.0 * 1e-8 / (r * 450.0);
            let (mut tl, mut tg, mut twl, mut tw) = (400.0, 450.0, 400.0, 450.0);
            let (mut qg, mut qw) = (0.0, 0.0);
            for _ in 0..4 {
                let pressure = 100000.0 * tg / 450.0;
                let kg = 0.026 * conductivity * (0.5 + (pressure - 1000.0) / 249000.0);
                let g = 1e-6 / (0.0005 / (tl - 200.0) + 0.0005 / kg);
                let q = 0.00025 * g * (tg - tl);
                tl += q / 0.002;
                tg = (((2.0_f64 + 0.001 * tg).powi(2) - 0.002 * q / (mg * r)).sqrt()
                    - 2.0) / 0.001;
                qg += q;
                let g = 1e-6 / (0.0005 / (twl - 200.0) + 0.0005 / (0.5 * (tw - 200.0)));
                let q = 0.00025 * g * (tw - twl);
                twl += q / 0.002;
                tw -= q / 0.001;
                qw += q;
            }
            assert_eq!(report.thermal.substeps, 4);
            assert_eq!(finite.substeps, 4);
            near(environment.thermal_endpoint[0].liquid_temperature_k, tl, 2e-11);
            near(temperature(&wall.thermal_endpoint[0].liquid), twl, 2e-11);
            near(gas.state().state().temperature, tg, 2e-11);
            near(gas.state().state().pressure, 100000.0 * tg / 450.0, 2e-11);
            near(environment.heat_to_liquid_j, qg, 2e-11);
            near(wall.heat_to_liquid_j, qw, 2e-11);
            near(report.flow.external_heat_j, qg, 2e-11);
            assert_eq!(report.flow.viscous_heat_j, 0.0);
            let after = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j()
                + gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
            assert!((after - before).abs() < 1e-9);
            assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
            near(environment.thermal_endpoint[0].gas.state().pressure,
                gas.state().state().pressure, 2e-13);

            let material = thermal.material_at(EXPOSED).unwrap();
            let sample = thermal.rheology_at(EXPOSED).unwrap();
            near(material.dynamic_viscosity_pa_s(),
                0.1 - 0.0002 * (temperature(material) - 400.0)
                    + 0.002 * sample.shear_rate_per_s, 2e-12);
            near(sample.collision_tau, 0.5 + 3.0 * material.dynamic_viscosity_pa_s(), 2e-12);
            assert_eq!(material.material().property(SURFACE_TENSION_PROPERTY).unwrap().value_si(), sigma);
            assert_eq!(material.material().properties().len(), 5);

            // The first native step establishes nonequilibrium stress. Let the
            // heated constitutive state act on that stress without reusing an
            // exposed-face geometry after VOF has advanced.
            thermal.step(gate).unwrap();
            (thermal.flow().lattice().grid().state_bits(), qg)
        })
    };
    let slow = run(0.001, 10.0);
    let fast = run(0.001, 100.0);
    let stronger = run(0.004, 10.0);
    assert!(fast.1 > slow.1);
    assert_eq!(slow.0.len(), fast.0.len());
    let delta = slow.0.iter().zip(&fast.0)
        .map(|(&a, &b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
        .fold(0.0_f64, f64::max);
    assert!(delta > 1e-13, "gas data must alter resolved native motion, not only receipts");
    assert_ne!(slow.0, stronger.0, "sigma data must retain its native mechanical effect");
}

#[test]
fn g0_full_bundle_and_constant_source_admission_refuse_missing_or_inferred_capillarity() {
    let original = source_card(0.001, false, 500.0);
    let full = bundle(&original, false, 400.0);
    let law = LiquidTemperatureLaw::try_from_material(
        &original, full.clone(), "T", LiquidViscosityBasis::Dynamic,
    ).unwrap();
    assert_eq!(law.reference_state().material(), &full);
    let warm = law.at_temperature(Temperature::new(450.0)).unwrap();
    for (a, b) in full.properties().iter().zip(warm.material().properties()) {
        assert_eq!(a.requirement(), b.requirement());
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    let other = source_card(0.004, false, 500.0);
    assert!(LiquidTemperatureLaw::try_from_material(
        &other, full.clone(), "T", LiquidViscosityBasis::Dynamic,
    ).is_err());
    assert!(LiquidTemperatureLaw::try_from_material(
        &original, full, "T", LiquidViscosityBasis::Kinematic,
    ).is_err());

    for curve in [true, false] {
        // Equal-valued curve samples are not a source-declared constant law.
        let card = if curve { super::super::card(0.1, 1000.0) } else {
            source_card(0.001, false, 500.0)
        };
        let source = liquid_source(&card, false);
        let law = if curve {
            LiquidTemperatureLaw::try_from_material(
                &card, source.material().clone(), "T", LiquidViscosityBasis::Dynamic,
            )
        } else {
            // The old four-property resolver cannot silently omit active sigma.
            LiquidTemperatureLaw::try_new(
                &card, &point(400.0), "T", LiquidViscosityBasis::Dynamic,
                MaterialPropertySelection::SingleClaimOnly,
            )
        }.unwrap();
        let chart = caloric(&card);
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap()
                .free_surface3_from_cells(domain(), &cells(), cx).unwrap();
            assert!(flow.with_enthalpy(
                &law, &chart, FreeSurfaceEnthalpyControls::default(), |_| 200000.0, &gate,
            ).is_err());
        });
    }
}

#[test]
fn g4_late_sigma_domain_exit_rolls_back_liquid_finite_wall_and_atmosphere() {
    let card = source_card(0.001, true, 401.2);
    let fresh = bundle(&card, true, 450.0);
    assert_eq!(fresh.property(SURFACE_TENSION_PROPERTY).unwrap().value_si(), 0.002);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_capillary(&card, true, |thermal, gate| {
        let mut gas = gas(100.0);
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let bits = thermal.flow().lattice().grid().state_bits();
        let mass = thermal.flow().ledger_mass_kg().unwrap();
        let energy = thermal.total_enthalpy_j().unwrap();
        let material = thermal.material_at([0, 3, 3]).unwrap().clone();
        let wall_h = bodies[0].specific_enthalpy_j_kg();
        let wall_source = bodies[0].material().clone();
        let gas_state = gas.state().clone();
        let gas_point = gas.point().clone();
        // The independent wall reference admits T=400.961538... after step 1,
        // then reaches 401.864956... after step 2. A fresh fallback sigma exists,
        // but the retained source must refuse instead of switching mid-trajectory.
        let error = thermal.step_with_finite_gas_environment(
            &[], &[], controls(), &mut bodies, &[body::contact()], &[],
            &mut gas, &[gas_face()], gas_controls(), false, gate,
        ).unwrap_err();
        assert!(matches!(error, FreeSurfaceEnthalpyError::Material(_)));
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.flow().ledger_mass_kg().unwrap(), mass);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
        assert_eq!(thermal.material_at([0, 3, 3]).unwrap(), &material);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_h);
        assert_eq!(bodies[0].material(), &wall_source);
        assert_eq!(gas.state(), &gas_state);
        assert_eq!(gas.point(), &gas_point);

        // A smaller, valid gas-only transfer remains retryable from the exact
        // original histories. Cancellation then preserves that accepted state.
        thermal.step_with_finite_gas_environment(
            &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[gas_face()],
            gas_controls(), false, gate,
        ).unwrap();
        let accepted = thermal.flow().lattice().grid().state_bits();
        let accepted_energy = thermal.total_enthalpy_j().unwrap();
        let accepted_gas = gas.state().clone();
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(thermal.step_with_finite_gas_environment(
            &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[gas_face()],
            gas_controls(), false, &cancelled,
        ).is_err());
        assert_eq!(thermal.flow().lattice().grid().state_bits(), accepted);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), accepted_energy);
        assert_eq!(gas.state(), &accepted_gas);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_h);
    });
}
