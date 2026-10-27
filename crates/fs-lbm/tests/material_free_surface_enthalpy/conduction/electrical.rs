//! Independent DC/thermal references on the actual moving-liquid timeline.

use super::*;
use fs_lbm::conjugate::{ElectricSetup, Electrode, ElectrodeDrive};
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceElectricalInput;
use fs_material::conductor::ElectricalPropertyBasis;

const CORE: [usize; 3] = [1, 3, 3];

fn electrical_card(
    basis: ElectricalPropertyBasis,
    scale: f64,
    varying: bool,
    shear: bool,
    upper_k: f64,
) -> MaterialCard {
    let base = source_card(0.001, shear, 500.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let reference = match basis {
        ElectricalPropertyBasis::Resistivity => 0.001 / scale,
        ElectricalPropertyBasis::Conductivity => 1000.0 * scale,
    };
    let mut claim = PropertyClaim {
        key: PropertyKey::new(basis.property_name(), basis.dimensions()),
        value: if varying {
            PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, reference * 0.75), (500.0, reference * 1.25)],
                dims: basis.dimensions(),
            }
        } else {
            PropertyValue::Scalar { value: reference, dims: basis.dimensions() }
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, upper_k)
            .with("pressure", 100000.0, 100000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: if varying { InterpolationPolicy::LinearInside }
            else { InterpolationPolicy::ConstantWithinValidity },
        observations: vec![],
        provenance: Provenance {
            source: "synthetic electrical/thermal liquid coupling fixture".into(),
            license: "test-only".into(), artifact: None,
        },
    };
    claims.insert_claim(claim.clone()).unwrap();
    if upper_k < 500.0 {
        // An available fresh source cannot replace the original selected claim.
        claim.validity = ValidityDomain::unconstrained()
            .with("T", upper_k + 0.000001, 500.0)
            .with("pressure", 100000.0, 100000.0);
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn electrical_bundle(
    card: &MaterialCard,
    basis: ElectricalPropertyBasis,
    shear: bool,
    t: f64,
) -> ResolvedMaterialStatePoint {
    let source = bundle(card, shear, t);
    let mut requirements = source.properties().iter()
        .map(|p| p.requirement().clone()).collect::<Vec<_>>();
    requirements.push(basis.requirement().unwrap());
    resolve_material_state_point(
        card, &source.query(), &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

fn current(voltage: f64) -> FreeSurfaceElectricalInput {
    FreeSurfaceElectricalInput {
        origin: [0, 3, 3], dimensions: [3, 1, 1],
        setup: ElectricSetup {
            conductors: vec![],
            electrodes: vec![
                Electrode { lo: [0; 3], hi: [1; 3], drive: ElectrodeDrive::Voltage(0.0) },
                Electrode { lo: [2, 0, 0], hi: [3, 1, 1], drive: ElectrodeDrive::Voltage(voltage) },
            ],
        },
        contacts: vec![], maximum_cells: 3, maximum_solves: 32,
        power_tolerance_w: 1e-9,
    }
}

fn with_unbound<R>(
    card: &MaterialCard,
    basis: ElectricalPropertyBasis,
    shear: bool,
    run: impl FnOnce(MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let material = electrical_bundle(card, basis, shear, 400.0);
    let source = if shear {
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity(
            material, LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected, "T", "rate",
        )
    } else {
        ResolvedLiquidSurfaceState::try_from_material(
            material, LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
        )
    }.unwrap();
    let law = LiquidTemperatureLaw::try_from_material(
        card, source.material().clone(), "T", LiquidViscosityBasis::Dynamic,
    ).unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut initial = cells();
    if !shear {
        for cell in &mut initial { cell.velocity_m_s = [0.0; 3]; }
    }
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        let mut liquid = if shear {
            scaling.free_surface3_from_cells_with_rheology(
                domain(), &initial, &law, &chart, FreeSurfaceEnthalpyControls::default(),
                RheologySolve {
                    maximum_iterations: 128, relaxation: 0.5,
                    absolute_relaxation_tolerance: 1e-13, relative_relaxation_tolerance: 1e-13,
                }, |_| 200000.0, cx, &gate,
            )
        } else {
            scaling.free_surface3_from_cells(domain(), &initial, cx).unwrap()
                .with_enthalpy(&law, &chart, FreeSurfaceEnthalpyControls::default(),
                    |_| 200000.0, &gate)
        }.unwrap();
        if shear {
            // Establish actual non-equilibrium stress before changing electrical
            // data. Equilibrium populations cannot reveal a viscosity change.
            liquid.step(&gate).unwrap();
        }
        run(liquid, &gate)
    })
}

#[test]
fn g1_current_and_deposited_heat_follow_independent_temperature_dependent_resistance() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = electrical_card(basis, 1.0, true, false, 500.0);
    with_unbound(&card, basis, false, |liquid, gate| {
        let mut liquid = liquid.with_electrical_core(current(4.0), gate).unwrap();
        let before = liquid.total_enthalpy_j().unwrap();
        let report = liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
        let electric = report.electrical.as_ref().unwrap();
        // One bulk cell, two zero-resistance electrodes: R=rho*dx/A=T/400 ohm.
        // Its real storage is rho_mass*dx^3*cp = 0.002 J/K. Compute the explicit
        // native thermal recurrence independently, without calling either solver.
        let (mut t, mut q, mut last_t) = (400.0, 0.0, 400.0);
        for _ in 0..4 {
            last_t = t;
            let heat = 0.00025 * 16.0 / (t / 400.0);
            t += heat / 0.002;
            q += heat;
        }
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(electric.substeps, 4);
        assert_eq!(electric.current_solves, 8);
        near(electric.last_applied.temperature_k[1], last_t, 2e-12);
        near(electric.thermal_endpoint.temperature_k[1], t, 2e-12);
        near(electric.thermal_endpoint.solution.electrodes[1].current_a, 4.0 / (t / 400.0), 2e-12);
        near(electric.joule_heat_j, q, 2e-12);
        near(electric.electrode_energy_j, q, 2e-12);
        near(report.flow.external_heat_j, q, 2e-12);
        near(liquid.total_enthalpy_j().unwrap() - before, q, 2e-9);
        near(temperature(liquid.material_at(CORE).unwrap()), t, 2e-11);
        assert!(electric.energy_defect_j.abs() < 1e-12);
        assert!(report.flow.balance_defect_j.abs() < 1e-9);
        assert_eq!(report.flow.viscous_heat_j, 0.0);
        let retained = &electric.thermal_endpoint.materials[1];
        assert_eq!(retained.material().properties().len(), 6);
        assert_eq!(retained.source_basis(), basis);
        assert_eq!(retained.source_property().answer().receipt.selected,
            electrical_bundle(&card, basis, false, 400.0)
                .property(basis.property_name()).unwrap().answer().receipt.selected);
    });
}

#[test]
fn g3_electrical_material_data_changes_native_motion_and_reuses_the_accepted_state() {
    let run = |basis, scale| {
        let card = electrical_card(basis, scale, false, true, 500.0);
        with_unbound(&card, basis, true, |liquid, gate| {
            let mut liquid = liquid.with_electrical_core(current(4.0), gate).unwrap();
            let first = liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
            let first_t = temperature(liquid.material_at(CORE).unwrap());
            let first_h = liquid.specific_enthalpy_at(CORE).unwrap();
            let first_bits = liquid.flow().lattice().grid().state_bits();
            let second = liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
            let electric = second.electrical.as_ref().unwrap();
            assert!(electric.last_applied.specific_enthalpy_j_kg[1] > first_h);
            assert!(electric.thermal_endpoint.temperature_k[1] > first_t);
            near(first.electrical.as_ref().unwrap().joule_heat_j, 0.016 * scale, 2e-12);
            near(electric.joule_heat_j, 0.016 * scale, 2e-12);
            let state = liquid.material_at(CORE).unwrap();
            let rate = liquid.rheology_at(CORE).unwrap().shear_rate_per_s;
            near(state.dynamic_viscosity_pa_s(),
                0.1 - 0.0002 * (temperature(state) - 400.0) + 0.002 * rate, 2e-12);
            assert_eq!(electric.thermal_endpoint.materials[1].source_basis(), basis);
            (first_bits, first_t, liquid.flow().lattice().grid().state_bits())
        })
    };
    let resistive = run(ElectricalPropertyBasis::Resistivity, 1.0);
    let conductive = run(ElectricalPropertyBasis::Conductivity, 1.0);
    assert_eq!(resistive, conductive, "the equivalent electrical source basis is a data choice");
    let stronger = run(ElectricalPropertyBasis::Conductivity, 4.0);
    assert!(stronger.1 > conductive.1);
    let delta = stronger.0.iter().zip(&conductive.0)
        .map(|(&a, &b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
        .fold(0.0_f64, f64::max);
    assert!(delta > 1e-12, "electrical material data must reach actual native collision populations");
    assert_ne!(stronger.2, conductive.2);
}

#[test]
fn g1_joule_heat_wall_and_finite_atmosphere_share_one_original_thermal_trajectory() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = electrical_card(basis, 1.0, true, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_unbound(&card, basis, false, |liquid, gate| {
        let mut liquid = liquid.with_electrical_core(current(4.0), gate).unwrap();
        let mut atmosphere = gas(100.0);
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let faces = (0..3).map(|x| FreeSurfaceConductionFace {
            cells: [[x, 3, 3], [x + 1, 3, 3]], area_m2: 1e-6, distances_m: [0.0005; 2],
        }).collect::<Vec<_>>();
        let before = liquid.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j()
            + atmosphere.mass_kg() * atmosphere.state().specific_internal_energy_j_kg();
        let report = liquid.step_with_finite_gas_environment(
            &faces, &[], controls(), &mut bodies, &[body::contact()], &[],
            &mut atmosphere, &[gas_face()], gas_controls(), false, gate,
        ).unwrap().conduction;
        // Six physical stores: four liquid cells, one wall cell and NASA9 gas.
        // Resolve every flux from the common old state before applying any heat.
        let (mut t, mut wall_t, mut gas_t) = ([400.0; 4], 450.0, 450.0);
        let (mut qe, mut qw, mut qg) = (0.0, 0.0, 0.0);
        let r = 8.314_462_618_153_24 / 0.028;
        let mg = 100000.0 * 1e-8 / (r * 450.0);
        for _ in 0..4 {
            let mut q = [0.0; 4];
            let joule = 0.00025 * 16.0 / (t[1] / 400.0);
            q[1] += joule;
            qe += joule;
            for i in 0..3 {
                let g = 1e-6 / (0.0005 / (t[i] - 200.0) + 0.0005 / (t[i + 1] - 200.0));
                let heat = 0.00025 * g * (t[i + 1] - t[i]);
                q[i] += heat;
                q[i + 1] -= heat;
            }
            let g = 1e-6 / (0.0005 / (t[0] - 200.0) + 0.0005 / (0.5 * (wall_t - 200.0)));
            let wall_heat = 0.00025 * g * (wall_t - t[0]);
            let pressure = 100000.0 * gas_t / 450.0;
            let kg = 2.6 * (0.5 + (pressure - 1000.0) / 249000.0);
            let g = 1e-6 / (0.0005 / (t[3] - 200.0) + 0.0005 / kg);
            let gas_heat = 0.00025 * g * (gas_t - t[3]);
            q[0] += wall_heat;
            q[3] += gas_heat;
            for i in 0..4 { t[i] += q[i] / 0.002; }
            wall_t -= wall_heat / 0.001;
            gas_t = (((2.0_f64 + 0.001 * gas_t).powi(2) - 0.002 * gas_heat / (mg * r)).sqrt()
                - 2.0) / 0.001;
            qw += wall_heat;
            qg += gas_heat;
        }
        assert_eq!(report.thermal.substeps, 4);
        let electric = report.electrical.as_ref().unwrap();
        near(electric.joule_heat_j, qe, 2e-11);
        near(electric.thermal_endpoint.temperature_k[1], t[1], 2e-11);
        near(report.bodies.as_ref().unwrap().heat_to_liquid_j, qw, 2e-11);
        near(report.environment.as_ref().unwrap().heat_to_liquid_j, qg, 2e-11);
        near(temperature(&report.thermal_endpoint[2].materials[1]), t[3], 2e-11);
        near(atmosphere.state().state().temperature, gas_t, 2e-11);
        near(bodies[0].specific_enthalpy_j_kg(), (wall_t - 300.0) * 500.0, 2e-11);
        near(report.flow.external_heat_j, qe + qg, 2e-11);
        let after = liquid.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j()
            + atmosphere.mass_kg() * atmosphere.state().specific_internal_energy_j_kg();
        assert!((after - before - electric.electrode_energy_j).abs() < 1e-9);
        assert!(report.finite_gas.as_ref().unwrap().combined_thermal_balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g4_late_electrical_source_and_solve_budget_failures_keep_all_participants_atomic() {
    for (upper_k, maximum_solves) in [(403.0, 32), (500.0, 3)] {
        let basis = ElectricalPropertyBasis::Resistivity;
        let card = electrical_card(basis, 1.0, false, false, upper_k);
        let solid = body::wall_card(0.5, 500.0);
        let chart = body::wall_chart(&solid);
        let context = point(450.0);
        with_unbound(&card, basis, false, |liquid, gate| {
            let mut input = current(4.0);
            input.maximum_solves = maximum_solves;
            let mut liquid = liquid.with_electrical_core(input, gate).unwrap();
            let mut atmosphere = gas(100.0);
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
            let bits = liquid.flow().lattice().grid().state_bits();
            let h = liquid.total_enthalpy_j().unwrap();
            let source = liquid.material_at(CORE).unwrap().clone();
            let wall_h = bodies[0].specific_enthalpy_j_kg();
            let wall_source = bodies[0].material().clone();
            let gas_source = atmosphere.state().clone();
            let gas_point = atmosphere.point().clone();
            // 16 W gives 2 K per substep: the source admits the FIRST 402 K
            // endpoint and refuses the second 404 K endpoint. The solve-budget
            // case also passes a full substep before refusing the next endpoint.
            for _ in 0..2 {
                let error = liquid.step_with_finite_gas_environment(
                    &[], &[], controls(), &mut bodies, &[body::contact()], &[],
                    &mut atmosphere, &[gas_face()], gas_controls(), false, gate,
                ).unwrap_err();
                if upper_k < 500.0 {
                    assert!(matches!(error, FreeSurfaceEnthalpyError::Material(_)));
                } else {
                    assert!(matches!(error, FreeSurfaceEnthalpyError::Invalid("electrical current-solve budget exhausted")));
                }
                assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
                assert_eq!(liquid.total_enthalpy_j().unwrap(), h);
                assert_eq!(liquid.material_at(CORE).unwrap(), &source);
                assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_h);
                assert_eq!(bodies[0].material(), &wall_source);
                assert_eq!(atmosphere.state(), &gas_source);
                assert_eq!(atmosphere.point(), &gas_point);
            }
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(liquid.step_with_conduction(&[], &[], controls(), &cancelled).is_err());
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), h);
        });
    }
}

#[test]
fn g0_partial_current_geometry_and_drive_skipping_are_explicit_refusals() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = electrical_card(basis, 1.0, false, false, 500.0);
    with_unbound(&card, basis, false, |liquid, gate| {
        let mut input = current(4.0);
        input.origin = [1, 3, 3]; // Includes the native x=3 interface cell.
        assert!(liquid.with_electrical_core(input, gate).is_err());
    });
    with_unbound(&card, basis, false, |liquid, gate| {
        let mut liquid = liquid.with_electrical_core(current(4.0), gate).unwrap();
        let bits = liquid.flow().lattice().grid().state_bits();
        let energy = liquid.total_enthalpy_j().unwrap();
        assert!(liquid.step(gate).is_err());
        assert!(liquid.step_with_viscous_heating(&[], gate).is_err());
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), energy);
        liquid.step_with_conduction(&[], &[], controls(), gate).unwrap();
    });
}
