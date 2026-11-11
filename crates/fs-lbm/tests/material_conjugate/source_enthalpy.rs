//! Source Cp(T), complete h(T), and actual steady native-channel heat transport.
use super::*;
use fs_evidence::ValidityDomain;
use fs_lbm::conjugate::ConvectionScheme;
use fs_lbm::material::conjugate::buoyant::enthalpy::MaterialEnthalpyConfig;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};

#[path = "source_enthalpy/transient.rs"]
mod transient;

#[test]
fn g1_source_caloric_endpoint_checks_actual_heat_without_changing_temperature() {
    use fs_lbm::conjugate::{
        AdvectiveEnthalpyField, AffineSpecificEnthalpy, ChtError, FlowFace, FlowField,
        FluidProperties,
    };
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Fluid).unwrap();
    let flow = FlowField::from_face_velocity(
        &domain,
        [
            FlowFace::Fixed,
            FlowFace::Fixed,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
        ],
        |axis, _| if axis == 0 { 1.0 } else { 0.0 },
    )
    .unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    let mut setup = ThermalSetup::new(faces);
    // h(theta)=theta+theta^2, Cp=1+2*theta. For theta=[1,3,4], the
    // independent conduction/upwind balances require sources [2,11,9] W.
    // Heat leaving is 2*theta[0] + h(theta[2]) = 22 W.
    setup.power_w = vec![2.0, 11.0, 9.0];
    let temperature = [301.0, 303.0, 304.0];
    let cells = temperature.map(|t| {
        let theta = t - 300.0;
        Some(AffineSpecificEnthalpy {
            temperature_k: t,
            enthalpy_j_kg: theta + theta * theta,
            heat_capacity_j_kg_k: 1.0 + 2.0 * theta,
        })
    });
    let field = AdvectiveEnthalpyField {
        cells: &cells,
        boundary_j_kg: [Some(0.0), Some(0.0), None, None, None, None],
    };
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 999.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    let config = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let evaluate = |values: &[f64], field, gate: &CancelGate| {
        config.evaluate_with_enthalpy(
            &domain,
            &fluid,
            &[],
            &flow,
            &setup,
            &[[1.0; 3]; 3],
            field,
            values,
            gate,
        )
    };
    let result = evaluate(&temperature, field, &CancelGate::new()).unwrap();
    assert_eq!(result.temperature, temperature);
    assert_eq!(result.report.iterations, 0);
    near(result.report.balance.source_w, 22.0, 1e-14);
    near(result.report.balance.advective_outflow_w, 20.0, 1e-14);
    near(result.report.balance.boundary_outflow_w, 22.0, 1e-14);
    assert!(result.report.balance.residual_w.abs() < 1e-13);
    // Re-anchoring the same tangent preserves every endpoint enthalpy and
    // heat prediction, but an endpoint audit must still refuse extrapolation.
    let stale = cells.map(|h| {
        h.map(|h| AffineSpecificEnthalpy {
            temperature_k: h.temperature_k - 1.0,
            enthalpy_j_kg: h.enthalpy_j_kg - h.heat_capacity_j_kg_k,
            ..h
        })
    });
    for (cell, &t) in temperature.iter().enumerate() {
        assert_eq!(
            stale[cell].unwrap().at_temperature(t),
            cells[cell].unwrap().enthalpy_j_kg,
        );
    }
    assert!(matches!(
        evaluate(
            &temperature,
            AdvectiveEnthalpyField {
                cells: &stale,
                ..field
            },
            &CancelGate::new()
        ),
        Err(ChtError::InvalidInput {
            field: "thermal.material_endpoint",
            ..
        })
    ));
    // Matching anchors do not replace the independent native heat residual:
    // h(theta)=theta, Cp=1 cannot supply the required nonlinear heat fluxes.
    let wrong_law = temperature.map(|t| {
        Some(AffineSpecificEnthalpy {
            temperature_k: t,
            enthalpy_j_kg: t - 300.0,
            heat_capacity_j_kg_k: 1.0,
        })
    });
    assert!(matches!(
        evaluate(
            &temperature,
            AdvectiveEnthalpyField {
                cells: &wrong_law,
                ..field
            },
            &CancelGate::new()
        ),
        Err(ChtError::SolverNotConverged { .. })
    ));
    assert!(evaluate(&temperature[..2], field, &CancelGate::new()).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        evaluate(&temperature, field, &gate),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(temperature, [301.0, 303.0, 304.0]);
}

#[test]
fn g0_source_caloric_endpoint_preserves_imposed_and_compact_temperatures() {
    use fs_lbm::Face3;
    use fs_lbm::conjugate::{
        AdvectiveEnthalpyField, ChtError, CompactComponent, FlowField, FluidProperties,
    };
    let domain = VoxelDomain::from_fn(1, 1, 4, 1.0, |_| Voxel::Solid(0)).unwrap();
    let flow = FlowField::quiescent(&domain);
    let solids = [SolidMaterial::new("synthetic compact surroundings", 1.0)];
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[4] = ThermalFace::Temperature(300.0);
    faces[5] = ThermalFace::Temperature(300.0);
    let mut setup = ThermalSetup::new(faces);
    setup.fixed_temperature.push((0, 300.0));
    setup.compact_components.push(CompactComponent {
        lo: [0, 0, 1],
        hi: [1, 1, 3],
        board_face: Face3::ZMin,
        power_w: 0.0,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: None,
    });
    let before = setup.clone();
    let config = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let evaluate = |temperature: &[f64]| {
        config.evaluate_with_enthalpy(
            &domain,
            &fluid,
            &solids,
            &flow,
            &setup,
            &[[1.0; 3]; 4],
            AdvectiveEnthalpyField {
                cells: &[None; 4],
                boundary_j_kg: [None; 6],
            },
            temperature,
            &CancelGate::new(),
        )
    };
    let result = evaluate(&[300.0; 4]).unwrap();
    assert_eq!(result.temperature, [300.0; 4]);
    assert_eq!(result.junctions[0].temperature_k, 300.0);
    assert_eq!(result.report.iterations, 0);
    // Both perturbations fit inside the linear residual tolerance, but an
    // endpoint audit may neither relax a prescribed value nor rewrite aliases.
    for cell in [0, 2] {
        let mut temperature = [300.0; 4];
        temperature[cell] = f64::from_bits(300.0_f64.to_bits() + 1);
        let unchanged = temperature;
        assert!(matches!(
            evaluate(&temperature),
            Err(ChtError::InvalidInput {
                field: "thermal.material_endpoint",
                ..
            })
        ));
        assert_eq!(temperature, unchanged);
    }
    assert_eq!(setup, before);
}

#[test]
fn g1_source_caloric_step_endpoint_closes_exact_storage_and_limited_fluxes() {
    use fs_lbm::conjugate::{
        AdvectiveEnthalpyField, AffineSpecificEnthalpy, ChtError, FlowFace, FlowField,
        FluidProperties,
    };
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Fluid).unwrap();
    let flow = FlowField::from_face_velocity(
        &domain,
        [
            FlowFace::Fixed,
            FlowFace::Fixed,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
        ],
        |axis, _| if axis == 0 { 1.0 } else { 0.0 },
    )
    .unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    let previous = [300.0; 3];
    let temperature = [301.0, 303.0, 304.0];
    let cells = temperature.map(|t| {
        let theta = t - 300.0;
        Some(AffineSpecificEnthalpy {
            temperature_k: t,
            enthalpy_j_kg: theta + theta * theta,
            heat_capacity_j_kg_k: 1.0 + 2.0 * theta,
        })
    });
    let field = AdvectiveEnthalpyField {
        cells: &cells,
        boundary_j_kg: [Some(0.0), Some(0.0), None, None, None, None],
    };
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 999.0,
        conductivity_w_m_k: 19.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    for scheme in [ConvectionScheme::Upwind, ConvectionScheme::EnthalpyMinmod] {
        let config = EnergyConfig { scheme, ..EnergyConfig::default() };
        let mut setup = ThermalSetup::new(faces);
        // h(theta)=theta+theta^2, unit masses, dt=1/2: the exact secant
        // capacities [2,4,5] store [2,12,20] J. The independently computed
        // upwind conduction/advection terms are [2,11,9] W. Minmod raises
        // the middle donor's outgoing h by min(10,8)/2=4 J/kg; this adds
        // +4/-4 W to its two incident cell balances, with zero net heat.
        let correction = if scheme == ConvectionScheme::EnthalpyMinmod { 4.0 } else { 0.0 };
        setup.power_w = vec![6.0, 35.0 + correction, 49.0 - correction];
        let before = setup.clone();
        let evaluate = |capacity: &[f64], values: &[f64], field, gate: &CancelGate| {
            config.evaluate_step_with_enthalpy(
                &domain, &fluid, &[], &flow, &setup, &[[1.0; 3]; 3],
                field, capacity, &previous, 0.5, values, gate,
            )
        };
        let (result, record) = evaluate(
            &[2.0, 4.0, 5.0], &temperature, field, &CancelGate::new(),
        ).unwrap();
        assert_eq!(result.temperature, temperature);
        assert_eq!(result.report.iterations, 0);
        assert_eq!(record.iterations, 0);
        near(record.stored_energy_change_j, 34.0, 1e-13);
        near(record.source_j, 45.0, 1e-13);
        near(record.boundary_outflow_j, 11.0, 1e-13);
        assert!(record.closure_j.abs() < 1e-11);
        near(result.report.balance.advective_outflow_w, 20.0, 1e-13);
        near(result.report.balance.residual_w, 68.0, 1e-13);
        // Point Cp at the new temperature is NOT a backward-Euler energy
        // secant. The wrong storage image must fail the original heat check.
        assert!(matches!(
            evaluate(&[3.0, 7.0, 9.0], &temperature, field, &CancelGate::new()),
            Err(ChtError::SolverNotConverged { .. })
        ));
        let stale = cells.map(|h| h.map(|h| AffineSpecificEnthalpy {
            temperature_k: h.temperature_k - 1.0,
            enthalpy_j_kg: h.enthalpy_j_kg - h.heat_capacity_j_kg_k,
            ..h
        }));
        assert!(matches!(
            evaluate(
                &[2.0, 4.0, 5.0], &temperature,
                AdvectiveEnthalpyField { cells: &stale, ..field }, &CancelGate::new(),
            ),
            Err(ChtError::InvalidInput { field: "thermal.material_endpoint", .. })
        ));
        assert!(evaluate(&[2.0, 4.0], &temperature, field, &CancelGate::new()).is_err());
        assert!(evaluate(&[2.0, 4.0, 5.0], &temperature[..2], field, &CancelGate::new()).is_err());
        let gate = CancelGate::new();
        gate.request();
        assert!(matches!(
            evaluate(&[2.0, 4.0, 5.0], &temperature, field, &gate),
            Err(ChtError::Cancelled)
        ));
        assert_eq!(setup, before);
        assert_eq!(previous, [300.0; 3]);
        assert_eq!(temperature, [301.0, 303.0, 304.0]);
    }
}

#[test]
fn g1_source_caloric_step_endpoint_uses_compact_energy_and_fixed_reservoir_jumps() {
    use fs_lbm::Face3;
    use fs_lbm::conjugate::{
        ChtError, CompactComponent, CompactHeatStorage, CompactHeatStorageLaw,
        FlowField, FluidProperties,
    };
    use fs_material::thermal::transport::storage::TabulatedHeatStorage;
    let domain = VoxelDomain::from_fn(1, 1, 4, 1.0, |_| Voxel::Solid(0)).unwrap();
    let flow = FlowField::quiescent(&domain);
    let solids = [SolidMaterial::new("synthetic compact surroundings", 123.0)];
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 999.0,
        conductivity_w_m_k: 19.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    let capacity = [2.0, 0.0, 0.0, 2.0];
    let previous = [300.0; 4];
    let temperature = [304.0, 310.0, 310.0, 304.0];
    let config = EnergyConfig::default();
    for imposed in [false, true] {
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        setup.power_w = vec![3.2, 0.0, 0.0, 3.2];
        if imposed {
            setup.power_w[0] = 0.0;
            setup.fixed_temperature.push((0, 304.0));
        }
        setup.compact_components.push(CompactComponent {
            lo: [0, 0, 1],
            hi: [1, 1, 3],
            board_face: Face3::ZMin,
            power_w: 79.6,
            junction_to_case_k_w: 1.0,
            junction_to_board_k_w: 1.0,
            heat_capacity_j_k: None,
        });
        setup.compact_heat_storage.push(CompactHeatStorage {
            component: 0,
            law: CompactHeatStorageLaw::Tabulated(
                TabulatedHeatStorage::try_new(&[[300.0, 2.0], [320.0, 22.0]]).unwrap(),
            ),
        });
        let before = setup.clone();
        let evaluate = |capacity: &[f64], previous: &[f64], temperature: &[f64], dt| {
            config.evaluate_step_with_conductivity(
                &domain, &fluid, &solids, &flow, &setup, &[[2.0; 3]; 4],
                capacity, previous, dt, temperature, &CancelGate::new(),
            )
        };
        let (result, record) = evaluate(&capacity, &previous, &temperature, 1.0).unwrap();
        // Each branch is G=1/(1+0.5/2)=0.8 W/K and carries 4.8 W.
        // Junction C(T)=2+(T-300) stores 20+50=70 J, not startup-Cp*10.
        // Ordinary cells store 2*4=8 J each; aliases add no second storage.
        assert_eq!(result.temperature, temperature);
        assert_eq!(result.report.iterations, 0);
        assert_eq!(record.iterations, 0);
        near(result.junctions[0].board_w, 4.8, 1e-12);
        near(result.junctions[0].case_w, 4.8, 1e-12);
        near(record.stored_energy_change_j, 86.0, 1e-12);
        near(record.source_j, if imposed { 82.8 } else { 86.0 }, 1e-12);
        // The prescribed reservoir supplies its 8 J temperature jump while
        // receiving 4.8 J from the junction: net injection is 3.2 J.
        near(record.boundary_outflow_j, if imposed { -3.2 } else { 0.0 }, 1e-12);
        assert!(record.closure_j.abs() < 1e-10);
        assert!(evaluate(&[2.0, 1.0, 0.0, 2.0], &previous, &temperature, 1.0).is_err());
        assert!(evaluate(&capacity, &previous[..3], &temperature, 1.0).is_err());
        assert!(evaluate(&capacity, &previous, &temperature, 0.0).is_err());
        let mut changed = temperature;
        let perturbed = if imposed { 0 } else { 2 };
        changed[perturbed] = f64::from_bits(changed[perturbed].to_bits() + 1);
        assert!(matches!(
            evaluate(&capacity, &previous, &changed, 1.0),
            Err(ChtError::InvalidInput { field: "thermal.material_endpoint", .. })
        ));
        assert!(evaluate(&capacity, &previous, &[304.0, 321.0, 321.0, 304.0], 1.0).is_err());
        assert_eq!(setup, before);
        assert_eq!(previous, [300.0; 4]);
        assert_eq!(temperature, [304.0, 310.0, 310.0, 304.0]);
    }
}

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), temperature)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap()
        .with("solute_mass_fraction", 0.2)
        .unwrap()
}

fn card(cp_slope: f64) -> MaterialCard {
    variable_card(cp_slope, 0.0, 300.0, 340.0)
}

fn variable_card(cp_slope: f64, density_slope: f64, lower: f64, upper: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (
            DENSITY_PROPERTY,
            Density::DIMS,
            [
                1.0 - density_slope * (lower - 300.0),
                1.0 - density_slope * (upper - 300.0),
            ],
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [
                2.0 + cp_slope * (lower - 300.0),
                2.0 + cp_slope * (upper - 300.0),
            ],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            [1.0, 1.0],
        ),
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            [1.0, 1.0],
        ),
    ] {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(lower, values[0]), (upper, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute(), lower, upper)
                    .with("pressure", 100000.0, 100000.0)
                    .with("solute_mass_fraction", 0.2, 0.2),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic caloric-flow oracle, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic coolant".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn caloric(card: &MaterialCard, datum: f64, upper: f64) -> ResolvedHeatCapacityCurve {
    caloric_span(card, datum, 300.0, upper)
}

fn caloric_span(
    card: &MaterialCard,
    datum: f64,
    lower: f64,
    upper: f64,
) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(lower),
            upper: &point(upper),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        datum,
        1e-5,
    )
    .unwrap()
}

fn liquid(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn controls() -> MaterialBuoyancyConfig<'static> {
    MaterialBuoyancyConfig {
        boussinesq: BoussinesqConfig {
            gravity_m_s2: [0.0; 3],
            max_relative_density_departure: 0.0,
            max_relative_specific_heat_departure: 0.9,
        },
        temperature_tolerance_k: 1e-9,
        ..natural_material_controls()
    }
}

fn enthalpy_controls() -> MaterialEnthalpyConfig {
    MaterialEnthalpyConfig {
        reference_temperature_k: 300.0,
        maximum_flux_defect_w: 1e-10,
    }
}

fn channel(reverse: bool) -> (VoxelDomain, SimpleConfig, ThermalSetup, EnergyConfig) {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut flow_faces = [FvBoundary::Symmetry; 6];
    flow_faces[usize::from(reverse)] = FvBoundary::Inlet {
        velocity: [if reverse { -0.1 } else { 0.1 }, 0.0, 0.0],
    };
    flow_faces[usize::from(!reverse)] = FvBoundary::Outlet;
    let mut thermal_faces = [ThermalFace::Adiabatic; 6];
    // An incoming pressure opening retains the actual reservoir h(310). No
    // artificial inlet diffusion contaminates this independent caloric balance.
    thermal_faces[0] = ThermalFace::Outflow {
        backflow_temperature: 310.0,
    };
    thermal_faces[1] = ThermalFace::Outflow {
        backflow_temperature: 310.0,
    };
    let mut thermal = ThermalSetup::new(thermal_faces);
    thermal.power_w = vec![2.6];
    (
        domain,
        SimpleConfig::new(flow_faces),
        thermal,
        EnergyConfig {
            scheme: ConvectionScheme::Upwind,
            ..EnergyConfig::default()
        },
    )
}

#[test]
fn source_enthalpy_changes_actual_channel_heat_and_preserves_datum_direction_and_sources() {
    let mut temperatures = Vec::new();
    for cp_slope in [0.0, 0.04] {
        let card = card(cp_slope);
        let law = liquid(&card);
        let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
        let mut datum_temperatures = Vec::new();
        for datum in [0.0, 1e9] {
            let caloric = caloric(&card, datum, 340.0);
            for reverse in [false, true] {
                let (domain, flow, thermal, energy) = channel(reverse);
                let solve = |_, t, _: &FluidTransportState| {
                    law.at_temperature(Temperature::new(t))
                        .map_err(|error| MaterialError::State {
                            what: error.to_string(),
                        })
                };
                // Exercise both entries; an empty solid source table is an
                // ordinary all-fluid domain on the SAME shared fixed point.
                let result = if reverse {
                    bound.solve_natural_convection_with_material_solids_and_enthalpy(
                        &domain,
                        &[],
                        &flow,
                        &thermal,
                        &energy,
                        &[310.0],
                        controls(),
                        &caloric,
                        enthalpy_controls(),
                        &CancelGate::new(),
                        solve,
                    )
                } else {
                    bound.solve_natural_convection_with_enthalpy(
                        &domain,
                        &[],
                        &flow,
                        &thermal,
                        &energy,
                        &[310.0],
                        controls(),
                        &caloric,
                        enthalpy_controls(),
                        &CancelGate::new(),
                        solve,
                    )
                }
                .unwrap();
                let convection = &result.convection;
                let t = convection.energy.temperature[0];
                // Cp=2+s*(T-300); mdot=.1, Qin=2.6 W and Tin=310 K.
                // Integrating source Cp gives Tout=323 K (s=0) or 320 K (s=.04).
                let expected = if cp_slope == 0.0 { 323.0 } else { 320.0 };
                assert!((t - expected).abs() < 1e-5, "{t} != {expected}");
                near(
                    convection
                        .flow
                        .field
                        .outward(&domain, 0, usize::from(!reverse)),
                    0.1,
                    1e-12,
                );
                near(
                    0.1 * caloric.specific_enthalpy_change(310.0, t).unwrap(),
                    2.6,
                    1e-9,
                );
                near(
                    convection.energy.report.balance.advective_outflow_w,
                    2.6,
                    1e-9,
                );
                assert!(convection.energy.report.balance.residual_w.abs() < 1e-9);
                assert!(result.enthalpy.report.maximum_face_flux_defect_w <= 1e-10);
                assert!(result.enthalpy.report.maximum_cell_flux_defect_w <= 1e-10);
                assert!(core::ptr::eq(result.caloric, &caloric));
                near(
                    result.enthalpy.boundary_j_kg[0].unwrap(),
                    20.0 + 50.0 * cp_slope,
                    1e-6,
                );
                let source = convection.endpoint[0].material.material().material();
                assert_eq!(source.card_identity(), card.content_hash());
                let receipt = &source
                    .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt;
                assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
                assert!(
                    receipt
                        .query_point
                        .contains(&("solute_mass_fraction".into(), 0.2))
                );
                datum_temperatures.push(t);
            }
        }
        for t in &datum_temperatures {
            near(*t, datum_temperatures[0], 1e-9);
        }
        temperatures.push(datum_temperatures[0]);
    }
    assert!(temperatures[0] - temperatures[1] > 2.9);
}

#[test]
fn g1_radiating_source_solid_and_caloric_channel_share_one_steady_heat_balance() {
    use fs_lbm::conjugate::{RadiationConfig, STEFAN_BOLTZMANN};

    // At the independent target Ts=325 K, Tf=320 K, the original source
    // solid has k=1.5 and the fluid k=1. Their half cells conduct 6 W.
    // Cp=2+0.04*(T-300), so h(320)-h(310)=26 J/kg; mdot=6/26
    // carries exactly those 6 W out. The remaining solid power radiates.
    let fluid_card = card(0.04);
    let fluid = liquid(&fluid_card);
    let caloric = caloric(&fluid_card, 1e9, 340.0);
    let solid_card = super::solid_feedback::card([0.02; 3], 340.0, false);
    let source_solid = super::solid_feedback::law(&solid_card);
    let bound = MaterialChtFluid::try_new(fluid.reference_state()).unwrap();
    let mut domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    let mut boundaries = [FvBoundary::Symmetry; 6];
    boundaries[0] = FvBoundary::Inlet {
        velocity: [6.0 / 26.0, 0.0, 0.0],
    };
    boundaries[1] = FvBoundary::Outlet;
    let mut flow = SimpleConfig::new(boundaries);
    flow.tolerance = 1e-11;
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Outflow {
        backflow_temperature: 310.0,
    };
    faces[1] = faces[0];
    let mut thermal = ThermalSetup::new(faces);
    let radiated = 0.8 * STEFAN_BOLTZMANN * (325.0_f64.powi(4) - 300.0_f64.powi(4));
    thermal.power_w = vec![6.0 + radiated, 0.0];
    let radiation = RadiationConfig::new(vec![0.8], [Some(300.0); 6]);
    let controls = MaterialBuoyancyConfig {
        radiation: Some(&radiation),
        ..controls()
    };
    let energy = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let gate = CancelGate::new();
    for source_defined_solid in [false, true] {
        let resolve = |_, t, _: &FluidTransportState| {
            fluid
                .at_temperature(Temperature::new(t))
                .map_err(|error| MaterialError::State {
                    what: error.to_string(),
                })
        };
        let run = if source_defined_solid {
            bound.solve_natural_convection_with_material_solids_and_enthalpy(
                &domain,
                std::slice::from_ref(&source_solid),
                &flow,
                &thermal,
                &energy,
                &[320.0; 2],
                controls,
                &caloric,
                enthalpy_controls(),
                &gate,
                resolve,
            )
        } else {
            bound.solve_natural_convection_with_enthalpy(
                &domain,
                &[SolidMaterial::new("reference plate", 1.5)],
                &flow,
                &thermal,
                &energy,
                &[320.0; 2],
                controls,
                &caloric,
                enthalpy_controls(),
                &gate,
                resolve,
            )
        }
        .unwrap();
        let result = &run.convection;
        near(result.energy.temperature[0], 325.0, 1e-9);
        near(result.energy.temperature[1], 320.0, 1e-9);
        near(result.energy.solid_to_fluid_heat_w(&domain), 6.0, 1e-8);
        near(result.energy.report.balance.advective_outflow_w, 6.0, 1e-8);
        near(result.report.radiation.unwrap().radiated_w, radiated, 1e-9);
        assert!(
            result
                .report
                .radiation
                .unwrap()
                .physical_energy_residual_w
                .abs()
                < 1e-8
        );
        assert!(
            run.enthalpy.report.maximum_cell_flux_defect_w
                <= enthalpy_controls().maximum_flux_defect_w
        );
        assert!(core::ptr::eq(run.caloric, &caloric));
        if source_defined_solid {
            assert_eq!(
                result.solid_endpoint[0].material.card_identity(),
                solid_card.content_hash()
            );
            near(result.solid_endpoint[0].conductivity_w_m_k[1], 1.5, 1e-9);
        }
    }
}

#[test]
fn source_enthalpy_requires_local_heat_closure_even_when_other_material_defects_pass() {
    let card = card(0.04);
    let law = liquid(&card);
    let caloric = caloric(&card, 0.0, 340.0);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    let (domain, flow, thermal, energy) = channel(false);
    let mut limits = controls();
    limits.maximum_passes = 1;
    limits.temperature_tolerance_k = 100.0;
    let result = bound.solve_natural_convection_with_enthalpy(
        &domain,
        &[],
        &flow,
        &thermal,
        &energy,
        &[310.0],
        limits,
        &caloric,
        enthalpy_controls(),
        &CancelGate::new(),
        |_, t, _| {
            law.at_temperature(Temperature::new(t))
                .map_err(|error| MaterialError::State {
                    what: error.to_string(),
                })
        },
    );
    match result {
        Err(MaterialBuoyancyError::EnthalpyNotConverged { material, enthalpy }) => {
            assert!(material.temperature_defect_k < 100.0);
            assert_eq!(material.transport_relative_defect, 0.0);
            assert_eq!(material.force_defect_n_m3, 0.0);
            assert!(enthalpy.maximum_face_flux_defect_w > 0.1);
            assert!(enthalpy.maximum_cell_flux_defect_w > 0.1);
        }
        _ => panic!("a small material defect cannot accept an unresolved h(T) flux"),
    }
}

#[test]
fn source_enthalpy_endpoint_domain_cancellation_and_original_cp_bound_refuse_atomically() {
    let card = card(0.04);
    let law = liquid(&card);
    let short_caloric = caloric(&card, 0.0, 315.0);
    let full_caloric = caloric(&card, 0.0, 340.0);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    let (domain, flow, thermal, energy) = channel(false);
    let saved = thermal.clone();
    let resolve = |_, t, _: &FluidTransportState| {
        law.at_temperature(Temperature::new(t))
            .map_err(|error| MaterialError::State {
                what: error.to_string(),
            })
    };
    assert!(matches!(
        bound.solve_natural_convection_with_enthalpy(
            &domain,
            &[],
            &flow,
            &thermal,
            &energy,
            &[310.0],
            controls(),
            &short_caloric,
            enthalpy_controls(),
            &CancelGate::new(),
            resolve
        ),
        Err(MaterialBuoyancyError::Caloric { .. })
    ));
    let mut strict = controls();
    strict.boussinesq.max_relative_specific_heat_departure = 0.0;
    assert!(
        bound
            .solve_natural_convection_with_enthalpy(
                &domain,
                &[],
                &flow,
                &thermal,
                &energy,
                &[310.0],
                strict,
                &full_caloric,
                enthalpy_controls(),
                &CancelGate::new(),
                resolve
            )
            .is_err()
    );
    let gate = CancelGate::new();
    gate.request();
    let calls = std::cell::Cell::new(0usize);
    assert!(matches!(
        bound.solve_natural_convection_with_enthalpy(
            &domain,
            &[],
            &flow,
            &thermal,
            &energy,
            &[310.0],
            controls(),
            &full_caloric,
            enthalpy_controls(),
            &gate,
            |cell, t, reference| {
                calls.set(calls.get() + 1);
                resolve(cell, t, reference)
            }
        ),
        Err(MaterialBuoyancyError::Solver(
            fs_lbm::conjugate::ChtError::Cancelled
        ))
    ));
    assert_eq!(calls.get(), 0);
    assert_eq!(thermal, saved);
}

#[test]
fn source_enthalpy_admits_transport_at_inlet_and_potential_backflow_temperatures() {
    // rho(305)=.975 is inside the 5% bound, but rho(280)=1.10 is outside.
    // The caloric chart supports BOTH temperatures: h/Cp admission alone
    // would otherwise allow this inadmissible source density into the flow.
    let card = variable_card(0.0, 0.005, 280.0, 340.0);
    let law = liquid(&card);
    let caloric = caloric_span(&card, 0.0, 280.0, 340.0);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    let mut limits = controls();
    limits.boussinesq.max_relative_density_departure = 0.05;
    limits.boussinesq.max_relative_specific_heat_departure = 0.0;
    for cold_face in [0, 1] {
        let (domain, flow, mut thermal, energy) = channel(false);
        thermal.faces[0] = ThermalFace::Outflow {
            backflow_temperature: 305.0,
        };
        thermal.faces[1] = ThermalFace::Outflow {
            backflow_temperature: 305.0,
        };
        thermal.faces[cold_face] = ThermalFace::Outflow {
            backflow_temperature: 280.0,
        };
        // With the unsafe admission, .1*Cp*(305-280)=5 W would exactly
        // sustain the supported cell temperature; the downstream state is
        // not the reason for refusing the cold inlet. Face 1 tests declared
        // backflow support even while the current flow points outward.
        thermal.power_w[0] = if cold_face == 0 { 5.0 } else { 0.0 };
        let saved = thermal.clone();
        let queries = std::cell::RefCell::new(Vec::new());
        let result = bound.solve_natural_convection_with_enthalpy(
            &domain,
            &[],
            &flow,
            &thermal,
            &energy,
            &[305.0],
            limits,
            &caloric,
            enthalpy_controls(),
            &CancelGate::new(),
            |cell, temperature, _| {
                queries.borrow_mut().push((cell, temperature));
                law.at_temperature(Temperature::new(temperature))
                    .map_err(|error| MaterialError::State {
                        what: error.to_string(),
                    })
            },
        );
        match result {
            Err(MaterialBuoyancyError::Material {
                cell: 0,
                temperature_k: 280.0,
                source,
            }) => {
                assert!(source.to_string().contains("density relative departure"));
            }
            _ => panic!("the original source-density bound must also admit every inlet/backflow"),
        }
        assert!(queries.borrow().contains(&(0, 305.0)));
        assert!(queries.borrow().contains(&(0, 280.0)));
        assert_eq!(thermal, saved);
    }
}

#[test]
fn source_enthalpy_graded_heated_buoyant_channel_closes_each_cell_and_boundary_heat() {
    let card = variable_card(0.04, 0.001, 300.0, 340.0);
    let law = liquid(&card);
    let caloric = caloric(&card, 1e9, 340.0);
    let bound = MaterialChtFluid::try_new(law.reference_state()).unwrap();
    let domain =
        VoxelDomain::graded([vec![0.4, 0.8, 1.2, 0.6], vec![0.5, 1.0, 0.75], vec![0.4]]).unwrap();
    let mut flow_faces = [FvBoundary::Symmetry; 6];
    flow_faces[0] = FvBoundary::Inlet {
        velocity: [0.1, 0.0, 0.0],
    };
    flow_faces[1] = FvBoundary::Outlet;
    let mut flow = SimpleConfig::new(flow_faces);
    flow.tolerance = 1e-9;
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Outflow {
        backflow_temperature: 310.0,
    };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 310.0,
    };
    let mut thermal = ThermalSetup::new(faces);
    thermal.power_w = vec![0.0; domain.cell_count()];
    thermal.power_w[domain.index(1, 0, 0)] = 1.5;
    let mut limits = controls();
    limits.boussinesq.gravity_m_s2 = [0.0, -0.2, 0.0];
    limits.boussinesq.max_relative_density_departure = 0.05;
    let energy = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };
    let result = bound
        .solve_natural_convection_with_enthalpy(
            &domain,
            &[],
            &flow,
            &thermal,
            &energy,
            &vec![310.0; domain.cell_count()],
            limits,
            &caloric,
            enthalpy_controls(),
            &CancelGate::new(),
            |_, temperature, _| {
                law.at_temperature(Temperature::new(temperature))
                    .map_err(|error| MaterialError::State {
                        what: error.to_string(),
                    })
            },
        )
        .unwrap();
    let solved = &result.convection;
    let temperatures = &solved.energy.temperature;
    assert!(
        solved
            .applied
            .iter()
            .any(|sample| sample.force_n_m3[1] > 0.001)
    );
    assert!(solved.flow.velocity_m_s.iter().any(|v| v[1].abs() > 1e-8));
    let low = temperatures.iter().copied().fold(f64::INFINITY, f64::min);
    let high = temperatures
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(high - low > 1e-3);
    let density = bound.transport().density_kg_m3();
    let h: Vec<_> = temperatures
        .iter()
        .map(|&t| caloric.specific_enthalpy_change(300.0, t).unwrap())
        .collect();
    let inlet_h = caloric.specific_enthalpy_change(300.0, 310.0).unwrap();
    let mut boundary_heat = 0.0;
    for cell in 0..domain.cell_count() {
        let mut local_heat = 0.0;
        let mut magnitude = 0.0;
        for face in 0..6 {
            let q = solved.flow.field.outward(&domain, cell, face);
            magnitude += q.abs();
            if let Some(next) = domain.neighbor(cell, face) {
                let axis = face / 2;
                // Independent finite-volume heat balance: source k=1 in
                // both cells, using their actual graded half-cell distances.
                let conductance = domain.face_area(cell, axis)
                    / (0.5 * domain.widths(cell)[axis] + 0.5 * domain.widths(next)[axis]);
                let donor_h = if q >= 0.0 { h[cell] } else { h[next] };
                local_heat +=
                    conductance * (temperatures[cell] - temperatures[next]) + density * q * donor_h;
            } else {
                let donor_h = if q >= 0.0 { h[cell] } else { inlet_h };
                let advected = density * q * donor_h;
                local_heat += advected;
                boundary_heat += advected;
            }
        }
        assert!(
            solved.flow.field.net_outflow(&domain, cell).abs() <= 128.0 * f64::EPSILON * magnitude
        );
        assert!(
            (local_heat - thermal.power_w[cell]).abs() < 1e-8,
            "cell {cell}: physical heat residual {}",
            local_heat - thermal.power_w[cell]
        );
    }
    assert!((boundary_heat - 1.5).abs() < 1e-8);
    assert!((solved.energy.report.balance.advective_outflow_w - boundary_heat).abs() < 1e-8);
    assert!(result.enthalpy.report.maximum_face_flux_defect_w <= 1e-10);
    assert!(result.enthalpy.report.maximum_cell_flux_defect_w <= 1e-10);
}
