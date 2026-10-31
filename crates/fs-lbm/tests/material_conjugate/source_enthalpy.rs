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

fn controls() -> MaterialBuoyancyConfig {
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
