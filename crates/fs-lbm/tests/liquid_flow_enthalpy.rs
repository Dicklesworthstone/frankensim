#![cfg(feature = "enthalpy-coupling")]
//! Synthetic source cards exercise native momentum/enthalpy unification.
//! These are numerical coupling fixtures, not measured liquid data.

use std::collections::BTreeMap;

use fs_conduction::transient::deposition::*;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    FacePatch, FanCurve, FvBoundary, InternalFan, SimpleConfig, VoxelDomain,
};
use fs_lbm::conjugate::enthalpy::advection::liquid::{
    LiquidEnthalpyError, LiquidEnthalpyObserverError, LiquidFlowEnthalpyInput,
};
use fs_lbm::material::conjugate::MaterialChtFluid;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const MAT: MaterialId = MaterialId(3);

#[path = "liquid_flow_enthalpy/feedback.rs"]
mod feedback;

#[path = "liquid_flow_enthalpy/mushy.rs"]
mod mushy;

#[path = "liquid_flow_enthalpy/viscous.rs"]
mod viscous;

#[path = "liquid_flow_enthalpy/mhd.rs"]
mod mhd;

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ))
}

fn point(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", absolute(), t).unwrap()
        .with("pressure", pressure).unwrap().with("lot", 7.0).unwrap()
}

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity, ValueForm::Static,
        ))
    } else { QuantitySpec::dimensional(dims) };
    PropertyClaim {
        key: PropertyKey::with_quantity(name, quantity),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", absolute(), 300.0, 810.0)
            .with("pressure", 100000.0, 200000.0).with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic liquid flow fixture, not measurements".into(),
            license: "test-only".into(), artifact: None,
        },
    }
}

fn claims(rho: f64, mu: f64, split_viscosity: bool) -> ClaimSet {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, rho),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 2.0),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 0.01),
    ] {
        claims.insert_claim(scalar(name, dims, value)).unwrap();
    }
    let mut viscosity = scalar(DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu);
    if split_viscosity {
        viscosity.validity = viscosity.validity.with("T", 300.0, 310.0);
        claims.insert_claim(viscosity.clone()).unwrap();
        viscosity.validity = viscosity.validity.with("T", 311.0, 810.0);
        viscosity.value = PropertyValue::Scalar { value: 4.0 * mu, dims: DYNAMIC_VISCOSITY_DIMS };
    }
    claims.insert_claim(viscosity).unwrap();
    claims
}

fn card(claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic liquid".into(), phase: "liquid".into(),
        process: "numerical fixture".into(), revision: 0,
    }, claims, vec![]).unwrap()
}

fn liquid(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(card, &point(300.0, 100000.0), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap()
}

fn source_caloric(
    card: &MaterialCard, pressure: f64, cp_name: &str, rho_name: &str, phase: SolidLiquidPhase,
) -> ResolvedHeatCapacityCurve {
    let cp = PropertyKey::new(cp_name, SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new(rho_name, Density::DIMS);
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T",
        lower: &point(300.0, pressure), upper: &point(810.0, pressure),
        selection: SelectionPolicy::SingleClaimOnly,
    }, phase, 0.0, 1e-8).unwrap()
}

fn liquid_caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    source_caloric(card, 100000.0, SPECIFIC_HEAT_CAPACITY_PROPERTY, DENSITY_PROPERTY,
        SolidLiquidPhase::Liquid)
}

fn config() -> StepConfig { StepConfig { max_dt_s: 0.125, energy_tolerance_j: 1e-9 } }

fn initial<'a>(
    domain: &VoxelDomain, rho: f64, caloric: &'a ResolvedHeatCapacityCurve,
) -> DepositionSimulation<'a> {
    let events = (0..domain.cell_count()).map(|c| {
        let [x, y, z] = domain.coords(c);
        let centre = domain.center(x, y, z);
        let widths = domain.widths(c);
        let mut start = centre;
        let mut end = centre;
        start[0] -= 0.5 * widths[0];
        end[0] += 0.5 * widths[0];
        TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(c as u64), segment_id: c as u64, body_id: 1,
            material_id: MAT, mass_kg: rho * domain.volume(c), specific_h_j_kg: 0.0,
            geometry: CellGeometry { start_m: start, end_m: end,
                volume_m3: domain.volume(c), surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0,
                reservoir: None, source_w: 0.0 },
        }) }
    }).collect();
    DepositionSimulation::new(events, &[(MAT, caloric.curve() as &dyn ThermalLaw)], DepositionBudget {
        max_cells: domain.cell_count(), max_links: 64,
        max_events: domain.cell_count(), max_substeps: 10000,
    }).unwrap().advance_to(0.0, config(), || Ok(())).unwrap().0
}

fn input<'a>(
    ids: &'a [Option<CellId>], liquid: &'a LiquidTemperatureLaw<'a>,
    caloric: &'a ResolvedHeatCapacityCurve,
) -> LiquidFlowEnthalpyInput<'a> {
    LiquidFlowEnthalpyInput { cell_ids: ids, material_id: MAT, liquid, caloric,
        inlet_specific_h_j_kg: [Some(0.0); 6], max_cells: ids.len(), max_connections: 64 }
}

fn throughflow() -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet { velocity: [0.25, 0.0, 0.0] };
    faces[1] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    // Same native positive fixture as flow_enthalpy.rs. The shared graph still
    // admits LOCAL mass balance to 64eps, independently of these solve criteria.
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config
}

#[derive(Clone, Default, Debug, PartialEq)]
struct Heater {
    durations: Vec<f64>,
    supplied_j: f64,
    reject_step: Option<usize>,
}
impl DepositionObserver for Heater {
    type Error = DepositionError;
    fn sources(&mut self, _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error>
    {
        Ok(BTreeMap::from([(CellId(1), 100.0)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, dt: f64, _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>, _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error>
    {
        self.durations.push(dt);
        self.supplied_j += dt * 100.0;
        if self.reject_step == Some(self.durations.len()) {
            return Err(DepositionError::InvalidInput("independent history refused endpoint"));
        }
        Ok(())
    }
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-8 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn sourced_native_duct_matches_upwind_energy_and_retains_endpoint_receipts() {
    let card = card(claims(2.0, 0.1, false));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let old = initial(&domain, 2.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let observer = solved.advect_liquid_enthalpy(
        &old, input(&ids, &liquid, &caloric), Heater::default(), &gate,
    ).unwrap();
    let (next, accepted, audit) = old.advance_observed(&observer, 1.0, config(), || Ok(())).unwrap();
    assert!(core::ptr::eq(accepted.material(), liquid.reference_state()));
    assert!(core::ptr::eq(accepted.field(), &solved.flow().field));
    assert_eq!(accepted.transport(), solved.transport());
    for flow in accepted.advection().flows() { near(flow.kg_s, 0.5); }
    let source_history = accepted.advection().inner();
    assert!(core::ptr::eq(source_history.material_law(), &liquid));
    assert!(core::ptr::eq(source_history.caloric_source(), &caloric));
    // Independent donor-cell recurrence at the actual accepted native substeps.
    let mut h = [0.0; 4];
    let mut exported_j = 0.0;
    for &dt in &source_history.inner().durations {
        exported_j += dt * 0.5 * h[3];
        let old_h = h;
        for c in 0..4 {
            let upstream = if c == 0 { 0.0 } else { old_h[c - 1] };
            let heat = if c == 1 { 100.0 } else { 0.0 };
            h[c] += dt * (heat + 0.5 * (upstream - old_h[c])) / 2.0;
        }
    }
    for (c, expected) in h.into_iter().enumerate() {
        near(next.cells()[&CellId(c as u64)].specific_h_j_kg, expected);
        near(next.temperature_k(CellId(c as u64)).unwrap(), 300.0 + expected / 2.0);
    }
    assert_eq!(source_history.samples().len(), 4);
    for sample in source_history.samples() {
        assert_eq!(sample.cell, CellId(sample.voxel as u64));
        assert_eq!(sample.phase.phase(), SolidLiquidPhase::Liquid);
        near(sample.phase.specific_enthalpy_j_kg(), h[sample.voxel]);
        for name in [DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY,
            THERMAL_CONDUCTIVITY_PROPERTY, DYNAMIC_VISCOSITY_PROPERTY]
        {
            let receipt = &sample.liquid.material().property(name).unwrap().answer().receipt;
            assert_eq!(receipt.selected, liquid.reference_state().material()
                .property(name).unwrap().answer().receipt.selected);
            assert_eq!(receipt.axis_quantities, *liquid.reference_point().axis_quantities());
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
            assert!(receipt.query_point.contains(&("lot".into(), 7.0)));
        }
    }
    near(accepted.advection().audit().external_enthalpy_j, -exported_j);
    near(accepted.advection().audit().inlet_mass_kg, 0.5);
    near(accepted.advection().audit().outlet_mass_kg, 0.5);
    near(audit.external_heat_j, source_history.inner().supplied_j - exported_j);
    near(next.energy_j() - old.energy_j(), audit.external_heat_j);
    near(next.mass_kg(), old.mass_kg());
    assert_eq!(observer.advection().audit().substeps, 0);
    assert_eq!(observer.advection().inner().inner(), &Heater::default());
    assert!(observer.advection().inner().samples().is_empty());
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
}

#[test]
fn changing_only_source_viscosity_changes_native_fan_flow_and_heat_transport() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut flow_config = SimpleConfig::new(faces);
    flow_config.tolerance = 1e-13;
    flow_config.pressure_tolerance = 1e-12;
    flow_config.internal_fans.push(InternalFan {
        patch: FacePatch { axis: 0, index: 2, lo: [0, 0, 0], hi: [4, 1, 1] },
        blows_positive: true,
        curve: FanCurve::new(&[(0.0, 8.0), (100.0, 8.0)]).unwrap(),
    });
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let mut outputs = Vec::new();
    for mu in [1.0, 4.0] {
        let card = card(claims(1.0, mu, false));
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
        let solved = fluid.solve_flow(&domain, &flow_config, &gate).unwrap();
        let old = initial(&domain, 1.0, &caloric);
        let observer = solved.advect_liquid_enthalpy(
            &old, input(&ids, &liquid, &caloric), Heater::default(), &gate,
        ).unwrap();
        let (next, accepted, audit) = old.advance_observed(&observer, 1.0, config(), || Ok(())).unwrap();
        let q = solved.flow().field.flux_x(2, 0, 0);
        assert!(q > 0.0);
        for edge in accepted.advection().flows() { near(edge.kg_s, q); }
        for sample in accepted.advection().inner().samples() {
            near(sample.liquid.dynamic_viscosity_pa_s(), mu);
        }
        near(next.energy_j() - old.energy_j(), audit.external_heat_j);
        near(audit.external_heat_j, 100.0 + accepted.advection().audit().external_enthalpy_j);
        near(next.mass_kg(), old.mass_kg());
        outputs.push((q, next.cells()[&CellId(1)].specific_h_j_kg));
    }
    assert!(outputs[1].0 < outputs[0].0);
    assert!(outputs[1].1 > outputs[0].1);
}

#[test]
fn source_boundary_refuses_a_claim_switch_and_rolls_back_the_complete_candidate() {
    let card = card(claims(2.0, 0.1, true));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    // A different same-name viscosity claim really is available above the boundary.
    let other = LiquidTemperatureLaw::try_new(&card, &point(330.0, 100000.0), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap();
    near(other.reference_state().dynamic_viscosity_pa_s(), 0.4);
    assert!(liquid.at_temperature(Temperature::new(330.0)).is_err());
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let old = initial(&domain, 2.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let observer = solved.advect_liquid_enthalpy(
        &old, input(&ids, &liquid, &caloric), Heater::default(), &gate,
    ).unwrap();
    let (_, short, _) = old.advance_observed(&observer, 0.25, config(), || Ok(())).unwrap();
    let durations = &short.advection().inner().inner().durations;
    assert!(!durations.is_empty());
    near(durations.iter().sum(), 0.25);
    assert!(matches!(
        old.advance_observed(&observer, 1.0, config(), || Ok(())),
        Err(LiquidEnthalpyObserverError::Binding(LiquidEnthalpyError::Material(_)))
    ));
    assert_eq!(observer.advection().audit().substeps, 0);
    assert_eq!(observer.advection().inner().inner(), &Heater::default());
    assert!(observer.advection().inner().samples().is_empty());
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
}

#[test]
fn inner_history_refusal_and_cancellation_preserve_sources_flow_and_thermal_state() {
    let card = card(claims(2.0, 0.1, false));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let old = initial(&domain, 2.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let history = Heater { reject_step: Some(2), ..Heater::default() };
    let observer = solved.advect_liquid_enthalpy(
        &old, input(&ids, &liquid, &caloric), history.clone(), &gate,
    ).unwrap();
    assert!(matches!(
        old.advance_observed(&observer, 1.0, config(), || Ok(())),
        Err(LiquidEnthalpyObserverError::Observer(DepositionError::InvalidInput(_)))
    ));
    assert!(matches!(
        old.advance_observed(&observer, 1.0, config(), || Err(DepositionError::Cancelled)),
        Err(LiquidEnthalpyObserverError::Binding(
            LiquidEnthalpyError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(observer.advection().inner().inner(), &history);
    assert!(observer.advection().inner().samples().is_empty());
    assert_eq!(observer.advection().audit().substeps, 0);
    assert!(core::ptr::eq(observer.material(), liquid.reference_state()));
    near(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
}

#[test]
fn same_card_cannot_hide_different_caloric_context_pins_or_borrowed_law() {
    let mut claims = claims(2.0, 0.1, false);
    claims.insert_claim(scalar("alternate_cp", SPECIFIC_HEAT_CAPACITY_DIMS, 4.0)).unwrap();
    let card = card(claims);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let high_pressure = source_caloric(&card, 200000.0, SPECIFIC_HEAT_CAPACITY_PROPERTY,
        DENSITY_PROPERTY, SolidLiquidPhase::Liquid);
    let alternate_cp = source_caloric(&card, 100000.0, "alternate_cp",
        DENSITY_PROPERTY, SolidLiquidPhase::Liquid);
    let equal_copy = caloric.clone();
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let old = initial(&domain, 2.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    assert!(matches!(solved.advect_liquid_enthalpy(
        &old, input(&ids, &liquid, &high_pressure), (), &gate,
    ), Err(LiquidEnthalpyError::Invalid(
        "caloric source context or temperature axis differs from the liquid law"
    ))));
    assert!(matches!(solved.advect_liquid_enthalpy(
        &old, input(&ids, &liquid, &alternate_cp), (), &gate,
    ), Err(LiquidEnthalpyError::Invalid(
        "caloric Cp/density source pins differ from the liquid law"
    ))));
    // Even an equal chart is not the actual borrowed caloric object in storage.
    assert_eq!(equal_copy, caloric);
    assert!(matches!(solved.advect_liquid_enthalpy(
        &old, input(&ids, &liquid, &equal_copy), (), &gate,
    ), Err(LiquidEnthalpyError::Thermal(_))));
    // Same numerical coefficients and card, but a different momentum context.
    let other_liquid = LiquidTemperatureLaw::try_new(&card, &point(300.0, 200000.0), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert!(matches!(solved.advect_liquid_enthalpy(
        &old, input(&ids, &other_liquid, &caloric), (), &gate,
    ), Err(LiquidEnthalpyError::Invalid(
        "momentum and thermal binding must retain the exact liquid reference state"
    ))));
}

#[test]
fn complete_caloric_chart_rejects_hidden_density_excursion_and_solid_phase() {
    let mut claims = claims(2.0, 0.1, false);
    let mut varying = scalar("alternate_density", Density::DIMS, 2.0);
    varying.value = PropertyValue::Curve {
        abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 2.0), (500.0, 1.6), (810.0, 2.0)], dims: Density::DIMS,
    };
    varying.interpolation = InterpolationPolicy::LinearInside;
    claims.insert_claim(varying).unwrap();
    let card = card(claims);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let varying_density = source_caloric(&card, 100000.0, SPECIFIC_HEAT_CAPACITY_PROPERTY,
        "alternate_density", SolidLiquidPhase::Liquid);
    let solid = source_caloric(&card, 100000.0, SPECIFIC_HEAT_CAPACITY_PROPERTY,
        DENSITY_PROPERTY, SolidLiquidPhase::Solid);
    let knots = varying_density.curve().knots();
    near(knots.first().unwrap().bulk_density_kg_m3, 2.0);
    near(knots.last().unwrap().bulk_density_kg_m3, 2.0);
    assert!(knots.iter().any(|knot| knot.bulk_density_kg_m3 < 1.9));
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let old = initial(&domain, 2.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    for incompatible in [&varying_density, &solid] {
        assert!(matches!(solved.advect_liquid_enthalpy(
            &old, input(&ids, &liquid, incompatible), (), &gate,
        ), Err(LiquidEnthalpyError::Invalid(
            "the complete caloric span must be fully liquid at the fixed source density"
        ))));
    }
}

#[test]
fn frozen_liquid_flow_refuses_rate_bounded_thermal_sources() {
    for bounded in [SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY] {
        let mut source_claims = ClaimSet::new();
        for (name, dims, value) in [
            (DENSITY_PROPERTY, Density::DIMS, 2.0),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 2.0),
            (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 0.01),
            (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 0.1),
        ] {
            let mut claim = scalar(name, dims, value);
            if name == bounded {
                claim.validity = claim.validity.with("shear_rate", 0.0, 10.0);
            }
            source_claims.insert_claim(claim).unwrap();
        }
        let card = card(source_claims);
        let lower = point(300.0, 100000.0).with("shear_rate", 1.0).unwrap();
        let upper = point(810.0, 100000.0).with("shear_rate", 1.0).unwrap();
        let liquid = LiquidTemperatureLaw::try_new(
            &card, &lower, "T", LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
        let rho = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &cp,
                density: &rho,
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Liquid, 0.0, 1e-8,
        ).unwrap();

        // These are valid instantaneous source queries. The coupled owner must
        // still refuse to freeze a thermal coefficient bounded by shear rate,
        // because this reduction evolves temperature and does not evolve rate.
        let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
        let gate = CancelGate::new();
        let fluid = MaterialChtFluid::try_new(liquid.reference_state()).unwrap();
        let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
        let old = initial(&domain, 2.0, &caloric);
        let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
        let result = solved.advect_liquid_enthalpy(
            &old, input(&ids, &liquid, &caloric), Heater::default(), &gate,
        );
        assert!(
            matches!(result, Err(LiquidEnthalpyError::Material(
                fs_material::state_point::MaterialStatePointError::InvalidRequirement {
                    ref property, ..
                }
            )) if property == bounded),
            "{bounded}: the frozen-flow adapter must refuse a shear-rate-bound thermal source",
        );
    }
}

