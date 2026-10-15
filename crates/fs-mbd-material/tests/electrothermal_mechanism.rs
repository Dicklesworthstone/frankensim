#![cfg(feature = "electrothermal-coupling")]
//! Actual spatial current, shared material heat history, and articulated dynamics.
#[path = "heating/common.rs"]
mod common;
use common::{budget, caloric, close, config, point};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{Se3, Vec3, Wrench};
use fs_lbm::conjugate::{ChtError, ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::enthalpy::transport::{EnthalpyTransport, EnthalpyTransportSample,
    advance_joule_enthalpy_with_material_transport};
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point};
use fs_mbd::articulated::{BaseState, JointModel, forward_dynamics};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::{HeatedArticulated, HeatedBody, HeatingError, LinkHeating,
    ThermalRegion};
use fs_qty::Temperature;
use fs_thermal::enthalpy::*;
use fs_thermal::observers::ObserverCouplingError;

const IDS: [CellId; 3] = [CellId(0), CellId(1), CellId(2)];
const LINKS: [u64; 2] = [70, 71];
const KEYS: [&str; 3] = ["left", "middle", "right"];

fn card(cp: f64, resistivity: f64, alpha_upper: f64) -> MaterialCard {
    // All five laws belong to ONE immutable card, not separate solver presets.
    // Synthetic coefficients deliberately expose the coupling, not named metals.
    let base = common::card(cp, 2e-3, alpha_upper);
    let mut claims = base.claims().clone();
    for (name, dims, low, high) in [
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, resistivity, 3. * resistivity),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 2., 6.),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300., low), (500., high)], dims },
            validity: ValidityDomain::unconstrained().with("T", 300., 500.).with("pressure", 1e5, 2e5),
            interpolation: InterpolationPolicy::LinearInside, uncertainty: UncertaintyModel::Unstated,
            observations: vec![], provenance: Provenance { source: "synthetic electrothermal fixture".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn geometry(slot: usize) -> VolumeMoments {
    let cube = VolumeMoments::cuboid(Vec3::new(1., 1., 1.)).unwrap();
    VolumeMoments::new(1., Vec3::new(slot as f64 + 0.5, 0.5, 0.5), cube.inertia_m5()).unwrap()
}
fn mechanism(card: &MaterialCard) -> ResolvedArticulated {
    let p = point(300., 1e5);
    let mass = |slot| RegionBinding { key: KEYS[slot], geometry: geometry(slot),
        material: card, point: &p, selection: MaterialPropertySelection::SingleClaimOnly };
    let mount = [mass(0), mass(2)];
    let rotor = [mass(1)];
    resolve_articulated(&[
        MaterialLink { name: "mount", parent: None, parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED, regions: &mount },
        MaterialLink { name: "rotor", parent: Some(0), parent_from_child_zero: Se3::identity(),
            joint: JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(), regions: &rotor },
    ], 3, || false).unwrap()
}
fn region<'a>(card: &'a MaterialCard, chart: &'a ResolvedHeatCapacityCurve, slot: usize) -> ThermalRegion<'a> {
    let mut r = common::region(card, chart, KEYS[slot], slot as u64, 0.);
    // This is the electrical reference grid; each body uses these coordinates
    // at q=0. Changing q during the electrical solve is outside this fixture.
    r.geometry.start_m = [slot as f64, 0.5, 0.5];
    r.geometry.end_m = [slot as f64 + 1., 0.5, 0.5];
    r
}
fn contacts(conductance: f64) -> Vec<TimedEvent> {
    (0..2).map(|i| TimedEvent { time_s: 0., event: Event::Contact(ThermalLink {
        id: LINKS[i], a: IDS[i], b: IDS[i+1], area_m2: 1., conductance_w_k: conductance,
    }) }).collect()
}
fn heated<'a>(r: &'a ResolvedArticulated, c: &'a MaterialCard, h: &'a ResolvedHeatCapacityCurve) -> HeatedArticulated<'a> {
    HeatedArticulated::new(r, vec![
        LinkHeating { link: "rotor", regions: vec![region(c, h, 1)] },
        LinkHeating { link: "mount", regions: vec![region(c, h, 2), region(c, h, 0)] },
    ], contacts(0.), budget(), config(), || Ok(())).unwrap()
}
fn setup(voltage: f64) -> ElectricSetup {
    ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo: [0,0,0], hi: [1,1,1], drive: ElectrodeDrive::Voltage(voltage) },
        Electrode { lo: [2,0,0], hi: [3,1,1], drive: ElectrodeDrive::Voltage(0.) },
    ] }
}
fn input<'a>(domain: &'a VoxelDomain, setup: &'a ElectricSetup) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput { domain, setup, contacts: &[], cell_ids: &IDS, max_cells: 3,
        max_current_solves: 1001, power_tolerance_w: 1e-8, energy_tolerance_j: 1e-8 }
}
fn resolve(card: &MaterialCard, t: f64) -> Result<EnthalpyTransport<ResolvedMaterialStatePoint>, JouleEnthalpyError> {
    let required: Vec<_> = [(ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS)].into_iter()
        .map(|(key, dims)| ScalarPropertyRequirement::try_new(key, dims,
            ScalarAdmissibility::StrictlyPositive).unwrap()).collect();
    let material = resolve_material_state_point(card, &point(t, 1e5), &required,
        MaterialPropertySelection::SingleClaimOnly).map_err(|source|
            ChtError::InvalidInput { field: "fixture material", reason: source.to_string() })?;
    let rho_e = material.property(ELECTRICAL_RESISTIVITY_PROPERTY).unwrap().value_si();
    let k = material.property(THERMAL_CONDUCTIVITY_PROPERTY).unwrap().value_si();
    EnthalpyTransport::try_new(material, Some(rho_e), k)
}
fn link(patch: &ContactPatch, sides: [&EnthalpyTransportSample<ResolvedMaterialStatePoint>; 2])
    -> Result<(f64, (f64, f64)), JouleEnthalpyError> {
    let ka = sides[0].response.conductivity_w_m_k();
    let kb = sides[1].response.conductivity_w_m_k();
    // Explicit metre-wide reference voxels, perfect thermal interface. Keep
    // BOTH same-epoch coefficients as the concrete link payload in the report.
    Ok((patch.link.area_m2 / (0.5 / ka + 0.5 / kb), (ka, kb)))
}
fn acceleration(body: &HeatedArticulated<'_>) -> f64 {
    forward_dynamics(body.configuration().model(), BaseState::stationary(Se3::identity()),
        &[0.], &[0.], &[1.], Vec3::default(), &[Wrench::default(); 2]).unwrap().generalized_acceleration[0]
}

#[test]
fn actual_current_heat_and_transport_use_the_same_material_and_existing_dynamics() {
    let c = card(10., 5., 500.);
    let h = caloric(&c, 1e5, SolidLiquidPhase::Solid);
    let r = mechanism(&c);
    let body = heated(&r, &c, &h);
    let domain = VoxelDomain::new(3,1,1,1.).unwrap();
    let drive = setup(20.);
    let gate = CancelGate::new();
    let (plain, report) = advance_joule_enthalpy_with_material_transport(body.thermal(), 2., config(),
        input(&domain, &drive), &LINKS, &gate, |_,_,t| resolve(&c,t), link).unwrap();
    let (next, coupled) = body.advance_joule(2., config(), input(&domain, &drive), &LINKS,
        &gate, |_,_,t| resolve(&c,t), link).unwrap();
    assert_eq!(next.thermal().cells(), plain.cells());
    assert_eq!(coupled.joule, report.joule);
    assert_eq!(coupled.endpoint, report.endpoint);
    assert_eq!(coupled.last_applied, report.last_applied);
    close(coupled.joule.thermal.external_heat_j, coupled.joule.deposited_joule_j);
    close(next.thermal().energy_j(), coupled.joule.electrode_energy_j);
    assert_eq!(coupled.joule.current_solves, coupled.joule.thermal.substeps + 1);
    assert!(next.thermal().temperature_k(IDS[0]).unwrap() > 300.); // cross-link heat, not electrode bulk heat
    let t = next.thermal().temperature_k(IDS[1]).unwrap();
    let s = 1. + 2e-3 * (t - 300.);
    // rho=2, unit cube COM (1.5,.5,.5): Izz at the hinge = 16/3 kg m2.
    close(acceleration(&next), 3. / (16. * s * s));
    assert!(acceleration(&next) < acceleration(&body));
    assert!(coupled.joule.endpoint.electric.electrodes[0].current_a < 4.);
    for sample in &coupled.endpoint.cells {
        assert_eq!(sample.response.material().card_identity(), c.content_hash());
        let alpha = 2e-3 * (sample.temperature_k - 300.);
        let body_index = if sample.cell.id == IDS[1] { 1 } else { 0 };
        let region_index = if sample.cell.id == IDS[2] { 1 } else { 0 };
        let region = &next.configuration().bodies()[body_index].regions()[region_index];
        close(region.expansion().free_linear_strain(), alpha);
        assert_eq!(sample.cell.mass_kg.to_bits(), region.inertia().mass().to_bits());
    }
    close(next.thermal().mass_kg(), 6.);
    assert_eq!(body.thermal().time_s(), 0.);
}

#[test]
fn cp_only_data_swap_and_source_shutoff_propagate_without_resetting_heat_history() {
    let domain = VoxelDomain::new(3,1,1,1.).unwrap();
    let drive = setup(20.);
    let off = setup(0.); // shorted ideal electrodes, NOT an open-circuit claim
    let a = card(10., 5., 500.);
    let b = card(20., 5., 500.);
    let ca = caloric(&a, 1e5, SolidLiquidPhase::Solid);
    let cb = caloric(&b, 1e5, SolidLiquidPhase::Solid);
    let ra = mechanism(&a);
    let rb = mechanism(&b);
    assert_eq!(ra.model(), rb.model());
    let old_a = heated(&ra, &a, &ca);
    let old_b = heated(&rb, &b, &cb);
    let gate = CancelGate::new();
    let (hot_a, report_a) = old_a.advance_joule(2., config(), input(&domain, &drive), &LINKS,
        &gate, |_,_,t| resolve(&a,t), link).unwrap();
    let (hot_b, report_b) = old_b.advance_joule(2., config(), input(&domain, &drive), &LINKS,
        &gate, |_,_,t| resolve(&b,t), link).unwrap();
    assert!(hot_a.thermal().temperature_k(IDS[1]).unwrap() > hot_b.thermal().temperature_k(IDS[1]).unwrap());
    assert!(acceleration(&hot_a) < acceleration(&hot_b));
    assert!(report_a.joule.endpoint.electric.electrodes[0].current_a
        < report_b.joule.endpoint.electric.electrodes[0].current_a);
    let (cool, stopped) = hot_a.advance_joule(4., config(), input(&domain, &off), &LINKS,
        &gate, |_,_,t| resolve(&a,t), link).unwrap();
    close(stopped.joule.deposited_joule_j, 0.);
    close(stopped.joule.electrode_energy_j, 0.);
    close(cool.thermal().energy_j(), hot_a.thermal().energy_j());
    assert!(cool.thermal().temperature_k(IDS[1]).unwrap() < hot_a.thermal().temperature_k(IDS[1]).unwrap());
    assert!(acceleration(&cool) > acceleration(&hot_a));
    assert!(cool.thermal().temperature_k(IDS[0]).unwrap() > hot_a.thermal().temperature_k(IDS[0]).unwrap());
}

#[test]
fn late_electrical_support_mechanical_support_and_cancellation_publish_nothing() {
    let c = card(10., 5., 500.);
    let h = caloric(&c, 1e5, SolidLiquidPhase::Solid);
    let r = mechanism(&c);
    let body = heated(&r, &c, &h);
    let domain = VoxelDomain::new(3,1,1,1.).unwrap();
    let drive = setup(20.);
    let gate = CancelGate::new();
    let cfg = StepConfig { max_dt_s: 0.25, energy_tolerance_j: 1e-9 };
    // No transport: exactly one start solve succeeds and heats the centre.
    // Its geometry is configured before the final current-support check fails.
    let mut request = input(&domain, &drive);
    request.max_current_solves = 1;
    let result = body.advance_joule(0.25, cfg, request, &[], &gate, |_,_,t| resolve(&c,t), link);
    assert!(matches!(result, Err(ObserverCouplingError::First(JouleEnthalpyError::CurrentBudget))));
    let result = body.advance_joule(0.25, cfg, input(&domain,&drive), &[], &gate,
        |_,_,t| if t > 300. { Err(JouleEnthalpyError::Invalid("endpoint source unavailable")) }
            else { resolve(&c,t) }, link);
    assert!(matches!(result, Err(ObserverCouplingError::First(JouleEnthalpyError::Invalid(_)))));
    let limited = card(10., 5., 300.5);
    let curve = caloric(&limited, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&limited);
    let narrow = heated(&reference, &limited, &curve);
    assert!(matches!(narrow.advance_joule(0.25, cfg, input(&domain,&drive), &[], &gate,
        |_,_,t| resolve(&limited,t), link), Err(ObserverCouplingError::Second(HeatingError::Mechanical(_)))));
    gate.request();
    assert!(body.advance_joule(0.25, cfg, input(&domain,&drive), &[], &gate,
        |_,_,t| resolve(&c,t), link).is_err());
    for unchanged in [&body, &narrow] {
        assert_eq!(unchanged.thermal().time_s(), 0.);
        assert_eq!(unchanged.thermal().energy_j(), 0.);
        assert_eq!(unchanged.configuration().model(), unchanged.reference().model());
    }
}

#[test]
fn one_body_electrothermal_adapter_retains_the_original_thermal_solution() {
    let c = card(10., 5., 500.);
    let h = caloric(&c, 1e5, SolidLiquidPhase::Solid);
    let p = point(300., 1e5);
    let regions: Vec<_> = (0..3).map(|slot| RegionBinding { key: KEYS[slot], geometry: geometry(slot),
        material: &c, point: &p, selection: MaterialPropertySelection::SingleClaimOnly }).collect();
    let reference: ResolvedBody = resolve_body(&regions, 3, || false).unwrap();
    let body = HeatedBody::new(&reference, (0..3).map(|slot| region(&c, &h, slot)).collect(),
        contacts(0.), budget(), config(), || Ok(())).unwrap();
    let domain = VoxelDomain::new(3,1,1,1.).unwrap();
    let drive = setup(20.);
    let gate = CancelGate::new();
    let (plain, expected) = advance_joule_enthalpy_with_material_transport(body.thermal(), 1., config(),
        input(&domain,&drive), &LINKS, &gate, |_,_,t| resolve(&c,t), link).unwrap();
    let (heated, actual) = body.advance_joule(1., config(), input(&domain,&drive), &LINKS,
        &gate, |_,_,t| resolve(&c,t), link).unwrap();
    assert_eq!(plain.cells(), heated.thermal().cells());
    assert_eq!(actual.joule, expected.joule);
    assert_eq!(actual.endpoint, expected.endpoint);
    assert_eq!(heated.configuration().inertia().mass().to_bits(), reference.inertia().mass().to_bits());
    assert!(heated.configuration().regions()[1].geometry().volume_m3() > 1.);
}
