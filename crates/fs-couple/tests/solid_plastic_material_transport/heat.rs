//! G0/G3: source-backed radiation/contact data reach loaded native J2 tetrahedra.
#[path = "heat/solid.rs"]
mod solid;
#[path = "heat/electrical.rs"]
mod electrical;
#[path = "heat/support.rs"]
mod support;

use super::*;
use fs_conduction::{GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity,
    ViewFactorMatrix, SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as CONTACT_R,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as CONTACT_DIMS};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_couple::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, MaterialHeatNetwork, PreparedHeatNetwork};
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_thermal::observers::ObserverCouplingError;

fn fixture(emissivity: f64) -> Fixture {
    fixture_with_property(emissivity, None)
}

fn fixture_with_property(emissivity: f64, extra: Option<PropertyClaim>) -> Fixture {
    let old = Fixture::new([5.0; 2], 0.01, false);
    let cards = old.cards.map(|card| {
        let mut claims = card.claims().clone();
        claims.insert_claim(scalar(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, emissivity)).unwrap();
        if let Some(property) = &extra { claims.insert_claim(property.clone()).unwrap(); }
        MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap()
    });
    // The augmented card is the authority for BOTH caloric and mechanical data.
    let lower = old.point.clone().with("T", 280.0).unwrap();
    let upper = old.point.clone().with("T", 320.0).unwrap();
    let cp = PropertyKey::new(CP, CP_DIMS);
    let rho = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    let caloric = [0, 1].map(|i| ResolvedHeatCapacityCurve::try_new(&cards[i],
        HeatCapacitySourceQuery { heat_capacity: &cp, density: &rho,
            temperature_axis: "T", lower: &lower, upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly },
        SolidLiquidPhase::Solid, -100.0, 1e-4).unwrap());
    Fixture { cards, caloric, point: old.point, conductivity: old.conductivity }
}

fn with_heat<R>(emissivity: f64, resistance: f64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &Fixture, &DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>, &PreparedHeatNetwork<'_, '_>) -> R) -> R {
    with_heat_fixture(fixture(emissivity), resistance, |cx, gate, f, thermal, solid, network| {
        let prepared = network.prepare(cx, thermal).unwrap();
        run(cx, gate, f, thermal, solid, &prepared)
    })
}

fn with_heat_fixture<R>(f: Fixture, resistance: f64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &Fixture, &DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>, &MaterialHeatNetwork<'_>) -> R) -> R {
    with_cx(|gate, cx| {
        let thermal = f.thermal();
        let solid = f.solid(&thermal, cx);
        let keys = vec![Some((f.conductivity.clone(), MaterialPropertySelection::SingleClaimOnly)); 2];
        let cells = solid.thermal_material_cells(&keys).unwrap();
        let links = links();
        let bulk = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
        let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &[] };
        let mut claims = ClaimSet::new();
        claims.insert_claim(scalar(CONTACT_R, CONTACT_DIMS, resistance)).unwrap();
        let interface = InterfaceSystemCard::assemble(
            SurfaceSpec { material: f.cards[0].id().clone(), texture_frame: "a".into() },
            SurfaceSpec { material: f.cards[1].id().clone(), texture_frame: "b".into() },
            SystemContext { medium: "fixed synthetic interface".into(), third_body: None,
                environment: "test".into(), history: "fixed".into() }, claims, vec![]).unwrap();
        let point = QueryPoint::new().with("Ta", 300.0).unwrap().with("Tb", 300.0).unwrap();
        let contacts = [MaterialContactLink { link: 91,
            sides: [0, 1].map(|i| ContactCell { cell: CELLS[i], material_id: MATS[i], curve: f.caloric[i].curve() }),
            point: ThermalInterfacePoint { interface: &interface, materials: [&f.cards[0], &f.cards[1]],
                context: &point, temperature_axes: ["Ta", "Tb"], liquid_fraction_axes: [None, None],
                mean_temperature_axis: None }, selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
        let surfaces = [0, 1].map(|i| RadiatingCell { surface: i, cell: CELLS[i], material_id: MATS[i],
            curve: f.caloric[i].curve(), material: cells[i].point, selection: MaterialPropertySelection::SingleClaimOnly });
        let (complex, positions) = box_grid([1, 1, 1], [L, L, L]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let traces = [0, 1].map(|i| {
            let phase = f.caloric[i].state_at_temperature(300.0).unwrap();
            let e = SurfaceEmissivity::from_thermal_state("synthetic", cells[i].point, phase,
                MaterialPropertySelection::SingleClaimOnly).unwrap().0;
            RadiationSurface::new(&mesh, format!("surface-{i}"),
                |face| on_box_face(face.centroid[0], i as f64 * L), e).unwrap()
        });
        // Algebraic parallel-plate fixture, not a reconstructed tetrahedral enclosure.
        let enclosure = GrayDiffuseEnclosure::new(traces.into_iter().collect(),
            ViewFactorMatrix::infinite_parallel_plates(L * L).unwrap()).unwrap();
        let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
            power_tolerance_w: 1e-11, radiosity_tolerance_w_m2: 1e-8 };
        let network = MaterialHeatNetwork { electrical: &electrical,
            thermal_contacts: &contacts, radiation: Some(&radiation) };
        run(cx, gate, &f, &thermal, &solid, &network)
    })
}

// Independent composition of the original evaluators through the same public
// observer protocol. It does not implement an alternative physical model.
#[derive(Clone)]
struct Manual<'p, 'n, 'a, 'scope> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    cx: &'p Cx<'scope>,
}
impl DepositionObserver for Manual<'_, '_, '_, '_> {
    type Error = HeatNetworkError;
    fn sources(&mut self, t: &DepositionSimulation<'_>, check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?; Ok(self.prepared.resolve(self.cx, t)?.bulk.sources_w)
    }
    fn conductances(&mut self, t: &DepositionSimulation<'_>, check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> {
        check()?; Ok(self.prepared.resolve(self.cx, t)?.bulk.conductances_w_k)
    }
    fn limit_timestep_s(&mut self, t: &DepositionSimulation<'_>, proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<f64, Self::Error> {
        check()?;
        let source = self.prepared.resolve(self.cx, t)?.bulk.source;
        let Some(radiation) = source.radiation else { return Ok(proposed_s); };
        let bounds = radiation.conductance_bounds_w_k(&mut *check)?;
        Ok(t.limit_timestep_for_additional_conductances(proposed_s, &bounds, check)?)
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>, _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> { check()?; Ok(()) }
}

#[test]
fn complete_heat_observer_replays_native_loaded_solid_and_keeps_endpoint_sources() {
    with_heat(0.7, 0.1, |cx, _, f, thermal, solid, prepared| {
        let observer = prepared.observe(cx, Recorder::default(), 1e-12).unwrap();
        let actual = solid.advance(thermal, &observer, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        let manual = solid.advance(thermal, &Manual { prepared, cx }, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        assert_eq!(actual.thermal().cells(), manual.thermal().cells());
        assert_eq!(actual.solid().solution().displacement_m, manual.solid().solution().displacement_m);
        assert_eq!(actual.solid().solution().element_states, manual.solid().solution().element_states);
        let sample = actual.observer().first.sample().unwrap();
        let exact = prepared.resolve(cx, actual.thermal()).unwrap().bulk.source;
        assert_eq!(sample.thermal_contacts, exact.thermal_contacts);
        assert_eq!(sample.radiation, exact.radiation);
        for surface in sample.radiation.as_ref().unwrap().surfaces() {
            let side = CELLS.iter().position(|id| *id == surface.cell).unwrap();
            assert_eq!(surface.material.card_identity(), f.cards[side].content_hash());
            assert_eq!(surface.phase.temperature_k(), actual.thermal().temperature_k(surface.cell).unwrap());
        }
        let net = actual.observer().first.closed_radiation_net_j();
        close(actual.report().external_heat_j - net, 0.0, 1e-10);
        assert!((actual.report().combined_residual_j + net).abs() < 2.0 * control().energy_tolerance_j);
        close(actual.observer().second.durations.iter().sum(), 1.0, 1e-12);
        assert!(actual.solid().solution().element_states[0].alpha > 0.0);
        assert_eq!(actual.solid().solution().element_states[1].alpha, 0.0);
        assert!(actual.report().numerical_loss_j > actual.report().physical_heat_j);
        assert!(observer.first.sample().is_none());
        assert!(observer.second.durations.is_empty());
        assert_eq!(solid.time_s(), 0.0); assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn emissivity_and_interface_data_change_yield_history_and_free_3d_displacement() {
    let run = |e, r| with_heat(e, r, |cx, _, _, t, solid, p| {
        let observer = p.observe(cx, (), 1e-12).unwrap();
        let next = solid.advance(t, &observer, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx).unwrap();
        [next.thermal().temperature_k(CELLS[0]).unwrap(), next.thermal().temperature_k(CELLS[1]).unwrap(),
            next.solid().solution().displacement_m[3][0], next.solid().solution().element_states[0].alpha]
    });
    let baseline = run(0.2, 0.01);
    for changed in [run(0.9, 0.01), run(0.2, 2.0)] {
        assert!((baseline[0] - changed[0]).abs() > 1e-4);
        assert!((baseline[1] - changed[1]).abs() > 1e-4);
        assert!((baseline[2] - changed[2]).abs() > 1e-10);
        assert!((baseline[3] - changed[3]).abs() > 1e-8);
    }
}

#[derive(Clone)]
struct CompetingConductance;
impl DepositionObserver for CompetingConductance {
    type Error = DepositionError;
    fn conductances(&mut self, _: &DepositionSimulation<'_>, _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> { Ok(BTreeMap::from([(91, 0.0)])) }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>, _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> { Ok(()) }
}

#[test]
fn competing_conductance_and_scope_cancellation_preserve_the_original_solid_and_heat() {
    with_heat(0.7, 0.1, |cx, gate, _, t, solid, p| {
        let original = t.cells().clone();
        let observer = p.observe(cx, CompetingConductance, 1e-12).unwrap();
        assert!(matches!(solid.advance(t, &observer, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx),
            Err(ThermoplasticError::Observer(ObserverCouplingError::ConflictingConductance(91)))));
        let observer = p.observe(cx, (), 1e-12).unwrap();
        gate.request();
        assert!(matches!(solid.advance(t, &observer, 1.0, &supports(SHEAR), &applied(),
            IncrementSettings::default(), thermal_config(), control(), cx), Err(ThermoplasticError::Cancelled)));
        assert_eq!(t.cells(), &original); assert_eq!(solid.time_s(), 0.0);
        assert!(observer.first.sample().is_none());
        assert_eq!(observer.first.closed_radiation_net_j(), 0.0);
    });
}

#[test]
fn richer_thermal_checkpoint_errors_and_invalid_energy_allowances_are_not_erased() {
    with_heat(0.7, 0.1, |cx, _, _, t, _, p| {
        for allowance in [-1.0, f64::NAN, f64::INFINITY] { assert!(p.observe(cx, (), allowance).is_err()); }
        let mut observer = p.observe(cx, (), 1e-12).unwrap();
        assert!(matches!(observer.first.sources(t, &mut || Err(DepositionError::Budget("test checkpoint"))),
            Err(HeatNetworkError::Thermal(DepositionError::Budget("test checkpoint")))));
        assert!(observer.first.sample().is_none());
    });
}

#[path = "../bridge_heat_network/common.rs"]
mod powered;
#[test]
fn thermal_observation_refuses_to_drop_an_actual_electrical_network() {
    powered::with_case(0.7, 0.1, true, false, vec![], |_, cx, _, network, t, _, _| {
        let prepared = network.prepare(cx, t).unwrap();
        assert!(matches!(prepared.observe(cx, (), 1e-12), Err(HeatNetworkError::Invalid(_))));
    });
}
