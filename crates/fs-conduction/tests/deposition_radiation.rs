//! Actual shared thermal evolution under source-resolved enclosure radiation.
mod support;
#[path="radiation_material/common.rs"] mod common;
use common::*;
use fs_conduction::transient::deposition::*;
use fs_conduction::transient::deposition::radiation::*;
use fs_conduction::STEFAN_BOLTZMANN_W_M2_K4;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_solver::NewtonKrylovConfig;
use fs_matdb::{MaterialCard, QueryPoint};
use std::collections::BTreeMap;

use fs_blake3::hash_bytes;
use fs_conduction::{EMISSIVITY_DIMS, GrayDiffuseEnclosure, SURFACE_EMISSIVITY_PROPERTY};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{Dims, Pressure};
use fs_thermal::enthalpy::conduction::implicit::coupled_source::ImplicitCoupledHeatSourceLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductivityLaw};

fn id(i:usize)->CellId { CellId(if i==0 {30} else {10}) }
fn initial<'a>(curves:&'a [EquilibriumEnthalpyPhaseCurve;2],h:[f64;2],power:f64)->DepositionSimulation<'a> {
    let events=[0,1].map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:id(i),segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),
        mass_kg:1.0,specific_h_j_kg:h[i],geometry:CellGeometry {
            start_m:[i as f64,0.0,0.0],end_m:[i as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:1.5},
        boundary:Boundary {ambient_k:350.0,htc_w_m2_k:0.0,reservoir:None,source_w:if i==0 {power} else {0.0}},
    })}).to_vec();
    DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),
        (MaterialId(1),&curves[1] as &dyn ThermalLaw)],DepositionBudget {
        max_cells:2,max_links:0,max_events:2,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-8},||Ok(())).unwrap().0
}
fn bindings<'a>(cards:&'a [MaterialCard;2],curves:&'a [EquilibriumEnthalpyPhaseCurve;2],context:&'a QueryPoint)
    ->[RadiatingCell<'a>;2] {
    // Deliberately reversed input order; matrix order is not CellId order either.
    [1,0].map(|i|RadiatingCell {surface:i,cell:id(i),material_id:MaterialId(i as u64),curve:&curves[i],
        material:binding(&cards[i],context),selection:MaterialPropertySelection::SingleClaimOnly})
}
fn controls()->RadiationStepSolve {
    RadiationStepSolve {maximum_cells:2,maximum_evaluations:300,maximum_iterations:32,difference_step_j:1e-4,
        cell_tolerance_j:1e-7,total_tolerance_j:1e-7,newton:NewtonKrylovConfig {
            absolute_tolerance:1e-9,relative_tolerance:0.0,forcing_minimum:1e-4,linear_restart:4,
            ..NewtonKrylovConfig::default()}}
}
#[test]
fn actual_implicit_cell_heat_matches_independent_nonlinear_parallel_plate_root() {
    let cards=[card(0.4,0.4),card(0.6,0.6)]; let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context=context();let surfaces=bindings(&cards,&curves,&context);
    let enclosure=enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
    let model=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
    let old=initial(&curves,[150.0,-50.0],0.0);
    let (next,report)=support::with_cx(|cx|model.step_candidate(cx,&old,0.01,controls(),&[0.0;2]).unwrap());
    let (mut lo,mut hi)=(0.0,50.0);
    for _ in 0..60 {
        let q=0.5*(lo+hi);
        let radiation=STEFAN_BOLTZMANN_W_M2_K4*((450.0_f64-q).powi(4)-(350.0_f64+q).powi(4))/(1.0/0.4+1.0/0.6-1.0);
        if q-0.01*radiation<0.0 {lo=q;} else {hi=q;}
    }
    let q=0.5*(lo+hi);
    near(next.cells()[&id(0)].specific_h_j_kg,150.0-q);
    near(next.cells()[&id(1)].specific_h_j_kg,-50.0+q);
    near(next.energy_j(),old.energy_j());near(report.external_heat_j,0.0);
    assert!(report.maximum_cell_residual_j<=controls().cell_tolerance_j);
    assert_eq!(old.time_s(),0.0);assert_eq!(old.cells()[&id(0)].specific_h_j_kg,150.0);
}
#[test]
fn implicit_radiation_retains_latent_heat_with_unresolvable_solver_coordinates() {
    let cards=[card(0.4,0.4),card(0.6,0.6)];
    let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context=context();let surfaces=bindings(&cards,&curves,&context);
    let enclosure=enclosure([emissivity(&cards[0],&curves[0],50.0),emissivity(&cards[1],&curves[1],50.0)]);
    let model=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,
        power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
    let old=initial(&curves,[50.0,50.0],5.0);
    let (next,report)=support::with_cx(|cx|
        model.step_candidate(cx,&old,0.01,controls(),&[1e-20,-1e-20]).unwrap());
    // Both cells remain on the 400 K plateau. The actual source changes the
    // hot cell's phase inventory while the enclosure stays in thermal balance.
    near(next.cells()[&id(0)].specific_h_j_kg,50.05);
    near(next.cells()[&id(1)].specific_h_j_kg,50.0);
    assert_eq!(next.temperature_k(id(0)).unwrap(),400.0);
    assert_eq!(next.temperature_k(id(1)).unwrap(),400.0);
    near(report.thermal.stored_energy_change_j,0.05);
    near(report.external_heat_j,0.05);
    assert!(report.maximum_cell_residual_j<=controls().cell_tolerance_j);
    assert_eq!(old.time_s(),0.0);
    assert_eq!(old.cells()[&id(0)].specific_h_j_kg,50.0);
}
#[test]
fn emissivity_data_changes_phase_evolution_at_fixed_kelvin_and_preserves_independent_heat() {
    let run=|scale:f64| {
        let cards=[card(0.2*scale,0.8*scale),card(0.6,0.6)];
        let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
        let context=context();let surfaces=bindings(&cards,&curves,&context);
        let enclosure=enclosure([emissivity(&cards[0],&curves[0],50.0),emissivity(&cards[1],&curves[1],-50.0)]);
        let model=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
        let old=initial(&curves,[50.0,-50.0],5.0);let mut state=old.clone();
        support::with_cx(|cx| {
            for i in 1..=20 {
                let (next,r)=model.step_candidate(cx,&state,i as f64*0.002,controls(),&[0.0;2]).unwrap();
                near(r.external_heat_j,0.002*5.0);
                let hot=&r.radiation.surfaces()[0];
                assert_eq!(hot.phase.temperature_k(),400.0);
                near(r.radiation.enclosure().surfaces()[0].emissivity().value(),scale*(0.2+0.6*hot.phase.liquid_mass_fraction()));
                near(next.energy_j()-old.energy_j(),i as f64*0.002*5.0);
                state=next;
            }
        });
        (curves[0].state_at_specific_enthalpy(state.cells()[&id(0)].specific_h_j_kg).unwrap().liquid_mass_fraction(),
            state.temperature_k(id(1)).unwrap())
    };
    let low=run(0.5);let high=run(1.0);
    assert!(high.0<low.0 && high.1>low.1);
}
#[test]
fn complete_binding_summed_coverage_and_transactional_source_addition_are_enforced() {
    let cards=[card(0.4,0.4),card(0.6,0.6)];let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context=context();let mut surfaces=bindings(&cards,&curves,&context);
    let enclosure=enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
    let old=initial(&curves,[150.0,-50.0],0.0);
    fn make<'a>(enclosure: &'a fs_conduction::GrayDiffuseEnclosure, cells: &'a [RadiatingCell<'a>]) -> MaterialRadiation<'a> {
        MaterialRadiation {enclosure,cells,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8}
    }
    support::with_cx(|cx| {
        let source=make(&enclosure,&surfaces).resolve(cx,&old).unwrap();
        let mut map=BTreeMap::from([(id(0),5.0),(id(1),-2.0)]);
        source.add_to(&mut map).unwrap();near(map.values().copied().sum(),3.0);
        let mut bad=BTreeMap::from([(id(0),f64::INFINITY),(id(1),7.0)]);let saved=bad.clone();
        assert!(source.add_to(&mut bad).is_err());assert_eq!(bad,saved);
        assert!(make(&enclosure,&surfaces[..1]).resolve(cx,&old).is_err());
    });
    // Each trace is individually 1 m2; jointly 2 m2 exceeds this cell's 1.5 m2.
    surfaces[0].cell=id(0);surfaces[0].material_id=MaterialId(0);surfaces[0].curve=&curves[0];
    surfaces[0].material=binding(&cards[0],&context);
    support::with_cx(|cx|assert!(make(&enclosure,&surfaces).resolve(cx,&old).is_err()));
}
#[test]
fn late_work_budget_and_cancellation_leave_the_original_timeline_unchanged() {
    let cards=[card(0.4,0.4),card(0.6,0.6)];let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context=context();let surfaces=bindings(&cards,&curves,&context);
    let enclosure=enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
    let model=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
    let old=initial(&curves,[150.0,-50.0],0.0);let saved=old.cells().clone();
    support::with_cx(|cx| {
        let (_,report)=model.step_candidate(cx,&old,0.01,controls(),&[0.0;2]).unwrap();
        assert!(matches!(model.step_candidate(cx,&old,0.01,
            RadiationStepSolve {maximum_evaluations:report.evaluations-1,..controls()},&[0.0;2]),
            Err(RadiationCellError::EvaluationBudget)));
    });
    support::with_cancelled_cx(|cx|assert!(matches!(model.step_candidate(cx,&old,0.01,controls(),&[0.0;2]),
        Err(RadiationCellError::Cancelled))));
    assert_eq!(old.cells(),&saved);assert_eq!(old.time_s(),0.0);
}

/// The original curve with zero conductivity: the cells exchange heat only
/// through the enclosure.
struct Insulated<'a>(&'a EquilibriumEnthalpyPhaseCurve);

impl ConductivityLaw for Insulated<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        ThermalLaw::temperature_k(self.0, h)?;
        Ok(0.0)
    }
}

impl ImplicitConductionLaw for Insulated<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        let knots = self.0.knots();
        [
            knots[0].specific_enthalpy_j_kg,
            knots[knots.len() - 1].specific_enthalpy_j_kg,
        ]
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: ThermalLaw::temperature_k(self.0, h)?,
            temperature_slope_kg_k_j: self
                .0
                .temperature_derivative_at_specific_enthalpy(h)
                .map_err(|_| DepositionError::MaterialDomain)?,
            conductivity_w_m_k: 0.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// One synthetic emissivity claim, linear in liquid fraction and valid from
/// `low_k` to 500 K, in the shape of `common::card`'s claim.
fn emissivity_claim(low_k: f64, e0: f64, e1: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: "melt".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.0, e0), (1.0, e1)],
            dims: EMISSIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with("T", low_k, 500.0)
            .with_quantity(
                "pressure",
                fs_qty::QuantitySpec::dimensional(Pressure::DIMS),
                90000.0,
                110000.0,
            )
            .with("lot", 7.0, 7.0)
            .with("melt", 0.0, 1.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic phase emissivity".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

/// A card whose second, observation-backed emissivity claim enters its
/// domain at 420 K.
fn switching_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(emissivity_claim(300.0, 0.3, 0.8))
        .unwrap();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic emissivity specimen".into(),
            method: "synthetic selection fixture".into(),
            artifact: hash_bytes(b"synthetic observation, not experimental validation"),
            caveats: "test fixture".into(),
            provenance: Provenance {
                source: "synthetic phase emissivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let mut observed = emissivity_claim(420.0, 0.5, 0.5);
    observed.observations = vec![observation];
    claims.insert_claim(observed).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "chart-defined".into(),
            process: "radiation fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn exchange_model<'a>(
    enclosure: &'a GrayDiffuseEnclosure,
    cells: &'a [RadiatingCell<'a>],
) -> MaterialRadiation<'a> {
    MaterialRadiation {
        enclosure,
        cells,
        power_tolerance_w: 1e-8,
        radiosity_tolerance_w_m2: 1e-8,
    }
}

#[test]
fn the_coupled_exchange_carries_resolve_powers_and_an_exact_enthalpy_tangent() {
    // The emissivity follows the liquid fraction. On the 400 K plateau only
    // the emissivity moves with h; on the solid and liquid ramps only the
    // temperature does. Each state puts the cells on different segments.
    let cards = [card(0.3, 0.8), card(0.7, 0.2)];
    let curves = [curve(&cards[0], 200.0), curve(&cards[1], 200.0)];
    let context = context();
    let surfaces = bindings(&cards, &curves, &context);
    let enclosure = enclosure([
        emissivity(&cards[0], &curves[0], 0.0),
        emissivity(&cards[1], &curves[1], 0.0),
    ]);
    let model = exchange_model(&enclosure, &surfaces);
    for h in [[-50.0, 250.0], [100.0, -50.0], [250.0, 120.0]] {
        let thermal = initial(&curves, h, 0.0);
        let exchange = model.exchange(&thermal, || Ok(())).unwrap();
        // CellId order: id(1) first, so exchange index k is cell 1 - k.
        assert_eq!(exchange.cells(), [id(1), id(0)]);
        let query = |h: [f64; 2]| {
            exchange
                .state_at_specific_enthalpies_and_time(&[h[1], h[0]], 0.0, &mut || Ok(()))
                .unwrap()
        };
        let state = query(h);
        let resolved = support::with_cx(|cx| model.resolve(cx, &thermal).unwrap());
        for (power, cell) in state.powers_w.iter().zip(exchange.cells()) {
            assert_eq!(power.to_bits(), resolved.sources_w()[cell].to_bits());
        }
        let step = 1e-3;
        for column in 0..2 {
            let (mut plus, mut minus) = (h, h);
            plus[1 - column] += step;
            minus[1 - column] -= step;
            let (high, low) = (query(plus), query(minus));
            for row in 0..2 {
                let difference = (high.powers_w[row] - low.powers_w[row]) / (2.0 * step);
                let exact = state.power_enthalpy_derivatives_w_kg_j[2 * row + column];
                assert!(
                    (exact - difference).abs() <= 1e-6 * exact.abs().max(1.0),
                    "h {h:?}, dP{row}/dh{column}: exact {exact:e}, central {difference:e}"
                );
            }
        }
    }
}

#[test]
fn the_implicit_graph_steps_the_coupled_exchange_as_the_radiation_step_does() {
    // A 450 K liquid cell radiates into a 350 K solid one, and both stay on
    // their ramps. The implicit graph, with the enclosure attached as its
    // coupled source, and the enclosure's own implicit step solve the same
    // backward-Euler equations.
    let cards = [card(0.3, 0.8), card(0.7, 0.2)];
    let curves = [curve(&cards[0], 200.0), curve(&cards[1], 200.0)];
    let context = context();
    let surfaces = bindings(&cards, &curves, &context);
    let enclosure = enclosure([
        emissivity(&cards[0], &curves[0], 0.0),
        emissivity(&cards[1], &curves[1], 0.0),
    ]);
    let model = exchange_model(&enclosure, &surfaces);
    let old = initial(&curves, [250.0, -50.0], 0.0);
    let laws = [Insulated(&curves[0]), Insulated(&curves[1])];
    let observer = ConductionObserver::new(
        &old,
        &[
            (MaterialId(0), &laws[0] as &dyn ConductivityLaw),
            (MaterialId(1), &laws[1] as &dyn ConductivityLaw),
        ],
        &[],
        &[],
        || Ok(()),
    )
    .unwrap();
    let exchange = model.exchange(&old, || Ok(())).unwrap();
    let graph = ImplicitConductionGraph::new(
        &old,
        &observer,
        &[
            (MaterialId(0), &laws[0] as &dyn ImplicitConductionLaw),
            (MaterialId(1), &laws[1] as &dyn ImplicitConductionLaw),
        ],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
    .with_coupled_heat_source(&old, exchange.cells(), &exchange, 4, || Ok(()))
    .unwrap();
    let dt = 0.05;
    let implicit = ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-9,
        total_tolerance_j: 1e-8,
        ..Default::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, report) = graph
        .advance_to(&old, dt, implicit, &mut workspace, || Ok(()))
        .unwrap();
    let (expected, _) = support::with_cx(|cx| {
        model
            .step_candidate(cx, &old, dt, controls(), &[0.0; 2])
            .unwrap()
    });
    for i in 0..2 {
        let [a, b] = [&next, &expected].map(|state| state.cells()[&id(i)].specific_h_j_kg);
        assert!(
            (a - b).abs() <= 1e-6,
            "cell {i}: graph {a}, radiation step {b}"
        );
    }
    assert!(next.cells()[&id(0)].specific_h_j_kg < 240.0);
    near(next.energy_j(), old.energy_j());
    // The enclosure only moves heat, hot to cold: its production is positive.
    assert!(report.entropy.passive_exchange_generation_j_k > 0.0);
}

#[test]
fn the_coupled_exchange_refuses_a_switching_claim_and_keeps_the_checkpoint_error() {
    // Cell 0's card gains a second, observation-backed emissivity claim at
    // 420 K. The attachment pins the claim the selection picks at 400 K. At
    // 419 K the exchange equals resolve under either selection. At 421 K the
    // single-claim selection refuses the tie and the preferring one moves to
    // the observed claim; the exchange refuses both.
    let cards = [switching_card(), card(0.6, 0.6)];
    let curves = [curve(&cards[0], 200.0), curve(&cards[1], 200.0)];
    let context = context();
    let enclosure = enclosure([
        emissivity(&cards[0], &curves[0], 0.0),
        emissivity(&cards[1], &curves[1], 0.0),
    ]);
    let below = initial(&curves, [219.0, -50.0], 0.0);
    let above = initial(&curves, [221.0, -50.0], 0.0);
    for prefer in [false, true] {
        let mut surfaces = bindings(&cards, &curves, &context);
        if prefer {
            // `bindings` lists cell 1 first: surfaces[1] binds cell 0.
            surfaces[1].selection = MaterialPropertySelection::PreferObservationBacked;
        }
        let model = exchange_model(&enclosure, &surfaces);
        let exchange = model.exchange(&below, || Ok(())).unwrap();
        let query = |hot: f64| {
            exchange.state_at_specific_enthalpies_and_time(&[-50.0, hot], 0.0, &mut || Ok(()))
        };
        let state = query(219.0).unwrap();
        let resolved = support::with_cx(|cx| model.resolve(cx, &below).unwrap());
        for (power, cell) in state.powers_w.iter().zip(exchange.cells()) {
            assert_eq!(power.to_bits(), resolved.sources_w()[cell].to_bits());
        }
        let moved = support::with_cx(|cx| model.resolve(cx, &above));
        assert_eq!(moved.is_ok(), prefer);
        match query(221.0) {
            Err(DepositionError::InvalidInput(message)) => {
                assert!(message.contains("selection"), "{message}");
            }
            other => panic!(
                "expected the switching refusal, got {:?}",
                other.map(|s| s.powers_w)
            ),
        }
    }
    let surfaces = bindings(&cards, &curves, &context);
    let model = exchange_model(&enclosure, &surfaces);
    let exchange = model.exchange(&below, || Ok(())).unwrap();
    // A checkpoint's own error comes back unchanged, here from inside the
    // radiosity assembly (the fourth poll).
    let mut polls = 0;
    let stopped = exchange.state_at_specific_enthalpies_and_time(&[-50.0, 219.0], 0.0, &mut || {
        polls += 1;
        if polls < 4 {
            Ok(())
        } else {
            Err(DepositionError::Budget("caller budget"))
        }
    });
    assert_eq!(
        stopped.err(),
        Some(DepositionError::Budget("caller budget"))
    );
    assert!(matches!(
        model.exchange(&below, || Err(DepositionError::Cancelled)),
        Err(RadiationCellError::Cancelled)
    ));
    assert!(matches!(
        exchange_model(&enclosure, &surfaces[..1]).exchange(&below, || Ok(())),
        Err(RadiationCellError::Invalid(_))
    ));
    assert!(matches!(
        exchange.state_at_specific_enthalpies_and_time(&[0.0], 0.0, &mut || Ok(())),
        Err(DepositionError::InvalidInput(_))
    ));
}
