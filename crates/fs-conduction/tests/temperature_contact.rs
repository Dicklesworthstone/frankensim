//! Source-backed nonlinear contact: independent eight-node P1 thermal references.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, ScalarField, SolveConfig, ThermalBoundary, ThermalBoundaryBuilder,
    ThermalInterfaces};
use fs_conduction::interface::temperature::{
    ContactIterationConfig, TemperatureContact, TemperatureContacts,
};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel};
use fs_rep_mesh::TetComplex;

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 19, kernel_id: 820, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn mesh() -> ConductionMesh {
    let positions = vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.],
        [0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    ConductionMesh::new(TetComplex::from_tets(8, vec![[0,1,2,3],[4,5,6,7]]), positions).unwrap()
}
fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh).adiabatic_remainder().finish().unwrap()
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("normal_pressure", 1e6).unwrap()
}
fn card(split: bool) -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    for (lo, hi) in if split { vec![(280.,340.),(341.,380.)] } else { vec![(280.,380.)] } {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: fs_conduction::TEMPERATURE_DIMS,
                knots: vec![(280.,0.1),(380.,0.7)], dims },
            validity: ValidityDomain::unconstrained().with("T",lo,hi)
                .with("normal_pressure", 999_999.,1_000_001.),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance {
                source: "synthetic R(T)=0.1+0.006(T-280), not measurements".into(),
                license: "internal-test-use".into(), artifact: None,
            },
        }).unwrap();
    }
    let side = |name: &str| SurfaceSpec { material: MaterialStateId {
        chemistry: name.into(),phase: "solid".into(),process: "fixture".into(),revision: 0,
    }, texture_frame: "fixed reference normal".into() };
    InterfaceSystemCard::assemble(side("A"),side("B"),SystemContext {
        medium: "dry".into(),third_body: None,environment: "vacuum".into(),history: "unaged".into(),
    },claims,vec![]).unwrap()
}
fn contacts<'m, 'c>(cx: &Cx<'_>, mesh: &'m ConductionMesh, bc: &ThermalBoundary,
    card: &'c InterfaceSystemCard) -> TemperatureContacts<'m, 'c> {
    let law = TemperatureContact::new("joint", ThermalInterfaces::coincident_face_pairs(mesh).unwrap(),
        card, &point(300.), "T", SelectionPolicy::SingleClaimOnly).unwrap();
    TemperatureContacts::new(cx,mesh,bc,vec![],vec![law],1).unwrap()
}
fn step_config() -> StepConfig {
    let mut linear = SolveConfig::default().linear;
    linear.tolerance = 1e-12; linear.max_iterations = 2048; linear.restart = 8;
    StepConfig { linear, energy_tolerance_j: 1e-8 }
}
fn coupling() -> ContactIterationConfig {
    ContactIterationConfig { max_iterations: 64,relaxation: 0.5,residual_tolerance_j: 1e-8 }
}
fn nonlinear() -> NonlinearStepConfig {
    NonlinearStepConfig { residual_rtol: 1e-12,residual_atol_j: 1e-11,..Default::default() }
}
fn material(variable: bool) -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(if variable {
        vec![(280.,1.6),(380.,3.6)]
    } else { vec![(280.,2.),(380.,2.)] }).unwrap())
}
fn old() -> Vec<f64> { vec![300.,300.,300.,300.,360.,360.,360.,360.] }
fn source() -> ScalarField { ScalarField::Nodal(vec![0.,0.,0.,0.,30.,30.,30.,30.]) }
fn near(a: f64,b: f64,tolerance: f64) { assert!((a-b).abs()<=tolerance,"{a} versus {b}"); }

#[test]
fn local_source_queries_keep_both_traces_pressure_and_one_pinned_claim() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=boundary(&mesh);let card=card(false);
        let set=contacts(cx,&mesh,&bc,&card);
        let evaluated=set.evaluate(cx,&bc,&old()).unwrap();
        let state=&evaluated.states()[0];
        near(state.side_a_temperature_k.min(state.side_b_temperature_k),300.,1e-12);
        near(state.side_a_temperature_k.max(state.side_b_temperature_k),360.,1e-12);
        near(state.interface_temperature_k,330.,1e-12);
        near(state.resistance.value_m2_k_per_w(),0.4,1e-14);
        assert_eq!(state.resistance.card_identity(),card.content_hash());
        assert!(state.resistance.receipt().query_point.contains(&("normal_pressure".into(),1e6)));
        let warmer=set.evaluate(cx,&bc,&[350.;8]).unwrap();
        near(warmer.states()[0].resistance.value_m2_k_per_w(),0.52,1e-14);
        assert_eq!(state.resistance.receipt().selected,warmer.states()[0].resistance.receipt().selected);
        assert_ne!(state.resistance.receipt().content_hash(),warmer.states()[0].resistance.receipt().content_hash());
        assert_eq!(evaluated.states(),set.evaluate(cx,&bc,&old()).unwrap().states());
    });
}

#[test]
fn source_support_and_quantity_refusals_never_choose_another_claim() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=boundary(&mesh);let card=card(true);
        let set=contacts(cx,&mesh,&bc,&card);
        assert!(set.evaluate(cx,&bc,&[350.;8]).is_err(),"original pin ends at340 K");
        // A newly selected source can answer 350 K, proving that the refusal
        // above is a source-continuity gate, not absence of another candidate.
        assert!(fs_conduction::InterfaceResistance::from_card("joint",&card,&point(350.),
            SelectionPolicy::SingleClaimOnly).is_ok());
        let pairs=ThermalInterfaces::coincident_face_pairs(&mesh).unwrap();
        let wrong=point(300.).with_quantity("T",fs_qty::QuantitySpec::dimensional(
            fs_qty::Dims([1,0,0,0,0,0])),300.).unwrap();
        assert!(TemperatureContact::new("joint",pairs.clone(),&card,&wrong,"T",
            SelectionPolicy::SingleClaimOnly).is_err());
        let wrong=point(300.).with("normal_pressure",2e6).unwrap();
        assert!(TemperatureContact::new("joint",pairs,&card,&wrong,"T",
            SelectionPolicy::SingleClaimOnly).is_err());
        assert!(set.evaluate(cx,&bc,&[300.;7]).is_err());
        assert!(set.evaluate(cx,&bc,&[f64::NAN;8]).is_err());
        assert!(set.evaluate(cx,&bc,&[0.;8]).is_err());
        assert!(set.evaluate(cx,&bc,&old()).is_ok(),"a failed query never poisons the law");
    });
}

// Independent P1 roots: diagonal capacities [10/24;4],[20/24;4],
// source loads [0;4],[30/24;4], K_e=grad(N) K(Tmean) grad(N)^T/6,
// and contact blocks +/- (1+delta_ij)/(24 R(mean of the six trace nodes)).
// Newton roots include dR/dT, unlike the production partitioned iteration.
const EXPECTED: [[[f64;8];3];2] = [
    [
        [306.6716930999188,307.1795003067839,307.1795003067839,300.7148242607056,
         356.76960697200474,356.6141111627513,356.6141111627513,360.0294117153965],
        [315.9320470201926,317.8966953274104,317.8966953274104,304.04359173934336,
         352.350202397157,351.6092028289733,351.6092028289733,359.546877237718],
        [323.48193256858633,326.66021099925763,326.66021099925763,309.5974034048414,
         349.09982742515854,347.63479565995596,347.63479565995596,358.4307022689581],
    ],[
        [306.655302613181,307.1775761701527,307.1775761701527,300.74749895490237,
         356.8408390295034,356.62369203650985,356.62369203650985,359.9327999432825],
        [315.8429951814067,317.8536064359657,317.8536064359657,304.3997743808933,
         352.59094341183174,351.6848607659552,351.6848607659552,359.0643438391421],
        [323.3907563852738,326.51753510115856,326.51753510115856,310.6198087064496,
         349.4222794244419,347.823168417243,347.823168417243,357.40856609405176],
    ],
];

#[test]
fn nonlinear_contact_and_bulk_conductivity_match_spatial_roots_and_window_energy() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=boundary(&mesh);let card=card(false);
        let set=contacts(cx,&mesh,&bc,&card);
        let capacities=[VolumetricHeatCapacity::declared(10.).unwrap(),
            VolumetricHeatCapacity::declared(20.).unwrap()];
        let engine=BackwardEuler::per_element(cx,&mesh,&capacities).unwrap();
        let source=source();
        for (case,variable) in [false,true].into_iter().enumerate() {
            let material=material(variable);
            let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
                element_materials:None,source:&source};
            let initial=old();let mut history=initial.clone();let mut stored=0.;
            for (step,dt) in [0.15,0.35,0.5].into_iter().enumerate() {
                let before=history.clone();
                let next=set.advance_capacity(cx,&engine,problem,&history,dt,step_config(),
                    Some(nonlinear()),coupling()).unwrap();
                assert_eq!(history,before,"trials cannot advance physical history");
                assert!(next.iterations>1 && next.iterations<=coupling().max_iterations);
                assert!(next.krylov_iterations>=next.conduction.krylov_iterations);
                assert!(next.physical.residual_norm_j<=coupling().residual_tolerance_j);
                near(next.physical.stored_energy_change_j,5.*dt,1e-8);
                near(next.physical.energy_residual_j,0.,1e-8);
                for (&actual,&expected) in next.conduction.temperature.iter().zip(&EXPECTED[case][step]) {
                    near(actual,expected,2e-7);
                }
                let face=&next.physical_contacts.states()[0];
                near(face.resistance.value_m2_k_per_w(),0.1+0.006*(face.interface_temperature_k-280.),1e-14);
                stored+=next.physical.stored_energy_change_j;
                history=next.conduction.temperature;
            }
            near(stored,5.,3e-8);
            let independent_energy: f64=history.iter().zip(initial).enumerate()
                .map(|(i,(&t,old))|if i<4 {10.*(t-old)/24.}else{20.*(t-old)/24.}).sum();
            near(independent_energy,5.,3e-8);
        }
    });
}

#[test]
fn total_energy_cannot_admit_a_frozen_contact_endpoint_and_refusals_leave_history_intact() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=boundary(&mesh);let card=card(false);
        let set=contacts(cx,&mesh,&bc,&card);let source=source();let material=material(false);
        let engine=BackwardEuler::per_element(cx,&mesh,&[
            VolumetricHeatCapacity::declared(10.).unwrap(),VolumetricHeatCapacity::declared(20.).unwrap(),
        ]).unwrap();
        let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let history=old();let before=history.clone();let snapshot=set.evaluate(cx,&bc,&history).unwrap();
        let frozen=engine.advance_nonlinear_prescribed(cx,problem,Some(snapshot.frozen_interfaces()),
            &history,0.15,step_config(),nonlinear()).unwrap().step;
        let actual=set.evaluate(cx,&bc,&frozen.temperature).unwrap();
        let audit=engine.evaluate_prescribed_endpoint(cx,problem,Some(actual.frozen_interfaces()),
            &history,&frozen.temperature,0.15,step_config(),None).unwrap();
        assert!(audit.energy_residual_j.abs()<1e-8);
        assert!(audit.residual_norm_j>0.2,"opposing contact errors are locally large");
        let limited=ContactIterationConfig {max_iterations:1,..coupling()};
        assert!(set.advance_capacity(cx,&engine,problem,&history,0.15,step_config(),
            Some(nonlinear()),limited).is_err());
        assert_eq!(history,before);
        let first=set.advance_capacity(cx,&engine,problem,&history,0.15,step_config(),Some(nonlinear()),coupling()).unwrap();
        let second=set.advance_capacity(cx,&engine,problem,&history,0.15,step_config(),Some(nonlinear()),coupling()).unwrap();
        assert_eq!(first.conduction.temperature,second.conduction.temperature);
        assert_eq!(first.physical_contacts.states(),second.physical_contacts.states());
        assert_eq!(first.iterations,second.iterations);
    });
}

#[test]
fn cancelled_contact_evaluation_and_step_publish_no_candidate() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=boundary(&mesh);let card=card(false);
        let set=contacts(cx,&mesh,&bc,&card);let source=source();let material=material(false);
        let engine=BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(10.).unwrap()).unwrap();
        let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let history=old();let before=history.clone();let gate=CancelGate::new_clock_free();gate.request();
        context(&gate,|stopped| {
            assert!(matches!(set.evaluate(stopped,&bc,&history),Err(ConductionError::Cancelled {..})));
            assert!(matches!(set.advance_capacity(stopped,&engine,problem,&history,0.15,
                step_config(),Some(nonlinear()),coupling()),Err(ConductionError::Cancelled {..})));
        });
        assert_eq!(history,before);
    });
}
