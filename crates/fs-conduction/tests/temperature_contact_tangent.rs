//! Full contact-law derivatives, independently checked against fresh source queries.
#[path = "temperature_contact_tangent/capacity.rs"]
mod capacity;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionError, ConductionMesh, ThermalBoundary, ThermalBoundaryBuilder,
    ThermalInterfaces, InterfaceFacePair};
use fs_conduction::interface::temperature::{TemperatureContact, TemperatureContacts};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel};
use fs_rep_mesh::TetComplex;

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 20, kernel_id: 821, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(8, vec![[0,1,2,3],[4,5,6,7]]),
        vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.],
            [0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).unwrap()
}
fn card_with_knots(knots: Vec<(f64, f64)>) -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let lo = knots[0].0;
    let hi = knots[knots.len()-1].0;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: fs_conduction::TEMPERATURE_DIMS,
            knots, dims },
        validity: ValidityDomain::unconstrained().with("T",lo,hi).with("normal_pressure",1e6,1e6),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "synthetic contact tangent fixture, not measurements".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let side = |name: &str| SurfaceSpec { material: MaterialStateId {
        chemistry: name.into(), phase: "solid".into(), process: "fixture".into(), revision: 0,
    }, texture_frame: "reference normal".into() };
    InterfaceSystemCard::assemble(side("A"),side("B"),SystemContext {
        medium: "dry".into(), third_body: None, environment: "vacuum".into(), history: "unaged".into(),
    },claims,vec![]).unwrap()
}
fn card(scale: f64) -> InterfaceSystemCard {
    card_with_knots(vec![(280.,0.1*scale),(380.,0.7*scale)])
}
fn contacts<'m, 'c>(cx: &Cx<'_>, mesh: &'m ConductionMesh, bc: &ThermalBoundary,
    card: &'c InterfaceSystemCard, reverse: bool) -> TemperatureContacts<'m, 'c> {
    let mut pairs = ThermalInterfaces::coincident_face_pairs(mesh).unwrap();
    if reverse {
        pairs = pairs.into_iter().map(|p| InterfaceFacePair {side_a:p.side_b,side_b:p.side_a}).collect();
    }
    let point=QueryPoint::new().with("T",300.).unwrap().with("normal_pressure",1e6).unwrap();
    let law=TemperatureContact::new("joint",pairs,card,&point,"T",SelectionPolicy::SingleClaimOnly).unwrap();
    TemperatureContacts::new(cx,mesh,bc,vec![],vec![law],1).unwrap()
}
fn near(a: f64,b: f64,tol: f64) { assert!((a-b).abs()<=tol,"{a} versus {b}; allowance {tol}"); }
fn dot(a: &[f64],b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }

#[test]
fn g3_contact_tangent_matches_fresh_queries_and_its_true_transpose() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card=card(1.);let set=contacts(cx,&mesh,&bc,&card,false);
        let temperature=[301.,303.,307.,309.,353.,349.,345.,341.];
        let direction=[0.2,-0.1,0.7,0.5,-0.3,0.9,0.4,-0.6];
        let weights=[0.9,-0.4,0.2,0.3,0.8,0.1,-0.5,0.7];
        let linear=set.linearize_temperature(cx,&bc,&temperature,1).unwrap();
        let tangent=linear.apply(cx,&direction).unwrap();
        let adjoint=linear.apply_transpose(cx,&weights).unwrap();
        near(dot(&weights,&tangent),dot(&direction,&adjoint),1e-13);
        near(tangent.iter().sum(),0.,1e-13);
        assert_eq!(linear.temperature(),temperature);
        let h=1e-3;
        let plus:Vec<_>=temperature.iter().zip(direction).map(|(t,d)|t+h*d).collect();
        let minus:Vec<_>=temperature.iter().zip(direction).map(|(t,d)|t-h*d).collect();
        let a=set.linearize_temperature(cx,&bc,&plus,1).unwrap();
        let b=set.linearize_temperature(cx,&bc,&minus,1).unwrap();
        for i in 0..8 {
            near(tangent[i],(a.nodal_outflow_w()[i]-b.nodal_outflow_w()[i])/(2.*h),2e-8);
        }
        let reversed=contacts(cx,&mesh,&bc,&card,true).linearize_temperature(cx,&bc,&temperature,1).unwrap();
        for (&a,&b) in tangent.iter().zip(reversed.apply(cx,&direction).unwrap().iter()) { near(a,b,1e-13); }
        let repeat=set.linearize_temperature(cx,&bc,&temperature,1).unwrap();
        assert_eq!(linear.apply(cx,&direction).unwrap(),repeat.apply(cx,&direction).unwrap());
    });
}

#[test]
fn g1_zero_mean_jump_keeps_nonuniform_contact_feedback() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card=card(1.);let set=contacts(cx,&mesh,&bc,&card,false);
        let t=[331.,329.,330.,330.,330.,330.,330.,330.];
        let linear=set.linearize_temperature(cx,&bc,&t,1).unwrap();
        let derivative=linear.apply(cx,&[1.;8]).unwrap();
        // R=.4, R'=.006; B*jump=[1,-1,0]/24. A uniform shift
        // changes R but not either jump. A mean-jump-only model loses this.
        let value=-0.006/(24.*0.4*0.4);
        for (actual,expected) in derivative.iter().zip([value,-value,0.,0.,-value,value,0.,0.]) {
            near(*actual,expected,1e-13);
        }
        assert!(derivative.iter().any(|v|v.abs()>1e-4));
        near(linear.nodal_outflow_w()[..4].iter().sum(),0.,1e-12);
        near(derivative.iter().sum(),0.,1e-13);
    });
}

#[test]
fn g0_contact_derivatives_refuse_source_corners_domains_shapes_and_budget() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card=card_with_knots(vec![(256.,0.25),(320.,0.5),(384.,0.75)]);
        let set=contacts(cx,&mesh,&bc,&card,false);
        assert!(set.linearize_temperature(cx,&bc,&[320.;8],1).is_ok(),"redundant equal-slope knot");
        for t in [256.,384.] { assert!(set.linearize_temperature(cx,&bc,&[t;8],1).is_err()); }
        let corner=card_with_knots(vec![(256.,0.25),(320.,0.5),(384.,1.)]);
        let kink=contacts(cx,&mesh,&bc,&corner,false);
        assert!(kink.evaluate(cx,&bc,&[320.;8]).is_ok(),"valid primal despite nondifferentiability");
        assert!(kink.linearize_temperature(cx,&bc,&[320.;8],1).is_err());
        assert!(set.linearize_temperature(cx,&bc,&[320.;8],0).is_err());
        assert!(set.linearize_temperature(cx,&bc,&[320.;7],1).is_err());
        let linear=set.linearize_temperature(cx,&bc,&[321.;8],1).unwrap();
        assert!(linear.apply(cx,&[0.;7]).is_err());
        assert!(linear.apply_transpose(cx,&[f64::NAN;8]).is_err());
        let gate=CancelGate::new_clock_free();gate.request();
        context(&gate,|stopped| {
            assert!(matches!(set.linearize_temperature(stopped,&bc,&[320.;8],1),Err(ConductionError::Cancelled {..})));
            assert!(matches!(linear.apply(stopped,&[0.;8]),Err(ConductionError::Cancelled {..})));
        });
    });
}
