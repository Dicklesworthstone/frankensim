//! Matching-contact FEM exercises the inverse fallback through the real Robin
//! feedback binding; a dense elimination oracle does not call production CG.
mod support;
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalAnalyzer, RobinFeedbackAnalysisConfig};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InterfaceFacePair,
    InterfaceResistance, InterfaceSurface, LinearConfig, ScalarField, ThermalBc,
    ThermalBoundary, ThermalBoundaryBuilder, ThermalInterfaces,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec, SystemContext,
    UncertaintyModel};
use fs_rep_mesh::TetComplex;
use fs_solver::goal::GoalResidualLimits;
use fs_solver::goal::feedback::{FeedbackInverseMethod, FeedbackResidualLimits};
use fs_sparse::Csr;
use support::{with_cx, with_cancelled_cx};

struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel,
    source: ScalarField, interfaces: ThermalInterfaces,
}
impl Fixture {
    fn new() -> Self {
        let (mut positions, mut tets) = (Vec::new(), Vec::new());
        for side in 0..2 {
            let (complex, points) = unit_cube(1);
            let offset = positions.len() as u32;
            tets.extend(complex.tets.into_iter().map(|tet| tet.map(|v| v+offset)));
            positions.extend(points.into_iter().map(|[x,y,z]| [x+side as f64,y,z]));
        }
        let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh).region("air",
            |f| on_box_face(f.centroid[0],0.) || on_box_face(f.centroid[0],2.),
            ThermalBc::robin(2.,300.).unwrap()).unwrap().adiabatic_remainder().finish().unwrap();
        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(RP,RD),
            value: PropertyValue::Scalar { value: 0.001, dims: RD }, validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(), provenance: Provenance { source: "synthetic strong-contact test".into(),
                license: "internal-test-use".into(), artifact: None } }).unwrap();
        let surface = |name: &str| SurfaceSpec { material: MaterialStateId { chemistry: name.into(),
            phase: "solid".into(), process: "fixture".into(), revision: 0 }, texture_frame: "declared".into() };
        let card = InterfaceSystemCard::assemble(surface("left"), surface("right"), SystemContext {
            medium: "declared".into(), third_body: None, environment: "declared".into(), history: "declared".into(),
        }, claims, Vec::new()).unwrap();
        let resistance = InterfaceResistance::from_card("bond", &card, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly).unwrap();
        let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap().into_iter().map(|pair| {
            if mesh.boundary()[pair.side_a].outward_normal[0] > 0. { pair }
            else { InterfaceFacePair { side_a: pair.side_b, side_b: pair.side_a } }
        }).collect();
        let interfaces = ThermalInterfaces::new(&mesh,&boundary,
            vec![InterfaceSurface::new("bond",pairs,resistance).unwrap()]).unwrap();
        Self { mesh, boundary, interfaces, material: ConductivityModel::isotropic_declared(10.).unwrap(), source: ScalarField::Uniform(0.) }
    }
    fn base(&self, cx: &fs_exec::Cx<'_>) -> LinearGoalAnalyzer<'_> {
        LinearGoalAnalyzer::new_for_maximum(cx, ConductionProblem {
            mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            source: &self.source, element_materials: None,
        }, Some(&self.interfaces), LinearConfig { tolerance: 1e-10, max_iterations: 200, restart: 32 },
            &vec![300.; self.mesh.vertex_count()], LinearGoalAnalysisConfig {
                residual_limits: limits().residual.solid, max_stability_iterations: 5000,
            }).unwrap()
    }
}
fn limits() -> RobinFeedbackAnalysisConfig {
    RobinFeedbackAnalysisConfig {
        residual: FeedbackResidualLimits { solid: GoalResidualLimits { max_rows: 256, max_nonzeros: 100_000 },
            max_ports: 4, max_transfer_nonzeros: 1000, max_response_entries: 1000,
            max_verification_entries: 100_000 }, max_response_iterations: 200, max_lowering_entries: 100_000,
    }
}
fn dense(a: &Csr, rhs: &[f64], b: &Csr, c: &Csr, d: &[f64]) -> Vec<f64> {
    let n=a.nrows();let mut m=vec![vec![0.;n+1];n];
    for i in 0..n {
        m[i][n]=rhs[i];
        for k in 0..d.len() {m[i][n]+=b.get(i,k)*d[k];}
        for j in 0..n {m[i][j]=a.get(i,j);for k in 0..d.len(){m[i][j]-=b.get(i,k)*c.get(k,j);}}
    }
    for k in 0..n {
        let pivot=(k..n).max_by(|&i,&j|m[i][k].abs().total_cmp(&m[j][k].abs())).unwrap();m.swap(k,pivot);
        let diagonal=m[k][k];assert!(diagonal.abs()>1e-12);
        for j in k..=n {m[k][j]/=diagonal;}
        for i in 0..n {if i!=k {let factor=m[i][k];for j in k..=n{let v=m[k][j];m[i][j]-=factor*v;}}}
    }
    m.iter().map(|row|row[n]).collect()
}

#[test]
fn strong_contact_keeps_its_verified_inverse_when_air_feedback_is_bound() {
    let f=Fixture::new();
    with_cx(|cx| {
        for slope in [0.1,-2.] {
            let base=f.base(cx);
            assert!(base.inverse_columns().is_some(),"fixture must exercise the non-comparison inverse path");
            let saved=base.inverse_columns().unwrap().to_vec();
            let coupled=base.with_robin_feedback(cx,&["air"],&[330.*(1.-slope)],&[slope],limits()).unwrap();
            assert_eq!(coupled.solid_inverse_columns().unwrap(),saved.as_slice());
            assert!(coupled.stability_iterations()<=5000);
            let t=vec![300.;f.mesh.vertex_count()];let vertices:Vec<_>=(0..t.len()).collect();
            let report=coupled.analyze_maximum(cx,&t,&vertices).unwrap();
            let (a,rhs,b,c,d)=coupled.stored_system();let solved=dense(a,rhs,b,c,d);
            let max=solved.iter().copied().fold(f64::NEG_INFINITY,f64::max);
            assert!((max-330.).abs()<1e-7);
            let [lo,hi]=report.interval_k().expect("checked contact inverse reaches the coupled checker");
            assert!(lo<=max && max<=hi,"{report:?}");
            assert!(report.algebraic_half_width_k().unwrap()>=29.99999);
            if slope<0. {assert_eq!(report.coupled().inverse_method(),Some(FeedbackInverseMethod::PortSchurDominance));}
            let mut reversed=vertices.clone();reversed.reverse();
            assert_eq!(report,coupled.analyze_maximum(cx,&t,&reversed).unwrap());
            let before=coupled.response_iterations();
            let close=coupled.analyze_maximum(cx,&vec![329.;t.len()],&vertices).unwrap();
            assert!(close.algebraic_half_width_k().unwrap()<report.algebraic_half_width_k().unwrap());
            assert_eq!(before,coupled.response_iterations());
        }
    });
}

#[test]
fn inverse_verification_is_admitted_before_response_solves() {
    let f=Fixture::new();
    with_cx(|cx| {
        let mut config=limits();config.residual.max_verification_entries=1;
        assert!(f.base(cx).with_robin_feedback(cx,&["air"],&[297.],&[0.1],config).is_err());
        let mut config=limits();config.residual.solid.max_nonzeros=1000;
        // A sparse matrix can fit while the required inverse pass cannot.
        assert!(f.base(cx).with_robin_feedback(cx,&["air"],&[297.],&[0.1],config).is_err());
    });
}

#[test]
fn contact_maximum_selection_and_cancellation_never_publish_partial_results() {
    let f=Fixture::new();
    with_cx(|cx| {
        let coupled=f.base(cx).with_robin_feedback(cx,&["air"],&[297.],&[0.1],limits()).unwrap();
        let t=vec![300.;f.mesh.vertex_count()];let all:Vec<_>=(0..t.len()).collect();
        let first=coupled.analyze_maximum(cx,&t,&all).unwrap();
        for selection in [vec![],vec![0,0],vec![t.len()]] {assert!(coupled.analyze_maximum(cx,&t,&selection).is_err());}
        let mut bad=t.clone();bad[0]=f64::NAN;assert!(coupled.analyze_maximum(cx,&bad,&all).is_err());
        with_cancelled_cx(|cancelled| assert!(matches!(coupled.analyze_maximum(cancelled,&t,&all),
            Err(fs_conduction::ConductionError::Cancelled { .. }))));
        assert_eq!(first,coupled.analyze_maximum(cx,&t,&all).unwrap());
    });
}
