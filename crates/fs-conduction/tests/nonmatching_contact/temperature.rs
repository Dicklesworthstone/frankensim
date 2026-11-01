//! Physical source laws on the retained intersection quadrature.
use super::*;
use fs_conduction::interface::temperature::{
    ContactIterationConfig, TemperatureContact, TemperatureContactCoordinates, TemperatureContacts,
};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};

fn source_value(a: f64, b: f64) -> f64 {
    0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0
}

pub(super) fn source_point() -> QueryPoint {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    let kelvin = fs_qty::QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    QueryPoint::new().with_quantity("TA", kelvin, 320.0).unwrap()
        .with_quantity("TB", kelvin, 320.0).unwrap()
        .with_quantity("normal_pressure",
            fs_qty::QuantitySpec::dimensional(fs_qty::Pressure::DIMS), 1e6).unwrap()
}

pub(super) fn source_card(scale: f64, reverse: bool) -> InterfaceSystemCard {
    let dims = AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Surface2D {
            axes: ["TA".into(), "TB".into()],
            axis_dims: [fs_qty::Temperature::DIMS; 2],
            samples: [[256.0,256.0], [384.0,256.0], [256.0,384.0], [384.0,384.0]]
                .into_iter().map(|p| (p, scale * source_value(p[0], p[1]))).collect(),
            triangles: vec![[0,1,2], [1,3,2]], dims,
        },
        validity: ValidityDomain::unconstrained().with("TA",256.0,384.0)
            .with("TB",256.0,384.0).with("normal_pressure",1e6,1e6),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside, observations: vec![],
        provenance: Provenance {
            source: "manufactured ordered surface law, not measured contact data".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId { chemistry:name.into(), phase:"solid".into(),
            process:"fixture".into(), revision:0 },
        texture_frame: "fixed reference normal".into(),
    };
    let names = if reverse { ["B","A"] } else { ["A","B"] };
    InterfaceSystemCard::assemble(side(names[0]), side(names[1]), SystemContext {
        medium:"dry".into(), third_body:None, environment:"test".into(), history:"unaged".into(),
    }, claims, vec![]).unwrap()
}

fn law<'a>(card: &'a InterfaceSystemCard, a: Vec<usize>, b: Vec<usize>,
    reverse: bool) -> TemperatureContact<'a> {
    let (a,b,axes) = if reverse { (b,a,["TB","TA"]) } else { (a,b,["TA","TB"]) };
    TemperatureContact::new_nonmatching("bond", a, b, card, &source_point(),
        TemperatureContactCoordinates::TwoSided(axes.map(String::from)),
        SelectionPolicy::SingleClaimOnly, options()).unwrap()
}

pub(super) fn bind<'m,'c>(cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh, bc: &ThermalBoundary,
    card: &'c InterfaceSystemCard, a: Vec<usize>, b: Vec<usize>, reverse: bool,
) -> TemperatureContacts<'m,'c> {
    TemperatureContacts::with_nonmatching(cx,mesh,bc,vec![],vec![],
        vec![law(card,a,b,reverse)],10_000,10_000).unwrap()
}

fn independent_load(mesh: &ConductionMesh, a: &[usize], b: &[usize], flux: f64) -> Vec<f64> {
    let mut load = vec![0.0;mesh.vertex_count()];
    // Constant physical heat flux integrates against each original trace's
    // P1 basis to area/3, independently of the opposite tessellation.
    for (faces,sign) in [(a,1.0),(b,-1.0)] {
        for &slot in faces {
            let face = &mesh.boundary()[slot];
            for vertex in face.vertices {
                load[vertex as usize] += sign * flux * face.area / 3.0;
            }
        }
    }
    load
}

fn affine_field(mesh: &ConductionMesh, on_a: &[bool]) -> Vec<f64> {
    mesh.positions().iter().zip(on_a).map(|(p,&a)| {
        let y = if p[1] > 1.5 { p[1]-2.0 } else { p[1] };
        if a {344.0+19.0*y} else {312.0+15.0*y}
    }).collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a,b)|a*b).sum()
}

#[test]
fn g1_pointwise_nonmatching_sources_recover_exact_heat_on_mixed_tessellations() {
    with_cx(|cx| {
        for right in [[[1,2,3],[1,3,2]], [[1,1,1],[1,2,3]]] {
            let (mesh,on_a) = parallel_contact_mesh(right);
            let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
            let (a,b) = parallel_sides(&mesh);
            let card = source_card(1.0,false);
            let set = bind(cx,&mesh,&bc,&card,a.clone(),b.clone(),false);
            let t = affine_field(&mesh,&on_a);
            // R=.5+y/16 and TA-TB=32+4y: q=64 W/m2 exactly.
            let linear = set.linearize_temperature(cx,&bc,&t,10_000).unwrap();
            let expected = independent_load(&mesh,&a,&b,64.0);
            for (&actual,&expected) in linear.nodal_outflow_w().iter().zip(&expected) {
                near(actual,expected,2e-10);
            }
            near(linear.nodal_outflow_w().iter().sum(),0.0,2e-10);
            let physical = linear.physical_contacts();
            assert!(physical.states().is_empty());
            let states = physical.quadrature_states();
            assert!(!states.is_empty());
            assert!(states.iter().all(|s| s.weight_m2>0.0));
            near(states.iter().map(|s|s.weight_m2).sum(),2.0,2e-12);
            let first = states[0].resistance.receipt().selected;
            for state in states {
                let ta = state.side_a_temperature_k;
                let tb = state.side_b_temperature_k;
                near(ta,344.0+19.0*(tb-312.0)/15.0,2e-12);
                near(state.resistance.value_m2_k_per_w(),source_value(ta,tb),2e-14);
                assert_eq!(state.resistance.card_identity(),card.content_hash());
                assert_eq!(state.resistance.receipt().selected,first);
                for coordinate in [("TA",ta),("TB",tb),("normal_pressure",1e6)] {
                    assert!(state.resistance.receipt().query_point
                        .contains(&(coordinate.0.into(),coordinate.1)));
                }
            }
            let flux = physical.frozen_interfaces().fluxes(&t).unwrap();
            assert_eq!(flux.len(),1);
            near(flux[0].heat_rate_a_to_b_w,128.0,2e-10);
            near(flux[0].area_m2,2.0,2e-12);
        }
    });
}

#[test]
fn g3_nonmatching_source_jacobian_and_transpose_match_fresh_physical_queries() {
    with_cx(|cx| {
        let (mesh,split) = mesh([1,2,3],0.0);
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let (a,b) = sides(&mesh,split,0.0);
        let card = source_card(1.0,false);
        let set = bind(cx,&mesh,&bc,&card,a.clone(),b.clone(),false);
        let on_a:Vec<_> = (0..mesh.vertex_count()).map(|i|i<split).collect();
        let t = affine_field(&mesh,&on_a);
        let direction:Vec<_> = (0..t.len()).map(|i|(i%7) as f64/5.0-0.4).collect();
        let weights:Vec<_> = (0..t.len()).map(|i|(i%11) as f64/7.0-0.7).collect();
        let linear = set.linearize_temperature(cx,&bc,&t,10_000).unwrap();
        let tangent = linear.apply(cx,&direction).unwrap();
        let dual = linear.apply_transpose(cx,&weights).unwrap();
        near(dot(&weights,&tangent),dot(&direction,&dual),2e-12);
        near(tangent.iter().sum(),0.0,2e-12);
        let h = 1e-3;
        let plus:Vec<_> = t.iter().zip(&direction).map(|(t,d)|t+h*d).collect();
        let minus:Vec<_> = t.iter().zip(&direction).map(|(t,d)|t-h*d).collect();
        let plus = set.linearize_temperature(cx,&bc,&plus,10_000).unwrap();
        let minus = set.linearize_temperature(cx,&bc,&minus,10_000).unwrap();
        for ((&actual,&plus),&minus) in tangent.iter()
            .zip(plus.nodal_outflow_w()).zip(minus.nodal_outflow_w()) {
            near(actual,(plus-minus)/(2.0*h),2e-8);
        }
        let reverse_card = source_card(1.0,true);
        let reverse = bind(cx,&mesh,&bc,&reverse_card,a,b,true)
            .linearize_temperature(cx,&bc,&t,10_000).unwrap();
        for (&forward,&backward) in linear.nodal_outflow_w().iter().zip(reverse.nodal_outflow_w()) {
            near(forward,backward,2e-10);
        }
        for (&forward,&backward) in tangent.iter().zip(reverse.apply(cx,&direction).unwrap().iter()) {
            near(forward,backward,2e-12);
        }
        near(reverse.physical_contacts().frozen_interfaces().fluxes(&t).unwrap()[0]
            .heat_rate_a_to_b_w,-64.0,2e-10);
    });
}

fn step_config() -> StepConfig {
    let mut linear = SolveConfig::default().linear;
    linear.tolerance=1e-12;linear.max_iterations=10_000;
    StepConfig {linear,energy_tolerance_j:1e-8}
}
fn coupling() -> ContactIterationConfig {
    ContactIterationConfig {max_iterations:64,relaxation:1.0,residual_tolerance_j:1e-9}
}

#[test]
fn g1_g3_nonmatching_contact_advances_manufactured_energy_and_differentiates_its_source() {
    with_cx(|cx| {
        let (mesh,split) = mesh([1,2,3],0.0);
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let (a,b) = sides(&mesh,split,0.0);
        let card = source_card(1.0,false);
        let set = bind(cx,&mesh,&bc,&card,a.clone(),b.clone(),false);
        let target:Vec<_> = (0..mesh.vertex_count()).map(|i|if i<split {344.0}else{312.0}).collect();
        let load = independent_load(&mesh,&a,&b,64.0);
        let mut mass = vec![0.0;target.len()];
        for (e,tet) in mesh.complex().tets.iter().enumerate() {
            for &vertex in tet { mass[vertex as usize] += 10.0*mesh.element_volume(e)/4.0; }
        }
        let dt = 0.01;
        // At the chosen endpoint both bulk gradients vanish. This old state
        // satisfies M(T-old)+dt*contact_load=0 independently of the solver.
        let old:Vec<_> = target.iter().zip(&load).zip(&mass).map(|((t,q),m)|t+dt*q/m).collect();
        let before = old.clone();
        let engine = BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(10.0).unwrap()).unwrap();
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let next = set.advance_capacity(cx,&engine,problem,&old,dt,step_config(),None,coupling()).unwrap();
        for (&actual,&expected) in next.conduction.temperature.iter().zip(&target) {
            near(actual,expected,3e-8);
        }
        near(next.physical.stored_energy_change_j,0.0,1e-8);
        near(next.physical.energy_residual_j,0.0,1e-8);
        near(mass.iter().zip(&next.conduction.temperature).zip(&old)
            .map(|((m,t),old)|m*(t-old)).sum(),0.0,1e-8);
        assert!(next.physical.residual_norm_j<=coupling().residual_tolerance_j);
        assert!(next.iterations>1);
        assert_eq!(old,before);
        let vertices:std::collections::BTreeSet<_> = a.iter()
            .flat_map(|&slot|mesh.boundary()[slot].vertices).collect();
        let mut weights = vec![0.0;target.len()];
        for &vertex in &vertices {weights[vertex as usize]=1.0/vertices.len() as f64;}
        let linear = engine.linearize_temperature_contact_endpoint(cx,problem,&set,&old,
            &next.conduction.temperature,dt,step_config(),coupling().residual_tolerance_j,&[],10_000).unwrap();
        let gradient = linear.pullback_with_resistance_controls(cx,&weights,&[],&[],10_000).unwrap();
        assert_eq!(gradient.resistance_multipliers.len(),1);
        assert_eq!(gradient.resistance_multipliers[0].interface,"bond");
        let objective = |scale| {
            let card = source_card(scale,false);
            let contacts = bind(cx,&mesh,&bc,&card,a.clone(),b.clone(),false);
            let solution = contacts.advance_capacity(cx,&engine,problem,&old,dt,
                step_config(),None,coupling()).unwrap();
            dot(&weights,&solution.conduction.temperature)
        };
        let h = 5e-4;
        let difference = (objective(1.0+h)-objective(1.0-h))/(2.0*h);
        near(gradient.resistance_multipliers[0].derivative,difference,2e-5);
        assert_eq!(old,before);
    });
}

#[test]
fn g0_nonmatching_source_support_and_quadrature_limits_refuse_transactionally() {
    with_cx(|cx| {
        let (mesh,split) = mesh([1,2,3],0.0);
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let (a,b) = sides(&mesh,split,0.0);
        let card = source_card(1.0,false);
        assert!(TemperatureContacts::with_nonmatching(cx,&mesh,&bc,vec![],vec![],
            vec![law(&card,a.clone(),b.clone(),false)],10_000,1).is_err());
        let mut missing = b.clone();missing.pop();
        assert!(TemperatureContacts::with_nonmatching(cx,&mesh,&bc,vec![],vec![],
            vec![law(&card,a.clone(),missing,false)],10_000,10_000).is_err());
        let set = bind(cx,&mesh,&bc,&card,a,b,false);
        let good = vec![320.0;mesh.vertex_count()];
        let before = set.evaluate(cx,&bc,&good).unwrap();
        assert!(set.evaluate(cx,&bc,&vec![400.0;mesh.vertex_count()]).is_err());
        assert!(set.linearize_temperature(cx,&bc,&good,0).is_err());
        assert_eq!(before.quadrature_states(),set.evaluate(cx,&bc,&good).unwrap().quadrature_states());
    });
}
