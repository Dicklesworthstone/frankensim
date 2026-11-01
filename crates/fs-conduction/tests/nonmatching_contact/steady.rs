//! Steady contact roots and reservoir accounting, independent of iteration.
use super::*;
use super::temperature::{bind, source_card, source_point};
use fs_conduction::interface::temperature::{TemperatureContact, TemperatureContacts};
use fs_conduction::interface::temperature::steady::SteadyContactConfig;

fn solver_config() -> SolveConfig {
    let mut config = SolveConfig::default();
    config.nonlinearity = Nonlinearity::FixedPoint {relaxation:1.0,max_backtracks:16};
    config.linear.tolerance=1e-12;config.linear.max_iterations=20_000;
    config.stop.residual_rtol=1e-13;config.stop.step_atol=0.0;
    config
}
fn coupling() -> SteadyContactConfig {
    SteadyContactConfig {max_iterations:64,relaxation:1.0,
        residual_tolerance_w:1e-8,energy_tolerance_w:1e-8}
}
fn reservoir_boundary(mesh:&ConductionMesh,robin:bool)->ThermalBoundary {
    if !robin {return boundary(mesh);}
    ThermalBoundaryBuilder::new(mesh)
        .region("cold",|f|on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .region("hot reservoir",|f|on_box_face(f.centroid[0],2.0),ThermalBc::robin(10.0,330.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}
fn contacts<'m,'c>(cx:&fs_exec::Cx<'_>,mesh:&'m ConductionMesh,bc:&ThermalBoundary,
    card:&'c InterfaceSystemCard,split:usize,matching:bool)->TemperatureContacts<'m,'c> {
    let (a,b)=sides(mesh,split,0.0);
    if !matching {return bind(cx,mesh,bc,card,a,b,false);}
    let pairs=ThermalInterfaces::coincident_face_pairs(mesh).unwrap().into_iter().map(|p|
        if a.contains(&p.side_a) {p} else {InterfaceFacePair {side_a:p.side_b,side_b:p.side_a}}
    ).collect();
    let law=TemperatureContact::new_two_sided("bond",pairs,card,&source_point(),
        ["TA","TB"],SelectionPolicy::SingleClaimOnly).unwrap();
    TemperatureContacts::new(cx,mesh,bc,vec![],vec![law],10_000).unwrap()
}
fn scalar_heat(robin:bool)->f64 {
    // Independent uniform one-dimensional series law: each bulk slab is
    // 0.1 K/W and the optional reservoir film contributes another 0.1 K/W.
    let film=if robin {0.1}else{0.0};
    let (mut lo,mut hi)=(0.0,30.0/(0.2+film));
    for _ in 0..80 {
        let heat=0.5*(lo+hi);
        let ta=300.0+0.1*heat;
        let tb=330.0-(0.1+film)*heat;
        let resistance=0.5+(ta-320.0)/1024.0+3.0*(tb-320.0)/1024.0;
        if heat*(0.2+film+resistance)>30.0 {hi=heat;}else{lo=heat;}
    }
    0.5*(lo+hi)
}

#[test]
fn g1_steady_contact_solves_matching_and_nonmatching_slabs_with_physical_reservoirs() {
    with_cx(|cx| {
        for matching in [true,false] {for robin in [false,true] {
            let (mesh,split)=mesh(if matching {[1,1,1]}else{[1,2,3]},0.0);
            let bc=reservoir_boundary(&mesh,robin);
            let card=source_card(1.0,false);
            let set=contacts(cx,&mesh,&bc,&card,split,matching);
            let material=ConductivityModel::isotropic_declared(10.0).unwrap();
            let source=ScalarField::Uniform(0.0);
            let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
                element_materials:None,source:&source};
            let initial=vec![315.0;mesh.vertex_count()];let before=initial.clone();
            let solved=set.solve_steady(cx,problem,&initial,solver_config(),coupling()).unwrap();
            let heat=scalar_heat(robin);
            let film=if robin {0.1}else{0.0};
            for (i,(p,&temperature)) in mesh.positions().iter().zip(&solved.conduction.temperature).enumerate() {
                let expected=if i<split {300.0+0.1*heat*p[0]}
                    else {330.0-heat*(film+0.1*(2.0-p[0]))};
                near(temperature,expected,2e-7);
            }
            assert!(solved.iterations>1);
            assert!(solved.residual_norm_w<=coupling().residual_tolerance_w);
            near(solved.energy.closure_w,0.0,coupling().energy_tolerance_w);
            assert!(solved.free_residual_w.windows(2).all(|p|p[0].0<p[1].0));
            let flux=solved.physical_contacts.frozen_interfaces()
                .fluxes(&solved.conduction.temperature).unwrap();
            near(flux[0].heat_rate_a_to_b_w,-heat,2e-7);
            for (x,expected) in if robin {vec![(0.0,-heat)]}else{vec![(0.0,-heat),(2.0,heat)]} {
                near(solved.dirichlet_reactions_w.iter().filter(|(i,_)|
                    on_box_face(mesh.positions()[*i][0],x)).map(|(_,q)|q).sum(),expected,2e-7);
            }
            if robin {
                assert_eq!(solved.robin_fluxes.len(),1);
                near(solved.energy.robin_out_w,-heat,2e-7);
                near(solved.robin_fluxes[0].heat_rate_w,-heat,2e-7);
                near(solved.energy.dirichlet_in_w,-heat,2e-7);
            }else{assert!(solved.robin_fluxes.is_empty());}
            assert_eq!(initial,before);
        }}
    });
}

#[test]
fn g0_g1_global_heat_cancellation_cannot_accept_a_frozen_steady_contact() {
    with_cx(|cx| {
        let (mesh,split)=mesh([1,2,3],0.0);let bc=boundary(&mesh);
        let card=source_card(1.0,false);
        let set=contacts(cx,&mesh,&bc,&card,split,false);
        let material=ConductivityModel::isotropic_declared(10.0).unwrap();
        let source=ScalarField::Uniform(0.0);
        let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let initial=vec![315.0;mesh.vertex_count()];let before=initial.clone();
        let applied=set.evaluate(cx,&bc,&initial).unwrap();
        let frozen=solve_with_interfaces(cx,problem,applied.frozen_interfaces(),solver_config()).unwrap();
        let actual=set.evaluate(cx,&bc,&frozen.temperature).unwrap();
        let system=assemble_operator_with_interfaces(cx,&mesh,&bc,&material,&source,
            &frozen.temperature,actual.frozen_interfaces()).unwrap();
        let mut residual=vec![0.0;mesh.vertex_count()];
        system.operator.spmv(&frozen.temperature,&mut residual);
        for (r,b) in residual.iter_mut().zip(&system.load) {*r-=b;}
        near(residual.iter().sum(),0.0,1e-8);
        let dofs=DofMap::new(&bc,mesh.vertex_count()).unwrap();
        assert!(dofs.free().iter().any(|&i|residual[i].abs()>1e-3));
        assert!(set.solve_steady(cx,problem,&initial,solver_config(),
            SteadyContactConfig {max_iterations:1,..coupling()}).is_err());
        assert!(set.solve_steady(cx,problem,&vec![400.0;mesh.vertex_count()],
            solver_config(),coupling()).is_err());
        assert!(set.solve_steady(cx,problem,&initial[..initial.len()-1],
            solver_config(),coupling()).is_err());
        assert!(set.solve_steady(cx,problem,&initial,solver_config(),
            SteadyContactConfig {residual_tolerance_w:0.0,..coupling()}).is_err());
        assert_eq!(initial,before);
        assert!(set.solve_steady(cx,problem,&initial,solver_config(),coupling()).is_ok());
    });
}
