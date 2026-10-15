//! A compressed slender solid must not be converted into stable audible modes.
use super::*;

#[test]
fn a_small_strain_thermal_column_retains_negative_modes_across_assembly_batches() {
    let sections=12usize;
    let mut nodes=Vec::new();
    for i in 0..=sections {
        for (y,z) in [(0.,0.),(0.1,0.),(0.,0.1),(0.1,0.1)] {
            nodes.push([4.*i as f64/sections as f64,y,z]);
        }
    }
    let mut cells=Vec::new();
    for i in 0..sections {
        let a=4*i;
        let corners=[a,a+4,a+1,a+5,a+2,a+6,a+3,a+7];
        for tet in [[0,1,3,7],[0,3,2,7],[0,2,6,7],[0,6,4,7],[0,4,5,7],[0,5,1,7]] {
            cells.push(tet.map(|v|corners[v]));
        }
    }
    assert_eq!(cells.len(),72); // both the complete 64-cell batch and its tail
    let fixed:Vec<_>=(0..4).chain(4*sections..4*(sections+1))
        .flat_map(|v|(0..3).map(move |a|3*v+a)).collect();
    let mat=material();
    let problem=TetLinearElasticProblem {nodes_m:&nodes,tetrahedra:&cells,
        materials:TetMaterialField::Uniform(&mat),fixed_dofs:&fixed,
        budget:TetAssemblyBudget {maximum_nodes:52,maximum_tetrahedra:72,
            maximum_free_dofs:156,minimum_scaled_jacobian:1e-10}};
    with_cx(|cx| {
        let original=problem.assemble(cx).unwrap();
        let thermal=strain(0.032);
        let field=TetThermalStrainField::Uniform(&thermal);
        let solution=problem.solve_thermal_displacement(field,TetStaticSolveConfig::default(),cx).unwrap();
        let hot=problem.assemble_thermal_tangent(&solution,field,ContentHash([1;32]),0.05,cx).unwrap();
        assert_eq!(hot.pencil().mass.to_dense(),original.mass.to_dense());
        assert_eq!(hot.elastic_stiffness().to_dense(),original.stiffness.to_dense());
        let stable=eigh_gen_dense(&original.stiffness.to_dense(),&original.mass.to_dense(),original.free_dofs.len()).unwrap();
        assert!(stable.iter().all(|mode|mode.lambda>0.));
        let heated=eigen(&hot);
        assert_eq!(heated.iter().filter(|&&lambda|lambda<0.).count(),2);
        // Independent P1 assembly: [-0.50288707, -0.22087328, 12.32365743,...].
        // This checks the discrete model, not continuum buckling convergence.
        assert!((heated[0]+0.50288707).abs()<1e-5);
        assert!((heated[1]+0.22087328).abs()<1e-5);
        assert!(hot.stress().elements().iter().all(|e|
            e.elastic_strain_mandel().iter().all(|v|v.abs()<0.05)));
    });
}
