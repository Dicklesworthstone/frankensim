use super::*;

fn diagonal(value: f64) -> Tangent6 { core::array::from_fn(|i| core::array::from_fn(|j| if i==j {value} else {0.0})) }
fn model() -> TensorProny {
    let mut a = [[0.0;6];6]; a[0] = [2.0,-1.0,0.0,3.0,0.5,0.0];
    let mut b = [[0.0;6];6]; b[2] = [0.0,1.0,2.0,0.0,-1.0,0.5];
    TensorProny::new_mandel(diagonal(12.0), vec![
        TensorPronyBranch::from_gram_rows(a,0.01).unwrap(),
        TensorPronyBranch::from_gram_rows(b,0.3).unwrap()],0.1).unwrap()
}
fn columns() -> Vec<Voigt> { vec![[1.0,0.0,-0.2,0.3,0.0,0.0],[-0.1,0.5,0.0,-0.4,0.2,0.0]] }
fn budget() -> TensorProjectionBudget { TensorProjectionBudget {max_coordinates:4,max_scalar_terms:10_000} }
fn close(a:f64,b:f64) { assert!((a-b).abs() < 2e-12*a.abs().max(b.abs()).max(1e-12),"{a} vs {b}"); }
fn projected_stiffness(p:&ProjectedTensorProny,omega:f64)->(Vec<f64>,Vec<f64>) {
    let n=p.coordinates(); let mut re=p.equilibrium_stiffness().to_vec(); let mut im=vec![0.0;n*n];
    for b in p.branches() {
        let z=omega*b.relaxation_time_s(); let den=1.0+z*z;
        for i in 0..n {for j in 0..n {
            let value=b.rows().iter().fold(0.0,|s,row|s+row[i]*row[j]);
            re[i*n+j]+=value*z*z/den; im[i*n+j]+=value*z/den;
        }}
    }
    (re,im)
}

#[test]
fn projected_tensor_preserves_full_complex_work_and_rank_deficient_branches() {
    let m=model(); let c=columns(); let p=m.project_strain(&c,0.25,budget(),||false).unwrap();
    assert_eq!(p.coordinates(),2); assert_eq!(p.branches().len(),2);
    assert!(p.branches()[0].rows()[1..].iter().flatten().all(|x|*x==0.0));
    for omega in [0.0,0.1,10.0,100.0,10_000.0] {
        let (storage,loss)=m.modulus_mandel(omega).unwrap();
        let (re,im)=projected_stiffness(&p,omega);
        for i in 0..2 {for j in 0..2 {
            close(re[2*i+j],0.25*dot(&mandel(&c[i]),&apply(&storage,&mandel(&c[j]))));
            close(im[2*i+j],0.25*dot(&mandel(&c[i]),&apply(&loss,&mandel(&c[j]))));
        }}
    }
    let x=[0.002,-0.003]; let strain=p.strain(&x).unwrap();
    let mut history=m.history(); let report=history.advance(strain,0.0,||false).unwrap();
    let force=p.conjugate_force(&report.stress).unwrap();
    close(0.5*(force[0]*x[0]+force[1]*x[1]),0.25*report.stored_energy);
    assert_ne!(p.equilibrium_stiffness()[1],0.0,"cross terms must survive");
}

#[test]
fn projected_factor_histories_match_the_original_material_on_reversing_strains() {
    let m=model(); let p=m.project_strain(&columns(),0.25,budget(),||false).unwrap();
    let mut history=m.history(); let mut z=vec![[0.0;6];p.branches().len()];
    let mut old=[0.0;2];
    for step in 0..100 {
        let x=[0.004*det::sin(step as f64*0.37),0.002*det::cos(step as f64*0.29)];
        let dt=if step==0 {0.0} else {0.003};
        let native=history.advance(p.strain(&x).unwrap(),dt,||false).unwrap();
        let mut force:Vec<_>=(0..2).map(|i|(0..2).map(|j|p.equilibrium_stiffness()[2*i+j]*x[j]).sum::<f64>()).collect();
        let mut stored=0.5*(force[0]*x[0]+force[1]*x[1]); let mut heat=0.0;
        for (z,branch) in z.iter_mut().zip(p.branches()) {
            for (k,row) in branch.rows().iter().enumerate() {
                let delta=row.iter().enumerate().fold(0.0,|s,(i,a)|s+a*(x[i]-old[i]));
                let (next,_,loss)=ramp(z[k],delta,dt/branch.relaxation_time_s());
                z[k]=next;heat+=loss;stored+=0.5*next*next;
                for i in 0..2 {force[i]+=row[i]*next;}
            }
        }
        let expected=p.conjugate_force(&native.stress).unwrap();
        for i in 0..2 {close(force[i],expected[i]);}
        close(stored,0.25*native.stored_energy); close(heat,0.25*native.dissipated);
        old=x;
    }
}

#[test]
fn projection_preserves_shear_work_and_volume_scaling() {
    let m=TensorProny::new_mandel(diagonal(2.0),vec![],0.1).unwrap();
    let c=[[0.0,0.0,0.0,1.0,0.0,0.0]];
    let p=m.project_strain(&c,0.5,budget(),||false).unwrap();
    close(p.equilibrium_stiffness()[0],2.0); // 2*V*C_xyxy, not V*C_xyxy.
    let large=m.project_strain(&c,2.0,budget(),||false).unwrap();
    close(large.equilibrium_stiffness()[0],4.0*p.equilibrium_stiffness()[0]);
    let stress=[0.0,0.0,0.0,7.0,0.0,0.0];
    close(p.conjugate_force(&stress).unwrap()[0],7.0);
    assert!(p.strain(&[0.08]).is_err(),"Frobenius shear norm exceeds 0.1");
}

#[test]
fn projection_refuses_domain_work_and_cancellation_without_mutating_the_law() {
    let m=model(); let saved=m.clone(); let c=columns();
    for v in [0.0,-1.0,f64::NAN,f64::INFINITY] {assert!(m.project_strain(&c,v,budget(),||false).is_err());}
    assert!(m.project_strain(&[],1.0,budget(),||false).is_err());
    assert!(m.project_strain(&c,1.0,TensorProjectionBudget {max_scalar_terms:1,..budget()},||false).is_err());
    let mut bad=c.clone();bad[0][0]=f64::NAN;assert!(m.project_strain(&bad,1.0,budget(),||false).is_err());
    let mut polls=0;let baseline=m.project_strain(&c,1.0,budget(),||{polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut at=0;assert!(matches!(m.project_strain(&c,1.0,budget(),||{at+=1;at==stop}),Err(TensorPronyError::Cancelled)));
    }
    assert_eq!(m,saved);
    assert_eq!(baseline,m.project_strain(&c,1.0,budget(),||false).unwrap());
    assert!(baseline.strain(&[f64::NAN,0.0]).is_err());
    assert!(baseline.strain(&[0.0]).is_err());
    assert!(baseline.conjugate_force(&[f64::NAN;6]).is_err());
}

#[test]
fn thermal_projection_retains_signed_full_tensor_derivatives_and_original_strain_work() {
    let m = model();
    let columns = columns();
    let volume = 0.25;
    let p = m.project_strain(&columns,volume,budget(),||false).unwrap();
    let mut equilibrium = diagonal(-0.7);
    equilibrium[0][3] = 1.3; equilibrium[3][0] = 1.3;
    let mut first = diagonal(0.2);
    first[0][0] = -1.1; first[1][4] = -0.8; first[4][1] = -0.8;
    let mut second = diagonal(-0.3);
    second[2][5] = 0.9; second[5][2] = 0.9;
    let derivative = source::thermal::TensorTemperatureDerivative::from_affine_parts(
        equilibrium, vec![first,second], [280.0,320.0],
    );
    let projected = p.thermal_stiffness_derivatives(&derivative,budget(),||false).unwrap();
    let pairs = [(equilibrium,projected.equilibrium()),
        (first,projected.branches()[0].as_slice()),
        (second,projected.branches()[1].as_slice())];
    // Independent physical Voigt contraction: both shear factors sqrt(2)
    // remain, including every normal/shear cross term and both tangent axes.
    let weights = [1.0,1.0,1.0,2.0_f64.sqrt(),2.0_f64.sqrt(),2.0_f64.sqrt()];
    for (matrix, reduced) in pairs {
        for x in [[0.004,-0.003],[-0.002,0.005]] {
            let strain: [f64;6] = std::array::from_fn(|i|
                weights[i]*(columns[0][i]*x[0]+columns[1][i]*x[1]));
            let direct = volume*(0..6).map(|i| (0..6)
                .map(|j| strain[i]*matrix[i][j]*strain[j]).sum::<f64>()).sum::<f64>();
            let reduced_work = reduced[0]*x[0]*x[0] + 2.0*reduced[1]*x[0]*x[1]
                + reduced[3]*x[1]*x[1];
            close(direct,reduced_work);
            assert_eq!(reduced[1],reduced[2]);
        }
    }
    let viscous = m.with_kelvin_viscosity_mandel(diagonal(0.2)).unwrap()
        .project_strain(&columns,volume,budget(),||false).unwrap();
    assert_eq!(projected,viscous.thermal_stiffness_derivatives(&derivative,budget(),||false).unwrap());
    assert!(p.thermal_stiffness_derivatives(&derivative,
        TensorProjectionBudget {max_scalar_terms:503,..budget()},||false).is_err());
    assert_eq!(projected,p.thermal_stiffness_derivatives(&derivative,
        TensorProjectionBudget {max_scalar_terms:504,..budget()},||false).unwrap());
    assert!(matches!(p.thermal_stiffness_derivatives(&derivative,budget(),||true),
        Err(TensorPronyError::Cancelled)));
    let mut incomplete = derivative;
    incomplete.branches_pa_k = incomplete.branches_pa_k[..1].to_vec().into();
    assert!(p.thermal_stiffness_derivatives(&incomplete,budget(),||false).is_err());
}

fn deviatoric_rows(scale: f64) -> Tangent6 {
    let mut rows = [[0.0;6];6];
    let a = scale/det::sqrt(2.0);
    let b = scale/det::sqrt(6.0);
    rows[0] = [a,-a,0.0,0.0,0.0,0.0];
    rows[1] = [b,b,-2.0*b,0.0,0.0,0.0];
    for i in 2..5 { rows[i][i+1] = scale; }
    rows
}
fn factor_projection(rows: Tangent6, free: Voigt) -> ProjectedTensorProny {
    TensorProny::new_mandel(diagonal(12.0), vec![
        TensorPronyBranch::from_gram_rows(rows,0.1).unwrap()],0.1).unwrap()
        .project_strain(&columns(),0.25,budget(),||false).unwrap()
        .with_eigenstrain_mandel(free).unwrap()
}

#[test]
fn changing_rank_five_factors_preserve_memory_and_thermal_eigenstrain() {
    let rows = deviatoric_rows(2.0);
    let mut next_rows = rows;
    next_rows[0] = core::array::from_fn(|i| 1.2*rows[0][i]+0.4*rows[1][i]);
    next_rows[1] = core::array::from_fn(|i| -0.2*rows[0][i]+0.8*rows[1][i]);
    next_rows[2] = rows[2].map(|v| 1.4*v);
    next_rows[3] = rows[3].map(|v| 0.7*v);
    let free = [1e-4,2e-4,-1e-4,0.0,3e-5,-2e-5];
    let old = factor_projection(rows,[0.0;6]);
    let next = factor_projection(next_rows,free);
    // The hydrostatic part is invisible to this shear arm. All five
    // observable directions, including signed normal/shear terms, remain.
    let q = [0.018,0.006,0.012,0.003,-0.002,0.001];
    let memory = apply(&rows,&q).map(|x| 0.5*x);
    let changed = old.transfer_factor_memory(&next,&memory,1e-11,||false).unwrap();
    let expected = apply(&next_rows,&core::array::from_fn(|i| q[i]-free[i])).map(|x| 0.5*x);
    for (actual,want) in changed.iter().zip(expected) {
        assert!((actual-want).abs() < 1e-15,"{actual} vs {want}");
    }
    let returned = next.transfer_factor_memory(&old,&changed,1e-11,||false).unwrap();
    for (actual,want) in returned.iter().zip(memory) { assert!((actual-want).abs() < 1e-15); }
    assert_eq!(memory[5],0.0);
    assert_eq!(changed[5],0.0);
    assert!(matches!(old.transfer_factor_memory(&next,&memory,1e-11,||true),
        Err(TensorPronyError::Cancelled)));
}

#[test]
fn singular_affine_thermal_memory_retains_entropy_and_expansion_capacity() {
    let rows = deviatoric_rows(2.0);
    let point = factor_projection(rows,[0.0;6]);
    let c = gram(&rows).unwrap();
    let slope = c.map(|row| row.map(|v| -0.005*v));
    let derivative = source::thermal::TensorTemperatureDerivative::from_affine_parts(
        [[0.0;6];6],vec![slope],[280.0,320.0]);
    let alpha = [1e-4,2e-4,-1e-4,3e-5,-2e-5,4e-5];
    let beta = [1e-6,-2e-6,1e-6,3e-7,2e-7,-4e-7];
    let expansion = TensorThermalExpansion::try_new(300.0,alpha,beta,[280.0,320.0]).unwrap();
    let x = [0.003,-0.001];
    let e = mandel(&point.strain(&x).unwrap());
    let q = [0.018,0.006,0.012,0.003,-0.002,0.001];
    let memory = apply(&rows,&q).map(|v| 0.5*v);
    let state = point.thermal_state_with_expansion_derivatives(
        &derivative,&x,&memory,300.0,&expansion,1e-11,||false).unwrap();
    // Independent physical-strain derivatives of .5 V(e:Ceq:e+q:C:q).
    let equilibrium = diagonal(12.0);
    let stress = core::array::from_fn(|i| apply(&equilibrium,&e)[i]+apply(&c,&q)[i]);
    let free = 0.125*(dot(&e,&apply(&equilibrium,&e))+dot(&q,&apply(&c,&q)));
    let entropy = 0.25*(dot(&alpha,&stress)-0.5*dot(&q,&apply(&slope,&q)));
    let curvature = 0.25*(dot(&alpha,&apply(&equilibrium,&alpha))
        +dot(&alpha,&apply(&c,&alpha))-2.0*dot(&alpha,&apply(&slope,&q))-dot(&beta,&stress));
    close(state.free_energy_j,free);
    close(state.entropy_j_k,entropy);
    close(state.internal_energy_j,free+300.0*entropy);
    close(state.heat_capacity_j_k,-300.0*curvature);
}

#[test]
fn thermal_transfer_retains_weak_positive_rows_and_refuses_rank_or_range_changes() {
    let mut rows = [[0.0;6];6]; rows[0][0] = 100.0; rows[1][1] = 1e-100;
    let mut changed_rows = rows; changed_rows[0][0] = 200.0; changed_rows[1][1] = 3e-100;
    let old = factor_projection(rows,[0.0;6]);
    let next = factor_projection(changed_rows,[0.0;6]);
    let q = [0.01,0.02,0.0,0.0,0.0,0.0];
    let memory = apply(&rows,&q).map(|v| 0.5*v);
    let result = old.transfer_factor_memory(&next,&memory,1e-11,||false).unwrap();
    for (actual,expected) in result.iter().zip(apply(&changed_rows,&q).map(|v| 0.5*v)) {
        if expected == 0.0 { assert_eq!(*actual,0.0); }
        else { assert!((actual/expected-1.0).abs() < 1e-12); }
    }
    assert!(result[1] > 0.0,"weak physical memory must survive");
    changed_rows[1][1] = 0.0;
    let dropped = factor_projection(changed_rows,[0.0;6]);
    assert!(old.transfer_factor_memory(&dropped,&[0.0;6],1e-11,||false).is_err());
    assert!(dropped.transfer_factor_memory(&old,&[0.0;6],1e-11,||false).is_err());
    changed_rows[1][2] = 3e-100;
    let rotated_nullspace = factor_projection(changed_rows,[0.0;6]);
    assert!(old.transfer_factor_memory(&rotated_nullspace,&[0.0;6],1e-11,||false).is_err());
}

#[test]
fn hidden_thermal_directions_and_unrecoverable_factor_memory_refuse_transactionally() {
    let rows = deviatoric_rows(2.0);
    let old = factor_projection(rows,[0.0;6]);
    let saved = old.clone();
    let next = factor_projection(rows.map(|row| row.map(|v| 2.0*v)),[0.0;6]);
    let mut memory = [0.0;6]; memory[5] = 1e-100;
    assert!(old.transfer_factor_memory(&next,&memory,1e-11,||false).is_err());
    let mut redundant = rows; redundant[5] = rows[0];
    assert!(old.transfer_factor_memory(&factor_projection(redundant,[0.0;6]),
        &[0.0;6],1e-11,||false).is_err());
    let hidden = core::array::from_fn(|i| core::array::from_fn(|j|
        if i<3 && j<3 {1e-100} else {0.0}));
    let derivative = source::thermal::TensorTemperatureDerivative::from_affine_parts(
        [[0.0;6];6],vec![hidden],[280.0,320.0]);
    assert!(old.thermal_state(&derivative,&[0.0;2],&[0.0;6],300.0,1e-11,||false).is_err());
    assert_eq!(memory[5],1e-100);
    assert_eq!(old,saved);
}

#[test]
fn expansion_capacity_matches_internal_energy_at_fixed_total_and_viscous_strain() {
    let column = [1.0, 0.2, -0.4, 0.3, -0.1, 0.15];
    let displacement = [0.002];
    let volume = 0.25;
    let alpha = [1e-5, -2e-5, 0.5e-5, 0.0, 0.0, 0.0];
    let viscous = [0.001, -0.0002, 0.0003, -0.0001, 0.0, 0.0002];
    let total = mandel(&column).map(|v| v * displacement[0]);
    let derivative = source::thermal::TensorTemperatureDerivative::from_affine_parts(
        diagonal(-2.0), vec![diagonal(3.0)], [280.0, 320.0],
    );
    // Independent scalar continuum potential. Both the total strain and the
    // physical viscous strain stay fixed as temperature changes; factor
    // coordinates themselves must change with branch stiffness and free strain.
    let reference = |temperature: f64| {
        let dt = temperature - 300.0;
        let equilibrium = 1200.0 - 2.0 * dt;
        let branch = 400.0 + 3.0 * dt;
        let elastic: [f64; 6] = core::array::from_fn(|i| total[i] - alpha[i] * dt);
        let spring: [f64; 6] = core::array::from_fn(|i| elastic[i] - viscous[i]);
        let free = 0.5 * volume * (0..6).map(|i|
            equilibrium * elastic[i] * elastic[i] + branch * spring[i] * spring[i]
        ).sum::<f64>();
        let entropy = volume * (0..6).map(|i|
            alpha[i] * (equilibrium * elastic[i] + branch * spring[i])
                + elastic[i] * elastic[i] - 1.5 * spring[i] * spring[i]
        ).sum::<f64>();
        let capacity = -temperature * volume * (0..6).map(|i|
            (equilibrium + branch) * alpha[i] * alpha[i]
                - 2.0 * alpha[i] * (-2.0 * elastic[i] + 3.0 * spring[i])
        ).sum::<f64>();
        (free, entropy, free + temperature * entropy, capacity, spring)
    };
    for temperature in [285.0, 300.0, 315.0] {
        let dt = temperature - 300.0;
        let branch = 400.0 + 3.0 * dt;
        let law = TensorProny::new_mandel(
            diagonal(1200.0 - 2.0 * dt),
            vec![TensorPronyBranch::from_gram_rows(diagonal(det::sqrt(branch)), 0.2).unwrap()],
            0.1,
        ).unwrap();
        let projected = law.project_strain(&[column], volume, budget(), || false).unwrap()
            .with_eigenstrain_mandel(alpha.map(|a| a * dt)).unwrap();
        let (free, entropy, internal, capacity, spring) = reference(temperature);
        let memory = spring.map(|q| det::sqrt(volume) * det::sqrt(branch) * q);
        let state = projected.thermal_state_with_expansion(
            &derivative, &displacement, &memory, temperature, alpha, 1e-10, || false,
        ).unwrap();
        close(state.free_energy_j, free);
        close(state.entropy_j_k, entropy);
        close(state.internal_energy_j, internal);
        close(state.heat_capacity_j_k, capacity);
        let step = 0.001;
        let numerical = (reference(temperature + step).2 - reference(temperature - step).2)
            / (2.0 * step);
        assert!((state.heat_capacity_j_k - numerical).abs()
            < 5e-8 * capacity.abs().max(1e-8));
        let bound = projected.thermal_capacity_bound_j_k2(&derivative, temperature, alpha).unwrap();
        assert!(bound >= -state.heat_capacity_j_k / temperature);
        close(bound, projected.affine_expansion_capacity_bound_j_k2(
            &derivative, temperature, alpha,
        ).unwrap());
    }
}


#[test]
fn variable_expansion_curvature_matches_energy_derivatives_with_singular_memory() {
    let volume: f64 = 0.25;
    let beta = [2e-7, -1e-7, 1.5e-7, 0.0, 0.0, 0.0];
    let expansion = TensorThermalExpansion::try_new(
        300.0, [0.0; 6], beta, [280.0, 320.0],
    ).unwrap();
    let derivative = source::thermal::TensorTemperatureDerivative::from_affine_parts(
        diagonal(0.0), vec![diagonal(0.0)], [280.0, 320.0],
    );
    let row = [20.0, -10.0, 5.0, 0.0, 0.0, 0.0];
    let mut rows = [[0.0; 6]; 6];
    rows[0] = row;
    let native = TensorProny::new_mandel(
        diagonal(1200.0),
        vec![TensorPronyBranch::from_gram_rows(rows, 0.2).unwrap()], 0.1,
    ).unwrap();
    let column = [1.0, 0.2, -0.4, 0.3, -0.1, 0.15];
    let total = mandel(&column).map(|x| 0.002 * x);
    let viscous = [0.001, -0.0002, 0.0003, -0.0001, 0.0, 0.0002];
    // Direct continuum potential at fixed total and viscous strains. This
    // reference uses the original rank-one spring, never a factor inverse.
    let reference = |t: f64| {
        let dt = t - 300.0;
        let alpha = beta.map(|b| b * dt);
        let free = beta.map(|b| 0.5 * b * dt * dt);
        let elastic: Voigt = core::array::from_fn(|i| total[i] - free[i]);
        let spring: Voigt = core::array::from_fn(|i| elastic[i] - viscous[i]);
        let spring_projection = (0..6).map(|i| row[i] * spring[i]).sum::<f64>();
        let stress: Voigt = core::array::from_fn(|i|
            1200.0 * elastic[i] + row[i] * spring_projection);
        let psi = 0.5 * volume * (1200.0 * elastic.iter().map(|x| x*x).sum::<f64>()
            + spring_projection * spring_projection);
        let entropy = volume * (0..6).map(|i| alpha[i] * stress[i]).sum::<f64>();
        let row_alpha = (0..6).map(|i| row[i] * alpha[i]).sum::<f64>();
        let curvature = volume * (
            1200.0 * alpha.iter().map(|a| a*a).sum::<f64>() + row_alpha * row_alpha
            - (0..6).map(|i| beta[i] * stress[i]).sum::<f64>());
        (psi, entropy, psi + t * entropy, -t * curvature, free, spring_projection)
    };
    for temperature in [285.0, 300.0, 315.0] {
        let (psi, entropy, internal, capacity, free, spring_projection) = reference(temperature);
        let projected = native.project_strain(&[column], volume, budget(), || false).unwrap()
            .with_eigenstrain_mandel(free).unwrap();
        let mut memory = [0.0; 6];
        memory[0] = det::sqrt(volume) * spring_projection;
        let original_memory = memory;
        let snapshot = expansion.at_temperature(temperature).unwrap();
        let actual = projected.thermal_state_with_expansion_derivatives(
            &derivative, &[0.002], &memory, temperature, &snapshot, 1e-10, || false,
        ).unwrap();
        close(actual.free_energy_j, psi);
        if entropy == 0.0 { assert_eq!(actual.entropy_j_k, 0.0); }
        else { close(actual.entropy_j_k, entropy); }
        close(actual.internal_energy_j, internal);
        close(actual.heat_capacity_j_k, capacity);
        let numerical = (reference(temperature + 0.001).2 - reference(temperature - 0.001).2) / 0.002;
        assert!((numerical-capacity).abs() < 5e-8 * capacity.abs().max(1e-8));
        let bound = projected.thermal_capacity_bound_with_expansion_j_k2(
            &derivative, temperature, &snapshot,
        ).unwrap();
        assert!(bound >= -actual.heat_capacity_j_k / temperature);
        assert_eq!(memory, original_memory);
        if temperature == 300.0 {
            assert_eq!(snapshot.alpha_mandel_per_k(), [0.0; 6]);
            assert_ne!(actual.heat_capacity_j_k, 0.0,
                "alpha crosses zero but its slope still contributes mechanical heat capacity");
        } else {
            assert!(projected.thermal_state_with_expansion_derivatives(
                &derivative, &[0.002], &memory, temperature, &expansion, 1e-10, || false,
            ).is_err(), "stale expansion derivative must not be reused");
        }
    }
    assert!(expansion.at_temperature(330.0).is_err());
}

#[test]
fn variable_expansion_admits_only_the_observable_singular_memory_envelope() {
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = 10.0;
    rows[1][1] = 20.0;
    let native = TensorProny::new_mandel(
        diagonal(1200.0),
        vec![TensorPronyBranch::from_gram_rows(rows, 0.2).unwrap()], 0.1,
    ).unwrap();
    let point = native.project_strain(
        &[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]], 1.0, budget(), || false,
    ).unwrap();
    let derivative = source::thermal::TensorTemperatureDerivative::from_affine_parts(
        diagonal(0.0), vec![diagonal(0.0)], [280.0, 320.0],
    );
    let expansion = TensorThermalExpansion::try_new(
        300.0, [0.0; 6], [1e-6, 0.0, 0.0, 0.0, 0.0, 0.0], [280.0, 320.0],
    ).unwrap();
    // A retained spring force orthogonal to beta has no curvature work.
    // No hidden strain is reconstructed from this singular material tensor.
    let orthogonal_memory = [0.0, 1000.0, 0.0, 0.0, 0.0, 0.0];
    let state = point.thermal_state_with_expansion_derivatives(
        &derivative, &[0.0], &orthogonal_memory, 300.0, &expansion, 1e-10, || false,
    ).unwrap();
    assert_eq!(state.free_energy_j, 500_000.0);
    assert_eq!(state.heat_capacity_j_k, 0.0);
    let outside = [1000.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(point.thermal_state_with_expansion_derivatives(
        &derivative, &[0.0], &outside, 300.0, &expansion, 1e-10, || false,
    ).is_err());
    // Constant alpha retains its original unrestricted constant-branch path.
    assert!(point.thermal_state(
        &derivative, &[0.0], &outside, 300.0, 1e-10, || false,
    ).is_ok());
    let wrong_interval = TensorThermalExpansion::try_new(
        300.0, [0.0; 6], [1e-6, 0.0, 0.0, 0.0, 0.0, 0.0], [290.0, 310.0],
    ).unwrap();
    assert!(point.thermal_capacity_bound_with_expansion_j_k2(
        &derivative, 300.0, &wrong_interval,
    ).is_err());
}
