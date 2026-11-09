//! Original positive constrained relaxation factor; no formed squared operator.
use super::*;
use fs_la::factor::svd_jacobi;
use fs_math::det;

impl PreparedTensorKelvinPlaneStress {
    /// Bind an original accepted history without advancing it. The transverse
    /// viscosity block may be SPD or have a structurally resolved original
    /// Gram rank below three. Redundant projected columns refuse.
    /// Original dense factor budgets and cooperative cancellation apply.
    pub fn new(history:&TensorPronyHistory,plane:MaterialPlane,
        budget:TensorKelvinStressControlBudget,mut cancelled:impl FnMut()->bool)->Result<Self>
    {
        poll(&mut cancelled)?;
        let model=&history.model;
        let n=3+6*model.branches().len();
        let storage=n.checked_mul(n).and_then(|v|v.checked_mul(16)).and_then(|v|v.checked_add(128*n+2048));
        let work=n.checked_pow(3).and_then(|v|v.checked_mul(1024)).and_then(|v|v.checked_add(100*n*n+20_000));
        if n>budget.max_state_coordinates || storage.is_none_or(|v|v>budget.max_dense_scalars)
            || work.is_none_or(|v|v>budget.max_preparation_work || v>fs_la::eigen::MAX_EIGEN_UNPOLLED_SCALAR_WORK) {
            return Err(invalid("Kelvin plane-stress preparation exceeds its explicit budget"));
        }
        let (inside,outside)=match plane {
            MaterialPlane::Axes12=>([0,1,3],[2,4,5]),MaterialPlane::Axes21=>([1,0,3],[2,4,5]),
            MaterialPlane::Axes13=>([0,2,5],[1,3,4]),MaterialPlane::Axes23=>([1,2,4],[0,3,5]),
        };
        let factor=model.kv_factor.as_ref().ok_or_else(||invalid("Kelvin plane stress requires transverse viscosity"))?;
        let columns:Vec<_>=(0..6).filter(|&column|outside.iter().any(|&axis|factor[axis][column]!=0.0)).collect();
        if columns.len()<3 {
            return Self::prepare_semidefinite(history,plane,inside,outside,&columns,&mut cancelled);
        }
        let gram=core::array::from_fn(|i|core::array::from_fn(|j|factor[j][i]));
        let equilibrium=Root3::new(&model.equilibrium_gram,outside)?;
        let viscosity=Root3::new(&gram,outside)?;
        let mut map=[[0.0;6];3];
        for j in 0..3 {
            map[j][inside[j]]=ENGINEERING[j];
            let rhs=core::array::from_fn(|i|-Dd::from_f64(model.c_eq[outside[i]][inside[j]])*Dd::from_f64(ENGINEERING[j]));
            let extension=equilibrium.solve(equilibrium.solve(rhs,true)?,false)?;
            for i in 0..3 {map[j][outside[i]]=extension[i].to_f64();}
        }
        let mut direct_gram=[[0.0;6];6];
        let mut direct_work=[[0.0;6];6];
        let mut w=[[0.0;3];3];
        let mut dissipative_map=map;
        for j in 0..3 {
            let mut column=apply(&gram,&map[j]);
            viscosity.qr.apply_qt(&mut column);
            let correction=viscosity.solve(core::array::from_fn(|i|Dd::from_f64(column[i])),false)?;
            for i in 0..3 {
                w[j][i]=correction[i].to_f64();
                dissipative_map[j][outside[i]]-=w[j][i];
                direct_gram[i][j]=column[i+3];
            }
        }
        // Independent original-tensor reaction, not the residual Gram square.
        let original_d=model.kv_pa_s().expect("admitted original viscosity");
        for i in 0..3 {for j in 0..3 {
            let mut sum=Sum::default();
            for a in 0..6 {for b in 0..6 {sum.add([map[i][a],original_d[a][b],dissipative_map[j][b],1.0]);}}
            direct_work[i][j]=number(sum.retained_value())?.to_f64();
        }}
        let mut state_rows=vec![[0.0;3];n];
        let mut force=vec![[0.0;3];n];
        let mut factor=zeros(n*n)?;
        for row in 0..3 {
            for i in 0..3 {state_rows[row][i]=equilibrium.scale*equilibrium.qr.r(row,i);}
            for j in 0..3 {force[row][j]=-equilibrium.apply(w[j].map(Dd::from_f64))?[row].to_f64();}
        }
        for (j,branch) in model.branches().iter().enumerate() {
            poll(&mut cancelled)?;
            let root=checked(1.0/det::sqrt(branch.tau_s),"unrepresentable original Maxwell rate root")?;
            if root<=0.0 {return Err(invalid("original Maxwell rate root underflowed"));}
            let b=branch.resolved_gram_rows(&model.factors[j]);
            for row in 0..6 {
                let index=3+6*j+row;
                state_rows[index]=outside.map(|i|b[row][i]);
                for k in 0..3 {
                    force[index][k]=(0..6).fold(Dd::ZERO,|s,i|s+Dd::from_f64(b[row][i])*Dd::from_f64(dissipative_map[k][i])).to_f64();
                }
                factor[index*n+index]=root;
            }
        }
        for column in 0..n {
            let response=viscosity.solve(state_rows[column].map(Dd::from_f64),true)?;
            for row in 0..3 {factor[row*n+column]=response[row].to_f64();}
        }
        let prepared=Self {model:Arc::clone(model),plane,inside,outside,equilibrium,
            viscosity:Some(viscosity),constraint:None,map,direct_gram,direct_work,
            rates:Vec::new(),rate_roots:Vec::new(),modes:Vec::new(),
            forcing:Vec::new(),strain_modes:Vec::new()};
        Self::finish_factor(prepared,factor,&force,n,None,&mut cancelled)
    }

    // One positive native factorization for both transverse ranks. The
    // definite path is the original identity basis with unchanged arithmetic.
    pub(super) fn finish_factor(mut prepared:Self,mut factor:Vec<f64>,force:&[[f64;3]],
        n:usize,basis:Option<&[f64]>,cancelled:&mut impl FnMut()->bool)->Result<Self>
    {
        if n==0 {poll(cancelled)?;return Ok(prepared);}
        let original=prepared.state_coordinates();
        let scale=factor.iter().fold(0.0_f64,|s,v|s.max(v.abs()));
        if !scale.is_finite() || scale<=0.0 || force.iter().flatten().any(|v|!v.is_finite()) {
            return Err(invalid("unrepresentable constrained Kelvin relaxation factor"));
        }
        for v in &mut factor {*v/=scale;}
        poll(cancelled)?;
        let decomposition=svd_jacobi(&factor,n,n);
        poll(cancelled)?;
        let tolerance=8192.0*f64::EPSILON*original as f64;
        let mut modes=zeros(original*n)?;
        for mode in 0..n {
            poll(cancelled)?;
            let singular=decomposition.sigma[mode];
            if !singular.is_finite() || singular<=0.0 {return Err(invalid("constrained Kelvin mode is unresolved"));}
            for row in 0..n {
                let actual=(0..n).fold(Dd::ZERO,|s,j|s+Dd::from_f64(factor[row*n+j])*Dd::from_f64(decomposition.v[j*n+mode]));
                let residual=(actual-Dd::from_f64(singular)*Dd::from_f64(decomposition.u[row*n+mode])).to_f64();
                if !residual.is_finite() || (residual/singular).abs()>tolerance {return Err(invalid("Kelvin plane-stress factor residual is unresolved"));}
            }
            for other in 0..n {for vectors in [&decomposition.u,&decomposition.v] {
                let value=(0..n).fold(Dd::ZERO,|s,row|s+Dd::from_f64(vectors[row*n+mode])*Dd::from_f64(vectors[row*n+other]));
                if (value.to_f64()-if mode==other {1.0}else{0.0}).abs()>tolerance {return Err(invalid("Kelvin plane-stress modes are not resolved orthogonal"));}
            }}
            for row in 0..original {
                modes[row*n+mode]=if let Some(basis)=basis {
                    number((0..n).fold(Dd::ZERO,|s,j|s+Dd::from_f64(basis[row*n+j])*Dd::from_f64(decomposition.v[j*n+mode])))?.to_f64()
                }else{decomposition.v[row*n+mode]};
            }
            let root=checked(scale*singular,"unrepresentable Kelvin plane-stress rate root")?;
            let rate=checked(root*root,"unrepresentable Kelvin plane-stress rate")?;
            if rate<=0.0 {return Err(invalid("Kelvin plane-stress rate underflowed"));}
            prepared.rate_roots.push(root);prepared.rates.push(rate);
            prepared.forcing.push(core::array::from_fn(|j|(0..original).fold(Dd::ZERO,|s,row|s+Dd::from_f64(modes[row*n+mode])*Dd::from_f64(force[row][j])).to_f64()));
            let free=prepared.equilibrium.solve(core::array::from_fn(|i|Dd::from_f64(modes[i*n+mode])),false)?;
            let mut direction=[0.0;6];for i in 0..3 {direction[prepared.outside[i]]=free[i].to_f64();}
            prepared.strain_modes.push(direction);
        }
        prepared.modes=modes;
        poll(cancelled)?;
        Ok(prepared)
    }
}
