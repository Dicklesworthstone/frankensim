//! Algebraic transverse constraints of structurally resolved Kelvin factors.
//!
//! Only exactly zero projected original Gram columns are removed. The remaining
//! zero/one/two independent columns establish the admitted rank; a redundant or
//! unresolved projected rank refuses. In z=[G_F*v;q], A=C*N imposes A^T*z=0.
//! Native QR supplies E spanning ker(A^T); F*E is the positive dynamic factor.
//! This eliminates algebraic equilibrium coordinates, never original branch
//! memories. Endpoint rates use the differentiated original algebraic constraint.

use super::*;
use fs_math::det;

pub(super) struct ConstrainedPlaneStress {
    pub(super) rank:usize,
    range:Vec<[f64;3]>,
    null:Vec<[f64;3]>,
    compliance:Vec<[f64;3]>,
    tolerance:f64,
}

fn dot3(a:&[f64;3],b:&[Dd;3])->Dd {
    (0..3).fold(Dd::ZERO,|s,i|s+Dd::from_f64(a[i])*b[i])
}
fn triangular(qr:&Qr,size:usize,rhs:[Dd;3],transpose:bool,scale:f64)->Result<[Dd;3]> {
    let mut answer=[Dd::ZERO;3];
    for index in 0..size {
        let i=if transpose {index}else{size-1-index};
        let mut value=rhs[i]/Dd::from_f64(scale);
        for j in 0..size {
            if (transpose && j<i)||(!transpose && j>i) {
                value=value-Dd::from_f64(if transpose {qr.r(j,i)}else{qr.r(i,j)})*answer[j];
            }
        }
        answer[i]=number(value/Dd::from_f64(qr.r(i,i)))?;
    }
    Ok(answer)
}
fn pivots(qr:&Qr,size:usize)->Result<()> {
    if (0..size).any(|i|!qr.r(i,i).is_finite()||qr.r(i,i).abs()<=64.0*f64::EPSILON) {
        return Err(invalid("original projected Kelvin rank is redundant or unresolved"));
    }
    Ok(())
}

impl ConstrainedPlaneStress {
    /// Compatibility is checked on the original spring stress. No projection
    /// or instantaneous strain jump is inserted into a continuous observation.
    pub(super) fn admit_history(&self,control:&PreparedTensorKelvinPlaneStress,history:&TensorPronyHistory)->Result<()> {
        let model=&control.model;
        let equilibrium=apply(&model.equilibrium_gram,&history.strain);
        let mut residual=[Dd::ZERO;3];let mut magnitude=[0.0;3];
        for i in 0..3 {
            let axis=control.outside[i];
            for row in 0..6 {
                let term=Dd::from_f64(model.equilibrium_gram[row][axis])*Dd::from_f64(equilibrium[row]);
                residual[i]=residual[i]+term;magnitude[i]+=term.abs().to_f64();
            }
            for (j,branch) in model.branches().iter().enumerate() {
                let gram=branch.resolved_gram_rows(&model.factors[j]);
                for row in 0..6 {
                    let term=Dd::from_f64(gram[row][axis])*Dd::from_f64(history.q[j][row]);
                    residual[i]=residual[i]+term;magnitude[i]+=term.abs().to_f64();
                }
            }
        }
        for direction in &self.null {
            let scale=(0..3).map(|i|direction[i].abs()*magnitude[i]).sum::<f64>();
            if !scale.is_finite()||dot3(direction,&residual).abs().to_f64()>self.tolerance*scale {
                return Err(invalid("initial history violates original Kelvin algebraic plane stress"));
            }
        }
        Ok(())
    }

    pub(super) fn endpoint_rate(&self,control:&PreparedTensorKelvinPlaneStress,
        spring:&[Dd;6],known_rate:&[Dd;6],q:&[Voigt])->Result<[Dd;3]>
    {
        let model=&control.model;let viscosity=model.kv_pa_s().expect("original Kelvin tensor");
        let mut residual=control.outside.map(|i|spring[i]);
        for i in 0..3 {for &axis in &control.inside {
            residual[i]=residual[i]+Dd::from_f64(viscosity[control.outside[i]][axis])*known_rate[axis];
        }}
        let mut rate=[Dd::ZERO;3];
        for direction in &self.range {
            let value=dot3(direction,&residual);
            for i in 0..3 {rate[i]=rate[i]-Dd::from_f64(direction[i])*value;}
        }
        // In ker(D_FF), differentiate the actual zero spring traction:
        // N^T*(C0*epsilon' - sum B_j^T*q_j/tau_j)=0.
        let mut c0=model.c_eq;
        for branch in model.branches() {for i in 0..6 {for j in 0..6 {c0[i][j]+=branch.c_pa[i][j];}}}
        let mut full=*known_rate;for i in 0..3 {full[control.outside[i]]=rate[i];}
        let mut null_rhs=[Dd::ZERO;3];
        for i in 0..3 {
            let axis=control.outside[i];let mut sum=Sum::default();
            for (j,branch) in model.branches().iter().enumerate() {
                let b=branch.resolved_gram_rows(&model.factors[j]);
                for row in 0..6 {sum.add_ratio(&[b[row][axis],q[j][row]],&[branch.tau_s]);}
            }
            for j in 0..6 {add(&mut sum,full[j],-c0[axis][j]);}
            null_rhs[i]=number(sum.retained_value())?;
        }
        for root in &self.compliance {
            let value=dot3(root,&null_rhs);
            for i in 0..3 {rate[i]=rate[i]+Dd::from_f64(root[i])*value;}
        }
        for i in 0..3 {full[control.outside[i]]=number(rate[i])?;}
        let mut differentiated=[Dd::ZERO;3];let mut magnitude=[0.0;3];
        for i in 0..3 {
            let axis=control.outside[i];
            for j in 0..6 {let term=Dd::from_f64(c0[axis][j])*full[j];
                differentiated[i]=differentiated[i]+term;magnitude[i]+=term.abs().to_f64();}
            for (j,branch) in model.branches().iter().enumerate() {
                let b=branch.resolved_gram_rows(&model.factors[j]);
                for row in 0..6 {
                    let mut term=Sum::default();term.add_ratio(&[b[row][axis],q[j][row]],&[branch.tau_s]);
                    let term=number(term.retained_value())?;
                    differentiated[i]=differentiated[i]-term;magnitude[i]+=term.abs().to_f64();
                }
            }
        }
        for direction in &self.null {
            let scale=(0..3).map(|i|direction[i].abs()*magnitude[i]).sum::<f64>();
            if !scale.is_finite()||dot3(direction,&differentiated).abs().to_f64()>self.tolerance*scale {
                return Err(invalid("original Kelvin algebraic strain rate is unresolved"));
            }
        }
        Ok(rate)
    }
}

impl PreparedTensorKelvinPlaneStress {
    // The public constructor has already admitted the original full-state
    // storage/work caps. This route never guesses a smaller positive rank.
    pub(super) fn prepare_semidefinite(history:&TensorPronyHistory,plane:MaterialPlane,
        inside:[usize;3],outside:[usize;3],columns:&[usize],cancelled:&mut impl FnMut()->bool)->Result<Self>
    {
        poll(cancelled)?;
        let model=&history.model;let original=3+6*model.branches().len();
        let rank=columns.len();let algebraic=3-rank;let dynamic=original-algebraic;
        let source=model.kv_factor.as_ref().expect("admitted original Kelvin factor");
        if source.iter().flatten().all(|v|*v==0.0) {
            return Err(invalid("zero Kelvin viscosity uses the original Maxwell controller"));
        }
        let tolerance=8192.0*f64::EPSILON*(original+6) as f64;
        let mut basis3=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        let mut range=Vec::with_capacity(rank);
        if rank>0 {
            let scale=outside.iter().flat_map(|&i|columns.iter().map(move |&j|source[i][j].abs())).fold(0.0_f64,f64::max);
            if !scale.is_finite()||scale<=0.0 {return Err(invalid("unrepresentable original transverse Kelvin scale"));}
            let thin:Vec<_>=outside.iter().flat_map(|&i|columns.iter().map(move |&j|source[i][j]/scale)).collect();
            let factor=qr(&thin,3,rank);pivots(&factor,rank)?;
            for vector in &mut basis3 {factor.apply_q(vector);}
            for column in 0..rank {
                let mut rhs=[Dd::ZERO;3];rhs[column]=Dd::ONE;
                let coefficients=triangular(&factor,rank,rhs,true,scale)?;
                let direction=core::array::from_fn(|axis|
                    (0..rank).fold(Dd::ZERO,|s,j|s+Dd::from_f64(basis3[j][axis])*coefficients[j]).to_f64());
                if direction.iter().any(|v|!v.is_finite()) {return Err(invalid("unrepresentable original Kelvin inverse root"));}
                range.push(direction);
            }
        }
        let null=basis3[rank..].to_vec();
        for i in 0..3 {for j in 0..3 {
            if (dot3(&basis3[i],&basis3[j].map(Dd::from_f64)).to_f64()-if i==j {1.0}else{0.0}).abs()>tolerance {
                return Err(invalid("original transverse Kelvin basis is unresolved"));
            }
        }}
        for (j,&column) in columns.iter().enumerate() {
            let vector=outside.map(|i|Dd::from_f64(source[i][column]));
            for (i,direction) in range.iter().enumerate() {
                if (dot3(direction,&vector).to_f64()-if i==j {1.0}else{0.0}).abs()>tolerance {
                    return Err(invalid("original transverse Kelvin rank is unresolved"));
                }
            }
            for direction in &null {
                let scale=(0..3).map(|i|direction[i].abs()*vector[i].abs().to_f64()).sum::<f64>();
                if !scale.is_finite()||dot3(direction,&vector).abs().to_f64()>tolerance*scale {
                    return Err(invalid("original transverse Kelvin nullspace is unresolved"));
                }
            }
        }
        let equilibrium=Root3::new(&model.equilibrium_gram,outside)?;
        let mut map=[[0.0;6];3];
        for j in 0..3 {
            map[j][inside[j]]=ENGINEERING[j];
            let rhs=core::array::from_fn(|i|-Dd::from_f64(model.c_eq[outside[i]][inside[j]])*Dd::from_f64(ENGINEERING[j]));
            let extension=equilibrium.solve(equilibrium.solve(rhs,true)?,false)?;
            for i in 0..3 {map[j][outside[i]]=extension[i].to_f64();}
        }
        let mut w=[[0.0;3];3];let mut dissipative_map=map;
        let mut direct_gram=[[0.0;6];6];let mut direct_work=[[0.0;6];6];
        for j in 0..3 {
            for row in 0..6 {
                let value=(0..6).fold(Dd::ZERO,|s,i|s+Dd::from_f64(source[i][row])*Dd::from_f64(map[j][i]));
                if let Some(index)=columns.iter().position(|&c|c==row) {
                    for i in 0..3 {w[j][i]=(Dd::from_f64(w[j][i])+Dd::from_f64(range[index][i])*value).to_f64();}
                }else{direct_gram[row][j]=number(value)?.to_f64();}
            }
            for i in 0..3 {dissipative_map[j][outside[i]]-=w[j][i];}
        }
        let original_d=model.kv_pa_s().expect("original Kelvin tensor");
        for i in 0..3 {for j in 0..3 {
            let mut sum=Sum::default();for a in 0..6 {for b in 0..6 {sum.add([map[i][a],original_d[a][b],dissipative_map[j][b],1.0]);}}
            direct_work[i][j]=number(sum.retained_value())?.to_f64();
        }}
        let mut rows=vec![[0.0;3];original];let mut force=vec![[0.0;3];original];
        for row in 0..3 {
            for i in 0..3 {rows[row][i]=equilibrium.scale*equilibrium.qr.r(row,i);}
            for j in 0..3 {force[row][j]=-equilibrium.apply(w[j].map(Dd::from_f64))?[row].to_f64();}
        }
        let mut roots=Vec::with_capacity(model.branches().len());
        for (j,branch) in model.branches().iter().enumerate() {
            poll(cancelled)?;let root=checked(1.0/det::sqrt(branch.tau_s),"unrepresentable original Maxwell rate root")?;
            if root<=0.0 {return Err(invalid("original Maxwell rate root underflowed"));}
            roots.push(root);let b=branch.resolved_gram_rows(&model.factors[j]);
            for row in 0..6 {
                let index=3+6*j+row;rows[index]=outside.map(|i|b[row][i]);
                for k in 0..3 {force[index][k]=(0..6).fold(Dd::ZERO,|s,i|s+Dd::from_f64(b[row][i])*Dd::from_f64(dissipative_map[k][i])).to_f64();}
            }
        }
        let mut a=zeros(original*algebraic)?;
        for row in 0..original {for column in 0..algebraic {
            a[row*algebraic+column]=dot3(&rows[row],&null[column].map(Dd::from_f64)).to_f64();
        }}
        let scale=a.iter().fold(0.0_f64,|s,v|s.max(v.abs()));
        if !scale.is_finite()||scale<=0.0 {return Err(invalid("unrepresentable original transverse elastic constraint"));}
        for value in &mut a {*value/=scale;}
        poll(cancelled)?;let elastic=qr(&a,original,algebraic);poll(cancelled)?;pivots(&elastic,algebraic)?;
        let mut compliance=Vec::with_capacity(algebraic);
        for column in 0..algebraic {
            let mut rhs=[Dd::ZERO;3];rhs[column]=Dd::ONE;
            let coefficients=triangular(&elastic,algebraic,rhs,false,scale)?;
            let root=core::array::from_fn(|i|(0..algebraic).fold(Dd::ZERO,|s,j|s+Dd::from_f64(null[j][i])*coefficients[j]).to_f64());
            if root.iter().any(|v|!v.is_finite()) {return Err(invalid("unrepresentable original transverse elastic compliance"));}
            compliance.push(root);
        }
        let mut basis=zeros(original*dynamic)?;let mut factor=zeros(dynamic*dynamic)?;
        for column in 0..dynamic {
            poll(cancelled)?;let mut vector=zeros(original)?;vector[algebraic+column]=1.0;elastic.apply_q(&mut vector);
            let mut spring=[Dd::ZERO;3];
            for row in 0..original {
                basis[row*dynamic+column]=vector[row];
                for i in 0..3 {spring[i]=spring[i]+Dd::from_f64(rows[row][i])*Dd::from_f64(vector[row]);}
            }
            for row in 0..rank {factor[row*dynamic+column]=number(dot3(&range[row],&spring))?.to_f64();}
            for row in 3..original {factor[(rank+row-3)*dynamic+column]=checked(roots[(row-3)/6]*vector[row],"unrepresentable constrained Maxwell factor")?;}
        }
        let prepared=Self {model:Arc::clone(model),plane,inside,outside,equilibrium,viscosity:None,
            constraint:Some(ConstrainedPlaneStress {rank,range,null,compliance,tolerance}),
            map,direct_gram,direct_work,rates:Vec::new(),rate_roots:Vec::new(),
            modes:Vec::new(),forcing:Vec::new(),strain_modes:Vec::new()};
        Self::finish_factor(prepared,factor,&force,dynamic,Some(&basis),cancelled)
    }
}
