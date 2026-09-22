//! Conservative, direction-updating tethers at shell point ports.
//!
//! A point moves by sum_a w_a (u_a + theta_a cross arm). Each tether pulls
//! toward a fixed Cartesian anchor. T = max(0, T_ref + k (length-length_ref));
//! k=0 explicitly selects prescribed tension, not a fixed-rest-length wire.
//! The force and BOTH material/geometric tangent terms derive from one radial
//! potential. No compression, follower-force approximation, or added mass.
use super::{ShellMesh, ShellModel, dot, norm};
use fs_sparse::Coo;

/// One massless tensile connection to a fixed anchor. The three nodes and
/// weights identify an existing shell facet; no new mesh or support is inferred.
#[derive(Clone, Debug)]
pub struct ShellTether {
    /// Nodes of a shell facet, in any order.
    pub nodes: [usize;3],
    /// Nonnegative barycentric weights, summing to one.
    pub weights: [f64;3],
    /// Reference bearing offset [m]; uses the shell's linear rotation map.
    pub arm_m: [f64;3],
    /// Fixed anchor position in the shell Cartesian frame [m].
    pub anchor_m: [f64;3],
    /// Tension at the undeformed reference [N], or target tension when k=0.
    pub reference_tension_n: f64,
    /// dT/d(length) on the taut branch [N/m]. Zero means tension-controlled.
    pub axial_stiffness_n_per_m: f64,
}

/// Final physical tether state, in the same order as the supplied connections.
#[derive(Clone, Copy, Debug)]
pub struct TetherResponse {
    /// Force ON the shell bearing [N]; its moment is arm cross force.
    pub force_n: [f64;3],
    /// Nonnegative current tension [N].
    pub tension_n: f64,
    /// Current anchor-to-bearing length [m].
    pub length_m: f64,
    /// Current minus reference length [m], evaluated without cancellation.
    pub length_change_m: f64,
    /// Potential relative to the reference [J]; can be negative on unloading.
    pub potential_change_j: f64,
}

pub(super) struct PreparedTether {
    dofs: [Option<usize>;18],
    columns: [[f64;3];18],
    chord: [f64;3],
    length: f64,
    tension: f64,
    stiffness: f64,
}
struct State {
    response: TetherResponse,
    gradient: [f64;3],
    hessian: [f64;9],
}
fn cross(a:[f64;3],b:[f64;3])->[f64;3] {
    [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
}
impl PreparedTether {
    pub(super) fn new(t:&ShellTether,mesh:&ShellMesh,model:&ShellModel)->Result<Self,String> {
        if t.nodes.iter().any(|n|*n>=mesh.nodes.len())
            || t.nodes[0]==t.nodes[1] || t.nodes[0]==t.nodes[2] || t.nodes[1]==t.nodes[2]
            || t.weights.iter().any(|w|!w.is_finite() || *w<0. || *w>1.)
            || (t.weights.iter().sum::<f64>()-1.).abs()>1e-12
            || t.arm_m.iter().chain(&t.anchor_m).any(|v|!v.is_finite())
            || !t.reference_tension_n.is_finite() || t.reference_tension_n<0.
            || !t.axial_stiffness_n_per_m.is_finite() || t.axial_stiffness_n_per_m<0.
            || !mesh.tris.iter().any(|tri|t.nodes.iter().all(|n|tri.contains(n))) {
            return Err("invalid shell tether facet, barycentric port, anchor or tension law".into());
        }
        let chord=std::array::from_fn(|c|t.arm_m[c]-t.anchor_m[c]
            +(0..3).map(|a|t.weights[a]*mesh.nodes[t.nodes[a]][c]).sum::<f64>());
        let length=norm(&chord);
        if !length.is_finite() || length<1e-9 {return Err("tether anchor coincides with its bearing".into());}
        let columns=std::array::from_fn(|i| {
            let mut axis=[0.;3];axis[i%3]=1.;
            let column=if i%6<3 {axis}else{cross(axis,t.arm_m)};
            column.map(|x|x*t.weights[i/6])
        });
        Ok(Self {dofs:std::array::from_fn(|i|model.dof_map[6*t.nodes[i/6]+i%6]),
            columns,chord,length,tension:t.reference_tension_n,stiffness:t.axial_stiffness_n_per_m})
    }
    fn movement(&self,u:&[f64])->[f64;3] {
        std::array::from_fn(|c|(0..18).map(|i|
            self.dofs[i].map_or(0.,|d|self.columns[i][c]*u[d])).sum())
    }
    fn state(&self,du:[f64;3])->Result<State,String> {
        let chord: [f64;3]=std::array::from_fn(|i|self.chord[i]+du[i]);
        let length=norm(&chord);
        if !length.is_finite() || length<1e-9 {return Err("tether trial reaches its anchor".into());}
        // Rationalize length-length_ref: exact zero at the reference and no
        // loss of a small extension in the subtraction of metre-scale lengths.
        let change=(2.*dot(self.chord,du)+dot(du,du))/(length+self.length);
        let trial=self.tension+self.stiffness*change;
        if !trial.is_finite() || !change.is_finite() {return Err("tether tension overflow".into());}
        let tension=trial.max(0.);
        let axial=if trial>0. {self.stiffness}else{0.};
        let energy=if trial<0. && self.stiffness>0. {
            -0.5*self.tension*(self.tension/self.stiffness)
        }else{change*(self.tension+0.5*self.stiffness*change)};
        let direction=chord.map(|x|x/length);
        let gradient=direction.map(|x|tension*x);
        let geometric=tension/length;
        let hessian=std::array::from_fn(|i| {
            let (r,c)=(i/3,i%3);
            (if r==c {geometric}else{0.})+(axial-geometric)*direction[r]*direction[c]
        });
        if !energy.is_finite() || gradient.iter().chain(&hessian).any(|v|!v.is_finite()) {
            return Err("tether energy or tangent overflow".into());
        }
        Ok(State {response:TetherResponse {force_n:gradient.map(|x|-x),tension_n:tension,
            length_m:length,length_change_m:change,potential_change_j:energy},gradient,hessian})
    }
    pub(super) fn response(&self,u:&[f64])->Result<TetherResponse,String> {
        Ok(self.state(self.movement(u))?.response)
    }
    pub(super) fn add_gradient(&self,u:&[f64],out:&mut[f64],factor:f64)->Result<(),String> {
        let state=self.state(self.movement(u))?;
        for i in 0..18 {if let Some(d)=self.dofs[i] {
            out[d]+=factor*dot(self.columns[i],state.gradient);
        }}
        Ok(())
    }
    pub(super) fn add_action(&self,u:&[f64],v:&[f64],out:&mut[f64],factor:f64)->Result<(),String> {
        let state=self.state(self.movement(u))?;
        let movement=self.movement(v);
        let hv=super::mul(&state.hessian,movement);
        for i in 0..18 {if let Some(d)=self.dofs[i] {
            out[d]+=factor*dot(self.columns[i],hv);
        }}
        Ok(())
    }
    pub(super) fn add_tangent(&self,u:&[f64],out:&mut Coo,factor:f64)->Result<(),String> {
        let state=self.state(self.movement(u))?;
        for i in 0..18 {if let Some(a)=self.dofs[i] {
            for j in 0..=i {if let Some(b)=self.dofs[j] {
                let value=factor*dot(self.columns[i],super::mul(&state.hessian,self.columns[j]));
                out.push(a,b,value);if a!=b {out.push(b,a,value);}
            }}
        }}
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn radial(stiffness:f64)->PreparedTether {
        let chord=[0.7,-0.2,0.3];
        PreparedTether {dofs:[None;18],columns:[[0.;3];18],chord,length:norm(&chord),
            tension:200.,stiffness}
    }
    #[test]
    fn energy_force_and_material_geometric_tangent_are_consistent() {
        for stiffness in [0.,30_000.] {
            let p=radial(stiffness);let u=[0.001,-0.002,0.003];
            let state=p.state(u).unwrap();let h=1e-6;
            for i in 0..3 {
                let mut a=u;let mut b=u;a[i]+=h;b[i]-=h;
                let a=p.state(a).unwrap();let b=p.state(b).unwrap();
                let fd=(a.response.potential_change_j-b.response.potential_change_j)/(2.*h);
                assert!((fd-state.gradient[i]).abs()<1e-6*(1.+state.gradient[i].abs()));
                for j in 0..3 {
                    let fd=(a.gradient[j]-b.gradient[j])/(2.*h);
                    assert!((fd-state.hessian[3*j+i]).abs()<1e-6*(1.+fd.abs()));
                }
            }
            assert!(state.response.force_n[0]<0.);
            assert!(state.response.length_change_m>0.);
            if stiffness==0. {assert_eq!(state.response.tension_n,200.);}
            else {assert!(state.response.tension_n>200.);}
        }
    }
    #[test]
    fn slack_tether_never_pushes_and_reference_energy_is_exact_zero() {
        let mut p=radial(10_000.);p.chord=[1.,0.,0.];p.length=1.;p.tension=10.;
        assert_eq!(p.state([0.;3]).unwrap().response.potential_change_j,0.);
        let slack=p.state([-0.01,0.,0.]).unwrap();
        assert_eq!(slack.response.tension_n,0.);assert_eq!(slack.gradient,[0.;3]);
        assert_eq!(slack.hessian,[0.;9]);
        assert!((slack.response.potential_change_j+0.005).abs()<1e-15);
        assert!(p.state([-1.,0.,0.]).is_err());
        assert!(p.state([f64::NAN,0.,0.]).is_err());
    }
    #[test]
    fn offset_port_force_action_and_sparse_tangent_use_the_same_map() {
        let mesh=ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],vec![[0,1,2]]).unwrap();
        let model=ShellModel {k:Coo::new(18,18).assemble(),m:Coo::new(18,18).assemble(),
            dof_map:(0..18).map(Some).collect(),free:18};
        let tether=ShellTether {nodes:[0,1,2],weights:[0.2,0.3,0.5],arm_m:[0.02,-0.03,0.04],
            anchor_m:[-0.4,0.1,-0.2],reference_tension_n:200.,axial_stiffness_n_per_m:30_000.};
        let p=PreparedTether::new(&tether,&mesh,&model).unwrap();
        let u:Vec<_>=(0..18).map(|i|(i as f64*0.71).sin()*0.001).collect();
        let v:Vec<_>=(0..18).map(|i|(i as f64*0.23).cos()).collect();
        let mut g=vec![0.;18];p.add_gradient(&u,&mut g,1.).unwrap();
        let mut k=Coo::new(18,18);p.add_tangent(&u,&mut k,1.).unwrap();
        let k=k.assemble();let mut kv=vec![0.;18];k.spmv(&v,&mut kv);
        let mut action=vec![0.;18];p.add_action(&u,&v,&mut action,1.).unwrap();
        for (a,b) in action.iter().zip(&kv) {assert!((a-b).abs()<1e-9*(1.+a.abs()));}
        let h=1e-7;
        let a:Vec<_>=u.iter().zip(&v).map(|(u,v)|u+h*v).collect();
        let b:Vec<_>=u.iter().zip(&v).map(|(u,v)|u-h*v).collect();
        let work:f64=g.iter().zip(&v).map(|(a,b)|a*b).sum();
        let fd=(p.response(&a).unwrap().potential_change_j-p.response(&b).unwrap().potential_change_j)/(2.*h);
        assert!((fd-work).abs()<1e-6*(1.+work.abs()));
        let mut invalid=tether.clone();invalid.weights=[0.;3];
        assert!(PreparedTether::new(&invalid,&mesh,&model).is_err());
        invalid=tether;invalid.nodes[2]=9;
        assert!(PreparedTether::new(&invalid,&mesh,&model).is_err());
    }
}
