//! Exact algebraic plane-stress elimination of the COMPLETE Maxwell energy.
//!
//! Let e=(e_parallel,e_normal) and a_j be the viscous Gram coordinates.
//! H=V e^T C_eq e/2 + sum |sqrt(V) B_j e-a_j|^2/2.
//! Zero normal traction minimizes this energy over the three normal strains.
//! It is one coupled Schur complement, NOT one condensation per Prony branch.
//! The retained a_j obey a_j' = -(1/tau_j) dH/da_j. Coupling this energy to
//! physical inertia therefore uses the original fs-phs J/R/Q/G integrator.
//! All six coordinates per branch survive, including zero/dependent rows.
use super::*;
use fs_material::visco_tensor::{TensorProny, TensorProjectionBudget};
use fs_material::Voigt;
use fs_math::det;

/// One homogeneous material quadrature contribution under local z-plane stress.
/// In-plane columns are [epsilon_xx,epsilon_yy,gamma_xy] per displacement [1/m].
/// Normal epsilon_zz,epsilon_yz,epsilon_zx are reconstructed at EVERY state.
/// A fixed material orientation must be applied before construction. Nonzero
/// Kelvin viscosity refuses: its normal strain-rate constraint is a different
/// differential-algebraic problem, not an algebraic Maxwell condensation.
#[derive(Clone, Debug)]
pub struct PlaneStressPoint {
    source: ProjectedTensorProny,
    coordinates: usize,
    energy: Vec<f64>,
    normal_map: Vec<Vec<f64>>,
    pub(super) relaxed_memory: Vec<Vec<f64>>,
    rates: Vec<f64>,
}

// Scaled Cholesky solve of the existing material's 3x3 normal block. No inverse,
// penalty stiffness, symmetrization, nullspace clipping or new numerical solver.
fn normal_solve(matrix:&[f64],d:usize,n:usize,rhs:&[f64])->Result<[f64;3]> {
    let scale:[f64;3]=core::array::from_fn(|i|det::sqrt(matrix[(n+i)*d+n+i]));
    if scale.iter().any(|v|!v.is_finite() || *v<=0.0) {return Err(bad("unresolved plane-stress normal stiffness"));}
    let normalized:Vec<_>=(0..9).map(|k|matrix[(n+k/3)*d+n+k%3]/scale[k/3]/scale[k%3]).collect();
    let factor=fs_la::factor::cholesky(&normalized,3)
        .map_err(|_|bad("plane-stress normal stiffness is not positive definite"))?;
    let mut x:[f64;3]=core::array::from_fn(|i|rhs[i]/scale[i]);
    factor.solve(&mut x);
    for i in 0..3 {
        let mut residual=-rhs[i]/scale[i];let mut magnitude=(rhs[i]/scale[i]).abs();
        for k in 0..3 {residual+=normalized[3*i+k]*x[k];magnitude+=(normalized[3*i+k]*x[k]).abs();}
        if !residual.is_finite() || !magnitude.is_finite()
            || residual.abs()>1e-10*magnitude.max(f64::MIN_POSITIVE) {
            return Err(bad("plane-stress normal solve residual exceeds tolerance"));
        }
    }
    for i in 0..3 {x[i]/=scale[i];}
    if !finite(&x) {return Err(bad("plane-stress normal solve overflow"));}
    Ok(x)
}

impl PlaneStressPoint {
    /// Prepare a work-conjugate plane-stress energy and strain reconstruction.
    /// Positive volume includes the quadrature weight. The projection budget's
    /// coordinate ceiling includes the THREE eliminated normal strain columns;
    /// its work ceiling covers both projection and full energy condensation.
    /// Nothing changes the material's original source, band, branches or domain.
    pub fn new(material:&TensorProny,columns:&[[f64;3]],volume_m3:f64,
        budget:TensorProjectionBudget,gate:&CancelGate)->Result<Self>
    {
        poll(gate)?;
        let n=columns.len();let d=n.checked_add(3).ok_or_else(||bad("plane-stress coordinate extent"))?;
        let memory=material.branches().len().checked_mul(6).ok_or_else(||bad("plane-stress memory extent"))?;
        let m=n.checked_add(memory).ok_or_else(||bad("plane-stress energy extent"))?;
        if n==0 || d>budget.max_coordinates || m>MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || columns.iter().flatten().any(|v|!v.is_finite()) {
            return Err(bad("plane stress requires bounded finite in-plane strain columns and full memory"));
        }
        if material.kv_pa_s().is_some_and(|c|c.iter().flatten().any(|v|*v!=0.0)) {
            return Err(bad("plane-stress Maxwell elimination does not omit Kelvin normal-rate terms"));
        }
        // Include the original strain projection, total K, normal solves, full
        // equilibrium congruence and every retained branch residual Gram product.
        let work=(material.branches().len() as u128+2)*36*d as u128
            +12*(d as u128).pow(2)+memory as u128*(d as u128).pow(2)
            +128*m as u128+(m as u128).pow(2)*((d as u128).pow(2)+memory as u128)
            +memory as u128*m as u128*d as u128;
        if work>budget.max_scalar_terms as u128 {return Err(bad("plane-stress condensation exceeds scalar-work budget"));}
        let mut full:Vec<Voigt>=columns.iter().map(|s|[s[0],s[1],0.0,0.5*s[2],0.0,0.0]).collect();
        full.extend([[0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0]]);
        let source=material.project_strain(&full,volume_m3,budget,||gate.is_requested())
            .map_err(TensorMaterialError::Material)?;
        let rows:Vec<_>=source.branches().iter().flat_map(|b|b.rows()).collect();
        let rates:Vec<_>=source.branches().iter().flat_map(|b|[b.relaxation_time_s().recip();6]).collect();
        if rates.iter().any(|v|!v.is_finite() || *v<=0.0) {return Err(bad("plane-stress relaxation rate is unrepresentable"));}
        let mut instantaneous=source.equilibrium_stiffness().to_vec();
        for row in &rows {
            poll(gate)?;
            for i in 0..d {for j in i..d {
                instantaneous[i*d+j]+=row[i]*row[j];instantaneous[j*d+i]=instantaneous[i*d+j];
            }}
        }
        if !finite(&instantaneous) {return Err(bad("plane-stress instantaneous stiffness overflow"));}
        let mut normal_map=vec![vec![0.0;m];3];
        for col in 0..m {
            poll(gate)?;
            let rhs:[f64;3]=core::array::from_fn(|i|if col<n {-instantaneous[(n+i)*d+col]} else {rows[col-n][n+i]});
            let x=normal_solve(&instantaneous,d,n,&rhs)?;
            for i in 0..3 {normal_map[i][col]=x[i];}
        }
        // e_full = T [u,a]. Compute the equilibrium congruence followed by
        // residual-factor Grams, rather than subtracting two large energies.
        let mut t=vec![vec![0.0;m];d];
        for i in 0..n {t[i][i]=1.0;}
        for i in 0..3 {t[n+i].copy_from_slice(&normal_map[i]);}
        let mut residuals=vec![vec![0.0;m];memory];
        for (r,row) in rows.iter().enumerate() {
            poll(gate)?;
            for col in 0..m {
                residuals[r][col]=(0..d).fold(0.0,|s,k|s+row[k]*t[k][col]);
                if col==n+r {residuals[r][col]-=1.0;}
            }
        }
        let mut energy=vec![0.0;m*m];
        for i in 0..m {
            poll(gate)?;
            for j in i..m {
                let mut value=0.0;
                for a in 0..d {for b in 0..d {value+=t[a][i]*source.equilibrium_stiffness()[a*d+b]*t[b][j];}}
                for row in &residuals {value+=row[i]*row[j];}
                energy[i*m+j]=value;energy[j*m+i]=value;
            }
        }
        let mut relaxed_memory=vec![vec![0.0;n];memory];
        for col in 0..n {
            poll(gate)?;
            let rhs:[f64;3]=core::array::from_fn(|i|-source.equilibrium_stiffness()[(n+i)*d+col]);
            let normal=normal_solve(source.equilibrium_stiffness(),d,n,&rhs)?;
            for (r,row) in rows.iter().enumerate() {
                relaxed_memory[r][col]=row[col]+(0..3).fold(0.0,|s,k|s+row[n+k]*normal[k]);
            }
        }
        if !finite(&energy) || relaxed_memory.iter().flatten().any(|v|!v.is_finite()) {
            return Err(bad("plane-stress energy or initialization overflow"));
        }
        poll(gate)?;
        Ok(Self {source,coordinates:n,energy,normal_map,relaxed_memory,rates})
    }
    /// Original in-plane displacement coordinates.
    pub fn coordinates(&self)->usize {self.coordinates}
    /// All original six-per-branch viscous factor coordinates a_j [sqrt(J)].
    pub fn memory_coordinates(&self)->usize {self.rates.len()}
    /// Complete local energy Q for [displacement, viscous factor coordinates].
    pub fn energy_matrix(&self)->&[f64] {&self.energy}
    /// Original immutable constitutive law; no fitted reduced law is substituted.
    pub fn material(&self)->&TensorProny {self.source.material()}

    fn full_strain_coordinates(&self,u:&[f64],a:&[f64])->Result<Vec<f64>> {
        if u.len()!=self.coordinates || a.len()!=self.memory_coordinates() || !finite(u) || !finite(a) {
            return Err(bad("plane-stress displacement and memory must match their point"));
        }
        let mut full=u.to_vec();
        for row in &self.normal_map {
            full.push(row.iter().zip(u.iter().chain(a)).fold(0.0,|s,(b,x)|s+b*x));
        }
        if !finite(&full) {return Err(bad("plane-stress strain reconstruction overflow"));}
        Ok(full)
    }
    /// Full tensor-shear strain, including the history-dependent normal strains.
    /// Enforces the original 3D small-strain limit at the supplied state.
    pub fn strain(&self,u:&[f64],a:&[f64])->Result<Voigt> {
        self.source.strain(&self.full_strain_coordinates(u,a)?).map_err(TensorMaterialError::Material)
    }
    /// Original 3D normal tractions [sigma_zz,sigma_yz,sigma_zx] [Pa].
    /// Recomputed from full material factors, not set to zero by definition.
    pub fn normal_stress_pa(&self,u:&[f64],a:&[f64])->Result<[f64;3]> {
        let full=self.full_strain_coordinates(u,a)?;
        self.source.strain(&full).map_err(TensorMaterialError::Material)?;
        let n=self.coordinates;let d=n+3;
        let mut force:[f64;3]=core::array::from_fn(|i|(0..d).fold(0.0,|s,k|s+self.source.equilibrium_stiffness()[(n+i)*d+k]*full[k]));
        for (row,&memory) in self.source.branches().iter().flat_map(|b|b.rows()).zip(a) {
            let z=row.iter().zip(&full).fold(0.0,|s,(b,x)|s+b*x)-memory;
            for i in 0..3 {force[i]+=row[n+i]*z;}
        }
        let volume=self.source.volume_m3();
        let stress=[force[0]/volume,0.5*force[1]/volume,0.5*force[2]/volume];
        if !finite(&stress) {return Err(bad("plane-stress traction overflow"));}
        Ok(stress)
    }
}

impl TensorMaterialBody {
    /// Assemble force-driven mechanics with exact local plane-stress constraints.
    /// The SAME accepted-history and exact-ZOH paths as full 3D points are used.
    /// No normal inertia, transverse stress or Kelvin rate term is invented.
    pub fn from_plane_stress(masses_kg:Vec<f64>,points:Vec<PlaneStressPoint>,
        budget:TensorMaterialBudget,gate:&CancelGate)->Result<Self>
    {
        poll(gate)?;
        let states=plane_stress_extent(masses_kg.len(),&points,budget)?;
        let inertia=Inertia::diagonal(masses_kg)?;
        Self::from_plane_stress_inertia(inertia,points,states,gate)
    }

    /// Assemble exact plane-stress memory with a full symmetric positive-definite
    /// physical mass matrix [kg], in row-major displacement-coordinate order.
    /// Momentum is p=M*v; the kinetic energy uses the complete inverse M^-1.
    /// Capacity is admitted before factorization and no off-diagonal inertia is
    /// dropped. Material condensation and force ports match `from_plane_stress`.
    pub fn from_plane_stress_mass_matrix(coordinates:usize,mass_kg:Vec<f64>,
        points:Vec<PlaneStressPoint>,budget:TensorMaterialBudget,gate:&CancelGate)->Result<Self>
    {
        poll(gate)?;
        let states=plane_stress_extent(coordinates,&points,budget)?;
        let inertia=Inertia::consistent(coordinates,mass_kg,gate)?;
        Self::from_plane_stress_inertia(inertia,points,states,gate)
    }

    fn from_plane_stress_inertia(inertia:Inertia,points:Vec<PlaneStressPoint>,
        states:usize,gate:&CancelGate)->Result<Self>
    {
        let n=inertia.dimension();
        let mut q=vec![0.0;states*states];let mut j=vec![0.0;states*states];
        let mut r=vec![0.0;states*states];let mut g=vec![0.0;states*n];
        for i in 0..n {
            for k in 0..n {q[(n+i)*states+n+k]=inertia.inverse_entry(i,k);}
            j[i*states+n+i]=1.0;j[(n+i)*states+i]=-1.0;g[(n+i)*n+i]=1.0;
        }
        let mut offset=2*n;let mut offsets=Vec::with_capacity(points.len());
        for point in &points {
            poll(gate)?;offsets.push(offset);
            let m=n+point.memory_coordinates();
            let address=|i|if i<n {i} else {offset+i-n};
            for i in 0..m {for k in 0..m {q[address(i)*states+address(k)]+=point.energy[i*m+k];}}
            for (i,rate) in point.rates.iter().enumerate() {r[(offset+i)*states+offset+i]=*rate;}
            offset+=point.memory_coordinates();
        }
        if !finite(&q) {return Err(bad("assembled plane-stress energy overflow"));}
        let storage=QuadraticStorage::new(q,states).map_err(TensorMaterialError::Propagation)?;
        let system=PortHamiltonian::new(states,n,j,r,g,Box::new(storage)).map_err(TensorMaterialError::Propagation)?;
        poll(gate)?;
        let split=split_quadratic(&system,LinearSplitBudget {max_modes:states,max_omega_dt:core::f64::consts::PI})
            .map_err(TensorMaterialError::Propagation)?;
        poll(gate)?;
        Ok(Self {system,split,inertia,points:Vec::new(),plane_points:points,plane_offsets:offsets})
    }
    /// Original plane-stress contributions; empty for a full-3D body.
    pub fn plane_stress_points(&self)->&[PlaneStressPoint] {&self.plane_points}
}

// Admit the complete constitutive state before any dense mass factorization.
fn plane_stress_extent(n:usize,points:&[PlaneStressPoint],budget:TensorMaterialBudget)->Result<usize> {
    if n==0 || n>128 || budget.max_state_coordinates==0
        || budget.max_state_coordinates>MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        || budget.max_material_points==0 || budget.max_material_points>4096
        || points.is_empty() || points.len()>budget.max_material_points
        || points.iter().any(|p|p.coordinates()!=n) {
        return Err(bad("bounded plane-stress points and matching inertia coordinates required"));
    }
    let states=points.iter().try_fold(2*n,|s,p|s.checked_add(p.memory_coordinates()))
        .ok_or_else(||bad("plane-stress body extent overflow"))?;
    let block=states.checked_add(n+1).and_then(|v|v.checked_mul(2))
        .ok_or_else(||bad("plane-stress exponential extent overflow"))?;
    if states>budget.max_state_coordinates || block>MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
        return Err(bad("complete plane-stress memory exceeds dense capacity; no coordinates removed"));
    }
    Ok(states)
}

impl TensorMaterialHistory<'_> {
    /// Reconstruct one accepted quadrature strain without exposing mutable memory.
    pub fn plane_stress_strain(&self,point:usize)->Result<Voigt> {
        let p=self.body.plane_points.get(point).ok_or_else(||bad("unknown plane-stress point"))?;
        let offset=self.body.plane_offsets[point];
        p.strain(self.displacement_m(),&self.state[offset..offset+p.memory_coordinates()])
    }
    /// Recompute accepted normal tractions from the original constitutive factors.
    pub fn plane_stress_normal_stress_pa(&self,point:usize)->Result<[f64;3]> {
        let p=self.body.plane_points.get(point).ok_or_else(||bad("unknown plane-stress point"))?;
        let offset=self.body.plane_offsets[point];
        p.normal_stress_pa(self.displacement_m(),&self.state[offset..offset+p.memory_coordinates()])
    }
}

#[cfg(test)]
mod tests;
