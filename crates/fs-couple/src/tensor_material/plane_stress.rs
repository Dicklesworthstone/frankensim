//! Exact plane-stress evolution of the COMPLETE Kelvin-Maxwell material.
//!
//! Let e=(e_parallel,e_normal) and a_j be the viscous Gram coordinates.
//! H=V e^T C_eq e/2 + sum |sqrt(V) B_j e-a_j|^2/2.
//! Without normal viscosity, zero normal traction minimizes this energy over
//! the three normal strains in one coupled Schur complement. Damped normal
//! directions remain massless dynamic coordinates; exactly undamped coordinate
//! directions are condensed through the same full material energy.
//! The retained a_j obey a_j' = -(1/tau_j) dH/da_j. Coupling this energy to
//! physical inertia therefore uses the original fs-phs J/R/Q/G integrator.
//! All six coordinates per branch survive, including zero/dependent rows.
use super::*;
use fs_material::visco_tensor::{TensorProny, TensorProjectionBudget};
use fs_material::Voigt;
use fs_math::det;

mod relaxation;
pub use relaxation::PlaneStressRelaxation;
mod descriptor;
pub use descriptor::{DescriptorPlaneStressBody, DescriptorPlaneStressHistory, DescriptorPlaneStressPoint};

/// One homogeneous material quadrature contribution under local z-plane stress.
/// In-plane columns are [epsilon_xx,epsilon_yy,gamma_xy] per displacement [1/m].
/// Normal epsilon_zz,epsilon_yz,epsilon_zx are reconstructed at EVERY state.
/// A fixed material orientation must be applied before construction. Kelvin
/// viscosity may act on any subset of the three local normal directions.
/// Eliminated directions require exactly zero Kelvin rows; the retained normal
/// block must be positive definite. Unresolved oblique singular blocks refuse.
#[derive(Clone, Debug)]
pub struct PlaneStressPoint {
    source: ProjectedTensorProny,
    coordinates: usize,
    energy: Vec<f64>,
    normal_map: Vec<Vec<f64>>,
    pub(super) relaxed_memory: Vec<Vec<f64>>,
    rates: Vec<f64>,
    viscosity: Vec<f64>,
    normal_viscosity: Option<NormalViscosity>,
}

#[derive(Clone, Debug)]
struct NormalViscosity {
    /// Retained local normal coordinates, in [zz,yz,zx] order.
    indices: Vec<usize>,
    /// D_ss^-1 on the retained normal coordinates.
    mobility: Vec<f64>,
    /// D_ss^-1 D_su, row-major (retained normals x displacement coordinates).
    coupling: Vec<f64>,
}

// Only exact structural zero rows may be eliminated. No threshold decides
// material rank, so arbitrarily small resolved positive viscosity survives.
fn retained_normals(material:&TensorProny)->Vec<usize> {
    material.kv_pa_s().map_or_else(Vec::new,|d|[2,4,5].iter().enumerate()
        .filter_map(|(k,&i)|(0..6).any(|j|d[i][j]!=0.0).then_some(k)).collect())
}

pub(super) fn internal_coordinate_count(material:&TensorProny)->Result<usize> {
    material.branches().len().checked_mul(6).and_then(|n|n.checked_add(retained_normals(material).len()))
        .ok_or_else(||bad("plane-stress internal-coordinate extent"))
}

// Scaled Cholesky solve of the existing material's normal principal block. No
// penalty stiffness, symmetrization, nullspace clipping or new numerical solver.
fn normal_solve(matrix:&[f64],d:usize,n:usize,rhs:&[f64])->Result<[f64;3]> {
    let x=normal_subset_solve(matrix,d,n,&[0,1,2],rhs)?;
    Ok([x[0],x[1],x[2]])
}

fn normal_subset_solve(matrix:&[f64],d:usize,n:usize,indices:&[usize],rhs:&[f64])->Result<Vec<f64>> {
    let count=indices.len();
    let scale:Vec<_>=indices.iter().map(|&i|det::sqrt(matrix[(n+i)*d+n+i])).collect();
    if scale.iter().any(|v|!v.is_finite() || *v<=0.0) {return Err(bad("unresolved plane-stress normal stiffness"));}
    let normalized:Vec<_>=(0..count*count).map(|k|
        matrix[(n+indices[k/count])*d+n+indices[k%count]]/scale[k/count]/scale[k%count]).collect();
    let factor=fs_la::factor::cholesky(&normalized,count)
        .map_err(|_|bad("plane-stress normal stiffness is not positive definite"))?;
    let mut x:Vec<_>=(0..count).map(|i|rhs[i]/scale[i]).collect();
    factor.solve(&mut x);
    for i in 0..count {
        let mut residual=-rhs[i]/scale[i];let mut magnitude=(rhs[i]/scale[i]).abs();
        for k in 0..count {residual+=normalized[count*i+k]*x[k];magnitude+=(normalized[count*i+k]*x[k]).abs();}
        if !residual.is_finite() || !magnitude.is_finite()
            || residual.abs()>1e-10*magnitude.max(f64::MIN_POSITIVE) {
            return Err(bad("plane-stress normal solve residual exceeds tolerance"));
        }
    }
    for i in 0..count {x[i]/=scale[i];}
    if !finite(&x) {return Err(bad("plane-stress normal solve overflow"));}
    Ok(x)
}

impl PlaneStressPoint {
    /// Prepare work-conjugate plane-stress energy, rates and reconstruction.
    /// Positive volume includes the quadrature weight. The projection budget's
    /// coordinate ceiling includes the THREE local normal strain columns;
    /// its work ceiling covers both projection and full energy condensation.
    /// Nothing changes the material's original source, band, branches or domain.
    pub fn new(material:&TensorProny,columns:&[[f64;3]],volume_m3:f64,
        budget:TensorProjectionBudget,gate:&CancelGate)->Result<Self>
    {
        poll(gate)?;
        let n=columns.len();let d=n.checked_add(3).ok_or_else(||bad("plane-stress coordinate extent"))?;
        let memory=material.branches().len().checked_mul(6).ok_or_else(||bad("plane-stress memory extent"))?;
        let internal=internal_coordinate_count(material)?;
        let normal_indices=retained_normals(material);
        let normal_count=normal_indices.len();
        let algebraic_indices:Vec<_>=(0..3).filter(|i|!normal_indices.contains(i)).collect();
        let m=n.checked_add(internal).ok_or_else(||bad("plane-stress energy extent"))?;
        if n==0 || d>budget.max_coordinates || m>MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || columns.iter().flatten().any(|v|!v.is_finite()) {
            return Err(bad("plane stress requires bounded finite in-plane strain columns and full memory"));
        }
        // Include the original strain projection, total K, normal solves, full
        // equilibrium congruence and every retained branch residual Gram product.
        let work=(material.branches().len() as u128+2)*36*d as u128
            +12*(d as u128).pow(2)+memory as u128*(d as u128).pow(2)
            +128*m as u128+(m as u128).pow(2)*((d as u128).pow(2)+memory as u128)
            +memory as u128*m as u128*d as u128
            +(material.branches().len() as u128+2)*36*n as u128+24*(d as u128).pow(2);
        if work>budget.max_scalar_terms as u128 {return Err(bad("plane-stress condensation exceeds scalar-work budget"));}
        let mut full:Vec<Voigt>=columns.iter().map(|s|[s[0],s[1],0.0,0.5*s[2],0.0,0.0]).collect();
        full.extend([[0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0]]);
        let source=material.project_strain(&full,volume_m3,budget,||gate.is_requested())
            .map_err(TensorMaterialError::Material)?;
        let damping=source.viscosity_matrix();
        // The original structural zero must survive projection exactly.
        // Lost nonzero rate contributions are never condensed away.
        if algebraic_indices.iter().any(|&i|(0..d).any(|j|damping[(n+i)*d+j]!=0.0)) {
            return Err(bad("undamped normal constraint has nonzero projected Kelvin viscosity"));
        }
        let mut viscosity=vec![0.0;n*n];
        let normal_viscosity=if normal_count>0 {
            let mut mobility=vec![0.0;normal_count*normal_count];let mut coupling=vec![0.0;normal_count*n];
            for col in 0..normal_count {
                poll(gate)?;
                let mut rhs=vec![0.0;normal_count];rhs[col]=1.0;
                let x=normal_subset_solve(damping,d,n,&normal_indices,&rhs)
                    .map_err(|_|bad("retained normal Kelvin block must be positive definite and resolved; oblique singular constraints are unsupported"))?;
                for i in 0..=col {mobility[i*normal_count+col]=x[i];mobility[col*normal_count+i]=x[i];}
            }
            for i in 0..normal_count {for k in 0..normal_count {
                let value=(0..normal_count).fold(0.0,|s,j|
                    s+damping[(n+normal_indices[i])*d+n+normal_indices[j]]*mobility[j*normal_count+k]);
                if !value.is_finite() || (value-if i==k {1.0} else {0.0}).abs()>192.0*f64::EPSILON {
                    return Err(bad("normal Kelvin inverse does not resolve the original viscosity"));
                }
            }}
            for col in 0..n {
                poll(gate)?;
                let rhs:Vec<_>=normal_indices.iter().map(|&i|damping[(n+i)*d+col]).collect();
                let x=normal_subset_solve(damping,d,n,&normal_indices,&rhs)?;
                for i in 0..normal_count {coupling[i*n+col]=x[i];}
            }
            // D_uu-D_un D_nn^-1 D_nu is evaluated as an original-factor
            // congruence by the material owner, avoiding subtraction of losses.
            let rate_columns:Vec<Voigt>=(0..n).map(|col|core::array::from_fn(|i|
                full[col][i]-(0..normal_count).fold(0.0,|s,k|
                    s+full[n+normal_indices[k]][i]*coupling[k*n+col]))).collect();
            let reduced=material.project_strain(&rate_columns,volume_m3,budget,||gate.is_requested())
                .map_err(TensorMaterialError::Material)?;
            viscosity.copy_from_slice(reduced.viscosity_matrix());
            Some(NormalViscosity {indices:normal_indices.clone(),mobility,coupling})
        } else {
            for i in 0..n {for j in 0..n {viscosity[i*n+j]=damping[i*d+j];}}
            None
        };
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
        for (k,&i) in normal_indices.iter().enumerate() {normal_map[i][n+memory+k]=1.0;}
        if !algebraic_indices.is_empty() {
            for col in 0..m {
                poll(gate)?;
                let rhs:Vec<_>=algebraic_indices.iter().map(|&i|if col<n {
                    -instantaneous[(n+i)*d+col]
                } else if col<n+memory {rows[col-n][n+i]} else {
                    -instantaneous[(n+i)*d+n+normal_indices[col-n-memory]]
                }).collect();
                let x=normal_subset_solve(&instantaneous,d,n,&algebraic_indices,&rhs)?;
                for (k,&i) in algebraic_indices.iter().enumerate() {normal_map[i][col]=x[k];}
            }
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
        let mut relaxed_memory=vec![vec![0.0;n];internal];
        for col in 0..n {
            poll(gate)?;
            let rhs:[f64;3]=core::array::from_fn(|i|-source.equilibrium_stiffness()[(n+i)*d+col]);
            let normal=normal_solve(source.equilibrium_stiffness(),d,n,&rhs)?;
            for (r,row) in rows.iter().enumerate() {
                relaxed_memory[r][col]=row[col]+(0..3).fold(0.0,|s,k|s+row[n+k]*normal[k]);
            }
            for (k,&i) in normal_indices.iter().enumerate() {relaxed_memory[memory+k][col]=normal[i];}
        }
        if !finite(&energy) || relaxed_memory.iter().flatten().any(|v|!v.is_finite()) {
            return Err(bad("plane-stress energy or initialization overflow"));
        }
        poll(gate)?;
        Ok(Self {source,coordinates:n,energy,normal_map,relaxed_memory,rates,viscosity,normal_viscosity})
    }
    /// Original in-plane displacement coordinates.
    pub fn coordinates(&self)->usize {self.coordinates}
    /// All original six-per-branch viscous coordinates a_j [sqrt(J)], followed
    /// by retained normal tensor strains, in [zz,yz,zx] order with exactly
    /// undamped directions omitted. No original Maxwell coordinate is removed.
    pub fn memory_coordinates(&self)->usize {self.rates.len()+self.normal_strain_coordinates()}
    /// Number of massless normal-strain history coordinates, from zero to three.
    pub fn normal_strain_coordinates(&self)->usize {self.normal_viscosity.as_ref().map_or(0,|v|v.indices.len())}
    /// Complete local energy Q for [displacement, material/history coordinates].
    pub fn energy_matrix(&self)->&[f64] {&self.energy}
    /// Original immutable constitutive law; no fitted reduced law is substituted.
    pub fn material(&self)->&TensorProny {self.source.material()}

    pub(super) fn admit_initial_memory(&self,u:&[f64],memory:InitialMemory)->Result<()> {
        if memory==InitialMemory::Instantaneous && self.normal_viscosity.is_some()
            && u.iter().any(|v|*v!=0.0) {
            return Err(bad("normal Kelvin viscosity requires relaxed initial displacement; an instantaneous jump has no finite-work prehistory"));
        }
        Ok(())
    }

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
    /// Normal Kelvin viscosity requires the velocity-aware variant.
    pub fn normal_stress_pa(&self,u:&[f64],a:&[f64])->Result<[f64;3]> {
        if self.normal_viscosity.is_some() {
            return Err(bad("normal Kelvin traction requires the actual in-plane velocity"));
        }
        self.normal_stress_pa_with_velocity(u,&vec![0.0;self.coordinates],a)
    }

    /// Recompute total elastic, Maxwell and Kelvin normal traction [Pa].
    /// Retained rates solve D_ss e_s'=-dH/de_s-D_su u'; remaining tractions
    /// vanish through the algebraic energy constraint, including branch memory.
    /// This exposes the original constitutive traction, not a literal zero.
    pub fn normal_stress_pa_with_velocity(&self,u:&[f64],velocity:&[f64],a:&[f64])->Result<[f64;3]> {
        if velocity.len()!=self.coordinates || !finite(velocity) {
            return Err(bad("plane-stress traction requires matching finite velocity"));
        }
        let full=self.full_strain_coordinates(u,a)?;
        self.source.strain(&full).map_err(TensorMaterialError::Material)?;
        let n=self.coordinates;let d=n+3;
        let mut force:[f64;3]=core::array::from_fn(|i|(0..d).fold(0.0,|s,k|s+self.source.equilibrium_stiffness()[(n+i)*d+k]*full[k]));
        for (row,&memory) in self.source.branches().iter().flat_map(|b|b.rows()).zip(a) {
            let z=row.iter().zip(&full).fold(0.0,|s,(b,x)|s+b*x)-memory;
            for i in 0..3 {force[i]+=row[n+i]*z;}
        }
        if let Some(normal)=&self.normal_viscosity {
            let damping=self.source.viscosity_matrix();
            let rhs:Vec<_>=normal.indices.iter().map(|&i|-force[i]
                -(0..n).fold(0.0,|s,j|s+damping[(n+i)*d+j]*velocity[j])).collect();
            let normal_rate=normal_subset_solve(damping,d,n,&normal.indices,&rhs)?;
            for i in 0..3 {
                force[i]+=(0..n).fold(0.0,|s,j|s+damping[(n+i)*d+j]*velocity[j]);
                force[i]+=(0..normal.indices.len()).fold(0.0,|s,j|
                    s+damping[(n+i)*d+n+normal.indices[j]]*normal_rate[j]);
            }
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
            for i in 0..n {for k in 0..n {r[(n+i)*states+n+k]+=point.viscosity[i*n+k];}}
            if let Some(normal)=&point.normal_viscosity {
                let z=offset+point.rates.len();
                // p'=-H_u-S*v+D_un*D_nn^-1*H_n;
                // e_n'=-D_nn^-1*D_nu*v-D_nn^-1*H_n.
                // The cross coupling is skew, so only the two original
                // Kelvin loss terms and Maxwell dashpots dissipate energy.
                let count=normal.indices.len();
                for i in 0..count {
                    for k in 0..count {r[(z+i)*states+z+k]=normal.mobility[count*i+k];}
                    for k in 0..n {
                        j[(n+k)*states+z+i]=normal.coupling[i*n+k];
                        j[(z+i)*states+n+k]=-normal.coupling[i*n+k];
                    }
                }
            }
            offset+=point.memory_coordinates();
        }
        if !finite(&q) || !finite(&j) || !finite(&r) {return Err(bad("assembled plane-stress operator overflow"));}
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
        p.normal_stress_pa_with_velocity(self.displacement_m(),&self.velocity_m_s(),
            &self.state[offset..offset+p.memory_coordinates()])
    }
}

#[cfg(test)]
mod tests;
