//! Geometry-owned modal reduction of a curved, variable-thickness shell.
//!
//! Membrane strain uses the Green--Lagrange metric of the P1 displacement
//! field: epsilon = sym(grad u) + grad(u)^T grad(u)/2. Bending and relative
//! drilling stabilization remain the existing linear DKT shell. Curvature
//! therefore produces quadratic AND cubic restoring forces; this is not a
//! rectangular von Karman tensor grafted onto a cymbal. All membrane terms
//! form a positive strain-energy density. The original linear membrane is
//! replaced, not counted twice. Rigid finite rotations are not exact in the
//! retained linear bending law: this is a moderate-rotation reduced model.
//!
//! Cold work projects geometry once. Evaluation uses rank-three displacement
//! gradient factors in O(facets*modes + modes^2), with no temporary allocation
//! and no dense fourth-order modal tensor. It is NOT a real-time benchmark or
//! a claim of adequate truncation: missing in-plane modes require enrichment.
use super::{ShellMesh, ShellModel, local_bending_stiffness};
use crate::{ModePair, PlateError, PlateSection};

/// Two-sided finite-thickness acoustic boundary and reciprocal modal projection.
pub mod radiation;

/// Explicit cold reduction ceilings. These are not acoustic error estimates.
#[derive(Debug, Clone, Copy)]
pub struct ReductionBudget {
    /// Maximum retained modes; all supplied modes must fit, none are dropped.
    pub max_modes: usize,
    /// Maximum facet-by-mode records evaluated each potential/gradient call.
    pub max_facet_modes: usize,
    /// Relative allowance on mass orthogonality and projection consistency.
    pub relative_tolerance: f64,
}
#[derive(Debug, Clone, Copy)]
struct ModeStrain { dx:[f64;3], dy:[f64;3], linear:[f64;3] }
#[derive(Debug, Clone)]
struct Facet { membrane:[f64;9], strains:Vec<ModeStrain> }

/// Mass-normalized potential and reciprocal geometric force/motion projections.
/// It stores no dynamic state or damping; a coupling/time owner advances those.
#[derive(Debug, Clone)]
pub struct ShellReduction {
    omegas:Vec<f64>,
    remainder:Vec<f64>,
    facets:Vec<Facet>,
    // Mode-major physical translations and axial rotations in the same basis.
    translations:Vec<[f64;3]>,
    rotations:Vec<[f64;3]>,
    reference_positions:Vec<[f64;3]>,
    section_thicknesses:Vec<f64>,
    nodes:usize,
    triangles:Vec<[usize;3]>,
    area_normals:Vec<[f64;3]>,
}
fn bad(what:&'static str)->PlateError {PlateError::BadSection{what}}
// Use the existing shell rigid-motion conformance scale, not the very small
// retained elastic eigenvalues, for cancellation in an assembled K projection.
// This is a floating-point admission allowance, not a physical soft spring.
const RIGID_PROJECTION_ROUNDOFF: f64 = 2.0e-12;

// Recognize only an EXACT supplied translation with zero physical rotations.
// Near-rigid eigenvectors and rigid rotations keep the ordinary spectral gate.
fn exact_translation(mesh:&ShellMesh, model:&ShellModel, mode:&ModePair)->bool {
    if mode.lambda!=0.0 {return false;}
    let at=|node:usize,c:usize|model.dof_map[6*node+c].map_or(0.0,|i|mode.phi[i]);
    let reference:[f64;3]=core::array::from_fn(|c|at(0,c));
    reference.iter().any(|v|*v!=0.0) && (0..mesh.nodes.len()).all(|node|
        (0..3).all(|c|at(node,c)==reference[c]) && (3..6).all(|c|at(node,c)==0.0))
}
fn dot(a:[f64;3],b:[f64;3])->f64 {a[0]*b[0]+a[1]*b[1]+a[2]*b[2]}
fn mul(a:&[f64;9],x:[f64;3])->[f64;3] {
    [a[0]*x[0]+a[1]*x[1]+a[2]*x[2],a[3]*x[0]+a[4]*x[1]+a[5]*x[2],a[6]*x[0]+a[7]*x[1]+a[8]*x[2]]
}
impl ShellReduction {
    /// Reduce the same admitted mesh, sections and pencil used for the supplied
    /// fs-modal modes. Verify M-orthonormality and the full small-signal stiffness
    /// against that pencil, so unrelated shape/frequency lists cannot masquerade
    /// as this geometry. No modal stiffness is guessed from an instrument name.
    ///
    /// # Errors
    /// Dimension, budget, finite-set, basis or section/pencil disagreement.
    pub fn new(mesh:&ShellMesh, sections:&[PlateSection], model:&ShellModel,
        modes:&[ModePair], budget:ReductionBudget)->Result<Self,PlateError> {
        mesh.validate()?;
        let n=modes.len();
        if n==0 || n>budget.max_modes || !budget.relative_tolerance.is_finite()
            || budget.relative_tolerance<=0.0 || budget.relative_tolerance>=1.0
            || n.checked_mul(mesh.tris.len()).is_none_or(|k|k>budget.max_facet_modes) {
            return Err(bad("shell reduction exceeds declared mode/work limits or tolerance is invalid"));
        }
        if sections.len()!=mesh.tris.len() {
            return Err(PlateError::SectionCount{expected:mesh.tris.len(),actual:sections.len()});
        }
        if model.dof_map.len()!=6*mesh.nodes.len() || model.k.nrows()!=model.free
            || model.k.ncols()!=model.free || model.m.nrows()!=model.free || model.m.ncols()!=model.free
            || model.dof_map.iter().flatten().any(|&i|i>=model.free)
            || modes.iter().any(|m|m.phi.len()!=model.free || !m.lambda.is_finite() || m.lambda<0.0
                || m.phi.iter().any(|v|!v.is_finite())) {
            return Err(bad("shell modal basis and pencil dimensions or frequencies disagree"));
        }
        let nn=n.checked_mul(n).ok_or_else(||bad("shell projection size overflow"))?;
        let translation:Vec<_>=modes.iter().map(|m|exact_translation(mesh,model,m)).collect();
        let mut original=vec![0.0;nn];
        let mut roundoff=vec![0.0;nn];
        let mut work=vec![0.0;model.free];
        let mut absolute_work=vec![0.0;model.free];
        for j in 0..n {
            model.m.spmv(&modes[j].phi,&mut work);
            for i in 0..n {
                let mass=modes[i].phi.iter().zip(&work).map(|(a,b)|a*b).sum::<f64>();
                if !mass.is_finite() || (mass-if i==j {1.0}else{0.0}).abs()>budget.relative_tolerance {
                    return Err(bad("shell modal basis is not mass orthonormal"));
                }
            }
            model.k.spmv(&modes[j].phi,&mut work);
            if translation.iter().any(|v|*v) {
                for (row,out) in absolute_work.iter_mut().enumerate() {
                    let (columns,values)=model.k.row(row);
                    *out=columns.iter().zip(values).map(|(&column,value)|
                        value.abs()*modes[j].phi[column].abs()).sum();
                    if !out.is_finite() {return Err(bad("shell rigid projection scale overflow"));}
                    // A geometrically constant displacement must also be a
                    // numerical null vector of the supplied pencil. A grounded
                    // or mismatched model must not get a free mode by declaration.
                    if translation[j] && (!work[row].is_finite()
                        || work[row].abs()>RIGID_PROJECTION_ROUNDOFF*(*out).max(1.0)) {
                        return Err(bad("shell supplied translation is not free in its pencil"));
                    }
                }
            }
            for i in 0..n {
                original[i*n+j]=modes[i].phi.iter().zip(&work).map(|(a,b)|a*b).sum();
                if translation[i] || translation[j] {
                    roundoff[i*n+j]=RIGID_PROJECTION_ROUNDOFF*modes[i].phi.iter()
                        .zip(&absolute_work).map(|(a,b)|a.abs()*b).sum::<f64>();
                    if !roundoff[i*n+j].is_finite() {return Err(bad("shell rigid projection allowance overflow"));}
                }
            }
        }
        // The ordinary spectral criterion is unchanged for elastic modes.
        // Exact translations additionally use |Phi|^T |K| |Phi| to account for
        // cancellation of large assembled entries, not an arbitrary floor on
        // the physical eigenvalues. No spring or frequency is inserted.
        let spectral_scale=modes.iter().map(|m|m.lambda).fold(0.0_f64,f64::max);
        if spectral_scale<=0.0 {return Err(bad("shell basis must retain at least one elastic mode"));}
        let scale_for=|i:usize| modes[i].lambda.max(1e-10*spectral_scale);
        for i in 0..n { for j in 0..n {
            let expected=if i==j { modes[i].lambda } else { 0.0 };
            let scale=scale_for(i).sqrt()*scale_for(j).sqrt();
            if !original[i*n+j].is_finite() || !scale.is_finite()
                || (original[i*n+j]-expected).abs()>budget.relative_tolerance*scale+roundoff[i*n+j] {
                return Err(bad("shell frequencies/shapes are not eigenpairs of the supplied pencil"));
            }
        }}
        let mut translations=Vec::with_capacity(n*mesh.nodes.len());
        let mut rotations=Vec::with_capacity(n*mesh.nodes.len());
        for mode in modes {for node in 0..mesh.nodes.len() {
            translations.push(core::array::from_fn(|c|model.dof_map[6*node+c].map_or(0.0,|i|mode.phi[i])));
            rotations.push(core::array::from_fn(|c|model.dof_map[6*node+3+c].map_or(0.0,|i|mode.phi[i])));
        }}
        let mut remainder=vec![0.0;nn];let mut membrane_linear=vec![0.0;nn];
        let mut facets=Vec::with_capacity(mesh.tris.len());let mut area_normals=Vec::with_capacity(mesh.tris.len());
        for (element,tri) in mesh.tris.iter().enumerate() {
            let section=&sections[element];section.validate()?;
            let g=mesh.facet(element)?;
            let kb=local_bending_stiffness(&g,section,element)?;
            let mut local=vec![[0.0;18];n];let mut strains=Vec::with_capacity(n);
            for (k,mode) in modes.iter().enumerate() {
                let(mut dx,mut dy)=([0.0;3],[0.0;3]);
                let reference=translations[k*mesh.nodes.len()+tri[0]];
                for a in 0..3 {
                    let u=translations[k*mesh.nodes.len()+tri[a]];
                    // Partition of unity: use differences BEFORE derivatives.
                    // This makes every constant translation exactly strain-free,
                    // including large rigid displacement beside elastic motion.
                    let difference=core::array::from_fn(|c|u[c]-reference[c]);
                    for c in 0..3 {dx[c]+=g.gradient[a][0]*difference[c];dy[c]+=g.gradient[a][1]*difference[c];}
                    let rotation=core::array::from_fn(|c|model.dof_map[6*tri[a]+3+c].map_or(0.0,|i|mode.phi[i]));
                    for c in 0..3 {
                        local[k][6*a+c]=dot(g.frame[c],difference);
                        local[k][6*a+3+c]=dot(g.frame[c],rotation);
                    }
                }
                strains.push(ModeStrain{dx,dy,linear:[dot(g.frame[0],dx),dot(g.frame[1],dy),
                    dot(g.frame[0],dy)+dot(g.frame[1],dx)]});
            }
            let membrane=section.d.map(|d|g.area_m2*(12.0/(section.thickness*section.thickness))*d);
            for i in 0..n {for j in 0..=i {
                let mut value=0.0;
                for a in 0..18 {let mut row=0.0;for b in 0..18 {row+=kb[18*a+b]*local[j][b];}value+=local[i][a]*row;}
                remainder[i*n+j]+=value;
                let value=dot(strains[i].linear,mul(&membrane,strains[j].linear));
                membrane_linear[i*n+j]+=value;
                if i!=j {remainder[j*n+i]=remainder[i*n+j];membrane_linear[j*n+i]=membrane_linear[i*n+j];}
            }}
            facets.push(Facet{membrane,strains});
            area_normals.push(g.frame[2].map(|v|v*g.area_m2));
        }
        for i in 0..n {for j in 0..n {
            let actual=remainder[i*n+j]+membrane_linear[i*n+j];
            let scale=scale_for(i).sqrt()*scale_for(j).sqrt();
            if !actual.is_finite() || !scale.is_finite() || scale<=0.0
                || (actual-original[i*n+j]).abs()>budget.relative_tolerance*scale+roundoff[i*n+j] {
                return Err(bad("shell reduction's small-signal stiffness does not match supplied pencil"));
            }
        }}
        Ok(Self{omegas:modes.iter().map(|m|m.lambda.sqrt()).collect(),remainder,facets,translations,
            nodes:mesh.nodes.len(),triangles:mesh.tris.clone(),area_normals,rotations,
            reference_positions:mesh.nodes.clone(),section_thicknesses:sections.iter().map(|s|s.thickness).collect()})
    }
    /// Linear frequencies of the supplied numerical basis [rad/s]; zero denotes
    /// an explicitly supplied rigid coordinate, not an invented soft spring.
    #[must_use]
    pub fn omegas(&self)->&[f64] {&self.omegas}
    /// Number of retained generalized coordinates, each with unit modal mass.
    #[must_use]
    pub fn mode_count(&self)->usize {self.omegas.len()}

    fn strain(f:&Facet,q:&[f64])->([f64;3],[f64;3],[f64;3]) {
        let(mut dx,mut dy,mut strain)=([0.0;3],[0.0;3],[0.0;3]);
        for (&q,s) in q.iter().zip(&f.strains) {for c in 0..3 {
            dx[c]+=s.dx[c]*q;dy[c]+=s.dy[c]*q;strain[c]+=s.linear[c]*q;
        }}
        strain[0]+=0.5*dot(dx,dx);strain[1]+=0.5*dot(dy,dy);strain[2]+=dot(dx,dy);
        (strain,dx,dy)
    }
    /// Mechanical potential [J] at mass-normalized displacements. Nonlinear
    /// membrane geometry couples the modes; no random excitation or pitch curve.
    /// Invalid cardinality/nonfinite input yields NaN for the time owner's
    /// finite-set gate; there is no hidden allocation, clamping or force repair.
    #[must_use]
    pub fn potential(&self,q:&[f64])->f64 {
        let n=self.mode_count();
        if q.len()!=n || q.iter().any(|v|!v.is_finite()) {return f64::NAN;}
        let mut energy=0.0;
        for i in 0..n {let mut row=0.0;for j in 0..n {row+=self.remainder[i*n+j]*q[j];}energy+=0.5*q[i]*row;}
        for f in &self.facets {let(s,_,_)=Self::strain(f,q);energy+=0.5*dot(s,mul(&f.membrane,s));}
        energy
    }
    /// Exact derivative of `potential`, with no heap allocation. A bad slice
    /// shape/nonfinite input fills output with NaN rather than indexing unsafely.
    pub fn gradient(&self,q:&[f64],out:&mut[f64]) {
        let n=self.mode_count();
        if q.len()!=n || out.len()!=n || q.iter().any(|v|!v.is_finite()) {out.fill(f64::NAN);return;}
        for i in 0..n {out[i]=0.0;for j in 0..n {out[i]+=self.remainder[i*n+j]*q[j];}}
        for f in &self.facets {
            let(strain,dx,dy)=Self::strain(f,q);let stress=mul(&f.membrane,strain);
            for (i,s) in f.strains.iter().enumerate() {
                let derivative=[s.linear[0]+dot(dx,s.dx),s.linear[1]+dot(dy,s.dy),
                    s.linear[2]+dot(dy,s.dx)+dot(dx,s.dy)];
                out[i]+=dot(derivative,stress);
            }
        }
    }
    /// Reciprocal point force/motion port at explicit triangle barycentrics.
    /// Generalized force is weights*F and physical velocity is weights dot v.
    /// The direction must be unit length. This is a point idealization, not
    /// an automatically resolved finite stick-tip contact patch.
    /// # Errors
    /// Invalid facet, barycentric location or direction.
    pub fn point_port(&self,triangle:usize,barycentric:[f64;3],direction:[f64;3])->Result<Vec<f64>,PlateError> {
        let tri=self.triangles.get(triangle).ok_or_else(||bad("point port facet out of range"))?;
        if barycentric.iter().any(|x|!x.is_finite() || *x<0.0)
            || (barycentric.iter().sum::<f64>()-1.0).abs()>1e-12
            || direction.iter().any(|v|!v.is_finite()) || (dot(direction,direction)-1.0).abs()>1e-12 {
            return Err(bad("point port needs barycentric coordinates and a finite unit direction"));
        }
        Ok((0..self.mode_count()).map(|k| {
            (0..3).map(|a|barycentric[a]*dot(direction,self.translations[k*self.nodes+tri[a]])).sum()
        }).collect())
    }
    /// Integrated oriented normal-displacement port [m^2/sqrt(kg)]. Multiplying
    /// by pressure gives generalized force; dotting modal velocity gives volume
    /// flow. This is a reciprocal coupling integral, NOT an acoustic radiation
    /// model or permission to treat an unbaffled cymbal as a monopole.
    #[must_use]
    pub fn volume_port(&self)->Vec<f64> {
        (0..self.mode_count()).map(|k|self.triangles.iter().zip(&self.area_normals).map(|(tri,&a)| {
            tri.iter().map(|&node|dot(a,self.translations[k*self.nodes+node])/3.0).sum::<f64>()
        }).sum()).collect()
    }
}
