//! Work-conjugate frame changes and harmonic plane-stress elimination.
//! Eliminate the TOTAL complex modulus, never each Prony branch separately:
//! Schur complements do not commute with adding relaxation branches.
use super::*;
use fs_math::c64::C64;

/// Complex in-plane matrix, in the explicitly documented coordinate convention.
pub type Complex3 = [[C64;3];3];
const IN_PLANE: [usize;3] = [0,1,3];
const OUT_OF_PLANE: [usize;3] = [2,4,5];
const ENGINEERING_TO_MANDEL: [f64;3] = [1.0,1.0,core::f64::consts::FRAC_1_SQRT_2];
fn cfinite(v:C64)->bool {v.re.is_finite() && v.im.is_finite()}
fn complex_apply(a:&Complex3,x:&[C64;3])->[C64;3] {
    core::array::from_fn(|i|(0..3).fold(C64::ZERO,|sum,j|sum+a[i][j]*x[j]))
}

// Three RHS columns, scaled together with their matrix before elimination.
// The common scale avoids squaring dimensional stiffness magnitudes in complex
// division. Partial-pivot ties retain the lowest row. No matrix is repaired.
fn solve3(mut a:Complex3,mut b:Complex3)->Result<Complex3> {
    let scale=a.iter().flatten().chain(b.iter().flatten())
        .fold(0.0_f64,|m,z|m.max(z.re.abs()).max(z.im.abs()));
    if !scale.is_finite() || scale<=0.0 {return Err(invalid("unresolved plane-stress elimination scale"));}
    for row in a.iter_mut().chain(b.iter_mut()) {for v in row {
        *v=C64::new(v.re/scale,v.im/scale);
    }}
    for k in 0..3 {
        let mut pivot=k;
        for i in k+1..3 {if a[i][k].abs()>a[pivot][k].abs() {pivot=i;}}
        if !cfinite(a[pivot][k]) || a[pivot][k].abs()==0.0 {
            return Err(invalid("unresolved plane-stress pivot"));
        }
        a.swap(k,pivot);b.swap(k,pivot);
        for i in k+1..3 {
            let factor=a[i][k]/a[k][k];a[i][k]=C64::ZERO;
            for j in k+1..3 {a[i][j]=a[i][j]-factor*a[k][j];}
            for j in 0..3 {b[i][j]=b[i][j]-factor*b[k][j];}
        }
    }
    let mut x=[[C64::ZERO;3];3];
    for i in (0..3).rev() {for col in 0..3 {
        let mut value=b[i][col];
        for j in i+1..3 {value=value-a[i][j]*x[j][col];}
        x[i][col]=value/a[i][i];
        if !cfinite(x[i][col]) {return Err(invalid("nonfinite plane-stress elimination"));}
    }}
    Ok(x)
}

impl TensorProny {
    /// Harmonic plane stress normal to local z, with exp(+i omega t).
    /// Eliminate zz, yz, zx strain from the whole frequency-dependent tensor.
    /// This is NOT a time-domain Prony fit or branch-by-branch condensation.
    pub fn plane_stress(&self,omega:f64)->Result<PlaneStressModulus> {
        let (storage,loss)=self.modulus_mandel(omega)?;
        let c: [[C64;6];6]=core::array::from_fn(|i|core::array::from_fn(|j|
            C64::new(storage[i][j],loss[i][j])));
        let normal=core::array::from_fn(|i|core::array::from_fn(|j|c[OUT_OF_PLANE[i]][OUT_OF_PLANE[j]]));
        let coupling=core::array::from_fn(|i|core::array::from_fn(|j|c[OUT_OF_PLANE[i]][IN_PLANE[j]]));
        let solve=solve3(normal,coupling)?;
        let mut reduced=[[C64::ZERO;3];3];
        for i in 0..3 {for j in 0..3 {
            let mut value=c[IN_PLANE[i]][IN_PLANE[j]];
            for k in 0..3 {value=value-c[IN_PLANE[i]][OUT_OF_PLANE[k]]*solve[k][j];}
            if !cfinite(value) {return Err(invalid("nonfinite plane-stress modulus"));}
            reduced[i][j]=value;
        }}
        Ok(PlaneStressModulus {mandel:reduced,normal_extension:solve,
            strain_limit:self.strain_limit})
    }
}

/// Condensed harmonic material with the full complex out-of-plane strain map.
#[derive(Clone,Debug,PartialEq)]
pub struct PlaneStressModulus {mandel:Complex3,normal_extension:Complex3,strain_limit:f64}
impl PlaneStressModulus {
    /// In-plane Mandel stiffness for [xx, yy, sqrt(2)xy] strain and stress [Pa].
    pub fn mandel(&self)->&Complex3 {&self.mandel}
    /// Symmetric work-conjugate map from [epsilon_xx, epsilon_yy, gamma_xy]
    /// to [sigma_xx, sigma_yy, sigma_xy], with gamma_xy = 2 epsilon_xy [Pa].
    pub fn engineering(&self)->Complex3 {
        core::array::from_fn(|i|core::array::from_fn(|j|
            self.mandel[i][j].scale(ENGINEERING_TO_MANDEL[i]*ENGINEERING_TO_MANDEL[j])))
    }
    /// Recover ALL six physical strain amplitudes for the imposed engineering
    /// in-plane strain. The eliminated stresses vanish to solve roundoff.
    /// The complex Frobenius amplitude bounds the peak physical strain and
    /// must remain inside the original material's declared strain limit.
    pub fn strain_tensor(&self,in_plane:[C64;3])->Result<[C64;6]> {
        if in_plane.iter().any(|v|!cfinite(*v)) {return Err(invalid("nonfinite in-plane strain"));}
        let a=core::array::from_fn(|i|in_plane[i].scale(ENGINEERING_TO_MANDEL[i]));
        let b=complex_apply(&self.normal_extension,&a);
        let mut strain=[C64::ZERO;6];
        for i in 0..3 {
            strain[IN_PLANE[i]]=a[i].scale(1.0/WEIGHTS[IN_PLANE[i]]);
            strain[OUT_OF_PLANE[i]]=(-b[i]).scale(1.0/WEIGHTS[OUT_OF_PLANE[i]]);
        }
        let norm2=(0..6).fold(0.0,|sum,i|sum+WEIGHTS[i]*WEIGHTS[i]
            *(strain[i].re*strain[i].re+strain[i].im*strain[i].im));
        if strain.iter().any(|v|!cfinite(*v)) || !norm2.is_finite()
            || det::sqrt(norm2)>self.strain_limit {
            return Err(invalid("plane-stress strain exceeds the material domain"));
        }
        Ok(strain)
    }
    /// Homogeneous Kirchhoff-Love section about its middle plane. Membrane
    /// and bending matrices are h*C and h^3*C/12; the coupling block is zero.
    /// No transverse-shear, through-thickness layering or shell geometry is inferred.
    pub fn section(&self,thickness_m:f64)->Result<HarmonicPlateSection> {
        let moment=thickness_m*thickness_m*thickness_m/12.0;
        if !thickness_m.is_finite() || thickness_m<=0.0 || !moment.is_finite() || moment<=0.0 {
            return Err(invalid("unrepresentable positive section thickness"));
        }
        let c=self.engineering();
        let membrane=core::array::from_fn(|i|core::array::from_fn(|j|c[i][j].scale(thickness_m)));
        let bending=core::array::from_fn(|i|core::array::from_fn(|j|c[i][j].scale(moment)));
        let section=HarmonicPlateSection {membrane,bending};
        if section.membrane.iter().flatten().chain(section.bending.iter().flatten()).any(|v|!cfinite(*v)) {
            return Err(invalid("section stiffness overflow"));
        }
        Ok(section)
    }
}

/// Frequency-specific homogeneous section, retaining every coupling coefficient.
#[derive(Clone,Debug,PartialEq)]
pub struct HarmonicPlateSection {
    /// Engineering membrane matrix [N/m].
    pub membrane:Complex3,
    /// Engineering bending matrix [N m].
    pub bending:Complex3,
}

#[cfg(test)]
mod tests;
