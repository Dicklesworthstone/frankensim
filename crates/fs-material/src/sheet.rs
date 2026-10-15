//! Effective sheets: source-owned areal mass and INDEPENDENT membrane/bending laws.
//!
//! FrankenTextiles' `ftx-textile-catalog/1` distinguishes a characterized sheet
//! from its constituent fibers. This adapter preserves that distinction in the
//! existing material state-point system. N/m is not Pa; measured N m bending is
//! not inferred from membrane stiffness, Poisson ratio or thickness cubed.
//! See docs/FRANKENTEXTILES_INTEGRATION.md for source and interchange semantics.

use crate::state_point::{MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec};

/// Reference mass per sheet area [kg/m2], NOT constituent bulk density.
pub const AREAL_DENSITY: &str = "areal_density";
/// Full sheet thickness at the characterized condition [m].
pub const THICKNESS: &str = "thickness";
/// Principal-axis effective membrane data, in N/m except dimensionless nu12.
pub const MEMBRANE_PROPERTIES: [&str; 4] = [
    "warp_stiffness", "weft_stiffness", "poisson_warp_weft", "shear_stiffness",
];
/// Independent D11, D22, D12, D66 in N m; NOT coupon flexural rigidity in N m2.
pub const BENDING_PROPERTIES: [&str; 4] = [
    "bending_warp", "bending_weft", "bending_coupling", "bending_twist",
];
/// Membrane resultant and energy-per-reference-area dimensions.
pub const MEMBRANE_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);
/// Moment-resultant / engineering-curvature stiffness dimensions.
pub const BENDING_DIMS: Dims = Dims([2, 1, -2, 0, 0, 0]);

/// A consumer chooses complete independent capabilities, never inferred ones.
#[derive(Debug, Clone, Copy)]
pub struct SheetCapabilities {
    /// Require all four equilibrium membrane constants.
    pub membrane: bool,
    /// Require all four independent directional bending constants.
    pub bending: bool,
}

/// Original source refusal or a physically unsupported sheet conversion.
#[derive(Debug, Clone, PartialEq)]
pub enum SheetError {
    /// The shared immutable material query failed.
    Material(MaterialStatePointError),
    /// Missing capability, incompatible meaning, nonfinite or unstable data.
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for SheetError {
    fn from(value: MaterialStatePointError) -> Self { Self::Material(value) }
}
impl core::fmt::Display for SheetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self { Self::Material(e) => e.fmt(f), Self::Invalid(s) => write!(f, "sheet: {s}") }
    }
}
impl std::error::Error for SheetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Material(e) => Some(e), _ => None }
    }
}

/// One whole source bundle: thermal/contact/other selected properties survive
/// either mechanical lowering. Refresh at changed material state; this is not
/// a live thermal/mechanical history or a rule for homogenizing fiber mixtures.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedSheet {
    material: ResolvedMaterialStatePoint,
    areal_density: f64,
    thickness: f64,
}
fn property(m: &ResolvedMaterialStatePoint, key: &'static str, dims: Dims)
    -> Result<f64, SheetError> {
    let p = m.property(key).ok_or(SheetError::Invalid(key))?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !p.is_constant_at_fixed_temperature() || !p.value_si().is_finite() {
        return Err(SheetError::Invalid("sheet needs static or temperature-only coefficients in the declared sheet units"));
    }
    Ok(p.value_si())
}
impl ResolvedSheet {
    /// Adopt the original source object, including every query/receipt/property.
    /// Only mass/thickness are required here. Each constitutive consumer below
    /// requires its own full data; absent bending does not become a t-cubed law.
    pub fn try_from_material(material: ResolvedMaterialStatePoint) -> Result<Self, SheetError> {
        let areal_density = property(&material, AREAL_DENSITY, Dims([-2, 1, 0, 0, 0, 0]))?;
        let thickness = property(&material, THICKNESS, Dims([1, 0, 0, 0, 0, 0]))?;
        if areal_density <= 0.0 || thickness <= 0.0 {
            return Err(SheetError::Invalid("positive reference areal density and full thickness required"));
        }
        Ok(Self { material, areal_density, thickness })
    }
    /// Complete original immutable material bundle.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint { &self.material }
    /// Reference mass per area [kg/m2]. Changing deformed area cannot create mass.
    pub const fn areal_density_kg_m2(&self) -> f64 { self.areal_density }
    /// Full characterized thickness [m], independent from the bending law.
    pub const fn thickness_m(&self) -> f64 { self.thickness }
    /// Equilibrium quadratic Green-strain law in material axes. Source values
    /// are effective membrane moduli, not instantaneous tangents of a nonlinear
    /// stress curve. This does not infer compression softening or Maxwell loss.
    pub fn membrane(&self) -> Result<QuadraticMembrane, SheetError> {
        let mut p = [0.0; 4];
        for (i, key) in MEMBRANE_PROPERTIES.iter().enumerate() {
            p[i] = property(&self.material, key, if i == 2 { Dims::NONE } else { MEMBRANE_DIMS })?;
        }
        QuadraticMembrane::new(p[0], p[1], p[2], p[3])
    }
    /// Characterized engineering-curvature stiffness, without any solid-plate
    /// conversion. Scalar seam/hinge bending and coupon N m2 cannot fill a D matrix.
    pub fn bending(&self) -> Result<SheetMatrix, SheetError> {
        let mut p = [0.0; 4];
        for (out, key) in p.iter_mut().zip(BENDING_PROPERTIES) {
            *out = property(&self.material, key, BENDING_DIMS)?;
        }
        SheetMatrix::orthotropic(p[0], p[1], p[2], p[3])
    }
}

/// Resolve only the selected independent sheet capabilities in one common query.
/// Extra capabilities can instead be requested with `resolve_material_state_point`
/// and adopted by `ResolvedSheet::try_from_material`, without re-query or data loss.
pub fn resolve_sheet(card: &MaterialCard, point: &QueryPoint,
    capabilities: SheetCapabilities, selection: MaterialPropertySelection)
    -> Result<ResolvedSheet, SheetError> {
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(AREAL_DENSITY, Dims([-2, 1, 0, 0, 0, 0]), ScalarAdmissibility::StrictlyPositive)?,
        ScalarPropertyRequirement::try_new(THICKNESS, Dims([1, 0, 0, 0, 0, 0]), ScalarAdmissibility::StrictlyPositive)?,
    ];
    if capabilities.membrane {
        for (i, key) in MEMBRANE_PROPERTIES.iter().enumerate() {
            requirements.push(ScalarPropertyRequirement::try_new(*key,
                if i == 2 { Dims::NONE } else { MEMBRANE_DIMS }, ScalarAdmissibility::Finite)?);
        }
    }
    if capabilities.bending {
        for key in BENDING_PROPERTIES {
            requirements.push(ScalarPropertyRequirement::try_new(key, BENDING_DIMS, ScalarAdmissibility::Finite)?);
        }
    }
    let sheet = ResolvedSheet::try_from_material(resolve_material_state_point(card, point, &requirements, selection)?)?;
    if capabilities.membrane { sheet.membrane()?; }
    if capabilities.bending { sheet.bending()?; }
    Ok(sheet)
}

/// Positive-definite symmetric engineering-strain/curvature matrix. The caller
/// determines N/m (membrane) or N m (bending); no thickness conversion occurs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SheetMatrix([f64; 9]);
impl SheetMatrix {
    /// Independent principal-axis direct, cross and engineering-shear stiffness.
    pub fn orthotropic(a: f64, b: f64, cross: f64, shear: f64) -> Result<Self, SheetError> {
        if [a,b,cross,shear].iter().any(|v| !v.is_finite()) || a <= 0.0 || b <= 0.0 || shear <= 0.0 {
            return Err(SheetError::Invalid("positive direct and shear stiffness and finite cross coupling required"));
        }
        let r = cross / a.sqrt() / b.sqrt();
        if !r.is_finite() || r.abs() >= 1.0 {
            return Err(SheetError::Invalid("sheet normal coupling must be strictly positive definite"));
        }
        Ok(Self([a,cross,0.0,cross,b,0.0,0.0,0.0,shear]))
    }
    /// Dense row-major coefficients, engineering shear/twist in the third slot.
    pub const fn coefficients(&self) -> &[f64; 9] { &self.0 }
    /// Material axis 1 is counterclockwise from mesh x. Transform by T^t D T,
    /// preserving energy and engineering factors of two; geometry is untouched.
    pub fn rotated(self, angle_rad: f64) -> Result<Self, SheetError> {
        if !angle_rad.is_finite() { return Err(SheetError::Invalid("finite material angle required")); }
        if angle_rad == 0.0 { return Ok(self); }
        let (s,c) = angle_rad.sin_cos();
        let (cc,ss,cs) = (c*c,s*s,c*s);
        let t = [[cc,ss,cs],[ss,cc,-cs],[-2.0*cs,2.0*cs,cc-ss]];
        let mut out = [0.0; 9];
        for i in 0..3 { for j in i..3 {
            let mut v = 0.0;
            for a in 0..3 { for b in 0..3 { v += t[a][i]*self.0[3*a+b]*t[b][j]; } }
            out[3*i+j] = v; out[3*j+i] = v;
        } }
        if out.iter().any(|v| !v.is_finite()) || [out[0],out[4],out[8]].iter().any(|v| *v <= 0.0) {
            return Err(SheetError::Invalid("unrepresentable rotated sheet stiffness"));
        }
        Ok(Self(out))
    }
}

/// Energy and its exact derivatives in `[E11,E22,2E12]` Green-strain coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SheetResponse {
    /// Energy per REFERENCE area [J/m2].
    pub energy_j_m2: f64,
    /// Work-conjugate second-Piola resultants `[S11,S22,S12]` [N/m].
    pub resultants_n_m: [f64; 3],
    /// Exact material tangent, dense row-major [N/m].
    pub tangent_n_m: [f64; 9],
}
/// Constitutive interface for sheet kinematics; assembly does not identify fibers
/// or choose material laws. Implementations must return consistent derivatives.
pub trait SheetLaw {
    /// Evaluate energy, stress and tangent at the same material Green strain.
    fn response(&self, strain: [f64; 3]) -> Result<SheetResponse, SheetError>;
}
/// Fixed-state effective orthotropic sheet, independent of bending and thickness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticMembrane(SheetMatrix);
impl QuadraticMembrane {
    /// Effective warp/weft moduli and shear [N/m], with reciprocal Poisson coupling.
    pub fn new(warp: f64, weft: f64, nu12: f64, shear: f64) -> Result<Self, SheetError> {
        if [warp,weft,nu12,shear].iter().any(|v| !v.is_finite()) || warp <= 0.0 || weft <= 0.0 || shear <= 0.0 {
            return Err(SheetError::Invalid("finite effective membrane constants with positive moduli required"));
        }
        let denominator = 1.0 - nu12*(nu12*weft/warp);
        if !denominator.is_finite() || denominator <= 0.0 { return Err(SheetError::Invalid("unstable reciprocal membrane compliance")); }
        Ok(Self(SheetMatrix::orthotropic(warp/denominator,weft/denominator,nu12*weft/denominator,shear)?))
    }
    /// Principal-axis matrix; use the physical element's grain rather than
    /// rotating both this law and its kinematics.
    pub const fn matrix(&self) -> &SheetMatrix { &self.0 }
}
impl SheetLaw for QuadraticMembrane {
    fn response(&self, strain: [f64; 3]) -> Result<SheetResponse, SheetError> {
        if strain.iter().any(|v| !v.is_finite()) { return Err(SheetError::Invalid("finite Green strain required")); }
        let mut stress = [0.0; 3];
        for i in 0..3 { for j in 0..3 { stress[i] += self.0.0[3*i+j]*strain[j]; } }
        let energy = 0.5*strain.iter().zip(stress).map(|(e,s)| e*s).sum::<f64>();
        if !energy.is_finite() || stress.iter().any(|v| !v.is_finite()) {
            return Err(SheetError::Invalid("unrepresentable sheet response"));
        }
        Ok(SheetResponse {energy_j_m2:energy,resultants_n_m:stress,tangent_n_m:self.0.0})
    }
}

/// Source-qualified sheet sensible heat on the shared material/thermal owners.
pub mod thermal;
