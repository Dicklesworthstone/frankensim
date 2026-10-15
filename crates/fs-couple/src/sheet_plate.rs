//! Source-bound effective-sheet bending on the EXISTING DKT and modal owners.
//!
//! Measured D is independent of membrane stiffness and thickness. This is flat
//! small-deflection plate dynamics, not a nonlinear garment-drape solver.
use fs_material::sheet::{ResolvedSheet, SheetError};
use fs_plate::{AssemblyOptions, PlateError, PlateMesh, PlateModel, PlateSection,
    PlateSectionField, assemble_with_sections};

/// Bounded topology before allocating per-triangle sections and DKT operators.
#[derive(Clone, Copy, Debug)]
pub struct SheetPlateBudget {
    /// Maximum source mesh nodes.
    pub max_nodes: usize,
    /// Maximum source triangles and material bindings.
    pub max_triangles: usize,
}
/// Original material/plate error or a refused source-to-geometry binding.
#[derive(Debug)]
pub enum SheetPlateError {
    /// Original shared material owner.
    Material(SheetError),
    /// Original DKT assembly/modal owner.
    Plate(PlateError),
    /// Incomplete binding, unrepresentable inertia, or exhausted topology budget.
    Invalid(&'static str),
    /// No candidate published after cancellation.
    Cancelled,
}
impl From<SheetError> for SheetPlateError { fn from(e:SheetError)->Self {Self::Material(e)} }
impl From<PlateError> for SheetPlateError { fn from(e:PlateError)->Self {Self::Plate(e)} }
impl core::fmt::Display for SheetPlateError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"sheet plate: {self:?}")}
}
impl std::error::Error for SheetPlateError {}

/// All original source bundles, grain directions and topology beside the actual
/// DKT pencil. Downstream modal/acoustic consumers use `model()`, not another
/// eigensolver. Fields are immutable so the source and operator cannot drift.
#[derive(Clone, Debug)]
pub struct SheetPlate {
    mesh: PlateMesh,
    sheets: Vec<ResolvedSheet>,
    grain_rad: Vec<f64>,
    model: PlateModel,
    reference_mass_kg: f64,
}
impl SheetPlate {
    /// One explicit whole-sheet source and material-axis angle per triangle.
    /// Mass is areal mass times reference area. `density = areal_mass/thickness`
    /// is ONLY a representation for the existing plate inertia; it is not a
    /// solid-fiber density or a constitutive law. Slope rotary inertia retains
    /// the DKT owner's uniform-through-thickness approximation, mu*t^2/12.
    /// Bending uses the independently characterized D, never E*t^3/12.
    ///
    /// Uniform prestress is an explicitly supplied load in `options`; it is not
    /// inferred from a material label or a fabricated membrane strain. Spatial
    /// source/angle ordering is exactly the original triangle ordering.
    /// Cancellation is polled per binding and around the bounded native assembly;
    /// this synchronous adapter cannot preempt an individual DKT assembly call.
    pub fn assemble(mesh:PlateMesh,sheets:Vec<ResolvedSheet>,grain_rad:Vec<f64>,
        boundary:&[usize],options:AssemblyOptions,budget:SheetPlateBudget,
        mut cancelled:impl FnMut()->bool)->Result<Self,SheetPlateError> {
        if cancelled() {return Err(SheetPlateError::Cancelled);}
        if mesh.nodes.len()>budget.max_nodes || mesh.tris.len()>budget.max_triangles
            || sheets.len()!=mesh.tris.len() || grain_rad.len()!=mesh.tris.len() {
            return Err(SheetPlateError::Invalid("bounded complete per-triangle sources and grain required"));
        }
        let mut sections=Vec::with_capacity(sheets.len());
        for (sheet,&angle) in sheets.iter().zip(&grain_rad) {
            if cancelled() {return Err(SheetPlateError::Cancelled);}
            let density=sheet.areal_density_kg_m2()/sheet.thickness_m();
            if !density.is_finite() || density<=0.0 {
                return Err(SheetPlateError::Invalid("unrepresentable sheet inertia lowering"));
            }
            sections.push(PlateSection {d:*sheet.bending()?.rotated(angle)?.coefficients(),
                thickness:sheet.thickness_m(),density});
        }
        let model=assemble_with_sections(&mesh,PlateSectionField::PerElement(&sections),
            boundary,&[],&options)?;
        let mut mass=0.0;
        // Native assembly admitted every node index and CCW triangle above.
        for (tri,sheet) in mesh.tris.iter().zip(&sheets) {
            if cancelled() {return Err(SheetPlateError::Cancelled);}
            let [a,b,c]=tri.map(|i|mesh.nodes[i]);
            let area=0.5*((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0));
            mass+=sheet.areal_density_kg_m2()*area;
        }
        if !mass.is_finite() || mass<=0.0 {return Err(SheetPlateError::Invalid("unrepresentable reference sheet mass"));}
        if cancelled() {return Err(SheetPlateError::Cancelled);}
        Ok(Self {mesh,sheets,grain_rad,model,reference_mass_kg:mass})
    }
    /// Original source mesh; no meshing or material averaging occurred.
    pub const fn mesh(&self)->&PlateMesh {&self.mesh}
    /// Original complete source bundles in triangle order.
    pub fn sheets(&self)->&[ResolvedSheet] {&self.sheets}
    /// Material axis 1 measured counterclockwise from reference x.
    pub fn grain_rad(&self)->&[f64] {&self.grain_rad}
    /// Original DKT model, directly accepted by fs-plate/fs-modal consumers.
    pub const fn model(&self)->&PlateModel {&self.model}
    /// Conserved reference mass [kg], not apparent density times deformed volume.
    pub const fn reference_mass_kg(&self)->f64 {self.reference_mass_kg}
}
