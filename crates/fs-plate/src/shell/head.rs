//! Circular prestressed films through the existing plate assembly, not a new
//! membrane stepper. Geometry, film elasticity/density, thickness and tension
//! are separate inputs. A drumhead brand does not specify the installed tension.
//! The fixed rim uses w=0 with free slope; finite bending stiffness is retained.
//! This baseline omits tuning-lug nonuniformity, coating, contact and enclosed
//! air. Those are physical coupled participants, not baked frequency corrections.
use super::profile::{ProfileBudget,ProfileStation,revolve};
use crate::{AssemblyOptions,EdgeSupport,PlateError,PlateMesh,PlateModel,PlateSection,assemble};

/// Explicit circular film and numerical sampling.
#[derive(Clone,Copy,Debug)]
pub struct TensionedDiskSpec {
    /// Vibrating radius [m], not an unexamined nominal shell diameter.
    pub radius_m:f64,
    /// Actual film thickness [m].
    pub thickness_m:f64,
    /// In-plane Young modulus [Pa].
    pub young_pa:f64,
    /// Poisson ratio.
    pub poisson:f64,
    /// Film density [kg/m^3].
    pub density_kg_m3:f64,
    /// Isotropic installed membrane tension [N/m], explicitly tuned.
    pub tension_n_m:f64,
    /// Number of radial intervals, including the center fan.
    pub radial_intervals:usize,
    /// Number of circumferential vertices per ring.
    pub azimuths:usize,
}
/// Actual film pencil and geometry, ready for the existing modal facility.
pub struct TensionedDisk {
    /// Supplied physical values.
    pub spec:TensionedDiskSpec,
    /// Generated planar triangles, with a unique center.
    pub mesh:PlateMesh,
    /// Actual section used in both stiffness and mass.
    pub section:PlateSection,
    /// Existing DKT + P1 membrane-prestress reduced pencil.
    pub model:PlateModel,
    /// Film mass in the sampled polygon [kg], before eliminating the rim DOFs.
    pub mass_kg:f64,
}
impl TensionedDisk {
    /// Construct a fixed-rim film without any authored resonant frequencies.
    /// # Errors
    /// Invalid physical values, discretization/budget or existing FEM refusal.
    pub fn new(spec:TensionedDiskSpec,budget:ProfileBudget)->Result<Self,PlateError> {
        let bad=||PlateError::BadSection{what:"disk requires positive finite radius/tension and bounded radial sampling"};
        if !spec.radius_m.is_finite() || spec.radius_m<=0.0 || !spec.tension_n_m.is_finite()
            || spec.tension_n_m<=0.0 || spec.radial_intervals==0
            || spec.radial_intervals.checked_add(1).is_none_or(|n|n>budget.max_nodes) {return Err(bad());}
        let stations:Vec<_>=(0..=spec.radial_intervals).map(|i|ProfileStation {
            radius_m:spec.radius_m*i as f64/spec.radial_intervals as f64,
            height_m:0.0,thickness_m:spec.thickness_m,
        }).collect();
        let generated=revolve(&stations,spec.azimuths,spec.young_pa,spec.poisson,spec.density_kg_m3,&[],&[],budget)?;
        let mesh=PlateMesh::from_unstructured(generated.mesh.nodes.iter().map(|p|(p[0],p[1])).collect(),generated.mesh.tris)?;
        let section=PlateSection::isotropic(spec.young_pa,spec.poisson,spec.thickness_m,spec.density_kg_m3)?;
        let model=assemble(&mesh,&section,&generated.outer_nodes,&[],&AssemblyOptions {
            pretension:spec.tension_n_m,support:EdgeSupport::SimplySupported,
        })?;
        Ok(Self{spec,mesh,section,model,mass_kg:generated.mass_kg})
    }
    /// Project one unit pressure into a supplied mass-normalized mode. Using
    /// the SAME weights for volume displacement gives reciprocal cavity work.
    /// # Errors
    /// A mode not in this pencil's finite reduced coordinate space.
    pub fn modal_area(&self,mode:&[f64])->Result<f64,PlateError> {
        if mode.len()!=self.model.free || mode.iter().any(|v|!v.is_finite()) {
            return Err(PlateError::BadSection{what:"disk mode does not match its finite pencil"});
        }
        let mut integral=0.0;
        for tri in &self.mesh.tris {
            let (x0,y0)=self.mesh.nodes[tri[0]];let (x1,y1)=self.mesh.nodes[tri[1]];let (x2,y2)=self.mesh.nodes[tri[2]];
            let area=0.5*((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0));
            for &node in tri {integral+=area/3.0*self.model.dof_map[3*node].map_or(0.0,|i|mode[i]);}
        }
        if !integral.is_finite() {return Err(PlateError::BadSection{what:"disk modal area overflow"});}
        Ok(integral)
    }
}
