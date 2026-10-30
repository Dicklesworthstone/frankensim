//! Original layer spectra on the full strain of the moving-support image.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::relaxation::{
    LongitudinalRelaxationSelection,elastic_section,
};
use crate::render::plate::impact::relaxation::InitialMemory;
use crate::render::plate::impact::supported::longitudinal::relaxation::SupportedRelaxationArm;
use fs_material::section::relaxation::ConcentricPronySection;

/// Original source selections and complete admitted layer spectra. Their
/// evolving normalized strain memories remain in ImpactSystem, never here.
#[derive(Clone, Debug)]
pub struct SupportedMaterialRelaxation {
    selections: Vec<LongitudinalRelaxationSelection>,
    sections: Vec<ConcentricPronySection>,
    initial: InitialMemory,
    maximum_states: usize,
}
impl SupportedMaterialRelaxation {
    pub fn selections(&self)->&[LongitudinalRelaxationSelection] {&self.selections}
    pub fn sections(&self)->&[ConcentricPronySection] {&self.sections}
    pub const fn initial_memory(&self)->InitialMemory {self.initial}
    pub const fn maximum_states(&self)->usize {self.maximum_states}
}
impl SupportedMaterialStrings {
    /// Bind the SAME source uniaxial spectra used by fixed-end strings to the
    /// original supported strain epsilon=u_x+(y_x^2+z_x^2)/2. Each layer/arm
    /// retains independent mean and cosine histories. Both transverse planes,
    /// axial waves, affine endpoint motion and their reciprocal support force
    /// enter the original impact owner's nonlinear H/R equation together.
    ///
    /// Data-only changes to Delta E, tau or layer area change the entire
    /// mechanics; equilibrium EA, solid inertia and maintained preload survive.
    /// This is the declared bonded, moderate-slope reduction, not core/wrap
    /// slip, evolving temperature or a geometrically exact rod. Initial memory
    /// is explicit, and all copied bodies retain this material law automatically.
    ///
    /// The source bands cover DC through the complete support/torsion bound,
    /// including the hereditary tangent's strain-history envelope. The clock
    /// must resolve every tau. The original time owner refuses a history outside
    /// its admitted envelope or total hereditary stress that loses local tension.
    pub fn with_material_longitudinal_relaxation(mut self,
        selections:Vec<LongitudinalRelaxationSelection>, initial:InitialMemory,
        maximum_states:usize, gate:&CancelGate) -> Result<Self,ImpactError>
    {
        poll(gate)?;
        let Some(attachments)=&self.longitudinal else {
            return Err(invalid("supported material memory requires original longitudinal waves and geometry"));
        };
        if self.axial_material_relaxation.is_some() || selections.is_empty()
            || selections.len()>self.sources.len() || !(1..=256).contains(&maximum_states) {
            return Err(invalid("supported material memory needs distinct cold selections and a 1..=256 state budget"));
        }
        let mut sections=Vec::with_capacity(selections.len());let mut arms=Vec::new();
        let dt=1.0/f64::from(self.sample_rate_hz);
        for (index,selection) in selections.iter().enumerate() {
            poll(gate)?;
            let source=self.sources.get(selection.string)
                .ok_or_else(||invalid("supported memory names an absent material string"))?;
            if selections[..index].iter().any(|s|s.string==selection.string) {
                return Err(invalid("supported material string cannot select duplicate memory owners"));
            }
            let section=elastic_section(source.source())?.with_prony_bending(&selection.layers)
                .map_err(|e|ImpactError::Owner(e.to_string()))?;
            let (solid,_)=attachments[selection.string].lower(source)?;
            let same=|a:f64,b:f64|(a-b).abs()<=16.0*f64::EPSILON*a.abs().max(b.abs());
            if !same(section.elastic().line_density_kg_m(),solid.solid_linear_density_kg_m)
                || !same(section.elastic().axial_stiffness_n(),source.axial_rigidity_n()) {
                return Err(invalid("supported memory must preserve original solid inertia and equilibrium rigidity"));
            }
            for (spectrum,geometry) in section.spectra().iter().zip(section.elastic().layer_integrals()) {
                for &(delta,tau) in &spectrum.terms {
                    if dt/tau>2.0 {return Err(invalid("source clock does not resolve supported material memory"));}
                    if delta==0.0 {continue;}
                    arms.push(SupportedRelaxationArm {string:selection.string,
                        axial_rigidity_n:delta*geometry.area_m2(),relaxation_time_s:tau});
                }
            }
            sections.push(section);
        }
        let supported=self.supported.clone().with_longitudinal_relaxation(arms,initial,maximum_states)?;
        let mut upper=supported.clock_bounds().iter().copied().fold(0.0_f64,f64::max);
        for body in self.bodies.iter().skip(1) {
            let BodyPotential::Linear(omegas)=&body.potential else {
                return Err(invalid("supported memory lost the original torsion body"));
            };
            upper=omegas.iter().copied().fold(upper,f64::max);
        }
        fs_phs::linear_exact::admit_midpoint_band(upper,dt).map_err(ImpactError::TimeResolution)?;
        if sections.iter().any(|section| {let(lo,hi)=section.omega_band_rad_s();lo>0.0 || hi<upper}) {
            return Err(invalid("supported material spectrum must cover DC through the complete hereditary tangent bound"));
        }
        // A new material arm stiffens the common object. Existing equilibrium
        // elasticity and independently selected viscosity must support it too.
        let mut ceiling=sections.iter().map(|s|s.omega_band_rad_s().1).fold(f64::INFINITY,f64::min);
        for (attachment,source) in attachments.iter().zip(&self.sources) {
            // Reuse the complete original source walker: includes bending
            // material/loss, each bonded layer and wrapped torsional spectra.
            // Axial lowering alone intentionally strips those other channels.
            let mut band=(0.0,f64::INFINITY);
            crate::bowed_string::runtime::reciprocal::bank::harmonic::string_band(source.source(),&mut band)
                .map_err(|e|ImpactError::Owner(e.to_string()))?;
            if band.0>0.0 || band.1<upper {
                return Err(invalid("supported hereditary band exceeds an original string material or loss channel"));
            }
            ceiling=ceiling.min(band.1);
            if let Some((lo,hi))=attachment.lower(source)?.1 {
                if lo>0.0 || hi<upper {
                    return Err(invalid("supported hereditary band exceeds original equilibrium material support"));
                }
                ceiling=ceiling.min(hi);
            }
        }
        if let Some(laws)=&self.axial_material_loss {for law in laws {
            let(lo,hi)=law.omega_band_rad_s();
            if lo>0.0 || hi<upper {return Err(invalid("supported hereditary band exceeds original axial viscosity support"));}
            ceiling=ceiling.min(hi);
        }}
        self.supported=supported.bind_longitudinal_relaxation_domain(ceiling,dt)?;
        self.bodies[0].potential=BodyPotential::SupportedStrings(self.supported.clone());
        self.axial_material_relaxation=Some(SupportedMaterialRelaxation {selections,sections,initial,maximum_states});
        poll(gate)?;Ok(self)
    }
    pub fn axial_material_relaxation(&self)->Option<&SupportedMaterialRelaxation> {
        self.axial_material_relaxation.as_ref()
    }
}
