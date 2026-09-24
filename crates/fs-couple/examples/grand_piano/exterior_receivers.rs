//! One explicit receiver policy for pressure, harmonic loading and played sound.
//! Moving or rotating a receiver never changes the source solve or impedance.
use super::{Boundary,Medium,Specification,RATE};
use fs_bem::{helmholtz::{self,RadiationSolution},near_field,panel3d::SpherePanels};
use fs_math::c64::C64;

pub struct ReceiverSet<'a> {
    surface:&'a SpherePanels,
    points:Vec<[f64;3]>,
    medium:Medium,
    geometry:Option<near_field::Geometry<'a>>,
    delays:Vec<f64>,
    patterns:Vec<near_field::FirstOrder>,
}
impl<'a> ReceiverSet<'a> {
    pub fn for_spec(boundary:&'a Boundary,spec:&Specification)->Result<Self,String> {
        if spec.receiver_patterns.keys().any(|i|*i>=spec.receivers.len())
            || (!spec.near_field_receivers && spec.receiver_patterns.values().any(|p|p.pressure_fraction()<1.)) {
            return Err("directional receiver needs its actual near-field position and matching pattern index".into());
        }
        let mut out=Self::new(boundary,&spec.receivers,spec.medium,spec.near_field_receivers)?;
        for (&index,&pattern) in &spec.receiver_patterns {out.patterns[index]=pattern;}
        Ok(out)
    }
    pub fn new(boundary:&'a Boundary,points:&[[f64;3]],medium:Medium,near:bool)->Result<Self,String> {
        if !(1..=2).contains(&points.len()) || !medium.sound_speed.is_finite() || medium.sound_speed<=0.
            || !medium.density.is_finite() || medium.density<=0. {
            return Err("invalid receiver count or physical medium".into());
        }
        let geometry=if near {
            Some(near_field::Geometry::new(&boundary.surface,points,2.*medium.sound_speed/f64::from(RATE))
                .map_err(|e|format!("near-field receiver admission: {e}"))?)
        } else {None};
        let delays=if let Some(geometry)=&geometry {
            // Shortest last leg from ANY panel, including rigid scattering:
            // conservative flight lower bound, not a nearest-centroid delay.
            geometry.clearances_m().iter().map(|d|d/medium.sound_speed).collect::<Vec<_>>()
        } else {
            points.iter().map(|p| {
                // Preserve the original arithmetic for unselected receivers.
                let radius=(0..3).fold(0.0_f64,|n,c|n.hypot(p[c]-boundary.center[c]));
                (radius-boundary.radius)/medium.sound_speed
            }).collect::<Vec<_>>()
        };
        if points.iter().flatten().any(|v|!v.is_finite()) || delays.iter().any(|d|
            !d.is_finite() || !(2./f64::from(RATE)..=0.5).contains(d)) {
            return Err("receiver flight lower bound must be 2 output samples to 0.5 s; near-field needs surface clearance, centroid mode needs the enclosing sphere".into());
        }
        Ok(Self {surface:&boundary.surface,points:points.to_vec(),medium,geometry,delays,
            patterns:vec![near_field::FirstOrder::default();points.len()]})
    }
    pub fn delays_s(&self)->&[f64] {&self.delays}
    pub fn prepare(&self,k:f64)->Result<Evaluation<'a>,String> {
        if let Some(geometry)=&self.geometry {
            if self.patterns.iter().any(|p|p.pressure_fraction()<1.) {
                let prepared=geometry.prepare_velocity(k,self.medium,near_field::Options::default())
                    .map_err(|e|format!("directional receiver quadrature: {e}"))?;
                Ok(Evaluation::Directional {prepared,patterns:self.patterns.clone(),medium:self.medium})
            } else {
                // Explicit alpha=1 and omitted patterns use the identical scalar
                // path, not a gradient calculation multiplied by zero afterward.
                Ok(Evaluation::Near(geometry.prepare(k,self.medium,near_field::Options::default())
                    .map_err(|e|format!("near-field receiver quadrature: {e}"))?))
            }
        } else {
            Ok(Evaluation::Centroid {surface:self.surface,points:self.points.clone(),medium:self.medium})
        }
    }
}
/// Rows are prepared once per frequency, shared by all solved modal fields.
pub enum Evaluation<'a> {
    Centroid {surface:&'a SpherePanels,points:Vec<[f64;3]>,medium:Medium},
    Near(near_field::Prepared<'a>),
    Directional {prepared:near_field::Prepared<'a>,patterns:Vec<near_field::FirstOrder>,medium:Medium},
}
impl Evaluation<'_> {
    /// Pressure-equivalent output for the selected first-order pattern. Raw
    /// pressure for an omni. Source boundary pressure is never changed here.
    pub fn pressure(&self,solution:&RadiationSolution)->Result<Vec<C64>,String> {
        match self {
            Self::Centroid {surface,points,medium}=>
                helmholtz::exterior_pressure_at_points(surface,solution,*medium,points).map_err(|e|e.to_string()),
            Self::Near(prepared)=>prepared.evaluate(solution).map(|result|result.pressure).map_err(|e|e.to_string()),
            Self::Directional {prepared,patterns,medium}=>{
                let field=prepared.evaluate_velocity(solution).map_err(|e|e.to_string())?;
                if patterns.len()!=field.scalar.pressure.len() || patterns.len()!=field.particle_velocity_m_s.len() {
                    return Err("directional receiver field cardinality differs".into());
                }
                patterns.iter().zip(field.scalar.pressure).zip(field.particle_velocity_m_s)
                    .map(|((pattern,p),v)|pattern.observe(p,v,*medium).map_err(|e|e.to_string())).collect()
            }
        }
    }
}
