//! One explicit receiver policy for pressure, harmonic loading and played sound.
//! Moving a receiver must never change the BEM boundary solve or impedance.
use super::{Boundary,Medium,Specification,RATE};
use fs_bem::{helmholtz::{self,RadiationSolution},near_field,panel3d::SpherePanels};
use fs_math::c64::C64;

pub struct ReceiverSet<'a> {
    surface:&'a SpherePanels,
    points:Vec<[f64;3]>,
    medium:Medium,
    geometry:Option<near_field::Geometry<'a>>,
    delays:Vec<f64>,
}
impl<'a> ReceiverSet<'a> {
    pub fn for_spec(boundary:&'a Boundary,spec:&Specification)->Result<Self,String> {
        Self::new(boundary,&spec.receivers,spec.medium,spec.near_field_receivers)
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
        Ok(Self {surface:&boundary.surface,points:points.to_vec(),medium,geometry,delays})
    }
    pub fn delays_s(&self)->&[f64] {&self.delays}
    pub fn prepare(&self,k:f64)->Result<Evaluation<'a>,String> {
        if let Some(geometry)=&self.geometry {
            Ok(Evaluation::Near(geometry.prepare(k,self.medium,near_field::Options::default())
                .map_err(|e|format!("near-field receiver quadrature: {e}"))?))
        } else {
            Ok(Evaluation::Centroid {surface:self.surface,points:self.points.clone(),medium:self.medium})
        }
    }
}
/// Rows are prepared once per frequency, shared by all solved modal fields.
pub enum Evaluation<'a> {
    Centroid {surface:&'a SpherePanels,points:Vec<[f64;3]>,medium:Medium},
    Near(near_field::Prepared<'a>),
}
impl Evaluation<'_> {
    pub fn pressure(&self,solution:&RadiationSolution)->Result<Vec<C64>,String> {
        match self {
            Self::Centroid {surface,points,medium}=>
                helmholtz::exterior_pressure_at_points(surface,solution,*medium,points).map_err(|e|e.to_string()),
            Self::Near(prepared)=>prepared.evaluate(solution).map(|result|result.pressure).map_err(|e|e.to_string()),
        }
    }
}
