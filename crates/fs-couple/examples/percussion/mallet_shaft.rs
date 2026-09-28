//! One loaded shaft/head inertia, four material histories, reciprocal moments.
use super::{Error, Spec, FeltStriker};
use super::super::shaft_playing::Built;
use fs_couple::render::plate::impact::striker::flexible::FlexibleStriker;

impl Spec {
    pub fn admit_shaft(&self, selected:bool) -> Result<(),Error> {
        match (self.head.is_some(),selected) {
            (true,false)=>Err("v2 physical mallet head requires its own --flexible-stick selection".into()),
            (false,true)=>Err("v1 mallet mass is already effective; use an explicit v2 physical head to couple shaft flexure without double-counting inertia".into()),
            _=>Ok(()),
        }
    }
    pub fn attach_inertia(&self, shaft:&FlexibleStriker) -> Result<FlexibleStriker,Error> {
        self.admit_shaft(true)?;
        let head=self.head.ok_or("missing physical head inertia")?;
        Ok(shaft.with_tip_inertia(self.jaw.mass_kg,head.rotary_kg_m2)?)
    }
    /// Reuse the original footprint, skin clearances, area split and all felt
    /// histories. The temporary scalar FeltStriker is NEVER time stepped; its
    /// inertia is replaced, not appended, by the already loaded rigid mode.
    /// Each site moves as u_tip + xi*theta_tip in the declared bending plane.
    pub fn bind_shaft(&self, mut tip:FeltStriker, hand:usize, shaft:&mut Built)
        -> Result<FeltStriker,Error>
    {
        if hand>=2 {return Err("invalid mallet hand".into());}
        self.admit_shaft(shaft.ports[hand].is_some())?;
        let Some(head)=self.head else {return Ok(tip);};
        if !shaft.loaded[hand] {return Err("felt head requires the loaded, not bare, shaft pencil".into());}
        let ports=shaft.ports[hand].as_ref().ok_or("missing loaded shaft ports")?;
        let rigid=ports.rigid_coordinate();
        let translation=ports.tip_row(shaft.total)?;
        let rotation=ports.tip_slope_row(shaft.total)?;
        if tip.port.coordinate!=rigid || tip.pads.len()!=4 {
            return Err("mallet footprint does not match its loaded shaft layout".into());
        }
        let (mut body,weight)=shaft.bodies[hand].take().ok_or("loaded shaft body was already consumed")?;
        // Launch the entire shaft/head in the exact rigid mode. The face is
        // initially horizontal with the supplied gap, so its affine rotation
        // datum is the initial shaft angle, not an extra strain or reset.
        let position=-self.jaw.initial_gap_m/weight;
        if !position.is_finite() || (self.jaw.initial_gap_m>0.0 && position==0.0) {
            return Err("loaded mallet gap is not representable".into());
        }
        body.initial[0].displacement_m_sqrt_kg=position;
        let initial_rotation=rotation[rigid]*position;
        let (s,c)=head.azimuth_rad.sin_cos();
        let a=self.radius_m/std::f64::consts::SQRT_2;
        let levers=[a*c,a*s,-a*c,-a*s];
        for (pad,lever) in tip.pads.iter_mut().zip(levers) {
            if pad.weights.len()!=shaft.total || pad.weights[rigid]!=tip.port.inverse_sqrt_mass {
                return Err("mallet must replace its own scalar contact coordinate".into());
            }
            for i in 0..shaft.total {
                if i!=rigid && (translation[i]!=0.0 || rotation[i]!=0.0) && pad.weights[i]!=0.0 {
                    return Err("mallet target aliases a shaft coordinate".into());
                }
            }
            pad.weights[rigid]=0.0;
            for i in 0..shaft.total {pad.weights[i]+=translation[i]+lever*rotation[i];}
            pad.precompression_m-=lever*initial_rotation;
            if !pad.precompression_m.is_finite() || pad.weights.iter().any(|v|!v.is_finite()) {
                return Err("loaded mallet contact projection overflow".into());
            }
        }
        tip.body=body;tip.port.inverse_sqrt_mass=weight;
        Ok(tip)
    }
}

#[cfg(test)]
#[path="mallet_shaft_tests.rs"]
mod tests;
