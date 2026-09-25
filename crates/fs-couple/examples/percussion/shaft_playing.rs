//! Flexible shaft composition shared by single-cymbal and drum performances.
//! Geometry owns beam modes; the existing impact owner owns contact and time.
use super::{Error, Experiment, ImpactBody, Stroke, flexible_sticks, mallets, sticks};
use super::mechanics::drive::{Input, Program, SpatialInput};
use fs_couple::render::plate::impact::striker::flexible::{FlexibleStriker, StrikerPorts};
use std::io::Write;

#[derive(Default)]
pub struct Selection {
    pub first: Option<FlexibleStriker>,
    pub second: Option<FlexibleStriker>,
}
impl Selection {
    pub fn options(args: &mut Vec<String>) -> Result<Self, Error> {
        Ok(Self { first: flexible_sticks::option(args, "--flexible-stick")?,
            second: flexible_sticks::option(args, "--second-flexible-stick")? })
    }
    pub fn enabled(&self) -> bool { self.first.is_some() || self.second.is_some() }
    pub fn admit(&self, command: &str, second: Option<Stroke>, mallets: &mallets::Selection) -> Result<(), Error> {
        sticks::admit_command(self.enabled(), command)?;
        if self.second.is_some() && second.is_none() {
            return Err("--second-flexible-stick requires --second-stick-position-m X Y".into());
        }
        for (shaft,mallet) in [self.first.as_ref(),self.second.as_ref()].into_iter()
            .zip([mallets.first.as_ref(),mallets.second.as_ref()]) {
            if let Some(mallet)=mallet {mallet.admit_shaft(shaft.is_some())?;}
        }
        Ok(())
    }
    pub fn build_with_mallets(&self, base:usize, second_coordinate:usize, first:Stroke,
        second:Option<Stroke>, dt_s:f64, mallets:&mallets::Selection) -> Result<Built,Error> {
        let load=|shaft:Option<&FlexibleStriker>,mallet:Option<&mallets::Spec>| -> Result<Option<FlexibleStriker>,Error> {
            match (shaft,mallet) {
                (Some(s),Some(m))=>Ok(Some(m.attach_inertia(s)?)),
                (Some(s),None)=>Ok(Some(s.clone())),
                (None,Some(m))=>{m.admit_shaft(false)?;Ok(None)},
                (None,None)=>Ok(None),
            }
        };
        let loaded=Self {first:load(self.first.as_ref(),mallets.first.as_ref())?,
            second:load(self.second.as_ref(),mallets.second.as_ref())?};
        loaded.build(base,second_coordinate,first,second,dt_s)
    }
    /// Append flexure AFTER every original solid coordinate, including wires
    /// and mute jaws. No head, shell, second rigid stick or carrier is moved.
    /// Cavity and material-private states are constructed afterwards as usual.
    pub fn build(&self, base: usize, second_coordinate: usize, first: Stroke,
        second: Option<Stroke>, dt_s: f64) -> Result<Built, Error>
    {
        if self.second.is_some() && second.is_none() {
            return Err("missing physical second-stick station for flexible shaft".into());
        }
        if base==0 || second.is_some() && (second_coordinate==0 || second_coordinate>=base) {
            return Err("rigid sticks must lie in the original solid prefix before shaft flexure".into());
        }
        let extra = [self.first.as_ref(), self.second.as_ref()].into_iter().flatten()
            .try_fold(0_usize, |n,s| n.checked_add(s.elastic_modes()))
            .ok_or("flexible shaft layout overflow")?;
        let total = base.checked_add(extra).ok_or("flexible shaft layout overflow")?;
        let mut built = Built { total, rigid: [0,second_coordinate], bodies: [None,None],
            elastic: Vec::new(), ports: [None,None], loaded:[false;2] };
        let mut offset = base;
        for (hand,shaft) in [self.first.as_ref(),self.second.as_ref()].into_iter().enumerate() {
            if let Some(shaft) = shaft {
                let stroke = if hand==0 { first } else { second.ok_or("missing second launch")? };
                let launch = flexible_sticks::Launch::new(Some(shaft),built.rigid[hand],offset,
                    total,dt_s,stroke.speed_m_s)?;
                built.bodies[hand] = Some((launch.body,launch.weight));
                built.elastic.push(launch.elastic.ok_or("missing admitted shaft flexure")?);
                built.ports[hand] = launch.ports;
                built.loaded[hand] = shaft.has_tip_inertia();
                offset += shaft.elastic_modes();
            }
        }
        Ok(built)
    }
}

pub struct Built {
    pub total: usize,
    rigid: [usize;2],
    pub bodies: [Option<(ImpactBody,f64)>;2],
    pub elastic: Vec<ImpactBody>,
    pub ports: [Option<StrikerPorts>;2],
    pub loaded: [bool;2],
}
impl Built {
    pub fn tip_row(&self, hand: usize, rigid_weight: f64) -> Result<Vec<f64>, Error> {
        match &self.ports[hand] {
            Some(p) => Ok(p.tip_row(self.total)?),
            None => {
                let mut row=vec![0.0;self.total];
                let value=row.get_mut(self.rigid[hand]).ok_or("stick outside mechanical basis")?;
                *value=rigid_weight; Ok(row)
            }
        }
    }
}

impl super::Mechanics {
    /// Every input is admitted before attaching one joint player clock. This
    /// wraps the selected time owner, not each shaft or each force component.
    pub fn with_player_drives(self, inputs: Vec<Input>, spatial: Vec<SpatialInput>,
        dt_s: f64, steps: u64, modes: usize) -> Result<Self, Error>
    {
        if matches!(&self, Self::Driven { .. }) || modes > self.state().len()/2 {
            return Err("player drive needs one unnested, dimensionally admitted mechanical image".into());
        }
        let drive=super::mechanics::drive::StickDrive::new_mixed(inputs,spatial,dt_s,steps,modes)?;
        Ok(Self::Driven{inner:Box::new(self),drive})
    }
}

/// Physical hand forces are not tip forces when a shaft bends. Keep their
/// signs and lever ratios, and let the caller add carrier/mute inputs before
/// admitting the complete shared performance.
pub fn player_inputs(e: &Experiment, programs: [Option<Program>;2])
    -> Result<(Vec<Input>,Vec<SpatialInput>),Error>
{
    let mut scalar=Vec::new();let mut spatial=Vec::new();
    for (hand,program) in programs.into_iter().enumerate() {
        if let Some(program)=program {
            if let Some(p)=&e.flexible_sticks[hand] {
                spatial.push(SpatialInput{program,weights:p.hand_row(e.force.len())?});
            } else {
                let p=if hand==0 {sticks::Port{coordinate:0,weight:e.stick_weight}}
                    else {e.second_stick.ok_or("missing physical second-stick force port")?};
                scalar.push(Input{program,coordinate:p.coordinate,tip_weight:p.weight});
            }
        }
    }
    Ok((scalar,spatial))
}

pub fn tip_motion(e: &Experiment, hand: usize) -> Result<(f64,f64),Error> {
    let state=e.system.state();
    if let Some(p)=&e.flexible_sticks[hand] {
        let o=p.observe(state)?;Ok((o.tip_displacement_m,o.tip_velocity_m_s))
    } else {
        let p=if hand==0 {sticks::Port{coordinate:0,weight:e.stick_weight}}
            else {e.second_stick.ok_or("missing second-stick observation")?};
        Ok((state[2*p.coordinate]*p.weight,state[2*p.coordinate+1]*p.weight))
    }
}
pub fn header(e: &Experiment, out: &mut impl Write) -> Result<(),Error> {
    for (hand,p) in e.flexible_sticks.iter().enumerate() {
        if p.is_some() {let n=hand+1;
            write!(out,",stick{n}_tip_down_m,stick{n}_tip_speed_m_s,stick{n}_hand_down_m,stick{n}_hand_speed_m_s,stick{n}_bending_j")?;}
    }
    Ok(())
}
pub fn row(e: &Experiment, out: &mut impl Write) -> Result<(),Error> {
    for p in e.flexible_sticks.iter().flatten() {
        let o=p.observe(e.system.state())?;
        write!(out,",{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",o.tip_displacement_m,o.tip_velocity_m_s,
            o.hand_displacement_m,o.hand_velocity_m_s,o.flexural_energy_j)?;
    }
    Ok(())
}

#[cfg(test)]
#[path="shaft_playing_tests.rs"]
mod tests;
