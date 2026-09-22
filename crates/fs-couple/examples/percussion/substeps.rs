//! Explicit hard-impact recovery inside the existing fixed output clock.
use super::{Mechanics, Error};
use fs_couple::render::plate::impact::ImpactSubstepConfig;

pub fn option(args: &mut Vec<String>) -> Result<Option<ImpactSubstepConfig>, Error> {
    let mut bounds = None;
    let mut i = 0;
    while i < args.len() {
        if args[i] != "--impact-substeps" { i += 1; continue; }
        if bounds.is_some() { return Err("--impact-substeps may be supplied only once".into()); }
        let values = args.get(i+1..i+3).ok_or("--impact-substeps needs maximum depth and solve attempts")?;
        let max_depth = values[0].parse::<u8>()?;
        let max_attempts = values[1].parse::<usize>()?;
        if max_depth > 10 || !(1..=2047).contains(&max_attempts) {
            return Err("--impact-substeps needs depth 0..10 and solve attempts 1..2047".into());
        }
        bounds = Some(ImpactSubstepConfig { max_depth, max_attempts });
        args.drain(i..i+3);
    }
    Ok(bounds)
}

impl Mechanics {
    /// Prepare this exact model (if needed), without resetting an accepted tick
    /// or force schedule. The caller explicitly chose this numerical realization.
    pub fn with_impact_substeps(self, bounds: ImpactSubstepConfig) -> Result<Self, Error> {
        match self {
            Self::Reference(s) => Ok(Self::Substepped(s.prepare()?.with_substeps(bounds)?)),
            Self::Nonlinear(s) => Ok(Self::Substepped(s.with_substeps(bounds)?)),
            Self::Substepped(_) => Err("impact refinement is already configured".into()),
            Self::Prepared(_) => Err("--impact-substeps requires nonlinear splash/drum/drum-stretch mechanics; modal/snare contact has a different time owner".into()),
            Self::Driven {inner, drive} => Ok(Self::Driven {
                inner:Box::new((*inner).with_impact_substeps(bounds)?), drive }),
        }
    }
}

#[cfg(test)]
mod tests;
