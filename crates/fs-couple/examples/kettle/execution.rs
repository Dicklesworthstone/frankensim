//! Bounded numerical choices around the SAME head/felt/air equations.
//! No force-envelope repair, dropped contact or observer-visible partial tick.
use super::Error;
use fs_couple::render::plate::impact::{ImpactSubstepConfig, ImpactSystem, SubsteppedImpactSystem};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Settings {
    pub substeps: ImpactSubstepConfig,
    pub backtracking_halvings: u8,
}
impl Default for Settings {
    fn default() -> Self {
        Self { substeps: ImpactSubstepConfig {max_depth:8,max_attempts:511}, backtracking_halvings:8 }
    }
}
impl Settings {
    // Cold CLI admission duplicates only the shared owners' integer ceilings.
    // Physical time/energy representability remains enforced by those owners.
    pub fn validate(self) -> Result<(), Error> {
        if self.substeps.max_depth>10 || !(1..=2047).contains(&self.substeps.max_attempts)
            || self.backtracking_halvings>16 {
            return Err("execution needs substep depth 0..10, 1..2047 solve attempts and 0..16 Newton backtracking halvings".into());
        }
        Ok(())
    }
    pub fn prepare(self, system: ImpactSystem) -> Result<SubsteppedImpactSystem, Error> {
        self.validate()?;
        let mut solver = system.prepare_analytic()?;
        solver.set_backtracking_limit(self.backtracking_halvings)?;
        Ok(solver.with_substeps(self.substeps)?)
    }
}

/// Consume both numerical options transactionally. A one-leaf, undamped
/// comparison is --impact-substeps 0 1 --newton-backtracking 0.
pub fn option(args: &mut Vec<String>) -> Result<Settings, Error> {
    let mut settings = Settings::default();
    let (mut subdivisions, mut backtracking) = (false, false);
    let mut rest = Vec::with_capacity(args.len()); let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--impact-substeps" => {
                if subdivisions || args.len()-i<3 {
                    return Err("supply --impact-substeps DEPTH ATTEMPTS at most once".into());
                }
                settings.substeps=ImpactSubstepConfig {max_depth:args[i+1].parse()?, max_attempts:args[i+2].parse()?};
                subdivisions=true; i+=3;
            }
            "--newton-backtracking" => {
                if backtracking || args.len()-i<2 {
                    return Err("supply --newton-backtracking HALVINGS at most once".into());
                }
                settings.backtracking_halvings=args[i+1].parse()?; backtracking=true; i+=2;
            }
            _ => {rest.push(args[i].clone());i+=1;}
        }
    }
    settings.validate()?;
    *args=rest; Ok(settings)
}

#[cfg(test)]
#[path="execution_tests.rs"]
mod tests;
