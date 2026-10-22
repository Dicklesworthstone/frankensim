//! Performance declarations select original material data, never a new solver.
use super::*;
pub use fs_couple::bowed_string::runtime::reciprocal::bank::relaxation::
    StringBendingRelaxationSelection as Selection;

#[derive(Clone, Debug)]
pub struct Input {
    max_states: usize,
    selections: Vec<Selection>,
}
impl Input {
    pub fn new(max_states: usize, selections: Vec<Selection>, strings: &[input::Wire]) -> Result<Self> {
        let value = Self { max_states, selections };
        value.validate(strings)?;
        Ok(value)
    }
    fn validate(&self, strings: &[input::Wire]) -> Result<()> {
        if !(1..=256).contains(&self.max_states) || self.selections.is_empty()
            || self.selections.len()>strings.len() {
            return Err("bending-memory needs 1..=256 states and nonempty original string selections".into());
        }
        for (index, selection) in self.selections.iter().enumerate() {
            if selection.string>=strings.len()
                || self.selections[..index].iter().any(|s|s.string==selection.string)
                || !(1..=8).contains(&selection.properties.len()) {
                return Err("string-memory needs distinct existing strings and 1..=8 source pairs".into());
            }
            if matches!(&strings[selection.string].bending_loss,input::BendingLoss::Prony { .. }) {
                return Err("one string cannot select both original-specimen Prony loss and additional bending-memory".into());
            }
            let mut roles=BTreeMap::new();
            let mut pairs=std::collections::BTreeSet::new();
            for (modulus,time) in &selection.properties {
                if !pairs.insert((modulus,time)) {
                    return Err("string-memory repeats a physical source pair".into());
                }
                for (name,dims) in [(modulus,Pressure::DIMS),(time,Dims([0,0,1,0,0,0]))] {
                    if name.is_empty() || name.len()>256
                        || roles.insert(name,dims).is_some_and(|previous|previous!=dims) {
                        return Err("string-memory property names must be bounded with consistent quantities".into());
                    }
                }
            }
        }
        Ok(())
    }
    pub fn preflight(&self, spec: &input::Spec) -> Result<()> {
        self.validate(&spec.strings)?;
        if spec.prepared_linear.is_some() {
            return Err("bending-memory requires the original exact owner; prepared-linear has no material strain states".into());
        }
        Ok(())
    }
    /// Resolve one shared pack once with the union of explicitly selected
    /// properties. Body packs and unrelated string packs gain no requirement.
    pub fn add_requirements(&self, path: &str, strings: &[input::Wire],
        requirements: &mut Vec<Requirement>) -> Result<()> {
        for selection in &self.selections {
            let wire=strings.get(selection.string).ok_or("absent bending-memory string")?;
            let uses_pack=wire.pack==path || match &wire.section {
                input::WireSection::Homogeneous=>false,
                input::WireSection::Bonded {sleeves,..}=>sleeves.iter().any(|s|s.pack==path),
            };
            if !uses_pack { continue; }
            for (modulus,time) in &selection.properties {
                for (name,dims,admit) in [(modulus,Pressure::DIMS,Admit::NonNegative),
                    (time,Dims([0,0,1,0,0,0]),Admit::StrictlyPositive)] {
                    if let Some(existing)=requirements.iter().find(|r|r.name()==name) {
                        if existing.dims()!=dims {
                            return Err("one string material property cannot supply incompatible quantities".into());
                        }
                    } else { requirements.push(Requirement::try_new(name.as_str(),dims,admit)?); }
                }
            }
        }
        Ok(())
    }
    /// This bound joins original-specimen and body histories in one setup budget.
    pub const fn max_memory_states(&self) -> usize { self.max_states }
    pub fn bind(&self, bank: ReciprocalStringBank, gate: &CancelGate) -> Result<ReciprocalStringBank> {
        poll(gate)?;
        Ok(bank.with_material_string_bending_relaxation(self.selections.clone(),None,self.max_states,gate)?)
    }
}
