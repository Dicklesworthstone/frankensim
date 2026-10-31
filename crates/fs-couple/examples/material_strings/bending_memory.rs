//! Performance declarations select original material data, never a new solver.
use super::*;
pub use fs_couple::bowed_string::runtime::reciprocal::bank::relaxation::
    LayeredStringBendingRelaxationSelection as Selection;
pub use fs_material::section::relaxation::LayerBendingRelaxationProperties as LayerProperties;

#[derive(Clone, Debug)]
pub struct Input {
    max_states: usize,
    selections: Vec<Selection>,
}
impl Input {
    pub fn uniform_selection(string: usize, properties: Vec<(String, String)>,
        strings: &[input::Wire]) -> Result<Selection> {
        if properties.is_empty() {
            return Err("string-memory requires 1..=8 source pairs".into());
        }
        let wire = strings.get(string).ok_or("absent bending-memory string")?;
        let count = match &wire.section {
            input::WireSection::Measured { .. } => return Err("whole-string measurements do not supply layer relaxation laws".into()),
            input::WireSection::Homogeneous | input::WireSection::CloseWound { .. } => 1,
            input::WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
        };
        Ok(Selection { string, layers: vec![LayerProperties {
            equilibrium_modulus: "young_modulus".into(), branches: properties,
        }; count] })
    }
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
            let wire = strings.get(selection.string).ok_or("absent bending-memory string")?;
            let layers = match &wire.section {
                input::WireSection::Measured { .. } => return Err("whole-string measurements do not supply layer relaxation laws".into()),
                input::WireSection::Homogeneous | input::WireSection::CloseWound { .. } => 1,
                input::WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
            };
            if self.selections[..index].iter().any(|s|s.string==selection.string)
                || selection.layers.len() != layers {
                return Err("string memory requires distinct strings and exactly one selection per physical layer".into());
            }
            if !matches!(&wire.bending_loss, input::BendingLoss::KelvinVoigt) {
                return Err("additional bending-memory requires the original Kelvin-Voigt bending selection".into());
            }
            for layer in &selection.layers {
                if layer.branches.len() > 8 || layer.equilibrium_modulus.is_empty()
                    || layer.equilibrium_modulus.len() > 256 {
                    return Err("layer-memory requires a bounded equilibrium property and at most eight source pairs".into());
                }
                let mut roles=BTreeMap::new();
                roles.insert(&layer.equilibrium_modulus, Pressure::DIMS);
                let mut pairs=std::collections::BTreeSet::new();
                for (modulus,time) in &layer.branches {
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
            let sleeves = match &wire.section {
                input::WireSection::Measured { .. } => return Err("whole-string measurements do not supply layer relaxation laws".into()),
                // Bending memory belongs to the core. A sliding helical wrap
                // is not another bonded bending layer and needs no such claims.
                input::WireSection::Homogeneous | input::WireSection::CloseWound { .. } => &[][..],
                input::WireSection::Bonded { sleeves, .. } => sleeves.as_slice(),
            };
            let packs = std::iter::once(wire.pack.as_str())
                .chain(sleeves.iter().map(|s| s.pack.as_str()));
            for (pack, layer) in packs.zip(&selection.layers) {
                if pack != path { continue; }
                let equilibrium = std::iter::once((&layer.equilibrium_modulus,
                    Pressure::DIMS, Admit::StrictlyPositive));
                let branches = layer.branches.iter().flat_map(|(modulus,time)| [
                    (modulus,Pressure::DIMS,Admit::NonNegative),
                    (time,Dims([0,0,1,0,0,0]),Admit::StrictlyPositive),
                ]);
                for (name,dims,admit) in equilibrium.chain(branches) {
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
        Ok(bank.with_material_layered_string_bending_relaxation(self.selections.clone(),None,self.max_states,gate)?)
    }
}
