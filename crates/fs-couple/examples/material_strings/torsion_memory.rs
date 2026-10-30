//! Per-layer source selection for the original torsional material-memory owner.
use super::*;
pub use fs_material::section::relaxation::LayerRelaxationProperties as LayerProperties;

impl input::Wire {
    /// Reserve every declared arm before source I/O; zero-valued arms can use
    /// fewer states after source admission. Torsion has one angular plane.
    pub(super) fn torsion_memory_states(&self) -> Result<usize> {
        let Some(layers) = &self.torsion_memory else { return Ok(0); };
        let physical_layers = match &self.section {
            input::WireSection::Measured { .. } => return Err("whole-string measurements do not supply torsional layer laws".into()),
            input::WireSection::Homogeneous => 1,
            input::WireSection::CloseWound { .. } => return Err("close-wound torsion requires a resolved winding law".into()),
            input::WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
        };
        if self.twist_modes == 0 || layers.len() != physical_layers {
            return Err("torsion-memory requires torsional modes and one selection per physical layer".into());
        }
        let mut arms = 0usize;
        for layer in layers {
            if layer.branches.len() > 8 || !bounded_name(&layer.equilibrium_modulus) {
                return Err("torsion-layer needs a bounded equilibrium property and at most eight source pairs".into());
            }
            let mut roles = BTreeMap::new();
            roles.insert(&layer.equilibrium_modulus, Pressure::DIMS);
            let mut pairs = std::collections::BTreeSet::new();
            for (modulus, time) in &layer.branches {
                if !pairs.insert((modulus, time)) {
                    return Err("torsion-layer repeats a physical source pair".into());
                }
                for (name, dims) in [(modulus, Pressure::DIMS), (time, Dims([0,0,1,0,0,0]))] {
                    if !bounded_name(name)
                        || roles.insert(name, dims).is_some_and(|prior| prior != dims) {
                        return Err("torsional memory property names require consistent physical quantities".into());
                    }
                }
            }
            arms = arms.checked_add(layer.branches.len())
                .ok_or("torsional material arm count overflow")?;
        }
        let states = arms.checked_mul(self.twist_modes)
            .ok_or("torsional material memory-state count overflow")?;
        if states > 256 {
            return Err("torsional material memory exceeds the original exact-state ceiling".into());
        }
        Ok(states)
    }

    /// Only the corresponding original layer pack supplies its selected law.
    /// A pack shared by several layers/strings resolves their union once.
    pub(super) fn add_torsion_memory_requirements(&self, path: &str,
        requirements: &mut Vec<Requirement>) -> Result<()> {
        let Some(layers) = &self.torsion_memory else { return Ok(()); };
        self.torsion_memory_states()?;
        let sleeves = match &self.section {
            input::WireSection::Measured { .. } => return Err("whole-string measurements do not supply torsional layer laws".into()),
            input::WireSection::Homogeneous => &[][..],
            input::WireSection::CloseWound { .. } => return Err("close-wound torsion requires a resolved winding law".into()),
            input::WireSection::Bonded { sleeves, .. } => sleeves.as_slice(),
        };
        let packs = std::iter::once(self.pack.as_str())
            .chain(sleeves.iter().map(|layer| layer.pack.as_str()));
        for (pack, layer) in packs.zip(layers) {
            if pack != path { continue; }
            let equilibrium = std::iter::once((&layer.equilibrium_modulus,
                Pressure::DIMS, Admit::StrictlyPositive));
            let branches = layer.branches.iter().flat_map(|(modulus,time)| [
                (modulus,Pressure::DIMS,Admit::NonNegative),
                (time,Dims([0,0,1,0,0,0]),Admit::StrictlyPositive),
            ]);
            for (name,dims,admit) in equilibrium.chain(branches) {
                if let Some(prior) = requirements.iter().find(|r| r.name() == name) {
                    if prior.dims() != dims {
                        return Err("one source property cannot supply incompatible torsional quantities".into());
                    }
                } else { requirements.push(Requirement::try_new(name.as_str(),dims,admit)?); }
            }
        }
        Ok(())
    }
}

fn bounded_name(name: &str) -> bool {
    !name.is_empty() && name.len() <= 256
        && !name.chars().any(char::is_whitespace)
}
