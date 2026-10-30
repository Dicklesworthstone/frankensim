//! Explicit axial material roles use the original source loader and wave owner.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::relaxation::
    LongitudinalRelaxationSelection;
use fs_material::section::relaxation::LayerRelaxationProperties;

#[derive(Clone, Debug)]
pub struct Input {
    pub max_states: usize,
    pub properties: LayerRelaxationProperties,
}
impl Input {
    pub fn validate(&self) -> Result<()> {
        if !(1..=256).contains(&self.max_states) || self.properties.branches.is_empty()
            || self.properties.branches.len() > 8 {
            return Err("longitudinal-memory needs 1..=256 states and 1..=8 source modulus/time pairs".into());
        }
        let mut requirements = Vec::new();
        self.add_requirements(&mut requirements)?;
        let mut pairs = std::collections::BTreeSet::new();
        if self.properties.branches.iter().any(|pair| !pairs.insert(pair)) {
            return Err("longitudinal-memory repeats one source-property pair".into());
        }
        Ok(())
    }
    pub fn add_requirements(&self, requirements: &mut Vec<Requirement>) -> Result<()> {
        let equilibrium = std::iter::once((&self.properties.equilibrium_modulus,
            Pressure::DIMS, Admit::StrictlyPositive));
        let branches = self.properties.branches.iter().flat_map(|(delta, tau)| [
            (delta, Pressure::DIMS, Admit::NonNegative),
            (tau, Dims([0,0,1,0,0,0]), Admit::StrictlyPositive),
        ]);
        for (name, dims, admit) in equilibrium.chain(branches) {
            if name.is_empty() || name.len() > 256 {
                return Err("axial-memory source property names must have 1..=256 bytes".into());
            }
            if let Some(prior) = requirements.iter().find(|r| r.name() == name) {
                if prior.dims() != dims {
                    return Err("one axial material source property has incompatible SI quantities".into());
                }
            } else { requirements.push(Requirement::try_new(name.as_str(), dims, admit)?); }
        }
        Ok(())
    }
    pub fn preflight(&self, spec: &input::Spec) -> Result<()> {
        self.validate()?;
        let Some((resolution, _)) = &spec.longitudinal_waves else {
            return Err("longitudinal-memory requires retained longitudinal-waves".into());
        };
        if spec.prepared_linear.is_some() {
            return Err("longitudinal-memory requires the exact owner with material history".into());
        }
        let mut states = 0usize;
        for wire in &spec.strings {
            let layers = match &wire.section {
                input::WireSection::Homogeneous => 1,
                input::WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
                input::WireSection::CloseWound { .. } => return Err(
                    "quasistatic wound sections have no resolved longitudinal material history".into()),
                input::WireSection::Measured { .. } => return Err(
                    "measured whole strings do not supply axial material history".into()),
            };
            states = layers.checked_mul(self.properties.branches.len())
                .and_then(|n| n.checked_mul(resolution.axial_modes.len()))
                .and_then(|n| n.checked_add(states)).ok_or("axial history count overflow")?;
        }
        if states > self.max_states {
            return Err("complete per-layer axial histories exceed longitudinal-memory MAX_STATES".into());
        }
        Ok(())
    }
    pub fn bind(&self, bank: ReciprocalStringBank, spec: &input::Spec, gate: &CancelGate)
        -> Result<ReciprocalStringBank>
    {
        self.preflight(spec)?;
        let selections = spec.strings.iter().enumerate().map(|(string, wire)| {
            let layers = match &wire.section {
                input::WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
                _ => 1, // Wound sections already refused by preflight.
            };
            LongitudinalRelaxationSelection { string, layers: vec![self.properties.clone(); layers] }
        }).collect();
        Ok(bank.with_material_longitudinal_relaxation(selections, None, self.max_states, gate)?)
    }
}
