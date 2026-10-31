//! Data-only axial/shear roles for the original mean-strain material owner.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::bank::relaxation::StretchingRelaxationSelection;
use fs_material::section::relaxation::LayerRelaxationProperties;

#[derive(Clone, Debug)]
pub struct Input {
    pub max_states: usize,
    pub core: LayerRelaxationProperties,
    pub winding: Option<LayerRelaxationProperties>,
}

fn requirements(properties: &LayerRelaxationProperties, out: &mut Vec<Requirement>) -> Result<()> {
    let equilibrium = std::iter::once((
        &properties.equilibrium_modulus,
        Pressure::DIMS,
        Admit::StrictlyPositive,
    ));
    let branches = properties.branches.iter().flat_map(|(delta, tau)| {
        [
            (delta, Pressure::DIMS, Admit::NonNegative),
            (tau, Dims([0, 0, 1, 0, 0, 0]), Admit::StrictlyPositive),
        ]
    });
    for (name, dims, admission) in equilibrium.chain(branches) {
        if name.is_empty() || name.len() > 256 {
            return Err("stretching-memory property names need 1..=256 bytes".into());
        }
        if let Some(previous) = out.iter().find(|p| p.name() == name) {
            if previous.dims() != dims {
                return Err(
                    "one stretching-memory source property has conflicting SI quantities".into(),
                );
            }
        } else {
            out.push(Requirement::try_new(name.as_str(), dims, admission)?);
        }
    }
    Ok(())
}

impl Input {
    pub fn validate(&self) -> Result<()> {
        if !(1..=256).contains(&self.max_states) {
            return Err("stretching-memory needs 1..=256 original history states".into());
        }
        for selection in std::iter::once(&self.core).chain(self.winding.as_ref()) {
            if selection.branches.len() > 8 {
                return Err("each stretching-memory source admits at most eight arms".into());
            }
            let mut pairs = std::collections::BTreeSet::new();
            if selection.branches.iter().any(|pair| !pairs.insert(pair)) {
                return Err("stretching-memory repeats a source-property pair".into());
            }
            requirements(selection, &mut Vec::new())?;
        }
        Ok(())
    }
    pub fn add_core_requirements(&self, out: &mut Vec<Requirement>) -> Result<()> {
        requirements(&self.core, out)
    }
    pub fn add_winding_requirements(&self, out: &mut Vec<Requirement>) -> Result<()> {
        requirements(
            self.winding
                .as_ref()
                .ok_or("wound stretching-memory needs explicit winding-memory roles")?,
            out,
        )
    }
    pub fn preflight(&self, spec: &input::Spec) -> Result<()> {
        self.validate()?;
        if spec.geometric_stretching.is_none()
            || spec.longitudinal_waves.is_some()
            || spec.prepared_linear.is_some()
        {
            return Err("stretching-memory requires the original exact geometric-stretching owner without longitudinal waves".into());
        }
        let mut count = 0usize;
        let mut has_winding = false;
        for string in &spec.strings {
            let layers = match &string.section {
                input::WireSection::Homogeneous => 1,
                input::WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
                input::WireSection::CloseWound { .. } => {
                    has_winding = true;
                    count += self
                        .winding
                        .as_ref()
                        .ok_or("every wound mean-strain source needs explicit winding-memory")?
                        .branches
                        .len();
                    1
                }
                input::WireSection::Measured { .. } => {
                    return Err(
                        "whole-string measurements do not determine axial relaxation spectra"
                            .into(),
                    );
                }
            };
            count = layers
                .checked_mul(self.core.branches.len())
                .and_then(|n| n.checked_add(count))
                .ok_or("complete stretching-memory extent overflow")?;
        }
        if !has_winding && self.winding.is_some() {
            return Err("winding-memory has no original wound member".into());
        }
        if count > self.max_states {
            return Err(
                "complete core and winding histories exceed stretching-memory MAX_STATES".into(),
            );
        }
        Ok(())
    }
    pub fn bind_after_preload(
        &self,
        performance: MultiBowPerformance,
        spec: &input::Spec,
        gate: &CancelGate,
    ) -> Result<MultiBowPerformance> {
        self.preflight(spec)?;
        let selections = spec
            .strings
            .iter()
            .enumerate()
            .map(|(string, source)| {
                let layers = match &source.section {
                    input::WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
                    _ => 1,
                };
                StretchingRelaxationSelection {
                    string,
                    layers: vec![self.core.clone(); layers],
                    winding: matches!(&source.section, input::WireSection::CloseWound { .. })
                        .then(|| self.winding.clone())
                        .flatten(),
                }
            })
            .collect();
        Ok(performance.with_material_stretching_relaxation(
            selections,
            None,
            self.max_states,
            gate,
        )?)
    }
}
