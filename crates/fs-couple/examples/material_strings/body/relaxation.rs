//! A loss selection resolves source data and invokes the shared physics owner.
//! No material history, force law or time stepper is implemented by the decoder.
use super::*;
use fs_couple::plate_viscosity::relaxation::{InitialPlateMemory, PlateRelaxationLimits};

#[derive(Clone, Debug)]
pub(super) enum PlateLoss {
    Tensor,
    RegionalMaxwell {
        properties: Vec<(String, String)>,
        initial: InitialPlateMemory,
        limits: PlateRelaxationLimits,
    },
}
impl PlateLoss {
    pub(super) fn parse(record: Option<&str>) -> Result<Self> {
        let Some(record) = record else { return Ok(Self::Tensor); };
        let mut fields = record.split_ascii_whitespace();
        if fields.next() != Some("loss") { return Err("expected optional body loss record".into()); }
        match fields.next() {
            Some("tensor") if fields.next().is_none() => Ok(Self::Tensor),
            Some("regional-maxwell") => {
                let initial = match fields.next() {
                    Some("relaxed") => InitialPlateMemory::Relaxed,
                    Some("unrelaxed") => InitialPlateMemory::Unrelaxed,
                    _ => return Err("body Maxwell memory must declare relaxed or unrelaxed initial strain".into()),
                };
                let max_memory_states: usize = fields.next().ok_or("missing body memory-state limit")?.parse()?;
                let max_projection_terms: usize = fields.next().ok_or("missing body projection-work limit")?.parse()?;
                let terms: usize = fields.next().ok_or("missing body Maxwell term count")?.parse()?;
                if max_memory_states > 256 || max_projection_terms == 0
                    || max_projection_terms as u64 > input::MAX_WORK || !(1..=8).contains(&terms) {
                    return Err("body Maxwell selection requires <=256 states, bounded positive work and 1..=8 terms".into());
                }
                let mut properties = Vec::with_capacity(terms);
                for _ in 0..terms {
                    let modulus = fields.next().ok_or("missing body excess-modulus property")?;
                    let time = fields.next().ok_or("missing body relaxation-time property")?;
                    if modulus.len() > 256 || time.len() > 256 || modulus == time
                        || properties.iter().any(|(a,b): &(String,String)| a == modulus && b == time) {
                        return Err("body Maxwell terms need distinct dimensioned property pairs within key limits".into());
                    }
                    properties.push((modulus.into(), time.into()));
                }
                if fields.next().is_some() { return Err("unconsumed body Maxwell loss fields".into()); }
                Ok(Self::RegionalMaxwell { properties, initial, limits: PlateRelaxationLimits {
                    max_modes: input::MAX_PLATE_BODY_MODES, max_projection_terms, max_memory_states,
                } })
            },
            _ => Err("body loss must be tensor or regional-maxwell".into()),
        }
    }
    pub(super) fn is_tensor(&self) -> bool { matches!(self, Self::Tensor) }
    pub(super) fn preflight(&self, spec: &input::Spec) -> Result<()> {
        if matches!(self, Self::RegionalMaxwell { .. }) && spec.prepared_linear.is_some() {
            return Err("body material memory requires the original exact owner; the selected prepared mechanical-only owner has no internal strain states".into());
        }
        Ok(())
    }
    pub(super) fn requirements(&self, model: PlateMaterialModel) -> Result<Vec<Requirement>> {
        match self {
            Self::Tensor => super::requirements(model),
            Self::RegionalMaxwell { properties, .. } => {
                let mut selected = elastic_requirements(model)?;
                for (modulus,time) in properties {
                    for (name,dims,admission) in [(modulus, Pressure::DIMS, Admit::NonNegative),
                        (time, Dims([0,0,1,0,0,0]), Admit::StrictlyPositive)] {
                        // A source time may serve several branches. Resolve the
                        // one original property once, not as duplicate demands.
                        if let Some(existing) = selected.iter().find(|r| r.name() == name) {
                            if existing.dims() != dims {
                                return Err("one body property key cannot have incompatible selected quantities".into());
                            }
                        } else {
                            selected.push(Requirement::try_new(name.as_str(),dims,admission)?);
                        }
                    }
                }
                Ok(selected)
            },
        }
    }
    pub(super) fn reserve_setup(&self, settings: &mut StringBankOptions) -> Result<()> {
        if let Self::RegionalMaxwell { limits, .. } = self {
            let n = settings.network.max_modes;
            let memory = limits.max_memory_states;
            let core = n.checked_mul(2).ok_or("body memory setup dimension overflow")?;
            let terms = core.checked_add(memory).and_then(|d| d.checked_mul(d))
                .and_then(|d| core.checked_mul(memory).and_then(|m| d.checked_add(m)))
                .and_then(|d| settings.network.max_setup_terms.checked_add(d))
                .ok_or("body memory setup work overflow")?;
            if terms as u64 > input::MAX_WORK { return Err("body memory exceeds the original instrument work budget".into()); }
            settings.network.max_setup_terms = terms;
        }
        Ok(())
    }
    pub(super) fn bind(&self, bank: ReciprocalStringBank, gate: &CancelGate)
        -> Result<ReciprocalStringBank>
    {
        poll(gate)?;
        Ok(match self {
            Self::Tensor => bank,
            Self::RegionalMaxwell { properties, initial, limits } =>
                bank.with_material_plate_relaxation(properties, *initial, *limits, gate)?,
        })
    }
}
