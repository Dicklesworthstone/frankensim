//! A loss selection resolves source data and invokes the shared physics owner.
//! No material history, force law or time stepper is implemented by the decoder.
use super::*;
use fs_couple::plate_viscosity::relaxation::{InitialPlateMemory, PlateRelaxationLimits};
use fs_couple::plate_viscosity::relaxation::tensor::TensorPlateRelaxationLimits;
use fs_material::visco_tensor::source::{TensorBranchProperties,TensorElasticity,TensorPronyProperties};

#[derive(Clone, Debug)]
pub(super) enum PlateLoss {
    Tensor,
    RegionalMaxwell {
        properties: Vec<(String, String)>,
        initial: InitialPlateMemory,
        limits: PlateRelaxationLimits,
    },
    TensorMaxwell {
        strain_limit: String,
        branches: Vec<TensorBranchProperties>,
        initial: InitialPlateMemory,
        limits: TensorPlateRelaxationLimits,
        sharing: Option<usize>,
    },
}
impl PlateLoss {
    pub(super) fn parse(record: Option<&str>, rows:&mut Records<'_>) -> Result<Self> {
        let Some(record) = record else { return Ok(Self::Tensor); };
        let mut fields = record.split_ascii_whitespace();
        if fields.next() != Some("loss") { return Err("expected optional body loss record".into()); }
        match fields.next() {
            Some("tensor") if fields.next().is_none() => Ok(Self::Tensor),
            Some(kind @ ("tensor-maxwell" | "tensor-maxwell-shared")) => {
                let initial=match fields.next() {
                    Some("relaxed")=>InitialPlateMemory::Relaxed,
                    Some("unrelaxed")=>InitialPlateMemory::Unrelaxed,
                    _=>return Err("tensor body memory must declare relaxed or unrelaxed initial strain".into()),
                };
                let strain_limit=key(fields.next().ok_or("missing tensor strain property")?)?;
                let count:usize=fields.next().ok_or("missing tensor branch count")?.parse()?;
                let max_material_points:usize=fields.next().ok_or("missing tensor material point limit")?.parse()?;
                let max_memory_states:usize=fields.next().ok_or("missing tensor history limit")?.parse()?;
                let max_projection_terms:usize=fields.next().ok_or("missing tensor projection limit")?.parse()?;
                let sharing=if kind=="tensor-maxwell-shared" {
                    let work:usize=fields.next().ok_or("missing tensor coordinate-sharing work limit")?.parse()?;
                    if work==0 || work>67_108_864 {
                        return Err("tensor coordinate sharing requires a bounded positive factor-work limit".into());
                    }
                    Some(work)
                } else {None};
                if !(1..=8).contains(&count) || !(1..=4096).contains(&max_material_points)
                    || max_memory_states>256 || max_projection_terms==0
                    || max_projection_terms>67_108_864 || fields.next().is_some() {
                    return Err("tensor body memory requires complete 1..=8 branches and bounded point/history/work limits".into());
                }
                let mut branches=Vec::with_capacity(count);
                for _ in 0..count {
                    let row=rows.row("tensor_branch",22)?;
                    let relaxation_time=key(row[0])?;
                    let mut stiffness_upper=core::array::from_fn(|_|String::new());
                    for (out,value) in stiffness_upper.iter_mut().zip(&row[1..]) { *out=key(value)?; }
                    branches.push(TensorBranchProperties {stiffness_upper,relaxation_time});
                }
                Ok(Self::TensorMaxwell {strain_limit,branches,initial,sharing,limits:TensorPlateRelaxationLimits {
                    max_modes:input::MAX_PLATE_BODY_MODES,max_material_points,max_memory_states,max_projection_terms,
                }})
            },
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
            _ => Err("body loss must be tensor, regional-maxwell, tensor-maxwell or tensor-maxwell-shared".into()),
        }
    }
    pub(super) fn is_tensor(&self) -> bool { matches!(self, Self::Tensor) }
    pub(super) fn preflight(&self, spec: &input::Spec) -> Result<()> {
        if !matches!(self, Self::Tensor) && spec.prepared_linear.is_some() {
            return Err("body material memory requires the original exact owner; the selected prepared mechanical-only owner has no internal strain states".into());
        }
        Ok(())
    }
    pub(super) fn requirements(&self, model: PlateMaterialModel) -> Result<Vec<Requirement>> {
        match self {
            Self::Tensor => super::requirements(model),
            Self::TensorMaxwell {..} => {
                let mut selected=vec![Requirement::try_new("density",Density::DIMS,Admit::StrictlyPositive)?];
                for property in self.tensor_properties(model).expect("tensor selection").requirements()? {
                    if let Some(existing)=selected.iter().find(|r|r.name()==property.name()) {
                        if existing.dims()!=property.dims() {
                            return Err("one tensor body property cannot select incompatible quantities".into());
                        }
                    } else {selected.push(property);}
                }
                Ok(selected)
            },
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
    pub(super) fn max_memory_states(&self) -> usize {
        match self {
            Self::Tensor => 0,
            Self::RegionalMaxwell { limits, .. } => limits.max_memory_states,
            Self::TensorMaxwell {limits,..}=>limits.max_memory_states,
        }
    }
    pub(super) fn bind(&self, bank: ReciprocalStringBank, gate: &CancelGate)
        -> Result<ReciprocalStringBank>
    {
        poll(gate)?;
        Ok(match self {
            Self::Tensor => bank,
            Self::RegionalMaxwell { properties, initial, limits } =>
                bank.with_material_plate_relaxation(properties, *initial, *limits, gate)?,
            Self::TensorMaxwell {initial,limits,sharing,..} => {
                let chart=bank.source().material_chart().ok_or("tensor body source requires its original regional plate")?;
                let selections=chart.specimen.regions().iter().map(|region|
                    self.tensor_properties(region.input.model).expect("tensor selection")).collect::<Vec<_>>();
                match sharing {
                    Some(work)=>bank.with_shared_tensor_material_plate_relaxation(&selections,*initial,*limits,*work,gate)?,
                    None=>bank.with_tensor_material_plate_relaxation(&selections,*initial,*limits,gate)?,
                }
            },
        })
    }

    fn tensor_properties(&self,model:PlateMaterialModel)->Option<TensorPronyProperties> {
        let Self::TensorMaxwell {strain_limit,branches,..}=self else {return None;};
        Some(TensorPronyProperties {elasticity:if model==PlateMaterialModel::Isotropic {
            TensorElasticity::Isotropic
        } else {TensorElasticity::Orthotropic},strain_limit:strain_limit.clone(),branches:branches.clone()})
    }
}

fn key(value:&str)->Result<String> {
    if value.is_empty() || value.len()>128 || value.chars().any(char::is_control) {
        return Err("tensor property names must be nonempty tokens within 128 bytes".into());
    }
    Ok(value.into())
}
