//! Version-3 contact records for the existing modal performance reader.
use super::{Reader, ModalPerformanceError, ModalForceVoice, input, read_attachment};
use crate::render::RenderError;
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::{ModalContact,ModalContactConfig};
use fs_dcontact::Obstacle;

pub(super) fn read(reader:&mut Reader<'_>,voices:&[ModalForceVoice],weights:&mut usize)
    -> Result<(ModalContact,ModalContactConfig),ModalPerformanceError>
{
    let mut row=reader.row("contact_limits")?;
    let config=ModalContactConfig {
        max_iterations:row.count(128)?,maximum_force_n:row.scalar()?,maximum_penetration_m:row.scalar()?,
        force_absolute_tolerance_n:row.scalar()?,force_relative_tolerance:row.scalar()?,
    };
    row.finish()?;
    let mut row=reader.row("contact")?;
    let stiffness=row.scalar()?;let alpha=row.scalar()?;let chi=row.scalar()?;
    let gap=row.scalar()?;let weight=row.scalar()?;
    let provenance=row.word()?;
    if provenance.len()>1024 {return Err(input(row.line,"contact source label exceeds 1024 bytes"));}
    let provenance=provenance.to_string();row.finish()?;
    let left=read_attachment(reader,"contact_left",voices,weights)?;
    let right=read_attachment(reader,"contact_right",voices,weights)?;
    let law=Obstacle::new(vec![-1.0],1,1,vec![gap],vec![weight],stiffness,alpha,provenance)
        .and_then(|law|law.with_internal_loss(chi))
        .map_err(|error|RenderError::Coupled(ModalCouplingError::ContactLaw(error)))?;
    Ok((ModalContact {left,right,law},config))
}
