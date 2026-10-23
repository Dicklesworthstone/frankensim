//! One exterior load on the actual shared body, retained through bank rebuilds.
use super::*;
use crate::render::plate::impact::radiation::Model;
use fs_blake3::{ContentHash, DomainHasher};

/// Geometry-derived planar Rayleigh producer on the retained body basis.
pub mod rayleigh;

/// Closed curved material surfaces through the shared exterior Helmholtz owner.
pub mod bem;

/// Original medium, complete body basis and passive impedance realization.
/// A supplied model keeps caller authority; the geometry producer owns the fit
/// and its continuum/finite-band accuracy. Nothing is inferred from a timbre.
#[derive(Clone, Debug)]
pub struct MaterialBodyRadiation {
    model: Model,
    gas: GasState,
    basis_identity: ContentHash,
    maximum_poles: usize,
}
impl MaterialBodyRadiation {
    pub const fn model(&self) -> &Model { &self.model }
    pub const fn gas(&self) -> &GasState { &self.gas }
    pub const fn basis_identity(&self) -> ContentHash { self.basis_identity }
}

impl MaterialStringBankSource {
    /// Complete exterior source retained through scheduled and simultaneous bows.
    pub const fn body_radiation(&self) -> Option<&MaterialBodyRadiation> { self.body_radiation.as_ref() }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn restore_body_radiation(
        &self, system: FrictionModalSystem, gate: &CancelGate,
    ) -> Result<FrictionModalSystem, ReciprocalStringError> {
        let Some(air) = &self.body_radiation else { return Ok(system); };
        let (identity, modes) = self.body_radiation_basis(&system, gate)?;
        if air.gas != self.gas || air.basis_identity != identity {
            return Err(invalid("exterior load no longer matches the original material body basis and gas"));
        }
        let system = system.with_radiation_load(&air.model, &modes, air.maximum_poles, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.validate_body_memory_frequencies(&system, gate)?;
        self.validate_material_bridge_frequencies(&system, gate)?;
        Ok(system)
    }

    fn body_radiation_basis(&self, system: &FrictionModalSystem, gate: &CancelGate)
        -> Result<(ContentHash, Vec<usize>), ReciprocalStringError>
    {
        poll(gate)?;
        let body = self.strings.len();
        let component = system.components().get(body)
            .ok_or_else(|| invalid("exterior radiation requires the original material body"))?;
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.material-body-radiation-basis.v1");
        if let Some(chart) = &self.material_chart {
            hash.update(&[0]);
            hash.update(chart.specimen.specimen_identity().as_bytes());
            let basis = self.material_chart_basis.as_ref()
                .ok_or_else(|| invalid("exterior load requires its retained complete body basis"))?;
            for dof in &basis.dof_map {
                hash.update(&(dof.map_or(u64::MAX, |i| i as u64)).to_le_bytes());
            }
            hash_rows(&mut hash, &basis.mass_normalized_modes, gate)?;
        } else if let Some(plate) = &self.material_plate {
            hash.update(&[1]);
            hash.update(plate.specimen.specimen_identity().as_bytes());
        } else if let Some(image) = &self.material_shell_image {
            hash.update(&[2]);
            hash_shell(&mut hash, image.source(), gate)?;
            hash_rows(&mut hash, image.physical_modes(), gate)?;
        } else if let Some(image) = &self.material_shell_assembly_image {
            hash.update(&[3]);
            for (component, rows) in image.source().components.iter().zip(image.physical_modes()) {
                hash_shell(&mut hash, &component.specimen, gate)?;
                hash_rows(&mut hash, rows, gate)?;
            }
        } else {
            return Err(invalid("raw modal bodies have no material/geometry basis for an exterior source binding"));
        }
        for mode in component.modes() {
            hash.update(&mode.angular_frequency_rad_s.to_bits().to_le_bytes());
        }
        let start = system.components()[..body].iter().map(|c| c.modes().len()).sum::<usize>();
        Ok((hash.finalize(), (start..start+component.modes().len()).collect()))
    }
}

impl ReciprocalStringBank {
    /// Couple a producer's complete exterior model to ALL actual body velocities
    /// in original mass-normalized mode order. Air reacts on the body and hence
    /// the bridges, every string, bow traction and accepted pressure acceleration.
    /// The current material/geometry basis and coherent gas are pinned together;
    /// later cold gas-inertia/torque/axial additions retain the original load.
    ///
    /// The model's producer must use this exact basis and exterior face ownership,
    /// independently of the interior cavity. This raw composition API cannot
    /// certify an authored model as geometry-derived. It adds no microphone
    /// filter: the existing compact observation retains its original scope.
    pub fn with_body_radiation_load(mut self, model: Model, maximum_poles: usize,
        gate: &CancelGate) -> Result<Self, ReciprocalStringError>
    {
        poll(gate)?;
        if self.source.body_radiation.is_some() || self.system.samples_rendered() != 0 {
            return Err(invalid("body exterior loading requires one original cold source binding"));
        }
        let (basis_identity, modes) = self.source.body_radiation_basis(&self.system, gate)?;
        if model.ports != modes.len() {
            return Err(invalid("body radiation requires every original body mode in its exact order"));
        }
        self.system = self.system.with_radiation_load(&model, &modes, maximum_poles, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.source.validate_body_memory_frequencies(&self.system, gate)?;
        self.source.validate_material_bridge_frequencies(&self.system, gate)?;
        self.source.body_radiation = Some(MaterialBodyRadiation { model, gas: self.source.gas,
            basis_identity, maximum_poles });
        poll(gate)?;
        Ok(self)
    }
}

fn hash_rows(hash: &mut DomainHasher, rows: &[Vec<f64>], gate: &CancelGate)
    -> Result<(), ReciprocalStringError>
{
    hash.update(&(rows.len() as u64).to_le_bytes());
    for row in rows {
        hash.update(&(row.len() as u64).to_le_bytes());
        for chunk in row.chunks(256) {
            poll(gate)?;
            for value in chunk { hash.update(&value.to_bits().to_le_bytes()); }
        }
    }
    Ok(())
}
fn hash_shell(hash: &mut DomainHasher, shell: &crate::shear_shell::ResolvedShearShell,
    gate: &CancelGate) -> Result<(), ReciprocalStringError>
{
    use fs_plate::shell::pressure::ShellFace;
    hash.update(&(shell.mesh().nodes.len() as u64).to_le_bytes());
    for node in &shell.mesh().nodes {
        poll(gate)?;
        for value in node { hash.update(&value.to_bits().to_le_bytes()); }
    }
    hash.update(&(shell.mesh().tris.len() as u64).to_le_bytes());
    for (i, tri) in shell.mesh().tris.iter().enumerate() {
        poll(gate)?;
        for index in tri { hash.update(&(*index as u64).to_le_bytes()); }
        let material = shell.material_at_element(i)
            .ok_or_else(|| invalid("exterior shell element has no original material source"))?;
        hash.update(material.identity().as_bytes());
    }
    for region in shell.regions() {
        poll(gate)?;
        hash.update(&region.thickness_m().to_bits().to_le_bytes());
        hash.update(&region.input().region.material_angle_rad.to_bits().to_le_bytes());
        for value in region.input().axis1_global { hash.update(&value.to_bits().to_le_bytes()); }
    }
    if let Some(surfaces) = shell.pressure_surfaces() {
        hash.update(&[1]);
        for selection in [&surfaces.exterior, &surfaces.cavity] {
            hash.update(&(selection.len() as u64).to_le_bytes());
            for facet in selection {
                poll(gate)?;
                hash.update(&(facet.triangle as u64).to_le_bytes());
                hash.update(&[match facet.side { ShellFace::Positive => 0, ShellFace::Negative => 1 }]);
            }
        }
    } else { hash.update(&[0]); }
    Ok(())
}
