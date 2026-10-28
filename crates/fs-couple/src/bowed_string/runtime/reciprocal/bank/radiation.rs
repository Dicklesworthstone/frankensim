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
    source_body_coordinates: Vec<usize>,
    // Present only when the original material-shell producer attached this load.
    shell_bem: Option<bem::ShellBemReport>,
    // Present only when the original retained plate produced this load.
    plate_rayleigh: Option<rayleigh::PlateRayleighReport>,
}
impl MaterialBodyRadiation {
    pub const fn model(&self) -> &Model { &self.model }
    pub const fn gas(&self) -> &GasState { &self.gas }
    pub const fn basis_identity(&self) -> ContentHash { self.basis_identity }
    /// Complete source order relative to the original body component: all
    /// structural modes, then explicitly declared cavity-neck coordinates.
    /// Interior standing-wave coordinates are never radiation sources.
    pub fn source_body_coordinates(&self) -> &[usize] { &self.source_body_coordinates }

    /// Preparation report of the source-derived closed-shell load actually
    /// attached to this body. Raw and planar loads do not carry this binding.
    /// It survives cold bank rebuilds with the original basis and gas checks.
    pub const fn material_shell_bem(&self) -> Option<&bem::ShellBemReport> {
        self.shell_bem.as_ref()
    }

    /// Original planar surface, retained basis and gas of the computed load.
    /// A raw supplied impedance cannot stand in for this geometry binding.
    pub const fn material_plate_rayleigh(&self) -> Option<&rayleigh::PlateRayleighReport> {
        self.plate_rayleigh.as_ref()
    }
}

impl MaterialStringBankSource {
    /// Complete exterior source retained through scheduled and simultaneous bows.
    pub const fn body_radiation(&self) -> Option<&MaterialBodyRadiation> { self.body_radiation.as_ref() }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn restore_body_radiation(
        &self, system: FrictionModalSystem, gate: &CancelGate,
    ) -> Result<FrictionModalSystem, ReciprocalStringError> {
        let Some(air) = &self.body_radiation else { return Ok(system); };
        let (identity, body_modes) = self.body_radiation_basis(&system, gate)?;
        let coordinates = self.radiation_body_coordinates(body_modes.len())?;
        if air.gas != self.gas || air.basis_identity != identity {
            return Err(invalid("exterior load no longer matches the original material body basis and gas"));
        }
        if coordinates != air.source_body_coordinates
            && !self.skin_prefix_coordinates(air, body_modes.len())
        {
            return Err(invalid("exterior load lost its original complete body/neck source map"));
        }
        let start = body_modes[0];
        let modes: Vec<_> = air.source_body_coordinates.iter().map(|coordinate|start+coordinate).collect();
        let system = system.with_radiation_load(&air.model, &modes, air.maximum_poles, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.validate_body_memory_frequencies(&system, gate)?;
        self.validate_material_bridge_frequencies(&system, gate)?;
        Ok(system)
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn radiation_body_coordinates(
        &self, structural_modes: usize,
    ) -> Result<Vec<usize>, ReciprocalStringError> {
        if structural_modes == 0 { return Err(invalid("exterior source requires a nonempty original body")); }
        if let Some(cavity) = &self.spatial_cavity {
            if cavity.coupling().structural_modes() != structural_modes {
                return Err(invalid("exterior source lost its original cavity/body basis"));
            }
            if cavity.opening_count() != 0 {
                return cavity.exterior_body_coordinates().map_err(|_|
                    invalid("exterior source requires every original internal-only cavity neck"));
            }
        }
        // An exterior-owned lumped mouth radiates through the one neck coordinate,
        // the first of the component after the body.
        if self.material_cavity.as_ref().is_some_and(|c|
            c.cavity.loss() == crate::cavity_modal::CavityNeckLoss::ExteriorMouth)
        {
            return Ok((0..=structural_modes).collect());
        }
        Ok((0..structural_modes).collect())
    }

    /// A skin load's coordinates: an ascending prefix of the structural modes,
    /// then the exterior-owned neck when there is one. Nothing else is admitted.
    fn skin_prefix_coordinates(&self, air: &MaterialBodyRadiation, structural_modes: usize) -> bool {
        if !air.shell_bem.as_ref().is_some_and(|report| report.skin_identity.is_some()) {
            return false;
        }
        let neck = self.material_cavity.as_ref().is_some_and(|c|
            c.cavity.loss() == crate::cavity_modal::CavityNeckLoss::ExteriorMouth);
        let stored = &air.source_body_coordinates;
        let loaded = stored.len() - usize::from(neck);
        (1..=structural_modes).contains(&loaded)
            && stored[..loaded].iter().copied().eq(0..loaded)
            && (!neck || stored[loaded] == structural_modes)
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn body_radiation_basis(&self, system: &FrictionModalSystem, gate: &CancelGate)
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
    /// followed by every explicitly declared compact cavity-neck velocity in
    /// source order. Neck columns use their original mass-normalized coordinate;
    /// a physical mouth trace therefore equals volume_weight / mouth_area.
    /// Air reacts on the body and hence
    /// the bridges, every string, bow traction and accepted pressure acceleration.
    /// The current material/geometry basis and coherent gas are pinned together;
    /// later cold gas-inertia/torque/axial additions retain the original load.
    ///
    /// The model's producer must use this exact basis and exterior face ownership,
    /// independently of the interior cavity. This raw composition API cannot
    /// certify an authored model as geometry-derived. It adds no microphone
    /// filter: the existing compact observation retains its original scope.
    pub fn with_body_radiation_load(self, model: Model, maximum_poles: usize,
        gate: &CancelGate) -> Result<Self, ReciprocalStringError>
    {
        self.attach_body_radiation(model, maximum_poles, None, gate)
    }

    /// The shared attachment. `coordinates` is a skin load's admitted prefix
    /// map; `None` is the complete body/neck map.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn attach_body_radiation(
        mut self, model: Model, maximum_poles: usize, coordinates: Option<Vec<usize>>,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if self.source.body_radiation.is_some() || self.system.samples_rendered() != 0 {
            return Err(invalid("body exterior loading requires one original cold source binding"));
        }
        let (basis_identity, body_modes) = self.source.body_radiation_basis(&self.system, gate)?;
        let source_body_coordinates = match coordinates {
            Some(coordinates) => coordinates,
            None => self.source.radiation_body_coordinates(body_modes.len())?,
        };
        let modes: Vec<_> = source_body_coordinates.iter().map(|coordinate|body_modes[0]+coordinate).collect();
        if model.ports != modes.len() {
            return Err(invalid("radiation requires every original body mode and actual mouth in source order"));
        }
        self.system = self.system.with_radiation_load(&model, &modes, maximum_poles, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.source.validate_body_memory_frequencies(&self.system, gate)?;
        self.source.validate_material_bridge_frequencies(&self.system, gate)?;
        self.source.body_radiation = Some(MaterialBodyRadiation { model, gas: self.source.gas,
            basis_identity, maximum_poles, source_body_coordinates, shell_bem: None, plate_rayleigh: None });
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
