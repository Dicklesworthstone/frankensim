//! Native molecular-bore input, using the original material parser and field.
//!
//! Three records precede `compile_limits`:
//! `molecular <fs-material relaxation declaration>`,
//! `molecular_resolution CONE_SLICES MAX_BLOCKS`, and
//! `molecular_wall none|diffusive POLES OMEGA_MIN OMEGA_MAX MAX_REL_ERROR`.
//! Angular-frequency bounds are in rad/s, including the line's near-DC probe
//! and Nyquist. Wall and cone resolution are independent of the FIR ladder.
//! All source/domain refusals remain with their original owners. A declared
//! reservoir is authored constitutive data, not a measurement certificate.
//! A preceding continuous gas_field or gas_steady block supplies the complete
//! molecular medium. Its exact inlet, rather than an unrelated ambient record,
//! admits the law; the retained declaration is then evaluated along the field.

use super::{GasRecord, Reader, ReedPerformanceError, bad};
use crate::acoustic_realize::map_drive;
use crate::driving_point::molecular::MolecularDrivingPoint;
use fs_duct::compose::{ComposeOptions, MAX_COMPOSITION_BLOCKS};
use fs_duct::molecular::MolecularWallResolution;
use fs_duct::{ConeSlices, Duct};
use fs_material::gas::GasState;
use fs_material::gas_field::AxialGasField;
use fs_material::relaxation::{
    RelaxationAbsorption, RelaxationState, RelaxingMedium,
    input::{medium_from_record, read_relaxation},
};

pub(super) struct MolecularInput {
    law: RelaxationAbsorption,
    state: RelaxationState,
    cones: ConeSlices,
    maximum_blocks: usize,
    wall: Option<MolecularWallResolution>,
}

pub(super) fn read(
    reader: &mut Reader<'_>,
    ambient: &GasRecord,
    field: Option<&AxialGasField>,
) -> Result<MolecularInput, ReedPerformanceError> {
    let mut row = reader.row("molecular")?;
    let line = row.line;
    let law = read_relaxation(&mut row.words)
        .map_err(|source| ReedPerformanceError::Gas { line, source })?;
    row.finish()?;
    let medium = match field {
        Some(field) => field
            .state_at(0.0)
            .and_then(|state| RelaxingMedium::from_thermochemical(&state)),
        None => medium_from_record(ambient),
    }
    .map_err(|source| ReedPerformanceError::Gas { line, source })?;
    let state = law
        .evaluate(&medium)
        .map_err(|source| ReedPerformanceError::Gas { line, source })?;

    let mut row = reader.row("molecular_resolution")?;
    let slices = row.count(fs_duct::profile::MAX_CONE_SLICES)?;
    let maximum_blocks = row.count(MAX_COMPOSITION_BLOCKS)?;
    if slices == 0 || maximum_blocks == 0 {
        return Err(bad(
            row.line,
            "molecular cone slices and composition cap must be positive",
        ));
    }
    row.finish()?;

    let mut row = reader.row("molecular_wall")?;
    let wall = match row.word()? {
        "none" => None,
        "diffusive" => {
            let poles_per_field = row.count(fs_phs::MAX_ZK_DIFFUSIVE_POLES)?;
            let band_rad_s = [row.positive()?, row.positive()?];
            let maximum_relative_error = row.positive()?;
            if poles_per_field == 0 {
                return Err(bad(
                    row.line,
                    "molecular diffusion requires a positive pole count",
                ));
            }
            Some(MolecularWallResolution {
                poles_per_field,
                band_rad_s,
                maximum_relative_error,
            })
        }
        _ => return Err(bad(row.line, "molecular_wall requires none or diffusive")),
    };
    row.finish()?;
    Ok(MolecularInput {
        law,
        state,
        cones: ConeSlices::Fixed(slices),
        maximum_blocks,
        wall,
    })
}

impl MolecularInput {
    pub(super) fn bind(
        self,
        geometry: &Duct,
        exterior: GasState,
    ) -> Result<MolecularDrivingPoint, ReedPerformanceError> {
        MolecularDrivingPoint::new(
            geometry,
            self.state,
            self.wall,
            self.cones,
            ComposeOptions {
                maximum_blocks: self.maximum_blocks,
                ..ComposeOptions::default()
            },
            exterior,
        )
        .map_err(|error| ReedPerformanceError::Model(map_drive(error)))
    }

    pub(super) fn bind_field(
        self,
        geometry: &Duct,
        field: &AxialGasField,
        steps_per_span: fs_duct::RiccatiSteps,
        exterior: GasState,
    ) -> Result<MolecularDrivingPoint, ReedPerformanceError> {
        MolecularDrivingPoint::from_gas_field(
            geometry,
            field,
            &self.law,
            self.wall,
            steps_per_span,
            self.cones,
            ComposeOptions {
                maximum_blocks: self.maximum_blocks,
                ..ComposeOptions::default()
            },
            exterior,
        )
        .map_err(|error| ReedPerformanceError::Model(map_drive(error)))
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod field_tests;
