//! Original shell tables enter the same body/profile and cavity owners as the
//! standalone assembler. File selection changes data, never a geometry model.
use super::{Result, body, input, native_input::SourceFiles};
use fs_couple::cavity_modal::HelmholtzModalCavity;
use std::path::Path;

pub(super) struct LoadedBody {
    pub profile: Option<body::Profile>,
    pub geometry: Option<ShellGeometry>,
}

/// Physical dimensions returned by the existing shell-table converter. Its
/// profile uses normalized x/y coordinates and thickness ratios, so these must
/// travel together; retaining old performance dimensions would rescale a swap.
pub(super) struct ShellGeometry {
    length_m: f64,
    width_m: f64,
    thickness_m: f64,
    cavity: Option<[f64; 3]>,
}

impl ShellGeometry {
    pub fn bind(&self, spec: &mut input::Spec) -> Result<()> {
        // Rebuild through the original cavity owner rather than copying old
        // inertia/loss coefficients onto a new volume or neck. An explicit
        // cavity-none comparison remains cavity-free.
        let cavity = spec
            .cavity
            .as_ref()
            .map(|original| {
                let [volume_m3, neck_radius_m, neck_length_m] = self.cavity.ok_or(
                    "the selected shell table has no cavity geometry for this performance",
                )?;
                let mut geometry = *original.cavity.geometry();
                geometry.volume_m3 = volume_m3;
                geometry.neck_radius_m = neck_radius_m;
                geometry.neck_length_m = neck_length_m;
                let mut rebound = original.clone();
                rebound.cavity = HelmholtzModalCavity::new(
                    geometry,
                    *original.cavity.gas(),
                    original.cavity.loss(),
                )?;
                Ok::<_, Box<dyn std::error::Error>>(rebound)
            })
            .transpose()?;
        spec.body.length = self.length_m;
        spec.body.width = self.width_m;
        spec.body.thickness = self.thickness_m;
        spec.cavity = cavity;
        Ok(())
    }
}

pub(super) fn load(
    profile: Option<&Path>,
    table: Option<&Path>,
    options: Option<&Path>,
    sources: &SourceFiles,
) -> Result<LoadedBody> {
    if table.is_some() != options.is_some() {
        return Err("--shell-table and --assembly-options must be supplied together".into());
    }
    if profile.is_some() && table.is_some() {
        return Err(
            "select --body-profile or --shell-table with --assembly-options, not both".into(),
        );
    }
    if let (Some(table), Some(options)) = (table, options) {
        let options = body::ShellTableOptions::parse(std::str::from_utf8(
            &sources.read_option(options, body::MAX_BYTES)?,
        )?)?;
        let converted = options.convert(std::str::from_utf8(
            &sources.read_option(table, body::MAX_SHELL_BYTES)?,
        )?)?;
        let profile = body::Profile::parse(converted.profile.as_bytes())?;
        return Ok(LoadedBody {
            profile: Some(profile),
            geometry: Some(ShellGeometry {
                length_m: converted.length_m,
                width_m: converted.width_m,
                thickness_m: options.reference_thickness_m(),
                cavity: converted.cavity,
            }),
        });
    }
    Ok(LoadedBody {
        profile: profile
            .map(|path| body::Profile::parse(&sources.read_option(path, body::MAX_SHELL_BYTES)?))
            .transpose()?,
        geometry: None,
    })
}
