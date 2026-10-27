//! One material-card declaration shared by DG jet and moving-valve sources.
//! The caller owns artifact I/O/cache; this module resolves the original
//! isotropic or principal-aligned orthotropic wall and retains its complete source.
//! No acoustic coefficient is inferred from a name and no fitted loss is introduced.

use super::{DgSourceError, Row, bad, checkpoint};
use crate::bernoulli_aperture::wall::material::source::WallMaterialLossSource;
use crate::bernoulli_aperture::wall::material::{AxialConstraint, LocalMaterialWall};
use fs_exec::CancelGate;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use std::{str::SplitAsciiWhitespace, sync::Arc};

/// Physical axial source partitions for the same uniform/regional DG wall owner.
pub mod regions;

/// The CLI bounds artifact bytes before decoding; direct callers share the cap.
pub const MAX_DG_WALL_PACK_BYTES: usize = 8 * 1024 * 1024;

enum WallConstitutiveInput {
    Isotropic(WallMaterialLossSource),
    Orthotropic {
        axes: [usize; 2],
        proportional_viscosity_time_property: Option<String>,
    },
}

/// One frozen, uniform solid around the original cylindrical main bore.
pub struct DgMaterialWallInput {
    line: usize,
    path: String,
    temperature_k: f64,
    angular_frequency_rad_s: f64,
    thickness_m: f64,
    axial: AxialConstraint,
    band_hz: [f64; 2],
    constitutive: WallConstitutiveInput,
}

/// Original source and resolved constitutive law, independent of DG resolution.
#[derive(Clone, Debug)]
pub struct DgMaterialWallRecord {
    path: String,
    pack: Arc<NormalizedMaterialCardPack>,
    wall: LocalMaterialWall,
}

impl DgMaterialWallRecord {
    /// Caller-relative source artifact, retained exactly as authored.
    pub fn path(&self) -> &str {
        &self.path
    }
    /// Full normalized source, shared with other users of the same artifact.
    pub fn pack(&self) -> &NormalizedMaterialCardPack {
        &self.pack
    }
    /// Original resolved material, thickness, axial restraint and loss spectrum.
    pub const fn wall(&self) -> &LocalMaterialWall {
        &self.wall
    }
}

impl DgMaterialWallInput {
    /// Read after the `dg_wall` keyword:
    /// `rigid`, `isotropic PACK T OMEGA H AXIAL FMIN FMAX LOSS`, or
    /// `orthotropic PACK T OMEGA H AXIAL FMIN FMAX HOOP_AXIS AXIAL_AXIS LOSS`.
    /// Axes are distinct one-based principal material directions in 1..=3,
    /// following the local cylindrical hoop/axial directions, not Cartesian grain.
    /// AXIAL is stress-free/strain-fixed. Isotropic LOSS retains elastic,
    /// kelvin-voigt or maxwell COUNT MAX_DT_OVER_TAU MAX_ANGULAR_STEP E TAU ... .
    /// Orthotropic LOSS is elastic or proportional-kelvin-voigt TIME_PROPERTY;
    /// the source time declares C_viscous=tau*C_elastic, never a fitted loss.
    /// General directional relaxation and off-axis normal/shear coupling refuse.
    pub fn read(
        line: usize,
        fields: &mut SplitAsciiWhitespace<'_>,
    ) -> Result<Option<Self>, DgSourceError> {
        let mut row = Row { line, fields };
        let orthotropic = match row.word()? {
            "rigid" => {
                row.finish()?;
                return Ok(None);
            }
            "isotropic" => false,
            "orthotropic" => true,
            _ => return Err(bad(line, "dg_wall requires rigid, isotropic or orthotropic material input")),
        };
        let path = row.word()?.to_owned();
        if path.len() > 512
            || path.contains('\\')
            || path.contains(':')
            || path.split('/').any(|part| part.is_empty() || part == "." || part == "..")
        {
            return Err(bad(line, "DG wall pack path must stay relative to the performance directory"));
        }
        let temperature_k = row.positive()?;
        let angular_frequency_rad_s = row.scalar()?;
        let thickness_m = row.positive()?;
        let axial = match row.word()? {
            "stress-free" => AxialConstraint::StressFree,
            "strain-fixed" => AxialConstraint::StrainFixed,
            _ => return Err(bad(line, "DG wall axial restraint must be stress-free or strain-fixed")),
        };
        let band_hz = [row.positive()?, row.positive()?];
        let constitutive = if orthotropic {
            let mut axes = [0; 2];
            for axis in &mut axes {
                *axis = row.word()?.parse::<usize>().ok().filter(|n| (1..=3).contains(n))
                    .ok_or_else(|| bad(line, "orthotropic wall axes must be distinct integers in 1..=3"))? - 1;
            }
            if axes[0] == axes[1] {
                return Err(bad(line, "orthotropic wall hoop and axial directions must differ"));
            }
            let proportional_viscosity_time_property = match row.word()? {
                "elastic" => None,
                "proportional-kelvin-voigt" => Some(row.word()?.to_owned()),
                _ => return Err(bad(line, "orthotropic wall loss requires elastic or proportional-kelvin-voigt TIME_PROPERTY")),
            };
            LocalMaterialWall::orthotropic_requirements(axes, proportional_viscosity_time_property.as_deref())
                .map_err(DgSourceError::Physics)?;
            WallConstitutiveInput::Orthotropic { axes, proportional_viscosity_time_property }
        } else {
            WallConstitutiveInput::Isotropic(
                WallMaterialLossSource::read(row.fields).map_err(|what| bad(line, what))?,
            )
        };
        row.finish()?;
        if angular_frequency_rad_s < 0.0 || band_hz[0] >= band_hz[1] {
            return Err(bad(line, "DG wall needs nonnegative source frequency and an ordered positive band"));
        }
        Ok(Some(Self {
            line, path, temperature_k, angular_frequency_rad_s, thickness_m,
            axial, band_hz, constitutive,
        }))
    }

    pub fn path(&self) -> &str {
        &self.path
    }
    pub const fn line(&self) -> usize {
        self.line
    }

    /// Resolve original source claims at one common point without defaults,
    /// geometry remapping, an effective isotropic card or second material history.
    pub fn resolve(
        self,
        pack: Arc<NormalizedMaterialCardPack>,
        gate: &CancelGate,
    ) -> Result<DgMaterialWallRecord, DgSourceError> {
        checkpoint(gate, self.line)?;
        if pack.to_bytes().len() > MAX_DG_WALL_PACK_BYTES {
            return Err(bad(self.line, "DG wall material artifact exceeds 8 MiB"));
        }
        let point = QueryPoint::new()
            .with("T", self.temperature_k)
            .and_then(|point| point.with("omega", self.angular_frequency_rad_s))
            .map_err(|_| bad(self.line, "invalid DG wall material query"))?;
        let requirements = match &self.constitutive {
            WallConstitutiveInput::Isotropic(loss) => loss.requirements()
                .map_err(|source| DgSourceError::Material { line: self.line, source })?,
            WallConstitutiveInput::Orthotropic { axes, proportional_viscosity_time_property } =>
                LocalMaterialWall::orthotropic_requirements(*axes, proportional_viscosity_time_property.as_deref())
                    .map_err(DgSourceError::Physics)?,
        };
        let material = resolve_material_state_point(pack.card(), &point, &requirements,
            MaterialPropertySelection::SingleClaimOnly)
            .map_err(|source| DgSourceError::Material { line: self.line, source })?;
        let wall = match self.constitutive {
            WallConstitutiveInput::Isotropic(loss) => {
                let loss = loss.resolve(&material, self.band_hz).map_err(DgSourceError::Physics)?;
                LocalMaterialWall::new(material, self.thickness_m, self.axial, loss, self.band_hz)
            }
            WallConstitutiveInput::Orthotropic { axes, proportional_viscosity_time_property } =>
                LocalMaterialWall::new_orthotropic(material, self.thickness_m, self.axial,
                    axes, proportional_viscosity_time_property.as_deref(), self.band_hz),
        }.map_err(DgSourceError::Physics)?;
        checkpoint(gate, self.line)?;
        Ok(DgMaterialWallRecord { path: self.path, pack, wall })
    }
}
