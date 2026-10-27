//! One material-card declaration shared by DG jet and moving-valve sources.
//! The caller owns artifact I/O/cache; this module resolves the original
//! isotropic wall and retains its complete source. No acoustic coefficients
//! are inferred from a material name and no fitted wall loss is introduced.

use super::{DgSourceError, Row, bad, checkpoint};
use crate::bernoulli_aperture::wall::material::source::WallMaterialLossSource;
use crate::bernoulli_aperture::wall::material::{AxialConstraint, LocalMaterialWall};
use fs_exec::CancelGate;
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{MaterialPropertySelection, resolve_material_state_point};
use std::{str::SplitAsciiWhitespace, sync::Arc};

/// The CLI bounds artifact bytes before decoding; direct callers share the cap.
pub const MAX_DG_WALL_PACK_BYTES: usize = 8 * 1024 * 1024;

/// One frozen, uniform solid around the original cylindrical main bore.
pub struct DgMaterialWallInput {
    line: usize,
    path: String,
    temperature_k: f64,
    angular_frequency_rad_s: f64,
    thickness_m: f64,
    axial: AxialConstraint,
    band_hz: [f64; 2],
    loss: WallMaterialLossSource,
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
    /// `rigid`, or `isotropic PACK T OMEGA H AXIAL FMIN FMAX LOSS`.
    /// AXIAL is stress-free/strain-fixed. LOSS is the existing elastic,
    /// kelvin-voigt or maxwell COUNT MAX_DT_OVER_TAU MAX_ANGULAR_STEP E TAU ...
    /// source syntax. The last two limits retain source identity; the exact DG
    /// propagator does not impose a midpoint wall clock.
    pub fn read(
        line: usize,
        fields: &mut SplitAsciiWhitespace<'_>,
    ) -> Result<Option<Self>, DgSourceError> {
        let mut row = Row { line, fields };
        match row.word()? {
            "rigid" => {
                row.finish()?;
                return Ok(None);
            }
            "isotropic" => {}
            _ => {
                return Err(bad(
                    line,
                    "dg_wall requires rigid or an isotropic material source",
                ));
            }
        }
        let path = row.word()?.to_owned();
        if path.len() > 512
            || path.contains('\\')
            || path.contains(':')
            || path
                .split('/')
                .any(|part| part.is_empty() || part == "." || part == "..")
        {
            return Err(bad(
                line,
                "DG wall pack path must stay relative to the performance directory",
            ));
        }
        let temperature_k = row.positive()?;
        let angular_frequency_rad_s = row.scalar()?;
        let thickness_m = row.positive()?;
        let axial = match row.word()? {
            "stress-free" => AxialConstraint::StressFree,
            "strain-fixed" => AxialConstraint::StrainFixed,
            _ => {
                return Err(bad(
                    line,
                    "DG wall axial restraint must be stress-free or strain-fixed",
                ));
            }
        };
        let band_hz = [row.positive()?, row.positive()?];
        let loss = WallMaterialLossSource::read(row.fields).map_err(|what| bad(line, what))?;
        row.finish()?;
        if angular_frequency_rad_s < 0.0 || band_hz[0] >= band_hz[1] {
            return Err(bad(
                line,
                "DG wall needs nonnegative source frequency and an ordered positive band",
            ));
        }
        Ok(Some(Self {
            line,
            path,
            temperature_k,
            angular_frequency_rad_s,
            thickness_m,
            axial,
            band_hz,
            loss,
        }))
    }

    pub fn path(&self) -> &str {
        &self.path
    }
    pub const fn line(&self) -> usize {
        self.line
    }

    /// Resolve the original source claims, with no defaults for missing material
    /// properties, no geometry remapping and no second material evolution.
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
        let material = (|| {
            let requirements = self.loss.requirements()?;
            resolve_material_state_point(
                pack.card(),
                &point,
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
        })()
        .map_err(|source| DgSourceError::Material {
            line: self.line,
            source,
        })?;
        let loss = self
            .loss
            .resolve(&material, self.band_hz)
            .map_err(DgSourceError::Physics)?;
        let wall =
            LocalMaterialWall::new(material, self.thickness_m, self.axial, loss, self.band_hz)
                .map_err(DgSourceError::Physics)?;
        checkpoint(gate, self.line)?;
        Ok(DgMaterialWallRecord {
            path: self.path,
            pack,
            wall,
        })
    }
}
