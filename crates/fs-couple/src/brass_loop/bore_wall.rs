//! The `bore_wall` row of a lip performance file: the breathing wall of a brass bore.
//!
//! `bore_wall rigid` is the row absent: the bore's walls are rigid and the audio is bit-identical (the input bytes, and so
//! the input hash in the receipt, differ).
//! `bore_wall thin-brass <thickness_m> <quality>` gives every cell of the nonlinear bore a radial wall
//! oscillator ([`fs_duct::nonlinear::bore::NonlinearBore::bind_shell`]): the mass and stiffness of the wall at
//! each cell's radius come from the landed membrane/Lame reduction
//! ([`LocalMaterialWall::base_patch`], membrane at thickness/radius up to 0.05, Lame's thick cylinder above it) for
//! UNS C26000 cartridge brass, the cone factor `cos phi` of the classical thin cone membrane is applied to the
//! membrane stiffness by the bore, and `quality` is the quality factor of the breathing mode of every cell.
//! The law of the wall is the classical thin cone membrane (hoop resultant `p r / cos phi`, so
//! `s = p r^2 / (E t cos phi)`); Kausel, Zietlow and Moore (JASA 128 (2010) 3161) print `1 / cos^3 phi` for the
//! static radial response, a factor `1 / cos^2 phi` more, which is not used. Bending and axial shell modes, a
//! rim bead, the outside air's added mass and a non-uniform thickness are not carried.
//!
//! The thickness and the quality factor are the card's numbers. A card that carries them names, in its
//! `bore_wall_source` row, whose instrument and whose measurement they are; nothing here supplies either.

use crate::bernoulli_aperture::wall::material::source::WallMaterialLossSource;
use crate::bernoulli_aperture::wall::material::{AxialConstraint, LocalMaterialWall, WallMaterialLoss};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure};

/// Cartridge brass 70 percent (UNS C26000 class), density at 68 F: `data/matdb/seed-v1/brass-cartridge-c26000-mil-hdbk-698a`.
pub const BRASS_DENSITY_KG_M3: f64 = 8525.4;
/// Its modulus of elasticity [Pa] (110316 MPa in the same row).
pub const BRASS_YOUNG_PA: f64 = 110_316.0e6;
/// Its Poisson ratio as that seed row prints it (derived from the modulus of elasticity and the shear modulus).
pub const BRASS_POISSON_RATIO: f64 = 0.333;
/// The band over which the wall's reduction is evaluated (the membrane law has no frequency dependence).
const BAND_HZ: [f64; 2] = [20.0, 20_000.0];
/// The seed row's own temperature, 68 F.
const KELVIN: f64 = 293.0;

/// The parsed `bore_wall thin-brass` row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoreWallRecord {
    /// Wall thickness [m].
    pub thickness_m: f64,
    /// Quality factor of every cell's breathing mode.
    pub quality: f64,
}

/// The parsed `bore_wall shell-brass` row: the bell as a shell of brass ([`super::bell_shell`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BellShellRecord {
    /// Wall thickness [m].
    pub thickness_m: f64,
    /// Loss factor of every mode (modal damping `eta omega`).
    pub loss_factor: f64,
    /// Modes up to this frequency [Hz] are retained.
    pub max_frequency_hz: f64,
    /// Points round the bell's mesh (even, `BELL_MIN_AZIMUTHS` to `BELL_MAX_AZIMUTHS`). The mesh is a polygon of flat facets, which is
    /// stiffer than the circle for a breathing that varies along the axis by an error that falls as the square of this count: a
    /// convergence quantity of the card, like the cut-off.
    pub azimuths: usize,
}

/// The fewest points round the bell's mesh: the shell owner admits a director only within 0.9 of the facet normals at its node, which
/// the 16-gon (neighbouring normals at cos(2 pi / 16) = 0.924) meets and the 12-gon (0.866) does not.
pub const BELL_MIN_AZIMUTHS: usize = 16;
/// The most points round the bell's mesh.
pub const BELL_MAX_AZIMUTHS: usize = 128;

/// The brass wall's reduction for the shell: surface density and stiffness at a cell radius.
#[derive(Clone, Debug)]
pub struct ThinBrassWall {
    #[cfg_attr(not(feature = "nonlinear-brass"), allow(dead_code))]
    wall: LocalMaterialWall,
}

impl ThinBrassWall {
    /// The wall of thickness `thickness_m` of the seed's cartridge brass, stress-free in the axial direction.
    ///
    /// # Errors
    /// A thickness the wall owner refuses, or the card assembly (which cannot fail for these constants).
    pub fn new(thickness_m: f64) -> Result<Self, String> {
        if !(thickness_m.is_finite() && thickness_m > 0.0) {
            return Err("bore_wall needs a positive finite thickness".into());
        }
        let material = resolved_brass()?;
        let wall = LocalMaterialWall::new(
            material,
            thickness_m,
            AxialConstraint::StressFree,
            WallMaterialLoss::Elastic,
            BAND_HZ,
        )
        .map_err(|e| format!("bore_wall wall refused: {e}"))?;
        Ok(Self { wall })
    }
}

fn resolved_brass() -> Result<ResolvedMaterialStatePoint, String> {
    let refuse = |what: &str, e: &dyn core::fmt::Debug| format!("bore_wall brass card refused at {what}: {e:?}");
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, BRASS_DENSITY_KG_M3),
        ("young_modulus", Pressure::DIMS, BRASS_YOUNG_PA),
        ("poisson_ratio", Dims::NONE, BRASS_POISSON_RATIO),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", KELVIN, KELVIN),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "MIL-HDBK-698A (1971) Table XV row 260, cartridge brass 70 percent, via data/matdb/seed-v1/brass-cartridge-c26000-mil-hdbk-698a".into(),
                    license: "Work-of-the-US-Government-Public-Use-Permitted".into(),
                    artifact: None,
                },
            })
            .map_err(|e| refuse(name, &e))?;
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "brass-cartridge-c26000".into(),
            phase: "solid".into(),
            process: "mil-hdbk-698a-table-xv".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .map_err(|e| refuse("assembly", &e))?;
    let requirements = WallMaterialLossSource::Elastic
        .requirements()
        .map_err(|e| refuse("requirements", &e))?;
    let point = QueryPoint::new()
        .with("T", KELVIN)
        .and_then(|p| p.with("omega", 1000.0))
        .map_err(|e| refuse("query", &e))?;
    resolve_material_state_point(&card, &point, &requirements, MaterialPropertySelection::SingleClaimOnly)
        .map_err(|e| refuse("resolution", &e))
}

#[cfg(feature = "nonlinear-brass")]
impl fs_duct::nonlinear::bore::ShellLaw for ThinBrassWall {
    fn strip(
        &self,
        radius_m: f64,
    ) -> Result<fs_duct::nonlinear::bore::ShellStrip, fs_duct::nonlinear::NonlinearDuctError> {
        let patch = self.wall.base_patch(radius_m, 1.0).map_err(|_| {
            fs_duct::nonlinear::NonlinearDuctError::UnsupportedGeometry(
                "the brass wall law refused this bore radius",
            )
        })?;
        Ok(fs_duct::nonlinear::bore::ShellStrip {
            surface_density_kg_m2: patch.wall.surface_density,
            stiffness_per_area_pa_m: patch.wall.stiffness_per_area,
            membrane: self.wall.runs_membrane_at(radius_m),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The three constants are the committed seed row's numbers; a change of the data file fails here.
    #[test]
    fn the_brass_constants_are_the_seed_rows_numbers() {
        let seed = include_str!(
            "../../../../data/matdb/seed-v1/brass-cartridge-c26000-mil-hdbk-698a/properties.tsv"
        );
        let scalar = |name: &str| -> (f64, String) {
            let line = seed
                .lines()
                .find(|l| l.starts_with("scalar\t") && l.split('\t').nth(1) == Some(name))
                .unwrap_or_else(|| panic!("seed row {name}"));
            let fields: Vec<&str> = line.split('\t').collect();
            (fields[4].parse().unwrap(), fields[5].to_string())
        };
        assert_eq!(scalar("density"), (BRASS_DENSITY_KG_M3, "kg/m3".to_string()));
        assert_eq!(scalar("young-modulus"), (BRASS_YOUNG_PA / 1.0e6, "MPa".to_string()));
        assert_eq!(scalar("poisson-ratio-derived"), (BRASS_POISSON_RATIO, "1".to_string()));
    }

    #[test]
    fn the_strip_is_the_membrane_expression_at_the_gate_and_lames_above_it() {
        let t = 0.5e-3;
        let wall = ThinBrassWall::new(t).unwrap();
        // Thin: sigma = rho t, K = E t / r^2 (stress-free axial restraint).
        let r = 30.0e-3;
        let patch = wall.wall.base_patch(r, 1.0).unwrap();
        assert!(wall.wall.runs_membrane_at(r));
        assert!((patch.wall.surface_density / (BRASS_DENSITY_KG_M3 * t) - 1.0).abs() <= 1e-12);
        assert!((patch.wall.stiffness_per_area / (BRASS_YOUNG_PA * t / (r * r)) - 1.0).abs() <= 1e-12);
        // A lead pipe, t / r = 0.083: the thick-wall branch. Lame's cylinder is softer than the membrane by about
        // (1/2 + nu) t / r = 0.069 here (the landed owner prints 4.06 % at the gate, t / r = 0.05).
        let r = 6.0e-3;
        assert!(!wall.wall.runs_membrane_at(r));
        let thick = wall.wall.base_patch(r, 1.0).unwrap();
        let membrane = BRASS_YOUNG_PA * t / (r * r);
        let ratio = thick.wall.stiffness_per_area / membrane;
        assert!((0.90..0.96).contains(&ratio), "{ratio}");
        // Refusals by name.
        assert!(ThinBrassWall::new(0.0).unwrap_err().contains("positive finite thickness"));
        assert!(ThinBrassWall::new(f64::NAN).unwrap_err().contains("positive finite thickness"));
    }
}
