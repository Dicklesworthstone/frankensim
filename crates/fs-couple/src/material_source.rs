//! Pack-backed material records for file-driven render inputs (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.32`, part 3).
//!
//! A render input names its material either by a compiled material-card pack
//! and the state that card was compiled at, or by bare numbers carrying an
//! explicit `fixture LABEL`. This module parses the shared record
//! `PACK STATE [AXIS=VALUE ...]`, checks the named state against the card, and
//! resolves the consumer's properties through `fs_material::state_point`. A
//! missing, ambiguous or out-of-domain property refuses by name; nothing
//! defaults.
//!
//! It never reads files: fs-couple is wasm-reachable, so the caller supplies
//! each named pack (the `music_render` CLI reads it relative to the
//! performance file). Naming the state is mandatory because one material can
//! be played in several states: a cane reed is played water-soaked, and its
//! dry modulus is about 1.6 times stiffer.

use fs_matdb::{MaterialCard, NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ORTHOTROPIC_POISSON_RATIO_PROPERTIES,
    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_plate::PlateSection;
use fs_qty::{Density, Dims, Pressure};

use crate::thin_plate::{PlateMaterialModel, plate_elasticity, plate_property};

/// Upper bound on a record's query axes.
pub const MAX_RECORD_AXES: usize = 8;
/// Upper bound on a pack path, state or fixture label.
pub const MAX_RECORD_TOKEN_BYTES: usize = 256;

/// A malformed record, or a pack that cannot supply what the record names.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MaterialIngressError(pub String);

impl core::fmt::Display for MaterialIngressError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for MaterialIngressError {}

fn refuse(what: impl Into<String>) -> MaterialIngressError {
    MaterialIngressError(what.into())
}

/// A token safe to echo into JSON provenance: printable ASCII without
/// quotes, backslashes or `=`, within the byte cap.
fn plain_token(token: &str, what: &str) -> Result<(), MaterialIngressError> {
    if token.is_empty()
        || token.len() > MAX_RECORD_TOKEN_BYTES
        || !token
            .bytes()
            .all(|b| b.is_ascii_graphic() && !matches!(b, b'"' | b'\\' | b'='))
    {
        return Err(refuse(format!(
            "{what} must be 1..={MAX_RECORD_TOKEN_BYTES} printable ASCII bytes without quotes, backslashes or '='"
        )));
    }
    Ok(())
}

/// `PACK STATE [AXIS=VALUE ...]`: a compiled material-card pack, the state its
/// card was compiled at, and the physical query point.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialRecord {
    /// Pack path, relative to the input that names it.
    pub pack: String,
    /// The card's material state (`MaterialStateId::process`), named explicitly.
    pub state: String,
    /// Query coordinates in the card's own normalized units, in record order.
    pub axes: Vec<(String, f64)>,
}

impl MaterialRecord {
    /// Parse the record's tokens.
    ///
    /// # Errors
    /// A missing pack or state, a path leaving its directory, an unsafe token,
    /// a malformed, repeated or nonfinite axis, or too many axes refuse.
    pub fn parse(tokens: &[&str]) -> Result<Self, MaterialIngressError> {
        let [pack, state, axes @ ..] = tokens else {
            return Err(refuse(
                "a material record needs PACK STATE [AXIS=VALUE ...]",
            ));
        };
        plain_token(pack, "a material pack path")?;
        if pack
            .split('/')
            .any(|part| part.is_empty() || part == "." || part == "..")
            || pack.contains(':')
        {
            return Err(refuse(
                "a material pack path must stay relative to its input's directory",
            ));
        }
        plain_token(state, "a material state")?;
        if axes.len() > MAX_RECORD_AXES {
            return Err(refuse(format!(
                "a material record takes at most {MAX_RECORD_AXES} query axes"
            )));
        }
        let mut parsed: Vec<(String, f64)> = Vec::with_capacity(axes.len());
        for token in axes {
            let Some((axis, value)) = token.split_once('=') else {
                return Err(refuse(format!("query axis '{token}' must be AXIS=VALUE")));
            };
            plain_token(axis, "a query axis")?;
            let value: f64 = value
                .parse()
                .map_err(|_| refuse(format!("query axis '{axis}' needs a numeric value")))?;
            if !value.is_finite() || parsed.iter().any(|(seen, _)| seen == axis) {
                return Err(refuse(format!(
                    "query axis '{axis}' must be finite and appear once"
                )));
            }
            parsed.push((axis.to_string(), value));
        }
        Ok(Self {
            pack: (*pack).to_string(),
            state: (*state).to_string(),
            axes: parsed,
        })
    }

    /// The pack's card, provided it was compiled at the named state.
    ///
    /// # Errors
    /// A card compiled at any other state refuses, naming both.
    pub fn card<'a>(
        &self,
        pack: &'a NormalizedMaterialCardPack,
    ) -> Result<&'a MaterialCard, MaterialIngressError> {
        let card = pack.card();
        if card.id().process != self.state {
            return Err(refuse(format!(
                "pack {} was compiled at state '{}', not the named state '{}'",
                self.pack,
                card.id().process,
                self.state
            )));
        }
        Ok(card)
    }

    /// The query point, each axis in the quantity schema the card's claims
    /// declare for it (untyped where no claim types it).
    ///
    /// # Errors
    /// A coordinate the point refuses.
    pub fn query_point(&self, card: &MaterialCard) -> Result<QueryPoint, MaterialIngressError> {
        let mut point = QueryPoint::new();
        for (axis, value) in &self.axes {
            let quantity = card.claims().claims_ordered().find_map(|(_, claim)| {
                claim.validity.axis_quantities().get(axis.as_str()).copied()
            });
            point = match quantity {
                Some(quantity) => point.with_quantity(axis.clone(), quantity, *value),
                None => point.with(axis.clone(), *value),
            }
            .map_err(|e| refuse(format!("query axis '{axis}': {e:?}")))?;
        }
        Ok(point)
    }

    /// `pack PATH HASH state STATE at AXIS=VALUE,...` for render provenance.
    #[must_use]
    pub fn provenance(&self, pack: &NormalizedMaterialCardPack) -> String {
        let axes = self
            .axes
            .iter()
            .map(|(axis, value)| format!("{axis}={value:?}"))
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "pack {} {} state {} at {}",
            self.pack,
            pack.content_hash().to_hex(),
            self.state,
            if axes.is_empty() { "-" } else { &axes }
        )
    }
}

/// Optional held environment after a native material record's existing fields:
/// `environment COUNT AXIS VALUE ...`. The source card owns quantity semantics;
/// this only admits bounded coordinates before any material artifact is loaded.
pub(crate) fn read_fixed_environment(
    fields: &mut std::str::SplitAsciiWhitespace<'_>,
) -> Result<Vec<(String, f64)>, MaterialIngressError> {
    let Some(selector) = fields.next() else { return Ok(Vec::new()); };
    if selector != "environment" {
        return Err(refuse("material state suffix requires environment COUNT AXIS VALUE ..."));
    }
    let count = fields.next().and_then(|value| value.parse::<usize>().ok())
        .filter(|count| (1..=fs_material::state_point::MAX_FIXED_ENVIRONMENT_COORDINATES).contains(count))
        .ok_or_else(|| refuse("material environment needs 1..=8 explicit coordinates"))?;
    let mut coordinates: Vec<(String, f64)> = Vec::with_capacity(count);
    for _ in 0..count {
        let axis = fields.next().ok_or_else(|| refuse("missing material environment axis"))?;
        plain_token(axis, "a material environment axis")?;
        let value = fields.next().and_then(|value| value.parse::<f64>().ok())
            .filter(|value| value.is_finite())
            .ok_or_else(|| refuse("material environment coordinates need finite numeric values"))?;
        if axis == "T" || axis == "omega" || coordinates.iter().any(|(prior, _)| prior == axis) {
            return Err(refuse("material environment coordinates must be distinct and cannot replace T or omega"));
        }
        coordinates.push((axis.to_owned(), value));
    }
    Ok(coordinates)
}

/// The label of a `fixture LABEL` suffix that follows bare material numbers.
///
/// # Errors
/// Anything but exactly `fixture LABEL`, or an unsafe label, refuses.
pub fn fixture_label<'a>(tokens: &[&'a str]) -> Result<&'a str, MaterialIngressError> {
    match tokens {
        ["fixture", label] => {
            plain_token(label, "a fixture label")?;
            Ok(label)
        }
        _ => Err(refuse(
            "bare material numbers need an explicit 'fixture LABEL'; name a compiled pack with a material record instead",
        )),
    }
}

/// A plate section from a record's card at `point`, through the same
/// requirement-driven resolution the plate material models use: isotropic
/// `density`, `young_modulus`, `poisson_ratio`; orthotropic `density`,
/// `young_modulus_1`, `young_modulus_2`, `shear_modulus_12`,
/// `poisson_ratio_12` (`Orthotropic21` exchanges the card's axes 1 and 2 with
/// Poisson reciprocity). An isotropic section is built exactly as bare
/// isotropic numbers are, so equal values give a bitwise equal section.
///
/// # Errors
/// A missing, ambiguous or out-of-domain property refuses by name, as does a
/// section the plate owner refuses.
pub fn plate_section(
    card: &MaterialCard,
    point: &QueryPoint,
    model: PlateMaterialModel,
    thickness_m: f64,
    angle_rad: f64,
) -> Result<(PlateSection, ResolvedMaterialStatePoint), MaterialIngressError> {
    let requirement = |name: &str, dims, admissibility| {
        ScalarPropertyRequirement::try_new(name, dims, admissibility)
            .map_err(|e| refuse(format!("{e:?}")))
    };
    let mut required = vec![requirement(
        DENSITY_PROPERTY,
        Density::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?];
    match model {
        PlateMaterialModel::Isotropic => {
            required.push(requirement(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
            required.push(requirement(
                POISSON_RATIO_PROPERTY,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: -1.0,
                    upper: 0.5,
                },
            )?);
        }
        PlateMaterialModel::Orthotropic12 | PlateMaterialModel::Orthotropic21 => {
            for name in [
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[0],
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[1],
            ] {
                required.push(requirement(
                    name,
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )?);
            }
            required.push(requirement(
                ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[0],
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
            required.push(requirement(
                ORTHOTROPIC_POISSON_RATIO_PROPERTIES[0],
                Dims::NONE,
                ScalarAdmissibility::Finite,
            )?);
        }
    }
    let state = resolve_material_state_point(
        card,
        point,
        &required,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|e| refuse(format!("{e:?}")))?;
    let owner = |e: crate::acoustic_realize::AcousticRealizeError| refuse(e.to_string());
    let density = plate_property(&state, DENSITY_PROPERTY, Density::DIMS).map_err(owner)?;
    let section = match model {
        PlateMaterialModel::Isotropic => {
            if angle_rad != 0.0 {
                return Err(refuse("an isotropic section takes no material angle"));
            }
            PlateSection::isotropic(
                plate_property(&state, YOUNG_MODULUS_PROPERTY, Pressure::DIMS).map_err(owner)?,
                plate_property(&state, POISSON_RATIO_PROPERTY, Dims::NONE).map_err(owner)?,
                thickness_m,
                density,
            )
        }
        PlateMaterialModel::Orthotropic12 | PlateMaterialModel::Orthotropic21 => {
            let (e1, e2, nu12, g12) = plate_elasticity(&state, model).map_err(owner)?;
            PlateSection::orthotropic_plane_stress_at_angle(
                e1,
                e2,
                nu12,
                g12,
                thickness_m,
                density,
                angle_rad,
            )
        }
    }
    .map_err(|e| refuse(format!("{e:?}")))?;
    Ok((section, state))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn records_parse_exactly_and_refuse_by_name() {
        let record = MaterialRecord::parse(&[
            "packs/cane.fsmatpk",
            "water-soaked",
            "T=293.15",
            "moisture_content=0.3",
        ])
        .unwrap();
        assert_eq!(record.pack, "packs/cane.fsmatpk");
        assert_eq!(record.state, "water-soaked");
        assert_eq!(
            record.axes,
            vec![
                ("T".to_string(), 293.15),
                ("moisture_content".to_string(), 0.3)
            ]
        );
        let axes = MaterialRecord::parse(&[
            "p", "s", "a=1", "b=2", "c=3", "d=4", "e=5", "f=6", "g=7", "h=8",
        ]);
        assert!(axes.is_ok(), "8 axes is the cap");
        for (tokens, needle) in [
            (vec!["p"], "PACK STATE"),
            (vec!["../p", "s"], "relative"),
            (vec!["a//p", "s"], "relative"),
            (vec!["p", "s", "T"], "AXIS=VALUE"),
            (vec!["p", "s", "T=x"], "numeric"),
            (vec!["p", "s", "T=inf"], "finite"),
            (vec!["p", "s", "T=1", "T=2"], "appear once"),
            (vec!["p", "s\"q"], "quotes"),
            (
                vec![
                    "p", "s", "a=1", "b=2", "c=3", "d=4", "e=5", "f=6", "g=7", "h=8", "i=9",
                ],
                "at most 8",
            ),
        ] {
            let error = MaterialRecord::parse(&tokens).unwrap_err().0;
            assert!(error.contains(needle), "{tokens:?}: {error}");
        }
        assert_eq!(
            fixture_label(&["fixture", "authored-blade"]).unwrap(),
            "authored-blade"
        );
        for tokens in [
            vec![],
            vec!["fixture"],
            vec!["label", "x"],
            vec!["fixture", "a", "b"],
            vec!["fixture", "a\\b"],
        ] {
            assert!(fixture_label(&tokens).is_err(), "{tokens:?}");
        }
    }
}
