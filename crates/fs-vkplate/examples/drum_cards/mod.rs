//! The drums' default build in the offline examples (.370): heads and stick
//! from sourced material cards, no recording-fitted head T60.

use fs_vkplate::kit::{self, DrumSpec, Piece};
use fs_vkplate::materials::KitMaterials;

/// The PET head-film card and the hickory stick card, both compiled from
/// data/matdb/seed-v1 (see the crate README).
const HEAD_FILM: &[u8] =
    include_bytes!("../../packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const HICKORY: &[u8] = include_bytes!("../../packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");

/// The density at which kit::STICK_TIP_KG, the Z5A profile's contact-effective
/// mass, was integrated (crates/fs-couple/examples/percussion). The mass over
/// this density is the profile's own effective volume, so the hickory card's
/// density sets the stick's moving mass without rescaling its geometry.
const Z5A_DENSITY_ESTIMATE_KG_M3: f64 = 800.0;

/// The presets' declared stick tip (kit::hertz_stick): radius and
/// Hunt-Crossley loss. Neither is a card property.
const TIP_RADIUS_M: f64 = 0.005;
const TIP_LOSS_S_M: f64 = 0.05;

pub const PRESET_HEADS_FLAG: &str = "--preset-heads";
pub const PRESET_HEADS_RECEIPT: &str =
    "heads: PRESET, recording-fitted T60 (opt-in --preset-heads); not a sourced head loss";

/// Builds a drum. By default its heads come from the film card and its stick
/// from the hickory card, at the standard gas's temperature. The film card
/// claims no loss, so the film loss is a declared stated unknown, whose
/// receipt joins the notes. The kick keeps its preset felt beater and says so
/// (ruling 6433 (b)). With `preset_heads` the drum is the recording-fitted
/// preset, and PRESET_HEADS_RECEIPT joins its notes.
pub fn drum(
    spec: &DrumSpec,
    hits: &[(f64, f64)],
    kick: bool,
    preset_heads: bool,
) -> Result<Piece, String> {
    if preset_heads {
        let mut piece = kit::drum(spec, hits)?;
        piece.notes.push(PRESET_HEADS_RECEIPT.into());
        return Ok(piece);
    }
    let gas = fs_vkplate::standard_gas();
    let cards = KitMaterials::from_packs_at_temperature(
        HEAD_FILM,
        &[],
        HICKORY,
        "T",
        gas.temperature,
        kit::STICK_TIP_KG / Z5A_DENSITY_ESTIMATE_KG_M3,
        TIP_RADIUS_M,
        TIP_LOSS_S_M,
    )
    .and_then(KitMaterials::with_stated_unknown_film_loss)?;
    if kick {
        kit::kick_with_materials(spec, hits, gas, &cards)
    } else {
        kit::drum_with_materials(spec, hits, gas, &cards)
    }
}
