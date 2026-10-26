//! The drums' default build in the offline examples (.370): heads and stick
//! from sourced material cards, no recording-fitted head T60, and snare wires
//! that are coils of sourced steel (ruling 7373).

use fs_matdb::QueryPoint;
use fs_vkplate::kit::{self, DrumSpec, Piece};
use fs_vkplate::materials::{KitMaterials, temperature_frequency_point};

/// The PET head-film card and the hickory stick card, both compiled from
/// data/matdb/seed-v1 (see the crate README).
const HEAD_FILM: &[u8] =
    include_bytes!("../../packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const HICKORY: &[u8] = include_bytes!("../../packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");

/// The snare wires' steel: the music-wire string card (Heptonstall, Waller and
/// Robertson, LIGO P1400229) and, for its shear modulus, the music-wire spring
/// moduli card (Fuchs 1968). Both compiled from data/matdb/seed-v1.
const WIRE_STEEL: &[u8] =
    include_bytes!("../../packs/steel-music-wire-string-heptonstall-ligo-p1400229.fsmatpk");
const WIRE_SHEAR: &[u8] = include_bytes!("../../packs/music-wire-spring-moduli-fuchs-1968.fsmatpk");

/// The steel card's loss tangent holds from 8 to 3000 Hz on its typed
/// `frequency` axis, so the card is read at 1 kHz inside that band, as the
/// piano reads it; all its claims are constant there.
const WIRE_STEEL_FREQUENCY: (&str, f64) = ("frequency", 1000.0);

/// The shear card's only state: its untyped `temperature` axis at 293 K, the
/// one temperature it is valid at (the kit's gas is at 293.15 K).
const WIRE_SHEAR_STATE: (&str, f64) = ("temperature", 293.0);

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
/// (ruling 6433 (b)). A drum with snare wires binds the wire steel cards, so
/// its strands are coils. With `preset_heads` the drum is the recording-fitted
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
    let cards = if spec.wires.is_some() {
        let shear = QueryPoint::new()
            .with(WIRE_SHEAR_STATE.0, WIRE_SHEAR_STATE.1)
            .map_err(|e| e.to_string())?;
        cards.with_snare_wire(
            WIRE_STEEL,
            &temperature_frequency_point(
                "T",
                gas.temperature,
                WIRE_STEEL_FREQUENCY.0,
                WIRE_STEEL_FREQUENCY.1,
            )?,
            WIRE_SHEAR,
            &shear,
        )?
    } else {
        cards
    };
    if kick {
        kit::kick_with_materials(spec, hits, gas, &cards)
    } else {
        kit::drum_with_materials(spec, hits, gas, &cards)
    }
}
