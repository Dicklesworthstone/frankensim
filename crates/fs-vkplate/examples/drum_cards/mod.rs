//! The kit's default build in the offline examples (.370): heads and stick
//! from sourced material cards, no recording-fitted head T60, a head film loss
//! from the film card's derived loss-tangent curve (ruling 7468), snare wires
//! that are coils of sourced steel (ruling 7373), and in place of every value
//! adopted because it graded closer to a recording (96ccf8f1b), the
//! recording's pitch for the tensions and 88cccbb78's declared estimates for
//! the rest (ruling 7486). The crash struck with the shaft takes its stick,
//! and the wood's crush across the grain, from the same hickory card (ruling
//! 7476).

use fs_matdb::QueryPoint;
use fs_qty::{Dims, QuantitySpec};
use fs_vkplate::kit::{self, CymbalSpec, DrumSpec, Piece};
use fs_vkplate::materials::{KitMaterials, temperature_frequency_point, temperature_point};
use fs_vkplate::shaft::{SHOULDER_STATION_M, ShaftStrike};

/// The PET head-film card and the hickory stick card, both compiled from
/// data/matdb/seed-v1 (see the crate README).
const HEAD_FILM: &[u8] =
    include_bytes!("../../packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const HICKORY: &[u8] = include_bytes!("../../packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");

fn hickory_point(temperature_k: f64) -> Result<QueryPoint, String> {
    temperature_point("T", temperature_k)?
        .with_quantity("moisture_content", QuantitySpec::dimensional(Dims::NONE), 0.12)
        .map_err(|error| error.to_string())
}

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

/// The film card's loss-tangent curve lies on its typed `frequency` axis.
const FILM_FREQUENCY_AXIS: &str = "frequency";

/// What the film card's curve is, for the drums' notes (the card's
/// hitachi-dma9 observation states it in full).
pub const FILM_LOSS_SOURCE_NOTE: &str = "head film loss source: DERIVED on the PET film card from Hitachi High-Tech DMA brief No. 9 (a 10 um annealed biaxially oriented film, not head gauge): its 1 Hz tan delta trace, shifted to 20 C by the beta relaxation's Arrhenius law (Ea 76 kJ/mol), 20 Hz to 20 kHz, about +/- 30 %";

/// The density at which kit::STICK_TIP_KG, the Z5A profile's contact-effective
/// mass, was integrated (crates/fs-couple/examples/percussion). The mass over
/// this density is the profile's own effective volume, so the hickory card's
/// density sets the stick's moving mass without rescaling its geometry.
const Z5A_DENSITY_ESTIMATE_KG_M3: f64 = 800.0;

/// The stick tip's radius (as the presets' kit::hertz_stick) and its
/// Hunt-Crossley loss, 88cccbb78's declared estimate (ruling 7486). Neither is
/// a card property.
const TIP_RADIUS_M: f64 = 0.005;
const TIP_LOSS_S_M: f64 = kit::DECLARED_STICK_LOSS_S_M;

pub const PRESET_HEADS_FLAG: &str = "--preset-heads";
pub const PRESET_HEADS_RECEIPT: &str =
    "heads: PRESET, recording-fitted T60 (opt-in --preset-heads); not a sourced head loss";

/// Builds a drum from a preset spec, with both head tensions scaled by `tune`
/// (a drummer's tuning turn; 1 is the drum's own tuning). By default the drum
/// is `kit::card_drum(spec)`, whose tensions, strike point, wires and kick
/// beater replace the values fitted to a recording (ruling 7486); its heads
/// come from the film card and its stick from the hickory card, at the
/// standard gas's temperature. The heads' film loss is the film card's
/// loss-tangent curve, read at each head mode's own frequency; its receipt,
/// FILM_LOSS_SOURCE_NOTE and the edge loss's stated unknown join the notes.
/// The kick keeps the spec's felt beater, card_drum's declared estimate, and
/// says so (ruling 6433 (b)); its notes also say that no pillow is modelled.
/// A drum with snare wires binds the wire steel cards, so its strands are
/// coils. With `preset_heads` the drum is the recording-fitted preset, and
/// PRESET_HEADS_RECEIPT joins its notes.
pub fn drum(
    spec: &DrumSpec,
    tune: f64,
    hits: &[(f64, f64)],
    kick: bool,
    preset_heads: bool,
) -> Result<Piece, String> {
    if preset_heads {
        let mut piece = kit::drum(&kit::retune(*spec, tune), hits)?;
        piece.notes.push(PRESET_HEADS_RECEIPT.into());
        return Ok(piece);
    }
    let spec = &kit::retune(kit::card_drum(spec)?, tune);
    let gas = fs_vkplate::standard_gas();
    let point = hickory_point(gas.temperature)?;
    let cards = KitMaterials::from_card_packs(
        HEAD_FILM,
        &[],
        HICKORY,
        &point,
        kit::STICK_TIP_KG / Z5A_DENSITY_ESTIMATE_KG_M3,
        TIP_RADIUS_M,
        TIP_LOSS_S_M,
    )
    .and_then(|m| m.with_film_loss_curve("T", FILM_FREQUENCY_AXIS, gas.temperature))?;
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
    let mut piece = if kick {
        kit::kick_with_materials(spec, hits, gas, &cards)
    } else {
        kit::drum_with_materials(spec, hits, gas, &cards)
    }?;
    piece.notes.push(FILM_LOSS_SOURCE_NOTE.into());
    Ok(piece)
}

/// The ride the examples play: by default `kit::card_ride()`, 88cccbb78's
/// declared shell (ruling 7486); with `preset` the recording-fitted preset.
pub fn ride(preset: bool) -> CymbalSpec {
    if preset {
        kit::ride()
    } else {
        kit::card_ride()
    }
}

/// A cymbal struck with 88cccbb78's declared stick loss
/// (kit::DECLARED_CYMBAL_INTERFACE), or with `preset` the presets' fitted one.
pub fn cymbal(
    spec: &CymbalSpec,
    hits: &[(f64, f64, f64, f64)],
    preset: bool,
) -> Result<Piece, String> {
    if preset {
        kit::cymbal(spec, hits)
    } else {
        kit::cymbal_with_interface(spec, hits, &kit::DECLARED_CYMBAL_INTERFACE)
    }
}

/// The hi-hat with 88cccbb78's declared rim law, flat rims and 5 mm stick tip
/// (kit::DECLARED_HIHAT_RIM_LAW, DECLARED_CYMBAL_INTERFACE), or with `preset`
/// the presets' fitted ones.
pub fn hihat(gap: f64, pedal_n: f64, hits: &[(f64, f64)], preset: bool) -> Result<Piece, String> {
    if preset {
        kit::hihat(gap, pedal_n, hits)
    } else {
        kit::hihat_with_interface(
            gap,
            pedal_n,
            hits,
            kit::DECLARED_HIHAT_RIM_LAW,
            &kit::DECLARED_CYMBAL_INTERFACE,
        )
    }
}

/// A preset-bronze cymbal struck once at its edge, at `speed` m/s, by the
/// hickory card's shaft at the declared shoulder station, the wood crushing
/// across its grain. The edge is taken fully rounded (half its thickness).
/// There is no default speed: a shaft strike is rendered as a ladder of speeds
/// (ruling 7373).
#[allow(dead_code)] // perform.rs shares this module and strikes no cymbal
pub fn shaft_cymbal(spec: &kit::CymbalSpec, speed: Option<f64>) -> Result<Piece, String> {
    let speed = speed.ok_or(
        "crash-shaft has no default speed: give SPEED_M_S, the shaft's speed at the edge (render a ladder of speeds, ruling 7373)",
    )?;
    let gas = fs_vkplate::standard_gas();
    let shaft = ShaftStrike::from_card(
        HICKORY,
        &hickory_point(gas.temperature)?,
        SHOULDER_STATION_M,
        0.5 * kit::edge_thickness(spec),
        (kit::BRONZE.0, kit::BRONZE.1),
    )?;
    kit::cymbal_struck_by(
        spec,
        &[(0.005, speed, 1.0)],
        gas,
        (shaft.mass_kg, shaft.law, Some(shaft.crush)),
        &shaft.notes,
    )
}

/// The open hat's rim opening on the card route: 88cccbb78's declared 8 mm (ruling 7520 G2).
pub const OPEN_HAT_OPENING_M: f64 = 0.008;
pub const OPEN_HAT_OPENING_NOTE: &str = "open hi-hat: rim opening 8 mm, the declared opening (88cccbb78); the recording's pedal position is unknown";
/// The presets' open-hat opening, picked by grade (bb1e5abf3).
pub const PRESET_OPEN_HAT_OPENING_M: f64 = 0.002;
pub const PRESET_OPEN_HAT_OPENING_NOTE: &str = "open hi-hat: rim opening 2 mm, PRESET, picked by grade (bb1e5abf3); not a declared or measured opening";

/// The open hi-hat, no pedal force: `hihat` at the declared opening, or with `preset` at the
/// presets' graded one. Its notes say which.
pub fn open_hihat(hits: &[(f64, f64)], preset: bool) -> Result<Piece, String> {
    let (gap, note) = if preset {
        (PRESET_OPEN_HAT_OPENING_M, PRESET_OPEN_HAT_OPENING_NOTE)
    } else {
        (OPEN_HAT_OPENING_M, OPEN_HAT_OPENING_NOTE)
    };
    let mut piece = hihat(gap, 0.0, hits, preset)?;
    piece.notes.push(note.into());
    Ok(piece)
}
