//! One source-law selection for characteristic networks and spatial DG bores.
//!
//! The bore owner changes neither the windway's material law nor the jet's
//! admissible profile. These records attach the original owners in the same
//! order, retaining their incompatibility, gas and validity checks.

use super::{FlueInputError, Rows, bad, poll};
use crate::bernoulli_aperture::contraction::ContractionSpec;
use crate::flue_loop::forced_jet::ForcedJetSpec;
use crate::flue_loop::nolle_jet::{ChannelClass, NolleSpec};
use crate::flue_loop::{FlueError, FlueVoice};
use fs_exec::CancelGate;

enum WindwayInput {
    None,
    QuasiSteady {
        length_m: f64,
    },
    Developing {
        length_m: f64,
    },
    Unsteady {
        length_m: f64,
        transverse_cells: usize,
        convection_cells: usize,
    },
}

/// Optional records between the physical flue geometry and the jet card.
pub(crate) struct FlueFlowInput {
    windway: WindwayInput,
    mouth_loss: Option<ContractionSpec>,
}

impl FlueFlowInput {
    pub(crate) fn read(rows: &mut Rows<'_>) -> Result<Self, FlueInputError> {
        let windway = if rows.next_is("windway") {
            let mut row = rows.row("windway")?;
            let selection = match row.fields.next() {
                Some("none") => WindwayInput::None,
                Some("slit") => WindwayInput::QuasiSteady {
                    length_m: row.scalar()?,
                },
                Some("slit-developing") => WindwayInput::Developing {
                    length_m: row.scalar()?,
                },
                Some("unsteady") => WindwayInput::Unsteady {
                    length_m: row.scalar()?,
                    transverse_cells: row.number()?,
                    convection_cells: row.number()?,
                },
                _ => {
                    return Err(bad(
                        row.line,
                        "windway requires none, slit LENGTH, slit-developing LENGTH, or unsteady LENGTH TRANSVERSE_CELLS CONVECTION_CELLS",
                    ));
                }
            };
            row.end()?;
            selection
        } else {
            WindwayInput::None
        };
        let mouth_loss = if rows.next_is("mouth_loss") {
            let mut row = rows.row("mouth_loss")?;
            let selection = match row.fields.next() {
                Some("none") => None,
                Some("contraction") => Some(ContractionSpec {
                    area_ratio: row.scalar()?,
                    maximum_contracted_mach: row.scalar()?,
                    provenance: row
                        .fields
                        .next()
                        .ok_or_else(|| {
                            bad(row.line, "mouth contraction needs an explicit source token")
                        })?
                        .to_string(),
                }),
                _ => {
                    return Err(bad(
                        row.line,
                        "mouth_loss requires none or contraction AREA_RATIO MAXIMUM_MACH SOURCE_TOKEN",
                    ));
                }
            };
            row.end()?;
            selection
        } else {
            None
        };
        Ok(Self {
            windway,
            mouth_loss,
        })
    }

    pub(crate) fn attach(
        self,
        voice: FlueVoice,
        held_pressure_pa: f64,
    ) -> Result<FlueVoice, FlueInputError> {
        let voice = match self.mouth_loss {
            Some(spec) => voice
                .with_mouth_contraction(spec)
                .map_err(FlueInputError::Physics)?,
            None => voice,
        };
        match self.windway {
            WindwayInput::None => Ok(voice),
            WindwayInput::QuasiSteady { length_m } => voice.with_viscous_windway(length_m),
            WindwayInput::Developing { length_m } => {
                voice.with_developing_windway(length_m, held_pressure_pa)
            }
            WindwayInput::Unsteady {
                length_m,
                transverse_cells,
                convection_cells,
            } => voice.with_unsteady_windway(length_m, transverse_cells, convection_cells),
        }
        .map_err(FlueInputError::Physics)
    }
}

/// Optional source selections after `jet`, before `jet_limits`.
pub(crate) struct FlueJetInput {
    card_centreline: bool,
    transfer: Option<NolleSpec>,
    response: Option<ForcedJetSpec>,
}

impl FlueJetInput {
    pub(crate) fn read(rows: &mut Rows<'_>) -> Result<Self, FlueInputError> {
        let card_centreline = if rows.next_is("jet_speed") {
            let mut row = rows.row("jet_speed")?;
            if row.fields.next() != Some("card-centreline") {
                return Err(bad(row.line, "jet_speed requires card-centreline"));
            }
            row.end()?;
            true
        } else {
            false
        };
        let transfer = if rows.next_is("jet_transfer") {
            let mut row = rows.row("jet_transfer")?;
            if row.fields.next() != Some("nolle-price-2015") {
                return Err(bad(
                    row.line,
                    "jet_transfer requires nolle-price-2015 CHANNEL CHAMFER_DEPTH_M",
                ));
            }
            let channel = match row.fields.next() {
                Some("long-channel") => ChannelClass::Long,
                Some("short-channel") => ChannelClass::Short,
                _ => {
                    return Err(bad(
                        row.line,
                        "jet_transfer channel must be long-channel or short-channel",
                    ));
                }
            };
            let chamfer_depth_m = row.scalar()?;
            row.end()?;
            Some(
                NolleSpec::price_2015(channel, chamfer_depth_m)
                    .map_err(|what| FlueInputError::Physics(FlueError::Invalid { what }))?,
            )
        } else {
            None
        };
        let response = if rows.next_is("jet_response") {
            let mut row = rows.row("jet_response")?;
            if row.fields.next() != Some("forced-bickley") {
                return Err(bad(
                    row.line,
                    "jet_response requires forced-bickley SPACING WINDOW",
                ));
            }
            let spec = ForcedJetSpec {
                spacing: row.scalar()?,
                window: row.scalar()?,
            };
            row.end()?;
            Some(spec)
        } else {
            None
        };
        if transfer.is_some() && response.is_some() {
            return Err(bad(
                rows.line,
                "jet_transfer and jet_response exclude each other",
            ));
        }
        Ok(Self {
            card_centreline,
            transfer,
            response,
        })
    }

    pub(crate) fn attach(
        self,
        voice: FlueVoice,
        held_pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<FlueVoice, FlueInputError> {
        let voice = if self.card_centreline {
            voice
                .with_card_centreline_speed(held_pressure_pa)
                .map_err(FlueInputError::Physics)?
        } else {
            voice
        };
        let voice = match self.transfer {
            Some(spec) => {
                poll(gate)?;
                let voice = voice
                    .with_nolle_jet(spec)
                    .map_err(FlueInputError::Physics)?;
                poll(gate)?;
                voice
            }
            None => voice,
        };
        match self.response {
            Some(spec) => {
                poll(gate)?;
                let voice = voice
                    .with_forced_jet(spec, held_pressure_pa)
                    .map_err(FlueInputError::Physics)?;
                poll(gate)?;
                Ok(voice)
            }
            None => Ok(voice),
        }
    }
}
